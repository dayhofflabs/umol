//! Tree-shaped constraint DSLs.
//!
//! Boundary types between the AST `Constraint` tree and its EDN form. Refs in
//! the tree carry either an integer index or a symbolic id; resolution to /
//! from the `AtomId` / `BondId` / ... on the AST is a separate fallible
//! step that consults the surrounding `Metadata`.

use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::aromatic::AromaticSystemConstraintDsl;
use super::atom::{AromaticValenceDsl, AtomConstraintDsl, MulticenterValenceDsl};
use super::bond::BondConstraintDsl;
use super::dative::DativeBondConstraintDsl;
use super::error::ParseError;
use super::molecule::{Metadata, MoleculeDsl};
use super::multicenter::MulticenterBondConstraintDsl;
use super::noncovalent::NoncovalentBondConstraintDsl;
use super::relational::{RelationalConstraintDsl, RELATIONAL_KEYS};
use super::stereo::{
    coset_lit, parse_stereo_coset, StereoAtomConstraintDsl, StereoBondConstraintDsl, StereoCosetDsl,
};
use super::value::{parse_value, ValueDsl};
use crate::ast::constraint::{
    AromaticValenceAst, AtomConstraint, BondConstraint, Constraint, Constraints,
    MoleculeConstraint, MulticenterValenceAst, RingMembershipAst, RingScope, SubPatternAnchor,
};
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::spin::SpinStateAst;
use crate::ast::stereo::{CisTransStereoAst, StereoCosetAst, StereoKind, TetrahedralStereoAst};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Per-entity counts for numeric-index bounds checking during constraint
/// resolution (DSL → AST). `from_ast` (AST → DSL) does not read counts.
///
/// Crate-internal. `Copy` (48 B of primitives) — callers pass by reference
/// for consistency with `Metadata`, not for cost reasons.
#[derive(Debug, Clone, Copy)]
pub(crate) struct EntityCounts {
    pub(crate) atom_count: usize,
    pub(crate) bond_count: usize,
    pub(crate) dative_bond_count: usize,
    pub(crate) aromatic_system_count: usize,
    pub(crate) multicenter_bond_count: usize,
    pub(crate) noncovalent_bond_count: usize,
    pub(crate) stereo_atom_count: usize,
    pub(crate) stereo_bond_count: usize,
}

impl EntityCounts {
    pub(crate) fn from_ast(ast: &MoleculeAst) -> Self {
        Self {
            atom_count: ast.atoms().count(),
            bond_count: ast.bonds().count(),
            dative_bond_count: ast.dative_bonds().count(),
            aromatic_system_count: ast.aromatic_systems().count(),
            multicenter_bond_count: ast.multicenter_bonds().count(),
            noncovalent_bond_count: ast.noncovalent_bonds().count(),
            stereo_atom_count: ast.stereo_atoms().count(),
            stereo_bond_count: ast.stereo_bonds().count(),
        }
    }
}

// Shared across the constraint-tree readers here and the molecule-map reader
// in `super::molecule`.

pub(super) fn eof_err() -> EdnError {
    DeError::Custom("unexpected end of input".into()).into()
}

pub(super) fn missing(key: &str, context: &'static str) -> EdnError {
    DeError::MissingField {
        key: key.to_string(),
        path: vec![context.into()],
    }
    .into()
}

pub(super) fn unexpected_byte_kind(b: u8) -> &'static str {
    match b {
        b'"' => "string",
        b':' => "keyword",
        b'[' => "vector",
        b'{' => "map",
        b'0'..=b'9' | b'-' | b'+' => "number",
        _ => "token",
    }
}

pub(super) fn read_vec<T>(
    de: &mut EdnStreamDeserializer<'_>,
    mut read_element: impl FnMut(&mut EdnStreamDeserializer<'_>) -> Result<T, EdnError>,
) -> Result<Vec<T>, EdnError> {
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            break;
        }
        out.push(read_element(de)?);
    }
    Ok(out)
}

pub(super) fn read_map(
    de: &mut EdnStreamDeserializer<'_>,
    mut on_entry: impl FnMut(&mut EdnStreamDeserializer<'_>, &str) -> Result<(), EdnError>,
) -> Result<(), EdnError> {
    de.consume_byte(b'{')?;
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?.into_owned();
        on_entry(de, key.as_str())?;
    }
    Ok(())
}

/// Consume `{:key value}` as a single-key map, returning the key and
/// leaving the stream positioned at the opening-map byte (caller has already
/// read the value). Errors if the map contains more than one key.
pub(super) fn read_single_key_map_header(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<String, EdnError> {
    de.consume_byte(b'{')?;
    Ok(de.read_keyword_name()?.into_owned())
}

pub(super) fn consume_single_key_map_close(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<(), EdnError> {
    if !de.try_consume_byte(b'}')? {
        return Err(DeError::Custom(format!("{} must have exactly one key", context)).into());
    }
    Ok(())
}

macro_rules! define_ref {
    ($name:ident, $id:ident, $accessor:ident, $kind:literal, $reader:ident) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
        }

        impl $name {
            /// Build a ref from an AST index, preferring an id from `metadata`
            /// if one is recorded for this index.
            pub fn from_ast(id: $id, metadata: &Metadata) -> Self {
                if let Some(name) = metadata.$accessor(id) {
                    Self::Id(name.to_string())
                } else {
                    Self::Index(id.index())
                }
            }

            /// Resolve this ref to an AST index against `metadata`. Fails on
            /// unknown id or out-of-range numeric index.
            pub fn into_ast(self, count: usize, metadata: &Metadata) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => {
                        for i in 0..count {
                            let id = $id::from(i);
                            if metadata.$accessor(id) == Some(name.as_str()) {
                                return Ok(id);
                            }
                        }
                        Err(ParseError::InvalidRef {
                            kind: $kind,
                            value: name,
                        })
                    }
                }
            }

            /// Resolve this ref against a pre-built id → index map. O(1) id
            /// lookup; intended for entity-loop resolution where cloning the
            /// full `Metadata` per call is wasteful.
            pub fn resolve(
                self,
                count: usize,
                id_to_idx: &IndexMap<String, $id>,
            ) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => id_to_idx.get(&name).copied().ok_or(ParseError::InvalidRef {
                        kind: $kind,
                        value: name,
                    }),
                }
            }
        }

        impl<'de> FromEdn<'de> for $name {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                match edn {
                    Edn::Int(n) => {
                        let i = usize::try_from(*n).map_err(|_| DeError::OutOfRange {
                            value: n.to_string(),
                            target: "usize",
                            path: Vec::new(),
                        })?;
                        Ok(Self::Index(i))
                    }
                    Edn::Keyword(k) => Ok(Self::Id(k.name().to_string())),
                    other => Err(DeError::TypeMismatch {
                        expected: concat!($kind, " ref (int or keyword)"),
                        got: other.kind(),
                        path: Vec::new(),
                    }),
                }
            }
        }

        impl ToEdn for $name {
            fn to_edn(&self) -> Edn<'static> {
                match self {
                    Self::Index(i) => Edn::Int(*i as i64),
                    Self::Id(name) => Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone())),
                }
            }
        }

        pub(super) fn $reader(de: &mut EdnStreamDeserializer<'_>) -> Result<$name, EdnError> {
            match de.peek_byte()?.ok_or_else(eof_err)? {
                b':' => Ok($name::Id(de.read_keyword_name()?.into_owned())),
                _ => {
                    let n = de.read_i64()?;
                    let i = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                        value: n.to_string(),
                        target: "usize",
                        path: Vec::new(),
                    })?;
                    Ok($name::Index(i))
                }
            }
        }
    };
}

define_ref!(AtomRef, AtomId, atom_id, "atom", read_atom_ref);
define_ref!(BondRef, BondId, bond_id, "bond", read_bond_ref);
define_ref!(
    DativeBondRef,
    DativeBondId,
    dative_bond_id,
    "dative-bond",
    read_dative_bond_ref
);
define_ref!(
    AromaticSystemRef,
    AromaticSystemId,
    aromatic_system_id,
    "aromatic-system",
    read_aromatic_system_ref
);
define_ref!(
    MulticenterBondRef,
    MulticenterBondId,
    multicenter_bond_id,
    "multicenter-bond",
    read_multicenter_bond_ref
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondId,
    noncovalent_bond_id,
    "noncovalent-bond",
    read_noncovalent_bond_ref
);
define_ref!(
    StereoAtomRef,
    StereoAtomId,
    stereo_atom_id,
    "stereo-atom",
    read_stereo_atom_ref
);
define_ref!(
    StereoBondRef,
    StereoBondId,
    stereo_bond_id,
    "stereo-bond",
    read_stereo_bond_ref
);

pub(super) fn read_value_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<ValueDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'"' => {
            let s = de.read_string()?;
            let v: ValueAst =
                parse_value(s.as_ref()).map_err(|e| DeError::subgrammar("value", e))?;
            Ok(ValueDsl(v))
        }
        b'[' => {
            let items = read_vec(de, |d| Ok(d.read_i64()?))?;
            Ok(ValueDsl(ValueAst::lit_set(items)))
        }
        b':' => {
            let name = de.read_keyword_name()?;
            if name.as_ref() == "undetermined" {
                Ok(ValueDsl(ValueAst::Undetermined))
            } else {
                Err(
                    DeError::Custom(format!("unexpected keyword :{} in value position", name))
                        .into(),
                )
            }
        }
        _ => Ok(ValueDsl(ValueAst::Lit(de.read_i64()?))),
    }
}

pub(super) fn read_spin_state(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<SpinStateAst, EdnError> {
    let mut unpaired = None;
    let mut multiplicity = None;
    read_map(de, |d, key| {
        match key {
            "unpaired" => unpaired = Some(read_value_dsl(d)?.into_ast(&())),
            "multiplicity" => multiplicity = Some(read_value_dsl(d)?.into_ast(&())),
            _ => d.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(SpinStateAst {
        unpaired: unpaired.ok_or_else(|| missing("unpaired", "spin"))?,
        multiplicity: multiplicity.ok_or_else(|| missing("multiplicity", "spin"))?,
    })
}

pub(super) fn read_aromatic_valence_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticValenceDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b':' => {
            let name = de.read_keyword_name()?;
            match name.as_ref() {
                "undetermined" => Ok(AromaticValenceDsl(AromaticValenceAst::Undetermined)),
                "not-aromatic" => Ok(AromaticValenceDsl(AromaticValenceAst::NotAromatic)),
                other => Err(DeError::Custom(format!(
                    "unknown aromatic-valence keyword :{}",
                    other
                ))
                .into()),
            }
        }
        b'{' => {
            let key = read_single_key_map_header(de)?;
            match key.as_str() {
                "aromatic" => {
                    let v = read_value_dsl(de)?.into_ast(&());
                    consume_single_key_map_close(de, "aromatic-valence")?;
                    Ok(AromaticValenceDsl(AromaticValenceAst::Aromatic(v)))
                }
                other => Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["aromatic-valence".into()],
                }
                .into()),
            }
        }
        b => Err(DeError::TypeMismatch {
            expected: ":undetermined / :not-aromatic / {:aromatic <value>}",
            got: unexpected_byte_kind(b),
            path: vec!["aromatic-valence".into()],
        }
        .into()),
    }
}

pub(super) fn read_multicenter_valence_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterValenceDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b':' => {
            let name = de.read_keyword_name()?;
            match name.as_ref() {
                "undetermined" => Ok(MulticenterValenceDsl(MulticenterValenceAst::Undetermined)),
                "not-multicenter" => {
                    Ok(MulticenterValenceDsl(MulticenterValenceAst::NotMulticenter))
                }
                other => Err(DeError::Custom(format!(
                    "unknown multicenter-valence keyword :{}",
                    other
                ))
                .into()),
            }
        }
        b'{' => {
            let key = read_single_key_map_header(de)?;
            match key.as_str() {
                "multicenter" => {
                    let v = read_value_dsl(de)?.into_ast(&());
                    consume_single_key_map_close(de, "multicenter-valence")?;
                    Ok(MulticenterValenceDsl(MulticenterValenceAst::Multicenter(v)))
                }
                other => Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["multicenter-valence".into()],
                }
                .into()),
            }
        }
        b => Err(DeError::TypeMismatch {
            expected: ":undetermined / :not-multicenter / {:multicenter <value>}",
            got: unexpected_byte_kind(b),
            path: vec!["multicenter-valence".into()],
        }
        .into()),
    }
}

/// Streaming-read the `{:size? :count}` map (value of a `:ring-membership` key).
fn read_ring_membership_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<RingMembershipAst, EdnError> {
    let mut size: Option<u8> = None;
    let mut count: Option<ValueAst> = None;
    read_map(de, |de, key| {
        match key {
            "size" => size = Some(de.read_i64()? as u8),
            "count" => count = Some(read_value_dsl(de)?.into_ast(&())),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["ring-membership".into()],
                }
                .into())
            }
        }
        Ok(())
    })?;
    let count =
        count.ok_or_else(|| DeError::Custom("ring-membership missing :count".to_string()))?;
    Ok(RingMembershipAst::new(
        size.map_or(RingScope::All, RingScope::Size),
        count,
    ))
}

/// EDN boundary for a ring-membership fact: `{:size? <int> :count <value>}`.
pub struct RingMembershipDsl(pub RingMembershipAst);

impl ToEdn for RingMembershipDsl {
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(2);
        if let RingScope::Size(s) = self.0.scope {
            m.insert(
                Edn::Keyword(EdnKeyword::owned("size".into())),
                Edn::Int(s as i64),
            );
        }
        m.insert(
            Edn::Keyword(EdnKeyword::owned("count".into())),
            ValueDsl::from_ast(&self.0.count, &()).to_edn(),
        );
        Edn::Map(m)
    }
}

impl<'de> FromEdn<'de> for RingMembershipDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(map) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "map",
                got: edn.kind(),
                path: vec!["ring-membership".into()],
            });
        };
        let mut scope = RingScope::All;
        let mut count = None;
        for (k, v) in map.iter() {
            let Edn::Keyword(key) = k else {
                return Err(DeError::TypeMismatch {
                    expected: "keyword",
                    got: k.kind(),
                    path: vec!["ring-membership".into()],
                });
            };
            match key.name() {
                "size" => {
                    let Edn::Int(n) = v else {
                        return Err(DeError::TypeMismatch {
                            expected: "int",
                            got: v.kind(),
                            path: vec!["ring-membership".into(), "size".into()],
                        });
                    };
                    scope = RingScope::Size(*n as u8);
                }
                "count" => count = Some(ValueDsl::from_edn(v)?.into_ast(&())),
                other => {
                    return Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["ring-membership".into()],
                    })
                }
            }
        }
        let count = count.ok_or_else(|| DeError::Custom("ring-membership missing :count".into()))?;
        Ok(Self(RingMembershipAst::new(scope, count)))
    }
}

pub(super) fn read_stereo_coset_dsl(
    de: &mut EdnStreamDeserializer<'_>,
    degree: usize,
) -> Result<StereoCosetDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'"' => {
            let s = de.read_string()?;
            let coset = parse_stereo_coset(s.as_ref(), degree)
                .map_err(|e| DeError::subgrammar("stereo coset", e))?;
            Ok(StereoCosetDsl(coset))
        }
        b'[' => {
            let items = read_vec(de, |d| Ok(d.read_i64()?))?;
            let set = items
                .into_iter()
                .map(coset_lit)
                .collect::<Result<Vec<_>, _>>()?;
            Ok(StereoCosetDsl(StereoCosetAst::LitSet(
                set.into_iter().collect(),
            )))
        }
        b':' => {
            let name = de.read_keyword_name()?;
            if name.as_ref() == "undetermined" {
                Ok(StereoCosetDsl(StereoCosetAst::Undetermined))
            } else {
                Err(
                    DeError::Custom(format!("unexpected keyword :{} in coset position", name))
                        .into(),
                )
            }
        }
        _ => Ok(StereoCosetDsl(StereoCosetAst::Lit(coset_lit(
            de.read_i64()?,
        )?))),
    }
}

/// Streaming counterpart of `stereo_site_dsl!`'s `FromEdn`: reads a fixed-kind
/// site value (`:undetermined`, `:not-stereo`, or `{:stereo <coset>}`) straight
/// from the deserializer. `$kind` fixes the coset degree.
macro_rules! read_stereo_site_dsl {
    ($name:ident, $ast:ident, $kind:expr) => {
        pub(super) fn $name(
            de: &mut EdnStreamDeserializer<'_>,
        ) -> Result<$ast, EdnError> {
            match de.peek_byte()?.ok_or_else(eof_err)? {
                b':' => {
                    let name = de.read_keyword_name()?;
                    match name.as_ref() {
                        "undetermined" => Ok($ast::Undetermined),
                        "not-stereo" => Ok($ast::NotStereo),
                        other => Err(DeError::Custom(format!(
                            "unknown stereo-configuration keyword :{}",
                            other
                        ))
                        .into()),
                    }
                }
                b'{' => {
                    let key = read_single_key_map_header(de)?;
                    match key.as_str() {
                        "stereo" => {
                            let coset = read_stereo_coset_dsl(de, $kind.degree())?.into_ast(&());
                            consume_single_key_map_close(de, "stereo-configuration")?;
                            Ok($ast::Stereo(coset))
                        }
                        other => Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["stereo-configuration".into()],
                        }
                        .into()),
                    }
                }
                b => Err(DeError::TypeMismatch {
                    expected: ":undetermined / :not-stereo / {:stereo <coset>}",
                    got: unexpected_byte_kind(b),
                    path: vec!["stereo-configuration".into()],
                }
                .into()),
            }
        }
    };
}

read_stereo_site_dsl! { read_tetrahedral_stereo_dsl, TetrahedralStereoAst, StereoKind::Tetrahedral }
read_stereo_site_dsl! { read_cis_trans_stereo_dsl, CisTransStereoAst, StereoKind::CisTrans }

pub(super) fn read_atom_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AtomConstraintDsl, EdnError> {
    let key = read_single_key_map_header(de)?;
    let c = match key.as_str() {
        "valence" => AtomConstraint::Valence(read_value_dsl(de)?.into_ast(&())),
        "total-valence" => AtomConstraint::TotalValence(read_value_dsl(de)?.into_ast(&())),
        "aromatic-valence" => {
            AtomConstraint::AromaticValence(read_aromatic_valence_dsl(de)?.into_ast(&()))
        }
        "multicenter-valence" => {
            AtomConstraint::MulticenterValence(read_multicenter_valence_dsl(de)?.into_ast(&()))
        }
        "donated-pairs" => AtomConstraint::DonatedPairs(read_value_dsl(de)?.into_ast(&())),
        "accepted-pairs" => AtomConstraint::AcceptedPairs(read_value_dsl(de)?.into_ast(&())),
        "degree" => AtomConstraint::Degree(read_value_dsl(de)?.into_ast(&())),
        "total-degree" => AtomConstraint::TotalDegree(read_value_dsl(de)?.into_ast(&())),
        "ring-degree" => AtomConstraint::RingDegree(read_value_dsl(de)?.into_ast(&())),
        "ring-valence" => AtomConstraint::RingValence(read_value_dsl(de)?.into_ast(&())),
        "total-hydrogens" => AtomConstraint::TotalHydrogens(read_value_dsl(de)?.into_ast(&())),
        "ring-membership" => AtomConstraint::RingMembership(read_ring_membership_dsl(de)?),
        "tetrahedral-stereo" => {
            AtomConstraint::TetrahedralStereo(read_tetrahedral_stereo_dsl(de)?)
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["atom-constraint".into()],
            }
            .into());
        }
    };
    consume_single_key_map_close(de, "atom-constraint")?;
    Ok(AtomConstraintDsl(c))
}

pub(super) fn read_bond_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<BondConstraintDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b':' => {
            let name = de.read_keyword_name()?;
            match name.as_ref() {
                "aromatic" => Ok(BondConstraintDsl(BondConstraint::Aromatic)),
                other => Err(DeError::Custom(format!(
                    "unknown bond-constraint keyword :{}",
                    other
                ))
                .into()),
            }
        }
        b'{' => {
            let key = read_single_key_map_header(de)?;
            let c = match key.as_str() {
                "ring-membership" => {
                    BondConstraint::RingMembership(read_ring_membership_dsl(de)?)
                }
                "cis-trans-stereo" => {
                    BondConstraint::CisTransStereo(read_cis_trans_stereo_dsl(de)?)
                }
                other => {
                    return Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["bond-constraint".into()],
                    }
                    .into());
                }
            };
            consume_single_key_map_close(de, "bond-constraint")?;
            Ok(BondConstraintDsl(c))
        }
        b => Err(DeError::TypeMismatch {
            expected: ":aromatic / {:ring-membership …}",
            got: unexpected_byte_kind(b),
            path: vec!["bond-constraint".into()],
        }
        .into()),
    }
}

pub(super) fn read_dative_bond_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DativeBondConstraintDsl, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b':' => {
            let name = de.read_keyword_name()?;
            match name.as_ref() {
                "aromatic" => Ok(DativeBondConstraintDsl::Aromatic),
                other => Err(DeError::Custom(format!(
                    "unknown dative-bond-constraint keyword :{}",
                    other
                ))
                .into()),
            }
        }
        b'{' => {
            let key = read_single_key_map_header(de)?;
            let c = match key.as_str() {
                "ring-membership" => {
                    DativeBondConstraintDsl::RingMembership(read_ring_membership_dsl(de)?)
                }
                other => {
                    return Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["dative-bond-constraint".into()],
                    }
                    .into());
                }
            };
            consume_single_key_map_close(de, "dative-bond-constraint")?;
            Ok(c)
        }
        b => Err(DeError::TypeMismatch {
            expected: ":aromatic / {:ring-membership …}",
            got: unexpected_byte_kind(b),
            path: vec!["dative-bond-constraint".into()],
        }
        .into()),
    }
}

pub(super) fn read_aromatic_system_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystemConstraintDsl, EdnError> {
    let key = read_single_key_map_header(de)?;
    let c = match key.as_str() {
        "electron-count" => {
            AromaticSystemConstraintDsl::ElectronCount(read_value_dsl(de)?.into_ast(&()))
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["aromatic-system-constraint".into()],
            }
            .into());
        }
    };
    consume_single_key_map_close(de, "aromatic-system-constraint")?;
    Ok(c)
}

pub(super) fn read_multicenter_bond_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBondConstraintDsl, EdnError> {
    let key = read_single_key_map_header(de)?;
    let c = match key.as_str() {
        "electron-count" => {
            MulticenterBondConstraintDsl::ElectronCount(read_value_dsl(de)?.into_ast(&()))
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["multicenter-bond-constraint".into()],
            }
            .into());
        }
    };
    consume_single_key_map_close(de, "multicenter-bond-constraint")?;
    Ok(c)
}

pub(super) fn read_noncovalent_bond_constraint_dsl(
    _de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondConstraintDsl, EdnError> {
    Err(DeError::Custom("no value-only noncovalent-bond constraints exist yet".to_string()).into())
}

/// Read a stereo constraint payload (the kind-bearing map) by capturing its value
/// slice and parsing it via `FromEdn` — the kind-aware build lives in the
/// `StereoAtomConstraintDsl::from_edn` impl, so the streaming path bridges rather
/// than reimplementing the map parse incrementally.
/// TODO: FIX THIS TO USE streaming parser
fn read_stereo_atom_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoAtomConstraintDsl, EdnError> {
    let slice = de.read_value_slice()?;
    let edn = umol_edn::read_string(slice)?;
    Ok(StereoAtomConstraintDsl::from_edn(&edn)?)
}

/// TODO: FIX THIS TO USE streaming parser
fn read_stereo_bond_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoBondConstraintDsl, EdnError> {
    let slice = de.read_value_slice()?;
    let edn = umol_edn::read_string(slice)?;
    Ok(StereoBondConstraintDsl::from_edn(&edn)?)
}

fn read_atom_ref_vec(de: &mut EdnStreamDeserializer<'_>) -> Result<Vec<AtomRef>, EdnError> {
    read_vec(de, read_atom_ref)
}

fn read_atom_ref_pair(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<[AtomRef; 2], EdnError> {
    de.consume_byte(b'[')?;
    let a = read_atom_ref(de)?;
    let b = read_atom_ref(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom(format!("{}: expected 2 elements", context)).into());
    }
    Ok([a, b])
}

fn read_atom_constraint_pair(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<[AtomConstraintDsl; 2], EdnError> {
    de.consume_byte(b'[')?;
    let a = read_atom_constraint_dsl(de)?;
    let b = read_atom_constraint_dsl(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom(format!("{}: expected 2 elements", context)).into());
    }
    Ok([a, b])
}

pub(super) fn read_relational_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
    key: &str,
) -> Result<RelationalConstraintDsl, EdnError> {
    use RelationalConstraintDsl as R;
    de.consume_byte(b'[')?;
    let c = match key {
        "dative-bond-donor" => R::DativeBondDonor {
            bond: read_dative_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "dative-bond-acceptor" => R::DativeBondAcceptor {
            bond: read_dative_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "dative-bond-parallels" => R::DativeBondParallels {
            dative: read_dative_bond_ref(de)?,
            parallel: read_bond_ref(de)?,
        },
        "dative-bond-donor-satisfies" => R::DativeBondDonorSatisfies {
            bond: read_dative_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "dative-bond-acceptor-satisfies" => R::DativeBondAcceptorSatisfies {
            bond: read_dative_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "aromatic-system-atoms" => R::AromaticSystemAtoms {
            system: read_aromatic_system_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "aromatic-system-contains" => R::AromaticSystemContains {
            system: read_aromatic_system_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "aromatic-system-contains-all" => R::AromaticSystemContainsAll {
            system: read_aromatic_system_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "aromatic-system-all-atoms" => R::AromaticSystemAllAtoms {
            system: read_aromatic_system_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "aromatic-system-any-atom" => R::AromaticSystemAnyAtom {
            system: read_aromatic_system_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "multicenter-bond-atoms" => R::MulticenterBondAtoms {
            bond: read_multicenter_bond_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "multicenter-bond-contains" => R::MulticenterBondContains {
            bond: read_multicenter_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "multicenter-bond-contains-all" => R::MulticenterBondContainsAll {
            bond: read_multicenter_bond_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "multicenter-bond-all-atoms" => R::MulticenterBondAllAtoms {
            bond: read_multicenter_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "multicenter-bond-any-atom" => R::MulticenterBondAnyAtom {
            bond: read_multicenter_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "noncovalent-bond-ends" => R::NoncovalentBondEnds {
            bond: read_noncovalent_bond_ref(de)?,
            atoms: read_atom_ref_pair(de, "noncovalent-bond-ends")?,
        },
        "noncovalent-bond-contains" => R::NoncovalentBondContains {
            bond: read_noncovalent_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "noncovalent-bond-ends-satisfy" => {
            let bond = read_noncovalent_bond_ref(de)?;
            let [a, b] = read_atom_constraint_pair(de, "noncovalent-bond-ends-satisfy")?;
            R::NoncovalentBondEndsSatisfy {
                bond,
                predicates: [Box::new(a), Box::new(b)],
            }
        }
        "stereo-atom-site" => R::StereoAtomSite {
            stereo_atom: read_stereo_atom_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "stereo-atom-contains" => R::StereoAtomContains {
            stereo_atom: read_stereo_atom_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "stereo-atom-ligands" => R::StereoAtomLigands {
            stereo_atom: read_stereo_atom_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "stereo-atom-all-ligands" => R::StereoAtomAllLigands {
            stereo_atom: read_stereo_atom_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "stereo-atom-any-ligand" => R::StereoAtomAnyLigand {
            stereo_atom: read_stereo_atom_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "stereo-bond-site" => R::StereoBondSite {
            stereo_bond: read_stereo_bond_ref(de)?,
            bond: read_bond_ref(de)?,
        },
        "stereo-bond-contains" => R::StereoBondContains {
            stereo_bond: read_stereo_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "stereo-bond-ligands" => R::StereoBondLigands {
            stereo_bond: read_stereo_bond_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "stereo-bond-all-ligands" => R::StereoBondAllLigands {
            stereo_bond: read_stereo_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "stereo-bond-any-ligand" => R::StereoBondAnyLigand {
            stereo_bond: read_stereo_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        other => {
            unreachable!("read_relational_constraint_dsl called with non-relational key {other}")
        }
    };
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom(format!("{}: expected 2 elements", key)).into());
    }
    Ok(c)
}

pub(super) fn read_sub_pattern_anchor_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<SubPatternAnchorDsl, EdnError> {
    let mut out = SubPatternAnchorDsl::default();
    read_map(de, |d, key| {
        match key {
            "atoms" => out.atoms = read_vec(d, |d| read_ref_pair(d, read_atom_ref, read_atom_ref))?,
            "bonds" => out.bonds = read_vec(d, |d| read_ref_pair(d, read_bond_ref, read_bond_ref))?,
            "dative-bonds" => {
                out.dative_bonds = read_vec(d, |d| {
                    read_ref_pair(d, read_dative_bond_ref, read_dative_bond_ref)
                })?
            }
            "aromatic-systems" => {
                out.aromatic_systems = read_vec(d, |d| {
                    read_ref_pair(d, read_aromatic_system_ref, read_aromatic_system_ref)
                })?
            }
            "multicenter-bonds" => {
                out.multicenter_bonds = read_vec(d, |d| {
                    read_ref_pair(d, read_multicenter_bond_ref, read_multicenter_bond_ref)
                })?
            }
            "noncovalent-bonds" => {
                out.noncovalent_bonds = read_vec(d, |d| {
                    read_ref_pair(d, read_noncovalent_bond_ref, read_noncovalent_bond_ref)
                })?
            }
            "stereo-atoms" => {
                out.stereo_atoms = read_vec(d, |d| {
                    read_ref_pair(d, read_stereo_atom_ref, read_stereo_atom_ref)
                })?
            }
            "stereo-bonds" => {
                out.stereo_bonds = read_vec(d, |d| {
                    read_ref_pair(d, read_stereo_bond_ref, read_stereo_bond_ref)
                })?
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["sub-pattern-anchor".into()],
                }
                .into());
            }
        }
        Ok(())
    })?;
    Ok(out)
}

fn read_ref_pair<A, B>(
    de: &mut EdnStreamDeserializer<'_>,
    read_a: fn(&mut EdnStreamDeserializer<'_>) -> Result<A, EdnError>,
    read_b: fn(&mut EdnStreamDeserializer<'_>) -> Result<B, EdnError>,
) -> Result<(A, B), EdnError> {
    de.consume_byte(b'[')?;
    let a = read_a(de)?;
    let b = read_b(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom("anchor pair must have 2 elements".into()).into());
    }
    Ok((a, b))
}

pub(super) fn read_molecule_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
    key: &str,
) -> Result<MoleculeConstraintDsl, EdnError> {
    // Caller has already consumed the outer `{` and the dispatch key.
    let c = match key {
        "charge-sum" => {
            let mut atoms = None;
            let mut sum = None;
            read_map(de, |d, k| {
                match k {
                    "atoms" => atoms = Some(read_vec(d, read_atom_ref)?),
                    "sum" => sum = Some(read_value_dsl(d)?),
                    _ => d.read_skip_value()?,
                }
                Ok(())
            })?;
            MoleculeConstraintDsl::ChargeSum {
                atoms,
                sum: sum.ok_or_else(|| missing("sum", "charge-sum"))?,
            }
        }
        "spin-sum" => {
            let mut atoms = None;
            let mut spin = None;
            read_map(de, |d, k| {
                match k {
                    "atoms" => atoms = Some(read_vec(d, read_atom_ref)?),
                    "spin" => spin = Some(read_spin_state(d)?),
                    _ => d.read_skip_value()?,
                }
                Ok(())
            })?;
            MoleculeConstraintDsl::SpinSum {
                atoms,
                spin: spin.ok_or_else(|| missing("spin", "spin-sum"))?,
            }
        }
        "bond-order-sum" => {
            let mut bonds = None;
            let mut sum = None;
            read_map(de, |d, k| {
                match k {
                    "bonds" => bonds = Some(read_vec(d, read_bond_ref)?),
                    "sum" => sum = Some(read_value_dsl(d)?),
                    _ => d.read_skip_value()?,
                }
                Ok(())
            })?;
            MoleculeConstraintDsl::BondOrderSum {
                bonds,
                sum: sum.ok_or_else(|| missing("sum", "bond-order-sum"))?,
            }
        }
        "connected" => {
            let mut atoms = None;
            read_map(de, |d, k| {
                match k {
                    "atoms" => atoms = Some(read_vec(d, read_atom_ref)?),
                    _ => d.read_skip_value()?,
                }
                Ok(())
            })?;
            MoleculeConstraintDsl::Connected { atoms }
        }
        "sub-pattern" => {
            let mut anchor = None;
            let mut pattern = None;
            read_map(de, |d, k| {
                match k {
                    "anchor" => anchor = Some(read_sub_pattern_anchor_dsl(d)?),
                    "pattern" => {
                        let input = super::molecule::read_molecule_input(d)?;
                        let (ast, _metadata) = input
                            .into_ast()
                            .map_err(|e| DeError::Custom(e.to_string()))?;
                        pattern = Some(Box::new(ast));
                    }
                    _ => d.read_skip_value()?,
                }
                Ok(())
            })?;
            MoleculeConstraintDsl::SubPattern {
                anchor: anchor.ok_or_else(|| missing("anchor", "sub-pattern"))?,
                pattern: pattern.ok_or_else(|| missing("pattern", "sub-pattern"))?,
            }
        }
        other => unreachable!("read_molecule_constraint_dsl called with non-molecule key {other}"),
    };
    Ok(c)
}

pub(super) fn read_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<ConstraintDsl, EdnError> {
    let key = read_single_key_map_header(de)?;
    let c = match key.as_str() {
        "atom" => {
            let (r, inner) = read_entity_leaf(de, read_atom_ref, read_atom_constraint_dsl, "atom")?;
            ConstraintDsl::Atom(r, inner)
        }
        "bond" => {
            let (r, inner) = read_entity_leaf(de, read_bond_ref, read_bond_constraint_dsl, "bond")?;
            ConstraintDsl::Bond(r, inner)
        }
        "dative-bond" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_dative_bond_ref,
                read_dative_bond_constraint_dsl,
                "dative-bond",
            )?;
            ConstraintDsl::DativeBond(r, inner)
        }
        "aromatic-system" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_aromatic_system_ref,
                read_aromatic_system_constraint_dsl,
                "aromatic-system",
            )?;
            ConstraintDsl::AromaticSystem(r, inner)
        }
        "multicenter-bond" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_multicenter_bond_ref,
                read_multicenter_bond_constraint_dsl,
                "multicenter-bond",
            )?;
            ConstraintDsl::MulticenterBond(r, inner)
        }
        "noncovalent-bond" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_noncovalent_bond_ref,
                read_noncovalent_bond_constraint_dsl,
                "noncovalent-bond",
            )?;
            ConstraintDsl::NoncovalentBond(r, inner)
        }
        "stereo-atom" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_stereo_atom_ref,
                read_stereo_atom_constraint_dsl,
                "stereo-atom",
            )?;
            ConstraintDsl::StereoAtom(r, inner)
        }
        "stereo-bond" => {
            let (r, inner) = read_entity_leaf(
                de,
                read_stereo_bond_ref,
                read_stereo_bond_constraint_dsl,
                "stereo-bond",
            )?;
            ConstraintDsl::StereoBond(r, inner)
        }
        "and" => ConstraintDsl::And(read_vec(de, read_constraint_dsl)?),
        "or" => ConstraintDsl::Or(read_vec(de, read_constraint_dsl)?),
        "not" => ConstraintDsl::Not(Box::new(read_constraint_dsl(de)?)),
        "charge-sum" | "spin-sum" | "bond-order-sum" | "connected" | "sub-pattern" => {
            ConstraintDsl::Molecule(read_molecule_constraint_dsl(de, key.as_str())?)
        }
        k if RELATIONAL_KEYS.contains(&k) => {
            ConstraintDsl::Relational(read_relational_constraint_dsl(de, k)?)
        }
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["constraint".into()],
            }
            .into());
        }
    };
    consume_single_key_map_close(de, "constraint")?;
    Ok(c)
}

fn read_entity_leaf<R, C>(
    de: &mut EdnStreamDeserializer<'_>,
    read_ref: fn(&mut EdnStreamDeserializer<'_>) -> Result<R, EdnError>,
    read_inner: fn(&mut EdnStreamDeserializer<'_>) -> Result<C, EdnError>,
    context: &'static str,
) -> Result<(R, C), EdnError> {
    de.consume_byte(b'[')?;
    let r = read_ref(de)?;
    let c = read_inner(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(
            DeError::Custom(format!("{} entity leaf must have 2 elements", context)).into(),
        );
    }
    Ok((r, c))
}

pub(super) fn read_constraints_dsl(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<Vec<ConstraintDsl>, EdnError> {
    read_vec(de, read_constraint_dsl)
}

/// Surface DSL wrapper around `MoleculeConstraint`. EDN form is a single-key
/// map keyed by the variant: `{:charge-sum {...}}`, `{:spin-sum {...}}`,
/// `{:bond-order-sum {...}}`, `{:connected {...}}`, `{:sub-pattern {...}}`.
///
/// For `ChargeSum` / `SpinSum` / `BondOrderSum` / `Connected`, the
/// `atoms` (or `bonds`) field is `None` to denote the entire molecule's
/// atoms (or bonds). Empty subset must be expressed explicitly as
/// `Some(vec![])`.
///
/// `SubPattern` carries a `Box<MoleculeAst>` directly: defaults are a
/// ground-input convenience that has no meaning for patterns, where
/// `Undetermined` is a wildcard. The AST↔DSL bridge for the pattern is
/// the identity; the EDN bridge wraps with empty `Metadata` so refs render
/// as numeric indices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraintDsl {
    ChargeSum {
        atoms: Option<Vec<AtomRef>>,
        sum: ValueDsl,
    },
    SpinSum {
        atoms: Option<Vec<AtomRef>>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Option<Vec<BondRef>>,
        sum: ValueDsl,
    },
    Connected {
        atoms: Option<Vec<AtomRef>>,
    },
    SubPattern {
        anchor: SubPatternAnchorDsl,
        pattern: Box<MoleculeAst>,
    },
}

fn atom_subset_from_ast(atoms: &Option<Vec<AtomId>>, meta: &Metadata) -> Option<Vec<AtomRef>> {
    atoms
        .as_ref()
        .map(|v| v.iter().map(|&a| AtomRef::from_ast(a, meta)).collect())
}

fn bond_subset_from_ast(bonds: &Option<Vec<BondId>>, meta: &Metadata) -> Option<Vec<BondRef>> {
    bonds
        .as_ref()
        .map(|v| v.iter().map(|&b| BondRef::from_ast(b, meta)).collect())
}

fn atom_subset_into_ast(
    atoms: Option<Vec<AtomRef>>,
    count: usize,
    meta: &Metadata,
) -> Result<Option<Vec<AtomId>>, ParseError> {
    atoms
        .map(|v| {
            v.into_iter()
                .map(|r| r.into_ast(count, meta))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
}

fn bond_subset_into_ast(
    bonds: Option<Vec<BondRef>>,
    count: usize,
    meta: &Metadata,
) -> Result<Option<Vec<BondId>>, ParseError> {
    bonds
        .map(|v| {
            v.into_iter()
                .map(|r| r.into_ast(count, meta))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
}

impl MoleculeConstraintDsl {
    pub(crate) fn from_ast(c: &MoleculeConstraint, meta: &Metadata) -> Result<Self, ParseError> {
        Ok(match c {
            MoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum {
                atoms: atom_subset_from_ast(atoms, meta),
                sum: ValueDsl::from_ast(sum, &()),
            },
            MoleculeConstraint::SpinSum { atoms, spin } => Self::SpinSum {
                atoms: atom_subset_from_ast(atoms, meta),
                spin: spin.clone(),
            },
            MoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum {
                bonds: bond_subset_from_ast(bonds, meta),
                sum: ValueDsl::from_ast(sum, &()),
            },
            MoleculeConstraint::Connected { atoms } => Self::Connected {
                atoms: atom_subset_from_ast(atoms, meta),
            },
            MoleculeConstraint::SubPattern { anchor, pattern } => {
                let pattern_meta = Metadata::default();
                let anchor_dsl = SubPatternAnchorDsl::from_ast_pair(anchor, meta, &pattern_meta);
                Self::SubPattern {
                    anchor: anchor_dsl,
                    pattern: pattern.clone(),
                }
            }
        })
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<MoleculeConstraint, ParseError> {
        Ok(match self {
            Self::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms: atom_subset_into_ast(atoms, counts.atom_count, meta)?,
                sum: sum.into_ast(&()),
            },
            Self::SpinSum { atoms, spin } => MoleculeConstraint::SpinSum {
                atoms: atom_subset_into_ast(atoms, counts.atom_count, meta)?,
                spin,
            },
            Self::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds: bond_subset_into_ast(bonds, counts.bond_count, meta)?,
                sum: sum.into_ast(&()),
            },
            Self::Connected { atoms } => MoleculeConstraint::Connected {
                atoms: atom_subset_into_ast(atoms, counts.atom_count, meta)?,
            },
            Self::SubPattern { anchor, pattern } => {
                let pattern_counts = EntityCounts::from_ast(&pattern);
                let pattern_meta = Metadata::default();
                let anchor_ast =
                    anchor.into_ast_pair(counts, meta, &pattern_counts, &pattern_meta)?;
                MoleculeConstraint::SubPattern {
                    anchor: anchor_ast,
                    pattern,
                }
            }
        })
    }
}

impl<'de> FromEdn<'de> for MoleculeConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "molecule-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "molecule-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["molecule-constraint".into()],
            });
        };
        Ok(match key.name() {
            "charge-sum" => {
                let (atoms, sum) = parse_sum_map::<AtomRef>(v, "charge-sum", "atoms")?;
                Self::ChargeSum { atoms, sum }
            }
            "spin-sum" => {
                let m = expect_map(v, "spin-sum")?;
                let atoms = parse_optional_refs::<AtomRef>(m, "atoms")?;
                let spin_edn = m.get_keyword("spin").ok_or_else(|| DeError::MissingField {
                    key: "spin".into(),
                    path: vec!["spin-sum".into()],
                })?;
                Self::SpinSum {
                    atoms,
                    spin: parse_spin(spin_edn)?,
                }
            }
            "bond-order-sum" => {
                let (bonds, sum) = parse_sum_map::<BondRef>(v, "bond-order-sum", "bonds")?;
                Self::BondOrderSum { bonds, sum }
            }
            "connected" => {
                let m = expect_map(v, "connected")?;
                Self::Connected {
                    atoms: parse_optional_refs::<AtomRef>(m, "atoms")?,
                }
            }
            "sub-pattern" => {
                let m = expect_map(v, "sub-pattern")?;
                let anchor_edn = m
                    .get_keyword("anchor")
                    .ok_or_else(|| DeError::MissingField {
                        key: "anchor".into(),
                        path: vec!["sub-pattern".into()],
                    })?;
                let pattern_edn =
                    m.get_keyword("pattern")
                        .ok_or_else(|| DeError::MissingField {
                            key: "pattern".into(),
                            path: vec!["sub-pattern".into()],
                        })?;
                let pattern_dsl = MoleculeDsl::from_edn(pattern_edn)?;
                let (pattern_ast, _) = pattern_dsl.into_parts();
                Self::SubPattern {
                    anchor: SubPatternAnchorDsl::from_edn(anchor_edn)?,
                    pattern: Box::new(pattern_ast),
                }
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for MoleculeConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::ChargeSum { atoms, sum } => ("charge-sum", render_sum_map("atoms", atoms, sum)),
            Self::SpinSum { atoms, spin } => {
                let mut m = EdnMap::with_capacity(2);
                if let Some(refs) = atoms {
                    m.insert(Edn::keyword("atoms"), render_refs(refs));
                }
                m.insert(Edn::keyword("spin"), render_spin(spin));
                ("spin-sum", Edn::Map(m))
            }
            Self::BondOrderSum { bonds, sum } => {
                ("bond-order-sum", render_sum_map("bonds", bonds, sum))
            }
            Self::Connected { atoms } => {
                let mut m = EdnMap::with_capacity(1);
                if let Some(refs) = atoms {
                    m.insert(Edn::keyword("atoms"), render_refs(refs));
                }
                ("connected", Edn::Map(m))
            }
            Self::SubPattern { anchor, pattern } => {
                let pattern_dsl = MoleculeDsl::from_parts((**pattern).clone(), Metadata::default());
                let mut m = EdnMap::with_capacity(2);
                m.insert(Edn::keyword("anchor"), anchor.to_edn());
                m.insert(Edn::keyword("pattern"), pattern_dsl.to_edn());
                ("sub-pattern", Edn::Map(m))
            }
        };
        let mut outer = EdnMap::with_capacity(1);
        outer.insert(Edn::Keyword(EdnKeyword::owned(key.into())), value);
        Edn::Map(outer)
    }
}

/// Surface DSL wrapper around `SubPatternAnchor`. Each vector carries
/// `(target, pattern)` ref pairs. Target-side refs resolve against the outer
/// molecule's `Metadata`; pattern-side refs resolve against the pattern
/// molecule's `Metadata`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SubPatternAnchorDsl {
    pub atoms: Vec<(AtomRef, AtomRef)>,
    pub bonds: Vec<(BondRef, BondRef)>,
    pub dative_bonds: Vec<(DativeBondRef, DativeBondRef)>,
    pub aromatic_systems: Vec<(AromaticSystemRef, AromaticSystemRef)>,
    pub multicenter_bonds: Vec<(MulticenterBondRef, MulticenterBondRef)>,
    pub noncovalent_bonds: Vec<(NoncovalentBondRef, NoncovalentBondRef)>,
    pub stereo_atoms: Vec<(StereoAtomRef, StereoAtomRef)>,
    pub stereo_bonds: Vec<(StereoBondRef, StereoBondRef)>,
}

impl SubPatternAnchorDsl {
    /// Build from an AST anchor. `target_meta` is the outer molecule's
    /// metadata; `pattern_meta` is the pattern molecule's metadata.
    pub fn from_ast_pair(
        anchor: &SubPatternAnchor,
        target_meta: &Metadata,
        pattern_meta: &Metadata,
    ) -> Self {
        Self {
            atoms: anchor
                .atoms()
                .iter()
                .map(|&(t, p)| {
                    (
                        AtomRef::from_ast(t, target_meta),
                        AtomRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            bonds: anchor
                .bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        BondRef::from_ast(t, target_meta),
                        BondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            dative_bonds: anchor
                .dative_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        DativeBondRef::from_ast(t, target_meta),
                        DativeBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            aromatic_systems: anchor
                .aromatic_systems()
                .iter()
                .map(|&(t, p)| {
                    (
                        AromaticSystemRef::from_ast(t, target_meta),
                        AromaticSystemRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            multicenter_bonds: anchor
                .multicenter_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        MulticenterBondRef::from_ast(t, target_meta),
                        MulticenterBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            noncovalent_bonds: anchor
                .noncovalent_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        NoncovalentBondRef::from_ast(t, target_meta),
                        NoncovalentBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            stereo_atoms: anchor
                .stereo_atoms()
                .iter()
                .map(|&(t, p)| {
                    (
                        StereoAtomRef::from_ast(t, target_meta),
                        StereoAtomRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
            stereo_bonds: anchor
                .stereo_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        StereoBondRef::from_ast(t, target_meta),
                        StereoBondRef::from_ast(p, pattern_meta),
                    )
                })
                .collect(),
        }
    }

    /// Resolve to an AST anchor. `target_*` carry outer-molecule counts +
    /// metadata; `pattern_*` carry pattern-molecule counts + metadata.
    pub(crate) fn into_ast_pair(
        self,
        target_counts: &EntityCounts,
        target_meta: &Metadata,
        pattern_counts: &EntityCounts,
        pattern_meta: &Metadata,
    ) -> Result<SubPatternAnchor, ParseError> {
        let mut anchor = SubPatternAnchor::new();
        for (t, p) in self.atoms {
            anchor.push_atom(
                t.into_ast(target_counts.atom_count, target_meta)?,
                p.into_ast(pattern_counts.atom_count, pattern_meta)?,
            );
        }
        for (t, p) in self.bonds {
            anchor.push_bond(
                t.into_ast(target_counts.bond_count, target_meta)?,
                p.into_ast(pattern_counts.bond_count, pattern_meta)?,
            );
        }
        for (t, p) in self.dative_bonds {
            anchor.push_dative_bond(
                t.into_ast(target_counts.dative_bond_count, target_meta)?,
                p.into_ast(pattern_counts.dative_bond_count, pattern_meta)?,
            );
        }
        for (t, p) in self.aromatic_systems {
            anchor.push_aromatic_system(
                t.into_ast(target_counts.aromatic_system_count, target_meta)?,
                p.into_ast(pattern_counts.aromatic_system_count, pattern_meta)?,
            );
        }
        for (t, p) in self.multicenter_bonds {
            anchor.push_multicenter_bond(
                t.into_ast(target_counts.multicenter_bond_count, target_meta)?,
                p.into_ast(pattern_counts.multicenter_bond_count, pattern_meta)?,
            );
        }
        for (t, p) in self.noncovalent_bonds {
            anchor.push_noncovalent_bond(
                t.into_ast(target_counts.noncovalent_bond_count, target_meta)?,
                p.into_ast(pattern_counts.noncovalent_bond_count, pattern_meta)?,
            );
        }
        for (t, p) in self.stereo_atoms {
            anchor.push_stereo_atom(
                t.into_ast(target_counts.stereo_atom_count, target_meta)?,
                p.into_ast(pattern_counts.stereo_atom_count, pattern_meta)?,
            );
        }
        for (t, p) in self.stereo_bonds {
            anchor.push_stereo_bond(
                t.into_ast(target_counts.stereo_bond_count, target_meta)?,
                p.into_ast(pattern_counts.stereo_bond_count, pattern_meta)?,
            );
        }
        Ok(anchor)
    }
}

impl<'de> FromEdn<'de> for SubPatternAnchorDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let m = expect_map(edn, "sub-pattern-anchor")?;
        let mut out = Self::default();
        for (k, v) in m.iter() {
            let Edn::Keyword(key) = k else {
                return Err(DeError::TypeMismatch {
                    expected: "keyword key",
                    got: k.kind(),
                    path: vec!["sub-pattern-anchor".into()],
                });
            };
            match key.name() {
                "atoms" => out.atoms = parse_ref_pairs::<AtomRef, AtomRef>(v)?,
                "bonds" => out.bonds = parse_ref_pairs::<BondRef, BondRef>(v)?,
                "dative-bonds" => {
                    out.dative_bonds = parse_ref_pairs::<DativeBondRef, DativeBondRef>(v)?
                }
                "aromatic-systems" => {
                    out.aromatic_systems =
                        parse_ref_pairs::<AromaticSystemRef, AromaticSystemRef>(v)?
                }
                "multicenter-bonds" => {
                    out.multicenter_bonds =
                        parse_ref_pairs::<MulticenterBondRef, MulticenterBondRef>(v)?
                }
                "noncovalent-bonds" => {
                    out.noncovalent_bonds =
                        parse_ref_pairs::<NoncovalentBondRef, NoncovalentBondRef>(v)?
                }
                "stereo-atoms" => {
                    out.stereo_atoms = parse_ref_pairs::<StereoAtomRef, StereoAtomRef>(v)?
                }
                "stereo-bonds" => {
                    out.stereo_bonds = parse_ref_pairs::<StereoBondRef, StereoBondRef>(v)?
                }
                other => {
                    return Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["sub-pattern-anchor".into()],
                    });
                }
            }
        }
        Ok(out)
    }
}

impl ToEdn for SubPatternAnchorDsl {
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(8);
        if !self.atoms.is_empty() {
            m.insert(Edn::keyword("atoms"), render_ref_pairs(&self.atoms));
        }
        if !self.bonds.is_empty() {
            m.insert(Edn::keyword("bonds"), render_ref_pairs(&self.bonds));
        }
        if !self.dative_bonds.is_empty() {
            m.insert(
                Edn::keyword("dative-bonds"),
                render_ref_pairs(&self.dative_bonds),
            );
        }
        if !self.aromatic_systems.is_empty() {
            m.insert(
                Edn::keyword("aromatic-systems"),
                render_ref_pairs(&self.aromatic_systems),
            );
        }
        if !self.multicenter_bonds.is_empty() {
            m.insert(
                Edn::keyword("multicenter-bonds"),
                render_ref_pairs(&self.multicenter_bonds),
            );
        }
        if !self.noncovalent_bonds.is_empty() {
            m.insert(
                Edn::keyword("noncovalent-bonds"),
                render_ref_pairs(&self.noncovalent_bonds),
            );
        }
        if !self.stereo_atoms.is_empty() {
            m.insert(
                Edn::keyword("stereo-atoms"),
                render_ref_pairs(&self.stereo_atoms),
            );
        }
        if !self.stereo_bonds.is_empty() {
            m.insert(
                Edn::keyword("stereo-bonds"),
                render_ref_pairs(&self.stereo_bonds),
            );
        }
        Edn::Map(m)
    }
}

fn expect_map<'e>(edn: &'e Edn<'e>, context: &'static str) -> Result<&'e EdnMap<'e>, DeError> {
    match edn {
        Edn::Map(m) => Ok(m),
        other => Err(DeError::TypeMismatch {
            expected: "map",
            got: other.kind(),
            path: vec![context.into()],
        }),
    }
}

fn parse_refs<R>(edn: &Edn<'_>) -> Result<Vec<R>, DeError>
where
    R: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of refs",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    v.iter().map(R::from_edn).collect()
}

fn render_refs<R: ToEdn>(refs: &[R]) -> Edn<'static> {
    Edn::Vector(refs.iter().map(R::to_edn).collect::<Vec<_>>().into())
}

fn parse_ref_pairs<A, B>(edn: &Edn<'_>) -> Result<Vec<(A, B)>, DeError>
where
    A: for<'de> FromEdn<'de>,
    B: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of [target pattern] pairs",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    v.iter()
        .map(|e| {
            let Edn::Vector(pair) = e else {
                return Err(DeError::TypeMismatch {
                    expected: "2-element vector [target pattern]",
                    got: e.kind(),
                    path: Vec::new(),
                });
            };
            if pair.len() != 2 {
                return Err(DeError::Custom(format!(
                    "anchor pair must have 2 elements, got {}",
                    pair.len()
                )));
            }
            Ok((A::from_edn(&pair[0])?, B::from_edn(&pair[1])?))
        })
        .collect()
}

fn render_ref_pairs<A: ToEdn, B: ToEdn>(pairs: &[(A, B)]) -> Edn<'static> {
    Edn::Vector(
        pairs
            .iter()
            .map(|(a, b)| Edn::Vector(vec![a.to_edn(), b.to_edn()].into()))
            .collect::<Vec<_>>()
            .into(),
    )
}

fn parse_optional_refs<R>(m: &EdnMap<'_>, refs_key: &'static str) -> Result<Option<Vec<R>>, DeError>
where
    R: for<'de> FromEdn<'de>,
{
    match m.get_keyword(refs_key) {
        Some(refs_edn) => parse_refs::<R>(refs_edn).map(Some),
        None => Ok(None),
    }
}

fn parse_sum_map<R>(
    edn: &Edn<'_>,
    context: &'static str,
    refs_key: &'static str,
) -> Result<(Option<Vec<R>>, ValueDsl), DeError>
where
    R: for<'de> FromEdn<'de>,
{
    let m = expect_map(edn, context)?;
    let refs = parse_optional_refs::<R>(m, refs_key)?;
    let sum_edn = m.get_keyword("sum").ok_or_else(|| DeError::MissingField {
        key: "sum".into(),
        path: vec![context.into()],
    })?;
    Ok((refs, ValueDsl::from_edn(sum_edn)?))
}

fn render_sum_map<R: ToEdn>(refs_key: &str, refs: &Option<Vec<R>>, sum: &ValueDsl) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    if let Some(v) = refs {
        m.insert(
            Edn::Keyword(EdnKeyword::owned(refs_key.into())),
            render_refs(v),
        );
    }
    m.insert(Edn::keyword("sum"), sum.to_edn());
    Edn::Map(m)
}

fn parse_spin(edn: &Edn<'_>) -> Result<SpinStateAst, DeError> {
    let m = expect_map(edn, "spin")?;
    let unpaired = m
        .get_keyword("unpaired")
        .ok_or_else(|| DeError::MissingField {
            key: "unpaired".into(),
            path: vec!["spin".into()],
        })?;
    let multiplicity = m
        .get_keyword("multiplicity")
        .ok_or_else(|| DeError::MissingField {
            key: "multiplicity".into(),
            path: vec!["spin".into()],
        })?;
    Ok(SpinStateAst {
        unpaired: ValueDsl::from_edn(unpaired)?.into_ast(&()),
        multiplicity: ValueDsl::from_edn(multiplicity)?.into_ast(&()),
    })
}

fn render_spin(spin: &SpinStateAst) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    m.insert(
        Edn::keyword("unpaired"),
        ValueDsl::from_ast(&spin.unpaired, &()).to_edn(),
    );
    m.insert(
        Edn::keyword("multiplicity"),
        ValueDsl::from_ast(&spin.multiplicity, &()).to_edn(),
    );
    Edn::Map(m)
}

/// Surface DSL wrapper around `Constraint`. Single-key-map EDN form:
/// `{:atom [<ref> <atom-constraint>]}` for narrow entity leaves;
/// `{:dative-bond-donor [<bond_ref> <atom_ref>]}` etc. for relational
/// (cross-entity) leaves; `{:charge-sum {...}}` etc. for molecule-scope
/// leaves (keys flattened from `MoleculeConstraintDsl`); `{:and [...]}` /
/// `{:or [...]}` / `{:not <c>}` for combinators.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstraintDsl {
    Atom(AtomRef, AtomConstraintDsl),
    Bond(BondRef, BondConstraintDsl),
    DativeBond(DativeBondRef, DativeBondConstraintDsl),
    AromaticSystem(AromaticSystemRef, AromaticSystemConstraintDsl),
    MulticenterBond(MulticenterBondRef, MulticenterBondConstraintDsl),
    NoncovalentBond(NoncovalentBondRef, NoncovalentBondConstraintDsl),
    StereoAtom(StereoAtomRef, StereoAtomConstraintDsl),
    StereoBond(StereoBondRef, StereoBondConstraintDsl),
    Relational(RelationalConstraintDsl),
    Molecule(MoleculeConstraintDsl),
    And(Vec<ConstraintDsl>),
    Or(Vec<ConstraintDsl>),
    Not(Box<ConstraintDsl>),
}

impl<'de> FromEdn<'de> for ConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let m = expect_map(edn, "constraint")?;
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["constraint".into()],
            });
        };
        Ok(match key.name() {
            "atom" => {
                let (r, c) = parse_entity_leaf::<AtomRef, AtomConstraintDsl>(v, "atom")?;
                Self::Atom(r, c)
            }
            "bond" => {
                let (r, c) = parse_entity_leaf::<BondRef, BondConstraintDsl>(v, "bond")?;
                Self::Bond(r, c)
            }
            "dative-bond" => {
                let (r, c) =
                    parse_entity_leaf::<DativeBondRef, DativeBondConstraintDsl>(v, "dative-bond")?;
                Self::DativeBond(r, c)
            }
            "aromatic-system" => {
                let (r, c) = parse_entity_leaf::<AromaticSystemRef, AromaticSystemConstraintDsl>(
                    v,
                    "aromatic-system",
                )?;
                Self::AromaticSystem(r, c)
            }
            "multicenter-bond" => {
                let (r, c) = parse_entity_leaf::<MulticenterBondRef, MulticenterBondConstraintDsl>(
                    v,
                    "multicenter-bond",
                )?;
                Self::MulticenterBond(r, c)
            }
            "noncovalent-bond" => {
                let (r, c) = parse_entity_leaf::<NoncovalentBondRef, NoncovalentBondConstraintDsl>(
                    v,
                    "noncovalent-bond",
                )?;
                Self::NoncovalentBond(r, c)
            }
            "stereo-atom" => {
                let (r, c) =
                    parse_entity_leaf::<StereoAtomRef, StereoAtomConstraintDsl>(v, "stereo-atom")?;
                Self::StereoAtom(r, c)
            }
            "stereo-bond" => {
                let (r, c) =
                    parse_entity_leaf::<StereoBondRef, StereoBondConstraintDsl>(v, "stereo-bond")?;
                Self::StereoBond(r, c)
            }
            "and" => Self::And(parse_constraint_vec(v, "and")?),
            "or" => Self::Or(parse_constraint_vec(v, "or")?),
            "not" => Self::Not(Box::new(ConstraintDsl::from_edn(v)?)),
            // Molecule-scope keys: delegate to MoleculeConstraintDsl.
            "charge-sum" | "spin-sum" | "bond-order-sum" | "connected" | "sub-pattern" => {
                Self::Molecule(MoleculeConstraintDsl::from_edn(edn)?)
            }
            k if RELATIONAL_KEYS.contains(&k) => {
                Self::Relational(RelationalConstraintDsl::from_edn(edn)?)
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for ConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Atom(r, c) => entity_leaf_edn("atom", r, c),
            Self::Bond(r, c) => entity_leaf_edn("bond", r, c),
            Self::DativeBond(r, c) => entity_leaf_edn("dative-bond", r, c),
            Self::AromaticSystem(r, c) => entity_leaf_edn("aromatic-system", r, c),
            Self::MulticenterBond(r, c) => entity_leaf_edn("multicenter-bond", r, c),
            Self::NoncovalentBond(_, c) => match *c {},
            Self::StereoAtom(r, c) => entity_leaf_edn("stereo-atom", r, c),
            Self::StereoBond(r, c) => entity_leaf_edn("stereo-bond", r, c),
            Self::Relational(r) => r.to_edn(),
            Self::Molecule(m) => m.to_edn(),
            Self::And(xs) => combinator_edn("and", xs),
            Self::Or(xs) => combinator_edn("or", xs),
            Self::Not(c) => {
                let mut m = EdnMap::with_capacity(1);
                m.insert(Edn::keyword("not"), c.to_edn());
                Edn::Map(m)
            }
        }
    }
}

impl ConstraintDsl {
    pub(crate) fn from_ast(c: &Constraint, meta: &Metadata) -> Result<Self, ParseError> {
        Ok(match c {
            Constraint::Atom(id, c) => Self::Atom(
                AtomRef::from_ast(*id, meta),
                AtomConstraintDsl::from_ast(c, &()),
            ),
            Constraint::Bond(id, c) => Self::Bond(
                BondRef::from_ast(*id, meta),
                BondConstraintDsl::from_ast(c, &()),
            ),
            Constraint::DativeBond(id, c) => Self::DativeBond(
                DativeBondRef::from_ast(*id, meta),
                DativeBondConstraintDsl::from_ast(c),
            ),
            Constraint::AromaticSystem(id, c) => Self::AromaticSystem(
                AromaticSystemRef::from_ast(*id, meta),
                AromaticSystemConstraintDsl::from_ast(c),
            ),
            Constraint::MulticenterBond(id, c) => Self::MulticenterBond(
                MulticenterBondRef::from_ast(*id, meta),
                MulticenterBondConstraintDsl::from_ast(c),
            ),
            Constraint::NoncovalentBond(_, c) => match *c {},
            Constraint::StereoAtom(id, kind, c) => Self::StereoAtom(
                StereoAtomRef::from_ast(*id, meta),
                StereoAtomConstraintDsl(*kind, c.clone()),
            ),
            Constraint::StereoBond(id, kind, c) => Self::StereoBond(
                StereoBondRef::from_ast(*id, meta),
                StereoBondConstraintDsl(*kind, c.clone()),
            ),
            Constraint::Relational(rel) => {
                Self::Relational(RelationalConstraintDsl::from_ast(rel, meta))
            }
            Constraint::Molecule(m) => Self::Molecule(MoleculeConstraintDsl::from_ast(m, meta)?),
            Constraint::And(xs) => Self::And(
                xs.iter()
                    .map(|c| ConstraintDsl::from_ast(c, meta))
                    .collect::<Result<_, _>>()?,
            ),
            Constraint::Or(xs) => Self::Or(
                xs.iter()
                    .map(|c| ConstraintDsl::from_ast(c, meta))
                    .collect::<Result<_, _>>()?,
            ),
            Constraint::Not(c) => Self::Not(Box::new(ConstraintDsl::from_ast(c, meta)?)),
        })
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<Constraint, ParseError> {
        Ok(match self {
            Self::Atom(r, c) => {
                Constraint::Atom(r.into_ast(counts.atom_count, meta)?, c.into_ast(&()))
            }
            Self::Bond(r, c) => {
                Constraint::Bond(r.into_ast(counts.bond_count, meta)?, c.into_ast(&()))
            }
            Self::DativeBond(r, c) => {
                Constraint::DativeBond(r.into_ast(counts.dative_bond_count, meta)?, c.into_ast())
            }
            Self::AromaticSystem(r, c) => Constraint::AromaticSystem(
                r.into_ast(counts.aromatic_system_count, meta)?,
                c.into_ast(),
            ),
            Self::MulticenterBond(r, c) => Constraint::MulticenterBond(
                r.into_ast(counts.multicenter_bond_count, meta)?,
                c.into_ast(),
            ),
            Self::NoncovalentBond(_, c) => match c {},
            Self::StereoAtom(r, StereoAtomConstraintDsl(kind, c)) => {
                Constraint::StereoAtom(r.into_ast(counts.stereo_atom_count, meta)?, kind, c)
            }
            Self::StereoBond(r, StereoBondConstraintDsl(kind, c)) => {
                Constraint::StereoBond(r.into_ast(counts.stereo_bond_count, meta)?, kind, c)
            }
            Self::Relational(r) => Constraint::Relational(r.into_ast(counts, meta)?),
            Self::Molecule(m) => Constraint::Molecule(m.into_ast(counts, meta)?),
            Self::And(xs) => Constraint::And(
                xs.into_iter()
                    .map(|c| c.into_ast(counts, meta))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Or(xs) => Constraint::Or(
                xs.into_iter()
                    .map(|c| c.into_ast(counts, meta))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Not(c) => Constraint::Not(Box::new(c.into_ast(counts, meta)?)),
        })
    }
}

/// Surface DSL wrapper around `Constraints` (a flat vec of `Constraint`).
/// EDN form: a vector of `ConstraintDsl` EDN forms.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConstraintsDsl(pub Vec<ConstraintDsl>);

impl<'de> FromEdn<'de> for ConstraintsDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Vector(v) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "vector of constraints",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        Ok(Self(
            v.iter()
                .map(ConstraintDsl::from_edn)
                .collect::<Result<_, _>>()?,
        ))
    }
}

impl ToEdn for ConstraintsDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Vector(self.0.iter().map(|c| c.to_edn()).collect::<Vec<_>>().into())
    }
}

impl ConstraintsDsl {
    /// Vacuous constraints (per `Constraint::is_vacuous`) are dropped
    /// during the AST → DSL lowering, matching the canonical-rendering
    /// rule: a constraint that asserts nothing does not appear in the
    /// canonical surface form.
    pub(crate) fn from_ast(cs: &Constraints, meta: &Metadata) -> Result<Self, ParseError> {
        Ok(Self(
            cs.iter()
                .filter(|c| !c.is_vacuous())
                .map(|c| ConstraintDsl::from_ast(c, meta))
                .collect::<Result<_, _>>()?,
        ))
    }

    pub(crate) fn into_ast(
        self,
        counts: &EntityCounts,
        meta: &Metadata,
    ) -> Result<Constraints, ParseError> {
        let mut out = Constraints::new();
        for c in self.0 {
            out.push(c.into_ast(counts, meta)?);
        }
        Ok(out)
    }
}

fn parse_entity_leaf<R, C>(edn: &Edn<'_>, context: &'static str) -> Result<(R, C), DeError>
where
    R: for<'de> FromEdn<'de>,
    C: for<'de> FromEdn<'de>,
{
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "2-element vector [ref constraint]",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    if v.len() != 2 {
        return Err(DeError::Custom(format!(
            "{} entity leaf must have 2 elements, got {}",
            context,
            v.len()
        )));
    }
    Ok((R::from_edn(&v[0])?, C::from_edn(&v[1])?))
}

fn entity_leaf_edn<R: ToEdn, C: ToEdn>(key: &str, r: &R, c: &C) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(key.into())),
        Edn::Vector(vec![r.to_edn(), c.to_edn()].into()),
    );
    Edn::Map(m)
}

fn parse_constraint_vec(
    edn: &Edn<'_>,
    context: &'static str,
) -> Result<Vec<ConstraintDsl>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of constraints",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    v.iter().map(ConstraintDsl::from_edn).collect()
}

fn combinator_edn(key: &str, xs: &[ConstraintDsl]) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(key.into())),
        Edn::Vector(xs.iter().map(|c| c.to_edn()).collect::<Vec<_>>().into()),
    );
    Edn::Map(m)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnKeyword};
    use umol_perm::{Orientation, Permutation};

    use super::*;
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraint, BondConstraint, DativeBondConstraint, FluxionalityAst,
        StereoLigandPair, LigandSymmetryAst, MulticenterValenceAst, OrientedLigandPermutation,
        LigandPermutation, RelationalConstraint, StereoAtomConstraint, StereogenicityAst, TopicityAst, TopicityRelationAst,
    };
    use crate::ast::id::StereoLigandId;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::operators::MemOp;
    use crate::ast::stereo::{CisTransStereoAst, StereoCosetAst, StereoKind, Stereogenicity, TetrahedralStereoAst, Topicity};
    use crate::ast::value::ValueAst;

    #[fixture]
    fn meta_with_atom_id() -> Metadata {
        Metadata::new().with_atom_id(AtomId(2), "c1")
    }

    #[rstest]
    #[case::int(Edn::Int(3), AtomRef::Index(3))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("c1".into())), AtomRef::Id("c1".into()))]
    fn test_atom_ref_from_edn(#[case] input: Edn<'static>, #[case] expected: AtomRef) {
        assert_eq!(AtomRef::from_edn(&input).unwrap(), expected);
    }

    #[rstest]
    fn test_atom_ref_from_edn_rejects_other_kinds() {
        let err = AtomRef::from_edn(&Edn::Str("x".into())).unwrap_err();
        assert!(matches!(
            err,
            DeError::TypeMismatch {
                expected: "atom ref (int or keyword)",
                ..
            }
        ));
    }

    #[rstest]
    #[case::index(AtomRef::Index(5), Edn::Int(5))]
    #[case::id(AtomRef::Id("c1".into()), Edn::Keyword(EdnKeyword::owned("c1".into())))]
    fn test_atom_ref_to_edn(#[case] input: AtomRef, #[case] expected: Edn<'static>) {
        assert_eq!(input.to_edn(), expected);
    }

    #[rstest]
    #[case::int("3", AtomRef::Index(3))]
    #[case::keyword(":c1", AtomRef::Id("c1".into()))]
    fn test_atom_ref_roundtrip_edn_string(#[case] input: &str, #[case] expected: AtomRef) {
        let tree = read_string(input).unwrap();
        let parsed = AtomRef::from_edn(&tree).unwrap();
        assert_eq!(parsed, expected);
        let rendered = parsed.to_edn();
        let reparsed = AtomRef::from_edn(&rendered).unwrap();
        assert_eq!(reparsed, expected);
    }

    #[rstest]
    fn test_atom_ref_from_ast_uses_id_when_present(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomId(2), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Id("c1".into()));
    }

    #[rstest]
    fn test_atom_ref_from_ast_falls_back_to_index_without_id(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomId(4), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Index(4));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_id(meta_with_atom_id: Metadata) {
        let id = AtomRef::Id("c1".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(id, AtomId(2));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_index(meta_with_atom_id: Metadata) {
        let id = AtomRef::Index(3).into_ast(5, &meta_with_atom_id).unwrap();
        assert_eq!(id, AtomId(3));
    }

    #[rstest]
    fn test_atom_ref_into_ast_out_of_range_index(meta_with_atom_id: Metadata) {
        let err = AtomRef::Index(9)
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "9".into(),
            }
        );
    }

    #[rstest]
    fn test_atom_ref_into_ast_unknown_id(meta_with_atom_id: Metadata) {
        let err = AtomRef::Id("nope".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "nope".into(),
            }
        );
    }

    #[fixture]
    fn full_counts() -> EntityCounts {
        EntityCounts {
            atom_count: 10,
            bond_count: 10,
            dative_bond_count: 10,
            aromatic_system_count: 10,
            multicenter_bond_count: 10,
            noncovalent_bond_count: 10,
            stereo_atom_count: 10,
            stereo_bond_count: 10,
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: ValueAst::Lit(0) }, "{:charge-sum {:atoms [0 1] :sum 0}}")]
    #[case::charge_sum_all(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }, "{:charge-sum {:sum 0}}")]
    #[case::spin_sum(MoleculeConstraint::SpinSum { atoms: Some(vec![AtomId(0)]), spin: (1_u8, 2_u8).into() },
        "{:spin-sum {:atoms [0] :spin {:unpaired 1 :multiplicity 2}}}")]
    #[case::spin_sum_all(MoleculeConstraint::SpinSum { atoms: None, spin: (0_u8, 1_u8).into() }, "{:spin-sum {:spin {:unpaired 0 :multiplicity 1}}}")]
    #[case::valence(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(1)]), sum: ValueAst::Lit(4) },
        "{:bond-order-sum {:bonds [0 1] :sum 4}}")]
    #[case::bond_order_sum_all(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(0) }, "{:bond-order-sum {:sum 0}}")]
    #[case::connected(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2)]) }, "{:connected {:atoms [0 1 2]}}")]
    #[case::connected_all(MoleculeConstraint::Connected { atoms: None }, "{:connected {}}")]
    fn test_molecule_constraint_dsl_roundtrip(
        #[from(full_counts)] counts: EntityCounts,
        #[case] input: MoleculeConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = MoleculeConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = MoleculeConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts, &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraint::Fluxionality(FluxionalityAst { perm: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])) })),
        "{:stereo-atom [0 {:kind :tetrahedral :fluxionality [[0 1]]}]}")]
    #[case::ligand_symmetry(
        Constraint::StereoAtom(StereoAtomId(1), StereoKind::Tetrahedral,
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
                perm: OrientedLigandPermutation { perm: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Improper },
                mem: MemOp::NotIn })),
        "{:stereo-atom [1 {:kind :tetrahedral :ligand-symmetry {:perm [[0 1]] :orientation :improper :member :not-in}}]}")]
    #[case::ligand_symmetry_defaults(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
                perm: OrientedLigandPermutation { perm: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])), orientation: Orientation::Proper },
                mem: MemOp::In })),
        "{:stereo-atom [0 {:kind :tetrahedral :ligand-symmetry {:perm [[0 1]]}}]}")]
    #[case::topicity(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Octahedral,
            StereoAtomConstraint::Topicity(TopicityAst {
                pair: StereoLigandPair::new(StereoLigandId(0), StereoLigandId(1)),
                rel: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        "{:stereo-atom [0 {:kind :octahedral :topicity {:pair [0 1] :relation :enantiotopic}}]}")]
    #[case::stereogenicity(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        "{:stereo-atom [0 {:kind :tetrahedral :stereogenicity {:relation :stereogenic}}]}")]
    fn test_constraint_dsl_stereo_atom_roundtrip(
        #[from(full_counts)] counts: EntityCounts,
        #[case] input: Constraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = ConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string(edn_source).unwrap(), "render mismatch");
        let back = ConstraintDsl::from_edn(&edn)
            .unwrap()
            .into_ast(&counts, &meta)
            .unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_molecule_constraint_dsl_rejects_wrong_shape() {
        let err = MoleculeConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_molecule_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = MoleculeConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_molecule_constraint_dsl_charge_sum_rejects_missing_sum() {
        let edn = read_string("{:charge-sum {:atoms [0 1]}}").unwrap();
        let err = MoleculeConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::MissingField { .. }));
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_empty_roundtrip(#[from(full_counts)] counts: EntityCounts) {
        let meta = Metadata::default();
        let anchor = SubPatternAnchor::new();
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        // Empty anchor renders as an empty map.
        assert_eq!(edn, read_string("{}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let back = parsed
            .into_ast_pair(&counts, &meta, &counts, &meta)
            .unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_atoms_roundtrip(#[from(full_counts)] counts: EntityCounts) {
        let meta = Metadata::default();
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomId(3), AtomId(0));
        anchor.push_atom(AtomId(5), AtomId(1));
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("{:atoms [[3 0] [5 1]]}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let back = parsed
            .into_ast_pair(&counts, &meta, &counts, &meta)
            .unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_stereo_roundtrip(#[from(full_counts)] counts: EntityCounts) {
        let meta = Metadata::default();
        let mut anchor = SubPatternAnchor::new();
        anchor.push_stereo_atom(StereoAtomId(2), StereoAtomId(0));
        anchor.push_stereo_bond(StereoBondId(4), StereoBondId(1));
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string("{:stereo-atoms [[2 0]] :stereo-bonds [[4 1]]}").unwrap()
        );
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let back = parsed
            .into_ast_pair(&counts, &meta, &counts, &meta)
            .unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus [[0 0]]}").unwrap();
        let err = SubPatternAnchorDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_rejects_wrong_pair_length() {
        let edn = read_string("{:atoms [[0]]}").unwrap();
        let err = SubPatternAnchorDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_leaf(Constraint::Atom(AtomId(0), AtomConstraint::Valence(ValueAst::Lit(4))), "{:atom [0 {:valence 4}]}")]
    #[case::bond_leaf(Constraint::Bond(BondId(1), BondConstraint::Aromatic), "{:bond [1 :aromatic]}")]
    #[case::dative_bond_leaf_ring_count(Constraint::DativeBond(DativeBondId(0), DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Lit(1))),
        "{:dative-bond [0 {:ring-membership {:count 1}}]}")]
    #[case::dative_bond_leaf_donor(Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(2) }),
        "{:dative-bond-donor [0 2]}")]
    #[case::dative_bond_leaf_acceptor(Constraint::Relational(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(0), atom: AtomId(3) }),
        "{:dative-bond-acceptor [0 3]}")]
    #[case::dative_bond_leaf_parallels(Constraint::Relational(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(2) }),
        "{:dative-bond-parallels [0 2]}")]
    #[case::dative_bond_leaf_donor_satisfies(Constraint::Relational(RelationalConstraint::DativeBondDonorSatisfies { bond: DativeBondId(0),
        predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(3))) }), "{:dative-bond-donor-satisfies [0 {:valence 3}]}")]
    #[case::aromatic_system_leaf_atoms(Constraint::Relational(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }), "{:aromatic-system-atoms [0 [0 1]]}")]
    #[case::aromatic_system_leaf_contains(Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(2) }),
        "{:aromatic-system-contains [0 2]}")]
    #[case::aromatic_system_leaf_all_atoms(Constraint::Relational(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(0),
        predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) }), "{:aromatic-system-all-atoms [0 {:valence 4}]}")]
    #[case::multicenter_leaf_atoms(Constraint::Relational(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)] }), "{:multicenter-bond-atoms [0 [0 1 2]]}")]
    #[case::multicenter_leaf_contains_all(Constraint::Relational(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }), "{:multicenter-bond-contains-all [0 [0 1]]}")]
    #[case::multicenter_leaf_any_atom(Constraint::Relational(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraint::Degree(ValueAst::Lit(3))) }), "{:multicenter-bond-any-atom [0 {:degree 3}]}")]
    #[case::noncovalent_leaf_ends(Constraint::Relational(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(3)] }),
        "{:noncovalent-bond-ends [0 [0 3]]}")]
    #[case::noncovalent_leaf_contains(Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(2) }),
        "{:noncovalent-bond-contains [0 2]}")]
    #[case::noncovalent_leaf_ends_satisfy(Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
        predicates: [Box::new(AtomConstraint::Valence(ValueAst::Lit(2))), Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))] }),
        "{:noncovalent-bond-ends-satisfy [0 [{:valence 2} {:valence 3}]]}")]
    #[case::molecule_connected(Constraint::Molecule(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(1)]) }), "{:connected {:atoms [0 1]}}")]
    #[case::molecule_charge_sum(Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: ValueAst::Lit(0) }),
        "{:charge-sum {:atoms [0 1] :sum 0}}")]
    #[case::molecule_spin_sum(Constraint::Molecule(MoleculeConstraint::SpinSum { atoms: Some(vec![AtomId(0)]), spin: (1_u8, 2_u8).into() }),
        "{:spin-sum {:atoms [0] :spin {:unpaired 1 :multiplicity 2}}}")]
    #[case::molecule_bond_order_sum(Constraint::Molecule(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(1)]), sum: ValueAst::Lit(4) }),
        "{:bond-order-sum {:bonds [0 1] :sum 4}}")]
    #[case::molecule_sub_pattern(Constraint::Molecule(MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) }),
        "{:sub-pattern {:anchor {} :pattern {:atoms [] :bonds []}}}")]
    #[case::not(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraint::Valence(ValueAst::Lit(3))))), "{:not {:atom [0 {:valence 3}]}}")]
    #[case::and(Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraint::Valence(ValueAst::Lit(4))), Constraint::Bond(BondId(0), BondConstraint::Aromatic)]),
        "{:and [{:atom [0 {:valence 4}]} {:bond [0 :aromatic]}]}")]
    #[case::or(Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraint::Degree(ValueAst::Lit(3))), Constraint::Atom(AtomId(0), AtomConstraint::Degree(ValueAst::Lit(4)))]),
        "{:or [{:atom [0 {:degree 3}]} {:atom [0 {:degree 4}]}]}")]
    fn test_constraint_dsl_roundtrip(
        #[from(full_counts)] counts: EntityCounts,
        #[case] input: Constraint,
        #[case] edn_source: &str,
    ) {
        let meta = Metadata::default();
        let dsl = ConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = ConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts, &meta).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::Valence(ValueAst::Lit(4)), "{:valence 4}")]
    #[case::aromatic_valence_not_aromatic(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic), "{:aromatic-valence :not-aromatic}")]
    #[case::aromatic_valence_aromatic(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(6))), "{:aromatic-valence {:aromatic 6}}")]
    #[case::aromatic_valence_undetermined(AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined), "{:aromatic-valence :undetermined}")]
    #[case::multicenter_valence_not_multicenter(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter), "{:multicenter-valence :not-multicenter}")]
    #[case::multicenter_valence_multicenter(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(3))), "{:multicenter-valence {:multicenter 3}}")]
    #[case::multicenter_valence_undetermined(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined), "{:multicenter-valence :undetermined}")]
    #[case::donated_pairs(AtomConstraint::DonatedPairs(ValueAst::Lit(1)), "{:donated-pairs 1}")]
    #[case::accepted_pairs(AtomConstraint::AcceptedPairs(ValueAst::Lit(2)), "{:accepted-pairs 2}")]
    #[case::degree(AtomConstraint::Degree(ValueAst::Lit(3)), "{:degree 3}")]
    #[case::total_degree(AtomConstraint::TotalDegree(ValueAst::Lit(4)), "{:total-degree 4}")]
    #[case::ring_degree(AtomConstraint::RingDegree(ValueAst::Lit(2)), "{:ring-degree 2}")]
    #[case::ring_valence(AtomConstraint::RingValence(ValueAst::Lit(3)), "{:ring-valence 3}")]
    #[case::total_valence(AtomConstraint::TotalValence(ValueAst::Lit(5)), "{:total-valence 5}")]
    #[case::total_hydrogens(AtomConstraint::TotalHydrogens(ValueAst::Lit(3)), "{:total-hydrogens 3}")]
    #[case::ring_membership_all(AtomConstraint::ring_membership(RingScope::All, ValueAst::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_size(AtomConstraint::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::tetrahedral_stereo_not_stereo(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::NotStereo), "{:tetrahedral-stereo :not-stereo}")]
    #[case::tetrahedral_stereo_lit(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1))), "{:tetrahedral-stereo {:stereo 1}}")]
    #[case::tetrahedral_stereo_set(AtomConstraint::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::lit_set([1, 2]))), "{:tetrahedral-stereo {:stereo [1 2]}}")]
    fn test_atom_constraint_dsl_roundtrip(
        #[case] input: AtomConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = AtomConstraintDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = AtomConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&());
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, "{:bond [0 :aromatic]}")]
    #[case::ring_membership_all(BondConstraint::ring_membership(RingScope::All, ValueAst::Lit(1)), "{:bond [0 {:ring-membership {:count 1}}]}")]
    #[case::ring_membership_size(BondConstraint::ring_membership(RingScope::Size(6), 1), "{:bond [0 {:ring-membership {:size 6 :count 1}}]}")]
    #[case::cis_trans_stereo_not_stereo(BondConstraint::CisTransStereo(CisTransStereoAst::NotStereo), "{:bond [0 {:cis-trans-stereo :not-stereo}]}")]
    #[case::cis_trans_stereo_lit(BondConstraint::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::Lit(1))), "{:bond [0 {:cis-trans-stereo {:stereo 1}}]}")]
    fn test_bond_constraint_dsl_roundtrip(
        #[from(full_counts)] counts: EntityCounts,
        #[case] input: BondConstraint,
        #[case] edn_source: &str,
    ) {
        let wrapped = Constraint::Bond(BondId(0), input.clone());
        let meta = Metadata::default();
        let dsl = ConstraintDsl::from_ast(&wrapped, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = ConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts, &meta).unwrap();
        assert_eq!(back, wrapped);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_unknown_key("{:atom [0 {:bogus 1}]}")]
    #[case::bond_unknown_key("{:bond [0 {:bogus 1}]}")]
    #[case::bond_unknown_keyword("{:bond [0 :bogus]}")]
    #[case::bond_wrong_type("{:bond [0 42]}")]
    #[case::dative_unknown_key("{:dative-bond [0 {:bogus 1}]}")]
    #[case::aromatic_system_unknown_key("{:aromatic-system [0 {:bogus 1}]}")]
    #[case::multicenter_unknown_key("{:multicenter-bond [0 {:bogus 1}]}")]
    #[case::noncovalent_unknown_key("{:noncovalent-bond [0 {:bogus 1}]}")]
    #[case::aromatic_valence_unknown_keyword("{:atom [0 {:aromatic-valence :bogus}]}")]
    #[case::aromatic_valence_unknown_key("{:atom [0 {:aromatic-valence {:bogus 1}}]}")]
    #[case::multicenter_valence_unknown_keyword("{:atom [0 {:multicenter-valence :bogus}]}")]
    #[case::multicenter_valence_unknown_key("{:atom [0 {:multicenter-valence {:bogus 1}}]}")]
    fn test_constraint_dsl_rejects_invalid_subvariant(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let err = ConstraintDsl::from_edn(&edn).unwrap_err();
        // Every case must fail; the specific error variant varies.
        assert!(
            matches!(
                err,
                DeError::UnknownField { .. }
                    | DeError::TypeMismatch { .. }
                    | DeError::Custom(_)
                    | DeError::MissingField { .. },
            ),
            "unexpected error: {err:?}",
        );
    }

    #[rstest]
    fn test_constraint_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = ConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_constraints_dsl_empty_roundtrip(#[from(full_counts)] counts: EntityCounts) {
        let meta = Metadata::default();
        let cs = Constraints::new();
        let dsl = ConstraintsDsl::from_ast(&cs, &meta).unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("[]").unwrap());
        let parsed = ConstraintsDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts, &meta).unwrap();
        assert_eq!(back, cs);
    }

    #[rstest]
    fn test_constraints_dsl_multi_roundtrip(#[from(full_counts)] counts: EntityCounts) {
        let meta = Metadata::default();
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(
            AtomId(0),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        ));
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        }));
        let dsl = ConstraintsDsl::from_ast(&cs, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected =
            read_string("[{:atom [0 {:valence 4}]} {:connected {:atoms [0 1]}}]").unwrap();
        assert_eq!(edn, expected);
        let parsed = ConstraintsDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&counts, &meta).unwrap();
        assert_eq!(back, cs);
    }
}
