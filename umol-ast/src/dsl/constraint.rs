//! Tree-shaped constraint DSLs.
//!
//! Boundary types between the AST `Constraint` tree and its EDN form. Refs in
//! the tree carry either an integer index or a symbolic id; resolution to /
//! from the `AtomId` / `BondId` / ... on the AST is a separate fallible
//! step that consults the surrounding `MoleculeMetadata`.

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_perm::{Orientation, Permutation};

use super::aromatic::AromaticSystemConstraintDsl;
use super::atom::{AromaticValenceDsl, AtomConstraintDsl, MulticenterValenceDsl};
use super::bond::BondConstraintDsl;
use super::boolean::read_boolean_dsl;
use super::dative::DativeBondConstraintDsl;
use super::edn_utils::{
    consume_single_key_map_close, eof_err, missing, read_map, read_single_key_map_header, read_vec,
    unexpected_byte_kind,
};
use super::error::ParseError;
use super::molecule::{read_molecule_input, Metadata, MoleculeDsl, MoleculeMetadata};
use super::multicenter::MulticenterBondConstraintDsl;
use super::namespace::{MoleculeNamespace, Namespace};
use super::noncovalent::NoncovalentBondConstraintDsl;
use super::refs::{
    read_aromatic_system_ref, read_atom_ref, read_bond_ref, read_dative_bond_ref,
    read_multicenter_bond_ref, read_noncovalent_bond_ref, read_stereo_atom_ref,
    read_stereo_bond_ref, AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef,
    NoncovalentBondRef, StereoAtomRef, StereoBondRef,
};
use super::relational::{RelationalConstraintDsl, RELATIONAL_KEYS};
use super::stereo::{
    coset_lit, parse_stereo_coset, read_stereogenicity_relation, read_topicity_relation,
    stereo_kind_from_name, RelationValue, StereoAtomConstraintDsl, StereoBondConstraintDsl,
    StereoCosetDsl,
};
use super::value::{parse_value, ValueDsl};
use crate::ast::boolean::BooleanAst;
use crate::ast::constraint::{
    AromaticValenceAst, AtomConstraintAst, BondConstraintAst, Constraint, Constraints,
    FluxionalityAst, LigandPermutation, LigandSymmetryAst, MoleculeConstraint,
    MulticenterValenceAst, OrientedLigandPermutation, RingMembershipAst, RingScope,
    StereoAtomConstraintAst, StereoBondConstraintAst, StereoLigandPair, StereogenicityAst,
    SubPatternAnchor, TopicityAst,
};
use crate::ast::id::{AtomId, BondId, StereoLigandPosition};
use crate::ast::molecule::MoleculeAst;
use crate::ast::spin::SpinStateAst;
use crate::ast::stereo::{CisTransStereoAst, StereoCosetAst, StereoKind, TetrahedralStereoAst};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

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
        let count =
            count.ok_or_else(|| DeError::Custom("ring-membership missing :count".into()))?;
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
        pub(super) fn $name(de: &mut EdnStreamDeserializer<'_>) -> Result<$ast, EdnError> {
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
        "valence" => AtomConstraintAst::Valence(read_value_dsl(de)?.into_ast(&())),
        "total-valence" => AtomConstraintAst::TotalValence(read_value_dsl(de)?.into_ast(&())),
        "aromatic-valence" => {
            AtomConstraintAst::AromaticValence(read_aromatic_valence_dsl(de)?.into_ast(&()))
        }
        "multicenter-valence" => {
            AtomConstraintAst::MulticenterValence(read_multicenter_valence_dsl(de)?.into_ast(&()))
        }
        "donated-pairs" => AtomConstraintAst::DonatedPairs(read_value_dsl(de)?.into_ast(&())),
        "accepted-pairs" => AtomConstraintAst::AcceptedPairs(read_value_dsl(de)?.into_ast(&())),
        "degree" => AtomConstraintAst::Degree(read_value_dsl(de)?.into_ast(&())),
        "total-degree" => AtomConstraintAst::TotalDegree(read_value_dsl(de)?.into_ast(&())),
        "ring-degree" => AtomConstraintAst::RingDegree(read_value_dsl(de)?.into_ast(&())),
        "ring-valence" => AtomConstraintAst::RingValence(read_value_dsl(de)?.into_ast(&())),
        "total-hydrogens" => AtomConstraintAst::TotalHydrogens(read_value_dsl(de)?.into_ast(&())),
        "ring-membership" => AtomConstraintAst::RingMembership(read_ring_membership_dsl(de)?),
        "tetrahedral-stereo" => {
            AtomConstraintAst::TetrahedralStereo(read_tetrahedral_stereo_dsl(de)?)
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
        b'{' => {
            let key = read_single_key_map_header(de)?;
            let c = match key.as_str() {
                "aromatic" => BondConstraintAst::Aromatic(read_boolean_dsl(de)?.0),
                "ring-membership" => {
                    BondConstraintAst::RingMembership(read_ring_membership_dsl(de)?)
                }
                "cis-trans-stereo" => {
                    BondConstraintAst::CisTransStereo(read_cis_trans_stereo_dsl(de)?)
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
            expected: "{:aromatic …} / {:ring-membership …} / {:cis-trans-stereo …}",
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
        b'{' => {
            let key = read_single_key_map_header(de)?;
            let c = match key.as_str() {
                "aromatic" => DativeBondConstraintDsl::Aromatic(read_boolean_dsl(de)?.0),
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
            expected: "{:aromatic …} / {:ring-membership …}",
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
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondConstraintDsl, EdnError> {
    let key = read_single_key_map_header(de)?;
    let c = match key.as_str() {
        "intramolecular" => NoncovalentBondConstraintDsl::Intramolecular(read_boolean_dsl(de)?.0),
        other => {
            return Err(DeError::UnknownField {
                key: other.to_string(),
                path: vec!["noncovalent-bond-constraint".into()],
            }
            .into());
        }
    };
    consume_single_key_map_close(de, "noncovalent-bond-constraint")?;
    Ok(c)
}

/// Membership polarity `:in` / `:not-in`. Shared by `#p` and the `#o`/`#g` relations.
/// A permutation as a vector of disjoint cycles `[[0 1 2] [3 4]]`; degree from the
/// stereo kind.
fn read_permutation(
    de: &mut EdnStreamDeserializer<'_>,
    degree: usize,
) -> Result<Permutation, EdnError> {
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    de.consume_byte(b'[')?;
    while !de.try_consume_byte(b']')? {
        de.consume_byte(b'[')?;
        let mut cycle = Vec::new();
        while !de.try_consume_byte(b']')? {
            let n = de.read_i64()?;
            let point = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                value: n.to_string(),
                target: "ligand position",
                path: Vec::new(),
            })?;
            cycle.push(point);
        }
        cycles.push(cycle);
    }
    Permutation::from_cycles(degree, &cycles)
        .map_err(|error| DeError::Custom(error.to_string()).into())
}

/// The `:relation` value: `:undetermined`, one keyword, a keyword vector (`LitSet`),
/// or a `{:not-in [members]}` complement map (`NotSet`).
fn read_relation_value(de: &mut EdnStreamDeserializer<'_>) -> Result<RelationValue, EdnError> {
    match de.peek_byte()? {
        Some(b'[') => Ok(RelationValue::Many(read_vec(de, |de| {
            Ok(de.read_keyword_name()?.into_owned())
        })?)),
        Some(b'{') => {
            let mut complement = None;
            read_map(de, |de, key| {
                match key {
                    "not-in" => {
                        complement =
                            Some(read_vec(de, |de| Ok(de.read_keyword_name()?.into_owned()))?);
                    }
                    other => {
                        return Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["relation".into()],
                        }
                        .into())
                    }
                }
                Ok(())
            })?;
            let complement = complement.ok_or_else(|| {
                DeError::Custom("relation complement missing :not-in".to_string())
            })?;
            Ok(RelationValue::NotIn(complement))
        }
        _ => {
            let kw = de.read_keyword_name()?.into_owned();
            Ok(if kw == "undetermined" {
                RelationValue::Undetermined
            } else {
                RelationValue::One(kw)
            })
        }
    }
}

fn read_ligand_symmetry(
    de: &mut EdnStreamDeserializer<'_>,
    kind: StereoKind,
) -> Result<LigandSymmetryAst, EdnError> {
    let mut permutation = None;
    let mut orientation = Orientation::Proper;
    let mut invariant = BooleanAst::Lit(true);
    read_map(de, |de, key| {
        match key {
            "permutation" => permutation = Some(read_permutation(de, kind.degree())?),
            "orientation" => {
                orientation = match de.read_keyword_name()?.as_ref() {
                    "proper" => Orientation::Proper,
                    "improper" => Orientation::Improper,
                    other => {
                        return Err(DeError::Custom(format!(
                            "expected :proper | :improper, got :{other}"
                        ))
                        .into())
                    }
                }
            }
            "invariant" => invariant = read_boolean_dsl(de)?.0,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["ligand-symmetry".into()],
                }
                .into())
            }
        }
        Ok(())
    })?;
    let permutation = permutation
        .ok_or_else(|| DeError::Custom("ligand-symmetry missing :permutation".to_string()))?;
    Ok(LigandSymmetryAst {
        permutation: OrientedLigandPermutation {
            permutation: LigandPermutation(permutation),
            orientation,
        },
        invariant,
    })
}

fn read_fluxionality(
    de: &mut EdnStreamDeserializer<'_>,
    kind: StereoKind,
) -> Result<FluxionalityAst, EdnError> {
    let mut permutation = None;
    let mut active = BooleanAst::Lit(true);
    read_map(de, |de, key| {
        match key {
            "permutation" => permutation = Some(read_permutation(de, kind.degree())?),
            "active" => active = read_boolean_dsl(de)?.0,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["fluxionality".into()],
                }
                .into())
            }
        }
        Ok(())
    })?;
    let permutation = permutation
        .ok_or_else(|| DeError::Custom("fluxionality missing :permutation".to_string()))?;
    Ok(FluxionalityAst {
        permutation: LigandPermutation(permutation),
        active,
    })
}

fn read_topicity(de: &mut EdnStreamDeserializer<'_>) -> Result<TopicityAst, EdnError> {
    let mut pair = None;
    let mut value = None;
    read_map(de, |de, key| {
        match key {
            "pair" => {
                let v = read_vec(de, |de| Ok(de.read_i64()?))?;
                let [a, b]: [i64; 2] = v[..].try_into().map_err(|_| {
                    DeError::Custom("topicity :pair must have 2 positions".to_string())
                })?;
                pair = Some(StereoLigandPair::new(
                    StereoLigandPosition(a as u32),
                    StereoLigandPosition(b as u32),
                ));
            }
            "relation" => value = Some(read_relation_value(de)?),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["topicity".into()],
                }
                .into())
            }
        }
        Ok(())
    })?;
    let pair = pair.ok_or_else(|| DeError::Custom("topicity missing :pair".to_string()))?;
    let value = value.ok_or_else(|| DeError::Custom("topicity missing :relation".to_string()))?;
    Ok(TopicityAst {
        pair,
        relation: read_topicity_relation(value)?,
    })
}

fn read_stereogenicity(de: &mut EdnStreamDeserializer<'_>) -> Result<StereogenicityAst, EdnError> {
    let mut value = None;
    read_map(de, |de, key| {
        match key {
            "relation" => value = Some(read_relation_value(de)?),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["stereogenicity".into()],
                }
                .into())
            }
        }
        Ok(())
    })?;
    let value =
        value.ok_or_else(|| DeError::Custom("stereogenicity missing :relation".to_string()))?;
    Ok(read_stereogenicity_relation(value)?)
}

/// Stream a stereo constraint as the positional 2-vector `[<kind> {<key> <value>}]`:
/// kind first (container-fixed) → degree known → the single-key payload streamed.
/// The per-key dispatch is the only per-collection part.
macro_rules! read_stereo_constraint_dsl {
    ($name:ident, $constraint:ident, $dsl:ident, $context:literal) => {
        fn $name(de: &mut EdnStreamDeserializer<'_>) -> Result<$dsl, EdnError> {
            de.consume_byte(b'[')?;
            let kind = stereo_kind_from_name(de.read_keyword_name()?.as_ref())?;
            let key = read_single_key_map_header(de)?;
            let constraint = match key.as_str() {
                "ligand-symmetry" => $constraint::LigandSymmetry(read_ligand_symmetry(de, kind)?),
                "fluxionality" => $constraint::Fluxionality(read_fluxionality(de, kind)?),
                "topicity" => $constraint::Topicity(read_topicity(de)?),
                "stereogenicity" => $constraint::Stereogenicity(read_stereogenicity(de)?),
                other => {
                    return Err(DeError::Custom(format!(
                        "unknown stereo constraint keyword :{other}"
                    ))
                    .into())
                }
            };
            consume_single_key_map_close(de, $context)?;
            de.consume_byte(b']')?;
            Ok($dsl(kind, constraint))
        }
    };
}

read_stereo_constraint_dsl! {
    read_stereo_atom_constraint_dsl, StereoAtomConstraintAst, StereoAtomConstraintDsl,
    "stereo-atom-constraint"
}
read_stereo_constraint_dsl! {
    read_stereo_bond_constraint_dsl, StereoBondConstraintAst, StereoBondConstraintDsl,
    "stereo-bond-constraint"
}

fn read_atom_ref_vec(de: &mut EdnStreamDeserializer<'_>) -> Result<Vec<AtomRef>, EdnError> {
    read_vec(de, read_atom_ref)
}

fn read_atom_ref_pair(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<[AtomRef; 2], EdnError> {
    de.consume_byte(b'[')?;
    let first = read_atom_ref(de)?;
    let second = read_atom_ref(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom(format!("{}: expected 2 elements", context)).into());
    }
    Ok([first, second])
}

fn read_atom_constraint_pair(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<[AtomConstraintDsl; 2], EdnError> {
    de.consume_byte(b'[')?;
    let first = read_atom_constraint_dsl(de)?;
    let second = read_atom_constraint_dsl(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom(format!("{}: expected 2 elements", context)).into());
    }
    Ok([first, second])
}

pub(super) fn read_relational_constraint_dsl(
    de: &mut EdnStreamDeserializer<'_>,
    key: &str,
) -> Result<RelationalConstraintDsl, EdnError> {
    use RelationalConstraintDsl as R;
    de.consume_byte(b'[')?;
    let c = match key {
        "dative-bond-donors" => R::DativeBondDonors {
            bond: read_dative_bond_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "dative-bond-donor" => R::DativeBondDonor {
            bond: read_dative_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "dative-bond-contains-all-donors" => R::DativeBondContainsAllDonors {
            bond: read_dative_bond_ref(de)?,
            atoms: read_atom_ref_vec(de)?,
        },
        "dative-bond-all-donors" => R::DativeBondAllDonors {
            bond: read_dative_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "dative-bond-any-donor" => R::DativeBondAnyDonor {
            bond: read_dative_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "dative-bond-acceptor" => R::DativeBondAcceptor {
            bond: read_dative_bond_ref(de)?,
            atom: read_atom_ref(de)?,
        },
        "dative-bond-acceptor-satisfies" => R::DativeBondAcceptorSatisfies {
            bond: read_dative_bond_ref(de)?,
            predicate: Box::new(read_atom_constraint_dsl(de)?),
        },
        "dative-bond-parallels" => R::DativeBondParallels {
            dative: read_dative_bond_ref(de)?,
            parallel: read_bond_ref(de)?,
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

fn read_ref_pair<R1, R2>(
    de: &mut EdnStreamDeserializer<'_>,
    read_first: fn(&mut EdnStreamDeserializer<'_>) -> Result<R1, EdnError>,
    read_second: fn(&mut EdnStreamDeserializer<'_>) -> Result<R2, EdnError>,
) -> Result<(R1, R2), EdnError> {
    de.consume_byte(b'[')?;
    let first = read_first(de)?;
    let second = read_second(de)?;
    if !de.try_consume_byte(b']')? {
        return Err(DeError::Custom("anchor pair must have 2 elements".into()).into());
    }
    Ok((first, second))
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
                        let input = read_molecule_input(d)?;
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
/// the identity; the EDN bridge wraps with empty `MoleculeMetadata` so refs render
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

fn denote_atom_subset<M: Metadata>(atoms: &Option<Vec<AtomId>>, meta: &M) -> Option<Vec<AtomRef>> {
    atoms
        .as_ref()
        .map(|v| v.iter().map(|&a| AtomRef::denote(a, meta)).collect())
}

fn denote_bond_subset<M: Metadata>(bonds: &Option<Vec<BondId>>, meta: &M) -> Option<Vec<BondRef>> {
    bonds
        .as_ref()
        .map(|v| v.iter().map(|&b| BondRef::denote(b, meta)).collect())
}

fn resolve_atom_subset<N: Namespace>(
    atoms: Option<Vec<AtomRef>>,
    namespace: &N,
) -> Result<Option<Vec<AtomId>>, ParseError> {
    atoms
        .map(|v| {
            v.into_iter()
                .map(|r| r.resolve(namespace))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
}

fn resolve_bond_subset<N: Namespace>(
    bonds: Option<Vec<BondRef>>,
    namespace: &N,
) -> Result<Option<Vec<BondId>>, ParseError> {
    bonds
        .map(|v| {
            v.into_iter()
                .map(|r| r.resolve(namespace))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
}

impl MoleculeConstraintDsl {
    pub(crate) fn from_ast<M: Metadata>(
        c: &MoleculeConstraint,
        meta: &M,
    ) -> Result<Self, ParseError> {
        Ok(match c {
            MoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum {
                atoms: denote_atom_subset(atoms, meta),
                sum: ValueDsl::from_ast(sum, &()),
            },
            MoleculeConstraint::SpinSum { atoms, spin } => Self::SpinSum {
                atoms: denote_atom_subset(atoms, meta),
                spin: spin.clone(),
            },
            MoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum {
                bonds: denote_bond_subset(bonds, meta),
                sum: ValueDsl::from_ast(sum, &()),
            },
            MoleculeConstraint::Connected { atoms } => Self::Connected {
                atoms: denote_atom_subset(atoms, meta),
            },
            MoleculeConstraint::SubPattern { anchor, pattern } => {
                let pattern_meta = MoleculeMetadata::default();
                let anchor_dsl = SubPatternAnchorDsl::from_ast_pair(anchor, meta, &pattern_meta);
                Self::SubPattern {
                    anchor: anchor_dsl,
                    pattern: pattern.clone(),
                }
            }
        })
    }

    pub(crate) fn into_ast<N: Namespace>(
        self,
        namespace: &N,
    ) -> Result<MoleculeConstraint, ParseError> {
        Ok(match self {
            Self::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms: resolve_atom_subset(atoms, namespace)?,
                sum: sum.into_ast(&()),
            },
            Self::SpinSum { atoms, spin } => MoleculeConstraint::SpinSum {
                atoms: resolve_atom_subset(atoms, namespace)?,
                spin,
            },
            Self::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds: resolve_bond_subset(bonds, namespace)?,
                sum: sum.into_ast(&()),
            },
            Self::Connected { atoms } => MoleculeConstraint::Connected {
                atoms: resolve_atom_subset(atoms, namespace)?,
            },
            Self::SubPattern { anchor, pattern } => {
                let pattern_namespace = MoleculeNamespace::from_ast(&pattern);
                let anchor_ast = anchor.into_ast_pair(namespace, &pattern_namespace)?;
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
                let (pattern_ast, pattern_meta) = pattern_dsl.into_parts();
                if !pattern_meta.is_empty() {
                    return Err(DeError::Custom(
                        "sub-pattern molecule must be anonymous: no :id or :atom-aliases".into(),
                    ));
                }
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
                let pattern_dsl =
                    MoleculeDsl::from_parts((**pattern).clone(), MoleculeMetadata::default());
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
/// molecule's `MoleculeMetadata`; pattern-side refs resolve against the pattern
/// molecule's `MoleculeMetadata`.
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
    pub fn from_ast_pair<M: Metadata>(
        anchor: &SubPatternAnchor,
        target_meta: &M,
        pattern_meta: &MoleculeMetadata,
    ) -> Self {
        Self {
            atoms: anchor
                .atoms()
                .iter()
                .map(|&(t, p)| {
                    (
                        AtomRef::denote(t, target_meta),
                        AtomRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            bonds: anchor
                .bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        BondRef::denote(t, target_meta),
                        BondRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            dative_bonds: anchor
                .dative_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        DativeBondRef::denote(t, target_meta),
                        DativeBondRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            aromatic_systems: anchor
                .aromatic_systems()
                .iter()
                .map(|&(t, p)| {
                    (
                        AromaticSystemRef::denote(t, target_meta),
                        AromaticSystemRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            multicenter_bonds: anchor
                .multicenter_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        MulticenterBondRef::denote(t, target_meta),
                        MulticenterBondRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            noncovalent_bonds: anchor
                .noncovalent_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        NoncovalentBondRef::denote(t, target_meta),
                        NoncovalentBondRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            stereo_atoms: anchor
                .stereo_atoms()
                .iter()
                .map(|&(t, p)| {
                    (
                        StereoAtomRef::denote(t, target_meta),
                        StereoAtomRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
            stereo_bonds: anchor
                .stereo_bonds()
                .iter()
                .map(|&(t, p)| {
                    (
                        StereoBondRef::denote(t, target_meta),
                        StereoBondRef::denote(p, pattern_meta),
                    )
                })
                .collect(),
        }
    }

    /// Resolve to an AST anchor. `target_*` carry outer-molecule counts +
    /// Resolve to an AST anchor. `host` is the enclosing molecule/reaction namespace (target side);
    /// `pattern_*` are the pattern molecule's counts + metadata — a stopgap until the pattern namespace
    /// is built (S3f).
    pub(crate) fn into_ast_pair(
        self,
        host: &impl Namespace,
        pattern: &impl Namespace,
    ) -> Result<SubPatternAnchor, ParseError> {
        let mut anchor = SubPatternAnchor::new();
        for (t, p) in self.atoms {
            anchor.push_atom(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.bonds {
            anchor.push_bond(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.dative_bonds {
            anchor.push_dative_bond(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.aromatic_systems {
            anchor.push_aromatic_system(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.multicenter_bonds {
            anchor.push_multicenter_bond(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.noncovalent_bonds {
            anchor.push_noncovalent_bond(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.stereo_atoms {
            anchor.push_stereo_atom(t.resolve(host)?, p.resolve(pattern)?);
        }
        for (t, p) in self.stereo_bonds {
            anchor.push_stereo_bond(t.resolve(host)?, p.resolve(pattern)?);
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
            Self::NoncovalentBond(r, c) => entity_leaf_edn("noncovalent-bond", r, c),
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
    pub(crate) fn from_ast<M: Metadata>(c: &Constraint, meta: &M) -> Result<Self, ParseError> {
        Ok(match c {
            Constraint::Atom(id, c) => Self::Atom(
                AtomRef::denote(*id, meta),
                AtomConstraintDsl::from_ast(c, &()),
            ),
            Constraint::Bond(id, c) => Self::Bond(
                BondRef::denote(*id, meta),
                BondConstraintDsl::from_ast(c, &()),
            ),
            Constraint::DativeBond(id, c) => Self::DativeBond(
                DativeBondRef::denote(*id, meta),
                DativeBondConstraintDsl::from_ast(c),
            ),
            Constraint::AromaticSystem(id, c) => Self::AromaticSystem(
                AromaticSystemRef::denote(*id, meta),
                AromaticSystemConstraintDsl::from_ast(c),
            ),
            Constraint::MulticenterBond(id, c) => Self::MulticenterBond(
                MulticenterBondRef::denote(*id, meta),
                MulticenterBondConstraintDsl::from_ast(c),
            ),
            Constraint::NoncovalentBond(id, c) => Self::NoncovalentBond(
                NoncovalentBondRef::denote(*id, meta),
                NoncovalentBondConstraintDsl::from_ast(c),
            ),
            Constraint::StereoAtom(id, kind, c) => Self::StereoAtom(
                StereoAtomRef::denote(*id, meta),
                StereoAtomConstraintDsl::from_ast(c, kind),
            ),
            Constraint::StereoBond(id, kind, c) => Self::StereoBond(
                StereoBondRef::denote(*id, meta),
                StereoBondConstraintDsl::from_ast(c, kind),
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

    pub(crate) fn into_ast<N: Namespace>(self, namespace: &N) -> Result<Constraint, ParseError> {
        Ok(match self {
            Self::Atom(r, c) => Constraint::Atom(r.resolve(namespace)?, c.into_ast(&())),
            Self::Bond(r, c) => Constraint::Bond(r.resolve(namespace)?, c.into_ast(&())),
            Self::DativeBond(r, c) => Constraint::DativeBond(r.resolve(namespace)?, c.into_ast()),
            Self::AromaticSystem(r, c) => {
                Constraint::AromaticSystem(r.resolve(namespace)?, c.into_ast())
            }
            Self::MulticenterBond(r, c) => {
                Constraint::MulticenterBond(r.resolve(namespace)?, c.into_ast())
            }
            Self::NoncovalentBond(r, c) => {
                Constraint::NoncovalentBond(r.resolve(namespace)?, c.into_ast())
            }
            Self::StereoAtom(r, StereoAtomConstraintDsl(kind, c)) => {
                Constraint::StereoAtom(r.resolve(namespace)?, kind, c)
            }
            Self::StereoBond(r, StereoBondConstraintDsl(kind, c)) => {
                Constraint::StereoBond(r.resolve(namespace)?, kind, c)
            }
            Self::Relational(r) => Constraint::Relational(r.into_ast(namespace)?),
            Self::Molecule(m) => Constraint::Molecule(m.into_ast(namespace)?),
            Self::And(xs) => Constraint::And(
                xs.into_iter()
                    .map(|c| c.into_ast(namespace))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Or(xs) => Constraint::Or(
                xs.into_iter()
                    .map(|c| c.into_ast(namespace))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Not(c) => Constraint::Not(Box::new(c.into_ast(namespace)?)),
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
    pub(crate) fn from_ast<M: Metadata>(cs: &Constraints, meta: &M) -> Result<Self, ParseError> {
        Ok(Self(
            cs.iter()
                .filter(|c| !c.is_vacuous())
                .map(|c| ConstraintDsl::from_ast(c, meta))
                .collect::<Result<_, _>>()?,
        ))
    }

    pub(crate) fn into_ast<N: Namespace>(self, namespace: &N) -> Result<Constraints, ParseError> {
        let mut out = Constraints::new();
        for c in self.0 {
            out.push(c.into_ast(namespace)?);
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
    use std::fs;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;
    use umol_perm::{Orientation, Permutation};

    use super::super::namespace::MoleculeNamespace;
    use super::*;
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraintAst, BondConstraintAst, DativeBondConstraintAst,
        FluxionalityAst, LigandPermutation, LigandSymmetryAst, MulticenterValenceAst,
        OrientedLigandPermutation, RelationalConstraint, StereoAtomConstraintAst, StereoLigandPair,
        StereogenicityAst, TopicityAst, TopicityRelationAst,
    };
    use crate::ast::id::{
        AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
        StereoBondId,
    };
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::stereo::{
        CisTransStereoAst, StereoCosetAst, StereoKind, Stereogenicity, TetrahedralStereoAst,
        Topicity,
    };
    use crate::ast::value::ValueAst;
    use crate::ast::BooleanAst;

    /// Every `fuzz_constraints` seed must parse as a `ConstraintDsl` or a `ConstraintsDsl` (tree) —
    /// guards the seed corpus against rot as the constraint DSL evolves.
    #[rstest]
    fn test_fuzz_constraints_seeds_valid() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fuzz/seeds/fuzz_constraints");
        let mut failures: Vec<String> = Vec::new();
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            let data = fs::read_to_string(&path).unwrap();
            let Ok(edn) = read_string(&data) else {
                failures.push(format!("{name}: not readable EDN"));
                continue;
            };
            if ConstraintDsl::from_edn(&edn).is_err() && ConstraintsDsl::from_edn(&edn).is_err() {
                failures.push(format!(
                    "{name}: neither ConstraintDsl nor ConstraintsDsl parses"
                ));
            }
        }
        assert!(
            failures.is_empty(),
            "invalid seeds:\n{}",
            failures.join("\n")
        );
    }

    #[rstest]
    #[case::empty_cycle("[[]]", Permutation::identity(4))]
    #[case::single_cycle("[[0 1 2]]", Permutation::from_image(&[1, 2, 0, 3]))]
    #[case::disjoint_cycles("[[0 1] [2 3]]", Permutation::from_image(&[1, 0, 3, 2]))]
    fn test_read_permutation(#[case] input: &str, #[case] expected: Permutation) {
        let mut de = EdnStreamDeserializer::new(input);
        assert_eq!(read_permutation(&mut de, 4), Ok(expected));
        assert_eq!(de.expect_eof(), Ok(()));
    }

    #[rstest]
    #[case::overlap(
        "[[0 1] [1 2]]",
        EdnError::De(DeError::Custom("cycle point 1 occurs more than once".to_string())),
    )]
    #[case::repeated(
        "[[0 1 0]]",
        EdnError::De(DeError::Custom("cycle point 0 occurs more than once".to_string())),
    )]
    #[case::out_of_range(
        "[[0 4]]",
        EdnError::De(DeError::Custom(
            "cycle point 4 at cycle 0, position 1 is outside 0..4".to_string(),
        )),
    )]
    #[case::negative(
        "[[0 -1]]",
        EdnError::De(DeError::OutOfRange {
            value: "-1".to_string(),
            target: "ligand position",
            path: Vec::new(),
        }),
    )]
    fn test_read_permutation_error(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(
            read_permutation(&mut EdnStreamDeserializer::new(input), 4),
            Err(expected),
        );
    }

    /// A namespace with ten entities of each kind, so index refs up to 9 resolve.
    #[fixture]
    fn full_namespace() -> MoleculeNamespace {
        let mut namespace = MoleculeNamespace::default();
        for _ in 0..10 {
            namespace.register_atom(None).unwrap();
        }
        for _ in 0..10 {
            namespace.register_bond(None, AtomId(0), AtomId(1)).unwrap();
            namespace
                .register_dative_bond(None, &[AtomId(0)], AtomId(1))
                .unwrap();
            namespace
                .register_aromatic_system(None, &[AtomId(0)])
                .unwrap();
            namespace
                .register_multicenter_bond(None, &[AtomId(0)])
                .unwrap();
            namespace
                .register_noncovalent_bond(None, AtomId(0), AtomId(1))
                .unwrap();
            namespace
                .register_stereo_atom(None, AtomId(0), &[])
                .unwrap();
            namespace
                .register_stereo_bond(None, BondId(0), &[])
                .unwrap();
        }
        namespace
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
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[case] input: MoleculeConstraint,
        #[case] edn_source: &str,
    ) {
        let meta = MoleculeMetadata::default();
        let dsl = MoleculeConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = MoleculeConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&namespace).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanAst::Lit(true) })),
        "{:stereo-atom [0 [:tetrahedral {:fluxionality {:permutation [[0 1]]}}]]}")]
    #[case::fluxionality_absent(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Fluxionality(FluxionalityAst { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanAst::Lit(false) })),
        "{:stereo-atom [0 [:tetrahedral {:fluxionality {:permutation [[0 1]] :active false}}]]}")]
    #[case::ligand_symmetry(
        Constraint::StereoAtom(StereoAtomId(1), StereoKind::Tetrahedral,
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst {
                permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Improper },
                invariant: BooleanAst::Lit(false) })),
        "{:stereo-atom [1 [:tetrahedral {:ligand-symmetry {:permutation [[0 1]] :orientation :improper :invariant false}}]]}")]
    #[case::ligand_symmetry_defaults(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraintAst::LigandSymmetry(LigandSymmetryAst {
                permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper },
                invariant: BooleanAst::Lit(true) })),
        "{:stereo-atom [0 [:tetrahedral {:ligand-symmetry {:permutation [[0 1]]}}]]}")]
    #[case::topicity(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Octahedral,
            StereoAtomConstraintAst::Topicity(TopicityAst {
                pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                relation: TopicityRelationAst::Lit(Topicity::Enantiotopic) })),
        "{:stereo-atom [0 [:octahedral {:topicity {:pair [0 1] :relation :enantiotopic}}]]}")]
    #[case::stereogenicity(
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
        "{:stereo-atom [0 [:tetrahedral {:stereogenicity {:relation :stereogenic}}]]}")]
    fn test_constraint_dsl_stereo_atom_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[case] input: Constraint,
        #[case] edn_source: &str,
    ) {
        let meta = MoleculeMetadata::default();
        let dsl = ConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string(edn_source).unwrap(), "render mismatch");
        let back = ConstraintDsl::from_edn(&edn)
            .unwrap()
            .into_ast(&namespace)
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
    #[case::named_entity(r#"{:sub-pattern {:anchor {} :pattern {:atoms [[:a "C"]]}}}"#)]
    #[case::aliased(r#"{:sub-pattern {:anchor {} :pattern {:atoms [:x] :atom-aliases [:x "C"]}}}"#)]
    fn test_molecule_constraint_dsl_sub_pattern_rejects_named_pattern(#[case] input: &str) {
        let DeError::Custom(msg) =
            MoleculeConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap_err()
        else {
            panic!("expected a Custom error");
        };
        assert_eq!(
            msg,
            "sub-pattern molecule must be anonymous: no :id or :atom-aliases"
        );
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_empty_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[from(full_namespace)] pattern: MoleculeNamespace,
    ) {
        let meta = MoleculeMetadata::default();
        let anchor = SubPatternAnchor::new();
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        // Empty anchor renders as an empty map.
        assert_eq!(edn, read_string("{}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast_pair(&namespace, &pattern).unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_atoms_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[from(full_namespace)] pattern: MoleculeNamespace,
    ) {
        let meta = MoleculeMetadata::default();
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomId(3), AtomId(0));
        anchor.push_atom(AtomId(5), AtomId(1));
        let dsl = SubPatternAnchorDsl::from_ast_pair(&anchor, &meta, &meta);
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("{:atoms [[3 0] [5 1]]}").unwrap());
        let parsed = SubPatternAnchorDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast_pair(&namespace, &pattern).unwrap();
        assert_eq!(back, anchor);
    }

    #[rstest]
    fn test_sub_pattern_anchor_dsl_stereo_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[from(full_namespace)] pattern: MoleculeNamespace,
    ) {
        let meta = MoleculeMetadata::default();
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
        let back = parsed.into_ast_pair(&namespace, &pattern).unwrap();
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
    #[case::atom_leaf(Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::Lit(4))), "{:atom [0 {:valence 4}]}")]
    #[case::bond_leaf(Constraint::Bond(BondId(1), BondConstraintAst::Aromatic(BooleanAst::Lit(true))), "{:bond [1 {:aromatic true}]}")]
    #[case::dative_bond_leaf_ring_count(Constraint::DativeBond(DativeBondId(0), DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1))),
        "{:dative-bond [0 {:ring-membership {:count 1}}]}")]
    #[case::dative_bond_leaf_donor(Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(2) }),
        "{:dative-bond-donor [0 2]}")]
    #[case::dative_bond_leaf_acceptor(Constraint::Relational(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(0), atom: AtomId(3) }),
        "{:dative-bond-acceptor [0 3]}")]
    #[case::dative_bond_leaf_parallels(Constraint::Relational(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(2) }),
        "{:dative-bond-parallels [0 2]}")]
    #[case::dative_bond_leaf_all_donors(Constraint::Relational(RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::Valence(ValueAst::Lit(3))) }), "{:dative-bond-all-donors [0 {:valence 3}]}")]
    #[case::aromatic_system_leaf_atoms(Constraint::Relational(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }), "{:aromatic-system-atoms [0 [0 1]]}")]
    #[case::aromatic_system_leaf_contains(Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(2) }),
        "{:aromatic-system-contains [0 2]}")]
    #[case::aromatic_system_leaf_all_atoms(Constraint::Relational(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(0),
        predicate: Box::new(AtomConstraintAst::Valence(ValueAst::Lit(4))) }), "{:aromatic-system-all-atoms [0 {:valence 4}]}")]
    #[case::multicenter_leaf_atoms(Constraint::Relational(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)] }), "{:multicenter-bond-atoms [0 [0 1 2]]}")]
    #[case::multicenter_leaf_contains_all(Constraint::Relational(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }), "{:multicenter-bond-contains-all [0 [0 1]]}")]
    #[case::multicenter_leaf_any_atom(Constraint::Relational(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraintAst::Degree(ValueAst::Lit(3))) }), "{:multicenter-bond-any-atom [0 {:degree 3}]}")]
    #[case::noncovalent_leaf_ends(Constraint::Relational(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(3)] }),
        "{:noncovalent-bond-ends [0 [0 3]]}")]
    #[case::noncovalent_leaf_contains(Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(2) }),
        "{:noncovalent-bond-contains [0 2]}")]
    #[case::noncovalent_leaf_ends_satisfy(Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
        predicates: [Box::new(AtomConstraintAst::Valence(ValueAst::Lit(2))), Box::new(AtomConstraintAst::Valence(ValueAst::Lit(3)))] }),
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
    #[case::not(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::Lit(3))))), "{:not {:atom [0 {:valence 3}]}}")]
    #[case::and(Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::Lit(4))), Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))]),
        "{:and [{:atom [0 {:valence 4}]} {:bond [0 {:aromatic true}]}]}")]
    #[case::or(Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::Degree(ValueAst::Lit(3))), Constraint::Atom(AtomId(0), AtomConstraintAst::Degree(ValueAst::Lit(4)))]),
        "{:or [{:atom [0 {:degree 3}]} {:atom [0 {:degree 4}]}]}")]
    fn test_constraint_dsl_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[case] input: Constraint,
        #[case] edn_source: &str,
    ) {
        let meta = MoleculeMetadata::default();
        let dsl = ConstraintDsl::from_ast(&input, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = ConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&namespace).unwrap();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintAst::Valence(ValueAst::Lit(4)), "{:valence 4}")]
    #[case::aromatic_valence_not_aromatic(AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic), "{:aromatic-valence :not-aromatic}")]
    #[case::aromatic_valence_aromatic(AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(6))), "{:aromatic-valence {:aromatic 6}}")]
    #[case::aromatic_valence_undetermined(AtomConstraintAst::AromaticValence(AromaticValenceAst::Undetermined), "{:aromatic-valence :undetermined}")]
    #[case::multicenter_valence_not_multicenter(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter), "{:multicenter-valence :not-multicenter}")]
    #[case::multicenter_valence_multicenter(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(3))), "{:multicenter-valence {:multicenter 3}}")]
    #[case::multicenter_valence_undetermined(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Undetermined), "{:multicenter-valence :undetermined}")]
    #[case::donated_pairs(AtomConstraintAst::DonatedPairs(ValueAst::Lit(1)), "{:donated-pairs 1}")]
    #[case::accepted_pairs(AtomConstraintAst::AcceptedPairs(ValueAst::Lit(2)), "{:accepted-pairs 2}")]
    #[case::degree(AtomConstraintAst::Degree(ValueAst::Lit(3)), "{:degree 3}")]
    #[case::total_degree(AtomConstraintAst::TotalDegree(ValueAst::Lit(4)), "{:total-degree 4}")]
    #[case::ring_degree(AtomConstraintAst::RingDegree(ValueAst::Lit(2)), "{:ring-degree 2}")]
    #[case::ring_valence(AtomConstraintAst::RingValence(ValueAst::Lit(3)), "{:ring-valence 3}")]
    #[case::total_valence(AtomConstraintAst::TotalValence(ValueAst::Lit(5)), "{:total-valence 5}")]
    #[case::total_hydrogens(AtomConstraintAst::TotalHydrogens(ValueAst::Lit(3)), "{:total-hydrogens 3}")]
    #[case::ring_membership_all(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_size(AtomConstraintAst::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::tetrahedral_stereo_not_stereo(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo), "{:tetrahedral-stereo :not-stereo}")]
    #[case::tetrahedral_stereo_lit(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1))), "{:tetrahedral-stereo {:stereo 1}}")]
    #[case::tetrahedral_stereo_set(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::lit_set([1, 2]))), "{:tetrahedral-stereo {:stereo [1 2]}}")]
    fn test_atom_constraint_dsl_roundtrip(
        #[case] input: AtomConstraintAst,
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
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanAst::Lit(true)), "{:bond [0 {:aromatic true}]}")]
    #[case::ring_membership_all(BondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)), "{:bond [0 {:ring-membership {:count 1}}]}")]
    #[case::ring_membership_size(BondConstraintAst::ring_membership(RingScope::Size(6), 1), "{:bond [0 {:ring-membership {:size 6 :count 1}}]}")]
    #[case::cis_trans_stereo_not_stereo(BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo), "{:bond [0 {:cis-trans-stereo :not-stereo}]}")]
    #[case::cis_trans_stereo_lit(BondConstraintAst::CisTransStereo(CisTransStereoAst::Stereo(StereoCosetAst::Lit(1))), "{:bond [0 {:cis-trans-stereo {:stereo 1}}]}")]
    fn test_bond_constraint_dsl_roundtrip(
        #[from(full_namespace)] namespace: MoleculeNamespace,
        #[case] input: BondConstraintAst,
        #[case] edn_source: &str,
    ) {
        let wrapped = Constraint::Bond(BondId(0), input.clone());
        let meta = MoleculeMetadata::default();
        let dsl = ConstraintDsl::from_ast(&wrapped, &meta).unwrap();
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = ConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&namespace).unwrap();
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
    fn test_constraints_dsl_empty_roundtrip(#[from(full_namespace)] namespace: MoleculeNamespace) {
        let meta = MoleculeMetadata::default();
        let cs = Constraints::new();
        let dsl = ConstraintsDsl::from_ast(&cs, &meta).unwrap();
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("[]").unwrap());
        let parsed = ConstraintsDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast(&namespace).unwrap();
        assert_eq!(back, cs);
    }

    #[rstest]
    fn test_constraints_dsl_multi_roundtrip(#[from(full_namespace)] namespace: MoleculeNamespace) {
        let meta = MoleculeMetadata::default();
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(
            AtomId(0),
            AtomConstraintAst::Valence(ValueAst::Lit(4)),
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
        let back = parsed.into_ast(&namespace).unwrap();
        assert_eq!(back, cs);
    }
}
