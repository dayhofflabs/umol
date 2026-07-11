//! Molecule DSL.
//!
//! `MoleculeDsl` wraps a `MoleculeAst` together with the `MoleculeMetadata` that records
//! the surface-form id/alias bindings (atom ids, bond ids, etc.). The EDN
//! form is a map keyed by `:atoms`, `:bonds`, `:dative-bonds`, `:aromatic-systems`,
//! `:multicenter-bonds`, `:noncovalent-bonds`, `:atom-aliases`/`:aliases`, and
//! `:constraints`. Each entity delegates to its own entity DSL. Constraints
//! parse directly into the typed `Constraint` tree.

// Closures like `|e| T::from_edn(e)` passed to `parse_vec` can't be replaced
// by bare `T::from_edn` — type-erasing the fn item loses the `for<'a>` HRTB
// on the `FromEdn<'a>` impl.
#![allow(clippy::redundant_closure)]

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::mem;
use std::str::FromStr;

use bimap::BiBTreeMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::{expand_bond_keyword, BondDsl};
use super::config::MoleculeDefaults;
use super::constraint::{read_constraints_dsl, ConstraintDsl, ConstraintsDsl};
use super::dative::DativeBondDsl;
use super::edn_utils::{
    atoms_pair, atoms_vec, eof_err, missing, optional_id, parse_vec, read_map, read_vec,
    required_key, two_atom_refs, unexpected_byte_kind,
};
use super::error::ParseError;
use super::multicenter::MulticenterBondDsl;
use super::namespace::{MoleculeNamespace, Namespace};
use super::noncovalent::NoncovalentBondDsl;
use super::refs::{
    parse_stereo_ligand, read_atom_ref, read_bond_ref, read_stereo_ligand, AtomRef, BondRef,
    StereoLigandRef,
};
use super::stereo::{
    expand_stereo_atom_keyword, expand_stereo_bond_keyword, StereoAtomDsl, StereoBondDsl,
};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::dative::DativeBondAst;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::{StereoLigand, StereoLigandKind};
use crate::ast::molecule::{MoleculeAst, MoleculeParts};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::NoncovalentBondAst;
use crate::ast::stereo::{StereoAtomAst, StereoBondAst};
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL for a whole molecule. Pairs `MoleculeAst` with `MoleculeMetadata`;
/// fields are private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    ast: MoleculeAst,
    metadata: MoleculeMetadata,
}

impl MoleculeDsl {
    pub fn from_parts(ast: MoleculeAst, metadata: MoleculeMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (MoleculeAst, MoleculeMetadata) {
        (self.ast, self.metadata)
    }
}

impl PartialEq for MoleculeDsl {
    fn eq(&self, other: &Self) -> bool {
        self.ast == other.ast && self.metadata == other.metadata
    }
}

impl Eq for MoleculeDsl {}

impl FromStr for MoleculeDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        MoleculeDsl::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for MoleculeDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let input = parse_molecule_input(edn)?;
        let (ast, namespace) = input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(
            ast,
            MoleculeMetadata::from(&namespace),
        ))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let mi = read_molecule_input(&mut de)?;
        de.expect_eof()?;
        let (ast, namespace) = mi.into_ast().map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(
            ast,
            MoleculeMetadata::from(&namespace),
        ))
    }
}

/// Surface-form metadata paired with a `MoleculeAst`. Records entity ids, atom aliases,
/// entity ids, `MoleculeDsl` keeps both fields private and rewraps atomically
/// through `from_parts`.
/// TODO: Review API, harmonize with `MoleculeNamespace`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeMetadata {
    atom_ids: IndexMap<AtomId, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondId, String>,
    dative_bond_ids: IndexMap<DativeBondId, String>,
    aromatic_system_ids: IndexMap<AromaticSystemId, String>,
    multicenter_bond_ids: IndexMap<MulticenterBondId, String>,
    noncovalent_bond_ids: IndexMap<NoncovalentBondId, String>,
    stereo_atom_ids: IndexMap<StereoAtomId, String>,
    stereo_bond_ids: IndexMap<StereoBondId, String>,
}

impl MoleculeMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether this metadata binds no `:id` keywords and no atom aliases — the shape an anonymous
    /// molecule (e.g. a sub-pattern) projects.
    pub fn is_empty(&self) -> bool {
        self.atom_ids.is_empty()
            && self.bond_ids.is_empty()
            && self.dative_bond_ids.is_empty()
            && self.aromatic_system_ids.is_empty()
            && self.multicenter_bond_ids.is_empty()
            && self.noncovalent_bond_ids.is_empty()
            && self.stereo_atom_ids.is_empty()
            && self.stereo_bond_ids.is_empty()
            && self.atom_aliases.is_empty()
    }

    pub fn atom_keyword(&self, id: AtomId) -> Option<&str> {
        self.atom_ids.get(&id).map(String::as_str)
    }

    pub fn bond_keyword(&self, id: BondId) -> Option<&str> {
        self.bond_ids.get(&id).map(String::as_str)
    }

    pub fn dative_bond_keyword(&self, id: DativeBondId) -> Option<&str> {
        self.dative_bond_ids.get(&id).map(String::as_str)
    }

    pub fn aromatic_system_keyword(&self, id: AromaticSystemId) -> Option<&str> {
        self.aromatic_system_ids.get(&id).map(String::as_str)
    }

    pub fn multicenter_bond_keyword(&self, id: MulticenterBondId) -> Option<&str> {
        self.multicenter_bond_ids.get(&id).map(String::as_str)
    }

    pub fn noncovalent_bond_keyword(&self, id: NoncovalentBondId) -> Option<&str> {
        self.noncovalent_bond_ids.get(&id).map(String::as_str)
    }

    pub fn stereo_atom_keyword(&self, id: StereoAtomId) -> Option<&str> {
        self.stereo_atom_ids.get(&id).map(String::as_str)
    }

    pub fn stereo_bond_keyword(&self, id: StereoBondId) -> Option<&str> {
        self.stereo_bond_ids.get(&id).map(String::as_str)
    }

    /// Whether `name` is already bound to any entity id (across all kinds). Linear scan.
    pub fn contains_id(&self, name: &str) -> bool {
        self.atom_ids.values().any(|n| n == name)
            || self.bond_ids.values().any(|n| n == name)
            || self.dative_bond_ids.values().any(|n| n == name)
            || self.aromatic_system_ids.values().any(|n| n == name)
            || self.multicenter_bond_ids.values().any(|n| n == name)
            || self.noncovalent_bond_ids.values().any(|n| n == name)
            || self.stereo_atom_ids.values().any(|n| n == name)
            || self.stereo_bond_ids.values().any(|n| n == name)
    }

    /// Name of the alias bound to this atom DSL, if any.
    pub fn atom_alias_for(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn has_atom_alias(&self, name: &str) -> bool {
        self.atom_aliases.contains_left(name)
    }

    pub fn has_atom_aliases(&self) -> bool {
        !self.atom_aliases.is_empty()
    }

    pub fn atom_aliases_len(&self) -> usize {
        self.atom_aliases.len()
    }

    pub fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_atom_keyword(&mut self, id: AtomId, name: impl Into<String>) {
        self.atom_ids.insert(id, name.into());
    }

    pub fn set_bond_keyword(&mut self, id: BondId, name: impl Into<String>) {
        self.bond_ids.insert(id, name.into());
    }

    pub fn set_dative_bond_keyword(&mut self, id: DativeBondId, name: impl Into<String>) {
        self.dative_bond_ids.insert(id, name.into());
    }

    pub fn set_aromatic_system_keyword(&mut self, id: AromaticSystemId, name: impl Into<String>) {
        self.aromatic_system_ids.insert(id, name.into());
    }

    pub fn set_multicenter_bond_keyword(&mut self, id: MulticenterBondId, name: impl Into<String>) {
        self.multicenter_bond_ids.insert(id, name.into());
    }

    pub fn set_noncovalent_bond_keyword(&mut self, id: NoncovalentBondId, name: impl Into<String>) {
        self.noncovalent_bond_ids.insert(id, name.into());
    }

    pub fn set_stereo_atom_keyword(&mut self, id: StereoAtomId, name: impl Into<String>) {
        self.stereo_atom_ids.insert(id, name.into());
    }

    pub fn set_stereo_bond_keyword(&mut self, id: StereoBondId, name: impl Into<String>) {
        self.stereo_bond_ids.insert(id, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases.insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_keyword(mut self, id: AtomId, name: impl Into<String>) -> Self {
        self.set_atom_keyword(id, name);
        self
    }

    pub fn with_bond_keyword(mut self, id: BondId, name: impl Into<String>) -> Self {
        self.set_bond_keyword(id, name);
        self
    }

    pub fn with_dative_bond_keyword(mut self, id: DativeBondId, name: impl Into<String>) -> Self {
        self.set_dative_bond_keyword(id, name);
        self
    }

    pub fn with_aromatic_system_keyword(
        mut self,
        id: AromaticSystemId,
        name: impl Into<String>,
    ) -> Self {
        self.set_aromatic_system_keyword(id, name);
        self
    }

    pub fn with_multicenter_bond_keyword(
        mut self,
        id: MulticenterBondId,
        name: impl Into<String>,
    ) -> Self {
        self.set_multicenter_bond_keyword(id, name);
        self
    }

    pub fn with_noncovalent_bond_keyword(
        mut self,
        id: NoncovalentBondId,
        name: impl Into<String>,
    ) -> Self {
        self.set_noncovalent_bond_keyword(id, name);
        self
    }

    pub fn with_stereo_atom_keyword(mut self, id: StereoAtomId, name: impl Into<String>) -> Self {
        self.set_stereo_atom_keyword(id, name);
        self
    }

    pub fn with_stereo_bond_keyword(mut self, id: StereoBondId, name: impl Into<String>) -> Self {
        self.set_stereo_bond_keyword(id, name);
        self
    }

    pub fn with_atom_alias(mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

/// The **render** query surface — `ref::denote` reads this to turn an AST id back into a surface
/// ref (id → keyword, else index). The mirror of `Namespace`: `MoleculeMetadata` implements it, and
/// later a reaction's `ReactionMetadata`. `denote` never emits the structural form, so `Metadata`
/// carries no participant index — it is a subset of the namespace surface.
pub trait Metadata {
    fn atom_keyword(&self, id: AtomId) -> Option<&str>;
    fn bond_keyword(&self, id: BondId) -> Option<&str>;
    fn dative_bond_keyword(&self, id: DativeBondId) -> Option<&str>;
    fn aromatic_system_keyword(&self, id: AromaticSystemId) -> Option<&str>;
    fn multicenter_bond_keyword(&self, id: MulticenterBondId) -> Option<&str>;
    fn noncovalent_bond_keyword(&self, id: NoncovalentBondId) -> Option<&str>;
    fn stereo_atom_keyword(&self, id: StereoAtomId) -> Option<&str>;
    fn stereo_bond_keyword(&self, id: StereoBondId) -> Option<&str>;
}

impl Metadata for MoleculeMetadata {
    fn atom_keyword(&self, id: AtomId) -> Option<&str> {
        self.atom_keyword(id)
    }
    fn bond_keyword(&self, id: BondId) -> Option<&str> {
        self.bond_keyword(id)
    }
    fn dative_bond_keyword(&self, id: DativeBondId) -> Option<&str> {
        self.dative_bond_keyword(id)
    }
    fn aromatic_system_keyword(&self, id: AromaticSystemId) -> Option<&str> {
        self.aromatic_system_keyword(id)
    }
    fn multicenter_bond_keyword(&self, id: MulticenterBondId) -> Option<&str> {
        self.multicenter_bond_keyword(id)
    }
    fn noncovalent_bond_keyword(&self, id: NoncovalentBondId) -> Option<&str> {
        self.noncovalent_bond_keyword(id)
    }
    fn stereo_atom_keyword(&self, id: StereoAtomId) -> Option<&str> {
        self.stereo_atom_keyword(id)
    }
    fn stereo_bond_keyword(&self, id: StereoBondId) -> Option<&str> {
        self.stereo_bond_keyword(id)
    }
}

impl From<&MoleculeNamespace> for MoleculeMetadata {
    /// Project the namespace to its roundtrip subset: the eight `id → keyword` maps (the inverse of
    /// the namespace's `find_by_keyword`) plus the atom aliases. The namespace is the source of truth;
    /// this is the derived view — parse-only data (participant indexes, counts) is dropped.
    fn from(namespace: &MoleculeNamespace) -> Self {
        let mut metadata = MoleculeMetadata::new();
        for (id, keyword) in namespace.atom_keywords() {
            metadata.set_atom_keyword(id, keyword);
        }
        for (id, keyword) in namespace.bond_keywords() {
            metadata.set_bond_keyword(id, keyword);
        }
        for (id, keyword) in namespace.dative_bond_keywords() {
            metadata.set_dative_bond_keyword(id, keyword);
        }
        for (id, keyword) in namespace.aromatic_system_keywords() {
            metadata.set_aromatic_system_keyword(id, keyword);
        }
        for (id, keyword) in namespace.multicenter_bond_keywords() {
            metadata.set_multicenter_bond_keyword(id, keyword);
        }
        for (id, keyword) in namespace.noncovalent_bond_keywords() {
            metadata.set_noncovalent_bond_keyword(id, keyword);
        }
        for (id, keyword) in namespace.stereo_atom_keywords() {
            metadata.set_stereo_atom_keyword(id, keyword);
        }
        for (id, keyword) in namespace.stereo_bond_keywords() {
            metadata.set_stereo_bond_keyword(id, keyword);
        }
        for (name, dsl) in namespace.atom_aliases() {
            metadata.add_atom_alias(name, dsl.clone());
        }
        metadata
    }
}

/// Direct EDN parsing for `MoleculeAst`. Accepts the same molecule-map
/// surface as [`MoleculeDsl::from_edn`]; any id keywords or aliases in the
/// input resolve to positional indices, then the metadata is discarded —
/// the result is metadata-free.
impl<'de> FromEdn<'de> for MoleculeAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        MoleculeDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        MoleculeDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for MoleculeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for MoleculeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `MoleculeAst`. Always emits canonical positional
/// refs (no id keywords, no aliases) since the AST carries no metadata.
/// For id-bearing surface output, wrap in [`MoleculeDsl`] with appropriate
/// [`MoleculeMetadata`] and call [`MoleculeDsl::to_edn`].
impl ToEdn for MoleculeAst {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(self, &MoleculeMetadata::default())
    }
}

// Streaming parse of the molecule map.
pub(super) fn read_molecule_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MoleculeInput, EdnError> {
    de.consume_byte(b'{')?;
    let mut mi = MoleculeInput::default();
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "atoms" => mi.atoms = read_vec(de, read_atom_entry)?,
            "bonds" => mi.bonds = read_vec(de, read_bond_entry)?,
            "dative-bonds" => mi.dative_bonds = read_vec(de, read_dative_bond_entry)?,
            "aromatic-systems" => mi.aromatic_systems = read_vec(de, read_aromatic_system_entry)?,
            "multicenter-bonds" => {
                mi.multicenter_bonds = read_vec(de, read_multicenter_bond_entry)?
            }
            "noncovalent-bonds" => {
                mi.noncovalent_bonds = read_vec(de, read_noncovalent_bond_entry)?
            }
            "stereo-atoms" => mi.stereo_atoms = read_vec(de, read_stereo_atom_entry)?,
            "stereo-bonds" => mi.stereo_bonds = read_vec(de, read_stereo_bond_entry)?,
            "atom-aliases" => mi.atom_aliases = read_atom_aliases(de)?,
            "constraints" => mi.constraints = read_constraints_dsl(de)?,
            "guards" => {
                de.read_skip_value()?;
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule".into()],
                }
                .into());
            }
        }
    }
    Ok(mi)
}

pub(super) fn read_atom_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AtomEntryInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'[' => {
            de.consume_byte(b'[')?;
            let id = de.read_keyword_name()?.into_owned();
            let spec = read_atom_spec(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryInput { id: Some(id), spec })
        }
        _ => Ok(AtomEntryInput {
            id: None,
            spec: read_atom_spec(de)?,
        }),
    }
}

fn read_atom_spec(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomSpecInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'"' => {
            let s = de.read_string()?;
            let dsl: AtomDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomSpecInput::Bare(Box::new(dsl)))
        }
        b':' => {
            let name = de.read_keyword_name()?;
            Ok(AtomSpecInput::Alias(name.into_owned()))
        }
        b => Err(DeError::TypeMismatch {
            expected: "atom-string or :alias",
            got: unexpected_byte_kind(b),
            path: vec!["atom-spec".into()],
        }
        .into()),
    }
}

fn read_bond_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<BondDsl, EdnError> {
    let byte = de.peek_byte()?.ok_or_else(eof_err)?;
    let text: Cow<'_, str> = match byte {
        b':' => {
            let name = de.read_keyword_name()?;
            let expanded = expand_bond_keyword(name.as_ref())
                .ok_or_else(|| DeError::Custom(format!("unknown bond keyword :{}", name)))?;
            Cow::Borrowed(expanded)
        }
        _ => de.read_string()?,
    };
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("bond", e).into())
}

fn read_dative_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<DativeBondDsl, EdnError> {
    let byte = de.peek_byte()?.ok_or_else(eof_err)?;
    let text: Cow<'_, str> = match byte {
        b':' => {
            let name = de.read_keyword_name()?;
            let expanded = super::dative::expand_dative_keyword(name.as_ref())
                .ok_or_else(|| DeError::Custom(format!("unknown dative keyword :{}", name)))?;
            Cow::Borrowed(expanded)
        }
        _ => de.read_string()?,
    };
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("dative", e).into())
}

/// A two-endpoint `:atoms` vector for a binary relation: exactly two refs.
pub(super) fn read_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<BondEntryInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'[' => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_dsl(de)?;
            de.consume_byte(b']')?;
            Ok(BondEntryInput {
                id: None,
                first: a,
                second: b,
                bond,
            })
        }
        b'{' => {
            let mut id = None;
            let mut atoms = None;
            let mut bond = None;
            read_map(de, |de, key| {
                match key {
                    "id" => id = Some(de.read_keyword_name()?.into_owned()),
                    "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
                    "type" => bond = Some(read_bond_dsl(de)?),
                    _ => de.read_skip_value()?,
                }
                Ok(())
            })?;
            let atoms = atoms.ok_or_else(|| missing("atoms", "bond-entry"))?;
            let [a, b] = two_atom_refs(atoms, "bond-entry")?;
            Ok(BondEntryInput {
                id,
                first: a,
                second: b,
                bond: bond.ok_or_else(|| missing("type", "bond-entry"))?,
            })
        }
        bb => Err(DeError::TypeMismatch {
            expected: "bond-entry map or 3-vec",
            got: unexpected_byte_kind(bb),
            path: vec!["bond-entry".into()],
        }
        .into()),
    }
}

pub(super) fn read_dative_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DativeBondEntryInput, EdnError> {
    let mut id = None;
    let mut donors = None;
    let mut acceptor = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "donors" => donors = Some(read_vec(de, read_atom_ref)?),
            "acceptor" => acceptor = Some(read_atom_ref(de)?),
            "type" => bond = Some(read_dative_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(DativeBondEntryInput {
        id,
        donors: donors.ok_or_else(|| missing("donors", "dative-bond-entry"))?,
        acceptor: acceptor.ok_or_else(|| missing("acceptor", "dative-bond-entry"))?,
        bond: bond.ok_or_else(|| missing("type", "dative-bond-entry"))?,
    })
}

pub(super) fn read_aromatic_system_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystemEntryInput, EdnError> {
    let mut id = None;
    let mut atoms = None;
    let mut system = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let s = de.read_string()?;
                system = Some(
                    s.as_ref()
                        .parse::<AromaticSystemDsl>()
                        .map_err(|e| DeError::subgrammar("aromatic", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let system = system.ok_or_else(|| missing("type", "aromatic-system-entry"))?;
    Ok(AromaticSystemEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "aromatic-system-entry"))?,
        system,
    })
}

pub(super) fn read_multicenter_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBondEntryInput, EdnError> {
    let mut id = None;
    let mut atoms = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let s = de.read_string()?;
                bond = Some(
                    s.as_ref()
                        .parse::<MulticenterBondDsl>()
                        .map_err(|e| DeError::subgrammar("multicenter", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let bond = bond.ok_or_else(|| missing("type", "multicenter-bond-entry"))?;
    Ok(MulticenterBondEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "multicenter-bond-entry"))?,
        bond,
    })
}

pub(super) fn read_noncovalent_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondEntryInput, EdnError> {
    let mut id = None;
    let mut atoms = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let text = de.read_string_or_keyword()?;
                bond = Some(
                    text.as_ref()
                        .parse::<NoncovalentBondDsl>()
                        .map_err(|e| DeError::subgrammar("noncovalent", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let atoms = atoms.ok_or_else(|| missing("atoms", "noncovalent-bond-entry"))?;
    let [a, b] = two_atom_refs(atoms, "noncovalent-bond-entry")?;
    Ok(NoncovalentBondEntryInput {
        id,
        first: a,
        second: b,
        bond: bond.ok_or_else(|| missing("type", "noncovalent-bond-entry"))?,
    })
}

fn read_stereo_atom_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<StereoAtomDsl, EdnError> {
    if de.peek_byte()?.ok_or_else(eof_err)? == b':' {
        let kw = de.read_keyword_name()?;
        let expanded = expand_stereo_atom_keyword(kw.as_ref())
            .ok_or_else(|| DeError::Custom(format!("unknown stereo atom keyword :{kw}")))?;
        expanded
            .parse::<StereoAtomDsl>()
            .map_err(|e| DeError::subgrammar("stereo atom", e).into())
    } else {
        de.read_string()?
            .as_ref()
            .parse::<StereoAtomDsl>()
            .map_err(|e| DeError::subgrammar("stereo atom", e).into())
    }
}

fn read_stereo_bond_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<StereoBondDsl, EdnError> {
    if de.peek_byte()?.ok_or_else(eof_err)? == b':' {
        let kw = de.read_keyword_name()?;
        let expanded = expand_stereo_bond_keyword(kw.as_ref())
            .ok_or_else(|| DeError::Custom(format!("unknown stereo bond keyword :{kw}")))?;
        expanded
            .parse::<StereoBondDsl>()
            .map_err(|e| DeError::subgrammar("stereo bond", e).into())
    } else {
        de.read_string()?
            .as_ref()
            .parse::<StereoBondDsl>()
            .map_err(|e| DeError::subgrammar("stereo bond", e).into())
    }
}

pub(super) fn read_stereo_atom_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoAtomEntryInput, EdnError> {
    let mut id = None;
    let mut site = None;
    let mut ligands = None;
    let mut stereo = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "site" => site = Some(read_atom_ref(de)?),
            "ligands" => ligands = Some(read_vec(de, read_stereo_ligand)?),
            "type" => stereo = Some(read_stereo_atom_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(StereoAtomEntryInput {
        id,
        site: site.ok_or_else(|| missing("site", "stereo-atom-entry"))?,
        ligands: ligands.ok_or_else(|| missing("ligands", "stereo-atom-entry"))?,
        stereo: stereo.ok_or_else(|| missing("type", "stereo-atom-entry"))?,
    })
}

pub(super) fn read_stereo_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoBondEntryInput, EdnError> {
    let mut id = None;
    let mut site = None;
    let mut ligands = None;
    let mut stereo = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "site" => site = Some(read_bond_ref(de)?),
            "ligands" => ligands = Some(read_vec(de, read_stereo_ligand)?),
            "type" => stereo = Some(read_stereo_bond_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(StereoBondEntryInput {
        id,
        site: site.ok_or_else(|| missing("site", "stereo-bond-entry"))?,
        ligands: ligands.ok_or_else(|| missing("ligands", "stereo-bond-entry"))?,
        stereo: stereo.ok_or_else(|| missing("type", "stereo-bond-entry"))?,
    })
}

pub(super) fn read_atom_aliases(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<Vec<(String, Box<AtomDsl>)>, EdnError> {
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            break;
        }
        let name = de.read_keyword_name()?.into_owned();
        if de.try_consume_byte(b']')? {
            return Err(DeError::Custom(
                ":atom-aliases must have even length (keyword/atom-string pairs)".into(),
            )
            .into());
        }
        let s = de.read_string()?;
        let dsl: AtomDsl = s
            .as_ref()
            .parse()
            .map_err(|e| DeError::subgrammar("atom", e))?;
        out.push((name, Box::new(dsl)));
    }
    Ok(out)
}

impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(&self.ast, &self.metadata)
    }
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    type Ctx = MoleculeDefaults;

    fn from_ast(ast: &MoleculeAst, cfg: &Self::Ctx) -> Self {
        let mut ast_out = ast.clone();
        for atom in ast_out.atoms_mut() {
            *atom = AtomDsl::from_ast(atom, &cfg.atom).0;
        }
        for bond in ast_out.bonds_mut() {
            *bond = BondDsl::from_ast(bond, &cfg.bond).0;
        }
        for system in ast_out.aromatic_systems_mut() {
            *system = AromaticSystemDsl::from_ast(system, &cfg.aromatic_system).0;
        }
        for bond in ast_out.multicenter_bonds_mut() {
            *bond = MulticenterBondDsl::from_ast(bond, &cfg.multicenter_bond).0;
        }
        for bond in ast_out.dative_bonds_mut() {
            *bond = DativeBondDsl::from_ast(bond, &cfg.dative_bond).0;
        }
        for bond in ast_out.noncovalent_bonds_mut() {
            *bond = NoncovalentBondDsl::from_ast(bond, &cfg.noncovalent_bond).0;
        }
        for stereo_atom in ast_out.stereo_atoms_mut() {
            *stereo_atom = StereoAtomDsl::from_ast(stereo_atom, &cfg.stereo_atom).0;
        }
        for stereo_bond in ast_out.stereo_bonds_mut() {
            *stereo_bond = StereoBondDsl::from_ast(stereo_bond, &cfg.stereo_bond).0;
        }
        MoleculeDsl {
            ast: ast_out,
            metadata: MoleculeMetadata::default(),
        }
    }
}

impl IntoAst<MoleculeAst> for MoleculeDsl {
    type Ctx = MoleculeDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> MoleculeAst {
        let mut ast = self.ast;
        for atom in ast.atoms_mut() {
            *atom = AtomDsl(mem::take(atom)).into_ast(&cfg.atom);
        }
        for bond in ast.bonds_mut() {
            *bond = BondDsl(mem::take(bond)).into_ast(&cfg.bond);
        }
        for bond in ast.dative_bonds_mut() {
            *bond = DativeBondDsl(mem::take(bond)).into_ast(&cfg.dative_bond);
        }
        for system in ast.aromatic_systems_mut() {
            *system = AromaticSystemDsl(mem::take(system)).into_ast(&cfg.aromatic_system);
        }
        for bond in ast.multicenter_bonds_mut() {
            *bond = MulticenterBondDsl(mem::take(bond)).into_ast(&cfg.multicenter_bond);
        }
        for bond in ast.noncovalent_bonds_mut() {
            *bond = NoncovalentBondDsl(mem::take(bond)).into_ast(&cfg.noncovalent_bond);
        }
        for stereo_atom in ast.stereo_atoms_mut() {
            *stereo_atom = StereoAtomDsl(mem::take(stereo_atom)).into_ast(&cfg.stereo_atom);
        }
        for stereo_bond in ast.stereo_bonds_mut() {
            *stereo_bond = StereoBondDsl(mem::take(stereo_bond)).into_ast(&cfg.stereo_bond);
        }
        ast
    }
}

pub(super) fn render_molecule_edn(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(8);
    map.insert(Edn::keyword("atoms"), render_atoms(ast, meta));
    map.insert(Edn::keyword("bonds"), render_bonds(ast, meta));
    if ast.dative_bonds().count() > 0 {
        map.insert(Edn::keyword("dative-bonds"), render_dative(ast, meta));
    }
    if ast.aromatic_systems().count() > 0 {
        map.insert(Edn::keyword("aromatic-systems"), render_aromatic(ast, meta));
    }
    if ast.multicenter_bonds().count() > 0 {
        map.insert(
            Edn::keyword("multicenter-bonds"),
            render_multicenter(ast, meta),
        );
    }
    if ast.noncovalent_bonds().count() > 0 {
        map.insert(
            Edn::keyword("noncovalent-bonds"),
            render_noncovalent(ast, meta),
        );
    }
    if ast.stereo_atoms().count() > 0 {
        map.insert(Edn::keyword("stereo-atoms"), render_stereo_atoms(ast, meta));
    }
    if ast.stereo_bonds().count() > 0 {
        map.insert(Edn::keyword("stereo-bonds"), render_stereo_bonds(ast, meta));
    }
    if meta.has_atom_aliases() {
        map.insert(Edn::keyword("atom-aliases"), render_atom_aliases(meta));
    }
    let constraints_dsl = ConstraintsDsl::from_ast(ast.constraints(), meta)
        .expect("ConstraintsDsl::from_ast is infallible for a well-formed AST");
    if !constraints_dsl.0.is_empty() {
        map.insert(Edn::keyword("constraints"), constraints_dsl.to_edn());
    }
    Edn::Map(map)
}

fn render_atoms(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .atoms()
        .iter()
        .map(|view| render_atom_entry(view.id, view.ast, meta))
        .collect();
    Edn::Vector(entries.into())
}

/// An atom value: its alias keyword if one is bound, else the atom-string.
pub(super) fn render_atom_value(atom: &AtomAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let dsl = AtomDsl::from_ref(atom);
    match meta.atom_alias_for(dsl) {
        Some(alias) => Edn::Keyword(EdnKeyword::owned(alias.to_string())),
        None => dsl.to_edn(),
    }
}

fn render_atom_entry(id: AtomId, atom: &AtomAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let spec = render_atom_value(atom, meta);
    match meta.atom_keyword(id) {
        Some(id) => Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(id.to_string())), spec].into()),
        None => spec,
    }
}

fn render_atom_ref(id: AtomId, meta: &MoleculeMetadata) -> Edn<'static> {
    match meta.atom_keyword(id) {
        Some(id) => Edn::Keyword(EdnKeyword::owned(id.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A bond entry — `[a b type]`, or `{:id … :atoms [a b] :type type}` when the bond has an id.
/// `type_edn` is the already-rendered `:type` Edn — one bond-dsl for a molecule, or a `[left right]`
/// vector / op-wrapped map for a span entry; it is not an ast.
pub(super) fn render_bond_entry(
    id: BondId,
    [a, b]: [AtomId; 2],
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let first = render_atom_ref(a, meta);
    let second = render_atom_ref(b, meta);
    match meta.bond_keyword(id) {
        Some(name) => {
            let mut m = EdnMap::with_capacity(3);
            m.insert(
                Edn::keyword("id"),
                Edn::Keyword(EdnKeyword::owned(name.to_string())),
            );
            m.insert(
                Edn::keyword("atoms"),
                Edn::Vector(vec![first, second].into()),
            );
            m.insert(Edn::keyword("type"), type_edn);
            Edn::Map(m)
        }
        None => Edn::Vector(vec![first, second, type_edn].into()),
    }
}

fn render_bonds(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .map(|view| {
            render_bond_entry(
                view.id,
                view.atom_ids(),
                BondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

// Overlay entries: `render_<entity>_entry` builds one entry map (`:id`? + participants + `:type`),
// with `:type` = the caller-supplied `type_edn`. `render_<entity>` passes the ast's rendered type;
// the span renderers (reaction_span) pass a `{:add|:modify|:remove}`-wrapped type over the same entry.

pub(super) fn render_dative_entry(
    id: DativeBondId,
    donors: impl Iterator<Item = AtomId>,
    acceptor: AtomId,
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = meta.dative_bond_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(
        Edn::keyword("donors"),
        Edn::Vector(
            donors
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("acceptor"), render_atom_ref(acceptor, meta));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_dative(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .dative_bonds()
        .iter()
        .map(|view| {
            render_dative_entry(
                view.id,
                view.donor_ids(),
                view.acceptor_id(),
                DativeBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_aromatic_entry(
    id: AromaticSystemId,
    atoms: impl Iterator<Item = AtomId>,
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(id) = meta.aromatic_system_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(
            atoms
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_aromatic(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .aromatic_systems()
        .iter()
        .map(|view| {
            render_aromatic_entry(
                view.id,
                view.atom_ids(),
                Edn::Str(Cow::Owned(
                    AromaticSystemDsl::from_ref(view.ast).to_string(),
                )),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_multicenter_entry(
    id: MulticenterBondId,
    atoms: impl Iterator<Item = AtomId>,
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(id) = meta.multicenter_bond_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(
            atoms
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_multicenter(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .multicenter_bonds()
        .iter()
        .map(|view| {
            render_multicenter_entry(
                view.id,
                view.atom_ids(),
                Edn::Str(Cow::Owned(
                    MulticenterBondDsl::from_ref(view.ast).to_string(),
                )),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_noncovalent_entry(
    id: NoncovalentBondId,
    [a, b]: [AtomId; 2],
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(id) = meta.noncovalent_bond_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(vec![render_atom_ref(a, meta), render_atom_ref(b, meta)].into()),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_noncovalent(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .map(|view| {
            render_noncovalent_entry(
                view.id,
                view.atom_ids(),
                NoncovalentBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_atom_entry(
    id: StereoAtomId,
    site: AtomId,
    ligands: Vec<Edn<'static>>,
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = meta.stereo_atom_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(Edn::keyword("site"), render_atom_ref(site, meta));
    m.insert(Edn::keyword("ligands"), Edn::Vector(ligands.into()));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_stereo_atoms(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .stereo_atoms()
        .iter()
        .map(|view| {
            render_stereo_atom_entry(
                view.id,
                view.site_id(),
                view.ligand_frame()
                    .into_iter()
                    .map(|l| render_stereo_ligand(l, meta))
                    .collect(),
                StereoAtomDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_bond_entry(
    id: StereoBondId,
    site: BondId,
    ligands: Vec<Edn<'static>>,
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = meta.stereo_bond_keyword(id) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.to_string())),
        );
    }
    m.insert(Edn::keyword("site"), render_bond_ref(site, meta));
    m.insert(Edn::keyword("ligands"), Edn::Vector(ligands.into()));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_stereo_bonds(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .stereo_bonds()
        .iter()
        .map(|view| {
            render_stereo_bond_entry(
                view.id,
                view.site_id(),
                view.ligand_frame()
                    .into_iter()
                    .map(|l| render_stereo_ligand(l, meta))
                    .collect(),
                StereoBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_ligand(ligand: StereoLigand, meta: &MoleculeMetadata) -> Edn<'static> {
    let atom = render_atom_ref(ligand.atom_id, meta);
    match ligand.kind {
        StereoLigandKind::Atom => atom,
        StereoLigandKind::ImplicitHydrogen => Edn::Vector(vec![Edn::keyword("h"), atom].into()),
        StereoLigandKind::LonePair => Edn::Vector(vec![Edn::keyword("lp"), atom].into()),
    }
}

fn render_bond_ref(id: BondId, meta: &MoleculeMetadata) -> Edn<'static> {
    match meta.bond_keyword(id) {
        Some(id) => Edn::Keyword(EdnKeyword::owned(id.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

fn render_atom_aliases(meta: &MoleculeMetadata) -> Edn<'static> {
    let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases_len() * 2);
    for (name, dsl) in meta.iter_atom_aliases() {
        pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
        pairs.push(dsl.to_edn());
    }
    Edn::Vector(pairs.into())
}

// Unresolved, owned-by-value tree that mirrors the EDN shape. Atom entries and
// per-bond endpoints carry `AtomRef` (index or id); constraint leaves carry
// typed per-entity `Constraint*` variants already parsed from their EDN form.
// Lowered destructively via `into_ast(self, cfg)` so that allocations move
// into the final `MoleculeAst`.

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct MoleculeInput {
    pub(crate) atoms: Vec<AtomEntryInput>,
    pub(crate) bonds: Vec<BondEntryInput>,
    pub(crate) dative_bonds: Vec<DativeBondEntryInput>,
    pub(crate) aromatic_systems: Vec<AromaticSystemEntryInput>,
    pub(crate) multicenter_bonds: Vec<MulticenterBondEntryInput>,
    pub(crate) noncovalent_bonds: Vec<NoncovalentBondEntryInput>,
    pub(crate) stereo_atoms: Vec<StereoAtomEntryInput>,
    pub(crate) stereo_bonds: Vec<StereoBondEntryInput>,
    pub(crate) atom_aliases: Vec<(String, Box<AtomDsl>)>,
    pub(crate) constraints: Vec<ConstraintDsl>,
}

/// Atom entry in a parsed molecule map. Mirrors the DSL spec §4 grammar
/// `atom-entry ::= atom-spec | [ keyword atom-spec ]`.
/// TODO: Fix pub(crate) visibility markers on the struct fields.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AtomEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) spec: AtomSpecInput,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum AtomSpecInput {
    Bare(Box<AtomDsl>),
    Alias(String),
}

/// Resolve an atom spec to its `AtomAst`: a bare value is its own atom; an alias is looked up in the
/// table (unknown → error). Shared by the molecule, reaction, and span `into_ast` paths.
pub(super) fn resolve_atom_spec(
    spec: AtomSpecInput,
    namespace: &impl Namespace,
) -> Result<AtomAst, ParseError> {
    match spec {
        AtomSpecInput::Bare(dsl) => Ok(dsl.0),
        AtomSpecInput::Alias(name) => match namespace.find_atom_alias(&name) {
            Some(dsl) => Ok(dsl.0.clone()),
            None => Err(ParseError::InvalidValue(format!(
                "unknown atom alias :{name}"
            ))),
        },
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) first: AtomRef,
    pub(crate) second: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) donors: Vec<AtomRef>,
    pub(crate) acceptor: AtomRef,
    pub(crate) bond: DativeBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AromaticSystemEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) system: AromaticSystemDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MulticenterBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) bond: MulticenterBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct NoncovalentBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) first: AtomRef,
    pub(crate) second: AtomRef,
    pub(crate) bond: NoncovalentBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct StereoAtomEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) site: AtomRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
    pub(crate) stereo: StereoAtomDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct StereoBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) site: BondRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
    pub(crate) stereo: StereoBondDsl,
}

impl MoleculeInput {
    /// Destructive lowering: consumes the input, resolves refs against the
    /// built id scopes, and produces the final `MoleculeAst` with its
    /// `MoleculeMetadata`. Called from `FromEdn::from_edn` and the streaming path.
    pub(crate) fn into_ast(self) -> Result<(MoleculeAst, MoleculeNamespace), ParseError> {
        let MoleculeInput {
            atoms: atom_entries,
            bonds: bond_entries,
            dative_bonds: dative_entries,
            aromatic_systems: aromatic_entries,
            multicenter_bonds: multicenter_entries,
            noncovalent_bonds: noncovalent_entries,
            stereo_atoms: stereo_atom_entries,
            stereo_bonds: stereo_bond_entries,
            atom_aliases: alias_entries,
            constraints: constraint_dsls,
        } = self;

        // Register atoms (positions + keywords), then the bijective aliases. `register_*` enforces
        // id-disjointness (across every entity kind + aliases) and alias bijectivity as it goes, so
        // there are no side tables; atom specs resolve against the namespace once it is complete.
        let mut namespace = MoleculeNamespace::default();
        for entry in &atom_entries {
            namespace.register_atom(entry.id.clone())?;
        }
        for (name, dsl) in alias_entries {
            namespace.register_atom_alias(name, dsl)?;
        }

        // Resolve atom aliases.
        let atoms: Vec<AtomAst> = atom_entries
            .into_iter()
            .map(|entry| resolve_atom_spec(entry.spec, &namespace))
            .collect::<Result<_, _>>()?;

        // Bonds.
        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::with_capacity(bond_entries.len());
        for entry in bond_entries {
            let a = entry.first.resolve(&namespace)?;
            let b = entry.second.resolve(&namespace)?;
            namespace.register_bond(entry.id, a, b)?;
            bonds.push((a, b, entry.bond.0));
        }

        // Dative bonds.
        let mut dative_list: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> =
            Vec::with_capacity(dative_entries.len());
        for entry in dative_entries {
            let donors = entry
                .donors
                .into_iter()
                .map(|d| d.resolve(&namespace))
                .collect::<Result<Vec<_>, _>>()?;
            if donors.is_empty() {
                return Err(ParseError::InvalidValue(
                    "dative bond requires at least one donor".to_string(),
                ));
            }
            let acceptor = entry.acceptor.resolve(&namespace)?;
            namespace.register_dative_bond(entry.id, &donors, acceptor)?;
            dative_list.push((donors, acceptor, entry.bond.0));
        }

        // Aromatic systems.
        let mut aromatic_list: Vec<(Vec<AtomId>, AromaticSystemAst)> =
            Vec::with_capacity(aromatic_entries.len());
        for entry in aromatic_entries {
            let atoms_resolved: Vec<AtomId> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(&namespace))
                .collect::<Result<_, _>>()?;
            namespace.register_aromatic_system(entry.id, &atoms_resolved)?;
            aromatic_list.push((atoms_resolved, entry.system.0));
        }

        // Multicenter bonds.
        let mut multicenter_list: Vec<(Vec<AtomId>, MulticenterBondAst)> =
            Vec::with_capacity(multicenter_entries.len());
        for entry in multicenter_entries {
            let atoms_resolved: Vec<AtomId> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(&namespace))
                .collect::<Result<_, _>>()?;
            namespace.register_multicenter_bond(entry.id, &atoms_resolved)?;
            multicenter_list.push((atoms_resolved, entry.bond.0));
        }

        // Noncovalent bonds.
        let mut noncovalent_list: Vec<(AtomId, AtomId, NoncovalentBondAst)> =
            Vec::with_capacity(noncovalent_entries.len());
        for entry in noncovalent_entries {
            let first = entry.first.resolve(&namespace)?;
            let second = entry.second.resolve(&namespace)?;
            namespace.register_noncovalent_bond(entry.id, first, second)?;
            noncovalent_list.push((first, second, entry.bond.0));
        }

        // Stereo atoms.
        let mut stereo_atom_list: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)> =
            Vec::with_capacity(stereo_atom_entries.len());
        for entry in stereo_atom_entries {
            let site = entry.site.resolve(&namespace)?;
            let ligands: Vec<StereoLigand> = entry
                .ligands
                .into_iter()
                .map(|l| Ok(StereoLigand::new(l.atom.resolve(&namespace)?, l.kind)))
                .collect::<Result<_, ParseError>>()?;
            namespace.register_stereo_atom(entry.id, site, &ligands)?;
            stereo_atom_list.push((site, ligands, entry.stereo.0));
        }

        // Stereo bonds.
        let mut stereo_bond_list: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> =
            Vec::with_capacity(stereo_bond_entries.len());
        for entry in stereo_bond_entries {
            let site = entry.site.resolve(&namespace)?;
            let ligands: Vec<StereoLigand> = entry
                .ligands
                .into_iter()
                .map(|l| Ok(StereoLigand::new(l.atom.resolve(&namespace)?, l.kind)))
                .collect::<Result<_, ParseError>>()?;
            namespace.register_stereo_bond(entry.id, site, &ligands)?;
            stereo_bond_list.push((site, ligands, entry.stereo.0));
        }

        // The namespace is complete; constraints resolve against it directly. `MoleculeMetadata` is
        // projected only at the DSL boundary, not here.
        let constraints = ConstraintsDsl(constraint_dsls).into_ast(&namespace)?;

        let ast = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            dative: dative_list,
            aromatic: aromatic_list,
            multicenter: multicenter_list,
            noncovalent: noncovalent_list,
            stereo_atoms: stereo_atom_list,
            stereo_bonds: stereo_bond_list,
            constraints,
        });
        Ok((ast, namespace))
    }
}

pub(super) fn parse_molecule_input(edn: &Edn<'_>) -> Result<MoleculeInput, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "molecule map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut input = MoleculeInput::default();
    for (k, v) in m.iter() {
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["molecule".into()],
            });
        };
        match key.name() {
            "atoms" => input.atoms = parse_vec(v, ":atoms", parse_atom_entry)?,
            "bonds" => input.bonds = parse_vec(v, ":bonds", parse_bond_entry)?,
            "dative-bonds" => {
                input.dative_bonds = parse_vec(v, ":dative-bonds", parse_dative_bond_entry)?
            }
            "aromatic-systems" => {
                input.aromatic_systems =
                    parse_vec(v, ":aromatic-systems", parse_aromatic_system_entry)?
            }
            "multicenter-bonds" => {
                input.multicenter_bonds =
                    parse_vec(v, ":multicenter-bonds", parse_multicenter_bond_entry)?
            }
            "noncovalent-bonds" => {
                input.noncovalent_bonds =
                    parse_vec(v, ":noncovalent-bonds", parse_noncovalent_bond_entry)?
            }
            "stereo-atoms" => {
                input.stereo_atoms = parse_vec(v, ":stereo-atoms", parse_stereo_atom_entry)?
            }
            "stereo-bonds" => {
                input.stereo_bonds = parse_vec(v, ":stereo-bonds", parse_stereo_bond_entry)?
            }
            "atom-aliases" => input.atom_aliases = parse_atom_aliases(v)?,
            "constraints" => {
                input.constraints = parse_vec(v, ":constraints", |e| ConstraintDsl::from_edn(e))?
            }
            "guards" => {
                // Spec §4 lists :guards as a future-reserved key; ignore for now.
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule".into()],
                });
            }
        }
    }
    Ok(input)
}

pub(super) fn parse_atom_entry(edn: &Edn<'_>) -> Result<AtomEntryInput, DeError> {
    match edn {
        Edn::Str(s) => {
            let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomEntryInput {
                id: None,
                spec: AtomSpecInput::Bare(Box::new(dsl)),
            })
        }
        Edn::Keyword(k) => Ok(AtomEntryInput {
            id: None,
            spec: AtomSpecInput::Alias(k.name().to_string()),
        }),
        Edn::Vector(v) if v.len() == 2 => {
            let Edn::Keyword(id_kw) = &v[0] else {
                return Err(DeError::TypeMismatch {
                    expected: "keyword id",
                    got: v[0].kind(),
                    path: vec!["atom-entry".into()],
                });
            };
            let spec = parse_atom_spec(&v[1])?;
            Ok(AtomEntryInput {
                id: Some(id_kw.name().to_string()),
                spec,
            })
        }
        other => Err(DeError::TypeMismatch {
            expected: "atom-string / keyword / [keyword atom-spec]",
            got: other.kind(),
            path: vec!["atom-entry".into()],
        }),
    }
}

fn parse_atom_spec(edn: &Edn<'_>) -> Result<AtomSpecInput, DeError> {
    match edn {
        Edn::Str(s) => {
            let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomSpecInput::Bare(Box::new(dsl)))
        }
        Edn::Keyword(k) => Ok(AtomSpecInput::Alias(k.name().to_string())),
        other => Err(DeError::TypeMismatch {
            expected: "atom-string or keyword alias",
            got: other.kind(),
            path: vec!["atom-spec".into()],
        }),
    }
}

pub(super) fn parse_bond_entry(edn: &Edn<'_>) -> Result<BondEntryInput, DeError> {
    match edn {
        Edn::Vector(v) if v.len() == 3 => Ok(BondEntryInput {
            id: None,
            first: AtomRef::from_edn(&v[0])?,
            second: AtomRef::from_edn(&v[1])?,
            bond: BondDsl::from_edn(&v[2])?,
        }),
        Edn::Map(m) => {
            let [a, b] = atoms_pair(m, "bond-entry")?;
            Ok(BondEntryInput {
                id: optional_id(m)?,
                first: a,
                second: b,
                bond: BondDsl::from_edn(required_key(m, "type", "bond-entry")?)?,
            })
        }
        other => Err(DeError::TypeMismatch {
            expected: "bond-entry map or 3-vec",
            got: other.kind(),
            path: vec!["bond-entry".into()],
        }),
    }
}

pub(super) fn parse_dative_bond_entry(edn: &Edn<'_>) -> Result<DativeBondEntryInput, DeError> {
    let m = expect_map(edn, "dative-bond-entry")?;
    let donors = parse_vec(
        required_key(m, "donors", "dative-bond-entry")?,
        ":donors",
        |e| AtomRef::from_edn(e),
    )?;
    Ok(DativeBondEntryInput {
        id: optional_id(m)?,
        donors,
        acceptor: AtomRef::from_edn(required_key(m, "acceptor", "dative-bond-entry")?)?,
        bond: DativeBondDsl::from_edn(required_key(m, "type", "dative-bond-entry")?)?,
    })
}

pub(super) fn parse_aromatic_system_entry(
    edn: &Edn<'_>,
) -> Result<AromaticSystemEntryInput, DeError> {
    let m = expect_map(edn, "aromatic-system-entry")?;
    let system = AromaticSystemDsl::from_edn(required_key(m, "type", "aromatic-system-entry")?)?;
    Ok(AromaticSystemEntryInput {
        id: optional_id(m)?,
        atoms: atoms_vec(m, "aromatic-system-entry")?,
        system,
    })
}

pub(super) fn parse_multicenter_bond_entry(
    edn: &Edn<'_>,
) -> Result<MulticenterBondEntryInput, DeError> {
    let m = expect_map(edn, "multicenter-bond-entry")?;
    let bond = MulticenterBondDsl::from_edn(required_key(m, "type", "multicenter-bond-entry")?)?;
    Ok(MulticenterBondEntryInput {
        id: optional_id(m)?,
        atoms: atoms_vec(m, "multicenter-bond-entry")?,
        bond,
    })
}

pub(super) fn parse_noncovalent_bond_entry(
    edn: &Edn<'_>,
) -> Result<NoncovalentBondEntryInput, DeError> {
    let m = expect_map(edn, "noncovalent-bond-entry")?;
    let [a, b] = atoms_pair(m, "noncovalent-bond-entry")?;
    Ok(NoncovalentBondEntryInput {
        id: optional_id(m)?,
        first: a,
        second: b,
        bond: NoncovalentBondDsl::from_edn(required_key(m, "type", "noncovalent-bond-entry")?)?,
    })
}

pub(super) fn parse_stereo_atom_entry(edn: &Edn<'_>) -> Result<StereoAtomEntryInput, DeError> {
    let m = expect_map(edn, "stereo-atom-entry")?;
    Ok(StereoAtomEntryInput {
        id: optional_id(m)?,
        site: AtomRef::from_edn(required_key(m, "site", "stereo-atom-entry")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-atom-entry")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
        stereo: StereoAtomDsl::from_edn(required_key(m, "type", "stereo-atom-entry")?)?,
    })
}

pub(super) fn parse_stereo_bond_entry(edn: &Edn<'_>) -> Result<StereoBondEntryInput, DeError> {
    let m = expect_map(edn, "stereo-bond-entry")?;
    Ok(StereoBondEntryInput {
        id: optional_id(m)?,
        site: BondRef::from_edn(required_key(m, "site", "stereo-bond-entry")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-bond-entry")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
        stereo: StereoBondDsl::from_edn(required_key(m, "type", "stereo-bond-entry")?)?,
    })
}

pub(super) fn parse_atom_aliases(edn: &Edn<'_>) -> Result<Vec<(String, Box<AtomDsl>)>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of keyword/atom-string pairs",
            got: edn.kind(),
            path: vec![":atom-aliases".into()],
        });
    };
    if !v.len().is_multiple_of(2) {
        return Err(DeError::Custom(
            ":atom-aliases must have even length (keyword/atom-string pairs)".into(),
        ));
    }
    let mut out = Vec::with_capacity(v.len() / 2);
    for pair in v.chunks(2) {
        let Edn::Keyword(name) = &pair[0] else {
            return Err(DeError::TypeMismatch {
                expected: "keyword (alias name)",
                got: pair[0].kind(),
                path: vec![":atom-aliases".into()],
            });
        };
        let Edn::Str(s) = &pair[1] else {
            return Err(DeError::TypeMismatch {
                expected: "atom-string",
                got: pair[1].kind(),
                path: vec![":atom-aliases".into()],
            });
        };
        let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
        out.push((name.name().to_string(), Box::new(dsl)));
    }
    Ok(out)
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

#[cfg(test)]
mod tests;
