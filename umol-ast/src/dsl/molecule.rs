//! Molecule DSL.
//!
//! `MoleculeDsl` wraps a `MoleculeAst` together with the `Metadata` that records
//! the surface-form id/alias bindings (atom ids, bond ids, etc.). The EDN
//! form is a map keyed by `:atoms`, `:bonds`, `:dative`, `:aromatic`,
//! `:multicenter`, `:noncovalent`, `:atom-aliases`/`:aliases`, and
//! `:constraints`. Each entity delegates to its own entity DSL. Constraints
//! parse directly into the typed `Constraint` tree.

// Closures like `|e| T::from_edn(e)` passed to `parse_vec` can't be replaced
// by bare `T::from_edn` — type-erasing the fn item loses the `for<'a>` HRTB
// on the `FromEdn<'a>` impl.
#![allow(clippy::redundant_closure)]

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::mem::take;
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::{BondDsl, expand_bond_keyword};
use super::config::MoleculeDefaults;
use super::constraint::{
    eof_err, missing, read_atom_ref, read_constraints_dsl, read_map, read_vec,
    unexpected_byte_kind, AtomRef, ConstraintDsl, ConstraintsDsl, EntityCounts,
};
use super::dative::DativeBondDsl;
use super::error::ParseError;
use super::multicenter::MulticenterBondDsl;
use super::noncovalent::NoncovalentBondDsl;
use super::value::ValueDsl;
use crate::ast::value::ValueAst;
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::dative::DativeBondAst;
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::NoncovalentBondAst;
use crate::ast::traits::{FromAst, IntoAst};

/// Surface-form metadata paired with a `MoleculeAst`. Records atom ids,
/// per-entity ids, and the atom-alias table. Never drifts onto a different
/// AST: `MoleculeDsl` keeps both fields private and rewraps atomically
/// through `from_parts`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    atom_ids: IndexMap<AtomIdx, String>,
    atom_aliases: BiMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondIdx, String>,
    dative_bond_ids: IndexMap<DativeBondIdx, String>,
    aromatic_system_ids: IndexMap<AromaticSystemIdx, String>,
    multicenter_bond_ids: IndexMap<MulticenterBondIdx, String>,
    noncovalent_bond_ids: IndexMap<NoncovalentBondIdx, String>,
}

impl Metadata {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom_id(&self, idx: AtomIdx) -> Option<&str> {
        self.atom_ids.get(&idx).map(String::as_str)
    }

    pub fn bond_id(&self, idx: BondIdx) -> Option<&str> {
        self.bond_ids.get(&idx).map(String::as_str)
    }

    pub fn dative_bond_id(&self, idx: DativeBondIdx) -> Option<&str> {
        self.dative_bond_ids.get(&idx).map(String::as_str)
    }

    pub fn aromatic_system_id(&self, idx: AromaticSystemIdx) -> Option<&str> {
        self.aromatic_system_ids.get(&idx).map(String::as_str)
    }

    pub fn multicenter_bond_id(&self, idx: MulticenterBondIdx) -> Option<&str> {
        self.multicenter_bond_ids.get(&idx).map(String::as_str)
    }

    pub fn noncovalent_bond_id(&self, idx: NoncovalentBondIdx) -> Option<&str> {
        self.noncovalent_bond_ids.get(&idx).map(String::as_str)
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

    pub fn set_atom_id(&mut self, idx: AtomIdx, name: impl Into<String>) {
        self.atom_ids.insert(idx, name.into());
    }

    pub fn set_bond_id(&mut self, idx: BondIdx, name: impl Into<String>) {
        self.bond_ids.insert(idx, name.into());
    }

    pub fn set_dative_bond_id(&mut self, idx: DativeBondIdx, name: impl Into<String>) {
        self.dative_bond_ids.insert(idx, name.into());
    }

    pub fn set_aromatic_system_id(&mut self, idx: AromaticSystemIdx, name: impl Into<String>) {
        self.aromatic_system_ids.insert(idx, name.into());
    }

    pub fn set_multicenter_bond_id(&mut self, idx: MulticenterBondIdx, name: impl Into<String>) {
        self.multicenter_bond_ids.insert(idx, name.into());
    }

    pub fn set_noncovalent_bond_id(&mut self, idx: NoncovalentBondIdx, name: impl Into<String>) {
        self.noncovalent_bond_ids.insert(idx, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases
            .insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_id(mut self, idx: AtomIdx, name: impl Into<String>) -> Self {
        self.set_atom_id(idx, name);
        self
    }

    pub fn with_bond_id(mut self, idx: BondIdx, name: impl Into<String>) -> Self {
        self.set_bond_id(idx, name);
        self
    }

    pub fn with_dative_bond_id(mut self, idx: DativeBondIdx, name: impl Into<String>) -> Self {
        self.set_dative_bond_id(idx, name);
        self
    }

    pub fn with_aromatic_system_id(
        mut self,
        idx: AromaticSystemIdx,
        name: impl Into<String>,
    ) -> Self {
        self.set_aromatic_system_id(idx, name);
        self
    }

    pub fn with_multicenter_bond_id(
        mut self,
        idx: MulticenterBondIdx,
        name: impl Into<String>,
    ) -> Self {
        self.set_multicenter_bond_id(idx, name);
        self
    }

    pub fn with_noncovalent_bond_id(
        mut self,
        idx: NoncovalentBondIdx,
        name: impl Into<String>,
    ) -> Self {
        self.set_noncovalent_bond_id(idx, name);
        self
    }

    pub fn with_atom_alias(
        mut self,
        name: impl Into<String>,
        atom: impl Into<AtomDsl>,
    ) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

/// Surface DSL for a whole molecule. Pairs `MoleculeAst` with `Metadata`;
/// fields are private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    ast: MoleculeAst,
    metadata: Metadata,
}

impl MoleculeDsl {
    pub fn from_parts(ast: MoleculeAst, metadata: Metadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (MoleculeAst, Metadata) {
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
        let (ast, metadata) = input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(ast, metadata))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let mi = read_molecule_input(&mut de)?;
        de.expect_eof()?;
        let (ast, metadata) = mi.into_ast().map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(ast, metadata))
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
/// [`Metadata`] and call [`MoleculeDsl::to_edn`].
impl ToEdn for MoleculeAst {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(self, &Metadata::default())
    }
}

// region: Streaming parser
//
// Single-pass parse of the molecule map directly off the source text, using
// only the byte-level / typed primitives on `EdnStreamDeserializer`. Each
// non-terminal in the grammar has one `read_*` function; no detour through
// an intermediate `Edn` tree.

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
            "dative" => mi.dative_bonds = read_vec(de, read_dative_bond_entry)?,
            "aromatic" => mi.aromatic_systems = read_vec(de, read_aromatic_system_entry)?,
            "multicenter" => mi.multicenter_bonds = read_vec(de, read_multicenter_bond_entry)?,
            "noncovalent" => mi.noncovalent_bonds = read_vec(de, read_noncovalent_bond_entry)?,
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

fn read_atom_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomEntryInput, EdnError> {
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
            let expanded =
                expand_bond_keyword(name.as_ref()).ok_or_else(|| {
                    DeError::Custom(format!("unknown bond keyword :{}", name))
                })?;
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
            let expanded = super::dative::expand_dative_keyword(name.as_ref()).ok_or_else(
                || DeError::Custom(format!("unknown dative keyword :{}", name)),
            )?;
            Cow::Borrowed(expanded)
        }
        _ => de.read_string()?,
    };
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("dative", e).into())
}

fn read_bond_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<BondEntryInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'[' => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_dsl(de)?;
            de.consume_byte(b']')?;
            Ok(BondEntryInput {
                id: None,
                a,
                b,
                bond,
            })
        }
        b'{' => {
            let mut id = None;
            let mut a = None;
            let mut b = None;
            let mut bond = None;
            read_map(de, |de, key| {
                match key {
                    "id" => id = Some(de.read_keyword_name()?.into_owned()),
                    "a" => a = Some(read_atom_ref(de)?),
                    "b" => b = Some(read_atom_ref(de)?),
                    "type" => bond = Some(read_bond_dsl(de)?),
                    _ => de.read_skip_value()?,
                }
                Ok(())
            })?;
            Ok(BondEntryInput {
                id,
                a: a.ok_or_else(|| missing("a", "bond-entry"))?,
                b: b.ok_or_else(|| missing("b", "bond-entry"))?,
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

fn read_dative_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DativeBondEntryInput, EdnError> {
    let mut id = None;
    let mut donor = None;
    let mut acceptor = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "donor" => donor = Some(read_atom_ref(de)?),
            "acceptor" => acceptor = Some(read_atom_ref(de)?),
            "type" => bond = Some(read_dative_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(DativeBondEntryInput {
        id,
        donor: donor.ok_or_else(|| missing("donor", "dative-bond-entry"))?,
        acceptor: acceptor.ok_or_else(|| missing("acceptor", "dative-bond-entry"))?,
        bond: bond.ok_or_else(|| missing("type", "dative-bond-entry"))?,
    })
}

fn read_aromatic_system_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystemEntryInput, EdnError> {
    let mut id = None;
    let mut atoms = None;
    let mut system = None;
    let mut electrons = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "electrons" => electrons = Some(read_value_vec(de)?),
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
    let mut system = system.ok_or_else(|| missing("type", "aromatic-system-entry"))?;
    if let Some(es) = electrons {
        system.0.electrons = es;
    }
    Ok(AromaticSystemEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "aromatic-system-entry"))?,
        system,
    })
}

fn read_value_vec(de: &mut EdnStreamDeserializer<'_>) -> Result<Vec<ValueAst>, EdnError> {
    read_vec(de, |de| {
        let slice = de.read_value_slice()?;
        let edn = umol_edn::read_string(slice)?;
        ValueDsl::from_edn(&edn)
            .map(|v| v.0)
            .map_err(EdnError::from)
    })
}

fn read_multicenter_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBondEntryInput, EdnError> {
    let mut id = None;
    let mut atoms = None;
    let mut bond = None;
    let mut electrons = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "electrons" => electrons = Some(read_value_vec(de)?),
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
    let mut bond = bond.ok_or_else(|| missing("type", "multicenter-bond-entry"))?;
    if let Some(es) = electrons {
        bond.0.electrons = es;
    }
    Ok(MulticenterBondEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "multicenter-bond-entry"))?,
        bond,
    })
}

fn read_noncovalent_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondEntryInput, EdnError> {
    let mut id = None;
    let mut a = None;
    let mut b = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "a" => a = Some(read_atom_ref(de)?),
            "b" => b = Some(read_atom_ref(de)?),
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
    Ok(NoncovalentBondEntryInput {
        id,
        a: a.ok_or_else(|| missing("a", "noncovalent-bond-entry"))?,
        b: b.ok_or_else(|| missing("b", "noncovalent-bond-entry"))?,
        bond: bond.ok_or_else(|| missing("type", "noncovalent-bond-entry"))?,
    })
}

fn read_atom_aliases(
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
        // `DativeBondDsl` and `NoncovalentBondDsl` use unit-shaped defaults
        // (empty struct), so there is nothing to strip here.
        MoleculeDsl {
            ast: ast_out,
            metadata: Metadata::default(),
        }
    }
}

impl IntoAst<MoleculeAst> for MoleculeDsl {
    type Ctx = MoleculeDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> MoleculeAst {
        let mut ast = self.ast;
        for atom in ast.atoms_mut() {
            *atom = AtomDsl(take(atom)).into_ast(&cfg.atom);
        }
        for bond in ast.bonds_mut() {
            *bond = BondDsl(take(bond)).into_ast(&cfg.bond);
        }
        for system in ast.aromatic_systems_mut() {
            *system = AromaticSystemDsl(take(system)).into_ast(&cfg.aromatic_system);
        }
        for bond in ast.multicenter_bonds_mut() {
            *bond = MulticenterBondDsl(take(bond)).into_ast(&cfg.multicenter_bond);
        }
        ast
    }
}

// endregion: Streaming parser

// region: Render

fn render_molecule_edn(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(8);
    map.insert(Edn::keyword("atoms"), render_atoms(ast, meta));
    map.insert(Edn::keyword("bonds"), render_bonds(ast, meta));
    if ast.dative_bonds().count() > 0 {
        map.insert(Edn::keyword("dative"), render_dative(ast, meta));
    }
    if ast.aromatic_systems().count() > 0 {
        map.insert(Edn::keyword("aromatic"), render_aromatic(ast, meta));
    }
    if ast.multicenter_bonds().count() > 0 {
        map.insert(Edn::keyword("multicenter"), render_multicenter(ast, meta));
    }
    if ast.noncovalent_bonds().count() > 0 {
        map.insert(Edn::keyword("noncovalent"), render_noncovalent(ast, meta));
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

fn render_atoms(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .atoms()
        .iter()
        .map(|view| render_atom_entry(view.idx, view.data, meta))
        .collect();
    Edn::Vector(entries.into())
}

fn render_atom_entry(idx: AtomIdx, atom: &AtomAst, meta: &Metadata) -> Edn<'static> {
    let dsl = AtomDsl::from_ref(atom);
    let spec = if let Some(alias) = meta.atom_alias_for(dsl) {
        Edn::Keyword(EdnKeyword::owned(alias.to_string()))
    } else {
        dsl.to_edn()
    };
    match meta.atom_id(idx) {
        Some(id) => Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(id.to_string())), spec].into()),
        None => spec,
    }
}

fn render_atom_ref(idx: AtomIdx, meta: &Metadata) -> Edn<'static> {
    match meta.atom_id(idx) {
        Some(id) => Edn::Keyword(EdnKeyword::owned(id.to_string())),
        None => Edn::Int(idx.index() as i64),
    }
}

fn render_bonds(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .map(|view| {
            let bond_edn = BondDsl::from_ref(view.data).to_edn();
            let a = render_atom_ref(view.atoms()[0], meta);
            let b = render_atom_ref(view.atoms()[1], meta);
            match meta.bond_id(view.idx) {
                Some(id) => {
                    let mut m = EdnMap::with_capacity(4);
                    m.insert(
                        Edn::keyword("id"),
                        Edn::Keyword(EdnKeyword::owned(id.to_string())),
                    );
                    m.insert(Edn::keyword("a"), a);
                    m.insert(Edn::keyword("b"), b);
                    m.insert(Edn::keyword("type"), bond_edn);
                    Edn::Map(m)
                }
                None => Edn::Vector(vec![a, b, bond_edn].into()),
            }
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_dative(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .dative_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.dative_bond_id(view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.to_string())),
                );
            }
            let donors: Vec<AtomIdx> = view.donors().collect();
            let donor_edn = if donors.len() == 1 {
                render_atom_ref(donors[0], meta)
            } else {
                Edn::Vector(
                    donors
                        .into_iter()
                        .map(|a| render_atom_ref(a, meta))
                        .collect::<Vec<_>>()
                        .into(),
                )
            };
            m.insert(Edn::keyword("donor"), donor_edn);
            m.insert(
                Edn::keyword("acceptor"),
                render_atom_ref(view.acceptor, meta),
            );
            m.insert(
                Edn::keyword("type"),
                DativeBondDsl::from_ref(view.data).to_edn(),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_aromatic(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .aromatic_systems()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.aromatic_system_id(view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.to_string())),
                );
            }
            let atoms: Vec<Edn<'static>> = view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            if !view.data.electrons.is_empty() {
                m.insert(
                    Edn::keyword("electrons"),
                    render_value_vec(&view.data.electrons),
                );
            }
            m.insert(
                Edn::keyword("type"),
                Edn::Str(Cow::Owned(AromaticSystemDsl::from_ref(view.data).to_string())),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_multicenter(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .multicenter_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.multicenter_bond_id(view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.to_string())),
                );
            }
            let atoms: Vec<Edn<'static>> = view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            if !view.data.electrons.is_empty() {
                m.insert(
                    Edn::keyword("electrons"),
                    render_value_vec(&view.data.electrons),
                );
            }
            m.insert(
                Edn::keyword("type"),
                Edn::Str(Cow::Owned(MulticenterBondDsl::from_ref(view.data).to_string())),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_value_vec(values: &[ValueAst]) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = values
        .iter()
        .map(|v| ValueDsl(v.clone()).to_edn())
        .collect();
    Edn::Vector(entries.into())
}

fn render_noncovalent(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.noncovalent_bond_id(view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.to_string())),
                );
            }
            m.insert(Edn::keyword("a"), render_atom_ref(view.atoms()[0], meta));
            m.insert(Edn::keyword("b"), render_atom_ref(view.atoms()[1], meta));
            m.insert(
                Edn::keyword("type"),
                NoncovalentBondDsl::from_ref(view.data).to_edn(),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_atom_aliases(meta: &Metadata) -> Edn<'static> {
    let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases_len() * 2);
    for (name, dsl) in meta.iter_atom_aliases() {
        pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
        pairs.push(dsl.to_edn());
    }
    Edn::Vector(pairs.into())
}

// endregion: Render

// region: Private parse intermediate
//
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
    pub(crate) atom_aliases: Vec<(String, Box<AtomDsl>)>,
    pub(crate) constraints: Vec<ConstraintDsl>,
}

/// Atom entry in a parsed molecule map. Mirrors the DSL spec §4 grammar
/// `atom-entry ::= atom-spec | [ keyword atom-spec ]`.
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) donor: AtomRef,
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
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: NoncovalentBondDsl,
}

impl MoleculeInput {
    /// Destructive lowering: consumes the input, resolves refs against the
    /// built id scopes, and produces the final `MoleculeAst` with its
    /// `Metadata`. Called from `FromEdn::from_edn` and the streaming path.
    pub(crate) fn into_ast(self) -> Result<(MoleculeAst, Metadata), ParseError> {
        let MoleculeInput {
            atoms: atom_entries,
            bonds: bond_entries,
            dative_bonds: dative_entries,
            aromatic_systems: aromatic_entries,
            multicenter_bonds: multicenter_entries,
            noncovalent_bonds: noncovalent_entries,
            atom_aliases: alias_entries,
            constraints: constraint_dsls,
        } = self;

        // Alias table: bijective. The parser enforces both directions —
        // duplicate names and duplicate atom-dsls are rejected at parse
        // time. Programmatic `Metadata::add_atom_alias` is last-wins.
        let mut alias_table: IndexMap<String, Box<AtomDsl>> = IndexMap::new();
        for (name, dsl) in alias_entries {
            if alias_table.contains_key(&name) {
                return Err(ParseError::DuplicateId(name));
            }
            if alias_table.values().any(|existing| existing == &dsl) {
                return Err(ParseError::InvalidValue(
                    "atom-aliases must be bijective: two names map to the same atom".into(),
                ));
            }
            alias_table.insert(name, dsl);
        }

        // Atoms: materialize AtomAst from each entry; collect ids.
        let mut atoms: Vec<AtomAst> = Vec::with_capacity(atom_entries.len());
        let mut metadata = Metadata::new();
        let mut atom_id_to_idx: IndexMap<String, AtomIdx> = IndexMap::new();
        for (pos, entry) in atom_entries.into_iter().enumerate() {
            let idx = AtomIdx(pos as u32);
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                atom_id_to_idx.insert(id.clone(), idx);
                metadata.set_atom_id(idx, id);
            }
            let ast = match entry.spec {
                AtomSpecInput::Bare(dsl) => dsl.0,
                AtomSpecInput::Alias(name) => match alias_table.get(&name) {
                    Some(dsl) => dsl.0.clone(),
                    None => {
                        return Err(ParseError::InvalidValue(format!(
                            "unknown atom alias :{}",
                            name
                        )));
                    }
                },
            };
            atoms.push(ast);
        }

        let atom_count = atoms.len();

        // Bonds.
        let mut bonds: Vec<(AtomIdx, AtomIdx, BondAst)> = Vec::with_capacity(bond_entries.len());
        let mut entry_ids: IndexMap<String, ()> = IndexMap::new();
        for (pos, entry) in bond_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                metadata.set_bond_id(BondIdx(pos as u32), id);
            }
            let a = entry.a.resolve(atom_count, &atom_id_to_idx)?;
            let b = entry.b.resolve(atom_count, &atom_id_to_idx)?;
            bonds.push((a, b, entry.bond.0));
        }

        // Dative bonds.
        let mut dative_list: Vec<(Vec<AtomIdx>, AtomIdx, DativeBondAst)> =
            Vec::with_capacity(dative_entries.len());
        for (pos, entry) in dative_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                metadata.set_dative_bond_id(DativeBondIdx(pos as u32), id);
            }
            let donor = entry.donor.resolve(atom_count, &atom_id_to_idx)?;
            let acceptor = entry.acceptor.resolve(atom_count, &atom_id_to_idx)?;
            dative_list.push((vec![donor], acceptor, entry.bond.0));
        }

        // Aromatic systems.
        let mut aromatic_list: Vec<(Vec<AtomIdx>, AromaticSystemAst)> =
            Vec::with_capacity(aromatic_entries.len());
        for (pos, entry) in aromatic_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                metadata.set_aromatic_system_id(AromaticSystemIdx(pos as u32), id);
            }
            let atoms_resolved: Vec<AtomIdx> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(atom_count, &atom_id_to_idx))
                .collect::<Result<_, _>>()?;
            aromatic_list.push((atoms_resolved, entry.system.0));
        }

        // Multicenter bonds.
        let mut multicenter_list: Vec<(Vec<AtomIdx>, MulticenterBondAst)> =
            Vec::with_capacity(multicenter_entries.len());
        for (pos, entry) in multicenter_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                metadata.set_multicenter_bond_id(MulticenterBondIdx(pos as u32), id);
            }
            let atoms_resolved: Vec<AtomIdx> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(atom_count, &atom_id_to_idx))
                .collect::<Result<_, _>>()?;
            multicenter_list.push((atoms_resolved, entry.bond.0));
        }

        // Noncovalent bonds.
        let mut noncovalent_list: Vec<(AtomIdx, AtomIdx, NoncovalentBondAst)> =
            Vec::with_capacity(noncovalent_entries.len());
        for (pos, entry) in noncovalent_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                metadata.set_noncovalent_bond_id(NoncovalentBondIdx(pos as u32), id);
            }
            let a = entry.a.resolve(atom_count, &atom_id_to_idx)?;
            let b = entry.b.resolve(atom_count, &atom_id_to_idx)?;
            noncovalent_list.push((a, b, entry.bond.0));
        }

        // Atom aliases. Names are guaranteed unique by the upstream
        // `parse_aliases` dedup; `add_atom_alias` is last-wins on
        // duplicate atom-dsl, which can't fire here.
        for (name, dsl) in alias_table {
            metadata.add_atom_alias(name, *dsl);
        }

        // Resolve constraint refs against the final metadata + counts.
        let counts = EntityCounts {
            atom_count,
            bond_count: bonds.len(),
            dative_bond_count: dative_list.len(),
            aromatic_system_count: aromatic_list.len(),
            multicenter_bond_count: multicenter_list.len(),
            noncovalent_bond_count: noncovalent_list.len(),
        };
        let constraints = ConstraintsDsl(constraint_dsls).into_ast(&counts, &metadata)?;

        let ast = MoleculeAst::from_parts(
            atoms,
            bonds,
            dative_list,
            aromatic_list,
            multicenter_list,
            noncovalent_list,
            constraints,
        );
        Ok((ast, metadata))
    }
}

// endregion: Private parse intermediate

// region: Parse

fn parse_molecule_input(edn: &Edn<'_>) -> Result<MoleculeInput, DeError> {
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
            "dative" => input.dative_bonds = parse_vec(v, ":dative", parse_dative_bond_entry)?,
            "aromatic" => {
                input.aromatic_systems = parse_vec(v, ":aromatic", parse_aromatic_system_entry)?
            }
            "multicenter" => {
                input.multicenter_bonds =
                    parse_vec(v, ":multicenter", parse_multicenter_bond_entry)?
            }
            "noncovalent" => {
                input.noncovalent_bonds =
                    parse_vec(v, ":noncovalent", parse_noncovalent_bond_entry)?
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

fn parse_vec<T>(
    edn: &Edn<'_>,
    context: &'static str,
    mut f: impl FnMut(&Edn<'_>) -> Result<T, DeError>,
) -> Result<Vec<T>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    v.iter().map(|e| f(e)).collect()
}

fn parse_atom_entry(edn: &Edn<'_>) -> Result<AtomEntryInput, DeError> {
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

fn parse_bond_entry(edn: &Edn<'_>) -> Result<BondEntryInput, DeError> {
    match edn {
        Edn::Vector(v) if v.len() == 3 => Ok(BondEntryInput {
            id: None,
            a: AtomRef::from_edn(&v[0])?,
            b: AtomRef::from_edn(&v[1])?,
            bond: BondDsl::from_edn(&v[2])?,
        }),
        Edn::Map(m) => Ok(BondEntryInput {
            id: optional_id(m)?,
            a: AtomRef::from_edn(required_key(m, "a", "bond-entry")?)?,
            b: AtomRef::from_edn(required_key(m, "b", "bond-entry")?)?,
            bond: BondDsl::from_edn(required_key(m, "type", "bond-entry")?)?,
        }),
        other => Err(DeError::TypeMismatch {
            expected: "bond-entry map or 3-vec",
            got: other.kind(),
            path: vec!["bond-entry".into()],
        }),
    }
}

fn parse_dative_bond_entry(edn: &Edn<'_>) -> Result<DativeBondEntryInput, DeError> {
    let m = expect_map(edn, "dative-bond-entry")?;
    Ok(DativeBondEntryInput {
        id: optional_id(m)?,
        donor: AtomRef::from_edn(required_key(m, "donor", "dative-bond-entry")?)?,
        acceptor: AtomRef::from_edn(required_key(m, "acceptor", "dative-bond-entry")?)?,
        bond: DativeBondDsl::from_edn(required_key(m, "type", "dative-bond-entry")?)?,
    })
}

fn parse_aromatic_system_entry(edn: &Edn<'_>) -> Result<AromaticSystemEntryInput, DeError> {
    let m = expect_map(edn, "aromatic-system-entry")?;
    let mut system =
        AromaticSystemDsl::from_edn(required_key(m, "type", "aromatic-system-entry")?)?;
    if let Some(edn) = m.get_keyword("electrons") {
        system.0.electrons = parse_value_vec(edn, ":electrons")?;
    }
    Ok(AromaticSystemEntryInput {
        id: optional_id(m)?,
        atoms: parse_vec(
            required_key(m, "atoms", "aromatic-system-entry")?,
            ":atoms",
            |e| AtomRef::from_edn(e),
        )?,
        system,
    })
}

fn parse_multicenter_bond_entry(edn: &Edn<'_>) -> Result<MulticenterBondEntryInput, DeError> {
    let m = expect_map(edn, "multicenter-bond-entry")?;
    let mut bond =
        MulticenterBondDsl::from_edn(required_key(m, "type", "multicenter-bond-entry")?)?;
    if let Some(edn) = m.get_keyword("electrons") {
        bond.0.electrons = parse_value_vec(edn, ":electrons")?;
    }
    Ok(MulticenterBondEntryInput {
        id: optional_id(m)?,
        atoms: parse_vec(
            required_key(m, "atoms", "multicenter-bond-entry")?,
            ":atoms",
            |e| AtomRef::from_edn(e),
        )?,
        bond,
    })
}

fn parse_value_vec(edn: &Edn<'_>, label: &'static str) -> Result<Vec<ValueAst>, DeError> {
    parse_vec(edn, label, |e| ValueDsl::from_edn(e).map(|v| v.0))
}

fn parse_noncovalent_bond_entry(edn: &Edn<'_>) -> Result<NoncovalentBondEntryInput, DeError> {
    let m = expect_map(edn, "noncovalent-bond-entry")?;
    Ok(NoncovalentBondEntryInput {
        id: optional_id(m)?,
        a: AtomRef::from_edn(required_key(m, "a", "noncovalent-bond-entry")?)?,
        b: AtomRef::from_edn(required_key(m, "b", "noncovalent-bond-entry")?)?,
        bond: NoncovalentBondDsl::from_edn(required_key(m, "type", "noncovalent-bond-entry")?)?,
    })
}

fn parse_atom_aliases(edn: &Edn<'_>) -> Result<Vec<(String, Box<AtomDsl>)>, DeError> {
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

fn required_key<'e>(
    m: &'e EdnMap<'e>,
    key: &'static str,
    context: &'static str,
) -> Result<&'e Edn<'e>, DeError> {
    m.get_keyword(key).ok_or_else(|| DeError::MissingField {
        key: key.to_string(),
        path: vec![context.into()],
    })
}

fn optional_id(m: &EdnMap<'_>) -> Result<Option<String>, DeError> {
    match m.get_keyword("id") {
        Some(Edn::Keyword(k)) => Ok(Some(k.name().to_string())),
        Some(other) => Err(DeError::TypeMismatch {
            expected: "keyword id",
            got: other.kind(),
            path: vec![":id".into()],
        }),
        None => Ok(None),
    }
}

/// Check that `id` is not already claimed by an atom id or alias name.
fn check_id_disjoint(
    id: &str,
    atom_id_to_idx: &IndexMap<String, AtomIdx>,
    alias_table: &IndexMap<String, Box<AtomDsl>>,
) -> Result<(), ParseError> {
    if atom_id_to_idx.contains_key(id) || alias_table.contains_key(id) {
        return Err(ParseError::DuplicateId(id.to_string()));
    }
    Ok(())
}

// endregion: Parse

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;
    use umol_shared::element::Element;

    use super::*;
    use crate::{dsl, mol};
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{Constraint, Constraints, MoleculeConstraint};
    use crate::ast::spin::SpinStateAst;
    use crate::ast::value::ValueAst;

    #[rstest]
    fn test_metadata_new() {
        let m = Metadata::new();
        assert_eq!(m, Metadata::default());
        assert!(m.atom_id(AtomIdx(0)).is_none());
        assert!(m.bond_id(BondIdx(0)).is_none());
        assert!(m.dative_bond_id(DativeBondIdx(0)).is_none());
        assert!(m.aromatic_system_id(AromaticSystemIdx(0)).is_none());
        assert!(m.multicenter_bond_id(MulticenterBondIdx(0)).is_none());
        assert!(m.noncovalent_bond_id(NoncovalentBondIdx(0)).is_none());
        assert!(!m.has_atom_aliases());
        assert_eq!(m.atom_aliases_len(), 0);
    }

    #[rstest]
    fn test_metadata_set_atom_id() {
        let mut m = Metadata::new();
        m.set_atom_id(AtomIdx(0), "c1");
        assert_eq!(m.atom_id(AtomIdx(0)), Some("c1"));
    }

    #[rstest]
    fn test_metadata_set_atom_id_last_wins() {
        let mut m = Metadata::new();
        m.set_atom_id(AtomIdx(0), "old");
        m.set_atom_id(AtomIdx(0), "new");
        assert_eq!(m.atom_id(AtomIdx(0)), Some("new"));
    }

    #[rstest]
    fn test_metadata_set_bond_id() {
        let mut m = Metadata::new();
        m.set_bond_id(BondIdx(2), "b1");
        assert_eq!(m.bond_id(BondIdx(2)), Some("b1"));
    }

    #[rstest]
    fn test_metadata_set_dative_bond_id() {
        let mut m = Metadata::new();
        m.set_dative_bond_id(DativeBondIdx(1), "d1");
        assert_eq!(m.dative_bond_id(DativeBondIdx(1)), Some("d1"));
    }

    #[rstest]
    fn test_metadata_set_aromatic_system_id() {
        let mut m = Metadata::new();
        m.set_aromatic_system_id(AromaticSystemIdx(0), "ring1");
        assert_eq!(m.aromatic_system_id(AromaticSystemIdx(0)), Some("ring1"));
    }

    #[rstest]
    fn test_metadata_set_multicenter_bond_id() {
        let mut m = Metadata::new();
        m.set_multicenter_bond_id(MulticenterBondIdx(0), "mc1");
        assert_eq!(m.multicenter_bond_id(MulticenterBondIdx(0)), Some("mc1"));
    }

    #[rstest]
    fn test_metadata_set_noncovalent_bond_id() {
        let mut m = Metadata::new();
        m.set_noncovalent_bond_id(NoncovalentBondIdx(0), "h1");
        assert_eq!(m.noncovalent_bond_id(NoncovalentBondIdx(0)), Some("h1"));
    }

    #[rstest]
    fn test_metadata_add_atom_alias() {
        let mut m = Metadata::new();
        let atom = AtomAst::from_element(Element::C).with_implicit_hydrogens(2_i64);
        m.add_atom_alias("HC2", atom.clone());
        assert!(m.has_atom_alias("HC2"));
        assert_eq!(m.atom_aliases_len(), 1);
        assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("HC2"));
    }

    #[rstest]
    fn test_metadata_add_atom_alias_duplicate_name_replaces_atom() {
        let mut m = Metadata::new();
        let first = AtomAst::from_element(Element::C);
        let second = AtomAst::from_element(Element::N);
        m.add_atom_alias("X", first.clone());
        m.add_atom_alias("X", second.clone());
        assert_eq!(m.atom_aliases_len(), 1);
        assert_eq!(m.atom_alias_for(&AtomDsl(second)), Some("X"));
        assert_eq!(m.atom_alias_for(&AtomDsl(first)), None);
    }

    #[rstest]
    fn test_metadata_add_atom_alias_duplicate_atom_replaces_name() {
        let mut m = Metadata::new();
        let atom = AtomAst::from_element(Element::C);
        m.add_atom_alias("first", atom.clone());
        m.add_atom_alias("second", atom.clone());
        assert_eq!(m.atom_aliases_len(), 1);
        assert!(!m.has_atom_alias("first"));
        assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("second"));
    }

    #[rstest]
    fn test_metadata_iter_atom_aliases() {
        let m = Metadata::new()
            .with_atom_alias("a", AtomAst::from_element(Element::C))
            .with_atom_alias("b", AtomAst::from_element(Element::N));
        let collected: Vec<(&str, AtomAst)> = m
            .iter_atom_aliases()
            .map(|(name, dsl)| (name, dsl.0.clone()))
            .collect();
        assert_eq!(collected.len(), 2);
        assert!(collected.contains(&("a", AtomAst::from_element(Element::C))));
        assert!(collected.contains(&("b", AtomAst::from_element(Element::N))));
    }

    #[rstest]
    fn test_metadata_with_atom_id_chains() {
        let m = Metadata::new()
            .with_atom_id(AtomIdx(0), "a")
            .with_atom_id(AtomIdx(1), "b");
        assert_eq!(m.atom_id(AtomIdx(0)), Some("a"));
        assert_eq!(m.atom_id(AtomIdx(1)), Some("b"));
    }

    #[rstest]
    fn test_metadata_with_bond_id() {
        let m = Metadata::new().with_bond_id(BondIdx(0), "b");
        assert_eq!(m.bond_id(BondIdx(0)), Some("b"));
    }

    #[rstest]
    fn test_metadata_with_dative_bond_id() {
        let m = Metadata::new().with_dative_bond_id(DativeBondIdx(0), "d");
        assert_eq!(m.dative_bond_id(DativeBondIdx(0)), Some("d"));
    }

    #[rstest]
    fn test_metadata_with_aromatic_system_id() {
        let m = Metadata::new().with_aromatic_system_id(AromaticSystemIdx(0), "r");
        assert_eq!(m.aromatic_system_id(AromaticSystemIdx(0)), Some("r"));
    }

    #[rstest]
    fn test_metadata_with_multicenter_bond_id() {
        let m = Metadata::new().with_multicenter_bond_id(MulticenterBondIdx(0), "mc");
        assert_eq!(m.multicenter_bond_id(MulticenterBondIdx(0)), Some("mc"));
    }

    #[rstest]
    fn test_metadata_with_noncovalent_bond_id() {
        let m = Metadata::new().with_noncovalent_bond_id(NoncovalentBondIdx(0), "h");
        assert_eq!(m.noncovalent_bond_id(NoncovalentBondIdx(0)), Some("h"));
    }

    #[rstest]
    fn test_metadata_with_atom_alias() {
        let atom = AtomAst::from_element(Element::C);
        let m = Metadata::new().with_atom_alias("c", atom.clone());
        assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("c"));
    }

    #[rstest]
    fn test_metadata_mixed_chain() {
        let m = Metadata::new()
            .with_atom_id(AtomIdx(0), "c1")
            .with_bond_id(BondIdx(0), "b1")
            .with_atom_alias("X", AtomAst::from_element(Element::C));
        assert_eq!(m.atom_id(AtomIdx(0)), Some("c1"));
        assert_eq!(m.bond_id(BondIdx(0)), Some("b1"));
        assert!(m.has_atom_alias("X"));
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_empty() {
        let ast = MoleculeAst::default();
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("{:atoms [] :bonds []}").unwrap());
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_two_atoms_one_bond() {
        let dsl = dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let edn = dsl.to_edn();
        // Canonical render: order-1 default bond becomes the `:single` keyword.
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 :single]]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_with_id() {
        let dsl = dsl!(r#"{:atoms [[:c1 "C"] "C"] :bonds []}"#);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_bond_with_id_uses_map_form() {
        let dsl = dsl!(r#"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"}]}"#);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type :single}]}"##)
                .unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_alias_substituted() {
        let dsl = dsl!(r#"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"#);
        let edn = dsl.to_edn();
        // Both atoms match the alias — rendered as :x keyword references; the
        // alias table emits the :atom-aliases key.
        assert_eq!(
            edn,
            read_string(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_display_matches_edn() {
        let dsl = dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        assert_eq!(dsl.to_string(), dsl.to_edn().to_string());
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_omits_empty_optional_sections() {
        let dsl = dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else {
            panic!("expected map");
        };
        assert!(m.get_keyword("dative").is_none());
        assert!(m.get_keyword("aromatic").is_none());
        assert!(m.get_keyword("multicenter").is_none());
        assert!(m.get_keyword("noncovalent").is_none());
        assert!(m.get_keyword("atom-aliases").is_none());
        assert!(m.get_keyword("constraints").is_none());
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_empty() {
        let edn = read_string("{:atoms [] :bonds []}").unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().atoms().count(), 0);
        assert_eq!(dsl.ast().bonds().count(), 0);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_two_atoms_one_bond() {
        let edn = read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().atoms().count(), 2);
        assert_eq!(dsl.ast().bonds().count(), 1);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_atom_with_id() {
        let edn = read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.metadata().atom_id(AtomIdx(0)), Some("c1"));
        assert_eq!(dsl.metadata().atom_id(AtomIdx(1)), None);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_bond_map_form_with_id() {
        let edn =
            read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"}]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().bonds().count(), 1);
        assert_eq!(dsl.metadata().bond_id(BondIdx(0)), Some("b1"));
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_atom_aliases() {
        let edn = read_string(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().atoms().count(), 2);
        assert!(dsl.metadata().has_atom_alias("x"));
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_unknown_alias_errors() {
        let edn = read_string(r##"{:atoms [:x] :bonds []}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_duplicate_atom_id_errors() {
        let edn = read_string(r##"{:atoms [[:a "C"] [:a "N"]] :bonds []}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_unknown_top_level_key_errors() {
        let edn = read_string(r##"{:atoms [] :bonds [] :bogus 1}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip() {
        let source = r##"{:atoms ["C" "C" "O"] :bonds [[0 1 :single] [1 2 :single]]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        let rendered = dsl.to_edn();
        assert_eq!(rendered, edn);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_with_ids_and_aliases() {
        let source = r##"{:atoms [[:a "C"] [:b "C"] :x] :bonds [{:id :b1 :a :a :b :b :type :single} [:b 2 :double]] :atom-aliases [:x "N"]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        let rendered = dsl.to_edn();
        assert_eq!(rendered, edn);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##"{:atoms [] :bonds []}"##)]
    #[case::small(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##)]
    #[case::with_ids(r##"{:atoms [[:a "C"] [:b "N"]] :bonds [{:id :b1 :a :a :b :b :type "1"}]}"##)]
    #[case::inline_atom_constraints(r##"{:atoms ["C#v4" "N#R+"] :bonds []}"##)]
    #[case::inline_bond_constraint(r##"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"##)]
    #[case::aromatic_section(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic [{:id :ar1 :atoms [0 1 2 3 4 5] :type "#e6"}]}"##)]
    #[case::multicenter_section(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type "#e2"}]}"##)]
    #[case::dative_section(r##"{:atoms ["C" "N"] :bonds [] :dative [{:id :d1 :donor 0 :acceptor 1 :type "1#R"}]}"##)]
    #[case::noncovalent_section(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}]}"##)]
    #[case::atom_aliases(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##)]
    #[case::constraints_connected(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"##)]
    #[case::constraints_bond_order_sum(r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"} {:id :b2 :a 1 :b 2 :type "1"}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##)]
    #[case::constraints_atom_leaf_in_not(r##"{:atoms [[:c1 "C"]] :bonds [] :constraints [{:not {:atom [:c1 {:valence 3}]}}]}"##)]
    #[case::constraints_nested_combinators(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:and [{:or [{:atom [0 {:valence 3}]} {:atom [0 {:valence 4}]}]} {:not {:connected {:atoms [0 1]}}}]}]}"##)]
    #[case::constraints_sub_pattern(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}]}"##)]
    #[case::constraints_atom_degree(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:atom [0 {:degree 3}]}]}"##)]
    #[case::constraints_atom_connectivity(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:connectivity 4}]}]}"##)]
    #[case::constraints_atom_total_hydrogens(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-hydrogens 2}]}]}"##)]
    #[case::constraints_atom_ring_count(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-count 1}]}]}"##)]
    #[case::constraints_atom_ring_size(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-size 6}]}]}"##)]
    #[case::constraints_atom_ring_connectivity(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-connectivity 2}]}]}"##)]
    #[case::constraints_atom_donated_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:donated-pairs 1}]}]}"##)]
    #[case::constraints_atom_accepted_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:accepted-pairs 1}]}]}"##)]
    #[case::constraints_atom_aromatic_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence :not-aromatic}]}]}"##)]
    #[case::constraints_atom_aromatic_valence_with_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence {:aromatic 6}}]}]}"##)]
    #[case::constraints_atom_multicenter_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence :not-multicenter}]}]}"##)]
    #[case::constraints_atom_multicenter_valence_with_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence {:multicenter 3}}]}]}"##)]
    #[case::constraints_bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 :aromatic]}]}"##)]
    #[case::constraints_bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-count 1}]}]}"##)]
    #[case::constraints_bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-size 6}]}]}"##)]
    #[case::constraints_dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:ring-count 1}]}]}"##)]
    #[case::constraints_dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
    #[case::constraints_dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
    #[case::constraints_dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
    #[case::constraints_dative_donor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-donor-satisfies [0 {:valence 3}]}]}"##)]
    #[case::constraints_dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor-satisfies [0 {:valence 3}]}]}"##)]
    #[case::constraints_aromatic_system_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
    #[case::constraints_aromatic_system_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
    #[case::constraints_aromatic_system_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
    #[case::constraints_aromatic_system_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-any-atom [0 {:valence 4}]}]}"##)]
    #[case::constraints_multicenter_contains(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter [{:atoms [0 1 2] :type "#e3"}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
    #[case::constraints_multicenter_contains_all(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter [{:atoms [0 1 2] :type "#e3"}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
    #[case::constraints_multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type "#e2"}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
    #[case::constraints_multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type "#e2"}] :constraints [{:multicenter-bond-any-atom [0 {:valence 4}]}]}"##)]
    #[case::constraints_noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
    #[case::constraints_noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
    #[case::constraints_noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 3} {:valence 1}]]}]}"##)]
    #[case::constraints_sub_pattern_multi_entity_anchor(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]] :bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [[0 1 "1"]]}}}]}"##)]
    #[case::constraints_sub_pattern_dative_anchor(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:sub-pattern {:anchor {:dative-bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}]}}}]}"##)]
    #[case::constraints_sub_pattern_aromatic_system_anchor(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}] :constraints [{:sub-pattern {:anchor {:aromatic-systems [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type "#e2"}]}}}]}"##)]
    #[case::constraints_sub_pattern_multicenter_anchor(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type "#e2"}] :constraints [{:sub-pattern {:anchor {:multicenter-bonds [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type "#e2"}]}}}]}"##)]
    #[case::constraints_sub_pattern_noncovalent_anchor(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:sub-pattern {:anchor {:noncovalent-bonds [[0 0]]} :pattern {:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}]}}}]}"##)]
    fn test_molecule_dsl_from_edn_str_matches_from_edn(#[case] source: &str) {
        let via_str = MoleculeDsl::from_edn_str(source).unwrap();
        let tree = read_string(source).unwrap();
        let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
        assert_eq!(via_str, via_tree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_constraint_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:bogus 1}]}]}"##)]
    #[case::atom_aromatic_valence_unknown_keyword(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence :bogus}]}]}"##)]
    #[case::atom_aromatic_valence_unknown_inner_key(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence {:bogus 1}}]}]}"##)]
    #[case::atom_aromatic_valence_wrong_type(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence 42}]}]}"##)]
    #[case::atom_multicenter_valence_unknown_keyword(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence :bogus}]}]}"##)]
    #[case::atom_multicenter_valence_unknown_inner_key(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence {:bogus 1}}]}]}"##)]
    #[case::atom_multicenter_valence_wrong_type(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence 42}]}]}"##)]
    #[case::bond_constraint_unknown_keyword(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 :bogus]}]}"##)]
    #[case::bond_constraint_unknown_inner_key(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:bogus 1}]}]}"##)]
    #[case::bond_constraint_wrong_type(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 42]}]}"##)]
    #[case::dative_constraint_unknown_key(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:bogus 1}]}]}"##)]
    #[case::sub_pattern_anchor_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:bogus [[0 0]]} :pattern {:atoms ["C"] :bonds []}}}]}"##)]
    #[case::constraint_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:bogus 1}]}"##)]
    #[case::noncovalent_ends_satisfy_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 2}]]}]}"##)]
    #[case::noncovalent_ends_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0]]}]}"##)]
    fn test_molecule_dsl_from_edn_str_rejects_invalid_constraints(#[case] source: &str) {
        let result = MoleculeDsl::from_edn_str(source);
        assert!(
            result.is_err(),
            "expected parse failure, but got: {:?}",
            result,
        );
    }

    #[rstest]
    fn test_molecule_dsl_from_str_parses_edn_source() {
        let source = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;
        let dsl: MoleculeDsl = source.parse().unwrap();
        assert_eq!(dsl.ast().atoms().count(), 2);
        assert_eq!(dsl.ast().bonds().count(), 1);
    }

    #[rstest]
    fn test_molecule_dsl_from_str_rejects_invalid() {
        let err = "not a map".parse::<MoleculeDsl>().unwrap_err();
        assert!(matches!(err, ParseError::EdnParse(_)));
    }

    #[rstest]
    fn test_molecule_dsl_dsl_to_ast_to_dsl_roundtrip_zeroed() {
        // Round-trip direction: DSL → AST (raise) → DSL (lower) is the
        // identity. AST → DSL → AST isn't, since raising `Undetermined`
        // fields to `Lit(0)` is one-way under `zeroed()`.
        let ast = mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let cfg = MoleculeDefaults::zeroed();
        let raised = dsl.clone().into_ast(&cfg);
        let lowered = MoleculeDsl::from_ast(&raised, &cfg);
        assert_eq!(lowered.ast(), dsl.ast());
    }

    #[rstest]
    fn test_molecule_dsl_from_ast_has_empty_metadata() {
        let ast = mol!(r#"{:atoms ["C"] :bonds []}"#);
        let cfg = MoleculeDefaults::zeroed();
        let dsl = MoleculeDsl::from_ast(&ast, &cfg);
        assert_eq!(dsl.metadata(), &Metadata::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}]}"##)]
    #[case::dative_with_id_and_type(r##"{:atoms ["C" "N"] :bonds [] :dative [{:id :d1 :donor 0 :acceptor 1 :type "1#R"}]}"##)]
    #[case::aromatic_minimal(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic [{:atoms [0 1 2 3 4 5] :type ""}]}"##)]
    #[case::aromatic_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:id :a1 :atoms [0 1] :type "#e6"}]}"##)]
    #[case::aromatic_with_electrons_lits(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic [{:atoms [0 1 2 3 4 5] :electrons [1 1 1 1 1 1] :type ""}]}"##)]
    #[case::aromatic_with_electrons_undetermined(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :electrons [:undetermined :undetermined] :type ""}]}"##)]
    #[case::aromatic_with_electrons_litset(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :electrons [[1 2] 1] :type ""}]}"##)]
    #[case::aromatic_with_electrons_expr(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :electrons ["?n + 1" 1] :type ""}]}"##)]
    #[case::aromatic_with_electrons_and_total(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic [{:atoms [0 1 2 3 4 5] :electrons [1 1 1 1 1 1] :type "#e6"}]}"##)]
    #[case::multicenter_minimal(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter [{:atoms [0 1 2] :type ""}]}"##)]
    #[case::multicenter_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:id :m1 :atoms [0 1] :type "#e2"}]}"##)]
    #[case::multicenter_with_electrons_lits(r##"{:atoms ["B" "H" "B"] :bonds [] :multicenter [{:atoms [0 1 2] :electrons [1 0 1] :type ""}]}"##)]
    #[case::multicenter_with_electrons_undetermined(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :electrons [:undetermined :undetermined] :type ""}]}"##)]
    #[case::multicenter_with_electrons_expr(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :electrons ["?n - 1" 1] :type ""}]}"##)]
    #[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}]}"##)]
    #[case::noncovalent_with_id(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:id :n1 :a 0 :b 1 :type "Hbd"}]}"##)]
    fn test_molecule_dsl_edn_roundtrip_non_localized_entities(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_connected_constraint() {
        let source =
            r##"{:atoms ["C" "C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1 2]}}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
        assert_eq!(dsl.ast().constraints().len(), 1);
    }

    /// Vacuous molecule-level constraints — `ChargeSum` / `BondOrderSum`
    /// with `Undetermined` sum, `SpinSum` with both spin fields
    /// `Undetermined` — are dropped during AST → DSL lowering. The
    /// canonical EDN omits the entire `:constraints` key when the only
    /// entries are vacuous.
    #[rstest]
    fn test_molecule_dsl_render_elides_vacuous_charge_sum() {
        let mut ast = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        ast.constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: ValueAst::Undetermined,
            }));
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else { panic!("expected map") };
        assert!(
            m.get_keyword("constraints").is_none(),
            "vacuous ChargeSum should not surface as :constraints, got {:?}",
            m.get_keyword("constraints"),
        );
    }

    #[rstest]
    fn test_molecule_dsl_render_elides_vacuous_bond_order_sum() {
        let mut ast = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        ast.constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::BondOrderSum {
                bonds: None,
                sum: ValueAst::Undetermined,
            }));
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else { panic!("expected map") };
        assert!(m.get_keyword("constraints").is_none());
    }

    #[rstest]
    fn test_molecule_dsl_render_elides_vacuous_spin_sum() {
        let mut ast = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        ast.constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::SpinSum {
                atoms: None,
                spin: SpinStateAst::default(),
            }));
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else { panic!("expected map") };
        assert!(m.get_keyword("constraints").is_none());
    }

    /// Non-vacuous molecule-level constraints survive the lowering and
    /// vacuous neighbors in the same constraint vec are dropped while
    /// the surviving entry is rendered.
    #[rstest]
    fn test_molecule_dsl_render_keeps_non_vacuous_drops_vacuous() {
        let mut ast = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        ast.constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: ValueAst::Undetermined,
            }));
        ast.constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: ValueAst::Lit(0),
            }));
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else { panic!("expected map") };
        let cs = m.get_keyword("constraints").expect("constraints key present");
        let Edn::Vector(v) = cs else { panic!("expected vec") };
        assert_eq!(v.len(), 1, "only the non-vacuous ChargeSum should survive");
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_connected_all_atoms() {
        let source = r##"{:atoms ["C" "C" "C"] :bonds [] :constraints [{:connected {}}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
        assert_eq!(dsl.ast().constraints().len(), 1);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_bond_order_sum_by_id() {
        let source = r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type :single} {:id :b2 :a 1 :b 2 :type :single}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_atom_leaf_constraint_by_id() {
        let source =
            r##"{:atoms [[:c1 "C"]] :bonds [] :constraints [{:not {:atom [:c1 {:valence 3}]}}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    // region: MoleculeAst direct EDN

    /// `MoleculeAst::to_edn` emits canonical EDN with positional refs only,
    /// regardless of any id keywords on the input. Parsing the canonical
    /// output back yields the same AST.
    #[rstest]
    fn test_molecule_ast_to_edn_canonical_positional_refs() {
        // Input has id keywords on atoms, bonds, and a constraint anchor.
        let source = r##"{:atoms [[:c1 "C"] [:c2 "C"]]
                          :bonds [{:id :b1 :a :c1 :b :c2 :type "1"}]
                          :constraints [{:atom [:c1 {:valence 4}]}
                                        {:bond [:b1 :aromatic]}]}"##;
        let dsl = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
        let (ast, _meta) = dsl.into_parts();

        // Canonical render: positional refs only.
        let canonical_source =
            r##"{:atoms ["C" "C"] :bonds [[0 1 :single]]
                 :constraints [{:atom [0 {:valence 4}]}
                               {:bond [0 :aromatic]}]}"##;
        assert_eq!(ast.to_edn(), read_string(canonical_source).unwrap());
    }

    #[rstest]
    fn test_molecule_ast_from_edn_tree_roundtrip() {
        let source = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;
        let edn = read_string(source).unwrap();
        let ast = MoleculeAst::from_edn(&edn).unwrap();
        assert_eq!(ast.atoms().count(), 2);
        assert_eq!(ast.bonds().count(), 1);
        // Render → parse → equal AST.
        let rendered = ast.to_edn();
        let reparsed = MoleculeAst::from_edn(&rendered).unwrap();
        assert_eq!(ast, reparsed);
    }

    #[rstest]
    fn test_molecule_ast_from_edn_str_fast_path() {
        let source = r##"{:atoms ["C" "O" "H"] :bonds [[0 1 "1"] [1 2 "1"]]}"##;
        let ast = MoleculeAst::from_edn_str(source).unwrap();
        assert_eq!(ast.atoms().count(), 3);
        assert_eq!(ast.bonds().count(), 2);
    }

    #[rstest]
    fn test_molecule_ast_from_edn_drops_id_metadata() {
        // Input carries ids; AST is metadata-free, so reparsing the rendered
        // form (which has no ids) should match the AST from the original parse.
        let source = r##"{:atoms [[:carbon "C"] [:oxygen "O"]]
                          :bonds [{:id :myb :a :carbon :b :oxygen :type "1"}]}"##;
        let ast = MoleculeAst::from_edn(&read_string(source).unwrap()).unwrap();
        let rendered = ast.to_edn().to_string();
        // No user-defined id keywords leaked through. (`:a` / `:b` / `:type`
        // are bond-entry field names, not ids.)
        assert!(!rendered.contains(":carbon"));
        assert!(!rendered.contains(":oxygen"));
        assert!(!rendered.contains(":myb"));
        let reparsed = MoleculeAst::from_edn_str(&rendered).unwrap();
        assert_eq!(ast, reparsed);
    }

    // endregion: MoleculeAst direct EDN

    #[rstest]
    fn test_molecule_dsl_constraint_unknown_ref_errors() {
        let source =
            r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [:nope 0]}}]}"##;
        let edn = read_string(source).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_sub_pattern() {
        let source = r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    fn test_molecule_dsl_sub_pattern_pattern_side_out_of_range_errors() {
        let source = r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 5]]} :pattern {:atoms ["N"] :bonds []}}}]}"##;
        let edn = read_string(source).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    #[case::valence(r##"{:atoms ["C#v4"] :bonds []}"##)]
    #[case::ring_count(r##"{:atoms ["N#R2"] :bonds []}"##)]
    #[case::atom_multiple(r##"{:atoms ["C#v4#R+"] :bonds []}"##)]
    fn test_molecule_dsl_edn_roundtrip_inline_constraints(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    #[case::single(r##"{:atoms ["C" "C"] :bonds [[0 1 :single]]}"##)]
    #[case::double(r##"{:atoms ["C" "C"] :bonds [[0 1 :double]]}"##)]
    #[case::triple(r##"{:atoms ["C" "C"] :bonds [[0 1 :triple]]}"##)]
    #[case::quadruple(r##"{:atoms ["C" "C"] :bonds [[0 1 :quadruple]]}"##)]
    #[case::aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 :aromatic]]}"##)]
    fn test_molecule_dsl_edn_roundtrip_bond_keyword_shorthands(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    // region: Entity endpoint ref errors

    #[rstest]
    fn test_molecule_dsl_bond_endpoint_out_of_range_errors() {
        let edn = read_string(r##"{:atoms ["C" "C"] :bonds [[0 5 "1"]]}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_bond_endpoint_unknown_id_errors() {
        let edn = read_string(r##"{:atoms ["C" "C"] :bonds [[:nope 0 "1"]]}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_noncovalent_endpoint_out_of_range_errors() {
        let edn = read_string(
            r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 99 :type "Hbd"}]}"##,
        )
        .unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_aromatic_atom_out_of_range_errors() {
        let edn =
            read_string(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 5] :type ""}]}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_dative_unknown_donor_id_errors() {
        let edn =
            read_string(
                r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor :nope :acceptor 1 :type :single}]}"##,
            )
            .unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    // endregion: Entity endpoint ref errors

    // region: :type required

    /// `:type` is required on every entry kind that has a DSL surface (bond,
    /// dative, aromatic, multicenter, noncovalent). Missing `:type` is a
    /// `MissingField` error in both the streaming and tree paths.
    #[rustfmt::skip]
    #[rstest]
    #[case::bond_without_type(r##"{:atoms ["C" "C"] :bonds [{:a 0 :b 1}]}"##)]
    #[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1}]}"##)]
    #[case::aromatic_without_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1]}]}"##)]
    #[case::multicenter_without_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1]}]}"##)]
    #[case::noncovalent_without_type(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1}]}"##)]
    fn test_molecule_dsl_type_required_tree(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(
            matches!(err, DeError::MissingField { .. }),
            "expected MissingField, got {:?}",
            err,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::bond_without_type(r##"{:atoms ["C" "C"] :bonds [{:a 0 :b 1}]}"##)]
    #[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1}]}"##)]
    #[case::aromatic_without_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1]}]}"##)]
    #[case::multicenter_without_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1]}]}"##)]
    #[case::noncovalent_without_type(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1}]}"##)]
    fn test_molecule_dsl_type_required_streaming(#[case] source: &str) {
        let err = MoleculeDsl::from_edn_str(source).unwrap_err();
        let de = match err {
            EdnError::De(de) => de,
            other => panic!("expected DeError, got {:?}", other),
        };
        assert!(
            matches!(de, DeError::MissingField { .. }),
            "expected MissingField, got {:?}",
            de,
        );
    }

    // endregion: :type required

    // region: :guards reserved-future key

    // endregion: :guards reserved-future key

    // region: Per-variant streaming parity

    #[rustfmt::skip]
    #[rstest]
    // AtomConstraint variants via :constraints [{:atom [0 <form>]}]
    #[case::atom_valence_lit(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence 4}]}]}"##)]
    #[case::atom_valence_set(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence [3 4]}]}]}"##)]
    #[case::atom_valence_undetermined(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence :undetermined}]}]}"##)]
    #[case::atom_valence_expr(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence "?h >= 1"}]}]}"##)]
    #[case::atom_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:degree 3}]}]}"##)]
    #[case::atom_connectivity(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:connectivity 4}]}]}"##)]
    #[case::atom_ring_connectivity(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-connectivity 2}]}]}"##)]
    #[case::atom_total_hydrogens(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-hydrogens 3}]}]}"##)]
    #[case::atom_ring_count(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-count 1}]}]}"##)]
    #[case::atom_ring_size(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-size 6}]}]}"##)]
    #[case::atom_donated_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:donated-pairs 1}]}]}"##)]
    #[case::atom_accepted_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:accepted-pairs 2}]}]}"##)]
    #[case::atom_aromatic_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence :not-aromatic}]}]}"##)]
    #[case::atom_aromatic_valence_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence {:aromatic 6}}]}]}"##)]
    #[case::atom_multicenter_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence :not-multicenter}]}]}"##)]
    #[case::atom_multicenter_valence_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence {:multicenter 3}}]}]}"##)]
    // BondConstraint variants
    #[case::bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 :aromatic]}]}"##)]
    #[case::bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-count 1}]}]}"##)]
    #[case::bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-size 6}]}]}"##)]
    // DativeBondConstraint variants
    #[case::dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-count 1}]}]}"##)]
    #[case::dative_ring_size(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-size 5}]}]}"##)]
    #[case::dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
    #[case::dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
    #[case::dative_donor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-donor-satisfies [0 {:valence 4}]}]}"##)]
    #[case::dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor-satisfies [0 {:degree 2}]}]}"##)]
    #[case::dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
    // RelationalConstraint variants for aromatic system
    #[case::aromatic_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-atoms [0 [0 1]]}]}"##)]
    #[case::aromatic_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
    #[case::aromatic_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
    #[case::aromatic_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
    #[case::aromatic_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-any-atom [0 {:degree 3}]}]}"##)]
    // RelationalConstraint variants for multicenter bond
    #[case::multicenter_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-atoms [0 [0 1]]}]}"##)]
    #[case::multicenter_contains(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
    #[case::multicenter_contains_all(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
    #[case::multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
    #[case::multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-any-atom [0 {:degree 3}]}]}"##)]
    // RelationalConstraint variants for noncovalent bond
    #[case::noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
    #[case::noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
    #[case::noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 4} {:valence 1}]]}]}"##)]
    // MoleculeConstraint variants (via flattened keys)
    #[case::molecule_charge_sum(r##"{:atoms ["C" "N"] :bonds [] :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"##)]
    #[case::molecule_spin_sum(r##"{:atoms ["C"] :bonds [] :constraints [{:spin-sum {:atoms [0] :spin {:unpaired 1 :multiplicity 2}}}]}"##)]
    // Anchor with multiple entity kinds (exercises all 6 ref-pair readers)
    #[case::sub_pattern_anchor_bonds_and_atoms(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]] :bonds [[0 0]]} :pattern {:atoms ["N" "N"] :bonds [[0 1 "1"]]}}}]}"##)]
    fn test_molecule_dsl_streaming_per_variant_parity(#[case] source: &str) {
        let via_str = MoleculeDsl::from_edn_str(source).unwrap();
        let tree = read_string(source).unwrap();
        let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
        assert_eq!(via_str, via_tree);
    }

    // endregion: Per-variant streaming parity

    // region: Streaming vs tree error parity

    #[rstest]
    #[case::missing_key_value(r##"{:atoms}"##)]
    #[case::string_key(r##"{"atoms" []}"##)]
    #[case::truncated_map(r##"{:atoms ["C""##)]
    #[case::truncated_outer(r##"{:atoms []"##)]
    #[case::unknown_top_key(r##"{:atoms [] :bonds [] :bogus 1}"##)]
    #[case::atom_out_of_range_in_bond(r##"{:atoms ["C"] :bonds [[0 99 "1"]]}"##)]
    #[case::unknown_constraint_key(r##"{:atoms ["C"] :bonds [] :constraints [{:bogus 1}]}"##)]
    #[case::unknown_atom_constraint_kind(
        r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:bogus 1}]}]}"##
    )]
    fn test_molecule_dsl_streaming_error_parity(#[case] source: &str) {
        let via_str_err = MoleculeDsl::from_edn_str(source).is_err();
        let via_tree_err = read_string(source)
            .map_err(|_| ())
            .and_then(|edn| MoleculeDsl::from_edn(&edn).map_err(|_| ()))
            .is_err();
        assert!(
            via_str_err,
            "{source:?}: streaming path should have errored"
        );
        assert!(via_tree_err, "{source:?}: tree path should have errored");
    }

    // endregion: Streaming vs tree error parity

    // region: Cross-entity id disjointness

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_vs_bond(r##"{:atoms [[:x "C"] [:y "C"]] :bonds [{:id :x :a 0 :b 1 :type "1"}]}"##)]
    #[case::atom_vs_alias(r##"{:atoms [[:x "C"]] :bonds [] :atom-aliases [:x "N"]}"##)]
    #[case::bond_vs_dative(r##"{:atoms ["C" "N"] :bonds [{:id :x :a 0 :b 1 :type "1"}] :dative [{:id :x :donor 0 :acceptor 1 :type :single}]}"##)]
    #[case::atom_vs_aromatic(r##"{:atoms [[:x "C"] [:y "C"]] :bonds [] :aromatic [{:id :x :atoms [0 1] :type ""}]}"##)]
    #[case::bond_vs_noncovalent(r##"{:atoms ["C" "C"] :bonds [{:id :x :a 0 :b 1 :type "1"}] :noncovalent [{:id :x :a 0 :b 1 :type "Hbd"}]}"##)]
    fn test_molecule_dsl_cross_entity_id_collision_errors(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    // endregion: Cross-entity id disjointness

    // region: Alias bijectivity

    #[rstest]
    fn test_molecule_dsl_duplicate_alias_name_errors() {
        let edn = read_string(r##"{:atoms [] :bonds [] :atom-aliases [:a "C" :a "N"]}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_aliases_must_be_bijective() {
        let edn = read_string(r##"{:atoms [] :bonds [] :atom-aliases [:a "C" :b "C"]}"##).unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_guards_key_accepted_and_ignored() {
        let source = r##"{:atoms ["C"] :bonds [] :guards [[:placeholder]]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        // :guards is silently accepted; the rendered form drops it since the
        // AST has no slot for it yet.
        let rendered = dsl.to_edn();
        let Edn::Map(m) = &rendered else {
            panic!("expected map");
        };
        assert!(m.get_keyword("guards").is_none());
        assert_eq!(dsl.ast().atoms().count(), 1);
    }
    // endregion: Alias bijectivity

    // region: MoleculeAst symmetric I/O

    #[rstest]
    fn test_molecule_ast_from_str_to_string_roundtrip() {
        let s = r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##;
        let ast: MoleculeAst = s.parse().unwrap();
        let rendered = ast.to_string();
        let back: MoleculeAst = rendered.parse().unwrap();
        assert_eq!(back, ast);
    }

    #[rstest]
    fn test_molecule_ast_to_edn_roundtrip() {
        let s = r##"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"##;
        let ast: MoleculeAst = s.parse().unwrap();
        let edn = ast.to_edn();
        let back = MoleculeAst::from_edn(&edn).unwrap();
        assert_eq!(back, ast);
    }

    // endregion: MoleculeAst symmetric I/O
}
