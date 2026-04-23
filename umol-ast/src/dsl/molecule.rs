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
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{
    DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn,
};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::constraint::{
    eof_err, missing, read_atom_ref, read_constraints_dsl, read_map, read_vec,
    unexpected_byte_kind, AtomRef, ConstraintDsl, ConstraintsDsl, ResolveContext,
};
use super::dative::DativeBondDsl;
use super::error::ParseError;
use super::multicenter::MulticenterBondDsl;
use super::noncovalent::NoncovalentBondDsl;
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
use crate::dsl::config::MoleculeDefaults;

/// Surface-form metadata paired with a `MoleculeAst`. Records atom ids,
/// per-entity ids, and the atom-alias table. Never drifts: rewrapped
/// atomically through `MoleculeDsl::from_parts`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_ids: IndexMap<AtomIdx, String>,
    pub atom_aliases: BiMap<String, AtomDsl>,
    pub bond_ids: IndexMap<BondIdx, String>,
    pub dative_bond_ids: IndexMap<DativeBondIdx, String>,
    pub aromatic_system_ids: IndexMap<AromaticSystemIdx, String>,
    pub multicenter_bond_ids: IndexMap<MulticenterBondIdx, String>,
    pub noncovalent_bond_ids: IndexMap<NoncovalentBondIdx, String>,
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
        let (ast, metadata) = mi
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(ast, metadata))
    }
}

// -- Streaming parser -------------------
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
            Ok(AtomEntryInput {
                id: Some(id),
                spec,
            })
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
    let text = de.read_string_or_keyword()?;
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("bond", e).into())
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
            "type" => {
                let s = de.read_string()?;
                bond = Some(
                    s.as_ref()
                        .parse::<DativeBondDsl>()
                        .map_err(|e| DeError::subgrammar("dative", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(DativeBondEntryInput {
        id,
        donor: donor.ok_or_else(|| missing("donor", "dative-bond-entry"))?,
        acceptor: acceptor.ok_or_else(|| missing("acceptor", "dative-bond-entry"))?,
        bond: bond.unwrap_or_default(),
    })
}

fn read_aromatic_system_entry(
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
    Ok(AromaticSystemEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "aromatic-system-entry"))?,
        system: system.unwrap_or_default(),
    })
}

fn read_multicenter_bond_entry(
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
    Ok(MulticenterBondEntryInput {
        id,
        atoms: atoms.ok_or_else(|| missing("atoms", "multicenter-bond-entry"))?,
        bond: bond.unwrap_or_default(),
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
) -> Result<Vec<(String, AtomDsl)>, EdnError> {
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
        out.push((name, dsl));
    }
    Ok(out)
}


impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(&self.ast, &self.metadata)
    }
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    type Ctx<'a> = MoleculeDefaults;
    type Error = ParseError;

    fn from_ast<'a>(ast: &MoleculeAst, cfg: &Self::Ctx<'a>) -> Result<Self, ParseError> {
        let mut ast_out = ast.clone();
        for atom in ast_out.atoms_mut() {
            *atom = AtomDsl::from_ast(atom, &cfg.atom)?.0;
        }
        for bond in ast_out.bonds_mut() {
            *bond = BondDsl::from_ast(bond, &cfg.bond)?.0;
        }
        for system in ast_out.aromatic_systems_mut() {
            *system = AromaticSystemDsl::from_ast(system, &cfg.aromatic_system)?.0;
        }
        for bond in ast_out.multicenter_bonds_mut() {
            *bond = MulticenterBondDsl::from_ast(bond, &cfg.multicenter_bond)?.0;
        }
        // `DativeBondDsl` and `NoncovalentBondDsl` use unit-shaped defaults
        // (empty struct), so there is nothing to strip here.
        Ok(MoleculeDsl {
            ast: ast_out,
            metadata: Metadata::default(),
        })
    }
}

impl IntoAst<MoleculeAst> for MoleculeDsl {
    type Ctx<'a> = MoleculeDefaults;
    type Error = ParseError;

    fn into_ast<'a>(self, cfg: &Self::Ctx<'a>) -> Result<MoleculeAst, ParseError> {
        let mut ast = self.ast;
        for atom in ast.atoms_mut() {
            *atom = AtomDsl(atom.clone()).into_ast(&cfg.atom)?;
        }
        for bond in ast.bonds_mut() {
            *bond = BondDsl(bond.clone()).into_ast(&cfg.bond)?;
        }
        for system in ast.aromatic_systems_mut() {
            *system = AromaticSystemDsl(system.clone()).into_ast(&cfg.aromatic_system)?;
        }
        for bond in ast.multicenter_bonds_mut() {
            *bond = MulticenterBondDsl(bond.clone()).into_ast(&cfg.multicenter_bond)?;
        }
        Ok(ast)
    }
}

// -- Render --------------------

fn render_molecule_edn(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(8);
    map.insert(Edn::keyword("atoms"), render_atoms(ast, meta));
    map.insert(Edn::keyword("bonds"), render_bonds(ast, meta));
    if ast.dative_bond_count() > 0 {
        map.insert(Edn::keyword("dative"), render_dative(ast, meta));
    }
    if ast.aromatic_system_count() > 0 {
        map.insert(Edn::keyword("aromatic"), render_aromatic(ast, meta));
    }
    if ast.multicenter_bond_count() > 0 {
        map.insert(Edn::keyword("multicenter"), render_multicenter(ast, meta));
    }
    if ast.noncovalent_bond_count() > 0 {
        map.insert(Edn::keyword("noncovalent"), render_noncovalent(ast, meta));
    }
    if !meta.atom_aliases.is_empty() {
        map.insert(Edn::keyword("atom-aliases"), render_atom_aliases(meta));
    }
    if !ast.constraints().is_empty() {
        let ctx = ResolveContext::for_rendering(meta);
        let dsl = ConstraintsDsl::from_ast(ast.constraints(), &ctx)
            .expect("ConstraintsDsl::from_ast is infallible for a well-formed AST");
        map.insert(Edn::keyword("constraints"), dsl.to_edn());
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
    let dsl = AtomDsl(atom.clone());
    let spec = if let Some(alias) = meta.atom_aliases.get_by_right(&dsl) {
        Edn::Keyword(EdnKeyword::owned(alias.clone()))
    } else {
        dsl.to_edn()
    };
    match meta.atom_ids.get(&idx) {
        Some(id) => Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(id.clone())), spec].into()),
        None => spec,
    }
}

fn render_atom_ref(idx: AtomIdx, meta: &Metadata) -> Edn<'static> {
    match meta.atom_ids.get(&idx) {
        Some(id) => Edn::Keyword(EdnKeyword::owned(id.clone())),
        None => Edn::Int(idx.index() as i64),
    }
}

fn render_bonds(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .map(|view| {
            let bond_edn = BondDsl(view.data.clone()).to_edn();
            let a = render_atom_ref(view.src, meta);
            let b = render_atom_ref(view.tgt, meta);
            match meta.bond_ids.get(&view.idx) {
                Some(id) => {
                    let mut m = EdnMap::with_capacity(4);
                    m.insert(
                        Edn::keyword("id"),
                        Edn::Keyword(EdnKeyword::owned(id.clone())),
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
            if let Some(id) = meta.dative_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            m.insert(Edn::keyword("donor"), render_atom_ref(view.donor, meta));
            m.insert(
                Edn::keyword("acceptor"),
                render_atom_ref(view.acceptor, meta),
            );
            let type_str = DativeBondDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
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
            let mut m = EdnMap::with_capacity(3);
            if let Some(id) = meta.aromatic_system_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            let atoms: Vec<Edn<'static>> = view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            let type_str = AromaticSystemDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
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
            let mut m = EdnMap::with_capacity(3);
            if let Some(id) = meta.multicenter_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            let atoms: Vec<Edn<'static>> = view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            let type_str = MulticenterBondDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_noncovalent(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.noncovalent_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            m.insert(Edn::keyword("a"), render_atom_ref(view.atoms[0], meta));
            m.insert(Edn::keyword("b"), render_atom_ref(view.atoms[1], meta));
            m.insert(
                Edn::keyword("type"),
                NoncovalentBondDsl(view.data.clone()).to_edn(),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_atom_aliases(meta: &Metadata) -> Edn<'static> {
    let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases.len() * 2);
    for (name, dsl) in meta.atom_aliases.iter() {
        pairs.push(Edn::Keyword(EdnKeyword::owned(name.clone())));
        pairs.push(dsl.to_edn());
    }
    Edn::Vector(pairs.into())
}

// -- Private parse intermediate ---------------------------------------------
//
// Unresolved, owned-by-value tree that mirrors the EDN shape. Atom entries and
// per-bond endpoints carry `AtomRef` (index or id); constraint leaves carry
// typed per-entity `Constraint*` variants already parsed from their EDN form.
// Lowered destructively via `into_ast(self, cfg)` so that allocations move
// into the final `MoleculeAst`.

#[derive(Clone, Debug, Default, PartialEq)]
#[allow(dead_code)]
pub(crate) struct MoleculeInput {
    pub(crate) atoms: Vec<AtomEntryInput>,
    pub(crate) bonds: Vec<BondEntryInput>,
    pub(crate) dative_bonds: Vec<DativeBondEntryInput>,
    pub(crate) aromatic_systems: Vec<AromaticSystemEntryInput>,
    pub(crate) multicenter_bonds: Vec<MulticenterBondEntryInput>,
    pub(crate) noncovalent_bonds: Vec<NoncovalentBondEntryInput>,
    pub(crate) atom_aliases: Vec<(String, AtomDsl)>,
    pub(crate) constraints: Vec<ConstraintDsl>,
}

/// Atom entry in a parsed molecule map. Mirrors the DSL spec §4 grammar
/// `atom-entry ::= atom-spec | [ keyword atom-spec ]`.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AtomEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) spec: AtomSpecInput,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomSpecInput {
    Bare(Box<AtomDsl>),
    Alias(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) donor: AtomRef,
    pub(crate) acceptor: AtomRef,
    pub(crate) bond: DativeBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AromaticSystemEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) system: AromaticSystemDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct MulticenterBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) bond: MulticenterBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
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

        // Alias table: check no duplicate alias names.
        let mut alias_table: IndexMap<String, AtomDsl> = IndexMap::new();
        for (name, dsl) in alias_entries {
            if alias_table.contains_key(&name) {
                return Err(ParseError::DuplicateId(name));
            }
            alias_table.insert(name, dsl);
        }

        // Atoms: materialize AtomAst from each entry; collect ids.
        let mut atoms: Vec<AtomAst> = Vec::with_capacity(atom_entries.len());
        let mut atom_ids: IndexMap<AtomIdx, String> = IndexMap::new();
        let mut atom_id_to_idx: IndexMap<String, AtomIdx> = IndexMap::new();
        for (pos, entry) in atom_entries.into_iter().enumerate() {
            let idx = AtomIdx(pos as u32);
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                atom_id_to_idx.insert(id.clone(), idx);
                atom_ids.insert(idx, id);
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
        let mut bond_ids: IndexMap<BondIdx, String> = IndexMap::new();
        let mut entry_ids: IndexMap<String, ()> = IndexMap::new();
        for (pos, entry) in bond_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                bond_ids.insert(BondIdx(pos as u32), id);
            }
            let a = entry
                .a
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            let b = entry
                .b
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            bonds.push((a, b, entry.bond.0));
        }

        // Dative bonds.
        let mut dative_list: Vec<(AtomIdx, AtomIdx, DativeBondAst)> =
            Vec::with_capacity(dative_entries.len());
        let mut dative_bond_ids: IndexMap<DativeBondIdx, String> = IndexMap::new();
        for (pos, entry) in dative_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                dative_bond_ids.insert(DativeBondIdx(pos as u32), id);
            }
            let donor = entry
                .donor
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            let acceptor = entry
                .acceptor
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            dative_list.push((donor, acceptor, entry.bond.0));
        }

        // Aromatic systems.
        let mut aromatic_list: Vec<(Vec<AtomIdx>, AromaticSystemAst)> =
            Vec::with_capacity(aromatic_entries.len());
        let mut aromatic_system_ids: IndexMap<AromaticSystemIdx, String> = IndexMap::new();
        for (pos, entry) in aromatic_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                aromatic_system_ids.insert(AromaticSystemIdx(pos as u32), id);
            }
            let atoms_resolved: Vec<AtomIdx> = entry
                .atoms
                .into_iter()
                .map(|r| r.into_ast(atom_count, &atom_only_metadata(&atom_ids)))
                .collect::<Result<_, _>>()?;
            aromatic_list.push((atoms_resolved, entry.system.0));
        }

        // Multicenter bonds.
        let mut multicenter_list: Vec<(Vec<AtomIdx>, MulticenterBondAst)> =
            Vec::with_capacity(multicenter_entries.len());
        let mut multicenter_bond_ids: IndexMap<MulticenterBondIdx, String> = IndexMap::new();
        for (pos, entry) in multicenter_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                multicenter_bond_ids.insert(MulticenterBondIdx(pos as u32), id);
            }
            let atoms_resolved: Vec<AtomIdx> = entry
                .atoms
                .into_iter()
                .map(|r| r.into_ast(atom_count, &atom_only_metadata(&atom_ids)))
                .collect::<Result<_, _>>()?;
            multicenter_list.push((atoms_resolved, entry.bond.0));
        }

        // Noncovalent bonds.
        let mut noncovalent_list: Vec<(AtomIdx, AtomIdx, NoncovalentBondAst)> =
            Vec::with_capacity(noncovalent_entries.len());
        let mut noncovalent_bond_ids: IndexMap<NoncovalentBondIdx, String> = IndexMap::new();
        for (pos, entry) in noncovalent_entries.into_iter().enumerate() {
            if let Some(id) = entry.id {
                check_id_disjoint(&id, &atom_id_to_idx, &alias_table)?;
                if entry_ids.insert(id.clone(), ()).is_some() {
                    return Err(ParseError::DuplicateId(id));
                }
                noncovalent_bond_ids.insert(NoncovalentBondIdx(pos as u32), id);
            }
            let a = entry
                .a
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            let b = entry
                .b
                .into_ast(atom_count, &atom_only_metadata(&atom_ids))?;
            noncovalent_list.push((a, b, entry.bond.0));
        }

        // Assemble alias bimap.
        let mut atom_aliases: BiMap<String, AtomDsl> = BiMap::new();
        for (name, dsl) in alias_table {
            if atom_aliases.insert_no_overwrite(name, dsl).is_err() {
                return Err(ParseError::InvalidValue(
                    "atom-aliases must be bijective".into(),
                ));
            }
        }

        let metadata = Metadata {
            atom_ids,
            atom_aliases,
            bond_ids,
            dative_bond_ids,
            aromatic_system_ids,
            multicenter_bond_ids,
            noncovalent_bond_ids,
        };

        // Resolve constraint refs against the final metadata + counts.
        let ctx = ResolveContext {
            atom_count,
            bond_count: bonds.len(),
            dative_bond_count: dative_list.len(),
            aromatic_system_count: aromatic_list.len(),
            multicenter_bond_count: multicenter_list.len(),
            noncovalent_bond_count: noncovalent_list.len(),
            metadata: &metadata,
        };
        let constraints = ConstraintsDsl(constraint_dsls).into_ast(&ctx)?;

        let ast = MoleculeAst::new(
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

// -- Parse --------------------

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
                input.constraints =
                    parse_vec(v, ":constraints", |e| ConstraintDsl::from_edn(e))?
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
        bond: match m.get_keyword("type") {
            Some(edn) => DativeBondDsl::from_edn(edn)?,
            None => DativeBondDsl::default(),
        },
    })
}

fn parse_aromatic_system_entry(edn: &Edn<'_>) -> Result<AromaticSystemEntryInput, DeError> {
    let m = expect_map(edn, "aromatic-system-entry")?;
    Ok(AromaticSystemEntryInput {
        id: optional_id(m)?,
        atoms: parse_vec(
            required_key(m, "atoms", "aromatic-system-entry")?,
            ":atoms",
            |e| AtomRef::from_edn(e),
        )?,
        system: match m.get_keyword("type") {
            Some(edn) => AromaticSystemDsl::from_edn(edn)?,
            None => AromaticSystemDsl::default(),
        },
    })
}

fn parse_multicenter_bond_entry(edn: &Edn<'_>) -> Result<MulticenterBondEntryInput, DeError> {
    let m = expect_map(edn, "multicenter-bond-entry")?;
    Ok(MulticenterBondEntryInput {
        id: optional_id(m)?,
        atoms: parse_vec(
            required_key(m, "atoms", "multicenter-bond-entry")?,
            ":atoms",
            |e| AtomRef::from_edn(e),
        )?,
        bond: match m.get_keyword("type") {
            Some(edn) => MulticenterBondDsl::from_edn(edn)?,
            None => MulticenterBondDsl::default(),
        },
    })
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

fn parse_atom_aliases(edn: &Edn<'_>) -> Result<Vec<(String, AtomDsl)>, DeError> {
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
        out.push((name.name().to_string(), dsl));
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
    alias_table: &IndexMap<String, AtomDsl>,
) -> Result<(), ParseError> {
    if atom_id_to_idx.contains_key(id) || alias_table.contains_key(id) {
        return Err(ParseError::DuplicateId(id.to_string()));
    }
    Ok(())
}

/// Minimal metadata carrying only the atom ids, used for per-entity ref
/// resolution during the first pass (before the full metadata exists). All
/// other ref kinds are zero-sized at this point.
fn atom_only_metadata(atom_ids: &IndexMap<AtomIdx, String>) -> Metadata {
    Metadata {
        atom_ids: atom_ids.clone(),
        atom_aliases: BiMap::new(),
        bond_ids: IndexMap::new(),
        dative_bond_ids: IndexMap::new(),
        aromatic_system_ids: IndexMap::new(),
        multicenter_bond_ids: IndexMap::new(),
        noncovalent_bond_ids: IndexMap::new(),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::{AtomAst, ElementAst};
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::value::ValueAst;

    fn c_atom() -> AtomAst {
        AtomAst::new(ElementAst::Lit(Element::C))
    }

    fn single_bond() -> BondAst {
        BondAst::new(ValueAst::Lit(1))
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
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_with_id() {
        let mut atom_ids = IndexMap::new();
        atom_ids.insert(AtomIdx(0), "c1".to_string());
        let meta = Metadata {
            atom_ids,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_bond_with_id_uses_map_form() {
        let mut bond_ids = IndexMap::new();
        bond_ids.insert(BondIdx(0), "b1".to_string());
        let meta = Metadata {
            bond_ids,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"}]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_alias_substituted() {
        let mut atom_aliases: BiMap<String, AtomDsl> = BiMap::new();
        atom_aliases.insert("x".into(), AtomDsl(c_atom()));
        let meta = Metadata {
            atom_aliases,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
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
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        assert_eq!(dsl.to_string(), dsl.to_edn().to_string());
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_omits_empty_optional_sections() {
        let ast = MoleculeAst::new(
            vec![c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
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
        assert_eq!(dsl.ast().atom_count(), 0);
        assert_eq!(dsl.ast().bond_count(), 0);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_two_atoms_one_bond() {
        let edn = read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().atom_count(), 2);
        assert_eq!(dsl.ast().bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_atom_with_id() {
        let edn = read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.metadata().atom_ids.get(&AtomIdx(0)), Some(&"c1".into()));
        assert_eq!(dsl.metadata().atom_ids.get(&AtomIdx(1)), None);
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_bond_map_form_with_id() {
        let edn =
            read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"}]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().bond_count(), 1);
        assert_eq!(dsl.metadata().bond_ids.get(&BondIdx(0)), Some(&"b1".into()));
    }

    #[rstest]
    fn test_molecule_dsl_from_edn_atom_aliases() {
        let edn = read_string(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.ast().atom_count(), 2);
        assert!(dsl.metadata().atom_aliases.contains_left("x"));
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
        let source = r##"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        let rendered = dsl.to_edn();
        assert_eq!(rendered, edn);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_with_ids_and_aliases() {
        let source = r##"{:atoms [[:a "C"] [:b "C"] :x] :bonds [{:id :b1 :a :a :b :b :type "1"} [:b 2 "2"]] :atom-aliases [:x "N"]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        let rendered = dsl.to_edn();
        assert_eq!(rendered, edn);
    }

    #[rstest]
    #[case::empty(r##"{:atoms [] :bonds []}"##)]
    #[case::small(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##)]
    #[case::with_ids(r##"{:atoms [[:a "C"] [:b "N"]] :bonds [{:id :b1 :a :a :b :b :type "1"}]}"##)]
    fn test_molecule_dsl_from_edn_str_matches_from_edn(#[case] source: &str) {
        let via_str = MoleculeDsl::from_edn_str(source).unwrap();
        let tree = read_string(source).unwrap();
        let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
        assert_eq!(via_str, via_tree);
    }

    #[rstest]
    fn test_molecule_dsl_from_str_parses_edn_source() {
        let source = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;
        let dsl: MoleculeDsl = source.parse().unwrap();
        assert_eq!(dsl.ast().atom_count(), 2);
        assert_eq!(dsl.ast().bond_count(), 1);
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
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let cfg = MoleculeDefaults::zeroed();
        let raised = dsl.clone().into_ast(&cfg).unwrap();
        let lowered = MoleculeDsl::from_ast(&raised, &cfg).unwrap();
        assert_eq!(lowered.ast(), dsl.ast());
    }

    #[rstest]
    fn test_molecule_dsl_from_ast_has_empty_metadata() {
        let ast = MoleculeAst::new(
            vec![c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let cfg = MoleculeDefaults::zeroed();
        let dsl = MoleculeDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.metadata(), &Metadata::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1}]}"##)]
    #[case::dative_with_id_and_type(r##"{:atoms ["C" "N"] :bonds [] :dative [{:id :d1 :donor 0 :acceptor 1 :type "#R"}]}"##)]
    #[case::aromatic_minimal(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic [{:atoms [0 1 2 3 4 5]}]}"##)]
    #[case::aromatic_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:id :a1 :atoms [0 1] :type "#e6"}]}"##)]
    #[case::multicenter_minimal(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter [{:atoms [0 1 2]}]}"##)]
    #[case::multicenter_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:id :m1 :atoms [0 1] :type "#e2"}]}"##)]
    #[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1 :type "Hbd"}]}"##)]
    #[case::noncovalent_with_id(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:id :n1 :a 0 :b 1 :type "Hbd"}]}"##)]
    fn test_molecule_dsl_edn_roundtrip_non_localized_entities(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_connected_constraint() {
        let source = r##"{:atoms ["C" "C" "C"] :bonds [] :constraints [{:connected [0 1 2]}]}"##;
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
        assert_eq!(dsl.ast().constraints().len(), 1);
    }

    #[rstest]
    fn test_molecule_dsl_edn_roundtrip_bond_order_sum_by_id() {
        let source = r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"} {:id :b2 :a 1 :b 2 :type "1"}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##;
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

    #[rstest]
    fn test_molecule_dsl_constraint_unknown_ref_errors() {
        let source = r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected [:nope 0]}]}"##;
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
    #[case::bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"##)]
    fn test_molecule_dsl_edn_roundtrip_inline_constraints(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    // -- Entity endpoint ref errors ----------------

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
        let edn = read_string(r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 5]}]}"##)
            .unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_molecule_dsl_dative_unknown_donor_id_errors() {
        let edn = read_string(
            r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor :nope :acceptor 1}]}"##,
        )
        .unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    // -- :type optionality ----------------

    #[rstest]
    #[case::dative_without_type(
        r##"{:atoms ["C" "N"] :bonds [] :dative [{:donor 0 :acceptor 1}]}"##
    )]
    #[case::aromatic_without_type(
        r##"{:atoms ["C" "C"] :bonds [] :aromatic [{:atoms [0 1]}]}"##
    )]
    #[case::multicenter_without_type(
        r##"{:atoms ["C" "C"] :bonds [] :multicenter [{:atoms [0 1]}]}"##
    )]
    fn test_molecule_dsl_type_field_optional(#[case] source: &str) {
        let edn = read_string(source).unwrap();
        let dsl = MoleculeDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.to_edn(), edn);
    }

    #[rstest]
    fn test_molecule_dsl_noncovalent_type_is_required() {
        let edn = read_string(r##"{:atoms ["N" "H"] :bonds [] :noncovalent [{:a 0 :b 1}]}"##)
            .unwrap();
        let err = MoleculeDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::MissingField { .. }));
    }

    // -- :guards reserved-future key ----------------

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
        assert_eq!(dsl.ast().atom_count(), 1);
    }
}
