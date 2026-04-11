//! Molecule DSL definitions

use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use umol_shared::SpinState;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::atom::parse_atom_dsl;
use super::error::ParseError;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::molecule::{
    AromaticSystem, DativeBond, LocalizedBond, MoleculeAst, MulticenterBond, NoncovalentBond,
};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_tags: IndexMap<usize, String>,
    pub atom_aliases: BiMap<String, AtomAst>,
    pub bond_ids: IndexMap<usize, String>,
    pub dative_bond_ids: IndexMap<usize, String>,
    pub aromatic_system_ids: IndexMap<usize, String>,
    pub multicenter_bond_ids: IndexMap<usize, String>,
    pub noncovalent_bond_ids: IndexMap<usize, String>,
}

/// Owns a `MoleculeAst` together with the `Metadata` bound to it. Fields are
/// private so that metadata cannot drift onto a different AST: once paired,
/// the pair is either rewrapped atomically or taken apart via `into_parts`.
/// This is the only type with `FromEdn`/`ToEdn` impls for the molecule DSL.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeAstWrapper {
    ast: MoleculeAst,
    metadata: Metadata,
}

impl MoleculeAstWrapper {
    pub fn new(ast: MoleculeAst, metadata: Metadata) -> Self {
        Self { ast, metadata }
    }

    pub fn from_ast(ast: MoleculeAst) -> Self {
        Self {
            ast,
            metadata: Metadata::default(),
        }
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

/// Parse a molecule DSL EDN string via the fused single-pass parser.
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeAstWrapper, ParseError> {
    let mut de = EdnStreamDeserializer::new(input);
    let mol_input =
        read_molecule_input(&mut de).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    de.expect_eof()
        .map_err(|e| ParseError::EdnParse(e.to_string()))?;
    mol_input.into_dsl()
}

impl FromStr for MoleculeAstWrapper {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_molecule_dsl(s)
    }
}

#[derive(Clone, Debug)]
enum AtomEntryInput {
    Str(String),
    Tagged(String, Box<AtomEntryInput>),
}

#[derive(Clone, Debug)]
enum AtomRefInput {
    Index(usize),
    Tag(String),
}

#[derive(Clone, Debug)]
struct LocalizedBondInput {
    id: Option<String>,
    a: AtomRefInput,
    b: AtomRefInput,
    bond: BondAst,
}

#[derive(Clone, Debug)]
struct DativeBondInput {
    id: Option<String>,
    donor: AtomRefInput,
    acceptor: AtomRefInput,
    bond: BondAst,
}

#[derive(Clone, Debug)]
struct AromaticSystemInput {
    id: Option<String>,
    atoms: Vec<AtomRefInput>,
}

#[derive(Clone, Debug)]
struct MulticenterBondInput {
    id: Option<String>,
    atoms: Vec<AtomRefInput>,
}

#[derive(Clone, Debug)]
struct NoncovalentBondInput {
    id: Option<String>,
    a: AtomRefInput,
    b: AtomRefInput,
    bond: BondAst,
}

struct RawMoleculeAst {
    atoms: Vec<AtomEntryInput>,
    bonds: Vec<LocalizedBondInput>,
    dative_bonds: Vec<DativeBondInput>,
    aromatic_systems: Vec<AromaticSystemInput>,
    multicenter_bonds: Vec<MulticenterBondInput>,
    noncovalent_bonds: Vec<NoncovalentBondInput>,
    atom_aliases: Vec<String>,
    charge: Option<i64>,
    spin: Option<SpinState>,
}

impl RawMoleculeAst {
    fn into_dsl(self) -> Result<MoleculeAstWrapper, ParseError> {
        let alias_table = Self::build_alias_table(&self.atom_aliases)?;

        let mut atoms: Vec<AtomAst> = Vec::with_capacity(self.atoms.len());
        let mut atom_tags: IndexMap<usize, String> = IndexMap::new();
        let mut tag_to_index: IndexMap<String, usize> = IndexMap::new();
        let mut atom_aliases: BiMap<String, AtomAst> = BiMap::new();

        for entry in self.atoms {
            let pos = atoms.len();
            let (tag, atom_str) = Self::resolve_entry(entry, &alias_table)?;
            let atom_ast = parse_atom_dsl(&atom_str)?;
            if let Some(tag_name) = tag {
                if tag_to_index.contains_key(&tag_name) || alias_table.contains_key(&tag_name) {
                    return Err(ParseError::DuplicateId(tag_name));
                }
                tag_to_index.insert(tag_name.clone(), pos);
                atom_tags.insert(pos, tag_name);
            }
            atoms.push(atom_ast);
        }

        for (name, def) in &alias_table {
            let atom_ast = parse_atom_dsl(def)?;
            atom_aliases.insert(name.clone(), atom_ast);
        }

        let resolve = |r: &AtomRefInput| -> Result<usize, ParseError> {
            match r {
                AtomRefInput::Index(i) => {
                    if *i < atoms.len() {
                        Ok(*i)
                    } else {
                        Err(ParseError::InvalidAtomIndex(i.to_string()))
                    }
                }
                AtomRefInput::Tag(name) => tag_to_index
                    .get(name)
                    .copied()
                    .ok_or_else(|| ParseError::InvalidAtomIndex(name.clone())),
            }
        };

        let mut seen_ids: HashSet<String> = HashSet::new();
        let mut check_id = |id: Option<String>| -> Result<Option<String>, ParseError> {
            match id {
                Some(s) => {
                    if !seen_ids.insert(s.clone()) {
                        return Err(ParseError::DuplicateId(s));
                    }
                    Ok(Some(s))
                }
                None => Ok(None),
            }
        };

        let mut bonds = Vec::with_capacity(self.bonds.len());
        let mut bond_ids = IndexMap::new();
        for (i, b) in self.bonds.into_iter().enumerate() {
            let a = resolve(&b.a)?;
            let bb = resolve(&b.b)?;
            if let Some(id) = check_id(b.id)? {
                bond_ids.insert(i, id);
            }
            bonds.push(LocalizedBond {
                a,
                b: bb,
                bond: b.bond,
            });
        }

        let mut dative_bonds = Vec::with_capacity(self.dative_bonds.len());
        let mut dative_bond_ids = IndexMap::new();
        for (i, db) in self.dative_bonds.into_iter().enumerate() {
            let donor = resolve(&db.donor)?;
            let acceptor = resolve(&db.acceptor)?;
            if let Some(id) = check_id(db.id)? {
                dative_bond_ids.insert(i, id);
            }
            dative_bonds.push(DativeBond {
                donor,
                acceptor,
                bond: db.bond,
            });
        }

        let mut aromatic_systems = Vec::with_capacity(self.aromatic_systems.len());
        let mut aromatic_system_ids = IndexMap::new();
        for (i, sys) in self.aromatic_systems.into_iter().enumerate() {
            let atom_indices: Vec<usize> =
                sys.atoms.iter().map(&resolve).collect::<Result<_, _>>()?;
            if let Some(id) = check_id(sys.id)? {
                aromatic_system_ids.insert(i, id);
            }
            aromatic_systems.push(AromaticSystem {
                atoms: atom_indices,
            });
        }

        let mut multicenter_bonds = Vec::with_capacity(self.multicenter_bonds.len());
        let mut multicenter_bond_ids = IndexMap::new();
        for (i, mc) in self.multicenter_bonds.into_iter().enumerate() {
            let atom_indices: Vec<usize> =
                mc.atoms.iter().map(&resolve).collect::<Result<_, _>>()?;
            if let Some(id) = check_id(mc.id)? {
                multicenter_bond_ids.insert(i, id);
            }
            multicenter_bonds.push(MulticenterBond {
                atoms: atom_indices,
            });
        }

        let mut noncovalent_bonds = Vec::with_capacity(self.noncovalent_bonds.len());
        let mut noncovalent_bond_ids = IndexMap::new();
        for (i, nc) in self.noncovalent_bonds.into_iter().enumerate() {
            let a = resolve(&nc.a)?;
            let bb = resolve(&nc.b)?;
            if let Some(id) = check_id(nc.id)? {
                noncovalent_bond_ids.insert(i, id);
            }
            noncovalent_bonds.push(NoncovalentBond {
                a,
                b: bb,
                bond: nc.bond,
            });
        }

        let ast = MoleculeAst {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            charge: self.charge,
            spin: self.spin,
            constraints: Vec::new(),
        };
        let metadata = Metadata {
            atom_tags,
            atom_aliases,
            bond_ids,
            dative_bond_ids,
            aromatic_system_ids,
            multicenter_bond_ids,
            noncovalent_bond_ids,
        };
        Ok(MoleculeAstWrapper { ast, metadata })
    }

    fn build_alias_table(raw: &[String]) -> Result<IndexMap<String, String>, ParseError> {
        if !raw.len().is_multiple_of(2) {
            return Err(ParseError::WrongFieldType {
                field: "atom-aliases".to_string(),
                expected: "flat vector of keyword/atom-spec pairs (even length)".to_string(),
            });
        }
        let mut table = IndexMap::new();
        for pair in raw.chunks(2) {
            let name = &pair[0];
            let def = &pair[1];
            if table.contains_key(name) {
                return Err(ParseError::DuplicateId(name.clone()));
            }
            table.insert(name.clone(), def.clone());
        }
        Ok(table)
    }

    fn resolve_entry(
        entry: AtomEntryInput,
        alias_table: &IndexMap<String, String>,
    ) -> Result<(Option<String>, String), ParseError> {
        match entry {
            AtomEntryInput::Str(s) => {
                if let Some(def) = alias_table.get(&s) {
                    Ok((None, def.clone()))
                } else {
                    Ok((None, s))
                }
            }
            AtomEntryInput::Tagged(tag, inner) => {
                let (_, atom_str) = Self::resolve_entry(*inner, alias_table)?;
                Ok((Some(tag), atom_str))
            }
        }
    }
}

fn read_molecule_input(de: &mut EdnStreamDeserializer<'_>) -> Result<RawMoleculeAst, EdnError> {
    de.consume_byte(b'{')?;

    let mut atoms: Option<Vec<AtomEntryInput>> = None;
    let mut bonds: Option<Vec<LocalizedBondInput>> = None;
    let mut dative_bonds: Vec<DativeBondInput> = Vec::new();
    let mut aromatic_systems: Vec<AromaticSystemInput> = Vec::new();
    let mut multicenter_bonds: Vec<MulticenterBondInput> = Vec::new();
    let mut noncovalent_bonds: Vec<NoncovalentBondInput> = Vec::new();
    let mut charge: Option<i64> = None;
    let mut spin: Option<SpinState> = None;
    let mut atom_aliases: Vec<String> = Vec::new();

    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "atoms" => atoms = Some(read_seq(de, read_atom_entry)?),
            "bonds" => bonds = Some(read_seq(de, read_localized_bond)?),
            "dative" => dative_bonds = read_seq(de, read_dative_bond)?,
            "aromatic" => aromatic_systems = read_seq(de, read_aromatic_system)?,
            "multicenter" => multicenter_bonds = read_seq(de, read_multicenter_bond)?,
            "noncovalent" => noncovalent_bonds = read_seq(de, read_noncovalent_bond)?,
            "charge" => charge = Some(de.read_i64()?),
            "spin" => {
                let s = de.read_string_or_keyword()?;
                spin = Some(
                    SpinState::from_str(s.as_ref()).map_err(|e| DeError::subgrammar("spin", e))?,
                );
            }
            "atom-aliases" | "aliases" => {
                atom_aliases = read_seq(de, |d| {
                    d.read_string_or_keyword()
                        .map(|c| c.into_owned())
                        .map_err(Into::into)
                })?;
            }
            _ => de.read_skip_value()?,
        }
    }

    let atoms = atoms.ok_or_else(|| DeError::MissingField {
        key: "atoms".to_string(),
        path: Vec::new(),
    })?;
    let bonds = bonds.ok_or_else(|| DeError::MissingField {
        key: "bonds".to_string(),
        path: Vec::new(),
    })?;

    Ok(RawMoleculeAst {
        atoms,
        bonds,
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        atom_aliases,
        charge,
        spin,
    })
}

fn read_seq<T, F>(de: &mut EdnStreamDeserializer<'_>, mut element: F) -> Result<Vec<T>, EdnError>
where
    F: FnMut(&mut EdnStreamDeserializer<'_>) -> Result<T, EdnError>,
{
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            return Ok(out);
        }
        out.push(element(de)?);
    }
}

fn read_atom_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomEntryInput, EdnError> {
    match de.peek_byte()? {
        Some(b'"') | Some(b':') => {
            let s = de.read_string_or_keyword()?;
            Ok(AtomEntryInput::Str(s.into_owned()))
        }
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let tag = de.read_string_or_keyword()?.into_owned();
            let inner = read_atom_entry(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryInput::Tagged(tag, Box::new(inner)))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_atom_ref(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomRefInput, EdnError> {
    match de.peek_byte()? {
        Some(b':') => Ok(AtomRefInput::Tag(de.read_keyword_name()?.into_owned())),
        Some(b'"') => Ok(AtomRefInput::Tag(de.read_string()?.into_owned())),
        Some(b) if b.is_ascii_digit() || b == b'-' || b == b'+' => {
            let n = de.read_i64()?;
            let idx = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                value: n.to_string(),
                target: "atom index",
                path: Vec::new(),
            })?;
            Ok(AtomRefInput::Index(idx))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_bond_spec(de: &mut EdnStreamDeserializer<'_>) -> Result<BondAst, EdnError> {
    let s = de.read_string_or_keyword()?;
    let aliases = super::bond::builtin_bond_aliases();
    if let Some(ast) = aliases.get_by_left(s.as_ref()) {
        return Ok(ast.clone());
    }
    BondAst::from_str(s.as_ref()).map_err(|e| DeError::subgrammar("bond", e).into())
}

fn read_localized_bond(de: &mut EdnStreamDeserializer<'_>) -> Result<LocalizedBondInput, EdnError> {
    match de.peek_byte()? {
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_spec(de)?;
            de.consume_byte(b']')?;
            Ok(LocalizedBondInput {
                id: None,
                a,
                b,
                bond,
            })
        }
        Some(b'{') => {
            let (id, a, b, bond) = read_endpoint_bond_map(de, EndpointBondKind::Localized)?;
            Ok(LocalizedBondInput { id, a, b, bond })
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_dative_bond(de: &mut EdnStreamDeserializer<'_>) -> Result<DativeBondInput, EdnError> {
    let (id, donor, acceptor, bond) = read_endpoint_bond_map(de, EndpointBondKind::Dative)?;
    Ok(DativeBondInput {
        id,
        donor,
        acceptor,
        bond,
    })
}

fn read_noncovalent_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondInput, EdnError> {
    let (id, a, b, bond) = read_endpoint_bond_map(de, EndpointBondKind::Noncovalent)?;
    Ok(NoncovalentBondInput { id, a, b, bond })
}

fn read_aromatic_system(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystemInput, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(AromaticSystemInput { id, atoms })
}

fn read_multicenter_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBondInput, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(MulticenterBondInput { id, atoms })
}

#[derive(Copy, Clone)]
enum EndpointBondKind {
    Localized,
    Dative,
    Noncovalent,
}

impl EndpointBondKind {
    fn first_key(self) -> &'static str {
        match self {
            EndpointBondKind::Dative => "donor",
            _ => "a",
        }
    }
    fn second_key(self) -> &'static str {
        match self {
            EndpointBondKind::Dative => "acceptor",
            _ => "b",
        }
    }
}

fn read_endpoint_bond_map(
    de: &mut EdnStreamDeserializer<'_>,
    kind: EndpointBondKind,
) -> Result<(Option<String>, AtomRefInput, AtomRefInput, BondAst), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<String> = None;
    let mut a: Option<AtomRefInput> = None;
    let mut b: Option<AtomRefInput> = None;
    let mut bond: Option<BondAst> = None;
    let first_key = kind.first_key();
    let second_key = kind.second_key();

    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        let key_ref = key.as_ref();
        if key_ref == "id" {
            id = Some(de.read_keyword_name()?.into_owned());
        } else if key_ref == first_key {
            a = Some(read_atom_ref(de)?);
        } else if key_ref == second_key {
            b = Some(read_atom_ref(de)?);
        } else if key_ref == "bond" {
            bond = Some(read_bond_spec(de)?);
        } else {
            de.read_skip_value()?;
        }
    }

    let a = a.ok_or_else(|| DeError::MissingField {
        key: first_key.to_string(),
        path: Vec::new(),
    })?;
    let b = b.ok_or_else(|| DeError::MissingField {
        key: second_key.to_string(),
        path: Vec::new(),
    })?;
    let bond = bond.ok_or_else(|| DeError::MissingField {
        key: "bond".to_string(),
        path: Vec::new(),
    })?;
    Ok((id, a, b, bond))
}

fn read_atoms_bond_map(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<(Option<String>, Vec<AtomRefInput>), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<String> = None;
    let mut atoms: Option<Vec<AtomRefInput>> = None;
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "id" => id = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_seq(de, read_atom_ref)?),
            _ => de.read_skip_value()?,
        }
    }
    let atoms = atoms.ok_or_else(|| DeError::MissingField {
        key: "atoms".to_string(),
        path: Vec::new(),
    })?;
    Ok((id, atoms))
}

fn unexpected(offset: usize, b: Option<u8>) -> EdnError {
    match b {
        Some(b) => umol_edn::ParseError::UnexpectedToken {
            offset,
            found: b as char,
        }
        .into(),
        None => umol_edn::ParseError::UnexpectedEof { offset }.into(),
    }
}

impl<'de> FromEdn<'de> for MoleculeAstWrapper {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Self::from_edn_str(&edn.to_string()).map_err(|e| DeError::subgrammar("molecule", e))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        parse_molecule_dsl(input).map_err(|e| DeError::subgrammar("molecule", e).into())
    }
}

impl ToEdn for MoleculeAstWrapper {
    fn to_edn(&self) -> Edn<'static> {
        let mut atom_elems = Vec::with_capacity(self.ast.atoms.len());
        for (i, atom) in self.ast.atoms.iter().enumerate() {
            let alias_name = self.metadata.atom_aliases.get_by_right(atom);
            let tag = self.metadata.atom_tags.get(&i);
            let atom_edn = if let Some(alias) = alias_name {
                Edn::Keyword(EdnKeyword::owned(alias.clone()))
            } else {
                atom.to_edn()
            };
            let entry = if let Some(tag_name) = tag {
                Edn::Vector(
                    vec![Edn::Keyword(EdnKeyword::owned(tag_name.clone())), atom_edn].into(),
                )
            } else {
                atom_edn
            };
            atom_elems.push(entry);
        }

        let render_endpoint = |idx: usize| -> Edn<'static> {
            if let Some(tag) = self.metadata.atom_tags.get(&idx) {
                Edn::Keyword(EdnKeyword::owned(tag.clone()))
            } else {
                Edn::Int(idx as i64)
            }
        };

        let bonds_edn: Vec<Edn<'static>> = self
            .ast
            .bonds
            .iter()
            .enumerate()
            .map(|(i, b)| render_localized(b, i, &self.metadata.bond_ids, &render_endpoint))
            .collect();

        let dative_edn: Vec<Edn<'static>> = self
            .ast
            .dative_bonds
            .iter()
            .enumerate()
            .map(|(i, b)| render_dative(b, i, &self.metadata.dative_bond_ids, &render_endpoint))
            .collect();

        let aromatic_edn: Vec<Edn<'static>> = self
            .ast
            .aromatic_systems
            .iter()
            .enumerate()
            .map(|(i, sys)| {
                render_atoms_map(
                    &sys.atoms,
                    i,
                    &self.metadata.aromatic_system_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let multicenter_edn: Vec<Edn<'static>> = self
            .ast
            .multicenter_bonds
            .iter()
            .enumerate()
            .map(|(i, mc)| {
                render_atoms_map(
                    &mc.atoms,
                    i,
                    &self.metadata.multicenter_bond_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let noncovalent_edn: Vec<Edn<'static>> = self
            .ast
            .noncovalent_bonds
            .iter()
            .enumerate()
            .map(|(i, nc)| {
                render_noncovalent(nc, i, &self.metadata.noncovalent_bond_ids, &render_endpoint)
            })
            .collect();

        let has_aliases = !self.metadata.atom_aliases.is_empty();
        let mut m = EdnMap::with_capacity(10);
        m.insert(Edn::keyword("atoms"), Edn::Vector(atom_elems.into()));
        m.insert(Edn::keyword("bonds"), Edn::Vector(bonds_edn.into()));
        if !self.ast.dative_bonds.is_empty() {
            m.insert(Edn::keyword("dative"), Edn::Vector(dative_edn.into()));
        }
        if !self.ast.aromatic_systems.is_empty() {
            m.insert(Edn::keyword("aromatic"), Edn::Vector(aromatic_edn.into()));
        }
        if !self.ast.multicenter_bonds.is_empty() {
            m.insert(
                Edn::keyword("multicenter"),
                Edn::Vector(multicenter_edn.into()),
            );
        }
        if !self.ast.noncovalent_bonds.is_empty() {
            m.insert(
                Edn::keyword("noncovalent"),
                Edn::Vector(noncovalent_edn.into()),
            );
        }
        if let Some(charge) = self.ast.charge {
            m.insert(Edn::keyword("charge"), Edn::Int(charge));
        }
        if let Some(spin) = &self.ast.spin {
            m.insert(
                Edn::keyword("spin"),
                Edn::Str(std::borrow::Cow::Owned(spin.to_string())),
            );
        }
        if has_aliases {
            let mut alias_elems = Vec::with_capacity(self.metadata.atom_aliases.len() * 2);
            for (name, atom) in self.metadata.atom_aliases.iter() {
                alias_elems.push(Edn::Keyword(EdnKeyword::owned(name.clone())));
                alias_elems.push(atom.to_edn());
            }
            m.insert(
                Edn::keyword("atom-aliases"),
                Edn::Vector(alias_elems.into()),
            );
        }
        Edn::Map(m)
    }
}

fn render_localized(
    b: &LocalizedBond,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let a = render_endpoint(b.a);
    let bb = render_endpoint(b.b);
    let bond = b.bond.to_edn();
    if let Some(id) = ids.get(&i) {
        let mut m = EdnMap::with_capacity(4);
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.clone())),
        );
        m.insert(Edn::keyword("a"), a);
        m.insert(Edn::keyword("b"), bb);
        m.insert(Edn::keyword("bond"), bond);
        Edn::Map(m)
    } else {
        Edn::Vector(vec![a, bb, bond].into())
    }
}

fn render_dative(
    b: &DativeBond,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let donor = render_endpoint(b.donor);
    let acceptor = render_endpoint(b.acceptor);
    let bond = b.bond.to_edn();
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = ids.get(&i) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.clone())),
        );
    }
    m.insert(Edn::keyword("donor"), donor);
    m.insert(Edn::keyword("acceptor"), acceptor);
    m.insert(Edn::keyword("bond"), bond);
    Edn::Map(m)
}

fn render_noncovalent(
    b: &NoncovalentBond,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let a = render_endpoint(b.a);
    let bb = render_endpoint(b.b);
    let bond = b.bond.to_edn();
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = ids.get(&i) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.clone())),
        );
    }
    m.insert(Edn::keyword("a"), a);
    m.insert(Edn::keyword("b"), bb);
    m.insert(Edn::keyword("bond"), bond);
    Edn::Map(m)
}

fn render_atoms_map(
    atoms: &[usize],
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let atom_vec: Vec<Edn<'static>> = atoms.iter().map(|&a| render_endpoint(a)).collect();
    let mut m = EdnMap::with_capacity(2);
    if let Some(id) = ids.get(&i) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(id.clone())),
        );
    }
    m.insert(Edn::keyword("atoms"), Edn::Vector(atom_vec.into()));
    Edn::Map(m)
}

impl fmt::Display for MoleculeAstWrapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::{e, Element};

    use super::*;
    use crate::ast::atom::{ElementAst, HydrogenAst};
    use crate::ast::config::MoleculeAstConfig;
    use crate::ast::value::ValueAst;
    use crate::ast::{FromAst, ToAst};
    use crate::graph_ir::molecule_builder::MoleculeBuilder;

    fn atoms(a: Vec<AtomAst>) -> MoleculeAst {
        MoleculeAst {
            atoms: a,
            ..Default::default()
        }
    }

    #[rstest]
    #[case::empty(
        r#"{:atoms [] :bonds []}"#,
        MoleculeAst::default(),
        Metadata::default()
    )]
    #[case::atom(
        r#"{:atoms ["C"] :bonds []}"#,
        atoms(vec![AtomAst::from_element(e!(C))]),
        Metadata::default()
    )]
    #[case::atom_tagged(
        r#"{:atoms [[:C "C"]] :bonds []}"#,
        atoms(vec![AtomAst::from_element(e!(C))]),
        Metadata {
            atom_tags: IndexMap::from([(0, "C".to_string())]),
            ..Default::default()
        }
    )]
    #[case::bond(
        r#"{:atoms ["N" "N"] :bonds [[0 1 :triple]]}"#,
        MoleculeAst {
            atoms: vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))],
            bonds: vec![LocalizedBond { a: 0, b: 1, bond: BondAst::from_order(3) }],
            ..Default::default()
        },
        Metadata::default()
    )]
    #[case::bond_with_tags(
        r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#,
        MoleculeAst {
            atoms: vec![AtomAst::from_element(e!(C)), AtomAst::from_element(e!(O))],
            bonds: vec![LocalizedBond { a: 0, b: 1, bond: BondAst::from_order(1) }],
            ..Default::default()
        },
        Metadata {
            atom_tags: IndexMap::from([(0, "C".to_string()), (1, "O".to_string())]),
            ..Default::default()
        }
    )]
    #[case::bond_id(
        r#"{:atoms ["H" "F"] :bonds [{:id :b1 :a 0 :b 1 :bond :single}]}"#,
        MoleculeAst {
            atoms: vec![AtomAst::from_element(e!(H)), AtomAst::from_element(e!(F))],
            bonds: vec![LocalizedBond { a: 0, b: 1, bond: BondAst::from_order(1) }],
            ..Default::default()
        },
        Metadata {
            bond_ids: IndexMap::from([(0, "b1".to_string())]),
            ..Default::default()
        }
    )]
    #[case::charge(
        r#"{:atoms [[:F "F#c-"]] :bonds [] :charge -1}"#,
        MoleculeAst {
            atoms: vec![AtomAst {
                element: ElementAst::Lit(Element::F),
                isotope_mass: None,
                implicit_hydrogens: None,
                charge: Some(ValueAst::Lit(-1)),
                lone_pairs: None,
                unpaired_electrons: None,
                multiplicity: None,
                valence: None,
                donated_pairs: None,
                accepted_pairs: None,
                aromatic_valence: None,
                multicenter_valence: None,
            }],
            charge: Some(-1),
            ..Default::default()
        },
        Metadata {
            atom_tags: IndexMap::from([(0, "F".to_string())]),
            ..Default::default()
        }
    )]
    #[case::alias_indexed(
        r#"{:atoms [:ch] :bonds [] :aliases [:ch "C #h1"]}"#,
        atoms(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: None,
            implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(1))),
            charge: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }]),
        Metadata {
            atom_aliases: BiMap::from_iter([(
                "ch".to_string(),
                AtomAst {
                    element: ElementAst::Lit(Element::C),
                    isotope_mass: None,
                    implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(1))),
                    charge: None,
                    lone_pairs: None,
                    unpaired_electrons: None,
                    multiplicity: None,
                    valence: None,
                    donated_pairs: None,
                    accepted_pairs: None,
                    aromatic_valence: None,
                    multicenter_valence: None,
                },
            )]),
            ..Default::default()
        }
    )]
    #[case::alias_reused(
        r#"{:atoms [:n :n] :bonds [[0 1 :single]] :aliases [:n "N"]}"#,
        MoleculeAst {
            atoms: vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))],
            bonds: vec![LocalizedBond { a: 0, b: 1, bond: BondAst::from_order(1) }],
            ..Default::default()
        },
        Metadata {
            atom_aliases: BiMap::from_iter([("n".to_string(), AtomAst::from_element(e!(N)))]),
            ..Default::default()
        }
    )]
    fn test_parse_molecule_dsl(
        #[case] input: &str,
        #[case] expected_ast: MoleculeAst,
        #[case] expected_meta: Metadata,
    ) {
        let dsl = parse_molecule_dsl(input).unwrap();
        assert_eq!(dsl.ast, expected_ast);
        assert_eq!(dsl.metadata, expected_meta);
    }

    #[test]
    fn test_parse_molecule_dsl_dative() {
        let dsl = parse_molecule_dsl(
            r#"{:atoms [[:B "B #h3"] [:N "N #h3"]]
                :bonds []
                :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#,
        )
        .unwrap();
        assert_eq!(dsl.ast.dative_bonds.len(), 1);
        assert_eq!(dsl.ast.dative_bonds[0].donor, 1);
        assert_eq!(dsl.ast.dative_bonds[0].acceptor, 0);
        assert_eq!(
            dsl.metadata.dative_bond_ids.get(&0),
            Some(&"d1".to_string())
        );
    }

    #[test]
    fn test_parse_molecule_dsl_aromatic() {
        let dsl = parse_molecule_dsl(
            r#"{:atoms [[:C1 :ch] [:C2 :ch] [:C3 :ch] [:C4 :ch] [:C5 :ch] [:C6 :ch]]
                :bonds [[:C1 :C2 :single] [:C2 :C3 :single] [:C3 :C4 :single] [:C4 :C5 :single] [:C5 :C6 :single] [:C6 :C1 :single]]
                :aromatic [{:id :ar1 :atoms [:C1 :C2 :C3 :C4 :C5 :C6]}]
                :aliases [:ch "C #h1 #v2 #a1"]}"#,
        )
        .unwrap();
        assert_eq!(dsl.ast.aromatic_systems.len(), 1);
        assert_eq!(
            dsl.metadata.aromatic_system_ids.get(&0),
            Some(&"ar1".to_string())
        );
        assert_eq!(dsl.ast.aromatic_systems[0].atoms, vec![0, 1, 2, 3, 4, 5]);
    }

    #[rstest]
    #[case::non_map("3")]
    #[case::missing_atoms(r#"{:bonds []}"#)]
    #[case::missing_bonds(r#"{:atoms ["C"]}"#)]
    #[case::unknown_endpoint(r#"{:atoms [[:C "C"]] :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#)]
    #[case::bad_atom_string(r##"{:atoms ["#h3"] :bonds []}"##)]
    #[case::trailing_content(r#"{:atoms ["C"] :bonds []} :extra :junk"#)]
    #[case::duplicate_tag(r#"{:atoms [[:C "C"] [:C "O"]] :bonds []}"#)]
    #[case::duplicate_alias(r#"{:aliases [:ch "C #h1" :ch "C #h2"] :atoms [] :bonds []}"#)]
    fn test_parse_molecule_dsl_invalid(#[case] input: &str) {
        assert!(
            parse_molecule_dsl(input).is_err(),
            "{input:?} should fail but succeeded"
        );
    }

    #[rstest]
    #[case::plain_vector(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::tagged_vector(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
    #[case::alias_vector(r#"{:atoms [:ch :ch "O"] :bonds [[0 2 :single]] :aliases [:ch "C #h1"]}"#)]
    #[case::mixed(
        r#"{:atoms [:ch [:ring-O :ch] "O"] :bonds [[0 :ring-O :single]] :aliases [:ch "C #h1"]}"#
    )]
    fn test_edn_roundtrip(#[case] input: &str) {
        let dsl1 = parse_molecule_dsl(input).unwrap();
        let edn = dsl1.to_string();
        let dsl2 = parse_molecule_dsl(&edn).unwrap();
        assert_eq!(dsl1.ast, dsl2.ast);
        assert_eq!(dsl1.metadata, dsl2.metadata);
    }

    #[test]
    fn test_builder_metadata_roundtrip() {
        let dsl =
            parse_molecule_dsl(r#"{:atoms [[:C "C #h3"] [:O "O #h1"]] :bonds [[:C :O :single]]}"#)
                .unwrap();

        let cfg = MoleculeAstConfig::zeroed();
        let builder = MoleculeBuilder::from_ast(&dsl.ast, &cfg).unwrap();
        let dsl2 = MoleculeAstWrapper::new(builder.to_ast(&cfg), dsl.metadata.clone());

        assert_eq!(dsl.ast.atoms.len(), dsl2.ast.atoms.len());
        assert!(dsl2.metadata.atom_tags.values().any(|t| t == "C"));
        assert!(dsl2.metadata.atom_tags.values().any(|t| t == "O"));
        let edn = dsl2.to_string();
        assert!(edn.contains(":C"), "serialized form should contain :C tag");
        assert!(edn.contains(":O"), "serialized form should contain :O tag");
    }

    #[test]
    fn test_builder_no_metadata_produces_indexed() {
        let dsl =
            parse_molecule_dsl(r#"{:atoms [[:C "C #h3"] [:O "O #h1"]] :bonds [[:C :O :single]]}"#)
                .unwrap();

        let cfg = MoleculeAstConfig::zeroed();
        let builder = MoleculeBuilder::from_ast(&dsl.ast, &cfg).unwrap();
        let ast2 = builder.to_ast(&cfg);

        assert_eq!(ast2.bonds[0].a, 0);
        assert_eq!(ast2.bonds[0].b, 1);
    }

    #[test]
    fn test_alias_tag_disjointness_error() {
        let result = parse_molecule_dsl(r#"{:atoms [[:ch "C"]] :bonds [] :aliases [:ch "C #h1"]}"#);
        assert!(result.is_err(), "alias and tag with same name should fail");
    }

    #[rstest]
    #[case::plain(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::tagged(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
    #[case::aliased(r#"{:atoms [:ch :ch "O"] :bonds [[0 2 :single]] :aliases [:ch "C #h1"]}"#)]
    #[case::charge_spin(r##"{:atoms ["C"] :bonds [] :charge -1 :spin "#u1"}"##)]
    #[case::dative(
        r#"{:atoms [[:N "N"] [:B "B"]] :bonds [[:N :B :single]] :dative [{:donor :N :acceptor :B :bond :single}]}"#
    )]
    #[case::aromatic(
        r#"{:atoms [[:C1 "C"] [:C2 "C"]] :bonds [[:C1 :C2 :single]] :aromatic [{:id :ar1 :atoms [:C1 :C2]}]}"#
    )]
    fn test_molecule_dsl_to_edn_roundtrip(#[case] input: &str) {
        let dsl = parse_molecule_dsl(input).unwrap();
        let edn = dsl.to_edn();
        let back = MoleculeAstWrapper::from_edn(&edn).unwrap();
        assert_eq!(dsl, back);
    }
}
