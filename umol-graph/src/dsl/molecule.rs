//! Molecule DSL definitions

use std::borrow::Cow;
use std::collections::HashSet;
use std::str::FromStr;
use std::{fmt, mem};

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_shared::spin::SpinState;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use super::atom::parse_atom_dsl;
use super::error::ParseError;
use crate::api::pattern::AtomPattern;
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::{
    AromaticValenceConstraint, AtomConstraint, BondConstraint, MoleculeConstraint,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::{AromaticSystemIdx, AtomIdx, BondIdx, MulticenterBondIdx};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_ids: IndexMap<usize, String>,
    pub atom_aliases: BiMap<String, AtomPattern>,
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
        // Inline here -- no need for free fn + trait impl
        parse_molecule_dsl(s)
    }
}

#[derive(Clone, Debug)]
enum AtomEntryInput {
    Str(String),
    WithId(String, Box<AtomEntryInput>),
}

#[derive(Clone, Debug)]
enum AtomRefInput {
    Index(usize),
    Id(String),
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
    constraints: Vec<Edn<'static>>,
}

impl RawMoleculeAst {
    fn into_dsl(self) -> Result<MoleculeAstWrapper, ParseError> {
        let alias_table = Self::build_alias_table(&self.atom_aliases)?;

        let mut atoms: Vec<AtomAst> = Vec::with_capacity(self.atoms.len());
        let mut lifted_constraints: Vec<MoleculeConstraint> = Vec::new();
        let mut atom_ids: IndexMap<usize, String> = IndexMap::new();
        let mut id_to_index: IndexMap<String, usize> = IndexMap::new();
        let mut atom_aliases: BiMap<String, AtomPattern> = BiMap::new();

        for entry in self.atoms {
            let pos = atoms.len();
            let (id, atom_str) = Self::resolve_entry(entry, &alias_table)?;
            let atom_pattern = parse_atom_dsl(&atom_str)?;
            if let Some(id_name) = id {
                if id_to_index.contains_key(&id_name) || alias_table.contains_key(&id_name) {
                    return Err(ParseError::DuplicateId(id_name));
                }
                id_to_index.insert(id_name.clone(), pos);
                atom_ids.insert(pos, id_name);
            }
            for c in atom_pattern.constraints {
                lifted_constraints.push(MoleculeConstraint::AtomPred(AtomIdx(pos as u32), c));
            }
            atoms.push(atom_pattern.ast);
        }

        for (name, def) in &alias_table {
            let atom_pattern = parse_atom_dsl(def)?;
            atom_aliases.insert(name.clone(), atom_pattern);
        }

        let atom_count = atoms.len();
        let resolve = |r: &AtomRefInput| -> Result<AtomIdx, ParseError> {
            match r {
                AtomRefInput::Index(i) => {
                    if *i < atom_count {
                        Ok(AtomIdx(*i as u32))
                    } else {
                        Err(ParseError::InvalidAtomIndex(i.to_string()))
                    }
                }
                AtomRefInput::Id(name) => id_to_index
                    .get(name)
                    .copied()
                    .map(|i| AtomIdx(i as u32))
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

        let mut bond_list = Vec::new();
        let mut bond_ids = IndexMap::new();
        for (i, b) in self.bonds.into_iter().enumerate() {
            let a = resolve(&b.a)?;
            let bb = resolve(&b.b)?;
            if let Some(id) = check_id(b.id)? {
                bond_ids.insert(i, id);
            }
            bond_list.push((a, bb, b.bond));
        }

        let mut dative_list = Vec::new();
        let mut dative_bond_ids = IndexMap::new();
        for (i, db) in self.dative_bonds.into_iter().enumerate() {
            let donor = resolve(&db.donor)?;
            let acceptor = resolve(&db.acceptor)?;
            if let Some(id) = check_id(db.id)? {
                dative_bond_ids.insert(i, id);
            }
            dative_list.push((donor, acceptor, db.bond));
        }

        let mut aromatic_list = Vec::new();
        let mut aromatic_system_ids = IndexMap::new();
        for (i, sys) in self.aromatic_systems.into_iter().enumerate() {
            let atom_indices: Vec<AtomIdx> =
                sys.atoms.iter().map(&resolve).collect::<Result<_, _>>()?;
            if let Some(id) = check_id(sys.id)? {
                aromatic_system_ids.insert(i, id);
            }
            aromatic_list.push((atom_indices, AromaticSystemAst::default()));
        }

        let mut multicenter_list = Vec::new();
        let mut multicenter_bond_ids = IndexMap::new();
        for (i, mc) in self.multicenter_bonds.into_iter().enumerate() {
            let atom_indices: Vec<AtomIdx> =
                mc.atoms.iter().map(&resolve).collect::<Result<_, _>>()?;
            if let Some(id) = check_id(mc.id)? {
                multicenter_bond_ids.insert(i, id);
            }
            multicenter_list.push((atom_indices, MulticenterBondAst {}));
        }

        let mut noncovalent_list = Vec::new();
        let mut noncovalent_bond_ids = IndexMap::new();
        for (i, nc) in self.noncovalent_bonds.into_iter().enumerate() {
            let a = resolve(&nc.a)?;
            let bb = resolve(&nc.b)?;
            if let Some(id) = check_id(nc.id)? {
                noncovalent_bond_ids.insert(i, id);
            }
            noncovalent_list.push((a, bb, nc.bond));
        }

        let mut constraints = lifted_constraints;
        if let Some(charge) = self.charge {
            constraints.push(MoleculeConstraint::TotalCharge(ValueAst::Lit(charge)));
        }
        if let Some(spin) = self.spin {
            constraints.push(MoleculeConstraint::TotalSpin(SpinStateAst::Lit(spin)));
        }
        if !self.constraints.is_empty() {
            let resolver = ConstraintResolver {
                atom_count,
                atom_ids: &atom_ids,
                bond_count: bond_list.len(),
                bond_ids: &bond_ids,
                aromatic_count: aromatic_list.len(),
                aromatic_ids: &aromatic_system_ids,
                multicenter_count: multicenter_list.len(),
                multicenter_ids: &multicenter_bond_ids,
            };
            for entry in &self.constraints {
                constraints.push(parse_molecule_constraint(entry, &resolver)?);
            }
        }

        let ast = MoleculeAst::new(
            atoms,
            bond_list,
            dative_list,
            noncovalent_list,
            aromatic_list,
            multicenter_list,
            constraints,
        );
        let metadata = Metadata {
            atom_ids,
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
            AtomEntryInput::WithId(id, inner) => {
                let (_, atom_str) = Self::resolve_entry(*inner, alias_table)?;
                Ok((Some(id), atom_str))
            }
        }
    }
}

struct ConstraintResolver<'a> {
    atom_count: usize,
    atom_ids: &'a IndexMap<usize, String>,
    bond_count: usize,
    bond_ids: &'a IndexMap<usize, String>,
    aromatic_count: usize,
    aromatic_ids: &'a IndexMap<usize, String>,
    multicenter_count: usize,
    multicenter_ids: &'a IndexMap<usize, String>,
}

fn reverse_lookup_id(ids: &IndexMap<usize, String>, name: &str) -> Option<usize> {
    ids.iter()
        .find(|(_, n)| n.as_str() == name)
        .map(|(i, _)| *i)
}

fn resolve_ref<T: From<usize>>(
    edn: &Edn<'_>,
    count: usize,
    ids: &IndexMap<usize, String>,
    kind: &'static str,
) -> Result<T, ParseError> {
    let i = match edn {
        Edn::Int(n) => {
            let i = usize::try_from(*n)
                .map_err(|_| ParseError::InvalidValue(format!("{kind} index {n}")))?;
            if i >= count {
                return Err(ParseError::InvalidValue(format!(
                    "{kind} index {i} out of range"
                )));
            }
            i
        }
        Edn::Keyword(k) => reverse_lookup_id(ids, k.name()).ok_or_else(|| {
            ParseError::InvalidValue(format!("unknown {kind} id :{}", k.name()))
        })?,
        other => {
            return Err(ParseError::InvalidValue(format!(
                "{kind} ref: {other}"
            )));
        }
    };
    Ok(T::from(i))
}

fn parse_value_ast(edn: &Edn<'_>) -> Result<ValueAst, ParseError> {
    match edn {
        Edn::Int(n) => Ok(ValueAst::Lit(*n)),
        Edn::Nil => Ok(ValueAst::Undetermined),
        Edn::Keyword(k) if k.name() == "undetermined" => Ok(ValueAst::Undetermined),
        Edn::Vector(v) => {
            let mut out = Vec::with_capacity(v.len());
            for e in v.iter() {
                let Edn::Int(n) = e else {
                    return Err(ParseError::InvalidValue(format!("value-set entry: {e}")));
                };
                out.push(*n);
            }
            Ok(ValueAst::LitSet(out))
        }
        other => Err(ParseError::InvalidValue(format!("value: {other}"))),
    }
}

fn parse_molecule_constraint(
    entry: &Edn<'_>,
    r: &ConstraintResolver<'_>,
) -> Result<MoleculeConstraint, ParseError> {
    let map = match entry {
        Edn::Map(m) => m,
        other => {
            return Err(ParseError::InvalidValue(format!(
                "constraint entry must be a map, got: {other}"
            )));
        }
    };
    if map.len() != 1 {
        return Err(ParseError::InvalidValue(
            "constraint entry must have exactly one key".to_string(),
        ));
    }
    let (key, value) = map.iter().next().unwrap();
    let key_name = match key {
        Edn::Keyword(k) => k.name(),
        other => {
            return Err(ParseError::InvalidValue(format!(
                "constraint key must be a keyword, got: {other}"
            )));
        }
    };
    match key_name {
        "atom-pred" => {
            let (atom_ref, form) = expect_pair(value, "atom-pred")?;
            let idx: AtomIdx = resolve_ref(atom_ref, r.atom_count, r.atom_ids, "atom")?;
            let c = parse_atom_constraint_form(form)?;
            Ok(MoleculeConstraint::AtomPred(idx, c))
        }
        "bond-pred" => {
            let (bond_ref, form) = expect_pair(value, "bond-pred")?;
            let idx: BondIdx = resolve_ref(bond_ref, r.bond_count, r.bond_ids, "bond")?;
            let c = parse_bond_constraint_form(form)?;
            Ok(MoleculeConstraint::BondPred(idx, c))
        }
        "total-charge" => Ok(MoleculeConstraint::TotalCharge(parse_value_ast(value)?)),
        "total-spin" => {
            let s = match value {
                Edn::Str(s) => SpinState::from_str(s.as_ref())?,
                Edn::Keyword(k) => SpinState::from_str(k.name())?,
                other => {
                    return Err(ParseError::InvalidValue(format!("total-spin: {other}")));
                }
            };
            Ok(MoleculeConstraint::TotalSpin(SpinStateAst::Lit(s)))
        }
        "aromatic-electron-count" => {
            let (sys_ref, val) = expect_pair(value, "aromatic-electron-count")?;
            let idx: AromaticSystemIdx =
                resolve_ref(sys_ref, r.aromatic_count, r.aromatic_ids, "aromatic system")?;
            Ok(MoleculeConstraint::AromaticElectronCount(
                idx,
                parse_value_ast(val)?,
            ))
        }
        "multicenter-electron-count" => {
            let (mc_ref, val) = expect_pair(value, "multicenter-electron-count")?;
            let idx: MulticenterBondIdx = resolve_ref(
                mc_ref,
                r.multicenter_count,
                r.multicenter_ids,
                "multicenter bond",
            )?;
            Ok(MoleculeConstraint::MulticenterElectronCount(
                idx,
                parse_value_ast(val)?,
            ))
        }
        "bond-order-sum" => {
            let m = expect_map(value, "bond-order-sum")?;
            let bonds_edn = m
                .get_keyword("bonds")
                .ok_or_else(|| ParseError::MissingKey(":bonds".to_string()))?;
            let equals_edn = m
                .get_keyword("equals")
                .ok_or_else(|| ParseError::MissingKey(":equals".to_string()))?;
            let bonds = match bonds_edn {
                Edn::Vector(v) => {
                    let mut out = Vec::with_capacity(v.len());
                    for e in v.iter() {
                        out.push(resolve_ref::<BondIdx>(e, r.bond_count, r.bond_ids, "bond")?);
                    }
                    out
                }
                other => {
                    return Err(ParseError::InvalidValue(format!(
                        "bond-order-sum :bonds must be a vector, got {other}"
                    )));
                }
            };
            Ok(MoleculeConstraint::BondOrderSum(
                bonds,
                parse_value_ast(equals_edn)?,
            ))
        }
        "connected" => {
            let v = match value {
                Edn::Vector(v) => v,
                other => {
                    return Err(ParseError::InvalidValue(format!(
                        "connected: expected vector, got {other}"
                    )));
                }
            };
            let mut atoms = Vec::with_capacity(v.len());
            for e in v.iter() {
                atoms.push(resolve_ref::<AtomIdx>(e, r.atom_count, r.atom_ids, "atom")?);
            }
            Ok(MoleculeConstraint::Connected(atoms))
        }
        "sub-pattern" => {
            let m = expect_map(value, "sub-pattern")?;
            let anchor_edn = m
                .get_keyword("anchor")
                .ok_or_else(|| ParseError::MissingKey(":anchor".to_string()))?;
            let pattern_edn = m
                .get_keyword("pattern")
                .ok_or_else(|| ParseError::MissingKey(":pattern".to_string()))?;
            let anchor: AtomIdx = resolve_ref(anchor_edn, r.atom_count, r.atom_ids, "atom")?;
            let wrapper = MoleculeAstWrapper::from_edn(pattern_edn)
                .map_err(|e| ParseError::InvalidValue(format!("sub-pattern: {e}")))?;
            Ok(MoleculeConstraint::SubPattern {
                anchor,
                pattern: Box::new(wrapper.ast),
            })
        }
        "and" | "or" => {
            let v = match value {
                Edn::Vector(v) => v,
                other => {
                    return Err(ParseError::InvalidValue(format!(
                        "{key_name}: expected vector, got {other}"
                    )));
                }
            };
            let mut children = Vec::with_capacity(v.len());
            for e in v.iter() {
                children.push(parse_molecule_constraint(e, r)?);
            }
            Ok(if key_name == "and" {
                MoleculeConstraint::And(children)
            } else {
                MoleculeConstraint::Or(children)
            })
        }
        "not" => Ok(MoleculeConstraint::Not(Box::new(parse_molecule_constraint(
            value, r,
        )?))),
        other => Err(ParseError::InvalidValue(format!(
            "unknown constraint key :{other}"
        ))),
    }
}

fn expect_pair<'e>(
    edn: &'e Edn<'e>,
    context: &'static str,
) -> Result<(&'e Edn<'e>, &'e Edn<'e>), ParseError> {
    match edn {
        Edn::Vector(v) if v.len() == 2 => Ok((&v[0], &v[1])),
        other => Err(ParseError::InvalidValue(format!(
            "{context}: expected [ref form] pair, got {other}"
        ))),
    }
}

fn expect_map<'e>(
    edn: &'e Edn<'e>,
    context: &'static str,
) -> Result<&'e EdnMap<'e>, ParseError> {
    match edn {
        Edn::Map(m) => Ok(m),
        other => Err(ParseError::InvalidValue(format!(
            "{context}: expected map, got {other}"
        ))),
    }
}

fn parse_atom_constraint_form(edn: &Edn<'_>) -> Result<AtomConstraint, ParseError> {
    match edn {
        Edn::Keyword(k) if k.name() == "in-ring" => Ok(AtomConstraint::InRing),
        Edn::Map(m) if m.len() == 1 => {
            let (key, value) = m.iter().next().unwrap();
            let key_name = match key {
                Edn::Keyword(k) => k.name(),
                other => {
                    return Err(ParseError::InvalidValue(format!(
                        "atom-constraint key: {other}"
                    )));
                }
            };
            match key_name {
                "valence" => Ok(AtomConstraint::Valence(parse_value_ast(value)?)),
                "aromatic-valence" => Ok(AtomConstraint::AromaticValence(
                    parse_aromatic_valence_form(value)?,
                )),
                "multicenter-valence" => {
                    Ok(AtomConstraint::MulticenterValence(parse_value_ast(value)?))
                }
                "donated-pairs" => Ok(AtomConstraint::DonatedPairs(parse_value_ast(value)?)),
                "accepted-pairs" => Ok(AtomConstraint::AcceptedPairs(parse_value_ast(value)?)),
                "degree" => Ok(AtomConstraint::Degree(parse_value_ast(value)?)),
                "connectivity" => Ok(AtomConstraint::Connectivity(parse_value_ast(value)?)),
                "total-h-count" => Ok(AtomConstraint::TotalHCount(parse_value_ast(value)?)),
                "ring-count" => Ok(AtomConstraint::RingCount(parse_value_ast(value)?)),
                "ring-size" => Ok(AtomConstraint::RingSize(parse_value_ast(value)?)),
                other => Err(ParseError::InvalidValue(format!(
                    "unknown atom-constraint key :{other}"
                ))),
            }
        }
        other => Err(ParseError::InvalidValue(format!(
            "atom-constraint form: {other}"
        ))),
    }
}

fn parse_aromatic_valence_form(
    edn: &Edn<'_>,
) -> Result<AromaticValenceConstraint, ParseError> {
    match edn {
        Edn::Keyword(k) if k.name() == "not-aromatic" => {
            Ok(AromaticValenceConstraint::NotAromatic)
        }
        _ => Ok(AromaticValenceConstraint::Value(parse_value_ast(edn)?)),
    }
}

fn parse_bond_constraint_form(edn: &Edn<'_>) -> Result<BondConstraint, ParseError> {
    match edn {
        Edn::Keyword(k) => match k.name() {
            "ring-bond" => Ok(BondConstraint::RingBond),
            "aromatic" => Ok(BondConstraint::Aromatic),
            other => Err(ParseError::InvalidValue(format!(
                "unknown bond-constraint :{other}"
            ))),
        },
        other => Err(ParseError::InvalidValue(format!(
            "bond-constraint form: {other}"
        ))),
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
    let mut constraints: Vec<Edn<'static>> = Vec::new();

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
            "constraints" => constraints = read_seq(de, read_constraint_entry)?,
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
        constraints,
    })
}

fn read_constraint_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<Edn<'static>, EdnError> {
    let slice = de.read_value_slice()?;
    Edn::from_str(slice).map_err(Into::into)
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
            let id = de.read_string_or_keyword()?.into_owned();
            let inner = read_atom_entry(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryInput::WithId(id, Box::new(inner)))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_atom_ref(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomRefInput, EdnError> {
    match de.peek_byte()? {
        Some(b':') => Ok(AtomRefInput::Id(de.read_keyword_name()?.into_owned())),
        Some(b'"') => Ok(AtomRefInput::Id(de.read_string()?.into_owned())),
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

fn constraint_sort_key(c: &AtomConstraint) -> u8 {
    match c {
        AtomConstraint::Valence(_) => 0,
        AtomConstraint::AromaticValence(_) => 1,
        AtomConstraint::DonatedPairs(_) => 2,
        AtomConstraint::AcceptedPairs(_) => 3,
        AtomConstraint::MulticenterValence(_) => 4,
        AtomConstraint::Degree(_) => 5,
        AtomConstraint::Connectivity(_) => 6,
        AtomConstraint::TotalHCount(_) => 7,
        AtomConstraint::InRing => 8,
        AtomConstraint::RingCount(_) => 9,
        AtomConstraint::RingSize(_) => 10,
    }
}

impl ToEdn for MoleculeAstWrapper {
    fn to_edn(&self) -> Edn<'static> {
        let mut per_atom_derived: Vec<Vec<AtomConstraint>> =
            vec![Vec::new(); self.ast.atoms().count()];
        for (idx, set) in self.ast.constraints().atoms() {
            let i = idx.index();
            if i < per_atom_derived.len() {
                per_atom_derived[i].extend(set.iter().cloned());
            }
        }

        let mut atom_elems = Vec::with_capacity(self.ast.atoms().count());
        for view in self.ast.atoms().iter() {
            let i = view.idx.index();
            let mut derived = mem::take(&mut per_atom_derived[i]);
            derived.sort_by_key(constraint_sort_key);
            let pattern = AtomPattern::with_constraints(view.data.clone(), derived);
            let alias_name = self.metadata.atom_aliases.get_by_right(&pattern);
            let id = self.metadata.atom_ids.get(&i);
            let atom_edn = if let Some(alias) = alias_name {
                Edn::Keyword(EdnKeyword::owned(alias.clone()))
            } else {
                pattern.to_edn()
            };
            let entry = if let Some(id_name) = id {
                Edn::Vector(
                    vec![Edn::Keyword(EdnKeyword::owned(id_name.clone())), atom_edn].into(),
                )
            } else {
                atom_edn
            };
            atom_elems.push(entry);
        }

        let render_endpoint = |idx: usize| -> Edn<'static> {
            if let Some(id) = self.metadata.atom_ids.get(&idx) {
                Edn::Keyword(EdnKeyword::owned(id.clone()))
            } else {
                Edn::Int(idx as i64)
            }
        };

        let bonds_edn: Vec<Edn<'static>> = self
            .ast
            .bonds()
            .iter()
            .enumerate()
            .map(|(i, b)| {
                render_localized(
                    b.src.index(),
                    b.tgt.index(),
                    b.data,
                    i,
                    &self.metadata.bond_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let dative_edn: Vec<Edn<'static>> = self
            .ast
            .dative_bonds()
            .iter()
            .enumerate()
            .map(|(i, v)| {
                render_dative(
                    v.donor.index(),
                    v.acceptor.index(),
                    v.data,
                    i,
                    &self.metadata.dative_bond_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let aromatic_edn: Vec<Edn<'static>> = self
            .ast
            .aromatic_systems()
            .iter()
            .enumerate()
            .map(|(i, v)| {
                render_atoms_map(
                    v.atoms(),
                    i,
                    &self.metadata.aromatic_system_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let multicenter_edn: Vec<Edn<'static>> = self
            .ast
            .multicenter_bonds()
            .iter()
            .enumerate()
            .map(|(i, v)| {
                render_atoms_map(
                    v.atoms(),
                    i,
                    &self.metadata.multicenter_bond_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let noncovalent_edn: Vec<Edn<'static>> = self
            .ast
            .noncovalent_bonds()
            .iter()
            .enumerate()
            .map(|(i, v)| {
                render_noncovalent(
                    v.atoms[0].index(),
                    v.atoms[1].index(),
                    v.data,
                    i,
                    &self.metadata.noncovalent_bond_ids,
                    &render_endpoint,
                )
            })
            .collect();

        let has_aliases = !self.metadata.atom_aliases.is_empty();
        let mut m = EdnMap::with_capacity(10);
        m.insert(Edn::keyword("atoms"), Edn::Vector(atom_elems.into()));
        m.insert(Edn::keyword("bonds"), Edn::Vector(bonds_edn.into()));
        if self.ast.dative_bonds().count() > 0 {
            m.insert(Edn::keyword("dative"), Edn::Vector(dative_edn.into()));
        }
        if self.ast.aromatic_systems().count() > 0 {
            m.insert(Edn::keyword("aromatic"), Edn::Vector(aromatic_edn.into()));
        }
        if self.ast.multicenter_bonds().count() > 0 {
            m.insert(
                Edn::keyword("multicenter"),
                Edn::Vector(multicenter_edn.into()),
            );
        }
        if self.ast.noncovalent_bonds().count() > 0 {
            m.insert(
                Edn::keyword("noncovalent"),
                Edn::Vector(noncovalent_edn.into()),
            );
        }
        for constraint in self.ast.constraints().global() {
            match constraint {
                MoleculeConstraint::TotalCharge(ValueAst::Lit(n)) => {
                    m.insert(Edn::keyword("charge"), Edn::Int(*n));
                }
                MoleculeConstraint::TotalSpin(SpinStateAst::Lit(s)) => {
                    m.insert(Edn::keyword("spin"), Edn::Str(Cow::Owned(s.to_string())));
                }
                _ => {}
            }
        }
        let constraint_entries = render_constraint_list(&self.ast, &self.metadata);
        if !constraint_entries.is_empty() {
            m.insert(
                Edn::keyword("constraints"),
                Edn::Vector(constraint_entries.into()),
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
    source: usize,
    target: usize,
    bond: &BondAst,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let a = render_endpoint(source);
    let bb = render_endpoint(target);
    let bond = bond.to_edn();
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
    donor: usize,
    acceptor: usize,
    bond: &BondAst,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let donor = render_endpoint(donor);
    let acceptor = render_endpoint(acceptor);
    let bond = bond.to_edn();
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
    source: usize,
    target: usize,
    bond: &BondAst,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let a = render_endpoint(source);
    let bb = render_endpoint(target);
    let bond = bond.to_edn();
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
    participants: impl Iterator<Item = AtomIdx>,
    i: usize,
    ids: &IndexMap<usize, String>,
    render_endpoint: &impl Fn(usize) -> Edn<'static>,
) -> Edn<'static> {
    let atom_vec: Vec<Edn<'static>> = participants.map(|a| render_endpoint(a.index())).collect();
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

fn render_id_or_int(ids: &IndexMap<usize, String>, i: usize) -> Edn<'static> {
    match ids.get(&i) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.clone())),
        None => Edn::Int(i as i64),
    }
}

fn render_value_ast(v: &ValueAst) -> Edn<'static> {
    match v {
        ValueAst::Lit(n) => Edn::Int(*n),
        ValueAst::Undetermined => Edn::Keyword(EdnKeyword::owned("undetermined".to_string())),
        ValueAst::LitSet(s) => {
            let elems: Vec<Edn<'static>> = s.iter().map(|n| Edn::Int(*n)).collect();
            Edn::Vector(elems.into())
        }
        ValueAst::Expr(_) => Edn::Nil,
    }
}

fn atom_constraint_has_packed_sugar(c: &AtomConstraint) -> bool {
    matches!(
        c,
        AtomConstraint::Valence(_)
            | AtomConstraint::DonatedPairs(_)
            | AtomConstraint::AcceptedPairs(_)
            | AtomConstraint::MulticenterValence(_)
            | AtomConstraint::AromaticValence(_)
    )
}

fn single_key_map(key: &str, value: Edn<'static>) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(Edn::Keyword(EdnKeyword::owned(key.to_string())), value);
    Edn::Map(m)
}

fn render_atom_constraint_form(c: &AtomConstraint) -> Edn<'static> {
    match c {
        AtomConstraint::Valence(v) => single_key_map("valence", render_value_ast(v)),
        AtomConstraint::AromaticValence(c) => {
            single_key_map("aromatic-valence", render_aromatic_valence_form(c))
        }
        AtomConstraint::MulticenterValence(v) => {
            single_key_map("multicenter-valence", render_value_ast(v))
        }
        AtomConstraint::DonatedPairs(v) => single_key_map("donated-pairs", render_value_ast(v)),
        AtomConstraint::AcceptedPairs(v) => single_key_map("accepted-pairs", render_value_ast(v)),
        AtomConstraint::Degree(v) => single_key_map("degree", render_value_ast(v)),
        AtomConstraint::Connectivity(v) => single_key_map("connectivity", render_value_ast(v)),
        AtomConstraint::TotalHCount(v) => single_key_map("total-h-count", render_value_ast(v)),
        AtomConstraint::InRing => Edn::Keyword(EdnKeyword::owned("in-ring".to_string())),
        AtomConstraint::RingCount(v) => single_key_map("ring-count", render_value_ast(v)),
        AtomConstraint::RingSize(v) => single_key_map("ring-size", render_value_ast(v)),
    }
}

fn render_aromatic_valence_form(c: &AromaticValenceConstraint) -> Edn<'static> {
    match c {
        AromaticValenceConstraint::NotAromatic => {
            Edn::Keyword(EdnKeyword::owned("not-aromatic".to_string()))
        }
        AromaticValenceConstraint::Value(v) => render_value_ast(v),
    }
}

fn render_bond_constraint_form(c: &BondConstraint) -> Edn<'static> {
    match c {
        BondConstraint::RingBond => Edn::Keyword(EdnKeyword::owned("ring-bond".to_string())),
        BondConstraint::Aromatic => Edn::Keyword(EdnKeyword::owned("aromatic".to_string())),
    }
}

fn render_atom_pred(
    idx: AtomIdx,
    c: &AtomConstraint,
    atom_ids: &IndexMap<usize, String>,
) -> Edn<'static> {
    let atom_ref = render_id_or_int(atom_ids, idx.index());
    let form = render_atom_constraint_form(c);
    single_key_map("atom-pred", Edn::Vector(vec![atom_ref, form].into()))
}

fn render_bond_pred(
    idx: BondIdx,
    c: &BondConstraint,
    bond_ids: &IndexMap<usize, String>,
) -> Edn<'static> {
    let bond_ref = render_id_or_int(bond_ids, idx.index());
    let form = render_bond_constraint_form(c);
    single_key_map("bond-pred", Edn::Vector(vec![bond_ref, form].into()))
}

fn render_molecule_constraint(
    c: &MoleculeConstraint,
    metadata: &Metadata,
) -> Edn<'static> {
    match c {
        MoleculeConstraint::AtomPred(idx, inner) => {
            render_atom_pred(*idx, inner, &metadata.atom_ids)
        }
        MoleculeConstraint::BondPred(idx, inner) => {
            render_bond_pred(*idx, inner, &metadata.bond_ids)
        }
        MoleculeConstraint::TotalCharge(v) => single_key_map("total-charge", render_value_ast(v)),
        MoleculeConstraint::TotalSpin(SpinStateAst::Lit(s)) => {
            single_key_map("total-spin", Edn::Str(Cow::Owned(s.to_string())))
        }
        MoleculeConstraint::TotalSpin(_) => single_key_map("total-spin", Edn::Nil),
        MoleculeConstraint::AromaticElectronCount(idx, v) => single_key_map(
            "aromatic-electron-count",
            Edn::Vector(
                vec![
                    render_id_or_int(&metadata.aromatic_system_ids, idx.index()),
                    render_value_ast(v),
                ]
                .into(),
            ),
        ),
        MoleculeConstraint::MulticenterElectronCount(idx, v) => single_key_map(
            "multicenter-electron-count",
            Edn::Vector(
                vec![
                    render_id_or_int(&metadata.multicenter_bond_ids, idx.index()),
                    render_value_ast(v),
                ]
                .into(),
            ),
        ),
        MoleculeConstraint::BondOrderSum(bonds, v) => {
            let bonds_edn: Vec<Edn<'static>> = bonds
                .iter()
                .map(|b| render_id_or_int(&metadata.bond_ids, b.index()))
                .collect();
            let mut inner = EdnMap::with_capacity(2);
            inner.insert(Edn::keyword("bonds"), Edn::Vector(bonds_edn.into()));
            inner.insert(Edn::keyword("equals"), render_value_ast(v));
            single_key_map("bond-order-sum", Edn::Map(inner))
        }
        MoleculeConstraint::Connected(atoms) => {
            let v: Vec<Edn<'static>> = atoms
                .iter()
                .map(|a| render_id_or_int(&metadata.atom_ids, a.index()))
                .collect();
            single_key_map("connected", Edn::Vector(v.into()))
        }
        MoleculeConstraint::SubPattern { anchor, pattern } => {
            let wrapper = MoleculeAstWrapper::from_ast((**pattern).clone());
            let mut inner = EdnMap::with_capacity(2);
            inner.insert(
                Edn::keyword("anchor"),
                render_id_or_int(&metadata.atom_ids, anchor.index()),
            );
            inner.insert(Edn::keyword("pattern"), wrapper.to_edn());
            single_key_map("sub-pattern", Edn::Map(inner))
        }
        MoleculeConstraint::And(xs) => single_key_map(
            "and",
            Edn::Vector(
                xs.iter()
                    .map(|c| render_molecule_constraint(c, metadata))
                    .collect::<Vec<_>>()
                    .into(),
            ),
        ),
        MoleculeConstraint::Or(xs) => single_key_map(
            "or",
            Edn::Vector(
                xs.iter()
                    .map(|c| render_molecule_constraint(c, metadata))
                    .collect::<Vec<_>>()
                    .into(),
            ),
        ),
        MoleculeConstraint::Not(inner) => {
            single_key_map("not", render_molecule_constraint(inner, metadata))
        }
    }
}

fn render_constraint_list(ast: &MoleculeAst, metadata: &Metadata) -> Vec<Edn<'static>> {
    let mut out = Vec::new();
    let constraints = ast.constraints();

    for (idx, set) in constraints.atoms() {
        for c in set.iter() {
            if !atom_constraint_has_packed_sugar(c) {
                out.push(render_atom_pred(*idx, c, &metadata.atom_ids));
            }
        }
    }

    for (idx, set) in constraints.bonds() {
        for c in set.iter() {
            out.push(render_bond_pred(*idx, c, &metadata.bond_ids));
        }
    }

    for (idx, v) in constraints.aromatic_systems() {
        out.push(single_key_map(
            "aromatic-electron-count",
            Edn::Vector(
                vec![
                    render_id_or_int(&metadata.aromatic_system_ids, idx.index()),
                    render_value_ast(v),
                ]
                .into(),
            ),
        ));
    }

    for (idx, v) in constraints.multicenter_bonds() {
        out.push(single_key_map(
            "multicenter-electron-count",
            Edn::Vector(
                vec![
                    render_id_or_int(&metadata.multicenter_bond_ids, idx.index()),
                    render_value_ast(v),
                ]
                .into(),
            ),
        ));
    }

    for c in constraints.global() {
        if matches!(
            c,
            MoleculeConstraint::TotalCharge(_) | MoleculeConstraint::TotalSpin(_)
        ) {
            continue;
        }
        out.push(render_molecule_constraint(c, metadata));
    }

    out
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
    use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::e;
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;

    fn mol_atoms(a: Vec<AtomAst>) -> MoleculeAst {
        MoleculeAst::new(a, vec![], vec![], vec![], vec![], vec![], vec![])
    }

    fn mol_with_bonds(a: Vec<AtomAst>, bonds: Vec<(usize, usize, BondAst)>) -> MoleculeAst {
        let bond_list: Vec<(AtomIdx, AtomIdx, BondAst)> = bonds
            .into_iter()
            .map(|(s, t, b)| (AtomIdx(s as u32), AtomIdx(t as u32), b))
            .collect();
        MoleculeAst::new(a, bond_list, vec![], vec![], vec![], vec![], vec![])
    }

    #[rstest]
    #[case::empty(
        r#"{:atoms [] :bonds []}"#,
        MoleculeAst::default(),
        Metadata::default()
    )]
    #[case::atom(
        r#"{:atoms ["C"] :bonds []}"#,
        mol_atoms(vec![AtomAst::from_element(e!(C))]),
        Metadata::default()
    )]
    #[case::atom_with_id(
        r#"{:atoms [[:C "C"]] :bonds []}"#,
        mol_atoms(vec![AtomAst::from_element(e!(C))]),
        Metadata {
            atom_ids: IndexMap::from([(0, "C".to_string())]),
            ..Default::default()
        }
    )]
    #[case::bond(
        r#"{:atoms ["N" "N"] :bonds [[0 1 :triple]]}"#,
        mol_with_bonds(
            vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))],
            vec![(0, 1, BondAst {
                order: ValueAst::Lit(3),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            })],
        ),
        Metadata::default()
    )]
    #[case::bond_with_ids(
        r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#,
        mol_with_bonds(
            vec![AtomAst::from_element(e!(C)), AtomAst::from_element(e!(O))],
            vec![(0, 1, BondAst {
                order: ValueAst::Lit(1),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            })],
        ),
        Metadata {
            atom_ids: IndexMap::from([(0, "C".to_string()), (1, "O".to_string())]),
            ..Default::default()
        }
    )]
    #[case::bond_id(
        r#"{:atoms ["H" "F"] :bonds [{:id :b1 :a 0 :b 1 :bond :single}]}"#,
        mol_with_bonds(
            vec![AtomAst::from_element(e!(H)), AtomAst::from_element(e!(F))],
            vec![(0, 1, BondAst {
                order: ValueAst::Lit(1),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            })],
        ),
        Metadata {
            bond_ids: IndexMap::from([(0, "b1".to_string())]),
            ..Default::default()
        }
    )]
    #[case::charge(
        r#"{:atoms [[:F "F#c-"]] :bonds [] :charge -1}"#,
        {
            let mut ast = mol_atoms(vec![AtomAst {
                element: ElementAst::Lit(Element::F),
                isotope_mass: IsotopeAst::Undetermined,
                implicit_hydrogens: HydrogenAst::Undetermined,
                charge: ValueAst::Lit(-1),
                lone_pairs: ValueAst::Undetermined,
                spin: SpinStateAst::default(),
            }]);
            ast.constraints_mut().insert(MoleculeConstraint::TotalCharge(ValueAst::Lit(-1)));
            ast
        },
        Metadata {
            atom_ids: IndexMap::from([(0, "F".to_string())]),
            ..Default::default()
        }
    )]
    #[case::alias_indexed(
        r#"{:atoms [:ch] :bonds [] :aliases [:ch "C #h1"]}"#,
        mol_atoms(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Undetermined,
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(1)),
            charge: ValueAst::Undetermined,
            lone_pairs: ValueAst::Undetermined,
            spin: SpinStateAst::default(),
        }]),
        Metadata {
            atom_aliases: BiMap::from_iter([(
                "ch".to_string(),
                AtomPattern::new(AtomAst {
                    element: ElementAst::Lit(Element::C),
                    isotope_mass: IsotopeAst::Undetermined,
                    implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(1)),
                    charge: ValueAst::Undetermined,
                    lone_pairs: ValueAst::Undetermined,
                    spin: SpinStateAst::default(),
                }),
            )]),
            ..Default::default()
        }
    )]
    #[case::alias_reused(
        r#"{:atoms [:n :n] :bonds [[0 1 :single]] :aliases [:n "N"]}"#,
        mol_with_bonds(
            vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))],
            vec![(0, 1, BondAst {
                order: ValueAst::Lit(1),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            })],
        ),
        Metadata {
            atom_aliases: BiMap::from_iter([(
                "n".to_string(),
                AtomPattern::new(AtomAst::from_element(e!(N))),
            )]),
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
        assert_eq!(dsl.ast.dative_bonds().count(), 1);
        let view = dsl.ast.dative_bonds().iter().next().unwrap();
        assert_eq!(view.donor, AtomIdx(1)); // donor = N
        assert_eq!(view.acceptor, AtomIdx(0)); // acceptor = B
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
        assert_eq!(dsl.ast.aromatic_systems().count(), 1);
        assert_eq!(
            dsl.metadata.aromatic_system_ids.get(&0),
            Some(&"ar1".to_string())
        );
        let view = dsl.ast.aromatic_systems().iter().next().unwrap();
        let p: Vec<AtomIdx> = view.atoms().collect();
        assert_eq!(
            p,
            vec![
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2),
                AtomIdx(3),
                AtomIdx(4),
                AtomIdx(5)
            ]
        );
    }

    #[rstest]
    #[case::non_map("3")]
    #[case::missing_atoms(r#"{:bonds []}"#)]
    #[case::missing_bonds(r#"{:atoms ["C"]}"#)]
    #[case::unknown_endpoint(r#"{:atoms [[:C "C"]] :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#)]
    #[case::bad_atom_string(r##"{:atoms ["#h3"] :bonds []}"##)]
    #[case::trailing_content(r#"{:atoms ["C"] :bonds []} :extra :junk"#)]
    #[case::duplicate_id(r#"{:atoms [[:C "C"] [:C "O"]] :bonds []}"#)]
    #[case::duplicate_alias(r#"{:aliases [:ch "C #h1" :ch "C #h2"] :atoms [] :bonds []}"#)]
    fn test_parse_molecule_dsl_invalid(#[case] input: &str) {
        assert!(
            parse_molecule_dsl(input).is_err(),
            "{input:?} should fail but succeeded"
        );
    }

    #[rstest]
    #[case::plain_vector(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::with_id_vector(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
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
    fn test_alias_id_disjointness_error() {
        let result = parse_molecule_dsl(r#"{:atoms [[:ch "C"]] :bonds [] :aliases [:ch "C #h1"]}"#);
        assert!(result.is_err(), "alias and id with same name should fail");
    }

    #[rstest]
    #[case::plain(r#"{:atoms ["C" "O"] :bonds [[0 1 :single]]}"#)]
    #[case::with_id(r#"{:atoms [[:C "C"] [:O "O"]] :bonds [[:C :O :single]]}"#)]
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

    #[rstest]
    #[case::total_charge(
        r#"{:atoms ["C"] :bonds [] :constraints [{:total-charge 0}]}"#
    )]
    #[case::total_spin(
        r##"{:atoms ["C"] :bonds [] :constraints [{:total-spin "#u1"}]}"##
    )]
    #[case::atom_pred_degree(
        r#"{:atoms ["C" "C"] :bonds [[0 1 :single]] :constraints [{:atom-pred [0 {:degree 3}]}]}"#
    )]
    #[case::atom_pred_degree_by_id(
        r#"{:atoms [[:c1 "C"] "C"] :bonds [[0 1 :single]] :constraints [{:atom-pred [:c1 {:degree 3}]}]}"#
    )]
    #[case::atom_pred_connectivity(
        r#"{:atoms ["C"] :bonds [] :constraints [{:atom-pred [0 {:connectivity 4}]}]}"#
    )]
    #[case::atom_pred_total_h_count(
        r#"{:atoms ["C"] :bonds [] :constraints [{:atom-pred [0 {:total-h-count 2}]}]}"#
    )]
    #[case::atom_pred_in_ring(
        r#"{:atoms ["C"] :bonds [] :constraints [{:atom-pred [0 :in-ring]}]}"#
    )]
    #[case::atom_pred_ring_count(
        r#"{:atoms ["C"] :bonds [] :constraints [{:atom-pred [0 {:ring-count 1}]}]}"#
    )]
    #[case::atom_pred_ring_size(
        r#"{:atoms ["C"] :bonds [] :constraints [{:atom-pred [0 {:ring-size 6}]}]}"#
    )]
    #[case::bond_pred_ring(
        r#"{:atoms ["C" "C"] :bonds [[0 1 :single]] :constraints [{:bond-pred [0 :ring-bond]}]}"#
    )]
    #[case::bond_pred_aromatic(
        r#"{:atoms ["C" "C"] :bonds [[0 1 :single]] :constraints [{:bond-pred [0 :aromatic]}]}"#
    )]
    #[case::aromatic_electron_count(
        r#"{:atoms ["C" "C"] :bonds [[0 1 :single]] :aromatic [{:id :ar1 :atoms [0 1]}] :constraints [{:aromatic-electron-count [:ar1 6]}]}"#
    )]
    #[case::multicenter_electron_count(
        r#"{:atoms ["C" "C"] :bonds [] :multicenter [{:id :mc1 :atoms [0 1]}] :constraints [{:multicenter-electron-count [:mc1 2]}]}"#
    )]
    #[case::bond_order_sum(
        r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single]] :constraints [{:bond-order-sum {:bonds [0 1] :equals 2}}]}"#
    )]
    #[case::connected(
        r#"{:atoms [[:a1 "C"] [:a2 "C"]] :bonds [[0 1 :single]] :constraints [{:connected [:a1 :a2]}]}"#
    )]
    #[case::sub_pattern(
        r#"{:atoms [[:a "C"]] :bonds [] :constraints [{:sub-pattern {:anchor :a :pattern {:atoms ["C"] :bonds []}}}]}"#
    )]
    #[case::and_combinator(
        r#"{:atoms ["C"] :bonds [] :constraints [{:and [{:total-charge 0} {:atom-pred [0 :in-ring]}]}]}"#
    )]
    #[case::or_combinator(
        r#"{:atoms ["C"] :bonds [] :constraints [{:or [{:atom-pred [0 {:degree 3}]} {:atom-pred [0 {:degree 4}]}]}]}"#
    )]
    #[case::not_combinator(
        r#"{:atoms ["C"] :bonds [] :constraints [{:not {:atom-pred [0 :in-ring]}}]}"#
    )]
    fn test_constraints_canonical_roundtrip(#[case] input: &str) {
        let dsl1 = parse_molecule_dsl(input).unwrap();
        let edn = dsl1.to_string();
        let dsl2 = parse_molecule_dsl(&edn).unwrap();
        assert_eq!(dsl1, dsl2);
    }
}
