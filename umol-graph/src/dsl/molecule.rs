//! Molecule DSL definitions

use std::borrow::Cow;
use std::collections::HashSet;
use std::str::FromStr;
use std::{fmt, mem};

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{
    DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, FromEdnMap,
    ToEdn, ToEdnMap,
};
use umol_shared::spin::SpinState;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::error::ParseError;
use crate::api::pattern::{AtomPattern, BondPattern};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::{
    AromaticValenceConstraint, AtomConstraint, BondConstraint, MoleculeConstraint,
};
use crate::ast::config::MoleculeAstConfig;
use crate::ast::error::LoweringError;
use crate::ast::molecule::MoleculeAst;
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::{AromaticSystemIdx, AtomIdx, BondIdx, FromAst, MulticenterBondIdx, ToAst};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_ids: IndexMap<usize, String>,
    pub atom_aliases: BiMap<String, AtomDsl>,
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
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    map: MoleculeMapDsl,
}

impl MoleculeDsl {
    pub fn new(ast: MoleculeAst, metadata: Metadata) -> Self {
        Self::from_parts(ast, metadata)
    }

    pub fn from_ast(ast: MoleculeAst) -> Self {
        Self::from_parts(ast, Metadata::default())
    }

    pub fn from_parts(ast: MoleculeAst, metadata: Metadata) -> Self {
        let Edn::Map(map) = render_molecule_edn(&ast, &metadata) else {
            unreachable!("molecule DSL rendering always returns a map")
        };
        let map = parse_molecule_map_edn(&map)
            .expect("rendered molecule DSL should always parse back to MoleculeMapDsl");
        Self { map }
    }

    pub fn lower_parts(&self) -> Result<(MoleculeAst, Metadata), LoweringError> {
        self.map
            .clone()
            .lower_parts()
            .map_err(|e| LoweringError::Custom(e.to_string()))
    }

    pub fn lower_ast(&self) -> Result<MoleculeAst, LoweringError> {
        self.lower_parts().map(|(ast, _)| ast)
    }
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    fn from_ast(ast: &MoleculeAst, _cfg: &MoleculeAstConfig) -> Result<Self, LoweringError> {
        Ok(Self::from_ast(ast.clone()))
    }
}

impl ToAst<MoleculeAst> for MoleculeDsl {
    fn to_ast(&self, _cfg: &MoleculeAstConfig) -> Result<MoleculeAst, LoweringError> {
        self.lower_ast()
    }
}

impl PartialEq for MoleculeDsl {
    fn eq(&self, other: &Self) -> bool {
        self.lower_parts().ok() == other.lower_parts().ok()
    }
}

impl Eq for MoleculeDsl {}

/// Parse a molecule DSL EDN string via the fused single-pass parser.
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeDsl, ParseError> {
    let mut de = EdnStreamDeserializer::new(input);
    let mol_input =
        read_molecule_input(&mut de).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    de.expect_eof()
        .map_err(|e| ParseError::EdnParse(e.to_string()))?;
    let dsl = MoleculeDsl { map: mol_input };
    dsl.lower_parts()
        .map_err(|e| ParseError::InvalidValue(e.to_string()))?;
    Ok(dsl)
}

impl FromStr for MoleculeDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Inline here -- no need for free fn + trait impl
        parse_molecule_dsl(s)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AtomEntryDsl {
    Str(String),
    WithId(String, Box<AtomEntryDsl>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AtomRefDsl {
    Index(usize),
    Id(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CovalentBondEntryDsl {
    id: Option<String>,
    a: AtomRefDsl,
    b: AtomRefDsl,
    bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct DativeBondEntryDsl {
    id: Option<String>,
    donor: AtomRefDsl,
    acceptor: AtomRefDsl,
    bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct AromaticEntryDsl {
    id: Option<String>,
    atoms: Vec<AtomRefDsl>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MulticenterEntryDsl {
    id: Option<String>,
    atoms: Vec<AtomRefDsl>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct NoncovalentBondEntryDsl {
    id: Option<String>,
    a: AtomRefDsl,
    b: AtomRefDsl,
    bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MoleculeConstraintDsl(Edn<'static>);

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct MoleculeMapDsl {
    atoms: Vec<AtomEntryDsl>,
    bonds: Vec<CovalentBondEntryDsl>,
    dative_bonds: Vec<DativeBondEntryDsl>,
    aromatic_systems: Vec<AromaticEntryDsl>,
    multicenter_bonds: Vec<MulticenterEntryDsl>,
    noncovalent_bonds: Vec<NoncovalentBondEntryDsl>,
    atom_aliases: Vec<String>,
    constraints: Vec<MoleculeConstraintDsl>,
}

impl MoleculeMapDsl {
    fn lower_parts(self) -> Result<(MoleculeAst, Metadata), ParseError> {
        let alias_table = Self::build_alias_table(&self.atom_aliases)?;

        let mut atoms: Vec<AtomAst> = Vec::with_capacity(self.atoms.len());
        let mut lifted_constraints: Vec<MoleculeConstraint> = Vec::new();
        let mut atom_ids: IndexMap<usize, String> = IndexMap::new();
        let mut id_to_index: IndexMap<String, usize> = IndexMap::new();
        let mut atom_aliases: BiMap<String, AtomDsl> = BiMap::new();

        for entry in self.atoms {
            let pos = atoms.len();
            let (id, atom_dsl) = Self::resolve_entry(entry, &alias_table)?;
            let (atom_ast, atom_constraints) = atom_dsl.lower_parts()?;
            if let Some(id_name) = id {
                if id_to_index.contains_key(&id_name) || alias_table.contains_key(&id_name) {
                    return Err(ParseError::DuplicateId(id_name));
                }
                id_to_index.insert(id_name.clone(), pos);
                atom_ids.insert(pos, id_name);
            }
            for c in atom_constraints {
                lifted_constraints.push(MoleculeConstraint::AtomPred(AtomIdx(pos as u32), c));
            }
            atoms.push(atom_ast);
        }

        for (name, def) in &alias_table {
            atom_aliases.insert(name.clone(), def.clone());
        }

        let atom_count = atoms.len();
        let resolve = |r: &AtomRefDsl| -> Result<AtomIdx, ParseError> {
            match r {
                AtomRefDsl::Index(i) => {
                    if *i < atom_count {
                        Ok(AtomIdx(*i as u32))
                    } else {
                        Err(ParseError::InvalidAtomIndex(i.to_string()))
                    }
                }
                AtomRefDsl::Id(name) => id_to_index
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
            let (bond_ast, bond_constraints) = b.bond.lower_parts()?;
            for c in bond_constraints {
                lifted_constraints.push(MoleculeConstraint::BondPred(BondIdx(i as u32), c));
            }
            bond_list.push((a, bb, bond_ast));
        }

        let mut dative_list = Vec::new();
        let mut dative_bond_ids = IndexMap::new();
        for (i, db) in self.dative_bonds.into_iter().enumerate() {
            let donor = resolve(&db.donor)?;
            let acceptor = resolve(&db.acceptor)?;
            if let Some(id) = check_id(db.id)? {
                dative_bond_ids.insert(i, id);
            }
            dative_list.push((donor, acceptor, db.bond.into_ast()?));
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
            noncovalent_list.push((a, bb, nc.bond.into_ast()?));
        }

        let mut constraints = lifted_constraints;
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
                constraints.push(entry.lower(&resolver)?);
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
        Ok((ast, metadata))
    }

    fn build_alias_table(raw: &[String]) -> Result<IndexMap<String, AtomDsl>, ParseError> {
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
            table.insert(name.clone(), def.parse()?);
        }
        Ok(table)
    }

    fn resolve_entry(
        entry: AtomEntryDsl,
        alias_table: &IndexMap<String, AtomDsl>,
    ) -> Result<(Option<String>, AtomDsl), ParseError> {
        match entry {
            AtomEntryDsl::Str(s) => {
                if let Some(def) = alias_table.get(&s) {
                    Ok((None, def.clone()))
                } else {
                    Ok((None, s.parse()?))
                }
            }
            AtomEntryDsl::WithId(id, inner) => {
                let (_, atom_dsl) = Self::resolve_entry(*inner, alias_table)?;
                Ok((Some(id), atom_dsl))
            }
        }
    }
}

impl MoleculeConstraintDsl {
    fn from_edn(edn: &Edn<'_>) -> Self {
        Self(edn.to_edn())
    }

    fn lower(&self, resolver: &ConstraintResolver<'_>) -> Result<MoleculeConstraint, ParseError> {
        parse_molecule_constraint(&self.0, resolver)
    }

    fn to_edn(&self) -> Edn<'static> {
        self.0.clone()
    }
}

impl ToAst<MoleculeAst> for MoleculeMapDsl {
    fn to_ast(&self, _cfg: &MoleculeAstConfig) -> Result<MoleculeAst, LoweringError> {
        self.clone()
            .lower_parts()
            .map(|(ast, _)| ast)
            .map_err(|e| LoweringError::Custom(e.to_string()))
    }
}

fn parse_keyword_or_string(edn: &Edn<'_>, context: &'static str) -> Result<String, ParseError> {
    match edn {
        Edn::Keyword(k) => Ok(k.name().to_string()),
        Edn::Str(s) => Ok(s.to_string()),
        other => Err(ParseError::InvalidValue(format!(
            "{context}: expected keyword or string, got {other}"
        ))),
    }
}

fn parse_keyword_id(edn: &Edn<'_>, context: &'static str) -> Result<String, ParseError> {
    match edn {
        Edn::Keyword(k) => Ok(k.name().to_string()),
        other => Err(ParseError::InvalidValue(format!(
            "{context}: expected keyword id, got {other}"
        ))),
    }
}

fn parse_vector<'e, T, F>(
    edn: &'e Edn<'e>,
    context: &'static str,
    mut element: F,
) -> Result<Vec<T>, ParseError>
where
    F: FnMut(&'e Edn<'e>) -> Result<T, ParseError>,
{
    let Edn::Vector(v) = edn else {
        return Err(ParseError::InvalidValue(format!(
            "{context}: expected vector, got {edn}"
        )));
    };
    v.iter().map(&mut element).collect()
}

fn parse_atom_entry_edn(edn: &Edn<'_>) -> Result<AtomEntryDsl, ParseError> {
    match edn {
        Edn::Str(_) | Edn::Keyword(_) => Ok(AtomEntryDsl::Str(parse_keyword_or_string(
            edn,
            "atom entry",
        )?)),
        Edn::Vector(v) if v.len() == 2 => Ok(AtomEntryDsl::WithId(
            parse_keyword_or_string(&v[0], "atom entry id")?,
            Box::new(parse_atom_entry_edn(&v[1])?),
        )),
        other => Err(ParseError::InvalidValue(format!("atom entry: {other}"))),
    }
}

fn parse_atom_ref_edn(edn: &Edn<'_>) -> Result<AtomRefDsl, ParseError> {
    match edn {
        Edn::Keyword(k) => Ok(AtomRefDsl::Id(k.name().to_string())),
        Edn::Str(s) => Ok(AtomRefDsl::Id(s.to_string())),
        Edn::Int(n) => {
            let idx = usize::try_from(*n).map_err(|_| ParseError::InvalidValue(format!(
                "atom ref index {n} out of range"
            )))?;
            Ok(AtomRefDsl::Index(idx))
        }
        other => Err(ParseError::InvalidValue(format!("atom ref: {other}"))),
    }
}

fn parse_bond_spec_edn(edn: &Edn<'_>) -> Result<BondDsl, ParseError> {
    match edn {
        Edn::Keyword(_) | Edn::Str(_) => BondDsl::from_edn(edn).map_err(|e| match e {
            DeError::Subgrammar { message, .. } => ParseError::InvalidBondSpec(message),
            other => ParseError::InvalidValue(format!("bond spec: {other}")),
        }),
        other => Err(ParseError::InvalidValue(format!("bond spec: {other}"))),
    }
}

fn parse_endpoint_bond_map_edn(
    map: &EdnMap<'_>,
    kind: EndpointBondKind,
) -> Result<(Option<String>, AtomRefDsl, AtomRefDsl, BondDsl), ParseError> {
    let first_key = kind.first_key();
    let second_key = kind.second_key();
    let mut id = None;
    let mut a = None;
    let mut b = None;
    let mut bond = None;

    for (key, value) in map.iter() {
        let Edn::Keyword(k) = key else {
            return Err(ParseError::InvalidValue(format!(
                "bond entry key must be keyword, got {key}"
            )));
        };
        match k.name() {
            "id" => id = Some(parse_keyword_id(value, "bond entry :id")?),
            name if name == first_key => a = Some(parse_atom_ref_edn(value)?),
            name if name == second_key => b = Some(parse_atom_ref_edn(value)?),
            "bond" => bond = Some(parse_bond_spec_edn(value)?),
            _ => {}
        }
    }

    let a = a.ok_or_else(|| ParseError::MissingKey(format!(":{first_key}")))?;
    let b = b.ok_or_else(|| ParseError::MissingKey(format!(":{second_key}")))?;
    let bond = bond.ok_or_else(|| ParseError::MissingKey(":bond".to_string()))?;
    Ok((id, a, b, bond))
}

fn parse_atoms_map_edn(map: &EdnMap<'_>) -> Result<(Option<String>, Vec<AtomRefDsl>), ParseError> {
    let mut id = None;
    let mut atoms = None;

    for (key, value) in map.iter() {
        let Edn::Keyword(k) = key else {
            return Err(ParseError::InvalidValue(format!(
                "relation entry key must be keyword, got {key}"
            )));
        };
        match k.name() {
            "id" => id = Some(parse_keyword_id(value, "relation entry :id")?),
            "atoms" => atoms = Some(parse_vector(value, "relation entry :atoms", parse_atom_ref_edn)?),
            _ => {}
        }
    }

    let atoms = atoms.ok_or_else(|| ParseError::MissingKey(":atoms".to_string()))?;
    Ok((id, atoms))
}

fn parse_constraint_entry_edn(edn: &Edn<'_>) -> MoleculeConstraintDsl {
    MoleculeConstraintDsl::from_edn(edn)
}

fn parse_molecule_map_edn(map: &EdnMap<'_>) -> Result<MoleculeMapDsl, ParseError> {
    let mut atoms: Option<Vec<AtomEntryDsl>> = None;
    let mut bonds: Option<Vec<CovalentBondEntryDsl>> = None;
    let mut dative_bonds: Vec<DativeBondEntryDsl> = Vec::new();
    let mut aromatic_systems: Vec<AromaticEntryDsl> = Vec::new();
    let mut multicenter_bonds: Vec<MulticenterEntryDsl> = Vec::new();
    let mut noncovalent_bonds: Vec<NoncovalentBondEntryDsl> = Vec::new();
    let mut atom_aliases: Vec<String> = Vec::new();
    let mut constraints: Vec<MoleculeConstraintDsl> = Vec::new();

    for (key, value) in map.iter() {
        let Edn::Keyword(k) = key else {
            return Err(ParseError::InvalidValue(format!(
                "molecule key must be keyword, got {key}"
            )));
        };
        match k.name() {
            "atoms" => atoms = Some(parse_vector(value, "atoms", parse_atom_entry_edn)?),
            "bonds" => {
                bonds = Some(parse_vector(value, "bonds", |entry| match entry {
                    Edn::Vector(v) if v.len() == 3 => Ok(CovalentBondEntryDsl {
                        id: None,
                        a: parse_atom_ref_edn(&v[0])?,
                        b: parse_atom_ref_edn(&v[1])?,
                        bond: parse_bond_spec_edn(&v[2])?,
                    }),
                    Edn::Map(m) => {
                        let (id, a, b, bond) =
                            parse_endpoint_bond_map_edn(m, EndpointBondKind::Localized)?;
                        Ok(CovalentBondEntryDsl { id, a, b, bond })
                    }
                    other => Err(ParseError::InvalidValue(format!("bond entry: {other}"))),
                })?)
            }
            "dative" => {
                dative_bonds = parse_vector(value, "dative", |entry| {
                    let Edn::Map(m) = entry else {
                        return Err(ParseError::InvalidValue(format!(
                            "dative bond entry: {entry}"
                        )));
                    };
                    let (id, donor, acceptor, bond) =
                        parse_endpoint_bond_map_edn(m, EndpointBondKind::Dative)?;
                    Ok(DativeBondEntryDsl {
                        id,
                        donor,
                        acceptor,
                        bond: bond,
                    })
                })?
            }
            "aromatic" => {
                aromatic_systems = parse_vector(value, "aromatic", |entry| {
                    let Edn::Map(m) = entry else {
                        return Err(ParseError::InvalidValue(format!("aromatic entry: {entry}")));
                    };
                    let (id, atoms) = parse_atoms_map_edn(m)?;
                    Ok(AromaticEntryDsl { id, atoms })
                })?
            }
            "multicenter" => {
                multicenter_bonds = parse_vector(value, "multicenter", |entry| {
                    let Edn::Map(m) = entry else {
                        return Err(ParseError::InvalidValue(format!(
                            "multicenter entry: {entry}"
                        )));
                    };
                    let (id, atoms) = parse_atoms_map_edn(m)?;
                    Ok(MulticenterEntryDsl { id, atoms })
                })?
            }
            "noncovalent" => {
                noncovalent_bonds = parse_vector(value, "noncovalent", |entry| {
                    let Edn::Map(m) = entry else {
                        return Err(ParseError::InvalidValue(format!(
                            "noncovalent bond entry: {entry}"
                        )));
                    };
                    let (id, a, b, bond) =
                        parse_endpoint_bond_map_edn(m, EndpointBondKind::Noncovalent)?;
                    Ok(NoncovalentBondEntryDsl {
                        id,
                        a,
                        b,
                        bond,
                    })
                })?
            }
            "atom-aliases" | "aliases" => {
                atom_aliases = parse_vector(value, "atom-aliases", |entry| {
                    parse_keyword_or_string(entry, "atom-aliases entry")
                })?
            }
            "constraints" => {
                constraints =
                    parse_vector(value, "constraints", |entry| Ok(parse_constraint_entry_edn(entry)))?
            }
            _ => {}
        }
    }

    let atoms = atoms.ok_or_else(|| ParseError::MissingKey(":atoms".to_string()))?;
    let bonds = bonds.ok_or_else(|| ParseError::MissingKey(":bonds".to_string()))?;

    Ok(MoleculeMapDsl {
        atoms,
        bonds,
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        atom_aliases,
        constraints,
    })
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
            let target_anchor_edn = m
                .get_keyword("target-anchor")
                .ok_or_else(|| ParseError::MissingKey(":target-anchor".to_string()))?;
            let pattern_anchor_edn = m
                .get_keyword("pattern-anchor")
                .ok_or_else(|| ParseError::MissingKey(":pattern-anchor".to_string()))?;
            let pattern_edn = m
                .get_keyword("pattern")
                .ok_or_else(|| ParseError::MissingKey(":pattern".to_string()))?;
            let target_anchor: AtomIdx =
                resolve_ref(target_anchor_edn, r.atom_count, r.atom_ids, "atom")?;
            let wrapper = MoleculeDsl::from_edn(pattern_edn)
                .map_err(|e| ParseError::InvalidValue(format!("sub-pattern: {e}")))?;
            let (pattern_ast, pattern_metadata) = wrapper
                .lower_parts()
                .map_err(|e| ParseError::InvalidValue(format!("sub-pattern: {e}")))?;
            let pattern_anchor: AtomIdx = resolve_ref(
                pattern_anchor_edn,
                pattern_ast.atoms().count(),
                &pattern_metadata.atom_ids,
                "atom",
            )?;
            Ok(MoleculeConstraint::SubPattern {
                target_anchor,
                pattern_anchor,
                pattern: Box::new(pattern_ast),
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

fn read_molecule_input(de: &mut EdnStreamDeserializer<'_>) -> Result<MoleculeMapDsl, EdnError> {
    de.consume_byte(b'{')?;

    let mut atoms: Option<Vec<AtomEntryDsl>> = None;
    let mut bonds: Option<Vec<CovalentBondEntryDsl>> = None;
    let mut dative_bonds: Vec<DativeBondEntryDsl> = Vec::new();
    let mut aromatic_systems: Vec<AromaticEntryDsl> = Vec::new();
    let mut multicenter_bonds: Vec<MulticenterEntryDsl> = Vec::new();
    let mut noncovalent_bonds: Vec<NoncovalentBondEntryDsl> = Vec::new();
    let mut atom_aliases: Vec<String> = Vec::new();
    let mut constraints: Vec<MoleculeConstraintDsl> = Vec::new();

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

    Ok(MoleculeMapDsl {
        atoms,
        bonds,
        dative_bonds,
        aromatic_systems,
        multicenter_bonds,
        noncovalent_bonds,
        atom_aliases,
        constraints,
    })
}

fn read_constraint_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MoleculeConstraintDsl, EdnError> {
    let slice = de.read_value_slice()?;
    Edn::from_str(slice).map(MoleculeConstraintDsl).map_err(Into::into)
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

fn read_atom_entry(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomEntryDsl, EdnError> {
    match de.peek_byte()? {
        Some(b'"') | Some(b':') => {
            let s = de.read_string_or_keyword()?;
            Ok(AtomEntryDsl::Str(s.into_owned()))
        }
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let id = de.read_string_or_keyword()?.into_owned();
            let inner = read_atom_entry(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryDsl::WithId(id, Box::new(inner)))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_atom_ref(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomRefDsl, EdnError> {
    match de.peek_byte()? {
        Some(b':') => Ok(AtomRefDsl::Id(de.read_keyword_name()?.into_owned())),
        Some(b'"') => Ok(AtomRefDsl::Id(de.read_string()?.into_owned())),
        Some(b) if b.is_ascii_digit() || b == b'-' || b == b'+' => {
            let n = de.read_i64()?;
            let idx = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                value: n.to_string(),
                target: "atom index",
                path: Vec::new(),
            })?;
            Ok(AtomRefDsl::Index(idx))
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_bond_spec(de: &mut EdnStreamDeserializer<'_>) -> Result<BondDsl, EdnError> {
    let s = de.read_string_or_keyword()?;
    let aliases = super::bond::builtin_bond_aliases();
    if let Some(ast) = aliases.get_by_left(s.as_ref()) {
        return Ok(BondDsl::from_pattern(BondPattern::new(ast.clone())));
    }
    s.as_ref()
        .parse::<BondDsl>()
        .map_err(|e| DeError::subgrammar("bond", e).into())
}

fn read_localized_bond(de: &mut EdnStreamDeserializer<'_>) -> Result<CovalentBondEntryDsl, EdnError> {
    match de.peek_byte()? {
        Some(b'[') => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_spec(de)?;
            de.consume_byte(b']')?;
            Ok(CovalentBondEntryDsl {
                id: None,
                a,
                b,
                bond,
            })
        }
        Some(b'{') => {
            let (id, a, b, bond) = read_endpoint_bond_map(de, EndpointBondKind::Localized)?;
            Ok(CovalentBondEntryDsl { id, a, b, bond })
        }
        other => Err(unexpected(de.position(), other)),
    }
}

fn read_dative_bond(de: &mut EdnStreamDeserializer<'_>) -> Result<DativeBondEntryDsl, EdnError> {
    let (id, donor, acceptor, bond) = read_endpoint_bond_map(de, EndpointBondKind::Dative)?;
    Ok(DativeBondEntryDsl {
        id,
        donor,
        acceptor,
        bond,
    })
}

fn read_noncovalent_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondEntryDsl, EdnError> {
    let (id, a, b, bond) = read_endpoint_bond_map(de, EndpointBondKind::Noncovalent)?;
    Ok(NoncovalentBondEntryDsl {
        id,
        a,
        b,
        bond,
    })
}

fn read_aromatic_system(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticEntryDsl, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(AromaticEntryDsl { id, atoms })
}

fn read_multicenter_bond(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterEntryDsl, EdnError> {
    let (id, atoms) = read_atoms_bond_map(de)?;
    Ok(MulticenterEntryDsl { id, atoms })
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
) -> Result<(Option<String>, AtomRefDsl, AtomRefDsl, BondDsl), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<String> = None;
    let mut a: Option<AtomRefDsl> = None;
    let mut b: Option<AtomRefDsl> = None;
    let mut bond: Option<BondDsl> = None;
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
) -> Result<(Option<String>, Vec<AtomRefDsl>), EdnError> {
    de.consume_byte(b'{')?;
    let mut id: Option<String> = None;
    let mut atoms: Option<Vec<AtomRefDsl>> = None;
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

impl<'de> FromEdnMap<'de> for MoleculeDsl {
    fn from_edn_map(map: &EdnMap<'de>) -> Result<Self, DeError> {
        let raw = parse_molecule_map_edn(map).map_err(|e| DeError::subgrammar("molecule", e))?;
        let dsl = Self { map: raw };
        dsl.lower_parts().map_err(|e| DeError::subgrammar("molecule", e))?;
        Ok(dsl)
    }
}

impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(map) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        Self::from_edn_map(map)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        parse_molecule_dsl(input).map_err(|e| DeError::subgrammar("molecule", e).into())
    }
}

impl ToEdnMap for MoleculeDsl {
    fn to_edn_map(&self) -> EdnMap<'static> {
        render_molecule_map_dsl(&self.map)
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

fn render_atom_ref_dsl(atom_ref: &AtomRefDsl) -> Edn<'static> {
    match atom_ref {
        AtomRefDsl::Index(i) => Edn::Int(*i as i64),
        AtomRefDsl::Id(id) => Edn::Keyword(EdnKeyword::owned(id.clone())),
    }
}

fn render_atom_entry_dsl(entry: &AtomEntryDsl, alias_names: &HashSet<&str>) -> Edn<'static> {
    match entry {
        AtomEntryDsl::Str(s) => {
            if alias_names.contains(s.as_str()) {
                Edn::Keyword(EdnKeyword::owned(s.clone()))
            } else {
                Edn::Str(Cow::Owned(s.clone()))
            }
        }
        AtomEntryDsl::WithId(id, inner) => Edn::Vector(
            vec![
                Edn::Keyword(EdnKeyword::owned(id.clone())),
                render_atom_entry_dsl(inner, alias_names),
            ]
            .into(),
        ),
    }
}

fn render_covalent_bond_entry_dsl(entry: &CovalentBondEntryDsl) -> Edn<'static> {
    let a = render_atom_ref_dsl(&entry.a);
    let b = render_atom_ref_dsl(&entry.b);
    let bond = entry.bond.to_edn();
    if let Some(id) = &entry.id {
        let mut m = EdnMap::with_capacity(4);
        m.insert(Edn::keyword("id"), Edn::Keyword(EdnKeyword::owned(id.clone())));
        m.insert(Edn::keyword("a"), a);
        m.insert(Edn::keyword("b"), b);
        m.insert(Edn::keyword("bond"), bond);
        Edn::Map(m)
    } else {
        Edn::Vector(vec![a, b, bond].into())
    }
}

fn render_dative_bond_entry_dsl(entry: &DativeBondEntryDsl) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = &entry.id {
        m.insert(Edn::keyword("id"), Edn::Keyword(EdnKeyword::owned(id.clone())));
    }
    m.insert(Edn::keyword("donor"), render_atom_ref_dsl(&entry.donor));
    m.insert(
        Edn::keyword("acceptor"),
        render_atom_ref_dsl(&entry.acceptor),
    );
    m.insert(Edn::keyword("bond"), entry.bond.to_edn());
    Edn::Map(m)
}

fn render_atoms_entry_dsl(id: &Option<String>, atoms: &[AtomRefDsl]) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    if let Some(id) = id {
        m.insert(Edn::keyword("id"), Edn::Keyword(EdnKeyword::owned(id.clone())));
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(atoms.iter().map(render_atom_ref_dsl).collect::<Vec<_>>().into()),
    );
    Edn::Map(m)
}

fn render_noncovalent_bond_entry_dsl(entry: &NoncovalentBondEntryDsl) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(id) = &entry.id {
        m.insert(Edn::keyword("id"), Edn::Keyword(EdnKeyword::owned(id.clone())));
    }
    m.insert(Edn::keyword("a"), render_atom_ref_dsl(&entry.a));
    m.insert(Edn::keyword("b"), render_atom_ref_dsl(&entry.b));
    m.insert(Edn::keyword("bond"), entry.bond.to_edn());
    Edn::Map(m)
}

fn render_molecule_map_dsl(map_dsl: &MoleculeMapDsl) -> EdnMap<'static> {
    let alias_names: HashSet<&str> = map_dsl
        .atom_aliases
        .chunks_exact(2)
        .map(|pair| pair[0].as_str())
        .collect();
    let mut m = EdnMap::with_capacity(10);
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(
            map_dsl
                .atoms
                .iter()
                .map(|entry| render_atom_entry_dsl(entry, &alias_names))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(
        Edn::keyword("bonds"),
        Edn::Vector(
            map_dsl
                .bonds
                .iter()
                .map(render_covalent_bond_entry_dsl)
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    if !map_dsl.dative_bonds.is_empty() {
        m.insert(
            Edn::keyword("dative"),
            Edn::Vector(
                map_dsl
                    .dative_bonds
                    .iter()
                    .map(render_dative_bond_entry_dsl)
                    .collect::<Vec<_>>()
                    .into(),
            ),
        );
    }
    if !map_dsl.aromatic_systems.is_empty() {
        m.insert(
            Edn::keyword("aromatic"),
            Edn::Vector(
                map_dsl
                    .aromatic_systems
                    .iter()
                    .map(|entry| render_atoms_entry_dsl(&entry.id, &entry.atoms))
                    .collect::<Vec<_>>()
                    .into(),
            ),
        );
    }
    if !map_dsl.multicenter_bonds.is_empty() {
        m.insert(
            Edn::keyword("multicenter"),
            Edn::Vector(
                map_dsl
                    .multicenter_bonds
                    .iter()
                    .map(|entry| render_atoms_entry_dsl(&entry.id, &entry.atoms))
                    .collect::<Vec<_>>()
                    .into(),
            ),
        );
    }
    if !map_dsl.noncovalent_bonds.is_empty() {
        m.insert(
            Edn::keyword("noncovalent"),
            Edn::Vector(
                map_dsl
                    .noncovalent_bonds
                    .iter()
                    .map(render_noncovalent_bond_entry_dsl)
                    .collect::<Vec<_>>()
                    .into(),
            ),
        );
    }
    if !map_dsl.atom_aliases.is_empty() {
        let alias_elems = map_dsl
            .atom_aliases
            .chunks_exact(2)
            .flat_map(|pair| {
                [
                    Edn::Keyword(EdnKeyword::owned(pair[0].clone())),
                    Edn::Str(Cow::Owned(pair[1].clone())),
                ]
            })
            .collect::<Vec<_>>();
        m.insert(Edn::keyword("atom-aliases"), Edn::Vector(alias_elems.into()));
    }
    if !map_dsl.constraints.is_empty() {
        m.insert(
            Edn::keyword("constraints"),
            Edn::Vector(
                map_dsl
                    .constraints
                    .iter()
                    .map(MoleculeConstraintDsl::to_edn)
                    .collect::<Vec<_>>()
                    .into(),
            ),
        );
    }
    m
}

fn render_molecule_edn(ast: &MoleculeAst, metadata: &Metadata) -> Edn<'static> {
    let mut per_atom_derived: Vec<Vec<AtomConstraint>> = vec![Vec::new(); ast.atoms().count()];
    for (idx, set) in ast.constraints().atoms() {
        let i = idx.index();
        if i < per_atom_derived.len() {
            per_atom_derived[i].extend(set.iter().cloned());
        }
    }

    let mut per_bond_derived: Vec<Vec<BondConstraint>> = vec![Vec::new(); ast.bonds().count()];
    for (idx, set) in ast.constraints().bonds() {
        let i = idx.index();
        if i < per_bond_derived.len() {
            per_bond_derived[i].extend(set.iter().cloned());
        }
    }

    let mut atom_elems = Vec::with_capacity(ast.atoms().count());
    for view in ast.atoms().iter() {
        let i = view.idx.index();
        let mut derived = mem::take(&mut per_atom_derived[i]);
        derived.sort_by_key(constraint_sort_key);
        let atom_dsl = AtomDsl::from_parts(view.data.clone(), derived);
        let alias_name = metadata.atom_aliases.get_by_right(&atom_dsl);
        let id = metadata.atom_ids.get(&i);
        let atom_edn = if let Some(alias) = alias_name {
            Edn::Keyword(EdnKeyword::owned(alias.clone()))
        } else {
            atom_dsl.to_edn()
        };
        let entry = if let Some(id_name) = id {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(id_name.clone())), atom_edn].into())
        } else {
            atom_edn
        };
        atom_elems.push(entry);
    }

    let render_endpoint = |idx: usize| -> Edn<'static> {
        if let Some(id) = metadata.atom_ids.get(&idx) {
            Edn::Keyword(EdnKeyword::owned(id.clone()))
        } else {
            Edn::Int(idx as i64)
        }
    };

    let bond_aliases = super::bond::builtin_bond_aliases();
    let bonds_edn: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let derived = mem::take(&mut per_bond_derived[i]);
            let bond_dsl = if bond_aliases.get_by_right(b.data).is_some() {
                BondDsl::from_parts(b.data.clone(), Vec::new())
            } else {
                BondDsl::from_parts(b.data.clone(), derived)
            };
            render_localized(
                b.src.index(),
                b.tgt.index(),
                &bond_dsl,
                i,
                &metadata.bond_ids,
                &render_endpoint,
            )
        })
        .collect();

    let dative_edn: Vec<Edn<'static>> = ast
        .dative_bonds()
        .iter()
        .enumerate()
        .map(|(i, v)| {
            render_dative(
                v.donor.index(),
                v.acceptor.index(),
                v.data,
                i,
                &metadata.dative_bond_ids,
                &render_endpoint,
            )
        })
        .collect();

    let aromatic_edn: Vec<Edn<'static>> = ast
        .aromatic_systems()
        .iter()
        .enumerate()
        .map(|(i, v)| render_atoms_map(v.atoms(), i, &metadata.aromatic_system_ids, &render_endpoint))
        .collect();

    let multicenter_edn: Vec<Edn<'static>> = ast
        .multicenter_bonds()
        .iter()
        .enumerate()
        .map(|(i, v)| {
            render_atoms_map(
                v.atoms(),
                i,
                &metadata.multicenter_bond_ids,
                &render_endpoint,
            )
        })
        .collect();

    let noncovalent_edn: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .enumerate()
        .map(|(i, v)| {
            render_noncovalent(
                v.atoms[0].index(),
                v.atoms[1].index(),
                v.data,
                i,
                &metadata.noncovalent_bond_ids,
                &render_endpoint,
            )
        })
        .collect();

    let has_aliases = !metadata.atom_aliases.is_empty();
    let mut m = EdnMap::with_capacity(10);
    m.insert(Edn::keyword("atoms"), Edn::Vector(atom_elems.into()));
    m.insert(Edn::keyword("bonds"), Edn::Vector(bonds_edn.into()));
    if ast.dative_bonds().count() > 0 {
        m.insert(Edn::keyword("dative"), Edn::Vector(dative_edn.into()));
    }
    if ast.aromatic_systems().count() > 0 {
        m.insert(Edn::keyword("aromatic"), Edn::Vector(aromatic_edn.into()));
    }
    if ast.multicenter_bonds().count() > 0 {
        m.insert(Edn::keyword("multicenter"), Edn::Vector(multicenter_edn.into()));
    }
    if ast.noncovalent_bonds().count() > 0 {
        m.insert(Edn::keyword("noncovalent"), Edn::Vector(noncovalent_edn.into()));
    }
    let constraint_entries = render_constraint_list(ast, metadata);
    if !constraint_entries.is_empty() {
        m.insert(
            Edn::keyword("constraints"),
            Edn::Vector(constraint_entries.into()),
        );
    }
    if has_aliases {
        let mut alias_elems = Vec::with_capacity(metadata.atom_aliases.len() * 2);
        for (name, atom) in metadata.atom_aliases.iter() {
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

impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Map(self.to_edn_map())
    }
}

fn render_localized(
    source: usize,
    target: usize,
    bond: &BondDsl,
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
        MoleculeConstraint::SubPattern {
            target_anchor,
            pattern_anchor,
            pattern,
        } => {
            let wrapper = MoleculeDsl::from_ast((**pattern).clone());
            let (_, wrapper_metadata) = wrapper
                .lower_parts()
                .expect("MoleculeDsl::from_ast should always lower");
            let mut inner = EdnMap::with_capacity(3);
            inner.insert(
                Edn::keyword("target-anchor"),
                render_id_or_int(&metadata.atom_ids, target_anchor.index()),
            );
            inner.insert(
                Edn::keyword("pattern-anchor"),
                render_id_or_int(&wrapper_metadata.atom_ids, pattern_anchor.index()),
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
            if !AtomDsl::packs_constraint(c) {
                out.push(render_atom_pred(*idx, c, &metadata.atom_ids));
            }
        }
    }

    let bond_aliases = super::bond::builtin_bond_aliases();
    for (idx, set) in constraints.bonds() {
        let bond_has_alias = bond_aliases.get_by_right(&ast[*idx]).is_some();
        for c in set.iter() {
            if bond_has_alias || !BondDsl::packs_constraint(c) {
                out.push(render_bond_pred(*idx, c, &metadata.bond_ids));
            }
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
        out.push(render_molecule_constraint(c, metadata));
    }

    out
}

impl fmt::Display for MoleculeDsl {
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
        r#"{:atoms [[:F "F#c-"]] :bonds [] :constraints [{:total-charge -1}]}"#,
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
                "C #h1".parse().unwrap(),
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
                AtomDsl::from_pattern(AtomPattern::new(AtomAst::from_element(e!(N)))),
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
        let (ast, metadata) = dsl.lower_parts().unwrap();
        assert_eq!(ast, expected_ast);
        assert_eq!(metadata, expected_meta);
    }

    #[test]
    fn test_parse_molecule_dsl_dative() {
        let dsl = parse_molecule_dsl(
            r#"{:atoms [[:B "B #h3"] [:N "N #h3"]]
                :bonds []
                :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#,
        )
        .unwrap();
        let (ast, metadata) = dsl.lower_parts().unwrap();
        assert_eq!(ast.dative_bonds().count(), 1);
        let view = ast.dative_bonds().iter().next().unwrap();
        assert_eq!(view.donor, AtomIdx(1)); // donor = N
        assert_eq!(view.acceptor, AtomIdx(0)); // acceptor = B
        assert_eq!(
            metadata.dative_bond_ids.get(&0),
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
        let (ast, metadata) = dsl.lower_parts().unwrap();
        assert_eq!(ast.aromatic_systems().count(), 1);
        assert_eq!(
            metadata.aromatic_system_ids.get(&0),
            Some(&"ar1".to_string())
        );
        let view = ast.aromatic_systems().iter().next().unwrap();
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
        assert_eq!(dsl1.lower_parts().unwrap(), dsl2.lower_parts().unwrap());
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
    #[case::charge_spin(r##"{:atoms ["C"] :bonds [] :constraints [{:total-charge -1} {:total-spin "#u1"}]}"##)]
    #[case::dative(
        r#"{:atoms [[:N "N"] [:B "B"]] :bonds [[:N :B :single]] :dative [{:donor :N :acceptor :B :bond :single}]}"#
    )]
    #[case::aromatic(
        r#"{:atoms [[:C1 "C"] [:C2 "C"]] :bonds [[:C1 :C2 :single]] :aromatic [{:id :ar1 :atoms [:C1 :C2]}]}"#
    )]
    fn test_molecule_dsl_to_edn_roundtrip(#[case] input: &str) {
        let dsl = parse_molecule_dsl(input).unwrap();
        let edn = dsl.to_edn();
        let back = MoleculeDsl::from_edn(&edn).unwrap();
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
        r#"{:atoms [[:a "C"]] :bonds [] :constraints [{:sub-pattern {:target-anchor :a :pattern-anchor 0 :pattern {:atoms ["C"] :bonds []}}}]}"#
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
