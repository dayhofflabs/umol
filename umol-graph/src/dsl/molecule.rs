//! Molecule map DSL: parser and AST

use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;

use indexmap::IndexMap;
use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::{self, SerializeMap, SerializeSeq, Serializer};
use serde::{Deserialize, Serialize};
use umol_data::SpinState;
use umol_edn::{Edn, EdnDeserializer};

use super::ast::DslAst;
use super::atom::{parse_atom_dsl, AtomAst};
use super::bond::{parse_bond_dsl, BondAst};
use super::config::MoleculeDslConfig;
use super::error::ParseError;

/// `:atoms` - either a named map or an indexed vector
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Atoms {
    Named(IndexMap<String, AtomAst>),
    Indexed(Vec<AtomAst>),
}

impl Default for Atoms {
    fn default() -> Self {
        Self::Indexed(vec![])
    }
}

/// `:bond` value on a bond entry: parsed bond-string or keyword shorthand
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondSpec {
    Literal(BondAst),
    Single,
    Double,
    Triple,
    Quadruple,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CovalentBond {
    pub id: Option<String>,
    pub a: String,
    pub b: String,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DativeBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub donor: String,
    pub acceptor: String,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AromaticSystem {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub atoms: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MulticenterBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub atoms: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NoncovalentBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    pub a: String,
    pub b: String,
    pub bond: BondSpec,
}

/// Parsed molecule map AST
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct MoleculeAst {
    pub atoms: Atoms,
    pub bonds: Vec<CovalentBond>,
    pub dative_bonds: Vec<DativeBond>,
    pub aromatic_systems: Vec<AromaticSystem>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub noncovalent_bonds: Vec<NoncovalentBond>,
    pub charge: Option<i64>,
    pub spin: Option<SpinState>,
}

impl DslAst for MoleculeAst {
    type Config = MoleculeDslConfig;
}

// Serde: BondSpec
impl Serialize for BondSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            BondSpec::Single => serializer.serialize_unit_variant("BondSpec", 0, "single"),
            BondSpec::Double => serializer.serialize_unit_variant("BondSpec", 1, "double"),
            BondSpec::Triple => serializer.serialize_unit_variant("BondSpec", 2, "triple"),
            BondSpec::Quadruple => serializer.serialize_unit_variant("BondSpec", 3, "quadruple"),
            BondSpec::Literal(ast) => ast.serialize(serializer),
        }
    }
}

impl<'de> Deserialize<'de> for BondSpec {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "single" => Ok(BondSpec::Single),
            "double" => Ok(BondSpec::Double),
            "triple" => Ok(BondSpec::Triple),
            "quadruple" => Ok(BondSpec::Quadruple),
            other => {
                let ast = parse_bond_dsl(other).map_err(de::Error::custom)?;
                Ok(BondSpec::Literal(ast))
            }
        }
    }
}

// Serde: Atoms
impl Serialize for Atoms {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Atoms::Named(map) => {
                let mut m = serializer.serialize_map(Some(map.len()))?;
                for (label, atom) in map {
                    m.serialize_entry(label, atom)?;
                }
                m.end()
            }
            Atoms::Indexed(vec) => {
                let mut s = serializer.serialize_seq(Some(vec.len()))?;
                for atom in vec {
                    s.serialize_element(atom)?;
                }
                s.end()
            }
        }
    }
}

impl<'de> Deserialize<'de> for Atoms {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(AtomsVisitor)
    }
}

struct AtomsVisitor;

impl<'de> Visitor<'de> for AtomsVisitor {
    type Value = Atoms;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a map or vector of atoms")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Atoms, A::Error> {
        let mut named = IndexMap::new();
        while let Some((label, atom)) = map.next_entry::<String, AtomAst>()? {
            named.insert(label, atom);
        }
        Ok(Atoms::Named(named))
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Atoms, A::Error> {
        let mut atoms = Vec::new();
        while let Some(atom) = seq.next_element::<AtomAst>()? {
            atoms.push(atom);
        }
        Ok(Atoms::Indexed(atoms))
    }

    fn visit_unit<E: de::Error>(self) -> Result<Atoms, E> {
        Ok(Atoms::Indexed(vec![]))
    }
}

// Serde: CovalentBond
impl Serialize for CovalentBond {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if self.id.is_none() {
            // Shorthand: [:a :b :bond-spec]
            let mut t = serializer.serialize_tuple(3)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.a)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.b)?;
            ser::SerializeTuple::serialize_element(&mut t, &self.bond)?;
            ser::SerializeTuple::end(t)
        } else {
            // Map form: {:id :b1 :a :0 :b :1 :bond :single}
            let mut m = serializer.serialize_struct("CovalentBond", 4)?;
            if let Some(id) = &self.id {
                ser::SerializeStruct::serialize_field(&mut m, "id", id)?;
            }
            ser::SerializeStruct::serialize_field(&mut m, "a", &self.a)?;
            ser::SerializeStruct::serialize_field(&mut m, "b", &self.b)?;
            ser::SerializeStruct::serialize_field(&mut m, "bond", &self.bond)?;
            ser::SerializeStruct::end(m)
        }
    }
}

impl<'de> Deserialize<'de> for CovalentBond {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(CovalentBondVisitor)
    }
}

struct CovalentBondVisitor;

impl<'de> Visitor<'de> for CovalentBondVisitor {
    type Value = CovalentBond;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a bond vector [:a :b :spec] or map {:a :A :b :B :bond :spec}")
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<CovalentBond, A::Error> {
        let a: String = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &"3-element vector"))?;
        let b: String = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &"3-element vector"))?;
        let bond: BondSpec = seq
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(2, &"3-element vector"))?;
        Ok(CovalentBond {
            id: None,
            a,
            b,
            bond,
        })
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<CovalentBond, A::Error> {
        let mut id = None;
        let mut a = None;
        let mut b = None;
        let mut bond = None;
        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "id" => id = Some(map.next_value::<String>()?),
                "a" => a = Some(map.next_value::<String>()?),
                "b" => b = Some(map.next_value::<String>()?),
                "bond" => bond = Some(map.next_value::<BondSpec>()?),
                _ => {
                    let _ = map.next_value::<de::IgnoredAny>()?;
                }
            }
        }
        Ok(CovalentBond {
            id,
            a: a.ok_or_else(|| de::Error::missing_field("a"))?,
            b: b.ok_or_else(|| de::Error::missing_field("b"))?,
            bond: bond.ok_or_else(|| de::Error::missing_field("bond"))?,
        })
    }
}

// Serde: MoleculeAst
impl Serialize for MoleculeAst {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        // Count fields: atoms + bonds are always present
        let mut count = 2;
        if !self.dative_bonds.is_empty() {
            count += 1;
        }
        if !self.aromatic_systems.is_empty() {
            count += 1;
        }
        if !self.multicenter_bonds.is_empty() {
            count += 1;
        }
        if !self.noncovalent_bonds.is_empty() {
            count += 1;
        }
        if self.charge.is_some() {
            count += 1;
        }
        if self.spin.is_some() {
            count += 1;
        }

        let mut m = serializer.serialize_struct("MoleculeAst", count)?;
        ser::SerializeStruct::serialize_field(&mut m, "atoms", &self.atoms)?;
        ser::SerializeStruct::serialize_field(&mut m, "bonds", &self.bonds)?;
        if !self.dative_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "dative", &self.dative_bonds)?;
        }
        if !self.aromatic_systems.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "aromatic", &self.aromatic_systems)?;
        }
        if !self.multicenter_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "multicenter", &self.multicenter_bonds)?;
        }
        if !self.noncovalent_bonds.is_empty() {
            ser::SerializeStruct::serialize_field(&mut m, "noncovalent", &self.noncovalent_bonds)?;
        }
        if let Some(charge) = self.charge {
            ser::SerializeStruct::serialize_field(&mut m, "charge", &charge)?;
        }
        if let Some(spin) = &self.spin {
            ser::SerializeStruct::serialize_field(&mut m, "spin", &spin.to_string())?;
        }
        ser::SerializeStruct::end(m)
    }
}

impl<'de> Deserialize<'de> for MoleculeAst {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(MoleculeAstVisitor)
    }
}

struct MoleculeAstVisitor;

impl<'de> Visitor<'de> for MoleculeAstVisitor {
    type Value = MoleculeAst;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a molecule map with :atoms and :bonds")
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<MoleculeAst, A::Error> {
        let mut atoms = None;
        let mut bonds = None;
        let mut dative_bonds = None;
        let mut aromatic_systems = None;
        let mut multicenter_bonds = None;
        let mut noncovalent_bonds = None;
        let mut charge = None;
        let mut spin = None;
        // Aliases are resolved inline during atom deserialization (not stored).
        // We collect them as raw strings here, then resolve after all fields are read.
        let mut aliases: Option<Vec<(String, String)>> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "atoms" => atoms = Some(map.next_value::<Atoms>()?),
                "bonds" => bonds = Some(map.next_value::<Vec<CovalentBond>>()?),
                "dative" => dative_bonds = Some(map.next_value::<Vec<DativeBond>>()?),
                "aromatic" => aromatic_systems = Some(map.next_value::<Vec<AromaticSystem>>()?),
                "multicenter" => {
                    multicenter_bonds = Some(map.next_value::<Vec<MulticenterBond>>()?)
                }
                "noncovalent" => {
                    noncovalent_bonds = Some(map.next_value::<Vec<NoncovalentBond>>()?)
                }
                "charge" => charge = map.next_value::<Option<i64>>()?,
                "spin" => {
                    let s: Option<String> = map.next_value()?;
                    spin = match s {
                        Some(ref s) => Some(
                            SpinState::from_str(s).map_err(|e| de::Error::custom(e.to_string()))?,
                        ),
                        None => None,
                    };
                }
                "aliases" => aliases = Some(map.next_value::<Vec<(String, String)>>()?),
                _ => {
                    let _ = map.next_value::<de::IgnoredAny>()?;
                }
            }
        }

        let atoms = atoms.ok_or_else(|| de::Error::missing_field("atoms"))?;
        let bonds = bonds.ok_or_else(|| de::Error::missing_field("bonds"))?;
        let _ = aliases;

        Ok(MoleculeAst {
            atoms,
            bonds,
            dative_bonds: dative_bonds.unwrap_or_default(),
            aromatic_systems: aromatic_systems.unwrap_or_default(),
            multicenter_bonds: multicenter_bonds.unwrap_or_default(),
            noncovalent_bonds: noncovalent_bonds.unwrap_or_default(),
            charge,
            spin,
        })
    }
}

// Alias resolution (post-deserialization)

/// Pre-validate atom strings at the EDN level so parse errors get proper types.
fn validate_edn_atom_strings(top: &Edn<'_>) -> Result<(), ParseError> {
    let Edn::Map(map) = top else { return Ok(()) };
    let Some(atoms) = map.get(&Edn::keyword("atoms")) else {
        return Ok(());
    };
    let values: Box<dyn Iterator<Item = &Edn<'_>>> = match atoms {
        Edn::Map(m) => Box::new(m.values()),
        Edn::Vector(v) => Box::new(v.iter()),
        _ => return Ok(()),
    };
    for v in values {
        if let Edn::Str(s) = v {
            parse_atom_dsl(s)?;
        }
    }
    Ok(())
}

fn validate(ast: &MoleculeAst) -> Result<(), ParseError> {
    let atom_labels: HashSet<String> = match &ast.atoms {
        Atoms::Named(m) => m.keys().cloned().collect(),
        Atoms::Indexed(v) => (0..v.len()).map(|i| i.to_string()).collect(),
    };

    let mut seen_ids: HashSet<&str> = HashSet::new();
    if let Atoms::Named(m) = &ast.atoms {
        for label in m.keys() {
            seen_ids.insert(label.as_str());
        }
    }
    for id in ast
        .bonds
        .iter()
        .filter_map(|e| e.id.as_deref())
        .chain(ast.dative_bonds.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.aromatic_systems.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.multicenter_bonds.iter().filter_map(|e| e.id.as_deref()))
        .chain(ast.noncovalent_bonds.iter().filter_map(|e| e.id.as_deref()))
    {
        if !seen_ids.insert(id) {
            return Err(ParseError::DuplicateId(id.to_string()));
        }
    }

    let check = |label: &str| -> Result<(), ParseError> {
        if atom_labels.contains(label) {
            Ok(())
        } else {
            Err(ParseError::InvalidAtomIndex(label.to_string()))
        }
    };

    for e in &ast.bonds {
        check(&e.a)?;
        check(&e.b)?;
    }
    for e in &ast.dative_bonds {
        check(&e.donor)?;
        check(&e.acceptor)?;
    }
    for e in &ast.aromatic_systems {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.multicenter_bonds {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.noncovalent_bonds {
        check(&e.a)?;
        check(&e.b)?;
    }

    Ok(())
}

// FromStr / Display / parse_molecule_dsl

/// Parse a molecule AST from an EDN string.
///
/// Alias resolution happens at the EDN level before serde deserialization:
/// keyword references (`:alias`) in the atoms section are replaced with their
/// atom string definitions from the `:aliases` vector.
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeAst, ParseError> {
    let top = umol_edn::read_string(input)?;

    // Pre-check: must be a map
    if !matches!(&top, Edn::Map(_)) {
        return Err(ParseError::EdnParse(
            "expected EDN map for top level".to_string(),
        ));
    }

    let top = resolve_edn_aliases(top)?;

    // Pre-check: required keys
    if let Edn::Map(ref m) = top {
        if !m.contains_key(&Edn::keyword("atoms")) {
            return Err(ParseError::MissingKey(":atoms".to_string()));
        }
        if !m.contains_key(&Edn::keyword("bonds")) {
            return Err(ParseError::MissingKey(":bonds".to_string()));
        }
    }

    // Pre-validate atom strings so parse errors surface with proper types
    validate_edn_atom_strings(&top)?;

    let ast = MoleculeAst::deserialize(EdnDeserializer(top)).map_err(ParseError::from)?;

    validate(&ast)?;
    Ok(ast)
}

/// Resolve `:aliases` at the EDN level. Replaces keyword references in `:atoms`
/// with their atom string definitions, then removes the `:aliases` key.
fn resolve_edn_aliases(top: Edn<'_>) -> Result<Edn<'_>, ParseError> {
    let Edn::Map(mut map) = top else {
        return Ok(top);
    };

    let aliases_key = Edn::keyword("aliases");
    let alias_map = if let Some(aliases_edn) = map.remove(&aliases_key) {
        parse_edn_aliases(&aliases_edn)?
    } else {
        IndexMap::new()
    };

    let atoms_key = Edn::keyword("atoms");
    if let Some(atoms) = map.remove(&atoms_key) {
        let resolved = resolve_atoms_aliases(atoms, &alias_map)?;
        map.insert(atoms_key, resolved);
    }

    Ok(Edn::Map(map))
}

fn parse_edn_aliases<'e>(edn: &Edn<'e>) -> Result<IndexMap<String, Edn<'e>>, ParseError> {
    let Edn::Vector(v) = edn else {
        return Err(ParseError::WrongFieldType {
            field: "aliases".to_string(),
            expected: "flat vector of keyword/atom-spec pairs".to_string(),
        });
    };
    if v.len() % 2 != 0 {
        return Err(ParseError::WrongFieldType {
            field: "aliases".to_string(),
            expected: "flat vector of keyword/atom-spec pairs (even length)".to_string(),
        });
    }
    let mut aliases = IndexMap::new();
    for pair in v.chunks(2) {
        let Edn::Keyword(k) = &pair[0] else {
            return Err(ParseError::EdnParse(
                "expected keyword as alias name".to_string(),
            ));
        };
        let name = k.as_str().to_string();
        if aliases.contains_key(&name) {
            return Err(ParseError::DuplicateId(name));
        }
        aliases.insert(name, pair[1].clone());
    }
    Ok(aliases)
}

fn resolve_atoms_aliases<'e>(
    atoms: Edn<'e>,
    aliases: &IndexMap<String, Edn<'e>>,
) -> Result<Edn<'e>, ParseError> {
    match atoms {
        Edn::Map(m) => {
            let mut resolved = std::collections::BTreeMap::new();
            for (k, v) in m {
                let v = resolve_one_atom(v, aliases)?;
                if let Edn::Keyword(ref kw) = k {
                    if aliases.contains_key(kw.as_str()) {
                        return Err(ParseError::DuplicateId(kw.as_str().to_string()));
                    }
                }
                resolved.insert(k, v);
            }
            Ok(Edn::Map(resolved))
        }
        Edn::Vector(items) => {
            let resolved = items
                .into_iter()
                .map(|e| resolve_one_atom(e, aliases))
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Edn::Vector(resolved))
        }
        other => Ok(other),
    }
}

fn resolve_one_atom<'e>(
    edn: Edn<'e>,
    aliases: &IndexMap<String, Edn<'e>>,
) -> Result<Edn<'e>, ParseError> {
    match &edn {
        Edn::Keyword(k) => aliases
            .get(k.as_str())
            .cloned()
            .ok_or_else(|| ParseError::UnknownAlias(k.as_str().to_string())),
        _ => Ok(edn),
    }
}

impl FromStr for MoleculeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_molecule_dsl(s)
    }
}

impl fmt::Display for MoleculeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = umol_edn::to_string(self).map_err(|_| fmt::Error)?;
        f.write_str(&s)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{e, spin, Element};

    use super::super::predicates::{ElementExpr, HydrogenExpr};
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    #[case::empty(r#"{:atoms [] :bonds []}"#, MoleculeAst::default())]
    #[case::atom(r#"{:atoms ["C"] :bonds []}"#, MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(C))]), ..Default::default() })]
    #[case::atom_id(r#"{:atoms {:C "C"} :bonds []}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C)))])), ..Default::default() })]
    #[case::atom_dsl(r#"{:atoms {:C "C #h4"} :bonds []}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None, implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(4))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), ..Default::default() })]
    #[case::bond(r#"{:atoms ["N" "N"] :bonds [[:0 :1 :triple]]}"#, MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![CovalentBond { id: None, a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Triple }], ..Default::default() })]
    #[case::bond_atom_ids(r#"{:atoms {:C "C" :O "O"} :bonds [[:C :O :single]]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: None, a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_id(r#"{:atoms ["H" "F"] :bonds [{:id :b1 :a :0 :b :1 :bond :single}]}"#,
        MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(H)), AtomAst::from_element(e!(F))]),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_id_atom_ids(r#"{:atoms {:C "C" :O "O"} :bonds [{:id :b1 :a :C :b :O :bond :single}]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_dsl(r#"{:atoms {:C "C" :O "O"} :bonds [{:id :b1 :a :C :b :O :bond "2"}]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Literal(BondAst { order: ValueAst::Lit(2), charge: None,
        unpaired_electrons: None, multiplicity: None }) }], ..Default::default() })]
    #[case::charge(r#"{:atoms {:F "F#c-"} :bonds [] :charge -1}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("F".to_string(), 
        AtomAst { element: ElementExpr::Lit(Element::F), isotope_mass: None, implicit_hydrogens: None, charge: Some(ValueAst::Lit(-1)), lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), charge: Some(-1), ..Default::default() })]
    #[case::spin(r##"{:atoms {:N "N #u3"} :bonds [] :spin "#u3"}"##, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("N".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::N), isotope_mass: None, implicit_hydrogens: None, charge: None, lone_pairs: None, unpaired_electrons: Some(ValueAst::Lit(3)),
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), spin: Some(spin!("#u3 #s4")), ..Default::default() })]
    #[case::alias_named(r#"{:atoms {:C :ch} :bonds [] :aliases [:ch "C #h1"]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None, implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None,
        unpaired_electrons: None, multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None })])), ..Default::default() })]
    #[case::alias_indexed(r#"{:atoms [:ch] :bonds [] :aliases [:ch "C #h1"]}"#,
        MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None }]), ..Default::default() })]
    #[case::alias_reused(r#"{:atoms [:n :n] :bonds [[:0 :1 :single]] :aliases [:n "N"]}"#,
        MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![CovalentBond { id: None, a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    fn test_parse_molecule_dsl(#[case] input: &str, #[case] expected: MoleculeAst) {
        let result = parse_molecule_dsl(input);
        assert!(
            result.is_ok(),
            "{input:?} should succeed, got {:?}",
            result.unwrap_err()
        );
        let ast = result.unwrap();
        assert_eq!(ast, expected);
    }

    #[test]
    fn test_parse_molecule_dsl_dative() {
        let result = parse_molecule_dsl(
            r#"{:atoms {:B "B #h3" :N "N #h3"}
                :bonds []
                :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#,
        )
        .unwrap();
        assert_eq!(
            result.dative_bonds,
            vec![DativeBond {
                id: Some("d1".to_string()),
                donor: "N".to_string(),
                acceptor: "B".to_string(),
                bond: BondSpec::Single,
            }]
        );
    }

    #[test]
    fn test_parse_molecule_dsl_aromatic() {
        let result = parse_molecule_dsl(
            r#"{:atoms {:C1 :ch :C2 :ch :C3 :ch :C4 :ch :C5 :ch :C6 :ch}
                :bonds [[:C1 :C2 :single] [:C2 :C3 :single] [:C3 :C4 :single] [:C4 :C5 :single] [:C5 :C6 :single] [:C6 :C1 :single]]
                :aromatic [{:id :ar1 :atoms [:C1 :C2 :C3 :C4 :C5 :C6]}]
                :aliases [:ch "C #h1 #v2 #a1"]}"#,
        )
        .unwrap();
        assert_eq!(result.aromatic_systems.len(), 1);
        assert_eq!(result.aromatic_systems[0].id, Some("ar1".to_string()));
        assert_eq!(
            result.aromatic_systems[0].atoms,
            vec!["C1", "C2", "C3", "C4", "C5", "C6"]
        );
    }

    #[rstest]
    #[case::empty(MoleculeAst::default())]
    #[case::named_atoms(MoleculeAst {
        atoms: Atoms::Named(IndexMap::from([
            ("C".to_string(), AtomAst::from_element(e!(C))),
            ("O".to_string(), AtomAst::from_element(e!(O))),
        ])),
        bonds: vec![CovalentBond { id: None, a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Single }],
        ..Default::default()
    })]
    #[case::indexed_atoms(MoleculeAst {
        atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![CovalentBond { id: None, a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Triple }],
        ..Default::default()
    })]
    #[case::bond_with_id(MoleculeAst {
        atoms: Atoms::Named(IndexMap::from([
            ("C".to_string(), AtomAst::from_element(e!(C))),
            ("O".to_string(), AtomAst::from_element(e!(O))),
        ])),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Double }],
        ..Default::default()
    })]
    #[case::dative(MoleculeAst {
        atoms: Atoms::Named(IndexMap::from([
            ("B".to_string(), AtomAst::from_element(e!(B))),
            ("N".to_string(), AtomAst::from_element(e!(N))),
        ])),
        bonds: vec![],
        dative_bonds: vec![DativeBond { id: Some("d1".to_string()), donor: "N".to_string(), acceptor: "B".to_string(), bond: BondSpec::Single }],
        ..Default::default()
    })]
    fn test_molecule_ast_json_roundtrip(#[case] ast: MoleculeAst) {
        let json = serde_json::to_string(&ast).expect("serialize to JSON");
        let back: MoleculeAst =
            serde_json::from_str(&json).expect(&format!("deserialize from JSON: {json}"));
        assert_eq!(ast, back);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::non_map("3", ParseError::EdnParse("expected EDN map for top level".to_string()))]
    #[case::missing_atoms(r#"{:bonds []}"#, ParseError::MissingKey(":atoms".to_string()))]
    #[case::missing_bonds(r#"{:atoms {:C "C"}}"#, ParseError::MissingKey(":bonds".to_string()))]
    #[case::unknown_endpoint(r#"{:atoms {:C "C"} :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#, ParseError::InvalidAtomIndex("X".to_string()))]
    #[case::duplicate_id(r#"{:atoms {:C "C" :O "O" :N "N"} :bonds [{:id :b1 :a :C :b :O :bond :single} {:id :b1 :a :O :b :N :bond :single}]}"#,
        ParseError::DuplicateId("b1".to_string()))]
    #[case::bad_atom_string(r##"{:atoms {:X "#h3"} :bonds []}"##, ParseError::InvalidElement("#h3".to_string()))]
    #[case::unknown_alias(r#"{:atoms {:C :ch} :bonds []}"#, ParseError::UnknownAlias("ch".to_string()))]
    #[case::trailing_content(r#"{:atoms {:C "C"} :bonds []} :extra :junk"#, ParseError::EdnParse("trailing content at byte 28".to_string()))]
    #[case::duplicate_atom_bond_id(r#"{:atoms {:b1 "C" :O "O"} :bonds [{:id :b1 :a :b1 :b :O :bond :single}]}"#, ParseError::DuplicateId("b1".to_string()))]
    #[case::duplicate_bond_dative_id(r#"{:atoms {:C "C" :O "O"} :bonds [{:id :b1 :a :C :b :O :bond :single}] :dative [{:id :b1 :donor :C :acceptor :O :bond :single}]}"#, ParseError::DuplicateId("b1".to_string()))]
    #[case::duplicate_id_alias(r#"{:aliases [:C "N"] :atoms {:C "C"} :bonds []}"#, ParseError::DuplicateId("C".to_string()))]
    #[case::duplicate_alias(r#"{:aliases [:ch "C #h1" :ch "C #h2"] :atoms [] :bonds []}"#, ParseError::DuplicateId("ch".to_string()))]
    fn test_parse_molecule_map_invalid(#[case] input: &str, #[case] expected: ParseError) {
        let result = parse_molecule_dsl(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        assert_eq!(result.unwrap_err(), expected, "for input {input:?}");
    }
}
