//! Molecule map DSL parser

use std::collections::HashSet;
use std::str::FromStr;

use clojure_reader::edn::{self, Edn};
use indexmap::IndexMap;
use umol_data::SpinState;

use super::atom::{parse_atom_dsl, AtomAst};
use super::bond::{parse_bond_dsl, BondAst};
use super::error::ParseError;

mod utils;
use utils::{extract_label, extract_list, extract_map, extract_tagged_str, map_get};

/// `:atoms` — either a named map or an indexed vector
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DativeBond {
    pub id: Option<String>,
    pub donor: String,
    pub acceptor: String,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    pub id: Option<String>,
    pub atoms: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MulticenterBond {
    pub id: Option<String>,
    pub atoms: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NoncovalentBond {
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

/// Parse a molecule AST from a EDN string
pub fn parse_molecule_dsl(input: &str) -> Result<MoleculeAst, ParseError> {
    let (top, rest) = edn::read(input).map_err(|e| ParseError::EdnParse(e.to_string()))?;
    let rest = rest.trim();
    if !rest.is_empty() {
        return Err(ParseError::EdnParse(format!("unexpected trailing content: {rest}")));
    }
    let map = extract_map(&top, "top level")?;

    let aliases = match map_get(map, "aliases") {
        Some(edn) => parse_aliases(edn)?,
        None => IndexMap::new(),
    };
    let atoms = parse_atoms(
        map_get(map, "atoms").ok_or_else(|| ParseError::MissingKey(":atoms".to_string()))?,
        &aliases,
    )?;
    if let Atoms::Named(m) = &atoms {
        for label in m.keys() {
            if aliases.contains_key(label) {
                return Err(ParseError::DuplicateId(label.clone()));
            }
        }
    }
    let bonds = parse_covalent_bonds(
        map_get(map, "bonds").ok_or_else(|| ParseError::MissingKey(":bonds".to_string()))?,
    )?;
    let dative = extract_list(map, "dative", parse_dative_bond)?;
    let aromatic = extract_list(map, "aromatic", parse_aromatic_system)?;
    let multicenter = extract_list(map, "multicenter", parse_multicenter_bond)?;
    let noncovalent = extract_list(map, "noncovalent", parse_noncovalent_bond)?;
    let charge = match map_get(map, "charge") {
        Some(Edn::Int(n)) => Some(*n),
        Some(Edn::Nil) | None => None,
        Some(_) => {
            return Err(ParseError::WrongFieldType {
                field: "charge".to_string(),
                expected: "integer or nil".to_string(),
            })
        }
    };
    let spin = match map_get(map, "spin") {
        Some(edn @ Edn::Str(_)) => Some(parse_spin_state(edn)?),
        Some(Edn::Nil) | None => None,
        Some(_) => {
            return Err(ParseError::WrongFieldType {
                field: "spin".to_string(),
                expected: "string or nil".to_string(),
            })
        }
    };

    let ast = MoleculeAst {
        atoms,
        bonds,
        dative_bonds: dative,
        aromatic_systems: aromatic,
        multicenter_bonds: multicenter,
        noncovalent_bonds: noncovalent,
        charge,
        spin,
    };
    validate(&ast)?;
    Ok(ast)
}

fn parse_aliases(edn: &Edn<'_>) -> Result<IndexMap<String, AtomAst>, ParseError> {
    let v = match edn {
        Edn::Vector(v) => v,
        _ => {
            return Err(ParseError::WrongFieldType {
                field: "aliases".to_string(),
                expected: "flat vector of keyword/atom-spec pairs".to_string(),
            })
        }
    };
    if v.len() % 2 != 0 {
        return Err(ParseError::WrongFieldType {
            field: "aliases".to_string(),
            expected: "flat vector of keyword/atom-spec pairs (even length)".to_string(),
        });
    }
    let mut aliases = IndexMap::new();
    for pair in v.chunks(2) {
        let name = extract_label(&pair[0])?;
        if aliases.contains_key(&name) {
            return Err(ParseError::DuplicateId(name));
        }
        let atom = parse_atom_dsl(extract_tagged_str(&pair[1], "atom")?)?;
        aliases.insert(name, atom);
    }
    Ok(aliases)
}

fn resolve_atom<'e>(
    edn: &'e Edn<'e>,
    aliases: &IndexMap<String, AtomAst>,
) -> Result<AtomAst, ParseError> {
    match edn {
        Edn::Key(name) => aliases
            .get(*name)
            .cloned()
            .ok_or_else(|| ParseError::UnknownAlias((*name).to_string())),
        _ => Ok(parse_atom_dsl(extract_tagged_str(edn, "atom")?)?),
    }
}

fn parse_atoms(edn: &Edn<'_>, aliases: &IndexMap<String, AtomAst>) -> Result<Atoms, ParseError> {
    match edn {
        Edn::Map(m) => {
            let mut named = IndexMap::new();
            for (k, v) in m {
                let label = extract_label(k)?;
                let atom = resolve_atom(v, aliases)?;
                named.insert(label, atom);
            }
            Ok(Atoms::Named(named))
        }
        Edn::Vector(items) => {
            let atoms = items
                .iter()
                .map(|e| resolve_atom(e, aliases))
                .collect::<Result<_, _>>()?;
            Ok(Atoms::Indexed(atoms))
        }
        _ => Err(ParseError::WrongFieldType {
            field: "atoms".to_string(),
            expected: "map or vector".to_string(),
        }),
    }
}

fn parse_bond_spec(edn: &Edn<'_>) -> Result<BondSpec, ParseError> {
    match edn {
        Edn::Tagged(t, v) if *t == "bond" => match v.as_ref() {
            Edn::Str(s) => Ok(BondSpec::Literal(parse_bond_dsl(s)?)),
            _ => Err(ParseError::InvalidBondDsl(
                "#bond must be followed by a string".to_string(),
            )),
        },
        Edn::Key("single") => Ok(BondSpec::Single),
        Edn::Key("double") => Ok(BondSpec::Double),
        Edn::Key("triple") => Ok(BondSpec::Triple),
        Edn::Key("quadruple") => Ok(BondSpec::Quadruple),
        _ => Err(ParseError::InvalidBondDsl(
            "bond spec must be #bond \"...\" or a keyword shorthand".to_string(),
        )),
    }
}

fn parse_covalent_bonds(edn: &Edn<'_>) -> Result<Vec<CovalentBond>, ParseError> {
    match edn {
        Edn::Vector(v) => v.iter().map(parse_covalent_bond).collect(),
        _ => Err(ParseError::WrongFieldType {
            field: "bonds".to_string(),
            expected: "vector".to_string(),
        }),
    }
}

fn parse_covalent_bond(edn: &Edn<'_>) -> Result<CovalentBond, ParseError> {
    match edn {
        Edn::Vector(v) => {
            if v.len() != 3 {
                return Err(ParseError::InvalidBond);
            }
            Ok(CovalentBond {
                id: None,
                a: extract_label(&v[0])?,
                b: extract_label(&v[1])?,
                bond: parse_bond_spec(&v[2])?,
            })
        }
        Edn::Map(_) => {
            let m = extract_map(edn, "covalent bond entry")?;
            Ok(CovalentBond {
                id: map_get(m, "id").map(extract_label).transpose()?,
                a: extract_label(
                    map_get(m, "a")
                        .ok_or_else(|| ParseError::MissingKey(":a in bond entry".to_string()))?,
                )?,
                b: extract_label(
                    map_get(m, "b")
                        .ok_or_else(|| ParseError::MissingKey(":b in bond entry".to_string()))?,
                )?,
                bond: parse_bond_spec(
                    map_get(m, "bond")
                        .ok_or_else(|| ParseError::MissingKey(":bond in bond entry".to_string()))?,
                )?,
            })
        }
        _ => Err(ParseError::InvalidBond),
    }
}

fn parse_dative_bond(edn: &Edn<'_>) -> Result<DativeBond, ParseError> {
    let m = extract_map(edn, "dative bond entry")?;
    Ok(DativeBond {
        id: map_get(m, "id").map(extract_label).transpose()?,
        donor: extract_label(
            map_get(m, "donor")
                .ok_or_else(|| ParseError::MissingKey(":donor in dative entry".to_string()))?,
        )?,
        acceptor: extract_label(
            map_get(m, "acceptor")
                .ok_or_else(|| ParseError::MissingKey(":acceptor in dative entry".to_string()))?,
        )?,
        bond: parse_bond_spec(
            map_get(m, "bond")
                .ok_or_else(|| ParseError::MissingKey(":bond in dative entry".to_string()))?,
        )?,
    })
}

fn parse_aromatic_system(edn: &Edn<'_>) -> Result<AromaticSystem, ParseError> {
    let m = extract_map(edn, "aromatic entry")?;
    let id = map_get(m, "id").map(extract_label).transpose()?;
    let atoms = match map_get(m, "atoms")
        .ok_or_else(|| ParseError::MissingKey(":atoms in aromatic entry".to_string()))?
    {
        Edn::Vector(v) => v.iter().map(extract_label).collect::<Result<_, _>>()?,
        _ => {
            return Err(ParseError::WrongFieldType {
                field: "atoms".to_string(),
                expected: "vector of keywords".to_string(),
            })
        }
    };
    Ok(AromaticSystem { id, atoms })
}

fn parse_multicenter_bond(edn: &Edn<'_>) -> Result<MulticenterBond, ParseError> {
    let m = extract_map(edn, "multicenter entry")?;
    let id = map_get(m, "id").map(extract_label).transpose()?;
    let atoms = match map_get(m, "atoms")
        .ok_or_else(|| ParseError::MissingKey(":atoms in multicenter entry".to_string()))?
    {
        Edn::Vector(v) => v.iter().map(extract_label).collect::<Result<_, _>>()?,
        _ => {
            return Err(ParseError::WrongFieldType {
                field: "atoms".to_string(),
                expected: "vector of keywords".to_string(),
            })
        }
    };
    Ok(MulticenterBond { id, atoms })
}

fn parse_noncovalent_bond(edn: &Edn<'_>) -> Result<NoncovalentBond, ParseError> {
    let m = extract_map(edn, "noncovalent bond entry")?;
    Ok(NoncovalentBond {
        id: map_get(m, "id").map(extract_label).transpose()?,
        a: extract_label(
            map_get(m, "a")
                .ok_or_else(|| ParseError::MissingKey(":a in noncovalent entry".to_string()))?,
        )?,
        b: extract_label(
            map_get(m, "b")
                .ok_or_else(|| ParseError::MissingKey(":b in noncovalent entry".to_string()))?,
        )?,
        bond: parse_bond_spec(
            map_get(m, "bond")
                .ok_or_else(|| ParseError::MissingKey(":bond in noncovalent entry".to_string()))?,
        )?,
    })
}

fn parse_spin_state(edn: &Edn<'_>) -> Result<SpinState, ParseError> {
    match edn {
        Edn::Str(s) => Ok(SpinState::from_str(s)?),
        _ => Err(ParseError::WrongFieldType {
            field: "spin".to_string(),
            expected: "string".to_string(),
        }),
    }
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

    // TODO: Add charge and spin validation

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{e, spin, Element};

    use super::super::atom::{ElementExpr, HydrogenExpr};
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    #[case::empty(r#"{:atoms [] :bonds []}"#, MoleculeAst::default())]
    #[case::atom(r#"{:atoms [#atom "C"] :bonds []}"#, MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(C))]), ..Default::default() })]
    #[case::atom_id(r#"{:atoms {:C #atom "C"} :bonds []}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C)))])), ..Default::default() })]
    #[case::atom_dsl(r#"{:atoms {:C #atom "C #h4"} :bonds []}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None, implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(4))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), ..Default::default() })]
    #[case::bond(r#"{:atoms [#atom "N" #atom "N"] :bonds [[:0 :1 :triple]]}"#, MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(N)), AtomAst::from_element(e!(N))]),
        bonds: vec![CovalentBond { id: None, a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Triple }], ..Default::default() })]
    #[case::bond_atom_ids(r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [[:C :O :single]]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: None, a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_id(r#"{:atoms [#atom "H" #atom "F"] :bonds [{:id :b1 :a :0 :b :1 :bond :single}]}"#,
        MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst::from_element(e!(H)), AtomAst::from_element(e!(F))]),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "0".to_string(), b: "1".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_id_atom_ids(r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :C :b :O :bond :single}]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Single }], ..Default::default() })]
    #[case::bond_dsl(r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :C :b :O :bond #bond "2"}]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(), AtomAst::from_element(e!(C))), ("O".to_string(), AtomAst::from_element(e!(O)))])),
        bonds: vec![CovalentBond { id: Some("b1".to_string()), a: "C".to_string(), b: "O".to_string(), bond: BondSpec::Literal(BondAst { order: ValueAst::Lit(2), charge: None,
        unpaired_electrons: None, multiplicity: None }) }], ..Default::default() })]
    #[case::charge(r#"{:atoms {:F #atom "F#c-"} :bonds [] :charge -1}"#, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("F".to_string(), 
        AtomAst { element: ElementExpr::Lit(Element::F), isotope_mass: None, implicit_hydrogens: None, charge: Some(ValueAst::Lit(-1)), lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), charge: Some(-1), ..Default::default() })]
    #[case::spin(r##"{:atoms {:N #atom "N #u3"} :bonds [] :spin "#u3"}"##, MoleculeAst { atoms: Atoms::Named(IndexMap::from([("N".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::N), isotope_mass: None, implicit_hydrogens: None, charge: None, lone_pairs: None, unpaired_electrons: Some(ValueAst::Lit(3)),
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None, })])), spin: Some(spin!("#u3 #s4")), ..Default::default() })]
    #[case::alias_named(r#"{:atoms {:C :ch} :bonds [] :aliases [:ch #atom "C #h1"]}"#,
        MoleculeAst { atoms: Atoms::Named(IndexMap::from([("C".to_string(),
        AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None, implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None,
        unpaired_electrons: None, multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None })])), ..Default::default() })]
    #[case::alias_indexed(r#"{:atoms [:ch] :bonds [] :aliases [:ch #atom "C #h1"]}"#,
        MoleculeAst { atoms: Atoms::Indexed(vec![AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: None,
        implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), charge: None, lone_pairs: None, unpaired_electrons: None,
        multiplicity: None, valence: None, donated_pairs: None, accepted_pairs: None, aromatic_valence: None, multicenter_valence: None }]), ..Default::default() })]
    #[case::alias_reused(r#"{:atoms [:n :n] :bonds [[:0 :1 :single]] :aliases [:n #atom "N"]}"#,
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
            r#"{:atoms {:B #atom "B #h3" :N #atom "N #h3"}
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
                :aliases [:ch #atom "C #h1 #v2 #a1"]}"#,
        )
        .unwrap();
        assert_eq!(result.aromatic_systems.len(), 1);
        assert_eq!(result.aromatic_systems[0].id, Some("ar1".to_string()));
        assert_eq!(
            result.aromatic_systems[0].atoms,
            vec!["C1", "C2", "C3", "C4", "C5", "C6"]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::non_map("3", ParseError::EdnParse("expected EDN map for top level".to_string()))]
    #[case::missing_atoms(r#"{:bonds []}"#, ParseError::MissingKey(":atoms".to_string()))]
    #[case::missing_bonds(r#"{:atoms {:C #atom "C"}}"#, ParseError::MissingKey(":bonds".to_string()))]
    #[case::unknown_endpoint(r#"{:atoms {:C #atom "C"} :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#, ParseError::InvalidAtomIndex("X".to_string()))]
    #[case::duplicate_id(r#"{:atoms {:C #atom "C" :O #atom "O" :N #atom "N"} :bonds [{:id :b1 :a :C :b :O :bond :single} {:id :b1 :a :O :b :N :bond :single}]}"#,
        ParseError::DuplicateId("b1".to_string()))]
    #[case::bad_atom_string(r##"{:atoms {:X #atom "#h3"} :bonds []}"##, ParseError::InvalidElement("#h3".to_string()))]
    #[case::unknown_alias(r#"{:atoms {:C :ch} :bonds []}"#, ParseError::UnknownAlias("ch".to_string()))]
    #[case::trailing_content(r#"{:atoms {:C #atom "C"} :bonds []} :extra :junk"#, ParseError::EdnParse("unexpected trailing content: :extra :junk".to_string()))]
    #[case::duplicate_atom_bond_id(r#"{:atoms {:b1 #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :b1 :b :O :bond :single}]}"#, ParseError::DuplicateId("b1".to_string()))]
    #[case::duplicate_bond_ids_cross_section(r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :C :b :O :bond :single}] :dative [{:id :b1 :donor :C :acceptor :O :bond :single}]}"#, ParseError::DuplicateId("b1".to_string()))]
    #[case::duplicate_atom_id_and_alias(r#"{:aliases [:C #atom "N"] :atoms {:C #atom "C"} :bonds []}"#, ParseError::DuplicateId("C".to_string()))]
    #[case::duplicate_alias(r#"{:aliases [:ch #atom "C #h1" :ch #atom "C #h2"] :atoms [] :bonds []}"#, ParseError::DuplicateId("ch".to_string()))]
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
