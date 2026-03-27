//! Molecule map DSL parser — `spec/umol-dsl-spec.md` §2–§4 and §7.7.

use std::collections::{BTreeMap, HashSet};

use clojure_reader::edn::{self, Edn};
use indexmap::IndexMap;

use super::atom::{parse_atom_dsl, AtomAst};
use super::bond::{parse_bond_dsl, BondAst};
use super::error::ParseError;

/// Keyword name label (EDN keyword name part, e.g. `"C"` from `:C`).
pub type Label = String;

/// `:atoms` — either a named map or an indexed vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomCollection {
    Named(IndexMap<Label, AtomAst>),
    Indexed(Vec<AtomAst>),
}

/// `:bond` value on a bond entry: parsed bond-string or keyword shorthand.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondSpec {
    Literal(BondAst),
    Single,
    Double,
    Triple,
    Quadruple,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CovalentBondEntry {
    pub id: Label,
    pub a: Label,
    pub b: Label,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DativeBondEntry {
    pub id: Label,
    pub donor: Label,
    pub acceptor: Label,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticEntry {
    pub id: Label,
    pub atoms: Vec<Label>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MulticenterEntry {
    pub id: Label,
    pub atoms: Vec<Label>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NoncovalentEntry {
    pub id: Label,
    pub a: Label,
    pub b: Label,
    pub bond: BondSpec,
}

/// Parsed molecule map AST.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeMapAst {
    pub atoms: AtomCollection,
    pub bonds: Vec<CovalentBondEntry>,
    pub dative: Vec<DativeBondEntry>,
    pub aromatic: Vec<AromaticEntry>,
    pub multicenter: Vec<MulticenterEntry>,
    pub noncovalent: Vec<NoncovalentEntry>,
    pub charge: Option<i64>,
}

/// Parse a molecule map from an EDN string.
pub fn parse_molecule_map(input: &str) -> Result<MoleculeMapAst, ParseError> {
    let top = edn::read_string(input)
        .map_err(|e| ParseError::EdnParse(e.to_string()))?;
    let map = require_map(&top, "top level")?;

    let atoms = parse_atom_collection(
        map_get(map, "atoms").ok_or_else(|| ParseError::MissingKey(":atoms".into()))?,
    )?;
    let bonds = parse_covalent_bonds(
        map_get(map, "bonds").ok_or_else(|| ParseError::MissingKey(":bonds".into()))?,
    )?;
    let dative = opt_list(map, "dative", parse_dative_entry)?;
    let aromatic = opt_list(map, "aromatic", parse_aromatic_entry)?;
    let multicenter = opt_list(map, "multicenter", parse_multicenter_entry)?;
    let noncovalent = opt_list(map, "noncovalent", parse_noncovalent_entry)?;
    let charge = match map_get(map, "charge") {
        Some(Edn::Int(n)) => Some(*n),
        Some(Edn::Nil) | None => None,
        Some(_) => {
            return Err(ParseError::InvalidMoleculeMap(
                ":charge must be an integer or nil".into(),
            ))
        }
    };

    let ast = MoleculeMapAst { atoms, bonds, dative, aromatic, multicenter, noncovalent, charge };
    validate(&ast)?;
    Ok(ast)
}

// ── low-level helpers ─────────────────────────────────────────────────────────

fn map_get<'e>(map: &'e BTreeMap<Edn<'e>, Edn<'e>>, key: &str) -> Option<&'e Edn<'e>> {
    map.iter()
        .find(|(k, _)| matches!(k, Edn::Key(s) if *s == key))
        .map(|(_, v)| v)
}

fn require_map<'e>(edn: &'e Edn<'e>, ctx: &str) -> Result<&'e BTreeMap<Edn<'e>, Edn<'e>>, ParseError> {
    match edn {
        Edn::Map(m) => Ok(m),
        _ => Err(ParseError::InvalidMoleculeMap(format!("expected EDN map for {ctx}"))),
    }
}

fn require_label(edn: &Edn<'_>) -> Result<Label, ParseError> {
    match edn {
        Edn::Key(s) => Ok((*s).to_string()),
        _ => Err(ParseError::InvalidMoleculeMap("expected EDN keyword as label".into())),
    }
}

fn require_tagged_str<'e>(edn: &'e Edn<'e>, tag: &str) -> Result<&'e str, ParseError> {
    match edn {
        Edn::Tagged(t, v) if *t == tag => match v.as_ref() {
            Edn::Str(s) => Ok(s),
            _ => Err(ParseError::InvalidMoleculeMap(format!("#{tag} value must be a string"))),
        },
        _ => Err(ParseError::InvalidMoleculeMap(format!("expected #{tag} tagged literal"))),
    }
}

fn parse_bond_spec(edn: &Edn<'_>) -> Result<BondSpec, ParseError> {
    match edn {
        Edn::Tagged(t, v) if *t == "bond" => match v.as_ref() {
            Edn::Str(s) => Ok(BondSpec::Literal(parse_bond_dsl(s)?)),
            _ => Err(ParseError::InvalidMoleculeMap("#bond must be followed by a string".into())),
        },
        Edn::Key("single") => Ok(BondSpec::Single),
        Edn::Key("double") => Ok(BondSpec::Double),
        Edn::Key("triple") => Ok(BondSpec::Triple),
        Edn::Key("quadruple") => Ok(BondSpec::Quadruple),
        _ => Err(ParseError::InvalidMoleculeMap(
            "bond spec must be #bond \"...\" or a keyword shorthand".into(),
        )),
    }
}

fn opt_list<'e, T>(
    map: &'e BTreeMap<Edn<'e>, Edn<'e>>,
    key: &str,
    f: impl Fn(&'e Edn<'e>) -> Result<T, ParseError>,
) -> Result<Vec<T>, ParseError> {
    match map_get(map, key) {
        None => Ok(Vec::new()),
        Some(Edn::Vector(v)) => v.iter().map(f).collect(),
        Some(_) => Err(ParseError::InvalidMoleculeMap(format!(":{key} must be a vector"))),
    }
}

// ── section parsers ───────────────────────────────────────────────────────────

fn parse_atom_collection(edn: &Edn<'_>) -> Result<AtomCollection, ParseError> {
    match edn {
        Edn::Map(m) => {
            let mut named = IndexMap::new();
            for (k, v) in m {
                let label = require_label(k)?;
                let atom = parse_atom_dsl(require_tagged_str(v, "atom")?)?;
                named.insert(label, atom);
            }
            Ok(AtomCollection::Named(named))
        }
        Edn::Vector(items) => {
            let atoms = items
                .iter()
                .map(|e| parse_atom_dsl(require_tagged_str(e, "atom")?))
                .collect::<Result<_, _>>()?;
            Ok(AtomCollection::Indexed(atoms))
        }
        _ => Err(ParseError::InvalidMoleculeMap(":atoms must be a map or vector".into())),
    }
}

fn parse_covalent_bonds(edn: &Edn<'_>) -> Result<Vec<CovalentBondEntry>, ParseError> {
    match edn {
        Edn::Vector(v) => v.iter().map(parse_covalent_entry).collect(),
        _ => Err(ParseError::InvalidMoleculeMap(":bonds must be a vector".into())),
    }
}

fn parse_covalent_entry(edn: &Edn<'_>) -> Result<CovalentBondEntry, ParseError> {
    let m = require_map(edn, "covalent bond entry")?;
    Ok(CovalentBondEntry {
        id: require_label(
            map_get(m, "id").ok_or_else(|| ParseError::MissingKey(":id in bond entry".into()))?,
        )?,
        a: require_label(
            map_get(m, "a").ok_or_else(|| ParseError::MissingKey(":a in bond entry".into()))?,
        )?,
        b: require_label(
            map_get(m, "b").ok_or_else(|| ParseError::MissingKey(":b in bond entry".into()))?,
        )?,
        bond: parse_bond_spec(
            map_get(m, "bond")
                .ok_or_else(|| ParseError::MissingKey(":bond in bond entry".into()))?,
        )?,
    })
}

fn parse_dative_entry(edn: &Edn<'_>) -> Result<DativeBondEntry, ParseError> {
    let m = require_map(edn, "dative bond entry")?;
    Ok(DativeBondEntry {
        id: require_label(
            map_get(m, "id")
                .ok_or_else(|| ParseError::MissingKey(":id in dative entry".into()))?,
        )?,
        donor: require_label(
            map_get(m, "donor")
                .ok_or_else(|| ParseError::MissingKey(":donor in dative entry".into()))?,
        )?,
        acceptor: require_label(
            map_get(m, "acceptor")
                .ok_or_else(|| ParseError::MissingKey(":acceptor in dative entry".into()))?,
        )?,
        bond: parse_bond_spec(
            map_get(m, "bond")
                .ok_or_else(|| ParseError::MissingKey(":bond in dative entry".into()))?,
        )?,
    })
}

fn parse_aromatic_entry(edn: &Edn<'_>) -> Result<AromaticEntry, ParseError> {
    let m = require_map(edn, "aromatic entry")?;
    let id = require_label(
        map_get(m, "id")
            .ok_or_else(|| ParseError::MissingKey(":id in aromatic entry".into()))?,
    )?;
    let atoms = match map_get(m, "atoms")
        .ok_or_else(|| ParseError::MissingKey(":atoms in aromatic entry".into()))?
    {
        Edn::Vector(v) => v.iter().map(require_label).collect::<Result<_, _>>()?,
        _ => {
            return Err(ParseError::InvalidMoleculeMap(
                ":atoms in aromatic entry must be a vector".into(),
            ))
        }
    };
    Ok(AromaticEntry { id, atoms })
}

fn parse_multicenter_entry(edn: &Edn<'_>) -> Result<MulticenterEntry, ParseError> {
    let m = require_map(edn, "multicenter entry")?;
    let id = require_label(
        map_get(m, "id")
            .ok_or_else(|| ParseError::MissingKey(":id in multicenter entry".into()))?,
    )?;
    let atoms = match map_get(m, "atoms")
        .ok_or_else(|| ParseError::MissingKey(":atoms in multicenter entry".into()))?
    {
        Edn::Vector(v) => v.iter().map(require_label).collect::<Result<_, _>>()?,
        _ => {
            return Err(ParseError::InvalidMoleculeMap(
                ":atoms in multicenter entry must be a vector".into(),
            ))
        }
    };
    Ok(MulticenterEntry { id, atoms })
}

fn parse_noncovalent_entry(edn: &Edn<'_>) -> Result<NoncovalentEntry, ParseError> {
    let m = require_map(edn, "noncovalent bond entry")?;
    Ok(NoncovalentEntry {
        id: require_label(
            map_get(m, "id")
                .ok_or_else(|| ParseError::MissingKey(":id in noncovalent entry".into()))?,
        )?,
        a: require_label(
            map_get(m, "a")
                .ok_or_else(|| ParseError::MissingKey(":a in noncovalent entry".into()))?,
        )?,
        b: require_label(
            map_get(m, "b")
                .ok_or_else(|| ParseError::MissingKey(":b in noncovalent entry".into()))?,
        )?,
        bond: parse_bond_spec(
            map_get(m, "bond")
                .ok_or_else(|| ParseError::MissingKey(":bond in noncovalent entry".into()))?,
        )?,
    })
}

// ── validation ────────────────────────────────────────────────────────────────

fn validate(ast: &MoleculeMapAst) -> Result<(), ParseError> {
    let atom_labels: HashSet<String> = match &ast.atoms {
        AtomCollection::Named(m) => m.keys().cloned().collect(),
        AtomCollection::Indexed(v) => (0..v.len()).map(|i| i.to_string()).collect(),
    };

    let mut seen_ids: HashSet<&str> = HashSet::new();
    for id in ast
        .bonds
        .iter()
        .map(|e| e.id.as_str())
        .chain(ast.dative.iter().map(|e| e.id.as_str()))
        .chain(ast.aromatic.iter().map(|e| e.id.as_str()))
        .chain(ast.multicenter.iter().map(|e| e.id.as_str()))
        .chain(ast.noncovalent.iter().map(|e| e.id.as_str()))
    {
        if !seen_ids.insert(id) {
            return Err(ParseError::DuplicateId(id.to_string()));
        }
    }

    let check = |label: &str| -> Result<(), ParseError> {
        if atom_labels.contains(label) {
            Ok(())
        } else {
            Err(ParseError::UnknownEndpoint(label.to_string()))
        }
    };

    for e in &ast.bonds {
        check(&e.a)?;
        check(&e.b)?;
    }
    for e in &ast.dative {
        check(&e.donor)?;
        check(&e.acceptor)?;
    }
    for e in &ast.aromatic {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.multicenter {
        for a in &e.atoms {
            check(a)?;
        }
    }
    for e in &ast.noncovalent {
        check(&e.a)?;
        check(&e.b)?;
    }

    Ok(())
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::super::atom::ElementExpr;
    use super::super::value::ValueAst;
    use super::*;

    fn atom_bare(element: Element) -> AtomAst {
        use super::super::atom::AtomAst;
        AtomAst {
            element: ElementExpr::Lit(element),
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    #[test]
    fn test_minimal_named() {
        let result = parse_molecule_map(r#"{:atoms {:C #atom "C"} :bonds []}"#).unwrap();
        assert!(matches!(result.atoms, AtomCollection::Named(_)));
        if let AtomCollection::Named(ref m) = result.atoms {
            assert_eq!(m.len(), 1);
            assert_eq!(m["C"], atom_bare(Element::C));
        }
        assert!(result.bonds.is_empty());
        assert!(result.dative.is_empty());
        assert!(result.aromatic.is_empty());
        assert_eq!(result.charge, None);
    }

    #[test]
    fn test_single_bond_keyword() {
        let result = parse_molecule_map(
            r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :C :b :O :bond :single}]}"#,
        )
        .unwrap();
        assert_eq!(
            result.bonds,
            vec![CovalentBondEntry {
                id: "b1".into(),
                a: "C".into(),
                b: "O".into(),
                bond: BondSpec::Single,
            }]
        );
    }

    #[test]
    fn test_bond_literal() {
        use super::super::bond::BondAst;
        let result = parse_molecule_map(
            r#"{:atoms {:C #atom "C" :O #atom "O"} :bonds [{:id :b1 :a :C :b :O :bond #bond "2"}]}"#,
        )
        .unwrap();
        assert_eq!(
            result.bonds[0].bond,
            BondSpec::Literal(BondAst {
                order: ValueAst::Lit(2),
                charge: None,
                unpaired_electrons: None,
                multiplicity: None,
            })
        );
    }

    #[test]
    fn test_all_bond_keywords() {
        let result = parse_molecule_map(
            r#"{:atoms {:A #atom "C" :B #atom "C" :C2 #atom "C" :D #atom "C"
                        :E #atom "C" :F #atom "C" :G #atom "C" :H2 #atom "C"}
                :bonds [{:id :b1 :a :A :b :B :bond :single}
                        {:id :b2 :a :C2 :b :D :bond :double}
                        {:id :b3 :a :E :b :F :bond :triple}
                        {:id :b4 :a :G :b :H2 :bond :quadruple}]}"#,
        )
        .unwrap();
        assert_eq!(result.bonds[0].bond, BondSpec::Single);
        assert_eq!(result.bonds[1].bond, BondSpec::Double);
        assert_eq!(result.bonds[2].bond, BondSpec::Triple);
        assert_eq!(result.bonds[3].bond, BondSpec::Quadruple);
    }

    #[test]
    fn test_charge_field() {
        let result = parse_molecule_map(
            r#"{:atoms {:N #atom "N"} :bonds [] :charge -1}"#,
        )
        .unwrap();
        assert_eq!(result.charge, Some(-1));
    }

    #[test]
    fn test_indexed_atoms() {
        let result = parse_molecule_map(
            r#"{:atoms [#atom "C" #atom "O"] :bonds [{:id :b1 :a :0 :b :1 :bond :single}]}"#,
        )
        .unwrap();
        assert!(matches!(result.atoms, AtomCollection::Indexed(_)));
        if let AtomCollection::Indexed(ref v) = result.atoms {
            assert_eq!(v.len(), 2);
            assert_eq!(v[0], atom_bare(Element::C));
            assert_eq!(v[1], atom_bare(Element::O));
        }
        assert_eq!(result.bonds[0].a, "0");
        assert_eq!(result.bonds[0].b, "1");
    }

    #[test]
    fn test_dative_section() {
        let result = parse_molecule_map(
            r#"{:atoms {:B #atom "B" :N #atom "N"}
                :bonds []
                :dative [{:id :d1 :donor :N :acceptor :B :bond :single}]}"#,
        )
        .unwrap();
        assert_eq!(
            result.dative,
            vec![DativeBondEntry {
                id: "d1".into(),
                donor: "N".into(),
                acceptor: "B".into(),
                bond: BondSpec::Single,
            }]
        );
    }

    #[test]
    fn test_aromatic_section() {
        let result = parse_molecule_map(
            r#"{:atoms {:C1 #atom "C" :C2 #atom "C" :C3 #atom "C"
                        :C4 #atom "C" :C5 #atom "C" :C6 #atom "C"}
                :bonds [{:id :b1 :a :C1 :b :C2 :bond :single}
                        {:id :b2 :a :C2 :b :C3 :bond :double}
                        {:id :b3 :a :C3 :b :C4 :bond :single}
                        {:id :b4 :a :C4 :b :C5 :bond :double}
                        {:id :b5 :a :C5 :b :C6 :bond :single}
                        {:id :b6 :a :C6 :b :C1 :bond :double}]
                :aromatic [{:id :ar1 :atoms [:C1 :C2 :C3 :C4 :C5 :C6]}]}"#,
        )
        .unwrap();
        assert_eq!(result.aromatic.len(), 1);
        assert_eq!(result.aromatic[0].id, "ar1");
        assert_eq!(
            result.aromatic[0].atoms,
            vec!["C1", "C2", "C3", "C4", "C5", "C6"]
        );
    }

    #[test]
    fn test_methanol() {
        use super::super::atom::HydrogenExpr;
        let result = parse_molecule_map(
            r#"{:atoms {:C #atom "C#h3" :O #atom "O#h1" :H #atom "H"}
                :bonds [{:id :b1 :a :C :b :O :bond :single}
                        {:id :b2 :a :O :b :H :bond :single}]}"#,
        )
        .unwrap();
        if let AtomCollection::Named(ref m) = result.atoms {
            assert_eq!(
                m["C"].implicit_hydrogens,
                Some(HydrogenExpr::Value(ValueAst::Lit(3)))
            );
            assert_eq!(
                m["O"].implicit_hydrogens,
                Some(HydrogenExpr::Value(ValueAst::Lit(1)))
            );
            assert_eq!(m["H"].element, ElementExpr::Lit(Element::H));
        } else {
            panic!("expected named atom collection");
        }
        assert_eq!(result.bonds.len(), 2);
    }

    #[rstest]
    #[case::non_map("42", ParseError::InvalidMoleculeMap("expected EDN map for top level".into()))]
    #[case::missing_atoms(r#"{:bonds []}"#, ParseError::MissingKey(":atoms".into()))]
    #[case::missing_bonds(r#"{:atoms {:C #atom "C"}}"#, ParseError::MissingKey(":bonds".into()))]
    #[case::unknown_endpoint(
        r#"{:atoms {:C #atom "C"} :bonds [{:id :b1 :a :C :b :X :bond :single}]}"#,
        ParseError::UnknownEndpoint("X".into())
    )]
    #[case::duplicate_id(
        r#"{:atoms {:C #atom "C" :O #atom "O" :N #atom "N"}
            :bonds [{:id :b1 :a :C :b :O :bond :single}
                    {:id :b1 :a :O :b :N :bond :single}]}"#,
        ParseError::DuplicateId("b1".into())
    )]
    #[case::bad_atom_string(
        "{:atoms {:X #atom \"#h3\"} :bonds []}",
        ParseError::InvalidAtomElement("#h3".into())
    )]
    fn test_parse_molecule_map_invalid(#[case] input: &str, #[case] expected: ParseError) {
        let result = parse_molecule_map(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        assert_eq!(result.unwrap_err(), expected, "for input {input:?}");
    }
}
