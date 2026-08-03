use std::fs;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_edn::read_string;

use super::*;
use crate::ast::atom::AtomAst;
use crate::ast::boolean::BooleanAst;
use crate::ast::constraint::{BondConstraintAst, Constraint, MoleculeConstraint};
use crate::ast::electrons::ElectronCountsAst;
use crate::ast::spin::UnpairedElectronsAst;
use crate::ast::value::ValueAst;
use crate::mol_dsl;

#[fixture]
fn populated_molecule_dsl() -> MoleculeDsl {
    r#"{
        :atoms [[:a "C"] "F" "Cl" "Br" "I"]
        :bonds [
            {:id :b :atoms [0 1] :type "2"}
            [0 2 "1"]
            [0 3 "1"]
            [0 4 "1"]
        ]
        :dative-bonds [{:id :d :donors [1] :acceptor 0 :type "1#R"}]
        :aromatic-systems [{:id :ar :atoms [0 1] :type "*#e2"}]
        :multicenter-bonds [{:id :m :atoms [0 1] :type "*#e2"}]
        :noncovalent-bonds [{:id :n :atoms [0 1] :type "Hbd"}]
        :stereo-atoms [{:id :sa :site 0 :ligands [1 2 3 4] :type "Th1"}]
        :stereo-bonds [{:id :sb :site 0 :ligands [2 3] :type "Ct1"}]
        :atom-aliases [:x "O"]
    }"#
    .parse()
    .unwrap()
}

/// Every `fuzz_molecule` seed must parse (tree and streaming) — guards the seed corpus against
/// rot as the molecule DSL evolves.
#[rstest]
fn test_fuzz_molecule_seeds_valid() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fuzz/seeds/fuzz_molecule");
    let mut count = 0;
    let mut failures: Vec<String> = Vec::new();
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        let data = fs::read_to_string(&path).unwrap();
        if let Err(e) = MoleculeDsl::from_edn_str(&data) {
            failures.push(format!("{name}: stream: {e:?}"));
        }
        let tree = read_string(&data)
            .ok()
            .and_then(|edn| MoleculeDsl::from_edn(&edn).ok());
        if tree.is_none() {
            failures.push(format!("{name}: tree parse failed"));
        }
        count += 1;
    }
    assert!(
        failures.is_empty(),
        "invalid seeds:\n{}",
        failures.join("\n")
    );
    assert_eq!(count, 28);
}

#[rstest]
#[case::empty(None)]
#[case::atom(Some(Entity::Atom(AtomId(4))))]
#[case::bond(Some(Entity::Bond(BondId(3))))]
#[case::dative_bond(Some(Entity::DativeBond(DativeBondId(0))))]
#[case::aromatic_system(Some(Entity::AromaticSystem(AromaticSystemId(0))))]
#[case::multicenter_bond(Some(Entity::MulticenterBond(MulticenterBondId(0))))]
#[case::noncovalent_bond(Some(Entity::NoncovalentBond(NoncovalentBondId(0))))]
#[case::stereo_atom(Some(Entity::StereoAtom(StereoAtomId(0))))]
#[case::stereo_bond(Some(Entity::StereoBond(StereoBondId(0))))]
fn test_molecule_dsl_new(populated_molecule_dsl: MoleculeDsl, #[case] entity: Option<Entity>) {
    let ast = populated_molecule_dsl.into_parts().0;
    let mut metadata = MoleculeMetadata::new();
    if let Some(entity) = entity {
        metadata.set_keyword(entity, "key").unwrap();
    }

    let actual = MoleculeDsl::new(ast.clone(), metadata.clone()).unwrap();

    assert_eq!(actual.into_parts(), (ast, metadata));
}

#[rstest]
#[case::atom(Entity::Atom(AtomId(5)))]
#[case::bond(Entity::Bond(BondId(4)))]
#[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
#[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
#[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
#[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
fn test_molecule_dsl_new_error(populated_molecule_dsl: MoleculeDsl, #[case] entity: Entity) {
    let ast = populated_molecule_dsl.into_parts().0;
    let mut metadata = MoleculeMetadata::new();
    metadata.set_keyword(entity, "key").unwrap();

    assert_eq!(
        MoleculeDsl::new(ast, metadata),
        Err(MetadataError::EntityOutOfRange(entity))
    );
}

#[rstest]
fn test_molecule_dsl_new_parsed(populated_molecule_dsl: MoleculeDsl) {
    let expected = populated_molecule_dsl.clone();
    let (ast, metadata) = populated_molecule_dsl.into_parts();

    assert_eq!(MoleculeDsl::new(ast, metadata).unwrap(), expected);
}

#[rstest]
#[case::empty("{:atoms [] :bonds []}", MoleculeAst::default())]
#[case::two_atoms_one_bond(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#, MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C); 2], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }))]
#[case::atom_with_keyword(r#"{:atoms [[:c1 "C"] "C"] :bonds []}"#, MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C); 2], bonds: vec![], ..Default::default() }))]
#[case::bond_with_id_field(r#"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type :single}]}"#, MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C); 2], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }))]
#[case::atom_alias(r#"{:atoms [:x :x] :bonds [[0 1 "1"]] :atom-aliases [:x "C"]}"#, MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C); 2], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }))]
fn test_mol_dsl_to_edn(#[case] input: &str, #[case] expected: MoleculeAst) {
    let dsl = input.parse::<MoleculeDsl>().unwrap();
    assert_eq!(dsl.into_ast(&MoleculeDefaults::default()), expected);
}

#[rstest]
fn test_molecule_dsl_display_to_edn_parity() {
    let dsl = r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#.parse::<MoleculeDsl>().unwrap();
    assert_eq!(dsl.to_string(), dsl.to_edn().to_string());
}

#[rstest]
fn test_molecule_dsl_to_edn_omits_empty_optional_sections() {
    let dsl = r#"{:atoms ["C"] :bonds []}"#.parse::<MoleculeDsl>().unwrap();
    let edn = dsl.to_edn();
    let Edn::Map(m) = &edn else {
        panic!("expected map");
    };
    assert!(m.get_keyword("dative-bonds").is_none());
    assert!(m.get_keyword("aromatic-systems").is_none());
    assert!(m.get_keyword("multicenter-bonds").is_none());
    assert!(m.get_keyword("noncovalent-bonds").is_none());
    assert!(m.get_keyword("atom-aliases").is_none());
    assert!(m.get_keyword("constraints").is_none());
}

#[rstest]
fn test_molecule_dsl_from_edn_empty() {
    let edn = read_string("{:atoms [] :bonds []}").unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(*dsl.ast(), MoleculeAst::default());
}

#[rstest]
fn test_molecule_dsl_from_edn_two_atoms_one_bond() {
    let edn = read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(MoleculeDsl::from_edn(&dsl.to_edn()).unwrap(), dsl);
}

#[rstest]
fn test_molecule_dsl_from_edn_atom_with_keyword() {
    let edn = read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.metadata().keyword(Entity::Atom(AtomId(0))), Some("c1"));
    assert_eq!(dsl.metadata().keyword(Entity::Atom(AtomId(1))), None);
}

#[rstest]
fn test_molecule_dsl_from_edn_bond_map_form_with_id_field() {
    let edn =
        read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]}"##).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.ast().bonds().count(), 1);
    assert_eq!(dsl.metadata().keyword(Entity::Bond(BondId(0))), Some("b1"));
}

#[rstest]
fn test_molecule_dsl_from_edn_atom_aliases() {
    let edn = read_string(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.ast().atoms().count(), 2);
    assert_eq!(
        dsl.metadata().atom_alias("x"),
        Some(&AtomDsl(AtomAst::from_element(Element::C)))
    );
}

#[rstest]
fn test_molecule_dsl_tree_streaming_metadata_equivalence() {
    let source = r##"{:atoms [[:a "C"] :x] :bonds [{:id :b :atoms [:a 1] :type :single}] :atom-aliases [:x "N"]}"##;
    let tree = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
    let streaming = MoleculeDsl::from_edn_str(source).unwrap();

    assert_eq!(streaming.metadata(), tree.metadata());
    assert_eq!(streaming.ast(), tree.ast());
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
fn test_molecule_dsl_edn_roundtrip_with_keywords_and_aliases() {
    let source = r##"{:atoms [[:a "C"] [:b "C"] :x] :bonds [{:id :b1 :atoms [:a :b] :type :single} [:b 2 :double]] :atom-aliases [:x "N"]}"##;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    let rendered = dsl.to_edn();
    assert_eq!(rendered, edn);
}

#[rustfmt::skip]
#[rstest]
#[case::empty(r##"{:atoms [] :bonds []}"##)]
#[case::small(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##)]
#[case::with_keywords(r##"{:atoms [[:a "C"] [:b "N"]] :bonds [{:id :b1 :atoms [:a :b] :type "1"}]}"##)]
#[case::inline_atom_constraints(r##"{:atoms ["C#v4" "N#R+"] :bonds []}"##)]
#[case::inline_bond_constraint(r##"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"##)]
#[case::aromatic_section(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:id :ar1 :atoms [0 1 2 3 4 5] :type "*#e6"}]}"##)]
#[case::multicenter_section(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]}"##)]
#[case::dative_section(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1#R"}]}"##)]
#[case::dative_multi_donor(r##"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1#R"}]}"##)]
#[case::noncovalent_section(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::stereo_atom_section(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_id_field_and_keyword_type(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :type :ccw}]}"##)]
#[case::stereo_atom_virtual_ligand(r##"{:atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 [:lp 0] [:h 0]] :type "Th1"}]}"##)]
#[case::stereo_bond_section(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_keyword_ref(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:id :db :atoms [1 2] :type "2"} [2 3 "1"]] :stereo-bonds [{:site :db :ligands [0 3] :type :e}]}"##)]
#[case::atom_aliases(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##)]
#[case::constraints_connected(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"##)]
#[case::constraints_bond_order_sum(r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"} {:id :b2 :atoms [1 2] :type "1"}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##)]
#[case::constraints_atom_leaf_in_not(r##"{:atoms [[:c1 "C"]] :bonds [] :constraints [{:not {:atom [:c1 {:valence 3}]}}]}"##)]
#[case::constraints_nested_combinators(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:and [{:or [{:atom [0 {:valence 3}]} {:atom [0 {:valence 4}]}]} {:not {:connected {:atoms [0 1]}}}]}]}"##)]
#[case::constraints_sub_pattern(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}]}"##)]
#[case::constraints_atom_degree(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:atom [0 {:degree 3}]}]}"##)]
#[case::constraints_atom_total_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-degree 4}]}]}"##)]
#[case::constraints_atom_total_hydrogens(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-hydrogens 2}]}]}"##)]
#[case::constraints_atom_ring_count(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::constraints_atom_ring_size(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::constraints_atom_ring_size_set(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 6 :count [5 6]}}]}]}"##)]
#[case::constraints_atom_ring_size_conj(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 5 :count 1}}]} {:atom [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::constraints_atom_total_valence(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-valence 4}]}]}"##)]
#[case::constraints_atom_ring_valence(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-valence 2}]}]}"##)]
#[case::constraints_atom_ring_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-degree 2}]}]}"##)]
#[case::constraints_atom_donated_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:donated-pairs 1}]}]}"##)]
#[case::constraints_atom_accepted_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:accepted-pairs 1}]}]}"##)]
#[case::constraints_atom_aromatic_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence :not-aromatic}]}]}"##)]
#[case::constraints_atom_aromatic_valence_with_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence {:aromatic 6}}]}]}"##)]
#[case::constraints_atom_multicenter_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence :not-multicenter}]}]}"##)]
#[case::constraints_atom_multicenter_valence_with_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence {:multicenter 3}}]}]}"##)]
#[case::constraints_bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:aromatic true}]}]}"##)]
#[case::constraints_bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::constraints_bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::constraints_dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::constraints_dative_donors(r##"{:atoms ["C" "C" "N"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1#R"}] :constraints [{:dative-bond-donors [0 [0 1]]}]}"##)]
#[case::constraints_dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
#[case::constraints_dative_contains_all_donors(r##"{:atoms ["C" "C" "N"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1#R"}] :constraints [{:dative-bond-contains-all-donors [0 [0 1]]}]}"##)]
#[case::constraints_dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
#[case::constraints_dative_any_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-any-donor [0 {:valence 3}]}]}"##)]
#[case::constraints_dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
#[case::constraints_dative_all_donors(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-all-donors [0 {:valence 3}]}]}"##)]
#[case::constraints_dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor-satisfies [0 {:valence 3}]}]}"##)]
#[case::constraints_aromatic_system_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
#[case::constraints_aromatic_system_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
#[case::constraints_aromatic_system_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
#[case::constraints_aromatic_system_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}] :constraints [{:aromatic-system-any-atom [0 {:valence 4}]}]}"##)]
#[case::constraints_multicenter_contains(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*#e3"}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
#[case::constraints_multicenter_contains_all(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*#e3"}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
#[case::constraints_multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
#[case::constraints_multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}] :constraints [{:multicenter-bond-any-atom [0 {:valence 4}]}]}"##)]
#[case::constraints_noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
#[case::constraints_noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
#[case::constraints_noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 3} {:valence 1}]]}]}"##)]
#[case::constraints_sub_pattern_multi_entity_anchor(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]] :bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [[0 1 "1"]]}}}]}"##)]
#[case::constraints_sub_pattern_dative_anchor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:sub-pattern {:anchor {:dative-bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]}}}]}"##)]
#[case::constraints_sub_pattern_aromatic_system_anchor(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}] :constraints [{:sub-pattern {:anchor {:aromatic-systems [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}]}}}]}"##)]
#[case::constraints_sub_pattern_multicenter_anchor(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}] :constraints [{:sub-pattern {:anchor {:multicenter-bonds [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]}}}]}"##)]
#[case::constraints_sub_pattern_noncovalent_anchor(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:sub-pattern {:anchor {:noncovalent-bonds [[0 0]]} :pattern {:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}}}]}"##)]
#[case::constraints_sub_pattern_stereo_atom_anchor(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:sub-pattern {:anchor {:stereo-atoms [[0 0]]} :pattern {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}}}]}"##)]
#[case::constraints_sub_pattern_stereo_bond_anchor(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}] :constraints [{:sub-pattern {:anchor {:stereo-bonds [[0 0]]} :pattern {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}}}]}"##)]
#[case::constraints_atom_tetrahedral_stereo_keyword(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:tetrahedral-stereo :not-stereo}]}]}"##)]
#[case::constraints_atom_tetrahedral_stereo_lit(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:tetrahedral-stereo {:stereo 1}}]}]}"##)]
#[case::constraints_atom_tetrahedral_stereo_coset_undetermined(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:tetrahedral-stereo {:stereo :undetermined}}]}]}"##)]
#[case::constraints_atom_tetrahedral_stereo_set(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:tetrahedral-stereo {:stereo [1 2]}}]}]}"##)]
#[case::constraints_atom_tetrahedral_stereo_expr(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:tetrahedral-stereo {:stereo "~1"}}]}]}"##)]
#[case::constraints_bond_cis_trans_stereo(r##"{:atoms ["C" "C"] :bonds [[0 1 "2"]] :constraints [{:bond [0 {:cis-trans-stereo {:stereo 1}}]}]}"##)]
#[case::stereo_atom_inline_constraints(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1#f(0,1,2)#g/"}]}"##)]
#[case::stereo_atom_molecule_constraint(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [0 [:tetrahedral {:ligand-symmetry {:permutation [[0 1]] :orientation :improper :invariant false}}]]}]}"##)]
#[case::stereo_atom_fluxionality_constraint(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [0 [:tetrahedral {:fluxionality {:permutation [[0 1 2]]}}]]}]}"##)]
#[case::stereo_bond_molecule_constraint(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}] :constraints [{:stereo-bond [0 [:cis-trans {:topicity {:pair [0 1] :relation :diastereotopic}}]]}]}"##)]
fn test_molecule_dsl_from_edn_str_from_edn_parity(#[case] source: &str) {
    let via_str = MoleculeDsl::from_edn_str(source).unwrap();
    let tree = read_string(source).unwrap();
    let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
    assert_eq!(via_str, via_tree);
}

#[rustfmt::skip]
#[rstest]
#[case::stereo_atom(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_id_field_virtual_ligands(r##"{:atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 [:lp 0] [:h 0]] :type "Th1"}]}"##)]
#[case::stereo_bond_keyword_ref(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:id :db :atoms [1 2] :type "2"} [2 3 "1"]] :stereo-bonds [{:site :db :ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_id_field(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :sb1 :site 1 :ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_atom_inline_constraints(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1#f(0,1,2)#g/"}]}"##)]
#[case::stereo_atom_molecule_constraint(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [0 [:tetrahedral {:stereogenicity {:relation :stereogenic}}]]}]}"##)]
#[case::stereo_bond_inline_constraints(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1#g/"}]}"##)]
#[case::stereo_bond_molecule_constraint(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}] :constraints [{:stereo-bond [0 [:cis-trans {:stereogenicity {:relation :stereogenic}}]]}]}"##)]
fn test_molecule_dsl_stereo_edn_roundtrip(#[case] source: &str) {
    let dsl = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
    let reparsed = MoleculeDsl::from_edn(&dsl.to_edn()).unwrap();
    assert_eq!(reparsed, dsl);
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
#[case::dative_constraint_unknown_key(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:bogus 1}]}]}"##)]
#[case::sub_pattern_anchor_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:bogus [[0 0]]} :pattern {:atoms ["C"] :bonds []}}}]}"##)]
#[case::constraint_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:bogus 1}]}"##)]
#[case::noncovalent_ends_satisfy_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 2}]]}]}"##)]
#[case::noncovalent_ends_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0]]}]}"##)]
fn test_molecule_dsl_from_edn_str_invalid_constraints(#[case] source: &str) {
    let result = MoleculeDsl::from_edn_str(source);
    assert!(
        result.is_err(),
        "expected parse failure, but got: {:?}",
        result,
    );
}

#[rstest]
fn test_molecule_dsl_from_str() {
    let source = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;
    let dsl: MoleculeDsl = source.parse().unwrap();
    assert_eq!(
        dsl,
        MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap()
    );
}

#[rstest]
fn test_molecule_dsl_from_str_error() {
    let err = "not a map".parse::<MoleculeDsl>().unwrap_err();
    assert!(matches!(err, ParseError::EdnParse(_)));
}

// Round-trip direction: DSL → AST (raise) → DSL (lower) is the
// identity. AST → DSL → AST isn't, since raising `Undetermined`
// fields to `Lit(0)` is one-way under `zeroed()`. One case per overlay
// kind so the `into_ast` / `from_ast` per-relation loops are exercised.
#[rustfmt::skip]
#[rstest]
#[case::atoms_bonds(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##)]
#[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}]}"##)]
#[case::aromatic(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}]}"##)]
#[case::multicenter(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*"}]}"##)]
#[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::stereo_atom(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_bond(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"##)]
fn test_molecule_dsl_dsl_to_ast_to_dsl_roundtrip_zeroed(#[case] source: &str) {
    let ast = mol_dsl!(source);
    let dsl = MoleculeDsl::new(ast, MoleculeMetadata::default()).unwrap();
    let cfg = MoleculeDefaults::zeroed();
    let raised = dsl.clone().into_ast(&cfg);
    let lowered = MoleculeDsl::from_ast(&raised, &cfg);
    assert_eq!(lowered.ast(), dsl.ast());
}

#[rstest]
fn test_molecule_dsl_from_ast_has_empty_metadata() {
    let ast = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
    let cfg = MoleculeDefaults::zeroed();
    let dsl = MoleculeDsl::from_ast(&ast, &cfg);
    assert_eq!(dsl.metadata(), &MoleculeMetadata::default());
}

#[rustfmt::skip]
#[rstest]
#[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}]}"##)]
#[case::dative_with_id_and_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1#R"}]}"##)]
#[case::dative_multi_donor(r##"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type :single}]}"##)]
#[case::aromatic_minimal(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*"}]}"##)]
#[case::aromatic_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:id :a1 :atoms [0 1] :type "*#e6"}]}"##)]
#[case::aromatic_with_electrons_literals(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]}"##)]
#[case::aromatic_with_electrons_and_total(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]#e6"}]}"##)]
#[case::multicenter_minimal(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*"}]}"##)]
#[case::multicenter_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:id :m1 :atoms [0 1] :type "*#e2"}]}"##)]
#[case::multicenter_with_electrons_literals(r##"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[1,0,1]"}]}"##)]
#[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::noncovalent_with_id(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:id :n1 :atoms [0 1] :type "Hbd"}]}"##)]
fn test_molecule_dsl_edn_roundtrip_non_localized_entities(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

#[rustfmt::skip]
#[rstest]
#[case::aromatic(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}]}"#)]
#[case::multicenter(r#"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}]}"#)]
fn test_molecule_dsl_edn_parse_electrons_undetermined(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    let ast = dsl.ast();
    let electrons = ast
        .aromatic_systems()
        .iter()
        .map(|v| v.ast.electrons.clone())
        .chain(
            ast.multicenter_bonds()
                .iter()
                .map(|v| v.ast.electrons.clone()),
        )
        .next()
        .unwrap();
    assert_eq!(electrons, ElectronCountsAst::Undetermined);
}

#[rstest]
fn test_molecule_dsl_edn_roundtrip_connected_constraint() {
    let source =
        r##"{:atoms ["C" "C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1 2]}}]}"##;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

/// `to_edn` drops vacuous molecule-level constraints — `ChargeSum`/`BondOrderSum`
/// with `Undetermined` sum, `UnpairedElectronCoupling` with undetermined
/// unpaired electrons — during lowering;
/// the surviving molecule constraints are exactly those left after render → reparse.
#[rustfmt::skip]
#[rstest]
#[case::charge_sum_vacuous(vec![MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Undetermined }], vec![])]
#[case::bond_order_sum_vacuous(vec![MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Undetermined }], vec![])]
#[case::unpaired_electron_coupling_vacuous(vec![MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::default() }], vec![])]
#[case::vacuous_dropped_concrete_kept(
    vec![
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Undetermined },
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
    ],
    vec![MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }])]
fn test_molecule_dsl_to_edn_vacuous_constraints(
    #[case] pushed: Vec<MoleculeConstraint>,
    #[case] expected: Vec<MoleculeConstraint>,
) {
    let mut ast = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
    for c in pushed {
        ast.constraints_mut().push(Constraint::Molecule(c));
    }
    let dsl = MoleculeDsl::new(ast, MoleculeMetadata::default()).unwrap();
    let reparsed = MoleculeAst::from_edn(&dsl.to_edn()).unwrap();
    let surviving: Vec<MoleculeConstraint> = reparsed
        .constraints()
        .iter()
        .filter_map(|c| match c {
            Constraint::Molecule(m) => Some(m.clone()),
            _ => None,
        })
        .collect();
    assert_eq!(surviving, expected);
}

#[rstest]
fn test_molecule_dsl_edn_roundtrip_connected_all_atoms() {
    let source = r#"{:atoms ["C" "C" "C"] :bonds [] :constraints [{:connected {}}]}"#;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

#[rstest]
fn test_molecule_ast_from_edn_structural_bond_ref() {
    // A structural bond ref ({:atoms [0 1]}) names the bond by its endpoints, resolved against the
    // namespace's participant lookup.
    let source = r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#;
    let ast = MoleculeAst::from_edn(&read_string(source).unwrap()).unwrap();
    let constraints: Vec<Constraint> = ast.constraints().iter().cloned().collect();
    assert_eq!(
        constraints,
        vec![Constraint::Bond(
            BondId(0),
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
        )]
    );
}

// A structural ref is input-only: parsing one and rendering yields the keyword-else-positional form
// (never structural), so it renders identically to the equivalent keyword / positional input.
#[rstest]
#[case::unnamed_bond_renders_positional(
    r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#,
    r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:aromatic true}]}]}"#,
)]
#[case::named_bond_renders_keyword(
    r#"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"}] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#,
    r#"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"}] :constraints [{:bond [:b1 {:aromatic true}]}]}"#,
)]
fn test_molecule_dsl_to_edn_structural_ref(#[case] structural: &str, #[case] canonical: &str) {
    let via_structural = MoleculeDsl::from_edn(&read_string(structural).unwrap()).unwrap();
    let via_canonical = MoleculeDsl::from_edn(&read_string(canonical).unwrap()).unwrap();
    assert_eq!(
        via_structural.to_edn().to_string(),
        via_canonical.to_edn().to_string(),
    );
}

// The streaming parser accepts structural refs identically to the tree parser, across the three
// structural map forms (`:atoms`, `:donors`/`:acceptor`, `:site`/`:ligands`).
#[rstest]
#[case::bond(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [{:atoms [0 1]} {:aromatic true}]}]}"#)]
#[case::dative(r#"{:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [{:donors [0] :acceptor 1} {:ring-membership {:count 1}}]}]}"#)]
#[case::stereo(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [{:site 0 :ligands [1 2 3 4]} [:tetrahedral {:ligand-symmetry {:permutation [[0 1]] :orientation :improper :invariant false}}]]}]}"#)]
fn test_molecule_dsl_from_edn_str_structural_ref(#[case] input: &str) {
    let via_stream = MoleculeDsl::from_edn_str(input).unwrap();
    let via_tree = MoleculeDsl::from_edn(&read_string(input).unwrap()).unwrap();
    assert_eq!(via_stream, via_tree);
}

#[rstest]
fn test_molecule_dsl_edn_roundtrip_bond_order_sum_by_keyword() {
    let source = r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :atoms [0 1] :type :single} {:id :b2 :atoms [1 2] :type :single}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

#[rstest]
fn test_molecule_dsl_edn_roundtrip_atom_leaf_constraint_by_keyword() {
    let source =
        r##"{:atoms [[:c1 "C"]] :bonds [] :constraints [{:not {:atom [:c1 {:valence 3}]}}]}"##;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

#[rstest]
fn test_molecule_dsl_edn_roundtrip_sub_pattern() {
    let source = r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]]} :pattern {:atoms ["N"] :bonds []}}}]}"##;
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

#[rstest]
#[case::valence(r##"{:atoms ["C#v4"] :bonds []}"##)]
#[case::ring_membership_all(r##"{:atoms ["N#R2"] :bonds []}"##)]
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

#[rustfmt::skip]
#[rstest]
#[case::atom_unknown_site(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site :nope :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::atom_out_of_range_site(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 99 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::atom_unknown_ligand(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [:nope 2 3 4] :type "Th1"}]}"##)]
#[case::atom_out_of_range_ligand(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [99 2 3 4] :type "Th1"}]}"##)]
#[case::bond_unknown_site(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site :nope :ligands [0 3] :type "Ct1"}]}"##)]
#[case::bond_out_of_range_ligand(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [99 3] :type "Ct1"}]}"##)]
fn test_molecule_dsl_stereo_ref_errors(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(
        matches!(err, DeError::Custom(_)),
        "expected Custom ref error, got {:?}",
        err,
    );
}

/// `:type` is required on every entry kind that has a DSL surface (bond,
/// dative, aromatic, multicenter, noncovalent). Missing `:type` is a
/// `MissingField` error in both the streaming and tree paths.
#[rustfmt::skip]
#[rstest]
#[case::bond_without_type(r##"{:atoms ["C" "C"] :bonds [{:atoms [0 1]}]}"##)]
#[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1}]}"##)]
#[case::aromatic_without_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1]}]}"##)]
#[case::multicenter_without_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1]}]}"##)]
#[case::noncovalent_without_type(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1]}]}"##)]
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
#[case::bond_without_type(r##"{:atoms ["C" "C"] :bonds [{:atoms [0 1]}]}"##)]
#[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1}]}"##)]
#[case::aromatic_without_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1]}]}"##)]
#[case::multicenter_without_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1]}]}"##)]
#[case::noncovalent_without_type(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1]}]}"##)]
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

#[rustfmt::skip]
#[rstest]
#[case::bond_missing_atoms(r##"{:atoms ["C" "C"] :bonds [{:type "1"}]}"##)]
#[case::aromatic_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:type "*"}]}"##)]
#[case::multicenter_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:type "*"}]}"##)]
#[case::noncovalent_missing_atoms(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:type "Hbd"}]}"##)]
#[case::dative_missing_donors(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:acceptor 1 :type :single}]}"##)]
#[case::dative_missing_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :type :single}]}"##)]
#[case::stereo_atom_missing_site(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_missing_ligands(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :type "Th1"}]}"##)]
#[case::stereo_atom_missing_type(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4]}]}"##)]
#[case::stereo_bond_missing_site(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_missing_ligands(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :type "Ct1"}]}"##)]
#[case::stereo_bond_missing_type(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3]}]}"##)]
fn test_molecule_dsl_required_field_missing_tree(#[case] source: &str) {
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
#[case::bond_missing_atoms(r##"{:atoms ["C" "C"] :bonds [{:type "1"}]}"##)]
#[case::aromatic_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:type "*"}]}"##)]
#[case::multicenter_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:type "*"}]}"##)]
#[case::noncovalent_missing_atoms(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:type "Hbd"}]}"##)]
#[case::dative_missing_donors(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:acceptor 1 :type :single}]}"##)]
#[case::dative_missing_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :type :single}]}"##)]
#[case::stereo_atom_missing_site(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_missing_ligands(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :type "Th1"}]}"##)]
#[case::stereo_atom_missing_type(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4]}]}"##)]
#[case::stereo_bond_missing_site(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_missing_ligands(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :type "Ct1"}]}"##)]
#[case::stereo_bond_missing_type(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3]}]}"##)]
fn test_molecule_dsl_required_field_missing_streaming(#[case] source: &str) {
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

#[rustfmt::skip]
#[rstest]
// AtomConstraintAst variants via :constraints [{:atom [0 <form>]}]
#[case::atom_valence_lit(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence 4}]}]}"##)]
#[case::atom_valence_set(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence [3 4]}]}]}"##)]
#[case::atom_valence_undetermined(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence :undetermined}]}]}"##)]
#[case::atom_valence_expr(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:valence "?h >= 1"}]}]}"##)]
#[case::atom_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:degree 3}]}]}"##)]
#[case::atom_total_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-degree 4}]}]}"##)]
#[case::atom_ring_degree(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-degree 2}]}]}"##)]
#[case::atom_total_hydrogens(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-hydrogens 3}]}]}"##)]
#[case::atom_ring_count(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::atom_ring_size(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::atom_ring_size_set(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 6 :count [5 6]}}]}]}"##)]
#[case::atom_ring_size_conj(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-membership {:size 5 :count 1}}]} {:atom [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::atom_total_valence(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:total-valence 4}]}]}"##)]
#[case::atom_ring_valence(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:ring-valence 2}]}]}"##)]
#[case::atom_donated_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:donated-pairs 1}]}]}"##)]
#[case::atom_accepted_pairs(r##"{:atoms ["N"] :bonds [] :constraints [{:atom [0 {:accepted-pairs 2}]}]}"##)]
#[case::atom_aromatic_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence :not-aromatic}]}]}"##)]
#[case::atom_aromatic_valence_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:aromatic-valence {:aromatic 6}}]}]}"##)]
#[case::atom_multicenter_valence_not(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence :not-multicenter}]}]}"##)]
#[case::atom_multicenter_valence_value(r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:multicenter-valence {:multicenter 3}}]}]}"##)]
// BondConstraintAst variants
#[case::bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:aromatic true}]}]}"##)]
#[case::bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
// DativeBondConstraintAst variants
#[case::dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::dative_ring_size(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-membership {:size 5 :count 1}}]}]}"##)]
#[case::dative_donors(r##"{:atoms ["C" "C" "N"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type :single}] :constraints [{:dative-bond-donors [0 [0 1]]}]}"##)]
#[case::dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
#[case::dative_contains_all_donors(r##"{:atoms ["C" "C" "N"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type :single}] :constraints [{:dative-bond-contains-all-donors [0 [0 1]]}]}"##)]
#[case::dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
#[case::dative_all_donors(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-all-donors [0 {:valence 4}]}]}"##)]
#[case::dative_any_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-any-donor [0 {:degree 1}]}]}"##)]
#[case::dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor-satisfies [0 {:degree 2}]}]}"##)]
#[case::dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [0] :acceptor 1 :type :single}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
// RelationalConstraint variants for aromatic system
#[case::aromatic_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}] :constraints [{:aromatic-system-atoms [0 [0 1]]}]}"##)]
#[case::aromatic_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
#[case::aromatic_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
#[case::aromatic_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
#[case::aromatic_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*"}] :constraints [{:aromatic-system-any-atom [0 {:degree 3}]}]}"##)]
// RelationalConstraint variants for multicenter bond
#[case::multicenter_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}] :constraints [{:multicenter-bond-atoms [0 [0 1]]}]}"##)]
#[case::multicenter_contains(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
#[case::multicenter_contains_all(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
#[case::multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
#[case::multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*"}] :constraints [{:multicenter-bond-any-atom [0 {:degree 3}]}]}"##)]
// RelationalConstraint variants for noncovalent bond
#[case::noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
#[case::noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
#[case::noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 4} {:valence 1}]]}]}"##)]
// MoleculeConstraint variants (via flattened keys)
#[case::molecule_charge_sum(r##"{:atoms ["C" "N"] :bonds [] :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"##)]
#[case::molecule_unpaired_electron_coupling(r##"{:atoms ["C"] :bonds [] :constraints [{:unpaired-electron-coupling {:atoms [0] :unpaired-electrons {:count 1 :multiplicity 2}}}]}"##)]
// Anchor with multiple entity kinds (exercises all 6 ref-pair readers)
#[case::sub_pattern_anchor_bonds_and_atoms(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]] :bonds [[0 0]]} :pattern {:atoms ["N" "N"] :bonds [[0 1 "1"]]}}}]}"##)]
fn test_molecule_dsl_streaming_per_variant_parity(#[case] source: &str) {
    let via_str = MoleculeDsl::from_edn_str(source).unwrap();
    let tree = read_string(source).unwrap();
    let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
    assert_eq!(via_str, via_tree);
}

#[rstest]
#[case::missing_key_value(r##"{:atoms}"##)]
#[case::string_key(r##"{"atoms" []}"##)]
#[case::truncated_map(r##"{:atoms ["C""##)]
#[case::truncated_outer(r##"{:atoms []"##)]
#[case::unknown_top_key(r##"{:atoms [] :bonds [] :bogus 1}"##)]
#[case::atom_out_of_range_in_bond(r##"{:atoms ["C"] :bonds [[0 99 "1"]]}"##)]
#[case::positional_bond_two_vec(r##"{:atoms ["C" "C"] :bonds [[0 1]]}"##)]
#[case::positional_bond_four_vec(r##"{:atoms ["C" "C"] :bonds [[0 1 "1" 1]]}"##)]
#[case::unknown_constraint_key(r##"{:atoms ["C"] :bonds [] :constraints [{:bogus 1}]}"##)]
#[case::unknown_atom_constraint_kind(
    r##"{:atoms ["C"] :bonds [] :constraints [{:atom [0 {:bogus 1}]}]}"##
)]
#[case::dative_donors_not_vector(
    r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors 0 :acceptor 1 :type :single}]}"##
)]
#[case::dative_donors_empty(
    r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [] :acceptor 1 :type :single}]}"##
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

/// Top-level EDN maps are unordered: any permutation of the section keys parses
/// to the same `MoleculeDsl`, on both the tree and streaming paths.
#[fixture]
fn permutation_reference() -> MoleculeDsl {
    MoleculeDsl::from_edn_str(
        r##"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [2] :acceptor 0 :type :single}] :constraints [{:atom [0 {:degree 1}]}]}"##,
    )
    .unwrap()
}

#[rustfmt::skip]
#[rstest]
#[case::canonical(r##"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [2] :acceptor 0 :type :single}] :constraints [{:atom [0 {:degree 1}]}]}"##)]
#[case::reversed(r##"{:constraints [{:atom [0 {:degree 1}]}] :dative-bonds [{:donors [2] :acceptor 0 :type :single}] :bonds [[0 1 "1"]] :atoms ["C" "C" "N"]}"##)]
#[case::shuffled(r##"{:dative-bonds [{:donors [2] :acceptor 0 :type :single}] :atoms ["C" "C" "N"] :constraints [{:atom [0 {:degree 1}]}] :bonds [[0 1 "1"]]}"##)]
fn test_molecule_dsl_top_level_key_permutation(
    permutation_reference: MoleculeDsl,
    #[case] source: &str,
) {
    let via_tree = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
    let via_str = MoleculeDsl::from_edn_str(source).unwrap();
    assert_eq!(via_tree, permutation_reference);
    assert_eq!(via_str, permutation_reference);
}

#[rustfmt::skip]
#[rstest]
#[case::atom_vs_bond(r##"{:atoms [[:x "C"] [:y "C"]] :bonds [{:id :x :atoms [0 1] :type "1"}]}"##)]
#[case::atom_vs_alias(r##"{:atoms [[:x "C"]] :bonds [] :atom-aliases [:x "N"]}"##)]
#[case::bond_vs_dative(r##"{:atoms ["C" "N"] :bonds [{:id :x :atoms [0 1] :type "1"}] :dative-bonds [{:id :x :donors [0] :acceptor 1 :type :single}]}"##)]
#[case::atom_vs_aromatic(r##"{:atoms [[:x "C"] [:y "C"]] :bonds [] :aromatic-systems [{:id :x :atoms [0 1] :type "*"}]}"##)]
#[case::bond_vs_noncovalent(r##"{:atoms ["C" "C"] :bonds [{:id :x :atoms [0 1] :type "1"}] :noncovalent-bonds [{:id :x :atoms [0 1] :type "Hbd"}]}"##)]
fn test_molecule_dsl_from_edn_cross_entity_keyword_collision_error(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)));
}

/// `from_edn` reference, alias, and keyword failures all surface as `DeError::Custom`.
#[rustfmt::skip]
#[rstest]
#[case::unknown_alias(r##"{:atoms [:x] :bonds []}"##)]
#[case::duplicate_atom_keyword(r##"{:atoms [[:a "C"] [:a "N"]] :bonds []}"##)]
#[case::duplicate_alias_name(r##"{:atoms [] :bonds [] :atom-aliases [:a "C" :a "N"]}"##)]
#[case::aliases_not_bijective(r##"{:atoms [] :bonds [] :atom-aliases [:a "C" :b "C"]}"##)]
#[case::bond_endpoint_out_of_range(r##"{:atoms ["C" "C"] :bonds [[0 5 "1"]]}"##)]
#[case::bond_endpoint_unknown_keyword(r##"{:atoms ["C" "C"] :bonds [[:nope 0 "1"]]}"##)]
#[case::aromatic_atom_out_of_range(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 5] :type "*"}]}"##)]
#[case::noncovalent_endpoint_out_of_range(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 99] :type "Hbd"}]}"##)]
#[case::dative_unknown_donor_keyword(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [:nope] :acceptor 1 :type :single}]}"##)]
#[case::constraint_unknown_ref(r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [:nope 0]}}]}"##)]
#[case::sub_pattern_ref_out_of_range(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:atoms [[0 5]]} :pattern {:atoms ["N"] :bonds []}}}]}"##)]
fn test_molecule_dsl_from_edn_error(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)), "expected Custom, got {:?}", err);
}

/// `MoleculeAst::to_edn` emits canonical EDN with positional refs only,
/// regardless of any entity keywords on the input. Parsing the canonical
/// output back yields the same AST.
#[rstest]
fn test_molecule_ast_to_edn_canonical_positional_refs() {
    // Input has entity keywords on atoms, bonds, and a constraint anchor.
    let source = r##"{:atoms [[:c1 "C"] [:c2 "C"]]
                      :bonds [{:id :b1 :atoms [:c1 :c2] :type "1"}]
                      :constraints [{:atom [:c1 {:valence 4}]}
                                    {:bond [:b1 {:aromatic true}]}]}"##;
    let dsl = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
    let (ast, _meta) = dsl.into_parts();

    // Canonical render: positional refs only.
    let canonical_source = r##"{:atoms ["C" "C"] :bonds [[0 1 :single]]
             :constraints [{:atom [0 {:valence 4}]}
                           {:bond [0 {:aromatic true}]}]}"##;
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
    assert_eq!(
        ast,
        MoleculeAst::from_edn(&read_string(source).unwrap()).unwrap()
    );
}

#[rstest]
fn test_molecule_ast_from_edn_keyword_metadata() {
    // Input carries entity keywords; AST is metadata-free, so reparsing the rendered
    // form (which has no keywords) should match the AST from the original parse.
    let source = r##"{:atoms [[:carbon "C"] [:oxygen "O"]]
                      :bonds [{:id :myb :atoms [:carbon :oxygen] :type "1"}]}"##;
    let ast = MoleculeAst::from_edn(&read_string(source).unwrap()).unwrap();
    let rendered = ast.to_edn().to_string();
    // No user-defined entity keywords leaked through. (`:a` / `:b` / `:type`
    // are bond-entry field names, not entity keywords.)
    assert!(!rendered.contains(":carbon"));
    assert!(!rendered.contains(":oxygen"));
    assert!(!rendered.contains(":myb"));
    let reparsed = MoleculeAst::from_edn_str(&rendered).unwrap();
    assert_eq!(ast, reparsed);
}

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

// Atom value + the atom / bond entry renderers (analogs of the overlay `render_<entity>_entry`):
// the atom value is its alias keyword or atom-string; the atom entry adds the `[keyword …]` framing;
// the bond entry places the caller's `type_edn` under `:type` (or in the `[a b type]` vector,
// without an `:id` field).

#[rstest]
#[case::no_alias(false, r#""C""#)]
#[case::alias(true, r#":x"#)]
fn test_render_atom_value(#[case] alias: bool, #[case] expected: &str) {
    let atom = AtomAst::from_element(Element::C);
    let mut meta = MoleculeMetadata::new();
    if alias {
        meta.add_atom_alias("x", atom.clone()).unwrap();
    }
    assert_eq!(
        render_atom_value(&atom, &meta),
        read_string(expected).unwrap()
    );
}

#[rstest]
#[case::no_keyword(None, r#""C""#)]
#[case::with_keyword(Some("c0"), r#"[:c0 "C"]"#)]
fn test_render_atom_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::Atom(AtomId(0)), keyword).unwrap();
    }
    let entry = render_atom_entry(AtomId(0), &AtomAst::from_element(Element::C), &meta);
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::positional(MoleculeMetadata::new(), AtomId(2), "2")]
#[case::keyword(
    {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(2)), "carbon").unwrap();
        metadata
    },
    AtomId(2),
    ":carbon",
)]
fn test_render_atom_ref(
    #[case] metadata: MoleculeMetadata,
    #[case] id: AtomId,
    #[case] expected: &str,
) {
    assert_eq!(
        render_atom_ref(id, &metadata),
        read_string(expected).unwrap()
    );
}

#[rstest]
#[case::no_id(None, r#"[0 1 "1"]"#)]
#[case::with_id(Some("b0"), r#"{:id :b0 :atoms [0 1] :type "1"}"#)]
fn test_render_bond_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::Bond(BondId(0)), keyword).unwrap();
    }
    let entry = render_bond_entry(BondId(0), [AtomId(0), AtomId(1)], Edn::string("1"), &meta);
    assert_eq!(entry, read_string(expected).unwrap());
}

// I4a — the overlay `render_<entity>_entry` renderers: `:id` present iff the metadata binds one,
// participants rendered as refs (positional without keywords), and the caller's rendered `:type` placed
// under `:type` verbatim (realistic values: a bond/dative order `"1"`, an aromatic/multicenter
// electron string `"[1,1,0]"`, a noncovalent `"Hbd"`, a stereo `"Th0"` / `"Ct0"`).

#[rstest]
#[case::no_id(None, r#"{:donors [0 2] :acceptor 1 :type "1"}"#)]
#[case::with_id(Some("d0"), r#"{:id :d0 :donors [0 2] :acceptor 1 :type "1"}"#)]
fn test_render_dative_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::DativeBond(DativeBondId(0)), keyword)
            .unwrap();
    }
    let entry = render_dative_entry(
        DativeBondId(0),
        [AtomId(0), AtomId(2)].into_iter(),
        AtomId(1),
        Edn::string("1"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::no_id(None, r#"{:atoms [0 1 2] :type "[1,1,0]"}"#)]
#[case::with_id(Some("r0"), r#"{:id :r0 :atoms [0 1 2] :type "[1,1,0]"}"#)]
fn test_render_aromatic_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::AromaticSystem(AromaticSystemId(0)), keyword)
            .unwrap();
    }
    let entry = render_aromatic_entry(
        AromaticSystemId(0),
        [AtomId(0), AtomId(1), AtomId(2)].into_iter(),
        Edn::string("[1,1,0]"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::no_id(None, r#"{:atoms [0 1 2] :type "[1,1,0]"}"#)]
#[case::with_id(Some("m0"), r#"{:id :m0 :atoms [0 1 2] :type "[1,1,0]"}"#)]
fn test_render_multicenter_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::MulticenterBond(MulticenterBondId(0)), keyword)
            .unwrap();
    }
    let entry = render_multicenter_entry(
        MulticenterBondId(0),
        [AtomId(0), AtomId(1), AtomId(2)].into_iter(),
        Edn::string("[1,1,0]"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::no_id(None, r#"{:atoms [0 1] :type "Hbd"}"#)]
#[case::with_id(Some("n0"), r#"{:id :n0 :atoms [0 1] :type "Hbd"}"#)]
fn test_render_noncovalent_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::NoncovalentBond(NoncovalentBondId(0)), keyword)
            .unwrap();
    }
    let entry = render_noncovalent_entry(
        NoncovalentBondId(0),
        [AtomId(0), AtomId(1)],
        Edn::string("Hbd"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::no_id(None, r#"{:site 0 :ligands [1 2] :type "Th0"}"#)]
#[case::with_id(Some("s0"), r#"{:id :s0 :site 0 :ligands [1 2] :type "Th0"}"#)]
fn test_render_stereo_atom_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::StereoAtom(StereoAtomId(0)), keyword)
            .unwrap();
    }
    let entry = render_stereo_atom_entry(
        StereoAtomId(0),
        AtomId(0),
        vec![Edn::Int(1), Edn::Int(2)],
        Edn::string("Th0"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::no_id(None, r#"{:site 0 :ligands [1] :type "Ct0"}"#)]
#[case::with_id(Some("s0"), r#"{:id :s0 :site 0 :ligands [1] :type "Ct0"}"#)]
fn test_render_stereo_bond_entry(#[case] keyword: Option<&str>, #[case] expected: &str) {
    let mut meta = MoleculeMetadata::new();
    if let Some(keyword) = keyword {
        meta.set_keyword(Entity::StereoBond(StereoBondId(0)), keyword)
            .unwrap();
    }
    let entry = render_stereo_bond_entry(
        StereoBondId(0),
        BondId(0),
        vec![Edn::Int(1)],
        Edn::string("Ct0"),
        &meta,
    );
    assert_eq!(entry, read_string(expected).unwrap());
}

#[rstest]
#[case::atom_positional(
    MoleculeMetadata::new(),
    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
    "2"
)]
#[case::atom_keyword(
    {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(2)), "carbon").unwrap();
        metadata
    },
    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
    ":carbon",
)]
#[case::implicit_hydrogen_positional(
    MoleculeMetadata::new(),
    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
    "[:h 2]"
)]
#[case::implicit_hydrogen_keyword(
    {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(2)), "carbon").unwrap();
        metadata
    },
    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
    "[:h :carbon]",
)]
#[case::lone_pair_positional(
    MoleculeMetadata::new(),
    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    "[:lp 2]"
)]
#[case::lone_pair_keyword(
    {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Atom(AtomId(2)), "carbon").unwrap();
        metadata
    },
    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    "[:lp :carbon]",
)]
fn test_render_stereo_ligand(
    #[case] metadata: MoleculeMetadata,
    #[case] ligand: StereoLigand,
    #[case] expected: &str,
) {
    assert_eq!(
        render_stereo_ligand(ligand, &metadata),
        read_string(expected).unwrap()
    );
}

#[rstest]
#[case::positional(MoleculeMetadata::new(), BondId(2), "2")]
#[case::keyword(
    {
        let mut metadata = MoleculeMetadata::new();
        metadata.set_keyword(Entity::Bond(BondId(2)), "bond").unwrap();
        metadata
    },
    BondId(2),
    ":bond",
)]
fn test_render_bond_ref(
    #[case] metadata: MoleculeMetadata,
    #[case] id: BondId,
    #[case] expected: &str,
) {
    assert_eq!(
        render_bond_ref(id, &metadata),
        read_string(expected).unwrap()
    );
}
