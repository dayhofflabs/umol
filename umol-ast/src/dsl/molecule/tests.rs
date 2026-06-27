use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_edn::read_string;

use super::*;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::{Constraint, Constraints, MoleculeConstraint};
use crate::ast::electrons::ElectronCountsAst;
use crate::ast::spin::SpinStateAst;
use crate::ast::value::ValueAst;
use crate::{dsl, mol};

#[rstest]
fn test_metadata_new() {
    let m = MoleculeMetadata::new();
    assert_eq!(m, MoleculeMetadata::default());
    assert!(m.atom_id(AtomId(0)).is_none());
    assert!(m.bond_id(BondId(0)).is_none());
    assert!(m.dative_bond_id(DativeBondId(0)).is_none());
    assert!(m.aromatic_system_id(AromaticSystemId(0)).is_none());
    assert!(m.multicenter_bond_id(MulticenterBondId(0)).is_none());
    assert!(m.noncovalent_bond_id(NoncovalentBondId(0)).is_none());
    assert!(!m.has_atom_aliases());
    assert_eq!(m.atom_aliases_len(), 0);
}

#[rstest]
#[case::set(&["c1"], "c1")]
#[case::last_wins(&["old", "new"], "new")]
fn test_metadata_set_atom_id(#[case] names: &[&str], #[case] expected: &str) {
    let mut m = MoleculeMetadata::new();
    for name in names {
        m.set_atom_id(AtomId(0), *name);
    }
    assert_eq!(m.atom_id(AtomId(0)), Some(expected));
}

#[rstest]
fn test_metadata_set_bond_id() {
    let mut m = MoleculeMetadata::new();
    m.set_bond_id(BondId(2), "b1");
    assert_eq!(m.bond_id(BondId(2)), Some("b1"));
}

#[rstest]
fn test_metadata_set_dative_bond_id() {
    let mut m = MoleculeMetadata::new();
    m.set_dative_bond_id(DativeBondId(1), "d1");
    assert_eq!(m.dative_bond_id(DativeBondId(1)), Some("d1"));
}

#[rstest]
fn test_metadata_set_aromatic_system_id() {
    let mut m = MoleculeMetadata::new();
    m.set_aromatic_system_id(AromaticSystemId(0), "ring1");
    assert_eq!(m.aromatic_system_id(AromaticSystemId(0)), Some("ring1"));
}

#[rstest]
fn test_metadata_set_multicenter_bond_id() {
    let mut m = MoleculeMetadata::new();
    m.set_multicenter_bond_id(MulticenterBondId(0), "mc1");
    assert_eq!(m.multicenter_bond_id(MulticenterBondId(0)), Some("mc1"));
}

#[rstest]
fn test_metadata_set_noncovalent_bond_id() {
    let mut m = MoleculeMetadata::new();
    m.set_noncovalent_bond_id(NoncovalentBondId(0), "h1");
    assert_eq!(m.noncovalent_bond_id(NoncovalentBondId(0)), Some("h1"));
}

#[rstest]
fn test_metadata_add_atom_alias() {
    let mut m = MoleculeMetadata::new();
    let atom = AtomAst::from_element(Element::C).with_implicit_hydrogens(2_i64);
    m.add_atom_alias("HC2", atom.clone());
    assert!(m.has_atom_alias("HC2"));
    assert_eq!(m.atom_aliases_len(), 1);
    assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("HC2"));
}

#[rstest]
fn test_metadata_add_atom_alias_duplicate_name_replaces_atom() {
    let mut m = MoleculeMetadata::new();
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
    let mut m = MoleculeMetadata::new();
    let atom = AtomAst::from_element(Element::C);
    m.add_atom_alias("first", atom.clone());
    m.add_atom_alias("second", atom.clone());
    assert_eq!(m.atom_aliases_len(), 1);
    assert!(!m.has_atom_alias("first"));
    assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("second"));
}

#[rstest]
fn test_metadata_iter_atom_aliases() {
    let m = MoleculeMetadata::new()
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
    let m = MoleculeMetadata::new()
        .with_atom_id(AtomId(0), "a")
        .with_atom_id(AtomId(1), "b");
    assert_eq!(m.atom_id(AtomId(0)), Some("a"));
    assert_eq!(m.atom_id(AtomId(1)), Some("b"));
}

#[rstest]
fn test_metadata_with_bond_id() {
    let m = MoleculeMetadata::new().with_bond_id(BondId(0), "b");
    assert_eq!(m.bond_id(BondId(0)), Some("b"));
}

#[rstest]
fn test_metadata_with_dative_bond_id() {
    let m = MoleculeMetadata::new().with_dative_bond_id(DativeBondId(0), "d");
    assert_eq!(m.dative_bond_id(DativeBondId(0)), Some("d"));
}

#[rstest]
fn test_metadata_with_aromatic_system_id() {
    let m = MoleculeMetadata::new().with_aromatic_system_id(AromaticSystemId(0), "r");
    assert_eq!(m.aromatic_system_id(AromaticSystemId(0)), Some("r"));
}

#[rstest]
fn test_metadata_with_multicenter_bond_id() {
    let m = MoleculeMetadata::new().with_multicenter_bond_id(MulticenterBondId(0), "mc");
    assert_eq!(m.multicenter_bond_id(MulticenterBondId(0)), Some("mc"));
}

#[rstest]
fn test_metadata_with_noncovalent_bond_id() {
    let m = MoleculeMetadata::new().with_noncovalent_bond_id(NoncovalentBondId(0), "h");
    assert_eq!(m.noncovalent_bond_id(NoncovalentBondId(0)), Some("h"));
}

#[rstest]
fn test_metadata_with_stereo_atom_id() {
    let m = MoleculeMetadata::new().with_stereo_atom_id(StereoAtomId(0), "s");
    assert_eq!(m.stereo_atom_id(StereoAtomId(0)), Some("s"));
}

#[rstest]
fn test_metadata_with_stereo_bond_id() {
    let m = MoleculeMetadata::new().with_stereo_bond_id(StereoBondId(0), "sb");
    assert_eq!(m.stereo_bond_id(StereoBondId(0)), Some("sb"));
}

#[rstest]
fn test_metadata_with_atom_alias() {
    let atom = AtomAst::from_element(Element::C);
    let m = MoleculeMetadata::new().with_atom_alias("c", atom.clone());
    assert_eq!(m.atom_alias_for(&AtomDsl(atom)), Some("c"));
}

#[rstest]
fn test_metadata_mixed_chain() {
    let m = MoleculeMetadata::new()
        .with_atom_id(AtomId(0), "c1")
        .with_bond_id(BondId(0), "b1")
        .with_atom_alias("X", AtomAst::from_element(Element::C));
    assert_eq!(m.atom_id(AtomId(0)), Some("c1"));
    assert_eq!(m.bond_id(BondId(0)), Some("b1"));
    assert!(m.has_atom_alias("X"));
}

#[rstest]
fn test_molecule_dsl_to_edn_empty() {
    let ast = MoleculeAst::default();
    let dsl = MoleculeDsl::from_parts(ast, MoleculeMetadata::default());
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
    let dsl = dsl!(r#"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]}"#);
    let edn = dsl.to_edn();
    assert_eq!(
        edn,
        read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type :single}]}"##)
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
    assert_eq!(dsl.metadata().atom_id(AtomId(0)), Some("c1"));
    assert_eq!(dsl.metadata().atom_id(AtomId(1)), None);
}

#[rstest]
fn test_molecule_dsl_from_edn_bond_map_form_with_id() {
    let edn =
        read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]}"##).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.ast().bonds().count(), 1);
    assert_eq!(dsl.metadata().bond_id(BondId(0)), Some("b1"));
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
#[case::with_ids(r##"{:atoms [[:a "C"] [:b "N"]] :bonds [{:id :b1 :atoms [:a :b] :type "1"}]}"##)]
#[case::inline_atom_constraints(r##"{:atoms ["C#v4" "N#R+"] :bonds []}"##)]
#[case::inline_bond_constraint(r##"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"##)]
#[case::aromatic_section(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:id :ar1 :atoms [0 1 2 3 4 5] :type "#e6"}]}"##)]
#[case::multicenter_section(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "#e2"}]}"##)]
#[case::dative_section(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:id :d1 :donor 0 :acceptor 1 :type "1#R"}]}"##)]
#[case::dative_multi_donor(r##"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donor [0 1] :acceptor 2 :type "1#R"}]}"##)]
#[case::noncovalent_section(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::stereo_atom_section(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_id_and_keyword(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :type :ccw}]}"##)]
#[case::stereo_atom_virtual_ligand(r##"{:atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 [:lp 0] [:h 0]] :type "Th1"}]}"##)]
#[case::stereo_bond_section(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_id_ref(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:id :db :atoms [1 2] :type "2"} [2 3 "1"]] :stereo-bonds [{:site :db :ligands [0 3] :type :e}]}"##)]
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
#[case::constraints_bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 :aromatic]}]}"##)]
#[case::constraints_bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::constraints_bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
#[case::constraints_dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::constraints_dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
#[case::constraints_dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
#[case::constraints_dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
#[case::constraints_dative_donor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-donor-satisfies [0 {:valence 3}]}]}"##)]
#[case::constraints_dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond-acceptor-satisfies [0 {:valence 3}]}]}"##)]
#[case::constraints_aromatic_system_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
#[case::constraints_aromatic_system_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
#[case::constraints_aromatic_system_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
#[case::constraints_aromatic_system_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}] :constraints [{:aromatic-system-any-atom [0 {:valence 4}]}]}"##)]
#[case::constraints_multicenter_contains(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "#e3"}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
#[case::constraints_multicenter_contains_all(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "#e3"}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
#[case::constraints_multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "#e2"}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
#[case::constraints_multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "#e2"}] :constraints [{:multicenter-bond-any-atom [0 {:valence 4}]}]}"##)]
#[case::constraints_noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
#[case::constraints_noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
#[case::constraints_noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 3} {:valence 1}]]}]}"##)]
#[case::constraints_sub_pattern_multi_entity_anchor(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :constraints [{:sub-pattern {:anchor {:atoms [[0 0]] :bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [[0 1 "1"]]}}}]}"##)]
#[case::constraints_sub_pattern_dative_anchor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:sub-pattern {:anchor {:dative-bonds [[0 0]]} :pattern {:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}]}}}]}"##)]
#[case::constraints_sub_pattern_aromatic_system_anchor(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}] :constraints [{:sub-pattern {:anchor {:aromatic-systems [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "#e2"}]}}}]}"##)]
#[case::constraints_sub_pattern_multicenter_anchor(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "#e2"}] :constraints [{:sub-pattern {:anchor {:multicenter-bonds [[0 0]]} :pattern {:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "#e2"}]}}}]}"##)]
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
#[case::stereo_atom_molecule_constraint(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [0 [:tetrahedral {:ligand-symmetry {:permutation [[0 1]] :orientation :improper :member :not-in}}]]}]}"##)]
#[case::stereo_atom_fluxionality_constraint(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}] :constraints [{:stereo-atom [0 [:tetrahedral {:fluxionality [[0 1 2]]}]]}]}"##)]
#[case::stereo_bond_molecule_constraint(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}] :constraints [{:stereo-bond [0 [:cis-trans {:topicity {:pair [0 1] :relation :diastereotopic}}]]}]}"##)]
fn test_molecule_dsl_from_edn_str_matches_from_edn(#[case] source: &str) {
    let via_str = MoleculeDsl::from_edn_str(source).unwrap();
    let tree = read_string(source).unwrap();
    let via_tree = MoleculeDsl::from_edn(&tree).unwrap();
    assert_eq!(via_str, via_tree);
}

#[rustfmt::skip]
#[rstest]
#[case::stereo_atom(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_atom_id_virtual(r##"{:atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 [:lp 0] [:h 0]] :type "Th1"}]}"##)]
#[case::stereo_bond_id_ref(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] {:id :db :atoms [1 2] :type "2"} [2 3 "1"]] :stereo-bonds [{:site :db :ligands [0 3] :type "Ct1"}]}"##)]
#[case::stereo_bond_own_id(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :sb1 :site 1 :ligands [0 3] :type "Ct1"}]}"##)]
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
#[case::dative_constraint_unknown_key(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1#R"}] :constraints [{:dative-bond [0 {:bogus 1}]}]}"##)]
#[case::sub_pattern_anchor_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:sub-pattern {:anchor {:bogus [[0 0]]} :pattern {:atoms ["C"] :bonds []}}}]}"##)]
#[case::constraint_unknown_key(r##"{:atoms ["C"] :bonds [] :constraints [{:bogus 1}]}"##)]
#[case::noncovalent_ends_satisfy_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 2}]]}]}"##)]
#[case::noncovalent_ends_wrong_pair_length(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0]]}]}"##)]
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

// Round-trip direction: DSL → AST (raise) → DSL (lower) is the
// identity. AST → DSL → AST isn't, since raising `Undetermined`
// fields to `Lit(0)` is one-way under `zeroed()`. One case per overlay
// kind so the `into_ast` / `from_ast` per-relation loops are exercised.
#[rustfmt::skip]
#[rstest]
#[case::atoms_bonds(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##)]
#[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}]}"##)]
#[case::aromatic(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}]}"##)]
#[case::multicenter(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type ""}]}"##)]
#[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::stereo_atom(r##"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"##)]
#[case::stereo_bond(r##"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"##)]
fn test_molecule_dsl_dsl_to_ast_to_dsl_roundtrip_zeroed(#[case] source: &str) {
    let ast = mol!(source);
    let dsl = MoleculeDsl::from_parts(ast, MoleculeMetadata::default());
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
    assert_eq!(dsl.metadata(), &MoleculeMetadata::default());
}

#[rustfmt::skip]
#[rstest]
#[case::dative(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}]}"##)]
#[case::dative_with_id_and_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:id :d1 :donor 0 :acceptor 1 :type "1#R"}]}"##)]
#[case::dative_multi_donor(r##"{:atoms ["C" "C" "C"] :bonds [] :dative-bonds [{:donor [0 1] :acceptor 2 :type :single}]}"##)]
#[case::aromatic_minimal(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type ""}]}"##)]
#[case::aromatic_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:id :a1 :atoms [0 1] :type "#e6"}]}"##)]
#[case::aromatic_with_electrons_literals(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :electrons [1 1 1 1 1 1] :type ""}]}"##)]
#[case::aromatic_with_electrons_and_total(r##"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3 4 5] :electrons [1 1 1 1 1 1] :type "#e6"}]}"##)]
#[case::multicenter_minimal(r##"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type ""}]}"##)]
#[case::multicenter_with_id_and_type(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:id :m1 :atoms [0 1] :type "#e2"}]}"##)]
#[case::multicenter_with_electrons_literals(r##"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :electrons [1 0 1] :type ""}]}"##)]
#[case::noncovalent(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"##)]
#[case::noncovalent_with_id(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:id :n1 :atoms [0 1] :type "Hbd"}]}"##)]
fn test_molecule_dsl_edn_roundtrip_non_localized_entities(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let dsl = MoleculeDsl::from_edn(&edn).unwrap();
    assert_eq!(dsl.to_edn(), edn);
}

/// Explicit `:electrons :undetermined` parses to `Undetermined` (the same
/// state as an omitted key); it is not a stable-roundtrip form because
/// rendering `Undetermined` omits the key.
#[rstest]
#[case::aromatic(
    r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :electrons :undetermined :type ""}]}"##
)]
#[case::multicenter(
    r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :electrons :undetermined :type ""}]}"##
)]
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
    assert_eq!(dsl.ast().constraints().len(), 1);
}

/// Vacuous molecule-level constraints — `ChargeSum` / `BondOrderSum`
/// with `Undetermined` sum, `SpinSum` with both spin fields
/// `Undetermined` — are dropped during AST → DSL lowering. The
/// canonical EDN omits the entire `:constraints` key when the only
/// entries are vacuous.
#[rstest]
#[case::charge_sum(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Undetermined })]
#[case::bond_order_sum(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Undetermined })]
#[case::spin_sum(MoleculeConstraint::SpinSum { atoms: None, spin: SpinStateAst::default() })]
fn test_molecule_dsl_render_elides_vacuous_molecule_constraint(
    #[case] constraint: MoleculeConstraint,
) {
    let mut ast = MoleculeAst::from_parts(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        vec![],
        vec![],
        vec![],
        vec![],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    );
    ast.constraints_mut().push(Constraint::Molecule(constraint));
    let dsl = MoleculeDsl::from_parts(ast, MoleculeMetadata::default());
    let Edn::Map(m) = &dsl.to_edn() else {
        panic!("expected map")
    };
    assert!(
        m.get_keyword("constraints").is_none(),
        "vacuous molecule constraint should not surface as :constraints",
    );
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
        Vec::new(),
        Vec::new(),
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
    let dsl = MoleculeDsl::from_parts(ast, MoleculeMetadata::default());
    let edn = dsl.to_edn();
    let Edn::Map(m) = &edn else {
        panic!("expected map")
    };
    let cs = m
        .get_keyword("constraints")
        .expect("constraints key present");
    let Edn::Vector(v) = cs else {
        panic!("expected vec")
    };
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
    let source = r##"{:atoms ["C" "C" "C"] :bonds [{:id :b1 :atoms [0 1] :type :single} {:id :b2 :atoms [1 2] :type :single}] :constraints [{:bond-order-sum {:bonds [:b1 :b2] :sum 2}}]}"##;
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
    let source = r##"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [:nope 0]}}]}"##;
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
        r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 99] :type "Hbd"}]}"##,
    )
    .unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)));
}

#[rstest]
fn test_molecule_dsl_aromatic_atom_out_of_range_errors() {
    let edn = read_string(
        r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 5] :type ""}]}"##,
    )
    .unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)));
}

#[rstest]
fn test_molecule_dsl_dative_unknown_donor_id_errors() {
    let edn = read_string(
        r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor :nope :acceptor 1 :type :single}]}"##,
    )
    .unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)));
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
#[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1}]}"##)]
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
#[case::dative_without_type(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1}]}"##)]
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

/// Every non-`:type` required field on an entry kind: `:atoms` (bond, aromatic,
/// multicenter, noncovalent), `:donor` / `:acceptor` (dative), and
/// `:site` / `:ligands` / `:type` (stereo atom and stereo bond). Each omission
/// is a `MissingField` error.
#[rustfmt::skip]
#[rstest]
#[case::bond_missing_atoms(r##"{:atoms ["C" "C"] :bonds [{:type "1"}]}"##)]
#[case::aromatic_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:type ""}]}"##)]
#[case::multicenter_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:type ""}]}"##)]
#[case::noncovalent_missing_atoms(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:type "Hbd"}]}"##)]
#[case::dative_missing_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:acceptor 1 :type :single}]}"##)]
#[case::dative_missing_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :type :single}]}"##)]
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
#[case::aromatic_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:type ""}]}"##)]
#[case::multicenter_missing_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:type ""}]}"##)]
#[case::noncovalent_missing_atoms(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:type "Hbd"}]}"##)]
#[case::dative_missing_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:acceptor 1 :type :single}]}"##)]
#[case::dative_missing_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :type :single}]}"##)]
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
// AtomConstraint variants via :constraints [{:atom [0 <form>]}]
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
// BondConstraint variants
#[case::bond_aromatic(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 :aromatic]}]}"##)]
#[case::bond_ring_count(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::bond_ring_size(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :constraints [{:bond [0 {:ring-membership {:size 6 :count 1}}]}]}"##)]
// DativeBondConstraint variants
#[case::dative_ring_count(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-membership {:count 1}}]}]}"##)]
#[case::dative_ring_size(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond [0 {:ring-membership {:size 5 :count 1}}]}]}"##)]
#[case::dative_donor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-donor [0 0]}]}"##)]
#[case::dative_acceptor(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor [0 1]}]}"##)]
#[case::dative_donor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-donor-satisfies [0 {:valence 4}]}]}"##)]
#[case::dative_acceptor_satisfies(r##"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-acceptor-satisfies [0 {:degree 2}]}]}"##)]
#[case::dative_parallels(r##"{:atoms ["C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donor 0 :acceptor 1 :type :single}] :constraints [{:dative-bond-parallels [0 0]}]}"##)]
// RelationalConstraint variants for aromatic system
#[case::aromatic_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-atoms [0 [0 1]]}]}"##)]
#[case::aromatic_contains(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-contains [0 0]}]}"##)]
#[case::aromatic_contains_all(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-contains-all [0 [0 1]]}]}"##)]
#[case::aromatic_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-all-atoms [0 {:valence 4}]}]}"##)]
#[case::aromatic_any_atom(r##"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type ""}] :constraints [{:aromatic-system-any-atom [0 {:degree 3}]}]}"##)]
// RelationalConstraint variants for multicenter bond
#[case::multicenter_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-atoms [0 [0 1]]}]}"##)]
#[case::multicenter_contains(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-contains [0 0]}]}"##)]
#[case::multicenter_contains_all(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-contains-all [0 [0 1]]}]}"##)]
#[case::multicenter_all_atoms(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-all-atoms [0 {:valence 4}]}]}"##)]
#[case::multicenter_any_atom(r##"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type ""}] :constraints [{:multicenter-bond-any-atom [0 {:degree 3}]}]}"##)]
// RelationalConstraint variants for noncovalent bond
#[case::noncovalent_ends(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends [0 [0 1]]}]}"##)]
#[case::noncovalent_contains(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-contains [0 0]}]}"##)]
#[case::noncovalent_ends_satisfy(r##"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}] :constraints [{:noncovalent-bond-ends-satisfy [0 [{:valence 4} {:valence 1}]]}]}"##)]
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
        r##"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donor 2 :acceptor 0 :type :single}] :constraints [{:atom [0 {:degree 1}]}]}"##,
    )
    .unwrap()
}

#[rustfmt::skip]
#[rstest]
#[case::canonical(r##"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"]] :dative-bonds [{:donor 2 :acceptor 0 :type :single}] :constraints [{:atom [0 {:degree 1}]}]}"##)]
#[case::reversed(r##"{:constraints [{:atom [0 {:degree 1}]}] :dative-bonds [{:donor 2 :acceptor 0 :type :single}] :bonds [[0 1 "1"]] :atoms ["C" "C" "N"]}"##)]
#[case::shuffled(r##"{:dative-bonds [{:donor 2 :acceptor 0 :type :single}] :atoms ["C" "C" "N"] :constraints [{:atom [0 {:degree 1}]}] :bonds [[0 1 "1"]]}"##)]
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
#[case::bond_vs_dative(r##"{:atoms ["C" "N"] :bonds [{:id :x :atoms [0 1] :type "1"}] :dative-bonds [{:id :x :donor 0 :acceptor 1 :type :single}]}"##)]
#[case::atom_vs_aromatic(r##"{:atoms [[:x "C"] [:y "C"]] :bonds [] :aromatic-systems [{:id :x :atoms [0 1] :type ""}]}"##)]
#[case::bond_vs_noncovalent(r##"{:atoms ["C" "C"] :bonds [{:id :x :atoms [0 1] :type "1"}] :noncovalent-bonds [{:id :x :atoms [0 1] :type "Hbd"}]}"##)]
fn test_molecule_dsl_cross_entity_id_collision_errors(#[case] source: &str) {
    let edn = read_string(source).unwrap();
    let err = MoleculeDsl::from_edn(&edn).unwrap_err();
    assert!(matches!(err, DeError::Custom(_)));
}

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

/// `MoleculeAst::to_edn` emits canonical EDN with positional refs only,
/// regardless of any id keywords on the input. Parsing the canonical
/// output back yields the same AST.
#[rstest]
fn test_molecule_ast_to_edn_canonical_positional_refs() {
    // Input has id keywords on atoms, bonds, and a constraint anchor.
    let source = r##"{:atoms [[:c1 "C"] [:c2 "C"]]
                      :bonds [{:id :b1 :atoms [:c1 :c2] :type "1"}]
                      :constraints [{:atom [:c1 {:valence 4}]}
                                    {:bond [:b1 :aromatic]}]}"##;
    let dsl = MoleculeDsl::from_edn(&read_string(source).unwrap()).unwrap();
    let (ast, _meta) = dsl.into_parts();

    // Canonical render: positional refs only.
    let canonical_source = r##"{:atoms ["C" "C"] :bonds [[0 1 :single]]
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
                      :bonds [{:id :myb :atoms [:carbon :oxygen] :type "1"}]}"##;
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
