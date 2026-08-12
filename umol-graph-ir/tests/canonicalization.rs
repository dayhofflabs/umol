//! Exact compatibility cases for aggregate canonicalization.
//!
//! The CDK and RDKit cases preserve published source-test partitions. The
//! InChI case preserves the molecule used by its permutation utility test.
//! External programs are not test dependencies: the expected partitions and
//! explicitly renumbered inputs are checked in here. The remaining cases cover
//! the graph-IR-specific dative, aromatic, multicenter, noncovalent, stereo-atom,
//! and stereo-bond entities.

#[path = "canonicalization/reaction_span.rs"]
mod reaction_span;

use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomId, ConstitutionColoring, GraphSymmetry, GraphSymmetryConfig, Molecule,
    MoleculeColoringFeatures,
};
use umol_graph_ir::mol_dsl;

// CDK's AtomDiscretePartitionRefinerTest applies aromaticity to
// TestMoleculeFactory.makeBiphenyl() and records
// `0,6|1,5,7,11|2,4,8,10|3,9`. The graph-IR translation therefore uses
// aromatic systems rather than treating the factory's Kekule bond orders as
// fixed distinguishing attributes.
const CDK_BIPHENYL: &str = r#"{
    :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
    :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
            [0 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"]
            [11 6 "1"]]
    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#e6"}
                       {:atoms [6 7 8 9 10 11] :attrs "*#e6"}]
}"#;
const CDK_BIPHENYL_RENUMBERED: &str = r#"{
    :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
    :bonds [[6 11 "1"] [11 10 "1"] [10 9 "1"] [9 8 "1"] [8 7 "1"] [7 6 "1"]
            [6 0 "1"] [0 5 "1"] [5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"]
            [1 0 "1"]]
    :aromatic-systems [{:atoms [6 11 10 9 8 7] :attrs "*#e6"}
                       {:atoms [0 5 4 3 2 1] :attrs "*#e6"}]
}"#;
const CDK_BIPHENYL_MAPPING: &[u32] = &[6, 11, 10, 9, 8, 7, 0, 5, 4, 3, 2, 1];
const CDK_BIPHENYL_ORBITS: &[&[u32]] = &[&[0, 6], &[1, 5, 7, 11], &[2, 4, 8, 10], &[3, 9]];

// RDKit's `ensure unused features are not used/isotopes` case ranks the two
// terminal carbons of [13CH3]OC separately when isotopes are included and
// together when they are excluded.
const RDKIT_ISOTOPE: &str = r#"{:atoms ["C#i13" "O" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#;
const RDKIT_ISOTOPE_RENUMBERED: &str = r#"{:atoms ["C" "C#i13" "O"] :bonds [[1 2 "1"] [2 0 "1"]]}"#;
const RDKIT_ISOTOPE_MAPPING: &[u32] = &[1, 2, 0];
const RDKIT_ISOTOPE_ORBITS: &[&[u32]] = &[&[0], &[1], &[2]];
const RDKIT_NO_ISOTOPE_ORBITS: &[&[u32]] = &[&[0, 2], &[1]];

// InChI's test_permutation_util.cpp uses this five-atom graph to verify that
// permuting atom indices leaves the generated identifier invariant. Here it
// supplies the corresponding exact orbit and renumbering fixture without
// invoking InChI during the test.
const INCHI_PERMUTATION: &str = r#"{
    :atoms ["O" "C" "C" "O" "Y"]
    :bonds [[0 1 "1"] [0 4 "1"] [1 2 "1"] [1 3 "2"]]
}"#;
const INCHI_PERMUTATION_RENUMBERED: &str = r#"{
    :atoms ["C" "O" "Y" "O" "C"]
    :bonds [[3 0 "1"] [3 2 "1"] [0 4 "1"] [0 1 "2"]]
}"#;
const INCHI_PERMUTATION_MAPPING: &[u32] = &[3, 0, 4, 1, 2];
const INCHI_PERMUTATION_ORBITS: &[&[u32]] = &[&[0], &[1], &[2], &[3], &[4]];

const DATIVE: &str = r#"{
    :atoms ["N" "N" "B"] :bonds []
    :dative-bonds [{:donors [0 1] :acceptor 2 :attrs "1#R"}]
}"#;
const DATIVE_RENUMBERED: &str = r#"{
    :atoms ["N" "B" "N"] :bonds []
    :dative-bonds [{:donors [2 0] :acceptor 1 :attrs "1#R"}]
}"#;
const DATIVE_MAPPING: &[u32] = &[2, 0, 1];
const DATIVE_ORBITS: &[&[u32]] = &[&[0, 1], &[2]];

const AROMATIC: &str = r#"{
    :atoms ["C" "C" "C" "C" "C" "C"]
    :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#e6"}]
}"#;
const AROMATIC_RENUMBERED: &str = r#"{
    :atoms ["C" "C" "C" "C" "C" "C"]
    :bonds [[3 0 "1"] [0 5 "1"] [5 2 "1"] [2 4 "1"] [4 1 "1"] [1 3 "1"]]
    :aromatic-systems [{:atoms [3 0 5 2 4 1] :attrs "*#e6"}]
}"#;
const AROMATIC_MAPPING: &[u32] = &[3, 0, 5, 2, 4, 1];
const AROMATIC_ORBITS: &[&[u32]] = &[&[0, 1, 2, 3, 4, 5]];

const MULTICENTER: &str = r#"{
    :atoms ["B" "H" "B"] :bonds []
    :multicenter-bonds [{:atoms [0 1 2] :attrs "*#e2"}]
}"#;
const MULTICENTER_RENUMBERED: &str = r#"{
    :atoms ["H" "B" "B"] :bonds []
    :multicenter-bonds [{:atoms [2 0 1] :attrs "*#e2"}]
}"#;
const MULTICENTER_MAPPING: &[u32] = &[2, 0, 1];
const MULTICENTER_ORBITS: &[&[u32]] = &[&[0, 2], &[1]];

const NONCOVALENT: &str = r#"{
    :atoms ["O" "O" "N"] :bonds []
    :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]
}"#;
const NONCOVALENT_RENUMBERED: &str = r#"{
    :atoms ["O" "N" "O"] :bonds []
    :noncovalent-bonds [{:atoms [2 0] :attrs "Hbd"}]
}"#;
const NONCOVALENT_MAPPING: &[u32] = &[2, 0, 1];
const NONCOVALENT_ORBITS: &[&[u32]] = &[&[0, 1], &[2]];

const STEREO_ATOM: &str = r#"{
    :atoms ["C" "Cl" "Cl" "F" "Br"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
}"#;
const STEREO_ATOM_RENUMBERED: &str = r#"{
    :atoms ["F" "Cl" "C" "Br" "Cl"]
    :bonds [[2 4 "1"] [2 1 "1"] [2 0 "1"] [2 3 "1"]]
    :stereo-atoms [{:site 2 :ligands [4 1 0 3] :attrs "Th1"}]
}"#;
const STEREO_ATOM_MAPPING: &[u32] = &[2, 4, 1, 0, 3];
const STEREO_ATOM_PROPER_ORBITS: &[&[u32]] = &[&[0], &[1], &[2], &[3], &[4]];
const STEREO_ATOM_STAR_ORBITS: &[&[u32]] = &[&[0], &[1, 2], &[3], &[4]];

const STEREO_BOND: &str = r#"{
    :atoms ["C" "C" "F" "Cl" "F" "Cl"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]
}"#;
const STEREO_BOND_RENUMBERED: &str = r#"{
    :atoms ["Cl" "C" "F" "Cl" "C" "F"]
    :bonds [[4 1 "2"] [4 5 "1"] [4 0 "1"] [1 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [5 0 2 3] :attrs "Ct1"}]
}"#;
const STEREO_BOND_MAPPING: &[u32] = &[4, 1, 5, 0, 2, 3];
const STEREO_BOND_PROPER_ORBITS: &[&[u32]] = &[&[0, 1], &[2, 4], &[3, 5]];
const STEREO_BOND_STAR_ORBITS: &[&[u32]] = &[&[0, 1], &[2, 4], &[3, 5]];

fn graph_symmetry(molecule: &Molecule, include_isotopes: bool) -> GraphSymmetry {
    let features = if include_isotopes {
        MoleculeColoringFeatures::all()
    } else {
        MoleculeColoringFeatures::all().difference(MoleculeColoringFeatures::ISOTOPE)
    };
    molecule.graph_symmetry(&GraphSymmetryConfig {
        coloring: ConstitutionColoring::new(features),
        iterate_to_fixpoint: true,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    })
}

fn proper_orbits(molecule: &Molecule, include_isotopes: bool) -> Vec<Vec<AtomId>> {
    let symmetry = graph_symmetry(molecule, include_isotopes);
    orbit_partition(molecule, |id| symmetry.proper_orbit_of(id))
}

fn star_orbits(molecule: &Molecule) -> Vec<Vec<AtomId>> {
    let symmetry = graph_symmetry(molecule, true);
    orbit_partition(molecule, |id| symmetry.star_orbit_of(id))
}

fn orbit_partition(
    molecule: &Molecule,
    orbit_of: impl Fn(AtomId) -> Vec<AtomId>,
) -> Vec<Vec<AtomId>> {
    let mut visited = vec![false; molecule.atoms().count()];
    molecule
        .atoms()
        .ids()
        .filter_map(|id| {
            if visited[id.index()] {
                return None;
            }
            let orbit = orbit_of(id);
            for member in &orbit {
                visited[member.index()] = true;
            }
            Some(orbit)
        })
        .collect()
}

fn atom_partition(partition: &[&[u32]]) -> Vec<Vec<AtomId>> {
    partition
        .iter()
        .map(|orbit| orbit.iter().copied().map(AtomId).collect())
        .collect()
}

fn remap_partition(partition: &[&[u32]], mapping: &[u32]) -> Vec<Vec<AtomId>> {
    let mut remapped: Vec<Vec<AtomId>> = partition
        .iter()
        .map(|orbit| {
            let mut remapped_orbit: Vec<AtomId> = orbit
                .iter()
                .map(|&id| AtomId(mapping[id as usize]))
                .collect();
            remapped_orbit.sort();
            remapped_orbit
        })
        .collect();
    remapped.sort();
    remapped
}

#[rstest]
#[case::cdk_biphenyl(
    CDK_BIPHENYL,
    CDK_BIPHENYL_RENUMBERED,
    CDK_BIPHENYL_MAPPING,
    CDK_BIPHENYL_ORBITS,
    true
)]
#[case::rdkit_isotope_included(
    RDKIT_ISOTOPE,
    RDKIT_ISOTOPE_RENUMBERED,
    RDKIT_ISOTOPE_MAPPING,
    RDKIT_ISOTOPE_ORBITS,
    true
)]
#[case::rdkit_isotope_excluded(
    RDKIT_ISOTOPE,
    RDKIT_ISOTOPE_RENUMBERED,
    RDKIT_ISOTOPE_MAPPING,
    RDKIT_NO_ISOTOPE_ORBITS,
    false
)]
#[case::inchi_permutation(
    INCHI_PERMUTATION,
    INCHI_PERMUTATION_RENUMBERED,
    INCHI_PERMUTATION_MAPPING,
    INCHI_PERMUTATION_ORBITS,
    true
)]
#[case::dative(DATIVE, DATIVE_RENUMBERED, DATIVE_MAPPING, DATIVE_ORBITS, true)]
#[case::aromatic(AROMATIC, AROMATIC_RENUMBERED, AROMATIC_MAPPING, AROMATIC_ORBITS, true)]
#[case::multicenter(
    MULTICENTER,
    MULTICENTER_RENUMBERED,
    MULTICENTER_MAPPING,
    MULTICENTER_ORBITS,
    true
)]
#[case::noncovalent(
    NONCOVALENT,
    NONCOVALENT_RENUMBERED,
    NONCOVALENT_MAPPING,
    NONCOVALENT_ORBITS,
    true
)]
fn test_molecule_graph_symmetry_corpus(
    #[case] input: &str,
    #[case] renumbered_input: &str,
    #[case] mapping: &[u32],
    #[case] expected_partition: &[&[u32]],
    #[case] include_isotopes: bool,
) {
    let molecule = mol_dsl!(input);
    let renumbered = mol_dsl!(renumbered_input);

    assert_eq!(
        proper_orbits(&molecule, include_isotopes),
        atom_partition(expected_partition),
    );
    assert_eq!(
        proper_orbits(&renumbered, include_isotopes),
        remap_partition(expected_partition, mapping),
    );
}

#[rstest]
#[case::stereo_atom(
    STEREO_ATOM,
    STEREO_ATOM_RENUMBERED,
    STEREO_ATOM_MAPPING,
    STEREO_ATOM_PROPER_ORBITS,
    STEREO_ATOM_STAR_ORBITS
)]
#[case::stereo_bond(
    STEREO_BOND,
    STEREO_BOND_RENUMBERED,
    STEREO_BOND_MAPPING,
    STEREO_BOND_PROPER_ORBITS,
    STEREO_BOND_STAR_ORBITS
)]
fn test_molecule_graph_symmetry_stereo_corpus(
    #[case] input: &str,
    #[case] renumbered_input: &str,
    #[case] mapping: &[u32],
    #[case] expected_proper_partition: &[&[u32]],
    #[case] expected_star_partition: &[&[u32]],
) {
    let molecule = mol_dsl!(input);
    let renumbered = mol_dsl!(renumbered_input);

    assert_eq!(
        proper_orbits(&molecule, true),
        atom_partition(expected_proper_partition),
    );
    assert_eq!(
        proper_orbits(&renumbered, true),
        remap_partition(expected_proper_partition, mapping),
    );
    assert_eq!(
        star_orbits(&molecule),
        atom_partition(expected_star_partition),
    );
    assert_eq!(
        star_orbits(&renumbered),
        remap_partition(expected_star_partition, mapping),
    );
}
