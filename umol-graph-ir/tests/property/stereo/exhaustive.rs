//! Binary coset actions and framed molecule equality from parity/endpoint blocks,
//! and one-carrier classifications from ligand equality. Molecules are
//! constructed directly, without production frame transport or symmetry calculation.
//! Literature cases exhaust binary assignments and group complete canonical forms; their counts
//! and unit sets are from Razinger et al. (1993), Table I and Figure 1.

use std::array;
use std::collections::BTreeSet;

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, Canonicalize, CanonicalizeContext, GraphSymmetryConfig,
    Molecule, MoleculeEntries, Reframe, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    Stereogenicity,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_perm::PermutationGroup;
use StereoSite::{Atom as A, Bond as B};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum StereoSite {
    Atom(u32),
    Bond(u32, u32),
}

fn binary_assignment(molecule: &Molecule, assignment: usize) -> Molecule {
    let mut result = molecule.clone();
    let atom_ids = molecule.stereo_atoms().ids().collect::<Vec<_>>();
    let bond_ids = molecule.stereo_bonds().ids().collect::<Vec<_>>();
    for (position, id) in atom_ids.iter().enumerate() {
        assert_eq!(
            molecule.stereo_atom(*id).attributes().configuration.kind(),
            Some(StereoKind::Tetrahedral)
        );
        result.stereo_atom_mut(*id).attributes_mut().configuration =
            StereoConfigurationForm::kinded(
                StereoKind::Tetrahedral,
                ((assignment >> position) & 1) as u32,
            );
    }
    for (position, id) in bond_ids.iter().enumerate() {
        assert_eq!(
            molecule.stereo_bond(*id).attributes().configuration.kind(),
            Some(StereoKind::CisTrans)
        );
        result.stereo_bond_mut(*id).attributes_mut().configuration =
            StereoConfigurationForm::kinded(
                StereoKind::CisTrans,
                ((assignment >> (atom_ids.len() + position)) & 1) as u32,
            );
    }
    result
}

fn exhaustive_stereoisomers(molecule: &Molecule) -> ([usize; 3], BTreeSet<StereoSite>) {
    let context = CanonicalizeContext {
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let sites = molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| A(stereo.site_id().0))
        .chain(molecule.stereo_bonds().iter().map(|stereo| {
            let mut atoms = molecule.bond(stereo.site_id()).atom_ids();
            atoms.sort();
            B(atoms[0].0, atoms[1].0)
        }))
        .collect::<Vec<_>>();
    let mut representatives = Vec::new();
    let mut classes = Vec::new();
    for assignment in 0..1usize << sites.len() {
        let canonical = binary_assignment(molecule, assignment)
            .canonicalize(&context)
            .unwrap();
        let class = match representatives.iter().position(|other| *other == canonical) {
            Some(class) => class,
            None => {
                representatives.push(canonical);
                representatives.len() - 1
            }
        };
        classes.push(class);
    }
    let mirror_mask = (1usize << molecule.stereo_atoms().count()) - 1;
    let chiral = classes
        .iter()
        .enumerate()
        .filter(|&(assignment, class)| *class != classes[assignment ^ mirror_mask])
        .map(|(_, &class)| class)
        .collect::<BTreeSet<_>>();
    let units = sites
        .into_iter()
        .enumerate()
        .filter(|(position, _)| {
            (0..classes.len())
                .any(|assignment| classes[assignment] != classes[assignment ^ (1 << position)])
        })
        .map(|(_, site)| site)
        .collect();
    (
        [
            representatives.len(),
            chiral.len(),
            representatives.len() - chiral.len(),
        ],
        units,
    )
}

fn exhaustive_coset_classes(molecule: &Molecule, kind: StereoKind) -> Vec<u32> {
    let context = CanonicalizeContext {
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    assert_eq!(molecule.stereo_atoms().count(), 1);
    let mut representatives = Vec::new();
    (0..kind.count() as u32)
        .map(|coset| {
            let mut assigned = molecule.clone();
            assigned
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::kinded(kind, coset);
            let canonical = assigned.canonicalize(&context).unwrap();
            match representatives.iter().position(|other| *other == canonical) {
                Some(class) => class as u32,
                None => {
                    representatives.push(canonical);
                    representatives.len() as u32 - 1
                }
            }
        })
        .collect()
}

pub(super) fn binary_coset_action(kind: StereoKind, coset: u32, image: &[usize; 4]) -> Option<u32> {
    let flip = match kind {
        StereoKind::Tetrahedral => {
            let inversions = (0..4)
                .flat_map(|a| (a + 1..4).map(move |b| (a, b)))
                .filter(|&(a, b)| image[a] > image[b])
                .count();
            (inversions % 2) as u32
        }
        StereoKind::CisTrans => {
            if image[0] / 2 != image[1] / 2 || image[2] / 2 != image[3] / 2 {
                return None;
            }
            ((image[0] % 2) ^ (image[2] % 2)) as u32
        }
        _ => return None,
    };
    Some(coset ^ flip)
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
fn test_stereo_kind_act(#[case] kind: StereoKind) {
    for &permutation in PermutationGroup::symmetric(4).elements() {
        let image = array::from_fn(|position| permutation.apply(position));
        for coset in 0..2 {
            assert_eq!(
                kind.act(coset, permutation),
                binary_coset_action(kind, coset, &image),
                "kind {kind:?}, coset {coset}, image {image:?}"
            );
        }
    }
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
fn test_molecule_framed_eq(#[case] kind: StereoKind, #[values(0, 1)] coset: u32) {
    let offset = if kind == StereoKind::Tetrahedral {
        1
    } else {
        2
    };
    let mut entries = MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); offset],
        ..Default::default()
    };
    entries.atoms.extend(
        [Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element),
    );
    if kind == StereoKind::CisTrans {
        entries
            .bonds
            .push((AtomId(0), AtomId(1), BondForm::from_order(2)));
    }
    for position in 0..4 {
        entries.bonds.push((
            AtomId::from(if kind == StereoKind::Tetrahedral {
                0
            } else {
                position / 2
            }),
            AtomId::from(offset + position),
            BondForm::from_order(1),
        ));
    }
    let ligands = (offset..offset + 4)
        .map(|atom| StereoLigand::new(AtomId::from(atom), StereoLigandKind::Atom))
        .collect();
    if kind == StereoKind::Tetrahedral {
        entries
            .stereo_atoms
            .push((AtomId(0), ligands, StereoAtomForm::new(kind, coset)));
    } else {
        entries
            .stereo_bonds
            .push((BondId(0), ligands, StereoBondForm::new(kind, coset)));
    }
    let source = Molecule::from_entries(entries.clone());

    for &permutation in PermutationGroup::symmetric(4).elements() {
        let image = array::from_fn(|position| permutation.apply(position));
        let Some(transported_coset) = binary_coset_action(kind, coset, &image) else {
            continue;
        };
        for other_coset in 0..2 {
            let ligands = image
                .into_iter()
                .map(|position| {
                    StereoLigand::new(AtomId::from(offset + position), StereoLigandKind::Atom)
                })
                .collect();
            let mut other_entries = entries.clone();
            if kind == StereoKind::Tetrahedral {
                other_entries.stereo_atoms[0] =
                    (AtomId(0), ligands, StereoAtomForm::new(kind, other_coset));
            } else {
                other_entries.stereo_bonds[0] =
                    (BondId(0), ligands, StereoBondForm::new(kind, other_coset));
            }
            let other = Molecule::from_entries(other_entries);
            let expected = other_coset == transported_coset;

            assert_eq!(
                source.framed_eq(&other),
                expected,
                "kind {kind:?}, coset {coset}, image {image:?}, other coset {other_coset}"
            );
            assert_eq!(
                other.framed_eq(&source),
                expected,
                "kind {kind:?}, coset {coset}, image {image:?}, other coset {other_coset}"
            );
        }
    }
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
fn test_stereo_view_stereogenicity(
    #[case] kind: StereoKind,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
) {
    let labels = [Element::H, Element::F, Element::Cl, Element::Br];
    for pattern in 0..256 {
        let elements: [Element; 4] =
            array::from_fn(|position| labels[(pattern >> (2 * position)) & 3]);
        let offset = if kind == StereoKind::Tetrahedral {
            1
        } else {
            2
        };
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); offset],
            ..Default::default()
        };
        entries
            .atoms
            .extend(elements.into_iter().map(AtomForm::from_element));
        if kind == StereoKind::CisTrans {
            entries
                .bonds
                .push((AtomId(0), AtomId(1), BondForm::from_order(2)));
        }
        for position in 0..4 {
            entries.bonds.push((
                AtomId::from(if kind == StereoKind::Tetrahedral {
                    0
                } else {
                    position / 2
                }),
                AtomId::from(offset + position),
                BondForm::from_order(1),
            ));
        }
        let ligands = (offset..offset + 4)
            .map(|atom| StereoLigand::new(AtomId::from(atom), StereoLigandKind::Atom))
            .collect();
        let expected = if kind == StereoKind::Tetrahedral {
            entries.stereo_atoms.push((
                AtomId(0),
                ligands,
                StereoAtomForm::new(kind, coset.clone()),
            ));
            let equal_pairs = (0..4)
                .flat_map(|first| (first + 1..4).map(move |second| (first, second)))
                .filter(|&(first, second)| elements[first] == elements[second])
                .count();
            match equal_pairs {
                0 => Stereogenicity::Stereogenic,
                1 => Stereogenicity::Prostereogenic,
                _ => Stereogenicity::Nonstereogenic,
            }
        } else {
            entries.stereo_bonds.push((
                BondId(0),
                ligands,
                StereoBondForm::new(kind, coset.clone()),
            ));
            match (elements[0] == elements[1], elements[2] == elements[3]) {
                (false, false) => Stereogenicity::Stereogenic,
                (true, true) => Stereogenicity::Nonstereogenic,
                _ => Stereogenicity::Prostereogenic,
            }
        };
        let molecule = Molecule::from_entries(entries);
        let symmetry = molecule
            .graph_symmetry(&GraphSymmetryConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                ..Default::default()
            })
            .into_graph_symmetry();
        let (actual, stored_coset) = if kind == StereoKind::Tetrahedral {
            let stereo = molecule
                .stereo_atom(StereoAtomId(0))
                .with_symmetry(&symmetry);
            (stereo.stereogenicity(), stereo.coset())
        } else {
            let stereo = molecule
                .stereo_bond(StereoBondId(0))
                .with_symmetry(&symmetry);
            (stereo.stereogenicity(), stereo.coset())
        };
        let action_is_stereogenic = PermutationGroup::symmetric(4)
            .elements()
            .iter()
            .filter(|permutation| {
                let image = array::from_fn(|position| permutation.apply(position));
                (0..4).all(|position| elements[position] == elements[image[position]])
                    && binary_coset_action(kind, 0, &image).is_some()
            })
            .all(|permutation| {
                let image = array::from_fn(|position| permutation.apply(position));
                (0..2).all(|index| binary_coset_action(kind, index, &image) == Some(index))
            });
        assert_eq!(
            action_is_stereogenic,
            expected == Stereogenicity::Stereogenic,
            "kind {kind:?}, coset {coset:?}, elements {elements:?}"
        );
        assert_eq!(actual == Stereogenicity::Stereogenic, action_is_stereogenic);
        assert_eq!(stored_coset, &coset);
        assert_eq!(
            actual, expected,
            "kind {kind:?}, coset {coset:?}, elements {elements:?}"
        );
    }
}

// Razinger et al., J. Chem. Inf. Comput. Sci. 1993, 33, 812–825, Table I and Figure 1.

/// Razinger et al. (1993), Figure 7 no. 1.
#[fixture]
fn bromochloroethanol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h3" "Br" "Cl" "O#h"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 13.
#[fixture]
fn bis_fluoroethyl_difluorohexane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C#h3" "C#h3" "C#h3" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                "F" "F" "F" "F"]
        :bonds [[8 9 "1"] [8 4 "1"] [8 5 "1"] [9 6 "1"] [9 7 "1"] [4 0 "1"] [5 1 "1"]
                [6 2 "1"] [7 3 "1"] [4 10 "1"] [5 11 "1"] [6 12 "1"] [7 13 "1"]]
        :stereo-atoms [{:site 4 :ligands [0 8 10 [:h 4]] :attrs "Th*"}
                       {:site 5 :ligands [1 8 11 [:h 5]] :attrs "Th*"}
                       {:site 6 :ligands [2 9 12 [:h 6]] :attrs "Th*"}
                       {:site 7 :ligands [3 9 13 [:h 7]] :attrs "Th*"}
                       {:site 8 :ligands [4 5 9 [:h 8]] :attrs "Th*"}
                       {:site 9 :ligands [6 7 8 [:h 9]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 14.
#[fixture]
fn dimethylcyclopentane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"] [0 5 "1"] [1 6 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 4 5 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 6 [:h 1]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 16.
#[fixture]
fn dimethylcyclohexane_1_2() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 17.
#[fixture]
fn dimethylcyclohexane_1_3() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 7 [:h 2]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 18.
#[fixture]
fn dimethylcyclohexane_1_4() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [3 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 7 [:h 3]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 19.
#[fixture]
fn trimethylcyclohexane_1_2_3() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h3" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"] [2 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 8 [:h 2]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 22.
#[fixture]
fn tetramethylcyclobutane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]
                [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 5 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 6 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [0 2 7 [:h 3]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 28, Figure 1 d.
#[fixture]
fn cyclohexanetriol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h" "C#h2" "C#h" "C#h2" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 7 [:h 2]] :attrs "Th*"}
                       {:site 4 :ligands [3 5 8 [:h 4]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 30.
#[fixture]
fn butene() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C#h" "C#h" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 38.
#[fixture]
fn bis_fluoroethyl_difluorohexene() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h3" "C#h3" "C#h3" "C#h3" "C#h" "C#h" "C#h" "C#h" "C" "C"
                "F" "F" "F" "F"]
        :bonds [[8 9 "2"] [8 4 "1"] [8 5 "1"] [9 6 "1"] [9 7 "1"] [4 0 "1"] [5 1 "1"]
                [6 2 "1"] [7 3 "1"] [4 10 "1"] [5 11 "1"] [6 12 "1"] [7 13 "1"]]
        :stereo-atoms [{:site 4 :ligands [0 8 10 [:h 4]] :attrs "Th*"}
                       {:site 5 :ligands [1 8 11 [:h 5]] :attrs "Th*"}
                       {:site 6 :ligands [2 9 12 [:h 6]] :attrs "Th*"}
                       {:site 7 :ligands [3 9 13 [:h 7]] :attrs "Th*"}]
        :stereo-bonds [{:site 0 :ligands [4 5 6 7] :attrs "Ct*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 48, Figure 1 h.
#[fixture]
fn decalin() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 0 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 9 [:h 0]] :attrs "Th*"}
                       {:site 5 :ligands [0 4 6 [:h 5]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 7 no. 69, Figure 1 e.
#[fixture]
fn inositol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#h" "O#h" "O#h" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"] [2 8 "1"] [3 9 "1"] [4 10 "1"] [5 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 8 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 9 [:h 3]] :attrs "Th*"}
                       {:site 4 :ligands [3 5 10 [:h 4]] :attrs "Th*"}
                       {:site 5 :ligands [0 4 11 [:h 5]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 a.
#[fixture]
fn diethylidenecyclohexane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C#h2" "C" "C#h2" "C#h2" "C#h" "C#h3" "C#h" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "2"] [6 7 "1"] [3 8 "2"] [8 9 "1"]]
        :stereo-bonds [{:site 6 :ligands [1 5 7 [:h 6]] :attrs "Ct*"}
                       {:site 8 :ligands [2 4 9 [:h 8]] :attrs "Ct*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 c.
#[fixture]
fn cyclohexanediol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "O#h" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [3 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 7 [:h 3]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 f.
#[fixture]
fn ethylidenemethylspiroundecane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C#h2" "C" "C#h2" "C#h2" "C#h2" "C#h2" "C#h" "C#h2" "C#h2"
                "C#h" "C#h3" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [3 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 3 "1"]
                [0 11 "2"] [11 12 "1"] [8 13 "1"]]
        :stereo-atoms [{:site 3 :ligands [2 4 6 10] :attrs "Th*"}
                       {:site 8 :ligands [7 9 13 [:h 8]] :attrs "Th*"}]
        :stereo-bonds [{:site 12 :ligands [1 5 12 [:h 11]] :attrs "Ct*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 g.
#[fixture]
fn ethylidenespiroundecane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h2" "C#h2" "C" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"
                "C#h" "C#h3"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [3 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 3 "1"]
                [0 11 "2"] [11 12 "1"]]
        :stereo-atoms [{:site 3 :ligands [2 4 6 10] :attrs "Th*"}]
        :stereo-bonds [{:site 12 :ligands [1 5 12 [:h 11]] :attrs "Ct*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 i.
#[fixture]
fn methylcyclohexanediol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h" "C#h2" "C#h" "C#h2" "C#h3" "O#h" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 7 [:h 2]] :attrs "Th*"}
                       {:site 4 :ligands [3 5 8 [:h 4]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 j.
#[fixture]
fn tris_hydroxybutyl_methane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "C#h" "C#h" "C#h"
                "C#h3" "C#h3" "C#h3" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [2 5 "1"] [3 6 "1"]
                [1 7 "1"] [2 8 "1"] [3 9 "1"] [7 10 "1"] [8 11 "1"] [9 12 "1"]
                [7 13 "1"] [8 14 "1"] [9 15 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 4 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [0 5 8 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [0 6 9 [:h 3]] :attrs "Th*"}
                       {:site 7 :ligands [1 10 13 [:h 7]] :attrs "Th*"}
                       {:site 8 :ligands [2 11 14 [:h 8]] :attrs "Th*"}
                       {:site 9 :ligands [3 12 15 [:h 9]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 k.
#[fixture]
fn tris_hydroxyethyl_methane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [2 5 "1"] [3 6 "1"]
                [1 7 "1"] [2 8 "1"] [3 9 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 4 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [0 5 8 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [0 6 9 [:h 3]] :attrs "Th*"}]
    }"#
    )
}

/// Razinger et al. (1993), Figure 1 n.
#[fixture]
fn bis_methylcyclohexyl_methanol() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2"
                "C#h" "C#h3" "C#h3" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]
                [0 12 "1"] [6 12 "1"] [3 13 "1"] [9 14 "1"] [12 15 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 12 [:h 0]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 13 [:h 3]] :attrs "Th*"}
                       {:site 6 :ligands [7 11 12 [:h 6]] :attrs "Th*"}
                       {:site 9 :ligands [8 10 14 [:h 9]] :attrs "Th*"}
                       {:site 12 :ligands [0 6 15 [:h 12]] :attrs "Th*"}]
    }"#
    )
}

#[fixture]
fn diamminedichloroplatinum_cis() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
    }"#
    )
}

#[fixture]
fn diamminedichloroplatinum_trans() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp1"}]
    }"#
    )
}

#[fixture]
fn chlorotetrafluorophosphorane_axial() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["P" "Cl" "F" "F" "F" "F"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb0"}]
    }"#
    )
}

#[fixture]
fn chlorotetrafluorophosphorane_equatorial() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["P" "Cl" "F" "F" "F" "F"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb8"}]
    }"#
    )
}

#[rstest]
#[case::bromochloroethanol(bromochloroethanol(), Some([2, 2, 0]), vec![A(0)], vec![])]
#[case::bis_fluoroethyl_difluorohexane(bis_fluoroethyl_difluorohexane(), Some([10, 6, 4]), vec![A(4), A(5), A(6), A(7)], vec![A(8), A(9)])]
#[case::dimethylcyclopentane(dimethylcyclopentane(), Some([3, 2, 1]), vec![A(0), A(1)], vec![])]
#[case::dimethylcyclohexane_1_2(dimethylcyclohexane_1_2(), Some([3, 2, 1]), vec![A(0), A(1)], vec![])]
#[case::dimethylcyclohexane_1_3(dimethylcyclohexane_1_3(), Some([3, 2, 1]), vec![A(0), A(2)], vec![])]
#[case::dimethylcyclohexane_1_4(dimethylcyclohexane_1_4(), Some([2, 0, 2]), vec![], vec![A(0), A(3)])]
#[case::trimethylcyclohexane_1_2_3(trimethylcyclohexane_1_2_3(), Some([4, 2, 2]), vec![A(0), A(2)], vec![A(1)])]
#[case::tetramethylcyclobutane(tetramethylcyclobutane(), Some([4, 0, 4]), vec![], vec![A(0), A(1), A(2), A(3)])]
#[case::cyclohexanetriol(cyclohexanetriol(), Some([2, 0, 2]), vec![], vec![A(0), A(2), A(4)])]
#[case::butene(butene(), Some([2, 0, 2]), vec![B(1, 2)], vec![])]
#[case::bis_fluoroethyl_difluorohexene(bis_fluoroethyl_difluorohexene(), Some([7, 4, 3]), vec![A(4), A(5), A(6), A(7)], vec![B(8, 9)])]
#[case::decalin(decalin(), Some([2, 0, 2]), vec![], vec![A(0), A(5)])]
#[case::inositol(inositol(), Some([9, 2, 7]), vec![], vec![A(0), A(1), A(2), A(3), A(4), A(5)])]
#[case::diethylidenecyclohexane(diethylidenecyclohexane(), None, vec![], vec![B(0, 6), B(3, 8)])]
#[case::cyclohexanediol(cyclohexanediol(), None, vec![], vec![A(0), A(3)])]
#[case::ethylidenemethylspiroundecane(ethylidenemethylspiroundecane(), None, vec![], vec![A(3), A(8), B(0, 11)])]
#[case::ethylidenespiroundecane(ethylidenespiroundecane(), None, vec![], vec![])]
#[case::methylcyclohexanediol(methylcyclohexanediol(), None, vec![A(2), A(4)], vec![A(0)])]
#[case::tris_hydroxybutyl_methane(tris_hydroxybutyl_methane(), None, vec![A(1), A(2), A(3), A(7), A(8), A(9)], vec![A(0)])]
#[case::tris_hydroxyethyl_methane(tris_hydroxyethyl_methane(), None, vec![A(1), A(2), A(3)], vec![])]
#[case::bis_methylcyclohexyl_methanol(bis_methylcyclohexyl_methanol(), None, vec![], vec![A(0), A(3), A(6), A(9), A(12)])]
fn test_exhaustive_stereoisomers(
    #[case] molecule: Molecule,
    #[case] expected_counts: Option<[usize; 3]>,
    #[case] constitutional: Vec<StereoSite>,
    #[case] configurational: Vec<StereoSite>,
) {
    let original = molecule.clone();
    assert_eq!(molecule.clone().edit().finish(), Ok(molecule.clone()));
    let (counts, units) = exhaustive_stereoisomers(&molecule);
    if let Some(expected) = expected_counts {
        assert_eq!(counts, expected);
    }
    assert_eq!(
        units,
        constitutional.into_iter().chain(configurational).collect()
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::square_planar_cis(diamminedichloroplatinum_cis(), StereoKind::SquarePlanar, vec![0, 1, 0], 0)]
#[case::square_planar_trans(diamminedichloroplatinum_trans(), StereoKind::SquarePlanar, vec![0, 1, 0], 1)]
#[case::trigonal_bipyramidal_axial(chlorotetrafluorophosphorane_axial(), StereoKind::TrigonalBipyramidal, (0..20).map(|coset| u32::from(coset >= 8)).collect(), 0)]
#[case::trigonal_bipyramidal_equatorial(chlorotetrafluorophosphorane_equatorial(), StereoKind::TrigonalBipyramidal, (0..20).map(|coset| u32::from(coset >= 8)).collect(), 1)]
fn test_exhaustive_coset_classes(
    #[case] molecule: Molecule,
    #[case] kind: StereoKind,
    #[case] expected: Vec<u32>,
    #[case] expected_stored: u32,
) {
    assert_eq!(molecule.clone().edit().finish(), Ok(molecule.clone()));
    let actual = exhaustive_coset_classes(&molecule, kind);
    let stored = molecule.stereo_atom(StereoAtomId(0));
    let StereoConfigurationForm::Kinded(stored_kind, StereoCoset::Lit(coset)) =
        &stored.attributes().configuration
    else {
        panic!("the named nonbinary fixtures require literal configurations");
    };
    assert_eq!(*stored_kind, kind);
    assert_eq!(actual[*coset as usize], expected_stored);
    assert_eq!(actual, expected);
}
