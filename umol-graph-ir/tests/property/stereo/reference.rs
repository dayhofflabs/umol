//! Binary coset actions from parity/endpoint blocks and one-carrier classifications
//! from ligand equality. These references use no production symmetry calculation.

use std::array;

use rstest::rstest;
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, ConstitutionColoring, GraphSymmetryConfig, Molecule,
    MoleculeEntries, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, Stereogenicity,
};
use umol_perm::PermutationGroup;

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
fn test_stereo_symmetry_stereogenicity(
    #[case] kind: StereoKind,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
    #[values(false, true)] iterate_to_fixpoint: bool,
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
                1 => Stereogenicity::Prochiral,
                _ => Stereogenicity::Symmetric,
            }
        } else {
            entries.stereo_bonds.push((
                BondId(0),
                ligands,
                StereoBondForm::new(kind, coset.clone()),
            ));
            if elements[0] != elements[1] && elements[2] != elements[3] {
                Stereogenicity::Stereogenic
            } else {
                Stereogenicity::Symmetric
            }
        };
        let molecule = Molecule::from_entries(entries);
        let graph = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint,
            max_iterations: 16,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = if kind == StereoKind::Tetrahedral {
            molecule.stereo_atom_symmetry(&graph, StereoAtomId(0))
        } else {
            molecule.stereo_bond_symmetry(&graph, StereoBondId(0))
        };
        let action_is_stereogenic = symmetry.group().elements().iter().all(|operation| {
            let image = array::from_fn(|position| operation.permutation().apply(position));
            (0..2).all(|index| binary_coset_action(kind, index, &image) == Some(index))
        });
        assert_eq!(
            action_is_stereogenic,
            expected == Stereogenicity::Stereogenic,
            "kind {kind:?}, coset {coset:?}, elements {elements:?}"
        );
        assert_eq!(symmetry.is_stereogenic(), action_is_stereogenic);
        assert_eq!(symmetry.coset(), &coset);
        assert_eq!(
            symmetry.stereogenicity(),
            expected,
            "kind {kind:?}, coset {coset:?}, elements {elements:?}"
        );
    }
}
