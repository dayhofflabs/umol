//! Substructure matching across algorithms on generated molecules, with exact
//! stereo matches specified independently by tetrahedral parity and endpoint swaps.

use proptest::prelude::*;
use rstest::rstest;
use umol_chem::element::Element;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};
use umol_graph_ir::ir::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, Molecule, MoleculeEntries, StereoAtomForm, StereoBondForm,
    StereoKind, StereoLigand, StereoLigandKind, SubstructureMatchAlgorithm,
    SubstructureMatchConfig,
};

use crate::strategies::molecule_strategy;

const SUBISO: [SubgraphIsomorphismAlgorithm; 6] = [
    Vf2,
    Ullmann,
    Ri,
    ArcMatch {
        path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
    },
    Vf2Rdkit,
    RayKirsch,
];
const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

fn sorted_matches(
    pattern: &Molecule,
    host: &Molecule,
    strategy: SubstructureMatchAlgorithm,
    subiso: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<AtomId>> {
    let mut occurrences: Vec<Vec<AtomId>> = pattern
        .substructure_matches(
            host,
            SubstructureMatchConfig {
                match_algorithm: strategy,
                subgraph_isomorphism_algorithm: subiso,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            },
        )
        .expect("generated molecules carry no molecule-scope constraints")
        .iter()
        .map(|c| {
            c.atoms()
                .matched_pairs()
                .iter()
                .map(|&(_, host)| host)
                .collect()
        })
        .collect();
    occurrences.sort();
    occurrences
}

#[rstest]
#[case::tetrahedral_same_frame(StereoKind::Tetrahedral, [0, 1, 2, 3], 0)]
#[case::tetrahedral_odd_frame(StereoKind::Tetrahedral, [1, 0, 2, 3], 1)]
#[case::tetrahedral_even_frame(StereoKind::Tetrahedral, [1, 2, 0, 3], 0)]
#[case::cis_trans_same_frame(StereoKind::CisTrans, [0, 1, 2, 3], 0)]
#[case::cis_trans_one_end(StereoKind::CisTrans, [1, 0, 2, 3], 1)]
#[case::cis_trans_both_ends(StereoKind::CisTrans, [1, 0, 3, 2], 0)]
#[case::cis_trans_exchanged_ends(StereoKind::CisTrans, [2, 3, 0, 1], 0)]
fn test_molecule_substructure_matches_stereo(
    #[case] kind: StereoKind,
    #[case] ligand_order: [usize; 4],
    #[case] transported_coset: u32,
    #[values(false, true)] opposite: bool,
    #[values(false, true)] relabel: bool,
) {
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
            .push((AtomId(0), ligands, StereoAtomForm::new(kind, 0u32)));
    } else {
        entries
            .stereo_bonds
            .push((BondId(0), ligands, StereoBondForm::new(kind, 0u32)));
    }
    let pattern = Molecule::from_entries(entries.clone());
    if relabel {
        entries.atoms.swap(offset, offset + 1);
    }
    let ligands = ligand_order
        .into_iter()
        .map(|position| {
            let position = match (relabel, position) {
                (true, 0) => 1,
                (true, 1) => 0,
                _ => position,
            };
            StereoLigand::new(AtomId::from(offset + position), StereoLigandKind::Atom)
        })
        .collect();
    let coset = transported_coset ^ u32::from(opposite);
    if kind == StereoKind::Tetrahedral {
        entries.stereo_atoms[0] = (AtomId(0), ligands, StereoAtomForm::new(kind, coset));
    } else {
        entries.stereo_bonds[0] = (BondId(0), ligands, StereoBondForm::new(kind, coset));
    }
    let host = Molecule::from_entries(entries);
    let mut mapping: Vec<_> = (0..offset + 4).map(AtomId::from).collect();
    if relabel {
        mapping.swap(offset, offset + 1);
    }
    let expected = if opposite { Vec::new() } else { vec![mapping] };

    for strategy in STRATEGIES {
        for subiso in SUBISO {
            assert_eq!(
                sorted_matches(&pattern, &host, strategy, subiso),
                expected,
                "{strategy:?}/{subiso:?}"
            );
            assert_eq!(
                sorted_matches(&host, &pattern, strategy, subiso),
                expected,
                "{strategy:?}/{subiso:?}"
            );
        }
    }
}

#[rstest]
#[case::atom(StereoLigand::new(AtomId(4), StereoLigandKind::Atom), false)]
#[case::hydrogen(
    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
    false
)]
#[case::lone_pair(StereoLigand::new(AtomId(0), StereoLigandKind::LonePair), false)]
#[case::virtual_ligand(StereoLigand::new(AtomId(0), StereoLigandKind::Virtual), true)]
fn test_molecule_substructure_matches_virtual_ligand(
    #[case] ligand: StereoLigand,
    #[case] matches: bool,
) {
    let mut entries = MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: (1..=4)
            .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    };
    let pattern = Molecule::from_entries(entries.clone());
    entries.stereo_atoms[0].1[3] = ligand;
    let host = Molecule::from_entries(entries);
    let expected = if matches {
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)]]
    } else {
        vec![]
    };
    for strategy in STRATEGIES {
        for subiso in SUBISO {
            assert_eq!(
                sorted_matches(&pattern, &host, strategy, subiso),
                expected,
                "{strategy:?}/{subiso:?}"
            );
            assert_eq!(
                sorted_matches(&host, &pattern, strategy, subiso),
                expected,
                "{strategy:?}/{subiso:?}"
            );
        }
    }
}

proptest! {
    #[test]
    fn test_substructure_cross_validation(
        host in molecule_strategy(),
        pattern in molecule_strategy(),
    ) {
        let reference = sorted_matches(&pattern, &host, GraphAndOverlays, Vf2);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                prop_assert_eq!(
                    sorted_matches(&pattern, &host, strategy, subiso),
                    reference.clone(),
                    "{:?}/{:?}", strategy, subiso
                );
            }
        }
    }

    #[test]
    fn test_substructure_cross_validation_planted(
        (host, subset) in molecule_strategy()
            .prop_filter("non-empty", |m| m.atoms().count() > 0)
            .prop_flat_map(|m| {
                let n = m.atoms().count();
                let ids: Vec<AtomId> = (0..n as u32).map(AtomId).collect();
                (Just(m), prop::sample::subsequence(ids, 1..=n))
            }),
    ) {
        // A pattern induced from `host` shares its structure, so it exercises the
        // match path far more than independent random pairs. (No `!is_empty` check:
        // the generator may emit constraints inconsistent with the topology — e.g. a
        // stored valence differing from the derived one — so a self-match is not
        // guaranteed; only cross-strategy / cross-algorithm agreement is.)
        let pattern = host.extract(&host.induced_subgraph(&subset));
        let reference = sorted_matches(&pattern, &host, GraphAndOverlays, Vf2);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                prop_assert_eq!(
                    sorted_matches(&pattern, &host, strategy, subiso),
                    reference.clone(),
                    "{:?}/{:?}", strategy, subiso
                );
            }
        }
    }
}
