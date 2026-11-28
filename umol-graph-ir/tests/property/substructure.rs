//! Substructure matching properties across strategies, subgraph-isomorphism algorithms,
//! automorphism algorithms, prepared carriers, and symmetry-reduction modes.
//!
//! # Validated properties and operational domains
//!
//! 1. **Exact stereo reframing (`test_molecule_substructure_matches_stereo`,
//!    `test_molecule_substructure_matches_virtual_ligand`):** tetrahedral parity, cis/trans
//!    endpoint swaps, and virtual-ligand kinds checked against independently computed coset
//!    actions across both strategies and all four `SubgraphIsomorphismAlgorithm` variants.
//! 2. **Stereo-preserving symmetry fixtures (`test_substructure_stereo_symmetry`):** named
//!    achiral, chiral, tetrahedral `A4`, `(R,R)` vs. `meso`, and cis/trans molecules across all
//!    four `(reduce_pattern_symmetry, reduce_host_symmetry)` modes, both strategies, all four
//!    `SubgraphIsomorphismAlgorithm` variants, and both `AutomorphismAlgorithm` variants.
//! 3. **Cross-validation on generated molecules (`test_substructure_cross_validation`,
//!    `test_substructure_cross_validation_planted`):** arbitrary and induced-subgraph molecules
//!    from `molecule_strategy()` validating `GraphAndOverlays` == `Incidence`, all four
//!    algorithms, mode splits (`visit_substructure_embeddings`, `has_substructure_match`,
//!    `first_substructure_match`, `count_substructure_matches`), and direct vs.
//!    `PreparedSubstructurePattern`/`PreparedSubstructureHost` parity.
//! 4. **Symmetry-reduced cross-validation (`test_substructure_symmetry_modes_planted`):**
//!    planted pairs from `molecule_strategy()` across all four symmetry modes verifying
//!    cross-algorithm and `TrailIr` == `Nauty` canonical-representative invariance within each
//!    strategy, `GraphAndOverlays` == `Incidence` count/existence parity, and direct vs. prepared
//!    parity.
//! 5. **Independent definition-level orbit validation (`test_substructure_symmetry_exact_orbits`):**
//!    small symmetric molecules (`n <= 5`) with elements, bond orders, aromatic overlays, and
//!    stereo atoms/bonds checked against an independent exhaustive atom-permutation enumeration
//!    of proper automorphisms and exact `(G_Q x G_H)`-orbit partitions.

use std::collections::BTreeSet;
use std::ops::ControlFlow;

use proptest::prelude::*;
use rstest::rstest;
use umol_chem::element::Element;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
use umol_graph_core::{
    AutomorphismAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomForm, AtomId, BondForm, BondId, ElectronCountsForm, FrameTransport,
    Molecule, MoleculeCorrespondence, MoleculeEntries, PreparedSubstructureHost,
    PreparedSubstructurePattern, StereoAtomForm, StereoBondForm, StereoKind, StereoLigand,
    StereoLigandKind, SubstructureMatchAlgorithm, SubstructureMatchConfig,
};
use umol_graph_ir::mol_dsl;
use umol_perm::Permutation;

use crate::strategies::molecule_strategy;

const SUBISO: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];
const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];
const AUTO_ALGS: [AutomorphismAlgorithm; 2] =
    [AutomorphismAlgorithm::TrailIr, AutomorphismAlgorithm::Nauty];
const SYMMETRY_MODES: [(bool, bool); 4] =
    [(false, false), (true, false), (false, true), (true, true)];

fn sorted_atom_embeddings(matches: &[MoleculeCorrespondence]) -> Vec<Vec<AtomId>> {
    let mut occurrences: Vec<Vec<AtomId>> = matches
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

fn sorted_matches(
    pattern: &Molecule,
    host: &Molecule,
    strategy: SubstructureMatchAlgorithm,
    subiso: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<AtomId>> {
    let matches = pattern.substructure_matches(
        host,
        SubstructureMatchConfig {
            match_algorithm: strategy,
            subgraph_isomorphism_algorithm: subiso,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        },
    );
    sorted_atom_embeddings(&matches)
}

/// Independent definition-level enumeration of every proper (`Orientation::Proper`)
/// atom automorphism of a small normalized molecule (`n <= 6`) by testing all `n!`
/// bijections directly against topology, atom/bond forms, aromatic overlays, and
/// stereo records.
fn exhaustive_proper_automorphisms(molecule: &Molecule) -> Vec<Vec<usize>> {
    let n = molecule.atoms().count();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut result = Vec::new();
    let mut stack = vec![0usize; n];
    if is_proper_atom_automorphism(molecule, &perm) {
        result.push(perm.clone());
    }
    let mut i = 0usize;
    while i < n {
        if stack[i] < i {
            if i.is_multiple_of(2) {
                perm.swap(0, i);
            } else {
                perm.swap(stack[i], i);
            }
            if is_proper_atom_automorphism(molecule, &perm) {
                result.push(perm.clone());
            }
            stack[i] += 1;
            i = 0;
        } else {
            stack[i] = 0;
            i += 1;
        }
    }
    result
}

fn is_proper_atom_automorphism(molecule: &Molecule, perm: &[usize]) -> bool {
    for atom in molecule.atoms().iter() {
        let image = AtomId(perm[atom.id().index()] as u32);
        if atom.attributes() != molecule.atom(image).attributes() {
            return false;
        }
    }

    for bond in molecule.bonds().iter() {
        let [a, b] = bond.atom_ids();
        let ma = AtomId(perm[a.index()] as u32);
        let mb = AtomId(perm[b.index()] as u32);
        let Some(mapped_bond) = molecule.bonds().iter().find(|other| {
            let [x, y] = other.atom_ids();
            (x == ma && y == mb) || (x == mb && y == ma)
        }) else {
            return false;
        };
        if bond.attributes() != mapped_bond.attributes() {
            return false;
        }
    }

    for system in molecule.aromatic_systems().iter() {
        let mut mapped_set: Vec<AtomId> = system
            .atom_ids()
            .map(|a| AtomId(perm[a.index()] as u32))
            .collect();
        mapped_set.sort();
        let Some(target) = molecule.aromatic_systems().iter().find(|other| {
            let mut other_set: Vec<AtomId> = other.atom_ids().collect();
            other_set.sort();
            other_set == mapped_set
        }) else {
            return false;
        };
        if system.attributes().charge != target.attributes().charge
            || system.attributes().unpaired_electrons != target.attributes().unpaired_electrons
            || system.attributes().constraints != target.attributes().constraints
        {
            return false;
        }
        match (system.electrons(), target.electrons()) {
            (ElectronCountsForm::Undetermined, ElectronCountsForm::Undetermined) => {}
            (ElectronCountsForm::Lit(src), ElectronCountsForm::Lit(dst))
                if src.len() == system.atom_ids().len() && dst.len() == target.atom_ids().len() =>
            {
                let dst_atoms: Vec<AtomId> = target.atom_ids().collect();
                for (pos, atom) in system.atom_ids().enumerate() {
                    let mapped_atom = AtomId(perm[atom.index()] as u32);
                    let Some(dst_pos) = dst_atoms.iter().position(|&a| a == mapped_atom) else {
                        return false;
                    };
                    if src[pos] != dst[dst_pos] {
                        return false;
                    }
                }
            }
            _ => return false,
        }
    }

    for stereo in molecule.stereo_atoms().iter() {
        let mapped_site = AtomId(perm[stereo.site_id().index()] as u32);
        let Some(target) = molecule.atom(mapped_site).stereo_atom() else {
            return false;
        };
        let mapped_ligands: Vec<StereoLigand> = stereo
            .ligand_ids()
            .iter()
            .map(|lig| StereoLigand::new(AtomId(perm[lig.atom_id.index()] as u32), lig.kind))
            .collect();
        let Some(action) = Permutation::between(&mapped_ligands, target.ligand_ids()) else {
            return false;
        };
        if stereo.attributes().clone().reframe_by(&action).as_ref() != Some(target.attributes()) {
            return false;
        }
    }

    for stereo in molecule.stereo_bonds().iter() {
        let [a, b] = molecule.bond(stereo.site_id()).atom_ids();
        let ma = AtomId(perm[a.index()] as u32);
        let mb = AtomId(perm[b.index()] as u32);
        let Some(mapped_bond) = molecule.bonds().iter().find(|other| {
            let [x, y] = other.atom_ids();
            (x == ma && y == mb) || (x == mb && y == ma)
        }) else {
            return false;
        };
        let Some(target) = mapped_bond.stereo_bond() else {
            return false;
        };
        let mapped_ligands: Vec<StereoLigand> = stereo
            .ligand_ids()
            .iter()
            .map(|lig| StereoLigand::new(AtomId(perm[lig.atom_id.index()] as u32), lig.kind))
            .collect();
        let Some(action) = Permutation::between(&mapped_ligands, target.ligand_ids()) else {
            return false;
        };
        if stereo.attributes().clone().reframe_by(&action).as_ref() != Some(target.attributes()) {
            return false;
        }
    }

    true
}

fn exact_orbit_partition(
    embeddings: &[Vec<AtomId>],
    q_group: &[Vec<usize>],
    h_group: &[Vec<usize>],
) -> BTreeSet<Vec<Vec<AtomId>>> {
    let mut orbits = BTreeSet::new();
    for embedding in embeddings {
        let mut orbit = BTreeSet::new();
        for q_perm in q_group {
            for h_perm in h_group {
                let mapped: Vec<AtomId> = q_perm
                    .iter()
                    .map(|&q| AtomId(h_perm[embedding[q].index()] as u32))
                    .collect();
                orbit.insert(mapped);
            }
        }
        orbits.insert(orbit.into_iter().collect());
    }
    orbits
}

/// Strategy emitting small (`n = 3..=5`) normalized molecules with repeated elements, optional
/// aromatic rings, and optional tetrahedral/cis-trans stereo so symmetric scaffolds occur
/// frequently enough for exhaustive `(G_Q x G_H)`-orbit verification.
fn small_symmetric_molecule_strategy() -> impl Strategy<Value = Molecule> {
    (3usize..=5usize).prop_flat_map(|n| {
        let elems = prop::collection::vec(prop_oneof![Just(Element::C), Just(Element::N)], n);
        let pair_count = n * (n - 1) / 2;
        let edge_orders =
            prop::collection::vec(prop_oneof![Just(0u8), Just(1u8), Just(2u8)], pair_count);
        let add_aromatic = any::<bool>();
        let aromatic_uniform = any::<bool>();
        let add_tetra = any::<bool>();
        let tetra_coset = 0u32..=1u32;
        (
            elems,
            edge_orders,
            add_aromatic,
            aromatic_uniform,
            add_tetra,
            tetra_coset,
        )
            .prop_map(
                move |(
                    elems,
                    edge_orders,
                    add_aromatic,
                    aromatic_uniform,
                    add_tetra,
                    tetra_coset,
                )| {
                    let atoms: Vec<AtomForm> =
                        elems.into_iter().map(AtomForm::from_element).collect();
                    let mut bonds = Vec::new();
                    let mut deg = vec![0usize; n];
                    let mut adj = vec![Vec::new(); n];
                    let mut idx = 0usize;
                    for u in 0..n {
                        for v in (u + 1)..n {
                            let order = edge_orders[idx];
                            idx += 1;
                            if order > 0 {
                                bonds.push((
                                    AtomId(u as u32),
                                    AtomId(v as u32),
                                    BondForm::from_order(order),
                                ));
                                deg[u] += 1;
                                deg[v] += 1;
                                adj[u].push(AtomId(v as u32));
                                adj[v].push(AtomId(u as u32));
                            }
                        }
                    }

                    let mut aromatic = Vec::new();
                    if add_aromatic {
                        let electrons = if aromatic_uniform {
                            ElectronCountsForm::Lit(vec![1; 3])
                        } else {
                            ElectronCountsForm::Lit(vec![2, 1, 1])
                        };
                        aromatic.push((
                            vec![AtomId(0), AtomId(1), AtomId(2)],
                            AromaticSystemForm::new(electrons),
                        ));
                    }

                    let mut stereo_atoms = Vec::new();
                    if add_tetra {
                        if let Some(center) = (0..n).find(|&v| deg[v] == 4) {
                            let ligands = adj[center]
                                .iter()
                                .map(|&nb| StereoLigand::new(nb, StereoLigandKind::Atom))
                                .collect();
                            stereo_atoms.push((
                                AtomId(center as u32),
                                ligands,
                                StereoAtomForm::new(StereoKind::Tetrahedral, tetra_coset),
                            ));
                        }
                    }

                    Molecule::from_entries(MoleculeEntries {
                        atoms,
                        bonds,
                        aromatic,
                        stereo_atoms,
                        ..Default::default()
                    })
                },
            )
    })
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

#[rstest]
#[case::tetrahedral_methane_a4(
    mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
    mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
    (12, 1, 1, 1)
)]
#[case::chiral_tetrahedral_reframed(
    mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
    mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
    (1, 1, 1, 1)
)]
#[case::rr_tartaric_half_pattern(
    mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0" "C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [3 4 "1"] [3 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"} {:site 3 :ligands [4 5 0 [:h 3]] :attrs "Th0"}]}"#),
    mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    (2, 2, 1, 1)
)]
#[case::meso_tartaric_half_pattern(
    mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0" "C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [3 4 "1"] [3 5 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"} {:site 3 :ligands [4 5 0 [:h 3]] :attrs "Th1"}]}"#),
    mol_dsl!(r#"{:atoms ["C #h1" "O #h1" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    (2, 2, 2, 2)
)]
#[case::cis_trans_symmetric_alkene_self(
    mol_dsl!(r#"{:atoms ["C" "C" "F" "Cl" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "F" "Cl" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
    (2, 1, 1, 1)
)]
#[case::cis_trans_single_end_swap_excluded(
    mol_dsl!(r#"{:atoms ["C" "C" "F" "F" "Cl" "Br"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]] :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
    mol_dsl!(r#"{:atoms ["C" "F"] :bonds [[0 1 "1"]]}"#),
    (2, 2, 2, 2)
)]
fn test_substructure_stereo_symmetry(
    #[case] host: Molecule,
    #[case] pattern: Molecule,
    #[case] expected_counts: (u64, u64, u64, u64),
) {
    let prepared_host = PreparedSubstructureHost::new(&host);
    let q_group = exhaustive_proper_automorphisms(&pattern);
    let h_group = exhaustive_proper_automorphisms(&host);
    let trivial_q = vec![(0..pattern.atoms().count()).collect::<Vec<_>>()];
    let trivial_h = vec![(0..host.atoms().count()).collect::<Vec<_>>()];
    let full = sorted_matches(&pattern, &host, GraphAndOverlays, BitsetRi);

    let cases = [
        (false, false, expected_counts.0, &trivial_q, &trivial_h),
        (true, false, expected_counts.1, &q_group, &trivial_h),
        (false, true, expected_counts.2, &trivial_q, &h_group),
        (true, true, expected_counts.3, &q_group, &h_group),
    ];

    for (red_q, red_h, expected, active_q, active_h) in cases {
        let expected_orbits = exact_orbit_partition(&full, active_q, active_h);
        assert_eq!(expected_orbits.len() as u64, expected);
        for strategy in STRATEGIES {
            let mut baseline_reps: Option<Vec<Vec<AtomId>>> = None;
            for auto_alg in AUTO_ALGS {
                for subiso in SUBISO {
                    let config = SubstructureMatchConfig {
                        match_algorithm: strategy,
                        subgraph_isomorphism_algorithm: subiso,
                        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                        automorphism_algorithm: auto_alg,
                        reduce_pattern_symmetry: red_q,
                        reduce_host_symmetry: red_h,
                    };
                    let direct =
                        sorted_atom_embeddings(&pattern.substructure_matches(&host, config));
                    let prepared = PreparedSubstructurePattern::new(&pattern, config);
                    let via_host =
                        sorted_atom_embeddings(&prepared.matches_in_host(&prepared_host));

                    assert_eq!(direct, via_host);
                    assert_eq!(direct.len() as u64, expected);
                    assert_eq!(pattern.count_substructure_matches(&host, config), expected);
                    assert_eq!(prepared.count_matches_in_host(&prepared_host), expected);
                    assert_eq!(pattern.has_substructure_match(&host, config), expected > 0);
                    assert_eq!(
                        exact_orbit_partition(&direct, active_q, active_h),
                        expected_orbits
                    );
                    if let Some(expected_reps) = &baseline_reps {
                        assert_eq!(&direct, expected_reps);
                    } else {
                        baseline_reps = Some(direct);
                    }
                }
            }
        }
    }
}

proptest! {
    #[test]
    fn test_substructure_cross_validation(
        host in molecule_strategy(),
        pattern in molecule_strategy(),
    ) {
        let prepared_host = PreparedSubstructureHost::new(&host);
        let baseline = sorted_matches(&pattern, &host, GraphAndOverlays, BitsetRi);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let direct = sorted_atom_embeddings(&pattern.substructure_matches(&host, config));
                let prepared = PreparedSubstructurePattern::new(&pattern, config);
                let via_host = sorted_atom_embeddings(&prepared.matches_in_host(&prepared_host));
                let mut embedded = Vec::new();
                let _: ControlFlow<()> =
                    pattern.visit_substructure_embeddings(&host, config, |emb| {
                        embedded.push(emb.to_vec());
                        ControlFlow::Continue(())
                    });
                embedded.sort();

                prop_assert_eq!(&direct, &baseline, "{:?}/{:?}", strategy, subiso);
                prop_assert_eq!(&via_host, &baseline, "{:?}/{:?}", strategy, subiso);
                prop_assert_eq!(&embedded, &baseline, "{:?}/{:?}", strategy, subiso);
                prop_assert_eq!(
                    pattern.has_substructure_match(&host, config),
                    !baseline.is_empty()
                );
                prop_assert_eq!(
                    prepared.has_match_in_host(&prepared_host),
                    !baseline.is_empty()
                );
                prop_assert_eq!(
                    pattern.count_substructure_matches(&host, config),
                    baseline.len() as u64
                );
                prop_assert_eq!(
                    prepared.count_matches_in_host(&prepared_host),
                    baseline.len() as u64
                );
                let first = pattern.first_substructure_match(&host, config);
                prop_assert_eq!(first.is_some(), !baseline.is_empty());
                if let Some(first_corr) = first {
                    let atoms: Vec<AtomId> = first_corr
                        .atoms()
                        .matched_pairs()
                        .iter()
                        .map(|&(_, h)| h)
                        .collect();
                    prop_assert!(baseline.contains(&atoms));
                }
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
        let prepared_host = PreparedSubstructureHost::new(&host);
        let baseline = sorted_matches(&pattern, &host, GraphAndOverlays, BitsetRi);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let direct = sorted_atom_embeddings(&pattern.substructure_matches(&host, config));
                let prepared = PreparedSubstructurePattern::new(&pattern, config);
                let via_host = sorted_atom_embeddings(&prepared.matches_in_host(&prepared_host));
                prop_assert_eq!(&direct, &baseline, "{:?}/{:?}", strategy, subiso);
                prop_assert_eq!(&via_host, &baseline, "{:?}/{:?}", strategy, subiso);
                prop_assert_eq!(
                    pattern.count_substructure_matches(&host, config),
                    baseline.len() as u64
                );
                prop_assert_eq!(
                    prepared.count_matches_in_host(&prepared_host),
                    baseline.len() as u64
                );
                prop_assert_eq!(
                    pattern.has_substructure_match(&host, config),
                    !baseline.is_empty()
                );
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_substructure_symmetry_modes_planted(
        (host, subset) in molecule_strategy()
            .prop_filter("non-empty", |m| m.atoms().count() > 0)
            .prop_flat_map(|m| {
                let n = m.atoms().count();
                let ids: Vec<AtomId> = (0..n as u32).map(AtomId).collect();
                (Just(m), prop::sample::subsequence(ids, 1..=n))
            }),
    ) {
        let pattern = host.extract(&host.induced_subgraph(&subset));
        let prepared_host = PreparedSubstructureHost::new(&host);
        let full_set: BTreeSet<Vec<AtomId>> =
            sorted_matches(&pattern, &host, GraphAndOverlays, BitsetRi)
                .into_iter()
                .collect();

        for (red_q, red_h) in SYMMETRY_MODES {
            let mut strategy_count: Option<u64> = None;
            for strategy in STRATEGIES {
                let mut baseline_reps: Option<Vec<Vec<AtomId>>> = None;
                for auto_alg in AUTO_ALGS {
                    for subiso in SUBISO {
                        let config = SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm: auto_alg,
                            reduce_pattern_symmetry: red_q,
                            reduce_host_symmetry: red_h,
                        };
                        let prepared = PreparedSubstructurePattern::new(&pattern, config);
                        let reps = sorted_atom_embeddings(&prepared.matches_in_host(&prepared_host));
                        let count = prepared.count_matches_in_host(&prepared_host);
                        let direct_count = pattern.count_substructure_matches(&host, config);
                        let has = prepared.has_match_in_host(&prepared_host);

                        prop_assert_eq!(count, reps.len() as u64);
                        prop_assert_eq!(direct_count, count);
                        prop_assert_eq!(has, !full_set.is_empty());
                        for rep in &reps {
                            prop_assert!(full_set.contains(rep));
                        }
                        if let Some(expected_reps) = &baseline_reps {
                            prop_assert_eq!(
                                &reps,
                                expected_reps,
                                "{:?}/{:?}/{:?} ({},{})",
                                strategy,
                                auto_alg,
                                subiso,
                                red_q,
                                red_h
                            );
                        } else {
                            baseline_reps = Some(reps);
                        }
                        if let Some(expected_count) = strategy_count {
                            prop_assert_eq!(
                                count,
                                expected_count,
                                "{:?}/{:?}/{:?} ({},{})",
                                strategy,
                                auto_alg,
                                subiso,
                                red_q,
                                red_h
                            );
                        } else {
                            strategy_count = Some(count);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_substructure_symmetry_exact_orbits(
        (host, subset) in small_symmetric_molecule_strategy().prop_flat_map(|m| {
            let n = m.atoms().count();
            let ids: Vec<AtomId> = (0..n as u32).map(AtomId).collect();
            (Just(m), prop::sample::subsequence(ids, 1..=n))
        }),
    ) {
        let pattern = host.extract(&host.induced_subgraph(&subset));
        let prepared_host = PreparedSubstructureHost::new(&host);
        let full = sorted_matches(&pattern, &host, GraphAndOverlays, BitsetRi);
        let q_group = exhaustive_proper_automorphisms(&pattern);
        let h_group = exhaustive_proper_automorphisms(&host);
        let trivial_q = vec![(0..pattern.atoms().count()).collect::<Vec<_>>()];
        let trivial_h = vec![(0..host.atoms().count()).collect::<Vec<_>>()];

        let modes = [
            (false, false, &trivial_q, &trivial_h),
            (true, false, &q_group, &trivial_h),
            (false, true, &trivial_q, &h_group),
            (true, true, &q_group, &h_group),
        ];

        for (red_q, red_h, active_q, active_h) in modes {
            let expected_orbits = exact_orbit_partition(&full, active_q, active_h);
            let expected_count = expected_orbits.len() as u64;

            for strategy in STRATEGIES {
                let mut baseline_reps: Option<Vec<Vec<AtomId>>> = None;
                for auto_alg in AUTO_ALGS {
                    for subiso in SUBISO {
                        let config = SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm: auto_alg,
                            reduce_pattern_symmetry: red_q,
                            reduce_host_symmetry: red_h,
                        };
                        let prepared = PreparedSubstructurePattern::new(&pattern, config);
                        let reps = sorted_atom_embeddings(&prepared.matches_in_host(&prepared_host));
                        let count = prepared.count_matches_in_host(&prepared_host);

                        prop_assert_eq!(
                            count,
                            expected_count,
                            "{:?}/{:?}/{:?} ({},{})",
                            strategy,
                            auto_alg,
                            subiso,
                            red_q,
                            red_h
                        );
                        prop_assert_eq!(reps.len() as u64, expected_count);
                        prop_assert_eq!(
                            &exact_orbit_partition(&reps, active_q, active_h),
                            &expected_orbits,
                            "{:?}/{:?}/{:?} ({},{})",
                            strategy,
                            auto_alg,
                            subiso,
                            red_q,
                            red_h
                        );
                        if let Some(expected_reps) = &baseline_reps {
                            prop_assert_eq!(&reps, expected_reps);
                        } else {
                            baseline_reps = Some(reps);
                        }
                    }
                }
            }
        }
    }
}
