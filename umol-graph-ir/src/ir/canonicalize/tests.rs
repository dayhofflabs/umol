//! Unit and bounded-exhaustive tests for aggregate canonicalization.
//!
//! These tests exercise exact typed-key ordering, adapter projection, search invariance,
//! stereo-frame transport, canonical representatives, and aggregate failure boundaries.

use std::cmp::Ordering;
use std::{array, iter};

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_perm::{Orientation, Permutation};

use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomFieldChange, AtomForm, AtomId,
    BondForm, BondId, BooleanForm, Constraint, ConstraintDelta, Constraints, DativeBondForm,
    DativeBondId, Entity, FluxionalityForm, IncidenceLevel, LigandPermutation, LigandSymmetryForm,
    MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, NoncovalentBondId, OrientedLigandPermutation, ReactionSpanEntries,
    StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondConstraintForm,
    StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand,
    StereoLigandPair, StereoTerm, Stereogenicity, StereogenicityForm, Topicity, TopicityForm,
    TopicityRelationForm,
};

fn node_branch_order(
    _adapter: &AutomorphismAdapter,
    _partition: &OrderedPartition,
    _algorithm: AutomorphismAlgorithm,
    _automorphisms: Option<&ProjectedAutomorphismOutput>,
    candidates: &mut [NodeId],
) -> bool {
    candidates.sort_unstable();
    false
}

fn reverse_node_branch_order(
    _adapter: &AutomorphismAdapter,
    _partition: &OrderedPartition,
    _algorithm: AutomorphismAlgorithm,
    _automorphisms: Option<&ProjectedAutomorphismOutput>,
    candidates: &mut [NodeId],
) -> bool {
    candidates.sort_unstable_by(|lhs, rhs| rhs.cmp(lhs));
    false
}

impl OrderedPartition {
    fn fixed_entity_prefix(&self, entity_count: usize) -> Vec<NodeId> {
        self.cells
            .iter()
            .take_while(|cell| cell.len() == 1)
            .flatten()
            .copied()
            .take_while(|node| node.index() < entity_count)
            .collect()
    }
}

impl AutomorphismAdapter {
    fn automorphisms(&self, algorithm: AutomorphismAlgorithm) -> ProjectedAutomorphismOutput {
        let output = self
            .graph()
            .automorphisms(|node| self.class(node), algorithm);
        self.project_automorphisms(&output)
    }
}

fn adapter_entity_blocks(incidence_graph: &IncidenceGraph) -> Vec<Vec<NodeId>> {
    let mut blocks = Vec::<Vec<NodeId>>::new();
    let mut previous_kind = None;
    for node in incidence_graph.graph().node_ids() {
        let kind = incidence_graph.entity(node).kind();
        if previous_kind != Some(kind) {
            blocks.push(Vec::new());
            previous_kind = Some(kind);
        }
        blocks
            .last_mut()
            .expect("current entity block is present")
            .push(node);
    }
    blocks
}

fn exhaustive_minimum<K, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    mut cells: Vec<Vec<NodeId>>,
    leaf_candidate: &LeafCandidate,
) -> CanonicalCandidate<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    fn visit_cells<K, LeafCandidate>(
        cells: &mut [Vec<NodeId>],
        cell_index: usize,
        order: &mut Vec<NodeId>,
        leaf_candidate: &LeafCandidate,
        best: &mut Option<CanonicalCandidate<K>>,
    ) where
        K: Ord,
        LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    {
        fn visit_permutations<K, LeafCandidate>(
            cells: &mut [Vec<NodeId>],
            cell_index: usize,
            position: usize,
            order: &mut Vec<NodeId>,
            leaf_candidate: &LeafCandidate,
            best: &mut Option<CanonicalCandidate<K>>,
        ) where
            K: Ord,
            LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
        {
            if position == cells[cell_index].len() {
                let old_len = order.len();
                order.extend_from_slice(&cells[cell_index]);
                visit_cells(cells, cell_index + 1, order, leaf_candidate, best);
                order.truncate(old_len);
                return;
            }

            for next in position..cells[cell_index].len() {
                cells[cell_index].swap(position, next);
                visit_permutations(cells, cell_index, position + 1, order, leaf_candidate, best);
                cells[cell_index].swap(position, next);
            }
        }

        if cell_index == cells.len() {
            let candidate = leaf_candidate(order);
            if best.as_ref().is_none_or(|best| candidate.key < best.key) {
                *best = Some(candidate);
            }
            return;
        }

        visit_permutations(cells, cell_index, 0, order, leaf_candidate, best);
    }

    let mut best = None;
    visit_cells(
        &mut cells,
        0,
        &mut Vec::with_capacity(adapter.source_node_count),
        leaf_candidate,
        &mut best,
    );
    best.expect("every finite partition has an entity ordering")
}

fn initial_classes(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<InitialClasses, Contradiction> {
    let (entity_keys, incidence_keys) = initial_class_keys(molecule, incidence_graph)?;
    Ok(rank_initial_classes(&entity_keys, &incidence_keys))
}

fn topology_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(topology_candidate(molecule, incidence_graph, order)?.key)
}

fn constitution_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(constitution_candidate(molecule, incidence_graph, order)?.key)
}

fn structure_comparison_key(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    order: &[NodeId],
) -> Result<CanonicalComparisonKey, Contradiction> {
    Ok(structure_candidate(molecule, incidence_graph, order)?.key)
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MinimumStereoFrames {
    configuration: StereoConfigurationForm,
    permutations: Vec<Permutation>,
}

fn minimum_kinded_stereo_frames(
    configuration: &StereoConfigurationForm,
    before: &[StereoLigand],
    after: &[StereoLigand],
) -> Result<Option<MinimumStereoFrames>, Contradiction> {
    let Some(kind) = configuration.kind() else {
        return Ok(None);
    };
    if before.len() != kind.degree() || after.len() != kind.degree() {
        return Ok(None);
    }

    let mut minimum: Option<StereoConfigurationForm> = None;
    let mut permutations = Vec::new();
    for permutation in Permutation::between_all(before, after)
        .into_iter()
        .filter(|permutation| kind.class_key().space().reindex(0, *permutation).is_some())
    {
        let candidate = configuration.apply(permutation).normalize()?;
        match minimum.as_ref().map(|value| candidate.cmp(value)) {
            None | Some(Ordering::Less) => {
                minimum = Some(candidate);
                permutations.clear();
                permutations.push(permutation);
            }
            Some(Ordering::Equal) => permutations.push(permutation),
            Some(Ordering::Greater) => {}
        }
    }

    Ok(minimum.map(|configuration| MinimumStereoFrames {
        configuration,
        permutations,
    }))
}

#[fixture]
fn initial_class_molecule() -> Molecule {
    let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(1),
        ArithExpr::Lit(2),
    ])));
    let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(0),
        ArithExpr::Lit(1),
    ])));
    let bond_ligands = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ];

    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::C)
                .with_charge(normalized_three)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
            AtomForm::from_element(Element::O).with_charge(3_i64),
            AtomForm::from_element(Element::C).with_charge(4_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
            (
                vec![AtomId(1)],
                AtomId(4),
                DativeBondForm::new(normalized_one),
            ),
            (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                AtomId(0),
                AtomId(5),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                AtomId(1),
                AtomId(5),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                AtomId(2),
                AtomId(5),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                bond_ligands.clone(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                bond_ligands,
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

#[fixture]
fn canonicalize_context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

#[fixture]
fn stereo_atom_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn stereo_bond_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn symmetric_stereo_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn meso_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
        ],
        ..Default::default()
    })
}

#[fixture]
fn repeated_ligand_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn para_stereo_canonicalization_molecule() -> Molecule {
    let outer_ligands = vec![
        StereoLigand::new(AtomId(10), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(11), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(12), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(13), StereoLigandKind::Atom),
    ];

    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(AtomForm::from_element)
        .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [2, 3, 4, 5]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                [6, 8, 7, 9]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(2),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(3),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                AtomId(4),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(5),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
            (
                AtomId(6),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(7),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                AtomId(8),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(9),
                outer_ligands,
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
        ],
        ..Default::default()
    })
}

fn selected_structure_key(molecule: &Molecule) -> CanonicalComparisonKey {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    structure_comparison_key(molecule, &incidence_graph, &order)
        .expect("canonical molecule has a structure comparison key")
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::axial(StereoKind::Axial)]
#[case::square_planar(StereoKind::SquarePlanar)]
#[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal)]
#[case::octahedral(StereoKind::Octahedral)]
fn test_reframe_stereo_atom(#[case] kind: StereoKind) {
    let degree = kind.degree();
    let swap = |first, second| {
        let mut image = (0..degree).collect::<Vec<_>>();
        image.swap(first, second);
        Permutation::from_image(&image)
    };
    let (frame, next_frame, expected_permutation, expected_pair) = if kind == StereoKind::Axial {
        (
            Permutation::from_image(&[2, 3, 0, 1]),
            Permutation::from_image(&[1, 0, 2, 3]),
            LigandPermutation(swap(2, 3)),
            StereoLigandPair::new(2usize.into(), 3usize.into()),
        )
    } else {
        (
            Permutation::from_image(&(1..degree).chain(iter::once(0)).collect::<Vec<_>>()),
            swap(1, 2),
            LigandPermutation(swap(0, degree - 1)),
            StereoLigandPair::new(0usize.into(), (degree - 1).into()),
        )
    };
    let source_permutation = LigandPermutation(swap(0, 1));
    let source_pair = StereoLigandPair::new(0usize.into(), 1usize.into());
    let source_constraints = vec![
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: source_permutation,
                orientation: Orientation::Improper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: source_permutation,
            active: BooleanForm::Lit(true),
        }),
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: source_pair,
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        }),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ];
    let expected_constraints = vec![
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: expected_permutation,
                orientation: Orientation::Improper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
            permutation: expected_permutation,
            active: BooleanForm::Lit(true),
        }),
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: expected_pair,
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        }),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ];
    let global_constraints = |constraints: &[StereoAtomConstraintForm]| {
        Constraints::from(Constraint::And(vec![
            Constraint::Not(Box::new(Constraint::StereoAtom(
                StereoAtomId(0),
                kind,
                constraints[0].clone(),
            ))),
            Constraint::Or(
                constraints[1..]
                    .iter()
                    .cloned()
                    .map(|constraint| Constraint::StereoAtom(StereoAtomId(0), kind, constraint))
                    .collect(),
            ),
        ]))
    };
    let atoms = (0..=degree)
        .map(|_| AtomForm::from_element(Element::C))
        .collect::<Vec<_>>();
    let ligands = (1..=degree)
        .map(|atom| StereoLigand::new(AtomId(atom as u32), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let source_form = StereoAtomForm {
        configuration: StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0)),
        constraints: source_constraints.clone().into(),
    };
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: atoms.clone(),
        stereo_atoms: vec![(AtomId(0), ligands.clone(), source_form.clone())],
        constraints: global_constraints(&source_constraints),
        ..Default::default()
    });
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms,
        stereo_atoms: vec![(
            AtomId(0),
            frame.act(&ligands),
            StereoAtomForm {
                configuration: source_form.configuration.apply(frame),
                constraints: expected_constraints.clone().into(),
            },
        )],
        constraints: global_constraints(&expected_constraints),
        ..Default::default()
    });

    let reframed = reframe_stereo_atom(&source, StereoAtomId(0), frame);
    assert_eq!(reframed, expected);
    assert_eq!(
        reframe_stereo_atom(&reframed, StereoAtomId(0), frame.inverse()),
        source
    );
    assert_eq!(
        reframe_stereo_atom(&reframed, StereoAtomId(0), next_frame),
        reframe_stereo_atom(&source, StereoAtomId(0), frame.compose(next_frame))
    );
}

#[rstest]
fn test_reframe_stereo_bond() {
    let degree = StereoKind::CisTrans.degree();
    let swap = |first, second| {
        let mut image = (0..degree).collect::<Vec<_>>();
        image.swap(first, second);
        Permutation::from_image(&image)
    };
    let frame = Permutation::from_image(&[2, 3, 0, 1]);
    let next_frame = Permutation::from_image(&[1, 0, 2, 3]);
    let source_permutation = LigandPermutation(swap(0, 1));
    let expected_permutation = LigandPermutation(swap(2, 3));
    let source_pair = StereoLigandPair::new(0usize.into(), 1usize.into());
    let expected_pair = StereoLigandPair::new(2usize.into(), 3usize.into());
    let source_constraints = vec![
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: source_permutation,
                orientation: Orientation::Proper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        StereoBondConstraintForm::Fluxionality(FluxionalityForm {
            permutation: source_permutation,
            active: BooleanForm::Lit(true),
        }),
        StereoBondConstraintForm::Topicity(TopicityForm {
            pair: source_pair,
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        }),
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ];
    let expected_constraints = vec![
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: expected_permutation,
                orientation: Orientation::Proper,
            },
            invariant: BooleanForm::Lit(true),
        }),
        StereoBondConstraintForm::Fluxionality(FluxionalityForm {
            permutation: expected_permutation,
            active: BooleanForm::Lit(true),
        }),
        StereoBondConstraintForm::Topicity(TopicityForm {
            pair: expected_pair,
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        }),
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ];
    let global_constraints = |constraints: &[StereoBondConstraintForm]| {
        Constraints::from(Constraint::And(vec![
            Constraint::Not(Box::new(Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                constraints[0].clone(),
            ))),
            Constraint::Or(
                constraints[1..]
                    .iter()
                    .cloned()
                    .map(|constraint| {
                        Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, constraint)
                    })
                    .collect(),
            ),
        ]))
    };
    let atoms = (0..6)
        .map(|_| AtomForm::from_element(Element::C))
        .collect::<Vec<_>>();
    let ligands = (2..6)
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let source_form = StereoBondForm {
        configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        constraints: source_constraints.clone().into(),
    };
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: atoms.clone(),
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        stereo_bonds: vec![(BondId(0), ligands.clone(), source_form.clone())],
        constraints: global_constraints(&source_constraints),
        ..Default::default()
    });
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        stereo_bonds: vec![(
            BondId(0),
            frame.act(&ligands),
            StereoBondForm {
                configuration: source_form.configuration.apply(frame),
                constraints: expected_constraints.clone().into(),
            },
        )],
        constraints: global_constraints(&expected_constraints),
        ..Default::default()
    });

    let reframed = reframe_stereo_bond(&source, StereoBondId(0), frame);
    assert_eq!(reframed, expected);
    assert_eq!(
        reframe_stereo_bond(&reframed, StereoBondId(0), frame.inverse()),
        source
    );
    assert_eq!(
        reframe_stereo_bond(&reframed, StereoBondId(0), next_frame),
        reframe_stereo_bond(&source, StereoBondId(0), frame.compose(next_frame))
    );
}

#[rstest]
#[case::literal(
        StereoCoset::Lit(0),
        StereoCoset::Lit(0),
        vec![Permutation::from_image(&[1, 2, 0, 3])]
    )]
#[case::undetermined(
        StereoCoset::Undetermined,
        StereoCoset::Undetermined,
        vec![
            Permutation::from_image(&[1, 2, 0, 3]),
            Permutation::from_image(&[2, 1, 0, 3]),
        ]
    )]
#[case::set_valued(
        StereoCoset::lit_set([0, 1]),
        StereoCoset::lit_set([0, 1]),
        vec![
            Permutation::from_image(&[1, 2, 0, 3]),
            Permutation::from_image(&[2, 1, 0, 3]),
        ]
    )]
#[case::symbolic(
        StereoCoset::term(StereoTerm::var("x")),
        StereoCoset::term(StereoTerm::var("x")),
        vec![Permutation::from_image(&[1, 2, 0, 3])]
    )]
fn test_minimum_kinded_stereo_frames(
    #[case] coset: StereoCoset,
    #[case] expected_coset: StereoCoset,
    #[case] expected_permutations: Vec<Permutation>,
) {
    let repeated = StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen);
    let before = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        repeated,
        repeated,
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
    ];
    let (after, _) = sort_ligand_frame(&before);
    let result = minimum_kinded_stereo_frames(
        &StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset),
        &before,
        &after,
    )
    .expect("fixed configurations normalize")
    .expect("the frames contain the same ligand multiset");

    assert_eq!(
        result,
        MinimumStereoFrames {
            configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, expected_coset,),
            permutations: expected_permutations,
        }
    );
}

#[rstest]
fn test_kindless_stereo_atom_frame_order() {
    let ligands = (0..7)
        .rev()
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let (sorted, order) = sort_ligand_frame(&ligands);
    let source = StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        })
        .into(),
    };
    let expected = StereoAtomForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(4usize.into(), 6usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        })
        .into(),
    };

    assert_eq!(
        sorted,
        (0..7)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>()
    );
    assert_eq!(
        reframe_stereo_atom_form_by_order(&source, &order),
        Some(expected)
    );
}

#[rstest]
fn test_kindless_stereo_bond_frame_order() {
    let ligands = (0..4)
        .rev()
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let (_, order) = sort_ligand_frame(&ligands);
    let source = StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: vec![
            StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                    orientation: Orientation::Proper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
            }),
        ]
        .into(),
    };
    let expected = StereoBondForm {
        configuration: StereoConfigurationForm::Undetermined,
        constraints: vec![
            StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: LigandPermutation(Permutation::from_image(&[0, 1, 3, 2])),
                    orientation: Orientation::Proper,
                },
                invariant: BooleanForm::Lit(true),
            }),
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(2usize.into(), 3usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
            }),
        ]
        .into(),
    };

    assert_eq!(
        reframe_stereo_bond_form_by_order(&source, &order),
        Some(expected)
    );
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
#[case::axial(StereoKind::Axial)]
#[case::square_planar(StereoKind::SquarePlanar)]
#[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal)]
#[case::octahedral(StereoKind::Octahedral)]
fn test_stereo_refinement_descriptor_frame_invariant(#[case] kind: StereoKind) {
    let degree = kind.degree();
    let frame = if matches!(kind, StereoKind::CisTrans | StereoKind::Axial) {
        Permutation::from_image(&[2, 3, 0, 1])
    } else {
        Permutation::from_image(&(1..degree).chain(iter::once(0)).collect::<Vec<_>>())
    };
    let ligands = (0..degree)
        .map(|class| (class as u32, StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let configuration = StereoConfigurationForm::kinded(kind, 0u32);

    assert_eq!(
        stereo_refinement_descriptor(7, &ligands, &configuration),
        stereo_refinement_descriptor(7, &frame.act(&ligands), &configuration.apply(frame),),
    );
}

#[rstest]
#[case::one_pass(false, 1)]
#[case::fixpoint(true, 2)]
fn test_structure_partition(
    para_stereo_canonicalization_molecule: Molecule,
    #[case] para_stereo: bool,
    #[case] expected_rounds: usize,
) {
    let molecule = para_stereo_canonicalization_molecule;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph)
        .expect("fixed molecule has initial classes");
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let (_, rounds) = structure_partition(
        &molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        para_stereo,
    )
    .expect("fixed molecule has a structure partition");

    assert_eq!(rounds, expected_rounds);
}

#[rstest]
fn test_structure_partition_no_stereo() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph)
        .expect("fixed molecule has initial classes");
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let (_, rounds) =
        structure_partition(&molecule, &incidence_graph, &adapter, &entity_keys, true)
            .expect("fixed molecule has a structure partition");

    assert_eq!(rounds, 1);
}

#[rstest]
fn test_structure_partition_distinct_stereo(stereo_atom_canonicalization_molecule: Molecule) {
    let molecule = stereo_atom_canonicalization_molecule;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph)
        .expect("fixed molecule has initial classes");
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let (_, rounds) =
        structure_partition(&molecule, &incidence_graph, &adapter, &entity_keys, true)
            .expect("fixed molecule has a structure partition");

    assert_eq!(rounds, 1);
}

#[rstest]
fn test_canonicalize_structure_para_stereo(
    para_stereo_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let correspondence = molecule_correspondence(&[
        vec![1, 0, 6, 8, 7, 9, 2, 4, 3, 5, 13, 11, 10, 12],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![1, 0, 6, 8, 7, 9, 2, 4, 3, 5],
        Vec::new(),
    ]);
    let renumbered = para_stereo_canonicalization_molecule.remap(&correspondence);
    let context = CanonicalizeContext {
        para_stereo: true,
        ..canonicalize_context
    };
    let canonical = canonicalize_structure(&para_stereo_canonicalization_molecule, &context)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        canonicalize_structure(&renumbered, &context),
        Ok(canonical.clone()),
    );
    assert_eq!(canonicalize_structure(&canonical, &context), Ok(canonical));
}

#[rstest]
fn test_structure_comparison_key(stereo_atom_canonicalization_molecule: Molecule) {
    let mut constrained_entries = molecule_entries(&stereo_atom_canonicalization_molecule);
    constrained_entries.stereo_atoms[0].2.constraints = StereoAtomConstraintForm::Stereogenicity(
        StereogenicityForm::Lit(Stereogenicity::Stereogenic),
    )
    .into();
    let constrained = Molecule::from_entries(constrained_entries);
    let incidence_graph =
        stereo_atom_canonicalization_molecule.incidence_graph(IncidenceLevel::Full);
    let constrained_incidence_graph = constrained.incidence_graph(IncidenceLevel::Full);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();

    assert_eq!(
        structure_comparison_key(
            &stereo_atom_canonicalization_molecule,
            &incidence_graph,
            &order,
        ),
        structure_comparison_key(&constrained, &constrained_incidence_graph, &order),
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_atom(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_correspondence(&[
        vec![4, 2, 0, 3, 1],
        vec![2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
        Vec::new(),
    ]);
    let renumbered = stereo_atom_canonicalization_molecule.remap(&renumbering);
    let canonical = canonicalize_structure(
        &stereo_atom_canonicalization_molecule,
        &canonicalize_context,
    )
    .expect("fixed molecule canonicalizes");

    assert_eq!(
        canonicalize_structure(&renumbered, &canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonicalize_structure(&canonical, &canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical.stereo_atoms().count(), 1);
}

#[rstest]
fn test_canonicalize_structure_configuration(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let mut opposite_entries = molecule_entries(&stereo_atom_canonicalization_molecule);
    opposite_entries.stereo_atoms[0].2.configuration =
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32);
    let opposite = Molecule::from_entries(opposite_entries);

    assert_ne!(
        canonicalize_structure(
            &stereo_atom_canonicalization_molecule,
            &canonicalize_context,
        ),
        canonicalize_structure(&opposite, &canonicalize_context),
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_bond(
    canonicalize_context: CanonicalizeContext,
    stereo_bond_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_correspondence(&[
        vec![5, 3, 1, 4, 2, 0],
        vec![4, 2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
    ]);
    let renumbered = stereo_bond_canonicalization_molecule.remap(&renumbering);
    let canonical = canonicalize_structure(
        &stereo_bond_canonicalization_molecule,
        &canonicalize_context,
    )
    .expect("fixed molecule canonicalizes");

    assert_eq!(
        canonicalize_structure(&renumbered, &canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonicalize_structure(&canonical, &canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical.stereo_bonds().count(), 1);
}

#[rstest]
fn test_canonicalize_structure_stereo_atom_constraints(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    let mut entries = molecule_entries(&stereo_atom_canonicalization_molecule);
    entries.stereo_atoms[0].2.constraints = constraint.clone().into();
    entries.constraints =
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
    let source = Molecule::from_entries(entries);
    let reframed = reframe_stereo_atom(
        &source,
        StereoAtomId(0),
        Permutation::from_image(&[1, 0, 2, 3]),
    );

    assert_eq!(
        canonicalize_structure(&reframed, &canonicalize_context),
        canonicalize_structure(&source, &canonicalize_context),
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_bond_constraints(
    canonicalize_context: CanonicalizeContext,
    stereo_bond_canonicalization_molecule: Molecule,
) {
    let constraint = StereoBondConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
    });
    let mut entries = molecule_entries(&stereo_bond_canonicalization_molecule);
    entries.stereo_bonds[0].2.constraints = constraint.clone().into();
    entries.constraints =
        Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, constraint).into();
    let source = Molecule::from_entries(entries);
    let reframed = reframe_stereo_bond(
        &source,
        StereoBondId(0),
        Permutation::from_image(&[2, 3, 0, 1]),
    );
    let left = canonicalize_structure(&source, &canonicalize_context)
        .expect("fixed molecule canonicalizes");
    let right = canonicalize_structure(&reframed, &canonicalize_context)
        .expect("reframed molecule canonicalizes");

    assert_eq!(
        selected_structure_key(&right),
        selected_structure_key(&left)
    );
}

#[rstest]
fn test_molecule_canonicalize(canonicalize_context: CanonicalizeContext) {
    let plain = AtomForm::from_element(Element::C);
    let mut constrained = plain.clone();
    constrained.constraints = AtomConstraintForm::valence(4).into();
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![plain.clone(), constrained.clone()],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(1), AtomId(1)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    });
    let renumbered = source.remap(&molecule_correspondence(&[
        vec![1, 0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ]));
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![constrained, plain],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    });

    assert_eq!(
        source.clone().canonicalize(&canonicalize_context),
        Ok(expected.clone()),
    );
    assert_eq!(
        renumbered.canonicalize(&canonicalize_context),
        Ok(expected.clone()),
    );
    assert_eq!(
        expected.clone().canonicalize(&canonicalize_context),
        Ok(expected.clone()),
    );
}

#[rstest]
#[case::feature_free_connected(
    r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 9 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [8 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#
)]
#[case::feature_free_disconnected(
    r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 8 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#
)]
#[case::symmetry_heavy_radicals(
    r#"
{:atoms ["H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[6 7 "3#c0#u0#s"]]}
"#
)]
fn test_molecule_canonicalize_retained(
    canonicalize_context: CanonicalizeContext,
    #[case] source: &str,
) {
    let source = Molecule::from_edn_str(source).expect("retained molecule parses");
    let expected = source.clone();

    for level in [
        CanonicalizeLevel::Topology,
        CanonicalizeLevel::Constitution,
        CanonicalizeLevel::Structure,
        CanonicalizeLevel::Full,
    ] {
        let canonical = source
            .clone()
            .canonicalize_by(level, &canonicalize_context)
            .expect("retained molecule canonicalizes");

        assert_eq!(canonical, expected);
        assert_eq!(
            canonical_key_by(&canonical, level, &canonicalize_context),
            canonical_key_by(&expected, level, &canonicalize_context),
        );
    }

    let (canonical, correspondence) = source
        .clone()
        .canonicalize_with_correspondence(&canonicalize_context)
        .expect("retained molecule canonicalizes with a correspondence");
    let transported = source.remap(&correspondence);

    assert_eq!(canonical, expected);
    assert_eq!(transported, expected);
    assert!(source.equiv_under(&canonical, &correspondence));
}

#[rstest]
fn test_molecule_canonicalize_stereo_frame(
    repeated_ligand_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    let mut entries = molecule_entries(&repeated_ligand_canonicalization_molecule);
    entries.stereo_atoms[0].2.constraints = constraint.clone().into();
    entries.constraints =
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
    let source = Molecule::from_entries(entries);
    let reframed = reframe_stereo_atom(
        &source,
        StereoAtomId(0),
        Permutation::from_image(&[0, 1, 3, 2]),
    );
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        reframed.canonicalize(&canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical.clone().canonicalize(&canonicalize_context),
        Ok(canonical),
    );
}

#[rstest]
fn test_molecule_canonicalize_contradiction(canonicalize_context: CanonicalizeContext) {
    let mut atom = AtomForm::from_element(Element::C);
    atom.constraints = AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())).into();
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom],
        ..Default::default()
    });

    assert_eq!(
        molecule.canonicalize(&canonicalize_context),
        Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
fn test_molecule_canonicalize_integrity_error(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let mut malformed = stereo_atom_canonicalization_molecule;
    malformed
        .stereo_atom_mut(StereoAtomId(0))
        .attributes
        .configuration = StereoConfigurationForm::kinded(StereoKind::Octahedral, 0u32);

    assert_eq!(
        malformed.canonicalize(&canonicalize_context),
        Err(MoleculeCanonicalizeError::Integrity(
            MoleculeIntegrityError::StereoLigandArity {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                kind: StereoKind::Octahedral,
                expected: 6,
                actual: 4,
            },
        )),
    );
}

#[rstest]
#[case::topology(CanonicalizeLevel::Topology)]
#[case::constitution(CanonicalizeLevel::Constitution)]
#[case::structure(CanonicalizeLevel::Structure)]
#[case::full(CanonicalizeLevel::Full)]
fn test_molecule_canonicalize_by(
    initial_class_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] level: CanonicalizeLevel,
) {
    let canonical = initial_class_molecule
        .clone()
        .canonicalize_by(level, &canonicalize_context)
        .expect("fixed molecule canonicalizes at every level");

    assert_eq!(
        canonical
            .clone()
            .canonicalize_by(level, &canonicalize_context),
        Ok(canonical),
    );
    if level == CanonicalizeLevel::Full {
        assert_eq!(
            initial_class_molecule
                .clone()
                .canonicalize_by(level, &canonicalize_context),
            initial_class_molecule.canonicalize(&canonicalize_context),
        );
    }
}

#[rstest]
fn test_molecule_canonical_eq(
    initial_class_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let renumbered = initial_class_molecule.remap(&reverse_correspondence(&initial_class_molecule));

    assert!(initial_class_molecule.canonical_eq(&renumbered, &canonicalize_context));
}

#[rstest]
fn test_molecule_canonical_eq_contradiction(canonicalize_context: CanonicalizeContext) {
    let mut left_contradiction = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    left_contradiction.atom_mut(AtomId(0)).attributes.charge = NumForm::lit_set(Vec::<i64>::new());
    let mut right_contradiction = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::N)],
        ..Default::default()
    });
    right_contradiction.atom_mut(AtomId(0)).attributes.charge = NumForm::lit_set(Vec::<i64>::new());
    let valid = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });

    assert!(left_contradiction.canonical_eq(&right_contradiction, &canonicalize_context));
    assert!(!left_contradiction.canonical_eq(&valid, &canonicalize_context));
}

#[rstest]
fn test_molecule_canonical_eq_integrity_error(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let mut malformed_left = stereo_atom_canonicalization_molecule.clone();
    malformed_left
        .stereo_atom_mut(StereoAtomId(0))
        .attributes
        .configuration = StereoConfigurationForm::kinded(StereoKind::Octahedral, 0u32);
    let mut malformed_right = stereo_atom_canonicalization_molecule;
    malformed_right
        .stereo_atom_mut(StereoAtomId(0))
        .attributes
        .configuration = StereoConfigurationForm::kinded(StereoKind::TrigonalBipyramidal, 0u32);

    assert!(malformed_left.canonical_eq(&malformed_left, &canonicalize_context));
    assert!(!malformed_left.canonical_eq(&malformed_right, &canonicalize_context));
}

#[rstest]
fn test_molecule_canonical_eq_by_topology(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    });

    assert!(left.canonical_eq_by(&right, CanonicalizeLevel::Topology, &canonicalize_context,));
    assert!(!left.canonical_eq_by(
        &right,
        CanonicalizeLevel::Constitution,
        &canonicalize_context,
    ));
}

#[rstest]
fn test_molecule_canonical_hash_by_topology(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    });

    assert_eq!(
        left.canonical_hash_by(CanonicalizeLevel::Topology, &canonicalize_context),
        right.canonical_hash_by(CanonicalizeLevel::Topology, &canonicalize_context),
    );
}

#[rstest]
fn test_molecule_canonical_eq_by_constitution(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let mut entries = molecule_entries(&stereo_atom_canonicalization_molecule);
    entries.stereo_atoms.clear();
    let constitution = Molecule::from_entries(entries);

    assert!(constitution.canonical_eq_by(
        &stereo_atom_canonicalization_molecule,
        CanonicalizeLevel::Constitution,
        &canonicalize_context,
    ));
    assert!(!constitution.canonical_eq_by(
        &stereo_atom_canonicalization_molecule,
        CanonicalizeLevel::Structure,
        &canonicalize_context,
    ));
}

#[rstest]
fn test_molecule_canonical_hash_by_constitution(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let mut entries = molecule_entries(&stereo_atom_canonicalization_molecule);
    entries.stereo_atoms.clear();
    let constitution = Molecule::from_entries(entries);

    assert_eq!(
        constitution.canonical_hash_by(CanonicalizeLevel::Constitution, &canonicalize_context),
        stereo_atom_canonicalization_molecule
            .canonical_hash_by(CanonicalizeLevel::Constitution, &canonicalize_context),
    );
}

#[rstest]
fn test_molecule_canonical_eq_by_structure(canonicalize_context: CanonicalizeContext) {
    let plain = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let constrained = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4))
        ],
        ..Default::default()
    });

    assert!(plain.canonical_eq_by(
        &constrained,
        CanonicalizeLevel::Structure,
        &canonicalize_context,
    ));
    assert!(!plain.canonical_eq_by(&constrained, CanonicalizeLevel::Full, &canonicalize_context,));
    assert_eq!(
        plain.canonical_eq_by(&constrained, CanonicalizeLevel::Full, &canonicalize_context,),
        plain.canonical_eq(&constrained, &canonicalize_context),
    );
}

#[rstest]
fn test_molecule_canonical_hash_by_structure(canonicalize_context: CanonicalizeContext) {
    let plain = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let constrained = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4))
        ],
        ..Default::default()
    });

    assert_eq!(
        plain.canonical_hash_by(CanonicalizeLevel::Structure, &canonicalize_context),
        constrained.canonical_hash_by(CanonicalizeLevel::Structure, &canonicalize_context),
    );
}

fn reaction_canonicalization_fixture() -> Reaction {
    let atom = AtomForm::from_element(Element::C).with_charge(0_i64);
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom.clone(), atom],
            ..Default::default()
        }),
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    )
}

#[rstest]
#[case::topology(CanonicalizeLevel::Topology)]
#[case::constitution(CanonicalizeLevel::Constitution)]
#[case::structure(CanonicalizeLevel::Structure)]
#[case::full(CanonicalizeLevel::Full)]
fn test_reaction_canonicalize_by(
    canonicalize_context: CanonicalizeContext,
    #[case] level: CanonicalizeLevel,
) {
    let source = reaction_canonicalization_fixture();
    let expected = source
        .to_reaction_span()
        .expect("fixed reaction materializes")
        .canonicalize_by(level, &canonicalize_context)
        .expect("fixed span canonicalizes")
        .to_reaction();

    assert_eq!(
        source.clone().canonicalize_by(level, &canonicalize_context),
        Ok(expected.clone()),
    );
    if level == CanonicalizeLevel::Full {
        assert_eq!(source.canonicalize(&canonicalize_context), Ok(expected),);
    }
}

#[rstest]
#[case::integrity(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                ..Default::default()
            }),
            [Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(1),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(0),
                    new: NumForm::Lit(1),
                },
            })].into_iter().collect(),
        ),
        ReactionCanonicalizeError::Integrity(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
#[case::contradiction(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
                ..Default::default()
            }),
            [Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })].into_iter().collect(),
        ),
        ReactionCanonicalizeError::Contradiction(Contradiction),
    )]
fn test_reaction_canonicalize_error(
    canonicalize_context: CanonicalizeContext,
    #[case] reaction: Reaction,
    #[case] expected: ReactionCanonicalizeError,
) {
    assert_eq!(reaction.canonicalize(&canonicalize_context), Err(expected));
}

#[rstest]
fn test_reaction_canonical_eq(canonicalize_context: CanonicalizeContext) {
    let source = reaction_canonicalization_fixture();
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed reaction canonicalizes");

    assert!(source.canonical_eq(&canonical, &canonicalize_context));
    assert_eq!(
        source.canonical_eq_by(&canonical, CanonicalizeLevel::Full, &canonicalize_context,),
        source.canonical_eq(&canonical, &canonicalize_context),
    );
}

#[rstest]
#[case::topology(CanonicalizeLevel::Topology)]
#[case::constitution(CanonicalizeLevel::Constitution)]
#[case::structure(CanonicalizeLevel::Structure)]
fn test_reaction_canonical_eq_by(
    canonicalize_context: CanonicalizeContext,
    #[case] level: CanonicalizeLevel,
) {
    let lhs = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let identity = Reaction::new(lhs.clone(), Deltas::new());
    let excluded_contradiction = Reaction::new(
        lhs,
        [Delta::Constraint(ConstraintDelta::Remove(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]
        .into_iter()
        .collect(),
    );

    assert!(excluded_contradiction.canonical_eq_by(&identity, level, &canonicalize_context,));
    assert!(!excluded_contradiction.canonical_eq(&identity, &canonicalize_context));
}

#[rstest]
#[case::topology(CanonicalizeLevel::Topology)]
#[case::constitution(CanonicalizeLevel::Constitution)]
#[case::structure(CanonicalizeLevel::Structure)]
fn test_reaction_canonical_hash_by(
    canonicalize_context: CanonicalizeContext,
    #[case] level: CanonicalizeLevel,
) {
    let lhs = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let identity = Reaction::new(lhs.clone(), Deltas::new());
    let excluded_contradiction = Reaction::new(
        lhs,
        [Delta::Constraint(ConstraintDelta::Remove(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]
        .into_iter()
        .collect(),
    );

    assert_eq!(
        excluded_contradiction.canonical_hash_by(level, &canonicalize_context),
        identity.canonical_hash_by(level, &canonicalize_context),
    );
}

#[rstest]
fn test_reaction_canonical_eq_error(canonicalize_context: CanonicalizeContext) {
    let lhs = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
        ..Default::default()
    });
    let contradiction = |new| {
        Reaction::new(
            lhs.clone(),
            [Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(new),
                },
            })]
            .into_iter()
            .collect(),
        )
    };
    let left_contradiction = contradiction(2);
    let right_contradiction = contradiction(3);
    let malformed = Reaction::new(
        lhs,
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    );

    assert!(left_contradiction.canonical_eq(&right_contradiction, &canonicalize_context));
    assert!(malformed.canonical_eq(&malformed, &canonicalize_context));
    assert!(!malformed.canonical_eq_by(
        &left_contradiction,
        CanonicalizeLevel::Topology,
        &canonicalize_context,
    ));
}

#[rstest]
fn test_molecule_canonical_eq_by_contradiction(canonicalize_context: CanonicalizeContext) {
    let mut left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    left.atom_mut(AtomId(0)).attributes.constraints =
        AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())).into();
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });

    assert!(left.canonical_eq_by(&right, CanonicalizeLevel::Structure, &canonicalize_context,));
    assert!(!left.canonical_eq(&right, &canonicalize_context));
}

#[rstest]
fn test_constraint_key() {
    let atom = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4));
    let bond = Constraint::Bond(
        BondId(0),
        BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
    );
    let left = Constraint::And(vec![atom.clone(), bond.clone(), atom.clone()]);
    let right = Constraint::And(vec![bond, atom]);

    assert_eq!(constraint_key(&left), constraint_key(&right));
}

#[rstest]
fn test_canonicalize_structure_selected_layer(
    initial_class_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let canonical = canonicalize_structure(&initial_class_molecule, &canonicalize_context)
        .expect("fixed molecule canonicalizes");
    let canonical_again = canonicalize_structure(&canonical, &canonicalize_context)
        .expect("canonical molecule canonicalizes");

    assert_eq!(
        selected_structure_key(&canonical_again),
        selected_structure_key(&canonical),
    );
    assert_eq!(canonical_again.check_integrity(), Ok(()));
}

#[rstest]
fn test_canonicalize_structure_renumbering(
    symmetric_stereo_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let canonical = canonicalize_structure(
        &symmetric_stereo_canonicalization_molecule,
        &canonicalize_context,
    )
    .expect("fixed molecule canonicalizes");
    let expected = selected_structure_key(&canonical);

    for rank in 0..(1..=5).product() {
        let permutation = Permutation::unrank(5, rank);
        let correspondence = molecule_correspondence(&[
            (0..5).map(|index| permutation.apply(index)).collect(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![0],
            Vec::new(),
        ]);
        let renumbered = symmetric_stereo_canonicalization_molecule.remap(&correspondence);
        let actual = canonicalize_structure(&renumbered, &canonicalize_context)
            .expect("renumbered molecule canonicalizes");

        assert_eq!(selected_structure_key(&actual), expected, "rank {rank}");
    }
}

#[rstest]
#[case::nauty(AutomorphismAlgorithm::Nauty)]
fn test_canonicalize_structure_minimum(
    symmetric_stereo_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] algorithm: AutomorphismAlgorithm,
) {
    let incidence_graph =
        symmetric_stereo_canonicalization_molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) = initial_class_keys(
        &symmetric_stereo_canonicalization_molecule,
        &incidence_graph,
    )
    .expect("fixed molecule has initial classes");
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let constitution_classes =
        constitution_entity_classes(&symmetric_stereo_canonicalization_molecule)
            .expect("fixed molecule has constitution classes");
    let descriptors = structure_partition_descriptors(
        &symmetric_stereo_canonicalization_molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        &constitution_classes,
    )
    .expect("fixed molecule has structure descriptors");
    let leaf_candidate = |order: &[NodeId]| {
        structure_candidate(
            &symmetric_stereo_canonicalization_molecule,
            &incidence_graph,
            order,
        )
        .expect("structure descriptors establish normalization")
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let expected = exhaustive_minimum(
        &adapter,
        adapter_entity_blocks(&incidence_graph),
        &leaf_candidate,
    );

    for options in [
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: node_branch_order,
        },
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: reverse_node_branch_order,
        },
    ] {
        let actual = canonical_search(
            &adapter,
            &descriptors,
            algorithm,
            options,
            &leaf_candidate,
            &no_prefix,
        );
        assert_eq!(actual.candidate.key, expected.key, "{options:?}");
    }

    for options in [
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: reverse_node_branch_order,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    ] {
        let (canonical, _) = canonicalize_structure_with_options(
            &symmetric_stereo_canonicalization_molecule,
            &CanonicalizeContext {
                automorphism_algorithm: algorithm,
                ..canonicalize_context
            },
            options,
        )
        .expect("fixed molecule canonicalizes");
        assert_eq!(
            selected_structure_key(&canonical),
            expected.key,
            "{options:?}"
        );
    }
}

#[rstest]
fn test_canonicalize_structure_meso(
    meso_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let correspondence = molecule_correspondence(&[
        vec![1, 0, 3, 2, 5, 4],
        vec![0, 3, 4, 1, 2],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![1, 0],
        Vec::new(),
    ]);
    let renumbered = meso_canonicalization_molecule.remap(&correspondence);

    assert_eq!(
        canonicalize_structure(&renumbered, &canonicalize_context),
        canonicalize_structure(&meso_canonicalization_molecule, &canonicalize_context,),
    );
}

#[rstest]
#[case::kinded(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))]
#[case::undetermined(StereoConfigurationForm::Undetermined)]
fn test_canonicalize_structure_repeated_ligands(
    repeated_ligand_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] configuration: StereoConfigurationForm,
) {
    let mut entries = molecule_entries(&repeated_ligand_canonicalization_molecule);
    entries.stereo_atoms[0].2.configuration = configuration;
    let source = Molecule::from_entries(entries);
    let reframed = reframe_stereo_atom(
        &source,
        StereoAtomId(0),
        Permutation::from_image(&[0, 1, 3, 2]),
    );

    assert_eq!(
        canonicalize_structure(&reframed, &canonicalize_context),
        canonicalize_structure(&source, &canonicalize_context),
    );
}

fn rank_paired_initial_classes(
    left: (&[InitialClassKey], &[InitialClassKey]),
    right: (&[InitialClassKey], &[InitialClassKey]),
) -> (InitialClasses, InitialClasses) {
    let ranks = left
        .0
        .iter()
        .chain(left.1)
        .chain(right.0)
        .chain(right.1)
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .enumerate()
        .map(|(rank, key)| (key, rank as u32))
        .collect::<BTreeMap<_, _>>();
    let rank =
        |entity_keys: &[InitialClassKey], incidence_keys: &[InitialClassKey]| InitialClasses {
            entities: entity_keys.iter().map(|key| ranks[key]).collect(),
            incidences: incidence_keys.iter().map(|key| ranks[key]).collect(),
        };

    (rank(left.0, left.1), rank(right.0, right.1))
}

fn colored_encoding_equivalent(left: &Molecule, right: &Molecule, level: IncidenceLevel) -> bool {
    fn key(adapter: &AutomorphismAdapter) -> Vec<u8> {
        adapter.graph().canonical_key(
            |node| {
                let (domain, rank) = match adapter.class(node) {
                    AutomorphismClass::Entity(rank) => (0, rank),
                    AutomorphismClass::Incidence(rank) => (1, rank),
                };
                let mut color = vec![domain];
                color.extend_from_slice(&rank.to_be_bytes());
                color
            },
            |_| Vec::new(),
            AutomorphismAlgorithm::Nauty,
        )
    }

    let left_incidence = left.incidence_graph(level);
    let right_incidence = right.incidence_graph(level);
    let (left_entity_keys, left_incidence_keys) =
        initial_class_keys(left, &left_incidence).unwrap();
    let (right_entity_keys, right_incidence_keys) =
        initial_class_keys(right, &right_incidence).unwrap();
    let (left_classes, right_classes) = rank_paired_initial_classes(
        (&left_entity_keys, &left_incidence_keys),
        (&right_entity_keys, &right_incidence_keys),
    );
    let left_adapter = AutomorphismAdapter::new(&left_incidence, &left_classes);
    let right_adapter = AutomorphismAdapter::new(&right_incidence, &right_classes);

    key(&left_adapter) == key(&right_adapter)
}

fn permutations(count: usize) -> Vec<Vec<usize>> {
    fn visit(values: &mut [usize], position: usize, output: &mut Vec<Vec<usize>>) {
        if position == values.len() {
            output.push(values.to_vec());
            return;
        }
        for next in position..values.len() {
            values.swap(position, next);
            visit(values, position + 1, output);
            values.swap(position, next);
        }
    }

    let mut values = (0..count).collect::<Vec<_>>();
    let mut output = Vec::new();
    visit(&mut values, 0, &mut output);
    output
}

fn explicitly_dense_equivalent(left: &Molecule, right: &Molecule) -> bool {
    fn visit(
        family: usize,
        permutations: &[Vec<Vec<usize>>; 8],
        images: &mut [Vec<usize>; 8],
        left: &Molecule,
        right: &Molecule,
    ) -> bool {
        if family == images.len() {
            return left.equiv_under(right, &molecule_correspondence(images));
        }

        permutations[family].iter().any(|permutation| {
            images[family].clone_from(permutation);
            visit(family + 1, permutations, images, left, right)
        })
    }

    let left_counts = molecule_counts(left);
    if left_counts != molecule_counts(right) {
        return false;
    }
    let permutations = left_counts.map(permutations);
    visit(
        0,
        &permutations,
        &mut array::from_fn(|_| Vec::new()),
        left,
        right,
    )
}

fn reverse_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
    let images = molecule_counts(molecule).map(|count| (0..count).rev().collect::<Vec<_>>());
    molecule_correspondence(&images)
}

fn direct_graph_adapter(source: &Graph) -> AutomorphismAdapter {
    AutomorphismAdapter {
        graph: source.clone(),
        classes: vec![AutomorphismClass::Entity(0); source.node_count()],
        node_sources: source.node_ids().map(SubdivisionNodeSource::Node).collect(),
        incidence_nodes: vec![None; source.edge_count()],
        source_node_count: source.node_count(),
    }
}

fn project_entries(mut entries: MoleculeEntries, level: IncidenceLevel) -> MoleculeEntries {
    entries.constraints = Constraints::new();
    match level {
        IncidenceLevel::Topology => {
            entries.dative.clear();
            entries.aromatic.clear();
            entries.multicenter.clear();
            entries.noncovalent.clear();
            entries.stereo_atoms.clear();
            entries.stereo_bonds.clear();
        }
        IncidenceLevel::Constitution => {
            entries.stereo_atoms.clear();
            entries.stereo_bonds.clear();
        }
        IncidenceLevel::Full => {}
    }
    entries
}

fn encoding_entries() -> MoleculeEntries {
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![(
            vec![AtomId(0), AtomId(1)],
            AtomId(2),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![1, 2, 3]),
        )],
        multicenter: vec![(
            vec![AtomId(1), AtomId(2), AtomId(3)],
            MulticenterBondForm::from_electrons(vec![2, 0, 1]),
        )],
        noncovalent: vec![(
            AtomId(2),
            AtomId(3),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    }
}

#[rstest]
fn test_incidence_cmp() {
    let incidences = [
        Incidence::BondEndpoint,
        Incidence::DativeDonor,
        Incidence::DativeAcceptor,
        Incidence::AromaticParticipant(NumForm::Undetermined),
        Incidence::AromaticParticipant(NumForm::Lit(-1)),
        Incidence::AromaticParticipant(NumForm::Lit(1)),
        Incidence::AromaticParticipant(NumForm::lit_set([0, 1])),
        Incidence::AromaticParticipant(NumForm::RangeFrom(0)),
        Incidence::AromaticParticipant(NumForm::RangeTo(0)),
        Incidence::AromaticParticipant(NumForm::var("x")),
        Incidence::AromaticParticipant(NumForm::pred_expr(PredExpr::Rel(
            ArithExpr::Var("x".into()),
            RelOp::Eq,
            ArithExpr::Lit(0),
        ))),
        Incidence::MulticenterParticipant(NumForm::Undetermined),
        Incidence::NoncovalentEndpoint,
        Incidence::StereoSite,
        Incidence::StereoLigand(StereoLigandKind::Atom),
        Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
        Incidence::StereoLigand(StereoLigandKind::LonePair),
    ];

    for pair in incidences.windows(2) {
        assert_eq!(pair[0].cmp(&pair[1]), Ordering::Less);
    }
    for lhs in &incidences {
        for rhs in &incidences {
            assert_eq!(
                lhs.cmp(rhs),
                incidence_key(lhs)
                    .unwrap()
                    .cmp(&incidence_key(rhs).unwrap()),
            );
        }
    }
}

#[rstest]
#[case::atom(
        Entity::Atom(AtomId(0)),
        vec![
            FieldPosition(0),
            FieldPosition(1),
            FieldPosition(2),
            FieldPosition(3),
            FieldPosition(4),
            FieldPosition(5),
        ]
    )]
#[case::bond(
        Entity::Bond(BondId(0)),
        vec![FieldPosition(1), FieldPosition(2), FieldPosition(3)]
    )]
#[case::dative_bond(
        Entity::DativeBond(DativeBondId(0)),
        vec![FieldPosition(2)]
    )]
#[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(0)),
        vec![FieldPosition(2), FieldPosition(3)]
    )]
#[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(0)),
        vec![FieldPosition(2), FieldPosition(3)]
    )]
#[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        vec![FieldPosition(1)]
    )]
#[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(0)),
        vec![FieldPosition(2)]
    )]
#[case::stereo_bond(
        Entity::StereoBond(StereoBondId(0)),
        vec![FieldPosition(2)]
    )]
fn test_entity_class_key_field_positions(
    initial_class_molecule: Molecule,
    #[case] entity: Entity,
    #[case] expected: Vec<FieldPosition>,
) {
    let InitialClassKey::Entity {
        value: CanonicalKeyValue::Product(fields),
        ..
    } = entity_class_key(&initial_class_molecule, entity).unwrap()
    else {
        panic!("entity class key must be a product");
    };

    assert_eq!(
        fields
            .into_iter()
            .map(|field| field.position)
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::normalized_atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), true)]
#[case::atom_element(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
#[case::atom_charge(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), false)]
#[case::normalized_bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), true)]
#[case::bond_order(Entity::Bond(BondId(0)), Entity::Bond(BondId(2)), false)]
#[case::normalized_dative(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(1)),
    true
)]
#[case::dative_order(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(2)),
    false
)]
#[case::aromatic_electrons_excluded(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    true
)]
#[case::aromatic_charge(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(2)),
    false
)]
#[case::multicenter_electrons_excluded(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(1)),
    true
)]
#[case::noncovalent_kind_equal(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(1)),
    true
)]
#[case::noncovalent_kind_distinct(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(2)),
    false
)]
#[case::stereo_atom_configuration_excluded(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    true
)]
#[case::stereo_atom_kind(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(2)),
    false
)]
#[case::stereo_bond_configuration_excluded(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    true
)]
#[case::entity_kind(Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
fn test_initial_classes(
    initial_class_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected_equal: bool,
) {
    let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Full);
    let classes = initial_classes(&initial_class_molecule, &incidence_graph).unwrap();
    let lhs_class = classes.entities[incidence_graph.node_of(lhs).index()];
    let rhs_class = classes.entities[incidence_graph.node_of(rhs).index()];

    assert_eq!(lhs_class == rhs_class, expected_equal);
}

#[rstest]
fn test_initial_classes_incidence(initial_class_molecule: Molecule) {
    let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Full);
    let classes = initial_classes(&initial_class_molecule, &incidence_graph).unwrap();
    let incidences = incidence_graph
        .incidences()
        .map(|(edge, incidence)| (incidence, classes.incidences[edge.index()]))
        .collect::<Vec<_>>();

    for (lhs, lhs_class) in &incidences {
        for (rhs, rhs_class) in &incidences {
            assert_eq!(lhs_class == rhs_class, lhs == rhs);
        }
    }
    for entity_class in &classes.entities {
        for (_, incidence_class) in &incidences {
            assert_ne!(entity_class, incidence_class);
        }
    }
}

#[rstest]
fn test_topology_comparison_key() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C)
                .with_charge(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(1),
                    ArithExpr::Lit(2),
                ]))))
                .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([]))),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![(
            AtomId(0),
            AtomId(1),
            BondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                ArithExpr::Lit(0),
                ArithExpr::Lit(1),
            ]))))
            .with_charge(-1_i64),
        )],
        dative: vec![(
            vec![AtomId(0)],
            AtomId(1),
            DativeBondForm::new(NumForm::lit_set([])),
        )],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let order = [
        incidence_graph.node_of(Entity::Atom(AtomId(0))),
        incidence_graph.node_of(Entity::Atom(AtomId(1))),
        incidence_graph.node_of(Entity::Bond(BondId(0))),
    ];
    let undetermined_num = NumForm::Undetermined;
    let undetermined_isotope = IsotopeMassForm::Undetermined;
    let undetermined_spin = UnpairedElectronsForm::default();

    assert_eq!(
        topology_comparison_key(&molecule, &incidence_graph, &order),
        Ok(CanonicalComparisonKey {
            entity_blocks: vec![
                PositionedKey {
                    position: EntityBlockPosition::ATOM,
                    value: sequence([
                        product([
                            element_form_key(&ElementForm::Lit(Element::C)),
                            isotope_mass_form_key(&undetermined_isotope),
                            num_form_key(&NumForm::Lit(3)),
                            num_form_key(&undetermined_num),
                            num_form_key(&undetermined_num),
                            unpaired_electrons_form_key(&undetermined_spin),
                        ]),
                        product([
                            element_form_key(&ElementForm::Lit(Element::O)),
                            isotope_mass_form_key(&undetermined_isotope),
                            num_form_key(&undetermined_num),
                            num_form_key(&undetermined_num),
                            num_form_key(&undetermined_num),
                            unpaired_electrons_form_key(&undetermined_spin),
                        ]),
                    ]),
                },
                PositionedKey {
                    position: EntityBlockPosition::BOND,
                    value: sequence([positioned_product([
                        (
                            0,
                            product([
                                CanonicalKeyValue::Unsigned(0),
                                CanonicalKeyValue::Unsigned(1),
                            ]),
                        ),
                        (1, num_form_key(&NumForm::Lit(1))),
                        (2, num_form_key(&NumForm::Lit(-1))),
                        (3, unpaired_electrons_form_key(&undetermined_spin)),
                    ])]),
                },
            ],
            constraints: Vec::new(),
        }),
    );
}

#[rstest]
fn test_topology_comparison_key_dense_remapping() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (
                AtomId(1),
                AtomId(2),
                BondForm::from_order(1).with_charge(-1_i64),
            ),
        ],
        ..Default::default()
    });
    let remapped = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (
                AtomId(0),
                AtomId(2),
                BondForm::from_order(1).with_charge(-1_i64),
            ),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
        ],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let remapped_incidence_graph = remapped.incidence_graph(IncidenceLevel::Topology);
    let order = [
        incidence_graph.node_of(Entity::Atom(AtomId(2))),
        incidence_graph.node_of(Entity::Atom(AtomId(0))),
        incidence_graph.node_of(Entity::Atom(AtomId(1))),
        incidence_graph.node_of(Entity::Bond(BondId(1))),
        incidence_graph.node_of(Entity::Bond(BondId(0))),
    ];
    let remapped_order = remapped_incidence_graph
        .graph()
        .node_ids()
        .collect::<Vec<_>>();

    assert_eq!(
        topology_comparison_key(&molecule, &incidence_graph, &order),
        topology_comparison_key(&remapped, &remapped_incidence_graph, &remapped_order),
    );
}

#[rstest]
#[case::dative_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            dative: vec![(
                vec![AtomId(0), AtomId(1)],
                AtomId(2),
                DativeBondForm::new(NumForm::RangeFrom(-1)),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::DATIVE_BOND,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(1),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (1, CanonicalKeyValue::Unsigned(0)),
            (2, num_form_key(&NumForm::RangeFrom(-1))),
        ]),
    )]
#[case::aromatic_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2])
                    .with_charge(NumForm::var("q")),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::AROMATIC_SYSTEM,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                electron_counts_form_key(&ElectronCountsForm::Lit(vec![2, 1])),
            ),
            (2, num_form_key(&NumForm::var("q"))),
            (
                3,
                unpaired_electrons_form_key(&UnpairedElectronsForm::default()),
            ),
        ]),
    )]
#[case::multicenter_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            multicenter: vec![(
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![2, 0])
                    .with_charge(NumForm::RangeTo(1)),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::MULTICENTER_BOND,
        positioned_product([
            (
                0,
                sequence([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                electron_counts_form_key(&ElectronCountsForm::Lit(vec![0, 2])),
            ),
            (2, num_form_key(&NumForm::RangeTo(1))),
            (
                3,
                unpaired_electrons_form_key(&UnpairedElectronsForm::default()),
            ),
        ]),
    )]
#[case::noncovalent_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            noncovalent: vec![(
                AtomId(0),
                AtomId(2),
                NoncovalentBondForm::default(),
            )],
            ..Default::default()
        }),
        EntityBlockPosition::NONCOVALENT_BOND,
        positioned_product([
            (
                0,
                product([
                    CanonicalKeyValue::Unsigned(0),
                    CanonicalKeyValue::Unsigned(2),
                ]),
            ),
            (
                1,
                noncovalent_bond_kind_form_key(&NoncovalentBondKindForm::Undetermined),
            ),
        ]),
    )]
fn test_constitution_comparison_key(
    #[case] molecule: Molecule,
    #[case] position: EntityBlockPosition,
    #[case] expected: CanonicalKeyValue,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let mut atom_ids = molecule.atoms().ids().collect::<Vec<_>>();
    atom_ids.reverse();
    let mut order = atom_ids
        .into_iter()
        .map(|id| incidence_graph.node_of(Entity::Atom(id)))
        .collect::<Vec<_>>();
    order.extend(
        incidence_graph
            .graph()
            .node_ids()
            .filter(|&node| !matches!(incidence_graph.entity(node), Entity::Atom(_))),
    );
    let key = constitution_comparison_key(&molecule, &incidence_graph, &order).unwrap();

    assert_eq!(
        key.entity_blocks
            .into_iter()
            .find(|block| block.position == position),
        Some(PositionedKey {
            position,
            value: sequence([expected]),
        }),
    );
}

#[rstest]
fn test_constitution_comparison_key_dense_remapping() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::F),
        ],
        dative: vec![(
            vec![AtomId(0), AtomId(2)],
            AtomId(1),
            DativeBondForm::new(NumForm::RangeFrom(1)),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]).with_charge(NumForm::var("q")),
        )],
        multicenter: vec![(
            vec![AtomId(1), AtomId(3)],
            MulticenterBondForm::new(ElectronCountsForm::Undetermined)
                .with_charge(NumForm::RangeTo(2)),
        )],
        noncovalent: vec![(AtomId(0), AtomId(3), NoncovalentBondForm::default())],
        ..Default::default()
    });
    let correspondence = molecule_correspondence(&[
        vec![3, 1, 0, 2],
        Vec::new(),
        vec![0],
        vec![0],
        vec![0],
        vec![0],
        Vec::new(),
        Vec::new(),
    ]);
    let remapped = molecule.remap(&correspondence);
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let remapped_incidence_graph = remapped.incidence_graph(IncidenceLevel::Constitution);
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    order.reverse();
    let remapped_order = order
        .iter()
        .map(|&node| incidence_graph.entity(node))
        .map(|entity| {
            correspondence
                .right_of(entity)
                .expect("dense correspondence maps every entity")
        })
        .map(|entity| remapped_incidence_graph.node_of(entity))
        .collect::<Vec<_>>();

    assert_eq!(
        constitution_comparison_key(&molecule, &incidence_graph, &order),
        constitution_comparison_key(&remapped, &remapped_incidence_graph, &remapped_order),
    );
}

#[rstest]
fn test_constitution_comparison_key_excluded_data() {
    let molecule = Molecule::from_entries(project_entries(
        encoding_entries(),
        IncidenceLevel::Constitution,
    ));
    let mut excluded = Molecule::from_entries(encoding_entries());
    excluded
        .modify_atoms(|atom| atom.with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))));
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let excluded_incidence_graph = excluded.incidence_graph(IncidenceLevel::Constitution);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let excluded_order = excluded_incidence_graph
        .graph()
        .node_ids()
        .collect::<Vec<_>>();

    assert_eq!(
        constitution_comparison_key(&molecule, &incidence_graph, &order),
        constitution_comparison_key(&excluded, &excluded_incidence_graph, &excluded_order),
    );
}

#[rstest]
#[case::dative(Entity::DativeBond(DativeBondId(0)))]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
fn test_correspondence_from_order(initial_class_molecule: Molecule, #[case] excluded: Entity) {
    let incidence_graph = initial_class_molecule.incidence_graph(IncidenceLevel::Topology);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let correspondence =
        correspondence_from_order(&initial_class_molecule, &incidence_graph, &order);

    assert_eq!(correspondence.right_of(excluded), Some(excluded));
}

#[rstest]
#[case::localized_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        vec![Incidence::BondEndpoint, Incidence::BondEndpoint],
        0,
    )]
#[case::repeated_virtual_ligand_anchor(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            Incidence::StereoSite,
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
        ],
        5,
    )]
fn test_automorphism_adapter_new(
    #[case] molecule: Molecule,
    #[case] expected_incidences: Vec<Incidence>,
    #[case] expected_occurrence_nodes: usize,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let classes = initial_classes(&molecule, &incidence_graph).unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let source = incidence_graph.graph();

    assert_eq!(
        incidence_graph
            .incidences()
            .map(|(_, incidence)| incidence.clone())
            .collect::<Vec<_>>(),
        expected_incidences,
    );
    assert_eq!(
        adapter.graph().node_count(),
        source.node_count() + expected_occurrence_nodes,
    );
    assert_eq!(
        adapter.graph().edge_count(),
        source.edge_count() + expected_occurrence_nodes,
    );
    assert!(adapter.graph().is_simple());

    for node in source.node_ids() {
        let adapter_node = adapter
            .node_of(SubdivisionNodeSource::Node(node))
            .expect("every source entity remains an adapter node");
        assert_eq!(
            adapter.node_source(adapter_node),
            SubdivisionNodeSource::Node(node),
        );
        assert_eq!(
            adapter.class(adapter_node),
            AutomorphismClass::Entity(classes.entities[node.index()]),
        );
    }
    for edge in source.edge_ids() {
        if let Some(adapter_node) = adapter.node_of(SubdivisionNodeSource::Edge(edge)) {
            assert_eq!(
                adapter.node_source(adapter_node),
                SubdivisionNodeSource::Edge(edge),
            );
            assert_eq!(
                adapter.class(adapter_node),
                AutomorphismClass::Incidence(classes.incidences[edge.index()]),
            );
            assert_eq!(adapter.graph().degree(adapter_node), 2);
        }
    }
}

#[rstest]
fn test_automorphism_adapter_automorphisms() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let classes = initial_classes(&molecule, &incidence_graph).unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);

    assert_eq!(
        adapter.automorphisms(AutomorphismAlgorithm::Nauty),
        ProjectedAutomorphismOutput {
            orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
            canonical_labels: vec![NodeId(0), NodeId(1), NodeId(2)],
            generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
        },
    );
}

#[rstest]
#[case::ordered_classes(
        vec![2_u32, 0, 2, 1, 0],
        OrderedPartition {
            cells: vec![
                vec![NodeId(1), NodeId(4)],
                vec![NodeId(3)],
                vec![NodeId(0), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_from_descriptors(
    #[case] descriptors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(OrderedPartition::from_descriptors(&descriptors), expected,);
}

#[rstest]
#[case::path(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        OrderedPartition {
            cells: vec![vec![NodeId(1), NodeId(2)], vec![NodeId(0), NodeId(3)]],
        },
    )]
fn test_ordered_partition_refine(#[case] graph: Graph, #[case] expected: OrderedPartition) {
    assert_eq!(
        OrderedPartition::from_descriptors(&[0_u32; 4]).refine(&graph),
        expected,
    );
}

#[rstest]
#[case::first_cell(
        OrderedPartition {
            cells: vec![vec![NodeId(0), NodeId(3)], vec![NodeId(1), NodeId(2)]],
        },
        0,
        NodeId(3),
        OrderedPartition {
            cells: vec![
                vec![NodeId(3)],
                vec![NodeId(0)],
                vec![NodeId(1), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_individualize(
    #[case] partition: OrderedPartition,
    #[case] cell_index: usize,
    #[case] node: NodeId,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(partition.individualize(cell_index, node), expected);
}

#[rstest]
#[case::path(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
#[case::branched(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
#[case::distinct_attributes(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
    )]
fn test_canonical_search(#[case] molecule: Molecule) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph).unwrap();
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        topology_candidate(&molecule, &incidence_graph, order)
            .expect("selected topology values normalize")
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let expected = exhaustive_minimum(
        &adapter,
        adapter_entity_blocks(&incidence_graph),
        &leaf_candidate,
    );
    let unpruned = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: node_branch_order,
        },
        &leaf_candidate,
        &no_prefix,
    );
    let reversed = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: reverse_node_branch_order,
        },
        &leaf_candidate,
        &no_prefix,
    );
    let pruned = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
        &leaf_candidate,
        &no_prefix,
    );

    assert_eq!(unpruned.candidate.key, expected.key);
    assert_eq!(reversed.candidate.key, expected.key);
    assert_eq!(pruned.candidate.key, expected.key);
    assert!(pruned.stats.visited_leaves <= unpruned.stats.visited_leaves);
}

#[rstest]
#[case::singleton(
    1,
    CanonicalSearchOptions {
        automorphism_pruning: false,
        prefix_pruning: false,
        branch_order: node_branch_order,
    },
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![],
        refinement_calls: 1,
        branch_order_calls: 0,
        backend_calls: 0,
        visited_leaves: 1,
        leaf_comparisons: 0,
        prefix_pruned_branches: 0,
        orbit_pruned_branches: 0,
    },
)]
#[case::symmetric(
    2,
    CanonicalSearchOptions {
        automorphism_pruning: false,
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    },
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![2],
        refinement_calls: 3,
        branch_order_calls: 1,
        backend_calls: 1,
        visited_leaves: 2,
        leaf_comparisons: 1,
        prefix_pruned_branches: 0,
        orbit_pruned_branches: 0,
    },
)]
#[case::orbit_pruned(
    2,
    CanonicalSearchOptions {
        automorphism_pruning: true,
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    },
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![2],
        refinement_calls: 2,
        branch_order_calls: 1,
        backend_calls: 1,
        visited_leaves: 1,
        leaf_comparisons: 0,
        prefix_pruned_branches: 0,
        orbit_pruned_branches: 1,
    },
)]
fn test_canonical_search_stats(
    #[case] node_count: usize,
    #[case] options: CanonicalSearchOptions,
    #[case] expected: CanonicalSearchStats,
) {
    let source = Graph::new(node_count, &[]);
    let adapter = direct_graph_adapter(&source);
    let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
        key: order.to_vec(),
        entity_order: order.to_vec(),
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<Vec<NodeId>>| false;
    let actual = canonical_search(
        &adapter,
        &adapter.classes,
        AutomorphismAlgorithm::Nauty,
        options,
        &leaf_candidate,
        &no_prefix,
    );

    assert_eq!(actual.stats, expected);
}

#[rstest]
fn test_canonical_search_color_classes() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, incidence_keys) = initial_class_keys(&molecule, &incidence_graph).unwrap();
    let classes = rank_initial_classes(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
    let mut relabeled = adapter.clone();
    relabeled.classes.iter_mut().for_each(|class| {
        *class = match *class {
            AutomorphismClass::Entity(value) => AutomorphismClass::Entity(u32::MAX - value),
            AutomorphismClass::Incidence(value) => AutomorphismClass::Incidence(u32::MAX - value),
        }
    });
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        topology_candidate(&molecule, &incidence_graph, order)
            .expect("selected topology values normalize")
    };
    let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
    let options = CanonicalSearchOptions {
        automorphism_pruning: true,
        prefix_pruning: false,
        branch_order: backend_canonical_branch_order,
    };

    let expected = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        options,
        &leaf_candidate,
        &no_prefix,
    );
    let actual = canonical_search(
        &relabeled,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        options,
        &leaf_candidate,
        &no_prefix,
    );

    assert_eq!(actual.candidate.key, expected.candidate.key);
}

#[rstest]
fn test_canonicalize_topology(canonicalize_context: CanonicalizeContext) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::C).with_charge(NumForm::ArithExpr(Box::new(
                ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(2)]),
            ))),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (
                AtomId(1),
                AtomId(2),
                BondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(0),
                    ArithExpr::Lit(1),
                ])))),
            ),
        ],
        dative: vec![(
            vec![AtomId(0)],
            AtomId(2),
            DativeBondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                ArithExpr::Lit(0),
                ArithExpr::Lit(1),
            ])))),
        )],
        ..Default::default()
    });
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(2)),
        ],
        dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    });

    assert_eq!(
        canonicalize_topology(&molecule, &canonicalize_context),
        Ok(expected),
    );
}

#[rstest]
#[case::disconnected(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::H),
            ],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(2)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        }),
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(2)),
            ],
            ..Default::default()
        }),
    )]
fn test_canonicalize_topology_components(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
    #[case] expected: Molecule,
) {
    assert_eq!(
        canonicalize_topology(&molecule, &canonicalize_context),
        Ok(expected),
    );
}

#[rstest]
#[case::selected_atom(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)
                .with_charge(NumForm::LitSet(Box::default()))],
            ..Default::default()
        }),
    )]
#[case::excluded_dative(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::B),
            ],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                DativeBondForm::new(NumForm::LitSet(Box::default())),
            )],
            ..Default::default()
        }),
    )]
fn test_canonicalize_topology_error(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
) {
    assert_eq!(
        canonicalize_topology(&molecule, &canonicalize_context),
        Err(MoleculeCanonicalizeError::Contradiction(Contradiction,)),
    );
}

#[rstest]
fn test_canonicalize_topology_excluded_data(canonicalize_context: CanonicalizeContext) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    });
    let remapping = molecule_correspondence(&[
        vec![1, 0],
        Vec::new(),
        vec![0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ]);
    let remapped = molecule.remap(&remapping);

    let (canonical, correspondence) = canonicalize_topology_with_options(
        &molecule,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let (canonical_remapped, remapped_correspondence) = canonicalize_topology_with_options(
        &remapped,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Topology);
    let canonical_remapped_incidence = canonical_remapped.incidence_graph(IncidenceLevel::Topology);
    let canonical_again = canonicalize_topology(&canonical, &canonicalize_context).unwrap();
    let canonical_again_incidence = canonical_again.incidence_graph(IncidenceLevel::Topology);

    assert!(molecule.equiv_under(&canonical, &correspondence));
    assert!(remapped.equiv_under(&canonical_remapped, &remapped_correspondence));
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical_remapped.check_integrity(), Ok(()));
    assert_eq!(
        topology_comparison_key(
            &canonical,
            &canonical_incidence,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
        topology_comparison_key(
            &canonical_again,
            &canonical_again_incidence,
            &canonical_again_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>(),
        ),
    );
    assert_eq!(
        topology_comparison_key(
            &canonical,
            &canonical_incidence,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
        topology_comparison_key(
            &canonical_remapped,
            &canonical_remapped_incidence,
            &canonical_remapped_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>(),
        ),
    );
}

#[rstest]
#[case::order_four(4)]
fn test_canonicalize_topology_exhaustive_domain(
    canonicalize_context: CanonicalizeContext,
    #[case] atom_count: usize,
) {
    let endpoint_pairs = (0..atom_count as u32)
        .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
        let bonds = endpoint_pairs
            .iter()
            .enumerate()
            .filter_map(|(position, &[first, second])| {
                ((edge_mask >> position) & 1 == 1).then_some((
                    AtomId(first),
                    AtomId(second),
                    BondForm::from_order(1),
                ))
            })
            .collect::<Vec<_>>();
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); atom_count],
            bonds,
            ..Default::default()
        });
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
        let (entity_keys, incidence_keys) =
            initial_class_keys(&molecule, &incidence_graph).unwrap();
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let leaf_candidate =
            |order: &[NodeId]| topology_candidate(&molecule, &incidence_graph, order).unwrap();
        let expected = exhaustive_minimum(
            &adapter,
            adapter_entity_blocks(&incidence_graph),
            &leaf_candidate,
        );

        let (canonical, correspondence) = canonicalize_topology_with_options(
            &molecule,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: backend_canonical_branch_order,
            },
        )
        .unwrap();
        let (unpruned, _) = canonicalize_topology_with_options(
            &molecule,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: reverse_node_branch_order,
            },
        )
        .unwrap();
        let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Topology);
        let canonical_order = canonical_incidence.graph().node_ids().collect::<Vec<_>>();

        assert_eq!(
            topology_comparison_key(&canonical, &canonical_incidence, &canonical_order),
            Ok(expected.key),
            "edge mask {edge_mask:#08b}",
        );
        assert_eq!(unpruned, canonical, "edge mask {edge_mask:#08b}");
        assert!(
            molecule.equiv_under(&canonical, &correspondence),
            "edge mask {edge_mask:#08b}",
        );
        assert_eq!(canonical.check_integrity(), Ok(()));
        assert_eq!(
            canonicalize_topology(&canonical, &canonicalize_context),
            Ok(canonical.clone()),
            "edge mask {edge_mask:#08b}",
        );

        for (index, atom_images) in permutations(atom_count).into_iter().enumerate() {
            let bond_count = molecule.bonds().count();
            let bond_images = if index % 2 == 0 {
                (0..bond_count).collect::<Vec<_>>()
            } else {
                (0..bond_count).rev().collect::<Vec<_>>()
            };
            let renumbering = molecule_correspondence(&[
                atom_images,
                bond_images,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ]);
            let renumbered = molecule.remap(&renumbering);

            assert_eq!(
                canonicalize_topology(&renumbered, &canonicalize_context),
                Ok(canonical.clone()),
                "edge mask {edge_mask:#08b}, renumbering {index}",
            );
        }
    }
}

#[rstest]
fn test_canonicalize_constitution(canonicalize_context: CanonicalizeContext) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O).with_constraint(AtomConstraintForm::Valence(
                NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(1),
                    ArithExpr::Lit(2),
                ]))),
            )),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::B),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::S),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        dative: vec![
            (
                vec![AtomId(0)],
                AtomId(1),
                DativeBondForm::new(NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
                    ArithExpr::Lit(0),
                    ArithExpr::Lit(1),
                ])))),
            ),
            (vec![AtomId(3)], AtomId(2), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                AtomId(1),
                AtomId(2),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
            (
                AtomId(0),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    });
    let expected_correspondence = molecule_correspondence(&[
        vec![3, 1, 2, 0, 4, 5],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![0, 1],
        vec![0, 1],
    ]);
    let expected = normalize_molecule(molecule.remap(&expected_correspondence)).unwrap();

    assert_eq!(
        canonicalize_constitution_with_options(
            &molecule,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: backend_canonical_branch_order,
            },
        ),
        Ok((expected, expected_correspondence)),
    );
}

#[rstest]
fn test_canonicalize_constitution_excluded_data(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(3))),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(3))),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    });

    let (_, left_correspondence) = canonicalize_constitution_with_options(
        &left,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let (_, right_correspondence) = canonicalize_constitution_with_options(
        &right,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();

    assert_eq!(right_correspondence, left_correspondence);
}

#[rstest]
fn test_canonicalize_constitution_properties(
    initial_class_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let normalized_source = normalize_molecule(initial_class_molecule.clone()).unwrap();
    let (canonical, correspondence) = canonicalize_constitution_with_options(
        &initial_class_molecule,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let acted = normalize_molecule(initial_class_molecule.remap(&correspondence)).unwrap();
    let inverse = correspondence.reverse();

    assert_eq!(acted, canonical);
    assert!(initial_class_molecule.equiv_under(&canonical, &correspondence));
    assert!(canonical.equiv_under(&normalized_source, &inverse));
    assert_eq!(canonical.check_integrity(), Ok(()));

    let (canonical_again, _) = canonicalize_constitution_with_options(
        &canonical,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
    let canonical_again_incidence = canonical_again.incidence_graph(IncidenceLevel::Constitution);
    assert_eq!(
        constitution_comparison_key(
            &canonical,
            &canonical_incidence,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
        constitution_comparison_key(
            &canonical_again,
            &canonical_again_incidence,
            &canonical_again_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>(),
        ),
    );

    let renumbering = reverse_correspondence(&initial_class_molecule);
    let renumbered = initial_class_molecule.remap(&renumbering);
    let (canonical_renumbered, renumbered_correspondence) = canonicalize_constitution_with_options(
        &renumbered,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            prefix_pruning: false,
            branch_order: backend_canonical_branch_order,
        },
    )
    .unwrap();
    let composed = renumbering.compose(&renumbered_correspondence);
    let composed_action = normalize_molecule(initial_class_molecule.remap(&composed)).unwrap();
    let canonical_renumbered_incidence =
        canonical_renumbered.incidence_graph(IncidenceLevel::Constitution);

    assert_eq!(composed_action, canonical_renumbered);
    assert!(initial_class_molecule.equiv_under(&canonical_renumbered, &composed));
    assert_eq!(
        constitution_comparison_key(
            &canonical,
            &canonical_incidence,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
        constitution_comparison_key(
            &canonical_renumbered,
            &canonical_renumbered_incidence,
            &canonical_renumbered_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>(),
        ),
    );

    let (unpruned, _) = canonicalize_constitution_with_options(
        &initial_class_molecule,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: false,
            branch_order: reverse_node_branch_order,
        },
    )
    .unwrap();
    let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
    assert_eq!(
        constitution_comparison_key(
            &canonical,
            &canonical_incidence,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
        constitution_comparison_key(
            &unpruned,
            &unpruned_incidence,
            &unpruned_incidence.graph().node_ids().collect::<Vec<_>>(),
        ),
    );
}

#[rstest]
fn test_canonicalize_constitution_family_minimum(canonicalize_context: CanonicalizeContext) {
    let atoms = vec![AtomForm::from_element(Element::C); 4];
    let cases = [
        (
            "dative",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                dative: vec![
                    (
                        vec![AtomId(0)],
                        AtomId(2),
                        DativeBondForm::new(NumForm::RangeFrom(1)),
                    ),
                    (
                        vec![AtomId(1)],
                        AtomId(3),
                        DativeBondForm::new(NumForm::RangeFrom(1)),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "aromatic",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                aromatic: vec![
                    (
                        vec![AtomId(0), AtomId(1)],
                        AromaticSystemForm::from_electrons(vec![1, 2])
                            .with_charge(NumForm::var("q")),
                    ),
                    (
                        vec![AtomId(2), AtomId(3)],
                        AromaticSystemForm::from_electrons(vec![1, 2])
                            .with_charge(NumForm::var("q")),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "multicenter",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                multicenter: vec![
                    (
                        vec![AtomId(0), AtomId(1)],
                        MulticenterBondForm::from_electrons(vec![2, 1]),
                    ),
                    (
                        vec![AtomId(2), AtomId(3)],
                        MulticenterBondForm::from_electrons(vec![2, 1]),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "noncovalent",
            Molecule::from_entries(MoleculeEntries {
                atoms,
                noncovalent: vec![
                    (
                        AtomId(0),
                        AtomId(1),
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                    (
                        AtomId(2),
                        AtomId(3),
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                ],
                ..Default::default()
            }),
        ),
    ];

    for (family, molecule) in cases {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
        let (entity_keys, incidence_keys) =
            initial_class_keys(&molecule, &incidence_graph).unwrap();
        let classes = rank_initial_classes(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &classes);
        let leaf_candidate =
            |order: &[NodeId]| constitution_candidate(&molecule, &incidence_graph, order).unwrap();
        let expected = exhaustive_minimum(
            &adapter,
            adapter_entity_blocks(&incidence_graph),
            &leaf_candidate,
        );
        let (canonical, correspondence) = canonicalize_constitution_with_options(
            &molecule,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: backend_canonical_branch_order,
            },
        )
        .unwrap();
        let (unpruned, _) = canonicalize_constitution_with_options(
            &molecule,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: reverse_node_branch_order,
            },
        )
        .unwrap();
        let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
        let canonical_order = canonical_incidence.graph().node_ids().collect::<Vec<_>>();
        let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
        let unpruned_order = unpruned_incidence.graph().node_ids().collect::<Vec<_>>();

        assert_eq!(
            constitution_comparison_key(&unpruned, &unpruned_incidence, &unpruned_order),
            Ok(expected.key.clone()),
            "unpruned {family}",
        );
        assert_eq!(
            constitution_comparison_key(&canonical, &canonical_incidence, &canonical_order),
            Ok(expected.key),
            "pruned {family}",
        );
        assert_eq!(unpruned, canonical, "{family}");
        assert!(
            molecule.equiv_under(&canonical, &correspondence),
            "{family}"
        );
        assert_eq!(canonical.check_integrity(), Ok(()), "{family}");

        for (index, atom_images) in permutations(molecule.atoms().count())
            .into_iter()
            .enumerate()
        {
            let mut images = molecule_counts(&molecule).map(|count| (0..count).collect::<Vec<_>>());
            images[0] = atom_images;
            if index % 2 == 1 {
                for family_images in &mut images[1..6] {
                    family_images.reverse();
                }
            }
            let renumbered = molecule.remap(&molecule_correspondence(&images));

            assert_eq!(
                canonicalize_constitution(&renumbered, &canonicalize_context),
                Ok(canonical.clone()),
                "{family}, renumbering {index}",
            );
        }
    }
}

#[rstest]
fn test_canonicalize_constitution_participant_order(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(3), AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![3, 2, 1]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(3), AtomId(5)],
                MulticenterBondForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(1), AtomId(2), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![3, 2, 1]),
            ),
        ],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![
            (
                vec![AtomId(2), AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![3, 1, 2]),
            ),
            (
                vec![AtomId(4), AtomId(5), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1, 3]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(5), AtomId(0), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![3, 1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(4), AtomId(1)],
                MulticenterBondForm::from_electrons(vec![2, 1, 3]),
            ),
        ],
        ..Default::default()
    });

    assert_eq!(
        canonicalize_constitution(&right, &canonicalize_context),
        canonicalize_constitution(&left, &canonicalize_context),
    );
}

#[rstest]
fn test_canonicalize_constitution_contradiction(canonicalize_context: CanonicalizeContext) {
    let selected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        dative: vec![(
            vec![AtomId(0)],
            AtomId(1),
            DativeBondForm::new(NumForm::lit_set([])),
        )],
        ..Default::default()
    });
    let excluded_constraint = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)
            .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([])))],
        ..Default::default()
    });
    let excluded_stereo = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
        )],
        ..Default::default()
    });

    for (location, molecule) in [
        ("selected", selected),
        ("excluded constraint", excluded_constraint),
        ("excluded stereo", excluded_stereo),
    ] {
        assert_eq!(
            canonicalize_constitution(&molecule, &canonicalize_context),
            Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
            "{location}",
        );
    }
}

#[rstest]
fn test_canonical_search_prefix() {
    let source = Graph::new(4, &[]);
    let adapter = direct_graph_adapter(&source);
    let entity_blocks = vec![source.node_ids().collect()];
    let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
        key: order.to_vec(),
        entity_order: order.to_vec(),
    };
    let prefix_worse = |partition: &OrderedPartition, best: &CanonicalCandidate<Vec<NodeId>>| {
        let prefix = partition.fixed_entity_prefix(4);
        prefix.as_slice() > &best.key[..prefix.len()]
    };
    let expected = exhaustive_minimum(&adapter, entity_blocks, &leaf_candidate);
    let actual = canonical_search(
        &adapter,
        &adapter.classes,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            prefix_pruning: true,
            branch_order: node_branch_order,
        },
        &leaf_candidate,
        &prefix_worse,
    );

    assert_eq!(actual.candidate.key, expected.key);
    assert_ne!(actual.stats.prefix_pruned_branches, 0);
}

#[rstest]
#[case::order_four(4)]
fn test_canonical_search_exhaustive(#[case] node_count: usize) {
    let endpoint_pairs = (0..node_count as u32)
        .flat_map(|first| ((first + 1)..node_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
        let edges = endpoint_pairs
            .iter()
            .enumerate()
            .filter_map(|(position, &edge)| ((edge_mask >> position) & 1 == 1).then_some(edge))
            .collect::<Vec<_>>();
        let source = Graph::new(node_count, &edges);
        let adapter = direct_graph_adapter(&source);
        let entity_blocks = vec![source.node_ids().collect()];
        let leaf_candidate = |order: &[NodeId]| {
            let mut positions = vec![0_u32; node_count];
            for (position, node) in order.iter().enumerate() {
                positions[node.index()] = position as u32;
            }
            let mut mapped_edges = source
                .edge_ids()
                .map(|edge| {
                    let [first, second] = source.edge_endpoints(edge);
                    let first = positions[first.index()];
                    let second = positions[second.index()];
                    [first.min(second), first.max(second)]
                })
                .collect::<Vec<_>>();
            mapped_edges.sort_unstable();
            CanonicalCandidate {
                key: mapped_edges,
                entity_order: order.to_vec(),
            }
        };
        let no_prefix = |_: &OrderedPartition, _: &CanonicalCandidate<_>| false;
        let expected = exhaustive_minimum(&adapter, entity_blocks, &leaf_candidate);

        for options in [
            CanonicalSearchOptions {
                automorphism_pruning: false,
                prefix_pruning: false,
                branch_order: reverse_node_branch_order,
            },
            CanonicalSearchOptions {
                automorphism_pruning: true,
                prefix_pruning: false,
                branch_order: backend_canonical_branch_order,
            },
        ] {
            assert_eq!(
                canonical_search(
                    &adapter,
                    &adapter.classes,
                    AutomorphismAlgorithm::Nauty,
                    options,
                    &leaf_candidate,
                    &no_prefix,
                )
                .candidate
                .key,
                expected.key,
                "edge mask {edge_mask:#08b}",
            );
        }
    }
}

#[rstest]
#[case::topology(IncidenceLevel::Topology)]
#[case::constitution(IncidenceLevel::Constitution)]
#[case::full(IncidenceLevel::Full)]
fn test_colored_encoding_dense_remapping_equivalence(#[case] level: IncidenceLevel) {
    let entries = encoding_entries();
    let complete = Molecule::from_entries(entries.clone());
    let molecule = Molecule::from_entries(project_entries(entries, level));
    let remapped = molecule.remap(&reverse_correspondence(&molecule));

    assert!(colored_encoding_equivalent(&complete, &molecule, level));
    assert_eq!(
        colored_encoding_equivalent(&molecule, &remapped, level),
        explicitly_dense_equivalent(&molecule, &remapped),
    );
    assert!(explicitly_dense_equivalent(&molecule, &remapped));

    let mut distinguished = remapped;
    distinguished.atom_mut(AtomId(0)).attributes.element = ElementForm::Lit(Element::O);
    assert_eq!(
        colored_encoding_equivalent(&molecule, &distinguished, level),
        explicitly_dense_equivalent(&molecule, &distinguished),
    );
    assert!(!explicitly_dense_equivalent(&molecule, &distinguished));
}

#[rstest]
#[case::order_four(4)]
fn test_colored_encoding_exhaustive_graph_domain(#[case] atom_count: usize) {
    let endpoint_pairs = (0..atom_count as u32)
        .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
        let bonds = endpoint_pairs
            .iter()
            .enumerate()
            .filter(|(position, _)| (edge_mask >> position) & 1 == 1)
            .map(|(_, &[first, second])| (AtomId(first), AtomId(second), BondForm::from_order(1)))
            .collect();
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); atom_count],
            bonds,
            ..Default::default()
        });
        let remapped = molecule.remap(&reverse_correspondence(&molecule));
        assert_eq!(
            colored_encoding_equivalent(&molecule, &remapped, IncidenceLevel::Topology),
            explicitly_dense_equivalent(&molecule, &remapped),
            "edge mask {edge_mask:#08b}",
        );

        let mut distinguished = remapped;
        distinguished.atom_mut(AtomId(0)).attributes.element = ElementForm::Lit(Element::O);
        assert_eq!(
            colored_encoding_equivalent(&molecule, &distinguished, IncidenceLevel::Topology,),
            explicitly_dense_equivalent(&molecule, &distinguished),
            "edge mask {edge_mask:#08b}",
        );
    }
}

#[rstest]
fn test_initial_classes_error() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C).with_charge(NumForm::LitSet(Box::default()))],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);

    assert_eq!(
        initial_classes(&molecule, &incidence_graph),
        Err(Contradiction)
    );
}

#[rstest]
#[case::integrity(
        MoleculeCanonicalizeError::from(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        MoleculeCanonicalizeError::Integrity(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
#[case::contradiction(
    MoleculeCanonicalizeError::from(Contradiction),
    MoleculeCanonicalizeError::Contradiction(Contradiction)
)]
fn test_molecule_canonicalize_error_from(
    #[case] actual: MoleculeCanonicalizeError,
    #[case] expected: MoleculeCanonicalizeError,
) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case::topology(CanonicalizeLevel::Topology)]
#[case::constitution(CanonicalizeLevel::Constitution)]
#[case::structure(CanonicalizeLevel::Structure)]
#[case::full(CanonicalizeLevel::Full)]
fn test_canonicalize_checked_reaction_span_by(
    canonicalize_context: CanonicalizeContext,
    #[case] level: CanonicalizeLevel,
) {
    let source = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
        ],
        bonds: vec![(
            AtomId(1),
            AtomId(2),
            EntitySpan::Removed(BondForm::from_order(1)),
        )],
        ..Default::default()
    });
    let expected = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![(
            AtomId(0),
            AtomId(1),
            EntitySpan::Removed(BondForm::from_order(1)),
        )],
        ..Default::default()
    });

    assert_eq!(
        canonicalize_checked_reaction_span_by(&source, level, &canonicalize_context,),
        Ok(expected.clone()),
    );
    assert_eq!(
        canonicalize_reaction_span_by(&source, level, &canonicalize_context),
        Ok(expected),
    );
}

#[rstest]
#[case::integrity(
        ReactionSpanCanonicalizeError::from(ReactionSpanIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionSpanCanonicalizeError::Integrity(
            ReactionSpanIntegrityError::InvalidReference {
                entity: Entity::Atom(AtomId(1)),
            },
        ),
    )]
#[case::contradiction(
    ReactionSpanCanonicalizeError::from(Contradiction),
    ReactionSpanCanonicalizeError::Contradiction(Contradiction)
)]
fn test_reaction_span_canonicalize_error_from(
    #[case] actual: ReactionSpanCanonicalizeError,
    #[case] expected: ReactionSpanCanonicalizeError,
) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case::integrity(
        ReactionCanonicalizeError::from(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionCanonicalizeError::Integrity(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
#[case::contradiction(
    ReactionCanonicalizeError::from(Contradiction),
    ReactionCanonicalizeError::Contradiction(Contradiction)
)]
fn test_reaction_canonicalize_error_from(
    #[case] actual: ReactionCanonicalizeError,
    #[case] expected: ReactionCanonicalizeError,
) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case::topology(StructuralDomainPosition::TOPOLOGY, 0)]
#[case::non_stereo(StructuralDomainPosition::NON_STEREO, 1)]
#[case::stereo(StructuralDomainPosition::STEREO, 2)]
fn test_structural_domain_position(
    #[case] position: StructuralDomainPosition,
    #[case] expected: u16,
) {
    assert_eq!(position.0, expected);
}

#[rstest]
#[case::atom(
    EntityBlockPosition::ATOM,
    EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 0)
)]
#[case::bond(
    EntityBlockPosition::BOND,
    EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 1)
)]
#[case::dative_bond(
    EntityBlockPosition::DATIVE_BOND,
    EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 0)
)]
#[case::aromatic_system(
    EntityBlockPosition::AROMATIC_SYSTEM,
    EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 1)
)]
#[case::multicenter_bond(
    EntityBlockPosition::MULTICENTER_BOND,
    EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 2)
)]
#[case::noncovalent_bond(
    EntityBlockPosition::NONCOVALENT_BOND,
    EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 3)
)]
#[case::stereo_atom(
    EntityBlockPosition::STEREO_ATOM,
    EntityBlockPosition::new(StructuralDomainPosition::STEREO, 0)
)]
#[case::stereo_bond(
    EntityBlockPosition::STEREO_BOND,
    EntityBlockPosition::new(StructuralDomainPosition::STEREO, 1)
)]
fn test_entity_block_position(
    #[case] position: EntityBlockPosition,
    #[case] expected: EntityBlockPosition,
) {
    assert_eq!(position, expected);
}

#[rstest]
#[case::topology_slots(EntityBlockPosition::ATOM, EntityBlockPosition::BOND)]
#[case::topology_before_non_stereo(
    EntityBlockPosition::new(StructuralDomainPosition::TOPOLOGY, 2),
    EntityBlockPosition::DATIVE_BOND
)]
#[case::non_stereo_slots(EntityBlockPosition::DATIVE_BOND, EntityBlockPosition::AROMATIC_SYSTEM)]
#[case::non_stereo_before_stereo(
    EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 4),
    EntityBlockPosition::STEREO_ATOM
)]
#[case::stereo_slots(EntityBlockPosition::STEREO_ATOM, EntityBlockPosition::STEREO_BOND)]
fn test_entity_block_position_cmp(
    #[case] lhs: EntityBlockPosition,
    #[case] rhs: EntityBlockPosition,
) {
    assert_eq!(lhs.cmp(&rhs), Ordering::Less);
}

#[rstest]
#[case::atom(
    ConstraintBlockPosition::ATOM,
    ConstraintBlockPosition::Inline(EntityBlockPosition::ATOM)
)]
#[case::bond(
    ConstraintBlockPosition::BOND,
    ConstraintBlockPosition::Inline(EntityBlockPosition::BOND)
)]
#[case::dative_bond(
    ConstraintBlockPosition::DATIVE_BOND,
    ConstraintBlockPosition::Inline(EntityBlockPosition::DATIVE_BOND)
)]
#[case::aromatic_system(
    ConstraintBlockPosition::AROMATIC_SYSTEM,
    ConstraintBlockPosition::Inline(EntityBlockPosition::AROMATIC_SYSTEM)
)]
#[case::multicenter_bond(
    ConstraintBlockPosition::MULTICENTER_BOND,
    ConstraintBlockPosition::Inline(EntityBlockPosition::MULTICENTER_BOND)
)]
#[case::noncovalent_bond(
    ConstraintBlockPosition::NONCOVALENT_BOND,
    ConstraintBlockPosition::Inline(EntityBlockPosition::NONCOVALENT_BOND)
)]
#[case::stereo_atom(
    ConstraintBlockPosition::STEREO_ATOM,
    ConstraintBlockPosition::Inline(EntityBlockPosition::STEREO_ATOM)
)]
#[case::stereo_bond(
    ConstraintBlockPosition::STEREO_BOND,
    ConstraintBlockPosition::Inline(EntityBlockPosition::STEREO_BOND)
)]
#[case::molecule(ConstraintBlockPosition::MOLECULE, ConstraintBlockPosition::Molecule)]
fn test_constraint_block_position(
    #[case] position: ConstraintBlockPosition,
    #[case] expected: ConstraintBlockPosition,
) {
    assert_eq!(position, expected);
}

#[rstest]
#[case::unchanged(SpanTagPosition::UNCHANGED, 0)]
#[case::added(SpanTagPosition::ADDED, 1)]
#[case::removed(SpanTagPosition::REMOVED, 2)]
#[case::modified(SpanTagPosition::MODIFIED, 3)]
fn test_span_tag_position(#[case] position: SpanTagPosition, #[case] expected: u16) {
    assert_eq!(position.0, expected);
}

#[rstest]
#[case::boolean(CanonicalKeyValue::Bool(false), 1)]
#[case::unsigned(CanonicalKeyValue::Unsigned(0), 2)]
#[case::signed(CanonicalKeyValue::Signed(0), 3)]
#[case::text(CanonicalKeyValue::Text(String::new()), 4)]
#[case::sequence(CanonicalKeyValue::Sequence(Vec::new()), 5)]
#[case::product(CanonicalKeyValue::Product(Vec::new()), 6)]
#[case::variant(CanonicalKeyValue::Variant(VariantKey {
        position: VariantPosition(0),
        fields: Vec::new(),
    }), 7)]
#[case::span(CanonicalKeyValue::Span(SpanKey {
        position: SpanTagPosition::UNCHANGED,
        values: Vec::new(),
    }), 8)]
fn test_canonical_key_value_position(#[case] value: CanonicalKeyValue, #[case] expected: u16) {
    assert_eq!(value.position(), expected);
}

#[rstest]
#[case::position_precedes_payload(
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(10),
        },
        PositionedKey {
            position: FieldPosition(1),
            value: CanonicalKeyValue::Signed(-10),
        },
        Ordering::Less,
    )]
#[case::payload_breaks_position_tie(
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(10),
        },
        PositionedKey {
            position: FieldPosition(0),
            value: CanonicalKeyValue::Signed(-10),
        },
        Ordering::Greater,
    )]
fn test_positioned_key_cmp(
    #[case] lhs: FieldKey,
    #[case] rhs: FieldKey,
    #[case] expected: Ordering,
) {
    assert_eq!(lhs.cmp(&rhs), expected);
}

#[rstest]
#[case::span_tag_precedes_value(
        SpanKey {
            position: SpanTagPosition::UNCHANGED,
            values: vec![CanonicalKeyValue::Signed(10)],
        },
        SpanKey {
            position: SpanTagPosition::ADDED,
            values: vec![CanonicalKeyValue::Signed(-10)],
        },
        Ordering::Less,
    )]
#[case::lhs_precedes_rhs(
        SpanKey {
            position: SpanTagPosition::MODIFIED,
            values: vec![CanonicalKeyValue::Signed(0), CanonicalKeyValue::Signed(10)],
        },
        SpanKey {
            position: SpanTagPosition::MODIFIED,
            values: vec![CanonicalKeyValue::Signed(1), CanonicalKeyValue::Signed(-10)],
        },
        Ordering::Less,
    )]
fn test_span_key_cmp(#[case] lhs: SpanKey, #[case] rhs: SpanKey, #[case] expected: Ordering) {
    assert_eq!(lhs.cmp(&rhs), expected);
}

#[rstest]
#[case::unchanged(
        EntitySpan::Unchanged(1),
        SpanTagPosition::UNCHANGED,
        vec![CanonicalKeyValue::Signed(1)],
    )]
#[case::added(
        EntitySpan::Added(1),
        SpanTagPosition::ADDED,
        vec![CanonicalKeyValue::Signed(1)],
    )]
#[case::removed(
        EntitySpan::Removed(1),
        SpanTagPosition::REMOVED,
        vec![CanonicalKeyValue::Signed(1)],
    )]
#[case::modified(
        EntitySpan::Modified { lhs: 1, rhs: 2 },
        SpanTagPosition::MODIFIED,
        vec![CanonicalKeyValue::Signed(1), CanonicalKeyValue::Signed(2)],
    )]
fn test_entity_span_key(
    #[case] span: EntitySpan<i64>,
    #[case] position: SpanTagPosition,
    #[case] values: Vec<CanonicalKeyValue>,
) {
    assert_eq!(
        entity_span_key(&span, |value| CanonicalKeyValue::Signed(*value)),
        CanonicalKeyValue::Span(SpanKey { position, values }),
    );
}

#[rstest]
#[case::absent_extension_preserves_key(None, Ordering::Equal)]
#[case::present_extension_appends_field(
        Some(PositionedKey {
            position: FieldPosition(8),
            value: CanonicalKeyValue::Unsigned(1),
        }),
        Ordering::Less,
    )]
fn test_canonical_key_value_append_only_extension(
    #[case] extension: Option<FieldKey>,
    #[case] expected: Ordering,
) {
    let original = CanonicalKeyValue::Product(vec![PositionedKey {
        position: FieldPosition(0),
        value: CanonicalKeyValue::Signed(1),
    }]);
    let mut extended_fields = match &original {
        CanonicalKeyValue::Product(fields) => fields.clone(),
        _ => unreachable!(),
    };
    extended_fields.extend(extension);
    let extended = CanonicalKeyValue::Product(extended_fields);

    assert_eq!(original.cmp(&extended), expected);
}

#[rstest]
#[case::entity_blocks_precede_constraints(
        CanonicalComparisonKey {
            entity_blocks: vec![PositionedKey {
                position: EntityBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(10),
            }],
        },
        CanonicalComparisonKey {
            entity_blocks: vec![PositionedKey {
                position: EntityBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(1),
            }],
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        Ordering::Less,
    )]
#[case::constraints_break_entity_tie(
        CanonicalComparisonKey {
            entity_blocks: Vec::new(),
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::ATOM,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        CanonicalComparisonKey {
            entity_blocks: Vec::new(),
            constraints: vec![PositionedKey {
                position: ConstraintBlockPosition::BOND,
                value: CanonicalKeyValue::Unsigned(0),
            }],
        },
        Ordering::Less,
    )]
fn test_canonical_comparison_key_cmp(
    #[case] lhs: CanonicalComparisonKey,
    #[case] rhs: CanonicalComparisonKey,
    #[case] expected: Ordering,
) {
    assert_eq!(lhs.cmp(&rhs), expected);
}

#[rstest]
#[case::non_stereo_before_stereo(
    ConstraintBlockPosition::Inline(EntityBlockPosition::new(
        StructuralDomainPosition::NON_STEREO,
        4,
    )),
    ConstraintBlockPosition::STEREO_ATOM
)]
#[case::stereo_before_molecule(
    ConstraintBlockPosition::STEREO_BOND,
    ConstraintBlockPosition::MOLECULE
)]
fn test_constraint_block_position_cmp(
    #[case] lhs: ConstraintBlockPosition,
    #[case] rhs: ConstraintBlockPosition,
) {
    assert_eq!(lhs.cmp(&rhs), Ordering::Less);
}

#[rstest]
#[case::local_slot(
    RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 0),
    RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 1)
)]
#[case::entity_slot(
    RelationalConstraintPosition::new(EntityBlockPosition::DATIVE_BOND, 7),
    RelationalConstraintPosition::new(EntityBlockPosition::AROMATIC_SYSTEM, 0)
)]
#[case::domain(
    RelationalConstraintPosition::new(
        EntityBlockPosition::new(StructuralDomainPosition::NON_STEREO, 4),
        0,
    ),
    RelationalConstraintPosition::new(EntityBlockPosition::STEREO_ATOM, 0)
)]
fn test_relational_constraint_position_cmp(
    #[case] lhs: RelationalConstraintPosition,
    #[case] rhs: RelationalConstraintPosition,
) {
    assert_eq!(lhs.cmp(&rhs), Ordering::Less);
}
