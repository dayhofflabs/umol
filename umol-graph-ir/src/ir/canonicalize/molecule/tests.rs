//! Unit and bounded-exhaustive tests for molecule canonicalization.

use std::collections::BTreeSet;
use std::iter;

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_perm::Permutation;

use super::super::tests::{
    canonicalize_context, exhaustive_minimum, explicitly_dense_equivalent, initial_color_molecule,
    para_stereo_canonicalization_molecule, permutations, reverse_remapping,
    stereo_atom_canonicalization_molecule, stereo_bond_canonicalization_molecule,
    stereo_constraint_canonicalization_molecule,
};
use super::super::*;
use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId,
    Constraint, Constraints, DativeBondForm, DativeBondId, Entity, IncidenceLevel,
    MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, NoncovalentBondId, StereoAtomConstraintForm, StereoAtomForm, StereoAtomId,
    StereoBondConstraintForm, StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandPair, Stereogenicity, StereogenicityForm, Topicity,
    TopicityForm, TopicityRelationForm,
};

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

fn initial_colors(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
) -> Result<InitialColors, Contradiction> {
    let (entity_keys, incidence_keys) = initial_color_keys(molecule, incidence_graph)?;
    Ok(rank_initial_colors(&entity_keys, &incidence_keys))
}

fn leaf_candidate_at(
    molecule: &Molecule,
    incidence_graph: &IncidenceGraph,
    level: DescriptionLevel,
    order: &[NodeId],
) -> CanonicalCandidate<CanonicalComparisonKey> {
    let entity_order = entity_order(molecule, incidence_graph, level, order);
    CanonicalCandidate {
        key: leaf_key(molecule, incidence_graph, level, &entity_order),
        entity_order,
    }
}

#[fixture]
fn symmetric_stereo_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
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
fn stereo_site_exchange_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(AtomForm::from_element)
        .collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
        ],
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

fn selected_structure_key(molecule: &Molecule) -> CanonicalComparisonKey {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    leaf_candidate_at(
        molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        &order,
    )
    .key
}

#[rstest]
#[case::topology(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
)]
#[case::constitution(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::B),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
)]
#[case::structure(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
)]
#[case::full(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
            AtomForm::from_element(Element::C),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
)]
fn test_molecule_canonicalize_description_level(
    canonicalize_context: CanonicalizeContext,
    #[case] source: Molecule,
    #[case] expected: DescriptionLevel,
) {
    let forced_full = source
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Full, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("representative molecule canonicalizes at Full");
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("representative molecule canonicalizes by its effective level");
    let (with_remapping, remapping) = source
        .clone()
        .tracked_canonicalize(&canonicalize_context)
        .expect("representative molecule canonicalizes with a remapping");
    let remapped = source.clone().remap(&reverse_remapping(&source));

    assert_eq!(molecule_canonicalize_level(&source), expected);
    assert_eq!(molecule_canonicalize_level(&remapped), expected);
    assert_eq!(canonical, forced_full);
    assert_eq!(with_remapping, canonical);
    assert_eq!(
        source.clone().remap(&remapping).reframe(),
        Ok(canonical.clone())
    );
    assert_eq!(
        source.canonical_hash(&canonicalize_context),
        canonical.canonical_hash(&canonicalize_context),
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
        stereo_refinement_descriptor(
            7,
            &frame.act(&ligands),
            &configuration
                .apply(frame)
                .expect("the permutation is a parent-group action of the form's kind"),
        ),
    );
}

#[rstest]
#[case::tetrahedral_short(StereoKind::Tetrahedral, 3)]
#[case::tetrahedral_long(StereoKind::Tetrahedral, 5)]
#[case::cis_trans_short(StereoKind::CisTrans, 3)]
#[case::cis_trans_long(StereoKind::CisTrans, 5)]
#[case::octahedral_short(StereoKind::Octahedral, 5)]
#[case::octahedral_long(StereoKind::Octahedral, 7)]
fn test_stereo_refinement_descriptor_error(#[case] kind: StereoKind, #[case] degree: usize) {
    let ligands = (0..degree)
        .map(|class| (class as u32, StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let configuration = StereoConfigurationForm::kinded(kind, 0u32);

    assert_eq!(
        stereo_refinement_descriptor(0, &ligands, &configuration),
        Err(Contradiction),
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
    let (entity_keys, incidence_keys) =
        initial_color_keys(&molecule, &incidence_graph).expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
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
    let (entity_keys, incidence_keys) =
        initial_color_keys(&molecule, &incidence_graph).expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (_, rounds) =
        structure_partition(&molecule, &incidence_graph, &adapter, &entity_keys, true)
            .expect("fixed molecule has a structure partition");

    assert_eq!(rounds, 1);
}

#[rstest]
fn test_structure_partition_distinct_stereo(stereo_atom_canonicalization_molecule: Molecule) {
    let molecule = stereo_atom_canonicalization_molecule;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) =
        initial_color_keys(&molecule, &incidence_graph).expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
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
    let atom_images = [1, 0, 6, 8, 7, 9, 2, 4, 3, 5, 13, 11, 10, 12].map(AtomId);
    let mut images = molecule_counts(&para_stereo_canonicalization_molecule)
        .map(|count| (0..count).collect::<Vec<_>>());
    images[0] = atom_images.iter().map(|id| id.index()).collect();
    let remapping = molecule_remapping(&images);
    let renumbered = para_stereo_canonicalization_molecule
        .clone()
        .remap(&remapping);
    let context = CanonicalizeContext {
        para_stereo: true,
        ..canonicalize_context
    };
    let (canonical, remapping) = canonicalize_structure(
        para_stereo_canonicalization_molecule.clone(),
        &context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .expect("fixed molecule canonicalizes");
    let (unpruned, unpruned_remapping) = canonicalize_structure(
        para_stereo_canonicalization_molecule.clone(),
        &context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
    )
    .expect("fixed molecule canonicalizes without orbit pruning");

    assert_eq!(unpruned, canonical);
    assert_eq!(
        para_stereo_canonicalization_molecule
            .clone()
            .remap(&remapping)
            .reframe(),
        Ok(canonical.clone()),
    );
    assert_eq!(
        para_stereo_canonicalization_molecule
            .remap(&unpruned_remapping)
            .reframe(),
        Ok(canonical.clone()),
    );
    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &context)
            .map(|(canonical, _)| canonical),
        Ok(canonical)
    );
}

#[rstest]
fn test_structure_comparison_key(stereo_atom_canonicalization_molecule: Molecule) {
    let mut constrained = stereo_atom_canonicalization_molecule.clone();
    constrained
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints = StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    ))
    .into();
    let incidence_graph =
        stereo_atom_canonicalization_molecule.incidence_graph(IncidenceLevel::Full);
    let constrained_incidence_graph = constrained.incidence_graph(IncidenceLevel::Full);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();

    assert_eq!(
        leaf_candidate_at(
            &stereo_atom_canonicalization_molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            &order
        )
        .key,
        leaf_candidate_at(
            &constrained,
            &constrained_incidence_graph,
            DescriptionLevel::Structure,
            &order
        )
        .key,
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_atom(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_remapping(&[
        vec![4, 2, 0, 3, 1],
        vec![2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
        Vec::new(),
    ]);
    let renumbered = stereo_atom_canonicalization_molecule
        .clone()
        .remap(&renumbering);
    let canonical = stereo_atom_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
    let mut opposite = stereo_atom_canonicalization_molecule.clone();
    opposite
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32);

    assert_ne!(
        stereo_atom_canonicalization_molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        opposite
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_bond(
    canonicalize_context: CanonicalizeContext,
    stereo_bond_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_remapping(&[
        vec![5, 3, 1, 4, 2, 0],
        vec![4, 2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
    ]);
    let renumbered = stereo_bond_canonicalization_molecule
        .clone()
        .remap(&renumbering);
    let canonical = stereo_bond_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
    let mut editor = stereo_atom_canonicalization_molecule.edit();
    editor
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints = constraint.clone().into();
    *editor.constraints_mut() =
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
    let source = editor.finish().unwrap();
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");

    assert_eq!(
        reframed
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        source
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
    let mut editor = stereo_bond_canonicalization_molecule.edit();
    editor
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .constraints = constraint.clone().into();
    *editor.constraints_mut() =
        Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, constraint).into();
    let source = editor.finish().unwrap();
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");
    let left = source
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");
    let right = reframed
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("reframed molecule canonicalizes");

    assert_eq!(
        selected_structure_key(&right),
        selected_structure_key(&left)
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
fn test_canonicalize_molecule_by_effective_regression(
    canonicalize_context: CanonicalizeContext,
    #[case] source: &str,
) {
    let source = Molecule::from_edn_str(source).expect("retained molecule parses");
    let expected = source.clone();

    for level in [
        DescriptionLevel::Topology,
        DescriptionLevel::Constitution,
        DescriptionLevel::Structure,
        DescriptionLevel::Full,
    ] {
        let canonical = source
            .clone()
            .tracked_canonicalize_at(level, &canonicalize_context)
            .map(|(canonical, _)| canonical)
            .expect("retained molecule canonicalizes");

        assert_eq!(canonical, expected);
        assert_eq!(
            canonical.canonical_key_at(level, &canonicalize_context),
            expected.canonical_key_at(level, &canonicalize_context),
        );
    }

    let (canonical, remapping) = source
        .clone()
        .tracked_canonicalize(&canonicalize_context)
        .expect("retained molecule canonicalizes with a remapping");
    let transported = source.clone().remap(&remapping);

    assert_eq!(canonical, expected);
    assert_eq!(transported, expected);
    assert!(source.framed_eq_under(&canonical, &remapping));
}

#[rstest]
fn test_canonical_key_by_topology(canonicalize_context: CanonicalizeContext) {
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
        left.canonical_key_at(DescriptionLevel::Topology, &canonicalize_context),
        right.canonical_key_at(DescriptionLevel::Topology, &canonicalize_context),
    );
    assert_ne!(
        left.canonical_key_at(DescriptionLevel::Constitution, &canonicalize_context),
        right.canonical_key_at(DescriptionLevel::Constitution, &canonicalize_context),
    );
}

#[rstest]
fn test_canonical_key_by_constitution(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let constitution = Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        ..Default::default()
    });

    assert_eq!(
        constitution.canonical_key_at(DescriptionLevel::Constitution, &canonicalize_context),
        stereo_atom_canonicalization_molecule
            .canonical_key_at(DescriptionLevel::Constitution, &canonicalize_context),
    );
    assert_ne!(
        constitution.canonical_key_at(DescriptionLevel::Structure, &canonicalize_context),
        stereo_atom_canonicalization_molecule
            .canonical_key_at(DescriptionLevel::Structure, &canonicalize_context),
    );
}

#[rstest]
fn test_canonical_key_by_structure(canonicalize_context: CanonicalizeContext) {
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
        plain.canonical_key_at(DescriptionLevel::Structure, &canonicalize_context),
        constrained.canonical_key_at(DescriptionLevel::Structure, &canonicalize_context),
    );
    assert_ne!(
        plain.canonical_key_at(DescriptionLevel::Full, &canonicalize_context),
        constrained.canonical_key_at(DescriptionLevel::Full, &canonicalize_context),
    );
}

#[rstest]
#[case::stereo_atom(stereo_atom_canonicalization_molecule(), [1, 1])]
#[case::stereo_bond(stereo_bond_canonicalization_molecule(), [1, 1])]
#[case::symmetric_stereo(symmetric_stereo_canonicalization_molecule(), [24, 24])]
#[case::para_stereo(para_stereo_canonicalization_molecule(), [8, 8])]
fn test_molecule_canonical_key_at(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
    #[case] expected_visited_leaves: [usize; 2],
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) =
        initial_color_keys(&molecule, &incidence_graph).expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = constitution_partition_descriptors(&adapter, &entity_keys, &incidence_graph);
    let leaf_candidate = |order: &[NodeId]| {
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        )
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let unpruned = canonical_search(
        &adapter,
        &descriptors,
        canonicalize_context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::NodeId,
        },
        &leaf_candidate,
    );
    let pruned = canonical_search_with_generator_filter(
        &adapter,
        &descriptors,
        canonicalize_context.automorphism_algorithm,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
        &leaf_candidate,
        &filter_generators,
    );

    assert_eq!(
        molecule.canonical_key_at(DescriptionLevel::Structure, &canonicalize_context),
        Ok(unpruned.candidate.key.clone()),
    );
    assert_eq!(pruned.candidate.key, unpruned.candidate.key);
    assert_eq!(
        [unpruned.stats.visited_leaves, pruned.stats.visited_leaves],
        expected_visited_leaves,
    );
}

#[rstest]
fn test_canonicalize_structure_selected_layer(
    initial_color_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let canonical = initial_color_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");
    let canonical_again = canonical
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
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
    let canonical = symmetric_stereo_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    for rank in 0..(1..=5).product() {
        let permutation = Permutation::unrank(5, rank);
        let atom_images = (0..5)
            .map(|index| AtomId(permutation.apply(index) as u32))
            .collect::<Vec<_>>();
        let remapping = molecule_remapping(&[
            atom_images.iter().map(|id| id.index()).collect(),
            vec![0, 1, 2, 3],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![0],
            Vec::new(),
        ]);
        let renumbered = symmetric_stereo_canonicalization_molecule
            .clone()
            .remap(&remapping);
        let (actual, canonical_remapping) = canonicalize_structure(
            renumbered.clone(),
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                branch_order: BranchOrdering::BackendCanonical,
            },
        )
        .expect("renumbered molecule canonicalizes");

        assert_eq!(actual, canonical, "rank {rank}");
        assert!(
            renumbered.framed_eq_under(&actual, &canonical_remapping),
            "rank {rank}",
        );
        assert_eq!(
            renumbered.remap(&canonical_remapping).reframe(),
            Ok(actual),
            "rank {rank}",
        );
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
    let (entity_keys, incidence_keys) = initial_color_keys(
        &symmetric_stereo_canonicalization_molecule,
        &incidence_graph,
    )
    .expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
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
        leaf_candidate_at(
            &symmetric_stereo_canonicalization_molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        )
    };
    let expected = exhaustive_minimum(
        &adapter,
        adapter_entity_blocks(&incidence_graph),
        &leaf_candidate,
    );
    let expected_remapping = remapping_from_order(
        &symmetric_stereo_canonicalization_molecule,
        &incidence_graph,
        &expected.entity_order,
    );
    let expected_molecule = symmetric_stereo_canonicalization_molecule
        .clone()
        .remap(&expected_remapping)
        .reframe()
        .expect("exhaustive minimum preserves compatible frames");

    for options in [
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::NodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
    ] {
        let actual = canonical_search(&adapter, &descriptors, algorithm, options, &leaf_candidate);
        assert_eq!(actual.candidate.key, expected.key, "{options:?}");
    }

    for options in [
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::NodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::NodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::ReverseNodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    ] {
        let (canonical, remapping) = canonicalize_structure(
            symmetric_stereo_canonicalization_molecule.clone(),
            &CanonicalizeContext {
                automorphism_algorithm: algorithm,
                ..canonicalize_context
            },
            options,
        )
        .expect("fixed molecule canonicalizes");

        assert_eq!(canonical, expected_molecule, "{options:?}");
        assert_eq!(
            selected_structure_key(&canonical),
            expected.key,
            "{options:?}"
        );
        assert!(
            symmetric_stereo_canonicalization_molecule.framed_eq_under(&canonical, &remapping),
            "{options:?}",
        );
        assert_eq!(
            symmetric_stereo_canonicalization_molecule
                .clone()
                .remap(&remapping)
                .reframe(),
            Ok(canonical),
            "{options:?}",
        );
    }
}

#[rstest]
#[case::nauty(AutomorphismAlgorithm::Nauty)]
fn test_canonicalize_structure_orbit_pruning(
    stereo_site_exchange_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] algorithm: AutomorphismAlgorithm,
) {
    let molecule = stereo_site_exchange_canonicalization_molecule;
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let (entity_keys, incidence_keys) =
        initial_color_keys(&molecule, &incidence_graph).expect("fixed molecule has initial colors");
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let (partition, _) = structure_partition(
        &molecule,
        &incidence_graph,
        &adapter,
        &entity_keys,
        canonicalize_context.para_stereo,
    )
    .expect("fixed molecule has a structure partition");
    let descriptors = partition.cell_indices(adapter.graph().node_count());
    let leaf_candidate = |order: &[NodeId]| {
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        )
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let entity_cells = partition
        .cells
        .iter()
        .filter(|cell| {
            cell.first()
                .is_some_and(|node| node.index() < adapter.source_node_count)
        })
        .cloned()
        .collect();
    let expected = exhaustive_minimum(&adapter, entity_cells, &leaf_candidate);
    let expected_remapping =
        remapping_from_order(&molecule, &incidence_graph, &expected.entity_order);
    let expected_molecule = molecule
        .clone()
        .remap(&expected_remapping)
        .reframe()
        .expect("exhaustive minimum preserves compatible frames");
    let mut stats = Vec::new();

    for options in [
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::NodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::NodeId,
        },
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::ReverseNodeId,
        },
    ] {
        let actual = canonical_search_with_generator_filter(
            &adapter,
            &descriptors,
            algorithm,
            options,
            &leaf_candidate,
            &filter_generators,
        );
        let remapping =
            remapping_from_order(&molecule, &incidence_graph, &actual.candidate.entity_order);
        let canonical = molecule
            .clone()
            .remap(&remapping)
            .reframe()
            .expect("selected minimum preserves compatible frames");

        assert_eq!(actual.candidate.key, expected.key, "{options:?}");
        assert_eq!(canonical, expected_molecule, "{options:?}");
        assert!(
            molecule.framed_eq_under(&canonical, &remapping),
            "{options:?}",
        );
        stats.push(actual.stats);
    }

    assert_eq!(
        stats,
        vec![
            CanonicalSearchStats {
                initial_residual_cell_sizes: vec![2; 10],
                refinement_calls: 3,
                branch_order_calls: 1,
                backend_calls: 1,
                visited_leaves: 2,
                leaf_comparisons: 1,
                orbit_pruned_branches: 0,
            },
            CanonicalSearchStats {
                initial_residual_cell_sizes: vec![2; 10],
                refinement_calls: 3,
                branch_order_calls: 1,
                backend_calls: 1,
                visited_leaves: 2,
                leaf_comparisons: 1,
                orbit_pruned_branches: 0,
            },
            CanonicalSearchStats {
                initial_residual_cell_sizes: vec![2; 10],
                refinement_calls: 2,
                branch_order_calls: 1,
                backend_calls: 1,
                visited_leaves: 1,
                leaf_comparisons: 0,
                orbit_pruned_branches: 1,
            },
            CanonicalSearchStats {
                initial_residual_cell_sizes: vec![2; 10],
                refinement_calls: 2,
                branch_order_calls: 1,
                backend_calls: 1,
                visited_leaves: 1,
                leaf_comparisons: 0,
                orbit_pruned_branches: 1,
            },
        ],
    );
}

#[rstest]
fn test_canonicalize_structure_meso(
    meso_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let remapping = molecule_remapping(&[
        vec![1, 0, 3, 2, 5, 4],
        vec![0, 3, 4, 1, 2],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![1, 0],
        Vec::new(),
    ]);
    let renumbered = meso_canonicalization_molecule.clone().remap(&remapping);

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        meso_canonicalization_molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
#[case::kinded(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))]
#[case::undetermined(StereoConfigurationForm::Undetermined)]
fn test_canonicalize_structure_stereo_frame(
    stereo_constraint_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] configuration: StereoConfigurationForm,
) {
    let mut source = stereo_constraint_canonicalization_molecule;
    source
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");

    assert_eq!(
        reframed
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        source
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

fn rank_paired_initial_colors(
    left: (&[InitialColorKey], &[InitialColorKey]),
    right: (&[InitialColorKey], &[InitialColorKey]),
) -> (InitialColors, InitialColors) {
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
        |entity_keys: &[InitialColorKey], incidence_keys: &[InitialColorKey]| InitialColors {
            entities: entity_keys.iter().map(|key| ranks[key]).collect(),
            incidences: incidence_keys.iter().map(|key| ranks[key]).collect(),
        };

    (rank(left.0, left.1), rank(right.0, right.1))
}

fn colored_encoding_equivalent(left: &Molecule, right: &Molecule, level: IncidenceLevel) -> bool {
    fn key(adapter: &AutomorphismAdapter) -> Vec<u8> {
        adapter.graph().canonical_key(
            |node| {
                let (domain, rank) = match adapter.color(node) {
                    AutomorphismColor::Entity(rank) => (0, rank),
                    AutomorphismColor::Incidence(rank) => (1, rank),
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
        initial_color_keys(left, &left_incidence).unwrap();
    let (right_entity_keys, right_incidence_keys) =
        initial_color_keys(right, &right_incidence).unwrap();
    let (left_colors, right_colors) = rank_paired_initial_colors(
        (&left_entity_keys, &left_incidence_keys),
        (&right_entity_keys, &right_incidence_keys),
    );
    let left_adapter = AutomorphismAdapter::new(&left_incidence, &left_colors);
    let right_adapter = AutomorphismAdapter::new(&right_incidence, &right_colors);

    key(&left_adapter) == key(&right_adapter)
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
            [AtomId(2), AtomId(3)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )],
        ..Default::default()
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
fn test_entity_color_key(
    initial_color_molecule: Molecule,
    #[case] entity: Entity,
    #[case] expected: Vec<FieldPosition>,
) {
    let InitialColorKey::Entity {
        value: CanonicalKeyValue::Product(fields),
        ..
    } = entity_color_key(&initial_color_molecule, entity).unwrap()
    else {
        panic!("entity color key must be a product");
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
fn test_initial_colors(
    initial_color_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected_equal: bool,
) {
    let incidence_graph = initial_color_molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&initial_color_molecule, &incidence_graph).unwrap();
    let lhs_color = colors.entities[incidence_graph.node_of(lhs).index()];
    let rhs_color = colors.entities[incidence_graph.node_of(rhs).index()];

    assert_eq!(lhs_color == rhs_color, expected_equal);
}

#[rstest]
fn test_initial_colors_incidence(initial_color_molecule: Molecule) {
    let incidence_graph = initial_color_molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&initial_color_molecule, &incidence_graph).unwrap();
    let contributions = [1, 2, 2, 1, 1, 2, 1, 2, 2, 1];
    let incidences = incidence_graph
        .incidences()
        .map(|(edge, incidence)| {
            let contribution = edge
                .index()
                .checked_sub(14)
                .and_then(|position| contributions.get(position));
            ((incidence, contribution), colors.incidences[edge.index()])
        })
        .collect::<Vec<_>>();

    for (lhs, lhs_color) in &incidences {
        for (rhs, rhs_color) in &incidences {
            assert_eq!(lhs_color == rhs_color, lhs == rhs);
        }
    }
    for entity_color in &colors.entities {
        for (_, incidence_color) in &incidences {
            assert_ne!(entity_color, incidence_color);
        }
    }
}

#[rstest]
#[case::incomplete_order_one(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    Graph::new(2, &[[0, 1]]),
    vec![AutomorphismColor::Entity(0); 2],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
    ],
    vec![vec![NodeId(0), NodeId(1)]],
)]
#[case::modal_majority(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2).into_concrete()),
            (AtomId(1), AtomId(2), BondForm::from_order(2).into_concrete()),
            (AtomId(0), AtomId(2), BondForm::from_order(1).into_concrete()),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 1], [1, 2], [0, 3], [3, 2]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(1),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
        SubdivisionNodeSource::Node(NodeId(3)),
    ],
    vec![vec![NodeId(0), NodeId(1), NodeId(2)], vec![NodeId(3)]],
)]
#[case::modal_tie(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2).into_concrete()),
            (AtomId(1), AtomId(2), BondForm::from_order(1).into_concrete()),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 3], [3, 1], [1, 2]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(2),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
        SubdivisionNodeSource::Node(NodeId(3)),
    ],
    vec![vec![NodeId(0), NodeId(1), NodeId(2)], vec![NodeId(3)]],
)]
fn test_topology_search_input(
    #[case] molecule: Molecule,
    #[case] expected_graph: Graph,
    #[case] expected_colors: Vec<AutomorphismColor>,
    #[case] expected_node_sources: Vec<SubdivisionNodeSource>,
    #[case] expected_partition: Vec<Vec<NodeId>>,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let (entity_keys, _) = initial_color_keys(&molecule, &incidence_graph).unwrap();
    let entity_colors = rank_initial_colors(&entity_keys, &[]).entities;
    let input = topology_search_input(&molecule, &incidence_graph, &entity_colors);

    assert_eq!(input.adapter.graph, expected_graph);
    assert_eq!(input.adapter.colors, expected_colors);
    assert_eq!(input.adapter.node_sources, expected_node_sources);
    assert_eq!(input.adapter.incidence_nodes, Vec::new());
    assert_eq!(
        input.adapter.source_node_count,
        incidence_graph.entity_count(EntityKind::Atom),
    );
    assert_eq!(input.initial_partition.cells, expected_partition);
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
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Topology,
            &order
        )
        .key,
        CanonicalComparisonKey {
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
        },
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
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Topology,
            &order
        )
        .key,
        leaf_candidate_at(
            &remapped,
            &remapped_incidence_graph,
            DescriptionLevel::Topology,
            &remapped_order
        )
        .key,
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
            noncovalent: vec![([AtomId(0), AtomId(2)],
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
    let key = leaf_candidate_at(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Constitution,
        &order,
    )
    .key;

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
        noncovalent: vec![([AtomId(0), AtomId(3)], NoncovalentBondForm::default())],
        ..Default::default()
    });
    let remapping = molecule_remapping(&[
        vec![3, 1, 0, 2],
        Vec::new(),
        vec![0],
        vec![0],
        vec![0],
        vec![0],
        Vec::new(),
        Vec::new(),
    ]);
    let remapped = molecule.clone().remap(&remapping);
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let remapped_incidence_graph = remapped.incidence_graph(IncidenceLevel::Constitution);
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    order.reverse();
    let remapped_order = order
        .iter()
        .map(|&node| incidence_graph.entity(node))
        .map(|entity| {
            MoleculeCorrespondence::from(&remapping)
                .right_of(entity)
                .expect("dense remapping maps every entity")
        })
        .map(|entity| remapped_incidence_graph.node_of(entity))
        .collect::<Vec<_>>();

    assert_eq!(
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Constitution,
            &order
        )
        .key,
        leaf_candidate_at(
            &remapped,
            &remapped_incidence_graph,
            DescriptionLevel::Constitution,
            &remapped_order
        )
        .key,
    );
}

#[rstest]
fn test_constitution_comparison_key_excluded_data() {
    let molecule = Molecule::from_entries(project_entries(
        encoding_entries(),
        IncidenceLevel::Constitution,
    ));
    let mut excluded = Molecule::from_entries(encoding_entries());
    for index in 0..excluded.atoms().count() {
        excluded
            .atom_mut(AtomId(index as u32))
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::Valence(NumForm::Lit(4)));
    }
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let excluded_incidence_graph = excluded.incidence_graph(IncidenceLevel::Constitution);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let excluded_order = excluded_incidence_graph
        .graph()
        .node_ids()
        .collect::<Vec<_>>();

    assert_eq!(
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Constitution,
            &order
        )
        .key,
        leaf_candidate_at(
            &excluded,
            &excluded_incidence_graph,
            DescriptionLevel::Constitution,
            &excluded_order
        )
        .key,
    );
}

#[rstest]
#[case::dative(Entity::DativeBond(DativeBondId(0)))]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
fn test_remapping_from_order(initial_color_molecule: Molecule, #[case] excluded: Entity) {
    let incidence_graph = initial_color_molecule.incidence_graph(IncidenceLevel::Topology);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let remapping = remapping_from_order(&initial_color_molecule, &incidence_graph, &order);

    assert_eq!(
        MoleculeCorrespondence::from(&remapping).right_of(excluded),
        Some(excluded)
    );
}

#[rstest]
#[case::bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    vec![
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
    ],
)]
#[case::dative_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        dative: vec![
            (vec![AtomId(0)], AtomId(2), DativeBondForm::from_order(1)),
            (vec![AtomId(1)], AtomId(3), DativeBondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::DativeBond(DativeBondId(1)),
        Entity::DativeBond(DativeBondId(0)),
    ],
)]
#[case::aromatic_system(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        aromatic: vec![
            (vec![AtomId(0), AtomId(1)], AromaticSystemForm::default()),
            (vec![AtomId(2), AtomId(3)], AromaticSystemForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::AromaticSystem(AromaticSystemId(1)),
        Entity::AromaticSystem(AromaticSystemId(0)),
    ],
)]
#[case::multicenter_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        multicenter: vec![
            (vec![AtomId(0), AtomId(1)], MulticenterBondForm::default()),
            (vec![AtomId(2), AtomId(3)], MulticenterBondForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::MulticenterBond(MulticenterBondId(1)),
        Entity::MulticenterBond(MulticenterBondId(0)),
    ],
)]
#[case::noncovalent_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        noncovalent: vec![
            ([AtomId(0), AtomId(1)], NoncovalentBondForm::default()),
            ([AtomId(2), AtomId(3)], NoncovalentBondForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1)),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    ],
)]
#[case::stereo_atom(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        stereo_atoms: vec![
            (AtomId(0), Vec::new(), StereoAtomForm::default()),
            (AtomId(1), Vec::new(), StereoAtomForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::StereoAtom(StereoAtomId(1)),
        Entity::StereoAtom(StereoAtomId(0)),
    ],
)]
#[case::stereo_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            ),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
        Entity::StereoBond(StereoBondId(1)),
        Entity::StereoBond(StereoBondId(0)),
    ],
)]
#[case::constitution_level(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::default())],
        stereo_atoms: vec![(AtomId(0), Vec::new(), StereoAtomForm::default())],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    ],
)]
fn test_entity_order(
    #[case] molecule: Molecule,
    #[case] level: DescriptionLevel,
    #[case] expected: Vec<Entity>,
) {
    let incidence_graph = molecule.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    order.reverse();

    assert_eq!(
        entity_order(&molecule, &incidence_graph, level, &order)
            .into_iter()
            .map(|node| incidence_graph.entity(node))
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(1), AtomId(0), BondForm::from_order(2)),
            (AtomId(1), AtomId(2), BondForm::from_order(1).with_charge(-1_i64)),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    EntityBlockPosition::BOND,
    Some(vec![
        positioned_product([
            (0, product([CanonicalKeyValue::Unsigned(0), CanonicalKeyValue::Unsigned(1)])),
            (1, num_form_key(&NumForm::Lit(2))),
            (2, num_form_key(&NumForm::Undetermined)),
            (3, unpaired_electrons_form_key(&UnpairedElectronsForm::default())),
        ]),
        positioned_product([
            (0, product([CanonicalKeyValue::Unsigned(1), CanonicalKeyValue::Unsigned(2)])),
            (1, num_form_key(&NumForm::Lit(1))),
            (2, num_form_key(&NumForm::Lit(-1))),
            (3, unpaired_electrons_form_key(&UnpairedElectronsForm::default())),
        ]),
    ]),
)]
#[case::dative_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        dative: vec![(
            vec![AtomId(1), AtomId(0)],
            AtomId(2),
            DativeBondForm::new(NumForm::RangeFrom(-1)),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    EntityBlockPosition::DATIVE_BOND,
    Some(vec![positioned_product([
        (
            0,
            sequence([CanonicalKeyValue::Unsigned(0), CanonicalKeyValue::Unsigned(1)]),
        ),
        (1, CanonicalKeyValue::Unsigned(2)),
        (2, num_form_key(&NumForm::RangeFrom(-1))),
    ])]),
)]
#[case::aromatic_system(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(NumForm::var("q")),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    EntityBlockPosition::AROMATIC_SYSTEM,
    Some(vec![positioned_product([
        (
            0,
            sequence([CanonicalKeyValue::Unsigned(0), CanonicalKeyValue::Unsigned(2)]),
        ),
        (1, electron_counts_form_key(&ElectronCountsForm::Lit(vec![1, 2]))),
        (2, num_form_key(&NumForm::var("q"))),
        (3, unpaired_electrons_form_key(&UnpairedElectronsForm::default())),
    ])]),
)]
#[case::multicenter_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        multicenter: vec![(
            vec![AtomId(0), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    EntityBlockPosition::MULTICENTER_BOND,
    Some(vec![positioned_product([
        (
            0,
            sequence([CanonicalKeyValue::Unsigned(0), CanonicalKeyValue::Unsigned(2)]),
        ),
        (1, electron_counts_form_key(&ElectronCountsForm::Lit(vec![2, 1]))),
        (2, num_form_key(&NumForm::Undetermined)),
        (3, unpaired_electrons_form_key(&UnpairedElectronsForm::default())),
    ])]),
)]
#[case::noncovalent_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        noncovalent: vec![([AtomId(2), AtomId(0)], NoncovalentBondForm::default())],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    EntityBlockPosition::NONCOVALENT_BOND,
    Some(vec![positioned_product([
        (0, product([CanonicalKeyValue::Unsigned(0), CanonicalKeyValue::Unsigned(2)])),
        (1, noncovalent_bond_kind_form_key(&NoncovalentBondKindForm::Undetermined)),
    ])]),
)]
#[case::stereo_atom(
    stereo_atom_canonicalization_molecule(),
    DescriptionLevel::Structure,
    EntityBlockPosition::STEREO_ATOM,
    Some(vec![positioned_product([
        (0, CanonicalKeyValue::Unsigned(0)),
        (
            1,
            sequence([
                stereo_ligand_key(4, StereoLigandKind::Atom),
                stereo_ligand_key(2, StereoLigandKind::Atom),
                stereo_ligand_key(1, StereoLigandKind::Atom),
                stereo_ligand_key(3, StereoLigandKind::Atom),
            ]),
        ),
        (
            2,
            stereo_configuration_form_key(&StereoConfigurationForm::kinded(
                StereoKind::Tetrahedral,
                0u32,
            )),
        ),
    ])]),
)]
#[case::stereo_atom_below_structure(
    stereo_atom_canonicalization_molecule(),
    DescriptionLevel::Constitution,
    EntityBlockPosition::STEREO_ATOM,
    None
)]
fn test_comparison_key(
    #[case] molecule: Molecule,
    #[case] level: DescriptionLevel,
    #[case] position: EntityBlockPosition,
    #[case] expected: Option<Vec<CanonicalKeyValue>>,
) {
    let key = comparison_key(&molecule, level).expect("fixed molecule has a comparison key");

    assert_eq!(
        key.entity_blocks
            .into_iter()
            .find(|block| block.position == position),
        expected.map(|rows| PositionedKey {
            position,
            value: sequence(rows),
        }),
    );
}

#[rstest]
#[case::topology(initial_color_molecule(), DescriptionLevel::Topology)]
#[case::constitution(initial_color_molecule(), DescriptionLevel::Constitution)]
#[case::structure(initial_color_molecule(), DescriptionLevel::Structure)]
#[case::stereo_atom(stereo_atom_canonicalization_molecule(), DescriptionLevel::Structure)]
#[case::stereo_bond(stereo_bond_canonicalization_molecule(), DescriptionLevel::Structure)]
#[case::full(stereo_constraint_canonicalization_molecule(), DescriptionLevel::Full)]
fn test_leaf_key(#[case] molecule: Molecule, #[case] level: DescriptionLevel) {
    let molecule = molecule.normalize().expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();

    for _ in 0..2 {
        let entity_order = entity_order(&molecule, &incidence_graph, level, &order);
        let remapping = remapping_from_order(&molecule, &incidence_graph, &entity_order);
        let canonical = molecule
            .clone()
            .remap(&remapping)
            .reframe()
            .expect("fixed molecule reframes");

        assert_eq!(
            leaf_key(&molecule, &incidence_graph, level, &entity_order),
            comparison_key(&canonical, level).expect("fixed molecule has a comparison key"),
        );
        order.reverse();
    }
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
#[case::mixed_virtual_ligand_anchor(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::StereoSite,
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            Incidence::StereoLigand(StereoLigandKind::LonePair),
        ],
        5,
    )]
fn test_automorphism_adapter_new(
    #[case] molecule: Molecule,
    #[case] expected_incidences: Vec<Incidence>,
    #[case] expected_occurrence_nodes: usize,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &incidence_graph).unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
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
            adapter.color(adapter_node),
            AutomorphismColor::Entity(colors.entities[node.index()]),
        );
    }
    for edge in source.edge_ids() {
        if let Some(adapter_node) = adapter.node_of(SubdivisionNodeSource::Edge(edge)) {
            assert_eq!(
                adapter.node_source(adapter_node),
                SubdivisionNodeSource::Edge(edge),
            );
            assert_eq!(
                adapter.color(adapter_node),
                AutomorphismColor::Incidence(colors.incidences[edge.index()]),
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
    let colors = initial_colors(&molecule, &incidence_graph).unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);

    assert_eq!(
        adapter.automorphisms(AutomorphismAlgorithm::Nauty),
        AutomorphismAdapterOutput {
            source_orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
            source_canonical_labels: vec![NodeId(0), NodeId(1), NodeId(2)],
            source_generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
        },
    );
}

#[rstest]
#[case::feature_free(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    vec![NodeId(1), NodeId(0), NodeId(2)],
    true,
)]
#[case::distinct_ligands(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    (0..10).map(NodeId::from).collect(),
    true,
)]
#[case::prochiral_explicit_hydrogens(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::stereo_site_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
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
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
#[case::stereo_bond_endpoint_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
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
            [2, 3, 4, 5]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(1), NodeId(0), NodeId(4), NodeId(5), NodeId(2), NodeId(3),
        NodeId(6), NodeId(9), NodeId(10), NodeId(7), NodeId(8), NodeId(11),
    ],
    true,
)]
fn test_retain_stereo_preserving_generators(
    #[case] molecule: Molecule,
    #[case] generator: Vec<NodeId>,
    #[case] expected_retained: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let mut generators = vec![generator.clone()];
    let expected = if expected_retained {
        vec![generator]
    } else {
        Vec::new()
    };

    retain_stereo_preserving_generators(&molecule, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
}

#[rstest]
#[case::order_four(4)]
fn test_ordered_partition_refine_exhaustive_domain(
    canonicalize_context: CanonicalizeContext,
    #[case] atom_count: usize,
) {
    fn compare_refinement(
        partition: OrderedPartition,
        renumbered_partition: OrderedPartition,
        graph: &Graph,
        renumbered_graph: &Graph,
    ) {
        let reference = partition.clone().refine_by_all_cells(graph);
        let active = partition.clone().refine_by_active_cells(graph);
        let actual = partition.refine(graph);
        assert_eq!(active, reference);
        assert_eq!(actual, reference);
        let renumbered = renumbered_partition.refine(renumbered_graph);
        let node_count = graph.node_count();
        let transported = OrderedPartition {
            cells: renumbered
                .cells
                .iter()
                .map(|cell| {
                    let mut cell = cell
                        .iter()
                        .map(|node| NodeId((node_count - 1 - node.index()) as u32))
                        .collect::<Vec<_>>();
                    cell.sort_unstable();
                    cell
                })
                .collect(),
        };
        assert_eq!(actual, transported);

        let Some(cell_index) = actual.first_non_singleton_entity_cell(graph.node_count()) else {
            return;
        };
        for node in actual.cells[cell_index].clone() {
            let image = NodeId((node_count - 1 - node.index()) as u32);
            compare_refinement(
                actual.individualize(cell_index, node),
                renumbered.individualize(cell_index, image),
                graph,
                renumbered_graph,
            );
        }
    }

    let endpoint_pairs = (0..atom_count as u32)
        .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for complete_bond_form in [false, true] {
        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let selected_pairs = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &pair)| ((edge_mask >> position) & 1 == 1).then_some(pair))
                .collect::<Vec<_>>();
            let graph = Graph::new(atom_count, &selected_pairs);
            let renumbered_pairs = selected_pairs
                .iter()
                .map(|&[first, second]| {
                    [
                        atom_count as u32 - 1 - first,
                        atom_count as u32 - 1 - second,
                    ]
                })
                .collect::<Vec<_>>();
            let renumbered_graph = Graph::new(atom_count, &renumbered_pairs);
            compare_refinement(
                OrderedPartition::from_descriptors(&vec![0_u32; atom_count]),
                OrderedPartition::from_descriptors(&vec![0_u32; atom_count]),
                &graph,
                &renumbered_graph,
            );

            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds: selected_pairs
                    .iter()
                    .map(|&[first, second]| {
                        let bond = BondForm::from_order(1);
                        (
                            AtomId(first),
                            AtomId(second),
                            if complete_bond_form {
                                bond.into_concrete()
                            } else {
                                bond
                            },
                        )
                    })
                    .collect(),
                ..Default::default()
            });
            let renumbered = molecule.clone().remap(&reverse_remapping(&molecule));
            let (expected, expected_remapping) = canonicalize_topology(
                molecule.clone(),
                &canonicalize_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    branch_order: BranchOrdering::BackendCanonical,
                },
            )
            .unwrap();
            let (actual, actual_remapping) = canonicalize_topology(
                renumbered.clone(),
                &canonicalize_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    branch_order: BranchOrdering::BackendCanonical,
                },
            )
            .unwrap();

            assert_eq!(actual, expected);
            assert!(molecule.framed_eq_under(&expected, &expected_remapping));
            assert!(renumbered.framed_eq_under(&actual, &actual_remapping));
        }
    }
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
    let (entity_keys, incidence_keys) = initial_color_keys(&molecule, &incidence_graph).unwrap();
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Topology,
            order,
        )
    };
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
            branch_order: BranchOrdering::NodeId,
        },
        &leaf_candidate,
    );
    let reversed = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
        &leaf_candidate,
    );
    let pruned = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
        &leaf_candidate,
    );

    assert_eq!(unpruned.candidate.key, expected.key);
    assert_eq!(reversed.candidate.key, expected.key);
    assert_eq!(pruned.candidate.key, expected.key);
    assert!(pruned.stats.visited_leaves <= unpruned.stats.visited_leaves);
}

#[rstest]
fn test_canonical_search_color_labels() {
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
    let (entity_keys, incidence_keys) = initial_color_keys(&molecule, &incidence_graph).unwrap();
    let colors = rank_initial_colors(&entity_keys, &incidence_keys);
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let mut relabeled = adapter.clone();
    relabeled.colors.iter_mut().for_each(|color| {
        *color = match *color {
            AutomorphismColor::Entity(value) => AutomorphismColor::Entity(u32::MAX - value),
            AutomorphismColor::Incidence(value) => AutomorphismColor::Incidence(u32::MAX - value),
        }
    });
    let descriptors = partition_descriptors(&adapter, &entity_keys, &incidence_keys);
    let leaf_candidate = |order: &[NodeId]| {
        leaf_candidate_at(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Topology,
            order,
        )
    };
    let options = CanonicalSearchOptions {
        automorphism_pruning: true,
        branch_order: BranchOrdering::BackendCanonical,
    };

    let expected = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        options,
        &leaf_candidate,
    );
    let actual = canonical_search(
        &relabeled,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        options,
        &leaf_candidate,
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
        molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(expected),
    );
}

#[rstest]
#[case::mixed_bond_orders(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 4],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (AtomId(1), AtomId(2), BondForm::from_order(2).into_concrete()),
        (AtomId(2), AtomId(3), BondForm::from_order(1).into_concrete()),
        (AtomId(3), AtomId(0), BondForm::from_order(3).into_concrete()),
    ],
    ..Default::default()
}))]
#[case::electronic_bond_forms(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 4],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (
            AtomId(1),
            AtomId(2),
            BondForm::from_order(1).with_charge(1_i64).into_concrete(),
        ),
        (
            AtomId(2),
            AtomId(3),
            BondForm::from_order(1)
                .with_charge(0_i64)
                .with_unpaired_electrons((2_u8, 3_u8)),
        ),
    ],
    ..Default::default()
}))]
#[case::disconnected(Molecule::from_entries(MoleculeEntries {
    atoms: vec![
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::O),
    ],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (AtomId(2), AtomId(3), BondForm::from_order(2).into_concrete()),
    ],
    ..Default::default()
}))]
#[case::topology_path(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 32],
    bonds: (0_u32..31)
        .map(|first| (
            AtomId(first),
            AtomId(first + 1),
            BondForm::from_order(1),
        ))
        .collect(),
    ..Default::default()
}))]
fn test_canonicalize_topology_relabeling(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
) {
    let options = CanonicalSearchOptions {
        automorphism_pruning: true,
        branch_order: BranchOrdering::BackendCanonical,
    };
    let expected = molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .unwrap();
    let renumbered = molecule.clone().remap(&reverse_remapping(&molecule));

    for source in [&molecule, &renumbered, &expected] {
        let (actual, remapping) =
            canonicalize_topology(source.clone(), &canonicalize_context, options).unwrap();

        assert_eq!(actual, expected);
        assert!(source.framed_eq_under(&actual, &remapping));
    }
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
        molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
        molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
    let remapping = molecule_remapping(&[
        vec![1, 0],
        Vec::new(),
        vec![0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ]);
    let remapped = molecule.clone().remap(&remapping);

    let (canonical, remapping) = canonicalize_topology(
        molecule.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let (canonical_remapped, remapped_remapping) = canonicalize_topology(
        remapped.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Topology);
    let canonical_remapped_incidence = canonical_remapped.incidence_graph(IncidenceLevel::Topology);
    let canonical_again = canonical
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .unwrap();
    let canonical_again_incidence = canonical_again.incidence_graph(IncidenceLevel::Topology);

    assert!(molecule.framed_eq_under(&canonical, &remapping));
    assert!(remapped.framed_eq_under(&canonical_remapped, &remapped_remapping));
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical_remapped.check_integrity(), Ok(()));
    assert_eq!(
        leaf_candidate_at(
            &canonical,
            &canonical_incidence,
            DescriptionLevel::Topology,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
        leaf_candidate_at(
            &canonical_again,
            &canonical_again_incidence,
            DescriptionLevel::Topology,
            &canonical_again_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>()
        )
        .key,
    );
    assert_eq!(
        leaf_candidate_at(
            &canonical,
            &canonical_incidence,
            DescriptionLevel::Topology,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
        leaf_candidate_at(
            &canonical_remapped,
            &canonical_remapped_incidence,
            DescriptionLevel::Topology,
            &canonical_remapped_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>()
        )
        .key,
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

    for complete_bond_form in [false, true] {
        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let bonds = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &[first, second])| {
                    if (edge_mask >> position) & 1 == 0 {
                        return None;
                    }
                    let bond = BondForm::from_order(1);
                    Some((
                        AtomId(first),
                        AtomId(second),
                        if complete_bond_form {
                            bond.into_concrete()
                        } else {
                            bond
                        },
                    ))
                })
                .collect::<Vec<_>>();
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds,
                ..Default::default()
            });

            let (canonical, remapping) = canonicalize_topology(
                molecule.clone(),
                &canonicalize_context,
                CanonicalSearchOptions {
                    automorphism_pruning: true,
                    branch_order: BranchOrdering::BackendCanonical,
                },
            )
            .unwrap();
            let (unpruned, _) = canonicalize_topology(
                molecule.clone(),
                &canonicalize_context,
                CanonicalSearchOptions {
                    automorphism_pruning: false,
                    branch_order: BranchOrdering::ReverseNodeId,
                },
            )
            .unwrap();

            assert_eq!(
                unpruned, canonical,
                "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}",
            );
            assert!(
                molecule.framed_eq_under(&canonical, &remapping),
                "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}",
            );
            assert_eq!(canonical.check_integrity(), Ok(()));
            assert_eq!(
                canonical
                    .clone()
                    .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
                    .map(|(canonical, _)| canonical),
                Ok(canonical.clone()),
                "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}",
            );

            for (index, atom_images) in permutations(atom_count).into_iter().enumerate() {
                let bond_count = molecule.bonds().count();
                let bond_images = if index % 2 == 0 {
                    (0..bond_count).collect::<Vec<_>>()
                } else {
                    (0..bond_count).rev().collect::<Vec<_>>()
                };
                let renumbering = molecule_remapping(&[
                    atom_images,
                    bond_images,
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                ]);
                let renumbered = molecule.clone().remap(&renumbering);

                assert_eq!(
                    renumbered.clone().tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context).map(|(canonical, _)| canonical),
                    Ok(canonical.clone()),
                    "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}, renumbering {index}",
                );
            }
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
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
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
                [AtomId(1), AtomId(2)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
            (
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    });
    let expected_remapping = molecule_remapping(&[
        vec![3, 1, 2, 0, 4, 5],
        vec![1, 0, 3, 2],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![1, 0],
        vec![0, 1],
        vec![0, 1],
    ]);
    let expected = molecule
        .clone()
        .remap(&expected_remapping)
        .reframe()
        .unwrap();

    assert_eq!(
        canonicalize_constitution(
            molecule.clone(),
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                branch_order: BranchOrdering::BackendCanonical,
            },
        ),
        Ok((expected, expected_remapping)),
    );
}

#[rstest]
fn test_canonicalize_constitution_excluded_data(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(3))),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
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
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    });

    let (_, left_remapping) = canonicalize_constitution(
        left.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let (_, right_remapping) = canonicalize_constitution(
        right.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();

    assert_eq!(right_remapping, left_remapping);
}

#[rstest]
fn test_canonicalize_constitution_properties(
    initial_color_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let normalized_source = initial_color_molecule.clone().normalize().unwrap();
    let (canonical, remapping) = canonicalize_constitution(
        initial_color_molecule.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let acted = initial_color_molecule
        .clone()
        .remap(&remapping)
        .reframe()
        .unwrap();
    fn inverse_images<Id: Copy + Into<usize> + From<usize>>(
        remapping: &Remapping<Id>,
    ) -> Vec<usize> {
        let mut sources = (0..remapping.len()).collect::<Vec<_>>();
        sources
            .sort_unstable_by_key(|&source| Into::<usize>::into(remapping.map(Id::from(source))));
        sources
    }
    let inverse = molecule_remapping(&[
        inverse_images(remapping.graph().nodes()),
        inverse_images(remapping.graph().edges()),
        inverse_images(remapping.dative_bonds()),
        inverse_images(remapping.aromatic_systems()),
        inverse_images(remapping.multicenter_bonds()),
        inverse_images(remapping.noncovalent_bonds()),
        inverse_images(remapping.stereo_atoms()),
        inverse_images(remapping.stereo_bonds()),
    ]);

    assert_eq!(acted, canonical);
    assert!(initial_color_molecule.framed_eq_under(&canonical, &remapping));
    assert!(canonical.framed_eq_under(&normalized_source, &inverse));
    assert_eq!(canonical.check_integrity(), Ok(()));

    let (canonical_again, _) = canonicalize_constitution(
        canonical.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
    let canonical_again_incidence = canonical_again.incidence_graph(IncidenceLevel::Constitution);
    assert_eq!(
        leaf_candidate_at(
            &canonical,
            &canonical_incidence,
            DescriptionLevel::Constitution,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
        leaf_candidate_at(
            &canonical_again,
            &canonical_again_incidence,
            DescriptionLevel::Constitution,
            &canonical_again_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>()
        )
        .key,
    );

    let renumbering = reverse_remapping(&initial_color_molecule);
    let renumbered = initial_color_molecule.clone().remap(&renumbering);
    let (canonical_renumbered, renumbered_remapping) = canonicalize_constitution(
        renumbered.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .unwrap();
    let composed = molecule_remapping(&[
        (0..renumbering.graph().nodes().len())
            .map(|idx| {
                renumbered_remapping
                    .graph()
                    .nodes()
                    .map(renumbering.graph().nodes().map(NodeId::from(idx)))
                    .index()
            })
            .collect(),
        (0..renumbering.graph().edges().len())
            .map(|idx| {
                renumbered_remapping
                    .graph()
                    .edges()
                    .map(renumbering.graph().edges().map(EdgeId::from(idx)))
                    .index()
            })
            .collect(),
        (0..renumbering.dative_bonds().len())
            .map(|idx| {
                renumbered_remapping
                    .dative_bonds()
                    .map(renumbering.dative_bonds().map(DativeBondId::from(idx)))
                    .index()
            })
            .collect(),
        (0..renumbering.aromatic_systems().len())
            .map(|idx| {
                renumbered_remapping
                    .aromatic_systems()
                    .map(
                        renumbering
                            .aromatic_systems()
                            .map(AromaticSystemId::from(idx)),
                    )
                    .index()
            })
            .collect(),
        (0..renumbering.multicenter_bonds().len())
            .map(|idx| {
                renumbered_remapping
                    .multicenter_bonds()
                    .map(
                        renumbering
                            .multicenter_bonds()
                            .map(MulticenterBondId::from(idx)),
                    )
                    .index()
            })
            .collect(),
        (0..renumbering.noncovalent_bonds().len())
            .map(|idx| {
                renumbered_remapping
                    .noncovalent_bonds()
                    .map(
                        renumbering
                            .noncovalent_bonds()
                            .map(NoncovalentBondId::from(idx)),
                    )
                    .index()
            })
            .collect(),
        (0..renumbering.stereo_atoms().len())
            .map(|idx| {
                renumbered_remapping
                    .stereo_atoms()
                    .map(renumbering.stereo_atoms().map(StereoAtomId::from(idx)))
                    .index()
            })
            .collect(),
        (0..renumbering.stereo_bonds().len())
            .map(|idx| {
                renumbered_remapping
                    .stereo_bonds()
                    .map(renumbering.stereo_bonds().map(StereoBondId::from(idx)))
                    .index()
            })
            .collect(),
    ]);
    let composed_action = initial_color_molecule
        .clone()
        .remap(&composed)
        .reframe()
        .unwrap();
    let canonical_renumbered_incidence =
        canonical_renumbered.incidence_graph(IncidenceLevel::Constitution);

    assert_eq!(composed_action, canonical_renumbered);
    assert!(initial_color_molecule.framed_eq_under(&canonical_renumbered, &composed));
    assert_eq!(
        leaf_candidate_at(
            &canonical,
            &canonical_incidence,
            DescriptionLevel::Constitution,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
        leaf_candidate_at(
            &canonical_renumbered,
            &canonical_renumbered_incidence,
            DescriptionLevel::Constitution,
            &canonical_renumbered_incidence
                .graph()
                .node_ids()
                .collect::<Vec<_>>()
        )
        .key,
    );

    let (unpruned, _) = canonicalize_constitution(
        initial_color_molecule.clone(),
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::ReverseNodeId,
        },
    )
    .unwrap();
    let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
    assert_eq!(
        leaf_candidate_at(
            &canonical,
            &canonical_incidence,
            DescriptionLevel::Constitution,
            &canonical_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
        leaf_candidate_at(
            &unpruned,
            &unpruned_incidence,
            DescriptionLevel::Constitution,
            &unpruned_incidence.graph().node_ids().collect::<Vec<_>>()
        )
        .key,
    );
}

#[rstest]
fn test_canonicalize_constitution_entity_kind_minimum(canonicalize_context: CanonicalizeContext) {
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
                        [AtomId(0), AtomId(1)],
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                    (
                        [AtomId(2), AtomId(3)],
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                ],
                ..Default::default()
            }),
        ),
    ];

    for (entity_kind, molecule) in cases {
        let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
        let (entity_keys, incidence_keys) =
            initial_color_keys(&molecule, &incidence_graph).unwrap();
        let colors = rank_initial_colors(&entity_keys, &incidence_keys);
        let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
        let leaf_candidate = |order: &[NodeId]| {
            leaf_candidate_at(
                &molecule,
                &incidence_graph,
                DescriptionLevel::Constitution,
                order,
            )
        };
        let expected = exhaustive_minimum(
            &adapter,
            adapter_entity_blocks(&incidence_graph),
            &leaf_candidate,
        );
        let (canonical, remapping) = canonicalize_constitution(
            molecule.clone(),
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: true,
                branch_order: BranchOrdering::BackendCanonical,
            },
        )
        .unwrap();
        let (unpruned, _) = canonicalize_constitution(
            molecule.clone(),
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning: false,
                branch_order: BranchOrdering::ReverseNodeId,
            },
        )
        .unwrap();
        let canonical_incidence = canonical.incidence_graph(IncidenceLevel::Constitution);
        let canonical_order = canonical_incidence.graph().node_ids().collect::<Vec<_>>();
        let unpruned_incidence = unpruned.incidence_graph(IncidenceLevel::Constitution);
        let unpruned_order = unpruned_incidence.graph().node_ids().collect::<Vec<_>>();

        assert_eq!(
            leaf_candidate_at(
                &unpruned,
                &unpruned_incidence,
                DescriptionLevel::Constitution,
                &unpruned_order
            )
            .key,
            expected.key.clone(),
            "unpruned {entity_kind}",
        );
        assert_eq!(
            leaf_candidate_at(
                &canonical,
                &canonical_incidence,
                DescriptionLevel::Constitution,
                &canonical_order
            )
            .key,
            expected.key,
            "pruned {entity_kind}",
        );
        assert_eq!(unpruned, canonical, "{entity_kind}");
        assert!(
            molecule.framed_eq_under(&canonical, &remapping),
            "{entity_kind}"
        );
        assert_eq!(canonical.check_integrity(), Ok(()), "{entity_kind}");

        for (index, atom_images) in permutations(molecule.atoms().count())
            .into_iter()
            .enumerate()
        {
            let mut images = molecule_counts(&molecule).map(|count| (0..count).collect::<Vec<_>>());
            images[0] = atom_images;
            if index % 2 == 1 {
                for entity_kind_images in &mut images[1..6] {
                    entity_kind_images.reverse();
                }
            }
            let renumbered = molecule.clone().remap(&molecule_remapping(&images));

            assert_eq!(
                renumbered
                    .clone()
                    .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
                    .map(|(canonical, _)| canonical),
                Ok(canonical.clone()),
                "{entity_kind}, renumbering {index}",
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
        right
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        left.clone()
            .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
            .map(|(canonical, _)| canonical),
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
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
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
            molecule
                .clone()
                .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
                .map(|(canonical, _)| canonical),
            Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
            "{location}",
        );
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
    let remapped = molecule.clone().remap(&reverse_remapping(&molecule));

    assert!(colored_encoding_equivalent(&complete, &molecule, level));
    assert_eq!(
        colored_encoding_equivalent(&molecule, &remapped, level),
        explicitly_dense_equivalent(&molecule, &remapped),
    );
    assert!(explicitly_dense_equivalent(&molecule, &remapped));

    let mut distinguished = remapped;
    distinguished.atom_mut(AtomId(0)).attributes_mut().element = ElementForm::Lit(Element::O);
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
        let remapped = molecule.clone().remap(&reverse_remapping(&molecule));
        assert_eq!(
            colored_encoding_equivalent(&molecule, &remapped, IncidenceLevel::Topology),
            explicitly_dense_equivalent(&molecule, &remapped),
            "edge mask {edge_mask:#08b}",
        );

        let mut distinguished = remapped;
        distinguished.atom_mut(AtomId(0)).attributes_mut().element = ElementForm::Lit(Element::O);
        assert_eq!(
            colored_encoding_equivalent(&molecule, &distinguished, IncidenceLevel::Topology,),
            explicitly_dense_equivalent(&molecule, &distinguished),
            "edge mask {edge_mask:#08b}",
        );
    }
}

#[rstest]
#[case::charge(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C).with_charge(NumForm::LitSet(Box::default()))],
        ..Default::default()
    }),
    IncidenceLevel::Topology,
)]
#[case::stereo_atom_configuration(
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(AtomForm::from_element)
            .into(),
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
        )],
        ..Default::default()
    }),
    IncidenceLevel::Full,
)]
#[case::stereo_bond_configuration(
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(AtomForm::from_element)
            .into(),
        bonds: [(0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1)]
            .map(|(first, second, order)| (AtomId(first), AtomId(second), BondForm::from_order(order)))
            .into(),
        stereo_bonds: vec![(
            BondId(0),
            [3, 2, 5, 4]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([])),
        )],
        ..Default::default()
    }),
    IncidenceLevel::Full,
)]
fn test_initial_colors_error(#[case] molecule: Molecule, #[case] incidence_level: IncidenceLevel) {
    let incidence_graph = molecule.incidence_graph(incidence_level);

    assert_eq!(
        initial_colors(&molecule, &incidence_graph),
        Err(Contradiction)
    );
}
