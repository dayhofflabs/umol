use std::collections::HashSet;
use std::iter;
use std::ops::ControlFlow;
use std::sync::Arc;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, BipartiteMaximumMatchingAlgorithm,
    Compaction, ConnectedComponentsAlgorithm, Correspondence, EdgeId,
    GeneralMaximumMatchingAlgorithm, Graph, GraphCompaction, GraphRemapping,
    MatchingEnumerationAlgorithm, MaximumIndependentSetAlgorithm, NodeId, NonBipartiteGraphError,
    RelevantCycleEnumerationAlgorithm, Remapping, ShortestCycleAlgorithm,
    SimpleCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_perm::{DynPermutation, Permutation, MAX_DEGREE};

use super::super::aromatic::AromaticSystemForm;
use super::super::atom::{AtomForm, ElementForm, IsotopeMassForm};
use super::super::bond::BondForm;
use super::super::boolean::BooleanForm;
use super::super::compact::MoleculeCompaction;
use super::super::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, AtomConstraintsForm, BondConstraintForm,
    BondConstraintsForm, Constraint, Constraints, DativeBondConstraintForm,
    DativeBondConstraintsForm, FluxionalityForm, LigandPermutation, MoleculeConstraint,
    MulticenterBondConstraintForm, NoncovalentBondConstraintForm, RelationalConstraint, RingScope,
    StereoAtomConstraintForm, StereoBondConstraintForm, StereoLigandPair, StereogenicityForm,
    TopicityForm, TopicityRelationForm,
};
use super::super::correspondence::MoleculeCorrespondence;
use super::super::dative::DativeBondForm;
use super::super::edit::{
    AddBond, AromaticSystemHandle, AtomFieldChange, AtomHandle, BondHandle, DativeBondHandle, Edit,
    Edits, MulticenterBondHandle, NoncovalentBondHandle, StereoAtomHandle, StereoBondHandle,
};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::{Entity, EntityKind};
use super::super::error::Contradiction;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::multicenter::MulticenterBondForm;
use super::super::noncovalent::{
    NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
};
use super::super::num::NumForm;
use super::super::ring::{RingConfig, RingModel, RingSetKind};
use super::super::spin::UnpairedElectronsForm;
use super::super::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind, Topicity,
};
use super::super::traits::{FrameTransport, Normalize, Reframe};
use super::transact::TransactionError;
use super::{
    AromaticSystems, DativeBonds, Molecule, MoleculeApplyError, MoleculeEntries,
    MoleculeIntegrityError, MulticenterBonds, NoncovalentBonds, StereoAtoms, StereoBonds,
};
use crate::ir::MoleculeRemapping;
use crate::{mol_dsl, mol_dsl_concrete};

fn ground_atom() -> AtomForm {
    let mut a = AtomForm::from_element(Element::C);
    a.isotope_mass = IsotopeMassForm::Natural;
    a.charge = NumForm::Lit(0);
    a.implicit_hydrogens = NumForm::Lit(4);
    a.lone_pairs = NumForm::Lit(0);
    a.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
    a
}

fn constraints_with_molecule(c: Constraint) -> Constraints {
    let mut out = Constraints::new();
    out.push(c);
    out
}

#[rstest]
fn test_molecule_new() {
    let m = Molecule::new();
    assert_eq!(m.atoms().count(), 0);
    assert_eq!(m.bonds().count(), 0);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
    assert_eq!(m.constraints().len(), 0);
}

#[rstest]
fn test_molecule_default_equals_new() {
    assert_eq!(Molecule::default(), Molecule::new());
}

#[rstest]
fn test_molecule_from_entries() {
    let atoms = vec![
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::O),
    ];
    let bonds = vec![(AtomId(0), AtomId(1), BondForm::from_order(1))];
    let m = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    });
    assert_eq!(m.atoms().count(), 2);
    assert_eq!(m.bonds().count(), 1);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
    assert_eq!(
        m.atom(AtomId(0)).attributes.element,
        ElementForm::Lit(Element::C)
    );
    assert_eq!(
        m.atom(AtomId(1)).attributes.element,
        ElementForm::Lit(Element::O)
    );
    assert_eq!(m.bond(BondId(0)).attributes.order, NumForm::Lit(1));
}

#[rstest]
fn test_molecule_builder() {
    assert_eq!(Molecule::builder().build(), Molecule::new());
}

#[rstest]
#[case::empty(Molecule::default(), true)]
#[case::ground_atom(
    mol_dsl_concrete!(r#"{:atoms ["C #h4"] :bonds []}"#),
    true,
)]
#[case::wildcard_element(
    mol_dsl!(r#"{:atoms ["*"] :bonds []}"#),
    false,
)]
#[case::wildcard_bond(
    mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "*"]]}"#),
    false,
)]
#[case::ground_atom_with_undetermined_constraint(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![ground_atom()],
        constraints: constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![]),
            sum: NumForm::Undetermined,
        })),
        ..Default::default()
    }),
    true,
)]
#[case::stereo_atom_ground_coset(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![ground_atom(); 5],
        bonds: (1..=4)
            .map(|id| {
                (
                    AtomId(0),
                    AtomId(id),
                    BondForm::from_order(1).into_concrete(),
                )
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::new(),
        ..Default::default()
    }),
    true,
)]
#[case::stereo_atom_undetermined_coset(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![ground_atom(); 5],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        )],
        constraints: Constraints::new(),
        ..Default::default()
    }),
    false,
)]
fn test_molecule_is_concrete(#[case] molecule: Molecule, #[case] expected: bool) {
    assert_eq!(molecule.is_concrete(), expected);
}

#[rstest]
#[case::hub(AtomId(0), vec![(AtomId(1), BondId(0)), (AtomId(2), BondId(1))])]
#[case::leaf_o(AtomId(1), vec![(AtomId(0), BondId(0))])]
#[case::leaf_n(AtomId(2), vec![(AtomId(0), BondId(1))])]
#[case::isolated(AtomId(3), vec![])]
fn test_molecule_neighbors(#[case] atom: AtomId, #[case] expected: Vec<(AtomId, BondId)>) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(2)),
        ],
        ..Default::default()
    });
    let mut neighbors = molecule.neighbors(atom);
    assert_eq!(neighbors.len(), expected.len());
    assert_eq!(
        neighbors.size_hint(),
        (expected.len(), Some(expected.len())),
    );
    while let Some(expected_neighbor) = expected.get(expected.len() - neighbors.len()) {
        let previous = neighbors.len();
        assert_eq!(
            neighbors.next().map(|n| (n.atom_id(), n.bond_id())),
            Some(*expected_neighbor),
        );
        let remaining = neighbors.len();
        assert_eq!(remaining, previous - 1);
        assert_eq!(neighbors.size_hint(), (remaining, Some(remaining)));
    }
    assert_eq!(neighbors.next().map(|n| (n.atom_id(), n.bond_id())), None,);
    assert_eq!(neighbors.len(), 0);
}

#[rstest]
fn test_molecule_editor_add_aromatic_system() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let mut b = molecule.edit();
    let id = b.add_aromatic_system(vec![AtomId(0), AtomId(1)], AromaticSystemForm::default());
    let new_molecule = b.build();
    assert_eq!(id, AromaticSystemId(0));
    let new_atoms: Vec<AtomId> = new_molecule
        .aromatic_system(AromaticSystemId(0))
        .atom_ids()
        .collect();
    assert_eq!(new_atoms, vec![AtomId(0), AtomId(1)]);
    assert_eq!(
        new_molecule.aromatic_systems().ids().collect::<Vec<_>>(),
        vec![AromaticSystemId(0)]
    );
    assert_eq!(
        molecule.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
}

#[fixture]
fn rich_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondForm::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            [AtomId(0), AtomId(3)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        ..Default::default()
    })
}

#[fixture]
fn equiv_molecule_entries() -> MoleculeEntries {
    let mut carbon = AtomForm::from_element(Element::C);
    carbon.charge = NumForm::Lit(1);

    MoleculeEntries {
        atoms: vec![
            carbon,
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![(
            vec![AtomId(1), AtomId(2)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![1, 2, 0]),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![2, 1, 0]),
        )],
        noncovalent: vec![(
            [AtomId(0), AtomId(3)],
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
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(2)]),
            },
        )),
    }
}

#[fixture]
fn molecule_reframe_source(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
) -> Molecule {
    entries.dative[0].0 = vec![AtomId(2), AtomId(1)];
    entries.aromatic[0] = (
        vec![AtomId(2), AtomId(0), AtomId(1)],
        AromaticSystemForm::from_electrons(vec![20, 10, 15]),
    );
    entries.multicenter[0] = (
        vec![AtomId(1), AtomId(2), AtomId(0)],
        MulticenterBondForm::from_electrons(vec![12, 14, 10]),
    );
    entries.noncovalent[0].0 = [AtomId(3), AtomId(0)];

    entries.stereo_atoms[0].1 = vec![
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ];
    entries.stereo_atoms[0]
        .2
        .constraints
        .set(StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        }));

    entries.stereo_bonds[0].1 = vec![
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    ];
    entries.stereo_bonds[0]
        .2
        .constraints
        .set(StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        }));

    entries.constraints = Constraints::new();
    entries.constraints.push(Constraint::Relational(
        RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::valence(1)),
                Box::new(AtomConstraintForm::valence(2)),
            ],
        },
    ));
    entries.constraints.push(Constraint::StereoAtom(
        StereoAtomId(0),
        StereoKind::Tetrahedral,
        StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }),
    ));
    entries.constraints.push(Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
        }),
    ));

    Molecule::from_entries(entries)
}

#[rstest]
fn test_molecule_try_from_entries(#[from(equiv_molecule_entries)] mut entries: MoleculeEntries) {
    entries.dative.push((
        vec![AtomId(0), AtomId(1)],
        AtomId(3),
        DativeBondForm::from_order(2),
    ));
    let expected = Molecule {
        graph: Graph::new(
            entries.atoms.len(),
            &entries
                .bonds
                .iter()
                .map(|(first, second, _)| [first.0, second.0])
                .collect::<Vec<_>>(),
        ),
        atoms: Arc::new(entries.atoms.clone()),
        bonds: Arc::new(
            entries
                .bonds
                .iter()
                .map(|(_, _, attributes)| attributes.clone())
                .collect(),
        ),
        dative_bonds: DativeBonds::new(entries.dative.clone()),
        aromatic_systems: AromaticSystems::new(entries.aromatic.clone()),
        multicenter_bonds: MulticenterBonds::new(entries.multicenter.clone()),
        noncovalent_bonds: NoncovalentBonds::new(entries.noncovalent.clone()),
        stereo_atoms: StereoAtoms::new(entries.stereo_atoms.clone()),
        stereo_bonds: StereoBonds::new(entries.stereo_bonds.clone()),
        constraints: entries.constraints.clone(),
    };

    assert_eq!(Molecule::try_from_entries(entries), Ok(expected));
}

#[rstest]
#[case::bond_endpoint(
    |entries: &mut MoleculeEntries| entries.bonds[0].0 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::dative_donor(
    |entries: &mut MoleculeEntries| entries.dative[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::dative_acceptor(
    |entries: &mut MoleculeEntries| entries.dative[0].1 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::aromatic_participant(
    |entries: &mut MoleculeEntries| entries.aromatic[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::multicenter_participant(
    |entries: &mut MoleculeEntries| entries.multicenter[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::noncovalent_endpoint(
    |entries: &mut MoleculeEntries| entries.noncovalent[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_atom_site(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].0 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_atom_ligand(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[0].atom_id = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_bond_site(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].0 = BondId(3),
    Entity::Bond(BondId(3)),
)]
#[case::stereo_bond_ligand(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].1[0].atom_id = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
fn test_molecule_try_from_entries_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] entity: Entity,
) {
    invalidate(&mut entries);

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::InvalidReference { entity }),
    );
}

#[rstest]
#[case::aromatic_electron_count(
    |entries: &mut MoleculeEntries| {
        entries.aromatic[0].1.electrons = ElectronCountsForm::Lit(vec![2]);
    },
    MoleculeIntegrityError::ElectronCountLengthMismatch {
        entity: Entity::AromaticSystem(AromaticSystemId(0)),
        participants: 3,
        electron_counts: 1,
    },
)]
#[case::multicenter_electron_count(
    |entries: &mut MoleculeEntries| {
        entries.multicenter[0].1.electrons = ElectronCountsForm::Lit(vec![2]);
    },
    MoleculeIntegrityError::ElectronCountLengthMismatch {
        entity: Entity::MulticenterBond(MulticenterBondId(0)),
        participants: 3,
        electron_counts: 1,
    },
)]
#[case::aromatic_duplicate_participant(
    |entries: &mut MoleculeEntries| entries.aromatic[0].0[2] = AtomId(1),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::AromaticSystem(AromaticSystemId(0)),
        atom: AtomId(1),
    },
)]
#[case::multicenter_duplicate_participant(
    |entries: &mut MoleculeEntries| entries.multicenter[0].0[2] = AtomId(1),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::MulticenterBond(MulticenterBondId(0)),
        atom: AtomId(1),
    },
)]
#[case::bond_self_loop(
    |entries: &mut MoleculeEntries| entries.bonds[0].1 = AtomId(0),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::Bond(BondId(0)),
        atom: AtomId(0),
    },
)]
#[case::bonds_parallel(
    |entries: &mut MoleculeEntries| {
        entries.bonds[1].0 = AtomId(0);
        entries.bonds[1].1 = AtomId(1);
    },
    MoleculeIntegrityError::BondsParallel {
        atoms: [AtomId(0), AtomId(1)],
    },
)]
#[case::dative_donor_duplicate(
    |entries: &mut MoleculeEntries| entries.dative[0].0[1] = AtomId(1),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::DativeBond(DativeBondId(0)),
        atom: AtomId(1),
    },
)]
#[case::dative_donor_acceptor_duplicate(
    |entries: &mut MoleculeEntries| entries.dative[0].0[1] = AtomId(3),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::DativeBond(DativeBondId(0)),
        atom: AtomId(3),
    },
)]
#[case::dative_bonds_identical(
    |entries: &mut MoleculeEntries| entries.dative.push((
        vec![AtomId(2), AtomId(1)],
        AtomId(3),
        DativeBondForm::from_order(2),
    )),
    MoleculeIntegrityError::DativeBondsIdentical {
        acceptor: AtomId(3),
        donors: vec![AtomId(1), AtomId(2)],
    },
)]
#[case::aromatic_systems_overlap(
    |entries: &mut MoleculeEntries| entries.aromatic.push((
        vec![AtomId(2), AtomId(3)],
        AromaticSystemForm::default(),
    )),
    MoleculeIntegrityError::AromaticSystemsOverlap { atom: AtomId(2) },
)]
#[case::multicenter_bonds_identical(
    |entries: &mut MoleculeEntries| entries.multicenter.push((
        vec![AtomId(0), AtomId(1), AtomId(2)],
        MulticenterBondForm::default(),
    )),
    MoleculeIntegrityError::MulticenterBondsIdentical {
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
    },
)]
#[case::noncovalent_self_loop(
    |entries: &mut MoleculeEntries| entries.noncovalent[0].0[1] = AtomId(0),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::NoncovalentBond(NoncovalentBondId(0)),
        atom: AtomId(0),
    },
)]
#[case::noncovalent_bonds_parallel_distinct_kinds(
    |entries: &mut MoleculeEntries| entries.noncovalent.push((
        [AtomId(3), AtomId(0)],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::VanDerWaals),
    )),
    MoleculeIntegrityError::NoncovalentBondsParallel {
        atoms: [AtomId(0), AtomId(3)],
    },
)]
#[case::stereo_atom_site_is_atom_ligand(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[0] =
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        atom: AtomId(1),
    },
)]
#[case::stereo_atom_ligand_duplicate(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[1] =
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    MoleculeIntegrityError::DuplicateStereoLigand {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        ligand: StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    },
)]
#[case::stereo_atom_sites_duplicate(
    |entries: &mut MoleculeEntries| {
        let stereo_atom = entries.stereo_atoms[0].clone();
        entries.stereo_atoms.push(stereo_atom);
    },
    MoleculeIntegrityError::StereoAtomSitesDuplicate { atom: AtomId(1) },
)]
#[case::stereo_bond_ligand_duplicate(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].1[1] =
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    MoleculeIntegrityError::DuplicateStereoLigand {
        entity: Entity::StereoBond(StereoBondId(0)),
        ligand: StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    },
)]
#[case::stereo_bond_sites_duplicate(
    |entries: &mut MoleculeEntries| {
        let stereo_bond = entries.stereo_bonds[0].clone();
        entries.stereo_bonds.push(stereo_bond);
    },
    MoleculeIntegrityError::StereoBondSitesDuplicate { bond: BondId(1) },
)]
#[case::stereo_atom_actual_ligand_incidence(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[1] =
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    MoleculeIntegrityError::StereoLigandIncidenceMismatch {
        entity: Entity::StereoAtom(StereoAtomId(0)),
    },
)]
#[case::stereo_atom_virtual_ligand_incidence(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[2] =
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
    MoleculeIntegrityError::StereoLigandIncidenceMismatch {
        entity: Entity::StereoAtom(StereoAtomId(0)),
    },
)]
#[case::stereo_bond_ligand_block_incidence(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].1.swap(1, 3),
    MoleculeIntegrityError::StereoLigandIncidenceMismatch {
        entity: Entity::StereoBond(StereoBondId(0)),
    },
)]
#[case::molecule_stereo_atom_arity(
    |entries: &mut MoleculeEntries| {
        entries.constraints = Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::TrigonalBipyramidal,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        )
        .into();
    },
    MoleculeIntegrityError::StereoLigandArity {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        kind: StereoKind::TrigonalBipyramidal,
        expected: 5,
        actual: 4,
    },
)]
#[case::molecule_stereo_atom_permutation(
    |entries: &mut MoleculeEntries| {
        entries.constraints = Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Fluxionality(FluxionalityForm {
                permutation: LigandPermutation(Permutation::identity(3)),
                active: BooleanForm::Lit(true),
            }),
        )
        .into();
    },
    MoleculeIntegrityError::StereoPermutationDegree {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        expected: 4,
        actual: 3,
    },
)]
#[case::molecule_stereo_bond_position(
    |entries: &mut MoleculeEntries| {
        entries.constraints = Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(0usize.into(), 4usize.into()),
                relation: TopicityRelationForm::Undetermined,
            }),
        )
        .into();
    },
    MoleculeIntegrityError::StereoLigandPositionOutOfRange {
        entity: Entity::StereoBond(StereoBondId(0)),
        position: 4,
        degree: 4,
    },
)]
#[case::stereo_ligand_arity(
    |entries: &mut MoleculeEntries| {
        entries.stereo_atoms[0].1.pop();
    },
    MoleculeIntegrityError::StereoLigandArity {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        kind: StereoKind::Tetrahedral,
        expected: 4,
        actual: 3,
    },
)]
#[case::stereo_coset(
    |entries: &mut MoleculeEntries| {
        entries.stereo_atoms[0].2.configuration = StereoConfigurationForm::kinded(
            StereoKind::Tetrahedral,
            StereoCoset::Lit(2),
        );
    },
    MoleculeIntegrityError::StereoCosetOutOfRange {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        kind: StereoKind::Tetrahedral,
        coset: 2,
        count: 2,
    },
)]
fn test_molecule_try_from_entries_integrity_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] expected: MoleculeIntegrityError,
) {
    invalidate(&mut entries);

    assert_eq!(Molecule::try_from_entries(entries), Err(expected));
}

#[rstest]
#[case::actual_and_virtual_ligands(vec![
    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
])]
fn test_molecule_try_from_entries_stereo_atom(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ligands: Vec<StereoLigand>,
) {
    entries.stereo_atoms[0].1 = ligands.clone();
    let molecule = Molecule::try_from_entries(entries).expect("entries satisfy molecule integrity");

    assert_eq!(
        molecule.stereo_atom(StereoAtomId(0)).ligand_frame(),
        ligands
    );
}

#[rstest]
#[case::site_endpoint_order(vec![
    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
])]
#[case::exchanged_endpoint_blocks(vec![
    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
])]
#[case::same_virtual_kind_on_opposite_endpoints(vec![
    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
])]
fn test_molecule_try_from_entries_stereo_bond(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ligands: Vec<StereoLigand>,
) {
    entries.stereo_bonds[0].1 = ligands.clone();
    let molecule = Molecule::try_from_entries(entries).expect("entries satisfy molecule integrity");

    assert_eq!(
        molecule.stereo_bond(StereoBondId(0)).ligand_frame(),
        ligands
    );
}

#[rstest]
#[case::atom_implicit_hydrogen(
    |entries: &mut MoleculeEntries| {
        entries.stereo_atoms[0].1[3] = entries.stereo_atoms[0].1[2];
    },
    Entity::StereoAtom(StereoAtomId(0)),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
)]
#[case::atom_lone_pair(
    |entries: &mut MoleculeEntries| {
        entries.stereo_atoms[0].1[2] = entries.stereo_atoms[0].1[3];
    },
    Entity::StereoAtom(StereoAtomId(0)),
    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
)]
#[case::bond_implicit_hydrogen(
    |entries: &mut MoleculeEntries| {
        entries.stereo_bonds[0].1[0] = entries.stereo_bonds[0].1[1];
    },
    Entity::StereoBond(StereoBondId(0)),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
)]
#[case::bond_lone_pair(
    |entries: &mut MoleculeEntries| {
        let ligand = StereoLigand::new(AtomId(1), StereoLigandKind::LonePair);
        entries.stereo_bonds[0].1[0] = ligand;
        entries.stereo_bonds[0].1[1] = ligand;
    },
    Entity::StereoBond(StereoBondId(0)),
    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
)]
fn test_molecule_try_from_entries_duplicate_stereo_ligand(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] entity: Entity,
    #[case] ligand: StereoLigand,
) {
    invalidate(&mut entries);

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::DuplicateStereoLigand { entity, ligand }),
    );
}

#[rstest]
#[case::atom(|entries: &mut MoleculeEntries| {
    entries.stereo_atoms[0].1 = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
    ];
    entries.stereo_atoms[0].2.configuration = StereoConfigurationForm::default();
}, Entity::StereoAtom(StereoAtomId(0)))]
#[case::bond(|entries: &mut MoleculeEntries| {
    entries.stereo_bonds[0].1 = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    ];
    entries.stereo_bonds[0].2.configuration = StereoConfigurationForm::default();
}, Entity::StereoBond(StereoBondId(0)))]
fn test_molecule_try_from_entries_stereo_frame_degree_too_large(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] entity: Entity,
) {
    invalidate(&mut entries);

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::StereoFrameDegreeTooLarge {
            entity,
            degree: MAX_DEGREE + 1,
            maximum: MAX_DEGREE,
        }),
    );
}

#[rstest]
#[case::reference_before_maximum_and_duplicate(
    |entries: &mut MoleculeEntries| {
        let invalid = StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen);
        entries.stereo_atoms[0].1.resize(MAX_DEGREE + 1, invalid);
    },
    MoleculeIntegrityError::InvalidReference {
        entity: Entity::Atom(AtomId(4)),
    },
)]
#[case::maximum_before_duplicate_and_arity(
    |entries: &mut MoleculeEntries| {
        let duplicate = entries.stereo_atoms[0].1[0];
        entries.stereo_atoms[0].1.resize(MAX_DEGREE + 1, duplicate);
    },
    MoleculeIntegrityError::StereoFrameDegreeTooLarge {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        degree: MAX_DEGREE + 1,
        maximum: MAX_DEGREE,
    },
)]
#[case::duplicate_before_arity(
    |entries: &mut MoleculeEntries| {
        entries.stereo_atoms[0].1.pop();
        entries.stereo_atoms[0].1[1] = entries.stereo_atoms[0].1[0];
    },
    MoleculeIntegrityError::DuplicateStereoLigand {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        ligand: StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
    },
)]
fn test_molecule_try_from_entries_stereo_frame_error_precedence(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] expected: MoleculeIntegrityError,
) {
    invalidate(&mut entries);

    assert_eq!(Molecule::try_from_entries(entries), Err(expected));
}

#[rstest]
#[case::tetrahedral_atom(StereoKind::Tetrahedral, 4)]
#[case::square_planar_atom(StereoKind::SquarePlanar, 4)]
#[case::axial_atom(StereoKind::Axial, 4)]
#[case::trigonal_bipyramidal_atom(StereoKind::TrigonalBipyramidal, 5)]
#[case::octahedral_atom(StereoKind::Octahedral, 6)]
fn test_molecule_try_from_entries_stereo_atom_kind(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] kind: StereoKind,
    #[case] degree: usize,
) {
    entries.atoms.extend([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::C),
    ]);
    entries.bonds.extend([
        (AtomId(1), AtomId(4), BondForm::from_order(1)),
        (AtomId(1), AtomId(5), BondForm::from_order(1)),
    ]);
    let mut ligands = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ];
    ligands.truncate(degree);
    entries.stereo_atoms[0].1 = ligands;
    entries.stereo_atoms[0].2.configuration =
        StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0));

    let molecule = Molecule::try_from_entries(entries).expect("entries satisfy molecule integrity");

    assert_eq!(molecule.stereo_atom(StereoAtomId(0)).kind(), kind);
}

#[rstest]
#[case::cis_trans_bond(StereoKind::CisTrans)]
#[case::axial_bond(StereoKind::Axial)]
fn test_molecule_try_from_entries_stereo_bond_kind(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] kind: StereoKind,
) {
    entries.stereo_bonds[0].2.configuration =
        StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0));

    let molecule = Molecule::try_from_entries(entries).expect("entries satisfy molecule integrity");

    assert_eq!(molecule.stereo_bond(StereoBondId(0)).kind(), kind);
}

/// A stereo kind names a coordination geometry, and a geometry belongs to an atom or to a bond.
/// Arity cannot separate them: every kind below shares degree 4 with the fixture's frame, so the
/// arity check passes and only the site-kind rule rejects the pairing.
#[rstest]
#[case::cis_trans_on_atom(Entity::StereoAtom(StereoAtomId(0)), StereoKind::CisTrans)]
#[case::tetrahedral_on_bond(Entity::StereoBond(StereoBondId(0)), StereoKind::Tetrahedral)]
#[case::square_planar_on_bond(Entity::StereoBond(StereoBondId(0)), StereoKind::SquarePlanar)]
#[case::trigonal_bipyramidal_on_bond(
    Entity::StereoBond(StereoBondId(0)),
    StereoKind::TrigonalBipyramidal
)]
#[case::octahedral_on_bond(Entity::StereoBond(StereoBondId(0)), StereoKind::Octahedral)]
fn test_molecule_try_from_entries_stereo_kind_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] entity: Entity,
    #[case] kind: StereoKind,
) {
    match entity {
        Entity::StereoAtom(id) => {
            entries.stereo_atoms[id.index()].2.configuration =
                StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0));
        }
        Entity::StereoBond(id) => {
            entries.stereo_bonds[id.index()].2.configuration =
                StereoConfigurationForm::kinded(kind, StereoCoset::Lit(0));
        }
        _ => unreachable!("test cases contain only stereo entities"),
    }

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch { entity, kind }),
    );
}

#[rstest]
#[case::cis_trans_on_atom(
    Constraint::StereoAtom(
        StereoAtomId(0),
        StereoKind::CisTrans,
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ),
    Entity::StereoAtom(StereoAtomId(0)),
    StereoKind::CisTrans
)]
#[case::tetrahedral_on_bond(
    Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::Tetrahedral,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    ),
    Entity::StereoBond(StereoBondId(0)),
    StereoKind::Tetrahedral
)]
fn test_molecule_try_from_entries_stereo_wrapper_kind_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] constraint: Constraint,
    #[case] entity: Entity,
    #[case] kind: StereoKind,
) {
    entries.constraints = constraint.into();

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch { entity, kind }),
    );
}

#[rstest]
#[case::atom(Entity::Atom(AtomId(4)))]
#[case::bond(Entity::Bond(BondId(3)))]
#[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
#[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
#[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
#[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
fn test_molecule_try_from_entries_constraint_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] entity: Entity,
) {
    let constraint = match entity {
        Entity::Atom(id) => Constraint::Atom(id, AtomConstraintForm::valence(NumForm::Lit(4))),
        Entity::Bond(id) => Constraint::Bond(id, BondConstraintForm::aromatic(false)),
        Entity::DativeBond(id) => {
            Constraint::DativeBond(id, DativeBondConstraintForm::aromatic(false))
        }
        Entity::AromaticSystem(id) => Constraint::AromaticSystem(
            id,
            AromaticSystemConstraintForm::electron_count(NumForm::Lit(6)),
        ),
        Entity::MulticenterBond(id) => Constraint::MulticenterBond(
            id,
            MulticenterBondConstraintForm::electron_count(NumForm::Lit(2)),
        ),
        Entity::NoncovalentBond(id) => {
            Constraint::NoncovalentBond(id, NoncovalentBondConstraintForm::intramolecular(true))
        }
        Entity::StereoAtom(id) => Constraint::StereoAtom(
            id,
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Entity::StereoBond(id) => Constraint::StereoBond(
            id,
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
    };
    entries.constraints = Constraint::Not(Box::new(constraint)).into();

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::InvalidReference { entity }),
    );
}

#[rstest]
#[case::relational_atom(
    Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }),
    Entity::Atom(AtomId(4)),
)]
#[case::relational_bond(
    Constraint::Relational(RelationalConstraint::DativeBondParallels {
        dative: DativeBondId(0),
        parallel: BondId(3),
    }),
    Entity::Bond(BondId(3)),
)]
#[case::atom_subset(
    Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(4)]),
    }),
    Entity::Atom(AtomId(4)),
)]
#[case::bond_subset(
    Constraint::Molecule(MoleculeConstraint::BondOrderSum {
        bonds: Some(vec![BondId(3)]),
        sum: NumForm::Lit(1),
    }),
    Entity::Bond(BondId(3)),
)]
fn test_molecule_try_from_entries_molecule_constraint_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] constraint: Constraint,
    #[case] entity: Entity,
) {
    entries.constraints = constraint.into();

    assert_eq!(
        Molecule::try_from_entries(entries),
        Err(MoleculeIntegrityError::InvalidReference { entity }),
    );
}

#[rstest]
#[should_panic(expected = "invalid molecule entries: molecule references unavailable atom 1")]
fn test_molecule_from_entries_error() {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default()],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
        ..Default::default()
    });
}

#[fixture]
fn framed_eq_under_molecules(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) -> (Molecule, Molecule, MoleculeRemapping) {
    let atom_images = [AtomId(2), AtomId(3), AtomId(0), AtomId(1)];
    let map_atom = |id: AtomId| atom_images[id.index()];

    let mut right_atoms = vec![AtomForm::default(); entries.atoms.len()];
    for (index, atom) in entries.atoms.iter().cloned().enumerate() {
        right_atoms[map_atom(AtomId(index as u32)).index()] = atom;
    }
    let right_bonds = entries
        .bonds
        .iter()
        .cloned()
        .map(|(first, second, attributes)| (map_atom(first), map_atom(second), attributes))
        .collect();
    let right_dative = entries
        .dative
        .iter()
        .cloned()
        .map(|(donors, acceptor, attributes)| {
            (
                donors.into_iter().map(map_atom).collect(),
                map_atom(acceptor),
                attributes,
            )
        })
        .collect();
    let right_aromatic = entries
        .aromatic
        .iter()
        .cloned()
        .map(|(atoms, attributes)| (atoms.into_iter().map(map_atom).collect(), attributes))
        .collect();
    let right_multicenter = entries
        .multicenter
        .iter()
        .cloned()
        .map(|(atoms, attributes)| (atoms.into_iter().map(map_atom).collect(), attributes))
        .collect();
    let right_noncovalent = entries
        .noncovalent
        .iter()
        .cloned()
        .map(|(atoms, attributes)| (atoms.map(map_atom), attributes))
        .collect();
    let right_stereo_atoms = entries
        .stereo_atoms
        .iter()
        .cloned()
        .map(|(site, ligands, attributes)| {
            (
                map_atom(site),
                ligands
                    .into_iter()
                    .map(|ligand| StereoLigand::new(map_atom(ligand.atom_id), ligand.kind))
                    .collect(),
                attributes,
            )
        })
        .collect();
    let right_stereo_bonds = entries
        .stereo_bonds
        .iter()
        .cloned()
        .map(|(site, ligands, attributes)| {
            (
                site,
                ligands
                    .into_iter()
                    .map(|ligand| StereoLigand::new(map_atom(ligand.atom_id), ligand.kind))
                    .collect(),
                attributes,
            )
        })
        .collect();

    let left = Molecule::from_entries(entries);
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: right_atoms,
        bonds: right_bonds,
        dative: right_dative,
        aromatic: right_aromatic,
        multicenter: right_multicenter,
        noncovalent: right_noncovalent,
        stereo_atoms: right_stereo_atoms,
        stereo_bonds: right_stereo_bonds,
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(2), AtomId(0)]),
            },
        )),
    });
    let correspondence = MoleculeRemapping::new(
        GraphRemapping::new(
            Remapping::new(atom_images.iter().copied().map(NodeId::from).collect()).unwrap(),
            Remapping::identity(left.bonds().count()),
        ),
        Remapping::identity(left.dative_bonds().count()),
        Remapping::identity(left.aromatic_systems().count()),
        Remapping::identity(left.multicenter_bonds().count()),
        Remapping::identity(left.noncovalent_bonds().count()),
        Remapping::identity(left.stereo_atoms().count()),
        Remapping::identity(left.stereo_bonds().count()),
    );

    (left, right, correspondence)
}

#[rstest]
fn test_molecule_normalize(#[from(equiv_molecule_entries)] mut entries: MoleculeEntries) {
    entries.dative[0].0.swap(0, 1);
    entries.aromatic[0].0.swap(0, 2);
    entries.multicenter[0].0.swap(0, 2);
    entries.noncovalent[0].0 = [AtomId(3), AtomId(0)];
    entries.stereo_atoms[0].1.swap(0, 1);
    entries.stereo_bonds[0].1.swap(0, 1);

    let mut expected_entries = entries.clone();
    expected_entries.aromatic[0].1.charge = NumForm::Lit(0);
    expected_entries.multicenter[0].1.charge = NumForm::Lit(0);
    let expected = Molecule::from_entries(expected_entries);

    entries.atoms[0].charge = NumForm::lit_set([1_i64]);
    entries.bonds[0].2.order = NumForm::lit_set([1_i64]);
    entries.dative[0].2.order = NumForm::lit_set([1_i64]);
    entries.aromatic[0].1.charge = NumForm::lit_set([0_i64]);
    entries.multicenter[0].1.charge = NumForm::lit_set([0_i64]);
    entries.stereo_atoms[0].2 =
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1]));
    entries.stereo_bonds[0].2 =
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1]));
    let duplicate = entries
        .constraints
        .iter()
        .next()
        .expect("fixture contains one molecule constraint")
        .clone();
    entries.constraints.push(duplicate);

    assert_eq!(Molecule::from_entries(entries).normalize(), Ok(expected));
}

#[rstest]
fn test_molecule_normalize_identity() {
    let molecule = Molecule::default();

    assert_eq!(molecule.clone().normalize(), Ok(molecule));
}

#[rstest]
fn test_molecule_normalize_shared_storage(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
) {
    let mut expected_entries = entries.clone();
    expected_entries.aromatic[0].1.charge = NumForm::Lit(0);
    let expected = Molecule::from_entries(expected_entries);

    entries.atoms[0].charge = NumForm::lit_set([1_i64]);
    entries.aromatic[0].1.charge = NumForm::lit_set([0_i64]);
    let source = Molecule::from_entries(entries);
    let snapshot = source.clone();
    let shared = source.clone();
    assert!(Arc::ptr_eq(&source.atoms, &shared.atoms));

    let normalized = shared.normalize().expect("the molecule is satisfiable");

    assert_eq!(source, snapshot);
    assert_eq!(normalized, expected);
}

#[rstest]
fn test_molecule_normalize_idempotence(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
) {
    entries.atoms[0].charge = NumForm::lit_set([1_i64]);
    let once = Molecule::from_entries(entries)
        .normalize()
        .expect("the molecule is satisfiable");

    assert_eq!(once.clone().normalize(), Ok(once));
}

#[rstest]
#[case::atom(|entries: &mut MoleculeEntries| {
    entries.atoms[0].charge = NumForm::lit_set(Vec::<i64>::new());
})]
#[case::aromatic_system(|entries: &mut MoleculeEntries| {
    entries.aromatic[0].1.charge = NumForm::lit_set(Vec::<i64>::new());
})]
#[case::stereo_atom(|entries: &mut MoleculeEntries| {
    entries.stereo_atoms[0].2 = StereoAtomForm::new(
        StereoKind::Tetrahedral,
        StereoCoset::lit_set(Vec::<u32>::new()),
    );
})]
fn test_molecule_normalize_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] contradict: fn(&mut MoleculeEntries),
) {
    contradict(&mut entries);

    assert_eq!(
        Molecule::from_entries(entries).normalize(),
        Err(Contradiction),
    );
}

#[rstest]
fn test_molecule_reframe_by(#[from(molecule_reframe_source)] source: Molecule) {
    let actions = source.representative_action();
    let reframed = source
        .clone()
        .reframe_by(&actions)
        .expect("the representative action covers the molecule");

    assert_eq!(
        reframed
            .dative_bonds
            .donors(DativeBondId(0))
            .collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2)],
    );
    assert_eq!(
        reframed
            .aromatic_systems
            .atoms(AromaticSystemId(0))
            .collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)],
    );
    assert_eq!(
        reframed.aromatic_systems.attributes(AromaticSystemId(0)),
        &AromaticSystemForm::from_electrons(vec![10, 15, 20]),
    );
    assert_eq!(
        reframed
            .multicenter_bonds
            .atoms(MulticenterBondId(0))
            .collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)],
    );
    assert_eq!(
        reframed.multicenter_bonds.attributes(MulticenterBondId(0)),
        &MulticenterBondForm::from_electrons(vec![10, 12, 14]),
    );
    assert_eq!(
        reframed.noncovalent_bonds.atoms(NoncovalentBondId(0)),
        [AtomId(0), AtomId(3)],
    );
    assert_eq!(
        reframed.stereo_atoms.ligands(StereoAtomId(0)),
        [
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        ],
    );
    assert_eq!(
        reframed.stereo_bonds.ligands(StereoBondId(0)),
        [
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ],
    );
    assert_eq!(
        reframed.constraints.as_slice(),
        &[
            Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::valence(2)),
                    Box::new(AtomConstraintForm::valence(1)),
                ],
            }),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(2usize.into(), 3usize.into()),
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                }),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(1usize.into(), 2usize.into()),
                    relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                }),
            ),
        ],
    );
}

#[rstest]
fn test_molecule_reframe_by_compatible_action(#[from(molecule_reframe_source)] source: Molecule) {
    let actions = source.representative_action();
    let mut compatible = source;
    compatible.atom_mut(AtomId(0)).attributes.charge = NumForm::Lit(2);

    let reframed = compatible
        .reframe_by(&actions)
        .expect("the action domain and participant degrees are compatible");

    assert_eq!(reframed.atom(AtomId(0)).attributes.charge, NumForm::Lit(2));
    assert_eq!(
        reframed.representative_action(),
        reframed.representative_action().identity(),
    );
}

#[rstest]
#[case::missing(|entries: &mut MoleculeEntries| entries.aromatic.clear())]
#[case::degree(|entries: &mut MoleculeEntries| {
    entries.aromatic[0] = (
        vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
        AromaticSystemForm::from_electrons(vec![1, 2, 3, 4]),
    );
})]
fn test_molecule_reframe_by_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] change_provider: fn(&mut MoleculeEntries),
) {
    let source = Molecule::from_entries(entries.clone());
    let mut provider = entries;
    change_provider(&mut provider);
    let actions = Molecule::from_entries(provider).representative_action();

    assert_eq!(source.reframe_by(&actions), None);
}

#[rstest]
fn test_molecule_representative_action(#[from(molecule_reframe_source)] source: Molecule) {
    let actions = source.representative_action();

    assert_eq!(
        actions
            .dative_bonds()
            .action(DativeBondId(0))
            .expect("the dense action covers the dative bond")
            .image(),
        [1, 0],
    );
    assert_eq!(
        actions
            .aromatic_systems()
            .action(AromaticSystemId(0))
            .expect("the dense action covers the aromatic system")
            .image(),
        [1, 2, 0],
    );
    assert_eq!(
        actions
            .multicenter_bonds()
            .action(MulticenterBondId(0))
            .expect("the dense action covers the multicenter bond")
            .image(),
        [2, 0, 1],
    );
    assert_eq!(
        actions.noncovalent_bonds().action(NoncovalentBondId(0)),
        Some(&DynPermutation::try_from(vec![1, 0]).expect("the expected action is valid")),
    );
    assert_eq!(
        actions.stereo_atoms().action(StereoAtomId(0)),
        Some(&Permutation::from_image(&[2, 3, 1, 0])),
    );
    assert_eq!(
        actions.stereo_bonds().action(StereoBondId(0)),
        Some(&Permutation::from_image(&[3, 2, 0, 1])),
    );
}

#[rstest]
fn test_molecule_reframe_with_action(#[from(molecule_reframe_source)] source: Molecule) {
    let (reframed, actions) = source
        .clone()
        .reframe_with_action()
        .expect("the molecule is satisfiable");

    let transported = source
        .clone()
        .normalize()
        .expect("the molecule is satisfiable")
        .reframe_by(&actions)
        .expect("the representative action covers the molecule")
        .normalize()
        .expect("the transported molecule is satisfiable");

    assert_eq!(transported, reframed);
    assert_eq!(source.reframe(), Ok(reframed));
}

#[rstest]
fn test_molecule_reframe_identity() {
    let molecule = Molecule::default();

    assert_eq!(molecule.clone().reframe(), Ok(molecule));
}

#[rstest]
fn test_molecule_reframe_idempotence(#[from(molecule_reframe_source)] source: Molecule) {
    let once = source.reframe().expect("the molecule is satisfiable");

    assert_eq!(once.clone().reframe(), Ok(once));
}

#[rstest]
fn test_molecule_reframe_shared_storage(#[from(molecule_reframe_source)] source: Molecule) {
    let snapshot = source.clone();
    let shared = source.clone();
    assert!(Arc::ptr_eq(&source.atoms, &shared.atoms));

    let reframed = shared.reframe().expect("the molecule is satisfiable");

    assert_eq!(source, snapshot);
    assert!(source.framed_eq(&reframed));
}

#[rstest]
fn test_molecule_reframe_error(#[from(equiv_molecule_entries)] mut entries: MoleculeEntries) {
    entries.stereo_atoms[0].2 = StereoAtomForm::new(
        StereoKind::Tetrahedral,
        StereoCoset::lit_set(Vec::<u32>::new()),
    );

    assert_eq!(
        Molecule::from_entries(entries).reframe(),
        Err(Contradiction)
    );
}

#[rstest]
fn test_molecule_framed_eq(#[from(molecule_reframe_source)] source: Molecule) {
    let representative = source
        .clone()
        .reframe()
        .expect("the molecule is satisfiable");
    let mut different = source.clone();
    different
        .stereo_atoms
        .attributes_mut(StereoAtomId(0))
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32);

    assert!(!source.normalized_eq(&representative));
    assert!(source.framed_eq(&representative));
    assert!(!source.framed_eq(&different));
}

#[rstest]
fn test_molecule_normalized_eq_entity_data(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = Molecule::from_entries(entries.clone());

    let mut canonical_encoding = entries.clone();
    canonical_encoding.atoms[0].charge = NumForm::lit_set([1]);
    let canonical_encoding = Molecule::from_entries(canonical_encoding);
    assert_ne!(base, canonical_encoding);
    assert!(base.normalized_eq(&canonical_encoding));

    let mut differences = Vec::new();

    let mut atom = entries.clone();
    atom.atoms[0].charge = NumForm::Lit(2);
    differences.push(Molecule::from_entries(atom));

    let mut bond = entries.clone();
    bond.bonds[0].2.order = NumForm::Lit(2);
    differences.push(Molecule::from_entries(bond));

    let mut dative = entries.clone();
    dative.dative[0].2.order = NumForm::Lit(2);
    differences.push(Molecule::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic[0].1.electrons = ElectronCountsForm::Lit(vec![2, 0, 1]);
    differences.push(Molecule::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter[0].1.electrons = ElectronCountsForm::Lit(vec![2, 0, 0]);
    differences.push(Molecule::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent[0].1.kind = NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic);
    differences.push(Molecule::from_entries(noncovalent));

    let mut stereo_atom = entries.clone();
    stereo_atom.stereo_atoms[0].2 = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
    differences.push(Molecule::from_entries(stereo_atom));

    let mut stereo_bond = entries.clone();
    stereo_bond.stereo_bonds[0].2 = StereoBondForm::new(StereoKind::CisTrans, 0u32);
    differences.push(Molecule::from_entries(stereo_bond));

    let mut constraint = entries;
    constraint.constraints =
        constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2)]),
        }));
    differences.push(Molecule::from_entries(constraint));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.normalized_eq(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_normalized_eq_relation_frames(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = Molecule::from_entries(entries.clone());
    let mut differences = Vec::new();

    let mut dative = entries.clone();
    dative.dative[0].0 = vec![AtomId(0), AtomId(2)];
    differences.push(Molecule::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic[0].0 = vec![AtomId(0), AtomId(1), AtomId(3)];
    differences.push(Molecule::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter[0].0 = vec![AtomId(0), AtomId(1), AtomId(3)];
    differences.push(Molecule::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent[0].0[1] = AtomId(2);
    differences.push(Molecule::from_entries(noncovalent));

    let mut stereo_atom_site = entries.clone();
    stereo_atom_site.stereo_atoms[0].0 = AtomId(2);
    stereo_atom_site.stereo_atoms[0].1 = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
    ];
    differences.push(Molecule::from_entries(stereo_atom_site));

    let mut stereo_atom_ligand = entries.clone();
    stereo_atom_ligand.stereo_atoms[0].1.swap(2, 3);
    differences.push(Molecule::from_entries(stereo_atom_ligand));

    let mut stereo_bond_site = entries.clone();
    stereo_bond_site.stereo_bonds[0].0 = BondId(2);
    stereo_bond_site.stereo_bonds[0].1 = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
    ];
    differences.push(Molecule::from_entries(stereo_bond_site));

    let mut stereo_bond_ligand = entries;
    stereo_bond_ligand.stereo_bonds[0].1[1] =
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair);
    differences.push(Molecule::from_entries(stereo_bond_ligand));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.normalized_eq(other))
            .collect::<Vec<_>>(),
        vec![false; 8],
    );
}

#[rstest]
fn test_molecule_normalized_eq_structure_and_counts(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = Molecule::from_entries(entries.clone());
    let mut differences = Vec::new();

    let mut topology = entries.clone();
    topology.bonds[2].1 = AtomId(0);
    topology.stereo_bonds[0].1[2] = StereoLigand::new(AtomId(2), StereoLigandKind::LonePair);
    differences.push(Molecule::from_entries(topology));

    let mut atoms = entries.clone();
    atoms.atoms.push(AtomForm::from_element(Element::F));
    differences.push(Molecule::from_entries(atoms));

    let mut bonds = entries.clone();
    bonds
        .bonds
        .push((AtomId(0), AtomId(3), BondForm::from_order(1)));
    differences.push(Molecule::from_entries(bonds));

    let mut dative = entries.clone();
    dative.dative.pop();
    differences.push(Molecule::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic.pop();
    differences.push(Molecule::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter.pop();
    differences.push(Molecule::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent.pop();
    differences.push(Molecule::from_entries(noncovalent));

    let mut stereo_atom = entries.clone();
    stereo_atom.stereo_atoms.pop();
    differences.push(Molecule::from_entries(stereo_atom));

    let mut stereo_bond = entries;
    stereo_bond.stereo_bonds.pop();
    differences.push(Molecule::from_entries(stereo_bond));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.normalized_eq(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_framed_eq_under_entity_ids(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
) {
    let (left, right, correspondence) = case;

    assert!(!left.normalized_eq(&right));
    assert!(left.framed_eq_under(&right, &correspondence));
    assert!(right.framed_eq_under(&left, &correspondence));
}

#[rstest]
fn test_molecule_framed_eq_under_aromatic_system_frame() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 4, 6]),
        )],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(2), AtomId(0), AtomId(1)],
            AromaticSystemForm::from_electrons(vec![6, 2, 4]),
        )],
        ..Default::default()
    });
    let correspondence = MoleculeRemapping::new(
        GraphRemapping::identity(left.atoms().count(), left.bonds().count()),
        Remapping::identity(left.dative_bonds().count()),
        Remapping::identity(left.aromatic_systems().count()),
        Remapping::identity(left.multicenter_bonds().count()),
        Remapping::identity(left.noncovalent_bonds().count()),
        Remapping::identity(left.stereo_atoms().count()),
        Remapping::identity(left.stereo_bonds().count()),
    );

    assert!(!left.normalized_eq(&right));
    assert!(left.framed_eq_under(&right, &correspondence));
}

#[rstest]
fn test_molecule_framed_eq_under_stereo_atom_constraint() {
    let atoms = vec![AtomForm::from_element(Element::C); 5];
    let bonds = (1..=4)
        .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
        .collect::<Vec<_>>();
    let left_ligands = (1..=4)
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let mut right_ligands = left_ligands.clone();
    right_ligands.swap(0, 1);
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: atoms.clone(),
        bonds: bonds.clone(),
        stereo_atoms: vec![(AtomId(0), left_ligands, StereoAtomForm::default())],
        constraints: Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            }),
        )
        .into(),
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        stereo_atoms: vec![(AtomId(0), right_ligands, StereoAtomForm::default())],
        constraints: Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(1usize.into(), 2usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            }),
        )
        .into(),
        ..Default::default()
    });
    let correspondence = MoleculeRemapping::new(
        GraphRemapping::identity(left.atoms().count(), left.bonds().count()),
        Remapping::identity(left.dative_bonds().count()),
        Remapping::identity(left.aromatic_systems().count()),
        Remapping::identity(left.multicenter_bonds().count()),
        Remapping::identity(left.noncovalent_bonds().count()),
        Remapping::identity(left.stereo_atoms().count()),
        Remapping::identity(left.stereo_bonds().count()),
    );

    assert!(left.framed_eq_under(&right, &correspondence));
}

#[rstest]
fn test_molecule_framed_eq_under_stereo_bond_block() {
    let atoms = vec![AtomForm::from_element(Element::C); 4];
    let bonds = vec![
        (AtomId(0), AtomId(1), BondForm::from_order(2)),
        (AtomId(0), AtomId(2), BondForm::from_order(1)),
        (AtomId(1), AtomId(3), BondForm::from_order(1)),
    ];
    let atom = |id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom);
    let lone_pair = |id| StereoLigand::new(AtomId(id), StereoLigandKind::LonePair);
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: atoms.clone(),
        bonds: bonds.clone(),
        stereo_bonds: vec![(
            BondId(0),
            vec![atom(2), lone_pair(0), atom(3), lone_pair(1)],
            StereoBondForm::default(),
        )],
        constraints: Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            }),
        )
        .into(),
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        stereo_bonds: vec![(
            BondId(0),
            vec![atom(3), lone_pair(1), atom(2), lone_pair(0)],
            StereoBondForm::default(),
        )],
        constraints: Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(2usize.into(), 3usize.into()),
                relation: TopicityRelationForm::Lit(Topicity::Homotopic),
            }),
        )
        .into(),
        ..Default::default()
    });
    let correspondence = MoleculeRemapping::new(
        GraphRemapping::identity(left.atoms().count(), left.bonds().count()),
        Remapping::identity(left.dative_bonds().count()),
        Remapping::identity(left.aromatic_systems().count()),
        Remapping::identity(left.multicenter_bonds().count()),
        Remapping::identity(left.noncovalent_bonds().count()),
        Remapping::identity(left.stereo_atoms().count()),
        Remapping::identity(left.stereo_bonds().count()),
    );

    assert!(!left.normalized_eq(&right));
    assert!(left.framed_eq_under(&right, &correspondence));
}

/// A correspondence must map each matched entity's participants onto its counterpart's. That is a
/// property of the correspondence, checked once, not of any per-entity-kind payload comparison.
///
/// Covers all six overlay kinds. Each case moves one participant of its entity onto a different
/// atom that is still legal for that entity kind, so both sides satisfy integrity and only the
/// correspondence is at fault.
#[rustfmt::skip]
#[rstest]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
#[case::dative(Entity::DativeBond(DativeBondId(0)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
fn test_molecule_framed_eq_under_participant_mismatch_error(#[case] entity: Entity) {
    // Atom 0 is bonded to 1, 4 and 5; atom 1 to 0, 2 and 5; atom 2 to 1 and 3. That gives every
    // entity kind a legal participant to move to.
    let build = |shifted: bool| {
        let atom = |kept: u32, moved: u32| AtomId(if shifted { moved } else { kept });
        let ligand = |kept: u32, moved: u32| StereoLigand::new(atom(kept, moved), StereoLigandKind::Atom);
        let virtual_ligand =
            |site: u32, kind| StereoLigand::new(AtomId(site), kind);
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5), (1, 5)]
                .into_iter()
                .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        };
        match entity {
            Entity::AromaticSystem(_) => {
                entries.aromatic = vec![(
                    vec![AtomId(1), AtomId(2), atom(3, 4)],
                    AromaticSystemForm::default(),
                )]
            }
            Entity::MulticenterBond(_) => {
                entries.multicenter = vec![(
                    vec![AtomId(1), AtomId(2), atom(3, 4)],
                    MulticenterBondForm::default(),
                )]
            }
            Entity::NoncovalentBond(_) => {
                entries.noncovalent =
                    vec![([AtomId(2), atom(4, 5)], NoncovalentBondForm::default())]
            }
            Entity::DativeBond(_) => {
                entries.dative = vec![(
                    vec![AtomId(1), atom(4, 5)],
                    AtomId(2),
                    DativeBondForm::default(),
                )]
            }
            Entity::StereoAtom(_) => {
                // Site 0, ligands drawn from its neighbours 1, 4, 5 plus two virtuals.
                entries.stereo_atoms = vec![(
                    AtomId(0),
                    vec![
                        ligand(1, 5),
                        ligand(4, 4),
                        virtual_ligand(0, StereoLigandKind::ImplicitHydrogen),
                        virtual_ligand(0, StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::default(),
                )]
            }
            Entity::StereoBond(_) => {
                // Site is bond 1 (atoms 1-2); endpoint 1's ligands come from 0 or 5, endpoint 2's
                // from 3.
                entries.stereo_bonds = vec![(
                    BondId(1),
                    vec![
                        ligand(0, 5),
                        virtual_ligand(1, StereoLigandKind::LonePair),
                        ligand(3, 3),
                        virtual_ligand(2, StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::default(),
                )]
            }
            _ => unreachable!("test cases contain only overlay kinds"),
        }
        Molecule::try_from_entries(entries)
    };

    let (Ok(left), Ok(right)) = (build(false), build(true)) else {
        panic!("{entity}: both molecules must satisfy integrity");
    };
    assert_ne!(left, right, "{entity}: the shift must actually move a participant");

    let correspondence = MoleculeRemapping::new(GraphRemapping::identity(left.atoms().count(), left.bonds().count()), Remapping::identity(left.dative_bonds().count()), Remapping::identity(left.aromatic_systems().count()), Remapping::identity(left.multicenter_bonds().count()), Remapping::identity(left.noncovalent_bonds().count()), Remapping::identity(left.stereo_atoms().count()), Remapping::identity(left.stereo_bonds().count()));

    assert!(
        !left.framed_eq_under(&right, &correspondence),
        "{entity}: participants that do not correspond must not compare equivalent",
    );
}

#[rstest]
fn test_molecule_framed_eq_under_count_error(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
) {
    let (left, right, _) = case;
    assert!(!left.framed_eq_under(&right, &MoleculeRemapping::default()));
}

#[rstest]
fn test_molecule_framed_eq_under_entity_id_mismatch_error(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
) {
    let (left, right, _) = case;
    let inconsistent = MoleculeRemapping::new(
        GraphRemapping::new(
            Remapping::new(vec![NodeId(2), NodeId(3), NodeId(0), NodeId(1)]).unwrap(),
            Remapping::new(vec![EdgeId(1), EdgeId(0), EdgeId(2)]).unwrap(),
        ),
        Remapping::identity(left.dative_bonds().count()),
        Remapping::identity(left.aromatic_systems().count()),
        Remapping::identity(left.multicenter_bonds().count()),
        Remapping::identity(left.noncovalent_bonds().count()),
        Remapping::identity(left.stereo_atoms().count()),
        Remapping::identity(left.stereo_bonds().count()),
    );
    assert!(!left.framed_eq_under(&right, &inconsistent));
}

#[rstest]
fn test_molecule_remap(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
) {
    let (left, right, remapping) = case;
    assert_eq!(left.remap(&remapping), right);
}

#[rstest]
#[case::atoms(0)]
#[case::bonds(1)]
#[case::dative_bonds(2)]
#[case::aromatic_systems(3)]
#[case::multicenter_bonds(4)]
#[case::noncovalent_bonds(5)]
#[case::stereo_atoms(6)]
#[case::stereo_bonds(7)]
fn test_molecule_try_remap_error(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
    #[case] kind: usize,
) {
    let (left, _, _) = case;
    let remapping = MoleculeRemapping::new(
        GraphRemapping::new(
            Remapping::identity(left.atoms().count() + usize::from(kind == 0)),
            Remapping::identity(left.bonds().count() + usize::from(kind == 1)),
        ),
        Remapping::identity(left.dative_bonds().count() + usize::from(kind == 2)),
        Remapping::identity(left.aromatic_systems().count() + usize::from(kind == 3)),
        Remapping::identity(left.multicenter_bonds().count() + usize::from(kind == 4)),
        Remapping::identity(left.noncovalent_bonds().count() + usize::from(kind == 5)),
        Remapping::identity(left.stereo_atoms().count() + usize::from(kind == 6)),
        Remapping::identity(left.stereo_bonds().count() + usize::from(kind == 7)),
    );
    assert_eq!(left.try_remap(&remapping), None);
}

#[rstest]
#[should_panic(expected = "molecule remapping requires matching entity counts")]
fn test_molecule_remap_error(
    #[from(framed_eq_under_molecules)] case: (Molecule, Molecule, MoleculeRemapping),
) {
    let (left, _, _) = case;
    left.remap(&MoleculeRemapping::default());
}

#[rstest]
#[case::c_c(BondId(0), AtomId(0), AtomId(1), NumForm::Lit(1))]
#[case::c_n(BondId(1), AtomId(1), AtomId(2), NumForm::Lit(2))]
#[case::n_o(BondId(2), AtomId(2), AtomId(3), NumForm::Lit(1))]
fn test_molecule_bond(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] id: BondId,
    #[case] first: AtomId,
    #[case] second: AtomId,
    #[case] order: NumForm,
) {
    let bv = molecule.bond(id);
    assert_eq!(bv.id, id);
    assert_eq!(bv.atom_ids()[0], first);
    assert_eq!(bv.atom_ids()[1], second);
    assert_eq!(bv.attributes.order, order);
}

#[rstest]
fn test_molecule_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(BondId, AtomId, AtomId, NumForm)> = molecule
        .bonds()
        .iter()
        .map(|v| {
            (
                v.id,
                v.atom_ids()[0],
                v.atom_ids()[1],
                v.attributes.order.clone(),
            )
        })
        .collect();
    assert_eq!(
        projected,
        vec![
            (BondId(0), AtomId(0), AtomId(1), NumForm::Lit(1)),
            (BondId(1), AtomId(1), AtomId(2), NumForm::Lit(2)),
            (BondId(2), AtomId(2), AtomId(3), NumForm::Lit(1)),
        ]
    );
}

#[rstest]
fn test_molecule_dative_bond(#[from(rich_molecule)] molecule: Molecule) {
    let dv = molecule.dative_bond(DativeBondId(0));
    assert_eq!(dv.id, DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(3));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(2)]);
    assert_eq!(
        dv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(2), AtomId(3)]
    );
    assert_eq!(dv.attributes.order, NumForm::Lit(1));
}

#[rstest]
fn test_molecule_dative_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(DativeBondId, Vec<AtomId>, AtomId)> = molecule
        .dative_bonds()
        .iter()
        .map(|v| (v.id, v.donor_ids().collect(), v.acceptor_id()))
        .collect();
    assert_eq!(
        projected,
        vec![(DativeBondId(0), vec![AtomId(2)], AtomId(3))]
    );
}

#[rstest]
fn test_molecule_aromatic_system(#[from(rich_molecule)] molecule: Molecule) {
    let av = molecule.aromatic_system(AromaticSystemId(0));
    assert_eq!(av.id, AromaticSystemId(0));
    assert_eq!(
        av.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
    assert_eq!(
        av.bond_ids().collect::<Vec<_>>(),
        vec![BondId(0), BondId(1)]
    );
}

#[rstest]
fn test_molecule_aromatic_systems(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(AromaticSystemId, Vec<AtomId>, Vec<BondId>)> = molecule
        .aromatic_systems()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect(), v.bond_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(
            AromaticSystemId(0),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![BondId(0), BondId(1)],
        )]
    );
}

#[rstest]
fn test_molecule_multicenter_bond(#[from(rich_molecule)] molecule: Molecule) {
    let mv = molecule.multicenter_bond(MulticenterBondId(0));
    assert_eq!(mv.id, MulticenterBondId(0));
    assert_eq!(
        mv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
}

#[rstest]
fn test_molecule_multicenter_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(MulticenterBondId, Vec<AtomId>)> = molecule
        .multicenter_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)],)]
    );
}

#[rstest]
fn test_molecule_noncovalent_bond(#[from(rich_molecule)] molecule: Molecule) {
    let nv = molecule.noncovalent_bond(NoncovalentBondId(0));
    assert_eq!(nv.id, NoncovalentBondId(0));
    assert_eq!(nv.atom_ids(), [AtomId(0), AtomId(3)]);
}

#[rstest]
fn test_molecule_noncovalent_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(NoncovalentBondId, [AtomId; 2])> = molecule
        .noncovalent_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids()))
        .collect();
    assert_eq!(
        projected,
        vec![(NoncovalentBondId(0), [AtomId(0), AtomId(3)])]
    );
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_of_id(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(molecule.bonds().of_id(a, b), expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
#[case::wrong_donor(AtomId(3), vec![AtomId(1)], None)]
fn test_dative_bond_views_of_id(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(molecule.dative_bonds().of_id(acceptor, &donors), expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_of_id(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(molecule.noncovalent_bonds().of_id(a, b), expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let inc: Vec<_> = molecule.dative_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let inc: Vec<_> = molecule.aromatic_systems().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let inc: Vec<_> = molecule.multicenter_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let inc: Vec<_> = molecule.noncovalent_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
#[case::disjoint(vec![AtomId(0), AtomId(1)], vec![])]
fn test_dative_bond_views_induced_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    assert_eq!(molecule.dative_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_aromatic_system_views_induced_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    assert_eq!(molecule.aromatic_systems().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_multicenter_bond_views_induced_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    assert_eq!(molecule.multicenter_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(1), AtomId(2)], vec![])]
fn test_noncovalent_bond_views_induced_ids(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    assert_eq!(molecule.noncovalent_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_of(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(molecule.bonds().of(a, b).map(|v| v.id), expected);
}

#[rstest]
#[case::pair(vec![AtomId(0), AtomId(1)], vec![BondId(0)])]
#[case::triangle(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1)])]
#[case::singleton(vec![AtomId(0)], vec![])]
fn test_bond_views_induced(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<BondId>,
) {
    let mut got: Vec<BondId> = molecule
        .bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    got.sort_unstable();
    assert_eq!(got, expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = molecule
        .dative_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
fn test_dative_bond_views_of(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(
        molecule.dative_bonds().of(acceptor, &donors).map(|v| v.id),
        expected
    );
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
fn test_dative_bond_views_induced(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = molecule
        .dative_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = molecule
        .aromatic_systems()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_aromatic_system_views_of(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(
        molecule.aromatic_systems().of(atoms).map(|v| v.id),
        expected
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_aromatic_system_views_induced(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = molecule
        .aromatic_systems()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = molecule
        .multicenter_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_multicenter_bond_views_of(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(
        molecule.multicenter_bonds().of(atoms).map(|v| v.id),
        expected,
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_multicenter_bond_views_induced(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = molecule
        .multicenter_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = molecule
        .noncovalent_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_of(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(
        molecule.noncovalent_bonds().of(a, b).map(|v| v.id),
        expected,
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_noncovalent_bond_views_induced(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = molecule
        .noncovalent_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::atom_0(AtomId(0), Element::C)]
#[case::atom_1(AtomId(1), Element::C)]
#[case::atom_2(AtomId(2), Element::N)]
#[case::atom_3(AtomId(3), Element::O)]
fn test_molecule_atom(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] id: AtomId,
    #[case] element: Element,
) {
    let av = molecule.atom(id);
    assert_eq!(av.id, id);
    assert_eq!(av.attributes.element, ElementForm::Lit(element));
}

#[rstest]
fn test_molecule_is_empty() {
    assert!(Molecule::default().is_empty());
}

#[rstest]
fn test_molecule_is_empty_rich(#[from(rich_molecule)] molecule: Molecule) {
    assert!(!molecule.is_empty());
}

#[rstest]
fn test_molecule_has_constraints_empty() {
    assert!(!Molecule::default().has_constraints());
}

#[rstest]
fn test_molecule_has_constraints_rich(#[from(rich_molecule)] molecule: Molecule) {
    assert!(!molecule.has_constraints());
}

#[rstest]
fn test_molecule_has_dative_bonds(#[from(rich_molecule)] molecule: Molecule) {
    assert!(molecule.has_dative_bonds());
}

#[rstest]
fn test_molecule_has_aromatic_systems(#[from(rich_molecule)] molecule: Molecule) {
    assert!(molecule.has_aromatic_systems());
}

#[rstest]
fn test_molecule_has_multicenter_bonds(#[from(rich_molecule)] molecule: Molecule) {
    assert!(molecule.has_multicenter_bonds());
}

#[rstest]
fn test_molecule_has_noncovalent_bonds(#[from(rich_molecule)] molecule: Molecule) {
    assert!(molecule.has_noncovalent_bonds());
}

#[rstest]
fn test_molecule_has_stereo_atoms() {
    let molecule = mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    );
    assert!(molecule.has_stereo_atoms());
    assert!(!molecule.has_stereo_bonds());
}

#[rstest]
fn test_molecule_has_stereo_bonds() {
    let molecule = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    );
    assert!(molecule.has_stereo_bonds());
    assert!(!molecule.has_stereo_atoms());
}

#[rstest]
fn test_molecule_atoms(#[from(rich_molecule)] molecule: Molecule) {
    let projected: Vec<(AtomId, ElementForm)> = molecule
        .atoms()
        .iter()
        .map(|v| (v.id, v.attributes.element.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (AtomId(0), ElementForm::Lit(Element::C)),
            (AtomId(1), ElementForm::Lit(Element::C)),
            (AtomId(2), ElementForm::Lit(Element::N)),
            (AtomId(3), ElementForm::Lit(Element::O)),
        ]
    );
}

#[test]
fn test_bond_views_induced_ids() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    });
    let bonds = molecule.bonds().induced_ids(&[AtomId(0), AtomId(1)]);
    assert_eq!(bonds, vec![BondId(0)]);

    let mut all = molecule
        .bonds()
        .induced_ids(&[AtomId(0), AtomId(1), AtomId(2)]);
    all.sort_unstable();
    assert_eq!(all, vec![BondId(0), BondId(1), BondId(2)]);
}

fn chain(n: usize) -> Molecule {
    let atoms = vec![AtomForm::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n.saturating_sub(1))
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId((i + 1) as u32),
                BondForm::from_order(1),
            )
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn ring(n: usize) -> Molecule {
    let atoms = vec![AtomForm::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n)
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId(((i + 1) % n) as u32),
                BondForm::from_order(1),
            )
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn two_components() -> Molecule {
    let atoms = vec![AtomForm::from_element(Element::C); 4];
    let bonds = vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1)),
        (AtomId(2), AtomId(3), BondForm::from_order(1)),
    ];
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

#[rstest]
#[case::isolated(chain(1), AtomId(0), 0)]
#[case::chain_end(chain(3), AtomId(0), 1)]
#[case::chain_mid(chain(3), AtomId(1), 2)]
#[case::ring_vertex(ring(6), AtomId(0), 2)]
fn test_molecule_degree(#[case] molecule: Molecule, #[case] atom: AtomId, #[case] expected: usize) {
    assert_eq!(molecule.graph().degree(atom), expected);
}

#[rstest]
#[case::single(chain(3), 1)]
#[case::two(two_components(), 2)]
#[case::empty(Molecule::default(), 0)]
fn test_molecule_enumerate_connected_components(
    #[case] molecule: Molecule,
    #[case] expected: usize,
) {
    let cc = molecule
        .graph()
        .enumerate_connected_components(ConnectedComponentsAlgorithm::Bfs);
    assert_eq!(cc.len(), expected);
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain(chain(5), 0)]
fn test_molecule_enumerate_biconnected_components(
    #[case] molecule: Molecule,
    #[case] expected: usize,
) {
    let bcc = molecule
        .graph()
        .enumerate_biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
    assert_eq!(bcc.len(), expected);
}

#[rstest]
#[case::ring_bond(ring(6), BondId(0), Some(6))]
#[case::chain_bond(chain(3), BondId(0), None)]
fn test_molecule_shortest_cycle_through_bond(
    #[case] molecule: Molecule,
    #[case] bond: BondId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        molecule
            .graph()
            .shortest_cycle_through_bond(bond, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::ring_atom(ring(6), AtomId(0), Some(6))]
#[case::chain_atom(chain(3), AtomId(1), None)]
fn test_molecule_shortest_cycle_through_atom(
    #[case] molecule: Molecule,
    #[case] atom: AtomId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        molecule
            .graph()
            .shortest_cycle_through_atom(atom, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
fn test_graph_view_visit_simple_cycles(
    #[case] molecule: Molecule,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    let mut cycles: Vec<Vec<AtomId>> = Vec::new();
    let flow: ControlFlow<()> = molecule.graph().visit_simple_cycles(
        max_size,
        SimpleCycleEnumerationAlgorithm::ReadTarjan,
        |cycle| {
            cycles.push(cycle.to_vec());
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    assert_eq!(cycles, expected);
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
#[case::chain(chain(5), 10, vec![])]
#[case::empty(Molecule::default(), 10, vec![])]
fn test_graph_view_enumerate_simple_cycles(
    #[case] molecule: Molecule,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    assert_eq!(
        molecule
            .graph()
            .enumerate_simple_cycles(max_size, SimpleCycleEnumerationAlgorithm::ReadTarjan),
        expected
    );
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
fn test_graph_view_visit_relevant_cycles(
    #[case] molecule: Molecule,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    let mut cycles: Vec<Vec<AtomId>> = Vec::new();
    let flow: ControlFlow<()> = molecule.graph().visit_relevant_cycles(
        max_size,
        RelevantCycleEnumerationAlgorithm::Vismara,
        |cycle| {
            cycles.push(cycle.to_vec());
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    assert_eq!(cycles, expected);
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
#[case::chain(chain(5), 10, vec![])]
#[case::empty(Molecule::default(), 10, vec![])]
fn test_graph_view_enumerate_relevant_cycles(
    #[case] molecule: Molecule,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    assert_eq!(
        molecule
            .graph()
            .enumerate_relevant_cycles(max_size, RelevantCycleEnumerationAlgorithm::Vismara),
        expected
    );
}

#[rstest]
#[case::triangle(ring(3), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_maximum_independent_set(#[case] molecule: Molecule, #[case] expected: usize) {
    let mis = molecule
        .graph()
        .maximum_independent_set(MaximumIndependentSetAlgorithm::BranchAndBound);
    assert_eq!(mis.len(), expected);
}

#[rstest]
#[case::ring_6(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 4 :single] [4 5 :single] [5 0 :single]]}"#),
    vec![BondId(0), BondId(2), BondId(4)],
)]
fn test_graph_view_bipartite_maximum_matching(
    #[case] molecule: Molecule,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = molecule.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        molecule
            .graph()
            .bipartite_maximum_matching(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
            )
            .expect("case graph is bipartite")
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::triangle(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#))]
fn test_graph_view_bipartite_maximum_matching_error(#[case] molecule: Molecule) {
    let node_order: Vec<AtomId> = molecule.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        molecule
            .graph()
            .bipartite_maximum_matching(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
            )
            .unwrap_err(),
        NonBipartiteGraphError,
    );
}

#[rstest]
#[case::chain_4(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single]]}"#),
    vec![BondId(0), BondId(2)],
)]
#[case::triangle(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#),
    vec![BondId(0)],
)]
#[case::single(mol_dsl!(r#"{:atoms ["C"]}"#), vec![])]
fn test_graph_view_general_maximum_matching(
    #[case] molecule: Molecule,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = molecule.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        molecule
            .graph()
            .general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds)
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::bipartite(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 0 :single]]}"#),
    vec![BondId(0), BondId(2)],
)]
#[case::general(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#),
    vec![BondId(0)],
)]
fn test_graph_view_bipartite_maximum_matching_or_general(
    #[case] molecule: Molecule,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = molecule.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        molecule
            .graph()
            .bipartite_maximum_matching_or_general(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                GeneralMaximumMatchingAlgorithm::Edmonds,
            )
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::ring_6(
    ring(6),
    vec![
        vec![BondId(0), BondId(2), BondId(4)],
        vec![BondId(1), BondId(3), BondId(5)],
    ],
)]
fn test_graph_view_visit_perfect_matchings(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Vec<BondId>>,
) {
    let mut matchings: Vec<Vec<BondId>> = Vec::new();
    let flow: ControlFlow<()> = molecule.graph().visit_perfect_matchings(
        MatchingEnumerationAlgorithm::BranchAndBound,
        |matching| {
            let mut bonds: Vec<BondId> = matching.bonds().collect();
            bonds.sort_unstable();
            matchings.push(bonds);
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    matchings.sort_unstable();
    assert_eq!(matchings, expected);
}

#[rstest]
#[case::ring_6(ring(6), 2)]
fn test_molecule_enumerate_perfect_matchings(#[case] molecule: Molecule, #[case] expected: usize) {
    let ms = molecule
        .graph()
        .enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound);
    assert_eq!(ms.len(), expected);
    for m in &ms {
        assert!(m.is_perfect(molecule.atoms().count()));
    }
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_automorphisms(#[case] molecule: Molecule, #[case] expected_orbits: usize) {
    let auto = molecule
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert_eq!(auto.orbit_count(), expected_orbits);
    assert_eq!(auto.atom_count(), molecule.atoms().count());
}

#[test]
fn test_atom_automorphism_same_orbit() {
    let molecule = ring(6);
    let auto = molecule
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert!(auto.same_orbit(AtomId(0), AtomId(3)));
}

#[rstest]
fn test_graph_view_visit_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let mut matches: Vec<Vec<AtomId>> = Vec::new();
    let flow: ControlFlow<()> = target.graph().visit_subgraph_isomorphisms(
        &query.graph(),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
        |embedding| {
            matches.push(embedding.to_vec());
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(0), AtomId(5)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
            vec![AtomId(3), AtomId(2)],
            vec![AtomId(3), AtomId(4)],
            vec![AtomId(4), AtomId(3)],
            vec![AtomId(4), AtomId(5)],
            vec![AtomId(5), AtomId(0)],
            vec![AtomId(5), AtomId(4)],
        ]
    );
}

#[rstest]
fn test_molecule_enumerate_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().enumerate_subgraph_isomorphisms(
        &query.graph(),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(0), AtomId(5)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
            vec![AtomId(3), AtomId(2)],
            vec![AtomId(3), AtomId(4)],
            vec![AtomId(4), AtomId(3)],
            vec![AtomId(4), AtomId(5)],
            vec![AtomId(5), AtomId(0)],
            vec![AtomId(5), AtomId(4)],
        ]
    );
}

#[rstest]
fn test_graph_view_visit_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let mut matches: Vec<Vec<AtomId>> = Vec::new();
    let flow: ControlFlow<()> = target.graph().visit_subgraph_isomorphisms_at(
        &query.graph(),
        (AtomId(0), AtomId(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
        |embedding| {
            matches.push(embedding.to_vec());
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(5)]]
    );
}

#[rstest]
fn test_molecule_enumerate_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().enumerate_subgraph_isomorphisms_at(
        &query.graph(),
        (AtomId(0), AtomId(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(5)],]
    );
}

#[rstest]
fn test_molecule_induced_subgraph(#[from(rich_molecule)] molecule: Molecule) {
    let sub = molecule.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    let extracted = molecule.extract(&sub);
    let atom_elements: Vec<_> = extracted
        .atoms()
        .iter()
        .map(|v| v.attributes.element.clone())
        .collect();
    assert_eq!(
        atom_elements,
        vec![
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::N),
        ]
    );
    let bonds: Vec<(AtomId, AtomId, NumForm)> = extracted
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.attributes.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), NumForm::Lit(1)),
            (AtomId(1), AtomId(2), NumForm::Lit(2)),
        ]
    );
    assert_eq!(
        sub.atoms().matched_pairs(),
        &[
            (AtomId(0), AtomId(0)),
            (AtomId(1), AtomId(1)),
            (AtomId(2), AtomId(2))
        ]
    );
    assert_eq!(
        sub.bonds().matched_pairs(),
        &[(BondId(0), BondId(0)), (BondId(1), BondId(1))]
    );
    assert_eq!(
        sub.aromatic_systems().matched_pairs(),
        &[(AromaticSystemId(0), AromaticSystemId(0))]
    );
    assert_eq!(
        sub.multicenter_bonds().matched_pairs(),
        &[(MulticenterBondId(0), MulticenterBondId(0))]
    );
    assert_eq!(
        sub.dative_bonds().matched_pairs(),
        &[] as &[(DativeBondId, DativeBondId)]
    );
    assert_eq!(
        sub.noncovalent_bonds().matched_pairs(),
        &[] as &[(NoncovalentBondId, NoncovalentBondId)]
    );
}

#[rstest]
fn test_molecule_induced_subgraph_preserves_dative(#[from(rich_molecule)] molecule: Molecule) {
    let sub = molecule.induced_subgraph(&[AtomId(2), AtomId(3)]);
    assert_eq!(
        sub.atoms().matched_pairs(),
        &[(AtomId(0), AtomId(2)), (AtomId(1), AtomId(3))]
    );
    assert_eq!(
        sub.dative_bonds().matched_pairs(),
        &[(DativeBondId(0), DativeBondId(0))]
    );
    let extracted = molecule.extract(&sub);
    let dv = extracted.dative_bond(DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(1));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(0)]);
    assert_eq!(dv.attributes.order, NumForm::Lit(1));
}

#[rstest]
fn test_molecule_edits(#[from(rich_molecule)] molecule: Molecule) {
    let sub = molecule.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    assert_eq!(
        molecule.edits(&sub),
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(3))],
            bonds: vec![BondHandle::Id(BondId(2))],
        }])
    );
}

#[rstest]
fn test_molecule_edits_identity(#[from(rich_molecule)] molecule: Molecule) {
    let atom_ids: Vec<AtomId> = molecule.atoms().iter().map(|v| v.id).collect();
    let sub = molecule.induced_subgraph(&atom_ids);
    assert_eq!(molecule.edits(&sub), Edits::new());
}

#[rstest]
#[case::add_atom(
    mol_dsl!(r#"{:atoms ["C"]}"#),
    Edits::from_iter([Edit::AddAtoms {
        atoms: vec![AtomForm::from_element(Element::N)],
    }]),
    mol_dsl!(r#"{:atoms ["C" "N"]}"#),
)]
fn test_molecule_apply(
    #[case] molecule: Molecule,
    #[case] edits: Edits,
    #[case] expected: Molecule,
) {
    let original = molecule.clone();

    assert_eq!(molecule.apply(edits), Ok(expected));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::transaction(
    mol_dsl!(r#"{:atoms ["C"]}"#),
    Edits::from_iter([
        Edit::AddAtoms {
            atoms: vec![AtomForm::from_element(Element::N)],
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        },
    ]),
    MoleculeApplyError::Transaction(TransactionError::OldStateMismatch),
)]
#[case::integrity(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
    Edits::from_iter([Edit::AddBonds {
        bonds: vec![AddBond {
            endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: BondForm::from_order(1),
        }],
    }]),
    MoleculeApplyError::Integrity(MoleculeIntegrityError::BondsParallel {
        atoms: [AtomId(0), AtomId(1)],
    }),
)]
fn test_molecule_apply_error(
    #[case] molecule: Molecule,
    #[case] edits: Edits,
    #[case] expected: MoleculeApplyError,
) {
    let original = molecule.clone();

    assert_eq!(molecule.tracked_apply(edits.clone()), Err(expected.clone()));
    assert_eq!(molecule.apply(edits), Err(expected));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::identity(false)]
#[case::remove_all(true)]
fn test_molecule_tracked_apply(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] remove_all: bool,
) {
    let source = Molecule::from_entries(entries);
    let mut edits = Edits::new();
    if remove_all {
        edits.remove_topology(
            (0..4).map(|idx| AtomHandle::Id(AtomId(idx))).collect(),
            vec![],
        );
    }
    let expected = if remove_all {
        Molecule::new()
    } else {
        source.clone()
    };
    let witness = if remove_all {
        MoleculeCorrespondence::new(
            Correspondence::new(vec![], 4, 0).unwrap(),
            Correspondence::new(vec![], 3, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
            Correspondence::new(vec![], 1, 0).unwrap(),
        )
    } else {
        MoleculeCorrespondence::new(
            Correspondence::identity(4),
            Correspondence::identity(3),
            Correspondence::identity(1),
            Correspondence::identity(1),
            Correspondence::identity(1),
            Correspondence::identity(1),
            Correspondence::identity(1),
            Correspondence::identity(1),
        )
    };
    assert_eq!(source.apply(edits.clone()), Ok(expected.clone()));
    assert_eq!(source.tracked_apply(edits), Ok((expected, witness)));
}

#[rstest]
fn test_molecule_tracked_extract(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![entries.atoms[0].clone(), entries.atoms[3].clone()],
        noncovalent: vec![(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        ..Default::default()
    });
    let source = Molecule::from_entries(entries);
    let selection = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(3)), (AtomId(1), AtomId(0))], 2, 4).unwrap(),
        Correspondence::new(vec![], 0, 3).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
    );
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(
            Compaction::new(4, vec![NodeId(1), NodeId(2)]).unwrap(),
            Compaction::new(3, vec![EdgeId(0), EdgeId(1), EdgeId(2)]).unwrap(),
        ),
        Compaction::new(1, vec![DativeBondId(0)]).unwrap(),
        Compaction::new(1, vec![AromaticSystemId(0)]).unwrap(),
        Compaction::new(1, vec![MulticenterBondId(0)]).unwrap(),
        Compaction::identity(1),
        Compaction::new(1, vec![StereoAtomId(0)]).unwrap(),
        Compaction::new(1, vec![StereoBondId(0)]).unwrap(),
    );

    assert_eq!(source.extract(&selection), expected);
    assert_eq!(
        source.tracked_extract(&selection),
        (expected, expected_compaction)
    );
}

#[rstest]
fn test_molecule_tracked_extract_empty() {
    let source = Molecule::new();
    let selection = MoleculeCorrespondence::empty();

    assert_eq!(source.extract(&selection), Molecule::new());
    assert_eq!(
        source.tracked_extract(&selection),
        (Molecule::new(), MoleculeCompaction::empty())
    );
}

#[rstest]
fn test_molecule_tracked_extract_identity(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let selection = source.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2), AtomId(3)]);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(4), Compaction::identity(3)),
        Compaction::identity(1),
        Compaction::identity(1),
        Compaction::identity(1),
        Compaction::identity(1),
        Compaction::identity(1),
        Compaction::identity(1),
    );

    assert_eq!(source.extract(&selection), source);
    assert_eq!(
        source.tracked_extract(&selection),
        (source, expected_compaction)
    );
}

#[rstest]
fn test_molecule_extract(#[from(rich_molecule)] molecule: Molecule) {
    let sub = molecule.induced_subgraph(&[AtomId(0), AtomId(1)]);
    let extracted = molecule.extract(&sub);
    assert_eq!(extracted.atoms().count(), 2);
}

#[rstest]
fn test_transaction_tracked_rollback_error() {
    let source = mol_dsl!(r#"{:atoms ["C"]}"#);
    let mut editor = source.edit();
    let mut edits = Edits::new();
    edits.add_atom(AtomForm::from_element(Element::N));
    let (transaction, _) = editor.tracked_transact(edits).unwrap();
    editor.remove(&[AtomId(1)], &[]);
    let before = editor.tracked_snapshot().unwrap();
    let mut plain = editor.clone();
    assert_eq!(
        transaction.clone().rollback(&mut plain),
        Err(TransactionError::RollbackStateMismatch)
    );
    assert_eq!(
        transaction.tracked_rollback(&mut editor),
        Err(TransactionError::RollbackStateMismatch)
    );
    assert_eq!(plain.tracked_snapshot(), Ok(before.clone()));
    assert_eq!(editor.tracked_snapshot(), Ok(before));
}

#[rstest]
fn test_molecule_editor_tracked_transact(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let source = Molecule::from_entries(entries);
    let mut editor = source.edit();
    editor.add_atom(AtomForm::from_element(Element::F));
    let before = editor.tracked_snapshot().unwrap();
    let mut plain = editor.clone();
    let mut edits = Edits::new();
    edits.remove_topology(
        (0..4).map(|idx| AtomHandle::Id(AtomId(idx))).collect(),
        vec![],
    );
    let added = edits.add_atom(AtomForm::from_element(Element::Cl));
    edits.remove_atom(added);
    edits.add_atom(AtomForm::from_element(Element::N));
    let expected = mol_dsl!(r#"{:atoms ["F" "N"]}"#);
    let expected_witness = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(4), AtomId(0))], 5, 2).unwrap(),
        Correspondence::new(vec![], 3, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
        Correspondence::new(vec![], 1, 0).unwrap(),
    );
    let transaction = plain.transact(edits.clone()).unwrap();
    let (tracked_transaction, witness) = editor.tracked_transact(edits).unwrap();
    assert_eq!(tracked_transaction, transaction);
    assert_eq!(witness, expected_witness);
    assert_eq!(editor.snapshot(), Ok(expected));
    assert_eq!(editor.tracked_snapshot(), plain.tracked_snapshot());
    assert_eq!(
        editor.tracked_snapshot().unwrap().1,
        before.1.compose(&witness).unwrap()
    );

    let reverse = tracked_transaction.tracked_rollback(&mut editor).unwrap();
    transaction.rollback(&mut plain).unwrap();
    assert_eq!(reverse, witness.reverse());
    assert_eq!(editor.snapshot(), Ok(before.0));
    assert_eq!(editor.tracked_snapshot(), plain.tracked_snapshot());
}

#[rstest]
fn test_molecule_editor_tracked_transact_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let mut editor = source.edit();
    editor.add_atom(AtomForm::from_element(Element::F));
    let before = editor.tracked_snapshot().unwrap();
    let mut plain = editor.clone();
    let mut edits = Edits::new();
    edits.remove_topology(
        (0..4).map(|idx| AtomHandle::Id(AtomId(idx))).collect(),
        vec![],
    );
    edits.remove_atom(AtomHandle::Id(AtomId(5)));
    let expected = TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom,
        index: 5,
        count: 5,
    };
    assert_eq!(
        editor.clone().tracked_apply(edits.clone()).err(),
        Some(expected.clone())
    );
    assert_eq!(plain.transact(edits.clone()), Err(expected.clone()));
    assert_eq!(editor.tracked_transact(edits), Err(expected));
    assert_eq!(plain.tracked_snapshot(), Ok(before.clone()));
    assert_eq!(editor.tracked_snapshot(), Ok(before));
}

#[rstest]
fn test_molecule_editor_tracked_apply() {
    let source = mol_dsl!(r#"{:atoms ["C" "N" "O"]}"#);
    let mut editor = source.edit();
    editor.remove(&[AtomId(0)], &[]);
    editor.add_atom(AtomForm::from_element(Element::F));
    let mut edits = Edits::new();
    edits.remove_atom(AtomHandle::Id(AtomId(0)));
    let added = edits.add_atom(AtomForm::from_element(Element::Cl));
    edits.remove_atom(added);
    edits.add_atom(AtomForm::from_element(Element::N));
    let plain = editor.clone().apply(edits.clone()).unwrap();
    let (tracked, witness) = editor.tracked_apply(edits).unwrap();
    let expected = mol_dsl!(r#"{:atoms ["O" "F" "N"]}"#);
    let local = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(1), AtomId(0)), (AtomId(2), AtomId(1))], 3, 3).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    let session = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(2), AtomId(0))], 3, 3).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    assert_eq!(witness, local);
    assert_eq!(plain.tracked_build(), (expected.clone(), session.clone()));
    assert_eq!(tracked.tracked_build(), (expected, session));
}

#[rstest]
fn test_molecule_editor_tracked_apply_transient() {
    let source = mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#);
    let mut editor = source.edit();
    editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    let mut plain = editor.clone();
    let (editor, applied) = editor.tracked_apply(Edits::new()).unwrap();
    let mut editor = editor;
    let (transaction, transacted) = editor.tracked_transact(Edits::new()).unwrap();
    let rolled_back = transaction.tracked_rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(2),
        Correspondence::identity(2),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    assert_eq!(applied, expected);
    assert_eq!(transacted, expected);
    assert_eq!(rolled_back, expected);
    assert_eq!(
        editor.snapshot(),
        Err(MoleculeIntegrityError::BondsParallel {
            atoms: [AtomId(0), AtomId(1)]
        })
    );
    plain.remove(&[], &[BondId(1)]);
    editor.remove(&[], &[BondId(1)]);
    assert_eq!(editor.tracked_snapshot(), plain.tracked_snapshot());
}

#[rstest]
fn test_molecule_editor_tracked_snapshot_identity(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let editor = source.edit();
    let witness = MoleculeCorrespondence::new(
        Correspondence::new(
            vec![
                (AtomId(0), AtomId(0)),
                (AtomId(1), AtomId(1)),
                (AtomId(2), AtomId(2)),
                (AtomId(3), AtomId(3)),
            ],
            4,
            4,
        )
        .unwrap(),
        Correspondence::new(
            vec![
                (BondId(0), BondId(0)),
                (BondId(1), BondId(1)),
                (BondId(2), BondId(2)),
            ],
            3,
            3,
        )
        .unwrap(),
        Correspondence::new(vec![(DativeBondId(0), DativeBondId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(0))], 1, 1).unwrap(),
        Correspondence::new(vec![(StereoBondId(0), StereoBondId(0))], 1, 1).unwrap(),
    );
    assert_eq!(
        editor.tracked_snapshot(),
        Ok((source.clone(), witness.clone()))
    );
    assert_eq!(
        editor.tracked_snapshot(),
        Ok((source.clone(), witness.clone()))
    );
    assert_eq!(
        editor.clone().try_tracked_build(),
        Ok((source.clone(), witness.clone()))
    );
    assert_eq!(editor.tracked_build(), (source, witness));
}

#[rstest]
fn test_molecule_editor_tracked_build_additions(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let expected = Molecule::from_entries(entries.clone());
    let mut editor = Molecule::new().edit();
    for atom in entries.atoms {
        editor.add_atom(atom);
    }
    for (a, b, bond) in entries.bonds {
        editor.add_bond(a, b, bond);
    }
    for (donors, acceptor, data) in entries.dative {
        editor.add_dative_bond(donors, acceptor, data);
    }
    for (atoms, data) in entries.aromatic {
        editor.add_aromatic_system(atoms, data);
    }
    for (atoms, data) in entries.multicenter {
        editor.add_multicenter_bond(atoms, data);
    }
    for (atoms, data) in entries.noncovalent {
        editor.add_noncovalent_bond(atoms, data);
    }
    for (site, ligands, data) in entries.stereo_atoms {
        editor.add_stereo_atom(site, ligands, data);
    }
    for (site, ligands, data) in entries.stereo_bonds {
        editor.add_stereo_bond(site, ligands, data);
    }
    for constraint in entries.constraints.as_slice() {
        editor.push_constraint(constraint.clone());
    }
    let witness = MoleculeCorrespondence::new(
        Correspondence::new(vec![], 0, 4).unwrap(),
        Correspondence::new(vec![], 0, 3).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::new(vec![], 0, 1).unwrap(),
    );
    assert_eq!(editor.snapshot(), Ok(expected.clone()));
    assert_eq!(
        editor.tracked_snapshot(),
        Ok((expected.clone(), witness.clone()))
    );
    assert_eq!(editor.tracked_build(), (expected, witness));
}

#[rstest]
fn test_molecule_editor_tracked_build_session() {
    let source = mol_dsl!(r#"{:atoms ["C" "N" "O" "F"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]}"#);
    let mut editor = source.edit();
    let first = editor.tracked_remove(&[AtomId(1)], &[]);
    let snapshot = editor.tracked_snapshot().unwrap();
    let added = editor.add_atom(AtomForm::from_element(Element::Cl));
    editor.add_bond(AtomId(2), added, BondForm::from_order(1));
    let second = editor.tracked_remove(&[AtomId(0)], &[]);
    let expected = mol_dsl!(r#"{:atoms ["O" "F" "Cl"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
    let witness = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(2), AtomId(0)), (AtomId(3), AtomId(1))], 4, 3).unwrap(),
        Correspondence::new(vec![(BondId(2), BondId(0))], 3, 2).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    let composed = MoleculeCorrespondence::from(&first)
        .extend_right(EntityKind::Atom, 1)
        .extend_right(EntityKind::Bond, 1)
        .compose(&MoleculeCorrespondence::from(&second))
        .unwrap();
    assert_eq!(composed, witness);
    assert_eq!(
        snapshot,
        (
            mol_dsl!(r#"{:atoms ["C" "O" "F"] :bonds [[1 2 "1"]]}"#),
            MoleculeCorrespondence::from(&first),
        )
    );
    assert_eq!(
        editor.tracked_snapshot(),
        Ok((expected.clone(), witness.clone()))
    );
    assert_eq!(editor.tracked_build(), (expected, witness));
}

#[rstest]
#[case::parallel_bond(MoleculeIntegrityError::BondsParallel { atoms: [AtomId(0), AtomId(1)] })]
fn test_molecule_editor_tracked_snapshot_error(#[case] expected: MoleculeIntegrityError) {
    let source = mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#);
    let mut editor = source.edit();
    let added = editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    assert_eq!(editor.snapshot(), Err(expected.clone()));
    assert_eq!(editor.tracked_snapshot(), Err(expected.clone()));
    assert_eq!(editor.clone().try_build(), Err(expected.clone()));
    assert_eq!(editor.clone().try_tracked_build(), Err(expected));
    editor.remove(&[], &[added]);
    let witness = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 2).unwrap(),
        Correspondence::new(vec![(BondId(0), BondId(0))], 1, 1).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    assert_eq!(editor.tracked_snapshot(), Ok((source, witness)));
}

#[rstest]
#[should_panic(expected = "invalid molecule editor state")]
fn test_molecule_editor_tracked_build_error() {
    let source = mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#);
    let mut editor = source.edit();
    editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    editor.tracked_build();
}

#[rstest]
fn test_molecule_editor_tracked_build_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let mut editor = source.edit();
    let mut edits = Edits::new();
    edits.remove_topology(
        (0..4).map(|idx| AtomHandle::Id(AtomId(idx))).collect(),
        vec![],
    );
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let witness = MoleculeCorrespondence::new(
        Correspondence::new(vec![], 4, 4).unwrap(),
        Correspondence::new(vec![], 3, 3).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::new(vec![], 1, 1).unwrap(),
    );
    assert_eq!(editor.tracked_build(), (source, witness));
}

#[rstest]
fn test_molecule_editor_tracked_snapshot_transaction_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let mut editor = source.edit();
    editor.add_atom(AtomForm::from_element(Element::F));
    let before = editor.tracked_snapshot().unwrap();
    let mut edits = Edits::new();
    edits.remove_topology(
        (0..4).map(|idx| AtomHandle::Id(AtomId(idx))).collect(),
        vec![],
    );
    edits.remove_atom(AtomHandle::Id(AtomId(5)));
    let result = editor.transact(edits);
    assert_eq!(
        result,
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 5,
            count: 5
        })
    );
    assert_eq!(editor.tracked_snapshot(), Ok(before));
}

#[rstest]
fn test_molecule_editor_tracked_build_dative_bonds_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (donors, acceptor, data) = &entries.dative[0];
    edits.remove_dative_bonds(vec![(
        DativeBondHandle::Id(DativeBondId(0)),
        donors
            .iter()
            .copied()
            .chain([*acceptor])
            .map(AtomHandle::Id)
            .collect(),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_aromatic_systems_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (atoms, data) = &entries.aromatic[0];
    edits.remove_aromatic_systems(vec![(
        AromaticSystemHandle::Id(AromaticSystemId(0)),
        atoms.iter().copied().map(AtomHandle::Id).collect(),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::identity(1),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_multicenter_bonds_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (atoms, data) = &entries.multicenter[0];
    edits.remove_multicenter_bonds(vec![(
        MulticenterBondHandle::Id(MulticenterBondId(0)),
        atoms.iter().copied().map(AtomHandle::Id).collect(),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_noncovalent_bonds_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (atoms, data) = &entries.noncovalent[0];
    edits.remove_noncovalent_bonds(vec![(
        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
        atoms.map(AtomHandle::Id),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::identity(1),
        Correspondence::identity(1),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_stereo_atoms_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (site, ligands, data) = &entries.stereo_atoms[0];
    edits.remove_stereo_atoms(vec![(
        StereoAtomHandle::Id(StereoAtomId(0)),
        AtomHandle::Id(*site),
        ligands
            .iter()
            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
            .collect(),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::new(vec![], 1, 1).unwrap(),
        Correspondence::identity(1),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_stereo_bonds_restoration(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let mut edits = Edits::new();
    let (site, ligands, data) = &entries.stereo_bonds[0];
    edits.remove_stereo_bonds(vec![(
        StereoBondHandle::Id(StereoBondId(0)),
        BondHandle::Id(*site),
        ligands
            .iter()
            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
            .collect(),
        data.clone(),
    )]);
    let transaction = editor.transact(edits).unwrap();
    transaction.rollback(&mut editor).unwrap();
    let expected = MoleculeCorrespondence::new(
        Correspondence::identity(4),
        Correspondence::identity(3),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::identity(1),
        Correspondence::new(vec![], 1, 1).unwrap(),
    );
    assert_eq!(editor.tracked_build(), (source, expected));
}

#[rstest]
fn test_molecule_editor_tracked_build_attributes(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries.clone());
    let mut editor = source.edit();
    let witness = editor.tracked_snapshot().unwrap().1;
    entries.atoms[0] = AtomForm::from_element(Element::S);
    entries.bonds[0].2 = BondForm::from_order(2);
    entries.dative[0].2 = DativeBondForm::from_order(2);
    entries.aromatic[0].1 = AromaticSystemForm::from_electrons(vec![0, 1, 2]);
    entries.multicenter[0].1 = MulticenterBondForm::from_electrons(vec![0, 1, 2]);
    entries.noncovalent[0].1 = NoncovalentBondForm::default();
    entries.stereo_atoms[0].2 = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
    entries.stereo_bonds[0].2 = StereoBondForm::new(StereoKind::CisTrans, 0u32);
    *editor.atom_mut(AtomId(0)).attributes = entries.atoms[0].clone();
    *editor.bond_mut(BondId(0)).attributes = entries.bonds[0].2.clone();
    *editor.dative_bond_mut(DativeBondId(0)).attributes = entries.dative[0].2.clone();
    *editor.aromatic_system_mut(AromaticSystemId(0)).attributes = entries.aromatic[0].1.clone();
    *editor.multicenter_bond_mut(MulticenterBondId(0)).attributes =
        entries.multicenter[0].1.clone();
    *editor.noncovalent_bond_mut(NoncovalentBondId(0)).attributes =
        entries.noncovalent[0].1.clone();
    *editor.stereo_atom_mut(StereoAtomId(0)).attributes = entries.stereo_atoms[0].2.clone();
    *editor.stereo_bond_mut(StereoBondId(0)).attributes = entries.stereo_bonds[0].2.clone();
    assert_eq!(
        editor.tracked_build(),
        (Molecule::from_entries(entries), witness)
    );
}

#[rstest]
fn test_molecule_editor_try_tracked_build_allocation(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let source = Molecule::from_entries(entries);
    let editor = source.edit();
    let atoms_ptr = editor.correspondence.atoms().matched_pairs().as_ptr();
    let bonds_ptr = editor.correspondence.bonds().matched_pairs().as_ptr();
    let dative_bonds_ptr = editor
        .correspondence
        .dative_bonds()
        .matched_pairs()
        .as_ptr();
    let aromatic_systems_ptr = editor
        .correspondence
        .aromatic_systems()
        .matched_pairs()
        .as_ptr();
    let multicenter_bonds_ptr = editor
        .correspondence
        .multicenter_bonds()
        .matched_pairs()
        .as_ptr();
    let noncovalent_bonds_ptr = editor
        .correspondence
        .noncovalent_bonds()
        .matched_pairs()
        .as_ptr();
    let stereo_atoms_ptr = editor
        .correspondence
        .stereo_atoms()
        .matched_pairs()
        .as_ptr();
    let stereo_bonds_ptr = editor
        .correspondence
        .stereo_bonds()
        .matched_pairs()
        .as_ptr();
    let (result, witness) = editor.try_tracked_build().unwrap();
    assert_eq!(result, source);
    assert_eq!(witness.atoms().matched_pairs().as_ptr(), atoms_ptr);
    assert_eq!(witness.bonds().matched_pairs().as_ptr(), bonds_ptr);
    assert_eq!(
        witness.dative_bonds().matched_pairs().as_ptr(),
        dative_bonds_ptr
    );
    assert_eq!(
        witness.aromatic_systems().matched_pairs().as_ptr(),
        aromatic_systems_ptr
    );
    assert_eq!(
        witness.multicenter_bonds().matched_pairs().as_ptr(),
        multicenter_bonds_ptr
    );
    assert_eq!(
        witness.noncovalent_bonds().matched_pairs().as_ptr(),
        noncovalent_bonds_ptr
    );
    assert_eq!(
        witness.stereo_atoms().matched_pairs().as_ptr(),
        stereo_atoms_ptr
    );
    assert_eq!(
        witness.stereo_bonds().matched_pairs().as_ptr(),
        stereo_bonds_ptr
    );
}

#[rstest]
fn test_molecule_editor_tracked_remove_dative_bonds_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_dative_bonds(&[]);
    assert_eq!(
        tracked.tracked_remove_dative_bonds(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_aromatic_systems_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_aromatic_systems(&[]);
    assert_eq!(
        tracked.tracked_remove_aromatic_systems(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_multicenter_bonds_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_multicenter_bonds(&[]);
    assert_eq!(
        tracked.tracked_remove_multicenter_bonds(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_noncovalent_bonds_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_noncovalent_bonds(&[]);
    assert_eq!(
        tracked.tracked_remove_noncovalent_bonds(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_stereo_atoms_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_stereo_atoms(&[]);
    assert_eq!(
        tracked.tracked_remove_stereo_atoms(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_stereo_bonds_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove_stereo_bonds(&[]);
    assert_eq!(
        tracked.tracked_remove_stereo_bonds(&[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
fn test_molecule_editor_tracked_remove_empty() {
    let mut plain = Molecule::new().edit();
    let mut tracked = Molecule::new().edit();

    plain.remove(&[], &[]);
    assert_eq!(
        tracked.tracked_remove(&[], &[]),
        MoleculeCompaction::empty()
    );
    assert_eq!(plain.build(), Molecule::new());
    assert_eq!(tracked.build(), Molecule::new());
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![DativeBondId(0)], vec![DativeBondId(0)])]
#[case::middle(vec![DativeBondId(1)], vec![DativeBondId(1)])]
#[case::unsorted_repeated(vec![DativeBondId(2), DativeBondId(0), DativeBondId(2)], vec![DativeBondId(0), DativeBondId(2)])]
#[case::all(vec![DativeBondId(0), DativeBondId(1), DativeBondId(2)], vec![DativeBondId(0), DativeBondId(1), DativeBondId(2)])]
fn test_molecule_editor_tracked_remove_dative_bonds(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<DativeBondId>,
    #[case] removed: Vec<DativeBondId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraintForm::aromatic(false),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&DativeBondId(idx)) {
                entries.dative.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::new(3, removed).unwrap(),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_dative_bonds(&ids);
    let compaction = tracked.tracked_remove_dative_bonds(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(DativeBondId(1))]
#[case::outside(DativeBondId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_dative_bonds_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: DativeBondId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_dative_bonds(&[id]);
    } else {
        editor.remove_dative_bonds(&[id]);
    }
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![AromaticSystemId(0)], vec![AromaticSystemId(0)])]
#[case::middle(vec![AromaticSystemId(1)], vec![AromaticSystemId(1)])]
#[case::unsorted_repeated(vec![AromaticSystemId(2), AromaticSystemId(0), AromaticSystemId(2)], vec![AromaticSystemId(0), AromaticSystemId(2)])]
#[case::all(vec![AromaticSystemId(0), AromaticSystemId(1), AromaticSystemId(2)], vec![AromaticSystemId(0), AromaticSystemId(1), AromaticSystemId(2)])]
fn test_molecule_editor_tracked_remove_aromatic_systems(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<AromaticSystemId>,
    #[case] removed: Vec<AromaticSystemId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::AromaticSystem(
        AromaticSystemId(0),
        AromaticSystemConstraintForm::electron_count(NumForm::Lit(3)),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&AromaticSystemId(idx)) {
                entries.aromatic.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::identity(3),
        Compaction::new(3, removed).unwrap(),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_aromatic_systems(&ids);
    let compaction = tracked.tracked_remove_aromatic_systems(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(AromaticSystemId(1))]
#[case::outside(AromaticSystemId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_aromatic_systems_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: AromaticSystemId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_aromatic_systems(&[id]);
    } else {
        editor.remove_aromatic_systems(&[id]);
    }
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![MulticenterBondId(0)], vec![MulticenterBondId(0)])]
#[case::middle(vec![MulticenterBondId(1)], vec![MulticenterBondId(1)])]
#[case::unsorted_repeated(vec![MulticenterBondId(2), MulticenterBondId(0), MulticenterBondId(2)], vec![MulticenterBondId(0), MulticenterBondId(2)])]
#[case::all(vec![MulticenterBondId(0), MulticenterBondId(1), MulticenterBondId(2)], vec![MulticenterBondId(0), MulticenterBondId(1), MulticenterBondId(2)])]
fn test_molecule_editor_tracked_remove_multicenter_bonds(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<MulticenterBondId>,
    #[case] removed: Vec<MulticenterBondId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::MulticenterBond(
        MulticenterBondId(0),
        MulticenterBondConstraintForm::electron_count(NumForm::Lit(3)),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&MulticenterBondId(idx)) {
                entries.multicenter.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::new(3, removed).unwrap(),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_multicenter_bonds(&ids);
    let compaction = tracked.tracked_remove_multicenter_bonds(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(MulticenterBondId(1))]
#[case::outside(MulticenterBondId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_multicenter_bonds_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: MulticenterBondId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_multicenter_bonds(&[id]);
    } else {
        editor.remove_multicenter_bonds(&[id]);
    }
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![NoncovalentBondId(0)], vec![NoncovalentBondId(0)])]
#[case::middle(vec![NoncovalentBondId(1)], vec![NoncovalentBondId(1)])]
#[case::unsorted_repeated(vec![NoncovalentBondId(2), NoncovalentBondId(0), NoncovalentBondId(2)], vec![NoncovalentBondId(0), NoncovalentBondId(2)])]
#[case::all(vec![NoncovalentBondId(0), NoncovalentBondId(1), NoncovalentBondId(2)], vec![NoncovalentBondId(0), NoncovalentBondId(1), NoncovalentBondId(2)])]
fn test_molecule_editor_tracked_remove_noncovalent_bonds(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<NoncovalentBondId>,
    #[case] removed: Vec<NoncovalentBondId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::NoncovalentBond(
        NoncovalentBondId(0),
        NoncovalentBondConstraintForm::intramolecular(true),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&NoncovalentBondId(idx)) {
                entries.noncovalent.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::new(3, removed).unwrap(),
        Compaction::identity(3),
        Compaction::identity(3),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_noncovalent_bonds(&ids);
    let compaction = tracked.tracked_remove_noncovalent_bonds(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(NoncovalentBondId(1))]
#[case::outside(NoncovalentBondId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_noncovalent_bonds_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: NoncovalentBondId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_noncovalent_bonds(&[id]);
    } else {
        editor.remove_noncovalent_bonds(&[id]);
    }
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![StereoAtomId(0)], vec![StereoAtomId(0)])]
#[case::middle(vec![StereoAtomId(1)], vec![StereoAtomId(1)])]
#[case::unsorted_repeated(vec![StereoAtomId(2), StereoAtomId(0), StereoAtomId(2)], vec![StereoAtomId(0), StereoAtomId(2)])]
#[case::all(vec![StereoAtomId(0), StereoAtomId(1), StereoAtomId(2)], vec![StereoAtomId(0), StereoAtomId(1), StereoAtomId(2)])]
fn test_molecule_editor_tracked_remove_stereo_atoms(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<StereoAtomId>,
    #[case] removed: Vec<StereoAtomId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::StereoAtom(
        StereoAtomId(0),
        StereoKind::Tetrahedral,
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&StereoAtomId(idx)) {
                entries.stereo_atoms.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::new(3, removed).unwrap(),
        Compaction::identity(3),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_stereo_atoms(&ids);
    let compaction = tracked.tracked_remove_stereo_atoms(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(StereoAtomId(1))]
#[case::outside(StereoAtomId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_stereo_atoms_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: StereoAtomId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_stereo_atoms(&[id]);
    } else {
        editor.remove_stereo_atoms(&[id]);
    }
}

#[rstest]
#[case::none(vec![], vec![])]
#[case::first(vec![StereoBondId(0)], vec![StereoBondId(0)])]
#[case::middle(vec![StereoBondId(1)], vec![StereoBondId(1)])]
#[case::unsorted_repeated(vec![StereoBondId(2), StereoBondId(0), StereoBondId(2)], vec![StereoBondId(0), StereoBondId(2)])]
#[case::all(vec![StereoBondId(0), StereoBondId(1), StereoBondId(2)], vec![StereoBondId(0), StereoBondId(1), StereoBondId(2)])]
fn test_molecule_editor_tracked_remove_stereo_bonds(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] ids: Vec<StereoBondId>,
    #[case] removed: Vec<StereoBondId>,
) {
    entries.constraints = Constraints::from_iter([Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    )]);
    let component = Molecule::from_entries(entries.clone());
    let source = Molecule::combine_all([&component, &component, &component]);
    let expected_components = (0..3)
        .map(|idx| {
            let mut entries = entries.clone();
            if removed.contains(&StereoBondId(idx)) {
                entries.stereo_bonds.clear();
                entries.constraints = Constraints::new();
            }
            Molecule::from_entries(entries)
        })
        .collect::<Vec<_>>();
    let expected = Molecule::combine_all(&expected_components);
    let expected_compaction = MoleculeCompaction::new(
        GraphCompaction::new(Compaction::identity(12), Compaction::identity(9)),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::identity(3),
        Compaction::new(3, removed).unwrap(),
    );
    let mut plain = source.edit();
    let mut tracked = source.edit();

    plain.remove_stereo_bonds(&ids);
    let compaction = tracked.tracked_remove_stereo_bonds(&ids);

    assert_eq!(compaction, expected_compaction);
    assert_eq!(plain.build(), expected);
    assert_eq!(
        tracked.tracked_build(),
        (expected, MoleculeCorrespondence::from(&expected_compaction))
    );
}

#[rstest]
#[case::boundary(StereoBondId(1))]
#[case::outside(StereoBondId(2))]
#[should_panic(expected = "removed entities belong to the source table")]
fn test_molecule_editor_tracked_remove_stereo_bonds_error(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
    #[case] id: StereoBondId,
    #[values(false, true)] tracked: bool,
) {
    let mut editor = Molecule::from_entries(entries).edit();
    if tracked {
        editor.tracked_remove_stereo_bonds(&[id]);
    } else {
        editor.remove_stereo_bonds(&[id]);
    }
}

#[rstest]
fn test_molecule_editor_remove_aromatic_systems(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    let result = b.build();
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.atoms().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]
    );
    assert_eq!(
        result.bonds().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![BondId(0), BondId(1), BondId(2)]
    );
}

#[rstest]
fn test_molecule_editor_remove_dative_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.remove_dative_bonds(&[DativeBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_remove_multicenter_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.remove_multicenter_bonds(&[MulticenterBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bonds().ids().collect::<Vec<_>>(),
        Vec::<MulticenterBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_remove_noncovalent_bonds(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_atom_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.atom_mut(AtomId(0)).attributes.element = ElementForm::Lit(Element::N);
    let result = b.build();
    assert_eq!(
        result.atom(AtomId(0)).attributes.element,
        ElementForm::Lit(Element::N)
    );
    assert_eq!(
        molecule.atom(AtomId(0)).attributes.element,
        ElementForm::Lit(Element::C)
    );
}

#[rstest]
fn test_molecule_editor_bond_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.bond_mut(BondId(0)).attributes.order = NumForm::Lit(3);
    let result = b.build();
    assert_eq!(result.bond(BondId(0)).attributes.order, NumForm::Lit(3));
    assert_eq!(molecule.bond(BondId(0)).attributes.order, NumForm::Lit(1));
}

#[rstest]
fn test_molecule_editor_atom_constraint_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.atom_mut(AtomId(0))
        .attributes
        .constraints
        .set(AtomConstraintForm::Degree(NumForm::Lit(2)));
    let result = b.build();
    assert_eq!(
        result.atom(AtomId(0)).attributes.constraints,
        AtomConstraintsForm::from_iter([AtomConstraintForm::Degree(NumForm::Lit(2))])
    );
    assert!(molecule.atom(AtomId(0)).attributes.constraints.is_empty());
}

#[rstest]
fn test_molecule_editor_add_dative_bond(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    let id = b.add_dative_bond(vec![AtomId(1)], AtomId(0), DativeBondForm::from_order(1));
    let result = b.build();
    assert_eq!(id, DativeBondId(1));
    let view = result.dative_bond(id);
    assert_eq!(view.acceptor_id(), AtomId(0));
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![AtomId(1)]);
}

#[rstest]
fn test_molecule_editor_add_multicenter_bond(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    let id = b.add_multicenter_bond(
        vec![AtomId(1), AtomId(2), AtomId(3)],
        MulticenterBondForm::default(),
    );
    let result = b.build();
    assert_eq!(id, MulticenterBondId(1));
    let atoms: Vec<AtomId> = result.multicenter_bond(id).atom_ids().collect();
    assert_eq!(atoms, vec![AtomId(1), AtomId(2), AtomId(3)]);
}

#[rstest]
fn test_molecule_editor_add_noncovalent_bond(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    let id = b.add_noncovalent_bond(
        [AtomId(1), AtomId(2)],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let result = b.build();
    assert_eq!(id, NoncovalentBondId(1));
    let view = result.noncovalent_bond(id);
    assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
}

#[rstest]
fn test_molecule_editor_push_constraint_and_constraints_mut(
    #[from(rich_molecule)] molecule: Molecule,
) {
    let mut b = molecule.edit();
    b.push_constraint(Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    }));
    b.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        }));
    let result = b.build();
    assert_eq!(result.constraints().len(), 2);
}

#[rstest]
fn test_molecule_editor_dative_bond_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.dative_bond_mut(DativeBondId(0))
        .attributes
        .constraints
        .set(DativeBondConstraintForm::ring_membership(
            RingScope::Size(5),
            1,
        ));
    let result = b.build();
    assert!(!result
        .dative_bond(DativeBondId(0))
        .attributes
        .constraints
        .is_empty());
    assert!(molecule
        .dative_bond(DativeBondId(0))
        .attributes
        .constraints
        .is_empty());
}

#[rstest]
fn test_molecule_editor_aromatic_system_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.aromatic_system_mut(AromaticSystemId(0)).attributes.charge = NumForm::Lit(0);
    let result = b.build();
    assert_eq!(
        result
            .aromatic_system(AromaticSystemId(0))
            .attributes
            .charge,
        NumForm::Lit(0)
    );
}

#[rstest]
fn test_molecule_editor_multicenter_bond_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.multicenter_bond_mut(MulticenterBondId(0))
        .attributes
        .electrons = ElectronCountsForm::Lit(vec![1, 1, 0]);
    let result = b.build();
    assert_eq!(
        result
            .multicenter_bond(MulticenterBondId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_editor_noncovalent_bond_mut(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.noncovalent_bond_mut(NoncovalentBondId(0)).attributes.kind =
        NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic);
    let result = b.build();
    assert_eq!(
        result
            .noncovalent_bond(NoncovalentBondId(0))
            .attributes
            .kind,
        NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
    );
}

#[rstest]
fn test_molecule_editor_remove_empty_is_noop(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    b.remove_dative_bonds(&[]);
    b.remove_aromatic_systems(&[]);
    b.remove_multicenter_bonds(&[]);
    b.remove_noncovalent_bonds(&[]);
    let result = b.build();
    assert_eq!(result.dative_bonds().count(), 1);
    assert_eq!(result.aromatic_systems().count(), 1);
    assert_eq!(result.multicenter_bonds().count(), 1);
    assert_eq!(result.noncovalent_bonds().count(), 1);
}

#[rstest]
#[case::hexagon(
    ring(6),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 6,
    },
    vec![(
        vec![
            AtomId(0),
            AtomId(1),
            AtomId(2),
            AtomId(3),
            AtomId(4),
            AtomId(5),
        ],
        vec![
            BondId(0),
            BondId(1),
            BondId(2),
            BondId(3),
            BondId(4),
            BondId(5),
        ],
    )],
)]
#[case::hexagon_cutoff(
    ring(6),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 5,
    },
    vec![],
)]
#[case::chain(
    chain(5),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 10,
    },
    vec![],
)]
#[case::empty(Molecule::default(), RingModel::default(), vec![])]
fn test_molecule_rings(
    #[case] molecule: Molecule,
    #[case] model: RingModel,
    #[case] expected: Vec<(Vec<AtomId>, Vec<BondId>)>,
) {
    let rings = molecule
        .rings(model, RingConfig::default())
        .iter()
        .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
        .collect::<Vec<_>>();
    assert_eq!(rings, expected);
}

#[rstest]
#[case::simple(
    RingSetKind::Simple,
    vec![
        vec![BondId(0), BondId(3), BondId(1)],
        vec![BondId(0), BondId(3), BondId(5), BondId(2)],
        vec![BondId(0), BondId(4), BondId(2)],
        vec![BondId(0), BondId(4), BondId(5), BondId(1)],
        vec![BondId(1), BondId(3), BondId(4), BondId(2)],
        vec![BondId(1), BondId(5), BondId(2)],
        vec![BondId(3), BondId(5), BondId(4)],
    ],
)]
#[case::relevant(
    RingSetKind::Relevant,
    vec![
        vec![BondId(0), BondId(3), BondId(1)],
        vec![BondId(0), BondId(4), BondId(2)],
        vec![BondId(1), BondId(5), BondId(2)],
        vec![BondId(3), BondId(5), BondId(4)],
    ],
)]
fn test_molecule_rings_kind(#[case] kind: RingSetKind, #[case] mut expected: Vec<Vec<BondId>>) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        ..Default::default()
    });
    let mut actual = molecule
        .rings(
            RingModel {
                kind,
                max_ring_size: 4,
            },
            RingConfig::default(),
        )
        .iter()
        .map(|ring| ring.bonds().to_vec())
        .collect::<Vec<_>>();
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected);
}

#[rstest]
fn test_molecule_editor_add_and_remove(#[from(rich_molecule)] molecule: Molecule) {
    let mut b = molecule.edit();
    let new_a = b.add_atom(AtomForm::from_element(Element::Br));
    b.add_bond(AtomId(0), new_a, BondForm::from_order(1));
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    b.remove(&[AtomId(3)], &[BondId(2)]);
    let result = b.build();
    let atoms: Vec<Element> = result
        .atoms()
        .iter()
        .map(|v| match v.attributes.element {
            ElementForm::Lit(e) => e,
            _ => panic!("non-ground element in editor result"),
        })
        .collect();
    assert_eq!(atoms, vec![Element::C, Element::C, Element::N, Element::Br]);
    let bonds: Vec<(AtomId, AtomId, NumForm)> = result
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.attributes.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), NumForm::Lit(1)),
            (AtomId(1), AtomId(2), NumForm::Lit(2)),
            (AtomId(0), AtomId(3), NumForm::Lit(1)),
        ]
    );
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
#[case::donor_below_acceptor(AtomId(0), AtomId(1))]
#[case::donor_above_acceptor(AtomId(1), AtomId(0))]
fn test_molecule_dative_acceptor_donor(#[case] donor: AtomId, #[case] acceptor: AtomId) {
    let atoms = vec![ground_atom(), ground_atom()];
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms,
        dative: vec![(vec![donor], acceptor, DativeBondForm::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let view = molecule.dative_bond(DativeBondId(0));
    assert_eq!(view.acceptor_id(), acceptor);
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![donor]);
}

#[rstest]
fn test_molecule_eq_canonical_across_bond_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let bond = BondForm {
        order: NumForm::Lit(1),
        charge: NumForm::Lit(0),
        unpaired_electrons: UnpairedElectronsForm::closed_shell(),
        constraints: BondConstraintsForm::new(),
    };
    let forward = Molecule::from_entries(MoleculeEntries {
        atoms: atoms_a,
        bonds: vec![(AtomId(0), AtomId(1), bond.clone())],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let reverse = Molecule::from_entries(MoleculeEntries {
        atoms: atoms_b,
        bonds: vec![(AtomId(1), AtomId(0), bond)],
        constraints: Constraints::new(),
        ..Default::default()
    });
    assert_eq!(forward, reverse);
}

#[rstest]
fn test_molecule_eq_canonical_across_dative_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let forward = Molecule::from_entries(MoleculeEntries {
        atoms: atoms_a,
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let reverse = Molecule::from_entries(MoleculeEntries {
        atoms: atoms_b,
        dative: vec![(vec![AtomId(1)], AtomId(0), DativeBondForm::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    assert_ne!(
        forward, reverse,
        "acceptor identity is part of dative bond; swapping donor/acceptor should differ"
    );
}

#[rstest]
fn test_molecule_raw_graph(#[from(rich_molecule)] molecule: Molecule) {
    let g = molecule.raw_graph();
    assert_eq!(g.node_count(), 4);
    assert_eq!(g.edge_count(), 3);
    assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
#[case::disjoint(
    HashSet::from([AtomId(3)]),
    None,
)]
fn test_aromatic_system_views_of_id(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(molecule.aromatic_systems().of_id(atoms), expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
fn test_multicenter_bond_views_of_id(
    #[from(rich_molecule)] molecule: Molecule,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(molecule.multicenter_bonds().of_id(atoms), expected);
}

#[rstest]
#[case::ring_4(
    ring(4),
    vec![
        vec![BondId(0), BondId(2)],
        vec![BondId(1), BondId(3)],
    ],
)]
fn test_graph_view_visit_maximum_matchings(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Vec<BondId>>,
) {
    let mut matchings: Vec<Vec<BondId>> = Vec::new();
    let flow: ControlFlow<()> = molecule.graph().visit_maximum_matchings(
        MatchingEnumerationAlgorithm::BranchAndBound,
        |matching| {
            let mut bonds: Vec<BondId> = matching.bonds().collect();
            bonds.sort_unstable();
            matchings.push(bonds);
            ControlFlow::Continue(())
        },
    );
    assert_eq!(flow, ControlFlow::Continue(()));
    matchings.sort_unstable();
    assert_eq!(matchings, expected);
}

#[rstest]
fn test_molecule_enumerate_maximum_matchings() {
    let molecule = ring(4);
    let mut ms: Vec<Vec<(AtomId, AtomId)>> = molecule
        .graph()
        .enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound)
        .into_iter()
        .map(|m| {
            let mut pairs: Vec<_> = (0..molecule.atoms().count())
                .map(AtomId::from)
                .filter_map(|a| m.mate(a).filter(|b| a < *b).map(|b| (a, b)))
                .collect();
            pairs.sort_unstable();
            pairs
        })
        .collect();
    ms.sort_unstable();
    assert_eq!(
        ms,
        vec![
            vec![(AtomId(0), AtomId(1)), (AtomId(2), AtomId(3))],
            vec![(AtomId(0), AtomId(3)), (AtomId(1), AtomId(2))],
        ]
    );
}

#[rstest]
fn test_molecule_index_atom(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(
        molecule.atom(AtomId(2)).attributes.element,
        ElementForm::Lit(Element::N)
    );
}

#[rstest]
fn test_molecule_index_bond(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(molecule.bond(BondId(1)).attributes.order, NumForm::Lit(2));
}

#[rstest]
fn test_molecule_index_dative_bond(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(
        molecule.dative_bond(DativeBondId(0)).attributes.order,
        NumForm::Lit(1)
    );
}

#[rstest]
fn test_molecule_index_aromatic_system(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(
        molecule
            .aromatic_system(AromaticSystemId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Undetermined
    );
}

#[rstest]
fn test_molecule_index_multicenter_bond(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(
        molecule
            .multicenter_bond(MulticenterBondId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Undetermined
    );
}

#[rstest]
fn test_molecule_index_noncovalent_bond(#[from(rich_molecule)] molecule: Molecule) {
    assert_eq!(
        molecule
            .noncovalent_bond(NoncovalentBondId(0))
            .attributes
            .kind,
        NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)
    );
}

#[rstest]
fn test_molecule_modify_atoms(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule.modify_atoms(|mut a| {
        a.charge = NumForm::Lit(1);
        a
    });
    let charges: Vec<NumForm> = molecule
        .atoms()
        .iter()
        .map(|v| v.attributes.charge.clone())
        .collect();
    assert_eq!(
        charges,
        vec![
            NumForm::Lit(1),
            NumForm::Lit(1),
            NumForm::Lit(1),
            NumForm::Lit(1),
        ]
    );
}

#[rstest]
fn test_molecule_modify_bonds(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule.modify_bonds(|mut b| {
        b.order = NumForm::Lit(1);
        b
    });
    let orders: Vec<NumForm> = molecule
        .bonds()
        .iter()
        .map(|v| v.attributes.order.clone())
        .collect();
    assert_eq!(
        orders,
        vec![NumForm::Lit(1), NumForm::Lit(1), NumForm::Lit(1)]
    );
}

#[rstest]
fn test_molecule_dative_bond_mut(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule
        .dative_bond_mut(DativeBondId(0))
        .attributes
        .constraints
        .set(DativeBondConstraintForm::ring_membership(
            RingScope::Size(6),
            1,
        ));
    assert_eq!(
        molecule.dative_bond(DativeBondId(0)).attributes.constraints,
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(
            RingScope::Size(6),
            1
        )])
    );
}

#[rstest]
fn test_molecule_aromatic_system_mut(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule
        .aromatic_system_mut(AromaticSystemId(0))
        .attributes
        .electrons = ElectronCountsForm::Lit(vec![1; 3]);
    assert_eq!(
        molecule
            .aromatic_system(AromaticSystemId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Lit(vec![1, 1, 1]),
    );
}

#[rstest]
fn test_molecule_modify_aromatic_systems(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule.modify_aromatic_systems(|mut a| {
        a.electrons = ElectronCountsForm::Lit(vec![1; 3]);
        a
    });
    let electrons: Vec<ElectronCountsForm> = molecule
        .aromatic_systems()
        .iter()
        .map(|v| v.attributes.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsForm::Lit(vec![1; 3])]);
}

#[rstest]
fn test_molecule_try_modify_aromatic_system(#[from(rich_molecule)] mut molecule: Molecule) {
    assert_eq!(
        molecule.try_modify_aromatic_system(AromaticSystemId(0), |form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 1, 0]);
        }),
        Ok(()),
    );
    assert_eq!(
        molecule
            .aromatic_system(AromaticSystemId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Lit(vec![2, 1, 0]),
    );
}

#[rstest]
#[case::invalid_reference(
    AromaticSystemId(1),
    ElectronCountsForm::Lit(vec![2, 1, 0]),
    MoleculeIntegrityError::InvalidReference {
        entity: Entity::AromaticSystem(AromaticSystemId(1)),
    },
)]
#[case::electron_count_length(
    AromaticSystemId(0),
    ElectronCountsForm::Lit(vec![2, 1]),
    MoleculeIntegrityError::ElectronCountLengthMismatch {
        entity: Entity::AromaticSystem(AromaticSystemId(0)),
        participants: 3,
        electron_counts: 2,
    },
)]
fn test_molecule_try_modify_aromatic_system_error(
    #[from(rich_molecule)] mut molecule: Molecule,
    #[case] id: AromaticSystemId,
    #[case] electrons: ElectronCountsForm,
    #[case] expected: MoleculeIntegrityError,
) {
    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_aromatic_system(id, |form| form.electrons = electrons),
        Err(expected),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_try_modify_aromatic_systems(#[from(rich_molecule)] mut molecule: Molecule) {
    assert_eq!(
        molecule.try_modify_aromatic_systems(|form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 1, 0]);
        }),
        Ok(()),
    );
    assert_eq!(
        molecule
            .aromatic_systems()
            .iter()
            .map(|view| view.attributes.electrons.clone())
            .collect::<Vec<_>>(),
        vec![ElectronCountsForm::Lit(vec![2, 1, 0])],
    );
}

#[rstest]
fn test_molecule_try_modify_aromatic_systems_error(#[from(rich_molecule)] mut molecule: Molecule) {
    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_aromatic_systems(|form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 1]);
        }),
        Err(MoleculeIntegrityError::ElectronCountLengthMismatch {
            entity: Entity::AromaticSystem(AromaticSystemId(0)),
            participants: 3,
            electron_counts: 2,
        }),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_multicenter_bond_mut(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule
        .multicenter_bond_mut(MulticenterBondId(0))
        .attributes
        .electrons = ElectronCountsForm::Lit(vec![1, 1, 0]);
    assert_eq!(
        molecule
            .multicenter_bond(MulticenterBondId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_modify_multicenter_bonds(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule.modify_multicenter_bonds(|mut m| {
        m.electrons = ElectronCountsForm::Lit(vec![1, 1, 0]);
        m
    });
    let electrons: Vec<ElectronCountsForm> = molecule
        .multicenter_bonds()
        .iter()
        .map(|v| v.attributes.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsForm::Lit(vec![1, 1, 0])],);
}

#[rstest]
fn test_molecule_try_modify_multicenter_bond(#[from(rich_molecule)] mut molecule: Molecule) {
    assert_eq!(
        molecule.try_modify_multicenter_bond(MulticenterBondId(0), |form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 0, 0]);
        }),
        Ok(()),
    );
    assert_eq!(
        molecule
            .multicenter_bond(MulticenterBondId(0))
            .attributes
            .electrons,
        ElectronCountsForm::Lit(vec![2, 0, 0]),
    );
}

#[rstest]
#[case::invalid_reference(
    MulticenterBondId(1),
    ElectronCountsForm::Lit(vec![2, 0, 0]),
    MoleculeIntegrityError::InvalidReference {
        entity: Entity::MulticenterBond(MulticenterBondId(1)),
    },
)]
#[case::electron_count_length(
    MulticenterBondId(0),
    ElectronCountsForm::Lit(vec![2, 0]),
    MoleculeIntegrityError::ElectronCountLengthMismatch {
        entity: Entity::MulticenterBond(MulticenterBondId(0)),
        participants: 3,
        electron_counts: 2,
    },
)]
fn test_molecule_try_modify_multicenter_bond_error(
    #[from(rich_molecule)] mut molecule: Molecule,
    #[case] id: MulticenterBondId,
    #[case] electrons: ElectronCountsForm,
    #[case] expected: MoleculeIntegrityError,
) {
    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_multicenter_bond(id, |form| form.electrons = electrons),
        Err(expected),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_try_modify_multicenter_bonds(#[from(rich_molecule)] mut molecule: Molecule) {
    assert_eq!(
        molecule.try_modify_multicenter_bonds(|form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 0, 0]);
        }),
        Ok(()),
    );
    assert_eq!(
        molecule
            .multicenter_bonds()
            .iter()
            .map(|view| view.attributes.electrons.clone())
            .collect::<Vec<_>>(),
        vec![ElectronCountsForm::Lit(vec![2, 0, 0])],
    );
}

#[rstest]
fn test_molecule_try_modify_multicenter_bonds_error(#[from(rich_molecule)] mut molecule: Molecule) {
    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_multicenter_bonds(|form| {
            form.electrons = ElectronCountsForm::Lit(vec![2, 0]);
        }),
        Err(MoleculeIntegrityError::ElectronCountLengthMismatch {
            entity: Entity::MulticenterBond(MulticenterBondId(0)),
            participants: 3,
            electron_counts: 2,
        }),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_noncovalent_bond_mut(#[from(rich_molecule)] mut molecule: Molecule) {
    molecule
        .noncovalent_bond_mut(NoncovalentBondId(0))
        .attributes
        .kind = NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic);
    assert_eq!(
        molecule
            .noncovalent_bond(NoncovalentBondId(0))
            .attributes
            .kind,
        NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)
    );
}

#[rstest]
fn test_molecule_try_modify_stereo_atom(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let mut molecule = Molecule::from_entries(entries);
    molecule
        .try_modify_stereo_atom(StereoAtomId(0), |form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        })
        .expect("tetrahedral configuration satisfies atom-site integrity");
    assert_eq!(
        molecule
            .stereo_atom(StereoAtomId(0))
            .attributes
            .configuration,
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
    );

    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_stereo_atom(StereoAtomId(0), |form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0));
        }),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            kind: StereoKind::CisTrans,
        }),
    );
    assert_eq!(molecule, before);
    assert_eq!(
        molecule.try_modify_stereo_atom(StereoAtomId(1), |_| {}),
        Err(MoleculeIntegrityError::InvalidReference {
            entity: Entity::StereoAtom(StereoAtomId(1)),
        }),
    );
}

#[rstest]
fn test_molecule_try_modify_stereo_atoms(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let mut molecule = Molecule::from_entries(entries);
    molecule
        .try_modify_stereo_atoms(|form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        })
        .expect("tetrahedral configurations satisfy atom-site integrity");

    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_stereo_atoms(|form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0));
        }),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            kind: StereoKind::CisTrans,
        }),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_try_modify_stereo_bond(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let mut molecule = Molecule::from_entries(entries);
    molecule
        .try_modify_stereo_bond(StereoBondId(0), |form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0));
        })
        .expect("cis/trans configuration satisfies bond-site integrity");
    assert_eq!(
        molecule
            .stereo_bond(StereoBondId(0))
            .attributes
            .configuration,
        StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
    );

    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_stereo_bond(StereoBondId(0), |form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        }),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch {
            entity: Entity::StereoBond(StereoBondId(0)),
            kind: StereoKind::Tetrahedral,
        }),
    );
    assert_eq!(molecule, before);
    assert_eq!(
        molecule.try_modify_stereo_bond(StereoBondId(1), |_| {}),
        Err(MoleculeIntegrityError::InvalidReference {
            entity: Entity::StereoBond(StereoBondId(1)),
        }),
    );
}

#[rstest]
fn test_molecule_try_modify_stereo_bonds(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let mut molecule = Molecule::from_entries(entries);
    molecule
        .try_modify_stereo_bonds(|form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0));
        })
        .expect("cis/trans configurations satisfy bond-site integrity");

    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_stereo_bonds(|form| {
            form.configuration =
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        }),
        Err(MoleculeIntegrityError::StereoKindSiteMismatch {
            entity: Entity::StereoBond(StereoBondId(0)),
            kind: StereoKind::Tetrahedral,
        }),
    );
    assert_eq!(molecule, before);
}

#[rstest]
fn test_molecule_try_modify_constraints(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let mut molecule = Molecule::from_entries(entries);
    let valid = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(NumForm::Lit(4)));
    molecule
        .try_modify_constraints(|constraints| constraints.push(valid.clone()))
        .expect("constraint references an available atom");
    assert!(molecule.constraints().iter().any(|entry| entry == &valid));

    let before = molecule.clone();
    assert_eq!(
        molecule.try_modify_constraints(|constraints| {
            constraints.push(Constraint::Atom(
                AtomId(4),
                AtomConstraintForm::degree(NumForm::Lit(4)),
            ));
        }),
        Err(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(4)),
        }),
    );
    assert_eq!(molecule, before);
}

// -- lift_constraints / inline_constraints ---------------------

/// Set-equality assertion for constraint vecs: order is unspecified, so the
/// test compares as multisets via sort+eq.
fn assert_same_constraints(a: &Constraints, b: &Constraints) {
    let mut x: Vec<&Constraint> = a.iter().collect();
    let mut y: Vec<&Constraint> = b.iter().collect();
    x.sort_by_key(|c| format!("{c:?}"));
    y.sort_by_key(|c| format!("{c:?}"));
    assert_eq!(x, y);
}

#[rstest]
fn test_molecule_lift_constraints_empty() {
    let mut molecule = Molecule::default();
    molecule.lift_constraints();
    assert!(molecule.constraints().is_empty());
}

#[rstest]
fn test_molecule_lift_constraints_drains_inline_stores(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    molecule
        .atom_mut(AtomId(0))
        .attributes
        .constraints
        .set(AtomConstraintForm::Valence(NumForm::Lit(4)));
    molecule
        .atom_mut(AtomId(2))
        .attributes
        .constraints
        .set(AtomConstraintForm::Degree(NumForm::Lit(3)));
    molecule
        .bond_mut(BondId(0))
        .attributes
        .constraints
        .set(BondConstraintForm::Aromatic(BooleanForm::Lit(true)));
    molecule
        .dative_bond_mut(DativeBondId(0))
        .attributes
        .constraints
        .set(DativeBondConstraintForm::ring_membership(
            RingScope::All,
            NumForm::Lit(1),
        ));

    molecule.lift_constraints();

    assert!(molecule.atom(AtomId(0)).attributes.constraints.is_empty());
    assert!(molecule.atom(AtomId(2)).attributes.constraints.is_empty());
    assert!(molecule.bond(BondId(0)).attributes.constraints.is_empty());
    assert!(molecule
        .dative_bond(DativeBondId(0))
        .attributes
        .constraints
        .is_empty());

    let mut expected = Constraints::new();
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));
    expected.push(Constraint::Atom(
        AtomId(2),
        AtomConstraintForm::Degree(NumForm::Lit(3)),
    ));
    expected.push(Constraint::Bond(
        BondId(0),
        BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
    ));
    expected.push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1)),
    ));
    assert_same_constraints(molecule.constraints(), &expected);
}

#[rstest]
fn test_molecule_lift_constraints_appends_to_existing(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    let prior = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    molecule.constraints_mut().push(prior.clone());
    molecule
        .atom_mut(AtomId(0))
        .attributes
        .constraints
        .set(AtomConstraintForm::Valence(NumForm::Lit(4)));

    molecule.lift_constraints();

    let mut expected = Constraints::new();
    expected.push(prior);
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));
    assert_same_constraints(molecule.constraints(), &expected);
}

#[rstest]
fn test_molecule_inline_constraints_drains_top_level_leaves(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    molecule.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));
    molecule.constraints_mut().push(Constraint::Bond(
        BondId(0),
        BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
    ));
    molecule.constraints_mut().push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1),
    ));

    molecule.inline_constraints().unwrap();

    assert!(molecule.constraints().is_empty());
    assert_eq!(
        molecule.atom(AtomId(0)).attributes.constraints,
        AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::Lit(4))])
    );
    assert_eq!(
        molecule.bond(BondId(0)).attributes.constraints,
        BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))])
    );
    assert_eq!(
        molecule.dative_bond(DativeBondId(0)).attributes.constraints,
        DativeBondConstraintsForm::from_iter([DativeBondConstraintForm::ring_membership(
            RingScope::Size(5),
            1
        )])
    );
}

#[rstest]
fn test_molecule_inline_constraints_last_wins_on_collision(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    molecule.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(3)),
    ));
    molecule.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));

    molecule.inline_constraints().unwrap();

    // Only one Valence survives; with two competing inserts of the same kind,
    // exactly one wins (which one is unspecified). Verify count and kind.
    assert_eq!(molecule.atom(AtomId(0)).attributes.constraints.len(), 1);
    let v = molecule
        .atom(AtomId(0))
        .attributes
        .constraints
        .iter()
        .next()
        .unwrap()
        .clone();
    assert!(matches!(v, AtomConstraintForm::Valence(_)));
}

#[rstest]
fn test_molecule_inline_constraints_skips_combinator_nested(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    let leaf = Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(4)));
    let nested = Constraint::And(vec![
        leaf.clone(),
        Constraint::Bond(
            BondId(0),
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
        ),
    ]);
    molecule.constraints_mut().push(nested.clone());

    molecule.inline_constraints().unwrap();

    let mut expected = Constraints::new();
    expected.push(nested);
    assert_same_constraints(molecule.constraints(), &expected);
    assert!(molecule.atom(AtomId(0)).attributes.constraints.is_empty());
    assert!(molecule.bond(BondId(0)).attributes.constraints.is_empty());
}

#[rstest]
fn test_molecule_inline_constraints_skips_relational_and_molecule(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    let rel = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    let mol = Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    });
    molecule.constraints_mut().push(rel.clone());
    molecule.constraints_mut().push(mol.clone());
    molecule.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(4)),
    ));

    molecule.inline_constraints().unwrap();

    let mut expected = Constraints::new();
    expected.push(rel);
    expected.push(mol);
    assert_same_constraints(molecule.constraints(), &expected);
    assert_eq!(
        molecule.atom(AtomId(0)).attributes.constraints,
        AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::Lit(4))])
    );
}

#[rstest]
fn test_molecule_lift_then_inline_roundtrips_inline_state(
    #[from(rich_molecule)] mut molecule: Molecule,
) {
    molecule
        .atom_mut(AtomId(0))
        .attributes
        .constraints
        .set(AtomConstraintForm::Valence(NumForm::Lit(4)));
    molecule
        .atom_mut(AtomId(0))
        .attributes
        .constraints
        .set(AtomConstraintForm::Degree(NumForm::Lit(3)));
    molecule
        .bond_mut(BondId(0))
        .attributes
        .constraints
        .set(BondConstraintForm::Aromatic(BooleanForm::Lit(true)));
    molecule
        .dative_bond_mut(DativeBondId(0))
        .attributes
        .constraints
        .set(DativeBondConstraintForm::ring_membership(
            RingScope::All,
            NumForm::Lit(1),
        ));

    let original = molecule.clone();

    molecule.lift_constraints();
    assert!(molecule.atom(AtomId(0)).attributes.constraints.is_empty());
    molecule.inline_constraints().unwrap();

    assert_eq!(molecule, original);
}

#[rstest]
#[case::empty(Vec::new(), Molecule::new())]
#[case::singleton(
    vec![Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    })],
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    }),
)]
#[case::multiple(
    vec![
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        Molecule::new(),
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            ..Default::default()
        }),
    ],
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(2))],
        ..Default::default()
    }),
)]
fn test_molecule_combine_all(#[case] molecules: Vec<Molecule>, #[case] expected: Molecule) {
    let combined = Molecule::combine_all(&molecules);

    assert_eq!(combined, expected);
}

#[rstest]
fn test_molecule_combine() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
        ..Default::default()
    });
    let union = left.combine(&right);

    assert_eq!(union.atoms().count(), 4);
    assert_eq!(union.bonds().count(), 2);
    assert_eq!(union.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    assert_eq!(union.bond(BondId(1)).atom_ids(), [AtomId(2), AtomId(3)]);
    assert_eq!(union.bond(BondId(1)).attributes, &BondForm::from_order(2));
}

#[rstest]
fn test_molecule_combine_from() {
    let mut left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    left.combine_from(&right);

    assert_eq!(left.atoms().count(), 3);
    assert_eq!(left.bond(BondId(0)).atom_ids(), [AtomId(1), AtomId(2)]);
}

#[rstest]
fn test_molecule_combine_from_storage() {
    let mut left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    Arc::get_mut(&mut left.atoms).unwrap().reserve(2);
    Arc::get_mut(&mut left.bonds).unwrap().reserve(1);
    let atom_storage = left.atoms.as_ptr();
    let bond_storage = left.bonds.as_ptr();
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::F),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
        ..Default::default()
    });

    left.combine_from(&right);

    assert_eq!(left.atoms.as_ptr(), atom_storage);
    assert_eq!(left.bonds.as_ptr(), bond_storage);
}

#[rstest]
fn test_molecule_combine_overlay() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1)],
            AromaticSystemForm::from_electrons(vec![1, 1]),
        )],
        ..Default::default()
    });
    let union = left.combine(&right);

    assert_eq!(union.aromatic_systems().count(), 1);
    // right's overlay over its atoms [0, 1] shifts by left's one atom
    assert_eq!(
        union
            .aromatic_system(AromaticSystemId(0))
            .atom_ids()
            .collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2)]
    );
}

#[rstest]
fn test_molecule_combine_from_stereo() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
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
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(4)]),
            },
        )]),
        ..Default::default()
    });
    let mut union = left;
    union.combine_from(&right);
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(5)]),
            },
        )]),
        ..Default::default()
    });
    assert_eq!(union, expected);
}

#[rstest]
fn test_molecule_combine_constraint() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
                sum: NumForm::Lit(0),
            },
        )),
        ..Default::default()
    });
    let union = left.combine(&right);

    // right's constraint over atoms [0, 1] is remapped to [1, 2] in the union
    let expected = Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(1), AtomId(2)]),
        sum: NumForm::Lit(0),
    });
    assert_eq!(
        union.constraints.iter().collect::<Vec<_>>(),
        vec![&expected]
    );
}

#[rstest]
fn test_molecule_split_empty() {
    let molecule = Molecule::new();

    assert!(molecule.split().is_empty());
    assert!(molecule.tracked_split().is_empty());
}

#[rstest]
fn test_molecule_split() {
    // two disconnected bonds → two components
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        ..Default::default()
    });
    let components = mol.tracked_split();
    assert_eq!(
        mol.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );

    assert_eq!(components.len(), 2);
    let (first, first_corr) = &components[0];
    assert_eq!(first.atoms().count(), 2);
    assert_eq!(first.bond(BondId(0)).attributes, &BondForm::from_order(1));
    assert_eq!(first_corr.atoms().left_of(AtomId(0)), Some(AtomId(0)));
    assert_eq!(first_corr.atoms().left_of(AtomId(1)), Some(AtomId(1)));
    let (second, second_corr) = &components[1];
    assert_eq!(second.bond(BondId(0)).attributes, &BondForm::from_order(2));
    assert_eq!(second_corr.atoms().left_of(AtomId(0)), Some(AtomId(2)));
    assert_eq!(second_corr.atoms().left_of(AtomId(1)), Some(AtomId(3)));
}

#[rstest]
fn test_molecule_split_interleaved() {
    let input = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(2)),
        ],
        constraints: Constraints::from_iter([
            Constraint::Molecule(MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(3), AtomId(0)]),
            }),
            Constraint::Molecule(MoleculeConstraint::BondOrderSum {
                bonds: Some(vec![BondId(1)]),
                sum: NumForm::Lit(2),
            }),
        ]),
        ..Default::default()
    });
    let first = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(0)]),
            },
        )]),
        ..Default::default()
    });
    let second = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::BondOrderSum {
                bonds: Some(vec![BondId(0)]),
                sum: NumForm::Lit(2),
            },
        )]),
        ..Default::default()
    });
    let third = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let first_correspondence = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(3))], 2, 5).unwrap(),
        Correspondence::new(vec![(BondId(0), BondId(0))], 1, 2).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    let second_correspondence = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(1)), (AtomId(1), AtomId(4))], 2, 5).unwrap(),
        Correspondence::new(vec![(BondId(0), BondId(1))], 1, 2).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    let third_correspondence = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(2))], 1, 5).unwrap(),
        Correspondence::new(vec![], 0, 2).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
    );
    assert_eq!(
        input.split(),
        vec![first.clone(), second.clone(), third.clone()]
    );
    assert_eq!(
        input.tracked_split(),
        vec![
            (first, first_correspondence.reverse()),
            (second, second_correspondence.reverse()),
            (third, third_correspondence.reverse()),
        ]
    );
}

#[rstest]
fn test_molecule_split_overlay_binds() {
    // two disconnected bonds, but an aromatic system over {1, 2} keeps all four atoms in one component
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![1, 1]),
        )],
        ..Default::default()
    });
    let components = mol.tracked_split();
    assert_eq!(
        mol.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );

    assert_eq!(components.len(), 1);
    assert_eq!(components[0].0.atoms().count(), 4);
}

#[rstest]
fn test_molecule_combine_split_roundtrip() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
        ..Default::default()
    });
    let union = left.combine(&right);
    let components = union.tracked_split();
    assert_eq!(
        union.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );

    assert_eq!(components.len(), 2);
    assert_eq!(components[0].0, left);
    assert_eq!(components[1].0, right);
}

#[rstest]
fn test_molecule_split_stereo() {
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: (0..7).map(|_| AtomForm::from_element(Element::C)).collect(),
        bonds: iter::once((AtomId(0), AtomId(1), BondForm::from_order(1)))
            .chain((3..=6).map(|id| (AtomId(2), AtomId(id), BondForm::from_order(1))))
            .collect(),
        stereo_atoms: vec![(
            AtomId(2),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(2), AtomId(6)]),
            },
        )]),
        ..Default::default()
    });
    let components = mol.tracked_split();
    assert_eq!(
        mol.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );
    let lone = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let bound = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::from_iter([Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(4)]),
            },
        )]),
        ..Default::default()
    });
    let lone_correspondence = MoleculeCorrespondence::new(
        Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 7).unwrap(),
        Correspondence::new(vec![(BondId(0), BondId(0))], 1, 5).unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::new(vec![], 0, 1).unwrap(),
        Correspondence::empty(),
    );
    let bound_correspondence = MoleculeCorrespondence::new(
        Correspondence::new(
            vec![
                (AtomId(0), AtomId(2)),
                (AtomId(1), AtomId(3)),
                (AtomId(2), AtomId(4)),
                (AtomId(3), AtomId(5)),
                (AtomId(4), AtomId(6)),
            ],
            5,
            7,
        )
        .unwrap(),
        Correspondence::new(
            vec![
                (BondId(0), BondId(1)),
                (BondId(1), BondId(2)),
                (BondId(2), BondId(3)),
                (BondId(3), BondId(4)),
            ],
            4,
            5,
        )
        .unwrap(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::empty(),
        Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(0))], 1, 1).unwrap(),
        Correspondence::empty(),
    );
    assert_eq!(
        components,
        vec![
            (lone, lone_correspondence.reverse()),
            (bound, bound_correspondence.reverse())
        ]
    );
}

#[rstest]
fn test_molecule_split_constraint_binds() {
    // two disconnected bonds, but a ChargeSum over {1, 2} binds all four atoms into one component
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
                sum: NumForm::Lit(0),
            },
        )),
        ..Default::default()
    });
    let components = mol.tracked_split();
    assert_eq!(
        mol.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );

    assert_eq!(components.len(), 1);
    assert_eq!(
        components[0].0.constraints.iter().collect::<Vec<_>>(),
        vec![&Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(1), AtomId(2)]),
            sum: NumForm::Lit(0),
        })]
    );
}

#[rstest]
fn test_molecule_split_constraint_routed() {
    // a constraint over the second component's atoms routes there, remapped to compact ids
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(2), AtomId(3)]),
                sum: NumForm::Lit(0),
            },
        )),
        ..Default::default()
    });
    let components = mol.tracked_split();
    assert_eq!(
        mol.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );

    assert_eq!(components.len(), 2);
    assert!(components[0].0.constraints.is_empty());
    assert_eq!(
        components[1].0.constraints.iter().collect::<Vec<_>>(),
        vec![&Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: NumForm::Lit(0),
        })]
    );
}

#[rstest]
fn test_molecule_split_constraint_entity_kinds(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let mut left_entries = entries.clone();
    left_entries.constraints = Constraints::new();
    let left = Molecule::from_entries(left_entries);

    let mut right_entries = entries;
    right_entries.constraints = vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(NumForm::Lit(4))),
        Constraint::Bond(BondId(0), BondConstraintForm::aromatic(false)),
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(false)),
        Constraint::AromaticSystem(
            AromaticSystemId(0),
            AromaticSystemConstraintForm::electron_count(NumForm::Lit(6)),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(0),
            MulticenterBondConstraintForm::electron_count(NumForm::Lit(2)),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(0),
            NoncovalentBondConstraintForm::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::Relational(RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(0),
            parallel: BondId(1),
        }),
        Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(3)]),
            sum: NumForm::Lit(1),
        }),
        Constraint::Molecule(MoleculeConstraint::BondOrderSum {
            bonds: Some(vec![BondId(0), BondId(2)]),
            sum: NumForm::Lit(2),
        }),
    ]
    .into();
    let right = Molecule::from_entries(right_entries);

    let combined = left.combine(&right);
    let components = combined.tracked_split();
    assert_eq!(
        combined.split(),
        components
            .iter()
            .map(|(component, _)| component.clone())
            .collect::<Vec<_>>()
    );
    let expected_left = MoleculeCorrespondence::new(
        Correspondence::from_images(&[AtomId(0), AtomId(1), AtomId(2), AtomId(3)], 8),
        Correspondence::from_images(&[BondId(0), BondId(1), BondId(2)], 6),
        Correspondence::from_images(&[DativeBondId(0)], 2),
        Correspondence::from_images(&[AromaticSystemId(0)], 2),
        Correspondence::from_images(&[MulticenterBondId(0)], 2),
        Correspondence::from_images(&[NoncovalentBondId(0)], 2),
        Correspondence::from_images(&[StereoAtomId(0)], 2),
        Correspondence::from_images(&[StereoBondId(0)], 2),
    );
    let expected_right = MoleculeCorrespondence::new(
        Correspondence::from_images(&[AtomId(4), AtomId(5), AtomId(6), AtomId(7)], 8),
        Correspondence::from_images(&[BondId(3), BondId(4), BondId(5)], 6),
        Correspondence::from_images(&[DativeBondId(1)], 2),
        Correspondence::from_images(&[AromaticSystemId(1)], 2),
        Correspondence::from_images(&[MulticenterBondId(1)], 2),
        Correspondence::from_images(&[NoncovalentBondId(1)], 2),
        Correspondence::from_images(&[StereoAtomId(1)], 2),
        Correspondence::from_images(&[StereoBondId(1)], 2),
    );

    assert_eq!(
        components,
        vec![
            (left, expected_left.reverse()),
            (right, expected_right.reverse())
        ]
    );
}
