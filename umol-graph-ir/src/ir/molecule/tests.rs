use std::collections::HashSet;
use std::ops::ControlFlow;
use std::sync::Arc;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, BipartiteMaximumMatchingAlgorithm,
    ConnectedComponentsAlgorithm, Correspondence, EdgeId, GeneralMaximumMatchingAlgorithm,
    MatchingEnumerationAlgorithm, MaximumIndependentSetAlgorithm, NodeId, NonBipartiteGraphError,
    RelevantCycleEnumerationAlgorithm, ShortestCycleAlgorithm, SimpleCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};
use umol_perm::Permutation;

use super::super::aromatic::AromaticSystemForm;
use super::super::atom::{AtomForm, ElementForm, IsotopeMassForm};
use super::super::bond::BondForm;
use super::super::boolean::BooleanForm;
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
use super::super::edit::{AtomFieldChange, AtomHandle, BondHandle, Edit, Edits};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
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
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
};
use super::{Molecule, MoleculeEntries, MoleculeIntegrityError, TransactionError};
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
            AtomId(0),
            AtomId(3),
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
            AtomId(0),
            AtomId(3),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
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
    |entries: &mut MoleculeEntries| entries.noncovalent[0].0 = AtomId(4),
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
#[case::dative_acceptor_is_donor(
    |entries: &mut MoleculeEntries| entries.dative[0].0[1] = AtomId(3),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::DativeBond(DativeBondId(0)),
        atom: AtomId(3),
    },
)]
#[case::dative_bonds_parallel(
    |entries: &mut MoleculeEntries| entries.dative.push((
        vec![AtomId(1)],
        AtomId(3),
        DativeBondForm::from_order(2),
    )),
    MoleculeIntegrityError::DativeBondsParallel {
        acceptor: AtomId(3),
        shared_donor: AtomId(1),
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
    |entries: &mut MoleculeEntries| entries.noncovalent[0].1 = AtomId(0),
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::NoncovalentBond(NoncovalentBondId(0)),
        atom: AtomId(0),
    },
)]
#[case::noncovalent_bonds_parallel_distinct_kinds(
    |entries: &mut MoleculeEntries| entries.noncovalent.push((
        AtomId(3),
        AtomId(0),
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
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::StereoAtom(StereoAtomId(0)),
        atom: AtomId(0),
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
    MoleculeIntegrityError::DuplicateParticipant {
        entity: Entity::StereoBond(StereoBondId(0)),
        atom: AtomId(0),
    },
)]
#[case::stereo_bond_sites_duplicate(
    |entries: &mut MoleculeEntries| {
        let stereo_bond = entries.stereo_bonds[0].clone();
        entries.stereo_bonds.push(stereo_bond);
    },
    MoleculeIntegrityError::StereoBondSitesDuplicate { bond: BondId(1) },
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
#[case::repeated_virtual_ligand_anchors(vec![
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
])]
fn test_molecule_try_from_entries(
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
#[case::axial_atom(Entity::StereoAtom(StereoAtomId(0)), StereoKind::Axial)]
#[case::cis_trans_atom(Entity::StereoAtom(StereoAtomId(0)), StereoKind::CisTrans)]
#[case::tetrahedral_bond(Entity::StereoBond(StereoBondId(0)), StereoKind::Tetrahedral)]
fn test_molecule_try_from_entries_stereo_kind(
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

    let molecule = Molecule::try_from_entries(entries).expect("entries satisfy molecule integrity");
    let actual = match entity {
        Entity::StereoAtom(id) => molecule.stereo_atom(id).kind(),
        Entity::StereoBond(id) => molecule.stereo_bond(id).kind(),
        _ => unreachable!("test cases contain only stereo entities"),
    };

    assert_eq!(actual, kind);
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
fn equiv_under_molecules(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) -> (Molecule, Molecule, MoleculeCorrespondence) {
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
        .map(|(first, second, attributes)| (map_atom(first), map_atom(second), attributes))
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
    let atom_correspondence = Correspondence::from_images(&atom_images, atom_images.len());
    let correspondence = MoleculeCorrespondence::induce(&left, &right, atom_correspondence)
        .expect("the atom correspondence describes the molecule pair");

    (left, right, correspondence)
}

#[rstest]
fn test_molecule_equiv_entity_data(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let base = Molecule::from_entries(entries.clone());

    let mut canonical_encoding = entries.clone();
    canonical_encoding.atoms[0].charge = NumForm::lit_set([1]);
    let canonical_encoding = Molecule::from_entries(canonical_encoding);
    assert_ne!(base, canonical_encoding);
    assert!(base.equiv(&canonical_encoding));

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
    noncovalent.noncovalent[0].2.kind = NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic);
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
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_equiv_relation_frames(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
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
    noncovalent.noncovalent[0].1 = AtomId(2);
    differences.push(Molecule::from_entries(noncovalent));

    let mut stereo_atom_site = entries.clone();
    stereo_atom_site.stereo_atoms[0].0 = AtomId(2);
    stereo_atom_site.stereo_atoms[0].1[1] = StereoLigand::new(AtomId(1), StereoLigandKind::Atom);
    differences.push(Molecule::from_entries(stereo_atom_site));

    let mut stereo_atom_ligand = entries.clone();
    stereo_atom_ligand.stereo_atoms[0].1[2] =
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair);
    differences.push(Molecule::from_entries(stereo_atom_ligand));

    let mut stereo_bond_site = entries.clone();
    stereo_bond_site.stereo_bonds[0].0 = BondId(2);
    differences.push(Molecule::from_entries(stereo_bond_site));

    let mut stereo_bond_ligand = entries;
    stereo_bond_ligand.stereo_bonds[0].1[1] = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);
    differences.push(Molecule::from_entries(stereo_bond_ligand));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 8],
    );
}

#[rstest]
fn test_molecule_equiv_structure_and_counts(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = Molecule::from_entries(entries.clone());
    let mut differences = Vec::new();

    let mut topology = entries.clone();
    topology.bonds[2].1 = AtomId(0);
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
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_equiv_under_non_identity(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;

    assert!(correspondence.is_total());
    assert!(!left.equiv(&right));
    assert!(left.equiv_under(&right, &correspondence));
    assert!(right.equiv_under(&left, &correspondence.reverse()));
}

#[rstest]
fn test_molecule_equiv_under_rejects_partial_correspondence(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;
    let partial = MoleculeCorrespondence::new(
        Correspondence::new(
            vec![
                (AtomId(0), AtomId(2)),
                (AtomId(1), AtomId(3)),
                (AtomId(2), AtomId(0)),
            ],
            4,
            4,
        )
        .expect("correspondence producer preserves partial-bijection invariants"),
        correspondence.bonds().clone(),
        correspondence.dative_bonds().clone(),
        correspondence.aromatic_systems().clone(),
        correspondence.multicenter_bonds().clone(),
        correspondence.noncovalent_bonds().clone(),
        correspondence.stereo_atoms().clone(),
        correspondence.stereo_bonds().clone(),
    );

    assert!(!left.equiv_under(&right, &partial));
}

#[rstest]
fn test_molecule_equiv_under_rejects_inconsistent_correspondence(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;
    let inconsistent = MoleculeCorrespondence::new(
        correspondence.atoms().clone(),
        Correspondence::from_images(&[BondId(1), BondId(0), BondId(2)], 3),
        correspondence.dative_bonds().clone(),
        correspondence.aromatic_systems().clone(),
        correspondence.multicenter_bonds().clone(),
        correspondence.noncovalent_bonds().clone(),
        correspondence.stereo_atoms().clone(),
        correspondence.stereo_bonds().clone(),
    );

    assert!(inconsistent.is_total());
    assert!(!left.equiv_under(&right, &inconsistent));
}

#[rstest]
fn test_molecule_remap(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;

    assert_eq!(left.remap(&correspondence), right);
}

#[rstest]
#[case::partial_correspondence(|molecule: Molecule, correspondence: &MoleculeCorrespondence| {
    let atoms = Correspondence::new(
        correspondence.atoms().matched_pairs()[..3].to_vec(),
        correspondence.atoms().left_count(),
        correspondence.atoms().right_count(),
    )
    .expect("the subset remains a partial bijection");
    (
        molecule,
        MoleculeCorrespondence::new(
            atoms,
            correspondence.bonds().clone(),
            correspondence.dative_bonds().clone(),
            correspondence.aromatic_systems().clone(),
            correspondence.multicenter_bonds().clone(),
            correspondence.noncovalent_bonds().clone(),
            correspondence.stereo_atoms().clone(),
            correspondence.stereo_bonds().clone(),
        ),
    )
})]
#[case::source_count(|molecule: Molecule, correspondence: &MoleculeCorrespondence| {
    (
        molecule,
        MoleculeCorrespondence::new(
            Correspondence::from_images(
                &[AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
                5,
            ),
            correspondence.bonds().clone(),
            correspondence.dative_bonds().clone(),
            correspondence.aromatic_systems().clone(),
            correspondence.multicenter_bonds().clone(),
            correspondence.noncovalent_bonds().clone(),
            correspondence.stereo_atoms().clone(),
            correspondence.stereo_bonds().clone(),
        ),
    )
})]
#[case::source_integrity(|mut molecule: Molecule, correspondence: &MoleculeCorrespondence| {
    molecule.constraints_mut().push(Constraint::Molecule(
        MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(4)]),
        },
    ));
    (molecule, correspondence.clone())
})]
fn test_molecule_try_remap_error(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
    #[case] prepare: fn(Molecule, &MoleculeCorrespondence) -> (Molecule, MoleculeCorrespondence),
) {
    let (left, _, correspondence) = case;
    let (left, correspondence) = prepare(left, &correspondence);

    assert_eq!(left.try_remap(&correspondence), None);
}

#[rstest]
#[should_panic(
    expected = "molecule remapping requires an integrity-valid source and a complete dense correspondence"
)]
fn test_molecule_remap_error(
    #[from(equiv_under_molecules)] case: (Molecule, Molecule, MoleculeCorrespondence),
) {
    let (left, _, correspondence) = case;
    let partial = MoleculeCorrespondence::new(
        Correspondence::new(
            correspondence.atoms().matched_pairs()[..3].to_vec(),
            correspondence.atoms().left_count(),
            correspondence.atoms().right_count(),
        )
        .expect("the subset remains a partial bijection"),
        correspondence.bonds().clone(),
        correspondence.dative_bonds().clone(),
        correspondence.aromatic_systems().clone(),
        correspondence.multicenter_bonds().clone(),
        correspondence.noncovalent_bonds().clone(),
        correspondence.stereo_atoms().clone(),
        correspondence.stereo_bonds().clone(),
    );

    left.remap(&partial);
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
#[case::stale_after_add(
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
    TransactionError::OldStateMismatch,
)]
fn test_molecule_apply_error(
    #[case] molecule: Molecule,
    #[case] edits: Edits,
    #[case] expected: TransactionError,
) {
    let original = molecule.clone();

    assert_eq!(molecule.apply(edits), Err(expected));
    assert_eq!(molecule, original);
}

#[rstest]
fn test_molecule_extract(#[from(rich_molecule)] molecule: Molecule) {
    let sub = molecule.induced_subgraph(&[AtomId(0), AtomId(1)]);
    let extracted = molecule.extract(&sub);
    assert_eq!(extracted.atoms().count(), 2);
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
    let _compaction = b.remove(&[AtomId(3)], &[BondId(2)]);
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
#[case::empty(Vec::new(), Molecule::new(), Vec::new())]
#[case::singleton(
    vec![Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    })],
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    }),
    vec![vec![(AtomId(0), AtomId(0))]],
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
    vec![
        vec![(AtomId(0), AtomId(0))],
        vec![],
        vec![(AtomId(0), AtomId(1)), (AtomId(1), AtomId(2))],
    ],
)]
fn test_molecule_combine_all(
    #[case] molecules: Vec<Molecule>,
    #[case] expected: Molecule,
    #[case] expected_atom_matched_pairs: Vec<Vec<(AtomId, AtomId)>>,
) {
    let (combined, correspondences) = Molecule::combine_all(&molecules);

    assert_eq!(combined, expected);
    assert_eq!(
        correspondences
            .iter()
            .map(|correspondence| correspondence.atoms().matched_pairs().to_vec())
            .collect::<Vec<_>>(),
        expected_atom_matched_pairs,
    );
    for (molecule, correspondence) in molecules.iter().zip(&correspondences) {
        assert_eq!(combined.extract(correspondence), *molecule);
        assert_eq!(
            correspondence.atoms().right_count(),
            combined.atoms().count()
        );
        assert_eq!(
            correspondence.bonds().right_count(),
            combined.bonds().count()
        );
        assert_eq!(
            correspondence.dative_bonds().right_count(),
            combined.dative_bonds().count()
        );
        assert_eq!(
            correspondence.aromatic_systems().right_count(),
            combined.aromatic_systems().count()
        );
        assert_eq!(
            correspondence.multicenter_bonds().right_count(),
            combined.multicenter_bonds().count()
        );
        assert_eq!(
            correspondence.noncovalent_bonds().right_count(),
            combined.noncovalent_bonds().count()
        );
        assert_eq!(
            correspondence.stereo_atoms().right_count(),
            combined.stereo_atoms().count()
        );
        assert_eq!(
            correspondence.stereo_bonds().right_count(),
            combined.stereo_bonds().count()
        );
    }
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
    let (union, correspondence) = left.combine(&right);

    assert_eq!(union.atoms().count(), 4);
    assert_eq!(union.bonds().count(), 2);
    assert_eq!(union.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    assert_eq!(union.bond(BondId(1)).atom_ids(), [AtomId(2), AtomId(3)]);
    assert_eq!(union.bond(BondId(1)).attributes, &BondForm::from_order(2));
    // right's ids map to their offset union ids; left's are the prefix (unchanged)
    assert_eq!(correspondence.atoms().right_of(AtomId(0)), Some(AtomId(2)));
    assert_eq!(correspondence.atoms().right_of(AtomId(1)), Some(AtomId(3)));
    assert_eq!(correspondence.bonds().right_of(BondId(0)), Some(BondId(1)));
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
    let correspondence = left.combine_from(&right);

    assert_eq!(left.atoms().count(), 3);
    assert_eq!(left.bond(BondId(0)).atom_ids(), [AtomId(1), AtomId(2)]);
    assert_eq!(correspondence.atoms().right_of(AtomId(0)), Some(AtomId(1)));
    assert_eq!(correspondence.atoms().right_of(AtomId(1)), Some(AtomId(2)));
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
    let (union, correspondence) = left.combine(&right);

    assert_eq!(union.aromatic_systems().count(), 1);
    // right's overlay over its atoms [0, 1] shifts by left's one atom
    assert_eq!(
        union
            .aromatic_system(AromaticSystemId(0))
            .atom_ids()
            .collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2)]
    );
    assert_eq!(
        correspondence
            .aromatic_systems()
            .right_of(AromaticSystemId(0)),
        Some(AromaticSystemId(0))
    );
}

#[rstest]
fn test_molecule_combine_stereo() {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
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
        constraints: Constraints::new(),
        ..Default::default()
    });
    let (union, _) = left.combine(&right);

    assert_eq!(union.stereo_atoms().count(), 1);
    let stereo = union.stereo_atoms().iter().next().unwrap();
    // right's site (atom 0) and ligands (atoms 1..=4) shift by left's one atom
    assert_eq!(stereo.site_id(), AtomId(1));
    assert_eq!(
        stereo.ligands().map(|l| l.atom_id()).collect::<Vec<_>>(),
        vec![AtomId(2), AtomId(3), AtomId(4), AtomId(5)]
    );
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
    let (union, _) = left.combine(&right);

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
    let components = mol.split();

    assert_eq!(components.len(), 2);
    let (first, first_corr) = &components[0];
    assert_eq!(first.atoms().count(), 2);
    assert_eq!(first.bond(BondId(0)).attributes, &BondForm::from_order(1));
    assert_eq!(first_corr.atoms().right_of(AtomId(0)), Some(AtomId(0)));
    assert_eq!(first_corr.atoms().right_of(AtomId(1)), Some(AtomId(1)));
    let (second, second_corr) = &components[1];
    assert_eq!(second.bond(BondId(0)).attributes, &BondForm::from_order(2));
    assert_eq!(second_corr.atoms().right_of(AtomId(0)), Some(AtomId(2)));
    assert_eq!(second_corr.atoms().right_of(AtomId(1)), Some(AtomId(3)));
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
    let components = mol.split();

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
    let (union, _) = left.combine(&right);
    let components = union.split();

    assert_eq!(components.len(), 2);
    assert_eq!(components[0].0, left);
    assert_eq!(components[1].0, right);
}

#[rstest]
fn test_molecule_split_stereo() {
    // a stereo atom binds its site + ligands into one component, separate from a lone bond
    let mol = Molecule::from_entries(MoleculeEntries {
        atoms: (0..7).map(|_| AtomForm::from_element(Element::C)).collect(),
        bonds: vec![(AtomId(5), AtomId(6), BondForm::from_order(1))],
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
        constraints: Constraints::new(),
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 2);
    let (bound, _) = &components[0];
    assert_eq!(bound.atoms().count(), 5);
    assert_eq!(bound.stereo_atoms().count(), 1);
    let stereo = bound.stereo_atoms().iter().next().unwrap();
    assert_eq!(stereo.site_id(), AtomId(0));
    assert_eq!(
        stereo.ligands().map(|l| l.atom_id()).collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)]
    );
    let (lone, _) = &components[1];
    assert_eq!(lone.atoms().count(), 2);
    assert_eq!(lone.stereo_atoms().count(), 0);
    assert_eq!(lone.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
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
    let components = mol.split();

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
    let components = mol.split();

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
