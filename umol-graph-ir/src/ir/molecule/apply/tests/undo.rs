use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{Compaction, GraphCompaction};

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::{AtomForm, ElementForm, IsotopeMassForm};
use crate::ir::bond::BondForm;
use crate::ir::compact::MoleculeCompaction;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AtomConstraintForm,
    AtomConstraintKey, BondConstraintForm, BondConstraintKey, Constraint, DativeBondConstraintForm,
    DativeBondConstraintKey, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, StereoAtomConstraintKey,
    StereoBondConstraintKey,
};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AddAromaticSystem, AddBond, AddMulticenterBond, AddStereoAtom, AddStereoBond,
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemFieldChange,
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    CascadedConstraints, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit, Edits,
    EntityHandle, ModifiedConstraint, MulticenterBondFieldChange, MulticenterBondHandle,
    NoncovalentBondFieldChange, RemovedAromaticSystem, RemovedAromaticSystems, RemovedAtom,
    RemovedConstraint, RemovedOverlays, RemovedTopology, StereoAtomFieldChange,
    StereoBondFieldChange, Undo,
};
use crate::ir::electrons::ElectronCountsForm;
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
};
use crate::ir::{Molecule, RingScope};

#[fixture]
fn empty() -> Molecule {
    Molecule::new()
}

#[fixture]
fn one_atom() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    })
}

#[fixture]
fn batched_overlays() -> Molecule {
    let mut editor = Molecule::default().edit();
    for _ in 0..6 {
        editor.add_atom(AtomForm::from_element(Element::C));
    }
    for index in 0..3_u32 {
        let first = AtomId(index * 2);
        let second = AtomId(index * 2 + 1);
        let bond = editor.add_bond(first, second, BondForm::from_order(1));
        editor.add_dative_bond(&[first], second, DativeBondForm::from_order(1));
        editor.add_aromatic_system(&[first, second], AromaticSystemForm::default());
        editor.add_multicenter_bond(&[first, second], MulticenterBondForm::default());
        editor.add_noncovalent_bond(
            [first, second],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        editor.add_stereo_atom(
            first,
            &[
                StereoLigand::new(second, StereoLigandKind::Atom),
                StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(first, StereoLigandKind::LonePair),
            ],
            StereoAtomForm::default(),
        );
        editor.add_stereo_bond(
            bond,
            &[
                StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(first, StereoLigandKind::LonePair),
                StereoLigand::new(second, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(second, StereoLigandKind::LonePair),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
    }
    editor.finish().unwrap()
}

#[rstest]
fn test_molecule_apply_undo_additions(removal_entries: MoleculeEntries) {
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let mut journal = Vec::new();
    for edit in [
        Edit::AddAtoms {
            atoms: vec![AtomForm::from_element(Element::F)],
        },
        Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::Id(AtomId(3)), AtomHandle::New(0)],
                attributes: BondForm::from_order(1),
            }],
        },
        Edit::AddDativeBond {
            donors: vec![AtomHandle::New(0)],
            acceptor: AtomHandle::Id(AtomId(3)),
            attributes: DativeBondForm::from_order(1),
        },
        Edit::AddAromaticSystem(Box::new(AddAromaticSystem {
            atoms: vec![AtomHandle::Id(AtomId(3)), AtomHandle::New(0)],
            attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
        })),
        Edit::AddMulticenterBond(Box::new(AddMulticenterBond {
            atoms: vec![AtomHandle::Id(AtomId(3)), AtomHandle::New(0)],
            attributes: MulticenterBondForm::from_electrons(vec![1, 2]),
        })),
        Edit::AddNoncovalentBond {
            atoms: [AtomHandle::Id(AtomId(2)), AtomHandle::New(0)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        },
        Edit::AddStereoAtom(Box::new(AddStereoAtom {
            site: AtomHandle::New(0),
            ligands: vec![(AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom)],
            attributes: StereoAtomForm::default(),
        })),
        Edit::AddStereoBond(Box::new(AddStereoBond {
            site: BondHandle::New(0),
            ligands: vec![
                (
                    AtomHandle::Id(AtomId(3)),
                    StereoLigandKind::ImplicitHydrogen,
                ),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::LonePair),
                (AtomHandle::New(0), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::New(0), StereoLigandKind::LonePair),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        })),
        Edit::AddMoleculeConstraint {
            constraint: Box::new(
                ConstraintEdit::new(
                    Constraint::Atom(AtomId(4), AtomConstraintForm::degree(1)),
                    |_| Some(EntityHandle::Atom(AtomHandle::New(0))),
                )
                .unwrap(),
            ),
        },
    ] {
        journal.push(
            molecule
                .apply_edit_with_undo(edit, &mut state)
                .unwrap()
                .unwrap(),
        );
    }
    let before_removal = molecule.clone();
    let removed = molecule
        .apply_edit_with_undo(
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            &mut state,
        )
        .unwrap()
        .unwrap();

    molecule.apply_undo(removed);
    assert_eq!(molecule, before_removal);
    for undo in journal.into_iter().rev() {
        molecule.apply_undo(undo);
    }
    assert_eq!(molecule, initial);
}

#[rstest]
fn test_molecule_apply_undo_added_topology() {
    let initial = MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
        ],
        ..Default::default()
    };
    let mut extended = initial.clone();
    extended.atoms.push(AtomForm::from_element(Element::F));
    extended
        .bonds
        .push((AtomId(2), AtomId(3), BondForm::from_order(1)));
    let mut molecule = Molecule::from_entries(extended);

    molecule.apply_undo(Undo::RemoveAddedTopology {
        atoms: vec![AddedAtom {
            id: AtomId(3),
            attributes: AtomForm::from_element(Element::F),
        }],
        bonds: vec![AddedBond {
            id: BondId(2),
            endpoints: [AtomId(2), AtomId(3)],
            attributes: BondForm::from_order(1),
        }],
    });

    assert_eq!(molecule, Molecule::from_entries(initial));
}

#[rstest]
fn test_molecule_apply_undo_added_topology_duplicate(mut one_atom: Molecule) {
    let added = AddedAtom {
        id: AtomId(0),
        attributes: AtomForm::from_element(Element::C),
    };
    one_atom.apply_undo(Undo::RemoveAddedTopology {
        atoms: vec![added.clone(), added],
        bonds: Vec::new(),
    });
}

#[rstest]
#[case::topology(Undo::RemoveAddedTopology {
    atoms: vec![AddedAtom { id: AtomId(0), attributes: AtomForm::default() }],
    bonds: vec![AddedBond { id: BondId(0), endpoints: [AtomId(0), AtomId(1)], attributes: BondForm::default() }],
})]
#[case::dative_bond(Undo::RemoveAddedDativeBond(AddedDativeBond { id: DativeBondId(0), donors: vec![], acceptor: AtomId(0), attributes: DativeBondForm::default() }))]
#[case::aromatic_system(Undo::RemoveAddedAromaticSystem(Box::new(AddedAromaticSystem { id: AromaticSystemId(0), atoms: vec![], attributes: AromaticSystemForm::default() })))]
#[case::multicenter_bond(Undo::RemoveAddedMulticenterBond(Box::new(AddedMulticenterBond { id: MulticenterBondId(0), atoms: vec![], attributes: MulticenterBondForm::default() })))]
#[case::noncovalent_bond(Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond { id: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)], attributes: NoncovalentBondForm::default() }))]
#[case::stereo_atom(Undo::RemoveAddedStereoAtom(Box::new(AddedStereoAtom { id: StereoAtomId(0), site: AtomId(0), ligands: vec![], attributes: StereoAtomForm::default() })))]
#[case::stereo_bond(Undo::RemoveAddedStereoBond(Box::new(AddedStereoBond { id: StereoBondId(0), site: BondId(0), ligands: vec![], attributes: StereoBondForm::default() })))]
fn test_molecule_apply_undo_added_entry_manipulated(#[case] undo: Undo) {
    Molecule::default().apply_undo(undo);
}

#[rstest]
fn test_molecule_apply_undo_removed_addition() {
    let mut molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C)],
        ..Default::default()
    });
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(
            Edit::AddAtoms {
                atoms: vec![AtomForm::from_element(Element::N)],
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    molecule.remove_topology(&[AtomId(1)], &[]);
    molecule.apply_undo(undo);
}

#[rstest]
#[case::atoms_short([0, 0, 0, 0, 0, 0, 0, 0])]
#[case::atoms_long([2, 0, 0, 0, 0, 0, 0, 0])]
#[case::bonds([1, 1, 0, 0, 0, 0, 0, 0])]
#[case::dative([1, 0, 1, 0, 0, 0, 0, 0])]
#[case::aromatic([1, 0, 0, 1, 0, 0, 0, 0])]
#[case::multicenter([1, 0, 0, 0, 1, 0, 0, 0])]
#[case::noncovalent([1, 0, 0, 0, 0, 1, 0, 0])]
#[case::stereo_atoms([1, 0, 0, 0, 0, 0, 1, 0])]
#[case::stereo_bonds([1, 0, 0, 0, 0, 0, 0, 1])]
fn test_molecule_apply_undo_compaction_counts(mut one_atom: Molecule, #[case] counts: [usize; 8]) {
    let compaction = MoleculeCompaction::new(
        GraphCompaction::new(
            Compaction::identity(counts[0]),
            Compaction::identity(counts[1]),
        ),
        Compaction::identity(counts[2]),
        Compaction::identity(counts[3]),
        Compaction::identity(counts[4]),
        Compaction::identity(counts[5]),
        Compaction::identity(counts[6]),
        Compaction::identity(counts[7]),
    );
    one_atom.apply_undo(Undo::RestoreRemovedTopology(Box::new(RemovedTopology {
        atoms: Vec::new(),
        bonds: Vec::new(),
        overlays: RemovedOverlays::default(),
        undo_compaction: compaction.undo_compaction(),
        compaction,
        cascade: CascadedConstraints::default(),
    })));
}

#[rstest]
fn test_molecule_apply_undo_compaction_dimension(mut one_atom: Molecule) {
    let compaction = MoleculeCompaction::empty();
    one_atom.apply_undo(Undo::RestoreRemovedTopology(Box::new(RemovedTopology {
        atoms: vec![RemovedAtom {
            id: AtomId(0),
            attributes: AtomForm::from_element(Element::N),
        }],
        bonds: Vec::new(),
        overlays: RemovedOverlays::default(),
        undo_compaction: compaction.undo_compaction(),
        compaction,
        cascade: CascadedConstraints::default(),
    })));
}

#[rstest]
fn test_molecule_apply_undo_reconstruction_entry(mut empty: Molecule) {
    empty.apply_undo(Undo::RestoreRemovedAromaticSystems(Box::new(
        RemovedAromaticSystems {
            removed: vec![RemovedAromaticSystem {
                id: AromaticSystemId(1),
                atoms: Vec::new(),
                attributes: AromaticSystemForm::default(),
            }],
            undo_compaction: MoleculeCompaction::empty().undo_compaction(),
            cascade: CascadedConstraints::default(),
        },
    )));
}

#[rstest]
#[case::dative_donors(Undo::RestoreDativeBondDonors {
        id: DativeBondId(9), donors: vec![AtomId(0)],
    })]
#[case::dative_acceptor(Undo::RestoreDativeBondAcceptor {
        id: DativeBondId(9), acceptor: AtomId(0),
    })]
#[case::aromatic_atoms(Undo::RestoreAromaticSystemAtoms {
        id: AromaticSystemId(9), atoms: vec![AtomId(0)],
    })]
#[case::multicenter_atoms(Undo::RestoreMulticenterBondAtoms {
        id: MulticenterBondId(9), atoms: vec![AtomId(0)],
    })]
#[case::noncovalent_atoms(Undo::RestoreNoncovalentBondAtoms {
        id: NoncovalentBondId(9), atoms: [AtomId(0), AtomId(1)],
    })]
#[case::stereo_atom_site(Undo::RestoreStereoAtomSite {
        id: StereoAtomId(9), site: AtomId(0),
    })]
#[case::stereo_atom_ligands(Undo::RestoreStereoAtomLigands {
        id: StereoAtomId(9), ligands: vec![],
    })]
#[case::stereo_bond_site(Undo::RestoreStereoBondSite {
        id: StereoBondId(9), site: BondId(0),
    })]
#[case::stereo_bond_ligands(Undo::RestoreStereoBondLigands {
        id: StereoBondId(9), ligands: vec![],
    })]
fn test_molecule_apply_undo_restore_target_manipulated(
    mut batched_overlays: Molecule,
    #[case] undo: Undo,
) {
    batched_overlays.apply_undo(undo);
}

#[rstest]
#[case::atom(Undo::ModifyAtomField {
    id: AtomId(0),
    change: AtomFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
})]
#[case::bond(Undo::ModifyBondField {
    id: BondId(0),
    change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(1) },
})]
#[case::dative_bond(Undo::ModifyDativeBondField {
    id: DativeBondId(0),
    change: DativeBondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(1) },
})]
#[case::aromatic_system(Undo::ModifyAromaticSystemField {
    id: AromaticSystemId(0),
    change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
})]
#[case::multicenter_bond(Undo::ModifyMulticenterBondField {
    id: MulticenterBondId(0),
    change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
})]
#[case::noncovalent_bond(Undo::ModifyNoncovalentBondField {
    id: NoncovalentBondId(0),
    change: NoncovalentBondFieldChange::Kind {
        old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
        new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
    },
})]
#[case::stereo_atom(Undo::ModifyStereoAtomField {
    id: StereoAtomId(0),
    change: StereoAtomFieldChange::Configuration {
        old: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        new: StereoConfigurationForm::Undetermined,
    },
})]
#[case::stereo_bond(Undo::ModifyStereoBondField {
    id: StereoBondId(0),
    change: StereoBondFieldChange::Configuration {
        old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
        new: StereoConfigurationForm::Undetermined,
    },
})]
fn test_molecule_apply_undo_field_manipulated(#[case] undo: Undo) {
    Molecule::default().apply_undo(undo);
}

#[rstest]
#[case(Undo::ModifyAtomElements { changes: vec![(AtomId(7), ElementForm::Lit(Element::C), ElementForm::Undetermined)] })]
#[case(Undo::ModifyAtomIsotopeMasses { changes: vec![(AtomId(7), IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined)] })]
#[case(Undo::ModifyAtomCharges { changes: vec![(AtomId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyAtomImplicitHydrogens { changes: vec![(AtomId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyAtomLonePairs { changes: vec![(AtomId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyAtomUnpairedElectrons { changes: vec![(AtomId(7), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::default())] })]
#[case(Undo::RestoreAtomConstraints { key: AtomConstraintKey::Valence, constraints: vec![(AtomId(7), Some(AtomConstraintForm::degree(2)))] })]
#[case(Undo::ModifyBondOrders { changes: vec![(BondId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyBondCharges { changes: vec![(BondId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyBondUnpairedElectrons { changes: vec![(BondId(7), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::default())] })]
#[case(Undo::RestoreBondConstraints { key: BondConstraintKey::Aromatic, constraints: vec![(BondId(7), Some(BondConstraintForm::aromatic(true)))] })]
#[case(Undo::ModifyDativeBondOrders { changes: vec![(DativeBondId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::RestoreDativeBondConstraints { key: DativeBondConstraintKey::Aromatic, constraints: vec![(DativeBondId(7), Some(DativeBondConstraintForm::undetermined(DativeBondConstraintKey::Aromatic)))] })]
#[case(Undo::ModifyAromaticSystemElectrons { changes: vec![(AromaticSystemId(7), ElectronCountsForm::Lit(vec![1, 2]), ElectronCountsForm::Undetermined)] })]
#[case(Undo::ModifyAromaticSystemCharges { changes: vec![(AromaticSystemId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyAromaticSystemUnpairedElectrons { changes: vec![(AromaticSystemId(7), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::default())] })]
#[case(Undo::RestoreAromaticSystemConstraints { key: AromaticSystemConstraintKey::ElectronCount, constraints: vec![(AromaticSystemId(7), Some(AromaticSystemConstraintForm::undetermined(AromaticSystemConstraintKey::ElectronCount)))] })]
#[case(Undo::ModifyMulticenterBondElectrons { changes: vec![(MulticenterBondId(7), ElectronCountsForm::Lit(vec![1, 2]), ElectronCountsForm::Undetermined)] })]
#[case(Undo::ModifyMulticenterBondCharges { changes: vec![(MulticenterBondId(7), NumForm::Lit(1), NumForm::Undetermined)] })]
#[case(Undo::ModifyMulticenterBondUnpairedElectrons { changes: vec![(MulticenterBondId(7), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm::default())] })]
#[case(Undo::RestoreMulticenterBondConstraints { key: MulticenterBondConstraintKey::ElectronCount, constraints: vec![(MulticenterBondId(7), Some(MulticenterBondConstraintForm::undetermined(MulticenterBondConstraintKey::ElectronCount)))] })]
#[case(Undo::ModifyNoncovalentBondKinds { changes: vec![(NoncovalentBondId(7), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Undetermined)] })]
#[case(Undo::RestoreNoncovalentBondConstraints { key: NoncovalentBondConstraintKey::Intramolecular, constraints: vec![(NoncovalentBondId(7), Some(NoncovalentBondConstraintForm::undetermined(NoncovalentBondConstraintKey::Intramolecular)))] })]
#[case(Undo::ModifyStereoAtomConfigurations { changes: vec![(StereoAtomId(7), StereoConfigurationForm::Undetermined, StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))] })]
#[case(Undo::RestoreStereoAtomConstraints { key: StereoAtomConstraintKey::Stereogenicity, constraints: vec![(StereoAtomId(7), None)] })]
#[case(Undo::ModifyStereoBondConfigurations { changes: vec![(StereoBondId(7), StereoConfigurationForm::Undetermined, StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))] })]
#[case(Undo::RestoreStereoBondConstraints { key: StereoBondConstraintKey::Stereogenicity, constraints: vec![(StereoBondId(7), None)] })]
fn test_molecule_apply_undo_batches(#[case] undo: Undo) {
    let mut molecule = Molecule::new();
    molecule.apply_undo(undo);
    assert_eq!(molecule, Molecule::new());
}

#[rstest]
#[case::order(
    BondForm::from_order(1),
    Edit::ModifyBondOrders { changes: vec![(BondHandle::Id(BondId(0)), NumForm::lit_set([1]), NumForm::Lit(2)), (BondHandle::Id(BondId(0)), NumForm::Lit(2), NumForm::Lit(3))] },
    Undo::ModifyBondOrders { changes: vec![(BondId(0), NumForm::Lit(2), NumForm::lit_set([1])), (BondId(0), NumForm::Lit(3), NumForm::Lit(2))] },
    BondForm { order: NumForm::lit_set([1]), ..Default::default() },
)]
#[case::constraint(
    BondForm { constraints: vec![BondConstraintForm::ring_membership(RingScope::All, 1)].into(), ..Default::default() },
    Edit::ModifyBondConstraints { key: BondConstraintKey::RingMembership(RingScope::All), changes: vec![(BondHandle::Id(BondId(0)), Some(BondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))), None)] },
    Undo::RestoreBondConstraints { key: BondConstraintKey::RingMembership(RingScope::All), constraints: vec![(BondId(0), Some(BondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))))] },
    BondForm { constraints: vec![BondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))].into(), ..Default::default() },
)]
fn test_molecule_apply_undo_bond_batches_normalized(
    #[case] initial: BondForm,
    #[case] edit: Edit,
    #[case] expected_undo: Undo,
    #[case] restored: BondForm,
) {
    let mut molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        bonds: vec![(AtomId(0), AtomId(1), initial)],
        ..Default::default()
    });
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    assert_eq!(undo, expected_undo);
    molecule.apply_undo(undo);
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 2],
            bonds: vec![(AtomId(0), AtomId(1), restored)],
            ..Default::default()
        })
    );
}

#[rstest]
#[case::dative_bond_order(
    Edit::ModifyDativeBondOrders { changes: vec![(DativeBondHandle::Id(DativeBondId(0)), NumForm::lit_set([1]), NumForm::Lit(2)), (DativeBondHandle::Id(DativeBondId(0)), NumForm::Lit(2), NumForm::Lit(3))] },
    Undo::ModifyDativeBondOrders { changes: vec![(DativeBondId(0), NumForm::Lit(2), NumForm::lit_set([1])), (DativeBondId(0), NumForm::Lit(3), NumForm::Lit(2))] },
    |entries: &mut MoleculeEntries| entries.dative[0].2.order = NumForm::lit_set([1]),
)]
#[case::aromatic_system_charge(
    Edit::ModifyAromaticSystemCharges { changes: vec![(AromaticSystemHandle::Id(AromaticSystemId(0)), NumForm::lit_set([1]), NumForm::Lit(2)), (AromaticSystemHandle::Id(AromaticSystemId(0)), NumForm::Lit(2), NumForm::Lit(3))] },
    Undo::ModifyAromaticSystemCharges { changes: vec![(AromaticSystemId(0), NumForm::Lit(2), NumForm::lit_set([1])), (AromaticSystemId(0), NumForm::Lit(3), NumForm::Lit(2))] },
    |entries: &mut MoleculeEntries| entries.aromatic[0].1.charge = NumForm::lit_set([1]),
)]
#[case::multicenter_bond_charge(
    Edit::ModifyMulticenterBondCharges { changes: vec![(MulticenterBondHandle::Id(MulticenterBondId(0)), NumForm::lit_set([1]), NumForm::Lit(2)), (MulticenterBondHandle::Id(MulticenterBondId(0)), NumForm::Lit(2), NumForm::Lit(3))] },
    Undo::ModifyMulticenterBondCharges { changes: vec![(MulticenterBondId(0), NumForm::Lit(2), NumForm::lit_set([1])), (MulticenterBondId(0), NumForm::Lit(3), NumForm::Lit(2))] },
    |entries: &mut MoleculeEntries| entries.multicenter[0].1.charge = NumForm::lit_set([1]),
)]
#[case::dative_bond_constraint(
    Edit::ModifyDativeBondConstraints { key: DativeBondConstraintKey::RingMembership(RingScope::All), changes: vec![(DativeBondHandle::Id(DativeBondId(0)), Some(DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))), None)] },
    Undo::RestoreDativeBondConstraints { key: DativeBondConstraintKey::RingMembership(RingScope::All), constraints: vec![(DativeBondId(0), Some(DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))))] },
    |entries: &mut MoleculeEntries| entries.dative[0].2.constraints = vec![DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::lit_set([1]))].into(),
)]
#[case::aromatic_system_constraint(
    Edit::ModifyAromaticSystemConstraints { key: AromaticSystemConstraintKey::ElectronCount, changes: vec![(AromaticSystemHandle::Id(AromaticSystemId(0)), Some(AromaticSystemConstraintForm::electron_count(NumForm::lit_set([1]))), None)] },
    Undo::RestoreAromaticSystemConstraints { key: AromaticSystemConstraintKey::ElectronCount, constraints: vec![(AromaticSystemId(0), Some(AromaticSystemConstraintForm::electron_count(NumForm::lit_set([1]))))] },
    |entries: &mut MoleculeEntries| entries.aromatic[0].1.constraints = vec![AromaticSystemConstraintForm::electron_count(NumForm::lit_set([1]))].into(),
)]
#[case::multicenter_bond_constraint(
    Edit::ModifyMulticenterBondConstraints { key: MulticenterBondConstraintKey::ElectronCount, changes: vec![(MulticenterBondHandle::Id(MulticenterBondId(0)), Some(MulticenterBondConstraintForm::electron_count(NumForm::lit_set([1]))), None)] },
    Undo::RestoreMulticenterBondConstraints { key: MulticenterBondConstraintKey::ElectronCount, constraints: vec![(MulticenterBondId(0), Some(MulticenterBondConstraintForm::electron_count(NumForm::lit_set([1]))))] },
    |entries: &mut MoleculeEntries| entries.multicenter[0].1.constraints = vec![MulticenterBondConstraintForm::electron_count(NumForm::lit_set([1]))].into(),
)]
fn test_molecule_apply_undo_relation_batches_normalized(
    mut removal_entries: MoleculeEntries,
    #[case] edit: Edit,
    #[case] expected_undo: Undo,
    #[case] restore: fn(&mut MoleculeEntries),
) {
    removal_entries.aromatic[0].1.charge = NumForm::Lit(1);
    removal_entries.multicenter[0].1.charge = NumForm::Lit(1);
    removal_entries.dative[0].2.constraints =
        vec![DativeBondConstraintForm::ring_membership(RingScope::All, 1)].into();
    removal_entries.aromatic[0].1.constraints =
        vec![AromaticSystemConstraintForm::electron_count(1)].into();
    removal_entries.multicenter[0].1.constraints =
        vec![MulticenterBondConstraintForm::electron_count(1)].into();
    let mut molecule = Molecule::from_entries(removal_entries.clone());
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    assert_eq!(undo, expected_undo);
    molecule.apply_undo(undo);
    restore(&mut removal_entries);
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
}

#[rstest]
#[case::atom(EntityKind::Atom)]
#[case::bond(EntityKind::Bond)]
#[case::dative_bond(EntityKind::DativeBond)]
#[case::aromatic_system(EntityKind::AromaticSystem)]
#[case::multicenter_bond(EntityKind::MulticenterBond)]
#[case::noncovalent_bond(EntityKind::NoncovalentBond)]
#[case::stereo_atom(EntityKind::StereoAtom)]
#[case::stereo_bond(EntityKind::StereoBond)]
fn test_molecule_apply_undo_field_receiver(mut empty: Molecule, #[case] kind: EntityKind) {
    let undo = match kind {
        EntityKind::Atom => Undo::ModifyAtomField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::default(),
            },
        },
        EntityKind::Bond => Undo::ModifyBondField {
            id: BondId(0),
            change: BondFieldChange::Order {
                old: NumForm::Lit(2),
                new: NumForm::Lit(1),
            },
        },
        EntityKind::DativeBond => Undo::ModifyDativeBondField {
            id: DativeBondId(0),
            change: DativeBondFieldChange::Order {
                old: NumForm::Lit(2),
                new: NumForm::Lit(1),
            },
        },
        EntityKind::AromaticSystem => Undo::ModifyAromaticSystemField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::default(),
            },
        },
        EntityKind::MulticenterBond => Undo::ModifyMulticenterBondField {
            id: MulticenterBondId(0),
            change: MulticenterBondFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::default(),
            },
        },
        EntityKind::NoncovalentBond => Undo::ModifyNoncovalentBondField {
            id: NoncovalentBondId(0),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            },
        },
        EntityKind::StereoAtom => Undo::ModifyStereoAtomField {
            id: StereoAtomId(0),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        },
        EntityKind::StereoBond => Undo::ModifyStereoBondField {
            id: StereoBondId(0),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
            },
        },
    };

    empty.apply_undo(undo);
}

#[rstest]
#[case::atom(Undo::RestoreAtomConstraint {
    id: AtomId(0), key: AtomConstraintKey::Degree, constraint: None,
})]
#[case::bond(Undo::RestoreBondConstraint {
    id: BondId(0), key: BondConstraintKey::Aromatic, constraint: None,
})]
#[case::dative_bond(Undo::RestoreDativeBondConstraint {
    id: DativeBondId(0), key: DativeBondConstraintKey::Aromatic, constraint: None,
})]
#[case::aromatic_system(Undo::RestoreAromaticSystemConstraint {
    id: AromaticSystemId(0), key: AromaticSystemConstraintKey::ElectronCount, constraint: None,
})]
#[case::multicenter_bond(Undo::RestoreMulticenterBondConstraint {
    id: MulticenterBondId(0), key: MulticenterBondConstraintKey::ElectronCount, constraint: None,
})]
#[case::noncovalent_bond(Undo::RestoreNoncovalentBondConstraint {
    id: NoncovalentBondId(0), key: NoncovalentBondConstraintKey::Intramolecular, constraint: None,
})]
#[case::stereo_atom(Undo::RestoreStereoAtomConstraint {
    id: StereoAtomId(0), key: StereoAtomConstraintKey::Stereogenicity, constraint: None,
})]
#[case::stereo_bond(Undo::RestoreStereoBondConstraint {
    id: StereoBondId(0), key: StereoBondConstraintKey::Stereogenicity, constraint: None,
})]
#[case::molecule(Undo::RemoveAddedMoleculeConstraint { position: 0 })]
fn test_molecule_apply_undo_constraint_manipulated(#[case] undo: Undo) {
    Molecule::default().apply_undo(undo);
}

#[rstest]
fn test_molecule_apply_undo_constraint_position(mut removal_entries: MoleculeEntries) {
    let removed = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    let repeated = Constraint::Atom(AtomId(3), AtomConstraintForm::degree(2));
    removal_entries.constraints = vec![removed, repeated.clone()].into();
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let added = molecule
        .apply_edit_with_undo(
            Edit::AddMoleculeConstraint {
                constraint: Box::new(repeated.clone().into()),
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(added, Undo::RemoveAddedMoleculeConstraint { position: 2 });
    let before_removal = molecule.clone();
    let topology = molecule
        .apply_edit_with_undo(
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    let after_removal = molecule.clone();
    let constraint = molecule
        .apply_edit_with_undo(
            Edit::RemoveMoleculeConstraint {
                constraint: Box::new(repeated.into()),
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        constraint,
        Undo::RestoreMoleculeConstraints(CascadedConstraints {
            removed: vec![
                RemovedConstraint {
                    position: 0,
                    constraint: Constraint::Atom(AtomId(2), AtomConstraintForm::degree(2)),
                },
                RemovedConstraint {
                    position: 1,
                    constraint: Constraint::Atom(AtomId(2), AtomConstraintForm::degree(2)),
                },
            ],
            modified: Vec::new(),
        })
    );

    molecule.apply_undo(constraint);
    assert_eq!(molecule, after_removal);
    molecule.apply_undo(topology);
    assert_eq!(molecule, before_removal);
    molecule.apply_undo(added);
    assert_eq!(molecule, initial);
}

#[rstest]
#[case::missing(7)]
#[case::overlapping(0)]
fn test_molecule_apply_undo_constraint_history(mut one_atom: Molecule, #[case] position: usize) {
    let constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    one_atom.extend_constraints(vec![constraint.clone()]);
    one_atom.apply_undo(Undo::RestoreMoleculeConstraints(CascadedConstraints {
        removed: vec![
            RemovedConstraint {
                position,
                constraint: constraint.clone(),
            },
            RemovedConstraint {
                position,
                constraint: constraint.clone(),
            },
        ],
        modified: vec![ModifiedConstraint {
            position,
            old: constraint,
            new: Constraint::Or(vec![]),
        }],
    }));
}

#[rstest]
#[case::empty(0)]
#[case::smaller(1)]
#[case::larger(4)]
fn test_molecule_apply_undo_unrelated(#[case] atom_count: usize) {
    let mut source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        ..Default::default()
    });
    let edits = Edits::from_iter([
        Edit::AddAtoms {
            atoms: vec![AtomForm::default()],
        },
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(1))],
            bonds: vec![],
        },
    ]);
    let mut state = ApplicationState::new(&source);
    let mut journal = Vec::new();
    for edit in edits {
        if let Some(undo) = source.apply_edit_with_undo(edit, &mut state).unwrap() {
            journal.push(undo);
        }
    }
    let mut unrelated = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); atom_count],
        ..Default::default()
    });
    for undo in journal.into_iter().rev() {
        unrelated.apply_undo(undo);
    }
}
