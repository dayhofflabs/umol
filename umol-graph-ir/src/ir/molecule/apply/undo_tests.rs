//! Undo replay with manipulated entries and unrelated receiver states.

use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{Compaction, GraphCompaction};

use super::ApplicationState;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::compact::MoleculeCompaction;
use crate::ir::constraint::{AtomConstraintForm, Constraint};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AddedAtom, AromaticSystemFieldChange, AtomFieldChange, AtomHandle, BondFieldChange,
    CascadedConstraints, DativeBondFieldChange, Edit, Edits, MulticenterBondFieldChange,
    NoncovalentBondFieldChange, RemovedAromaticSystem, RemovedAtom, RemovedConstraint,
    RemovedOverlays, StereoAtomFieldChange, StereoBondFieldChange, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::{Molecule, MoleculeEntries};
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
};
use crate::ir::ModifiedConstraint;

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
    editor.build()
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
fn test_molecule_apply_undo_reconstruction_entry(mut empty: Molecule) {
    empty.apply_undo(Undo::RestoreRemovedAromaticSystems {
        removed: vec![RemovedAromaticSystem {
            id: AromaticSystemId(1),
            atoms: Vec::new(),
            attributes: AromaticSystemForm::default(),
        }],
        undo_compaction: MoleculeCompaction::empty().undo_compaction(),
        cascade: CascadedConstraints::default(),
    });
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
    one_atom.apply_undo(Undo::RestoreRemovedTopology {
        atoms: Vec::new(),
        bonds: Vec::new(),
        overlays: RemovedOverlays::default(),
        undo_compaction: compaction.undo_compaction(),
        compaction,
        cascade: CascadedConstraints::default(),
    });
}

#[rstest]
fn test_molecule_apply_undo_compaction_dimension(mut one_atom: Molecule) {
    let compaction = MoleculeCompaction::empty();
    one_atom.apply_undo(Undo::RestoreRemovedTopology {
        atoms: vec![RemovedAtom {
            id: AtomId(0),
            attributes: AtomForm::from_element(Element::N),
        }],
        bonds: Vec::new(),
        overlays: RemovedOverlays::default(),
        undo_compaction: compaction.undo_compaction(),
        compaction,
        cascade: CascadedConstraints::default(),
    });
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
