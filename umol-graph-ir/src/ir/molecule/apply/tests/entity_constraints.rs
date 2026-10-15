use rstest::*;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AtomConstraintForm,
    AtomConstraintKey, BondConstraintForm, BondConstraintKey, DativeBondConstraintForm,
    DativeBondConstraintKey, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoBondConstraintForm, StereoBondConstraintKey, StereogenicityForm,
};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondHandle, DativeBondHandle, Edit, Edits,
    MulticenterBondHandle, NoncovalentBondHandle, StereoAtomHandle, StereoBondHandle, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::NoncovalentBondForm;
use crate::ir::num::NumForm;
use crate::ir::stereo::{CisTransStereo, StereoKind, Stereogenicity};
use crate::ir::traits::Normalize;
use crate::ir::{
    Molecule, MoleculeApplyError, RingScope, RollbackMethod, Transaction, TransactionError,
};

mod aromatic;
mod atom;
mod bond;
mod dative;
mod multicenter;
mod noncovalent;
mod stereo;

#[rstest]
#[case::atom(Edit::ModifyAtomConstraint {
    id: AtomHandle::Id(AtomId(0)), old: None, new: None,
}, EntityKind::Atom)]
#[case::bond(Edit::ModifyBondConstraint {
    id: BondHandle::Id(BondId(0)), old: None, new: None,
}, EntityKind::Bond)]
#[case::dative_bond(Edit::ModifyDativeBondConstraint {
    id: DativeBondHandle::Id(DativeBondId(0)), old: None, new: None,
}, EntityKind::DativeBond)]
#[case::aromatic_system(Edit::ModifyAromaticSystemConstraint {
    id: AromaticSystemHandle::Id(AromaticSystemId(0)), old: None, new: None,
}, EntityKind::AromaticSystem)]
#[case::multicenter_bond(Edit::ModifyMulticenterBondConstraint {
    id: MulticenterBondHandle::Id(MulticenterBondId(0)), old: None, new: None,
}, EntityKind::MulticenterBond)]
#[case::noncovalent_bond(Edit::ModifyNoncovalentBondConstraint {
    id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), old: None, new: None,
}, EntityKind::NoncovalentBond)]
#[case::stereo_atom(Edit::ModifyStereoAtomConstraint {
    id: StereoAtomHandle::Id(StereoAtomId(0)), kind: None, old: None, new: None,
}, EntityKind::StereoAtom)]
#[case::stereo_bond(Edit::ModifyStereoBondConstraint {
    id: StereoBondHandle::Id(StereoBondId(0)), kind: None, old: None, new: None,
}, EntityKind::StereoBond)]
fn test_molecule_apply_edit_constraint_target_error(
    #[case] edit: Edit,
    #[case] kind: EntityKind,
    #[values(false, true)] journaled: bool,
) {
    let mut molecule = Molecule::default();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(
        result,
        Err(TransactionError::HandleOutOfRange {
            kind,
            index: 0,
            count: 0
        })
    );
    assert_eq!(molecule, Molecule::default());
}
