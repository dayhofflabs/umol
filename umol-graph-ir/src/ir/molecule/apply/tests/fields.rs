use rstest::*;
use umol_chem::element::Element;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::{AtomForm, ElementForm, IsotopeMassForm};
use crate::ir::bond::BondForm;
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange,
    BondHandle, DativeBondFieldChange, DativeBondHandle, Edit, Edits, MulticenterBondFieldChange,
    MulticenterBondHandle, NoncovalentBondFieldChange, NoncovalentBondHandle,
    StereoAtomFieldChange, StereoAtomHandle, StereoBondFieldChange, StereoBondHandle, Undo,
};
use crate::ir::electrons::ElectronCountsForm;
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
};
use crate::ir::traits::Normalize;
use crate::ir::{Molecule, MoleculeApplyError, RollbackMethod, Transaction, TransactionError};

mod aromatic;
mod atom;
mod bond;
mod dative;
mod multicenter;
mod noncovalent;
mod stereo;

#[rstest]
#[case(Edit::ModifyAtomElements { changes: vec![] })]
#[case(Edit::ModifyAtomIsotopeMasses { changes: vec![] })]
#[case(Edit::ModifyAtomCharges { changes: vec![] })]
#[case(Edit::ModifyAtomImplicitHydrogens { changes: vec![] })]
#[case(Edit::ModifyAtomLonePairs { changes: vec![] })]
#[case(Edit::ModifyAtomUnpairedElectrons { changes: vec![] })]
#[case(Edit::ModifyBondOrders { changes: vec![] })]
#[case(Edit::ModifyBondCharges { changes: vec![] })]
#[case(Edit::ModifyBondUnpairedElectrons { changes: vec![] })]
#[case(Edit::ModifyDativeBondOrders { changes: vec![] })]
#[case(Edit::ModifyAromaticSystemElectrons { changes: vec![] })]
#[case(Edit::ModifyAromaticSystemCharges { changes: vec![] })]
#[case(Edit::ModifyAromaticSystemUnpairedElectrons { changes: vec![] })]
#[case(Edit::ModifyMulticenterBondElectrons { changes: vec![] })]
#[case(Edit::ModifyMulticenterBondCharges { changes: vec![] })]
#[case(Edit::ModifyMulticenterBondUnpairedElectrons { changes: vec![] })]
#[case(Edit::ModifyNoncovalentBondKinds { changes: vec![] })]
#[case(Edit::ModifyStereoAtomConfigurations { changes: vec![] })]
#[case(Edit::ModifyStereoBondConfigurations { changes: vec![] })]
fn test_molecule_apply_field_batches_identity(#[case] edit: Edit) {
    let source = Molecule::new();
    assert_eq!(
        source
            .clone()
            .apply(Edits::from_iter([edit.clone()]))
            .unwrap(),
        source
    );
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    molecule.apply_undo(undo);
    assert_eq!(molecule, source);
}
