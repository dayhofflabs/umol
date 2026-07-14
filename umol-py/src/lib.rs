//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types mirror the Rust API (see `umol-ast`).

use pyo3::prelude::*;

#[cfg(feature = "graph")]
use crate::{
    aromatic::{
        AromaticSystemAst, AromaticSystemConstraintAst, AromaticSystemConstraintKey,
        AromaticSystemConstraintsAst, AromaticSystemConstraintsView, AromaticSystemView,
        AromaticSystemViews,
    },
    atom::{
        AtomAst, AtomConstraintAst, AtomConstraintKey, AtomConstraintsAst, AtomConstraintsView,
        AtomRingSizeCounts, AtomView, AtomViews, ElementAst, IsotopeMassAst, SpinStateAst,
    },
    bond::{
        BondAst, BondConstraintAst, BondConstraintKey, BondConstraintsAst, BondConstraintsView,
        BondRingSizeCounts, BondView, BondViews,
    },
    boolean::BooleanAst,
    constraint::{AromaticValenceAst, MulticenterValenceAst, RingMembershipAst, RingScope},
    dative::{
        DativeBondAst, DativeBondConstraintAst, DativeBondConstraintKey, DativeBondConstraintsAst,
        DativeBondConstraintsView, DativeBondRingSizeCounts, DativeBondView, DativeBondViews,
    },
    electrons::ElectronCountsAst,
    element::Element,
    error::ParseError,
    molecule::MoleculeAst,
    multicenter::{
        MulticenterBondAst, MulticenterBondConstraintAst, MulticenterBondConstraintKey,
        MulticenterBondConstraintsAst, MulticenterBondConstraintsView, MulticenterBondView,
        MulticenterBondViews,
    },
    noncovalent::{
        NoncovalentBondAst, NoncovalentBondConstraintAst, NoncovalentBondConstraintKey,
        NoncovalentBondConstraintsAst, NoncovalentBondConstraintsView, NoncovalentBondKind,
        NoncovalentBondKindAst, NoncovalentBondView, NoncovalentBondViews,
    },
    stereo::{
        CisTransStereo, CisTransStereoAst, FluxionalityAst, LigandPermutation, LigandSymmetryAst,
        Orientation, OrientedLigandPermutation, Permutation, StereoAtomAst,
        StereoAtomConstraintAst, StereoAtomConstraintKey, StereoAtomConstraintsAst,
        StereoAtomConstraintsView, StereoAtomView, StereoAtomViews, StereoBondAst,
        StereoBondConstraintAst, StereoBondConstraintKey, StereoBondConstraintsAst,
        StereoBondConstraintsView, StereoBondView, StereoBondViews, StereoConfigurationAst,
        StereoCosetAst, StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, StereoTerm,
        Stereogenicity, StereogenicityAst, TetrahedralStereo, TetrahedralStereoAst, Topicity,
        TopicityAst, TopicityRelationAst,
    },
    value::{MemOp, RelOp, ValueAst, ValuePredicate, ValueTerm},
};

#[cfg(feature = "graph")]
mod aromatic;
#[cfg(feature = "graph")]
mod atom;
#[cfg(feature = "graph")]
mod bond;
#[cfg(feature = "graph")]
mod boolean;
#[cfg(feature = "graph")]
mod constraint;
#[cfg(feature = "graph")]
mod convert;
#[cfg(feature = "graph")]
mod dative;
#[cfg(feature = "graph")]
mod electrons;
#[cfg(feature = "graph")]
mod element;
#[cfg(feature = "graph")]
mod error;
#[cfg(feature = "graph")]
mod molecule;
#[cfg(feature = "graph")]
mod molecule_constraint;
#[cfg(feature = "graph")]
mod multicenter;
#[cfg(feature = "graph")]
mod noncovalent;
#[cfg(feature = "graph")]
mod stereo;
#[cfg(feature = "graph")]
mod value;

/// The native extension module. Wrapper types are registered here as the
/// binding grows; the graph domain is gated behind the `graph` feature.
#[pymodule]
#[cfg_attr(not(feature = "graph"), allow(unused_variables))]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "graph")]
    {
        module.add_class::<Element>()?;
        module.add("ParseError", module.py().get_type::<ParseError>())?;
        module.add_class::<MoleculeAst>()?;
        module.add_class::<RelOp>()?;
        module.add_class::<MemOp>()?;
        module.add_class::<ValueTerm>()?;
        module.add_class::<ValuePredicate>()?;
        module.add_class::<ValueAst>()?;
        module.add_class::<ElementAst>()?;
        module.add_class::<IsotopeMassAst>()?;
        module.add_class::<SpinStateAst>()?;
        module.add_class::<AtomAst>()?;
        module.add_class::<AtomView>()?;
        module.add_class::<AtomViews>()?;
        module.add_class::<AromaticValenceAst>()?;
        module.add_class::<MulticenterValenceAst>()?;
        module.add_class::<RingScope>()?;
        module.add_class::<RingMembershipAst>()?;
        module.add_class::<AtomConstraintAst>()?;
        module.add_class::<AtomConstraintKey>()?;
        module.add_class::<AtomConstraintsAst>()?;
        module.add_class::<AtomConstraintsView>()?;
        module.add_class::<AtomRingSizeCounts>()?;
        module.add_class::<BondAst>()?;
        module.add_class::<BondConstraintAst>()?;
        module.add_class::<BondConstraintKey>()?;
        module.add_class::<BondConstraintsAst>()?;
        module.add_class::<BondConstraintsView>()?;
        module.add_class::<BondRingSizeCounts>()?;
        module.add_class::<BondView>()?;
        module.add_class::<BondViews>()?;
        module.add_class::<DativeBondAst>()?;
        module.add_class::<DativeBondConstraintAst>()?;
        module.add_class::<DativeBondConstraintKey>()?;
        module.add_class::<DativeBondConstraintsAst>()?;
        module.add_class::<DativeBondConstraintsView>()?;
        module.add_class::<DativeBondRingSizeCounts>()?;
        module.add_class::<DativeBondView>()?;
        module.add_class::<DativeBondViews>()?;
        module.add_class::<ElectronCountsAst>()?;
        module.add_class::<AromaticSystemAst>()?;
        module.add_class::<AromaticSystemConstraintAst>()?;
        module.add_class::<AromaticSystemConstraintKey>()?;
        module.add_class::<AromaticSystemConstraintsAst>()?;
        module.add_class::<AromaticSystemConstraintsView>()?;
        module.add_class::<AromaticSystemView>()?;
        module.add_class::<AromaticSystemViews>()?;
        module.add_class::<MulticenterBondAst>()?;
        module.add_class::<MulticenterBondConstraintAst>()?;
        module.add_class::<MulticenterBondConstraintKey>()?;
        module.add_class::<MulticenterBondConstraintsAst>()?;
        module.add_class::<MulticenterBondConstraintsView>()?;
        module.add_class::<MulticenterBondView>()?;
        module.add_class::<MulticenterBondViews>()?;
        module.add_class::<NoncovalentBondKind>()?;
        module.add_class::<NoncovalentBondKindAst>()?;
        module.add_class::<NoncovalentBondAst>()?;
        module.add_class::<NoncovalentBondConstraintAst>()?;
        module.add_class::<NoncovalentBondConstraintKey>()?;
        module.add_class::<NoncovalentBondConstraintsAst>()?;
        module.add_class::<NoncovalentBondConstraintsView>()?;
        module.add_class::<NoncovalentBondView>()?;
        module.add_class::<NoncovalentBondViews>()?;
        module.add_class::<BooleanAst>()?;
        module.add_class::<Permutation>()?;
        module.add_class::<StereoTerm>()?;
        module.add_class::<StereoCosetAst>()?;
        module.add_class::<TetrahedralStereoAst>()?;
        module.add_class::<TetrahedralStereo>()?;
        module.add_class::<CisTransStereoAst>()?;
        module.add_class::<CisTransStereo>()?;
        module.add_class::<StereoKind>()?;
        module.add_class::<StereoLigandKind>()?;
        module.add_class::<StereoLigand>()?;
        module.add_class::<Topicity>()?;
        module.add_class::<Stereogenicity>()?;
        module.add_class::<StereoConfigurationAst>()?;
        module.add_class::<Orientation>()?;
        module.add_class::<LigandPermutation>()?;
        module.add_class::<OrientedLigandPermutation>()?;
        module.add_class::<StereoLigandPair>()?;
        module.add_class::<TopicityRelationAst>()?;
        module.add_class::<StereogenicityAst>()?;
        module.add_class::<LigandSymmetryAst>()?;
        module.add_class::<FluxionalityAst>()?;
        module.add_class::<TopicityAst>()?;
        module.add_class::<StereoAtomConstraintKey>()?;
        module.add_class::<StereoAtomConstraintAst>()?;
        module.add_class::<StereoAtomConstraintsAst>()?;
        module.add_class::<StereoBondConstraintKey>()?;
        module.add_class::<StereoBondConstraintAst>()?;
        module.add_class::<StereoBondConstraintsAst>()?;
        module.add_class::<StereoAtomAst>()?;
        module.add_class::<StereoBondAst>()?;
        module.add_class::<StereoAtomConstraintsView>()?;
        module.add_class::<StereoBondConstraintsView>()?;
        module.add_class::<StereoAtomView>()?;
        module.add_class::<StereoBondView>()?;
        module.add_class::<StereoAtomViews>()?;
        module.add_class::<StereoBondViews>()?;
    }
    Ok(())
}
