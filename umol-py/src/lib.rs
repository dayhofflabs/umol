//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types correspond to the Rust API (see `umol-ast`).

use pyo3::prelude::*;

#[cfg(feature = "graph")]
use crate::{
    algorithm::{
        AutomorphismAlgorithm, CommonSubgraphEnumerationAlgorithm, ConnectedComponentsAlgorithm,
        MaximumIndependentSetAlgorithm, RelevantCycleEnumerationAlgorithm,
        SimpleCycleEnumerationAlgorithm, SubgraphEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm, SubstructureMatchAlgorithm,
    },
    aromatic::{AromaticSystemAst, AromaticSystemView, AromaticSystemViews},
    atom::{AtomAst, AtomView, AtomViews, ElementAst, IsotopeMass, IsotopeMassAst},
    bond::{BondAst, BondView, BondViews},
    boolean::BooleanAst,
    constraint::{
        aromatic::{
            AromaticSystemConstraintAst, AromaticSystemConstraintKey, AromaticSystemConstraintsAst,
            AromaticSystemConstraintsView,
        },
        atom::{
            AromaticValence, AromaticValenceAst, AtomConstraintAst, AtomConstraintKey,
            AtomConstraintsAst, AtomConstraintsView, AtomRingSizeCounts, MulticenterValence,
            MulticenterValenceAst,
        },
        bond::{
            BondConstraintAst, BondConstraintKey, BondConstraintsAst, BondConstraintsView,
            BondRingSizeCounts,
        },
        dative::{
            DativeBondConstraintAst, DativeBondConstraintKey, DativeBondConstraintsAst,
            DativeBondConstraintsView, DativeBondRingSizeCounts,
        },
        molecule::{
            Constraint, Constraints, ConstraintsView, MoleculeConstraint, RelationalConstraint,
            SubPatternAnchor,
        },
        multicenter::{
            MulticenterBondConstraintAst, MulticenterBondConstraintKey,
            MulticenterBondConstraintsAst, MulticenterBondConstraintsView,
        },
        noncovalent::{
            NoncovalentBondConstraintAst, NoncovalentBondConstraintKey,
            NoncovalentBondConstraintsAst, NoncovalentBondConstraintsView,
        },
        ring::{RingMembershipAst, RingScope},
        stereo::{
            FluxionalityAst, LigandSymmetryAst, StereoAtomConstraintAst, StereoAtomConstraintKey,
            StereoAtomConstraintsAst, StereoAtomConstraintsView, StereoBondConstraintAst,
            StereoBondConstraintKey, StereoBondConstraintsAst, StereoBondConstraintsView,
            StereogenicityAst, TopicityAst, TopicityRelationAst,
        },
    },
    correspondence::{Correspondence, MoleculeCorrespondence},
    dative::{DativeBondAst, DativeBondView, DativeBondViews},
    defaults::{MoleculeDefaults, ReactionDefaults},
    delta::{
        AromaticSystemDelta, AromaticSystemFieldChange, AtomDelta, AtomFieldChange, BondDelta,
        BondFieldChange, ConstraintDelta, DativeBondDelta, DativeBondFieldChange, Delta, Deltas,
        MulticenterBondDelta, MulticenterBondFieldChange, NoncovalentBondDelta,
        NoncovalentBondFieldChange, StereoAtomDelta, StereoAtomFieldChange, StereoBondDelta,
        StereoBondFieldChange,
    },
    electrons::ElectronCountsAst,
    element::Element,
    error::{
        ContradictionError, InvalidStructureError, MetadataError, ModelConversionError, ParseError,
        UnderdeterminedError,
    },
    fingerprint::config::{
        EcfpHashScheme, HashedFingerprintConfig, PatternFingerprintConfig,
        ReactionCombinedFingerprintConfig, RefinementRounds, StructuralFingerprintConfig,
        WlHashScheme,
    },
    fingerprint::reaction::{
        ReactionCombinedFingerprint, ReactionSide, RoleTaggedHashedFeatureSet,
        SignedHashedFeatureSet,
    },
    fingerprint::value::{BitFp, CountedHashedFeatureSet, HashedFeatureSet, StructuralFeatureSet},
    metadata::{Entity, MoleculeMetadata, ReactionMetadata},
    model::{
        aromaticity::{AromaticityConfig, AromaticityModel, RingLimits},
        stereo::{StereoKindModel, StereoModel},
        valence::{AtomTypeRegistry, ValenceEntry, ValenceModel, ValenceTable},
        ChemistryModel, ElementScope,
    },
    molecule::MoleculeAst,
    multicenter::{MulticenterBondAst, MulticenterBondView, MulticenterBondViews},
    noncovalent::{
        NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst, NoncovalentBondView,
        NoncovalentBondViews,
    },
    reaction::{
        ReactionApplicationConfig, ReactionAst, ReactionCompositionConfig, ReactionDerivation,
    },
    resolve::{
        AromaticBondConstraintMismatchPolicy, AromaticityFailurePolicy, AromaticityMismatchPolicy,
        AromaticityResolveConfig, ResolveConfig, StereoFailurePolicy, StereoMismatchPolicy,
        StereoResolveConfig,
    },
    ring::RingConfig,
    smiles::{SmilesIoConfig, SmilesSyntaxFlags},
    spin::{SpinState, UnpairedElectrons, UnpairedElectronsAst},
    stereo::{
        CisTransConfiguration, CisTransStereo, CisTransStereoAst, LigandPermutation, Orientation,
        OrientedLigandPermutation, Permutation, StereoAtomAst, StereoAtomView, StereoAtomViews,
        StereoBondAst, StereoBondView, StereoBondViews, StereoConfigurationAst, StereoCoset,
        StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, StereoTerm, Stereogenicity,
        TetrahedralConfiguration, TetrahedralStereo, TetrahedralStereoAst, Topicity,
    },
    substructure::SubstructureSearchConfig,
    value::{MemOp, RelOp, ValueAst, ValuePredicate, ValueTerm},
};

#[cfg(feature = "graph")]
mod algorithm;
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
mod correspondence;
#[cfg(feature = "graph")]
mod dative;
#[cfg(feature = "graph")]
mod defaults;
#[cfg(feature = "graph")]
mod delta;
#[cfg(feature = "graph")]
mod electrons;
#[cfg(feature = "graph")]
mod element;
#[cfg(feature = "graph")]
mod error;
#[cfg(feature = "graph")]
mod fingerprint;
#[cfg(feature = "graph")]
mod lattice;
#[cfg(feature = "graph")]
mod metadata;
#[cfg(feature = "graph")]
mod model;
#[cfg(feature = "graph")]
mod molecule;
#[cfg(feature = "graph")]
mod multicenter;
#[cfg(feature = "graph")]
mod noncovalent;
#[cfg(feature = "graph")]
mod reaction;
#[cfg(feature = "graph")]
mod resolve;
#[cfg(feature = "graph")]
mod ring;
#[cfg(feature = "graph")]
mod smiles;
#[cfg(feature = "graph")]
mod spin;
#[cfg(feature = "graph")]
mod stereo;
#[cfg(feature = "graph")]
mod substructure;
#[cfg(feature = "graph")]
mod value;

/// The native extension module. Wrapper types are registered here as the
/// binding grows; the graph domain is gated behind the `graph` feature.
#[pymodule]
#[cfg_attr(not(feature = "graph"), allow(unused_variables))]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "graph")]
    {
        module.add_class::<AutomorphismAlgorithm>()?;
        module.add_class::<CommonSubgraphEnumerationAlgorithm>()?;
        module.add_class::<ConnectedComponentsAlgorithm>()?;
        module.add_class::<SimpleCycleEnumerationAlgorithm>()?;
        module.add_class::<RelevantCycleEnumerationAlgorithm>()?;
        module.add_class::<MaximumIndependentSetAlgorithm>()?;
        module.add_class::<SubgraphEnumerationAlgorithm>()?;
        module.add_class::<SubgraphIsomorphismAlgorithm>()?;
        module.add_class::<SubstructureMatchAlgorithm>()?;
        module.add_class::<Element>()?;
        module.add(
            "ContradictionError",
            module.py().get_type::<ContradictionError>(),
        )?;
        module.add(
            "InvalidStructureError",
            module.py().get_type::<InvalidStructureError>(),
        )?;
        module.add(
            "ModelConversionError",
            module.py().get_type::<ModelConversionError>(),
        )?;
        module.add("MetadataError", module.py().get_type::<MetadataError>())?;
        module.add("ParseError", module.py().get_type::<ParseError>())?;
        module.add(
            "UnderdeterminedError",
            module.py().get_type::<UnderdeterminedError>(),
        )?;
        module.add_class::<MoleculeAst>()?;
        module.add_class::<MoleculeDefaults>()?;
        module.add_class::<ReactionDefaults>()?;
        module.add_class::<AtomTypeRegistry>()?;
        module.add_class::<ValenceEntry>()?;
        module.add_class::<ValenceModel>()?;
        module.add_class::<ValenceTable>()?;
        module.add_class::<ElementScope>()?;
        module.add_class::<RingLimits>()?;
        module.add_class::<AromaticityConfig>()?;
        module.add_class::<AromaticityModel>()?;
        module.add_class::<AromaticBondConstraintMismatchPolicy>()?;
        module.add_class::<AromaticityFailurePolicy>()?;
        module.add_class::<AromaticityMismatchPolicy>()?;
        module.add_class::<StereoFailurePolicy>()?;
        module.add_class::<StereoMismatchPolicy>()?;
        module.add_class::<StereoKindModel>()?;
        module.add_class::<StereoModel>()?;
        module.add_class::<ChemistryModel>()?;
        module.add_class::<AromaticityResolveConfig>()?;
        module.add_class::<StereoResolveConfig>()?;
        module.add_class::<ResolveConfig>()?;
        module.add_class::<AromaticSystemDelta>()?;
        module.add_class::<AtomDelta>()?;
        module.add_class::<AtomFieldChange>()?;
        module.add_class::<BondDelta>()?;
        module.add_class::<BondFieldChange>()?;
        module.add_class::<DativeBondDelta>()?;
        module.add_class::<DativeBondFieldChange>()?;
        module.add_class::<AromaticSystemFieldChange>()?;
        module.add_class::<MulticenterBondDelta>()?;
        module.add_class::<MulticenterBondFieldChange>()?;
        module.add_class::<NoncovalentBondDelta>()?;
        module.add_class::<NoncovalentBondFieldChange>()?;
        module.add_class::<StereoAtomDelta>()?;
        module.add_class::<StereoAtomFieldChange>()?;
        module.add_class::<StereoBondDelta>()?;
        module.add_class::<StereoBondFieldChange>()?;
        module.add_class::<ConstraintDelta>()?;
        module.add_class::<Delta>()?;
        module.add_class::<Deltas>()?;
        module.add_class::<ReactionCompositionConfig>()?;
        module.add_class::<ReactionApplicationConfig>()?;
        module.add_class::<ReactionAst>()?;
        module.add_class::<ReactionDerivation>()?;
        module.add_class::<RingConfig>()?;
        module.add_class::<SmilesIoConfig>()?;
        module.add_class::<SmilesSyntaxFlags>()?;
        module.add_class::<RefinementRounds>()?;
        module.add_class::<WlHashScheme>()?;
        module.add_class::<EcfpHashScheme>()?;
        module.add_class::<HashedFingerprintConfig>()?;
        module.add_class::<PatternFingerprintConfig>()?;
        module.add_class::<StructuralFingerprintConfig>()?;
        module.add_class::<ReactionCombinedFingerprintConfig>()?;
        module.add_class::<HashedFeatureSet>()?;
        module.add_class::<CountedHashedFeatureSet>()?;
        module.add_class::<BitFp>()?;
        module.add_class::<StructuralFeatureSet>()?;
        module.add_class::<ReactionSide>()?;
        module.add_class::<SignedHashedFeatureSet>()?;
        module.add_class::<RoleTaggedHashedFeatureSet>()?;
        module.add_class::<ReactionCombinedFingerprint>()?;
        module.add_class::<Correspondence>()?;
        module.add_class::<MoleculeCorrespondence>()?;
        module.add_class::<Entity>()?;
        module.add_class::<MoleculeMetadata>()?;
        module.add_class::<ReactionMetadata>()?;
        module.add_class::<SubstructureSearchConfig>()?;
        module.add_class::<SubPatternAnchor>()?;
        module.add_class::<RelationalConstraint>()?;
        module.add_class::<MoleculeConstraint>()?;
        module.add_class::<Constraint>()?;
        module.add_class::<Constraints>()?;
        module.add_class::<ConstraintsView>()?;
        module.add_class::<RelOp>()?;
        module.add_class::<MemOp>()?;
        module.add_class::<ValueTerm>()?;
        module.add_class::<ValuePredicate>()?;
        module.add_class::<ValueAst>()?;
        module.add_class::<ElementAst>()?;
        module.add_class::<IsotopeMass>()?;
        module.add_class::<IsotopeMassAst>()?;
        module.add_class::<UnpairedElectrons>()?;
        module.add_class::<SpinState>()?;
        module.add_class::<UnpairedElectronsAst>()?;
        module.add_class::<AtomAst>()?;
        module.add_class::<AtomView>()?;
        module.add_class::<AtomViews>()?;
        module.add_class::<AromaticValence>()?;
        module.add_class::<AromaticValenceAst>()?;
        module.add_class::<MulticenterValence>()?;
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
        module.add_class::<StereoCoset>()?;
        module.add_class::<TetrahedralStereoAst>()?;
        module.add_class::<TetrahedralStereo>()?;
        module.add_class::<TetrahedralConfiguration>()?;
        module.add_class::<CisTransStereoAst>()?;
        module.add_class::<CisTransStereo>()?;
        module.add_class::<CisTransConfiguration>()?;
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
