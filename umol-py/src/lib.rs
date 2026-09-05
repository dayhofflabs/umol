//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types correspond to the Rust API (see `umol-graph-ir`).

use pyo3::prelude::*;

#[cfg(feature = "depiction")]
use crate::depict::{DepictConfig, Depiction, MoleculeLayoutAlgorithm};
#[cfg(feature = "graph")]
use crate::{
    algorithm::{
        AutomorphismAlgorithm, CommonSubgraphEnumerationAlgorithm, ConnectedComponentsAlgorithm,
        MaximumIndependentSetAlgorithm, RelevantCycleEnumerationAlgorithm,
        SimpleCycleEnumerationAlgorithm, SubgraphEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm, SubstructureMatchAlgorithm,
    },
    aromatic::{AromaticSystemForm, AromaticSystemUpdate, AromaticSystemView, AromaticSystemViews},
    atom::{AtomForm, AtomUpdate, AtomView, AtomViews, ElementForm, IsotopeMass, IsotopeMassForm},
    bond::{BondForm, BondUpdate, BondView, BondViews},
    boolean::BooleanForm,
    canonicalize::CanonicalizeConfig,
    constraint::{
        aromatic::{
            AromaticSystemConstraintForm, AromaticSystemConstraintKey,
            AromaticSystemConstraintsForm, AromaticSystemConstraintsView,
        },
        atom::{
            AromaticValence, AromaticValenceForm, AtomConstraintForm, AtomConstraintKey,
            AtomConstraintsForm, AtomConstraintsView, AtomRingSizeCounts, MulticenterValence,
            MulticenterValenceForm,
        },
        bond::{
            BondConstraintForm, BondConstraintKey, BondConstraintsForm, BondConstraintsView,
            BondRingSizeCounts,
        },
        dative::{
            DativeBondConstraintForm, DativeBondConstraintKey, DativeBondConstraintsForm,
            DativeBondConstraintsView, DativeBondRingSizeCounts,
        },
        molecule::{
            Constraint, Constraints, ConstraintsView, MoleculeConstraint, RelationalConstraint,
        },
        multicenter::{
            MulticenterBondConstraintForm, MulticenterBondConstraintKey,
            MulticenterBondConstraintsForm, MulticenterBondConstraintsView,
        },
        noncovalent::{
            NoncovalentBondConstraintForm, NoncovalentBondConstraintKey,
            NoncovalentBondConstraintsForm, NoncovalentBondConstraintsView,
        },
        ring::{RingMembershipForm, RingScope},
        stereo::{
            FluxionalityForm, LigandSymmetryForm, StereoAtomConstraintForm,
            StereoAtomConstraintKey, StereoAtomConstraintsForm, StereoAtomConstraintsView,
            StereoBondConstraintForm, StereoBondConstraintKey, StereoBondConstraintsForm,
            StereoBondConstraintsView, StereogenicityForm, TopicityForm, TopicityRelationForm,
        },
    },
    correspondence::{Correspondence, MoleculeCorrespondence},
    dative::{DativeBondForm, DativeBondUpdate, DativeBondView, DativeBondViews},
    defaults::{MoleculeDefaults, ReactionDefaults},
    delta::{
        AromaticSystemDelta, AromaticSystemFieldChange, AtomDelta, AtomFieldChange, BondDelta,
        BondFieldChange, ConstraintDelta, DativeBondDelta, DativeBondFieldChange, Delta, Deltas,
        MulticenterBondDelta, MulticenterBondFieldChange, NoncovalentBondDelta,
        NoncovalentBondFieldChange, StereoAtomDelta, StereoAtomFieldChange, StereoBondDelta,
        StereoBondFieldChange,
    },
    edit::{ConstraintEdit, Edit, Edits, New},
    electrons::ElectronCountsForm,
    element::Element,
    error::{
        ContradictionError, InvalidStructureError, MetadataError, ModelConversionError, ParseError,
        TransactionError, UnderdeterminedError,
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
        aromaticity::{
            AromaticityConfig, AromaticityModel, AromaticityRule, AromaticityTieBreak, RingLimits,
        },
        connectivity::ConnectivityModel,
        stereo::{StereoKindModel, StereoModel},
        valence::{
            AtomTypeRegistry, ValenceCandidateSource, ValenceEntry, ValenceModel, ValenceTable,
            ValenceTieBreak,
        },
        ChemistryModel, ElementScope,
    },
    molecule::Molecule,
    multicenter::{
        MulticenterBondForm, MulticenterBondUpdate, MulticenterBondView, MulticenterBondViews,
    },
    noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm, NoncovalentBondUpdate,
        NoncovalentBondView, NoncovalentBondViews,
    },
    num::{ArithExpr, MemOp, NumForm, PredExpr, RelOp},
    reaction::{Reaction, ReactionApplicationConfig, ReactionCompositionConfig},
    reaction_span::ReactionSpan,
    remap::{MoleculeRemapping, Remapping},
    resolve::{
        AromaticBondConstraintMismatchPolicy, AromaticityFailurePolicy, AromaticityMismatchPolicy,
        AromaticityResolveConfig, AtomCompletions, ResolveConfig, ResolveContradiction,
        ResolveReport, Solution, StereoFailurePolicy, StereoMismatchPolicy, StereoResolveConfig,
    },
    ring::RingConfig,
    smiles::{SmilesIoConfig, SmilesSyntaxFlags},
    spin::{SpinState, UnpairedElectrons, UnpairedElectronsForm, UnpairedElectronsUpdate},
    stereo::{
        CisTransConfiguration, CisTransStereo, CisTransStereoForm, LigandPermutation, Orientation,
        OrientedLigandPermutation, Permutation, StereoAtomForm, StereoAtomUpdate, StereoAtomView,
        StereoAtomViews, StereoBondForm, StereoBondUpdate, StereoBondView, StereoBondViews,
        StereoConfigurationForm, StereoConfigurationUpdate, StereoCoset, StereoKind, StereoLigand,
        StereoLigandKind, StereoLigandPair, StereoTerm, Stereogenicity, TetrahedralConfiguration,
        TetrahedralStereo, TetrahedralStereoForm, Topicity,
    },
    substructure::SubstructureSearchConfig,
    transaction::{MoleculeEditor, Transaction},
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
mod canonicalize;
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
#[cfg(feature = "depiction")]
mod depict;
#[cfg(feature = "graph")]
mod edit;
#[cfg(feature = "graph")]
mod electrons;
#[cfg(feature = "graph")]
mod element;
#[cfg(feature = "graph")]
mod entity;
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
mod num;
#[cfg(feature = "graph")]
mod reaction;
#[cfg(feature = "graph")]
mod reaction_span;
mod remap;
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
mod transaction;

/// The native extension module. Wrapper types are registered here as the
/// binding grows; the graph domain is gated behind the `graph` feature.
#[pymodule]
#[cfg_attr(not(feature = "graph"), allow(unused_variables))]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "graph")]
    {
        #[cfg(feature = "depiction")]
        module.add_class::<MoleculeLayoutAlgorithm>()?;
        #[cfg(feature = "depiction")]
        module.add_class::<DepictConfig>()?;
        #[cfg(feature = "depiction")]
        module.add_class::<Depiction>()?;
        module.add_class::<AutomorphismAlgorithm>()?;
        module.add_class::<CanonicalizeConfig>()?;
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
            "TransactionError",
            module.py().get_type::<TransactionError>(),
        )?;
        module.add(
            "UnderdeterminedError",
            module.py().get_type::<UnderdeterminedError>(),
        )?;
        module.add_class::<Molecule>()?;
        module.add_class::<MoleculeEditor>()?;
        module.add_class::<Transaction>()?;
        module.add_class::<MoleculeDefaults>()?;
        module.add_class::<ReactionDefaults>()?;
        module.add_class::<AtomTypeRegistry>()?;
        module.add_class::<ValenceEntry>()?;
        module.add_class::<ValenceCandidateSource>()?;
        module.add_class::<ValenceModel>()?;
        module.add_class::<ValenceTieBreak>()?;
        module.add_class::<ValenceTable>()?;
        module.add_class::<ElementScope>()?;
        module.add_class::<RingLimits>()?;
        module.add_class::<AromaticityConfig>()?;
        module.add_class::<AromaticityModel>()?;
        module.add_class::<AromaticityRule>()?;
        module.add_class::<AromaticityTieBreak>()?;
        module.add_class::<AromaticBondConstraintMismatchPolicy>()?;
        module.add_class::<AromaticityFailurePolicy>()?;
        module.add_class::<AromaticityMismatchPolicy>()?;
        module.add_class::<StereoFailurePolicy>()?;
        module.add_class::<StereoMismatchPolicy>()?;
        module.add_class::<StereoKindModel>()?;
        module.add_class::<StereoModel>()?;
        module.add_class::<ConnectivityModel>()?;
        module.add_class::<ChemistryModel>()?;
        module.add_class::<AromaticityResolveConfig>()?;
        module.add_class::<StereoResolveConfig>()?;
        module.add_class::<ResolveConfig>()?;
        module.add_class::<ResolveContradiction>()?;
        module.add_class::<ResolveReport>()?;
        module.add_class::<Solution>()?;
        module.add_class::<AtomCompletions>()?;
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
        module.add_class::<New>()?;
        module.add_class::<ConstraintEdit>()?;
        module.add_class::<Edit>()?;
        module.add_class::<Edits>()?;
        module.add_class::<ReactionCompositionConfig>()?;
        module.add_class::<ReactionApplicationConfig>()?;
        module.add_class::<Reaction>()?;
        module.add_class::<ReactionSpan>()?;
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
        module.add_class::<Remapping>()?;
        module.add_class::<MoleculeRemapping>()?;
        module.add_class::<MoleculeCorrespondence>()?;
        module.add_class::<Entity>()?;
        module.add_class::<MoleculeMetadata>()?;
        module.add_class::<ReactionMetadata>()?;
        module.add_class::<SubstructureSearchConfig>()?;
        module.add_class::<RelationalConstraint>()?;
        module.add_class::<MoleculeConstraint>()?;
        module.add_class::<Constraint>()?;
        module.add_class::<Constraints>()?;
        module.add_class::<ConstraintsView>()?;
        module.add_class::<RelOp>()?;
        module.add_class::<MemOp>()?;
        module.add_class::<ArithExpr>()?;
        module.add_class::<PredExpr>()?;
        module.add_class::<NumForm>()?;
        module.add_class::<ElementForm>()?;
        module.add_class::<IsotopeMass>()?;
        module.add_class::<IsotopeMassForm>()?;
        module.add_class::<UnpairedElectrons>()?;
        module.add_class::<SpinState>()?;
        module.add_class::<UnpairedElectronsForm>()?;
        module.add_class::<UnpairedElectronsUpdate>()?;
        module.add_class::<AtomForm>()?;
        module.add_class::<AtomUpdate>()?;
        module.add_class::<AtomView>()?;
        module.add_class::<AtomViews>()?;
        module.add_class::<AromaticValence>()?;
        module.add_class::<AromaticValenceForm>()?;
        module.add_class::<MulticenterValence>()?;
        module.add_class::<MulticenterValenceForm>()?;
        module.add_class::<RingScope>()?;
        module.add_class::<RingMembershipForm>()?;
        module.add_class::<AtomConstraintForm>()?;
        module.add_class::<AtomConstraintKey>()?;
        module.add_class::<AtomConstraintsForm>()?;
        module.add_class::<AtomConstraintsView>()?;
        module.add_class::<AtomRingSizeCounts>()?;
        module.add_class::<BondForm>()?;
        module.add_class::<BondUpdate>()?;
        module.add_class::<BondConstraintForm>()?;
        module.add_class::<BondConstraintKey>()?;
        module.add_class::<BondConstraintsForm>()?;
        module.add_class::<BondConstraintsView>()?;
        module.add_class::<BondRingSizeCounts>()?;
        module.add_class::<BondView>()?;
        module.add_class::<BondViews>()?;
        module.add_class::<DativeBondForm>()?;
        module.add_class::<DativeBondUpdate>()?;
        module.add_class::<DativeBondConstraintForm>()?;
        module.add_class::<DativeBondConstraintKey>()?;
        module.add_class::<DativeBondConstraintsForm>()?;
        module.add_class::<DativeBondConstraintsView>()?;
        module.add_class::<DativeBondRingSizeCounts>()?;
        module.add_class::<DativeBondView>()?;
        module.add_class::<DativeBondViews>()?;
        module.add_class::<ElectronCountsForm>()?;
        module.add_class::<AromaticSystemForm>()?;
        module.add_class::<AromaticSystemUpdate>()?;
        module.add_class::<AromaticSystemConstraintForm>()?;
        module.add_class::<AromaticSystemConstraintKey>()?;
        module.add_class::<AromaticSystemConstraintsForm>()?;
        module.add_class::<AromaticSystemConstraintsView>()?;
        module.add_class::<AromaticSystemView>()?;
        module.add_class::<AromaticSystemViews>()?;
        module.add_class::<MulticenterBondForm>()?;
        module.add_class::<MulticenterBondUpdate>()?;
        module.add_class::<MulticenterBondConstraintForm>()?;
        module.add_class::<MulticenterBondConstraintKey>()?;
        module.add_class::<MulticenterBondConstraintsForm>()?;
        module.add_class::<MulticenterBondConstraintsView>()?;
        module.add_class::<MulticenterBondView>()?;
        module.add_class::<MulticenterBondViews>()?;
        module.add_class::<NoncovalentBondKind>()?;
        module.add_class::<NoncovalentBondKindForm>()?;
        module.add_class::<NoncovalentBondForm>()?;
        module.add_class::<NoncovalentBondUpdate>()?;
        module.add_class::<NoncovalentBondConstraintForm>()?;
        module.add_class::<NoncovalentBondConstraintKey>()?;
        module.add_class::<NoncovalentBondConstraintsForm>()?;
        module.add_class::<NoncovalentBondConstraintsView>()?;
        module.add_class::<NoncovalentBondView>()?;
        module.add_class::<NoncovalentBondViews>()?;
        module.add_class::<BooleanForm>()?;
        module.add_class::<Permutation>()?;
        module.add_class::<StereoTerm>()?;
        module.add_class::<StereoCoset>()?;
        module.add_class::<TetrahedralStereoForm>()?;
        module.add_class::<TetrahedralStereo>()?;
        module.add_class::<TetrahedralConfiguration>()?;
        module.add_class::<CisTransStereoForm>()?;
        module.add_class::<CisTransStereo>()?;
        module.add_class::<CisTransConfiguration>()?;
        module.add_class::<StereoKind>()?;
        module.add_class::<StereoLigandKind>()?;
        module.add_class::<StereoLigand>()?;
        module.add_class::<Topicity>()?;
        module.add_class::<Stereogenicity>()?;
        module.add_class::<StereoConfigurationForm>()?;
        module.add_class::<StereoConfigurationUpdate>()?;
        module.add_class::<Orientation>()?;
        module.add_class::<LigandPermutation>()?;
        module.add_class::<OrientedLigandPermutation>()?;
        module.add_class::<StereoLigandPair>()?;
        module.add_class::<TopicityRelationForm>()?;
        module.add_class::<StereogenicityForm>()?;
        module.add_class::<LigandSymmetryForm>()?;
        module.add_class::<FluxionalityForm>()?;
        module.add_class::<TopicityForm>()?;
        module.add_class::<StereoAtomConstraintKey>()?;
        module.add_class::<StereoAtomConstraintForm>()?;
        module.add_class::<StereoAtomConstraintsForm>()?;
        module.add_class::<StereoBondConstraintKey>()?;
        module.add_class::<StereoBondConstraintForm>()?;
        module.add_class::<StereoBondConstraintsForm>()?;
        module.add_class::<StereoAtomForm>()?;
        module.add_class::<StereoAtomUpdate>()?;
        module.add_class::<StereoBondForm>()?;
        module.add_class::<StereoBondUpdate>()?;
        module.add_class::<StereoAtomConstraintsView>()?;
        module.add_class::<StereoBondConstraintsView>()?;
        module.add_class::<StereoAtomView>()?;
        module.add_class::<StereoBondView>()?;
        module.add_class::<StereoAtomViews>()?;
        module.add_class::<StereoBondViews>()?;
    }
    Ok(())
}
