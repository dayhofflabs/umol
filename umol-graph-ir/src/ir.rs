//! Semantic graph IR layer.

pub(crate) mod aromatic;
pub(crate) mod atom;
pub(crate) mod bond;
pub(crate) mod boolean;
pub(crate) mod canonicalization;
pub(crate) mod coloring;
pub(crate) mod compose;
pub(crate) mod constraint;
pub(crate) mod correspondence;
pub(crate) mod dative;
pub(crate) mod delta;
pub(crate) mod edit;
pub(crate) mod electrons;
pub(crate) mod entity;
pub(crate) mod error;
pub(crate) mod id;
pub(crate) mod incidence;
pub(crate) mod ligand;
pub(crate) mod matching;
pub(crate) mod molecule;
pub(crate) mod multicenter;
pub(crate) mod noncovalent;
pub(crate) mod num;
pub(crate) mod operators;
pub(crate) mod reaction;
pub(crate) mod reaction_derivation;
pub(crate) mod reaction_span;
pub(crate) mod remap;
pub(crate) mod ring;
pub(crate) mod spin;
pub(crate) mod stereo;
pub(crate) mod substructure;
pub(crate) mod symmetry;
pub(crate) mod traits;
pub(crate) mod validate;
pub(crate) mod view;

pub use aromatic::{AromaticSystemForm, AromaticSystemUpdate};
pub use atom::{AtomForm, AtomUpdate, ElementForm, IsotopeMass, IsotopeMassForm};
pub use bond::{BondForm, BondUpdate};
pub use boolean::BooleanForm;
pub use canonicalization::{
    CanonicalizationContext, CanonicalizationLevel, MoleculeCanonicalizationError,
    ReactionCanonicalizationError, ReactionSpanCanonicalizationError,
};
pub use coloring::{ConstitutionColoring, MoleculeColoring, MoleculeColoringFeatures};
pub use constraint::{
    aromatic_covalence, AromaticSystemConstraintForm, AromaticSystemConstraintKey,
    AromaticSystemConstraintsForm, AromaticValence, AromaticValenceForm, AtomConstraintForm,
    AtomConstraintKey, AtomConstraintsForm, BondConstraintForm, BondConstraintKey,
    BondConstraintsForm, Constraint, Constraints, DativeBondConstraintForm,
    DativeBondConstraintKey, DativeBondConstraintsForm, FluxionalityForm, LigandPermutation,
    LigandSymmetryForm, MoleculeConstraint, MulticenterBondConstraintForm,
    MulticenterBondConstraintKey, MulticenterBondConstraintsForm, MulticenterValence,
    MulticenterValenceForm, NoncovalentBondConstraintForm, NoncovalentBondConstraintKey,
    NoncovalentBondConstraintsForm, OrientedLigandPermutation, RelationalConstraint,
    RingMembershipForm, RingScope, StereoAtomConstraintForm, StereoAtomConstraintKey,
    StereoAtomConstraintsForm, StereoBondConstraintForm, StereoBondConstraintKey,
    StereoBondConstraintsForm, StereoLigandPair, StereogenicityForm, TopicityForm,
    TopicityRelationForm,
};
pub use correspondence::MoleculeCorrespondence;
pub use dative::{DativeBondForm, DativeBondUpdate};
pub use delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, ConstraintSpan, DativeBondDelta,
    Delta, Deltas, EntitySpan, MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta,
    StereoBondDelta,
};
pub use edit::{
    AddBond, AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemFieldChange,
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    CascadedConstraints, ConstraintEdit, ConstraintEditError, DativeBondFieldChange,
    DativeBondHandle, Edit, Edits, EntityHandle, ModifiedConstraint, MulticenterBondFieldChange,
    MulticenterBondHandle, NoncovalentBondFieldChange, NoncovalentBondHandle,
    RemovedAromaticSystem, RemovedAtom, RemovedBond, RemovedConstraint, RemovedDativeBond,
    RemovedMulticenterBond, RemovedNoncovalentBond, RemovedOverlays, RemovedStereoAtom,
    RemovedStereoBond, StereoAtomFieldChange, StereoAtomHandle, StereoAtomRemoval,
    StereoBondFieldChange, StereoBondHandle, StereoBondRemoval, Undo,
};
pub use electrons::ElectronCountsForm;
pub use entity::{Entity, EntityKind};
pub use error::{ApplyError, ApplyPreconditionError, Contradiction, NoJoin};
pub use id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId, StereoLigandPosition,
};
pub use incidence::{IncidenceGraph, IncidenceLevel};
pub use ligand::{StereoLigand, StereoLigandKind};
pub use matching::BondMatching;
pub use molecule::transact::{Transaction, TransactionError};
pub use molecule::{
    spec, AtomArg, Fragment, Molecule, MoleculeBuilder, MoleculeEditor, MoleculeEntries,
    MoleculeIntegrityError, MoleculeSpec, MoleculeSpecTerm, Port, PortArg,
};
pub use multicenter::{MulticenterBondForm, MulticenterBondUpdate};
pub use noncovalent::{
    NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm, NoncovalentBondUpdate,
};
pub use num::{ArithExpr, NumForm, PredExpr};
pub use operators::{MemOp, RelOp};
pub use reaction::Reaction;
pub use reaction_derivation::ReactionDerivation;
pub use reaction_span::{ReactionSpan, ReactionSpanEntries, ReactionSpanIntegrityError};
pub use remap::{IdCompaction, IdRemapping, UndoCompaction};
pub use ring::{
    RingConfig, RingConnection, RingGraph, RingId, RingModel, RingRelation, RingSet, RingSetKind,
};
pub use spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
pub use stereo::{
    CisTransConfiguration, CisTransStereo, CisTransStereoForm, StereoAtomForm, StereoAtomUpdate,
    StereoBondForm, StereoBondUpdate, StereoConfiguration, StereoConfigurationForm,
    StereoConfigurationUpdate, StereoCoset, StereoKind, StereoTerm, Stereogenicity,
    TetrahedralConfiguration, TetrahedralStereo, TetrahedralStereoForm, Topicity,
};
pub use substructure::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
pub use symmetry::{GraphSymmetry, GraphSymmetryConfig, StereoSymmetry};
pub use traits::{
    AsLit, BiRelationEquiv, EntityPatch, Equiv, FromIr, IntoIr, Lattice, Normalize, Normalized,
    RelationEquiv, TryFromIr, TryIntoIr,
};
pub use validate::{DpoContradiction, ReactionIntegrityError};
pub use view::{
    AromaticSystemView, AromaticSystemViewMut, AromaticSystemViews, AtomAutomorphism, AtomView,
    AtomViewMut, AtomViews, BondView, BondViewMut, BondViews, DativeBondView, DativeBondViewMut,
    DativeBondViews, MulticenterBondView, MulticenterBondViewMut, MulticenterBondViews,
    NeighborView, NoncovalentBondView, NoncovalentBondViewMut, NoncovalentBondViews, RingAtomView,
    RingBondView, RingView, RingViews, StereoAtomView, StereoAtomViewMut, StereoAtomViews,
    StereoBondView, StereoBondViewMut, StereoBondViews, StereoLigandView,
};
