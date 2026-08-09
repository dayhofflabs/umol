//! Semantic graph IR layer.

pub(crate) mod aromatic;
pub(crate) mod atom;
pub(crate) mod bond;
pub(crate) mod boolean;
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
pub(crate) mod value;
pub(crate) mod view;

pub use aromatic::{AromaticSystemAst, AromaticSystemUpdate};
pub use atom::{AtomAst, AtomUpdate, ElementForm, IsotopeMass, IsotopeMassForm};
pub use bond::{BondAst, BondUpdate};
pub use boolean::BooleanForm;
pub use coloring::{ConstitutionColoring, ConstitutionFeatures, MoleculeColoring};
pub use constraint::{
    aromatic_covalence, AromaticSystemConstraintAst, AromaticSystemConstraintKey,
    AromaticSystemConstraintsAst, AromaticValence, AromaticValenceAst, AtomConstraintAst,
    AtomConstraintKey, AtomConstraintsAst, BondConstraintAst, BondConstraintKey,
    BondConstraintsAst, Constraint, Constraints, DativeBondConstraintAst, DativeBondConstraintKey,
    DativeBondConstraintsAst, FluxionalityAst, LigandPermutation, LigandSymmetryAst,
    MoleculeConstraint, MulticenterBondConstraintAst, MulticenterBondConstraintKey,
    MulticenterBondConstraintsAst, MulticenterValence, MulticenterValenceAst,
    NoncovalentBondConstraintAst, NoncovalentBondConstraintKey, NoncovalentBondConstraintsAst,
    OrientedLigandPermutation, RelationalConstraint, RingMembershipAst, RingScope,
    StereoAtomConstraintAst, StereoAtomConstraintKey, StereoAtomConstraintsAst,
    StereoBondConstraintAst, StereoBondConstraintKey, StereoBondConstraintsAst, StereoLigandPair,
    StereogenicityAst, SubPatternAnchor, TopicityAst, TopicityRelationAst,
};
pub use correspondence::MoleculeCorrespondence;
pub use dative::{DativeBondAst, DativeBondUpdate};
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
pub use incidence::{IncidenceGraph, IncidenceNodeSelection};
pub use ligand::{StereoLigand, StereoLigandKind};
pub use matching::BondMatching;
pub use molecule::transact::{Transaction, TransactionError};
pub use molecule::{
    spec, AtomArg, Fragment, MoleculeAst, MoleculeBuilder, MoleculeEditor, MoleculeEntries,
    MoleculeEntriesError, MoleculeSpec, MoleculeSpecTerm, Port, PortArg,
};
pub use multicenter::{MulticenterBondAst, MulticenterBondUpdate};
pub use noncovalent::{
    NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindForm, NoncovalentBondUpdate,
};
pub use operators::{MemOp, RelOp};
pub use reaction::ReactionAst;
pub use reaction_derivation::ReactionDerivation;
pub use reaction_span::{ReactionSpanAst, ReactionSpanEntries, ReactionSpanEntriesError};
pub use remap::{IdCompaction, IdRemapping, UndoCompaction};
pub use ring::{
    RingConfig, RingConnection, RingGraph, RingId, RingModel, RingRelation, RingSet, RingSetKind,
};
pub use spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
pub use stereo::{
    CisTransConfiguration, CisTransStereo, CisTransStereoForm, StereoAtomAst, StereoAtomUpdate,
    StereoBondAst, StereoBondUpdate, StereoConfiguration, StereoConfigurationForm,
    StereoConfigurationUpdate, StereoCoset, StereoKind, StereoTerm, Stereogenicity,
    TetrahedralConfiguration, TetrahedralStereo, TetrahedralStereoForm, Topicity,
};
pub use substructure::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
pub use symmetry::{GraphSymmetry, GraphSymmetryConfig, StereoSymmetry};
pub use traits::{
    AsLit, BiEquiv, Canonical, Canonicalize, EntityPatch, Equiv, FromIr, IntoIr, Lattice,
    TryFromIr, TryIntoIr,
};
pub use validate::{
    ConnectivityContradiction, ConnectivityError, ConnectivityModel, ConnectivityValidator,
    ConstraintContradiction, ConstraintError, ConstraintValidateConfig, ConstraintValidator,
    DpoContradiction, DpoError, DpoValidator, EntityStructureContradiction, EntityStructureError,
    EntityStructureValidator, IncidenceConstraintContradiction, IncidenceConstraintValidator,
    MoleculeConstraintContradiction, MoleculeConstraintValidator, ReactionIntegrityContradiction,
    ReactionIntegrityError, ReactionIntegrityValidator, RelationalConstraintContradiction,
    RelationalConstraintValidator, RingConstraintContradiction, RingConstraintValidator,
};
pub use value::{ArithExpr, NumForm, PredExpr};
pub use view::{
    AromaticSystemView, AromaticSystemViewMut, AromaticSystemViews, AtomAutomorphism, AtomView,
    AtomViewMut, AtomViews, BondView, BondViewMut, BondViews, DativeBondView, DativeBondViewMut,
    DativeBondViews, MulticenterBondView, MulticenterBondViewMut, MulticenterBondViews,
    NeighborView, NoncovalentBondView, NoncovalentBondViewMut, NoncovalentBondViews, RingAtomView,
    RingBondView, RingView, RingViews, StereoAtomView, StereoAtomViewMut, StereoAtomViews,
    StereoBondView, StereoBondViewMut, StereoBondViews, StereoLigandView,
};
