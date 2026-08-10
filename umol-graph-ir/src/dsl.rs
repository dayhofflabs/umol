//! Surface DSLs for molecules, reactions, and reaction spans.
//!
//! Whole-value DSL wrappers pair a graph-IR value with persistent surface metadata:
//! entity keywords and atom aliases used for lossless DSL rendering. Parsing
//! uses short-lived molecule and reaction contexts that additionally maintain
//! entity counts and participant indexes while resolving positional, keyword,
//! and structural references.

pub(crate) mod aromatic;
pub(crate) mod atom;
pub(crate) mod bond;
pub(crate) mod boolean;
pub(crate) mod config;
pub(crate) mod constraint;
pub(crate) mod dative;
pub(crate) mod edit;
pub(crate) mod edn_utils;
pub(crate) mod electrons;
pub(crate) mod error;
mod metadata;
pub(crate) mod molecule;
pub(crate) mod multicenter;
pub(crate) mod namespace;
pub(crate) mod noncovalent;
pub(crate) mod num;
pub(crate) mod operators;
mod predicate;
pub(crate) mod reaction;
pub(crate) mod reaction_span;
pub(crate) mod refs;
pub(crate) mod relational;
pub(crate) mod stereo;

pub use aromatic::{
    parse_aromatic_system, parse_aromatic_system_update, AromaticSystemConstraintDsl,
    AromaticSystemDsl, AromaticSystemPredicate, AromaticSystemUpdateDsl,
};
pub use atom::{
    parse_atom, parse_atom_update, AromaticValenceDsl, AtomConstraintDsl, AtomDsl, AtomPredicate,
    AtomUpdateDsl, MulticenterValenceDsl,
};
pub use bond::{
    parse_bond, parse_bond_update, BondConstraintDsl, BondDsl, BondPredicate, BondUpdateDsl,
};
pub use boolean::{parse_boolean, BooleanDsl};
pub use config::{
    AromaticSystemDefaults, AromaticSystemOverrides, AromaticValenceDefault, AtomDefaults,
    AtomOverrides, BondDefaults, BondOverrides, DativeBondDefaults, DativeBondOverrides,
    DeltaDefaults, IsotopeDefault, MoleculeDefaults, MoleculeOverrides, MulticenterBondDefaults,
    MulticenterBondOverrides, MulticenterValenceDefault, MultiplicityDefault,
    NoncovalentBondDefaults, NoncovalentBondOverrides, NumericDefault, ReactionDefaults,
    ReactionOverrides, StereoAtomDefaults, StereoAtomOverrides, StereoBondDefaults,
    StereoBondOverrides, StereoDefault, UnpairedElectronsDefault,
};
pub use constraint::{ConstraintDsl, ConstraintsDsl, MoleculeConstraintDsl};
pub use dative::{
    parse_dative_bond, parse_dative_bond_update, DativeBondConstraintDsl, DativeBondDsl,
    DativeBondPredicate, DativeBondUpdateDsl,
};
pub use edit::EditsDsl;
pub use error::ParseError;
pub use metadata::{Metadata, MetadataError, MoleculeMetadata, ReactionMetadata};
pub use molecule::MoleculeDsl;
pub use multicenter::{
    parse_multicenter_bond, parse_multicenter_bond_update, MulticenterBondConstraintDsl,
    MulticenterBondDsl, MulticenterBondPredicate, MulticenterBondUpdateDsl,
};
pub use namespace::{MoleculeContext, Namespace};
pub use noncovalent::{
    parse_noncovalent_bond, parse_noncovalent_bond_update, NoncovalentBondConstraintDsl,
    NoncovalentBondDsl, NoncovalentBondUpdateDsl,
};
pub use num::{parse_num, NumDsl};
pub use reaction::{ReactionContext, ReactionDsl};
pub use reaction_span::ReactionSpanDsl;
pub use refs::{
    AromaticSystemRef, AtomRef, BondRef, DativeBondParticipants, DativeBondRef, MulticenterBondRef,
    NoncovalentBondRef, StereoAtomParticipants, StereoAtomRef, StereoBondParticipants,
    StereoBondRef, StereoLigandRef,
};
pub use relational::RelationalConstraintDsl;
pub use stereo::{
    parse_stereo_atom, parse_stereo_atom_update, parse_stereo_bond, parse_stereo_bond_update,
    StereoAtomConstraintDsl, StereoAtomDsl, StereoAtomUpdateDsl, StereoBondConstraintDsl,
    StereoBondDsl, StereoBondUpdateDsl,
};
