//! Surface DSL layer.
//!
//! `umol_ast::dsl` is the canonical import location for every DSL type.
//! Sub-modules are `pub(crate)` implementation detail; external users access
//! everything through this facade.

pub(crate) mod aromatic;
pub(crate) mod atom;
pub(crate) mod bond;
pub(crate) mod config;
pub(crate) mod constraint;
pub(crate) mod dative;
pub(crate) mod error;
pub(crate) mod molecule;
pub(crate) mod multicenter;
pub(crate) mod noncovalent;
mod predicates;
pub(crate) mod refs;
pub(crate) mod value;

pub use aromatic::{AromaticPredicate, AromaticSystemConstraintDsl, AromaticSystemDsl, parse_aromatic};
pub use atom::{
    AromaticValenceDsl, AtomConstraintDsl, AtomDsl, AtomPredicate, MulticenterValenceDsl,
    parse_atom,
};
pub use bond::{BondConstraintDsl, BondDsl, BondPredicate, parse_bond};
pub use config::{
    AromaticSystemDefaults, AromaticSystemOverrides, AromaticValenceDefault, AtomDefaults,
    AtomOverrides, BondDefaults, BondOverrides, DativeBondDefaults, DativeBondOverrides,
    ImplicitDefault, IsotopeDefault, MoleculeDefaults, MoleculeOverrides, MulticenterBondDefaults,
    MulticenterBondOverrides, MulticenterValenceDefault, MultiplicityDefault,
    NoncovalentBondDefaults, NoncovalentBondOverrides, NumericDefault, UnpairedElectronsDefault,
};
pub use constraint::{
    AromaticSystemRef, AtomRef, BondRef, ConstraintDsl, ConstraintsDsl, DativeBondRef,
    MoleculeConstraintDsl, MulticenterBondRef, NoncovalentBondRef, SubPatternAnchorDsl,
};
pub use dative::{DativeBondConstraintDsl, DativeBondDsl, DativePredicate, parse_dative};
pub use error::ParseError;
pub use molecule::{Metadata, MoleculeDsl};
pub use multicenter::{
    MulticenterBondConstraintDsl, MulticenterBondDsl, MulticenterPredicate, parse_multicenter,
};
pub use noncovalent::{NoncovalentBondConstraintDsl, NoncovalentBondDsl, parse_noncovalent};
pub use refs::{
    AromaticSystemRefDsl, AtomRefDsl, BondRefDsl, DativeBondRefDsl, MulticenterBondRefDsl,
    NoncovalentBondRefDsl,
};
pub use value::{ValueDsl, parse_value};
