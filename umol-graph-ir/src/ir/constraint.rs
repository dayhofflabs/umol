//! Graph-IR constraints: per-scope predicates and their containers.
//!
//! Per-scope enums (`AtomConstraintForm`, `BondConstraintForm`, `DativeBondConstraintForm`,
//! `AromaticSystemConstraintForm`, `MulticenterBondConstraintForm`,
//! `NoncovalentBondConstraintForm`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Per-entity constraints live inline on entity forms via the typed
//! containers (`AtomConstraintsForm`, `BondConstraintsForm`, `DativeBondConstraintsForm`,
//! `AromaticSystemConstraintsForm`, `MulticenterBondConstraintsForm`,
//! `NoncovalentBondConstraintsForm`, `StereoAtomConstraintsForm`, and
//! `StereoBondConstraintsForm`). Each exposes a uniform `new`/`len`/`iter`/
//! `set`/`retain`/`clear` surface; `set` keeps keys unique and overwrites an
//! existing value at the same key.
//! The molecule-level `Constraints` (under `molecule`) is a flat
//! `Vec<Constraint>` for molecule-scope predicates and combinator forms.
//!
//! Each entity container's `asserted` and `asserted_complete` readings return
//! normal forms with meet's `Option` semantics: `None` is bottom and `Some`
//! includes top. Absence contributes top; completion closes top only where the
//! key defines a negative. These readings require no molecule or ring context.
//!
//! # Semantic properties
//!
//! Own-assertion readings interpret only the requested key and leave storage
//! unchanged. Successful container normalization preserves both readings.
//! A molecule-backed view gives the same own-assertion readings as its container.

mod aromatic;
mod atom;
mod bond;
mod dative;
mod molecule;
mod multicenter;
mod noncovalent;
mod relational;
mod ring;
mod stereo;

pub use aromatic::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemConstraintsForm,
};
pub(crate) use atom::read_atom_assertion;
pub use atom::{
    aromatic_covalence, AromaticValence, AromaticValenceForm, AtomConstraintForm,
    AtomConstraintKey, AtomConstraintsForm, MulticenterValence, MulticenterValenceForm,
};
pub use bond::{BondConstraintForm, BondConstraintKey, BondConstraintsForm};
pub use dative::{DativeBondConstraintForm, DativeBondConstraintKey, DativeBondConstraintsForm};
pub use molecule::{Constraint, Constraints, MoleculeConstraint};
pub(crate) use molecule::{
    ConstraintFrameActionDomain, ConstraintFrameActionMap, ConstraintFrameActions,
};
pub use multicenter::{
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondConstraintsForm,
};
pub use noncovalent::{
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, NoncovalentBondConstraintsForm,
};
pub use relational::RelationalConstraint;
pub use ring::{RingMembershipForm, RingScope};
pub use stereo::{
    FluxionalityForm, LigandPermutation, LigandSymmetryForm, OrientedLigandPermutation,
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomConstraintsForm,
    StereoBondConstraintForm, StereoBondConstraintKey, StereoBondConstraintsForm, StereoLigandPair,
    StereogenicityForm, TopicityForm, TopicityRelationForm,
};
