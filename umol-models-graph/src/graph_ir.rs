//! Graph-based molecular intermediate representation.

pub mod atom;
pub mod atom_matcher;
pub mod atom_spec;
pub mod atom_spec_registry;
pub mod atom_validator;
pub mod bond;
pub mod bond_matcher;
pub mod bond_spec;
pub mod bond_spec_registry;
pub mod convert;
pub mod diagnostics;
pub mod error;
pub mod molecule;

pub use atom::{Atom, AtomBuilder};
pub use atom_matcher::{
    AtomMatcher, ALWAYS_ATOM_MATCHER, DEFAULT_ATOM_MATCHER, LENIENT_ATOM_MATCHER,
    STRICT_ATOM_MATCHER,
};
pub use atom_spec::AtomSpec;
pub use atom_spec_registry::AtomSpecRegistry;
pub use atom_validator::{
    AtomValidator, ALWAYS_ATOM_VALIDATOR, DEFAULT_ATOM_VALIDATOR, LENIENT_ATOM_VALIDATOR,
    STRICT_ATOM_VALIDATOR,
};
pub use bond::{Bond, BondBuilder};
pub use bond_matcher::{
    BondMatcher, ALWAYS_BOND_MATCHER, DEFAULT_BOND_MATCHER, LENIENT_BOND_MATCHER,
    STRICT_BOND_MATCHER,
};
pub use bond_spec::{BondDonation, BondOrder, BondSpec};
pub use bond_spec_registry::BondSpecRegistry;
pub use convert::sir_to_gir;
pub use error::GraphError;
pub use molecule::{AtomIndex, BondIndex, Molecule, MoleculeBuilder};
