// Core types and traits

use std::fmt::Debug;

/// A trait for types that can link atoms together
pub trait AtomLink {
    type Site: AtomSite;
    type SiteRef: Debug + Clone + Copy + PartialEq + Eq + std::hash::Hash;
}

/// A trait for types that represent atom sites
pub trait AtomSite {
    fn element(&self) -> Option<crate::Element>;
}

/// Index type for atoms in a molecule
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomIndex(pub usize);

/// Index type for bonds in a molecule
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BondIndex(pub usize);
