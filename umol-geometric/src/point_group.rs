//! Point group data for Born-Oppenheimer model

/// Marker trait for point group symmetry types:w
pub trait PointGroup {}

/// Trivial point group (no symmetry).
#[derive(Debug, Clone, Copy, Default)]
pub struct C1;
impl PointGroup for C1 {}

