//! Geometric primitives for umol.

pub(crate) mod orientation;
pub(crate) mod plane;
pub(crate) mod point;

pub use orientation::signed_volume;
pub use plane::complementary_direction;
pub use point::{Point2D, Point3D};
