//! Geometric primitives for umol.

pub(crate) mod orientation;
pub(crate) mod plane;
pub(crate) mod point;

pub use orientation::{same_side_of_axis, signed_volume, AXIS_SIDE_TOLERANCE};
pub use plane::complementary_direction;
pub use point::{Point2D, Point3D};
