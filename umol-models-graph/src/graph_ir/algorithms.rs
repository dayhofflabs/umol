//! Index-agnostic graph algorithms over contiguous integer node ids.
//!
//! Kernels in this module operate on `usize` node indices (`0..node_count`).
//! Domain layers are responsible for converting typed indices to/from this form.

pub mod bcc;
pub mod cycles;
pub mod mis;
