//! Safe data boundary for umol's vendored CoordGen integration.
//!
//! Enabling the `native` feature compiles the vendored C++ source and exposes
//! [`generate_coordinates`]. The default feature set does not invoke a C++ compiler.

/// Vendored upstream CoordGen version.
pub const COORDGEN_VERSION: &str = "3.0.2";

#[cfg(feature = "native")]
mod native {
    use std::error::Error;
    use std::fmt::{self, Display, Formatter};

    const OK: i32 = 0;
    const NULL_POINTER: i32 = 1;
    const ATOM_OUT_OF_BOUNDS: i32 = 2;
    const ALLOCATION_FAILED: i32 = 3;
    const BACKEND_EXCEPTION: i32 = 4;

    unsafe extern "C" {
        fn umol_coordgen_generate(
            atom_count: usize,
            atomic_numbers: *const u16,
            bond_count: usize,
            bonds: *const Bond,
            points: *mut Point,
        ) -> i32;
    }

    /// Bond input for native coordinate generation.
    ///
    /// Atom positions are zero-based indices into the atomic-number slice supplied to
    /// [`generate_coordinates`].
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct Bond {
        pub atom_0: usize,
        pub atom_1: usize,
        pub order: u8,
    }

    /// Two-dimensional coordinate returned by CoordGen.
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Point {
        pub x: f64,
        pub y: f64,
    }

    /// Failure while validating coordinate-generation input or running CoordGen.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub enum CoordgenError {
        BondAtomOutOfBounds {
            bond_index: usize,
            atom_index: usize,
            atom_count: usize,
        },
        AllocationFailed,
        BackendException,
        NonFinitePoint {
            atom_index: usize,
        },
        UnexpectedStatus {
            status: i32,
        },
    }

    impl Display for CoordgenError {
        fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
            match self {
                Self::BondAtomOutOfBounds {
                    bond_index,
                    atom_index,
                    atom_count,
                } => write!(
                    formatter,
                    "bond {bond_index} references atom {atom_index} outside frame of size {atom_count}"
                ),
                Self::AllocationFailed => formatter.write_str("CoordGen allocation failed"),
                Self::BackendException => formatter.write_str("CoordGen raised an exception"),
                Self::NonFinitePoint { atom_index } => {
                    write!(formatter, "CoordGen returned a non-finite point for atom {atom_index}")
                }
                Self::UnexpectedStatus { status } => {
                    write!(formatter, "CoordGen returned unexpected status {status}")
                }
            }
        }
    }

    impl Error for CoordgenError {}

    /// Generates one point per input atom while preserving the input atom order.
    ///
    /// Atomic number zero denotes a generic layout atom. Empty input is valid. The backend's
    /// native bond-length scale is preserved; the graph-IR adapter owns normalization.
    ///
    /// # Errors
    ///
    /// Returns [`CoordgenError::BondAtomOutOfBounds`] for a bond endpoint outside the supplied
    /// atom frame. Native allocation failures, exceptions, and non-finite output are reported by
    /// their corresponding variants.
    pub fn generate_coordinates(
        atomic_numbers: &[u16],
        bonds: &[Bond],
    ) -> Result<Vec<Point>, CoordgenError> {
        for (bond_index, bond) in bonds.iter().enumerate() {
            for atom_index in [bond.atom_0, bond.atom_1] {
                if atom_index >= atomic_numbers.len() {
                    return Err(CoordgenError::BondAtomOutOfBounds {
                        bond_index,
                        atom_index,
                        atom_count: atomic_numbers.len(),
                    });
                }
            }
        }

        let mut points = vec![Point { x: 0.0, y: 0.0 }; atomic_numbers.len()];
        let status = unsafe {
            umol_coordgen_generate(
                atomic_numbers.len(),
                atomic_numbers.as_ptr(),
                bonds.len(),
                bonds.as_ptr(),
                points.as_mut_ptr(),
            )
        };
        match status {
            OK => {}
            NULL_POINTER | ATOM_OUT_OF_BOUNDS => {
                return Err(CoordgenError::UnexpectedStatus { status });
            }
            ALLOCATION_FAILED => return Err(CoordgenError::AllocationFailed),
            BACKEND_EXCEPTION => return Err(CoordgenError::BackendException),
            _ => return Err(CoordgenError::UnexpectedStatus { status }),
        }

        if let Some(atom_index) = points
            .iter()
            .position(|point| !point.x.is_finite() || !point.y.is_finite())
        {
            return Err(CoordgenError::NonFinitePoint { atom_index });
        }
        Ok(points)
    }
}

#[cfg(feature = "native")]
pub use native::{generate_coordinates, Bond, CoordgenError, Point};
