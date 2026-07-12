//! Safe data boundary for umol's vendored nauty integration.
//!
//! Only the stable umol-owned C shim is declared here. Upstream nauty structs,
//! options, statistics, and allocation macros remain private to the C side.

use std::error::Error;
use std::ffi::c_void;
use std::fmt::{self, Display, Formatter};

/// Vendored upstream nauty version.
pub const NAUTY_VERSION: &str = "2.9.3";

type GeneratorCallback = unsafe extern "C" fn(*mut c_void, *const u32, u32);

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(not(test), allow(dead_code))]
enum RawError {
    Ok = 0,
    NullPointer = 1,
    InvalidVertexCount = 2,
    InvalidCsr = 3,
    IntegerOverflow = 4,
    AllocationFailed = 5,
    ReentrantCall = 6,
    NTooBig = 7,
    MTooBig = 8,
    CanonGraphMissing = 9,
    Aborted = 10,
    Killed = 11,
    Unknown = 12,
}

unsafe extern "C" {
    #[allow(dead_code)]
    fn umol_nauty_run(
        vertex_count: u32,
        offsets: *const usize,
        neighbors: *const u32,
        colors: *const u32,
        canonical_labels: *mut u32,
        orbits: *mut u32,
        report_generator: Option<GeneratorCallback>,
        callback_context: *mut c_void,
        group_mantissa: *mut f64,
        group_exponent: *mut i32,
    ) -> RawError;
}

/// Failure while validating input or running the native shim.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NautyError {
    OffsetCount { expected: usize, actual: usize },
    FirstOffset { actual: usize },
    NonmonotonicOffsets { vertex: usize },
    TerminalOffset { expected: usize, actual: usize },
    ColorCount { expected: usize, actual: usize },
    NeighborOutOfBounds { position: usize, neighbor: u32 },
    VertexCountOverflow { count: usize },
    DegreeOverflow { vertex: usize, degree: usize },
    NullPointer,
    InvalidVertexCount,
    InvalidCsr,
    IntegerOverflow,
    AllocationFailed,
    ReentrantCall,
    NTooBig,
    MTooBig,
    CanonGraphMissing,
    Aborted,
    Killed,
    Unknown,
}

impl Display for NautyError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl Error for NautyError {}

#[cfg_attr(not(test), allow(dead_code))]
impl RawError {
    fn into_result(self) -> Result<(), NautyError> {
        match self {
            Self::Ok => Ok(()),
            Self::NullPointer => Err(NautyError::NullPointer),
            Self::InvalidVertexCount => Err(NautyError::InvalidVertexCount),
            Self::InvalidCsr => Err(NautyError::InvalidCsr),
            Self::IntegerOverflow => Err(NautyError::IntegerOverflow),
            Self::AllocationFailed => Err(NautyError::AllocationFailed),
            Self::ReentrantCall => Err(NautyError::ReentrantCall),
            Self::NTooBig => Err(NautyError::NTooBig),
            Self::MTooBig => Err(NautyError::MTooBig),
            Self::CanonGraphMissing => Err(NautyError::CanonGraphMissing),
            Self::Aborted => Err(NautyError::Aborted),
            Self::Killed => Err(NautyError::Killed),
            Self::Unknown => Err(NautyError::Unknown),
        }
    }
}

/// Owned CSR topology and ranked vertex colors accepted by nauty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NautyInput {
    offsets: Vec<usize>,
    neighbors: Vec<u32>,
    colors: Vec<u32>,
}

impl NautyInput {
    pub fn try_new(
        vertex_count: usize,
        offsets: Vec<usize>,
        neighbors: Vec<u32>,
        colors: Vec<u32>,
    ) -> Result<Self, NautyError> {
        validate_vertex_count(vertex_count)?;

        let expected_offsets =
            vertex_count
                .checked_add(1)
                .ok_or(NautyError::VertexCountOverflow {
                    count: vertex_count,
                })?;
        if offsets.len() != expected_offsets {
            return Err(NautyError::OffsetCount {
                expected: expected_offsets,
                actual: offsets.len(),
            });
        }
        if colors.len() != vertex_count {
            return Err(NautyError::ColorCount {
                expected: vertex_count,
                actual: colors.len(),
            });
        }
        if offsets[0] != 0 {
            return Err(NautyError::FirstOffset { actual: offsets[0] });
        }

        for vertex in 0..vertex_count {
            let begin = offsets[vertex];
            let end = offsets[vertex + 1];
            if end < begin {
                return Err(NautyError::NonmonotonicOffsets { vertex });
            }
            validate_degree(vertex, end - begin)?;
        }

        let terminal = offsets[vertex_count];
        if terminal != neighbors.len() {
            return Err(NautyError::TerminalOffset {
                expected: neighbors.len(),
                actual: terminal,
            });
        }
        for (position, &neighbor) in neighbors.iter().enumerate() {
            if neighbor as usize >= vertex_count {
                return Err(NautyError::NeighborOutOfBounds { position, neighbor });
            }
        }

        Ok(Self {
            offsets,
            neighbors,
            colors,
        })
    }

    pub fn vertex_count(&self) -> usize {
        self.colors.len()
    }

    pub fn directed_edge_count(&self) -> usize {
        self.neighbors.len()
    }

    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    pub fn neighbors(&self) -> &[u32] {
        &self.neighbors
    }

    pub fn colors(&self) -> &[u32] {
        &self.colors
    }
}

fn validate_vertex_count(count: usize) -> Result<(), NautyError> {
    if count > i32::MAX as usize {
        Err(NautyError::VertexCountOverflow { count })
    } else {
        Ok(())
    }
}

fn validate_degree(vertex: usize, degree: usize) -> Result<(), NautyError> {
    if degree > i32::MAX as usize {
        Err(NautyError::DegreeOverflow { vertex, degree })
    } else {
        Ok(())
    }
}

/// Group-size representation returned by nauty: `mantissa × 10^exponent`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NautyGroupOrder {
    pub mantissa: f64,
    pub exponent: i32,
}

/// Owned output from one nauty invocation.
#[derive(Clone, Debug, PartialEq)]
pub struct NautyOutput {
    pub canonical_labels: Vec<u32>,
    pub orbits: Vec<u32>,
    pub group_order: NautyGroupOrder,
    pub generators: Vec<Vec<u32>>,
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::empty(0, vec![0], vec![], vec![])]
    #[case::isolated(1, vec![0, 0], vec![], vec![7])]
    #[case::edge(2, vec![0, 1, 2], vec![1, 0], vec![3, 3])]
    fn test_nauty_input_try_new(
        #[case] vertex_count: usize,
        #[case] offsets: Vec<usize>,
        #[case] neighbors: Vec<u32>,
        #[case] colors: Vec<u32>,
    ) {
        let expected = NautyInput {
            offsets: offsets.clone(),
            neighbors: neighbors.clone(),
            colors: colors.clone(),
        };
        assert_eq!(
            NautyInput::try_new(vertex_count, offsets, neighbors, colors),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::offset_count(0, vec![], vec![], vec![], NautyError::OffsetCount { expected: 1, actual: 0 })]
    #[case::color_count(1, vec![0, 0], vec![], vec![], NautyError::ColorCount { expected: 1, actual: 0 })]
    #[case::first_offset(0, vec![1], vec![], vec![], NautyError::FirstOffset { actual: 1 })]
    #[case::nonmonotonic(2, vec![0, 2, 1], vec![0], vec![0, 0], NautyError::NonmonotonicOffsets { vertex: 1 })]
    #[case::terminal(1, vec![0, 1], vec![], vec![0], NautyError::TerminalOffset { expected: 0, actual: 1 })]
    #[case::neighbor(1, vec![0, 1], vec![1], vec![0], NautyError::NeighborOutOfBounds { position: 0, neighbor: 1 })]
    fn test_nauty_input_try_new_error(
        #[case] vertex_count: usize,
        #[case] offsets: Vec<usize>,
        #[case] neighbors: Vec<u32>,
        #[case] colors: Vec<u32>,
        #[case] expected: NautyError,
    ) {
        assert_eq!(
            NautyInput::try_new(vertex_count, offsets, neighbors, colors),
            Err(expected)
        );
    }

    #[rstest]
    #[case::maximum(i32::MAX as usize, Ok(()))]
    #[case::overflow(
        i32::MAX as usize + 1,
        Err(NautyError::VertexCountOverflow { count: i32::MAX as usize + 1 })
    )]
    fn test_validate_vertex_count(#[case] count: usize, #[case] expected: Result<(), NautyError>) {
        assert_eq!(validate_vertex_count(count), expected);
    }

    #[rstest]
    #[case::maximum(4, i32::MAX as usize, Ok(()))]
    #[case::overflow(
        4,
        i32::MAX as usize + 1,
        Err(NautyError::DegreeOverflow { vertex: 4, degree: i32::MAX as usize + 1 })
    )]
    fn test_validate_degree(
        #[case] vertex: usize,
        #[case] degree: usize,
        #[case] expected: Result<(), NautyError>,
    ) {
        assert_eq!(validate_degree(vertex, degree), expected);
    }

    #[rstest]
    #[case::ok(RawError::Ok, Ok(()))]
    #[case::null_pointer(RawError::NullPointer, Err(NautyError::NullPointer))]
    #[case::invalid_vertex_count(RawError::InvalidVertexCount, Err(NautyError::InvalidVertexCount))]
    #[case::invalid_csr(RawError::InvalidCsr, Err(NautyError::InvalidCsr))]
    #[case::integer_overflow(RawError::IntegerOverflow, Err(NautyError::IntegerOverflow))]
    #[case::allocation_failed(RawError::AllocationFailed, Err(NautyError::AllocationFailed))]
    #[case::reentrant_call(RawError::ReentrantCall, Err(NautyError::ReentrantCall))]
    #[case::n_too_big(RawError::NTooBig, Err(NautyError::NTooBig))]
    #[case::m_too_big(RawError::MTooBig, Err(NautyError::MTooBig))]
    #[case::canon_graph_missing(RawError::CanonGraphMissing, Err(NautyError::CanonGraphMissing))]
    #[case::aborted(RawError::Aborted, Err(NautyError::Aborted))]
    #[case::killed(RawError::Killed, Err(NautyError::Killed))]
    #[case::unknown(RawError::Unknown, Err(NautyError::Unknown))]
    fn test_raw_error_into_result(
        #[case] input: RawError,
        #[case] expected: Result<(), NautyError>,
    ) {
        assert_eq!(input.into_result(), expected);
    }
}
