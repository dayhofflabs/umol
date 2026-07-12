//! Safe data boundary for umol's vendored nauty integration.
//!
//! Only the stable umol-owned C shim is declared here. Upstream nauty structs,
//! options, statistics, and allocation macros remain private to the C side.

use std::error::Error;
use std::ffi::c_void;
use std::fmt::{self, Display, Formatter};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

/// Vendored upstream nauty version.
pub const NAUTY_VERSION: &str = "2.9.3";

type GeneratorCallback = unsafe extern "C" fn(*mut c_void, *const u32, u32);

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)] // Variants are constructed by the C return ABI.
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
    InvalidPartition = 13,
}

unsafe extern "C" {
    fn umol_nauty_run(
        vertex_count: u32,
        offsets: *const usize,
        neighbors: *const u32,
        colors: *const u32,
        partition: *const u32,
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
    PartitionCount { expected: usize, actual: usize },
    PartitionVertexOutOfBounds { position: usize, vertex: u32 },
    DuplicatePartitionVertex { vertex: u32 },
    NonmonotonicPartitionColors { position: usize },
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
    InvalidPartition,
    GeneratorCallbackPanicked,
}

impl Display for NautyError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl Error for NautyError {}

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
            Self::InvalidPartition => Err(NautyError::InvalidPartition),
        }
    }
}

/// Owned CSR topology and ranked vertex colors accepted by nauty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NautyInput {
    offsets: Vec<usize>,
    neighbors: Vec<u32>,
    colors: Vec<u32>,
    partition: Vec<u32>,
}

impl NautyInput {
    pub fn try_new(
        vertex_count: usize,
        offsets: Vec<usize>,
        neighbors: Vec<u32>,
        colors: Vec<u32>,
        partition: Vec<u32>,
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
        if partition.len() != vertex_count {
            return Err(NautyError::PartitionCount {
                expected: vertex_count,
                actual: partition.len(),
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
        let mut seen = vec![false; vertex_count];
        let mut previous_color = None;
        for (position, &vertex) in partition.iter().enumerate() {
            let vertex = vertex as usize;
            if vertex >= vertex_count {
                return Err(NautyError::PartitionVertexOutOfBounds {
                    position,
                    vertex: partition[position],
                });
            }
            if seen[vertex] {
                return Err(NautyError::DuplicatePartitionVertex {
                    vertex: partition[position],
                });
            }
            let color = colors[vertex];
            if previous_color.is_some_and(|previous| color < previous) {
                return Err(NautyError::NonmonotonicPartitionColors { position });
            }
            seen[vertex] = true;
            previous_color = Some(color);
        }

        Ok(Self {
            offsets,
            neighbors,
            colors,
            partition,
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

    pub fn partition(&self) -> &[u32] {
        &self.partition
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

#[derive(Default)]
struct GeneratorCollector {
    generators: Vec<Vec<u32>>,
    panicked: bool,
}

unsafe extern "C" fn collect_generator(
    context: *mut c_void,
    permutation: *const u32,
    vertex_count: u32,
) {
    let collector = unsafe { &mut *context.cast::<GeneratorCollector>() };
    let result = catch_unwind(AssertUnwindSafe(|| {
        let permutation = unsafe { slice::from_raw_parts(permutation, vertex_count as usize) };
        collector.generators.push(permutation.to_vec());
    }));
    if result.is_err() {
        collector.panicked = true;
    }
}

/// Run sparse nauty and return owned canonical-label, orbit, group-order, and
/// generator output.
pub fn run(input: &NautyInput) -> Result<NautyOutput, NautyError> {
    let vertex_count = input.vertex_count();
    if vertex_count == 0 {
        return Ok(NautyOutput {
            canonical_labels: Vec::new(),
            orbits: Vec::new(),
            group_order: NautyGroupOrder {
                mantissa: 1.0,
                exponent: 0,
            },
            generators: Vec::new(),
        });
    }

    let vertex_count =
        u32::try_from(vertex_count).map_err(|_| NautyError::VertexCountOverflow {
            count: input.vertex_count(),
        })?;
    let mut canonical_labels = vec![0; vertex_count as usize];
    let mut orbits = vec![0; vertex_count as usize];
    let mut group_mantissa = 0.0;
    let mut group_exponent = 0;
    let mut collector = GeneratorCollector::default();

    let status = unsafe {
        umol_nauty_run(
            vertex_count,
            input.offsets.as_ptr(),
            input.neighbors.as_ptr(),
            input.colors.as_ptr(),
            input.partition.as_ptr(),
            canonical_labels.as_mut_ptr(),
            orbits.as_mut_ptr(),
            Some(collect_generator),
            (&mut collector as *mut GeneratorCollector).cast(),
            &mut group_mantissa,
            &mut group_exponent,
        )
    };
    status.into_result()?;
    if collector.panicked {
        return Err(NautyError::GeneratorCallbackPanicked);
    }

    Ok(NautyOutput {
        canonical_labels,
        orbits,
        group_order: NautyGroupOrder {
            mantissa: group_mantissa,
            exponent: group_exponent,
        },
        generators: collector.generators,
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::iter::once;
    use std::thread;

    use rstest::*;

    use super::*;

    const GROUP_ORDER_MANTISSA_TOLERANCE: f64 = 1.0e-9;

    fn assert_output_semantics(input: &NautyInput, output: &NautyOutput) {
        let vertex_count = input.vertex_count();
        let expected_vertices: Vec<u32> = (0..vertex_count as u32).collect();
        let edges: HashSet<(u32, u32)> = (0..vertex_count)
            .flat_map(|vertex| {
                input.neighbors[input.offsets[vertex]..input.offsets[vertex + 1]]
                    .iter()
                    .map(move |&neighbor| (vertex as u32, neighbor))
            })
            .collect();

        let mut canonical_labels = output.canonical_labels.clone();
        canonical_labels.sort_unstable();
        assert_eq!(canonical_labels, expected_vertices);
        assert_eq!(output.orbits.len(), vertex_count);

        for generator in &output.generators {
            let mut image = generator.clone();
            image.sort_unstable();
            assert_eq!(image, expected_vertices);
            for (source, &target) in generator.iter().enumerate() {
                assert_eq!(input.colors[source], input.colors[target as usize]);
            }
            for &(source, target) in &edges {
                assert!(edges.contains(&(generator[source as usize], generator[target as usize])));
            }
        }
    }

    #[rstest]
    #[case::empty(0, vec![0], vec![], vec![], vec![])]
    #[case::isolated(1, vec![0, 0], vec![], vec![7], vec![0])]
    #[case::edge(2, vec![0, 1, 2], vec![1, 0], vec![3, 3], vec![0, 1])]
    fn test_nauty_input_try_new(
        #[case] vertex_count: usize,
        #[case] offsets: Vec<usize>,
        #[case] neighbors: Vec<u32>,
        #[case] colors: Vec<u32>,
        #[case] partition: Vec<u32>,
    ) {
        let expected = NautyInput {
            offsets: offsets.clone(),
            neighbors: neighbors.clone(),
            colors: colors.clone(),
            partition: partition.clone(),
        };
        assert_eq!(
            NautyInput::try_new(vertex_count, offsets, neighbors, colors, partition),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::offset_count(0, vec![], vec![], vec![], vec![], NautyError::OffsetCount { expected: 1, actual: 0 })]
    #[case::color_count(1, vec![0, 0], vec![], vec![], vec![0], NautyError::ColorCount { expected: 1, actual: 0 })]
    #[case::partition_count(1, vec![0, 0], vec![], vec![0], vec![], NautyError::PartitionCount { expected: 1, actual: 0 })]
    #[case::first_offset(0, vec![1], vec![], vec![], vec![], NautyError::FirstOffset { actual: 1 })]
    #[case::nonmonotonic(2, vec![0, 2, 1], vec![0], vec![0, 0], vec![0, 1], NautyError::NonmonotonicOffsets { vertex: 1 })]
    #[case::terminal(1, vec![0, 1], vec![], vec![0], vec![0], NautyError::TerminalOffset { expected: 0, actual: 1 })]
    #[case::neighbor(1, vec![0, 1], vec![1], vec![0], vec![0], NautyError::NeighborOutOfBounds { position: 0, neighbor: 1 })]
    #[case::partition_vertex(2, vec![0, 0, 0], vec![], vec![0, 0], vec![0, 2], NautyError::PartitionVertexOutOfBounds { position: 1, vertex: 2 })]
    #[case::partition_duplicate(2, vec![0, 0, 0], vec![], vec![0, 0], vec![0, 0], NautyError::DuplicatePartitionVertex { vertex: 0 })]
    #[case::partition_colors(2, vec![0, 0, 0], vec![], vec![1, 0], vec![0, 1], NautyError::NonmonotonicPartitionColors { position: 1 })]
    fn test_nauty_input_try_new_error(
        #[case] vertex_count: usize,
        #[case] offsets: Vec<usize>,
        #[case] neighbors: Vec<u32>,
        #[case] colors: Vec<u32>,
        #[case] partition: Vec<u32>,
        #[case] expected: NautyError,
    ) {
        assert_eq!(
            NautyInput::try_new(vertex_count, offsets, neighbors, colors, partition),
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
    #[case::invalid_partition(RawError::InvalidPartition, Err(NautyError::InvalidPartition))]
    fn test_raw_error_into_result(
        #[case] input: RawError,
        #[case] expected: Result<(), NautyError>,
    ) {
        assert_eq!(input.into_result(), expected);
    }

    #[rstest]
    #[case::empty(
        NautyInput::try_new(0, vec![0], vec![], vec![], vec![]).unwrap(),
        NautyGroupOrder { mantissa: 1.0, exponent: 0 },
        vec![],
        vec![]
    )]
    #[case::singleton(
        NautyInput::try_new(1, vec![0, 0], vec![], vec![0], vec![0]).unwrap(),
        NautyGroupOrder { mantissa: 1.0, exponent: 0 },
        vec![0],
        vec![0]
    )]
    #[case::same_color_edge(
        NautyInput::try_new(2, vec![0, 1, 2], vec![1, 0], vec![0, 0], vec![0, 1]).unwrap(),
        NautyGroupOrder { mantissa: 2.0, exponent: 0 },
        vec![0, 1],
        vec![0, 0]
    )]
    #[case::different_color_edge(
        NautyInput::try_new(2, vec![0, 1, 2], vec![1, 0], vec![0, 1], vec![0, 1]).unwrap(),
        NautyGroupOrder { mantissa: 1.0, exponent: 0 },
        vec![0, 1],
        vec![0, 1]
    )]
    #[case::path(
        NautyInput::try_new(
            3,
            vec![0, 1, 3, 4],
            vec![1, 0, 2, 1],
            vec![0, 1, 0],
            vec![0, 2, 1]
        ).unwrap(),
        NautyGroupOrder { mantissa: 2.0, exponent: 0 },
        vec![0, 2, 1],
        vec![0, 1, 0]
    )]
    #[case::cycle(
        NautyInput::try_new(
            4,
            vec![0, 2, 4, 6, 8],
            vec![1, 3, 0, 2, 1, 3, 0, 2],
            vec![0, 0, 0, 0],
            vec![0, 1, 2, 3]
        ).unwrap(),
        NautyGroupOrder { mantissa: 8.0, exponent: 0 },
        vec![0, 2, 1, 3],
        vec![0, 0, 0, 0]
    )]
    #[case::disconnected_edges(
        NautyInput::try_new(
            4,
            vec![0, 1, 2, 3, 4],
            vec![1, 0, 3, 2],
            vec![0, 0, 0, 0],
            vec![0, 1, 2, 3]
        ).unwrap(),
        NautyGroupOrder { mantissa: 8.0, exponent: 0 },
        vec![0, 2, 3, 1],
        vec![0, 0, 0, 0]
    )]
    fn test_run(
        #[case] input: NautyInput,
        #[case] expected_order: NautyGroupOrder,
        #[case] expected_labels: Vec<u32>,
        #[case] expected_orbits: Vec<u32>,
    ) {
        let output = run(&input).unwrap();
        assert_eq!(output.group_order, expected_order);
        assert_eq!(output.canonical_labels, expected_labels);
        assert_eq!(output.orbits, expected_orbits);
        assert_output_semantics(&input, &output);
    }

    #[rstest]
    #[case::same_color_edge(
        NautyInput::try_new(2, vec![0, 1, 2], vec![1, 0], vec![0, 0], vec![0, 1]).unwrap(),
        vec![vec![1, 0]]
    )]
    #[case::different_color_edge(
        NautyInput::try_new(2, vec![0, 1, 2], vec![1, 0], vec![0, 1], vec![0, 1]).unwrap(),
        vec![]
    )]
    fn test_run_generators(#[case] input: NautyInput, #[case] expected: Vec<Vec<u32>>) {
        assert_eq!(run(&input).unwrap().generators, expected);
    }

    #[rstest]
    #[case::complete_15(15, 130.7674368, 10)]
    fn test_run_group_order(
        #[case] vertex_count: usize,
        #[case] expected_mantissa: f64,
        #[case] expected_exponent: i32,
    ) {
        let mut offsets = Vec::with_capacity(vertex_count + 1);
        let mut neighbors = Vec::with_capacity(vertex_count * (vertex_count - 1));
        offsets.push(0);
        for vertex in 0..vertex_count {
            neighbors.extend(
                (0..vertex_count)
                    .filter(|&neighbor| neighbor != vertex)
                    .map(|neighbor| neighbor as u32),
            );
            offsets.push(neighbors.len());
        }
        let input = NautyInput::try_new(
            vertex_count,
            offsets,
            neighbors,
            vec![0; vertex_count],
            (0..vertex_count as u32).collect(),
        )
        .unwrap();
        let order = run(&input).unwrap().group_order;
        assert!((order.mantissa - expected_mantissa).abs() < GROUP_ORDER_MANTISSA_TOLERANCE);
        assert_eq!(order.exponent, expected_exponent);
    }

    #[rstest]
    #[case::cycle_6(
        NautyInput::try_new(
            6,
            vec![0, 2, 4, 6, 8, 10, 12],
            vec![1, 5, 0, 2, 1, 3, 2, 4, 3, 5, 0, 4],
            vec![0, 0, 1, 0, 0, 0],
            vec![0, 1, 3, 4, 5, 2]
        ).unwrap(),
        2
    )]
    fn test_run_stabilizer(#[case] input: NautyInput, #[case] site: usize) {
        let output = run(&input).unwrap();
        assert_eq!(
            output.group_order,
            NautyGroupOrder {
                mantissa: 2.0,
                exponent: 0
            }
        );
        assert!(output
            .generators
            .iter()
            .all(|generator| generator[site] == site as u32));
        assert_output_semantics(&input, &output);
    }

    #[rstest]
    #[case::parallel_cycles(8)]
    fn test_run_concurrency(#[case] thread_count: usize) {
        let handles: Vec<_> = (0..thread_count)
            .map(|site| {
                thread::spawn(move || {
                    let site = site % 6;
                    let mut colors = vec![0; 6];
                    colors[site] = 1;
                    let input = NautyInput::try_new(
                        6,
                        vec![0, 2, 4, 6, 8, 10, 12],
                        vec![1, 5, 0, 2, 1, 3, 2, 4, 3, 5, 0, 4],
                        colors,
                        (0..6)
                            .filter(|&vertex| vertex != site as u32)
                            .chain(once(site as u32))
                            .collect(),
                    )
                    .unwrap();
                    let output = run(&input).unwrap();
                    assert_output_semantics(&input, &output);
                    (site, output)
                })
            })
            .collect();

        for handle in handles {
            let (site, output) = handle.join().expect("nauty worker succeeds");
            assert_eq!(
                output.group_order,
                NautyGroupOrder {
                    mantissa: 2.0,
                    exponent: 0
                }
            );
            assert!(output
                .generators
                .iter()
                .all(|generator| generator[site] == site as u32));
        }
    }
}
