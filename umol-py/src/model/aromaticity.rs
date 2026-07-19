//! Python bindings for aromaticity-model values.

use pyo3::prelude::*;
use umol_graph::ops::model::RingLimits as GraphRingLimits;

/// Ring-size and fused-ring search bounds for aromaticity perception.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingLimits(GraphRingLimits);

#[pymethods]
impl RingLimits {
    #[new]
    #[pyo3(signature = (*, min_ring_size=3, max_ring_size=22, include_fused=true, max_fused_combination=6, max_fused_search=10_000))]
    fn new(
        min_ring_size: usize,
        max_ring_size: usize,
        include_fused: bool,
        max_fused_combination: usize,
        max_fused_search: usize,
    ) -> Self {
        Self(GraphRingLimits {
            min_ring_size,
            max_ring_size,
            include_fused,
            max_fused_combination,
            max_fused_search,
        })
    }

    #[getter]
    fn min_ring_size(&self) -> usize {
        self.0.min_ring_size
    }

    #[getter]
    fn max_ring_size(&self) -> usize {
        self.0.max_ring_size
    }

    #[getter]
    fn include_fused(&self) -> bool {
        self.0.include_fused
    }

    #[getter]
    fn max_fused_combination(&self) -> usize {
        self.0.max_fused_combination
    }

    #[getter]
    fn max_fused_search(&self) -> usize {
        self.0.max_fused_search
    }

    fn __repr__(&self) -> String {
        format!(
            "RingLimits(min_ring_size={}, max_ring_size={}, include_fused={}, max_fused_combination={}, max_fused_search={})",
            self.0.min_ring_size,
            self.0.max_ring_size,
            if self.0.include_fused { "True" } else { "False" },
            self.0.max_fused_combination,
            self.0.max_fused_search,
        )
    }
}

impl RingLimits {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for AromaticityModel configuration"
    )]
    pub(crate) fn from_rust(limits: &GraphRingLimits) -> Self {
        Self(limits.clone())
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for AromaticityModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphRingLimits {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::default(3, 22, true, 6, 10_000, GraphRingLimits::default())]
    #[case::zero(0, 0, false, 0, 0, GraphRingLimits {
        min_ring_size: 0,
        max_ring_size: 0,
        include_fused: false,
        max_fused_combination: 0,
        max_fused_search: 0,
    })]
    #[case::nondefault(5, 18, false, 4, 2_500, GraphRingLimits {
        min_ring_size: 5,
        max_ring_size: 18,
        include_fused: false,
        max_fused_combination: 4,
        max_fused_search: 2_500,
    })]
    fn test_ring_limits_new(
        #[case] min_ring_size: usize,
        #[case] max_ring_size: usize,
        #[case] include_fused: bool,
        #[case] max_fused_combination: usize,
        #[case] max_fused_search: usize,
        #[case] expected: GraphRingLimits,
    ) {
        assert_eq!(
            RingLimits::new(
                min_ring_size,
                max_ring_size,
                include_fused,
                max_fused_combination,
                max_fused_search,
            )
            .0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        RingLimits::new(3, 22, true, 6, 10_000),
        "RingLimits(min_ring_size=3, max_ring_size=22, include_fused=True, max_fused_combination=6, max_fused_search=10000)"
    )]
    #[case::nondefault(
        RingLimits::new(5, 18, false, 4, 2_500),
        "RingLimits(min_ring_size=5, max_ring_size=18, include_fused=False, max_fused_combination=4, max_fused_search=2500)"
    )]
    fn test_ring_limits_repr(#[case] limits: RingLimits, #[case] expected: &str) {
        assert_eq!(limits.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphRingLimits::default())]
    #[case::nondefault(GraphRingLimits {
        min_ring_size: 5,
        max_ring_size: 18,
        include_fused: false,
        max_fused_combination: 4,
        max_fused_search: 2_500,
    })]
    fn test_ring_limits_from_rust(#[case] limits: GraphRingLimits) {
        assert_eq!(RingLimits::from_rust(&limits).0, limits);
    }

    #[rstest]
    #[case::default(RingLimits::new(3, 22, true, 6, 10_000))]
    #[case::nondefault(RingLimits::new(5, 18, false, 4, 2_500))]
    fn test_ring_limits_to_rust(#[case] limits: RingLimits) {
        assert_eq!(limits.to_rust(), limits.0);
    }
}
