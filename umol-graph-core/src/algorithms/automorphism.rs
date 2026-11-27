//! Graph automorphisms, canonical labeling, and canonical keys for edge-colored
//! simple graphs.
//!
//! Two backends implement the same edge-colored symmetry semantics:
//! - [`AutomorphismAlgorithm::Nauty`]: vendored sparse nauty extended with
//!   edge-weighted dispatch over `sparsegraph.w` (see [McKay and Piperno,
//!   *Practical Graph Isomorphism, II*
//!   (2014)](https://doi.org/10.1016/j.jsc.2013.09.003) and the
//!   [nauty and Traces site](https://pallini.di.uniroma1.it/)).
//! - [`AutomorphismAlgorithm::TrailIr`]: pure-Rust trail-backed
//!   individualization--refinement combining backtracking undo trails and gated
//!   recursive non-uniform component recursion (Bliss; [Junttila and Kaski,
//!   ALENEX 2007](https://doi.org/10.1137/1.9781611972870.13),
//!   [TAPAS 2011](https://doi.org/10.1007/978-3-642-19754-3_16)) with exact
//!   false-twin and monochromatic true-twin cell collapse (sassy; [Anders and
//!   Schweitzer, ALENEX 2021](https://doi.org/10.1137/1.9781611976472.12)).

use std::fmt::{self, Display, Formatter};

use umol_nauty_sys::{run_weighted as run_weighted_nauty, WeightedNautyInput};

use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

mod trail_ir;

use trail_ir::TrailIrWorkspace;

/// Algorithm selecting how graph automorphisms and canonical labelings are
/// computed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AutomorphismAlgorithm {
    /// Vendored sparse nauty extended with native edge-weighted dispatch over
    /// `sparsegraph.w`.
    ///
    /// See [McKay and Piperno, *Practical Graph Isomorphism, II*
    /// (2014)](https://doi.org/10.1016/j.jsc.2013.09.003) and the
    /// [nauty and Traces site](https://pallini.di.uniroma1.it/).
    Nauty,
    /// Pure-Rust trail-backed individualization--refinement solver for
    /// edge-colored simple graphs.
    ///
    /// Combines a backtracking undo-trail partition and gated recursive
    /// non-uniform component recursion (Bliss; [Junttila and Kaski,
    /// ALENEX 2007](https://doi.org/10.1137/1.9781611972870.13),
    /// [TAPAS 2011](https://doi.org/10.1007/978-3-642-19754-3_16)) with exact
    /// false-twin and monochromatic true-twin cell collapse (sassy; [Anders and
    /// Schweitzer, ALENEX 2021](https://doi.org/10.1137/1.9781611976472.12))
    /// and edge-colored Hopcroft refinement.
    TrailIr,
}

/// Size of an automorphism group, exact when representable without loss and
/// otherwise retained in the solver's scientific representation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AutomorphismGroupOrder {
    Exact(u128),
    Scientific { mantissa: f64, exponent: i32 },
}

impl AutomorphismGroupOrder {
    /// Construct from `mantissa × 10^exponent`, promoting to [`Self::Exact`]
    /// only when the floating-point representation proves the integer exactly.
    pub fn from_scientific(mantissa: f64, exponent: i32) -> Self {
        match exact_scientific_value(mantissa, exponent) {
            Some(value) => Self::Exact(value),
            None => Self::Scientific { mantissa, exponent },
        }
    }

    /// The exact group order when it can be recovered without loss.
    pub fn exact_value(self) -> Option<u128> {
        match self {
            Self::Exact(value) => Some(value),
            Self::Scientific { mantissa, exponent } => exact_scientific_value(mantissa, exponent),
        }
    }
}

impl Display for AutomorphismGroupOrder {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Exact(value) => Display::fmt(value, formatter),
            Self::Scientific { mantissa, exponent } => {
                write!(formatter, "{mantissa}e{exponent}")
            }
        }
    }
}

fn exact_scientific_value(mantissa: f64, exponent: i32) -> Option<u128> {
    const MAX_EXACT_F64_INTEGER: f64 = 9_007_199_254_740_992.0;
    const EXACT_INTEGER_ROUNDING_ULPS: f64 = 2.0;

    let exponent = u32::try_from(exponent).ok()?;
    let factor = 10_u128.checked_pow(exponent)?;
    let value = mantissa * factor as f64;
    let rounded = value.round();
    let rounding_tolerance = value.abs().max(1.0) * f64::EPSILON * EXACT_INTEGER_ROUNDING_ULPS;
    if !value.is_finite()
        || !(0.0..=MAX_EXACT_F64_INTEGER).contains(&value)
        || (value - rounded).abs() > rounding_tolerance
    {
        return None;
    }
    let exact = rounded as u128;
    (exact as f64 == rounded).then_some(exact)
}

/// Symmetry group of an edge-colored graph, issued by [`Graph::automorphisms`],
/// [`Graph::automorphisms_with_scratch`], or [`CanonicalLabeling`].
///
/// Every instance preserves these invariants:
/// - `orbits` has length [`Self::node_count`], and each vertex's orbit
///   representative is the minimum [`NodeId`] in its orbit.
/// - [`Self::orbit_count`] equals the number of distinct orbit representatives.
/// - Every generator in [`Self::generators`] has length [`Self::node_count`]
///   and is a vertex bijection preserving adjacency, node colors, and edge
///   colors.
/// - The subgroup generated by [`Self::generators`] has order
///   [`Self::group_order`] and induces the orbit partition returned by
///   [`Self::orbit_of`] and [`Self::same_orbit`].
///
/// Queries taking a [`NodeId`] assume the identifier belongs to the graph that
/// produced this value.
///
/// # Semantic properties
///
/// - Orbits, [`Self::orbit_count`], [`Self::group_order`], and the generated
///   subgroup are invariant across [`AutomorphismAlgorithm`] variants.
/// - The specific generating basis returned by [`Self::generators`] is
///   algorithm-dependent.
#[derive(Debug, Clone, PartialEq)]
pub struct AutomorphismGroup {
    orbits: Vec<NodeId>,
    node_count: usize,
    orbit_count: usize,
    group_order: AutomorphismGroupOrder,
    generators: Vec<Vec<NodeId>>,
}

impl AutomorphismGroup {
    /// Number of nodes in the graph that produced this group.
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Number of distinct vertex orbits.
    pub fn orbit_count(&self) -> usize {
        self.orbit_count
    }

    /// Minimum [`NodeId`] in the orbit containing `node`.
    ///
    /// # Panics
    ///
    /// Panics if `node.index() >= self.node_count()`.
    pub fn orbit_of(&self, node: NodeId) -> NodeId {
        self.orbits[node.index()]
    }

    /// Whether `first` and `second` lie in the same automorphism orbit.
    ///
    /// # Panics
    ///
    /// Panics if `first.index() >= self.node_count()` or
    /// `second.index() >= self.node_count()`.
    pub fn same_orbit(&self, first: NodeId, second: NodeId) -> bool {
        self.orbits[first.index()] == self.orbits[second.index()]
    }

    /// Order of the automorphism group.
    pub fn group_order(&self) -> AutomorphismGroupOrder {
        self.group_order
    }

    /// A generating set of the automorphism group, each generator stored as a
    /// permutation image over `0..node_count`.
    pub fn generators(&self) -> &[Vec<NodeId>] {
        &self.generators
    }
}

/// Canonical vertex ordering paired with the graph's [`AutomorphismGroup`],
/// issued by [`Graph::canonical_labeling`] or
/// [`Graph::canonical_labeling_with_scratch`].
///
/// [`Self::labels`] has length `self.automorphisms().node_count()` and is a
/// bijection over `0..node_count`, where `labels()[rank]` is the source
/// [`NodeId`] assigned canonical position `rank`.
///
/// # Semantic properties
///
/// - For a fixed [`AutomorphismAlgorithm`] and color ordering, [`Self::labels`]
///   induces an identical canonical colored adjacency matrix across all vertex
///   renumberings of isomorphic graphs.
/// - Tie-breaking between equivalent canonical orderings is
///   algorithm-dependent.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalLabeling {
    labels: Vec<NodeId>,
    automorphisms: AutomorphismGroup,
}

impl CanonicalLabeling {
    /// Canonical vertex order: `labels()[rank]` is the source [`NodeId`]
    /// assigned canonical position `rank`.
    pub fn labels(&self) -> &[NodeId] {
        &self.labels
    }

    /// Automorphism group computed alongside this canonical labeling.
    pub fn automorphisms(&self) -> &AutomorphismGroup {
        &self.automorphisms
    }

    /// Consume this labeling and return its [`AutomorphismGroup`].
    pub fn into_automorphisms(self) -> AutomorphismGroup {
        self.automorphisms
    }
}

/// Opaque reusable scratch storage for [`Graph::automorphisms_with_scratch`]
/// and [`Graph::canonical_labeling_with_scratch`].
///
/// Safe to reuse across arbitrary graphs, sizes, colorings, and
/// [`AutomorphismAlgorithm`] variants.
///
/// # Semantic properties
///
/// Prior calls on a scratch instance never affect the output of subsequent
/// calls; `*_with_scratch` methods produce results identical to their
/// no-scratch counterparts.
#[derive(Clone, Debug, Default)]
pub struct AutomorphismScratch {
    node_ranks: Vec<u32>,
    edge_ranks: Vec<u32>,
    node_order: Vec<u32>,
    edge_order: Vec<u32>,
    trail_ir: TrailIrWorkspace,
}

impl AutomorphismScratch {
    /// Allocate empty reusable scratch storage.
    pub fn new() -> Self {
        Self::default()
    }
}

impl Graph {
    /// Compute the automorphism group of `self` under `node_color` and
    /// `edge_color` using `alg`.
    ///
    /// Skips canonical-labeling construction when only symmetry orbits, group
    /// order, and generators are needed.
    ///
    /// # Semantic properties
    ///
    /// Orbits, orbit count, group order, and the generated subgroup depend only
    /// on the graph and the equivalence classes induced by `node_color` and
    /// `edge_color`, and agree across [`AutomorphismAlgorithm`] variants.
    pub fn automorphisms<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> AutomorphismGroup {
        let mut scratch = AutomorphismScratch::new();
        self.automorphisms_with_scratch(node_color, edge_color, alg, &mut scratch)
    }

    /// Compute the automorphism group of `self`, reusing buffers in `scratch`.
    ///
    /// # Semantic properties
    ///
    /// Produces output identical to [`Self::automorphisms`] regardless of any
    /// prior calls performed with `scratch`.
    pub fn automorphisms_with_scratch<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> AutomorphismGroup {
        self.solve_automorphisms(node_color, edge_color, alg, scratch, false)
            .0
    }

    /// Compute a canonical vertex labeling and automorphism group of `self`
    /// under `node_color` and `edge_color` using `alg`.
    ///
    /// # Semantic properties
    ///
    /// For a fixed [`AutomorphismAlgorithm`] and color ordering, the returned
    /// labeling induces an identical canonical colored adjacency matrix across
    /// all vertex renumberings of isomorphic graphs.
    pub fn canonical_labeling<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> CanonicalLabeling {
        let mut scratch = AutomorphismScratch::new();
        self.canonical_labeling_with_scratch(node_color, edge_color, alg, &mut scratch)
    }

    /// Compute a canonical vertex labeling and automorphism group of `self`,
    /// reusing buffers in `scratch`.
    ///
    /// # Semantic properties
    ///
    /// Produces output identical to [`Self::canonical_labeling`] regardless of
    /// any prior calls performed with `scratch`.
    pub fn canonical_labeling_with_scratch<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> CanonicalLabeling {
        let (automorphisms, labels) =
            self.solve_automorphisms(node_color, edge_color, alg, scratch, true);
        CanonicalLabeling {
            labels,
            automorphisms,
        }
    }

    fn solve_automorphisms<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
        get_canon: bool,
    ) -> (AutomorphismGroup, Vec<NodeId>) {
        let node_count = self.node_count();
        if node_count == 0 {
            return (
                AutomorphismGroup {
                    orbits: Vec::new(),
                    node_count: 0,
                    orbit_count: 0,
                    group_order: AutomorphismGroupOrder::Exact(1),
                    generators: Vec::new(),
                },
                Vec::new(),
            );
        }

        let node_colors: Vec<NC> = self.node_ids().map(node_color).collect();
        scratch.node_order.clear();
        scratch.node_order.extend(0..node_count as u32);
        scratch
            .node_order
            .sort_by(|&a, &b| node_colors[a as usize].cmp(&node_colors[b as usize]));
        scratch.node_ranks.resize(node_count, 0);
        let mut node_rank = 0_u32;
        for position in 0..node_count {
            let vertex = scratch.node_order[position] as usize;
            if position > 0
                && node_colors[vertex] != node_colors[scratch.node_order[position - 1] as usize]
            {
                node_rank += 1;
            }
            scratch.node_ranks[vertex] = node_rank;
        }

        let edge_count = self.edge_count();
        scratch.edge_ranks.resize(edge_count, 0);
        if edge_count > 0 {
            let edge_colors: Vec<EC> = self.edge_ids().map(edge_color).collect();
            if edge_colors[1..]
                .iter()
                .all(|color| color == &edge_colors[0])
            {
                scratch.edge_ranks.fill(0);
            } else {
                scratch.edge_order.clear();
                scratch.edge_order.extend(0..edge_count as u32);
                scratch
                    .edge_order
                    .sort_by(|&a, &b| edge_colors[a as usize].cmp(&edge_colors[b as usize]));
                let mut edge_rank = 0_u32;
                for position in 0..edge_count {
                    let edge = scratch.edge_order[position] as usize;
                    if position > 0
                        && edge_colors[edge]
                            != edge_colors[scratch.edge_order[position - 1] as usize]
                    {
                        edge_rank += 1;
                    }
                    scratch.edge_ranks[edge] = edge_rank;
                }
            }
        }

        match alg {
            AutomorphismAlgorithm::Nauty => {
                let mut offsets = Vec::with_capacity(node_count + 1);
                let mut neighbors = Vec::with_capacity(2 * edge_count);
                let mut weights = Vec::with_capacity(2 * edge_count);
                offsets.push(0);
                for node in self.node_ids() {
                    for neighbor in self.neighbors(node) {
                        neighbors.push(neighbor.node.0);
                        weights.push(scratch.edge_ranks[neighbor.edge.index()]);
                    }
                    offsets.push(neighbors.len());
                }

                let input = WeightedNautyInput::try_new(
                    node_count,
                    offsets,
                    neighbors,
                    weights,
                    scratch.node_ranks.clone(),
                    scratch.node_order.clone(),
                )
                .expect("Graph produces valid weighted nauty input");
                let output = run_weighted_nauty(&input, get_canon)
                    .expect("weighted vendored nauty succeeds");
                let orbit_count = output
                    .orbits
                    .iter()
                    .enumerate()
                    .filter(|&(vertex, &rep)| vertex == rep as usize)
                    .count();
                let group = AutomorphismGroup {
                    orbits: output.orbits.into_iter().map(NodeId).collect(),
                    node_count,
                    orbit_count,
                    group_order: AutomorphismGroupOrder::from_scientific(
                        output.group_order.mantissa,
                        output.group_order.exponent,
                    ),
                    generators: output
                        .generators
                        .into_iter()
                        .map(|generator| generator.into_iter().map(NodeId).collect())
                        .collect(),
                };
                let labels = output.canonical_labels.into_iter().map(NodeId).collect();
                (group, labels)
            }
            AutomorphismAlgorithm::TrailIr => {
                scratch
                    .trail_ir
                    .solve(self, &scratch.node_ranks, &scratch.edge_ranks, get_canon)
            }
        }
    }

    /// Numbering-invariant canonical key of `self` as an edge-colored graph.
    ///
    /// Encodes the canonical node colors followed by the sorted canonical
    /// `(u, v, edge_color)` edge list directly from
    /// [`Self::canonical_labeling`]. Topology only — overlays and N-ary
    /// relations belong to a richer incidence built by the caller.
    ///
    /// # Semantic properties
    ///
    /// For a fixed [`AutomorphismAlgorithm`], two graphs yield equal keys if
    /// and only if they are isomorphic under `node_color` and `edge_color`.
    /// Key bytes may differ across [`AutomorphismAlgorithm`] variants when
    /// canonical tie-breaking differs.
    pub fn canonical_key(
        &self,
        node_color: impl Fn(NodeId) -> Vec<u8>,
        edge_color: impl Fn(EdgeId) -> Vec<u8>,
        alg: AutomorphismAlgorithm,
    ) -> Vec<u8> {
        let node_colors: Vec<Vec<u8>> = self.node_ids().map(node_color).collect();
        let edge_colors: Vec<Vec<u8>> = self.edge_ids().map(edge_color).collect();
        let labeling = self.canonical_labeling(
            |node| node_colors[node.index()].as_slice(),
            |edge| edge_colors[edge.index()].as_slice(),
            alg,
        );
        let canonical = labeling.labels();

        let mut position = vec![0u32; self.node_count()];
        for (rank, &node) in canonical.iter().enumerate() {
            position[node.index()] = rank as u32;
        }

        let mut key = Vec::new();
        key.extend_from_slice(&(self.node_count() as u32).to_le_bytes());
        for &node in canonical {
            let color = &node_colors[node.index()];
            key.extend_from_slice(&(color.len() as u32).to_le_bytes());
            key.extend_from_slice(color);
        }

        let mut canonical_edges: Vec<(u32, u32, &[u8])> = self
            .edge_ids()
            .map(|edge| {
                let [u, v] = self.edge_endpoints(edge);
                let (u, v) = (position[u.index()], position[v.index()]);
                (u.min(v), u.max(v), edge_colors[edge.index()].as_slice())
            })
            .collect();
        canonical_edges.sort_unstable();

        key.extend_from_slice(&(canonical_edges.len() as u32).to_le_bytes());
        for (u, v, color) in canonical_edges {
            key.extend_from_slice(&u.to_le_bytes());
            key.extend_from_slice(&v.to_le_bytes());
            key.extend_from_slice(&(color.len() as u32).to_le_bytes());
            key.extend_from_slice(color);
        }
        key
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::fmt::Debug;
    use std::thread;

    use rstest::*;

    use super::AutomorphismAlgorithm::{Nauty, TrailIr};
    use super::*;
    use crate::union_find::UnionFind;

    fn assert_automorphism_semantics<NC: Debug + Eq, EC: Debug + Eq>(
        graph: &Graph,
        node_colors: &[NC],
        edge_colors: &[EC],
        aut: &AutomorphismGroup,
    ) {
        let expected_nodes: Vec<NodeId> = graph.node_ids().collect();
        let edges: HashMap<(usize, usize), &EC> = graph
            .edge_ids()
            .map(|edge| {
                let [a, b] = graph.edge_endpoints(edge);
                (
                    (a.index().min(b.index()), a.index().max(b.index())),
                    &edge_colors[edge.index()],
                )
            })
            .collect();
        let mut generated_orbits = UnionFind::new(graph.node_count());

        for generator in aut.generators() {
            let mut image = generator.clone();
            image.sort_unstable();
            assert_eq!(image, expected_nodes);
            for (source, &target) in generator.iter().enumerate() {
                assert_eq!(node_colors[source], node_colors[target.index()]);
                generated_orbits.union(source, target.index());
            }
            for (&(a, b), &color) in &edges {
                let a_image = generator[a].index();
                let b_image = generator[b].index();
                assert_eq!(
                    edges.get(&(a_image.min(b_image), a_image.max(b_image))),
                    Some(&color)
                );
            }
        }

        let mut min_rep = vec![usize::MAX; graph.node_count()];
        for v in 0..graph.node_count() {
            let root = generated_orbits.find(v);
            min_rep[root] = min_rep[root].min(v);
        }
        for a in 0..graph.node_count() {
            assert_eq!(
                aut.orbit_of(NodeId(a as u32)).index(),
                min_rep[generated_orbits.find(a)]
            );
            for b in 0..graph.node_count() {
                assert_eq!(
                    generated_orbits.find(a) == generated_orbits.find(b),
                    aut.same_orbit(NodeId(a as u32), NodeId(b as u32))
                );
            }
        }
    }

    fn canonical_colored_adjacency<NC: Clone, EC: Clone + Ord>(
        graph: &Graph,
        node_colors: &[NC],
        edge_colors: &[EC],
        labels: &[NodeId],
    ) -> (Vec<NC>, Vec<(u32, u32, EC)>) {
        let mut position = vec![0u32; graph.node_count()];
        for (rank, &node) in labels.iter().enumerate() {
            position[node.index()] = rank as u32;
        }
        let ordered_nodes = labels
            .iter()
            .map(|&node| node_colors[node.index()].clone())
            .collect();
        let mut ordered_edges: Vec<(u32, u32, EC)> = graph
            .edge_ids()
            .map(|edge| {
                let [u, v] = graph.edge_endpoints(edge);
                let (pu, pv) = (position[u.index()], position[v.index()]);
                (pu.min(pv), pu.max(pv), edge_colors[edge.index()].clone())
            })
            .collect();
        ordered_edges.sort();
        (ordered_nodes, ordered_edges)
    }

    #[rstest]
    #[case::integer(8.0, 0, AutomorphismGroupOrder::Exact(8))]
    #[case::positive_exponent(130.7674368, 10, AutomorphismGroupOrder::Exact(1_307_674_368_000))]
    #[case::negative_exponent(
        1.5,
        -1,
        AutomorphismGroupOrder::Scientific { mantissa: 1.5, exponent: -1 }
    )]
    #[case::fractional(
        1.5,
        0,
        AutomorphismGroupOrder::Scientific { mantissa: 1.5, exponent: 0 }
    )]
    #[case::overflow(
        1.0,
        39,
        AutomorphismGroupOrder::Scientific { mantissa: 1.0, exponent: 39 }
    )]
    #[case::non_finite(
        f64::INFINITY,
        0,
        AutomorphismGroupOrder::Scientific { mantissa: f64::INFINITY, exponent: 0 }
    )]
    fn test_automorphism_group_order_from_scientific(
        #[case] mantissa: f64,
        #[case] exponent: i32,
        #[case] expected: AutomorphismGroupOrder,
    ) {
        assert_eq!(
            AutomorphismGroupOrder::from_scientific(mantissa, exponent),
            expected
        );
    }

    #[rstest]
    #[case::exact(AutomorphismGroupOrder::Exact(12), Some(12))]
    #[case::recoverable(
        AutomorphismGroupOrder::Scientific { mantissa: 2.5, exponent: 2 },
        Some(250)
    )]
    #[case::fractional(
        AutomorphismGroupOrder::Scientific { mantissa: 2.5, exponent: 0 },
        None
    )]
    fn test_automorphism_group_order_exact_value(
        #[case] order: AutomorphismGroupOrder,
        #[case] expected: Option<u128>,
    ) {
        assert_eq!(order.exact_value(), expected);
    }

    #[rstest]
    #[case::exact(AutomorphismGroupOrder::Exact(12), "12")]
    #[case::scientific(
        AutomorphismGroupOrder::Scientific { mantissa: 1.25, exponent: 20 },
        "1.25e20"
    )]
    fn test_automorphism_group_order_display(
        #[case] order: AutomorphismGroupOrder,
        #[case] expected: &str,
    ) {
        assert_eq!(order.to_string(), expected);
    }

    #[rstest]
    #[case::two_orbits(AutomorphismGroup {
        orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
        node_count: 3,
        orbit_count: 2,
        group_order: AutomorphismGroupOrder::Exact(2),
        generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
    })]
    fn test_automorphism_group_queries(#[case] group: AutomorphismGroup) {
        assert_eq!(group.node_count(), 3);
        assert_eq!(group.orbit_count(), 2);
        assert_eq!(group.orbit_of(NodeId(1)), NodeId(0));
        assert!(group.same_orbit(NodeId(0), NodeId(1)));
        assert!(!group.same_orbit(NodeId(0), NodeId(2)));
        assert_eq!(group.group_order(), AutomorphismGroupOrder::Exact(2));
        assert_eq!(group.generators(), &[vec![NodeId(1), NodeId(0), NodeId(2)]]);
    }

    #[rstest]
    #[should_panic]
    #[case::orbit_of_oob(NodeId(2), None)]
    #[should_panic]
    #[case::same_orbit_first_oob(NodeId(2), Some(NodeId(0)))]
    #[should_panic]
    #[case::same_orbit_second_oob(NodeId(0), Some(NodeId(2)))]
    fn test_automorphism_group_out_of_bounds_panic(
        #[case] first: NodeId,
        #[case] second: Option<NodeId>,
    ) {
        let group = Graph::new(2, &[[0, 1]]).automorphisms(|_| 0u8, |_| 0u8, TrailIr);
        match second {
            None => {
                let _ = group.orbit_of(first);
            }
            Some(second) => {
                let _ = group.same_orbit(first, second);
            }
        }
    }

    #[rstest]
    #[case::two_orbits(CanonicalLabeling {
        labels: vec![NodeId(2), NodeId(0), NodeId(1)],
        automorphisms: AutomorphismGroup {
            orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
            node_count: 3,
            orbit_count: 2,
            group_order: AutomorphismGroupOrder::Exact(2),
            generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
        },
    })]
    fn test_canonical_labeling_queries(#[case] labeling: CanonicalLabeling) {
        assert_eq!(labeling.labels(), &[NodeId(2), NodeId(0), NodeId(1)]);
        let expected_group = labeling.automorphisms().clone();
        assert_eq!(expected_group.orbit_count(), 2);
        assert_eq!(labeling.into_automorphisms(), expected_group);
    }

    #[rstest]
    #[case::empty(Graph::default(), vec![], vec![], vec![], AutomorphismGroupOrder::Exact(1))]
    #[case::singleton(
        Graph::new(1, &[]),
        vec![0],
        vec![],
        vec![0],
        AutomorphismGroupOrder::Exact(1)
    )]
    #[case::same_color_edge(
        Graph::new(2, &[[0, 1]]),
        vec![0, 0],
        vec![0],
        vec![0, 0],
        AutomorphismGroupOrder::Exact(2)
    )]
    #[case::different_color_edge(
        Graph::new(2, &[[0, 1]]),
        vec![0, 1],
        vec![0],
        vec![0, 1],
        AutomorphismGroupOrder::Exact(1)
    )]
    #[case::uniform_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0, 0, 0, 0],
        vec![0, 0, 0, 0],
        vec![0, 0, 0, 0],
        AutomorphismGroupOrder::Exact(8)
    )]
    #[case::alternating_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0, 0, 0, 0],
        vec![1, 0, 1, 0],
        vec![0, 0, 0, 0],
        AutomorphismGroupOrder::Exact(4)
    )]
    #[case::one_double_bond_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0, 0, 0, 0],
        vec![1, 0, 0, 0],
        vec![0, 0, 2, 2],
        AutomorphismGroupOrder::Exact(2)
    )]
    #[case::colored_path(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0, 1, 0],
        vec![0, 0],
        vec![0, 1, 0],
        AutomorphismGroupOrder::Exact(2)
    )]
    #[case::three_color_path_4(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        vec![0, 0, 0, 0],
        vec![0, 1, 2],
        vec![0, 1, 2, 3],
        AutomorphismGroupOrder::Exact(1)
    )]
    #[case::disconnected_edges(
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0, 0, 0, 0],
        vec![0, 0],
        vec![0, 0, 0, 0],
        AutomorphismGroupOrder::Exact(8)
    )]
    fn test_graph_automorphisms(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] graph: Graph,
        #[case] node_colors: Vec<u8>,
        #[case] edge_colors: Vec<u8>,
        #[case] expected_orbits: Vec<u32>,
        #[case] expected_order: AutomorphismGroupOrder,
    ) {
        let aut = graph.automorphisms(
            |node| node_colors[node.index()],
            |edge| edge_colors[edge.index()],
            alg,
        );
        let expected_orbit_count = expected_orbits
            .iter()
            .enumerate()
            .filter(|&(vertex, &rep)| vertex == rep as usize)
            .count();
        assert_eq!(aut.node_count(), graph.node_count());
        assert_eq!(aut.orbit_count(), expected_orbit_count);
        assert_eq!(aut.group_order(), expected_order);
        for (a, &rep_a) in expected_orbits.iter().enumerate() {
            assert_eq!(aut.orbit_of(NodeId(a as u32)), NodeId(rep_a));
            for (b, &rep_b) in expected_orbits.iter().enumerate() {
                assert_eq!(
                    aut.same_orbit(NodeId(a as u32), NodeId(b as u32)),
                    rep_a == rep_b
                );
            }
        }
        assert_automorphism_semantics(&graph, &node_colors, &edge_colors, &aut);
    }

    #[rstest]
    #[case::swap(
        Graph::new(2, &[[0, 1]]),
        vec![0u8, 0],
        vec![0u8],
        vec![vec![NodeId(1), NodeId(0)]]
    )]
    #[case::colored_path(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 0],
        vec![0u8, 0],
        vec![vec![NodeId(2), NodeId(1), NodeId(0)]]
    )]
    #[case::trivial_node_color(Graph::new(2, &[[0, 1]]), vec![0u8, 1], vec![0u8], vec![])]
    #[case::trivial_edge_color(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0u8, 1, 0],
        vec![0u8, 1],
        vec![]
    )]
    fn test_graph_automorphisms_generators(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] graph: Graph,
        #[case] node_colors: Vec<u8>,
        #[case] edge_colors: Vec<u8>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let aut = graph.automorphisms(
            |node| node_colors[node.index()],
            |edge| edge_colors[edge.index()],
            alg,
        );
        assert_eq!(aut.generators(), expected);
    }

    #[rstest]
    #[case::complete_13(13, 6_227_020_800)]
    fn test_graph_automorphisms_large_order(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] nodes: usize,
        #[case] expected: u128,
    ) {
        let mut edges = Vec::new();
        for a in 0..nodes as u32 {
            for b in a + 1..nodes as u32 {
                edges.push([a, b]);
            }
        }
        let graph = Graph::new(nodes, &edges);
        let aut = graph.automorphisms(|_| 0u8, |_| 0u8, alg);
        assert_eq!(aut.group_order(), AutomorphismGroupOrder::Exact(expected));
        assert_automorphism_semantics(
            &graph,
            &vec![0u8; nodes],
            &vec![0u8; graph.edge_count()],
            &aut,
        );
    }

    #[rstest]
    #[case::uniform_cycle_6(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        vec![0u8; 6],
        NodeId(2),
        AutomorphismGroupOrder::Exact(2)
    )]
    #[case::opposite_double_bonds_fix_2(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        vec![1u8, 0, 0, 1, 0, 0],
        NodeId(2),
        AutomorphismGroupOrder::Exact(2)
    )]
    #[case::opposite_double_bonds_fix_0(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        vec![1u8, 0, 0, 1, 0, 0],
        NodeId(0),
        AutomorphismGroupOrder::Exact(1)
    )]
    fn test_graph_automorphisms_stabilizer(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] graph: Graph,
        #[case] edge_colors: Vec<u8>,
        #[case] site: NodeId,
        #[case] expected_order: AutomorphismGroupOrder,
    ) {
        let node_colors: Vec<bool> = graph.node_ids().map(|node| node == site).collect();
        let aut = graph.automorphisms(
            |node| node_colors[node.index()],
            |edge| edge_colors[edge.index()],
            alg,
        );
        assert_eq!(aut.group_order(), expected_order);
        assert!(aut
            .generators()
            .iter()
            .all(|generator| generator[site.index()] == site));
        assert_automorphism_semantics(&graph, &node_colors, &edge_colors, &aut);
    }

    #[rstest]
    #[case::non_copy_strings(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec!["C".to_owned(), "C".to_owned(), "C".to_owned(), "C".to_owned()],
        vec![b"double".to_vec(), b"single".to_vec(), b"double".to_vec(), b"single".to_vec()],
        1,
        AutomorphismGroupOrder::Exact(4)
    )]
    #[case::non_copy_asymmetric_edges(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec!["C".to_owned(), "C".to_owned(), "C".to_owned(), "C".to_owned()],
        vec![b"double".to_vec(), b"single".to_vec(), b"single".to_vec(), b"single".to_vec()],
        2,
        AutomorphismGroupOrder::Exact(2)
    )]
    fn test_graph_automorphisms_non_copy_colors(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] graph: Graph,
        #[case] node_colors: Vec<String>,
        #[case] edge_colors: Vec<Vec<u8>>,
        #[case] expected_orbits: usize,
        #[case] expected_order: AutomorphismGroupOrder,
    ) {
        let aut = graph.automorphisms(
            |node| node_colors[node.index()].clone(),
            |edge| edge_colors[edge.index()].clone(),
            alg,
        );
        assert_eq!(aut.orbit_count(), expected_orbits);
        assert_eq!(aut.group_order(), expected_order);
        assert_automorphism_semantics(&graph, &node_colors, &edge_colors, &aut);
    }

    #[rstest]
    fn test_graph_automorphisms_with_scratch() {
        let sequence: Vec<(Graph, Vec<u8>, Vec<u8>, AutomorphismAlgorithm)> = vec![
            (Graph::default(), vec![], vec![], TrailIr),
            (Graph::new(1, &[]), vec![0], vec![], Nauty),
            (
                Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
                vec![0, 0, 0, 0],
                vec![1, 0, 1, 0],
                TrailIr,
            ),
            (
                Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
                vec![0, 0, 0, 0],
                vec![1, 0, 0, 0],
                Nauty,
            ),
            (
                Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
                vec![0, 0, 1, 0, 0, 0],
                vec![1, 0, 0, 1, 0, 0],
                TrailIr,
            ),
            (Graph::new(2, &[[0, 1]]), vec![0, 0], vec![3], Nauty),
            (Graph::default(), vec![], vec![], Nauty),
        ];

        let mut scratch = AutomorphismScratch::new();
        for (graph, node_colors, edge_colors, alg) in &sequence {
            let expected = graph.automorphisms(
                |node| node_colors[node.index()],
                |edge| edge_colors[edge.index()],
                *alg,
            );
            let actual = graph.automorphisms_with_scratch(
                |node| node_colors[node.index()],
                |edge| edge_colors[edge.index()],
                *alg,
                &mut scratch,
            );
            assert_eq!(actual, expected);
        }
    }

    #[rstest]
    #[case::parallel_cycles(8)]
    fn test_graph_automorphisms_concurrency(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] thread_count: usize,
    ) {
        let handles: Vec<_> = (0..thread_count)
            .map(|worker| {
                thread::spawn(move || {
                    let graph = Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]);
                    let site = NodeId((worker % graph.node_count()) as u32);
                    let node_colors: Vec<bool> =
                        graph.node_ids().map(|node| node == site).collect();
                    let edge_a = ((site.index() + 1) % 6, (site.index() + 2) % 6);
                    let edge_b = ((site.index() + 4) % 6, (site.index() + 5) % 6);
                    let edge_colors: Vec<u8> = graph
                        .edge_ids()
                        .map(|edge| {
                            let [u, v] = graph.edge_endpoints(edge);
                            let pair = (u.index(), v.index());
                            let rev = (v.index(), u.index());
                            u8::from(
                                pair == edge_a || rev == edge_a || pair == edge_b || rev == edge_b,
                            )
                        })
                        .collect();
                    let aut = graph.automorphisms(
                        |node| node_colors[node.index()],
                        |edge| edge_colors[edge.index()],
                        alg,
                    );
                    assert_automorphism_semantics(&graph, &node_colors, &edge_colors, &aut);
                    (site, aut)
                })
            })
            .collect();

        for handle in handles {
            let (site, aut) = handle.join().expect("automorphism worker succeeds");
            assert_eq!(aut.group_order(), AutomorphismGroupOrder::Exact(2));
            assert!(aut
                .generators()
                .iter()
                .all(|generator| generator[site.index()] == site));
        }
    }

    #[rstest]
    #[case::empty(Graph::default(), vec![], vec![])]
    #[case::singleton(Graph::new(1, &[]), vec![0u8], vec![])]
    #[case::alternating_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0u8; 4],
        vec![1u8, 0, 1, 0]
    )]
    #[case::one_double_bond_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0u8; 4],
        vec![1u8, 0, 0, 0]
    )]
    #[case::three_color_path_4(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        vec![0u8; 4],
        vec![0u8, 1, 2]
    )]
    fn test_graph_canonical_labeling(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] graph: Graph,
        #[case] node_colors: Vec<u8>,
        #[case] edge_colors: Vec<u8>,
    ) {
        let labeling = graph.canonical_labeling(
            |node| node_colors[node.index()],
            |edge| edge_colors[edge.index()],
            alg,
        );
        let expected_nodes: Vec<NodeId> = graph.node_ids().collect();
        let mut sorted = labeling.labels().to_vec();
        sorted.sort_unstable();
        assert_eq!(sorted, expected_nodes);

        let group = graph.automorphisms(
            |node| node_colors[node.index()],
            |edge| edge_colors[edge.index()],
            alg,
        );
        assert_eq!(labeling.automorphisms().orbit_count(), group.orbit_count());
        assert_eq!(labeling.automorphisms().group_order(), group.group_order());
        for node in graph.node_ids() {
            assert_eq!(
                labeling.automorphisms().orbit_of(node),
                group.orbit_of(node)
            );
        }
        assert_automorphism_semantics(&graph, &node_colors, &edge_colors, labeling.automorphisms());

        let expected_group = labeling.automorphisms().clone();
        assert_eq!(labeling.into_automorphisms(), expected_group);
    }

    #[rstest]
    #[case::one_double_bond_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0u8; 4],
        vec![1u8, 0, 0, 0],
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0u8; 4],
        vec![0u8, 0, 1, 0]
    )]
    #[case::scrambled_3color_path(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        vec![0u8, 1, 1, 0],
        vec![0u8, 1, 2],
        Graph::new(4, &[[2, 0], [0, 3], [3, 1]]),
        vec![1u8, 0, 0, 1],
        vec![0u8, 1, 2]
    )]
    fn test_graph_canonical_labeling_relabeling(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] a: Graph,
        #[case] a_nodes: Vec<u8>,
        #[case] a_edges: Vec<u8>,
        #[case] b: Graph,
        #[case] b_nodes: Vec<u8>,
        #[case] b_edges: Vec<u8>,
    ) {
        let la = a.canonical_labeling(|n| a_nodes[n.index()], |e| a_edges[e.index()], alg);
        let lb = b.canonical_labeling(|n| b_nodes[n.index()], |e| b_edges[e.index()], alg);
        assert_eq!(
            canonical_colored_adjacency(&a, &a_nodes, &a_edges, la.labels()),
            canonical_colored_adjacency(&b, &b_nodes, &b_edges, lb.labels())
        );
    }

    #[rstest]
    fn test_graph_canonical_labeling_with_scratch() {
        let sequence: Vec<(Graph, Vec<u8>, Vec<u8>, AutomorphismAlgorithm)> = vec![
            (Graph::default(), vec![], vec![], Nauty),
            (
                Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
                vec![0, 0, 0, 0],
                vec![1, 0, 1, 0],
                TrailIr,
            ),
            (
                Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
                vec![0, 0, 0, 0],
                vec![1, 0, 0, 0],
                Nauty,
            ),
            (
                Graph::new(3, &[[0, 1], [1, 2]]),
                vec![0, 1, 0],
                vec![2, 0],
                TrailIr,
            ),
            (Graph::new(1, &[]), vec![5], vec![], TrailIr),
        ];

        let mut scratch = AutomorphismScratch::default();
        for (graph, node_colors, edge_colors, alg) in &sequence {
            let expected = graph.canonical_labeling(
                |node| node_colors[node.index()],
                |edge| edge_colors[edge.index()],
                *alg,
            );
            let actual = graph.canonical_labeling_with_scratch(
                |node| node_colors[node.index()],
                |edge| edge_colors[edge.index()],
                *alg,
                &mut scratch,
            );
            assert_eq!(actual, expected);
        }
    }

    #[rstest]
    #[case::empty(Graph::default(), vec![], vec![], Graph::default(), vec![], vec![])]
    #[case::singleton(
        Graph::new(1, &[]),
        vec![vec![3u8]],
        vec![],
        Graph::new(1, &[]),
        vec![vec![3u8]],
        vec![]
    )]
    #[case::reversed_path(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![0u8], vec![1], vec![2]],
        vec![vec![9u8], vec![9]],
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![2u8], vec![1], vec![0]],
        vec![vec![9u8], vec![9]]
    )]
    #[case::swapped_edge(
        Graph::new(2, &[[0, 1]]),
        vec![vec![0u8], vec![1]],
        vec![vec![7u8]],
        Graph::new(2, &[[0, 1]]),
        vec![vec![1u8], vec![0]],
        vec![vec![7u8]]
    )]
    #[case::rotated_colored_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![vec![0u8]; 4],
        vec![vec![1u8], vec![2], vec![3], vec![4]],
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![vec![0u8]; 4],
        vec![vec![3u8], vec![4], vec![1], vec![2]]
    )]
    fn test_graph_canonical_key_isomorphic(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] a: Graph,
        #[case] a_nodes: Vec<Vec<u8>>,
        #[case] a_edges: Vec<Vec<u8>>,
        #[case] b: Graph,
        #[case] b_nodes: Vec<Vec<u8>>,
        #[case] b_edges: Vec<Vec<u8>>,
    ) {
        let key_a = a.canonical_key(
            |node| a_nodes[node.index()].clone(),
            |edge| a_edges[edge.index()].clone(),
            alg,
        );
        let key_b = b.canonical_key(
            |node| b_nodes[node.index()].clone(),
            |edge| b_edges[edge.index()].clone(),
            alg,
        );
        assert_eq!(key_a, key_b);
    }

    #[rstest]
    #[case::node_color(
        Graph::new(2, &[[0, 1]]),
        vec![vec![0u8], vec![1]],
        vec![vec![7u8]],
        Graph::new(2, &[[0, 1]]),
        vec![vec![0u8], vec![2]],
        vec![vec![7u8]]
    )]
    #[case::edge_color(
        Graph::new(2, &[[0, 1]]),
        vec![vec![0u8], vec![1]],
        vec![vec![7u8]],
        Graph::new(2, &[[0, 1]]),
        vec![vec![0u8], vec![1]],
        vec![vec![8u8]]
    )]
    #[case::topology(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![0u8], vec![0], vec![0]],
        vec![vec![7u8], vec![7]],
        Graph::new(3, &[[0, 1], [1, 2], [0, 2]]),
        vec![vec![0u8], vec![0], vec![0]],
        vec![vec![7u8], vec![7], vec![7]]
    )]
    #[case::colored_square_adjacent_vs_opposite(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![vec![0u8]; 4],
        vec![vec![1u8], vec![1], vec![0], vec![0]],
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![vec![0u8]; 4],
        vec![vec![1u8], vec![0], vec![1], vec![0]]
    )]
    fn test_graph_canonical_key_distinct(
        #[values(Nauty, TrailIr)] alg: AutomorphismAlgorithm,
        #[case] a: Graph,
        #[case] a_nodes: Vec<Vec<u8>>,
        #[case] a_edges: Vec<Vec<u8>>,
        #[case] b: Graph,
        #[case] b_nodes: Vec<Vec<u8>>,
        #[case] b_edges: Vec<Vec<u8>>,
    ) {
        let key_a = a.canonical_key(
            |node| a_nodes[node.index()].clone(),
            |edge| a_edges[edge.index()].clone(),
            alg,
        );
        let key_b = b.canonical_key(
            |node| b_nodes[node.index()].clone(),
            |edge| b_edges[edge.index()].clone(),
            alg,
        );
        assert_ne!(key_a, key_b);
    }
}
