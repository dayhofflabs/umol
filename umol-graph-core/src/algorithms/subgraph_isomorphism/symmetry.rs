//! Symmetry-breaking constraints for subgraph-isomorphism search.
//!
//! Provides explicit symmetry carriers for both the query graph ([`QuerySymmetryConstraints`])
//! and the target graph ([`TargetSymmetry`]), plus crate-internal stabilizer-chain and combined
//! double-orbit canonicality helpers consumed by the matching engines.
//!
//! Both carriers use the deterministic query vertex priority `(Reverse(degree(q)), q)` and
//! compute exact orbits and point stabilizers via Schreier--Sims stabilization, so all orbit
//! partitions, Grochow--Kellis inequalities, and canonical orbit representatives are invariant
//! across [`AutomorphismAlgorithm`] variants and generating bases.

use std::cmp::Reverse;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use crate::algorithms::automorphism::{
    AutomorphismAlgorithm, AutomorphismGroup, AutomorphismScratch,
};
use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

/// Error returned when externally supplied symmetry generators are invalid for a graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymmetryGeneratorError {
    /// A generator's image slice does not have length `node_count`.
    LengthMismatch {
        generator_index: usize,
        len: usize,
        node_count: usize,
    },
    /// A generator maps a node outside `0..node_count`.
    NodeOutOfBounds {
        generator_index: usize,
        node: NodeId,
        node_count: usize,
    },
    /// A generator is not a bijection: two distinct nodes map to `image`.
    DuplicateImage {
        generator_index: usize,
        image: NodeId,
    },
    /// A generator does not map a graph edge to an edge of the graph.
    EdgeNotPreserved {
        generator_index: usize,
        endpoints: [NodeId; 2],
    },
}

impl Display for SymmetryGeneratorError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match *self {
            Self::LengthMismatch {
                generator_index,
                len,
                node_count,
            } => write!(
                f,
                "generator {generator_index} has length {len}, expected {node_count}"
            ),
            Self::NodeOutOfBounds {
                generator_index,
                node,
                node_count,
            } => write!(
                f,
                "generator {generator_index} maps a node to {node:?}, out of bounds for {node_count} nodes"
            ),
            Self::DuplicateImage {
                generator_index,
                image,
            } => write!(
                f,
                "generator {generator_index} is not a bijection: duplicate image {image:?}"
            ),
            Self::EdgeNotPreserved {
                generator_index,
                endpoints: [first, second],
            } => write!(
                f,
                "generator {generator_index} does not preserve edge [{first:?}, {second:?}]"
            ),
        }
    }
}

impl Error for SymmetryGeneratorError {}

/// Deterministic query vertex priority `(Reverse(degree(q)), q)` shared across
/// Grochow--Kellis query constraints, host stabilizer push-down, and double-orbit canonicality.
pub(crate) fn query_priority_order(query: &Graph) -> Vec<u32> {
    let mut order: Vec<u32> = (0..query.node_count() as u32).collect();
    order.sort_unstable_by_key(|&q| (Reverse(query.degree(NodeId(q))), q));
    order
}

/// Compact support bitset (`v` set iff `perm[v] != v`) for exact disjoint-support pruning.
#[derive(Clone, Debug)]
enum Support {
    Inline(u64),
    Wide(Vec<u64>),
}

impl Support {
    fn from_perm(perm: &[u32]) -> Self {
        let n = perm.len();
        if n <= 64 {
            let mut mask = 0u64;
            for (i, &img) in perm.iter().enumerate() {
                if img as usize != i {
                    mask |= 1u64 << i;
                }
            }
            Self::Inline(mask)
        } else {
            let mut words = vec![0u64; n.div_ceil(64)];
            for (i, &img) in perm.iter().enumerate() {
                if img as usize != i {
                    words[i / 64] |= 1u64 << (i % 64);
                }
            }
            Self::Wide(words)
        }
    }

    fn is_disjoint(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Inline(a), Self::Inline(b)) => (a & b) == 0,
            (Self::Wide(a), Self::Wide(b)) => a.iter().zip(b).all(|(&x, &y)| (x & y) == 0),
            _ => unreachable!("supports in the same group have matching width"),
        }
    }
}

/// Schreier tree for one base point in a stabilizer chain.
#[derive(Clone, Debug)]
struct SchreierTree {
    base: u32,
    orbit: Vec<u32>,
    parent_pt: Vec<u32>,
    parent_gen: Vec<u32>,
    generators: Vec<Vec<u32>>,
    inv_generators: Vec<Vec<u32>>,
}

const UNVISITED: u32 = u32::MAX;
const ROOT_GEN: u32 = u32::MAX;

impl SchreierTree {
    /// Build the exact orbit of `base` and its BFS Schreier tree under `generators`.
    fn build(base: usize, node_count: usize, generators: &[Vec<u32>]) -> Self {
        let mut parent_pt = vec![UNVISITED; node_count];
        let mut parent_gen = vec![ROOT_GEN; node_count];
        let mut orbit = Vec::new();
        parent_pt[base] = base as u32;
        orbit.push(base as u32);

        let mut head = 0;
        while head < orbit.len() {
            let x = orbit[head] as usize;
            head += 1;
            for (gen_idx, gen) in generators.iter().enumerate() {
                let y = gen[x] as usize;
                if parent_pt[y] == UNVISITED {
                    parent_pt[y] = x as u32;
                    parent_gen[y] = gen_idx as u32;
                    orbit.push(y as u32);
                }
            }
        }

        let inv_generators = generators.iter().map(|g| invert_perm(g)).collect();
        Self {
            base: base as u32,
            orbit,
            parent_pt,
            parent_gen,
            generators: generators.to_vec(),
            inv_generators,
        }
    }

    /// Reconstruct `u_x` (mapping `base -> x`) in `out`.
    fn coset_rep_into(&self, target: usize, out: &mut [u32], tmp: &mut [u32]) {
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = i as u32;
        }
        let mut cur = target;
        while cur != self.base as usize {
            let g = &self.generators[self.parent_gen[cur] as usize];
            for (i, slot) in tmp.iter_mut().enumerate() {
                *slot = out[g[i] as usize];
            }
            out.copy_from_slice(tmp);
            cur = self.parent_pt[cur] as usize;
        }
    }

    /// Compute `out = g_cur ∘ u_x` where `u_x(base) = x`.
    fn compose_right_coset_rep(
        &self,
        g_cur: &[u32],
        target: usize,
        out: &mut [u32],
        tmp: &mut [u32],
    ) {
        out.copy_from_slice(g_cur);
        let mut cur = target;
        while cur != self.base as usize {
            let g = &self.generators[self.parent_gen[cur] as usize];
            for (i, slot) in tmp.iter_mut().enumerate() {
                *slot = out[g[i] as usize];
            }
            out.copy_from_slice(tmp);
            cur = self.parent_pt[cur] as usize;
        }
    }

    /// Apply `u_x^{-1}` (mapping `x -> base`) on the left to `perm` in place: `perm <- u_x^{-1} ∘ perm`.
    fn apply_left_inv_coset_rep(&self, target: usize, perm: &mut [u32]) {
        let mut cur = target;
        while cur != self.base as usize {
            let g_inv = &self.inv_generators[self.parent_gen[cur] as usize];
            for slot in perm.iter_mut() {
                *slot = g_inv[*slot as usize];
            }
            cur = self.parent_pt[cur] as usize;
        }
    }
}

fn invert_perm(perm: &[u32]) -> Vec<u32> {
    let mut inv = vec![0u32; perm.len()];
    for (i, &img) in perm.iter().enumerate() {
        inv[img as usize] = i as u32;
    }
    inv
}

fn first_moved(perm: &[u32]) -> Option<usize> {
    perm.iter()
        .enumerate()
        .find_map(|(i, &img)| (img as usize != i).then_some(i))
}

/// Sims's filter: maintains a small generating set with identical generated subgroup.
struct SimsFilter {
    /// For each `i`, entries `(j, gen_idx)` where generator `gen_idx` fixes `0..i` and maps `i -> j`.
    buckets: Vec<Vec<(u32, usize)>>,
    generators: Vec<Vec<u32>>,
    inv_generators: Vec<Vec<u32>>,
}

impl SimsFilter {
    fn new(node_count: usize) -> Self {
        Self {
            buckets: vec![Vec::new(); node_count],
            generators: Vec::new(),
            inv_generators: Vec::new(),
        }
    }

    fn insert(&mut self, mut perm: Vec<u32>) {
        while let Some(i) = first_moved(&perm) {
            let j = perm[i];
            if let Some(&(_, gen_idx)) = self.buckets[i].iter().find(|&&(img, _)| img == j) {
                let g_inv = &self.inv_generators[gen_idx];
                for slot in &mut perm {
                    *slot = g_inv[*slot as usize];
                }
            } else {
                let gen_idx = self.generators.len();
                let inv = invert_perm(&perm);
                self.buckets[i].push((j, gen_idx));
                self.generators.push(perm);
                self.inv_generators.push(inv);
                break;
            }
        }
    }

    fn into_generators(self) -> Vec<Vec<u32>> {
        self.generators
    }
}

/// Compute an exact Sims-filtered generating set for `Stab_<generators>(tree.base)`.
fn stabilize_point(tree: &SchreierTree, node_count: usize) -> Vec<Vec<u32>> {
    if tree.generators.is_empty() {
        return Vec::new();
    }
    let base = tree.base as usize;
    if tree.orbit.len() <= 1 {
        return tree.generators.clone();
    }

    let mut filter = SimsFilter::new(node_count);
    let gen_supports: Vec<Support> = tree
        .generators
        .iter()
        .map(|g| Support::from_perm(g))
        .collect();

    // Seed existing generators that already fix `base` (the `x == base` Schreier generators).
    for g in &tree.generators {
        if g[base] as usize == base {
            filter.insert(g.clone());
        }
    }

    let mut u_x = vec![0u32; node_count];
    let mut tmp = vec![0u32; node_count];
    let mut candidate = vec![0u32; node_count];

    for &x_u32 in &tree.orbit {
        let x = x_u32 as usize;
        tree.coset_rep_into(x, &mut u_x, &mut tmp);
        let u_x_supp = (x != base).then(|| Support::from_perm(&u_x));

        for (gen_idx, gen) in tree.generators.iter().enumerate() {
            let y = gen[x] as usize;
            // Skip `x == base` when `gen` fixes `base` (already seeded above).
            if x == base && y == base {
                continue;
            }
            // Skip Schreier tree edges (`u_y^{-1} ∘ gen ∘ u_x == id`).
            if tree.parent_pt[y] as usize == x && tree.parent_gen[y] as usize == gen_idx {
                continue;
            }
            // Disjoint supports commute (`u_x^{-1} ∘ gen ∘ u_x == gen`, already seeded).
            if u_x_supp
                .as_ref()
                .is_some_and(|u_supp| gen_supports[gen_idx].is_disjoint(u_supp))
            {
                continue;
            }

            for v in 0..node_count {
                candidate[v] = gen[u_x[v] as usize];
            }
            tree.apply_left_inv_coset_rep(y, &mut candidate);
            if first_moved(&candidate).is_some() {
                filter.insert(candidate.clone());
            }
        }
    }

    filter.into_generators()
}

/// Compute the minimum-`NodeId` orbit partition induced by `generators` on `0..node_count`.
fn compute_orbits(node_count: usize, generators: &[Vec<u32>]) -> (Vec<NodeId>, usize) {
    let mut parent: Vec<usize> = (0..node_count).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for gen in generators {
        for (v, &img) in gen.iter().enumerate() {
            let rv = find(&mut parent, v);
            let ri = find(&mut parent, img as usize);
            if rv != ri {
                if rv < ri {
                    parent[ri] = rv;
                } else {
                    parent[rv] = ri;
                }
            }
        }
    }
    let mut orbits = Vec::with_capacity(node_count);
    let mut orbit_count = 0;
    for v in 0..node_count {
        let rep = find(&mut parent, v);
        if rep == v {
            orbit_count += 1;
        }
        orbits.push(NodeId(rep as u32));
    }
    (orbits, orbit_count)
}

/// Validate externally supplied generators against `graph` and convert to non-identity `u32` images.
fn validate_generators(
    graph: &Graph,
    generators: &[Vec<NodeId>],
) -> Result<Vec<Vec<u32>>, SymmetryGeneratorError> {
    let node_count = graph.node_count();
    let mut seen = vec![false; node_count];
    let mut result = Vec::new();

    for (generator_index, gen) in generators.iter().enumerate() {
        if gen.len() != node_count {
            return Err(SymmetryGeneratorError::LengthMismatch {
                generator_index,
                len: gen.len(),
                node_count,
            });
        }
        seen.fill(false);
        let mut perm = Vec::with_capacity(node_count);
        let mut is_identity = true;
        for (idx, &node) in gen.iter().enumerate() {
            let target = node.index();
            if target >= node_count {
                return Err(SymmetryGeneratorError::NodeOutOfBounds {
                    generator_index,
                    node,
                    node_count,
                });
            }
            if seen[target] {
                return Err(SymmetryGeneratorError::DuplicateImage {
                    generator_index,
                    image: node,
                });
            }
            seen[target] = true;
            if target != idx {
                is_identity = false;
            }
            perm.push(node.0);
        }
        if is_identity {
            continue;
        }
        for edge in graph.edge_ids() {
            let [first, second] = graph.edge_endpoints(edge);
            let mapped_first = NodeId(perm[first.index()]);
            let mapped_second = NodeId(perm[second.index()]);
            if graph.find_edge(mapped_first, mapped_second).is_none() {
                return Err(SymmetryGeneratorError::EdgeNotPreserved {
                    generator_index,
                    endpoints: [first, second],
                });
            }
        }
        result.push(perm);
    }

    Ok(result)
}

fn generators_to_u32(generators: &[Vec<NodeId>]) -> Vec<Vec<u32>> {
    generators
        .iter()
        .filter_map(|gen| {
            let perm: Vec<u32> = gen.iter().map(|node| node.0).collect();
            first_moved(&perm).map(|_| perm)
        })
        .collect()
}

/// Grochow--Kellis strict-order symmetry-breaking constraints for a query graph.
///
/// Each stored constraint `(lo, hi)` enforces `phi(lo) < phi(hi)` on every emitted embedding `phi`
/// of the query. Constraints are constructed along the deterministic base priority
/// `(Reverse(query.degree(q)), q)`, forming a strict DAG over `0..node_count`.
///
/// # Semantic properties
///
/// - For a fixed `(query, node_color, edge_color)`, [`Self::from_graph`] and
///   [`Self::from_graph_with_scratch`] produce identical constraints under every
///   [`AutomorphismAlgorithm`] variant.
/// - Enforcing every constraint `(lo, hi)` emits exactly one canonical embedding per query
///   symmetry orbit when paired with an injective embedding relation whose predicates are
///   invariant under the represented symmetry group.
#[derive(Clone, Debug)]
pub struct QuerySymmetryConstraints {
    node_count: usize,
    constraints: Vec<(NodeId, NodeId)>,
    priority_order: Vec<u32>,
    /// Non-trivial stabilizer-chain levels indexed by rank in `priority_order`.
    levels: Vec<Option<SchreierTree>>,
    /// Highest rank `r` with `levels[r].is_some()`, or `0` when trivial.
    last_active_rank: usize,
}

impl QuerySymmetryConstraints {
    /// Trivial symmetry constraints (`{1}`, zero inequalities) for a query with `node_count` nodes.
    pub fn none(node_count: usize) -> Self {
        Self {
            node_count,
            constraints: Vec::new(),
            priority_order: (0..node_count as u32).collect(),
            levels: Vec::new(),
            last_active_rank: 0,
        }
    }

    /// Compute Grochow--Kellis symmetry constraints for `query` under `node_color` and
    /// `edge_color` using `alg`.
    pub fn from_graph<NC: Ord, EC: Ord>(
        query: &Graph,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> Self {
        let mut scratch = AutomorphismScratch::new();
        Self::from_graph_with_scratch(query, node_color, edge_color, alg, &mut scratch)
    }

    /// Compute Grochow--Kellis symmetry constraints for `query`, reusing `scratch`.
    pub fn from_graph_with_scratch<NC: Ord, EC: Ord>(
        query: &Graph,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> Self {
        let group = query.automorphisms_with_scratch(node_color, edge_color, alg, scratch);
        let generators = generators_to_u32(group.generators());
        Self::from_validated_generators(query, generators)
    }

    /// Build Grochow--Kellis symmetry constraints for an externally supplied subgroup of
    /// `Aut(query)`, such as a stereo- or reaction-preserving proper subgroup.
    ///
    /// # Errors
    ///
    /// Returns [`SymmetryGeneratorError`] if any generator has wrong length, maps a node out of
    /// bounds, is not a bijection, or fails to preserve an edge of `query`.
    pub fn try_from_generators(
        query: &Graph,
        generators: &[Vec<NodeId>],
    ) -> Result<Self, SymmetryGeneratorError> {
        let validated = validate_generators(query, generators)?;
        Ok(Self::from_validated_generators(query, validated))
    }

    fn from_validated_generators(query: &Graph, initial_generators: Vec<Vec<u32>>) -> Self {
        let node_count = query.node_count();
        let priority_order = query_priority_order(query);
        if initial_generators.is_empty() || node_count <= 1 {
            return Self {
                node_count,
                constraints: Vec::new(),
                priority_order,
                levels: Vec::new(),
                last_active_rank: 0,
            };
        }

        // Sims-filter initial generators so redundant/large generating sets stay compact.
        let mut filter = SimsFilter::new(node_count);
        for gen in initial_generators {
            filter.insert(gen);
        }
        let mut current_gens = filter.into_generators();

        let mut constraints = Vec::new();
        let mut levels: Vec<Option<SchreierTree>> = vec![None; node_count];
        let mut last_active_rank = 0;

        for (rank, &base_u32) in priority_order.iter().enumerate() {
            if current_gens.is_empty() {
                break;
            }
            let base = base_u32 as usize;
            let mut tree = SchreierTree::build(base, node_count, &current_gens);
            if tree.orbit.len() <= 1 {
                continue;
            }
            let mut others: Vec<u32> = tree
                .orbit
                .iter()
                .copied()
                .filter(|&u| u != base_u32)
                .collect();
            others.sort_unstable();
            for u in others {
                constraints.push((NodeId(base_u32), NodeId(u)));
            }
            let next_gens = stabilize_point(&tree, node_count);
            // Sort orbit deterministically (with `base` first) for reproducible traversal.
            tree.orbit.sort_unstable();
            if let Some(pos) = tree.orbit.iter().position(|&u| u == base_u32) {
                tree.orbit.swap(0, pos);
            }
            levels[rank] = Some(tree);
            last_active_rank = rank;
            current_gens = next_gens;
        }

        if constraints.is_empty() {
            levels.clear();
        }

        Self {
            node_count,
            constraints,
            priority_order,
            levels,
            last_active_rank,
        }
    }

    /// Number of nodes in the query graph domain.
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Whether the represented symmetry group is trivial (zero constraints).
    pub fn is_trivial(&self) -> bool {
        self.constraints.is_empty()
    }

    /// Strict precedence pairs `(lo, hi)` enforcing `phi(lo) < phi(hi)`.
    pub fn constraints(&self) -> &[(NodeId, NodeId)] {
        &self.constraints
    }

    #[allow(
        dead_code,
        reason = "consumed by PreparedSubgraphQuery and solvers in S1c/S2"
    )]
    pub(crate) fn priority_order(&self) -> &[u32] {
        &self.priority_order
    }
}

/// Operation-issued or generator-validated symmetry carrier for a target (host) graph.
///
/// Stores the target's orbit partition, bitmasks of minimum orbit representatives, and a compact
/// generating set for exact point-stabilizer pruning during search.
///
/// # Semantic properties
///
/// - [`Self::orbit_of`], [`Self::same_orbit`], [`Self::orbit_count`],
///   [`Self::is_orbit_representative`], and emitted host-reduced canonical embedding sets are
///   invariant across [`AutomorphismAlgorithm`] variants.
/// - Asymmetric targets (`orbit_count == node_count`) satisfy [`Self::is_trivial`] and incur zero
///   search-time stabilizer work.
#[derive(Clone, Debug)]
pub struct TargetSymmetry {
    node_count: usize,
    orbit_count: usize,
    orbits: Vec<NodeId>,
    orbit_reps_u64: u64,
    orbit_reps_u128: u128,
    orbit_reps_words: Vec<u64>,
    generators: Vec<Vec<u32>>,
}

impl TargetSymmetry {
    /// Trivial symmetry (`{1}`, every vertex in its own singleton orbit) for a target with
    /// `node_count` nodes.
    pub fn none(node_count: usize) -> Self {
        let orbits: Vec<NodeId> = (0..node_count as u32).map(NodeId).collect();
        let (orbit_reps_u64, orbit_reps_u128, orbit_reps_words) = build_orbit_rep_masks(&orbits);
        Self {
            node_count,
            orbit_count: node_count,
            orbits,
            orbit_reps_u64,
            orbit_reps_u128,
            orbit_reps_words,
            generators: Vec::new(),
        }
    }

    /// Compute target symmetry for `target` under `node_color` and `edge_color` using `alg`.
    pub fn from_graph<NC: Ord, EC: Ord>(
        target: &Graph,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> Self {
        let mut scratch = AutomorphismScratch::new();
        Self::from_graph_with_scratch(target, node_color, edge_color, alg, &mut scratch)
    }

    /// Compute target symmetry for `target`, reusing `scratch`.
    pub fn from_graph_with_scratch<NC: Ord, EC: Ord>(
        target: &Graph,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> Self {
        let group = target.automorphisms_with_scratch(node_color, edge_color, alg, scratch);
        Self::from_automorphism_group(&group)
    }

    /// Construct target symmetry from an operation-issued [`AutomorphismGroup`].
    pub fn from_automorphism_group(group: &AutomorphismGroup) -> Self {
        let node_count = group.node_count();
        let orbit_count = group.orbit_count();
        let orbits: Vec<NodeId> = (0..node_count as u32)
            .map(|v| group.orbit_of(NodeId(v)))
            .collect();
        let (orbit_reps_u64, orbit_reps_u128, orbit_reps_words) = build_orbit_rep_masks(&orbits);
        let generators = if orbit_count == node_count {
            Vec::new()
        } else {
            let mut filter = SimsFilter::new(node_count);
            for gen in generators_to_u32(group.generators()) {
                filter.insert(gen);
            }
            filter.into_generators()
        };
        Self {
            node_count,
            orbit_count,
            orbits,
            orbit_reps_u64,
            orbit_reps_u128,
            orbit_reps_words,
            generators,
        }
    }

    /// Construct target symmetry from an externally supplied subgroup of `Aut(target)`, such as a
    /// stereo-preserving proper molecular symmetry group.
    ///
    /// # Errors
    ///
    /// Returns [`SymmetryGeneratorError`] if any generator has wrong length, maps a node out of
    /// bounds, is not a bijection, or fails to preserve an edge of `target`.
    pub fn try_from_generators(
        target: &Graph,
        generators: &[Vec<NodeId>],
    ) -> Result<Self, SymmetryGeneratorError> {
        let node_count = target.node_count();
        let validated = validate_generators(target, generators)?;
        if validated.is_empty() {
            return Ok(Self::none(node_count));
        }
        let mut filter = SimsFilter::new(node_count);
        for gen in validated {
            filter.insert(gen);
        }
        let filtered = filter.into_generators();
        let (orbits, orbit_count) = compute_orbits(node_count, &filtered);
        let (orbit_reps_u64, orbit_reps_u128, orbit_reps_words) = build_orbit_rep_masks(&orbits);
        Ok(Self {
            node_count,
            orbit_count,
            orbits,
            orbit_reps_u64,
            orbit_reps_u128,
            orbit_reps_words,
            generators: if orbit_count == node_count {
                Vec::new()
            } else {
                filtered
            },
        })
    }

    /// Number of nodes in the target graph domain.
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Number of distinct vertex orbits under the target symmetry group.
    pub fn orbit_count(&self) -> usize {
        self.orbit_count
    }

    /// Whether the target symmetry group is trivial (`orbit_count == node_count`).
    pub fn is_trivial(&self) -> bool {
        self.orbit_count == self.node_count
    }

    /// Minimum [`NodeId`] in the orbit containing `node`.
    ///
    /// # Panics
    ///
    /// Panics if `node.index() >= self.node_count()`.
    pub fn orbit_of(&self, node: NodeId) -> NodeId {
        self.orbits[node.index()]
    }

    /// Whether `first` and `second` lie in the same target symmetry orbit.
    ///
    /// # Panics
    ///
    /// Panics if either node is outside `0..self.node_count()`.
    pub fn same_orbit(&self, first: NodeId, second: NodeId) -> bool {
        self.orbits[first.index()] == self.orbits[second.index()]
    }

    /// Whether `node` is the minimum [`NodeId`] in its target symmetry orbit.
    ///
    /// # Panics
    ///
    /// Panics if `node.index() >= self.node_count()`.
    pub fn is_orbit_representative(&self, node: NodeId) -> bool {
        self.orbits[node.index()] == node
    }

    #[allow(dead_code, reason = "consumed by BitsetRi and Glasgow solvers in S2")]
    pub(crate) fn orbit_reps_u64(&self) -> u64 {
        self.orbit_reps_u64
    }

    #[allow(dead_code, reason = "consumed by BitsetRi and Glasgow solvers in S2")]
    pub(crate) fn orbit_reps_u128(&self) -> u128 {
        self.orbit_reps_u128
    }

    #[allow(dead_code, reason = "consumed by BitsetRi and Glasgow solvers in S2")]
    pub(crate) fn orbit_reps_words(&self) -> &[u64] {
        &self.orbit_reps_words
    }
}

fn build_orbit_rep_masks(orbits: &[NodeId]) -> (u64, u128, Vec<u64>) {
    let n = orbits.len();
    let mut u64_mask = 0u64;
    let mut u128_mask = 0u128;
    let mut words = vec![0u64; n.div_ceil(64)];
    for (v, &rep) in orbits.iter().enumerate() {
        if rep.index() == v {
            if v < 64 {
                u64_mask |= 1u64 << v;
            }
            if v < 128 {
                u128_mask |= 1u128 << v;
            }
            words[v / 64] |= 1u64 << (v % 64);
        }
    }
    (u64_mask, u128_mask, words)
}

/// One level of a target point-stabilizer chain `S_d = Stab_{G_H}(t_0, ..., t_{d-1})`.
#[derive(Clone, Debug)]
struct TargetStabilizerLevel {
    orbits: Vec<NodeId>,
    generators: Vec<Vec<u32>>,
    /// Schreier tree rooted at `t_d` (populated when `S_d` is non-trivial and pushed past `t_d`).
    tree: Option<SchreierTree>,
}

/// Incremental target point-stabilizer chain supporting both push-down host-symmetry pruning
/// and exact combined `G_Q × G_H` leaf canonicality checks.
#[derive(Clone, Debug, Default)]
pub(crate) struct TargetStabilizerChain {
    node_count: usize,
    levels: Vec<TargetStabilizerLevel>,
}

impl TargetStabilizerChain {
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Initialize level 0 from `target_sym`.
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn reset(&mut self, target_sym: &TargetSymmetry) {
        self.node_count = target_sym.node_count;
        self.levels.clear();
        if !target_sym.is_trivial() {
            self.levels.push(TargetStabilizerLevel {
                orbits: target_sym.orbits.clone(),
                generators: target_sym.generators.clone(),
                tree: None,
            });
        }
    }

    /// Truncate the chain back to `depth` stabilized points (leaving `levels.len() == depth + 1`,
    /// or `0` when initialized from a trivial symmetry group).
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn truncate(&mut self, depth: usize) {
        if self.levels.is_empty() {
            return;
        }
        let keep = (depth + 1).min(self.levels.len());
        self.levels.truncate(keep);
        if let Some(last) = self.levels.last_mut() {
            last.tree = None;
        }
    }

    /// Current stabilizer depth (`d` means `S_d = Stab(t_0, ..., t_{d-1})` is at the top).
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn depth(&self) -> usize {
        self.levels.len().saturating_sub(1)
    }

    /// Whether the current top stabilizer `S_d` is trivial (`{1}`).
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn is_current_trivial(&self) -> bool {
        self.levels
            .last()
            .is_none_or(|level| level.generators.is_empty())
    }

    /// Whether `target_node` is the minimum vertex in its orbit under the current stabilizer `S_d`.
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn is_orbit_min_at_current(&self, target_node: usize) -> bool {
        match self.levels.last() {
            None => true,
            Some(level) => level.orbits[target_node].index() == target_node,
        }
    }

    /// Verify that `target_node` is minimal in its orbit under the current stabilizer `S_d`,
    /// and if so push `S_{d+1} = Stab_{S_d}(target_node)` onto the chain.
    ///
    /// Returns `false` (leaving the chain unchanged) if `target_node` is not an `S_d`-orbit
    /// minimum.
    #[allow(dead_code, reason = "consumed by solvers in S2")]
    pub(crate) fn try_push_orbit_min(&mut self, target_node: usize) -> bool {
        let Some(top) = self.levels.last_mut() else {
            return true;
        };
        if top.generators.is_empty() {
            return true;
        }
        if top.orbits[target_node].index() != target_node {
            return false;
        }
        let tree = SchreierTree::build(target_node, self.node_count, &top.generators);
        let (next_gens, next_orbits) = if tree.orbit.len() <= 1 {
            (top.generators.clone(), top.orbits.clone())
        } else {
            let next_gens = stabilize_point(&tree, self.node_count);
            let (next_orbits, _) = compute_orbits(self.node_count, &next_gens);
            (next_gens, next_orbits)
        };
        top.tree = Some(tree);
        self.levels.push(TargetStabilizerLevel {
            orbits: next_orbits,
            generators: next_gens,
            tree: None,
        });
        true
    }
}

/// Exact leaf canonicality guard for combined `G_Q × G_H` symmetry reduction.
///
/// Returns `true` if and only if `embedding` is the unique lexicographic minimum of its
/// `(G_Q × G_H)`-orbit `{h ∘ embedding ∘ g : g ∈ G_Q, h ∈ G_H}` ordered along `priority_order`.
#[allow(dead_code, reason = "consumed by solvers in S2")]
pub(crate) fn is_canonical_double_orbit_embedding(
    embedding: &[usize],
    priority_order: &[u32],
    query_sym: &QuerySymmetryConstraints,
    target_sym: &TargetSymmetry,
    chain: &mut TargetStabilizerChain,
) -> bool {
    let n1 = embedding.len();
    debug_assert_eq!(n1, query_sym.node_count);
    debug_assert_eq!(n1, priority_order.len());
    debug_assert!(query_sym.is_trivial() || query_sym.priority_order == priority_order);
    if n1 == 0 || target_sym.is_trivial() {
        return true;
    }

    // Precompute `suffix_is_global_min[r]`: true iff `embedding[p_r..]` already consists of the
    // `n1 - r` smallest unused target ids in ascending order, so no injective completion can
    // beat `embedding` from rank `r` onward.
    let mut suffix_is_global_min = vec![false; n1 + 1];
    suffix_is_global_min[n1] = true;
    for r in (0..n1).rev() {
        let t_r = embedding[priority_order[r] as usize];
        let smaller_in_prefix = (0..r)
            .filter(|&i| embedding[priority_order[i] as usize] < t_r)
            .count();
        suffix_is_global_min[r] = suffix_is_global_min[r + 1] && smaller_in_prefix == t_r;
    }
    if suffix_is_global_min[0] {
        return true;
    }

    // Build exact target stabilizer chain S_0 >= S_1 >= ... along `priority_order`.
    chain.reset(target_sym);
    for &p_r in priority_order {
        if chain.is_current_trivial() {
            break;
        }
        if !chain.try_push_orbit_min(embedding[p_r as usize]) {
            return false;
        }
    }

    if query_sym.is_trivial() {
        return true;
    }

    let n2 = target_sym.node_count;
    let g_init: Vec<u32> = (0..n1 as u32).collect();
    let h_init: Vec<u32> = (0..n2 as u32).collect();
    let mut g_tmp = vec![0u32; n1];

    check_double_orbit_recursive(
        0,
        &g_init,
        &h_init,
        embedding,
        query_sym,
        &chain.levels,
        &suffix_is_global_min,
        &mut g_tmp,
    )
}

#[allow(clippy::too_many_arguments)]
fn check_double_orbit_recursive(
    rank: usize,
    g_cur: &[u32],
    h_cur: &[u32],
    embedding: &[usize],
    query_sym: &QuerySymmetryConstraints,
    target_levels: &[TargetStabilizerLevel],
    suffix_is_global_min: &[bool],
    g_tmp: &mut [u32],
) -> bool {
    let n1 = embedding.len();
    if rank == n1 || suffix_is_global_min[rank] {
        return true;
    }

    let priority_order = &query_sym.priority_order;
    let target_level = target_levels.get(rank).filter(|l| !l.generators.is_empty());
    let query_done = rank > query_sym.last_active_rank;

    if query_done && target_level.is_none() {
        for &p_j in &priority_order[rank..] {
            let val = h_cur[embedding[g_cur[p_j as usize] as usize]] as usize;
            let expected = embedding[p_j as usize];
            if val < expected {
                return false;
            }
            if val > expected {
                return true;
            }
        }
        return true;
    }

    let p_r = priority_order[rank] as usize;
    let t_r = embedding[p_r];

    let check_branch = |g_next: &[u32], g_tmp: &mut [u32]| -> bool {
        let w = h_cur[embedding[g_next[p_r] as usize]] as usize;
        match target_level {
            None => {
                if w < t_r {
                    return false;
                }
                if w > t_r {
                    return true;
                }
                check_double_orbit_recursive(
                    rank + 1,
                    g_next,
                    h_cur,
                    embedding,
                    query_sym,
                    target_levels,
                    suffix_is_global_min,
                    g_tmp,
                )
            }
            Some(level) => {
                let w_min = level.orbits[w].index();
                if w_min < t_r {
                    return false;
                }
                if w_min > t_r {
                    return true;
                }
                if w == t_r {
                    check_double_orbit_recursive(
                        rank + 1,
                        g_next,
                        h_cur,
                        embedding,
                        query_sym,
                        target_levels,
                        suffix_is_global_min,
                        g_tmp,
                    )
                } else {
                    let mut h_next = h_cur.to_vec();
                    level
                        .tree
                        .as_ref()
                        .expect("non-trivial target level carries Schreier tree")
                        .apply_left_inv_coset_rep(w, &mut h_next);
                    check_double_orbit_recursive(
                        rank + 1,
                        g_next,
                        &h_next,
                        embedding,
                        query_sym,
                        target_levels,
                        suffix_is_global_min,
                        g_tmp,
                    )
                }
            }
        }
    };

    if let Some(Some(q_tree)) = query_sym.levels.get(rank) {
        let mut g_next = vec![0u32; n1];
        for &u in &q_tree.orbit {
            if u as usize == p_r {
                if !check_branch(g_cur, g_tmp) {
                    return false;
                }
            } else {
                q_tree.compose_right_coset_rep(g_cur, u as usize, &mut g_next, g_tmp);
                if !check_branch(&g_next, g_tmp) {
                    return false;
                }
            }
        }
        true
    } else {
        check_branch(g_cur, g_tmp)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::ops::ControlFlow;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::algorithms::subgraph_isomorphism::SubgraphIsomorphismAlgorithm::Vf2Rdkit;
    use crate::union_find::UnionFind;

    const SIMPLE_GRAPHS_GRAPH6: &str = include_str!("../../../tests/data/simple-through-8.g6");

    fn parse_graph6(source: &str) -> Graph {
        let bytes = source.as_bytes();
        let node_count = (bytes[0] - 63) as usize;
        let mut edges = Vec::new();
        let mut bit = 0;
        for second in 1..node_count {
            for first in 0..second {
                let value = bytes[1 + bit / 6] - 63;
                if value & (1u8 << (5 - bit % 6)) != 0 {
                    edges.push([first as u32, second as u32]);
                }
                bit += 1;
            }
        }
        Graph::new(node_count, &edges)
    }

    fn all_embeddings(target: &Graph, query: &Graph) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        let _: ControlFlow<()> = target.visit_subgraph_isomorphisms(
            query,
            &mut |_, _| true,
            &mut |_, _| true,
            Vf2Rdkit,
            |emb| {
                out.push(emb.iter().map(|n| n.index()).collect());
                ControlFlow::Continue(())
            },
        );
        out.sort_unstable();
        out
    }

    fn count_orbits_by_generators(
        embeddings: &[Vec<usize>],
        query_gens: &[Vec<NodeId>],
        target_gens: &[Vec<NodeId>],
    ) -> (usize, Vec<Vec<Vec<usize>>>) {
        let index: HashMap<&[usize], usize> = embeddings
            .iter()
            .enumerate()
            .map(|(i, e)| (e.as_slice(), i))
            .collect();
        let mut uf = UnionFind::new(embeddings.len());
        let mut buf = Vec::new();

        for (i, emb) in embeddings.iter().enumerate() {
            for g in query_gens {
                buf.clear();
                buf.extend(g.iter().map(|&q| emb[q.index()]));
                let j = index[buf.as_slice()];
                uf.union(i, j);
            }
            for h in target_gens {
                buf.clear();
                buf.extend(emb.iter().map(|&t| h[t].index()));
                let j = index[buf.as_slice()];
                uf.union(i, j);
            }
        }

        let mut groups: HashMap<usize, Vec<Vec<usize>>> = HashMap::new();
        for (i, emb) in embeddings.iter().enumerate() {
            groups.entry(uf.find(i)).or_default().push(emb.clone());
        }
        let mut orbits: Vec<Vec<Vec<usize>>> = groups.into_values().collect();
        for orbit in &mut orbits {
            orbit.sort_unstable();
        }
        orbits.sort_unstable();
        (orbits.len(), orbits)
    }

    fn satisfies_query_constraints(emb: &[usize], sym: &QuerySymmetryConstraints) -> bool {
        sym.constraints()
            .iter()
            .all(|&(lo, hi)| emb[lo.index()] < emb[hi.index()])
    }

    fn satisfies_target_pushdown(
        emb: &[usize],
        priority_order: &[u32],
        target_sym: &TargetSymmetry,
        chain: &mut TargetStabilizerChain,
    ) -> bool {
        if target_sym.is_trivial() {
            return true;
        }
        chain.reset(target_sym);
        for &p_r in priority_order {
            if chain.is_current_trivial() {
                return true;
            }
            if !chain.try_push_orbit_min(emb[p_r as usize]) {
                return false;
            }
        }
        true
    }

    fn priority_key(emb: &[usize], priority_order: &[u32]) -> Vec<usize> {
        priority_order.iter().map(|&q| emb[q as usize]).collect()
    }

    #[rstest]
    #[case::empty(0)]
    #[case::singleton(1)]
    #[case::five(5)]
    fn test_query_symmetry_constraints_none(#[case] node_count: usize) {
        let sym = QuerySymmetryConstraints::none(node_count);
        assert_eq!(sym.node_count(), node_count);
        assert!(sym.is_trivial());
        assert_eq!(sym.constraints(), &[]);
    }

    #[rstest]
    #[case::empty(Graph::new(0, &[]), vec![], vec![], true)]
    #[case::singleton(Graph::new(1, &[]), vec![0], vec![], true)]
    #[case::single_edge(Graph::new(2, &[[0, 1]]), vec![0, 0], vec![0], false)]
    #[case::colored_edge(Graph::new(2, &[[0, 1]]), vec![0, 1], vec![0], true)]
    #[case::path_4(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), vec![0; 4], vec![0; 3], false)]
    #[case::cycle_4(Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), vec![0; 4], vec![0; 4], false)]
    #[case::edge_colored_cycle_4(Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), vec![0; 4], vec![1, 2, 1, 2], false)]
    #[case::cycle_6(Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]), vec![0; 6], vec![0; 6], false)]
    #[case::k4(Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]), vec![0; 4], vec![0; 6], false)]
    #[case::disconnected_edges(Graph::new(4, &[[0, 1], [2, 3]]), vec![0; 4], vec![0, 0], false)]
    fn test_query_symmetry_constraints_from_graph(
        #[case] query: Graph,
        #[case] node_colors: Vec<u8>,
        #[case] edge_colors: Vec<u8>,
        #[case] expected_trivial: bool,
    ) {
        let mut scratch = AutomorphismScratch::new();
        let trail = QuerySymmetryConstraints::from_graph(
            &query,
            |n| node_colors[n.index()],
            |e| edge_colors[e.index()],
            AutomorphismAlgorithm::TrailIr,
        );
        let trail_scratch = QuerySymmetryConstraints::from_graph_with_scratch(
            &query,
            |n| node_colors[n.index()],
            |e| edge_colors[e.index()],
            AutomorphismAlgorithm::TrailIr,
            &mut scratch,
        );
        let nauty = QuerySymmetryConstraints::from_graph(
            &query,
            |n| node_colors[n.index()],
            |e| edge_colors[e.index()],
            AutomorphismAlgorithm::Nauty,
        );

        assert_eq!(trail.node_count(), query.node_count());
        assert_eq!(trail.is_trivial(), expected_trivial);
        assert_eq!(trail.constraints(), trail_scratch.constraints());
        assert_eq!(trail.constraints(), nauty.constraints());

        // Every constraint obeys strict priority DAG order.
        let order = query_priority_order(&query);
        let mut rank = vec![0usize; query.node_count()];
        for (r, &v) in order.iter().enumerate() {
            rank[v as usize] = r;
        }
        for &(lo, hi) in trail.constraints() {
            assert_ne!(lo, hi);
            assert!(lo.index() < query.node_count());
            assert!(hi.index() < query.node_count());
            assert!(rank[lo.index()] < rank[hi.index()]);
        }
    }

    #[rstest]
    fn test_query_symmetry_constraints_try_from_generators() {
        // Rotational subgroup C_4 < D_4 on 4-cycle: 4 elements, 1 constraint orbit {0,1,2,3},
        // stabilizer of 0 in C_4 is trivial (unlike D_4 which also constrains (1, 3)).
        let cycle4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let rot_gen = vec![vec![NodeId(1), NodeId(2), NodeId(3), NodeId(0)]];
        let c4 = QuerySymmetryConstraints::try_from_generators(&cycle4, &rot_gen).unwrap();
        assert_eq!(
            c4.constraints(),
            &[
                (NodeId(0), NodeId(1)),
                (NodeId(0), NodeId(2)),
                (NodeId(0), NodeId(3)),
            ]
        );

        let d4 = QuerySymmetryConstraints::from_graph(
            &cycle4,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        assert_eq!(
            d4.constraints(),
            &[
                (NodeId(0), NodeId(1)),
                (NodeId(0), NodeId(2)),
                (NodeId(0), NodeId(3)),
                (NodeId(1), NodeId(3)),
            ]
        );

        // Empty or identity-only generators yield trivial constraints.
        let id = vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]];
        assert!(QuerySymmetryConstraints::try_from_generators(&cycle4, &id)
            .unwrap()
            .is_trivial());
    }

    #[rstest]
    #[case::length_mismatch(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![NodeId(0), NodeId(1)]],
        SymmetryGeneratorError::LengthMismatch { generator_index: 0, len: 2, node_count: 3 },
    )]
    #[case::out_of_bounds(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![NodeId(0), NodeId(1), NodeId(3)]],
        SymmetryGeneratorError::NodeOutOfBounds { generator_index: 0, node: NodeId(3), node_count: 3 },
    )]
    #[case::duplicate_image(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![NodeId(0), NodeId(1), NodeId(1)]],
        SymmetryGeneratorError::DuplicateImage { generator_index: 0, image: NodeId(1) },
    )]
    #[case::edge_not_preserved(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
        SymmetryGeneratorError::EdgeNotPreserved { generator_index: 0, endpoints: [NodeId(1), NodeId(2)] },
    )]
    fn test_symmetry_try_from_generators_error(
        #[case] graph: Graph,
        #[case] generators: Vec<Vec<NodeId>>,
        #[case] expected: SymmetryGeneratorError,
    ) {
        assert_eq!(
            QuerySymmetryConstraints::try_from_generators(&graph, &generators).unwrap_err(),
            expected,
        );
        assert_eq!(
            TargetSymmetry::try_from_generators(&graph, &generators).unwrap_err(),
            expected,
        );
    }

    #[rstest]
    fn test_target_symmetry_from_graph() {
        let empty = TargetSymmetry::none(0);
        assert!(empty.is_trivial());
        assert_eq!(empty.node_count(), 0);
        assert_eq!(empty.orbit_count(), 0);
        assert_eq!(empty.orbit_reps_u64(), 0);

        let path4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let trail =
            TargetSymmetry::from_graph(&path4, |_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let nauty =
            TargetSymmetry::from_graph(&path4, |_| 0u8, |_| 0u8, AutomorphismAlgorithm::Nauty);
        let group = path4.automorphisms(|_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let from_group = TargetSymmetry::from_automorphism_group(&group);

        for sym in [&trail, &nauty, &from_group] {
            assert!(!sym.is_trivial());
            assert_eq!(sym.node_count(), 4);
            assert_eq!(sym.orbit_count(), 2);
            assert_eq!(sym.orbit_of(NodeId(0)), NodeId(0));
            assert_eq!(sym.orbit_of(NodeId(3)), NodeId(0));
            assert_eq!(sym.orbit_of(NodeId(1)), NodeId(1));
            assert_eq!(sym.orbit_of(NodeId(2)), NodeId(1));
            assert!(sym.same_orbit(NodeId(0), NodeId(3)));
            assert!(!sym.same_orbit(NodeId(0), NodeId(1)));
            assert!(sym.is_orbit_representative(NodeId(0)));
            assert!(sym.is_orbit_representative(NodeId(1)));
            assert!(!sym.is_orbit_representative(NodeId(2)));
            assert!(!sym.is_orbit_representative(NodeId(3)));
            assert_eq!(sym.orbit_reps_u64(), 0b0011);
            assert_eq!(sym.orbit_reps_u128(), 0b0011);
            assert_eq!(sym.orbit_reps_words(), &[0b0011]);
        }

        // Proper subgroup C_4 < D_4 on 4-cycle: transitive (1 orbit), stabilizer of 0 is trivial.
        let cycle4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let c4 = TargetSymmetry::try_from_generators(
            &cycle4,
            &[vec![NodeId(1), NodeId(2), NodeId(3), NodeId(0)]],
        )
        .unwrap();
        assert_eq!(c4.orbit_count(), 1);
        assert_eq!(c4.orbit_reps_u64(), 0b0001);
        let mut chain = TargetStabilizerChain::new();
        chain.reset(&c4);
        assert!(!chain.is_current_trivial());
        assert!(chain.try_push_orbit_min(0));
        assert!(chain.is_current_trivial());

        // Wide target (`n = 144 > 128`): 72 colored disjoint edges -> even nodes are orbit reps.
        let edges: Vec<[u32; 2]> = (0..72).map(|i| [2 * i, 2 * i + 1]).collect();
        let wide = Graph::new(144, &edges);
        let wide_sym = TargetSymmetry::from_graph(
            &wide,
            |n| n.index() / 2,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        assert_eq!(wide_sym.orbit_count(), 72);
        assert_eq!(wide_sym.orbit_reps_u64(), 0x5555_5555_5555_5555u64);
        assert_eq!(
            wide_sym.orbit_reps_u128(),
            0x5555_5555_5555_5555_5555_5555_5555_5555u128
        );
        assert_eq!(
            wide_sym.orbit_reps_words(),
            &[
                0x5555_5555_5555_5555u64,
                0x5555_5555_5555_5555u64,
                0x0000_0000_0000_5555u64,
            ]
        );
    }

    #[rstest]
    fn test_is_canonical_double_orbit_embedding() {
        // H = {0-2, 1-3} with proper subgroup G_H = <(0 3)(1 2)> (orbits {0, 3} and {1, 2});
        // Q = {0-1} with G_Q = <(0 1)>. Both (0, 2) and (1, 3) pass simultaneous push-down,
        // but only (0, 2) is canonical under G_Q × G_H.
        let target = Graph::new(4, &[[0, 2], [1, 3]]);
        let query = Graph::new(2, &[[0, 1]]);
        let p_order = query_priority_order(&query);
        let q_sym = QuerySymmetryConstraints::from_graph(
            &query,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        let h_gen = vec![vec![NodeId(3), NodeId(2), NodeId(1), NodeId(0)]];
        let t_sym = TargetSymmetry::try_from_generators(&target, &h_gen).unwrap();

        let embs = all_embeddings(&target, &query);
        assert_eq!(embs, vec![vec![0, 2], vec![1, 3], vec![2, 0], vec![3, 1]]);

        let mut chain = TargetStabilizerChain::new();
        let pushdown: Vec<Vec<usize>> = embs
            .iter()
            .filter(|e| {
                satisfies_query_constraints(e, &q_sym)
                    && satisfies_target_pushdown(e, &p_order, &t_sym, &mut chain)
            })
            .cloned()
            .collect();
        assert_eq!(pushdown, vec![vec![0, 2], vec![1, 3]]);

        let canonical: Vec<Vec<usize>> = pushdown
            .into_iter()
            .filter(|e| {
                is_canonical_double_orbit_embedding(e, &p_order, &q_sym, &t_sym, &mut chain)
            })
            .collect();
        assert_eq!(canonical, vec![vec![0, 2]]);
    }

    #[rstest]
    fn test_symmetry_large_clique_and_direct_product() {
        // K_10 in K_12: |G_Q| = 10!, |G_H| = 12!. Suffix global-min shortcut handles (0..10)
        // without enumerating 10! branches.
        let mut k10_edges = Vec::new();
        for i in 0..10u32 {
            for j in (i + 1)..10 {
                k10_edges.push([i, j]);
            }
        }
        let mut k12_edges = Vec::new();
        for i in 0..12u32 {
            for j in (i + 1)..12 {
                k12_edges.push([i, j]);
            }
        }
        let k10 = Graph::new(10, &k10_edges);
        let k12 = Graph::new(12, &k12_edges);
        let p_order = query_priority_order(&k10);
        let q_sym = QuerySymmetryConstraints::from_graph(
            &k10,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        let t_sym =
            TargetSymmetry::from_graph(&k12, |_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let mut chain = TargetStabilizerChain::new();
        let emb: Vec<usize> = (0..10).collect();
        assert!(is_canonical_double_orbit_embedding(
            &emb, &p_order, &q_sym, &t_sym, &mut chain
        ));
        let non_canon: Vec<usize> = (1..11).collect();
        assert!(!is_canonical_double_orbit_embedding(
            &non_canon, &p_order, &q_sym, &t_sym, &mut chain
        ));

        // 12 disjoint triangles colored by component: |G| = 6^12, pairwise disjoint supports.
        let mut tri_edges = Vec::new();
        for c in 0..12u32 {
            let b = 3 * c;
            tri_edges.push([b, b + 1]);
            tri_edges.push([b + 1, b + 2]);
            tri_edges.push([b, b + 2]);
        }
        let tris = Graph::new(36, &tri_edges);
        let tri_sym = QuerySymmetryConstraints::from_graph(
            &tris,
            |n| n.index() / 3,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        // Each triangle has 3 constraints: (b, b+1), (b, b+2), (b+1, b+2).
        assert_eq!(tri_sym.constraints().len(), 12 * 3);
    }

    #[rstest]
    fn test_symmetry_simple_through_8_invariance() {
        let mut trail_scratch = AutomorphismScratch::new();
        let mut nauty_scratch = AutomorphismScratch::new();
        let mut checked = 0usize;

        for line in SIMPLE_GRAPHS_GRAPH6.lines() {
            let graph = parse_graph6(line);
            checked += 1;

            let q_trail = QuerySymmetryConstraints::from_graph_with_scratch(
                &graph,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::TrailIr,
                &mut trail_scratch,
            );
            let q_nauty = QuerySymmetryConstraints::from_graph_with_scratch(
                &graph,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::Nauty,
                &mut nauty_scratch,
            );
            assert_eq!(q_trail.constraints(), q_nauty.constraints());
            assert_eq!(q_trail.is_trivial(), q_nauty.is_trivial());

            let t_trail = TargetSymmetry::from_graph_with_scratch(
                &graph,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::TrailIr,
                &mut trail_scratch,
            );
            let t_nauty = TargetSymmetry::from_graph_with_scratch(
                &graph,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::Nauty,
                &mut nauty_scratch,
            );
            assert_eq!(t_trail.orbit_count(), t_nauty.orbit_count());
            assert_eq!(t_trail.orbit_reps_u64(), t_nauty.orbit_reps_u64());
            for v in graph.node_ids() {
                assert_eq!(t_trail.orbit_of(v), t_nauty.orbit_of(v));
            }
        }
        // All 1 + 2 + 4 + 11 + 34 + 156 + 1,044 + 12,346 = 13,598 simple graphs through order 8.
        assert_eq!(checked, 13_598);
    }

    #[rstest]
    fn test_symmetry_orbit_reduction_simple_through_5() {
        let graphs: Vec<Graph> = SIMPLE_GRAPHS_GRAPH6
            .lines()
            .map(parse_graph6)
            .filter(|g| g.node_count() <= 5)
            .collect();
        // 1 + 2 + 4 + 11 + 34 = 52 simple graphs through order 5.
        assert_eq!(graphs.len(), 52);

        let mut scratch = AutomorphismScratch::new();
        let groups: Vec<AutomorphismGroup> = graphs
            .iter()
            .map(|g| {
                g.automorphisms_with_scratch(
                    |_| 0u8,
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                    &mut scratch,
                )
            })
            .collect();
        let q_syms: Vec<QuerySymmetryConstraints> = graphs
            .iter()
            .map(|g| {
                QuerySymmetryConstraints::from_graph_with_scratch(
                    g,
                    |_| 0u8,
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                    &mut scratch,
                )
            })
            .collect();
        let t_syms: Vec<TargetSymmetry> = groups
            .iter()
            .map(TargetSymmetry::from_automorphism_group)
            .collect();

        let mut chain = TargetStabilizerChain::new();
        for (qi, query) in graphs.iter().enumerate() {
            let p_order = q_syms[qi].priority_order();
            for (hi, target) in graphs.iter().enumerate() {
                if query.node_count() > target.node_count()
                    || query.edge_count() > target.edge_count()
                {
                    continue;
                }
                let embs = all_embeddings(target, query);
                if embs.is_empty() {
                    continue;
                }

                // Mode 1: G_Q alone — push-down selects exact 1 rep per G_Q-orbit.
                let (q_orbit_count, q_orbits) =
                    count_orbits_by_generators(&embs, groups[qi].generators(), &[]);
                let q_kept: Vec<Vec<usize>> = embs
                    .iter()
                    .filter(|e| satisfies_query_constraints(e, &q_syms[qi]))
                    .cloned()
                    .collect();
                assert_eq!(q_kept.len(), q_orbit_count, "G_Q count ({qi}, {hi})");
                for orbit in &q_orbits {
                    let reps = orbit
                        .iter()
                        .filter(|e| satisfies_query_constraints(e, &q_syms[qi]))
                        .count();
                    assert_eq!(reps, 1, "G_Q unique rep ({qi}, {hi})");
                }

                // Mode 2: G_H alone — push-down selects exact 1 lex-min rep per G_H-orbit.
                let (h_orbit_count, h_orbits) =
                    count_orbits_by_generators(&embs, &[], groups[hi].generators());
                let h_kept: Vec<Vec<usize>> = embs
                    .iter()
                    .filter(|e| satisfies_target_pushdown(e, p_order, &t_syms[hi], &mut chain))
                    .cloned()
                    .collect();
                assert_eq!(h_kept.len(), h_orbit_count, "G_H count ({qi}, {hi})");
                for orbit in &h_orbits {
                    let expected_min = orbit
                        .iter()
                        .min_by_key(|e| priority_key(e, p_order))
                        .unwrap();
                    let reps: Vec<&Vec<usize>> = orbit
                        .iter()
                        .filter(|e| satisfies_target_pushdown(e, p_order, &t_syms[hi], &mut chain))
                        .collect();
                    assert_eq!(reps, vec![expected_min], "G_H lex-min rep ({qi}, {hi})");
                }

                // Mode 3: G_Q × G_H — simultaneous push-down + leaf guard selects exact 1 lex-min
                // rep per double-orbit.
                let (qh_orbit_count, qh_orbits) = count_orbits_by_generators(
                    &embs,
                    groups[qi].generators(),
                    groups[hi].generators(),
                );
                let qh_kept: Vec<Vec<usize>> = embs
                    .iter()
                    .filter(|e| {
                        satisfies_query_constraints(e, &q_syms[qi])
                            && satisfies_target_pushdown(e, p_order, &t_syms[hi], &mut chain)
                            && is_canonical_double_orbit_embedding(
                                e,
                                p_order,
                                &q_syms[qi],
                                &t_syms[hi],
                                &mut chain,
                            )
                    })
                    .cloned()
                    .collect();
                assert_eq!(qh_kept.len(), qh_orbit_count, "G_Q×G_H count ({qi}, {hi})");
                for orbit in &qh_orbits {
                    let expected_min = orbit
                        .iter()
                        .min_by_key(|e| priority_key(e, p_order))
                        .unwrap();
                    assert!(
                        satisfies_query_constraints(expected_min, &q_syms[qi])
                            && satisfies_target_pushdown(
                                expected_min,
                                p_order,
                                &t_syms[hi],
                                &mut chain
                            ),
                        "push-down dropped canonical rep ({qi}, {hi})"
                    );
                    let reps: Vec<&Vec<usize>> = orbit
                        .iter()
                        .filter(|e| {
                            satisfies_query_constraints(e, &q_syms[qi])
                                && satisfies_target_pushdown(e, p_order, &t_syms[hi], &mut chain)
                                && is_canonical_double_orbit_embedding(
                                    e,
                                    p_order,
                                    &q_syms[qi],
                                    &t_syms[hi],
                                    &mut chain,
                                )
                        })
                        .collect();
                    assert_eq!(reps, vec![expected_min], "G_Q×G_H lex-min rep ({qi}, {hi})");
                }
            }
        }
    }
}
