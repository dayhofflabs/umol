//! Aggregate canonicalization inputs and failures.

use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::{iter, mem};

use thiserror::Error;
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismOutput, EdgeId, Graph, GraphRemapping, NodeId, Remapping,
    SubdivisionNodeSource, UnionFind,
};

use super::delta::Delta;
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
    StereoBondId,
};
use super::incidence::{Incidence, IncidenceGraph};
use super::molecule::Molecule;
use super::reaction::Reaction;
use super::reaction_span::ReactionSpan;
use super::remap::MoleculeRemapping;
use super::traits::Reframe;

mod molecule;
mod reaction_span;

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeContext {
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum DescriptionLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}

fn delta_canonicalize_level(delta: &Delta) -> DescriptionLevel {
    match delta {
        Delta::Atom(_) | Delta::Bond(_) => DescriptionLevel::Topology,
        Delta::DativeBond(_)
        | Delta::AromaticSystem(_)
        | Delta::MulticenterBond(_)
        | Delta::NoncovalentBond(_) => DescriptionLevel::Constitution,
        Delta::StereoAtom(_) | Delta::StereoBond(_) => DescriptionLevel::Structure,
        Delta::Constraint(_) => DescriptionLevel::Full,
    }
}

fn molecule_canonicalize_level(molecule: &Molecule) -> DescriptionLevel {
    if molecule.has_constraints() {
        DescriptionLevel::Full
    } else if molecule.has_stereo_atoms() || molecule.has_stereo_bonds() {
        DescriptionLevel::Structure
    } else if molecule.has_dative_bonds()
        || molecule.has_aromatic_systems()
        || molecule.has_multicenter_bonds()
        || molecule.has_noncovalent_bonds()
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
}

fn reaction_canonicalize_level(reaction: &Reaction) -> DescriptionLevel {
    reaction.deltas().iter().fold(
        molecule_canonicalize_level(reaction.lhs()),
        |level, delta| level.max(delta_canonicalize_level(delta)),
    )
}

fn reaction_span_canonicalize_level(span: &ReactionSpan) -> DescriptionLevel {
    if !span.constraints().is_empty() {
        DescriptionLevel::Full
    } else if span.stereo_atoms().count() > 0 || span.stereo_bonds().count() > 0 {
        DescriptionLevel::Structure
    } else if span.dative_bonds().count() > 0
        || span.aromatic_systems().count() > 0
        || span.multicenter_bonds().count() > 0
        || span.noncovalent_bonds().count() > 0
    {
        DescriptionLevel::Constitution
    } else {
        DescriptionLevel::Topology
    }
}

/// Canonical entity-frame selection for complete indexed graph-IR aggregates.
///
/// Canonicalization extends [`Reframe`]: normalization is the fixed-frame prefix, reframing adds
/// participant-frame selection, and canonicalization additionally selects entity ids. It preserves
/// the represented aggregate, transports every reference and position-sensitive value, and
/// normalizes carried forms in the selected frame. Construction of the canonical form is fallible;
/// canonical equality totalizes those failures according to the aggregate's semantic contract.
///
/// # Semantic properties
///
/// For every integrity-valid aggregate and fixed context:
///
/// - successful complete canonicalization is exactly idempotent and invariant under valid dense
///   entity remapping;
/// - `canonical_eq` is reflexive, symmetric, and transitive under its documented failure
///   totalization;
/// - successful canonical hashes are invariant under valid dense entity remapping, and canonical
///   equality implies equal canonical hashes when both hash operations succeed.
///
/// For a fixed umol release and context, canonicalization is deterministic. During the 0.x series,
/// the order of keys and the resulting representatives may change between releases. Returned
/// canonical forms are ordinary IR values without version provenance and must not be used as
/// persistent cross-release identifiers.
pub trait Canonicalize: Reframe {
    type Error;

    /// Construct the complete canonical form.
    ///
    /// # Errors
    ///
    /// Returns the aggregate-specific canonicalization error when intrinsic normalization is
    /// unsatisfiable or a carried attribute cannot be transported in its stored frame.
    fn canonicalize(self, context: &CanonicalizeContext) -> Result<Self, Self::Error>;

    /// Construct the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping maps every entity id in the input frame to its id in the returned
    /// canonical frame. Each of its eight entity-kind components is total on both sides and
    /// therefore represents a dense bijection. For [`Reaction`], the two frames are the complete
    /// union frames of the materialized input and returned reaction spans.
    ///
    /// # Errors
    ///
    /// Returns the same aggregate-specific canonicalization error as
    /// [`Self::canonicalize`]. A [`Reaction`] also reports failure to materialize its reaction span.
    ///
    /// # Semantic properties
    ///
    /// Discarding the remapping yields exactly [`Self::canonicalize`] under the same context.
    /// Remapping the input through the remapping and then applying [`Reframe::reframe`]
    /// reconstructs the returned canonical value.
    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error>;

    /// Hash the complete canonical form.
    ///
    /// The returned value uses Rust's [`DefaultHasher`]. It may collide and is not a persistent
    /// identifier: both the hasher and umol's canonical representatives may change between
    /// releases.
    ///
    /// # Errors
    ///
    /// Returns the same error as [`Self::canonicalize`].
    ///
    /// # Semantic properties
    ///
    /// Successful hashes are invariant under valid dense entity remapping. Canonically equal
    /// values have equal canonical hashes when both hash operations succeed.
    fn canonical_hash(&self, context: &CanonicalizeContext) -> Result<u64, Self::Error>;

    /// Compare complete canonical forms.
    ///
    /// Structural identity short-circuits the operation. Otherwise, two intrinsic contradictions
    /// compare equal.
    fn canonical_eq(&self, other: &Self, context: &CanonicalizeContext) -> bool;
}

fn hash_value<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct InitialColors {
    entities: Vec<u32>,
    incidences: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AutomorphismColor {
    Entity(u32),
    Incidence(u32),
}

/// Colored simple graph and source-domain projection used by canonical search.
///
/// Source nodes occupy the leading `source_node_count` positions. Remaining vertices encode
/// subdivided incidences or bonds and are excluded when backend automorphisms are projected onto
/// the search domain.
#[derive(Clone, Debug)]
struct AutomorphismAdapter {
    graph: Graph,
    colors: Vec<AutomorphismColor>,
    node_sources: Vec<SubdivisionNodeSource>,
    source_node_count: usize,
}

/// Automorphism data over both the source-entity and complete adapter-node domains.
#[derive(Clone, Debug, PartialEq, Eq)]
struct AutomorphismAdapterOutput {
    source_orbits: Vec<NodeId>,
    /// Source entity nodes in backend canonical-label order, used only to order search branches.
    source_canonical_labels: Vec<NodeId>,
    source_generators: Vec<Vec<NodeId>>,
}

impl AutomorphismAdapter {
    /// Construct the exact graph adapter used by canonicalization.
    fn new(incidence_graph: &IncidenceGraph, initial_colors: &InitialColors) -> Self {
        let source = incidence_graph.graph();
        debug_assert_eq!(initial_colors.entities.len(), source.node_count());
        debug_assert_eq!(initial_colors.incidences.len(), source.edge_count());

        let mut colors = initial_colors
            .entities
            .iter()
            .copied()
            .map(AutomorphismColor::Entity)
            .collect::<Vec<_>>();
        let mut node_sources = source
            .node_ids()
            .map(SubdivisionNodeSource::Node)
            .collect::<Vec<_>>();
        let mut edges = Vec::new();

        let direct_pair_counts = source
            .edge_ids()
            .filter(|&edge| {
                matches!(
                    incidence_graph.incidence(edge),
                    Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
                )
            })
            .fold(BTreeMap::<[NodeId; 2], usize>::new(), |mut counts, edge| {
                *counts.entry(source.edge_endpoints(edge)).or_default() += 1;
                counts
            });

        let mut push_edge = |endpoints: [NodeId; 2]| {
            edges.push([endpoints[0].0, endpoints[1].0]);
        };
        for edge in source.edge_ids() {
            let endpoints = source.edge_endpoints(edge);
            let direct = matches!(
                incidence_graph.incidence(edge),
                Incidence::BondEndpoint | Incidence::NoncovalentEndpoint
            ) && direct_pair_counts[&endpoints] == 1;
            if direct {
                push_edge(endpoints);
                continue;
            }

            let occurrence = NodeId(node_sources.len() as u32);
            node_sources.push(SubdivisionNodeSource::Edge(edge));
            colors.push(AutomorphismColor::Incidence(
                initial_colors.incidences[edge.index()],
            ));
            push_edge([endpoints[0], occurrence]);
            push_edge([occurrence, endpoints[1]]);
        }
        let graph = Graph::new(node_sources.len(), &edges);
        debug_assert!(graph.is_simple());

        Self {
            graph,
            colors,
            node_sources,
            source_node_count: source.node_count(),
        }
    }

    fn graph(&self) -> &Graph {
        &self.graph
    }

    fn color(&self, node: NodeId) -> AutomorphismColor {
        self.colors[node.index()]
    }

    fn node_source(&self, node: NodeId) -> SubdivisionNodeSource {
        self.node_sources[node.index()]
    }

    fn automorphisms_for_partition(
        &self,
        partition: &OrderedPartition,
        algorithm: AutomorphismAlgorithm,
    ) -> AutomorphismAdapterOutput {
        let cell_indices = partition.cell_indices(self.graph().node_count());
        // Search partitions may deliberately coarsen covariant occurrence data. Retaining the
        // exact adapter color here keeps orbit pruning restricted to true automorphisms.
        let output = self.graph().automorphisms(
            |node| (self.color(node), cell_indices[node.index()]),
            algorithm,
        );

        self.project_automorphisms(&output)
    }

    fn project_automorphisms(&self, output: &AutomorphismOutput) -> AutomorphismAdapterOutput {
        let source_node = |node| match self.node_source(node) {
            SubdivisionNodeSource::Node(source) if source.index() < self.source_node_count => {
                source
            }
            SubdivisionNodeSource::Node(_) => {
                unreachable!("disjoint colors preserve the source-node domain")
            }
            SubdivisionNodeSource::Edge(_) => {
                unreachable!("disjoint colors preserve the adapter node domain")
            }
        };
        let source_orbits = (0..self.source_node_count)
            .map(|index| source_node(output.orbit_of(NodeId(index as u32))))
            .collect();
        let source_canonical_labels = output
            .canonical_labels()
            .iter()
            .filter_map(|&node| match self.node_source(node) {
                SubdivisionNodeSource::Node(source) if source.index() < self.source_node_count => {
                    Some(source)
                }
                SubdivisionNodeSource::Node(_) => None,
                SubdivisionNodeSource::Edge(_) => None,
            })
            .collect();
        let source_generators = output
            .generators()
            .iter()
            .map(|generator| {
                generator[..self.source_node_count]
                    .iter()
                    .copied()
                    .map(source_node)
                    .collect()
            })
            .collect();
        AutomorphismAdapterOutput {
            source_orbits,
            source_canonical_labels,
            source_generators,
        }
    }
}

fn source_orbits_from_generators(
    source_node_count: usize,
    generators: &[Vec<NodeId>],
) -> Vec<NodeId> {
    let mut sets = UnionFind::new(source_node_count);
    for generator in generators {
        debug_assert_eq!(generator.len(), source_node_count);
        for (source, image) in generator.iter().enumerate() {
            sets.union(source, image.index());
        }
    }

    let roots = (0..source_node_count)
        .map(|source| sets.find(source))
        .collect::<Vec<_>>();
    let mut representatives = vec![usize::MAX; source_node_count];
    for (source, root) in roots.iter().copied().enumerate() {
        representatives[root] = representatives[root].min(source);
    }
    roots
        .into_iter()
        .map(|root| NodeId(representatives[root] as u32))
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OrderedPartition {
    cells: Vec<Vec<NodeId>>,
}

/// Minimum carrier size at which refinement tracks only cells created by the preceding split.
const ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT: usize = 32;

impl OrderedPartition {
    fn from_descriptors<C: Clone + Ord>(descriptors: &[C]) -> Self {
        let mut cells = BTreeMap::<C, Vec<NodeId>>::new();
        for (index, descriptor) in descriptors.iter().enumerate() {
            cells
                .entry(descriptor.clone())
                .or_default()
                .push(NodeId(index as u32));
        }

        Self {
            cells: cells.into_values().collect(),
        }
    }

    fn from_color_ranks(colors: impl IntoIterator<Item = u32>) -> Self {
        let mut cells = Vec::<Vec<NodeId>>::new();
        for (index, color) in colors.into_iter().enumerate() {
            let color = color as usize;
            if cells.len() <= color {
                cells.resize_with(color + 1, Vec::new);
            }
            cells[color].push(NodeId(index as u32));
        }
        cells.retain(|cell| !cell.is_empty());
        Self { cells }
    }

    fn refine(self, graph: &Graph) -> Self {
        if graph.node_count() < ACTIVE_CELL_REFINEMENT_MIN_NODE_COUNT {
            self.refine_by_all_cells(graph)
        } else {
            self.refine_by_active_cells(graph)
        }
    }

    fn refine_by_active_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut active_cells = vec![true; self.cells.len()];
        let mut node_active_columns = vec![usize::MAX; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());
        let mut next_active_cells = Vec::<bool>::with_capacity(self.cells.len());

        loop {
            node_active_columns.fill(usize::MAX);
            let mut active_count = 0;
            for (cell_index, cell) in self.cells.iter().enumerate() {
                if active_cells[cell_index] {
                    for node in cell {
                        node_active_columns[node.index()] = active_count;
                    }
                    active_count += 1;
                }
            }
            signatures.clear();
            signatures.resize(node_count * active_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * active_count;
                for neighbor in graph.neighbors(node) {
                    let column = node_active_columns[neighbor.node.index()];
                    if column != usize::MAX {
                        signatures[offset + column] += 1;
                    }
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * active_count;
                    &signatures[offset..offset + active_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                    next_active_cells.push(false);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    let new_cell_count = next_cells.len() - next_active_cells.len();
                    next_active_cells.extend(iter::repeat_n(true, new_cell_count));
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
            if !changed {
                return self;
            }
            mem::swap(&mut active_cells, &mut next_active_cells);
            next_active_cells.clear();
        }
    }

    fn refine_by_all_cells(mut self, graph: &Graph) -> Self {
        let node_count = graph.node_count();
        let mut cell_indices = vec![0; node_count];
        let mut signatures = Vec::<u32>::new();
        let mut next_cells = Vec::<Vec<NodeId>>::with_capacity(self.cells.len());

        loop {
            for (cell_index, cell) in self.cells.iter().enumerate() {
                for node in cell {
                    cell_indices[node.index()] = cell_index as u32;
                }
            }
            let cell_count = self.cells.len();
            signatures.clear();
            signatures.resize(node_count * cell_count, 0);
            for node in graph.node_ids() {
                let offset = node.index() * cell_count;
                for neighbor in graph.neighbors(node) {
                    signatures[offset + cell_indices[neighbor.node.index()] as usize] += 1;
                }
            }

            let mut changed = false;
            for mut cell in self.cells.drain(..) {
                let signature = |node: NodeId| {
                    let offset = node.index() * cell_count;
                    &signatures[offset..offset + cell_count]
                };
                cell.sort_unstable_by(|lhs, rhs| {
                    signature(*rhs)
                        .cmp(signature(*lhs))
                        .then_with(|| lhs.cmp(rhs))
                });

                let mut group_start = 0;
                for index in 1..cell.len() {
                    if signature(cell[index]) != signature(cell[group_start]) {
                        next_cells.push(cell[group_start..index].to_vec());
                        group_start = index;
                    }
                }
                if group_start == 0 {
                    next_cells.push(cell);
                } else {
                    next_cells.push(cell[group_start..].to_vec());
                    changed = true;
                }
            }

            mem::swap(&mut self.cells, &mut next_cells);
            if !changed {
                return self;
            }
        }
    }

    fn individualize(&self, cell_index: usize, node: NodeId) -> Self {
        let mut cells = Vec::with_capacity(self.cells.len() + 1);
        for (index, cell) in self.cells.iter().enumerate() {
            if index != cell_index {
                cells.push(cell.clone());
                continue;
            }

            cells.push(vec![node]);
            let remainder = cell
                .iter()
                .copied()
                .filter(|&candidate| candidate != node)
                .collect::<Vec<_>>();
            if !remainder.is_empty() {
                cells.push(remainder);
            }
        }

        Self { cells }
    }

    fn first_non_singleton_entity_cell(&self, entity_count: usize) -> Option<usize> {
        self.cells.iter().position(|cell| {
            cell.len() > 1 && cell.first().is_some_and(|node| node.index() < entity_count)
        })
    }

    fn entity_order(&self, entity_count: usize) -> Vec<NodeId> {
        self.cells
            .iter()
            .flatten()
            .copied()
            .filter(|node| node.index() < entity_count)
            .collect()
    }

    fn cell_indices(&self, node_count: usize) -> Vec<u32> {
        let mut indices = vec![0; node_count];
        for (cell_index, cell) in self.cells.iter().enumerate() {
            for node in cell {
                indices[node.index()] = cell_index as u32;
            }
        }
        indices
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BranchOrdering {
    BackendCanonical,
    #[cfg(test)]
    NodeId,
    #[cfg(test)]
    ReverseNodeId,
}

impl BranchOrdering {
    fn requires_automorphisms(self) -> bool {
        self == Self::BackendCanonical
    }

    fn order(
        self,
        adapter: &AutomorphismAdapter,
        automorphisms: Option<&AutomorphismAdapterOutput>,
        candidates: &mut [NodeId],
    ) {
        match self {
            Self::BackendCanonical => {
                let labels = &automorphisms
                    .expect("backend branch order requested automorphisms")
                    .source_canonical_labels;
                let mut ranks = vec![0; adapter.source_node_count];
                for (rank, node) in labels.iter().enumerate() {
                    ranks[node.index()] = rank;
                }
                candidates.sort_unstable_by_key(|node| ranks[node.index()]);
            }
            #[cfg(test)]
            Self::NodeId => candidates.sort_unstable(),
            #[cfg(test)]
            Self::ReverseNodeId => candidates.sort_unstable_by(|lhs, rhs| rhs.cmp(lhs)),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct CanonicalSearchOptions {
    automorphism_pruning: bool,
    branch_order: BranchOrdering,
}

#[cfg(test)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct CanonicalSearchStats {
    initial_residual_cell_sizes: Vec<usize>,
    refinement_calls: usize,
    branch_order_calls: usize,
    backend_calls: usize,
    visited_leaves: usize,
    leaf_comparisons: usize,
    orbit_pruned_branches: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalCandidate<K> {
    key: K,
    entity_order: Vec<NodeId>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalSearchResult<K> {
    candidate: CanonicalCandidate<K>,
    #[cfg(test)]
    stats: CanonicalSearchStats,
}

/// Minimize a typed leaf candidate over an exact partition ordered by semantic descriptors.
///
/// Adapter colors are opaque equality labels and do not order the partition. Automorphism pruning
/// requires the leaf key to be invariant under adapter automorphisms.
fn canonical_search<K, Descriptor, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        options,
        leaf_candidate,
        false,
        &|_| {},
    )
}

fn canonical_search_from_partition<K, LeafCandidate>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
{
    canonical_search_impl(
        adapter,
        partition,
        algorithm,
        options,
        leaf_candidate,
        false,
        &|_| {},
    )
}

fn canonical_search_with_generator_filter<K, Descriptor, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition_descriptors: &[Descriptor],
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    Descriptor: Clone + Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    canonical_search_impl(
        adapter,
        OrderedPartition::from_descriptors(partition_descriptors),
        algorithm,
        options,
        leaf_candidate,
        true,
        filter_generators,
    )
}

fn canonical_search_impl<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    initial: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
) -> CanonicalSearchResult<K>
where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let initial = initial.refine(adapter.graph());
    let mut best = None;
    #[cfg(test)]
    let mut stats = CanonicalSearchStats {
        initial_residual_cell_sizes: initial
            .cells
            .iter()
            .filter(|cell| {
                cell.len() > 1
                    && cell
                        .first()
                        .is_some_and(|node| node.index() < adapter.source_node_count)
            })
            .map(Vec::len)
            .collect(),
        refinement_calls: 1,
        ..Default::default()
    };

    search_partition(
        adapter,
        initial,
        algorithm,
        options,
        leaf_candidate,
        apply_generator_filter,
        filter_generators,
        &mut best,
        #[cfg(test)]
        &mut stats,
    );

    CanonicalSearchResult {
        candidate: best.expect("every finite partition has a discrete entity refinement"),
        #[cfg(test)]
        stats,
    }
}

#[allow(clippy::too_many_arguments)]
fn search_partition<K, LeafCandidate, FilterGenerators>(
    adapter: &AutomorphismAdapter,
    partition: OrderedPartition,
    algorithm: AutomorphismAlgorithm,
    options: CanonicalSearchOptions,
    leaf_candidate: &LeafCandidate,
    apply_generator_filter: bool,
    filter_generators: &FilterGenerators,
    best: &mut Option<CanonicalCandidate<K>>,
    #[cfg(test)] stats: &mut CanonicalSearchStats,
) where
    K: Ord,
    LeafCandidate: Fn(&[NodeId]) -> CanonicalCandidate<K>,
    FilterGenerators: Fn(&mut Vec<Vec<NodeId>>),
{
    let Some(cell_index) = partition.first_non_singleton_entity_cell(adapter.source_node_count)
    else {
        #[cfg(test)]
        {
            stats.visited_leaves += 1;
        }
        let entity_order = partition.entity_order(adapter.source_node_count);
        let candidate = leaf_candidate(&entity_order);
        let improves = if let Some(best) = best.as_ref() {
            #[cfg(test)]
            {
                stats.leaf_comparisons += 1;
            }
            candidate.key < best.key
        } else {
            true
        };
        if improves {
            *best = Some(candidate);
        }
        return;
    };

    let mut candidates = partition.cells[cell_index].clone();
    let automorphisms = (options.automorphism_pruning
        || apply_generator_filter
        || options.branch_order.requires_automorphisms())
    .then(|| {
        #[cfg(test)]
        {
            stats.backend_calls += 1;
        }
        let mut output = adapter.automorphisms_for_partition(&partition, algorithm);
        if apply_generator_filter {
            filter_generators(&mut output.source_generators);
            output.source_orbits =
                source_orbits_from_generators(adapter.source_node_count, &output.source_generators);
        }
        output
    });

    if options.automorphism_pruning {
        let orbits = &automorphisms
            .as_ref()
            .expect("automorphisms requested for orbit pruning")
            .source_orbits;
        let mut representatives = BTreeMap::<NodeId, NodeId>::new();
        for candidate in candidates {
            representatives
                .entry(orbits[candidate.index()])
                .and_modify(|representative| *representative = (*representative).min(candidate))
                .or_insert(candidate);
        }
        #[cfg(test)]
        {
            stats.orbit_pruned_branches +=
                partition.cells[cell_index].len() - representatives.len();
        }
        candidates = representatives.into_values().collect();
    }

    #[cfg(test)]
    {
        stats.branch_order_calls += 1;
    }
    options
        .branch_order
        .order(adapter, automorphisms.as_ref(), &mut candidates);

    for candidate in candidates {
        #[cfg(test)]
        {
            stats.refinement_calls += 1;
        }
        let child = partition
            .individualize(cell_index, candidate)
            .refine(adapter.graph());
        search_partition(
            adapter,
            child,
            algorithm,
            options,
            leaf_candidate,
            apply_generator_filter,
            filter_generators,
            best,
            #[cfg(test)]
            stats,
        );
    }
}

/// A leaf key that is built when it is first compared or read.
struct DeferredKey<'a, K, Build> {
    order: Vec<NodeId>,
    build: &'a Build,
    key: OnceCell<K>,
}

impl<'a, K, Build: Fn(&[NodeId]) -> K> DeferredKey<'a, K, Build> {
    fn new(order: &[NodeId], build: &'a Build) -> Self {
        Self {
            order: order.to_vec(),
            build,
            key: OnceCell::new(),
        }
    }

    fn get(&self) -> &K {
        self.key.get_or_init(|| (self.build)(&self.order))
    }

    /// The key, if a comparison built it.
    fn into_built(self) -> Option<K> {
        self.key.into_inner()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialEq for DeferredKey<'_, K, Build> {
    fn eq(&self, other: &Self) -> bool {
        self.get() == other.get()
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Eq for DeferredKey<'_, K, Build> {}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> PartialOrd for DeferredKey<'_, K, Build> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: Ord, Build: Fn(&[NodeId]) -> K> Ord for DeferredKey<'_, K, Build> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.get().cmp(other.get())
    }
}

/// Topology search input with atoms and non-modal bond subdivisions as vertices.
///
/// Bonds in the largest normalized bond-form class are represented directly as edges. Every other
/// bond remains a colored subdivision vertex. The complete typed leaf key and remapping still
/// order and map every bond.
#[derive(Clone, Debug)]
struct TopologySearchInput {
    adapter: AutomorphismAdapter,
    initial_partition: OrderedPartition,
}

/// Assign dense integers to the colors of one entity kind or of the incidences, in the order of the
/// colors and starting at `first`. Each color comes with the position of its integer in `ranks`.
/// The integer after the last one used is returned.
fn rank_colors<T: Ord>(mut colors: Vec<(T, usize)>, first: u32, ranks: &mut [u32]) -> u32 {
    colors.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    let mut rank = first;
    let mut previous = None;
    for (color, position) in &colors {
        if previous.is_some_and(|previous| previous != color) {
            rank += 1;
        }
        ranks[*position] = rank;
        previous = Some(color);
    }
    if colors.is_empty() {
        first
    } else {
        rank + 1
    }
}

fn molecule_counts(molecule: &Molecule) -> [usize; 8] {
    [
        molecule.atoms().count(),
        molecule.bonds().count(),
        molecule.dative_bonds().count(),
        molecule.aromatic_systems().count(),
        molecule.multicenter_bonds().count(),
        molecule.noncovalent_bonds().count(),
        molecule.stereo_atoms().count(),
        molecule.stereo_bonds().count(),
    ]
}

fn reaction_span_counts(span: &ReactionSpan) -> [usize; 8] {
    [
        span.atoms().len(),
        span.bonds().len(),
        span.dative_bonds().count(),
        span.aromatic_systems().count(),
        span.multicenter_bonds().count(),
        span.noncovalent_bonds().count(),
        span.stereo_atoms().count(),
        span.stereo_bonds().count(),
    ]
}

fn molecule_remapping(images: &[Vec<usize>; 8]) -> MoleculeRemapping {
    fn remapping<Id: Copy + Debug + Into<usize> + From<usize>>(images: &[usize]) -> Remapping<Id> {
        Remapping::new(images.iter().copied().map(Id::from).collect())
            .expect("selected entity order is a permutation")
    }
    MoleculeRemapping::new(
        GraphRemapping::new(
            remapping::<NodeId>(&images[0]),
            remapping::<EdgeId>(&images[1]),
        ),
        remapping::<DativeBondId>(&images[2]),
        remapping::<AromaticSystemId>(&images[3]),
        remapping::<MulticenterBondId>(&images[4]),
        remapping::<NoncovalentBondId>(&images[5]),
        remapping::<StereoAtomId>(&images[6]),
        remapping::<StereoBondId>(&images[7]),
    )
}

/// Failure to construct a canonical [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizeError {
    /// A carried value is contradictory or incompatible with its stored atom or ligand frame.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizeError {
    /// Normalization or span materialization failed, or literal electron counts
    /// do not match their atom list.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

impl From<ReactionSpanCanonicalizeError> for ReactionCanonicalizeError {
    fn from(error: ReactionSpanCanonicalizeError) -> Self {
        match error {
            ReactionSpanCanonicalizeError::Contradiction(contradiction) => {
                Self::Contradiction(contradiction)
            }
        }
    }
}

#[cfg(test)]
mod tests;
