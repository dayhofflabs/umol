//! Ring detection and enumeration.

use std::collections::{BTreeMap, HashSet, VecDeque};

use umol_graph_core::{Graph, RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm};

use super::id::{AtomId, BondId};
use super::view::RingView;

/// Index of a ring within a `RingSet`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RingId(pub u32);

impl RingId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Ring {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

impl Ring {
    fn new(atoms: Vec<AtomId>, bonds: Vec<BondId>) -> Option<Self> {
        if atoms.len() < 3 || atoms.len() != bonds.len() {
            return None;
        }
        Some(Self { atoms, bonds })
    }
}

pub(crate) fn intersection<T: Copy + Eq>(a: &[T], b: &[T]) -> Vec<T> {
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    small
        .iter()
        .copied()
        .filter(|x| large.contains(x))
        .collect()
}

/// Topological relation between two rings in a `RingSet`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingRelation {
    Identical,
    Disjoint,
    Spiro,
    Fused,
    Bridged,
    MultiSpiro,
    Noncontiguous,
}

/// Selection of which ring set to enumerate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingSetKind {
    /// Every elementary cycle, up to the configured maximum ring size.
    Simple,
    /// The union of all minimum cycle bases, up to the configured maximum ring
    /// size.
    Relevant,
}

/// Semantic definition of the rings returned for a molecule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingModel {
    /// The cycle-set semantics used to define rings.
    pub kind: RingSetKind,
    /// Largest number of bonds admitted in a ring.
    pub max_ring_size: usize,
}

impl Default for RingModel {
    fn default() -> Self {
        Self {
            kind: RingSetKind::Relevant,
            max_ring_size: 22,
        }
    }
}

/// Algorithms used to compute each supported [`RingSetKind`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingConfig {
    /// Algorithm used when [`RingModel::kind`] is [`RingSetKind::Simple`].
    pub simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm,
    /// Algorithm used when [`RingModel::kind`] is [`RingSetKind::Relevant`].
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}

impl Default for RingConfig {
    fn default() -> Self {
        Self {
            simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        }
    }
}

/// Collection of rings enumerated from a parent molecule, with atom/bond →
/// ring reverse indices and a `RingGraph` describing pairwise relations
/// between rings.
#[derive(Debug, Clone)]
pub struct RingSet {
    kind: RingSetKind,
    max_ring_size: usize,
    rings: Vec<Ring>,
    atom_to_rings: BTreeMap<AtomId, Vec<RingId>>,
    bond_to_rings: BTreeMap<BondId, Vec<RingId>>,
    ring_graph: RingGraph,
}

impl RingSet {
    /// Enumerate the rings selected by `model` using `config`.
    pub(super) fn enumerate(graph: &Graph, model: RingModel, config: RingConfig) -> Self {
        let raw_cycles = match model.kind {
            RingSetKind::Simple => {
                graph.enumerate_simple_cycles(model.max_ring_size, config.simple_cycle_algorithm)
            }
            RingSetKind::Relevant => graph
                .enumerate_relevant_cycles(model.max_ring_size, config.relevant_cycle_algorithm),
        };

        let all_rings: Vec<Ring> = raw_cycles
            .into_iter()
            .filter_map(|cycle| {
                let ring_atoms = cycle.nodes().iter().copied().map(AtomId::from).collect();
                let ring_bonds = cycle.edges().iter().copied().map(BondId::from).collect();
                Ring::new(ring_atoms, ring_bonds)
            })
            .collect();

        Self::from_parts(model.kind, model.max_ring_size, all_rings)
    }

    fn from_parts(kind: RingSetKind, max_ring_size: usize, rings: Vec<Ring>) -> Self {
        let mut atom_to_rings: BTreeMap<AtomId, Vec<RingId>> = BTreeMap::new();
        let mut bond_to_rings: BTreeMap<BondId, Vec<RingId>> = BTreeMap::new();
        for (id, ring) in rings.iter().enumerate() {
            let ring_idx = RingId(id as u32);
            for &atom in &ring.atoms {
                atom_to_rings.entry(atom).or_default().push(ring_idx);
            }
            for &bond in &ring.bonds {
                bond_to_rings.entry(bond).or_default().push(ring_idx);
            }
        }

        let ring_graph = RingGraph::new(&rings);

        Self {
            kind,
            max_ring_size,
            rings,
            atom_to_rings,
            bond_to_rings,
            ring_graph,
        }
    }

    pub fn kind(&self) -> RingSetKind {
        self.kind
    }

    pub fn max_ring_size(&self) -> usize {
        self.max_ring_size
    }

    pub fn count(&self) -> usize {
        self.rings.len()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = RingId> {
        (0..self.rings.len()).map(|i| RingId(i as u32))
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = RingView<'_>> + '_ {
        self.rings
            .iter()
            .enumerate()
            .map(|(i, r)| RingView::new(RingId(i as u32), &r.atoms, &r.bonds))
    }

    pub fn get(&self, id: RingId) -> Option<RingView<'_>> {
        let r = self.rings.get(id.index())?;
        Some(RingView::new(id, &r.atoms, &r.bonds))
    }

    pub fn shared_atoms(&self, a: RingId, b: RingId) -> Vec<AtomId> {
        let (Some(ra), Some(rb)) = (self.rings.get(a.index()), self.rings.get(b.index())) else {
            return Vec::new();
        };
        intersection(&ra.atoms, &rb.atoms)
    }

    pub fn shared_bonds(&self, a: RingId, b: RingId) -> Vec<BondId> {
        let (Some(ra), Some(rb)) = (self.rings.get(a.index()), self.rings.get(b.index())) else {
            return Vec::new();
        };
        intersection(&ra.bonds, &rb.bonds)
    }

    pub fn relation(&self, a: RingId, b: RingId) -> RingRelation {
        self.ring_graph.relation(a, b)
    }

    pub fn are_spiro(&self, a: RingId, b: RingId) -> bool {
        self.relation(a, b) == RingRelation::Spiro
    }

    pub fn are_fused(&self, a: RingId, b: RingId) -> bool {
        self.relation(a, b) == RingRelation::Fused
    }

    pub fn are_bridged(&self, a: RingId, b: RingId) -> bool {
        self.relation(a, b) == RingRelation::Bridged
    }

    pub fn spiro_neighbors(&self, i: RingId) -> Vec<RingId> {
        let mut result: Vec<RingId> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Spiro).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn fused_neighbors(&self, i: RingId) -> Vec<RingId> {
        let mut result: Vec<RingId> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Fused).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn bridged_neighbors(&self, i: RingId) -> Vec<RingId> {
        let mut result: Vec<RingId> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Bridged).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn fused_components(&self) -> Vec<Vec<RingId>> {
        let mut visited: HashSet<RingId> = HashSet::new();
        let mut components: Vec<Vec<RingId>> = Vec::new();

        for ring in self.ids() {
            if visited.contains(&ring) {
                continue;
            }
            let component = self.fused_component(ring);
            for &r in &component {
                visited.insert(r);
            }
            components.push(component);
        }

        components.sort_by_key(|component| component.first().copied().map(RingId::index));
        components
    }

    pub fn fused_component(&self, root: RingId) -> Vec<RingId> {
        let mut visited: HashSet<RingId> = HashSet::new();
        let mut queue: VecDeque<RingId> = VecDeque::new();
        queue.push_back(root);
        visited.insert(root);

        while let Some(current) = queue.pop_front() {
            for neighbor in self.fused_neighbors(current) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        let mut result: Vec<RingId> = visited.into_iter().collect();
        result.sort_unstable();
        result
    }

    pub fn contains_atom(&self, atom: AtomId) -> bool {
        self.atom_to_rings.contains_key(&atom)
    }

    pub fn atom_smallest_ring_size(&self, atom: AtomId) -> Option<usize> {
        self.atom_to_rings.get(&atom).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].atoms.len())
                .min()
        })
    }

    pub fn contains_bond(&self, bond: BondId) -> bool {
        self.bond_to_rings.contains_key(&bond)
    }

    pub fn bond_smallest_ring_size(&self, bond: BondId) -> Option<usize> {
        self.bond_to_rings.get(&bond).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].atoms.len())
                .min()
        })
    }

    pub fn graph(&self) -> RingGraph {
        self.ring_graph.clone()
    }
}

/// Edge in a `RingGraph`: the two rings it connects and their relation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RingConnection {
    pub first: RingId,
    pub second: RingId,
    pub relation: RingRelation,
}

/// Ring graph: `RingConnection`s plus adjacency lists. Connected
/// components of the fused/bridged/spiro subgraph correspond to ring
/// systems.
#[derive(Debug, Clone)]
pub struct RingGraph {
    edges: Vec<RingConnection>,
    neighbors: Vec<Vec<(RingId, RingRelation)>>,
}

impl RingGraph {
    fn new(rings: &[Ring]) -> Self {
        let mut edges = Vec::new();
        let mut neighbors = vec![Vec::new(); rings.len()];
        let ids: Vec<RingId> = (0..rings.len()).map(|i| RingId(i as u32)).collect();
        for (i, &a) in ids.iter().enumerate() {
            for &b in &ids[i + 1..] {
                let relation = classify_ring_relation(&rings[a.index()], &rings[b.index()]);
                if relation == RingRelation::Disjoint || relation == RingRelation::Identical {
                    continue;
                }
                edges.push(RingConnection {
                    first: a,
                    second: b,
                    relation,
                });
                neighbors[a.index()].push((b, relation));
                neighbors[b.index()].push((a, relation));
            }
        }
        edges.sort_by_key(|e| (e.first, e.second, e.relation as u8));
        for n in &mut neighbors {
            n.sort_by_key(|(id, rel)| (*id, *rel as u8));
        }
        Self { edges, neighbors }
    }

    pub fn edges(&self) -> &[RingConnection] {
        &self.edges
    }

    pub fn neighbors(&self, ring: RingId) -> Vec<(RingId, RingRelation)> {
        self.neighbors
            .get(ring.index())
            .cloned()
            .unwrap_or_default()
    }

    pub fn relation(&self, a: RingId, b: RingId) -> RingRelation {
        if a == b {
            return RingRelation::Identical;
        }
        self.neighbors
            .get(a.index())
            .and_then(|neighbors| {
                neighbors
                    .iter()
                    .find_map(|(id, rel)| (*id == b).then_some(*rel))
            })
            .unwrap_or(RingRelation::Disjoint)
    }
}

fn classify_ring_relation(a: &Ring, b: &Ring) -> RingRelation {
    let shared_bonds = intersection(&a.bonds, &b.bonds);
    if shared_bonds.is_empty() {
        return match intersection(&a.atoms, &b.atoms).len() {
            0 => RingRelation::Disjoint,
            1 => RingRelation::Spiro,
            _ => RingRelation::MultiSpiro,
        };
    }

    if shared_bonds.len() == 1 {
        return RingRelation::Fused;
    }

    let bonds_a = &a.bonds;
    let n = bonds_a.len();
    let mut runs = 0usize;
    for i in 0..n {
        let curr_shared = shared_bonds.contains(&bonds_a[i]);
        let prev_shared = shared_bonds.contains(&bonds_a[(i + n - 1) % n]);
        if curr_shared && !prev_shared {
            runs += 1;
        }
    }

    if runs <= 1 {
        RingRelation::Bridged
    } else {
        RingRelation::Noncontiguous
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[fixture]
    fn triangle_set() -> RingSet {
        RingSet::from_parts(
            RingSetKind::Simple,
            6,
            vec![Ring::new(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                vec![BondId(0), BondId(1), BondId(2)],
            )
            .unwrap()],
        )
    }

    #[fixture]
    fn fused_pair() -> RingSet {
        RingSet::from_parts(
            RingSetKind::Simple,
            10,
            vec![
                Ring::new(
                    vec![
                        AtomId(0),
                        AtomId(1),
                        AtomId(2),
                        AtomId(3),
                        AtomId(4),
                        AtomId(5),
                    ],
                    vec![
                        BondId(0),
                        BondId(1),
                        BondId(2),
                        BondId(3),
                        BondId(4),
                        BondId(5),
                    ],
                )
                .unwrap(),
                Ring::new(
                    vec![
                        AtomId(1),
                        AtomId(2),
                        AtomId(6),
                        AtomId(7),
                        AtomId(8),
                        AtomId(9),
                    ],
                    vec![
                        BondId(1),
                        BondId(6),
                        BondId(7),
                        BondId(8),
                        BondId(9),
                        BondId(10),
                    ],
                )
                .unwrap(),
            ],
        )
    }

    #[fixture]
    fn spiro_pair() -> RingSet {
        RingSet::from_parts(
            RingSetKind::Simple,
            6,
            vec![
                Ring::new(
                    vec![AtomId(0), AtomId(1), AtomId(2)],
                    vec![BondId(0), BondId(1), BondId(2)],
                )
                .unwrap(),
                Ring::new(
                    vec![AtomId(0), AtomId(3), AtomId(4)],
                    vec![BondId(3), BondId(4), BondId(5)],
                )
                .unwrap(),
            ],
        )
    }

    #[fixture]
    fn bridged_pair() -> RingSet {
        RingSet::from_parts(
            RingSetKind::Simple,
            6,
            vec![
                Ring::new(
                    vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
                    vec![BondId(0), BondId(1), BondId(2), BondId(3)],
                )
                .unwrap(),
                Ring::new(
                    vec![AtomId(0), AtomId(1), AtomId(2), AtomId(4)],
                    vec![BondId(0), BondId(1), BondId(4), BondId(5)],
                )
                .unwrap(),
            ],
        )
    }

    /// Atom 1 sits in both a triangle and a 6-ring.
    #[fixture]
    fn shared_atom_set() -> RingSet {
        RingSet::from_parts(
            RingSetKind::Simple,
            10,
            vec![
                Ring::new(
                    vec![AtomId(1), AtomId(2), AtomId(3)],
                    vec![BondId(0), BondId(1), BondId(2)],
                )
                .unwrap(),
                Ring::new(
                    vec![
                        AtomId(1),
                        AtomId(4),
                        AtomId(5),
                        AtomId(6),
                        AtomId(7),
                        AtomId(8),
                    ],
                    vec![
                        BondId(3),
                        BondId(4),
                        BondId(5),
                        BondId(6),
                        BondId(7),
                        BondId(8),
                    ],
                )
                .unwrap(),
            ],
        )
    }

    #[rstest]
    #[case::zero(RingId(0), 0)]
    #[case::seven(RingId(7), 7)]
    fn test_ring_id_index(#[case] id: RingId, #[case] expected: usize) {
        assert_eq!(id.index(), expected);
    }

    #[rstest]
    #[case::valid_triangle(
        vec![AtomId(0), AtomId(1), AtomId(2)],
        vec![BondId(0), BondId(1), BondId(2)],
        true,
    )]
    #[case::too_small(
        vec![AtomId(0), AtomId(1)],
        vec![BondId(0), BondId(1)],
        false,
    )]
    #[case::empty(vec![], vec![], false)]
    #[case::atom_bond_len_mismatch(
        vec![AtomId(0), AtomId(1), AtomId(2)],
        vec![BondId(0), BondId(1)],
        false,
    )]
    fn test_ring_new(#[case] atoms: Vec<AtomId>, #[case] bonds: Vec<BondId>, #[case] valid: bool) {
        assert_eq!(Ring::new(atoms, bonds).is_some(), valid);
    }

    #[rstest]
    #[case::both_empty(vec![], vec![], vec![])]
    #[case::first_empty(vec![], vec![1, 2], vec![])]
    #[case::disjoint(vec![1, 2, 3], vec![4, 5, 6], vec![])]
    #[case::full_overlap(vec![1, 2, 3], vec![1, 2, 3], vec![1, 2, 3])]
    #[case::partial(vec![1, 2, 3], vec![2, 3, 4], vec![2, 3])]
    #[case::first_shorter(vec![1, 2], vec![1, 2, 3, 4], vec![1, 2])]
    #[case::second_shorter(vec![1, 2, 3, 4], vec![2, 3], vec![2, 3])]
    fn test_intersection(#[case] a: Vec<u32>, #[case] b: Vec<u32>, #[case] expected: Vec<u32>) {
        assert_eq!(intersection(&a, &b), expected);
    }

    #[rstest]
    fn test_ring_model_default() {
        assert_eq!(
            RingModel::default(),
            RingModel {
                kind: RingSetKind::Relevant,
                max_ring_size: 22,
            }
        );
    }

    #[rstest]
    #[case::equal(
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: 8,
        },
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: 8,
        },
        true,
    )]
    #[case::different_kind(
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: 8,
        },
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size: 8,
        },
        false,
    )]
    #[case::different_max_ring_size(
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size: 8,
        },
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size: 9,
        },
        false,
    )]
    fn test_ring_model_eq(
        #[case] left: RingModel,
        #[case] right: RingModel,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }

    #[rstest]
    fn test_ring_config_default() {
        assert_eq!(
            RingConfig::default(),
            RingConfig {
                simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            }
        );
    }

    #[rstest]
    fn test_ring_view_atoms(triangle_set: RingSet) {
        let view = triangle_set.get(RingId(0)).unwrap();
        assert_eq!(view.atoms(), &[AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_ring_view_bonds(triangle_set: RingSet) {
        let view = triangle_set.get(RingId(0)).unwrap();
        assert_eq!(view.bonds(), &[BondId(0), BondId(1), BondId(2)]);
    }

    #[rstest]
    fn test_ring_view_len(triangle_set: RingSet) {
        assert_eq!(triangle_set.get(RingId(0)).unwrap().len(), 3);
    }

    #[rstest]
    fn test_ring_view_is_empty(triangle_set: RingSet) {
        assert!(!triangle_set.get(RingId(0)).unwrap().is_empty());
    }

    #[rstest]
    fn test_ring_view_shared_atoms(fused_pair: RingSet) {
        let va = fused_pair.get(RingId(0)).unwrap();
        let vb = fused_pair.get(RingId(1)).unwrap();
        assert_eq!(va.shared_atoms(&vb), vec![AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_ring_view_shared_bonds(fused_pair: RingSet) {
        let va = fused_pair.get(RingId(0)).unwrap();
        let vb = fused_pair.get(RingId(1)).unwrap();
        assert_eq!(va.shared_bonds(&vb), vec![BondId(1)]);
    }

    #[rstest]
    #[case::simple_hexagon(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: 6,
        },
        vec![(
            vec![
                AtomId(0),
                AtomId(1),
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5),
            ],
            vec![
                BondId(0),
                BondId(1),
                BondId(2),
                BondId(3),
                BondId(4),
                BondId(5),
            ],
        )],
    )]
    #[case::simple_hexagon_cutoff(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: 5,
        },
        vec![],
    )]
    #[case::relevant_k4(
        Graph::new(
            4,
            &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
        ),
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size: 4,
        },
        vec![
            (
                vec![AtomId(0), AtomId(1), AtomId(2)],
                vec![BondId(0), BondId(3), BondId(1)],
            ),
            (
                vec![AtomId(0), AtomId(1), AtomId(3)],
                vec![BondId(0), BondId(4), BondId(2)],
            ),
            (
                vec![AtomId(0), AtomId(2), AtomId(3)],
                vec![BondId(1), BondId(5), BondId(2)],
            ),
            (
                vec![AtomId(1), AtomId(2), AtomId(3)],
                vec![BondId(3), BondId(5), BondId(4)],
            ),
        ],
    )]
    fn test_ring_set_enumerate(
        #[case] graph: Graph,
        #[case] model: RingModel,
        #[case] expected: Vec<(Vec<AtomId>, Vec<BondId>)>,
    ) {
        let actual = RingSet::enumerate(&graph, model, RingConfig::default())
            .iter()
            .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::relevant(RingSetKind::Relevant, 5)]
    #[case::simple(RingSetKind::Simple, 5)]
    fn test_ring_set_from_parts(#[case] kind: RingSetKind, #[case] max_ring_size: usize) {
        let set = RingSet::from_parts(kind, max_ring_size, vec![]);
        assert_eq!(set.kind(), kind);
        assert_eq!(set.max_ring_size(), max_ring_size);
        assert_eq!(set.count(), 0);
        assert_eq!(set.ids().collect::<Vec<_>>(), Vec::<RingId>::new());
        assert_eq!(set.iter().count(), 0);
        assert!(set.get(RingId(0)).is_none());
        assert!(!set.contains_atom(AtomId(0)));
        assert!(!set.contains_bond(BondId(0)));
        assert_eq!(set.atom_smallest_ring_size(AtomId(0)), None);
        assert_eq!(set.bond_smallest_ring_size(BondId(0)), None);
        assert_eq!(set.shared_atoms(RingId(0), RingId(1)), Vec::<AtomId>::new());
        assert_eq!(set.shared_bonds(RingId(0), RingId(1)), Vec::<BondId>::new());
        assert_eq!(set.graph().edges(), &[]);
    }

    #[rstest]
    fn test_ring_set_kind(triangle_set: RingSet) {
        assert_eq!(triangle_set.kind(), RingSetKind::Simple);
    }

    #[rstest]
    fn test_ring_set_max_ring_size(triangle_set: RingSet) {
        assert_eq!(triangle_set.max_ring_size(), 6);
    }

    #[rstest]
    fn test_ring_set_count(triangle_set: RingSet) {
        assert_eq!(triangle_set.count(), 1);
    }

    #[rstest]
    fn test_ring_set_ids(triangle_set: RingSet) {
        let empty = RingSet::from_parts(RingSetKind::Simple, 6, vec![]);
        let mut empty_ids = empty.ids();
        assert_eq!(empty_ids.len(), 0);
        assert_eq!(empty_ids.size_hint(), (0, Some(0)));
        assert_eq!(empty_ids.next(), None);

        let mut ids = triangle_set.ids();
        assert_eq!(ids.len(), 1);
        assert_eq!(ids.size_hint(), (1, Some(1)));
        assert_eq!(ids.next(), Some(RingId(0)));
        assert_eq!(ids.len(), 0);
        assert_eq!(ids.size_hint(), (0, Some(0)));
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_ring_set_iter(triangle_set: RingSet) {
        let empty = RingSet::from_parts(RingSetKind::Simple, 6, vec![]);
        let mut empty_iter = empty.iter();
        assert_eq!(empty_iter.len(), 0);
        assert_eq!(empty_iter.size_hint(), (0, Some(0)));
        assert_eq!(empty_iter.next().map(|view| view.id), None);

        let mut iter = triangle_set.iter();
        assert_eq!(iter.len(), 1);
        assert_eq!(iter.size_hint(), (1, Some(1)));
        assert_eq!(iter.next().map(|view| view.id), Some(RingId(0)));
        assert_eq!(iter.len(), 0);
        assert_eq!(iter.size_hint(), (0, Some(0)));
        assert_eq!(iter.next().map(|view| view.id), None);
    }

    #[rstest]
    #[case::in_range(RingId(0), Some(RingId(0)))]
    #[case::out_of_range(RingId(99), None)]
    fn test_ring_set_get(
        triangle_set: RingSet,
        #[case] id: RingId,
        #[case] expected: Option<RingId>,
    ) {
        assert_eq!(triangle_set.get(id).map(|v| v.id), expected);
    }

    #[rstest]
    #[case::ring_to_ring(RingId(0), RingId(1), vec![AtomId(1), AtomId(2)])]
    #[case::oob_second(RingId(0), RingId(99), vec![])]
    #[case::oob_first(RingId(99), RingId(0), vec![])]
    fn test_ring_set_shared_atoms(
        fused_pair: RingSet,
        #[case] a: RingId,
        #[case] b: RingId,
        #[case] expected: Vec<AtomId>,
    ) {
        assert_eq!(fused_pair.shared_atoms(a, b), expected);
    }

    #[rstest]
    #[case::ring_to_ring(RingId(0), RingId(1), vec![BondId(1)])]
    #[case::oob_second(RingId(0), RingId(99), vec![])]
    #[case::oob_first(RingId(99), RingId(0), vec![])]
    fn test_ring_set_shared_bonds(
        fused_pair: RingSet,
        #[case] a: RingId,
        #[case] b: RingId,
        #[case] expected: Vec<BondId>,
    ) {
        assert_eq!(fused_pair.shared_bonds(a, b), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), RingRelation::Fused)]
    #[case::spiro(spiro_pair(), RingRelation::Spiro)]
    #[case::bridged(bridged_pair(), RingRelation::Bridged)]
    fn test_ring_set_relation(#[case] set: RingSet, #[case] expected: RingRelation) {
        assert_eq!(set.relation(RingId(0), RingId(1)), expected);
        assert_eq!(set.relation(RingId(0), RingId(0)), RingRelation::Identical);
    }

    #[rstest]
    #[case::fused(fused_pair(), true)]
    #[case::spiro(spiro_pair(), false)]
    #[case::bridged(bridged_pair(), false)]
    fn test_ring_set_are_fused(#[case] set: RingSet, #[case] expected: bool) {
        assert_eq!(set.are_fused(RingId(0), RingId(1)), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), false)]
    #[case::spiro(spiro_pair(), true)]
    #[case::bridged(bridged_pair(), false)]
    fn test_ring_set_are_spiro(#[case] set: RingSet, #[case] expected: bool) {
        assert_eq!(set.are_spiro(RingId(0), RingId(1)), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), false)]
    #[case::spiro(spiro_pair(), false)]
    #[case::bridged(bridged_pair(), true)]
    fn test_ring_set_are_bridged(#[case] set: RingSet, #[case] expected: bool) {
        assert_eq!(set.are_bridged(RingId(0), RingId(1)), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), vec![])]
    #[case::spiro(spiro_pair(), vec![RingId(1)])]
    #[case::bridged(bridged_pair(), vec![])]
    fn test_ring_set_spiro_neighbors(#[case] set: RingSet, #[case] expected: Vec<RingId>) {
        assert_eq!(set.spiro_neighbors(RingId(0)), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), vec![RingId(1)])]
    #[case::spiro(spiro_pair(), vec![])]
    #[case::bridged(bridged_pair(), vec![])]
    fn test_ring_set_fused_neighbors(#[case] set: RingSet, #[case] expected: Vec<RingId>) {
        assert_eq!(set.fused_neighbors(RingId(0)), expected);
    }

    #[rstest]
    #[case::fused(fused_pair(), vec![])]
    #[case::spiro(spiro_pair(), vec![])]
    #[case::bridged(bridged_pair(), vec![RingId(1)])]
    fn test_ring_set_bridged_neighbors(#[case] set: RingSet, #[case] expected: Vec<RingId>) {
        assert_eq!(set.bridged_neighbors(RingId(0)), expected);
    }

    #[rstest]
    fn test_ring_set_fused_component(triangle_set: RingSet) {
        assert_eq!(triangle_set.fused_component(RingId(0)), vec![RingId(0)]);
    }

    #[rstest]
    fn test_ring_set_fused_components() {
        let set = RingSet::from_parts(
            RingSetKind::Simple,
            10,
            vec![
                Ring::new(
                    vec![
                        AtomId(0),
                        AtomId(1),
                        AtomId(2),
                        AtomId(3),
                        AtomId(4),
                        AtomId(5),
                    ],
                    vec![
                        BondId(0),
                        BondId(1),
                        BondId(2),
                        BondId(3),
                        BondId(4),
                        BondId(5),
                    ],
                )
                .unwrap(),
                Ring::new(
                    vec![
                        AtomId(1),
                        AtomId(2),
                        AtomId(6),
                        AtomId(7),
                        AtomId(8),
                        AtomId(9),
                    ],
                    vec![
                        BondId(1),
                        BondId(6),
                        BondId(7),
                        BondId(8),
                        BondId(9),
                        BondId(10),
                    ],
                )
                .unwrap(),
                Ring::new(
                    vec![AtomId(20), AtomId(21), AtomId(22)],
                    vec![BondId(20), BondId(21), BondId(22)],
                )
                .unwrap(),
            ],
        );
        assert_eq!(
            set.fused_components(),
            vec![vec![RingId(0), RingId(1)], vec![RingId(2)]],
        );
    }

    #[rstest]
    #[case::present(AtomId(0), true)]
    #[case::absent(AtomId(99), false)]
    fn test_ring_set_contains_atom(
        triangle_set: RingSet,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(triangle_set.contains_atom(atom), expected);
    }

    #[rstest]
    #[case::triangle_only(triangle_set(), AtomId(0), Some(3))]
    #[case::triangle_absent(triangle_set(), AtomId(99), None)]
    #[case::shared_picks_min(shared_atom_set(), AtomId(1), Some(3))]
    #[case::six_only(shared_atom_set(), AtomId(4), Some(6))]
    fn test_ring_set_atom_smallest_ring_size(
        #[case] set: RingSet,
        #[case] atom: AtomId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(set.atom_smallest_ring_size(atom), expected);
    }

    #[rstest]
    #[case::present(BondId(1), true)]
    #[case::absent(BondId(99), false)]
    fn test_ring_set_contains_bond(
        triangle_set: RingSet,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(triangle_set.contains_bond(bond), expected);
    }

    #[rstest]
    #[case::in_triangle(BondId(0), Some(3))]
    #[case::absent(BondId(99), None)]
    fn test_ring_set_bond_smallest_ring_size(
        triangle_set: RingSet,
        #[case] bond: BondId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(triangle_set.bond_smallest_ring_size(bond), expected);
    }

    #[rstest]
    fn test_ring_graph_new(fused_pair: RingSet) {
        let graph = fused_pair.graph();
        assert_eq!(
            graph.edges(),
            &[RingConnection {
                first: RingId(0),
                second: RingId(1),
                relation: RingRelation::Fused,
            }],
        );
    }

    #[rstest]
    #[case::ring_0(RingId(0), vec![(RingId(1), RingRelation::Fused)])]
    #[case::ring_1(RingId(1), vec![(RingId(0), RingRelation::Fused)])]
    #[case::oob(RingId(99), vec![])]
    fn test_ring_graph_neighbors(
        fused_pair: RingSet,
        #[case] ring: RingId,
        #[case] expected: Vec<(RingId, RingRelation)>,
    ) {
        assert_eq!(fused_pair.graph().neighbors(ring), expected);
    }

    #[rstest]
    #[case::self_is_identical(RingId(0), RingId(0), RingRelation::Identical)]
    #[case::oob_first(RingId(99), RingId(0), RingRelation::Disjoint)]
    #[case::oob_second(RingId(0), RingId(99), RingRelation::Disjoint)]
    fn test_ring_graph_relation(
        triangle_set: RingSet,
        #[case] a: RingId,
        #[case] b: RingId,
        #[case] expected: RingRelation,
    ) {
        assert_eq!(triangle_set.graph().relation(a, b), expected);
    }

    #[rstest]
    #[case::disjoint(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1), BondId(2)]).unwrap(),
        Ring::new(vec![AtomId(10), AtomId(11), AtomId(12)], vec![BondId(10), BondId(11), BondId(12)]).unwrap(),
        RingRelation::Disjoint,
    )]
    #[case::spiro(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1), BondId(2)]).unwrap(),
        Ring::new(vec![AtomId(0), AtomId(3), AtomId(4)], vec![BondId(3), BondId(4), BondId(5)]).unwrap(),
        RingRelation::Spiro,
    )]
    #[case::multispiro(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1), BondId(2)]).unwrap(),
        Ring::new(vec![AtomId(0), AtomId(2), AtomId(3)], vec![BondId(3), BondId(4), BondId(5)]).unwrap(),
        RingRelation::MultiSpiro,
    )]
    #[case::fused(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1), BondId(2)]).unwrap(),
        Ring::new(vec![AtomId(1), AtomId(2), AtomId(3)], vec![BondId(1), BondId(3), BondId(4)]).unwrap(),
        RingRelation::Fused,
    )]
    #[case::bridged(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)], vec![BondId(0), BondId(1), BondId(2), BondId(3)]).unwrap(),
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(4)], vec![BondId(0), BondId(1), BondId(4), BondId(5)]).unwrap(),
        RingRelation::Bridged,
    )]
    #[case::noncontiguous(
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)], vec![BondId(0), BondId(1), BondId(2), BondId(3), BondId(4), BondId(5)]).unwrap(),
        Ring::new(vec![AtomId(0), AtomId(1), AtomId(9), AtomId(3), AtomId(4), AtomId(8)], vec![BondId(0), BondId(10), BondId(11), BondId(3), BondId(12), BondId(13)]).unwrap(),
        RingRelation::Noncontiguous,
    )]
    fn test_classify_ring_relation(
        #[case] a: Ring,
        #[case] b: Ring,
        #[case] expected: RingRelation,
    ) {
        assert_eq!(classify_ring_relation(&a, &b), expected);
    }
}
