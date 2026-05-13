//! Ring detection and ring-system analysis.
//!
//! Ring enumeration uses Vismara's relevant-cycle algorithm from
//! `umol_graph_core`, decomposed over biconnected components.

use std::collections::{BTreeMap, HashSet, VecDeque};

use umol_graph_core::{BiconnectedComponentsAlgorithm, CycleEnumerationAlgorithm, Graph, NodeId};

use super::idx::{AtomId, BondId};

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

fn intersection<T: Copy + Eq>(a: &[T], b: &[T]) -> Vec<T> {
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    small
        .iter()
        .copied()
        .filter(|x| large.contains(x))
        .collect()
}

/// Borrowed view of a ring: its index plus atom and bond membership.
#[derive(Debug, Clone, Copy)]
pub struct RingView<'a> {
    pub id: RingId,
    atoms: &'a [AtomId],
    bonds: &'a [BondId],
}

impl<'a> RingView<'a> {
    pub fn atoms(&self) -> &'a [AtomId] {
        self.atoms
    }

    pub fn bonds(&self) -> &'a [BondId] {
        self.bonds
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn shared_atoms(&self, other: &RingView<'_>) -> Vec<AtomId> {
        intersection(self.atoms, other.atoms)
    }

    pub fn shared_bonds(&self, other: &RingView<'_>) -> Vec<BondId> {
        intersection(self.bonds, other.bonds)
    }
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

/// Selection of which cycle family to enumerate: `Simple` for the minimum
/// cycle basis, `Relevant` for the Vismara relevant-cycle set (the union of
/// all minimum cycle bases — unique and ordering-independent).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingFamily {
    Simple,
    Relevant,
}

/// Collection of rings for a `MoleculeAst`, with atom/bond → ring reverse
/// indices and a `RingGraph` describing pairwise relations between rings.
#[derive(Debug, Clone)]
pub struct RingSet {
    family: RingFamily,
    max_ring_size: usize,
    rings: Vec<Ring>,
    atom_to_rings: BTreeMap<AtomId, Vec<RingId>>,
    bond_to_rings: BTreeMap<BondId, Vec<RingId>>,
    ring_graph: RingGraph,
}

impl RingSet {
    pub fn empty() -> Self {
        Self {
            family: RingFamily::Simple,
            max_ring_size: 0,
            rings: Vec::new(),
            atom_to_rings: BTreeMap::new(),
            bond_to_rings: BTreeMap::new(),
            ring_graph: RingGraph {
                edges: Vec::new(),
                neighbors: Vec::new(),
            },
        }
    }

    fn from_rings(family: RingFamily, max_ring_size: usize, rings: Vec<Ring>) -> Self {
        if rings.is_empty() {
            let mut empty = Self::empty();
            empty.family = family;
            empty.max_ring_size = max_ring_size;
            return empty;
        }

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

        let ring_graph = RingGraph::from_ring_list(&rings);

        Self {
            family,
            max_ring_size,
            rings,
            atom_to_rings,
            bond_to_rings,
            ring_graph,
        }
    }

    pub fn family(&self) -> RingFamily {
        self.family
    }

    pub fn max_ring_size(&self) -> usize {
        self.max_ring_size
    }

    pub fn count(&self) -> usize {
        self.rings.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = RingId> {
        (0..self.rings.len()).map(|i| RingId(i as u32))
    }

    pub fn iter(&self) -> impl Iterator<Item = RingView<'_>> {
        self.rings.iter().enumerate().map(|(i, r)| RingView {
            id: RingId(i as u32),
            atoms: &r.atoms,
            bonds: &r.bonds,
        })
    }

    pub fn get(&self, id: RingId) -> Option<RingView<'_>> {
        let r = self.rings.get(id.index())?;
        Some(RingView {
            id,
            atoms: &r.atoms,
            bonds: &r.bonds,
        })
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
pub struct RingGraphEdge {
    pub source: RingId,
    pub target: RingId,
    pub relation: RingRelation,
}

/// Graph over rings, with `RingRelation`-labeled edges. Connected
/// components of the fused/bridged/spiro subgraph correspond to ring
/// systems.
#[derive(Debug, Clone)]
pub struct RingGraph {
    edges: Vec<RingGraphEdge>,
    neighbors: Vec<Vec<(RingId, RingRelation)>>,
}

impl RingGraph {
    fn from_ring_list(rings: &[Ring]) -> Self {
        let mut edges = Vec::new();
        let mut neighbors = vec![Vec::new(); rings.len()];
        let indices: Vec<RingId> = (0..rings.len()).map(|i| RingId(i as u32)).collect();
        for (i, &a) in indices.iter().enumerate() {
            for &b in &indices[i + 1..] {
                let relation = classify_ring_relation(&rings[a.index()], &rings[b.index()]);
                if relation == RingRelation::Disjoint || relation == RingRelation::Identical {
                    continue;
                }
                edges.push(RingGraphEdge {
                    source: a,
                    target: b,
                    relation,
                });
                neighbors[a.index()].push((b, relation));
                neighbors[b.index()].push((a, relation));
            }
        }
        edges.sort_by_key(|e| (e.source, e.target, e.relation as u8));
        for n in &mut neighbors {
            n.sort_by_key(|(id, rel)| (*id, *rel as u8));
        }
        Self { edges, neighbors }
    }

    pub fn edges(&self) -> &[RingGraphEdge] {
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

fn is_induced_cycle(graph: &Graph, cycle: &[NodeId]) -> bool {
    let n = cycle.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            if graph.find_edge(cycle[i], cycle[j]).is_some() {
                return false;
            }
        }
    }
    true
}

pub(crate) fn enumerate_rings(
    graph: &Graph,
    family: RingFamily,
    max_ring_size: usize,
    atom_filter: impl Fn(AtomId) -> bool,
) -> RingSet {
    let filtered_nodes: Vec<NodeId> = graph
        .node_ids()
        .filter(|&n| atom_filter(AtomId::from(n)))
        .collect();

    let use_subgraph = filtered_nodes.len() < graph.node_count();

    let (sub, node_map) = if use_subgraph {
        let sub = graph.induced_subgraph(&filtered_nodes);
        (sub.graph, sub.node_map)
    } else {
        let node_map: Vec<NodeId> = graph.node_ids().collect();
        (graph.clone(), node_map)
    };

    let bcc = sub.biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);

    let mut all_rings: Vec<Ring> = Vec::new();
    for component in &bcc {
        let comp_sub = sub.induced_subgraph(component);
        let raw_cycles = comp_sub
            .graph
            .enumerate_cycles(max_ring_size, CycleEnumerationAlgorithm::Vismara);

        let component_rings: Vec<Ring> = raw_cycles
            .into_iter()
            .filter(|cycle| match family {
                RingFamily::Relevant => is_induced_cycle(&comp_sub.graph, cycle),
                RingFamily::Simple => true,
            })
            .filter_map(|cycle| {
                let ring_atoms: Vec<AtomId> = cycle
                    .iter()
                    .map(|&local| {
                        let sub_node = comp_sub.node_map[local.index()];
                        let orig_node = node_map[sub_node.index()];
                        AtomId::from(orig_node)
                    })
                    .collect();
                let n = ring_atoms.len();
                let mut ring_bonds = Vec::with_capacity(n);
                for i in 0..n {
                    let a_orig = NodeId::from(ring_atoms[i]);
                    let b_orig = NodeId::from(ring_atoms[(i + 1) % n]);
                    let edge = graph.find_edge(a_orig, b_orig)?;
                    ring_bonds.push(BondId::from(edge));
                }
                Ring::new(ring_atoms, ring_bonds)
            })
            .collect();

        all_rings.extend(component_rings);
    }

    RingSet::from_rings(family, max_ring_size, all_rings)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn ring_of(atoms: &[u32], bonds: &[u32]) -> Ring {
        Ring::new(
            atoms.iter().copied().map(AtomId).collect(),
            bonds.iter().copied().map(BondId).collect(),
        )
        .expect("valid ring")
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
    fn test_ring_new(
        #[case] atoms: Vec<AtomId>,
        #[case] bonds: Vec<BondId>,
        #[case] valid: bool,
    ) {
        assert_eq!(Ring::new(atoms, bonds).is_some(), valid);
    }

    #[rstest]
    #[case(RingId(0), 0)]
    #[case(RingId(7), 7)]
    fn test_ring_idx_index(#[case] id: RingId, #[case] expected: usize) {
        assert_eq!(id.index(), expected);
    }

    #[rstest]
    #[case::disjoint(
        ring_of(&[0, 1, 2], &[0, 1, 2]),
        ring_of(&[10, 11, 12], &[10, 11, 12]),
        RingRelation::Disjoint,
    )]
    #[case::spiro(
        ring_of(&[0, 1, 2], &[0, 1, 2]),
        ring_of(&[0, 3, 4], &[3, 4, 5]),
        RingRelation::Spiro,
    )]
    #[case::multispiro(
        ring_of(&[0, 1, 2], &[0, 1, 2]),
        ring_of(&[0, 2, 3], &[3, 4, 5]),
        RingRelation::MultiSpiro,
    )]
    #[case::fused(
        ring_of(&[0, 1, 2], &[0, 1, 2]),
        ring_of(&[1, 2, 3], &[1, 3, 4]),
        RingRelation::Fused,
    )]
    #[case::bridged(
        ring_of(&[0, 1, 2, 3], &[0, 1, 2, 3]),
        ring_of(&[0, 1, 2, 4], &[0, 1, 4, 5]),
        RingRelation::Bridged,
    )]
    #[case::noncontiguous(
        ring_of(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4, 5]),
        ring_of(&[0, 1, 9, 3, 4, 8], &[0, 10, 11, 3, 12, 13]),
        RingRelation::Noncontiguous,
    )]
    fn test_classify_ring_relation(
        #[case] a: Ring,
        #[case] b: Ring,
        #[case] expected: RingRelation,
    ) {
        assert_eq!(classify_ring_relation(&a, &b), expected);
    }

    #[rstest]
    #[case::too_short(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![NodeId(0), NodeId(1)],
        false,
    )]
    #[case::induced_triangle(
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        vec![NodeId(0), NodeId(1), NodeId(2)],
        true,
    )]
    #[case::square_with_chord(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0], [0, 2]]),
        vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
        false,
    )]
    fn test_is_induced_cycle(
        #[case] graph: Graph,
        #[case] cycle: Vec<NodeId>,
        #[case] expected: bool,
    ) {
        assert_eq!(is_induced_cycle(&graph, &cycle), expected);
    }

    #[fixture]
    fn triangle_set() -> RingSet {
        RingSet::from_rings(RingFamily::Simple, 6, vec![ring_of(&[0, 1, 2], &[0, 1, 2])])
    }

    #[fixture]
    fn fused_pair() -> RingSet {
        let r1 = ring_of(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4, 5]);
        let r2 = ring_of(&[1, 2, 6, 7, 8, 9], &[1, 6, 7, 8, 9, 10]);
        RingSet::from_rings(RingFamily::Simple, 10, vec![r1, r2])
    }

    #[fixture]
    fn spiro_pair() -> RingSet {
        let r1 = ring_of(&[0, 1, 2], &[0, 1, 2]);
        let r2 = ring_of(&[0, 3, 4], &[3, 4, 5]);
        RingSet::from_rings(RingFamily::Simple, 6, vec![r1, r2])
    }

    #[fixture]
    fn bridged_pair() -> RingSet {
        let r1 = ring_of(&[0, 1, 2, 3], &[0, 1, 2, 3]);
        let r2 = ring_of(&[0, 1, 2, 4], &[0, 1, 4, 5]);
        RingSet::from_rings(RingFamily::Simple, 6, vec![r1, r2])
    }

    #[rstest]
    fn test_ring_set_empty() {
        let set = RingSet::empty();
        assert_eq!(set.family(), RingFamily::Simple);
        assert_eq!(set.max_ring_size(), 0);
        assert_eq!(set.count(), 0);
        assert_eq!(set.ids().collect::<Vec<_>>(), Vec::<RingId>::new());
        assert_eq!(set.iter().count(), 0);
        assert!(set.get(RingId(0)).is_none());
        assert!(!set.contains_atom(AtomId(0)));
        assert!(!set.contains_bond(BondId(0)));
        assert_eq!(set.atom_smallest_ring_size(AtomId(0)), None);
        assert_eq!(set.bond_smallest_ring_size(BondId(0)), None);
        assert_eq!(
            set.shared_atoms(RingId(0), RingId(1)),
            Vec::<AtomId>::new()
        );
        assert_eq!(
            set.shared_bonds(RingId(0), RingId(1)),
            Vec::<BondId>::new()
        );
        assert_eq!(set.graph().edges(), &[]);
    }

    #[rstest]
    #[case(RingFamily::Relevant)]
    #[case(RingFamily::Simple)]
    fn test_ring_set_from_rings_empty_preserves_family(#[case] family: RingFamily) {
        let set = RingSet::from_rings(family, 5, vec![]);
        assert_eq!(set.family(), family);
        assert_eq!(set.max_ring_size(), 5);
        assert_eq!(set.count(), 0);
    }

    #[rstest]
    fn test_ring_set_from_rings_accessors(triangle_set: RingSet) {
        assert_eq!(triangle_set.family(), RingFamily::Simple);
        assert_eq!(triangle_set.max_ring_size(), 6);
        assert_eq!(triangle_set.count(), 1);
        assert_eq!(triangle_set.ids().collect::<Vec<_>>(), vec![RingId(0)]);
        let views: Vec<RingId> = triangle_set.iter().map(|v| v.id).collect();
        assert_eq!(views, vec![RingId(0)]);
    }

    #[rstest]
    fn test_ring_set_get_out_of_range(triangle_set: RingSet) {
        assert!(triangle_set.get(RingId(99)).is_none());
    }

    #[rstest]
    fn test_ring_view_accessors(triangle_set: RingSet) {
        let view = triangle_set.get(RingId(0)).unwrap();
        assert_eq!(view.atoms(), &[AtomId(0), AtomId(1), AtomId(2)]);
        assert_eq!(view.bonds(), &[BondId(0), BondId(1), BondId(2)]);
        assert_eq!(view.len(), 3);
        assert!(!view.is_empty());
    }

    #[rstest]
    fn test_ring_view_shared_atoms_and_bonds(fused_pair: RingSet) {
        let va = fused_pair.get(RingId(0)).unwrap();
        let vb = fused_pair.get(RingId(1)).unwrap();
        assert_eq!(va.shared_atoms(&vb), vec![AtomId(1), AtomId(2)]);
        assert_eq!(va.shared_bonds(&vb), vec![BondId(1)]);
    }

    #[rstest]
    fn test_ring_set_membership(triangle_set: RingSet) {
        assert!(triangle_set.contains_atom(AtomId(0)));
        assert!(!triangle_set.contains_atom(AtomId(99)));
        assert!(triangle_set.contains_bond(BondId(1)));
        assert!(!triangle_set.contains_bond(BondId(99)));
        assert_eq!(triangle_set.atom_smallest_ring_size(AtomId(0)), Some(3));
        assert_eq!(triangle_set.atom_smallest_ring_size(AtomId(99)), None);
        assert_eq!(triangle_set.bond_smallest_ring_size(BondId(0)), Some(3));
        assert_eq!(triangle_set.bond_smallest_ring_size(BondId(99)), None);
    }

    #[rstest]
    fn test_ring_set_smallest_picks_minimum_of_multiple() {
        // Atom 1 is in both a triangle and a 6-ring.
        let r_small = ring_of(&[1, 2, 3], &[0, 1, 2]);
        let r_large = ring_of(&[1, 4, 5, 6, 7, 8], &[3, 4, 5, 6, 7, 8]);
        let set = RingSet::from_rings(RingFamily::Simple, 10, vec![r_small, r_large]);
        assert_eq!(set.atom_smallest_ring_size(AtomId(1)), Some(3));
        assert_eq!(set.atom_smallest_ring_size(AtomId(4)), Some(6));
    }

    #[rstest]
    fn test_ring_set_shared_oob_returns_empty(triangle_set: RingSet) {
        assert_eq!(
            triangle_set.shared_atoms(RingId(0), RingId(99)),
            Vec::<AtomId>::new(),
        );
        assert_eq!(
            triangle_set.shared_bonds(RingId(0), RingId(99)),
            Vec::<BondId>::new(),
        );
        assert_eq!(
            triangle_set.shared_atoms(RingId(99), RingId(0)),
            Vec::<AtomId>::new(),
        );
    }

    #[rstest]
    fn test_ring_set_shared_hits(fused_pair: RingSet) {
        assert_eq!(
            fused_pair.shared_atoms(RingId(0), RingId(1)),
            vec![AtomId(1), AtomId(2)],
        );
        assert_eq!(
            fused_pair.shared_bonds(RingId(0), RingId(1)),
            vec![BondId(1)],
        );
    }

    #[rstest]
    #[case::fused(fused_pair(), RingRelation::Fused)]
    #[case::spiro(spiro_pair(), RingRelation::Spiro)]
    #[case::bridged(bridged_pair(), RingRelation::Bridged)]
    fn test_ring_set_relation_by_kind(#[case] set: RingSet, #[case] expected: RingRelation) {
        assert_eq!(set.relation(RingId(0), RingId(1)), expected);
        assert_eq!(
            set.relation(RingId(0), RingId(0)),
            RingRelation::Identical
        );
        assert_eq!(
            set.are_fused(RingId(0), RingId(1)),
            expected == RingRelation::Fused,
        );
        assert_eq!(
            set.are_spiro(RingId(0), RingId(1)),
            expected == RingRelation::Spiro,
        );
        assert_eq!(
            set.are_bridged(RingId(0), RingId(1)),
            expected == RingRelation::Bridged,
        );
    }

    #[rstest]
    #[case::fused(fused_pair(), RingRelation::Fused)]
    #[case::spiro(spiro_pair(), RingRelation::Spiro)]
    #[case::bridged(bridged_pair(), RingRelation::Bridged)]
    fn test_ring_set_neighbors_by_kind(#[case] set: RingSet, #[case] expected_kind: RingRelation) {
        let hit = vec![RingId(1)];
        let miss = Vec::<RingId>::new();
        assert_eq!(
            set.fused_neighbors(RingId(0)),
            if expected_kind == RingRelation::Fused {
                hit.clone()
            } else {
                miss.clone()
            },
        );
        assert_eq!(
            set.spiro_neighbors(RingId(0)),
            if expected_kind == RingRelation::Spiro {
                hit.clone()
            } else {
                miss.clone()
            },
        );
        assert_eq!(
            set.bridged_neighbors(RingId(0)),
            if expected_kind == RingRelation::Bridged {
                hit
            } else {
                miss
            },
        );
    }

    #[rstest]
    fn test_ring_set_fused_component_single(triangle_set: RingSet) {
        assert_eq!(triangle_set.fused_component(RingId(0)), vec![RingId(0)]);
    }

    #[rstest]
    fn test_ring_set_fused_components_mixed() {
        let r1 = ring_of(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4, 5]);
        let r2 = ring_of(&[1, 2, 6, 7, 8, 9], &[1, 6, 7, 8, 9, 10]);
        let r3 = ring_of(&[20, 21, 22], &[20, 21, 22]);
        let set = RingSet::from_rings(RingFamily::Simple, 10, vec![r1, r2, r3]);
        assert_eq!(
            set.fused_components(),
            vec![vec![RingId(0), RingId(1)], vec![RingId(2)]],
        );
    }

    #[rstest]
    fn test_ring_graph_edges_and_neighbors(fused_pair: RingSet) {
        let graph = fused_pair.graph();
        assert_eq!(
            graph.edges(),
            &[RingGraphEdge {
                source: RingId(0),
                target: RingId(1),
                relation: RingRelation::Fused,
            }],
        );
        assert_eq!(
            graph.neighbors(RingId(0)),
            vec![(RingId(1), RingRelation::Fused)],
        );
        assert_eq!(
            graph.neighbors(RingId(1)),
            vec![(RingId(0), RingRelation::Fused)],
        );
        assert_eq!(
            graph.neighbors(RingId(99)),
            Vec::<(RingId, RingRelation)>::new(),
        );
    }

    #[rstest]
    #[case::self_is_identical(RingId(0), RingId(0), RingRelation::Identical)]
    #[case::oob_first(RingId(99), RingId(0), RingRelation::Disjoint)]
    #[case::oob_second(RingId(0), RingId(99), RingRelation::Disjoint)]
    fn test_ring_graph_relation_edges(
        triangle_set: RingSet,
        #[case] a: RingId,
        #[case] b: RingId,
        #[case] expected: RingRelation,
    ) {
        assert_eq!(triangle_set.graph().relation(a, b), expected);
    }

    #[rstest]
    #[case::full_hexagon(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        RingFamily::Simple,
        10,
        |_: AtomId| true,
        1,
    )]
    #[case::filter_breaks_cycle(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        RingFamily::Simple,
        10,
        |a: AtomId| a.0 != 5,
        0,
    )]
    fn test_enumerate_rings_count(
        #[case] graph: Graph,
        #[case] family: RingFamily,
        #[case] max_ring_size: usize,
        #[case] atom_filter: fn(AtomId) -> bool,
        #[case] expected_count: usize,
    ) {
        let set = enumerate_rings(&graph, family, max_ring_size, atom_filter);
        assert_eq!(set.count(), expected_count);
    }

    #[rstest]
    fn test_enumerate_rings_induced_keeps_only_chord_free_cycles() {
        // K4: 4 nodes fully connected. Simple enumeration yields the three
        // 4-cycles plus four 3-cycles; Induced keeps only the 3-cycles.
        let graph = Graph::new(4, &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
        let simple = enumerate_rings(&graph, RingFamily::Simple, 4, |_| true);
        let induced = enumerate_rings(&graph, RingFamily::Relevant, 4, |_| true);
        assert!(simple.count() >= induced.count());
        for view in induced.iter() {
            assert_eq!(view.len(), 3);
        }
    }
}
