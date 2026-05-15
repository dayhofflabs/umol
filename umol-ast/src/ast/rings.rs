//! Ring detection and ring-system analysis.
//!
//! Ring enumeration uses Vismara's relevant-cycle algorithm from
//! `umol_graph_core`, decomposed over biconnected components.

use std::collections::{BTreeMap, HashSet, VecDeque};

use umol_graph_core::{CycleEnumerationAlgorithm, Graph, NodeId};

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

/// Borrowed view of a ring: its index plus its atom and bond membership.
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

/// Collection of rings enumerated from a parent molecule, with atom/bond →
/// ring reverse indices and a `RingGraph` describing pairwise relations
/// between rings.
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
    /// Enumerate rings of `family` up to `max_ring_size` on the atoms of
    /// `graph` passing `atom_filter`. The actual cycle algorithm is
    /// Vismara on biconnected components; `RingFamily::Relevant` keeps only
    /// chordless cycles.
    pub(super) fn enumerate(
        family: RingFamily,
        max_ring_size: usize,
        atom_filter: impl Fn(AtomId) -> bool,
        graph: &Graph,
    ) -> Self {
        let filtered_nodes: Vec<NodeId> = graph
            .node_ids()
            .filter(|&n| atom_filter(AtomId::from(n)))
            .collect();

        let use_subgraph = filtered_nodes.len() < graph.node_count();

        let (sub, host_nodes) = if use_subgraph {
            let embedding = graph.induced_subgraph(&filtered_nodes);
            (embedding.extract(), embedding.host_nodes().to_vec())
        } else {
            let host_nodes: Vec<NodeId> = graph.node_ids().collect();
            (graph.clone(), host_nodes)
        };

        let raw_cycles =
            sub.enumerate_cycles(max_ring_size, CycleEnumerationAlgorithm::Vismara);

        let all_rings: Vec<Ring> = raw_cycles
            .into_iter()
            .filter(|cycle| match family {
                RingFamily::Relevant => is_induced_cycle(&sub, cycle),
                RingFamily::Simple => true,
            })
            .filter_map(|cycle| {
                let ring_atoms: Vec<AtomId> = cycle
                    .iter()
                    .map(|&sub_node| AtomId::from(host_nodes[sub_node.index()]))
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

        Self::from_parts(family, max_ring_size, all_rings)
    }

    fn from_parts(family: RingFamily, max_ring_size: usize, rings: Vec<Ring>) -> Self {
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

    pub fn iter(&self) -> impl Iterator<Item = RingView<'_>> + '_ {
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
    fn new(rings: &[Ring]) -> Self {
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeAst;

    #[fixture]
    fn triangle_set() -> RingSet {
        RingSet::from_parts(
            RingFamily::Simple,
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
            RingFamily::Simple,
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
            RingFamily::Simple,
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
            RingFamily::Simple,
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
            RingFamily::Simple,
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
    #[case::full_hexagon(
        6,
        &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]][..],
        RingFamily::Simple,
        10,
        |_: AtomId| true,
        1,
    )]
    #[case::filter_breaks_cycle(
        6,
        &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]][..],
        RingFamily::Simple,
        10,
        |a: AtomId| a.0 != 5,
        0,
    )]
    fn test_ring_set_enumerate(
        #[case] node_count: usize,
        #[case] edges: &[[u32; 2]],
        #[case] family: RingFamily,
        #[case] max_ring_size: usize,
        #[case] atom_filter: fn(AtomId) -> bool,
        #[case] expected_count: usize,
    ) {
        let atoms = vec![AtomAst::default(); node_count];
        let bonds: Vec<_> = edges
            .iter()
            .map(|[a, b]| (AtomId(*a), AtomId(*b), BondAst::default()))
            .collect();
        let mol = MoleculeAst::from_atoms_and_bonds(atoms, bonds);
        let set = mol.rings_with(family, max_ring_size, atom_filter);
        assert_eq!(set.count(), expected_count);
    }

    #[rstest]
    fn test_ring_set_enumerate_relevant() {
        // K4 (4 nodes fully connected): Simple includes 4-cycles spanning
        // chords; Relevant keeps only the chordless 3-cycles.
        let atoms = vec![AtomAst::default(); 4];
        let edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
        let bonds: Vec<_> = edges
            .iter()
            .map(|[a, b]| (AtomId(*a), AtomId(*b), BondAst::default()))
            .collect();
        let mol = MoleculeAst::from_atoms_and_bonds(atoms, bonds);
        let simple = mol.rings_with(RingFamily::Simple, 4, |_| true);
        let relevant = mol.rings_with(RingFamily::Relevant, 4, |_| true);
        assert!(simple.count() >= relevant.count());
        for view in relevant.iter() {
            assert_eq!(view.len(), 3);
        }
    }

    #[rstest]
    #[case::relevant(RingFamily::Relevant, 5)]
    #[case::simple(RingFamily::Simple, 5)]
    fn test_ring_set_from_parts(#[case] family: RingFamily, #[case] max_ring_size: usize) {
        let set = RingSet::from_parts(family, max_ring_size, vec![]);
        assert_eq!(set.family(), family);
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
    fn test_ring_set_family(triangle_set: RingSet) {
        assert_eq!(triangle_set.family(), RingFamily::Simple);
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
        assert_eq!(triangle_set.ids().collect::<Vec<_>>(), vec![RingId(0)]);
    }

    #[rstest]
    fn test_ring_set_iter(triangle_set: RingSet) {
        let views: Vec<RingId> = triangle_set.iter().map(|v| v.id).collect();
        assert_eq!(views, vec![RingId(0)]);
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
            RingFamily::Simple,
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
            &[RingGraphEdge {
                source: RingId(0),
                target: RingId(1),
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
}
