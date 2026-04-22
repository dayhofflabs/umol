//! Ring detection and ring-system analysis.
//!
//! Ring enumeration uses Vismara's relevant-cycle algorithm from
//! `umol_graph_core`, decomposed over biconnected components.

use std::collections::{BTreeMap, HashSet, VecDeque};

use umol_graph_core::{
    BiconnectedComponentsAlgorithm, CycleEnumerationAlgorithm, Graph, NodeId,
};

use super::idx::{AtomIdx, BondIdx};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RingIdx(pub u32);

impl RingIdx {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Ring {
    atoms: Vec<AtomIdx>,
    bonds: Vec<BondIdx>,
}

impl Ring {
    fn new(atoms: Vec<AtomIdx>, bonds: Vec<BondIdx>) -> Option<Self> {
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

#[derive(Debug, Clone, Copy)]
pub struct RingView<'a> {
    pub idx: RingIdx,
    atoms: &'a [AtomIdx],
    bonds: &'a [BondIdx],
}

impl<'a> RingView<'a> {
    pub fn atoms(&self) -> &'a [AtomIdx] {
        self.atoms
    }

    pub fn bonds(&self) -> &'a [BondIdx] {
        self.bonds
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn shared_atoms(&self, other: &RingView<'_>) -> Vec<AtomIdx> {
        intersection(self.atoms, other.atoms)
    }

    pub fn shared_bonds(&self, other: &RingView<'_>) -> Vec<BondIdx> {
        intersection(self.bonds, other.bonds)
    }
}

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingFamily {
    Simple,
    Induced,
}

#[derive(Debug, Clone)]
pub struct RingSet {
    family: RingFamily,
    max_ring_size: usize,
    rings: Vec<Ring>,
    atom_to_rings: BTreeMap<AtomIdx, Vec<RingIdx>>,
    bond_to_rings: BTreeMap<BondIdx, Vec<RingIdx>>,
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

        let mut atom_to_rings: BTreeMap<AtomIdx, Vec<RingIdx>> = BTreeMap::new();
        let mut bond_to_rings: BTreeMap<BondIdx, Vec<RingIdx>> = BTreeMap::new();
        for (idx, ring) in rings.iter().enumerate() {
            let ring_idx = RingIdx(idx as u32);
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

    pub fn ids(&self) -> impl Iterator<Item = RingIdx> {
        (0..self.rings.len()).map(|i| RingIdx(i as u32))
    }

    pub fn iter(&self) -> impl Iterator<Item = RingView<'_>> {
        self.rings.iter().enumerate().map(|(i, r)| RingView {
            idx: RingIdx(i as u32),
            atoms: &r.atoms,
            bonds: &r.bonds,
        })
    }

    pub fn get(&self, idx: RingIdx) -> Option<RingView<'_>> {
        let r = self.rings.get(idx.index())?;
        Some(RingView {
            idx,
            atoms: &r.atoms,
            bonds: &r.bonds,
        })
    }

    pub fn shared_atoms(&self, a: RingIdx, b: RingIdx) -> Vec<AtomIdx> {
        let (Some(ra), Some(rb)) = (self.rings.get(a.index()), self.rings.get(b.index())) else {
            return Vec::new();
        };
        intersection(&ra.atoms, &rb.atoms)
    }

    pub fn shared_bonds(&self, a: RingIdx, b: RingIdx) -> Vec<BondIdx> {
        let (Some(ra), Some(rb)) = (self.rings.get(a.index()), self.rings.get(b.index())) else {
            return Vec::new();
        };
        intersection(&ra.bonds, &rb.bonds)
    }

    pub fn relation(&self, a: RingIdx, b: RingIdx) -> RingRelation {
        self.ring_graph.relation(a, b)
    }

    pub fn are_spiro(&self, a: RingIdx, b: RingIdx) -> bool {
        self.relation(a, b) == RingRelation::Spiro
    }

    pub fn are_fused(&self, a: RingIdx, b: RingIdx) -> bool {
        self.relation(a, b) == RingRelation::Fused
    }

    pub fn are_bridged(&self, a: RingIdx, b: RingIdx) -> bool {
        self.relation(a, b) == RingRelation::Bridged
    }

    pub fn spiro_neighbors(&self, i: RingIdx) -> Vec<RingIdx> {
        let mut result: Vec<RingIdx> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Spiro).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn fused_neighbors(&self, i: RingIdx) -> Vec<RingIdx> {
        let mut result: Vec<RingIdx> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Fused).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn bridged_neighbors(&self, i: RingIdx) -> Vec<RingIdx> {
        let mut result: Vec<RingIdx> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Bridged).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn fused_components(&self) -> Vec<Vec<RingIdx>> {
        let mut visited: HashSet<RingIdx> = HashSet::new();
        let mut components: Vec<Vec<RingIdx>> = Vec::new();

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

        components.sort_by_key(|component| component.first().copied().map(RingIdx::index));
        components
    }

    pub fn fused_component(&self, root: RingIdx) -> Vec<RingIdx> {
        let mut visited: HashSet<RingIdx> = HashSet::new();
        let mut queue: VecDeque<RingIdx> = VecDeque::new();
        queue.push_back(root);
        visited.insert(root);

        while let Some(current) = queue.pop_front() {
            for neighbor in self.fused_neighbors(current) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        let mut result: Vec<RingIdx> = visited.into_iter().collect();
        result.sort_unstable();
        result
    }

    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.atom_to_rings.contains_key(&atom)
    }

    pub fn atom_smallest_ring_size(&self, atom: AtomIdx) -> Option<usize> {
        self.atom_to_rings.get(&atom).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].atoms.len())
                .min()
        })
    }

    pub fn contains_bond(&self, bond: BondIdx) -> bool {
        self.bond_to_rings.contains_key(&bond)
    }

    pub fn bond_smallest_ring_size(&self, bond: BondIdx) -> Option<usize> {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RingGraphEdge {
    pub source: RingIdx,
    pub target: RingIdx,
    pub relation: RingRelation,
}

#[derive(Debug, Clone)]
pub struct RingGraph {
    edges: Vec<RingGraphEdge>,
    neighbors: Vec<Vec<(RingIdx, RingRelation)>>,
}

impl RingGraph {
    fn from_ring_list(rings: &[Ring]) -> Self {
        let mut edges = Vec::new();
        let mut neighbors = vec![Vec::new(); rings.len()];
        let indices: Vec<RingIdx> = (0..rings.len()).map(|i| RingIdx(i as u32)).collect();
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
            n.sort_by_key(|(idx, rel)| (*idx, *rel as u8));
        }
        Self { edges, neighbors }
    }

    pub fn edges(&self) -> &[RingGraphEdge] {
        &self.edges
    }

    pub fn neighbors(&self, ring: RingIdx) -> Vec<(RingIdx, RingRelation)> {
        self.neighbors
            .get(ring.index())
            .cloned()
            .unwrap_or_default()
    }

    pub fn relation(&self, a: RingIdx, b: RingIdx) -> RingRelation {
        if a == b {
            return RingRelation::Identical;
        }
        self.neighbors
            .get(a.index())
            .and_then(|neighbors| {
                neighbors
                    .iter()
                    .find_map(|(idx, rel)| (*idx == b).then_some(*rel))
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
    atom_filter: impl Fn(AtomIdx) -> bool,
) -> RingSet {
    let filtered_nodes: Vec<NodeId> = graph
        .node_ids()
        .filter(|&n| atom_filter(AtomIdx::from(n)))
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
        let raw_cycles =
            comp_sub
                .graph
                .enumerate_cycles(max_ring_size, CycleEnumerationAlgorithm::Vismara);

        let component_rings: Vec<Ring> = raw_cycles
            .into_iter()
            .filter(|cycle| match family {
                RingFamily::Induced => is_induced_cycle(&comp_sub.graph, cycle),
                RingFamily::Simple => true,
            })
            .filter_map(|cycle| {
                let ring_atoms: Vec<AtomIdx> = cycle
                    .iter()
                    .map(|&local| {
                        let sub_node = comp_sub.node_map[local.index()];
                        let orig_node = node_map[sub_node.index()];
                        AtomIdx::from(orig_node)
                    })
                    .collect();
                let n = ring_atoms.len();
                let mut ring_bonds = Vec::with_capacity(n);
                for i in 0..n {
                    let a_orig = NodeId::from(ring_atoms[i]);
                    let b_orig = NodeId::from(ring_atoms[(i + 1) % n]);
                    let edge = graph.find_edge(a_orig, b_orig)?;
                    ring_bonds.push(BondIdx::from(edge));
                }
                Ring::new(ring_atoms, ring_bonds)
            })
            .collect();

        all_rings.extend(component_rings);
    }

    RingSet::from_rings(family, max_ring_size, all_rings)
}
