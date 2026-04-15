//! Ring detection primitives for GraphIR.
//!
//! Used for ring size queries and bounded ring enumeration.

use std::collections::{BTreeMap, HashSet, VecDeque};

use umol_graph_core::NodeId;
use umol_shared::atom_ast::AromaticValenceAst;
use umol_shared::element::Element;

use super::config::RingEnumerationStrategy;
use crate::ast::{AtomIdx, BondIdx};
use crate::ast::molecule::MoleculeAst;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RingIndex(pub u32);

impl RingIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    atoms: Vec<AtomIdx>,
    bonds: Vec<BondIdx>,
}

impl Ring {
    pub fn new(atoms: Vec<AtomIdx>, bonds: Vec<BondIdx>) -> Result<Self, String> {
        if atoms.len() < 3 {
            return Err("ring must contain at least 3 atoms".to_string());
        }
        if atoms.len() != bonds.len() {
            return Err("ring atoms/bonds length mismatch".to_string());
        }
        Ok(Self { atoms, bonds })
    }

    pub fn atoms(&self) -> &[AtomIdx] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[BondIdx] {
        &self.bonds
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn shared_atoms(&self, other: &Ring) -> Vec<AtomIdx> {
        let (small, large) = if self.atoms.len() <= other.atoms.len() {
            (&self.atoms, &other.atoms)
        } else {
            (&other.atoms, &self.atoms)
        };
        small
            .iter()
            .copied()
            .filter(|atom| large.contains(atom))
            .collect()
    }

    pub fn shared_bonds(&self, other: &Ring) -> Vec<BondIdx> {
        let (small, large) = if self.bonds.len() <= other.bonds.len() {
            (&self.bonds, &other.bonds)
        } else {
            (&other.bonds, &self.bonds)
        };
        small
            .iter()
            .copied()
            .filter(|bond| large.contains(bond))
            .collect()
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
    InducedBenzenoid,
    Mcb,
    Relevant,
    Essential,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingScope {
    All,
    AromaticSubgraph,
    AtomSubset,
}

#[derive(Debug, Clone)]
pub struct RingSet {
    pub family: RingFamily,
    pub scope: RingScope,
    pub max_ring_size: usize,
    pub rings: Vec<Ring>,
    pub atom_to_rings: BTreeMap<AtomIdx, Vec<RingIndex>>,
    pub bond_to_rings: BTreeMap<BondIdx, Vec<RingIndex>>,
    ring_graph: RingGraph,
}

impl RingSet {
    fn empty_with(family: RingFamily, scope: RingScope) -> Self {
        Self {
            family,
            scope,
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

    pub fn empty() -> Self {
        Self::empty_with(RingFamily::Simple, RingScope::All)
    }

    pub fn from_rings(
        family: RingFamily,
        scope: RingScope,
        max_ring_size: usize,
        rings: Vec<Ring>,
    ) -> Self {
        if rings.is_empty() {
            let mut empty = Self::empty_with(family, scope);
            empty.max_ring_size = max_ring_size;
            return empty;
        }

        let mut atom_to_rings: BTreeMap<AtomIdx, Vec<RingIndex>> = BTreeMap::new();
        let mut bond_to_rings: BTreeMap<BondIdx, Vec<RingIndex>> = BTreeMap::new();
        for (idx, ring) in rings.iter().enumerate() {
            let ring_idx = RingIndex(idx as u32);
            for &atom in ring.atoms() {
                atom_to_rings.entry(atom).or_default().push(ring_idx);
            }
            for &bond in ring.bonds() {
                bond_to_rings.entry(bond).or_default().push(ring_idx);
            }
        }

        let ring_graph = RingGraph::from_ring_list(&rings);

        Self {
            family,
            scope,
            max_ring_size,
            rings,
            atom_to_rings,
            bond_to_rings,
            ring_graph,
        }
    }

    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    pub fn ring_indices(&self) -> impl Iterator<Item = RingIndex> {
        (0..self.rings.len()).map(|i| RingIndex(i as u32))
    }

    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn ring(&self, idx: RingIndex) -> Option<&Ring> {
        self.rings.get(idx.index())
    }

    pub fn shared_atoms(&self, a: RingIndex, b: RingIndex) -> Vec<AtomIdx> {
        let (Some(ra), Some(rb)) = (self.ring(a), self.ring(b)) else {
            return Vec::new();
        };
        ra.shared_atoms(rb)
    }

    pub fn shared_bonds(&self, a: RingIndex, b: RingIndex) -> Vec<BondIdx> {
        let (Some(ra), Some(rb)) = (self.ring(a), self.ring(b)) else {
            return Vec::new();
        };
        ra.shared_bonds(rb)
    }

    pub fn ring_relation(&self, a: RingIndex, b: RingIndex) -> RingRelation {
        self.ring_graph.relation(a, b)
    }

    pub fn are_spiro(&self, a: RingIndex, b: RingIndex) -> bool {
        self.ring_relation(a, b) == RingRelation::Spiro
    }

    pub fn are_fused(&self, a: RingIndex, b: RingIndex) -> bool {
        self.ring_relation(a, b) == RingRelation::Fused
    }

    pub fn are_bridged(&self, a: RingIndex, b: RingIndex) -> bool {
        self.ring_relation(a, b) == RingRelation::Bridged
    }

    pub fn ring_spiro_neighbors(&self, i: RingIndex) -> Vec<RingIndex> {
        let mut result: Vec<RingIndex> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Spiro).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn ring_fused_neighbors(&self, i: RingIndex) -> Vec<RingIndex> {
        let mut result: Vec<RingIndex> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Fused).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn ring_bridged_neighbors(&self, i: RingIndex) -> Vec<RingIndex> {
        let mut result: Vec<RingIndex> = self
            .ring_graph
            .neighbors(i)
            .into_iter()
            .filter_map(|(j, relation)| (relation == RingRelation::Bridged).then_some(j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn fused_components(&self) -> Vec<Vec<RingIndex>> {
        let mut visited: HashSet<RingIndex> = HashSet::new();
        let mut components: Vec<Vec<RingIndex>> = Vec::new();

        for ring in self.ring_indices() {
            if visited.contains(&ring) {
                continue;
            }
            let component = self.ring_fused_component(ring);
            for &r in &component {
                visited.insert(r);
            }
            components.push(component);
        }

        components.sort_by_key(|component| component.first().copied().map(RingIndex::index));
        components
    }

    pub fn ring_fused_component(&self, root: RingIndex) -> Vec<RingIndex> {
        let mut visited: HashSet<RingIndex> = HashSet::new();
        let mut queue: VecDeque<RingIndex> = VecDeque::new();
        queue.push_back(root);
        visited.insert(root);

        while let Some(current) = queue.pop_front() {
            for neighbor in self.ring_fused_neighbors(current) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        let mut result: Vec<RingIndex> = visited.into_iter().collect();
        result.sort_unstable();
        result
    }

    pub fn is_ring_atom(&self, atom: AtomIdx) -> bool {
        self.atom_to_rings.contains_key(&atom)
    }

    pub fn atom_smallest_ring_size(&self, atom: AtomIdx) -> Option<usize> {
        self.atom_to_rings.get(&atom).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].len())
                .min()
        })
    }

    pub fn is_ring_bond(&self, bond: BondIdx) -> bool {
        self.bond_to_rings.contains_key(&bond)
    }

    pub fn bond_smallest_ring_size(&self, bond: BondIdx) -> Option<usize> {
        self.bond_to_rings.get(&bond).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].len())
                .min()
        })
    }

    pub fn ring_graph(&self) -> RingGraph {
        self.ring_graph.clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RingGraphEdge {
    pub source: RingIndex,
    pub target: RingIndex,
    pub relation: RingRelation,
}

#[derive(Debug, Clone)]
pub struct RingGraph {
    edges: Vec<RingGraphEdge>,
    neighbors: Vec<Vec<(RingIndex, RingRelation)>>,
}

impl RingGraph {
    pub fn from_ring_list(rings: &[Ring]) -> Self {
        let mut edges = Vec::new();
        let mut neighbors = vec![Vec::new(); rings.len()];
        let indices: Vec<RingIndex> = (0..rings.len()).map(|i| RingIndex(i as u32)).collect();
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

    pub fn neighbors(&self, ring: RingIndex) -> Vec<(RingIndex, RingRelation)> {
        self.neighbors
            .get(ring.index())
            .cloned()
            .unwrap_or_default()
    }

    pub fn relation(&self, a: RingIndex, b: RingIndex) -> RingRelation {
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

pub struct RingEnumerator {
    family: RingFamily,
    aromatic_only: bool,
    max_ring_size: usize,
    max_rings_per_component: usize,
}

impl RingEnumerator {
    pub fn new(family: RingFamily, strategy: &RingEnumerationStrategy) -> Self {
        Self {
            family,
            aromatic_only: strategy.aromatic_only,
            max_ring_size: strategy.max_ring_size,
            max_rings_per_component: strategy.max_rings_per_component,
        }
    }

    pub fn enumerate(&self, ast: &MoleculeAst) -> RingSet {
        let graph = ast.graph();

        let atom_filter: Option<Vec<NodeId>> = match self.family {
            RingFamily::Simple | RingFamily::Induced if self.aromatic_only => {
                Some(
                    graph
                        .node_ids()
                        .filter(|&n| {
                            let atom = ast.atom(AtomIdx(n.0));
                            !matches!(
                                atom.aromatic_valence,
                                AromaticValenceAst::NotAromatic | AromaticValenceAst::Undetermined
                            )
                        })
                        .collect(),
                )
            }
            RingFamily::InducedBenzenoid => {
                Some(
                    graph
                        .node_ids()
                        .filter(|&n| {
                            let atom = ast.atom(AtomIdx(n.0));
                            matches!(atom.element, umol_shared::atom_ast::ElementAst::Lit(Element::C))
                                && !matches!(
                                    atom.aromatic_valence,
                                    AromaticValenceAst::NotAromatic
                                        | AromaticValenceAst::Undetermined
                                )
                        })
                        .collect(),
                )
            }
            _ => None,
        };

        let max_cycle = if self.family == RingFamily::InducedBenzenoid {
            6
        } else {
            self.max_ring_size
        };

        let (sub, node_map, _edge_map) = if let Some(ref nodes) = atom_filter {
            let sub = graph.induced_subgraph(nodes);
            (sub.graph, sub.node_map, sub.edge_map)
        } else {
            let node_map: Vec<NodeId> = graph.node_ids().collect();
            let edge_map: Vec<umol_graph_core::EdgeId> = graph.edge_ids().collect();
            (graph.clone(), node_map, edge_map)
        };

        let bcc = sub.biconnected_components();

        let mut all_rings: Vec<Ring> = Vec::new();
        for component in &bcc {
            let comp_sub = sub.induced_subgraph(component);
            let raw_cycles = comp_sub.graph.enumerate_simple_cycles(max_cycle);

            let mut component_rings: Vec<Ring> = raw_cycles
                .into_iter()
                .filter(|cycle| {
                    if self.family == RingFamily::Induced {
                        is_induced_cycle_graph(&comp_sub.graph, cycle)
                    } else if self.family == RingFamily::InducedBenzenoid {
                        cycle.len() == 6
                    } else {
                        true
                    }
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
                    Ring::new(ring_atoms, ring_bonds).ok()
                })
                .collect();

            component_rings.truncate(self.max_rings_per_component);
            all_rings.extend(component_rings);
        }

        let scope = if atom_filter.is_some() {
            RingScope::AromaticSubgraph
        } else {
            RingScope::All
        };

        RingSet::from_rings(self.family, scope, max_cycle, all_rings)
    }
}

fn is_induced_cycle_graph(graph: &umol_graph_core::Graph, cycle: &[NodeId]) -> bool {
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

fn classify_ring_relation(a: &Ring, b: &Ring) -> RingRelation {
    let shared_bonds = a.shared_bonds(b);
    if shared_bonds.is_empty() {
        return match a.shared_atoms(b).len() {
            0 => RingRelation::Disjoint,
            1 => RingRelation::Spiro,
            _ => RingRelation::MultiSpiro,
        };
    }

    if shared_bonds.len() == 1 {
        return RingRelation::Fused;
    }

    let bonds_a = a.bonds();
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
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst};
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::graph_ir::config::RingEnumerationStrategy;

    fn enumerate_simple(ast: &MoleculeAst, max_ring_size: usize) -> RingSet {
        RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: false,
                max_ring_size,
                max_rings_per_component: usize::MAX,
            },
        )
        .enumerate(ast)
    }

    fn mol(n: usize, edges: &[(usize, usize)]) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    #[fixture]
    fn carbon_ring(#[default(6)] n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    #[fixture]
    fn carbon_chain(#[default(5)] n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n - 1)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx((i + 1) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    #[fixture]
    fn disjoint() -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); 12];
        let mut edges = Vec::new();
        for i in 0..6u32 {
            edges.push((AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)));
        }
        for i in 6..12u32 {
            edges.push((
                AtomIdx(i),
                AtomIdx(6 + ((i + 1 - 6) % 6)),
                BondAst::from_order(1),
            ));
        }
        MoleculeAst::new(atoms, edges, vec![], vec![], vec![], vec![], vec![])
    }

    #[rustfmt::skip]
    #[fixture]
    fn spiro() -> MoleculeAst {
        mol(5, &[(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeAst {
        mol(10, &[
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
        ])
    }

    #[rustfmt::skip]
    #[fixture]
    fn bridged() -> MoleculeAst {
        mol(5, &[(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1)])
    }

    #[rustfmt::skip]
    #[fixture]
    fn multi_spiro() -> MoleculeAst {
        mol(6, &[(0, 1), (1, 2), (0, 3), (3, 2), (0, 4), (4, 2), (0, 5), (5, 2)])
    }

    #[rustfmt::skip]
    #[fixture]
    fn cubane() -> MoleculeAst {
        mol(8, &[
            (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
            (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
        ])
    }

    #[fixture]
    fn substituted() -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); 7];
        let mut edges: Vec<_> = (0..6)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % 6) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        edges.push((AtomIdx(0), AtomIdx(6), BondAst::from_order(1)));
        MoleculeAst::new(atoms, edges, vec![], vec![], vec![], vec![], vec![])
    }

    #[fixture]
    fn aromatic_plus_saturated() -> MoleculeAst {
        let mut atoms: Vec<AtomAst> = (0..6)
            .map(|_| AtomAst {
                element: ElementAst::Lit(Element::C),
                aromatic_valence: AromaticValenceAst::Value(ValueAst::Lit(1)),
                ..Default::default()
            })
            .collect();
        atoms.extend((0..6).map(|_| AtomAst::from_element(Element::C)));
        let mut edges = Vec::new();
        for i in 0..6u32 {
            edges.push((AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)));
        }
        for i in 6..12u32 {
            edges.push((
                AtomIdx(i),
                AtomIdx(6 + ((i + 1 - 6) % 6)),
                BondAst::from_order(1),
            ));
        }
        edges.push((AtomIdx(0), AtomIdx(6), BondAst::from_order(1)));
        MoleculeAst::new(atoms, edges, vec![], vec![], vec![], vec![], vec![])
    }

    #[fixture]
    fn ring_rings(#[default(6)] n: usize) -> RingSet {
        enumerate_simple(&carbon_ring(n), n)
    }

    #[fixture]
    fn disjoint_rings() -> RingSet {
        enumerate_simple(&disjoint(), 10)
    }

    #[fixture]
    fn spiro_rings() -> RingSet {
        enumerate_simple(&spiro(), 10)
    }

    #[fixture]
    fn fused_rings() -> RingSet {
        enumerate_simple(&naphthalene(), 10)
    }

    #[fixture]
    fn bridged_rings() -> RingSet {
        enumerate_simple(&bridged(), 5)
    }

    #[fixture]
    fn cubane_rings() -> RingSet {
        enumerate_simple(&cubane(), 8)
    }

    #[fixture]
    fn multi_spiro_rings() -> RingSet {
        enumerate_simple(&multi_spiro(), 4)
    }

    #[rstest]
    #[case::empty_0(MoleculeAst::default(), 0, 0)]
    #[case::empty_3(MoleculeAst::default(), 3, 0)]
    #[case::single_atom_0(mol(1, &[]), 0, 0)]
    #[case::single_atom_3(mol(1, &[]), 3, 0)]
    #[case::pentane_0(carbon_chain(5), 0, 0)]
    #[case::pentane_3(carbon_chain(5), 3, 0)]
    #[case::cyclohexane_0(carbon_ring(6), 0, 0)]
    #[case::cyclohexane_3(carbon_ring(6), 3, 0)]
    #[case::cyclohexane_6(carbon_ring(6), 6, 1)]
    #[case::cyclohexane_8(carbon_ring(6), 8, 1)]
    #[case::disjoint_0(disjoint(), 0, 0)]
    #[case::disjoint_3(disjoint(), 3, 0)]
    #[case::disjoint_6(disjoint(), 6, 2)]
    #[case::disjoint_8(disjoint(), 8, 2)]
    #[case::spiro_3(spiro(), 3, 2)]
    #[case::spiro_6(spiro(), 6, 2)]
    #[case::fused_3(naphthalene(), 3, 0)]
    #[case::fused_6(naphthalene(), 6, 2)]
    #[case::fused_10(naphthalene(), 10, 3)]
    #[case::fused_20(naphthalene(), 20, 3)]
    #[case::bridged_3(bridged(), 3, 0)]
    #[case::bridged_4(bridged(), 4, 3)]
    #[case::multi_spiro_3(multi_spiro(), 3, 0)]
    #[case::multi_spiro_6(multi_spiro(), 6, 6)]
    #[case::multi_spiro_20(multi_spiro(), 20, 6)]
    #[case::cubane_4(cubane(), 4, 6)]
    #[case::cubane_6(cubane(), 6, 22)]
    #[case::cubane_8(cubane(), 8, 28)]
    #[case::cubane_20(cubane(), 20, 28)]
    fn test_ring_set_enumerate(
        #[case] ast: MoleculeAst,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rings = enumerate_simple(&ast, max_ring_size);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::non_aromatic(carbon_ring(6), 0)]
    #[case::mixed(aromatic_plus_saturated(), 1)]
    fn test_ring_set_enumerate_aromatic(#[case] ast: MoleculeAst, #[case] expected: usize) {
        let rings = RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: true,
                max_ring_size: 22,
                max_rings_per_component: usize::MAX,
            },
        )
        .enumerate(&ast);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::capped_3(8, 3, 3)]
    #[case::uncapped(8, 28, 28)]
    fn test_ring_set_enumerate_capped(
        cubane: MoleculeAst,
        #[case] max_ring_size: usize,
        #[case] max_rings_per_component: usize,
        #[case] expected: usize,
    ) {
        let rings = RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: false,
                max_ring_size,
                max_rings_per_component,
            },
        )
        .enumerate(&cubane);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::ring(ring_rings(6), vec![RingIndex(0)])]
    #[case::fused(fused_rings(), vec![RingIndex(0), RingIndex(1), RingIndex(2)])]
    fn test_ring_set_ring_indices(#[case] rings: RingSet, #[case] expected: Vec<RingIndex>) {
        let indices: Vec<RingIndex> = rings.ring_indices().collect();
        assert_eq!(indices, expected);
    }

    #[rstest]
    #[case::ring(ring_rings(6), RingIndex(0), Some(6))]
    #[case::non_existent(ring_rings(6), RingIndex(1), None)]
    fn test_ring_set_ring(
        #[case] rings: RingSet,
        #[case] idx: RingIndex,
        #[case] expected_len: Option<usize>,
    ) {
        assert_eq!(rings.ring(idx).map(|r| r.len()), expected_len);
    }

    #[rstest]
    #[case::identical(ring_rings(6), RingIndex(0), RingIndex(0), 6)]
    #[case::disjoint(disjoint_rings(), RingIndex(0), RingIndex(1), 0)]
    #[case::spiro(spiro_rings(), RingIndex(0), RingIndex(1), 1)]
    #[case::fused(fused_rings(), RingIndex(0), RingIndex(1), 2)]
    #[case::bridged(bridged_rings(), RingIndex(0), RingIndex(1), 3)]
    #[case::multi_spiro(multi_spiro_rings(), RingIndex(0), RingIndex(5), 2)]
    #[case::cubane(cubane_rings(), RingIndex(0), RingIndex(25), 4)]
    #[case::non_existent(ring_rings(6), RingIndex(0), RingIndex(2), 0)]
    fn test_ring_set_shared_atoms(
        #[case] rings: RingSet,
        #[case] a: RingIndex,
        #[case] b: RingIndex,
        #[case] expected: usize,
    ) {
        assert_eq!(rings.shared_atoms(a, b).len(), expected);
    }

    #[rstest]
    #[case::identical(ring_rings(6), RingIndex(0), RingIndex(0), 6)]
    #[case::disjoint(disjoint_rings(), RingIndex(0), RingIndex(1), 0)]
    #[case::spiro(spiro_rings(), RingIndex(0), RingIndex(1), 0)]
    #[case::fused(fused_rings(), RingIndex(0), RingIndex(1), 1)]
    #[case::bridged(bridged_rings(), RingIndex(0), RingIndex(1), 2)]
    #[case::multi_spiro(multi_spiro_rings(), RingIndex(0), RingIndex(5), 0)]
    #[case::cubane(cubane_rings(), RingIndex(0), RingIndex(25), 2)]
    #[case::non_existent(ring_rings(6), RingIndex(0), RingIndex(2), 0)]
    fn test_ring_set_shared_bonds(
        #[case] rings: RingSet,
        #[case] a: RingIndex,
        #[case] b: RingIndex,
        #[case] expected: usize,
    ) {
        assert_eq!(rings.shared_bonds(a, b).len(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::identical(ring_rings(6), RingIndex(0), RingIndex(0), RingRelation::Identical)]
    #[case::disjoint(disjoint_rings(), RingIndex(0), RingIndex(1), RingRelation::Disjoint)]
    #[case::spiro(spiro_rings(), RingIndex(0), RingIndex(1), RingRelation::Spiro)]
    #[case::fused(fused_rings(), RingIndex(0), RingIndex(1), RingRelation::Fused)]
    #[case::bridged(bridged_rings(), RingIndex(0), RingIndex(1), RingRelation::Bridged)]
    #[case::multi_spiro(multi_spiro_rings(), RingIndex(0), RingIndex(5), RingRelation::MultiSpiro)]
    #[case::cubane(cubane_rings(), RingIndex(0), RingIndex(25), RingRelation::Noncontiguous)]
    #[case::non_existent(ring_rings(6), RingIndex(0), RingIndex(2), RingRelation::Disjoint)]
    fn test_ring_set_ring_relation(
        #[case] rings: RingSet,
        #[case] a: RingIndex,
        #[case] b: RingIndex,
        #[case] expected: RingRelation,
    ) {
        assert_eq!(rings.ring_relation(a, b), expected);
    }

    #[rstest]
    #[case::single_ring(ring_rings(6), RingIndex(0), 0)]
    #[case::spiro_0(spiro_rings(), RingIndex(0), 1)]
    #[case::spiro_1(spiro_rings(), RingIndex(1), 1)]
    #[case::fused_0(fused_rings(), RingIndex(0), 0)]
    #[case::fused_1(fused_rings(), RingIndex(1), 0)]
    #[case::cubane_0(cubane_rings(), RingIndex(0), 0)]
    fn test_ring_set_spiro_neighbors(
        #[case] rings: RingSet,
        #[case] idx: RingIndex,
        #[case] expected: usize,
    ) {
        assert_eq!(rings.ring_spiro_neighbors(idx).len(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single_ring(ring_rings(6), RingIndex(0), 0)]
    #[case::spiro_0(spiro_rings(), RingIndex(0), 0)]
    #[case::spiro_1(spiro_rings(), RingIndex(1), 0)]
    #[case::fused_0(fused_rings(), RingIndex(0), 1)]
    #[case::fused_1(fused_rings(), RingIndex(1), 1)]
    #[case::cubane_0(cubane_rings(), RingIndex(0), 8)]
    fn test_ring_set_fused_neighbors(
        #[case] rings: RingSet,
        #[case] idx: RingIndex,
        #[case] expected: usize,
    ) {
        assert_eq!(rings.ring_fused_neighbors(idx).len(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single_ring(ring_rings(6), RingIndex(0), 0)]
    #[case::fused_0(fused_rings(), RingIndex(0), 1)]
    #[case::fused_1(fused_rings(), RingIndex(1), 1)]
    #[case::bridged_0(bridged_rings(), RingIndex(0), 2)]
    #[case::bridged_1(bridged_rings(), RingIndex(1), 2)]
    #[case::cubane_0(cubane_rings(), RingIndex(0), 16)]
    fn test_ring_set_bridged_neighbors(
        #[case] rings: RingSet,
        #[case] idx: RingIndex,
        #[case] expected: usize,
    ) {
        assert_eq!(rings.ring_bridged_neighbors(idx).len(), expected);
    }

    #[rstest]
    #[case::cyclohexane(ring_rings(6), vec![1])]
    #[case::naphthalene(fused_rings(), vec![1, 2])]
    #[case::cubane(cubane_rings(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 18])]
    fn test_ring_set_fused_components(
        #[case] rings: RingSet,
        #[case] mut expected_sizes: Vec<usize>,
    ) {
        let mut sizes: Vec<usize> = rings.fused_components().iter().map(|c| c.len()).collect();
        sizes.sort_unstable();
        expected_sizes.sort_unstable();
        assert_eq!(sizes, expected_sizes);
    }

    #[rstest]
    #[case::cyclohexane(ring_rings(6), RingIndex(0), 1)]
    #[case::naphthalene_0(fused_rings(), RingIndex(0), 2)]
    #[case::naphthalene_1(fused_rings(), RingIndex(1), 2)]
    #[case::cubane_0(cubane_rings(), RingIndex(0), 18)]
    fn test_ring_set_ring_fused_component(
        #[case] rings: RingSet,
        #[case] idx: RingIndex,
        #[case] expected_size: usize,
    ) {
        let component = rings.ring_fused_component(idx);
        assert!(component.contains(&idx));
        assert_eq!(component.len(), expected_size);
    }

    #[rstest]
    #[case::empty(MoleculeAst::default(), 0, false)]
    #[case::single_atom(mol(1, &[]), 0, false)]
    #[case::cyclohexane_in(carbon_ring(6), 0, true)]
    #[case::cyclohexane_in_3(carbon_ring(6), 3, true)]
    #[case::cyclohexane_out(carbon_ring(6), 6, false)]
    #[case::naphthalene_in(naphthalene(), 0, true)]
    #[case::naphthalene_shared(naphthalene(), 3, true)]
    #[case::naphthalene_second(naphthalene(), 6, true)]
    #[case::naphthalene_out(naphthalene(), 10, false)]
    fn test_ring_set_is_ring_atom(
        #[case] ast: MoleculeAst,
        #[case] atom_index: u32,
        #[case] expected: bool,
    ) {
        let rings = enumerate_simple(&ast, 10);
        assert_eq!(rings.is_ring_atom(AtomIdx(atom_index)), expected);
    }

    #[rstest]
    fn test_ring_set_is_ring_atom_substituted(substituted: MoleculeAst) {
        let rings = enumerate_simple(&substituted, 10);
        for i in 0..6 {
            assert!(rings.is_ring_atom(AtomIdx(i)));
        }
        assert!(!rings.is_ring_atom(AtomIdx(6)));
    }

    #[rstest]
    #[case::empty(MoleculeAst::default(), 10, 0, None)]
    #[case::single_atom(mol(1, &[]), 10, 0, None)]
    #[case::pentane(carbon_chain(5), 10, 0, None)]
    #[case::cyclohexane(carbon_ring(6), 10, 0, Some(6))]
    #[case::cyclohexane_3(carbon_ring(6), 10, 3, Some(6))]
    #[case::naphthalene(naphthalene(), 10, 0, Some(6))]
    #[case::naphthalene_shared(naphthalene(), 10, 3, Some(6))]
    #[case::cubane(cubane(), 10, 0, Some(4))]
    #[case::spiro(spiro(), 10, 0, Some(3))]
    #[case::bridged(bridged(), 10, 0, Some(4))]
    fn test_ring_set_atom_smallest_ring_size(
        #[case] ast: MoleculeAst,
        #[case] max_ring_size: usize,
        #[case] atom_index: u32,
        #[case] expected: Option<usize>,
    ) {
        let rings = enumerate_simple(&ast, max_ring_size);
        assert_eq!(rings.atom_smallest_ring_size(AtomIdx(atom_index)), expected);
    }

    #[rstest]
    #[case::empty(MoleculeAst::default(), 0, false)]
    #[case::single_atom(mol(1, &[]), 0, false)]
    #[case::cyclohexane_in(carbon_ring(6), 0, true)]
    #[case::cyclohexane_in_3(carbon_ring(6), 3, true)]
    fn test_ring_set_is_ring_bond(
        #[case] ast: MoleculeAst,
        #[case] bond_index: u32,
        #[case] expected: bool,
    ) {
        let rings = enumerate_simple(&ast, 10);
        assert_eq!(rings.is_ring_bond(BondIdx(bond_index)), expected);
    }

    #[rstest]
    fn test_ring_set_is_ring_bond_naphthalene(naphthalene: MoleculeAst) {
        let rings = enumerate_simple(&naphthalene, 10);
        for i in 0..11 {
            assert!(rings.is_ring_bond(BondIdx(i)));
        }
    }

    #[rstest]
    fn test_ring_set_is_ring_bond_substituted(substituted: MoleculeAst) {
        let rings = enumerate_simple(&substituted, 10);
        for i in 0..6 {
            assert!(rings.is_ring_bond(BondIdx(i)));
        }
        assert!(!rings.is_ring_bond(BondIdx(6)));
    }

    #[rstest]
    #[case::empty(MoleculeAst::default(), 10, 0, None)]
    #[case::single_atom(mol(1, &[]), 10, 0, None)]
    #[case::pentane(carbon_chain(5), 10, 0, None)]
    #[case::cyclohexane(carbon_ring(6), 10, 0, Some(6))]
    #[case::cyclohexane_3(carbon_ring(6), 10, 3, Some(6))]
    #[case::cubane(cubane(), 10, 0, Some(4))]
    #[case::spiro(spiro(), 10, 0, Some(3))]
    #[case::bridged(bridged(), 10, 0, Some(4))]
    fn test_ring_set_bond_smallest_ring_size(
        #[case] ast: MoleculeAst,
        #[case] max_ring_size: usize,
        #[case] bond_index: u32,
        #[case] expected: Option<usize>,
    ) {
        let rings = enumerate_simple(&ast, max_ring_size);
        assert_eq!(rings.bond_smallest_ring_size(BondIdx(bond_index)), expected);
    }

    #[rstest]
    fn test_ring_set_ring_graph(fused_rings: RingSet) {
        let ring_graph = fused_rings.ring_graph();
        assert!(!ring_graph.edges().is_empty());
        assert!(ring_graph
            .edges()
            .iter()
            .any(|e| matches!(e.relation, RingRelation::Fused)));
    }
}
