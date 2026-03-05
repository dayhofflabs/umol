//! Topology projections for GraphIR molecules.
//!
//! `TopologyGraph` stores a projected undirected graph and explicit node/edge
//! mappings back to GraphIR entities.

use std::collections::{HashMap, HashSet, VecDeque};

use nalgebra::DMatrix;
use petgraph::graph::{EdgeIndex, Graph, NodeIndex};
use petgraph::prelude::*;
use petgraph::visit::EdgeRef;
use thiserror::Error;

use super::builder::MoleculeBuilder;
use super::{
    AtomIndex, BondIndex, DativeBondIndex, Molecule, MulticenterBondIndex, NoncovalentBondIndex,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionKind {
    Ordinary,
    LineGraph,
    BipartiteIncidence,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DativeProjection {
    Skip,
    Undirected,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoncovalentProjection {
    Skip,
    Undirected,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MulticenterProjection {
    Skip,
    CliqueExpansion,
    IncidenceNode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopologyProjection {
    pub kind: ProjectionKind,
    pub dative: DativeProjection,
    pub noncovalent: NoncovalentProjection,
    pub multicenter: MulticenterProjection,
}

impl TopologyProjection {
    pub fn ordinary() -> Self {
        Self {
            kind: ProjectionKind::Ordinary,
            dative: DativeProjection::Skip,
            noncovalent: NoncovalentProjection::Skip,
            multicenter: MulticenterProjection::Skip,
        }
    }

    pub fn line_graph() -> Self {
        Self {
            kind: ProjectionKind::LineGraph,
            ..Self::ordinary()
        }
    }

    pub fn bipartite_incidence() -> Self {
        Self {
            kind: ProjectionKind::BipartiteIncidence,
            ..Self::ordinary()
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TopologyNodeRef {
    Atom(AtomIndex),
    Bond(BondIndex),
    DativeBond(DativeBondIndex),
    NoncovalentBond(NoncovalentBondIndex),
    MulticenterBond(MulticenterBondIndex),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TopologyEdgeRef {
    BondLikeConnection,
    SharedAtom(AtomIndex),
    Incidence,
}

#[derive(Debug, Clone)]
pub(crate) enum TopologyEdge {
    Edge {
        node_ref: TopologyNodeRef,
        a: AtomIndex,
        b: AtomIndex,
    },
    Hyperedge {
        node_ref: TopologyNodeRef,
        atoms: Vec<AtomIndex>,
    },
}

#[derive(Debug, Clone)]
pub struct TopologyGraph {
    graph: Graph<(), (), Undirected>,
    node_map: Vec<TopologyNodeRef>,
    edge_map: Vec<TopologyEdgeRef>,
}

#[derive(Debug, Error, Clone, Copy, PartialEq, Eq)]
pub enum TopologyExportError {
    #[error("graph6 export requires a simple graph without incidence edges")]
    HasIncidenceEdges,
    #[error("graph6 export requires a simple graph without self-loops")]
    HasSelfLoops,
    #[error("graph6 export requires a simple graph without parallel edges")]
    HasParallelEdges,
    #[error("graph has too many nodes for graph6: {0}")]
    TooManyNodes(usize),
}

impl TopologyGraph {
    pub fn graph(&self) -> &Graph<(), (), Undirected> {
        &self.graph
    }

    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn node_map(&self) -> &[TopologyNodeRef] {
        self.node_map.as_slice()
    }

    pub fn edge_map(&self) -> &[TopologyEdgeRef] {
        self.edge_map.as_slice()
    }

    pub fn node_ref(&self, index: NodeIndex) -> Option<TopologyNodeRef> {
        self.node_map.get(index.index()).copied()
    }

    pub fn edge_ref(&self, index: EdgeIndex) -> Option<TopologyEdgeRef> {
        self.edge_map.get(index.index()).copied()
    }

    pub fn from_molecule(mol: &Molecule, projection: TopologyProjection) -> Self {
        Self::build(
            projection,
            mol.topology_nodes(),
            mol.topology_edges(projection),
        )
    }

    pub fn from_builder(builder: &MoleculeBuilder, projection: TopologyProjection) -> Self {
        Self::build(
            projection,
            builder.topology_nodes(),
            builder.topology_edges(projection),
        )
    }

    pub fn adjacency_matrix(&self) -> DMatrix<u8> {
        let n = self.node_count();
        let mut matrix = DMatrix::<u8>::zeros(n, n);
        for edge in self.graph.edge_references() {
            let u = edge.source().index();
            let v = edge.target().index();
            matrix[(u, v)] = matrix[(u, v)].saturating_add(1);
            if u != v {
                matrix[(v, u)] = matrix[(v, u)].saturating_add(1);
            }
        }
        matrix
    }

    pub fn incidence_matrix(&self) -> DMatrix<u8> {
        let n = self.node_count();
        let m = self.edge_count();
        let mut matrix = DMatrix::<u8>::zeros(n, m);
        for (col, edge) in self.graph.edge_references().enumerate() {
            let u = edge.source().index();
            let v = edge.target().index();
            matrix[(u, col)] = matrix[(u, col)].saturating_add(1);
            matrix[(v, col)] = matrix[(v, col)].saturating_add(1);
        }
        matrix
    }

    pub fn canonical_bfs(&self) -> Vec<NodeIndex> {
        self.canonical_bfs_with_rank(|_| usize::MAX)
    }

    pub fn canonical_bfs_with_rank<F>(&self, rank: F) -> Vec<NodeIndex>
    where
        F: Fn(TopologyNodeRef) -> usize,
    {
        let n = self.node_count();
        if n == 0 {
            return Vec::new();
        }
        let mut order = Vec::with_capacity(n);
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();

        let mut all_nodes: Vec<NodeIndex> = self.graph.node_indices().collect();
        all_nodes.sort_by_key(|&nidx| self.bfs_key(nidx, &rank));

        for start in all_nodes {
            if visited[start.index()] {
                continue;
            }
            visited[start.index()] = true;
            queue.push_back(start);

            while let Some(node) = queue.pop_front() {
                order.push(node);
                let mut neighbors: Vec<NodeIndex> = self
                    .graph
                    .neighbors(node)
                    .filter(|nidx| !visited[nidx.index()])
                    .collect();
                neighbors.sort_by_key(|&nidx| self.bfs_key(nidx, &rank));
                neighbors.dedup_by_key(|nidx| nidx.index());
                for neigh in neighbors {
                    if !visited[neigh.index()] {
                        visited[neigh.index()] = true;
                        queue.push_back(neigh);
                    }
                }
            }
        }

        order
    }

    pub fn to_graph6(&self) -> Result<String, TopologyExportError> {
        self.validate_graph6_export()?;
        let order: Vec<NodeIndex> = (0..self.node_count()).map(NodeIndex::new).collect();
        Ok(self.encode_graph6_with_order(&order))
    }

    pub fn to_graph6_canonical(
        &self,
    ) -> Result<(String, Vec<NodeIndex>), TopologyExportError> {
        self.validate_graph6_export()?;
        let order = self.canonical_bfs();
        let g6 = self.encode_graph6_with_order(&order);
        Ok((g6, order))
    }

    pub fn to_graph6_canonical_with_rank<F>(
        &self,
        rank: F,
    ) -> Result<(String, Vec<NodeIndex>), TopologyExportError>
    where
        F: Fn(TopologyNodeRef) -> usize,
    {
        self.validate_graph6_export()?;
        let order = self.canonical_bfs_with_rank(rank);
        let g6 = self.encode_graph6_with_order(&order);
        Ok((g6, order))
    }

    fn build(
        projection: TopologyProjection,
        atoms: impl Iterator<Item = AtomIndex>,
        edges: impl Iterator<Item = TopologyEdge>,
    ) -> TopologyGraph {
        match projection.kind {
            ProjectionKind::Ordinary => Self::build_ordinary(atoms, edges),
            ProjectionKind::LineGraph => Self::build_line_graph(edges),
            ProjectionKind::BipartiteIncidence => Self::build_bipartite(atoms, edges),
        }
    }

    fn build_ordinary(
        atoms: impl Iterator<Item = AtomIndex>,
        edges: impl Iterator<Item = TopologyEdge>,
    ) -> TopologyGraph {
        let mut graph = Graph::<(), (), Undirected>::default();
        let mut node_map = Vec::new();
        let mut edge_map = Vec::new();
        let mut atom_node = HashMap::<AtomIndex, NodeIndex>::new();
        for atom in atoms {
            let n = graph.add_node(());
            atom_node.insert(atom, n);
            node_map.push(TopologyNodeRef::Atom(atom));
        }

        for edge in edges {
            match edge {
                TopologyEdge::Edge { a, b, .. } => {
                    if let (Some(&na), Some(&nb)) = (atom_node.get(&a), atom_node.get(&b)) {
                        graph.add_edge(na, nb, ());
                        edge_map.push(TopologyEdgeRef::BondLikeConnection);
                    }
                }
                TopologyEdge::Hyperedge { node_ref, atoms } => {
                    let hn = graph.add_node(());
                    node_map.push(node_ref);
                    for atom in atoms {
                        if let Some(&an) = atom_node.get(&atom) {
                            graph.add_edge(an, hn, ());
                            edge_map.push(TopologyEdgeRef::Incidence);
                        }
                    }
                }
            }
        }

        TopologyGraph {
            graph,
            node_map,
            edge_map,
        }
    }

    fn build_line_graph(edges: impl Iterator<Item = TopologyEdge>) -> TopologyGraph {
        let mut graph = Graph::<(), (), Undirected>::default();
        let mut node_map = Vec::new();
        let mut edge_map = Vec::new();
        let mut line_nodes = Vec::<(NodeIndex, Vec<AtomIndex>)>::new();

        for edge in edges {
            let (node_ref, atoms) = match edge {
                TopologyEdge::Edge { node_ref, a, b } => (node_ref, vec![a, b]),
                TopologyEdge::Hyperedge { node_ref, atoms } => (node_ref, atoms),
            };
            let n = graph.add_node(());
            node_map.push(node_ref);
            line_nodes.push((n, atoms));
        }

        let mut incidence = HashMap::<AtomIndex, Vec<NodeIndex>>::new();
        for (node, atoms) in &line_nodes {
            for atom in atoms {
                incidence.entry(*atom).or_default().push(*node);
            }
        }

        let mut by_atom: Vec<(AtomIndex, Vec<NodeIndex>)> = incidence.into_iter().collect();
        by_atom.sort_unstable_by_key(|(a, _)| a.index());

        let mut seen = HashSet::<(usize, usize)>::new();
        for (atom, nodes) in by_atom {
            for i in 0..nodes.len() {
                for j in (i + 1)..nodes.len() {
                    let u = nodes[i];
                    let v = nodes[j];
                    let key = if u.index() <= v.index() {
                        (u.index(), v.index())
                    } else {
                        (v.index(), u.index())
                    };
                    if seen.insert(key) {
                        graph.add_edge(u, v, ());
                        edge_map.push(TopologyEdgeRef::SharedAtom(atom));
                    }
                }
            }
        }

        TopologyGraph {
            graph,
            node_map,
            edge_map,
        }
    }

    fn build_bipartite(
        atoms: impl Iterator<Item = AtomIndex>,
        edges: impl Iterator<Item = TopologyEdge>,
    ) -> TopologyGraph {
        let mut graph = Graph::<(), (), Undirected>::default();
        let mut node_map = Vec::new();
        let mut edge_map = Vec::new();
        let mut atom_node = HashMap::<AtomIndex, NodeIndex>::new();
        for atom in atoms {
            let n = graph.add_node(());
            atom_node.insert(atom, n);
            node_map.push(TopologyNodeRef::Atom(atom));
        }

        for edge in edges {
            match edge {
                TopologyEdge::Edge { node_ref, a, b } => {
                    let en = graph.add_node(());
                    node_map.push(node_ref);
                    if let Some(&an) = atom_node.get(&a) {
                        graph.add_edge(an, en, ());
                        edge_map.push(TopologyEdgeRef::Incidence);
                    }
                    if let Some(&bn) = atom_node.get(&b) {
                        graph.add_edge(bn, en, ());
                        edge_map.push(TopologyEdgeRef::Incidence);
                    }
                }
                TopologyEdge::Hyperedge { node_ref, atoms } => {
                    let en = graph.add_node(());
                    node_map.push(node_ref);
                    for atom in atoms {
                        if let Some(&an) = atom_node.get(&atom) {
                            graph.add_edge(an, en, ());
                            edge_map.push(TopologyEdgeRef::Incidence);
                        }
                    }
                }
            }
        }

        TopologyGraph {
            graph,
            node_map,
            edge_map,
        }
    }

    fn bfs_key<F>(&self, nidx: NodeIndex, rank: &F) -> (usize, usize, usize)
    where
        F: Fn(TopologyNodeRef) -> usize,
    {
        let node_ref = self.node_map[nidx.index()];
        (rank(node_ref), node_ref_key(node_ref), nidx.index())
    }

    fn append_graph6_n(out: &mut String, n: u64) {
        if n <= 62 {
            out.push((n as u8 + 63) as char);
            return;
        }
        if n <= 258_047 {
            out.push('~');
            out.push((((n >> 12) & 63) as u8 + 63) as char);
            out.push((((n >> 6) & 63) as u8 + 63) as char);
            out.push(((n & 63) as u8 + 63) as char);
            return;
        }
        out.push('~');
        out.push('~');
        for shift in [30, 24, 18, 12, 6, 0] {
            out.push((((n >> shift) & 63) as u8 + 63) as char);
        }
    }

    fn validate_graph6_export(&self) -> Result<(), TopologyExportError> {
        if self
            .edge_map
            .iter()
            .any(|e| *e == TopologyEdgeRef::Incidence)
        {
            return Err(TopologyExportError::HasIncidenceEdges);
        }
        let n = self.node_count();
        let n_u64 = u64::try_from(n).map_err(|_| TopologyExportError::TooManyNodes(n))?;
        if n_u64 > 68_719_476_735 {
            return Err(TopologyExportError::TooManyNodes(n));
        }
        let mut seen = HashSet::<(usize, usize)>::new();
        for edge in self.graph.edge_references() {
            let u = edge.source().index();
            let v = edge.target().index();
            if u == v {
                return Err(TopologyExportError::HasSelfLoops);
            }
            let key = if u <= v { (u, v) } else { (v, u) };
            if !seen.insert(key) {
                return Err(TopologyExportError::HasParallelEdges);
            }
        }
        Ok(())
    }

    fn encode_graph6_with_order(&self, order: &[NodeIndex]) -> String {
        let n = order.len();
        let mut out = String::new();
        Self::append_graph6_n(&mut out, n as u64);

        let mut bits = Vec::<bool>::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                bits.push(self.graph.find_edge(order[i], order[j]).is_some());
            }
        }
        while bits.len() % 6 != 0 {
            bits.push(false);
        }
        for chunk in bits.chunks(6) {
            let mut value = 0u8;
            for (k, bit) in chunk.iter().enumerate() {
                if *bit {
                    value |= 1 << (5 - k);
                }
            }
            out.push((value + 63) as char);
        }
        out
    }
}

fn node_ref_key(node_ref: TopologyNodeRef) -> usize {
    match node_ref {
        TopologyNodeRef::Atom(i) => i.index(),
        TopologyNodeRef::Bond(i) => 1_000_000 + i.index(),
        TopologyNodeRef::DativeBond(i) => 2_000_000 + i.index(),
        TopologyNodeRef::NoncovalentBond(i) => 3_000_000 + i.index(),
        TopologyNodeRef::MulticenterBond(i) => 4_000_000 + i.index(),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::bond::BondNoncovalent;
    use crate::graph_ir::atom::AtomBuilder;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::dative::DativeBond;
    use crate::graph_ir::molecule::MoleculeBuilder;
    use crate::graph_ir::multicenter::{MulticenterBond, MulticenterSet};
    use crate::graph_ir::noncovalent::NoncovalentBond;

    #[fixture]
    fn c3_molecule() -> MoleculeBuilder {
        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(AtomBuilder::new(Element::C));
        let a1 = b.add_atom(AtomBuilder::new(Element::C));
        let a2 = b.add_atom(AtomBuilder::new(Element::C));
        b.add_bond_unchecked(a0, a1, BondBuilder::new(1, Some(false)));
        b.add_bond_unchecked(a1, a2, BondBuilder::new(1, Some(false)));
        b
    }

    #[fixture]
    fn c3_molecule_reordered() -> MoleculeBuilder {
        let mut b = MoleculeBuilder::new();
        let a2 = b.add_atom(AtomBuilder::new(Element::C));
        let a1 = b.add_atom(AtomBuilder::new(Element::C));
        let a0 = b.add_atom(AtomBuilder::new(Element::C));
        b.add_bond_unchecked(a0, a1, BondBuilder::new(1, Some(false)));
        b.add_bond_unchecked(a1, a2, BondBuilder::new(1, Some(false)));
        b
    }

    #[fixture]
    fn n_c_dative_molecule() -> MoleculeBuilder {
        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(AtomBuilder::new(Element::N));
        let a1 = b.add_atom(AtomBuilder::new(Element::C));
        b.add_dative_bond(DativeBond::new(a0, a1, 1));
        b
    }

    #[fixture]
    fn bh2_multicenter_molecule() -> MoleculeBuilder {
        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(AtomBuilder::new(Element::B));
        let a1 = b.add_atom(AtomBuilder::new(Element::H));
        let a2 = b.add_atom(AtomBuilder::new(Element::H));
        b.add_multicenter_bond(MulticenterBond::new([MulticenterSet::topology_only([
            a0, a1, a2,
        ])]));
        b
    }

    #[fixture]
    fn o_h_noncovalent_molecule() -> MoleculeBuilder {
        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(AtomBuilder::new(Element::O));
        let a1 = b.add_atom(AtomBuilder::new(Element::H));
        b.add_noncovalent_bond(NoncovalentBond::new(a0, a1, BondNoncovalent::Hydrogen));
        b
    }

    #[rstest]
    fn test_ordinary(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        assert_eq!(tg.node_count(), 3);
        assert_eq!(tg.edge_count(), 2);
        assert_eq!(tg.edge_map()[0], TopologyEdgeRef::BondLikeConnection);
    }

    #[rstest]
    fn test_ordinary_dative(n_c_dative_molecule: MoleculeBuilder) {
        let mut p = TopologyProjection::ordinary();
        p.dative = DativeProjection::Undirected;
        let tg = n_c_dative_molecule.topology_graph(p);
        assert_eq!(tg.node_count(), 2);
        assert_eq!(tg.edge_count(), 1);
    }

    #[rstest]
    fn test_line_graph(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection {
            kind: ProjectionKind::LineGraph,
            ..TopologyProjection::ordinary()
        });
        assert_eq!(tg.node_count(), 2);
        assert_eq!(tg.edge_count(), 1);
        assert_eq!(
            tg.edge_map()[0],
            TopologyEdgeRef::SharedAtom(AtomIndex::new(1))
        );
    }

    #[rstest]
    fn test_bipartite_incidence(o_h_noncovalent_molecule: MoleculeBuilder) {
        let mut p = TopologyProjection::bipartite_incidence();
        p.noncovalent = NoncovalentProjection::Undirected;
        let tg = o_h_noncovalent_molecule.topology_graph(p);
        assert_eq!(tg.node_count(), 3);
        assert_eq!(tg.edge_count(), 2);
        assert!(tg
            .edge_map()
            .iter()
            .all(|e| *e == TopologyEdgeRef::Incidence));
    }

    #[rstest]
    fn test_multicenter_incidence_node(bh2_multicenter_molecule: MoleculeBuilder) {
        let mut p = TopologyProjection::ordinary();
        p.multicenter = MulticenterProjection::IncidenceNode;
        let tg = bh2_multicenter_molecule.topology_graph(p);
        assert_eq!(tg.node_count(), 4);
        assert_eq!(tg.edge_count(), 3);
        assert!(tg
            .edge_map()
            .iter()
            .all(|e| *e == TopologyEdgeRef::Incidence));
    }

    #[rstest]
    fn test_multicenter_clique_expansion(bh2_multicenter_molecule: MoleculeBuilder) {
        let mut p = TopologyProjection::ordinary();
        p.multicenter = MulticenterProjection::CliqueExpansion;
        let edges: Vec<TopologyEdge> = bh2_multicenter_molecule.topology_edges(p).collect();
        assert_eq!(edges.len(), 3);
        assert!(edges.iter().all(|e| matches!(
            e,
            TopologyEdge::Edge {
                node_ref: TopologyNodeRef::MulticenterBond(MulticenterBondIndex(0)),
                ..
            }
        )));
    }

    #[rstest]
    fn test_adjacency_matrix(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        assert_eq!(
            tg.adjacency_matrix(),
            DMatrix::from_row_slice(3, 3, &[0, 1, 0, 1, 0, 1, 0, 1, 0])
        );
    }

    #[rstest]
    fn test_incidence_matrix(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        assert_eq!(
            tg.incidence_matrix(),
            DMatrix::from_row_slice(3, 2, &[1, 0, 1, 1, 0, 1])
        );
    }

    #[rstest]
    fn test_canonical_bfs_default(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        let bfs = tg.canonical_bfs();
        assert_eq!(
            bfs.into_iter().map(|n| n.index()).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
    }

    #[rstest]
    fn test_canonical_bfs_with_rank(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        let bfs = tg.canonical_bfs_with_rank(|node_ref| match node_ref {
            TopologyNodeRef::Atom(a) => usize::MAX - a.index(),
            _ => usize::MAX,
        });
        assert_eq!(
            bfs.into_iter().map(|n| n.index()).collect::<Vec<_>>(),
            vec![2, 1, 0]
        );
    }

    #[rstest]
    fn test_to_graph6(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        assert_eq!(tg.to_graph6().unwrap(), "Bg");
    }

    #[rstest]
    fn test_to_graph6_canonical(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        let (g6, _order) = tg.to_graph6_canonical().unwrap();
        assert_eq!(g6, "Bg");
    }

    #[rstest]
    fn test_to_graph6_canonical_with_rank(c3_molecule: MoleculeBuilder) {
        let tg = c3_molecule.topology_graph(TopologyProjection::ordinary());
        let (g6, _order) = tg
            .to_graph6_canonical_with_rank(|node_ref| match node_ref {
                TopologyNodeRef::Atom(a) => usize::MAX - a.index(),
                _ => usize::MAX,
            })
            .unwrap();
        assert_eq!(g6, "Bg");
    }

    #[rstest]
    fn test_builder_topology_graph6_canonical(
        c3_molecule: MoleculeBuilder,
        c3_molecule_reordered: MoleculeBuilder,
    ) {
        let p = TopologyProjection::ordinary();
        let a = c3_molecule.topology_graph6_canonical(p).unwrap();
        let b = c3_molecule_reordered.topology_graph6_canonical(p).unwrap();
        assert_eq!(a, b);
    }

    #[rstest]
    fn test_to_graph6_error(bh2_multicenter_molecule: MoleculeBuilder) {
        let mut p = TopologyProjection::ordinary();
        p.multicenter = MulticenterProjection::IncidenceNode;
        let tg = bh2_multicenter_molecule.topology_graph(p);
        assert_eq!(tg.to_graph6(), Err(TopologyExportError::HasIncidenceEdges));
    }
}
