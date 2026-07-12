//! Graph automorphism and canonical labeling.

use std::cell::RefCell;
use std::collections::HashSet;
use std::mem;
use std::os::raw::c_int;

use nauty_Traces_sys::*;

use crate::graph::{EdgeId, Graph, NodeId};

thread_local! {
    /// Accumulates the generators nauty emits via `userautomproc` during one
    /// `sparsenauty` call; cleared before the call and drained right after. The
    /// crate enables nauty's `tls` feature, so each thread runs independently.
    static GENERATORS: RefCell<Vec<Vec<NodeId>>> = const { RefCell::new(Vec::new()) };
}

/// nauty `userautomproc`: invoked once per generator with its permutation image
/// `perm` over the `n` vertices (`perm[i]` is the image of vertex `i`).
unsafe extern "C" fn capture_generator(
    _count: c_int,
    perm: *mut c_int,
    _orbits: *mut c_int,
    _numorbits: c_int,
    _stabvertex: c_int,
    n: c_int,
) {
    let image: Vec<NodeId> = (0..n as usize)
        .map(|i| NodeId(unsafe { *perm.add(i) } as u32))
        .collect();
    GENERATORS.with(|g| g.borrow_mut().push(image));
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AutomorphismAlgorithm {
    Nauty,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AutoGroupOrder {
    Exact(u32),
    Approx(f64),
}

#[derive(Debug, Clone)]
pub struct Automorphism {
    orbits: Vec<NodeId>,
    canonical_lab: Vec<NodeId>,
    node_count: usize,
    orbit_count: usize,
    group_order: AutoGroupOrder,
    generators: Vec<Vec<NodeId>>,
}

impl Graph {
    pub fn automorphisms<C: Ord + Copy>(
        &self,
        node_color: impl Fn(NodeId) -> C,
        alg: AutomorphismAlgorithm,
    ) -> Automorphism {
        match alg {
            AutomorphismAlgorithm::Nauty => self.automorphisms_nauty(node_color),
        }
    }

    // McKay & Piperno 2014 "Practical graph isomorphism, II". Impl: nauty-Traces-sys FFI.
    fn automorphisms_nauty<C: Ord + Copy>(&self, node_color: impl Fn(NodeId) -> C) -> Automorphism {
        let n = self.node_count();

        if n == 0 {
            return Automorphism {
                orbits: vec![],
                canonical_lab: vec![],
                node_count: 0,
                orbit_count: 0,
                group_order: AutoGroupOrder::Exact(1),
                generators: vec![],
            };
        }

        let mut indexed: Vec<(usize, C)> = self
            .node_ids()
            .map(|id| (id.index(), node_color(id)))
            .collect();
        indexed.sort_by_key(|&(_, c)| c);

        let mut lab = vec![0 as c_int; n];
        let mut ptn = vec![0 as c_int; n];
        for (pos, &(v, _)) in indexed.iter().enumerate() {
            lab[pos] = v as c_int;
        }
        for pos in 0..n.saturating_sub(1) {
            ptn[pos] = if indexed[pos].1 == indexed[pos + 1].1 {
                1
            } else {
                0
            };
        }

        let edge_count = self.edge_ids().count();
        let n_dir_edges = 2 * edge_count;
        let mut degree = vec![0usize; n];
        for eid in self.edge_ids() {
            let [a, b] = self.edge_endpoints(eid);
            degree[a.index()] += 1;
            degree[b.index()] += 1;
        }

        let mut sg = SparseGraph::new(n, n_dir_edges);
        let mut pos = 0usize;
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            sg.v[i] = pos;
            sg.d[i] = degree[i] as c_int;
            pos += degree[i];
        }

        let mut offset = vec![0usize; n];
        for eid in self.edge_ids() {
            let [a, b] = self.edge_endpoints(eid);
            let ai = a.index();
            let bi = b.index();
            sg.e[sg.v[ai] + offset[ai]] = bi as c_int;
            offset[ai] += 1;
            sg.e[sg.v[bi] + offset[bi]] = ai as c_int;
            offset[bi] += 1;
        }

        let mut orbits = vec![0 as c_int; n];
        let mut options = optionblk::default_sparse();
        options.getcanon = TRUE;
        options.defaultptn = FALSE;
        options.userautomproc = Some(capture_generator);
        let mut stats = statsblk::default();
        let mut cg = sparsegraph::default();

        GENERATORS.with(|g| g.borrow_mut().clear());
        let m = SETWORDSNEEDED(n);
        unsafe {
            nauty_check(
                WORDSIZE as c_int,
                m as c_int,
                n as c_int,
                NAUTYVERSIONID as c_int,
            );
            sparsenauty(
                &mut (&mut sg).into(),
                lab.as_mut_ptr(),
                ptn.as_mut_ptr(),
                orbits.as_mut_ptr(),
                &mut options,
                &mut stats,
                &mut cg,
            );
            SG_FREE(&mut cg);
        }

        let generators = GENERATORS.with(|g| mem::take(&mut *g.borrow_mut()));
        let orbits: Vec<NodeId> = orbits.iter().map(|&o| NodeId(o as u32)).collect();
        let canonical_lab: Vec<NodeId> = lab.iter().map(|&v| NodeId(v as u32)).collect();

        let orbit_count = {
            let mut reps = HashSet::new();
            for &o in &orbits {
                reps.insert(o);
            }
            reps.len()
        };

        let group_order = {
            let g1 = stats.grpsize1;
            let g2 = stats.grpsize2;
            if g2 == 0 && g1 >= 0.0 && g1 <= u32::MAX as f64 && g1.fract() == 0.0 {
                AutoGroupOrder::Exact(g1 as u32)
            } else if g2 == 0 {
                AutoGroupOrder::Approx(g1)
            } else {
                AutoGroupOrder::Approx(g1 * 10.0_f64.powi(g2))
            }
        };

        Automorphism {
            orbits,
            canonical_lab,
            node_count: n,
            orbit_count,
            group_order,
            generators,
        }
    }

    /// Numbering-invariant canonical key of `self` as an edge-colored graph: two
    /// graphs yield equal keys iff isomorphic under `node_color` and `edge_color`.
    /// Each edge is subdivided into a class-disjoint colored vertex so nauty (vertex
    /// colors only) canonicalizes edge colors; the key is the canonical node colors
    /// followed by the canonical edge list. Topology only — overlays/N-ary relations
    /// belong to a richer incidence built by the caller.
    pub fn canonical_key(
        &self,
        node_color: impl Fn(NodeId) -> Vec<u8>,
        edge_color: impl Fn(EdgeId) -> Vec<u8>,
        alg: AutomorphismAlgorithm,
    ) -> Vec<u8> {
        let nodes: Vec<NodeId> = self.node_ids().collect();
        let edges: Vec<EdgeId> = self.edge_ids().collect();
        let node_total = nodes.len();
        let subdivided_total = node_total + edges.len();

        let mut dense = vec![0u32; self.node_bound()];
        for (pos, &id) in nodes.iter().enumerate() {
            dense[id.index()] = pos as u32;
        }

        // Subdivide: original nodes 0..node_total, then one vertex per edge.
        let mut subdivided_edges: Vec<[u32; 2]> = Vec::with_capacity(2 * edges.len());
        for (k, &eid) in edges.iter().enumerate() {
            let [a, b] = self.edge_endpoints(eid);
            let edge_vertex = (node_total + k) as u32;
            subdivided_edges.push([dense[a.index()], edge_vertex]);
            subdivided_edges.push([edge_vertex, dense[b.index()]]);
        }
        let subdivided = Graph::new(subdivided_total, &subdivided_edges);

        // Class-prefixed colors keep node (0) and edge (1) vertices disjoint; ranks
        // feed nauty's vertex partition while the bytes go into the key.
        let mut colors: Vec<Vec<u8>> = Vec::with_capacity(subdivided_total);
        for &id in &nodes {
            let mut color = vec![0u8];
            color.extend_from_slice(&node_color(id));
            colors.push(color);
        }
        for &eid in &edges {
            let mut color = vec![1u8];
            color.extend_from_slice(&edge_color(eid));
            colors.push(color);
        }

        let mut distinct = colors.clone();
        distinct.sort();
        distinct.dedup();
        let ranks: Vec<u32> = colors
            .iter()
            .map(|color| distinct.binary_search(color).expect("color is present") as u32)
            .collect();

        let canonical = subdivided
            .automorphisms(|node| ranks[node.index()], alg)
            .canonical_labeling()
            .to_vec();
        let mut position = vec![0u32; subdivided_total];
        for (rank, node) in canonical.iter().enumerate() {
            position[node.index()] = rank as u32;
        }

        let mut key: Vec<u8> = Vec::new();
        key.extend_from_slice(&(subdivided_total as u32).to_le_bytes());
        for &node in &canonical {
            let color = &colors[node.index()];
            key.extend_from_slice(&(color.len() as u32).to_le_bytes());
            key.extend_from_slice(color);
        }
        let mut canonical_edges: Vec<(u32, u32)> = subdivided
            .edge_ids()
            .map(|edge| {
                let [u, v] = subdivided.edge_endpoints(edge);
                let (u, v) = (position[u.index()], position[v.index()]);
                (u.min(v), u.max(v))
            })
            .collect();
        canonical_edges.sort_unstable();
        key.extend_from_slice(&(canonical_edges.len() as u32).to_le_bytes());
        for (u, v) in canonical_edges {
            key.extend_from_slice(&u.to_le_bytes());
            key.extend_from_slice(&v.to_le_bytes());
        }
        key
    }
}

impl Automorphism {
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn orbit_count(&self) -> usize {
        self.orbit_count
    }

    pub fn orbit_of(&self, v: NodeId) -> NodeId {
        self.orbits[v.index()]
    }

    pub fn same_orbit(&self, a: NodeId, b: NodeId) -> bool {
        self.orbits[a.index()] == self.orbits[b.index()]
    }

    pub fn canonical_labeling(&self) -> &[NodeId] {
        &self.canonical_lab
    }

    pub fn auto_group_order(&self) -> AutoGroupOrder {
        self.group_order
    }

    /// A generating set of the automorphism group, each generator a permutation
    /// image over `0..node_count` (`generators()[k][i]` is the image of node `i`).
    /// Empty iff the group is trivial.
    pub fn generators(&self) -> &[Vec<NodeId>] {
        &self.generators
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::fmt::Debug;
    use std::thread;

    use rstest::*;

    use super::AutomorphismAlgorithm::Nauty;
    use super::*;
    use crate::union_find::UnionFind;

    fn assert_automorphism_semantics<C: Copy + Debug + Eq>(
        graph: &Graph,
        colors: &[C],
        aut: &Automorphism,
    ) {
        let expected_nodes: Vec<NodeId> = graph.node_ids().collect();
        let edges: HashSet<(usize, usize)> = graph
            .edge_ids()
            .map(|edge| {
                let [a, b] = graph.edge_endpoints(edge);
                (a.index().min(b.index()), a.index().max(b.index()))
            })
            .collect();
        let mut generated_orbits = UnionFind::new(graph.node_count());

        for generator in aut.generators() {
            let mut image = generator.clone();
            image.sort_unstable();
            assert_eq!(image, expected_nodes);
            for (source, &target) in generator.iter().enumerate() {
                assert_eq!(colors[source], colors[target.index()]);
                generated_orbits.union(source, target.index());
            }
            for &(a, b) in &edges {
                let a_image = generator[a].index();
                let b_image = generator[b].index();
                assert!(edges.contains(&(a_image.min(b_image), a_image.max(b_image))));
            }
        }

        for a in 0..graph.node_count() {
            for b in 0..graph.node_count() {
                assert_eq!(
                    generated_orbits.find(a) == generated_orbits.find(b),
                    aut.same_orbit(NodeId(a as u32), NodeId(b as u32))
                );
            }
        }
    }

    #[rstest]
    #[case::empty(Graph::default(), vec![], 0, AutoGroupOrder::Exact(1), vec![], vec![])]
    #[case::singleton(
        Graph::new(1, &[]),
        vec![0],
        1,
        AutoGroupOrder::Exact(1),
        vec![(0, 0)],
        vec![]
    )]
    #[case::same_color_edge(
        Graph::new(2, &[[0, 1]]),
        vec![0, 0],
        1,
        AutoGroupOrder::Exact(2),
        vec![(0, 1)],
        vec![]
    )]
    #[case::different_color_edge(
        Graph::new(2, &[[0, 1]]),
        vec![0, 1],
        2,
        AutoGroupOrder::Exact(1),
        vec![],
        vec![(0, 1)]
    )]
    #[case::uniform_square(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        vec![0, 0, 0, 0],
        1,
        AutoGroupOrder::Exact(8),
        vec![(0, 1), (0, 2), (0, 3)],
        vec![]
    )]
    #[case::colored_path(
        Graph::new(3, &[[0, 1], [1, 2]]),
        vec![0, 1, 0],
        2,
        AutoGroupOrder::Exact(2),
        vec![(0, 2)],
        vec![(0, 1)]
    )]
    #[case::disconnected_edges(
        Graph::new(4, &[[0, 1], [2, 3]]),
        vec![0, 0, 0, 0],
        1,
        AutoGroupOrder::Exact(8),
        vec![(0, 1), (0, 2), (0, 3)],
        vec![]
    )]
    fn test_graph_automorphisms(
        #[case] graph: Graph,
        #[case] colors: Vec<u8>,
        #[case] expected_orbits: usize,
        #[case] expected_order: AutoGroupOrder,
        #[case] same_orbit: Vec<(u32, u32)>,
        #[case] different_orbit: Vec<(u32, u32)>,
    ) {
        let aut = graph.automorphisms(|node| colors[node.index()], Nauty);
        assert_eq!(aut.node_count(), graph.node_count());
        assert_eq!(aut.orbit_count(), expected_orbits);
        assert_eq!(aut.auto_group_order(), expected_order);
        for (a, b) in same_orbit {
            assert!(aut.same_orbit(NodeId(a), NodeId(b)));
        }
        for (a, b) in different_orbit {
            assert!(!aut.same_orbit(NodeId(a), NodeId(b)));
        }
        assert_automorphism_semantics(&graph, &colors, &aut);
    }

    #[rstest]
    #[case::swap(Graph::new(2, &[[0, 1]]), vec![0u8, 0], vec![vec![NodeId(1), NodeId(0)]])]
    #[case::colored_path(Graph::new(3, &[[0, 1], [1, 2]]), vec![0u8, 1, 0], vec![vec![NodeId(2), NodeId(1), NodeId(0)]])]
    #[case::trivial(Graph::new(2, &[[0, 1]]), vec![0u8, 1], vec![])]
    fn test_automorphism_generators(
        #[case] g: Graph,
        #[case] colors: Vec<u8>,
        #[case] expected: Vec<Vec<NodeId>>,
    ) {
        let aut = g.automorphisms(|n| colors[n.index()], Nauty);
        assert_eq!(aut.generators(), expected);
    }

    #[rstest]
    #[case::complete_13(13, 6_227_020_800.0)]
    fn test_graph_automorphisms_large_order(#[case] nodes: usize, #[case] expected: f64) {
        let mut edges = Vec::new();
        for a in 0..nodes as u32 {
            for b in a + 1..nodes as u32 {
                edges.push([a, b]);
            }
        }
        let graph = Graph::new(nodes, &edges);
        let aut = graph.automorphisms(|_| 0u8, Nauty);
        let AutoGroupOrder::Approx(order) = aut.auto_group_order() else {
            panic!("order above u32 must use the approximate representation");
        };
        assert!((order - expected).abs() < 0.5);
        assert_automorphism_semantics(&graph, &vec![0; nodes], &aut);
    }

    #[rstest]
    #[case::cycle_6(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        NodeId(2),
        AutoGroupOrder::Exact(2)
    )]
    fn test_graph_automorphisms_stabilizer(
        #[case] graph: Graph,
        #[case] site: NodeId,
        #[case] expected_order: AutoGroupOrder,
    ) {
        let colors: Vec<bool> = graph.node_ids().map(|node| node == site).collect();
        let aut = graph.automorphisms(|node| colors[node.index()], Nauty);
        assert_eq!(aut.auto_group_order(), expected_order);
        assert!(aut
            .generators()
            .iter()
            .all(|generator| generator[site.index()] == site));
    }

    #[rstest]
    #[case::parallel_cycles(8)]
    fn test_graph_automorphisms_concurrency(#[case] thread_count: usize) {
        let handles: Vec<_> = (0..thread_count)
            .map(|site| {
                thread::spawn(move || {
                    let graph = Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]);
                    let site = NodeId((site % graph.node_count()) as u32);
                    let colors: Vec<bool> = graph.node_ids().map(|node| node == site).collect();
                    let aut = graph.automorphisms(|node| colors[node.index()], Nauty);
                    assert_automorphism_semantics(&graph, &colors, &aut);
                    (site, aut)
                })
            })
            .collect();

        for handle in handles {
            let (site, aut) = handle.join().expect("automorphism worker succeeds");
            assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(2));
            assert!(aut
                .generators()
                .iter()
                .all(|generator| generator[site.index()] == site));
        }
    }

    #[rstest]
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
    fn test_graph_canonical_key_isomorphic(
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
            Nauty,
        );
        let key_b = b.canonical_key(
            |node| b_nodes[node.index()].clone(),
            |edge| b_edges[edge.index()].clone(),
            Nauty,
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
    fn test_graph_canonical_key_distinct(
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
            Nauty,
        );
        let key_b = b.canonical_key(
            |node| b_nodes[node.index()].clone(),
            |edge| b_edges[edge.index()].clone(),
            Nauty,
        );
        assert_ne!(key_a, key_b);
    }
}
