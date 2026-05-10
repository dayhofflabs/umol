//! Graph automorphism and canonical labeling.

use std::collections::HashSet;
use std::os::raw::c_int;

use nauty_Traces_sys::*;

use crate::graph::{Graph, NodeId};

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
        let mut stats = statsblk::default();
        let mut cg = sparsegraph::default();

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

        let orbits: Vec<NodeId> = orbits.iter().map(|&o| NodeId(o as u32)).collect();
        let canonical_lab: Vec<NodeId> = lab.iter().map(|&v| NodeId(v as u32)).collect();

        let num_orbits = {
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
            orbit_count: num_orbits,
            group_order,
        }
    }
}

impl Automorphism {
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn num_orbits(&self) -> usize {
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
}

#[cfg(test)]
mod tests {
    use super::AutomorphismAlgorithm::Nauty;
    use super::*;

    #[test]
    fn test_automorphisms_empty() {
        let g = Graph::default();
        let aut = g.automorphisms(|_: NodeId| 0u8, Nauty);
        assert_eq!(aut.num_orbits(), 0);
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphisms_single_vertex() {
        let g = Graph::new(1, &[]);
        let aut = g.automorphisms(|_| 0u8, Nauty);
        assert_eq!(aut.num_orbits(), 1);
        assert_eq!(aut.orbit_of(NodeId(0)), NodeId(0));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphisms_two_same_color() {
        let g = Graph::new(2, &[[0, 1]]);
        let aut = g.automorphisms(|_| 0u8, Nauty);
        assert_eq!(aut.num_orbits(), 1);
        assert!(aut.same_orbit(NodeId(0), NodeId(1)));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(2));
    }

    #[test]
    fn test_automorphisms_two_different_color() {
        let g = Graph::new(2, &[[0, 1]]);
        let aut = g.automorphisms(|n| n.index() as u8, Nauty);
        assert_eq!(aut.num_orbits(), 2);
        assert!(!aut.same_orbit(NodeId(0), NodeId(1)));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphisms_square_uniform() {
        let g = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let aut = g.automorphisms(|_| 0u8, Nauty);
        assert_eq!(aut.num_orbits(), 1);
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(8));
    }

    #[test]
    fn test_automorphisms_path_colored() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let colors = [0u8, 1, 0];
        let aut = g.automorphisms(|n| colors[n.index()], Nauty);
        assert_eq!(aut.num_orbits(), 2);
        assert!(aut.same_orbit(NodeId(0), NodeId(2)));
        assert!(!aut.same_orbit(NodeId(0), NodeId(1)));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(2));
    }
}
