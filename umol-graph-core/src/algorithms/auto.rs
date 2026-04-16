//! Graph automorphism and canonical labeling via nauty.

use std::collections::HashSet;
use std::os::raw::c_int;

use nauty_Traces_sys::*;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AutoGroupOrder {
    Exact(u32),
    Approx(f64),
}

#[derive(Debug, Clone)]
pub struct Automorphism {
    orbits: Vec<c_int>,
    canonical_lab: Vec<c_int>,
    n: usize,
    num_orbits: usize,
    grpsize1: f64,
    grpsize2: c_int,
}

impl Automorphism {
    /// Compute automorphism orbits and canonical labeling for a colored graph.
    ///
    /// `n` is the total vertex count. `edges` lists undirected edges as `(u, v)`
    /// pairs with `u, v < n`. `colors` assigns an integer color to each vertex;
    /// vertices with different colors are never in the same orbit.
    pub fn compute<C: Ord + Copy>(
        n: usize,
        edges: &[(usize, usize)],
        colors: &[C],
    ) -> Self {
        assert_eq!(colors.len(), n);

        if n == 0 {
            return Self {
                orbits: vec![],
                canonical_lab: vec![],
                n: 0,
                num_orbits: 0,
                grpsize1: 1.0,
                grpsize2: 0,
            };
        }

        let mut indexed: Vec<(usize, C)> = colors.iter().copied().enumerate().collect();
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

        let n_dir_edges = 2 * edges.len();
        let mut degree = vec![0usize; n];
        for &(a, b) in edges {
            degree[a] += 1;
            degree[b] += 1;
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
        for &(a, b) in edges {
            sg.e[sg.v[a] + offset[a]] = b as c_int;
            offset[a] += 1;
            sg.e[sg.v[b] + offset[b]] = a as c_int;
            offset[b] += 1;
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

        let num_orbits = {
            let mut reps = HashSet::new();
            for &o in &orbits {
                reps.insert(o);
            }
            reps.len()
        };

        Self {
            orbits,
            canonical_lab: lab,
            n,
            num_orbits,
            grpsize1: stats.grpsize1,
            grpsize2: stats.grpsize2,
        }
    }

    pub fn vertex_count(&self) -> usize {
        self.n
    }

    pub fn num_orbits(&self) -> usize {
        self.num_orbits
    }

    pub fn orbit_of(&self, v: usize) -> usize {
        self.orbits[v] as usize
    }

    pub fn same_orbit(&self, a: usize, b: usize) -> bool {
        self.orbits[a] == self.orbits[b]
    }

    pub fn canonical_labeling(&self) -> &[c_int] {
        &self.canonical_lab
    }

    pub fn auto_group_order(&self) -> AutoGroupOrder {
        if self.grpsize2 == 0 {
            let g = self.grpsize1;
            if g >= 0.0 && g <= u32::MAX as f64 && g.fract() == 0.0 {
                return AutoGroupOrder::Exact(g as u32);
            }
        }
        let approx = if self.grpsize2 == 0 {
            self.grpsize1
        } else {
            self.grpsize1 * 10.0_f64.powi(self.grpsize2)
        };
        AutoGroupOrder::Approx(approx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_automorphism_empty() {
        let aut = Automorphism::compute::<u8>(0, &[], &[]);
        assert_eq!(aut.num_orbits(), 0);
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphism_single_vertex() {
        let aut = Automorphism::compute(1, &[], &[0u8]);
        assert_eq!(aut.num_orbits(), 1);
        assert_eq!(aut.orbit_of(0), 0);
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphism_two_same_color() {
        let aut = Automorphism::compute(2, &[(0, 1)], &[0u8, 0]);
        assert_eq!(aut.num_orbits(), 1);
        assert!(aut.same_orbit(0, 1));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(2));
    }

    #[test]
    fn test_automorphism_two_different_color() {
        let aut = Automorphism::compute(2, &[(0, 1)], &[0u8, 1]);
        assert_eq!(aut.num_orbits(), 2);
        assert!(!aut.same_orbit(0, 1));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(1));
    }

    #[test]
    fn test_automorphism_square_uniform() {
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];
        let aut = Automorphism::compute(4, &edges, &[0u8, 0, 0, 0]);
        assert_eq!(aut.num_orbits(), 1);
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(8));
    }

    #[test]
    fn test_automorphism_path_colored() {
        let edges = vec![(0, 1), (1, 2)];
        let aut = Automorphism::compute(3, &edges, &[0u8, 1, 0]);
        assert_eq!(aut.num_orbits(), 2);
        assert!(aut.same_orbit(0, 2));
        assert!(!aut.same_orbit(0, 1));
        assert_eq!(aut.auto_group_order(), AutoGroupOrder::Exact(2));
    }
}
