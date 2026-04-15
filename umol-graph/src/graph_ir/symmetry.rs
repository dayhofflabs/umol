//! Graph automorphism and canonical labeling via nauty.

use std::collections::{BTreeMap, HashSet};
use std::os::raw::c_int;

use index_vec::Idx;
use nauty_Traces_sys::*;
use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin::SpinMultiplicity;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::AtomIdx;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::molecule::MoleculeAst;

/// Vertex color for nauty partitioning.
/// Atom and Bond variants are in separate cells, so they can never share an orbit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum VertexColor {
    Atom {
        element: Element,
        isotope_mass: Option<u32>,
        charge: Option<i8>,
        hydrogens: Option<u8>,
        unpaired_electrons: Option<u8>,
        multiplicity: Option<SpinMultiplicity>,
    },
    Bond {
        order: u8,
        donation: u8,
    },
}

fn atom_color(atom: &AtomAst) -> VertexColor {
    VertexColor::Atom {
        element: match atom.element {
            ElementAst::Lit(e) => e,
            _ => Element::Og,
        },
        isotope_mass: match &atom.isotope_mass {
            IsotopeAst::Lit(n) => Some(*n),
            _ => None,
        },
        charge: match &atom.charge {
            ValueAst::Lit(n) => Some(*n as i8),
            _ => None,
        },
        hydrogens: match &atom.implicit_hydrogens {
            HydrogenAst::Value(ValueAst::Lit(n)) => Some(*n as u8),
            _ => None,
        },
        unpaired_electrons: match &atom.spin {
            SpinStateAst::Lit(s) => Some(s.unpaired_electrons()),
            _ => None,
        },
        multiplicity: match &atom.spin {
            SpinStateAst::Lit(s) => Some(s.multiplicity()),
            _ => None,
        },
    }
}

fn bond_color(bond: &BondAst) -> VertexColor {
    VertexColor::Bond {
        order: match &bond.order {
            ValueAst::Lit(n) => *n as u8,
            _ => 0,
        },
        donation: 0,
    }
}

/// Automorphism group order: exact when it fits in `u32`, else approximate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AutoGroupOrder {
    Exact(u32),
    Approx(f64),
}

/// Result of graph automorphism computation on a molecular graph.
#[derive(Debug, Clone)]
pub struct GraphSymmetry {
    nauty_to_atom: Vec<AtomIdx>,
    atom_to_nauty: Vec<usize>,
    orbits: Vec<c_int>,
    canonical_lab: Vec<c_int>,
    n_atoms: usize,
    num_orbits: usize,
    grpsize1: f64,
    grpsize2: c_int,
}

impl GraphSymmetry {
    pub fn num_orbits(&self) -> usize {
        self.num_orbits
    }

    pub fn orbit_representative(&self, atom: AtomIdx) -> AtomIdx {
        let ni = self.atom_to_nauty[atom.index()];
        debug_assert!(ni < self.n_atoms);
        debug_assert!((self.orbits[ni] as usize) < self.n_atoms);
        self.nauty_to_atom[self.orbits[ni] as usize]
    }

    pub fn same_orbit(&self, a: AtomIdx, b: AtomIdx) -> bool {
        debug_assert!(a.index() < self.n_atoms);
        debug_assert!(b.index() < self.n_atoms);
        debug_assert!(self.atom_to_nauty[a.index()] < self.n_atoms);
        debug_assert!(self.atom_to_nauty[b.index()] < self.n_atoms);
        self.orbits[self.atom_to_nauty[a.index()]] == self.orbits[self.atom_to_nauty[b.index()]]
    }

    pub fn orbit_partition(&self) -> Vec<Vec<AtomIdx>> {
        let mut groups: BTreeMap<c_int, Vec<AtomIdx>> = BTreeMap::new();
        for i in 0..self.n_atoms {
            groups
                .entry(self.orbits[i])
                .or_default()
                .push(self.nauty_to_atom[i]);
        }
        groups.into_values().collect()
    }

    pub fn canonical_order(&self) -> Vec<AtomIdx> {
        self.canonical_lab
            .iter()
            .filter_map(|&v| {
                let v = v as usize;
                (v < self.n_atoms).then(|| self.nauty_to_atom[v])
            })
            .collect()
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

    pub fn auto_group_order_exact(&self) -> Option<u32> {
        match self.auto_group_order() {
            AutoGroupOrder::Exact(n) => Some(n),
            AutoGroupOrder::Approx(_) => None,
        }
    }
}

/// Compute graph automorphism orbits and canonical labeling for a molecule.
///
/// Bond types are encoded via edge subdivision: each bond becomes an auxiliary
/// vertex (colored by bond order and donation) connecting its two endpoints.
pub fn compute_symmetry(ast: &MoleculeAst) -> GraphSymmetry {
    let n_atoms = ast.atom_count();

    if n_atoms == 0 {
        return GraphSymmetry {
            nauty_to_atom: vec![],
            atom_to_nauty: vec![],
            orbits: vec![],
            canonical_lab: vec![],
            n_atoms: 0,
            num_orbits: 0,
            grpsize1: 1.0,
            grpsize2: 0,
        };
    }

    let mut atom_to_nauty = vec![usize::MAX; n_atoms];
    let mut nauty_to_atom = Vec::with_capacity(n_atoms);
    for (ni, (idx, _)) in ast.atoms().enumerate() {
        atom_to_nauty[idx.index()] = ni;
        nauty_to_atom.push(idx);
    }

    let bonds: Vec<_> = ast.bonds().collect();
    let n_bonds = bonds.len();
    let n_total = n_atoms + n_bonds;

    let mut colored: Vec<(usize, VertexColor)> = Vec::with_capacity(n_total);
    for (ni, (_, atom)) in ast.atoms().enumerate() {
        colored.push((ni, atom_color(atom)));
    }

    let mut bond_endpoints: Vec<(usize, usize, usize)> = Vec::with_capacity(n_bonds);
    for (i, (_, src, tgt, bond)) in bonds.iter().enumerate() {
        let aux = n_atoms + i;
        colored.push((aux, bond_color(bond)));
        bond_endpoints.push((atom_to_nauty[src.index()], atom_to_nauty[tgt.index()], aux));
    }

    colored.sort_by_key(|&(_, c)| c);
    let mut lab = vec![0 as c_int; n_total];
    let mut ptn = vec![0 as c_int; n_total];
    for (pos, &(v, _)) in colored.iter().enumerate() {
        lab[pos] = v as c_int;
    }
    for pos in 0..n_total.saturating_sub(1) {
        ptn[pos] = if colored[pos].1 == colored[pos + 1].1 {
            1
        } else {
            0
        };
    }

    let n_dir_edges = 4 * n_bonds;
    let mut degree = vec![0usize; n_total];
    for &(a, b, aux) in &bond_endpoints {
        degree[a] += 1;
        degree[b] += 1;
        degree[aux] = 2;
    }

    let mut sg = SparseGraph::new(n_total, n_dir_edges);
    let mut pos = 0usize;
    #[allow(clippy::needless_range_loop)]
    for i in 0..n_total {
        sg.v[i] = pos;
        sg.d[i] = degree[i] as c_int;
        pos += degree[i];
    }

    let mut offset = vec![0usize; n_total];
    for &(a, b, aux) in &bond_endpoints {
        sg.e[sg.v[a] + offset[a]] = aux as c_int;
        offset[a] += 1;
        sg.e[sg.v[aux] + offset[aux]] = a as c_int;
        offset[aux] += 1;
        sg.e[sg.v[b] + offset[b]] = aux as c_int;
        offset[b] += 1;
        sg.e[sg.v[aux] + offset[aux]] = b as c_int;
        offset[aux] += 1;
    }

    let mut orbits = vec![0 as c_int; n_total];
    let mut options = optionblk::default_sparse();
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    let mut stats = statsblk::default();
    let mut cg = sparsegraph::default();

    let m = SETWORDSNEEDED(n_total);
    unsafe {
        nauty_check(
            WORDSIZE as c_int,
            m as c_int,
            n_total as c_int,
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
        for &orbit in &orbits[..n_atoms] {
            reps.insert(orbit);
        }
        reps.len()
    };

    GraphSymmetry {
        nauty_to_atom,
        atom_to_nauty,
        orbits,
        canonical_lab: lab,
        n_atoms,
        num_orbits,
        grpsize1: stats.grpsize1,
        grpsize2: stats.grpsize2,
    }
}

#[cfg(test)]
mod tests {
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::AtomIdx;
    use crate::ast::bond::BondAst;

    fn atom(element: Element) -> AtomAst {
        AtomAst::from_element(element)
    }

    fn mol(atoms: Vec<AtomAst>, edges: &[(usize, usize, u8)]) -> MoleculeAst {
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b, order)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(order)))
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    #[test]
    fn test_compute_symmetry_empty() {
        let ast = MoleculeAst::default();
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 0);
        assert_eq!(sym.orbit_partition().len(), 0);
        assert_eq!(sym.canonical_order().len(), 0);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(1)));
    }

    #[test]
    fn test_compute_symmetry_single_atom() {
        let ast = mol(vec![atom(Element::C)], &[]);
        let a = AtomIdx(0);
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 1);
        assert_eq!(sym.orbit_representative(a), a);
        assert_eq!(sym.orbit_partition(), vec![vec![a]]);
        assert_eq!(sym.canonical_order(), vec![a]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
    }

    #[test]
    fn test_compute_symmetry_h2_equivalent() {
        let ast = mol(vec![atom(Element::H), atom(Element::H)], &[(0, 1, 1)]);
        let h1 = AtomIdx(0);
        let h2 = AtomIdx(1);
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(h1, h2));
        assert_eq!(sym.orbit_representative(h1), h1);
        assert_eq!(sym.orbit_representative(h2), h1);
        assert_eq!(sym.orbit_partition(), vec![vec![h1, h2]]);
        assert_eq!(sym.canonical_order(), vec![h1, h2]);
        assert_eq!(sym.auto_group_order_exact(), Some(2));
    }

    #[test]
    fn test_compute_symmetry_hf_distinct() {
        let ast = mol(vec![atom(Element::H), atom(Element::F)], &[(0, 1, 1)]);
        let h = AtomIdx(0);
        let f = AtomIdx(1);
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 2);
        assert!(!sym.same_orbit(h, f));
        assert_eq!(sym.orbit_representative(h), h);
        assert_eq!(sym.orbit_representative(f), f);
        assert_eq!(sym.orbit_partition(), vec![vec![h], vec![f]]);
        assert_eq!(sym.canonical_order(), vec![h, f]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
    }

    #[test]
    fn test_compute_symmetry_square_uniform_bonds() {
        let ast = mol(
            vec![atom(Element::C); 4],
            &[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1)],
        );
        let c: Vec<_> = (0..4).map(|i| AtomIdx(i)).collect();
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(c[0], c[1]));
        assert!(sym.same_orbit(c[1], c[2]));
        assert!(sym.same_orbit(c[2], c[3]));
        assert_eq!(sym.orbit_partition(), vec![vec![c[0], c[1], c[2], c[3]]]);
        assert_eq!(sym.canonical_order(), vec![c[0], c[1], c[3], c[2]]);
        assert_eq!(sym.auto_group_order_exact(), Some(8));
    }

    #[test]
    fn test_compute_symmetry_linear_mixed_bond_orders() {
        let ast = mol(
            vec![atom(Element::C); 3],
            &[(0, 1, 2), (1, 2, 1)],
        );
        let c: Vec<_> = (0..3).map(|i| AtomIdx(i)).collect();
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 3);
        assert!(!sym.same_orbit(c[0], c[1]));
        assert!(!sym.same_orbit(c[1], c[2]));
        assert!(!sym.same_orbit(c[2], c[0]));
        assert_eq!(sym.canonical_order(), vec![c[0], c[2], c[1]]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
    }

    #[test]
    fn test_compute_symmetry_alternating_cycle() {
        let ast = mol(
            vec![atom(Element::C); 4],
            &[(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 0, 1)],
        );
        let c: Vec<_> = (0..4).map(|i| AtomIdx(i)).collect();
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(c[0], c[1]));
        assert!(sym.same_orbit(c[1], c[2]));
        assert!(sym.same_orbit(c[2], c[3]));
        assert_eq!(sym.canonical_order(), vec![c[0], c[3], c[1], c[2]]);
        assert_eq!(sym.auto_group_order_exact(), Some(4));
    }

    #[test]
    fn test_compute_symmetry_water() {
        let ast = mol(
            vec![atom(Element::H), atom(Element::O), atom(Element::H)],
            &[(0, 1, 1), (1, 2, 1)],
        );
        let h1 = AtomIdx(0);
        let o = AtomIdx(1);
        let h2 = AtomIdx(2);
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 2);
        assert!(sym.same_orbit(h1, h2));
        assert!(!sym.same_orbit(h1, o));
        assert_eq!(sym.orbit_partition(), vec![vec![h1, h2], vec![o]]);
        assert_eq!(sym.canonical_order(), vec![h1, h2, o]);
        assert_eq!(sym.auto_group_order_exact(), Some(2));
    }

    #[test]
    fn test_compute_symmetry_benzene_ring() {
        let ast = mol(
            vec![atom(Element::C); 6],
            &[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)],
        );
        let c: Vec<_> = (0..6).map(|i| AtomIdx(i)).collect();
        let sym = compute_symmetry(&ast);
        assert_eq!(sym.num_orbits(), 1);
        for i in 0..6 {
            assert!(sym.same_orbit(c[0], c[i]));
        }
        assert_eq!(sym.orbit_partition(), vec![vec![c[0], c[1], c[2], c[3], c[4], c[5]]]);
        assert_eq!(sym.canonical_order(), vec![c[0], c[2], c[4], c[3], c[1], c[5]]);
        assert_eq!(sym.auto_group_order_exact(), Some(12));
    }

    #[test]
    fn test_compute_symmetry_canonical_order_deterministic() {
        let ast = mol(
            vec![atom(Element::C), atom(Element::N), atom(Element::O)],
            &[],
        );
        let order1 = compute_symmetry(&ast).canonical_order();
        let order2 = compute_symmetry(&ast).canonical_order();
        assert_eq!(order1, order2);
        assert_eq!(order1.len(), 3);
    }
}
