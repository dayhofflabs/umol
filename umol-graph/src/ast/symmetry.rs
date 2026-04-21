//! Graph automorphism and canonical labeling on molecular ASTs.
//!
//! This is a thin adapter over the generic automorphism engine in
//! `umol_graph_core::algorithms::auto`. It encodes atom attributes as vertex
//! colors, represents bond types via edge subdivision (each bond becomes an
//! auxiliary vertex colored by its type), and wraps results with `AtomIdx`
//! handles.

use std::collections::{BTreeMap, HashSet};

use umol_graph_core::algorithms::auto::Automorphism;
use umol_ast::ast::atom::{ElementAst, ImplicitHydrogensAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin::SpinMultiplicity;
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

use super::AtomIdx;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::molecule::MoleculeAst;

pub use umol_graph_core::algorithms::auto::AutoGroupOrder;

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
            ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => Some(*n as u8),
            _ => None,
        },
        unpaired_electrons: atom
            .spin
            .try_into_ground()
            .ok()
            .flatten()
            .map(|s| s.unpaired_electrons()),
        multiplicity: atom
            .spin
            .try_into_ground()
            .ok()
            .flatten()
            .map(|s| s.multiplicity()),
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

/// Result of graph automorphism computation on a molecular graph.
#[derive(Debug, Clone)]
pub struct GraphSymmetry {
    nauty_to_atom: Vec<AtomIdx>,
    atom_to_nauty: Vec<usize>,
    n_atoms: usize,
    aut: Automorphism,
}

impl GraphSymmetry {
    pub fn num_orbits(&self) -> usize {
        let mut reps = HashSet::new();
        for i in 0..self.n_atoms {
            reps.insert(self.aut.orbit_of(i));
        }
        reps.len()
    }

    pub fn orbit_representative(&self, atom: AtomIdx) -> AtomIdx {
        let ni = self.atom_to_nauty[atom.index()];
        self.nauty_to_atom[self.aut.orbit_of(ni)]
    }

    pub fn same_orbit(&self, a: AtomIdx, b: AtomIdx) -> bool {
        self.aut
            .same_orbit(self.atom_to_nauty[a.index()], self.atom_to_nauty[b.index()])
    }

    pub fn orbit_partition(&self) -> Vec<Vec<AtomIdx>> {
        let mut groups: BTreeMap<usize, Vec<AtomIdx>> = BTreeMap::new();
        for i in 0..self.n_atoms {
            groups
                .entry(self.aut.orbit_of(i))
                .or_default()
                .push(self.nauty_to_atom[i]);
        }
        groups.into_values().collect()
    }

    pub fn canonical_order(&self) -> Vec<AtomIdx> {
        self.aut
            .canonical_labeling()
            .iter()
            .filter_map(|&v| {
                let v = v as usize;
                (v < self.n_atoms).then(|| self.nauty_to_atom[v])
            })
            .collect()
    }

    pub fn auto_group_order(&self) -> AutoGroupOrder {
        self.aut.auto_group_order()
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
    let n_atoms = ast.atoms().count();

    let mut atom_to_nauty = vec![usize::MAX; n_atoms];
    let mut nauty_to_atom = Vec::with_capacity(n_atoms);
    for (ni, view) in ast.atoms().iter().enumerate() {
        atom_to_nauty[view.idx.index()] = ni;
        nauty_to_atom.push(view.idx);
    }

    let bonds: Vec<_> = ast.bonds().iter().collect();
    let n_bonds = bonds.len();
    let n_total = n_atoms + n_bonds;

    let mut colors: Vec<VertexColor> = Vec::with_capacity(n_total);
    for view in ast.atoms().iter() {
        colors.push(atom_color(view.data));
    }
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(2 * n_bonds);
    for (i, b) in bonds.iter().enumerate() {
        let aux = n_atoms + i;
        colors.push(bond_color(b.data));
        edges.push((atom_to_nauty[b.src.index()], aux));
        edges.push((atom_to_nauty[b.tgt.index()], aux));
    }

    let aut = Automorphism::compute(n_total, &edges, &colors);

    GraphSymmetry {
        nauty_to_atom,
        atom_to_nauty,
        n_atoms,
        aut,
    }
}

#[cfg(test)]
mod tests {
    use umol_shared::element::Element;

    use super::*;
    use super::AtomIdx;
    use super::super::bond::BondAst;

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
        let c: Vec<_> = (0..4).map(AtomIdx).collect();
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
        let c: Vec<_> = (0..3).map(AtomIdx).collect();
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
        let c: Vec<_> = (0..4).map(AtomIdx).collect();
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
        let c: Vec<_> = (0..6).map(AtomIdx).collect();
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
