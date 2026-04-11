//! Graph automorphism and canonical labeling via nauty.

use std::collections::{BTreeMap, HashSet};
use std::os::raw::c_int;

use nauty_Traces_sys::*;
use umol_shared::element::Element;
use umol_shared::spin::SpinMultiplicity;

use super::atom_pattern::AtomPattern;
use super::bond_pattern::BondPattern;
use super::molecule::AtomIndex;
use super::molecule_builder::MoleculeBuilder;

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

fn atom_color(ab: &AtomPattern) -> VertexColor {
    use crate::graph_ir::atom_pattern::{HydrogenPattern, IsotopePattern};
    VertexColor::Atom {
        element: ab.element(),
        isotope_mass: match ab.isotope_mass {
            IsotopePattern::Is(n) => Some(n),
            _ => None,
        },
        charge: ab.charge.into_option(),
        hydrogens: match ab.implicit_hydrogens {
            HydrogenPattern::Is(h) => Some(h),
            _ => None,
        },
        unpaired_electrons: ab.unpaired_electrons.into_option(),
        multiplicity: ab.multiplicity.into_option(),
    }
}

fn bond_color(b: &BondPattern) -> VertexColor {
    VertexColor::Bond {
        order: b.order(),
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
    nauty_to_atom: Vec<AtomIndex>,
    atom_to_nauty: Vec<usize>,
    orbits: Vec<c_int>,
    canonical_lab: Vec<c_int>,
    n_atoms: usize,
    num_orbits: usize,
    grpsize1: f64,
    grpsize2: c_int,
}

impl GraphSymmetry {
    /// Number of distinct orbits among atoms.
    pub fn num_orbits(&self) -> usize {
        self.num_orbits
    }

    /// Canonical orbit representative for an atom.
    pub fn orbit_representative(&self, atom: AtomIndex) -> AtomIndex {
        let ni = self.atom_to_nauty[atom.index()];
        debug_assert!(ni < self.n_atoms);
        debug_assert!((self.orbits[ni] as usize) < self.n_atoms);
        self.nauty_to_atom[self.orbits[ni] as usize]
    }

    /// Whether two atoms belong to the same orbit.
    pub fn same_orbit(&self, a: AtomIndex, b: AtomIndex) -> bool {
        debug_assert!(a.index() < self.n_atoms);
        debug_assert!(b.index() < self.n_atoms);
        debug_assert!(self.atom_to_nauty[a.index()] < self.n_atoms);
        debug_assert!(self.atom_to_nauty[b.index()] < self.n_atoms);
        self.orbits[self.atom_to_nauty[a.index()]] == self.orbits[self.atom_to_nauty[b.index()]]
    }

    /// Atoms grouped by orbit.
    pub fn orbit_partition(&self) -> Vec<Vec<AtomIndex>> {
        let mut groups: BTreeMap<c_int, Vec<AtomIndex>> = BTreeMap::new();
        for i in 0..self.n_atoms {
            groups
                .entry(self.orbits[i])
                .or_default()
                .push(self.nauty_to_atom[i]);
        }
        groups.into_values().collect()
    }

    /// Atoms in canonical order (deterministic across identical graphs).
    pub fn canonical_order(&self) -> Vec<AtomIndex> {
        self.canonical_lab
            .iter()
            .filter_map(|&v| {
                let v = v as usize;
                (v < self.n_atoms).then(|| self.nauty_to_atom[v])
            })
            .collect()
    }

    /// Automorphism group order: exact if it fits in `u32`, else approximate.
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

    /// Exact automorphism group order when representable as `u32`.
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
pub fn compute_symmetry(builder: &MoleculeBuilder) -> GraphSymmetry {
    let atom_indices: Vec<AtomIndex> = builder.atom_indices().collect();
    let n_atoms = atom_indices.len();

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

    // AtomIndex <-> dense nauty index mapping
    let max_atom_idx = atom_indices.iter().map(|a| a.index()).max().unwrap();
    let mut atom_to_nauty = vec![usize::MAX; max_atom_idx + 1];
    let mut nauty_to_atom = Vec::with_capacity(n_atoms);
    for (ni, &ai) in atom_indices.iter().enumerate() {
        atom_to_nauty[ai.index()] = ni;
        nauty_to_atom.push(ai);
    }

    let bond_indices: Vec<_> = builder.bond_indices().collect();
    let n_bonds = bond_indices.len();
    let n_total = n_atoms + n_bonds;

    // Vertex colors and bond endpoint info
    let mut colored: Vec<(usize, VertexColor)> = Vec::with_capacity(n_total);
    for (ni, &ai) in atom_indices.iter().enumerate() {
        colored.push((ni, atom_color(builder.atom(ai).unwrap())));
    }

    let mut bond_endpoints: Vec<(usize, usize, usize)> = Vec::with_capacity(n_bonds);
    for (i, &bi) in bond_indices.iter().enumerate() {
        let aux = n_atoms + i;
        let (a, b) = builder.bond_atom_indices(bi).unwrap();
        colored.push((aux, bond_color(builder.bond(bi).unwrap())));
        bond_endpoints.push((atom_to_nauty[a.index()], atom_to_nauty[b.index()], aux));
    }

    // Build partition (lab/ptn) from vertex colors
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

    // Build nauty sparse graph (each bond → auxiliary vertex with 2 edges)
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

    // Run nauty
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
    use crate::graph_ir::atom_pattern::AtomPattern;

    #[test]
    fn empty() {
        let b = MoleculeBuilder::new();
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 0);
        assert_eq!(sym.orbit_partition().len(), 0);
        assert_eq!(sym.canonical_order().len(), 0);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(1)));
    }

    #[test]
    fn single_atom() {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(AtomPattern::new(Element::C));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert_eq!(sym.orbit_representative(a), a);
        assert_eq!(sym.orbit_partition(), vec![vec![a]]);
        assert_eq!(sym.canonical_order(), vec![a]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(1)));
    }

    #[test]
    fn h2_equivalent() {
        let mut b = MoleculeBuilder::new();
        let h1 = b.add_atom(AtomPattern::new(Element::H));
        let h2 = b.add_atom(AtomPattern::new(Element::H));
        b.add_bond(h1, h2, BondPattern::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(h1, h2));
        assert_eq!(sym.orbit_representative(h1), h1);
        assert_eq!(sym.orbit_representative(h2), h1);
        assert_eq!(sym.orbit_partition(), vec![vec![h1, h2]]);
        assert_eq!(sym.canonical_order(), vec![h1, h2]);
        assert_eq!(sym.auto_group_order_exact(), Some(2));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(2)));
    }

    #[test]
    fn hf_distinct() {
        let mut b = MoleculeBuilder::new();
        let h = b.add_atom(AtomPattern::new(Element::H));
        let f = b.add_atom(AtomPattern::new(Element::F));
        b.add_bond(h, f, BondPattern::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 2);
        assert!(!sym.same_orbit(h, f));
        assert_eq!(sym.orbit_representative(h), h);
        assert_eq!(sym.orbit_representative(f), f);
        assert_eq!(sym.orbit_partition(), vec![vec![h], vec![f]]);
        assert_eq!(sym.canonical_order(), vec![h, f]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(1)));
    }

    #[test]
    fn square_uniform_bonds() {
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..4)
            .map(|_| b.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..4 {
            b.add_bond(c[i], c[(i + 1) % 4], BondPattern::new(1));
        }
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(c[0], c[1]));
        assert!(sym.same_orbit(c[1], c[2]));
        assert!(sym.same_orbit(c[2], c[3]));
        assert!(sym.same_orbit(c[3], c[0]));
        assert_eq!(sym.orbit_representative(c[0]), c[0]);
        assert_eq!(sym.orbit_representative(c[1]), c[0]);
        assert_eq!(sym.orbit_representative(c[2]), c[0]);
        assert_eq!(sym.orbit_representative(c[3]), c[0]);
        let partition = sym.orbit_partition();
        assert_eq!(partition, vec![vec![c[0], c[1], c[2], c[3]]]);
        assert_eq!(sym.canonical_order(), vec![c[0], c[1], c[3], c[2]]);
        assert_eq!(sym.auto_group_order_exact(), Some(8));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(8)));
    }

    #[test]
    fn linear_mixed_bond_orders() {
        // C=C-C: all three atoms in distinct orbits
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(AtomPattern::new(Element::C));
        let c2 = b.add_atom(AtomPattern::new(Element::C));
        let c3 = b.add_atom(AtomPattern::new(Element::C));
        b.add_bond(c1, c2, BondPattern::new(2));
        b.add_bond(c2, c3, BondPattern::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 3);
        assert!(!sym.same_orbit(c1, c2));
        assert!(!sym.same_orbit(c2, c3));
        assert!(!sym.same_orbit(c3, c1));
        assert_eq!(sym.orbit_representative(c1), c1);
        assert_eq!(sym.orbit_representative(c2), c2);
        assert_eq!(sym.orbit_representative(c3), c3);
        let partition = sym.orbit_partition();
        assert_eq!(partition, vec![vec![c1], vec![c2], vec![c3]]);
        assert_eq!(sym.canonical_order(), vec![c1, c3, c2]);
        assert_eq!(sym.auto_group_order_exact(), Some(1));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(1)));
    }

    #[test]
    fn alternating_cycle() {
        // C=C-C=C cycle: all atoms equivalent
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..4)
            .map(|_| b.add_atom(AtomPattern::new(Element::C)))
            .collect();
        b.add_bond(c[0], c[1], BondPattern::new(2));
        b.add_bond(c[1], c[2], BondPattern::new(1));
        b.add_bond(c[2], c[3], BondPattern::new(2));
        b.add_bond(c[3], c[0], BondPattern::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(c[0], c[1]));
        assert!(sym.same_orbit(c[1], c[2]));
        assert!(sym.same_orbit(c[2], c[3]));
        assert!(sym.same_orbit(c[3], c[0]));
        assert_eq!(sym.orbit_representative(c[0]), c[0]);
        assert_eq!(sym.orbit_representative(c[1]), c[0]);
        assert_eq!(sym.orbit_representative(c[2]), c[0]);
        assert_eq!(sym.orbit_representative(c[3]), c[0]);
        let partition = sym.orbit_partition();
        assert_eq!(partition, vec![vec![c[0], c[1], c[2], c[3]]]);
        assert_eq!(sym.canonical_order(), vec![c[0], c[3], c[1], c[2]]);
        assert_eq!(sym.auto_group_order_exact(), Some(4));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(4)));
    }

    #[test]
    fn water() {
        // H-O-H: two orbits (O and the two H's)
        let mut b = MoleculeBuilder::new();
        let h1 = b.add_atom(AtomPattern::new(Element::H));
        let o = b.add_atom(AtomPattern::new(Element::O));
        let h2 = b.add_atom(AtomPattern::new(Element::H));
        b.add_bond(h1, o, BondPattern::new(1));
        b.add_bond(o, h2, BondPattern::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 2);
        assert!(sym.same_orbit(h1, h2));
        assert!(!sym.same_orbit(h1, o));
        assert!(!sym.same_orbit(h2, o));
        assert_eq!(sym.orbit_representative(h1), h1);
        assert_eq!(sym.orbit_representative(o), o);
        assert_eq!(sym.orbit_representative(h2), h1);
        let partition = sym.orbit_partition();
        assert_eq!(partition, vec![vec![h1, h2], vec![o]]);
        assert_eq!(sym.canonical_order(), vec![h1, h2, o]);
        assert_eq!(sym.auto_group_order_exact(), Some(2));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(2)));
    }

    #[test]
    fn benzene_ring_uniform_bonds() {
        // 6 C atoms in a ring, all single bonds: all equivalent
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..6)
            .map(|_| b.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(c[i], c[(i + 1) % 6], BondPattern::new(1));
        }
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(c[0], c[1]));
        assert!(sym.same_orbit(c[1], c[2]));
        assert!(sym.same_orbit(c[2], c[3]));
        assert!(sym.same_orbit(c[3], c[4]));
        assert!(sym.same_orbit(c[4], c[5]));
        assert!(sym.same_orbit(c[5], c[0]));
        assert_eq!(sym.orbit_representative(c[0]), c[0]);
        assert_eq!(sym.orbit_representative(c[1]), c[0]);
        assert_eq!(sym.orbit_representative(c[2]), c[0]);
        assert_eq!(sym.orbit_representative(c[3]), c[0]);
        assert_eq!(sym.orbit_representative(c[4]), c[0]);
        assert_eq!(sym.orbit_representative(c[5]), c[0]);
        let partition = sym.orbit_partition();
        assert_eq!(partition, vec![vec![c[0], c[1], c[2], c[3], c[4], c[5]]]);
        assert_eq!(
            sym.canonical_order(),
            vec![c[0], c[2], c[4], c[3], c[1], c[5]]
        );
        assert_eq!(sym.auto_group_order_exact(), Some(12));
        assert!(matches!(sym.auto_group_order(), AutoGroupOrder::Exact(12)));
    }

    #[test]
    fn canonical_order_deterministic() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(AtomPattern::new(Element::C));
        b.add_atom(AtomPattern::new(Element::N));
        b.add_atom(AtomPattern::new(Element::O));
        let order1 = compute_symmetry(&b).canonical_order();
        let order2 = compute_symmetry(&b).canonical_order();
        assert_eq!(order1, order2);
        assert_eq!(order1.len(), 3);
    }
}
