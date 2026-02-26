//! Graph automorphism and canonical labeling via nauty.

use std::collections::{BTreeMap, HashSet};
use std::os::raw::c_int;

use nauty_Traces_sys::*;
use umol_data::{Element, SpinMultiplicity};

use super::atom::AtomBuilder;
use super::bond::Bond;
use super::molecule::{AtomIndex, MoleculeBuilder};
use crate::bond::BondDonation;

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

fn atom_color(ab: &AtomBuilder) -> VertexColor {
    VertexColor::Atom {
        element: ab.element(),
        isotope_mass: ab.isotope_mass(),
        charge: ab.charge(),
        hydrogens: ab.hydrogens(),
        unpaired_electrons: ab.unpaired_electrons(),
        multiplicity: ab.multiplicity(),
    }
}

fn bond_color(b: &Bond) -> VertexColor {
    VertexColor::Bond {
        order: b.order(),
        donation: match b.donation() {
            None => 0,
            Some(BondDonation::Shared) => 1,
            Some(BondDonation::Donating) => 2,
            Some(BondDonation::Accepting) => 3,
        },
    }
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
}

impl GraphSymmetry {
    /// Canonical orbit representative for an atom.
    pub fn orbit_representative(&self, atom: AtomIndex) -> AtomIndex {
        let ni = self.atom_to_nauty[atom.index()];
        debug_assert!(ni < self.n_atoms);
        self.nauty_to_atom[self.orbits[ni] as usize]
    }

    /// Whether two atoms belong to the same orbit.
    pub fn same_orbit(&self, a: AtomIndex, b: AtomIndex) -> bool {
        self.orbits[self.atom_to_nauty[a.index()]] == self.orbits[self.atom_to_nauty[b.index()]]
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

    /// Number of distinct orbits among atoms.
    pub fn num_orbits(&self) -> usize {
        self.num_orbits
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
        for i in 0..n_atoms {
            reps.insert(orbits[i]);
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
    }
}

#[cfg(test)]
mod tests {
    use umol_data::Element;

    use super::*;

    #[test]
    fn single_atom() {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(AtomBuilder::new(Element::C));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert_eq!(sym.canonical_order(), vec![a]);
    }

    #[test]
    fn h2_equivalent() {
        let mut b = MoleculeBuilder::new();
        let h1 = b.add_atom(AtomBuilder::new(Element::H));
        let h2 = b.add_atom(AtomBuilder::new(Element::H));
        b.add_bond(h1, h2, Bond::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert!(sym.same_orbit(h1, h2));
    }

    #[test]
    fn hf_distinct() {
        let mut b = MoleculeBuilder::new();
        let h = b.add_atom(AtomBuilder::new(Element::H));
        let f = b.add_atom(AtomBuilder::new(Element::F));
        b.add_bond(h, f, Bond::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 2);
        assert!(!sym.same_orbit(h, f));
    }

    #[test]
    fn square_uniform_bonds() {
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..4)
            .map(|_| b.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..4 {
            b.add_bond(c[i], c[(i + 1) % 4], Bond::new(1));
        }
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
    }

    #[test]
    fn linear_mixed_bond_orders() {
        // C=C-C: all three atoms in distinct orbits
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(AtomBuilder::new(Element::C));
        let c2 = b.add_atom(AtomBuilder::new(Element::C));
        let c3 = b.add_atom(AtomBuilder::new(Element::C));
        b.add_bond(c1, c2, Bond::new(2));
        b.add_bond(c2, c3, Bond::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 3);
    }

    #[test]
    fn alternating_cycle() {
        // C=C-C=C cycle: all atoms equivalent
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..4)
            .map(|_| b.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        b.add_bond(c[0], c[1], Bond::new(2));
        b.add_bond(c[1], c[2], Bond::new(1));
        b.add_bond(c[2], c[3], Bond::new(2));
        b.add_bond(c[3], c[0], Bond::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
    }

    #[test]
    fn water_like() {
        // H-O-H: two orbits (O and the two H's)
        let mut b = MoleculeBuilder::new();
        let h1 = b.add_atom(AtomBuilder::new(Element::H));
        let o = b.add_atom(AtomBuilder::new(Element::O));
        let h2 = b.add_atom(AtomBuilder::new(Element::H));
        b.add_bond(h1, o, Bond::new(1));
        b.add_bond(o, h2, Bond::new(1));
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 2);
        assert!(sym.same_orbit(h1, h2));
        assert!(!sym.same_orbit(h1, o));
        let partition = sym.orbit_partition();
        assert_eq!(partition.len(), 2);
    }

    #[test]
    fn canonical_order_deterministic() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(AtomBuilder::new(Element::C));
        b.add_atom(AtomBuilder::new(Element::N));
        b.add_atom(AtomBuilder::new(Element::O));
        let order1 = compute_symmetry(&b).canonical_order();
        let order2 = compute_symmetry(&b).canonical_order();
        assert_eq!(order1, order2);
        assert_eq!(order1.len(), 3);
    }

    #[test]
    fn empty_graph() {
        let b = MoleculeBuilder::new();
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 0);
        assert!(sym.canonical_order().is_empty());
    }

    #[test]
    fn benzene_ring_uniform_bonds() {
        // 6 C atoms in a ring, all single bonds: all equivalent
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..6)
            .map(|_| b.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(c[i], c[(i + 1) % 6], Bond::new(1));
        }
        let sym = compute_symmetry(&b);
        assert_eq!(sym.num_orbits(), 1);
        assert_eq!(sym.canonical_order().len(), 6);
    }
}
