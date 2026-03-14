//! Ring detection primitives for GraphIR.
//!
//! Used for ring size queries and bounded ring enumeration.

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use super::molecule::builder::MoleculeBuilder;
use super::molecule::{enumerate_rings, AtomIndex, BondIndex, Molecule};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RingIndex(pub u32);

impl RingIndex {
    pub fn index(self) -> usize {
        self.0 as usize
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

#[derive(Debug, Clone)]
/// Set of simple cycles up to max_ring_size. Not a minimal cycle basis.
pub struct MoleculeRings {
    pub max_ring_size: usize,
    pub ring_atoms: Vec<Vec<AtomIndex>>,
    pub ring_bonds: Vec<Vec<BondIndex>>,
    pub atom_rings: BTreeMap<AtomIndex, Vec<RingIndex>>,
    pub bond_rings: BTreeMap<BondIndex, Vec<RingIndex>>,
}

impl MoleculeRings {
    pub fn empty() -> Self {
        Self {
            max_ring_size: 0,
            ring_atoms: Vec::new(),
            ring_bonds: Vec::new(),
            atom_rings: BTreeMap::new(),
            bond_rings: BTreeMap::new(),
        }
    }

    pub fn from_builder(builder: &MoleculeBuilder, max_ring_size: usize) -> Self {
        let bcc = builder.biconnected_components();
        let full_adj = builder.adjacency_list();
        Self::build(&bcc, &full_adj, max_ring_size, |a, b| {
            builder.connecting_bond_index(a, b)
        })
    }

    pub fn from_molecule(molecule: &Molecule, max_ring_size: usize) -> Self {
        let bcc = molecule.biconnected_components();
        let full_adj = molecule.adjacency_list();
        Self::build(&bcc, &full_adj, max_ring_size, |a, b| {
            molecule.connecting_bond_index(a, b)
        })
    }

    fn build(
        bcc: &[Vec<AtomIndex>],
        full_adj: &HashMap<AtomIndex, Vec<AtomIndex>>,
        max_ring_size: usize,
        bond_lookup: impl Fn(AtomIndex, AtomIndex) -> Option<BondIndex>,
    ) -> Self {
        let mut all_rings: Vec<Vec<AtomIndex>> = Vec::new();

        for component in bcc {
            let component_set: HashSet<AtomIndex> = component.iter().copied().collect();
            let mut sub_adj: HashMap<AtomIndex, Vec<AtomIndex>> = HashMap::new();
            for &atom in component {
                let neighbors: Vec<AtomIndex> = full_adj
                    .get(&atom)
                    .map(|ns| {
                        ns.iter()
                            .copied()
                            .filter(|n| component_set.contains(n))
                            .collect()
                    })
                    .unwrap_or_default();
                sub_adj.insert(atom, neighbors);
            }
            all_rings.extend(enumerate_rings(&sub_adj, max_ring_size));
        }

        let mut ring_atoms: BTreeMap<AtomIndex, Vec<RingIndex>> = BTreeMap::new();
        let mut ring_bonds: BTreeMap<BondIndex, Vec<RingIndex>> = BTreeMap::new();
        let mut cycle_bonds: Vec<Vec<BondIndex>> = Vec::with_capacity(all_rings.len());

        for (idx, ring) in all_rings.iter().enumerate() {
            let ring_idx = RingIndex(idx as u32);
            for &atom in ring {
                ring_atoms.entry(atom).or_default().push(ring_idx);
            }
            let n = ring.len();
            let mut bonds = Vec::with_capacity(n);
            for i in 0..n {
                let a = ring[i];
                let b = ring[(i + 1) % n];
                if let Some(bond) = bond_lookup(a, b) {
                    ring_bonds.entry(bond).or_default().push(ring_idx);
                    bonds.push(bond);
                }
            }
            cycle_bonds.push(bonds);
        }

        Self {
            max_ring_size,
            ring_atoms: all_rings,
            ring_bonds: cycle_bonds,
            atom_rings: ring_atoms,
            bond_rings: ring_bonds,
        }
    }

    pub fn ring_count(&self) -> usize {
        self.ring_atoms.len()
    }

    pub fn ring_indices(&self) -> impl Iterator<Item = RingIndex> {
        (0..self.ring_atoms.len()).map(|i| RingIndex(i as u32))
    }

    pub fn ring(&self, idx: RingIndex) -> Option<&[AtomIndex]> {
        self.ring_atoms.get(idx.index()).map(|v| v.as_slice())
    }

    pub fn shared_atoms(&self, a: RingIndex, b: RingIndex) -> Vec<AtomIndex> {
        let (Some(ra), Some(rb)) = (self.ring(a), self.ring(b)) else {
            return Vec::new();
        };
        let atoms_a: HashSet<AtomIndex> = ra.iter().copied().collect();
        rb.iter()
            .copied()
            .filter(|atom| atoms_a.contains(atom))
            .collect()
    }

    pub fn shared_bonds(&self, a: RingIndex, b: RingIndex) -> Vec<BondIndex> {
        let (Some(ba), Some(bb)) = (
            self.ring_bonds.get(a.index()),
            self.ring_bonds.get(b.index()),
        ) else {
            return Vec::new();
        };
        let bonds_b: HashSet<BondIndex> = bb.iter().copied().collect();
        ba.iter()
            .copied()
            .filter(|bond| bonds_b.contains(bond))
            .collect()
    }

    pub fn ring_relation(&self, a: RingIndex, b: RingIndex) -> RingRelation {
        if a == b {
            return RingRelation::Identical;
        }
        let shared = self.shared_bonds(a, b);
        if shared.is_empty() {
            return match self.shared_atoms(a, b).len() {
                0 => RingRelation::Disjoint,
                1 => RingRelation::Spiro,
                _ => RingRelation::MultiSpiro,
            };
        }
        if shared.len() == 1 {
            return RingRelation::Fused;
        }
        let Some(bonds_a) = self.ring_bonds.get(a.index()) else {
            return RingRelation::Disjoint;
        };
        let n = bonds_a.len();
        let shared_set: HashSet<BondIndex> = shared.into_iter().collect();
        let mut runs = 0usize;
        for i in 0..n {
            if shared_set.contains(&bonds_a[i]) && !shared_set.contains(&bonds_a[(i + n - 1) % n]) {
                runs += 1;
            }
        }
        if runs <= 1 {
            RingRelation::Bridged
        } else {
            RingRelation::Noncontiguous
        }
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
            .ring_indices()
            .filter(|&j| j != i && self.are_spiro(i, j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn ring_fused_neighbors(&self, i: RingIndex) -> Vec<RingIndex> {
        let mut result: Vec<RingIndex> = self
            .ring_indices()
            .filter(|&j| j != i && self.are_fused(i, j))
            .collect();
        result.sort_unstable();
        result
    }

    pub fn ring_bridged_neighbors(&self, i: RingIndex) -> Vec<RingIndex> {
        let mut result: Vec<RingIndex> = self
            .ring_indices()
            .filter(|&j| j != i && self.are_bridged(i, j))
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

    pub fn is_ring_atom(&self, atom: AtomIndex) -> bool {
        self.atom_rings.contains_key(&atom)
    }

    pub fn atom_smallest_ring_size(&self, atom: AtomIndex) -> Option<usize> {
        self.atom_rings.get(&atom).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.ring_atoms[i.index()].len())
                .min()
        })
    }

    pub fn is_ring_bond(&self, bond: BondIndex) -> bool {
        self.bond_rings.contains_key(&bond)
    }

    pub fn bond_smallest_ring_size(&self, bond: BondIndex) -> Option<usize> {
        self.bond_rings.get(&bond).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.ring_atoms[i.index()].len())
                .min()
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::atom;
    use crate::graph_ir::atom::AtomBuilder;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::config::ResolveConfig;

    #[fixture]
    fn empty_builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    #[fixture]
    fn single_atom_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomBuilder::new(Element::C));
        builder
    }

    #[fixture]
    fn chain_builder(#[default(5)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..n - 1 {
            builder.add_bond_unchecked(atoms[i], atoms[i + 1], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn ring_builder(#[default(6)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    // bicyclohexyl
    fn disjoint_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..12)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        for i in 6..12 {
            builder.add_bond_unchecked(
                atoms[i],
                atoms[6 + ((i + 1 - 6) % 6)],
                BondBuilder::new(1, None),
            );
        }
        builder
    }

    #[fixture]
    // spiro[3.4]pentane
    fn spiro_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..5)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    // naphthalene
    fn fused_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..10)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let ring1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (a, b) in ring1_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let ring2_edges = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
        for (a, b) in ring2_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    // bicyclo[1.1.1]pentane
    fn bridged_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..5)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn multi_spiro_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (0, 3),
            (3, 2),
            (0, 4),
            (4, 2),
            (0, 5),
            (5, 2),
        ];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    fn cubane_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..8)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
    }

    #[fixture]
    // methylcyclohexane
    fn substituted_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomBuilder::new(Element::C)))
            .collect();
        let methyl = builder.add_atom(AtomBuilder::new(Element::C));
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        // Attach methyl group to C1 (atom 0)
        builder.add_bond_unchecked(atoms[0], methyl, BondBuilder::new(1, None));
        builder
    }

    #[fixture]
    fn ring_molecule(#[default(6)] n: usize) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n).map(|_| builder.add_atom(atom!("{CH2v2}"))).collect();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondBuilder::new(1, None));
        }
        builder
            .build(&ResolveConfig::default())
            .expect("test molecule should build")
    }

    #[fixture]
    fn ring_molecule_rings(#[default(6)] n: usize) -> MoleculeRings {
        let molecule = ring_builder(n);
        MoleculeRings::from_builder(&molecule, n)
    }

    #[fixture]
    fn disjoint_molecule_rings() -> MoleculeRings {
        let molecule = disjoint_builder();
        MoleculeRings::from_builder(&molecule, 10)
    }

    #[fixture]
    fn spiro_molecule_rings() -> MoleculeRings {
        let molecule = spiro_builder();
        MoleculeRings::from_builder(&molecule, 10)
    }

    #[fixture]
    fn fused_molecule_rings() -> MoleculeRings {
        let molecule = fused_builder();
        MoleculeRings::from_builder(&molecule, 10)
    }

    #[fixture]
    fn bridged_molecule_rings() -> MoleculeRings {
        let molecule = bridged_builder();
        MoleculeRings::from_builder(&molecule, 5)
    }

    #[fixture]
    fn cubane_molecule_rings() -> MoleculeRings {
        let molecule = cubane_builder();
        MoleculeRings::from_builder(&molecule, 8)
    }

    #[fixture]
    fn multi_spiro_molecule_rings() -> MoleculeRings {
        let molecule = multi_spiro_builder();
        MoleculeRings::from_builder(&molecule, 4)
    }

    #[rstest]
    #[case::empty_max_size_0(empty_builder(), 0, 0)]
    #[case::empty_max_size_3(empty_builder(), 3, 0)]
    #[case::single_atom_max_size_0(single_atom_builder(), 0, 0)]
    #[case::single_atom_max_size_3(single_atom_builder(), 3, 0)]
    #[case::pentane_max_size_0(chain_builder(5), 0, 0)]
    #[case::pentane_max_size_3(chain_builder(5), 3, 0)]
    #[case::cyclohexane_max_size_0(ring_builder(6), 0, 0)]
    #[case::cyclohexane_max_size_3(ring_builder(6), 3, 0)]
    #[case::cyclohexane_max_size_6(ring_builder(6), 6, 1)]
    #[case::cyclohexane_max_size_8(ring_builder(6), 8, 1)]
    #[case::disjoint_max_size_0(disjoint_builder(), 0, 0)]
    #[case::disjoint_max_size_3(disjoint_builder(), 3, 0)]
    #[case::disjoint_max_size_6(disjoint_builder(), 6, 2)]
    #[case::disjoint_max_size_8(disjoint_builder(), 8, 2)]
    #[case::spiro_max_size_3(spiro_builder(), 3, 2)]
    #[case::spiro_max_size_6(spiro_builder(), 6, 2)]
    #[case::fused_max_size_3(fused_builder(), 3, 0)]
    #[case::fused_max_size_6(fused_builder(), 6, 2)]
    #[case::fused_max_size_10(fused_builder(), 10, 3)]
    #[case::fused_max_size_20(fused_builder(), 20, 3)]
    #[case::bridged_max_size_3(bridged_builder(), 3, 0)]
    #[case::bridged_max_size_4(bridged_builder(), 4, 3)]
    #[case::multi_spiro_max_size_3(multi_spiro_builder(), 3, 0)]
    #[case::multi_spiro_max_size_4(multi_spiro_builder(), 6, 6)]
    #[case::multi_spiro_max_size_20(multi_spiro_builder(), 20, 6)]
    #[case::cubane_max_size_4(cubane_builder(), 4, 6)]
    #[case::cubane_max_size_6(cubane_builder(), 6, 22)]
    #[case::cubane_max_size_8(cubane_builder(), 8, 28)]
    #[case::cubane_max_size_20(cubane_builder(), 20, 28)]
    fn test_molecule_rings_from_builder(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rings = MoleculeRings::from_builder(&builder, max_ring_size);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::cyclohexane_max_size_0(ring_molecule(6), 0, 0)]
    #[case::cyclohexane_max_size_3(ring_molecule(6), 3, 0)]
    #[case::cyclohexane_max_size_6(ring_molecule(6), 6, 1)]
    #[case::cyclohexane_max_size_8(ring_molecule(6), 8, 1)]
    fn test_molecule_rings_from_molecule(
        #[case] molecule: Molecule,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rings = MoleculeRings::from_molecule(&molecule, max_ring_size);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::ring(ring_molecule_rings(6), vec![RingIndex(0)])]
    #[case::fused(fused_molecule_rings(), vec![RingIndex(0), RingIndex(1), RingIndex(2)])]
    fn test_molecule_rings_ring_indices(
        #[case] rings: MoleculeRings,
        #[case] expected: Vec<RingIndex>,
    ) {
        let ring_indices: Vec<RingIndex> = rings.ring_indices().collect();
        assert_eq!(ring_indices, expected);
    }

    #[rstest]
    #[case::ring(ring_molecule_rings(6), RingIndex(0),
         Some(vec![AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2), AtomIndex::new(3),
              AtomIndex::new(4), AtomIndex::new(5)]))]
    #[case::non_existent(ring_molecule_rings(6), RingIndex(1), None)]
    #[case::fused(fused_molecule_rings(), RingIndex(0),
        Some(vec![AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2),
             AtomIndex::new(3), AtomIndex::new(4), AtomIndex::new(5)]))]
    fn test_molecule_rings_ring(
        #[case] rings: MoleculeRings,
        #[case] ring_index: RingIndex,
        #[case] expected: Option<Vec<AtomIndex>>,
    ) {
        let ring = rings.ring(ring_index).map(|v| v.to_vec());
        assert_eq!(ring, expected);
    }

    #[rstest]
    #[case::identical(ring_molecule_rings(6), RingIndex(0), RingIndex(0),
        vec![AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2), AtomIndex::new(3), AtomIndex::new(4), AtomIndex::new(5)])]
    #[case::disjoint(disjoint_molecule_rings(), RingIndex(0), RingIndex(1), vec![])]
    #[case::spiro(spiro_molecule_rings(), RingIndex(0), RingIndex(1), vec![AtomIndex::new(0)])]
    #[case::fused(fused_molecule_rings(), RingIndex(0), RingIndex(1), vec![AtomIndex::new(3), AtomIndex::new(4)])]
    #[case::bridged(bridged_molecule_rings(), RingIndex(0), RingIndex(1), vec![AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2)])]
    #[case::multi_spiro(multi_spiro_molecule_rings(), RingIndex(0), RingIndex(5), vec![AtomIndex::new(0), AtomIndex::new(2)])]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0), RingIndex(25), vec![AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2), AtomIndex::new(3)])]
    #[case::non_existent(ring_molecule_rings(6), RingIndex(0), RingIndex(2), vec![])]
    fn test_molecule_rings_shared_atoms(
        #[case] rings: MoleculeRings,
        #[case] ring_index_a: RingIndex,
        #[case] ring_index_b: RingIndex,
        #[case] expected: Vec<AtomIndex>,
    ) {
        let mut shared_atoms = rings.shared_atoms(ring_index_a, ring_index_b);
        shared_atoms.sort_unstable();
        assert_eq!(shared_atoms, expected);
    }

    #[rstest]
    #[case::identical(ring_molecule_rings(6), RingIndex(0), RingIndex(0),
        vec![BondIndex::new(0), BondIndex::new(1), BondIndex::new(2), BondIndex::new(3), BondIndex::new(4), BondIndex::new(5)])]
    #[case::disjoint(disjoint_molecule_rings(), RingIndex(0), RingIndex(1), vec![])]
    #[case::spiro(spiro_molecule_rings(), RingIndex(0), RingIndex(1), vec![])]
    #[case::fused(fused_molecule_rings(), RingIndex(0), RingIndex(1), vec![BondIndex::new(3)])]
    #[case::bridged(bridged_molecule_rings(), RingIndex(0), RingIndex(1), vec![BondIndex::new(0), BondIndex::new(1)])]
    #[case::multi_spiro(multi_spiro_molecule_rings(), RingIndex(0), RingIndex(5), vec![])]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0), RingIndex(25), vec![BondIndex::new(0), BondIndex::new(2)])]
    #[case::non_existent(ring_molecule_rings(6), RingIndex(0), RingIndex(2), vec![])]
    fn test_molecule_rings_shared_bonds(
        #[case] rings: MoleculeRings,
        #[case] ring_index_a: RingIndex,
        #[case] ring_index_b: RingIndex,
        #[case] expected: Vec<BondIndex>,
    ) {
        let mut shared_bonds = rings.shared_bonds(ring_index_a, ring_index_b);
        shared_bonds.sort_unstable();
        assert_eq!(shared_bonds, expected);
    }

    #[rstest]
    #[case::identical(
        ring_molecule_rings(6),
        RingIndex(0),
        RingIndex(0),
        RingRelation::Identical
    )]
    #[case::disjoint(
        disjoint_molecule_rings(),
        RingIndex(0),
        RingIndex(1),
        RingRelation::Disjoint
    )]
    #[case::spiro(
        spiro_molecule_rings(),
        RingIndex(0),
        RingIndex(1),
        RingRelation::Spiro
    )]
    #[case::fused(
        fused_molecule_rings(),
        RingIndex(0),
        RingIndex(1),
        RingRelation::Fused
    )]
    #[case::bridged(
        bridged_molecule_rings(),
        RingIndex(0),
        RingIndex(1),
        RingRelation::Bridged
    )]
    #[case::multi_spiro(
        multi_spiro_molecule_rings(),
        RingIndex(0),
        RingIndex(5),
        RingRelation::MultiSpiro
    )]
    #[case::cubane(
        cubane_molecule_rings(),
        RingIndex(0),
        RingIndex(25),
        RingRelation::Noncontiguous
    )]
    #[case::non_existent(
        ring_molecule_rings(6),
        RingIndex(0),
        RingIndex(2),
        RingRelation::Disjoint
    )]
    fn test_molecule_rings_ring_relation(
        #[case] rings: MoleculeRings,
        #[case] ring_index_a: RingIndex,
        #[case] ring_index_b: RingIndex,
        #[case] expected: RingRelation,
    ) {
        let relation = rings.ring_relation(ring_index_a, ring_index_b);
        assert_eq!(relation, expected);
    }

    #[rstest]
    #[case::single_ring(ring_molecule_rings(6), RingIndex(0), vec![])]
    #[case::spiro_1(spiro_molecule_rings(), RingIndex(0), vec![RingIndex(1)])]
    #[case::spiro_2(spiro_molecule_rings(), RingIndex(1), vec![RingIndex(0)])]
    #[case::fused_1(fused_molecule_rings(), RingIndex(0), vec![])]
    #[case::fused_2(fused_molecule_rings(), RingIndex(1), vec![])]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0), vec![])]
    fn test_molecule_rings_spiro_neighbors(
        #[case] rings: MoleculeRings,
        #[case] ring_index: RingIndex,
        #[case] expected: Vec<RingIndex>,
    ) {
        let mut neighbors = rings.ring_spiro_neighbors(ring_index);
        neighbors.sort_unstable();
        assert_eq!(neighbors, expected);
    }

    #[rstest]
    #[case::single_ring(ring_molecule_rings(6), RingIndex(0), vec![])]
    #[case::spiro_1(spiro_molecule_rings(), RingIndex(0), vec![])]
    #[case::spiro_2(spiro_molecule_rings(), RingIndex(1), vec![])]
    #[case::fused_1(fused_molecule_rings(), RingIndex(0), vec![RingIndex(1)])]
    #[case::fused_2(fused_molecule_rings(), RingIndex(1), vec![RingIndex(0)])]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0),
        vec![RingIndex(1), RingIndex(2), RingIndex(3), RingIndex(4), RingIndex(13), RingIndex(17), RingIndex(20), RingIndex(21)])]
    fn test_molecule_rings_fused_neighbors(
        #[case] rings: MoleculeRings,
        #[case] ring_index: RingIndex,
        #[case] expected: Vec<RingIndex>,
    ) {
        let mut neighbors = rings.ring_fused_neighbors(ring_index);
        neighbors.sort_unstable();
        assert_eq!(neighbors, expected);
    }

    #[rstest]
    #[case::single_ring(ring_molecule_rings(6), RingIndex(0), vec![])]
    #[case::fused_1(fused_molecule_rings(), RingIndex(0), vec![RingIndex(2)])]
    #[case::fused_2(fused_molecule_rings(), RingIndex(1), vec![RingIndex(2)])]
    #[case::bridged_1(bridged_molecule_rings(), RingIndex(0), vec![RingIndex(1), RingIndex(2)])]
    #[case::bridged_2(bridged_molecule_rings(), RingIndex(1), vec![RingIndex(0), RingIndex(2)])]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0), vec![
        RingIndex(6), RingIndex(7), RingIndex(8), RingIndex(9), RingIndex(10), RingIndex(11), RingIndex(12), RingIndex(14), RingIndex(15),
        RingIndex(16), RingIndex(18), RingIndex(19), RingIndex(22), RingIndex(23), RingIndex(24), RingIndex(26)])]
    fn test_molecule_rings_bridged_neighbors(
        #[case] rings: MoleculeRings,
        #[case] ring_index: RingIndex,
        #[case] expected: Vec<RingIndex>,
    ) {
        let mut neighbors = rings.ring_bridged_neighbors(ring_index);
        neighbors.sort_unstable();
        assert_eq!(neighbors, expected);
    }

    #[rstest]
    #[case::cyclohexane(ring_molecule_rings(6), vec![1])]
    #[case::naphthalene(fused_molecule_rings(), vec![1, 2])]
    #[case::cubane(cubane_molecule_rings(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 18])]
    fn test_molecule_rings_fused_components(
        #[case] rings: MoleculeRings,
        #[case] mut expected_sizes: Vec<usize>,
    ) {
        let components = rings.fused_components();
        let mut sizes: Vec<usize> = components.iter().map(|c| c.len()).collect();
        sizes.sort_unstable();
        expected_sizes.sort_unstable();
        assert_eq!(sizes, expected_sizes);
    }

    #[rstest]
    #[case::cyclohexane(ring_molecule_rings(6), RingIndex(0), 1)]
    #[case::naphthalene_0(fused_molecule_rings(), RingIndex(0), 2)]
    #[case::naphthalene_1(fused_molecule_rings(), RingIndex(1), 2)]
    #[case::cubane_0(cubane_molecule_rings(), RingIndex(0), 18)]
    fn test_molecule_rings_ring_fused_component(
        #[case] rings: MoleculeRings,
        #[case] ring_index: RingIndex,
        #[case] expected_size: usize,
    ) {
        let component = rings.ring_fused_component(ring_index);
        assert!(component.contains(&ring_index));
        assert_eq!(component.len(), expected_size);
    }

    #[rstest]
    #[case::empty(empty_builder(), 0, false)]
    #[case::empty_atom_1(empty_builder(), 1, false)]
    #[case::single_atom(single_atom_builder(), 0, false)]
    #[case::cyclohexane_ring(ring_builder(6), 0, true)]
    #[case::cyclohexane_ring_atom_3(ring_builder(6), 3, true)]
    #[case::cyclohexane_not_in(ring_builder(6), 6, false)]
    #[case::naphthalene_ring_0_atom_0(fused_builder(), 0, true)]
    #[case::naphthalene_ring_0_atom_3(fused_builder(), 3, true)]
    #[case::naphthalene_ring_0_atom_6(fused_builder(), 6, true)]
    #[case::naphthalene_not_in(fused_builder(), 10, false)]
    #[case::methylcyclohexane_ring(substituted_builder(), 0, true)]
    #[case::methylcyclohexane_ring_atom_5(substituted_builder(), 5, true)]
    #[case::methylcyclohexane_methyl(substituted_builder(), 6, false)]
    #[case::nonexistent(ring_builder(6), 6, false)]
    fn test_molecule_rings_is_ring_atom(
        #[case] builder: MoleculeBuilder,
        #[case] atom_index: usize,
        #[case] expected: bool,
    ) {
        let rings = MoleculeRings::from_builder(&builder, 10);
        let atom = AtomIndex::new(atom_index);
        assert_eq!(rings.is_ring_atom(atom), expected);
    }

    #[rstest]
    #[case::first_empty(empty_builder(), 10, 0, None)]
    #[case::first_single_atom(single_atom_builder(), 10, 0, None)]
    #[case::first_pentane(chain_builder(5), 10, 0, None)]
    #[case::first_cyclohexane(ring_builder(6), 10, 0, Some(6))]
    #[case::shared_hexagon(ring_builder(6), 10, 3, Some(6))]
    #[case::first_naphthalene(fused_builder(), 10, 0, Some(6))]
    #[case::shared_naphthalene(fused_builder(), 10, 3, Some(6))]
    #[case::first_cubane(cubane_builder(), 10, 0, Some(4))]
    #[case::first_spiro(spiro_builder(), 10, 0, Some(3))]
    #[case::first_bridged(bridged_builder(), 10, 0, Some(4))]
    fn test_molecule_rings_atom_smallest_ring_size(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] atom_index: usize,
        #[case] expected_ring_size: Option<usize>,
    ) {
        let rings = MoleculeRings::from_builder(&builder, max_ring_size);
        let atom = AtomIndex::new(atom_index);
        assert_eq!(rings.atom_smallest_ring_size(atom), expected_ring_size);
    }

    #[rstest]
    #[case::empty(empty_builder(), 0, false)]
    #[case::empty_bond_1(empty_builder(), 1, false)]
    #[case::single_atom(single_atom_builder(), 0, false)]
    #[case::cyclohexane_bond_0(ring_builder(6), 0, true)]
    #[case::cyclohexane_bond_3(ring_builder(6), 3, true)]
    #[case::naphthalene_ring_0_bond_0(fused_builder(), 0, true)]
    #[case::naphthalene_ring_0_bond_5(fused_builder(), 5, true)]
    #[case::naphthalene_ring_0_bond_6(fused_builder(), 6, true)]
    #[case::methylcyclohexane_ring(substituted_builder(), 0, true)]
    #[case::methylcyclohexane_ring_bond_5(substituted_builder(), 5, true)]
    #[case::methylcyclohexane_methyl(substituted_builder(), 6, false)]
    #[case::nonexistent(ring_builder(6), 6, false)]
    fn test_molecule_rings_is_ring_bond(
        #[case] builder: MoleculeBuilder,
        #[case] bond_index: usize,
        #[case] expected: bool,
    ) {
        let rings = MoleculeRings::from_builder(&builder, 10);
        let bond = BondIndex::new(bond_index);
        assert_eq!(rings.is_ring_bond(bond), expected);
    }

    #[rstest]
    #[case::empty(empty_builder(), 10, 0, None)]
    #[case::single_atom(single_atom_builder(), 10, 0, None)]
    #[case::pentane(chain_builder(5), 10, 0, None)]
    #[case::cyclohexane(ring_builder(6), 10, 0, Some(6))]
    #[case::cyclohexane_shared(ring_builder(6), 10, 3, Some(6))]
    #[case::naphthalene(fused_builder(), 10, 0, Some(6))]
    #[case::naphthalene_shared(fused_builder(), 10, 3, Some(6))]
    #[case::cubane(cubane_builder(), 10, 0, Some(4))]
    #[case::spiro(spiro_builder(), 10, 0, Some(3))]
    #[case::bridged(bridged_builder(), 10, 0, Some(4))]
    fn test_molecule_rings_bond_smallest_ring_size(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] bond_index: usize,
        #[case] expected_ring_size: Option<usize>,
    ) {
        let rings = MoleculeRings::from_builder(&builder, max_ring_size);
        let bond = BondIndex::new(bond_index);
        assert_eq!(rings.bond_smallest_ring_size(bond), expected_ring_size);
    }
}
