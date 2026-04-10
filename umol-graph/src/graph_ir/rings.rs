//! Ring detection primitives for GraphIR.
//!
//! Used for ring size queries and bounded ring enumeration.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use umol_data::Element;

use super::config::RingEnumerationStrategy;
use super::molecule::{AtomIndex, BondIndex, Molecule};
use super::molecule_builder::MoleculeBuilder;
use crate::algorithms::{biconnected_components, enumerate_simple_cycles};

#[derive(Debug, Clone)]
struct AtomAdjacency {
    neighbors: BTreeMap<AtomIndex, Vec<AtomIndex>>,
}

#[derive(Debug, Clone)]
struct DenseProjection {
    atoms: Vec<AtomIndex>,
    adj: Vec<Vec<usize>>,
}

impl AtomAdjacency {
    fn from_builder(builder: &MoleculeBuilder) -> Self {
        Self::from_map(builder.adjacency_list())
    }

    fn from_molecule(molecule: &Molecule) -> Self {
        Self::from_map(molecule.adjacency_list())
    }

    fn from_map(neighbors: HashMap<AtomIndex, Vec<AtomIndex>>) -> Self {
        Self {
            neighbors: neighbors.into_iter().collect(),
        }
    }

    fn atoms(&self) -> Vec<AtomIndex> {
        self.neighbors.keys().copied().collect()
    }

    fn induced(&self, atoms: &HashSet<AtomIndex>) -> Self {
        let mut induced_neighbors: BTreeMap<AtomIndex, Vec<AtomIndex>> = BTreeMap::new();
        for &atom in atoms {
            let mut neighbors: Vec<AtomIndex> = self
                .neighbors
                .get(&atom)
                .map(|ns| ns.iter().copied().filter(|n| atoms.contains(n)).collect())
                .unwrap_or_default();
            neighbors.sort_unstable();
            neighbors.dedup();
            induced_neighbors.insert(atom, neighbors);
        }
        Self {
            neighbors: induced_neighbors,
        }
    }

    fn to_dense(&self) -> DenseProjection {
        self.to_dense_for_atoms(&self.atoms())
    }

    fn to_dense_for_atoms(&self, atoms: &[AtomIndex]) -> DenseProjection {
        let atoms_sorted: Vec<AtomIndex> = atoms
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let atom_to_id: HashMap<AtomIndex, usize> = atoms_sorted
            .iter()
            .copied()
            .enumerate()
            .map(|(i, a)| (a, i))
            .collect();

        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); atoms_sorted.len()];
        for &atom in &atoms_sorted {
            let mut neighbors = self
                .neighbors
                .get(&atom)
                .map(|ns| {
                    ns.iter()
                        .filter_map(|&n| atom_to_id.get(&n).copied())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            neighbors.sort_unstable();
            neighbors.dedup();
            adj[atom_to_id[&atom]] = neighbors;
        }

        DenseProjection {
            atoms: atoms_sorted,
            adj,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RingIndex(pub u32);

impl RingIndex {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    atoms: Vec<AtomIndex>,
    bonds: Vec<BondIndex>,
}

impl Ring {
    pub fn new(atoms: Vec<AtomIndex>, bonds: Vec<BondIndex>) -> Result<Self, String> {
        if atoms.len() < 3 {
            return Err("ring must contain at least 3 atoms".to_string());
        }
        if atoms.len() != bonds.len() {
            return Err("ring atoms/bonds length mismatch".to_string());
        }
        Ok(Self { atoms, bonds })
    }

    pub fn atoms(&self) -> &[AtomIndex] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[BondIndex] {
        &self.bonds
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn shared_atoms(&self, other: &Ring) -> Vec<AtomIndex> {
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

    pub fn shared_bonds(&self, other: &Ring) -> Vec<BondIndex> {
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
/// Molecule-indexed ring set.
pub struct RingSet {
    pub family: RingFamily,
    pub scope: RingScope,
    pub max_ring_size: usize,
    pub rings: Vec<Ring>,
    pub atom_to_rings: BTreeMap<AtomIndex, Vec<RingIndex>>,
    pub bond_to_rings: BTreeMap<BondIndex, Vec<RingIndex>>,
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

        let mut atom_to_rings: BTreeMap<AtomIndex, Vec<RingIndex>> = BTreeMap::new();
        let mut bond_to_rings: BTreeMap<BondIndex, Vec<RingIndex>> = BTreeMap::new();
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

    pub fn induced_from_molecule_atoms(molecule: &Molecule, atoms: &[AtomIndex]) -> Self {
        if atoms.len() < 3 {
            return Self::empty_with(RingFamily::Induced, RingScope::AtomSubset);
        }

        let atom_set: HashSet<AtomIndex> = atoms.iter().copied().collect();
        if atom_set.len() < 3 {
            return Self::empty_with(RingFamily::Induced, RingScope::AtomSubset);
        }

        let full_adj = AtomAdjacency::from_molecule(molecule);
        let sub_adj = full_adj.induced(&atom_set);
        let dense = sub_adj.to_dense();
        let component_atoms = dense.atoms.clone();

        let mut bond_map: HashMap<(AtomIndex, AtomIndex), BondIndex> = HashMap::new();
        for bond in molecule.bond_indices() {
            if let Some((a, b)) = molecule.bond_atom_indices(bond) {
                bond_map.insert((a, b), bond);
                bond_map.insert((b, a), bond);
            }
        }

        let mut rings: Vec<Ring> =
            enumerate_simple_cycles(component_atoms.len(), &dense.adj, component_atoms.len())
                .into_iter()
                .filter(|cycle| is_induced_cycle(cycle, &dense.adj))
                .filter_map(|cycle| {
                    let ring_atoms: Vec<AtomIndex> =
                        cycle.into_iter().map(|i| component_atoms[i]).collect();
                    let n = ring_atoms.len();
                    let mut ring_bonds = Vec::with_capacity(n);
                    for i in 0..n {
                        let a = ring_atoms[i];
                        let b = ring_atoms[(i + 1) % n];
                        let bond = *bond_map.get(&(a, b))?;
                        ring_bonds.push(bond);
                    }
                    Ring::new(ring_atoms, ring_bonds).ok()
                })
                .collect();
        rings.sort_by_key(|ring| {
            let mut atoms: Vec<usize> = ring.atoms().iter().map(|a| a.index()).collect();
            atoms.sort_unstable();
            (ring.len(), atoms)
        });

        Self::from_rings(
            RingFamily::Induced,
            RingScope::AtomSubset,
            component_atoms.len(),
            rings,
        )
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

    pub fn shared_atoms(&self, a: RingIndex, b: RingIndex) -> Vec<AtomIndex> {
        let (Some(ra), Some(rb)) = (self.ring(a), self.ring(b)) else {
            return Vec::new();
        };
        ra.shared_atoms(rb)
    }

    pub fn shared_bonds(&self, a: RingIndex, b: RingIndex) -> Vec<BondIndex> {
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

    pub fn is_ring_atom(&self, atom: AtomIndex) -> bool {
        self.atom_to_rings.contains_key(&atom)
    }

    pub fn atom_smallest_ring_size(&self, atom: AtomIndex) -> Option<usize> {
        self.atom_to_rings.get(&atom).and_then(|ring_indices| {
            ring_indices
                .iter()
                .map(|i| self.rings[i.index()].len())
                .min()
        })
    }

    pub fn is_ring_bond(&self, bond: BondIndex) -> bool {
        self.bond_to_rings.contains_key(&bond)
    }

    pub fn bond_smallest_ring_size(&self, bond: BondIndex) -> Option<usize> {
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

    pub fn enumerate_builder(&self, builder: &MoleculeBuilder) -> RingSet {
        let mut bond_map: HashMap<(AtomIndex, AtomIndex), BondIndex> = HashMap::new();
        for bond in builder.bond_indices() {
            if let Some((a, b)) = builder.bond_atom_indices(bond) {
                bond_map.insert((a, b), bond);
                bond_map.insert((b, a), bond);
            }
        }
        let full_adj = AtomAdjacency::from_builder(builder);

        match self.family {
            RingFamily::Simple | RingFamily::Induced => {
                let (adj, bcc) = if self.aromatic_only {
                    let pi_atoms: HashSet<AtomIndex> = builder
                        .atom_indices()
                        .filter(|&atom| builder.atom_aromatic_hint(atom))
                        .collect();
                    let adj = full_adj.induced(&pi_atoms);
                    let atoms = adj.atoms();
                    let bcc = molecule_biconnected_components(&atoms, &adj);
                    (adj, bcc)
                } else {
                    let adj = full_adj;
                    let atoms = adj.atoms();
                    let bcc = molecule_biconnected_components(&atoms, &adj);
                    (adj, bcc)
                };
                self.build(&bcc, &adj, &bond_map)
            }
            RingFamily::InducedBenzenoid => {
                let aromatic_carbons: HashSet<AtomIndex> = builder
                    .atom_indices()
                    .filter(|&atom| {
                        builder.atom(atom).is_some_and(|a| {
                            a.element() == Element::C && builder.atom_has_aromatic_candidate(atom)
                        })
                    })
                    .collect();
                let adj = full_adj.induced(&aromatic_carbons);
                let atoms = adj.atoms();
                let bcc = molecule_biconnected_components(&atoms, &adj);
                self.build_induced_benzenoid(&bcc, &adj, &bond_map)
            }
            RingFamily::Mcb | RingFamily::Relevant | RingFamily::Essential => {
                todo!()
            }
        }
    }

    pub fn enumerate_molecule(&self, molecule: &Molecule) -> RingSet {
        let mut bond_map: HashMap<(AtomIndex, AtomIndex), BondIndex> = HashMap::new();
        for bond in molecule.bond_indices() {
            if let Some((a, b)) = molecule.bond_atom_indices(bond) {
                bond_map.insert((a, b), bond);
                bond_map.insert((b, a), bond);
            }
        }
        let full_adj = AtomAdjacency::from_molecule(molecule);

        match self.family {
            RingFamily::Simple | RingFamily::Induced => {
                let (adj, bcc) = if self.aromatic_only {
                    let pi_atoms: HashSet<AtomIndex> = molecule
                        .atom_indices()
                        .filter(|&atom| molecule.atom(atom).is_some_and(|a| a.is_aromatic()))
                        .collect();
                    let adj = full_adj.induced(&pi_atoms);
                    let atoms = adj.atoms();
                    let bcc = molecule_biconnected_components(&atoms, &adj);
                    (adj, bcc)
                } else {
                    let adj = full_adj;
                    let atoms = adj.atoms();
                    let bcc = molecule_biconnected_components(&atoms, &adj);
                    (adj, bcc)
                };
                self.build(&bcc, &adj, &bond_map)
            }
            RingFamily::InducedBenzenoid => {
                let aromatic_carbons: HashSet<AtomIndex> = molecule
                    .atom_indices()
                    .filter(|&atom| {
                        molecule.atom(atom).is_some_and(|a| {
                            a.element() == Element::C && molecule.atom_aromatic_valence(atom) > 0
                        })
                    })
                    .collect();
                let adj = full_adj.induced(&aromatic_carbons);
                let atoms = adj.atoms();
                let bcc = molecule_biconnected_components(&atoms, &adj);
                self.build_induced_benzenoid(&bcc, &adj, &bond_map)
            }
            RingFamily::Mcb | RingFamily::Relevant | RingFamily::Essential => {
                todo!()
            }
        }
    }

    fn build(
        &self,
        bcc: &[Vec<AtomIndex>],
        adj: &AtomAdjacency,
        bond_map: &HashMap<(AtomIndex, AtomIndex), BondIndex>,
    ) -> RingSet {
        let mut all_ring_atoms: Vec<Vec<AtomIndex>> = Vec::new();
        for component in bcc {
            // Domain adapter invariant: each ring-component projection uses
            // deterministic contiguous ids before invoking integer kernels.
            let dense = adj.to_dense_for_atoms(component);
            let mut rings =
                enumerate_simple_cycles(dense.atoms.len(), &dense.adj, self.max_ring_size)
                    .into_iter()
                    .filter(|cycle| {
                        self.family != RingFamily::Induced || is_induced_cycle(cycle, &dense.adj)
                    })
                    .map(|cycle| cycle.into_iter().map(|i| dense.atoms[i]).collect())
                    .collect::<Vec<Vec<AtomIndex>>>();
            rings.truncate(self.max_rings_per_component);
            all_ring_atoms.extend(rings);
        }

        let mut rings: Vec<Ring> = Vec::with_capacity(all_ring_atoms.len());
        for ring_atoms in &all_ring_atoms {
            let n = ring_atoms.len();
            let mut bonds = Vec::with_capacity(n);
            for i in 0..n {
                let a = ring_atoms[i];
                let b = ring_atoms[(i + 1) % n];
                if let Some(&bond) = bond_map.get(&(a, b)) {
                    bonds.push(bond);
                }
            }
            rings.push(
                Ring::new(ring_atoms.clone(), bonds)
                    .expect("enumerated ring must have aligned atom/bond cycle"),
            );
        }

        let is_aromatic_scope = self.aromatic_only || self.family == RingFamily::InducedBenzenoid;

        RingSet::from_rings(
            self.family,
            if is_aromatic_scope {
                RingScope::AromaticSubgraph
            } else {
                RingScope::All
            },
            self.max_ring_size,
            rings,
        )
    }

    fn build_induced_benzenoid(
        &self,
        bcc: &[Vec<AtomIndex>],
        adj: &AtomAdjacency,
        bond_map: &HashMap<(AtomIndex, AtomIndex), BondIndex>,
    ) -> RingSet {
        let mut rings: Vec<Ring> = Vec::new();
        for component in bcc {
            let dense = adj.to_dense_for_atoms(component);
            if dense.atoms.len() < 6 {
                continue;
            }

            let mut component_rings: Vec<Ring> =
                enumerate_simple_cycles(dense.atoms.len(), &dense.adj, 6)
                    .into_iter()
                    .filter(|cycle| cycle.len() == 6)
                    .filter_map(|cycle| {
                        let ring_atoms: Vec<AtomIndex> =
                            cycle.into_iter().map(|i| dense.atoms[i]).collect();
                        let mut ring_bonds = Vec::with_capacity(6);
                        for i in 0..6 {
                            let a = ring_atoms[i];
                            let b = ring_atoms[(i + 1) % 6];
                            let bond = *bond_map.get(&(a, b))?;
                            ring_bonds.push(bond);
                        }
                        Ring::new(ring_atoms, ring_bonds).ok()
                    })
                    .collect();

            component_rings.sort_by_key(|ring| {
                let mut atoms: Vec<usize> = ring.atoms().iter().map(|a| a.index()).collect();
                atoms.sort_unstable();
                atoms
            });
            component_rings.truncate(self.max_rings_per_component);
            rings.extend(component_rings);
        }

        RingSet::from_rings(
            RingFamily::InducedBenzenoid,
            RingScope::AromaticSubgraph,
            6,
            rings,
        )
    }
}

/// Check if a cycle is induced (chordless)
fn is_induced_cycle(cycle: &[usize], adj: &[Vec<usize>]) -> bool {
    let n = cycle.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let adjacent_in_cycle = j == i + 1 || (i == 0 && j == n - 1);
            if adjacent_in_cycle {
                continue;
            }
            if adj[cycle[i]].contains(&cycle[j]) {
                return false;
            }
        }
    }
    true
}

/// Compute the relation between two rings based on shared atoms and bonds.
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

/// Compute the biconnected components of the molecular adjacency list.
fn molecule_biconnected_components(
    atoms: &[AtomIndex],
    adj: &AtomAdjacency,
) -> Vec<Vec<AtomIndex>> {
    let dense = adj.to_dense_for_atoms(atoms);

    biconnected_components(dense.atoms.len(), &dense.adj)
        .into_iter()
        .map(|component| component.into_iter().map(|i| dense.atoms[i]).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::graph_ir::atom::Atom;
    use crate::graph_ir::atom_pattern::AtomPattern;
    use crate::graph_ir::bond_pattern::BondPattern;
    use crate::graph_ir::config::{ResolveConfig, RingEnumerationStrategy};

    fn enumerate(builder: &MoleculeBuilder, max_ring_size: usize) -> RingSet {
        RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: false,
                max_ring_size,
                max_rings_per_component: usize::MAX,
            },
        )
        .enumerate_builder(builder)
    }

    fn enumerate_capped(
        builder: &MoleculeBuilder,
        max_ring_size: usize,
        max_rings_per_component: usize,
    ) -> RingSet {
        RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: false,
                max_ring_size,
                max_rings_per_component,
            },
        )
        .enumerate_builder(builder)
    }

    fn enumerate_molecule(molecule: &Molecule, max_ring_size: usize) -> RingSet {
        RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: false,
                max_ring_size,
                max_rings_per_component: usize::MAX,
            },
        )
        .enumerate_molecule(molecule)
    }

    fn enumerate_aromatic(builder: &MoleculeBuilder, max_ring_size: usize) -> RingSet {
        RingEnumerator::new(
            RingFamily::Simple,
            &RingEnumerationStrategy {
                aromatic_only: true,
                max_ring_size,
                max_rings_per_component: usize::MAX,
            },
        )
        .enumerate_builder(builder)
    }

    #[fixture]
    fn empty_builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    #[fixture]
    fn single_atom_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomPattern::new(Element::C));
        builder
    }

    #[fixture]
    fn chain_builder(#[default(5)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..n - 1 {
            builder.add_bond_unchecked(atoms[i], atoms[i + 1], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    fn ring_builder(#[default(6)] n: usize) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    // bicyclohexyl
    fn disjoint_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..12)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondPattern::new(1));
        }
        for i in 6..12 {
            builder.add_bond_unchecked(atoms[i], atoms[6 + ((i + 1 - 6) % 6)], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    // spiro[3.4]pentane
    fn spiro_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..5)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    // naphthalene
    fn fused_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..10)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let ring1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (a, b) in ring1_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        let ring2_edges = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
        for (a, b) in ring2_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    // bicyclo[1.1.1]pentane
    fn bridged_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..5)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let edges = [(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[rustfmt::skip]
    #[fixture]
    fn multi_spiro_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let edges = [(0, 1), (1, 2), (0, 3), (3, 2), (0, 4), (4, 2), (0, 5), (5, 2)];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[rustfmt::skip]
    #[fixture]
    fn cubane_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..8)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let edges = [
            (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
            (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    // methylcyclohexane
    fn substituted_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let methyl = builder.add_atom(AtomPattern::new(Element::C));
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondPattern::new(1));
        }
        // Attach methyl group to C1 (atom 0)
        builder.add_bond_unchecked(atoms[0], methyl, BondPattern::new(1));
        builder
    }

    #[fixture]
    fn ring_molecule(#[default(6)] n: usize) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..n)
            .map(|_| builder.add_resolved_atom("C#h2#v2".parse::<Atom>().unwrap()))
            .collect();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondPattern::new(1));
        }
        builder
            .build(&ResolveConfig::default())
            .expect("test molecule should build")
    }

    #[fixture]
    fn ring_molecule_rings(#[default(6)] n: usize) -> RingSet {
        let molecule = ring_builder(n);
        enumerate(&molecule, n)
    }

    #[fixture]
    fn disjoint_molecule_rings() -> RingSet {
        let molecule = disjoint_builder();
        enumerate(&molecule, 10)
    }

    #[fixture]
    fn spiro_molecule_rings() -> RingSet {
        let molecule = spiro_builder();
        enumerate(&molecule, 10)
    }

    #[fixture]
    fn fused_molecule_rings() -> RingSet {
        let molecule = fused_builder();
        enumerate(&molecule, 10)
    }

    #[fixture]
    fn bridged_molecule_rings() -> RingSet {
        let molecule = bridged_builder();
        enumerate(&molecule, 5)
    }

    #[fixture]
    fn cubane_molecule_rings() -> RingSet {
        let molecule = cubane_builder();
        enumerate(&molecule, 8)
    }

    #[fixture]
    fn multi_spiro_molecule_rings() -> RingSet {
        let molecule = multi_spiro_builder();
        enumerate(&molecule, 4)
    }

    /// Benzene with aromatic hints, linked to a saturated cyclohexane.
    #[fixture]
    fn aromatic_plus_saturated_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let aromatic: Vec<AtomIndex> = (0..6)
            .map(|_| {
                let idx = builder.add_atom(AtomPattern::new(Element::C));
                builder
                    .set_atom_aromatic_hint(idx, true)
                    .expect("newly added atom index must be valid");
                idx
            })
            .collect();
        for i in 0..6 {
            let bond_idx =
                builder.add_bond_unchecked(aromatic[i], aromatic[(i + 1) % 6], BondPattern::new(1));
            builder.set_bond_aromatic_hint(bond_idx, true);
        }
        let sat: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(sat[i], sat[(i + 1) % 6], BondPattern::new(1));
        }
        builder.add_bond_unchecked(aromatic[0], sat[0], BondPattern::new(1));
        builder
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
    fn test_ring_set_enumerate(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rings = enumerate(&builder, max_ring_size);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::cyclohexane_max_size_0(ring_molecule(6), 0, 0)]
    #[case::cyclohexane_max_size_3(ring_molecule(6), 3, 0)]
    #[case::cyclohexane_max_size_6(ring_molecule(6), 6, 1)]
    #[case::cyclohexane_max_size_8(ring_molecule(6), 8, 1)]
    fn test_ring_set_enumerate_molecule(
        #[case] molecule: Molecule,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rings = enumerate_molecule(&molecule, max_ring_size);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::non_aromatic_only(ring_builder(6), 0)]
    #[case::non_aromatic_and_aromatic(aromatic_plus_saturated_builder(), 1)]
    fn test_ring_set_enumerate_aromatic(#[case] builder: MoleculeBuilder, #[case] expected: usize) {
        let rings = enumerate_aromatic(&builder, 22);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::cubane(cubane_builder(), 8, 3, 3)]
    #[case::cubane(cubane_builder(), 8, 28, 28)]
    fn test_ring_set_enumerate_capped(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] max_rings_per_component: usize,
        #[case] expected: usize,
    ) {
        let rings = enumerate_capped(&builder, max_ring_size, max_rings_per_component);
        assert_eq!(rings.ring_count(), expected);
    }

    #[rstest]
    #[case::ring(ring_molecule_rings(6), vec![RingIndex(0)])]
    #[case::fused(fused_molecule_rings(), vec![RingIndex(0), RingIndex(1), RingIndex(2)])]
    fn test_ring_set_ring_indices(#[case] rings: RingSet, #[case] expected: Vec<RingIndex>) {
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
    fn test_ring_set_ring(
        #[case] rings: RingSet,
        #[case] ring_index: RingIndex,
        #[case] expected: Option<Vec<AtomIndex>>,
    ) {
        let ring = rings.ring(ring_index).map(|v| v.atoms().to_vec());
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
    fn test_ring_set_shared_atoms(
        #[case] rings: RingSet,
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
    fn test_ring_set_shared_bonds(
        #[case] rings: RingSet,
        #[case] ring_index_a: RingIndex,
        #[case] ring_index_b: RingIndex,
        #[case] expected: Vec<BondIndex>,
    ) {
        let mut shared_bonds = rings.shared_bonds(ring_index_a, ring_index_b);
        shared_bonds.sort_unstable();
        assert_eq!(shared_bonds, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::identical(ring_molecule_rings(6), RingIndex(0), RingIndex(0), RingRelation::Identical)]
    #[case::disjoint(disjoint_molecule_rings(), RingIndex(0), RingIndex(1), RingRelation::Disjoint)]
    #[case::spiro(spiro_molecule_rings(), RingIndex(0), RingIndex(1), RingRelation::Spiro)]
    #[case::fused(fused_molecule_rings(), RingIndex(0), RingIndex(1), RingRelation::Fused)]
    #[case::bridged(bridged_molecule_rings(), RingIndex(0), RingIndex(1), RingRelation::Bridged)]
    #[case::multi_spiro(multi_spiro_molecule_rings(), RingIndex(0), RingIndex(5), RingRelation::MultiSpiro)]
    #[case::cubane(cubane_molecule_rings(), RingIndex(0), RingIndex(25), RingRelation::Noncontiguous)]
    #[case::non_existent(ring_molecule_rings(6), RingIndex(0), RingIndex(2), RingRelation::Disjoint)]
    fn test_ring_set_ring_relation(
        #[case] rings: RingSet,
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
    fn test_ring_set_spiro_neighbors(
        #[case] rings: RingSet,
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
    fn test_ring_set_fused_neighbors(
        #[case] rings: RingSet,
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
    fn test_ring_set_bridged_neighbors(
        #[case] rings: RingSet,
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
    fn test_ring_set_fused_components(
        #[case] rings: RingSet,
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
    fn test_ring_set_ring_fused_component(
        #[case] rings: RingSet,
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
    fn test_ring_set_is_ring_atom(
        #[case] builder: MoleculeBuilder,
        #[case] atom_index: usize,
        #[case] expected: bool,
    ) {
        let rings = enumerate(&builder, 10);
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
    fn test_ring_set_atom_smallest_ring_size(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] atom_index: usize,
        #[case] expected_ring_size: Option<usize>,
    ) {
        let rings = enumerate(&builder, max_ring_size);
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
    fn test_ring_set_is_ring_bond(
        #[case] builder: MoleculeBuilder,
        #[case] bond_index: usize,
        #[case] expected: bool,
    ) {
        let rings = enumerate(&builder, 10);
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
    fn test_ring_set_bond_smallest_ring_size(
        #[case] builder: MoleculeBuilder,
        #[case] max_ring_size: usize,
        #[case] bond_index: usize,
        #[case] expected_ring_size: Option<usize>,
    ) {
        let rings = enumerate(&builder, max_ring_size);
        let bond = BondIndex::new(bond_index);
        assert_eq!(rings.bond_smallest_ring_size(bond), expected_ring_size);
    }

    #[test]
    fn test_ring_set_ring_graph() {
        let rings = fused_molecule_rings();
        let ring_graph = rings.ring_graph();
        assert!(!ring_graph.edges().is_empty());
        assert!(ring_graph
            .edges()
            .iter()
            .any(|e| matches!(e.relation, RingRelation::Fused)));
    }
}
