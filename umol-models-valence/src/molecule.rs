//! Molecular representation as valence graph

use std::collections::{HashMap, HashSet};
use std::fmt::{self, Display};

use indexmap::IndexMap;
use petgraph::graph::NodeIndex;
use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use umol::error::DataError;
use umol::Result;

use crate::{Atom, AtomBuilder, Bond, BondBuilder};

/// Internal atom and bond indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

/// Graph model of atoms and bonds, with valence constraints
#[derive(Debug, Clone)]
pub struct Molecule {
    data: StableGraph<Atom, Bond, Undirected, usize>,
}

impl Molecule {
    pub fn atom_count(&self) -> usize {
        self.data.node_count()
    }

    pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.data.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.data.node_weight(index)
    }

    pub fn bond_count(&self) -> usize {
        self.data.edge_count()
    }

    pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edge_weights()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.data.edge_weight(index)
    }

    // TODO: Review naming
    pub fn bond_atoms(
        &self,
        index: BondIndex,
    ) -> Option<(&Atom, &Atom)> {
        self.data.edge_endpoints(index).map(|(a, b)| {
            (
                self.data.node_weight(a).unwrap(),
                self.data.node_weight(b).unwrap(),
            )
        })
    }

    // TODO: Review naming
    pub fn bond_atom_indices(
        &self,
        index: BondIndex,
    ) -> Option<(AtomIndex, AtomIndex)> {
        self.data.edge_endpoints(index)
    }

    // TODO: Review naming
    pub fn atom_bonds<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edges(index).map(|e| e.weight())
    }

    // TODO: Review naming
    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.data.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbors<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.data
            .neighbors(index)
            .map(|n| self.data.node_weight(n).unwrap())
    }

    pub fn atom_neighbor_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AtomIndex> + 'graph {
        self.data.neighbors(index)
    }
    // TODO: Add methods for converting some atoms/bonds to builders
}

impl Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "Molecule with {} atoms and {} bonds:",
            self.atom_count(),
            self.bond_count()
        )?;
        for (i, atom) in self.atoms().enumerate() {
            writeln!(f, "  Atom {}: {:?}", i, atom)?;
        }
        for (i, bond) in self.bonds().enumerate() {
            if let Some((a, b)) = self.bond_atoms(BondIndex::new(i)) {
                writeln!(
                    f,
                    "  Bond {}: {} between atoms {:?} and {:?}",
                    i, bond, a, b
                )?;
            }
        }
        Ok(())
    }
}

/// Builder type for ValenceGraphs, allowing incremental construction and strict validation.
pub struct MoleculeBuilder {
    atom_builders: HashMap<usize, AtomBuilder>,
    bond_builders: HashMap<(usize, usize), BondBuilder>,
}

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            atom_builders: HashMap::new(),
            bond_builders: HashMap::new(),
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            atom_builders: HashMap::with_capacity(atom_capacity),
            bond_builders: HashMap::with_capacity(bond_capacity),
        }
    }

    pub fn create_atom<A: Into<AtomBuilder>>(&mut self, atom: A) -> (usize, &mut AtomBuilder) {
        let builder = atom.into();
        let idx = self.atom_builders.len();
        self.atom_builders.insert(idx, builder);
        (idx, self.atom_builders.get_mut(&idx).unwrap())
    }

    pub fn create_atoms<A: Into<AtomBuilder>>(
        &mut self,
        atoms: impl IntoIterator<Item = A>,
    ) -> impl Iterator<Item = usize> {
        let builders_iter = atoms.into_iter().map(|atom| atom.into());
        let (lbound, _ubound) = builders_iter.size_hint();
        let offset = self.atom_builders.len();
        let indices = builders_iter.enumerate().fold(
            Vec::with_capacity(lbound),
            |mut acc, (idx, builder)| {
                self.atom_builders.insert(offset + idx, builder);
                acc.push(offset + idx);
                acc
            },
        );
        indices.into_iter()
    }

    pub fn add_atom<A: Into<AtomBuilder>>(
        &mut self,
        idx: usize,
        atom: A,
    ) -> Result<(usize, &mut AtomBuilder)> {
        let builder = atom.into();
        if self.atom_builders.contains_key(&idx) {
            return Err(DataError::DuplicateAtomIndex(idx).into());
        }
        self.atom_builders.insert(idx, builder);
        Ok((idx, self.atom_builders.get_mut(&idx).unwrap()))
    }

    pub fn add_atoms<A: Into<AtomBuilder>>(
        &mut self,
        atoms: impl IntoIterator<Item = (usize, A)>,
    ) -> Result<impl Iterator<Item = usize>> {
        // Collect input into a temporary Vec and validate
        let staged_atoms: Vec<(usize, AtomBuilder)> = atoms
            .into_iter()
            .map(|(idx, atom)| (idx, atom.into()))
            .collect();
        let mut seen_indices = HashSet::with_capacity(staged_atoms.len());
        for (idx, _) in &staged_atoms {
            if !seen_indices.insert(*idx) {
                return Err(DataError::DuplicateAtomIndex(*idx).into());
            }
            if self.atom_builders.contains_key(idx) {
                return Err(DataError::DuplicateAtomIndex(*idx).into());
            }
        }

        // Commit changes
        let mut indices = Vec::with_capacity(staged_atoms.len());
        for (idx, builder) in staged_atoms {
            self.atom_builders.insert(idx, builder);
            indices.push(idx);
        }

        Ok(indices.into_iter())
    }

    pub fn add_bond<B: Into<BondBuilder>>(
        &mut self,
        idx1: usize,
        idx2: usize,
        bond: B,
    ) -> Result<(usize, usize, &mut BondBuilder)> {
        if idx1 == idx2 {
            return Err(DataError::LoopBond(idx1).into());
        }
        if !self.atom_builders.contains_key(&idx1) {
            return Err(DataError::MissingAtomIndex(idx1).into());
        }
        if !self.atom_builders.contains_key(&idx2) {
            return Err(DataError::MissingAtomIndex(idx2).into());
        }
        let builder = bond.into();
        if self.bond_builders.contains_key(&(idx1, idx2)) {
            return Err(DataError::DuplicateBondIndex(idx1, idx2).into());
        }
        self.bond_builders.insert((idx1, idx2), builder);
        Ok((
            idx1,
            idx2,
            self.bond_builders.get_mut(&(idx1, idx2)).unwrap(),
        ))
    }

    pub fn add_bonds<B: Into<BondBuilder>>(
        &mut self,
        bonds: impl IntoIterator<Item = (usize, usize, B)>,
    ) -> Result<impl Iterator<Item = (usize, usize)>> {
        // Helper for canonical keys
        let canonical_bond_key = |idx1: usize, idx2: usize| (idx1.min(idx2), idx1.max(idx2));

        // Collect input into a temporary Vec and validate
        let staged_bonds: Vec<(usize, usize, BondBuilder)> = bonds
            .into_iter()
            .map(|(idx1, idx2, bond)| (idx1, idx2, bond.into()))
            .collect();
        let mut seen_keys = HashSet::with_capacity(staged_bonds.len());
        for (idx1, idx2, _builder) in &staged_bonds {
            if idx1 == idx2 {
                return Err(DataError::LoopBond(*idx1).into());
            }
            if !self.atom_builders.contains_key(idx1) {
                return Err(DataError::MissingAtomIndex(*idx1).into());
            }
            if !self.atom_builders.contains_key(idx2) {
                return Err(DataError::MissingAtomIndex(*idx2).into());
            }

            let key = canonical_bond_key(*idx1, *idx2);
            if !seen_keys.insert(key) {
                return Err(DataError::DuplicateBondIndex(*idx1, *idx2).into());
            }
            if self.bond_builders.contains_key(&key) {
                return Err(DataError::DuplicateBondIndex(*idx1, *idx2).into());
            }
        }

        // Commit changes
        let mut indices = Vec::with_capacity(staged_bonds.len());
        for (idx1, idx2, builder) in staged_bonds {
            let key = canonical_bond_key(idx1, idx2);
            self.bond_builders.insert(key, builder);
            indices.push((idx1, idx2));
        }

        Ok(indices.into_iter())
    }

    pub fn build(self) -> Result<Molecule> {
        let mut atom_builders = self.atom_builders;
        let bond_builders = self.bond_builders;
        let mut built_bonds = HashMap::with_capacity(bond_builders.len());

        // Build bonds and update atom valences
        for (key @ (idx1, idx2), bond_builder) in bond_builders {
            let bond = bond_builder.build()?;
            let valence = bond.order().value();

            let atom1_builder = atom_builders
                .get_mut(&idx1)
                .expect("Atom builder missing for idx1 during build");
            atom1_builder.update_valence(|v| v + valence);

            let atom2_builder = atom_builders
                .get_mut(&idx2)
                .expect("Atom builder missing for idx2 during build");
            atom2_builder.update_valence(|v| v + valence);
            built_bonds.insert(key, bond);
        }

        // Build atoms
        let mut built_atoms = IndexMap::with_capacity(atom_builders.len());
        for (idx, atom_builder) in atom_builders {
            let atom = atom_builder.build()?;
            built_atoms.insert(idx, atom);
        }

        // Initialize graph and add atoms
        let mut graph = StableGraph::with_capacity(built_atoms.len(), built_bonds.len());
        let mut atom_indices = HashMap::with_capacity(built_atoms.len());
        for (idx, atom) in built_atoms {
            let node_index = graph.add_node(atom);
            atom_indices.insert(idx, node_index);
        }

        // Add bonds to graph
        for ((idx1, idx2), bond) in built_bonds {
            let node1 = *atom_indices
                .get(&idx1)
                .expect("Node index map missing mapping for idx1 in step 5");
            let node2 = *atom_indices
                .get(&idx2)
                .expect("Node index map missing mapping for idx2 in step 5");
            graph.add_edge(node1, node2, bond);
        }

        Ok(Molecule { data: graph })
    }
}

#[cfg(test)]
mod tests {
    use umol::Error;
    use umol_data::Element;

    use super::*;
    use crate::{AtomBuilder, BondBuilder, BondOrder};

    #[test]
    fn test_builder_new() {
        let builder = MoleculeBuilder::new();
        assert!(builder.atom_builders.is_empty());
        assert!(builder.bond_builders.is_empty());
    }

    #[test]
    fn test_builder_with_capacity() {
        let builder = MoleculeBuilder::with_capacity(10, 20);
        assert!(builder.atom_builders.is_empty());
        assert!(builder.bond_builders.is_empty());
    }

    #[test]
    fn test_builder_create_atom() {
        let mut builder = MoleculeBuilder::new();
        let (idx0, atom0_builder) = builder.create_atom(Element::C);
        atom0_builder.set_charge(1);
        assert_eq!(idx0, 0);
        assert!(builder.atom_builders.contains_key(&0));
        assert_eq!(builder.atom_builders.len(), 1);
        assert_eq!(builder.atom_builders[&0].charge(), Some(1));

        let (idx1, _) = builder.create_atom(Element::O);
        assert_eq!(idx1, 1);
        assert_eq!(builder.atom_builders.len(), 2);
        assert!(builder.atom_builders.contains_key(&1));
    }

    #[test]
    fn test_builder_create_atoms() {
        let mut builder = MoleculeBuilder::new();
        let atoms_to_add = vec![AtomBuilder::new(Element::H), AtomBuilder::new(Element::H)];
        let indices: Vec<_> = builder.create_atoms(atoms_to_add).collect();

        assert_eq!(indices, vec![0, 1]);
        assert_eq!(builder.atom_builders.len(), 2);
        assert!(builder.atom_builders.contains_key(&0));
        assert_eq!(builder.atom_builders[&0].element(), Element::H);
        assert!(builder.atom_builders.contains_key(&1));
        assert_eq!(builder.atom_builders[&1].element(), Element::H);

        // Add more atoms to existing builder
        let more_atoms = vec![AtomBuilder::new(Element::O)];
        let indices2: Vec<_> = builder.create_atoms(more_atoms).collect();
        assert_eq!(indices2, vec![2]);
        assert_eq!(builder.atom_builders.len(), 3);
        assert!(builder.atom_builders.contains_key(&2));
        assert_eq!(builder.atom_builders[&2].element(), Element::O);
    }

    #[test]
    fn test_builder_add_atom() {
        let mut builder = MoleculeBuilder::new();
        let (idx5, atom5_builder) = builder.add_atom(5, Element::N).unwrap();
        atom5_builder.set_lone_pairs(1);
        assert_eq!(idx5, 5);
        assert_eq!(builder.atom_builders.len(), 1);
        assert!(builder.atom_builders.contains_key(&5));
        assert_eq!(builder.atom_builders[&5].lone_pairs(), Some(1));

        let result = builder.add_atom(5, Element::C);
        assert!(result.is_err());
        match result {
            Err(Error::Data(DataError::DuplicateAtomIndex(idx))) => assert_eq!(idx, 5),
            _ => panic!("Expected DuplicateAtomIndex error"),
        }
    }

    #[test]
    fn test_builder_add_atoms() {
        let mut builder = MoleculeBuilder::new();
        builder.create_atom(Element::C); // index 0

        let atoms_to_add = vec![
            (2, AtomBuilder::new(Element::O)),
            (1, AtomBuilder::new(Element::N)),
        ];
        let result = builder.add_atoms(atoms_to_add);
        assert!(result.is_ok());
        let indices: Vec<_> = result.unwrap().collect();
        assert_eq!(indices.len(), 2);
        assert!(indices.contains(&1));
        assert!(indices.contains(&2));

        assert_eq!(builder.atom_builders.len(), 3);
        assert!(builder.atom_builders.contains_key(&1));
        assert_eq!(builder.atom_builders[&1].element(), Element::N);
        assert!(builder.atom_builders.contains_key(&2));
        assert_eq!(builder.atom_builders[&2].element(), Element::O);

        // Test duplicate index within input
        let current_atom_count = builder.atom_builders.len(); // Should be 3
        let duplicate_indices = vec![
            (3, AtomBuilder::new(Element::H)),
            (3, AtomBuilder::new(Element::H)),
        ];
        let result_dup = builder.add_atoms(duplicate_indices);
        assert!(result_dup.is_err());
        match result_dup {
            Err(Error::Data(DataError::DuplicateAtomIndex(idx))) => assert_eq!(idx, 3),
            _ => panic!("Expected DuplicateAtomIndex error"),
        }
        // Atomicity check: length should be unchanged after failed batch add
        assert_eq!(builder.atom_builders.len(), current_atom_count);

        // Test duplicate index conflicting with existing
        let current_atom_count_before_conflict = builder.atom_builders.len(); // Should still be 3
        let conflict_indices = vec![
            (4, AtomBuilder::new(Element::F)),
            (0, AtomBuilder::new(Element::P)), // Conflicts with existing index 0
        ];
        let result_conf = builder.add_atoms(conflict_indices);
        assert!(result_conf.is_err());
        match result_conf {
            Err(Error::Data(DataError::DuplicateAtomIndex(idx))) => assert_eq!(idx, 0),
            _ => panic!("Expected DuplicateAtomIndex error for existing index"),
        }
        // Atomicity check: length should be unchanged after failed batch add
        assert_eq!(
            builder.atom_builders.len(),
            current_atom_count_before_conflict
        );
    }

    #[test]
    fn test_builder_add_bond() {
        let mut builder = MoleculeBuilder::new();
        builder.create_atom(Element::C); // 0
        builder.create_atom(Element::O); // 1
        builder.create_atom(Element::N); // 2

        let (idx1, idx2, bond_builder) = builder.add_bond(0, 1, BondOrder::Single).unwrap();
        bond_builder.set_order(BondOrder::Double);
        assert_eq!((idx1, idx2), (0, 1));
        assert_eq!(builder.bond_builders.len(), 1);
        assert!(builder.bond_builders.contains_key(&(0, 1)));
        assert_eq!(builder.bond_builders[&(0, 1)].order(), BondOrder::Double);

        // Test loop bond
        let res_loop = builder.add_bond(2, 2, BondOrder::Single);
        assert!(matches!(res_loop, Err(Error::Data(DataError::LoopBond(2)))));

        // Test missing atom 1
        let res_missing1 = builder.add_bond(99, 1, BondOrder::Single);
        assert!(matches!(
            res_missing1,
            Err(Error::Data(DataError::MissingAtomIndex(99)))
        ));

        // Test missing atom 2
        let res_missing2 = builder.add_bond(0, 98, BondOrder::Single);
        assert!(matches!(
            res_missing2,
            Err(Error::Data(DataError::MissingAtomIndex(98)))
        ));

        // Test duplicate bond
        let res_dup = builder.add_bond(0, 1, BondOrder::Triple);
        assert!(matches!(
            res_dup,
            Err(Error::Data(DataError::DuplicateBondIndex(0, 1)))
        ));
    }

    #[test]
    fn test_builder_add_bonds() {
        let mut builder = MoleculeBuilder::new();
        builder.create_atom(Element::C); // 0
        builder.create_atom(Element::O); // 1
        builder.create_atom(Element::N); // 2
        builder.create_atom(Element::H); // 3

        let bonds_to_add = vec![
            (0, 1, BondBuilder::new(BondOrder::Single)),
            (1, 2, BondBuilder::new(BondOrder::Double)),
        ];
        let result = builder.add_bonds(bonds_to_add);
        assert!(result.is_ok());
        let indices: Vec<_> = result.unwrap().collect();
        assert_eq!(indices.len(), 2);
        assert!(indices.contains(&(0, 1)));
        assert!(indices.contains(&(1, 2)));

        assert_eq!(builder.bond_builders.len(), 2);
        assert!(builder.bond_builders.contains_key(&(0, 1)));
        assert!(builder.bond_builders.contains_key(&(1, 2)));
        assert_eq!(builder.bond_builders[&(0, 1)].order(), BondOrder::Single);
        assert_eq!(builder.bond_builders[&(1, 2)].order(), BondOrder::Double);

        // Test duplicate within input
        let current_bond_count = builder.bond_builders.len(); // Should be 2
        let bonds_dup = vec![
            (2, 3, BondBuilder::new(BondOrder::Single)),
            (2, 3, BondBuilder::new(BondOrder::Triple)), // Duplicate key (canonical)
        ];
        let res_dup = builder.add_bonds(bonds_dup);
        assert!(matches!(
            res_dup,
            Err(Error::Data(DataError::DuplicateBondIndex(2, 3)))
        ));
        // Atomicity check: length should be unchanged after failed batch add
        assert_eq!(builder.bond_builders.len(), current_bond_count);

        // Test conflict with existing
        let current_bond_count_before_conflict = builder.bond_builders.len(); // Should still be 2
        let bonds_conf = vec![
            (3, 0, BondBuilder::new(BondOrder::Single)),
            (0, 1, BondBuilder::new(BondOrder::Single)), // Conflicts with existing (0,1)
        ];
        let res_conf = builder.add_bonds(bonds_conf);
        assert!(matches!(
            res_conf,
            Err(Error::Data(DataError::DuplicateBondIndex(0, 1)))
        ));
        // Atomicity check: length should be unchanged after failed batch add
        assert_eq!(
            builder.bond_builders.len(),
            current_bond_count_before_conflict
        );

        // Test missing atom
        let bonds_missing = vec![(0, 99, BondBuilder::new(BondOrder::Single))];
        let res_missing = builder.add_bonds(bonds_missing);
        assert!(matches!(
            res_missing,
            Err(Error::Data(DataError::MissingAtomIndex(99)))
        ));

        // Test loop bond
        let bonds_loop = vec![(3, 3, BondBuilder::new(BondOrder::Single))];
        let res_loop = builder.add_bonds(bonds_loop);
        assert!(matches!(res_loop, Err(Error::Data(DataError::LoopBond(3)))));

        // Test Canonical Key Handling
        let current_bond_count_before_canonical_test = builder.bond_builders.len(); // Should be 2

        // Test conflict with existing bond (0,1) but adding (1,0)
        let bond_conf_reversed = vec![
            (1, 0, BondBuilder::new(BondOrder::Single)), // Conflicts with existing (0,1) via canonical key
        ];
        let res_conf_rev = builder.add_bonds(bond_conf_reversed);
        assert!(matches!(
            res_conf_rev,
            Err(Error::Data(DataError::DuplicateBondIndex(1, 0)))
                | Err(Error::Data(DataError::DuplicateBondIndex(0, 1)))
        ));
        assert_eq!(
            builder.bond_builders.len(),
            current_bond_count_before_canonical_test
        );

        // Test duplicate within batch using reversed indices
        builder.create_atom(Element::F); // Atom 4
        builder.create_atom(Element::Cl); // Atom 5
        let bonds_dup_reversed = vec![
            (4, 5, BondBuilder::new(BondOrder::Single)),
            (5, 4, BondBuilder::new(BondOrder::Double)), // Duplicate via canonical key
        ];
        let res_dup_rev = builder.add_bonds(bonds_dup_reversed);
        assert!(matches!(
            res_dup_rev,
            Err(Error::Data(DataError::DuplicateBondIndex(5, 4)))
                | Err(Error::Data(DataError::DuplicateBondIndex(4, 5)))
        ));
        assert_eq!(
            builder.bond_builders.len(),
            current_bond_count_before_canonical_test
        ); // Should still be 2, as (4,5)/(5,4) batch failed
    }

    #[test]
    fn test_builder_build_success() {
        let mut builder = MoleculeBuilder::new();
        let h1 = builder.create_atom(Element::H).0;
        let o = builder.create_atom(Element::O).0;
        let h2 = builder.create_atom(Element::H).0;

        builder.add_bond(h1, o, BondOrder::Single).unwrap();
        builder.add_bond(o, h2, BondOrder::Single).unwrap();

        let result = builder.build();
        assert!(result.is_ok());

        let molecule = result.unwrap();
        assert_eq!(molecule.atom_count(), 3);
        assert_eq!(molecule.bond_count(), 2);

        let oxygen = molecule
            .atoms()
            .find(|a| a.element() == Element::O)
            .unwrap();
        assert_eq!(oxygen.implicit_hydrogens(), 0);
        assert_eq!(oxygen.valence(), 2);
    }

    #[test]
    fn test_builder_build_fail_atom_valence() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.create_atom(Element::C).0;
        let h1 = builder.create_atom(Element::H).0;
        let h2 = builder.create_atom(Element::H).0;
        let h3 = builder.create_atom(Element::H).0;
        let h4 = builder.create_atom(Element::H).0;
        let h5 = builder.create_atom(Element::H).0; // The problematic 5th hydrogen

        builder.add_bond(c, h1, BondOrder::Single).unwrap();
        builder.add_bond(c, h2, BondOrder::Single).unwrap();
        builder.add_bond(c, h3, BondOrder::Single).unwrap();
        builder.add_bond(c, h4, BondOrder::Single).unwrap();
        builder.add_bond(c, h5, BondOrder::Single).unwrap(); // 5th bond

        let result = builder.build();
        assert!(result.is_err());

        // Check for an appropriate error, likely stemming from AtomBuilder::build validation
        match result {
            Err(Error::Data(DataError::NoAtomSpec(msg))) => {
                // Check if the error message relates to the Carbon atom
                assert!(msg.contains("element: C"));
                assert!(msg.contains("valence: Some(5)")); // AtomBuilder accumulated valence 5
            }
            Err(e) => panic!(
                "Expected NoAtomSpec error due to invalid valence, got {:?}",
                e
            ),
            Ok(_) => panic!("Build succeeded unexpectedly for invalid molecule"),
        }
    }
}
