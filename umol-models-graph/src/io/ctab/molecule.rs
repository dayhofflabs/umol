//! Molecule type for CTab format.

use crate::io::ctab::atom::{Atom, AtomStandard, AtomSymbol};
use crate::io::ctab::bond::{Bond, BondStandard};
use crate::io::ctab::sgroup::SGroup;
use umol_data::{e, Element};

use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::graph6::get_graph6_representation;
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// Type aliases for the node and edge indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

fn element_symbol_key(element: Element) -> [u8; 2] {
    let symbol = element.symbol();
    let bytes = symbol.as_bytes();
    [
        bytes.get(0).copied().unwrap_or(0),
        bytes.get(1).copied().unwrap_or(0),
    ]
}

/// Format sum formula according to Hill notation
/// 
/// Hill notation: C first, H second, then other elements alphabetically by symbol
fn format_sum_formula(
    c_count: usize,
    h_count: usize,
    atom_counts: BTreeMap<[u8; 2], (Element, usize)>,
    charge: i8,
) -> String {
    let mut sum_formula = String::new();
    
    // Carbon first
    if c_count > 1 {
        sum_formula.push_str(&format!("C{}", c_count));
    } else if c_count == 1 {
        sum_formula.push_str("C");
    }
    
    // Hydrogen second
    if h_count > 1 {
        sum_formula.push_str(&format!("H{}", h_count));
    } else if h_count == 1 {
        sum_formula.push_str("H");
    }
    
    // Other elements alphabetically by symbol (BTreeMap with [u8; 2] keys maintains order)
    for (_, (element, count)) in atom_counts {
        if count > 1 {
            sum_formula.push_str(&format!("{}{}", element, count));
        } else {
            sum_formula.push_str(&element.to_string());
        }
    }
    
    // Charge at the end
    if charge != 0 {
        sum_formula.push_str(&format!("{:+}", charge));
    }
    
    sum_formula
}

/// Graph-based molecule representation with full MOL file semantics (including queries)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
}

/// Graph-based molecule representation for standard (non-query) molecules only
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MoleculeStandard {
    pub graph: StableGraph<AtomStandard, BondStandard, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
}

impl Molecule {
    /// Create empty molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<Atom, Bond, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
        }
    }

    /// Get number of atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get number of bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Get molecule-level property by key
    pub fn property(&self, key: &str) -> Option<&String> {
        self.properties.get(key)
    }

    /// Set molecule-level property by key
    pub fn set_property(&mut self, key: String, value: String) {
        self.properties.insert(key, value);
    }

    /// Get molecule-level properties as hashmap
    pub fn properties(&self) -> &HashMap<String, String> {
        &self.properties
    }

    /// Get mutable reference to molecule-level properties map
    pub fn properties_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.properties
    }

    /// Add atom to the molecule and update index mappings
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: Bond) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get iterator over atom indices
    pub fn atom_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.node_indices().map(|i| i.index())
    }

    /// Get iterator over atoms
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> + '_ {
        self.graph.node_weights()
    }

    /// Get immutable reference to atom by index
    pub fn atom(&self, idx: usize) -> Option<&Atom> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut Atom> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get iterator over bond indices
    pub fn bond_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.edge_indices().map(|i| i.index())
    }

    /// Get iterator over bonds
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> + '_ {
        self.graph.edge_weights()
    }

    /// Get immutable reference to bond by index
    pub fn bond(&self, idx: usize) -> Option<&Bond> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut Bond> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }

    /// Get iterator over neighbor atom indices for atom index
    pub fn neighbors(&self, idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.graph.neighbors(AtomIndex::new(idx)).map(|i| i.index())
    }

    /// Get iterator over sgroup indices
    pub fn sgroup_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.sgroups.keys().copied()
    }

    /// Get iterator over sgroups
    pub fn sgroups(&self) -> impl Iterator<Item = &SGroup> + '_ {
        self.sgroups.values()
    }

    /// Get immutable reference to sgroup by index
    pub fn sgroup(&self, idx: usize) -> Option<&SGroup> {
        self.sgroups.get(&idx)
    }

    /// Get mutable reference to sgroup by index
    pub fn sgroup_mut(&mut self, idx: usize) -> Option<&mut SGroup> {
        self.sgroups.get_mut(&idx)
    }

    /// Add sgroup to the molecule
    pub fn add_sgroup(&mut self, sgroup_index: usize, sgroup: SGroup) {
        self.sgroups.insert(sgroup_index, sgroup);
    }

    /// Get sum formula
    pub fn sum_formula(&self) -> String {
        let mut atom_counts = BTreeMap::new();
        let mut c_count = 0;
        let mut h_count = 0;
        let mut charge = 0;
        for atom in self.atoms() {
            match atom.symbol {
                AtomSymbol::Element(element) => {
                    if element == e!(C) {
                        c_count += 1;
                    } else if element == e!(H) {
                        h_count += 1;
                    } else {
                        let key = element_symbol_key(element);
                        atom_counts.entry(key).or_insert((element, 0)).1 += 1;
                    }
                }
                AtomSymbol::NamedIsotope(_) => {
                    // Only H is allowed in named isotopes
                    h_count += 1;
                }
                _ => {}
            }
            charge += atom.charge;
        }
        format_sum_formula(c_count, h_count, atom_counts, charge)
    }

    /// Get graph6 representation of the molecule (excluding atom and bond labels)
    pub fn graph6(&self) -> String {
        get_graph6_representation(&self.graph)
    }
}

impl MoleculeStandard {
    /// Create empty standard molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<AtomStandard, BondStandard, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
        }
    }

    /// Get number of atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get number of bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Add atom to the molecule and update index mappings
    pub fn add_atom(&mut self, atom: AtomStandard) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: BondStandard) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get iterator over atom indices
    pub fn atom_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.node_indices().map(|i| i.index())
    }

    /// Get iterator over atoms
    pub fn atoms(&self) -> impl Iterator<Item = &AtomStandard> + '_ {
        self.graph.node_weights()
    }

    /// Get immutable reference to atom by index
    pub fn atom(&self, idx: usize) -> Option<&AtomStandard> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut AtomStandard> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get iterator over bond indices
    pub fn bond_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.graph.edge_indices().map(|i| i.index())
    }

    /// Get iterator over bonds
    pub fn bonds(&self) -> impl Iterator<Item = &BondStandard> + '_ {
        self.graph.edge_weights()
    }

    /// Get immutable reference to bond by index
    pub fn bond(&self, idx: usize) -> Option<&BondStandard> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut BondStandard> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }

    /// Get sum formula
    pub fn sum_formula(&self) -> String {
        let mut atom_counts = BTreeMap::new();
        let mut charge = 0;
        let mut c_count = 0;
        let mut h_count = 0;
        for atom in self.atoms() {
            match atom.element {
                e!(C) => {
                    c_count += 1;
                }
                e!(H) => {
                    h_count += 1;
                }
                _ => {
                    let key = element_symbol_key(atom.element);
                    atom_counts.entry(key).or_insert((atom.element, 0)).1 += 1;
                }
            }
            charge += atom.charge;
        }
        format_sum_formula(c_count, h_count, atom_counts, charge)
    }

    /// Get graph6 representation of the molecule (excluding atom and bond labels)
    pub fn graph6(&self) -> String {
        get_graph6_representation(&self.graph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_standard_serialize() {
        let graph =
            StableGraph::<AtomStandard, BondStandard, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();

        let molecule = MoleculeStandard {
            graph,
            sgroups,
            properties,
        };

        let yaml =
            serde_yaml::to_string(&molecule).expect("Failed to serialize MoleculeStandard to YAML");
        let deserialized: MoleculeStandard =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize MoleculeStandard from YAML");
        assert_eq!(molecule.properties, deserialized.properties);
    }

    #[test]
    fn test_molecule_standard_sum_formula_glycine() {
        use crate::io::ctab::atom::AtomStandard;
        use crate::io::ctab::bond::{BondStandard, BondType};
        use umol_data::{e, Element};

        // Glycine: NH2-CH2-COOH (C2H5NO2)
        let mut molecule = MoleculeStandard::new();
        molecule.add_atom(AtomStandard::new(e!(N)));  // 0: N
        molecule.add_atom(AtomStandard::new(e!(H)));  // 1: H
        molecule.add_atom(AtomStandard::new(e!(H)));  // 2: H
        molecule.add_atom(AtomStandard::new(e!(C)));  // 3: C
        molecule.add_atom(AtomStandard::new(e!(H)));  // 4: H
        molecule.add_atom(AtomStandard::new(e!(H)));  // 5: H
        molecule.add_atom(AtomStandard::new(e!(C)));  // 6: C
        molecule.add_atom(AtomStandard::new(e!(O)));  // 7: O
        molecule.add_atom(AtomStandard::new(e!(O)));  // 8: O
        molecule.add_atom(AtomStandard::new(e!(H)));  // 9: H
        
        // Bonds
        molecule.add_bond(0, 1, BondStandard::new(BondType::Single));  // N-H
        molecule.add_bond(0, 2, BondStandard::new(BondType::Single));  // N-H
        molecule.add_bond(0, 3, BondStandard::new(BondType::Single));  // N-C
        molecule.add_bond(3, 4, BondStandard::new(BondType::Single));  // C-H
        molecule.add_bond(3, 5, BondStandard::new(BondType::Single));  // C-H
        molecule.add_bond(3, 6, BondStandard::new(BondType::Single));  // C-C
        molecule.add_bond(6, 7, BondStandard::new(BondType::Double)); // C=O
        molecule.add_bond(6, 8, BondStandard::new(BondType::Single)); // C-O
        molecule.add_bond(8, 9, BondStandard::new(BondType::Single)); // O-H
        
        assert_eq!(molecule.sum_formula(), "C2H5NO2");
    }

    #[test]
    fn test_molecule_serialize() {
        let graph = StableGraph::<Atom, Bond, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();
        let molecule = Molecule {
            graph,
            sgroups,
            properties,
        };

        let yaml = serde_yaml::to_string(&molecule).expect("Failed to serialize Molecule to YAML");
        let deserialized: Molecule =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Molecule from YAML");
        assert_eq!(molecule.properties, deserialized.properties);
    }

    #[test]
    fn test_molecule_sum_formula() {
        use crate::io::ctab::atom::AtomSymbol;
        use crate::io::ctab::bond::BondType;
        use umol_data::{e, Element};

        let mut molecule = Molecule::new();
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(O))));
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));
        molecule.add_bond(0, 1, Bond::new(BondType::Single));
        molecule.add_bond(0, 2, Bond::new(BondType::Single));
        assert_eq!(molecule.sum_formula(), "H2O");
    }

    #[test]
    fn test_molecule_sum_formula_glycine() {
        use crate::io::ctab::atom::AtomSymbol;
        use crate::io::ctab::bond::BondType;
        use umol_data::{e, Element};

        // Glycine: NH2-CH2-COOH (C2H5NO2)
        let mut molecule = Molecule::new();
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(N))));  // 0: N
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));  // 1: H
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));  // 2: H
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));  // 3: C
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));  // 4: H
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));  // 5: H
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(C))));  // 6: C
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(O))));  // 7: O
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(O))));  // 8: O
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));  // 9: H
        
        // Bonds
        molecule.add_bond(0, 1, Bond::new(BondType::Single));  // N-H
        molecule.add_bond(0, 2, Bond::new(BondType::Single));  // N-H
        molecule.add_bond(0, 3, Bond::new(BondType::Single));  // N-C
        molecule.add_bond(3, 4, Bond::new(BondType::Single));  // C-H
        molecule.add_bond(3, 5, Bond::new(BondType::Single));  // C-H
        molecule.add_bond(3, 6, Bond::new(BondType::Single));  // C-C
        molecule.add_bond(6, 7, Bond::new(BondType::Double)); // C=O
        molecule.add_bond(6, 8, Bond::new(BondType::Single)); // C-O
        molecule.add_bond(8, 9, Bond::new(BondType::Single)); // O-H
        
        assert_eq!(molecule.sum_formula(), "C2H5NO2");
    }

    #[test]
    fn test_format_sum_formula() {
        use umol_data::{e, Element};
        use std::collections::BTreeMap;

        // Empty molecule
        assert_eq!(format_sum_formula(0, 0, BTreeMap::new(), 0), "");

        // Water (H2O)
        let mut atom_counts = BTreeMap::new();
        atom_counts.insert(element_symbol_key(e!(O)), (e!(O), 1));
        assert_eq!(format_sum_formula(0, 2, atom_counts, 0), "H2O");

        // Methane (CH4)
        assert_eq!(format_sum_formula(1, 4, BTreeMap::new(), 0), "CH4");

        // Glycine (C2H5NO2) - alphabetical order N before O
        let mut atom_counts = BTreeMap::new();
        atom_counts.insert(element_symbol_key(e!(N)), (e!(N), 1));
        atom_counts.insert(element_symbol_key(e!(O)), (e!(O), 2));
        assert_eq!(format_sum_formula(2, 5, atom_counts, 0), "C2H5NO2");

        // Charged molecule (NH4+)
        let mut atom_counts = BTreeMap::new();
        atom_counts.insert(element_symbol_key(e!(N)), (e!(N), 1));
        assert_eq!(format_sum_formula(0, 4, atom_counts, 1), "H4N+1");

        // Multiple elements alphabetically sorted
        let mut atom_counts = BTreeMap::new();
        atom_counts.insert(element_symbol_key(e!(Br)), (e!(Br), 2));
        atom_counts.insert(element_symbol_key(e!(Cl)), (e!(Cl), 3));
        atom_counts.insert(element_symbol_key(e!(F)), (e!(F), 1));
        assert_eq!(format_sum_formula(1, 0, atom_counts, 0), "CBr2Cl3F");
    }

    #[test]
    fn test_element_symbol_key() {
        use umol_data::e;

        // Single character elements (padded with 0)
        assert_eq!(element_symbol_key(e!(C)), [b'C', 0]);
        assert_eq!(element_symbol_key(e!(H)), [b'H', 0]);
        assert_eq!(element_symbol_key(e!(N)), [b'N', 0]);
        assert_eq!(element_symbol_key(e!(O)), [b'O', 0]);
        assert_eq!(element_symbol_key(e!(F)), [b'F', 0]);

        // Two character elements
        assert_eq!(element_symbol_key(e!(Br)), [b'B', b'r']);
        assert_eq!(element_symbol_key(e!(Cl)), [b'C', b'l']);

        // Verify alphabetical ordering
        let mut keys = vec![
            element_symbol_key(e!(Cl)),
            element_symbol_key(e!(Br)),
            element_symbol_key(e!(F)),
            element_symbol_key(e!(N)),
            element_symbol_key(e!(O)),
        ];
        keys.sort();
        
        // Should be: Br, Cl, F, N, O (alphabetical)
        assert_eq!(keys, vec![
            [b'B', b'r'], // Br
            [b'C', b'l'], // Cl
            [b'F', 0],    // F
            [b'N', 0],    // N
            [b'O', 0],    // O
        ]);
    }

    #[test]
    fn test_molecule_graph6() {
        use crate::io::ctab::atom::AtomSymbol;
        use crate::io::ctab::bond::BondType;
        use umol_data::{e, Element};

        let mut molecule = Molecule::new();
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(O))));
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));
        molecule.add_atom(Atom::new(AtomSymbol::Element(e!(H))));
        molecule.add_bond(0, 1, Bond::new(BondType::Single));
        molecule.add_bond(0, 2, Bond::new(BondType::Single));
        assert_eq!(molecule.graph6(), "Bo");
    }
}
