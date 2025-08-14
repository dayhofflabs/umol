//! Molecule type for CTab format.

use crate::io::ctab::atom::{Atom, AtomStandard};
use crate::io::ctab::bond::{Bond, BondStandard};
use crate::io::ctab::sgroup::SGroup;
use crate::io::mol::parser::{mol_block, mol_block_standard};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use umol::error::DataError;
use umol::Result;

/// Type aliases for the node and edge indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

/// MOL file header information (3 lines)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Header {
    pub name: String,         // Molecule title / name (line 1)
    pub program_info: String, // Program info (line 2)
    pub comment: String,      // Comment (line 3)
}

impl Header {
    /// Create a new header
    pub fn new(name: String, program_info: String, comment: String) -> Self {
        Self {
            name,
            program_info,
            comment,
        }
    }

    /// Create an empty header
    pub fn empty() -> Self {
        Self {
            name: String::new(),
            program_info: String::new(),
            comment: String::new(),
        }
    }
}

impl Default for Header {
    fn default() -> Self {
        Self::empty()
    }
}

/// Graph-based molecule representation with full MOL file semantics (including queries)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
    pub header: Header,
}

/// Graph-based molecule representation for standard (non-query) molecules only
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MoleculeStandard {
    pub graph: StableGraph<AtomStandard, BondStandard, Undirected, usize>,
    pub sgroups: BTreeMap<usize, SGroup>,
    pub properties: HashMap<String, String>,
    pub header: Header,
}

/// Result of parsing a MOL file with information about query features
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParsedMol {
    molecule: Molecule,
    is_nonstandard_molecule: bool,
}

impl ParsedMol {
    /// Create a new ParsedMol
    pub fn new(molecule: Molecule, is_nonstandard_molecule: bool) -> Self {
        Self {
            molecule,
            is_nonstandard_molecule,
        }
    }

    /// Check if the parsed molecule contains non-standard features
    pub fn has_nonstandard_features(&self) -> bool {
        self.is_nonstandard_molecule
    }

    /// Extract the molecule (works regardless of query features)
    pub fn into_molecule(self) -> Molecule {
        self.molecule
    }

    /// Get a reference to the molecule
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    /// Try to convert to a standard molecule (fails if query features present)
    pub fn try_into_standard(self) -> Result<MoleculeStandard> {
        if self.is_nonstandard_molecule {
            Err(DataError::InvalidFeature(
                "Query features present in molecule marked as standard".to_string(),
            )
            .into())
        } else {
            // TODO: Implement conversion from Molecule to MoleculeStandard
            // For now, this is a placeholder
            Err(DataError::InvalidMolFormat("Conversion not yet implemented".to_string()).into())
        }
    }

    /// Try to convert to a standard molecule, returning a reference to the error
    pub fn as_standard(&self) -> Result<&MoleculeStandard> {
        if self.is_nonstandard_molecule {
            Err(DataError::InvalidFeature(
                "Query features present in molecule marked as standard".to_string(),
            )
            .into())
        } else {
            // TODO: This would need a different approach since we can't return a reference
            // to a converted value that doesn't exist yet
            Err(
                DataError::InvalidMolFormat("Reference conversion not yet implemented".to_string())
                    .into(),
            )
        }
    }
}

/// Parse a MOL string into a ParsedMol
///
/// This is the main parsing function that handles both standard and query molecules.
/// Users can then extract the appropriate type based on their needs.
pub fn parse_mol_str(input: &str) -> Result<ParsedMol> {
    parse_mol(input.as_bytes())
}

/// Parse a MOL string into a MoleculeStandard (optimized, standard molecules only)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_standard_str(input: &str) -> Result<MoleculeStandard> {
    parse_mol_standard(input.as_bytes())
}

/// Parse MOL bytes into a ParsedMol
///
/// This is the primary parsing function. It parses once and allows flexible
/// extraction of either Molecule or MoleculeStandard depending on content.
pub fn parse_mol(input: &[u8]) -> Result<ParsedMol> {
    match mol_block(input) {
        Ok((remaining, parsed_mol)) => {
            // Check if there's unexpected remaining data
            if !remaining.is_empty() && !remaining.iter().all(|&b| b.is_ascii_whitespace()) {
                return Err(DataError::InvalidMolFormat(format!(
                    "Unexpected data after MOL block: {} bytes remaining",
                    remaining.len()
                ))
                .into());
            }
            Ok(parsed_mol)
        }
        Err(nom::Err::Error(e)) | Err(nom::Err::Failure(e)) => {
            Err(DataError::InvalidMolFormat(format!("MOL parsing failed: {:?}", e)).into())
        }
        Err(nom::Err::Incomplete(_)) => {
            Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
        }
    }
}

/// Parse MOL bytes into a MoleculeStandard (optimized, standard molecules only)
///
/// This is the high-performance parsing function for standard molecules.
/// It will fail if the MOL file contains query features.
pub fn parse_mol_standard(input: &[u8]) -> Result<MoleculeStandard> {
    match mol_block_standard(input) {
        Ok((remaining, molecule)) => {
            // Check if there's unexpected remaining data
            if !remaining.is_empty() && !remaining.iter().all(|&b| b.is_ascii_whitespace()) {
                return Err(DataError::InvalidMolFormat(format!(
                    "Unexpected data after MOL block: {} bytes remaining",
                    remaining.len()
                ))
                .into());
            }
            Ok(molecule)
        }
        Err(nom::Err::Error(e)) | Err(nom::Err::Failure(e)) => {
            Err(DataError::InvalidMolFormat(format!("Standard MOL parsing failed: {:?}", e)).into())
        }
        Err(nom::Err::Incomplete(_)) => {
            Err(DataError::InvalidMolFormat("Incomplete MOL data".to_string()).into())
        }
    }
}

impl Molecule {
    /// Create empty molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<Atom, Bond, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
            header: Header::empty(),
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
}

impl MoleculeStandard {
    /// Create empty standard molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<AtomStandard, BondStandard, Undirected, usize>::default(),
            sgroups: BTreeMap::new(),
            properties: HashMap::new(),
            header: Header::empty(),
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

    /// Get immutable reference to atom by index
    pub fn atom(&self, idx: usize) -> Option<&AtomStandard> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut AtomStandard> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get immutable reference to bond by index
    pub fn bond(&self, idx: usize) -> Option<&BondStandard> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut BondStandard> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[test]
    fn test_header_serialize() {
        let header = Header::new(
            "Test Molecule".to_string(),
            "Generated by test".to_string(),
            "Test comment".to_string(),
        );

        // Test YAML serialization
        let yaml = serde_yaml::to_string(&header).expect("Failed to serialize Header to YAML");
        let deserialized: Header =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Header from YAML");
        assert_eq!(header, deserialized);

        // Test JSON serialization
        let json =
            serde_json::to_string_pretty(&header).expect("Failed to serialize Header to JSON");
        let deserialized: Header =
            serde_json::from_str(&json).expect("Failed to deserialize Header from JSON");
        assert_eq!(header, deserialized);
    }

    #[test]
    fn test_molecule_standard_serialize() {
        let header = Header::new("Test".to_string(), "Test".to_string(), "Test".to_string());
        let graph =
            StableGraph::<AtomStandard, BondStandard, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();

        let molecule = MoleculeStandard {
            graph,
            sgroups,
            properties,
            header,
        };

        let yaml =
            serde_yaml::to_string(&molecule).expect("Failed to serialize MoleculeStandard to YAML");
        let deserialized: MoleculeStandard =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize MoleculeStandard from YAML");
        assert_eq!(molecule.header, deserialized.header);
        assert_eq!(molecule.properties, deserialized.properties);
    }

    #[test]
    fn test_molecule_serialize() {
        let header = Header::new("Test".to_string(), "Test".to_string(), "Test".to_string());
        let graph = StableGraph::<Atom, Bond, Undirected, usize>::with_capacity(0, 0);
        let sgroups = BTreeMap::new();
        let properties = HashMap::new();
        let molecule = Molecule {
            graph,
            sgroups,
            properties,
            header,
        };

        let yaml = serde_yaml::to_string(&molecule).expect("Failed to serialize Molecule to YAML");
        let deserialized: Molecule =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Molecule from YAML");
        assert_eq!(molecule.header, deserialized.header);
        assert_eq!(molecule.properties, deserialized.properties);
    }

    #[rstest]
    #[case(b"Methane\nRDKit          3D\nGenerated by RDKit\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0\nM  END\n",
       "valid")]
    fn test_parse_mol(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol(mol_str);
        assert!(
            result.is_ok(),
            "{} should have succeeded: {:?}",
            desc,
            result.err()
        );

        let parsed_mol = result.unwrap();
        assert!(
            !parsed_mol.has_nonstandard_features(),
            "{} should not have query features",
            desc
        );

        let molecule = parsed_mol.into_molecule();
        assert_eq!(molecule.atom_count(), 1, "{} should have 1 atom", desc);
        assert_eq!(molecule.bond_count(), 0, "{} should have 0 bonds", desc);
        assert_eq!(
            molecule.header.name, "Methane",
            "{} header title should match",
            desc
        );
        assert_eq!(
            molecule.header.program_info, "RDKit          3D",
            "{} header program info should match",
            desc
        );
        assert_eq!(
            molecule.header.comment, "Generated by RDKit",
            "{} header comment should match",
            desc
        );
    }

    #[rstest]
    #[case(b"Ethane\nRDKit          3D\nGenerated by RDKit\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n",
      "standard ethane")]
    fn test_parse_mol_standard(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol_standard(mol_str);
        assert!(
            result.is_ok(),
            "{} should have succeeded: {:?}",
            desc,
            result.err()
        );

        let molecule = result.unwrap();
        assert_eq!(molecule.atom_count(), 2, "{} should have 2 atoms", desc);
        assert_eq!(molecule.bond_count(), 1, "{} should have 1 bond", desc);
        assert_eq!(
            molecule.header.name, "Ethane",
            "{} header title should match",
            desc
        );
        assert_eq!(
            molecule.header.program_info, "RDKit          3D",
            "{} header program info should match",
            desc
        );
        assert_eq!(
            molecule.header.comment, "Generated by RDKit",
            "{} header comment should match",
            desc
        );
    }

    #[test]
    fn test_parse_mol_standard_query() {
        // MOL with query atom 'A' (any atom)
        let query_mol = b"Query\nRDKit          3D\nWith query atom\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

        let result = parse_mol_standard(query_mol);
        assert!(
            result.is_err(),
            "Standard parser should fail on query atoms"
        );

        let error = result.unwrap_err();
        let error_string = format!("{}", error);
        assert!(
            error_string.contains("MOL parsing failed")
                || error_string.contains("Standard MOL parsing failed"),
            "Error should mention parsing failure: {}",
            error_string
        );
    }
}
