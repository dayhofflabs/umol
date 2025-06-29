//! Molecular graph model.

use crate::atom::{Atom, AtomStandard, AtomSymbol};
use crate::bond::{Bond, BondType};
use crate::conformer::{Conformer, Point3D};
use crate::io::ctab::mol::{mol_block, mol_block_standard};
use crate::sgroup::SGroup;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableGraph;
use petgraph::Undirected;
use std::collections::HashMap;
use umol::error::DataError;
use umol::Result;

/// Type aliases for the node and edge indices
pub type AtomIndex = NodeIndex<usize>;
pub type BondIndex = EdgeIndex<usize>;

/// MOL file header information (3 lines)
#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    /// Molecule title/name (line 1)
    pub title: String,
    /// Program info (line 2)
    pub program_info: String,
    /// Comment (line 3)
    pub comment: String,
}

impl Header {
    /// Create a new header
    pub fn new(title: String, program_info: String, comment: String) -> Self {
        Self {
            title,
            program_info,
            comment,
        }
    }

    /// Create an empty header
    pub fn empty() -> Self {
        Self {
            title: String::new(),
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
#[derive(Debug, Clone)]
pub struct Molecule {
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    pub conformers: Vec<Conformer>,
    pub sgroups: Vec<SGroup>,
    pub properties: HashMap<String, String>,
    pub header: Header,
}

/// Graph-based molecule representation for standard (non-query) molecules only
#[derive(Debug, Clone)]
pub struct MoleculeStandard {
    pub graph: StableGraph<AtomStandard, Bond, Undirected, usize>,
    pub conformers: Vec<Conformer>,
    pub sgroups: Vec<SGroup>,
    pub properties: HashMap<String, String>,
    pub header: Header,
}

/// Result of parsing a MOL file with information about query features
#[derive(Debug, Clone)]
pub struct ParsedMol {
    molecule: Molecule,
    is_query_molecule: bool,
}

impl ParsedMol {
    /// Create a new ParsedMol
    pub fn new(molecule: Molecule, is_query_molecule: bool) -> Self {
        Self {
            molecule,
            is_query_molecule,
        }
    }

    /// Check if the parsed molecule contains query features
    pub fn has_query_features(&self) -> bool {
        self.is_query_molecule
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
        if self.is_query_molecule {
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
    pub fn as_standard(&self) -> std::result::Result<&MoleculeStandard, umol::Error> {
        if self.is_query_molecule {
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
            conformers: Vec::new(),
            sgroups: Vec::new(),
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
    ///
    /// - `atom`: Atom to add (Molecule takes ownership)
    ///
    /// Return index of added atom.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Add bond between two atoms specified by external/MOL indices
    ///
    /// - `idx1`, `idx2`: Atom indices
    /// - `bond`: Bond to add (Molecule takes ownership)
    ///
    /// Return index of added bond.
    pub fn add_bond(&mut self, idx1: usize, idx2: usize, bond: Bond) -> usize {
        self.graph
            .add_edge(AtomIndex::new(idx1), AtomIndex::new(idx2), bond)
            .index()
    }

    /// Get immutable reference to atom by index
    ///
    /// - `idx`: Atom index
    ///
    /// Return immutable reference to atom.
    pub fn atom(&self, idx: usize) -> Option<&Atom> {
        self.graph.node_weight(AtomIndex::new(idx))
    }

    /// Get mutable reference to atom by index
    ///
    /// - `idx`: Atom index
    ///
    /// Return mutable reference to atom.
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut Atom> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }

    /// Get immutable reference to bond by index
    ///
    /// - `idx`: Bond index
    ///
    /// Return immutable reference to bond.
    pub fn bond(&self, idx: usize) -> Option<&Bond> {
        self.graph.edge_weight(BondIndex::new(idx))
    }

    /// Get mutable reference to bond by index
    ///
    /// - `idx`: Bond index
    ///
    /// Return mutable reference to bond.
    pub fn bond_mut(&mut self, idx: usize) -> Option<&mut Bond> {
        self.graph.edge_weight_mut(BondIndex::new(idx))
    }

    /// Get iterator over neighbor atom indices for atom index
    pub fn neighbors(&self, idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.graph.neighbors(AtomIndex::new(idx)).map(|i| i.index())
    }

    /// Get immutable slice of all conformers
    pub fn conformers(&self) -> &[Conformer] {
        &self.conformers
    }

    /// Get mutable reference to vector of conformers
    pub fn conformers_mut(&mut self) -> &mut Vec<Conformer> {
        &mut self.conformers
    }

    /// Add conformer to the molecule
    ///
    /// - `conformer`: Conformer to add
    ///
    /// Return error if the number of positions in the conformer does not match
    /// the number of atoms in the molecule.
    pub fn add_conformer(&mut self, conformer: Conformer) -> Result<()> {
        let num_atoms = self.atom_count();
        let num_positions = conformer.positions.len();

        if num_atoms != num_positions {
            return Err(DataError::InvalidConformer(format!(
                "Expected {} positions, found {}",
                num_atoms, num_positions
            ))
            .into());
        }

        self.conformers.push(conformer);
        Ok(())
    }
}

impl MoleculeStandard {
    /// Create empty standard molecule
    pub fn new() -> Self {
        Self {
            graph: StableGraph::<AtomStandard, Bond, Undirected, usize>::default(),
            conformers: Vec::new(),
            sgroups: Vec::new(),
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
    ///
    /// - `atom`: Atom to add (MoleculeStandard takes ownership)
    ///
    /// Return index of added atom.
    pub fn add_atom(&mut self, atom: AtomStandard) -> usize {
        self.graph.add_node(atom).index()
    }

    /// Get mutable reference to atom by index
    ///
    /// - `idx`: Atom index
    ///
    /// Return mutable reference to atom.
    pub fn atom_mut(&mut self, idx: usize) -> Option<&mut AtomStandard> {
        self.graph.node_weight_mut(AtomIndex::new(idx))
    }
}

impl Molecule {
    /// Serialize to MOL format string
    pub fn to_mol_string(&self) -> String {
        let mut output = String::new();
        self.write_mol_to_string(&mut output);
        output
    }

    /// Write MOL format to writer
    pub fn write_mol<W: std::io::Write>(&self, mut writer: W) -> Result<()> {
        let mol_string = self.to_mol_string();
        writer.write_all(mol_string.as_bytes()).map_err(|e| {
            let serialization_error: umol::error::SerializationError =
                umol::error::SerializationError::IoError(e);
            let umol_error: umol::Error = serialization_error.into();
            umol_error
        })?;
        Ok(())
    }

    /// Internal method to write MOL format to string
    fn write_mol_to_string(&self, output: &mut String) {
        // Write header (3 lines)
        output.push_str(&self.header.title);
        output.push('\n');
        output.push_str(&self.header.program_info);
        output.push('\n');
        output.push_str(&self.header.comment);
        output.push('\n');

        // Write counts line
        let atom_count = self.graph.node_count();
        let bond_count = self.graph.edge_count();
        output.push_str(&format!(
            "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
            atom_count, bond_count
        ));

        // Get first conformer for coordinates (or use zero coordinates)
        let empty_coords: Vec<Point3D> = Vec::new();
        let coordinates = self
            .conformers
            .first()
            .map(|c| &c.positions)
            .unwrap_or_else(|| &empty_coords);

        // Write atom block
        for (i, node_idx) in self.graph.node_indices().enumerate() {
            if let Some(atom) = self.graph.node_weight(node_idx) {
                let coord = coordinates
                    .get(i)
                    .copied()
                    .unwrap_or(Point3D::new(0.0, 0.0, 0.0));

                // Format: x10.4, y10.4, z10.4, symbol3, mass_diff2, charge3
                let symbol_str = match &atom.symbol {
                    AtomSymbol::Element(element) => element.symbol().to_string(),
                    AtomSymbol::NamedIsotope(isotope) => isotope.element().symbol().to_string(),
                    AtomSymbol::AtomList(atom_list) => {
                        atom_list.elements.first().unwrap_or(&umol_data::Element::C).symbol().to_string()
                    },
                    AtomSymbol::Unspecified(c) => c.to_string(),
                    AtomSymbol::LonePair => "LP".to_string(),
                    AtomSymbol::RGroup(n) => format!("R{}", n),
                };

                // Use precise F10.4 format: 10 characters wide, 4 decimal places, right-aligned
                // Symbol is exactly 3 characters, left-aligned after a single space
                output.push_str(&format!(
                    "{:10.4}{:10.4}{:10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                    coord.x, coord.y, coord.z, symbol_str
                ));
            }
        }

        // Write bond block
        for edge_idx in self.graph.edge_indices() {
            if let Some((idx1, idx2)) = self.graph.edge_endpoints(edge_idx) {
                if let Some(bond) = self.graph.edge_weight(edge_idx) {
                    let idx1_1based = idx1.index() + 1;
                    let idx2_1based = idx2.index() + 1;
                    let bond_type_code = match bond.bond_type {
                        BondType::Single => 1,
                        BondType::Double => 2,
                        BondType::Triple => 3,
                        BondType::Aromatic => 4,
                        _ => 1, // Fallback
                    };

                    output.push_str(&format!(
                        "{:3}{:3}{:3}  0  0  0  0\n",
                        idx1_1based, idx2_1based, bond_type_code
                    ));
                }
            }
        }

        // Write properties block (simplified for now)
        // TODO: Implement property serialization

        // Write M  END
        output.push_str("M  END\n");
    }
}

impl MoleculeStandard {
    /// Serialize to MOL format string
    pub fn to_mol_string(&self) -> String {
        let mut output = String::new();
        self.write_mol_to_string(&mut output);
        output
    }

    /// Write MOL format to writer
    pub fn write_mol<W: std::io::Write>(&self, mut writer: W) -> Result<()> {
        let mol_string = self.to_mol_string();
        writer.write_all(mol_string.as_bytes()).map_err(|e| {
            let serialization_error: umol::error::SerializationError =
                umol::error::SerializationError::IoError(e);
            let umol_error: umol::Error = serialization_error.into();
            umol_error
        })?;
        Ok(())
    }

    /// Internal method to write MOL format to string
    fn write_mol_to_string(&self, output: &mut String) {
        // Write header (3 lines)
        output.push_str(&self.header.title);
        output.push('\n');
        output.push_str(&self.header.program_info);
        output.push('\n');
        output.push_str(&self.header.comment);
        output.push('\n');

        // Write counts line
        let atom_count = self.graph.node_count();
        let bond_count = self.graph.edge_count();
        output.push_str(&format!(
            "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
            atom_count, bond_count
        ));

        // Get first conformer for coordinates (or use zero coordinates)
        let empty_coords: Vec<Point3D> = Vec::new();
        let coordinates = self
            .conformers
            .first()
            .map(|c| &c.positions)
            .unwrap_or(&empty_coords);

        // Write atom block
        for (i, node_idx) in self.graph.node_indices().enumerate() {
            if let Some(atom) = self.graph.node_weight(node_idx) {
                let coord = coordinates
                    .get(i)
                    .copied()
                    .unwrap_or(Point3D::new(0.0, 0.0, 0.0));

                // Format: x10.4, y10.4, z10.4, symbol3, mass_diff2, charge3
                let symbol_str = atom.element.symbol().to_string();

                // Use precise F10.4 format: 10 characters wide, 4 decimal places, right-aligned
                // Symbol is exactly 3 characters, left-aligned after a single space
                output.push_str(&format!(
                    "{:10.4}{:10.4}{:10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                    coord.x, coord.y, coord.z, symbol_str
                ));
            }
        }

        // Write bond block
        for edge_idx in self.graph.edge_indices() {
            if let Some((idx1, idx2)) = self.graph.edge_endpoints(edge_idx) {
                if let Some(bond) = self.graph.edge_weight(edge_idx) {
                    let idx1_1based = idx1.index() + 1;
                    let idx2_1based = idx2.index() + 1;
                    let bond_type_code = match bond.bond_type {
                        BondType::Single => 1,
                        BondType::Double => 2,
                        BondType::Triple => 3,
                        BondType::Aromatic => 4,
                        _ => 1, // Fallback
                    };

                    output.push_str(&format!(
                        "{:3}{:3}{:3}  0  0  0  0\n",
                        idx1_1based, idx2_1based, bond_type_code
                    ));
                }
            }
        }

        // Write properties block (simplified for now)
        // TODO: Implement property serialization

        // Write M  END
        output.push_str("M  END\n");
    }
}

// Standard library integration
impl std::fmt::Display for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_mol_string())
    }
}

impl std::fmt::Display for MoleculeStandard {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_mol_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

         #[rstest]
     #[case(b"Methane\nRDKit          3D\nGenerated by RDKit\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C  -2  3  0  0  0  4  0  0  0  0  0  0\nM  END\n",
       "valid")]
     fn test_parse_mol(#[case] mol_str: &[u8], #[case] desc: &str) {
         let result = parse_mol(mol_str);
         assert!(result.is_ok(), "{} should have succeeded: {:?}", desc, result.err());

        let parsed_mol = result.unwrap();
        assert!(
            !parsed_mol.has_query_features(),
            "{} should not have query features",
            desc
        );

        let molecule = parsed_mol.into_molecule();
        assert_eq!(molecule.atom_count(), 1, "{} should have 1 atom", desc);
        assert_eq!(molecule.bond_count(), 0, "{} should have 0 bonds", desc);
        assert_eq!(
            molecule.header.title, "Methane",
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
    fn test_write_mol() {
        // Create a simple molecule
        let mut molecule = Molecule::new();
        molecule.header = Header::new(
            "Test Molecule".to_string(),
            "umol-test      ".to_string(),
            "Test comment".to_string(),
        );

        // Add a carbon atom
        let atom = Atom::new(AtomSymbol::Element(umol_data::Element::C));
        molecule.add_atom(atom);

        // Test serialization
        let mol_string = molecule.to_mol_string();

        // Basic checks
        assert!(mol_string.contains("Test Molecule"));
        assert!(mol_string.contains("umol-test"));
        assert!(mol_string.contains("Test comment"));
        assert!(mol_string.contains("  1  0  0  0  0  0  0  0  0  0999 V2000"));
        assert!(mol_string.contains("M  END"));

        println!("Generated MOL:\n{}", mol_string);
    }

    #[rstest]
    #[case(b"Ethane\nRDKit          3D\nGenerated by RDKit\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n",
      "standard ethane")]
    fn test_parse_mol_standard(#[case] mol_str: &[u8], #[case] desc: &str) {
        let result = parse_mol_standard(mol_str);
        assert!(result.is_ok(), "{} should have succeeded: {:?}", desc, result.err());

        let molecule = result.unwrap();
        assert_eq!(molecule.atom_count(), 2, "{} should have 2 atoms", desc);
        assert_eq!(molecule.bond_count(), 1, "{} should have 1 bond", desc);
        assert_eq!(
            molecule.header.title, "Ethane",
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
        assert!(result.is_err(), "Standard parser should fail on query atoms");
        
        let error = result.unwrap_err();
        let error_string = format!("{}", error);
        assert!(error_string.contains("MOL parsing failed") || error_string.contains("Standard MOL parsing failed"), 
                "Error should mention parsing failure: {}", error_string);
    }

    #[test]
    fn test_write_mol_standard() {
        // Create a simple standard molecule
        let mut molecule = MoleculeStandard::new();
        molecule.header = Header::new(
            "Test Standard".to_string(),
            "umol-standard  ".to_string(),
            "Standard comment".to_string(),
        );

        // Add two carbon atoms
        let atom1 = AtomStandard::new(umol_data::Element::C);
        let atom2 = AtomStandard::new(umol_data::Element::N);
        molecule.add_atom(atom1);
        molecule.add_atom(atom2);

        // Add a bond directly to the graph (since MoleculeStandard doesn't have add_bond)
        let bond = Bond::new(BondType::Single);
        molecule.graph.add_edge(AtomIndex::new(0), AtomIndex::new(1), bond);

        // Test serialization
        let mol_string = molecule.to_mol_string();

        // Basic checks
        assert!(mol_string.contains("Test Standard"));
        assert!(mol_string.contains("umol-standard"));
        assert!(mol_string.contains("Standard comment"));
        assert!(mol_string.contains("  2  1  0  0  0  0  0  0  0  0999 V2000"));
        assert!(mol_string.contains("C  "));
        assert!(mol_string.contains("N  "));
        assert!(mol_string.contains("  1  2  1  0  0  0  0"));
        assert!(mol_string.contains("M  END"));

        println!("Generated Standard MOL:\n{}", mol_string);
    }

    #[test]
    fn test_round_trip_standard() {
        // Test that we can parse a standard MOL and write it back
        let original_mol = b"Methanol\nRDKit          3D\nSimple molecule\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n";
        
        // Parse with standard parser
        let parsed = parse_mol_standard(original_mol).unwrap();
        
        // Write back to string
        let regenerated = parsed.to_mol_string();
        println!("Original bytes: {:?}", original_mol);
        println!("Regenerated string: {:?}", regenerated);
        println!("Regenerated bytes: {:?}", regenerated.as_bytes());
        
        // Parse the regenerated string to verify it's valid
        let reparsed = parse_mol_standard(regenerated.as_bytes()).unwrap();
        
        // Verify structure is preserved
        assert_eq!(reparsed.atom_count(), 2);
        assert_eq!(reparsed.bond_count(), 1);
        assert_eq!(reparsed.header.title, "Methanol");
    }
}
