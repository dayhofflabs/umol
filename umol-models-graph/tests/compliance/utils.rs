//! Common utilities and test infrastructure for compliance testing

use std::path::{Path, PathBuf};
use serde::{Deserialize, Serialize};

/// Test file categories for automatic discovery and organization
#[derive(Debug, Clone)]
pub struct TestFileCategories {
    pub basic_molecules: Vec<PathBuf>,
    pub property_examples: Vec<PathBuf>,
    pub sgroup_examples: Vec<PathBuf>,
    pub rgroup_examples: Vec<PathBuf>,
    pub query_features: Vec<PathBuf>,
    pub real_world_complex: Vec<PathBuf>,
    pub edge_cases: Vec<PathBuf>,
    pub known_failures: Vec<PathBuf>,
}

/// Complexity level of test cases
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TestComplexity {
    Simple,
    Properties,
    SGroups,
    RGroups,
    Query,
    Complex,
    EdgeCase,
}

/// Expected test outcome
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TestExpectation {
    ShouldParse,
    ShouldError,
    ShouldWarn,
    KnownFailure,
}

/// Load a MOL test file from the compliance test data directory
pub fn load_test_mol_file(filename: &str) -> Result<String, std::io::Error> {
    let test_data_path = Path::new("tests/compliance/data/mol").join(filename);
    std::fs::read_to_string(&test_data_path)
}

/// Discover all MOL files in a directory
pub fn discover_mol_files(directory: &Path) -> Result<Vec<PathBuf>, std::io::Error> {
    let mut mol_files = Vec::new();
    
    if directory.exists() {
        for entry in std::fs::read_dir(directory)? {
            let path = entry?.path();
            if path.extension().map_or(false, |ext| ext == "mol") {
                mol_files.push(path);
            }
        }
    }
    
    mol_files.sort();
    Ok(mol_files)
}

/// Categorize test files based on filename patterns
pub fn categorize_test_files() -> Result<TestFileCategories, std::io::Error> {
    let test_data_dir = Path::new("tests/compliance/data/mol");
    let files = discover_mol_files(test_data_dir)?;
    
    let mut categories = TestFileCategories {
        basic_molecules: Vec::new(),
        property_examples: Vec::new(),
        sgroup_examples: Vec::new(),
        rgroup_examples: Vec::new(),
        query_features: Vec::new(),
        real_world_complex: Vec::new(),
        edge_cases: Vec::new(),
        known_failures: Vec::new(),
    };
    
    for file in files {
        let filename = file.file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        
        match categorize_by_filename(filename) {
            (TestComplexity::Simple, _) => categories.basic_molecules.push(file),
            (TestComplexity::Properties, _) => categories.property_examples.push(file),
            (TestComplexity::SGroups, _) => categories.sgroup_examples.push(file),
            (TestComplexity::RGroups, _) => categories.rgroup_examples.push(file),
            (TestComplexity::Query, _) => categories.query_features.push(file),
            (TestComplexity::Complex, _) => categories.real_world_complex.push(file),
            (TestComplexity::EdgeCase, _) => categories.edge_cases.push(file),
        }
    }
    
    Ok(categories)
}

/// Categorize a file based on its filename
pub fn categorize_by_filename(filename: &str) -> (TestComplexity, TestExpectation) {
    // Convert to lowercase for case-insensitive matching
    let name = filename.to_lowercase();
    
    // Known simple molecules
    if matches!(name.as_str(), "h3bnh3.mol" | "butanoic_acid.mol" | "seaborgium.mol" | "seaborgium_abs.mol") {
        return (TestComplexity::Simple, TestExpectation::ShouldParse);
    }
    
    // Property examples
    if name.contains("atom") && (name.contains("value") || name.contains("alias")) 
        || name.contains("isotope") || name.contains("radical") 
        || name.contains("charge") || name.contains("hyd") {
        return (TestComplexity::Properties, TestExpectation::ShouldParse);
    }
    
    // S-Group examples
    if name.starts_with("sgroup") || name.contains("sgroups") {
        return (TestComplexity::SGroups, TestExpectation::ShouldParse);
    }
    
    // R-Group examples  
    if name.starts_with("rg") || name.contains("rgroup") {
        return (TestComplexity::RGroups, TestExpectation::ShouldParse);
    }
    
    // Query features
    if name.starts_with("query") || name.contains("atomlist") || name.contains("link") {
        return (TestComplexity::Query, TestExpectation::ShouldParse);
    }
    
    // Edge cases
    if name.contains("short") || name.contains("bug") || name.contains("bad") {
        return (TestComplexity::EdgeCase, TestExpectation::ShouldParse);
    }
    
    // Complex real-world molecules (ChEBI, large structures)
    if name.contains("chebi") || name.contains("coordination") 
        || name.contains("peptide") || name.contains("phosphate") {
        return (TestComplexity::Complex, TestExpectation::ShouldParse);
    }
    
    // Default to simple for now
    (TestComplexity::Simple, TestExpectation::ShouldParse)
}

/// Molecule snapshot for serialization testing
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct MoleculeSnapshot {
    pub format_version: String,
    pub source_file: String,
    pub metadata: MoleculeMetadata,
    pub atoms: Vec<AtomSnapshot>,
    pub bonds: Vec<BondSnapshot>,
    pub properties: PropertySnapshot,
    pub validation: ValidationSnapshot,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct MoleculeMetadata {
    pub atom_count: usize,
    pub bond_count: usize,
    pub has_properties: bool,
    pub has_sgroups: bool,
    pub has_rgroups: bool,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct AtomSnapshot {
    pub index: usize,
    pub element: String,
    pub coordinates: [f64; 3],
    pub charge: i8,
    pub hydrogen_count: Option<u8>,
    pub properties: std::collections::HashMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct BondSnapshot {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_type: String,
    pub stereo: String,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct PropertySnapshot {
    pub charges: Vec<(usize, i8)>,
    pub isotopes: Vec<(usize, u32)>,
    pub radicals: Vec<(usize, String)>,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ValidationSnapshot {
    pub bond_indices_valid: bool,
    pub atom_properties_consistent: bool,
    pub molecular_formula: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_categorize_by_filename() {
        assert_eq!(
            categorize_by_filename("H3BNH3.mol"),
            (TestComplexity::Simple, TestExpectation::ShouldParse)
        );
        
        assert_eq!(
            categorize_by_filename("atomValueLines.mol"),
            (TestComplexity::Properties, TestExpectation::ShouldParse)
        );
        
        assert_eq!(
            categorize_by_filename("sgroup-peptide.mol"),
            (TestComplexity::SGroups, TestExpectation::ShouldParse)
        );
        
        assert_eq!(
            categorize_by_filename("rgfile.1.mol"),
            (TestComplexity::RGroups, TestExpectation::ShouldParse)
        );
    }
}
