//! Common utilities and test infrastructure for compliance testing

use std::path::{Path, PathBuf};
use umol_models_graph::io::ctab::molecule::{parse_mol_str, Molecule};

/// Load a MOL test file from compliance test data
/// filename can be relative to tests/compliance/data/mol/ (e.g., "basic/valid/seaborgium.mol")
pub fn load_test_mol_file(filename: &str) -> Result<String, std::io::Error> {
    let test_data_path = Path::new("tests/compliance/data/mol").join(filename);
    std::fs::read_to_string(&test_data_path)
}

/// Discover MOL files in a specific category directory
pub fn discover_mol_files(category: &str) -> Result<Vec<PathBuf>, std::io::Error> {
    let dir = Path::new("tests/compliance/data/mol").join(category);
    if !dir.exists() {
        return Ok(Vec::new());
    }
    
    let mut mol_files = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) == Some("mol") {
            mol_files.push(path);
        }
    }
    
    Ok(mol_files)
}

/// Parse MOL file with error handling
pub fn parse_mol_file(content: &str) -> Result<Molecule, String> {
    parse_mol_str(content)
        .map(|parsed| parsed.molecule().clone())
        .map_err(|e| format!("Parse error: {}", e))
}

/// Parse MOL file with lenient handling of trailing data after M END
pub fn try_parse_mol_file(content: &str) -> Result<Molecule, String> {
    match parse_mol_str(content) {
        Ok(parsed) => Ok(parsed.molecule().clone()),
        Err(e) => {
            let error_msg = format!("{}", e);
            if error_msg.contains("Unexpected data after MOL block") {
                if let Some(end_pos) = content.find("M  END") {
                    let truncated = &content[..end_pos + 6]; // Include "M  END"
                    match parse_mol_str(truncated) {
                        Ok(parsed) => return Ok(parsed.molecule().clone()),
                        Err(_) => return Err(format!("Parse error: {}", e)),
                    }
                }
            }
            Err(format!("Parse error: {}", e))
        }
    }
}

/// Basic assertion that a molecule has valid structure
pub fn assert_basic_molecule_structure(molecule: &Molecule) {
    assert!(molecule.graph.node_count() > 0, "Molecule should have atoms");
    
    // Check that all bonds reference valid atoms
    for edge_index in molecule.graph.edge_indices() {
        let endpoints = molecule.graph.edge_endpoints(edge_index);
        assert!(endpoints.is_some(), "Bond should have valid endpoints");
    }
}