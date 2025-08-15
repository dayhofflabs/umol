//! Testing utils for MOL files

use std::fs::{read_dir, read_to_string};
use std::io;
use std::path::{Path, PathBuf};

/// Load a MOL test file from test data
/// filename can be relative to tests/data/mol/ (e.g., "basic/valid/seaborgium.mol")
pub fn load_test_mol_file(filename: &str) -> Result<String, io::Error> {
    let test_data_path = Path::new("tests/data/mol").join(filename);
    read_to_string(&test_data_path)
}

/// Discover MOL files in a specific category directory
pub fn discover_mol_files(category: &str) -> Result<Vec<PathBuf>, io::Error> {
    let dir = Path::new("tests/data/mol").join(category);
    if !dir.exists() {
        return Ok(Vec::new());
    }

    let mut mol_files = Vec::new();
    for entry in read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) == Some("mol") {
            mol_files.push(path);
        }
    }
    Ok(mol_files)
}