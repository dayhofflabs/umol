//! Real-world MOL file parsing tests using categorized test data

use crate::common::*;
use insta::assert_yaml_snapshot;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_molecules_exist() {
        let mol_files = discover_mol_files("basic/valid")
            .expect("Should be able to read basic test data directory");
        assert!(!mol_files.is_empty(), "Should have basic MOL test files");
        println!("Found {} basic MOL files", mol_files.len());
    }
    
    #[test]
    fn test_parse_basic_molecules() {
        let mol_files = match discover_mol_files("basic/valid") {
            Ok(files) => files,
            Err(_) => {
                println!("Skipping test - basic test files not available");
                return;
            }
        };
        
        let mut successful_parses = 0;
        let mut failed_parses = 0;
        
        for mol_file in &mol_files {
            let filename = mol_file.file_name().unwrap().to_str().unwrap();
            
            match std::fs::read_to_string(&mol_file) {
                Ok(content) => {
                    match try_parse_mol_file(&content) {
                        Ok(molecule) => {
                            successful_parses += 1;
                            
                            println!("✓ Parsed basic molecule: {} ({} atoms, {} bonds)", 
                                     filename, molecule.graph.node_count(), molecule.graph.edge_count());
                        }
                        Err(e) => {
                            failed_parses += 1;
                            println!("✗ Failed to parse {}: {}", filename, e);
                        }
                    }
                }
                Err(e) => {
                    failed_parses += 1;
                    println!("✗ Failed to load {}: {}", filename, e);
                }
            }
        }
        
        println!("Basic molecules: {} successful, {} failed", successful_parses, failed_parses);
        
        if successful_parses + failed_parses > 0 {
            assert!(successful_parses > 0, "Should parse at least some basic molecules successfully");
            let success_rate = successful_parses as f64 / (successful_parses + failed_parses) as f64 * 100.0;
            println!("Success rate: {:.1}%", success_rate);
        }
    }

    #[test] 
    fn test_snapshot_seaborgium() {
        let content = match load_test_mol_file("basic/valid/seaborgium.mol") {
            Ok(content) => content,
            Err(_) => {
                println!("Skipping snapshot test - seaborgium.mol not found");
                return;
            }
        };
        
        let molecule = try_parse_mol_file(&content)
            .expect("Should be able to parse seaborgium.mol");
        
        // Create a snapshot test using direct serde serialization
        assert_yaml_snapshot!("seaborgium_molecule", molecule);
    }

    #[test]
    fn test_snapshot_element_with_alias() {
        let content = match load_test_mol_file("with_properties/valid/element-with-alias.mol") {
            Ok(content) => content,
            Err(_) => {
                println!("Skipping snapshot test - element-with-alias.mol not found");
                return;
            }
        };
        
        let molecule = try_parse_mol_file(&content)
            .expect("Should be able to parse element-with-alias.mol");
        
        // Create a snapshot test using direct serde serialization
        assert_yaml_snapshot!("element_with_alias_molecule", molecule);
    }
    
    #[test]
    fn test_parse_properties_molecules() {
        let mol_files = match discover_mol_files("with_properties/valid") {
            Ok(files) => files,
            Err(_) => {
                println!("Skipping test - properties test files not available");
                return;
            }
        };
        
        let mut successful_parses = 0;
        let mut failed_parses = 0;
        
        for mol_file in &mol_files {
            let filename = mol_file.file_name().unwrap().to_str().unwrap();
            
            match std::fs::read_to_string(&mol_file) {
                Ok(content) => {
                    match try_parse_mol_file(&content) {
                        Ok(_molecule) => {
                            successful_parses += 1;
                            println!("✓ Parsed properties molecule: {}", filename);
                        }
                        Err(e) => {
                            failed_parses += 1;
                            println!("✗ Failed to parse {}: {}", filename, e);
                        }
                    }
                }
                Err(e) => {
                    failed_parses += 1;
                    println!("✗ Failed to load {}: {}", filename, e);
                }
            }
        }
        
        println!("Properties molecules: {} successful, {} failed", successful_parses, failed_parses);
        
        if successful_parses + failed_parses > 0 {
            assert!(successful_parses > 0, "Should parse at least some properties molecules successfully");
        }
    }
}