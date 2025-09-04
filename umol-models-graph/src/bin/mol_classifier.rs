use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::os::unix::fs as unix_fs;
use umol_models_graph::io::mol::parser::{parse_mol, parse_mol_moleculelike};

#[derive(Debug)]
struct ClassificationStats {
    molecule: usize,
    moleculelike: usize,
    invalid: usize,
    total: usize,
}

impl ClassificationStats {
    fn new() -> Self {
        Self {
            molecule: 0,
            moleculelike: 0,
            invalid: 0,
            total: 0,
        }
    }

    fn add_molecule(&mut self) {
        self.molecule += 1;
        self.total += 1;
    }

    fn add_moleculelike(&mut self) {
        self.moleculelike += 1;
        self.total += 1;
    }

    fn add_invalid(&mut self) {
        self.invalid += 1;
        self.total += 1;
    }

    #[allow(dead_code)]
    fn molecule_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            (self.molecule as f64 / self.total as f64) * 100.0
        }
    }

    fn valid_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            ((self.molecule + self.moleculelike) as f64 / self.total as f64) * 100.0
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum FileClassification {
    Molecule,      // Works with parse_mol (basic parser)
    MoleculeLike,  // Works with parse_mol_moleculelike but not parse_mol  
    Invalid,       // Doesn't work with either
}

fn classify_mol_file(file_path: &Path) -> Result<FileClassification, Box<dyn std::error::Error>> {
    let mol_bytes = fs::read(file_path)?;
    
    // Reject files that don't look like MOL files (e.g. RXN files)
    if mol_bytes.starts_with(b"$RXN") {
        return Ok(FileClassification::Invalid);
    }
    
    // Try basic parser first
    match parse_mol(&mol_bytes) {
        Ok(_mol) => Ok(FileClassification::Molecule),
        Err(_) => {
            // Try extended parser if basic fails
            match parse_mol_moleculelike(&mol_bytes) {
                Ok(_mol) => Ok(FileClassification::MoleculeLike),
                Err(_) => Ok(FileClassification::Invalid),
            }
        }
    }
}

fn clean_existing_symlinks(data_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let categories = ["molecule", "moleculelike", "invalid"];
    
    for category in categories {
        let category_path = Path::new(data_path).join(category);
        if !category_path.exists() {
            continue;
        }
        
        // Remove all subdirectories (which contain the symlinks)
        for entry in fs::read_dir(&category_path)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                println!("Removing symlinks in {}", path.display());
                if let Err(e) = fs::remove_dir_all(&path) {
                    eprintln!("Warning: Failed to remove {}: {}", path.display(), e);
                }
            }
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let should_sort = args.len() > 1 && args[1] == "--sort";
    
    let data_raw_path = "tests/mol_parsing/data_raw";
    let data_path = "tests/mol_parsing/data";
    
    if !Path::new(data_raw_path).exists() {
        eprintln!("Error: {} directory not found", data_raw_path);
        std::process::exit(1);
    }

    // Clean existing symlinks before sorting
    if should_sort {
        println!("Cleaning existing symlinks...");
        clean_existing_symlinks(data_path)?;
    }

    let mut source_stats: HashMap<String, ClassificationStats> = HashMap::new();
    let mut total_files = 0;
    let mut processed_files = 0;
    let mut error_files = 0;

    // Walk through all source directories
    for entry in fs::read_dir(data_raw_path)? {
        let entry = entry?;
        let source_path = entry.path();
        
        if !source_path.is_dir() {
            continue;
        }

        let source_name = source_path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();

        let mut stats = ClassificationStats::new();

        // Find all .mol files in this source directory (recursive)
        fn collect_mol_files(dir: &Path) -> Result<Vec<std::path::PathBuf>, Box<dyn std::error::Error>> {
            let mut files = Vec::new();
            for entry in fs::read_dir(dir)? {
                let entry = entry?;
                let path = entry.path();
                
                if path.is_dir() {
                    files.extend(collect_mol_files(&path)?);
                } else if path.extension() == Some(std::ffi::OsStr::new("mol")) {
                    files.push(path);
                }
            }
            Ok(files)
        }

        let mol_files = collect_mol_files(&source_path)?;
        
        for path in mol_files {
            total_files += 1;
            
            match classify_mol_file(&path) {
                Ok(classification) => {
                    match classification {
                        FileClassification::Molecule => stats.add_molecule(),
                        FileClassification::MoleculeLike => stats.add_moleculelike(),
                        FileClassification::Invalid => stats.add_invalid(),
                    }
                    
                    // Sort files if requested
                    if should_sort {
                        let category = match classification {
                            FileClassification::Molecule => "molecule",
                            FileClassification::MoleculeLike => "moleculelike",
                            FileClassification::Invalid => "invalid",
                        };
                        
                        let dest_dir = Path::new(data_path).join(category).join(&source_name);
                        let dest_file = dest_dir.join(path.file_name().unwrap());
                        
                        // Create destination directory if it doesn't exist
                        if let Err(e) = fs::create_dir_all(&dest_dir) {
                            eprintln!("Failed to create directory {}: {}", dest_dir.display(), e);
                            continue;
                        }
                        
                        // Create relative path from dest_file to source file  
                        let relative_path = Path::new("../../../data_raw")
                            .join(&source_name)
                            .join(path.file_name().unwrap());
                        
                        if let Err(e) = unix_fs::symlink(&relative_path, &dest_file) {
                            eprintln!("Failed to create symlink for {}: {}", path.display(), e);
                        }
                    }
                    
                    processed_files += 1;
                }
                Err(e) => {
                    eprintln!("Error processing {}: {}", path.display(), e);
                    error_files += 1;
                }
            }
        }

        if stats.total > 0 {
            source_stats.insert(source_name, stats);
        }
    }

    // Generate markdown table
    println!("# MOL File Classification Results");
    println!();
    println!("Classification based on parser compatibility:");
    println!("- **Molecule**: Files that parse successfully with `parse_mol` (basic parser)");
    println!("- **MoleculeLike**: Files that require `parse_mol_moleculelike` (extended features)");
    println!("- **Invalid**: Files that fail both parsers");
    println!();
    println!("| Source | Total | Molecule | MoleculeLike | Invalid | Valid % |");
    println!("| --- | --- | --- | --- | --- | --- |");

    // Sort sources alphabetically for consistent output
    let mut sources: Vec<_> = source_stats.iter().collect();
    sources.sort_by_key(|(name, _)| name.as_str());

    let mut total_molecule = 0;
    let mut total_moleculelike = 0;
    let mut total_invalid = 0;
    let mut grand_total = 0;

    for (source, stats) in sources {
        println!(
            "| {} | {} | {} | {} | {} | {:.1}% |",
            source,
            stats.total,
            stats.molecule,
            stats.moleculelike,
            stats.invalid,
            stats.valid_percentage()
        );
        
        total_molecule += stats.molecule;
        total_moleculelike += stats.moleculelike;
        total_invalid += stats.invalid;
        grand_total += stats.total;
    }

    // Add totals row
    let overall_valid_percentage = if grand_total == 0 {
        0.0
    } else {
        ((total_molecule + total_moleculelike) as f64 / grand_total as f64) * 100.0
    };

    println!("| **Total** | **{}** | **{}** | **{}** | **{}** | **{:.1}%** |", 
             grand_total, total_molecule, total_moleculelike, total_invalid, overall_valid_percentage);
    
    println!();
    println!();
    println!("**Summary:**");
    println!("- Processed: {}/{} files", processed_files, total_files);
    if error_files > 0 {
        println!("- Errors: {} files could not be read", error_files);
    }
    println!("- {:.1}% of files are valid (parseable)", overall_valid_percentage);
    println!("- {:.1}% of files work with basic parser (`parse_mol`)", (total_molecule as f64 / grand_total as f64) * 100.0);
    
    if !should_sort {
        println!();
        println!("**Note:** Run with `--sort` flag to create symlinks in data/ directories:");
        println!("```");
        println!("cargo run --bin mol_classifier -- --sort");
        println!("```");
    } else {
        println!();
        println!("**Files sorted into data/ directories via symlinks**");
    }

    Ok(())
}