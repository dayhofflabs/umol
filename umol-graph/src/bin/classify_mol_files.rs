//! MOL file classifier for conformance test organization.
//!
//! Classifies MOL files into categories based on which parser configurations succeed:
//! - molecule: passes all 4 parsers
//! - molecule_lenient: needs lenient flags for basic parser  
//! - extended_molecule: needs extended parser
//! - extended_molecule_lenient: needs lenient flags for extended parser
//! - invalid: fails all parsers
//! - bug: hierarchy violation (indicates parser inconsistency)
//!
//! Usage:
//!   cargo run --bin classify_mol_files           # Show classification stats
//!   cargo run --bin classify_mol_files -- --sort # Copy files to data/ directories

use std::collections::HashMap;
use std::error::Error;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::{env, fs, process};

use umol_io::io::ctfile::config::CtfileIoConfig;
use umol_io::io::ctfile::parser::{
    parse_extended_mol_bytes_with, parse_mol_bytes_to_table_ir_with,
};

/// Results from running all 4 parsers on a file
#[derive(Debug, Clone, Copy)]
struct ParseResults {
    mol_basic: bool,
    mol_lenient: bool,
    ext_extended: bool,
    ext_lenient: bool,
}

impl ParseResults {
    fn pattern(&self) -> String {
        let chars = [
            if self.mol_basic {
                "+".to_string()
            } else {
                "-".to_string()
            },
            if self.mol_lenient {
                "+".to_string()
            } else {
                "-".to_string()
            },
            if self.ext_extended {
                "+".to_string()
            } else {
                "-".to_string()
            },
            if self.ext_lenient {
                "+".to_string()
            } else {
                "-".to_string()
            },
        ];
        chars.join("")
    }

    fn has_hierarchy_violation(&self) -> bool {
        // mol(basic) → mol(lenient)
        if self.mol_basic && !self.mol_lenient {
            println!(
                "mol_basic ({}) -> mol_lenient ({})",
                self.mol_basic, self.mol_lenient
            );
            return true;
        }
        // mol(basic) → extended(extended)
        if self.mol_basic && !self.ext_extended {
            println!(
                "mol_basic ({}) -> ext_extended ({})",
                self.mol_basic, self.ext_extended
            );
            return true;
        }
        // extended(extended) → extended(lenient)
        if self.ext_extended && !self.ext_lenient {
            println!(
                "ext_extended ({}) -> ext_lenient ({})",
                self.ext_extended, self.ext_lenient
            );
            return true;
        }
        false
    }

    fn violation_description(&self) -> Option<String> {
        let mut violations = Vec::new();
        if self.mol_basic && !self.mol_lenient {
            violations.push("mol(basic) succeeded but mol(lenient) failed");
        }
        if self.mol_basic && !self.ext_extended {
            violations.push("mol(basic) succeeded but extended(extended) failed");
        }
        if self.ext_extended && !self.ext_lenient {
            violations.push("extended(extended) succeeded but extended(lenient) failed");
        }
        if violations.is_empty() {
            None
        } else {
            Some(violations.join("; "))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Category {
    Molecule,
    MoleculeLenient,
    ExtendedMolecule,
    ExtendedMoleculeLenient,
    Invalid,
    Bug,
}

impl Category {
    fn from_results(results: &ParseResults) -> Self {
        if results.has_hierarchy_violation() {
            return Category::Bug;
        }
        match &results.pattern()[..] {
            "++++" => Category::Molecule,
            "-+-+" => Category::MoleculeLenient,
            "--++" => Category::ExtendedMolecule,
            "---+" => Category::ExtendedMoleculeLenient,
            "----" => Category::Invalid,
            _ => {
                println!("PATTERN {}", results.pattern());
                Category::Bug
            }
        }
    }

    fn dir_name(&self) -> &'static str {
        match self {
            Category::Molecule => "molecule",
            Category::MoleculeLenient => "molecule_lenient",
            Category::ExtendedMolecule => "extended_molecule",
            Category::ExtendedMoleculeLenient => "extended_molecule_lenient",
            Category::Invalid => "invalid",
            Category::Bug => "bug",
        }
    }
}

#[derive(Debug, Default)]
struct ClassificationStats {
    molecule: usize,
    molecule_lenient: usize,
    extended_molecule: usize,
    extended_molecule_lenient: usize,
    invalid: usize,
    bug: usize,
    total: usize,
}

impl ClassificationStats {
    fn add(&mut self, category: Category) {
        match category {
            Category::Molecule => self.molecule += 1,
            Category::MoleculeLenient => self.molecule_lenient += 1,
            Category::ExtendedMolecule => self.extended_molecule += 1,
            Category::ExtendedMoleculeLenient => self.extended_molecule_lenient += 1,
            Category::Invalid => self.invalid += 1,
            Category::Bug => self.bug += 1,
        }
        self.total += 1;
    }

    fn valid_count(&self) -> usize {
        self.molecule
            + self.molecule_lenient
            + self.extended_molecule
            + self.extended_molecule_lenient
    }

    fn valid_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            (self.valid_count() as f64 / self.total as f64) * 100.0
        }
    }
}

fn classify_mol_file(file_path: &Path) -> Result<(Category, ParseResults), Box<dyn Error>> {
    let mol_bytes = fs::read(file_path)?;

    // For basic parser: BASIC vs BASIC LENIENT (lenient parsing features, no extended atoms/bonds)
    let mol_basic_config = CtfileIoConfig::basic();
    let mol_lenient_config = CtfileIoConfig::basic_lenient();

    // For extended parser: EXTENDED vs LENIENT (includes all extended features)
    let ext_extended_config = CtfileIoConfig::extended();
    let ext_lenient_config = CtfileIoConfig::extended_lenient();

    let results = ParseResults {
        mol_basic: parse_mol_bytes_to_table_ir_with(&mol_bytes, &mol_basic_config).is_ok(),
        mol_lenient: parse_mol_bytes_to_table_ir_with(&mol_bytes, &mol_lenient_config).is_ok(),
        ext_extended: parse_extended_mol_bytes_with(&mol_bytes, &ext_extended_config).is_ok(),
        ext_lenient: parse_extended_mol_bytes_with(&mol_bytes, &ext_lenient_config).is_ok(),
    };

    let category = Category::from_results(&results);
    Ok((category, results))
}

fn clean_existing_files(data_path: &str) -> Result<(), Box<dyn Error>> {
    let categories = [
        "molecule",
        "molecule_lenient",
        "extended_molecule",
        "extended_molecule_lenient",
        "invalid",
        "bug",
    ];

    for category in categories {
        let category_path = Path::new(data_path).join(category);
        if !category_path.exists() {
            continue;
        }

        for entry in fs::read_dir(&category_path)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                println!("Removing {}", path.display());
                if let Err(e) = fs::remove_dir_all(&path) {
                    eprintln!("Warning: Failed to remove {}: {}", path.display(), e);
                }
            }
        }
    }
    Ok(())
}

fn collect_mol_files(dir: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() {
            files.extend(collect_mol_files(&path)?);
        } else if path.extension() == Some(OsStr::new("mol")) {
            files.push(path);
        }
    }
    Ok(files)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if !(args.len() == 1 || args.len() == 2 && args[1] == "--sort") {
        eprintln!(
            "Usage: cargo run --bin classify_mol_files           # Show classification stats"
        );
        eprintln!(
            "       cargo run --bin classify_mol_files -- --sort # Copy files to data/ directories"
        );
        process::exit(1);
    }
    let should_sort = args.iter().any(|a| a == "--sort");

    // Paths relative to workspace root, prefixed with package dir
    let base_path = "umol-models-graph";
    let data_raw_path = format!("{}/tests/mol_parsing/data_raw", base_path);
    let data_path = format!("{}/tests/mol_parsing/data", base_path);

    if !Path::new(&data_raw_path).exists() {
        eprintln!("Error: {} directory not found", data_raw_path);
        process::exit(1);
    }

    if should_sort {
        println!("Cleaning existing organized files...");
        clean_existing_files(&data_path)?;
    }

    let mut source_stats: HashMap<String, ClassificationStats> = HashMap::new();
    let mut bug_files: Vec<(String, String, ParseResults)> = Vec::new();
    let mut total_files = 0;
    let mut processed_files = 0;
    let mut error_files = 0;

    for entry in fs::read_dir(&data_raw_path)? {
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

        let mut stats = ClassificationStats::default();
        let mol_files = collect_mol_files(&source_path)?;

        for path in mol_files {
            total_files += 1;

            match classify_mol_file(&path) {
                Ok((category, results)) => {
                    stats.add(category);

                    if category == Category::Bug {
                        let filename = path.file_name().unwrap().to_string_lossy().to_string();
                        bug_files.push((source_name.clone(), filename, results));
                    }

                    if should_sort {
                        let dest_dir = Path::new(&data_path)
                            .join(category.dir_name())
                            .join(&source_name);
                        let dest_file = dest_dir.join(path.file_name().unwrap());

                        if let Err(e) = fs::create_dir_all(&dest_dir) {
                            eprintln!("Failed to create directory {}: {}", dest_dir.display(), e);
                            continue;
                        }

                        if let Err(e) = fs::copy(&path, &dest_file) {
                            eprintln!(
                                "Failed to copy file from {} to {}: {}",
                                path.display(),
                                dest_file.display(),
                                e
                            );
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

    // Print markdown report
    println!("# MOL File Classification Results\n");
    println!("Classification based on parser compatibility:");
    println!("- **molecule** (++++): All parsers succeed");
    println!("- **molecule_lenient** (-+-+): Needs lenient flags for basic parser");
    println!("- **extended_molecule** (--++): Needs extended parser");
    println!("- **extended_molecule_lenient** (---+): Needs lenient flags for extended parser");
    println!("- **invalid** (----): All parsers fail");
    println!("- **bug**: Parser hierarchy violation\n");

    println!("| Source | Total | molecule | mol_lenient | ext_mol | ext_lenient | invalid | bug | Valid % |");
    println!("| --- | --- | --- | --- | --- | --- | --- | --- | --- |");

    let mut sources: Vec<_> = source_stats.iter().collect();
    sources.sort_by_key(|(name, _)| name.as_str());

    let mut totals = ClassificationStats::default();

    for (source, stats) in &sources {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {:.1}% |",
            source,
            stats.total,
            stats.molecule,
            stats.molecule_lenient,
            stats.extended_molecule,
            stats.extended_molecule_lenient,
            stats.invalid,
            stats.bug,
            stats.valid_percentage()
        );

        totals.molecule += stats.molecule;
        totals.molecule_lenient += stats.molecule_lenient;
        totals.extended_molecule += stats.extended_molecule;
        totals.extended_molecule_lenient += stats.extended_molecule_lenient;
        totals.invalid += stats.invalid;
        totals.bug += stats.bug;
        totals.total += stats.total;
    }

    println!(
        "| **Total** | **{}** | **{}** | **{}** | **{}** | **{}** | **{}** | **{}** | **{:.1}%** |",
        totals.total,
        totals.molecule,
        totals.molecule_lenient,
        totals.extended_molecule,
        totals.extended_molecule_lenient,
        totals.invalid,
        totals.bug,
        totals.valid_percentage()
    );

    println!("\n\n**Summary:**");
    println!("- Processed: {}/{} files", processed_files, total_files);
    if error_files > 0 {
        println!("- Errors: {} files could not be read", error_files);
    }
    println!(
        "- {:.1}% of files are valid (parseable)",
        totals.valid_percentage()
    );

    if totals.bug > 0 {
        println!(
            "\n⚠️  **WARNING: {} files have parser hierarchy violations!**\n",
            totals.bug
        );
        println!("| Source | File | Pattern | Violation |");
        println!("| --- | --- | --- | --- |");
        for (source, filename, results) in &bug_files {
            println!(
                "| {} | {} | {} | {} |",
                source,
                filename,
                results.pattern(),
                results.violation_description().unwrap_or_default()
            );
        }
    }

    if !should_sort {
        println!("\n**Note:** Run with `--sort` flag to copy files to data/ directories:");
        println!("```");
        println!("cargo run --bin classify_mol_files -- --sort");
        println!("```");
    } else {
        println!("\n**Files sorted into data/ directories**");
    }

    Ok(())
}
