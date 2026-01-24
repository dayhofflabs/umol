//! SMILES string classifier for conformance test organization.
//!
//! Classifies SMILES strings into categories based on parser success:
//! - opensmiles_strict: passes strict OpenSMILES parser
//! - invalid: fails parser
//!
//! Usage:
//!   cargo run --bin classify_smiles_strings                    # Show classification stats
//!   cargo run --bin classify_smiles_strings -- --sort          # Sort into data/ directories
//!   cargo run --bin classify_smiles_strings -- --sort --max-samples 100

use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use clap::Parser;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use umol_models_graph::io::smiles::parse_smiles;

#[derive(Parser)]
#[command(name = "classify_smiles_strings")]
#[command(about = "Classify SMILES strings by parser compatibility")]
struct Args {
    /// Sort files into data/ directories
    #[arg(long)]
    sort: bool,

    /// Maximum samples per file (0 = no limit)
    #[arg(long, default_value = "200")]
    max_samples: usize,

    /// Random seed for sampling
    #[arg(long, default_value = "0")]
    seed: u64,
}

#[derive(Debug, Clone)]
struct SmilesEntry {
    smiles: String,
    source_file: String,
    line_number: usize,
}


fn read_smiles_file(path: &Path) -> io::Result<Vec<SmilesEntry>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();
    let mut line_number = 0;

    let source_file = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string_lossy().to_string());

    for line in reader.lines() {
        let line = line?;
        line_number += 1;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            continue;
        }
        if trimmed.starts_with("SMILES")
            || trimmed.starts_with("smiles")
            || trimmed.starts_with("Smiles")
        {
            continue;
        }

        let smiles = if let Some(tab_pos) = trimmed.find('\t') {
            &trimmed[..tab_pos]
        } else if let Some(space_pos) = trimmed.find(' ') {
            &trimmed[..space_pos]
        } else {
            trimmed
        };

        entries.push(SmilesEntry {
            smiles: smiles.to_string(),
            source_file: source_file.clone(),
            line_number,
        });
    }

    Ok(entries)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Category {
    OpensmilesStrict,
    Invalid,
}

impl Category {
    fn dir_name(&self) -> &'static str {
        match self {
            Category::OpensmilesStrict => "opensmiles_strict",
            Category::Invalid => "invalid",
        }
    }
}

#[derive(Debug, Default)]
struct ClassificationStats {
    opensmiles_strict: usize,
    invalid: usize,
    total: usize,
}

impl ClassificationStats {
    fn add(&mut self, category: Category) {
        match category {
            Category::OpensmilesStrict => self.opensmiles_strict += 1,
            Category::Invalid => self.invalid += 1,
        }
        self.total += 1;
    }

    fn valid_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            (self.opensmiles_strict as f64 / self.total as f64) * 100.0
        }
    }
}

fn classify_smiles(smiles: &str) -> Category {
    match parse_smiles(smiles.as_bytes()) {
        Ok(_) => Category::OpensmilesStrict,
        Err(_) => Category::Invalid,
    }
}

fn collect_smi_files(dir: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() {
            files.extend(collect_smi_files(&path)?);
        } else {
            let ext = path.extension().and_then(|e| e.to_str());
            if ext == Some("smi") || ext == Some("smiles") {
                files.push(path);
            }
        }
    }
    Ok(files)
}

fn clean_existing_files(data_path: &str) -> Result<(), Box<dyn Error>> {
    let categories = ["opensmiles_strict", "invalid"];

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

fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' | ' ' | '\t' | '\n' | '\r' => '_',
            c if c.is_ascii_control() => '_',
            c => c,
        })
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let base_path = "umol-models-graph";
    let data_raw_path = format!("{}/tests/smiles_parsing/data_raw", base_path);
    let data_path = format!("{}/tests/smiles_parsing/data", base_path);

    if !Path::new(&data_raw_path).exists() {
        eprintln!("Error: {} directory not found", data_raw_path);
        std::process::exit(1);
    }

    if args.sort {
        println!("Cleaning existing organized files...");
        clean_existing_files(&data_path)?;
    }

    let mut rng = ChaCha8Rng::seed_from_u64(args.seed);
    let mut source_stats: HashMap<String, ClassificationStats> = HashMap::new();
    let mut processed_smiles = 0usize;
    let mut error_files = 0usize;

    // Collect entries by category for sorting
    let mut categorized: HashMap<String, HashMap<Category, Vec<SmilesEntry>>> = HashMap::new();

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
        let smi_files = collect_smi_files(&source_path)?;
        let mut source_entries: HashMap<Category, Vec<SmilesEntry>> = HashMap::new();

        for file_path in &smi_files {
            let entries = match read_smiles_file(file_path) {
                Ok(e) => e,
                Err(e) => {
                    eprintln!("Error reading {}: {}", file_path.display(), e);
                    error_files += 1;
                    continue;
                }
            };

            // Sample per-file if needed
            let entries_to_process = if args.max_samples > 0 && entries.len() > args.max_samples {
                let mut sampled = entries;
                sampled.shuffle(&mut rng);
                sampled.truncate(args.max_samples);
                sampled
            } else {
                entries
            };

            for entry in entries_to_process {
                let category = classify_smiles(&entry.smiles);
                stats.add(category);

                if args.sort {
                    source_entries
                        .entry(category)
                        .or_default()
                        .push(entry.clone());
                }

                processed_smiles += 1;
            }
        }

        if stats.total > 0 {
            source_stats.insert(source_name.clone(), stats);
            if args.sort {
                categorized.insert(source_name, source_entries);
            }
        }
    }

    // Write sorted files - one .smiles file per SMILES
    if args.sort {
        for (source_name, entries_by_category) in &categorized {
            for (category, entries) in entries_by_category {
                let dest_dir = Path::new(&data_path)
                    .join(category.dir_name())
                    .join(source_name);

                fs::create_dir_all(&dest_dir)?;

                for (idx, entry) in entries.iter().enumerate() {
                    // Create filename from index and sanitized SMILES prefix (max 20 chars)
                    let truncated = entry.smiles.len() > 20;
                    let smiles_prefix: String = entry.smiles.chars().take(20).collect();
                    let suffix = if truncated { "_" } else { "" };
                    let filename = format!(
                        "{:04}_{}{}.smiles",
                        idx,
                        sanitize_filename(&smiles_prefix),
                        suffix
                    );
                    let dest_file = dest_dir.join(&filename);

                    let mut file = File::create(&dest_file)?;
                    writeln!(file, "# {}:{}", entry.source_file, entry.line_number)?;
                    writeln!(file, "{}", entry.smiles)?;
                }
            }
        }
    }

    // Print markdown report
    println!("\n# SMILES Classification Results\n");
    println!("Classification based on parser compatibility:");
    println!("- **opensmiles_strict**: Passes strict OpenSMILES parser");
    println!("- **invalid**: Fails parser\n");

    if args.max_samples > 0 {
        println!(
            "**Sampling**: max {} per file, seed {}\n",
            args.max_samples, args.seed
        );
    }

    println!("| Source | Total | opensmiles_strict | invalid | Valid % |");
    println!("| --- | --- | --- | --- | --- |");

    let mut sources: Vec<_> = source_stats.iter().collect();
    sources.sort_by_key(|(name, _)| name.as_str());

    let mut totals = ClassificationStats::default();

    for (source, stats) in &sources {
        println!(
            "| {} | {} | {} | {} | {:.1}% |",
            source,
            stats.total,
            stats.opensmiles_strict,
            stats.invalid,
            stats.valid_percentage()
        );

        totals.opensmiles_strict += stats.opensmiles_strict;
        totals.invalid += stats.invalid;
        totals.total += stats.total;
    }

    println!(
        "| **Total** | **{}** | **{}** | **{}** | **{:.1}%** |",
        totals.total,
        totals.opensmiles_strict,
        totals.invalid,
        totals.valid_percentage()
    );

    println!("\n\n**Summary:**");
    println!("- Processed: {} SMILES strings", processed_smiles);
    if error_files > 0 {
        println!("- Errors: {} files could not be read", error_files);
    }
    println!(
        "- {:.1}% of SMILES are valid (parseable with strict OpenSMILES)",
        totals.valid_percentage()
    );

    if !args.sort {
        println!("\n**Note:** Run with `--sort` flag to create individual .smiles files:");
        println!("```");
        println!("cargo run --bin classify_smiles_strings -- --sort");
        println!("cargo run --bin classify_smiles_strings -- --sort --max-samples 100 --seed 42");
        println!("```");
    } else {
        println!("\n**Files sorted into data/ directories (one .smiles file per SMILES)**");
    }

    Ok(())
}
