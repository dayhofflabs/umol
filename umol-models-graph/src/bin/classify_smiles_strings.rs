//! SMILES string classifier for conformance test organization.
//!
//! Classifies SMILES strings into categories based on parser success:
//! - basic_opensmiles: passes strict OpenSMILES parser (no wildcards)
//! - opensmiles: passes extended OpenSMILES parser (with wildcards)
//! - invalid: fails all parsers
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
use std::process;

use clap::Parser;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use umol_models_graph::io::smiles::config::SmilesIoConfig;
use umol_models_graph::io::smiles::{parse_extended_smiles_bytes_with, parse_smiles};

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

    /// Number of rows to check for column detection
    #[arg(long, default_value = "5")]
    probe_rows: usize,
}

#[derive(Debug, Clone)]
struct SmilesEntry {
    smiles: String,
    source_file: String,
    line_number: usize,
}

/// Detected file format
#[derive(Debug, Clone)]
struct FileFormat {
    delimiter: Option<char>,
    smiles_column: usize,
}

impl Default for FileFormat {
    fn default() -> Self {
        Self {
            delimiter: None,
            smiles_column: 0,
        }
    }
}

/// Check if a string looks like valid SMILES using lenient parser
fn is_likely_smiles(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    let config = SmilesIoConfig::lenient();
    parse_extended_smiles_bytes_with(s.as_bytes(), &config).is_ok()
}

/// Split a line by delimiter
fn split_line(line: &str, delimiter: Option<char>) -> Vec<&str> {
    match delimiter {
        Some(d) => line.split(d).collect(),
        None => vec![line],
    }
}

/// Extract SMILES from a line, handling CXSMILES format.
/// CXSMILES has format: `SMILES |...|` or `SMILES |...| other_data`
/// Returns the SMILES/CXSMILES portion.
fn extract_cxsmiles(line: &str) -> Option<&str> {
    // Look for ` |` which indicates start of CXSMILES extension
    if let Some(pipe_start) = line.find(" |") {
        // Find the closing `|`
        let extension_start = pipe_start + 2; // skip " |"
        if let Some(rel_end) = line[extension_start..].find('|') {
            let end = extension_start + rel_end + 1; // include closing |
            return Some(&line[..end]);
        }
    }
    None
}

/// Detect file format by probing first N non-empty, non-comment lines
fn detect_file_format(path: &Path, probe_rows: usize) -> io::Result<FileFormat> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut probe_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Skip header lines
        if trimmed.to_lowercase().starts_with("smiles") {
            continue;
        }

        probe_lines.push(trimmed.to_string());
        if probe_lines.len() >= probe_rows {
            break;
        }
    }

    if probe_lines.is_empty() {
        return Ok(FileFormat::default());
    }

    // Detect delimiter from first line (tab or comma only, not space - space is handled by CXSMILES extraction)
    let delimiter = if probe_lines[0].contains('\t') {
        Some('\t')
    } else if probe_lines[0].contains(',') {
        Some(',')
    } else {
        None
    };

    // Single column or no delimiter, SMILES is first
    if delimiter.is_none() {
        return Ok(FileFormat {
            delimiter: None,
            smiles_column: 0,
        });
    }

    // Determine number of columns
    let num_columns = split_line(&probe_lines[0], delimiter).len();

    if num_columns == 1 {
        return Ok(FileFormat {
            delimiter,
            smiles_column: 0,
        });
    }

    // Try each column to find SMILES
    let mut column_scores: Vec<usize> = vec![0; num_columns];

    for line in &probe_lines {
        let parts = split_line(line, delimiter);
        for (col_idx, part) in parts.iter().enumerate() {
            if col_idx < num_columns && is_likely_smiles(part.trim()) {
                column_scores[col_idx] += 1;
            }
        }
    }

    // Pick column with highest score
    let smiles_column = column_scores
        .iter()
        .enumerate()
        .max_by_key(|(_, &score)| score)
        .map(|(idx, _)| idx)
        .unwrap_or(0);

    Ok(FileFormat {
        delimiter,
        smiles_column,
    })
}

fn read_smiles_file(path: &Path, probe_rows: usize) -> io::Result<Vec<SmilesEntry>> {
    let format = detect_file_format(path, probe_rows)?;

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();

    let source_file = path
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string_lossy().to_string());

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Skip header lines
        if trimmed.to_lowercase().starts_with("smiles") {
            continue;
        }

        // First try CXSMILES extraction (handles space before |...|)
        let smiles = if let Some(cxsmiles) = extract_cxsmiles(trimmed) {
            cxsmiles
        } else if let Some(delimiter) = format.delimiter {
            // Use column extraction for tab/comma delimited
            let parts = split_line(trimmed, Some(delimiter));
            parts
                .get(format.smiles_column)
                .map(|s| s.trim())
                .unwrap_or(trimmed)
        } else {
            // No delimiter detected, take everything before first space (or whole line)
            trimmed.split_once(' ').map(|(s, _)| s).unwrap_or(trimmed)
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
    BasicOpensmiles,
    Opensmiles,
    Invalid,
    Bug,
}

impl Category {
    fn dir_name(&self) -> &'static str {
        match self {
            Category::BasicOpensmiles => "basic_opensmiles",
            Category::Opensmiles => "opensmiles",
            Category::Invalid => "invalid",
            Category::Bug => "bug",
        }
    }
}

#[derive(Debug, Default)]
struct ClassificationStats {
    basic_opensmiles: usize,
    opensmiles: usize,
    invalid: usize,
    bug: usize,
    total: usize,
}

impl ClassificationStats {
    fn add(&mut self, category: Category) {
        match category {
            Category::BasicOpensmiles => self.basic_opensmiles += 1,
            Category::Opensmiles => self.opensmiles += 1,
            Category::Invalid => self.invalid += 1,
            Category::Bug => self.bug += 1,
        }
        self.total += 1;
    }

    fn valid_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            ((self.basic_opensmiles + self.opensmiles) as f64 / self.total as f64) * 100.0
        }
    }
}

fn classify_smiles(smiles: &str) -> Category {
    // Run both parsers to enforce parser hierarchy
    // extended parser must be a superset of basic parser
    let basic_ok = parse_smiles(smiles).is_ok();
    let config = SmilesIoConfig::opensmiles();
    let extended_ok = parse_extended_smiles_bytes_with(smiles.as_bytes(), &config).is_ok();

    match (basic_ok, extended_ok) {
        (true, true) => Category::BasicOpensmiles,
        (false, true) => Category::Opensmiles,
        (false, false) => Category::Invalid,
        (true, false) => Category::Bug, // hierarchy violation
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
    let categories = ["basic_opensmiles", "opensmiles", "invalid", "bug"];

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
        process::exit(1);
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
            let entries = match read_smiles_file(file_path, args.probe_rows) {
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
    println!("- **basic_opensmiles**: Passes strict OpenSMILES parser (no wildcards)");
    println!("- **opensmiles**: Passes extended OpenSMILES parser (with wildcards)");
    println!("- **invalid**: Fails all parsers\n");

    if args.max_samples > 0 {
        println!(
            "**Sampling**: max {} per file, seed {}\n",
            args.max_samples, args.seed
        );
    }

    println!("| Source | Total | basic_opensmiles | opensmiles | invalid | bug | Valid % |");
    println!("| --- | --- | --- | --- | --- | --- | --- |");

    let mut sources: Vec<_> = source_stats.iter().collect();
    sources.sort_by_key(|(name, _)| name.as_str());

    let mut totals = ClassificationStats::default();

    for (source, stats) in &sources {
        println!(
            "| {} | {} | {} | {} | {} | {} | {:.1}% |",
            source,
            stats.total,
            stats.basic_opensmiles,
            stats.opensmiles,
            stats.invalid,
            stats.bug,
            stats.valid_percentage()
        );

        totals.basic_opensmiles += stats.basic_opensmiles;
        totals.opensmiles += stats.opensmiles;
        totals.invalid += stats.invalid;
        totals.bug += stats.bug;
        totals.total += stats.total;
    }

    println!(
        "| **Total** | **{}** | **{}** | **{}** | **{}** | **{}** | **{:.1}%** |",
        totals.total,
        totals.basic_opensmiles,
        totals.opensmiles,
        totals.invalid,
        totals.bug,
        totals.valid_percentage()
    );

    println!("\n\n**Summary:**");
    println!("- Processed: {} SMILES strings", processed_smiles);
    if error_files > 0 {
        println!("- Errors: {} files could not be read", error_files);
    }
    println!(
        "- {:.1}% of SMILES are valid (parseable with OpenSMILES)",
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
