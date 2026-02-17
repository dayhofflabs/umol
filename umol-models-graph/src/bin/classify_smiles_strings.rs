//! SMILES string classifier for conformance test organization.
//!
//! Classifies SMILES strings into categories based on parser success:
//! - basic_opensmiles: passes strict OpenSMILES parser (no wildcards)
//! - opensmiles: passes extended OpenSMILES parser (with wildcards)
//! - basic_chemaxon: passes basic Chemaxon parser (no wildcards)
//! - chemaxon: passes extended Chemaxon parser (with wildcards)
//! - chemaxon_invalid: SMILES part parses, but CX block is invalid/unhandled
//! - invalid: fails all parsers
//! - bug: unexpected parser outcome (indicates parser inconsistency)
//!
//! Usage:
//!   cargo run --bin classify_smiles_strings                    # Show classification stats
//!   cargo run --bin classify_smiles_strings -- --sort          # Sort into data/ directories
//!   cargo run --bin classify_smiles_strings -- --sort --max-samples 100

use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Cursor, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::LazyLock;

use clap::Parser;
use murmur3::murmur3_32;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use regex::Regex;
use umol_models_graph::io::smiles::config::SmilesIoConfig;
use umol_models_graph::io::smiles::{parse_extended_smiles_bytes_with, parse_smiles_bytes_with};

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
#[derive(Debug, Clone, Default)]
struct FileFormat {
    delimiter: Option<char>,
    smiles_column: usize,
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

struct ParseResults {
    has_cx_annotations: bool,
    basic_opensmiles: bool,
    opensmiles: bool,
    basic_chemaxon: bool,
    chemaxon: bool,
}

impl ParseResults {
    fn pattern(&self) -> String {
        format!(
            "{}{}{}{}{}",
            if self.has_cx_annotations { "+" } else { "-" },
            if self.basic_opensmiles { "+" } else { "-" },
            if self.basic_chemaxon { "+" } else { "-" },
            if self.opensmiles { "+" } else { "-" },
            if self.chemaxon { "+" } else { "-" },
        )
    }

    fn has_hierarchy_violation(&self) -> bool {
        // basic_opensmiles → opensmiles
        if self.basic_opensmiles && !self.opensmiles {
            println!(
                "basic_opensmiles ({}) -> opensmiles ({})",
                self.basic_opensmiles, self.opensmiles
            );
            return true;
        }
        // If there's no CX block, basic_chemaxon should behave like basic_opensmiles
        // and chemaxon should behave like opensmiles.
        if !self.has_cx_annotations && self.basic_chemaxon && !self.basic_opensmiles {
            println!(
                "basic_chemaxon ({}) -> basic_opensmiles ({})",
                self.basic_chemaxon, self.chemaxon
            );
            return true;
        }
        if !self.has_cx_annotations && self.opensmiles && !self.chemaxon {
            println!(
                "opensmiles ({}) -> chemaxon ({})",
                self.opensmiles, self.chemaxon
            );
            return true;
        }
        // basic_chemaxon → chemaxon
        if self.basic_chemaxon && !self.chemaxon {
            println!(
                "basic_chemaxon ({}) -> chemaxon ({})",
                self.basic_chemaxon, self.chemaxon
            );
            return true;
        }
        false
    }

    fn violation_description(&self) -> Option<String> {
        let mut violations = Vec::new();
        if self.basic_opensmiles && !self.opensmiles {
            violations.push("basic_opensmiles succeeded but opensmiles failed");
        }
        if !self.has_cx_annotations && self.opensmiles && !self.chemaxon {
            violations.push("opensmiles succeeded but chemaxon failed");
        }
        if self.basic_chemaxon && !self.chemaxon {
            violations.push("basic_chemaxon succeeded but chemaxon failed");
        }
        if violations.is_empty() {
            None
        } else {
            Some(violations.join("; "))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Category {
    BasicOpensmiles,
    Opensmiles,
    BasicChemaxon,
    Chemaxon,
    ChemaxonInvalid,
    Invalid,
    Bug,
}

impl Category {
    fn from_results(results: &ParseResults) -> Self {
        if results.has_hierarchy_violation() {
            return Category::Bug;
        }
        match &results.pattern()[..] {
            // Pattern: XABCD
            // X: has_cx
            // A: basic_opensmiles
            // B: basic_chemaxon
            // C: opensmiles
            // D: chemaxon
            "-++++" => Category::BasicOpensmiles,
            "---++" => Category::Opensmiles,
            "+++++" => Category::BasicChemaxon,
            "++-++" | "+--++" => Category::Chemaxon,
            "++-+-" | "+--+-" => Category::ChemaxonInvalid,
            "-----" | "+----" => Category::Invalid,
            _ => {
                println!("PATTERN {}", results.pattern());
                Category::Bug
            }
        }
    }

    fn dir_name(&self) -> &'static str {
        match self {
            Category::BasicOpensmiles => "basic_opensmiles",
            Category::Opensmiles => "opensmiles",
            Category::BasicChemaxon => "basic_chemaxon",
            Category::Chemaxon => "chemaxon",
            Category::Invalid => "invalid",
            Category::ChemaxonInvalid => "chemaxon_invalid",
            Category::Bug => "bug",
        }
    }
}

#[derive(Debug, Default)]
struct ClassificationStats {
    basic_opensmiles: usize,
    opensmiles: usize,
    basic_chemaxon: usize,
    chemaxon: usize,
    chemaxon_invalid: usize,
    invalid: usize,
    bug: usize,
    total: usize,
}

impl ClassificationStats {
    fn add(&mut self, category: Category) {
        match category {
            Category::BasicOpensmiles => self.basic_opensmiles += 1,
            Category::Opensmiles => self.opensmiles += 1,
            Category::BasicChemaxon => self.basic_chemaxon += 1,
            Category::Chemaxon => self.chemaxon += 1,
            Category::ChemaxonInvalid => self.chemaxon_invalid += 1,
            Category::Invalid => self.invalid += 1,
            Category::Bug => self.bug += 1,
        }
        self.total += 1;
    }

    fn valid_percentage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            let valid =
                self.basic_opensmiles + self.opensmiles + self.basic_chemaxon + self.chemaxon;
            (valid as f64 / self.total as f64) * 100.0
        }
    }
}

/// Check if SMILES has CX annotations (` |...|` block)
fn has_cx_annotations(smiles: &str) -> bool {
    static CX_ANNOTATIONS_RE: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^\S+\s+\|.*\|").expect("CX annotations regex"));
    CX_ANNOTATIONS_RE.is_match(smiles)
}

fn classify_smiles(smiles: &str) -> Result<(Category, ParseResults), Box<dyn Error>> {
    let basic_config = SmilesIoConfig::basic();
    let opensmiles_config = SmilesIoConfig::opensmiles();
    let basic_chemaxon_config = SmilesIoConfig::basic_chemaxon();
    let chemaxon_config = SmilesIoConfig::chemaxon();

    let results = ParseResults {
        has_cx_annotations: has_cx_annotations(smiles),
        basic_opensmiles: parse_smiles_bytes_with(smiles.as_bytes(), &basic_config).is_ok(),
        opensmiles: parse_extended_smiles_bytes_with(smiles.as_bytes(), &opensmiles_config).is_ok(),
        basic_chemaxon: parse_smiles_bytes_with(smiles.as_bytes(), &basic_chemaxon_config).is_ok(),
        chemaxon: parse_extended_smiles_bytes_with(smiles.as_bytes(), &chemaxon_config).is_ok(),
    };

    let category = Category::from_results(&results);
    Ok((category, results))
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
            if ext == Some("smi") || ext == Some("smiles") || ext == Some("cxsmi") {
                files.push(path);
            }
        }
    }
    Ok(files)
}

fn clean_existing_files(data_path: &str) -> Result<(), Box<dyn Error>> {
    let categories = [
        "basic_opensmiles",
        "opensmiles",
        "basic_chemaxon",
        "chemaxon",
        "chemaxon_invalid",
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
    let mut bug_smiles: Vec<(String, String, ParseResults)> = Vec::new();
    let mut processed_smiles = 0usize;
    let mut error_files = 0usize;
    let mut error_smiles = 0usize;

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
                match classify_smiles(&entry.smiles) {
                    Ok((category, results)) => {
                        stats.add(category);
                        if category == Category::Bug {
                            bug_smiles.push((source_name.clone(), entry.smiles.clone(), results));
                        }
                        if args.sort {
                            source_entries
                                .entry(category)
                                .or_default()
                                .push(entry.clone());
                        }
                    }
                    Err(e) => {
                        eprintln!("Error classifying {}: {}", entry.smiles, e);
                        error_smiles += 1;
                    }
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
                    // Create filename from index, sanitized SMILES prefix (max 20 chars),
                    // and a stable hash suffix to ensure uniqueness for different CX variants
                    let truncated = entry.smiles.len() > 20;
                    let smiles_prefix: String = entry.smiles.chars().take(20).collect();
                    let prefix_suffix = if truncated { "_" } else { "" };
                    let hash = murmur3_32(&mut Cursor::new(entry.smiles.as_bytes()), 0).unwrap();
                    let filename = format!(
                        "{:04}_{}{}_{:04x}.smiles",
                        idx,
                        sanitize_filename(&smiles_prefix),
                        prefix_suffix,
                        (hash & 0xFFFF) as u16
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
    println!("- **basic_chemaxon**: Requires CHEMAXON_EXTENSIONS (basic parser)");
    println!("- **chemaxon**: Requires CHEMAXON_EXTENSIONS (extended parser)");
    println!("- **invalid**: Fails all parsers\n");
    println!("- **chemaxon_invalid**: SMILES parses but CX block is invalid/unhandled");

    if args.max_samples > 0 {
        println!(
            "**Sampling**: max {} per file, seed {}\n",
            args.max_samples, args.seed
        );
    }

    println!("| Source | Total | basic_opensmiles | opensmiles | basic_chemaxon | chemaxon | invalid | chemaxon_invalid | bug | Valid % |");
    println!("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |");

    let mut sources: Vec<_> = source_stats.iter().collect();
    sources.sort_by_key(|(name, _)| name.as_str());

    let mut totals = ClassificationStats::default();

    for (source, stats) in &sources {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.1}% |",
            source,
            stats.total,
            stats.basic_opensmiles,
            stats.opensmiles,
            stats.basic_chemaxon,
            stats.chemaxon,
            stats.invalid,
            stats.chemaxon_invalid,
            stats.bug,
            stats.valid_percentage()
        );

        totals.basic_opensmiles += stats.basic_opensmiles;
        totals.opensmiles += stats.opensmiles;
        totals.basic_chemaxon += stats.basic_chemaxon;
        totals.chemaxon += stats.chemaxon;
        totals.chemaxon_invalid += stats.chemaxon_invalid;
        totals.invalid += stats.invalid;
        totals.bug += stats.bug;
        totals.total += stats.total;
    }

    println!(
        "| **Total** | **{}** | **{}** | **{}** | **{}** | **{}** | **{}** | **{}** | **{}** | **{:.1}%** |",
        totals.total,
        totals.basic_opensmiles,
        totals.opensmiles,
        totals.basic_chemaxon,
        totals.chemaxon,
        totals.chemaxon_invalid,
        totals.invalid,
        totals.bug,
        totals.valid_percentage()
    );

    println!("\n\n**Summary:**");
    println!("- Processed: {} SMILES strings", processed_smiles);
    if error_files > 0 {
        println!("- Errors: {} files could not be read", error_files);
    }
    if error_smiles > 0 {
        println!("- Errors: {} SMILES could not be classified", error_smiles);
    }
    println!("- {:.1}% of SMILES are valid", totals.valid_percentage());

    if !bug_smiles.is_empty() {
        println!(
            "\n**WARNING: {} SMILES have parser hierarchy violations!**",
            bug_smiles.len()
        );
        println!("| Source | SMILES | Pattern | Violation |");
        println!("| --- | --- | --- | --- |");
        for (source, smiles, results) in &bug_smiles {
            println!(
                "| {} | {} | {} | {} |",
                source,
                smiles,
                results.pattern(),
                results.violation_description().unwrap_or_default()
            );
        }
    }

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
