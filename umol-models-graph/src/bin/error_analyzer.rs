use std::collections::HashMap;
use std::path::Path;
use std::{env, fs};

use umol_models_graph::io::mol::parser::parse_mol_moleculelike;

fn main() {
    let args: Vec<String> = env::args().collect();
    let target_dir = if args.len() > 1 {
        &args[1]
    } else {
        "tests/mol_parsing/data/invalid"
    };

    let mut error_patterns = HashMap::new();
    let mut files_tested = 0;

    println!("Analyzing parsing errors in: {}", target_dir);
    println!("========================================");

    collect_mol_files(
        Path::new(target_dir),
        &mut error_patterns,
        &mut files_tested,
    );

    println!("\n=== ERROR SUMMARY ===");
    println!("Total files tested: {}", files_tested);
    println!();

    let mut sorted_errors: Vec<_> = error_patterns.iter().collect();
    sorted_errors.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    for (error_type, files) in sorted_errors {
        println!("{}: {} files", error_type, files.len());
        for file in files.iter().take(3) {
            println!("  - {}", file);
        }
        if files.len() > 3 {
            println!("  ... and {} more", files.len() - 3);
        }
        println!();
    }
}

fn collect_mol_files(
    dir: &Path,
    error_patterns: &mut HashMap<String, Vec<String>>,
    files_tested: &mut usize,
) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                collect_mol_files(&path, error_patterns, files_tested);
            } else if path.extension().map_or(false, |ext| ext == "mol") {
                *files_tested += 1;

                if let Ok(content) = fs::read(&path) {
                    match parse_mol_moleculelike(&content) {
                        Ok(_) => {
                            // This shouldn't happen for invalid files
                            let error_type = "UNEXPECTED_SUCCESS".to_string();
                            error_patterns
                                .entry(error_type)
                                .or_default()
                                .push(path.display().to_string());
                        }
                        Err(e) => {
                            let error_str = format!("{}", e);
                            let error_type = extract_error_pattern(&error_str);
                            error_patterns
                                .entry(error_type)
                                .or_default()
                                .push(path.display().to_string());
                        }
                    }
                }
            }
        }
    }
}

fn extract_error_pattern(error_str: &str) -> String {
    if error_str.contains("code: Eof") {
        "Eof (unexpected end of file)".to_string()
    } else if error_str.contains("code: MapRes") {
        "MapRes (mapping/conversion error)".to_string()
    } else if error_str.contains("code: Digit") {
        "Digit (expected digit, found other)".to_string()
    } else if error_str.contains("code: CrLf") {
        "CrLf (line ending issue)".to_string()
    } else if error_str.contains("code: IsNot") {
        "IsNot (format mismatch)".to_string()
    } else if error_str.contains("code: Tag") {
        "Tag (expected tag not found)".to_string()
    } else if error_str.contains("R-group label can only be applied") {
        "RGroup (R-group validation error)".to_string()
    } else {
        // Extract the actual error code if possible
        if let Some(start) = error_str.find("code: ") {
            let code_part = &error_str[start + 6..];
            if let Some(end) = code_part.find(" ") {
                format!("Other ({})", &code_part[..end])
            } else if let Some(end) = code_part.find("}") {
                format!("Other ({})", &code_part[..end])
            } else {
                "Other (unknown)".to_string()
            }
        } else {
            "Other (unknown)".to_string()
        }
    }
}
