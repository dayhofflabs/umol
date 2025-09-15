use std::fs;
use std::path::PathBuf;

use umol_models_graph::io::mol::parser::parse_mol_moleculelike;

#[derive(Debug, Clone)]
struct FileAnalysis {
    path: PathBuf,
    error_type: String,
    specific_issue: String,
    has_rgp: bool,
    has_rxn: bool,
    has_v2000: bool,
    file_size: usize,
}

fn analyze_file(path: &PathBuf) -> Option<FileAnalysis> {
    let content = fs::read(path).ok()?;
    let content_str = String::from_utf8_lossy(&content);

    let has_rgp = content_str.contains("$RGP");
    let has_rxn = content_str.contains("$RXN");
    let has_v2000 = content_str.contains("V2000");
    let file_size = content.len();

    let relative_path = path
        .strip_prefix("tests/mol_parsing/data/invalid/")
        .unwrap()
        .to_path_buf();

    match parse_mol_moleculelike(&content) {
        Err(e) => {
            let error_str = format!("{:?}", e);
            let error_type = error_str.split('(').next().unwrap_or("Unknown").to_string();

            // Analyze specific issues
            let specific_issue = match error_type.as_str() {
                "Data" => {
                    if error_str.contains("Eof") {
                        if !has_v2000 {
                            "Missing V2000 tag".to_string()
                        } else if file_size < 50 {
                            "Truncated/malformed file".to_string()
                        } else {
                            "Unexpected end of file".to_string()
                        }
                    } else if error_str.contains("Digit") {
                        if let Some(counts_line) = content_str.lines().nth(3) {
                            if counts_line.contains("999V2000") || !counts_line.contains(' ') {
                                "Missing space before V2000".to_string()
                            } else if !counts_line
                                .chars()
                                .take(30)
                                .all(|c| c.is_ascii_digit() || c.is_ascii_whitespace())
                            {
                                "Non-digit in counts line".to_string()
                            } else {
                                "Digit parsing error".to_string()
                            }
                        } else {
                            "Counts line missing".to_string()
                        }
                    } else if error_str.contains("MapRes") {
                        // Check for unusual bond orders
                        let mut issues = Vec::new();
                        for line in content_str.lines().skip(4) {
                            if line.trim().is_empty() || line.starts_with("M  END") {
                                break;
                            }

                            // Bond line check
                            let parts: Vec<&str> = line.split_whitespace().collect();
                            if parts.len() >= 3
                                && parts.iter().take(3).all(|p| p.parse::<i32>().is_ok())
                            {
                                if let Ok(bond_order) = parts[2].parse::<u8>() {
                                    if bond_order > 4 || bond_order == 0 {
                                        issues.push(format!("Bond order {}", bond_order));
                                    }
                                }
                            }

                            // Atom line check (look for unusual elements)
                            if parts.len() >= 4 && parts[3].len() <= 3 {
                                let element = parts[3];
                                if element.starts_with('R')
                                    || element == "L"
                                    || element == "A"
                                    || element == "Q"
                                    || element == "D"
                                    || element == "T"
                                    || element.len() > 2
                                {
                                    issues.push(format!("Element '{}'", element));
                                }
                            }
                        }
                        if issues.is_empty() {
                            "Unknown mapping error".to_string()
                        } else {
                            issues.join(", ")
                        }
                    } else if error_str.contains("IsNot") {
                        if has_rxn {
                            "RXN file (not MOL)".to_string()
                        } else if has_rgp {
                            "Contains $RGP blocks".to_string()
                        } else {
                            "Format mismatch".to_string()
                        }
                    } else {
                        "Unknown data error".to_string()
                    }
                }
                _ => error_type.clone(),
            };

            Some(FileAnalysis {
                path: relative_path,
                error_type,
                specific_issue,
                has_rgp,
                has_rxn,
                has_v2000,
                file_size,
            })
        }
        Ok(_) => None, // File parsed successfully
    }
}

fn main() {
    let invalid_dir = PathBuf::from("tests/mol_parsing/data/invalid");
    let mut analyses = Vec::new();

    // Collect all analyses
    for entry in fs::read_dir(&invalid_dir).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_dir() {
            let subdir = entry.path();
            for file_entry in fs::read_dir(&subdir).unwrap() {
                let file_entry = file_entry.unwrap();
                if file_entry
                    .path()
                    .extension()
                    .is_some_and(|ext| ext == "mol")
                {
                    if let Some(analysis) = analyze_file(&file_entry.path()) {
                        analyses.push(analysis);
                    }
                }
            }
        }
    }

    // Sort by directory then filename
    analyses.sort_by(|a, b| a.path.cmp(&b.path));

    println!(
        "=== DETAILED FILE-BY-FILE ANALYSIS ({} files) ===\n",
        analyses.len()
    );

    for analysis in &analyses {
        println!("📁 {}", analysis.path.display());
        println!(
            "   Error: {} - {}",
            analysis.error_type, analysis.specific_issue
        );
        println!(
            "   Features: $RGP={}, $RXN={}, V2000={}, Size={}b",
            analysis.has_rgp, analysis.has_rxn, analysis.has_v2000, analysis.file_size
        );
        println!();
    }

    // Summary statistics
    let rgp_count = analyses.iter().filter(|a| a.has_rgp).count();
    let rxn_count = analyses.iter().filter(|a| a.has_rxn).count();
    let missing_v2000 = analyses.iter().filter(|a| !a.has_v2000).count();
    let empty_files = analyses.iter().filter(|a| a.file_size < 10).count();

    println!("=== SUMMARY ===");
    println!("Files with $RGP blocks: {}", rgp_count);
    println!("Files with $RXN blocks: {}", rxn_count);
    println!("Files missing V2000: {}", missing_v2000);
    println!("Empty/tiny files: {}", empty_files);

    // Group by specific issues
    let mut issue_groups: std::collections::HashMap<String, Vec<&FileAnalysis>> =
        std::collections::HashMap::new();
    for analysis in &analyses {
        issue_groups
            .entry(analysis.specific_issue.clone())
            .or_default()
            .push(analysis);
    }

    println!("\n=== ISSUE BREAKDOWN ===");
    for (issue, files) in issue_groups {
        println!("{}: {} files", issue, files.len());
        for file in files.iter().take(3) {
            println!("  - {}", file.path.display());
        }
        if files.len() > 3 {
            println!("  ... and {} more", files.len() - 3);
        }
        println!();
    }
}
