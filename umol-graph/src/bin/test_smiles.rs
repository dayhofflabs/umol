use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;

use clap::{Parser, ValueEnum};
use umol_io::io::smiles::config::{SmilesIoConfig, SmilesParseFlags};
use umol_io::io::smiles::parse_extended_smiles_bytes_with;
use umol_io::io::smiles::parser::parse_smiles_bytes_to_table_ir_with;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ParserType {
    /// Basic OpenSMILES (strict, no extensions)
    BasicOpensmiles,
    /// OpenSMILES with wildcards (extended parser)
    Opensmiles,
    /// Basic parser with default flags
    Basic,
    /// Basic lenient parser
    BasicLenient,
    /// Extended parser with default flags
    Extended,
    /// Extended lenient parser
    Lenient,
    /// Basic parser with CXSMILES support
    BasicChemaxon,
    /// Extended parser with CXSMILES support
    Chemaxon,
}

#[derive(Parser)]
#[command(name = "test_smiles")]
#[command(about = "Test SMILES parsing on a file or stdin")]
struct Args {
    /// Parser type to use
    #[arg(short, long, default_value = "basic-opensmiles")]
    parser: ParserType,

    /// Path to .smi file (use "-" or omit for stdin)
    #[arg(default_value = "-")]
    path: String,
}

/// Extract CXSMILES portion from a line (SMILES + |...|).
fn extract_cxsmiles(line: &str) -> Option<&str> {
    if let Some(pipe_start) = line.find(" |") {
        let extension_start = pipe_start + 2;
        if let Some(rel_end) = line[extension_start..].find('|') {
            let end = extension_start + rel_end + 1;
            return Some(&line[..end]);
        }
    }
    None
}

fn main() {
    let args = Args::parse();

    let input: Box<dyn BufRead> = if args.path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let file = File::open(&args.path).unwrap_or_else(|e| {
            eprintln!("error: cannot open {}: {}", args.path, e);
            process::exit(2)
        });
        Box::new(BufReader::new(file))
    };

    let (use_extended, config, include_chemaxon) = match args.parser {
        ParserType::BasicOpensmiles => (false, SmilesIoConfig::basic_opensmiles(), false),
        ParserType::Opensmiles => (true, SmilesIoConfig::opensmiles(), false),
        ParserType::Basic => (false, SmilesIoConfig::basic(), false),
        ParserType::BasicLenient => (
            false,
            SmilesIoConfig::with_parse_flags(
                SmilesParseFlags::BASIC_MAX & SmilesParseFlags::LENIENT,
            ),
            false,
        ),
        ParserType::Extended => (true, SmilesIoConfig::extended(), false),
        ParserType::Lenient => (true, SmilesIoConfig::lenient(), false),
        ParserType::BasicChemaxon => (false, SmilesIoConfig::basic_chemaxon(), true),
        ParserType::Chemaxon => (true, SmilesIoConfig::chemaxon(), true),
    };

    let mut n = 0usize;
    let mut ok = 0usize;
    for line in input.lines() {
        let line = match line {
            Ok(s) => s,
            Err(_) => continue,
        };
        let smiles = if include_chemaxon {
            extract_cxsmiles(&line).unwrap_or_else(|| line.split_whitespace().next().unwrap_or(""))
        } else {
            line.split_whitespace().next().unwrap_or("")
        };
        if smiles.is_empty() {
            continue;
        }
        n += 1;
        let success = if use_extended {
            parse_extended_smiles_bytes_with(smiles.as_bytes(), &config).is_ok()
        } else {
            parse_smiles_bytes_to_table_ir_with(smiles.as_bytes(), &config).is_ok()
        };
        if success {
            ok += 1;
        }
    }
    eprintln!("parsed: {} lines (ok: {})", n, ok);
}
