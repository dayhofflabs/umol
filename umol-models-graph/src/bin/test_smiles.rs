use std::fs::File;
use std::io::{self, BufRead, BufReader};

use clap::{Parser, ValueEnum};
use umol_models_graph::io::smiles::config::{SmilesIoConfig, SmilesParseFlags};
use umol_models_graph::io::smiles::{parse_extended_smiles_bytes_with, parse_smiles_bytes_with};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ParserType {
    /// Basic OpenSMILES (strict, no extensions)
    BasicOpensmiles,
    /// OpenSMILES with wildcards (extended parser)
    Opensmiles,
    /// Basic parser with default flags
    Basic,
    /// Basic parser with BASIC_MAX flags
    BasicLenient,
    /// Extended parser with default flags
    Extended,
    /// Extended parser with LENIENT flags
    Lenient,
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

fn main() {
    let args = Args::parse();

    let input: Box<dyn BufRead> = if args.path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let file = File::open(&args.path).unwrap_or_else(|e| {
            eprintln!("error: cannot open {}: {}", args.path, e);
            std::process::exit(2)
        });
        Box::new(BufReader::new(file))
    };

    let (use_extended, config) = match args.parser {
        ParserType::BasicOpensmiles => (false, SmilesIoConfig::basic_opensmiles()),
        ParserType::Opensmiles => (true, SmilesIoConfig::opensmiles()),
        ParserType::Basic => (false, SmilesIoConfig::basic()),
        ParserType::BasicLenient => {
            (false, SmilesIoConfig::with_parse_flags(SmilesParseFlags::BASIC_MAX))
        }
        ParserType::Extended => (true, SmilesIoConfig::extended()),
        ParserType::Lenient => (true, SmilesIoConfig::lenient()),
    };

    let mut n = 0usize;
    let mut ok = 0usize;
    for line in input.lines() {
        let line = match line {
            Ok(s) => s,
            Err(_) => continue,
        };
        // Allow SMILES [tab or space] name
        let smiles = line.split_whitespace().next().unwrap_or("");
        if smiles.is_empty() {
            continue;
        }
        n += 1;
        let success = if use_extended {
            parse_extended_smiles_bytes_with(smiles.as_bytes(), &config).is_ok()
        } else {
            parse_smiles_bytes_with(smiles.as_bytes(), &config).is_ok()
        };
        if success {
            ok += 1;
        }
    }
    eprintln!("parsed: {} lines (ok: {})", n, ok);
}
