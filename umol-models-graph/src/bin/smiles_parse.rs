use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

use umol_models_graph::io::smiles::parse_smiles;

fn main() {
    let mut args = env::args().skip(1);
    let path = args.next().unwrap_or_else(|| {
        eprintln!("usage: smiles_parse <path-to-.smi>");
        std::process::exit(2)
    });
    let file = File::open(&path).unwrap_or_else(|e| {
        eprintln!("error: cannot open {}: {}", path, e);
        std::process::exit(2)
    });
    let reader = BufReader::new(file);

    let mut n = 0usize;
    let mut ok = 0usize;
    for line in reader.lines() {
        let line = match line {
            Ok(s) => s,
            Err(_) => continue,
        };
        // Allow SMILES [tab or space] name
        let smiles = line.split_whitespace().next().unwrap_or("");
        if smiles.is_empty() { continue; }
        n += 1;
        let _ = parse_smiles(smiles.as_bytes()).map(|_| { ok += 1; });
    }
    eprintln!("parsed: {} lines (ok: {})", n, ok);
}


