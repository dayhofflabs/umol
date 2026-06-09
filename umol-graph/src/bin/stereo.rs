//! Dump the raised molecule — with per-atom `#T` / per-bond `#C` stereo cosets — for a
//! SMILES or MOL input. Parse + raise only (no resolve, no perception), so the output is
//! exactly what Phase B asserts. For blind validation against reference structures.
//!
//! Usage:
//!   stereo smiles "<SMILES>"     (or `-` to read the SMILES from stdin)
//!   stereo mol <file.mol>        (or `-` to read the MOL from stdin)

use std::env;
use std::fs;
use std::io::{stdin, Read};
use std::path::Path;
use std::process;

use umol_ast::dsl::{Metadata, MoleculeDsl};
use umol_io::ctfile::parse_mol_bytes_to_ast;
use umol_io::smiles::parse_smiles_to_ast;

fn read_input(arg: &str) -> Vec<u8> {
    if arg == "-" {
        let mut buf = Vec::new();
        stdin().read_to_end(&mut buf).expect("read stdin");
        buf
    } else if Path::new(arg).is_file() {
        fs::read(arg).expect("read file")
    } else {
        arg.as_bytes().to_vec()
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: stereo <smiles|mol> <input | file | ->");
        process::exit(2);
    }
    let raw = read_input(&args[2]);
    let ast = match args[1].as_str() {
        "smiles" => {
            let s = String::from_utf8(raw).expect("utf-8 SMILES");
            parse_smiles_to_ast(s.trim()).unwrap_or_else(|e| {
                eprintln!("parse error: {e}");
                process::exit(1);
            })
        }
        "mol" => parse_mol_bytes_to_ast(&raw).unwrap_or_else(|e| {
            eprintln!("parse error: {e}");
            process::exit(1);
        }),
        other => {
            eprintln!("unknown format `{other}` (expected `smiles` or `mol`)");
            process::exit(2);
        }
    };
    println!("{}", MoleculeDsl::from_parts(ast, Metadata::new()));
}
