use std::{env, fs, process};

use umol_models_graph::io::mol::parser::{parse_mol, parse_mol_moleculelike};

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <mol_file>", args[0]);
        eprintln!(
            "Example: {} tests/mol_parsing/data_raw/cdk/hisotopes.mol",
            args[0]
        );
        process::exit(1);
    }

    let file_path = &args[1];

    // Read the MOL file
    let mol_bytes = match fs::read(file_path) {
        Ok(bytes) => bytes,
        Err(e) => {
            eprintln!("Error reading file '{}': {}", file_path, e);
            process::exit(1);
        }
    };

    println!("Testing file: {}", file_path);
    println!("Size: {} bytes", mol_bytes.len());
    println!();

    // Test with basic parser (parse_mol)
    println!("=== BASIC PARSER (parse_mol) ===");
    match parse_mol(&mol_bytes) {
        Ok(molecule) => {
            println!("✅ SUCCESS");
            println!("   Atoms: {}", molecule.atom_count());
            println!("   Bonds: {}", molecule.bond_count());
        }
        Err(e) => {
            println!("❌ FAILED");
            println!("   Error: {:?}", e);
        }
    }

    println!();

    // Test with extended parser (parse_mol_moleculelike)
    println!("=== EXTENDED PARSER (parse_mol_moleculelike) ===");
    match parse_mol_moleculelike(&mol_bytes) {
        Ok(molecule_like) => {
            println!("✅ SUCCESS");
            println!("   Atoms: {}", molecule_like.atom_count());
            println!("   Bonds: {}", molecule_like.bond_count());
        }
        Err(e) => {
            println!("❌ FAILED");
            println!("   Error: {:?}", e);
        }
    }

    println!();

    // Summary
    let basic_ok = parse_mol(&mol_bytes).is_ok();
    let extended_ok = parse_mol_moleculelike(&mol_bytes).is_ok();

    println!("=== CLASSIFICATION ===");
    match (basic_ok, extended_ok) {
        (true, true) => println!("📁 MOLECULE (basic parser works)"),
        (false, true) => println!("📁 MOLECULELIKE (needs extended features)"),
        (true, false) => println!("🐛 BUG: basic succeeds but extended fails!"),
        (false, false) => println!("❌ INVALID (both parsers fail)"),
    }
}
