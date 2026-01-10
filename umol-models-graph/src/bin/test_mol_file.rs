use std::{env, fs, process};

use umol_models_graph::io::ctfile::config::CtfileIoConfig;
use umol_models_graph::io::ctfile::{
    parse_extended_mol_bytes, parse_extended_mol_bytes_with, parse_mol_bytes, parse_mol_bytes_with,
};

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

    // Test with basic parser (parse_mol_bytes)
    println!("=== BASIC PARSER (parse_mol_bytes) ===");
    let basic_parsed = parse_mol_bytes(&mol_bytes);
    let basic_ok = basic_parsed.is_ok();
    match basic_parsed {
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

    // Test with basic lenient parser (parse_mol_bytes_with)
    println!("=== BASIC LENIENT PARSER (parse_mol_bytes_with) ===");
    let basic_lenient_parsed = parse_mol_bytes_with(&mol_bytes, &CtfileIoConfig::basic_lenient());
    let basic_lenient_ok = basic_lenient_parsed.is_ok();
    match basic_lenient_parsed {
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

    // Test with extended parser (parse_extended_mol_bytes)
    println!("=== EXTENDED PARSER (parse_extended_mol) ===");
    let extended_parsed = parse_extended_mol_bytes(&mol_bytes);
    let extended_ok = extended_parsed.is_ok();
    match extended_parsed {
        Ok(extended_molecule) => {
            println!("✅ SUCCESS");
            println!("   Atoms: {}", extended_molecule.atom_count());
            println!("   Bonds: {}", extended_molecule.bond_count());
        }
        Err(e) => {
            println!("❌ FAILED");
            println!("   Error: {:?}", e);
        }
    }

    println!();

    // Test with extended lenient parser (parse_extended_mol_bytes_with)
    println!("=== EXTENDED LENIENT PARSER (parse_extended_mol_bytes_with) ===");
    let extended_lenient_parsed =
        parse_extended_mol_bytes_with(&mol_bytes, &CtfileIoConfig::extended_lenient());
    let extended_lenient_ok = extended_lenient_parsed.is_ok();
    match extended_lenient_parsed {
        Ok(extended_molecule) => {
            println!("✅ SUCCESS");
            println!("   Atoms: {}", extended_molecule.atom_count());
            println!("   Bonds: {}", extended_molecule.bond_count());
        }
        Err(e) => {
            println!("❌ FAILED");
            println!("   Error: {:?}", e);
        }
    }

    println!();

    println!("=== CLASSIFICATION ===");
    match (basic_ok, basic_lenient_ok, extended_ok, extended_lenient_ok) {
        (true, true, true, true) => println!("📁 MOLECULE (basic parser works)"),
        (false, true, false, true) => println!("📁 NON-STANDARD BASIC MOLECULE (basic lenient parser works)"),
        (false, false, true, true) => println!("📁 EXTENDED MOLECULE (extended parser works)"),
        (false, false, false, true) => println!("📁 NON-STANDARD EXTENDED MOLECULE (extended lenient parser works)"),
        (false, false, false, false) => println!("❌ INVALID (all parsers fail)"),
       _ => println!("🐛 BUG: hierarchy violation (basic: {basic_ok}, basic lenient: {basic_lenient_ok}, extended: {extended_ok}, extended lenient: {extended_lenient_ok})"),
    }
}
