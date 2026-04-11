use std::process::exit;

use clap::Parser;
use umol_graph::graph_ir::molecule_builder::MoleculeBuilder;
use umol_graph::graph_ir::resolve_molecule;
use umol_graph::io::smiles::parse_smiles;

#[derive(Parser)]
#[command(name = "parse-smiles")]
#[command(about = "Parse a SMILES string and show per-atom resolution results")]
struct Args {
    smiles: String,
}

fn main() {
    let args = Args::parse();
    let table_mol = match parse_smiles(&args.smiles) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("parse error: {}", e);
            exit(1);
        }
    };

    println!("table atoms: {}", table_mol.atoms.len());
    for (i, atom) in table_mol.atoms.iter().enumerate() {
        println!(
            "  [{}] {:?} charge={:?} H={:?} aromatic={:?}",
            i, atom.element, atom.charge, atom.implicit_hydrogens, atom.aromatic
        );
    }
    println!("table bonds: {}", table_mol.bonds.len());
    for (i, bond) in table_mol.bonds.iter().enumerate() {
        println!(
            "  [{}] {}-{} order={:?}",
            i,
            bond.atoms.first(),
            bond.atoms.second(),
            bond.order
        );
    }

    let mut builder = MoleculeBuilder::from_table_molecule(&table_mol);
    let result = resolve_molecule(&mut builder);
    match result {
        Ok(()) => {
            println!("\nresolved atoms:");
            for ai in builder.atom_indices() {
                let atom = builder.atom(ai).unwrap();
                println!("  [{}] {}", ai.index(), atom);
            }
        }
        Err(e) => {
            println!("resolution error: {}", e);
        }
    }
}
