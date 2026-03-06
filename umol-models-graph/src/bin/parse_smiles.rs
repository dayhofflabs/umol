use clap::Parser;
use umol_models_graph::graph_ir::resolve_molecule;
use umol_models_graph::io::smiles::parse_smiles;

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
            std::process::exit(1);
        }
    };

    println!("table atoms: {}", table_mol.atoms.len());
    for (i, atom) in table_mol.atoms.iter().enumerate() {
        println!(
            "  [{}] {:?} charge={:?} H={:?} implicit_H={} aromatic={:?}",
            i, atom.element, atom.charge, atom.hydrogens, atom.implicit_hydrogens, atom.aromatic
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

    match resolve_molecule(&table_mol) {
        Ok(mol) => {
            println!("\nresolved atoms:");
            for ai in mol.atom_indices() {
                let atom = mol.atom(ai).unwrap();
                println!("  [{}] {}", ai.index(), atom.to_spec());
            }
        }
        Err(e) => {
            println!("\nresolution error: {}", e);
        }
    }
}
