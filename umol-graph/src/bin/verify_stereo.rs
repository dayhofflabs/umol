//! Parse a SMILES (or MOL) to a molecule AST and lower it to DSL — parse + lower
//! only (no resolve, no perception), so the output is exactly the inline `#T` /
//! `#C` stereo cosets Phase B produces. Used to generate the stereo resolution
//! conformance inputs (`umol-graph/tests/resolution/data/stereo_*`): pick a
//! SMILES, draw + verify it by name in ChemDraw, then commit the lowered DSL.
//!
//! Usage:
//!   verify_stereo smiles "<SMILES>"   (or `-` to read the SMILES from stdin)
//!   verify_stereo mol <file.mol>      (or `-` to read the MOL from stdin)
//!   verify_stereo table               (lower every entry in `EXAMPLES`, print DSL)
//!   verify_stereo write <data-dir>    (write `<data-dir>/<category>/<name>.edn`
//!                                      = `{:input <DSL>}` for every `EXAMPLES` entry)
//!
//! NOTE: the `EXAMPLES` SMILES are best-effort; the R/S and Z/E assignments must
//! be verified by name against a ChemDraw structure before the DSL is committed.

use std::io::{stdin, Read};
use std::path::Path;
use std::{env, fs, process};

use umol_ast::ast::{FromAst, MoleculeAst};
use umol_ast::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_io::ctfile::parse_mol_bytes_to_ast;
use umol_io::smiles::parse_smiles_to_ast;

/// `(category, name, smiles)` — the stereo conformance corpus. `category` is the
/// `data/<category>/` directory; `name` is the `<name>.edn` stem.
const EXAMPLES: &[(&str, &str, &str)] = &[
    // --- Tetrahedral (#T) -------------------------------------------------
    ("stereo_tetrahedral", "r-cfclbri", "[C@@](F)(Cl)(Br)I"),
    ("stereo_tetrahedral", "s-cfclbri", "[C@](F)(Cl)(Br)I"),
    ("stereo_tetrahedral", "cfclbri-unknown", "FC(Cl)(Br)I"),
    ("stereo_tetrahedral", "l-alanine", "N[C@@H](C)C(=O)O"),
    ("stereo_tetrahedral", "r-butan-2-ol", "CC[C@@H](C)O"),
    (
        "stereo_tetrahedral",
        "methyl-ethyl-sulfoxide",
        "C[S@@+]([O-])CC",
    ),
    (
        "stereo_tetrahedral",
        "methyl-ethyl-propyl-isopropyl-phosphonium",
        "CC[P@@+](C)(CCC)C(C)C",
    ),
    (
        "stereo_tetrahedral",
        "2r3r-dichlorobutane",
        "C[C@@H](Cl)[C@H](Cl)C",
    ),
    (
        "stereo_tetrahedral",
        "2s3s-dichlorobutane",
        "C[C@H](Cl)[C@@H](Cl)C",
    ),
    (
        "stereo_tetrahedral",
        "meso-dichlorobutane",
        "C[C@H](Cl)[C@H](Cl)C",
    ),
    ("stereo_tetrahedral", "r-methyloxirane", "C[C@@H]1CO1"),
    (
        "stereo_tetrahedral",
        "cis-1-2-dichlorocyclohexane",
        "Cl[C@H]1CCCC[C@H]1Cl",
    ),
    (
        "stereo_tetrahedral",
        "trans-1-2-dichlorocyclohexane",
        "Cl[C@H]1CCCC[C@@H]1Cl",
    ),
    (
        "stereo_tetrahedral",
        "2-3-dichloropentane-partial",
        "C[C@H](Cl)C(Cl)CC",
    ),
    (
        "stereo_tetrahedral",
        "l-ascorbic-acid",
        "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    ),
    (
        "stereo_tetrahedral",
        "r-1-phenylethanol",
        "C[C@@H](O)c1ccccc1",
    ),
    (
        "stereo_tetrahedral",
        "methyl-ethyl-propyl-isopropyl-ammonium",
        "CC[N@@+](C)(CCC)C(C)C",
    ),
    (
        "stereo_tetrahedral",
        "2-3-4-trichloropentane",
        "C[C@H](Cl)[C@H](Cl)[C@@H](Cl)C",
    ),
    (
        "stereo_tetrahedral",
        "trans-decalin",
        "[H][C@]12CCCC[C@]1([H])CCCC2",
    ),
    (
        "stereo_tetrahedral",
        "cis-decalin",
        "[H][C@]12CCCC[C@@]1([H])CCCC2"
    ),
    (
        "stereo_tetrahedral",
        "alpha-d-glucopyranose",
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    ),
    (
        "stereo_tetrahedral",
        "2r3e-pent-3-en-2-ol",
        "C[C@@H](O)/C=C/C",
    ),
    // --- Cis/trans (#C) ---------------------------------------------------
    (
        "stereo_cis_trans",
        "z-2-3-difluorobut-2-ene",
        r"C/C(F)=C(F)\C",
    ),
    (
        "stereo_cis_trans",
        "e-2-3-difluorobut-2-ene",
        "C/C(F)=C(F)/C",
    ),
    (
        "stereo_cis_trans",
        "difluorobut-2-ene-unknown",
        "CC(F)=C(F)C",
    ),
    ("stereo_cis_trans", "z-but-2-ene", r"C/C=C\C"),
    ("stereo_cis_trans", "z-2-fluorobut-2-ene", r"C/C(F)=C/C"),
    ("stereo_cis_trans", "e-azomethane", "C/N=N/C"),
    ("stereo_cis_trans", "z-azomethane", r"C/N=N\C"),
    ("stereo_cis_trans", "z-butan-2-one-oxime", r"C/C(=N/O)CC"),
    ("stereo_cis_trans", "2e4e-hexa-2-4-diene", "C/C=C/C=C/C"),
    (
        "stereo_cis_trans",
        "2e-hexa-2-4-diene-partial",
        "C/C=C/C=CC",
    ),
    ("stereo_cis_trans", "cyclohexene", "C1=CCCCC1"),
    ("stereo_cis_trans", "z-cyclooctene", r"C1/C=C\CCCCC1"),
    ("stereo_cis_trans", "e-cyclooctene", "C1/C=C/CCCCC1"),
];

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

/// Lower an AST to DSL with **zeroed** defaults — fully explicit (nothing elided),
/// matching how the resolution harness lowers its output.
fn lower(ast: &MoleculeAst) -> String {
    MoleculeDsl::from_ast(ast, &MoleculeDefaults::zeroed()).to_string()
}

/// Parse the SMILES to an AST, then lower the AST to DSL.
fn parse_and_lower(smiles: &str) -> Result<String, String> {
    let ast = parse_smiles_to_ast(smiles.trim()).map_err(|e| e.to_string())?;
    Ok(lower(&ast))
}

fn print_table() {
    for (category, name, smiles) in EXAMPLES {
        println!(";; ===== {category} / {name} =====");
        println!(";; smiles: {smiles}");
        match parse_and_lower(smiles) {
            Ok(dsl) => println!("{dsl}"),
            Err(e) => println!(";; PARSE ERROR: {e}"),
        }
        println!();
    }
}

/// Config-overrides paired with every zeroed-lowered input: pin charge and aromatic
/// valence (left undetermined by `MoleculeDefaults::default()`, which the harness
/// raises input with) so resolution is well-posed. Implicit hydrogens are left
/// undetermined on purpose — the resolver fills them. Atoms with a non-zero charge
/// or aromatic flag carry those inline (≠ zeroed → rendered) and override-unaffected.
const CONFIG_OVERRIDES: &str = "{:atom {:aromatic-valence :not-aromatic :charge :zero}}";

/// Write `<base>/<category>/<name>.edn` = `{:config-overrides … :input <DSL>}` for every entry.
fn write_inputs(base: &str) {
    for (category, name, smiles) in EXAMPLES {
        let dsl = match parse_and_lower(smiles) {
            Ok(dsl) => dsl,
            Err(e) => {
                eprintln!("{category}/{name}: PARSE ERROR: {e}");
                continue;
            }
        };
        let dir = Path::new(base).join(category);
        fs::create_dir_all(&dir).expect("create category dir");
        let path = dir.join(format!("{name}.edn"));
        let edn = format!("{{:config-overrides {CONFIG_OVERRIDES}\n :input {dsl}}}\n");
        fs::write(&path, edn).expect("write edn");
        println!("wrote {}", path.display());
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("table") if args.len() == 2 => print_table(),
        Some("write") if args.len() == 3 => write_inputs(&args[2]),
        Some("smiles") if args.len() == 3 => {
            let s = String::from_utf8(read_input(&args[2])).expect("utf-8 SMILES");
            match parse_and_lower(&s) {
                Ok(dsl) => println!("{dsl}"),
                Err(e) => {
                    eprintln!("parse error: {e}");
                    process::exit(1);
                }
            }
        }
        Some("mol") if args.len() == 3 => {
            let ast = parse_mol_bytes_to_ast(&read_input(&args[2])).unwrap_or_else(|e| {
                eprintln!("parse error: {e}");
                process::exit(1);
            });
            println!("{}", lower(&ast));
        }
        _ => {
            eprintln!("usage: verify_stereo <smiles <input|file|-> | mol <input|file|-> | table>");
            process::exit(2);
        }
    }
}
