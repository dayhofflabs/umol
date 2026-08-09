#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_graph_ir::ir::{Molecule, TryIntoIr};
use umol_io::smiles::Smiles;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        if let Ok(smiles) = Smiles::parse_bytes(data) {
            let _: Result<Molecule, _> = smiles.as_table_ir().try_into_ir(&());
        }
    });
});
