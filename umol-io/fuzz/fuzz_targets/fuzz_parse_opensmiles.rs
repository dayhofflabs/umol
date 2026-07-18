#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_io::smiles::Smiles;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        if let Ok(smiles) = Smiles::parse_bytes(data) {
            let _: Result<MoleculeAst, _> = smiles.as_table_ir().try_into_ast(&());
        }
    });
});
