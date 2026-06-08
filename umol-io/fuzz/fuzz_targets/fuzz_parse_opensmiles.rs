#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_io::smiles::parser::parse_smiles_bytes_to_table_ir;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        let _ = parse_smiles_bytes_to_table_ir(data);
    });
});
