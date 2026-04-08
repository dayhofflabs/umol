#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_graph::io::smiles::parse_smiles_bytes;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        let _ = parse_smiles_bytes(data);
    });
});
