#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_models_graph::io::smiles::parse_smiles;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        let _ = parse_smiles(data);
    });
});
