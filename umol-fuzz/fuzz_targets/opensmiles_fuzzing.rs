#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_models_graph::io::smiles::parse_smiles_m3;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        let _ = parse_smiles_m3(data);
    });
});
