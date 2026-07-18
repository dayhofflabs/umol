#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_io::smiles::Smiles;

fuzz_target!(|data: &[u8]| {
    let _ = std::panic::catch_unwind(|| {
        let _ = Smiles::parse_bytes(data);
    });
});
