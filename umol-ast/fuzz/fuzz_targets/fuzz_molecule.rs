#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::MoleculeDsl;
use umol_edn::{read_string, FromEdn};

fuzz_target!(|data: &str| {
    // Streaming path: single-pass over the bytes.
    let _ = MoleculeDsl::from_edn_str(data);

    // Tree path: parse EDN tree first, then lift to DSL.
    if let Ok(edn) = read_string(data) {
        let _ = MoleculeDsl::from_edn(&edn);
    }
});
