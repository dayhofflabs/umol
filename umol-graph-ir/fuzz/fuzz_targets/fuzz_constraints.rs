#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_graph_ir::dsl::{ConstraintDsl, ConstraintsDsl};
use umol_edn::{read_string, FromEdn};

fuzz_target!(|data: &str| {
    // Streaming paths.
    let _ = ConstraintDsl::from_edn_str(data);
    let _ = ConstraintsDsl::from_edn_str(data);

    // Tree paths.
    if let Ok(edn) = read_string(data) {
        let _ = ConstraintDsl::from_edn(&edn);
        let _ = ConstraintsDsl::from_edn(&edn);
    }
});
