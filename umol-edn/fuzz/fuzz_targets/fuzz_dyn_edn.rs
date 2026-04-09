#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_edn::serde::{from_str_with, to_string, DynEdn};
use umol_edn::ParseConfig;

fuzz_target!(|data: &str| {
    let mut cfg = ParseConfig::default();
    cfg.allow_unknown_tags = true;

    if let Ok(val) = from_str_with::<DynEdn>(data, &cfg) {
        // Serialize must not panic.
        if let Ok(rendered) = to_string(&val) {
            // Re-deserialize must not panic.
            let _ = from_str_with::<DynEdn>(&rendered, &cfg);
        }
    }
});
