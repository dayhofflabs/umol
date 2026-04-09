#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_edn::{read_string, read_string_with, DuplicateKeyPolicy, ParseConfig};

fuzz_target!(|data: &str| {
    // Default config: parse → display → re-parse must agree.
    if let Ok(val) = read_string(data) {
        let displayed = val.to_string();
        let reparsed = read_string(&displayed).expect("display produced unparseable EDN");
        assert_eq!(val, reparsed, "roundtrip mismatch: {displayed}");
    }

    // Permissive config: unknown tags allowed, duplicate keys last-wins.
    let permissive = ParseConfig {
        allow_unknown_tags: true,
        duplicate_keys: DuplicateKeyPolicy::LastWins,
        ..Default::default()
    };
    if let Ok(val) = read_string_with(data, &permissive) {
        let displayed = val.to_string();
        let reparsed =
            read_string_with(&displayed, &permissive).expect("display produced unparseable EDN");
        assert_eq!(val, reparsed, "roundtrip mismatch (permissive): {displayed}");
    }
});
