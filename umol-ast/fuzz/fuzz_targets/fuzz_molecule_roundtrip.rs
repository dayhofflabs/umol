#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::MoleculeDsl;
use umol_edn::{FromEdn, ToEdn};

fuzz_target!(|data: &str| {
    // If a string parses to a MoleculeDsl, rendering it back to EDN and
    // re-parsing must yield an equivalent MoleculeDsl.
    if let Ok(dsl) = MoleculeDsl::from_edn_str(data) {
        let rendered = dsl.to_edn().to_string();
        let reparsed = MoleculeDsl::from_edn_str(&rendered)
            .expect("render produced unparseable molecule DSL");
        assert_eq!(dsl, reparsed, "roundtrip mismatch:\n  in:  {data:?}\n  out: {rendered}");
    }
});
