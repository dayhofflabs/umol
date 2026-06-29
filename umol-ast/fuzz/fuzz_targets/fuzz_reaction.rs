#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::ReactionDsl;
use umol_edn::{read_string, FromEdn};

fuzz_target!(|data: &str| {
    // Streaming path: single-pass over the bytes.
    let stream = ReactionDsl::from_edn_str(data).ok();

    // Tree path: parse the EDN tree first, then lift to the DSL.
    let tree = read_string(data)
        .ok()
        .and_then(|edn| ReactionDsl::from_edn(&edn).ok());

    // Parse-or-error parity: neither path may panic, and when both produce a
    // value they must agree.
    if let (Some(s), Some(t)) = (&stream, &tree) {
        assert_eq!(s, t, "streaming and tree reaction parsers disagree");
    }
});
