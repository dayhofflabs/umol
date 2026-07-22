#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::ReactionSpanDsl;
use umol_edn::{read_string, FromEdn};

fuzz_target!(|data: &str| {
    // Streaming path: single-pass over the bytes.
    let stream = ReactionSpanDsl::from_edn_str(data).ok();

    // Tree path: parse the EDN tree first, then lift to the DSL.
    let tree = read_string(data)
        .ok()
        .and_then(|edn| ReactionSpanDsl::from_edn(&edn).ok());

    // Parse-or-error parity: neither path may panic, and both must either
    // reject or produce the same value.
    assert_eq!(stream, tree, "streaming and tree reaction-span parsers disagree");
});
