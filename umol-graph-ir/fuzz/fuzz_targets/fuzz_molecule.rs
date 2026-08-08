#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::MoleculeDsl;
use umol_edn::{read_string, FromEdn};

fuzz_target!(|data: &str| {
    // Streaming path: single-pass over the bytes.
    let stream = MoleculeDsl::from_edn_str(data).ok();

    // Tree path: parse EDN tree first, then lift to DSL.
    let tree = read_string(data)
        .ok()
        .and_then(|edn| MoleculeDsl::from_edn(&edn).ok());

    assert_eq!(stream, tree, "streaming and tree molecule parsers disagree");
});
