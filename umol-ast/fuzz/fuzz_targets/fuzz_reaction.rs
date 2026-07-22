#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::ast::{IntoAst, SubstructureMatchAlgorithm};
use umol_ast::dsl::{ReactionDefaults, ReactionDsl};
use umol_edn::{read_string, FromEdn};
use umol_graph_core::SubgraphIsomorphismAlgorithm;

fuzz_target!(|data: &str| {
    // Streaming path: single-pass over the bytes.
    let stream = ReactionDsl::from_edn_str(data).ok();

    // Tree path: parse the EDN tree first, then lift to the DSL.
    let tree = read_string(data)
        .ok()
        .and_then(|edn| ReactionDsl::from_edn(&edn).ok());

    // Parse-or-error parity: neither path may panic, and both must either
    // reject or produce the same value.
    assert_eq!(stream, tree, "streaming and tree reaction parsers disagree");

    if let Some(dsl) = stream {
        let reaction = dsl.into_ast(&ReactionDefaults::default());
        let _ = reaction.validate_application(&reaction.lhs);
        if let Ok(applications) = reaction.apply(
            &reaction.lhs,
            SubstructureMatchAlgorithm::GraphAndOverlays,
            SubgraphIsomorphismAlgorithm::Vf2,
        ) {
            for application in applications.take(16) {
                let _ = application;
            }
        };
    }
});
