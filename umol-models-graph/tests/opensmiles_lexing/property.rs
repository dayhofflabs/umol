//! Property-based lexing tests for OpenSMILES (UMOL)

use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;
use umol_models_graph::io::smiles::lexer::{Lexer, Token};

fn token_piece_strategy() -> impl Strategy<Value = String> {
    // Curated token-like pieces; include valid and some invalid forms
    let fixed = prop_oneof![
        // Atoms (single and multi-char)
        Just("C".to_string()),
        Just("N".to_string()),
        Just("c".to_string()),
        Just("*".to_string()),
        Just("Cl".to_string()),
        Just("Br".to_string()),
        Just("se".to_string()),
        Just("as".to_string()),
        // Brackets
        Just("[".to_string()),
        Just("]".to_string()),
        // Bonds and punctuation
        Just("-".to_string()),
        Just("=".to_string()),
        Just(":".to_string()),
        Just("#".to_string()),
        Just("$".to_string()),
        Just("/".to_string()),
        Just("\\".to_string()),
        Just(".".to_string()),
        // Percent forms
        Just("%12".to_string()),
        Just("%01".to_string()),
        Just("%0".to_string()),
        // Chirality heads
        Just("@".to_string()),
        Just("@@".to_string()),
        Just("@TH".to_string()),
        Just("@AL".to_string()),
        Just("@SP".to_string()),
        Just("@TB".to_string()),
        Just("@OH".to_string()),
        // Charges
        Just("+".to_string()),
        Just("-".to_string()),
        Just("++".to_string()),
        Just("--".to_string()),
        // Digits (0-9)
        Just("0".to_string()),
        Just("1".to_string()),
        Just("2".to_string()),
        Just("3".to_string()),
        Just("4".to_string()),
        Just("5".to_string()),
        Just("6".to_string()),
        Just("7".to_string()),
        Just("8".to_string()),
        Just("9".to_string()),
        // Whitespace
        Just(" ".to_string()),
        Just("\t".to_string()),
        Just("\n".to_string()),
    ];
    fixed
}

fn input_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(token_piece_strategy(), 0..32).prop_map(|pieces| pieces.concat())
}

fn collect_spanned(input: &str) -> Vec<(usize, Token, usize)> {
    Lexer::new(input.as_bytes())
        .map(|t| t.ok())
        .collect::<Option<Vec<_>>>()
        .expect("Lexer wrapper yields Ok with Token::Error for bad chars")
}

fn collect_tokens(input: &str) -> Vec<Token> {
    collect_spanned(input)
        .into_iter()
        .map(|(_, t, _)| t)
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 128, max_shrink_time: 20_000,
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct("umol-models-graph/tests/proptest-regressions"))),
        .. ProptestConfig::default() })]

    #[test]
    fn spans_cover_input_and_do_not_overlap(input in input_strategy()) {
        let spanned = collect_spanned(&input);
        let mut prev_end = 0usize;
        for (start, _tok, end) in &spanned {
            prop_assert!(*start == prev_end, "non-contiguous spans: prev_end={}, start={}", prev_end, start);
            prev_end = *end;
        }
        prop_assert!(prev_end == input.len(), "final end {} != input len {}", prev_end, input.len());
    }

    #[test]
    fn idempotent_relex_of_span_slices(input in input_strategy()) {
        let spanned1 = collect_spanned(&input);
        let tokens1: Vec<Token> = spanned1.iter().map(|(_, t, _)| t.clone()).collect();
        // Reconstruct by concatenating slices defined by spans
        let mut rebuilt = String::with_capacity(input.len());
        for (start, _t, end) in &spanned1 {
            rebuilt.push_str(&input[*start..*end]);
        }
        prop_assert_eq!(rebuilt.as_str(), input.as_str());
        let tokens2 = collect_tokens(&rebuilt);
        prop_assert_eq!(tokens1, tokens2);
    }
}
