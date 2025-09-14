//! Invalid tests for OpenSMILES (UMOL)

use rstest::rstest;
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::ParseState;

fn rejects(input: &str) -> bool {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    parser.parse(&mut state, lexer).is_err()
}

fn enumerated_invalid() -> Vec<&'static str> {
    vec![
        // rings
        "C1CCC",        // open ring
        "C11",          // self loop
        "C1C1",         // 2-member ring
        // skip conflicting ring dirs until generator adds rings reliably
        "%0", "%01",  // invalid percent indices
        // bonds
        "C-", "CC/", "C=", "C#", "C$",
        // whitespace mid-string
        "C C", "C \t C", "C\nC",
        // brackets errors
        "[CH3", "[CH3-+]", "[C@H", "[C-:0:1]",
        // branch ambiguity
        "(CO)N",
    ]
}

fn mutate_missing_bracket() -> Vec<String> {
    vec!["[CH3]", "[NH4+]", "[13C-]", "[C@H+]", "[C:0]"].into_iter()
        .map(|s| s.trim_end_matches(']')).map(|s| s.to_string()).collect()
}

fn mutate_trailing_bond() -> Vec<String> {
    vec!["C", "CC", "C=C", "c1ccccc1"].into_iter()
        .flat_map(|s| ["-", "=", "/", "\\", ":"].into_iter().map(move |b| format!("{}{}", s, b)))
        .collect()
}

fn mutate_whitespace() -> Vec<String> {
    vec!["CC", "CO", "C1CCCCC1"].into_iter()
        .flat_map(|_s| [" ", "\t", "\n"].into_iter().map(move |ws| format!("C{}C", ws)))
        .collect()
}

#[rstest]
fn invalid_enumerated_reject() {
    for s in enumerated_invalid() {
        assert!(rejects(s), "Expected rejection: {}", s);
    }
}

#[rstest]
fn invalid_mutations_reject() {
    for s in mutate_missing_bracket() { assert!(rejects(&s), "Expected rejection: {}", s); }
    for s in mutate_trailing_bond() { assert!(rejects(&s), "Expected rejection: {}", s); }
    for s in mutate_whitespace() { assert!(rejects(&s), "Expected rejection: {}", s); }
}
