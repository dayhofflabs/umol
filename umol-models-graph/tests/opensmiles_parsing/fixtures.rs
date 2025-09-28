//! Fixtures for OpenSMILES (UMOL) parsing tests
//! 
use rstest::fixture;
use umol_models_graph::io::smiles::lexer_old::Lexer;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::ParseState;
// removed IR import; parser returns unit in current setup

#[fixture]
pub fn seed() -> u64 {
    20250913
}

#[fixture]
pub fn rng(seed: u64) -> fastrand::Rng {
    let mut r = fastrand::Rng::new();
    r.seed(seed);
    r
}

#[fixture]
pub fn parser() -> MoleculeParser {
    MoleculeParser::new()
}

#[fixture]
pub fn parse_state() -> ParseState {
    ParseState::default()
}

#[allow(dead_code)]
pub fn count_tokens(input: &str) -> usize {
    let lexer = Lexer::new(input);
    lexer.map(|t| t.ok()).count()
}

pub fn parse_and_assert_invariants(input: &str) -> Result<(), String> {
    let mut state = ParseState::default();
    let parser = MoleculeParser::new();
    let lexer = Lexer::new(input);
    parser.parse(&mut state, lexer).map_err(|e| format!("{:?}", e))
}
