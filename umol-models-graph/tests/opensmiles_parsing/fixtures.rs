use rstest::fixture;
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::ParseState;

#[fixture]
pub fn seed() -> u64 { 20250913 }

#[fixture]
pub fn rng(seed: u64) -> fastrand::Rng {
    let mut r = fastrand::Rng::new();
    r.seed(seed);
    r
}

#[fixture]
pub fn parser() -> MoleculeParser { MoleculeParser::new() }

#[fixture]
pub fn parse_state() -> ParseState { ParseState::default() }

pub fn accepts(input: &str) -> bool {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    parser.parse(&mut state, lexer).is_ok()
}
