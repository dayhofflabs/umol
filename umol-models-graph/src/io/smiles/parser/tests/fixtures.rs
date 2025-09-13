use rstest::*;
use umol::logging::setup_logger;
use umol::with_logger;
use slog::{o, Level};

use crate::io::smiles::state::ParseState;

#[fixture]
pub fn parse_state() -> ParseState {
    let mut state = ParseState::default();
    let root = setup_logger(Level::Debug);
    state.log = Some(with_logger!(root, "io::smiles::parser"));
    state
}


