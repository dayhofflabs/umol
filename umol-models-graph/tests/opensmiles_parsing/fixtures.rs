//! Fixtures for OpenSMILES (UMOL) parsing tests
//! 
use rstest::fixture;
use umol_models_graph::io::smiles::lexer::Lexer;
use umol_models_graph::io::smiles::parser::grammar::MoleculeParser;
use umol_models_graph::io::smiles::state::ParseState;
use umol_models_graph::io::ir::Molecule as IRMolecule;

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

pub fn parse_and_assert_invariants(input: &str) -> Result<(), String> {
    let mut state = ParseState::default();
    let lexer = Lexer::new(input);
    let parser = MoleculeParser::new();
    parser
        .parse(&mut state, lexer)
        .map_err(|_| format!("parse failed: {}", input))?;

    if !state.rings.is_empty() {
        return Err(format!("open rings after parse: {}", input));
    }

    let mols: Vec<IRMolecule> = state.drain_molecules();
    for (mi, m) in mols.iter().enumerate() {
        if m.bonds.len() < m.atoms.len().saturating_sub(1) {
            return Err(format!(
                "edges < vertices-1 in molecule {}: {} bonds, {} atoms; {}",
                mi,
                m.bonds.len(),
                m.atoms.len(),
                input
            ));
        }
        for (bi, b) in m.bonds.iter().enumerate() {
            if let (Some(sa), Some(ea)) = (b.start_atom, b.end_atom) {
                if sa == ea {
                    return Err(format!("self-loop bond {} in molecule {}: {}", bi, mi, input));
                }
            }
        }
    }
    Ok(())
}
