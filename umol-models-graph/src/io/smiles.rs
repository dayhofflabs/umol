//! SMILES format parser and writer.

// Removed legacy logos-based lexer modules
// pub mod lexer;
// pub mod lexer_old;
pub mod linter;
pub mod parser;
// pub mod parser_old;
pub mod iterators;
pub mod state;
// Removed legacy FSM stages M0–M5 and old M6 alias; consolidated into parser
// pub mod fsm_m0;
// pub mod fsm_m1;
// pub mod fsm_m2;
// pub mod fsm_m3;
// pub mod fsm_m4;
// pub mod fsm_m5;
// pub mod fsm_m6;
#[cfg(test)]
pub mod test_support;

pub use parser::parse_smiles;
pub use parser::M6Error as ParseError;
