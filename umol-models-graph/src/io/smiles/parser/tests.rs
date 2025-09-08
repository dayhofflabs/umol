use pretty_assertions::assert_eq;

use super::*;
use crate::io::ir::{Molecule, SourceFormat};
use crate::io::smiles::lexer::Lexer;

#[test]
fn test_parser() {
    let lexer = Lexer::new("C");
    let mut errors = Vec::new();
    let parser = grammar::MoleculeParser::new();
    let result = parser.parse(&mut errors, lexer).unwrap();
    assert_eq!(result, Molecule::new(SourceFormat::SMILES));
}
