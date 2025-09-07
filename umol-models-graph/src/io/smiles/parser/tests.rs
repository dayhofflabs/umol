use super::*;
use umol_data::Element;

#[test]
fn test_grammar() {
    println!("test_grammar");
    let result = grammar::MoleculeParser::new().parse("C");
    assert_eq!(result, Ok(Element::C));
}