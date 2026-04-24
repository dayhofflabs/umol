#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_ast::dsl::{
    parse_aromatic_system, parse_atom, parse_bond, parse_dative_bond, parse_multicenter_bond,
    parse_noncovalent_bond, parse_value,
};

fuzz_target!(|data: &str| {
    // Entity-string parsers: must not panic, integer overflow, or hang on
    // arbitrary bytes. All return Result<_, ParseError>.
    let _ = parse_atom(data);
    let _ = parse_bond(data);
    let _ = parse_aromatic_system(data);
    let _ = parse_dative_bond(data);
    let _ = parse_multicenter_bond(data);
    let _ = parse_noncovalent_bond(data);
    let _ = parse_value(data);
});
