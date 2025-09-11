//! SMILES format parser.

use lalrpop_util::lalrpop_mod;

// SMILES grammar parser
lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub grammar
);

// SMILES chain-only grammar parser
lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub chain
);

// SMILES tree grammar parser
lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub tree
);

#[cfg(test)]
mod tests;