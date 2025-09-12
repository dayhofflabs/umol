//! SMILES format parser.

use lalrpop_util::lalrpop_mod;

// SMILES unbranched grammar parser
lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub unbranched
);

// SMILES branched grammar parser
lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub branched
);

#[cfg(test)]
mod tests;