//! SMILES format parser.

use lalrpop_util::lalrpop_mod;

lalrpop_mod!(
    #[allow(unused_imports)]
    #[rustfmt::skip]
    pub grammar
);

pub mod utils;

#[cfg(test)]
mod tests;