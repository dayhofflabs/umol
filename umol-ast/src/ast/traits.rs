//! AST marker and conversion traits.
//!
//! `Ast` pins each AST type to its lowering/raising `Config`. `FromAst` and
//! `ToAst` are the boundary between the AST and external forms (concrete
//! graph-backed types, DSL surface types). Each impl picks its own `Error`.

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::config::{
    AromaticSystemAstConfig, AtomAstConfig, BondAstConfig, MulticenterBondAstConfig,
};
use super::multicenter::MulticenterBondAst;

pub trait Ast {
    type Config;
}

pub trait FromAst<A: Ast>: Sized {
    type Error;
    fn from_ast(ast: &A, cfg: &A::Config) -> Result<Self, Self::Error>;
}

pub trait ToAst<A: Ast> {
    type Error;
    fn to_ast(&self, cfg: &A::Config) -> Result<A, Self::Error>;
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}

impl Ast for BondAst {
    type Config = BondAstConfig;
}

impl Ast for AromaticSystemAst {
    type Config = AromaticSystemAstConfig;
}

impl Ast for MulticenterBondAst {
    type Config = MulticenterBondAstConfig;
}
