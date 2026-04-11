//! Bond structural AST.

use crate::ast::config::BondAstConfig;
use crate::ast::value::ValueAst;
use crate::ast::Ast;

/// Bond AST: structural representation of a bond (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: Option<ValueAst>,
    pub unpaired_electrons: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self {
            order: ValueAst::Lit(order as i64),
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
        }
    }
}

impl Ast for BondAst {
    type Config = BondAstConfig;
}
