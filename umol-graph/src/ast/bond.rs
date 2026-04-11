//! Bond structural AST.

use umol_shared::{SpinStateAst, ValueAst};

use crate::ast::config::BondAstConfig;
use crate::ast::Ast;

/// Bond AST: structural representation of a bond (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: Option<ValueAst>,
    pub spin: Option<SpinStateAst>,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: None,
            spin: None,
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self {
            order: ValueAst::Lit(order as i64),
            charge: None,
            spin: None,
        }
    }
}

impl Ast for BondAst {
    type Config = BondAstConfig;
}
