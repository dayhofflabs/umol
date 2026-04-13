//! Bond structural AST.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

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

    pub fn is_ground(&self) -> bool {
        self.order.is_ground()
            && self.charge.as_ref().map_or(true, |v| v.is_ground())
            && self.spin.as_ref().map_or(true, |v| v.is_ground())
    }
}

impl Ast for BondAst {
    type Config = BondAstConfig;
}
