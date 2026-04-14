//! Bond structural AST.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::BondAstConfig;
use crate::ast::Ast;

/// Bond AST: structural representation of a bond (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: ValueAst::default(),
            spin: SpinStateAst::default(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(ValueAst::Lit(order as i64))
    }

    pub fn matches_ground(&self, target: &BondAst) -> bool {
        (match &target.order {
            ValueAst::Lit(n) => self.order.matches(*n),
            _ => false,
        }) && (match &target.charge {
            ValueAst::Lit(n) => self.charge.matches(*n),
            _ => matches!(self.charge, ValueAst::Undetermined),
        }) && (match &target.spin {
            SpinStateAst::Lit(s) => self.spin.matches(*s),
            _ => matches!(self.spin, SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined }),
        })
    }

    pub fn is_ground(&self) -> bool {
        self.order.is_ground()
            && self.charge.is_ground()
            && self.spin.is_ground()
    }
}

impl Ast for BondAst {
    type Config = BondAstConfig;
}
