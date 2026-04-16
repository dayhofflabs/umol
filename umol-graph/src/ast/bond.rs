//! Bond structural AST.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::table_ir::bond::{Bond as TableBond, BondOrder};

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

    /// Lift a `table_ir::Bond` to a `BondAst` field-by-field. `BondOrder::Aromatic`
    /// maps to order 1 (the aromatic hint is carried on atoms, not bonds).
    /// Query orders (`SingleOrDouble`, `Any`, ...) and absent fields map to
    /// `Undetermined`; callers that need a ground AST apply defaults explicitly.
    pub fn from_table_bond(bond: &TableBond) -> Self {
        let order = match bond.order {
            BondOrder::Aromatic => ValueAst::Lit(1),
            o => o
                .value()
                .map_or(ValueAst::Undetermined, |v| ValueAst::Lit(v as i64)),
        };
        let charge = bond
            .charge
            .map_or(ValueAst::Undetermined, |c| ValueAst::Lit(c as i64));
        let u = bond
            .unpaired_electrons
            .map_or(ValueAst::Undetermined, |u| ValueAst::Lit(u as i64));
        let m = bond
            .multiplicity
            .map_or(ValueAst::Undetermined, |m| ValueAst::Lit(m.multiplicity() as i64));
        let mut bond_ast = Self {
            order,
            charge,
            spin: SpinStateAst::from_pair(u, m),
        };
        if let Ok(Some(state)) = bond_ast.spin.try_into_ground() {
            bond_ast.spin = SpinStateAst::Lit(state);
        }
        bond_ast
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
