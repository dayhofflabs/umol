//! Bond AST data structures.

use std::mem;

use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

use super::config::{BondAstConfig, MultiplicityMode, NumericMode, UnpairedElectronsMode};
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
            spin: SpinStateAst::from_values(u, m),
        };
        if let Ok(Some(state)) = bond_ast.spin.try_into_ground() {
            bond_ast.spin = SpinStateAst::from_state(state);
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
        }) && (match target.spin.try_into_ground() {
            Ok(Some(s)) => self.spin.matches(s),
            _ => matches!(self.spin.unpaired, ValueAst::Undetermined)
                && matches!(self.spin.multiplicity, ValueAst::Undetermined),
        })
    }

    pub fn is_ground(&self) -> bool {
        self.order.is_ground()
            && self.charge.is_ground()
            && self.spin.is_ground()
    }

    pub fn coerce(&mut self, cfg: &BondAstConfig) {
        if matches!(self.charge, ValueAst::Undetermined) {
            self.charge = match cfg.charge_mode {
                NumericMode::Zero => ValueAst::Lit(0),
                NumericMode::Required => ValueAst::Undetermined,
            };
        }
        coerce_spin(&mut self.spin, cfg);
    }

    pub fn release(&mut self, cfg: &BondAstConfig) {
        if matches!(
            (&cfg.charge_mode, &self.charge),
            (NumericMode::Zero, ValueAst::Lit(0))
        ) {
            self.charge = ValueAst::Undetermined;
        }
        release_spin(&mut self.spin, cfg);
    }
}

pub(super) fn coerce_spin(spin: &mut SpinStateAst, cfg: &BondAstConfig) {
    match mem::take(spin) {
        SpinStateAst {
            unpaired,
            multiplicity,
        } => {
            let resolved_u = if matches!(unpaired, ValueAst::Undetermined) {
                match cfg.unpaired_electrons_mode {
                    UnpairedElectronsMode::Zero => ValueAst::Lit(0),
                    UnpairedElectronsMode::Required => ValueAst::Undetermined,
                    UnpairedElectronsMode::Derived => match &multiplicity {
                        ValueAst::Lit(m) => ValueAst::Lit(m - 1),
                        _ => ValueAst::Undetermined,
                    },
                }
            } else {
                unpaired
            };
            let resolved_m = if matches!(multiplicity, ValueAst::Undetermined) {
                match cfg.multiplicity_mode {
                    MultiplicityMode::Required => ValueAst::Undetermined,
                    MultiplicityMode::Derived => match &resolved_u {
                        ValueAst::Lit(u) => ValueAst::Lit(u + 1),
                        _ => ValueAst::Undetermined,
                    },
                }
            } else {
                multiplicity
            };
            *spin = SpinStateAst::from_values(resolved_u, resolved_m);
            if let Ok(Some(state)) = spin.try_into_ground() {
                *spin = SpinStateAst::from_state(state);
            }
        }
        lit => *spin = lit,
    }
}

pub(super) fn release_spin(spin: &mut SpinStateAst, cfg: &BondAstConfig) {
    match mem::take(spin) {
        SpinStateAst::from_state(state) => {
            let u_value = state.unpaired_electrons();
            let m_value = state.multiplicity();
            let derived = m_value.multiplicity() == u_value + 1;
            let u_ast = match cfg.unpaired_electrons_mode {
                UnpairedElectronsMode::Zero if u_value == 0 => ValueAst::Undetermined,
                UnpairedElectronsMode::Derived if derived => match cfg.multiplicity_mode {
                    MultiplicityMode::Required => ValueAst::Undetermined,
                    MultiplicityMode::Derived if u_value == 0 => ValueAst::Undetermined,
                    MultiplicityMode::Derived => ValueAst::Lit(u_value as i64),
                },
                _ => ValueAst::Lit(u_value as i64),
            };
            let m_ast = match cfg.multiplicity_mode {
                MultiplicityMode::Derived if derived => ValueAst::Undetermined,
                _ => ValueAst::Lit(m_value.multiplicity() as i64),
            };
            *spin = SpinStateAst::from_values(u_ast, m_ast);
        }
        pair => *spin = pair,
    }
}
