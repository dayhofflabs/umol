//! Atom structural AST.

use std::mem;

use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::{
    AtomAstConfig, ImplicitHydrogenMode, IsotopeMode, MultiplicityMode, NumericMode,
    UnpairedElectronsMode,
};
use crate::table_ir::atom::{Atom as TableAtom, ImplicitHydrogens};

/// Atom AST: structural representation of an atom (ground or pattern).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: HydrogenAst,
    pub lone_pairs: ValueAst,
    pub spin: SpinStateAst,
}

impl AtomAst {
    pub fn new(element: ElementAst) -> Self {
        Self {
            element,
            ..Default::default()
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementAst::Lit(element))
    }

    /// Lift the base fields of a `table_ir::Atom` to an `AtomAst`. Topology- and
    /// chemistry-derived fields (`valence`, `aromatic`) are lifted at the
    /// `AtomPattern` level so they land in the molecule constraint vec.
    pub fn from_table_atom(atom: &TableAtom) -> Self {
        let mut ast = Self::from_element(atom.element);
        if let Some(mass) = atom.isotope_mass {
            ast.isotope_mass = IsotopeAst::Lit(mass);
        }
        if let Some(charge) = atom.charge {
            ast.charge = ValueAst::Lit(charge as i64);
        }
        match atom.implicit_hydrogens {
            Some(ImplicitHydrogens::Hydrogens(h)) => {
                ast.implicit_hydrogens = HydrogenAst::Value(ValueAst::Lit(h as i64));
            }
            Some(ImplicitHydrogens::Normal) => {
                ast.implicit_hydrogens = HydrogenAst::Normal;
            }
            None => {}
        }
        if let Some(lp) = atom.lone_pairs {
            ast.lone_pairs = ValueAst::Lit(lp as i64);
        }
        let u = atom
            .unpaired_electrons
            .map_or(ValueAst::Undetermined, |u| ValueAst::Lit(u as i64));
        let m = atom
            .multiplicity
            .map_or(ValueAst::Undetermined, |m| ValueAst::Lit(m.multiplicity() as i64));
        if !matches!(u, ValueAst::Undetermined) || !matches!(m, ValueAst::Undetermined) {
            ast.spin = SpinStateAst::from_pair(u, m);
        }
        ast
    }

    pub fn matches_ground(&self, target: &AtomAst) -> bool {
        (match &target.element {
            ElementAst::Lit(e) => self.element.matches(e),
            _ => false,
        }) && self.isotope_mass.matches(&target.isotope_mass)
            && (match &target.charge {
                ValueAst::Lit(n) => self.charge.matches(*n),
                _ => matches!(self.charge, ValueAst::Undetermined),
            })
            && self.implicit_hydrogens.matches(&target.implicit_hydrogens)
            && (match &target.lone_pairs {
                ValueAst::Lit(n) => self.lone_pairs.matches(*n),
                _ => matches!(self.lone_pairs, ValueAst::Undetermined),
            })
            && (match &target.spin {
                SpinStateAst::Lit(s) => self.spin.matches(*s),
                _ => matches!(self.spin, SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined }),
            })
    }

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.is_ground()
            && self.charge.is_ground()
            && self.implicit_hydrogens.is_ground()
            && self.lone_pairs.is_ground()
            && self.spin.is_ground()
    }

    pub fn charge_or_zero(&self) -> i8 {
        match &self.charge {
            ValueAst::Lit(n) => *n as i8,
            _ => 0,
        }
    }
}

impl AtomAst {
    // TODO: Verify Aromatic valence mode: Aromatic
    pub fn coerce(&mut self, cfg: &AtomAstConfig) {
        if matches!(self.isotope_mass, IsotopeAst::Undetermined) {
            self.isotope_mass = match cfg.isotope_mode {
                IsotopeMode::Natural => IsotopeAst::Natural,
                IsotopeMode::Required => IsotopeAst::Undetermined,
            };
        }
        if matches!(self.charge, ValueAst::Undetermined) {
            self.charge = match cfg.charge_mode {
                NumericMode::Zero => ValueAst::Lit(0),
                NumericMode::Required => ValueAst::Undetermined,
            };
        }
        if matches!(self.implicit_hydrogens, HydrogenAst::Undetermined) {
            self.implicit_hydrogens = match cfg.implicit_h_mode {
                ImplicitHydrogenMode::Normal => HydrogenAst::Normal,
                ImplicitHydrogenMode::Zero => HydrogenAst::Value(ValueAst::Lit(0)),
                ImplicitHydrogenMode::Required => HydrogenAst::Undetermined,
            };
        }
        match mem::take(&mut self.spin) {
            SpinStateAst::Pair {
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
                self.spin = SpinStateAst::from_pair(resolved_u, resolved_m);
                if let Ok(Some(state)) = self.spin.try_into_ground() {
                    self.spin = SpinStateAst::Lit(state);
                }
            }
            lit => self.spin = lit,
        }
        let coerce_numeric = |v: &mut ValueAst, mode: &NumericMode| {
            if matches!(*v, ValueAst::Undetermined) {
                *v = match mode {
                    NumericMode::Zero => ValueAst::Lit(0),
                    NumericMode::Required => ValueAst::Undetermined,
                };
            }
        };
        coerce_numeric(&mut self.lone_pairs, &cfg.lone_pairs_mode);
    }

    /// Collapse fields back to `Undetermined` where the current value is what
    /// `coerce` would have produced. Call after solving to restore roundtrip
    /// fidelity with the DSL.
    pub fn release(&mut self, cfg: &AtomAstConfig) {
        if matches!(
            (&cfg.isotope_mode, &self.isotope_mass),
            (IsotopeMode::Natural | IsotopeMode::Required, IsotopeAst::Natural)
        ) {
            self.isotope_mass = IsotopeAst::Undetermined;
        }
        if matches!(
            (&cfg.charge_mode, &self.charge),
            (NumericMode::Zero, ValueAst::Lit(0))
        ) {
            self.charge = ValueAst::Undetermined;
        }
        match (&cfg.implicit_h_mode, &self.implicit_hydrogens) {
            (ImplicitHydrogenMode::Normal, HydrogenAst::Normal) => {
                self.implicit_hydrogens = HydrogenAst::Undetermined;
            }
            (ImplicitHydrogenMode::Zero, HydrogenAst::Value(ValueAst::Lit(0))) => {
                self.implicit_hydrogens = HydrogenAst::Undetermined;
            }
            _ => {}
        }
        match mem::take(&mut self.spin) {
            SpinStateAst::Lit(state) => {
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
                self.spin = SpinStateAst::from_pair(u_ast, m_ast);
            }
            pair => self.spin = pair,
        }
        let release_numeric = |v: &mut ValueAst, mode: &NumericMode| {
            if matches!((mode, &*v), (NumericMode::Zero, ValueAst::Lit(0))) {
                *v = ValueAst::Undetermined;
            }
        };
        release_numeric(&mut self.lone_pairs, &cfg.lone_pairs_mode);
    }
}
