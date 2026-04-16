//! Atom structural AST.

use std::mem;

use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::{
    AromaticValenceMode, AtomAstConfig, ImplicitHydrogenMode, IsotopeMode, MultiplicityMode,
    NumericMode, UnpairedElectronsMode,
};
use crate::ast::Ast;

/// Atom AST: structural representation of an atom (ground or pattern).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: HydrogenAst,
    pub lone_pairs: ValueAst,
    pub spin: SpinStateAst,
    pub valence: ValueAst,
    pub donated_pairs: ValueAst,
    pub accepted_pairs: ValueAst,
    pub aromatic_valence: AromaticValenceAst,
    pub multicenter_valence: ValueAst,
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
            && (match &target.valence {
                ValueAst::Lit(n) => self.valence.matches(*n),
                _ => matches!(self.valence, ValueAst::Undetermined),
            })
            && (match &target.donated_pairs {
                ValueAst::Lit(n) => self.donated_pairs.matches(*n),
                _ => matches!(self.donated_pairs, ValueAst::Undetermined),
            })
            && (match &target.accepted_pairs {
                ValueAst::Lit(n) => self.accepted_pairs.matches(*n),
                _ => matches!(self.accepted_pairs, ValueAst::Undetermined),
            })
            && self.aromatic_valence.matches(&target.aromatic_valence)
            && (match &target.multicenter_valence {
                ValueAst::Lit(n) => self.multicenter_valence.matches(*n),
                _ => matches!(self.multicenter_valence, ValueAst::Undetermined),
            })
    }

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.is_ground()
            && self.charge.is_ground()
            && self.implicit_hydrogens.is_ground()
            && self.lone_pairs.is_ground()
            && self.spin.is_ground()
            && self.valence.is_ground()
            && self.donated_pairs.is_ground()
            && self.accepted_pairs.is_ground()
            && self.aromatic_valence.is_ground()
            && self.multicenter_valence.is_ground()
    }

    pub fn charge_or_zero(&self) -> i8 {
        match &self.charge {
            ValueAst::Lit(n) => *n as i8,
            _ => 0,
        }
    }
}

impl AtomAst {
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
        if matches!(self.aromatic_valence, AromaticValenceAst::Undetermined) {
            self.aromatic_valence = match cfg.aromatic_valence_mode {
                AromaticValenceMode::NotAromatic => AromaticValenceAst::NotAromatic,
                AromaticValenceMode::Aromatic => AromaticValenceAst::Value(ValueAst::Undetermined),
                AromaticValenceMode::Required => AromaticValenceAst::Undetermined,
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
        coerce_numeric(&mut self.valence, &cfg.valence_mode);
        coerce_numeric(&mut self.donated_pairs, &cfg.donated_pairs_mode);
        coerce_numeric(&mut self.accepted_pairs, &cfg.accepted_pairs_mode);
        coerce_numeric(&mut self.multicenter_valence, &cfg.multicenter_valence_mode);
    }

    /// Collapse fields back to `Undetermined` where the current value is what
    /// `coerce` would have produced. Call after solving to restore roundtrip
    /// fidelity with the DSL.
    pub fn release(&mut self, cfg: &AtomAstConfig) {
        if matches!(
            (&cfg.isotope_mode, &self.isotope_mass),
            (IsotopeMode::Natural, IsotopeAst::Natural)
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
        match (&cfg.aromatic_valence_mode, &self.aromatic_valence) {
            (
                AromaticValenceMode::NotAromatic | AromaticValenceMode::Required,
                AromaticValenceAst::NotAromatic,
            ) => {
                self.aromatic_valence = AromaticValenceAst::Undetermined;
            }
            (AromaticValenceMode::Aromatic, AromaticValenceAst::Value(ValueAst::Undetermined)) => {
                self.aromatic_valence = AromaticValenceAst::Undetermined;
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
        release_numeric(&mut self.valence, &cfg.valence_mode);
        release_numeric(&mut self.donated_pairs, &cfg.donated_pairs_mode);
        release_numeric(&mut self.accepted_pairs, &cfg.accepted_pairs_mode);
        release_numeric(&mut self.multicenter_valence, &cfg.multicenter_valence_mode);
    }
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}
