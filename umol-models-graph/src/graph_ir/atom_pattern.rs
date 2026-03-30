//! Atom pattern representation and conversion into concrete ground atoms.

use umol_data::{Element, SpinMultiplicity, SpinState};

use std::str::FromStr;

use crate::atom::{AromaticValence, IsotopeMass};
use crate::dsl::atom::{AtomAst, AtomLowerConfig, AromaticMode, ChargeMode, ImplicitHydrogenMode, IsotopeMode, parse_atom_dsl};
use crate::dsl::error::LoweringError;
use crate::dsl::lowering::FromAst;
use crate::dsl::predicates::{AromaticExpr, ElementExpr, HydrogenExpr, IsotopeExpr};
use crate::dsl::value::ValueAst;
use crate::graph_ir::atom::Atom;
use crate::graph_ir::atom_type::AtomError;
use crate::graph_ir::error::ValidationError;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};
use crate::table_ir::atom::ImplicitHydrogens;
use crate::table_ir::atom::Atom as TableAtom;

/// Generic pattern for a scalar-valued field: unconstrained or an exact value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Pattern<T> {
    Any,
    Is(T),
}

impl<T: Copy> Copy for Pattern<T> {}

impl<T: PartialEq + Copy> Pattern<T> {
    pub fn matches(&self, value: T) -> bool {
        match self {
            Self::Any => true,
            Self::Is(v) => *v == value,
        }
    }

    /// Convert to `Option<T>`: `Any` → `None`, `Is(v)` → `Some(v)`.
    pub fn into_option(self) -> Option<T> {
        match self {
            Self::Any => None,
            Self::Is(v) => Some(v),
        }
    }
}

/// Pattern on a chemical element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ElementPattern {
    Any,
    Is(Element),
    OneOf(Vec<Element>),
}

impl ElementPattern {
    pub fn matches(&self, element: Element) -> bool {
        match self {
            Self::Any => true,
            Self::Is(e) => *e == element,
            Self::OneOf(set) => set.contains(&element),
        }
    }
}

/// Pattern on isotope mass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IsotopePattern {
    Any,
    Natural,
    Is(u32),
    OneOf(Vec<u32>),
}

impl IsotopePattern {
    pub fn matches(&self, mass: IsotopeMass) -> bool {
        match self {
            Self::Any => true,
            Self::Natural => mass == IsotopeMass::Natural,
            Self::Is(n) => mass == IsotopeMass::MassNumber(*n),
            Self::OneOf(set) => match mass {
                IsotopeMass::MassNumber(v) => set.contains(&v),
                _ => false,
            },
        }
    }
}

/// Pattern on implicit hydrogen count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HydrogenPattern {
    Any,
    Normal,
    Is(u8),
}

impl HydrogenPattern {
    pub fn matches(&self, n: u8) -> bool {
        match self {
            Self::Any => true,
            Self::Normal => false, // requires context to evaluate
            Self::Is(h) => *h == n,
        }
    }
}

/// Pattern on aromatic valence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticPattern {
    Any,
    NotAromatic,
    Aromatic,
    Is(u8),
}

impl AromaticPattern {
    pub fn matches(&self, av: AromaticValence) -> bool {
        match self {
            Self::Any => true,
            Self::NotAromatic => av == AromaticValence::NotAromatic,
            Self::Aromatic => av.is_aromatic(),
            Self::Is(n) => av == AromaticValence::Valence(*n),
        }
    }
}

/// Atom pattern for use in `MoleculeBuilder` graph nodes and queries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPattern {
    pub element: ElementPattern,
    pub isotope_mass: IsotopePattern,
    pub charge: Pattern<i8>,
    pub implicit_hydrogens: HydrogenPattern,
    pub lone_pairs: Pattern<u8>,
    pub unpaired_electrons: Pattern<u8>,
    pub multiplicity: Pattern<SpinMultiplicity>,
    pub valence: Pattern<u8>,
    pub donated_pairs: Pattern<u8>,
    pub accepted_pairs: Pattern<u8>,
    pub aromatic_valence: AromaticPattern,
    pub multicenter_valence: Pattern<u8>,
}

impl AtomPattern {
    pub fn new(element: Element) -> Self {
        Self {
            element: ElementPattern::Is(element),
            isotope_mass: IsotopePattern::Any,
            charge: Pattern::Any,
            implicit_hydrogens: HydrogenPattern::Any,
            lone_pairs: Pattern::Any,
            unpaired_electrons: Pattern::Any,
            multiplicity: Pattern::Any,
            valence: Pattern::Any,
            donated_pairs: Pattern::Any,
            accepted_pairs: Pattern::Any,
            aromatic_valence: AromaticPattern::Any,
            multicenter_valence: Pattern::Any,
        }
    }

    pub fn from_atom(atom: &Atom) -> Self {
        let isotope = match atom.isotope_mass() {
            IsotopeMass::Natural => IsotopePattern::Natural,
            IsotopeMass::MassNumber(n) => IsotopePattern::Is(n),
        };
        let aromatic = match atom.aromatic_valence() {
            AromaticValence::NotAromatic => AromaticPattern::NotAromatic,
            AromaticValence::Valence(n) => AromaticPattern::Is(n),
        };
        Self {
            element: ElementPattern::Is(atom.element()),
            isotope_mass: isotope,
            charge: Pattern::Is(atom.charge()),
            implicit_hydrogens: HydrogenPattern::Is(atom.implicit_hydrogens()),
            lone_pairs: Pattern::Is(atom.lone_pairs()),
            unpaired_electrons: Pattern::Is(atom.unpaired_electrons()),
            multiplicity: Pattern::Is(atom.multiplicity()),
            valence: Pattern::Is(atom.valence()),
            donated_pairs: Pattern::Is(atom.donated_pairs()),
            accepted_pairs: Pattern::Is(atom.accepted_pairs()),
            aromatic_valence: aromatic,
            multicenter_valence: Pattern::Is(atom.multicenter_valence()),
        }
    }

    /// Create a pattern from a table IR atom.
    ///
    /// Computed fields (valence, donated/accepted pairs, aromatic valence,
    /// multicenter valence) are left as `Any` — filled in during resolution.
    pub fn from_table_atom(atom: &TableAtom) -> Self {
        Self {
            element: ElementPattern::Is(atom.element),
            isotope_mass: atom
                .isotope_mass
                .map_or(IsotopePattern::Any, IsotopePattern::Is),
            charge: atom.charge.map_or(Pattern::Any, Pattern::Is),
            implicit_hydrogens: match atom.implicit_hydrogens {
                Some(ImplicitHydrogens::Normal) => HydrogenPattern::Normal,
                Some(ImplicitHydrogens::Hydrogens(h)) => HydrogenPattern::Is(h),
                None => HydrogenPattern::Any,
            },
            lone_pairs: atom.lone_pairs.map_or(Pattern::Any, Pattern::Is),
            unpaired_electrons: atom.unpaired_electrons.map_or(Pattern::Any, Pattern::Is),
            multiplicity: atom.multiplicity.map_or(Pattern::Any, Pattern::Is),
            valence: Pattern::Any,
            donated_pairs: Pattern::Any,
            accepted_pairs: Pattern::Any,
            aromatic_valence: AromaticPattern::Any,
            multicenter_valence: Pattern::Any,
        }
    }

    /// Create a pattern from a builder atom, incorporating bond-graph-derived valence.
    ///
    /// The resulting pattern is suitable for registry lookup via
    /// `AtomTypeRegistry::candidates_for`.
    pub fn from_builder_atom(builder: &MoleculeBuilder, atom_index: AtomIndex) -> Self {
        let atom = builder.atom(atom_index).expect("atom_index must be valid");
        let valence = builder.atom_bond_order_sum(atom_index);
        let (donated_pairs, accepted_pairs) = builder.atom_dative_bond_order_sums(atom_index);
        let implicit_hydrogens = if builder.atom_has_normal_implicit_hydrogens(atom_index) {
            HydrogenPattern::Normal
        } else {
            atom.implicit_hydrogens.clone()
        };
        let aromatic_valence = if builder.atom_aromatic_hint(atom_index) {
            AromaticPattern::Aromatic
        } else if builder.atom_explicit_aromatic_hint(atom_index) == Some(false) {
            AromaticPattern::NotAromatic
        } else {
            AromaticPattern::Any
        };
        let multicenter_valence = if builder.atom_has_multicenter_bonds(atom_index) {
            Pattern::Any
        } else {
            Pattern::Is(0)
        };
        Self {
            element: atom.element.clone(),
            isotope_mass: IsotopePattern::Any,
            charge: atom.charge,
            implicit_hydrogens,
            lone_pairs: atom.lone_pairs,
            unpaired_electrons: atom.unpaired_electrons,
            multiplicity: atom.multiplicity,
            valence: Pattern::Is(valence),
            donated_pairs: Pattern::Is(donated_pairs),
            accepted_pairs: Pattern::Is(accepted_pairs),
            aromatic_valence,
            multicenter_valence,
        }
    }

    /// Extract the literal element. Panics if the element pattern is not `Is`.
    ///
    /// All atoms stored in a `MoleculeBuilder` graph are created with `Is`
    /// elements, so this is safe for graph-node patterns.
    pub fn element(&self) -> Element {
        match &self.element {
            ElementPattern::Is(e) => *e,
            _ => panic!("atom pattern has non-literal element expression"),
        }
    }

    /// Validate ground fields; non-ground (`Any`, `Normal`, `Aromatic`) fields are skipped.
    ///
    /// The electron invariant is only checked when all relevant fields are ground.
    pub fn check_invariants(&self) -> Result<(), AtomError> {
        let element = match &self.element {
            ElementPattern::Is(e) => *e,
            _ => return Ok(()),
        };
        let charge = match self.charge {
            Pattern::Any => return Ok(()),
            Pattern::Is(c) => c,
        };
        let implicit_hydrogens = match &self.implicit_hydrogens {
            HydrogenPattern::Any | HydrogenPattern::Normal => return Ok(()),
            HydrogenPattern::Is(h) => *h,
        };
        let aromatic_valence: AromaticValence = match &self.aromatic_valence {
            AromaticPattern::Any | AromaticPattern::Aromatic => return Ok(()),
            AromaticPattern::NotAromatic => AromaticValence::NotAromatic,
            AromaticPattern::Is(n) => AromaticValence::Valence(*n),
        };
        let lone_pairs = match self.lone_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let unpaired_electrons = match self.unpaired_electrons {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let multiplicity = match self.multiplicity {
            Pattern::Any => {
                let m = unpaired_electrons.checked_add(1).ok_or_else(|| {
                    AtomError::InvalidMultiplicity((unpaired_electrons as u16 + 1).to_string())
                })?;
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or_else(|| AtomError::InvalidMultiplicity(m.to_string()))?
            }
            Pattern::Is(m) => m,
        };
        let valence = match self.valence {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let donated_pairs = match self.donated_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let accepted_pairs = match self.accepted_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let multicenter_valence = match self.multicenter_valence {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };

        let spin = SpinState::try_new(unpaired_electrons, multiplicity)?;

        let (min_charge, max_charge) = element.charge_bounds();
        if charge < min_charge || charge > max_charge {
            return Err(AtomError::ChargeOutOfBounds {
                element,
                charge,
                min_charge,
                max_charge,
            });
        }

        let max_valence = element.max_valence();
        if valence > max_valence {
            return Err(AtomError::OutOfRange {
                field: "valence",
                value: valence as i64,
                min: 0,
                max: max_valence as i64,
            });
        }

        let max_unpaired_electrons = element.max_unpaired_electrons();
        if spin.unpaired_electrons() > max_unpaired_electrons {
            return Err(AtomError::OutOfRange {
                field: "unpaired_electrons",
                value: spin.unpaired_electrons() as i64,
                min: 0,
                max: max_unpaired_electrons as i64,
            });
        }

        let max_implicit_hydrogens = element.max_implicit_hydrogens();
        if implicit_hydrogens > max_implicit_hydrogens {
            return Err(AtomError::OutOfRange {
                field: "implicit_hydrogens",
                value: implicit_hydrogens as i64,
                min: 0,
                max: max_implicit_hydrogens as i64,
            });
        }

        let aromatic_valence_i16 = aromatic_valence.valence() as i16;
        let aromatic_increment = aromatic_increment(aromatic_valence) as i16;
        let total_e_inv_o = spin.unpaired_electrons() as i16
            + (2 * lone_pairs as i16)
            + (2 * donated_pairs as i16)
            + (2 * accepted_pairs as i16)
            + (2 * implicit_hydrogens as i16)
            + (2 * valence as i16)
            + aromatic_valence_i16
            + aromatic_increment
            + (multicenter_valence as i16);

        let total_e_inv_e = (element.valence_electrons() as i16) - (charge as i16)
            + (implicit_hydrogens as i16)
            + (valence as i16)
            + aromatic_increment
            + (2 * accepted_pairs as i16);

        if total_e_inv_o != total_e_inv_e {
            return Err(AtomError::ElectronInvariantMismatch {
                element,
                orbital_invariant: total_e_inv_o,
                electron_invariant: total_e_inv_e,
            });
        }

        Ok(())
    }

    /// Convert a ground pattern into a concrete atom.
    ///
    /// Fields that are `Any` return 0 for numeric fields,
    /// `NotAromatic` for aromatic valence. Returns `ValidationError::NonGround`
    /// for fields with non-ground patterns that cannot be defaulted:
    /// `ElementPattern::{Any, OneOf}`, `HydrogenPattern::Normal`,
    /// `IsotopePattern::OneOf`, `AromaticPattern::Aromatic`.
    pub fn to_atom(&self) -> Result<Atom, ValidationError> {
        let element = match &self.element {
            ElementPattern::Is(e) => *e,
            _ => return Err(ValidationError::NonGround { field: "element" }),
        };

        let isotope_mass: Option<u32> = match &self.isotope_mass {
            IsotopePattern::Any | IsotopePattern::Natural => None,
            IsotopePattern::Is(n) => Some(*n),
            IsotopePattern::OneOf(_) => {
                return Err(ValidationError::NonGround {
                    field: "isotope_mass",
                })
            }
        };

        let charge = match self.charge {
            Pattern::Any => 0,
            Pattern::Is(c) => c,
        };

        let implicit_hydrogens: u8 = match &self.implicit_hydrogens {
            HydrogenPattern::Any => 0,
            HydrogenPattern::Normal => {
                return Err(ValidationError::NonGround {
                    field: "implicit_hydrogens",
                })
            }
            HydrogenPattern::Is(h) => *h,
        };

        let lone_pairs = match self.lone_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let unpaired_electrons = match self.unpaired_electrons {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };

        let multiplicity = match self.multiplicity {
            Pattern::Any => {
                let m = unpaired_electrons
                    .checked_add(1)
                    .ok_or(ValidationError::InvalidMultiplicity(u8::MAX))?;
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or(ValidationError::InvalidMultiplicity(m))?
            }
            Pattern::Is(m) => m,
        };

        let valence = match self.valence {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let donated_pairs = match self.donated_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };
        let accepted_pairs = match self.accepted_pairs {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };

        let aromatic_valence: AromaticValence = match &self.aromatic_valence {
            AromaticPattern::Any | AromaticPattern::NotAromatic => AromaticValence::NotAromatic,
            AromaticPattern::Is(n) => AromaticValence::Valence(*n),
            AromaticPattern::Aromatic => {
                return Err(ValidationError::NonGround {
                    field: "aromatic_valence",
                })
            }
        };

        let multicenter_valence = match self.multicenter_valence {
            Pattern::Any => 0,
            Pattern::Is(v) => v,
        };

        let atom = Atom::try_new(
            element,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            valence,
            donated_pairs,
            accepted_pairs,
            aromatic_valence,
            multicenter_valence,
        )
        .map_err(ValidationError::Atom)?;

        atom.check_invariants().map_err(ValidationError::Atom)?;

        Ok(atom)
    }

    /// Test whether this pattern matches a concrete atom.
    pub fn matches_atom(&self, atom: &Atom) -> bool {
        self.element.matches(atom.element())
            && self.isotope_mass.matches(atom.isotope_mass())
            && self.charge.matches(atom.charge())
            && self.implicit_hydrogens.matches(atom.implicit_hydrogens())
            && self.lone_pairs.matches(atom.lone_pairs())
            && self.unpaired_electrons.matches(atom.unpaired_electrons())
            && self.multiplicity.matches(atom.multiplicity())
            && self.valence.matches(atom.valence())
            && self.donated_pairs.matches(atom.donated_pairs())
            && self.accepted_pairs.matches(atom.accepted_pairs())
            && self.aromatic_valence.matches(atom.aromatic_valence())
            && self.multicenter_valence.matches(atom.multicenter_valence())
    }
}

impl FromAst<AtomAst> for AtomPattern {
    fn from_ast(ast: AtomAst, cfg: &AtomLowerConfig) -> Result<Self, LoweringError> {
        let element = match ast.element {
            ElementExpr::Lit(e) => ElementPattern::Is(e),
            ElementExpr::Wildcard => ElementPattern::Any,
            ElementExpr::Set(set) => ElementPattern::OneOf(set),
            ElementExpr::Bind { .. } | ElementExpr::Ref(_) => {
                return Err(LoweringError::NonGround { field: "element" })
            }
        };

        let isotope_mass = match ast.isotope_mass.or_else(|| match cfg.isotope_mode {
            IsotopeMode::Normal => Some(IsotopeExpr::Natural),
            IsotopeMode::Provided => None,
        }) {
            None => IsotopePattern::Any,
            Some(IsotopeExpr::Natural) => IsotopePattern::Natural,
            Some(IsotopeExpr::Wildcard) => IsotopePattern::Any,
            Some(IsotopeExpr::Lit(n)) => IsotopePattern::Is(n),
            Some(IsotopeExpr::Set(s)) => IsotopePattern::OneOf(s),
            Some(IsotopeExpr::Bind { .. }) | Some(IsotopeExpr::Ref(_)) => {
                return Err(LoweringError::NonGround {
                    field: "isotope_mass",
                })
            }
        };

        let charge = match ast.charge.or_else(|| match cfg.charge_mode {
            ChargeMode::Zero => Some(ValueAst::Lit(0)),
            ChargeMode::Provided => None,
        }) {
            None => Pattern::Any,
            Some(ValueAst::Wildcard) => Pattern::Any,
            Some(ValueAst::Lit(n)) => Pattern::Is(
                i8::try_from(n).map_err(|_| LoweringError::NonGround { field: "charge" })?,
            ),
            Some(_) => return Err(LoweringError::NonGround { field: "charge" }),
        };

        let implicit_hydrogens =
            match ast
                .implicit_hydrogens
                .or_else(|| match cfg.implicit_h_mode {
                    ImplicitHydrogenMode::Normal => Some(HydrogenExpr::Normal),
                    ImplicitHydrogenMode::Zero => Some(HydrogenExpr::Value(ValueAst::Lit(0))),
                    ImplicitHydrogenMode::Provided => None,
                }) {
                None => HydrogenPattern::Any,
                Some(HydrogenExpr::Normal) => HydrogenPattern::Normal,
                Some(HydrogenExpr::Value(ValueAst::Wildcard)) => HydrogenPattern::Any,
                Some(HydrogenExpr::Value(ValueAst::Lit(n))) => {
                    HydrogenPattern::Is(u8::try_from(n).map_err(|_| LoweringError::NonGround {
                        field: "implicit_hydrogens",
                    })?)
                }
                Some(HydrogenExpr::Value(_)) => {
                    return Err(LoweringError::NonGround {
                        field: "implicit_hydrogens",
                    })
                }
            };

        let aromatic_valence = match ast.aromatic_valence.or_else(|| match cfg.aromatic_mode {
            AromaticMode::None => Some(AromaticExpr::None),
            AromaticMode::Any => Some(AromaticExpr::Value(ValueAst::Wildcard)),
            AromaticMode::Provided => None,
        }) {
            None => AromaticPattern::Any,
            Some(AromaticExpr::None) => AromaticPattern::NotAromatic,
            Some(AromaticExpr::Value(ValueAst::Wildcard)) => AromaticPattern::Aromatic,
            Some(AromaticExpr::Value(ValueAst::Lit(n))) => {
                AromaticPattern::Is(u8::try_from(n).map_err(|_| LoweringError::NonGround {
                    field: "aromatic_valence",
                })?)
            }
            Some(AromaticExpr::Value(_)) => {
                return Err(LoweringError::NonGround {
                    field: "aromatic_valence",
                })
            }
        };

        let multiplicity = match ast.multiplicity {
            None => Pattern::Any,
            Some(ValueAst::Wildcard) => Pattern::Any,
            Some(ValueAst::Lit(n)) => {
                let m = u8::try_from(n).map_err(|_| LoweringError::NonGround {
                    field: "multiplicity",
                })?;
                Pattern::Is(
                    SpinMultiplicity::from_multiplicity(m)
                        .ok_or(LoweringError::InvalidMultiplicity(m))?,
                )
            }
            Some(_) => {
                return Err(LoweringError::NonGround {
                    field: "multiplicity",
                })
            }
        };

        let lower_u8_opt =
            |v: Option<ValueAst>, field: &'static str| -> Result<Pattern<u8>, LoweringError> {
                match v {
                    None => Ok(Pattern::Any),
                    Some(ValueAst::Wildcard) => Ok(Pattern::Any),
                    Some(ValueAst::Lit(n)) => u8::try_from(n)
                        .map(Pattern::Is)
                        .map_err(|_| LoweringError::NonGround { field }),
                    Some(_) => Err(LoweringError::NonGround { field }),
                }
            };

        Ok(AtomPattern {
            element,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs: lower_u8_opt(ast.lone_pairs, "lone_pairs")?,
            unpaired_electrons: lower_u8_opt(ast.unpaired_electrons, "unpaired_electrons")?,
            multiplicity,
            valence: lower_u8_opt(ast.valence, "valence")?,
            donated_pairs: lower_u8_opt(ast.donated_pairs, "donated_pairs")?,
            accepted_pairs: lower_u8_opt(ast.accepted_pairs, "accepted_pairs")?,
            aromatic_valence,
            multicenter_valence: lower_u8_opt(ast.multicenter_valence, "multicenter_valence")?,
        })
    }
}

impl FromAst<AtomAst> for Atom {
    fn from_ast(ast: AtomAst, cfg: &AtomLowerConfig) -> Result<Self, LoweringError> {
        let pattern = AtomPattern::from_ast(ast, cfg)?;
        pattern.to_atom().map_err(|e| match e {
            ValidationError::NonGround { field } => LoweringError::NonGround { field },
            ValidationError::InvalidMultiplicity(n) => LoweringError::InvalidMultiplicity(n),
            ValidationError::Atom(ae) => match ae {
                AtomError::SpinState(se) => LoweringError::SpinState(se),
                other => LoweringError::Atom(other.to_string()),
            },
            other => LoweringError::Atom(other.to_string()),
        })
    }
}

impl FromStr for AtomPattern {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_atom_dsl(s).map_err(|e| LoweringError::Atom(e.to_string()))?;
        Self::from_ast(ast, &AtomLowerConfig::default())
    }
}

fn aromatic_increment(aromatic_valence: AromaticValence) -> u8 {
    match aromatic_valence {
        AromaticValence::NotAromatic => 0,
        AromaticValence::Valence(0) => 0,
        AromaticValence::Valence(1) => 1,
        AromaticValence::Valence(_) => 0,
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{Element, SpinMultiplicity, SpinStateError};

    use super::*;
    use crate::graph_ir::atom_type::AtomError;
    use crate::graph_ir::error::ValidationError;

    #[rstest]
    #[case::defaults(AtomPattern { implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) }, Element::C, 0, 4, 0, SpinMultiplicity::Singlet)]
    fn test_atom_pattern_to_atom(
        #[case] pattern: AtomPattern,
        #[case] element: Element,
        #[case] charge: i8,
        #[case] implicit_hydrogens: u8,
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: SpinMultiplicity,
    ) {
        let atom = pattern.to_atom().unwrap();
        assert_eq!(atom.element(), element);
        assert_eq!(atom.charge(), charge);
        assert_eq!(atom.implicit_hydrogens(), implicit_hydrogens);
        assert_eq!(atom.unpaired_electrons(), unpaired_electrons);
        assert_eq!(atom.multiplicity(), multiplicity);
    }

    #[rstest]
    #[case::invalid_spin_state(
        AtomPattern { unpaired_electrons: Pattern::Is(2), multiplicity: Pattern::Is(SpinMultiplicity::Quartet), ..AtomPattern::new(Element::C) },
        ValidationError::Atom(AtomError::SpinState(SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::Quartet }))
    )]
    #[case::invariant_mismatch(
        AtomPattern { valence: Pattern::Is(2), ..AtomPattern::new(Element::O) },
        ValidationError::Atom(AtomError::ElectronInvariantMismatch { element: Element::O, orbital_invariant: 4, electron_invariant: 8 })
    )]
    fn test_atom_pattern_to_atom_error(
        #[case] pattern: AtomPattern,
        #[case] expected: ValidationError,
    ) {
        let result = pattern.to_atom();
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), expected);
    }
}
