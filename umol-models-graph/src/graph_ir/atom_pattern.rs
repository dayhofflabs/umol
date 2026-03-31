//! Atom pattern representation and conversion into concrete ground atoms.

use std::fmt;
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::{Element, SpinMultiplicity, SpinState};

use crate::atom::{AromaticValence, IsotopeMass};
use crate::dsl::ast::{FromAst, ToAst};
use crate::dsl::atom::{
    parse_atom_dsl, AromaticMode, AtomAst, AtomLowerConfig, ChargeMode, ImplicitHydrogenMode,
    IsotopeMode,
};
use crate::dsl::error::LoweringError;
use crate::dsl::predicates::{AromaticExpr, ElementExpr, HydrogenExpr, IsotopeExpr};
use crate::dsl::value::ValueAst;
use crate::graph_ir::atom::Atom;
use crate::graph_ir::atom_type::AtomError;
use crate::graph_ir::error::ValidationError;
use crate::graph_ir::molecule::AtomIndex;
use crate::graph_ir::molecule_builder::MoleculeBuilder;
use crate::table_ir::atom::{Atom as TableAtom, ImplicitHydrogens};

/// Generic pattern for a scalar-valued field: unconstrained or an exact value.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum HydrogenPattern {
    Any,
    Normal,
    Is(u8),
}

// Normal is a deferred constraint — must be resolved before matching
impl HydrogenPattern {
    pub fn matches(&self, n: u8) -> bool {
        match self {
            Self::Any => true,
            Self::Normal => false,
            Self::Is(h) => *h == n,
        }
    }
}

/// Pattern on aromatic valence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
            Self::Aromatic => av.is_aromatic(),
            Self::NotAromatic => av == AromaticValence::NotAromatic,
            Self::Is(n) => av == AromaticValence::Valence(*n),
        }
    }
}

/// Atom pattern term.
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
    /// `IsotopePattern::OneOf`, `AromaticPattern::Any`.
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
            AromaticPattern::Aromatic => {
                return Err(ValidationError::NonGround {
                    field: "aromatic_valence",
                })
            }
            AromaticPattern::Is(n) => AromaticValence::Valence(*n),
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

fn aromatic_increment(aromatic_valence: AromaticValence) -> u8 {
    match aromatic_valence {
        AromaticValence::NotAromatic => 0,
        AromaticValence::Valence(0) => 0,
        AromaticValence::Valence(1) => 1,
        AromaticValence::Valence(_) => 0,
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
            IsotopeMode::Natural => Some(IsotopeExpr::Natural),
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
            AromaticMode::None => Some(AromaticExpr::NotAromatic),
            AromaticMode::Any => Some(AromaticExpr::Unspecified),
            AromaticMode::Provided => None,
        }) {
            None => AromaticPattern::Any,
            Some(AromaticExpr::Unspecified) => AromaticPattern::Any,
            Some(AromaticExpr::NotAromatic) => AromaticPattern::NotAromatic,
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

impl ToAst<AtomAst> for AtomPattern {
    fn to_ast(&self) -> AtomAst {
        AtomAst {
            element: match &self.element {
                ElementPattern::Any => ElementExpr::Wildcard,
                ElementPattern::Is(e) => ElementExpr::Lit(*e),
                ElementPattern::OneOf(es) => ElementExpr::Set(es.clone()),
            },
            isotope_mass: match &self.isotope_mass {
                IsotopePattern::Any => Some(IsotopeExpr::Wildcard),
                IsotopePattern::Natural => Some(IsotopeExpr::Natural),
                IsotopePattern::Is(n) => Some(IsotopeExpr::Lit(*n)),
                IsotopePattern::OneOf(ns) => Some(IsotopeExpr::Set(ns.clone())),
            },
            charge: match self.charge {
                Pattern::Any => Some(ValueAst::Wildcard),
                Pattern::Is(n) => Some(ValueAst::Lit(n as i32)),
            },
            implicit_hydrogens: match &self.implicit_hydrogens {
                HydrogenPattern::Any => Some(HydrogenExpr::Value(ValueAst::Wildcard)),
                HydrogenPattern::Normal => Some(HydrogenExpr::Normal),
                HydrogenPattern::Is(n) => Some(HydrogenExpr::Value(ValueAst::Lit(*n as i32))),
            },
            lone_pairs: wildcard_u8(self.lone_pairs),
            unpaired_electrons: wildcard_u8(self.unpaired_electrons),
            multiplicity: match self.multiplicity {
                Pattern::Any => Some(ValueAst::Wildcard),
                Pattern::Is(m) => Some(ValueAst::Lit(m.multiplicity() as i32)),
            },
            valence: wildcard_u8(self.valence),
            donated_pairs: wildcard_u8(self.donated_pairs),
            accepted_pairs: wildcard_u8(self.accepted_pairs),
            aromatic_valence: match &self.aromatic_valence {
                AromaticPattern::Any => Some(AromaticExpr::Unspecified),
                AromaticPattern::NotAromatic => Some(AromaticExpr::NotAromatic),
                AromaticPattern::Aromatic => Some(AromaticExpr::Value(ValueAst::Wildcard)),
                AromaticPattern::Is(n) => Some(AromaticExpr::Value(ValueAst::Lit(*n as i32))),
            },
            multicenter_valence: wildcard_u8(self.multicenter_valence),
        }
    }
}

fn wildcard_u8(pat: Pattern<u8>) -> Option<ValueAst> {
    match pat {
        Pattern::Any => Some(ValueAst::Wildcard),
        Pattern::Is(n) => Some(ValueAst::Lit(n as i32)),
    }
}

impl FromStr for AtomPattern {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_atom_dsl(s).map_err(|e| LoweringError::Atom(e.to_string()))?;
        Self::from_ast(ast, &AtomLowerConfig::default())
    }
}

impl fmt::Display for AtomPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast().fmt(f)
    }
}

impl Serialize for AtomPattern {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for AtomPattern {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        let ast = parse_atom_dsl(&s).map_err(SerdeError::custom)?;
        AtomPattern::from_ast(
            ast,
            &AtomLowerConfig {
                isotope_mode: IsotopeMode::Provided,
                charge_mode: ChargeMode::Provided,
                implicit_h_mode: ImplicitHydrogenMode::Provided,
                aromatic_mode: AromaticMode::Provided,
            },
        )
        .map_err(SerdeError::custom)
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

    #[rustfmt::skip]
    #[rstest]
    #[case::defaults(
        AtomAst::from_element(Element::C),
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }
    )]
    #[case::charge_and_hydrogens(
        AtomAst { charge: Some(ValueAst::Lit(1)), implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), ..AtomPattern::new(Element::C) }
    )]
    #[case::wildcard_charge(
        AtomAst { charge: Some(ValueAst::Wildcard), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }
    )]
    #[case::aromatic_wildcard(
        AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Wildcard)), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::Aromatic, ..AtomPattern::new(Element::C) }
    )]
    #[case::aromatic_unspecified(
        AtomAst { aromatic_valence: Some(AromaticExpr::Unspecified), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }
    )]
    #[case::aromatic_not_aromatic(
        AtomAst { aromatic_valence: Some(AromaticExpr::NotAromatic), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) }
    )]
    #[case::aromatic_specific(
        AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(2))), ..AtomAst::from_element(Element::C) },
        AtomLowerConfig::default(),
        AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::Is(2), ..AtomPattern::new(Element::C) }
    )]
    #[case::absent_aromatic_mode_any(
        AtomAst::from_element(Element::C),
        AtomLowerConfig { aromatic_mode: AromaticMode::Any, ..AtomLowerConfig::default() },
        AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }
    )]
    #[case::absent_aromatic_mode_none(
        AtomAst::from_element(Element::C),
        AtomLowerConfig { aromatic_mode: AromaticMode::None, ..AtomLowerConfig::default() },
        AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) }
    )]
    fn test_atom_pattern_from_ast(
        #[case] ast: AtomAst,
        #[case] cfg: AtomLowerConfig,
        #[case] expected: AtomPattern,
    ) {
        assert_eq!(AtomPattern::from_ast(ast, &cfg).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_any(
        AtomPattern::new(Element::C),
        AtomAst {
            element: ElementExpr::Lit(Element::C),
            isotope_mass: Some(IsotopeExpr::Wildcard),
            charge: Some(ValueAst::Wildcard),
            implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Wildcard)),
            lone_pairs: Some(ValueAst::Wildcard),
            unpaired_electrons: Some(ValueAst::Wildcard),
            multiplicity: Some(ValueAst::Wildcard),
            valence: Some(ValueAst::Wildcard),
            donated_pairs: Some(ValueAst::Wildcard),
            accepted_pairs: Some(ValueAst::Wildcard),
            aromatic_valence: Some(AromaticExpr::Unspecified),
            multicenter_valence: Some(ValueAst::Wildcard),
        }
    )]
    #[case::ground(
        AtomPattern { isotope_mass: IsotopePattern::Natural, charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), aromatic_valence: AromaticPattern::Is(1), ..AtomPattern::new(Element::C) },
        AtomAst {
            charge: Some(ValueAst::Lit(1)),
            implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))),
            aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))),
            isotope_mass: Some(IsotopeExpr::Natural),
            .."C".parse::<AtomPattern>().unwrap().to_ast()
        }
    )]
    #[case::aromatic_variants(
        AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::Aromatic, ..AtomPattern::new(Element::C) },
        AtomAst {
            aromatic_valence: Some(AromaticExpr::Value(ValueAst::Wildcard)),
            isotope_mass: Some(IsotopeExpr::Natural),
            .."C".parse::<AtomPattern>().unwrap().to_ast()
        }
    )]
    fn test_atom_pattern_to_ast(#[case] pattern: AtomPattern, #[case] expected: AtomAst) {
        assert_eq!(pattern.to_ast(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::element("C", AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) })]
    #[case::hydrogens("C#h4", AtomPattern { isotope_mass: IsotopePattern::Natural, implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) })]
    #[case::charge_plus("C#c+#h3", AtomPattern { isotope_mass: IsotopePattern::Natural, charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), ..AtomPattern::new(Element::C) })]
    #[case::isotope("C#i13#h4", AtomPattern { isotope_mass: IsotopePattern::Is(13), implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) })]
    #[case::aromatic_none("C#a!", AtomPattern { isotope_mass: IsotopePattern::Natural, aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) })]
    #[case::wildcard("*", AtomPattern { element: ElementPattern::Any, isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) })]
    fn test_atom_pattern_from_str(#[case] input: &str, #[case] expected: AtomPattern) {
        assert_eq!(input.parse::<AtomPattern>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_any(AtomPattern::new(Element::C), "C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a?#m*")]
    #[case::isotope_natural(AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }, "C#i=#c*#h*#n*#u*#s*#v*#d*#r*#a?#m*")]
    #[case::hydrogens(AtomPattern { implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) }, "C#i*#c*#h4#n*#u*#s*#v*#d*#r*#a?#m*")]
    #[case::charge_plus(AtomPattern { charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), ..AtomPattern::new(Element::C) }, "C#i*#c+#h3#n*#u*#s*#v*#d*#r*#a?#m*")]
    #[case::isotope_mass(AtomPattern { isotope_mass: IsotopePattern::Is(13), implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) }, "C#i13#c*#h4#n*#u*#s*#v*#d*#r*#a?#m*")]
    #[case::aromatic_none(AtomPattern { aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) }, "C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a!#m*")]
    #[case::aromatic(AtomPattern { aromatic_valence: AromaticPattern::Aromatic, ..AtomPattern::new(Element::C) }, "C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a*#m*")]
    fn test_atom_pattern_display(#[case] pattern: AtomPattern, #[case] expected: &str) {
        assert_eq!(pattern.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_any(AtomPattern::new(Element::C), r#""C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a?#m*""#)]
    #[case::isotope_natural(AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) }, r#""C#i=#c*#h*#n*#u*#s*#v*#d*#r*#a?#m*""#)]
    #[case::hydrogens(AtomPattern { implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) }, r#""C#i*#c*#h4#n*#u*#s*#v*#d*#r*#a?#m*""#)]
    #[case::charge_plus(AtomPattern { charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), ..AtomPattern::new(Element::C) }, r#""C#i*#c+#h3#n*#u*#s*#v*#d*#r*#a?#m*""#)]
    #[case::isotope_mass(AtomPattern { isotope_mass: IsotopePattern::Is(13), implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) }, r#""C#i13#c*#h4#n*#u*#s*#v*#d*#r*#a?#m*""#)]
    #[case::aromatic_none(AtomPattern { aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) }, r#""C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a!#m*""#)]
    #[case::aromatic(AtomPattern { aromatic_valence: AromaticPattern::Aromatic, ..AtomPattern::new(Element::C) }, r#""C#i*#c*#h*#n*#u*#s*#v*#d*#r*#a*#m*""#)]
    fn test_atom_pattern_serialize(#[case] pattern: AtomPattern, #[case] expected: &str) {
        let json = serde_json::to_string(&pattern).unwrap();
        assert_eq!(json, expected);
    }

    // deserialize uses IsotopeMode::Provided: absent isotope → IsotopePattern::Any
    #[rustfmt::skip]
    #[rstest]
    #[case::all_any(r#""C""#, AtomPattern::new(Element::C))]
    #[case::isotope_natural(r#""C#i=""#, AtomPattern { isotope_mass: IsotopePattern::Natural, ..AtomPattern::new(Element::C) })]
    #[case::hydrogens(r#""C#h4""#, AtomPattern { implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) })]
    #[case::charge_plus(r#""C#c+#h3""#, AtomPattern { charge: Pattern::Is(1), implicit_hydrogens: HydrogenPattern::Is(3), ..AtomPattern::new(Element::C) })]
    #[case::isotope_mass(r#""C#i13#h4""#, AtomPattern { isotope_mass: IsotopePattern::Is(13), implicit_hydrogens: HydrogenPattern::Is(4), ..AtomPattern::new(Element::C) })]
    #[case::aromatic_dontcare(r#""C#a?""#, AtomPattern::new(Element::C))]
    #[case::aromatic_none(r#""C#a!""#, AtomPattern { aromatic_valence: AromaticPattern::NotAromatic, ..AtomPattern::new(Element::C) })]
    #[case::aromatic(r#""C#a*""#, AtomPattern { aromatic_valence: AromaticPattern::Aromatic, ..AtomPattern::new(Element::C) })]
    fn test_atom_pattern_deserialize(#[case] input: &str, #[case] expected: AtomPattern) {
        let pattern: AtomPattern = serde_json::from_str(input).unwrap();
        assert_eq!(pattern, expected);
    }
}
