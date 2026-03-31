//! Atom type for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::{Element, SpinMultiplicity, SpinState};

use super::ast_utils::{raise_i8, raise_spin_m_ground, raise_spin_u_ground, raise_u8};
use crate::atom::{AromaticValence, IsotopeMass};
use crate::dsl::ast::{FromAst, ToAst};
use crate::dsl::atom::{parse_atom_dsl, AtomAst};
use crate::dsl::config::{AromaticValenceMode, AtomDslConfig, ImplicitHydrogenMode, IsotopeMode};
use crate::dsl::error::LoweringError;
use crate::dsl::predicates::{AromaticExpr, ElementExpr, HydrogenExpr, IsotopeExpr};
use crate::dsl::value::ValueAst;
use crate::graph_ir::atom_pattern::AtomPattern;
use crate::graph_ir::atom_type::AtomError;
use crate::graph_ir::error::ValidationError;

/// Atom in GraphIR (ground term)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atom {
    element: Element,
    isotope_mass: IsotopeMass,
    charge: i8,
    implicit_hydrogens: u8,
    lone_pairs: u8,
    spin: SpinState,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    multicenter_valence: u8,
}

impl Atom {
    /// Construct a new `Atom`, validating only spin compatibility.
    /// Full electron invariant validation requires `check_invariants`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_new(
        element: Element,
        isotope_mass: Option<u32>,
        charge: i8,
        implicit_hydrogens: u8,
        lone_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
        valence: u8,
        donated_pairs: u8,
        accepted_pairs: u8,
        aromatic_valence: AromaticValence,
        multicenter_valence: u8,
    ) -> Result<Self, AtomError> {
        let spin = SpinState::try_new(unpaired_electrons, multiplicity)?;
        Ok(Self {
            element,
            isotope_mass: isotope_mass.map_or(IsotopeMass::Natural, IsotopeMass::MassNumber),
            charge,
            implicit_hydrogens,
            lone_pairs,
            spin,
            valence,
            donated_pairs,
            accepted_pairs,
            aromatic_valence,
            multicenter_valence,
        })
    }

    // TODO: Add constructors that accept constraints

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn isotope_mass(&self) -> IsotopeMass {
        self.isotope_mass
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }

    pub fn lone_pairs(&self) -> u8 {
        self.lone_pairs
    }

    pub fn spin(&self) -> SpinState {
        self.spin
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.spin.unpaired_electrons()
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.spin.multiplicity()
    }

    pub fn valence(&self) -> u8 {
        self.valence
    }

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
    }

    pub fn aromatic_valence(&self) -> AromaticValence {
        self.aromatic_valence
    }

    pub fn multicenter_valence(&self) -> u8 {
        self.multicenter_valence
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic_valence.is_aromatic()
    }

    /// Validate the electron invariant and field bounds of this atom.
    /// TODO: Integrate into constraint resolution framework.
    pub fn check_invariants(&self) -> Result<(), AtomError> {
        let (min_charge, max_charge) = self.element.charge_bounds();
        if self.charge < min_charge || self.charge > max_charge {
            return Err(AtomError::ChargeOutOfBounds {
                element: self.element,
                charge: self.charge,
                min_charge,
                max_charge,
            });
        }

        let max_valence = self.element.max_valence();
        if self.valence > max_valence {
            return Err(AtomError::OutOfRange {
                field: "valence",
                value: self.valence as i64,
                min: 0,
                max: max_valence as i64,
            });
        }

        let unpaired_electrons = self.spin.unpaired_electrons();
        let max_unpaired_electrons = self.element.max_unpaired_electrons();
        if unpaired_electrons > max_unpaired_electrons {
            return Err(AtomError::OutOfRange {
                field: "unpaired_electrons",
                value: unpaired_electrons as i64,
                min: 0,
                max: max_unpaired_electrons as i64,
            });
        }

        let max_implicit_hydrogens = self.element.max_implicit_hydrogens();
        if self.implicit_hydrogens > max_implicit_hydrogens {
            return Err(AtomError::OutOfRange {
                field: "implicit_hydrogens",
                value: self.implicit_hydrogens as i64,
                min: 0,
                max: max_implicit_hydrogens as i64,
            });
        }

        let aromatic_valence = self.aromatic_valence.valence() as i16;
        let aromatic_increment = aromatic_increment(self.aromatic_valence) as i16;
        let total_e_inv_o = unpaired_electrons as i16
            + (2 * self.lone_pairs as i16)
            + (2 * self.donated_pairs as i16)
            + (2 * self.accepted_pairs as i16)
            + (2 * self.implicit_hydrogens as i16)
            + (2 * self.valence as i16)
            + aromatic_valence
            + aromatic_increment
            + (self.multicenter_valence as i16);

        let total_e_inv_e = (self.element.valence_electrons() as i16) - (self.charge as i16)
            + (self.implicit_hydrogens as i16)
            + (self.valence as i16)
            + aromatic_increment
            + (2 * self.accepted_pairs as i16);

        if total_e_inv_o != total_e_inv_e {
            return Err(AtomError::ElectronInvariantMismatch {
                element: self.element,
                orbital_invariant: total_e_inv_o,
                electron_invariant: total_e_inv_e,
            });
        }

        Ok(())
    }
}

impl FromAst<AtomAst> for Atom {
    fn from_ast(ast: AtomAst, cfg: &AtomDslConfig) -> Result<Self, LoweringError> {
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

impl ToAst<AtomAst> for Atom {
    fn to_ast(&self, cfg: &AtomDslConfig) -> AtomAst {
        let u = self.unpaired_electrons();
        let m = self.multiplicity();
        AtomAst {
            element: ElementExpr::Lit(self.element()),
            isotope_mass: match (self.isotope_mass(), &cfg.isotope_mode) {
                (IsotopeMass::Natural, IsotopeMode::Natural) => None,
                (IsotopeMass::Natural, IsotopeMode::Required) => Some(IsotopeExpr::Natural),
                (IsotopeMass::MassNumber(n), _) => Some(IsotopeExpr::Lit(n)),
            },
            charge: raise_i8(self.charge(), &cfg.charge_mode),
            implicit_hydrogens: match (&cfg.implicit_h_mode, self.implicit_hydrogens()) {
                (ImplicitHydrogenMode::Zero, 0) => None,
                (_, n) => Some(HydrogenExpr::Value(ValueAst::Lit(n as i32))),
            },
            lone_pairs: raise_u8(self.lone_pairs(), &cfg.lone_pairs_mode),
            unpaired_electrons: raise_spin_u_ground(
                u,
                m,
                &cfg.unpaired_electrons_mode,
                &cfg.multiplicity_mode,
            ),
            multiplicity: raise_spin_m_ground(
                u,
                m,
                &cfg.unpaired_electrons_mode,
                &cfg.multiplicity_mode,
            ),
            valence: raise_u8(self.valence(), &cfg.valence_mode),
            donated_pairs: raise_u8(self.donated_pairs(), &cfg.donated_pairs_mode),
            accepted_pairs: raise_u8(self.accepted_pairs(), &cfg.accepted_pairs_mode),
            aromatic_valence: match (self.aromatic_valence(), &cfg.aromatic_valence_mode) {
                (AromaticValence::NotAromatic, AromaticValenceMode::NotAromatic) => None,
                // Unspecified → Any → NotAromatic in to_atom(); safe to suppress
                (AromaticValence::NotAromatic, AromaticValenceMode::Required) => None,
                (AromaticValence::NotAromatic, AromaticValenceMode::Aromatic) => {
                    Some(AromaticExpr::NotAromatic)
                }
                (AromaticValence::Valence(n), _) => {
                    Some(AromaticExpr::Value(ValueAst::Lit(n as i32)))
                }
            },
            multicenter_valence: raise_u8(
                self.multicenter_valence(),
                &cfg.multicenter_valence_mode,
            ),
        }
    }
}

impl FromStr for Atom {
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // TODO: Fix error handling
        let ast = parse_atom_dsl(s).map_err(|e| AtomError::InvalidTag(e.to_string()))?;
        Self::from_ast(ast, &AtomDslConfig::zeroed()).map_err(|e| match e {
            LoweringError::SpinState(se) => AtomError::SpinState(se),
            other => AtomError::InvalidTag(other.to_string()),
        })
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast(&AtomDslConfig::zeroed()).fmt(f)
    }
}

impl Serialize for Atom {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Atom {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(SerdeError::custom)
    }
}

// TODO: Combine with identical function in atom_pattern.rs
fn aromatic_increment(aromatic_valence: AromaticValence) -> u8 {
    match aromatic_valence {
        AromaticValence::NotAromatic => 0,
        AromaticValence::Valence(0) => 0,
        AromaticValence::Valence(1) => 1,
        AromaticValence::Valence(2) => 0,
        AromaticValence::Valence(_) => 0,
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{Element, SpinMultiplicity, SpinState};

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::helium(AtomAst::from_element(Element::He), AtomDslConfig::zeroed(), "He".parse::<Atom>().unwrap())]
    #[case::isotope(AtomAst { isotope_mass: Some(IsotopeExpr::Lit(13)), implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(4))),
        ..AtomAst::from_element(Element::C) }, AtomDslConfig::zeroed(), "C#i13#h4".parse::<Atom>().unwrap())]
    #[case::charge(AtomAst { charge: Some(ValueAst::Lit(1)), implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))),
        ..AtomAst::from_element(Element::C) }, AtomDslConfig::zeroed(), "C#c+#h3".parse::<Atom>().unwrap())]
    fn test_atom_from_ast(#[case] ast: AtomAst, #[case] cfg: AtomDslConfig, #[case] expected: Atom) {
        let atom = Atom::from_ast(ast, &cfg).unwrap();
        assert_eq!(atom, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::helium("He".parse::<Atom>().unwrap(), "He".parse::<Atom>().unwrap().to_ast(&AtomDslConfig::zeroed()))]
    #[case::isotope("C#i13#h4".parse::<Atom>().unwrap(), AtomAst { element: ElementExpr::Lit(Element::C), isotope_mass: Some(IsotopeExpr::Lit(13)),
            implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(4))), .."He".parse::<Atom>().unwrap().to_ast(&AtomDslConfig::zeroed()) })]
    #[case::aromatic("C#h#v2#a1".parse::<Atom>().unwrap(), AtomAst { element: ElementExpr::Lit(Element::C),
            implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), valence: Some(ValueAst::Lit(2)), aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))),
            .."He".parse::<Atom>().unwrap().to_ast(&AtomDslConfig::zeroed()) })]
    fn test_atom_to_ast(#[case] atom: Atom, #[case] expected: AtomAst) {
        assert_eq!(atom.to_ast(&AtomDslConfig::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::element("He", Atom { element: Element::He, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 0, lone_pairs: 0,
        spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::whitespace("  He  ", Atom { element: Element::He, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 0, lone_pairs: 0,
        spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::hydrogens("C#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::hydrogens_whitespace("C #h 4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::isotope_mass_number("C#i13#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::MassNumber(13), charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::isotope_natural("C#i=#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::charge_plus("C#c+#h3", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 1, implicit_hydrogens: 3, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::charge_minus("C#c-#h3#n1", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: -1, implicit_hydrogens: 3, lone_pairs: 1,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    #[case::aromatic_none("C#a!#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::NotAromatic, multicenter_valence: 0 })]
    fn test_atom_from_str(
        #[case] input: &str,
        #[case] expected: Atom,
    ) {
        let result = Atom::from_str(input);
        assert!(result.is_ok(), "from_str should succeed for {}", input);
        let atom = result.unwrap();
        assert_eq!(atom, expected, "from_str mismatch for {}", input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_element("X", AtomError::InvalidTag("Invalid atom element: X".to_string()))]
    #[case::wildcard("*", AtomError::InvalidTag("non-ground value for field 'element'".to_string()))]
    #[case::unknown_predicate("C#x1", AtomError::InvalidTag("Unknown atom predicate: #x".to_string()))]
    #[case::duplicate_charge("C#c+#c-", AtomError::InvalidTag("Duplicate #c atom predicate".to_string()))]
    #[case::duplicate_h("C#h3#h2", AtomError::InvalidTag("Duplicate #h atom predicate".to_string()))]
    #[case::non_ground_payload("C#h*", AtomError::InvalidTag("invalid atom spec: electron invariant mismatch for C: inv_o=0, inv_e=4".to_string()))]
    #[case::malformed_number("C#vabc", AtomError::InvalidTag("Trailing input: \"abc\"".to_string()))]
    #[case::trailing_input("C foo", AtomError::InvalidTag("Trailing input: \"foo\"".to_string()))]
    fn test_atom_from_str_invalid(#[case] input: &str, #[case] expected: AtomError) {
        let result = Atom::from_str(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        let err = result.unwrap_err();
        assert_eq!(err, expected, "{input:?} should fail with {expected:?}, got {err:?}"
        );
    }

    #[rstest]
    #[case::defaults("He".parse::<Atom>().unwrap(), "He")]
    #[case::charge_plus("C#c+#h3".parse::<Atom>().unwrap(), "C#c+#h3")]
    #[case::charge_minus("C#c-#h3#n1".parse::<Atom>().unwrap(), "C#c-#h3#n")]
    #[case::isotope_mass("C#i13#h4".parse::<Atom>().unwrap(), "C#i13#h4")]
    #[case::aromatic("C#h#v2#a1".parse::<Atom>().unwrap(), "C#h#v2#a")]
    #[case::multiplicity("C#n#u2#s".parse::<Atom>().unwrap(), "C#n#u2#s")]
    fn test_atom_display(#[case] atom: Atom, #[case] expected: &str) {
        assert_eq!(atom.to_string(), expected);
    }

    #[rstest]
    #[case::defaults("He".parse::<Atom>().unwrap(), r#""He""#)]
    #[case::charge_plus("C#c+#h3".parse::<Atom>().unwrap(), r#""C#c+#h3""#)]
    #[case::charge_minus("C#c-#h3#n1".parse::<Atom>().unwrap(), r#""C#c-#h3#n""#)]
    #[case::isotope_mass("C#i13#h4".parse::<Atom>().unwrap(), r#""C#i13#h4""#)]
    #[case::aromatic("C#h#v2#a1".parse::<Atom>().unwrap(), r#""C#h#v2#a""#)]
    #[case::multiplicity("C#n#u2#s".parse::<Atom>().unwrap(), r#""C#n#u2#s""#)]
    fn test_atom_serialize(#[case] atom: Atom, #[case] expected: &str) {
        let json = serde_json::to_string(&atom).unwrap();
        assert_eq!(json, expected);
    }

    #[rstest]
    #[case::defaults(r#""He""#, "He".parse::<Atom>().unwrap(), )]
    #[case::charge_plus(r#""C#c+#h3""#, "C#c+#h3".parse::<Atom>().unwrap())]
    #[case::charge_minus(r#""C#c-#h3#n1""#, "C#c-#h3#n".parse::<Atom>().unwrap())]
    #[case::isotope_mass(r#""C#i13#h4""#, "C#i13#h4".parse::<Atom>().unwrap())]
    #[case::aromatic(r#""C#h#v2#a1""#, "C#h#v2#a".parse::<Atom>().unwrap())]
    #[case::multiplicity(r#""C#n#u2#s""#, "C#n#u2#s".parse::<Atom>().unwrap())]
    fn test_atom_deserialize(#[case] input: &str, #[case] expected: Atom) {
        let atom: Atom = serde_json::from_str(input).unwrap();
        assert_eq!(atom, expected);
    }
}
