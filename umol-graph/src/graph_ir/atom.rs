//! Atom type for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use umol_shared::atom_ast::{AromaticAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::error::SpinStateError;
use umol_shared::spin::{SpinMultiplicity, SpinState};
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;
use umol_edn::{DeError, Edn, FromEdn, ToEdn};

use super::ast_utils::{raise_i8_ground, raise_spin_ground, raise_u8_ground};
use crate::atom::{AromaticValence, IsotopeMass};
use crate::ast::atom::AtomAst;
use crate::ast::config::{AromaticValenceMode, AtomAstConfig, ImplicitHydrogenMode, IsotopeMode};
use crate::ast::error::LoweringError;
use crate::ast::{FromAst, ToAst};
use crate::dsl::atom::parse_atom_dsl;
use crate::graph_ir::atom_pattern::AtomPattern;
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
    ) -> Result<Self, ValidationError> {
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
    pub fn check_invariants(&self) -> Result<(), ValidationError> {
        let (min_charge, max_charge) = self.element.charge_bounds();
        if self.charge < min_charge || self.charge > max_charge {
            return Err(ValidationError::ChargeOutOfBounds {
                element: self.element,
                charge: self.charge,
                min_charge,
                max_charge,
            });
        }

        let max_valence = self.element.max_valence();
        if self.valence > max_valence {
            return Err(ValidationError::OutOfRange {
                field: "valence",
                value: self.valence as i64,
                min: 0,
                max: max_valence as i64,
            });
        }

        let unpaired_electrons = self.spin.unpaired_electrons();
        let max_unpaired_electrons = self.element.max_unpaired_electrons();
        if unpaired_electrons > max_unpaired_electrons {
            return Err(ValidationError::OutOfRange {
                field: "unpaired_electrons",
                value: unpaired_electrons as i64,
                min: 0,
                max: max_unpaired_electrons as i64,
            });
        }

        let max_implicit_hydrogens = self.element.max_implicit_hydrogens();
        if self.implicit_hydrogens > max_implicit_hydrogens {
            return Err(ValidationError::OutOfRange {
                field: "implicit_hydrogens",
                value: self.implicit_hydrogens as i64,
                min: 0,
                max: max_implicit_hydrogens as i64,
            });
        }

        let aromatic_valence = self.aromatic_valence.valence() as i16;
        let aromatic_increment = self.aromatic_valence.valence_increment() as i16;
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
            return Err(ValidationError::ElectronInvariantMismatch {
                element: self.element,
                orbital_invariant: total_e_inv_o,
                electron_invariant: total_e_inv_e,
            });
        }

        Ok(())
    }
}

impl FromAst<AtomAst> for Atom {
    fn from_ast(ast: &AtomAst, cfg: &AtomAstConfig) -> Result<Self, LoweringError> {
        let pattern = AtomPattern::from_ast(ast, cfg)?;
        pattern.to_atom().map_err(|e| match e {
            ValidationError::NonGround { field } => LoweringError::NonGround { field },
            ValidationError::InvalidMultiplicity(n) => LoweringError::InvalidMultiplicity(n),
            ValidationError::SpinUnderdetermined => {
                LoweringError::SpinState(SpinStateError::Underdetermined)
            }
            ValidationError::SpinIncompatible {
                unpaired_electrons,
                multiplicity,
            } => LoweringError::SpinState(SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            }),
            other => LoweringError::Atom(other.to_string()),
        })
    }
}

impl ToAst<AtomAst> for Atom {
    fn to_ast(&self, cfg: &AtomAstConfig) -> AtomAst {
        let (spin_u, spin_m) = raise_spin_ground(
            self.unpaired_electrons(),
            self.multiplicity(),
            &cfg.unpaired_electrons_mode,
            &cfg.multiplicity_mode,
        );
        AtomAst {
            element: ElementAst::Lit(self.element()),
            isotope_mass: match (self.isotope_mass(), &cfg.isotope_mode) {
                (IsotopeMass::Natural, IsotopeMode::Natural) => None,
                (IsotopeMass::Natural, IsotopeMode::Required) => Some(IsotopeAst::Natural),
                (IsotopeMass::MassNumber(n), _) => Some(IsotopeAst::Lit(n)),
            },
            charge: raise_i8_ground(self.charge(), &cfg.charge_mode),
            implicit_hydrogens: match (&cfg.implicit_h_mode, self.implicit_hydrogens()) {
                (ImplicitHydrogenMode::Zero, 0) => None,
                (_, n) => Some(HydrogenAst::Value(ValueAst::Lit(n as i64))),
            },
            lone_pairs: raise_u8_ground(self.lone_pairs(), &cfg.lone_pairs_mode),
            spin: SpinStateAst::from_pair(spin_u, spin_m),
            valence: raise_u8_ground(self.valence(), &cfg.valence_mode),
            donated_pairs: raise_u8_ground(self.donated_pairs(), &cfg.donated_pairs_mode),
            accepted_pairs: raise_u8_ground(self.accepted_pairs(), &cfg.accepted_pairs_mode),
            aromatic_valence: match (self.aromatic_valence(), &cfg.aromatic_valence_mode) {
                (AromaticValence::NotAromatic, AromaticValenceMode::NotAromatic) => None,
                // Unspecified → Any → NotAromatic in to_atom(); safe to suppress
                (AromaticValence::NotAromatic, AromaticValenceMode::Required) => None,
                (AromaticValence::NotAromatic, AromaticValenceMode::Aromatic) => {
                    Some(AromaticAst::NotAromatic)
                }
                (AromaticValence::Valence(n), _) => {
                    Some(AromaticAst::Value(ValueAst::Lit(n as i64)))
                }
            },
            multicenter_valence: raise_u8_ground(
                self.multicenter_valence(),
                &cfg.multicenter_valence_mode,
            ),
        }
    }
}

impl FromStr for Atom {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_atom_dsl(s).map_err(|e| LoweringError::Atom(e.to_string()))?;
        Self::from_ast(&ast, &AtomAstConfig::zeroed())
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast(&AtomAstConfig::zeroed()).fmt(f)
    }
}

impl<'de> FromEdn<'de> for Atom {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let ast = AtomAst::from_edn(edn)?;
        Self::from_ast(&ast, &AtomAstConfig::zeroed()).map_err(|e| DeError::subgrammar("atom", e))
    }
}

impl ToEdn for Atom {
    fn to_edn(&self) -> Edn<'static> {
        self.to_ast(&AtomAstConfig::zeroed()).to_edn()
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;
    use umol_shared::spin::{SpinMultiplicity, SpinState};

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::helium(AtomAst::from_element(Element::He), AtomAstConfig::zeroed(), "He".parse::<Atom>().unwrap())]
    #[case::isotope(AtomAst { isotope_mass: Some(IsotopeAst::Lit(13)), implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(4))),
        ..AtomAst::from_element(Element::C) }, AtomAstConfig::zeroed(), "C#i13#h4".parse::<Atom>().unwrap())]
    #[case::charge(AtomAst { charge: Some(ValueAst::Lit(1)), implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(3))),
        ..AtomAst::from_element(Element::C) }, AtomAstConfig::zeroed(), "C#c+#h3".parse::<Atom>().unwrap())]
    fn test_atom_from_ast(#[case] ast: AtomAst, #[case] cfg: AtomAstConfig, #[case] expected: Atom) {
        let atom = Atom::from_ast(&ast, &cfg).unwrap();
        assert_eq!(atom, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::helium("He".parse::<Atom>().unwrap(), "He".parse::<Atom>().unwrap().to_ast(&AtomAstConfig::zeroed()))]
    #[case::isotope("C#i13#h4".parse::<Atom>().unwrap(), AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: Some(IsotopeAst::Lit(13)),
            implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(4))), .."He".parse::<Atom>().unwrap().to_ast(&AtomAstConfig::zeroed()) })]
    #[case::aromatic("C#h#v2#a1".parse::<Atom>().unwrap(), AtomAst { element: ElementAst::Lit(Element::C),
            implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(1))), valence: Some(ValueAst::Lit(2)), aromatic_valence: Some(AromaticAst::Value(ValueAst::Lit(1))),
            .."He".parse::<Atom>().unwrap().to_ast(&AtomAstConfig::zeroed()) })]
    fn test_atom_to_ast(#[case] atom: Atom, #[case] expected: AtomAst) {
        assert_eq!(atom.to_ast(&AtomAstConfig::zeroed()), expected);
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
    #[case::unknown_element("X", LoweringError::Atom("Invalid atom element: X".to_string()))]
    #[case::wildcard("*", LoweringError::NonGround { field: "element" })]
    #[case::unknown_predicate("C#x1", LoweringError::Atom("Unknown atom predicate: #x".to_string()))]
    #[case::duplicate_charge("C#c+#c-", LoweringError::Atom("Duplicate #c atom predicate".to_string()))]
    #[case::duplicate_h("C#h3#h2", LoweringError::Atom("Duplicate #h atom predicate".to_string()))]
    #[case::non_ground_payload("C#h*", LoweringError::Atom("electron invariant mismatch for C: inv_o=0, inv_e=4".to_string()))]
    #[case::malformed_number("C#vabc", LoweringError::Atom("Trailing input: \"abc\"".to_string()))]
    #[case::trailing_input("C foo", LoweringError::Atom("Trailing input: \"foo\"".to_string()))]
    fn test_atom_from_str_invalid(#[case] input: &str, #[case] expected: LoweringError) {
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
    fn test_atom_to_edn(#[case] atom: Atom, #[case] expected: &str) {
        let edn = atom.to_edn();
        assert_eq!(edn.to_string(), expected);
    }

    #[rstest]
    #[case::defaults(r#""He""#, "He".parse::<Atom>().unwrap())]
    #[case::charge_plus(r#""C#c+#h3""#, "C#c+#h3".parse::<Atom>().unwrap())]
    #[case::charge_minus(r#""C#c-#h3#n1""#, "C#c-#h3#n".parse::<Atom>().unwrap())]
    #[case::isotope_mass(r#""C#i13#h4""#, "C#i13#h4".parse::<Atom>().unwrap())]
    #[case::aromatic(r#""C#h#v2#a1""#, "C#h#v2#a".parse::<Atom>().unwrap())]
    #[case::multiplicity(r#""C#n#u2#s""#, "C#n#u2#s".parse::<Atom>().unwrap())]
    fn test_atom_from_edn(#[case] input: &str, #[case] expected: Atom) {
        let atom = Atom::from_edn_str(input).unwrap();
        assert_eq!(atom, expected);
    }
}
