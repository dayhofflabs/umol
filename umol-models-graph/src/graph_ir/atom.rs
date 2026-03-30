//! Atom type and builder for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use smallvec::SmallVec;
use umol_data::{Element, SpinMultiplicity, SpinState};

use self::parser::parse_ground_atom_dsl;
use crate::atom::{AromaticValence, Chirality, ImplicitHydrogens, IsotopeMass};
use crate::graph_ir::atom_pattern::AtomPattern;
use crate::graph_ir::atom_type::{AtomError, AtomTypeSpec};
use crate::graph_ir::error::ResolutionError;
use crate::table_ir::atom::Atom as TableAtom;
mod parser;

/// Atom in GraphIR (ground term)
#[derive(Debug, Clone, PartialEq)]
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

    pub fn to_spec(&self) -> AtomTypeSpec {
        AtomTypeSpec::new(
            self.element,
            self.isotope_mass.mass_number(),
            self.charge,
            self.implicit_hydrogens,
            self.lone_pairs,
            self.spin.unpaired_electrons(),
            self.spin.multiplicity(),
            self.valence,
            self.donated_pairs,
            self.accepted_pairs,
            self.aromatic_valence,
            self.multicenter_valence,
        )
        .expect("resolved Atom fields always form a valid AtomTypeSpec")
    }

    pub fn to_builder(&self) -> AtomBuilder {
        let spec = self.to_spec();
        AtomBuilder {
            element: self.element,
            isotope_mass: self.isotope_mass.mass_number(),
            charge: Some(self.charge),
            implicit_hydrogens: Some(ImplicitHydrogens::Hydrogens(self.implicit_hydrogens)),
            lone_pairs: Some(self.lone_pairs),
            unpaired_electrons: Some(self.spin.unpaired_electrons()),
            multiplicity: Some(self.spin.multiplicity()),
            aromatic_hint: None,
            chirality_hint: None,
            candidates: SmallVec::from_elem(spec, 1),
        }
    }

    pub(crate) fn from_spec(spec: AtomTypeSpec) -> Self {
        Self {
            element: spec.element(),
            isotope_mass: spec
                .isotope_mass()
                .map_or(IsotopeMass::Natural, IsotopeMass::MassNumber),
            charge: spec.charge(),
            implicit_hydrogens: spec.implicit_hydrogens(),
            lone_pairs: spec.lone_pairs(),
            spin: spec.spin(),
            valence: spec.valence(),
            donated_pairs: spec.donated_pairs(),
            accepted_pairs: spec.accepted_pairs(),
            aromatic_valence: spec.aromatic_valence(),
            multicenter_valence: spec.multicenter_valence(),
        }
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.element)?;

        if let IsotopeMass::MassNumber(m) = self.isotope_mass {
            write!(f, "#i{}", m)?;
        }

        match self.charge {
            0 => {}
            1 => write!(f, "#c+")?,
            -1 => write!(f, "#c-")?,
            c if c > 0 => write!(f, "#c+{}", c)?,
            c => write!(f, "#c{}", c)?,
        }

        if self.implicit_hydrogens > 0 {
            if self.implicit_hydrogens == 1 {
                write!(f, "#h")?;
            } else {
                write!(f, "#h{}", self.implicit_hydrogens)?;
            }
        }

        if self.lone_pairs > 0 {
            if self.lone_pairs == 1 {
                write!(f, "#n")?;
            } else {
                write!(f, "#n{}", self.lone_pairs)?;
            }
        }

        let unpaired = self.spin.unpaired_electrons();
        if unpaired > 0 {
            if unpaired == 1 {
                write!(f, "#u")?;
            } else {
                write!(f, "#u{}", unpaired)?;
            }
        }

        let multiplicity = self.spin.multiplicity().multiplicity();
        if multiplicity != unpaired + 1 {
            if multiplicity == 1 {
                write!(f, "#s")?;
            } else {
                write!(f, "#s{}", multiplicity)?;
            }
        }

        if self.valence > 0 {
            if self.valence == 1 {
                write!(f, "#v")?;
            } else {
                write!(f, "#v{}", self.valence)?;
            }
        }

        if self.donated_pairs > 0 {
            if self.donated_pairs == 1 {
                write!(f, "#d")?;
            } else {
                write!(f, "#d{}", self.donated_pairs)?;
            }
        }

        if self.accepted_pairs > 0 {
            if self.accepted_pairs == 1 {
                write!(f, "#r")?;
            } else {
                write!(f, "#r{}", self.accepted_pairs)?;
            }
        }

        match self.aromatic_valence {
            AromaticValence::None => {}
            AromaticValence::Valence(1) => write!(f, "#a")?,
            AromaticValence::Valence(n) => write!(f, "#a{}", n)?,
        }

        if self.multicenter_valence > 0 {
            if self.multicenter_valence == 1 {
                write!(f, "#m")?;
            } else {
                write!(f, "#m{}", self.multicenter_valence)?;
            }
        }

        Ok(())
    }
}

impl FromStr for Atom {
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_ground_atom_dsl(s)?.to_atom()
    }
}

/// Builder for constructing `Atom` values. Used as graph node weights during
/// resolution phases and for manual molecule construction. Fields progress
/// from `None` (unknown) to `Some` as resolution phases fill them in.
#[derive(Debug, Clone)]
pub struct AtomBuilder {
    element: Element,
    isotope_mass: Option<u32>,
    charge: Option<i8>,
    implicit_hydrogens: Option<ImplicitHydrogens>,
    lone_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<SpinMultiplicity>,
    aromatic_hint: Option<bool>,
    chirality_hint: Option<Chirality>,
    candidates: SmallVec<[AtomTypeSpec; 4]>,
}

impl AtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic_hint: None,
            chirality_hint: None,
            candidates: SmallVec::new(),
        }
    }

    pub fn from_table_atom(atom: &TableAtom) -> Self {
        Self {
            element: atom.element,
            isotope_mass: atom.isotope_mass,
            charge: atom.charge,
            implicit_hydrogens: atom.implicit_hydrogens,
            lone_pairs: atom.lone_pairs,
            unpaired_electrons: atom.unpaired_electrons,
            multiplicity: atom.multiplicity,
            aromatic_hint: atom.aromatic,
            chirality_hint: atom.chirality,
            candidates: SmallVec::new(),
        }
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn isotope_mass(&self) -> Option<u32> {
        self.isotope_mass
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn hydrogen_count(&self) -> Option<u8> {
        match self.implicit_hydrogens {
            Some(ImplicitHydrogens::Hydrogens(h)) => Some(h),
            _ => None,
        }
    }

    pub fn implicit_hydrogens(&self) -> Option<ImplicitHydrogens> {
        self.implicit_hydrogens
    }

    pub fn lone_pairs(&self) -> Option<u8> {
        self.lone_pairs
    }

    pub fn unpaired_electrons(&self) -> Option<u8> {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> Option<SpinMultiplicity> {
        self.multiplicity
    }

    pub fn aromatic_hint(&self) -> Option<bool> {
        self.aromatic_hint
    }

    pub fn chirality_hint(&self) -> Option<&Chirality> {
        self.chirality_hint.as_ref()
    }

    pub fn candidates(&self) -> &[AtomTypeSpec] {
        &self.candidates
    }

    pub fn set_isotope_mass(&mut self, isotope_mass: u32) -> &mut Self {
        self.isotope_mass = Some(isotope_mass);
        self
    }

    pub fn set_charge(&mut self, charge: i8) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn set_implicit_hydrogens(&mut self, hydrogens: u8) -> &mut Self {
        self.implicit_hydrogens = Some(ImplicitHydrogens::Hydrogens(hydrogens));
        self
    }

    pub fn normal_implicit_hydrogens(&mut self) -> &mut Self {
        self.implicit_hydrogens = Some(ImplicitHydrogens::Normal);
        self
    }

    pub fn set_lone_pairs(&mut self, lone_pairs: u8) -> &mut Self {
        self.lone_pairs = Some(lone_pairs);
        self
    }

    pub fn set_unpaired_electrons(
        &mut self,
        unpaired_electrons: u8,
    ) -> Result<&mut Self, ResolutionError> {
        if let Some(m) = self.multiplicity {
            SpinState::try_new(unpaired_electrons, m)?;
        }
        self.unpaired_electrons = Some(unpaired_electrons);
        Ok(self)
    }

    pub fn set_multiplicity(
        &mut self,
        multiplicity: SpinMultiplicity,
    ) -> Result<&mut Self, ResolutionError> {
        let count = self.unpaired_electrons.unwrap_or(0);
        SpinState::try_new(count, multiplicity)?;
        self.multiplicity = Some(multiplicity);
        Ok(self)
    }

    pub fn set_aromatic_hint(&mut self, aromatic: bool) -> &mut Self {
        self.aromatic_hint = Some(aromatic);
        self
    }

    pub fn set_chirality_hint(&mut self, chirality: Chirality) -> &mut Self {
        self.chirality_hint = Some(chirality);
        self
    }

    pub fn set_candidates(&mut self, candidates: SmallVec<[AtomTypeSpec; 4]>) -> &mut Self {
        self.candidates = candidates;
        self
    }

    pub fn add_candidate(&mut self, candidate: AtomTypeSpec) -> &mut Self {
        if !self.candidates.contains(&candidate) {
            self.candidates.push(candidate);
        }
        self
    }

    fn checked_candidate(&self) -> Result<&AtomTypeSpec, ResolutionError> {
        let candidate = match self.candidates.len() {
            0 => {
                return Err(ResolutionError::ValenceNoMatch(format!(
                    "no valence match for {:?}",
                    self.element
                )))
            }
            1 => &self.candidates[0],
            n => {
                let specs: Vec<String> = self.candidates.iter().map(|s| s.to_spec_str()).collect();
                return Err(ResolutionError::ValenceAmbiguous(format!(
                    "{} valence matches for {:?}: {}",
                    n,
                    self.element,
                    specs.join(", ")
                )));
            }
        };

        if self
            .charge
            .is_some_and(|charge| charge != candidate.charge())
            || self.implicit_hydrogens.is_some_and(
                |h| matches!(h, ImplicitHydrogens::Hydrogens(n) if n != candidate.implicit_hydrogens()),
            )
            || self
                .lone_pairs
                .is_some_and(|lone_pairs| lone_pairs != candidate.lone_pairs())
            || self
                .unpaired_electrons
                .is_some_and(|u| u != candidate.unpaired_electrons())
            || self
                .multiplicity
                .is_some_and(|m| m != candidate.multiplicity())
        {
            return Err(ResolutionError::ValenceViolation(
                self.element,
                format!("atom candidate mismatch for {}", candidate.to_spec_str()),
            ));
        }

        if let Err(error) = candidate.check_invariants() {
            return Err(ResolutionError::ValenceViolation(
                self.element,
                format!(
                    "atom invariant verification failed for {}: {}",
                    candidate.to_spec_str(),
                    error
                ),
            ));
        }

        Ok(candidate)
    }

    pub fn can_build(&self) -> bool {
        self.checked_candidate().is_ok()
    }

    /// Build the final `Atom` from the builder state.
    ///
    /// Requires exactly one valence candidate remaining (resolution phases
    /// must have narrowed the set). All fields on that candidate become the
    /// definite atom properties.
    pub fn build(&self) -> Result<Atom, ResolutionError> {
        let candidate = self.checked_candidate()?;

        Ok(Atom {
            element: self.element,
            isotope_mass: self
                .isotope_mass
                .map_or(IsotopeMass::Natural, IsotopeMass::MassNumber),
            charge: candidate.charge(),
            implicit_hydrogens: candidate.implicit_hydrogens(),
            lone_pairs: candidate.lone_pairs(),
            spin: candidate.spin(),
            valence: candidate.valence(),
            donated_pairs: candidate.donated_pairs(),
            accepted_pairs: candidate.accepted_pairs(),
            aromatic_valence: candidate.aromatic_valence(),
            multicenter_valence: candidate.multicenter_valence(),
        })
    }
}

impl FromStr for AtomBuilder {
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let spec = AtomTypeSpec::from_spec_str(s)?;
        Ok(AtomBuilder {
            element: spec.element(),
            isotope_mass: None,
            charge: Some(spec.charge()),
            implicit_hydrogens: Some(ImplicitHydrogens::Hydrogens(spec.implicit_hydrogens())),
            lone_pairs: Some(spec.lone_pairs()),
            unpaired_electrons: Some(spec.unpaired_electrons()),
            multiplicity: Some(spec.multiplicity()),
            aromatic_hint: None,
            chirality_hint: None,
            candidates: SmallVec::from_elem(spec, 1),
        })
    }
}

impl From<AtomBuilder> for AtomPattern {
    fn from(value: AtomBuilder) -> Self {
        AtomPattern {
            element: Some(value.element),
            isotope_mass: value.isotope_mass.map(IsotopeMass::MassNumber),
            charge: value.charge,
            implicit_hydrogens: match value.implicit_hydrogens {
                Some(ImplicitHydrogens::Hydrogens(h)) => Some(h),
                Some(ImplicitHydrogens::Normal) | None => None,
            },
            lone_pairs: value.lone_pairs,
            unpaired_electrons: value.unpaired_electrons,
            multiplicity: value.multiplicity,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }
}

#[macro_export]
macro_rules! atom {
    ($spec:expr) => {
        $spec
            .parse::<$crate::graph_ir::atom::AtomBuilder>()
            .expect("invalid atom spec")
    };
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_data::{Element, SpinMultiplicity, SpinState};

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::element("He", Atom { element: Element::He, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 0, lone_pairs: 0,
        spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::whitespace("  He  ", Atom { element: Element::He, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 0, lone_pairs: 0,
        spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::hydrogens("C#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::isotope_natural("C#i=#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::charge_plus("C#c+#h3", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 1, implicit_hydrogens: 3, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::charge_minus("C#c-#h3#n1", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: -1, implicit_hydrogens: 3, lone_pairs: 1,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::aromatic_none("C#a!#h4", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 4, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 0, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::None, multicenter_valence: 0 })]
    #[case::no_charg("C#a#h1#v2", Atom { element: Element::C, isotope_mass: IsotopeMass::Natural, charge: 0, implicit_hydrogens: 1, lone_pairs: 0,
            spin: SpinState::new(0, SpinMultiplicity::Singlet), valence: 2, donated_pairs: 0, accepted_pairs: 0, aromatic_valence: AromaticValence::Valence(1), multicenter_valence: 0 })]
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
    #[case::unknown_element("X", AtomError::InvalidElement("X".to_string()))]
    #[case::wildcard("*", AtomError::InvalidElement("*".to_string()))]
    #[case::unknown_predicate("C#x1", AtomError::UnexpectedTag("#x".to_string()))]
    #[case::duplicate_charge("C#c+#c-", AtomError::DuplicateTag("#c".to_string()))]
    #[case::duplicate_h("C#h3#h2", AtomError::DuplicateTag("#h".to_string()))]
    #[case::non_ground_payload("C#h*", AtomError::InvalidImplicitHydrogens("*".to_string()))]
    #[case::malformed_number("C#vabc", AtomError::InvalidValence("abc".to_string()))]
    #[case::trailing_input("C foo", AtomError::UnexpectedTag("f".to_string()))]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::defaults(Atom::from_spec(AtomTypeSpec::new(Element::He, None, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0).unwrap()), "He")]
    #[case::charge_plus(Atom::from_spec(AtomTypeSpec::new(Element::C, None, 1, 3, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0).unwrap()), "C#c+#h3")]
    #[case::charge_minus(Atom::from_spec(AtomTypeSpec::new(Element::C, None, -1, 3, 1, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0).unwrap()), "C#c-#h3#n")]
    #[case::isotope_mass(Atom::from_spec(AtomTypeSpec::new(Element::C, Some(13), 0, 4, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0).unwrap()), "C#i13#h4")]
    #[case::aromatic(Atom::from_spec(AtomTypeSpec::new(Element::C, None, 0, 1, 0, 0, SpinMultiplicity::Singlet, 2, 0, 0, AromaticValence::Valence(1), 0).unwrap()), "C#h#v2#a")]
    #[case::multiplicity(Atom::from_spec(AtomTypeSpec::new(Element::C, None, 0, 0, 1, 2, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0).unwrap()), "C#n#u2#s")]
    fn test_atom_display(#[case] atom: Atom, #[case] expected: &str) {
        assert_eq!(atom.to_string(), expected);
    }
}
