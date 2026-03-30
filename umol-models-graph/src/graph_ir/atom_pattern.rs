//! Atom pattern representation and conversion into concrete ground atoms.

use umol_data::{Element, SpinMultiplicity, SpinState};

use crate::atom::{AromaticValence, IsotopeMass};
use crate::graph_ir::atom::Atom;
use crate::graph_ir::atom_type::{AtomError, AtomTypeSpec};
use crate::table_ir::atom::Atom as TableAtom;

/// Transitional atom pattern type.
///
/// All fields are optional so this can represent partial atom assignments
/// during parsing and future query/pattern workflows.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AtomPattern {
    // TODO: Check if this needs to be Option<Element> or Element
    pub element: Option<Element>,
    pub isotope_mass: Option<IsotopeMass>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub lone_pairs: Option<u8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub valence: Option<u8>,
    pub donated_pairs: Option<u8>,
    pub accepted_pairs: Option<u8>,
    pub aromatic_valence: Option<AromaticValence>,
    pub multicenter_valence: Option<u8>,
}

impl AtomPattern {
    pub fn from_atom(atom: &Atom) -> Self {
        Self {
            element: Some(atom.element()),
            isotope_mass: Some(atom.isotope_mass()),
            charge: Some(atom.charge()),
            implicit_hydrogens: Some(atom.implicit_hydrogens()),
            lone_pairs: Some(atom.lone_pairs()),
            unpaired_electrons: Some(atom.unpaired_electrons()),
            multiplicity: Some(atom.multiplicity()),
            valence: Some(atom.valence()),
            donated_pairs: Some(atom.donated_pairs()),
            accepted_pairs: Some(atom.accepted_pairs()),
            aromatic_valence: Some(atom.aromatic_valence()),
            multicenter_valence: Some(atom.multicenter_valence()),
        }
    }

    pub fn new(element: Element) -> Self {
        Self {
            element: Some(element),
            ..Self::default()
        }
    }

    pub fn from_table_atom(atom: &TableAtom) -> Self {
        Self {
            element: Some(atom.element),
            // TODO: Check SMILES spec if natural isotope is default
            isotope_mass: atom.isotope_mass.map(IsotopeMass::MassNumber),
            charge: atom.charge,
            implicit_hydrogens: match atom.implicit_hydrogens {
                Some(crate::atom::ImplicitHydrogens::Hydrogens(h)) => Some(h),
                _ => None,
            },
            lone_pairs: atom.lone_pairs,
            unpaired_electrons: atom.unpaired_electrons,
            multiplicity: atom.multiplicity,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    pub fn element(&self) -> Element {
        self.element.expect("atom pattern element must be set")
    }

    pub fn isotope_mass(&self) -> Option<u32> {
        self.isotope_mass.and_then(|m| m.mass_number())
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn hydrogen_count(&self) -> Option<u8> {
        self.implicit_hydrogens
    }

    pub fn implicit_hydrogens(&self) -> Option<crate::atom::ImplicitHydrogens> {
        self.implicit_hydrogens
            .map(crate::atom::ImplicitHydrogens::Hydrogens)
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

    pub fn check_invariants(&self) -> Result<(), AtomError> {
        let element = self
            .element
            .ok_or_else(|| AtomError::InvalidElement("<missing>".to_string()))?;
        let charge = self.charge.unwrap_or(0);
        let implicit_hydrogens = self.implicit_hydrogens.unwrap_or(0);
        let lone_pairs = self.lone_pairs.unwrap_or(0);
        let unpaired_electrons = self.unpaired_electrons.unwrap_or(0);
        let multiplicity = match self.multiplicity {
            Some(m) => m,
            None => {
                let m = unpaired_electrons.checked_add(1).ok_or_else(|| {
                    AtomError::InvalidMultiplicity((unpaired_electrons as u16 + 1).to_string())
                })?;
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or_else(|| AtomError::InvalidMultiplicity(m.to_string()))?
            }
        };
        let valence = self.valence.unwrap_or(0);
        let donated_pairs = self.donated_pairs.unwrap_or(0);
        let accepted_pairs = self.accepted_pairs.unwrap_or(0);
        let aromatic_valence = self.aromatic_valence.unwrap_or(AromaticValence::None);
        let multicenter_valence = self.multicenter_valence.unwrap_or(0);

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

    /// Convert a pattern into a concrete atom by applying ground defaults and
    /// validating spin + atom invariants.
    pub fn to_atom(&self) -> Result<Atom, AtomError> {
        self.check_invariants()?;

        let element = self
            .element
            .ok_or_else(|| AtomError::InvalidElement("<missing>".to_string()))?;
        let isotope_mass = self.isotope_mass.unwrap_or(IsotopeMass::Natural);
        let charge = self.charge.unwrap_or(0);
        let implicit_hydrogens = self.implicit_hydrogens.unwrap_or(0);
        let lone_pairs = self.lone_pairs.unwrap_or(0);
        let unpaired_electrons = self.unpaired_electrons.unwrap_or(0);
        let multiplicity = match self.multiplicity {
            Some(m) => m,
            None => {
                let m = unpaired_electrons.checked_add(1).ok_or_else(|| {
                    AtomError::InvalidMultiplicity((unpaired_electrons as u16 + 1).to_string())
                })?;
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or_else(|| AtomError::InvalidMultiplicity(m.to_string()))?
            }
        };
        let valence = self.valence.unwrap_or(0);
        let donated_pairs = self.donated_pairs.unwrap_or(0);
        let accepted_pairs = self.accepted_pairs.unwrap_or(0);
        let aromatic_valence = self.aromatic_valence.unwrap_or(AromaticValence::None);
        let multicenter_valence = self.multicenter_valence.unwrap_or(0);

        let spec = AtomTypeSpec::new(
            element,
            isotope_mass.mass_number(),
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
        )?;
        Ok(Atom::from_spec(spec))
    }

    pub fn matches_atom(&self, atom: &Atom) -> bool {
        self.element.is_none_or(|v| v == atom.element())
            && self
                .isotope_mass
                .is_none_or(|v| v == atom.isotope_mass())
            && self.charge.is_none_or(|v| v == atom.charge())
            && self
                .implicit_hydrogens
                .is_none_or(|v| v == atom.implicit_hydrogens())
            && self.lone_pairs.is_none_or(|v| v == atom.lone_pairs())
            && self
                .unpaired_electrons
                .is_none_or(|v| v == atom.unpaired_electrons())
            && self.multiplicity.is_none_or(|v| v == atom.multiplicity())
            && self.valence.is_none_or(|v| v == atom.valence())
            && self.donated_pairs.is_none_or(|v| v == atom.donated_pairs())
            && self.accepted_pairs.is_none_or(|v| v == atom.accepted_pairs())
            && self
                .aromatic_valence
                .is_none_or(|v| v == atom.aromatic_valence())
            && self
                .multicenter_valence
                .is_none_or(|v| v == atom.multicenter_valence())
    }
}

fn aromatic_increment(aromatic_valence: AromaticValence) -> u8 {
    match aromatic_valence {
        AromaticValence::None => 0,
        AromaticValence::Valence(0) => 0,
        AromaticValence::Valence(1) => 1,
        AromaticValence::Valence(2) => 0,
        AromaticValence::Valence(_) => 0,
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_data::{Element, SpinMultiplicity, SpinStateError};

    use super::*;

    #[rstest]
    #[case::defaults(AtomPattern { element: Some(Element::C), implicit_hydrogens: Some(4), ..AtomPattern::default() }, Element::C, 0, 4, 0, SpinMultiplicity::Singlet)]
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
    #[case::invalid_spin_state(AtomPattern { element: Some(Element::C), unpaired_electrons: Some(2), multiplicity: Some(SpinMultiplicity::Quartet), ..AtomPattern::default() },
        AtomError::SpinState(SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::Quartet, }))]
    #[case::invariant_mismatch(AtomPattern { element: Some(Element::O), valence: Some(2), ..AtomPattern::default() },
        AtomError::ElectronInvariantMismatch { element: Element::O, orbital_invariant: 4, electron_invariant: 8 })]
    fn test_atom_pattern_to_atom_error(#[case] pattern: AtomPattern, #[case] expected: AtomError) {
        let result = pattern.to_atom();
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), expected);
    }
}
