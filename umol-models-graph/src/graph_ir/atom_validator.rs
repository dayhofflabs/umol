//! Atom validator infrastructure copied from `umol-models-valence`.

use std::sync::LazyLock;

use super::atom::Atom;
use super::error::ResolutionError;

pub struct AtomValidator {
    #[allow(clippy::type_complexity)]
    validators: Vec<Box<dyn Fn(&Atom) -> Result<(), ResolutionError> + Send + Sync>>,
}

impl AtomValidator {
    pub fn new(
        validators: Vec<Box<dyn Fn(&Atom) -> Result<(), ResolutionError> + Send + Sync>>,
    ) -> Self {
        Self { validators }
    }

    pub fn strict() -> Self {
        Self::new(vec![Box::new(|atom| {
            let element = atom.element();
            let (min_charge, max_charge) = element.charge_bounds();
            let min_charge = i32::from(min_charge);
            let max_charge = i32::from(max_charge);

            let charge = atom.charge();
            if charge < min_charge || charge > max_charge {
                return Err(ResolutionError::ValenceViolation(
                    element,
                    format!("Charge {} is out of bounds", charge),
                ));
            }

            let unpaired = atom.unpaired_electrons();
            if unpaired > u32::from(element.max_unpaired_electrons()) {
                return Err(ResolutionError::ValenceViolation(
                    element,
                    format!("Unpaired electrons {} exceed max", unpaired),
                ));
            }

            let multiplicity = atom.multiplicity();
            if multiplicity > u32::from(element.max_unpaired_electrons()) + 1 {
                return Err(ResolutionError::ValenceViolation(
                    element,
                    format!("Multiplicity {} exceeds max", multiplicity),
                ));
            }

            let implicit_hydrogens = atom.implicit_hydrogens();
            if implicit_hydrogens > u32::from(element.max_implicit_hydrogens()) {
                return Err(ResolutionError::ValenceViolation(
                    element,
                    format!("Implicit hydrogens {} exceed max", implicit_hydrogens),
                ));
            }

            let valence = atom.valence();
            if valence > u32::from(element.max_valence()) {
                return Err(ResolutionError::ValenceViolation(
                    element,
                    format!("Valence {} exceeds max", valence),
                ));
            }

            Ok(())
        })])
    }

    pub fn lenient() -> Self {
        Self::always()
    }

    pub fn always() -> Self {
        Self::new(vec![Box::new(|_| Ok(()))])
    }

    pub fn with_validator(
        mut self,
        validator: impl Fn(&Atom) -> Result<(), ResolutionError> + Send + Sync + 'static,
    ) -> Self {
        self.validators.push(Box::new(validator));
        self
    }

    pub fn validate(&self, atom: &Atom) -> Result<(), ResolutionError> {
        for validator in &self.validators {
            validator(atom)?;
        }
        Ok(())
    }
}

impl Default for AtomValidator {
    fn default() -> Self {
        Self::strict()
    }
}

pub static DEFAULT_ATOM_VALIDATOR: LazyLock<AtomValidator> = LazyLock::new(AtomValidator::default);
pub static STRICT_ATOM_VALIDATOR: LazyLock<AtomValidator> = LazyLock::new(AtomValidator::strict);
pub static LENIENT_ATOM_VALIDATOR: LazyLock<AtomValidator> = LazyLock::new(AtomValidator::lenient);
pub static ALWAYS_ATOM_VALIDATOR: LazyLock<AtomValidator> = LazyLock::new(AtomValidator::always);
