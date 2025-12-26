//! Atom validator infrastructure copied from `umol-models-valence`.

use once_cell::sync::Lazy;
use super::atom::Atom;
use super::error::GraphError;

type Result<T> = std::result::Result<T, GraphError>;

pub struct AtomValidator {
    #[allow(clippy::type_complexity)]
    validators: Vec<Box<dyn Fn(&Atom) -> Result<()> + Send + Sync>>,
}

impl AtomValidator {
    pub fn new(validators: Vec<Box<dyn Fn(&Atom) -> Result<()> + Send + Sync>>) -> Self {
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
                return Err(GraphError::ValenceViolation(
                    element,
                    format!("Charge {} is out of bounds", charge),
                ));
            }

            let unpaired = atom.unpaired_e();
            if unpaired > u32::from(element.max_unpaired_electrons()) {
                return Err(GraphError::ValenceViolation(
                    element,
                    format!("Unpaired electrons {} exceed max", unpaired),
                ));
            }

            let multiplicity = atom.multiplicity();
            if multiplicity > u32::from(element.max_unpaired_electrons()) + 1 {
                return Err(GraphError::ValenceViolation(
                    element,
                    format!("Multiplicity {} exceeds max", multiplicity),
                ));
            }

            let implicit_h = atom.implicit_h();
            if implicit_h > u32::from(element.max_implicit_hydrogens()) {
                return Err(GraphError::ValenceViolation(
                    element,
                    format!("Implicit hydrogens {} exceed max", implicit_h),
                ));
            }

            let valence = atom.valence();
            if valence > u32::from(element.max_valence()) {
                return Err(GraphError::ValenceViolation(
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
        validator: impl Fn(&Atom) -> Result<()> + Send + Sync + 'static,
    ) -> Self {
        self.validators.push(Box::new(validator));
        self
    }

    pub fn validate(&self, atom: &Atom) -> Result<()> {
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

pub static DEFAULT_ATOM_VALIDATOR: Lazy<AtomValidator> = Lazy::new(AtomValidator::default);
pub static STRICT_ATOM_VALIDATOR: Lazy<AtomValidator> = Lazy::new(AtomValidator::strict);
pub static LENIENT_ATOM_VALIDATOR: Lazy<AtomValidator> = Lazy::new(AtomValidator::lenient);
pub static ALWAYS_ATOM_VALIDATOR: Lazy<AtomValidator> = Lazy::new(AtomValidator::always);
