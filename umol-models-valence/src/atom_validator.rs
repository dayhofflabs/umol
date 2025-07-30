//! Atom validators

use crate::AtomBuilder;
use once_cell::sync::Lazy;
use umol::error::DataError;
use umol::Result;

/// Validators for `Atom` type. Checks atom properties against element-specific constraints.
/// The `AtomValidator` type is a collection of validators that are applied to an `AtomBuilder` type.
/// Default validator checks charge, lone pairs, donated pairs, unpaired electrons, multiplicity,
/// implicit hydrogens, and valence against element-specific constraints.
/// Interrelations between properties are not validated here, they are checked as part of the
/// `AtomType` matching.
pub struct AtomValidator {
    #[allow(clippy::type_complexity)]
    validators: Vec<Box<dyn Fn(&AtomBuilder) -> Result<()> + Send + Sync>>,
}

impl AtomValidator {
    pub fn new(validators: Vec<Box<dyn Fn(&AtomBuilder) -> Result<()> + Send + Sync>>) -> Self {
        Self { validators }
    }

    pub fn strict() -> Self {
        Self::new(vec![Box::new(|builder| {
            let element = builder.element();
            let (min_charge, max_charge) = element.charge_bounds();
            if let Some(charge) = builder.charge() {
                if charge < min_charge || charge > max_charge {
                    return Err(DataError::InvalidAtomCharge(format!(
                        "Charge {} is out of bounds for element {}",
                        charge, element
                    ))
                    .into());
                }
            }
            if let Some(unpaired_electrons) = builder.unpaired_electrons() {
                if unpaired_electrons > element.max_unpaired_electrons() {
                    return Err(DataError::InvalidAtomUnpairedElectrons(format!(
                        "Unpaired electrons {} exceed max for element {}",
                        unpaired_electrons, element
                    ))
                    .into());
                }
            }
            if let Some(multiplicity) = builder.multiplicity() {
                if multiplicity > element.max_unpaired_electrons() + 1 {
                    return Err(DataError::InvalidAtomMultiplicity(format!(
                        "Multiplicity {} exceeds max for element {}",
                        multiplicity, element
                    ))
                    .into());
                }
            }
            if let Some(implicit_hydrogens) = builder.implicit_hydrogens() {
                if implicit_hydrogens > element.max_implicit_hydrogens() {
                    return Err(DataError::InvalidAtomImplicitHydrogens(format!(
                        "Implicit hydrogens {} exceed max for element {}",
                        implicit_hydrogens, element
                    ))
                    .into());
                }
            }
            if let Some(valence) = builder.valence() {
                if valence > element.max_valence() {
                    return Err(DataError::InvalidAtomValence(format!(
                        "Valence {} exceeds max for element {}",
                        valence, element
                    ))
                    .into());
                }
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
        validator: impl Fn(&AtomBuilder) -> Result<()> + Send + Sync + 'static,
    ) -> Self {
        self.validators.push(Box::new(validator));
        self
    }

    pub fn validate(&self, builder: &AtomBuilder) -> Result<()> {
        for validator in &self.validators {
            validator(builder)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use umol_data::{e, Element};

    #[test]
    fn test_atom_validator() {
        let validator = AtomValidator::default();
        let builder = AtomBuilder::new(Element::H);
        let result = validator.validate(&builder);
        assert!(result.is_ok());
    }

    #[test]
    fn test_atom_validator_invalid_charge() {
        let validator = AtomValidator::default();
        let mut builder = AtomBuilder::new(Element::H);
        builder.set_charge(2);
        let result = validator.validate(&builder);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            DataError::InvalidAtomCharge(format!(
                "Charge {} is out of bounds for element {}",
                2,
                Element::H
            ))
            .to_string()
        );
    }

    #[test]
    fn test_atom_validator_invalid_unpaired_electrons() {
        let validator = AtomValidator::default();
        let mut builder = AtomBuilder::new(Element::C);
        builder.set_unpaired_electrons(5);
        let result = validator.validate(&builder);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            DataError::InvalidAtomUnpairedElectrons(format!(
                "Unpaired electrons {} exceed max for element {}",
                5,
                Element::C
            ))
            .to_string()
        );
    }

    #[test]
    fn test_atom_validator_invalid_multiplicity() {
        let validator = AtomValidator::default();
        let mut builder = AtomBuilder::new(Element::C);
        builder.set_multiplicity(6);
        let result = validator.validate(&builder);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            DataError::InvalidAtomMultiplicity(format!(
                "Multiplicity {} exceeds max for element {}",
                6,
                Element::C
            ))
            .to_string()
        );
    }

    #[test]
    fn test_atom_validator_invalid_implicit_hydrogens() {
        let validator = AtomValidator::default();
        let mut builder = AtomBuilder::new(Element::C);
        builder.set_implicit_hydrogens(5);
        let result = validator.validate(&builder);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            DataError::InvalidAtomImplicitHydrogens(format!(
                "Implicit hydrogens {} exceed max for element {}",
                5,
                Element::C
            ))
            .to_string()
        );
    }

    #[test]
    fn test_atom_validator_invalid_valence() {
        let validator = AtomValidator::default();
        let mut builder = AtomBuilder::new(Element::C);
        builder.set_valence(10);
        let result = validator.validate(&builder);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            DataError::InvalidAtomValence(format!(
                "Valence {} exceeds max for element {}",
                10,
                Element::C
            ))
            .to_string()
        );
    }

    #[test]
    fn test_atom_validator_custom() {
        let validator = AtomValidator::default().with_validator(|builder| {
            if builder.element() == Element::C {
                Ok(())
            } else {
                Err(DataError::InvalidElement(format!(
                    "Element {} is not valid for this validator",
                    builder.element()
                ))
                .into())
            }
        });
        let builder = AtomBuilder::new(e!(C));
        let result = validator.validate(&builder);
        assert!(result.is_ok());
        let builder = AtomBuilder::new(e!(H));
        let result = validator.validate(&builder);
        assert!(result.is_err());
    }

    #[test]
    fn test_atom_validator_always() {
        let validator = AtomValidator::always();
        let builder = AtomBuilder::new(e!(C));
        let result = validator.validate(&builder);
        assert!(result.is_ok());
    }

    #[test]
    fn test_default_atom_validator_lazy_static() {
        let validator = DEFAULT_ATOM_VALIDATOR.validate(&AtomBuilder::new(e!(C)));
        assert!(validator.is_ok());
    }

    #[test]
    fn test_always_atom_validator_lazy_static() {
        let validator = ALWAYS_ATOM_VALIDATOR.validate(&AtomBuilder::new(e!(C)));
        assert!(validator.is_ok());
    }
}
