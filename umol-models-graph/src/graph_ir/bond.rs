//! Bond types for GraphIR.

use std::str::FromStr;

use thiserror::Error;
use umol_data::{SpinMultiplicity, SpinState, SpinStateError};

use crate::graph_ir::error::ResolutionError;
use crate::table_ir::bond::{Bond as TableBond, BondOrder};

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum BondError {
    #[error("bond spec must use b{{...}} notation")]
    InvalidFormat,
    #[error("bond spec is empty")]
    EmptySpec,
    #[error("bond order token 'o<n>' is required")]
    MissingOrder,
    #[error("duplicate token '{token}'")]
    DuplicateToken { token: char },
    #[error("invalid bond order token")]
    InvalidOrder,
    #[error("invalid charge token")]
    InvalidCharge,
    #[error("invalid multiplicity token")]
    InvalidMultiplicity,
    #[error("invalid aromatic hint token; expected a0 or a1")]
    InvalidAromaticHint,
    #[error("unexpected token '{token}'")]
    UnexpectedToken { token: char },
    #[error(transparent)]
    SpinState(#[from] SpinStateError),
}

impl From<BondError> for ResolutionError {
    fn from(value: BondError) -> Self {
        ResolutionError::InvalidBond(value.to_string())
    }
}

/// Resolved shared (covalent) bond in GraphIR. Order is the localized (σ-skeleton)
/// bond order. Dative and non-covalent bonds are stored separately.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
    charge: i8,
    spin: SpinState,
}

impl Bond {
    pub fn new(order: u8) -> Self {
        Self {
            order,
            charge: 0,
            spin: SpinState::closed_shell(),
        }
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn spin(&self) -> SpinState {
        self.spin
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.spin.multiplicity()
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.spin.unpaired_electrons()
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn to_builder(&self) -> BondBuilder {
        BondBuilder {
            order: self.order,
            charge: Some(self.charge),
            unpaired_electrons: Some(self.spin.unpaired_electrons()),
            multiplicity: Some(self.spin.multiplicity()),
            aromatic_hint: None,
        }
    }
}

/// Mutable bond representation used during resolution phases.
/// Carries an aromaticity hint that is consumed by Kekulization;
/// `build()` produces the final `Bond` with a definite order.
#[derive(Debug, Clone, PartialEq)]
pub struct BondBuilder {
    order: u8,
    charge: Option<i8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<SpinMultiplicity>,
    aromatic_hint: Option<bool>,
}

impl BondBuilder {
    pub fn new(order: u8, aromatic_hint: Option<bool>) -> Self {
        Self {
            order,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic_hint,
        }
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn multiplicity(&self) -> Option<SpinMultiplicity> {
        self.multiplicity
    }

    pub fn unpaired_electrons(&self) -> Option<u8> {
        self.unpaired_electrons
    }

    pub fn aromatic_hint(&self) -> Option<bool> {
        self.aromatic_hint
    }

    pub fn set_order(&mut self, order: u8) {
        self.order = order;
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = Some(charge);
    }

    pub fn set_unpaired_electrons(&mut self, unpaired_electrons: u8) {
        self.unpaired_electrons = Some(unpaired_electrons);
    }

    pub fn set_multiplicity(&mut self, multiplicity: SpinMultiplicity) {
        self.multiplicity = Some(multiplicity);
    }

    pub fn set_aromatic_hint(&mut self, aromatic: Option<bool>) {
        self.aromatic_hint = aromatic;
    }

    pub fn from_table_bond(bond: &TableBond) -> Self {
        debug_assert!(
            !bond.order.is_query(),
            "query bond orders must be resolved before conversion to BondBuilder"
        );
        let (order, aromatic_hint) = match bond.order {
            BondOrder::Aromatic => (1, Some(true)),
            o => (
                o.value()
                    .expect("non-query, non-aromatic bond order must have a value"),
                Some(false),
            ),
        };
        Self {
            order,
            charge: bond.charge,
            unpaired_electrons: bond.unpaired_electrons,
            multiplicity: bond.multiplicity,
            aromatic_hint,
        }
    }

    fn checked_spin_and_charge(&self) -> Result<(i8, SpinState), ResolutionError> {
        let charge = self.charge.unwrap_or(0);
        let electrons: i16 = 2 * self.order as i16 - charge as i16;
        if electrons < 0 {
            return Err(ResolutionError::BondInvariantViolation(format!(
                "bond state is not buildable: order={}, charge={:?}, unpaired_electrons={:?}, multiplicity={:?}",
                self.order, self.charge, self.unpaired_electrons, self.multiplicity
            )));
        }
        let spin = match (self.unpaired_electrons, self.multiplicity) {
            (Some(unpaired), Some(multiplicity)) => SpinState::try_new(unpaired, multiplicity),
            (Some(unpaired), None) => SpinState::max_multiplicity(unpaired).ok_or(
                SpinStateError::UnpairedElectronsExceedMax {
                    unpaired_electrons: unpaired,
                },
            ),
            (None, Some(multiplicity)) => {
                SpinState::try_new(multiplicity.multiplicity() - 1, multiplicity)
            }
            (None, None) => Ok(SpinState::closed_shell()),
        }
        .map_err(|_e| {
            ResolutionError::BondInvariantViolation(format!(
                "bond state is not buildable: order={}, charge={:?}, unpaired_electrons={:?}, multiplicity={:?}",
                self.order, self.charge, self.unpaired_electrons, self.multiplicity
            ))
        })?;
        if !spin.is_compatible_with(electrons as u8) {
            return Err(ResolutionError::BondInvariantViolation(format!(
                "bond state is not buildable: order={}, charge={:?}, unpaired_electrons={:?}, multiplicity={:?}",
                self.order, self.charge, self.unpaired_electrons, self.multiplicity
            )));
        }
        Ok((charge, spin))
    }

    pub fn can_build(&self) -> bool {
        self.checked_spin_and_charge().is_ok()
    }

    pub fn build(&self) -> Result<Bond, ResolutionError> {
        let (charge, spin) = self.checked_spin_and_charge()?;
        Ok(Bond {
            order: self.order,
            charge,
            spin,
        })
    }
}

impl FromStr for BondBuilder {
    type Err = BondError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        let body = trimmed
            .strip_prefix("b{")
            .and_then(|rest| rest.strip_suffix('}'))
            .ok_or(BondError::InvalidFormat)?;
        if body.trim().is_empty() {
            return Err(BondError::EmptySpec);
        }

        let mut chars = body.chars().peekable();

        let mut order: Option<u8> = None;
        let mut charge: Option<i8> = None;
        let mut unpaired_electrons: Option<u8> = None;
        let mut multiplicity: Option<SpinMultiplicity> = None;
        let mut aromatic_hint: Option<bool> = None;

        while let Some(token) = chars.next() {
            if token.is_ascii_whitespace() {
                continue;
            }

            let mut number = String::new();
            while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
                number.push(chars.next().expect("peeked digit must exist"));
            }

            match token {
                'o' => {
                    if order.is_some() {
                        return Err(BondError::DuplicateToken { token });
                    }
                    if number.is_empty() {
                        return Err(BondError::InvalidOrder);
                    }
                    order = Some(number.parse::<u8>().map_err(|_| BondError::InvalidOrder)?);
                }
                '+' | '-' => {
                    if charge.is_some() {
                        return Err(BondError::DuplicateToken { token });
                    }
                    let magnitude = if number.is_empty() {
                        1
                    } else {
                        number.parse::<i8>().map_err(|_| BondError::InvalidCharge)?
                    };
                    charge = Some(if token == '-' { -magnitude } else { magnitude });
                }
                '^' => {
                    if unpaired_electrons.is_some() {
                        return Err(BondError::DuplicateToken { token });
                    }
                    let n = if number.is_empty() {
                        1
                    } else {
                        number
                            .parse::<u8>()
                            .map_err(|_| BondError::InvalidMultiplicity)?
                    };
                    unpaired_electrons = Some(n);
                }
                'x' => {
                    if multiplicity.is_some() {
                        return Err(BondError::DuplicateToken { token });
                    }
                    let m = if number.is_empty() {
                        1
                    } else {
                        number
                            .parse::<u8>()
                            .map_err(|_| BondError::InvalidMultiplicity)?
                    };
                    multiplicity = Some(
                        SpinMultiplicity::from_multiplicity(m)
                            .ok_or(BondError::InvalidMultiplicity)?,
                    );
                }
                'a' => {
                    if aromatic_hint.is_some() {
                        return Err(BondError::DuplicateToken { token });
                    }
                    let hint = match number.as_str() {
                        "0" => false,
                        "1" => true,
                        _ => return Err(BondError::InvalidAromaticHint),
                    };
                    aromatic_hint = Some(hint);
                }
                _ => return Err(BondError::UnexpectedToken { token }),
            }
        }

        let order = order.ok_or(BondError::MissingOrder)?;
        if let (Some(unpaired), Some(mult)) = (unpaired_electrons, multiplicity) {
            SpinState::try_new(unpaired, mult)?;
        }

        let mut builder = BondBuilder::new(order, aromatic_hint);
        if let Some(charge) = charge {
            builder.set_charge(charge);
        }
        if let Some(unpaired) = unpaired_electrons {
            builder.set_unpaired_electrons(unpaired);
        }
        if let Some(mult) = multiplicity {
            builder.set_multiplicity(mult);
        }
        Ok(builder)
    }
}

#[macro_export]
macro_rules! bond {
    ($spec:expr) => {
        $spec
            .parse::<$crate::graph_ir::bond::BondBuilder>()
            .expect("invalid bond spec")
    };
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::{SpinMultiplicity, SpinStateError};

    use super::*;

    fn bond_builder_with_state(
        order: u8,
        charge: Option<i8>,
        unpaired_electrons: Option<u8>,
        multiplicity: Option<SpinMultiplicity>,
    ) -> BondBuilder {
        let mut builder = BondBuilder::new(order, None);
        if let Some(charge) = charge {
            builder.set_charge(charge);
        }
        if let Some(unpaired_electrons) = unpaired_electrons {
            builder.set_unpaired_electrons(unpaired_electrons);
        }
        if let Some(multiplicity) = multiplicity {
            builder.set_multiplicity(multiplicity);
        }
        builder
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell(1, None, None, None, true)]
    #[case::high_spin(1, Some(1), Some(1), None, true)]
    #[case::from_multiplicity(2, Some(0), None, Some(SpinMultiplicity::Triplet), true)]
    #[case::negative_electrons(0, Some(1), None, None, false)]
    #[case::incompatible_spin_pair(1, None, Some(0), Some(SpinMultiplicity::Triplet), false)]
    #[case::electron_parity_mismatch(1, Some(0), Some(1), None, false)]
    #[case::max_unpaired_exceeded(1, Some(0), Some(10), None, false)]
    fn test_bond_builder_can_build(
        #[case] order: u8,
        #[case] charge: Option<i8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
        #[case] expected: bool,
    ) {
        let builder = bond_builder_with_state(order, charge, unpaired_electrons, multiplicity);
        assert_eq!(builder.can_build(), expected);
        assert_eq!(builder.build().is_ok(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell(1, None, None, None, 0, 0, SpinMultiplicity::Singlet)]
    #[case::high_spin(1, Some(1), Some(1), None, 1, 1, SpinMultiplicity::Doublet)]
    #[case::from_multiplicity(2, Some(0), None, Some(SpinMultiplicity::Triplet), 0, 2, SpinMultiplicity::Triplet)]
    #[case::complete(1, Some(0), Some(2), Some(SpinMultiplicity::Singlet), 0, 2, SpinMultiplicity::Singlet)]
    fn test_bond_builder_build(
        #[case] order: u8,
        #[case] charge: Option<i8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
        #[case] expected_charge: i8,
        #[case] expected_unpaired: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let builder = bond_builder_with_state(order, charge, unpaired_electrons, multiplicity);
        let bond = builder.build().expect("expected build success");
        assert_eq!(bond.order(), order);
        assert_eq!(bond.charge(), expected_charge);
        assert_eq!(bond.unpaired_electrons(), expected_unpaired);
        assert_eq!(bond.multiplicity(), expected_multiplicity);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(1, BondOrder::Aromatic, Some(1), Some(1), Some(SpinMultiplicity::Doublet), Some(true))]
    #[case::single(1, BondOrder::Single, None, None, None, Some(false))]
    #[case::double(2, BondOrder::Double, Some(-1), Some(1), Some(SpinMultiplicity::Doublet), Some(false))]
    fn test_bond_builder_from_table_bond(
        #[case] expected_order: u8,
        #[case] order: BondOrder,
        #[case] charge: Option<i8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
        #[case] expected_aromatic_hint: Option<bool>,
    ) {
        let mut bond = TableBond::new(0, 1, order);
        bond.charge = charge;
        bond.unpaired_electrons = unpaired_electrons;
        bond.multiplicity = multiplicity;

        let builder = BondBuilder::from_table_bond(&bond);
        assert_eq!(builder.order(), expected_order);
        assert_eq!(builder.charge(), charge);
        assert_eq!(builder.unpaired_electrons(), unpaired_electrons);
        assert_eq!(builder.multiplicity(), multiplicity);
        assert_eq!(builder.aromatic_hint(), expected_aromatic_hint);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order_only("b{o1}", 1, None, None, None, None)]
    #[case::full("b{o2+1^1x2a1}", 2, Some(1), Some(1), Some(SpinMultiplicity::Doublet), Some(true))]
    #[case::spaced(" b{ o1 -2 ^2 x1 a0 } ", 1, Some(-2), Some(2), Some(SpinMultiplicity::Singlet), Some(false))]
    fn test_bond_builder_from_str(
        #[case] input: &str,
        #[case] expected_order: u8,
        #[case] expected_charge: Option<i8>,
        #[case] expected_unpaired_electrons: Option<u8>,
        #[case] expected_multiplicity: Option<SpinMultiplicity>,
        #[case] expected_aromatic_hint: Option<bool>,
    ) {
        let builder: BondBuilder = input.parse().expect("expected parse success");
        assert_eq!(builder.order(), expected_order);
        assert_eq!(builder.charge(), expected_charge);
        assert_eq!(builder.unpaired_electrons(), expected_unpaired_electrons);
        assert_eq!(builder.multiplicity(), expected_multiplicity);
        assert_eq!(builder.aromatic_hint(), expected_aromatic_hint);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::missing_tag("{o1}", BondError::InvalidFormat)]
    #[case::empty("b{}", BondError::EmptySpec)]
    #[case::missing_order("b{+1}", BondError::MissingOrder)]
    #[case::duplicate_order("b{o1o2}", BondError::DuplicateToken { token: 'o' })]
    #[case::invalid_order("b{o}", BondError::InvalidOrder)]
    #[case::invalid_aromatic("b{o1a2}", BondError::InvalidAromaticHint)]
    #[case::unexpected_token("b{o1q2}", BondError::UnexpectedToken { token: 'q' })]
    #[case::invalid_spin_pair("b{o1^0x3}", BondError::SpinState(SpinStateError::Incompatible {
        unpaired_electrons: 0,
        multiplicity: SpinMultiplicity::Triplet,
    }))]
    fn test_bond_builder_from_str_error(#[case] input: &str, #[case] expected: BondError) {
        assert_eq!(input.parse::<BondBuilder>().unwrap_err(), expected);
    }
}
