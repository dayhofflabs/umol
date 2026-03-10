//! Bond types for GraphIR.

use umol_data::SpinMultiplicity;

use crate::table_ir::bond::{Bond as TableBond, BondOrder};

/// Resolved shared (covalent) bond in GraphIR. Order is the localized (σ-skeleton)
/// bond order. Dative and non-covalent bonds are stored separately.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
    charge: i8,
    multiplicity: SpinMultiplicity,
}

impl Bond {
    pub fn new(order: u8) -> Self {
        Self {
            order,
            charge: 0,
            multiplicity: SpinMultiplicity::Singlet,
        }
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn to_builder(&self) -> BondBuilder {
        BondBuilder {
            order: self.order,
            charge: Some(self.charge),
            multiplicity: Some(self.multiplicity),
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
    multiplicity: Option<SpinMultiplicity>,
    aromatic_hint: Option<bool>,
}

impl BondBuilder {
    pub fn new(order: u8, aromatic_hint: Option<bool>) -> Self {
        Self {
            order,
            charge: None,
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

    pub fn aromatic_hint(&self) -> Option<bool> {
        self.aromatic_hint
    }

    pub fn set_order(&mut self, order: u8) {
        self.order = order;
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = Some(charge);
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
            multiplicity: bond.multiplicity,
            aromatic_hint,
        }
    }

    pub fn build(&self) -> Bond {
        Bond {
            order: self.order,
            charge: self.charge.unwrap_or(0),
            multiplicity: self.multiplicity.unwrap_or(SpinMultiplicity::Singlet),
        }
    }
}
