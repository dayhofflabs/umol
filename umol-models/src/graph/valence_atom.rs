//! Valence atom implementation
//!
//! Valence atom is the node type of valence graphs and is
//! defined by its element, charge, unpaired electrons,
//! implicit hydrogens, and lone pairs.

use serde::{Deserialize, Serialize};
use std::{
    fmt::{self, Display},
    str::FromStr,
};
use umol::{Error, Result};
use umol_data::{Element, ValenceState};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ValenceAtom {
    element: Element,
    charge: i8,
    lone_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: u8,
    implicit_hydrogens: u8,
    bond_sum: u8,
}

impl ValenceAtom {
    pub fn new(
        element: Element,
        charge: i8,
        lone_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: u8,
        implicit_hydrogens: u8,
        bond_sum: u8,
    ) -> Self {
        ValenceAtom {
            element,
            charge,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            implicit_hydrogens,
            bond_sum,
        }
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn lone_pairs(&self) -> u8 {
        self.lone_pairs
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> u8 {
        self.multiplicity
    }

    pub fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }

    pub fn bond_sum(&self) -> u8 {
        self.bond_sum
    }

    pub fn valence(&self) -> u8 {
        self.bond_sum + self.implicit_hydrogens
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(element, 0, 0, 0, 1, 0, 0)
    }
}

pub struct ValenceAtomBuilder {
    element: Element,
    charge: i8,
    lone_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: u8,
    implicit_hydrogens: u8,
    bond_sum: u8,
}

impl ValenceAtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: 0,
            lone_pairs: 0,
            unpaired_electrons: 0,
            multiplicity: 1,
            implicit_hydrogens: 0,
            bond_sum: 0,
        }
    }

    pub fn charge(mut self, charge: i8) -> Self {
        self.charge = charge;
        self
    }

    pub fn lone_pairs(mut self, lone_pairs: u8) -> Self {
        self.lone_pairs = lone_pairs;
        self
    }

    pub fn unpaired_electrons(mut self, unpaired_electrons: u8) -> Self {
        self.unpaired_electrons = unpaired_electrons;
        self
    }

    pub fn multiplicity(mut self, multiplicity: u8) -> Self {
        self.multiplicity = multiplicity;
        self
    }

    pub fn implicit_hydrogens(mut self, implicit_hydrogens: u8) -> Self {
        self.implicit_hydrogens = implicit_hydrogens;
        self
    }

    pub fn bond_sum(mut self, bond_sum: u8) -> Self {
        self.bond_sum = bond_sum;
        self
    }

    pub fn build(self) -> Result<ValenceAtom> {
        Ok(ValenceAtom {
            element: self.element,
            charge: self.charge,
            lone_pairs: self.lone_pairs,
            unpaired_electrons: self.unpaired_electrons,
            multiplicity: self.multiplicity,
            implicit_hydrogens: self.implicit_hydrogens,
            bond_sum: self.bond_sum,
        })
    }
}

impl From<Element> for ValenceAtom {
    fn from(element: Element) -> Self {
        ValenceAtom::from_element(element)
    }
}

impl From<ValenceState> for ValenceAtom {
    fn from(state: ValenceState) -> Self {
        ValenceAtom::new(
            state.element(),
            state.charge(),
            state.lone_pairs(),
            state.unpaired_electrons(),
            state.multiplicity(),
            0,
            state.valence(),
        )
    }
}

impl TryFrom<&str> for ValenceAtom {
    type Error = Error;

    fn try_from(symbol: &str) -> Result<Self> {
        Element::try_from(symbol).map(|element| element.into())
    }
}

impl FromStr for ValenceAtom {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        ValenceAtom::try_from(s)
    }
}

impl Display for ValenceAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let symbol = self.element.symbol();

        // Use SMILES-compatible notation (except for unpaired electrons)
        // Always use brackets if:
        // - Symbol is 2 characters
        // - Has charge
        // - Has unpaired electrons
        let needs_brackets = symbol.len() == 2 || self.charge != 0 || self.unpaired_electrons > 0;

        if needs_brackets {
            write!(f, "[{}", symbol)?;

            if self.charge > 0 {
                write!(f, "+{}", self.charge)?;
            } else if self.charge < 0 {
                write!(f, "{}", self.charge)?; // Negative sign included automatically
            }

            if self.unpaired_electrons > 0 {
                write!(f, "^{}", self.unpaired_electrons)?;
            }

            write!(f, "]")
        } else {
            write!(f, "{}", symbol)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use umol_data::Element;

    #[test]
    fn test_valence_atom_new() {
        let atom = ValenceAtom::new(Element::C, 0, 1, 2, 3, 0, 4);
        assert_eq!(atom.element(), Element::C);
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.multiplicity(), 3);
        assert_eq!(atom.bond_sum(), 4);
    }

    #[test]
    fn test_valence_atom_from_valence_state() {
        let state = ValenceState::new(Element::C, 0, 1, 2, 3, 4);
        let atom = ValenceAtom::from(state);
        assert_eq!(atom, ValenceAtom::new(Element::C, 0, 1, 2, 3, 0, 4));
    }

    #[test]
    fn test_valence_atom_from_element() {
        let atom: ValenceAtom = Element::O.into();
        assert_eq!(atom.element(), Element::O);
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
    }

    #[test]
    fn test_valence_atom_builder() {
        let atom = ValenceAtomBuilder::new(Element::N)
            .charge(-1)
            .lone_pairs(1)
            .unpaired_electrons(2)
            .multiplicity(3)
            .implicit_hydrogens(0)
            .bond_sum(1)
            .build()
            .unwrap();

        assert_eq!(atom.element(), Element::N);
        assert_eq!(atom.charge(), -1);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.multiplicity(), 3);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.bond_sum(), 1);
    }

    #[test]
    fn test_valence_atom_builder_validation() {
        // Valid combination
        let result = ValenceAtomBuilder::new(Element::C)
            .charge(-1)
            .unpaired_electrons(3)
            .build();
        assert!(result.is_ok());
    }

    #[test]
    fn test_valence_atom_try_from_str() {
        let atom = ValenceAtom::try_from("C").unwrap();
        assert_eq!(atom, ValenceAtom::from_element(Element::C));

        let atom = ValenceAtom::try_from("N").unwrap();
        assert_eq!(atom, ValenceAtom::from_element(Element::N));

        assert!(ValenceAtom::try_from("Xx").is_err());
    }

    #[rstest]
    #[case(Element::C, 0, 0, "C")]
    #[case(Element::C, 1, 0, "[C+1]")]
    #[case(Element::C, -1, 0, "[C-1]")]
    #[case(Element::C, 0, 2, "[C^2]")]
    #[case(Element::C, 1, 2, "[C+1^2]")]
    #[case(Element::He, 0, 0, "[He]")] // Two-letter element always gets brackets
    fn test_valence_atom_display(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] unpaired: u8,
        #[case] expected: &str,
    ) {
        let atom = ValenceAtom::new(element, charge, 0, unpaired, 1, 0, 0);
        assert_eq!(format!("{}", atom), expected);
    }
}
