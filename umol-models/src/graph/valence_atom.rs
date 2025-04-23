//! Valence atom implementation
//!
//! Valence atom is the node type of valence graphs and is
//! defined by its element, charge, unpaired electrons,
//! implicit hydrogens, and lone pairs.

use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use umol::Result;
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
}

impl Display for ValenceAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        ValenceState::new(
            self.element,
            self.charge,
            self.lone_pairs,
            self.unpaired_electrons,
            self.multiplicity,
            self.implicit_hydrogens,
        )
        .fmt(f)
    }
}

#[derive(Debug)]
pub struct ValenceAtomBuilder {
    element: Option<Element>,
    charge: Option<i8>,
    lone_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<u8>,
    implicit_hydrogens: Option<u8>,
    bond_sum: Option<u8>,
}

impl ValenceAtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element: Some(element),
            charge: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            implicit_hydrogens: None,
            bond_sum: None,
        }
    }

    pub fn element(&mut self, element: Element) -> &mut Self {
        self.element = Some(element);
        self
    }

    pub fn charge(&mut self, charge: i8) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn lone_pairs(&mut self, lone_pairs: u8) -> &mut Self {
        self.lone_pairs = Some(lone_pairs);
        self
    }

    pub fn unpaired_electrons(&mut self, unpaired_electrons: u8) -> &mut Self {
        self.unpaired_electrons = Some(unpaired_electrons);
        self
    }

    pub fn multiplicity(&mut self, multiplicity: u8) -> &mut Self {
        self.multiplicity = Some(multiplicity);
        self
    }

    pub fn implicit_hydrogens(&mut self, implicit_hydrogens: u8) -> &mut Self {
        self.implicit_hydrogens = Some(implicit_hydrogens);
        self
    }

    pub fn bond_sum(&mut self, bond_sum: u8) -> &mut Self {
        self.bond_sum = Some(bond_sum);
        self
    }

    pub fn build(self) -> Result<ValenceAtom> {
        self.validate()?;
        self.infer_type()
    }

    fn validate(&self) -> Result<()> {
        Ok(())
    }

    fn infer_type(self) -> Result<ValenceAtom> {
        Ok(ValenceAtom {
            element: self.element.unwrap(),
            charge: self.charge.unwrap_or(0),
            lone_pairs: self.lone_pairs.unwrap_or(0),
            unpaired_electrons: self.unpaired_electrons.unwrap_or(0),
            multiplicity: self.multiplicity.unwrap_or(0),
            implicit_hydrogens: self.implicit_hydrogens.unwrap_or(0),
            bond_sum: self.bond_sum.unwrap_or(0),
        })
    }
}

impl From<ValenceState> for ValenceAtomBuilder {
    fn from(state: ValenceState) -> Self {
        ValenceAtomBuilder {
            element: Some(state.element()),
            charge: Some(state.charge()),
            lone_pairs: Some(state.lone_pairs()),
            unpaired_electrons: Some(state.unpaired_electrons()),
            multiplicity: Some(state.multiplicity()),
            implicit_hydrogens: None,
            bond_sum: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use umol_data::{Element, ValenceState};

    #[test]
    fn test_valence_atom_display() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder
            .charge(0)
            .lone_pairs(1)
            .unpaired_electrons(2)
            .multiplicity(3)
            .implicit_hydrogens(0)
            .bond_sum(0);

        let atom = builder.build().unwrap();
        assert_eq!(format!("{}", atom), "[C/1^2]");
    }

    #[test]
    fn test_valence_atom_serialize() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder
            .charge(0)
            .lone_pairs(1)
            .unpaired_electrons(2)
            .multiplicity(3)
            .implicit_hydrogens(0)
            .bond_sum(0);

        let atom = builder.build().unwrap();
        let serialized = serde_json::to_string(&atom).unwrap();
        assert_eq!(serialized, "{\"element\":\"C\",\"charge\":0,\"lone_pairs\":1,\"unpaired_electrons\":2,\"multiplicity\":3,\"implicit_hydrogens\":0,\"bond_sum\":0}");

        let deserialized: ValenceAtom = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, atom);
    }

    #[test]
    fn test_valence_atom_builder() {
        let mut builder = ValenceAtomBuilder::new(Element::N);

        builder
            .charge(-1)
            .lone_pairs(1)
            .unpaired_electrons(2)
            .multiplicity(3)
            .implicit_hydrogens(0)
            .bond_sum(1);

        let atom = builder.build().unwrap();

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
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder.charge(-1);
        builder.unpaired_electrons(3);
        builder.implicit_hydrogens(0);
        let result = builder.build();
        assert!(result.is_ok());
    }

    #[test]
    fn test_valence_atom_builder_from_valence_state() {
        let atom: ValenceAtomBuilder = ValenceState::new(Element::C, 0, 1, 2, 3, 4).into();
        assert_eq!(atom.element, Some(Element::C));
        assert_eq!(atom.charge, Some(0));
        assert_eq!(atom.lone_pairs, Some(1));
        assert_eq!(atom.unpaired_electrons, Some(2));
        assert_eq!(atom.multiplicity, Some(3));
        assert_eq!(atom.implicit_hydrogens, None);
        assert_eq!(atom.bond_sum, None);
    }
}
