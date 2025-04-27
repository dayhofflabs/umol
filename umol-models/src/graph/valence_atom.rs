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
    donated_pairs: u8,
    accepted_pairs: u8,
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

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
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

    pub fn to_builder(self) -> ValenceAtomBuilder {
        ValenceAtomBuilder::from(self)
    }
}

impl Display for ValenceAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        ValenceState::new(
            self.element,
            self.charge,
            self.lone_pairs,
            self.donated_pairs,
            self.accepted_pairs,
            self.unpaired_electrons,
            self.multiplicity,
            self.bond_sum + self.implicit_hydrogens,
        )
        .fmt(f)
    }
}

#[derive(Debug)]
pub struct ValenceAtomBuilder {
    element: Option<Element>,
    charge: Option<i8>,
    lone_pairs: Option<u8>,
    donated_pairs: Option<u8>,
    accepted_pairs: Option<u8>,
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
            donated_pairs: None,
            accepted_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            implicit_hydrogens: None,
            bond_sum: None,
        }
    }

    pub fn from_valence_state(state: ValenceState) -> Self {
        Self {
            element: Some(state.element()),
            charge: Some(state.charge()),
            lone_pairs: Some(state.lone_pairs()),
            donated_pairs: Some(state.donated_pairs()),
            accepted_pairs: Some(state.accepted_pairs()),
            unpaired_electrons: Some(state.unpaired_electrons()),
            multiplicity: Some(state.multiplicity()),
            implicit_hydrogens: None,
            bond_sum: None,
        }
    }

    pub fn element(&self) -> Option<Element> {
        self.element
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn lone_pairs(&self) -> Option<u8> {
        self.lone_pairs
    }

    pub fn donated_pairs(&self) -> Option<u8> {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> Option<u8> {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> Option<u8> {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> Option<u8> {
        self.multiplicity
    }

    pub fn implicit_hydrogens(&self) -> Option<u8> {
        self.implicit_hydrogens
    }

    pub fn bond_sum(&self) -> Option<u8> {
        self.bond_sum
    }

    pub fn set_element(&mut self, element: Element) -> &mut Self {
        self.element = Some(element);
        self
    }

    pub fn set_charge(&mut self, charge: i8) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn set_lone_pairs(&mut self, lone_pairs: u8) -> &mut Self {
        self.lone_pairs = Some(lone_pairs);
        self
    }

    pub fn set_donated_pairs(&mut self, donated_pairs: u8) -> &mut Self {
        self.donated_pairs = Some(donated_pairs);
        self
    }

    pub fn set_accepted_pairs(&mut self, accepted_pairs: u8) -> &mut Self {
        self.accepted_pairs = Some(accepted_pairs);
        self
    }

    pub fn set_unpaired_electrons(&mut self, unpaired_electrons: u8) -> &mut Self {
        self.unpaired_electrons = Some(unpaired_electrons);
        self
    }

    pub fn set_multiplicity(&mut self, multiplicity: u8) -> &mut Self {
        self.multiplicity = Some(multiplicity);
        self
    }

    pub fn set_implicit_hydrogens(&mut self, implicit_hydrogens: u8) -> &mut Self {
        self.implicit_hydrogens = Some(implicit_hydrogens);
        self
    }

    pub fn set_bond_sum(&mut self, bond_sum: u8) -> &mut Self {
        self.bond_sum = Some(bond_sum);
        self
    }

    pub fn update_element(&mut self, f: impl FnOnce(Element) -> Element) -> &mut Self {
        self.element = Some(f(self.element.unwrap()));
        self
    }

    pub fn update_charge(&mut self, f: impl FnOnce(i8) -> i8) -> &mut Self {
        self.charge = Some(f(self.charge.unwrap_or(0)));
        self
    }

    pub fn update_lone_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.lone_pairs = Some(f(self.lone_pairs.unwrap_or(0)));
        self
    }

    pub fn update_donated_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.donated_pairs = Some(f(self.donated_pairs.unwrap_or(0)));
        self
    }

    pub fn update_accepted_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.accepted_pairs = Some(f(self.accepted_pairs.unwrap_or(0)));
        self
    }

    pub fn update_unpaired_electrons(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.unpaired_electrons = Some(f(self.unpaired_electrons.unwrap_or(0)));
        self
    }

    pub fn update_multiplicity(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.multiplicity = Some(f(self.multiplicity.unwrap_or(0)));
        self
    }

    pub fn update_implicit_hydrogens(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.implicit_hydrogens = Some(f(self.implicit_hydrogens.unwrap_or(0)));
        self
    }

    pub fn update_bond_sum(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.bond_sum = Some(f(self.bond_sum.unwrap_or(0)));
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
            donated_pairs: self.donated_pairs.unwrap_or(0),
            accepted_pairs: self.accepted_pairs.unwrap_or(0),
            unpaired_electrons: self.unpaired_electrons.unwrap_or(0),
            multiplicity: self.multiplicity.unwrap_or(0),
            implicit_hydrogens: self.implicit_hydrogens.unwrap_or(0),
            bond_sum: self.bond_sum.unwrap_or(0),
        })
    }
}

impl From<ValenceState> for ValenceAtomBuilder {
    fn from(state: ValenceState) -> Self {
        ValenceAtomBuilder::from_valence_state(state)
    }
}

impl From<ValenceAtom> for ValenceAtomBuilder {
    fn from(atom: ValenceAtom) -> Self {
        ValenceAtomBuilder {
            element: Some(atom.element()),
            charge: Some(atom.charge()),
            lone_pairs: Some(atom.lone_pairs()),
            donated_pairs: Some(atom.donated_pairs()),
            accepted_pairs: Some(atom.accepted_pairs()),
            unpaired_electrons: Some(atom.unpaired_electrons()),
            multiplicity: Some(atom.multiplicity()),
            implicit_hydrogens: Some(atom.implicit_hydrogens()),
            bond_sum: Some(atom.bond_sum()),
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
            .set_charge(0)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_bond_sum(0);

        let atom = builder.build().unwrap();
        assert_eq!(format!("{}", atom), "[C/^2]");
    }

    #[test]
    fn test_valence_atom_serialize() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder
            .set_charge(0)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_bond_sum(0);

        let atom = builder.build().unwrap();
        let serialized = serde_json::to_string(&atom).unwrap();
        assert_eq!(
            serialized,
            "{\"element\":\"C\",\"charge\":0,\"lone_pairs\":1,\"donated_pairs\":0,\"accepted_pairs\":0,\"unpaired_electrons\":2,\"multiplicity\":3,\"implicit_hydrogens\":0,\"bond_sum\":0}"
        );

        let deserialized: ValenceAtom = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, atom);
    }

    #[test]
    fn test_valence_atom_to_builder() {
        let atom = ValenceAtomBuilder::new(Element::C).build().unwrap();
        let builder = atom.to_builder();
        assert_eq!(builder.element(), Some(Element::C));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(0));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(0));
        assert_eq!(builder.multiplicity(), Some(0));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.bond_sum(), Some(0));
    }

    #[test]
    fn test_valence_atom_builder_new() {
        let atom = ValenceAtomBuilder::new(Element::C).build().unwrap();
        assert_eq!(atom.element(), Element::C);
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.lone_pairs(), 0);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
        assert_eq!(atom.multiplicity(), 0);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.bond_sum(), 0);
    }

    #[test]
    fn test_valence_atom_builder_from_valence_state() {
        let state = ValenceState::new(Element::C, 0, 1, 0, 0, 2, 3, 4);
        let builder = ValenceAtomBuilder::from_valence_state(state);
        assert_eq!(builder.element(), Some(Element::C));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(1));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(2));
        assert_eq!(builder.multiplicity(), Some(3));
        assert_eq!(builder.implicit_hydrogens(), None);
        assert_eq!(builder.bond_sum(), None);
    }

    #[test]
    fn test_valence_atom_builder_properties() {
        let mut builder = ValenceAtomBuilder::new(Element::N);
        builder
            .set_charge(-1)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_bond_sum(1);

        assert_eq!(builder.element(), Some(Element::N));
        assert_eq!(builder.charge(), Some(-1));
        assert_eq!(builder.lone_pairs(), Some(1));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(2));
        assert_eq!(builder.multiplicity(), Some(3));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.bond_sum(), Some(1));
    }

    #[test]
    fn test_valence_atom_builder_set() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder.set_charge(0);
        builder.set_lone_pairs(1);
        builder.set_donated_pairs(0);
        builder.set_accepted_pairs(0);
        builder.set_unpaired_electrons(2);
        builder.set_multiplicity(3);
        builder.set_implicit_hydrogens(0);
        builder.set_bond_sum(2);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), Element::C);
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.multiplicity(), 3);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.bond_sum(), 2);
    }

    #[test]
    fn test_valence_atom_builder_update() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder.update_element(|x| x.next().unwrap());
        builder.update_charge(|x| x + 1);
        builder.update_lone_pairs(|x| x + 1);
        builder.update_donated_pairs(|x| x + 1);
        builder.update_implicit_hydrogens(|x| x + 1);
        builder.update_bond_sum(|x| x + 2);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), Element::N);
        assert_eq!(atom.charge(), 1);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.donated_pairs(), 1);
        assert_eq!(atom.implicit_hydrogens(), 1);
        assert_eq!(atom.bond_sum(), 2);
    }

    #[test]
    fn test_valence_atom_builder_build() {
        let mut builder = ValenceAtomBuilder::new(Element::N);
        builder
            .set_charge(-1)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_bond_sum(1);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), Element::N);
        assert_eq!(atom.charge(), -1);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.multiplicity(), 3);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.bond_sum(), 1);
    }

    #[test]
    fn test_valence_atom_builder_validation() {
        let mut builder = ValenceAtomBuilder::new(Element::C);
        builder.set_charge(-1);
        builder.set_unpaired_electrons(3);
        builder.set_implicit_hydrogens(0);
        let result = builder.build();
        assert!(result.is_ok());
    }

    #[test]
    fn test_valence_state_into_valence_atom_builder() {
        let atom: ValenceAtomBuilder = ValenceState::new(Element::C, 0, 1, 0, 0, 2, 3, 4).into();
        assert_eq!(atom.element(), Some(Element::C));
        assert_eq!(atom.charge(), Some(0));
        assert_eq!(atom.lone_pairs(), Some(1));
        assert_eq!(atom.donated_pairs(), Some(0));
        assert_eq!(atom.accepted_pairs(), Some(0));
        assert_eq!(atom.unpaired_electrons(), Some(2));
        assert_eq!(atom.multiplicity(), Some(3));
        assert_eq!(atom.implicit_hydrogens(), None);
        assert_eq!(atom.bond_sum(), None);
    }

    #[test]
    fn test_valence_atom_into_valence_atom_builder() {
        let atom = ValenceAtomBuilder::new(Element::C).build().unwrap();
        let builder: ValenceAtomBuilder = atom.into();
        assert_eq!(builder.element(), Some(Element::C));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(0));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(0));
        assert_eq!(builder.multiplicity(), Some(0));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.bond_sum(), Some(0));
    }
}
