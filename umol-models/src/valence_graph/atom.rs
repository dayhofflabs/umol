//! Valence atom type and builder
//!
//! Valence atom is the node type of valence graphs and is defined by its element and
//! properties. It's strictly typed, meaning that it has to match one of the predefined
//! atom specs. It cannot be created directly, but only through the `AtomBuilder` type,
//! which ensures correct typing.

use crate::AtomSpec;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use umol::{error::DataError, Result};
use umol_data::Element;

use super::{AtomMatcher, AtomValidator, DEFAULT_ATOM_MATCHER, DEFAULT_ATOM_VALIDATOR};

/// Valence atom type including strict typing. Cannot be created directly, but only through
/// the `AtomBuilder` type, which performs validation of the atom properties. Mutations are
/// possible by converting back to a builder using the `to_builder` method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Atom {
    element: Element,
    charge: i8,
    lone_pairs: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: u8,
    implicit_hydrogens: u8,
    valence: u8,
}

impl Atom {
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

    pub fn valence(&self) -> u8 {
        self.valence
    }

    pub fn to_builder(self) -> AtomBuilder {
        AtomBuilder {
            element: self.element,
            charge: Some(self.charge),
            lone_pairs: Some(self.lone_pairs),
            donated_pairs: Some(self.donated_pairs),
            accepted_pairs: Some(self.accepted_pairs),
            unpaired_electrons: Some(self.unpaired_electrons),
            multiplicity: Some(self.multiplicity),
            implicit_hydrogens: Some(self.implicit_hydrogens),
            valence: Some(self.valence),
        }
    }

    pub fn to_spec(self) -> AtomSpec {
        AtomSpec::new(
            self.element,
            self.charge,
            self.lone_pairs,
            self.donated_pairs,
            self.accepted_pairs,
            self.unpaired_electrons,
            self.multiplicity,
            self.implicit_hydrogens,
            self.valence,
        )
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_spec())
    }
}

impl From<Atom> for AtomBuilder {
    fn from(atom: Atom) -> Self {
        atom.to_builder()
    }
}

/// Builder type for creating and mutating `Atom` types including strict typing.
/// The resulting `Atom` objects must match the predefined `AtomType` types.
#[derive(Debug)]
pub struct AtomBuilder {
    element: Element,
    charge: Option<i8>,
    lone_pairs: Option<u8>,
    donated_pairs: Option<u8>,
    accepted_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<u8>,
    implicit_hydrogens: Option<u8>,
    valence: Option<u8>,
}

impl AtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: None,
            lone_pairs: None,
            donated_pairs: None,
            accepted_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            implicit_hydrogens: None,
            valence: None,
        }
    }

    pub fn from_spec(atom_spec: AtomSpec) -> Self {
        Self {
            element: atom_spec.element(),
            charge: Some(atom_spec.charge()),
            lone_pairs: Some(atom_spec.lone_pairs()),
            donated_pairs: Some(atom_spec.donated_pairs()),
            accepted_pairs: Some(atom_spec.accepted_pairs()),
            unpaired_electrons: Some(atom_spec.unpaired_electrons()),
            multiplicity: Some(atom_spec.multiplicity()),
            implicit_hydrogens: Some(atom_spec.implicit_hydrogens()),
            valence: Some(atom_spec.valence()),
        }
    }

    pub fn element(&self) -> Element {
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

    pub fn valence(&self) -> Option<u8> {
        self.valence
    }

    pub fn set_element(&mut self, element: Element) -> &mut Self {
        self.element = element;
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

    pub fn set_valence(&mut self, valence: u8) -> &mut Self {
        self.valence = Some(valence);
        self
    }

    pub fn update_element(&mut self, f: impl FnOnce(Element) -> Element) -> &mut Self {
        self.element = f(self.element);
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
        self.multiplicity = Some(f(self.multiplicity.unwrap_or(1)));
        self
    }

    pub fn update_implicit_hydrogens(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.implicit_hydrogens = Some(f(self.implicit_hydrogens.unwrap_or(0)));
        self
    }

    pub fn update_valence(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.valence = Some(f(self.valence.unwrap_or(0)));
        self
    }

    pub fn build(self) -> Result<Atom> {
        self.build_with(&DEFAULT_ATOM_VALIDATOR, &DEFAULT_ATOM_MATCHER)
    }

    pub fn build_with(self, validator: &AtomValidator, matcher: &AtomMatcher) -> Result<Atom> {
        validator.validate(&self)?;
        let atom_specs = matcher.find(&self)?;
        if atom_specs.is_empty() {
            return Err(DataError::NoAtomSpec(format!("{:?}", self)).into());
        } else if atom_specs.len() > 1 {
            return Err(DataError::MultipleAtomSpecs(format!(
                "{:?}: {}",
                self,
                atom_specs
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            ))
            .into());
        }
        let atom_spec = atom_specs.first().unwrap();
        Ok(Atom {
            element: atom_spec.element(),
            charge: atom_spec.charge(),
            lone_pairs: atom_spec.lone_pairs(),
            donated_pairs: atom_spec.donated_pairs(),
            accepted_pairs: atom_spec.accepted_pairs(),
            unpaired_electrons: atom_spec.unpaired_electrons(),
            multiplicity: atom_spec.multiplicity(),
            implicit_hydrogens: atom_spec.implicit_hydrogens(),
            valence: atom_spec.valence(),
        })
    }
}

impl From<AtomSpec> for AtomBuilder {
    fn from(atom_spec: AtomSpec) -> Self {
        AtomBuilder::from_spec(atom_spec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomSpec;
    use umol_data::e;

    #[test]
    fn test_atom_display() {
        let mut builder = AtomBuilder::new(e!(C));
        builder
            .set_charge(0)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_valence(0);

        let atom = builder.build().unwrap();
        assert_eq!(format!("{}", atom), "[C/^2]");
    }

    #[test]
    fn test_atom_serialize() {
        let mut builder = AtomBuilder::new(e!(C));
        builder
            .set_charge(0)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_valence(0);

        let atom = builder.build().unwrap();
        let serialized = serde_json::to_string(&atom).unwrap();
        assert_eq!(
            serialized,
            "{\"element\":\"C\",\"charge\":0,\"lone_pairs\":1,\"donated_pairs\":0,\"accepted_pairs\":0,\"unpaired_electrons\":2,\"multiplicity\":3,\"implicit_hydrogens\":0,\"valence\":0}"
        );

        let deserialized: Atom = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, atom);
    }

    #[test]
    fn test_atom_to_builder() {
        let atom = AtomBuilder::new(e!(Ne)).build().unwrap();
        let builder = atom.to_builder();
        assert_eq!(builder.element(), e!(Ne));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(4));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(0));
        assert_eq!(builder.multiplicity(), Some(1));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.valence(), Some(0));
    }

    #[test]
    fn test_atom_builder_new() {
        let atom = AtomBuilder::new(e!(Ne)).build().unwrap();
        assert_eq!(atom.element(), e!(Ne));
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.lone_pairs(), 4);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
        assert_eq!(atom.multiplicity(), 1);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.valence(), 0);
    }

    #[test]
    fn test_atom_builder_from_atom_type() {
        let state = AtomSpec::new(e!(C), 0, 1, 0, 0, 2, 3, 0, 4);
        let builder = AtomBuilder::from_spec(state);
        assert_eq!(builder.element(), e!(C));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(1));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(2));
        assert_eq!(builder.multiplicity(), Some(3));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.valence(), Some(4));
    }

    #[test]
    fn test_atom_builder_properties() {
        let mut builder = AtomBuilder::new(e!(N));
        builder
            .set_charge(-1)
            .set_lone_pairs(1)
            .set_donated_pairs(0)
            .set_accepted_pairs(0)
            .set_unpaired_electrons(2)
            .set_multiplicity(3)
            .set_implicit_hydrogens(0)
            .set_valence(1);

        assert_eq!(builder.element(), e!(N));
        assert_eq!(builder.charge(), Some(-1));
        assert_eq!(builder.lone_pairs(), Some(1));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(2));
        assert_eq!(builder.multiplicity(), Some(3));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.valence(), Some(1));
    }

    #[test]
    fn test_atom_builder_set() {
        let mut builder = AtomBuilder::new(e!(C));
        builder.set_charge(0);
        builder.set_lone_pairs(1);
        builder.set_donated_pairs(0);
        builder.set_accepted_pairs(0);
        builder.set_unpaired_electrons(2);
        builder.set_multiplicity(3);
        builder.set_implicit_hydrogens(0);
        builder.set_valence(2);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), e!(C));
        assert_eq!(atom.charge(), 0);
        assert_eq!(atom.lone_pairs(), 1);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 2);
        assert_eq!(atom.multiplicity(), 3);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.valence(), 2);
    }

    #[test]
    fn test_atom_builder_update() {
        let mut builder = AtomBuilder::new(e!(C));
        builder.update_element(|elem| elem.next().unwrap());
        builder.update_charge(|x| x + 1);
        builder.update_valence(|x| x + 4);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), e!(N));
        assert_eq!(atom.charge(), 1);
        assert_eq!(atom.lone_pairs(), 0);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
        assert_eq!(atom.multiplicity(), 1);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.valence(), 4);
    }

    #[test]
    fn test_atom_builder_build() {
        let mut builder = AtomBuilder::new(e!(N));
        builder
            .set_charge(1)
            .set_valence(4);

        let atom = builder.build().unwrap();
        assert_eq!(atom.element(), e!(N));
        assert_eq!(atom.charge(), 1);
        assert_eq!(atom.lone_pairs(), 0);
        assert_eq!(atom.donated_pairs(), 0);
        assert_eq!(atom.accepted_pairs(), 0);
        assert_eq!(atom.unpaired_electrons(), 0);
        assert_eq!(atom.multiplicity(), 1);
        assert_eq!(atom.implicit_hydrogens(), 0);
        assert_eq!(atom.valence(), 4);
    }

    #[test]
    fn test_atom_builder_validation() {
        let mut builder = AtomBuilder::new(e!(C));
        builder.set_unpaired_electrons(3);
        let result = builder.build();
        assert!(result.is_ok());
    }

    #[test]
    fn test_atom_type_into_atom_builder() {
        let atom: AtomBuilder = AtomSpec::new(e!(C), 0, 1, 0, 0, 2, 3, 0, 4).into();
        assert_eq!(atom.element(), e!(C));
        assert_eq!(atom.charge(), Some(0));
        assert_eq!(atom.lone_pairs(), Some(1));
        assert_eq!(atom.donated_pairs(), Some(0));
        assert_eq!(atom.accepted_pairs(), Some(0));
        assert_eq!(atom.unpaired_electrons(), Some(2));
        assert_eq!(atom.multiplicity(), Some(3));
        assert_eq!(atom.implicit_hydrogens(), Some(0));
        assert_eq!(atom.valence(), Some(4));
    }

    #[test]
    fn test_atom_into_atom_builder() {
        let mut atom = AtomBuilder::new(e!(C));
        atom.set_unpaired_electrons(2);
        atom.set_multiplicity(3);
        atom.set_valence(2);
        let atom = atom.build().unwrap();
        let builder: AtomBuilder = atom.into();
        assert_eq!(builder.element(), e!(C));
        assert_eq!(builder.charge(), Some(0));
        assert_eq!(builder.lone_pairs(), Some(1));
        assert_eq!(builder.donated_pairs(), Some(0));
        assert_eq!(builder.accepted_pairs(), Some(0));
        assert_eq!(builder.unpaired_electrons(), Some(2));
        assert_eq!(builder.multiplicity(), Some(3));
        assert_eq!(builder.implicit_hydrogens(), Some(0));
        assert_eq!(builder.valence(), Some(2));
    }
}
