//! Valence-aware Atom type copied from `umol-models-valence` for GraphIR usage.

use std::fmt::{self, Display};

use serde::{Deserialize, Serialize};
use umol_data::Element;

use super::atom_matcher::{AtomMatcher, DEFAULT_ATOM_MATCHER};
use super::atom_spec::AtomSpec;
use super::atom_validator::{AtomValidator, DEFAULT_ATOM_VALIDATOR};
use super::error::GraphError;
use crate::position::Point3D;
use crate::table_ir::Chirality;
use crate::span::Span;

type Result<T> = std::result::Result<T, GraphError>;

/// Valence atom type including strict typing.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    element: Element,
    charge: i32,
    valence: u32,
    multiplicity: u32,
    lone_pairs: u32,
    donated_pairs: u32,
    accepted_pairs: u32,
    unpaired_electrons: u32,
    implicit_hydrogens: u32,
    aromatic: Option<bool>,
    chirality: Option<Chirality>,
    class_num: Option<u32>,
    isotope_mass: Option<u32>,
    position: Option<Point3D>,
    span: Option<Span>,
}

impl Atom {
    pub fn element(&self) -> Element {
        self.element
    }

    pub fn charge(&self) -> i32 {
        self.charge
    }

    pub fn lone_pairs(&self) -> u32 {
        self.lone_pairs
    }

    pub fn donated_pairs(&self) -> u32 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u32 {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> u32 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    pub fn implicit_hydrogens(&self) -> u32 {
        self.implicit_hydrogens
    }

    pub fn valence(&self) -> u32 {
        self.valence
    }

    pub fn position(&self) -> Option<&Point3D> {
        self.position.as_ref()
    }

    pub fn isotope(&self) -> Option<u32> {
        self.isotope
    }

    pub fn aromatic(&self) -> Option<bool> {
        self.aromatic
    }

    pub fn chirality(&self) -> Option<Chirality> {
        self.chirality
    }

    pub fn class_num(&self) -> Option<u32> {
        self.class_num
    }

    pub fn span(&self) -> Option<Span> {
        self.span
    }

    pub fn span_bytes(&self) -> Option<(u32, u32)> {
        self.span.and_then(|s| s.bytes_range())
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
            valence: Some(self.valence),
            implicit_hydrogens: Some(self.implicit_hydrogens),
            position: self.position,
            isotope: self.isotope_mass,
            aromatic: self.aromatic,
            chirality: self.chirality,
            class_num: self.class_num,
            span: self.span,
        }
    }

    pub fn to_spec(&self) -> AtomSpec {
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
#[derive(Debug)]
pub struct AtomBuilder {
    element: Element,
    charge: Option<i32>,
    valence: Option<u32>,
    multiplicity: Option<u32>,
    lone_pairs: Option<u32>,
    donated_pairs: Option<u32>,
    accepted_pairs: Option<u32>,
    unpaired_electrons: Option<u32>,
    implicit_hydrogens: Option<u32>,
    position: Option<Point3D>,
    isotope_mass: Option<u32>,
    aromatic: Option<bool>,
    chirality: Option<Chirality>,
    class_num: Option<u32>,
    span: Option<Span>,
}

impl AtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: None,
            multiplicity: None,
            valence: None,
            lone_pairs: None,
            donated_pairs: None,
            accepted_pairs: None,
            unpaired_electrons: None,
            implicit_hydrogens: None,
            position: None,
            isotope: None,
            aromatic: None,
            chirality: None,
            class_num: None,
            span: None,
        }
    }

    pub fn from_spec(atom_spec: AtomSpec) -> Self {
        Self {
            element: atom_spec.element(),
            charge: Some(atom_spec.charge()),
            valence: Some(atom_spec.valence()),
            multiplicity: Some(atom_spec.multiplicity()),
            lone_pairs: Some(atom_spec.lone_pairs()),
            donated_pairs: Some(atom_spec.donated_pairs()),
            accepted_pairs: Some(atom_spec.accepted_pairs()),
            unpaired_electrons: Some(atom_spec.unpaired_electrons()),
            implicit_hydrogens: Some(atom_spec.implicit_hydrogens()),
            position: None,
            isotope: None,
            aromatic: None,
            chirality: None,
            class_num: None,
            span: None,
        }
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn charge(&self) -> Option<i32> {
        self.charge
    }

    pub fn lone_pairs(&self) -> Option<u32> {
        self.lone_pairs
    }

    pub fn donated_pairs(&self) -> Option<u32> {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> Option<u32> {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> Option<u32> {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> Option<u32> {
        self.multiplicity
    }

    pub fn valence(&self) -> Option<u32> {
        self.valence
    }

    pub fn implicit_hydrogens(&self) -> Option<u32> {
        self.implicit_hydrogens
    }

    pub fn position(&self) -> Option<&Point3D> {
        self.position.as_ref()
    }

    pub fn isotope(&self) -> Option<u32> {
        self.isotope
    }

    pub fn aromatic(&self) -> Option<bool> {
        self.aromatic
    }

    pub fn chirality(&self) -> Option<Chirality> {
        self.chirality
    }

    pub fn class_num(&self) -> Option<u32> {
        self.class_num
    }

    pub fn span(&self) -> Option<Span> {
        self.span
    }

    pub fn span_bytes(&self) -> Option<(u32, u32)> {
        self.span.and_then(|s| s.bytes_range())
    }

    pub fn set_element(&mut self, element: Element) -> &mut Self {
        self.element = element;
        self
    }

    pub fn set_charge(&mut self, charge: i32) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn set_lone_pairs(&mut self, count: u32) -> &mut Self {
        self.lone_pairs = Some(count);
        self
    }

    pub fn set_donated_pairs(&mut self, count: u32) -> &mut Self {
        self.donated_pairs = Some(count);
        self
    }

    pub fn set_accepted_pairs(&mut self, count: u32) -> &mut Self {
        self.accepted_pairs = Some(count);
        self
    }

    pub fn set_unpaired_electrons(&mut self, count: u32) -> &mut Self {
        self.unpaired_electrons = Some(count);
        self
    }

    pub fn set_multiplicity(&mut self, multiplicity: u32) -> &mut Self {
        self.multiplicity = Some(multiplicity);
        self
    }

    pub fn set_valence(&mut self, valence: u32) -> &mut Self {
        self.valence = Some(valence);
        self
    }

    pub fn set_implicit_hydrogens(&mut self, count: u32) -> &mut Self {
        self.implicit_hydrogens = Some(count);
        self
    }

    pub fn set_position(&mut self, position: Point3D) -> &mut Self {
        self.position = Some(position);
        self
    }

    pub fn set_isotope(&mut self, isotope: u32) -> &mut Self {
        self.isotope_mass = Some(isotope);
        self
    }

    pub fn set_aromatic(&mut self, value: bool) -> &mut Self {
        self.aromatic = Some(value);
        self
    }

    pub fn set_chirality(&mut self, chirality: Chirality) -> &mut Self {
        self.chirality = Some(chirality);
        self
    }

    pub fn set_class_num(&mut self, class_num: u32) -> &mut Self {
        self.class_num = Some(class_num);
        self
    }

    pub fn set_span(&mut self, start: Option<u32>, end: Option<u32>) -> &mut Self {
        self.span = Span::from_bytes_opt(start, end);
        self
    }

    pub fn set_span_opt(&mut self, span: Option<Span>) -> &mut Self {
        self.span = span;
        self
    }

    pub fn update_element(&mut self, f: impl FnOnce(Element) -> Element) -> &mut Self {
        self.element = f(self.element);
        self
    }

    pub fn update_charge(&mut self, f: impl FnOnce(i32) -> i32) -> &mut Self {
        self.charge = Some(f(self.charge.unwrap_or(0)));
        self
    }

    pub fn update_lone_pairs(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.lone_pairs = Some(f(self.lone_pairs.unwrap_or(0)));
        self
    }

    pub fn update_donated_pairs(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.donated_pairs = Some(f(self.donated_pairs.unwrap_or(0)));
        self
    }

    pub fn update_accepted_pairs(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.accepted_pairs = Some(f(self.accepted_pairs.unwrap_or(0)));
        self
    }

    pub fn update_unpaired_electrons(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.unpaired_electrons = Some(f(self.unpaired_electrons.unwrap_or(0)));
        self
    }

    pub fn update_multiplicity(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.multiplicity = Some(f(self.multiplicity.unwrap_or(1)));
        self
    }

    pub fn update_valence(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.valence = Some(f(self.valence.unwrap_or(0)));
        self
    }

    pub fn build(self) -> Result<Atom> {
        self.build_with(&DEFAULT_ATOM_VALIDATOR, &DEFAULT_ATOM_MATCHER)
    }

    pub fn build_with(self, validator: &AtomValidator, matcher: &AtomMatcher) -> Result<Atom> {
        let atom_specs = matcher.find(&self)?;
        if atom_specs.is_empty() {
            return Err(GraphError::InvalidAtomSpec(format!("{:?}", self)));
        } else if atom_specs.len() > 1 {
            return Err(GraphError::InvalidAtomSpec(format!(
                "{:?}: {}",
                self,
                atom_specs
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            )));
        }
        let atom_spec = atom_specs.first().unwrap();
        let atom = Atom {
            element: atom_spec.element(),
            charge: self.charge.unwrap_or(atom_spec.charge()),
            valence: self.valence.unwrap_or_else(|| atom_spec.valence()),
            multiplicity: self
                .multiplicity
                .unwrap_or_else(|| atom_spec.multiplicity()),
            lone_pairs: self.lone_pairs.unwrap_or_else(|| atom_spec.lone_pairs()),
            donated_pairs: self
                .donated_pairs
                .unwrap_or_else(|| atom_spec.donated_pairs()),
            accepted_pairs: self
                .accepted_pairs
                .unwrap_or_else(|| atom_spec.accepted_pairs()),
            unpaired_electrons: self.unpaired_electrons.unwrap_or_else(|| atom_spec.unpaired_electrons()),
            implicit_hydrogens: self.implicit_hydrogens.unwrap_or_else(|| atom_spec.implicit_hydrogens()),
            position: self.position,
            isotope: self.isotope_mass,
            aromatic: self.aromatic,
            chirality: self.chirality,
            class_num: self.class_num,
            span: self.span,
        };
        validator.validate(&atom)?;
        Ok(atom)
    }
}

impl From<AtomSpec> for AtomBuilder {
    fn from(atom_spec: AtomSpec) -> Self {
        AtomBuilder::from_spec(atom_spec)
    }
}

impl From<Element> for AtomBuilder {
    fn from(element: Element) -> Self {
        AtomBuilder::new(element)
    }
}
