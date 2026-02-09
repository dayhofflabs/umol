//! Reaction IR for TableIR.
//!
//! `Reaction` is TableIR reaction representation for reactions with fixed composition
//! Delaying the implementation for `ExtendedReaction` as long its semantics is unclear.
//! (ordinary reactions between extended molecules (variable composition, polymers, etc.)
//! and reaction templates are distinct use cases).

use std::collections::BTreeMap;

use indexmap::IndexMap;

use super::molecule::{ExtendedMolecule, Molecule};
use super::source::SourceFormat;

/// Basic reaction IR
#[derive(Clone, Debug, PartialEq)]
pub struct Reaction {
    pub reactants: Vec<Molecule>,
    pub products: Vec<Molecule>,
    pub agents: Vec<Molecule>,
    pub atom_mapping: BTreeMap<u32, (Vec<(usize, usize)>, Vec<(usize, usize)>)>,

    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl Reaction {
    pub fn empty() -> Self {
        Self {
            reactants: Vec::new(),
            products: Vec::new(),
            agents: Vec::new(),
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn from_molecules(
        reactants: Vec<Molecule>,
        products: Vec<Molecule>,
        agents: Vec<Molecule>,
    ) -> Self {
        Self {
            reactants,
            products,
            agents,
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn reactant_count(&self) -> usize {
        self.reactants.len()
    }

    pub fn product_count(&self) -> usize {
        self.products.len()
    }

    pub fn agent_count(&self) -> usize {
        self.agents.len()
    }

    pub fn add_reactant(&mut self, molecule: Molecule) {
        self.reactants.push(molecule);
    }

    pub fn add_product(&mut self, molecule: Molecule) {
        self.products.push(molecule);
    }

    pub fn add_agent(&mut self, molecule: Molecule) {
        self.agents.push(molecule);
    }
}

/// Extended reaction IR
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedReaction {
    pub reactants: Vec<ExtendedMolecule>,
    pub products: Vec<ExtendedMolecule>,
    pub agents: Vec<ExtendedMolecule>,
    pub atom_mapping: BTreeMap<u32, (Vec<(usize, usize)>, Vec<(usize, usize)>)>,

    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl ExtendedReaction {
    pub fn empty() -> Self {
        Self {
            reactants: Vec::new(),
            products: Vec::new(),
            agents: Vec::new(),
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn from_extended_molecules(
        reactants: Vec<ExtendedMolecule>,
        products: Vec<ExtendedMolecule>,
        agents: Vec<ExtendedMolecule>,
    ) -> Self {
        Self {
            reactants,
            products,
            agents,
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn reactant_count(&self) -> usize {
        self.reactants.len()
    }

    pub fn product_count(&self) -> usize {
        self.products.len()
    }

    pub fn agent_count(&self) -> usize {
        self.agents.len()
    }

    pub fn add_reactant(&mut self, molecule: ExtendedMolecule) {
        self.reactants.push(molecule);
    }

    pub fn add_product(&mut self, molecule: ExtendedMolecule) {
        self.products.push(molecule);
    }

    pub fn add_agent(&mut self, molecule: ExtendedMolecule) {
        self.agents.push(molecule);
    }
}
