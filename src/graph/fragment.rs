// Molecular fragments implementation

use std::marker::PhantomData;
use crate::Element;
use crate::error::Error;
use crate::graph::{GraphAtom, GraphBond, GraphMolecule, AtomIndex};
use crate::graph::bond::BondOrder;
use crate::graph::builder::MoleculeBuilder;

/// Trait for molecular fragments that can be added to a molecule
pub trait Fragment {
    /// Add this fragment to a molecule builder starting from the active atom
    fn add_to_builder(&self, builder: MoleculeBuilder) -> Result<MoleculeBuilder, Error>;
    
    /// Create a standalone molecule from this fragment
    fn to_molecule(&self) -> Result<GraphMolecule, Error> {
        let builder = MoleculeBuilder::new();
        self.add_to_builder(builder)?.build()
    }
}

/// Ring fragment with configurable size and element type
pub struct Ring<E> {
    size: usize,
    element: E,
    bond_order: BondOrder,
    aromatic: bool,
    _phantom: PhantomData<E>,
}

impl<E: Into<GraphAtom> + Clone> Ring<E> {
    pub fn new(size: usize, element: E, bond_order: BondOrder) -> Self {
        if size < 3 {
            panic!("Ring size must be at least 3");
        }
        
        Self { 
            size, 
            element, 
            bond_order, 
            aromatic: false,
            _phantom: PhantomData,
        }
    }
    
    pub fn aromatic(mut self) -> Self {
        self.aromatic = true;
        self
    }
}

impl<E: Into<GraphAtom> + Clone> Fragment for Ring<E> {
    fn add_to_builder(&self, builder: MoleculeBuilder) -> Result<MoleculeBuilder, Error> {
        let mut builder = builder;
        let start_atom = builder.active_atom();
        
        // If no active atom, create the first atom of the ring
        let first_idx = if start_atom.is_none() {
            let atom_context = builder.atom(self.element.clone());
            let idx = atom_context.atom_idx();
            builder = atom_context.done();
            idx
        } else {
            start_atom.unwrap()
        };
        
        let mut prev_idx = first_idx;
        let mut indices = vec![prev_idx];
        
        // Add the rest of the ring atoms
        for _ in 1..self.size {
            let atom_context = builder.atom(self.element.clone());
            let new_idx = atom_context.atom_idx();
            builder = atom_context.done();
            
            // Connect to previous atom
            let bond_context = builder.bond(self.bond_order)
                .from(prev_idx)?
                .to(new_idx)?;
            builder = bond_context.done();
            
            indices.push(new_idx);
            prev_idx = new_idx;
        }
        
        // Close the ring
        let bond_context = builder.bond(self.bond_order)
            .from(prev_idx)?
            .to(first_idx)?;
        builder = bond_context.done();
        
        // If aromatic, add aromaticity property
        if self.aromatic {
            // Future: set aromaticity property when implemented
        }
        
        Ok(builder)
    }
}

/// Chain fragment with configurable length and element type
pub struct Chain<E> {
    length: usize,
    element: E,
    bond_order: BondOrder,
    _phantom: PhantomData<E>,
}

impl<E: Into<GraphAtom> + Clone> Chain<E> {
    pub fn new(length: usize, element: E, bond_order: BondOrder) -> Self {
        if length < 1 {
            panic!("Chain length must be at least 1");
        }
        
        Self { 
            length, 
            element, 
            bond_order,
            _phantom: PhantomData,
        }
    }
}

impl<E: Into<GraphAtom> + Clone> Fragment for Chain<E> {
    fn add_to_builder(&self, builder: MoleculeBuilder) -> Result<MoleculeBuilder, Error> {
        let mut builder = builder;
        let start_atom = builder.active_atom();
        
        // If no active atom, create the first atom of the chain
        if start_atom.is_none() {
            let atom_context = builder.atom(self.element.clone());
            builder = atom_context.done();
        }
        
        // Add the rest of the chain
        for _ in 1..self.length {
            let atom_context = builder.atom(self.element.clone());
            let new_idx = atom_context.atom_idx();
            builder = atom_context.done();
            
            // Bond to previous atom (which is the active atom)
            let prev_idx = builder.active_atom().unwrap();
            let bond_context = builder.bond(self.bond_order)
                .from(prev_idx)?
                .to(new_idx)?;
            builder = bond_context.done();
            
            // Make the new atom active
            builder = builder.set_active_atom(new_idx);
        }
        
        Ok(builder)
    }
}

/// Functional group fragments
pub enum FunctionalGroup {
    Hydroxyl,
    Carbonyl,
    Amino,
    Carboxyl,
    Methyl,
    Ethyl,
    Phenyl,
}

impl Fragment for FunctionalGroup {
    fn add_to_builder(&self, builder: MoleculeBuilder) -> Result<MoleculeBuilder, Error> {
        let active = builder.active_atom().ok_or_else(|| 
            Error::InvalidOperation("Cannot add functional group without an active atom".to_string())
        )?;
        
        match self {
            FunctionalGroup::Hydroxyl => {
                // Add OH group to active atom
                let atom_context = builder.atom(Element::O);
                let o_idx = atom_context.atom_idx();
                let builder = atom_context.done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(active)?
                    .to(o_idx)?;
                let builder = bond_context.done();
                
                // Add implicit hydrogen to O
                let selection_context = builder.select(o_idx)?;
                let builder = selection_context.modify(|atom| {
                    Ok(*atom = atom.with_implicit_hydrogens(1))
                })?.done();
                
                Ok(builder)
            },
            FunctionalGroup::Carbonyl => {
                // Add C=O group
                let atom_context = builder.atom(Element::O);
                let o_idx = atom_context.atom_idx();
                let builder = atom_context.done();
                
                let bond_context = builder.bond(BondOrder::Double)
                    .from(active)?
                    .to(o_idx)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
            FunctionalGroup::Amino => {
                // Add NH2 group
                let atom_context = builder.atom(Element::N);
                let n_idx = atom_context.atom_idx();
                let builder = atom_context.with_implicit_hydrogens(2).done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(active)?
                    .to(n_idx)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
            FunctionalGroup::Carboxyl => {
                // Add COOH group
                let atom_context = builder.atom(Element::C);
                let c_idx = atom_context.atom_idx();
                let builder = atom_context.done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(active)?
                    .to(c_idx)?;
                let builder = bond_context.done();
                
                // Add C=O
                let atom_context = builder.atom(Element::O);
                let o1_idx = atom_context.atom_idx();
                let builder = atom_context.done();
                
                let bond_context = builder.bond(BondOrder::Double)
                    .from(c_idx)?
                    .to(o1_idx)?;
                let builder = bond_context.done();
                
                // Add C-OH
                let atom_context = builder.atom(Element::O);
                let o2_idx = atom_context.atom_idx();
                let builder = atom_context.with_implicit_hydrogens(1).done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(c_idx)?
                    .to(o2_idx)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
            FunctionalGroup::Methyl => {
                // Add CH3 group
                let atom_context = builder.atom(Element::C);
                let c_idx = atom_context.atom_idx();
                let builder = atom_context.with_implicit_hydrogens(3).done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(active)?
                    .to(c_idx)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
            FunctionalGroup::Ethyl => {
                // Add CH2CH3 group
                let builder = FunctionalGroup::Methyl.add_to_builder(builder)?;
                
                // Add another carbon - the methyl carbon should be active
                let c1_idx = builder.active_atom().unwrap();
                
                let atom_context = builder.atom(Element::C);
                let c2_idx = atom_context.atom_idx();
                let builder = atom_context.with_implicit_hydrogens(2).done();
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(c1_idx)?
                    .to(c2_idx)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
            FunctionalGroup::Phenyl => {
                // Reuse benzene ring implementation
                let benzene = Ring::new(6, Element::C, BondOrder::Single).aromatic();
                benzene.add_to_builder(builder)?;
                
                // Connect to active atom
                let active = builder.active_atom().unwrap();
                let ring_start = builder.active_atom().unwrap(); // Use last added atom from ring
                
                let bond_context = builder.bond(BondOrder::Single)
                    .from(active)?
                    .to(ring_start)?;
                let builder = bond_context.done();
                
                Ok(builder)
            },
        }
    }
}

/// Common fragment implementations as functions
pub fn benzene() -> Ring<Element> {
    Ring::new(6, Element::C, BondOrder::Single).aromatic()
}

pub fn cyclohexane() -> Ring<Element> {
    Ring::new(6, Element::C, BondOrder::Single)
}

pub fn cyclopentane() -> Ring<Element> {
    Ring::new(5, Element::C, BondOrder::Single)
}

pub fn methane() -> impl Fragment {
    FunctionalGroup::Methyl
}