// Fluent builder API for graph molecules

use crate::{Error, Element};
use crate::graph::{GraphAtom, GraphBond, GraphMolecule, AtomIndex, BondIndex};
use crate::graph::bond::BondOrder;
use crate::graph::fragment::Fragment;
use crate::validation::ValidationSet;
use petgraph::stable_graph::{StableGraph, NodeIndex, EdgeIndex};
use std::collections::{HashMap, HashSet};

/// The main molecule builder
pub struct MoleculeBuilder {
    graph: StableGraph<GraphAtom, GraphBond>,
    active_atom: Option<AtomIndex>,
    operations: Vec<BuildOperation>,
    validation_set: ValidationSet,
}

/// Internal tracking of operations for validation/undo
enum BuildOperation {
    AddAtom { idx: AtomIndex, atom: GraphAtom },
    AddBond { idx: BondIndex, from: AtomIndex, to: AtomIndex, bond: GraphBond },
    ModifyAtom { idx: AtomIndex, property: String },
    ModifyBond { idx: BondIndex, property: String },
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::new(),
            active_atom: None,
            operations: Vec::new(),
            validation_set: ValidationSet::standard(),
        }
    }
    
    /// Get the current active atom
    pub fn active_atom(&self) -> Option<AtomIndex> {
        self.active_atom
    }
    
    /// Set the active atom
    pub fn set_active_atom(mut self, idx: AtomIndex) -> Self {
        if self.graph.contains_node(idx.into()) {
            self.active_atom = Some(idx);
        }
        self
    }
    
    /// Start atom addition flow by specifying element
    pub fn atom<T: Into<GraphAtom>>(mut self, atom_spec: T) -> AtomContext {
        let atom = atom_spec.into();
        let idx = self.graph.add_node(atom).into();
        self.active_atom = Some(idx);
        self.operations.push(BuildOperation::AddAtom { idx, atom });
        
        AtomContext { builder: self, atom_idx: idx }
    }
    
    /// Start bond addition flow from active atom
    pub fn bond(self, order: BondOrder) -> BondContext {
        BondContext { builder: self, order, from: self.active_atom }
    }
    
    /// Start selection flow
    pub fn select(self, idx: AtomIndex) -> Result<SelectionContext, Error> {
        if !self.graph.contains_node(idx.into()) {
            return Err(Error::InvalidAtomIndex(idx));
        }
        Ok(SelectionContext { builder: self, selected: vec![idx] })
    }
    
    /// Add a molecular fragment at the active atom
    pub fn add_fragment<F: Fragment>(self, fragment: F) -> Result<Self, Error> {
        fragment.add_to_builder(self)
    }
    
    /// Validate the molecule under construction
    pub fn validate(&self) -> Result<(), Vec<Error>> {
        // Create a temporary molecule for validation
        let molecule = GraphMolecule {
            graph: self.graph.clone(),
            properties: HashMap::new(),
        };
        
        self.validation_set.validate(&molecule)
    }
    
    /// Finalize and build the molecule
    pub fn build(self) -> Result<GraphMolecule, Error> {
        // Optional validation on build
        if let Err(errors) = self.validate() {
            return Err(Error::ValidationFailed(errors));
        }
        
        Ok(GraphMolecule {
            graph: self.graph,
            properties: HashMap::new(),
        })
    }
}

/// Context after adding an atom - methods specific to a newly added atom
pub struct AtomContext {
    builder: MoleculeBuilder,
    atom_idx: AtomIndex,
}

impl AtomContext {
    /// Get the atom index for this context
    pub fn atom_idx(&self) -> AtomIndex {
        self.atom_idx
    }
    
    /// Set charge on the atom
    pub fn with_charge(mut self, charge: i8) -> Result<Self, Error> {
        if let Some(atom) = self.builder.graph.node_weight_mut(self.atom_idx.into()) {
            *atom = atom.with_charge(charge)?;
            self.builder.operations.push(BuildOperation::ModifyAtom { 
                idx: self.atom_idx, 
                property: format!("charge:{}", charge) 
            });
        }
        Ok(self)
    }
    
    /// Set unpaired electrons on the atom
    pub fn with_unpaired_electrons(mut self, unpaired: u8) -> Result<Self, Error> {
        if let Some(atom) = self.builder.graph.node_weight_mut(self.atom_idx.into()) {
            *atom = atom.with_unpaired_electrons(unpaired)?;
            self.builder.operations.push(BuildOperation::ModifyAtom { 
                idx: self.atom_idx, 
                property: format!("unpaired:{}", unpaired) 
            });
        }
        Ok(self)
    }
    
    /// Set implicit hydrogens on the atom
    pub fn with_implicit_hydrogens(mut self, count: u8) -> Self {
        if let Some(atom) = self.builder.graph.node_weight_mut(self.atom_idx.into()) {
            *atom = atom.with_implicit_hydrogens(count);
            self.builder.operations.push(BuildOperation::ModifyAtom { 
                idx: self.atom_idx, 
                property: format!("hydrogens:{}", count) 
            });
        }
        self
    }
    
    /// Bond this atom to another atom
    pub fn bond_to(self, target: AtomIndex, order: BondOrder) -> Result<BondResultContext, Error> {
        let bond = GraphBond::new(order);
        if !self.builder.graph.contains_node(target.into()) {
            return Err(Error::InvalidAtomIndex(target));
        }
        
        let bond_idx = self.builder.graph.add_edge(
            self.atom_idx.into(), 
            target.into(), 
            bond
        ).into();
        
        self.builder.operations.push(BuildOperation::AddBond { 
            idx: bond_idx, 
            from: self.atom_idx, 
            to: target, 
            bond 
        });
        
        Ok(BondResultContext { 
            builder: self.builder, 
            bond_idx 
        })
    }
    
    /// Continue chain of atoms
    pub fn attach<T: Into<GraphAtom>>(self, atom_spec: T, order: BondOrder) -> Result<AtomContext, Error> {
        let new_atom = atom_spec.into();
        let new_idx = self.builder.graph.add_node(new_atom).into();
        
        let bond = GraphBond::new(order);
        let bond_idx = self.builder.graph.add_edge(
            self.atom_idx.into(), 
            new_idx.into(), 
            bond
        ).into();
        
        let mut builder = self.builder;
        builder.operations.push(BuildOperation::AddAtom { 
            idx: new_idx, 
            atom: new_atom 
        });
        
        builder.operations.push(BuildOperation::AddBond { 
            idx: bond_idx, 
            from: self.atom_idx, 
            to: new_idx, 
            bond 
        });
        
        builder.active_atom = Some(new_idx);
        
        Ok(AtomContext { 
            builder, 
            atom_idx: new_idx 
        })
    }
    
    /// Return to the main builder with this atom as active
    pub fn done(mut self) -> MoleculeBuilder {
        self.builder.active_atom = Some(self.atom_idx);
        self.builder
    }
}

/// Context after specifying a bond
pub struct BondContext {
    builder: MoleculeBuilder,
    order: BondOrder,
    from: Option<AtomIndex>,
}

impl BondContext {
    /// Specify the source atom if not set
    pub fn from(mut self, idx: AtomIndex) -> Result<Self, Error> {
        if !self.builder.graph.contains_node(idx.into()) {
            return Err(Error::InvalidAtomIndex(idx));
        }
        self.from = Some(idx);
        Ok(self)
    }
    
    /// Specify target atom and create the bond
    pub fn to(self, idx: AtomIndex) -> Result<BondResultContext, Error> {
        if !self.builder.graph.contains_node(idx.into()) {
            return Err(Error::InvalidAtomIndex(idx));
        }
        
        if let Some(from_idx) = self.from {
            let bond = GraphBond::new(self.order);
            let bond_idx = self.builder.graph.add_edge(
                from_idx.into(), 
                idx.into(), 
                bond
            ).into();
            
            self.builder.operations.push(BuildOperation::AddBond { 
                idx: bond_idx, 
                from: from_idx, 
                to: idx, 
                bond 
            });
            
            Ok(BondResultContext { 
                builder: self.builder, 
                bond_idx 
            })
        } else {
            Err(Error::InvalidOperation("No source atom specified for bond".into()))
        }
    }
}

/// Context after creating a bond
pub struct BondResultContext {
    builder: MoleculeBuilder,
    bond_idx: BondIndex,
}

impl BondResultContext {
    /// Set bond order
    pub fn with_order(mut self, order: BondOrder) -> Self {
        if let Some(bond) = self.builder.graph.edge_weight_mut(self.bond_idx.into()) {
            *bond = bond.with_order(order);
            self.builder.operations.push(BuildOperation::ModifyBond { 
                idx: self.bond_idx, 
                property: format!("order:{}", order) 
            });
        }
        self
    }
    
    /// Return to main builder
    pub fn done(self) -> MoleculeBuilder {
        self.builder
    }
}

/// Context for selection operations
pub struct SelectionContext {
    builder: MoleculeBuilder,
    selected: Vec<AtomIndex>,
}

impl SelectionContext {
    /// Add another atom to selection
    pub fn and(mut self, idx: AtomIndex) -> Result<Self, Error> {
        if !self.builder.graph.contains_node(idx.into()) {
            return Err(Error::InvalidAtomIndex(idx));
        }
        self.selected.push(idx);
        Ok(self)
    }
    
    /// Select all atoms matching a predicate
    pub fn all<F>(mut self, predicate: F) -> Self 
    where 
        F: FnMut(&GraphAtom) -> bool 
    {
        let mut matching = self.builder.graph
            .node_indices()
            .filter_map(|idx| {
                let atom = self.builder.graph.node_weight(idx)?;
                if predicate(atom) {
                    Some(AtomIndex::from(idx))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        
        self.selected.append(&mut matching);
        self
    }
    
    /// Modify all selected atoms
    pub fn modify<F>(mut self, modifier: F) -> Result<Self, Error>
    where
        F: FnMut(&mut GraphAtom) -> Result<(), Error>
    {
        let mut modifier = modifier;
        
        for &idx in &self.selected {
            if let Some(atom) = self.builder.graph.node_weight_mut(idx.into()) {
                modifier(atom)?;
                self.builder.operations.push(BuildOperation::ModifyAtom { 
                    idx, 
                    property: "modified".to_string() 
                });
            }
        }
        
        Ok(self)
    }
    
    /// Return to builder with last selected atom active
    pub fn done(mut self) -> MoleculeBuilder {
        if let Some(&last) = self.selected.last() {
            self.builder.active_atom = Some(last);
        }
        self.builder
    }
}