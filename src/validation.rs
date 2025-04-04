// Validation rules and infrastructure for molecules

use crate::{Error, Element};
use crate::graph::{GraphMolecule, GraphAtom, GraphBond, AtomIndex, BondIndex};
use std::collections::HashSet;

/// Trait for any validation rule
pub trait ValidationRule {
    /// Check if the molecule satisfies this rule
    fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>>;
    
    /// Get a description of this rule
    fn description(&self) -> &str;
}

/// Validation for atom valence based on element properties
pub struct ValenceRule;

impl ValidationRule for ValenceRule {
    fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>> {
        let mut errors = Vec::new();
        
        for (idx, atom) in molecule.atoms() {
            if let Some(element) = atom.element() {
                let neighbors = molecule.neighbors(idx);
                let bond_sum: u8 = neighbors.iter()
                    .filter_map(|(_, bond_idx)| {
                        let bond = molecule.bond(*bond_idx)?;
                        Some(bond.order().value())
                    })
                    .sum();
                    
                let h_count = atom.implicit_hydrogens();
                let total_bonds = bond_sum + h_count;
                
                // Check against element's allowed valence
                // This is a simplified check - real validation would account for charge
                if total_bonds > element.max_unpaired_electrons() as u8 {
                    errors.push(Error::InvalidValence {
                        element,
                        atom_idx: idx,
                        expected: element.max_unpaired_electrons() as u8,
                        actual: total_bonds,
                    });
                }
            }
        }
        
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    fn description(&self) -> &str {
        "Checks that all atoms have valid valence counts"
    }
}

/// Validation for connectivity (no disconnected fragments)
pub struct ConnectivityRule;

impl ValidationRule for ConnectivityRule {
    fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>> {
        if molecule.atom_count() <= 1 {
            return Ok(());
        }
        
        // Use a breadth-first search to check connectivity
        let mut visited = HashSet::new();
        let start_idx = molecule.atoms().next().map(|(idx, _)| idx).unwrap();
        
        let mut queue = vec![start_idx];
        visited.insert(start_idx);
        
        while let Some(idx) = queue.pop() {
            for (neighbor, _) in molecule.neighbors(idx) {
                if !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }
        
        if visited.len() != molecule.atom_count() {
            Err(vec![Error::DisconnectedMolecule])
        } else {
            Ok(())
        }
    }
    
    fn description(&self) -> &str {
        "Checks that the molecule is fully connected"
    }
}

/// Validation for chemical reasonableness 
pub struct ChemicalReasonablenessRule;

impl ValidationRule for ChemicalReasonablenessRule {
    fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>> {
        let mut errors = Vec::new();
        
        // Check for atoms with multiple bonds to the same neighbor
        for (idx, _) in molecule.atoms() {
            let neighbors = molecule.neighbors(idx);
            let mut neighbor_counts = HashMap::new();
            
            for (neighbor, _) in neighbors {
                *neighbor_counts.entry(neighbor).or_insert(0) += 1;
            }
            
            for (neighbor, count) in neighbor_counts {
                if count > 1 {
                    errors.push(Error::UnreasonableStructure(
                        format!("Atom {:?} has {} bonds to atom {:?}", idx, count, neighbor)
                    ));
                }
            }
        }
        
        // Implement other reasonableness checks:
        // - No small rings (< 3 atoms)
        // - No strained geometries
        // - etc.
        
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    fn description(&self) -> &str {
        "Checks that the molecule follows basic chemical reasonableness rules"
    }
}

/// Collection of validation rules
pub struct ValidationSet {
    rules: Vec<Box<dyn ValidationRule>>,
}

impl ValidationSet {
    pub fn new() -> Self {
        Self { rules: Vec::new() }
    }
    
    pub fn add_rule<R: ValidationRule + 'static>(&mut self, rule: R) {
        self.rules.push(Box::new(rule));
    }
    
    pub fn validate(&self, molecule: &GraphMolecule) -> Result<(), Vec<Error>> {
        let mut all_errors = Vec::new();
        
        for rule in &self.rules {
            if let Err(errors) = rule.validate(molecule) {
                all_errors.extend(errors);
            }
        }
        
        if all_errors.is_empty() {
            Ok(())
        } else {
            let error_strings = all_errors.iter()
                .map(|e| format!("{}", e))
                .collect();
            Err(vec![Error::ValidationFailed(error_strings)])
        }
    }
    
    /// Create a standard validation set with common rules
    pub fn standard() -> Self {
        let mut set = Self::new();
        set.add_rule(ValenceRule);
        set.add_rule(ConnectivityRule);
        set.add_rule(ChemicalReasonablenessRule);
        set
    }
}