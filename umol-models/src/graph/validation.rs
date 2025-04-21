//! Validation types and functions for valence graphs

use crate::ValenceGraph;
use umol::{Error, Result};
use umol::error::ValidationError;

/// Trait for validation rules
pub trait ValidationRule {
    fn validate(&self, graph: &ValenceGraph) -> Result<(), Vec<ValidationError>>;

    fn name(&self) -> &str;

    fn description(&self) -> &str {
        self.name()
    }

}

/// Validation of atom valences
pub struct ValenceValidationRule;

impl ValidationRule for ValenceValidationRule {
    fn validate(&self, graph: &ValenceGraph) -> Result<(), Vec<ValidationError>> {
        let mut errors = Vec::new();

        for (idx, atom) in graph.atoms().enumerate() {
            let element = atom.element();
            let bond_sum = graph.atom_bonds(idx).map(|b| graph.order(b)).sum();
                
            }



// /// Trait for validation sets
// pub trait ValidationSet {
//     fn validate(&self, graph: &ValenceGraph) -> Result<(), Vec<Error>>;
// }