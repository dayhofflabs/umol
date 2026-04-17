//! Top-level theory combination consumed by resolution, validation, and matching engines.

use crate::unify::aromaticity::AromaticityTheory;
use crate::unify::propagate::ValenceTheory;
use crate::unify::valence::AtomTypeRegistry;

#[derive(Clone, Debug)]
pub struct Chemistry {
    pub valence: ValenceTheory,
    pub aromaticity: AromaticityTheory,
}

impl Default for Chemistry {
    fn default() -> Self {
        Self {
            valence: ValenceTheory::AtomTyping {
                registry: AtomTypeRegistry::default_registry().clone(),
            },
            aromaticity: AromaticityTheory::daylight(),
        }
    }
}
