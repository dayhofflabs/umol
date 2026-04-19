//! Top-level theory combination consumed by resolution, validation, and matching engines.

use crate::ops::aromaticity::AromaticityTheory;
use crate::ops::propagate::ValenceTheory;
use crate::ops::valence::AtomTypeRegistry;

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
