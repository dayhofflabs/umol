//! Aromaticity perception framework.
//!
//! `AromaticityModel` dispatches to concrete implementations in submodules:
//! `hueckel_rule`, `hmo`, `clar`.

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

use self::clar::ClarAromaticity;
use self::hmo::HmoAromaticity;
use self::hueckel_rule::HueckelRuleAromaticity;
use crate::graph_ir::aromatic::AromaticSystem;
use crate::graph_ir::config::AromaticityStrategy;
use crate::graph_ir::error::ResolutionError;
use crate::graph_ir::molecule::builder::MoleculeBuilder;
use crate::graph_ir::rings::MoleculeRings;

pub enum AromaticityModel {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityModel {
    pub fn from_strategy(strategy: &AromaticityStrategy) -> Self {
        match strategy {
            AromaticityStrategy::HueckelRule {
                element_scope,
                ring_limits,
            } => Self::HueckelRule(HueckelRuleAromaticity::new(
                element_scope.clone(),
                ring_limits.clone(),
            )),
            AromaticityStrategy::Hmo {
                element_scope,
                stabilization_threshold,
            } => Self::Hmo(HmoAromaticity::new(
                element_scope.clone(),
                *stabilization_threshold,
            )),
            AromaticityStrategy::Clar => Self::Clar(ClarAromaticity),
        }
    }

    pub fn max_ring_size(&self) -> usize {
        match self {
            Self::HueckelRule(m) => m.ring_limits.max_ring_size,
            Self::Hmo(_) => 22,
            Self::Clar(_) => 6,
        }
    }

    pub fn aromatic_systems(
        &self,
        builder: &MoleculeBuilder,
        rings: &MoleculeRings,
    ) -> Result<Vec<AromaticSystem>, ResolutionError> {
        match self {
            Self::HueckelRule(m) => Ok(m.find_from_rings(builder, rings)),
            Self::Hmo(m) => m.find_from_rings(builder, rings),
            Self::Clar(m) => Ok(m.find_from_rings(builder, rings)),
        }
    }
}
