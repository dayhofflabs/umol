//! Molecule transforms.
//!
//! Transforms are applied on `MoleculeBuilder` and composed explicitly with
//! conversion (`MoleculeBuilder::from_molecule`) and finalization (`build`).

use std::collections::BTreeSet;

use thiserror::Error;

use crate::graph_ir::atom_pattern::Pattern;
use crate::graph_ir::aromaticity::{AromaticityError, AromaticityModel};
use crate::graph_ir::config::{AromaticityStrategy, RingEnumerationStrategy};
use crate::graph_ir::kekule::{kekulize, KekuleConfig, KekulizationError};
use crate::graph_ir::molecule::MoleculeBuilder;
use crate::graph_ir::rings::{RingEnumerator, RingFamily};

#[derive(Debug, Error)]
pub enum TransformError {
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Kekulization(#[from] KekulizationError),
    #[error("kekulization failed for aromatic system {system_index}: {source}")]
    KekulizationAtSystem {
        system_index: usize,
        #[source]
        source: KekulizationError,
    },
    #[error("build failed after transform sequence: {0}")]
    Build(String),
}

pub trait Transform {
    fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError>;
}

#[derive(Clone, Debug)]
pub enum TransformConfig {
    Aromatize(AromatizeConfig),
    Kekulize(KekulizeConfig),
}

#[derive(Default)]
pub struct TransformSequence {
    steps: Vec<Box<dyn Transform>>,
}

impl TransformSequence {
    pub fn new() -> Self {
        Self { steps: Vec::new() }
    }

    pub fn from_configs(configs: Vec<TransformConfig>) -> Self {
        let mut sequence = Self::new();
        for config in configs {
            match config {
                TransformConfig::Aromatize(cfg) => sequence.add(Aromatize::new(cfg)),
                TransformConfig::Kekulize(cfg) => sequence.add(Kekulize::new(cfg)),
            }
        }
        sequence
    }

    pub fn add<T: Transform + 'static>(&mut self, transform: T) {
        self.steps.push(Box::new(transform));
    }

    pub fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError> {
        for step in &self.steps {
            step.apply(builder)?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct AromatizeConfig {
    pub aromaticity_strategy: AromaticityStrategy,
    pub enumeration_strategy: RingEnumerationStrategy,
}

impl Default for AromatizeConfig {
    fn default() -> Self {
        Self {
            aromaticity_strategy: AromaticityStrategy::daylight(),
            enumeration_strategy: RingEnumerationStrategy::default(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Aromatize {
    config: AromatizeConfig,
}

impl Aromatize {
    pub fn new(config: AromatizeConfig) -> Self {
        Self { config }
    }
}

impl Transform for Aromatize {
    fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError> {
        let ring_family = match self.config.aromaticity_strategy {
            AromaticityStrategy::Clar => RingFamily::InducedBenzenoid,
            AromaticityStrategy::HueckelRule { .. } | AromaticityStrategy::Hmo { .. } => {
                RingFamily::Simple
            }
        };
        let model = AromaticityModel::new(&self.config.aromaticity_strategy);
        let enumerator = RingEnumerator::new(ring_family, &self.config.enumeration_strategy);
        let rings = enumerator.enumerate_builder(builder);

        let systems = model.aromatic_systems(builder, &rings)?;
        let aromatic_bonds: BTreeSet<_> = systems
            .iter()
            .flat_map(|system| system.rings().iter())
            .flat_map(|ring| ring.bonds().iter().copied())
            .collect();
        for bond_idx in aromatic_bonds {
            builder.bond_mut(bond_idx).unwrap().order = Pattern::Is(1);
            builder.set_bond_aromatic_hint(bond_idx, true);
        }

        builder.clear_aromatic_systems();
        for system in systems {
            builder.add_aromatic_system(system);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KekulizeAlgorithm {
    Dfs,
    Matching,
}

#[derive(Clone, Debug)]
pub struct KekulizeConfig {
    pub algorithm: KekulizeAlgorithm,
    pub max_backtrack_steps: usize,
}

impl Default for KekulizeConfig {
    fn default() -> Self {
        Self {
            algorithm: KekulizeAlgorithm::Dfs,
            max_backtrack_steps: 100_000,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Kekulize {
    config: KekulizeConfig,
}

impl Kekulize {
    pub fn new(config: KekulizeConfig) -> Self {
        Self { config }
    }
}

impl Transform for Kekulize {
    fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError> {
        let systems: Vec<_> = builder.aromatic_systems().cloned().collect();
        for (idx, system) in systems.iter().enumerate() {
            let assignment = match self.config.algorithm {
                KekulizeAlgorithm::Dfs => kekulize(
                    builder,
                    system,
                    &KekuleConfig {
                        max_backtrack_steps: self.config.max_backtrack_steps,
                        bond_order_hints: None,
                    },
                ),
                KekulizeAlgorithm::Matching => Err(KekulizationError::UnsupportedAlgorithm(
                    "matching".to_string(),
                )),
            }
            .map_err(|e| TransformError::KekulizationAtSystem {
                system_index: idx,
                source: e,
            })?;

            for (bond_idx, order) in assignment.bond_orders {
                if let Some(bond) = builder.bond_mut(bond_idx) {
                    bond.order = Pattern::Is(order);
                }
            }
        }

        builder.clear_aromatic_systems();
        Ok(())
    }
}
