//! Molecule transforms.
//!
//! Transforms are applied on `MoleculeBuilder` and composed explicitly with
//! conversion (`MoleculeBuilder::from_molecule`) and finalization (`build`).

use thiserror::Error;

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

        builder.clear_aromatic_systems();
        let systems = model.aromatic_systems(builder, &rings)?;
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
                    bond.set_order(order);
                }
            }
        }

        builder.clear_aromatic_systems();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use smallvec::SmallVec;
    use umol_data::Element;

    use super::*;
    use crate::graph_ir::aromaticity::{AromaticContribution, AromaticSystem};
    use crate::graph_ir::atom::AtomBuilder;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::config::ResolveConfig;
    use crate::graph_ir::molecule::Molecule;
    use crate::spec;

    fn carbon_aromatic_1() -> AtomBuilder {
        let spec = spec!("{Cv2a1H}");
        let mut ab = AtomBuilder::new(Element::C);
        ab.set_candidates(SmallVec::from_elem(spec, 1));
        ab
    }

    fn benzene_kekule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        builder
            .build(&ResolveConfig::default())
            .expect("benzene should build")
    }

    fn naphthalene_kekule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        let ring1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (a, b) in ring1_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let ring2_edges = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
        for (a, b) in ring2_edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
            .build(&ResolveConfig::default())
            .expect("naphthalene should build")
    }

    #[test]
    fn no_op_sequence() {
        let mol = benzene_kekule();
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        let seq = TransformSequence::new();
        seq.apply(&mut builder).expect("no-op should succeed");
        let out = builder
            .build(&ResolveConfig::default())
            .expect("build should succeed");
        assert_eq!(out.atom_count(), mol.atom_count());
        assert_eq!(out.bond_count(), mol.bond_count());
        assert_eq!(out.aromatic_system_count(), mol.aromatic_system_count());
    }

    #[test]
    fn ordering_aromatize_then_kekulize() {
        let mol = benzene_kekule();
        let seq = TransformSequence::from_configs(vec![
            TransformConfig::Aromatize(AromatizeConfig {
                aromaticity_strategy: AromaticityStrategy::daylight(),
                enumeration_strategy: RingEnumerationStrategy::default(),
            }),
            TransformConfig::Kekulize(KekulizeConfig::default()),
        ]);
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        seq.apply(&mut builder).expect("sequence works");
        let out = builder
            .build(&ResolveConfig::default())
            .expect("build should succeed");
        assert_eq!(out.aromatic_system_count(), 0);
        let doubles = out.bonds().filter(|b| b.order() == 2).count();
        assert_eq!(doubles, 3);
    }

    #[test]
    fn ordering_kekulize_then_aromatize() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        builder.add_aromatic_system(AromaticSystem::new(
            atoms
                .iter()
                .copied()
                .map(|a| AromaticContribution::new(a, 1)),
        ));
        let mol = builder
            .build(&ResolveConfig::default())
            .expect("aromatic benzene should build");

        let seq = TransformSequence::from_configs(vec![
            TransformConfig::Kekulize(KekulizeConfig::default()),
            TransformConfig::Aromatize(AromatizeConfig {
                aromaticity_strategy: AromaticityStrategy::daylight(),
                enumeration_strategy: RingEnumerationStrategy::default(),
            }),
        ]);
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        seq.apply(&mut builder).expect("sequence works");
        let out = builder
            .build(&ResolveConfig::default())
            .expect("build should succeed");
        assert_eq!(out.aromatic_system_count(), 1);
    }

    #[test]
    fn aromatize_detects_naphthalene() {
        let mol = naphthalene_kekule();
        let seq = TransformSequence::from_configs(vec![TransformConfig::Aromatize(
            AromatizeConfig::default(),
        )]);
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        seq.apply(&mut builder).expect("aromatization works");
        let out = builder
            .build(&ResolveConfig::default())
            .expect("build should succeed");
        assert_eq!(out.aromatic_system_count(), 1);
        let system = out
            .aromatic_system(crate::graph_ir::molecule::AromaticSystemIndex(0))
            .expect("system exists");
        assert_eq!(system.atom_count(), 10);
    }

    #[test]
    fn kekulize_naphthalene_localizes_bonds() {
        let seq = TransformSequence::from_configs(vec![
            TransformConfig::Aromatize(AromatizeConfig::default()),
            TransformConfig::Kekulize(KekulizeConfig::default()),
        ]);
        let mut builder = MoleculeBuilder::from_molecule(&naphthalene_kekule());
        seq.apply(&mut builder).expect("kekulization works");
        let out = builder
            .build(&ResolveConfig::default())
            .expect("build should succeed");
        assert_eq!(out.aromatic_system_count(), 0);
        let doubles = out.bonds().filter(|b| b.order() == 2).count();
        assert_eq!(doubles, 5);
    }

    #[test]
    fn kekulize_failure_is_deterministic() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..3)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        for i in 0..3 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 3], BondBuilder::new(1, None));
        }
        builder.add_aromatic_system(AromaticSystem::new(
            atoms
                .iter()
                .copied()
                .map(|a| AromaticContribution::new(a, 1)),
        ));
        let mol = builder
            .build(&ResolveConfig::default())
            .expect("triangle should build");
        let seq = TransformSequence::from_configs(vec![TransformConfig::Kekulize(
            KekulizeConfig::default(),
        )]);
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        let err = seq.apply(&mut builder).expect_err("kekulize should fail");
        assert!(err
            .to_string()
            .contains("kekulization failed for aromatic system 0"));
    }

    #[test]
    fn add_executes_in_order() {
        struct AddSystem;
        impl Transform for AddSystem {
            fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError> {
                let atoms: Vec<_> = builder.atom_indices().collect();
                if !atoms.is_empty() {
                    builder.add_aromatic_system(AromaticSystem::new(
                        atoms.into_iter().map(|a| AromaticContribution::new(a, 1)),
                    ));
                }
                Ok(())
            }
        }

        struct ClearSystems;
        impl Transform for ClearSystems {
            fn apply(&self, builder: &mut MoleculeBuilder) -> Result<(), TransformError> {
                builder.clear_aromatic_systems();
                Ok(())
            }
        }

        let mol = benzene_kekule();
        let mut seq = TransformSequence::new();
        seq.add(AddSystem);
        seq.add(ClearSystems);

        let mut builder = MoleculeBuilder::from_molecule(&mol);
        seq.apply(&mut builder).expect("sequence should succeed");
        assert_eq!(builder.aromatic_system_count(), 0);
    }
}
