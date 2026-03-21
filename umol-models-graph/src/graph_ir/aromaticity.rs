//! Aromaticity types and perception models for GraphIR.
//!
//! Core types (`AromaticSystem`, `AromaticContribution`) describe delocalized
//! π systems. Perception models (`hueckel_rule`, `hmo`, `clar`) detect aromatic
//! systems from ring topology and atom properties. `AromaticityModel` dispatches
//! to the configured model.

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

pub use self::clar::*;
pub use self::hmo::*;
pub use self::hueckel_rule::*;
use thiserror::Error;
use crate::graph_ir::config::AromaticityStrategy;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};
use crate::graph_ir::rings::{Ring, RingSet};

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    #[error("ring enumeration error: {0}")]
    RingEnumeration(String),
    #[error("hueckel input error: {0}")]
    HueckelInputError(String),
    #[error("hmo missing atom: {0}")]
    HmoMissingAtom(String),
    #[error("hmo missing parameters: {0}")]
    HmoMissingParameters(String),
    #[error("hmo invalid input: {0}")]
    HmoInvalidInput(String),
    #[error("clar input error: {0}")]
    ClarInputError(String),
}

/// Per-atom contribution to an aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticContribution {
    atom: AtomIndex,
    aromatic_valence: u8,
}

impl AromaticContribution {
    pub fn new(atom: AtomIndex, aromatic_valence: u8) -> Self {
        Self {
            atom,
            aromatic_valence,
        }
    }

    pub fn atom(&self) -> AtomIndex {
        self.atom
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }
}

// TODO: Add multiplicity field
/// An aromatic system consisting of atoms and number of electrons contributed.
/// Charge is delocalized charge, not assignable to any individual atom.
/// Each atom can participate in at most one aromatic system, appears only once
/// in the contributions list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    contributions: Vec<AromaticContribution>,
    charge: i8,
    rings: Vec<Ring>,
}

impl AromaticSystem {
    fn normalized_contributions<I>(contributions: I) -> Vec<AromaticContribution>
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        let mut contributions: Vec<AromaticContribution> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        contributions
    }

    pub fn new<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        Self::with_rings(contributions, Vec::new())
    }

    pub fn with_rings<I>(contributions: I, rings: Vec<Ring>) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        Self {
            contributions: Self::normalized_contributions(contributions),
            charge: 0,
            rings,
        }
    }

    pub fn contributions(&self) -> &[AromaticContribution] {
        &self.contributions
    }

    pub fn atom_count(&self) -> usize {
        self.contributions.len()
    }

    pub fn electron_count(&self) -> u8 {
        self.contributions.iter().map(|c| c.aromatic_valence).sum()
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = charge;
    }

    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.contributions
            .binary_search_by_key(&atom, |c| c.atom)
            .is_ok()
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.contributions.iter().map(|c| c.atom)
    }
}

pub enum AromaticityModel {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityModel {
    pub fn new(strategy: &AromaticityStrategy) -> Self {
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

    pub fn aromatic_systems(
        &self,
        builder: &MoleculeBuilder,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        match self {
            Self::HueckelRule(m) => Ok(m.find_from_rings(builder, rings)),
            Self::Hmo(m) => m.find_from_rings(builder, rings),
            Self::Clar(m) => m.find_from_rings(builder, rings),
        }
    }
}
