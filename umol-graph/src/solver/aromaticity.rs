//! Aromaticity types, perception models, and configuration.
//!
//! Core types (`AromaticSystem`, `AromaticContribution`) describe delocalized
//! π systems. Perception models (`hueckel_rule`, `hmo`, `clar`) detect aromatic
//! systems from ring topology and atom properties. `AromaticityModel` dispatches
//! to the configured model.

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

pub use clar::*;
pub use hmo::*;
pub use hueckel_rule::*;
use thiserror::Error;
use umol_shared::element::Element;
use umol_shared::spin::SpinState;

use crate::ast::AtomIdx;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{Ring, RingSet};

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

/// Elements eligible for aromaticity perception.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementScope {
    Any,
    AllowList(Vec<Element>),
}

/// Ring size and fused-ring constraints for HueckelRule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingLimits {
    pub min_ring_size: usize,
    pub max_ring_size: usize,
    pub include_fused: bool,
    pub max_fused_combination: usize,
    pub max_fused_search: usize,
}

impl Default for RingLimits {
    fn default() -> Self {
        Self {
            min_ring_size: 3,
            max_ring_size: 22,
            include_fused: true,
            max_fused_combination: 6,
            max_fused_search: 10_000,
        }
    }
}

#[derive(Clone, Debug)]
pub enum AromaticityStrategy {
    HueckelRule {
        element_scope: ElementScope,
        ring_limits: RingLimits,
    },
    Hmo {
        element_scope: ElementScope,
        /// Delocalization energy per pi-electron (in units of |beta|) required
        /// for classification as aromatic. Benzene: dE/n ~ 0.33|beta|.
        stabilization_threshold: f64,
    },
    Clar,
}

/// Policy for mismatches between aromatic hints and detected aromatic systems.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityHintPolicy {
    Strict,
    Ignore,
}

impl AromaticityStrategy {
    /// Daylight (SMILES) aromaticity: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::AllowList(vec![
                Element::C,
                Element::N,
                Element::O,
                Element::S,
                Element::Se,
                Element::As,
            ]),
            ring_limits: RingLimits::default(),
        }
    }

    /// MDL (MOL/SDF) aromaticity: C and N only. Minimum ring size 6.
    pub fn mdl() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        }
    }

    /// Permissive aromaticity: any element with aromatic valence states.
    pub fn permissive() -> Self {
        Self::HueckelRule {
            element_scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        }
    }
}

/// Per-atom contribution to an aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticContribution {
    atom: AtomIdx,
    aromatic_valence: u8,
}

impl AromaticContribution {
    pub fn new(atom: AtomIdx, aromatic_valence: u8) -> Self {
        Self {
            atom,
            aromatic_valence,
        }
    }

    pub fn atom(&self) -> AtomIdx {
        self.atom
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }
}

/// An aromatic system consisting of atoms and number of electrons contributed.
/// Charge is delocalized charge, not assignable to any individual atom.
/// Each atom can participate in at most one aromatic system, appears only once
/// in the contributions list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AromaticSystem {
    contributions: Vec<AromaticContribution>,
    charge: i8,
    spin: SpinState,
    // TODO: Check if this computed property should be removed.
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
            spin: SpinState::closed_shell(),
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

    pub fn spin(&self) -> SpinState {
        self.spin
    }

    pub fn set_spin(&mut self, spin: SpinState) {
        self.spin = spin;
    }

    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.contributions
            .binary_search_by_key(&atom, |c| c.atom)
            .is_ok()
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
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
        ast: &MoleculeAst,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let mut systems = match self {
            Self::HueckelRule(m) => Ok(m.find_from_rings(ast, rings)),
            Self::Hmo(m) => m.find_from_rings(ast, rings),
            Self::Clar(m) => m.find_from_rings(ast, rings),
        }?;
        systems.sort_by(|a, b| {
            let min_a = a.atoms().min();
            let min_b = b.atoms().min();
            min_a.cmp(&min_b)
        });
        Ok(systems)
    }
}
