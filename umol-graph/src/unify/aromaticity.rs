//! Aromaticity types, perception theories, and configuration.
//!
//! Core types (`AromaticSystem`, `AromaticContribution`) describe delocalized
//! π systems. Perception implementations (`hueckel_rule`, `hmo`, `clar`) detect
//! aromatic systems from ring topology and atom properties. `AromaticityTheory`
//! dispatches to the configured perception.

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
use crate::ast::molecule::{AromaticSystemAst, MoleculeAst};
use crate::ast::rings::{RingEnumerationStrategy, RingEnumerator, RingFamily, RingSet};
use crate::unify::resolve::Progress;

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

/// Policy for mismatches between aromatic hints and detected aromatic systems.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityHintPolicy {
    Strict,
    Ignore,
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
}

impl AromaticSystem {
    pub fn new<I>(contributions: I) -> Self
    where
        I: IntoIterator<Item = AromaticContribution>,
    {
        let mut contributions: Vec<AromaticContribution> = contributions.into_iter().collect();
        contributions.sort_unstable();
        contributions.dedup_by_key(|c| c.atom);
        Self {
            contributions,
            charge: 0,
            spin: SpinState::closed_shell(),
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

    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.contributions
            .binary_search_by_key(&atom, |c| c.atom)
            .is_ok()
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.contributions.iter().map(|c| c.atom)
    }
}

#[derive(Clone, Debug)]
pub struct AromaticityTheory {
    pub kind: AromaticityKind,
    pub ring_enumeration: RingEnumerationStrategy,
}

#[derive(Clone, Debug)]
pub enum AromaticityKind {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityTheory {
    /// Daylight (SMILES) aromaticity: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self {
            kind: AromaticityKind::HueckelRule(HueckelRuleAromaticity::new(
                ElementScope::AllowList(vec![
                    Element::C,
                    Element::N,
                    Element::O,
                    Element::S,
                    Element::Se,
                    Element::As,
                ]),
                RingLimits::default(),
            )),
            ring_enumeration: RingEnumerationStrategy::default(),
        }
    }

    /// MDL (MOL/SDF) aromaticity: C and N only. Minimum ring size 6.
    pub fn mdl() -> Self {
        Self {
            kind: AromaticityKind::HueckelRule(HueckelRuleAromaticity::new(
                ElementScope::AllowList(vec![Element::C, Element::N]),
                RingLimits {
                    min_ring_size: 6,
                    ..RingLimits::default()
                },
            )),
            ring_enumeration: RingEnumerationStrategy::default(),
        }
    }

    /// Permissive aromaticity: any element with aromatic valence states.
    pub fn permissive() -> Self {
        Self {
            kind: AromaticityKind::HueckelRule(HueckelRuleAromaticity::new(
                ElementScope::Any,
                RingLimits::default(),
            )),
            ring_enumeration: RingEnumerationStrategy::default(),
        }
    }

    pub fn aromatic_systems(
        &self,
        ast: &MoleculeAst,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let mut systems = match &self.kind {
            AromaticityKind::HueckelRule(m) => Ok(m.find_from_rings(ast, rings)),
            AromaticityKind::Hmo(m) => m.find_from_rings(ast, rings),
            AromaticityKind::Clar(m) => m.find_from_rings(ast, rings),
        }?;
        systems.sort_by(|a, b| {
            let min_a = a.atoms().min();
            let min_b = b.atoms().min();
            min_a.cmp(&min_b)
        });
        Ok(systems)
    }

    pub fn refine(&self, ast: &mut MoleculeAst) -> Result<Progress, AromaticityError> {
        let ring_family = match &self.kind {
            AromaticityKind::Clar(_) => RingFamily::InducedBenzenoid,
            AromaticityKind::HueckelRule(_) | AromaticityKind::Hmo(_) => RingFamily::Simple,
        };
        let enumerator = RingEnumerator::new(ring_family, &self.ring_enumeration);
        let rings = enumerator.enumerate(ast);
        let systems = self.aromatic_systems(ast, &rings)?;
        if systems.is_empty() {
            return Ok(Progress::Fixpoint);
        }
        let mut builder = ast.edit();
        for sys in &systems {
            let atoms: Vec<AtomIdx> = sys.atoms().collect();
            builder.add_aromatic_system(atoms, AromaticSystemAst {});
        }
        *ast = builder.build();
        Ok(Progress::Advanced)
    }
}
