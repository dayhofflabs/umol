//! Aromaticity perception theories and configuration.
//!
//! Perception implementations (`hueckel_rule`, `hmo`, `clar`) detect delocalized
//! π systems from ring topology and atom properties, producing
//! [`AromaticSystem`](crate::ast::aromatic::AromaticSystem) values lifted into
//! the molecule AST. `AromaticityTheory` dispatches to the configured perception.

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

pub use clar::*;
pub use hmo::*;
pub use hueckel_rule::*;
use thiserror::Error;
use umol_shared::element::Element;

use crate::ast::constraint::MoleculeConstraint;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{RingEnumerationStrategy, RingSet};
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

#[derive(Clone, Debug)]
pub struct AromaticityTheory {
    pub strategy: AromaticityStrategy,
    pub ring_enumeration: RingEnumerationStrategy,
}

#[derive(Clone, Debug)]
pub enum AromaticityStrategy {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityTheory {
    /// Daylight (SMILES) aromaticity: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self {
            strategy: AromaticityStrategy::HueckelRule(HueckelRuleAromaticity::new(
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
    /// TODO: Verify
    pub fn mdl() -> Self {
        Self {
            strategy: AromaticityStrategy::HueckelRule(HueckelRuleAromaticity::new(
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
            strategy: AromaticityStrategy::HueckelRule(HueckelRuleAromaticity::new(
                ElementScope::Any,
                RingLimits::default(),
            )),
            ring_enumeration: RingEnumerationStrategy::default(),
        }
    }

    pub fn refine(
        &self,
        ast: &mut MoleculeAst,
        rings: &RingSet,
    ) -> Result<Progress, AromaticityError> {
        let mut systems = match &self.strategy {
            AromaticityStrategy::HueckelRule(m) => Ok(m.find_from_rings(ast, rings)),
            AromaticityStrategy::Hmo(m) => m.find_from_rings(ast, rings),
            AromaticityStrategy::Clar(m) => m.find_from_rings(ast, rings),
        }?;
        if systems.is_empty() {
            return Ok(Progress::Fixpoint);
        }
        systems.sort_by(|a, b| a.atoms().first().cmp(&b.atoms().first()));

        let mut builder = ast.edit();
        for sys in &systems {
            builder.add_aromatic_system(sys.atoms().to_vec(), sys.ast().clone());
            for (idx, c) in sys.atoms().iter().zip(sys.atom_constraints()) {
                builder.push_constraint(MoleculeConstraint::AtomPred(*idx, c.clone()));
            }
            for (idx, c) in sys.bonds().iter().zip(sys.bond_constraints()) {
                builder.push_constraint(MoleculeConstraint::BondPred(*idx, c.clone()));
            }
        }
        *ast = builder.build();
        Ok(Progress::Advanced)
    }
}
