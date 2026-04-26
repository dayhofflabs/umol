//! Chemistry model: top-level configuration consumed by the resolver and
//! validator engines.
//!
//! `ChemistryModel` wraps a `ValenceModel` (atom-typing or counts) and an
//! `AromaticityModel` (Hückel rule, HMO, or Clar). Engines read this; engines
//! and configs are kept as distinct types so multiple engine instances can
//! share one model.

use thiserror::Error;
use umol_shared::element::Element;

use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

#[derive(Debug, Clone)]
pub struct ChemistryModel {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
}

impl Default for ChemistryModel {
    fn default() -> Self {
        Self {
            valence: ValenceModel::AtomTyping {
                registry: AtomTypeRegistry::default_registry().clone(),
            },
            aromaticity: AromaticityModel::daylight(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum ValenceModel {
    AtomTyping {
        registry: AtomTypeRegistry,
    },
    Counts {
        table: ValenceTable,
        allow_implicit_hydrogens: bool,
    },
}

#[derive(Debug, Clone)]
pub enum AromaticityModel {
    HueckelRule {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
    Hmo {
        scope: ElementScope,
        stabilization_threshold: f64,
    },
    Clar {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
}

impl AromaticityModel {
    /// Daylight (SMILES) aromaticity scope: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self::HueckelRule {
            scope: ElementScope::AllowList(vec![
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

    /// MDL (MOL/SDF) aromaticity scope: C and N only, minimum ring size 6.
    pub fn mdl() -> Self {
        Self::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        }
    }

    /// Permissive aromaticity scope: any element.
    pub fn permissive() -> Self {
        Self::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        }
    }
}

/// Elements eligible for aromaticity perception.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementScope {
    Any,
    AllowList(Vec<Element>),
}

impl ElementScope {
    pub fn contains(&self, element: Element) -> bool {
        match self {
            Self::Any => true,
            Self::AllowList(list) => list.contains(&element),
        }
    }
}

/// Ring-size and fused-ring search bounds for ring-based aromaticity perception.
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

/// Setup-time errors loading model data (TOML registries / valence tables).
/// Distinct from the per-engine `*Error` types that surface at resolve time.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConfigError {
    #[error("invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("invalid valence table: {0}")]
    InvalidValenceTable(String),
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_shared::element::Element;

    use super::*;

    #[test]
    fn test_chemistry_model_default() {
        let model = ChemistryModel::default();
        assert!(matches!(model.valence, ValenceModel::AtomTyping { .. }));
        assert!(matches!(model.aromaticity, AromaticityModel::HueckelRule { .. }));
    }

    #[rstest]
    #[case::any(ElementScope::Any, Element::U, true)]
    #[case::allow_match(ElementScope::AllowList(vec![Element::C]), Element::C, true)]
    #[case::allow_miss(ElementScope::AllowList(vec![Element::C]), Element::N, false)]
    fn test_element_scope_contains(
        #[case] scope: ElementScope,
        #[case] element: Element,
        #[case] expected: bool,
    ) {
        assert_eq!(scope.contains(element), expected);
    }

    #[test]
    fn test_aromaticity_model_daylight_scope() {
        match AromaticityModel::daylight() {
            AromaticityModel::HueckelRule { scope, .. } => {
                assert!(scope.contains(Element::C));
                assert!(scope.contains(Element::N));
                assert!(!scope.contains(Element::B));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }

    #[test]
    fn test_aromaticity_model_mdl_min_ring_size() {
        match AromaticityModel::mdl() {
            AromaticityModel::HueckelRule {
                ring_limits, scope, ..
            } => {
                assert_eq!(ring_limits.min_ring_size, 6);
                assert!(scope.contains(Element::N));
                assert!(!scope.contains(Element::O));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }

    #[test]
    fn test_aromaticity_model_permissive_scope() {
        match AromaticityModel::permissive() {
            AromaticityModel::HueckelRule { scope, .. } => {
                assert!(matches!(scope, ElementScope::Any));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }
}
