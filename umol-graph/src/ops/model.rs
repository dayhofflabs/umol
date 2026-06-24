//! Chemistry model: top-level configuration consumed by the resolver and
//! validator engines.
//!
//! `ChemistryModel` wraps a `ValenceModel` (atom-typing or counts) and an
//! `AromaticityModel` (Hückel rule, HMO, or Clar). Engines read this; engines
//! and configs are kept as distinct types so multiple engine instances can
//! share one model. The model carries no resolution behavior of its own — the
//! resolver and validator engines do the work.

use std::array;
use std::borrow::Cow;

use strum::EnumCount;
use thiserror::Error;
use umol_ast::ast::{ConstitutionColoring, GraphSymmetryConfig, StereoKind};
use umol_chem::element::Element;

use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

#[derive(Debug, Clone)]
pub struct ChemistryModel {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
    pub stereo: StereoModel,
}

impl Default for ChemistryModel {
    fn default() -> Self {
        Self {
            valence: ValenceModel::AtomTyping(AtomTypingModel {
                registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
            }),
            aromaticity: AromaticityModel::daylight(),
            stereo: StereoModel::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum ValenceModel {
    AtomTyping(AtomTypingModel),
    Counts(CountsModel),
}

/// Atom-typing valence model: the registry of `AtomAst` patterns.
#[derive(Debug, Clone)]
pub struct AtomTypingModel {
    pub registry: Cow<'static, AtomTypeRegistry>,
}

/// Counts valence model: the per-element covalence table.
#[derive(Debug, Clone)]
pub struct CountsModel {
    pub table: Cow<'static, ValenceTable>,
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

/// Stereo perception model. `kind_models` is a per-`StereoKind` slot map (indexed
/// by the kind's discriminant); a `None` slot means that kind is not perceived.
/// `para_stereo` enables the graph-symmetry fixpoint iteration that resolves
/// para-stereocenters; `inconsistency` governs how the resolver handles a
/// `#T`/`#C` assertion it cannot realize.
#[derive(Debug, Clone)]
pub struct StereoModel {
    pub kind_models: [Option<StereoKindModel>; StereoKind::COUNT],
    pub para_stereo: bool,
    pub max_iterations: usize,
    pub inconsistency: InconsistencyPolicy,
}

/// Per-kind perception settings: the elements eligible to bear this kind and
/// whether the kind's sites are treated as fluxional.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoKindModel {
    pub scope: ElementScope,
    pub fluxionality: bool,
}

/// How the stereo resolver handles a `#T`/`#C` assertion it cannot (fully)
/// realize: keep what it can, strip the unrealizable element, or error. Never a
/// silent drop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InconsistencyPolicy {
    Keep,
    Strip,
    Error,
}

impl StereoModel {
    /// The per-kind model for `kind`, or `None` if the kind is not perceived.
    pub fn kind_model(&self, kind: StereoKind) -> Option<&StereoKindModel> {
        self.kind_models[kind as usize].as_ref()
    }

    /// Build the umol-ast graph-symmetry config the validator runs: full
    /// constitution coloring, fixpoint iteration gated by `para_stereo`.
    pub fn graph_symmetry_config(&self) -> GraphSymmetryConfig<ConstitutionColoring> {
        GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: self.para_stereo,
            max_iterations: self.max_iterations,
        }
    }
}

impl Default for StereoModel {
    /// Perceive the two realized binary kinds — tetrahedral atoms and cis/trans
    /// bonds — for any element; the higher geometries (square-planar,
    /// trigonal-bipyramidal, octahedral, axial) are staged off. No para-stereo
    /// fixpoint by default; inconsistency is an error.
    fn default() -> Self {
        let mut kind_models: [Option<StereoKindModel>; StereoKind::COUNT] =
            array::from_fn(|_| None);
        kind_models[StereoKind::Tetrahedral as usize] = Some(StereoKindModel {
            scope: ElementScope::Any,
            fluxionality: false,
        });
        kind_models[StereoKind::CisTrans as usize] = Some(StereoKindModel {
            scope: ElementScope::Any,
            fluxionality: false,
        });
        Self {
            kind_models,
            para_stereo: false,
            max_iterations: 16,
            inconsistency: InconsistencyPolicy::Error,
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
    use umol_chem::element::Element;

    use super::*;

    #[test]
    fn test_chemistry_model_default() {
        let model = ChemistryModel::default();
        assert!(matches!(model.valence, ValenceModel::AtomTyping(_)));
        assert!(matches!(
            model.aromaticity,
            AromaticityModel::HueckelRule { .. }
        ));
        assert!(!model.stereo.para_stereo);
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

    #[test]
    fn test_stereo_model_default() {
        let model = StereoModel::default();
        assert!(!model.para_stereo);
        assert_eq!(model.max_iterations, 16);
        assert_eq!(model.inconsistency, InconsistencyPolicy::Error);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, Some(StereoKindModel { scope: ElementScope::Any, fluxionality: false }))]
    #[case::cis_trans(StereoKind::CisTrans, Some(StereoKindModel { scope: ElementScope::Any, fluxionality: false }))]
    #[case::axial(StereoKind::Axial, None)]
    #[case::square_planar(StereoKind::SquarePlanar, None)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, None)]
    #[case::octahedral(StereoKind::Octahedral, None)]
    fn test_stereo_model_kind_model(
        #[case] kind: StereoKind,
        #[case] expected: Option<StereoKindModel>,
    ) {
        assert_eq!(StereoModel::default().kind_model(kind), expected.as_ref());
    }

    #[rstest]
    #[case::no_para(false, false)]
    #[case::para(true, true)]
    fn test_stereo_model_graph_symmetry_config(
        #[case] para_stereo: bool,
        #[case] expected_fixpoint: bool,
    ) {
        let model = StereoModel {
            para_stereo,
            ..StereoModel::default()
        };
        let cfg = model.graph_symmetry_config();
        assert_eq!(cfg.iterate_to_fixpoint, expected_fixpoint);
        assert_eq!(cfg.max_iterations, 16);
    }
}
