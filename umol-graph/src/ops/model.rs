//! Chemistry model: top-level configuration consumed by the resolver and
//! validator engines.
//!
//! `ChemistryModel` combines a `ValenceModel` (atom-typing or counts), an
//! `AromaticityModel` (Hückel rule, HMO, or Clar), and a `StereoModel`. Engines
//! read this; engines and configs are kept as distinct types so multiple engine
//! instances can share one model. The model carries no resolution behavior of
//! its own — the resolver and validator engines do the work.

use std::array;
use std::borrow::Cow;

use strum::EnumCount;
use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::StereoKind;

use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

#[derive(Debug, Clone, PartialEq)]
pub struct ChemistryModel {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
    pub stereo: StereoModel,
}

impl Default for ChemistryModel {
    fn default() -> Self {
        Self {
            valence: ValenceModel::AtomTyping {
                registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
            },
            aromaticity: AromaticityModel::daylight(),
            stereo: StereoModel::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValenceModel {
    /// Atom-typing valence model: the registry of `AtomAst` patterns.
    AtomTyping {
        registry: Cow<'static, AtomTypeRegistry>,
    },
    /// Counts valence model: the per-element covalence table.
    Counts { table: Cow<'static, ValenceTable> },
}

#[derive(Debug, Clone, PartialEq)]
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

/// Stereo perception model. `kind_models` is a per-`StereoKind` array (indexed
/// by the kind's discriminant); a `None` entry means that kind is not perceived.
/// `para_stereo` enables the graph-symmetry fixpoint iteration that resolves
/// para-stereocenters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoModel {
    pub kind_models: [Option<StereoKindModel>; StereoKind::COUNT],
    pub para_stereo: bool,
}

/// Per-kind perception settings: the elements eligible to bear this kind and
/// whether the kind's sites are treated as fluxional.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoKindModel {
    pub scope: ElementScope,
    pub fluxionality: bool,
}

impl StereoModel {
    /// The per-kind model for `kind`, or `None` if the kind is not perceived.
    pub fn kind_model(&self, kind: StereoKind) -> Option<&StereoKindModel> {
        self.kind_models[kind as usize].as_ref()
    }
}

impl Default for StereoModel {
    /// Perceive the two realized binary kinds — tetrahedral atoms and cis/trans
    /// bonds — for any element; the higher geometries (square-planar,
    /// trigonal-bipyramidal, octahedral, axial) are staged off. No para-stereo
    /// fixpoint by default.
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
    use std::borrow::Cow;
    use std::{array, ptr};

    use rstest::rstest;
    use umol_chem::element::Element;

    use super::*;
    use crate::{registry, valence_table};

    #[rstest]
    fn test_chemistry_model_default() {
        let model = ChemistryModel::default();
        assert_eq!(
            model,
            ChemistryModel {
                valence: ValenceModel::AtomTyping {
                    registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
                },
                aromaticity: AromaticityModel::daylight(),
                stereo: StereoModel::default(),
            },
        );
        match model.valence {
            ValenceModel::AtomTyping {
                registry: Cow::Borrowed(registry),
            } => assert!(ptr::eq(registry, AtomTypeRegistry::default_registry())),
            other => panic!("expected borrowed default atom-typing registry, got {other:?}"),
        }
    }

    #[rstest]
    #[case::valence(ChemistryModel {
        valence: ValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4]]),
        },
        ..ChemistryModel::default()
    })]
    #[case::aromaticity(ChemistryModel {
        aromaticity: AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.5,
        },
        ..ChemistryModel::default()
    })]
    #[case::stereo(ChemistryModel {
        stereo: StereoModel {
            para_stereo: true,
            ..StereoModel::default()
        },
        ..ChemistryModel::default()
    })]
    fn test_chemistry_model_eq_difference(#[case] other: ChemistryModel) {
        assert_ne!(ChemistryModel::default(), other);
    }

    #[rstest]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4"]),
        },
        ValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4"]),
        },
    )]
    #[case::counts(
        ValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4], O => [2]]),
        },
        ValenceModel::Counts {
            table: Cow::Owned(valence_table![O => [2], C => [4]]),
        },
    )]
    fn test_valence_model_eq(#[case] left: ValenceModel, #[case] right: ValenceModel) {
        assert_eq!(left, right);
    }

    #[rstest]
    #[case::variant(
        ValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4"]),
        },
        ValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4]]),
        },
    )]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4"]),
        },
        ValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v3"]),
        },
    )]
    #[case::counts(
        ValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4]]),
        },
        ValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [3]]),
        },
    )]
    fn test_valence_model_eq_difference(#[case] left: ValenceModel, #[case] right: ValenceModel) {
        assert_ne!(left, right);
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits::default(),
        },
    )]
    #[case::hmo(
        AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.5,
        },
        AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.5,
        },
    )]
    #[case::clar(
        AromaticityModel::Clar {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::Clar {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        },
    )]
    fn test_aromaticity_model_eq(#[case] left: AromaticityModel, #[case] right: AromaticityModel) {
        assert_eq!(left, right);
    }

    #[rstest]
    #[case::variant(
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
    )]
    #[case::hueckel_scope(
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        },
    )]
    #[case::hueckel_ring_limits(
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits {
                min_ring_size: 4,
                ..RingLimits::default()
            },
        },
    )]
    #[case::hmo_scope(
        AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.5,
        },
        AromaticityModel::Hmo {
            scope: ElementScope::AllowList(vec![Element::C]),
            stabilization_threshold: 0.5,
        },
    )]
    #[case::hmo_threshold(
        AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.5,
        },
        AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.6,
        },
    )]
    #[case::clar_scope(
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::Clar {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        },
    )]
    #[case::clar_ring_limits(
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits {
                include_fused: false,
                ..RingLimits::default()
            },
        },
    )]
    fn test_aromaticity_model_eq_difference(
        #[case] left: AromaticityModel,
        #[case] right: AromaticityModel,
    ) {
        assert_ne!(left, right);
    }

    #[rstest]
    fn test_aromaticity_model_daylight() {
        assert_eq!(
            AromaticityModel::daylight(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![
                    Element::C,
                    Element::N,
                    Element::O,
                    Element::S,
                    Element::Se,
                    Element::As,
                ]),
                ring_limits: RingLimits::default(),
            },
        );
    }

    #[rstest]
    fn test_aromaticity_model_mdl() {
        assert_eq!(
            AromaticityModel::mdl(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C, Element::N]),
                ring_limits: RingLimits {
                    min_ring_size: 6,
                    ..RingLimits::default()
                },
            },
        );
    }

    #[rstest]
    fn test_aromaticity_model_permissive() {
        assert_eq!(
            AromaticityModel::permissive(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
        );
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

    #[rstest]
    fn test_ring_limits_default() {
        assert_eq!(
            RingLimits::default(),
            RingLimits {
                min_ring_size: 3,
                max_ring_size: 22,
                include_fused: true,
                max_fused_combination: 6,
                max_fused_search: 10_000,
            },
        );
    }

    #[rstest]
    #[case::kind_models(StereoModel {
        kind_models: array::from_fn(|_| None),
        ..StereoModel::default()
    })]
    #[case::para_stereo(StereoModel {
        para_stereo: true,
        ..StereoModel::default()
    })]
    fn test_stereo_model_eq_difference(#[case] other: StereoModel) {
        assert_ne!(StereoModel::default(), other);
    }

    #[rstest]
    fn test_stereo_model_default() {
        let mut kind_models = array::from_fn(|_| None);
        kind_models[StereoKind::Tetrahedral as usize] = Some(StereoKindModel {
            scope: ElementScope::Any,
            fluxionality: false,
        });
        kind_models[StereoKind::CisTrans as usize] = Some(StereoKindModel {
            scope: ElementScope::Any,
            fluxionality: false,
        });

        assert_eq!(
            StereoModel::default(),
            StereoModel {
                kind_models,
                para_stereo: false,
            },
        );
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
}
