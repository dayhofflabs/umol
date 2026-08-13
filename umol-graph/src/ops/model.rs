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
use umol_graph_ir::ir::{AtomFieldKind, StereoKind};

use crate::ops::valence::{AtomTypeRegistry, ValenceTable};
use crate::ops::validate::ConnectivityModel;
use crate::utils::SortingDirection;

#[derive(Debug, Clone, PartialEq)]
pub struct ChemistryModel {
    pub connectivity: ConnectivityModel,
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
    pub stereo: StereoModel,
}

impl Default for ChemistryModel {
    fn default() -> Self {
        Self {
            connectivity: ConnectivityModel::default(),
            valence: ValenceModel::default(),
            aromaticity: AromaticityModel::daylight(),
            stereo: StereoModel::default(),
        }
    }
}

/// Valence model: where candidate states come from, and how plural survivors
/// are disposed of.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceModel {
    pub candidates: ValenceCandidateSource,
    pub tie_break: ValenceTieBreak,
}

impl ValenceModel {
    /// The atom-typing source with the `Strict` tie-break.
    pub fn atom_typing(registry: Cow<'static, AtomTypeRegistry>) -> Self {
        Self {
            candidates: ValenceCandidateSource::AtomTyping { registry },
            tie_break: ValenceTieBreak::Strict,
        }
    }

    /// The counts source with the `Strict` tie-break.
    pub fn counts(table: Cow<'static, ValenceTable>) -> Self {
        Self {
            candidates: ValenceCandidateSource::Counts { table },
            tie_break: ValenceTieBreak::Strict,
        }
    }

    /// The umol SMILES reading: the owned SMILES valence table with the
    /// `MostSaturated` tie-break.
    pub fn smiles() -> Self {
        Self {
            candidates: ValenceCandidateSource::Counts {
                table: Cow::Borrowed(ValenceTable::smiles_table()),
            },
            tie_break: ValenceTieBreak::MostSaturated,
        }
    }

    /// The MDL/CTfile reading: the frozen MDL valence table with the
    /// `MostSaturated` tie-break.
    pub fn mdl() -> Self {
        Self {
            candidates: ValenceCandidateSource::Counts {
                table: Cow::Borrowed(ValenceTable::mdl_table()),
            },
            tie_break: ValenceTieBreak::MostSaturated,
        }
    }
}

impl Default for ValenceModel {
    /// The default atom-typing registry with the `Strict` tie-break.
    fn default() -> Self {
        Self::atom_typing(Cow::Borrowed(AtomTypeRegistry::default_registry()))
    }
}

/// Where valence candidate states come from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValenceCandidateSource {
    /// The registry of `AtomForm` patterns.
    AtomTyping {
        registry: Cow<'static, AtomTypeRegistry>,
    },
    /// The per-element covalence table.
    Counts { table: Cow<'static, ValenceTable> },
}

/// Disposal policy for plural candidate survivors: named lexicographic keys
/// only, no open key construction.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ValenceTieBreak {
    /// No selection: plural survivors stay in the report.
    #[default]
    Strict,
    /// Max implicit hydrogens, then max lone pairs, then min unpaired
    /// electrons — the closed-shell most-saturated reading.
    MostSaturated,
}

impl ValenceTieBreak {
    /// The policy's lexicographic key: candidates are ordered by each pair in
    /// sequence and the greatest is selected. An empty key selects nothing;
    /// a tie surviving the full key stays plural.
    pub fn key(&self) -> &'static [(AtomFieldKind, SortingDirection)] {
        match self {
            Self::Strict => &[],
            Self::MostSaturated => &[
                (
                    AtomFieldKind::ImplicitHydrogens,
                    SortingDirection::Ascending,
                ),
                (AtomFieldKind::LonePairs, SortingDirection::Ascending),
                (
                    AtomFieldKind::UnpairedElectrons,
                    SortingDirection::Descending,
                ),
            ],
        }
    }
}

/// Aromaticity model: the participating elements and the perception rule.
#[derive(Debug, Clone, PartialEq)]
pub struct AromaticityModel {
    pub scope: ElementScope,
    pub rule: AromaticityRule,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AromaticityRule {
    Hueckel { ring_limits: RingLimits },
    Hmo { stabilization_threshold: f64 },
    Clar { ring_limits: RingLimits },
}

impl AromaticityModel {
    /// Daylight (SMILES) aromaticity scope: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self {
            scope: ElementScope::AllowList(vec![
                Element::C,
                Element::N,
                Element::O,
                Element::S,
                Element::Se,
                Element::As,
            ]),
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
        }
    }

    /// MDL (MOL/SDF) aromaticity scope: C and N only, minimum ring size 6.
    pub fn mdl() -> Self {
        Self {
            scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits {
                    min_ring_size: 6,
                    ..RingLimits::default()
                },
            },
        }
    }

    /// Permissive aromaticity scope: any element.
    pub fn permissive() -> Self {
        Self {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
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
                connectivity: ConnectivityModel::default(),
                valence: ValenceModel::atom_typing(Cow::Borrowed(
                    AtomTypeRegistry::default_registry()
                )),
                aromaticity: AromaticityModel::daylight(),
                stereo: StereoModel::default(),
            },
        );
        assert_eq!(model.valence.tie_break, ValenceTieBreak::Strict);
        match model.valence.candidates {
            ValenceCandidateSource::AtomTyping {
                registry: Cow::Borrowed(registry),
            } => assert!(ptr::eq(registry, AtomTypeRegistry::default_registry())),
            other => panic!("expected borrowed default atom-typing registry, got {other:?}"),
        }
    }

    #[rstest]
    #[case::atom_typing(
        ValenceModel::atom_typing(Cow::Borrowed(AtomTypeRegistry::default_registry())),
        ValenceCandidateSource::AtomTyping {
            registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
        },
        ValenceTieBreak::Strict,
    )]
    #[case::counts(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        ValenceCandidateSource::Counts {
            table: Cow::Borrowed(ValenceTable::default_table()),
        },
        ValenceTieBreak::Strict,
    )]
    #[case::smiles(
        ValenceModel::smiles(),
        ValenceCandidateSource::Counts {
            table: Cow::Borrowed(ValenceTable::smiles_table()),
        },
        ValenceTieBreak::MostSaturated,
    )]
    #[case::mdl(
        ValenceModel::mdl(),
        ValenceCandidateSource::Counts {
            table: Cow::Borrowed(ValenceTable::mdl_table()),
        },
        ValenceTieBreak::MostSaturated,
    )]
    fn test_valence_model_constructors(
        #[case] model: ValenceModel,
        #[case] candidates: ValenceCandidateSource,
        #[case] tie_break: ValenceTieBreak,
    ) {
        assert_eq!(
            model,
            ValenceModel {
                candidates,
                tie_break,
            },
        );
    }

    #[rstest]
    #[case::strict(ValenceTieBreak::Strict, &[])]
    #[case::most_saturated(
        ValenceTieBreak::MostSaturated,
        &[
            (AtomFieldKind::ImplicitHydrogens, SortingDirection::Ascending),
            (AtomFieldKind::LonePairs, SortingDirection::Ascending),
            (AtomFieldKind::UnpairedElectrons, SortingDirection::Descending),
        ],
    )]
    fn test_valence_tie_break_key(
        #[case] tie_break: ValenceTieBreak,
        #[case] expected: &[(AtomFieldKind, SortingDirection)],
    ) {
        assert_eq!(tie_break.key(), expected);
    }

    #[rstest]
    #[case::valence(ChemistryModel {
        valence: ValenceModel::counts(Cow::Owned(valence_table![C => [4]])),
        ..ChemistryModel::default()
    })]
    #[case::aromaticity(ChemistryModel {
        aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
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
        ValenceModel::atom_typing(Cow::Owned(registry!["C#c0#v4"])),
        ValenceModel::atom_typing(Cow::Owned(registry!["C#c0#v4"])),
    )]
    #[case::counts(
        ValenceModel::counts(Cow::Owned(valence_table![C => [4], O => [2]])),
        ValenceModel::counts(Cow::Owned(valence_table![O => [2], C => [4]])),
    )]
    fn test_valence_model_eq(#[case] left: ValenceModel, #[case] right: ValenceModel) {
        assert_eq!(left, right);
    }

    #[rstest]
    #[case::variant(
        ValenceModel::atom_typing(Cow::Owned(registry!["C#c0#v4"])),
        ValenceModel::counts(Cow::Owned(valence_table![C => [4]])),
    )]
    #[case::atom_typing(
        ValenceModel::atom_typing(Cow::Owned(registry!["C#c0#v4"])),
        ValenceModel::atom_typing(Cow::Owned(registry!["C#c0#v3"])),
    )]
    #[case::counts(
        ValenceModel::counts(Cow::Owned(valence_table![C => [4]])),
        ValenceModel::counts(Cow::Owned(valence_table![C => [3]])),
    )]
    #[case::tie_break(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        ValenceModel {
            tie_break: ValenceTieBreak::MostSaturated,
            ..ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()))
        },
    )]
    fn test_valence_model_eq_difference(#[case] left: ValenceModel, #[case] right: ValenceModel) {
        assert_ne!(left, right);
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C, Element::N]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C, Element::N]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
    )]
    #[case::hmo(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
    )]
    #[case::clar(
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
    )]
    fn test_aromaticity_model_eq(#[case] left: AromaticityModel, #[case] right: AromaticityModel) {
        assert_eq!(left, right);
    }

    #[rstest]
    #[case::variant(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
    )]
    #[case::hueckel_scope(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
    )]
    #[case::hueckel_ring_limits(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hueckel { ring_limits: RingLimits {
                min_ring_size: 4,
                ..RingLimits::default()
            } } },
    )]
    #[case::hmo_scope(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
    )]
    #[case::hmo_threshold(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.6 } },
    )]
    #[case::clar_scope(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
    )]
    #[case::clar_ring_limits(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits {
                include_fused: false,
                ..RingLimits::default()
            } } },
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
            AromaticityModel {
                scope: ElementScope::AllowList(vec![
                    Element::C,
                    Element::N,
                    Element::O,
                    Element::S,
                    Element::Se,
                    Element::As,
                ]),
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::default()
                }
            },
        );
    }

    #[rstest]
    fn test_aromaticity_model_mdl() {
        assert_eq!(
            AromaticityModel::mdl(),
            AromaticityModel {
                scope: ElementScope::AllowList(vec![Element::C, Element::N]),
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits {
                        min_ring_size: 6,
                        ..RingLimits::default()
                    }
                }
            },
        );
    }

    #[rstest]
    fn test_aromaticity_model_permissive() {
        assert_eq!(
            AromaticityModel::permissive(),
            AromaticityModel {
                scope: ElementScope::Any,
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::default()
                }
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
