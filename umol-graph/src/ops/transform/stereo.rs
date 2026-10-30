//! Configuration for stereo transformations.

use std::collections::BTreeSet;

use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    RingConfig, StereoConfigurationForm, StereoCoset, StereoKind, Stereogenicity,
};

use crate::ops::stereo::StereoSelection;

/// Site selection and algorithms for stereo addition.
///
/// Default selects tetrahedral and cis/trans kinds, stereogenic sites only,
/// and all topology sites. It uses Nauty and RingConfig::default().
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule. The executing operation checks ids before indexing or mutation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddStereoConfig {
    pub selection: StereoSelection,
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub ring_config: RingConfig,
}

impl Default for AddStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection {
                sites: None,
                kinds: Some(BTreeSet::from([
                    StereoKind::Tetrahedral,
                    StereoKind::CisTrans,
                ])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            ring_config: RingConfig::default(),
        }
    }
}

/// Which configuration states stereo removal accepts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RemoveStereoMode {
    /// Accept every configuration, including an undetermined kind.
    All,
    /// Accept only a known kind with an explicitly undetermined coset.
    Undetermined,
}

impl RemoveStereoMode {
    /// Whether the configuration satisfies this removal mode.
    ///
    /// Undetermined tests the coset state, independently of the known kind.
    /// An unknown kind has no coset to test; finite sets and symbolic cosets
    /// are not undetermined cosets. No classification or resolution is performed.
    pub fn matches(self, configuration: &StereoConfigurationForm) -> bool {
        match self {
            Self::All => true,
            Self::Undetermined => matches!(configuration.coset(), Some(StereoCoset::Undetermined)),
        }
    }
}

/// Site selection, removal mode, and algorithm for stereo removal.
///
/// Default accepts all configurations and leaves sites, kinds, and
/// stereogenicities unrestricted. It uses Nauty. One mode applies uniformly
/// to every selected site and kind.
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule or calculating symmetry. The executing operation checks ids
/// before indexing or mutation and classifies sites only when selection needs it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemoveStereoConfig {
    pub selection: StereoSelection,
    pub mode: RemoveStereoMode,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for RemoveStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection::default(),
            mode: RemoveStereoMode::All,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

/// Site selection and algorithm for stereo enumeration.
///
/// Default selects tetrahedral and cis/trans kinds, all topology sites, and
/// all stereogenicity values. It uses Nauty. Processing order is supplied
/// separately from this configuration.
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule. The executing operation checks ids before indexing or mutation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EnumerateStereoConfig {
    pub selection: StereoSelection,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for EnumerateStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection {
                sites: None,
                kinds: Some(BTreeSet::from([
                    StereoKind::Tetrahedral,
                    StereoKind::CisTrans,
                ])),
                stereogenicities: None,
            },
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AtomId, BondId, StereoTerm};

    use super::*;
    use crate::ops::stereo::StereoSite;

    #[rstest]
    fn test_add_stereo_config_default() {
        assert_eq!(
            AddStereoConfig::default(),
            AddStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: Some(BTreeSet::from([
                        StereoKind::Tetrahedral,
                        StereoKind::CisTrans,
                    ])),
                    stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                },
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                ring_config: RingConfig::default(),
            },
        );
    }

    #[rstest]
    #[case::default(
        AddStereoConfig::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::atom_site(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::bond_site(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::tetrahedral(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::cis_trans(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::intersection(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([
                    StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
                    StereoSite::Bond(BondId(0)),
                ])),
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::prochiral(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    #[case::symmetric(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
    )]
    #[case::all_stereogenicities(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: None,
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
            StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::empty_kinds(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::new()),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    #[case::unbound_sites(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([
                    StereoSite::Atom(AtomId(u32::MAX)), StereoSite::Bond(BondId(u32::MAX)),
                ])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    fn test_add_stereo_config_selection(
        #[case] config: AddStereoConfig,
        #[case] expected: Vec<StereoSite>,
    ) {
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoKind::Tetrahedral,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoKind::Tetrahedral,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoKind::Tetrahedral,
                Stereogenicity::Symmetric,
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoKind::SquarePlanar,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoKind::CisTrans,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoKind::CisTrans,
                Stereogenicity::Symmetric,
            ),
        ];
        let selected: Vec<_> = candidates
            .into_iter()
            .filter(|&(site, kind, stereogenicity)| {
                config
                    .selection
                    .matches(site, Some(kind), Some(stereogenicity))
            })
            .map(|(site, _, _)| site)
            .collect();
        assert_eq!(selected, expected);
    }

    #[rstest]
    #[case::undetermined(StereoCoset::Undetermined, true)]
    #[case::literal_zero(StereoCoset::Lit(0), false)]
    #[case::literal_one(StereoCoset::Lit(1), false)]
    #[case::finite_set(StereoCoset::lit_set([0, 1]), false)]
    #[case::singleton_set(StereoCoset::lit_set([1]), false)]
    #[case::empty_set(StereoCoset::lit_set([]), false)]
    #[case::symbolic(StereoCoset::term(StereoTerm::var("configuration")), false)]
    fn test_remove_stereo_mode_matches(
        #[case] coset: StereoCoset,
        #[case] undetermined: bool,
        #[values(
            StereoKind::Tetrahedral,
            StereoKind::CisTrans,
            StereoKind::SquarePlanar
        )]
        kind: StereoKind,
        #[values(RemoveStereoMode::All, RemoveStereoMode::Undetermined)] mode: RemoveStereoMode,
    ) {
        let configuration = StereoConfigurationForm::Kinded(kind, coset);
        let expected = match mode {
            RemoveStereoMode::All => true,
            RemoveStereoMode::Undetermined => undetermined,
        };
        assert_eq!(mode.matches(&configuration), expected);
    }

    #[rstest]
    #[case::all(RemoveStereoMode::All, true)]
    #[case::undetermined(RemoveStereoMode::Undetermined, false)]
    fn test_remove_stereo_mode_matches_kindless(
        #[case] mode: RemoveStereoMode,
        #[case] expected: bool,
    ) {
        assert_eq!(
            mode.matches(&StereoConfigurationForm::Undetermined),
            expected
        );
    }

    #[rstest]
    fn test_remove_stereo_config_default() {
        assert_eq!(
            RemoveStereoConfig::default(),
            RemoveStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: None,
                    stereogenicities: None,
                },
                mode: RemoveStereoMode::All,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            },
        );
    }

    #[rstest]
    #[case::unrestricted(
        StereoSelection::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
            StereoSite::Atom(AtomId(2)), StereoSite::Atom(AtomId(3)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1)),
            StereoSite::Atom(AtomId(4)), StereoSite::Bond(BondId(2))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(3)), StereoSite::Bond(BondId(0))],
    )]
    #[case::tetrahedral(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2))],
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::cis_trans(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::both_kinds(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::site_kind_intersection(
        StereoSelection {
            sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0)),
                StereoSite::Bond(BondId(1))])),
            kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
            stereogenicities: None,
        },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::stereogenic(
        StereoSelection {
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::non_stereogenic(
        StereoSelection {
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral])),
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
        vec![],
    )]
    #[case::empty_kinds(
        StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() },
        vec![], vec![],
    )]
    #[case::empty_sites(
        StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() },
        vec![], vec![],
    )]
    fn test_remove_stereo_config_selection(
        #[case] selection: StereoSelection,
        #[case] all: Vec<StereoSite>,
        #[case] undetermined: Vec<StereoSite>,
        #[values(RemoveStereoMode::All, RemoveStereoMode::Undetermined)] mode: RemoveStereoMode,
    ) {
        let config = RemoveStereoConfig {
            selection,
            mode,
            ..Default::default()
        };
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                Some(Stereogenicity::Stereogenic),
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                Some(Stereogenicity::Prochiral),
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoConfigurationForm::Kinded(
                    StereoKind::Tetrahedral,
                    StereoCoset::lit_set([0, 1]),
                ),
                Some(Stereogenicity::Symmetric),
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoConfigurationForm::Kinded(
                    StereoKind::SquarePlanar,
                    StereoCoset::Undetermined,
                ),
                None,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined),
                Some(Stereogenicity::Stereogenic),
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                Some(Stereogenicity::Symmetric),
            ),
            (
                StereoSite::Atom(AtomId(4)),
                StereoConfigurationForm::Undetermined,
                None,
            ),
            (
                StereoSite::Bond(BondId(2)),
                StereoConfigurationForm::Undetermined,
                None,
            ),
        ];
        let selected: Vec<_> = candidates
            .iter()
            .filter(|(site, configuration, stereogenicity)| {
                config.mode.matches(configuration)
                    && config
                        .selection
                        .matches(*site, configuration.kind(), *stereogenicity)
            })
            .map(|(site, _, _)| *site)
            .collect();
        let expected = match mode {
            RemoveStereoMode::All => all,
            RemoveStereoMode::Undetermined => undetermined,
        };
        assert_eq!(selected, expected);
    }

    #[rstest]
    fn test_enumerate_stereo_config_default() {
        assert_eq!(
            EnumerateStereoConfig::default(),
            EnumerateStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: Some(BTreeSet::from([
                        StereoKind::Tetrahedral,
                        StereoKind::CisTrans,
                    ])),
                    stereogenicities: None,
                },
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            },
        );
    }

    #[rstest]
    #[case::default(
        EnumerateStereoConfig::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::tetrahedral(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2))],
    )]
    #[case::cis_trans(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::sites(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(1))])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(1))],
    )]
    #[case::intersection(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0)),
                    StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))])),
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::prochiral(
        EnumerateStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    #[case::symmetric(
        EnumerateStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
    )]
    #[case::empty_kinds(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::new()),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    #[case::unbound_sites(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(u32::MAX)),
                    StereoSite::Bond(BondId(u32::MAX))])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    fn test_enumerate_stereo_config_selection(
        #[case] config: EnumerateStereoConfig,
        #[case] expected: Vec<StereoSite>,
    ) {
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoKind::Tetrahedral,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoKind::Tetrahedral,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoKind::Tetrahedral,
                Stereogenicity::Symmetric,
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoKind::SquarePlanar,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoKind::CisTrans,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoKind::CisTrans,
                Stereogenicity::Symmetric,
            ),
        ];
        let selected: Vec<_> = candidates
            .into_iter()
            .filter(|&(site, kind, stereogenicity)| {
                config
                    .selection
                    .matches(site, Some(kind), Some(stereogenicity))
            })
            .map(|(site, _, _)| site)
            .collect();
        assert_eq!(selected, expected);
    }
}
