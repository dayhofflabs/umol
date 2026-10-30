//! Configuration for stereo transformations.

use std::collections::BTreeSet;

use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{RingConfig, StereoKind, Stereogenicity};

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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AtomId, BondId};

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
}
