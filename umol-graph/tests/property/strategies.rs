//! Shared proptest generators for the umol-graph property suite.
//!
//! `select_scenario` emits the aromatic-selection operational domain: one- and
//! two-ring skeletons (fused or coupled) under the Hückel rule, every ring
//! atom carrying a carrier entry, at most three flexible atoms drawn from a
//! fixed completion pool, all contributions literal, no stored systems and no
//! atom-level assertions. Policies vary over both failure policies, both
//! aromaticity tie-breaks, and both valence tie-breaks.

use proptest::collection::vec;
use proptest::prelude::*;
use proptest::sample::subsequence;
use smallvec::{smallvec, SmallVec};
use umol_graph::ops::model::{
    AromaticityModel, AromaticityRule, AromaticityTieBreak, ElementScope, RingLimits,
    ValenceTieBreak,
};
use umol_graph::ops::resolve::{AromaticityFailurePolicy, AromaticityResolveConfig};
use umol_graph::ops::valence::AtomCompletions;
use umol_graph_ir::atom_dsl;
use umol_graph_ir::ir::{AtomForm, AtomId, Molecule};

#[derive(Clone, Debug)]
pub(crate) struct SelectScenario {
    pub(crate) molecule: Molecule,
    pub(crate) completions: AtomCompletions,
    pub(crate) model: AromaticityModel,
    pub(crate) config: AromaticityResolveConfig,
    pub(crate) tie_break: ValenceTieBreak,
}

#[derive(Clone, Debug)]
enum RingTopology {
    Single(usize),
    Fused(usize, usize),
    Coupled(usize, usize),
}

/// Atom count and bond list of the skeleton: ring A on `0..a`; a fused ring B
/// shares the edge `(a - 2, a - 1)`; a coupled ring B is bridged over
/// `(0, a)`.
fn skeleton(topology: &RingTopology) -> (usize, Vec<(usize, usize)>) {
    match *topology {
        RingTopology::Single(a) => {
            let bonds = (0..a).map(|i| (i, (i + 1) % a)).collect();
            (a, bonds)
        }
        RingTopology::Fused(a, b) => {
            let mut bonds: Vec<(usize, usize)> = (0..a).map(|i| (i, (i + 1) % a)).collect();
            bonds.push((a - 2, a));
            for i in 0..(b - 3) {
                bonds.push((a + i, a + i + 1));
            }
            bonds.push((a + b - 3, a - 1));
            (a + b - 2, bonds)
        }
        RingTopology::Coupled(a, b) => {
            let mut bonds: Vec<(usize, usize)> = (0..a).map(|i| (i, (i + 1) % a)).collect();
            bonds.extend((0..b).map(|i| (a + i, a + (i + 1) % b)));
            bonds.push((0, a));
            (a + b, bonds)
        }
    }
}

/// The completion pool: index 0 is the fixed unit contributor; the others are
/// flexible — aromatic-only pairs, a mixed aromatic/non-aromatic pair, and an
/// aromatic-only triple.
fn completion_pool(index: usize) -> SmallVec<[AtomForm; 1]> {
    match index {
        0 => smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")],
        1 => smallvec![
            atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
        ],
        2 => smallvec![
            atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a"),
            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
        ],
        3 => smallvec![
            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            atom_dsl!("C#i=#c0#h2#n0#u0#s#v2#a!"),
        ],
        _ => smallvec![
            atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
            atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a"),
            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
        ],
    }
}

pub(crate) fn select_scenario() -> impl Strategy<Value = SelectScenario> {
    let topology = prop_oneof![
        (5usize..=6).prop_map(RingTopology::Single),
        ((5usize..=6), (5usize..=6)).prop_map(|(a, b)| RingTopology::Fused(a, b)),
        ((5usize..=6), (5usize..=6)).prop_map(|(a, b)| RingTopology::Coupled(a, b)),
    ];
    topology.prop_flat_map(|topology| {
        let (atom_count, bonds) = skeleton(&topology);
        (
            subsequence((0..atom_count).collect::<Vec<usize>>(), 0..=3),
            vec(1usize..=4, 3),
            prop_oneof![
                Just(AromaticityFailurePolicy::Error),
                Just(AromaticityFailurePolicy::Keep),
            ],
            prop_oneof![
                Just(AromaticityTieBreak::Strict),
                Just(AromaticityTieBreak::MinElectronCount),
            ],
            prop_oneof![
                Just(ValenceTieBreak::Strict),
                Just(ValenceTieBreak::MostSaturated)
            ],
        )
            .prop_map(
                move |(flexible_atoms, pool_picks, failure, structural, value)| {
                    let atoms = vec!["\"C#c0\""; atom_count].join(" ");
                    let bond_list = bonds
                        .iter()
                        .map(|(from, to)| format!("[{from} {to} \"1\"]"))
                        .collect::<Vec<_>>()
                        .join(" ");
                    let molecule: Molecule = format!("{{:atoms [{atoms}] :bonds [{bond_list}]}}")
                        .parse()
                        .expect("generated skeleton parses");
                    let completions = AtomCompletions::from_iter((0..atom_count).map(|atom| {
                        let pool = flexible_atoms
                            .iter()
                            .position(|&flexible| flexible == atom)
                            .map_or(0, |slot| pool_picks[slot]);
                        (
                            AtomId(u32::try_from(atom).expect("small skeleton")),
                            completion_pool(pool),
                        )
                    }));
                    SelectScenario {
                        molecule,
                        completions,
                        model: AromaticityModel {
                            scope: ElementScope::Any,
                            rule: AromaticityRule::Hueckel {
                                ring_limits: RingLimits::default(),
                            },
                            tie_break: structural,
                        },
                        config: AromaticityResolveConfig {
                            aromatic_valence_failure: failure,
                            ..AromaticityResolveConfig::default()
                        },
                        tie_break: value,
                    }
                },
            )
    })
}
