//! Shared proptest generators for the umol-graph property suite.
//!
//! `select_scenario` emits the aromatic-selection operational domain: one- and
//! two-ring skeletons (fused or coupled) under the Hückel rule, every ring
//! atom carrying a carrier entry, at most three flexible atoms drawn from a
//! fixed completion pool, all contributions literal, no stored systems and no
//! atom-level assertions. One localized bond carries a positive aromaticity
//! constraint to activate selection. Policies vary over both failure policies,
//! both aromaticity tie-breaks, and both valence tie-breaks.
//!
//! `aromatic_skeleton` emits aromatic skeletons: one five- or six-membered
//! ring, two fused rings, or two rings coupled by a bond, over C, N, O, and S
//! with single bonds and atom charge 0. Each atom asserts no aromatic valence,
//! `#a+`, `#a0` to `#a2`, or `#a!`, and has an undetermined hydrogen count or
//! a literal one from 0 to 2; each bond asserts no aromaticity, `#a`, or
//! `#a!`. Optionally the first ring is a stored concrete closed-shell
//! aromatic system with literal contributions from 0 to 2. Each atom with
//! three ring neighbors optionally asserts `#T+`, `#T0`, `#T1`, or `#T!`.
//! The domain keeps contradictory and underdetermined skeletons. Each
//! skeleton also renders as SMILES and as a V2000 molfile, stating what the
//! format can express.

use std::collections::BTreeSet;

use proptest::collection::vec;
use proptest::prelude::*;
use proptest::sample::subsequence;
use proptest::strategy::ValueTree;
use proptest::test_runner::{Config, FileFailurePersistence, TestRunner};
use rstest::rstest;
use smallvec::{smallvec, SmallVec};
use umol_chem::element::Element;
use umol_graph::ops::model::{
    AromaticityModel, AromaticityRule, AromaticityTieBreak, ElementScope, RingLimits,
    ValenceTieBreak,
};
use umol_graph::ops::resolve::{AromaticityFailurePolicy, AromaticityResolveConfig};
use umol_graph::ops::valence::AtomCompletions;
use umol_graph_ir::atom_dsl;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AsLit, AtomConstraintForm, AtomForm, AtomId,
    BondConstraintForm, BondForm, BooleanForm, ElementForm, Lattice, Molecule, MoleculeEntries,
    NumForm, StereoCoset, TetrahedralStereoForm, TryIntoIr, UnpairedElectronsForm,
};
use umol_io::ctfile::parser::parse_mol_to_ir;
use umol_io::smiles::Smiles;

#[derive(Clone, Debug)]
pub(crate) struct SelectScenario {
    pub(crate) molecule: Molecule,
    pub(crate) completions: AtomCompletions,
    pub(crate) model: AromaticityModel,
    pub(crate) config: AromaticityResolveConfig,
    pub(crate) tie_break: ValenceTieBreak,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
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

fn ring_topology() -> impl Strategy<Value = RingTopology> {
    prop_oneof![
        (5usize..=6).prop_map(RingTopology::Single),
        ((5usize..=6), (5usize..=6)).prop_map(|(a, b)| RingTopology::Fused(a, b)),
        ((5usize..=6), (5usize..=6)).prop_map(|(a, b)| RingTopology::Coupled(a, b)),
    ]
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
    ring_topology().prop_flat_map(|topology| {
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
                        .enumerate()
                        .map(|(index, (from, to))| {
                            let attributes = if index == 0 { "1#a+" } else { "1" };
                            format!("[{from} {to} \"{attributes}\"]")
                        })
                        .collect::<Vec<_>>()
                        .join(" ");
                    let molecule: Molecule = format!("{{:atoms [{atoms}] :bonds [{bond_list}]}}")
                        .parse()
                        .expect("generated skeleton parses");
                    let completions = AtomCompletions::from_iter((0..atom_count).map(|atom| {
                        let pool = flexible_atoms
                            .iter()
                            .position(|&flexible| flexible == atom)
                            .map_or(0, |idx| pool_picks[idx]);
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

#[derive(Clone, Debug)]
struct SkeletonAtom {
    element: Element,
    implicit_hydrogens: NumForm,
    aromatic_valence: AromaticValenceForm,
    tetrahedral_stereo: TetrahedralStereoForm,
}

/// One skeleton of the `aromatic_skeleton` domain, held as the parameters its molecule and
/// renderings derive from.
#[derive(Clone, Debug)]
pub(crate) struct AromaticSkeleton {
    topology: RingTopology,
    atoms: Vec<SkeletonAtom>,
    bond_aromaticity: Vec<BooleanForm>,
    /// Electron contributions of the stored aromatic system on the first ring.
    stored_aromatic_system: Option<Vec<i64>>,
}

/// Depth-first spanning tree from atom 0, visiting neighbors in bond order: the tree
/// children of each atom with their bonds, and the ring-closure bonds at each atom.
struct SmilesTree {
    children: Vec<Vec<(usize, usize)>>,
    ring_closures: Vec<Vec<usize>>,
}

impl SmilesTree {
    fn new(atom_count: usize, bonds: &[(usize, usize)]) -> Self {
        let mut neighbors = vec![Vec::new(); atom_count];
        for (bond, &(from, to)) in bonds.iter().enumerate() {
            neighbors[from].push((to, bond));
            neighbors[to].push((from, bond));
        }
        let mut tree = Self {
            children: vec![Vec::new(); atom_count],
            ring_closures: vec![Vec::new(); atom_count],
        };
        tree.visit(
            0,
            &neighbors,
            &mut vec![false; atom_count],
            &mut vec![false; bonds.len()],
        );
        tree
    }

    fn visit(
        &mut self,
        atom: usize,
        neighbors: &[Vec<(usize, usize)>],
        visited: &mut [bool],
        traversed: &mut [bool],
    ) {
        visited[atom] = true;
        for &(next, bond) in &neighbors[atom] {
            if traversed[bond] {
                continue;
            }
            traversed[bond] = true;
            if visited[next] {
                self.ring_closures[atom].push(bond);
                self.ring_closures[next].push(bond);
            } else {
                self.children[atom].push((next, bond));
                self.visit(next, neighbors, visited, traversed);
            }
        }
    }
}

pub(crate) fn aromatic_skeleton() -> impl Strategy<Value = AromaticSkeleton> {
    ring_topology().prop_flat_map(|topology| {
        let (atom_count, bonds) = skeleton(&topology);
        let first_ring = match topology {
            RingTopology::Single(a) | RingTopology::Fused(a, _) | RingTopology::Coupled(a, _) => a,
        };
        let atoms: Vec<BoxedStrategy<SkeletonAtom>> = (0..atom_count)
            .map(|atom| {
                let degree = bonds
                    .iter()
                    .filter(|&&(from, to)| from == atom || to == atom)
                    .count();
                let tetrahedral_stereo = if degree == 3 {
                    prop_oneof![
                        Just(TetrahedralStereoForm::Undetermined),
                        Just(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined)),
                        (0u32..=1).prop_map(|coset| {
                            TetrahedralStereoForm::Stereo(StereoCoset::Lit(coset))
                        }),
                        Just(TetrahedralStereoForm::NotStereo),
                    ]
                    .boxed()
                } else {
                    Just(TetrahedralStereoForm::Undetermined).boxed()
                };
                (
                    prop::sample::select(vec![Element::C, Element::N, Element::O, Element::S]),
                    prop_oneof![
                        Just(NumForm::Undetermined),
                        (0i64..=2).prop_map(NumForm::Lit)
                    ],
                    prop_oneof![
                        Just(AromaticValenceForm::Undetermined),
                        Just(AromaticValenceForm::Aromatic(NumForm::Undetermined)),
                        (0i64..=2).prop_map(|electrons| {
                            AromaticValenceForm::Aromatic(NumForm::Lit(electrons))
                        }),
                        Just(AromaticValenceForm::NotAromatic),
                    ],
                    tetrahedral_stereo,
                )
                    .prop_map(
                        |(element, implicit_hydrogens, aromatic_valence, tetrahedral_stereo)| {
                            SkeletonAtom {
                                element,
                                implicit_hydrogens,
                                aromatic_valence,
                                tetrahedral_stereo,
                            }
                        },
                    )
                    .boxed()
            })
            .collect();
        (
            atoms,
            vec(
                prop_oneof![
                    Just(BooleanForm::Undetermined),
                    Just(BooleanForm::Lit(true)),
                    Just(BooleanForm::Lit(false)),
                ],
                bonds.len(),
            ),
            prop::option::of(vec(0i64..=2, first_ring)),
        )
            .prop_map(move |(atoms, bond_aromaticity, stored_aromatic_system)| {
                AromaticSkeleton {
                    topology: topology.clone(),
                    atoms,
                    bond_aromaticity,
                    stored_aromatic_system,
                }
            })
    })
}

impl AromaticSkeleton {
    /// The skeleton as a molecule. An undetermined parameter stores no assertion.
    pub(crate) fn molecule(&self) -> Molecule {
        let (_, bonds) = skeleton(&self.topology);
        let atom_id = |index: usize| AtomId(u32::try_from(index).expect("small skeleton"));
        Molecule::from_entries(MoleculeEntries {
            atoms: self
                .atoms
                .iter()
                .map(|atom| AtomForm {
                    element: ElementForm::Lit(atom.element),
                    charge: NumForm::Lit(0),
                    implicit_hydrogens: atom.implicit_hydrogens.clone(),
                    constraints: [
                        AtomConstraintForm::AromaticValence(atom.aromatic_valence.clone()),
                        AtomConstraintForm::TetrahedralStereo(atom.tetrahedral_stereo.clone()),
                    ]
                    .into_iter()
                    .filter(|constraint| !constraint.is_undetermined())
                    .collect(),
                    ..AtomForm::default()
                })
                .collect(),
            bonds: bonds
                .iter()
                .zip(&self.bond_aromaticity)
                .map(|(&(from, to), &aromatic)| {
                    (
                        atom_id(from),
                        atom_id(to),
                        BondForm {
                            order: NumForm::Lit(1),
                            constraints: [BondConstraintForm::Aromatic(aromatic)]
                                .into_iter()
                                .filter(|constraint| !constraint.is_undetermined())
                                .collect(),
                            ..BondForm::default()
                        },
                    )
                })
                .collect(),
            aromatic: self
                .stored_aromatic_system
                .iter()
                .map(|electrons| {
                    (
                        (0..electrons.len()).map(atom_id).collect(),
                        AromaticSystemForm::from_electrons(electrons.clone())
                            .with_charge(NumForm::Lit(0))
                            .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                    )
                })
                .collect(),
            ..MoleculeEntries::default()
        })
    }

    /// SMILES text for the skeleton. An atom asserting an aromatic valence is lowercase, an
    /// atom with a literal hydrogen count is a bracket atom with that count, a bond asserting
    /// aromaticity is `:`, and one asserting non-aromaticity is `-`. A bracket atom with at
    /// most one hydrogen writes coset 0 as `@` and coset 1 as `@@`; the parsed coset follows
    /// SMILES neighbor order. SMILES cannot state a contribution count, a non-aromatic atom,
    /// a stored aromatic system, `#T+`, or `#T!`, and the text omits them.
    pub(crate) fn smiles(&self) -> String {
        let (atom_count, bonds) = skeleton(&self.topology);
        let tree = SmilesTree::new(atom_count, &bonds);
        let mut labels = vec![None; bonds.len()];
        let mut text = String::new();
        self.write_smiles(0, &tree, &mut labels, &mut text);
        text
    }

    fn write_smiles(
        &self,
        atom: usize,
        tree: &SmilesTree,
        labels: &mut [Option<usize>],
        text: &mut String,
    ) {
        let skeleton_atom = &self.atoms[atom];
        let symbol = match skeleton_atom.aromatic_valence {
            AromaticValenceForm::Aromatic(_) => skeleton_atom.element.symbol().to_ascii_lowercase(),
            _ => skeleton_atom.element.symbol().to_owned(),
        };
        match skeleton_atom.implicit_hydrogens {
            NumForm::Lit(hydrogens) => {
                let chirality = match skeleton_atom.tetrahedral_stereo {
                    TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)) if hydrogens <= 1 => "@",
                    TetrahedralStereoForm::Stereo(StereoCoset::Lit(_)) if hydrogens <= 1 => "@@",
                    _ => "",
                };
                let hydrogen = match hydrogens {
                    0 => String::new(),
                    1 => "H".to_owned(),
                    count => format!("H{count}"),
                };
                text.push_str(&format!("[{symbol}{chirality}{hydrogen}]"));
            }
            _ => text.push_str(&symbol),
        }
        for &bond in &tree.ring_closures[atom] {
            let label = match labels[bond] {
                Some(label) => {
                    text.push_str(self.smiles_bond(bond));
                    label
                }
                None => {
                    let label = labels.iter().flatten().count() + 1;
                    labels[bond] = Some(label);
                    label
                }
            };
            text.push_str(&label.to_string());
        }
        let children = &tree.children[atom];
        for (index, &(child, bond)) in children.iter().enumerate() {
            let branch = index + 1 < children.len();
            if branch {
                text.push('(');
            }
            text.push_str(self.smiles_bond(bond));
            self.write_smiles(child, tree, labels, text);
            if branch {
                text.push(')');
            }
        }
    }

    fn smiles_bond(&self, bond: usize) -> &'static str {
        match self.bond_aromaticity[bond] {
            BooleanForm::Lit(true) => ":",
            BooleanForm::Lit(false) => "-",
            BooleanForm::Undetermined => "",
        }
    }

    /// V2000 molfile text for the skeleton, every atom at the origin. A literal hydrogen
    /// count fills the hydrogen-count field, a bond asserting aromaticity has type 4 and every
    /// other bond type 1, and an atom with a literal count of at most one hydrogen writes
    /// coset 0, coset 1, and `#T+` as parity 1, 2, and 3. The format has no atom aromaticity,
    /// contribution count, stored aromatic system, or `#T!`, and type 1 covers both a bond
    /// asserting non-aromaticity and one without an assertion.
    pub(crate) fn molfile(&self) -> String {
        let (_, bonds) = skeleton(&self.topology);
        let mut text = format!(
            "skeleton\n  umol\n\n{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
            self.atoms.len(),
            bonds.len()
        );
        for atom in &self.atoms {
            let (hydrogens, parity) = match atom.implicit_hydrogens {
                NumForm::Lit(count) => (
                    count + 1,
                    match atom.tetrahedral_stereo {
                        TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)) if count <= 1 => 1,
                        TetrahedralStereoForm::Stereo(StereoCoset::Lit(_)) if count <= 1 => 2,
                        TetrahedralStereoForm::Stereo(StereoCoset::Undetermined) if count <= 1 => 3,
                        _ => 0,
                    },
                ),
                _ => (0, 0),
            };
            text.push_str(&format!(
                "{:10.4}{:10.4}{:10.4} {:<3} 0  0{parity:3}{hydrogens:3}  0  0  0  0  0  0  0  0\n",
                0.0,
                0.0,
                0.0,
                atom.element.symbol()
            ));
        }
        for (&(from, to), &aromatic) in bonds.iter().zip(&self.bond_aromaticity) {
            let bond_type = if aromatic == BooleanForm::Lit(true) {
                4
            } else {
                1
            };
            text.push_str(&format!(
                "{:3}{:3}{bond_type:3}  0  0  0  0\n",
                from + 1,
                to + 1
            ));
        }
        text.push_str("M  END\n");
        text
    }
}

#[rstest]
fn test_aromatic_skeleton_coverage() {
    let strategy = aromatic_skeleton();
    let mut runner = TestRunner::deterministic();
    let skeletons: Vec<AromaticSkeleton> = (0..512)
        .map(|_| {
            strategy
                .new_tree(&mut runner)
                .expect("skeleton generation")
                .current()
        })
        .collect();
    let molecules: Vec<Molecule> = skeletons.iter().map(AromaticSkeleton::molecule).collect();
    let atoms: Vec<AtomForm> = molecules
        .iter()
        .flat_map(|molecule| {
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes().clone())
                .collect::<Vec<_>>()
        })
        .collect();
    assert_eq!(
        skeletons
            .iter()
            .map(|skeleton| skeleton.topology.clone())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            RingTopology::Single(5),
            RingTopology::Single(6),
            RingTopology::Fused(5, 5),
            RingTopology::Fused(5, 6),
            RingTopology::Fused(6, 5),
            RingTopology::Fused(6, 6),
            RingTopology::Coupled(5, 5),
            RingTopology::Coupled(5, 6),
            RingTopology::Coupled(6, 5),
            RingTopology::Coupled(6, 6),
        ])
    );
    assert_eq!(
        atoms
            .iter()
            .map(|atom| atom.element.clone())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([Element::C, Element::N, Element::O, Element::S].map(ElementForm::Lit))
    );
    assert_eq!(
        atoms
            .iter()
            .map(|atom| atom.implicit_hydrogens.clone())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            NumForm::Undetermined,
            NumForm::Lit(0),
            NumForm::Lit(1),
            NumForm::Lit(2),
        ])
    );
    assert_eq!(
        atoms
            .iter()
            .map(|atom| atom.constraints.aromatic_valence().cloned())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            None,
            Some(AromaticValenceForm::Aromatic(NumForm::Undetermined)),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(0))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(2))),
            Some(AromaticValenceForm::NotAromatic),
        ])
    );
    assert_eq!(
        atoms
            .iter()
            .map(|atom| atom.constraints.tetrahedral_stereo().cloned())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            None,
            Some(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined)),
            Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))),
            Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))),
            Some(TetrahedralStereoForm::NotStereo),
        ])
    );
    assert_eq!(
        molecules
            .iter()
            .flat_map(|molecule| {
                molecule
                    .bonds()
                    .iter()
                    .map(|bond| bond.attributes().constraints.aromatic())
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            BooleanForm::Undetermined,
            BooleanForm::Lit(true),
            BooleanForm::Lit(false),
        ])
    );
    assert_eq!(
        molecules
            .iter()
            .map(|molecule| molecule.aromatic_systems().count())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([0, 1])
    );
    assert_eq!(
        molecules
            .iter()
            .flat_map(|molecule| {
                molecule
                    .aromatic_systems()
                    .iter()
                    .flat_map(|aromatic_system| {
                        aromatic_system
                            .attributes()
                            .electrons
                            .as_lit()
                            .expect("stored contributions are literal")
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([0, 1, 2])
    );
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_aromatic_skeleton_smiles(skeleton in aromatic_skeleton()) {
        let raised: Molecule = Smiles::parse(&skeleton.smiles())
            .expect("SMILES rendering parses")
            .as_table_ir()
            .try_into_ir(&())
            .expect("SMILES rendering raises");
        let mut expected: Vec<(ElementForm, NumForm, bool)> = skeleton
            .atoms
            .iter()
            .map(|atom| (
                ElementForm::Lit(atom.element),
                atom.implicit_hydrogens.clone(),
                matches!(atom.aromatic_valence, AromaticValenceForm::Aromatic(_)),
            ))
            .collect();
        expected.sort();
        let mut actual: Vec<(ElementForm, NumForm, bool)> = raised
            .atoms()
            .iter()
            .map(|atom| {
                let attributes = atom.attributes();
                (
                    attributes.element.clone(),
                    attributes.implicit_hydrogens.clone(),
                    matches!(
                        attributes.constraints.aromatic_valence(),
                        Some(AromaticValenceForm::Aromatic(_))
                    ),
                )
            })
            .collect();
        actual.sort();
        prop_assert_eq!(actual, expected);
        prop_assert_eq!(raised.bonds().count(), skeleton.bond_aromaticity.len());
    }

    #[test]
    fn test_aromatic_skeleton_molfile(skeleton in aromatic_skeleton()) {
        let raised = parse_mol_to_ir(&skeleton.molfile()).expect("molfile rendering raises");
        prop_assert_eq!(
            raised
                .atoms()
                .iter()
                .map(|atom| (
                    atom.attributes().element.clone(),
                    atom.attributes().implicit_hydrogens.clone(),
                ))
                .collect::<Vec<_>>(),
            skeleton
                .atoms
                .iter()
                .map(|atom| (ElementForm::Lit(atom.element), atom.implicit_hydrogens.clone()))
                .collect::<Vec<_>>()
        );
        prop_assert_eq!(
            raised
                .bonds()
                .iter()
                .map(|bond| bond.attributes().constraints.aromatic())
                .collect::<Vec<_>>(),
            skeleton
                .bond_aromaticity
                .iter()
                .map(|&aromatic| match aromatic {
                    BooleanForm::Lit(true) => aromatic,
                    _ => BooleanForm::Undetermined,
                })
                .collect::<Vec<_>>()
        );
    }
}
