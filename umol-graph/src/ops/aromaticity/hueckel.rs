//! Hückel-rule (4n+2 electron counting) aromaticity perception.
//!
//! Filters candidate atoms by element scope and the per-atom aromatic-valence
//! constraint, enumerates rings within configured bounds, checks the Hueckel
//! 4n+2 rule on individual and fused ring combinations, and produces aromatic
//! system tuples `(Vec<AtomId>, AromaticSystemForm)` ready for `Molecule::edit`.

use std::collections::HashSet;

use umol_graph_core::UnionFind;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomId, ElementForm, Molecule, RingId, RingSet, UnpairedElectronsForm,
};

use crate::ops::model::{ElementScope, RingLimits};

#[derive(Clone, Debug)]
pub struct HueckelAromaticity {
    pub element_scope: ElementScope,
    pub ring_limits: RingLimits,
}

impl HueckelAromaticity {
    pub fn new(element_scope: ElementScope, ring_limits: RingLimits) -> Self {
        Self {
            element_scope,
            ring_limits,
        }
    }

    /// Whether some total contribution reachable from the members' ranges is
    /// perceived as aromatic by the rule: a candidate whose range admits no
    /// 4n+2 value is settled early.
    pub fn accepts_range(&self, members: &[(u32, u32)]) -> bool {
        let lower: u32 = members.iter().map(|&(lower, _)| lower).sum();
        let upper: u32 = members.iter().map(|&(_, upper)| upper).sum();
        let lower = lower.max(2);
        let first = lower + ((2 + 4 - lower % 4) % 4);
        first <= upper
    }

    pub fn find_from_rings<F>(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        electrons_at: &F,
    ) -> Vec<(Vec<AtomId>, AromaticSystemForm)>
    where
        F: Fn(AtomId) -> Option<u8>,
    {
        let eligible_cycles: Vec<RingId> = rings
            .ids()
            .filter(|&i| self.filter_ring(molecule, rings, i, electrons_at))
            .collect();

        let mut aromatic_atom_sets: Vec<HashSet<AtomId>> = Vec::new();

        for &cycle_idx in &eligible_cycles {
            let Some(ring) = rings.get(cycle_idx) else {
                continue;
            };
            let ring_atoms: Vec<AtomId> = ring.atoms().to_vec();
            if let Some(electrons) = ring_electron_count(&ring_atoms, electrons_at) {
                if check_4n_plus_2(electrons) {
                    aromatic_atom_sets.push(ring_atoms.into_iter().collect());
                }
            }
        }

        if self.ring_limits.include_unions {
            let fused_systems = self.enumerate_unions(rings, &eligible_cycles);
            for atoms in fused_systems {
                let atom_vec: Vec<AtomId> = atoms.iter().copied().collect();
                if let Some(electrons) = ring_electron_count(&atom_vec, electrons_at) {
                    if check_4n_plus_2(electrons) {
                        aromatic_atom_sets.push(atoms);
                    }
                }
            }
        }

        let merged = merge_overlapping_systems(aromatic_atom_sets);

        let mut candidates = Vec::new();
        for atom_set in merged {
            let mut atoms: Vec<AtomId> = atom_set.into_iter().collect();
            atoms.sort_unstable();

            let mut electrons: Vec<i64> = Vec::with_capacity(atoms.len());
            let mut valid = true;
            for &atom in &atoms {
                if let Some(e) = electrons_at(atom) {
                    electrons.push(i64::from(e));
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }

            candidates.push((
                atoms,
                AromaticSystemForm::from_electrons(electrons)
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            ));
        }

        candidates
    }

    fn is_atom_eligible<F>(&self, molecule: &Molecule, id: AtomId, electrons_at: &F) -> bool
    where
        F: Fn(AtomId) -> Option<u8>,
    {
        let view = molecule.atom(id);
        let element = match view.attributes.element {
            ElementForm::Lit(e) => e,
            _ => return false,
        };
        if !self.element_scope.contains(element) {
            return false;
        }
        electrons_at(id).is_some()
    }

    fn filter_ring<F>(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        ring: RingId,
        electrons_at: &F,
    ) -> bool
    where
        F: Fn(AtomId) -> Option<u8>,
    {
        let Some(ring) = rings.get(ring) else {
            return false;
        };
        let len = ring.len();
        if len < self.ring_limits.min_ring_size || len > self.ring_limits.max_ring_size {
            return false;
        }
        ring.atoms()
            .iter()
            .all(|&a| self.is_atom_eligible(molecule, a, electrons_at))
    }

    /// Connected ring unions for acceptance, walking rings that share at
    /// least one contiguous run of bonds — the `Fused` and `Bridged`
    /// relations. `Spiro` sharing breaks conjugation at the shared atom and
    /// `Noncontiguous` sharing is excluded until a need arises.
    pub(crate) fn enumerate_unions(
        &self,
        rings: &RingSet,
        eligible: &[RingId],
    ) -> Vec<HashSet<AtomId>> {
        let max_ring_count = self.ring_limits.max_ring_count;
        if max_ring_count < 2 {
            return Vec::new();
        }

        let eligible_set: HashSet<RingId> = eligible.iter().copied().collect();
        let mut results: Vec<HashSet<AtomId>> = Vec::new();
        let mut seen_unions: HashSet<Vec<RingId>> = HashSet::new();

        'outer: for &start in eligible {
            let mut stack: Vec<(Vec<RingId>, HashSet<AtomId>)> = Vec::new();
            let Some(start_ring) = rings.get(start) else {
                continue;
            };
            let start_atoms: HashSet<AtomId> = start_ring.atoms().iter().copied().collect();
            stack.push((vec![start], start_atoms));

            while let Some((union_rings, atoms)) = stack.pop() {
                if union_rings.len() >= 2 {
                    let mut key = union_rings.clone();
                    key.sort_unstable();
                    if seen_unions.insert(key) {
                        results.push(atoms.clone());
                        if results.len() >= self.ring_limits.max_unions {
                            break 'outer;
                        }
                    }
                }

                if union_rings.len() >= max_ring_count {
                    continue;
                }

                let last = *union_rings.last().unwrap();
                let mut neighbors = rings.fused_neighbors(last);
                neighbors.extend(rings.bridged_neighbors(last));
                neighbors.sort_unstable();
                for neighbor_idx in neighbors {
                    if !eligible_set.contains(&neighbor_idx) || union_rings.contains(&neighbor_idx)
                    {
                        continue;
                    }
                    if neighbor_idx <= union_rings[0] {
                        continue;
                    }
                    let mut extended_rings = union_rings.clone();
                    extended_rings.push(neighbor_idx);
                    let mut new_atoms = atoms.clone();
                    if let Some(nr) = rings.get(neighbor_idx) {
                        new_atoms.extend(nr.atoms().iter().copied());
                    }
                    stack.push((extended_rings, new_atoms));
                }
            }
        }

        results
    }
}

fn check_4n_plus_2(electron_count: u32) -> bool {
    if electron_count < 2 {
        return false;
    }
    (electron_count - 2).is_multiple_of(4)
}

fn ring_electron_count<F>(ids: &[AtomId], electrons_at: &F) -> Option<u32>
where
    F: Fn(AtomId) -> Option<u8>,
{
    let mut total: u32 = 0;
    for &id in ids {
        total += electrons_at(id)? as u32;
    }
    Some(total)
}

fn merge_overlapping_systems(aromatic_systems: Vec<HashSet<AtomId>>) -> Vec<HashSet<AtomId>> {
    if aromatic_systems.len() < 2 {
        return aromatic_systems;
    }

    let n = aromatic_systems.len();
    let mut components = UnionFind::new(n);
    let atom_capacity = aromatic_systems
        .iter()
        .flat_map(|system| system.iter())
        .map(|atom| atom.index() + 1)
        .max()
        .unwrap_or(0);
    let mut owner = vec![None; atom_capacity];

    for (system_index, system) in aromatic_systems.iter().enumerate() {
        for atom in system {
            match owner[atom.index()] {
                Some(previous) => components.union(system_index, previous),
                None => owner[atom.index()] = Some(system_index),
            }
        }
    }

    let mut result: Vec<HashSet<AtomId>> = Vec::new();
    let mut output_for_root: Vec<Option<usize>> = vec![None; n];
    for (system_index, system) in aromatic_systems.into_iter().enumerate() {
        let root = components.find(system_index);
        match output_for_root[root] {
            Some(output) => result[output].extend(system),
            None => {
                output_for_root[root] = Some(result.len());
                result.push(system);
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId, BondForm, ElectronCountsForm,
        ElementForm, Molecule, MoleculeEntries, NumForm, RingConfig, RingModel, RingSetKind,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    fn aromatic(element: Element, pi: i64) -> (AtomForm, Option<i64>) {
        (AtomForm::from_element(element), Some(pi))
    }

    fn aromatic_charged(element: Element, charge: i64, pi: i64) -> (AtomForm, Option<i64>) {
        (
            AtomForm {
                element: ElementForm::Lit(element),
                charge: NumForm::Lit(charge),
                ..Default::default()
            },
            Some(pi),
        )
    }

    fn plain(element: Element) -> (AtomForm, Option<i64>) {
        (AtomForm::from_element(element), None)
    }

    fn apply_pi(specs: Vec<(AtomForm, Option<i64>)>) -> Vec<AtomForm> {
        specs
            .into_iter()
            .map(|(mut atom, pi)| {
                if let Some(n) = pi {
                    atom.constraints.set(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(n)),
                    ));
                }
                atom
            })
            .collect()
    }

    fn make_ring(specs: Vec<(AtomForm, Option<i64>)>) -> Molecule {
        let n = specs.len();
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId(((i + 1) % n) as u32),
                    BondForm::from_order(1),
                )
            })
            .collect();
        let atoms = apply_pi(specs);
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn make_fused(specs: Vec<(AtomForm, Option<i64>)>, edges: &[(usize, usize)]) -> Molecule {
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomId(a as u32), AtomId(b as u32), BondForm::from_order(1)))
            .collect();
        let atoms = apply_pi(specs);
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn daylight_model() -> HueckelAromaticity {
        HueckelAromaticity::new(
            ElementScope::AllowList(vec![
                Element::C,
                Element::N,
                Element::O,
                Element::S,
                Element::Se,
                Element::As,
            ]),
            RingLimits::default(),
        )
    }

    fn mdl_model() -> HueckelAromaticity {
        HueckelAromaticity::new(
            ElementScope::AllowList(vec![Element::C, Element::N]),
            RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        )
    }

    fn permissive_model() -> HueckelAromaticity {
        HueckelAromaticity::new(ElementScope::Any, RingLimits::default())
    }

    #[fixture]
    fn benzene() -> Molecule {
        make_ring(vec![aromatic(Element::C, 1); 6])
    }

    #[fixture]
    fn pyridine() -> Molecule {
        make_ring(vec![
            aromatic(Element::N, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn pyrrole() -> Molecule {
        make_ring(vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn furan() -> Molecule {
        make_ring(vec![
            aromatic(Element::O, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn thiophene() -> Molecule {
        make_ring(vec![
            aromatic(Element::S, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn imidazole() -> Molecule {
        make_ring(vec![
            aromatic(Element::N, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn tropylium() -> Molecule {
        make_ring(vec![
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic_charged(Element::C, 1, 0),
        ])
    }

    #[fixture]
    fn cyclopentadienyl_anion() -> Molecule {
        make_ring(vec![
            aromatic_charged(Element::C, -1, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> Molecule {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn azulene() -> Molecule {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    fn phenanthrene() -> Molecule {
        make_fused(
            vec![aromatic(Element::C, 1); 14],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (3, 6), (6, 7),
                (7, 8), (8, 9), (9, 4), (8, 10), (10, 11), (11, 12), (12, 13), (13, 9),
            ],
        )
    }

    #[fixture]
    fn cyclobutadiene() -> Molecule {
        make_ring(vec![aromatic(Element::C, 1); 4])
    }

    #[fixture]
    fn cyclohexane() -> Molecule {
        make_ring(vec![plain(Element::C); 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn cubane() -> Molecule {
        make_fused(
            vec![plain(Element::C); 8],
            &[
                (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
            ],
        )
    }

    #[fixture]
    fn borazine() -> Molecule {
        make_ring(vec![
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
        ])
    }

    fn electron_total(system: &(Vec<AtomId>, AromaticSystemForm)) -> i64 {
        match &system.1.electrons {
            ElectronCountsForm::Lit(counts) => counts.iter().sum(),
            ElectronCountsForm::Undetermined => 0,
        }
    }

    #[rstest]
    #[case::point_accepted(&[(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)], true)]
    #[case::point_rejected(&[(1, 1), (1, 1), (1, 1), (1, 1), (1, 1)], false)]
    #[case::span_accepted(&[(0, 3), (1, 1), (1, 1), (1, 1), (1, 1)], true)]
    #[case::span_rejected(&[(3, 4), (1, 1), (1, 1), (1, 1), (1, 1)], false)]
    #[case::two_electrons(&[(2, 2)], true)]
    #[case::below_two(&[(0, 1)], false)]
    #[case::empty(&[], false)]
    fn test_hueckel_aromaticity_aromaticity_accepts_range(
        #[case] members: &[(u32, u32)],
        #[case] expected: bool,
    ) {
        assert_eq!(daylight_model().accepts_range(members), expected);
    }

    #[rstest]
    #[case::benzene(benzene(), 6, 6)]
    #[case::pyridine(pyridine(), 6, 6)]
    #[case::pyrrole(pyrrole(), 5, 6)]
    #[case::furan(furan(), 5, 6)]
    #[case::thiophene(thiophene(), 5, 6)]
    #[case::imidazole(imidazole(), 5, 6)]
    #[case::naphthalene(naphthalene(), 10, 10)]
    #[case::azulene(azulene(), 10, 10)]
    #[case::phenanthrene(phenanthrene(), 14, 14)]
    #[case::tropylium(tropylium(), 7, 6)]
    #[case::cyclopentadienyl_anion(cyclopentadienyl_anion(), 5, 6)]
    fn test_hueckel_aromaticity_find_from_rings_aromatic(
        #[case] molecule: Molecule,
        #[case] expected_atoms: usize,
        #[case] expected_electrons: i64,
    ) {
        let rings = molecule
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: RingLimits::default().max_ring_size,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let model = daylight_model();
        let systems = model.find_from_rings(&molecule, &rings, &|v| match molecule
            .atom(v)
            .attributes
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceForm::Undetermined)
        {
            AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
            _ => None,
        });
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].0.len(), expected_atoms);
        assert_eq!(electron_total(&systems[0]), expected_electrons);
    }

    #[rstest]
    #[case::cyclobutadiene(cyclobutadiene(), daylight_model())]
    #[case::cyclohexane(cyclohexane(), daylight_model())]
    #[case::cubane(cubane(), daylight_model())]
    #[case::borazine_daylight(borazine(), daylight_model())]
    #[case::pyrrole_mdl(pyrrole(), mdl_model())]
    fn test_hueckel_aromaticity_find_from_rings_non_aromatic(
        #[case] molecule: Molecule,
        #[case] model: HueckelAromaticity,
    ) {
        let rings = molecule
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: RingLimits::default().max_ring_size,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let systems = model.find_from_rings(&molecule, &rings, &|v| match molecule
            .atom(v)
            .attributes
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceForm::Undetermined)
        {
            AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
            _ => None,
        });
        assert!(systems.is_empty());
    }

    #[rstest]
    fn test_hueckel_aromaticity_find_from_rings_bridged() {
        // Two five-rings sharing the two-bond run 3-4-0: the `Bridged`
        // relation, not `Fused`. Neither ring passes alone (sum 4); their
        // union does (sum 6), so the combination walk must cross bridged
        // sharing.
        let molecule = mol_dsl!(
            r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a0" "C#a" "C#a"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]
                        [0 5 "1"] [5 6 "1"] [6 3 "1"]]}"#
        );
        let rings = molecule
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: RingLimits::default().max_ring_size,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let systems = daylight_model().find_from_rings(&molecule, &rings, &|v| match molecule
            .atom(v)
            .attributes
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceForm::Undetermined)
        {
            AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
            _ => None,
        });

        assert_eq!(
            systems,
            vec![(
                (0..7).map(AtomId).collect::<Vec<_>>(),
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 0, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )]
        );
    }

    #[rstest]
    fn test_hueckel_aromaticity_find_from_rings_borazine_permissive(borazine: Molecule) {
        let rings = borazine
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: RingLimits::default().max_ring_size,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let systems = permissive_model().find_from_rings(&borazine, &rings, &|v| match borazine
            .atom(v)
            .attributes
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceForm::Undetermined)
        {
            AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
            _ => None,
        });
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].0.len(), 6);
        assert_eq!(electron_total(&systems[0]), 6);
    }

    #[rstest]
    #[case::empty(vec![], vec![])]
    #[case::single(
        vec![HashSet::from([AtomId(0), AtomId(1)])],
        vec![vec![AtomId(0), AtomId(1)]]
    )]
    #[case::disjoint(
        vec![
            HashSet::from([AtomId(0), AtomId(1)]),
            HashSet::from([AtomId(3), AtomId(4)]),
        ],
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(3), AtomId(4)],
        ]
    )]
    #[case::overlapping(
        vec![
            HashSet::from([AtomId(0), AtomId(1)]),
            HashSet::from([AtomId(1), AtomId(2)]),
        ],
        vec![vec![AtomId(0), AtomId(1), AtomId(2)]]
    )]
    #[case::transitive(
        vec![
            HashSet::from([AtomId(0), AtomId(1)]),
            HashSet::from([AtomId(3), AtomId(4)]),
            HashSet::from([AtomId(1), AtomId(3)]),
        ],
        vec![vec![
            AtomId(0),
            AtomId(1),
            AtomId(3),
            AtomId(4),
        ]]
    )]
    #[case::two_components(
        vec![
            HashSet::from([AtomId(0), AtomId(1)]),
            HashSet::from([AtomId(1), AtomId(2)]),
            HashSet::from([AtomId(5), AtomId(6)]),
            HashSet::from([AtomId(6), AtomId(7)]),
        ],
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![AtomId(5), AtomId(6), AtomId(7)],
        ]
    )]
    fn test_merge_overlapping_systems(
        #[case] aromatic_systems: Vec<HashSet<AtomId>>,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let mut actual: Vec<Vec<AtomId>> = merge_overlapping_systems(aromatic_systems)
            .into_iter()
            .map(|system| {
                let mut atoms: Vec<AtomId> = system.into_iter().collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        actual.sort_unstable();

        assert_eq!(actual, expected);
    }
}
