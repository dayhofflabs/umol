//! HueckelRule (4n+2 electron counting) aromaticity perception.
//!
//! Filters candidate atoms by element scope and the per-atom aromatic-valence
//! constraint, enumerates rings within configured bounds, checks the Hueckel
//! 4n+2 rule on individual and fused ring combinations, and produces aromatic
//! system tuples `(Vec<AtomId>, AromaticSystemForm)` ready for `Molecule::edit`.

use std::collections::{HashMap, HashSet};

use umol_graph_core::UnionFind;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomId, AtomView, ElementForm, Molecule, RingId, RingSet, RingView,
    UnpairedElectronsForm,
};

use crate::ops::model::{ElementScope, RingLimits};

#[derive(Clone, Debug)]
pub struct HueckelRuleAromaticity {
    pub element_scope: ElementScope,
    pub ring_limits: RingLimits,
}

impl HueckelRuleAromaticity {
    pub fn new(element_scope: ElementScope, ring_limits: RingLimits) -> Self {
        Self {
            element_scope,
            ring_limits,
        }
    }

    pub fn find_from_rings<F>(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        electrons_at: &F,
    ) -> Vec<(Vec<AtomId>, AromaticSystemForm)>
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let eligible_cycles: Vec<RingId> = rings
            .ids()
            .filter(|&i| {
                rings
                    .get(i)
                    .is_some_and(|r| self.filter_ring(molecule, r, electrons_at))
            })
            .collect();

        let mut aromatic_atom_sets: Vec<HashSet<AtomId>> = Vec::new();

        for &cycle_idx in &eligible_cycles {
            let Some(ring) = rings.get(cycle_idx) else {
                continue;
            };
            let ring_atoms: Vec<AtomId> = ring.atoms().to_vec();
            if let Some(electrons) = ring_electron_count(molecule, &ring_atoms, electrons_at) {
                if check_4n_plus_2(electrons) {
                    aromatic_atom_sets.push(ring_atoms.into_iter().collect());
                }
            }
        }

        if self.ring_limits.include_fused {
            let fused_systems = self.enumerate_fused_combinations(rings, &eligible_cycles);
            for atoms in fused_systems {
                let atom_vec: Vec<AtomId> = atoms.iter().copied().collect();
                if let Some(electrons) = ring_electron_count(molecule, &atom_vec, electrons_at) {
                    if check_4n_plus_2(electrons) {
                        aromatic_atom_sets.push(atoms);
                    }
                }
            }
        }

        let merged = merge_overlapping_systems(&aromatic_atom_sets);

        let mut candidates = Vec::new();
        for atom_set in merged {
            let mut atoms: Vec<AtomId> = atom_set.into_iter().collect();
            atoms.sort_unstable();

            let mut electrons: Vec<i64> = Vec::with_capacity(atoms.len());
            let mut valid = true;
            for &atom in &atoms {
                if let Some(e) = electrons_at(&molecule.atom(atom)) {
                    electrons.push(e as i64);
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
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let view = molecule.atom(id);
        let element = match view.attributes.element {
            ElementForm::Lit(e) => e,
            _ => return false,
        };
        if !self.element_scope.contains(element) {
            return false;
        }
        electrons_at(&view).is_some()
    }

    fn filter_ring<F>(&self, molecule: &Molecule, ring: RingView<'_>, electrons_at: &F) -> bool
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let len = ring.len();
        if len < self.ring_limits.min_ring_size || len > self.ring_limits.max_ring_size {
            return false;
        }
        ring.atoms()
            .iter()
            .all(|&a| self.is_atom_eligible(molecule, a, electrons_at))
    }

    fn enumerate_fused_combinations(
        &self,
        rings: &RingSet,
        eligible: &[RingId],
    ) -> Vec<HashSet<AtomId>> {
        let max_combo = self.ring_limits.max_fused_combination;
        if max_combo < 2 {
            return Vec::new();
        }

        let eligible_set: HashSet<RingId> = eligible.iter().copied().collect();
        let mut results: Vec<HashSet<AtomId>> = Vec::new();
        let mut seen_combos: HashSet<Vec<RingId>> = HashSet::new();

        'outer: for &start in eligible {
            let mut stack: Vec<(Vec<RingId>, HashSet<AtomId>)> = Vec::new();
            let Some(start_ring) = rings.get(start) else {
                continue;
            };
            let start_atoms: HashSet<AtomId> = start_ring.atoms().iter().copied().collect();
            stack.push((vec![start], start_atoms));

            while let Some((combo, atoms)) = stack.pop() {
                if combo.len() >= 2 {
                    let mut key = combo.clone();
                    key.sort_unstable();
                    if seen_combos.insert(key) {
                        results.push(atoms.clone());
                        if results.len() >= self.ring_limits.max_fused_search {
                            break 'outer;
                        }
                    }
                }

                if combo.len() >= max_combo {
                    continue;
                }

                let last = *combo.last().unwrap();
                for neighbor_idx in rings.fused_neighbors(last) {
                    if !eligible_set.contains(&neighbor_idx) || combo.contains(&neighbor_idx) {
                        continue;
                    }
                    if neighbor_idx <= combo[0] {
                        continue;
                    }
                    let mut new_combo = combo.clone();
                    new_combo.push(neighbor_idx);
                    let mut new_atoms = atoms.clone();
                    if let Some(nr) = rings.get(neighbor_idx) {
                        new_atoms.extend(nr.atoms().iter().copied());
                    }
                    stack.push((new_combo, new_atoms));
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

fn ring_electron_count<F>(molecule: &Molecule, ids: &[AtomId], electrons_at: &F) -> Option<u32>
where
    F: Fn(&AtomView<'_>) -> Option<u8>,
{
    let mut total: u32 = 0;
    for &id in ids {
        total += electrons_at(&molecule.atom(id))? as u32;
    }
    Some(total)
}

fn merge_overlapping_systems(aromatic_systems: &[HashSet<AtomId>]) -> Vec<HashSet<AtomId>> {
    if aromatic_systems.is_empty() {
        return Vec::new();
    }

    let n = aromatic_systems.len();
    let mut uf = UnionFind::new(n);

    for i in 0..n {
        for j in (i + 1)..n {
            if !aromatic_systems[i].is_disjoint(&aromatic_systems[j]) {
                uf.union(i, j);
            }
        }
    }

    // TODO: Check if usize-typed indices are correct.
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        groups.entry(uf.find(i)).or_default().push(i);
    }

    let mut result = Vec::new();
    for (_, indices) in groups {
        let mut merged_atoms: HashSet<AtomId> = HashSet::new();
        for &idx in &indices {
            merged_atoms.extend(aromatic_systems[idx].iter());
        }
        result.push(merged_atoms);
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

    fn daylight_model() -> HueckelRuleAromaticity {
        HueckelRuleAromaticity::new(
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

    fn mdl_model() -> HueckelRuleAromaticity {
        HueckelRuleAromaticity::new(
            ElementScope::AllowList(vec![Element::C, Element::N]),
            RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        )
    }

    fn permissive_model() -> HueckelRuleAromaticity {
        HueckelRuleAromaticity::new(ElementScope::Any, RingLimits::default())
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
    fn test_hueckel_rule_find_from_rings_aromatic(
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
        let systems = model.find_from_rings(&molecule, &rings, &|v| match v
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
    fn test_hueckel_rule_find_from_rings_non_aromatic(
        #[case] molecule: Molecule,
        #[case] model: HueckelRuleAromaticity,
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
        let systems = model.find_from_rings(&molecule, &rings, &|v| match v
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
    fn test_hueckel_rule_find_from_rings_borazine_permissive(borazine: Molecule) {
        let rings = borazine
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: RingLimits::default().max_ring_size,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let systems = permissive_model().find_from_rings(&borazine, &rings, &|v| match v
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
}
