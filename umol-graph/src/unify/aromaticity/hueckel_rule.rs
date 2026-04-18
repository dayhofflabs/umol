//! HueckelRule (4n+2 electron counting) aromaticity model.
//!
//! Filters candidate atoms by element scope and aromatic valence, enumerates
//! rings within configured bounds, checks the Hueckel 4n+2 rule on individual
//! and fused ring combinations, and produces `AromaticSystem` objects.

use std::collections::{HashMap, HashSet};

use umol_graph_core::UnionFind;

use umol_shared::atom_ast::ElementAst;
use umol_shared::value_ast::ValueAst;

use super::{ElementScope, RingLimits};
use crate::ast::AtomIdx;
use crate::ast::aromatic::AromaticSystem;
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, BondConstraint};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{RingIdx, RingSet, RingView};

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

    pub fn find_from_rings(
        &self,
        ast: &MoleculeAst,
        rings: &RingSet,
    ) -> Vec<AromaticSystem> {
        let eligible_cycles: Vec<RingIdx> = rings
            .ids()
            .filter(|&i| rings.get(i).is_some_and(|r| self.filter_ring(ast, r)))
            .collect();

        let mut aromatic_atom_sets: Vec<HashSet<AtomIdx>> = Vec::new();

        for &cycle_idx in &eligible_cycles {
            let Some(ring) = rings.get(cycle_idx) else {
                continue;
            };
            if let Some(electrons) = self.ring_electron_count(ast, ring.atoms()) {
                if self.check_4n_plus_2(electrons) {
                    let atom_set: HashSet<AtomIdx> = ring.atoms().iter().copied().collect();
                    aromatic_atom_sets.push(atom_set);
                }
            }
        }

        if self.ring_limits.include_fused {
            let fused_systems = self.enumerate_fused_combinations(rings, &eligible_cycles);
            for atoms in fused_systems {
                let atom_vec: Vec<AtomIdx> = atoms.iter().copied().collect();
                if let Some(electrons) = self.ring_electron_count(ast, &atom_vec) {
                    if self.check_4n_plus_2(electrons) {
                        aromatic_atom_sets.push(atoms);
                    }
                }
            }
        }

        let merged = merge_overlapping_systems(&aromatic_atom_sets);

        let mut candidates = Vec::new();
        for atom_set in merged {
            let mut atoms: Vec<AtomIdx> = atom_set.into_iter().collect();
            atoms.sort_unstable();

            let mut atom_constraints: Vec<AtomConstraint> = Vec::with_capacity(atoms.len());
            let mut valid = true;
            for &atom in &atoms {
                if let Some(e) = self.aromatic_electron_count(ast, atom) {
                    atom_constraints.push(AtomConstraint::AromaticValence(
                        AromaticValenceConstraint::Value(ValueAst::Lit(e as i64)),
                    ));
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }

            let bonds = ast.induced_bonds(&atoms);
            let bond_constraints = vec![BondConstraint::Aromatic; bonds.len()];

            candidates.push(AromaticSystem::new(
                atoms,
                bonds,
                AromaticSystemAst::default(),
                atom_constraints,
                bond_constraints,
            ));
        }

        candidates
    }

    fn is_atom_eligible(&self, ast: &MoleculeAst, atom: AtomIdx) -> bool {
        let atom_ast = ast.atom(atom);
        let element = match atom_ast.data.element {
            ElementAst::Lit(e) => e,
            _ => return false,
        };
        match &self.element_scope {
            ElementScope::Any => {}
            ElementScope::AllowList(allowed) => {
                if !allowed.contains(&element) {
                    return false;
                }
            }
        }
        ast.atom_aromatic_pi_electrons(atom).is_some()
    }

    fn aromatic_electron_count(&self, ast: &MoleculeAst, atom: AtomIdx) -> Option<u8> {
        ast.atom_aromatic_pi_electrons(atom)
    }

    fn filter_ring(&self, ast: &MoleculeAst, ring: RingView<'_>) -> bool {
        let len = ring.len();
        if len < self.ring_limits.min_ring_size || len > self.ring_limits.max_ring_size {
            return false;
        }
        ring.atoms().iter().all(|&a| self.is_atom_eligible(ast, a))
    }

    fn check_4n_plus_2(&self, electron_count: u32) -> bool {
        if electron_count < 2 {
            return false;
        }
        (electron_count - 2).is_multiple_of(4)
    }

    fn ring_electron_count(&self, ast: &MoleculeAst, atoms: &[AtomIdx]) -> Option<u32> {
        let mut total: u32 = 0;
        for &atom in atoms {
            total += self.aromatic_electron_count(ast, atom)? as u32;
        }
        Some(total)
    }

    fn enumerate_fused_combinations(
        &self,
        rings: &RingSet,
        eligible: &[RingIdx],
    ) -> Vec<HashSet<AtomIdx>> {
        let max_combo = self.ring_limits.max_fused_combination;
        if max_combo < 2 {
            return Vec::new();
        }

        let eligible_set: HashSet<RingIdx> = eligible.iter().copied().collect();
        let mut results: Vec<HashSet<AtomIdx>> = Vec::new();
        let mut seen_combos: HashSet<Vec<RingIdx>> = HashSet::new();

        'outer: for &start in eligible {
            let mut stack: Vec<(Vec<RingIdx>, HashSet<AtomIdx>)> = Vec::new();
            let Some(start_ring) = rings.get(start) else {
                continue;
            };
            let start_atoms: HashSet<AtomIdx> = start_ring.atoms().iter().copied().collect();
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

fn merge_overlapping_systems(
    aromatic_systems: &[HashSet<AtomIdx>],
) -> Vec<HashSet<AtomIdx>> {
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

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        groups.entry(uf.find(i)).or_default().push(i);
    }

    let mut result = Vec::new();
    for (_, members) in groups {
        let mut merged_atoms: HashSet<AtomIdx> = HashSet::new();
        for &idx in &members {
            merged_atoms.extend(aromatic_systems[idx].iter());
        }
        result.push(merged_atoms);
    }

    result
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::atom_ast::ElementAst;
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::AtomIdx;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{
        AromaticValenceConstraint, AtomConstraint, MoleculeConstraint,
    };
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::rings::RingEnumerationStrategy;
    use crate::ast::rings::{RingEnumerator, RingFamily};

    fn aromatic(element: Element, pi: i64) -> (AtomAst, Option<i64>) {
        (AtomAst::from_element(element), Some(pi))
    }

    fn aromatic_charged(element: Element, charge: i64, pi: i64) -> (AtomAst, Option<i64>) {
        (
            AtomAst {
                element: ElementAst::Lit(element),
                charge: ValueAst::Lit(charge),
                ..Default::default()
            },
            Some(pi),
        )
    }

    fn plain(element: Element) -> (AtomAst, Option<i64>) {
        (AtomAst::from_element(element), None)
    }

    fn pi_constraints(specs: &[(AtomAst, Option<i64>)]) -> Vec<MoleculeConstraint> {
        specs
            .iter()
            .enumerate()
            .filter_map(|(i, (_, pi))| {
                pi.map(|n| {
                    MoleculeConstraint::AtomPred(
                        AtomIdx(i as u32),
                        AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                            ValueAst::Lit(n),
                        )),
                    )
                })
            })
            .collect()
    }

    fn make_ring(specs: Vec<(AtomAst, Option<i64>)>) -> MoleculeAst {
        let n = specs.len();
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        let constraints = pi_constraints(&specs);
        let atoms: Vec<AtomAst> = specs.into_iter().map(|(a, _)| a).collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], constraints)
    }

    fn make_fused(specs: Vec<(AtomAst, Option<i64>)>, edges: &[(usize, usize)]) -> MoleculeAst {
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(1)))
            .collect();
        let constraints = pi_constraints(&specs);
        let atoms: Vec<AtomAst> = specs.into_iter().map(|(a, _)| a).collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], constraints)
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
    fn benzene() -> MoleculeAst {
        make_ring(vec![aromatic(Element::C, 1); 6])
    }

    #[fixture]
    fn pyridine() -> MoleculeAst {
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
    fn pyrrole() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn furan() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::O, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn thiophene() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::S, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn imidazole() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::N, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn tropylium() -> MoleculeAst {
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
    fn cyclopentadienyl_anion() -> MoleculeAst {
        make_ring(vec![aromatic_charged(Element::C, -1, 2); 5])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeAst {
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
    fn azulene() -> MoleculeAst {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    fn phenanthrene() -> MoleculeAst {
        make_fused(
            vec![aromatic(Element::C, 1); 14],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (3, 6), (6, 7),
                (7, 8), (8, 9), (9, 4), (8, 10), (10, 11), (11, 12), (12, 13), (13, 9),
            ],
        )
    }

    #[fixture]
    fn cyclobutadiene() -> MoleculeAst {
        make_ring(vec![aromatic(Element::C, 1); 4])
    }

    #[fixture]
    fn cyclohexane() -> MoleculeAst {
        make_ring(vec![plain(Element::C); 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn cubane() -> MoleculeAst {
        make_fused(
            vec![plain(Element::C); 8],
            &[
                (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
            ],
        )
    }

    #[fixture]
    fn borazine() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
            aromatic(Element::B, 0),
            aromatic(Element::N, 2),
        ])
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
    #[case::cyclopentadienyl_anion(cyclopentadienyl_anion(), 5, 10)]
    fn test_find_from_rings_aromatic(
        #[case] ast: MoleculeAst,
        #[case] expected_atoms: usize,
        #[case] expected_electrons: u8,
    ) {
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate(&ast);
        let model = daylight_model();
        let systems = model.find_from_rings(&ast, &rings);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atom_count(), expected_atoms);
        assert_eq!(systems[0].electron_count(), expected_electrons);
    }

    #[rstest]
    #[case::cyclobutadiene(cyclobutadiene(), daylight_model())]
    #[case::cyclohexane(cyclohexane(), daylight_model())]
    #[case::cubane(cubane(), daylight_model())]
    #[case::borazine_daylight(borazine(), daylight_model())]
    #[case::pyrrole_mdl(pyrrole(), mdl_model())]
    fn test_find_from_rings_non_aromatic(
        #[case] ast: MoleculeAst,
        #[case] model: HueckelRuleAromaticity,
    ) {
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate(&ast);
        let systems = model.find_from_rings(&ast, &rings);
        assert!(systems.is_empty());
    }

    #[rstest]
    fn test_find_from_rings_borazine_permissive(borazine: MoleculeAst) {
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate(&borazine);
        let systems = permissive_model().find_from_rings(&borazine, &rings);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atom_count(), 6);
        assert_eq!(systems[0].electron_count(), 6);
    }
}
