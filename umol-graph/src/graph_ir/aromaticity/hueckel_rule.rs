//! HueckelRule (4n+2 electron counting) aromaticity model.
//!
//! Filters candidate atoms by element scope and aromatic valence, enumerates
//! rings within configured bounds, checks the Hueckel 4n+2 rule on individual
//! and fused ring combinations, and produces `AromaticSystem` objects.

use std::collections::{HashMap, HashSet};

use petgraph::unionfind::UnionFind;

use super::{AromaticContribution, AromaticSystem};
use crate::graph_ir::config::{ElementScope, RingLimits};
use crate::graph_ir::molecule::AtomIndex;
use crate::graph_ir::molecule_builder::MoleculeBuilder;
use crate::graph_ir::rings::{Ring, RingIndex, RingSet};

/// HueckelRule aromaticity model. Parameterized by element scope and ring scope.
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
        builder: &MoleculeBuilder,
        rings: &RingSet,
    ) -> Vec<AromaticSystem> {
        let eligible_cycles: Vec<RingIndex> = rings
            .ring_indices()
            .filter(|&i| rings.ring(i).is_some_and(|r| self.filter_ring(builder, r)))
            .collect();

        let mut aromatic_atom_sets: Vec<(HashSet<AtomIndex>, Vec<Ring>)> = Vec::new();

        for &cycle_idx in &eligible_cycles {
            let Some(ring) = rings.ring(cycle_idx) else {
                continue;
            };
            if let Some(electrons) = self.ring_electron_count(builder, ring.atoms()) {
                if self.check_4n_plus_2(electrons) {
                    let atom_set: HashSet<AtomIndex> = ring.atoms().iter().copied().collect();
                    aromatic_atom_sets.push((atom_set, vec![ring.clone()]));
                }
            }
        }

        // Check fused ring combinations if enabled.
        if self.ring_limits.include_fused {
            let fused_systems = self.enumerate_fused_combinations(rings, &eligible_cycles);
            for (atoms, rings) in fused_systems {
                let atom_vec: Vec<AtomIndex> = atoms.iter().copied().collect();
                if let Some(electrons) = self.ring_electron_count(builder, &atom_vec) {
                    if self.check_4n_plus_2(electrons) {
                        aromatic_atom_sets.push((atoms, rings));
                    }
                }
            }
        }

        // Merge overlapping aromatic systems: if two systems share atoms,
        // combine them into a single system.
        let merged = merge_overlapping_systems(&aromatic_atom_sets);

        // Build candidates from merged systems.
        let mut candidates = Vec::new();
        for (atom_set, rings) in merged {
            let mut contributions: Vec<AromaticContribution> = Vec::new();
            let mut valid = true;
            for &atom in &atom_set {
                if let Some(e) = self.aromatic_electron_count(builder, atom) {
                    contributions.push(AromaticContribution::new(atom, e));
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }
            candidates.push(AromaticSystem::with_rings(contributions, rings));
        }

        candidates
    }

    fn is_atom_eligible(&self, builder: &MoleculeBuilder, atom: AtomIndex) -> bool {
        let atom_data = match builder.atom(atom) {
            Some(a) => a,
            None => return false,
        };
        match &self.element_scope {
            ElementScope::Any => {}
            ElementScope::AllowList(allowed) => {
                if !allowed.contains(&atom_data.element()) {
                    return false;
                }
            }
        }
        builder.atom_has_aromatic_candidate(atom)
    }

    fn aromatic_electron_count(&self, builder: &MoleculeBuilder, atom: AtomIndex) -> Option<u8> {
        builder.atom(atom)?;
        Some(builder.atom_aromatic_valence(atom))
    }

    fn filter_ring(&self, builder: &MoleculeBuilder, ring: &Ring) -> bool {
        let len = ring.len();
        if len < self.ring_limits.min_ring_size || len > self.ring_limits.max_ring_size {
            return false;
        }
        ring.atoms()
            .iter()
            .all(|&a| self.is_atom_eligible(builder, a))
    }

    fn check_4n_plus_2(&self, electron_count: u32) -> bool {
        if electron_count < 2 {
            return false;
        }
        (electron_count - 2).is_multiple_of(4)
    }

    fn ring_electron_count(&self, builder: &MoleculeBuilder, atoms: &[AtomIndex]) -> Option<u32> {
        let mut total: u32 = 0;
        for &atom in atoms {
            total += self.aromatic_electron_count(builder, atom)? as u32;
        }
        Some(total)
    }

    /// Enumerate fused ring combinations by iteratively merging rings that
    /// share at least one bond. Uses the `fused_neighbors` relationship from
    /// `RingSet` and bounds exploration to `max_fused_combination`.
    fn enumerate_fused_combinations(
        &self,
        rings: &RingSet,
        eligible: &[RingIndex],
    ) -> Vec<(HashSet<AtomIndex>, Vec<Ring>)> {
        let max_combo = self.ring_limits.max_fused_combination;
        if max_combo < 2 {
            return Vec::new();
        }

        let eligible_set: HashSet<RingIndex> = eligible.iter().copied().collect();
        let mut results: Vec<(HashSet<AtomIndex>, Vec<Ring>)> = Vec::new();
        let mut seen_combos: HashSet<Vec<RingIndex>> = HashSet::new();

        'outer: for &start in eligible {
            let mut stack: Vec<(Vec<RingIndex>, HashSet<AtomIndex>)> = Vec::new();
            let Some(start_ring) = rings.ring(start) else {
                continue;
            };
            let start_atoms: HashSet<AtomIndex> = start_ring.atoms().iter().copied().collect();
            stack.push((vec![start], start_atoms));

            while let Some((combo, atoms)) = stack.pop() {
                if combo.len() >= 2 {
                    let mut key = combo.clone();
                    key.sort_unstable();
                    if seen_combos.insert(key) {
                        let rings: Vec<Ring> = combo
                            .iter()
                            .filter_map(|&i| rings.ring(i).cloned())
                            .collect();
                        results.push((atoms.clone(), rings));
                        if results.len() >= self.ring_limits.max_fused_search {
                            break 'outer;
                        }
                    }
                }

                if combo.len() >= max_combo {
                    continue;
                }

                let last = *combo.last().unwrap();
                for neighbor_idx in rings.ring_fused_neighbors(last) {
                    if !eligible_set.contains(&neighbor_idx) || combo.contains(&neighbor_idx) {
                        continue;
                    }
                    if neighbor_idx <= combo[0] {
                        continue;
                    }
                    let mut new_combo = combo.clone();
                    new_combo.push(neighbor_idx);
                    let mut new_atoms = atoms.clone();
                    if let Some(nr) = rings.ring(neighbor_idx) {
                        new_atoms.extend(nr.atoms().iter().copied());
                    }
                    stack.push((new_combo, new_atoms));
                }
            }
        }

        results
    }
}

/// Merge overlapping aromatic systems. Two systems overlap if they share
/// any atoms. Returns deduplicated, merged results.
fn merge_overlapping_systems(
    aromatic_systems: &[(HashSet<AtomIndex>, Vec<Ring>)],
) -> Vec<(HashSet<AtomIndex>, Vec<Ring>)> {
    if aromatic_systems.is_empty() {
        return Vec::new();
    }

    let n = aromatic_systems.len();
    let mut uf = UnionFind::<usize>::new(n);

    for i in 0..n {
        for j in (i + 1)..n {
            if !aromatic_systems[i].0.is_disjoint(&aromatic_systems[j].0) {
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
        let mut merged_atoms: HashSet<AtomIndex> = HashSet::new();
        let mut merged_rings: Vec<Ring> = Vec::new();
        let mut seen_rings: HashSet<Vec<AtomIndex>> = HashSet::new();
        for &idx in &members {
            merged_atoms.extend(aromatic_systems[idx].0.iter());
            for ring in &aromatic_systems[idx].1 {
                let mut normalized = ring.atoms().to_vec();
                normalized.sort_unstable();
                if seen_rings.insert(normalized) {
                    merged_rings.push(ring.clone());
                }
            }
        }
        result.push((merged_atoms, merged_rings));
    }

    result
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::graph_ir::bond_pattern::BondPattern;
    use crate::graph_ir::config::RingEnumerationStrategy;
    use crate::graph_ir::rings::{RingEnumerator, RingFamily};

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

    fn make_ring(atom_specs: &[&str]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_resolved_atom(s.parse().unwrap()))
            .collect();
        let n = atoms.len();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondPattern::new(1));
        }
        builder
    }

    fn make_fused(atom_specs: &[&str], edges: &[(usize, usize)]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_resolved_atom(s.parse().unwrap()))
            .collect();
        for &(a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    const C1: &str = "C#h#v2#a";
    const C0: &str = "C#v4";

    #[fixture]
    fn benzene() -> MoleculeBuilder {
        make_ring(&[C1; 6])
    }

    #[fixture]
    fn pyridine() -> MoleculeBuilder {
        make_ring(&["N#n#v2#a", C1, C1, C1, C1, C1])
    }

    #[fixture]
    fn pyrrole() -> MoleculeBuilder {
        make_ring(&["N#h#v2#a2", C1, C1, C1, C1])
    }

    #[fixture]
    fn furan() -> MoleculeBuilder {
        make_ring(&["O#n#v2#a2", C1, C1, C1, C1])
    }

    #[fixture]
    fn thiophene() -> MoleculeBuilder {
        make_ring(&["S#n#v2#a2", C1, C1, C1, C1])
    }

    #[fixture]
    fn imidazole() -> MoleculeBuilder {
        make_ring(&["N#n#v2#a", C1, C1, "N#h#v2#a2", C1])
    }

    #[fixture]
    // 6 neutral C (a=1) + 1 C+ (a=0), charge-separated representation
    fn tropylium() -> MoleculeBuilder {
        make_ring(&[C1, C1, C1, C1, C1, C1, "C#c+#h#v2#a0"])
    }

    #[fixture]
    fn cyclopentadienyl_anion() -> MoleculeBuilder {
        make_ring(&["C#c-#h#v2#a2"; 5])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeBuilder {
        make_fused(
            &[C1; 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn azulene() -> MoleculeBuilder {
        make_fused(
            &[C1; 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    fn phenanthrene() -> MoleculeBuilder {
        make_fused(
            &[C1; 14],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (3, 6), (6, 7),
                (7, 8), (8, 9), (9, 4), (8, 10), (10, 11), (11, 12), (12, 13), (13, 9),
            ],
        )
    }

    #[fixture]
    fn cyclobutadiene() -> MoleculeBuilder {
        make_ring(&[C1; 4])
    }

    #[fixture]
    fn cyclohexane() -> MoleculeBuilder {
        make_ring(&[C0; 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn cubane() -> MoleculeBuilder {
        make_fused(
            &[C0; 8],
            &[
                (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
            ],
        )
    }

    #[fixture]
    fn borazine() -> MoleculeBuilder {
        // B contributes 0π electrons (empty p orbital); N contributes 2π (lone pair).
        make_ring(&[
            "B#h#v2#a0", "N#h#v2#a2", "B#h#v2#a0", "N#h#v2#a2", "B#h#v2#a0", "N#h#v2#a2",
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
        #[case] builder: MoleculeBuilder,
        #[case] expected_atoms: usize,
        #[case] expected_electrons: u8,
    ) {
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate_builder(&builder);
        let model = daylight_model();
        let systems = model.find_from_rings(&builder, &rings);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].contributions().len(), expected_atoms);
        assert_eq!(systems[0].electron_count(), expected_electrons);
    }

    #[rstest]
    #[case::cyclobutadiene(cyclobutadiene(), daylight_model())]
    #[case::cyclohexane(cyclohexane(), daylight_model())]
    #[case::cubane(cubane(), daylight_model())]
    #[case::borazine_daylight(borazine(), daylight_model())]
    #[case::pyrrole_mdl(pyrrole(), mdl_model())]
    fn test_find_from_rings_non_aromatic(
        #[case] builder: MoleculeBuilder,
        #[case] model: HueckelRuleAromaticity,
    ) {
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate_builder(&builder);
        let systems = model.find_from_rings(&builder, &rings);
        assert!(systems.is_empty());
    }

    #[rstest]
    fn test_borazine_permissive_aromatic(borazine: MoleculeBuilder) {
        // B contributes 0π, N contributes 2π: 6 total → Hückel 4(1)+2 aromatic.
        // The daylight model excludes B, so only the permissive model finds this system.
        let rings = RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
            .enumerate_builder(&borazine);
        let systems = permissive_model().find_from_rings(&borazine, &rings);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].contributions().len(), 6);
        assert_eq!(systems[0].electron_count(), 6);
    }
}
