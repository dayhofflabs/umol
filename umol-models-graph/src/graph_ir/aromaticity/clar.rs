//! Clar aromaticity model (pi-sextet optimization).
//!
//! Maximizes the number of disjoint aromatic pi-sextets across a fused ring
//! system using a branch-and-bound solver over candidate rings. Applicable
//! only to benzenoid hydrocarbons (all-carbon, all 6-membered rings). Returns
//! an error for non-benzenoid or heterocyclic systems.

use std::collections::HashSet;

use umol_data::Element;

use super::{AromaticContribution, AromaticSystem, AromaticityError};
use crate::algorithms::maximum_independent_set;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};
use crate::graph_ir::rings::{Ring, RingIndex, RingSet};

#[derive(Clone, Debug)]
pub struct ClarAromaticity;

impl ClarAromaticity {
    pub fn find_from_rings(
        &self,
        builder: &MoleculeBuilder,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let has_non_benzenoid_aromatic = builder.atom_indices().any(|atom| {
            builder.atom(atom).is_some_and(|a| {
                builder.atom_has_aromatic_candidate(atom) && a.element() != Element::C
            })
        });
        if has_non_benzenoid_aromatic {
            return Err(AromaticityError::ClarInputError(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present"
                    .to_string(),
            ));
        }

        let sextet_indices: Vec<RingIndex> = rings
            .ring_indices()
            .filter(|&i| {
                let Some(cycle) = rings.ring(i) else {
                    return false;
                };
                cycle.len() == 6
                    && cycle.atoms().iter().all(|&atom| {
                        builder
                            .atom(atom)
                            .map(|a| {
                                a.element() == Element::C
                                    && builder.atom_has_aromatic_candidate(atom)
                            })
                            .unwrap_or(false)
                    })
            })
            .collect();

        if sextet_indices.is_empty() {
            return Ok(Vec::new());
        }

        let best_sextet_indices = select_disjoint_sextets(rings, &sextet_indices);
        if best_sextet_indices.is_empty() {
            return Ok(Vec::new());
        }

        let selected_atoms: HashSet<AtomIndex> = best_sextet_indices
            .iter()
            .filter_map(|&i| rings.ring(i))
            .flat_map(|r| r.atoms().iter().copied())
            .collect();

        let contributions: Vec<AromaticContribution> = selected_atoms
            .iter()
            .map(|&atom| AromaticContribution::new(atom, builder.atom_aromatic_valence(atom)))
            .collect();

        let selected_rings: Vec<Ring> = best_sextet_indices
            .iter()
            .filter_map(|&ring_idx| rings.ring(ring_idx).cloned())
            .collect();

        Ok(vec![AromaticSystem::with_rings(
            contributions,
            selected_rings,
        )])
    }
}

fn select_disjoint_sextets(rings: &RingSet, candidates: &[RingIndex]) -> Vec<RingIndex> {
    if candidates.is_empty() {
        return Vec::new();
    }

    // Domain adapter invariant: candidates are mapped to contiguous integer
    // ids in stable input order, then mapped back after MIS selection.
    let candidate_atoms: Vec<HashSet<AtomIndex>> = candidates
        .iter()
        .map(|&ring_idx| {
            rings
                .ring(ring_idx)
                .map(|ring| ring.atoms().iter().copied().collect())
                .unwrap_or_default()
        })
        .collect();

    let n = candidates.len();
    let mut conflict_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if !candidate_atoms[i].is_disjoint(&candidate_atoms[j]) {
                conflict_adj[i].push(j);
                conflict_adj[j].push(i);
            }
        }
    }

    let selected = maximum_independent_set(&conflict_adj);
    selected.into_iter().map(|i| candidates[i]).collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::config::RingEnumerationStrategy;
    use crate::graph_ir::rings::{RingEnumerator, RingFamily};
    const C1: &str = "C#h#v2#a";
    const C0: &str = "C#v4";

    fn make_ring(atom_specs: &[&str]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_resolved_atom(s.parse().unwrap()))
            .collect();
        let n = atoms.len();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondBuilder::new(1, None));
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
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        builder
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
    fn phenanthrene() -> MoleculeBuilder {
        make_fused(
            &[C1; 14],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
                (8, 10), (10, 11), (11, 12), (12, 13), (13, 9),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn coronene() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..24)
            .map(|_| builder.add_resolved_atom(C1.parse().unwrap()))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        for i in 0..6 {
            let a = i;
            let b = (i + 1) % 6;
            let c = 6 + i * 3;
            let d = 6 + i * 3 + 1;
            let e = 6 + i * 3 + 2;
            builder.add_bond_unchecked(atoms[a], atoms[c], BondBuilder::new(1, None));
            builder.add_bond_unchecked(atoms[c], atoms[d], BondBuilder::new(1, None));
            builder.add_bond_unchecked(atoms[d], atoms[e], BondBuilder::new(1, None));
            builder.add_bond_unchecked(atoms[e], atoms[b], BondBuilder::new(1, None));
        }
        for i in 0..6 {
            let this_last = 6 + i * 3 + 2;
            let next_first = 6 + ((i + 1) % 6) * 3;
            builder.add_bond_unchecked(
                atoms[this_last],
                atoms[next_first],
                BondBuilder::new(1, None),
            );
        }
        builder
    }

    fn hex_ring_indices(builder: &MoleculeBuilder, ring_info: &RingSet) -> Vec<RingIndex> {
        ring_info
            .ring_indices()
            .filter(|&i| {
                ring_info.ring(i).is_some_and(|cycle| {
                    cycle.len() == 6
                        && cycle.atoms().iter().all(|&atom| {
                            builder
                                .atom(atom)
                                .map(|a| a.element() == Element::C)
                                .unwrap_or(false)
                        })
                })
            })
            .collect()
    }

    #[rstest]
    #[case::benzene(make_ring(&[C1; 6]), 1)]
    #[case::naphthalene(naphthalene(), 1)]
    #[case::phenanthrene(phenanthrene(), 2)]
    #[case::coronene(coronene(), 3)]
    fn test_clar_model_sextet_count(
        #[case] builder: MoleculeBuilder,
        #[case] expected_sextets: usize,
    ) {
        let ring_info = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate_builder(&builder);
        let candidates = hex_ring_indices(&builder, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), expected_sextets);
    }

    #[rstest]
    #[case::benzene(make_ring(&[C1; 6]), 1, Some(6))]
    #[case::naphthalene(naphthalene(), 1, Some(6))]
    #[case::phenanthrene(phenanthrene(), 1, Some(12))]
    #[case::coronene(coronene(), 1, Some(18))]
    #[case::cyclohexane(make_ring(&[C0; 6]), 0, None)]
    fn test_clar_model_find_from_rings(
        #[case] builder: MoleculeBuilder,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let rings = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate_builder(&builder);
        let model = ClarAromaticity;
        let systems = model.find_from_rings(&builder, &rings).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(
            systems.get(0).map(|s| s.contributions().len()),
            expected_atoms
        );
        if let Some(system) = systems.first() {
            let ring_set: HashSet<RingIndex> = system
                .rings()
                .iter()
                .filter_map(|ring| {
                    rings
                        .ring_indices()
                        .find(|&idx| rings.ring(idx) == Some(ring))
                })
                .collect();
            let expected_selected: HashSet<RingIndex> =
                select_disjoint_sextets(&rings, &hex_ring_indices(&builder, &rings))
                    .into_iter()
                    .collect();
            assert_eq!(ring_set, expected_selected);
        }
    }

    #[rstest]
    #[case::pyrrole(make_ring(&["N#h#v2#a2", C1, C1, C1, C1]))]
    #[case::pyridine(make_ring(&["N#n#v2#a", C1, C1, C1, C1, C1]))]
    #[case::furan(make_ring(&["O#n#v2#a2", C1, C1, C1, C1]))]
    fn test_clar_model_find_from_rings_error(#[case] builder: MoleculeBuilder) {
        let rings = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate_builder(&builder);
        let model = ClarAromaticity;
        assert!(model.find_from_rings(&builder, &rings).is_err());
    }

    #[rstest]
    fn test_clar_solver(phenanthrene: MoleculeBuilder) {
        let ring_info = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate_builder(&phenanthrene);
        let candidates = hex_ring_indices(&phenanthrene, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), 2);
        for &sextet_idx in &sextets {
            let ring = ring_info.ring(sextet_idx).unwrap();
            assert_eq!(ring.len(), 6);
        }
    }
}
