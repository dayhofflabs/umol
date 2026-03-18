//! Clar aromaticity model (pi-sextet optimization).
//!
//! Maximizes the number of disjoint aromatic pi-sextets across a fused ring
//! system using a branch-and-bound solver over candidate rings. Applicable
//! only to benzenoid hydrocarbons (all-carbon, all 6-membered rings). Returns
//! an error for non-benzenoid or heterocyclic systems.

use std::collections::HashSet;

use umol_data::Element;

use super::{AromaticContribution, AromaticSystem};
use crate::atom::AromaticValence;
use crate::graph_ir::error::ResolutionError;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};
use crate::graph_ir::rings::{MoleculeRings, RingIndex};

#[derive(Clone, Debug)]
pub struct ClarAromaticity;

impl ClarAromaticity {
    pub fn find_from_rings(
        &self,
        builder: &MoleculeBuilder,
        rings: &MoleculeRings,
    ) -> Result<Vec<AromaticSystem>, ResolutionError> {
        let has_non_benzenoid_aromatic = rings.atom_rings.keys().any(|&atom| {
            builder.atom(atom).is_some_and(|a| {
                a.candidates()
                    .iter()
                    .any(|c| c.aromatic_valence().is_aromatic())
                    && a.element() != Element::C
            })
        });
        if has_non_benzenoid_aromatic {
            return Err(ResolutionError::AromaticityInconsistent(
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
                    && cycle.iter().all(|&atom| {
                        builder
                            .atom(atom)
                            .map(|a| {
                                a.element() == Element::C
                                    && a.candidates()
                                        .iter()
                                        .any(|c| c.aromatic_valence().is_aromatic())
                            })
                            .unwrap_or(false)
                    })
            })
            .collect();

        if sextet_indices.is_empty() {
            return Ok(Vec::new());
        }

        let all_atoms: HashSet<AtomIndex> = sextet_indices
            .iter()
            .filter_map(|&i| rings.ring(i))
            .flat_map(|r| r.iter().copied())
            .collect();

        let mut solver = ClarSolver::new(rings, &sextet_indices);
        let best_sextet_indices = solver.solve();

        if best_sextet_indices.is_empty() {
            return Ok(Vec::new());
        }

        let contributions: Vec<AromaticContribution> = all_atoms
            .iter()
            .map(|&atom| {
                let e = builder
                    .atom(atom)
                    .and_then(|a| {
                        a.candidates()
                            .iter()
                            .find_map(|c| match c.aromatic_valence() {
                                AromaticValence::Valence(e) => Some(e),
                                AromaticValence::None => None,
                            })
                    })
                    .unwrap_or(1);
                AromaticContribution::new(atom, e)
            })
            .collect();

        let all_rings: Vec<Vec<AtomIndex>> = sextet_indices
            .iter()
            .filter_map(|&ring_idx| rings.ring(ring_idx).map(|r| r.to_vec()))
            .collect();

        Ok(vec![AromaticSystem::with_rings(contributions, all_rings)])
    }
}

pub(crate) struct ClarSolver<'a> {
    rings: &'a MoleculeRings,
    candidates: &'a [RingIndex],
    current: Vec<RingIndex>,
    best: Vec<RingIndex>,
    used_atoms: HashSet<AtomIndex>,
}

impl<'a> ClarSolver<'a> {
    pub(crate) fn new(rings: &'a MoleculeRings, candidates: &'a [RingIndex]) -> Self {
        Self {
            rings,
            candidates,
            current: Vec::new(),
            best: Vec::new(),
            used_atoms: HashSet::new(),
        }
    }

    pub(crate) fn solve(&mut self) -> Vec<RingIndex> {
        self.branch(0);
        self.best.clone()
    }

    fn branch(&mut self, pos: usize) {
        let remaining = self.candidates.len().saturating_sub(pos);
        if self.current.len() + remaining <= self.best.len() {
            return;
        }
        if pos == self.candidates.len() {
            if self.current.len() > self.best.len() {
                self.best = self.current.clone();
            }
            return;
        }

        let ring = self.candidates[pos];
        if self.can_add(ring) {
            let added_atoms = self.add_ring(ring);
            self.branch(pos + 1);
            self.remove_ring(ring, &added_atoms);
        }
        self.branch(pos + 1);
    }

    fn can_add(&self, ring: RingIndex) -> bool {
        self.rings
            .ring(ring)
            .is_some_and(|r| r.iter().all(|a| !self.used_atoms.contains(a)))
    }

    fn add_ring(&mut self, ring: RingIndex) -> Vec<AtomIndex> {
        self.current.push(ring);
        let atoms: Vec<AtomIndex> = self
            .rings
            .ring(ring)
            .map(|r| r.to_vec())
            .unwrap_or_default();
        self.used_atoms.extend(atoms.iter().copied());
        atoms
    }

    fn remove_ring(&mut self, _ring: RingIndex, atoms: &[AtomIndex]) {
        self.current.pop();
        for &atom in atoms {
            self.used_atoms.remove(&atom);
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::atom;
    use crate::graph_ir::bond::BondBuilder;
    use crate::graph_ir::config::RingEnumerationStrategy;
    use crate::graph_ir::rings::RingEnumerator;

    const C1: &str = "{Cv2a1H}";
    const C0: &str = "{Cv4}";

    fn make_ring(atom_specs: &[&str]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_atom(atom!(s)))
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
            .map(|s| builder.add_atom(atom!(s)))
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
            .map(|_| builder.add_atom(atom!(C1)))
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

    fn hex_ring_indices(builder: &MoleculeBuilder, ring_info: &MoleculeRings) -> Vec<RingIndex> {
        ring_info
            .ring_indices()
            .filter(|&i| {
                ring_info.ring(i).is_some_and(|cycle| {
                    cycle.len() == 6
                        && cycle.iter().all(|&atom| {
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
        let ring_info =
            RingEnumerator::new(&RingEnumerationStrategy::default()).enumerate_builder(&builder);
        let candidates = hex_ring_indices(&builder, &ring_info);
        let mut solver = ClarSolver::new(&ring_info, &candidates);
        let sextets = solver.solve();
        assert_eq!(sextets.len(), expected_sextets);
    }

    #[rstest]
    #[case::benzene(make_ring(&[C1; 6]), 1, Some(6))]
    #[case::naphthalene(naphthalene(), 1, Some(10))]
    #[case::phenanthrene(phenanthrene(), 1, Some(14))]
    #[case::coronene(coronene(), 1, Some(24))]
    #[case::cyclohexane(make_ring(&[C0; 6]), 0, None)]
    fn test_clar_model_find_from_rings(
        #[case] builder: MoleculeBuilder,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let rings =
            RingEnumerator::new(&RingEnumerationStrategy::default()).enumerate_builder(&builder);
        let model = ClarAromaticity;
        let systems = model.find_from_rings(&builder, &rings).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(
            systems.get(0).map(|s| s.contributions().len()),
            expected_atoms
        );
    }

    #[rstest]
    #[case::pyrrole(make_ring(&["{Nv2a2H}", C1, C1, C1, C1]))]
    #[case::pyridine(make_ring(&["{N/1v2a1}", C1, C1, C1, C1, C1]))]
    #[case::furan(make_ring(&["{O/1v2a2}", C1, C1, C1, C1]))]
    fn test_clar_model_find_from_rings_error(#[case] builder: MoleculeBuilder) {
        let rings =
            RingEnumerator::new(&RingEnumerationStrategy::default()).enumerate_builder(&builder);
        let model = ClarAromaticity;
        assert!(model.find_from_rings(&builder, &rings).is_err());
    }

    #[rstest]
    fn test_clar_solver(phenanthrene: MoleculeBuilder) {
        let ring_info = RingEnumerator::new(&RingEnumerationStrategy::default())
            .enumerate_builder(&phenanthrene);
        let candidates = hex_ring_indices(&phenanthrene, &ring_info);
        let mut solver = ClarSolver::new(&ring_info, &candidates);
        let sextets = solver.solve();
        assert_eq!(sextets.len(), 2);
        for &sextet_idx in &sextets {
            let ring = ring_info.ring(sextet_idx).unwrap();
            assert_eq!(ring.len(), 6);
        }
    }
}
