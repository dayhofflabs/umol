//! Clar aromaticity model (pi-sextet optimization).
//!
//! Maximizes the number of disjoint aromatic pi-sextets across a fused ring
//! system using a branch-and-bound solver over candidate rings. Applicable
//! only to benzenoid hydrocarbons (all-carbon, all 6-membered rings). Returns
//! an error for non-benzenoid or heterocyclic systems.

use std::collections::HashSet;

use umol_shared::element::Element;

use umol_shared::atom_ast::{AromaticValenceAst, ElementAst};
use umol_shared::value_ast::ValueAst;

use umol_graph_core::Graph;

use super::{AromaticContribution, AromaticSystem, AromaticityError};
use crate::ast::AtomIdx;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{Ring, RingIndex, RingSet};

#[derive(Clone, Debug)]
pub struct ClarAromaticity;

impl ClarAromaticity {
    pub fn find_from_rings(
        &self,
        ast: &MoleculeAst,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let has_non_benzenoid = ast.atoms().any(|(_, atom)| {
            !matches!(atom.element, ElementAst::Lit(Element::C))
                && matches!(atom.aromatic_valence, AromaticValenceAst::Value(ValueAst::Lit(_)))
        });
        if has_non_benzenoid {
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
                        let a = ast.atom(atom);
                        matches!(a.element, ElementAst::Lit(Element::C))
                            && matches!(
                                a.aromatic_valence,
                                AromaticValenceAst::Value(ValueAst::Lit(_))
                            )
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

        let selected_atoms: HashSet<AtomIdx> = best_sextet_indices
            .iter()
            .filter_map(|&i| rings.ring(i))
            .flat_map(|r| r.atoms().iter().copied())
            .collect();

        let contributions: Vec<AromaticContribution> = selected_atoms
            .iter()
            .map(|&atom| {
                let a = ast.atom(atom);
                let valence = match a.aromatic_valence {
                    AromaticValenceAst::Value(ValueAst::Lit(n)) => n as u8,
                    _ => 0,
                };
                AromaticContribution::new(atom, valence)
            })
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

    let candidate_atoms: Vec<HashSet<AtomIdx>> = candidates
        .iter()
        .map(|&ring_idx| {
            rings
                .ring(ring_idx)
                .map(|ring| ring.atoms().iter().copied().collect())
                .unwrap_or_default()
        })
        .collect();

    let n = candidates.len();
    let mut edges: Vec<[u32; 2]> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if !candidate_atoms[i].is_disjoint(&candidate_atoms[j]) {
                edges.push([i as u32, j as u32]);
            }
        }
    }

    let conflict_graph = Graph::new(n, &edges);
    let selected = conflict_graph.maximum_independent_set();
    selected
        .into_iter()
        .map(|node_id| candidates[node_id.index()])
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst};
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::AtomIdx;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::rings::RingEnumerationStrategy;
    use crate::ast::rings::{RingEnumerator, RingFamily};

    fn aromatic_atom(element: Element, pi: i64) -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(element),
            aromatic_valence: AromaticValenceAst::Value(ValueAst::Lit(pi)),
            ..Default::default()
        }
    }

    fn make_ring(atoms: Vec<AtomAst>) -> MoleculeAst {
        let n = atoms.len();
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    fn make_fused(atoms: Vec<AtomAst>, edges: &[(usize, usize)]) -> MoleculeAst {
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    fn hex_ring_indices(ast: &MoleculeAst, ring_info: &RingSet) -> Vec<RingIndex> {
        ring_info
            .ring_indices()
            .filter(|&i| {
                ring_info.ring(i).is_some_and(|cycle| {
                    cycle.len() == 6
                        && cycle.atoms().iter().all(|&atom| {
                            matches!(ast.atom(atom).element, ElementAst::Lit(Element::C))
                        })
                })
            })
            .collect()
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeAst {
        make_fused(
            vec![aromatic_atom(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn phenanthrene() -> MoleculeAst {
        make_fused(
            vec![aromatic_atom(Element::C, 1); 14],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
                (8, 10), (10, 11), (11, 12), (12, 13), (13, 9),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn coronene() -> MoleculeAst {
        let atoms = vec![aromatic_atom(Element::C, 1); 24];
        let mut edges = Vec::new();
        for i in 0..6 {
            edges.push((i, (i + 1) % 6));
        }
        for i in 0..6 {
            let a = i;
            let b = (i + 1) % 6;
            let c = 6 + i * 3;
            let d = 6 + i * 3 + 1;
            let e = 6 + i * 3 + 2;
            edges.push((a, c));
            edges.push((c, d));
            edges.push((d, e));
            edges.push((e, b));
        }
        for i in 0..6 {
            let this_last = 6 + i * 3 + 2;
            let next_first = 6 + ((i + 1) % 6) * 3;
            edges.push((this_last, next_first));
        }
        make_fused(atoms, &edges)
    }

    #[rstest]
    #[case::benzene(make_ring(vec![aromatic_atom(Element::C, 1); 6]), 1)]
    #[case::naphthalene(naphthalene(), 1)]
    #[case::phenanthrene(phenanthrene(), 2)]
    #[case::coronene(coronene(), 3)]
    fn test_clar_aromaticity_sextet_count(
        #[case] ast: MoleculeAst,
        #[case] expected_sextets: usize,
    ) {
        let ring_info = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate(&ast);
        let candidates = hex_ring_indices(&ast, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), expected_sextets);
    }

    #[rstest]
    #[case::benzene(make_ring(vec![aromatic_atom(Element::C, 1); 6]), 1, Some(6))]
    #[case::naphthalene(naphthalene(), 1, Some(6))]
    #[case::phenanthrene(phenanthrene(), 1, Some(12))]
    #[case::coronene(coronene(), 1, Some(18))]
    #[case::cyclohexane(make_ring(vec![AtomAst::from_element(Element::C); 6]), 0, None)]
    fn test_clar_aromaticity_find_from_rings(
        #[case] ast: MoleculeAst,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let rings = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate(&ast);
        let model = ClarAromaticity;
        let systems = model.find_from_rings(&ast, &rings).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(
            systems.first().map(|s| s.contributions().len()),
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
                select_disjoint_sextets(&rings, &hex_ring_indices(&ast, &rings))
                    .into_iter()
                    .collect();
            assert_eq!(ring_set, expected_selected);
        }
    }

    #[rstest]
    #[case::pyrrole(make_ring(vec![
        aromatic_atom(Element::N, 2),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
    ]))]
    #[case::pyridine(make_ring(vec![
        aromatic_atom(Element::N, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
    ]))]
    #[case::furan(make_ring(vec![
        aromatic_atom(Element::O, 2),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
        aromatic_atom(Element::C, 1),
    ]))]
    fn test_clar_aromaticity_find_from_rings_error(#[case] ast: MoleculeAst) {
        let rings = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate(&ast);
        let model = ClarAromaticity;
        assert!(model.find_from_rings(&ast, &rings).is_err());
    }

    #[rstest]
    fn test_clar_aromaticity_solver(phenanthrene: MoleculeAst) {
        let ring_info = RingEnumerator::new(
            RingFamily::InducedBenzenoid,
            &RingEnumerationStrategy::default(),
        )
        .enumerate(&phenanthrene);
        let candidates = hex_ring_indices(&phenanthrene, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), 2);
        for &sextet_idx in &sextets {
            let ring = ring_info.ring(sextet_idx).unwrap();
            assert_eq!(ring.len(), 6);
        }
    }
}
