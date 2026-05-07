//! Clar aromaticity model (π-sextet optimization).
//!
//! Maximizes the number of disjoint aromatic π-sextets across a fused-ring
//! system using branch-and-bound on a sextet-conflict graph. Applies only to
//! benzenoid hydrocarbons (all-carbon, all 6-membered rings); rejects
//! non-benzenoid or heterocyclic input.

use std::collections::HashSet;

use umol_ast::ast::{
    AromaticSystemAst, AtomIdx, AtomView, ElementAst, MoleculeAst, RingIdx, RingSet, SpinStateAst,
    ValueAst,
};
use umol_graph_core::Graph;
use umol_shared::element::Element;

use thiserror::Error;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ClarError {
    #[error("clar: non-benzenoid input: {0}")]
    NonBenzenoid(String),
}

#[derive(Clone, Debug)]
pub struct ClarAromaticity;

impl ClarAromaticity {
    pub fn find_from_rings<F>(
        &self,
        ast: &MoleculeAst,
        rings: &RingSet,
        electrons_at: &F,
    ) -> Result<Vec<(Vec<AtomIdx>, AromaticSystemAst)>, ClarError>
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let has_non_benzenoid = ast.atoms().iter().any(|view| {
            !matches!(view.data.element, ElementAst::Lit(Element::C))
                && electrons_at(&view).is_some()
        });
        if has_non_benzenoid {
            return Err(ClarError::NonBenzenoid(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present"
                    .to_string(),
            ));
        }

        let sextet_indices: Vec<RingIdx> = rings
            .ids()
            .filter(|&i| {
                let Some(cycle) = rings.get(i) else {
                    return false;
                };
                cycle.len() == 6
                    && cycle.atoms().iter().all(|&atom| {
                        let a = ast.atom(atom);
                        matches!(a.data.element, ElementAst::Lit(Element::C))
                            && electrons_at(&a).is_some()
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
            .filter_map(|&i| rings.get(i))
            .flat_map(|r| r.atoms().iter().copied())
            .collect();

        let mut atoms: Vec<AtomIdx> = selected_atoms.into_iter().collect();
        atoms.sort_unstable();

        let electrons: Vec<ValueAst> = atoms
            .iter()
            .map(|&atom| {
                let pi = electrons_at(&ast.atom(atom)).unwrap_or(0);
                ValueAst::Lit(pi as i64)
            })
            .collect();

        Ok(vec![(
            atoms,
            AromaticSystemAst::new(electrons, ValueAst::Lit(0), SpinStateAst::closed_shell()),
        )])
    }
}

fn select_disjoint_sextets(rings: &RingSet, candidates: &[RingIdx]) -> Vec<RingIdx> {
    if candidates.is_empty() {
        return Vec::new();
    }

    let candidate_atoms: Vec<HashSet<AtomIdx>> = candidates
        .iter()
        .map(|&ring_idx| {
            rings
                .get(ring_idx)
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
    let selected = conflict_graph
        .maximum_independent_set(umol_graph_core::MaxIndependentSetAlgorithm::BranchAndBound);
    selected
        .into_iter()
        .map(|node_id| candidates[node_id.index()])
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst, Constraints, ElementAst,
        MoleculeAst, RingFamily, RingIdx, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;
    use crate::ops::aromaticity::electrons_from_aromatic_constraint;

    fn aromatic(element: Element, pi: i64) -> (AtomAst, Option<i64>) {
        (AtomAst::from_element(element), Some(pi))
    }

    fn plain(element: Element) -> (AtomAst, Option<i64>) {
        (AtomAst::from_element(element), None)
    }

    fn apply_pi(specs: Vec<(AtomAst, Option<i64>)>) -> Vec<AtomAst> {
        specs
            .into_iter()
            .map(|(mut atom, pi)| {
                if let Some(n) = pi {
                    atom.constraints
                        .add(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(
                            ValueAst::Lit(n),
                        )));
                }
                atom
            })
            .collect()
    }

    fn make_ring(specs: Vec<(AtomAst, Option<i64>)>) -> MoleculeAst {
        let n = specs.len();
        let atoms = apply_pi(specs);
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    fn make_fused(specs: Vec<(AtomAst, Option<i64>)>, edges: &[(usize, usize)]) -> MoleculeAst {
        let atoms = apply_pi(specs);
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    fn enumerate_induced(ast: &MoleculeAst) -> RingSet {
        ast.enumerate_rings(RingFamily::Simple, 6, |_| true)
    }

    fn hex_ring_indices(ast: &MoleculeAst, ring_info: &RingSet) -> Vec<RingIdx> {
        ring_info
            .ids()
            .filter(|&i| {
                ring_info.get(i).is_some_and(|cycle| {
                    cycle.len() == 6
                        && cycle.atoms().iter().all(|&atom| {
                            matches!(ast.atom(atom).data.element, ElementAst::Lit(Element::C))
                        })
                })
            })
            .collect()
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
    fn phenanthrene() -> MoleculeAst {
        make_fused(
            vec![aromatic(Element::C, 1); 14],
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
        // 24 atoms, 30 edges: inner hexagon (atoms 0..6), outer 18-cycle
        // (atoms 6..24), and 6 spokes from inner i to outer 6+3i. Real
        // coronene topology — 7 hexagonal faces in the planar embedding.
        let atoms: Vec<(AtomAst, Option<i64>)> =
            (0..24).map(|_| aromatic(Element::C, 1)).collect();
        let mut edges = Vec::new();
        for i in 0..6 {
            edges.push((i, (i + 1) % 6));
        }
        for i in 6..24 {
            edges.push((i, if i == 23 { 6 } else { i + 1 }));
        }
        for i in 0..6 {
            edges.push((i, 6 + 3 * i));
        }
        make_fused(atoms, &edges)
    }

    #[rstest]
    #[case::benzene(make_ring(vec![aromatic(Element::C, 1); 6]), 1)]
    #[case::naphthalene(naphthalene(), 1)]
    #[case::phenanthrene(phenanthrene(), 2)]
    #[case::coronene(coronene(), 3)]
    fn test_clar_aromaticity_sextet_count(
        #[case] ast: MoleculeAst,
        #[case] expected_sextets: usize,
    ) {
        let ring_info = enumerate_induced(&ast);
        let candidates = hex_ring_indices(&ast, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), expected_sextets);
    }

    #[rstest]
    #[case::benzene(make_ring(vec![aromatic(Element::C, 1); 6]), 1, Some(6))]
    #[case::naphthalene(naphthalene(), 1, Some(6))]
    #[case::phenanthrene(phenanthrene(), 1, Some(12))]
    #[case::coronene(coronene(), 1, Some(18))]
    #[case::cyclohexane(make_ring(vec![plain(Element::C); 6]), 0, None)]
    fn test_clar_aromaticity_find_from_rings(
        #[case] ast: MoleculeAst,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let rings = enumerate_induced(&ast);
        let model = ClarAromaticity;
        let systems = model.find_from_rings(&ast, &rings, &electrons_from_aromatic_constraint).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(systems.first().map(|s| s.0.len()), expected_atoms);
        if let Some((system_atoms_vec, _)) = systems.first() {
            let system_atoms: HashSet<AtomIdx> = system_atoms_vec.iter().copied().collect();
            let expected_atoms: HashSet<AtomIdx> =
                select_disjoint_sextets(&rings, &hex_ring_indices(&ast, &rings))
                    .into_iter()
                    .filter_map(|idx| rings.get(idx))
                    .flat_map(|r| r.atoms().iter().copied())
                    .collect();
            assert_eq!(system_atoms, expected_atoms);
        }
    }

    #[rstest]
    #[case::pyrrole(make_ring(vec![
        aromatic(Element::N, 2),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
    ]))]
    #[case::pyridine(make_ring(vec![
        aromatic(Element::N, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
    ]))]
    #[case::furan(make_ring(vec![
        aromatic(Element::O, 2),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
        aromatic(Element::C, 1),
    ]))]
    fn test_clar_aromaticity_find_from_rings_error(#[case] ast: MoleculeAst) {
        let rings = enumerate_induced(&ast);
        let model = ClarAromaticity;
        assert!(model.find_from_rings(&ast, &rings, &electrons_from_aromatic_constraint).is_err());
    }

    #[rstest]
    fn test_clar_aromaticity_solver(phenanthrene: MoleculeAst) {
        let ring_info = enumerate_induced(&phenanthrene);
        let candidates = hex_ring_indices(&phenanthrene, &ring_info);
        let sextets = select_disjoint_sextets(&ring_info, &candidates);
        assert_eq!(sextets.len(), 2);
        for &sextet_idx in &sextets {
            let ring = ring_info.get(sextet_idx).unwrap();
            assert_eq!(ring.len(), 6);
        }
    }
}
