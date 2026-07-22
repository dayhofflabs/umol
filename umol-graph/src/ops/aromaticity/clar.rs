//! Clar aromaticity model (π-sextet optimization).
//!
//! Maximizes the number of disjoint aromatic π-sextets across a fused-ring
//! system using branch-and-bound on a sextet-conflict graph. Applies only to
//! benzenoid hydrocarbons (all-carbon, all 6-membered rings); rejects
//! non-benzenoid or heterocyclic input.

use std::collections::HashSet;

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemAst, AtomId, AtomView, ElementAst, MoleculeAst, RingId, RingSet, SpinStateAst,
};
use umol_chem::element::Element;
use umol_graph_core::{Graph, MaximumIndependentSetAlgorithm};

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
    ) -> Result<Vec<(Vec<AtomId>, AromaticSystemAst)>, ClarError>
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let has_non_benzenoid = ast.atoms().iter().any(|view| {
            !matches!(view.ast.element, ElementAst::Lit(Element::C))
                && electrons_at(&view).is_some()
        });
        if has_non_benzenoid {
            return Err(ClarError::NonBenzenoid(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present"
                    .to_string(),
            ));
        }

        let sextet_indices: Vec<RingId> = rings
            .ids()
            .filter(|&i| {
                let Some(cycle) = rings.get(i) else {
                    return false;
                };
                cycle.len() == 6
                    && cycle.atoms().iter().all(|&atom| {
                        let a = ast.atom(atom);
                        matches!(a.ast.element, ElementAst::Lit(Element::C))
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

        let selected_atoms: HashSet<AtomId> = best_sextet_indices
            .iter()
            .filter_map(|&i| rings.get(i))
            .flat_map(|r| r.atoms().iter().copied())
            .collect();

        let mut atoms: Vec<AtomId> = selected_atoms.into_iter().collect();
        atoms.sort_unstable();

        let electrons: Vec<i64> = atoms
            .iter()
            .map(|&atom| electrons_at(&ast.atom(atom)).unwrap_or(0) as i64)
            .collect();

        Ok(vec![(
            atoms,
            AromaticSystemAst::from_electrons(electrons)
                .with_charge(0)
                .with_spin(SpinStateAst::closed_shell()),
        )])
    }
}

fn select_disjoint_sextets(rings: &RingSet, candidates: &[RingId]) -> Vec<RingId> {
    if candidates.is_empty() {
        return Vec::new();
    }

    let candidate_atoms: Vec<HashSet<AtomId>> = candidates
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
    let selected =
        conflict_graph.maximum_independent_set(MaximumIndependentSetAlgorithm::BranchAndBound);
    selected
        .into_iter()
        .map(|node_id| candidates[node_id.index()])
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticValenceAst, AtomAst, AtomConstraintAst, AtomId, BondAst, ElementAst, MoleculeAst,
        MoleculeParts, RingFamily, RingId, ValueAst,
    };
    use umol_chem::element::Element;

    use super::*;

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
                    atom.constraints.set(AtomConstraintAst::AromaticValence(
                        AromaticValenceAst::Aromatic(ValueAst::Lit(n)),
                    ));
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
                    AtomId(i as u32),
                    AtomId(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn make_fused(specs: Vec<(AtomAst, Option<i64>)>, edges: &[(usize, usize)]) -> MoleculeAst {
        let atoms = apply_pi(specs);
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomId(a as u32), AtomId(b as u32), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn enumerate_induced(ast: &MoleculeAst) -> RingSet {
        ast.rings_with(RingFamily::Simple, 6, |_| true)
            .into_ring_set()
    }

    fn hex_ring_indices(ast: &MoleculeAst, ring_info: &RingSet) -> Vec<RingId> {
        ring_info
            .ids()
            .filter(|&i| {
                ring_info.get(i).is_some_and(|cycle| {
                    cycle.len() == 6
                        && cycle.atoms().iter().all(|&atom| {
                            matches!(ast.atom(atom).ast.element, ElementAst::Lit(Element::C))
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

    #[rstest]
    #[case::benzene(make_ring(vec![aromatic(Element::C, 1); 6]), 1)]
    #[case::naphthalene(naphthalene(), 1)]
    #[case::phenanthrene(phenanthrene(), 2)]
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
    #[case::cyclohexane(make_ring(vec![plain(Element::C); 6]), 0, None)]
    fn test_clar_aromaticity_find_from_rings(
        #[case] ast: MoleculeAst,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let rings = enumerate_induced(&ast);
        let model = ClarAromaticity;
        let systems = model
            .find_from_rings(&ast, &rings, &|v| match v
                .ast
                .constraints
                .aromatic_valence()
                .unwrap_or(&AromaticValenceAst::Undetermined)
            {
                AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                _ => None,
            })
            .unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(systems.first().map(|s| s.0.len()), expected_atoms);
        if let Some((system_atoms_vec, _)) = systems.first() {
            let system_atoms: HashSet<AtomId> = system_atoms_vec.iter().copied().collect();
            let expected_atoms: HashSet<AtomId> =
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
        assert!(model
            .find_from_rings(&ast, &rings, &|v| {
                match v
                    .ast
                    .constraints
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined)
                {
                    AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                    _ => None,
                }
            })
            .is_err());
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
