//! Substructure matcher over MoleculeAst using VF2 subgraph isomorphism.

use petgraph::algo::subgraph_isomorphisms_iter;
use petgraph::graph::Graph;
use petgraph::Directed;

use index_vec::Idx;

use crate::ast::molecule::MoleculeAst;
use crate::solver::Solver;

/// Query atom index → target atom index mapping.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assignment(pub Vec<usize>);

/// Precomputed graph view of a target molecule for substructure matching.
pub struct MatchTarget<'a> {
    ast: &'a MoleculeAst,
    graph: Graph<usize, usize, Directed>,
}

/// Precomputed graph view of a query pattern for substructure matching.
pub struct MatchQuery<'a> {
    ast: &'a MoleculeAst,
    graph: Graph<usize, usize, Directed>,
}

fn build_graph(ast: &MoleculeAst) -> Graph<usize, usize, Directed> {
    let mut graph = Graph::with_capacity(ast.atom_count(), ast.bond_count() * 2);
    let nodes: Vec<_> = (0..ast.atom_count()).map(|i| graph.add_node(i)).collect();
    for (bond_idx, src, tgt, _) in ast.bonds() {
        let s = src.index();
        let t = tgt.index();
        graph.add_edge(nodes[s], nodes[t], bond_idx.index());
        graph.add_edge(nodes[t], nodes[s], bond_idx.index());
    }
    graph
}

impl<'a> MatchTarget<'a> {
    pub fn new(ast: &'a MoleculeAst) -> Self {
        Self {
            graph: build_graph(ast),
            ast,
        }
    }

    pub fn ast(&self) -> &MoleculeAst {
        self.ast
    }
}

impl<'a> MatchQuery<'a> {
    pub fn new(ast: &'a MoleculeAst) -> Self {
        Self {
            graph: build_graph(ast),
            ast,
        }
    }

    pub fn ast(&self) -> &MoleculeAst {
        self.ast
    }
}

/// Find all substructure matches of `query` in `target`.
pub fn find_matches(query: &MatchQuery, target: &MatchTarget) -> Vec<Assignment> {
    let q_ast = query.ast;
    let t_ast = target.ast;

    if q_ast.atom_count() == 0 {
        return if post_filter(q_ast, t_ast, &[]) {
            vec![Assignment(vec![])]
        } else {
            vec![]
        };
    }

    let mut node_match = |q_idx: &usize, t_idx: &usize| {
        q_ast
            .atom((*q_idx).into())
            .matches_ground(t_ast.atom((*t_idx).into()))
    };
    let mut edge_match = |q_idx: &usize, t_idx: &usize| {
        q_ast
            .bond((*q_idx).into())
            .matches_ground(t_ast.bond((*t_idx).into()))
    };

    let q_ref = &query.graph;
    let t_ref = &target.graph;
    let Some(iter) = subgraph_isomorphisms_iter(
        &q_ref,
        &t_ref,
        &mut node_match,
        &mut edge_match,
    ) else {
        return vec![];
    };

    iter.filter(|assignment| post_filter(q_ast, t_ast, assignment))
        .map(Assignment)
        .collect()
}

/// Find all substructure matches of `query` in `target`, then post-filter
/// through the solver's valence validation.
pub fn find_matches_with(query: &MatchQuery, target: &MatchTarget, solver: &Solver) -> Vec<Assignment> {
    let assignments = find_matches(query, target);
    solver.filter(query.ast, target.ast, assignments)
}

fn post_filter(query: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    check_dative_bonds(query, target, assignment)
        && check_noncovalent_bonds(query, target, assignment)
        && check_aromatic_systems(query, target, assignment)
        && check_multicenter_bonds(query, target, assignment)
}

fn check_dative_bonds(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.dative_bond_ids().all(|qid| {
        let q_parts = query.dative_bond_participants(qid);
        let q_data = query.dative_bond(qid);
        let mapped_source = assignment[q_parts[0].index()];
        let mapped_target = assignment[q_parts[1].index()];
        target.dative_bond_ids().any(|tid| {
            let t_parts = target.dative_bond_participants(tid);
            t_parts[0].index() == mapped_source
                && t_parts[1].index() == mapped_target
                && q_data.matches_ground(target.dative_bond(tid))
        })
    })
}

fn check_noncovalent_bonds(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.noncovalent_bond_ids().all(|qid| {
        let q_parts = query.noncovalent_bond_participants(qid);
        let q_data = query.noncovalent_bond(qid);
        let mapped_source = assignment[q_parts[0].index()];
        let mapped_target = assignment[q_parts[1].index()];
        target.noncovalent_bond_ids().any(|tid| {
            let t_parts = target.noncovalent_bond_participants(tid);
            t_parts[0].index() == mapped_source
                && t_parts[1].index() == mapped_target
                && q_data.matches_ground(target.noncovalent_bond(tid))
        })
    })
}

fn check_aromatic_systems(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.aromatic_system_ids().all(|qid| {
        let q_atoms = query.aromatic_system_participants(qid);
        let mapped: Vec<usize> = q_atoms.iter().map(|a| assignment[a.index()]).collect();
        target.aromatic_system_ids().any(|tid| {
            let t_atoms = target.aromatic_system_participants(tid);
            mapped.iter().all(|m| t_atoms.iter().any(|a| a.index() == *m))
        })
    })
}

fn check_multicenter_bonds(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.multicenter_bond_ids().all(|qid| {
        let q_atoms = query.multicenter_bond_participants(qid);
        let mapped: Vec<usize> = q_atoms.iter().map(|a| assignment[a.index()]).collect();
        target.multicenter_bond_ids().any(|tid| {
            let t_atoms = target.multicenter_bond_participants(tid);
            mapped.iter().all(|m| t_atoms.iter().any(|a| a.index() == *m))
        })
    })
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::atom_ast::ElementAst;
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::{AromaticSystemAst, MulticenterBondAst};
    use crate::ast::AtomIdx;

    fn mol_with_atoms(elements: &[Element]) -> MoleculeAst {
        let atoms = elements.iter().map(|&e| AtomAst::from_element(e)).collect();
        MoleculeAst::new(atoms, vec![], vec![], vec![], vec![], vec![], vec![])
    }

    fn find(query: &MoleculeAst, target: &MoleculeAst) -> Vec<Assignment> {
        find_matches(&MatchQuery::new(query), &MatchTarget::new(target))
    }

    fn sorted(mut assignments: Vec<Assignment>) -> Vec<Assignment> {
        assignments.sort_by(|a, b| a.0.cmp(&b.0));
        assignments
    }

    #[test]
    fn test_find_matches_single_atom() {
        let mol = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&mol, &mol), vec![Assignment(vec![0])]);
    }

    #[test]
    fn test_find_matches_chain_identity() {
        let mol = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(find(&mol, &mol), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_substructure() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            ],
            vec![], vec![], vec![], vec![], vec![],
        );

        let results = sorted(find(&query, &target));
        assert_eq!(
            results,
            vec![
                Assignment(vec![0, 1]),
                Assignment(vec![1, 0]),
                Assignment(vec![1, 2]),
                Assignment(vec![2, 1]),
            ]
        );
    }

    #[test]
    fn test_find_matches_element_mismatch() {
        let target = mol_with_atoms(&[Element::N]);
        let query = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_bond_order_mismatch() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![], vec![], vec![], vec![], vec![],
        );
        let target = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_wildcard_element() {
        let query = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let target = mol_with_atoms(&[Element::C, Element::N]);
        let results = sorted(find(&query, &target));
        assert_eq!(
            results,
            vec![Assignment(vec![0]), Assignment(vec![1])]
        );
    }

    #[test]
    fn test_find_matches_wildcard_bond_order() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
            vec![], vec![], vec![], vec![], vec![],
        );
        let target = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(find(&query, &target), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_optional_field_unconstrained() {
        let query = mol_with_atoms(&[Element::C]);
        let target = MoleculeAst::new(
            vec![AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::from_element(Element::C) }],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(find(&query, &target), vec![Assignment(vec![0])]);
    }

    #[test]
    fn test_find_matches_optional_field_present_vs_absent() {
        let query = MoleculeAst::new(
            vec![AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::from_element(Element::C) }],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let target = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_dative_bond_direction() {
        let atoms = vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::B)];
        let query = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![],
        );
        let target_correct = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![],
        );
        let target_reversed = MoleculeAst::new(
            atoms,
            vec![],
            vec![(AtomIdx(1), AtomIdx(0), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![],
        );

        assert_eq!(find(&query, &target_correct), vec![Assignment(vec![0, 1])]);
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_noncovalent_bond_direction() {
        let atoms = vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::O)];
        let query = MoleculeAst::new(
            atoms.clone(),
            vec![], vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![],
        );
        let target_correct = MoleculeAst::new(
            atoms.clone(),
            vec![], vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![],
        );
        let target_reversed = MoleculeAst::new(
            atoms,
            vec![], vec![],
            vec![(AtomIdx(1), AtomIdx(0), BondAst::from_order(1))],
            vec![], vec![], vec![],
        );

        assert_eq!(find(&query, &target_correct), vec![Assignment(vec![0, 1])]);
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_aromatic_system_subset() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst {})],
            vec![], vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
                (AtomIdx(2), AtomIdx(0), BondAst::from_order(1)),
            ],
            vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], AromaticSystemAst {})],
            vec![], vec![],
        );

        assert!(!find(&query, &target).is_empty());
    }

    #[test]
    fn test_find_matches_aromatic_system_mismatch() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst {})],
            vec![], vec![],
        );
        let target = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        );

        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_multicenter_subset() {
        let query = MoleculeAst::new(
            vec![AtomAst::from_element(Element::B), AtomAst::from_element(Element::H)],
            vec![], vec![], vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], MulticenterBondAst {})],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::B),
            ],
            vec![], vec![], vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], MulticenterBondAst {})],
            vec![],
        );

        assert!(!find(&query, &target).is_empty());
    }

    #[test]
    fn test_find_matches_multicenter_identity() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::B),
            ],
            vec![], vec![], vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], MulticenterBondAst {})],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::H),
            ],
            vec![], vec![], vec![], vec![],
            vec![
                (vec![AtomIdx(0), AtomIdx(2), AtomIdx(1)], MulticenterBondAst {}),
                (vec![AtomIdx(0), AtomIdx(3), AtomIdx(1)], MulticenterBondAst {}),
            ],
            vec![],
        );

        let results = find(&query, &target);
        assert!(!results.is_empty());
        for a in &results {
            let mapped: Vec<usize> = vec![a.0[0], a.0[1], a.0[2]];
            let in_single = target.multicenter_bond_ids().any(|tid| {
                let t_atoms = target.multicenter_bond_participants(tid);
                mapped.iter().all(|m| t_atoms.iter().any(|a| a.index() == *m))
            });
            assert!(in_single, "assignment {a:?} spans multiple multicenter bonds");
        }
    }

    #[test]
    fn test_find_matches_empty_query() {
        let query = MoleculeAst::default();
        let target = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(find(&query, &target), vec![Assignment(vec![])]);
    }

    #[test]
    fn test_find_matches_no_match() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            ],
            vec![], vec![], vec![], vec![], vec![],
        );
        let target = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_with_identity() {
        let query = mol_with_atoms(&[Element::C]);
        let target = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![], vec![], vec![], vec![], vec![],
        );

        let solver = Solver::default();
        let q = MatchQuery::new(&query);
        let t = MatchTarget::new(&target);
        assert_eq!(find_matches_with(&q, &t, &solver), find_matches(&q, &t));
    }
}
