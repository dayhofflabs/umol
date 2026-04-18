//! Substructure matcher over MoleculeAst using VF2 subgraph isomorphism.

use umol_graph_core::{subgraph_isomorphisms, EdgeId, NodeId};

use crate::ast::molecule::MoleculeAst;
use crate::unify::chemistry::Chemistry;
use crate::unify::propagate::ElectronInvariant;

/// Query atom index → target atom index mapping.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assignment(pub Vec<usize>);

/// Precomputed graph view of a target molecule for substructure matching.
pub struct MatchTarget<'a> {
    ast: &'a MoleculeAst,
}

/// Precomputed graph view of a query pattern for substructure matching.
pub struct MatchQuery<'a> {
    ast: &'a MoleculeAst,
}

impl<'a> MatchTarget<'a> {
    pub fn new(ast: &'a MoleculeAst) -> Self {
        Self { ast }
    }

    pub fn ast(&self) -> &MoleculeAst {
        self.ast
    }
}

impl<'a> MatchQuery<'a> {
    pub fn new(ast: &'a MoleculeAst) -> Self {
        Self { ast }
    }

    pub fn ast(&self) -> &MoleculeAst {
        self.ast
    }
}

/// Find all substructure matches of `query` in `target`.
pub fn find_matches(query: &MatchQuery, target: &MatchTarget) -> Vec<Assignment> {
    let q_ast = query.ast;
    let t_ast = target.ast;

    if q_ast.atoms().count() == 0 {
        return if post_filter(q_ast, t_ast, &[]) {
            vec![Assignment(vec![])]
        } else {
            vec![]
        };
    }

    let mut node_match = |q: NodeId, t: NodeId| {
        q_ast
            .atom(q.into())
            .data
            .matches_ground(t_ast.atom(t.into()).data)
    };
    let mut edge_match = |q: EdgeId, t: EdgeId| {
        q_ast
            .bond(q.into())
            .data
            .matches_ground(t_ast.bond(t.into()).data)
    };

    subgraph_isomorphisms(
        q_ast.graph(),
        t_ast.graph(),
        &mut node_match,
        &mut edge_match,
    )
    .into_iter()
    .filter(|assignment| post_filter(q_ast, t_ast, assignment))
    .map(Assignment)
    .collect()
}

/// Find all substructure matches of `query` in `target`, then post-filter
/// through the chemistry's valence validation.
pub fn find_matches_with(
    query: &MatchQuery,
    target: &MatchTarget,
    chemistry: &Chemistry,
) -> Vec<Assignment> {
    find_matches(query, target)
        .into_iter()
        .filter(|a| {
            a.0.iter().all(|&t_idx| {
                ElectronInvariant.validate(target.ast, t_idx)
                    && chemistry.valence.validate(target.ast, t_idx)
            })
        })
        .collect()
}

fn post_filter(query: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    check_dative_bonds(query, target, assignment)
        && check_noncovalent_bonds(query, target, assignment)
        && check_aromatic_systems(query, target, assignment)
        && check_multicenter_bonds(query, target, assignment)
}

fn check_dative_bonds(query: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    query.dative_bonds().iter().all(|q| {
        let mapped_donor = assignment[q.donor.index()];
        let mapped_acceptor = assignment[q.acceptor.index()];
        target.dative_bonds().iter().any(|t| {
            t.donor.index() == mapped_donor
                && t.acceptor.index() == mapped_acceptor
                && q.data.matches_ground(t.data)
        })
    })
}

fn check_noncovalent_bonds(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.noncovalent_bonds().iter().all(|q| {
        let mapped_a = assignment[q.atoms[0].index()];
        let mapped_b = assignment[q.atoms[1].index()];
        target.noncovalent_bonds().iter().any(|t| {
            t.atoms[0].index() == mapped_a
                && t.atoms[1].index() == mapped_b
                && q.data.matches_ground(t.data)
        })
    })
}

fn check_aromatic_systems(query: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    query.aromatic_systems().iter().all(|q| {
        let mapped: Vec<usize> = q.atoms().map(|a| assignment[a.index()]).collect();
        target.aromatic_systems().iter().any(|t| {
            let t_atoms: Vec<usize> = t.atoms().map(|a| a.index()).collect();
            mapped.iter().all(|m| t_atoms.contains(m))
        })
    })
}

fn check_multicenter_bonds(
    query: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    query.multicenter_bonds().iter().all(|q| {
        let mapped: Vec<usize> = q.atoms().map(|a| assignment[a.index()]).collect();
        target.multicenter_bonds().iter().any(|t| {
            let t_atoms: Vec<usize> = t.atoms().map(|a| a.index()).collect();
            mapped.iter().all(|m| t_atoms.contains(m))
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
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::multicenter::MulticenterBondAst;
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
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(find(&mol, &mol), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_substructure() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
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
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
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
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_wildcard_element() {
        let query = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target = mol_with_atoms(&[Element::C, Element::N]);
        let results = sorted(find(&query, &target));
        assert_eq!(results, vec![Assignment(vec![0]), Assignment(vec![1])]);
    }

    #[test]
    fn test_find_matches_wildcard_bond_order() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(find(&query, &target), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_optional_field_unconstrained() {
        let query = mol_with_atoms(&[Element::C]);
        let target = MoleculeAst::new(
            vec![AtomAst {
                charge: ValueAst::Lit(-1),
                ..AtomAst::from_element(Element::C)
            }],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(find(&query, &target), vec![Assignment(vec![0])]);
    }

    #[test]
    fn test_find_matches_optional_field_present_vs_absent() {
        let query = MoleculeAst::new(
            vec![AtomAst {
                charge: ValueAst::Lit(-1),
                ..AtomAst::from_element(Element::C)
            }],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_dative_bond_direction() {
        let atoms = vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::B),
        ];
        let query = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target_correct = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target_reversed = MoleculeAst::new(
            atoms,
            vec![],
            vec![(AtomIdx(1), AtomIdx(0), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
        );

        assert_eq!(find(&query, &target_correct), vec![Assignment(vec![0, 1])]);
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_noncovalent_bond_direction() {
        let atoms = vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
        ];
        let query = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
        );
        let target_correct = MoleculeAst::new(
            atoms.clone(),
            vec![],
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
        );
        let target_reversed = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![(AtomIdx(1), AtomIdx(0), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
        );

        assert_eq!(find(&query, &target_correct), vec![Assignment(vec![0, 1])]);
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_aromatic_system_subset() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
            vec![],
            vec![],
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
            vec![],
            vec![],
            vec![(
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
                AromaticSystemAst::default(),
            )],
            vec![],
            vec![],
        );

        assert!(!find(&query, &target).is_empty());
    }

    #[test]
    fn test_find_matches_aromatic_system_mismatch() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
            vec![],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );

        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_multicenter_subset() {
        let query = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], MulticenterBondAst {})],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::B),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
                MulticenterBondAst {},
            )],
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
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
                MulticenterBondAst {},
            )],
            vec![],
        );
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::H),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![
                (
                    vec![AtomIdx(0), AtomIdx(2), AtomIdx(1)],
                    MulticenterBondAst {},
                ),
                (
                    vec![AtomIdx(0), AtomIdx(3), AtomIdx(1)],
                    MulticenterBondAst {},
                ),
            ],
            vec![],
        );

        let results = find(&query, &target);
        assert!(!results.is_empty());
        for a in &results {
            let mapped: Vec<usize> = vec![a.0[0], a.0[1], a.0[2]];
            let in_single = target.multicenter_bonds().iter().any(|t| {
                let t_atoms: Vec<usize> = t.atoms().map(|a| a.index()).collect();
                mapped.iter().all(|m| t_atoms.contains(m))
            });
            assert!(
                in_single,
                "assignment {a:?} spans multiple multicenter bonds"
            );
        }
    }

    #[test]
    fn test_find_matches_empty_query() {
        let query = MoleculeAst::default();
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
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
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let target = mol_with_atoms(&[Element::C]);
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_with_identity() {
        let query = mol_with_atoms(&[Element::C]);
        let target = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );

        let chemistry = Chemistry::default();
        let q = MatchQuery::new(&query);
        let t = MatchTarget::new(&target);
        assert_eq!(find_matches_with(&q, &t, &chemistry), find_matches(&q, &t));
    }
}
