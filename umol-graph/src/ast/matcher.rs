//! Substructure matcher over MoleculeAst using VF2 subgraph isomorphism.

use petgraph::algo::subgraph_isomorphisms_iter;
use petgraph::graph::Graph;
use petgraph::Directed;

use crate::ast::molecule::{AromaticSystem, BondTuple, MoleculeAst, MulticenterBond};

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
    let mut graph = Graph::with_capacity(ast.atoms.len(), ast.bonds.len() * 2);
    let nodes: Vec<_> = (0..ast.atoms.len()).map(|i| graph.add_node(i)).collect();
    for (bond_idx, bond) in ast.bonds.iter().enumerate() {
        graph.add_edge(nodes[bond.source], nodes[bond.target], bond_idx);
        graph.add_edge(nodes[bond.target], nodes[bond.source], bond_idx);
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
}

impl<'a> MatchQuery<'a> {
    pub fn new(ast: &'a MoleculeAst) -> Self {
        Self {
            graph: build_graph(ast),
            ast,
        }
    }
}

/// Find all substructure matches of `query` in `target`.
pub fn find_matches(query: &MatchQuery, target: &MatchTarget) -> Vec<Assignment> {
    let q_ast = query.ast;
    let t_ast = target.ast;

    if q_ast.atoms.is_empty() {
        return if post_filter(q_ast, t_ast, &[]) {
            vec![Assignment(vec![])]
        } else {
            vec![]
        };
    }

    let mut node_match = |q_idx: &usize, t_idx: &usize| {
        q_ast.atoms[*q_idx].matches_ground(&t_ast.atoms[*t_idx])
    };
    let mut edge_match = |q_idx: &usize, t_idx: &usize| {
        q_ast.bonds[*q_idx].bond.matches_ground(&t_ast.bonds[*t_idx].bond)
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

fn post_filter(query: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    check_directed_bonds(&query.dative_bonds, &target.dative_bonds, assignment)
        && check_directed_bonds(&query.noncovalent_bonds, &target.noncovalent_bonds, assignment)
        && check_group_subset(&query.aromatic_systems, &target.aromatic_systems, assignment)
        && check_group_subset(&query.multicenter_bonds, &target.multicenter_bonds, assignment)
}

fn check_directed_bonds(
    query_bonds: &[BondTuple],
    target_bonds: &[BondTuple],
    assignment: &[usize],
) -> bool {
    query_bonds.iter().all(|qb| {
        let mapped_source = assignment[qb.source];
        let mapped_target = assignment[qb.target];
        target_bonds.iter().any(|tb| {
            tb.source == mapped_source
                && tb.target == mapped_target
                && qb.bond.matches_ground(&tb.bond)
        })
    })
}

fn check_group_subset<T: HasAtoms>(
    query_groups: &[T],
    target_groups: &[T],
    assignment: &[usize],
) -> bool {
    query_groups.iter().all(|qg| {
        let mapped: Vec<usize> = qg.atoms().iter().map(|&a| assignment[a]).collect();
        target_groups
            .iter()
            .any(|tg| mapped.iter().all(|m| tg.atoms().contains(m)))
    })
}

trait HasAtoms {
    fn atoms(&self) -> &[usize];
}

impl HasAtoms for AromaticSystem {
    fn atoms(&self) -> &[usize] {
        &self.atoms
    }
}

impl HasAtoms for MulticenterBond {
    fn atoms(&self) -> &[usize] {
        &self.atoms
    }
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
    use crate::ast::molecule::{AromaticSystem, BondTuple, MulticenterBond};

    fn atom(e: Element) -> AtomAst {
        AtomAst::from_element(e)
    }

    fn bond(source: usize, target: usize, order: u8) -> BondTuple {
        BondTuple {
            source,
            target,
            bond: BondAst::from_order(order),
        }
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
        let mol = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        assert_eq!(find(&mol, &mol), vec![Assignment(vec![0])]);
    }

    #[test]
    fn test_find_matches_chain_identity() {
        let mol = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        assert_eq!(find(&mol, &mol), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_substructure() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1), bond(1, 2, 1)],
            ..Default::default()
        };
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
        let query = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::N)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_bond_order_mismatch() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![bond(0, 1, 2)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_wildcard_element() {
        let query = MoleculeAst {
            atoms: vec![AtomAst::new(ElementAst::Wildcard)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::N)],
            ..Default::default()
        };
        let results = sorted(find(&query, &target));
        assert_eq!(
            results,
            vec![Assignment(vec![0]), Assignment(vec![1])]
        );
    }

    #[test]
    fn test_find_matches_wildcard_bond_order() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![BondTuple {
                source: 0,
                target: 1,
                bond: BondAst::new(ValueAst::Wildcard),
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![bond(0, 1, 2)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![Assignment(vec![0, 1])]);
    }

    #[test]
    fn test_find_matches_optional_field_unconstrained() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![AtomAst {
                charge: Some(ValueAst::Lit(-1)),
                ..atom(Element::C)
            }],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![Assignment(vec![0])]);
    }

    #[test]
    fn test_find_matches_optional_field_present_vs_absent() {
        let query = MoleculeAst {
            atoms: vec![AtomAst {
                charge: Some(ValueAst::Lit(-1)),
                ..atom(Element::C)
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_dative_bond_direction() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::B)],
            dative_bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        let target_correct = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::B)],
            dative_bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        let target_reversed = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::B)],
            dative_bonds: vec![bond(1, 0, 1)],
            ..Default::default()
        };
        assert_eq!(
            find(&query, &target_correct),
            vec![Assignment(vec![0, 1])]
        );
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_noncovalent_bond_direction() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::O)],
            noncovalent_bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        let target_correct = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::O)],
            noncovalent_bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        let target_reversed = MoleculeAst {
            atoms: vec![atom(Element::N), atom(Element::O)],
            noncovalent_bonds: vec![bond(1, 0, 1)],
            ..Default::default()
        };
        assert_eq!(
            find(&query, &target_correct),
            vec![Assignment(vec![0, 1])]
        );
        assert_eq!(find(&query, &target_reversed), vec![]);
    }

    #[test]
    fn test_find_matches_aromatic_system_subset() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1)],
            aromatic_systems: vec![AromaticSystem {
                atoms: vec![0, 1],
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1), bond(1, 2, 1), bond(2, 0, 1)],
            aromatic_systems: vec![AromaticSystem {
                atoms: vec![0, 1, 2],
            }],
            ..Default::default()
        };
        let results = find(&query, &target);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_find_matches_aromatic_system_mismatch() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1)],
            aromatic_systems: vec![AromaticSystem {
                atoms: vec![0, 1],
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![]);
    }

    #[test]
    fn test_find_matches_multicenter_subset() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::B), atom(Element::H)],
            multicenter_bonds: vec![MulticenterBond {
                atoms: vec![0, 1],
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::B), atom(Element::H), atom(Element::B)],
            multicenter_bonds: vec![MulticenterBond {
                atoms: vec![0, 1, 2],
            }],
            ..Default::default()
        };
        let results = find(&query, &target);
        assert!(!results.is_empty());
    }

    #[test]
    fn test_find_matches_multicenter_identity() {
        // B2H6: two 3c-2e bonds sharing the B-B pair
        let query = MoleculeAst {
            atoms: vec![atom(Element::B), atom(Element::H), atom(Element::B)],
            multicenter_bonds: vec![MulticenterBond {
                atoms: vec![0, 1, 2],
            }],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![
                atom(Element::B), // 0
                atom(Element::B), // 1
                atom(Element::H), // 2 (bridge 1)
                atom(Element::H), // 3 (bridge 2)
            ],
            multicenter_bonds: vec![
                MulticenterBond {
                    atoms: vec![0, 2, 1],
                },
                MulticenterBond {
                    atoms: vec![0, 3, 1],
                },
            ],
            ..Default::default()
        };
        let results = find(&query, &target);
        // Each 3c-2e bond should produce matches (B,H,B permutations)
        assert!(!results.is_empty());
        // All assignments must map into a single multicenter bond
        for a in &results {
            let mapped: Vec<usize> = vec![a.0[0], a.0[1], a.0[2]];
            let in_single = target.multicenter_bonds.iter().any(|mc| {
                mapped.iter().all(|m| mc.atoms.contains(m))
            });
            assert!(in_single, "assignment {a:?} spans multiple multicenter bonds");
        }
    }

    #[test]
    fn test_find_matches_empty_query() {
        let query = MoleculeAst::default();
        let target = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::N)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![Assignment(vec![])]);
    }

    #[test]
    fn test_find_matches_no_match() {
        let query = MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1), bond(1, 2, 1)],
            ..Default::default()
        };
        let target = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        assert_eq!(find(&query, &target), vec![]);
    }
}
