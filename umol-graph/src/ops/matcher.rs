//! Substructure matcher over `MoleculePattern` / `Molecule`.

use umol_graph_core::{subgraph_isomorphisms, EdgeId, NodeId};

use crate::api::molecule::Molecule;
use crate::api::pattern::MoleculePattern;
use crate::ast::constraint::MoleculeConstraint;
use crate::ast::molecule::MoleculeAst;
use crate::ast::{AtomIdx, BondIdx};

/// Query atom index → target atom index mapping.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assignment(pub Vec<usize>);

/// Substructure matcher.
pub struct Matcher;

impl Matcher {
    pub fn new() -> Self {
        Self
    }

    /// All substructure matches of `pattern` in `target`.
    pub fn find(&self, pattern: &MoleculePattern, target: &Molecule) -> Vec<Assignment> {
        let pattern_ast = pattern.ast();
        let target_ast = target.ast();

        if pattern_ast.atoms().count() == 0 {
            return if post_filter(pattern_ast, target, &[]) {
                vec![Assignment(vec![])]
            } else {
                vec![]
            };
        }

        let mut node_match = |q: NodeId, t: NodeId| {
            pattern_ast
                .atom(q.into())
                .data
                .matches_ground(target_ast.atom(t.into()).data)
        };
        let mut edge_match = |q: EdgeId, t: EdgeId| {
            pattern_ast
                .bond(q.into())
                .data
                .matches_ground(target_ast.bond(t.into()).data)
        };

        subgraph_isomorphisms(
            pattern_ast.graph(),
            target_ast.graph(),
            &mut node_match,
            &mut edge_match,
        )
        .into_iter()
        .filter(|a| post_filter(pattern_ast, target, a))
        .map(Assignment)
        .collect()
    }
}

impl Default for Matcher {
    fn default() -> Self {
        Self::new()
    }
}

fn post_filter(pattern: &MoleculeAst, target: &Molecule, assignment: &[usize]) -> bool {
    let target_ast = target.ast();
    check_dative_bonds(pattern, target_ast, assignment)
        && check_noncovalent_bonds(pattern, target_ast, assignment)
        && check_aromatic_systems(pattern, target_ast, assignment)
        && check_multicenter_bonds(pattern, target_ast, assignment)
        && check_constraints(pattern, target, assignment)
}

fn check_constraints(pattern: &MoleculeAst, target: &Molecule, assignment: &[usize]) -> bool {
    pattern
        .constraints()
        .iter()
        .all(|c| evaluate_remapped(&c, pattern, target, assignment))
}

fn evaluate_remapped(
    constraint: &MoleculeConstraint,
    pattern: &MoleculeAst,
    target: &Molecule,
    assignment: &[usize],
) -> bool {
    match constraint {
        MoleculeConstraint::AtomPred(idx, c) => {
            let mapped = AtomIdx(assignment[idx.index()] as u32);
            c.evaluate(target, mapped)
        }
        MoleculeConstraint::BondPred(idx, c) => {
            let Some(mapped) = remap_bond(pattern, target, *idx, assignment) else {
                return false;
            };
            c.evaluate(target, mapped)
        }
        MoleculeConstraint::Connected(atoms) => {
            let mapped: Vec<AtomIdx> = atoms
                .iter()
                .map(|a| AtomIdx(assignment[a.index()] as u32))
                .collect();
            MoleculeConstraint::Connected(mapped).evaluate(target)
        }
        MoleculeConstraint::BondOrderSum(bonds, v) => {
            let mapped: Option<Vec<BondIdx>> = bonds
                .iter()
                .map(|b| remap_bond(pattern, target, *b, assignment))
                .collect();
            match mapped {
                Some(mb) => MoleculeConstraint::BondOrderSum(mb, v.clone()).evaluate(target),
                None => false,
            }
        }
        MoleculeConstraint::TotalCharge(_) | MoleculeConstraint::TotalSpin(_) => {
            constraint.evaluate(target)
        }
        MoleculeConstraint::AromaticElectronCount(_, _)
        | MoleculeConstraint::MulticenterElectronCount(_, _) => true,
        MoleculeConstraint::SubPattern { .. } => false,
        MoleculeConstraint::And(xs) => xs
            .iter()
            .all(|x| evaluate_remapped(x, pattern, target, assignment)),
        MoleculeConstraint::Or(xs) => xs
            .iter()
            .any(|x| evaluate_remapped(x, pattern, target, assignment)),
        MoleculeConstraint::Not(inner) => !evaluate_remapped(inner, pattern, target, assignment),
    }
}

fn remap_bond(
    pattern: &MoleculeAst,
    target: &Molecule,
    bond: BondIdx,
    assignment: &[usize],
) -> Option<BondIdx> {
    let view = pattern.bond(bond);
    let src: NodeId = AtomIdx(assignment[view.src.index()] as u32).into();
    let tgt: NodeId = AtomIdx(assignment[view.tgt.index()] as u32).into();
    target.graph().find_edge(src, tgt).map(Into::into)
}

fn check_dative_bonds(pattern: &MoleculeAst, target: &MoleculeAst, assignment: &[usize]) -> bool {
    pattern.dative_bonds().iter().all(|q| {
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
    pattern: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    pattern.noncovalent_bonds().iter().all(|q| {
        let mapped_a = assignment[q.atoms[0].index()];
        let mapped_b = assignment[q.atoms[1].index()];
        target.noncovalent_bonds().iter().any(|t| {
            t.atoms[0].index() == mapped_a
                && t.atoms[1].index() == mapped_b
                && q.data.matches_ground(t.data)
        })
    })
}

fn check_aromatic_systems(
    pattern: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    pattern.aromatic_systems().iter().all(|q| {
        let mapped: Vec<usize> = q.atoms().map(|a| assignment[a.index()]).collect();
        target.aromatic_systems().iter().any(|t| {
            let t_atoms: Vec<usize> = t.atoms().map(|a| a.index()).collect();
            mapped.iter().all(|m| t_atoms.contains(m))
        })
    })
}

fn check_multicenter_bonds(
    pattern: &MoleculeAst,
    target: &MoleculeAst,
    assignment: &[usize],
) -> bool {
    pattern.multicenter_bonds().iter().all(|q| {
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
    use crate::api::pattern::MoleculePattern;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::config::MoleculeAstConfig;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::AtomIdx;

    fn mol_with_atoms(elements: &[Element]) -> MoleculeAst {
        let atoms = elements.iter().map(|&e| AtomAst::from_element(e)).collect();
        MoleculeAst::new(atoms, vec![], vec![], vec![], vec![], vec![], vec![])
    }

    fn find(pattern: &MoleculeAst, target: &MoleculeAst) -> Vec<Assignment> {
        let pattern = MoleculePattern::new(pattern.clone());
        let mut grounded = target.clone();
        grounded.coerce(&MoleculeAstConfig::zeroed());
        let target = Molecule::new(grounded).unwrap();
        Matcher::new().find(&pattern, &target)
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

    mod constraint_eval {
        use pretty_assertions::assert_eq;
        use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
        use umol_shared::element::Element;
        use umol_shared::spin::SpinState;
        use umol_shared::spin_ast::SpinStateAst;
        use umol_shared::value_ast::ValueAst;

        use super::super::*;
        use crate::api::molecule::Molecule;
        use crate::api::pattern::MoleculePattern;
        use crate::ast::atom::AtomAst;
        use crate::ast::bond::BondAst;
        use crate::ast::constraint::{AtomConstraint, BondConstraint, MoleculeConstraint};
        use crate::ast::{AtomIdx, BondIdx};

        fn ground_c(h: i64) -> AtomAst {
            AtomAst {
                element: ElementAst::Lit(Element::C),
                isotope_mass: IsotopeAst::Natural,
                charge: ValueAst::Lit(0),
                implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(h)),
                lone_pairs: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            }
        }

        fn ground_bond(order: i64) -> BondAst {
            BondAst {
                order: ValueAst::Lit(order),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::Lit(SpinState::closed_shell()),
            }
        }

        fn propane() -> Molecule {
            let ast = MoleculeAst::new(
                vec![ground_c(3), ground_c(2), ground_c(3)],
                vec![
                    (AtomIdx(0), AtomIdx(1), ground_bond(1)),
                    (AtomIdx(1), AtomIdx(2), ground_bond(1)),
                ],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            );
            Molecule::new(ast).unwrap()
        }

        #[test]
        fn test_find_matches_atom_predicate_filters() {
            let target = propane();
            let mut ast = MoleculeAst::new(
                vec![AtomAst::new(ElementAst::Undetermined)],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            );
            ast.constraints_mut().insert(MoleculeConstraint::AtomPred(
                AtomIdx(0),
                AtomConstraint::Degree(ValueAst::Lit(2)),
            ));
            let pattern = MoleculePattern::new(ast);
            let matches = Matcher::new().find(&pattern, &target);
            assert_eq!(matches, vec![Assignment(vec![1])]);
        }

        #[test]
        fn test_find_matches_bond_predicate_filters() {
            let cyclopropane = {
                let ast = MoleculeAst::new(
                    vec![ground_c(2), ground_c(2), ground_c(2)],
                    vec![
                        (AtomIdx(0), AtomIdx(1), ground_bond(1)),
                        (AtomIdx(1), AtomIdx(2), ground_bond(1)),
                        (AtomIdx(2), AtomIdx(0), ground_bond(1)),
                    ],
                    vec![],
                    vec![],
                    vec![],
                    vec![],
                    vec![],
                );
                Molecule::new(ast).unwrap()
            };
            let mut ast = MoleculeAst::new(
                vec![
                    AtomAst::new(ElementAst::Undetermined),
                    AtomAst::new(ElementAst::Undetermined),
                ],
                vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            );
            ast.constraints_mut().insert(MoleculeConstraint::BondPred(
                BondIdx(0),
                BondConstraint::RingBond,
            ));
            let pattern = MoleculePattern::new(ast);
            let matches = Matcher::new().find(&pattern, &cyclopropane);
            assert!(!matches.is_empty());
            assert!(matches.iter().all(|a| a.0.len() == 2));
        }

        #[test]
        fn test_find_matches_atom_predicate_no_match() {
            let target = propane();
            let mut ast = MoleculeAst::new(
                vec![AtomAst::new(ElementAst::Undetermined)],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            );
            ast.constraints_mut().insert(MoleculeConstraint::AtomPred(
                AtomIdx(0),
                AtomConstraint::Degree(ValueAst::Lit(5)),
            ));
            let pattern = MoleculePattern::new(ast);
            assert_eq!(Matcher::new().find(&pattern, &target), vec![]);
        }

        #[test]
        fn test_find_matches_and_combinator() {
            let target = propane();
            let mut ast = MoleculeAst::new(
                vec![AtomAst::new(ElementAst::Undetermined)],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            );
            ast.constraints_mut()
                .insert(MoleculeConstraint::And(vec![
                    MoleculeConstraint::AtomPred(
                        AtomIdx(0),
                        AtomConstraint::Degree(ValueAst::Lit(1)),
                    ),
                    MoleculeConstraint::AtomPred(
                        AtomIdx(0),
                        AtomConstraint::TotalHCount(ValueAst::Lit(3)),
                    ),
                ]));
            let pattern = MoleculePattern::new(ast);
            let matches = Matcher::new().find(&pattern, &target);
            let mut idx: Vec<usize> = matches.iter().map(|a| a.0[0]).collect();
            idx.sort();
            assert_eq!(idx, vec![0, 2]);
        }
    }
}
