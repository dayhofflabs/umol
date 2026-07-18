use proptest::prelude::*;
use umol_ast::ast::{
    AromaticSystemDelta, AromaticSystemId, AtomDelta, AtomFieldChange, AtomId, BondDelta,
    BondFieldChange, BondId, Canonicalize, Constraint, ConstraintDelta, DativeBondDelta,
    DativeBondId, Delta, Deltas, EntityPatch, MoleculeConstraint, ValueAst,
};

use crate::strategies::*;

/// Apply a `diff` (only `ModifyField` / `ModifyConstraint` deltas) to an atom state.
fn apply_atom_diff(mut ast: AtomAst, diff: Vec<AtomDelta>) -> AtomAst {
    for d in diff {
        match d {
            AtomDelta::ModifyField { change, .. } => {
                AtomDelta::apply_field(&mut ast, change).unwrap()
            }
            AtomDelta::ModifyConstraint { old, new, .. } => {
                AtomDelta::apply_constraint(&mut ast, old, new).unwrap()
            }
            other => unreachable!("diff yields only modify deltas, got {other:?}"),
        }
    }
    ast
}

/// Apply a `diff` (only `ModifyField` / `ModifyConstraint` deltas) to a bond state.
fn apply_bond_diff(mut ast: BondAst, diff: Vec<BondDelta>) -> BondAst {
    for d in diff {
        match d {
            BondDelta::ModifyField { change, .. } => {
                BondDelta::apply_field(&mut ast, change).unwrap()
            }
            BondDelta::ModifyConstraint { old, new, .. } => {
                BondDelta::apply_constraint(&mut ast, old, new).unwrap()
            }
            other => unreachable!("diff yields only modify deltas, got {other:?}"),
        }
    }
    ast
}

fn apply_dative_bond_diff(mut ast: DativeBondAst, diff: Vec<DativeBondDelta>) -> DativeBondAst {
    for delta in diff {
        match delta {
            DativeBondDelta::ModifyField { change, .. } => {
                DativeBondDelta::apply_field(&mut ast, change).unwrap()
            }
            DativeBondDelta::ModifyConstraint { old, new, .. } => {
                DativeBondDelta::apply_constraint(&mut ast, old, new).unwrap()
            }
            other => unreachable!("diff yields only modify deltas, got {other:?}"),
        }
    }
    ast
}

fn apply_aromatic_system_diff(
    mut ast: AromaticSystemAst,
    diff: Vec<AromaticSystemDelta>,
) -> AromaticSystemAst {
    for delta in diff {
        match delta {
            AromaticSystemDelta::ModifyField { change, .. } => {
                AromaticSystemDelta::apply_field(&mut ast, change).unwrap()
            }
            AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
                AromaticSystemDelta::apply_constraint(&mut ast, old, new).unwrap()
            }
            other => unreachable!("diff yields only modify deltas, got {other:?}"),
        }
    }
    ast
}

fn atom_id_strategy() -> impl Strategy<Value = AtomId> {
    (0u32..3).prop_map(AtomId)
}

fn bond_id_strategy() -> impl Strategy<Value = BondId> {
    (0u32..3).prop_map(BondId)
}

fn atoms_strategy() -> impl Strategy<Value = [AtomId; 2]> {
    (0u32..3, 0u32..3).prop_map(|(a, b)| [AtomId(a), AtomId(b)])
}

fn atom_delta_strategy() -> impl Strategy<Value = AtomDelta> {
    prop_oneof![
        (atom_id_strategy(), atom_ast_strategy()).prop_map(|(id, ast)| AtomDelta::Add { id, ast }),
        (atom_id_strategy(), atom_ast_strategy())
            .prop_map(|(id, ast)| AtomDelta::Remove { id, ast }),
        (atom_id_strategy(), value_basic(0..=3), value_basic(0..=3)).prop_map(|(id, old, new)| {
            AtomDelta::ModifyField {
                id,
                change: AtomFieldChange::Charge { old, new },
            }
        }),
        (
            atom_id_strategy(),
            prop::option::of(atom_constraint_strategy()),
            prop::option::of(atom_constraint_strategy()),
        )
            .prop_map(|(id, old, new)| AtomDelta::ModifyConstraint { id, old, new }),
    ]
}

fn bond_delta_strategy() -> impl Strategy<Value = BondDelta> {
    prop_oneof![
        (bond_id_strategy(), atoms_strategy(), bond_ast_strategy())
            .prop_map(|(id, atoms, ast)| BondDelta::Add { id, atoms, ast }),
        (bond_id_strategy(), atoms_strategy(), bond_ast_strategy())
            .prop_map(|(id, atoms, ast)| BondDelta::Remove { id, atoms, ast }),
        (bond_id_strategy(), value_basic(1..=3), value_basic(1..=3)).prop_map(|(id, old, new)| {
            BondDelta::ModifyField {
                id,
                change: BondFieldChange::Order { old, new },
            }
        }),
        (
            bond_id_strategy(),
            prop::option::of(bond_constraint_strategy()),
            prop::option::of(bond_constraint_strategy()),
        )
            .prop_map(|(id, old, new)| BondDelta::ModifyConstraint { id, old, new }),
    ]
}

fn constraint_delta_strategy() -> impl Strategy<Value = ConstraintDelta> {
    // A small set of distinct constraints so adds/removes collide and exercise the
    // multiset fold.
    let constraint = (0i64..3).prop_map(|sum| {
        Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: ValueAst::Lit(sum),
        })
    });
    prop_oneof![
        constraint.clone().prop_map(ConstraintDelta::Add),
        constraint.prop_map(ConstraintDelta::Remove),
    ]
}

fn deltas_strategy() -> impl Strategy<Value = Deltas> {
    let delta = prop_oneof![
        atom_delta_strategy().prop_map(Delta::Atom),
        bond_delta_strategy().prop_map(Delta::Bond),
        constraint_delta_strategy().prop_map(Delta::Constraint),
    ];
    prop::collection::vec(delta, 0..8).prop_map(Deltas::from_iter)
}

proptest! {
    /// Canonicalize is idempotent: re-canonicalizing a canonical `Deltas` is a fixpoint.
    /// This is the confluence check — the normal form is unique, independent of the
    /// order the deltas arrived in.
    #[test]
    fn test_deltas_canonicalize_idempotent(deltas in deltas_strategy()) {
        if let Ok(once) = deltas.canonicalize() {
            prop_assert_eq!(once.clone().canonicalize().unwrap(), once);
        }
    }

    /// `apply(lhs, diff(lhs, rhs)) == rhs` for atoms — the patch algebra law.
    #[test]
    fn test_atom_delta_diff_apply(lhs in atom_ast_strategy(), rhs in atom_ast_strategy()) {
        let diff = AtomDelta::diff(AtomId(0), &lhs, &rhs);
        prop_assert_eq!(apply_atom_diff(lhs, diff), rhs);
    }

    /// `diff(x, x)` is empty and applying it is the identity (atoms).
    #[test]
    fn test_atom_delta_diff_identity(atom in atom_ast_strategy()) {
        let diff = AtomDelta::diff(AtomId(0), &atom, &atom);
        prop_assert!(diff.is_empty());
        prop_assert_eq!(apply_atom_diff(atom.clone(), diff), atom);
    }

    /// Applying the directed atom update recovers the target up to canonical equality.
    #[test]
    fn test_atom_ast_difference_to(lhs in atom_ast_strategy(), rhs in atom_ast_strategy()) {
        let update = lhs.difference_to(&rhs);
        prop_assert!(lhs.update(&update).canonical_eq(&rhs));
    }

    /// `apply(lhs, diff(lhs, rhs)) == rhs` for bonds.
    #[test]
    fn test_bond_delta_diff_apply(lhs in bond_ast_strategy(), rhs in bond_ast_strategy()) {
        let diff = BondDelta::diff(BondId(0), &lhs, &rhs);
        prop_assert_eq!(apply_bond_diff(lhs, diff), rhs);
    }

    /// `diff(x, x)` is empty and applying it is the identity (bonds).
    #[test]
    fn test_bond_delta_diff_identity(bond in bond_ast_strategy()) {
        let diff = BondDelta::diff(BondId(0), &bond, &bond);
        prop_assert!(diff.is_empty());
        prop_assert_eq!(apply_bond_diff(bond.clone(), diff), bond);
    }

    /// Applying the directed bond update recovers the target up to canonical equality.
    #[test]
    fn test_bond_ast_difference_to(lhs in bond_ast_strategy(), rhs in bond_ast_strategy()) {
        let update = lhs.difference_to(&rhs);
        prop_assert!(lhs.update(&update).canonical_eq(&rhs));
    }

    #[test]
    fn test_dative_bond_delta_diff_apply(
        lhs in dative_bond_strategy(),
        rhs in dative_bond_strategy(),
    ) {
        let diff = DativeBondDelta::diff(DativeBondId(0), &lhs, &rhs);
        prop_assert_eq!(apply_dative_bond_diff(lhs, diff), rhs);
    }

    #[test]
    fn test_dative_bond_delta_diff_identity(bond in dative_bond_strategy()) {
        let diff = DativeBondDelta::diff(DativeBondId(0), &bond, &bond);
        prop_assert!(diff.is_empty());
        prop_assert_eq!(apply_dative_bond_diff(bond.clone(), diff), bond);
    }

    #[test]
    fn test_dative_bond_ast_difference_to(
        lhs in dative_bond_strategy(),
        rhs in dative_bond_strategy(),
    ) {
        let update = lhs.difference_to(&rhs);
        prop_assert!(lhs.update(&update).canonical_eq(&rhs));
    }

    #[test]
    fn test_aromatic_system_delta_diff_apply(
        lhs in aromatic_system_patch_ast_strategy(),
        rhs in aromatic_system_patch_ast_strategy(),
    ) {
        let diff = AromaticSystemDelta::diff(AromaticSystemId(0), &lhs, &rhs);
        prop_assert_eq!(apply_aromatic_system_diff(lhs, diff), rhs);
    }

    #[test]
    fn test_aromatic_system_delta_diff_identity(system in aromatic_system_patch_ast_strategy()) {
        let diff = AromaticSystemDelta::diff(AromaticSystemId(0), &system, &system);
        prop_assert!(diff.is_empty());
        prop_assert_eq!(apply_aromatic_system_diff(system.clone(), diff), system);
    }

    #[test]
    fn test_aromatic_system_ast_difference_to(
        lhs in aromatic_system_patch_ast_strategy(),
        rhs in aromatic_system_patch_ast_strategy(),
    ) {
        let update = lhs.difference_to(&rhs);
        prop_assert!(lhs.update(&update).canonical_eq(&rhs));
    }
}
