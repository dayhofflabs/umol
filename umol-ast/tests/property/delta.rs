use proptest::prelude::*;
use umol_ast::ast::{
    AtomDelta, AtomFieldChange, AtomId, BondDelta, BondFieldChange, BondId, Canonicalize,
    Constraint, ConstraintDelta, Delta, Deltas, MoleculeConstraint, ValueAst,
};

use crate::strategies::*;

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
}
