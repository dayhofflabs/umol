//! Property-based tests for MoleculeDsl roundtrip invariants.
//!
//! Generates structured `MoleculeAst` values, wraps them in `MoleculeDsl`,
//! and asserts the render → parse cycle is the identity.

use std::collections::HashSet;
use std::iter::repeat_with;
use std::ops::RangeInclusive;

use proptest::prelude::*;
use umol_ast::ast::{
    ArithOp, AromaticSystemAst, AromaticSystemConstraint, AromaticSystemConstraints,
    AromaticSystemIdx, AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraints, AtomIdx,
    BondAst, BondConstraint, BondConstraints, BondIdx, Constraint, Constraints, DativeBondAst,
    DativeBondConstraint, DativeBondConstraints, DativeBondIdx, ElementAst, Expr,
    ImplicitHydrogensAst, IsotopeAst, MoleculeAst, MoleculeConstraint, MulticenterBondAst,
    MulticenterBondConstraint, MulticenterBondConstraints, MulticenterBondIdx,
    MulticenterValenceAst, NoncovalentBondAst, NoncovalentBondIdx, NoncovalentBondKind,
    NoncovalentBondKindAst, RelOp, RelationalConstraint, SpinStateAst, SubPatternAnchor, ValueAst,
};
use umol_ast::dsl::{
    parse_value, AromaticSystemDsl, AtomDsl, BondDsl, DativeBondDsl, Metadata, MoleculeDsl,
    MulticenterBondDsl, NoncovalentBondDsl, ValueDsl,
};
use umol_edn::{read_string, Edn, FromEdn, ToEdn};
use umol_shared::element::Element;

const ELEMENTS: &[Element] = &[
    Element::H,
    Element::C,
    Element::N,
    Element::O,
    Element::F,
    Element::P,
    Element::S,
    Element::Cl,
    Element::Br,
];

const NONCOVALENT_KINDS: &[NoncovalentBondKind] = &[
    NoncovalentBondKind::HydrogenBond,
    NoncovalentBondKind::HalogenBond,
    NoncovalentBondKind::ChalcogenBond,
    NoncovalentBondKind::Ionic,
    NoncovalentBondKind::VanDerWaals,
];

fn element_strategy() -> impl Strategy<Value = Element> {
    prop::sample::select(ELEMENTS)
}

fn id_strategy() -> impl Strategy<Value = String> {
    "[a-zA-Z][a-zA-Z0-9_]{0,3}".prop_map(|s| s.to_string())
}

fn element_ast_strategy() -> impl Strategy<Value = ElementAst> {
    prop_oneof![
        6 => element_strategy().prop_map(ElementAst::Lit),
        2 => Just(ElementAst::Undetermined),
        2 => prop::collection::vec(element_strategy(), 1..=3).prop_map(|mut v| {
            // Deduplicate to keep shape canonical for roundtrip.
            v.dedup();
            ElementAst::Set(v)
        }),
        1 => (id_strategy(), prop::collection::vec(element_strategy(), 1..=3))
            .prop_map(|(id, mut set)| {
                set.dedup();
                ElementAst::Bind { id, set }
            }),
        1 => id_strategy().prop_map(ElementAst::Ref),
    ]
}

fn value_basic(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        4 => Just(ValueAst::Undetermined),
        4 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::LitSet(v)
        }),
        1 => top_expr_strategy().prop_map(ValueAst::Expr),
    ]
}

/// Arithmetic-typed Expr: produces only the arithmetic subset of `Expr`
/// (`Lit`, `Var`, `Neg(arith)`, `BinOp(arith, op, arith)`). Includes
/// negative `Lit` and `Neg(Neg(_))` shapes that the parser canonicalizes,
/// to be paired with `simplify()` for roundtrip testing.
fn arith_expr_strategy() -> BoxedStrategy<Expr> {
    let leaf = prop_oneof![
        (-10i64..=10).prop_map(Expr::Lit),
        id_strategy().prop_map(Expr::Var),
    ]
    .boxed();
    leaf.prop_recursive(3, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|e| Expr::Neg(Box::new(e))),
            (inner.clone(), arith_op_strategy(), inner).prop_map(|(l, op, r)| Expr::BinOp(
                Box::new(l),
                op,
                Box::new(r)
            )),
        ]
        .boxed()
    })
    .boxed()
}

/// Boolean-typed Expr: `Rel(arith, op, arith)`, `Mem(arith, set)`,
/// `Not(bool)`, `And(bool*)`, `Or(bool*)`. Each boolean recursion correctly
/// roots in arithmetic leaves so the parser accepts the rendered form.
fn bool_expr_strategy() -> BoxedStrategy<Expr> {
    let arith = arith_expr_strategy();
    let leaf = prop_oneof![
        (arith.clone(), rel_op_strategy(), arith.clone()).prop_map(|(l, op, r)| Expr::Rel(
            Box::new(l),
            op,
            Box::new(r)
        )),
        (arith, prop::collection::vec(-10i64..=10, 1..=3))
            .prop_map(|(e, set)| Expr::Mem(Box::new(e), set)),
    ]
    .boxed();
    leaf.prop_recursive(2, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|e| Expr::Not(Box::new(e))),
            prop::collection::vec(inner.clone(), 1..=3).prop_map(Expr::And),
            prop::collection::vec(inner, 1..=3).prop_map(Expr::Or),
        ]
        .boxed()
    })
    .boxed()
}

fn any_expr_strategy() -> BoxedStrategy<Expr> {
    prop_oneof![arith_expr_strategy(), bool_expr_strategy()].boxed()
}

fn any_value_ast_strategy() -> BoxedStrategy<ValueAst> {
    prop_oneof![
        Just(ValueAst::Undetermined),
        (-10i64..=10).prop_map(ValueAst::Lit),
        prop::collection::vec(-10i64..=10, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::LitSet(v)
        }),
        any_expr_strategy().prop_map(ValueAst::Expr),
    ]
    .boxed()
}

fn arith_op_strategy() -> impl Strategy<Value = ArithOp> {
    prop_oneof![
        Just(ArithOp::Add),
        Just(ArithOp::Sub),
        Just(ArithOp::Mul),
        Just(ArithOp::Div),
        Just(ArithOp::Rem),
    ]
}

fn rel_op_strategy() -> impl Strategy<Value = RelOp> {
    prop_oneof![
        Just(RelOp::Le),
        Just(RelOp::Ge),
        Just(RelOp::Eq),
        Just(RelOp::Lt),
        Just(RelOp::Gt),
    ]
}

/// Expr leaf: non-negative `Lit` or `Var`. Safe as a subexpression of any
/// operator; **not** safe as the outermost `ValueAst::Expr` wrapper (see
/// `top_expr_strategy`). Negative literals are excluded because the Expr
/// grammar has no `Lit(-n)` parse — `-n` inside an Expr always parses as
/// `Neg(Lit(n))`, so emitting `Expr::Lit(-n)` from the generator would fail
/// the structural roundtrip equality even though the two forms are
/// semantically identical under `is_ground` / `evaluate`. Negative integers
/// still appear elsewhere (top-level `ValueAst::Lit`, `LitSet`, and
/// `Expr::Mem` sets all route through `dec_int` which reads signed).
fn expr_leaf_strategy() -> impl Strategy<Value = Expr> {
    prop_oneof![
        (0i64..=10).prop_map(Expr::Lit),
        id_strategy().prop_map(Expr::Var),
    ]
}

/// Expr tree intended as the outermost `ValueAst::Expr(e)`. Must not render
/// to a pure integer literal with optional sign, since the `value` parser's
/// `dec_int` alt would match first (bare `Lit(n)` → `n` → `ValueAst::Lit(n)`
/// and `Neg(Lit(n))` → `-n` → `ValueAst::Lit(-n)`). Also avoids `Neg(Neg(_))`
/// anywhere (the parser folds consecutive signs) and `Or`/`And` children
/// that are themselves `Or`/`And` (the parser flattens consecutive same-op
/// tokens).
///
/// Returns a boxed strategy to keep the composed type size bounded when
/// plugged into every value field across the molecule tree.
fn top_expr_strategy() -> BoxedStrategy<Expr> {
    let set = prop::collection::vec(-10i64..=10, 1..=3);
    prop_oneof![
        // Var alone renders as `?id`, distinguishable from a bare integer.
        id_strategy().prop_map(Expr::Var),
        (
            expr_leaf_strategy(),
            arith_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| Expr::BinOp(Box::new(a), op, Box::new(b))),
        (
            expr_leaf_strategy(),
            rel_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| Expr::Rel(Box::new(a), op, Box::new(b))),
        (expr_leaf_strategy(), set).prop_map(|(e, s)| Expr::Mem(Box::new(e), s)),
        (
            expr_leaf_strategy(),
            rel_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| {
                Expr::Not(Box::new(Expr::Rel(Box::new(a), op, Box::new(b))))
            }),
        // `Neg` of `Var` renders `-?id`; safe (non-Lit inner, no sign folding).
        id_strategy().prop_map(|x| Expr::Neg(Box::new(Expr::Var(x)))),
        // `Or` / `And` with exactly leaf children so the parser can't flatten.
        prop::collection::vec(expr_leaf_strategy(), 2..=3).prop_map(Expr::Or),
        prop::collection::vec(expr_leaf_strategy(), 2..=3).prop_map(Expr::And),
    ]
    .boxed()
}

fn isotope_strategy() -> impl Strategy<Value = IsotopeAst> {
    prop_oneof![
        3 => Just(IsotopeAst::Natural),
        3 => Just(IsotopeAst::Undetermined),
        3 => (1i64..=250).prop_map(IsotopeAst::Lit),
        1 => prop::collection::vec(1i64..=250, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            IsotopeAst::LitSet(v)
        }),
        1 => top_expr_strategy().prop_map(IsotopeAst::Expr),
    ]
}

fn implicit_hydrogens_strategy() -> impl Strategy<Value = ImplicitHydrogensAst> {
    prop_oneof![
        3 => Just(ImplicitHydrogensAst::Normal),
        3 => Just(ImplicitHydrogensAst::Undetermined),
        3 => (0i64..=4).prop_map(ImplicitHydrogensAst::Lit),
        1 => prop::collection::vec(0i64..=4, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ImplicitHydrogensAst::LitSet(v)
        }),
        1 => top_expr_strategy().prop_map(ImplicitHydrogensAst::Expr),
    ]
}

fn spin_state_strategy() -> impl Strategy<Value = SpinStateAst> {
    // DSL preserves spin fields field-wise. Physical (u, m) parity is a
    // tier-2 solver invariant, not a parse-time check, so any independent
    // pair must roundtrip.
    (value_basic(0..=6), value_basic(1..=7)).prop_map(|(u, m)| SpinStateAst { unpaired: u, multiplicity: m })
}

/// `SpinStateAst` with at least one of `unpaired` / `multiplicity` not
/// `Undetermined`. Used inside `MoleculeConstraint::SpinSum` and similar
/// where a fully-vacuous spin state would elide on render.
fn non_vacuous_spin_state_strategy() -> impl Strategy<Value = SpinStateAst> {
    (value_basic(0..=6), value_basic(1..=7))
        .prop_map(|(u, m)| SpinStateAst { unpaired: u, multiplicity: m })
        .prop_filter("non-vacuous spin", |s| !s.is_undetermined())
}

/// Simple value strategy used inside constraint values: `Undetermined`,
/// `Lit`, and `LitSet`. No `Expr` — the constraint formatters route to
/// `fmt_value_field_required` / `fmt_ring_count` / the various `#r` blocks,
/// and `Expr(Lit(n))` or `Expr(Neg(Lit(n)))` would render to a pure integer
/// that the parser then re-reads as a plain `Lit`, breaking roundtrip. The
/// molecule-level EDN tests cover `Expr` on constraint values through the
/// tree-based path, so the gap is contained.
fn constraint_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => Just(ValueAst::Undetermined),
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::LitSet(v)
        }),
    ]
}

/// `Lit`/`LitSet` only — still used by the ring-size strategies where
/// `Undetermined` on the inner value collapses into a dropped constraint
/// in the entity-level formatter (see `BondConstraint::RingSize` /
/// `DativeBondConstraint::RingSize` — vacuous, intentionally dropped).
fn constraint_inner_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        range.clone().prop_map(ValueAst::Lit),
        prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::LitSet(v)
        }),
    ]
}

/// `AromaticValenceAst::Undetermined` is vacuous (renders empty per the
/// canonical-rendering rule). Excluded so the strategy stays inside the
/// render → reparse identity.
fn aromatic_valence_ast_strategy() -> impl Strategy<Value = AromaticValenceAst> {
    prop_oneof![
        Just(AromaticValenceAst::NotAromatic),
        constraint_value_strategy(0..=6).prop_map(AromaticValenceAst::Aromatic),
    ]
}

fn multicenter_valence_ast_strategy() -> impl Strategy<Value = MulticenterValenceAst> {
    prop_oneof![
        Just(MulticenterValenceAst::NotMulticenter),
        constraint_value_strategy(0..=6).prop_map(MulticenterValenceAst::Multicenter),
    ]
}

/// Atom constraints route through `fmt_value_field_required` (or
/// `fmt_ring_count` for `#R`), which elide vacuous (Undetermined) payloads
/// per the canonical-rendering rule. Generators excluding `Undetermined`
/// keep the render → reparse identity intact.
fn atom_constraint_strategy() -> BoxedStrategy<AtomConstraint> {
    prop_oneof![
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::Valence),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraint::DonatedPairs),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraint::AcceptedPairs),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::Degree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::Connectivity),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::RingConnectivity),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::TotalHydrogens),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(AtomConstraint::RingSize),
        aromatic_valence_ast_strategy().prop_map(AtomConstraint::AromaticValence),
        multicenter_valence_ast_strategy().prop_map(AtomConstraint::MulticenterValence),
    ]
    .boxed()
}

fn atom_constraints_strategy() -> impl Strategy<Value = AtomConstraints> {
    prop::collection::vec(atom_constraint_strategy(), 0..=3).prop_map(|list| {
        let mut cs = AtomConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

fn bond_constraint_strategy() -> BoxedStrategy<BondConstraint> {
    prop_oneof![
        Just(BondConstraint::Aromatic),
        constraint_inner_value_strategy(0..=6).prop_map(BondConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(BondConstraint::RingSize),
    ]
    .boxed()
}

fn bond_constraints_strategy() -> impl Strategy<Value = BondConstraints> {
    prop::collection::vec(bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = BondConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

fn dative_bond_constraint_strategy() -> BoxedStrategy<DativeBondConstraint> {
    // Only the ring-membership constraints render into the entity's string
    // form; Donor/Acceptor/Parallels/{Donor,Acceptor}Satisfies carry
    // molecule-scope refs and are skipped by `fmt_constraint`.
    prop_oneof![
        constraint_inner_value_strategy(0..=6).prop_map(DativeBondConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(DativeBondConstraint::RingSize),
    ]
    .boxed()
}

fn dative_bond_constraints_strategy() -> impl Strategy<Value = DativeBondConstraints> {
    prop::collection::vec(dative_bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = DativeBondConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

prop_compose! {
    fn atom_ast_strategy()
    (
        element in element_ast_strategy(),
        isotope in isotope_strategy(),
        charge in value_basic(-2..=2),
        implicit_hydrogens in implicit_hydrogens_strategy(),
        lone_pairs in value_basic(0..=4),
        spin in spin_state_strategy(),
        constraints in atom_constraints_strategy(),
    ) -> AtomAst {
        AtomAst {
            element,
            isotope_mass: isotope,
            charge,
            implicit_hydrogens,
            lone_pairs,
            spin,
            constraints,
        }
    }
}

prop_compose! {
    fn bond_ast_strategy()
    (
        order in value_basic(1..=4),
        charge in value_basic(-1..=1),
        constraints in bond_constraints_strategy(),
    ) -> BondAst {
        BondAst {
            order,
            charge,
            spin: SpinStateAst::default(),
            constraints,
        }
    }
}

/// `BondAst` shapes that render to bond keyword shorthands per spec §7.7:
/// `:single`, `:double`, `:triple`, `:quadruple`, plus `:aromatic` (an
/// order-1 bond with the inline `Aromatic` flag).
fn canonical_keyword_bond_strategy() -> impl Strategy<Value = BondAst> {
    prop_oneof![
        Just(BondAst::new(ValueAst::Lit(1))),
        Just(BondAst::new(ValueAst::Lit(2))),
        Just(BondAst::new(ValueAst::Lit(3))),
        Just(BondAst::new(ValueAst::Lit(4))),
        Just({
            let mut bond = BondAst::new(ValueAst::Lit(1));
            bond.constraints.add(BondConstraint::Aromatic);
            bond
        }),
    ]
}

/// Generate a random undirected, simple-graph edge set over `atom_count`
/// vertices: no self-loops, no duplicates.
fn edge_set_strategy(atom_count: usize) -> impl Strategy<Value = Vec<[u32; 2]>> {
    if atom_count < 2 {
        return Just(Vec::new()).boxed();
    }
    let max_edges = atom_count.min(8);
    prop::collection::vec((0..atom_count as u32, 0..atom_count as u32), 0..=max_edges)
        .prop_map(|pairs| {
            let mut seen: HashSet<(u32, u32)> = HashSet::new();
            let mut out = Vec::new();
            for (a, b) in pairs {
                if a == b {
                    continue;
                }
                let key = if a < b { (a, b) } else { (b, a) };
                if seen.insert(key) {
                    out.push([key.0, key.1]);
                }
            }
            out
        })
        .boxed()
}

fn dative_bond_strategy() -> impl Strategy<Value = DativeBondAst> {
    // `acceptor_slot` is set by the owning molecule from endpoint ordering.
    // Order is sampled from the small literal range that the DSL keyword
    // shorthands cover (`:single` / `:double` / `:triple`), keeping
    // canonical-form roundtrip exercised across haptic-pair counts.
    let order_strategy = prop_oneof![
        Just(ValueAst::Lit(1)),
        Just(ValueAst::Lit(2)),
        Just(ValueAst::Lit(3)),
        Just(ValueAst::Undetermined),
    ];
    (order_strategy, dative_bond_constraints_strategy()).prop_map(
        |(order, constraints)| DativeBondAst {
            acceptor_slot: 0,
            order,
            constraints,
        },
    )
}

/// Optional `ElectronCount` constraint (the asserted total). The strategy
/// emits `None` half the time, otherwise wraps a `ValueAst::Lit` or
/// `LitSet`. `Undetermined` is excluded because it has no canonical
/// surface form in the entity-string `#e<n>` slot — `#e*` is admitted on
/// parse but the renderer omits the predicate entirely, breaking
/// roundtrip.
fn optional_aromatic_electron_count() -> impl Strategy<Value = AromaticSystemConstraints> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=12)).prop_map(|opt| {
        let mut cs = AromaticSystemConstraints::new();
        if let Some(v) = opt {
            cs.add(AromaticSystemConstraint::ElectronCount(v));
        }
        cs
    })
}

fn optional_multicenter_electron_count() -> impl Strategy<Value = MulticenterBondConstraints> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=8)).prop_map(|opt| {
        let mut cs = MulticenterBondConstraints::new();
        if let Some(v) = opt {
            cs.add(MulticenterBondConstraint::ElectronCount(v));
        }
        cs
    })
}

fn electron_count_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::LitSet(v)
        }),
    ]
}

/// Stand-alone strategy for entity-string roundtrip tests. The per-atom
/// `electrons` vec is empty because the entity string carries no per-atom
/// data; the `ElectronCount` constraint is exercised here via `#e<n>`.
fn aromatic_system_ast_strategy() -> impl Strategy<Value = AromaticSystemAst> {
    (value_basic(-2..=2), optional_aromatic_electron_count()).prop_map(
        |(charge, constraints)| AromaticSystemAst {
            electrons: Vec::new(),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        },
    )
}

/// Atom-count-aware variant: generates an `AromaticSystemAst` whose
/// `electrons` vec has exactly `atom_count` entries. Includes an optional
/// `ElectronCount` constraint so the molecule-level prop tests exercise
/// both the per-atom vec and the asserted total in the same pass.
fn aromatic_system_ast_for(atom_count: usize) -> impl Strategy<Value = AromaticSystemAst> {
    (
        value_basic(-2..=2),
        prop::collection::vec(value_basic(0..=2), atom_count),
        optional_aromatic_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| AromaticSystemAst {
            electrons,
            charge,
            spin: SpinStateAst::default(),
            constraints,
        })
}

fn multicenter_bond_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (value_basic(-2..=2), optional_multicenter_electron_count()).prop_map(
        |(charge, constraints)| MulticenterBondAst {
            electrons: Vec::new(),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        },
    )
}

fn multicenter_bond_ast_for(atom_count: usize) -> impl Strategy<Value = MulticenterBondAst> {
    (
        value_basic(-2..=2),
        prop::collection::vec(value_basic(0..=2), atom_count),
        optional_multicenter_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| MulticenterBondAst {
            electrons,
            charge,
            spin: SpinStateAst::default(),
            constraints,
        })
}

fn noncovalent_bond_ast_strategy() -> impl Strategy<Value = NoncovalentBondAst> {
    prop::sample::select(NONCOVALENT_KINDS).prop_map(|kind| NoncovalentBondAst {
        kind: NoncovalentBondKindAst::Lit(kind),
        constraints: Default::default(),
    })
}

/// Generate k distinct atom indices in [0, atom_count).
fn distinct_atoms_strategy(
    atom_count: usize,
    min_k: usize,
    max_k: usize,
) -> BoxedStrategy<Vec<AtomIdx>> {
    if atom_count < min_k {
        return Just(Vec::new()).boxed();
    }
    let max_k = max_k.min(atom_count);
    (min_k..=max_k)
        .prop_flat_map(move |k| {
            prop::collection::vec(0..atom_count as u32, k).prop_map(move |mut v| {
                v.sort_unstable();
                v.dedup();
                // If dedup shrank the vec below k, pad from the start (always valid).
                let mut i = 0u32;
                while v.len() < k && (i as usize) < atom_count {
                    if !v.contains(&i) {
                        v.push(i);
                    }
                    i += 1;
                }
                v.sort_unstable();
                v.into_iter().map(AtomIdx).collect()
            })
        })
        .boxed()
}

fn molecule_ast_strategy() -> impl Strategy<Value = MoleculeAst> {
    (0usize..=5)
        .prop_flat_map(|atom_count| {
            let atoms = prop::collection::vec(atom_ast_strategy(), atom_count);
            let edges = edge_set_strategy(atom_count);
            let bond_data = prop::collection::vec(bond_ast_strategy(), 0..=8);
            (Just(atom_count), atoms, edges, bond_data)
        })
        .prop_flat_map(|(atom_count, atoms, edges, bond_pool)| {
            // Truncate bond pool to the number of edges generated.
            let bond_count = edges.len();
            let bonds: Vec<BondAst> = bond_pool
                .into_iter()
                .chain(repeat_with(|| BondAst::from_order(1)))
                .take(bond_count)
                .collect();
            let bonds_full: Vec<_> = edges
                .iter()
                .zip(bonds)
                .map(|(&[a, b], bond)| (AtomIdx(a), AtomIdx(b), bond))
                .collect();

            let dative_count_max = (atom_count / 2).min(2);
            let aromatic_count_max = (atom_count / 3).min(2);
            let multicenter_count_max = (atom_count / 3).min(2);
            let noncovalent_count_max = (atom_count / 2).min(2);

            let datives = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    dative_bond_strategy(),
                ),
                0..=dative_count_max,
            );
            let aromatics = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3)))
                    .prop_flat_map(|atoms| {
                        let n = atoms.len();
                        (Just(atoms), aromatic_system_ast_for(n))
                    }),
                0..=aromatic_count_max,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3)))
                    .prop_flat_map(|atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_ast_for(n))
                    }),
                0..=multicenter_count_max,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_ast_strategy(),
                ),
                0..=noncovalent_count_max,
            );

            (
                Just(atoms),
                Just(bonds_full),
                datives,
                aromatics,
                multicenters,
                noncovalents,
                Just(atom_count),
            )
        })
        .prop_map(
            |(atoms, bonds, datives, aromatics, multicenters, noncovalents, _n)| {
                let dative_triples: Vec<_> = datives
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((vec![*a], *b, data)),
                        _ => None,
                    })
                    .collect();
                let aromatic_entries: Vec<_> = aromatics
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let multicenter_entries: Vec<_> = multicenters
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let noncovalent_triples: Vec<_> = noncovalents
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((*a, *b, data)),
                        _ => None,
                    })
                    .collect();
                MoleculeAst::from_parts(
                    atoms,
                    bonds,
                    dative_triples,
                    aromatic_entries,
                    multicenter_entries,
                    noncovalent_triples,
                    Constraints::new(),
                )
            },
        )
}

/// Per-entity counts for a generated `MoleculeAst`. Carried into the
/// constraint generator so that ref-bearing variants (`Constraint::Atom`,
/// relational constraints, anchors) only emit valid in-bounds indices.
#[derive(Clone, Copy)]
struct ConstraintCounts {
    atom: usize,
    bond: usize,
    dative: usize,
    aromatic: usize,
    multicenter: usize,
    noncovalent: usize,
}

impl ConstraintCounts {
    fn from_ast(ast: &MoleculeAst) -> Self {
        Self {
            atom: ast.atom_count(),
            bond: ast.bond_count(),
            dative: ast.dative_bond_count(),
            aromatic: ast.aromatic_system_count(),
            multicenter: ast.multicenter_bond_count(),
            noncovalent: ast.noncovalent_bond_count(),
        }
    }
}

fn atom_idx_strategy(atom_count: usize) -> BoxedStrategy<AtomIdx> {
    (0u32..atom_count as u32).prop_map(AtomIdx).boxed()
}

fn bond_idx_strategy(bond_count: usize) -> BoxedStrategy<BondIdx> {
    (0u32..bond_count as u32).prop_map(BondIdx).boxed()
}

fn dative_bond_idx_strategy(count: usize) -> BoxedStrategy<DativeBondIdx> {
    (0u32..count as u32).prop_map(DativeBondIdx).boxed()
}

fn aromatic_system_idx_strategy(count: usize) -> BoxedStrategy<AromaticSystemIdx> {
    (0u32..count as u32).prop_map(AromaticSystemIdx).boxed()
}

fn multicenter_bond_idx_strategy(count: usize) -> BoxedStrategy<MulticenterBondIdx> {
    (0u32..count as u32).prop_map(MulticenterBondIdx).boxed()
}

fn noncovalent_bond_idx_strategy(count: usize) -> BoxedStrategy<NoncovalentBondIdx> {
    (0u32..count as u32).prop_map(NoncovalentBondIdx).boxed()
}

/// Non-recursive constraint leaves: every value-only and relational
/// variant. Combinators wrap these in `constraint_strategy` below.
fn constraint_leaf_strategy(counts: ConstraintCounts) -> BoxedStrategy<Constraint> {
    let mut choices: Vec<BoxedStrategy<Constraint>> = Vec::new();

    if counts.atom > 0 {
        let atom_idx = atom_idx_strategy(counts.atom);

        // Constraint::Atom carrying any AtomConstraint variant.
        let atom_leaf = (atom_idx.clone(), atom_constraint_strategy())
            .prop_map(|(idx, c)| Constraint::Atom(idx, c));
        choices.push(atom_leaf.boxed());

        // MoleculeConstraint variants over atom refs.
        let max_atoms = counts.atom.min(3);
        let atoms_vec = prop::collection::vec(atom_idx.clone(), 1..=max_atoms);
        let optional_atoms = prop::option::of(atoms_vec.clone()).boxed();
        let molecule_connected = optional_atoms
            .clone()
            .prop_map(|atoms| Constraint::Molecule(MoleculeConstraint::Connected { atoms }))
            .boxed();
        // Vacuous molecule-level constraints (Undetermined sum / fully
        // Undetermined spin) elide on render and would break round-trip;
        // restrict the value/spin strategies accordingly.
        let molecule_charge_sum = (optional_atoms.clone(), constraint_inner_value_strategy(-3..=3))
            .prop_map(|(atoms, sum)| {
                Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms, sum })
            })
            .boxed();
        let molecule_spin_sum = (optional_atoms, non_vacuous_spin_state_strategy())
            .prop_map(|(atoms, spin)| {
                Constraint::Molecule(MoleculeConstraint::SpinSum { atoms, spin })
            })
            .boxed();
        choices.push(molecule_connected);
        choices.push(molecule_charge_sum);
        choices.push(molecule_spin_sum);
    }

    if counts.bond > 0 {
        let bond_idx = bond_idx_strategy(counts.bond);
        let bond_leaf = (bond_idx.clone(), bond_constraint_strategy())
            .prop_map(|(idx, c)| Constraint::Bond(idx, c))
            .boxed();
        choices.push(bond_leaf);

        let max_bonds = counts.bond.min(3);
        let optional_bonds =
            prop::option::of(prop::collection::vec(bond_idx, 1..=max_bonds)).boxed();
        let molecule_bond_order_sum = (optional_bonds, constraint_inner_value_strategy(0..=8))
            .prop_map(|(bonds, sum)| {
                Constraint::Molecule(MoleculeConstraint::BondOrderSum { bonds, sum })
            })
            .boxed();
        choices.push(molecule_bond_order_sum);
    }

    if counts.dative > 0 {
        let dative_idx = dative_bond_idx_strategy(counts.dative);
        let dative_leaf = (dative_idx.clone(), dative_bond_constraint_strategy())
            .prop_map(|(idx, c)| Constraint::DativeBond(idx, c))
            .boxed();
        choices.push(dative_leaf);

        if counts.atom > 0 {
            let atom_idx = atom_idx_strategy(counts.atom);
            let donor = (dative_idx.clone(), atom_idx.clone())
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::DativeBondDonor { bond, atom })
                })
                .boxed();
            let acceptor = (dative_idx.clone(), atom_idx.clone())
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAcceptor { bond, atom })
                })
                .boxed();
            let donor_satisfies = (dative_idx.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::DativeBondDonorSatisfies {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            let acceptor_satisfies = (dative_idx.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAcceptorSatisfies {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            choices.push(donor);
            choices.push(acceptor);
            choices.push(donor_satisfies);
            choices.push(acceptor_satisfies);
        }
        if counts.bond > 0 {
            let parallels = (dative_idx, bond_idx_strategy(counts.bond))
                .prop_map(|(dative, parallel)| {
                    Constraint::Relational(RelationalConstraint::DativeBondParallels {
                        dative,
                        parallel,
                    })
                })
                .boxed();
            choices.push(parallels);
        }
    }

    if counts.aromatic > 0 && counts.atom > 0 {
        let system_idx = aromatic_system_idx_strategy(counts.aromatic);
        let atom_idx = atom_idx_strategy(counts.atom);
        let max_atoms = counts.atom.min(3);
        let atoms_vec = prop::collection::vec(atom_idx.clone(), 1..=max_atoms);

        let atoms = (system_idx.clone(), atoms_vec.clone())
            .prop_map(|(system, atoms)| {
                Constraint::Relational(RelationalConstraint::AromaticSystemAtoms { system, atoms })
            })
            .boxed();
        let contains = (system_idx.clone(), atom_idx)
            .prop_map(|(system, atom)| {
                Constraint::Relational(RelationalConstraint::AromaticSystemContains {
                    system,
                    atom,
                })
            })
            .boxed();
        let contains_all = (system_idx.clone(), atoms_vec)
            .prop_map(|(system, atoms)| {
                Constraint::Relational(RelationalConstraint::AromaticSystemContainsAll {
                    system,
                    atoms,
                })
            })
            .boxed();
        let all_atoms = (system_idx.clone(), atom_constraint_strategy())
            .prop_map(|(system, predicate)| {
                Constraint::Relational(RelationalConstraint::AromaticSystemAllAtoms {
                    system,
                    predicate: Box::new(predicate),
                })
            })
            .boxed();
        let any_atom = (system_idx, atom_constraint_strategy())
            .prop_map(|(system, predicate)| {
                Constraint::Relational(RelationalConstraint::AromaticSystemAnyAtom {
                    system,
                    predicate: Box::new(predicate),
                })
            })
            .boxed();
        choices.push(atoms);
        choices.push(contains);
        choices.push(contains_all);
        choices.push(all_atoms);
        choices.push(any_atom);
    }

    if counts.multicenter > 0 && counts.atom > 0 {
        let bond_idx = multicenter_bond_idx_strategy(counts.multicenter);
        let atom_idx = atom_idx_strategy(counts.atom);
        let max_atoms = counts.atom.min(3);
        let atoms_vec = prop::collection::vec(atom_idx.clone(), 1..=max_atoms);

        let atoms = (bond_idx.clone(), atoms_vec.clone())
            .prop_map(|(bond, atoms)| {
                Constraint::Relational(RelationalConstraint::MulticenterBondAtoms { bond, atoms })
            })
            .boxed();
        let contains = (bond_idx.clone(), atom_idx)
            .prop_map(|(bond, atom)| {
                Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond, atom })
            })
            .boxed();
        let contains_all = (bond_idx.clone(), atoms_vec)
            .prop_map(|(bond, atoms)| {
                Constraint::Relational(RelationalConstraint::MulticenterBondContainsAll {
                    bond,
                    atoms,
                })
            })
            .boxed();
        let all_atoms = (bond_idx.clone(), atom_constraint_strategy())
            .prop_map(|(bond, predicate)| {
                Constraint::Relational(RelationalConstraint::MulticenterBondAllAtoms {
                    bond,
                    predicate: Box::new(predicate),
                })
            })
            .boxed();
        let any_atom = (bond_idx, atom_constraint_strategy())
            .prop_map(|(bond, predicate)| {
                Constraint::Relational(RelationalConstraint::MulticenterBondAnyAtom {
                    bond,
                    predicate: Box::new(predicate),
                })
            })
            .boxed();
        choices.push(atoms);
        choices.push(contains);
        choices.push(contains_all);
        choices.push(all_atoms);
        choices.push(any_atom);
    }

    if counts.noncovalent > 0 && counts.atom > 0 {
        let bond_idx = noncovalent_bond_idx_strategy(counts.noncovalent);
        let atom_idx = atom_idx_strategy(counts.atom);

        let ends = (bond_idx.clone(), atom_idx.clone(), atom_idx.clone())
            .prop_map(|(bond, a, b)| {
                Constraint::Relational(RelationalConstraint::NoncovalentBondEnds {
                    bond,
                    atoms: [a, b],
                })
            })
            .boxed();
        let contains = (bond_idx.clone(), atom_idx)
            .prop_map(|(bond, atom)| {
                Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond, atom })
            })
            .boxed();
        let ends_satisfy = (
            bond_idx,
            atom_constraint_strategy(),
            atom_constraint_strategy(),
        )
            .prop_map(|(bond, p1, p2)| {
                Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
                    bond,
                    predicates: [Box::new(p1), Box::new(p2)],
                })
            })
            .boxed();
        choices.push(ends);
        choices.push(contains);
        choices.push(ends_satisfy);
    }

    // SubPattern: pattern molecule and a small anchor pinning the first few
    // entities to themselves on both sides (capped to keep refs valid).
    let target_counts = counts;
    let sub_pattern = molecule_ast_strategy()
        .prop_flat_map(move |pattern| {
            let pattern_counts = ConstraintCounts::from_ast(&pattern);
            (
                Just(pattern),
                sub_pattern_anchor_strategy(target_counts, pattern_counts),
            )
        })
        .prop_map(|(pattern, anchor)| {
            Constraint::Molecule(MoleculeConstraint::SubPattern {
                anchor,
                pattern: Box::new(pattern),
            })
        })
        .boxed();
    choices.push(sub_pattern);

    if choices.is_empty() {
        return Just(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: None,
        }))
        .boxed();
    }
    prop::strategy::Union::new(choices).boxed()
}

/// Sub-pattern anchor: link the first few entities of each kind pairwise,
/// capped at the minimum of the two molecules' counts on each side so all
/// refs are valid in their respective metadata scopes.
fn sub_pattern_anchor_strategy(
    target: ConstraintCounts,
    pattern: ConstraintCounts,
) -> BoxedStrategy<SubPatternAnchor> {
    let atom_pairs = target.atom.min(pattern.atom).min(2);
    let bond_pairs = target.bond.min(pattern.bond).min(2);
    let dative_pairs = target.dative.min(pattern.dative).min(1);
    let aromatic_pairs = target.aromatic.min(pattern.aromatic).min(1);
    let multicenter_pairs = target.multicenter.min(pattern.multicenter).min(1);
    let noncovalent_pairs = target.noncovalent.min(pattern.noncovalent).min(1);
    (
        0..=atom_pairs,
        0..=bond_pairs,
        0..=dative_pairs,
        0..=aromatic_pairs,
        0..=multicenter_pairs,
        0..=noncovalent_pairs,
    )
        .prop_map(|(a, b, d, ar, mc, nc)| {
            let mut anchor = SubPatternAnchor::new();
            for i in 0..a {
                anchor.push_atom(AtomIdx(i as u32), AtomIdx(i as u32));
            }
            for i in 0..b {
                anchor.push_bond(BondIdx(i as u32), BondIdx(i as u32));
            }
            for i in 0..d {
                anchor.push_dative_bond(DativeBondIdx(i as u32), DativeBondIdx(i as u32));
            }
            for i in 0..ar {
                anchor
                    .push_aromatic_system(AromaticSystemIdx(i as u32), AromaticSystemIdx(i as u32));
            }
            for i in 0..mc {
                anchor.push_multicenter_bond(
                    MulticenterBondIdx(i as u32),
                    MulticenterBondIdx(i as u32),
                );
            }
            for i in 0..nc {
                anchor.push_noncovalent_bond(
                    NoncovalentBondIdx(i as u32),
                    NoncovalentBondIdx(i as u32),
                );
            }
            anchor
        })
        .boxed()
}

/// Constraint tree: leaves wrapped in bounded-depth combinators (And/Or/Not).
fn constraint_strategy(counts: ConstraintCounts) -> BoxedStrategy<Constraint> {
    constraint_leaf_strategy(counts)
        .prop_recursive(2, 6, 3, |inner| {
            prop_oneof![
                prop::collection::vec(inner.clone(), 1..=3).prop_map(Constraint::And),
                prop::collection::vec(inner.clone(), 1..=3).prop_map(Constraint::Or),
                inner.prop_map(|c| Constraint::Not(Box::new(c))),
            ]
            .boxed()
        })
        .boxed()
}

fn molecule_ast_with_constraints_strategy() -> impl Strategy<Value = MoleculeAst> {
    molecule_ast_strategy().prop_flat_map(|ast| {
        let counts = ConstraintCounts::from_ast(&ast);
        let max_constraints = 4usize;
        (
            Just(ast),
            prop::collection::vec(constraint_strategy(counts), 0..=max_constraints),
        )
            .prop_map(|(ast, constraints)| {
                let mut cs = Constraints::new();
                for c in constraints {
                    cs.push(c);
                }
                let mut b = ast.edit();
                *b.constraints_mut() = cs;
                b.build()
            })
    })
}

fn molecule_dsl_strategy() -> impl Strategy<Value = MoleculeDsl> {
    molecule_ast_with_constraints_strategy()
        .prop_map(|ast| MoleculeDsl::from_parts(ast, Metadata::default()))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn test_molecule_dsl_to_edn_from_edn_tree_roundtrip(dsl in molecule_dsl_strategy()) {
        let edn = dsl.to_edn();
        let parsed = MoleculeDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_to_edn_from_edn_str_roundtrip(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_streaming_matches_tree(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let via_stream = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming: {e}")))?;
        let tree = read_string(&rendered)
            .map_err(|e| TestCaseError::fail(format!("edn parse: {e}")))?;
        let via_tree = MoleculeDsl::from_edn(&tree)
            .map_err(|e| TestCaseError::fail(format!("tree: {e}")))?;
        prop_assert_eq!(via_stream, via_tree);
    }

    /// Direct `MoleculeAst::ToEdn` / `FromEdn` round-trips are the identity.
    /// Refs render as positional integers (no id keywords); the AST carries
    /// no metadata, so canonical EDN parses back to an equal AST.
    #[test]
    fn test_molecule_ast_to_edn_from_edn_tree_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let edn = ast.to_edn();
        let parsed = MoleculeAst::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(ast, parsed);
    }

    #[test]
    fn test_molecule_ast_to_edn_from_edn_str_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let rendered = ast.to_edn().to_string();
        let parsed = MoleculeAst::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(ast, parsed);
    }

    #[test]
    fn test_molecule_dsl_double_render_is_stable(dsl in molecule_dsl_strategy()) {
        let s1 = dsl.to_edn().to_string();
        let d1 = MoleculeDsl::from_edn_str(&s1)
            .map_err(|e| TestCaseError::fail(format!("first parse: {e}")))?;
        let s2 = d1.to_edn().to_string();
        prop_assert_eq!(s1, s2);
    }

    /// `lift_constraints` followed by `inline_constraints` is idempotent:
    /// running the pair twice yields the same `MoleculeAst` as running it
    /// once. This holds even if the original AST has duplicate (entity, kind)
    /// entries across the inline + molecule scopes — the first pass collapses
    /// them via the entity store's last-wins policy and the second pass is
    /// a fixpoint.
    #[test]
    fn test_lift_inline_idempotent(ast in molecule_ast_with_constraints_strategy()) {
        let mut once = ast.clone();
        once.lift_constraints();
        once.inline_constraints();

        let mut twice = once.clone();
        twice.lift_constraints();
        twice.inline_constraints();

        prop_assert_eq!(once, twice);
    }

    /// `lift_constraints` drains every entity's inline `constraints` store.
    #[test]
    fn test_lift_drains_entity_stores(ast in molecule_ast_with_constraints_strategy()) {
        let mut a = ast;
        a.lift_constraints();
        for view in a.atoms().iter() {
            prop_assert!(view.data.constraints.is_empty());
        }
        for view in a.bonds().iter() {
            prop_assert!(view.data.constraints.is_empty());
        }
        for view in a.dative_bonds().iter() {
            prop_assert!(view.data.constraints.is_empty());
        }
    }

    /// `inline_constraints` removes every TOP-LEVEL inline-capable narrow
    /// leaf from the molecule list. Combinator-nested entries, relational
    /// leaves, molecule-scope leaves are preserved.
    #[test]
    fn test_inline_removes_top_level_leaves(ast in molecule_ast_with_constraints_strategy()) {
        let mut a = ast;
        a.inline_constraints();
        for c in a.constraints().iter() {
            prop_assert!(
                !matches!(
                    c,
                    Constraint::Atom(..)
                        | Constraint::Bond(..)
                        | Constraint::DativeBond(..)
                        | Constraint::AromaticSystem(..)
                        | Constraint::MulticenterBond(..)
                        | Constraint::NoncovalentBond(..)
                ),
                "inline-capable narrow leaf survived inline_constraints: {c:?}",
            );
        }
    }

    // -- Per-entity Display ↔ FromStr roundtrip ----------------------------
    //
    // The entity DSL types carry a compact string form (Display) parsed by
    // their own `FromStr`. The invariant is `parse(display(x)) == x` for any
    // generator-produced AST. Exercises the same render/parse pairing at the
    // entity layer that the molecule-level EDN tests cover for the full map.

    #[test]
    fn test_atom_dsl_display_from_str_roundtrip(atom in atom_ast_strategy()) {
        let dsl = AtomDsl(atom);
        let rendered = dsl.to_string();
        let parsed: AtomDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_bond_dsl_display_from_str_roundtrip(bond in bond_ast_strategy()) {
        let dsl = BondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: BondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// `BondDsl::ToEdn` ↔ `FromEdn` round-trips for any generated bond
    /// shape. Non-canonical bonds render as bond strings; canonical
    /// shapes (order-only, no charge / spin / non-aromatic constraints,
    /// or order-1 with the `Aromatic` flag) render as keyword shorthands.
    #[test]
    fn test_bond_dsl_to_edn_from_edn_roundtrip(bond in bond_ast_strategy()) {
        let dsl = BondDsl(bond);
        let edn = dsl.to_edn();
        let parsed = BondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical-shape bonds render as keyword shorthands, and the
    /// keyword form parses back to the same AST.
    #[test]
    fn test_bond_dsl_keyword_to_edn_from_edn_roundtrip(
        bond in canonical_keyword_bond_strategy(),
    ) {
        let dsl = BondDsl(bond);
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical bond, got {edn:?}",
        );
        let parsed = BondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_aromatic_system_dsl_display_from_str_roundtrip(
        system in aromatic_system_ast_strategy(),
    ) {
        let dsl = AromaticSystemDsl(system);
        let rendered = dsl.to_string();
        let parsed: AromaticSystemDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_multicenter_bond_dsl_display_from_str_roundtrip(
        bond in multicenter_bond_ast_strategy(),
    ) {
        let dsl = MulticenterBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: MulticenterBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_dative_bond_dsl_display_from_str_roundtrip(
        bond in dative_bond_strategy(),
    ) {
        let dsl = DativeBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: DativeBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_noncovalent_bond_dsl_display_from_str_roundtrip(
        bond in noncovalent_bond_ast_strategy(),
    ) {
        let dsl = NoncovalentBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: NoncovalentBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    // region: ValueAst::simplify

    /// `simplify` is idempotent: `x.simplify().simplify() == x.simplify()`.
    #[test]
    fn test_value_ast_simplify_idempotent(v in any_value_ast_strategy()) {
        let once = v.simplify();
        let twice = once.clone().simplify();
        prop_assert_eq!(once, twice);
    }

    /// `simplify()` is the canonical form: for any generated `ValueAst`,
    /// rendering and parsing yields a value that — once simplified —
    /// equals `simplify()` on the original. The parser produces a partly
    /// canonical form (it folds within `Expr` but doesn't always lift
    /// `Expr(Lit(n))` to `ValueAst::Lit(n)`); simplify completes the
    /// canonicalization on both sides.
    #[test]
    fn test_value_ast_render_parse_equals_simplify(v in any_value_ast_strategy()) {
        let dsl = ValueDsl(v.clone());
        let rendered = dsl.to_string();
        let parsed = parse_value(&rendered).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered:?}"))
        })?;
        prop_assert_eq!(parsed.simplify(), v.simplify());
    }

    // endregion: ValueAst::simplify
}
