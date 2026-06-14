//! Shared proptest generators for the umol-ast property suite. Domain imports
//! are re-exported (`pub(crate) use`) so the per-area test modules need only
//! `use proptest::prelude::*; use crate::strategies::*;`.

pub(crate) use std::collections::HashSet;
pub(crate) use std::fmt::Debug;
pub(crate) use std::iter::repeat_with;
pub(crate) use std::ops::RangeInclusive;

use proptest::prelude::*;
pub(crate) use umol_ast::ast::{
    AddBond, ArithOp, AromaticSystemAst, AromaticSystemConstraint, AromaticSystemConstraintKind,
    AromaticSystemConstraints, AromaticSystemId, AromaticValenceAst, AtomAst, AtomConstraint,
    AtomConstraintKind, AtomConstraints, AtomFieldChange, AtomId, AtomRef, BondAst, BondConstraint,
    BondConstraintKind, BondConstraints, BondFieldChange, BondId, BondRef, Constraint, Constraints,
    DativeBondAst, DativeBondConstraint, DativeBondConstraintKind, DativeBondConstraints,
    DativeBondId, Edit, ElementAst, FluxionalityAst, IsotopeMassAst, Lattice, LigandPairAst,
    LigandSymmetryAst, MemOp, MoleculeAst, MoleculeConstraint, MulticenterBondAst,
    MulticenterBondConstraint, MulticenterBondConstraintKind, MulticenterBondConstraints,
    MulticenterBondId, MulticenterValenceAst, NoncovalentBondAst, NoncovalentBondId,
    NoncovalentBondKind, NoncovalentBondKindAst, OrientedPermutationAst, PermutationAst, RelOp,
    RelationalConstraint, SpinStateAst, StereoAtomAst, StereoAtomConstraint, StereoAtomId,
    StereoBondAst, StereoBondConstraint, StereoBondId, StereoConfigurationAst, StereoCosetAst,
    StereoExpr, StereoKind, StereoLigand, StereoLigandId, StereoLigandKind, Stereogenicity,
    StereogenicityAst,
    StereogenicityRelationAst, SubPatternAnchor, Topicity, TopicityAst, TopicityRelationAst,
    ValueAst, ValueExpr,
};
pub(crate) use umol_ast::dsl::{
    parse_value, AromaticSystemDsl, AtomDsl, BondDsl, DativeBondDsl, Metadata, MoleculeDsl,
    MulticenterBondDsl, NoncovalentBondDsl, StereoAtomConstraintDsl, StereoAtomDsl,
    StereoBondConstraintDsl, StereoBondDsl, ValueDsl,
};
pub(crate) use umol_edn::{read_string, Edn, FromEdn, ToEdn};
pub(crate) use umol_perm::{Orientation, Permutation};
pub(crate) use umol_shared::element::Element;

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

pub(crate) fn element_strategy() -> impl Strategy<Value = Element> {
    prop::sample::select(ELEMENTS)
}

pub(crate) fn id_strategy() -> impl Strategy<Value = String> {
    "[a-zA-Z][a-zA-Z0-9_]{0,3}".prop_map(|s| s.to_string())
}

pub(crate) fn element_ast_strategy() -> impl Strategy<Value = ElementAst> {
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
                ElementAst::bind(id, set, MemOp::In)
            }),
        1 => id_strategy().prop_map(ElementAst::Ref),
    ]
}

pub(crate) fn value_basic(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        4 => Just(ValueAst::Undetermined),
        4 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range.clone(), 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::set(v)
        }),
        1 => id_strategy().prop_map(ValueAst::reference),
        1 => (id_strategy(), prop::collection::vec(range, 1..=3)).prop_map(|(id, mut v)| {
            v.sort_unstable();
            v.dedup();
            ValueAst::bind(id, v)
        }),
        1 => top_expr_strategy().prop_map(ValueAst::expr),
    ]
}

/// Arithmetic-typed ValueExpr: produces only the arithmetic subset of `ValueExpr`
/// (`Lit`, `Var`, `Neg(arith)`, `BinOp(arith, op, arith)`). Includes
/// negative `Lit` and `Neg(Neg(_))` shapes that the parser canonicalizes,
/// to be paired with `simplify()` for roundtrip testing.
pub(crate) fn arith_expr_strategy() -> BoxedStrategy<ValueExpr> {
    let leaf = prop_oneof![
        (-10i64..=10).prop_map(ValueExpr::Lit),
        id_strategy().prop_map(ValueExpr::Var),
    ]
    .boxed();
    leaf.prop_recursive(3, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|e| ValueExpr::Neg(Box::new(e))),
            (inner.clone(), arith_op_strategy(), inner).prop_map(|(l, op, r)| ValueExpr::BinOp(
                Box::new(l),
                op,
                Box::new(r)
            )),
        ]
        .boxed()
    })
    .boxed()
}

/// Boolean-typed ValueExpr: `Rel(arith, op, arith)`, `Mem(arith, set)`,
/// `Not(bool)`, `And(bool*)`, `Or(bool*)`. Each boolean recursion correctly
/// roots in arithmetic leaves so the parser accepts the rendered form.
pub(crate) fn bool_expr_strategy() -> BoxedStrategy<ValueExpr> {
    let arith = arith_expr_strategy();
    let leaf =
        prop_oneof![
            (arith.clone(), rel_op_strategy(), arith.clone())
                .prop_map(|(l, op, r)| ValueExpr::Rel(Box::new(l), op, Box::new(r))),
            (arith, prop::collection::vec(-10i64..=10, 1..=3))
                .prop_map(|(e, set)| ValueExpr::Mem(Box::new(e), set)),
        ]
        .boxed();
    leaf.prop_recursive(2, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|e| ValueExpr::Not(Box::new(e))),
            prop::collection::vec(inner.clone(), 1..=3).prop_map(ValueExpr::And),
            prop::collection::vec(inner, 1..=3).prop_map(ValueExpr::Or),
        ]
        .boxed()
    })
    .boxed()
}

pub(crate) fn any_expr_strategy() -> BoxedStrategy<ValueExpr> {
    prop_oneof![arith_expr_strategy(), bool_expr_strategy()].boxed()
}

pub(crate) fn any_value_ast_strategy() -> BoxedStrategy<ValueAst> {
    prop_oneof![
        Just(ValueAst::Undetermined),
        (-10i64..=10).prop_map(ValueAst::Lit),
        prop::collection::vec(-10i64..=10, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::set(v)
        }),
        id_strategy().prop_map(ValueAst::reference),
        (id_strategy(), prop::collection::vec(-10i64..=10, 1..=3)).prop_map(|(id, mut v)| {
            v.sort_unstable();
            v.dedup();
            ValueAst::bind(id, v)
        }),
        any_expr_strategy().prop_map(ValueAst::expr),
    ]
    .boxed()
}

pub(crate) fn arith_op_strategy() -> impl Strategy<Value = ArithOp> {
    prop_oneof![
        Just(ArithOp::Add),
        Just(ArithOp::Sub),
        Just(ArithOp::Mul),
        Just(ArithOp::Div),
        Just(ArithOp::Rem),
    ]
}

pub(crate) fn rel_op_strategy() -> impl Strategy<Value = RelOp> {
    prop_oneof![
        Just(RelOp::Le),
        Just(RelOp::Ge),
        Just(RelOp::Eq),
        Just(RelOp::Lt),
        Just(RelOp::Gt),
    ]
}

/// ValueExpr leaf: non-negative `Lit` or `Var`. Safe as a subexpression of any
/// operator; **not** safe as the outermost `ValueAst::ValueExpr` wrapper (see
/// `top_expr_strategy`). Negative literals are excluded because the ValueExpr
/// grammar has no `Lit(-n)` parse — `-n` inside an ValueExpr always parses as
/// `Neg(Lit(n))`, so emitting `ValueExpr::Lit(-n)` from the generator would fail
/// the structural roundtrip equality even though the two forms are
/// semantically identical under `is_ground` / `evaluate`. Negative integers
/// still appear elsewhere (top-level `ValueAst::Lit`, `Set`, and
/// `ValueExpr::Mem` sets all route through `dec_int` which reads signed).
pub(crate) fn expr_leaf_strategy() -> impl Strategy<Value = ValueExpr> {
    prop_oneof![
        (0i64..=10).prop_map(ValueExpr::Lit),
        id_strategy().prop_map(ValueExpr::Var),
    ]
}

/// ValueExpr tree intended as the outermost `ValueAst::ValueExpr(e)`. Constraints on
/// the outermost shape — inner compositions (BinOp / Rel / Neg / And / Or)
/// may freely contain `Var` leaves, so `?h + 1` and `(?h <= 1) | (?h >= 2)`
/// generate:
///
/// - Must not render to a pure integer literal with optional sign — the
///   `value` parser's `dec_int` alt would match first (bare `Lit(n)` → `n`
///   → `ValueAst::Lit(n)`, `Neg(Lit(n))` → `-n` → `ValueAst::Lit(-n)`).
/// - Must not be a bare `Var(_)` — renders to `?id`, which the parser now
///   intercepts as `ValueAst::Ref` (separate generator arm covers Ref).
/// - Must not be `Mem(Var(_), _)` — renders to `?id :: {set}`, which the
///   parser now intercepts as `ValueAst::Bind` (separate generator arm
///   covers Bind).
/// - Avoids `Neg(Neg(_))` anywhere (the parser folds consecutive signs).
/// - `Or` / `And` children are leaves so the parser can't flatten
///   consecutive same-op tokens.
///
/// Returns a boxed strategy to keep the composed type size bounded when
/// plugged into every value field across the molecule tree.
pub(crate) fn top_expr_strategy() -> BoxedStrategy<ValueExpr> {
    let set = prop::collection::vec(-10i64..=10, 1..=3);
    let non_var_leaf = (0i64..=10).prop_map(ValueExpr::Lit);
    prop_oneof![
        (
            expr_leaf_strategy(),
            arith_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| ValueExpr::BinOp(Box::new(a), op, Box::new(b))),
        (
            expr_leaf_strategy(),
            rel_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| ValueExpr::Rel(Box::new(a), op, Box::new(b))),
        (non_var_leaf, set).prop_map(|(e, s)| ValueExpr::Mem(Box::new(e), s)),
        (
            expr_leaf_strategy(),
            rel_op_strategy(),
            expr_leaf_strategy()
        )
            .prop_map(|(a, op, b)| {
                ValueExpr::Not(Box::new(ValueExpr::Rel(Box::new(a), op, Box::new(b))))
            }),
        // `Neg` of `Var` renders `-?id`; safe (non-Lit inner, no sign folding).
        id_strategy().prop_map(|x| ValueExpr::Neg(Box::new(ValueExpr::Var(x)))),
        // `Or` / `And` with exactly leaf children so the parser can't flatten.
        prop::collection::vec(expr_leaf_strategy(), 2..=3).prop_map(ValueExpr::Or),
        prop::collection::vec(expr_leaf_strategy(), 2..=3).prop_map(ValueExpr::And),
    ]
    .boxed()
}

pub(crate) fn isotope_strategy() -> impl Strategy<Value = IsotopeMassAst> {
    prop_oneof![
        3 => Just(IsotopeMassAst::Natural),
        3 => Just(IsotopeMassAst::Undetermined),
        3 => (1i64..=250).prop_map(IsotopeMassAst::Lit),
        1 => prop::collection::vec(1i64..=250, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            IsotopeMassAst::set(v)
        }),
        1 => (1i64..=250).prop_map(IsotopeMassAst::Not),
        1 => prop::collection::vec(1i64..=250, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            IsotopeMassAst::not_set(v)
        }),
        1 => id_strategy().prop_map(IsotopeMassAst::reference),
        1 => (id_strategy(), prop::collection::vec(1i64..=250, 1..=3), prop_oneof![
            Just(MemOp::In),
            Just(MemOp::NotIn),
        ]).prop_map(|(id, mut v, polarity)| {
            v.sort_unstable();
            v.dedup();
            IsotopeMassAst::bind(id, v, polarity)
        }),
    ]
}

pub(crate) fn spin_state_strategy() -> impl Strategy<Value = SpinStateAst> {
    // DSL preserves spin fields field-wise. Physical (u, m) parity is
    // not a parse-time check, so any independent pair must roundtrip.
    (value_basic(0..=6), value_basic(1..=7)).prop_map(|(u, m)| SpinStateAst {
        unpaired: u,
        multiplicity: m,
    })
}

/// `SpinStateAst` with at least one of `unpaired` / `multiplicity` not
/// `Undetermined`. Used inside `MoleculeConstraint::SpinSum` and similar
/// where a fully-vacuous spin state would elide on render.
pub(crate) fn non_vacuous_spin_state_strategy() -> impl Strategy<Value = SpinStateAst> {
    (value_basic(0..=6), value_basic(1..=7))
        .prop_map(|(u, m)| SpinStateAst {
            unpaired: u,
            multiplicity: m,
        })
        .prop_filter("non-vacuous spin", |s| !s.is_undetermined())
}

/// Simple value strategy used inside constraint values: `Undetermined`,
/// `Lit`, and `Set`. No `ValueExpr` — the constraint formatters route to
/// `fmt_value_field_required` / `fmt_ring_count` / the various `#r` blocks,
/// and `ValueExpr(Lit(n))` or `ValueExpr(Neg(Lit(n)))` would render to a pure integer
/// that the parser then re-reads as a plain `Lit`, breaking roundtrip. The
/// molecule-level EDN tests cover `ValueExpr` on constraint values through the
/// tree-based path, so the gap is contained.
pub(crate) fn constraint_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => Just(ValueAst::Undetermined),
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::set(v)
        }),
    ]
}

/// `Lit`/`Set` only — still used by the ring-size strategies where
/// `Undetermined` on the inner value collapses into a dropped constraint
/// in the entity-level formatter (see `BondConstraint::RingSize` /
/// `DativeBondConstraint::RingSize` — vacuous, intentionally dropped).
pub(crate) fn constraint_inner_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        range.clone().prop_map(ValueAst::Lit),
        prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::set(v)
        }),
    ]
}

/// `AromaticValenceAst::Undetermined` is vacuous (renders empty per the
/// canonical-rendering rule). Excluded so the strategy stays inside the
/// render → reparse identity.
pub(crate) fn aromatic_valence_ast_strategy() -> impl Strategy<Value = AromaticValenceAst> {
    prop_oneof![
        Just(AromaticValenceAst::NotAromatic),
        constraint_value_strategy(0..=6).prop_map(AromaticValenceAst::Aromatic),
    ]
}

pub(crate) fn multicenter_valence_ast_strategy() -> impl Strategy<Value = MulticenterValenceAst> {
    prop_oneof![
        Just(MulticenterValenceAst::NotMulticenter),
        constraint_value_strategy(0..=6).prop_map(MulticenterValenceAst::Multicenter),
    ]
}

/// Atom constraints route through `fmt_value_field_required` (or
/// `fmt_ring_count` for `#R`), which elide vacuous (Undetermined) payloads
/// per the canonical-rendering rule. Generators excluding `Undetermined`
/// keep the render → reparse identity intact.
pub(crate) fn atom_constraint_strategy() -> BoxedStrategy<AtomConstraint> {
    prop_oneof![
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::Valence),
        constraint_inner_value_strategy(0..=8).prop_map(AtomConstraint::TotalValence),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraint::DonatedPairs),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraint::AcceptedPairs),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::Degree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::TotalDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::RingDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::RingValence),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::TotalHydrogens),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(AtomConstraint::RingSize),
        aromatic_valence_ast_strategy().prop_map(AtomConstraint::AromaticValence),
        multicenter_valence_ast_strategy().prop_map(AtomConstraint::MulticenterValence),
        stereo_config_strategy().prop_map(AtomConstraint::TetrahedralStereo),
    ]
    .boxed()
}

pub(crate) fn atom_constraints_strategy() -> impl Strategy<Value = AtomConstraints> {
    prop::collection::vec(atom_constraint_strategy(), 0..=3).prop_map(|list| {
        let mut cs = AtomConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

pub(crate) fn bond_constraint_strategy() -> BoxedStrategy<BondConstraint> {
    prop_oneof![
        Just(BondConstraint::Aromatic),
        constraint_inner_value_strategy(0..=6).prop_map(BondConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(BondConstraint::RingSize),
        stereo_config_strategy().prop_map(BondConstraint::CisTransStereo),
    ]
    .boxed()
}

pub(crate) fn bond_constraints_strategy() -> impl Strategy<Value = BondConstraints> {
    prop::collection::vec(bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = BondConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

pub(crate) fn dative_bond_constraint_strategy() -> BoxedStrategy<DativeBondConstraint> {
    prop_oneof![
        Just(DativeBondConstraint::Aromatic),
        constraint_inner_value_strategy(0..=6).prop_map(DativeBondConstraint::RingCount),
        constraint_inner_value_strategy(3..=10).prop_map(DativeBondConstraint::RingSize),
    ]
    .boxed()
}

pub(crate) fn dative_bond_constraints_strategy() -> impl Strategy<Value = DativeBondConstraints> {
    prop::collection::vec(dative_bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = DativeBondConstraints::new();
        for c in list {
            cs.add(c);
        }
        cs
    })
}

prop_compose! {
    pub(crate) fn atom_ast_strategy()
    (
        element in element_ast_strategy(),
        isotope in isotope_strategy(),
        charge in value_basic(-2..=2),
        implicit_hydrogens in value_basic(0..=4),
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
    pub(crate) fn bond_ast_strategy()
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
pub(crate) fn canonical_keyword_bond_strategy() -> impl Strategy<Value = BondAst> {
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
pub(crate) fn edge_set_strategy(atom_count: usize) -> impl Strategy<Value = Vec<[u32; 2]>> {
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

pub(crate) fn dative_bond_strategy() -> impl Strategy<Value = DativeBondAst> {
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
    (order_strategy, dative_bond_constraints_strategy()).prop_map(|(order, constraints)| {
        DativeBondAst {
            acceptor_slot: 0,
            order,
            constraints,
        }
    })
}

/// Optional `ElectronCount` constraint (the asserted total). The strategy
/// emits `None` half the time, otherwise wraps a `ValueAst::Lit` or
/// `Set`. `Undetermined` is excluded because it has no canonical
/// surface form in the entity-string `#e<n>` slot — `#e*` is admitted on
/// parse but the renderer omits the predicate entirely, breaking
/// roundtrip.
pub(crate) fn optional_aromatic_electron_count() -> impl Strategy<Value = AromaticSystemConstraints> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=12)).prop_map(|opt| {
        let mut cs = AromaticSystemConstraints::new();
        if let Some(v) = opt {
            cs.add(AromaticSystemConstraint::ElectronCount(v));
        }
        cs
    })
}

pub(crate) fn optional_multicenter_electron_count() -> impl Strategy<Value = MulticenterBondConstraints> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=8)).prop_map(|opt| {
        let mut cs = MulticenterBondConstraints::new();
        if let Some(v) = opt {
            cs.add(MulticenterBondConstraint::ElectronCount(v));
        }
        cs
    })
}

pub(crate) fn electron_count_value_strategy(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::set(v)
        }),
    ]
}

/// Stand-alone strategy for entity-string roundtrip tests. The per-atom
/// `electrons` vec is empty because the entity string carries no per-atom
/// data; the `ElectronCount` constraint is exercised here via `#e<n>`.
pub(crate) fn aromatic_system_ast_strategy() -> impl Strategy<Value = AromaticSystemAst> {
    (value_basic(-2..=2), optional_aromatic_electron_count()).prop_map(|(charge, constraints)| {
        AromaticSystemAst {
            electrons: Vec::new(),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        }
    })
}

/// Atom-count-aware variant: generates an `AromaticSystemAst` whose
/// `electrons` vec has exactly `atom_count` entries. Includes an optional
/// `ElectronCount` constraint so the molecule-level prop tests exercise
/// both the per-atom vec and the asserted total in the same pass.
pub(crate) fn aromatic_system_ast_for(atom_count: usize) -> impl Strategy<Value = AromaticSystemAst> {
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

pub(crate) fn multicenter_bond_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (value_basic(-2..=2), optional_multicenter_electron_count()).prop_map(
        |(charge, constraints)| MulticenterBondAst {
            electrons: Vec::new(),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        },
    )
}

pub(crate) fn multicenter_bond_ast_for(atom_count: usize) -> impl Strategy<Value = MulticenterBondAst> {
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

pub(crate) fn noncovalent_bond_ast_strategy() -> impl Strategy<Value = NoncovalentBondAst> {
    prop::sample::select(NONCOVALENT_KINDS).prop_map(|kind| NoncovalentBondAst {
        kind: NoncovalentBondKindAst::Lit(kind),
        constraints: Default::default(),
    })
}

/// Coset forms that round-trip through both the entity `:type` string and the
/// EDN coset-form: `Undetermined` (`*`), `Lit`, and a literal set
/// (`{a,b,…}` ↔ EDN vector). The `~`/`^`/`?var` operator-exprs are reserved
/// (§7.14) and excluded.
pub(crate) fn stereo_coset_strategy() -> impl Strategy<Value = StereoCosetAst> {
    prop_oneof![
        Just(StereoCosetAst::Undetermined),
        (0u32..=6).prop_map(StereoCosetAst::Lit),
        prop::collection::vec(0u32..=6, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(v)))
        }),
    ]
}

pub(crate) fn stereo_ligand_kind_strategy() -> impl Strategy<Value = StereoLigandKind> {
    prop_oneof![
        Just(StereoLigandKind::Atom),
        Just(StereoLigandKind::ImplicitHydrogen),
        Just(StereoLigandKind::LonePair),
    ]
}

/// `StereoConfigurationAst` for `#T` / `#C` constraints, excluding the vacuous
/// `Undetermined` config (it renders empty per the canonical-rendering rule,
/// breaking render → reparse — mirrors `aromatic_valence_ast_strategy`).
/// `Stereo(Undetermined)` (the `+` form) is non-vacuous and kept.
pub(crate) fn stereo_config_strategy() -> impl Strategy<Value = StereoConfigurationAst> {
    prop_oneof![
        Just(StereoConfigurationAst::NotStereo),
        stereo_coset_strategy().prop_map(StereoConfigurationAst::Stereo),
    ]
}

pub(crate) fn stereo_atom_kind_strategy() -> impl Strategy<Value = StereoKind> {
    prop_oneof![
        Just(StereoKind::Tetrahedral),
        Just(StereoKind::SquarePlanar),
        Just(StereoKind::TrigonalBipyramidal),
        Just(StereoKind::Octahedral),
    ]
}

/// A permutation of the kind's `degree` positions, as a shuffled one-line image.
pub(crate) fn permutation_strategy(degree: usize) -> impl Strategy<Value = Permutation> {
    Just((0..degree as u8).collect::<Vec<u8>>())
        .prop_shuffle()
        .prop_map(move |image| Permutation::from_image(degree, &image))
}

pub(crate) fn orientation_strategy() -> impl Strategy<Value = Orientation> {
    prop_oneof![Just(Orientation::Proper), Just(Orientation::Improper)]
}

pub(crate) fn mem_op_strategy() -> impl Strategy<Value = MemOp> {
    prop_oneof![Just(MemOp::In), Just(MemOp::NotIn)]
}

pub(crate) fn ligand_pair_strategy(degree: usize) -> impl Strategy<Value = LigandPairAst> {
    (0..degree as u8, 0..degree as u8)
        .prop_map(|(a, b)| LigandPairAst::new(StereoLigandId(a), StereoLigandId(b)))
}

/// Non-vacuous topicity relations only (`Undetermined` elides on render, so it
/// would break the render → reparse roundtrip — mirrors `stereo_config_strategy`).
pub(crate) fn topicity_relation_strategy() -> impl Strategy<Value = TopicityRelationAst> {
    prop_oneof![
        Just(TopicityRelationAst::Lit(Topicity::Homotopic)),
        Just(TopicityRelationAst::Lit(Topicity::Enantiotopic)),
        Just(TopicityRelationAst::Lit(Topicity::Diastereotopic)),
        Just(TopicityRelationAst::NotSet(vec![Topicity::Homotopic])),
        Just(TopicityRelationAst::NotSet(vec![Topicity::Enantiotopic])),
        Just(TopicityRelationAst::NotSet(vec![Topicity::Diastereotopic])),
    ]
}

pub(crate) fn stereogenicity_relation_strategy() -> impl Strategy<Value = StereogenicityRelationAst> {
    prop_oneof![
        Just(StereogenicityRelationAst::Lit(Stereogenicity::Symmetric)),
        Just(StereogenicityRelationAst::Lit(Stereogenicity::Prochiral)),
        Just(StereogenicityRelationAst::Lit(Stereogenicity::Stereogenic)),
        Just(StereogenicityRelationAst::NotSet(vec![Stereogenicity::Symmetric])),
        Just(StereogenicityRelationAst::NotSet(vec![Stereogenicity::Prochiral])),
        Just(StereogenicityRelationAst::NotSet(vec![Stereogenicity::Stereogenic])),
    ]
}

/// Topicity relations spanning the full lattice: the non-vacuous singletons /
/// complements plus `Undetermined` (top).
pub(crate) fn topicity_relation_lattice_strategy() -> impl Strategy<Value = TopicityRelationAst> {
    prop_oneof![
        Just(TopicityRelationAst::Undetermined),
        topicity_relation_strategy(),
    ]
}

pub(crate) fn stereogenicity_relation_lattice_strategy() -> impl Strategy<Value = StereogenicityRelationAst> {
    prop_oneof![
        Just(StereogenicityRelationAst::Undetermined),
        stereogenicity_relation_strategy(),
    ]
}

pub(crate) fn ligand_symmetry_strategy(degree: usize) -> impl Strategy<Value = LigandSymmetryAst> {
    (
        permutation_strategy(degree),
        orientation_strategy(),
        mem_op_strategy(),
    )
        .prop_map(|(perm, orientation, mem)| LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(perm),
                orientation,
            },
            mem,
        })
}

/// Lattice laws for a stereo relation: meet / join commutativity and
/// associativity, absorption, idempotence, and `matches(t) ⇔ meet(t) == Some(t)`.
pub(crate) fn assert_relation_lattice_laws<L: Lattice + Debug>(
    a: &L,
    b: &L,
    c: &L,
) -> Result<(), TestCaseError> {
    prop_assert_eq!(a.meet(b), b.meet(a));
    prop_assert_eq!(a.join(b), b.join(a));
    prop_assert_eq!(
        a.meet(b).and_then(|ab| ab.meet(c)),
        b.meet(c).and_then(|bc| a.meet(&bc))
    );
    prop_assert_eq!(a.join(b).join(c), a.join(&b.join(c)));
    prop_assert_eq!(a.meet(&a.join(b)), Some(a.clone()));
    prop_assert_eq!(a.meet(a), Some(a.clone()));
    prop_assert_eq!(a.join(a), a.clone());
    prop_assert_eq!(a.matches(b), a.meet(b) == Some(b.clone()));
    Ok(())
}

/// One stereo constraint of each kind, with permutation degree = `kind.degree()`.
/// Atom- and bond-centered share the leaf types; only the enum wrapper differs.
macro_rules! stereo_constraint_strategy {
    ($name:ident, $constraint:ident) => {
        pub(crate) fn $name(kind: StereoKind) -> BoxedStrategy<$constraint> {
            let degree = kind.degree();
            prop_oneof![
                (
                    permutation_strategy(degree),
                    orientation_strategy(),
                    mem_op_strategy()
                )
                    .prop_map(|(perm, orientation, mem)| $constraint::LigandSymmetry(
                        LigandSymmetryAst {
                            perm: OrientedPermutationAst {
                                perm: PermutationAst(perm),
                                orientation,
                            },
                            mem,
                        }
                    )),
                permutation_strategy(degree).prop_map(|perm| $constraint::Fluxionality(
                    FluxionalityAst {
                        perm: PermutationAst(perm),
                    }
                )),
                (ligand_pair_strategy(degree), topicity_relation_strategy())
                    .prop_map(|(pair, rel)| $constraint::Topicity(TopicityAst { pair, rel })),
                stereogenicity_relation_strategy()
                    .prop_map(|rel| $constraint::Stereogenicity(StereogenicityAst(rel))),
            ]
            .boxed()
        }
    };
}

stereo_constraint_strategy! { stereo_atom_constraint_strategy, StereoAtomConstraint }
stereo_constraint_strategy! { stereo_bond_constraint_strategy, StereoBondConstraint }

pub(crate) fn stereo_atom_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = Vec<StereoAtomConstraint>> {
    prop::collection::vec(stereo_atom_constraint_strategy(kind), 0..=3)
}

pub(crate) fn stereo_bond_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = Vec<StereoBondConstraint>> {
    prop::collection::vec(stereo_bond_constraint_strategy(kind), 0..=3)
}

pub(crate) fn stereo_atom_ast_strategy() -> impl Strategy<Value = StereoAtomAst> {
    (stereo_atom_kind_strategy(), stereo_coset_strategy()).prop_flat_map(|(kind, coset)| {
        stereo_atom_constraints_strategy(kind)
            .prop_map(move |cs| StereoAtomAst::new(kind, coset.clone()).with_constraints(cs))
    })
}

pub(crate) fn stereo_bond_ast_strategy() -> impl Strategy<Value = StereoBondAst> {
    stereo_coset_strategy().prop_flat_map(|coset| {
        stereo_bond_constraints_strategy(StereoKind::CisTrans).prop_map(move |cs| {
            StereoBondAst::new(StereoKind::CisTrans, coset.clone()).with_constraints(cs)
        })
    })
}

/// Stereo-atom entries for a molecule of `atom_count` atoms. Sites are the
/// distinct indices `0..n` (one stereo element per site, §4.1). Plain-atom
/// ligands reference any atom; virtual ligands (implicit-H / lone-pair) carry
/// the site atom (their bearing atom).
pub(crate) fn stereo_atom_entries_strategy(
    atom_count: usize,
    max: usize,
) -> BoxedStrategy<Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>> {
    if atom_count == 0 || max == 0 {
        return Just(Vec::new()).boxed();
    }
    let entry = (
        prop::collection::vec((stereo_ligand_kind_strategy(), 0..atom_count as u32), 1..=4),
        stereo_atom_ast_strategy(),
    );
    prop::collection::vec(entry, 0..=max)
        .prop_map(|entries| {
            entries
                .into_iter()
                .enumerate()
                .map(|(i, (ligand_specs, ast))| {
                    let site = AtomId(i as u32);
                    let ligands = ligand_specs
                        .into_iter()
                        .map(|(kind, a)| match kind {
                            StereoLigandKind::Atom => StereoLigand::new(AtomId(a), kind),
                            _ => StereoLigand::new(site, kind),
                        })
                        .collect();
                    (site, ligands, ast)
                })
                .collect()
        })
        .boxed()
}

/// Stereo-bond entries. Sites are the distinct bond indices `0..n`; ligands
/// reference any atom (a virtual ligand's bearing atom is a double-bond
/// terminus, modeled here as any in-range atom — roundtrip is independent of
/// the choice).
pub(crate) fn stereo_bond_entries_strategy(
    atom_count: usize,
    bond_count: usize,
    max: usize,
) -> BoxedStrategy<Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>> {
    if atom_count == 0 || bond_count == 0 || max == 0 {
        return Just(Vec::new()).boxed();
    }
    let entry = (
        prop::collection::vec((stereo_ligand_kind_strategy(), 0..atom_count as u32), 1..=4),
        stereo_bond_ast_strategy(),
    );
    prop::collection::vec(entry, 0..=max)
        .prop_map(|entries| {
            entries
                .into_iter()
                .enumerate()
                .map(|(i, (ligand_specs, ast))| {
                    let site = BondId(i as u32);
                    let ligands = ligand_specs
                        .into_iter()
                        .map(|(kind, a)| StereoLigand::new(AtomId(a), kind))
                        .collect();
                    (site, ligands, ast)
                })
                .collect()
        })
        .boxed()
}

/// Generate k distinct atom indices in [0, atom_count).
pub(crate) fn distinct_atoms_strategy(
    atom_count: usize,
    min_k: usize,
    max_k: usize,
) -> BoxedStrategy<Vec<AtomId>> {
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
                v.into_iter().map(AtomId).collect()
            })
        })
        .boxed()
}

pub(crate) fn molecule_ast_strategy() -> impl Strategy<Value = MoleculeAst> {
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
                .map(|(&[a, b], bond)| (AtomId(a), AtomId(b), bond))
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
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), aromatic_system_ast_for(n))
                    },
                ),
                0..=aromatic_count_max,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_ast_for(n))
                    },
                ),
                0..=multicenter_count_max,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_ast_strategy(),
                ),
                0..=noncovalent_count_max,
            );

            let stereo_atoms = stereo_atom_entries_strategy(atom_count, atom_count.min(2));
            let stereo_bonds =
                stereo_bond_entries_strategy(atom_count, bond_count, bond_count.min(2));

            (
                Just(atoms),
                Just(bonds_full),
                datives,
                aromatics,
                multicenters,
                noncovalents,
                stereo_atoms,
                stereo_bonds,
            )
        })
        .prop_map(
            |(
                atoms,
                bonds,
                datives,
                aromatics,
                multicenters,
                noncovalents,
                stereo_atoms,
                stereo_bonds,
            )| {
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
                    stereo_atoms,
                    stereo_bonds,
                    Constraints::new(),
                )
            },
        )
}

/// Per-entity counts for a generated `MoleculeAst`. Carried into the
/// constraint generator so that ref-bearing variants (`Constraint::Atom`,
/// relational constraints, anchors) only emit valid in-bounds indices.
#[derive(Clone, Copy)]
pub(crate) struct ConstraintCounts {
    atom: usize,
    bond: usize,
    dative: usize,
    aromatic: usize,
    multicenter: usize,
    noncovalent: usize,
    stereo_atom: usize,
    stereo_bond: usize,
}

impl ConstraintCounts {
    fn from_ast(ast: &MoleculeAst) -> Self {
        Self {
            atom: ast.atoms().count(),
            bond: ast.bonds().count(),
            dative: ast.dative_bonds().count(),
            aromatic: ast.aromatic_systems().count(),
            multicenter: ast.multicenter_bonds().count(),
            noncovalent: ast.noncovalent_bonds().count(),
            stereo_atom: ast.stereo_atoms().count(),
            stereo_bond: ast.stereo_bonds().count(),
        }
    }
}

pub(crate) fn atom_idx_strategy(atom_count: usize) -> BoxedStrategy<AtomId> {
    (0u32..atom_count as u32).prop_map(AtomId).boxed()
}

pub(crate) fn bond_idx_strategy(bond_count: usize) -> BoxedStrategy<BondId> {
    (0u32..bond_count as u32).prop_map(BondId).boxed()
}

pub(crate) fn dative_bond_idx_strategy(count: usize) -> BoxedStrategy<DativeBondId> {
    (0u32..count as u32).prop_map(DativeBondId).boxed()
}

pub(crate) fn aromatic_system_idx_strategy(count: usize) -> BoxedStrategy<AromaticSystemId> {
    (0u32..count as u32).prop_map(AromaticSystemId).boxed()
}

pub(crate) fn multicenter_bond_idx_strategy(count: usize) -> BoxedStrategy<MulticenterBondId> {
    (0u32..count as u32).prop_map(MulticenterBondId).boxed()
}

pub(crate) fn noncovalent_bond_idx_strategy(count: usize) -> BoxedStrategy<NoncovalentBondId> {
    (0u32..count as u32).prop_map(NoncovalentBondId).boxed()
}

pub(crate) fn stereo_atom_idx_strategy(count: usize) -> BoxedStrategy<StereoAtomId> {
    (0u32..count as u32).prop_map(StereoAtomId).boxed()
}

pub(crate) fn stereo_bond_idx_strategy(count: usize) -> BoxedStrategy<StereoBondId> {
    (0u32..count as u32).prop_map(StereoBondId).boxed()
}

/// Non-recursive constraint leaves: every value-only and relational
/// variant. Combinators wrap these in `constraint_strategy` below.
pub(crate) fn constraint_leaf_strategy(counts: ConstraintCounts) -> BoxedStrategy<Constraint> {
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
        let molecule_charge_sum = (
            optional_atoms.clone(),
            constraint_inner_value_strategy(-3..=3),
        )
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

    if counts.aromatic > 0 {
        let system_idx = aromatic_system_idx_strategy(counts.aromatic);

        let aromatic_leaf = (system_idx.clone(), electron_count_value_strategy(0..=12))
            .prop_map(|(system, v)| {
                Constraint::AromaticSystem(system, AromaticSystemConstraint::ElectronCount(v))
            })
            .boxed();
        choices.push(aromatic_leaf);

        if counts.atom > 0 {
            let atom_idx = atom_idx_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);
            let atoms_vec = prop::collection::vec(atom_idx.clone(), 1..=max_atoms);

            let atoms = (system_idx.clone(), atoms_vec.clone())
                .prop_map(|(system, atoms)| {
                    Constraint::Relational(RelationalConstraint::AromaticSystemAtoms {
                        system,
                        atoms,
                    })
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
    }

    if counts.multicenter > 0 {
        let bond_idx = multicenter_bond_idx_strategy(counts.multicenter);

        let multicenter_leaf = (bond_idx.clone(), electron_count_value_strategy(0..=8))
            .prop_map(|(bond, v)| {
                Constraint::MulticenterBond(bond, MulticenterBondConstraint::ElectronCount(v))
            })
            .boxed();
        choices.push(multicenter_leaf);

        if counts.atom > 0 {
            let atom_idx = atom_idx_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);
            let atoms_vec = prop::collection::vec(atom_idx.clone(), 1..=max_atoms);

            let atoms = (bond_idx.clone(), atoms_vec.clone())
                .prop_map(|(bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondAtoms {
                        bond,
                        atoms,
                    })
                })
                .boxed();
            let contains = (bond_idx.clone(), atom_idx)
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondContains {
                        bond,
                        atom,
                    })
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

    if counts.stereo_atom > 0 && counts.atom > 0 {
        let sa_idx = stereo_atom_idx_strategy(counts.stereo_atom);
        let atom_idx = atom_idx_strategy(counts.atom);
        let max_atoms = counts.atom.min(3);

        let site = (sa_idx.clone(), atom_idx.clone())
            .prop_map(|(stereo_atom, atom)| {
                Constraint::Relational(RelationalConstraint::StereoAtomSite { stereo_atom, atom })
            })
            .boxed();
        let contains = (sa_idx.clone(), atom_idx.clone())
            .prop_map(|(stereo_atom, atom)| {
                Constraint::Relational(RelationalConstraint::StereoAtomContains {
                    stereo_atom,
                    atom,
                })
            })
            .boxed();
        let ligands = (
            sa_idx.clone(),
            prop::collection::vec(atom_idx, 1..=max_atoms),
        )
            .prop_map(|(stereo_atom, atoms)| {
                Constraint::Relational(RelationalConstraint::StereoAtomLigands {
                    stereo_atom,
                    atoms,
                })
            })
            .boxed();
        let all_ligands = (sa_idx.clone(), atom_constraint_strategy())
            .prop_map(|(stereo_atom, c)| {
                Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                    stereo_atom,
                    predicate: Box::new(c),
                })
            })
            .boxed();
        let any_ligand = (sa_idx, atom_constraint_strategy())
            .prop_map(|(stereo_atom, c)| {
                Constraint::Relational(RelationalConstraint::StereoAtomAnyLigand {
                    stereo_atom,
                    predicate: Box::new(c),
                })
            })
            .boxed();
        choices.push(site);
        choices.push(contains);
        choices.push(ligands);
        choices.push(all_ligands);
        choices.push(any_ligand);
    }

    if counts.stereo_bond > 0 {
        let sb_idx = stereo_bond_idx_strategy(counts.stereo_bond);

        if counts.bond > 0 {
            let bond_idx = bond_idx_strategy(counts.bond);
            let site = (sb_idx.clone(), bond_idx)
                .prop_map(|(stereo_bond, bond)| {
                    Constraint::Relational(RelationalConstraint::StereoBondSite {
                        stereo_bond,
                        bond,
                    })
                })
                .boxed();
            choices.push(site);
        }

        if counts.atom > 0 {
            let atom_idx = atom_idx_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);

            let contains = (sb_idx.clone(), atom_idx.clone())
                .prop_map(|(stereo_bond, atom)| {
                    Constraint::Relational(RelationalConstraint::StereoBondContains {
                        stereo_bond,
                        atom,
                    })
                })
                .boxed();
            let ligands = (
                sb_idx.clone(),
                prop::collection::vec(atom_idx, 1..=max_atoms),
            )
                .prop_map(|(stereo_bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::StereoBondLigands {
                        stereo_bond,
                        atoms,
                    })
                })
                .boxed();
            let all_ligands = (sb_idx.clone(), atom_constraint_strategy())
                .prop_map(|(stereo_bond, c)| {
                    Constraint::Relational(RelationalConstraint::StereoBondAllLigands {
                        stereo_bond,
                        predicate: Box::new(c),
                    })
                })
                .boxed();
            let any_ligand = (sb_idx, atom_constraint_strategy())
                .prop_map(|(stereo_bond, c)| {
                    Constraint::Relational(RelationalConstraint::StereoBondAnyLigand {
                        stereo_bond,
                        predicate: Box::new(c),
                    })
                })
                .boxed();
            choices.push(contains);
            choices.push(ligands);
            choices.push(all_ligands);
            choices.push(any_ligand);
        }
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
pub(crate) fn sub_pattern_anchor_strategy(
    target: ConstraintCounts,
    pattern: ConstraintCounts,
) -> BoxedStrategy<SubPatternAnchor> {
    let atom_pairs = target.atom.min(pattern.atom).min(2);
    let bond_pairs = target.bond.min(pattern.bond).min(2);
    let dative_pairs = target.dative.min(pattern.dative).min(1);
    let aromatic_pairs = target.aromatic.min(pattern.aromatic).min(1);
    let multicenter_pairs = target.multicenter.min(pattern.multicenter).min(1);
    let noncovalent_pairs = target.noncovalent.min(pattern.noncovalent).min(1);
    let stereo_atom_pairs = target.stereo_atom.min(pattern.stereo_atom).min(1);
    let stereo_bond_pairs = target.stereo_bond.min(pattern.stereo_bond).min(1);
    (
        0..=atom_pairs,
        0..=bond_pairs,
        0..=dative_pairs,
        0..=aromatic_pairs,
        0..=multicenter_pairs,
        0..=noncovalent_pairs,
        0..=stereo_atom_pairs,
        0..=stereo_bond_pairs,
    )
        .prop_map(|(a, b, d, ar, mc, nc, sa, sb)| {
            let mut anchor = SubPatternAnchor::new();
            for i in 0..a {
                anchor.push_atom(AtomId(i as u32), AtomId(i as u32));
            }
            for i in 0..b {
                anchor.push_bond(BondId(i as u32), BondId(i as u32));
            }
            for i in 0..d {
                anchor.push_dative_bond(DativeBondId(i as u32), DativeBondId(i as u32));
            }
            for i in 0..ar {
                anchor.push_aromatic_system(AromaticSystemId(i as u32), AromaticSystemId(i as u32));
            }
            for i in 0..mc {
                anchor.push_multicenter_bond(
                    MulticenterBondId(i as u32),
                    MulticenterBondId(i as u32),
                );
            }
            for i in 0..nc {
                anchor.push_noncovalent_bond(
                    NoncovalentBondId(i as u32),
                    NoncovalentBondId(i as u32),
                );
            }
            for i in 0..sa {
                anchor.push_stereo_atom(StereoAtomId(i as u32), StereoAtomId(i as u32));
            }
            for i in 0..sb {
                anchor.push_stereo_bond(StereoBondId(i as u32), StereoBondId(i as u32));
            }
            anchor
        })
        .boxed()
}

/// Constraint tree: leaves wrapped in bounded-depth combinators (And/Or/Not).
pub(crate) fn constraint_strategy(counts: ConstraintCounts) -> BoxedStrategy<Constraint> {
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

pub(crate) fn molecule_ast_with_constraints_strategy() -> impl Strategy<Value = MoleculeAst> {
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

/// Generate a `Metadata` populated for an AST of the given counts. Entity
/// ids use deterministic prefixed names (`atom0`, `bond1`, ...) so that
/// names are unique across kinds and disjoint from alias names. Atom
/// aliases are capped at 3 and use a 3-element pool (`C`, `N`, `O`) for
/// the alias atom-DSL values, keeping bijectivity (each alias name
/// distinct, each alias atom distinct).
pub(crate) fn metadata_for(counts: ConstraintCounts) -> BoxedStrategy<Metadata> {
    const ALIAS_ELEMENTS: [Element; 3] = [Element::C, Element::N, Element::O];
    let id_flag = || prop::option::weighted(0.4, Just(()));
    let atom_flags = prop::collection::vec(id_flag(), counts.atom);
    let bond_flags = prop::collection::vec(id_flag(), counts.bond);
    let dative_flags = prop::collection::vec(id_flag(), counts.dative);
    let aromatic_flags = prop::collection::vec(id_flag(), counts.aromatic);
    let multicenter_flags = prop::collection::vec(id_flag(), counts.multicenter);
    let noncovalent_flags = prop::collection::vec(id_flag(), counts.noncovalent);
    (
        atom_flags,
        bond_flags,
        dative_flags,
        aromatic_flags,
        multicenter_flags,
        noncovalent_flags,
    )
        .prop_map(
            |(atoms, bonds, datives, aromatics, multicenters, noncovalents)| {
                let mut meta = Metadata::new();
                for (i, slot) in atoms.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_atom_id(AtomId(i as u32), format!("atom{i}"));
                    }
                }
                for (i, slot) in bonds.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_bond_id(BondId(i as u32), format!("bond{i}"));
                    }
                }
                for (i, slot) in datives.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_dative_bond_id(DativeBondId(i as u32), format!("dative{i}"));
                    }
                }
                for (i, slot) in aromatics.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_aromatic_system_id(
                            AromaticSystemId(i as u32),
                            format!("aromatic{i}"),
                        );
                    }
                }
                for (i, slot) in multicenters.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_multicenter_bond_id(
                            MulticenterBondId(i as u32),
                            format!("multicenter{i}"),
                        );
                    }
                }
                for (i, slot) in noncovalents.iter().enumerate() {
                    if slot.is_some() {
                        meta.set_noncovalent_bond_id(
                            NoncovalentBondId(i as u32),
                            format!("noncovalent{i}"),
                        );
                    }
                }
                for (i, element) in ALIAS_ELEMENTS.iter().enumerate() {
                    meta.add_atom_alias(format!("al{i}"), AtomAst::from_element(*element));
                }
                meta
            },
        )
        .boxed()
}

pub(crate) fn molecule_dsl_strategy() -> impl Strategy<Value = MoleculeDsl> {
    molecule_ast_with_constraints_strategy().prop_flat_map(|ast| {
        let counts = ConstraintCounts::from_ast(&ast);
        metadata_for(counts)
            .prop_map(move |metadata| MoleculeDsl::from_parts(ast.clone(), metadata))
    })
}

pub(crate) fn transaction_atom_count_strategy() -> impl Strategy<Value = usize> {
    1usize..=6
}

pub(crate) fn transaction_atoms(count: usize) -> Vec<AtomAst> {
    (0..count)
        .map(|idx| {
            let element = ELEMENTS[idx % ELEMENTS.len()];
            AtomAst::from_element(element)
        })
        .collect()
}

pub(crate) fn transaction_path_bonds(count: usize) -> Vec<AddBond> {
    (0..count.saturating_sub(1))
        .map(|idx| {
            Edit::add_bond(
                AtomRef::New(idx),
                AtomRef::New(idx + 1),
                BondAst::from_order((idx % 3 + 1) as u8),
            )
        })
        .map(|edit| match edit {
            Edit::AddBonds { mut bonds } => bonds.remove(0),
            _ => unreachable!(),
        })
        .collect()
}

pub(crate) fn transaction_path_molecule(count: usize) -> MoleculeAst {
    let atoms = transaction_atoms(count);
    let bonds = (0..count.saturating_sub(1))
        .map(|idx| {
            (
                AtomId(idx as u32),
                AtomId((idx + 1) as u32),
                BondAst::from_order((idx % 3 + 1) as u8),
            )
        })
        .collect();
    MoleculeAst::from_atoms_and_bonds(atoms, bonds)
}

pub(crate) fn transaction_add_path_edits(count: usize) -> Vec<Edit> {
    vec![
        Edit::AddAtoms {
            atoms: transaction_atoms(count),
        },
        Edit::AddBonds {
            bonds: transaction_path_bonds(count),
        },
    ]
}

#[derive(Clone, Debug)]
pub(crate) enum TransactionCase {
    AddPath {
        count: usize,
    },
    RemoveAtom {
        count: usize,
        idx: usize,
    },
    RemoveBond {
        count: usize,
        idx: usize,
    },
    SetAtomCharge {
        count: usize,
        idx: usize,
        charge: i64,
    },
    SetBondOrder {
        count: usize,
        idx: usize,
        order: u8,
    },
    AddAtomConstraint {
        count: usize,
        idx: usize,
        size: i64,
    },
    AddDativeBond {
        count: usize,
        donor: usize,
        acceptor: usize,
    },
}

impl TransactionCase {
    pub(crate) fn base(&self) -> MoleculeAst {
        match self {
            Self::AddPath { .. } => MoleculeAst::default(),
            Self::RemoveAtom { count, .. }
            | Self::RemoveBond { count, .. }
            | Self::SetAtomCharge { count, .. }
            | Self::SetBondOrder { count, .. }
            | Self::AddAtomConstraint { count, .. }
            | Self::AddDativeBond { count, .. } => transaction_path_molecule(*count),
        }
    }

    pub(crate) fn edits(&self) -> Vec<Edit> {
        match self {
            Self::AddPath { count } => transaction_add_path_edits(*count),
            Self::RemoveAtom { count, idx } => {
                vec![Edit::remove_atom(AtomRef::Id(AtomId((idx % count) as u32)))]
            }
            Self::RemoveBond { count, idx } => vec![Edit::remove_bond(BondRef::Id(BondId(
                (idx % (count - 1)) as u32,
            )))],
            Self::SetAtomCharge { count, idx, charge } => {
                vec![Edit::SetAtomField {
                    id: AtomRef::Id(AtomId((idx % count) as u32)),
                    change: AtomFieldChange::Charge {
                        old: ValueAst::default(),
                        new: ValueAst::Lit(*charge),
                    },
                }]
            }
            Self::SetBondOrder { count, idx, order } => {
                let bond_idx = idx % (count - 1);
                vec![Edit::SetBondField {
                    id: BondRef::Id(BondId(bond_idx as u32)),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit((bond_idx % 3 + 1) as i64),
                        new: ValueAst::Lit(*order as i64),
                    },
                }]
            }
            Self::AddAtomConstraint { count, idx, size } => {
                vec![Edit::AddAtomConstraint {
                    id: AtomRef::Id(AtomId((idx % count) as u32)),
                    constraint: AtomConstraint::ring_size(*size),
                }]
            }
            Self::AddDativeBond {
                count,
                donor,
                acceptor,
            } => {
                let donor = donor % count;
                let mut acceptor = acceptor % count;
                if acceptor == donor {
                    acceptor = (acceptor + 1) % count;
                }
                vec![Edit::AddDativeBond {
                    atoms: vec![
                        AtomRef::Id(AtomId(donor as u32)),
                        AtomRef::Id(AtomId(acceptor as u32)),
                    ],
                    ast: DativeBondAst::from_order(1),
                }]
            }
        }
    }
}

pub(crate) fn transaction_case_strategy() -> impl Strategy<Value = TransactionCase> {
    prop_oneof![
        transaction_atom_count_strategy().prop_map(|count| TransactionCase::AddPath { count }),
        (1usize..=6, 0usize..6).prop_map(|(count, idx)| TransactionCase::RemoveAtom { count, idx }),
        (2usize..=6, 0usize..5).prop_map(|(count, idx)| TransactionCase::RemoveBond { count, idx }),
        (1usize..=6, 0usize..6, -3i64..=3).prop_map(|(count, idx, charge)| {
            TransactionCase::SetAtomCharge { count, idx, charge }
        }),
        (2usize..=6, 0usize..5, 1u8..=3).prop_map(|(count, idx, order)| {
            TransactionCase::SetBondOrder { count, idx, order }
        }),
        (1usize..=6, 0usize..6, 3i64..=8).prop_map(|(count, idx, size)| {
            TransactionCase::AddAtomConstraint { count, idx, size }
        }),
        (2usize..=6, 0usize..6, 0usize..6).prop_map(|(count, donor, acceptor)| {
            TransactionCase::AddDativeBond {
                count,
                donor,
                acceptor,
            }
        }),
    ]
}

