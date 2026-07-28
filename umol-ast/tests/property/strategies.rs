//! Shared proptest generators for the umol-ast property suite. Domain imports
//! are re-exported (`pub(crate) use`) so the per-area test modules need only
//! `use proptest::prelude::*; use crate::strategies::*;`.

use std::collections::BTreeSet;
pub(crate) use std::collections::HashSet;
pub(crate) use std::fmt::Debug;
use std::hash::Hash;
pub(crate) use std::iter::repeat_with;
pub(crate) use std::ops::RangeInclusive;

use proptest::bool::weighted;
use proptest::prelude::*;
pub(crate) use umol_ast::ast::{
    AddBond, AromaticSystemAst, AromaticSystemConstraintAst, AromaticSystemConstraintKey,
    AromaticSystemConstraintsAst, AromaticSystemDelta, AromaticSystemFieldChange,
    AromaticSystemHandle, AromaticSystemId, AromaticSystemUpdate, AromaticValenceAst, AtomAst,
    AtomConstraintAst, AtomConstraintKey, AtomConstraintsAst, AtomDelta, AtomFieldChange,
    AtomHandle, AtomId, AtomUpdate, BondAst, BondConstraintAst, BondConstraintKey,
    BondConstraintsAst, BondDelta, BondFieldChange, BondHandle, BondId, BondUpdate, BooleanAst,
    Canonicalize, CisTransStereoAst, Constraint, Constraints, DativeBondAst,
    DativeBondConstraintAst, DativeBondConstraintKey, DativeBondConstraintsAst, DativeBondDelta,
    DativeBondFieldChange, DativeBondHandle, DativeBondId, DativeBondUpdate, Delta, Deltas,
    DpoValidator, Edit, ElectronCountsAst, ElementAst, Entity, FluxionalityAst, FromAst, IntoAst,
    IsotopeMassAst, Lattice, LigandPermutation, LigandSymmetryAst, MemOp, MoleculeAst,
    MoleculeConstraint, MoleculeCorrespondence, MoleculeParts, MulticenterBondAst,
    MulticenterBondConstraintAst, MulticenterBondConstraintKey, MulticenterBondConstraintsAst,
    MulticenterBondDelta, MulticenterBondFieldChange, MulticenterBondHandle, MulticenterBondId,
    MulticenterBondUpdate, MulticenterValenceAst, NoncovalentBondAst, NoncovalentBondConstraintAst,
    NoncovalentBondConstraintsAst, NoncovalentBondDelta, NoncovalentBondHandle, NoncovalentBondId,
    NoncovalentBondKind, NoncovalentBondKindAst, NoncovalentBondUpdate, OrientedLigandPermutation,
    ReactionAst, ReactionSpanAst, RelOp, RelationalConstraint, RingMembershipAst, RingScope,
    SpinStateAst, SpinStateUpdate, StereoAtomAst, StereoAtomConstraintAst,
    StereoAtomConstraintsAst, StereoAtomDelta, StereoAtomFieldChange, StereoAtomId,
    StereoAtomUpdate, StereoBondAst, StereoBondConstraintAst, StereoBondConstraintsAst,
    StereoBondDelta, StereoBondFieldChange, StereoBondId, StereoBondUpdate, StereoConfigurationAst,
    StereoConfigurationUpdate, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    StereoLigandPair, StereoLigandPosition, Stereogenicity, StereogenicityAst, SubPatternAnchor,
    TetrahedralStereoAst, Topicity, TopicityAst, TopicityRelationAst, ValueAst, ValuePredicate,
    ValueTerm,
};
pub(crate) use umol_ast::dsl::{
    parse_value, AromaticSystemDsl, AromaticSystemUpdateDsl, AtomDsl, AtomUpdateDsl, BondDsl,
    BondUpdateDsl, DativeBondDsl, DativeBondParticipants, DativeBondUpdateDsl, MoleculeContext,
    MoleculeDefaults, MoleculeDsl, MoleculeMetadata, MulticenterBondDsl, MulticenterBondUpdateDsl,
    NoncovalentBondDsl, NoncovalentBondUpdateDsl, ParseError, ReactionDefaults, ReactionDsl,
    ReactionSpanDsl, StereoAtomConstraintDsl, StereoAtomDsl, StereoAtomParticipants,
    StereoAtomUpdateDsl, StereoBondConstraintDsl, StereoBondDsl, StereoBondParticipants,
    StereoBondUpdateDsl, StereoLigandRef, ValueDsl,
};
pub(crate) use umol_chem::element::Element;
pub(crate) use umol_edn::{read_string, Edn, FromEdn, ToEdn};
use umol_graph_core::{Correspondence, EdgeId, NodeId};
pub(crate) use umol_perm::{Orientation, Permutation};

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
        2 => prop::sample::subsequence(Element::all().to_vec(), 1..=118).prop_map(ElementAst::lit_set),
        1 => prop::sample::subsequence(Element::all().to_vec(), 1..=118).prop_map(ElementAst::not_set),
        1 => id_strategy().prop_map(ElementAst::var),
        1 => (id_strategy(), prop::sample::subsequence(Element::all().to_vec(), 1..=118))
            .prop_map(|(id, set)| ElementAst::var_in(id, set)),
        1 => (id_strategy(), prop::sample::subsequence(Element::all().to_vec(), 1..=118))
            .prop_map(|(id, set)| ElementAst::var_not_in(id, set)),
    ]
    .prop_map(|e| e.canonicalize().unwrap_or(ElementAst::Undetermined))
}

pub(crate) fn value_basic(range: RangeInclusive<i64>) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        4 => Just(ValueAst::Undetermined),
        4 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 2..=3).prop_map(ValueAst::lit_set),
        2 => value_term_strategy().prop_map(ValueAst::term),
        2 => value_predicate_strategy().prop_map(ValueAst::predicate),
    ]
    .prop_map(canonicalize_value)
}

/// Full `ValueTerm` grammar: `Lit`/`Var` leaves under `Neg`, n-ary `Sum`/
/// `Product`, and binary `Div`/`Rem`. Generated raw; `value_basic` canonicalizes.
fn value_term_strategy() -> BoxedStrategy<ValueTerm> {
    let leaf = prop_oneof![
        (-10i64..=10).prop_map(ValueTerm::Lit),
        id_strategy().prop_map(ValueTerm::Var),
    ]
    .boxed();
    leaf.prop_recursive(3, 8, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|t| ValueTerm::Neg(Box::new(t))),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(ValueTerm::Sum),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(ValueTerm::Product),
            (inner.clone(), inner.clone())
                .prop_map(|(a, b)| ValueTerm::Div(Box::new(a), Box::new(b))),
            (inner.clone(), inner).prop_map(|(a, b)| ValueTerm::Rem(Box::new(a), Box::new(b))),
        ]
        .boxed()
    })
    .boxed()
}

/// Full `ValuePredicate` grammar: `Rel`/`Mem` leaves over terms, under `Not`,
/// `And`, `Or`. Generated raw; `value_basic` canonicalizes (folding/NNF).
fn value_predicate_strategy() -> BoxedStrategy<ValuePredicate> {
    let term = value_term_strategy();
    let leaf = prop_oneof![
        (term.clone(), rel_op_strategy(), term.clone())
            .prop_map(|(a, op, b)| ValuePredicate::Rel(a, op, b)),
        (
            term,
            mem_op_strategy(),
            prop::collection::vec(-10i64..=10, 1..=3)
        )
            .prop_map(|(e, op, s)| ValuePredicate::Mem(e, op, s.into_iter().collect())),
    ]
    .boxed();
    leaf.prop_recursive(2, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|p| ValuePredicate::Not(Box::new(p))),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(ValuePredicate::And),
            prop::collection::vec(inner, 2..=3).prop_map(ValuePredicate::Or),
        ]
        .boxed()
    })
    .boxed()
}

fn rel_op_strategy() -> impl Strategy<Value = RelOp> {
    prop_oneof![
        Just(RelOp::Le),
        Just(RelOp::Ge),
        Just(RelOp::Eq),
        Just(RelOp::Lt),
        Just(RelOp::Gt),
        Just(RelOp::Ne),
    ]
}

/// Canonicalize a generated value so the property suite operates on canonical
/// forms: the lattice laws and the render/parse identity compare against the
/// generated value itself, so a non-canonical input would spuriously fail.
/// The unsatisfiable case is unreachable for these generators.
fn canonicalize_value(v: ValueAst) -> ValueAst {
    v.canonicalize().unwrap_or(ValueAst::Undetermined)
}

pub(crate) fn any_value_ast_strategy() -> BoxedStrategy<ValueAst> {
    value_basic(-10..=10).boxed()
}

/// Possibly **non-canonical** (but satisfiable) `ValueAst`: unlike `value_basic`
/// it does not canonicalize, so it exercises the input-canonicality-independent
/// lattice laws on raw `Term`/`Predicate` forms. Unsatisfiable draws are filtered
/// out — on an unsatisfiable target the `matches` law's meet-derived RHS only
/// agrees with the default for satisfiable targets.
pub(crate) fn raw_value_ast_strategy() -> BoxedStrategy<ValueAst> {
    prop_oneof![
        2 => Just(ValueAst::Undetermined),
        2 => (-10i64..=10).prop_map(ValueAst::Lit),
        3 => value_term_strategy().prop_map(ValueAst::term),
        3 => value_predicate_strategy().prop_map(ValueAst::predicate),
    ]
    .prop_filter("satisfiable", |v| v.clone().canonicalize().is_ok())
    .boxed()
}

// Raw (non-canonical, satisfiable) generators for the remaining
// `canonical()`-overriding leaves, to fuzz the *fold* path of the universal
// lattice laws (the canonicalizing strategies only ever reach the borrow path).
// Each mixes deliberately non-canonical draws with the canonical strategy, then
// filters to satisfiable values (the `matches` law's RHS only agrees with the
// default on satisfiable targets).

pub(crate) fn raw_element_ast_strategy() -> BoxedStrategy<ElementAst> {
    prop_oneof![
        3 => element_strategy().prop_map(|e| ElementAst::lit_set([e])),
        3 => id_strategy().prop_map(|id| ElementAst::var_in(id, Element::all().to_vec())),
        2 => element_ast_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_isotope_strategy() -> BoxedStrategy<IsotopeMassAst> {
    prop_oneof![
        3 => (1u32..=250).prop_map(|m| IsotopeMassAst::lit_set(vec![m])),
        2 => isotope_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_aromatic_valence_ast_strategy() -> BoxedStrategy<AromaticValenceAst> {
    prop_oneof![
        3 => raw_value_ast_strategy().prop_map(AromaticValenceAst::Aromatic),
        2 => aromatic_valence_ast_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_multicenter_valence_ast_strategy() -> BoxedStrategy<MulticenterValenceAst> {
    prop_oneof![
        3 => raw_value_ast_strategy().prop_map(MulticenterValenceAst::Multicenter),
        2 => multicenter_valence_ast_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_tetrahedral_stereo_strategy() -> BoxedStrategy<TetrahedralStereoAst> {
    prop_oneof![
        Just(TetrahedralStereoAst::Undetermined),
        Just(TetrahedralStereoAst::NotStereo),
        stereo_coset_strategy().prop_map(TetrahedralStereoAst::Stereo),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_cis_trans_stereo_strategy() -> BoxedStrategy<CisTransStereoAst> {
    prop_oneof![
        Just(CisTransStereoAst::Undetermined),
        Just(CisTransStereoAst::NotStereo),
        stereo_coset_strategy().prop_map(CisTransStereoAst::Stereo),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_stereo_configuration_strategy() -> BoxedStrategy<StereoConfigurationAst> {
    prop_oneof![
        Just(StereoConfigurationAst::Undetermined),
        (stereo_atom_kind_strategy(), stereo_coset_strategy())
            .prop_map(|(kind, coset)| StereoConfigurationAst::kinded(kind, coset)),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_topicity_relation_strategy() -> BoxedStrategy<TopicityRelationAst> {
    prop_oneof![
        2 => Just(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic]))),
        2 => Just(TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Diastereotopic]))),
        3 => topicity_relation_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn raw_stereogenicity_relation_strategy() -> BoxedStrategy<StereogenicityAst> {
    prop_oneof![
        2 => Just(StereogenicityAst::LitSet(BTreeSet::from([Stereogenicity::Symmetric]))),
        2 => Just(StereogenicityAst::LitSet(BTreeSet::from([Stereogenicity::Stereogenic]))),
        3 => stereogenicity_relation_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().canonicalize().is_ok())
    .boxed()
}

pub(crate) fn isotope_strategy() -> impl Strategy<Value = IsotopeMassAst> {
    prop_oneof![
        3 => Just(IsotopeMassAst::Natural),
        3 => Just(IsotopeMassAst::Undetermined),
        3 => (1u32..=250).prop_map(IsotopeMassAst::Lit),
        2 => prop::collection::vec(1u32..=250, 1..=3).prop_map(IsotopeMassAst::lit_set),
        1 => id_strategy().prop_map(IsotopeMassAst::var),
        1 => (id_strategy(), prop::collection::vec(1u32..=250, 1..=3))
            .prop_map(|(id, v)| IsotopeMassAst::var_in(id, v)),
    ]
    .prop_map(|i| i.canonicalize().unwrap_or(IsotopeMassAst::Undetermined))
}

pub(crate) fn spin_state_strategy() -> impl Strategy<Value = SpinStateAst> {
    // DSL preserves spin fields field-wise. Physical (u, m) parity is
    // not a parse-time check, so any independent pair must roundtrip.
    (value_basic(0..=6), value_basic(1..=7)).prop_map(|(u, m)| SpinStateAst {
        unpaired: u,
        multiplicity: m,
    })
}

pub(crate) fn spin_state_update_strategy() -> impl Strategy<Value = SpinStateUpdate> {
    (
        prop::option::of(value_basic(0..=6)),
        prop::option::of(value_basic(1..=7)),
    )
        .prop_map(|(unpaired, multiplicity)| SpinStateUpdate {
            unpaired,
            multiplicity,
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
/// `Lit`, and `LitSet`. No symbolic `Term`/`Predicate` — the constraint
/// formatters route to `fmt_value_field_required` / `fmt_ring_count` / the
/// various `#r` blocks, and a `Term(Lit(n))` would render to a pure integer
/// that the parser then re-reads as a plain `Lit`, breaking roundtrip. The
/// molecule-level EDN tests cover symbolic values on constraint values through
/// the tree-based path, so the gap is contained.
pub(crate) fn constraint_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => Just(ValueAst::Undetermined),
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::lit_set(v)
        }),
    ]
}

/// `Lit`/`Set` only — still used by the ring-size strategies where
/// `Undetermined` on the inner value collapses into a dropped constraint
/// in the entity-level formatter (see vacuous `RingMembership(_, Undetermined)`, intentionally dropped).
pub(crate) fn constraint_inner_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        range.clone().prop_map(ValueAst::Lit),
        prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::lit_set(v)
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
    .prop_map(|v| v.canonicalize().unwrap_or(AromaticValenceAst::Undetermined))
}

pub(crate) fn multicenter_valence_ast_strategy() -> impl Strategy<Value = MulticenterValenceAst> {
    prop_oneof![
        Just(MulticenterValenceAst::NotMulticenter),
        constraint_value_strategy(0..=6).prop_map(MulticenterValenceAst::Multicenter),
    ]
    .prop_map(|v| {
        v.canonicalize()
            .unwrap_or(MulticenterValenceAst::Undetermined)
    })
}

/// Atom constraints route through `fmt_value_field_required` (or
/// `fmt_ring_count` for `#R`), which elide vacuous (Undetermined) payloads
/// per the canonical-rendering rule. Generators excluding `Undetermined`
/// keep the render → reparse identity intact.
pub(crate) fn atom_constraint_strategy() -> BoxedStrategy<AtomConstraintAst> {
    prop_oneof![
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::Valence),
        constraint_inner_value_strategy(0..=8).prop_map(AtomConstraintAst::TotalValence),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraintAst::DonatedPairs),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraintAst::AcceptedPairs),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::Degree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::TotalDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::RingDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::RingValence),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintAst::TotalHydrogens),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| AtomConstraintAst::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6))
            .prop_map(|(s, count)| AtomConstraintAst::ring_membership(RingScope::Size(s), count)),
        aromatic_valence_ast_strategy().prop_map(AtomConstraintAst::AromaticValence),
        multicenter_valence_ast_strategy().prop_map(AtomConstraintAst::MulticenterValence),
        tetrahedral_stereo_strategy().prop_map(AtomConstraintAst::TetrahedralStereo),
    ]
    .boxed()
}

pub(crate) fn atom_constraints_strategy() -> impl Strategy<Value = AtomConstraintsAst> {
    prop::collection::vec(atom_constraint_strategy(), 0..=3).prop_map(|list| {
        let mut cs = AtomConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn atom_update_constraints_strategy() -> impl Strategy<Value = AtomConstraintsAst> {
    prop::collection::vec(
        prop_oneof![
            atom_constraint_strategy(),
            atom_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(AtomConstraintsAst::from_iter)
}

pub(crate) fn bond_constraint_strategy() -> BoxedStrategy<BondConstraintAst> {
    prop_oneof![
        any::<bool>().prop_map(|b| BondConstraintAst::Aromatic(BooleanAst::Lit(b))),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| BondConstraintAst::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6))
            .prop_map(|(s, count)| BondConstraintAst::ring_membership(RingScope::Size(s), count)),
        cis_trans_stereo_strategy().prop_map(BondConstraintAst::CisTransStereo),
    ]
    .boxed()
}

pub(crate) fn bond_constraints_strategy() -> impl Strategy<Value = BondConstraintsAst> {
    prop::collection::vec(bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = BondConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn bond_update_constraints_strategy() -> impl Strategy<Value = BondConstraintsAst> {
    prop::collection::vec(
        prop_oneof![
            bond_constraint_strategy(),
            bond_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=2,
    )
    .prop_map(BondConstraintsAst::from_iter)
}

pub(crate) fn dative_bond_constraint_strategy() -> BoxedStrategy<DativeBondConstraintAst> {
    prop_oneof![
        any::<bool>().prop_map(|b| DativeBondConstraintAst::Aromatic(BooleanAst::Lit(b))),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| DativeBondConstraintAst::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6)).prop_map(|(s, count)| {
            DativeBondConstraintAst::ring_membership(RingScope::Size(s), count)
        }),
    ]
    .boxed()
}

pub(crate) fn dative_bond_constraints_strategy() -> impl Strategy<Value = DativeBondConstraintsAst>
{
    prop::collection::vec(dative_bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = DativeBondConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn dative_bond_update_constraints_strategy(
) -> impl Strategy<Value = DativeBondConstraintsAst> {
    prop::collection::vec(
        prop_oneof![
            dative_bond_constraint_strategy(),
            dative_bond_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=2,
    )
    .prop_map(DativeBondConstraintsAst::from_iter)
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
    pub(crate) fn atom_update_strategy()
    (
        element in prop::option::of(element_ast_strategy()),
        isotope_mass in prop::option::of(isotope_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        implicit_hydrogens in prop::option::of(value_basic(0..=4)),
        lone_pairs in prop::option::of(value_basic(0..=4)),
        spin in spin_state_update_strategy(),
        constraints in atom_update_constraints_strategy(),
    ) -> AtomUpdate {
        AtomUpdate {
            element,
            isotope_mass,
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
        spin in spin_state_strategy(),
        constraints in bond_constraints_strategy(),
    ) -> BondAst {
        BondAst {
            order,
            charge,
            spin,
            constraints,
        }
    }
}

prop_compose! {
    pub(crate) fn bond_update_strategy()
    (
        order in prop::option::of(value_basic(1..=4)),
        charge in prop::option::of(value_basic(-1..=1)),
        spin in spin_state_update_strategy(),
        constraints in bond_update_constraints_strategy(),
    ) -> BondUpdate {
        BondUpdate {
            order,
            charge,
            spin,
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
        Just(
            BondAst::new(ValueAst::Lit(1))
                .with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ),
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
    // Order is sampled from the small literal range that the DSL keyword
    // shorthands cover (`:single` / `:double` / `:triple`), keeping
    // canonical-form roundtrip exercised across haptic-pair counts.
    let order_strategy = prop_oneof![
        Just(ValueAst::Lit(1)),
        Just(ValueAst::Lit(2)),
        Just(ValueAst::Lit(3)),
        Just(ValueAst::Undetermined),
    ];
    (order_strategy, dative_bond_constraints_strategy())
        .prop_map(|(order, constraints)| DativeBondAst { order, constraints })
}

prop_compose! {
    pub(crate) fn dative_bond_update_strategy()
    (
        order in prop::option::of(value_basic(1..=4)),
        constraints in dative_bond_update_constraints_strategy(),
    ) -> DativeBondUpdate {
        DativeBondUpdate { order, constraints }
    }
}

/// Optional `ElectronCount` constraint (the asserted total). The strategy
/// emits `None` half the time, otherwise wraps a `ValueAst::Lit` or
/// `Set`. `Undetermined` is excluded because it has no canonical
/// surface form in the entity-string `#e<n>` field — `#e*` is admitted on
/// parse but the renderer omits the predicate entirely, breaking
/// roundtrip.
pub(crate) fn optional_aromatic_electron_count(
) -> impl Strategy<Value = AromaticSystemConstraintsAst> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=12)).prop_map(|opt| {
        let mut cs = AromaticSystemConstraintsAst::new();
        if let Some(v) = opt {
            cs.set(AromaticSystemConstraintAst::ElectronCount(v));
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn aromatic_system_update_constraints_strategy(
) -> impl Strategy<Value = AromaticSystemConstraintsAst> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            electron_count_value_strategy(0..=12),
            Just(ValueAst::Undetermined),
        ],
    )
    .prop_map(|value| {
        value
            .map(AromaticSystemConstraintAst::ElectronCount)
            .map(AromaticSystemConstraintsAst::from)
            .unwrap_or_default()
    })
}

pub(crate) fn optional_multicenter_electron_count(
) -> impl Strategy<Value = MulticenterBondConstraintsAst> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=8)).prop_map(|opt| {
        let mut cs = MulticenterBondConstraintsAst::new();
        if let Some(v) = opt {
            cs.set(MulticenterBondConstraintAst::ElectronCount(v));
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn multicenter_bond_update_constraints_strategy(
) -> impl Strategy<Value = MulticenterBondConstraintsAst> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            electron_count_value_strategy(0..=8),
            Just(ValueAst::Undetermined),
        ],
    )
    .prop_map(|value| {
        value
            .map(MulticenterBondConstraintAst::ElectronCount)
            .map(MulticenterBondConstraintsAst::from)
            .unwrap_or_default()
    })
}

pub(crate) fn electron_count_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = ValueAst> {
    prop_oneof![
        3 => range.clone().prop_map(ValueAst::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            ValueAst::lit_set(v)
        }),
    ]
}

/// Leaf strategy: `Undetermined` or a concrete `Lit` count vector (length 1–4).
pub(crate) fn electron_counts_ast_strategy() -> impl Strategy<Value = ElectronCountsAst> {
    prop_oneof![
        Just(ElectronCountsAst::Undetermined),
        prop::collection::vec(0i64..=8, 1..=4).prop_map(ElectronCountsAst::Lit),
    ]
}

/// Stand-alone strategy for entity-string roundtrip tests. `electrons` is
/// `Undetermined` because the entity string carries no per-atom data; the
/// `ElectronCount` constraint is exercised here via `#e<n>`.
pub(crate) fn aromatic_system_ast_strategy() -> impl Strategy<Value = AromaticSystemAst> {
    (value_basic(-2..=2), optional_aromatic_electron_count()).prop_map(|(charge, constraints)| {
        AromaticSystemAst {
            electrons: ElectronCountsAst::Undetermined,
            charge,
            spin: SpinStateAst::default(),
            constraints,
        }
    })
}

pub(crate) fn aromatic_system_patch_ast_strategy() -> impl Strategy<Value = AromaticSystemAst> {
    (
        electron_counts_ast_strategy(),
        value_basic(-2..=2),
        spin_state_strategy(),
        optional_aromatic_electron_count(),
    )
        .prop_map(|(electrons, charge, spin, constraints)| AromaticSystemAst {
            electrons,
            charge,
            spin,
            constraints,
        })
}

prop_compose! {
    pub(crate) fn aromatic_system_update_strategy()
    (
        electrons in prop::option::of(electron_counts_ast_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        spin in spin_state_update_strategy(),
        constraints in aromatic_system_update_constraints_strategy(),
    ) -> AromaticSystemUpdate {
        AromaticSystemUpdate { electrons, charge, spin, constraints }
    }
}

/// Atom-count-aware variant: generates an `AromaticSystemAst` whose
/// `electrons` `Lit` vector has exactly `atom_count` entries. Includes an
/// optional `ElectronCount` constraint so the molecule-level prop tests
/// exercise both the per-atom counts and the asserted total in the same pass.
pub(crate) fn aromatic_system_ast_for(
    atom_count: usize,
) -> impl Strategy<Value = AromaticSystemAst> {
    (
        value_basic(-2..=2),
        prop::collection::vec(0i64..=2, atom_count),
        optional_aromatic_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| AromaticSystemAst {
            electrons: ElectronCountsAst::Lit(electrons),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        })
}

pub(crate) fn multicenter_bond_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (value_basic(-2..=2), optional_multicenter_electron_count()).prop_map(
        |(charge, constraints)| MulticenterBondAst {
            electrons: ElectronCountsAst::Undetermined,
            charge,
            spin: SpinStateAst::default(),
            constraints,
        },
    )
}

pub(crate) fn multicenter_bond_patch_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (
        electron_counts_ast_strategy(),
        value_basic(-2..=2),
        spin_state_strategy(),
        optional_multicenter_electron_count(),
    )
        .prop_map(
            |(electrons, charge, spin, constraints)| MulticenterBondAst {
                electrons,
                charge,
                spin,
                constraints,
            },
        )
}

prop_compose! {
    pub(crate) fn multicenter_bond_update_strategy()
    (
        electrons in prop::option::of(electron_counts_ast_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        spin in spin_state_update_strategy(),
        constraints in multicenter_bond_update_constraints_strategy(),
    ) -> MulticenterBondUpdate {
        MulticenterBondUpdate { electrons, charge, spin, constraints }
    }
}

pub(crate) fn multicenter_bond_ast_for(
    atom_count: usize,
) -> impl Strategy<Value = MulticenterBondAst> {
    (
        value_basic(-2..=2),
        prop::collection::vec(0i64..=2, atom_count),
        optional_multicenter_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(electrons),
            charge,
            spin: SpinStateAst::default(),
            constraints,
        })
}

pub(crate) fn noncovalent_bond_kind_ast_strategy() -> impl Strategy<Value = NoncovalentBondKindAst>
{
    prop_oneof![
        Just(NoncovalentBondKindAst::Undetermined),
        prop::sample::select(NONCOVALENT_KINDS).prop_map(NoncovalentBondKindAst::Lit),
    ]
}

pub(crate) fn noncovalent_bond_constraint_strategy() -> BoxedStrategy<NoncovalentBondConstraintAst>
{
    any::<bool>()
        .prop_map(|b| NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Lit(b)))
        .boxed()
}

pub(crate) fn noncovalent_bond_constraints_strategy(
) -> impl Strategy<Value = NoncovalentBondConstraintsAst> {
    prop::collection::vec(noncovalent_bond_constraint_strategy(), 0..=1).prop_map(|list| {
        let mut cs = NoncovalentBondConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn noncovalent_bond_update_constraints_strategy(
) -> impl Strategy<Value = NoncovalentBondConstraintsAst> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            any::<bool>().prop_map(|value| {
                NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Lit(value))
            }),
            Just(NoncovalentBondConstraintAst::Intramolecular(
                BooleanAst::Undetermined,
            )),
        ],
    )
    .prop_map(|constraint| {
        constraint
            .map(NoncovalentBondConstraintsAst::from)
            .unwrap_or_default()
    })
}

pub(crate) fn noncovalent_bond_ast_strategy() -> impl Strategy<Value = NoncovalentBondAst> {
    (
        prop::sample::select(NONCOVALENT_KINDS),
        noncovalent_bond_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondAst {
            kind: NoncovalentBondKindAst::Lit(kind),
            constraints,
        })
}

pub(crate) fn noncovalent_bond_patch_ast_strategy() -> impl Strategy<Value = NoncovalentBondAst> {
    (
        noncovalent_bond_kind_ast_strategy(),
        noncovalent_bond_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondAst { kind, constraints })
}

pub(crate) fn noncovalent_bond_update_strategy() -> impl Strategy<Value = NoncovalentBondUpdate> {
    (
        prop::option::of(noncovalent_bond_kind_ast_strategy()),
        noncovalent_bond_update_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondUpdate { kind, constraints })
}

/// Coset forms that round-trip through both the entity `:type` string and the
/// EDN coset-form: `Undetermined` (`*`), `Lit`, and a literal set
/// (`{a,b,…}` ↔ EDN vector). The `~`/`^`/`?var` operator-exprs are reserved
/// (§7.14) and excluded.
pub(crate) fn stereo_coset_strategy() -> impl Strategy<Value = StereoCoset> {
    prop_oneof![
        Just(StereoCoset::Undetermined),
        (0u32..=6).prop_map(StereoCoset::Lit),
        prop::collection::vec(0u32..=6, 1..=3).prop_map(StereoCoset::lit_set),
    ]
}

pub(crate) fn stereo_ligand_kind_strategy() -> impl Strategy<Value = StereoLigandKind> {
    prop_oneof![
        Just(StereoLigandKind::Atom),
        Just(StereoLigandKind::ImplicitHydrogen),
        Just(StereoLigandKind::LonePair),
    ]
}

/// Per-kind `#T` / `#C` constraint values, excluding the vacuous `Undetermined`
/// site (it is dropped on render, breaking render → reparse). `Stereo(coset)`
/// (incl. the `+` form `Stereo(Undetermined)`) and `NotStereo` are kept. The
/// `_lattice_strategy` variants add `Undetermined` (the lattice top).
pub(crate) fn tetrahedral_stereo_strategy() -> impl Strategy<Value = TetrahedralStereoAst> {
    prop_oneof![
        Just(TetrahedralStereoAst::NotStereo),
        stereo_coset_strategy().prop_map(TetrahedralStereoAst::Stereo),
    ]
    .prop_map(|s| {
        s.canonicalize()
            .unwrap_or(TetrahedralStereoAst::Undetermined)
    })
}

pub(crate) fn tetrahedral_stereo_lattice_strategy() -> impl Strategy<Value = TetrahedralStereoAst> {
    prop_oneof![
        Just(TetrahedralStereoAst::Undetermined),
        tetrahedral_stereo_strategy(),
    ]
}

pub(crate) fn cis_trans_stereo_strategy() -> impl Strategy<Value = CisTransStereoAst> {
    prop_oneof![
        Just(CisTransStereoAst::NotStereo),
        stereo_coset_strategy().prop_map(CisTransStereoAst::Stereo),
    ]
    .prop_map(|s| s.canonicalize().unwrap_or(CisTransStereoAst::Undetermined))
}

pub(crate) fn cis_trans_stereo_lattice_strategy() -> impl Strategy<Value = CisTransStereoAst> {
    prop_oneof![
        Just(CisTransStereoAst::Undetermined),
        cis_trans_stereo_strategy(),
    ]
}

/// `StereoConfigurationAst` over the atom geometry kinds, including the kindless
/// `Undetermined` top.
pub(crate) fn stereo_configuration_lattice_strategy(
) -> impl Strategy<Value = StereoConfigurationAst> {
    prop_oneof![
        Just(StereoConfigurationAst::Undetermined),
        (stereo_atom_kind_strategy(), stereo_coset_strategy())
            .prop_map(|(kind, coset)| StereoConfigurationAst::kinded(kind, coset)),
    ]
    .prop_map(|c| {
        c.canonicalize()
            .unwrap_or(StereoConfigurationAst::Undetermined)
    })
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
    Just((0..degree).collect::<Vec<usize>>())
        .prop_shuffle()
        .prop_map(|image| Permutation::from_image(&image))
}

/// A permutation realizable by the stereo kind's parent group. This differs
/// from an arbitrary degree-matched permutation for partitioned kinds such as
/// cis/trans, whose parent preserves the two ligand sides.
pub(crate) fn stereo_frame_permutation_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = Permutation> {
    let space = kind.class_key().space();
    (0..space.group().order(), 0..space.count()).prop_map(move |(group, coset)| {
        space.group().elements()[group].compose(
            space
                .unindex(coset as u32)
                .expect("generated coset index is in range"),
        )
    })
}

pub(crate) fn orientation_strategy() -> impl Strategy<Value = Orientation> {
    prop_oneof![Just(Orientation::Proper), Just(Orientation::Improper)]
}

pub(crate) fn mem_op_strategy() -> impl Strategy<Value = MemOp> {
    prop_oneof![Just(MemOp::In), Just(MemOp::NotIn)]
}

pub(crate) fn ligand_pair_strategy(degree: usize) -> impl Strategy<Value = StereoLigandPair> {
    (0..degree as u32, 0..degree as u32)
        .prop_map(|(a, b)| StereoLigandPair::new(StereoLigandPosition(a), StereoLigandPosition(b)))
}

/// Non-vacuous topicity relations only (`Undetermined` elides on render, so it
/// would break the render → reparse roundtrip — mirrors `tetrahedral_stereo_strategy`).
pub(crate) fn topicity_relation_strategy() -> impl Strategy<Value = TopicityRelationAst> {
    prop_oneof![
        Just(TopicityRelationAst::Lit(Topicity::Homotopic)),
        Just(TopicityRelationAst::Lit(Topicity::Enantiotopic)),
        Just(TopicityRelationAst::Lit(Topicity::Diastereotopic)),
        Just(TopicityRelationAst::NotSet(BTreeSet::from([
            Topicity::Homotopic
        ]))),
        Just(TopicityRelationAst::NotSet(BTreeSet::from([
            Topicity::Enantiotopic
        ]))),
        Just(TopicityRelationAst::NotSet(BTreeSet::from([
            Topicity::Diastereotopic
        ]))),
    ]
}

pub(crate) fn stereogenicity_relation_strategy() -> impl Strategy<Value = StereogenicityAst> {
    prop_oneof![
        Just(StereogenicityAst::Lit(Stereogenicity::Symmetric)),
        Just(StereogenicityAst::Lit(Stereogenicity::Prochiral)),
        Just(StereogenicityAst::Lit(Stereogenicity::Stereogenic)),
        Just(StereogenicityAst::NotSet(BTreeSet::from([
            Stereogenicity::Symmetric
        ]))),
        Just(StereogenicityAst::NotSet(BTreeSet::from([
            Stereogenicity::Prochiral
        ]))),
        Just(StereogenicityAst::NotSet(BTreeSet::from([
            Stereogenicity::Stereogenic
        ]))),
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

pub(crate) fn stereogenicity_relation_lattice_strategy() -> impl Strategy<Value = StereogenicityAst>
{
    prop_oneof![
        Just(StereogenicityAst::Undetermined),
        stereogenicity_relation_strategy(),
    ]
}

pub(crate) fn ligand_symmetry_strategy(degree: usize) -> impl Strategy<Value = LigandSymmetryAst> {
    (
        permutation_strategy(degree),
        orientation_strategy(),
        any::<bool>(),
    )
        .prop_map(|(permutation, orientation, invariant)| LigandSymmetryAst {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(permutation),
                orientation,
            },
            invariant: BooleanAst::Lit(invariant),
        })
}

pub(crate) fn fluxionality_strategy(degree: usize) -> impl Strategy<Value = FluxionalityAst> {
    (permutation_strategy(degree), any::<bool>()).prop_map(|(permutation, active)| {
        FluxionalityAst {
            permutation: LigandPermutation(permutation),
            active: BooleanAst::Lit(active),
        }
    })
}

pub(crate) fn topicity_strategy(degree: usize) -> impl Strategy<Value = TopicityAst> {
    (
        ligand_pair_strategy(degree),
        topicity_relation_lattice_strategy(),
    )
        .prop_map(|(pair, relation)| TopicityAst { pair, relation })
}

/// Canonical, fiber-spanning `RingMembershipAst`: the `scope` varies (`All` and
/// `Size(3..=10)`) so a value triple lands in different fibers, exercising the
/// cross-scope `meet` → `None` / `join` → `Err(NoJoin)` path.
pub(crate) fn ring_membership_lattice_strategy() -> impl Strategy<Value = RingMembershipAst> {
    prop_oneof![
        constraint_value_strategy(0..=6)
            .prop_map(|count| RingMembershipAst::new(RingScope::All, count)),
        (3u8..=10, constraint_value_strategy(0..=6))
            .prop_map(|(size, count)| RingMembershipAst::new(RingScope::Size(size), count)),
    ]
    .prop_map(|membership| {
        membership
            .canonicalize()
            .expect("non-empty count strategy never contradicts")
    })
}

/// Universal lattice laws — hold for **any** inputs (canonical or not): meet/join
/// commutativity and associativity, `matches` ⇔ meet-derived, and the
/// Lattice→Canonicalize correspondence that `meet`/`join` land in canonical form.
pub(crate) fn assert_lattice_laws<L: Lattice + Debug>(
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
    prop_assert_eq!(
        a.join(b).and_then(|ab| ab.join(c)),
        b.join(c).and_then(|bc| a.join(&bc))
    );
    prop_assert_eq!(a.matches(b), a.meet(b) == b.clone().canonicalize().ok());
    if let Some(m) = a.meet(b) {
        prop_assert_eq!(m.clone().canonicalize(), Ok(m));
    }
    if let Ok(j) = a.join(b) {
        prop_assert_eq!(j.clone().canonicalize(), Ok(j));
    }
    // `canonical()` (the borrow fast-path) agrees with `canonicalize()`.
    prop_assert_eq!(
        a.canonical().map(|c| c.into_owned()),
        a.clone().canonicalize()
    );
    // `canonical_eq` is canonical equality.
    prop_assert_eq!(
        a.canonical_eq(b),
        a.clone().canonicalize().ok() == b.clone().canonicalize().ok()
    );
    Ok(())
}

/// Lattice laws that assume **canonical** inputs: each input is a `canonicalize`
/// fixpoint, plus idempotence and absorption (whose RHS is the input verbatim,
/// which only holds when the input is already canonical).
pub(crate) fn assert_canonical_lattice_laws<L: Lattice + Debug>(
    a: &L,
    b: &L,
    c: &L,
) -> Result<(), TestCaseError> {
    prop_assert_eq!(a.clone().canonicalize(), Ok(a.clone()));
    prop_assert_eq!(b.clone().canonicalize(), Ok(b.clone()));
    prop_assert_eq!(c.clone().canonicalize(), Ok(c.clone()));
    prop_assert_eq!(a.meet(a), Some(a.clone()));
    prop_assert_eq!(a.join(a), Ok(a.clone()));
    if let Ok(ab) = a.join(b) {
        prop_assert_eq!(a.meet(&ab), Some(a.clone()));
    }
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
                    any::<bool>()
                )
                    .prop_map(|(permutation, orientation, invariant)| {
                        $constraint::LigandSymmetry(LigandSymmetryAst {
                            permutation: OrientedLigandPermutation {
                                permutation: LigandPermutation(permutation),
                                orientation,
                            },
                            invariant: BooleanAst::Lit(invariant),
                        })
                    }),
                (permutation_strategy(degree), any::<bool>()).prop_map(|(permutation, active)| {
                    $constraint::Fluxionality(FluxionalityAst {
                        permutation: LigandPermutation(permutation),
                        active: BooleanAst::Lit(active),
                    })
                }),
                (ligand_pair_strategy(degree), topicity_relation_strategy()).prop_map(
                    |(pair, relation)| $constraint::Topicity(TopicityAst { pair, relation })
                ),
                stereogenicity_relation_strategy().prop_map(|rel| $constraint::Stereogenicity(rel)),
            ]
            .boxed()
        }
    };
}

stereo_constraint_strategy! { stereo_atom_constraint_strategy, StereoAtomConstraintAst }
stereo_constraint_strategy! { stereo_bond_constraint_strategy, StereoBondConstraintAst }

pub(crate) fn stereo_atom_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoAtomConstraintsAst> {
    prop::collection::vec(stereo_atom_constraint_strategy(kind), 0..=3).prop_map(|list| {
        let mut cs = StereoAtomConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn stereo_bond_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoBondConstraintsAst> {
    prop::collection::vec(stereo_bond_constraint_strategy(kind), 0..=3).prop_map(|list| {
        let mut cs = StereoBondConstraintsAst::new();
        for c in list {
            cs.set(c);
        }
        cs.canonicalize().unwrap_or_default()
    })
}

pub(crate) fn stereo_atom_ast_strategy() -> impl Strategy<Value = StereoAtomAst> {
    (stereo_atom_kind_strategy(), stereo_coset_strategy()).prop_flat_map(|(kind, coset)| {
        stereo_atom_constraints_strategy(kind)
            .prop_map(move |cs| StereoAtomAst::new(kind, coset.clone()).with_constraints(cs))
    })
}

pub(crate) fn stereo_atom_update_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoAtomConstraintsAst> {
    prop::collection::vec(
        prop_oneof![
            stereo_atom_constraint_strategy(kind),
            stereo_atom_constraint_strategy(kind)
                .prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(StereoAtomConstraintsAst::from_iter)
}

pub(crate) fn stereo_atom_update_strategy() -> impl Strategy<Value = StereoAtomUpdate> {
    prop_oneof![
        Just(StereoAtomUpdate::default()),
        Just(StereoAtomUpdate {
            configuration: StereoConfigurationUpdate::Undetermined,
            ..Default::default()
        }),
        (
            stereo_atom_kind_strategy(),
            prop::option::of(stereo_coset_strategy())
        )
            .prop_flat_map(|(kind, coset)| {
                stereo_atom_update_constraints_strategy(kind).prop_map(move |constraints| {
                    StereoAtomUpdate {
                        configuration: StereoConfigurationUpdate::Kinded {
                            kind,
                            coset: coset.clone(),
                        },
                        constraints,
                    }
                })
            },),
    ]
}

pub(crate) fn stereo_bond_ast_strategy() -> impl Strategy<Value = StereoBondAst> {
    stereo_coset_strategy().prop_flat_map(|coset| {
        stereo_bond_constraints_strategy(StereoKind::CisTrans).prop_map(move |cs| {
            StereoBondAst::new(StereoKind::CisTrans, coset.clone()).with_constraints(cs)
        })
    })
}

pub(crate) fn stereo_bond_update_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoBondConstraintsAst> {
    prop::collection::vec(
        prop_oneof![
            stereo_bond_constraint_strategy(kind),
            stereo_bond_constraint_strategy(kind)
                .prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(StereoBondConstraintsAst::from_iter)
}

pub(crate) fn stereo_bond_update_strategy() -> impl Strategy<Value = StereoBondUpdate> {
    prop_oneof![
        Just(StereoBondUpdate::default()),
        Just(StereoBondUpdate {
            configuration: StereoConfigurationUpdate::Undetermined,
            ..Default::default()
        }),
        prop::option::of(stereo_coset_strategy()).prop_flat_map(|coset| {
            stereo_bond_update_constraints_strategy(StereoKind::CisTrans).prop_map(
                move |constraints| StereoBondUpdate {
                    configuration: StereoConfigurationUpdate::Kinded {
                        kind: StereoKind::CisTrans,
                        coset: coset.clone(),
                    },
                    constraints,
                },
            )
        }),
    ]
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
                MoleculeAst::from_parts(MoleculeParts {
                    atoms,
                    bonds,
                    dative: dative_triples,
                    aromatic: aromatic_entries,
                    multicenter: multicenter_entries,
                    noncovalent: noncovalent_triples,
                    stereo_atoms,
                    stereo_bonds,
                    constraints: Constraints::new(),
                })
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

pub(crate) fn atom_id_strategy(atom_count: usize) -> BoxedStrategy<AtomId> {
    (0u32..atom_count as u32).prop_map(AtomId).boxed()
}

pub(crate) fn bond_id_strategy(bond_count: usize) -> BoxedStrategy<BondId> {
    (0u32..bond_count as u32).prop_map(BondId).boxed()
}

pub(crate) fn dative_bond_id_strategy(count: usize) -> BoxedStrategy<DativeBondId> {
    (0u32..count as u32).prop_map(DativeBondId).boxed()
}

pub(crate) fn aromatic_system_id_strategy(count: usize) -> BoxedStrategy<AromaticSystemId> {
    (0u32..count as u32).prop_map(AromaticSystemId).boxed()
}

pub(crate) fn multicenter_bond_id_strategy(count: usize) -> BoxedStrategy<MulticenterBondId> {
    (0u32..count as u32).prop_map(MulticenterBondId).boxed()
}

pub(crate) fn noncovalent_bond_id_strategy(count: usize) -> BoxedStrategy<NoncovalentBondId> {
    (0u32..count as u32).prop_map(NoncovalentBondId).boxed()
}

pub(crate) fn stereo_atom_id_strategy(count: usize) -> BoxedStrategy<StereoAtomId> {
    (0u32..count as u32).prop_map(StereoAtomId).boxed()
}

pub(crate) fn stereo_bond_id_strategy(count: usize) -> BoxedStrategy<StereoBondId> {
    (0u32..count as u32).prop_map(StereoBondId).boxed()
}

/// Non-recursive constraint leaves: every value-only and relational
/// variant. Combinators wrap these in `constraint_strategy` below.
pub(crate) fn constraint_leaf_strategy(counts: ConstraintCounts) -> BoxedStrategy<Constraint> {
    let mut choices: Vec<BoxedStrategy<Constraint>> = Vec::new();

    if counts.atom > 0 {
        let atom_id = atom_id_strategy(counts.atom);

        // Constraint::Atom carrying any AtomConstraintAst variant.
        let atom_leaf = (atom_id.clone(), atom_constraint_strategy())
            .prop_map(|(id, c)| Constraint::Atom(id, c));
        choices.push(atom_leaf.boxed());

        // MoleculeConstraint variants over atom refs.
        let max_atoms = counts.atom.min(3);
        let atoms_vec = prop::collection::vec(atom_id.clone(), 1..=max_atoms);
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
        let bond_id = bond_id_strategy(counts.bond);
        let bond_leaf = (bond_id.clone(), bond_constraint_strategy())
            .prop_map(|(id, c)| Constraint::Bond(id, c))
            .boxed();
        choices.push(bond_leaf);

        let max_bonds = counts.bond.min(3);
        let optional_bonds =
            prop::option::of(prop::collection::vec(bond_id, 1..=max_bonds)).boxed();
        let molecule_bond_order_sum = (optional_bonds, constraint_inner_value_strategy(0..=8))
            .prop_map(|(bonds, sum)| {
                Constraint::Molecule(MoleculeConstraint::BondOrderSum { bonds, sum })
            })
            .boxed();
        choices.push(molecule_bond_order_sum);
    }

    if counts.dative > 0 {
        let dative_id = dative_bond_id_strategy(counts.dative);
        let dative_leaf = (dative_id.clone(), dative_bond_constraint_strategy())
            .prop_map(|(id, c)| Constraint::DativeBond(id, c))
            .boxed();
        choices.push(dative_leaf);

        if counts.atom > 0 {
            let atom_id = atom_id_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);
            let atoms_vec = prop::collection::vec(atom_id.clone(), 1..=max_atoms);
            let donors = (dative_id.clone(), atoms_vec.clone())
                .prop_map(|(bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::DativeBondDonors { bond, atoms })
                })
                .boxed();
            let donor = (dative_id.clone(), atom_id.clone())
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::DativeBondDonor { bond, atom })
                })
                .boxed();
            let contains_all_donors = (dative_id.clone(), atoms_vec)
                .prop_map(|(bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::DativeBondContainsAllDonors {
                        bond,
                        atoms,
                    })
                })
                .boxed();
            let all_donors = (dative_id.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAllDonors {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            let any_donor = (dative_id.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAnyDonor {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            let acceptor = (dative_id.clone(), atom_id)
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAcceptor { bond, atom })
                })
                .boxed();
            let acceptor_satisfies = (dative_id.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::DativeBondAcceptorSatisfies {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            choices.push(donors);
            choices.push(donor);
            choices.push(contains_all_donors);
            choices.push(all_donors);
            choices.push(any_donor);
            choices.push(acceptor);
            choices.push(acceptor_satisfies);
        }
        if counts.bond > 0 {
            let parallels = (dative_id, bond_id_strategy(counts.bond))
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
        let system_id = aromatic_system_id_strategy(counts.aromatic);

        let aromatic_leaf = (system_id.clone(), electron_count_value_strategy(0..=12))
            .prop_map(|(system, v)| {
                Constraint::AromaticSystem(system, AromaticSystemConstraintAst::ElectronCount(v))
            })
            .boxed();
        choices.push(aromatic_leaf);

        if counts.atom > 0 {
            let atom_id = atom_id_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);
            let atoms_vec = prop::collection::vec(atom_id.clone(), 1..=max_atoms);

            let atoms = (system_id.clone(), atoms_vec.clone())
                .prop_map(|(system, atoms)| {
                    Constraint::Relational(RelationalConstraint::AromaticSystemAtoms {
                        system,
                        atoms,
                    })
                })
                .boxed();
            let contains = (system_id.clone(), atom_id)
                .prop_map(|(system, atom)| {
                    Constraint::Relational(RelationalConstraint::AromaticSystemContains {
                        system,
                        atom,
                    })
                })
                .boxed();
            let contains_all = (system_id.clone(), atoms_vec)
                .prop_map(|(system, atoms)| {
                    Constraint::Relational(RelationalConstraint::AromaticSystemContainsAll {
                        system,
                        atoms,
                    })
                })
                .boxed();
            let all_atoms = (system_id.clone(), atom_constraint_strategy())
                .prop_map(|(system, predicate)| {
                    Constraint::Relational(RelationalConstraint::AromaticSystemAllAtoms {
                        system,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            let any_atom = (system_id, atom_constraint_strategy())
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
        let bond_id = multicenter_bond_id_strategy(counts.multicenter);

        let multicenter_leaf = (bond_id.clone(), electron_count_value_strategy(0..=8))
            .prop_map(|(bond, v)| {
                Constraint::MulticenterBond(bond, MulticenterBondConstraintAst::ElectronCount(v))
            })
            .boxed();
        choices.push(multicenter_leaf);

        if counts.atom > 0 {
            let atom_id = atom_id_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);
            let atoms_vec = prop::collection::vec(atom_id.clone(), 1..=max_atoms);

            let atoms = (bond_id.clone(), atoms_vec.clone())
                .prop_map(|(bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondAtoms {
                        bond,
                        atoms,
                    })
                })
                .boxed();
            let contains = (bond_id.clone(), atom_id)
                .prop_map(|(bond, atom)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondContains {
                        bond,
                        atom,
                    })
                })
                .boxed();
            let contains_all = (bond_id.clone(), atoms_vec)
                .prop_map(|(bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondContainsAll {
                        bond,
                        atoms,
                    })
                })
                .boxed();
            let all_atoms = (bond_id.clone(), atom_constraint_strategy())
                .prop_map(|(bond, predicate)| {
                    Constraint::Relational(RelationalConstraint::MulticenterBondAllAtoms {
                        bond,
                        predicate: Box::new(predicate),
                    })
                })
                .boxed();
            let any_atom = (bond_id, atom_constraint_strategy())
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
        let bond_id = noncovalent_bond_id_strategy(counts.noncovalent);
        let atom_id = atom_id_strategy(counts.atom);

        let ends = (bond_id.clone(), atom_id.clone(), atom_id.clone())
            .prop_map(|(bond, a, b)| {
                Constraint::Relational(RelationalConstraint::NoncovalentBondEnds {
                    bond,
                    atoms: [a, b],
                })
            })
            .boxed();
        let contains = (bond_id.clone(), atom_id)
            .prop_map(|(bond, atom)| {
                Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond, atom })
            })
            .boxed();
        let ends_satisfy = (
            bond_id,
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
        let sa_id = stereo_atom_id_strategy(counts.stereo_atom);
        let atom_id = atom_id_strategy(counts.atom);
        let max_atoms = counts.atom.min(3);

        let site = (sa_id.clone(), atom_id.clone())
            .prop_map(|(stereo_atom, atom)| {
                Constraint::Relational(RelationalConstraint::StereoAtomSite { stereo_atom, atom })
            })
            .boxed();
        let contains = (sa_id.clone(), atom_id.clone())
            .prop_map(|(stereo_atom, atom)| {
                Constraint::Relational(RelationalConstraint::StereoAtomContains {
                    stereo_atom,
                    atom,
                })
            })
            .boxed();
        let ligands = (sa_id.clone(), prop::collection::vec(atom_id, 1..=max_atoms))
            .prop_map(|(stereo_atom, atoms)| {
                Constraint::Relational(RelationalConstraint::StereoAtomLigands {
                    stereo_atom,
                    atoms,
                })
            })
            .boxed();
        let all_ligands = (sa_id.clone(), atom_constraint_strategy())
            .prop_map(|(stereo_atom, c)| {
                Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                    stereo_atom,
                    predicate: Box::new(c),
                })
            })
            .boxed();
        let any_ligand = (sa_id, atom_constraint_strategy())
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
        let sb_id = stereo_bond_id_strategy(counts.stereo_bond);

        if counts.bond > 0 {
            let bond_id = bond_id_strategy(counts.bond);
            let site = (sb_id.clone(), bond_id)
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
            let atom_id = atom_id_strategy(counts.atom);
            let max_atoms = counts.atom.min(3);

            let contains = (sb_id.clone(), atom_id.clone())
                .prop_map(|(stereo_bond, atom)| {
                    Constraint::Relational(RelationalConstraint::StereoBondContains {
                        stereo_bond,
                        atom,
                    })
                })
                .boxed();
            let ligands = (sb_id.clone(), prop::collection::vec(atom_id, 1..=max_atoms))
                .prop_map(|(stereo_bond, atoms)| {
                    Constraint::Relational(RelationalConstraint::StereoBondLigands {
                        stereo_bond,
                        atoms,
                    })
                })
                .boxed();
            let all_ligands = (sb_id.clone(), atom_constraint_strategy())
                .prop_map(|(stereo_bond, c)| {
                    Constraint::Relational(RelationalConstraint::StereoBondAllLigands {
                        stereo_bond,
                        predicate: Box::new(c),
                    })
                })
                .boxed();
            let any_ligand = (sb_id, atom_constraint_strategy())
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

pub(crate) fn molecule_ast_with_atom_subset_strategy(
) -> impl Strategy<Value = (MoleculeAst, Vec<AtomId>)> {
    molecule_ast_structurally_unambiguous_strategy().prop_flat_map(|ast| {
        let atom_count = ast.atoms().count();
        (Just(ast), prop::collection::vec(any::<bool>(), atom_count)).prop_map(|(ast, keep)| {
            let atoms = keep
                .into_iter()
                .enumerate()
                .filter_map(|(index, keep)| keep.then_some(AtomId(index as u32)))
                .collect();
            (ast, atoms)
        })
    })
}

pub(crate) fn molecule_ast_with_removals_strategy(
) -> impl Strategy<Value = (MoleculeAst, Vec<AtomId>, Vec<BondId>)> {
    molecule_ast_strategy().prop_flat_map(|ast| {
        let atom_count = ast.atoms().count();
        let bond_count = ast.bonds().count();
        (
            Just(ast),
            prop::collection::vec(any::<bool>(), atom_count),
            prop::collection::vec(any::<bool>(), bond_count),
        )
            .prop_map(|(ast, atom_mask, bond_mask)| {
                let atoms = atom_mask
                    .into_iter()
                    .enumerate()
                    .filter_map(|(index, remove)| remove.then_some(AtomId(index as u32)))
                    .collect();
                let bonds = bond_mask
                    .into_iter()
                    .enumerate()
                    .filter_map(|(index, remove)| remove.then_some(BondId(index as u32)))
                    .collect();
                (ast, atoms, bonds)
            })
    })
}

pub(crate) fn molecule_ast_structurally_unambiguous_strategy() -> impl Strategy<Value = MoleculeAst>
{
    molecule_ast_strategy().prop_filter(
        "entity incidence identifies at most one entity of each family",
        |ast| {
            let atom_images = (0..ast.atoms().count())
                .map(NodeId::from)
                .collect::<Vec<_>>();
            let correspondence = MoleculeCorrespondence::induce(
                ast,
                ast,
                Correspondence::from_images(&atom_images, ast.atoms().count()),
            );
            correspondence.is_total()
                && correspondence_rights_are_unique(correspondence.atoms())
                && correspondence_rights_are_unique(correspondence.bonds())
                && correspondence_rights_are_unique(correspondence.dative_bonds())
                && correspondence_rights_are_unique(correspondence.aromatic_systems())
                && correspondence_rights_are_unique(correspondence.multicenter_bonds())
                && correspondence_rights_are_unique(correspondence.noncovalent_bonds())
                && correspondence_rights_are_unique(correspondence.stereo_atoms())
                && correspondence_rights_are_unique(correspondence.stereo_bonds())
        },
    )
}

fn correspondence_rights_are_unique<Id>(correspondence: &Correspondence<Id>) -> bool
where
    Id: Copy + Eq + Ord + From<usize> + Hash,
{
    correspondence
        .matched_pairs()
        .iter()
        .map(|&(_, right)| right)
        .collect::<HashSet<_>>()
        .len()
        == correspondence.matched_pair_count()
}

/// Generate a `MoleculeMetadata` populated for an AST of the given counts. Entity
/// ids use deterministic prefixed names (`atom0`, `bond1`, ...) so that
/// names are unique across kinds and disjoint from alias names. Atom
/// aliases are capped at 3 and use a 3-element pool (`C`, `N`, `O`) for
/// the alias atom-DSL values, keeping bijectivity (each alias name
/// distinct, each alias atom distinct).
pub(crate) fn metadata_for(counts: ConstraintCounts) -> BoxedStrategy<MoleculeMetadata> {
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
                let mut meta = MoleculeMetadata::new();
                for (i, atom) in atoms.iter().enumerate() {
                    if atom.is_some() {
                        meta.set_keyword(Entity::Atom(AtomId(i as u32)), format!("atom{i}"))
                            .unwrap();
                    }
                }
                for (i, bond) in bonds.iter().enumerate() {
                    if bond.is_some() {
                        meta.set_keyword(Entity::Bond(BondId(i as u32)), format!("bond{i}"))
                            .unwrap();
                    }
                }
                for (i, dative) in datives.iter().enumerate() {
                    if dative.is_some() {
                        meta.set_keyword(
                            Entity::DativeBond(DativeBondId(i as u32)),
                            format!("dative{i}"),
                        )
                        .unwrap();
                    }
                }
                for (i, aromatic) in aromatics.iter().enumerate() {
                    if aromatic.is_some() {
                        meta.set_keyword(
                            Entity::AromaticSystem(AromaticSystemId(i as u32)),
                            format!("aromatic{i}"),
                        )
                        .unwrap();
                    }
                }
                for (i, multicenter) in multicenters.iter().enumerate() {
                    if multicenter.is_some() {
                        meta.set_keyword(
                            Entity::MulticenterBond(MulticenterBondId(i as u32)),
                            format!("multicenter{i}"),
                        )
                        .unwrap();
                    }
                }
                for (i, noncovalent) in noncovalents.iter().enumerate() {
                    if noncovalent.is_some() {
                        meta.set_keyword(
                            Entity::NoncovalentBond(NoncovalentBondId(i as u32)),
                            format!("noncovalent{i}"),
                        )
                        .unwrap();
                    }
                }
                for (i, element) in ALIAS_ELEMENTS.iter().enumerate() {
                    meta.add_atom_alias(format!("al{i}"), AtomAst::from_element(*element))
                        .unwrap();
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
        .map(|id| {
            let element = ELEMENTS[id % ELEMENTS.len()];
            AtomAst::from_element(element)
        })
        .collect()
}

pub(crate) fn transaction_path_bonds(count: usize) -> Vec<AddBond> {
    (0..count.saturating_sub(1))
        .map(|id| {
            Edit::add_bond(
                AtomHandle::New(id),
                AtomHandle::New(id + 1),
                BondAst::from_order((id % 3 + 1) as u8),
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
        .map(|id| {
            (
                AtomId(id as u32),
                AtomId((id + 1) as u32),
                BondAst::from_order((id % 3 + 1) as u8),
            )
        })
        .collect();
    MoleculeAst::from_parts(MoleculeParts {
        atoms,
        bonds,
        ..Default::default()
    })
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
        id: usize,
    },
    RemoveBond {
        count: usize,
        id: usize,
    },
    SetAtomCharge {
        count: usize,
        id: usize,
        charge: i64,
    },
    SetBondOrder {
        count: usize,
        id: usize,
        order: u8,
    },
    AddAtomConstraint {
        count: usize,
        id: usize,
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
            Self::RemoveAtom { count, id } => {
                vec![Edit::remove_atom(AtomHandle::Id(AtomId(
                    (id % count) as u32,
                )))]
            }
            Self::RemoveBond { count, id } => vec![Edit::remove_bond(BondHandle::Id(BondId(
                (id % (count - 1)) as u32,
            )))],
            Self::SetAtomCharge { count, id, charge } => {
                vec![Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId((id % count) as u32)),
                    change: AtomFieldChange::Charge {
                        old: ValueAst::default(),
                        new: ValueAst::Lit(*charge),
                    },
                }]
            }
            Self::SetBondOrder { count, id, order } => {
                let bond_id = id % (count - 1);
                vec![Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(bond_id as u32)),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit((bond_id % 3 + 1) as i64),
                        new: ValueAst::Lit(*order as i64),
                    },
                }]
            }
            Self::AddAtomConstraint { count, id, size } => {
                vec![Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId((id % count) as u32)),
                    old: None,
                    new: Some(AtomConstraintAst::ring_membership(
                        RingScope::Size(*size as u8),
                        1,
                    )),
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
                        AtomHandle::Id(AtomId(donor as u32)),
                        AtomHandle::Id(AtomId(acceptor as u32)),
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
        (1usize..=6, 0usize..6).prop_map(|(count, id)| TransactionCase::RemoveAtom { count, id }),
        (2usize..=6, 0usize..5).prop_map(|(count, id)| TransactionCase::RemoveBond { count, id }),
        (1usize..=6, 0usize..6, -3i64..=3).prop_map(|(count, id, charge)| {
            TransactionCase::SetAtomCharge { count, id, charge }
        }),
        (2usize..=6, 0usize..5, 1u8..=3)
            .prop_map(|(count, id, order)| { TransactionCase::SetBondOrder { count, id, order } }),
        (1usize..=6, 0usize..6, 3i64..=8).prop_map(|(count, id, size)| {
            TransactionCase::AddAtomConstraint { count, id, size }
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

// Overlay-transaction fixture: a fixed 6-carbon path carrying two overlays of each DAMN kind, so
// the edit generator can remove ≥2 of one kind (the batched-remove path) and mix overlay removes
// with atom appends and a topology removal in one transaction. Atom sets are the single source of
// truth for both `overlay_transaction_base` and the remove edits, so the `OldState` check matches.
const AROMATIC_SETS: [&[u32]; 2] = [&[0, 1, 2], &[3, 4, 5]];
const MULTICENTER_SETS: [&[u32]; 2] = [&[1, 2, 3], &[0, 4, 5]];
const DATIVE_DONORS: [&[u32]; 2] = [&[0], &[4]];
const DATIVE_ACCEPTORS: [u32; 2] = [1, 5];
const NONCOVALENT_PAIRS: [[u32; 2]; 2] = [[0, 2], [3, 5]];

pub(crate) fn overlay_transaction_base() -> MoleculeAst {
    let atoms: Vec<AtomAst> = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
    let bonds = (0..5)
        .map(|i| (AtomId(i), AtomId(i + 1), BondAst::from_order(1)))
        .collect();
    let dative = (0..2)
        .map(|i| {
            (
                DATIVE_DONORS[i]
                    .iter()
                    .map(|&a| AtomId(a))
                    .collect::<Vec<_>>(),
                AtomId(DATIVE_ACCEPTORS[i]),
                DativeBondAst::from_order(1),
            )
        })
        .collect();
    let aromatic = AROMATIC_SETS
        .iter()
        .map(|set| {
            (
                set.iter().map(|&a| AtomId(a)).collect::<Vec<_>>(),
                AromaticSystemAst::default(),
            )
        })
        .collect();
    let multicenter = MULTICENTER_SETS
        .iter()
        .map(|set| {
            (
                set.iter().map(|&a| AtomId(a)).collect::<Vec<_>>(),
                MulticenterBondAst::default(),
            )
        })
        .collect();
    let noncovalent = NONCOVALENT_PAIRS
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                AtomId(*b),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )
        })
        .collect();
    MoleculeAst::from_parts(MoleculeParts {
        atoms,
        bonds,
        dative,
        aromatic,
        multicenter,
        noncovalent,
        constraints: Constraints::new(),
        ..Default::default()
    })
}

/// A valid transaction over `overlay_transaction_base`: optional atom appends, then one batched
/// remove edit per DAMN kind (a chosen subset of that kind's two ids — exercises ≥2 same-kind
/// removal), then a topology removal of a chosen atom subset (cascade-removes overlays not removed
/// explicitly). Ordered adds → overlay removes → topology, mirroring `apply_at`; every id resolves
/// against the pre-removal base state, so `transact` succeeds and the round-trip properties apply.
pub(crate) fn overlay_transaction_strategy() -> impl Strategy<Value = (MoleculeAst, Vec<Edit>)> {
    (
        prop::collection::vec(any::<bool>(), 2),
        prop::collection::vec(any::<bool>(), 2),
        prop::collection::vec(any::<bool>(), 2),
        prop::collection::vec(any::<bool>(), 2),
        0usize..=2,
        prop::collection::vec(any::<bool>(), 6),
        prop::collection::vec(any::<bool>(), 6),
        prop::collection::vec(any::<bool>(), 2),
        prop::collection::vec(any::<bool>(), 2),
    )
        .prop_map(
            |(rm_ar, rm_mc, rm_dv, rm_nc, add, mod_at, rm_at, con_ar, con_mc)| {
                let mut edits: Vec<Edit> = Vec::new();
                if add > 0 {
                    edits.push(Edit::AddAtoms {
                        atoms: (0..add)
                            .map(|_| AtomAst::from_element(Element::C))
                            .collect(),
                    });
                }
                // Base carbons carry the default (`Undetermined`) charge, so that is the `old` value.
                for i in (0..6).filter(|&i| mod_at[i]) {
                    edits.push(Edit::ModifyAtomField {
                        id: AtomHandle::Id(AtomId(i as u32)),
                        change: AtomFieldChange::Charge {
                            old: ValueAst::default(),
                            new: ValueAst::Lit(1),
                        },
                    });
                }
                // Molecule-level constraints referencing overlays, added before the removals so that
                // removing a referenced overlay exercises constraint drop/remap + its rollback restore.
                for i in (0..2).filter(|&i| con_ar[i]) {
                    edits.push(Edit::AddMoleculeConstraint {
                        constraint: Constraint::AromaticSystem(
                            AromaticSystemId(i as u32),
                            AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6)),
                        ),
                    });
                }
                for i in (0..2).filter(|&i| con_mc[i]) {
                    edits.push(Edit::AddMoleculeConstraint {
                        constraint: Constraint::MulticenterBond(
                            MulticenterBondId(i as u32),
                            MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(4)),
                        ),
                    });
                }
                let dative: Vec<_> = (0..2)
                    .filter(|&i| rm_dv[i])
                    .map(|i| {
                        let mut atoms: Vec<AtomHandle> = DATIVE_DONORS[i]
                            .iter()
                            .map(|&a| AtomHandle::Id(AtomId(a)))
                            .collect();
                        atoms.push(AtomHandle::Id(AtomId(DATIVE_ACCEPTORS[i])));
                        (
                            DativeBondHandle::Id(DativeBondId(i as u32)),
                            atoms,
                            DativeBondAst::from_order(1),
                        )
                    })
                    .collect();
                if !dative.is_empty() {
                    edits.push(Edit::RemoveDativeBonds { removes: dative });
                }
                let aromatic: Vec<_> = (0..2)
                    .filter(|&i| rm_ar[i])
                    .map(|i| {
                        let atoms = AROMATIC_SETS[i]
                            .iter()
                            .map(|&a| AtomHandle::Id(AtomId(a)))
                            .collect();
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(i as u32)),
                            atoms,
                            AromaticSystemAst::default(),
                        )
                    })
                    .collect();
                if !aromatic.is_empty() {
                    edits.push(Edit::RemoveAromaticSystems { removes: aromatic });
                }
                let multicenter: Vec<_> = (0..2)
                    .filter(|&i| rm_mc[i])
                    .map(|i| {
                        let atoms = MULTICENTER_SETS[i]
                            .iter()
                            .map(|&a| AtomHandle::Id(AtomId(a)))
                            .collect();
                        (
                            MulticenterBondHandle::Id(MulticenterBondId(i as u32)),
                            atoms,
                            MulticenterBondAst::default(),
                        )
                    })
                    .collect();
                if !multicenter.is_empty() {
                    edits.push(Edit::RemoveMulticenterBonds {
                        removes: multicenter,
                    });
                }
                let noncovalent: Vec<_> = (0..2)
                    .filter(|&i| rm_nc[i])
                    .map(|i| {
                        (
                            NoncovalentBondHandle::Id(NoncovalentBondId(i as u32)),
                            [
                                AtomHandle::Id(AtomId(NONCOVALENT_PAIRS[i][0])),
                                AtomHandle::Id(AtomId(NONCOVALENT_PAIRS[i][1])),
                            ],
                            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                        )
                    })
                    .collect();
                if !noncovalent.is_empty() {
                    edits.push(Edit::RemoveNoncovalentBonds {
                        removes: noncovalent,
                    });
                }
                let atoms: Vec<AtomHandle> = (0..6)
                    .filter(|&i| rm_at[i])
                    .map(|i| AtomHandle::Id(AtomId(i as u32)))
                    .collect();
                if !atoms.is_empty() {
                    edits.push(Edit::RemoveTopology {
                        atoms,
                        bonds: vec![],
                    });
                }
                (overlay_transaction_base(), edits)
            },
        )
}

/// `(base, edits)` pairs for the transact round-trip properties: the single-edit `TransactionCase`
/// coverage plus the multi-edit overlay-removal sequences.
pub(crate) fn transaction_edits_strategy() -> impl Strategy<Value = (MoleculeAst, Vec<Edit>)> {
    prop_oneof![
        transaction_case_strategy().prop_map(|case| (case.base(), case.edits())),
        overlay_transaction_strategy(),
    ]
}

/// A small localized molecule: 1–4 element atoms over a simple edge set, bond orders 1–3.
fn simple_molecule_strategy() -> impl Strategy<Value = MoleculeAst> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                prop::collection::vec(
                    element_strategy().prop_map(AtomAst::from_element),
                    atom_count,
                ),
                edge_set_strategy(atom_count),
            )
        })
        .prop_flat_map(|(atoms, edges)| {
            let orders = prop::collection::vec(1u8..=3, edges.len());
            (Just(atoms), Just(edges), orders)
        })
        .prop_map(|(atoms, edges, orders)| {
            let bonds = edges
                .iter()
                .zip(orders)
                .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondAst::from_order(order)))
                .collect();
            MoleculeAst::from_parts(MoleculeParts {
                atoms,
                bonds,
                ..Default::default()
            })
        })
}

pub(crate) fn reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(simple_molecule_strategy())
}

pub(crate) fn replacement_reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    (molecule_ast_strategy(), molecule_ast_strategy()).prop_map(|(lhs, rhs)| {
        let correspondence =
            Correspondence::new(Vec::new(), lhs.atoms().count(), rhs.atoms().count());
        ReactionAst::from_sides(lhs, rhs, correspondence)
    })
}

pub(crate) fn comprehensive_reaction_strategy() -> BoxedStrategy<ReactionAst> {
    prop_oneof![
        2 => overlay_reaction_strategy(),
        1 => replacement_reaction_strategy(),
    ]
    .boxed()
}

/// A localized molecule with DAMN overlays (dative / aromatic / multicenter / noncovalent) plus
/// stereo (tetrahedral atoms / cis-trans bonds) and no molecule constraints (orthogonal). 1–4 atoms;
/// overlays generated as in `molecule_ast_strategy`, scoped.
fn overlay_molecule_strategy() -> impl Strategy<Value = MoleculeAst> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                Just(atom_count),
                prop::collection::vec(
                    element_strategy().prop_map(AtomAst::from_element),
                    atom_count,
                ),
                edge_set_strategy(atom_count),
            )
        })
        .prop_flat_map(|(atom_count, atoms, edges)| {
            let orders = prop::collection::vec(1u8..=3, edges.len());
            let datives = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    dative_bond_strategy(),
                ),
                0..=1,
            );
            let aromatics = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), aromatic_system_ast_for(n))
                    },
                ),
                0..=1,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_ast_for(n))
                    },
                ),
                0..=1,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_ast_strategy(),
                ),
                0..=1,
            );
            // A tetrahedral stereo atom: a site atom plus four ligands. Real atoms fill the first
            // slots (ids need not be graph neighbors — tier-1 only requires the kind's ligand
            // count); virtual implicit-H / lone-pair fills pad to `degree == 4`, all bearing the
            // site atom. 0..=1 so many molecules have none.
            let stereo_atoms = stereo_atom_overlay_strategy(atom_count);
            // A cis/trans stereo bond: a bond as site plus two ligand atoms (padded with virtual
            // fills to `degree == 4`). Requires a bond to name as site.
            let stereo_bonds = if edges.is_empty() {
                Just(Vec::new()).boxed()
            } else {
                stereo_bond_overlay_strategy(atom_count, edges.len())
            };
            (
                Just(atoms),
                Just(edges),
                orders,
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
                edges,
                orders,
                datives,
                aromatics,
                multicenters,
                noncovalents,
                stereo_atoms,
                stereo_bonds,
            )| {
                let bonds = edges
                    .iter()
                    .zip(orders)
                    .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondAst::from_order(order)))
                    .collect();
                let dative = datives
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((vec![*a], *b, data)),
                        _ => None,
                    })
                    .collect();
                let aromatic = aromatics
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let multicenter = multicenters
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let noncovalent = noncovalents
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((*a, *b, data)),
                        _ => None,
                    })
                    .collect();
                MoleculeAst::from_parts(MoleculeParts {
                    atoms,
                    bonds,
                    dative,
                    aromatic,
                    multicenter,
                    noncovalent,
                    stereo_atoms,
                    stereo_bonds,
                    constraints: Constraints::new(),
                })
            },
        )
}

/// Cosets valid for `kind`: `Undetermined` or an in-range `Lit` index (`0..kind.count()`). Relative
/// reaction ops (`Swap` / `Mirror` / `Apply`) act on the coset through the kind's algebra, which
/// panics on an out-of-range index — so unlike the generic `stereo_coset_strategy`, indices are
/// bounded by the kind's coset count.
pub(crate) fn stereo_coset_for_kind(kind: StereoKind) -> impl Strategy<Value = StereoCoset> {
    let count = kind.count() as u32;
    prop_oneof![
        Just(StereoCoset::Undetermined),
        (0..count).prop_map(StereoCoset::Lit),
    ]
}

pub(crate) fn aromatic_system_update_for(
    atom_count: usize,
) -> impl Strategy<Value = AromaticSystemUpdate> {
    (
        prop::option::of(prop_oneof![
            Just(ElectronCountsAst::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsAst::Lit),
        ]),
        prop::option::of(value_basic(-2..=2)),
        spin_state_update_strategy(),
        aromatic_system_update_constraints_strategy(),
    )
        .prop_map(
            |(electrons, charge, spin, constraints)| AromaticSystemUpdate {
                electrons,
                charge,
                spin,
                constraints,
            },
        )
}

pub(crate) fn multicenter_bond_update_for(
    atom_count: usize,
) -> impl Strategy<Value = MulticenterBondUpdate> {
    (
        prop::option::of(prop_oneof![
            Just(ElectronCountsAst::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsAst::Lit),
        ]),
        prop::option::of(value_basic(-2..=2)),
        spin_state_update_strategy(),
        multicenter_bond_update_constraints_strategy(),
    )
        .prop_map(
            |(electrons, charge, spin, constraints)| MulticenterBondUpdate {
                electrons,
                charge,
                spin,
                constraints,
            },
        )
}

pub(crate) fn stereo_atom_application_update_strategy() -> impl Strategy<Value = StereoAtomUpdate> {
    (
        prop_oneof![
            Just(StereoConfigurationUpdate::Unchanged),
            Just(StereoConfigurationUpdate::Undetermined),
            prop::option::of(stereo_coset_for_kind(StereoKind::Tetrahedral)).prop_map(|coset| {
                StereoConfigurationUpdate::Kinded {
                    kind: StereoKind::Tetrahedral,
                    coset,
                }
            }),
        ],
        stereo_atom_update_constraints_strategy(StereoKind::Tetrahedral),
    )
        .prop_map(|(configuration, constraints)| StereoAtomUpdate {
            configuration,
            constraints,
        })
}

pub(crate) fn stereo_bond_application_update_strategy() -> impl Strategy<Value = StereoBondUpdate> {
    (
        prop_oneof![
            Just(StereoConfigurationUpdate::Unchanged),
            Just(StereoConfigurationUpdate::Undetermined),
            prop::option::of(stereo_coset_for_kind(StereoKind::CisTrans)).prop_map(|coset| {
                StereoConfigurationUpdate::Kinded {
                    kind: StereoKind::CisTrans,
                    coset,
                }
            }),
        ],
        stereo_bond_update_constraints_strategy(StereoKind::CisTrans),
    )
        .prop_map(|(configuration, constraints)| StereoBondUpdate {
            configuration,
            constraints,
        })
}

/// A `degree`-length ligand frame of *distinct* `StereoLigand`s over `atom_count` atoms. The overlay
/// matcher (`permutation_for_ligands`) rejects a non-unique frame, so `apply` finds no identity
/// match — hence ligands must be unique. Real-atom ligands come first (distinct atoms); virtual
/// implicit-H / lone-pair fills pad by distinct `(atom, kind)` pairs. A frame of `degree` unique
/// ligands needs `atom_count * 3 >= degree`, so callers gate on `atom_count`.
fn unique_ligand_frame(
    atom_count: usize,
    degree: usize,
) -> impl Strategy<Value = Vec<StereoLigand>> {
    let pool: Vec<StereoLigand> = (0..atom_count as u32)
        .flat_map(|a| {
            [
                StereoLigandKind::Atom,
                StereoLigandKind::ImplicitHydrogen,
                StereoLigandKind::LonePair,
            ]
            .into_iter()
            .map(move |kind| StereoLigand::new(AtomId(a), kind))
        })
        .collect();
    Just(pool).prop_shuffle().prop_map(move |mut pool| {
        pool.truncate(degree);
        pool
    })
}

/// 0..=1 tetrahedral stereo atoms over an `atom_count`-atom molecule (needs `atom_count >= 2` for a
/// `degree`-length unique ligand frame). Site is any atom; ligands are distinct real/virtual ligands
/// whose atoms need not be graph neighbors (tier-1 requires only the ligand count for the kind).
fn stereo_atom_overlay_strategy(
    atom_count: usize,
) -> BoxedStrategy<Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>> {
    let degree = StereoKind::Tetrahedral.degree();
    if atom_count * 3 < degree {
        return Just(Vec::new()).boxed();
    }
    prop::collection::vec(
        (
            0..atom_count as u32,
            unique_ligand_frame(atom_count, degree),
            stereo_coset_for_kind(StereoKind::Tetrahedral),
        ),
        0..=1,
    )
    .prop_map(move |entries| {
        entries
            .into_iter()
            .map(|(site, ligands, coset)| {
                let ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
                (AtomId(site), ligands, ast)
            })
            .collect()
    })
    .boxed()
}

/// 0..=1 cis/trans stereo bonds (needs `atom_count >= 2` for a `degree`-length unique frame). Site is
/// any bond; ligands are distinct real/virtual ligands (their atoms need not be double-bond termini).
fn stereo_bond_overlay_strategy(
    atom_count: usize,
    bond_count: usize,
) -> BoxedStrategy<Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>> {
    let degree = StereoKind::CisTrans.degree();
    if bond_count == 0 || atom_count * 3 < degree {
        return Just(Vec::new()).boxed();
    }
    prop::collection::vec(
        (
            0..bond_count as u32,
            unique_ligand_frame(atom_count, degree),
            stereo_coset_for_kind(StereoKind::CisTrans),
        ),
        0..=1,
    )
    .prop_map(move |entries| {
        entries
            .into_iter()
            .map(|(site, ligands, coset)| {
                let ast = StereoBondAst::new(StereoKind::CisTrans, coset);
                (BondId(site), ligands, ast)
            })
            .collect()
    })
    .boxed()
}

/// A reaction whose `lhs` carries DAMN overlays — exercises overlay carry, correspondence, and
/// co-deletion through compose.
pub(crate) fn overlay_reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(overlay_molecule_strategy())
}

/// An optional edit to a surviving stereo entity. The relative ops (`Swap` / `Mirror` / `Apply`)
/// resolve `old` from the matched host coset at apply, carrying no pre-state; `SetCoset` becomes a
/// `ModifyField { Configuration }` whose `old` is read from `lhs`, so apply's precondition holds.
#[derive(Clone, Debug)]
enum StereoOp {
    Swap,
    Mirror,
    Apply(Permutation),
    SetCoset(StereoCoset),
}

/// Per-surviving-stereo-entity optional op: `Swap` / `Mirror` use the kind's in-group generators,
/// `SetCoset` is bounded to the kind's in-range cosets, and `Apply` a permutation in the kind's
/// parent group. The coset algebra rejects out-of-group permutations (`reindex` → `None`, which
/// `act` unwraps), so `Apply` only draws arbitrary permutations for kinds whose parent is the full
/// symmetric group (Tetrahedral); other kinds omit it and lean on the in-group `Swap` / `Mirror`.
fn stereo_op_strategy(kind: StereoKind) -> impl Strategy<Value = Option<StereoOp>> {
    let base = prop_oneof![
        Just(StereoOp::Swap),
        Just(StereoOp::Mirror),
        stereo_coset_for_kind(kind).prop_map(StereoOp::SetCoset),
    ]
    .boxed();
    let ops = if kind == StereoKind::Tetrahedral {
        prop_oneof![
            base,
            permutation_strategy(kind.degree()).prop_map(StereoOp::Apply),
        ]
        .boxed()
    } else {
        base
    };
    prop::option::weighted(0.5, ops)
}

/// A valid reaction over any generated `lhs`: DPO-valid atom deletions (each removed atom takes its
/// incident bonds, overlays, and stereo entities), per-surviving-entity optional field / relative
/// edits (the absolute `old` read from `lhs`, so apply's precondition holds), plus up to two new
/// atoms bonded to the lowest survivor. No dangling by construction.
fn reaction_over(
    molecule: impl Strategy<Value = MoleculeAst>,
) -> impl Strategy<Value = ReactionAst> {
    molecule
        .prop_flat_map(|lhs| {
            let atom_count = lhs.atoms().count();
            let bond_count = lhs.bonds().count();
            let dative_count = lhs.dative_bonds().count();
            let aromatic_count = lhs.aromatic_systems().count();
            let multicenter_count = lhs.multicenter_bonds().count();
            let stereo_atom_count = lhs.stereo_atoms().count();
            let stereo_bond_count = lhs.stereo_bonds().count();
            (
                Just(lhs),
                prop::collection::vec(weighted(0.25), atom_count),
                prop::collection::vec(prop::option::of(-2i64..=2), atom_count),
                prop::collection::vec(prop::option::of(1i64..=3), bond_count),
                prop::collection::vec(element_strategy(), 0..=2),
                (
                    // Overlay `ModifyField` on survivors: dative order, aromatic / multicenter charge.
                    prop::collection::vec(prop::option::of(1i64..=3), dative_count),
                    prop::collection::vec(prop::option::of(-2i64..=2), aromatic_count),
                    prop::collection::vec(prop::option::of(-2i64..=2), multicenter_count),
                    // Add an `Aromatic` constraint to a surviving dative (guarded on absence).
                    prop::collection::vec(weighted(0.3), dative_count),
                    // Add a noncovalent overlay between the two newly-added atoms.
                    weighted(0.4),
                ),
                (
                    prop::collection::vec(
                        stereo_op_strategy(StereoKind::Tetrahedral),
                        stereo_atom_count,
                    ),
                    prop::collection::vec(
                        stereo_op_strategy(StereoKind::CisTrans),
                        stereo_bond_count,
                    ),
                ),
            )
        })
        .prop_map(
            |(lhs, removals, charges, orders, additions, overlay_ops, stereo_ops)| {
                build_reaction(
                    lhs,
                    removals,
                    charges,
                    orders,
                    additions,
                    overlay_ops,
                    stereo_ops,
                )
            },
        )
}

/// Per-entity overlay `ModifyField` / `Add` / `ModifyConstraint` randomness: dative orders, aromatic
/// charges, multicenter charges, dative-Aromatic-constraint flags, and the add-noncovalent flag.
type OverlayOps = (
    Vec<Option<i64>>,
    Vec<Option<i64>>,
    Vec<Option<i64>>,
    Vec<bool>,
    bool,
);

/// Per-stereo-entity optional op: stereo atoms, then stereo bonds.
type StereoOps = (Vec<Option<StereoOp>>, Vec<Option<StereoOp>>);

fn build_reaction(
    lhs: MoleculeAst,
    removals: Vec<bool>,
    charges: Vec<Option<i64>>,
    orders: Vec<Option<i64>>,
    additions: Vec<Element>,
    overlay_ops: OverlayOps,
    stereo_ops: StereoOps,
) -> ReactionAst {
    let atom_count = lhs.atoms().count();
    let bond_count = lhs.bonds().count();
    let removed_atoms: HashSet<AtomId> = removals
        .iter()
        .enumerate()
        .filter(|&(_, &remove)| remove)
        .map(|(index, _)| AtomId(index as u32))
        .collect();
    // A removed atom takes all its incident bonds with it (DPO-valid; apply never dangles).
    let mut removed_bonds: HashSet<BondId> = HashSet::new();
    for j in 0..bond_count as u32 {
        let [x, y] = lhs.raw_graph().edge_endpoints(EdgeId(j));
        if removed_atoms.contains(&AtomId::from(x)) || removed_atoms.contains(&AtomId::from(y)) {
            removed_bonds.insert(BondId(j));
        }
    }

    let mut deltas: Vec<Delta> = Vec::new();
    for &id in &removed_atoms {
        deltas.push(Delta::Atom(AtomDelta::Remove {
            id,
            ast: lhs.atom(id).ast.clone(),
        }));
    }
    for &id in &removed_bonds {
        let [x, y] = lhs.raw_graph().edge_endpoints(EdgeId(id.0));
        deltas.push(Delta::Bond(BondDelta::Remove {
            id,
            atoms: [AtomId::from(x), AtomId::from(y)],
            ast: lhs.bond(id).ast.clone(),
        }));
    }
    // A removed atom also takes its incident overlays (DPO-valid; apply never dangles on overlays).
    let mut removed_dative: HashSet<DativeBondId> = HashSet::new();
    let mut removed_aromatic: HashSet<AromaticSystemId> = HashSet::new();
    let mut removed_multicenter: HashSet<MulticenterBondId> = HashSet::new();
    let mut removed_noncovalent: HashSet<NoncovalentBondId> = HashSet::new();
    for &id in &removed_atoms {
        let view = lhs.atom(id);
        removed_dative.extend(view.dative_bond_ids());
        if let Some(system) = view.aromatic_system_id() {
            removed_aromatic.insert(system);
        }
        removed_multicenter.extend(view.multicenter_bond_ids());
        removed_noncovalent.extend(view.noncovalent_bond_ids());
    }
    for &id in &removed_dative {
        let view = lhs.dative_bond(id);
        deltas.push(Delta::DativeBond(DativeBondDelta::Remove {
            id,
            donors: view.donor_ids().collect(),
            acceptor: view.acceptor_id(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_aromatic {
        let view = lhs.aromatic_system(id);
        deltas.push(Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_multicenter {
        let view = lhs.multicenter_bond(id);
        deltas.push(Delta::MulticenterBond(MulticenterBondDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_noncovalent {
        let view = lhs.noncovalent_bond(id);
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
            id,
            atoms: view.atom_ids(),
            ast: view.ast.clone(),
        }));
    }
    // A removed atom also takes its incident stereo entities (site OR ligand incidence), else
    // apply / span / DpoValidator dangle. `incident_ids` covers both.
    let mut removed_stereo_atom: HashSet<StereoAtomId> = HashSet::new();
    let mut removed_stereo_bond: HashSet<StereoBondId> = HashSet::new();
    for &id in &removed_atoms {
        removed_stereo_atom.extend(lhs.stereo_atoms().incident_ids(id));
        removed_stereo_bond.extend(lhs.stereo_bonds().incident_ids(id));
    }
    for &id in &removed_stereo_atom {
        let view = lhs.stereo_atom(id);
        deltas.push(Delta::StereoAtom(StereoAtomDelta::Remove {
            id,
            site: view.site_id(),
            ligands: view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_stereo_bond {
        let view = lhs.stereo_bond(id);
        deltas.push(Delta::StereoBond(StereoBondDelta::Remove {
            id,
            site: view.site_id(),
            ligands: view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect(),
            ast: view.ast.clone(),
        }));
    }
    for (index, new_charge) in charges.into_iter().enumerate() {
        let id = AtomId(index as u32);
        if removed_atoms.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.atom(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::Atom(AtomDelta::ModifyField {
                id,
                change: AtomFieldChange::Charge { old, new },
            }));
        }
    }
    for (index, new_order) in orders.into_iter().enumerate() {
        let id = BondId(index as u32);
        if removed_bonds.contains(&id) {
            continue;
        }
        let Some(order) = new_order else { continue };
        let old = lhs.bond(id).ast.order.clone();
        let new = ValueAst::Lit(order);
        if old != new {
            deltas.push(Delta::Bond(BondDelta::ModifyField {
                id,
                change: BondFieldChange::Order { old, new },
            }));
        }
    }
    // Part A — overlay `ModifyField` on survivors: read `old` from `lhs`, emit only when it changes.
    let (
        dative_orders,
        aromatic_charges,
        multicenter_charges,
        dative_aromatic_flags,
        add_noncovalent,
    ) = overlay_ops;
    for (index, new_order) in dative_orders.into_iter().enumerate() {
        let id = DativeBondId(index as u32);
        if removed_dative.contains(&id) {
            continue;
        }
        let Some(order) = new_order else { continue };
        let old = lhs.dative_bond(id).ast.order.clone();
        let new = ValueAst::Lit(order);
        if old != new {
            deltas.push(Delta::DativeBond(DativeBondDelta::ModifyField {
                id,
                change: DativeBondFieldChange::Order { old, new },
            }));
        }
    }
    for (index, new_charge) in aromatic_charges.into_iter().enumerate() {
        let id = AromaticSystemId(index as u32);
        if removed_aromatic.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.aromatic_system(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                id,
                change: AromaticSystemFieldChange::Charge { old, new },
            }));
        }
    }
    for (index, new_charge) in multicenter_charges.into_iter().enumerate() {
        let id = MulticenterBondId(index as u32);
        if removed_multicenter.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.multicenter_bond(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::MulticenterBond(MulticenterBondDelta::ModifyField {
                id,
                change: MulticenterBondFieldChange::Charge { old, new },
            }));
        }
    }
    // Part A — add an `Aromatic` constraint to a surviving dative, guarded on its absence (apply's
    // `old: None` precondition requires no existing constraint under that key).
    for (index, add) in dative_aromatic_flags.into_iter().enumerate() {
        let id = DativeBondId(index as u32);
        if !add || removed_dative.contains(&id) {
            continue;
        }
        let has_aromatic = lhs
            .dative_bond(id)
            .ast
            .constraints
            .iter()
            .any(|c| matches!(c, DativeBondConstraintAst::Aromatic(_)));
        if has_aromatic {
            continue;
        }
        deltas.push(Delta::DativeBond(DativeBondDelta::ModifyConstraint {
            id,
            old: None,
            new: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        }));
    }
    // Part B — stereo edits on survivors. Relative ops resolve `old` from the host at apply;
    // `SetCoset` reads `old` from `lhs`. Every op is emitted only when it *changes* the entity's
    // configuration value: a value no-op (a relative op on an `Undetermined` coset, a `Mirror` on an
    // achiral kind, a stabilizer permutation, or a `SetCoset` to the current value) would materialize
    // a spurious `Modified { X, X }` span state that `to_reaction` diffs back to empty, breaking the
    // span roundtrip.
    let (stereo_atom_ops, stereo_bond_ops) = stereo_ops;
    for (index, op) in stereo_atom_ops.into_iter().enumerate() {
        let id = StereoAtomId(index as u32);
        if removed_stereo_atom.contains(&id) {
            continue;
        }
        let Some(op) = op else { continue };
        let kind = lhs.stereo_atom(id).kind();
        let old = lhs.stereo_atom(id).ast.configuration.clone();
        let (new, delta) = match &op {
            StereoOp::Swap => (old.swap(), StereoAtomDelta::Swap { id, kind }),
            StereoOp::Mirror => (old.mirror(), StereoAtomDelta::Mirror { id, kind }),
            StereoOp::Apply(permutation) => (
                old.apply(*permutation),
                StereoAtomDelta::Apply {
                    id,
                    kind,
                    permutation: *permutation,
                },
            ),
            StereoOp::SetCoset(coset) => {
                let new = StereoConfigurationAst::kinded(kind, coset.clone());
                (
                    new.clone(),
                    StereoAtomDelta::ModifyField {
                        id,
                        change: StereoAtomFieldChange::Configuration {
                            old: old.clone(),
                            new,
                        },
                    },
                )
            }
        };
        if new != old {
            deltas.push(Delta::StereoAtom(delta));
        }
    }
    for (index, op) in stereo_bond_ops.into_iter().enumerate() {
        let id = StereoBondId(index as u32);
        if removed_stereo_bond.contains(&id) {
            continue;
        }
        let Some(op) = op else { continue };
        let kind = lhs.stereo_bond(id).kind();
        let old = lhs.stereo_bond(id).ast.configuration.clone();
        let (new, delta) = match &op {
            StereoOp::Swap => (old.swap(), StereoBondDelta::Swap { id, kind }),
            StereoOp::Mirror => (old.mirror(), StereoBondDelta::Mirror { id, kind }),
            StereoOp::Apply(permutation) => (
                old.apply(*permutation),
                StereoBondDelta::Apply {
                    id,
                    kind,
                    permutation: *permutation,
                },
            ),
            StereoOp::SetCoset(coset) => {
                let new = StereoConfigurationAst::kinded(kind, coset.clone());
                (
                    new.clone(),
                    StereoBondDelta::ModifyField {
                        id,
                        change: StereoBondFieldChange::Configuration {
                            old: old.clone(),
                            new,
                        },
                    },
                )
            }
        };
        if new != old {
            deltas.push(Delta::StereoBond(delta));
        }
    }
    // Append atoms bonded to the lowest surviving atom (isolated if every atom is removed).
    let anchor = (0..atom_count as u32)
        .map(AtomId)
        .find(|id| !removed_atoms.contains(id));
    let mut added_atom_ids: Vec<AtomId> = Vec::new();
    for (offset, element) in additions.into_iter().enumerate() {
        let atom = AtomId((atom_count + offset) as u32);
        added_atom_ids.push(atom);
        deltas.push(Delta::Atom(AtomDelta::Add {
            id: atom,
            ast: AtomAst::from_element(element),
        }));
        if let Some(anchor) = anchor {
            deltas.push(Delta::Bond(BondDelta::Add {
                id: BondId((bond_count + offset) as u32),
                atoms: [anchor, atom],
                ast: BondAst::from_order(1),
            }));
        }
    }
    // Part A — overlay `Add`: a noncovalent bond between the two newly-added atoms (both created in
    // this reaction, so no dangling). Ids append past the lhs noncovalent count.
    if add_noncovalent && added_atom_ids.len() >= 2 {
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(lhs.noncovalent_bonds().count() as u32),
            atoms: [added_atom_ids[0], added_atom_ids[1]],
            ast: NoncovalentBondAst {
                kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::VanDerWaals),
                constraints: Default::default(),
            },
        }));
    }
    ReactionAst::new(lhs, Deltas::from_iter(deltas))
}
