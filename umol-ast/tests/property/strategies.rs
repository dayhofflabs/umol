//! Shared proptest generators for the umol-ast property suite. Domain imports
//! are re-exported (`pub(crate) use`) so the per-area test modules need only
//! `use proptest::prelude::*; use crate::strategies::*;`.

pub(crate) use std::collections::HashSet;
use std::collections::{BTreeSet, HashMap};
pub(crate) use std::fmt::Debug;
use std::hash::Hash;
pub(crate) use std::iter::repeat_with;
pub(crate) use std::ops::RangeInclusive;

use proptest::bool::weighted;
use proptest::prelude::*;
pub(crate) use umol_ast::ast::{
    aromatic_covalence, AddBond, AromaticSystemAst, AromaticSystemConstraintAst,
    AromaticSystemConstraintKey, AromaticSystemConstraintsAst, AromaticSystemDelta,
    AromaticSystemFieldChange, AromaticSystemHandle, AromaticSystemId, AromaticSystemUpdate,
    AromaticValence, AromaticValenceAst, AsLit, AtomAst, AtomConstraintAst, AtomConstraintKey,
    AtomConstraintsAst, AtomDelta, AtomFieldChange, AtomHandle, AtomId, AtomUpdate, BondAst,
    BondConstraintAst, BondConstraintKey, BondConstraintsAst, BondDelta, BondFieldChange,
    BondHandle, BondId, BondUpdate, BooleanAst, Canonicalize, CisTransStereoAst, Constraint,
    ConstraintEdit, Constraints, DativeBondAst, DativeBondConstraintAst, DativeBondConstraintKey,
    DativeBondConstraintsAst, DativeBondDelta, DativeBondFieldChange, DativeBondHandle,
    DativeBondId, DativeBondUpdate, Delta, Deltas, DpoValidator, Edit, Edits, ElectronCountsAst,
    ElementAst, Entity, EntityHandle, EntityKind, FluxionalityAst, FromAst, IntoAst,
    IsotopeMassAst, Lattice, LigandPermutation, LigandSymmetryAst, MemOp, MoleculeAst,
    MoleculeConstraint, MoleculeCorrespondence, MoleculeEntries, MulticenterBondAst,
    MulticenterBondConstraintAst, MulticenterBondConstraintKey, MulticenterBondConstraintsAst,
    MulticenterBondDelta, MulticenterBondFieldChange, MulticenterBondHandle, MulticenterBondId,
    MulticenterBondUpdate, MulticenterValenceAst, NoncovalentBondAst, NoncovalentBondConstraintAst,
    NoncovalentBondConstraintsAst, NoncovalentBondDelta, NoncovalentBondFieldChange,
    NoncovalentBondHandle, NoncovalentBondId, NoncovalentBondKind, NoncovalentBondKindAst,
    NoncovalentBondUpdate, OrientedLigandPermutation, ReactionAst, ReactionSpanAst, RelOp,
    RelationalConstraint, RingMembershipAst, RingScope, StereoAtomAst, StereoAtomConstraintAst,
    StereoAtomConstraintsAst, StereoAtomDelta, StereoAtomFieldChange, StereoAtomHandle,
    StereoAtomId, StereoAtomUpdate, StereoBondAst, StereoBondConstraintAst,
    StereoBondConstraintsAst, StereoBondDelta, StereoBondFieldChange, StereoBondHandle,
    StereoBondId, StereoBondUpdate, StereoConfigurationAst, StereoConfigurationUpdate, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, StereoLigandPosition,
    Stereogenicity, StereogenicityAst, SubPatternAnchor, TetrahedralStereoAst, Topicity,
    TopicityAst, TopicityRelationAst, TransactionError, UnpairedElectronsAst,
    UnpairedElectronsUpdate, ValueAst, ValuePredicate, ValueTerm,
};
pub(crate) use umol_ast::dsl::{
    parse_value, AromaticSystemDsl, AromaticSystemUpdateDsl, AtomDsl, AtomUpdateDsl, BondDsl,
    BondUpdateDsl, DativeBondDsl, DativeBondParticipants, DativeBondUpdateDsl, EditsDsl,
    MetadataError, MoleculeContext, MoleculeDefaults, MoleculeDsl, MoleculeMetadata,
    MulticenterBondDsl, MulticenterBondUpdateDsl, NoncovalentBondDsl, NoncovalentBondUpdateDsl,
    ParseError, ReactionDefaults, ReactionDsl, ReactionMetadata, ReactionSpanDsl,
    StereoAtomConstraintDsl, StereoAtomDsl, StereoAtomParticipants, StereoAtomUpdateDsl,
    StereoBondConstraintDsl, StereoBondDsl, StereoBondParticipants, StereoBondUpdateDsl,
    StereoLigandRef, ValueDsl,
};
pub(crate) use umol_chem::element::Element;
pub(crate) use umol_edn::{read_string, Edn, FromEdn, ToEdn};
use umol_graph_core::{Correspondence, EdgeId};
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
        3 => (0u32..=250).prop_map(|m| IsotopeMassAst::lit_set(vec![m])),
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
        3 => (0u32..=250).prop_map(IsotopeMassAst::Lit),
        2 => prop::collection::vec(0u32..=250, 1..=3).prop_map(IsotopeMassAst::lit_set),
        1 => id_strategy().prop_map(IsotopeMassAst::var),
        1 => (id_strategy(), prop::collection::vec(0u32..=250, 1..=3))
            .prop_map(|(id, v)| IsotopeMassAst::var_in(id, v)),
    ]
    .prop_map(|i| i.canonicalize().unwrap_or(IsotopeMassAst::Undetermined))
}

pub(crate) fn unpaired_electrons_strategy() -> impl Strategy<Value = UnpairedElectronsAst> {
    // The components are structurally independent; physical compatibility is
    // validated only when converting to `SpinState`.
    (value_basic(0..=6), value_basic(1..=7)).prop_map(|(count, multiplicity)| {
        UnpairedElectronsAst {
            count,
            multiplicity,
        }
    })
}

pub(crate) fn raw_unpaired_electrons_strategy() -> impl Strategy<Value = UnpairedElectronsAst> {
    (raw_value_ast_strategy(), raw_value_ast_strategy()).prop_map(|(count, multiplicity)| {
        UnpairedElectronsAst {
            count,
            multiplicity,
        }
    })
}

pub(crate) fn unpaired_electrons_update_strategy() -> impl Strategy<Value = UnpairedElectronsUpdate>
{
    (
        prop::option::of(value_basic(0..=6)),
        prop::option::of(value_basic(1..=7)),
    )
        .prop_map(|(count, multiplicity)| UnpairedElectronsUpdate {
            count,
            multiplicity,
        })
}

pub(crate) fn partial_unpaired_electrons_update_strategy(
) -> impl Strategy<Value = UnpairedElectronsUpdate> {
    prop_oneof![
        value_basic(0..=6).prop_map(|count| UnpairedElectronsUpdate {
            count: Some(count),
            multiplicity: None,
        }),
        value_basic(1..=7).prop_map(|multiplicity| UnpairedElectronsUpdate {
            count: None,
            multiplicity: Some(multiplicity),
        }),
    ]
}

/// `UnpairedElectronsAst` with at least one of `count` / `multiplicity` not
/// `Undetermined`. Used inside `MoleculeConstraint::UnpairedElectronCoupling` and similar
/// where a fully vacuous unpaired-electron state would elide on render.
pub(crate) fn non_vacuous_unpaired_electrons_strategy(
) -> impl Strategy<Value = UnpairedElectronsAst> {
    (value_basic(0..=6), value_basic(1..=7))
        .prop_map(|(u, m)| UnpairedElectronsAst {
            count: u,
            multiplicity: m,
        })
        .prop_filter("non-vacuous unpaired-electron state", |s| {
            !s.is_undetermined()
        })
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
        unpaired_electrons in unpaired_electrons_strategy(),
        constraints in atom_constraints_strategy(),
    ) -> AtomAst {
        AtomAst {
            element,
            isotope_mass: isotope,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired_electrons,
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
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in atom_update_constraints_strategy(),
    ) -> AtomUpdate {
        AtomUpdate {
            element,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired_electrons,
            constraints,
        }
    }
}

prop_compose! {
    pub(crate) fn bond_ast_strategy()
    (
        order in value_basic(1..=4),
        charge in value_basic(-1..=1),
        unpaired_electrons in unpaired_electrons_strategy(),
        constraints in bond_constraints_strategy(),
    ) -> BondAst {
        BondAst {
            order,
            charge,
            unpaired_electrons,
            constraints,
        }
    }
}

prop_compose! {
    pub(crate) fn bond_update_strategy()
    (
        order in prop::option::of(value_basic(1..=4)),
        charge in prop::option::of(value_basic(-1..=1)),
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in bond_update_constraints_strategy(),
    ) -> BondUpdate {
        BondUpdate {
            order,
            charge,
            unpaired_electrons,
            constraints,
        }
    }
}

/// `BondAst` shapes that render to bond keyword shorthands per spec §7.6:
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
            unpaired_electrons: UnpairedElectronsAst::default(),
            constraints,
        }
    })
}

pub(crate) fn aromatic_system_patch_ast_strategy() -> impl Strategy<Value = AromaticSystemAst> {
    (
        electron_counts_ast_strategy(),
        value_basic(-2..=2),
        unpaired_electrons_strategy(),
        optional_aromatic_electron_count(),
    )
        .prop_map(
            |(electrons, charge, unpaired_electrons, constraints)| AromaticSystemAst {
                electrons,
                charge,
                unpaired_electrons,
                constraints,
            },
        )
}

prop_compose! {
    pub(crate) fn aromatic_system_update_strategy()
    (
        electrons in prop::option::of(electron_counts_ast_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in aromatic_system_update_constraints_strategy(),
    ) -> AromaticSystemUpdate {
        AromaticSystemUpdate { electrons, charge, unpaired_electrons, constraints }
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
            unpaired_electrons: UnpairedElectronsAst::default(),
            constraints,
        })
}

pub(crate) fn multicenter_bond_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (value_basic(-2..=2), optional_multicenter_electron_count()).prop_map(
        |(charge, constraints)| MulticenterBondAst {
            electrons: ElectronCountsAst::Undetermined,
            charge,
            unpaired_electrons: UnpairedElectronsAst::default(),
            constraints,
        },
    )
}

pub(crate) fn multicenter_bond_patch_ast_strategy() -> impl Strategy<Value = MulticenterBondAst> {
    (
        electron_counts_ast_strategy(),
        value_basic(-2..=2),
        unpaired_electrons_strategy(),
        optional_multicenter_electron_count(),
    )
        .prop_map(
            |(electrons, charge, unpaired_electrons, constraints)| MulticenterBondAst {
                electrons,
                charge,
                unpaired_electrons,
                constraints,
            },
        )
}

prop_compose! {
    pub(crate) fn multicenter_bond_update_strategy()
    (
        electrons in prop::option::of(electron_counts_ast_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in multicenter_bond_update_constraints_strategy(),
    ) -> MulticenterBondUpdate {
        MulticenterBondUpdate { electrons, charge, unpaired_electrons, constraints }
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
            unpaired_electrons: UnpairedElectronsAst::default(),
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
/// (§5.8) and excluded.
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

pub(crate) fn molecule_entries_strategy() -> impl Strategy<Value = MoleculeEntries> {
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
                MoleculeEntries {
                    atoms,
                    bonds,
                    dative: dative_triples,
                    aromatic: aromatic_entries,
                    multicenter: multicenter_entries,
                    noncovalent: noncovalent_triples,
                    stereo_atoms,
                    stereo_bonds,
                    constraints: Constraints::new(),
                }
            },
        )
}

pub(crate) fn molecule_ast_strategy() -> impl Strategy<Value = MoleculeAst> {
    molecule_entries_strategy().prop_map(MoleculeAst::from_entries)
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
    fn from_entries(entries: &MoleculeEntries) -> Self {
        Self {
            atom: entries.atoms.len(),
            bond: entries.bonds.len(),
            dative: entries.dative.len(),
            aromatic: entries.aromatic.len(),
            multicenter: entries.multicenter.len(),
            noncovalent: entries.noncovalent.len(),
            stereo_atom: entries.stereo_atoms.len(),
            stereo_bond: entries.stereo_bonds.len(),
        }
    }

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
        // undetermined unpaired-electron state) elide on render and would break round-trip;
        // restrict the value and unpaired-electron strategies accordingly.
        let molecule_charge_sum = (
            optional_atoms.clone(),
            constraint_inner_value_strategy(-3..=3),
        )
            .prop_map(|(atoms, sum)| {
                Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms, sum })
            })
            .boxed();
        let molecule_unpaired_electron_coupling =
            (optional_atoms, non_vacuous_unpaired_electrons_strategy())
                .prop_map(|(atoms, unpaired_electrons)| {
                    Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                        atoms,
                        unpaired_electrons,
                    })
                })
                .boxed();
        choices.push(molecule_connected);
        choices.push(molecule_charge_sum);
        choices.push(molecule_unpaired_electron_coupling);
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
    molecule_entries_with_constraints_strategy().prop_map(MoleculeAst::from_entries)
}

pub(crate) fn molecule_entries_with_constraints_strategy() -> impl Strategy<Value = MoleculeEntries>
{
    molecule_entries_strategy().prop_flat_map(|entries| {
        let counts = ConstraintCounts::from_entries(&entries);
        let max_constraints = 4usize;
        (
            Just(entries),
            prop::collection::vec(constraint_strategy(counts), 0..=max_constraints),
        )
            .prop_map(|(mut entries, constraints)| {
                let mut cs = Constraints::new();
                for c in constraints {
                    cs.push(c);
                }
                entries.constraints = cs;
                entries
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
        molecule_entity_incidence_is_unique,
    )
}

fn molecule_entity_incidence_is_unique(ast: &MoleculeAst) -> bool {
    all_unique(ast.bonds().iter().map(|bond| sorted_pair(bond.atom_ids())))
        && all_unique(
            ast.dative_bonds()
                .iter()
                .map(|dative| (dative.acceptor_id(), sorted(dative.donor_ids().collect()))),
        )
        && all_unique(
            ast.aromatic_systems()
                .iter()
                .map(|aromatic| sorted(aromatic.atom_ids().collect())),
        )
        && all_unique(
            ast.multicenter_bonds()
                .iter()
                .map(|multicenter| sorted(multicenter.atom_ids().collect())),
        )
        && all_unique(
            ast.noncovalent_bonds()
                .iter()
                .map(|noncovalent| sorted_pair(noncovalent.atom_ids())),
        )
        && all_unique(
            ast.stereo_atoms()
                .iter()
                .map(|stereo| (stereo.site_id(), sorted(stereo.ligand_frame()))),
        )
        && all_unique(
            ast.stereo_bonds()
                .iter()
                .map(|stereo| (stereo.site_id(), sorted(stereo.ligand_frame()))),
        )
}

fn all_unique<T>(values: impl IntoIterator<Item = T>) -> bool
where
    T: Eq + Hash,
{
    let mut seen = HashSet::new();
    values.into_iter().all(|value| seen.insert(value))
}

fn sorted<T: Ord>(mut values: Vec<T>) -> Vec<T> {
    values.sort_unstable();
    values
}

fn sorted_pair<T: Ord>([first, second]: [T; 2]) -> [T; 2] {
    if first <= second {
        [first, second]
    } else {
        [second, first]
    }
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
    let stereo_atom_flags = prop::collection::vec(id_flag(), counts.stereo_atom);
    let stereo_bond_flags = prop::collection::vec(id_flag(), counts.stereo_bond);
    (
        atom_flags,
        bond_flags,
        dative_flags,
        aromatic_flags,
        multicenter_flags,
        noncovalent_flags,
        stereo_atom_flags,
        stereo_bond_flags,
    )
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
                for (i, stereo_atom) in stereo_atoms.iter().enumerate() {
                    if stereo_atom.is_some() {
                        meta.set_keyword(
                            Entity::StereoAtom(StereoAtomId(i as u32)),
                            format!("stereo_atom{i}"),
                        )
                        .unwrap();
                    }
                }
                for (i, stereo_bond) in stereo_bonds.iter().enumerate() {
                    if stereo_bond.is_some() {
                        meta.set_keyword(
                            Entity::StereoBond(StereoBondId(i as u32)),
                            format!("stereo_bond{i}"),
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

pub(crate) fn invalid_metadata_for(
    counts: ConstraintCounts,
) -> BoxedStrategy<(MoleculeMetadata, Entity)> {
    (0u8..8)
        .prop_map(move |kind| {
            let kind = EntityKind::try_from(kind).unwrap();
            let count = match kind {
                EntityKind::Atom => counts.atom,
                EntityKind::Bond => counts.bond,
                EntityKind::DativeBond => counts.dative,
                EntityKind::AromaticSystem => counts.aromatic,
                EntityKind::MulticenterBond => counts.multicenter,
                EntityKind::NoncovalentBond => counts.noncovalent,
                EntityKind::StereoAtom => counts.stereo_atom,
                EntityKind::StereoBond => counts.stereo_bond,
            };
            let entity = kind.with_id(count as u32);
            let mut metadata = MoleculeMetadata::new();
            metadata.set_keyword(entity, "invalid".to_string()).unwrap();
            (metadata, entity)
        })
        .boxed()
}

pub(crate) fn molecule_dsl_strategy() -> impl Strategy<Value = MoleculeDsl> {
    molecule_ast_with_constraints_strategy().prop_flat_map(|ast| {
        let counts = ConstraintCounts::from_ast(&ast);
        metadata_for(counts).prop_map(move |metadata| {
            MoleculeDsl::new(ast.clone(), metadata).expect("generated metadata is coherent")
        })
    })
}

pub(crate) fn invalid_molecule_dsl_parts_strategy(
) -> impl Strategy<Value = (MoleculeAst, MoleculeMetadata, Entity)> {
    molecule_ast_with_constraints_strategy().prop_flat_map(|ast| {
        let counts = ConstraintCounts::from_ast(&ast);
        invalid_metadata_for(counts)
            .prop_map(move |(metadata, entity)| (ast.clone(), metadata, entity))
    })
}

pub(crate) fn molecule_metadata_with_atom_subset_strategy(
) -> impl Strategy<Value = (MoleculeAst, MoleculeMetadata, Vec<AtomId>)> {
    molecule_ast_with_atom_subset_strategy().prop_flat_map(|(ast, atoms)| {
        metadata_for(ConstraintCounts::from_ast(&ast))
            .prop_map(move |metadata| (ast.clone(), metadata, atoms.clone()))
    })
}

fn added_entities(reaction: &ReactionAst) -> Vec<Entity> {
    reaction
        .deltas
        .iter()
        .filter_map(|delta| match delta {
            Delta::Atom(AtomDelta::Add { id, .. }) => Some(Entity::Atom(*id)),
            Delta::Bond(BondDelta::Add { id, .. }) => Some(Entity::Bond(*id)),
            Delta::DativeBond(DativeBondDelta::Add { id, .. }) => Some(Entity::DativeBond(*id)),
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, .. }) => {
                Some(Entity::AromaticSystem(*id))
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, .. }) => {
                Some(Entity::MulticenterBond(*id))
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, .. }) => {
                Some(Entity::NoncovalentBond(*id))
            }
            Delta::StereoAtom(StereoAtomDelta::Add { id, .. }) => Some(Entity::StereoAtom(*id)),
            Delta::StereoBond(StereoBondDelta::Add { id, .. }) => Some(Entity::StereoBond(*id)),
            _ => None,
        })
        .collect()
}

fn delta_keyword(entity: Entity) -> String {
    let kind = match entity {
        Entity::Atom(_) => "atom",
        Entity::Bond(_) => "bond",
        Entity::DativeBond(_) => "dative",
        Entity::AromaticSystem(_) => "aromatic",
        Entity::MulticenterBond(_) => "multicenter",
        Entity::NoncovalentBond(_) => "noncovalent",
        Entity::StereoAtom(_) => "stereo_atom",
        Entity::StereoBond(_) => "stereo_bond",
    };
    format!("delta_{kind}_{}", entity.id_index())
}

pub(crate) fn reaction_dsl_strategy() -> impl Strategy<Value = ReactionDsl> {
    comprehensive_reaction_strategy().prop_flat_map(|reaction| {
        metadata_for(ConstraintCounts::from_ast(&reaction.lhs)).prop_map(move |lhs| {
            let mut metadata = ReactionMetadata::from(lhs);
            for entity in added_entities(&reaction) {
                metadata
                    .set_delta_keyword(entity, delta_keyword(entity))
                    .expect("generated delta keywords are disjoint");
            }
            metadata
                .add_atom_alias("reaction_alias", AtomAst::from_element(Element::F))
                .expect("generated reaction alias is disjoint and bijective");
            ReactionDsl::new(reaction.clone(), metadata)
                .expect("generated reaction metadata is coherent")
        })
    })
}

pub(crate) fn invalid_reaction_dsl_parts_strategy(
) -> impl Strategy<Value = (ReactionAst, ReactionMetadata, MetadataError)> {
    prop_oneof![
        comprehensive_reaction_strategy().prop_flat_map(|reaction| {
            invalid_metadata_for(ConstraintCounts::from_ast(&reaction.lhs)).prop_map(
                move |(lhs, entity)| {
                    (
                        reaction.clone(),
                        ReactionMetadata::from(lhs),
                        MetadataError::EntityOutOfRange(entity),
                    )
                },
            )
        }),
        (comprehensive_reaction_strategy(), 0u8..8).prop_map(|(reaction, kind)| {
            let entity = EntityKind::try_from(kind).unwrap().with_id(u32::MAX);
            let mut metadata = ReactionMetadata::default();
            metadata
                .set_delta_keyword(entity, "invalid")
                .expect("fresh invalid keyword is unique");
            (reaction, metadata, MetadataError::EntityNotAdded(entity))
        }),
    ]
}

pub(crate) fn reaction_span_dsl_strategy() -> impl Strategy<Value = ReactionSpanDsl> {
    comprehensive_reaction_strategy()
        .prop_filter_map("reaction must have a materializable span", |reaction| {
            reaction.to_reaction_span().ok()
        })
        .prop_flat_map(|span| {
            metadata_for(ConstraintCounts {
                atom: span.atoms().len(),
                bond: span.bonds().len(),
                dative: 0,
                aromatic: 0,
                multicenter: 0,
                noncovalent: 0,
                stereo_atom: 0,
                stereo_bond: 0,
            })
            .prop_map(move |metadata| {
                ReactionSpanDsl::new(span.clone(), metadata)
                    .expect("generated reaction-span metadata is coherent")
            })
        })
}

pub(crate) fn invalid_reaction_span_dsl_parts_strategy(
) -> impl Strategy<Value = (ReactionSpanAst, MoleculeMetadata, Entity)> {
    comprehensive_reaction_strategy()
        .prop_filter_map("reaction must have a materializable span", |reaction| {
            reaction.to_reaction_span().ok()
        })
        .prop_map(|span| {
            let entity = Entity::Atom(AtomId(span.atoms().len() as u32));
            let mut metadata = MoleculeMetadata::new();
            metadata
                .set_keyword(entity, "invalid")
                .expect("fresh invalid keyword is unique");
            (span, metadata, entity)
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
        .map(|id| AddBond {
            endpoints: [AtomHandle::New(id), AtomHandle::New(id + 1)],
            ast: BondAst::from_order((id % 3 + 1) as u8),
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
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

pub(crate) fn transaction_add_path_edits(count: usize) -> Edits {
    let mut edits = Edits::new();
    edits.add_atoms(transaction_atoms(count));
    edits.add_bonds(transaction_path_bonds(count));
    edits
}

const INITIAL_HANDLE_ELEMENTS: [Element; 4] = [Element::H, Element::C, Element::N, Element::O];
const CREATED_HANDLE_ELEMENTS: [Element; 4] = [Element::F, Element::P, Element::S, Element::Cl];
const SENTINEL_HANDLE_ELEMENT: Element = Element::Br;

/// A compacting atom transaction with independently labeled initial and created entities.
///
/// Removal masks and the target origin/index are generated independently, then the target's mask is
/// fixed to the requested liveness. Expected states filter these labels directly rather than using
/// transaction compaction machinery.
#[derive(Clone, Debug)]
pub(crate) struct StableAtomHandleTrace {
    pub(crate) initial_count: usize,
    pub(crate) created_count: usize,
    pub(crate) remove_initial: Vec<bool>,
    pub(crate) remove_created: Vec<bool>,
    pub(crate) target_created: bool,
    pub(crate) target_index: usize,
}

impl StableAtomHandleTrace {
    pub(crate) fn base(&self) -> MoleculeAst {
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: INITIAL_HANDLE_ELEMENTS[..self.initial_count]
                .iter()
                .copied()
                .map(AtomAst::from_element)
                .collect(),
            ..Default::default()
        })
    }

    pub(crate) fn edits(&self) -> Edits {
        let mut edits = Edits::new();
        let created = edits.add_atoms(
            CREATED_HANDLE_ELEMENTS[..self.created_count]
                .iter()
                .copied()
                .map(AtomAst::from_element),
        );
        edits.remove_topology(
            self.remove_initial
                .iter()
                .enumerate()
                .filter(|&(_, remove)| *remove)
                .map(|(index, _)| AtomHandle::Id(AtomId(index as u32)))
                .chain(
                    self.remove_created
                        .iter()
                        .enumerate()
                        .filter(|&(_, remove)| *remove)
                        .map(|(index, _)| created[index].clone()),
                )
                .collect(),
            Vec::new(),
        );
        let sentinel = edits.add_atom(AtomAst::from_element(SENTINEL_HANDLE_ELEMENT));
        edits.push(Edit::ModifyAtomField {
            id: if self.target_created {
                created[self.target_index].clone()
            } else {
                AtomHandle::Id(AtomId(self.target_index as u32))
            },
            change: AtomFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(7),
            },
        });
        edits.push(Edit::ModifyAtomField {
            id: sentinel,
            change: AtomFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(9),
            },
        });
        edits
    }

    pub(crate) fn expected(&self) -> MoleculeAst {
        let initial = INITIAL_HANDLE_ELEMENTS[..self.initial_count]
            .iter()
            .copied()
            .enumerate()
            .filter(|(index, _)| !self.remove_initial[*index])
            .map(|(index, element)| {
                if !self.target_created && index == self.target_index {
                    AtomAst::from_element(element).with_charge(7_i64)
                } else {
                    AtomAst::from_element(element)
                }
            });
        let created = CREATED_HANDLE_ELEMENTS[..self.created_count]
            .iter()
            .copied()
            .enumerate()
            .filter(|(index, _)| !self.remove_created[*index])
            .map(|(index, element)| {
                if self.target_created && index == self.target_index {
                    AtomAst::from_element(element).with_charge(7_i64)
                } else {
                    AtomAst::from_element(element)
                }
            });
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: initial
                .chain(created)
                .chain([AtomAst::from_element(SENTINEL_HANDLE_ELEMENT).with_charge(9_i64)])
                .collect(),
            ..Default::default()
        })
    }

    pub(crate) fn expected_removed_error(&self) -> TransactionError {
        TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: self.target_index,
        }
    }
}

pub(crate) fn stable_atom_handle_trace_strategy(
    target_removed: bool,
) -> impl Strategy<Value = StableAtomHandleTrace> {
    (1usize..=4, 1usize..=4, any::<bool>())
        .prop_flat_map(|(initial_count, created_count, target_created)| {
            let target_count = if target_created {
                created_count
            } else {
                initial_count
            };
            (
                Just((initial_count, created_count, target_created)),
                prop::collection::vec(any::<bool>(), initial_count),
                prop::collection::vec(any::<bool>(), created_count),
                0usize..target_count,
            )
        })
        .prop_map(
            move |(
                (initial_count, created_count, target_created),
                mut remove_initial,
                mut remove_created,
                target_index,
            )| {
                if target_created {
                    remove_created[target_index] = target_removed;
                } else {
                    remove_initial[target_index] = target_removed;
                }
                StableAtomHandleTrace {
                    initial_count,
                    created_count,
                    remove_initial,
                    remove_created,
                    target_created,
                    target_index,
                }
            },
        )
        .prop_filter("trace must compact at least one atom", |trace| {
            trace
                .remove_initial
                .iter()
                .chain(&trace.remove_created)
                .any(|remove| *remove)
        })
}

pub(crate) fn transaction_entity_kind_order_strategy() -> impl Strategy<Value = Vec<EntityKind>> {
    (
        1usize..=8,
        Just(vec![
            EntityKind::Atom,
            EntityKind::Bond,
            EntityKind::DativeBond,
            EntityKind::AromaticSystem,
            EntityKind::MulticenterBond,
            EntityKind::NoncovalentBond,
            EntityKind::StereoAtom,
            EntityKind::StereoBond,
        ])
        .prop_shuffle(),
    )
        .prop_map(|(count, mut kinds)| {
            kinds.truncate(count);
            kinds
        })
}

/// One generated batched edit with exactly one invalid handle at an arbitrary position.
#[derive(Clone, Debug)]
pub(crate) struct InvalidTransactionBatch {
    pub(crate) kind: EntityKind,
    pub(crate) count: usize,
    pub(crate) invalid_position: usize,
}

impl InvalidTransactionBatch {
    pub(crate) fn base(&self) -> MoleculeAst {
        let atoms = (0..self.count * 2)
            .map(|index| AtomAst::from_element(INITIAL_HANDLE_ELEMENTS[index % 4]))
            .collect();
        let bonds = (0..self.count)
            .map(|index| {
                (
                    AtomId((index * 2) as u32),
                    AtomId((index * 2 + 1) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        let dative = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32)],
                    AtomId((index * 2 + 1) as u32),
                    DativeBondAst::from_order(1),
                )
            })
            .collect();
        let aromatic = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32), AtomId((index * 2 + 1) as u32)],
                    AromaticSystemAst::default(),
                )
            })
            .collect();
        let multicenter = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32), AtomId((index * 2 + 1) as u32)],
                    MulticenterBondAst::default(),
                )
            })
            .collect();
        let noncovalent = (0..self.count)
            .map(|index| {
                (
                    AtomId((index * 2) as u32),
                    AtomId((index * 2 + 1) as u32),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )
            })
            .collect();
        let stereo_atoms = (0..self.count)
            .map(|index| {
                (
                    AtomId((index * 2) as u32),
                    vec![StereoLigand::new(
                        AtomId((index * 2 + 1) as u32),
                        StereoLigandKind::Atom,
                    )],
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )
            })
            .collect();
        let stereo_bonds = (0..self.count)
            .map(|index| {
                (
                    BondId(index as u32),
                    vec![
                        StereoLigand::new(AtomId((index * 2) as u32), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId((index * 2 + 1) as u32), StereoLigandKind::Atom),
                    ],
                    StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )
            })
            .collect();
        MoleculeAst::from_entries(MoleculeEntries {
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
    }

    pub(crate) fn edits(&self) -> Edits {
        let invalid = self.count as u32;
        let edit = match self.kind {
            EntityKind::Bond => Edit::AddBonds {
                bonds: (0..self.count)
                    .map(|position| AddBond {
                        endpoints: [
                            AtomHandle::Id(AtomId(0)),
                            AtomHandle::Id(AtomId(if position == self.invalid_position {
                                invalid * 2
                            } else {
                                1
                            })),
                        ],
                        ast: BondAst::from_order(1),
                    })
                    .collect(),
            },
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            DativeBondHandle::Id(DativeBondId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId((position * 2) as u32)),
                                AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                            ],
                            DativeBondAst::from_order(1),
                        )
                    })
                    .collect(),
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId((position * 2) as u32)),
                                AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                            ],
                            AromaticSystemAst::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            MulticenterBondHandle::Id(MulticenterBondId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId((position * 2) as u32)),
                                AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                            ],
                            MulticenterBondAst::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            NoncovalentBondHandle::Id(NoncovalentBondId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            [
                                AtomHandle::Id(AtomId((position * 2) as u32)),
                                AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                            ],
                            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            StereoAtomHandle::Id(StereoAtomId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            AtomHandle::Id(AtomId((position * 2) as u32)),
                            vec![(
                                AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                                StereoLigandKind::Atom,
                            )],
                            StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: (0..self.count)
                    .map(|position| {
                        (
                            StereoBondHandle::Id(StereoBondId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            BondHandle::Id(BondId(position as u32)),
                            vec![
                                (
                                    AtomHandle::Id(AtomId((position * 2) as u32)),
                                    StereoLigandKind::Atom,
                                ),
                                (
                                    AtomHandle::Id(AtomId((position * 2 + 1) as u32)),
                                    StereoLigandKind::Atom,
                                ),
                            ],
                            StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::Atom => unreachable!(),
        };
        Edits::from_iter([edit])
    }

    pub(crate) fn expected_error(&self) -> TransactionError {
        if self.kind == EntityKind::Bond {
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: self.count * 2,
                count: self.count * 2,
            }
        } else {
            TransactionError::HandleOutOfRange {
                kind: self.kind,
                index: self.count,
                count: self.count,
            }
        }
    }
}

pub(crate) fn invalid_transaction_batch_strategy() -> impl Strategy<Value = InvalidTransactionBatch>
{
    (
        prop::sample::select(vec![
            EntityKind::Bond,
            EntityKind::DativeBond,
            EntityKind::AromaticSystem,
            EntityKind::MulticenterBond,
            EntityKind::NoncovalentBond,
            EntityKind::StereoAtom,
            EntityKind::StereoBond,
        ]),
        1usize..=5,
    )
        .prop_flat_map(|(kind, count)| (Just((kind, count)), 0usize..count))
        .prop_map(
            |((kind, count), invalid_position)| InvalidTransactionBatch {
                kind,
                count,
                invalid_position,
            },
        )
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

    pub(crate) fn edits(&self) -> Edits {
        match self {
            Self::AddPath { count } => transaction_add_path_edits(*count),
            Self::RemoveAtom { count, id } => {
                let mut edits = Edits::new();
                edits.remove_atom(AtomHandle::Id(AtomId((id % count) as u32)));
                edits
            }
            Self::RemoveBond { count, id } => {
                let mut edits = Edits::new();
                edits.remove_bond(BondHandle::Id(BondId((id % (count - 1)) as u32)));
                edits
            }
            Self::SetAtomCharge { count, id, charge } => {
                Edits::from_iter([Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId((id % count) as u32)),
                    change: AtomFieldChange::Charge {
                        old: ValueAst::default(),
                        new: ValueAst::Lit(*charge),
                    },
                }])
            }
            Self::SetBondOrder { count, id, order } => {
                let bond_id = id % (count - 1);
                Edits::from_iter([Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(bond_id as u32)),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit((bond_id % 3 + 1) as i64),
                        new: ValueAst::Lit(*order as i64),
                    },
                }])
            }
            Self::AddAtomConstraint { count, id, size } => {
                Edits::from_iter([Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId((id % count) as u32)),
                    old: None,
                    new: Some(AtomConstraintAst::ring_membership(
                        RingScope::Size(*size as u8),
                        1,
                    )),
                }])
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
                let mut edits = Edits::new();
                edits.add_dative_bond(
                    vec![
                        AtomHandle::Id(AtomId(donor as u32)),
                        AtomHandle::Id(AtomId(acceptor as u32)),
                    ],
                    DativeBondAst::from_order(1),
                );
                edits
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

fn transaction_all_entities_molecule() -> MoleculeAst {
    let ligands = (0..4)
        .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: (0..4).map(|_| AtomAst::from_element(Element::C)).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondAst::default(),
        )],
        noncovalent: vec![(
            AtomId(0),
            AtomId(3),
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            ligands.clone(),
            StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )],
        stereo_bonds: vec![(
            BondId(0),
            ligands,
            StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    })
}

fn transaction_field_cases() -> Vec<(MoleculeAst, Edits)> {
    let base = transaction_all_entities_molecule();
    let value = |change| (base.clone(), Edits::from_iter([change]));
    vec![
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Element {
                old: ElementAst::Lit(Element::C),
                new: ElementAst::Lit(Element::N),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::IsotopeMass {
                old: IsotopeMassAst::default(),
                new: IsotopeMassAst::Lit(13),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(1),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::ImplicitHydrogens {
                old: ValueAst::default(),
                new: ValueAst::Lit(3),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::LonePairs {
                old: ValueAst::default(),
                new: ValueAst::Lit(1),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst::default(),
                new: UnpairedElectronsAst::from((2_u8, 1_u8)),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Order {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(-1),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst::default(),
                new: UnpairedElectronsAst::from((2_u8, 3_u8)),
            },
        }),
        value(Edit::ModifyDativeBondField {
            id: DativeBondHandle::Id(DativeBondId(0)),
            change: DativeBondFieldChange::Order {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsAst::default(),
                new: ElectronCountsAst::Lit(vec![1, 1, 1]),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(1),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst::default(),
                new: UnpairedElectronsAst::from((1_u8, 2_u8)),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::Electrons {
                old: ElectronCountsAst::default(),
                new: ElectronCountsAst::Lit(vec![1, 1, 1]),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::Charge {
                old: ValueAst::default(),
                new: ValueAst::Lit(-1),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst::default(),
                new: UnpairedElectronsAst::from((1_u8, 2_u8)),
            },
        }),
        value(Edit::ModifyNoncovalentBondField {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
            },
        }),
        value(Edit::ModifyStereoAtomField {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
        }),
        value(Edit::ModifyStereoBondField {
            id: StereoBondHandle::Id(StereoBondId(0)),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                new: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        }),
    ]
}

fn transaction_constraint_cases() -> Vec<(MoleculeAst, Edits)> {
    let base = transaction_all_entities_molecule();
    vec![
        Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintAst::degree(3)),
        },
        Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: None,
            new: Some(BondConstraintAst::aromatic(true)),
        },
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::Id(DativeBondId(0)),
            old: None,
            new: Some(DativeBondConstraintAst::aromatic(true)),
        },
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            old: None,
            new: Some(AromaticSystemConstraintAst::electron_count(6)),
        },
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            old: None,
            new: Some(MulticenterBondConstraintAst::electron_count(2)),
        },
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            old: None,
            new: Some(NoncovalentBondConstraintAst::intramolecular(true)),
        },
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoAtomConstraintAst::Stereogenicity(
                StereogenicityAst::Lit(Stereogenicity::Stereogenic),
            )),
        },
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::Id(StereoBondId(0)),
            kind: Some(StereoKind::CisTrans),
            old: None,
            new: Some(StereoBondConstraintAst::Stereogenicity(
                StereogenicityAst::Lit(Stereogenicity::Stereogenic),
            )),
        },
    ]
    .into_iter()
    .map(|edit| (base.clone(), Edits::from_iter([edit])))
    .collect()
}

fn transaction_removal_cases() -> Vec<(MoleculeAst, Edits)> {
    let base = transaction_all_entities_molecule();
    let atom_handles = |ids: &[u32]| {
        ids.iter()
            .map(|&id| AtomHandle::Id(AtomId(id)))
            .collect::<Vec<_>>()
    };
    vec![
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(3))],
            bonds: Vec::new(),
        },
        Edit::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![BondHandle::Id(BondId(1))],
        },
        Edit::RemoveDativeBonds {
            removes: vec![(
                DativeBondHandle::Id(DativeBondId(0)),
                atom_handles(&[0, 1]),
                DativeBondAst::from_order(1),
            )],
        },
        Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                atom_handles(&[0, 1, 2]),
                AromaticSystemAst::default(),
            )],
        },
        Edit::RemoveMulticenterBonds {
            removes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                atom_handles(&[0, 1, 2]),
                MulticenterBondAst::default(),
            )],
        },
        Edit::RemoveNoncovalentBonds {
            removes: vec![(
                NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(3))],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        },
        Edit::RemoveStereoAtoms {
            removes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(0)),
                atom_handles(&[0, 1, 2, 3])
                    .into_iter()
                    .map(|atom| (atom, StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
        },
        Edit::RemoveStereoBonds {
            removes: vec![(
                StereoBondHandle::Id(StereoBondId(0)),
                BondHandle::Id(BondId(0)),
                atom_handles(&[0, 1, 2, 3])
                    .into_iter()
                    .map(|atom| (atom, StereoLigandKind::Atom))
                    .collect(),
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
        },
    ]
    .into_iter()
    .map(|edit| (base.clone(), Edits::from_iter([edit])))
    .collect()
}

fn transaction_creation_case(include_created_constraint: bool) -> (MoleculeAst, Edits) {
    let base = transaction_all_entities_molecule();
    let mut edits = Edits::new();
    let atom = edits.add_atom(AtomAst::from_element(Element::N));
    let bond = edits.add_bond(
        AtomHandle::Id(AtomId(1)),
        AtomHandle::Id(AtomId(2)),
        BondAst::from_order(2),
    );
    let dative = edits.add_dative_bond(
        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        DativeBondAst::from_order(1),
    );
    let aromatic = edits.add_aromatic_system(
        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        AromaticSystemAst::default(),
    );
    let multicenter = edits.add_multicenter_bond(
        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        MulticenterBondAst::default(),
    );
    let noncovalent = edits.add_noncovalent_bond(
        [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic),
    );
    let ligands = (0..4)
        .map(|id| (AtomHandle::Id(AtomId(id)), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let stereo_atom = edits.add_stereo_atom(
        AtomHandle::Id(AtomId(1)),
        ligands.clone(),
        StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
    );
    let stereo_bond = edits.add_stereo_bond(
        BondHandle::Id(BondId(1)),
        ligands,
        StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
    );
    let source = Constraint::And(vec![
        Constraint::Atom(AtomId(7), AtomConstraintAst::degree(3)),
        Constraint::Bond(BondId(7), BondConstraintAst::aromatic(true)),
        Constraint::DativeBond(DativeBondId(7), DativeBondConstraintAst::aromatic(true)),
        Constraint::AromaticSystem(
            AromaticSystemId(7),
            AromaticSystemConstraintAst::electron_count(6),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(7),
            MulticenterBondConstraintAst::electron_count(2),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(7),
            NoncovalentBondConstraintAst::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(7),
            StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(7),
            StereoKind::CisTrans,
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
    ]);
    let mappings = HashMap::from([
        (Entity::Atom(AtomId(7)), EntityHandle::Atom(atom)),
        (Entity::Bond(BondId(7)), EntityHandle::Bond(bond)),
        (
            Entity::DativeBond(DativeBondId(7)),
            EntityHandle::DativeBond(dative),
        ),
        (
            Entity::AromaticSystem(AromaticSystemId(7)),
            EntityHandle::AromaticSystem(aromatic),
        ),
        (
            Entity::MulticenterBond(MulticenterBondId(7)),
            EntityHandle::MulticenterBond(multicenter),
        ),
        (
            Entity::NoncovalentBond(NoncovalentBondId(7)),
            EntityHandle::NoncovalentBond(noncovalent),
        ),
        (
            Entity::StereoAtom(StereoAtomId(7)),
            EntityHandle::StereoAtom(stereo_atom),
        ),
        (
            Entity::StereoBond(StereoBondId(7)),
            EntityHandle::StereoBond(stereo_bond),
        ),
    ]);
    if include_created_constraint {
        edits.add_molecule_constraint(
            ConstraintEdit::new(source, |entity| mappings.get(&entity).cloned()).unwrap(),
        );
    }
    (base, edits)
}

fn complete_transaction_cases(include_created_constraint: bool) -> Vec<(MoleculeAst, Edits)> {
    let mut cases = transaction_field_cases();
    cases.extend(transaction_constraint_cases());
    cases.extend(transaction_removal_cases());
    cases.push(transaction_creation_case(include_created_constraint));
    let constraint = Constraint::Atom(AtomId(0), AtomConstraintAst::degree(3));
    cases.push((
        transaction_all_entities_molecule(),
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: constraint.clone().into(),
        }]),
    ));
    cases.push((
        transaction_compaction_molecule(Constraints::from_iter([
            constraint.clone(),
            constraint.clone(),
        ])),
        Edits::from_iter([Edit::RemoveMoleculeConstraint {
            constraint: constraint.into(),
        }]),
    ));
    cases
}

pub(crate) fn complete_transaction_strategy() -> impl Strategy<Value = (MoleculeAst, Edits)> {
    prop::sample::select(complete_transaction_cases(true))
}

fn transaction_compaction_molecule(constraints: Constraints) -> MoleculeAst {
    let atoms = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
    let bonds = (0..3)
        .map(|index| {
            (
                AtomId(index * 2),
                AtomId(index * 2 + 1),
                BondAst::from_order(1),
            )
        })
        .collect();
    let pairs = [[0_u32, 1_u32], [2, 3], [4, 5]];
    let dative = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a)], AtomId(*b), DativeBondAst::from_order(1)))
        .collect();
    let aromatic = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a), AtomId(*b)], AromaticSystemAst::default()))
        .collect();
    let multicenter = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a), AtomId(*b)], MulticenterBondAst::default()))
        .collect();
    let noncovalent = pairs
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                AtomId(*b),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )
        })
        .collect();
    let stereo_atoms = pairs
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                vec![StereoLigand::new(AtomId(*b), StereoLigandKind::Atom)],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )
        })
        .collect();
    let stereo_bonds = pairs
        .iter()
        .enumerate()
        .map(|(index, [a, b])| {
            (
                BondId(index as u32),
                vec![
                    StereoLigand::new(AtomId(*a), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(*b), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )
        })
        .collect();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        dative,
        aromatic,
        multicenter,
        noncovalent,
        stereo_atoms,
        stereo_bonds,
        constraints,
    })
}

fn transaction_constraint(kind: EntityKind, id: u32, value: i64) -> Constraint {
    match kind {
        EntityKind::Atom => Constraint::Atom(AtomId(id), AtomConstraintAst::degree(value)),
        EntityKind::Bond => Constraint::Bond(
            BondId(id),
            BondConstraintAst::ring_membership(RingScope::All, value),
        ),
        EntityKind::DativeBond => Constraint::DativeBond(
            DativeBondId(id),
            DativeBondConstraintAst::aromatic(value % 2 == 0),
        ),
        EntityKind::AromaticSystem => Constraint::AromaticSystem(
            AromaticSystemId(id),
            AromaticSystemConstraintAst::electron_count(value),
        ),
        EntityKind::MulticenterBond => Constraint::MulticenterBond(
            MulticenterBondId(id),
            MulticenterBondConstraintAst::electron_count(value),
        ),
        EntityKind::NoncovalentBond => Constraint::NoncovalentBond(
            NoncovalentBondId(id),
            NoncovalentBondConstraintAst::intramolecular(value % 2 == 0),
        ),
        EntityKind::StereoAtom => Constraint::StereoAtom(
            StereoAtomId(id),
            StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(match value % 3 {
                0 => Stereogenicity::Symmetric,
                1 => Stereogenicity::Prochiral,
                _ => Stereogenicity::Stereogenic,
            })),
        ),
        EntityKind::StereoBond => Constraint::StereoBond(
            StereoBondId(id),
            StereoKind::CisTrans,
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(match value % 3 {
                0 => Stereogenicity::Symmetric,
                1 => Stereogenicity::Prochiral,
                _ => Stereogenicity::Stereogenic,
            })),
        ),
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ConstraintCompactionCase {
    kind: EntityKind,
    base: MoleculeAst,
    expected: Vec<Constraint>,
}

impl ConstraintCompactionCase {
    pub(crate) fn base(&self) -> MoleculeAst {
        self.base.clone()
    }

    pub(crate) fn edits(&self) -> Edits {
        let edit = match self.kind {
            EntityKind::Atom => Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            EntityKind::Bond => Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0))],
            },
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    DativeBondAst::from_order(1),
                )],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemAst::default(),
                )],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondAst::default(),
                )],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    ],
                    StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )],
            },
        };
        Edits::from_iter([edit])
    }

    pub(crate) fn expected(&self) -> &[Constraint] {
        &self.expected
    }
}

pub(crate) fn constraint_compaction_case_strategy(
) -> impl Strategy<Value = ConstraintCompactionCase> {
    let cases = [
        EntityKind::Atom,
        EntityKind::Bond,
        EntityKind::DativeBond,
        EntityKind::AromaticSystem,
        EntityKind::MulticenterBond,
        EntityKind::NoncovalentBond,
        EntityKind::StereoAtom,
        EntityKind::StereoBond,
    ]
    .into_iter()
    .map(|kind| {
        let constraints = vec![
            transaction_constraint(kind, 2, 30),
            transaction_constraint(kind, 0, 10),
            transaction_constraint(kind, 1, 20),
            transaction_constraint(kind, 1, 20),
            transaction_constraint(kind, 2, 31),
        ];
        let expected = vec![
            transaction_constraint(kind, 1, 30),
            transaction_constraint(kind, 0, 20),
            transaction_constraint(kind, 0, 20),
            transaction_constraint(kind, 1, 31),
        ];
        ConstraintCompactionCase {
            kind,
            base: transaction_compaction_molecule(Constraints::from_iter(constraints)),
            expected,
        }
    })
    .collect::<Vec<_>>();
    prop::sample::select(cases)
}

pub(crate) fn consecutive_transaction_strategy(
) -> impl Strategy<Value = (MoleculeAst, Edits, Edits)> {
    (-4_i64..=4, -4_i64..=4)
        .prop_filter("successive charges must differ", |(first, second)| {
            first != second
        })
        .prop_map(|(first, second)| {
            let base = MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::from_element(Element::C)],
                ..Default::default()
            });
            let first_edits = Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: ValueAst::default(),
                    new: ValueAst::Lit(first),
                },
            }]);
            let second_edits = Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(first),
                    new: ValueAst::Lit(second),
                },
            }]);
            (base, first_edits, second_edits)
        })
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
    MoleculeAst::from_entries(MoleculeEntries {
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
pub(crate) fn overlay_transaction_strategy() -> impl Strategy<Value = (MoleculeAst, Edits)> {
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
                let mut edits = Edits::new();
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
                        )
                        .into(),
                    });
                }
                for i in (0..2).filter(|&i| con_mc[i]) {
                    edits.push(Edit::AddMoleculeConstraint {
                        constraint: Constraint::MulticenterBond(
                            MulticenterBondId(i as u32),
                            MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(4)),
                        )
                        .into(),
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
pub(crate) fn transaction_edits_strategy() -> impl Strategy<Value = (MoleculeAst, Edits)> {
    prop_oneof![
        transaction_case_strategy().prop_map(|case| (case.base(), case.edits())),
        complete_transaction_strategy(),
        overlay_transaction_strategy(),
        stable_atom_handle_trace_strategy(false).prop_map(|trace| (trace.base(), trace.edits())),
    ]
}

/// Edit sequences in the standalone grammar's exact representational shape.
///
/// Transactions accept batched atom and bond additions, while the standalone grammar represents
/// each addition as one entry. Keep those batches in `transaction_edits_strategy` for transaction
/// coverage and split only the inputs to the syntax round-trip property.
pub(crate) fn edits_dsl_strategy() -> impl Strategy<Value = Edits> {
    prop_oneof![
        transaction_case_strategy().prop_map(|case| (case.base(), case.edits())),
        prop::sample::select(complete_transaction_cases(true)),
        overlay_transaction_strategy(),
        stable_atom_handle_trace_strategy(false).prop_map(|trace| (trace.base(), trace.edits())),
    ]
    .prop_map(|(_, edits)| {
        edits
            .into_iter()
            .flat_map(|edit| match edit {
                Edit::AddAtoms { atoms } => atoms
                    .into_iter()
                    .map(|atom| Edit::AddAtoms { atoms: vec![atom] })
                    .collect(),
                Edit::AddBonds { bonds } => bonds
                    .into_iter()
                    .map(|bond| Edit::AddBonds { bonds: vec![bond] })
                    .collect(),
                edit => vec![edit],
            })
            .collect()
    })
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
            MoleculeAst::from_entries(MoleculeEntries {
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
            Correspondence::new(Vec::new(), lhs.atoms().count(), rhs.atoms().count())
                .expect("correspondence producer preserves partial-bijection invariants");
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
                MoleculeAst::from_entries(MoleculeEntries {
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
        unpaired_electrons_update_strategy(),
        aromatic_system_update_constraints_strategy(),
    )
        .prop_map(|(electrons, charge, unpaired_electrons, constraints)| {
            AromaticSystemUpdate {
                electrons,
                charge,
                unpaired_electrons,
                constraints,
            }
        })
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
        unpaired_electrons_update_strategy(),
        multicenter_bond_update_constraints_strategy(),
    )
        .prop_map(|(electrons, charge, unpaired_electrons, constraints)| {
            MulticenterBondUpdate {
                electrons,
                charge,
                unpaired_electrons,
                constraints,
            }
        })
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

    let mut deltas = Deltas::new();
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
    ReactionAst::new(lhs, deltas)
}
