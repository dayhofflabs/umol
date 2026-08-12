//! Shared proptest generators for the umol-graph-ir property suite. Domain imports
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
pub(crate) use umol_chem::element::Element;
pub(crate) use umol_edn::{read_string, Edn, FromEdn, ToEdn};
use umol_graph_core::{Correspondence, EdgeId};
pub(crate) use umol_graph_ir::dsl::{
    parse_num, AromaticSystemDsl, AromaticSystemUpdateDsl, AtomDsl, AtomUpdateDsl, BondDsl,
    BondUpdateDsl, DativeBondDsl, DativeBondParticipants, DativeBondUpdateDsl, EditsDsl,
    MetadataError, MoleculeContext, MoleculeDefaults, MoleculeDsl, MoleculeMetadata,
    MulticenterBondDsl, MulticenterBondUpdateDsl, NoncovalentBondDsl, NoncovalentBondUpdateDsl,
    NumDsl, ParseError, ReactionDefaults, ReactionDsl, ReactionMetadata, ReactionSpanDsl,
    StereoAtomConstraintDsl, StereoAtomDsl, StereoAtomParticipants, StereoAtomUpdateDsl,
    StereoBondConstraintDsl, StereoBondDsl, StereoBondParticipants, StereoBondUpdateDsl,
    StereoLigandRef,
};
pub(crate) use umol_graph_ir::ir::{
    aromatic_covalence, AddBond, ArithExpr, AromaticSystemConstraintForm,
    AromaticSystemConstraintKey, AromaticSystemConstraintsForm, AromaticSystemDelta,
    AromaticSystemFieldChange, AromaticSystemForm, AromaticSystemHandle, AromaticSystemId,
    AromaticSystemUpdate, AromaticValence, AromaticValenceForm, AsLit, AtomConstraintForm,
    AtomConstraintKey, AtomConstraintsForm, AtomDelta, AtomFieldChange, AtomForm, AtomHandle,
    AtomId, AtomUpdate, BondConstraintForm, BondConstraintKey, BondConstraintsForm, BondDelta,
    BondFieldChange, BondForm, BondHandle, BondId, BondUpdate, BooleanForm, CisTransStereoForm,
    Constraint, ConstraintEdit, Constraints, DativeBondConstraintForm, DativeBondConstraintKey,
    DativeBondConstraintsForm, DativeBondDelta, DativeBondFieldChange, DativeBondForm,
    DativeBondHandle, DativeBondId, DativeBondUpdate, Delta, Deltas, Edit, Edits,
    ElectronCountsForm, ElementForm, Entity, EntityHandle, EntityKind, Equiv, FluxionalityForm,
    FromIr, IntoIr, IsotopeMassForm, Lattice, LigandPermutation, LigandSymmetryForm, MemOp,
    Molecule, MoleculeConstraint, MoleculeCorrespondence, MoleculeEntries, MoleculeIntegrityError,
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondConstraintsForm,
    MulticenterBondDelta, MulticenterBondFieldChange, MulticenterBondForm, MulticenterBondHandle,
    MulticenterBondId, MulticenterBondUpdate, MulticenterValenceForm,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintsForm, NoncovalentBondDelta,
    NoncovalentBondFieldChange, NoncovalentBondForm, NoncovalentBondHandle, NoncovalentBondId,
    NoncovalentBondKind, NoncovalentBondKindForm, NoncovalentBondUpdate, Normalize, NumForm,
    OrientedLigandPermutation, PredExpr, Reaction, ReactionSpan, RelOp, RelationalConstraint,
    RingMembershipForm, RingScope, StereoAtomConstraintForm, StereoAtomConstraintsForm,
    StereoAtomDelta, StereoAtomFieldChange, StereoAtomForm, StereoAtomHandle, StereoAtomId,
    StereoAtomUpdate, StereoBondConstraintForm, StereoBondConstraintsForm, StereoBondDelta,
    StereoBondFieldChange, StereoBondForm, StereoBondHandle, StereoBondId, StereoBondUpdate,
    StereoConfigurationForm, StereoConfigurationUpdate, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind, StereoLigandPair, StereoLigandPosition, Stereogenicity, StereogenicityForm,
    TetrahedralStereoForm, Topicity, TopicityForm, TopicityRelationForm, TransactionError,
    UnpairedElectronsForm, UnpairedElectronsUpdate,
};
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

pub(crate) fn element_form_strategy() -> impl Strategy<Value = ElementForm> {
    prop_oneof![
        6 => element_strategy().prop_map(ElementForm::Lit),
        2 => Just(ElementForm::Undetermined),
        2 => prop::sample::subsequence(Element::all().to_vec(), 1..=118).prop_map(ElementForm::lit_set),
        1 => prop::sample::subsequence(Element::all().to_vec(), 1..=118).prop_map(ElementForm::not_set),
        1 => id_strategy().prop_map(ElementForm::var),
        1 => (id_strategy(), prop::sample::subsequence(Element::all().to_vec(), 1..=118))
            .prop_map(|(id, set)| ElementForm::var_in(id, set)),
        1 => (id_strategy(), prop::sample::subsequence(Element::all().to_vec(), 1..=118))
            .prop_map(|(id, set)| ElementForm::var_not_in(id, set)),
    ]
    .prop_map(|e| e.normalize().unwrap_or(ElementForm::Undetermined))
}

pub(crate) fn value_basic(range: RangeInclusive<i64>) -> impl Strategy<Value = NumForm> {
    prop_oneof![
        4 => Just(NumForm::Undetermined),
        4 => range.clone().prop_map(NumForm::Lit),
        1 => prop::collection::vec(range, 2..=3).prop_map(NumForm::lit_set),
        2 => arith_expr_strategy().prop_map(NumForm::arith_expr),
        2 => pred_expr_strategy().prop_map(NumForm::pred_expr),
    ]
    .prop_map(normalize_value)
}

/// Full `ArithExpr` grammar: `Lit`/`Var` leaves under `Neg`, n-ary `Sum`/
/// `Product`, and binary `Div`/`Rem`. Generated raw; `value_basic` normalizes.
fn arith_expr_strategy() -> BoxedStrategy<ArithExpr> {
    let leaf = prop_oneof![
        (-10i64..=10).prop_map(ArithExpr::Lit),
        id_strategy().prop_map(ArithExpr::Var),
    ]
    .boxed();
    leaf.prop_recursive(3, 8, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|t| ArithExpr::Neg(Box::new(t))),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(ArithExpr::Sum),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(ArithExpr::Product),
            (inner.clone(), inner.clone())
                .prop_map(|(a, b)| ArithExpr::Div(Box::new(a), Box::new(b))),
            (inner.clone(), inner).prop_map(|(a, b)| ArithExpr::Rem(Box::new(a), Box::new(b))),
        ]
        .boxed()
    })
    .boxed()
}

/// Full `PredExpr` grammar: `Rel`/`Mem` leaves over arithmetic expressions, under `Not`,
/// `And`, `Or`. Generated raw; `value_basic` normalizes (folding/NNF).
fn pred_expr_strategy() -> BoxedStrategy<PredExpr> {
    let term = arith_expr_strategy();
    let leaf = prop_oneof![
        (term.clone(), rel_op_strategy(), term.clone())
            .prop_map(|(a, op, b)| PredExpr::Rel(a, op, b)),
        (
            term,
            mem_op_strategy(),
            prop::collection::vec(-10i64..=10, 1..=3)
        )
            .prop_map(|(e, op, s)| PredExpr::Mem(e, op, s.into_iter().collect())),
    ]
    .boxed();
    leaf.prop_recursive(2, 6, 3, |inner| {
        prop_oneof![
            inner.clone().prop_map(|p| PredExpr::Not(Box::new(p))),
            prop::collection::vec(inner.clone(), 2..=3).prop_map(PredExpr::And),
            prop::collection::vec(inner, 2..=3).prop_map(PredExpr::Or),
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

/// Normalize a generated value so the property suite operates on normal
/// forms: the lattice laws and the render/parse identity compare against the
/// generated value itself, so a non-normal input would spuriously fail.
/// The unsatisfiable case is unreachable for these generators.
fn normalize_value(v: NumForm) -> NumForm {
    v.normalize().unwrap_or(NumForm::Undetermined)
}

pub(crate) fn any_num_form_strategy() -> BoxedStrategy<NumForm> {
    value_basic(-10..=10).boxed()
}

/// Possibly **non-normal** (but satisfiable) `NumForm`: unlike `value_basic`
/// it does not normalize, so it exercises the normalization-independent
/// lattice laws on raw `ArithExpr`/`PredExpr` forms. Unsatisfiable draws are filtered
/// out — on an unsatisfiable target the `matches` law's meet-derived RHS only
/// agrees with the default for satisfiable targets.
pub(crate) fn raw_num_form_strategy() -> BoxedStrategy<NumForm> {
    prop_oneof![
        2 => Just(NumForm::Undetermined),
        2 => (-10i64..=10).prop_map(NumForm::Lit),
        3 => arith_expr_strategy().prop_map(NumForm::arith_expr),
        3 => pred_expr_strategy().prop_map(NumForm::pred_expr),
    ]
    .prop_filter("satisfiable", |v| v.clone().normalize().is_ok())
    .boxed()
}

// Raw (non-normal, satisfiable) generators for the remaining
// `normalized()`-overriding leaves, to fuzz the *fold* path of the universal
// lattice laws (the normalizing strategies only ever reach the borrow path).
// Each mixes deliberately non-normal draws with the normalized strategy, then
// filters to satisfiable values (the `matches` law's RHS only agrees with the
// default on satisfiable targets).

pub(crate) fn raw_element_form_strategy() -> BoxedStrategy<ElementForm> {
    prop_oneof![
        3 => element_strategy().prop_map(|e| ElementForm::lit_set([e])),
        3 => id_strategy().prop_map(|id| ElementForm::var_in(id, Element::all().to_vec())),
        2 => element_form_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_isotope_strategy() -> BoxedStrategy<IsotopeMassForm> {
    prop_oneof![
        3 => (0u32..=250).prop_map(|m| IsotopeMassForm::lit_set(vec![m])),
        2 => isotope_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_aromatic_valence_form_strategy() -> BoxedStrategy<AromaticValenceForm> {
    prop_oneof![
        3 => raw_num_form_strategy().prop_map(AromaticValenceForm::Aromatic),
        2 => aromatic_valence_form_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_multicenter_valence_form_strategy() -> BoxedStrategy<MulticenterValenceForm> {
    prop_oneof![
        3 => raw_num_form_strategy().prop_map(MulticenterValenceForm::Multicenter),
        2 => multicenter_valence_form_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_tetrahedral_stereo_strategy() -> BoxedStrategy<TetrahedralStereoForm> {
    prop_oneof![
        Just(TetrahedralStereoForm::Undetermined),
        Just(TetrahedralStereoForm::NotStereo),
        stereo_coset_strategy().prop_map(TetrahedralStereoForm::Stereo),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_cis_trans_stereo_strategy() -> BoxedStrategy<CisTransStereoForm> {
    prop_oneof![
        Just(CisTransStereoForm::Undetermined),
        Just(CisTransStereoForm::NotStereo),
        stereo_coset_strategy().prop_map(CisTransStereoForm::Stereo),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_stereo_configuration_strategy() -> BoxedStrategy<StereoConfigurationForm> {
    prop_oneof![
        Just(StereoConfigurationForm::Undetermined),
        (stereo_atom_kind_strategy(), stereo_coset_strategy())
            .prop_map(|(kind, coset)| StereoConfigurationForm::kinded(kind, coset)),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_topicity_relation_strategy() -> BoxedStrategy<TopicityRelationForm> {
    prop_oneof![
        2 => Just(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic]))),
        2 => Just(TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Diastereotopic]))),
        3 => topicity_relation_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn raw_stereogenicity_relation_strategy() -> BoxedStrategy<StereogenicityForm> {
    prop_oneof![
        2 => Just(StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Symmetric]))),
        2 => Just(StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Stereogenic]))),
        3 => stereogenicity_relation_strategy(),
    ]
    .prop_filter("satisfiable", |x| x.clone().normalize().is_ok())
    .boxed()
}

pub(crate) fn isotope_strategy() -> impl Strategy<Value = IsotopeMassForm> {
    prop_oneof![
        3 => Just(IsotopeMassForm::Natural),
        3 => Just(IsotopeMassForm::Undetermined),
        3 => (0u32..=250).prop_map(IsotopeMassForm::Lit),
        2 => prop::collection::vec(0u32..=250, 1..=3).prop_map(IsotopeMassForm::lit_set),
        1 => id_strategy().prop_map(IsotopeMassForm::var),
        1 => (id_strategy(), prop::collection::vec(0u32..=250, 1..=3))
            .prop_map(|(id, v)| IsotopeMassForm::var_in(id, v)),
    ]
    .prop_map(|i| i.normalize().unwrap_or(IsotopeMassForm::Undetermined))
}

pub(crate) fn unpaired_electrons_strategy() -> impl Strategy<Value = UnpairedElectronsForm> {
    // The components are structurally independent; physical compatibility is
    // validated only when converting to `SpinState`.
    (value_basic(0..=6), value_basic(1..=7)).prop_map(|(count, multiplicity)| {
        UnpairedElectronsForm {
            count,
            multiplicity,
        }
    })
}

pub(crate) fn raw_unpaired_electrons_strategy() -> impl Strategy<Value = UnpairedElectronsForm> {
    (raw_num_form_strategy(), raw_num_form_strategy()).prop_map(|(count, multiplicity)| {
        UnpairedElectronsForm {
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

/// `UnpairedElectronsForm` with at least one of `count` / `multiplicity` not
/// `Undetermined`. Used inside `MoleculeConstraint::UnpairedElectronCoupling` and similar
/// where a fully vacuous unpaired-electron state would elide on render.
pub(crate) fn non_vacuous_unpaired_electrons_strategy(
) -> impl Strategy<Value = UnpairedElectronsForm> {
    (value_basic(0..=6), value_basic(1..=7))
        .prop_map(|(u, m)| UnpairedElectronsForm {
            count: u,
            multiplicity: m,
        })
        .prop_filter("non-vacuous unpaired-electron state", |s| {
            !s.is_undetermined()
        })
}

/// Simple value strategy used inside constraint values: `Undetermined`,
/// `Lit`, and `LitSet`. No symbolic `ArithExpr`/`PredExpr` — the constraint
/// formatters route to `fmt_num_field_required` / `fmt_ring_count` / the
/// various `#r` blocks, and an `ArithExpr(Lit(n))` would render to a pure integer
/// that the parser then re-reads as a plain `Lit`, breaking roundtrip. The
/// molecule-level EDN tests cover symbolic values on constraint values through
/// the tree-based path, so the gap is contained.
pub(crate) fn constraint_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = NumForm> {
    prop_oneof![
        3 => Just(NumForm::Undetermined),
        3 => range.clone().prop_map(NumForm::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            NumForm::lit_set(v)
        }),
    ]
}

/// `Lit`/`Set` only — still used by the ring-size strategies where
/// `Undetermined` on the inner value collapses into a dropped constraint
/// in the entity-level formatter (see vacuous `RingMembership(_, Undetermined)`, intentionally dropped).
pub(crate) fn constraint_inner_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = NumForm> {
    prop_oneof![
        range.clone().prop_map(NumForm::Lit),
        prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            NumForm::lit_set(v)
        }),
    ]
}

/// `AromaticValenceForm::Undetermined` is vacuous (renders empty per the
/// canonical-rendering rule). Excluded so the strategy stays inside the
/// render → reparse identity.
pub(crate) fn aromatic_valence_form_strategy() -> impl Strategy<Value = AromaticValenceForm> {
    prop_oneof![
        Just(AromaticValenceForm::NotAromatic),
        constraint_value_strategy(0..=6).prop_map(AromaticValenceForm::Aromatic),
    ]
    .prop_map(|v| v.normalize().unwrap_or(AromaticValenceForm::Undetermined))
}

pub(crate) fn multicenter_valence_form_strategy() -> impl Strategy<Value = MulticenterValenceForm> {
    prop_oneof![
        Just(MulticenterValenceForm::NotMulticenter),
        constraint_value_strategy(0..=6).prop_map(MulticenterValenceForm::Multicenter),
    ]
    .prop_map(|v| {
        v.normalize()
            .unwrap_or(MulticenterValenceForm::Undetermined)
    })
}

/// Atom constraints route through `fmt_num_field_required` (or
/// `fmt_ring_count` for `#R`), which elide vacuous (Undetermined) payloads
/// per the canonical-rendering rule. Generators excluding `Undetermined`
/// keep the render → reparse identity intact.
pub(crate) fn atom_constraint_strategy() -> BoxedStrategy<AtomConstraintForm> {
    prop_oneof![
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::Valence),
        constraint_inner_value_strategy(0..=8).prop_map(AtomConstraintForm::TotalValence),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraintForm::DonatedPairs),
        constraint_inner_value_strategy(0..=4).prop_map(AtomConstraintForm::AcceptedPairs),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::Degree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::TotalDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::RingDegree),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::RingValence),
        constraint_inner_value_strategy(0..=6).prop_map(AtomConstraintForm::TotalHydrogens),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| AtomConstraintForm::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6))
            .prop_map(|(s, count)| AtomConstraintForm::ring_membership(RingScope::Size(s), count)),
        aromatic_valence_form_strategy().prop_map(AtomConstraintForm::AromaticValence),
        multicenter_valence_form_strategy().prop_map(AtomConstraintForm::MulticenterValence),
        tetrahedral_stereo_strategy().prop_map(AtomConstraintForm::TetrahedralStereo),
    ]
    .boxed()
}

pub(crate) fn atom_constraints_strategy() -> impl Strategy<Value = AtomConstraintsForm> {
    prop::collection::vec(atom_constraint_strategy(), 0..=3).prop_map(|list| {
        let mut cs = AtomConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn atom_update_constraints_strategy() -> impl Strategy<Value = AtomConstraintsForm> {
    prop::collection::vec(
        prop_oneof![
            atom_constraint_strategy(),
            atom_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(AtomConstraintsForm::from_iter)
}

pub(crate) fn bond_constraint_strategy() -> BoxedStrategy<BondConstraintForm> {
    prop_oneof![
        any::<bool>().prop_map(|b| BondConstraintForm::Aromatic(BooleanForm::Lit(b))),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| BondConstraintForm::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6))
            .prop_map(|(s, count)| BondConstraintForm::ring_membership(RingScope::Size(s), count)),
        cis_trans_stereo_strategy().prop_map(BondConstraintForm::CisTransStereo),
    ]
    .boxed()
}

pub(crate) fn bond_constraints_strategy() -> impl Strategy<Value = BondConstraintsForm> {
    prop::collection::vec(bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = BondConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn bond_update_constraints_strategy() -> impl Strategy<Value = BondConstraintsForm> {
    prop::collection::vec(
        prop_oneof![
            bond_constraint_strategy(),
            bond_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=2,
    )
    .prop_map(BondConstraintsForm::from_iter)
}

pub(crate) fn dative_bond_constraint_strategy() -> BoxedStrategy<DativeBondConstraintForm> {
    prop_oneof![
        any::<bool>().prop_map(|b| DativeBondConstraintForm::Aromatic(BooleanForm::Lit(b))),
        constraint_inner_value_strategy(0..=6)
            .prop_map(|v| DativeBondConstraintForm::ring_membership(RingScope::All, v)),
        (3u8..=10, constraint_inner_value_strategy(0..=6)).prop_map(|(s, count)| {
            DativeBondConstraintForm::ring_membership(RingScope::Size(s), count)
        }),
    ]
    .boxed()
}

pub(crate) fn dative_bond_constraints_strategy() -> impl Strategy<Value = DativeBondConstraintsForm>
{
    prop::collection::vec(dative_bond_constraint_strategy(), 0..=2).prop_map(|list| {
        let mut cs = DativeBondConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn dative_bond_update_constraints_strategy(
) -> impl Strategy<Value = DativeBondConstraintsForm> {
    prop::collection::vec(
        prop_oneof![
            dative_bond_constraint_strategy(),
            dative_bond_constraint_strategy().prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=2,
    )
    .prop_map(DativeBondConstraintsForm::from_iter)
}

prop_compose! {
    pub(crate) fn atom_form_strategy()
    (
        element in element_form_strategy(),
        isotope in isotope_strategy(),
        charge in value_basic(-2..=2),
        implicit_hydrogens in value_basic(0..=4),
        lone_pairs in value_basic(0..=4),
        unpaired_electrons in unpaired_electrons_strategy(),
        constraints in atom_constraints_strategy(),
    ) -> AtomForm {
        AtomForm {
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
        element in prop::option::of(element_form_strategy()),
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
    pub(crate) fn bond_form_strategy()
    (
        order in value_basic(1..=4),
        charge in value_basic(-1..=1),
        unpaired_electrons in unpaired_electrons_strategy(),
        constraints in bond_constraints_strategy(),
    ) -> BondForm {
        BondForm {
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

/// `BondForm` shapes that render to bond keyword shorthands per spec §7.6:
/// `:single`, `:double`, `:triple`, `:quadruple`, plus `:aromatic` (an
/// order-1 bond with the inline `Aromatic` flag).
pub(crate) fn canonical_keyword_bond_strategy() -> impl Strategy<Value = BondForm> {
    prop_oneof![
        Just(BondForm::new(NumForm::Lit(1))),
        Just(BondForm::new(NumForm::Lit(2))),
        Just(BondForm::new(NumForm::Lit(3))),
        Just(BondForm::new(NumForm::Lit(4))),
        Just(
            BondForm::new(NumForm::Lit(1))
                .with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
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

pub(crate) fn dative_bond_strategy() -> impl Strategy<Value = DativeBondForm> {
    // Order is sampled from the small literal range that the DSL keyword
    // shorthands cover (`:single` / `:double` / `:triple`), keeping
    // canonical-form roundtrip exercised across haptic-pair counts.
    let order_strategy = prop_oneof![
        Just(NumForm::Lit(1)),
        Just(NumForm::Lit(2)),
        Just(NumForm::Lit(3)),
        Just(NumForm::Undetermined),
    ];
    (order_strategy, dative_bond_constraints_strategy())
        .prop_map(|(order, constraints)| DativeBondForm { order, constraints })
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
/// emits `None` half the time, otherwise wraps a `NumForm::Lit` or
/// `Set`. `Undetermined` is excluded because it has no canonical
/// surface form in the entity-string `#e<n>` field — `#e*` is admitted on
/// parse but the renderer omits the predicate entirely, breaking
/// roundtrip.
pub(crate) fn optional_aromatic_electron_count(
) -> impl Strategy<Value = AromaticSystemConstraintsForm> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=12)).prop_map(|opt| {
        let mut cs = AromaticSystemConstraintsForm::new();
        if let Some(v) = opt {
            cs.set(AromaticSystemConstraintForm::ElectronCount(v));
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn aromatic_system_update_constraints_strategy(
) -> impl Strategy<Value = AromaticSystemConstraintsForm> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            electron_count_value_strategy(0..=12),
            Just(NumForm::Undetermined),
        ],
    )
    .prop_map(|value| {
        value
            .map(AromaticSystemConstraintForm::ElectronCount)
            .map(AromaticSystemConstraintsForm::from)
            .unwrap_or_default()
    })
}

pub(crate) fn optional_multicenter_electron_count(
) -> impl Strategy<Value = MulticenterBondConstraintsForm> {
    prop::option::weighted(0.5, electron_count_value_strategy(0..=8)).prop_map(|opt| {
        let mut cs = MulticenterBondConstraintsForm::new();
        if let Some(v) = opt {
            cs.set(MulticenterBondConstraintForm::ElectronCount(v));
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn multicenter_bond_update_constraints_strategy(
) -> impl Strategy<Value = MulticenterBondConstraintsForm> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            electron_count_value_strategy(0..=8),
            Just(NumForm::Undetermined),
        ],
    )
    .prop_map(|value| {
        value
            .map(MulticenterBondConstraintForm::ElectronCount)
            .map(MulticenterBondConstraintsForm::from)
            .unwrap_or_default()
    })
}

pub(crate) fn electron_count_value_strategy(
    range: RangeInclusive<i64>,
) -> impl Strategy<Value = NumForm> {
    prop_oneof![
        3 => range.clone().prop_map(NumForm::Lit),
        1 => prop::collection::vec(range, 1..=3).prop_map(|mut v| {
            v.sort_unstable();
            v.dedup();
            NumForm::lit_set(v)
        }),
    ]
}

/// Leaf strategy: `Undetermined` or a concrete `Lit` count vector (length 1–4).
pub(crate) fn electron_counts_form_strategy() -> impl Strategy<Value = ElectronCountsForm> {
    prop_oneof![
        Just(ElectronCountsForm::Undetermined),
        prop::collection::vec(0i64..=8, 1..=4).prop_map(ElectronCountsForm::Lit),
    ]
}

/// Stand-alone strategy for entity-string roundtrip tests. `electrons` is
/// `Undetermined` because the entity string carries no per-atom data; the
/// `ElectronCount` constraint is exercised here via `#e<n>`.
pub(crate) fn aromatic_system_form_strategy() -> impl Strategy<Value = AromaticSystemForm> {
    (value_basic(-2..=2), optional_aromatic_electron_count()).prop_map(|(charge, constraints)| {
        AromaticSystemForm {
            electrons: ElectronCountsForm::Undetermined,
            charge,
            unpaired_electrons: UnpairedElectronsForm::default(),
            constraints,
        }
    })
}

pub(crate) fn aromatic_system_patch_form_strategy() -> impl Strategy<Value = AromaticSystemForm> {
    (
        electron_counts_form_strategy(),
        value_basic(-2..=2),
        unpaired_electrons_strategy(),
        optional_aromatic_electron_count(),
    )
        .prop_map(
            |(electrons, charge, unpaired_electrons, constraints)| AromaticSystemForm {
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
        electrons in prop::option::of(electron_counts_form_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in aromatic_system_update_constraints_strategy(),
    ) -> AromaticSystemUpdate {
        AromaticSystemUpdate { electrons, charge, unpaired_electrons, constraints }
    }
}

/// Atom-count-aware variant: generates an `AromaticSystemForm` whose
/// `electrons` `Lit` vector has exactly `atom_count` entries. Includes an
/// optional `ElectronCount` constraint so the molecule-level prop tests
/// exercise both the per-atom counts and the asserted total in the same pass.
pub(crate) fn aromatic_system_form_for(
    atom_count: usize,
) -> impl Strategy<Value = AromaticSystemForm> {
    (
        value_basic(-2..=2),
        prop::collection::vec(0i64..=2, atom_count),
        optional_aromatic_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| AromaticSystemForm {
            electrons: ElectronCountsForm::Lit(electrons),
            charge,
            unpaired_electrons: UnpairedElectronsForm::default(),
            constraints,
        })
}

pub(crate) fn multicenter_bond_form_strategy() -> impl Strategy<Value = MulticenterBondForm> {
    (value_basic(-2..=2), optional_multicenter_electron_count()).prop_map(
        |(charge, constraints)| MulticenterBondForm {
            electrons: ElectronCountsForm::Undetermined,
            charge,
            unpaired_electrons: UnpairedElectronsForm::default(),
            constraints,
        },
    )
}

pub(crate) fn multicenter_bond_patch_form_strategy() -> impl Strategy<Value = MulticenterBondForm> {
    (
        electron_counts_form_strategy(),
        value_basic(-2..=2),
        unpaired_electrons_strategy(),
        optional_multicenter_electron_count(),
    )
        .prop_map(
            |(electrons, charge, unpaired_electrons, constraints)| MulticenterBondForm {
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
        electrons in prop::option::of(electron_counts_form_strategy()),
        charge in prop::option::of(value_basic(-2..=2)),
        unpaired_electrons in unpaired_electrons_update_strategy(),
        constraints in multicenter_bond_update_constraints_strategy(),
    ) -> MulticenterBondUpdate {
        MulticenterBondUpdate { electrons, charge, unpaired_electrons, constraints }
    }
}

pub(crate) fn multicenter_bond_form_for(
    atom_count: usize,
) -> impl Strategy<Value = MulticenterBondForm> {
    (
        value_basic(-2..=2),
        prop::collection::vec(0i64..=2, atom_count),
        optional_multicenter_electron_count(),
    )
        .prop_map(|(charge, electrons, constraints)| MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(electrons),
            charge,
            unpaired_electrons: UnpairedElectronsForm::default(),
            constraints,
        })
}

pub(crate) fn noncovalent_bond_kind_form_strategy() -> impl Strategy<Value = NoncovalentBondKindForm>
{
    prop_oneof![
        Just(NoncovalentBondKindForm::Undetermined),
        prop::sample::select(NONCOVALENT_KINDS).prop_map(NoncovalentBondKindForm::Lit),
    ]
}

pub(crate) fn noncovalent_bond_constraint_strategy() -> BoxedStrategy<NoncovalentBondConstraintForm>
{
    any::<bool>()
        .prop_map(|b| NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Lit(b)))
        .boxed()
}

pub(crate) fn noncovalent_bond_constraints_strategy(
) -> impl Strategy<Value = NoncovalentBondConstraintsForm> {
    prop::collection::vec(noncovalent_bond_constraint_strategy(), 0..=1).prop_map(|list| {
        let mut cs = NoncovalentBondConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn noncovalent_bond_update_constraints_strategy(
) -> impl Strategy<Value = NoncovalentBondConstraintsForm> {
    prop::option::weighted(
        0.5,
        prop_oneof![
            any::<bool>().prop_map(|value| {
                NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Lit(value))
            }),
            Just(NoncovalentBondConstraintForm::Intramolecular(
                BooleanForm::Undetermined,
            )),
        ],
    )
    .prop_map(|constraint| {
        constraint
            .map(NoncovalentBondConstraintsForm::from)
            .unwrap_or_default()
    })
}

pub(crate) fn noncovalent_bond_form_strategy() -> impl Strategy<Value = NoncovalentBondForm> {
    (
        prop::sample::select(NONCOVALENT_KINDS),
        noncovalent_bond_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondForm {
            kind: NoncovalentBondKindForm::Lit(kind),
            constraints,
        })
}

pub(crate) fn noncovalent_bond_patch_form_strategy() -> impl Strategy<Value = NoncovalentBondForm> {
    (
        noncovalent_bond_kind_form_strategy(),
        noncovalent_bond_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondForm { kind, constraints })
}

pub(crate) fn noncovalent_bond_update_strategy() -> impl Strategy<Value = NoncovalentBondUpdate> {
    (
        prop::option::of(noncovalent_bond_kind_form_strategy()),
        noncovalent_bond_update_constraints_strategy(),
    )
        .prop_map(|(kind, constraints)| NoncovalentBondUpdate { kind, constraints })
}

/// Coset forms that round-trip through both the entity attribute string and the
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
pub(crate) fn tetrahedral_stereo_strategy() -> impl Strategy<Value = TetrahedralStereoForm> {
    prop_oneof![
        Just(TetrahedralStereoForm::NotStereo),
        stereo_coset_strategy().prop_map(TetrahedralStereoForm::Stereo),
    ]
    .prop_map(|s| s.normalize().unwrap_or(TetrahedralStereoForm::Undetermined))
}

pub(crate) fn tetrahedral_stereo_lattice_strategy() -> impl Strategy<Value = TetrahedralStereoForm>
{
    prop_oneof![
        Just(TetrahedralStereoForm::Undetermined),
        tetrahedral_stereo_strategy(),
    ]
}

pub(crate) fn cis_trans_stereo_strategy() -> impl Strategy<Value = CisTransStereoForm> {
    prop_oneof![
        Just(CisTransStereoForm::NotStereo),
        stereo_coset_strategy().prop_map(CisTransStereoForm::Stereo),
    ]
    .prop_map(|s| s.normalize().unwrap_or(CisTransStereoForm::Undetermined))
}

pub(crate) fn cis_trans_stereo_lattice_strategy() -> impl Strategy<Value = CisTransStereoForm> {
    prop_oneof![
        Just(CisTransStereoForm::Undetermined),
        cis_trans_stereo_strategy(),
    ]
}

/// `StereoConfigurationForm` over the atom geometry kinds, including the kindless
/// `Undetermined` top.
pub(crate) fn stereo_configuration_lattice_strategy(
) -> impl Strategy<Value = StereoConfigurationForm> {
    prop_oneof![
        Just(StereoConfigurationForm::Undetermined),
        (stereo_atom_kind_strategy(), stereo_coset_strategy())
            .prop_map(|(kind, coset)| StereoConfigurationForm::kinded(kind, coset)),
    ]
    .prop_map(|c| {
        c.normalize()
            .unwrap_or(StereoConfigurationForm::Undetermined)
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
pub(crate) fn topicity_relation_strategy() -> impl Strategy<Value = TopicityRelationForm> {
    prop_oneof![
        Just(TopicityRelationForm::Lit(Topicity::Homotopic)),
        Just(TopicityRelationForm::Lit(Topicity::Enantiotopic)),
        Just(TopicityRelationForm::Lit(Topicity::Diastereotopic)),
        Just(TopicityRelationForm::NotSet(BTreeSet::from([
            Topicity::Homotopic
        ]))),
        Just(TopicityRelationForm::NotSet(BTreeSet::from([
            Topicity::Enantiotopic
        ]))),
        Just(TopicityRelationForm::NotSet(BTreeSet::from([
            Topicity::Diastereotopic
        ]))),
    ]
}

pub(crate) fn stereogenicity_relation_strategy() -> impl Strategy<Value = StereogenicityForm> {
    prop_oneof![
        Just(StereogenicityForm::Lit(Stereogenicity::Symmetric)),
        Just(StereogenicityForm::Lit(Stereogenicity::Prochiral)),
        Just(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
        Just(StereogenicityForm::NotSet(BTreeSet::from([
            Stereogenicity::Symmetric
        ]))),
        Just(StereogenicityForm::NotSet(BTreeSet::from([
            Stereogenicity::Prochiral
        ]))),
        Just(StereogenicityForm::NotSet(BTreeSet::from([
            Stereogenicity::Stereogenic
        ]))),
    ]
}

/// Topicity relations spanning the full lattice: the non-vacuous singletons /
/// complements plus `Undetermined` (top).
pub(crate) fn topicity_relation_lattice_strategy() -> impl Strategy<Value = TopicityRelationForm> {
    prop_oneof![
        Just(TopicityRelationForm::Undetermined),
        topicity_relation_strategy(),
    ]
}

pub(crate) fn stereogenicity_relation_lattice_strategy() -> impl Strategy<Value = StereogenicityForm>
{
    prop_oneof![
        Just(StereogenicityForm::Undetermined),
        stereogenicity_relation_strategy(),
    ]
}

pub(crate) fn ligand_symmetry_strategy(degree: usize) -> impl Strategy<Value = LigandSymmetryForm> {
    (
        permutation_strategy(degree),
        orientation_strategy(),
        any::<bool>(),
    )
        .prop_map(|(permutation, orientation, invariant)| LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(permutation),
                orientation,
            },
            invariant: BooleanForm::Lit(invariant),
        })
}

pub(crate) fn fluxionality_strategy(degree: usize) -> impl Strategy<Value = FluxionalityForm> {
    (permutation_strategy(degree), any::<bool>()).prop_map(|(permutation, active)| {
        FluxionalityForm {
            permutation: LigandPermutation(permutation),
            active: BooleanForm::Lit(active),
        }
    })
}

pub(crate) fn topicity_strategy(degree: usize) -> impl Strategy<Value = TopicityForm> {
    (
        ligand_pair_strategy(degree),
        topicity_relation_lattice_strategy(),
    )
        .prop_map(|(pair, relation)| TopicityForm { pair, relation })
}

/// Normalized, fiber-spanning `RingMembershipForm`: the `scope` varies (`All` and
/// `Size(3..=10)`) so a value triple lands in different fibers, exercising the
/// cross-scope `meet` → `None` / `join` → `Err(NoJoin)` path.
pub(crate) fn ring_membership_lattice_strategy() -> impl Strategy<Value = RingMembershipForm> {
    prop_oneof![
        constraint_value_strategy(0..=6)
            .prop_map(|count| RingMembershipForm::new(RingScope::All, count)),
        (3u8..=10, constraint_value_strategy(0..=6))
            .prop_map(|(size, count)| RingMembershipForm::new(RingScope::Size(size), count)),
    ]
    .prop_map(|membership| {
        membership
            .normalize()
            .expect("non-empty count strategy never contradicts")
    })
}

/// Universal lattice laws — hold for **any** inputs (normalized or not): meet/join
/// commutativity and associativity, `matches` ⇔ meet-derived, and the
/// Lattice→Normalize correspondence that `meet`/`join` land in normal form.
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
    prop_assert_eq!(a.matches(b), a.meet(b) == b.clone().normalize().ok());
    if let Some(m) = a.meet(b) {
        prop_assert_eq!(m.clone().normalize(), Ok(m));
    }
    if let Ok(j) = a.join(b) {
        prop_assert_eq!(j.clone().normalize(), Ok(j));
    }
    // `normalized()` (the borrow fast-path) agrees with `normalize()`.
    prop_assert_eq!(
        a.normalized().map(|c| c.into_owned()),
        a.clone().normalize()
    );
    // `equiv` is equality of normal forms.
    prop_assert_eq!(
        a.equiv(b),
        a.clone().normalize().ok() == b.clone().normalize().ok()
    );
    Ok(())
}

/// Lattice laws that assume **normalized** inputs: each input is a `normalize`
/// fixpoint, plus idempotence and absorption (whose RHS is the input verbatim,
/// which only holds when the input is already normalized).
pub(crate) fn assert_normalized_lattice_laws<L: Lattice + Debug>(
    a: &L,
    b: &L,
    c: &L,
) -> Result<(), TestCaseError> {
    prop_assert_eq!(a.clone().normalize(), Ok(a.clone()));
    prop_assert_eq!(b.clone().normalize(), Ok(b.clone()));
    prop_assert_eq!(c.clone().normalize(), Ok(c.clone()));
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
                        $constraint::LigandSymmetry(LigandSymmetryForm {
                            permutation: OrientedLigandPermutation {
                                permutation: LigandPermutation(permutation),
                                orientation,
                            },
                            invariant: BooleanForm::Lit(invariant),
                        })
                    }),
                (permutation_strategy(degree), any::<bool>()).prop_map(|(permutation, active)| {
                    $constraint::Fluxionality(FluxionalityForm {
                        permutation: LigandPermutation(permutation),
                        active: BooleanForm::Lit(active),
                    })
                }),
                (ligand_pair_strategy(degree), topicity_relation_strategy()).prop_map(
                    |(pair, relation)| $constraint::Topicity(TopicityForm { pair, relation })
                ),
                stereogenicity_relation_strategy().prop_map(|rel| $constraint::Stereogenicity(rel)),
            ]
            .boxed()
        }
    };
}

stereo_constraint_strategy! { stereo_atom_constraint_strategy, StereoAtomConstraintForm }
stereo_constraint_strategy! { stereo_bond_constraint_strategy, StereoBondConstraintForm }

pub(crate) fn stereo_atom_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoAtomConstraintsForm> {
    prop::collection::vec(stereo_atom_constraint_strategy(kind), 0..=3).prop_map(|list| {
        let mut cs = StereoAtomConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn stereo_bond_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoBondConstraintsForm> {
    prop::collection::vec(stereo_bond_constraint_strategy(kind), 0..=3).prop_map(|list| {
        let mut cs = StereoBondConstraintsForm::new();
        for c in list {
            cs.set(c);
        }
        cs.normalize().unwrap_or_default()
    })
}

pub(crate) fn stereo_atom_form_strategy() -> impl Strategy<Value = StereoAtomForm> {
    stereo_atom_kind_strategy().prop_flat_map(|kind| {
        (
            stereo_coset_for_kind(kind),
            stereo_atom_constraints_strategy(kind),
        )
            .prop_map(move |(coset, cs)| StereoAtomForm::new(kind, coset).with_constraints(cs))
    })
}

pub(crate) fn stereo_atom_update_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoAtomConstraintsForm> {
    prop::collection::vec(
        prop_oneof![
            stereo_atom_constraint_strategy(kind),
            stereo_atom_constraint_strategy(kind)
                .prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(StereoAtomConstraintsForm::from_iter)
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

pub(crate) fn stereo_bond_form_strategy() -> impl Strategy<Value = StereoBondForm> {
    (
        stereo_coset_for_kind(StereoKind::CisTrans),
        stereo_bond_constraints_strategy(StereoKind::CisTrans),
    )
        .prop_map(|(coset, cs)| {
            StereoBondForm::new(StereoKind::CisTrans, coset).with_constraints(cs)
        })
}

pub(crate) fn stereo_bond_update_constraints_strategy(
    kind: StereoKind,
) -> impl Strategy<Value = StereoBondConstraintsForm> {
    prop::collection::vec(
        prop_oneof![
            stereo_bond_constraint_strategy(kind),
            stereo_bond_constraint_strategy(kind)
                .prop_map(|constraint| constraint.as_undetermined()),
        ],
        0..=3,
    )
    .prop_map(StereoBondConstraintsForm::from_iter)
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
) -> BoxedStrategy<Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>> {
    if atom_count == 0 || max == 0 {
        return Just(Vec::new()).boxed();
    }
    let entry = stereo_atom_kind_strategy().prop_flat_map(move |kind| {
        (
            prop::collection::vec(
                (stereo_ligand_kind_strategy(), 0..atom_count as u32),
                kind.degree(),
            ),
            stereo_coset_for_kind(kind),
            stereo_atom_constraints_strategy(kind),
        )
            .prop_map(move |(ligands, coset, constraints)| {
                (
                    ligands,
                    StereoAtomForm::new(kind, coset).with_constraints(constraints),
                )
            })
    });
    prop::collection::vec(entry, 0..=max)
        .prop_map(|entries| {
            entries
                .into_iter()
                .enumerate()
                .map(|(i, (ligand_specs, attributes))| {
                    let site = AtomId(i as u32);
                    let mut actual_atoms = HashSet::new();
                    let ligands = ligand_specs
                        .into_iter()
                        .map(|(kind, a)| {
                            let atom = AtomId(a);
                            match kind {
                                StereoLigandKind::Atom
                                    if atom != site && actual_atoms.insert(atom) =>
                                {
                                    StereoLigand::new(atom, kind)
                                }
                                StereoLigandKind::Atom => {
                                    StereoLigand::new(site, StereoLigandKind::ImplicitHydrogen)
                                }
                                _ => StereoLigand::new(site, kind),
                            }
                        })
                        .collect();
                    (site, ligands, attributes)
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
) -> BoxedStrategy<Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>> {
    if atom_count == 0 || bond_count == 0 || max == 0 {
        return Just(Vec::new()).boxed();
    }
    let kind = StereoKind::CisTrans;
    let entry = (
        prop::collection::vec(
            (stereo_ligand_kind_strategy(), 0..atom_count as u32),
            kind.degree(),
        ),
        stereo_coset_for_kind(kind),
        stereo_bond_constraints_strategy(kind),
    )
        .prop_map(move |(ligands, coset, constraints)| {
            (
                ligands,
                StereoBondForm::new(kind, coset).with_constraints(constraints),
            )
        });
    prop::collection::vec(entry, 0..=max)
        .prop_map(|entries| {
            entries
                .into_iter()
                .enumerate()
                .map(|(i, (ligand_specs, attributes))| {
                    let site = BondId(i as u32);
                    let mut actual_atoms = HashSet::new();
                    let ligands = ligand_specs
                        .into_iter()
                        .map(|(kind, a)| {
                            let atom = AtomId(a);
                            match kind {
                                StereoLigandKind::Atom if actual_atoms.insert(atom) => {
                                    StereoLigand::new(atom, kind)
                                }
                                StereoLigandKind::Atom => {
                                    StereoLigand::new(atom, StereoLigandKind::ImplicitHydrogen)
                                }
                                _ => StereoLigand::new(atom, kind),
                            }
                        })
                        .collect();
                    (site, ligands, attributes)
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
            let atoms = prop::collection::vec(atom_form_strategy(), atom_count);
            let edges = edge_set_strategy(atom_count);
            let bond_data = prop::collection::vec(bond_form_strategy(), 0..=8);
            (Just(atom_count), atoms, edges, bond_data)
        })
        .prop_flat_map(|(atom_count, atoms, edges, bond_pool)| {
            // Truncate bond pool to the number of edges generated.
            let bond_count = edges.len();
            let bonds: Vec<BondForm> = bond_pool
                .into_iter()
                .chain(repeat_with(|| BondForm::from_order(1)))
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
                        (Just(atoms), aromatic_system_form_for(n))
                    },
                ),
                0..=aromatic_count_max,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_form_for(n))
                    },
                ),
                0..=multicenter_count_max,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_form_strategy(),
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
                let mut dative_incidences = HashSet::new();
                let dative_triples: Vec<_> = datives
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [donor, acceptor]
                            if donor != acceptor
                                && dative_incidences.insert((*acceptor, *donor)) =>
                        {
                            Some((vec![*donor], *acceptor, data))
                        }
                        _ => None,
                    })
                    .collect();
                let mut aromatic_atoms = HashSet::<AtomId>::new();
                let aromatic_entries: Vec<_> = aromatics
                    .into_iter()
                    .filter(|(atoms, _)| {
                        atoms.len() >= 3
                            && atoms.iter().all(|atom| !aromatic_atoms.contains(atom))
                            && {
                                aromatic_atoms.extend(atoms.iter().copied());
                                true
                            }
                    })
                    .collect();
                let mut multicenter_sets = HashSet::new();
                let multicenter_entries: Vec<_> = multicenters
                    .into_iter()
                    .filter(|(atoms, _)| {
                        atoms.len() >= 3
                            && multicenter_sets
                                .insert(atoms.iter().copied().collect::<BTreeSet<_>>())
                    })
                    .collect();
                let mut noncovalent_pairs = HashSet::new();
                let noncovalent_triples: Vec<_> = noncovalents
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b]
                            if a != b
                                && noncovalent_pairs.insert(if a < b {
                                    [*a, *b]
                                } else {
                                    [*b, *a]
                                }) =>
                        {
                            Some((*a, *b, data))
                        }
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

pub(crate) fn molecule_strategy() -> impl Strategy<Value = Molecule> {
    molecule_entries_strategy().prop_map(Molecule::from_entries)
}

/// Per-entity counts for a generated `Molecule`. Carried into the
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

    fn from_ir(molecule: &Molecule) -> Self {
        Self {
            atom: molecule.atoms().count(),
            bond: molecule.bonds().count(),
            dative: molecule.dative_bonds().count(),
            aromatic: molecule.aromatic_systems().count(),
            multicenter: molecule.multicenter_bonds().count(),
            noncovalent: molecule.noncovalent_bonds().count(),
            stereo_atom: molecule.stereo_atoms().count(),
            stereo_bond: molecule.stereo_bonds().count(),
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

        // Constraint::Atom carrying any AtomConstraintForm variant.
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
                Constraint::AromaticSystem(system, AromaticSystemConstraintForm::ElectronCount(v))
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
                Constraint::MulticenterBond(bond, MulticenterBondConstraintForm::ElectronCount(v))
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

    if choices.is_empty() {
        return Just(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: None,
        }))
        .boxed();
    }
    prop::strategy::Union::new(choices).boxed()
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

pub(crate) fn molecule_with_constraints_strategy() -> impl Strategy<Value = Molecule> {
    molecule_entries_with_constraints_strategy().prop_map(Molecule::from_entries)
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

pub(crate) fn molecule_with_atom_subset_strategy() -> impl Strategy<Value = (Molecule, Vec<AtomId>)>
{
    molecule_structurally_unambiguous_strategy().prop_flat_map(|molecule| {
        let atom_count = molecule.atoms().count();
        (
            Just(molecule),
            prop::collection::vec(any::<bool>(), atom_count),
        )
            .prop_map(|(molecule, keep)| {
                let atoms = keep
                    .into_iter()
                    .enumerate()
                    .filter_map(|(index, keep)| keep.then_some(AtomId(index as u32)))
                    .collect();
                (molecule, atoms)
            })
    })
}

pub(crate) fn molecule_with_removals_strategy(
) -> impl Strategy<Value = (Molecule, Vec<AtomId>, Vec<BondId>)> {
    molecule_strategy().prop_flat_map(|molecule| {
        let atom_count = molecule.atoms().count();
        let bond_count = molecule.bonds().count();
        (
            Just(molecule),
            prop::collection::vec(any::<bool>(), atom_count),
            prop::collection::vec(any::<bool>(), bond_count),
        )
            .prop_map(|(molecule, atom_mask, bond_mask)| {
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
                (molecule, atoms, bonds)
            })
    })
}

pub(crate) fn molecule_structurally_unambiguous_strategy() -> impl Strategy<Value = Molecule> {
    molecule_entries_structurally_unambiguous_strategy().prop_map(Molecule::from_entries)
}

pub(crate) fn molecule_entries_structurally_unambiguous_strategy(
) -> impl Strategy<Value = MoleculeEntries> {
    molecule_entries_strategy().prop_filter(
        "entity incidence identifies at most one entity of each family",
        |entries| molecule_entity_incidence_is_unique(&Molecule::from_entries(entries.clone())),
    )
}

fn molecule_entity_incidence_is_unique(molecule: &Molecule) -> bool {
    all_unique(
        molecule
            .bonds()
            .iter()
            .map(|bond| sorted_pair(bond.atom_ids())),
    ) && all_unique(
        molecule
            .dative_bonds()
            .iter()
            .map(|dative| (dative.acceptor_id(), sorted(dative.donor_ids().collect()))),
    ) && all_unique(
        molecule
            .aromatic_systems()
            .iter()
            .map(|aromatic| sorted(aromatic.atom_ids().collect())),
    ) && all_unique(
        molecule
            .multicenter_bonds()
            .iter()
            .map(|multicenter| sorted(multicenter.atom_ids().collect())),
    ) && all_unique(
        molecule
            .noncovalent_bonds()
            .iter()
            .map(|noncovalent| sorted_pair(noncovalent.atom_ids())),
    ) && all_unique(
        molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| (stereo.site_id(), sorted(stereo.ligand_frame()))),
    ) && all_unique(
        molecule
            .stereo_bonds()
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

/// Generate a `MoleculeMetadata` populated for a molecule of the given counts. Entity
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
                    meta.add_atom_alias(format!("al{i}"), AtomForm::from_element(*element))
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
    molecule_with_constraints_strategy().prop_flat_map(|molecule| {
        let counts = ConstraintCounts::from_ir(&molecule);
        metadata_for(counts).prop_map(move |metadata| {
            MoleculeDsl::new(molecule.clone(), metadata).expect("generated metadata is coherent")
        })
    })
}

pub(crate) fn invalid_molecule_dsl_parts_strategy(
) -> impl Strategy<Value = (Molecule, MoleculeMetadata, Entity)> {
    molecule_with_constraints_strategy().prop_flat_map(|molecule| {
        let counts = ConstraintCounts::from_ir(&molecule);
        invalid_metadata_for(counts)
            .prop_map(move |(metadata, entity)| (molecule.clone(), metadata, entity))
    })
}

pub(crate) fn molecule_metadata_with_atom_subset_strategy(
) -> impl Strategy<Value = (Molecule, MoleculeMetadata, Vec<AtomId>)> {
    molecule_with_atom_subset_strategy().prop_flat_map(|(molecule, atoms)| {
        metadata_for(ConstraintCounts::from_ir(&molecule))
            .prop_map(move |metadata| (molecule.clone(), metadata, atoms.clone()))
    })
}

fn added_entities(reaction: &Reaction) -> Vec<Entity> {
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
        metadata_for(ConstraintCounts::from_ir(&reaction.lhs)).prop_map(move |lhs| {
            let mut metadata = ReactionMetadata::from(lhs);
            for entity in added_entities(&reaction) {
                metadata
                    .set_delta_keyword(entity, delta_keyword(entity))
                    .expect("generated delta keywords are disjoint");
            }
            metadata
                .add_atom_alias("reaction_alias", AtomForm::from_element(Element::F))
                .expect("generated reaction alias is disjoint and bijective");
            ReactionDsl::new(reaction.clone(), metadata)
                .expect("generated reaction metadata is coherent")
        })
    })
}

pub(crate) fn invalid_reaction_dsl_parts_strategy(
) -> impl Strategy<Value = (Reaction, ReactionMetadata, MetadataError)> {
    prop_oneof![
        comprehensive_reaction_strategy().prop_flat_map(|reaction| {
            invalid_metadata_for(ConstraintCounts::from_ir(&reaction.lhs)).prop_map(
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
) -> impl Strategy<Value = (ReactionSpan, MoleculeMetadata, Entity)> {
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

pub(crate) fn transaction_atoms(count: usize) -> Vec<AtomForm> {
    (0..count)
        .map(|id| {
            let element = ELEMENTS[id % ELEMENTS.len()];
            AtomForm::from_element(element)
        })
        .collect()
}

pub(crate) fn transaction_path_bonds(count: usize) -> Vec<AddBond> {
    (0..count.saturating_sub(1))
        .map(|id| AddBond {
            endpoints: [AtomHandle::New(id), AtomHandle::New(id + 1)],
            attributes: BondForm::from_order((id % 3 + 1) as u8),
        })
        .collect()
}

pub(crate) fn transaction_path_molecule(count: usize) -> Molecule {
    let atoms = transaction_atoms(count);
    let bonds = (0..count.saturating_sub(1))
        .map(|id| {
            (
                AtomId(id as u32),
                AtomId((id + 1) as u32),
                BondForm::from_order((id % 3 + 1) as u8),
            )
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
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
    pub(crate) fn base(&self) -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: INITIAL_HANDLE_ELEMENTS[..self.initial_count]
                .iter()
                .copied()
                .map(AtomForm::from_element)
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
                .map(AtomForm::from_element),
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
        let sentinel = edits.add_atom(AtomForm::from_element(SENTINEL_HANDLE_ELEMENT));
        edits.push(Edit::ModifyAtomField {
            id: if self.target_created {
                created[self.target_index].clone()
            } else {
                AtomHandle::Id(AtomId(self.target_index as u32))
            },
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(7),
            },
        });
        edits.push(Edit::ModifyAtomField {
            id: sentinel,
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(9),
            },
        });
        edits
    }

    pub(crate) fn expected(&self) -> Molecule {
        let initial = INITIAL_HANDLE_ELEMENTS[..self.initial_count]
            .iter()
            .copied()
            .enumerate()
            .filter(|(index, _)| !self.remove_initial[*index])
            .map(|(index, element)| {
                if !self.target_created && index == self.target_index {
                    AtomForm::from_element(element).with_charge(7_i64)
                } else {
                    AtomForm::from_element(element)
                }
            });
        let created = CREATED_HANDLE_ELEMENTS[..self.created_count]
            .iter()
            .copied()
            .enumerate()
            .filter(|(index, _)| !self.remove_created[*index])
            .map(|(index, element)| {
                if self.target_created && index == self.target_index {
                    AtomForm::from_element(element).with_charge(7_i64)
                } else {
                    AtomForm::from_element(element)
                }
            });
        Molecule::from_entries(MoleculeEntries {
            atoms: initial
                .chain(created)
                .chain([AtomForm::from_element(SENTINEL_HANDLE_ELEMENT).with_charge(9_i64)])
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
    pub(crate) fn base(&self) -> Molecule {
        let atoms = (0..self.count * 2)
            .map(|index| AtomForm::from_element(INITIAL_HANDLE_ELEMENTS[index % 4]))
            .collect();
        let bonds = (0..self.count)
            .map(|index| {
                (
                    AtomId((index * 2) as u32),
                    AtomId((index * 2 + 1) as u32),
                    BondForm::from_order(1),
                )
            })
            .collect();
        let dative = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32)],
                    AtomId((index * 2 + 1) as u32),
                    DativeBondForm::from_order(1),
                )
            })
            .collect();
        let aromatic = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32), AtomId((index * 2 + 1) as u32)],
                    AromaticSystemForm::default(),
                )
            })
            .collect();
        let multicenter = (0..self.count)
            .map(|index| {
                (
                    vec![AtomId((index * 2) as u32), AtomId((index * 2 + 1) as u32)],
                    MulticenterBondForm::default(),
                )
            })
            .collect();
        let noncovalent = (0..self.count)
            .map(|index| {
                (
                    AtomId((index * 2) as u32),
                    AtomId((index * 2 + 1) as u32),
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )
            })
            .collect();
        let stereo_atoms = (0..self.count)
            .map(|index| {
                let site = AtomId((index * 2) as u32);
                let other = AtomId((index * 2 + 1) as u32);
                (
                    site,
                    vec![
                        StereoLigand::new(other, StereoLigandKind::Atom),
                        StereoLigand::new(site, StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(site, StereoLigandKind::LonePair),
                        StereoLigand::new(other, StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )
            })
            .collect();
        let stereo_bonds = (0..self.count)
            .map(|index| {
                let first = AtomId((index * 2) as u32);
                let second = AtomId((index * 2 + 1) as u32);
                (
                    BondId(index as u32),
                    vec![
                        StereoLigand::new(first, StereoLigandKind::Atom),
                        StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(second, StereoLigandKind::Atom),
                        StereoLigand::new(second, StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )
            })
            .collect();
        Molecule::from_entries(MoleculeEntries {
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
                        attributes: BondForm::from_order(1),
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
                            DativeBondForm::from_order(1),
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
                            AromaticSystemForm::default(),
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
                            MulticenterBondForm::default(),
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
                            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: (0..self.count)
                    .map(|position| {
                        let site = AtomId((position * 2) as u32);
                        let other = AtomId((position * 2 + 1) as u32);
                        (
                            StereoAtomHandle::Id(StereoAtomId(
                                if position == self.invalid_position {
                                    invalid
                                } else {
                                    position as u32
                                },
                            )),
                            AtomHandle::Id(site),
                            vec![
                                (AtomHandle::Id(other), StereoLigandKind::Atom),
                                (AtomHandle::Id(site), StereoLigandKind::ImplicitHydrogen),
                                (AtomHandle::Id(site), StereoLigandKind::LonePair),
                                (AtomHandle::Id(other), StereoLigandKind::ImplicitHydrogen),
                            ],
                            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: (0..self.count)
                    .map(|position| {
                        let first = AtomId((position * 2) as u32);
                        let second = AtomId((position * 2 + 1) as u32);
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
                                (AtomHandle::Id(first), StereoLigandKind::Atom),
                                (AtomHandle::Id(first), StereoLigandKind::ImplicitHydrogen),
                                (AtomHandle::Id(second), StereoLigandKind::Atom),
                                (AtomHandle::Id(second), StereoLigandKind::ImplicitHydrogen),
                            ],
                            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
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
    pub(crate) fn base(&self) -> Molecule {
        match self {
            Self::AddPath { .. } => Molecule::default(),
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
                        old: NumForm::default(),
                        new: NumForm::Lit(*charge),
                    },
                }])
            }
            Self::SetBondOrder { count, id, order } => {
                let bond_id = id % (count - 1);
                Edits::from_iter([Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(bond_id as u32)),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit((bond_id % 3 + 1) as i64),
                        new: NumForm::Lit(*order as i64),
                    },
                }])
            }
            Self::AddAtomConstraint { count, id, size } => {
                Edits::from_iter([Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId((id % count) as u32)),
                    old: None,
                    new: Some(AtomConstraintForm::ring_membership(
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
                    DativeBondForm::from_order(1),
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

fn transaction_all_entities_molecule() -> Molecule {
    let bond_ligands = (0..4)
        .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let atom_ligands = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
    ];
    Molecule::from_entries(MoleculeEntries {
        atoms: (0..4).map(|_| AtomForm::from_element(Element::C)).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            AtomId(0),
            AtomId(3),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            atom_ligands,
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )],
        stereo_bonds: vec![(
            BondId(0),
            bond_ligands,
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    })
}

fn transaction_field_cases() -> Vec<(Molecule, Edits)> {
    let base = transaction_all_entities_molecule();
    let value = |change| (base.clone(), Edits::from_iter([change]));
    vec![
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Element {
                old: ElementForm::Lit(Element::C),
                new: ElementForm::Lit(Element::N),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::IsotopeMass {
                old: IsotopeMassForm::default(),
                new: IsotopeMassForm::Lit(13),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::ImplicitHydrogens {
                old: NumForm::default(),
                new: NumForm::Lit(3),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::LonePairs {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }),
        value(Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm::default(),
                new: UnpairedElectronsForm::from((2_u8, 1_u8)),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        }),
        value(Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm::default(),
                new: UnpairedElectronsForm::from((2_u8, 3_u8)),
            },
        }),
        value(Edit::ModifyDativeBondField {
            id: DativeBondHandle::Id(DativeBondId(0)),
            change: DativeBondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::default(),
                new: ElectronCountsForm::Lit(vec![1, 1, 1]),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }),
        value(Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change: AromaticSystemFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm::default(),
                new: UnpairedElectronsForm::from((1_u8, 2_u8)),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::Electrons {
                old: ElectronCountsForm::default(),
                new: ElectronCountsForm::Lit(vec![1, 1, 1]),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        }),
        value(Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm::default(),
                new: UnpairedElectronsForm::from((1_u8, 2_u8)),
            },
        }),
        value(Edit::ModifyNoncovalentBondField {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            },
        }),
        value(Edit::ModifyStereoAtomField {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
        }),
        value(Edit::ModifyStereoBondField {
            id: StereoBondHandle::Id(StereoBondId(0)),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        }),
    ]
}

fn transaction_constraint_cases() -> Vec<(Molecule, Edits)> {
    let base = transaction_all_entities_molecule();
    vec![
        Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintForm::degree(3)),
        },
        Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: None,
            new: Some(BondConstraintForm::aromatic(true)),
        },
        Edit::ModifyDativeBondConstraint {
            id: DativeBondHandle::Id(DativeBondId(0)),
            old: None,
            new: Some(DativeBondConstraintForm::aromatic(true)),
        },
        Edit::ModifyAromaticSystemConstraint {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            old: None,
            new: Some(AromaticSystemConstraintForm::electron_count(6)),
        },
        Edit::ModifyMulticenterBondConstraint {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            old: None,
            new: Some(MulticenterBondConstraintForm::electron_count(2)),
        },
        Edit::ModifyNoncovalentBondConstraint {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            old: None,
            new: Some(NoncovalentBondConstraintForm::intramolecular(true)),
        },
        Edit::ModifyStereoAtomConstraint {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoAtomConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            )),
        },
        Edit::ModifyStereoBondConstraint {
            id: StereoBondHandle::Id(StereoBondId(0)),
            kind: Some(StereoKind::CisTrans),
            old: None,
            new: Some(StereoBondConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            )),
        },
    ]
    .into_iter()
    .map(|edit| (base.clone(), Edits::from_iter([edit])))
    .collect()
}

fn transaction_removal_cases() -> Vec<(Molecule, Edits)> {
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
                DativeBondForm::from_order(1),
            )],
        },
        Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                atom_handles(&[0, 1, 2]),
                AromaticSystemForm::default(),
            )],
        },
        Edit::RemoveMulticenterBonds {
            removes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                atom_handles(&[0, 1, 2]),
                MulticenterBondForm::default(),
            )],
        },
        Edit::RemoveNoncovalentBonds {
            removes: vec![(
                NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(3))],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        },
        Edit::RemoveStereoAtoms {
            removes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(0)),
                vec![
                    (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(0)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
        },
    ]
    .into_iter()
    .map(|edit| (base.clone(), Edits::from_iter([edit])))
    .collect()
}

fn transaction_creation_case(include_created_constraint: bool) -> (Molecule, Edits) {
    let base = transaction_all_entities_molecule();
    let mut edits = Edits::new();
    let atom = edits.add_atom(AtomForm::from_element(Element::N));
    let bond = edits.add_bond(
        AtomHandle::Id(AtomId(1)),
        AtomHandle::Id(AtomId(2)),
        BondForm::from_order(2),
    );
    let dative = edits.add_dative_bond(
        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        DativeBondForm::from_order(1),
    );
    let aromatic = edits.add_aromatic_system(
        vec![AtomHandle::Id(AtomId(3)), atom.clone()],
        AromaticSystemForm::default(),
    );
    let multicenter = edits.add_multicenter_bond(
        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        MulticenterBondForm::default(),
    );
    let noncovalent = edits.add_noncovalent_bond(
        [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
    );
    let bond_ligands = (0..4)
        .map(|id| (AtomHandle::Id(AtomId(id)), StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let atom_ligands = vec![
        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
        (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
        (
            AtomHandle::Id(AtomId(1)),
            StereoLigandKind::ImplicitHydrogen,
        ),
    ];
    let stereo_atom = edits.add_stereo_atom(
        AtomHandle::Id(AtomId(1)),
        atom_ligands,
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
    );
    let stereo_bond = edits.add_stereo_bond(
        BondHandle::Id(BondId(1)),
        bond_ligands,
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
    );
    let source = Constraint::And(vec![
        Constraint::Atom(AtomId(7), AtomConstraintForm::degree(3)),
        Constraint::Bond(BondId(7), BondConstraintForm::aromatic(true)),
        Constraint::DativeBond(DativeBondId(7), DativeBondConstraintForm::aromatic(true)),
        Constraint::AromaticSystem(
            AromaticSystemId(7),
            AromaticSystemConstraintForm::electron_count(6),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(7),
            MulticenterBondConstraintForm::electron_count(2),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(7),
            NoncovalentBondConstraintForm::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(7),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(7),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
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

fn complete_transaction_cases(include_created_constraint: bool) -> Vec<(Molecule, Edits)> {
    let mut cases = transaction_field_cases();
    cases.extend(transaction_constraint_cases());
    cases.extend(transaction_removal_cases());
    cases.push(transaction_creation_case(include_created_constraint));
    let constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(3));
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

pub(crate) fn complete_transaction_strategy() -> impl Strategy<Value = (Molecule, Edits)> {
    prop::sample::select(complete_transaction_cases(true))
}

fn transaction_compaction_molecule(constraints: Constraints) -> Molecule {
    let atoms = (0..6).map(|_| AtomForm::from_element(Element::C)).collect();
    let bonds = (0..3)
        .map(|index| {
            (
                AtomId(index * 2),
                AtomId(index * 2 + 1),
                BondForm::from_order(1),
            )
        })
        .collect();
    let pairs = [[0_u32, 1_u32], [2, 3], [4, 5]];
    let dative = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a)], AtomId(*b), DativeBondForm::from_order(1)))
        .collect();
    let aromatic = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a), AtomId(*b)], AromaticSystemForm::default()))
        .collect();
    let multicenter = pairs
        .iter()
        .map(|[a, b]| (vec![AtomId(*a), AtomId(*b)], MulticenterBondForm::default()))
        .collect();
    let noncovalent = pairs
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                AtomId(*b),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )
        })
        .collect();
    let stereo_atoms = pairs
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                vec![
                    StereoLigand::new(AtomId(*b), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(*a), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(*a), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(*b), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
                    StereoLigand::new(AtomId(*a), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(*b), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(*b), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
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
        EntityKind::Atom => Constraint::Atom(AtomId(id), AtomConstraintForm::degree(value)),
        EntityKind::Bond => Constraint::Bond(
            BondId(id),
            BondConstraintForm::ring_membership(RingScope::All, value),
        ),
        EntityKind::DativeBond => Constraint::DativeBond(
            DativeBondId(id),
            DativeBondConstraintForm::aromatic(value % 2 == 0),
        ),
        EntityKind::AromaticSystem => Constraint::AromaticSystem(
            AromaticSystemId(id),
            AromaticSystemConstraintForm::electron_count(value),
        ),
        EntityKind::MulticenterBond => Constraint::MulticenterBond(
            MulticenterBondId(id),
            MulticenterBondConstraintForm::electron_count(value),
        ),
        EntityKind::NoncovalentBond => Constraint::NoncovalentBond(
            NoncovalentBondId(id),
            NoncovalentBondConstraintForm::intramolecular(value % 2 == 0),
        ),
        EntityKind::StereoAtom => Constraint::StereoAtom(
            StereoAtomId(id),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(match value % 3 {
                0 => Stereogenicity::Symmetric,
                1 => Stereogenicity::Prochiral,
                _ => Stereogenicity::Stereogenic,
            })),
        ),
        EntityKind::StereoBond => Constraint::StereoBond(
            StereoBondId(id),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(match value % 3 {
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
    base: Molecule,
    expected: Vec<Constraint>,
}

impl ConstraintCompactionCase {
    pub(crate) fn base(&self) -> Molecule {
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
                    DativeBondForm::from_order(1),
                )],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                )],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                )],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
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

pub(crate) fn consecutive_transaction_strategy() -> impl Strategy<Value = (Molecule, Edits, Edits)>
{
    (-4_i64..=4, -4_i64..=4)
        .prop_filter("successive charges must differ", |(first, second)| {
            first != second
        })
        .prop_map(|(first, second)| {
            let base = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                ..Default::default()
            });
            let first_edits = Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(first),
                },
            }]);
            let second_edits = Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(first),
                    new: NumForm::Lit(second),
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

pub(crate) fn overlay_transaction_base() -> Molecule {
    let atoms: Vec<AtomForm> = (0..6).map(|_| AtomForm::from_element(Element::C)).collect();
    let bonds = (0..5)
        .map(|i| (AtomId(i), AtomId(i + 1), BondForm::from_order(1)))
        .collect();
    let dative = (0..2)
        .map(|i| {
            (
                DATIVE_DONORS[i]
                    .iter()
                    .map(|&a| AtomId(a))
                    .collect::<Vec<_>>(),
                AtomId(DATIVE_ACCEPTORS[i]),
                DativeBondForm::from_order(1),
            )
        })
        .collect();
    let aromatic = AROMATIC_SETS
        .iter()
        .map(|set| {
            (
                set.iter().map(|&a| AtomId(a)).collect::<Vec<_>>(),
                AromaticSystemForm::default(),
            )
        })
        .collect();
    let multicenter = MULTICENTER_SETS
        .iter()
        .map(|set| {
            (
                set.iter().map(|&a| AtomId(a)).collect::<Vec<_>>(),
                MulticenterBondForm::default(),
            )
        })
        .collect();
    let noncovalent = NONCOVALENT_PAIRS
        .iter()
        .map(|[a, b]| {
            (
                AtomId(*a),
                AtomId(*b),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
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
pub(crate) fn overlay_transaction_strategy() -> impl Strategy<Value = (Molecule, Edits)> {
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
                            .map(|_| AtomForm::from_element(Element::C))
                            .collect(),
                    });
                }
                // Base carbons carry the default (`Undetermined`) charge, so that is the `old` value.
                for i in (0..6).filter(|&i| mod_at[i]) {
                    edits.push(Edit::ModifyAtomField {
                        id: AtomHandle::Id(AtomId(i as u32)),
                        change: AtomFieldChange::Charge {
                            old: NumForm::default(),
                            new: NumForm::Lit(1),
                        },
                    });
                }
                // Molecule-level constraints referencing overlays, added before the removals so that
                // removing a referenced overlay exercises constraint drop/remap + its rollback restore.
                for i in (0..2).filter(|&i| con_ar[i]) {
                    edits.push(Edit::AddMoleculeConstraint {
                        constraint: Constraint::AromaticSystem(
                            AromaticSystemId(i as u32),
                            AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6)),
                        )
                        .into(),
                    });
                }
                for i in (0..2).filter(|&i| con_mc[i]) {
                    edits.push(Edit::AddMoleculeConstraint {
                        constraint: Constraint::MulticenterBond(
                            MulticenterBondId(i as u32),
                            MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(4)),
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
                            DativeBondForm::from_order(1),
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
                            AromaticSystemForm::default(),
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
                            MulticenterBondForm::default(),
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
                            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
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
pub(crate) fn transaction_edits_strategy() -> impl Strategy<Value = (Molecule, Edits)> {
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
fn simple_molecule_strategy() -> impl Strategy<Value = Molecule> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                prop::collection::vec(
                    element_strategy().prop_map(AtomForm::from_element),
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
                .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondForm::from_order(order)))
                .collect();
            Molecule::from_entries(MoleculeEntries {
                atoms,
                bonds,
                ..Default::default()
            })
        })
}

pub(crate) fn reaction_strategy() -> impl Strategy<Value = Reaction> {
    reaction_over(simple_molecule_strategy())
}

pub(crate) fn replacement_reaction_strategy() -> impl Strategy<Value = Reaction> {
    (molecule_strategy(), molecule_strategy()).prop_map(|(lhs, rhs)| {
        let correspondence =
            Correspondence::new(Vec::new(), lhs.atoms().count(), rhs.atoms().count())
                .expect("correspondence producer preserves partial-bijection invariants");
        Reaction::from_sides(lhs, rhs, correspondence)
            .expect("an empty atom correspondence uniquely relates replacement sides")
    })
}

pub(crate) fn comprehensive_reaction_strategy() -> BoxedStrategy<Reaction> {
    prop_oneof![
        2 => overlay_reaction_strategy(),
        1 => replacement_reaction_strategy(),
    ]
    .boxed()
}

/// A localized molecule with DAMN overlays (dative / aromatic / multicenter / noncovalent) plus
/// stereo (tetrahedral atoms / cis-trans bonds) and no molecule constraints (orthogonal). 1–4 atoms;
/// overlays generated as in `molecule_strategy`, scoped.
fn overlay_molecule_strategy() -> impl Strategy<Value = Molecule> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                Just(atom_count),
                prop::collection::vec(
                    element_strategy().prop_map(AtomForm::from_element),
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
                        (Just(atoms), aromatic_system_form_for(n))
                    },
                ),
                0..=1,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_form_for(n))
                    },
                ),
                0..=1,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_form_strategy(),
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
                    .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondForm::from_order(order)))
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
                Molecule::from_entries(MoleculeEntries {
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

/// Cosets valid for `kind`: `Undetermined` or an in-range `Lit` index (`0..kind.count()`). Unlike
/// the generic `stereo_coset_strategy`, indices are bounded by the kind's coset count.
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
            Just(ElectronCountsForm::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsForm::Lit),
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
            Just(ElectronCountsForm::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsForm::Lit),
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

/// 0..=1 tetrahedral stereo atoms over an `atom_count`-atom molecule. The frame contains distinct
/// actual atoms other than the site plus at most one implicit hydrogen and lone pair anchored at the
/// site. Ligand atoms need not be graph neighbors (tier-1 requires only the ligand count for the
/// kind).
fn stereo_atom_overlay_strategy(
    atom_count: usize,
) -> BoxedStrategy<Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>> {
    let degree = StereoKind::Tetrahedral.degree();
    if atom_count + 1 < degree {
        return Just(Vec::new()).boxed();
    }
    let entry = (0..atom_count as u32).prop_flat_map(move |site| {
        let site_id = AtomId(site);
        let pool: Vec<StereoLigand> = (0..atom_count as u32)
            .filter(|atom| *atom != site)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .chain([
                StereoLigand::new(site_id, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(site_id, StereoLigandKind::LonePair),
            ])
            .collect();
        (
            Just(site_id),
            Just(pool).prop_shuffle().prop_map(move |mut pool| {
                pool.truncate(degree);
                pool
            }),
            stereo_coset_for_kind(StereoKind::Tetrahedral),
        )
    });
    prop::collection::vec(entry, 0..=1)
        .prop_map(move |entries| {
            entries
                .into_iter()
                .map(|(site, ligands, coset)| {
                    let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset);
                    (site, ligands, attributes)
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
) -> BoxedStrategy<Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>> {
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
                let attributes = StereoBondForm::new(StereoKind::CisTrans, coset);
                (BondId(site), ligands, attributes)
            })
            .collect()
    })
    .boxed()
}

/// A reaction whose `lhs` carries DAMN overlays — exercises overlay carry, correspondence, and
/// co-deletion through compose.
pub(crate) fn overlay_reaction_strategy() -> impl Strategy<Value = Reaction> {
    reaction_over(overlay_molecule_strategy())
}

/// An optional source operation used to derive an absolute stereo configuration delta.
#[derive(Clone, Debug)]
enum StereoOp {
    Swap,
    Mirror,
    Apply(Permutation),
    SetCoset(StereoCoset),
}

/// Generate optional source operations over valid configurations. These operations are evaluated
/// while constructing the reaction and emitted as absolute before/after deltas.
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
/// incident bonds, overlays, and stereo entities), per-surviving-entity optional field
/// edits (the absolute `old` read from `lhs`, so apply's precondition holds), plus up to two new
/// atoms bonded to the lowest survivor. No dangling by construction.
fn reaction_over(molecule: impl Strategy<Value = Molecule>) -> impl Strategy<Value = Reaction> {
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
    lhs: Molecule,
    removals: Vec<bool>,
    charges: Vec<Option<i64>>,
    orders: Vec<Option<i64>>,
    additions: Vec<Element>,
    overlay_ops: OverlayOps,
    stereo_ops: StereoOps,
) -> Reaction {
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
            attributes: lhs.atom(id).attributes.clone(),
        }));
    }
    for &id in &removed_bonds {
        let [x, y] = lhs.raw_graph().edge_endpoints(EdgeId(id.0));
        deltas.push(Delta::Bond(BondDelta::Remove {
            id,
            atoms: [AtomId::from(x), AtomId::from(y)],
            attributes: lhs.bond(id).attributes.clone(),
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
            attributes: view.attributes.clone(),
        }));
    }
    for &id in &removed_aromatic {
        let view = lhs.aromatic_system(id);
        deltas.push(Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            attributes: view.attributes.clone(),
        }));
    }
    for &id in &removed_multicenter {
        let view = lhs.multicenter_bond(id);
        deltas.push(Delta::MulticenterBond(MulticenterBondDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            attributes: view.attributes.clone(),
        }));
    }
    for &id in &removed_noncovalent {
        let view = lhs.noncovalent_bond(id);
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
            id,
            atoms: view.atom_ids(),
            attributes: view.attributes.clone(),
        }));
    }
    // A removed atom also takes its incident stereo entities (site OR ligand incidence), else
    // apply and span materialization would otherwise leave dangling references. `incident_ids`
    // covers both site and ligand incidence.
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
            attributes: view.attributes.clone(),
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
            attributes: view.attributes.clone(),
        }));
    }
    for (index, new_charge) in charges.into_iter().enumerate() {
        let id = AtomId(index as u32);
        if removed_atoms.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.atom(id).attributes.charge.clone();
        let new = NumForm::Lit(charge);
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
        let old = lhs.bond(id).attributes.order.clone();
        let new = NumForm::Lit(order);
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
        let old = lhs.dative_bond(id).attributes.order.clone();
        let new = NumForm::Lit(order);
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
        let old = lhs.aromatic_system(id).attributes.charge.clone();
        let new = NumForm::Lit(charge);
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
        let old = lhs.multicenter_bond(id).attributes.charge.clone();
        let new = NumForm::Lit(charge);
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
            .attributes
            .constraints
            .iter()
            .any(|c| matches!(c, DativeBondConstraintForm::Aromatic(_)));
        if has_aromatic {
            continue;
        }
        deltas.push(Delta::DativeBond(DativeBondDelta::ModifyConstraint {
            id,
            old: None,
            new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        }));
    }
    // Part B — stereo edits on survivors. Source operations are evaluated against the lhs and
    // emitted as absolute before/after deltas. Value no-ops are omitted because they would
    // materialize a spurious `Modified { X, X }` span state.
    let (stereo_atom_ops, stereo_bond_ops) = stereo_ops;
    for (index, op) in stereo_atom_ops.into_iter().enumerate() {
        let id = StereoAtomId(index as u32);
        if removed_stereo_atom.contains(&id) {
            continue;
        }
        let Some(op) = op else { continue };
        let kind = lhs.stereo_atom(id).kind();
        let old = lhs.stereo_atom(id).attributes.configuration.clone();
        let new = match &op {
            StereoOp::Swap => old.swap(),
            StereoOp::Mirror => old.mirror(),
            StereoOp::Apply(permutation) => old.apply(*permutation),
            StereoOp::SetCoset(coset) => StereoConfigurationForm::kinded(kind, coset.clone()),
        };
        let delta = StereoAtomDelta::ModifyField {
            id,
            change: StereoAtomFieldChange::Configuration {
                old: old.clone(),
                new: new.clone(),
            },
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
        let old = lhs.stereo_bond(id).attributes.configuration.clone();
        let new = match &op {
            StereoOp::Swap => old.swap(),
            StereoOp::Mirror => old.mirror(),
            StereoOp::Apply(permutation) => old.apply(*permutation),
            StereoOp::SetCoset(coset) => StereoConfigurationForm::kinded(kind, coset.clone()),
        };
        let delta = StereoBondDelta::ModifyField {
            id,
            change: StereoBondFieldChange::Configuration {
                old: old.clone(),
                new: new.clone(),
            },
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
            attributes: AtomForm::from_element(element),
        }));
        if let Some(anchor) = anchor {
            deltas.push(Delta::Bond(BondDelta::Add {
                id: BondId((bond_count + offset) as u32),
                atoms: [anchor, atom],
                attributes: BondForm::from_order(1),
            }));
        }
    }
    // Part A — overlay `Add`: a noncovalent bond between the two newly-added atoms (both created in
    // this reaction, so no dangling). Ids append past the lhs noncovalent count.
    if add_noncovalent && added_atom_ids.len() >= 2 {
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(lhs.noncovalent_bonds().count() as u32),
            atoms: [added_atom_ids[0], added_atom_ids[1]],
            attributes: NoncovalentBondForm {
                kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::VanDerWaals),
                constraints: Default::default(),
            },
        }));
    }
    Reaction::new(lhs, deltas)
}
