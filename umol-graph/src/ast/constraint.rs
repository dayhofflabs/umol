//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use std::mem;

use indexmap::IndexMap;
use strum::EnumDiscriminants;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use super::molecule::MoleculeAst;
use super::{AromaticSystemIdx, AtomIdx, BondIdx, MulticenterBondIdx};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    AtomPred(AtomIdx, AtomConstraint),
    BondPred(BondIdx, BondConstraint),
    TotalCharge(ValueAst),
    TotalSpin(SpinStateAst),
    AromaticElectronCount(AromaticSystemIdx, ValueAst),
    MulticenterElectronCount(MulticenterBondIdx, ValueAst),
    BondOrderSum(Vec<BondIdx>, ValueAst),
    Connected(Vec<AtomIdx>),
    SubPattern {
        anchor: AtomIdx,
        pattern: Box<MoleculeAst>,
    },
    And(Vec<MoleculeConstraint>),
    Or(Vec<MoleculeConstraint>),
    Not(Box<MoleculeConstraint>),
}

impl MoleculeConstraint {
    /// A ground assertion carries only literal values (no wildcards, variables,
    /// or expressions) and is not a query combinator. These are facts about a
    /// resolved molecule.
    pub fn is_ground_assertion(&self) -> bool {
        match self {
            Self::AtomPred(_, c) => c.is_ground(),
            Self::BondPred(_, c) => c.is_ground(),
            Self::TotalCharge(v)
            | Self::AromaticElectronCount(_, v)
            | Self::MulticenterElectronCount(_, v)
            | Self::BondOrderSum(_, v) => v.is_ground(),
            Self::TotalSpin(s) => s.is_ground(),
            Self::Connected(_) => true,
            Self::SubPattern { .. } | Self::And(_) | Self::Or(_) | Self::Not(_) => false,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticValenceConstraint {
    NotAromatic,
    Value(ValueAst),
}

impl AromaticValenceConstraint {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::NotAromatic => true,
            Self::Value(v) => v.is_ground(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash))]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceConstraint),
    MulticenterValence(ValueAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    Connectivity(ValueAst),
    TotalHCount(ValueAst),
    InRing,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl AtomConstraint {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::MulticenterValence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::Connectivity(v)
            | Self::TotalHCount(v)
            | Self::RingCount(v)
            | Self::RingSize(v) => v.is_ground(),
            Self::AromaticValence(c) => c.is_ground(),
            Self::InRing => true,
        }
    }

    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash))]
pub enum BondConstraint {
    RingBond,
    Aromatic,
}

impl BondConstraint {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::RingBond | Self::Aromatic => true,
        }
    }

    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }
}

/// Set of `AtomConstraint`s attached to a single atom.
///
/// Uniqueness invariant: at most one constraint per `AtomConstraintKind`. Inserting
/// a constraint of a kind already present replaces the prior entry and returns it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomConstraintSet {
    items: Vec<AtomConstraint>,
}

impl AtomConstraintSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get(&self, kind: AtomConstraintKind) -> Option<&AtomConstraint> {
        self.items.iter().find(|c| c.kind() == kind)
    }

    pub fn insert(&mut self, c: AtomConstraint) -> Option<AtomConstraint> {
        let kind = c.kind();
        match self.items.iter().position(|e| e.kind() == kind) {
            Some(pos) => Some(mem::replace(&mut self.items[pos], c)),
            None => {
                self.items.push(c);
                None
            }
        }
    }

    pub fn remove(&mut self, kind: AtomConstraintKind) -> Option<AtomConstraint> {
        let pos = self.items.iter().position(|c| c.kind() == kind)?;
        Some(self.items.remove(pos))
    }

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.items.iter()
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AtomConstraint) -> bool) {
        self.items.retain(|c| f(c));
    }
}

impl FromIterator<AtomConstraint> for AtomConstraintSet {
    fn from_iter<I: IntoIterator<Item = AtomConstraint>>(iter: I) -> Self {
        let mut set = Self::new();
        for c in iter {
            set.insert(c);
        }
        set
    }
}

/// Set of `BondConstraint`s attached to a single bond. Uniqueness is per `BondConstraintKind`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondConstraintSet {
    items: Vec<BondConstraint>,
}

impl BondConstraintSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn get(&self, kind: BondConstraintKind) -> Option<&BondConstraint> {
        self.items.iter().find(|c| c.kind() == kind)
    }

    pub fn insert(&mut self, c: BondConstraint) -> Option<BondConstraint> {
        let kind = c.kind();
        match self.items.iter().position(|e| e.kind() == kind) {
            Some(pos) => Some(mem::replace(&mut self.items[pos], c)),
            None => {
                self.items.push(c);
                None
            }
        }
    }

    pub fn remove(&mut self, kind: BondConstraintKind) -> Option<BondConstraint> {
        let pos = self.items.iter().position(|c| c.kind() == kind)?;
        Some(self.items.remove(pos))
    }

    pub fn iter(&self) -> impl Iterator<Item = &BondConstraint> {
        self.items.iter()
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn retain(&mut self, mut f: impl FnMut(&BondConstraint) -> bool) {
        self.items.retain(|c| f(c));
    }
}

impl FromIterator<BondConstraint> for BondConstraintSet {
    fn from_iter<I: IntoIterator<Item = BondConstraint>>(iter: I) -> Self {
        let mut set = Self::new();
        for c in iter {
            set.insert(c);
        }
        set
    }
}

/// Partitioned constraint storage for `MoleculeAst`.
///
/// Per-target constraints are keyed by their participant index. Combinators,
/// global predicates, and structural queries live in `global`. Insertion
/// dispatches on the `MoleculeConstraint` variant; conflicting inserts on a
/// per-target slot replace the prior entry (the displaced value is dropped).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeConstraints {
    atoms: IndexMap<AtomIdx, AtomConstraintSet>,
    bonds: IndexMap<BondIdx, BondConstraintSet>,
    aromatic_systems: IndexMap<AromaticSystemIdx, ValueAst>,
    multicenter_bonds: IndexMap<MulticenterBondIdx, ValueAst>,
    global: Vec<MoleculeConstraint>,
}

impl MoleculeConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(&mut self, c: MoleculeConstraint) {
        match c {
            MoleculeConstraint::AtomPred(idx, inner) => {
                self.atoms.entry(idx).or_default().insert(inner);
            }
            MoleculeConstraint::BondPred(idx, inner) => {
                self.bonds.entry(idx).or_default().insert(inner);
            }
            MoleculeConstraint::AromaticElectronCount(idx, v) => {
                self.aromatic_systems.insert(idx, v);
            }
            MoleculeConstraint::MulticenterElectronCount(idx, v) => {
                self.multicenter_bonds.insert(idx, v);
            }
            other => self.global.push(other),
        }
    }

    pub fn atoms(&self) -> &IndexMap<AtomIdx, AtomConstraintSet> {
        &self.atoms
    }

    pub fn atoms_mut(&mut self) -> &mut IndexMap<AtomIdx, AtomConstraintSet> {
        &mut self.atoms
    }

    pub fn bonds(&self) -> &IndexMap<BondIdx, BondConstraintSet> {
        &self.bonds
    }

    pub fn bonds_mut(&mut self) -> &mut IndexMap<BondIdx, BondConstraintSet> {
        &mut self.bonds
    }

    pub fn aromatic_systems(&self) -> &IndexMap<AromaticSystemIdx, ValueAst> {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &IndexMap<MulticenterBondIdx, ValueAst> {
        &self.multicenter_bonds
    }

    pub fn global(&self) -> &[MoleculeConstraint] {
        &self.global
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
            && self.bonds.is_empty()
            && self.aromatic_systems.is_empty()
            && self.multicenter_bonds.is_empty()
            && self.global.is_empty()
    }

    pub fn len(&self) -> usize {
        self.atoms.values().map(|s| s.len()).sum::<usize>()
            + self.bonds.values().map(|s| s.len()).sum::<usize>()
            + self.aromatic_systems.len()
            + self.multicenter_bonds.len()
            + self.global.len()
    }

    /// Reconstitute a flat iterator of `MoleculeConstraint`s from the partitioned
    /// storage. Order: atom predicates, bond predicates, aromatic-system
    /// predicates, multicenter-bond predicates, global constraints.
    pub fn iter(&self) -> impl Iterator<Item = MoleculeConstraint> + '_ {
        let atom_iter = self.atoms.iter().flat_map(|(idx, set)| {
            set.iter()
                .map(move |c| MoleculeConstraint::AtomPred(*idx, c.clone()))
        });
        let bond_iter = self.bonds.iter().flat_map(|(idx, set)| {
            set.iter()
                .map(move |c| MoleculeConstraint::BondPred(*idx, c.clone()))
        });
        let aromatic_iter = self
            .aromatic_systems
            .iter()
            .map(|(idx, v)| MoleculeConstraint::AromaticElectronCount(*idx, v.clone()));
        let multicenter_iter = self
            .multicenter_bonds
            .iter()
            .map(|(idx, v)| MoleculeConstraint::MulticenterElectronCount(*idx, v.clone()));
        let global_iter = self.global.iter().cloned();
        atom_iter
            .chain(bond_iter)
            .chain(aromatic_iter)
            .chain(multicenter_iter)
            .chain(global_iter)
    }
}

impl FromIterator<MoleculeConstraint> for MoleculeConstraints {
    fn from_iter<I: IntoIterator<Item = MoleculeConstraint>>(iter: I) -> Self {
        let mut s = Self::new();
        for c in iter {
            s.insert(c);
        }
        s
    }
}

impl From<Vec<MoleculeConstraint>> for MoleculeConstraints {
    fn from(v: Vec<MoleculeConstraint>) -> Self {
        v.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use super::AtomIdx;

    #[rstest]
    #[case::and_pair(
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomPred(AtomIdx(1), AtomConstraint::InRing),
        ]),
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomPred(AtomIdx(1), AtomConstraint::InRing),
        ]),
        true,
    )]
    #[case::and_order_matters(
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomPred(AtomIdx(1), AtomConstraint::InRing),
        ]),
        MoleculeConstraint::And(vec![
            MoleculeConstraint::AtomPred(AtomIdx(1), AtomConstraint::InRing),
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
        ]),
        false,
    )]
    #[case::or_distinct_payload(
        MoleculeConstraint::Or(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::TotalCharge(ValueAst::Lit(1)),
        ]),
        MoleculeConstraint::Or(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::TotalCharge(ValueAst::Lit(2)),
        ]),
        false,
    )]
    #[case::not_idempotent_eq(
        MoleculeConstraint::Not(Box::new(MoleculeConstraint::TotalCharge(ValueAst::Lit(-1)))),
        MoleculeConstraint::Not(Box::new(MoleculeConstraint::TotalCharge(ValueAst::Lit(-1)))),
        true,
    )]
    #[case::and_or_distinct(
        MoleculeConstraint::And(vec![MoleculeConstraint::TotalCharge(ValueAst::Lit(0))]),
        MoleculeConstraint::Or(vec![MoleculeConstraint::TotalCharge(ValueAst::Lit(0))]),
        false,
    )]
    fn test_molecule_constraint_combinators_eq(
        #[case] left: MoleculeConstraint,
        #[case] right: MoleculeConstraint,
        #[case] equal: bool,
    ) {
        assert_eq!(left == right, equal);
    }

    #[test]
    fn test_molecule_constraint_combinators_nested() {
        let inner = MoleculeConstraint::Or(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(-1)),
            MoleculeConstraint::Not(Box::new(MoleculeConstraint::AtomPred(
                AtomIdx(0),
                AtomConstraint::InRing,
            ))),
        ]);
        let outer = MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            inner.clone(),
        ]);

        let MoleculeConstraint::And(children) = &outer else {
            panic!("expected And");
        };
        assert_eq!(children.len(), 2);
        assert_eq!(children[1], inner);
    }

    #[rstest]
    #[case::total_charge_lit(MoleculeConstraint::TotalCharge(ValueAst::Lit(0)), true)]
    #[case::total_charge_undetermined(MoleculeConstraint::TotalCharge(ValueAst::Undetermined), false)]
    #[case::atom_derived_ground(
        MoleculeConstraint::AtomPred(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
        true,
    )]
    #[case::atom_derived_undetermined(
        MoleculeConstraint::AtomPred(AtomIdx(0), AtomConstraint::Valence(ValueAst::Undetermined)),
        false,
    )]
    #[case::bond_derived_ring(
        MoleculeConstraint::BondPred(BondIdx(0), BondConstraint::RingBond),
        true,
    )]
    #[case::connected(MoleculeConstraint::Connected(vec![AtomIdx(0), AtomIdx(1)]), true)]
    #[case::sub_pattern(
        MoleculeConstraint::SubPattern {
            anchor: AtomIdx(0),
            pattern: Box::new(MoleculeAst::default()),
        },
        false,
    )]
    #[case::and_combinator(
        MoleculeConstraint::And(vec![MoleculeConstraint::TotalCharge(ValueAst::Lit(0))]),
        false,
    )]
    fn test_molecule_constraint_is_ground_assertion(
        #[case] constraint: MoleculeConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.is_ground_assertion(), expected);
    }

    #[rstest]
    #[case::valence_lit(AtomConstraint::Valence(ValueAst::Lit(4)), true)]
    #[case::valence_undetermined(AtomConstraint::Valence(ValueAst::Undetermined), false)]
    #[case::aromatic_valence_lit(
        AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(3))),
        true,
    )]
    #[case::aromatic_valence_set(
        AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::LitSet(vec![2, 3]))),
        false,
    )]
    #[case::aromatic_not_aromatic(
        AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic),
        true,
    )]
    #[case::multicenter_valence_lit(AtomConstraint::MulticenterValence(ValueAst::Lit(1)), true)]
    #[case::donated_pairs_lit(AtomConstraint::DonatedPairs(ValueAst::Lit(0)), true)]
    #[case::accepted_pairs_undetermined(AtomConstraint::AcceptedPairs(ValueAst::Undetermined), false)]
    #[case::degree_lit(AtomConstraint::Degree(ValueAst::Lit(3)), true)]
    #[case::connectivity_lit(AtomConstraint::Connectivity(ValueAst::Lit(4)), true)]
    #[case::total_h_count_lit(AtomConstraint::TotalHCount(ValueAst::Lit(2)), true)]
    #[case::in_ring(AtomConstraint::InRing, true)]
    #[case::ring_count_lit(AtomConstraint::RingCount(ValueAst::Lit(1)), true)]
    #[case::ring_size_lit(AtomConstraint::RingSize(ValueAst::Lit(6)), true)]
    #[case::ring_size_undetermined(AtomConstraint::RingSize(ValueAst::Undetermined), false)]
    fn test_atom_constraint_is_ground(#[case] constraint: AtomConstraint, #[case] expected: bool) {
        assert_eq!(constraint.is_ground(), expected);
    }

    #[rstest]
    #[case::valence_eq(
        AtomConstraint::Valence(ValueAst::Lit(4)),
        AtomConstraint::Valence(ValueAst::Lit(4)),
        true,
    )]
    #[case::valence_payload_diff(
        AtomConstraint::Valence(ValueAst::Lit(4)),
        AtomConstraint::Valence(ValueAst::Lit(3)),
        false,
    )]
    #[case::variant_diff(
        AtomConstraint::Valence(ValueAst::Lit(4)),
        AtomConstraint::MulticenterValence(ValueAst::Lit(4)),
        false,
    )]
    #[case::in_ring_eq(AtomConstraint::InRing, AtomConstraint::InRing, true)]
    #[case::in_ring_vs_ring_count(
        AtomConstraint::InRing,
        AtomConstraint::RingCount(ValueAst::Lit(1)),
        false,
    )]
    fn test_atom_constraint_eq(
        #[case] left: AtomConstraint,
        #[case] right: AtomConstraint,
        #[case] equal: bool,
    ) {
        assert_eq!(left == right, equal);
    }

    #[test]
    fn test_atom_constraint_clone() {
        let original = AtomConstraint::RingSize(ValueAst::Lit(6));
        let cloned = original.clone();
        assert_eq!(original, cloned);
    }

    #[rstest]
    #[case::ring_bond(BondConstraint::RingBond, true)]
    #[case::aromatic(BondConstraint::Aromatic, true)]
    fn test_bond_constraint_is_ground(#[case] constraint: BondConstraint, #[case] expected: bool) {
        assert_eq!(constraint.is_ground(), expected);
    }

    #[test]
    fn test_bond_constraint_clone() {
        let original = BondConstraint::RingBond;
        let cloned = original.clone();
        assert_eq!(original, cloned);
    }
}
