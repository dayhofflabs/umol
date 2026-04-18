//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::molecule::MoleculeAst;
use crate::ast::{AromaticSystemIdx, AtomIdx, BondIdx, MulticenterBondIdx};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    AtomDerived(AtomIdx, AtomConstraint),
    BondDerived(BondIdx, BondConstraint),
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
            Self::AtomDerived(_, c) => c.is_ground(),
            Self::BondDerived(_, c) => c.is_ground(),
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

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AtomConstraint {
    ValenceSum(ValueAst),
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
            Self::ValenceSum(v)
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
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum BondConstraint {
    RingBond,
}

impl BondConstraint {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::RingBond => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::AtomIdx;

    #[rstest]
    #[case::and_pair(
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomDerived(AtomIdx(1), AtomConstraint::InRing),
        ]),
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomDerived(AtomIdx(1), AtomConstraint::InRing),
        ]),
        true,
    )]
    #[case::and_order_matters(
        MoleculeConstraint::And(vec![
            MoleculeConstraint::TotalCharge(ValueAst::Lit(0)),
            MoleculeConstraint::AtomDerived(AtomIdx(1), AtomConstraint::InRing),
        ]),
        MoleculeConstraint::And(vec![
            MoleculeConstraint::AtomDerived(AtomIdx(1), AtomConstraint::InRing),
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
            MoleculeConstraint::Not(Box::new(MoleculeConstraint::AtomDerived(
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
        MoleculeConstraint::AtomDerived(AtomIdx(0), AtomConstraint::ValenceSum(ValueAst::Lit(4))),
        true,
    )]
    #[case::atom_derived_undetermined(
        MoleculeConstraint::AtomDerived(AtomIdx(0), AtomConstraint::ValenceSum(ValueAst::Undetermined)),
        false,
    )]
    #[case::bond_derived_ring(
        MoleculeConstraint::BondDerived(BondIdx(0), BondConstraint::RingBond),
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
    #[case::valence_sum_lit(AtomConstraint::ValenceSum(ValueAst::Lit(4)), true)]
    #[case::valence_sum_undetermined(AtomConstraint::ValenceSum(ValueAst::Undetermined), false)]
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
    #[case::valence_sum_eq(
        AtomConstraint::ValenceSum(ValueAst::Lit(4)),
        AtomConstraint::ValenceSum(ValueAst::Lit(4)),
        true,
    )]
    #[case::valence_sum_payload_diff(
        AtomConstraint::ValenceSum(ValueAst::Lit(4)),
        AtomConstraint::ValenceSum(ValueAst::Lit(3)),
        false,
    )]
    #[case::variant_diff(
        AtomConstraint::ValenceSum(ValueAst::Lit(4)),
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
