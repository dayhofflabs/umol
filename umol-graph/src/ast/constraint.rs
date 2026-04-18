//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::AtomIdx;
use crate::ast::molecule::MoleculeAst;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    SubPattern {
        anchor: AtomIdx,
        pattern: Box<MoleculeAst>,
    },
    Derived {
        predicate: DerivedPred,
        refs: RelationRefs,
    },
    Matcher(MatcherFlag),
    And(Vec<MoleculeConstraint>),
    Or(Vec<MoleculeConstraint>),
    Not(Box<MoleculeConstraint>),
}

impl MoleculeConstraint {
    /// A ground assertion is a `Derived` constraint whose predicate carries only
    /// literal values (no wildcards, variables, or expressions). These are facts
    /// about a resolved molecule, not queries.
    pub fn is_ground_assertion(&self) -> bool {
        match self {
            Self::Derived { predicate, .. } => predicate.is_ground(),
            _ => false,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RelationRefs {
    pub atoms: Vec<AtomIdx>,
    pub bonds: Vec<usize>,
    pub dative_bonds: Vec<usize>,
    pub aromatic_systems: Vec<usize>,
    pub multicenter_bonds: Vec<usize>,
    pub noncovalent_bonds: Vec<usize>,
}

impl RelationRefs {
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
            && self.bonds.is_empty()
            && self.dative_bonds.is_empty()
            && self.aromatic_systems.is_empty()
            && self.multicenter_bonds.is_empty()
            && self.noncovalent_bonds.is_empty()
    }

    pub fn atoms(atoms: Vec<AtomIdx>) -> Self {
        Self { atoms, ..Self::default() }
    }

    pub fn bonds(bonds: Vec<usize>) -> Self {
        Self { bonds, ..Self::default() }
    }

    pub fn dative_bonds(dative_bonds: Vec<usize>) -> Self {
        Self { dative_bonds, ..Self::default() }
    }

    pub fn aromatic_systems(aromatic_systems: Vec<usize>) -> Self {
        Self { aromatic_systems, ..Self::default() }
    }

    pub fn multicenter_bonds(multicenter_bonds: Vec<usize>) -> Self {
        Self { multicenter_bonds, ..Self::default() }
    }

    pub fn noncovalent_bonds(noncovalent_bonds: Vec<usize>) -> Self {
        Self { noncovalent_bonds, ..Self::default() }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DerivedPred {
    TotalCharge(ValueAst),
    TotalSpin(SpinStateAst),
    ValenceSum(ValueAst),
    AromaticElectronCount(ValueAst),
    RingSize(ValueAst),
    InRing,
    NotInRing,
    InRelation(RelationSym),
    NotInRelation(RelationSym),
}

impl DerivedPred {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::TotalCharge(v) | Self::ValenceSum(v) | Self::AromaticElectronCount(v) | Self::RingSize(v) => v.is_ground(),
            Self::TotalSpin(s) => s.is_ground(),
            Self::InRing | Self::NotInRing | Self::InRelation(_) | Self::NotInRelation(_) => true,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MatcherFlag {
    Injective,
    NonInjective,
    Induced,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RelationSym {
    Atoms,
    Bonds,
    DativeBonds,
    AromaticSystems,
    MulticenterBonds,
    NoncovalentBonds,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomConstraint {
    ValenceSum(ValueAst),
    AromaticValence(ValueAst),
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
            | Self::AromaticValence(v)
            | Self::MulticenterValence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::Connectivity(v)
            | Self::TotalHCount(v)
            | Self::RingCount(v)
            | Self::RingSize(v) => v.is_ground(),
            Self::InRing => true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
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

    fn charge(n: i64) -> MoleculeConstraint {
        MoleculeConstraint::Derived {
            predicate: DerivedPred::TotalCharge(ValueAst::Lit(n)),
            refs: RelationRefs::default(),
        }
    }

    fn in_ring_atom(atom: usize) -> MoleculeConstraint {
        MoleculeConstraint::Derived {
            predicate: DerivedPred::InRing,
            refs: RelationRefs::atoms(vec![AtomIdx(atom as u32)]),
        }
    }

    #[rstest]
    #[case::and_pair(
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        true,
    )]
    #[case::and_order_matters(
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        MoleculeConstraint::And(vec![in_ring_atom(1), charge(0)]),
        false,
    )]
    #[case::or_distinct_payload(
        MoleculeConstraint::Or(vec![charge(0), charge(1)]),
        MoleculeConstraint::Or(vec![charge(0), charge(2)]),
        false,
    )]
    #[case::not_idempotent_eq(
        MoleculeConstraint::Not(Box::new(charge(-1))),
        MoleculeConstraint::Not(Box::new(charge(-1))),
        true,
    )]
    #[case::and_or_distinct(
        MoleculeConstraint::And(vec![charge(0)]),
        MoleculeConstraint::Or(vec![charge(0)]),
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
            charge(-1),
            MoleculeConstraint::Not(Box::new(in_ring_atom(0))),
        ]);
        let outer = MoleculeConstraint::And(vec![charge(0), inner.clone()]);

        let MoleculeConstraint::And(children) = &outer else {
            panic!("expected And");
        };
        assert_eq!(children.len(), 2);
        assert_eq!(children[1], inner);
    }

    #[rstest]
    #[case::valence_sum_lit(AtomConstraint::ValenceSum(ValueAst::Lit(4)), true)]
    #[case::valence_sum_undetermined(AtomConstraint::ValenceSum(ValueAst::Undetermined), false)]
    #[case::aromatic_valence_lit(AtomConstraint::AromaticValence(ValueAst::Lit(3)), true)]
    #[case::aromatic_valence_set(AtomConstraint::AromaticValence(ValueAst::LitSet(vec![2, 3])), false)]
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
        AtomConstraint::AromaticValence(ValueAst::Lit(4)),
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
