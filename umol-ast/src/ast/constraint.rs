//! AST constraints: per-scope predicates and their routing.
//!
//! Per-scope enums (`AtomConstraint`, `BondConstraint`, `AromaticSystemConstraint`,
//! `MulticenterBondConstraint`, `MoleculeConstraint`) each carry the
//! predicates admissible at that scope. `Constraint` routes per-entity
//! predicates to their target index, wraps molecule-scope predicates, and
//! provides logical combinators.
//!
//! Each AST node stores its own locally-scoped predicates in an inherent
//! `constraints: Vec<…>` field (collocation). During this migration phase
//! `Constraints` carries only molecule-global entries; lift/unlift
//! between local and global storage is future work.

use strum::EnumDiscriminants;

use super::idx::{AromaticSystemIdx, AtomIdx, BondIdx, MulticenterBondIdx};
use super::molecule::MoleculeAst;
use super::spin::SpinStateAst;
use super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash))]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceConstraint),
    MulticenterValence(MulticenterValenceConstraint),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    Connectivity(ValueAst),
    RingConnectivity(ValueAst),
    TotalHydrogens(ValueAst),
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl AtomConstraint {
    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticValenceConstraint {
    NotAromatic,
    Value(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterValenceConstraint {
    NotMulticenter,
    Value(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash))]
pub enum BondConstraint {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl BondConstraint {
    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl DativeBondConstraint {
    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraint {
    Electrons(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondConstraint {
    Electrons(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraint {}

/// Molecule-scope predicates: non-logical, unanchored assertions whose scope
/// is the molecule as a whole rather than a single entity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    ChargeSum(ValueAst),
    SpinSum(SpinStateAst),
    BondOrderSum {
        bonds: Vec<BondIdx>,
        sum: ValueAst,
    },
    Connected(Vec<AtomIdx>),
    SubPattern {
        target_anchor: AtomIdx,
        pattern_anchor: AtomIdx,
        pattern: Box<MoleculeAst>,
    },
}

/// Constraint entry on a `MoleculeAst`. Per-entity variants route a
/// scope-specific predicate to its target index; `Molecule` wraps
/// molecule-scope predicates; `And`/`Or`/`Not` are logical combinators.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Constraint {
    Atom(AtomIdx, AtomConstraint),
    Bond(BondIdx, BondConstraint),
    AromaticSystem(AromaticSystemIdx, AromaticSystemConstraint),
    MulticenterBond(MulticenterBondIdx, MulticenterBondConstraint),
    Molecule(MoleculeConstraint),
    And(Vec<Constraint>),
    Or(Vec<Constraint>),
    Not(Box<Constraint>),
}

/// Constraint structure for `MoleculeAst`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Constraints {
    items: Vec<Constraint>,
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, c: Constraint) {
        self.items.push(c);
    }

    pub fn iter(&self) -> impl Iterator<Item = &Constraint> {
        self.items.iter()
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn retain(&mut self, mut f: impl FnMut(&Constraint) -> bool) {
        self.items.retain(|c| f(c));
    }
}

impl FromIterator<Constraint> for Constraints {
    fn from_iter<I: IntoIterator<Item = Constraint>>(iter: I) -> Self {
        Self {
            items: iter.into_iter().collect(),
        }
    }
}

impl From<Vec<Constraint>> for Constraints {
    fn from(items: Vec<Constraint>) -> Self {
        Self { items }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_constraints_push_and_iter() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum(
            ValueAst::Lit(0),
        )));
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
            AtomIdx(1),
        ])));
        assert_eq!(cs.len(), 2);
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                Constraint::Molecule(MoleculeConstraint::ChargeSum(ValueAst::Lit(0))),
                Constraint::Molecule(MoleculeConstraint::Connected(vec![AtomIdx(0), AtomIdx(1),])),
            ],
        );
    }

    #[test]
    fn test_constraints_from_vec_roundtrip() {
        let items = vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Bond(BondIdx(0), BondConstraint::Aromatic),
            Constraint::Molecule(MoleculeConstraint::ChargeSum(ValueAst::Lit(0))),
        ];
        let cs = Constraints::from(items.clone());
        let roundtrip: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(roundtrip, items);
    }

    #[test]
    fn test_constraints_retain() {
        let mut cs = Constraints::from(vec![
            Constraint::Molecule(MoleculeConstraint::ChargeSum(ValueAst::Lit(0))),
            Constraint::Molecule(MoleculeConstraint::Connected(vec![AtomIdx(0)])),
        ]);
        cs.retain(|c| matches!(c, Constraint::Molecule(MoleculeConstraint::Connected(_)),));
        assert_eq!(cs.len(), 1);
    }
}
