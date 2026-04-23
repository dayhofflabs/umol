//! Molecule-scope constraints, the `Constraint` combinator tree, and the
//! molecule-level `Constraints` store.

use std::mem;
use std::slice::Iter;

use super::aromatic::AromaticSystemConstraint;
use super::atom::AtomConstraint;
use super::bond::BondConstraint;
use super::dative::DativeBondConstraint;
use super::multicenter::MulticenterBondConstraint;
use super::noncovalent::NoncovalentBondConstraint;
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::remap::IdxRemapping;
use crate::ast::spin::SpinStateAst;
use crate::ast::value::ValueAst;

/// Tree node type: per-entity leaf, molecule-scope leaf, or combinator. The
/// bare entity-leaf forms appear only inside a combinator (e.g.
/// `And(Atom(..), Bond(..))`) or a molecule-scope predicate; unconditional
/// per-entity constraints live inline on the entity AST and are lifted there
/// at DSL → AST conversion time.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Constraint {
    Atom(AtomIdx, AtomConstraint),
    Bond(BondIdx, BondConstraint),
    DativeBond(DativeBondIdx, DativeBondConstraint),
    AromaticSystem(AromaticSystemIdx, AromaticSystemConstraint),
    MulticenterBond(MulticenterBondIdx, MulticenterBondConstraint),
    NoncovalentBond(NoncovalentBondIdx, NoncovalentBondConstraint),
    Molecule(MoleculeConstraint),
    And(Vec<Constraint>),
    Or(Vec<Constraint>),
    Not(Box<Constraint>),
}

impl Constraint {
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Constraint::Atom(idx, c) => remap.atom(idx).map(|i| Constraint::Atom(i, c)),
            Constraint::Bond(idx, c) => remap.bond(idx).map(|i| Constraint::Bond(i, c)),
            Constraint::DativeBond(idx, c) => {
                let i = remap.dative_bond(idx)?;
                c.remap(remap).map(|c| Constraint::DativeBond(i, c))
            }
            Constraint::AromaticSystem(idx, c) => {
                let i = remap.aromatic_system(idx)?;
                c.remap(remap).map(|c| Constraint::AromaticSystem(i, c))
            }
            Constraint::MulticenterBond(idx, c) => {
                let i = remap.multicenter_bond(idx)?;
                c.remap(remap).map(|c| Constraint::MulticenterBond(i, c))
            }
            Constraint::NoncovalentBond(idx, c) => {
                let i = remap.noncovalent_bond(idx)?;
                c.remap(remap).map(|c| Constraint::NoncovalentBond(i, c))
            }
            Constraint::Molecule(m) => m.remap(remap).map(Constraint::Molecule),
            Constraint::And(xs) => xs
                .into_iter()
                .map(|c| c.remap(remap))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::And),
            Constraint::Or(xs) => xs
                .into_iter()
                .map(|c| c.remap(remap))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::Or),
            Constraint::Not(x) => x.remap(remap).map(|c| Constraint::Not(Box::new(c))),
        }
    }
}

/// Molecule-level constraint store: a flat list of `Constraint` tree nodes
/// (molecule-scope predicates, combinators, and entity-leaves that appear
/// inside combinators). Unconditional per-entity constraints live on the
/// entity AST's own `constraints` field; the DSL parser lifts them there.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Constraints(Vec<Constraint>);

impl Constraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[Constraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, Constraint> {
        self.0.iter()
    }

    pub fn push(&mut self, c: Constraint) {
        self.0.push(c);
    }

    pub fn retain(&mut self, mut f: impl FnMut(&Constraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn take(&mut self) -> Vec<Constraint> {
        mem::take(&mut self.0)
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Remap entity indices. Entries that reference a removed entity (directly
    /// or via a combinator subtree) are dropped.
    pub fn remap(&mut self, remap: &IdxRemapping) {
        self.0 = mem::take(&mut self.0)
            .into_iter()
            .filter_map(|c| c.remap(remap))
            .collect();
    }
}

/// Molecule-scope predicates: non-logical, unanchored assertions whose scope
/// is the molecule as a whole or a declared subset of entities.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    ChargeSum {
        atoms: Vec<AtomIdx>,
        sum: ValueAst,
    },
    SpinSum {
        atoms: Vec<AtomIdx>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Vec<BondIdx>,
        sum: ValueAst,
    },
    Connected(Vec<AtomIdx>),
    SubPattern {
        anchor: SubPatternAnchor,
        pattern: Box<MoleculeAst>,
    },
}

impl MoleculeConstraint {
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(|atoms| MoleculeConstraint::ChargeSum { atoms, sum })
            }
            MoleculeConstraint::SpinSum { atoms, spin } => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(|atoms| MoleculeConstraint::SpinSum { atoms, spin })
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => {
                let bonds: Option<Vec<_>> = bonds.into_iter().map(|b| remap.bond(b)).collect();
                bonds.map(|bonds| MoleculeConstraint::BondOrderSum { bonds, sum })
            }
            MoleculeConstraint::Connected(atoms) => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(MoleculeConstraint::Connected)
            }
            MoleculeConstraint::SubPattern { anchor, pattern } => anchor
                .remap(remap)
                .map(|anchor| MoleculeConstraint::SubPattern { anchor, pattern }),
        }
    }
}

/// Multi-correspondence anchor for a `SubPattern` constraint. Each vec carries
/// `(target, pattern)` pairs pinning a target-molecule entity to a
/// pattern-molecule entity of the same kind. An empty anchor denotes an
/// unanchored match (pattern can embed anywhere).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SubPatternAnchor {
    atoms: Vec<(AtomIdx, AtomIdx)>,
    bonds: Vec<(BondIdx, BondIdx)>,
    dative_bonds: Vec<(DativeBondIdx, DativeBondIdx)>,
    aromatic_systems: Vec<(AromaticSystemIdx, AromaticSystemIdx)>,
    multicenter_bonds: Vec<(MulticenterBondIdx, MulticenterBondIdx)>,
    noncovalent_bonds: Vec<(NoncovalentBondIdx, NoncovalentBondIdx)>,
}

impl SubPatternAnchor {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
            && self.bonds.is_empty()
            && self.dative_bonds.is_empty()
            && self.aromatic_systems.is_empty()
            && self.multicenter_bonds.is_empty()
            && self.noncovalent_bonds.is_empty()
    }

    pub fn atoms(&self) -> &[(AtomIdx, AtomIdx)] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[(BondIdx, BondIdx)] {
        &self.bonds
    }

    pub fn dative_bonds(&self) -> &[(DativeBondIdx, DativeBondIdx)] {
        &self.dative_bonds
    }

    pub fn aromatic_systems(&self) -> &[(AromaticSystemIdx, AromaticSystemIdx)] {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &[(MulticenterBondIdx, MulticenterBondIdx)] {
        &self.multicenter_bonds
    }

    pub fn noncovalent_bonds(&self) -> &[(NoncovalentBondIdx, NoncovalentBondIdx)] {
        &self.noncovalent_bonds
    }

    pub fn push_atom(&mut self, target: AtomIdx, pattern: AtomIdx) {
        self.atoms.push((target, pattern));
    }

    pub fn push_bond(&mut self, target: BondIdx, pattern: BondIdx) {
        self.bonds.push((target, pattern));
    }

    pub fn push_dative_bond(&mut self, target: DativeBondIdx, pattern: DativeBondIdx) {
        self.dative_bonds.push((target, pattern));
    }

    pub fn push_aromatic_system(&mut self, target: AromaticSystemIdx, pattern: AromaticSystemIdx) {
        self.aromatic_systems.push((target, pattern));
    }

    pub fn push_multicenter_bond(
        &mut self,
        target: MulticenterBondIdx,
        pattern: MulticenterBondIdx,
    ) {
        self.multicenter_bonds.push((target, pattern));
    }

    pub fn push_noncovalent_bond(
        &mut self,
        target: NoncovalentBondIdx,
        pattern: NoncovalentBondIdx,
    ) {
        self.noncovalent_bonds.push((target, pattern));
    }

    /// Remap target-side indices per `remap`. Returns `None` if any target
    /// index in the anchor has been removed.
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .into_iter()
            .map(|(t, p)| remap.atom(t).map(|t| (t, p)))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .into_iter()
            .map(|(t, p)| remap.bond(t).map(|t| (t, p)))
            .collect();
        let dative_bonds: Option<Vec<_>> = self
            .dative_bonds
            .into_iter()
            .map(|(t, p)| remap.dative_bond(t).map(|t| (t, p)))
            .collect();
        let aromatic_systems: Option<Vec<_>> = self
            .aromatic_systems
            .into_iter()
            .map(|(t, p)| remap.aromatic_system(t).map(|t| (t, p)))
            .collect();
        let multicenter_bonds: Option<Vec<_>> = self
            .multicenter_bonds
            .into_iter()
            .map(|(t, p)| remap.multicenter_bond(t).map(|t| (t, p)))
            .collect();
        let noncovalent_bonds: Option<Vec<_>> = self
            .noncovalent_bonds
            .into_iter()
            .map(|(t, p)| remap.noncovalent_bond(t).map(|t| (t, p)))
            .collect();
        Some(Self {
            atoms: atoms?,
            bonds: bonds?,
            dative_bonds: dative_bonds?,
            aromatic_systems: aromatic_systems?,
            multicenter_bonds: multicenter_bonds?,
            noncovalent_bonds: noncovalent_bonds?,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::idx::{AromaticSystemIdx, AtomIdx, BondIdx};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::spin::SpinStateAst;
    use crate::ast::value::ValueAst;

    fn idx_remapping(removed_nodes: Vec<u32>, removed_edges: Vec<u32>) -> IdxRemapping {
        IdxRemapping::new(
            Remapping {
                removed_nodes,
                removed_edges,
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(vec![], 0)]
    #[case::molecule_leaves(vec![Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: vec![AtomIdx(0), AtomIdx(1)], sum: ValueAst::Lit(0) }),
            Constraint::Molecule(MoleculeConstraint::SpinSum { atoms: vec![AtomIdx(0)], spin: SpinStateAst::new(0, 1) })], 2)]
    #[case::combinator(vec![Constraint::And(vec![Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Bond(BondIdx(0), BondConstraint::Aromatic)])], 1)]
    fn test_constraints_push(
        #[case] items: Vec<Constraint>,
        #[case] expected_len: usize,
    ) {
        let mut cs = Constraints::new();
        for c in items {
            cs.push(c);
        }
        assert_eq!(cs.len(), expected_len);
        assert_eq!(cs.is_empty(), expected_len == 0);
    }

    #[rstest]
    fn test_constraints_retain() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0)],
            sum: ValueAst::Lit(0),
        }));
        cs.push(Constraint::And(vec![]));

        cs.retain(|c| matches!(c, Constraint::Molecule(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_constraints_take_drains() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
            AtomIdx(1),
        ])));
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0)],
            sum: ValueAst::Lit(0),
        }));

        let taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_constraints_clear() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
        ])));
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drops_entity_leaf_on_removed_atom(vec![Constraint::Atom(AtomIdx(1), AtomConstraint::Valence(ValueAst::Lit(4))),
        Constraint::Atom(AtomIdx(2), AtomConstraint::Valence(ValueAst::Lit(3)))],
        idx_remapping(vec![1], vec![]), vec![Constraint::Atom(AtomIdx(1), AtomConstraint::Valence(ValueAst::Lit(3)))])]
    #[case::shifts_remaining_leaves(vec![Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
        Constraint::Atom(AtomIdx(2), AtomConstraint::Degree(ValueAst::Lit(3)))],
        idx_remapping(vec![1], vec![]), vec![Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
        Constraint::Atom(AtomIdx(1), AtomConstraint::Degree(ValueAst::Lit(3)))])]
    #[case::drops_combinator_if_any_leaf_dropped(vec![Constraint::And(vec![Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
        Constraint::Atom(AtomIdx(1), AtomConstraint::Degree(ValueAst::Lit(3)))])], idx_remapping(vec![1], vec![]), vec![])]
    #[case::subpattern_shifts_anchor_atoms(vec![Constraint::Molecule(MoleculeConstraint::SubPattern {
        anchor: { let mut a = SubPatternAnchor::new(); a.push_atom(AtomIdx(3), AtomIdx(0)); a }, pattern: Box::new(MoleculeAst::default()) })],
        idx_remapping(vec![1], vec![]), vec![Constraint::Molecule(MoleculeConstraint::SubPattern { anchor: { let mut a = SubPatternAnchor::new(); a.push_atom(AtomIdx(2), AtomIdx(0)); a },
        pattern: Box::new(MoleculeAst::default()) })])]
    fn test_constraints_remap(
        #[case] items: Vec<Constraint>,
        #[case] remap: IdxRemapping,
        #[case] expected: Vec<Constraint>,
    ) {
        let mut cs = Constraints::new();
        for c in items {
            cs.push(c);
        }
        cs.remap(&remap);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_sub_pattern_anchor_default_is_empty() {
        let a = SubPatternAnchor::default();
        assert!(a.is_empty());
        assert!(a.atoms().is_empty());
        assert!(a.bonds().is_empty());
        assert!(a.dative_bonds().is_empty());
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_and_accessors() {
        let mut a = SubPatternAnchor::new();
        a.push_atom(AtomIdx(3), AtomIdx(0));
        a.push_bond(BondIdx(5), BondIdx(1));
        a.push_aromatic_system(AromaticSystemIdx(2), AromaticSystemIdx(0));

        assert!(!a.is_empty());
        assert_eq!(a.atoms(), &[(AtomIdx(3), AtomIdx(0))]);
        assert_eq!(a.bonds(), &[(BondIdx(5), BondIdx(1))]);
        assert_eq!(
            a.aromatic_systems(),
            &[(AromaticSystemIdx(2), AromaticSystemIdx(0))]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::shifts_target({ let mut a = SubPatternAnchor::new(); a.push_atom(AtomIdx(3), AtomIdx(0)); a.push_bond(BondIdx(5), BondIdx(1)); a },
        idx_remapping(vec![1], vec![2]), Some({ let mut a = SubPatternAnchor::new(); a.push_atom(AtomIdx(2), AtomIdx(0)); a.push_bond(BondIdx(4), BondIdx(1)); a }))]
    #[case::drops_on_removed_target_atom({ let mut a = SubPatternAnchor::new(); a.push_atom(AtomIdx(2), AtomIdx(0)); a }, idx_remapping(vec![2], vec![]), None)]
    fn test_sub_pattern_anchor_remap(
        #[case] anchor: SubPatternAnchor,
        #[case] remap: IdxRemapping,
        #[case] expected: Option<SubPatternAnchor>,
    ) {
        assert_eq!(anchor.remap(&remap), expected);
    }
}
