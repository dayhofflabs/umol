//! Molecule-scope constraints, the `Constraint` combinator tree, and the
//! molecule-level `Constraints` store.

use std::mem;
use std::slice::Iter;

use super::super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::super::molecule::MoleculeAst;
use super::super::remap::IdxRemapping;
use super::super::spin::SpinStateAst;
use super::super::value::ValueAst;
use super::aromatic::AromaticSystemConstraint;
use super::atom::AtomConstraint;
use super::bond::BondConstraint;
use super::dative::DativeBondConstraint;
use super::multicenter::MulticenterBondConstraint;
use super::noncovalent::NoncovalentBondConstraint;
use super::relational::RelationalConstraint;

/// Tree node type: per-entity leaf, molecule-scope leaf, relational leaf, or
/// combinator. The bare entity-leaf forms appear only inside a combinator
/// (e.g. `And(Atom(..), Bond(..))`) or a molecule-scope predicate;
/// unconditional per-entity value-only constraints live inline on the entity
/// AST and are lifted there at DSL → AST conversion time. Cross-entity
/// ref-bearing constraints (e.g. a dative-bond donor identity, aromatic
/// system membership, noncovalent endpoints) live only at molecule scope
/// via `Relational`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Constraint {
    Atom(AtomIdx, AtomConstraint),
    Bond(BondIdx, BondConstraint),
    DativeBond(DativeBondIdx, DativeBondConstraint),
    AromaticSystem(AromaticSystemIdx, AromaticSystemConstraint),
    MulticenterBond(MulticenterBondIdx, MulticenterBondConstraint),
    NoncovalentBond(NoncovalentBondIdx, NoncovalentBondConstraint),
    Relational(RelationalConstraint),
    Molecule(MoleculeConstraint),
    And(Vec<Constraint>),
    Or(Vec<Constraint>),
    Not(Box<Constraint>),
}

impl Constraint {
    /// A constraint is vacuous when it asserts nothing: the entity-leaf
    /// variants delegate to their inner `is_undetermined`; molecule-scope
    /// leaves delegate to `MoleculeConstraint::is_vacuous`. Combinators
    /// (`And`/`Or`) are vacuous only when empty; `Not(c)` is treated as
    /// non-vacuous (negating a vacuous claim is a meaningful unsat claim).
    /// `Relational` and `NoncovalentBond` are always non-vacuous (no
    /// `Undetermined` payload to elide).
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::Atom(_, c) => c.is_undetermined(),
            Self::Bond(_, c) => c.is_undetermined(),
            Self::DativeBond(_, c) => c.is_undetermined(),
            Self::AromaticSystem(_, c) => c.is_undetermined(),
            Self::MulticenterBond(_, c) => c.is_undetermined(),
            Self::NoncovalentBond(_, _) => false,
            Self::Relational(_) => false,
            Self::Molecule(c) => c.is_vacuous(),
            Self::And(xs) | Self::Or(xs) => xs.is_empty(),
            Self::Not(_) => false,
        }
    }

    /// Recursively simplify every contained `ValueAst`. Refs are unchanged;
    /// constraint kinds are preserved. SubPattern's inner `MoleculeAst` is
    /// recursively simplified via [`MoleculeAst::simplify_values`].
    pub fn simplify(self) -> Self {
        match self {
            Constraint::Atom(idx, c) => Constraint::Atom(idx, c.simplify()),
            Constraint::Bond(idx, c) => Constraint::Bond(idx, c.simplify()),
            Constraint::DativeBond(idx, c) => Constraint::DativeBond(idx, c.simplify()),
            Constraint::AromaticSystem(idx, c) => Constraint::AromaticSystem(idx, c.simplify()),
            Constraint::MulticenterBond(idx, c) => Constraint::MulticenterBond(idx, c.simplify()),
            Constraint::NoncovalentBond(_, c) => match c {},
            Constraint::Relational(r) => Constraint::Relational(r.simplify()),
            Constraint::Molecule(m) => Constraint::Molecule(m.simplify()),
            Constraint::And(xs) => Constraint::And(xs.into_iter().map(|c| c.simplify()).collect()),
            Constraint::Or(xs) => Constraint::Or(xs.into_iter().map(|c| c.simplify()).collect()),
            Constraint::Not(c) => Constraint::Not(Box::new((*c).simplify())),
        }
    }

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
            Constraint::Relational(r) => r.remap(remap).map(Constraint::Relational),
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

    /// Simplify every contained constraint's value(s) in place by
    /// recursively calling [`Constraint::simplify`].
    pub fn simplify_each(&mut self) {
        self.0 = mem::take(&mut self.0)
            .into_iter()
            .map(|c| c.simplify())
            .collect();
    }
}

/// Molecule-scope predicates: non-logical, unanchored assertions whose scope
/// is the molecule as a whole or a declared subset of entities.
///
/// For `ChargeSum` / `SpinSum` / `BondOrderSum` / `Connected`, an `atoms`
/// (or `bonds`) value of `None` denotes the entire molecule's atoms (or
/// bonds), making the predicate stable across structural growth. `Some(vec)`
/// denotes a fixed subset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    ChargeSum {
        atoms: Option<Vec<AtomIdx>>,
        sum: ValueAst,
    },
    SpinSum {
        atoms: Option<Vec<AtomIdx>>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Option<Vec<BondIdx>>,
        sum: ValueAst,
    },
    Connected {
        atoms: Option<Vec<AtomIdx>>,
    },
    SubPattern {
        anchor: SubPatternAnchor,
        pattern: Box<MoleculeAst>,
    },
}

impl MoleculeConstraint {
    /// A constraint is vacuous when its value-bearing payload is
    /// `Undetermined`: `ChargeSum`/`BondOrderSum` with `Undetermined` sum,
    /// `SpinSum` with both spin fields `Undetermined`. `Connected` and
    /// `SubPattern` are structural — never vacuous in this sense.
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::ChargeSum { sum, .. } => sum.is_undetermined(),
            Self::BondOrderSum { sum, .. } => sum.is_undetermined(),
            Self::SpinSum { spin, .. } => spin.is_undetermined(),
            Self::Connected { .. } => false,
            Self::SubPattern { .. } => false,
        }
    }

    /// Simplify every contained `ValueAst` and `SpinStateAst` in place.
    /// `Connected` carries no values to simplify; `SubPattern`'s pattern
    /// recurses via [`MoleculeAst::simplify_values`].
    pub fn simplify(self) -> Self {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms,
                sum: sum.simplify(),
            },
            MoleculeConstraint::SpinSum { atoms, mut spin } => {
                spin.simplify_values();
                MoleculeConstraint::SpinSum { atoms, spin }
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds,
                sum: sum.simplify(),
            },
            MoleculeConstraint::Connected { atoms } => MoleculeConstraint::Connected { atoms },
            MoleculeConstraint::SubPattern {
                anchor,
                mut pattern,
            } => {
                pattern.simplify_values();
                MoleculeConstraint::SubPattern { anchor, pattern }
            }
        }
    }

    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => {
                let atoms = remap_atom_subset(atoms, remap)?;
                Some(MoleculeConstraint::ChargeSum { atoms, sum })
            }
            MoleculeConstraint::SpinSum { atoms, spin } => {
                let atoms = remap_atom_subset(atoms, remap)?;
                Some(MoleculeConstraint::SpinSum { atoms, spin })
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => {
                let bonds = remap_bond_subset(bonds, remap)?;
                Some(MoleculeConstraint::BondOrderSum { bonds, sum })
            }
            MoleculeConstraint::Connected { atoms } => {
                let atoms = remap_atom_subset(atoms, remap)?;
                Some(MoleculeConstraint::Connected { atoms })
            }
            MoleculeConstraint::SubPattern { anchor, pattern } => anchor
                .remap(remap)
                .map(|anchor| MoleculeConstraint::SubPattern { anchor, pattern }),
        }
    }
}

/// Remap an `Option<Vec<AtomIdx>>`. `None` (all atoms) passes through.
/// `Some(vec)` remaps each element; if any atom was removed the whole
/// constraint is dropped (returns outer `None`).
fn remap_atom_subset(
    atoms: Option<Vec<AtomIdx>>,
    remap: &IdxRemapping,
) -> Option<Option<Vec<AtomIdx>>> {
    match atoms {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|a| remap.atom(a))
            .collect::<Option<Vec<_>>>()
            .map(Some),
    }
}

/// Remap an `Option<Vec<BondIdx>>`. Same semantics as `remap_atom_subset`.
fn remap_bond_subset(
    bonds: Option<Vec<BondIdx>>,
    remap: &IdxRemapping,
) -> Option<Option<Vec<BondIdx>>> {
    match bonds {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|b| remap.bond(b))
            .collect::<Option<Vec<_>>>()
            .map(Some),
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
    use crate::ast::idx::{
        AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
    };
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

    fn relation_remapping(
        removed_dative: Vec<u32>,
        removed_aromatic: Vec<u32>,
        removed_multicenter: Vec<u32>,
        removed_noncovalent: Vec<u32>,
    ) -> IdxRemapping {
        IdxRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            removed_dative,
            removed_aromatic,
            removed_multicenter,
            removed_noncovalent,
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(vec![], 0)]
    #[case::molecule_leaves(vec![Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomIdx(0), AtomIdx(1)]), sum: ValueAst::Lit(0) }),
            Constraint::Molecule(MoleculeConstraint::SpinSum { atoms: Some(vec![AtomIdx(0)]), spin: SpinStateAst::from((0_u8, 1_u8)) })], 2)]
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
            atoms: Some(vec![AtomIdx(0)]),
            sum: ValueAst::Lit(0),
        }));
        cs.push(Constraint::And(vec![]));

        cs.retain(|c| matches!(c, Constraint::Molecule(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_constraints_take_drains() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomIdx(0), AtomIdx(1)]),
        }));
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomIdx(0)]),
            sum: ValueAst::Lit(0),
        }));

        let taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_constraints_clear() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomIdx(0)]),
        }));
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
        a.push_dative_bond(DativeBondIdx(4), DativeBondIdx(0));
        a.push_aromatic_system(AromaticSystemIdx(2), AromaticSystemIdx(0));
        a.push_multicenter_bond(MulticenterBondIdx(7), MulticenterBondIdx(2));
        a.push_noncovalent_bond(NoncovalentBondIdx(1), NoncovalentBondIdx(3));

        assert!(!a.is_empty());
        assert_eq!(a.atoms(), &[(AtomIdx(3), AtomIdx(0))]);
        assert_eq!(a.bonds(), &[(BondIdx(5), BondIdx(1))]);
        assert_eq!(a.dative_bonds(), &[(DativeBondIdx(4), DativeBondIdx(0))]);
        assert_eq!(
            a.aromatic_systems(),
            &[(AromaticSystemIdx(2), AromaticSystemIdx(0))],
        );
        assert_eq!(
            a.multicenter_bonds(),
            &[(MulticenterBondIdx(7), MulticenterBondIdx(2))],
        );
        assert_eq!(
            a.noncovalent_bonds(),
            &[(NoncovalentBondIdx(1), NoncovalentBondIdx(3))],
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

    #[rstest]
    fn test_sub_pattern_anchor_remap_relations_shift() {
        let mut a = SubPatternAnchor::new();
        a.push_dative_bond(DativeBondIdx(2), DativeBondIdx(0));
        a.push_aromatic_system(AromaticSystemIdx(3), AromaticSystemIdx(1));
        a.push_multicenter_bond(MulticenterBondIdx(4), MulticenterBondIdx(2));
        a.push_noncovalent_bond(NoncovalentBondIdx(5), NoncovalentBondIdx(3));

        let remap = relation_remapping(vec![0], vec![1], vec![0], vec![2]);
        let mapped = a.remap(&remap).expect("all targets survive");

        assert_eq!(
            mapped.dative_bonds(),
            &[(DativeBondIdx(1), DativeBondIdx(0))]
        );
        assert_eq!(
            mapped.aromatic_systems(),
            &[(AromaticSystemIdx(2), AromaticSystemIdx(1))]
        );
        assert_eq!(
            mapped.multicenter_bonds(),
            &[(MulticenterBondIdx(3), MulticenterBondIdx(2))]
        );
        assert_eq!(
            mapped.noncovalent_bonds(),
            &[(NoncovalentBondIdx(4), NoncovalentBondIdx(3))]
        );
    }

    #[rstest]
    #[case::dative_dropped({
        let mut a = SubPatternAnchor::new();
        a.push_dative_bond(DativeBondIdx(1), DativeBondIdx(0));
        a
    }, relation_remapping(vec![1], vec![], vec![], vec![]))]
    #[case::aromatic_dropped({
        let mut a = SubPatternAnchor::new();
        a.push_aromatic_system(AromaticSystemIdx(2), AromaticSystemIdx(0));
        a
    }, relation_remapping(vec![], vec![2], vec![], vec![]))]
    #[case::multicenter_dropped({
        let mut a = SubPatternAnchor::new();
        a.push_multicenter_bond(MulticenterBondIdx(3), MulticenterBondIdx(0));
        a
    }, relation_remapping(vec![], vec![], vec![3], vec![]))]
    #[case::noncovalent_dropped({
        let mut a = SubPatternAnchor::new();
        a.push_noncovalent_bond(NoncovalentBondIdx(4), NoncovalentBondIdx(0));
        a
    }, relation_remapping(vec![], vec![], vec![], vec![4]))]
    fn test_sub_pattern_anchor_remap_relation_drops_return_none(
        #[case] anchor: SubPatternAnchor,
        #[case] remap: IdxRemapping,
    ) {
        assert_eq!(anchor.remap(&remap), None);
    }

    #[rstest]
    #[case::bond_shift(
        Constraint::Bond(BondIdx(3), BondConstraint::Aromatic),
        idx_remapping(vec![], vec![1]),
        Some(Constraint::Bond(BondIdx(2), BondConstraint::Aromatic)),
    )]
    #[case::bond_dropped(
        Constraint::Bond(BondIdx(1), BondConstraint::Aromatic),
        idx_remapping(vec![], vec![1]),
        None,
    )]
    #[case::or_all_survive(
        Constraint::Or(vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Atom(AtomIdx(2), AtomConstraint::Valence(ValueAst::Lit(2))),
        ]),
        idx_remapping(vec![1], vec![]),
        Some(Constraint::Or(vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Atom(AtomIdx(1), AtomConstraint::Valence(ValueAst::Lit(2))),
        ])),
    )]
    #[case::or_drops_if_any_leaf_drops(
        Constraint::Or(vec![
            Constraint::Atom(AtomIdx(1), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Atom(AtomIdx(2), AtomConstraint::Valence(ValueAst::Lit(2))),
        ]),
        idx_remapping(vec![1], vec![]),
        None,
    )]
    #[case::not_wraps_child(
        Constraint::Not(Box::new(Constraint::Atom(
            AtomIdx(2),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        ))),
        idx_remapping(vec![1], vec![]),
        Some(Constraint::Not(Box::new(Constraint::Atom(
            AtomIdx(1),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        )))),
    )]
    #[case::not_drops_child(
        Constraint::Not(Box::new(Constraint::Atom(
            AtomIdx(1),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        ))),
        idx_remapping(vec![1], vec![]),
        None,
    )]
    fn test_constraint_remap_logical_and_bond(
        #[case] c: Constraint,
        #[case] remap: IdxRemapping,
        #[case] expected: Option<Constraint>,
    ) {
        assert_eq!(c.remap(&remap), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor_shifts_inner_atom(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondIdx(2), atom: AtomIdx(3) }),
        idx_remapping(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondIdx(2), atom: AtomIdx(2) })),
    )]
    #[case::dative_donor_dropped_when_bond_removed(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondIdx(1), atom: AtomIdx(0) }),
        relation_remapping(vec![1], vec![], vec![], vec![]),
        None,
    )]
    #[case::aromatic_system_contains_shifts_inner_atom(
        Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemIdx(1), atom: AtomIdx(2) }),
        idx_remapping(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemIdx(1), atom: AtomIdx(1) })),
    )]
    #[case::multicenter_bond_contains_shifts_inner_atom(
        Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondIdx(0), atom: AtomIdx(2) }),
        idx_remapping(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondIdx(0), atom: AtomIdx(1) })),
    )]
    #[case::noncovalent_bond_contains_shifts_inner_atom(
        Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondIdx(0), atom: AtomIdx(3) }),
        idx_remapping(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondIdx(0), atom: AtomIdx(2) })),
    )]
    fn test_constraint_remap_relation_variants(
        #[case] c: Constraint,
        #[case] remap: IdxRemapping,
        #[case] expected: Option<Constraint>,
    ) {
        assert_eq!(c.remap(&remap), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_shifts(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomIdx(0), AtomIdx(2)]), sum: ValueAst::Lit(1) },
        idx_remapping(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomIdx(0), AtomIdx(1)]), sum: ValueAst::Lit(1) }),
    )]
    #[case::charge_sum_drops_when_atom_removed(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomIdx(1), AtomIdx(2)]), sum: ValueAst::Lit(0) },
        idx_remapping(vec![1], vec![]),
        None,
    )]
    #[case::charge_sum_all_atoms_passes_through(
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
        idx_remapping(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }),
    )]
    #[case::spin_sum_shifts(
        MoleculeConstraint::SpinSum { atoms: Some(vec![AtomIdx(0), AtomIdx(2)]), spin: SpinStateAst::from((0_u8, 1_u8)) },
        idx_remapping(vec![1], vec![]),
        Some(MoleculeConstraint::SpinSum { atoms: Some(vec![AtomIdx(0), AtomIdx(1)]), spin: SpinStateAst::from((0_u8, 1_u8)) }),
    )]
    #[case::bond_order_sum_shifts(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondIdx(0), BondIdx(2)]), sum: ValueAst::Lit(4) },
        idx_remapping(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondIdx(0), BondIdx(1)]), sum: ValueAst::Lit(4) }),
    )]
    #[case::bond_order_sum_drops_when_bond_removed(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondIdx(1)]), sum: ValueAst::Lit(2) },
        idx_remapping(vec![], vec![1]),
        None,
    )]
    #[case::bond_order_sum_all_bonds_passes_through(
        MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(0) },
        idx_remapping(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(0) }),
    )]
    #[case::connected_shifts(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomIdx(0), AtomIdx(2), AtomIdx(3)]) },
        idx_remapping(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: Some(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]) }),
    )]
    #[case::connected_drops_when_atom_removed(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomIdx(1)]) },
        idx_remapping(vec![1], vec![]),
        None,
    )]
    #[case::connected_all_atoms_passes_through(
        MoleculeConstraint::Connected { atoms: None },
        idx_remapping(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: None }),
    )]
    fn test_molecule_constraint_remap_variants(
        #[case] c: MoleculeConstraint,
        #[case] remap: IdxRemapping,
        #[case] expected: Option<MoleculeConstraint>,
    ) {
        assert_eq!(c.remap(&remap), expected);
    }
}
