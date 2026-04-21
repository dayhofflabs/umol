//! AST constraints: per-scope predicates and their routing.
//!
//! Per-scope enums (`AtomConstraint`, `BondConstraint`, `DativeBondConstraint`,
//! `AromaticSystemConstraint`, `MulticenterBondConstraint`,
//! `NoncovalentBondConstraint`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Storage is dual. Each entity AST carries its own inline `constraints`
//! vec; `Constraints` on `MoleculeAst` carries per-scope `IndexMap` buckets
//! plus a flat `molecule` vec for molecule-scope and combinator forms.
//! Consumers read the union of inline and molecule-level entries for any
//! given idx; there is no invariant between the two stores.

use std::hash::Hash;
use std::mem;

use indexmap::IndexMap;
use strum::EnumDiscriminants;

use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::molecule::MoleculeAst;
use super::remap::IdxRemapping;
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
        target_anchor: AtomIdx,
        pattern_anchor: AtomIdx,
        pattern: Box<MoleculeAst>,
    },
}

/// Tree node type: per-entity leaf, molecule-scope leaf, or combinator. Used
/// inside `Constraints::molecule`.
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

/// Molecule-level constraint storage. Independent of entity-inline
/// `constraints` vecs; consumers must read both.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Constraints {
    atom: IndexMap<AtomIdx, Vec<AtomConstraint>>,
    bond: IndexMap<BondIdx, Vec<BondConstraint>>,
    dative_bond: IndexMap<DativeBondIdx, Vec<DativeBondConstraint>>,
    aromatic_system: IndexMap<AromaticSystemIdx, Vec<AromaticSystemConstraint>>,
    multicenter_bond: IndexMap<MulticenterBondIdx, Vec<MulticenterBondConstraint>>,
    noncovalent_bond: IndexMap<NoncovalentBondIdx, Vec<NoncovalentBondConstraint>>,
    molecule: Vec<Constraint>,
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.atom.is_empty()
            && self.bond.is_empty()
            && self.dative_bond.is_empty()
            && self.aromatic_system.is_empty()
            && self.multicenter_bond.is_empty()
            && self.noncovalent_bond.is_empty()
            && self.molecule.is_empty()
    }

    pub fn len(&self) -> usize {
        self.atom.values().map(Vec::len).sum::<usize>()
            + self.bond.values().map(Vec::len).sum::<usize>()
            + self.dative_bond.values().map(Vec::len).sum::<usize>()
            + self.aromatic_system.values().map(Vec::len).sum::<usize>()
            + self.multicenter_bond.values().map(Vec::len).sum::<usize>()
            + self.noncovalent_bond.values().map(Vec::len).sum::<usize>()
            + self.molecule.len()
    }

    pub fn atom(&self, idx: AtomIdx) -> &[AtomConstraint] {
        self.atom.get(&idx).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn bond(&self, idx: BondIdx) -> &[BondConstraint] {
        self.bond.get(&idx).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn dative_bond(&self, idx: DativeBondIdx) -> &[DativeBondConstraint] {
        self.dative_bond.get(&idx).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn aromatic_system(&self, idx: AromaticSystemIdx) -> &[AromaticSystemConstraint] {
        self.aromatic_system
            .get(&idx)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondIdx) -> &[MulticenterBondConstraint] {
        self.multicenter_bond
            .get(&idx)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondIdx) -> &[NoncovalentBondConstraint] {
        self.noncovalent_bond
            .get(&idx)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn molecule(&self) -> &[Constraint] {
        &self.molecule
    }

    pub fn push_atom(&mut self, idx: AtomIdx, c: AtomConstraint) {
        self.atom.entry(idx).or_default().push(c);
    }

    pub fn push_bond(&mut self, idx: BondIdx, c: BondConstraint) {
        self.bond.entry(idx).or_default().push(c);
    }

    pub fn push_dative_bond(&mut self, idx: DativeBondIdx, c: DativeBondConstraint) {
        self.dative_bond.entry(idx).or_default().push(c);
    }

    pub fn push_aromatic_system(&mut self, idx: AromaticSystemIdx, c: AromaticSystemConstraint) {
        self.aromatic_system.entry(idx).or_default().push(c);
    }

    pub fn push_multicenter_bond(&mut self, idx: MulticenterBondIdx, c: MulticenterBondConstraint) {
        self.multicenter_bond.entry(idx).or_default().push(c);
    }

    pub fn push_noncovalent_bond(&mut self, idx: NoncovalentBondIdx, c: NoncovalentBondConstraint) {
        self.noncovalent_bond.entry(idx).or_default().push(c);
    }

    pub fn push_molecule(&mut self, c: Constraint) {
        self.molecule.push(c);
    }

    pub fn remove_atom(&mut self, idx: AtomIdx) -> Vec<AtomConstraint> {
        self.atom.shift_remove(&idx).unwrap_or_default()
    }

    pub fn remove_bond(&mut self, idx: BondIdx) -> Vec<BondConstraint> {
        self.bond.shift_remove(&idx).unwrap_or_default()
    }

    pub fn remove_dative_bond(&mut self, idx: DativeBondIdx) -> Vec<DativeBondConstraint> {
        self.dative_bond.shift_remove(&idx).unwrap_or_default()
    }

    pub fn remove_aromatic_system(
        &mut self,
        idx: AromaticSystemIdx,
    ) -> Vec<AromaticSystemConstraint> {
        self.aromatic_system.shift_remove(&idx).unwrap_or_default()
    }

    pub fn remove_multicenter_bond(
        &mut self,
        idx: MulticenterBondIdx,
    ) -> Vec<MulticenterBondConstraint> {
        self.multicenter_bond.shift_remove(&idx).unwrap_or_default()
    }

    pub fn remove_noncovalent_bond(
        &mut self,
        idx: NoncovalentBondIdx,
    ) -> Vec<NoncovalentBondConstraint> {
        self.noncovalent_bond.shift_remove(&idx).unwrap_or_default()
    }

    pub fn retain_atom(&mut self, idx: AtomIdx, mut f: impl FnMut(&AtomConstraint) -> bool) {
        if let Some(v) = self.atom.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.atom.shift_remove(&idx);
            }
        }
    }

    pub fn retain_bond(&mut self, idx: BondIdx, mut f: impl FnMut(&BondConstraint) -> bool) {
        if let Some(v) = self.bond.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.bond.shift_remove(&idx);
            }
        }
    }

    pub fn retain_dative_bond(
        &mut self,
        idx: DativeBondIdx,
        mut f: impl FnMut(&DativeBondConstraint) -> bool,
    ) {
        if let Some(v) = self.dative_bond.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.dative_bond.shift_remove(&idx);
            }
        }
    }

    pub fn retain_aromatic_system(
        &mut self,
        idx: AromaticSystemIdx,
        mut f: impl FnMut(&AromaticSystemConstraint) -> bool,
    ) {
        if let Some(v) = self.aromatic_system.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.aromatic_system.shift_remove(&idx);
            }
        }
    }

    pub fn retain_multicenter_bond(
        &mut self,
        idx: MulticenterBondIdx,
        mut f: impl FnMut(&MulticenterBondConstraint) -> bool,
    ) {
        if let Some(v) = self.multicenter_bond.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.multicenter_bond.shift_remove(&idx);
            }
        }
    }

    pub fn retain_noncovalent_bond(
        &mut self,
        idx: NoncovalentBondIdx,
        mut f: impl FnMut(&NoncovalentBondConstraint) -> bool,
    ) {
        if let Some(v) = self.noncovalent_bond.get_mut(&idx) {
            v.retain(|c| f(c));
            if v.is_empty() {
                self.noncovalent_bond.shift_remove(&idx);
            }
        }
    }

    pub fn remove_molecule(&mut self) -> Vec<Constraint> {
        mem::take(&mut self.molecule)
    }

    pub fn retain_molecule(&mut self, mut f: impl FnMut(&Constraint) -> bool) {
        self.molecule.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.atom.clear();
        self.bond.clear();
        self.dative_bond.clear();
        self.aromatic_system.clear();
        self.multicenter_bond.clear();
        self.noncovalent_bond.clear();
        self.molecule.clear();
    }

    pub fn remap(&mut self, remap: &IdxRemapping) {
        self.atom = remap_keys(mem::take(&mut self.atom), |k| remap.atom(k));
        self.bond = remap_keys(mem::take(&mut self.bond), |k| remap.bond(k));
        self.dative_bond = remap_keys(mem::take(&mut self.dative_bond), |k| remap.dative_bond(k));
        self.aromatic_system = remap_keys(mem::take(&mut self.aromatic_system), |k| {
            remap.aromatic_system(k)
        });
        self.multicenter_bond = remap_keys(mem::take(&mut self.multicenter_bond), |k| {
            remap.multicenter_bond(k)
        });
        self.noncovalent_bond = remap_keys(mem::take(&mut self.noncovalent_bond), |k| {
            remap.noncovalent_bond(k)
        });
        self.molecule = mem::take(&mut self.molecule)
            .into_iter()
            .filter_map(|c| c.remap(remap))
            .collect();
    }
}

impl Constraint {
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Constraint::Atom(idx, c) => remap.atom(idx).map(|i| Constraint::Atom(i, c)),
            Constraint::Bond(idx, c) => remap.bond(idx).map(|i| Constraint::Bond(i, c)),
            Constraint::DativeBond(idx, c) => {
                remap.dative_bond(idx).map(|i| Constraint::DativeBond(i, c))
            }
            Constraint::AromaticSystem(idx, c) => remap
                .aromatic_system(idx)
                .map(|i| Constraint::AromaticSystem(i, c)),
            Constraint::MulticenterBond(idx, c) => remap
                .multicenter_bond(idx)
                .map(|i| Constraint::MulticenterBond(i, c)),
            Constraint::NoncovalentBond(idx, c) => remap
                .noncovalent_bond(idx)
                .map(|i| Constraint::NoncovalentBond(i, c)),
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
            MoleculeConstraint::SubPattern {
                target_anchor,
                pattern_anchor,
                pattern,
            } => remap
                .atom(target_anchor)
                .map(|target_anchor| MoleculeConstraint::SubPattern {
                    target_anchor,
                    pattern_anchor,
                    pattern,
                }),
        }
    }
}

fn remap_keys<K: Hash + Eq + Copy, V>(
    map: IndexMap<K, V>,
    f: impl Fn(K) -> Option<K>,
) -> IndexMap<K, V> {
    map.into_iter()
        .filter_map(|(k, v)| f(k).map(|nk| (nk, v)))
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_constraints_push_per_scope() {
        let mut cs = Constraints::new();
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.push_atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(3)));
        cs.push_bond(BondIdx(0), BondConstraint::Aromatic);

        assert_eq!(cs.atom(AtomIdx(0)).len(), 2);
        assert_eq!(cs.bond(BondIdx(0)).len(), 1);
        assert_eq!(cs.dative_bond(DativeBondIdx(0)).len(), 0);
        assert_eq!(cs.len(), 3);
    }

    #[test]
    fn test_constraints_push_molecule() {
        let mut cs = Constraints::new();
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0), AtomIdx(1)],
            sum: ValueAst::Lit(0),
        }));
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::SpinSum {
            atoms: vec![AtomIdx(0)],
            spin: SpinStateAst::new(0, 1),
        }));
        assert_eq!(cs.molecule().len(), 2);
    }

    #[test]
    fn test_constraints_push_molecule_combinator() {
        let mut cs = Constraints::new();
        cs.push_molecule(Constraint::And(vec![
            Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4))),
            Constraint::Bond(BondIdx(0), BondConstraint::Aromatic),
        ]));
        assert_eq!(cs.molecule().len(), 1);
    }

    #[test]
    fn test_constraints_is_empty() {
        let mut cs = Constraints::new();
        assert!(cs.is_empty());
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
        assert!(!cs.is_empty());
    }

    #[test]
    fn test_constraints_remove_atom() {
        let mut cs = Constraints::new();
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.push_atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(3)));
        cs.push_atom(AtomIdx(1), AtomConstraint::Valence(ValueAst::Lit(2)));

        let removed = cs.remove_atom(AtomIdx(0));
        assert_eq!(
            removed,
            vec![
                AtomConstraint::Valence(ValueAst::Lit(4)),
                AtomConstraint::Degree(ValueAst::Lit(3)),
            ]
        );
        assert_eq!(cs.atom(AtomIdx(0)).len(), 0);
        assert_eq!(cs.atom(AtomIdx(1)).len(), 1);
    }

    #[test]
    fn test_constraints_remove_atom_missing() {
        let mut cs = Constraints::new();
        assert_eq!(cs.remove_atom(AtomIdx(0)), Vec::<AtomConstraint>::new());
    }

    #[test]
    fn test_constraints_remove_bond() {
        let mut cs = Constraints::new();
        cs.push_bond(BondIdx(0), BondConstraint::Aromatic);
        cs.push_bond(BondIdx(1), BondConstraint::RingCount(ValueAst::Lit(1)));

        let removed = cs.remove_bond(BondIdx(0));
        assert_eq!(removed, vec![BondConstraint::Aromatic]);
        assert_eq!(cs.bond(BondIdx(0)).len(), 0);
        assert_eq!(cs.bond(BondIdx(1)).len(), 1);
    }

    #[test]
    fn test_constraints_remove_molecule() {
        let mut cs = Constraints::new();
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0)],
            sum: ValueAst::Lit(0),
        }));
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
            AtomIdx(1),
        ])));

        let removed = cs.remove_molecule();
        assert_eq!(removed.len(), 2);
        assert_eq!(cs.molecule().len(), 0);
    }

    #[test]
    fn test_constraints_retain_atom() {
        let mut cs = Constraints::new();
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.push_atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(3)));

        cs.retain_atom(AtomIdx(0), |c| matches!(c, AtomConstraint::Valence(_)));
        assert_eq!(
            cs.atom(AtomIdx(0)),
            &[AtomConstraint::Valence(ValueAst::Lit(4))]
        );
    }

    #[test]
    fn test_constraints_retain_atom_empties_bucket() {
        let mut cs = Constraints::new();
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));

        cs.retain_atom(AtomIdx(0), |_| false);
        assert!(cs.is_empty());
    }

    #[test]
    fn test_constraints_retain_molecule() {
        let mut cs = Constraints::new();
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![AtomIdx(0)],
            sum: ValueAst::Lit(0),
        }));
        cs.push_molecule(Constraint::And(vec![]));

        cs.retain_molecule(|c| matches!(c, Constraint::Molecule(_)));
        assert_eq!(cs.molecule().len(), 1);
    }

    #[test]
    fn test_constraints_clear() {
        let mut cs = Constraints::new();
        cs.push_atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.push_bond(BondIdx(0), BondConstraint::Aromatic);
        cs.push_molecule(Constraint::Molecule(MoleculeConstraint::Connected(vec![
            AtomIdx(0),
        ])));

        cs.clear();
        assert!(cs.is_empty());
    }
}
