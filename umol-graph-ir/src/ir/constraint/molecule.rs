//! Molecule-scope constraints, the `Constraint` combinator tree, and the
//! molecule-level `Constraints` store.

use std::cmp::Ordering;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use super::super::edit::{CascadedConstraints, ModifiedConstraint, RemovedConstraint};
use super::super::error::Contradiction;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::molecule::MoleculeAst;
use super::super::remap::{IdCompaction, IdRemapping};
use super::super::spin::UnpairedElectronsAst;
use super::super::stereo::StereoKind;
use super::super::traits::{Canonicalize, Lattice};
use super::super::value::ValueAst;
use super::aromatic::AromaticSystemConstraintAst;
use super::atom::AtomConstraintAst;
use super::bond::BondConstraintAst;
use super::dative::DativeBondConstraintAst;
use super::multicenter::MulticenterBondConstraintAst;
use super::noncovalent::NoncovalentBondConstraintAst;
use super::relational::RelationalConstraint;
use super::stereo::{StereoAtomConstraintAst, StereoBondConstraintAst};

/// Tree node type: per-entity leaf, molecule-scope leaf, relational leaf, or
/// combinator. The bare entity-leaf forms appear only inside a combinator
/// (e.g. `And(Atom(..), Bond(..))`) or a molecule-scope predicate;
/// unconditional per-entity value-only constraints live inline on the entity
/// AST and are lifted there at DSL → AST conversion time. Cross-entity
/// ref-bearing constraints (e.g. a dative-bond donor identity, aromatic
/// system membership, noncovalent endpoints) live only at molecule scope
/// via `Relational`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Constraint {
    Atom(AtomId, AtomConstraintAst),
    Bond(BondId, BondConstraintAst),
    DativeBond(DativeBondId, DativeBondConstraintAst),
    AromaticSystem(AromaticSystemId, AromaticSystemConstraintAst),
    MulticenterBond(MulticenterBondId, MulticenterBondConstraintAst),
    NoncovalentBond(NoncovalentBondId, NoncovalentBondConstraintAst),
    StereoAtom(StereoAtomId, StereoKind, StereoAtomConstraintAst),
    StereoBond(StereoBondId, StereoKind, StereoBondConstraintAst),
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
    /// `Relational` is always non-vacuous (no `Undetermined` payload to
    /// elide); the stereo leaves delegate to the inner `is_undetermined` (a
    /// `#o`/`#g` with an `Undetermined` relation constrains nothing).
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::Atom(_, c) => c.is_undetermined(),
            Self::Bond(_, c) => c.is_undetermined(),
            Self::DativeBond(_, c) => c.is_undetermined(),
            Self::AromaticSystem(_, c) => c.is_undetermined(),
            Self::MulticenterBond(_, c) => c.is_undetermined(),
            Self::NoncovalentBond(_, c) => c.is_undetermined(),
            Self::StereoAtom(_, _, c) => c.is_undetermined(),
            Self::StereoBond(_, _, c) => c.is_undetermined(),
            Self::Relational(_) => false,
            Self::Molecule(c) => c.is_vacuous(),
            Self::And(xs) | Self::Or(xs) => xs.is_empty(),
            Self::Not(_) => false,
        }
    }

    pub fn compact(self, compaction: &IdCompaction) -> Option<Self> {
        match self {
            Constraint::Atom(id, c) => {
                let i = compaction.compact_atom(id)?;
                c.compact(compaction).map(|c| Constraint::Atom(i, c))
            }
            Constraint::Bond(id, c) => {
                let i = compaction.compact_bond(id)?;
                c.compact(compaction).map(|c| Constraint::Bond(i, c))
            }
            Constraint::DativeBond(id, c) => {
                let i = compaction.compact_dative_bond(id)?;
                c.compact(compaction).map(|c| Constraint::DativeBond(i, c))
            }
            Constraint::AromaticSystem(id, c) => {
                let i = compaction.compact_aromatic_system(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::AromaticSystem(i, c))
            }
            Constraint::MulticenterBond(id, c) => {
                let i = compaction.compact_multicenter_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::MulticenterBond(i, c))
            }
            Constraint::NoncovalentBond(id, c) => {
                let i = compaction.compact_noncovalent_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::NoncovalentBond(i, c))
            }
            Constraint::StereoAtom(id, kind, c) => {
                let i = compaction.compact_stereo_atom(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::StereoAtom(i, kind, c))
            }
            Constraint::StereoBond(id, kind, c) => {
                let i = compaction.compact_stereo_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::StereoBond(i, kind, c))
            }
            Constraint::Relational(r) => r.compact(compaction).map(Constraint::Relational),
            Constraint::Molecule(m) => m.compact(compaction).map(Constraint::Molecule),
            Constraint::And(xs) => xs
                .into_iter()
                .map(|c| c.compact(compaction))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::And),
            Constraint::Or(xs) => xs
                .into_iter()
                .map(|c| c.compact(compaction))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::Or),
            Constraint::Not(x) => x.compact(compaction).map(|c| Constraint::Not(Box::new(c))),
        }
    }

    /// Re-anchor every entity ref through a total id remapping (match-based: lhs → host,
    /// created → appended). Total — never drops (the parallel of `compact`, which compacts after
    /// removal). Stays a separate flow from `compact`; the two are not bridged.
    pub fn remap(self, map: &IdRemapping) -> Self {
        match self {
            Constraint::Atom(id, c) => {
                let i = map.map_atom(id);
                Constraint::Atom(i, c.remap(map))
            }
            Constraint::Bond(id, c) => {
                let i = map.map_bond(id);
                Constraint::Bond(i, c.remap(map))
            }
            Constraint::DativeBond(id, c) => {
                Constraint::DativeBond(map.map_dative(id), c.remap(map))
            }
            Constraint::AromaticSystem(id, c) => {
                Constraint::AromaticSystem(map.map_aromatic(id), c.remap(map))
            }
            Constraint::MulticenterBond(id, c) => {
                Constraint::MulticenterBond(map.map_multicenter(id), c.remap(map))
            }
            Constraint::NoncovalentBond(id, c) => {
                Constraint::NoncovalentBond(map.map_noncovalent(id), c.remap(map))
            }
            Constraint::StereoAtom(id, kind, c) => {
                Constraint::StereoAtom(map.map_stereo_atom(id), kind, c.remap(map))
            }
            Constraint::StereoBond(id, kind, c) => {
                Constraint::StereoBond(map.map_stereo_bond(id), kind, c.remap(map))
            }
            Constraint::Relational(r) => Constraint::Relational(r.remap(map)),
            Constraint::Molecule(m) => Constraint::Molecule(m.remap(map)),
            Constraint::And(xs) => Constraint::And(xs.into_iter().map(|c| c.remap(map)).collect()),
            Constraint::Or(xs) => Constraint::Or(xs.into_iter().map(|c| c.remap(map)).collect()),
            Constraint::Not(x) => Constraint::Not(Box::new(x.remap(map))),
        }
    }
}

impl Canonicalize for Constraint {
    /// Canonicalize the inner predicate of each leaf; for `And`/`Or`, recurse,
    /// flatten the same combinator, drop empty `And`/`Or`, then sort + dedup
    /// children by the `Constraint` order. `Not` canonicalizes its inner node.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Atom(id, c) => Self::Atom(id, c.canonicalize()?),
            Self::Bond(id, c) => Self::Bond(id, c.canonicalize()?),
            Self::DativeBond(id, c) => Self::DativeBond(id, c.canonicalize()?),
            Self::AromaticSystem(id, c) => Self::AromaticSystem(id, c.canonicalize()?),
            Self::MulticenterBond(id, c) => Self::MulticenterBond(id, c.canonicalize()?),
            Self::NoncovalentBond(id, c) => Self::NoncovalentBond(id, c.canonicalize()?),
            Self::StereoAtom(id, kind, c) => Self::StereoAtom(id, kind, c.canonicalize()?),
            Self::StereoBond(id, kind, c) => Self::StereoBond(id, kind, c.canonicalize()?),
            Self::Relational(r) => Self::Relational(r.canonicalize()?),
            Self::Molecule(m) => Self::Molecule(m.canonicalize()?),
            Self::And(xs) => Self::And(canonicalize_logical_constraints(xs, true)?),
            Self::Or(xs) => Self::Or(canonicalize_logical_constraints(xs, false)?),
            Self::Not(c) => Self::Not(Box::new((*c).canonicalize()?)),
        })
    }
}

impl Canonicalize for Constraints {
    /// The store is an implicit conjunction, so it canonicalizes like an `And`:
    /// flatten top-level `And` entries, drop empty `And`/`Or`, sort + dedup.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self(canonicalize_logical_constraints(self.0, true)?))
    }
}

/// Canonicalize each child, splice same-combinator children (flatten), drop
/// empty `And`/`Or`, then sort + dedup. `is_and` selects which combinator is the
/// parent: `true` flattens nested `And` (and the conjunctive top-level store),
/// `false` flattens nested `Or`.
fn canonicalize_logical_constraints(
    constraints: Vec<Constraint>,
    is_and: bool,
) -> Result<Vec<Constraint>, Contradiction> {
    let mut out = Vec::new();
    for child in constraints {
        match child.canonicalize()? {
            Constraint::And(inner) if is_and => out.extend(inner),
            Constraint::Or(inner) if !is_and => out.extend(inner),
            Constraint::And(inner) | Constraint::Or(inner) if inner.is_empty() => {}
            other => out.push(other),
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

/// Molecule-level constraint store: a flat list of `Constraint` tree nodes
/// (molecule-scope predicates, combinators, and entity-leaves that appear
/// inside combinators). Unconditional per-entity constraints live on the
/// entity AST's own `constraints` field; the DSL parser lifts them there.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
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

    pub fn remove_at(&mut self, position: usize) -> Constraint {
        self.0.remove(position)
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
    pub fn compact(&mut self, compaction: &IdCompaction) {
        self.0 = mem::take(&mut self.0)
            .into_iter()
            .filter_map(|c| c.compact(compaction))
            .collect();
    }

    /// Remap entity indices and return the patch needed to restore or inspect
    /// constraints that were dropped or rewritten by the compaction.
    pub fn compact_with_update(&mut self, compaction: &IdCompaction) -> CascadedConstraints {
        let mut update = CascadedConstraints::default();
        let mut next = Vec::new();
        for (position, constraint) in mem::take(&mut self.0).into_iter().enumerate() {
            match constraint.clone().compact(compaction) {
                Some(mapped) => {
                    if mapped != constraint {
                        update.modified.push(ModifiedConstraint {
                            position,
                            old: constraint,
                            new: mapped.clone(),
                        });
                    }
                    next.push(mapped);
                }
                None => update.removed.push(RemovedConstraint {
                    position,
                    constraint,
                }),
            }
        }
        self.0 = next;
        update
    }
}

impl FromIterator<Constraint> for Constraints {
    fn from_iter<I: IntoIterator<Item = Constraint>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl IntoIterator for Constraints {
    type Item = Constraint;
    type IntoIter = IntoIter<Constraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<Constraint> for Constraints {
    fn from(c: Constraint) -> Self {
        Self(vec![c])
    }
}

impl From<Vec<Constraint>> for Constraints {
    fn from(cs: Vec<Constraint>) -> Self {
        Self(cs)
    }
}

/// Molecule-scope predicates: non-logical, unanchored assertions whose scope
/// is the molecule as a whole or a declared subset of entities.
///
/// For `ChargeSum` / `UnpairedElectronCoupling` / `BondOrderSum` / `Connected`,
/// an `atoms` (or `bonds`) value of `None` denotes the entire molecule's atoms
/// (or bonds), making the predicate stable across structural growth.
/// `Some(vec)` denotes a fixed subset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    ChargeSum {
        atoms: Option<Vec<AtomId>>,
        sum: ValueAst,
    },
    UnpairedElectronCoupling {
        atoms: Option<Vec<AtomId>>,
        unpaired_electrons: UnpairedElectronsAst,
    },
    BondOrderSum {
        bonds: Option<Vec<BondId>>,
        sum: ValueAst,
    },
    Connected {
        atoms: Option<Vec<AtomId>>,
    },
    SubPattern {
        anchor: SubPatternAnchor,
        pattern: Box<MoleculeAst>,
    },
}

impl MoleculeConstraint {
    /// A constraint is vacuous when its value-bearing payload is
    /// `Undetermined`: `ChargeSum`/`BondOrderSum` with `Undetermined` sum,
    /// `UnpairedElectronCoupling` with both unpaired-electron fields
    /// `Undetermined`. `Connected` and `SubPattern` are structural — never
    /// vacuous in this sense.
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::ChargeSum { sum, .. } => sum.is_undetermined(),
            Self::BondOrderSum { sum, .. } => sum.is_undetermined(),
            Self::UnpairedElectronCoupling {
                unpaired_electrons, ..
            } => unpaired_electrons.is_undetermined(),
            Self::Connected { .. } => false,
            Self::SubPattern { .. } => false,
        }
    }

    pub fn compact(self, compaction: &IdCompaction) -> Option<Self> {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::ChargeSum { atoms, sum })
            }
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::UnpairedElectronCoupling {
                    atoms,
                    unpaired_electrons,
                })
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => {
                let bonds = compact_bond_subset(bonds, compaction)?;
                Some(MoleculeConstraint::BondOrderSum { bonds, sum })
            }
            MoleculeConstraint::Connected { atoms } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::Connected { atoms })
            }
            MoleculeConstraint::SubPattern { anchor, pattern } => anchor
                .compact(compaction)
                .map(|anchor| MoleculeConstraint::SubPattern { anchor, pattern }),
        }
    }

    pub fn remap(self, map: &IdRemapping) -> Self {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms: remap_atom_subset(atoms, map),
                sum,
            },
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => MoleculeConstraint::UnpairedElectronCoupling {
                atoms: remap_atom_subset(atoms, map),
                unpaired_electrons,
            },
            MoleculeConstraint::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds: remap_bond_subset(bonds, map),
                sum,
            },
            MoleculeConstraint::Connected { atoms } => MoleculeConstraint::Connected {
                atoms: remap_atom_subset(atoms, map),
            },
            MoleculeConstraint::SubPattern { anchor, pattern } => MoleculeConstraint::SubPattern {
                anchor: anchor.remap(map),
                pattern,
            },
        }
    }
}

impl Canonicalize for MoleculeConstraint {
    /// Canonicalize the value payload (`sum` / `unpaired_electrons`) and sort
    /// each atom/bond subset; refs are otherwise unchanged. `SubPattern` is a
    /// **no-op** — the inner pattern is not recursed into (a nested pattern
    /// normalizes at its own top level via lift/inline and entity
    /// canonicalization).
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::ChargeSum { atoms, sum } => Self::ChargeSum {
                atoms: atoms.map(|mut v| {
                    v.sort_unstable();
                    v
                }),
                sum: sum.canonicalize()?,
            },
            Self::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => Self::UnpairedElectronCoupling {
                atoms: atoms.map(|mut v| {
                    v.sort_unstable();
                    v
                }),
                unpaired_electrons: unpaired_electrons.canonicalize()?,
            },
            Self::BondOrderSum { bonds, sum } => Self::BondOrderSum {
                bonds: bonds.map(|mut v| {
                    v.sort_unstable();
                    v
                }),
                sum: sum.canonicalize()?,
            },
            Self::Connected { atoms } => Self::Connected {
                atoms: atoms.map(|mut v| {
                    v.sort_unstable();
                    v
                }),
            },
            other @ Self::SubPattern { .. } => other,
        })
    }
}

impl PartialOrd for MoleculeConstraint {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MoleculeConstraint {
    /// Variant declaration order, then payload. `SubPattern` orders by `anchor`
    /// only — the inner `MoleculeAst` has no total order (graph), so same-anchor
    /// patterns compare `Equal` here. This is intentionally weaker than `Eq`,
    /// which does compare the pattern; canonicalization only needs the order for
    /// a stable sort, and dedup falls back to `PartialEq`.
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::ChargeSum { atoms: a1, sum: s1 }, Self::ChargeSum { atoms: a2, sum: s2 }) => {
                (a1, s1).cmp(&(a2, s2))
            }
            (
                Self::UnpairedElectronCoupling {
                    atoms: a1,
                    unpaired_electrons: s1,
                },
                Self::UnpairedElectronCoupling {
                    atoms: a2,
                    unpaired_electrons: s2,
                },
            ) => (a1, s1).cmp(&(a2, s2)),
            (
                Self::BondOrderSum { bonds: b1, sum: s1 },
                Self::BondOrderSum { bonds: b2, sum: s2 },
            ) => (b1, s1).cmp(&(b2, s2)),
            (Self::Connected { atoms: a1 }, Self::Connected { atoms: a2 }) => a1.cmp(a2),
            (Self::SubPattern { anchor: a1, .. }, Self::SubPattern { anchor: a2, .. }) => {
                a1.cmp(a2)
            }
            _ => {
                let rank = |c: &Self| match c {
                    Self::ChargeSum { .. } => 0u8,
                    Self::UnpairedElectronCoupling { .. } => 1,
                    Self::BondOrderSum { .. } => 2,
                    Self::Connected { .. } => 3,
                    Self::SubPattern { .. } => 4,
                };
                rank(self).cmp(&rank(other))
            }
        }
    }
}

/// Remap an `Option<Vec<AtomId>>`. `None` (all atoms) passes through.
/// `Some(vec)` compacts each element; if any atom was removed the whole
/// constraint is dropped (returns outer `None`).
fn compact_atom_subset(
    atoms: Option<Vec<AtomId>>,
    compaction: &IdCompaction,
) -> Option<Option<Vec<AtomId>>> {
    match atoms {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|a| compaction.compact_atom(a))
            .collect::<Option<Vec<_>>>()
            .map(Some),
    }
}

/// Remap an `Option<Vec<BondId>>`. Same semantics as `compact_atom_subset`.
fn compact_bond_subset(
    bonds: Option<Vec<BondId>>,
    compaction: &IdCompaction,
) -> Option<Option<Vec<BondId>>> {
    match bonds {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|b| compaction.compact_bond(b))
            .collect::<Option<Vec<_>>>()
            .map(Some),
    }
}

/// Re-anchor an `Option<Vec<AtomId>>` through a total id remapping (the parallel of
/// `compact_atom_subset`). `None` (all atoms) passes through; total — never drops.
fn remap_atom_subset(atoms: Option<Vec<AtomId>>, map: &IdRemapping) -> Option<Vec<AtomId>> {
    atoms.map(|vec| vec.into_iter().map(|a| map.map_atom(a)).collect())
}

/// Re-anchor an `Option<Vec<BondId>>`. Same semantics as `remap_atom_subset`.
fn remap_bond_subset(bonds: Option<Vec<BondId>>, map: &IdRemapping) -> Option<Vec<BondId>> {
    bonds.map(|vec| vec.into_iter().map(|b| map.map_bond(b)).collect())
}

/// Multi-correspondence anchor for a `SubPattern` constraint. Each vec carries
/// `(target, pattern)` pairs constraining a target-molecule entity to a
/// pattern-molecule entity of the same kind. An empty anchor denotes an
/// unanchored match (pattern can embed anywhere).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct SubPatternAnchor {
    atoms: Vec<(AtomId, AtomId)>,
    bonds: Vec<(BondId, BondId)>,
    dative_bonds: Vec<(DativeBondId, DativeBondId)>,
    aromatic_systems: Vec<(AromaticSystemId, AromaticSystemId)>,
    multicenter_bonds: Vec<(MulticenterBondId, MulticenterBondId)>,
    noncovalent_bonds: Vec<(NoncovalentBondId, NoncovalentBondId)>,
    stereo_atoms: Vec<(StereoAtomId, StereoAtomId)>,
    stereo_bonds: Vec<(StereoBondId, StereoBondId)>,
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
            && self.stereo_atoms.is_empty()
            && self.stereo_bonds.is_empty()
    }

    pub fn atoms(&self) -> &[(AtomId, AtomId)] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[(BondId, BondId)] {
        &self.bonds
    }

    pub fn dative_bonds(&self) -> &[(DativeBondId, DativeBondId)] {
        &self.dative_bonds
    }

    pub fn aromatic_systems(&self) -> &[(AromaticSystemId, AromaticSystemId)] {
        &self.aromatic_systems
    }

    pub fn multicenter_bonds(&self) -> &[(MulticenterBondId, MulticenterBondId)] {
        &self.multicenter_bonds
    }

    pub fn noncovalent_bonds(&self) -> &[(NoncovalentBondId, NoncovalentBondId)] {
        &self.noncovalent_bonds
    }

    pub fn stereo_atoms(&self) -> &[(StereoAtomId, StereoAtomId)] {
        &self.stereo_atoms
    }

    pub fn stereo_bonds(&self) -> &[(StereoBondId, StereoBondId)] {
        &self.stereo_bonds
    }

    pub fn push_atom(&mut self, target: AtomId, pattern: AtomId) {
        self.atoms.push((target, pattern));
    }

    pub fn push_bond(&mut self, target: BondId, pattern: BondId) {
        self.bonds.push((target, pattern));
    }

    pub fn push_dative_bond(&mut self, target: DativeBondId, pattern: DativeBondId) {
        self.dative_bonds.push((target, pattern));
    }

    pub fn push_aromatic_system(&mut self, target: AromaticSystemId, pattern: AromaticSystemId) {
        self.aromatic_systems.push((target, pattern));
    }

    pub fn push_multicenter_bond(&mut self, target: MulticenterBondId, pattern: MulticenterBondId) {
        self.multicenter_bonds.push((target, pattern));
    }

    pub fn push_noncovalent_bond(&mut self, target: NoncovalentBondId, pattern: NoncovalentBondId) {
        self.noncovalent_bonds.push((target, pattern));
    }

    pub fn push_stereo_atom(&mut self, target: StereoAtomId, pattern: StereoAtomId) {
        self.stereo_atoms.push((target, pattern));
    }

    pub fn push_stereo_bond(&mut self, target: StereoBondId, pattern: StereoBondId) {
        self.stereo_bonds.push((target, pattern));
    }

    /// Compact target-side indices per `compaction`. Returns `None` if any target
    /// index in the anchor has been removed.
    pub fn compact(self, compaction: &IdCompaction) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .into_iter()
            .map(|(t, p)| compaction.compact_atom(t).map(|t| (t, p)))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .into_iter()
            .map(|(t, p)| compaction.compact_bond(t).map(|t| (t, p)))
            .collect();
        let dative_bonds: Option<Vec<_>> = self
            .dative_bonds
            .into_iter()
            .map(|(t, p)| compaction.compact_dative_bond(t).map(|t| (t, p)))
            .collect();
        let aromatic_systems: Option<Vec<_>> = self
            .aromatic_systems
            .into_iter()
            .map(|(t, p)| compaction.compact_aromatic_system(t).map(|t| (t, p)))
            .collect();
        let multicenter_bonds: Option<Vec<_>> = self
            .multicenter_bonds
            .into_iter()
            .map(|(t, p)| compaction.compact_multicenter_bond(t).map(|t| (t, p)))
            .collect();
        let noncovalent_bonds: Option<Vec<_>> = self
            .noncovalent_bonds
            .into_iter()
            .map(|(t, p)| compaction.compact_noncovalent_bond(t).map(|t| (t, p)))
            .collect();
        let stereo_atoms: Option<Vec<_>> = self
            .stereo_atoms
            .into_iter()
            .map(|(t, p)| compaction.compact_stereo_atom(t).map(|t| (t, p)))
            .collect();
        let stereo_bonds: Option<Vec<_>> = self
            .stereo_bonds
            .into_iter()
            .map(|(t, p)| compaction.compact_stereo_bond(t).map(|t| (t, p)))
            .collect();
        Some(Self {
            atoms: atoms?,
            bonds: bonds?,
            dative_bonds: dative_bonds?,
            aromatic_systems: aromatic_systems?,
            multicenter_bonds: multicenter_bonds?,
            noncovalent_bonds: noncovalent_bonds?,
            stereo_atoms: stereo_atoms?,
            stereo_bonds: stereo_bonds?,
        })
    }

    /// Re-anchor target-side indices through a total id remapping. Total — the parallel of
    /// `compact`, never drops.
    pub(crate) fn remap(self, map: &IdRemapping) -> Self {
        Self {
            atoms: self
                .atoms
                .into_iter()
                .map(|(t, p)| (map.map_atom(t), p))
                .collect(),
            bonds: self
                .bonds
                .into_iter()
                .map(|(t, p)| (map.map_bond(t), p))
                .collect(),
            dative_bonds: self
                .dative_bonds
                .into_iter()
                .map(|(t, p)| (map.map_dative(t), p))
                .collect(),
            aromatic_systems: self
                .aromatic_systems
                .into_iter()
                .map(|(t, p)| (map.map_aromatic(t), p))
                .collect(),
            multicenter_bonds: self
                .multicenter_bonds
                .into_iter()
                .map(|(t, p)| (map.map_multicenter(t), p))
                .collect(),
            noncovalent_bonds: self
                .noncovalent_bonds
                .into_iter()
                .map(|(t, p)| (map.map_noncovalent(t), p))
                .collect(),
            stereo_atoms: self
                .stereo_atoms
                .into_iter()
                .map(|(t, p)| (map.map_stereo_atom(t), p))
                .collect(),
            stereo_bonds: self
                .stereo_bonds
                .into_iter()
                .map(|(t, p)| (map.map_stereo_bond(t), p))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{Compaction, RelationId};

    use super::*;
    use crate::ir::atom::AtomAst;
    use crate::ir::constraint::RingScope;
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    };
    use crate::ir::molecule::{MoleculeAst, MoleculeEntries};
    use crate::ir::spin::UnpairedElectronsAst;
    use crate::ir::value::{ValueAst, ValueTerm};
    use crate::ir::BooleanAst;

    fn id_compaction(removed_nodes: Vec<u32>, removed_edges: Vec<u32>) -> IdCompaction {
        IdCompaction::new(
            Compaction::new(removed_nodes, removed_edges),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    fn relation_compaction(
        removed_dative: Vec<u32>,
        removed_aromatic: Vec<u32>,
        removed_multicenter: Vec<u32>,
        removed_noncovalent: Vec<u32>,
        removed_stereo_atoms: Vec<u32>,
        removed_stereo_bonds: Vec<u32>,
    ) -> IdCompaction {
        let rel = |v: Vec<u32>| v.into_iter().map(RelationId).collect::<Vec<_>>();
        IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            rel(removed_dative),
            rel(removed_aromatic),
            rel(removed_multicenter),
            rel(removed_noncovalent),
            rel(removed_stereo_atoms),
            rel(removed_stereo_bonds),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_lit(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)), false)]
    #[case::atom_undetermined(Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::Undetermined)), true)]
    #[case::bond_lit(Constraint::Bond(BondId(0), BondConstraintAst::ring_membership(RingScope::Size(6), 1)), false)]
    #[case::bond_undetermined(Constraint::Bond(BondId(0), BondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)), true)]
    #[case::bond_aromatic_flag(Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))), false)]
    #[case::dative_undetermined(Constraint::DativeBond(DativeBondId(0), DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)), true)]
    #[case::aromatic_system_undetermined(Constraint::AromaticSystem(AromaticSystemId(0),
        AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)), true)]
    #[case::multicenter_undetermined(Constraint::MulticenterBond(MulticenterBondId(0),
        MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined)), true)]
    #[case::relational(Constraint::Relational(RelationalConstraint::DativeBondDonor {
        bond: DativeBondId(0), atom: AtomId(0) }), false)]
    #[case::molecule_undetermined(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: ValueAst::Undetermined }), true)]
    #[case::molecule_lit(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: ValueAst::Lit(0) }), false)]
    #[case::and_empty(Constraint::And(vec![]), true)]
    #[case::and_nonempty(Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))]), false)]
    #[case::or_empty(Constraint::Or(vec![]), true)]
    #[case::or_nonempty(Constraint::Or(vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))]), false)]
    #[case::not(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)))), false)]
    fn test_constraint_is_vacuous(#[case] c: Constraint, #[case] expected: bool) {
        assert_eq!(c.is_vacuous(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::leaf_folds(
        Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::term(ValueTerm::Lit(4)))),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))),
    )]
    #[case::and_flattens_nested(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::And(vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))]),
        ]),
        Ok(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ])),
    )]
    #[case::and_drops_empty_or_child(
        Constraint::And(vec![Constraint::Or(vec![]), Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))]),
        Ok(Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))])),
    )]
    #[case::and_sorts_and_dedups(
        Constraint::And(vec![
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
        ]),
        Ok(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ])),
    )]
    #[case::or_flattens_nested(
        Constraint::Or(vec![
            Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))]),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ]),
        Ok(Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ])),
    )]
    #[case::or_drops_empty_and_child(
        Constraint::Or(vec![Constraint::And(vec![]), Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))]),
        Ok(Constraint::Or(vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))])),
    )]
    #[case::not_folds_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::term(ValueTerm::Lit(4)))))),
        Ok(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4))))),
    )]
    #[case::inner_contradiction_propagates(
        Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::lit_set(Vec::<i64>::new())))]),
        Err(Contradiction),
    )]
    fn test_constraint_canonicalize(
        #[case] input: Constraint,
        #[case] expected: Result<Constraint, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_shifts(
        Constraint::Atom(AtomId(2), AtomConstraintAst::valence(4)),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4))),
    )]
    #[case::atom_dropped(
        Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::bond_shifts(
        Constraint::Bond(BondId(3), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        id_compaction(vec![], vec![1]),
        Some(Constraint::Bond(BondId(2), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))),
    )]
    #[case::bond_dropped(
        Constraint::Bond(BondId(1), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        id_compaction(vec![], vec![1]),
        None,
    )]
    #[case::dative_shifts(
        Constraint::DativeBond(DativeBondId(2), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        relation_compaction(vec![0], vec![], vec![], vec![], vec![], vec![]),
        Some(Constraint::DativeBond(DativeBondId(1), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)))),
    )]
    #[case::dative_dropped(
        Constraint::DativeBond(DativeBondId(1), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        relation_compaction(vec![1], vec![], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::aromatic_system_shifts(
        Constraint::AromaticSystem(AromaticSystemId(2), AromaticSystemConstraintAst::electron_count(6)),
        relation_compaction(vec![], vec![0], vec![], vec![], vec![], vec![]),
        Some(Constraint::AromaticSystem(AromaticSystemId(1), AromaticSystemConstraintAst::electron_count(6))),
    )]
    #[case::aromatic_system_dropped(
        Constraint::AromaticSystem(AromaticSystemId(1), AromaticSystemConstraintAst::electron_count(6)),
        relation_compaction(vec![], vec![1], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::multicenter_shifts(
        Constraint::MulticenterBond(MulticenterBondId(2), MulticenterBondConstraintAst::electron_count(2)),
        relation_compaction(vec![], vec![], vec![0], vec![], vec![], vec![]),
        Some(Constraint::MulticenterBond(MulticenterBondId(1), MulticenterBondConstraintAst::electron_count(2))),
    )]
    #[case::multicenter_dropped(
        Constraint::MulticenterBond(MulticenterBondId(1), MulticenterBondConstraintAst::electron_count(2)),
        relation_compaction(vec![], vec![], vec![1], vec![], vec![], vec![]),
        None,
    )]
    #[case::relational_dative_donor_shifts_atom(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(2), atom: AtomId(3) }),
        id_compaction(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(2), atom: AtomId(2) })),
    )]
    #[case::relational_dative_donor_dropped_when_bond_removed(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(0) }),
        relation_compaction(vec![1], vec![], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::relational_aromatic_system_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(2) }),
        id_compaction(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(1) })),
    )]
    #[case::relational_multicenter_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(2) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(1) })),
    )]
    #[case::relational_noncovalent_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(3) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(2) })),
    )]
    #[case::molecule_charge_sum_shifts(
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: ValueAst::Lit(1) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: ValueAst::Lit(1) })),
    )]
    #[case::and_all_survive(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        Some(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(2)),
        ])),
    )]
    #[case::and_drops_if_any_leaf_drops(
        Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::or_all_survive(
        Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(2)),
        ])),
    )]
    #[case::or_drops_if_any_leaf_drops(
        Constraint::Or(vec![
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::not_wraps_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(2), AtomConstraintAst::valence(4)))),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Not(Box::new(Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4))))),
    )]
    #[case::not_drops_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)))),
        id_compaction(vec![1], vec![]),
        None,
    )]
    fn test_constraint_compact(
        #[case] c: Constraint,
        #[case] compaction: IdCompaction,
        #[case] expected: Option<Constraint>,
    ) {
        assert_eq!(c.compact(&compaction), expected);
    }

    fn id_remapping(
        atom: &[(u32, u32)],
        bond: &[(u32, u32)],
        dative: &[(u32, u32)],
    ) -> IdRemapping {
        IdRemapping::new(
            atom.iter().map(|&(a, b)| (AtomId(a), AtomId(b))).collect(),
            bond.iter().map(|&(a, b)| (BondId(a), BondId(b))).collect(),
            dative
                .iter()
                .map(|&(a, b)| (DativeBondId(a), DativeBondId(b)))
                .collect(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom(
        Constraint::Atom(AtomId(2), AtomConstraintAst::valence(4)),
        id_remapping(&[(2, 5)], &[], &[]),
        Constraint::Atom(AtomId(5), AtomConstraintAst::valence(4)),
    )]
    #[case::bond(
        Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        id_remapping(&[], &[(0, 3)], &[]),
        Constraint::Bond(BondId(3), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
    )]
    #[case::dative_leaf(
        Constraint::DativeBond(DativeBondId(1), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        id_remapping(&[], &[], &[(1, 0)]),
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
    )]
    #[case::molecule_charge_sum(
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: ValueAst::Lit(1) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(3), AtomId(4)]), sum: ValueAst::Lit(1) }),
    )]
    #[case::molecule_unpaired_electron_coupling_subset(
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(3), AtomId(4)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
    )]
    #[case::molecule_unpaired_electron_coupling_all_atoms(
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
    )]
    #[case::relational_dative_donor(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(2) }),
        id_remapping(&[(2, 5)], &[], &[(1, 0)]),
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(5) }),
    )]
    #[case::and(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(2)),
        ]),
        id_remapping(&[(0, 1), (2, 3)], &[], &[]),
        Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(3), AtomConstraintAst::valence(2)),
        ]),
    )]
    #[case::not(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(2), AtomConstraintAst::valence(4)))),
        id_remapping(&[(2, 0)], &[], &[]),
        Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)))),
    )]
    fn test_constraint_remap(
        #[case] c: Constraint,
        #[case] map: IdRemapping,
        #[case] expected: Constraint,
    ) {
        assert_eq!(c.remap(&map), expected);
    }

    #[rstest]
    fn test_constraints_new() {
        let cs = Constraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[Constraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(vec![], 0)]
    #[case::molecule_leaves(vec![
            Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: ValueAst::Lit(0) }),
            Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
        ], 2)]
    #[case::combinator(vec![Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ])], 1)]
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
            atoms: Some(vec![AtomId(0)]),
            sum: ValueAst::Lit(0),
        }));
        cs.push(Constraint::And(vec![]));

        cs.retain(|c| matches!(c, Constraint::Molecule(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_constraints_take() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        }));
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
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
            atoms: Some(vec![AtomId(0)]),
        }));
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_constraints_iter() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)));
        cs.push(Constraint::Bond(
            BondId(0),
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
        ));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintAst::Aromatic(BooleanAst::Lit(true))
                ),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drops_entity_leaf_on_removed_atom(
        vec![
            Constraint::Atom(AtomId(1), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::valence(3)),
        ],
        id_compaction(vec![1], vec![]),
        vec![Constraint::Atom(AtomId(1), AtomConstraintAst::valence(3))],
    )]
    #[case::shifts_remaining_leaves(
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintAst::degree(3)),
        ],
        id_compaction(vec![1], vec![]),
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintAst::degree(3)),
        ],
    )]
    #[case::drops_combinator_if_any_leaf_dropped(
        vec![Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintAst::degree(3)),
        ])],
        id_compaction(vec![1], vec![]),
        vec![],
    )]
    #[case::subpattern_shifts_anchor_atoms(
        vec![Constraint::Molecule(MoleculeConstraint::SubPattern {
            anchor: { let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(3), AtomId(0)); a },
            pattern: Box::new(MoleculeAst::default()),
        })],
        id_compaction(vec![1], vec![]),
        vec![Constraint::Molecule(MoleculeConstraint::SubPattern {
            anchor: { let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(2), AtomId(0)); a },
            pattern: Box::new(MoleculeAst::default()),
        })],
    )]
    fn test_constraints_compact(
        #[case] items: Vec<Constraint>,
        #[case] compaction: IdCompaction,
        #[case] expected: Vec<Constraint>,
    ) {
        let mut cs = Constraints::new();
        for c in items {
            cs.push(c);
        }
        cs.compact(&compaction);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_constraints_compact_with_update() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)));
        cs.push(Constraint::Atom(AtomId(1), AtomConstraintAst::degree(3)));
        cs.push(Constraint::Bond(
            BondId(2),
            BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        ));
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(2)]),
        }));

        let update = cs.compact_with_update(&id_compaction(vec![1], vec![1]));

        assert_eq!(
            cs.as_slice(),
            &[
                Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
                Constraint::Bond(
                    BondId(1),
                    BondConstraintAst::ring_membership(RingScope::Size(6), 1)
                ),
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                }),
            ],
        );
        assert_eq!(
            update,
            CascadedConstraints {
                removed: vec![RemovedConstraint {
                    position: 1,
                    constraint: Constraint::Atom(AtomId(1), AtomConstraintAst::degree(3)),
                }],
                modified: vec![
                    ModifiedConstraint {
                        position: 2,
                        old: Constraint::Bond(
                            BondId(2),
                            BondConstraintAst::ring_membership(RingScope::Size(6), 1)
                        ),
                        new: Constraint::Bond(
                            BondId(1),
                            BondConstraintAst::ring_membership(RingScope::Size(6), 1)
                        ),
                    },
                    ModifiedConstraint {
                        position: 3,
                        old: Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(2)]),
                        }),
                        new: Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(1)]),
                        }),
                    },
                ],
            },
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::flattens_top_level_and_then_sorts(
        Constraints::from(vec![
            Constraint::And(vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))]),
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
        ]),
        Ok(Constraints::from(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ])),
    )]
    #[case::drops_empty_or_and_dedups(
        Constraints::from(vec![
            Constraint::Or(vec![]),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
            Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        ]),
        Ok(Constraints::from(vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))])),
    )]
    #[case::inner_contradiction_propagates(
        Constraints::from(vec![Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::lit_set(Vec::<i64>::new())))]),
        Err(Contradiction),
    )]
    fn test_constraints_canonicalize(
        #[case] input: Constraints,
        #[case] expected: Result<Constraints, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_lit(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }, false)]
    #[case::charge_sum_undetermined(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Undetermined }, true)]
    #[case::unpaired_electron_coupling_ground(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }, false)]
    #[case::unpaired_electron_coupling_undetermined(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::default() }, true)]
    #[case::bond_order_sum_lit(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(4) }, false)]
    #[case::bond_order_sum_undetermined(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Undetermined }, true)]
    #[case::connected(MoleculeConstraint::Connected { atoms: None }, false)]
    #[case::sub_pattern(MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) }, false)]
    fn test_molecule_constraint_is_vacuous(
        #[case] c: MoleculeConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_vacuous(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_sorts_and_folds(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(2), AtomId(0)]), sum: ValueAst::term(ValueTerm::Lit(1)) },
        Ok(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: ValueAst::Lit(1) }),
    )]
    #[case::unpaired_electron_coupling_sorts_and_folds(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(2), AtomId(0)]),
            unpaired_electrons: UnpairedElectronsAst { count: ValueAst::term(ValueTerm::Lit(0)), multiplicity: ValueAst::term(ValueTerm::Lit(1)) } },
        Ok(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
    )]
    #[case::bond_order_sum_sorts_and_folds(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(2), BondId(0)]), sum: ValueAst::term(ValueTerm::Lit(4)) },
        Ok(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(2)]), sum: ValueAst::Lit(4) }),
    )]
    #[case::connected_sorts(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(3), AtomId(1), AtomId(2)]) },
        Ok(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(1), AtomId(2), AtomId(3)]) }),
    )]
    #[case::charge_sum_empty_litset_contradiction(
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::lit_set(Vec::<i64>::new()) },
        Err(Contradiction),
    )]
    fn test_molecule_constraint_canonicalize(
        #[case] input: MoleculeConstraint,
        #[case] expected: Result<MoleculeConstraint, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::connected_none(MoleculeConstraint::Connected { atoms: None })]
    #[case::sub_pattern_noop(MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) })]
    fn test_molecule_constraint_canonicalize_identity(#[case] input: MoleculeConstraint) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_before_connected(
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
        MoleculeConstraint::Connected { atoms: None },
        Ordering::Less,
    )]
    #[case::connected_before_sub_pattern(
        MoleculeConstraint::Connected { atoms: None },
        MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) },
        Ordering::Less,
    )]
    #[case::sub_pattern_ignores_pattern(
        MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) },
        MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::from_entries(MoleculeEntries { atoms: vec![AtomAst::default()], bonds: vec![], ..Default::default() })) },
        Ordering::Equal,
    )]
    #[case::sub_pattern_orders_by_anchor(
        MoleculeConstraint::SubPattern { anchor: SubPatternAnchor::new(), pattern: Box::new(MoleculeAst::default()) },
        MoleculeConstraint::SubPattern {
            anchor: { let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(0), AtomId(0)); a },
            pattern: Box::new(MoleculeAst::default()),
        },
        Ordering::Less,
    )]
    fn test_molecule_constraint_cmp(
        #[case] a: MoleculeConstraint,
        #[case] b: MoleculeConstraint,
        #[case] expected: Ordering,
    ) {
        assert_eq!(a.cmp(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_shifts(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: ValueAst::Lit(1) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: ValueAst::Lit(1) }),
    )]
    #[case::charge_sum_drops_when_atom_removed(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(1), AtomId(2)]), sum: ValueAst::Lit(0) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::charge_sum_all_atoms_passes_through(
        MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }),
    )]
    #[case::unpaired_electron_coupling_shifts(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(1)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
    )]
    #[case::unpaired_electron_coupling_drops_when_atom_removed(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(1), AtomId(2)]), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::unpaired_electron_coupling_all_atoms_passes_through(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)) }),
    )]
    #[case::bond_order_sum_shifts(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(2)]), sum: ValueAst::Lit(4) },
        id_compaction(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(1)]), sum: ValueAst::Lit(4) }),
    )]
    #[case::bond_order_sum_drops_when_bond_removed(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(1)]), sum: ValueAst::Lit(2) },
        id_compaction(vec![], vec![1]),
        None,
    )]
    #[case::bond_order_sum_all_bonds_passes_through(
        MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(0) },
        id_compaction(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: None, sum: ValueAst::Lit(0) }),
    )]
    #[case::connected_shifts(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(2), AtomId(3)]) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2)]) }),
    )]
    #[case::connected_drops_when_atom_removed(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(1)]) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::connected_all_atoms_passes_through(
        MoleculeConstraint::Connected { atoms: None },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: None }),
    )]
    fn test_molecule_constraint_compact(
        #[case] c: MoleculeConstraint,
        #[case] compaction: IdCompaction,
        #[case] expected: Option<MoleculeConstraint>,
    ) {
        assert_eq!(c.compact(&compaction), expected);
    }

    #[rstest]
    fn test_sub_pattern_anchor_new() {
        let a = SubPatternAnchor::new();
        assert!(a.is_empty());
        assert!(a.atoms().is_empty());
        assert!(a.bonds().is_empty());
        assert!(a.dative_bonds().is_empty());
        assert!(a.aromatic_systems().is_empty());
        assert!(a.multicenter_bonds().is_empty());
        assert!(a.noncovalent_bonds().is_empty());
        assert!(a.stereo_atoms().is_empty());
        assert!(a.stereo_bonds().is_empty());
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_atom() {
        let mut a = SubPatternAnchor::new();
        a.push_atom(AtomId(3), AtomId(0));
        assert!(!a.is_empty());
        assert_eq!(a.atoms(), &[(AtomId(3), AtomId(0))]);
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_bond() {
        let mut a = SubPatternAnchor::new();
        a.push_bond(BondId(5), BondId(1));
        assert_eq!(a.bonds(), &[(BondId(5), BondId(1))]);
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_dative_bond() {
        let mut a = SubPatternAnchor::new();
        a.push_dative_bond(DativeBondId(4), DativeBondId(0));
        assert_eq!(a.dative_bonds(), &[(DativeBondId(4), DativeBondId(0))]);
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_aromatic_system() {
        let mut a = SubPatternAnchor::new();
        a.push_aromatic_system(AromaticSystemId(2), AromaticSystemId(0));
        assert_eq!(
            a.aromatic_systems(),
            &[(AromaticSystemId(2), AromaticSystemId(0))],
        );
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_multicenter_bond() {
        let mut a = SubPatternAnchor::new();
        a.push_multicenter_bond(MulticenterBondId(7), MulticenterBondId(2));
        assert_eq!(
            a.multicenter_bonds(),
            &[(MulticenterBondId(7), MulticenterBondId(2))],
        );
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_noncovalent_bond() {
        let mut a = SubPatternAnchor::new();
        a.push_noncovalent_bond(NoncovalentBondId(1), NoncovalentBondId(3));
        assert_eq!(
            a.noncovalent_bonds(),
            &[(NoncovalentBondId(1), NoncovalentBondId(3))],
        );
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_stereo_atom() {
        let mut a = SubPatternAnchor::new();
        a.push_stereo_atom(StereoAtomId(2), StereoAtomId(0));
        assert_eq!(a.stereo_atoms(), &[(StereoAtomId(2), StereoAtomId(0))]);
    }

    #[rstest]
    fn test_sub_pattern_anchor_push_stereo_bond() {
        let mut a = SubPatternAnchor::new();
        a.push_stereo_bond(StereoBondId(4), StereoBondId(1));
        assert_eq!(a.stereo_bonds(), &[(StereoBondId(4), StereoBondId(1))]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::shifts_target(
        { let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(3), AtomId(0)); a.push_bond(BondId(5), BondId(1)); a },
        id_compaction(vec![1], vec![2]),
        Some({ let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(2), AtomId(0)); a.push_bond(BondId(4), BondId(1)); a }),
    )]
    #[case::drops_on_removed_target_atom(
        { let mut a = SubPatternAnchor::new(); a.push_atom(AtomId(2), AtomId(0)); a },
        id_compaction(vec![2], vec![]),
        None,
    )]
    #[case::dative_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_dative_bond(DativeBondId(1), DativeBondId(0)); a },
        relation_compaction(vec![1], vec![], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::aromatic_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_aromatic_system(AromaticSystemId(2), AromaticSystemId(0)); a },
        relation_compaction(vec![], vec![2], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::multicenter_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_multicenter_bond(MulticenterBondId(3), MulticenterBondId(0)); a },
        relation_compaction(vec![], vec![], vec![3], vec![], vec![], vec![]),
        None,
    )]
    #[case::noncovalent_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_noncovalent_bond(NoncovalentBondId(4), NoncovalentBondId(0)); a },
        relation_compaction(vec![], vec![], vec![], vec![4], vec![], vec![]),
        None,
    )]
    #[case::stereo_atom_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_stereo_atom(StereoAtomId(2), StereoAtomId(0)); a },
        relation_compaction(vec![], vec![], vec![], vec![], vec![2], vec![]),
        None,
    )]
    #[case::stereo_bond_dropped(
        { let mut a = SubPatternAnchor::new(); a.push_stereo_bond(StereoBondId(3), StereoBondId(0)); a },
        relation_compaction(vec![], vec![], vec![], vec![], vec![], vec![3]),
        None,
    )]
    fn test_sub_pattern_anchor_compact(
        #[case] anchor: SubPatternAnchor,
        #[case] compaction: IdCompaction,
        #[case] expected: Option<SubPatternAnchor>,
    ) {
        assert_eq!(anchor.compact(&compaction), expected);
    }

    #[rstest]
    fn test_sub_pattern_anchor_compact_relations_shift() {
        let mut a = SubPatternAnchor::new();
        a.push_dative_bond(DativeBondId(2), DativeBondId(0));
        a.push_aromatic_system(AromaticSystemId(3), AromaticSystemId(1));
        a.push_multicenter_bond(MulticenterBondId(4), MulticenterBondId(2));
        a.push_noncovalent_bond(NoncovalentBondId(5), NoncovalentBondId(3));
        a.push_stereo_atom(StereoAtomId(6), StereoAtomId(4));
        a.push_stereo_bond(StereoBondId(7), StereoBondId(5));

        let compaction = relation_compaction(vec![0], vec![1], vec![0], vec![2], vec![1], vec![3]);
        let mapped = a.compact(&compaction).expect("all targets survive");

        assert_eq!(mapped.dative_bonds(), &[(DativeBondId(1), DativeBondId(0))],);
        assert_eq!(
            mapped.aromatic_systems(),
            &[(AromaticSystemId(2), AromaticSystemId(1))],
        );
        assert_eq!(
            mapped.multicenter_bonds(),
            &[(MulticenterBondId(3), MulticenterBondId(2))],
        );
        assert_eq!(
            mapped.noncovalent_bonds(),
            &[(NoncovalentBondId(4), NoncovalentBondId(3))],
        );
        assert_eq!(mapped.stereo_atoms(), &[(StereoAtomId(5), StereoAtomId(4))],);
        assert_eq!(mapped.stereo_bonds(), &[(StereoBondId(6), StereoBondId(5))],);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(
        vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))],
        vec![Constraint::Bond(BondId(0), BondConstraintAst::Aromatic(BooleanAst::Lit(true)))],
    )]
    #[case::preserves_order_and_duplicates(
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(3)),
        ],
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(3)),
        ],
    )]
    #[case::empty(vec![], vec![])]
    fn test_constraints_from_iter(
        #[case] input: Vec<Constraint>,
        #[case] expected: Vec<Constraint>,
    ) {
        let cs = Constraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_constraints_into_iter() {
        let cs = Constraints::from_iter([
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(
                BondId(0),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            ),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintAst::Aromatic(BooleanAst::Lit(true))
                ),
            ],
        );
    }

    #[rstest]
    fn test_constraints_from_constraint() {
        let cs: Constraints = Constraint::Bond(
            BondId(0),
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
        )
        .into();
        assert_eq!(
            cs.as_slice(),
            &[Constraint::Bond(
                BondId(0),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true))
            )],
        );
    }

    #[rstest]
    fn test_constraints_from_vec() {
        let cs: Constraints = vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
            Constraint::Bond(
                BondId(0),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            ),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                Constraint::Atom(AtomId(0), AtomConstraintAst::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintAst::Aromatic(BooleanAst::Lit(true))
                ),
            ],
        );
    }
}
