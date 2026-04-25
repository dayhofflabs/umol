//! Cross-entity relational constraints.
//!
//! A `RelationalConstraint` relates one DAMN entity (dative bond, aromatic
//! system, multicenter bond, noncovalent bond) to one or more other entities
//! by reference. Every variant carries at least two indices (the outer entity
//! plus one or more inner atom/bond refs) or one index plus a delegated
//! predicate over a role slot.
//!
//! Relational constraints live **only** at molecule scope — as entries in
//! `MoleculeAst::constraints` (via `Constraint::Relational(...)`) or inside
//! `And`/`Or`/`Not` combinators. They cannot be inline on the entity AST:
//! the per-entity `XxxConstraints` containers are narrowed to value-only
//! variants so ref-bearing constraints are unrepresentable inline.
//!
//! Two sub-patterns share the enum:
//! - **Role identity**: `Donor`, `Acceptor`, `Parallels`, `Ends`, `Atoms`,
//!   `Contains`, `ContainsAll` — pin an atom/bond identity to a role or set
//!   membership.
//! - **Role predicate**: `DonorSatisfies`, `AcceptorSatisfies`, `AllAtoms`,
//!   `AnyAtom`, `EndsSatisfy` — delegate an `AtomConstraint` to a role slot,
//!   quantified over the matching participants.

use super::super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::super::remap::IdxRemapping;
use super::atom::AtomConstraint;

/// Cross-entity constraint relating one DAMN entity (dative bond, aromatic
/// system, multicenter bond, noncovalent bond) to atoms, bonds, or atom
/// predicates by reference. Lives only at molecule scope (in
/// [`Constraint::Relational`](super::Constraint::Relational) or inside
/// `And`/`Or`/`Not`); cannot appear inline on an entity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RelationalConstraint {
    // region: Dative bond
    /// The donor end of dative bond `bond` is `atom`.
    DativeBondDonor {
        bond: DativeBondIdx,
        atom: AtomIdx,
    },
    /// The acceptor end of dative bond `bond` is `atom`.
    DativeBondAcceptor {
        bond: DativeBondIdx,
        atom: AtomIdx,
    },
    /// Dative bond `dative` is parallel to a covalent bond `parallel`
    /// (same atom pair).
    DativeBondParallels {
        dative: DativeBondIdx,
        parallel: BondIdx,
    },
    /// The donor end of dative bond `bond` satisfies `predicate`.
    DativeBondDonorSatisfies {
        bond: DativeBondIdx,
        predicate: Box<AtomConstraint>,
    },
    /// The acceptor end of dative bond `bond` satisfies `predicate`.
    DativeBondAcceptorSatisfies {
        bond: DativeBondIdx,
        predicate: Box<AtomConstraint>,
    },

    // endregion: Dative bond

    // region: Aromatic system
    /// Aromatic system `system` consists of exactly `atoms` (as a set).
    AromaticSystemAtoms {
        system: AromaticSystemIdx,
        atoms: Vec<AtomIdx>,
    },
    /// `atom` is one of the atoms participating in aromatic system `system`.
    AromaticSystemContains {
        system: AromaticSystemIdx,
        atom: AtomIdx,
    },
    /// Every atom in `atoms` participates in aromatic system `system`
    /// (set inclusion, not equality).
    AromaticSystemContainsAll {
        system: AromaticSystemIdx,
        atoms: Vec<AtomIdx>,
    },
    /// Every atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAllAtoms {
        system: AromaticSystemIdx,
        predicate: Box<AtomConstraint>,
    },
    /// At least one atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAnyAtom {
        system: AromaticSystemIdx,
        predicate: Box<AtomConstraint>,
    },

    // endregion: Aromatic system

    // region: Multicenter bond
    /// Multicenter bond `bond` consists of exactly `atoms` (as a set).
    MulticenterBondAtoms {
        bond: MulticenterBondIdx,
        atoms: Vec<AtomIdx>,
    },
    /// `atom` is one of the participants in multicenter bond `bond`.
    MulticenterBondContains {
        bond: MulticenterBondIdx,
        atom: AtomIdx,
    },
    /// Every atom in `atoms` participates in multicenter bond `bond`.
    MulticenterBondContainsAll {
        bond: MulticenterBondIdx,
        atoms: Vec<AtomIdx>,
    },
    /// Every participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAllAtoms {
        bond: MulticenterBondIdx,
        predicate: Box<AtomConstraint>,
    },
    /// At least one participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAnyAtom {
        bond: MulticenterBondIdx,
        predicate: Box<AtomConstraint>,
    },

    // endregion: Multicenter bond

    // region: Noncovalent bond
    /// Noncovalent bond `bond` connects exactly the pair `atoms` (unordered).
    NoncovalentBondEnds {
        bond: NoncovalentBondIdx,
        atoms: [AtomIdx; 2],
    },
    /// `atom` is one of the two endpoints of noncovalent bond `bond`.
    NoncovalentBondContains {
        bond: NoncovalentBondIdx,
        atom: AtomIdx,
    },
    /// The two endpoints of noncovalent bond `bond` satisfy `predicates[0]`
    /// and `predicates[1]` respectively. Order is not symmetric: the bond
    /// stores its endpoints as an unordered pair, but each predicate is
    /// associated with one specific slot.
    NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondIdx,
        predicates: [Box<AtomConstraint>; 2],
    },
    // endregion: Noncovalent bond
}

impl RelationalConstraint {
    /// Remap all indices this constraint carries. Returns `None` if any
    /// referenced entity has been removed by the remapping.
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        use RelationalConstraint::*;
        Some(match self {
            DativeBondDonor { bond, atom } => DativeBondDonor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            DativeBondAcceptor { bond, atom } => DativeBondAcceptor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            DativeBondParallels { dative, parallel } => DativeBondParallels {
                dative: remap.dative_bond(dative)?,
                parallel: remap.bond(parallel)?,
            },
            DativeBondDonorSatisfies { bond, predicate } => DativeBondDonorSatisfies {
                bond: remap.dative_bond(bond)?,
                predicate,
            },
            DativeBondAcceptorSatisfies { bond, predicate } => DativeBondAcceptorSatisfies {
                bond: remap.dative_bond(bond)?,
                predicate,
            },
            AromaticSystemAtoms { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                AromaticSystemAtoms {
                    system,
                    atoms: atoms?,
                }
            }
            AromaticSystemContains { system, atom } => AromaticSystemContains {
                system: remap.aromatic_system(system)?,
                atom: remap.atom(atom)?,
            },
            AromaticSystemContainsAll { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                AromaticSystemContainsAll {
                    system,
                    atoms: atoms?,
                }
            }
            AromaticSystemAllAtoms { system, predicate } => AromaticSystemAllAtoms {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            AromaticSystemAnyAtom { system, predicate } => AromaticSystemAnyAtom {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            MulticenterBondAtoms { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                MulticenterBondAtoms {
                    bond,
                    atoms: atoms?,
                }
            }
            MulticenterBondContains { bond, atom } => MulticenterBondContains {
                bond: remap.multicenter_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            MulticenterBondContainsAll { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                MulticenterBondContainsAll {
                    bond,
                    atoms: atoms?,
                }
            }
            MulticenterBondAllAtoms { bond, predicate } => MulticenterBondAllAtoms {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            MulticenterBondAnyAtom { bond, predicate } => MulticenterBondAnyAtom {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            NoncovalentBondEnds { bond, atoms } => {
                let bond = remap.noncovalent_bond(bond)?;
                let [a, b] = atoms;
                NoncovalentBondEnds {
                    bond,
                    atoms: [remap.atom(a)?, remap.atom(b)?],
                }
            }
            NoncovalentBondContains { bond, atom } => NoncovalentBondContains {
                bond: remap.noncovalent_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            NoncovalentBondEndsSatisfy { bond, predicates } => NoncovalentBondEndsSatisfy {
                bond: remap.noncovalent_bond(bond)?,
                predicates,
            },
        })
    }
}
