//! Cross-entity relational constraints.
//!
//! A `RelationalConstraint` relates one entity kind to one or more other entity kinds
//! by reference. Every variant carries at least two indices (the outer entity
//! plus one or more inner atom/bond refs) or one index plus a delegated
//! predicate over a role slot.
//!
//! Relational constraints live **only** at molecule scope — as entries in
//! `MoleculeAst::constraints` (via `Constraint::Relational(...)`) or inside
//! `And`/`Or`/`Not` combinators. They cannot be inline on the entity AST:
//! the per-entity constraint containers are narrowed to value-only
//! variants so ref-bearing constraints are unrepresentable inline.
//!
//! Two sub-patterns share the enum:
//! - **Role identity / set membership**: `Donor`, `Donors`, `ContainsAllDonors`,
//!   `Acceptor`, `Parallels`, `Ends`, `Atoms`, `Contains`, `ContainsAll` —
//!   constrain an atom/bond identity to a role or set membership.
//! - **Role predicate**: `AllDonors`, `AnyDonor`, `AcceptorSatisfies`,
//!   `AllAtoms`, `AnyAtom`, `EndsSatisfy` — delegate an `AtomConstraint` to a
//!   role slot, quantified over the matching participants.

use super::super::error::Contradiction;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::remap::IdRemapping;
use super::super::traits::Canonicalize;
use super::atom::AtomConstraint;

/// Cross-entity constraint relating one overlay entity (dative bond, aromatic
/// system, multicenter bond, noncovalent bond, stereo atom, stereo bond) to
/// atoms, bonds, or atom predicates by reference. Stereo variants constrain the
/// site identity and the atom-kind ligands. Lives only at molecule scope (in
/// [`Constraint::Relational`](super::Constraint::Relational) or inside
/// `And`/`Or`/`Not`); cannot appear inline on an entity.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RelationalConstraint {
    /// Dative bond `bond` has exactly `atoms` as its donors (as a set).
    DativeBondDonors {
        bond: DativeBondId,
        atoms: Vec<AtomId>,
    },
    /// `atom` is one of the donors of dative bond `bond`.
    DativeBondDonor { bond: DativeBondId, atom: AtomId },
    /// Every atom in `atoms` is a donor of dative bond `bond` (set inclusion,
    /// not equality).
    DativeBondContainsAllDonors {
        bond: DativeBondId,
        atoms: Vec<AtomId>,
    },
    /// Every donor of dative bond `bond` satisfies `predicate`.
    DativeBondAllDonors {
        bond: DativeBondId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one donor of dative bond `bond` satisfies `predicate`.
    DativeBondAnyDonor {
        bond: DativeBondId,
        predicate: Box<AtomConstraint>,
    },
    /// The acceptor of dative bond `bond` is `atom`.
    DativeBondAcceptor { bond: DativeBondId, atom: AtomId },
    /// The acceptor of dative bond `bond` satisfies `predicate`.
    DativeBondAcceptorSatisfies {
        bond: DativeBondId,
        predicate: Box<AtomConstraint>,
    },
    /// Dative bond `dative` is parallel to a localized bond `parallel`
    /// (same atom pair).
    DativeBondParallels {
        dative: DativeBondId,
        parallel: BondId,
    },

    /// Aromatic system `system` consists of exactly `atoms` (as a set).
    AromaticSystemAtoms {
        system: AromaticSystemId,
        atoms: Vec<AtomId>,
    },
    /// `atom` is one of the atoms participating in aromatic system `system`.
    AromaticSystemContains {
        system: AromaticSystemId,
        atom: AtomId,
    },
    /// Every atom in `atoms` participates in aromatic system `system`
    /// (set inclusion, not equality).
    AromaticSystemContainsAll {
        system: AromaticSystemId,
        atoms: Vec<AtomId>,
    },
    /// Every atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAllAtoms {
        system: AromaticSystemId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAnyAtom {
        system: AromaticSystemId,
        predicate: Box<AtomConstraint>,
    },

    /// Multicenter bond `bond` consists of exactly `atoms` (as a set).
    MulticenterBondAtoms {
        bond: MulticenterBondId,
        atoms: Vec<AtomId>,
    },
    /// `atom` is one of the participants in multicenter bond `bond`.
    MulticenterBondContains {
        bond: MulticenterBondId,
        atom: AtomId,
    },
    /// Every atom in `atoms` participates in multicenter bond `bond`.
    MulticenterBondContainsAll {
        bond: MulticenterBondId,
        atoms: Vec<AtomId>,
    },
    /// Every participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAllAtoms {
        bond: MulticenterBondId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAnyAtom {
        bond: MulticenterBondId,
        predicate: Box<AtomConstraint>,
    },

    /// Noncovalent bond `bond` connects exactly the pair `atoms` (unordered).
    NoncovalentBondEnds {
        bond: NoncovalentBondId,
        atoms: [AtomId; 2],
    },
    /// `atom` is one of the two endpoints of noncovalent bond `bond`.
    NoncovalentBondContains {
        bond: NoncovalentBondId,
        atom: AtomId,
    },
    /// The two endpoints of noncovalent bond `bond` satisfy `predicates[0]`
    /// and `predicates[1]` respectively. Order is not symmetric: the bond
    /// stores its endpoints as an unordered pair, but each predicate is
    /// associated with one specific slot.
    NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondId,
        predicates: [Box<AtomConstraint>; 2],
    },

    /// The site of stereo atom `stereo_atom` is `atom`.
    StereoAtomSite {
        stereo_atom: StereoAtomId,
        atom: AtomId,
    },
    /// `atom` is one of the atom-kind ligands of stereo atom `stereo_atom`.
    StereoAtomContains {
        stereo_atom: StereoAtomId,
        atom: AtomId,
    },
    /// The atom-kind ligands of stereo atom `stereo_atom` are exactly `atoms`
    /// (as a set; ligand order/configuration is the inline coset's concern).
    StereoAtomLigands {
        stereo_atom: StereoAtomId,
        atoms: Vec<AtomId>,
    },
    /// Every atom-kind ligand of stereo atom `stereo_atom` satisfies `predicate`.
    StereoAtomAllLigands {
        stereo_atom: StereoAtomId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one atom-kind ligand of stereo atom `stereo_atom` satisfies
    /// `predicate`.
    StereoAtomAnyLigand {
        stereo_atom: StereoAtomId,
        predicate: Box<AtomConstraint>,
    },

    /// The site of stereo bond `stereo_bond` is `bond`.
    StereoBondSite {
        stereo_bond: StereoBondId,
        bond: BondId,
    },
    /// `atom` is one of the atom-kind ligands of stereo bond `stereo_bond`.
    StereoBondContains {
        stereo_bond: StereoBondId,
        atom: AtomId,
    },
    /// The atom-kind ligands of stereo bond `stereo_bond` are exactly `atoms`
    /// (as a set).
    StereoBondLigands {
        stereo_bond: StereoBondId,
        atoms: Vec<AtomId>,
    },
    /// Every atom-kind ligand of stereo bond `stereo_bond` satisfies `predicate`.
    StereoBondAllLigands {
        stereo_bond: StereoBondId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one atom-kind ligand of stereo bond `stereo_bond` satisfies
    /// `predicate`.
    StereoBondAnyLigand {
        stereo_bond: StereoBondId,
        predicate: Box<AtomConstraint>,
    },
}

impl RelationalConstraint {
    /// Remap all indices this constraint carries. Returns `None` if any
    /// referenced entity has been removed by the remapping.
    pub fn remap(self, remap: &IdRemapping) -> Option<Self> {
        Some(match self {
            Self::DativeBondDonors { bond, atoms } => {
                let bond = remap.dative_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::DativeBondDonors {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::DativeBondDonor { bond, atom } => Self::DativeBondDonor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::DativeBondContainsAllDonors { bond, atoms } => {
                let bond = remap.dative_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::DativeBondContainsAllDonors {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::DativeBondAllDonors { bond, predicate } => Self::DativeBondAllDonors {
                bond: remap.dative_bond(bond)?,
                predicate,
            },
            Self::DativeBondAnyDonor { bond, predicate } => Self::DativeBondAnyDonor {
                bond: remap.dative_bond(bond)?,
                predicate,
            },
            Self::DativeBondAcceptor { bond, atom } => Self::DativeBondAcceptor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies {
                    bond: remap.dative_bond(bond)?,
                    predicate,
                }
            }
            Self::DativeBondParallels { dative, parallel } => Self::DativeBondParallels {
                dative: remap.dative_bond(dative)?,
                parallel: remap.bond(parallel)?,
            },
            Self::AromaticSystemAtoms { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::AromaticSystemAtoms {
                    system,
                    atoms: atoms?,
                }
            }
            Self::AromaticSystemContains { system, atom } => Self::AromaticSystemContains {
                system: remap.aromatic_system(system)?,
                atom: remap.atom(atom)?,
            },
            Self::AromaticSystemContainsAll { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::AromaticSystemContainsAll {
                    system,
                    atoms: atoms?,
                }
            }
            Self::AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            Self::AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            Self::MulticenterBondAtoms { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::MulticenterBondAtoms {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::MulticenterBondContains { bond, atom } => Self::MulticenterBondContains {
                bond: remap.multicenter_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::MulticenterBondContainsAll { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::MulticenterBondContainsAll {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            Self::MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            Self::NoncovalentBondEnds { bond, atoms } => {
                let bond = remap.noncovalent_bond(bond)?;
                let [a, b] = atoms;
                Self::NoncovalentBondEnds {
                    bond,
                    atoms: [remap.atom(a)?, remap.atom(b)?],
                }
            }
            Self::NoncovalentBondContains { bond, atom } => Self::NoncovalentBondContains {
                bond: remap.noncovalent_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::NoncovalentBondEndsSatisfy { bond, predicates } => {
                Self::NoncovalentBondEndsSatisfy {
                    bond: remap.noncovalent_bond(bond)?,
                    predicates,
                }
            }
            Self::StereoAtomSite { stereo_atom, atom } => Self::StereoAtomSite {
                stereo_atom: remap.stereo_atom(stereo_atom)?,
                atom: remap.atom(atom)?,
            },
            Self::StereoAtomContains { stereo_atom, atom } => Self::StereoAtomContains {
                stereo_atom: remap.stereo_atom(stereo_atom)?,
                atom: remap.atom(atom)?,
            },
            Self::StereoAtomLigands { stereo_atom, atoms } => {
                let stereo_atom = remap.stereo_atom(stereo_atom)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::StereoAtomLigands {
                    stereo_atom,
                    atoms: atoms?,
                }
            }
            Self::StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAllLigands {
                stereo_atom: remap.stereo_atom(stereo_atom)?,
                predicate,
            },
            Self::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand {
                stereo_atom: remap.stereo_atom(stereo_atom)?,
                predicate,
            },
            Self::StereoBondSite { stereo_bond, bond } => Self::StereoBondSite {
                stereo_bond: remap.stereo_bond(stereo_bond)?,
                bond: remap.bond(bond)?,
            },
            Self::StereoBondContains { stereo_bond, atom } => Self::StereoBondContains {
                stereo_bond: remap.stereo_bond(stereo_bond)?,
                atom: remap.atom(atom)?,
            },
            Self::StereoBondLigands { stereo_bond, atoms } => {
                let stereo_bond = remap.stereo_bond(stereo_bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::StereoBondLigands {
                    stereo_bond,
                    atoms: atoms?,
                }
            }
            Self::StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => Self::StereoBondAllLigands {
                stereo_bond: remap.stereo_bond(stereo_bond)?,
                predicate,
            },
            Self::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand {
                stereo_bond: remap.stereo_bond(stereo_bond)?,
                predicate,
            },
        })
    }
}

impl Canonicalize for RelationalConstraint {
    /// Canonicalize any inner `AtomConstraint` predicate. Refs (`bond`, `atom`,
    /// `system`, `parallel`, atom sets) are unchanged.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::DativeBondAllDonors { bond, predicate } => Self::DativeBondAllDonors {
                bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::DativeBondAnyDonor { bond, predicate } => Self::DativeBondAnyDonor {
                bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies {
                    bond,
                    predicate: Box::new((*predicate).canonicalize()?),
                }
            }
            Self::AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::NoncovalentBondEndsSatisfy { bond, predicates } => {
                let [a, b] = predicates;
                Self::NoncovalentBondEndsSatisfy {
                    bond,
                    predicates: [
                        Box::new((*a).canonicalize()?),
                        Box::new((*b).canonicalize()?),
                    ],
                }
            }
            Self::StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAllLigands {
                stereo_atom,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand {
                stereo_atom,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => Self::StereoBondAllLigands {
                stereo_bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            Self::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand {
                stereo_bond,
                predicate: Box::new((*predicate).canonicalize()?),
            },
            other @ (Self::DativeBondDonors { .. }
            | Self::DativeBondDonor { .. }
            | Self::DativeBondContainsAllDonors { .. }
            | Self::DativeBondAcceptor { .. }
            | Self::DativeBondParallels { .. }
            | Self::AromaticSystemAtoms { .. }
            | Self::AromaticSystemContains { .. }
            | Self::AromaticSystemContainsAll { .. }
            | Self::MulticenterBondAtoms { .. }
            | Self::MulticenterBondContains { .. }
            | Self::MulticenterBondContainsAll { .. }
            | Self::NoncovalentBondEnds { .. }
            | Self::NoncovalentBondContains { .. }
            | Self::StereoAtomSite { .. }
            | Self::StereoAtomContains { .. }
            | Self::StereoAtomLigands { .. }
            | Self::StereoBondSite { .. }
            | Self::StereoBondContains { .. }
            | Self::StereoBondLigands { .. }) => other,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{RelationId, Remapping};

    use super::*;
    use crate::ast::value::ValueAst;

    #[allow(clippy::too_many_arguments)]
    fn remapping(
        removed_nodes: Vec<u32>,
        removed_edges: Vec<u32>,
        removed_dative: Vec<u32>,
        removed_aromatic: Vec<u32>,
        removed_multicenter: Vec<u32>,
        removed_noncovalent: Vec<u32>,
        removed_stereo_atom: Vec<u32>,
        removed_stereo_bond: Vec<u32>,
    ) -> IdRemapping {
        let rel = |v: Vec<u32>| v.into_iter().map(RelationId).collect::<Vec<_>>();
        IdRemapping::new(
            Remapping::new(removed_nodes, removed_edges),
            rel(removed_dative),
            rel(removed_aromatic),
            rel(removed_multicenter),
            rel(removed_noncovalent),
            rel(removed_stereo_atom),
            rel(removed_stereo_bond),
        )
    }

    fn val_pred() -> Box<AtomConstraint> {
        Box::new(AtomConstraint::Valence(ValueAst::Lit(4)))
    }

    /// Drop atom 1; drop dative 0; preserve other entities. Indices above
    /// the removed slot shift down by one.
    fn one_atom_one_dative() -> IdRemapping {
        remapping(
            vec![1],
            vec![],
            vec![0],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
    }

    /// Drop bond 0; preserve other entities.
    fn drop_bond0() -> IdRemapping {
        remapping(
            vec![],
            vec![0],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
    }

    /// Drop aromatic system 0.
    fn drop_aromatic0() -> IdRemapping {
        remapping(
            vec![],
            vec![],
            vec![],
            vec![0],
            vec![],
            vec![],
            vec![],
            vec![],
        )
    }

    /// Drop multicenter bond 0.
    fn drop_multicenter0() -> IdRemapping {
        remapping(
            vec![],
            vec![],
            vec![],
            vec![],
            vec![0],
            vec![],
            vec![],
            vec![],
        )
    }

    /// Drop noncovalent bond 0.
    fn drop_noncovalent0() -> IdRemapping {
        remapping(
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![0],
            vec![],
            vec![],
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor_shifts(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(1) }))]
    #[case::dative_donor_drops_when_atom_removed(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(1) },
        one_atom_one_dative(), None)]
    #[case::dative_donor_drops_when_bond_removed(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(2) },
        one_atom_one_dative(), None)]
    #[case::dative_acceptor_shifts(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(1), atom: AtomId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(0), atom: AtomId(1) }))]
    #[case::dative_parallels_shifts(RelationalConstraint::DativeBondParallels { dative: DativeBondId(1), parallel: BondId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(2) }))]
    #[case::dative_parallels_drops_when_bond_removed(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(0) },
        drop_bond0(), None)]
    #[case::dative_all_donors_shifts(RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(1), predicate: val_pred() },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(0), predicate: val_pred() }))]
    #[case::dative_donors_shifts(RelationalConstraint::DativeBondDonors { bond: DativeBondId(1), atoms: vec![AtomId(2)] },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondDonors { bond: DativeBondId(0), atoms: vec![AtomId(1)] }))]
    #[case::dative_acceptor_satisfies_shifts(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondId(1), predicate: val_pred() },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondId(0), predicate: val_pred() }))]
    #[case::aromatic_atoms_shifts(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(1), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![0], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::aromatic_atoms_drops_when_atom_removed(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0), atoms: vec![AtomId(1), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::aromatic_atoms_drops_when_system_removed(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0), atoms: vec![AtomId(0)] },
        drop_aromatic0(), None)]
    #[case::aromatic_contains_shifts(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![0], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(1) }))]
    #[case::aromatic_contains_all_shifts(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::aromatic_contains_all_drops_when_atom_removed(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::aromatic_all_atoms_shifts(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![0], vec![], vec![], vec![], vec![]),
        Some(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(0), predicate: val_pred() }))]
    #[case::aromatic_any_atom_shifts(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![0], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemId(0), predicate: val_pred() }))]
    #[case::multicenter_atoms_shifts(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(1), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![0], vec![], vec![], vec![]), Some(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::multicenter_atoms_drops_when_atom_removed(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::multicenter_atoms_drops_when_bond_removed(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0), atoms: vec![AtomId(0)] },
        drop_multicenter0(), None)]
    #[case::multicenter_contains_shifts(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![0], vec![], vec![], vec![]), Some(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(1) }))]
    #[case::multicenter_contains_all_shifts(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::multicenter_contains_all_drops_when_atom_removed(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::multicenter_all_atoms_shifts(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![0], vec![], vec![], vec![]),
        Some(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondId(0), predicate: val_pred() }))]
    #[case::multicenter_any_atom_shifts(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![0], vec![], vec![], vec![]),
        Some(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(0), predicate: val_pred() }))]
    #[case::noncovalent_ends_shifts(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(1), atoms: [AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![0], vec![], vec![]), Some(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)] }))]
    #[case::noncovalent_ends_drops_when_atom_removed(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::noncovalent_ends_drops_when_bond_removed(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)] },
        drop_noncovalent0(), None)]
    #[case::noncovalent_contains_shifts(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![0], vec![], vec![]),
        Some(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(1) }))]
    #[case::noncovalent_ends_satisfy_shifts(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(1), predicates: [val_pred(), val_pred()] },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![0], vec![], vec![]), Some(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
        predicates: [val_pred(), val_pred()] }))]
    #[case::stereo_atom_site_shifts(RelationalConstraint::StereoAtomSite { stereo_atom: StereoAtomId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::StereoAtomSite { stereo_atom: StereoAtomId(0), atom: AtomId(1) }))]
    #[case::stereo_atom_site_drops_when_stereo_removed(RelationalConstraint::StereoAtomSite { stereo_atom: StereoAtomId(0), atom: AtomId(2) },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![0], vec![]), None)]
    #[case::stereo_atom_site_drops_when_atom_removed(RelationalConstraint::StereoAtomSite { stereo_atom: StereoAtomId(0), atom: AtomId(1) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::stereo_atom_contains_shifts(RelationalConstraint::StereoAtomContains { stereo_atom: StereoAtomId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::StereoAtomContains { stereo_atom: StereoAtomId(0), atom: AtomId(1) }))]
    #[case::stereo_atom_ligands_shifts(RelationalConstraint::StereoAtomLigands { stereo_atom: StereoAtomId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]),
        Some(RelationalConstraint::StereoAtomLigands { stereo_atom: StereoAtomId(0), atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::stereo_atom_ligands_drops_when_atom_removed(RelationalConstraint::StereoAtomLigands { stereo_atom: StereoAtomId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::stereo_atom_all_ligands_shifts(RelationalConstraint::StereoAtomAllLigands { stereo_atom: StereoAtomId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::StereoAtomAllLigands { stereo_atom: StereoAtomId(0), predicate: val_pred() }))]
    #[case::stereo_atom_any_ligand_shifts(RelationalConstraint::StereoAtomAnyLigand { stereo_atom: StereoAtomId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::StereoAtomAnyLigand { stereo_atom: StereoAtomId(0), predicate: val_pred() }))]
    #[case::stereo_bond_site_shifts(RelationalConstraint::StereoBondSite { stereo_bond: StereoBondId(1), bond: BondId(2) },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![], vec![0]),
        Some(RelationalConstraint::StereoBondSite { stereo_bond: StereoBondId(0), bond: BondId(2) }))]
    #[case::stereo_bond_site_drops_when_bond_removed(RelationalConstraint::StereoBondSite { stereo_bond: StereoBondId(0), bond: BondId(0) },
        remapping(vec![], vec![0], vec![], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::stereo_bond_contains_shifts(RelationalConstraint::StereoBondContains { stereo_bond: StereoBondId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![0]),
        Some(RelationalConstraint::StereoBondContains { stereo_bond: StereoBondId(0), atom: AtomId(1) }))]
    #[case::stereo_bond_ligands_shifts(RelationalConstraint::StereoBondLigands { stereo_bond: StereoBondId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![], vec![], vec![]),
        Some(RelationalConstraint::StereoBondLigands { stereo_bond: StereoBondId(0), atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::stereo_bond_all_ligands_shifts(RelationalConstraint::StereoBondAllLigands { stereo_bond: StereoBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![], vec![0]),
        Some(RelationalConstraint::StereoBondAllLigands { stereo_bond: StereoBondId(0), predicate: val_pred() }))]
    #[case::stereo_bond_any_ligand_shifts(RelationalConstraint::StereoBondAnyLigand { stereo_bond: StereoBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![], vec![], vec![0]),
        Some(RelationalConstraint::StereoBondAnyLigand { stereo_bond: StereoBondId(0), predicate: val_pred() }))]
    fn test_relational_constraint_remap(
        #[case] input: RelationalConstraint,
        #[case] remap: IdRemapping,
        #[case] expected: Option<RelationalConstraint>,
    ) {
        assert_eq!(input.remap(&remap), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::predicate_litset_singleton(
        RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::lit_set([4]))) },
        Ok(RelationalConstraint::DativeBondAllDonors { bond: DativeBondId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::Lit(4))) }))]
    #[case::both_ends_canonicalize(
        RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
            predicates: [Box::new(AtomConstraint::Valence(ValueAst::lit_set([4]))), Box::new(AtomConstraint::Degree(ValueAst::lit_set([2])))] },
        Ok(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
            predicates: [Box::new(AtomConstraint::Valence(ValueAst::Lit(4))), Box::new(AtomConstraint::Degree(ValueAst::Lit(2)))] }))]
    #[case::predicate_empty_litset_contradiction(
        RelationalConstraint::StereoAtomAllLigands { stereo_atom: StereoAtomId(0), predicate: Box::new(AtomConstraint::Valence(ValueAst::lit_set(Vec::<i64>::new()))) },
        Err(Contradiction))]
    fn test_relational_constraint_canonicalize(
        #[case] input: RelationalConstraint,
        #[case] expected: Result<RelationalConstraint, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(2) })]
    #[case::dative_donors(RelationalConstraint::DativeBondDonors { bond: DativeBondId(0), atoms: vec![AtomId(1), AtomId(2)] })]
    #[case::aromatic_contains(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(1) })]
    fn test_relational_constraint_canonicalize_identity(#[case] input: RelationalConstraint) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }
}
