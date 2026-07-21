//! Reaction AST: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations live in `reaction_span.rs`).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::iter::from_fn;

use umol_graph_core::{Correspondence, NodeId, SubgraphIsomorphismAlgorithm};
use umol_perm::Permutation;

use super::aromatic::{AromaticSystemAst, AromaticSystemUpdate};
use super::atom::{AtomAst, AtomUpdate};
use super::bond::{BondAst, BondUpdate};
use super::correspondence::MoleculeCorrespondence;
use super::dative::{DativeBondAst, DativeBondUpdate};
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::edit::{
    AddBond, AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle,
    BondFieldChange, BondHandle, DativeBondFieldChange, DativeBondHandle, Edit,
    MulticenterBondFieldChange, MulticenterBondHandle, NoncovalentBondFieldChange,
    NoncovalentBondHandle, StereoAtomFieldChange, StereoAtomHandle, StereoAtomRemoval,
    StereoBondFieldChange, StereoBondHandle, StereoBondRemoval,
};
use super::error::{ApplyError, ApplyPreconditionError, Contradiction};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::MoleculeAst;
#[cfg(test)]
use super::molecule::MoleculeParts;
use super::multicenter::{MulticenterBondAst, MulticenterBondUpdate};
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondUpdate};
use super::reaction_derivation::ReactionDerivation;
use super::remap::IdRemapping;
use super::stereo::StereoConfigurationAst;
use super::substructure::SubstructureMatchAlgorithm;
use super::traits::Canonicalize;
use super::validate::{
    DpoValidator, EntityStructureValidator, ReactionIntegrityContradiction,
    ReactionIntegrityValidator,
};

/// A reaction as one full molecule state (`lhs`) plus one resolved delta (`deltas`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionAst {
    pub lhs: MoleculeAst,
    pub deltas: Deltas,
}

impl ReactionAst {
    pub fn new(lhs: MoleculeAst, deltas: Deltas) -> Self {
        Self { lhs, deltas }
    }

    /// The reaction transforming `lhs` into `rhs` under the atom correspondence `atom`: induce the
    /// full per-entity correspondence, diff the two sides into deltas, and pair them with `lhs`. The
    /// inverse of reading a reaction's two sides back off its span.
    pub fn from_sides(lhs: MoleculeAst, rhs: MoleculeAst, atom: Correspondence<NodeId>) -> Self {
        let correspondence = MoleculeCorrespondence::induce(&lhs, &rhs, atom);
        let deltas = lhs.difference_to(&rhs, &correspondence);
        Self::new(lhs, deltas)
    }

    /// Apply the reaction at one match of `lhs` into `host` — the injective pattern→host
    /// `correspondence` — producing the derivation `lhs ⇒ rhs` (the transformed host plus the
    /// lhs↔rhs comap). DPO: a deleted host atom must carry no localized bond the rule does not also
    /// delete (else `ApplyError::Dangling`). Created atoms/bonds are appended, preserved entities are
    /// mutated in place, deleted entities are removed (the host renumbers). Molecule-level constraints
    /// are added/removed with their entity refs re-anchored through the match (lhs → host, created →
    /// appended); transact's renumbering compacts them on removal.
    pub fn apply_at(
        &self,
        host: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> Result<ReactionDerivation, ApplyError> {
        let deltas = self.deltas.clone().canonicalize()?;
        self.apply_at_canonical(host, correspondence, deltas)
    }

    fn apply_at_canonical(
        &self,
        host: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
        mut deltas: Deltas,
    ) -> Result<ReactionDerivation, ApplyError> {
        // A stereo coset is stated relative to a ligand ordering; the rule writes its cosets in the
        // rule's frame, the host stores the matched center in its own. Restate the rule's absolute
        // stereo deltas into the host frame before lowering (identity when the frames agree).
        reframe_stereo(&mut deltas, &self.lhs, host, correspondence);
        // Host id of a matched pattern entity (total-on-pattern, so always `Some`).
        let host_atom = |id: AtomId| {
            AtomId::from(
                correspondence
                    .atoms()
                    .right_of(NodeId::from(id))
                    .expect("a matched pattern atom maps to a host atom"),
            )
        };
        let host_bond = |id: BondId| {
            correspondence
                .bonds()
                .right_of(id)
                .expect("matched bond maps to host")
        };
        let host_dative = |id: DativeBondId| {
            correspondence
                .dative_bonds()
                .right_of(id)
                .expect("matched dative bond maps to host")
        };
        let host_aromatic = |id: AromaticSystemId| {
            correspondence
                .aromatic_systems()
                .right_of(id)
                .expect("matched aromatic system maps to host")
        };
        let host_multicenter = |id: MulticenterBondId| {
            correspondence
                .multicenter_bonds()
                .right_of(id)
                .expect("matched multicenter bond maps to host")
        };
        let host_noncovalent = |id: NoncovalentBondId| {
            correspondence
                .noncovalent_bonds()
                .right_of(id)
                .expect("matched noncovalent bond maps to host")
        };
        let host_stereo_atom = |id: StereoAtomId| {
            correspondence
                .stereo_atoms()
                .right_of(id)
                .expect("matched stereo atom maps to host")
        };
        let host_stereo_bond = |id: StereoBondId| {
            correspondence
                .stereo_bonds()
                .right_of(id)
                .expect("matched stereo bond maps to host")
        };

        let mut created_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut created_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut sets: Vec<Edit> = Vec::new();
        let mut remove_atoms: Vec<AtomHandle> = Vec::new();
        let mut remove_bonds: Vec<BondHandle> = Vec::new();
        let mut removed_host_atoms: Vec<AtomId> = Vec::new();
        let mut removed_host_bonds: HashSet<BondId> = HashSet::new();
        let mut removed_host_dative: HashSet<DativeBondId> = HashSet::new();
        let mut removed_host_aromatic: HashSet<AromaticSystemId> = HashSet::new();
        let mut removed_host_multicenter: HashSet<MulticenterBondId> = HashSet::new();
        let mut removed_host_noncovalent: HashSet<NoncovalentBondId> = HashSet::new();
        let mut removed_host_stereo_atom: HashSet<StereoAtomId> = HashSet::new();
        let mut removed_host_stereo_bond: HashSet<StereoBondId> = HashSet::new();
        let mut constraint_deltas: Vec<ConstraintDelta> = Vec::new();

        for delta in deltas.iter() {
            match delta {
                Delta::Atom(AtomDelta::Add { id, ast }) => {
                    created_atoms.insert(*id, ast.clone());
                }
                Delta::Atom(AtomDelta::Remove { id, .. }) => {
                    let removed = host_atom(*id);
                    removed_host_atoms.push(removed);
                    remove_atoms.push(AtomHandle::Id(removed));
                }
                Delta::Atom(AtomDelta::ModifyField { id, change }) => {
                    let update = match change {
                        AtomFieldChange::Element { new, .. } => AtomUpdate {
                            element: Some(new.clone()),
                            ..Default::default()
                        },
                        AtomFieldChange::IsotopeMass { new, .. } => AtomUpdate {
                            isotope_mass: Some(new.clone()),
                            ..Default::default()
                        },
                        AtomFieldChange::Charge { new, .. } => AtomUpdate {
                            charge: Some(new.clone()),
                            ..Default::default()
                        },
                        AtomFieldChange::ImplicitHydrogens { new, .. } => AtomUpdate {
                            implicit_hydrogens: Some(new.clone()),
                            ..Default::default()
                        },
                        AtomFieldChange::LonePairs { new, .. } => AtomUpdate {
                            lone_pairs: Some(new.clone()),
                            ..Default::default()
                        },
                        AtomFieldChange::Spin { old, new } => AtomUpdate {
                            spin: old.difference_to(new),
                            ..Default::default()
                        },
                    };
                    let host_id = host_atom(*id);
                    sets.extend(Edit::for_atom_update(
                        AtomHandle::Id(host_id),
                        host.atom(host_id).ast,
                        &update,
                    ));
                }
                Delta::Atom(AtomDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_atom(*id);
                        sets.extend(Edit::for_atom_update(
                            AtomHandle::Id(host_id),
                            host.atom(host_id).ast,
                            &AtomUpdate {
                                constraints: constraint.into(),
                                ..Default::default()
                            },
                        ));
                    }
                }
                Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                    created_bonds.insert(*id, (*atoms, ast.clone()));
                }
                Delta::Bond(BondDelta::Remove { id, .. }) => {
                    let removed = host_bond(*id);
                    removed_host_bonds.insert(removed);
                    remove_bonds.push(BondHandle::Id(removed));
                }
                Delta::Bond(BondDelta::ModifyField { id, change }) => {
                    let update = match change {
                        BondFieldChange::Order { new, .. } => BondUpdate {
                            order: Some(new.clone()),
                            ..Default::default()
                        },
                        BondFieldChange::Charge { new, .. } => BondUpdate {
                            charge: Some(new.clone()),
                            ..Default::default()
                        },
                        BondFieldChange::Spin { old, new } => BondUpdate {
                            spin: old.difference_to(new),
                            ..Default::default()
                        },
                    };
                    let host_id = host_bond(*id);
                    sets.extend(Edit::for_bond_update(
                        BondHandle::Id(host_id),
                        host.bond(host_id).ast,
                        &update,
                    ));
                }
                Delta::Bond(BondDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_bond(*id);
                        sets.extend(Edit::for_bond_update(
                            BondHandle::Id(host_id),
                            host.bond(host_id).ast,
                            &BondUpdate {
                                constraints: constraint.into(),
                                ..Default::default()
                            },
                        ));
                    }
                }
                Delta::DativeBond(d) => match d {
                    DativeBondDelta::ModifyField { id, change } => {
                        let update = match change {
                            DativeBondFieldChange::Order { new, .. } => DativeBondUpdate {
                                order: Some(new.clone()),
                                ..Default::default()
                            },
                        };
                        let host_id = host_dative(*id);
                        sets.extend(Edit::for_dative_bond_update(
                            DativeBondHandle::Id(host_id),
                            host.dative_bond(host_id).ast,
                            &update,
                        ));
                    }
                    DativeBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_dative(*id);
                            sets.extend(Edit::for_dative_bond_update(
                                DativeBondHandle::Id(host_id),
                                host.dative_bond(host_id).ast,
                                &DativeBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            ));
                        }
                    }
                    DativeBondDelta::Add { .. } => {}
                    DativeBondDelta::Remove { id, .. } => {
                        removed_host_dative.insert(host_dative(*id));
                    }
                },
                Delta::AromaticSystem(a) => match a {
                    AromaticSystemDelta::ModifyField { id, change } => {
                        let update = match change {
                            AromaticSystemFieldChange::Electrons { new, .. } => {
                                AromaticSystemUpdate {
                                    electrons: Some(new.clone()),
                                    ..Default::default()
                                }
                            }
                            AromaticSystemFieldChange::Charge { new, .. } => AromaticSystemUpdate {
                                charge: Some(new.clone()),
                                ..Default::default()
                            },
                            AromaticSystemFieldChange::Spin { old, new } => AromaticSystemUpdate {
                                spin: old.difference_to(new),
                                ..Default::default()
                            },
                        };
                        let host_id = host_aromatic(*id);
                        sets.extend(Edit::for_aromatic_system_update(
                            AromaticSystemHandle::Id(host_id),
                            host.aromatic_system(host_id).ast,
                            &update,
                        ));
                    }
                    AromaticSystemDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_aromatic(*id);
                            sets.extend(Edit::for_aromatic_system_update(
                                AromaticSystemHandle::Id(host_id),
                                host.aromatic_system(host_id).ast,
                                &AromaticSystemUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            ));
                        }
                    }
                    AromaticSystemDelta::Add { .. } => {}
                    AromaticSystemDelta::Remove { id, .. } => {
                        removed_host_aromatic.insert(host_aromatic(*id));
                    }
                },
                Delta::MulticenterBond(mc) => match mc {
                    MulticenterBondDelta::ModifyField { id, change } => {
                        let update = match change {
                            MulticenterBondFieldChange::Electrons { new, .. } => {
                                MulticenterBondUpdate {
                                    electrons: Some(new.clone()),
                                    ..Default::default()
                                }
                            }
                            MulticenterBondFieldChange::Charge { new, .. } => {
                                MulticenterBondUpdate {
                                    charge: Some(new.clone()),
                                    ..Default::default()
                                }
                            }
                            MulticenterBondFieldChange::Spin { old, new } => {
                                MulticenterBondUpdate {
                                    spin: old.difference_to(new),
                                    ..Default::default()
                                }
                            }
                        };
                        let host_id = host_multicenter(*id);
                        sets.extend(Edit::for_multicenter_bond_update(
                            MulticenterBondHandle::Id(host_id),
                            host.multicenter_bond(host_id).ast,
                            &update,
                        ));
                    }
                    MulticenterBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_multicenter(*id);
                            sets.extend(Edit::for_multicenter_bond_update(
                                MulticenterBondHandle::Id(host_id),
                                host.multicenter_bond(host_id).ast,
                                &MulticenterBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            ));
                        }
                    }
                    MulticenterBondDelta::Add { .. } => {}
                    MulticenterBondDelta::Remove { id, .. } => {
                        removed_host_multicenter.insert(host_multicenter(*id));
                    }
                },
                Delta::NoncovalentBond(nc) => match nc {
                    NoncovalentBondDelta::ModifyField { id, change } => {
                        let update = match change {
                            NoncovalentBondFieldChange::Kind { new, .. } => NoncovalentBondUpdate {
                                kind: Some(new.clone()),
                                ..Default::default()
                            },
                        };
                        let host_id = host_noncovalent(*id);
                        sets.extend(Edit::for_noncovalent_bond_update(
                            NoncovalentBondHandle::Id(host_id),
                            host.noncovalent_bond(host_id).ast,
                            &update,
                        ));
                    }
                    NoncovalentBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_noncovalent(*id);
                            sets.extend(Edit::for_noncovalent_bond_update(
                                NoncovalentBondHandle::Id(host_id),
                                host.noncovalent_bond(host_id).ast,
                                &NoncovalentBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            ));
                        }
                    }
                    NoncovalentBondDelta::Add { .. } => {}
                    NoncovalentBondDelta::Remove { id, .. } => {
                        removed_host_noncovalent.insert(host_noncovalent(*id));
                    }
                },
                // Stereo: the four set-ops lower directly; the relative ops resolve against the
                // matched host config (same frame — no reindex, like the other overlays) and emit an
                // absolute `Configuration`. `Add` is lowered in the second pass; `Remove` tracks the
                // host id for the DPO dangling check.
                Delta::StereoAtom(s) => match s {
                    StereoAtomDelta::ModifyField { id, change } => {
                        let host_id = host_stereo_atom(*id);
                        let StereoAtomFieldChange::Configuration { new, .. } = change;
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration {
                                old: host.stereo_atom(host_id).ast.configuration.clone(),
                                new: new.clone(),
                            },
                        })
                    }
                    StereoAtomDelta::ModifyConstraint { id, old, new, .. } => {
                        if let Some(constraint) = new.as_ref().or(old.as_ref()) {
                            let host_id = host_stereo_atom(*id);
                            sets.push(Edit::ModifyStereoAtomConstraint {
                                id: StereoAtomHandle::Id(host_id),
                                old: host
                                    .stereo_atom(host_id)
                                    .ast
                                    .constraints
                                    .get(constraint.key())
                                    .cloned(),
                                new: new.clone(),
                            })
                        }
                    }
                    StereoAtomDelta::Apply {
                        id,
                        kind,
                        permutation,
                    } => {
                        let host_id = host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.apply(*permutation);
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Swap { id, kind } => {
                        let host_id = host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.swap();
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Mirror { id, kind } => {
                        let host_id = host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.mirror();
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Add { .. } => {}
                    StereoAtomDelta::Remove { id, .. } => {
                        removed_host_stereo_atom.insert(host_stereo_atom(*id));
                    }
                },
                Delta::StereoBond(s) => match s {
                    StereoBondDelta::ModifyField { id, change } => {
                        let host_id = host_stereo_bond(*id);
                        let StereoBondFieldChange::Configuration { new, .. } = change;
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondHandle::Id(host_id),
                            change: StereoBondFieldChange::Configuration {
                                old: host.stereo_bond(host_id).ast.configuration.clone(),
                                new: new.clone(),
                            },
                        })
                    }
                    StereoBondDelta::ModifyConstraint { id, old, new, .. } => {
                        if let Some(constraint) = new.as_ref().or(old.as_ref()) {
                            let host_id = host_stereo_bond(*id);
                            sets.push(Edit::ModifyStereoBondConstraint {
                                id: StereoBondHandle::Id(host_id),
                                old: host
                                    .stereo_bond(host_id)
                                    .ast
                                    .constraints
                                    .get(constraint.key())
                                    .cloned(),
                                new: new.clone(),
                            })
                        }
                    }
                    StereoBondDelta::Apply {
                        id,
                        kind,
                        permutation,
                    } => {
                        let host_id = host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.apply(*permutation);
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondHandle::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Swap { id, kind } => {
                        let host_id = host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.swap();
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondHandle::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Mirror { id, kind } => {
                        let host_id = host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.mirror();
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondHandle::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Add { .. } => {}
                    StereoBondDelta::Remove { id, .. } => {
                        removed_host_stereo_bond.insert(host_stereo_bond(*id));
                    }
                },
                Delta::Constraint(c) => constraint_deltas.push(c.clone()),
            }
        }

        // DPO gluing condition: a deleted host atom keeps no bond or overlay the rule does not
        // also delete.
        for &host_atom in &removed_host_atoms {
            let atom = host.atom(host_atom);
            for bond in atom.bond_ids() {
                if !removed_host_bonds.contains(&bond) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            for dative in atom.dative_bond_ids() {
                if !removed_host_dative.contains(&dative) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            if let Some(aromatic) = atom.aromatic_system_id() {
                if !removed_host_aromatic.contains(&aromatic) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            for multicenter in atom.multicenter_bond_ids() {
                if !removed_host_multicenter.contains(&multicenter) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            for noncovalent in atom.noncovalent_bond_ids() {
                if !removed_host_noncovalent.contains(&noncovalent) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            // Stereo incidence (site or ligand) via the stereo views; a stereo bond's site is a bond,
            // so a deleted atom touches a stereo bond only as a ligand — `incident_ids` covers both.
            for stereo_atom in host.stereo_atoms().incident_ids(host_atom) {
                if !removed_host_stereo_atom.contains(&stereo_atom) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
            for stereo_bond in host.stereo_bonds().incident_ids(host_atom) {
                if !removed_host_stereo_bond.contains(&stereo_bond) {
                    return Err(ApplyError::Dangling { host_atom });
                }
            }
        }

        // `AddAtoms` is the first edit, so created atoms take `New(0..k)` in ascending id order.
        let new_atom_index: HashMap<AtomId, usize> = created_atoms
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();
        let atom_handle = |id: AtomId| match new_atom_index.get(&id) {
            Some(&index) => AtomHandle::New(index),
            None => AtomHandle::Id(host_atom(id)),
        };
        // `AddBonds` follows `AddAtoms`, so a created bond continues the shared `New(..)` created-entity
        // numbering after the created atoms (a stereo bond's `New` site indexes this joint list).
        let new_bond_index: HashMap<BondId, usize> = created_bonds
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, created_atoms.len() + index))
            .collect();
        let bond_handle = |id: BondId| match new_bond_index.get(&id) {
            Some(&index) => BondHandle::New(index),
            None => BondHandle::Id(host_bond(id)),
        };

        // Overlay create/remove need `atom_handle` (created participants resolve to `New`), so they
        // are lowered in a second pass: adds after the topology adds, removes before
        // `RemoveTopology`. Removes are collected per kind and emitted as one batched edit each,
        // so each overlay id space is compacted once against the pre-removal state (a sequence of
        // single-id removes would stale the not-yet-processed ids). Dative `atoms` is
        // `[donors…, acceptor]` (acceptor last, per transact).
        let mut overlay_adds: Vec<Edit> = Vec::new();
        let mut remove_dative: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondAst)> = Vec::new();
        let mut remove_aromatic: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemAst)> =
            Vec::new();
        let mut remove_multicenter: Vec<(
            MulticenterBondHandle,
            Vec<AtomHandle>,
            MulticenterBondAst,
        )> = Vec::new();
        let mut remove_noncovalent: Vec<(
            NoncovalentBondHandle,
            [AtomHandle; 2],
            NoncovalentBondAst,
        )> = Vec::new();
        let mut remove_stereo_atom: Vec<StereoAtomRemoval> = Vec::new();
        let mut remove_stereo_bond: Vec<StereoBondRemoval> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::DativeBond(DativeBondDelta::Add {
                    donors,
                    acceptor,
                    ast,
                    ..
                }) => {
                    let mut atoms: Vec<AtomHandle> =
                        donors.iter().map(|a| atom_handle(*a)).collect();
                    atoms.push(atom_handle(*acceptor));
                    overlay_adds.push(Edit::AddDativeBond {
                        atoms,
                        ast: ast.clone(),
                    });
                }
                Delta::DativeBond(DativeBondDelta::Remove {
                    id,
                    donors,
                    acceptor,
                    ast,
                }) => {
                    let mut atoms: Vec<AtomHandle> =
                        donors.iter().map(|a| atom_handle(*a)).collect();
                    atoms.push(atom_handle(*acceptor));
                    remove_dative.push((
                        DativeBondHandle::Id(host_dative(*id)),
                        atoms,
                        ast.clone(),
                    ));
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddAromaticSystem {
                        atoms: atoms.iter().map(|a| atom_handle(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, ast }) => {
                    remove_aromatic.push((
                        AromaticSystemHandle::Id(host_aromatic(*id)),
                        atoms.iter().map(|a| atom_handle(*a)).collect(),
                        ast.clone(),
                    ));
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddMulticenterBond {
                        atoms: atoms.iter().map(|a| atom_handle(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, ast }) => {
                    remove_multicenter.push((
                        MulticenterBondHandle::Id(host_multicenter(*id)),
                        atoms.iter().map(|a| atom_handle(*a)).collect(),
                        ast.clone(),
                    ));
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddNoncovalentBond {
                        atoms: [atom_handle(atoms[0]), atom_handle(atoms[1])],
                        ast: ast.clone(),
                    });
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, ast }) => {
                    remove_noncovalent.push((
                        NoncovalentBondHandle::Id(host_noncovalent(*id)),
                        [atom_handle(atoms[0]), atom_handle(atoms[1])],
                        ast.clone(),
                    ));
                }
                Delta::StereoAtom(StereoAtomDelta::Add {
                    site, ligands, ast, ..
                }) => {
                    overlay_adds.push(Edit::AddStereoAtom {
                        site: atom_handle(*site),
                        ligands: ligands
                            .iter()
                            .map(|l| (atom_handle(l.atom_id), l.kind))
                            .collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::StereoAtom(StereoAtomDelta::Remove {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    remove_stereo_atom.push((
                        StereoAtomHandle::Id(host_stereo_atom(*id)),
                        atom_handle(*site),
                        ligands
                            .iter()
                            .map(|l| (atom_handle(l.atom_id), l.kind))
                            .collect(),
                        ast.clone(),
                    ));
                }
                Delta::StereoBond(StereoBondDelta::Add {
                    site, ligands, ast, ..
                }) => {
                    overlay_adds.push(Edit::AddStereoBond {
                        site: bond_handle(*site),
                        ligands: ligands
                            .iter()
                            .map(|l| (atom_handle(l.atom_id), l.kind))
                            .collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::StereoBond(StereoBondDelta::Remove {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    remove_stereo_bond.push((
                        StereoBondHandle::Id(host_stereo_bond(*id)),
                        bond_handle(*site),
                        ligands
                            .iter()
                            .map(|l| (atom_handle(l.atom_id), l.kind))
                            .collect(),
                        ast.clone(),
                    ));
                }
                _ => {}
            }
        }
        let mut overlay_removes: Vec<Edit> = Vec::new();
        if !remove_dative.is_empty() {
            overlay_removes.push(Edit::RemoveDativeBonds {
                removes: remove_dative,
            });
        }
        if !remove_aromatic.is_empty() {
            overlay_removes.push(Edit::RemoveAromaticSystems {
                removes: remove_aromatic,
            });
        }
        if !remove_multicenter.is_empty() {
            overlay_removes.push(Edit::RemoveMulticenterBonds {
                removes: remove_multicenter,
            });
        }
        if !remove_noncovalent.is_empty() {
            overlay_removes.push(Edit::RemoveNoncovalentBonds {
                removes: remove_noncovalent,
            });
        }
        if !remove_stereo_atom.is_empty() {
            overlay_removes.push(Edit::RemoveStereoAtoms {
                removes: remove_stereo_atom,
            });
        }
        if !remove_stereo_bond.is_empty() {
            overlay_removes.push(Edit::RemoveStereoBonds {
                removes: remove_stereo_bond,
            });
        }

        // Molecule-level constraints lower to `Edit::{Add,Remove}MoleculeConstraint`, refs
        // re-anchored through the match: lhs entities → host (via `m`), created atoms/bonds → their
        // appended host id. Emitted before all removals (overlay + topology) so each removal's
        // constraint compaction updates them (remapping surviving refs, dropping refs to a deleted
        // entity). The overlay/stereo maps cover only `lhs` overlays — a constraint referencing a
        // rule-created overlay is unsupported.
        let mut constraint_edits: Vec<Edit> = Vec::new();
        if !constraint_deltas.is_empty() {
            let host_atom_count = host.atoms().count();
            let host_bond_count = host.bonds().count();
            let mut atom: HashMap<AtomId, AtomId> = (0..self.lhs.atoms().count() as u32)
                .map(|i| (AtomId(i), host_atom(AtomId(i))))
                .collect();
            for (&created, &index) in &new_atom_index {
                atom.insert(created, AtomId((host_atom_count + index) as u32));
            }
            let mut bond: HashMap<BondId, BondId> = (0..self.lhs.bonds().count() as u32)
                .map(|i| (BondId(i), host_bond(BondId(i))))
                .collect();
            for (index, &created) in created_bonds.keys().enumerate() {
                bond.insert(created, BondId((host_bond_count + index) as u32));
            }
            let match_map = IdRemapping::new(
                atom,
                bond,
                (0..self.lhs.dative_bonds().count() as u32)
                    .map(|i| (DativeBondId(i), host_dative(DativeBondId(i))))
                    .collect(),
                (0..self.lhs.aromatic_systems().count() as u32)
                    .map(|i| (AromaticSystemId(i), host_aromatic(AromaticSystemId(i))))
                    .collect(),
                (0..self.lhs.multicenter_bonds().count() as u32)
                    .map(|i| (MulticenterBondId(i), host_multicenter(MulticenterBondId(i))))
                    .collect(),
                (0..self.lhs.noncovalent_bonds().count() as u32)
                    .map(|i| (NoncovalentBondId(i), host_noncovalent(NoncovalentBondId(i))))
                    .collect(),
                (0..self.lhs.stereo_atoms().count() as u32)
                    .map(|i| (StereoAtomId(i), host_stereo_atom(StereoAtomId(i))))
                    .collect(),
                (0..self.lhs.stereo_bonds().count() as u32)
                    .map(|i| (StereoBondId(i), host_stereo_bond(StereoBondId(i))))
                    .collect(),
            );
            for delta in constraint_deltas {
                match delta {
                    ConstraintDelta::Add(c) => constraint_edits.push(Edit::AddMoleculeConstraint {
                        constraint: c.remap(&match_map),
                    }),
                    ConstraintDelta::Remove(c) => {
                        constraint_edits.push(Edit::RemoveMoleculeConstraint {
                            constraint: c.remap(&match_map),
                        })
                    }
                }
            }
        }

        let mut edits: Vec<Edit> = Vec::new();
        if !created_atoms.is_empty() {
            edits.push(Edit::AddAtoms {
                atoms: created_atoms.values().cloned().collect(),
            });
        }
        if !created_bonds.is_empty() {
            edits.push(Edit::AddBonds {
                bonds: created_bonds
                    .values()
                    .map(|(atoms, ast)| AddBond {
                        endpoints: [atom_handle(atoms[0]), atom_handle(atoms[1])],
                        ast: ast.clone(),
                    })
                    .collect(),
            });
        }
        edits.extend(overlay_adds);
        edits.extend(sets);
        // Constraints precede all removals (overlay and topology) so each removal's constraint
        // compaction updates them — a constraint referencing a surviving overlay whose lower-id
        // sibling is removed would otherwise carry a stale id.
        edits.extend(constraint_edits);
        edits.extend(overlay_removes);
        if !remove_atoms.is_empty() || !remove_bonds.is_empty() {
            edits.push(Edit::RemoveTopology {
                atoms: remove_atoms,
                bonds: remove_bonds,
            });
        }

        let mut builder = host.edit();
        builder.transact(edits)?;
        let product = builder.build();

        // Emit-compliance: the product is a generated molecule, so it must satisfy every per-entity
        // structural invariant (a rule's adds can land a parallel bond, an overlapping system, or a
        // second stereo center on an occupied site). The per-entity `has_conflict` primitives are the
        // shared gates (also consulted by the validator and `meet_pushout`); enforced per generating op
        // pending a single central emit gate.
        if product.bonds().has_conflict()
            || product.dative_bonds().has_conflict()
            || product.aromatic_systems().has_conflict()
            || product.multicenter_bonds().has_conflict()
            || product.noncovalent_bonds().has_conflict()
            || product.stereo_atoms().has_conflict()
            || product.stereo_bonds().has_conflict()
        {
            return Err(ApplyError::StructuralConflict);
        }

        // The host↔product comap: preserved host atoms mate to their compacted product id (survivors
        // keep ascending order), removed atoms are left-exposed, created atoms right-exposed. `induce`
        // derives the bond and overlay correspondences from this atom map.
        let removed: HashSet<AtomId> = removed_host_atoms.iter().copied().collect();
        let mut atom_mates: Vec<(NodeId, NodeId)> = Vec::new();
        let mut product_atom = 0u32;
        for host_atom in 0..host.atoms().count() as u32 {
            if removed.contains(&AtomId(host_atom)) {
                continue;
            }
            atom_mates.push((NodeId(host_atom), NodeId(product_atom)));
            product_atom += 1;
        }
        let atom_map =
            Correspondence::new(atom_mates, host.atoms().count(), product.atoms().count());
        let comap = MoleculeCorrespondence::induce(host, &product, atom_map);
        Ok(ReactionDerivation::new(host.clone(), product, comap))
    }

    /// Validate the structural preconditions shared by every match against `host`.
    pub fn validate_application(&self, host: &MoleculeAst) -> Result<(), ApplyPreconditionError> {
        self.application_deltas(host).map(drop)
    }

    fn application_deltas(&self, host: &MoleculeAst) -> Result<Deltas, ApplyPreconditionError> {
        let deltas = self
            .deltas
            .clone()
            .canonicalize()
            .map_err(|_| ApplyPreconditionError::InconsistentReaction)?;

        let reaction_integrity = match ReactionIntegrityValidator.validate(&self.lhs, &deltas) {
            Ok(outcome) => outcome,
            Err(error) => match error {},
        };
        reaction_integrity
            .into_observation()
            .map_err(|contradiction| match contradiction {
                ReactionIntegrityContradiction::InvalidReference { entity } => {
                    ApplyPreconditionError::InvalidReactionReference { entity }
                }
                ReactionIntegrityContradiction::IncidenceMismatch { entity } => {
                    ApplyPreconditionError::ReactionIncidenceMismatch { entity }
                }
            })?;

        let lhs_structure = match EntityStructureValidator.validate(&self.lhs) {
            Ok(outcome) => outcome,
            Err(error) => match error {},
        };
        lhs_structure
            .into_observation()
            .map_err(ApplyPreconditionError::ReactionStructure)?;

        let dpo = match DpoValidator.validate_reaction(&self.lhs, &deltas) {
            Ok(outcome) => outcome,
            Err(error) => match error {},
        };
        dpo.into_observation()
            .map_err(ApplyPreconditionError::ReactionDpo)?;

        let host_structure = match EntityStructureValidator.validate(host) {
            Ok(outcome) => outcome,
            Err(error) => match error {},
        };
        host_structure
            .into_observation()
            .map_err(ApplyPreconditionError::HostStructure)?;

        Ok(deltas)
    }

    /// Every product of applying the reaction to `host`: one per injective match of `lhs` into
    /// `host` (via `subiso`) that satisfies the match-local DPO and structural conditions.
    /// Structural preconditions are checked before match enumeration. Match-local rejection is
    /// skipped; an internal application failure is yielded once and terminates the iterator.
    pub fn apply<'h>(
        &'h self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Result<
        impl Iterator<Item = Result<ReactionDerivation, ApplyError>> + 'h,
        ApplyPreconditionError,
    > {
        let deltas = self.application_deltas(host)?;
        let mut correspondences = self
            .lhs
            .substructure_matches(host, SubstructureMatchAlgorithm::GraphAndOverlays, subiso)
            .into_iter();
        let mut failed = false;

        Ok(from_fn(move || {
            while !failed {
                let correspondence = correspondences.next()?;
                match self.apply_at_canonical(host, &correspondence, deltas.clone()) {
                    Ok(derivation) => return Some(Ok(derivation)),
                    Err(error) if error.is_match_rejection() => {}
                    Err(error) => {
                        failed = true;
                        return Some(Err(error));
                    }
                }
            }
            None
        }))
    }
}

/// Restate `deltas`' absolute stereo cosets from the rule (`lhs`) frame into the matched `host` frame.
/// The coset is meaningful only per ligand ordering, so a `ModifyField`/`Remove` delta lowered onto a
/// host whose matching center is numbered differently must carry its cosets across — the delta-side
/// mirror of the matcher's `coset_for`. `before` is the rule's ligand order mapped into the host id
/// space, `after` the host's stored order; identity when they agree. The relative ops
/// (`Apply`/`Swap`/`Mirror`) resolve against the host coset, `Add` creates a fresh overlay, and stereo
/// constraints are positionless — none are reframed; a delta with no host correspondent is skipped.
fn reframe_stereo(
    deltas: &mut Deltas,
    lhs: &MoleculeAst,
    host: &MoleculeAst,
    correspondence: &MoleculeCorrespondence,
) {
    let into_host = |l: &StereoLigand| {
        StereoLigand::new(
            AtomId::from(
                correspondence
                    .atoms()
                    .right_of(NodeId::from(l.atom_id))
                    .expect("a matched rule ligand maps into the host"),
            ),
            l.kind,
        )
    };
    let from_host = |l: &StereoLigand| {
        StereoLigand::new(
            AtomId::from(
                correspondence
                    .atoms()
                    .left_of(NodeId::from(l.atom_id))
                    .expect("a matched host ligand maps back to the rule"),
            ),
            l.kind,
        )
    };
    for delta in deltas.iter_mut() {
        match delta {
            Delta::StereoAtom(s) => {
                let Some(host_id) = correspondence.stereo_atoms().right_of(s.id()) else {
                    continue;
                };
                let before: Vec<StereoLigand> = lhs
                    .stereo_atom(s.id())
                    .ligand_frame()
                    .iter()
                    .map(into_host)
                    .collect();
                let after = host.stereo_atom(host_id).ligand_frame();
                match s {
                    StereoAtomDelta::ModifyField {
                        change: StereoAtomFieldChange::Configuration { old, new },
                        ..
                    } => {
                        let sigma = Permutation::between(&before, &after);
                        *old = old.apply(sigma);
                        *new = new.apply(sigma);
                    }
                    StereoAtomDelta::Remove { ligands, ast, .. } => {
                        *ast = ast.transform_frame(&before, &after);
                        *ligands = after.iter().map(from_host).collect();
                    }
                    _ => {}
                }
            }
            Delta::StereoBond(s) => {
                let Some(host_id) = correspondence.stereo_bonds().right_of(s.id()) else {
                    continue;
                };
                let before: Vec<StereoLigand> = lhs
                    .stereo_bond(s.id())
                    .ligand_frame()
                    .iter()
                    .map(into_host)
                    .collect();
                let after = host.stereo_bond(host_id).ligand_frame();
                match s {
                    StereoBondDelta::ModifyField {
                        change: StereoBondFieldChange::Configuration { old, new },
                        ..
                    } => {
                        let sigma = Permutation::between(&before, &after);
                        *old = old.apply(sigma);
                        *new = new.apply(sigma);
                    }
                    StereoBondDelta::Remove { ligands, ast, .. } => {
                        *ast = ast.transform_frame(&before, &after);
                        *ligands = after.iter().map(from_host).collect();
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }
}

impl Canonicalize for ReactionAst {
    /// Value-level in a fixed atom id space: `deltas` are canonicalized;
    /// `lhs` is passed through (`MoleculeAst` has no whole-molecule canonical form — its
    /// equality is structural). Equality up to atom renumbering is a separate `umol-graph`
    /// operation.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            lhs: self.lhs,
            deltas: self.deltas.canonicalize()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::constraint::{Constraint, Constraints, MoleculeConstraint};
    use super::super::edit::{AtomFieldChange, BondFieldChange};
    use super::super::entity::Entity;
    use super::super::ligand::StereoLigandKind;
    use super::super::molecule::transact::TransactionError;
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use super::super::stereo::{StereoAtomAst, StereoBondAst, StereoCosetAst, StereoKind};
    use super::super::validate::{DpoContradiction, EntityStructureContradiction};
    use super::super::value::ValueAst;
    use super::*;

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(id),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(old),
                new: ValueAst::Lit(new),
            },
        })
    }

    #[rstest]
    fn test_reaction_ast_from_sides() {
        // C-C (order 1) → C-C (order 2) under the total atom correspondence: one bond-order modify.
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2);
        assert_eq!(
            ReactionAst::from_sides(left.clone(), right, atoms),
            ReactionAst::new(
                left,
                Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: ValueAst::Lit(1),
                        new: ValueAst::Lit(2),
                    },
                })]),
            ),
        );
    }

    #[rstest]
    fn test_reaction_ast_canonicalize() {
        // The delta chain fuses; the lhs is passed through unchanged.
        let reaction = ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)]),
        );
        assert_eq!(
            reaction.canonicalize().unwrap(),
            ReactionAst::new(
                MoleculeAst::default(),
                Deltas::from_iter([charge_set(0, 0, 2)])
            ),
        );
    }

    #[rstest]
    #[case::bond_order(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() }),
    )]
    #[case::overlay_removed(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::O) }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O)], bonds: vec![], ..Default::default() }),
    )]
    fn test_reaction_ast_apply_at(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] atom_map: Vec<AtomId>,
        #[case] expected: MoleculeAst,
    ) {
        let atom_images: Vec<NodeId> = atom_map.iter().map(|&a| NodeId::from(a)).collect();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&atom_images, host.atoms().count()),
        );
        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected
        );
    }

    // `dangling_*`: the rule deletes a host atom still carrying an undeleted bond/overlay (DPO gluing
    // condition). `structural_conflict`: the rule's add lands a second bond on an already-bonded atom
    // pair, so the product would carry parallel bonds — an emit-compliance invariant (`has_conflict`).
    #[rstest]
    #[case::dangling_bond(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
        vec![AtomId(0)],
        ApplyError::Dangling { host_atom: AtomId(0) },
    )]
    #[case::dangling_noncovalent(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::O),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        vec![AtomId(0)],
        ApplyError::Dangling { host_atom: AtomId(0) },
    )]
    #[case::structural_conflict(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::from_order(1),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        ApplyError::StructuralConflict,
    )]
    fn test_reaction_ast_apply_at_error(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] atom_map: Vec<AtomId>,
        #[case] expected: ApplyError,
    ) {
        let images: Vec<NodeId> = atom_map.iter().map(|&a| NodeId::from(a)).collect();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&images, host.atoms().count()),
        );
        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap_err(),
            expected
        );
    }

    #[rstest]
    fn test_reaction_ast_apply_at_molecule_constraint() {
        // A reaction adding a molecule-level `ChargeSum` over its lhs atoms; applied at a match
        // that maps lhs atoms 0,1 → host atoms 1,2, the constraint's refs re-anchor to the host.
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                    sum: ValueAst::Lit(0),
                }),
            ))]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&[NodeId(1), NodeId(2)], host.atoms().count()),
        );
        let result = reaction.apply_at(&host, &correspondence).unwrap();
        assert_eq!(
            result.rhs().constraints(),
            &Constraints::from(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
                sum: ValueAst::Lit(0),
            })),
        );
    }

    #[rstest]
    #[case::valid(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
            Deltas::new(),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
    )]
    #[case::canonical_add_remove_cancellation(
        ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
                Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
            ]),
        ),
        MoleculeAst::default(),
    )]
    #[case::unordered_bond_incidence(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(1), AtomId(0)],
                ast: BondAst::from_order(1),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
    )]
    fn test_reaction_ast_validate_application(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
    ) {
        assert_eq!(reaction.validate_application(&host), Ok(()));
    }

    #[rstest]
    #[case::inconsistent_reaction(
        ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
                Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::O) }),
            ]),
        ),
        MoleculeAst::default(),
        ApplyPreconditionError::InconsistentReaction,
    )]
    #[case::reaction_structure(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(0), AtomId(1), BondAst::from_order(2)),
                ],
                ..Default::default()
            }),
            Deltas::new(),
        ),
        MoleculeAst::default(),
        ApplyPreconditionError::ReactionStructure(EntityStructureContradiction::BondsParallel { atoms: [AtomId(0), AtomId(1)] }),
    )]
    #[case::reaction_dpo(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) })]),
        ),
        MoleculeAst::default(),
        ApplyPreconditionError::ReactionDpo(DpoContradiction::DanglingBond { atom: AtomId(0), bond: BondId(0) }),
    )]
    #[case::host_structure(
        ReactionAst::new(MoleculeAst::default(), Deltas::new()),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(1), BondAst::from_order(2)),
            ],
            ..Default::default()
        }),
        ApplyPreconditionError::HostStructure(EntityStructureContradiction::BondsParallel { atoms: [AtomId(0), AtomId(1)] }),
    )]
    fn test_reaction_ast_validate_application_error(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] expected: ApplyPreconditionError,
    ) {
        assert_eq!(reaction.validate_application(&host).unwrap_err(), expected);
    }

    #[rstest]
    #[case::atom(
        Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::default() }),
        Entity::Atom(AtomId(0)),
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Remove { id: BondId(0), atoms: [AtomId(0), AtomId(1)], ast: BondAst::default() }),
        Entity::Bond(BondId(0)),
    )]
    #[case::dative_bond(
        Delta::DativeBond(DativeBondDelta::Remove { id: DativeBondId(0), donors: vec![AtomId(0)], acceptor: AtomId(1), ast: DativeBondAst::default() }),
        Entity::DativeBond(DativeBondId(0)),
    )]
    #[case::aromatic_system(
        Delta::AromaticSystem(AromaticSystemDelta::Remove { id: AromaticSystemId(0), atoms: vec![AtomId(0)], ast: AromaticSystemAst::default() }),
        Entity::AromaticSystem(AromaticSystemId(0)),
    )]
    #[case::multicenter_bond(
        Delta::MulticenterBond(MulticenterBondDelta::Remove { id: MulticenterBondId(0), atoms: vec![AtomId(0)], ast: MulticenterBondAst::default() }),
        Entity::MulticenterBond(MulticenterBondId(0)),
    )]
    #[case::noncovalent_bond(
        Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)], ast: NoncovalentBondAst::default() }),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    )]
    #[case::stereo_atom(
        Delta::StereoAtom(StereoAtomDelta::Remove { id: StereoAtomId(0), site: AtomId(0), ligands: vec![], ast: StereoAtomAst::default() }),
        Entity::StereoAtom(StereoAtomId(0)),
    )]
    #[case::stereo_bond(
        Delta::StereoBond(StereoBondDelta::Remove { id: StereoBondId(0), site: BondId(0), ligands: vec![], ast: StereoBondAst::default() }),
        Entity::StereoBond(StereoBondId(0)),
    )]
    fn test_reaction_ast_validate_application_rejects_missing_delta_target(
        #[case] delta: Delta,
        #[case] entity: Entity,
    ) {
        let reaction = ReactionAst::new(MoleculeAst::default(), Deltas::from_iter([delta]));
        assert_eq!(
            reaction.validate_application(&MoleculeAst::default()),
            Err(ApplyPreconditionError::InvalidReactionReference { entity }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_created_id_collision() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::O),
            })]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::InvalidReactionReference {
                entity: Entity::Atom(AtomId(0)),
            }),
        );
    }

    #[rstest]
    #[case::bond_endpoint(Delta::Bond(BondDelta::Add {
        id: BondId(0),
        atoms: [AtomId(0), AtomId(1)],
        ast: BondAst::default(),
    }))]
    #[case::dative_participant(Delta::DativeBond(DativeBondDelta::Add {
        id: DativeBondId(0),
        donors: vec![AtomId(1)],
        acceptor: AtomId(0),
        ast: DativeBondAst::default(),
    }))]
    #[case::aromatic_participant(Delta::AromaticSystem(AromaticSystemDelta::Add {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        ast: AromaticSystemAst::default(),
    }))]
    #[case::multicenter_participant(Delta::MulticenterBond(MulticenterBondDelta::Add {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        ast: MulticenterBondAst::default(),
    }))]
    #[case::noncovalent_endpoint(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        ast: NoncovalentBondAst::default(),
    }))]
    #[case::stereo_atom_site(Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(1),
        ligands: vec![],
        ast: StereoAtomAst::default(),
    }))]
    #[case::stereo_atom_ligand(Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
        ast: StereoAtomAst::default(),
    }))]
    fn test_reaction_ast_validate_application_rejects_missing_structural_reference(
        #[case] delta: Delta,
    ) {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = ReactionAst::new(lhs.clone(), Deltas::from_iter([delta]));
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::InvalidReactionReference {
                entity: Entity::Atom(AtomId(1)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_missing_stereo_bond_site() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Add {
                id: StereoBondId(0),
                site: BondId(0),
                ligands: vec![],
                ast: StereoBondAst::default(),
            })]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::InvalidReactionReference {
                entity: Entity::Bond(BondId(0)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_bond_incidence_mismatch() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(2)],
                ast: BondAst::from_order(1),
            })]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::Bond(BondId(0)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_dative_incidence_mismatch() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::O),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::default())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                id: DativeBondId(0),
                donors: vec![AtomId(2)],
                acceptor: AtomId(1),
                ast: DativeBondAst::default(),
            })]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::DativeBond(DativeBondId(0)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_stereo_frame_incidence_mismatch() {
        let stored_ligands = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                stored_ligands.clone(),
                StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        });
        let mut removed_ligands = stored_ligands;
        removed_ligands.swap(0, 1);
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: removed_ligands,
                ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
            })]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::StereoAtom(StereoAtomId(0)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_validate_application_rejects_recursive_constraint_reference() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        });
        let constraint = Constraint::Not(Box::new(Constraint::And(vec![Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
            },
        )])));
        let reaction = ReactionAst::new(
            lhs.clone(),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(constraint))]),
        );
        assert_eq!(
            reaction.validate_application(&lhs),
            Err(ApplyPreconditionError::InvalidReactionReference {
                entity: Entity::Atom(AtomId(1)),
            }),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::bond_order(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) } })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
        vec![MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))], ..Default::default() })],
    )]
    #[case::match_rejection(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(1), AtomId(2), BondAst::from_order(1))], ..Default::default() }),
        vec![MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() })],
    )]
    #[case::host_relative_update(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge { old: ValueAst::Undetermined, new: ValueAst::Lit(1) },
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C).with_charge(0_i64)], ..Default::default() }),
        vec![MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C).with_charge(1_i64)], ..Default::default() })],
    )]
    fn test_reaction_ast_apply(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] expected: Vec<MoleculeAst>,
    ) {
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        assert_eq!(products, expected);
    }

    #[rstest]
    #[case::transaction(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C)],
                constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0)]),
                    sum: ValueAst::Lit(0),
                })),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0)]),
                    sum: ValueAst::Lit(0),
                }),
            ))]),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
        ApplyError::Transaction(TransactionError::MissingEntry),
    )]
    fn test_reaction_ast_apply_error(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] expected: ApplyError,
    ) {
        let mut applications = reaction
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .unwrap();

        assert_eq!(applications.next().unwrap().unwrap_err(), expected);
        assert_eq!(applications.next(), None);
    }

    #[rstest]
    #[case::host_structure(
        ReactionAst::new(MoleculeAst::default(), Deltas::new()),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(1), BondAst::from_order(2)),
            ],
            ..Default::default()
        }),
        ApplyPreconditionError::HostStructure(EntityStructureContradiction::BondsParallel { atoms: [AtomId(0), AtomId(1)] }),
    )]
    fn test_reaction_ast_apply_precondition_error(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] expected: ApplyPreconditionError,
    ) {
        match reaction.apply(&host, SubgraphIsomorphismAlgorithm::Vf2) {
            Err(error) => assert_eq!(error, expected),
            Ok(_) => panic!("invalid input passed application integrity validation"),
        }
    }

    #[fixture]
    fn tetrahedral_inversion() -> ReactionAst {
        // Invert a tetrahedral C(0) whose ligands F,Cl,Br,I are stated in ascending order: coset 0 → 1.
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::F),
                    AtomAst::from_element(Element::Cl),
                    AtomAst::from_element(Element::Br),
                    AtomAst::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(0), AtomId(2), BondAst::from_order(1)),
                    (AtomId(0), AtomId(3), BondAst::from_order(1)),
                    (AtomId(0), AtomId(4), BondAst::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(0),
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(1),
                    ),
                },
            })]),
        )
    }

    // Applying the ascending-frame inversion rule to a host that states the same center in a different
    // ligand order: the match succeeds (the matcher reframes), and `apply_at` now reframes the rule's
    // `ModifyField` coset into the host frame before lowering it, so the derivation inverts the host's
    // stored coset in the host's own frame. `same_frame` is the control; `swapped_frame` (ligands 1↔2,
    // its physically-equal coset 1) forces the reframe.
    #[rstest]
    #[case::same_frame([1, 2, 3, 4], 0, 1)]
    #[case::swapped_frame([2, 1, 3, 4], 1, 0)]
    fn test_reaction_ast_apply_stereo_cross_frame(
        tetrahedral_inversion: ReactionAst,
        #[case] host_ligands: [u32; 4],
        #[case] host_coset: u32,
        #[case] product_coset: u32,
    ) {
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
                (AtomId(0), AtomId(3), BondAst::from_order(1)),
                (AtomId(0), AtomId(4), BondAst::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                host_ligands
                    .iter()
                    .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, host_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
                (AtomId(0), AtomId(3), BondAst::from_order(1)),
                (AtomId(0), AtomId(4), BondAst::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                host_ligands
                    .iter()
                    .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, product_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let rhs = tetrahedral_inversion
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .unwrap()
            .next()
            .expect("the inversion rule matches the host")
            .unwrap()
            .rhs()
            .clone();
        assert_eq!(rhs, expected);
    }

    // Adding a stereo bond whose site is a rule-created bond: the site resolves through a `New` handle
    // into the shared created-entity list, which `AddBonds` fills after `AddAtoms` — so the created
    // bond's `New` index must clear the created atoms. Regression for the stereo-bond-only compose
    // failure (a created-bond `New` site aliasing a created atom → `RefTypeMismatch`).
    #[rstest]
    #[case::coset_0(0u32)]
    fn test_reaction_ast_apply_stereo_bond_created_site(#[case] coset: u32) {
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C)],
                bonds: vec![],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(2),
                }),
                Delta::StereoBond(StereoBondDelta::Add {
                    id: StereoBondId(0),
                    site: BondId(0),
                    ligands: vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    ],
                    ast: StereoBondAst::new(StereoKind::CisTrans, 0u32),
                }),
            ]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondAst::new(StereoKind::CisTrans, coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let rhs = reaction
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .unwrap()
            .next()
            .expect("the reaction applies to a lone carbon")
            .unwrap()
            .rhs()
            .clone();
        assert_eq!(rhs, expected);
    }

    // A molecule with two stereo centers — where one center's site is the other's ligand — must match
    // itself: `verify_overlays` selects the host stereo atom whose *site* is the mapped site, not the
    // first one merely incident to it. Regression for the two-distinct-site self-apply failure that the
    // stereo compose completeness surfaced.
    #[rstest]
    #[case::undetermined(StereoCosetAst::Undetermined)]
    fn test_reaction_ast_apply_two_stereo_centers(#[case] coset: StereoCosetAst) {
        let center = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, coset.clone()),
                ),
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, coset.clone()),
                ),
            ],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let rhs = ReactionAst::new(center.clone(), Deltas::new())
            .apply(&center, SubgraphIsomorphismAlgorithm::Vf2)
            .unwrap()
            .next()
            .expect("a two-stereo-center molecule matches itself")
            .unwrap()
            .rhs()
            .clone();
        assert_eq!(rhs, center);
    }

    #[rstest]
    fn test_reaction_ast_apply_at_comap() {
        // Remove atom O (id 1) and its bond: host C-O ⇒ product C. Atom 0 is preserved (mated), atom
        // 1 is deleted (left-exposed), so the comap's atom map records exactly that.
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
            ]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&[NodeId(0), NodeId(1)], host.atoms().count()),
        );
        let derivation = reaction.apply_at(&host, &correspondence).unwrap();
        assert_eq!(
            derivation.rhs(),
            &MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C)],
                bonds: vec![],
                ..Default::default()
            })
        );
        assert_eq!(derivation.atom_map().mates(), &[(NodeId(0), NodeId(0))]);
        assert_eq!(derivation.atom_map().left_exposed(), vec![NodeId(1)]);
    }
}
