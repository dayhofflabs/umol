//! Reaction AST: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations live in `reaction_span.rs`).

use std::collections::{BTreeMap, HashMap, HashSet};

use umol_graph_core::SubgraphIsomorphismAlgorithm;

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::dative::DativeBondAst;
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::stereo::StereoConfigurationAst;
use super::edit::{
    AddBond, AromaticSystemRef, AtomRef, BondRef, DativeBondRef, Edit, MulticenterBondRef,
    NoncovalentBondRef, StereoAtomFieldChange, StereoAtomRef, StereoAtomRemoval,
    StereoBondFieldChange, StereoBondRef, StereoBondRemoval,
};
use super::embedding::MoleculeEmbedding;
use super::error::{ApplyError, Contradiction};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::molecule::MoleculeAst;
use super::remap::IdRemapping;
use super::substructure::SubstructureMatchAlgorithm;
use super::traits::Canonicalize;

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

    /// Apply the reaction at one match `m` of `lhs` into a host (`m.ast()`), producing the
    /// transformed host. DPO: a deleted host atom must carry no localized bond the rule does not
    /// also delete (else `ApplyError::Dangling`). Created atoms/bonds are appended, preserved
    /// entities are mutated in place, deleted entities are removed (the host renumbers).
    /// Molecule-level constraints are added/removed with their entity refs re-anchored through the
    /// match (lhs → host, created → appended); transact's renumbering compacts them on removal.
    pub fn apply_at(&self, m: &MoleculeEmbedding) -> Result<MoleculeAst, ApplyError> {
        let deltas = self.deltas.clone().canonicalize()?;
        let host = m.ast();

        let mut created_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut created_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut sets: Vec<Edit> = Vec::new();
        let mut remove_atoms: Vec<AtomRef> = Vec::new();
        let mut remove_bonds: Vec<BondRef> = Vec::new();
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
                    let host_atom = m.host_atom(*id);
                    removed_host_atoms.push(host_atom);
                    remove_atoms.push(AtomRef::Id(host_atom));
                }
                Delta::Atom(AtomDelta::ModifyField { id, change }) => {
                    sets.push(Edit::ModifyAtomField {
                        id: AtomRef::Id(m.host_atom(*id)),
                        change: change.clone(),
                    })
                }
                Delta::Atom(AtomDelta::ModifyConstraint { id, old, new }) => {
                    sets.push(Edit::ModifyAtomConstraint {
                        id: AtomRef::Id(m.host_atom(*id)),
                        old: old.clone(),
                        new: new.clone(),
                    })
                }
                Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                    created_bonds.insert(*id, (*atoms, ast.clone()));
                }
                Delta::Bond(BondDelta::Remove { id, .. }) => {
                    let host_bond = m.host_bond(*id);
                    removed_host_bonds.insert(host_bond);
                    remove_bonds.push(BondRef::Id(host_bond));
                }
                Delta::Bond(BondDelta::ModifyField { id, change }) => {
                    sets.push(Edit::ModifyBondField {
                        id: BondRef::Id(m.host_bond(*id)),
                        change: change.clone(),
                    })
                }
                Delta::Bond(BondDelta::ModifyConstraint { id, old, new }) => {
                    sets.push(Edit::ModifyBondConstraint {
                        id: BondRef::Id(m.host_bond(*id)),
                        old: old.clone(),
                        new: new.clone(),
                    })
                }
                Delta::DativeBond(d) => match d {
                    DativeBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyDativeBondField {
                            id: DativeBondRef::Id(m.host_dative_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    DativeBondDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyDativeBondConstraint {
                            id: DativeBondRef::Id(m.host_dative_bond(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    DativeBondDelta::Add { .. } => {}
                    DativeBondDelta::Remove { id, .. } => {
                        removed_host_dative.insert(m.host_dative_bond(*id));
                    }
                },
                Delta::AromaticSystem(a) => match a {
                    AromaticSystemDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyAromaticSystemField {
                            id: AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                            change: change.clone(),
                        })
                    }
                    AromaticSystemDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyAromaticSystemConstraint {
                            id: AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    AromaticSystemDelta::Add { .. } => {}
                    AromaticSystemDelta::Remove { id, .. } => {
                        removed_host_aromatic.insert(m.host_aromatic_system(*id));
                    }
                },
                Delta::MulticenterBond(mc) => match mc {
                    MulticenterBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyMulticenterBondField {
                            id: MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    MulticenterBondDelta::ModifyConstraint { id, old, new } => {
                        sets.push(Edit::ModifyMulticenterBondConstraint {
                            id: MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    MulticenterBondDelta::Add { .. } => {}
                    MulticenterBondDelta::Remove { id, .. } => {
                        removed_host_multicenter.insert(m.host_multicenter_bond(*id));
                    }
                },
                Delta::NoncovalentBond(nc) => match nc {
                    NoncovalentBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyNoncovalentBondField {
                            id: NoncovalentBondRef::Id(m.host_noncovalent_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    // `NoncovalentBondConstraint` is uninhabited — no `Edit` variant, no-op.
                    NoncovalentBondDelta::ModifyConstraint { .. } => {}
                    NoncovalentBondDelta::Add { .. } => {}
                    NoncovalentBondDelta::Remove { id, .. } => {
                        removed_host_noncovalent.insert(m.host_noncovalent_bond(*id));
                    }
                },
                // Stereo: the four set-ops lower directly; the relative ops resolve against the
                // matched host config (same frame — no reindex, like the other overlays) and emit an
                // absolute `Configuration`. `Add` is lowered in the second pass; `Remove` tracks the
                // host id for the DPO dangling check.
                Delta::StereoAtom(s) => match s {
                    StereoAtomDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomRef::Id(m.host_stereo_atom(*id)),
                            change: change.clone(),
                        })
                    }
                    StereoAtomDelta::ModifyConstraint { id, old, new, .. } => {
                        sets.push(Edit::ModifyStereoAtomConstraint {
                            id: StereoAtomRef::Id(m.host_stereo_atom(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    StereoAtomDelta::Apply { id, kind, permutation } => {
                        let host_id = m.host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.apply(*permutation);
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomRef::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Swap { id, kind } => {
                        let host_id = m.host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.swap();
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomRef::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Mirror { id, kind } => {
                        let host_id = m.host_stereo_atom(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_atom(host_id).coset().clone(),
                        );
                        let new = old.mirror();
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomRef::Id(host_id),
                            change: StereoAtomFieldChange::Configuration { old, new },
                        })
                    }
                    StereoAtomDelta::Add { .. } => {}
                    StereoAtomDelta::Remove { id, .. } => {
                        removed_host_stereo_atom.insert(m.host_stereo_atom(*id));
                    }
                },
                Delta::StereoBond(s) => match s {
                    StereoBondDelta::ModifyField { id, change } => {
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondRef::Id(m.host_stereo_bond(*id)),
                            change: change.clone(),
                        })
                    }
                    StereoBondDelta::ModifyConstraint { id, old, new, .. } => {
                        sets.push(Edit::ModifyStereoBondConstraint {
                            id: StereoBondRef::Id(m.host_stereo_bond(*id)),
                            old: old.clone(),
                            new: new.clone(),
                        })
                    }
                    StereoBondDelta::Apply { id, kind, permutation } => {
                        let host_id = m.host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.apply(*permutation);
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondRef::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Swap { id, kind } => {
                        let host_id = m.host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.swap();
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondRef::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Mirror { id, kind } => {
                        let host_id = m.host_stereo_bond(*id);
                        let old = StereoConfigurationAst::Kinded(
                            *kind,
                            host.stereo_bond(host_id).coset().clone(),
                        );
                        let new = old.mirror();
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondRef::Id(host_id),
                            change: StereoBondFieldChange::Configuration { old, new },
                        })
                    }
                    StereoBondDelta::Add { .. } => {}
                    StereoBondDelta::Remove { id, .. } => {
                        removed_host_stereo_bond.insert(m.host_stereo_bond(*id));
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
        let new_index: HashMap<AtomId, usize> = created_atoms
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();
        let atom_ref = |id: AtomId| match new_index.get(&id) {
            Some(&index) => AtomRef::New(index),
            None => AtomRef::Id(m.host_atom(id)),
        };
        let new_bond_index: HashMap<BondId, usize> = created_bonds
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();
        let bond_ref = |id: BondId| match new_bond_index.get(&id) {
            Some(&index) => BondRef::New(index),
            None => BondRef::Id(m.host_bond(id)),
        };

        // Overlay create/remove need `atom_ref` (created participants resolve to `New`), so they
        // are lowered in a second pass: adds after the topology adds, removes before
        // `RemoveTopology`. Removes are collected per kind and emitted as one batched edit each,
        // so each overlay id space is compacted once against the pre-removal state (a sequence of
        // single-id removes would stale the not-yet-processed ids). Dative `atoms` is
        // `[donors…, acceptor]` (acceptor last, per transact).
        let mut overlay_adds: Vec<Edit> = Vec::new();
        let mut remove_dative: Vec<(DativeBondRef, Vec<AtomRef>, DativeBondAst)> = Vec::new();
        let mut remove_aromatic: Vec<(AromaticSystemRef, Vec<AtomRef>, AromaticSystemAst)> =
            Vec::new();
        let mut remove_multicenter: Vec<(MulticenterBondRef, Vec<AtomRef>, MulticenterBondAst)> =
            Vec::new();
        let mut remove_noncovalent: Vec<(NoncovalentBondRef, [AtomRef; 2], NoncovalentBondAst)> =
            Vec::new();
        let mut remove_stereo_atom: Vec<StereoAtomRemoval> = Vec::new();
        let mut remove_stereo_bond: Vec<StereoBondRemoval> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::DativeBond(DativeBondDelta::Add {
                    donors, acceptor, ast, ..
                }) => {
                    let mut atoms: Vec<AtomRef> = donors.iter().map(|a| atom_ref(*a)).collect();
                    atoms.push(atom_ref(*acceptor));
                    overlay_adds.push(Edit::AddDativeBond {
                        atoms,
                        ast: ast.clone(),
                    });
                }
                Delta::DativeBond(DativeBondDelta::Remove {
                    id, donors, acceptor, ast,
                }) => {
                    let mut atoms: Vec<AtomRef> = donors.iter().map(|a| atom_ref(*a)).collect();
                    atoms.push(atom_ref(*acceptor));
                    remove_dative.push((
                        DativeBondRef::Id(m.host_dative_bond(*id)),
                        atoms,
                        ast.clone(),
                    ));
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddAromaticSystem {
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, ast }) => {
                    remove_aromatic.push((
                        AromaticSystemRef::Id(m.host_aromatic_system(*id)),
                        atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast.clone(),
                    ));
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddMulticenterBond {
                        atoms: atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, ast }) => {
                    remove_multicenter.push((
                        MulticenterBondRef::Id(m.host_multicenter_bond(*id)),
                        atoms.iter().map(|a| atom_ref(*a)).collect(),
                        ast.clone(),
                    ));
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, ast, .. }) => {
                    overlay_adds.push(Edit::AddNoncovalentBond {
                        atoms: [atom_ref(atoms[0]), atom_ref(atoms[1])],
                        ast: ast.clone(),
                    });
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, ast }) => {
                    remove_noncovalent.push((
                        NoncovalentBondRef::Id(m.host_noncovalent_bond(*id)),
                        [atom_ref(atoms[0]), atom_ref(atoms[1])],
                        ast.clone(),
                    ));
                }
                Delta::StereoAtom(StereoAtomDelta::Add {
                    site, ligands, ast, ..
                }) => {
                    overlay_adds.push(Edit::AddStereoAtom {
                        site: atom_ref(*site),
                        ligands: ligands.iter().map(|l| (atom_ref(l.atom_id), l.kind)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::StereoAtom(StereoAtomDelta::Remove {
                    id, site, ligands, ast,
                }) => {
                    remove_stereo_atom.push((
                        StereoAtomRef::Id(m.host_stereo_atom(*id)),
                        atom_ref(*site),
                        ligands.iter().map(|l| (atom_ref(l.atom_id), l.kind)).collect(),
                        ast.clone(),
                    ));
                }
                Delta::StereoBond(StereoBondDelta::Add {
                    site, ligands, ast, ..
                }) => {
                    overlay_adds.push(Edit::AddStereoBond {
                        site: bond_ref(*site),
                        ligands: ligands.iter().map(|l| (atom_ref(l.atom_id), l.kind)).collect(),
                        ast: ast.clone(),
                    });
                }
                Delta::StereoBond(StereoBondDelta::Remove {
                    id, site, ligands, ast,
                }) => {
                    remove_stereo_bond.push((
                        StereoBondRef::Id(m.host_stereo_bond(*id)),
                        bond_ref(*site),
                        ligands.iter().map(|l| (atom_ref(l.atom_id), l.kind)).collect(),
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
                .map(|i| (AtomId(i), m.host_atom(AtomId(i))))
                .collect();
            for (&created, &index) in &new_index {
                atom.insert(created, AtomId((host_atom_count + index) as u32));
            }
            let mut bond: HashMap<BondId, BondId> = (0..self.lhs.bonds().count() as u32)
                .map(|i| (BondId(i), m.host_bond(BondId(i))))
                .collect();
            for (index, &created) in created_bonds.keys().enumerate() {
                bond.insert(created, BondId((host_bond_count + index) as u32));
            }
            let match_map = IdRemapping::new(
                atom,
                bond,
                (0..self.lhs.dative_bonds().count() as u32)
                    .map(|i| (DativeBondId(i), m.host_dative_bond(DativeBondId(i))))
                    .collect(),
                (0..self.lhs.aromatic_systems().count() as u32)
                    .map(|i| (AromaticSystemId(i), m.host_aromatic_system(AromaticSystemId(i))))
                    .collect(),
                (0..self.lhs.multicenter_bonds().count() as u32)
                    .map(|i| {
                        (
                            MulticenterBondId(i),
                            m.host_multicenter_bond(MulticenterBondId(i)),
                        )
                    })
                    .collect(),
                (0..self.lhs.noncovalent_bonds().count() as u32)
                    .map(|i| {
                        (
                            NoncovalentBondId(i),
                            m.host_noncovalent_bond(NoncovalentBondId(i)),
                        )
                    })
                    .collect(),
                (0..self.lhs.stereo_atoms().count() as u32)
                    .map(|i| (StereoAtomId(i), m.host_stereo_atom(StereoAtomId(i))))
                    .collect(),
                (0..self.lhs.stereo_bonds().count() as u32)
                    .map(|i| (StereoBondId(i), m.host_stereo_bond(StereoBondId(i))))
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
                        endpoints: [atom_ref(atoms[0]), atom_ref(atoms[1])],
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
        Ok(builder.build())
    }

    /// Every product of applying the reaction to `host`: one per injective match of `lhs` into
    /// `host` (via `subiso`) that satisfies the DPO gluing condition. Matches that dangle are
    /// skipped.
    pub fn apply<'h>(
        &'h self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> impl Iterator<Item = MoleculeAst> + 'h {
        self.lhs
            .substructure_matches(host, SubstructureMatchAlgorithm::GraphAndOverlays, subiso)
            .into_iter()
            .filter_map(move |m| self.apply_at(&m).ok())
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
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
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
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        vec![AtomId(0), AtomId(1)],
        vec![],
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
        ),
    )]
    #[case::overlay_removed(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
                vec![], vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::O) }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            vec![], vec![], vec![], vec![],
            vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            vec![], vec![],
            Constraints::new(),
        ),
        vec![AtomId(0), AtomId(1)],
        vec![NoncovalentBondId(0)],
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::O)], vec![]),
    )]
    fn test_reaction_ast_apply_at(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] atom_map: Vec<AtomId>,
        #[case] host_noncovalent: Vec<NoncovalentBondId>,
        #[case] expected: MoleculeAst,
    ) {
        let embedding = MoleculeEmbedding::from_match(
            &host,
            &reaction.lhs,
            atom_map,
            vec![],
            vec![],
            vec![],
            host_noncovalent,
            vec![],
            vec![],
        );
        assert_eq!(reaction.apply_at(&embedding).unwrap(), expected);
    }

    #[rstest]
    #[case::dangling_bond(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        ),
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
    )]
    #[case::dangling_noncovalent(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::O)], vec![]),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::O),
            })]),
        ),
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            vec![], vec![], vec![], vec![],
            vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            vec![], vec![],
            Constraints::new(),
        ),
    )]
    fn test_reaction_ast_apply_at_error(#[case] reaction: ReactionAst, #[case] host: MoleculeAst) {
        // The rule deletes a host atom that still carries an undeleted bond/overlay → dangling.
        let embedding = MoleculeEmbedding::from_match(
            &host,
            &reaction.lhs,
            vec![AtomId(0)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(
            reaction.apply_at(&embedding),
            Err(ApplyError::Dangling {
                host_atom: AtomId(0),
            }),
        );
    }

    #[rstest]
    fn test_reaction_ast_apply_at_molecule_constraint() {
        // A reaction adding a molecule-level `ChargeSum` over its lhs atoms; applied at a match
        // that maps lhs atoms 0,1 → host atoms 1,2, the constraint's refs re-anchor to the host.
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
                MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                    sum: ValueAst::Lit(0),
                },
            )))]),
        );
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
        );
        let embedding = MoleculeEmbedding::from_match(
            &host,
            &reaction.lhs,
            vec![AtomId(1), AtomId(2)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let result = reaction.apply_at(&embedding).unwrap();
        assert_eq!(
            result.constraints(),
            &Constraints::from(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
                sum: ValueAst::Lit(0),
            })),
        );
    }

    #[rstest]
    fn test_reaction_ast_apply() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .collect();
        assert_eq!(
            products,
            vec![MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            )],
        );
    }
}
