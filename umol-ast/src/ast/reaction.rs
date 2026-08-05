//! Reaction AST: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations live in `reaction_span.rs`).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::iter::from_fn;

use umol_graph_core::{Correspondence, NodeId};
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
    BondFieldChange, BondHandle, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit,
    Edits, EntityHandle, MulticenterBondFieldChange, MulticenterBondHandle,
    NoncovalentBondFieldChange, NoncovalentBondHandle, StereoAtomFieldChange, StereoAtomHandle,
    StereoAtomRemoval, StereoBondFieldChange, StereoBondHandle, StereoBondRemoval,
};
use super::entity::Entity;
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
use super::stereo::{StereoConfigurationAst, StereoCoset, StereoKind, StereoTerm};
use super::substructure::SubstructureMatchConfig;
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

fn stereo_delta_domains_are_valid(lhs: &MoleculeAst, deltas: &Deltas) -> bool {
    fn permutation_is_valid(kind: StereoKind, permutation: Permutation) -> bool {
        permutation.degree() == kind.degree()
            && (0..kind.count() as u32).all(|coset| {
                kind.class_key()
                    .space()
                    .reindex(coset, permutation)
                    .is_some()
            })
    }

    fn term_is_valid(kind: StereoKind, term: &StereoTerm) -> bool {
        match term {
            StereoTerm::Var(var) => var.1.as_ref().is_none_or(|domain| {
                !domain.is_empty() && domain.iter().all(|i| *i < kind.count() as u32)
            }),
            StereoTerm::Lit(index) => *index < kind.count() as u32,
            StereoTerm::LitSet(indices) => {
                !indices.is_empty() && indices.iter().all(|i| *i < kind.count() as u32)
            }
            StereoTerm::Swap(inner) | StereoTerm::Mirror(inner) => term_is_valid(kind, inner),
            StereoTerm::Apply(inner, permutation) => {
                term_is_valid(kind, inner) && permutation_is_valid(kind, *permutation)
            }
        }
    }

    fn configuration_is_valid(configuration: &StereoConfigurationAst) -> bool {
        match configuration {
            StereoConfigurationAst::Undetermined => true,
            StereoConfigurationAst::Kinded(kind, coset) => match coset {
                StereoCoset::Undetermined => true,
                StereoCoset::Lit(index) => *index < kind.count() as u32,
                StereoCoset::LitSet(indices) => {
                    !indices.is_empty() && indices.iter().all(|i| *i < kind.count() as u32)
                }
                StereoCoset::Term(term) => term_is_valid(*kind, term),
            },
        }
    }

    fn configurations_are_compatible(
        lhs: Option<&StereoConfigurationAst>,
        old: &StereoConfigurationAst,
        new: &StereoConfigurationAst,
    ) -> bool {
        let lhs_kind = lhs.and_then(StereoConfigurationAst::kind);
        let old_kind = old.kind();
        let new_kind = new.kind();
        configuration_is_valid(old)
            && configuration_is_valid(new)
            && old_kind.zip(new_kind).is_none_or(|(old, new)| old == new)
            && lhs_kind.zip(old_kind).is_none_or(|(lhs, old)| lhs == old)
            && lhs_kind.zip(new_kind).is_none_or(|(lhs, new)| lhs == new)
    }

    if lhs
        .stereo_atoms()
        .iter()
        .any(|view| !configuration_is_valid(&view.ast.configuration))
        || lhs
            .stereo_bonds()
            .iter()
            .any(|view| !configuration_is_valid(&view.ast.configuration))
    {
        return false;
    }

    deltas.iter().all(|delta| match delta {
        Delta::StereoAtom(StereoAtomDelta::Add { ast, .. }) => {
            configuration_is_valid(&ast.configuration)
        }
        Delta::StereoAtom(StereoAtomDelta::Remove { id, ast, .. }) => {
            let lhs_configuration = lhs
                .stereo_atoms()
                .get(*id)
                .map(|view| &view.ast.configuration);
            configurations_are_compatible(lhs_configuration, &ast.configuration, &ast.configuration)
        }
        Delta::StereoAtom(StereoAtomDelta::ModifyField {
            id,
            change: StereoAtomFieldChange::Configuration { old, new },
        }) => configurations_are_compatible(
            lhs.stereo_atoms()
                .get(*id)
                .map(|view| &view.ast.configuration),
            old,
            new,
        ),
        Delta::StereoAtom(StereoAtomDelta::Apply {
            id,
            kind,
            permutation,
        }) => {
            lhs.stereo_atoms()
                .get(*id)
                .and_then(|view| view.ast.configuration.kind())
                .is_none_or(|lhs_kind| lhs_kind == *kind)
                && permutation_is_valid(*kind, *permutation)
        }
        Delta::StereoAtom(
            StereoAtomDelta::Swap { id, kind } | StereoAtomDelta::Mirror { id, kind },
        ) => lhs
            .stereo_atoms()
            .get(*id)
            .and_then(|view| view.ast.configuration.kind())
            .is_none_or(|lhs_kind| lhs_kind == *kind),
        Delta::StereoBond(StereoBondDelta::Add { ast, .. }) => {
            configuration_is_valid(&ast.configuration)
        }
        Delta::StereoBond(StereoBondDelta::Remove { id, ast, .. }) => {
            let lhs_configuration = lhs
                .stereo_bonds()
                .get(*id)
                .map(|view| &view.ast.configuration);
            configurations_are_compatible(lhs_configuration, &ast.configuration, &ast.configuration)
        }
        Delta::StereoBond(StereoBondDelta::ModifyField {
            id,
            change: StereoBondFieldChange::Configuration { old, new },
        }) => configurations_are_compatible(
            lhs.stereo_bonds()
                .get(*id)
                .map(|view| &view.ast.configuration),
            old,
            new,
        ),
        Delta::StereoBond(StereoBondDelta::Apply {
            id,
            kind,
            permutation,
        }) => {
            lhs.stereo_bonds()
                .get(*id)
                .and_then(|view| view.ast.configuration.kind())
                .is_none_or(|lhs_kind| lhs_kind == *kind)
                && permutation_is_valid(*kind, *permutation)
        }
        Delta::StereoBond(
            StereoBondDelta::Swap { id, kind } | StereoBondDelta::Mirror { id, kind },
        ) => lhs
            .stereo_bonds()
            .get(*id)
            .and_then(|view| view.ast.configuration.kind())
            .is_none_or(|lhs_kind| lhs_kind == *kind),
        _ => true,
    })
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
    /// appended); transact's renumbering compacts them on removal. The supplied correspondence must
    /// be total on the pattern and agree with the mapped topology, overlay incidence, and stereo
    /// sites; incompatible stereo ligand frames are reported separately.
    pub fn apply_at(
        &self,
        host: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> Result<ReactionDerivation, ApplyError> {
        if correspondence.atoms().left_count() != self.lhs.atoms().count()
            || correspondence.atoms().right_count() != host.atoms().count()
            || correspondence.atoms().matched_pair_count() != self.lhs.atoms().count()
        {
            return Err(ApplyError::CorrespondenceMismatch {
                entity: Entity::Atom(AtomId(0)),
            });
        }
        let mut host_atoms = HashSet::new();
        for left in self.lhs.atoms().ids() {
            let entity = Entity::Atom(left);
            let Some(right) = correspondence.atoms().right_of(NodeId::from(left)) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            if right.index() >= host.atoms().count() || !host_atoms.insert(right) {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
        }

        let induced =
            MoleculeCorrespondence::induce(&self.lhs, host, correspondence.atoms().clone());
        macro_rules! require_induced_family {
            ($family:ident, $entity:ident, $fallback:expr) => {{
                let supplied = correspondence.$family();
                let expected = induced.$family();
                if supplied != expected {
                    let id = supplied
                        .matched_pairs()
                        .iter()
                        .find_map(|&(left, right)| {
                            (expected.right_of(left) != Some(right)).then_some(left)
                        })
                        .or_else(|| {
                            expected.matched_pairs().iter().find_map(|&(left, right)| {
                                (supplied.right_of(left) != Some(right)).then_some(left)
                            })
                        })
                        .unwrap_or($fallback);
                    return Err(ApplyError::CorrespondenceMismatch {
                        entity: Entity::$entity(id),
                    });
                }
            }};
        }
        require_induced_family!(bonds, Bond, BondId(0));
        require_induced_family!(dative_bonds, DativeBond, DativeBondId(0));
        require_induced_family!(aromatic_systems, AromaticSystem, AromaticSystemId(0));
        require_induced_family!(multicenter_bonds, MulticenterBond, MulticenterBondId(0));
        require_induced_family!(noncovalent_bonds, NoncovalentBond, NoncovalentBondId(0));

        if correspondence.stereo_atoms().left_count() != self.lhs.stereo_atoms().count()
            || correspondence.stereo_atoms().right_count() != host.stereo_atoms().count()
            || correspondence.stereo_atoms().matched_pair_count() != self.lhs.stereo_atoms().count()
        {
            return Err(ApplyError::CorrespondenceMismatch {
                entity: Entity::StereoAtom(StereoAtomId(0)),
            });
        }
        let mut host_stereo_atoms = HashSet::new();
        for left in self.lhs.stereo_atoms().ids() {
            let entity = Entity::StereoAtom(left);
            let Some(right) = correspondence.stereo_atoms().right_of(left) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            let Some(host_view) = host.stereo_atoms().get(right) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            if !host_stereo_atoms.insert(right) {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
            let rule_site = self.lhs.stereo_atom(left).site_id();
            if correspondence
                .atoms()
                .right_of(NodeId::from(rule_site))
                .map(AtomId::from)
                != Some(host_view.site_id())
            {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
        }
        if correspondence.stereo_bonds().left_count() != self.lhs.stereo_bonds().count()
            || correspondence.stereo_bonds().right_count() != host.stereo_bonds().count()
            || correspondence.stereo_bonds().matched_pair_count() != self.lhs.stereo_bonds().count()
        {
            return Err(ApplyError::CorrespondenceMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            });
        }
        let mut host_stereo_bonds = HashSet::new();
        for left in self.lhs.stereo_bonds().ids() {
            let entity = Entity::StereoBond(left);
            let Some(right) = correspondence.stereo_bonds().right_of(left) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            let Some(host_view) = host.stereo_bonds().get(right) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            if !host_stereo_bonds.insert(right) {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
            let rule_site = self.lhs.stereo_bond(left).site_id();
            if correspondence.bonds().right_of(rule_site) != Some(host_view.site_id()) {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
        }

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
        reframe_stereo(&mut deltas, &self.lhs, host, correspondence)?;
        let host_atom = |id: AtomId| {
            correspondence
                .atoms()
                .right_of(NodeId::from(id))
                .map(AtomId::from)
                .filter(|right| right.index() < host.atoms().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::Atom(id),
                })
        };
        let host_bond = |id: BondId| {
            correspondence
                .bonds()
                .right_of(id)
                .filter(|right| right.index() < host.bonds().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::Bond(id),
                })
        };
        let host_dative = |id: DativeBondId| {
            correspondence
                .dative_bonds()
                .right_of(id)
                .filter(|right| right.index() < host.dative_bonds().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::DativeBond(id),
                })
        };
        let host_aromatic = |id: AromaticSystemId| {
            correspondence
                .aromatic_systems()
                .right_of(id)
                .filter(|right| right.index() < host.aromatic_systems().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::AromaticSystem(id),
                })
        };
        let host_multicenter = |id: MulticenterBondId| {
            correspondence
                .multicenter_bonds()
                .right_of(id)
                .filter(|right| right.index() < host.multicenter_bonds().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::MulticenterBond(id),
                })
        };
        let host_noncovalent = |id: NoncovalentBondId| {
            correspondence
                .noncovalent_bonds()
                .right_of(id)
                .filter(|right| right.index() < host.noncovalent_bonds().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::NoncovalentBond(id),
                })
        };
        let host_stereo_atom = |id: StereoAtomId| {
            correspondence
                .stereo_atoms()
                .right_of(id)
                .filter(|right| right.index() < host.stereo_atoms().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::StereoAtom(id),
                })
        };
        let host_stereo_bond = |id: StereoBondId| {
            correspondence
                .stereo_bonds()
                .right_of(id)
                .filter(|right| right.index() < host.stereo_bonds().count())
                .ok_or(ApplyError::CorrespondenceMismatch {
                    entity: Entity::StereoBond(id),
                })
        };

        let mut created_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut created_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut sets = Edits::new();
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
                    let removed = host_atom(*id)?;
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
                        AtomFieldChange::UnpairedElectrons { old, new } => AtomUpdate {
                            unpaired_electrons: old.difference_to(new),
                            ..Default::default()
                        },
                    };
                    let host_id = host_atom(*id)?;
                    sets.update_atom(AtomHandle::Id(host_id), host.atom(host_id).ast, &update);
                }
                Delta::Atom(AtomDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_atom(*id)?;
                        sets.update_atom(
                            AtomHandle::Id(host_id),
                            host.atom(host_id).ast,
                            &AtomUpdate {
                                constraints: constraint.into(),
                                ..Default::default()
                            },
                        );
                    }
                }
                Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                    created_bonds.insert(*id, (*atoms, ast.clone()));
                }
                Delta::Bond(BondDelta::Remove { id, .. }) => {
                    let removed = host_bond(*id)?;
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
                        BondFieldChange::UnpairedElectrons { old, new } => BondUpdate {
                            unpaired_electrons: old.difference_to(new),
                            ..Default::default()
                        },
                    };
                    let host_id = host_bond(*id)?;
                    sets.update_bond(BondHandle::Id(host_id), host.bond(host_id).ast, &update);
                }
                Delta::Bond(BondDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_bond(*id)?;
                        sets.update_bond(
                            BondHandle::Id(host_id),
                            host.bond(host_id).ast,
                            &BondUpdate {
                                constraints: constraint.into(),
                                ..Default::default()
                            },
                        );
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
                        let host_id = host_dative(*id)?;
                        sets.update_dative_bond(
                            DativeBondHandle::Id(host_id),
                            host.dative_bond(host_id).ast,
                            &update,
                        );
                    }
                    DativeBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_dative(*id)?;
                            sets.update_dative_bond(
                                DativeBondHandle::Id(host_id),
                                host.dative_bond(host_id).ast,
                                &DativeBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    DativeBondDelta::Add { .. } => {}
                    DativeBondDelta::Remove { id, .. } => {
                        removed_host_dative.insert(host_dative(*id)?);
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
                            AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                                AromaticSystemUpdate {
                                    unpaired_electrons: old.difference_to(new),
                                    ..Default::default()
                                }
                            }
                        };
                        let host_id = host_aromatic(*id)?;
                        sets.update_aromatic_system(
                            AromaticSystemHandle::Id(host_id),
                            host.aromatic_system(host_id).ast,
                            &update,
                        );
                    }
                    AromaticSystemDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_aromatic(*id)?;
                            sets.update_aromatic_system(
                                AromaticSystemHandle::Id(host_id),
                                host.aromatic_system(host_id).ast,
                                &AromaticSystemUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    AromaticSystemDelta::Add { .. } => {}
                    AromaticSystemDelta::Remove { id, .. } => {
                        removed_host_aromatic.insert(host_aromatic(*id)?);
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
                            MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                                MulticenterBondUpdate {
                                    unpaired_electrons: old.difference_to(new),
                                    ..Default::default()
                                }
                            }
                        };
                        let host_id = host_multicenter(*id)?;
                        sets.update_multicenter_bond(
                            MulticenterBondHandle::Id(host_id),
                            host.multicenter_bond(host_id).ast,
                            &update,
                        );
                    }
                    MulticenterBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_multicenter(*id)?;
                            sets.update_multicenter_bond(
                                MulticenterBondHandle::Id(host_id),
                                host.multicenter_bond(host_id).ast,
                                &MulticenterBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    MulticenterBondDelta::Add { .. } => {}
                    MulticenterBondDelta::Remove { id, .. } => {
                        removed_host_multicenter.insert(host_multicenter(*id)?);
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
                        let host_id = host_noncovalent(*id)?;
                        sets.update_noncovalent_bond(
                            NoncovalentBondHandle::Id(host_id),
                            host.noncovalent_bond(host_id).ast,
                            &update,
                        );
                    }
                    NoncovalentBondDelta::ModifyConstraint { id, old, new } => {
                        let constraint = new.clone().or_else(|| {
                            old.as_ref().map(|constraint| constraint.as_undetermined())
                        });
                        if let Some(constraint) = constraint {
                            let host_id = host_noncovalent(*id)?;
                            sets.update_noncovalent_bond(
                                NoncovalentBondHandle::Id(host_id),
                                host.noncovalent_bond(host_id).ast,
                                &NoncovalentBondUpdate {
                                    constraints: constraint.into(),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    NoncovalentBondDelta::Add { .. } => {}
                    NoncovalentBondDelta::Remove { id, .. } => {
                        removed_host_noncovalent.insert(host_noncovalent(*id)?);
                    }
                },
                // Stereo: the four set-ops lower directly; the relative ops resolve against the
                // matched host config (same frame — no reindex, like the other overlays) and emit an
                // absolute `Configuration`. `Add` is lowered in the second pass; `Remove` tracks the
                // host id for the DPO dangling check.
                Delta::StereoAtom(s) => match s {
                    StereoAtomDelta::ModifyField { id, change } => {
                        let host_id = host_stereo_atom(*id)?;
                        let StereoAtomFieldChange::Configuration { new, .. } = change;
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration {
                                old: host.stereo_atom(host_id).ast.configuration.clone(),
                                new: new.clone(),
                            },
                        })
                    }
                    StereoAtomDelta::ModifyConstraint { id, kind, old, new } => {
                        if let Some(constraint) = new.as_ref().or(old.as_ref()) {
                            let host_id = host_stereo_atom(*id)?;
                            sets.push(Edit::ModifyStereoAtomConstraint {
                                id: StereoAtomHandle::Id(host_id),
                                kind: *kind,
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
                        let host_id = host_stereo_atom(*id)?;
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
                        let host_id = host_stereo_atom(*id)?;
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
                        let host_id = host_stereo_atom(*id)?;
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
                        removed_host_stereo_atom.insert(host_stereo_atom(*id)?);
                    }
                },
                Delta::StereoBond(s) => match s {
                    StereoBondDelta::ModifyField { id, change } => {
                        let host_id = host_stereo_bond(*id)?;
                        let StereoBondFieldChange::Configuration { new, .. } = change;
                        sets.push(Edit::ModifyStereoBondField {
                            id: StereoBondHandle::Id(host_id),
                            change: StereoBondFieldChange::Configuration {
                                old: host.stereo_bond(host_id).ast.configuration.clone(),
                                new: new.clone(),
                            },
                        })
                    }
                    StereoBondDelta::ModifyConstraint { id, kind, old, new } => {
                        if let Some(constraint) = new.as_ref().or(old.as_ref()) {
                            let host_id = host_stereo_bond(*id)?;
                            sets.push(Edit::ModifyStereoBondConstraint {
                                id: StereoBondHandle::Id(host_id),
                                kind: *kind,
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
                        let host_id = host_stereo_bond(*id)?;
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
                        let host_id = host_stereo_bond(*id)?;
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
                        let host_id = host_stereo_bond(*id)?;
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
                        removed_host_stereo_bond.insert(host_stereo_bond(*id)?);
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
            Some(&index) => Ok(AtomHandle::New(index)),
            None => host_atom(id).map(AtomHandle::Id),
        };
        // Created bonds use their own `New(0..k)` namespace in ascending id order.
        let new_bond_index: HashMap<BondId, usize> = created_bonds
            .keys()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();
        let bond_handle = |id: BondId| match new_bond_index.get(&id) {
            Some(&index) => Ok(BondHandle::New(index)),
            None => host_bond(id).map(BondHandle::Id),
        };

        // Overlay create/remove need `atom_handle` (created participants resolve to `New`), so they
        // are lowered in a second pass: adds after the topology adds, removes before
        // `RemoveTopology`. Removes are collected per kind and emitted as one batched edit each,
        // so each overlay id space is compacted once against the pre-removal state (a sequence of
        // single-id removes would stale the not-yet-processed ids). Dative `atoms` is
        // `[donors…, acceptor]` (acceptor last, per transact).
        let mut overlay_adds = Edits::new();
        let mut new_dative_handles = HashMap::new();
        let mut new_aromatic_handles = HashMap::new();
        let mut new_multicenter_handles = HashMap::new();
        let mut new_noncovalent_handles = HashMap::new();
        let mut new_stereo_atom_handles = HashMap::new();
        let mut new_stereo_bond_handles = HashMap::new();
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
                    id,
                    donors,
                    acceptor,
                    ast,
                }) => {
                    let mut atoms: Vec<AtomHandle> = donors
                        .iter()
                        .map(|a| atom_handle(*a))
                        .collect::<Result<_, _>>()?;
                    atoms.push(atom_handle(*acceptor)?);
                    let handle = overlay_adds.add_dative_bond(atoms, ast.clone());
                    new_dative_handles.insert(*id, handle);
                }
                Delta::DativeBond(DativeBondDelta::Remove {
                    id,
                    donors,
                    acceptor,
                    ast,
                }) => {
                    let mut atoms: Vec<AtomHandle> = donors
                        .iter()
                        .map(|a| atom_handle(*a))
                        .collect::<Result<_, _>>()?;
                    atoms.push(atom_handle(*acceptor)?);
                    remove_dative.push((
                        DativeBondHandle::Id(host_dative(*id)?),
                        atoms,
                        ast.clone(),
                    ));
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, ast }) => {
                    let handle = overlay_adds.add_aromatic_system(
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    );
                    new_aromatic_handles.insert(*id, handle);
                }
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, ast }) => {
                    remove_aromatic.push((
                        AromaticSystemHandle::Id(host_aromatic(*id)?),
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    ));
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, ast }) => {
                    let handle = overlay_adds.add_multicenter_bond(
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    );
                    new_multicenter_handles.insert(*id, handle);
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, ast }) => {
                    remove_multicenter.push((
                        MulticenterBondHandle::Id(host_multicenter(*id)?),
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    ));
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, ast }) => {
                    let handle = overlay_adds.add_noncovalent_bond(
                        [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                        ast.clone(),
                    );
                    new_noncovalent_handles.insert(*id, handle);
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, ast }) => {
                    remove_noncovalent.push((
                        NoncovalentBondHandle::Id(host_noncovalent(*id)?),
                        [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                        ast.clone(),
                    ));
                }
                Delta::StereoAtom(StereoAtomDelta::Add {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    let handle = overlay_adds.add_stereo_atom(
                        atom_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    );
                    new_stereo_atom_handles.insert(*id, handle);
                }
                Delta::StereoAtom(StereoAtomDelta::Remove {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    remove_stereo_atom.push((
                        StereoAtomHandle::Id(host_stereo_atom(*id)?),
                        atom_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    ));
                }
                Delta::StereoBond(StereoBondDelta::Add {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    let handle = overlay_adds.add_stereo_bond(
                        bond_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    );
                    new_stereo_bond_handles.insert(*id, handle);
                }
                Delta::StereoBond(StereoBondDelta::Remove {
                    id,
                    site,
                    ligands,
                    ast,
                }) => {
                    remove_stereo_bond.push((
                        StereoBondHandle::Id(host_stereo_bond(*id)?),
                        bond_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        ast.clone(),
                    ));
                }
                _ => {}
            }
        }
        let mut overlay_removes = Edits::new();
        if !remove_dative.is_empty() {
            overlay_removes.remove_dative_bonds(remove_dative);
        }
        if !remove_aromatic.is_empty() {
            overlay_removes.remove_aromatic_systems(remove_aromatic);
        }
        if !remove_multicenter.is_empty() {
            overlay_removes.remove_multicenter_bonds(remove_multicenter);
        }
        if !remove_noncovalent.is_empty() {
            overlay_removes.remove_noncovalent_bonds(remove_noncovalent);
        }
        if !remove_stereo_atom.is_empty() {
            overlay_removes.remove_stereo_atoms(remove_stereo_atom);
        }
        if !remove_stereo_bond.is_empty() {
            overlay_removes.remove_stereo_bonds(remove_stereo_bond);
        }

        // Molecule-level constraints retain stable handles until transaction application. Preserved
        // LHS entities map to host `Id` handles and every created entity maps to the `New` handle
        // issued by its scheduled addition. Constraints precede all removals so each removal's
        // compaction updates surviving references and drops references to deleted entities.
        let mut constraint_edits = Edits::new();
        if !constraint_deltas.is_empty() {
            let mut handles = HashMap::new();
            for left in self.lhs.atoms().ids() {
                handles.insert(
                    Entity::Atom(left),
                    EntityHandle::Atom(AtomHandle::Id(host_atom(left)?)),
                );
            }
            for (&created, &index) in &new_atom_index {
                handles.insert(
                    Entity::Atom(created),
                    EntityHandle::Atom(AtomHandle::New(index)),
                );
            }
            for left in self.lhs.bonds().ids() {
                handles.insert(
                    Entity::Bond(left),
                    EntityHandle::Bond(BondHandle::Id(host_bond(left)?)),
                );
            }
            for (&created, &index) in &new_bond_index {
                handles.insert(
                    Entity::Bond(created),
                    EntityHandle::Bond(BondHandle::New(index)),
                );
            }
            for left in self.lhs.dative_bonds().ids() {
                handles.insert(
                    Entity::DativeBond(left),
                    EntityHandle::DativeBond(DativeBondHandle::Id(host_dative(left)?)),
                );
            }
            for (created, handle) in new_dative_handles {
                handles.insert(
                    Entity::DativeBond(created),
                    EntityHandle::DativeBond(handle),
                );
            }
            for left in self.lhs.aromatic_systems().ids() {
                handles.insert(
                    Entity::AromaticSystem(left),
                    EntityHandle::AromaticSystem(AromaticSystemHandle::Id(host_aromatic(left)?)),
                );
            }
            for (created, handle) in new_aromatic_handles {
                handles.insert(
                    Entity::AromaticSystem(created),
                    EntityHandle::AromaticSystem(handle),
                );
            }
            for left in self.lhs.multicenter_bonds().ids() {
                handles.insert(
                    Entity::MulticenterBond(left),
                    EntityHandle::MulticenterBond(MulticenterBondHandle::Id(host_multicenter(
                        left,
                    )?)),
                );
            }
            for (created, handle) in new_multicenter_handles {
                handles.insert(
                    Entity::MulticenterBond(created),
                    EntityHandle::MulticenterBond(handle),
                );
            }
            for left in self.lhs.noncovalent_bonds().ids() {
                handles.insert(
                    Entity::NoncovalentBond(left),
                    EntityHandle::NoncovalentBond(NoncovalentBondHandle::Id(host_noncovalent(
                        left,
                    )?)),
                );
            }
            for (created, handle) in new_noncovalent_handles {
                handles.insert(
                    Entity::NoncovalentBond(created),
                    EntityHandle::NoncovalentBond(handle),
                );
            }
            for left in self.lhs.stereo_atoms().ids() {
                handles.insert(
                    Entity::StereoAtom(left),
                    EntityHandle::StereoAtom(StereoAtomHandle::Id(host_stereo_atom(left)?)),
                );
            }
            for (created, handle) in new_stereo_atom_handles {
                handles.insert(
                    Entity::StereoAtom(created),
                    EntityHandle::StereoAtom(handle),
                );
            }
            for left in self.lhs.stereo_bonds().ids() {
                handles.insert(
                    Entity::StereoBond(left),
                    EntityHandle::StereoBond(StereoBondHandle::Id(host_stereo_bond(left)?)),
                );
            }
            for (created, handle) in new_stereo_bond_handles {
                handles.insert(
                    Entity::StereoBond(created),
                    EntityHandle::StereoBond(handle),
                );
            }
            for delta in constraint_deltas {
                let (constraint, add) = match delta {
                    ConstraintDelta::Add(constraint) => (constraint, true),
                    ConstraintDelta::Remove(constraint) => (constraint, false),
                };
                let constraint =
                    ConstraintEdit::new(constraint, |entity| handles.get(&entity).cloned())
                        .map_err(|_| ApplyError::InternalInvariant)?;
                if add {
                    constraint_edits.add_molecule_constraint(constraint);
                } else {
                    constraint_edits.remove_molecule_constraint(constraint);
                }
            }
        }

        let mut edits = Edits::new();
        if !created_atoms.is_empty() {
            edits.add_atoms(created_atoms.values().cloned());
        }
        if !created_bonds.is_empty() {
            edits.add_bonds(
                created_bonds
                    .values()
                    .map(|(atoms, ast)| {
                        Ok(AddBond {
                            endpoints: [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                            ast: ast.clone(),
                        })
                    })
                    .collect::<Result<Vec<_>, ApplyError>>()?,
            );
        }
        for edit in overlay_adds {
            edits.push(edit);
        }
        for edit in sets {
            edits.push(edit);
        }
        // Constraints precede all removals (overlay and topology) so each removal's constraint
        // compaction updates them — a constraint referencing a surviving overlay whose lower-id
        // sibling is removed would otherwise carry a stale id.
        for edit in constraint_edits {
            edits.push(edit);
        }
        for edit in overlay_removes {
            edits.push(edit);
        }
        if !remove_atoms.is_empty() || !remove_bonds.is_empty() {
            edits.remove_topology(remove_atoms, remove_bonds);
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

        // The host↔product comap: preserved host atoms match their compacted product id (survivors
        // keep ascending order), removed atoms are left-unmatched, created atoms right-unmatched.
        // `induce` derives the bond and overlay correspondences from this atom map.
        let removed: HashSet<AtomId> = removed_host_atoms.iter().copied().collect();
        let mut atom_matched_pairs: Vec<(NodeId, NodeId)> = Vec::new();
        let mut product_atom = 0u32;
        for host_atom in 0..host.atoms().count() as u32 {
            if removed.contains(&AtomId(host_atom)) {
                continue;
            }
            atom_matched_pairs.push((NodeId(host_atom), NodeId(product_atom)));
            product_atom += 1;
        }
        let atom_map = Correspondence::new(
            atom_matched_pairs,
            host.atoms().count(),
            product.atoms().count(),
        )
        .expect("correspondence producer preserves partial-bijection invariants");
        let comap = MoleculeCorrespondence::induce(host, &product, atom_map);
        Ok(ReactionDerivation::new(host.clone(), product, comap))
    }

    /// Validate the structural preconditions shared by every match against `host`.
    pub fn validate_application(&self, host: &MoleculeAst) -> Result<(), ApplyPreconditionError> {
        self.application_deltas(host).map(drop)
    }

    fn application_deltas(&self, host: &MoleculeAst) -> Result<Deltas, ApplyPreconditionError> {
        if !stereo_delta_domains_are_valid(&self.lhs, &self.deltas) {
            return Err(ApplyPreconditionError::InconsistentReaction);
        }
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
    /// `host` (using `match_config`) that satisfies the
    /// match-local DPO and structural conditions.
    /// Structural preconditions are checked before match enumeration. Match-local rejection is
    /// skipped; an internal application failure is yielded once and terminates the iterator.
    pub fn apply<'h>(
        &'h self,
        host: &'h MoleculeAst,
        match_config: SubstructureMatchConfig,
    ) -> Result<
        impl Iterator<Item = Result<ReactionDerivation, ApplyError>> + 'h,
        ApplyPreconditionError,
    > {
        let deltas = self.application_deltas(host)?;
        let mut correspondences = self
            .lhs
            .substructure_matches(host, match_config)
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
) -> Result<(), ApplyError> {
    let into_host = |l: &StereoLigand| {
        correspondence
            .atoms()
            .right_of(NodeId::from(l.atom_id))
            .map(|atom| StereoLigand::new(AtomId::from(atom), l.kind))
            .ok_or(ApplyError::InternalInvariant)
    };
    let from_host = |l: &StereoLigand| {
        correspondence
            .atoms()
            .left_of(NodeId::from(l.atom_id))
            .map(|atom| StereoLigand::new(AtomId::from(atom), l.kind))
            .ok_or(ApplyError::InternalInvariant)
    };
    for delta in deltas.iter_mut() {
        match delta {
            Delta::StereoAtom(s) => {
                let Some(host_id) = correspondence.stereo_atoms().right_of(s.id()) else {
                    continue;
                };
                let entity = Entity::StereoAtom(s.id());
                let rule_view = lhs
                    .stereo_atoms()
                    .get(s.id())
                    .ok_or(ApplyError::InternalInvariant)?;
                let host_view = host
                    .stereo_atoms()
                    .get(host_id)
                    .ok_or(ApplyError::InternalInvariant)?;
                let before: Vec<StereoLigand> = rule_view
                    .ligand_frame()
                    .iter()
                    .map(into_host)
                    .collect::<Result<_, _>>()?;
                let after = host_view.ligand_frame();
                match s {
                    StereoAtomDelta::ModifyField {
                        change: StereoAtomFieldChange::Configuration { old, new },
                        ..
                    } => {
                        let sigma = Permutation::between(&before, &after)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        *old = old.apply(sigma);
                        *new = new.apply(sigma);
                    }
                    StereoAtomDelta::Remove { ligands, ast, .. } => {
                        *ast = ast
                            .transform_frame(&before, &after)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        *ligands = after.iter().map(from_host).collect::<Result<_, _>>()?;
                    }
                    _ => {}
                }
            }
            Delta::StereoBond(s) => {
                let Some(host_id) = correspondence.stereo_bonds().right_of(s.id()) else {
                    continue;
                };
                let entity = Entity::StereoBond(s.id());
                let rule_view = lhs
                    .stereo_bonds()
                    .get(s.id())
                    .ok_or(ApplyError::InternalInvariant)?;
                let host_view = host
                    .stereo_bonds()
                    .get(host_id)
                    .ok_or(ApplyError::InternalInvariant)?;
                let before: Vec<StereoLigand> = rule_view
                    .ligand_frame()
                    .iter()
                    .map(into_host)
                    .collect::<Result<_, _>>()?;
                let after = host_view.ligand_frame();
                match s {
                    StereoBondDelta::ModifyField {
                        change: StereoBondFieldChange::Configuration { old, new },
                        ..
                    } => {
                        let sigma = Permutation::between(&before, &after)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        *old = old.apply(sigma);
                        *new = new.apply(sigma);
                    }
                    StereoBondDelta::Remove { ligands, ast, .. } => {
                        *ast = ast
                            .transform_frame(&before, &after)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        *ligands = after.iter().map(from_host).collect::<Result<_, _>>()?;
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }
    Ok(())
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
    use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};

    use super::super::constraint::{
        AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, Constraint, Constraints,
        DativeBondConstraintAst, MoleculeConstraint, MulticenterBondConstraintAst,
        NoncovalentBondConstraintAst, RelationalConstraint, StereoAtomConstraintAst,
        StereoBondConstraintAst, StereogenicityAst,
    };
    use super::super::edit::{AtomFieldChange, BondFieldChange};
    use super::super::entity::Entity;
    use super::super::ligand::StereoLigandKind;
    use super::super::molecule::transact::TransactionError;
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use super::super::stereo::{StereoAtomAst, StereoBondAst, StereoCoset, StereoKind};
    use super::super::substructure::SubstructureMatchAlgorithm;
    use super::super::validate::{DpoContradiction, EntityStructureContradiction};
    use super::super::value::ValueAst;
    use super::*;

    const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
        match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    };

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
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2)
            .expect("correspondence producer preserves partial-bijection invariants");
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

    #[rustfmt::skip]
    #[rstest]
    #[case::missing_atom(
        ReactionAst::new(MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }), Deltas::new()),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C)], ..Default::default() }),
        MoleculeCorrespondence::new(
            Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        ),
        ApplyError::CorrespondenceMismatch { entity: Entity::Atom(AtomId(0)) },
    )]
    #[case::bond_incidence(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
            Deltas::new(),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))], ..Default::default() }),
        MoleculeCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1)], 2),
            Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        ),
        ApplyError::CorrespondenceMismatch { entity: Entity::Bond(BondId(0)) },
    )]
    #[case::noncovalent_incidence(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O); 3], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::new(),
        ),
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst::from_element(Element::O); 3], noncovalent: vec![(AtomId(0), AtomId(2), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        MoleculeCorrespondence::new(
            Correspondence::from_images(&[NodeId(0), NodeId(1), NodeId(2)], 3),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        ),
        ApplyError::CorrespondenceMismatch { entity: Entity::NoncovalentBond(NoncovalentBondId(0)) },
    )]
    fn test_reaction_ast_apply_at_correspondence_error(
        #[case] reaction: ReactionAst,
        #[case] host: MoleculeAst,
        #[case] correspondence: MoleculeCorrespondence,
        #[case] expected: ApplyError,
    ) {
        assert_eq!(reaction.apply_at(&host, &correspondence), Err(expected));
    }

    #[rstest]
    #[case::field(Delta::StereoAtom(StereoAtomDelta::ModifyField {
        id: StereoAtomId(0),
        change: StereoAtomFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0u32),
            new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 1u32),
        },
    }))]
    #[case::removal(Delta::StereoAtom(StereoAtomDelta::Remove {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ],
        ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
    }))]
    fn test_reaction_ast_apply_at_stereo_atom_error(#[case] delta: Delta) {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 6],
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
            ..Default::default()
        });
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 6],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new((0..6u32).map(|id| (NodeId(id), NodeId(id))).collect(), 6, 6)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let reaction = ReactionAst::new(lhs, Deltas::from_iter([delta]));

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap_err(),
            ApplyError::StereoFrameMismatch {
                entity: Entity::StereoAtom(StereoAtomId(0)),
            },
        );
    }

    #[rstest]
    #[case::field(Delta::StereoBond(StereoBondDelta::ModifyField {
        id: StereoBondId(0),
        change: StereoBondFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0u32),
            new: StereoConfigurationAst::kinded(StereoKind::CisTrans, 1u32),
        },
    }))]
    #[case::removal(Delta::StereoBond(StereoBondDelta::Remove {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: vec![
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
        ],
        ast: StereoBondAst::new(StereoKind::CisTrans, 0u32),
    }))]
    fn test_reaction_ast_apply_at_stereo_bond_error(#[case] delta: Delta) {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 7],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 7],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new((0..7u32).map(|id| (NodeId(id), NodeId(id))).collect(), 7, 7)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(BondId(0), BondId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(StereoBondId(0), StereoBondId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let reaction = ReactionAst::new(lhs, Deltas::from_iter([delta]));

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap_err(),
            ApplyError::StereoFrameMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            },
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
    fn test_reaction_ast_apply_at_molecule_constraint_created() {
        let constraint = Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintAst::valence(3_i64)),
            Constraint::Bond(BondId(0), BondConstraintAst::aromatic(true)),
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintAst::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintAst::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(0),
                MulticenterBondConstraintAst::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(0),
                NoncovalentBondConstraintAst::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0),
                parallel: BondId(0),
            }),
        ]);
        let reaction = ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(0),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::DativeBond(DativeBondDelta::Add {
                    id: DativeBondId(0),
                    donors: vec![AtomId(0)],
                    acceptor: AtomId(1),
                    ast: DativeBondAst::from_order(1),
                }),
                Delta::AromaticSystem(AromaticSystemDelta::Add {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    ast: AromaticSystemAst::default(),
                }),
                Delta::MulticenterBond(MulticenterBondDelta::Add {
                    id: MulticenterBondId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    ast: MulticenterBondAst::default(),
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
                Delta::StereoAtom(StereoAtomDelta::Add {
                    id: StereoAtomId(0),
                    site: AtomId(0),
                    ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                    ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                }),
                Delta::StereoBond(StereoBondDelta::Add {
                    id: StereoBondId(0),
                    site: BondId(0),
                    ligands: vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    ],
                    ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                }),
                Delta::Constraint(ConstraintDelta::Add(constraint.clone())),
            ]),
        );

        let host = MoleculeAst::default();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let result = reaction.apply_at(&host, &correspondence).unwrap();

        assert_eq!(result.rhs().constraints(), &Constraints::from(constraint));
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
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        assert_eq!(products, expected);
    }

    #[rstest]
    #[case::graph_and_overlays(SubstructureMatchAlgorithm::GraphAndOverlays)]
    #[case::incidence(SubstructureMatchAlgorithm::Incidence)]
    fn test_reaction_ast_apply_match_algorithm(
        #[case] match_algorithm: SubstructureMatchAlgorithm,
    ) {
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = ReactionAst::new(host.clone(), Deltas::new());
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                SubstructureMatchConfig {
                    match_algorithm,
                    ..MATCH_CONFIG
                },
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        assert_eq!(products, vec![host]);
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
        let mut applications = reaction.apply(&host, MATCH_CONFIG).unwrap();

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
        match reaction.apply(&host, MATCH_CONFIG) {
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
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
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
            .apply(&host, MATCH_CONFIG)
            .unwrap()
            .next()
            .expect("the inversion rule matches the host")
            .unwrap()
            .rhs()
            .clone();
        assert_eq!(rhs, expected);
    }

    // A stereo-bond addition can name a bond added by the same reaction. Its `BondHandle::New`
    // index is the bond creation ordinal and is independent of the atom creation namespace.
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
            .apply(&host, MATCH_CONFIG)
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
    #[case::undetermined(StereoCoset::Undetermined)]
    fn test_reaction_ast_apply_two_stereo_centers(#[case] coset: StereoCoset) {
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
            .apply(&center, MATCH_CONFIG)
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
        // Remove atom O (id 1) and its bond: host C-O ⇒ product C. Atom 0 is preserved (matched),
        // atom 1 is deleted (left-unmatched), so the comap's atom map records exactly that.
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
        assert_eq!(
            derivation.atom_map().matched_pairs(),
            &[(NodeId(0), NodeId(0))]
        );
        assert_eq!(derivation.atom_map().left_unmatched(), vec![NodeId(1)]);
    }
}
