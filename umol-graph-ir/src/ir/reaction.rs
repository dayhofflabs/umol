//! Reaction graph IR: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations live in `reaction_span.rs`).

mod dpo;
mod integrity;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::vec::IntoIter;

use dpo::check_reaction_dpo;
pub use dpo::DpoContradiction;
pub use integrity::ReactionIntegrityError;
use umol_graph_core::Correspondence;
use umol_perm::{DynPermutation, Permutation};

use super::aromatic::{
    aromatic_system_representative_action, AromaticSystemForm, AromaticSystemUpdate,
};
use super::atom::{AtomForm, AtomUpdate};
use super::bond::{BondForm, BondUpdate};
use super::constraint::{
    ConstraintFrameActionDomain, ConstraintFrameActionMap, ConstraintFrameActions,
};
use super::correspondence::MoleculeCorrespondence;
use super::dative::{dative_bond_representative_action, DativeBondForm, DativeBondUpdate};
use super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::derivation::ReactionDerivation;
use super::edit::{
    AddBond, AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle,
    BondFieldChange, BondHandle, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit,
    Edits, EntityHandle, MulticenterBondFieldChange, MulticenterBondHandle,
    NoncovalentBondFieldChange, NoncovalentBondHandle, StereoAtomFieldChange, StereoAtomHandle,
    StereoAtomRemoval, StereoBondFieldChange, StereoBondHandle, StereoBondRemoval,
};
use super::entity::Entity;
use super::error::{ApplyError, ApplyPreconditionError, Contradiction};
use super::frame::{
    AromaticSystemsFrameAction, DativeBondsFrameAction, MulticenterBondsFrameAction,
    NoncovalentBondsFrameAction, OverlaysFrameAction, StereoAtomsFrameAction,
    StereoBondsFrameAction,
};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
#[cfg(test)]
use super::molecule::MoleculeEntries;
use super::molecule::{Molecule, MoleculeIntegrityError};
use super::multicenter::{
    multicenter_bond_representative_action, MulticenterBondForm, MulticenterBondUpdate,
};
use super::noncovalent::{
    noncovalent_bond_representative_action, NoncovalentBondForm, NoncovalentBondUpdate,
};
use super::stereo::{
    stereo_atom_representative_action, stereo_bond_representative_action, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoTerm,
};
use super::substructure::SubstructureMatchConfig;
use super::traits::{FrameTransport, Lattice, Normalize, Reframe};

/// A reaction as one full molecule state (`lhs`) plus one resolved delta (`deltas`).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct Reaction {
    lhs: Molecule,
    deltas: Deltas,
}

/// One-shot reaction applications over an eagerly enumerated correspondence set.
///
/// This operation-issued iterator is created by [`Reaction::apply`] and has no independent public
/// constructor. It owns snapshots of the reaction and host plus the reaction's normalized deltas,
/// so later changes to the inputs cannot affect iteration. Matching is completed when the iterator
/// is created; derivations are constructed lazily in match order. Match-local rejection is skipped;
/// another application failure is yielded once and terminates the iterator.
#[derive(Debug)]
pub struct ReactionApplicationIter {
    reaction: Reaction,
    host: Molecule,
    deltas: Deltas,
    correspondences: IntoIter<MoleculeCorrespondence>,
    failed: bool,
}

impl ReactionApplicationIter {
    fn new(
        reaction: Reaction,
        host: Molecule,
        match_config: SubstructureMatchConfig,
    ) -> Result<Self, ApplyPreconditionError> {
        let deltas = reaction.application_deltas()?;
        let correspondences = reaction
            .lhs
            .substructure_matches(&host, match_config)?
            .into_iter();
        Ok(Self {
            reaction,
            host,
            deltas,
            correspondences,
            failed: false,
        })
    }
}

impl Iterator for ReactionApplicationIter {
    type Item = Result<ReactionDerivation, ApplyError>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.failed {
            let correspondence = self.correspondences.next()?;
            match self
                .reaction
                .apply_at_canonical(&self.host, &correspondence, self.deltas.clone())
            {
                Ok(derivation) => return Some(Ok(derivation)),
                Err(error) if error.is_match_rejection() => {}
                Err(error) => {
                    self.failed = true;
                    return Some(Err(error));
                }
            }
        }
        None
    }
}

/// One-shot product component collections derived lazily from reaction applications.
///
/// This operation-issued iterator is created through [`React`] and has no independent public
/// constructor. It owns the underlying [`ReactionApplicationIter`] and therefore the reaction and
/// reactant snapshots. Each successful application is replaced lazily by the conservative
/// connected-component split of its right-hand side. Component order is inherited from
/// [`Molecule::split`]. The split correspondences and the rest of the derivation are intentionally
/// discarded. Application errors pass through unchanged.
#[derive(Debug)]
pub struct ReactionProductsIter {
    applications: ReactionApplicationIter,
}

impl ReactionProductsIter {
    fn new(applications: ReactionApplicationIter) -> Self {
        Self { applications }
    }
}

impl Iterator for ReactionProductsIter {
    type Item = Result<Vec<Molecule>, ApplyError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.applications.next().map(|application| {
            application.map(|derivation| {
                derivation
                    .rhs()
                    .split()
                    .into_iter()
                    .map(|(component, _)| component)
                    .collect()
            })
        })
    }
}

/// Apply a reaction and emit only the connected product components of each successful match.
///
/// Implementations preserve reaction match order, [`Molecule::split`] component order, and both
/// application error channels. The operation borrows reusable reactants and takes an explicit
/// matching configuration; its returned iterator owns snapshots of every input.
///
/// # Semantic properties
///
/// For a molecule slice `reactants`, the result is identical to combining `reactants` in slice
/// order, applying `reaction`, and splitting every successful derivation's right-hand side while
/// discarding the split correspondences. An empty slice follows the same rule through the empty
/// combined molecule.
///
/// # Examples
///
/// ```
/// use umol_graph_core::{
///     RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
/// };
/// use umol_graph_ir::ir::{
///     Molecule, React, Reaction, SubstructureMatchAlgorithm, SubstructureMatchConfig,
/// };
///
/// let reactants = [Molecule::new(), Molecule::new()];
/// let reaction = Reaction::default();
/// let config = SubstructureMatchConfig {
///     match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
///     subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
///     relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
/// };
/// let product_sets = reactants.react(&reaction, config)?;
///
/// for products in product_sets {
///     for product in products? {
///         println!("{product:?}");
///     }
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait React {
    /// Apply `reaction` and lazily emit the connected product components for every successful
    /// match.
    ///
    /// The returned iterator owns snapshots of the reaction and reactants. Matching is eager;
    /// product construction and splitting are lazy.
    ///
    /// # Errors
    ///
    /// Returns [`ApplyPreconditionError`] before issuing an iterator when the reaction fails a
    /// reaction-wide application precondition. Once issued, the iterator yields [`ApplyError`] for
    /// a non-rejection failure while realizing a selected match.
    fn react(
        &self,
        reaction: &Reaction,
        match_config: SubstructureMatchConfig,
    ) -> Result<ReactionProductsIter, ApplyPreconditionError>;
}

impl React for Molecule {
    fn react(
        &self,
        reaction: &Reaction,
        match_config: SubstructureMatchConfig,
    ) -> Result<ReactionProductsIter, ApplyPreconditionError> {
        reaction
            .apply(self, match_config)
            .map(ReactionProductsIter::new)
    }
}

impl React for [Molecule] {
    fn react(
        &self,
        reaction: &Reaction,
        match_config: SubstructureMatchConfig,
    ) -> Result<ReactionProductsIter, ApplyPreconditionError> {
        let (host, _) = Molecule::combine_all(self);
        ReactionApplicationIter::new(reaction.clone(), host, match_config)
            .map(ReactionProductsIter::new)
    }
}

fn stereo_delta_domains_are_valid(lhs: &Molecule, deltas: &Deltas) -> bool {
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

    fn configuration_is_valid(configuration: &StereoConfigurationForm) -> bool {
        match configuration {
            StereoConfigurationForm::Undetermined => true,
            StereoConfigurationForm::Kinded(kind, coset) => match coset {
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
        lhs: Option<&StereoConfigurationForm>,
        old: &StereoConfigurationForm,
        new: &StereoConfigurationForm,
    ) -> bool {
        let lhs_kind = lhs.and_then(StereoConfigurationForm::kind);
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
        .any(|view| !configuration_is_valid(&view.attributes.configuration))
        || lhs
            .stereo_bonds()
            .iter()
            .any(|view| !configuration_is_valid(&view.attributes.configuration))
    {
        return false;
    }

    deltas.iter().all(|delta| match delta {
        Delta::StereoAtom(StereoAtomDelta::Add { attributes, .. }) => {
            configuration_is_valid(&attributes.configuration)
        }
        Delta::StereoAtom(StereoAtomDelta::Remove { id, attributes, .. }) => {
            let lhs_configuration = lhs
                .stereo_atoms()
                .get(*id)
                .map(|view| &view.attributes.configuration);
            configurations_are_compatible(
                lhs_configuration,
                &attributes.configuration,
                &attributes.configuration,
            )
        }
        Delta::StereoAtom(StereoAtomDelta::ModifyField {
            id,
            change: StereoAtomFieldChange::Configuration { old, new },
        }) => configurations_are_compatible(
            lhs.stereo_atoms()
                .get(*id)
                .map(|view| &view.attributes.configuration),
            old,
            new,
        ),
        Delta::StereoBond(StereoBondDelta::Add { attributes, .. }) => {
            configuration_is_valid(&attributes.configuration)
        }
        Delta::StereoBond(StereoBondDelta::Remove { id, attributes, .. }) => {
            let lhs_configuration = lhs
                .stereo_bonds()
                .get(*id)
                .map(|view| &view.attributes.configuration);
            configurations_are_compatible(
                lhs_configuration,
                &attributes.configuration,
                &attributes.configuration,
            )
        }
        Delta::StereoBond(StereoBondDelta::ModifyField {
            id,
            change: StereoBondFieldChange::Configuration { old, new },
        }) => configurations_are_compatible(
            lhs.stereo_bonds()
                .get(*id)
                .map(|view| &view.attributes.configuration),
            old,
            new,
        ),
        _ => true,
    })
}

/// Normalize reaction deltas after restating every overlay removal from its explicit local frame
/// into the entity's owning frame. Existing entities are owned by the lhs; created entities are
/// owned by their unique add delta. This contextual step must precede per-entity folding because
/// field and constraint deltas are stated in the owning frame.
pub(super) fn normalize_reaction_deltas(
    lhs: &Molecule,
    deltas: &Deltas,
) -> Result<Deltas, Contradiction> {
    let mut dative_adds = HashMap::new();
    let mut aromatic_adds = HashMap::new();
    let mut multicenter_adds = HashMap::new();
    let mut noncovalent_adds = HashMap::new();
    let mut stereo_atom_adds = HashMap::new();
    let mut stereo_bond_adds = HashMap::new();
    for delta in deltas.iter() {
        match delta {
            Delta::DativeBond(DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ..
            }) => {
                dative_adds.insert(*id, (donors.clone(), *acceptor));
            }
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. }) => {
                aromatic_adds.insert(*id, atoms.clone());
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. }) => {
                multicenter_adds.insert(*id, atoms.clone());
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, .. }) => {
                noncovalent_adds.insert(*id, *atoms);
            }
            Delta::StereoAtom(StereoAtomDelta::Add {
                id, site, ligands, ..
            }) => {
                stereo_atom_adds.insert(*id, (*site, ligands.clone()));
            }
            Delta::StereoBond(StereoBondDelta::Add {
                id, site, ligands, ..
            }) => {
                stereo_bond_adds.insert(*id, (*site, ligands.clone()));
            }
            _ => {}
        }
    }

    let mut restated = deltas.clone();
    for delta in restated.iter_mut() {
        match delta {
            Delta::DativeBond(DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                attributes,
            }) => {
                let owner = lhs
                    .dative_bonds()
                    .get(*id)
                    .map(|view| (view.donor_ids().collect(), view.acceptor_id()))
                    .or_else(|| dative_adds.get(id).cloned())
                    .ok_or(Contradiction)?;
                let action = DynPermutation::between(donors, &owner.0).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *donors = owner.0;
                *acceptor = owner.1;
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove {
                id,
                atoms,
                attributes,
            }) => {
                let owner = lhs
                    .aromatic_systems()
                    .get(*id)
                    .map(|view| view.atom_ids().collect())
                    .or_else(|| aromatic_adds.get(id).cloned())
                    .ok_or(Contradiction)?;
                let action = DynPermutation::between(atoms, &owner).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *atoms = owner;
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove {
                id,
                atoms,
                attributes,
            }) => {
                let owner = lhs
                    .multicenter_bonds()
                    .get(*id)
                    .map(|view| view.atom_ids().collect())
                    .or_else(|| multicenter_adds.get(id).cloned())
                    .ok_or(Contradiction)?;
                let action = DynPermutation::between(atoms, &owner).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *atoms = owner;
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                id,
                atoms,
                attributes,
            }) => {
                let owner = lhs
                    .noncovalent_bonds()
                    .get(*id)
                    .map(|view| view.atom_ids())
                    .or_else(|| noncovalent_adds.get(id).copied())
                    .ok_or(Contradiction)?;
                let action = DynPermutation::between(atoms, &owner).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *atoms = owner;
            }
            Delta::StereoAtom(StereoAtomDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            }) => {
                let owner = lhs
                    .stereo_atoms()
                    .get(*id)
                    .map(|view| (view.site_id(), view.ligand_frame()))
                    .or_else(|| stereo_atom_adds.get(id).cloned())
                    .ok_or(Contradiction)?;
                let action = Permutation::between(ligands, &owner.1).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *site = owner.0;
                *ligands = owner.1;
            }
            Delta::StereoBond(StereoBondDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            }) => {
                let owner = lhs
                    .stereo_bonds()
                    .get(*id)
                    .map(|view| (view.site_id(), view.ligand_frame()))
                    .or_else(|| stereo_bond_adds.get(id).cloned())
                    .ok_or(Contradiction)?;
                let action = Permutation::between(ligands, &owner.1).ok_or(Contradiction)?;
                *attributes = attributes
                    .clone()
                    .reframe_by(&action)
                    .ok_or(Contradiction)?;
                *site = owner.0;
                *ligands = owner.1;
            }
            _ => {}
        }
    }
    restated.normalize()
}

#[derive(Default)]
struct ReactionFrameActionDomain {
    deltas: ConstraintFrameActionDomain,
    constraints: ConstraintFrameActionDomain,
}

impl ReactionFrameActionDomain {
    fn from_deltas(deltas: &Deltas) -> Self {
        let mut domain = Self::default();
        for delta in deltas.iter() {
            match delta {
                Delta::Constraint(delta) => {
                    delta.collect_frame_action_domain(&mut domain.constraints)
                }
                delta => delta.collect_frame_action_domain(&mut domain.deltas),
            }
        }
        domain
    }

    fn contains_dative_bond(&self, id: DativeBondId) -> bool {
        self.deltas.contains_dative_bond(id) || self.constraints.contains_dative_bond(id)
    }

    fn contains_aromatic_system(&self, id: AromaticSystemId) -> bool {
        self.deltas.contains_aromatic_system(id) || self.constraints.contains_aromatic_system(id)
    }

    fn contains_multicenter_bond(&self, id: MulticenterBondId) -> bool {
        self.deltas.contains_multicenter_bond(id) || self.constraints.contains_multicenter_bond(id)
    }

    fn contains_noncovalent_bond(&self, id: NoncovalentBondId) -> bool {
        self.deltas.contains_noncovalent_bond(id) || self.constraints.contains_noncovalent_bond(id)
    }

    fn contains_stereo_atom(&self, id: StereoAtomId) -> bool {
        self.deltas.contains_stereo_atom(id) || self.constraints.contains_stereo_atom(id)
    }

    fn contains_stereo_bond(&self, id: StereoBondId) -> bool {
        self.deltas.contains_stereo_bond(id) || self.constraints.contains_stereo_bond(id)
    }
}

fn reaction_frame_action(
    lhs: &Molecule,
    deltas: &Deltas,
    domain: Option<&ReactionFrameActionDomain>,
) -> OverlaysFrameAction {
    let mut dative_bonds = BTreeMap::new();
    let mut aromatic_systems = BTreeMap::new();
    let mut multicenter_bonds = BTreeMap::new();
    let mut noncovalent_bonds = BTreeMap::new();
    let mut stereo_atoms = BTreeMap::new();
    let mut stereo_bonds = BTreeMap::new();

    for view in lhs.dative_bonds().iter() {
        if domain.is_none_or(|domain| domain.contains_dative_bond(view.id)) {
            dative_bonds.insert(
                view.id,
                dative_bond_representative_action(view.donor_ids().collect()),
            );
        }
    }
    for view in lhs.aromatic_systems().iter() {
        if domain.is_none_or(|domain| domain.contains_aromatic_system(view.id)) {
            aromatic_systems.insert(
                view.id,
                aromatic_system_representative_action(view.atom_ids().collect()),
            );
        }
    }
    for view in lhs.multicenter_bonds().iter() {
        if domain.is_none_or(|domain| domain.contains_multicenter_bond(view.id)) {
            multicenter_bonds.insert(
                view.id,
                multicenter_bond_representative_action(view.atom_ids().collect()),
            );
        }
    }
    for view in lhs.noncovalent_bonds().iter() {
        if domain.is_none_or(|domain| domain.contains_noncovalent_bond(view.id)) {
            noncovalent_bonds.insert(
                view.id,
                noncovalent_bond_representative_action(view.atom_ids()),
            );
        }
    }
    for view in lhs.stereo_atoms().iter() {
        if domain.is_none_or(|domain| domain.contains_stereo_atom(view.id)) {
            stereo_atoms.insert(
                view.id,
                stereo_atom_representative_action(&view.ligand_frame())
                    .expect("integrity-valid stereo-atom frames fit the bounded action"),
            );
        }
    }
    for view in lhs.stereo_bonds().iter() {
        if domain.is_none_or(|domain| domain.contains_stereo_bond(view.id)) {
            stereo_bonds.insert(
                view.id,
                stereo_bond_representative_action(&view.ligand_frame())
                    .expect("integrity-valid stereo-bond frames admit a standard-frame action"),
            );
        }
    }

    for delta in deltas.iter() {
        match delta {
            Delta::DativeBond(DativeBondDelta::Add { id, donors, .. })
                if domain.is_none_or(|domain| domain.contains_dative_bond(*id)) =>
            {
                dative_bonds.insert(*id, dative_bond_representative_action(donors.clone()));
            }
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. })
                if domain.is_none_or(|domain| domain.contains_aromatic_system(*id)) =>
            {
                aromatic_systems.insert(*id, aromatic_system_representative_action(atoms.clone()));
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. })
                if domain.is_none_or(|domain| domain.contains_multicenter_bond(*id)) =>
            {
                multicenter_bonds
                    .insert(*id, multicenter_bond_representative_action(atoms.clone()));
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, .. })
                if domain.is_none_or(|domain| domain.contains_noncovalent_bond(*id)) =>
            {
                noncovalent_bonds.insert(*id, noncovalent_bond_representative_action(*atoms));
            }
            Delta::StereoAtom(StereoAtomDelta::Add { id, ligands, .. })
                if domain.is_none_or(|domain| domain.contains_stereo_atom(*id)) =>
            {
                stereo_atoms.insert(
                    *id,
                    stereo_atom_representative_action(ligands)
                        .expect("integrity-valid stereo-atom frames fit the bounded action"),
                );
            }
            Delta::StereoBond(StereoBondDelta::Add { id, ligands, .. })
                if domain.is_none_or(|domain| domain.contains_stereo_bond(*id)) =>
            {
                stereo_bonds.insert(
                    *id,
                    stereo_bond_representative_action(ligands)
                        .expect("integrity-valid stereo-bond frames admit a standard-frame action"),
                );
            }
            _ => {}
        }
    }

    OverlaysFrameAction::new(
        DativeBondsFrameAction::from_action_map(dative_bonds)
            .expect("every dynamic permutation is a dative-bond action"),
        AromaticSystemsFrameAction::from_action_map(aromatic_systems)
            .expect("every dynamic permutation is an aromatic-system action"),
        MulticenterBondsFrameAction::from_action_map(multicenter_bonds)
            .expect("every dynamic permutation is a multicenter-bond action"),
        NoncovalentBondsFrameAction::from_action_map(noncovalent_bonds)
            .expect("every noncovalent-bond action has degree two"),
        StereoAtomsFrameAction::from_action_map(stereo_atoms)
            .expect("every bounded permutation is a stereo-atom action"),
        StereoBondsFrameAction::from_action_map(stereo_bonds)
            .expect("every selected stereo-bond action preserves endpoint blocks"),
    )
}

fn reframe_reaction_deltas(
    lhs: &Molecule,
    deltas: Deltas,
    actions: &OverlaysFrameAction,
) -> Option<Deltas> {
    let mut dative_adds = HashMap::new();
    let mut aromatic_adds = HashMap::new();
    let mut multicenter_adds = HashMap::new();
    let mut noncovalent_adds = HashMap::new();
    let mut stereo_atom_adds = HashMap::new();
    let mut stereo_bond_adds = HashMap::new();
    for delta in deltas.iter() {
        match delta {
            Delta::DativeBond(DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ..
            }) => {
                dative_adds.insert(*id, (donors.clone(), *acceptor));
            }
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. }) => {
                aromatic_adds.insert(*id, atoms.clone());
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. }) => {
                multicenter_adds.insert(*id, atoms.clone());
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, .. }) => {
                noncovalent_adds.insert(*id, *atoms);
            }
            Delta::StereoAtom(StereoAtomDelta::Add {
                id, site, ligands, ..
            }) => {
                stereo_atom_adds.insert(*id, (*site, ligands.clone()));
            }
            Delta::StereoBond(StereoBondDelta::Add {
                id, site, ligands, ..
            }) => {
                stereo_bond_adds.insert(*id, (*site, ligands.clone()));
            }
            _ => {}
        }
    }

    // A raw removal stays in its local coordinates while its owner moves, so transport uses the
    // conjugate of the owner action by the local-to-owner action.
    deltas
        .into_iter()
        .map(|delta| {
            Some(match delta {
                Delta::DativeBond(delta @ DativeBondDelta::Remove { .. }) => {
                    let DativeBondDelta::Remove { id, donors, .. } = &delta else {
                        unreachable!()
                    };
                    let id = *id;
                    let owner = lhs
                        .dative_bonds()
                        .get(id)
                        .map(|view| (view.donor_ids().collect(), view.acceptor_id()))
                        .or_else(|| dative_adds.get(&id).cloned())?;
                    let local_to_owner = DynPermutation::between(donors, &owner.0)?;
                    let local_action = local_to_owner
                        .compose(actions.dative_bonds().action(id)?)?
                        .compose(&local_to_owner.inverse())?;
                    Delta::DativeBond(delta.reframe_by(&local_action)?)
                }
                Delta::AromaticSystem(delta @ AromaticSystemDelta::Remove { .. }) => {
                    let AromaticSystemDelta::Remove { id, atoms, .. } = &delta else {
                        unreachable!()
                    };
                    let id = *id;
                    let owner = lhs
                        .aromatic_systems()
                        .get(id)
                        .map(|view| view.atom_ids().collect())
                        .or_else(|| aromatic_adds.get(&id).cloned())?;
                    let local_to_owner = DynPermutation::between(atoms, &owner)?;
                    let local_action = local_to_owner
                        .compose(actions.aromatic_systems().action(id)?)?
                        .compose(&local_to_owner.inverse())?;
                    Delta::AromaticSystem(delta.reframe_by(&local_action)?)
                }
                Delta::MulticenterBond(delta @ MulticenterBondDelta::Remove { .. }) => {
                    let MulticenterBondDelta::Remove { id, atoms, .. } = &delta else {
                        unreachable!()
                    };
                    let id = *id;
                    let owner = lhs
                        .multicenter_bonds()
                        .get(id)
                        .map(|view| view.atom_ids().collect())
                        .or_else(|| multicenter_adds.get(&id).cloned())?;
                    let local_to_owner = DynPermutation::between(atoms, &owner)?;
                    let local_action = local_to_owner
                        .compose(actions.multicenter_bonds().action(id)?)?
                        .compose(&local_to_owner.inverse())?;
                    Delta::MulticenterBond(delta.reframe_by(&local_action)?)
                }
                Delta::NoncovalentBond(delta @ NoncovalentBondDelta::Remove { .. }) => {
                    let NoncovalentBondDelta::Remove { id, atoms, .. } = &delta else {
                        unreachable!()
                    };
                    let id = *id;
                    let owner = lhs
                        .noncovalent_bonds()
                        .get(id)
                        .map(|view| view.atom_ids())
                        .or_else(|| noncovalent_adds.get(&id).copied())?;
                    let local_to_owner = DynPermutation::between(atoms, &owner)?;
                    let local_action = local_to_owner
                        .compose(actions.noncovalent_bonds().action(id)?)?
                        .compose(&local_to_owner.inverse())?;
                    Delta::NoncovalentBond(delta.reframe_by(&local_action)?)
                }
                Delta::StereoAtom(delta @ StereoAtomDelta::Remove { .. }) => {
                    let id = delta.id();
                    let owner = lhs
                        .stereo_atoms()
                        .get(id)
                        .map(|view| (view.site_id(), view.ligand_frame()))
                        .or_else(|| stereo_atom_adds.get(&id).cloned())?;
                    let StereoAtomDelta::Remove { ligands, .. } = &delta else {
                        unreachable!()
                    };
                    let local_to_owner = Permutation::between(ligands, &owner.1)?;
                    let owner_to_target = *actions.stereo_atoms().action(id)?;
                    (local_to_owner.degree() == owner_to_target.degree()).then_some(())?;
                    let local_action = local_to_owner
                        .compose(owner_to_target)
                        .compose(local_to_owner.inverse());
                    Delta::StereoAtom(delta.reframe_by(&local_action)?)
                }
                Delta::StereoBond(delta @ StereoBondDelta::Remove { .. }) => {
                    let id = delta.id();
                    let owner = lhs
                        .stereo_bonds()
                        .get(id)
                        .map(|view| (view.site_id(), view.ligand_frame()))
                        .or_else(|| stereo_bond_adds.get(&id).cloned())?;
                    let StereoBondDelta::Remove { ligands, .. } = &delta else {
                        unreachable!()
                    };
                    let local_to_owner = Permutation::between(ligands, &owner.1)?;
                    let owner_to_target = *actions.stereo_bonds().action(id)?;
                    (local_to_owner.degree() == owner_to_target.degree()).then_some(())?;
                    let local_action = local_to_owner
                        .compose(owner_to_target)
                        .compose(local_to_owner.inverse());
                    Delta::StereoBond(delta.reframe_by(&local_action)?)
                }
                delta => delta.reframe_by(actions)?,
            })
        })
        .collect()
}

impl Normalize for Reaction {
    fn normalize(self) -> Result<Self, Contradiction> {
        let deltas = normalize_reaction_deltas(&self.lhs, &self.deltas)?;
        let lhs = self.lhs.normalize()?;
        Ok(Self { lhs, deltas })
    }
}

impl FrameTransport for Reaction {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        let deltas = reframe_reaction_deltas(&self.lhs, self.deltas, actions)?;
        let lhs = self.lhs.reframe_by(actions)?;
        Some(Self { lhs, deltas })
    }
}

impl Reframe for Reaction {
    fn representative_action(&self) -> Self::Action {
        reaction_frame_action(&self.lhs, &self.deltas, None)
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        let normalized = self.normalize()?;
        let domain = ReactionFrameActionDomain::from_deltas(&normalized.deltas);
        let actions = reaction_frame_action(&normalized.lhs, &normalized.deltas, Some(&domain));
        let deltas = reframe_reaction_deltas(&normalized.lhs, normalized.deltas, &actions)
            .ok_or(Contradiction)?;
        let lhs = normalized.lhs.reframe()?;
        Self { lhs, deltas }.normalize()
    }
}

impl Reaction {
    /// Construct a reaction from an lhs and deltas, asserting representation integrity.
    ///
    /// # Panics
    ///
    /// Panics if the deltas violate reaction representation integrity relative to the closed lhs.
    /// Use [`Self::try_new`] for independently assembled deltas.
    pub fn new(lhs: Molecule, deltas: Deltas) -> Self {
        Self::try_new(lhs, deltas).expect("invalid reaction")
    }

    /// Construct a reaction after checking its representation integrity.
    ///
    /// The check covers delta references, added stereo entries, stereo constraint site kinds, and
    /// removal incidence compatible with each source entity's participant structure. The lhs is an
    /// already closed [`Molecule`]. Removal payloads are interpreted in their recorded local frame.
    /// The check does not require the deltas to materialize a reaction span or impose DPO or
    /// chemistry semantics.
    ///
    /// # Errors
    ///
    /// Returns [`ReactionIntegrityError`] when the supplied parts do not form an interpretable
    /// reaction representation.
    pub fn try_new(lhs: Molecule, deltas: Deltas) -> Result<Self, ReactionIntegrityError> {
        let reaction = Self { lhs, deltas };
        reaction.check_integrity()?;
        Ok(reaction)
    }

    /// Borrow the left-hand-side molecule.
    pub fn lhs(&self) -> &Molecule {
        &self.lhs
    }

    /// Borrow the resolved delta collection.
    pub fn deltas(&self) -> &Deltas {
        &self.deltas
    }

    /// Consume the reaction and return its left-hand side and deltas.
    pub fn into_parts(self) -> (Molecule, Deltas) {
        (self.lhs, self.deltas)
    }

    /// The reaction transforming `lhs` into `rhs` under `atom_correspondence`: induce the full
    /// per-entity correspondence, diff the two sides into deltas, and pair them with `lhs`. The
    /// inverse of reading a reaction's two sides back off its span. Returns `None` when the atom
    /// correspondence is not compatible with the supplied sides or their entity incidence does not
    /// induce unique partners.
    ///
    /// # Semantic properties
    ///
    /// The returned reaction retains `lhs` exactly. Its materialized rhs is the supplied `rhs`
    /// reindexed into the lhs-anchored reaction frame: preserved entities retain lhs ids and
    /// rhs-only entities are appended. The two right-hand molecules are equivalent under the
    /// induced total correspondence, but need not be structurally equal when matched pairs cross
    /// entity order.
    pub fn from_sides(
        lhs: Molecule,
        rhs: Molecule,
        atom_correspondence: Correspondence<AtomId>,
    ) -> Option<Self> {
        let correspondence = MoleculeCorrespondence::induce(&lhs, &rhs, atom_correspondence)?;
        let deltas = lhs
            .difference_to(&rhs, &correspondence)
            .expect("induced molecule correspondence is compatible with its source molecules");
        Some(Self::new(lhs, deltas))
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
    ///
    /// # Semantic properties
    ///
    /// Every matched overlay delta is transported from its normalized rule-owner frame into the
    /// host's stored participant frame before lowering. A rule-side `old` value is a pattern: it
    /// must match the transported host value, and the concrete host value becomes the realized
    /// transaction pre-state. Added overlays retain their rule-owned frame because they have no
    /// host counterpart.
    pub fn apply_at(
        &self,
        host: &Molecule,
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
            let Some(right) = correspondence.atoms().right_of(left) else {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            };
            if right.index() >= host.atoms().count() || !host_atoms.insert(right) {
                return Err(ApplyError::CorrespondenceMismatch { entity });
            }
        }

        let Some(induced) =
            MoleculeCorrespondence::induce(&self.lhs, host, correspondence.atoms().clone())
        else {
            return Err(ApplyError::CorrespondenceMismatch {
                entity: Entity::Atom(AtomId(0)),
            });
        };
        macro_rules! require_induced_entity_correspondence {
            ($entity_set:ident, $entity:ident, $fallback:expr) => {{
                let supplied = correspondence.$entity_set();
                let expected = induced.$entity_set();
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
        require_induced_entity_correspondence!(bonds, Bond, BondId(0));
        require_induced_entity_correspondence!(dative_bonds, DativeBond, DativeBondId(0));
        require_induced_entity_correspondence!(
            aromatic_systems,
            AromaticSystem,
            AromaticSystemId(0)
        );
        require_induced_entity_correspondence!(
            multicenter_bonds,
            MulticenterBond,
            MulticenterBondId(0)
        );
        require_induced_entity_correspondence!(
            noncovalent_bonds,
            NoncovalentBond,
            NoncovalentBondId(0)
        );

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
            if correspondence.atoms().right_of(rule_site) != Some(host_view.site_id()) {
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

        let deltas = normalize_reaction_deltas(&self.lhs, &self.deltas)?;
        self.apply_at_canonical(host, correspondence, deltas)
    }

    fn apply_at_canonical(
        &self,
        host: &Molecule,
        correspondence: &MoleculeCorrespondence,
        mut deltas: Deltas,
    ) -> Result<ReactionDerivation, ApplyError> {
        deltas = reframe_application_deltas(deltas, &self.lhs, host, correspondence)?;
        let host_atom = |id: AtomId| {
            correspondence
                .atoms()
                .right_of(id)
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

        let mut created_atoms: BTreeMap<AtomId, AtomForm> = BTreeMap::new();
        let mut created_bonds: BTreeMap<BondId, ([AtomId; 2], BondForm)> = BTreeMap::new();
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
                Delta::Atom(AtomDelta::Add { id, attributes }) => {
                    created_atoms.insert(*id, attributes.clone());
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
                    sets.update_atom(
                        AtomHandle::Id(host_id),
                        host.atom(host_id).attributes,
                        &update,
                    );
                }
                Delta::Atom(AtomDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_atom(*id)?;
                        sets.update_atom(
                            AtomHandle::Id(host_id),
                            host.atom(host_id).attributes,
                            &AtomUpdate {
                                constraints: constraint.into(),
                                ..Default::default()
                            },
                        );
                    }
                }
                Delta::Bond(BondDelta::Add {
                    id,
                    atoms,
                    attributes,
                }) => {
                    created_bonds.insert(*id, (*atoms, attributes.clone()));
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
                    sets.update_bond(
                        BondHandle::Id(host_id),
                        host.bond(host_id).attributes,
                        &update,
                    );
                }
                Delta::Bond(BondDelta::ModifyConstraint { id, old, new }) => {
                    let constraint = new
                        .clone()
                        .or_else(|| old.as_ref().map(|constraint| constraint.as_undetermined()));
                    if let Some(constraint) = constraint {
                        let host_id = host_bond(*id)?;
                        sets.update_bond(
                            BondHandle::Id(host_id),
                            host.bond(host_id).attributes,
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
                            host.dative_bond(host_id).attributes,
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
                                host.dative_bond(host_id).attributes,
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
                            host.aromatic_system(host_id).attributes,
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
                                host.aromatic_system(host_id).attributes,
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
                            host.multicenter_bond(host_id).attributes,
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
                                host.multicenter_bond(host_id).attributes,
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
                            host.noncovalent_bond(host_id).attributes,
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
                                host.noncovalent_bond(host_id).attributes,
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
                // Stereo modifications lower directly after reframing into the matched host frame.
                // `Add` is lowered in the second pass; `Remove` tracks the host id for the DPO
                // dangling check.
                Delta::StereoAtom(s) => match s {
                    StereoAtomDelta::ModifyField { id, change } => {
                        let host_id = host_stereo_atom(*id)?;
                        let StereoAtomFieldChange::Configuration { new, .. } = change;
                        sets.push(Edit::ModifyStereoAtomField {
                            id: StereoAtomHandle::Id(host_id),
                            change: StereoAtomFieldChange::Configuration {
                                old: host.stereo_atom(host_id).attributes.configuration.clone(),
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
                                    .attributes
                                    .constraints
                                    .get(constraint.key())
                                    .cloned(),
                                new: new.clone(),
                            })
                        }
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
                                old: host.stereo_bond(host_id).attributes.configuration.clone(),
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
                                    .attributes
                                    .constraints
                                    .get(constraint.key())
                                    .cloned(),
                                new: new.clone(),
                            })
                        }
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
        let mut remove_dative: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondForm)> =
            Vec::new();
        let mut remove_aromatic: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemForm)> =
            Vec::new();
        let mut remove_multicenter: Vec<(
            MulticenterBondHandle,
            Vec<AtomHandle>,
            MulticenterBondForm,
        )> = Vec::new();
        let mut remove_noncovalent: Vec<(
            NoncovalentBondHandle,
            [AtomHandle; 2],
            NoncovalentBondForm,
        )> = Vec::new();
        let mut remove_stereo_atom: Vec<StereoAtomRemoval> = Vec::new();
        let mut remove_stereo_bond: Vec<StereoBondRemoval> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::DativeBond(DativeBondDelta::Add {
                    id,
                    donors,
                    acceptor,
                    attributes,
                }) => {
                    let mut atoms: Vec<AtomHandle> = donors
                        .iter()
                        .map(|a| atom_handle(*a))
                        .collect::<Result<_, _>>()?;
                    atoms.push(atom_handle(*acceptor)?);
                    let handle = overlay_adds.add_dative_bond(atoms, attributes.clone());
                    new_dative_handles.insert(*id, handle);
                }
                Delta::DativeBond(DativeBondDelta::Remove {
                    id,
                    donors,
                    acceptor,
                    attributes,
                }) => {
                    let mut atoms: Vec<AtomHandle> = donors
                        .iter()
                        .map(|a| atom_handle(*a))
                        .collect::<Result<_, _>>()?;
                    atoms.push(atom_handle(*acceptor)?);
                    remove_dative.push((
                        DativeBondHandle::Id(host_dative(*id)?),
                        atoms,
                        attributes.clone(),
                    ));
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add {
                    id,
                    atoms,
                    attributes,
                }) => {
                    let handle = overlay_adds.add_aromatic_system(
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    );
                    new_aromatic_handles.insert(*id, handle);
                }
                Delta::AromaticSystem(AromaticSystemDelta::Remove {
                    id,
                    atoms,
                    attributes,
                }) => {
                    remove_aromatic.push((
                        AromaticSystemHandle::Id(host_aromatic(*id)?),
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    ));
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add {
                    id,
                    atoms,
                    attributes,
                }) => {
                    let handle = overlay_adds.add_multicenter_bond(
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    );
                    new_multicenter_handles.insert(*id, handle);
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove {
                    id,
                    atoms,
                    attributes,
                }) => {
                    remove_multicenter.push((
                        MulticenterBondHandle::Id(host_multicenter(*id)?),
                        atoms
                            .iter()
                            .map(|a| atom_handle(*a))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    ));
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id,
                    atoms,
                    attributes,
                }) => {
                    let handle = overlay_adds.add_noncovalent_bond(
                        [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                        attributes.clone(),
                    );
                    new_noncovalent_handles.insert(*id, handle);
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id,
                    atoms,
                    attributes,
                }) => {
                    remove_noncovalent.push((
                        NoncovalentBondHandle::Id(host_noncovalent(*id)?),
                        [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                        attributes.clone(),
                    ));
                }
                Delta::StereoAtom(StereoAtomDelta::Add {
                    id,
                    site,
                    ligands,
                    attributes,
                }) => {
                    let handle = overlay_adds.add_stereo_atom(
                        atom_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    );
                    new_stereo_atom_handles.insert(*id, handle);
                }
                Delta::StereoAtom(StereoAtomDelta::Remove {
                    id,
                    site,
                    ligands,
                    attributes,
                }) => {
                    remove_stereo_atom.push((
                        StereoAtomHandle::Id(host_stereo_atom(*id)?),
                        atom_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    ));
                }
                Delta::StereoBond(StereoBondDelta::Add {
                    id,
                    site,
                    ligands,
                    attributes,
                }) => {
                    let handle = overlay_adds.add_stereo_bond(
                        bond_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
                    );
                    new_stereo_bond_handles.insert(*id, handle);
                }
                Delta::StereoBond(StereoBondDelta::Remove {
                    id,
                    site,
                    ligands,
                    attributes,
                }) => {
                    remove_stereo_bond.push((
                        StereoBondHandle::Id(host_stereo_bond(*id)?),
                        bond_handle(*site)?,
                        ligands
                            .iter()
                            .map(|l| atom_handle(l.atom_id).map(|atom| (atom, l.kind)))
                            .collect::<Result<_, _>>()?,
                        attributes.clone(),
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
                    .map(|(atoms, attributes)| {
                        Ok(AddBond {
                            endpoints: [atom_handle(atoms[0])?, atom_handle(atoms[1])?],
                            attributes: attributes.clone(),
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
        let product = match builder.try_build() {
            Ok(product) => product,
            Err(
                MoleculeIntegrityError::DuplicateParticipant { .. }
                | MoleculeIntegrityError::BondsParallel { .. }
                | MoleculeIntegrityError::DativeBondsIdentical { .. }
                | MoleculeIntegrityError::NoncovalentBondsParallel { .. }
                | MoleculeIntegrityError::AromaticSystemsOverlap { .. }
                | MoleculeIntegrityError::MulticenterBondsIdentical { .. }
                | MoleculeIntegrityError::StereoAtomSitesDuplicate { .. }
                | MoleculeIntegrityError::StereoBondSitesDuplicate { .. }
                | MoleculeIntegrityError::StereoLigandIncidenceMismatch { .. },
            ) => {
                return Err(ApplyError::StructuralConflict);
            }
            Err(_) => return Err(ApplyError::InternalInvariant),
        };

        // The host↔product comap: preserved host atoms match their compacted product id (survivors
        // keep ascending order), removed atoms are left-unmatched, created atoms right-unmatched.
        // `induce` derives the bond and overlay correspondences from this atom map.
        let removed: HashSet<AtomId> = removed_host_atoms.iter().copied().collect();
        let mut atom_matched_pairs: Vec<(AtomId, AtomId)> = Vec::new();
        let mut product_atom = 0u32;
        for host_atom in 0..host.atoms().count() as u32 {
            if removed.contains(&AtomId(host_atom)) {
                continue;
            }
            atom_matched_pairs.push((AtomId(host_atom), AtomId(product_atom)));
            product_atom += 1;
        }
        let atom_map = Correspondence::new(
            atom_matched_pairs,
            host.atoms().count(),
            product.atoms().count(),
        )
        .expect("correspondence producer preserves partial-bijection invariants");
        let comap = MoleculeCorrespondence::induce(host, &product, atom_map)
            .expect("successful reaction application preserves unique entity incidence");
        Ok(ReactionDerivation::new(host.clone(), product, comap))
    }

    /// Check the reaction-local structural preconditions shared by every application.
    pub fn check_preconditions(&self) -> Result<(), ApplyPreconditionError> {
        self.application_deltas().map(drop)
    }

    fn application_deltas(&self) -> Result<Deltas, ApplyPreconditionError> {
        if !stereo_delta_domains_are_valid(&self.lhs, &self.deltas) {
            return Err(ApplyPreconditionError::InconsistentReaction);
        }
        let deltas = normalize_reaction_deltas(&self.lhs, &self.deltas)
            .map_err(|_| ApplyPreconditionError::InconsistentReaction)?;

        check_reaction_dpo(&self.lhs, &deltas).map_err(ApplyPreconditionError::ReactionDpo)?;

        Ok(deltas)
    }

    /// Every derivation of applying the reaction to `host`: one per injective match of `lhs` into
    /// `host` under `match_config` that satisfies the match-local DPO and structural conditions.
    ///
    /// The returned iterator owns snapshots of this reaction and `host`. Matching is eager;
    /// derivation construction is lazy and follows match order. Match-local rejection is skipped;
    /// another application failure is yielded once and terminates the iterator.
    ///
    /// # Errors
    ///
    /// Returns [`ApplyPreconditionError`] before match enumeration when this reaction fails a
    /// reaction-wide structural precondition. Failures arising while realizing a selected match
    /// remain [`ApplyError`] iterator items.
    pub fn apply(
        &self,
        host: &Molecule,
        match_config: SubstructureMatchConfig,
    ) -> Result<ReactionApplicationIter, ApplyPreconditionError> {
        ReactionApplicationIter::new(self.clone(), host.clone(), match_config)
    }
}

fn application_frame_actions(
    lhs: &Molecule,
    deltas: &Deltas,
    host: &Molecule,
    correspondence: &MoleculeCorrespondence,
) -> Result<ConstraintFrameActionMap, ApplyError> {
    let domain = ReactionFrameActionDomain::from_deltas(deltas);
    let mut actions = ConstraintFrameActionMap::default();
    let mut dative_adds = HashMap::new();
    let mut aromatic_adds = HashMap::new();
    let mut multicenter_adds = HashMap::new();
    let mut noncovalent_adds = HashSet::new();
    let mut stereo_atom_adds = HashMap::new();
    let mut stereo_bond_adds = HashMap::new();
    for delta in deltas.iter() {
        match delta {
            Delta::DativeBond(DativeBondDelta::Add { id, donors, .. })
                if domain.constraints.contains_dative_bond(*id) =>
            {
                dative_adds.insert(*id, donors.len());
            }
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. })
                if domain.constraints.contains_aromatic_system(*id) =>
            {
                aromatic_adds.insert(*id, atoms.len());
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. })
                if domain.constraints.contains_multicenter_bond(*id) =>
            {
                multicenter_adds.insert(*id, atoms.len());
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, .. })
                if domain.constraints.contains_noncovalent_bond(*id) =>
            {
                noncovalent_adds.insert(*id);
            }
            Delta::StereoAtom(StereoAtomDelta::Add { id, ligands, .. })
                if domain.constraints.contains_stereo_atom(*id) =>
            {
                stereo_atom_adds.insert(*id, ligands.len());
            }
            Delta::StereoBond(StereoBondDelta::Add { id, ligands, .. })
                if domain.constraints.contains_stereo_bond(*id) =>
            {
                stereo_bond_adds.insert(*id, ligands.len());
            }
            _ => {}
        }
    }

    let mapped_atom = |atom: AtomId, entity: Entity| {
        correspondence
            .atoms()
            .right_of(atom)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })
    };
    let mapped_ligands = |ligands: &[StereoLigand], entity: Entity| {
        ligands
            .iter()
            .map(|ligand| {
                mapped_atom(ligand.atom_id, entity)
                    .map(|atom_id| StereoLigand::new(atom_id, ligand.kind))
            })
            .collect::<Result<Vec<_>, _>>()
    };

    for rule_view in lhs
        .dative_bonds()
        .iter()
        .filter(|view| domain.contains_dative_bond(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::DativeBond(id);
        let host_id = correspondence
            .dative_bonds()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .dative_bonds()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = rule_view
            .donor_ids()
            .map(|atom| mapped_atom(atom, entity))
            .collect::<Result<Vec<_>, _>>()?;
        let target = host_view.donor_ids().collect::<Vec<_>>();
        let action = DynPermutation::between(&mapped, &target)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        actions.insert_dative_bond(id, action);
    }
    for (id, degree) in dative_adds {
        if domain.constraints.contains_dative_bond(id) && lhs.dative_bonds().get(id).is_none() {
            actions.insert_dative_bond(id, DynPermutation::identity(degree));
        }
    }

    for rule_view in lhs
        .aromatic_systems()
        .iter()
        .filter(|view| domain.contains_aromatic_system(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::AromaticSystem(id);
        let host_id = correspondence
            .aromatic_systems()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .aromatic_systems()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = rule_view
            .atom_ids()
            .map(|atom| mapped_atom(atom, entity))
            .collect::<Result<Vec<_>, _>>()?;
        let target = host_view.atom_ids().collect::<Vec<_>>();
        let action = DynPermutation::between(&mapped, &target)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        actions.insert_aromatic_system(id, action);
    }
    for (id, degree) in aromatic_adds {
        if domain.constraints.contains_aromatic_system(id)
            && lhs.aromatic_systems().get(id).is_none()
        {
            actions.insert_aromatic_system(id, DynPermutation::identity(degree));
        }
    }

    for rule_view in lhs
        .multicenter_bonds()
        .iter()
        .filter(|view| domain.contains_multicenter_bond(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::MulticenterBond(id);
        let host_id = correspondence
            .multicenter_bonds()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .multicenter_bonds()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = rule_view
            .atom_ids()
            .map(|atom| mapped_atom(atom, entity))
            .collect::<Result<Vec<_>, _>>()?;
        let target = host_view.atom_ids().collect::<Vec<_>>();
        let action = DynPermutation::between(&mapped, &target)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        actions.insert_multicenter_bond(id, action);
    }
    for (id, degree) in multicenter_adds {
        if domain.constraints.contains_multicenter_bond(id)
            && lhs.multicenter_bonds().get(id).is_none()
        {
            actions.insert_multicenter_bond(id, DynPermutation::identity(degree));
        }
    }

    for rule_view in lhs
        .noncovalent_bonds()
        .iter()
        .filter(|view| domain.contains_noncovalent_bond(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::NoncovalentBond(id);
        let host_id = correspondence
            .noncovalent_bonds()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .noncovalent_bonds()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = rule_view
            .atom_ids()
            .into_iter()
            .map(|atom| mapped_atom(atom, entity))
            .collect::<Result<Vec<_>, _>>()?;
        let action = DynPermutation::between(&mapped, &host_view.atom_ids())
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        actions.insert_noncovalent_bond(id, action);
    }
    for id in noncovalent_adds {
        if domain.constraints.contains_noncovalent_bond(id)
            && lhs.noncovalent_bonds().get(id).is_none()
        {
            actions.insert_noncovalent_bond(id, DynPermutation::identity(2));
        }
    }

    for rule_view in lhs
        .stereo_atoms()
        .iter()
        .filter(|view| domain.contains_stereo_atom(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::StereoAtom(id);
        let host_id = correspondence
            .stereo_atoms()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .stereo_atoms()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = mapped_ligands(&rule_view.ligand_frame(), entity)?;
        let action = Permutation::between(&mapped, &host_view.ligand_frame())
            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
        actions.insert_stereo_atom(id, action);
    }
    for (id, degree) in stereo_atom_adds {
        if domain.constraints.contains_stereo_atom(id) && lhs.stereo_atoms().get(id).is_none() {
            actions.insert_stereo_atom(id, Permutation::identity(degree));
        }
    }

    for rule_view in lhs
        .stereo_bonds()
        .iter()
        .filter(|view| domain.contains_stereo_bond(view.id))
    {
        let id = rule_view.id;
        let entity = Entity::StereoBond(id);
        let host_id = correspondence
            .stereo_bonds()
            .right_of(id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let host_view = host
            .stereo_bonds()
            .get(host_id)
            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
        let mapped = mapped_ligands(&rule_view.ligand_frame(), entity)?;
        let action = Permutation::between(&mapped, &host_view.ligand_frame())
            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
        actions.insert_stereo_bond(id, action);
    }
    for (id, degree) in stereo_bond_adds {
        if domain.constraints.contains_stereo_bond(id) && lhs.stereo_bonds().get(id).is_none() {
            actions.insert_stereo_bond(id, Permutation::identity(degree));
        }
    }

    Ok(actions)
}

fn optional_pattern_matches<T: Lattice>(pattern: &Option<T>, target: &Option<T>) -> bool {
    match (pattern, target) {
        (None, _) => true,
        (Some(pattern), Some(target)) => pattern.matches(target),
        (Some(pattern), None) => pattern.is_undetermined(),
    }
}

fn reframe_application_deltas(
    deltas: Deltas,
    lhs: &Molecule,
    host: &Molecule,
    correspondence: &MoleculeCorrespondence,
) -> Result<Deltas, ApplyError> {
    let actions = application_frame_actions(lhs, &deltas, host, correspondence)?;
    deltas
        .into_iter()
        .map(|delta| {
            Ok(match delta {
                Delta::DativeBond(mut delta) => {
                    let id = match &delta {
                        DativeBondDelta::Add { id, .. }
                        | DativeBondDelta::Remove { id, .. }
                        | DativeBondDelta::ModifyField { id, .. }
                        | DativeBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::DativeBond(id);
                    if matches!(delta, DativeBondDelta::Add { .. }) {
                        return Ok(Delta::DativeBond(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .dative_bond_action(id)
                            .ok_or(ApplyError::InternalInvariant)?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::InternalInvariant)?;
                    }
                    if let DativeBondDelta::Remove { attributes, .. } = &mut delta {
                        let host_id = correspondence
                            .dative_bonds()
                            .right_of(id)
                            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                        let host_attributes = host
                            .dative_bonds()
                            .get(host_id)
                            .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                            .attributes;
                        if !attributes.matches(host_attributes) {
                            return Err(ApplyError::CorrespondenceMismatch { entity });
                        }
                        *attributes = host_attributes.clone();
                    }
                    Delta::DativeBond(delta)
                }
                Delta::AromaticSystem(mut delta) => {
                    let id = match &delta {
                        AromaticSystemDelta::Add { id, .. }
                        | AromaticSystemDelta::Remove { id, .. }
                        | AromaticSystemDelta::ModifyField { id, .. }
                        | AromaticSystemDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::AromaticSystem(id);
                    if matches!(delta, AromaticSystemDelta::Add { .. }) {
                        return Ok(Delta::AromaticSystem(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .aromatic_system_action(id)
                            .ok_or(ApplyError::InternalInvariant)?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::InternalInvariant)?;
                    }
                    let host_id = correspondence
                        .aromatic_systems()
                        .right_of(id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                    let host_attributes = host
                        .aromatic_systems()
                        .get(host_id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                        .attributes;
                    match &mut delta {
                        AromaticSystemDelta::Remove { attributes, .. } => {
                            if !attributes.matches(host_attributes) {
                                return Err(ApplyError::CorrespondenceMismatch { entity });
                            }
                            *attributes = host_attributes.clone();
                        }
                        AromaticSystemDelta::ModifyField {
                            change: AromaticSystemFieldChange::Electrons { old, .. },
                            ..
                        } => {
                            if !old.matches(&host_attributes.electrons) {
                                return Err(ApplyError::CorrespondenceMismatch { entity });
                            }
                            *old = host_attributes.electrons.clone();
                        }
                        _ => {}
                    }
                    Delta::AromaticSystem(delta)
                }
                Delta::MulticenterBond(mut delta) => {
                    let id = match &delta {
                        MulticenterBondDelta::Add { id, .. }
                        | MulticenterBondDelta::Remove { id, .. }
                        | MulticenterBondDelta::ModifyField { id, .. }
                        | MulticenterBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::MulticenterBond(id);
                    if matches!(delta, MulticenterBondDelta::Add { .. }) {
                        return Ok(Delta::MulticenterBond(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .multicenter_bond_action(id)
                            .ok_or(ApplyError::InternalInvariant)?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::InternalInvariant)?;
                    }
                    let host_id = correspondence
                        .multicenter_bonds()
                        .right_of(id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                    let host_attributes = host
                        .multicenter_bonds()
                        .get(host_id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                        .attributes;
                    match &mut delta {
                        MulticenterBondDelta::Remove { attributes, .. } => {
                            if !attributes.matches(host_attributes) {
                                return Err(ApplyError::CorrespondenceMismatch { entity });
                            }
                            *attributes = host_attributes.clone();
                        }
                        MulticenterBondDelta::ModifyField {
                            change: MulticenterBondFieldChange::Electrons { old, .. },
                            ..
                        } => {
                            if !old.matches(&host_attributes.electrons) {
                                return Err(ApplyError::CorrespondenceMismatch { entity });
                            }
                            *old = host_attributes.electrons.clone();
                        }
                        _ => {}
                    }
                    Delta::MulticenterBond(delta)
                }
                Delta::NoncovalentBond(mut delta) => {
                    let id = match &delta {
                        NoncovalentBondDelta::Add { id, .. }
                        | NoncovalentBondDelta::Remove { id, .. }
                        | NoncovalentBondDelta::ModifyField { id, .. }
                        | NoncovalentBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::NoncovalentBond(id);
                    if matches!(delta, NoncovalentBondDelta::Add { .. }) {
                        return Ok(Delta::NoncovalentBond(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .noncovalent_bond_action(id)
                            .ok_or(ApplyError::InternalInvariant)?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::InternalInvariant)?;
                    }
                    if let NoncovalentBondDelta::Remove { attributes, .. } = &mut delta {
                        let host_id = correspondence
                            .noncovalent_bonds()
                            .right_of(id)
                            .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                        let host_attributes = host
                            .noncovalent_bonds()
                            .get(host_id)
                            .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                            .attributes;
                        if !attributes.matches(host_attributes) {
                            return Err(ApplyError::CorrespondenceMismatch { entity });
                        }
                        *attributes = host_attributes.clone();
                    }
                    Delta::NoncovalentBond(delta)
                }
                Delta::StereoAtom(mut delta) => {
                    let id = delta.id();
                    let entity = Entity::StereoAtom(id);
                    if matches!(delta, StereoAtomDelta::Add { .. }) {
                        return Ok(Delta::StereoAtom(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .stereo_atom_action(id)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                    }
                    let host_id = correspondence
                        .stereo_atoms()
                        .right_of(id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                    let host_attributes = host
                        .stereo_atoms()
                        .get(host_id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                        .attributes;
                    match &mut delta {
                        StereoAtomDelta::Remove { attributes, .. } => {
                            if !attributes.matches(host_attributes) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *attributes = host_attributes.clone();
                        }
                        StereoAtomDelta::ModifyField {
                            change: StereoAtomFieldChange::Configuration { old, .. },
                            ..
                        } => {
                            if !old.matches(&host_attributes.configuration) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *old = host_attributes.configuration.clone();
                        }
                        StereoAtomDelta::ModifyConstraint { old, new, .. } => {
                            let key = old.as_ref().or(new.as_ref()).map(|value| value.key());
                            let host_old =
                                key.and_then(|key| host_attributes.constraints.get(key).cloned());
                            if !optional_pattern_matches(old, &host_old) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *old = host_old;
                        }
                        StereoAtomDelta::Add { .. } => unreachable!(),
                    }
                    Delta::StereoAtom(delta)
                }
                Delta::StereoBond(mut delta) => {
                    let id = delta.id();
                    let entity = Entity::StereoBond(id);
                    if matches!(delta, StereoBondDelta::Add { .. }) {
                        return Ok(Delta::StereoBond(delta));
                    }
                    if delta.uses_participant_frame() {
                        let action = actions
                            .stereo_bond_action(id)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                        delta = delta
                            .reframe_by(action)
                            .ok_or(ApplyError::StereoFrameMismatch { entity })?;
                    }
                    let host_id = correspondence
                        .stereo_bonds()
                        .right_of(id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?;
                    let host_attributes = host
                        .stereo_bonds()
                        .get(host_id)
                        .ok_or(ApplyError::CorrespondenceMismatch { entity })?
                        .attributes;
                    match &mut delta {
                        StereoBondDelta::Remove { attributes, .. } => {
                            if !attributes.matches(host_attributes) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *attributes = host_attributes.clone();
                        }
                        StereoBondDelta::ModifyField {
                            change: StereoBondFieldChange::Configuration { old, .. },
                            ..
                        } => {
                            if !old.matches(&host_attributes.configuration) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *old = host_attributes.configuration.clone();
                        }
                        StereoBondDelta::ModifyConstraint { old, new, .. } => {
                            let key = old.as_ref().or(new.as_ref()).map(|value| value.key());
                            let host_old =
                                key.and_then(|key| host_attributes.constraints.get(key).cloned());
                            if !optional_pattern_matches(old, &host_old) {
                                return Err(ApplyError::StereoFrameMismatch { entity });
                            }
                            *old = host_old;
                        }
                        StereoBondDelta::Add { .. } => unreachable!(),
                    }
                    Delta::StereoBond(delta)
                }
                Delta::Constraint(delta) => {
                    let delta =
                        delta
                            .reframe_by_actions(&actions)
                            .map_err(|entity| match entity {
                                Entity::StereoAtom(_) | Entity::StereoBond(_) => {
                                    ApplyError::StereoFrameMismatch { entity }
                                }
                                _ => ApplyError::InternalInvariant,
                            })?;
                    Delta::Constraint(delta)
                }
                invariant => invariant,
            })
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Deltas::from_iter)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};
    use umol_perm::MAX_DEGREE;

    use super::super::constraint::{
        AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
        Constraints, DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
        NoncovalentBondConstraintForm, RelationalConstraint, StereoAtomConstraintForm,
        StereoBondConstraintForm, StereoLigandPair, StereogenicityForm, TopicityForm,
        TopicityRelationForm,
    };
    use super::super::edit::{AtomFieldChange, BondFieldChange};
    use super::super::electrons::ElectronCountsForm;
    use super::super::entity::{Entity, EntityKind};
    use super::super::id::StereoLigandPosition;
    use super::super::ligand::StereoLigandKind;
    use super::super::molecule::transact::TransactionError;
    use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use super::super::num::NumForm;
    use super::super::stereo::{StereoAtomForm, StereoBondForm, StereoCoset, StereoKind, Topicity};
    use super::super::substructure::SubstructureMatchAlgorithm;
    use super::*;

    const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
        match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    };

    fn removal_frame_reactions(entity_kind: EntityKind) -> (Reaction, Reaction) {
        let atom_ligands: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let bond_ligands: Vec<StereoLigand> = (0..=3)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(2)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(5), AtomId(1), BondForm::from_order(1)),
                (AtomId(6), AtomId(2), BondForm::from_order(1)),
                (AtomId(6), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        };
        let (owner_removal, local_removal) = match entity_kind {
            EntityKind::DativeBond => {
                let attributes = DativeBondForm::from_order(2);
                entries
                    .dative
                    .push((vec![AtomId(0), AtomId(1)], AtomId(2), attributes.clone()));
                (
                    Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: vec![AtomId(0), AtomId(1)],
                        acceptor: AtomId(2),
                        attributes: attributes.clone(),
                    }),
                    Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: vec![AtomId(1), AtomId(0)],
                        acceptor: AtomId(2),
                        attributes,
                    }),
                )
            }
            EntityKind::AromaticSystem => {
                let attributes = AromaticSystemForm::from_electrons(vec![1, 2]);
                entries
                    .aromatic
                    .push((vec![AtomId(0), AtomId(1)], attributes.clone()));
                (
                    Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: attributes.clone(),
                    }),
                    Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(1), AtomId(0)],
                        attributes: attributes
                            .reframe_by(
                                &DynPermutation::between(
                                    &[AtomId(0), AtomId(1)],
                                    &[AtomId(1), AtomId(0)],
                                )
                                .expect("the frames have the same atoms"),
                            )
                            .expect("the action has the form's degree"),
                    }),
                )
            }
            EntityKind::MulticenterBond => {
                let attributes = MulticenterBondForm::from_electrons(vec![3, 4]);
                entries
                    .multicenter
                    .push((vec![AtomId(0), AtomId(1)], attributes.clone()));
                (
                    Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: attributes.clone(),
                    }),
                    Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(1), AtomId(0)],
                        attributes: attributes
                            .reframe_by(
                                &DynPermutation::between(
                                    &[AtomId(0), AtomId(1)],
                                    &[AtomId(1), AtomId(0)],
                                )
                                .expect("the frames have the same atoms"),
                            )
                            .expect("the action has the form's degree"),
                    }),
                )
            }
            EntityKind::NoncovalentBond => {
                let attributes = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);
                entries
                    .noncovalent
                    .push((AtomId(0), AtomId(1), attributes.clone()));
                (
                    Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: attributes.clone(),
                    }),
                    Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(1), AtomId(0)],
                        attributes,
                    }),
                )
            }
            EntityKind::StereoAtom => {
                let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
                entries
                    .stereo_atoms
                    .push((AtomId(0), atom_ligands.clone(), attributes.clone()));
                let mut local_ligands = atom_ligands.clone();
                local_ligands.swap(0, 1);
                let local_attributes = attributes
                    .clone()
                    .reframe_by(
                        &Permutation::between(&atom_ligands, &local_ligands)
                            .expect("the frames have the same ligands"),
                    )
                    .expect("the local frame is a tetrahedral action");
                (
                    Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: atom_ligands,
                        attributes,
                    }),
                    Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: local_ligands,
                        attributes: local_attributes,
                    }),
                )
            }
            EntityKind::StereoBond => {
                let attributes = StereoBondForm::new(StereoKind::CisTrans, 0u32);
                entries
                    .stereo_bonds
                    .push((BondId(4), bond_ligands.clone(), attributes.clone()));
                let mut local_ligands = bond_ligands.clone();
                local_ligands.swap(0, 1);
                let local_attributes = attributes
                    .clone()
                    .reframe_by(
                        &Permutation::between(&bond_ligands, &local_ligands)
                            .expect("the frames have the same ligands"),
                    )
                    .expect("the local frame preserves the stereo-bond endpoint blocks");
                (
                    Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(4),
                        ligands: bond_ligands,
                        attributes,
                    }),
                    Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(4),
                        ligands: local_ligands,
                        attributes: local_attributes,
                    }),
                )
            }
            EntityKind::Atom | EntityKind::Bond => {
                unreachable!("only overlay entity kinds have participant frames")
            }
        };
        let lhs = Molecule::from_entries(entries);
        (
            Reaction::new(lhs.clone(), Deltas::from_iter([owner_removal])),
            Reaction::new(lhs, Deltas::from_iter([local_removal])),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::successful(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) } })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() })],
    )]
    #[case::match_rejection(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))], ..Default::default() }),
        vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() })],
    )]
    fn test_reaction_application_iter(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] expected: Vec<Molecule>,
    ) {
        let mut applications =
            ReactionApplicationIter::new(reaction, host, MATCH_CONFIG).unwrap();
        let products = applications
            .by_ref()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect::<Vec<_>>();

        assert_eq!(products, expected);
        assert_eq!(applications.next(), None);
    }

    #[rstest]
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    fn test_reaction_application_iter_error() {
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                constraints: Constraints::from(Constraint::Molecule(
                    MoleculeConstraint::ChargeSum {
                        atoms: Some(vec![AtomId(0)]),
                        sum: NumForm::Lit(0),
                    },
                )),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0)]),
                    sum: NumForm::Lit(0),
                }),
            ))]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let mut applications = ReactionApplicationIter::new(reaction, host, MATCH_CONFIG).unwrap();

        assert_eq!(
            applications.next(),
            Some(Err(ApplyError::Transaction(TransactionError::MissingEntry))),
        );
        assert_eq!(applications.next(), None);
    }

    #[rstest]
    fn test_reaction_application_iter_snapshot() {
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = Reaction::new(
            host.clone(),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Undetermined,
                    new: NumForm::Lit(1),
                },
            })]),
        );
        let mut applications =
            ReactionApplicationIter::new(reaction.clone(), host.clone(), MATCH_CONFIG).unwrap();
        let mut changed_reaction = reaction;
        changed_reaction.deltas = Deltas::new();
        let expected_changed_reaction = Reaction::new(changed_reaction.lhs.clone(), Deltas::new());
        let mut changed_host = host;
        changed_host.combine_from(&Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O)],
            ..Default::default()
        }));

        let derivation = applications.next().unwrap().unwrap();
        assert_eq!(
            derivation.rhs(),
            &Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C).with_charge(1_i64)],
                ..Default::default()
            }),
        );
        assert_eq!(applications.next(), None);
        assert_eq!(changed_reaction, expected_changed_reaction);
        assert_eq!(
            changed_host,
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                ..Default::default()
            }),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zero_components(
        Reaction::default(),
        Molecule::default(),
        vec![vec![]],
    )]
    #[case::one_component(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Deltas::new(),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
        vec![vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() })]],
    )]
    #[case::multiple_components(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], ..Default::default() }),
            Deltas::new(),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], ..Default::default() }),
        vec![vec![
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O)], ..Default::default() }),
        ]],
    )]
    fn test_reaction_products_iter(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] expected: Vec<Vec<Molecule>>,
    ) {
        let mut products = ReactionProductsIter {
            applications: ReactionApplicationIter::new(reaction, host, MATCH_CONFIG).unwrap(),
        };
        let actual = products
            .by_ref()
            .map(Result::unwrap)
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
        assert_eq!(products.next(), None);
    }

    #[rstest]
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    fn test_reaction_products_iter_error() {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        });
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                constraints: Constraints::from(constraint.clone()),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(constraint))]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let mut products = ReactionProductsIter {
            applications: ReactionApplicationIter::new(reaction, host, MATCH_CONFIG).unwrap(),
        };

        assert_eq!(
            products.next(),
            Some(Err(ApplyError::Transaction(TransactionError::MissingEntry))),
        );
        assert_eq!(products.next(), None);
    }

    #[rstest]
    fn test_react_react() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = Reaction::new(molecule.clone(), Deltas::new());

        let products = molecule
            .react(&reaction, MATCH_CONFIG)
            .unwrap()
            .map(Result::unwrap)
            .collect::<Vec<_>>();

        assert_eq!(products, vec![vec![molecule]]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::multiple(
        vec![
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O)], ..Default::default() }),
        ],
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], ..Default::default() }),
            Deltas::new(),
        ),
        vec![vec![
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O)], ..Default::default() }),
        ]],
    )]
    #[case::empty(
        vec![],
        Reaction::default(),
        vec![vec![]],
    )]
    fn test_react_react_slice(
        #[case] reactants: Vec<Molecule>,
        #[case] reaction: Reaction,
        #[case] expected: Vec<Vec<Molecule>>,
    ) {
        let products = reactants
            .react(&reaction, MATCH_CONFIG)
            .unwrap()
            .map(Result::unwrap)
            .collect::<Vec<_>>();

        assert_eq!(products, expected);
    }

    #[rstest]
    fn test_reaction_new() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let deltas = Deltas::new();

        let reaction = Reaction::new(lhs.clone(), deltas.clone());

        assert_eq!(reaction, Reaction { lhs, deltas });
    }

    #[rstest]
    #[should_panic(expected = "invalid reaction: InvalidReference { entity: Atom(AtomId(0)) }")]
    fn test_reaction_new_error() {
        Reaction::new(
            Molecule::default(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::default(),
            })]),
        );
    }

    #[rstest]
    fn test_reaction_try_new() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: [Element::C, Element::F, Element::Cl, Element::H, Element::H]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        });
        let deltas = Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: (1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        })]);

        let reaction = Reaction::try_new(lhs.clone(), deltas.clone()).unwrap();

        assert_eq!(reaction, Reaction { lhs, deltas });
    }

    #[rstest]
    #[case::atom(
        Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoAtomForm::default(),
        }),
        MoleculeIntegrityError::DuplicateStereoLigand {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            ligand: StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        },
    )]
    #[case::bond(
        Delta::StereoBond(StereoBondDelta::Add {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            ],
            attributes: StereoBondForm::default(),
        }),
        MoleculeIntegrityError::DuplicateStereoLigand {
            entity: Entity::StereoBond(StereoBondId(0)),
            ligand: StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        },
    )]
    #[case::atom_oversized(
        Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: (1..=(MAX_DEGREE as u32 + 1))
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoAtomForm::default(),
        }),
        MoleculeIntegrityError::StereoFrameDegreeTooLarge {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            degree: MAX_DEGREE + 1,
            maximum: MAX_DEGREE,
        },
    )]
    #[case::bond_oversized(
        Delta::StereoBond(StereoBondDelta::Add {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: (0..=(MAX_DEGREE as u32))
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoBondForm::default(),
        }),
        MoleculeIntegrityError::StereoFrameDegreeTooLarge {
            entity: Entity::StereoBond(StereoBondId(0)),
            degree: MAX_DEGREE + 1,
            maximum: MAX_DEGREE,
        },
    )]
    fn test_reaction_try_new_stereo_add_error(
        #[case] delta: Delta,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); MAX_DEGREE + 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });

        assert_eq!(
            Reaction::try_new(lhs, Deltas::from_iter([delta])),
            Err(ReactionIntegrityError::StereoIntegrityError(expected)),
        );
    }

    #[rstest]
    #[case::atom(
        Delta::StereoAtom(StereoAtomDelta::ModifyConstraint {
            id: StereoAtomId(0),
            kind: Some(StereoKind::CisTrans),
            old: None,
            new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        }),
        Entity::StereoAtom(StereoAtomId(0)),
        StereoKind::CisTrans,
    )]
    #[case::bond(
        Delta::StereoBond(StereoBondDelta::ModifyConstraint {
            id: StereoBondId(0),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        }),
        Entity::StereoBond(StereoBondId(0)),
        StereoKind::Tetrahedral,
    )]
    fn test_reaction_try_new_stereo_modify_constraint_error(
        #[case] delta: Delta,
        #[case] entity: Entity,
        #[case] kind: StereoKind,
    ) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(2)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(5), AtomId(1), BondForm::from_order(1)),
                (AtomId(6), AtomId(2), BondForm::from_order(1)),
                (AtomId(6), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::default(),
            )],
            stereo_bonds: vec![(
                BondId(4),
                (0..=3)
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoBondForm::default(),
            )],
            ..Default::default()
        });

        assert_eq!(
            Reaction::try_new(lhs, Deltas::from_iter([delta])),
            Err(ReactionIntegrityError::StereoIntegrityError(
                MoleculeIntegrityError::StereoKindSiteMismatch { entity, kind },
            )),
        );
    }

    #[rstest]
    #[case::atom_existing(EntityKind::StereoAtom, true)]
    #[case::atom_added(EntityKind::StereoAtom, false)]
    #[case::bond_existing(EntityKind::StereoBond, true)]
    #[case::bond_added(EntityKind::StereoBond, false)]
    fn test_reaction_try_new_stereo_constraint_wrapper_error(
        #[case] entity_kind: EntityKind,
        #[case] existing: bool,
    ) {
        let atom_ligands: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let bond_ligands: Vec<StereoLigand> = (0..=3)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(2)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(5), AtomId(1), BondForm::from_order(1)),
                (AtomId(6), AtomId(2), BondForm::from_order(1)),
                (AtomId(6), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        };
        let (entity, kind, addition, constraint) = match entity_kind {
            EntityKind::StereoAtom => {
                if existing {
                    entries.stereo_atoms.push((
                        AtomId(0),
                        atom_ligands.clone(),
                        StereoAtomForm::default(),
                    ));
                }
                (
                    Entity::StereoAtom(StereoAtomId(0)),
                    StereoKind::CisTrans,
                    Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: atom_ligands,
                        attributes: StereoAtomForm::default(),
                    }),
                    Constraint::StereoAtom(
                        StereoAtomId(0),
                        StereoKind::CisTrans,
                        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                    ),
                )
            }
            EntityKind::StereoBond => {
                if existing {
                    entries.stereo_bonds.push((
                        BondId(4),
                        bond_ligands.clone(),
                        StereoBondForm::default(),
                    ));
                }
                (
                    Entity::StereoBond(StereoBondId(0)),
                    StereoKind::Tetrahedral,
                    Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(0),
                        site: BondId(4),
                        ligands: bond_ligands,
                        attributes: StereoBondForm::default(),
                    }),
                    Constraint::StereoBond(
                        StereoBondId(0),
                        StereoKind::Tetrahedral,
                        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                    ),
                )
            }
            _ => unreachable!("test cases contain only stereo entity kinds"),
        };
        let lhs = Molecule::from_entries(entries);
        let mut deltas = Vec::new();
        if !existing {
            deltas.push(addition);
        }
        deltas.push(Delta::Constraint(ConstraintDelta::Add(constraint)));

        assert_eq!(
            Reaction::try_new(lhs, Deltas::from_iter(deltas)),
            Err(ReactionIntegrityError::StereoIntegrityError(
                MoleculeIntegrityError::StereoKindSiteMismatch { entity, kind },
            )),
        );
    }

    #[rstest]
    #[case::dative_existing(EntityKind::DativeBond, true)]
    #[case::dative_created(EntityKind::DativeBond, false)]
    #[case::aromatic_existing(EntityKind::AromaticSystem, true)]
    #[case::aromatic_created(EntityKind::AromaticSystem, false)]
    #[case::multicenter_existing(EntityKind::MulticenterBond, true)]
    #[case::multicenter_created(EntityKind::MulticenterBond, false)]
    #[case::noncovalent_existing(EntityKind::NoncovalentBond, true)]
    #[case::noncovalent_created(EntityKind::NoncovalentBond, false)]
    #[case::stereo_atom_existing(EntityKind::StereoAtom, true)]
    #[case::stereo_atom_created(EntityKind::StereoAtom, false)]
    #[case::stereo_bond_existing(EntityKind::StereoBond, true)]
    #[case::stereo_bond_created(EntityKind::StereoBond, false)]
    fn test_reaction_try_new_local_removal_frame(
        #[case] entity_kind: EntityKind,
        #[case] existing: bool,
    ) {
        let atom_ligands: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let bond_ligands: Vec<StereoLigand> = (0..=3)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(2)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(5), AtomId(1), BondForm::from_order(1)),
                (AtomId(6), AtomId(2), BondForm::from_order(1)),
                (AtomId(6), AtomId(3), BondForm::from_order(1)),
            ],
            ..Default::default()
        };
        let (_entity, addition, removal) = match entity_kind {
            EntityKind::DativeBond => {
                if existing {
                    entries.dative.push((
                        vec![AtomId(0), AtomId(1)],
                        AtomId(2),
                        DativeBondForm::default(),
                    ));
                }
                (
                    Entity::DativeBond(DativeBondId(0)),
                    Delta::DativeBond(DativeBondDelta::Add {
                        id: DativeBondId(0),
                        donors: vec![AtomId(0), AtomId(1)],
                        acceptor: AtomId(2),
                        attributes: DativeBondForm::default(),
                    }),
                    Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: vec![AtomId(1), AtomId(0)],
                        acceptor: AtomId(2),
                        attributes: DativeBondForm::default(),
                    }),
                )
            }
            EntityKind::AromaticSystem => {
                if existing {
                    entries
                        .aromatic
                        .push((vec![AtomId(0), AtomId(1)], AromaticSystemForm::default()));
                }
                (
                    Entity::AromaticSystem(AromaticSystemId(0)),
                    Delta::AromaticSystem(AromaticSystemDelta::Add {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: AromaticSystemForm::default(),
                    }),
                    Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(1), AtomId(0)],
                        attributes: AromaticSystemForm::default(),
                    }),
                )
            }
            EntityKind::MulticenterBond => {
                if existing {
                    entries
                        .multicenter
                        .push((vec![AtomId(0), AtomId(1)], MulticenterBondForm::default()));
                }
                (
                    Entity::MulticenterBond(MulticenterBondId(0)),
                    Delta::MulticenterBond(MulticenterBondDelta::Add {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: MulticenterBondForm::default(),
                    }),
                    Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(1), AtomId(0)],
                        attributes: MulticenterBondForm::default(),
                    }),
                )
            }
            EntityKind::NoncovalentBond => {
                if existing {
                    entries.noncovalent.push((
                        AtomId(0),
                        AtomId(1),
                        NoncovalentBondForm::default(),
                    ));
                }
                (
                    Entity::NoncovalentBond(NoncovalentBondId(0)),
                    Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: NoncovalentBondForm::default(),
                    }),
                    Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(1), AtomId(0)],
                        attributes: NoncovalentBondForm::default(),
                    }),
                )
            }
            EntityKind::StereoAtom => {
                if existing {
                    entries.stereo_atoms.push((
                        AtomId(0),
                        atom_ligands.clone(),
                        StereoAtomForm::default(),
                    ));
                }
                let mut removed_ligands = atom_ligands.clone();
                removed_ligands.swap(0, 1);
                (
                    Entity::StereoAtom(StereoAtomId(0)),
                    Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: atom_ligands,
                        attributes: StereoAtomForm::default(),
                    }),
                    Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: removed_ligands,
                        attributes: StereoAtomForm::default(),
                    }),
                )
            }
            EntityKind::StereoBond => {
                if existing {
                    entries.stereo_bonds.push((
                        BondId(4),
                        bond_ligands.clone(),
                        StereoBondForm::default(),
                    ));
                }
                let mut removed_ligands = bond_ligands.clone();
                removed_ligands.swap(0, 1);
                (
                    Entity::StereoBond(StereoBondId(0)),
                    Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(0),
                        site: BondId(4),
                        ligands: bond_ligands,
                        attributes: StereoBondForm::default(),
                    }),
                    Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(4),
                        ligands: removed_ligands,
                        attributes: StereoBondForm::default(),
                    }),
                )
            }
            EntityKind::Atom | EntityKind::Bond => {
                unreachable!("only overlay entity kinds have participant frames")
            }
        };
        let lhs = Molecule::from_entries(entries);
        let deltas = if existing {
            Deltas::from_iter([removal])
        } else {
            Deltas::from_iter([addition, removal])
        };

        Reaction::try_new(lhs, deltas).expect("the local removal frame is incidence-compatible");
    }

    #[rstest]
    #[case::dative(EntityKind::DativeBond)]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::multicenter(EntityKind::MulticenterBond)]
    #[case::noncovalent(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_reaction_local_removal_frame_semantics(#[case] entity_kind: EntityKind) {
        let (owner, local) = removal_frame_reactions(entity_kind);
        let host = owner.lhs().clone();
        let atom_ids: Vec<AtomId> = host.atoms().ids().collect();
        let correspondence = MoleculeCorrespondence::induce(
            owner.lhs(),
            &host,
            Correspondence::from_images(&atom_ids, host.atoms().count()),
        )
        .expect("identity atom images induce the identity molecule correspondence");

        assert_eq!(
            local.to_reaction_span().unwrap(),
            owner.to_reaction_span().unwrap(),
        );
        assert_eq!(
            local.apply_at(&host, &correspondence).unwrap(),
            owner.apply_at(&host, &correspondence).unwrap(),
        );
        assert_eq!(
            local.reverse().unwrap().to_reaction_span().unwrap(),
            owner.reverse().unwrap().to_reaction_span().unwrap(),
        );
    }

    #[rstest]
    #[case::dative(EntityKind::DativeBond)]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::multicenter(EntityKind::MulticenterBond)]
    #[case::noncovalent(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_reaction_normalize_removal(#[case] entity_kind: EntityKind) {
        let (owner, local) = removal_frame_reactions(entity_kind);

        assert_eq!(local.normalize(), owner.normalize());
    }

    #[rstest]
    #[case::dative(EntityKind::DativeBond)]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::multicenter(EntityKind::MulticenterBond)]
    #[case::noncovalent(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_reaction_reframe_by_identity(#[case] entity_kind: EntityKind) {
        let (_, local) = removal_frame_reactions(entity_kind);
        let identity = local.representative_action().identity();

        assert_eq!(local.clone().reframe_by(&identity), Some(local));
    }

    #[rstest]
    fn test_reaction_reframe_by_composition() {
        let owner_atoms = vec![AtomId(2), AtomId(0), AtomId(1)];
        let local_atoms = vec![AtomId(0), AtomId(2), AtomId(1)];
        let owner_attributes = AromaticSystemForm::from_electrons(vec![3, 1, 2]);
        let owner_to_local = DynPermutation::between(&owner_atoms, &local_atoms)
            .expect("the frames have the same atoms");
        let local_attributes = owner_attributes
            .clone()
            .reframe_by(&owner_to_local)
            .expect("the action has the form's degree");
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                aromatic: vec![(owner_atoms, owner_attributes)],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                id: AromaticSystemId(0),
                atoms: local_atoms,
                attributes: local_attributes,
            })]),
        );
        let action = reaction.representative_action();
        let inverse = action.inverse();
        let composite = action
            .compose(&inverse)
            .expect("an action and its inverse have the same domain");

        assert_eq!(
            reaction
                .clone()
                .reframe_by(&action)
                .and_then(|transported| transported.reframe_by(&inverse)),
            reaction.clone().reframe_by(&composite),
        );
        assert_eq!(reaction.clone().reframe_by(&composite), Some(reaction));
    }

    #[rstest]
    fn test_reaction_reframe_with_action_erased_entity() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            aromatic: vec![(
                vec![AtomId(1), AtomId(0)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            )],
            ..Default::default()
        });
        let attributes = MulticenterBondForm::from_electrons(vec![3, 4]);
        let id = MulticenterBondId(7);
        let reaction = Reaction::new(
            lhs.clone(),
            Deltas::from_iter([
                Delta::MulticenterBond(MulticenterBondDelta::Add {
                    id,
                    atoms: vec![AtomId(1), AtomId(0)],
                    attributes: attributes.clone(),
                }),
                Delta::MulticenterBond(MulticenterBondDelta::Remove {
                    id,
                    atoms: vec![AtomId(1), AtomId(0)],
                    attributes,
                }),
            ]),
        );

        let (reframed, action) = reaction
            .reframe_with_action()
            .expect("the created-then-removed entity is satisfiable");

        assert_eq!(
            reframed,
            Reaction::new(
                lhs.reframe().expect("the lhs is satisfiable"),
                Deltas::new()
            ),
        );
        assert_eq!(action.aromatic_systems().count(), 1);
        assert_eq!(
            action
                .multicenter_bonds()
                .action(id)
                .map(DynPermutation::image),
            Some([1, 0].as_slice()),
        );
    }

    #[rstest]
    fn test_reaction_reframe_change_chain() {
        let id = AromaticSystemId(0);
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                aromatic: vec![(
                    vec![AtomId(2), AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![30, 10, 20]),
                )],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                    id,
                    change: AromaticSystemFieldChange::Electrons {
                        old: ElectronCountsForm::Lit(vec![30, 10, 20]),
                        new: ElectronCountsForm::Lit(vec![31, 11, 21]),
                    },
                }),
                Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                    id,
                    change: AromaticSystemFieldChange::Electrons {
                        old: ElectronCountsForm::Lit(vec![31, 11, 21]),
                        new: ElectronCountsForm::Lit(vec![32, 12, 22]),
                    },
                }),
            ]),
        );
        let expected = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                aromatic: vec![(
                    vec![AtomId(0), AtomId(1), AtomId(2)],
                    AromaticSystemForm::from_electrons(vec![10, 20, 30]),
                )],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                id,
                change: AromaticSystemFieldChange::Electrons {
                    old: ElectronCountsForm::Lit(vec![10, 20, 30]),
                    new: ElectronCountsForm::Lit(vec![12, 22, 32]),
                },
            })]),
        );

        assert_eq!(reaction.reframe(), Ok(expected));
    }

    #[rstest]
    fn test_reaction_reframe_constraint_delta() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            noncovalent: vec![(
                AtomId(1),
                AtomId(0),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        let input = Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
            ],
        });
        let expected = Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
            ],
        });
        let reaction = Reaction::new(
            lhs.clone(),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(input))]),
        );

        assert_eq!(
            reaction.reframe(),
            Ok(Reaction::new(
                lhs.reframe().expect("the lhs is satisfiable"),
                Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(expected))]),
            )),
        );
    }

    #[rstest]
    #[case::atom(EntityKind::StereoAtom)]
    #[case::bond(EntityKind::StereoBond)]
    fn test_reaction_reframe_stereo_constraint(#[case] entity_kind: EntityKind) {
        let atom_ligands = [2, 1, 3, 4]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let bond_ligands = [1, 0, 3, 2]
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(2)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
                (AtomId(5), AtomId(1), BondForm::from_order(1)),
                (AtomId(6), AtomId(2), BondForm::from_order(1)),
                (AtomId(6), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                atom_ligands,
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            stereo_bonds: vec![(
                BondId(4),
                bond_ligands,
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let pair = StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2));
        let delta = match entity_kind {
            EntityKind::StereoAtom => Delta::StereoAtom(StereoAtomDelta::ModifyConstraint {
                id: StereoAtomId(0),
                kind: Some(StereoKind::Tetrahedral),
                old: None,
                new: Some(StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair,
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                })),
            }),
            EntityKind::StereoBond => Delta::StereoBond(StereoBondDelta::ModifyConstraint {
                id: StereoBondId(0),
                kind: Some(StereoKind::CisTrans),
                old: None,
                new: Some(StereoBondConstraintForm::Topicity(TopicityForm {
                    pair,
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                })),
            }),
            _ => unreachable!("test cases contain only stereo entity kinds"),
        };
        let reaction = Reaction::new(lhs, Deltas::from_iter([delta]));
        // This exact domain exercises sparse action discovery for frame-relative stereo
        // constraints; the comprehensive reaction strategy does not generate these delta arms.
        let expected = reaction
            .clone()
            .reframe_with_action()
            .map(|(reframed, _)| reframed);

        assert_eq!(reaction.reframe(), expected);
    }

    #[rstest]
    fn test_reaction_representative_action_contradiction() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            aromatic: vec![(
                vec![AtomId(1), AtomId(0)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            )],
            ..Default::default()
        });
        let reaction = Reaction::new(
            lhs,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Undetermined,
                        new: NumForm::Lit(1),
                    },
                }),
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(2),
                        new: NumForm::Lit(3),
                    },
                }),
            ]),
        );

        let action = reaction.representative_action();

        assert_eq!(
            action
                .aromatic_systems()
                .action(AromaticSystemId(0))
                .map(DynPermutation::image),
            Some([1, 0].as_slice()),
        );
        assert_eq!(reaction.normalize(), Err(Contradiction));
    }

    #[rstest]
    fn test_reaction_reframe_by_action_coverage() {
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                dative: vec![(
                    vec![AtomId(1), AtomId(0)],
                    AtomId(2),
                    DativeBondForm::from_order(1),
                )],
                ..Default::default()
            }),
            Deltas::new(),
        );
        let action = Reaction::default().representative_action();

        assert_eq!(reaction.reframe_by(&action), None);
    }

    #[rstest]
    fn test_reaction_reframe_by_action_degree() {
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                dative: vec![(
                    vec![AtomId(1), AtomId(0)],
                    AtomId(2),
                    DativeBondForm::from_order(1),
                )],
                ..Default::default()
            }),
            Deltas::new(),
        );
        let action = OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![DynPermutation::identity(1)])
                .expect("a dynamic permutation is a dative-bond action"),
            AromaticSystemsFrameAction::from_vec(vec![])
                .expect("the empty aromatic-system action is admissible"),
            MulticenterBondsFrameAction::from_vec(vec![])
                .expect("the empty multicenter-bond action is admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![])
                .expect("the empty noncovalent-bond action is admissible"),
            StereoAtomsFrameAction::from_vec(vec![])
                .expect("the empty stereo-atom action is admissible"),
            StereoBondsFrameAction::from_vec(vec![])
                .expect("the empty stereo-bond action is admissible"),
        );

        assert_eq!(reaction.reframe_by(&action), None);
    }

    #[rstest]
    #[case::invalid_reference(
        vec![AtomId(1), AtomId(7)],
        ReactionIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(7)) },
    )]
    #[case::incidence(
        vec![AtomId(0), AtomId(3)],
        ReactionIntegrityError::IncidenceMismatch {
            entity: Entity::DativeBond(DativeBondId(0)),
        },
    )]
    fn test_reaction_try_new_removal_error_precedence(
        #[case] donors: Vec<AtomId>,
        #[case] expected: ReactionIntegrityError,
    ) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            dative: vec![(
                vec![AtomId(0), AtomId(1)],
                AtomId(2),
                DativeBondForm::default(),
            )],
            ..Default::default()
        });
        let deltas = Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
            id: DativeBondId(0),
            donors,
            acceptor: AtomId(2),
            attributes: DativeBondForm::default(),
        })]);

        assert_eq!(Reaction::try_new(lhs, deltas), Err(expected));
    }

    #[rstest]
    fn test_reaction_try_new_stereo_bond_cross_block_error() {
        let ligands = (2..=5)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                ligands,
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let deltas = Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: [2, 4, 3, 5]
                .into_iter()
                .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
        })]);

        assert_eq!(
            Reaction::try_new(lhs, deltas),
            Err(ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            }),
        );
    }

    #[rstest]
    fn test_reaction_lhs() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = Reaction::new(lhs.clone(), Deltas::new());

        assert_eq!(reaction.lhs(), &lhs);
    }

    #[rstest]
    fn test_reaction_deltas() {
        let deltas = Deltas::from_iter([Delta::Atom(AtomDelta::Add {
            id: AtomId(0),
            attributes: AtomForm::from_element(Element::C),
        })]);
        let reaction = Reaction::new(Molecule::default(), deltas.clone());

        assert_eq!(reaction.deltas(), &deltas);
    }

    #[rstest]
    fn test_reaction_into_parts() {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let deltas = Deltas::new();
        let reaction = Reaction::new(lhs.clone(), deltas.clone());

        assert_eq!(reaction.into_parts(), (lhs, deltas));
    }

    #[rstest]
    #[case::unavailable_atom(
        Reaction {
            lhs: Molecule::default(),
            deltas: Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::default(),
            })]),
        },
        ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(0)),
        },
    )]
    fn test_reaction_try_new_reference_error(
        #[case] reaction: Reaction,
        #[case] expected: ReactionIntegrityError,
    ) {
        assert_eq!(
            Reaction::try_new(reaction.lhs, reaction.deltas),
            Err(expected)
        );
    }

    /// A stereo entity keeps its kind across a configuration change. An undetermined side asserts
    /// no geometry and restricts nothing, and the same-kind change is an ordinary modification.
    #[rustfmt::skip]
    #[rstest]
    #[case::kind_change(
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::kinded(StereoKind::Axial, StereoCoset::Lit(0)),
        Err(ReactionIntegrityError::StereoKindModified {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            old: StereoKind::Tetrahedral,
            new: StereoKind::Axial,
        }),
    )]
    #[case::same_kind(
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        Ok(()),
    )]
    #[case::old_undetermined(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::kinded(StereoKind::Axial, StereoCoset::Lit(0)),
        Ok(()),
    )]
    #[case::new_undetermined(
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Undetermined,
        Ok(()),
    )]
    fn test_reaction_try_new_stereo_kind(
        #[case] old: StereoConfigurationForm,
        #[case] new: StereoConfigurationForm,
        #[case] expected: Result<(), ReactionIntegrityError>,
    ) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        });
        let deltas = Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration { old, new },
            })]);

        assert_eq!(Reaction::try_new(lhs, deltas).map(drop), expected);
    }

    #[rstest]
    fn test_reaction_from_sides() {
        // C-C (order 1) → C-C (order 2) under the total atom correspondence: one bond-order modify.
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let right = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 2)
            .expect("correspondence producer preserves partial-bijection invariants");
        assert_eq!(
            Reaction::from_sides(left.clone(), right, atoms),
            Some(Reaction::new(
                left,
                Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                })]),
            )),
        );
    }

    #[rstest]
    #[case::bond_order(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }),
    )]
    #[case::overlay_removed(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::O) }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O)], bonds: vec![], ..Default::default() }),
    )]
    fn test_reaction_apply_at(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] atom_map: Vec<AtomId>,
        #[case] expected: Molecule,
    ) {
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&atom_map, host.atoms().count()),
        )
        .expect("the atom correspondence describes the molecule pair");
        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected
        );
    }

    #[rstest]
    fn test_reaction_apply_at_dative_bond_frame() {
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                dative: vec![(
                    vec![AtomId(0), AtomId(1)],
                    AtomId(2),
                    DativeBondForm::default(),
                )],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                id: DativeBondId(0),
                donors: vec![AtomId(0), AtomId(1)],
                acceptor: AtomId(2),
                attributes: DativeBondForm::default(),
            })]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            dative: vec![(
                vec![AtomId(1), AtomId(0)],
                AtomId(2),
                DativeBondForm::from_order(2),
            )],
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &host,
            Correspondence::from_images(&[AtomId(0), AtomId(1), AtomId(2)], 3),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected,
        );
    }

    #[rstest]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    fn test_reaction_apply_at_electron_frame(#[case] entity_kind: EntityKind) {
        let atoms = vec![AtomForm::from_element(Element::C); 3];
        let rule_frame = vec![AtomId(0), AtomId(1), AtomId(2)];
        let host_frame = vec![AtomId(2), AtomId(0), AtomId(1)];
        let (reaction, host, expected) = match entity_kind {
            EntityKind::AromaticSystem => (
                Reaction::new(
                    Molecule::from_entries(MoleculeEntries {
                        atoms: atoms.clone(),
                        aromatic: vec![(rule_frame.clone(), AromaticSystemForm::default())],
                        ..Default::default()
                    }),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                        id: AromaticSystemId(0),
                        change: AromaticSystemFieldChange::Electrons {
                            old: ElectronCountsForm::Undetermined,
                            new: ElectronCountsForm::Lit(vec![2, 3, 5]),
                        },
                    })]),
                ),
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    aromatic: vec![(
                        host_frame.clone(),
                        AromaticSystemForm {
                            electrons: ElectronCountsForm::Lit(vec![7, 11, 13]),
                            ..Default::default()
                        },
                    )],
                    ..Default::default()
                }),
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    aromatic: vec![(
                        host_frame.clone(),
                        AromaticSystemForm {
                            electrons: ElectronCountsForm::Lit(vec![5, 2, 3]),
                            ..Default::default()
                        },
                    )],
                    ..Default::default()
                }),
            ),
            EntityKind::MulticenterBond => (
                Reaction::new(
                    Molecule::from_entries(MoleculeEntries {
                        atoms: atoms.clone(),
                        multicenter: vec![(rule_frame, MulticenterBondForm::default())],
                        ..Default::default()
                    }),
                    Deltas::from_iter([Delta::MulticenterBond(
                        MulticenterBondDelta::ModifyField {
                            id: MulticenterBondId(0),
                            change: MulticenterBondFieldChange::Electrons {
                                old: ElectronCountsForm::Undetermined,
                                new: ElectronCountsForm::Lit(vec![2, 3, 5]),
                            },
                        },
                    )]),
                ),
                Molecule::from_entries(MoleculeEntries {
                    atoms: atoms.clone(),
                    multicenter: vec![(
                        host_frame.clone(),
                        MulticenterBondForm {
                            electrons: ElectronCountsForm::Lit(vec![7, 11, 13]),
                            ..Default::default()
                        },
                    )],
                    ..Default::default()
                }),
                Molecule::from_entries(MoleculeEntries {
                    atoms,
                    multicenter: vec![(
                        host_frame,
                        MulticenterBondForm {
                            electrons: ElectronCountsForm::Lit(vec![5, 2, 3]),
                            ..Default::default()
                        },
                    )],
                    ..Default::default()
                }),
            ),
            _ => unreachable!("test cases contain only electron-bearing entity kinds"),
        };
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &host,
            Correspondence::from_images(&[AtomId(0), AtomId(1), AtomId(2)], 3),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected,
        );
    }

    #[rstest]
    #[case::addition(false)]
    #[case::removal(true)]
    fn test_reaction_apply_at_constraint_frame(#[case] removal: bool) {
        let rule_constraint =
            Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(1))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(2))),
                ],
            });
        let host_constraint =
            Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(1))),
                ],
            });
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            constraints: if removal {
                Constraints::from(rule_constraint.clone())
            } else {
                Constraints::default()
            },
            ..Default::default()
        });
        let reaction = Reaction::new(
            lhs,
            Deltas::from_iter([Delta::Constraint(if removal {
                ConstraintDelta::Remove(rule_constraint)
            } else {
                ConstraintDelta::Add(rule_constraint)
            })]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            noncovalent: vec![(
                AtomId(1),
                AtomId(0),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            constraints: if removal {
                Constraints::from(host_constraint.clone())
            } else {
                Constraints::default()
            },
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            noncovalent: vec![(
                AtomId(1),
                AtomId(0),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            constraints: if removal {
                Constraints::default()
            } else {
                Constraints::from(host_constraint)
            },
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &host,
            Correspondence::from_images(&[AtomId(0), AtomId(1)], 2),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected,
        );
    }

    #[rstest]
    fn test_reaction_apply_at_stereo_constraint_frame() {
        let atom_rule_frame = [1, 2, 3, 4]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let atom_host_frame = [2, 1, 3, 4]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let bond_rule_frame = [0, 1, 2, 3]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let bond_host_frame = [2, 3, 0, 1]
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .to_vec();
        let atom_rule_constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        });
        let atom_host_constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        });
        let bond_rule_constraint = StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        });
        let bond_host_constraint = StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(2), StereoLigandPosition(3)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        });
        let atoms = vec![AtomForm::from_element(Element::C); 7];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(2)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
            (AtomId(5), AtomId(1), BondForm::from_order(1)),
            (AtomId(6), AtomId(2), BondForm::from_order(1)),
            (AtomId(6), AtomId(3), BondForm::from_order(1)),
        ];
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                atom_rule_frame,
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )],
            stereo_bonds: vec![(
                BondId(4),
                bond_rule_frame,
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )],
            ..Default::default()
        });
        let reaction = Reaction::new(
            lhs,
            Deltas::from_iter([
                Delta::StereoAtom(StereoAtomDelta::ModifyConstraint {
                    id: StereoAtomId(0),
                    kind: Some(StereoKind::Tetrahedral),
                    old: None,
                    new: Some(atom_rule_constraint.clone()),
                }),
                Delta::StereoBond(StereoBondDelta::ModifyConstraint {
                    id: StereoBondId(0),
                    kind: Some(StereoKind::CisTrans),
                    old: None,
                    new: Some(bond_rule_constraint.clone()),
                }),
                Delta::Constraint(ConstraintDelta::Add(Constraint::And(vec![
                    Constraint::StereoAtom(
                        StereoAtomId(0),
                        StereoKind::Tetrahedral,
                        atom_rule_constraint,
                    ),
                    Constraint::StereoBond(
                        StereoBondId(0),
                        StereoKind::CisTrans,
                        bond_rule_constraint,
                    ),
                ]))),
            ]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                atom_host_frame.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            stereo_bonds: vec![(
                BondId(4),
                bond_host_frame.clone(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                atom_host_frame,
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(atom_host_constraint.clone()),
            )],
            stereo_bonds: vec![(
                BondId(4),
                bond_host_frame,
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))
                    .with_constraint(bond_host_constraint.clone()),
            )],
            constraints: Constraints::from(Constraint::And(vec![
                Constraint::StereoAtom(
                    StereoAtomId(0),
                    StereoKind::Tetrahedral,
                    atom_host_constraint,
                ),
                Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, bond_host_constraint),
            ])),
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &host,
            Correspondence::from_images(
                &[
                    AtomId(0),
                    AtomId(1),
                    AtomId(2),
                    AtomId(3),
                    AtomId(4),
                    AtomId(5),
                    AtomId(6),
                ],
                7,
            ),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap().rhs(),
            &expected,
        );
    }

    #[rstest]
    fn test_reaction_apply_at_created_frame() {
        let constraint = Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(2))),
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(1))),
            ],
        });
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            ..Default::default()
        });
        let reaction = Reaction::new(
            lhs.clone(),
            Deltas::from_iter([
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(1), AtomId(0)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
                Delta::Constraint(ConstraintDelta::Add(constraint.clone())),
            ]),
        );
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            noncovalent: vec![(
                AtomId(1),
                AtomId(0),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            constraints: Constraints::from(constraint),
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &lhs,
            Correspondence::from_images(&[AtomId(0), AtomId(1)], 2),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&lhs, &correspondence).unwrap().rhs(),
            &expected,
        );
    }

    #[rstest]
    fn test_reaction_apply_at_frame_error() {
        let atoms = [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element)
            .collect::<Vec<_>>();
        let bonds = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect::<Vec<_>>();
        let ligands = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )],
            ..Default::default()
        });
        let reaction = Reaction::new(
            lhs,
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
                    new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
                },
            })]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                ligands,
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            )],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            &host,
            Correspondence::from_images(
                &[AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
                5,
            ),
        )
        .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            reaction.apply_at(&host, &correspondence),
            Err(ApplyError::StereoFrameMismatch {
                entity: Entity::StereoAtom(StereoAtomId(0)),
            }),
        );
    }

    // `dangling_*`: the rule deletes a host atom still carrying an undeleted bond/overlay (DPO gluing
    // condition). `structural_conflict`: the rule's add lands a second bond on an already-bonded atom
    // pair, so checked product publication rejects the parallel bonds.
    #[rstest]
    #[case::dangling_bond(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        vec![AtomId(0)],
        ApplyError::Dangling { host_atom: AtomId(0) },
    )]
    #[case::dangling_noncovalent(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O)], bonds: vec![], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::O),
            })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        vec![AtomId(0)],
        ApplyError::Dangling { host_atom: AtomId(0) },
    )]
    #[case::structural_conflict(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(1)],
                attributes: BondForm::from_order(1),
            })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        vec![AtomId(0), AtomId(1)],
        ApplyError::StructuralConflict,
    )]
    fn test_reaction_apply_at_error(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] atom_map: Vec<AtomId>,
        #[case] expected: ApplyError,
    ) {
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&atom_map, host.atoms().count()),
        )
        .expect("the atom correspondence describes the molecule pair");
        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap_err(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::missing_atom(
        Reaction::new(Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }), Deltas::new()),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
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
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::new(),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        MoleculeCorrespondence::new(
            Correspondence::from_images(&[AtomId(0), AtomId(1)], 2),
            Correspondence::new(vec![], 1, 1).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        ),
        ApplyError::CorrespondenceMismatch { entity: Entity::Bond(BondId(0)) },
    )]
    #[case::noncovalent_incidence(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O); 3], noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
            Deltas::new(),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::O); 3], noncovalent: vec![(AtomId(0), AtomId(2), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))], constraints: Constraints::new(), ..Default::default() }),
        MoleculeCorrespondence::new(
            Correspondence::from_images(&[AtomId(0), AtomId(1), AtomId(2)], 3),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1).expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), Correspondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        ),
        ApplyError::CorrespondenceMismatch { entity: Entity::NoncovalentBond(NoncovalentBondId(0)) },
    )]
    fn test_reaction_apply_at_correspondence_error(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] correspondence: MoleculeCorrespondence,
        #[case] expected: ApplyError,
    ) {
        assert_eq!(reaction.apply_at(&host, &correspondence), Err(expected));
    }

    #[rstest]
    #[case::field(Delta::StereoAtom(StereoAtomDelta::ModifyField {
        id: StereoAtomId(0),
        change: StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
            new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
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
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
    }))]
    fn test_reaction_apply_at_stereo_atom_error(#[case] delta: Delta) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: (1..=4)
                .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        });
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: [1, 2, 3, 5]
                .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
                .into(),
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new((0..6u32).map(|id| (AtomId(id), AtomId(id))).collect(), 6, 6)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new((0..3u32).map(|id| (BondId(id), BondId(id))).collect(), 4, 4)
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
        let reaction = Reaction::new(lhs, Deltas::from_iter([delta]));

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
            old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0u32),
            new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1u32),
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
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
    }))]
    #[case::constraint(Delta::Constraint(ConstraintDelta::Add(Constraint::StereoBond(
        StereoBondId(0),
        StereoKind::CisTrans,
        StereoBondConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(
                StereoLigandPosition(0),
                StereoLigandPosition(2),
            ),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }),
    ))))]
    fn test_reaction_apply_at_stereo_bond_error(#[case] delta: Delta) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(6), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            )],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new((0..7u32).map(|id| (AtomId(id), AtomId(id))).collect(), 7, 7)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new((0..4u32).map(|id| (BondId(id), BondId(id))).collect(), 5, 5)
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
        let reaction = Reaction::new(lhs, Deltas::from_iter([delta]));

        assert_eq!(
            reaction.apply_at(&host, &correspondence).unwrap_err(),
            ApplyError::StereoFrameMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            },
        );
    }

    #[rstest]
    fn test_reaction_apply_at_molecule_constraint() {
        // A reaction adding a molecule-level `ChargeSum` over its lhs atoms; applied at a match
        // that maps lhs atoms 0,1 → host atoms 1,2, the constraint's refs re-anchor to the host.
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                    sum: NumForm::Lit(0),
                }),
            ))]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&[AtomId(1), AtomId(2)], host.atoms().count()),
        )
        .expect("the atom correspondence describes the molecule pair");
        let result = reaction.apply_at(&host, &correspondence).unwrap();
        assert_eq!(
            result.rhs().constraints(),
            &Constraints::from(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
                sum: NumForm::Lit(0),
            })),
        );
    }

    #[rstest]
    fn test_reaction_apply_at_molecule_constraint_created() {
        let constraint = Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(0), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(0),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(0),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0),
                parallel: BondId(0),
            }),
        ]);
        let reaction = Reaction::new(
            Molecule::default(),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(0),
                    attributes: AtomForm::from_element(Element::C),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::N),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    attributes: AtomForm::from_element(Element::H),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(3),
                    attributes: AtomForm::from_element(Element::H),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(2),
                    atoms: [AtomId(0), AtomId(3)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::DativeBond(DativeBondDelta::Add {
                    id: DativeBondId(0),
                    donors: vec![AtomId(0)],
                    acceptor: AtomId(1),
                    attributes: DativeBondForm::from_order(1),
                }),
                Delta::AromaticSystem(AromaticSystemDelta::Add {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    attributes: AromaticSystemForm::default(),
                }),
                Delta::MulticenterBond(MulticenterBondDelta::Add {
                    id: MulticenterBondId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    attributes: MulticenterBondForm::default(),
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
                Delta::StereoAtom(StereoAtomDelta::Add {
                    id: StereoAtomId(0),
                    site: AtomId(0),
                    ligands: vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    ],
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                }),
                Delta::StereoBond(StereoBondDelta::Add {
                    id: StereoBondId(0),
                    site: BondId(0),
                    ligands: vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                }),
                Delta::Constraint(ConstraintDelta::Add(constraint.clone())),
            ]),
        );

        let host = Molecule::default();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::new(Vec::new(), 0, 0)
                .expect("correspondence producer preserves partial-bijection invariants"),
        )
        .expect("the atom correspondence describes the molecule pair");
        let result = reaction.apply_at(&host, &correspondence).unwrap();

        assert_eq!(result.rhs().constraints(), &Constraints::from(constraint));
    }

    #[rstest]
    #[case::valid(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Deltas::new(),
        ),
    )]
    #[case::normalized_add_remove_cancellation(
        Reaction::new(
            Molecule::default(),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add { id: AtomId(0), attributes: AtomForm::from_element(Element::C) }),
                Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) }),
            ]),
        ),
    )]
    #[case::unordered_bond_incidence(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(1), AtomId(0)],
                attributes: BondForm::from_order(1),
            })]),
        ),
    )]
    fn test_reaction_check_preconditions(#[case] reaction: Reaction) {
        assert_eq!(reaction.check_preconditions(), Ok(()));
    }

    #[rstest]
    #[case::inconsistent_reaction(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                }),
                Delta::Atom(AtomDelta::ModifyField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(2),
                        new: NumForm::Lit(3),
                    },
                }),
            ]),
        ),
        ApplyPreconditionError::InconsistentReaction,
    )]
    #[case::reaction_dpo(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) })]),
        ),
        ApplyPreconditionError::ReactionDpo(DpoContradiction::DanglingBond { atom: AtomId(0), bond: BondId(0) }),
    )]
    fn test_reaction_check_preconditions_error(
        #[case] reaction: Reaction,
        #[case] expected: ApplyPreconditionError,
    ) {
        assert_eq!(reaction.check_preconditions().unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::bond_order(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) } })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() }),
        vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() })],
    )]
    #[case::match_rejection(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))], ..Default::default() }),
        vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))], ..Default::default() })],
    )]
    #[case::host_relative_update(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
            })]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)], ..Default::default() }),
        vec![Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C).with_charge(1_i64)], ..Default::default() })],
    )]
    fn test_reaction_apply(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] expected: Vec<Molecule>,
    ) {
        let products: Vec<Molecule> = reaction
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
    fn test_reaction_apply_match_algorithm(#[case] match_algorithm: SubstructureMatchAlgorithm) {
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        let reaction = Reaction::new(host.clone(), Deltas::new());
        let products: Vec<Molecule> = reaction
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
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    #[case::transaction(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0)]),
                    sum: NumForm::Lit(0),
                })),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: Some(vec![AtomId(0)]),
                    sum: NumForm::Lit(0),
                }),
            ))]),
        ),
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C)], ..Default::default() }),
        ApplyError::Transaction(TransactionError::MissingEntry),
    )]
    fn test_reaction_apply_error(
        #[case] reaction: Reaction,
        #[case] host: Molecule,
        #[case] expected: ApplyError,
    ) {
        let mut applications = reaction.apply(&host, MATCH_CONFIG).unwrap();

        assert_eq!(applications.next().unwrap().unwrap_err(), expected);
        assert_eq!(applications.next(), None);
    }

    #[fixture]
    fn tetrahedral_inversion() -> Reaction {
        // Invert a tetrahedral C(0) whose ligands F,Cl,Br,I are stated in ascending order: coset 0 → 1.
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                    AtomForm::from_element(Element::I),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::Kinded(
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
    fn test_reaction_apply_stereo_cross_frame(
        tetrahedral_inversion: Reaction,
        #[case] host_ligands: [u32; 4],
        #[case] host_coset: u32,
        #[case] product_coset: u32,
    ) {
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                host_ligands
                    .iter()
                    .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, host_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                host_ligands
                    .iter()
                    .map(|&x| StereoLigand::new(AtomId(x), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, product_coset),
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

    #[fixture]
    fn explicit_hydrogen_removal() -> Reaction {
        let ligands = vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::H),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(AtomId(0), ligands.clone(), attributes.clone())],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands,
                attributes,
            })]),
        )
    }

    #[rstest]
    #[case::same_frame([0, 1, 2, 3], 0)]
    #[case::reordered_frame([1, 0, 2, 3], 1)]
    fn test_reaction_apply_stereo_atom_removal_distinct_ligands(
        explicit_hydrogen_removal: Reaction,
        #[case] host_order: [usize; 4],
        #[case] host_coset: u32,
    ) {
        let atoms = vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ];
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                host_order
                    .into_iter()
                    .map(|position| {
                        [
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                        ][position]
                    })
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, host_coset),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            constraints: Constraints::new(),
            ..Default::default()
        });
        let rhs = explicit_hydrogen_removal
            .apply(&host, MATCH_CONFIG)
            .unwrap()
            .next()
            .expect("the removal rule matches the host")
            .unwrap()
            .rhs()
            .clone();
        assert_eq!(rhs, expected);
    }

    #[fixture]
    fn square_planar_modification() -> Reaction {
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::Pt),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::Br),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::SquarePlanar,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::Kinded(
                        StereoKind::SquarePlanar,
                        StereoCoset::Lit(2),
                    ),
                },
            })]),
        )
    }

    #[rstest]
    #[case::matching_configuration(0, Some(2))]
    #[case::different_configuration(1, None)]
    #[case::product_configuration(2, None)]
    fn test_reaction_apply_stereo_atom_modification_distinct_ligands(
        square_planar_modification: Reaction,
        #[case] host_coset: u32,
        #[case] product_coset: Option<u32>,
    ) {
        let molecule = |coset: u32| {
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::Pt),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::Br),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(4), BondForm::from_order(1)),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::SquarePlanar, coset),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            })
        };
        let mut applications = square_planar_modification
            .apply(&molecule(host_coset), MATCH_CONFIG)
            .unwrap();
        match product_coset {
            Some(coset) => {
                let rhs = applications
                    .next()
                    .expect("the modification rule matches the host")
                    .unwrap()
                    .rhs()
                    .clone();
                assert_eq!(rhs, molecule(coset));
            }
            None => assert!(applications.next().is_none()),
        }
    }

    // A stereo-bond addition can name a bond added by the same reaction. Its `BondHandle::New`
    // index is the bond creation ordinal and is independent of the atom creation namespace.
    #[rstest]
    #[case::coset_0(0u32)]
    fn test_reaction_apply_stereo_bond_created_site(#[case] coset: u32) {
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                bonds: vec![],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(2),
                }),
                Delta::StereoBond(StereoBondDelta::Add {
                    id: StereoBondId(0),
                    site: BondId(0),
                    ligands: vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                }),
            ]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, coset),
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
    fn test_reaction_apply_two_stereo_centers(#[case] coset: StereoCoset) {
        let center = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, coset.clone()),
                ),
                (
                    AtomId(1),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, coset.clone()),
                ),
            ],
            constraints: Constraints::new(),
            ..Default::default()
        });
        let rhs = Reaction::new(center.clone(), Deltas::new())
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
    fn test_reaction_apply_at_comap() {
        // Remove atom O (id 1) and its bond: host C-O ⇒ product C. Atom 0 is preserved (matched),
        // atom 1 is deleted (left-unmatched), so the comap's atom map records exactly that.
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    attributes: BondForm::from_order(1),
                }),
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    attributes: AtomForm::from_element(Element::O),
                }),
            ]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&[AtomId(0), AtomId(1)], host.atoms().count()),
        )
        .expect("the atom correspondence describes the molecule pair");
        let derivation = reaction.apply_at(&host, &correspondence).unwrap();
        assert_eq!(
            derivation.rhs(),
            &Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                bonds: vec![],
                ..Default::default()
            })
        );
        assert_eq!(
            derivation.atom_correspondence().matched_pairs(),
            &[(AtomId(0), AtomId(0))]
        );
        assert_eq!(
            derivation.atom_correspondence().left_unmatched(),
            vec![AtomId(1)]
        );
    }
}
