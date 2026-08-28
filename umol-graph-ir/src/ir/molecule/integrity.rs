//! Representation-integrity checks for [`Molecule`].

use std::collections::{BTreeSet, HashSet};
use std::iter;

use thiserror::Error;
use umol_graph_core::NodeId;
use umol_perm::{Permutation, MAX_DEGREE};

use super::super::constraint::{
    Constraint, StereoAtomConstraintForm, StereoBondConstraintForm, StereoLigandPair,
};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::id::{AtomId, BondId};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::stereo::{StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm};
use super::Molecule;

/// Failure of the representation contract required to interpret a [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeIntegrityError {
    #[error("molecule references unavailable {entity}")]
    InvalidReference { entity: Entity },
    #[error("molecule stores {actual} {kind:?} values for {expected} graph entities")]
    EntityCountMismatch {
        kind: super::super::entity::EntityKind,
        expected: usize,
        actual: usize,
    },
    #[error(
        "{entity}: electron-count vector has length {electron_counts}, expected {participants}"
    )]
    ElectronCountLengthMismatch {
        entity: Entity,
        participants: usize,
        electron_counts: usize,
    },
    #[error("{entity}: participant atom {atom:?} is duplicated")]
    DuplicateParticipant { entity: Entity, atom: AtomId },
    #[error("bond: parallel bonds on atoms {atoms:?}")]
    BondsParallel { atoms: [AtomId; 2] },
    #[error(
        "dative bond: parallel datives to acceptor {acceptor:?} sharing donor {shared_donor:?}"
    )]
    DativeBondsParallel {
        acceptor: AtomId,
        shared_donor: AtomId,
    },
    #[error("noncovalent bond: parallel bonds on atoms {atoms:?}")]
    NoncovalentBondsParallel { atoms: [AtomId; 2] },
    #[error("aromatic systems: overlap on atom {atom:?}")]
    AromaticSystemsOverlap { atom: AtomId },
    #[error("multicenter bonds: identical participant set {atoms:?}")]
    MulticenterBondsIdentical { atoms: Vec<AtomId> },
    #[error("stereo atom: duplicate site {atom:?}")]
    StereoAtomSitesDuplicate { atom: AtomId },
    #[error("stereo bond: duplicate site {bond:?}")]
    StereoBondSitesDuplicate { bond: BondId },
    #[error("{entity}: stereo ligand {ligand:?} is duplicated in the frame")]
    DuplicateStereoLigand {
        entity: Entity,
        ligand: StereoLigand,
    },
    #[error(
        "{entity}: stereo frame has degree {degree}, exceeding the supported maximum {maximum}"
    )]
    StereoFrameDegreeTooLarge {
        entity: Entity,
        degree: usize,
        maximum: usize,
    },
    #[error("{entity}: ligand frame does not match stereo-site incidence")]
    StereoLigandIncidenceMismatch { entity: Entity },
    #[error("{entity}: stereo kind {kind:?} is not admissible for this site type")]
    StereoKindSiteMismatch { entity: Entity, kind: StereoKind },
    #[error("{entity}: stereo frame has {actual} ligands, expected {expected} for {kind:?}")]
    StereoLigandArity {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    #[error("{entity}: coset {coset} is outside 0..{count} for {kind:?}")]
    StereoCosetOutOfRange {
        entity: Entity,
        kind: StereoKind,
        coset: u32,
        count: usize,
    },
    #[error(
        "{entity}: permutation has degree {actual}, expected {expected} for the stored ligand frame"
    )]
    StereoPermutationDegree {
        entity: Entity,
        expected: usize,
        actual: usize,
    },
    #[error("{entity}: ligand position {position} is outside 0..{degree}")]
    StereoLigandPositionOutOfRange {
        entity: Entity,
        position: usize,
        degree: usize,
    },
}

impl Molecule {
    /// Check the representation invariants required to interpret this molecule.
    ///
    /// This checks stored references, parallel collection shapes, the fixed relation semantics of
    /// every entity family, and kind-dependent stereo domains. It does not check chemistry or
    /// constraint satisfaction.
    pub fn check_integrity(&self) -> Result<(), MoleculeIntegrityError> {
        use super::super::entity::EntityKind;

        if self.atoms.len() != self.graph.node_count() {
            return Err(MoleculeIntegrityError::EntityCountMismatch {
                kind: EntityKind::Atom,
                expected: self.graph.node_count(),
                actual: self.atoms.len(),
            });
        }
        if self.bonds.len() != self.graph.edge_count() {
            return Err(MoleculeIntegrityError::EntityCountMismatch {
                kind: EntityKind::Bond,
                expected: self.graph.edge_count(),
                actual: self.bonds.len(),
            });
        }

        let contains = |entity| match entity {
            Entity::Atom(id) => id.index() < self.atoms.len(),
            Entity::Bond(id) => id.index() < self.bonds.len(),
            Entity::DativeBond(id) => id.index() < self.dative_bonds.count(),
            Entity::AromaticSystem(id) => id.index() < self.aromatic_systems.count(),
            Entity::MulticenterBond(id) => id.index() < self.multicenter_bonds.count(),
            Entity::NoncovalentBond(id) => id.index() < self.noncovalent_bonds.count(),
            Entity::StereoAtom(id) => id.index() < self.stereo_atoms.count(),
            Entity::StereoBond(id) => id.index() < self.stereo_bonds.count(),
        };

        let mut bond_pairs = HashSet::new();
        for view in self.bonds().iter() {
            let entity = Entity::Bond(view.id);
            let atoms = view.atom_ids();
            require_references(&contains, atoms.into_iter().map(Entity::Atom))?;
            check_unique_pair(entity, atoms)?;
            let pair = unordered_pair(atoms);
            if !bond_pairs.insert(pair) {
                return Err(MoleculeIntegrityError::BondsParallel { atoms: pair });
            }
        }

        let mut dative_incidences = HashSet::new();
        for view in self.dative_bonds().iter() {
            let entity = Entity::DativeBond(view.id);
            let acceptor = view.acceptor_id();
            require_references(&contains, view.atom_ids().map(Entity::Atom))?;
            check_unique_participants(entity, view.atom_ids())?;
            for donor in view.donor_ids() {
                if !dative_incidences.insert((acceptor, donor)) {
                    return Err(MoleculeIntegrityError::DativeBondsParallel {
                        acceptor,
                        shared_donor: donor,
                    });
                }
            }
        }

        let mut aromatic_membership = HashSet::new();
        for view in self.aromatic_systems().iter() {
            let entity = Entity::AromaticSystem(view.id);
            require_references(&contains, view.atom_ids().map(Entity::Atom))?;
            check_unique_participants(entity, view.atom_ids())?;
            for atom in view.atom_ids() {
                if !aromatic_membership.insert(atom) {
                    return Err(MoleculeIntegrityError::AromaticSystemsOverlap { atom });
                }
            }
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }

        let mut multicenter_participant_sets = HashSet::new();
        for view in self.multicenter_bonds().iter() {
            let entity = Entity::MulticenterBond(view.id);
            require_references(&contains, view.atom_ids().map(Entity::Atom))?;
            check_unique_participants(entity, view.atom_ids())?;
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            if !multicenter_participant_sets.insert(atoms.iter().copied().collect::<BTreeSet<_>>())
            {
                return Err(MoleculeIntegrityError::MulticenterBondsIdentical { atoms });
            }
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }

        let mut noncovalent_pairs = HashSet::new();
        for view in self.noncovalent_bonds().iter() {
            let entity = Entity::NoncovalentBond(view.id);
            let atoms = view.atom_ids();
            require_references(&contains, atoms.into_iter().map(Entity::Atom))?;
            check_unique_pair(entity, atoms)?;
            let pair = unordered_pair(atoms);
            if !noncovalent_pairs.insert(pair) {
                return Err(MoleculeIntegrityError::NoncovalentBondsParallel { atoms: pair });
            }
        }

        let mut stereo_atom_sites = HashSet::new();
        for view in self.stereo_atoms().iter() {
            let entity = Entity::StereoAtom(view.id);
            let site = view.site_id();
            let ligand_frame = view.ligand_frame();
            require_reference(&contains, Entity::Atom(site))?;
            require_references(
                &contains,
                ligand_frame
                    .iter()
                    .map(|ligand| Entity::Atom(ligand.atom_id)),
            )?;
            check_stereo_atom_entry(entity, site, &ligand_frame, view.attributes)?;
            if !stereo_atom_sites.insert(site) {
                return Err(MoleculeIntegrityError::StereoAtomSitesDuplicate { atom: site });
            }
            if ligand_frame
                .iter()
                .any(|&ligand| !ligand_matches_site(self, ligand, site, None))
            {
                return Err(MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity });
            }
        }

        let mut stereo_bond_sites = HashSet::new();
        for view in self.stereo_bonds().iter() {
            let entity = Entity::StereoBond(view.id);
            let site = view.site_id();
            let ligand_frame = view.ligand_frame();
            require_reference(&contains, Entity::Bond(site))?;
            require_references(
                &contains,
                ligand_frame
                    .iter()
                    .map(|ligand| Entity::Atom(ligand.atom_id)),
            )?;
            check_stereo_bond_entry(entity, &ligand_frame, view.attributes)?;
            if !stereo_bond_sites.insert(site) {
                return Err(MoleculeIntegrityError::StereoBondSitesDuplicate { bond: site });
            }
            let [first, second] = view.site().atom_ids();
            let matches_endpoint_order = |first, second| {
                ligand_frame.len() == 4
                    && ligand_frame[..2]
                        .iter()
                        .all(|&ligand| ligand_matches_site(self, ligand, first, Some(second)))
                    && ligand_frame[2..]
                        .iter()
                        .all(|&ligand| ligand_matches_site(self, ligand, second, Some(first)))
            };
            if !matches_endpoint_order(first, second) && !matches_endpoint_order(second, first) {
                return Err(MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity });
            }
        }
        for constraint in self.constraints.iter() {
            super::validate_constraint_references(constraint, &contains)
                .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })?;
            check_molecule_constraint(self, constraint)?;
        }
        Ok(())
    }
}

fn ligand_matches_site(
    molecule: &Molecule,
    ligand: StereoLigand,
    site: AtomId,
    opposite_site: Option<AtomId>,
) -> bool {
    match ligand.kind {
        StereoLigandKind::Atom => {
            Some(ligand.atom_id) != opposite_site
                && molecule
                    .graph
                    .find_edge(NodeId::from(site), NodeId::from(ligand.atom_id))
                    .is_some()
        }
        StereoLigandKind::ImplicitHydrogen | StereoLigandKind::LonePair => ligand.atom_id == site,
    }
}

fn unordered_pair([first, second]: [AtomId; 2]) -> [AtomId; 2] {
    if first <= second {
        [first, second]
    } else {
        [second, first]
    }
}

fn check_unique_pair(
    entity: Entity,
    [first, second]: [AtomId; 2],
) -> Result<(), MoleculeIntegrityError> {
    if first == second {
        Err(MoleculeIntegrityError::DuplicateParticipant {
            entity,
            atom: first,
        })
    } else {
        Ok(())
    }
}

fn check_unique_participants(
    entity: Entity,
    participants: impl IntoIterator<Item = AtomId>,
) -> Result<(), MoleculeIntegrityError> {
    let mut seen = HashSet::new();
    for atom in participants {
        if !seen.insert(atom) {
            return Err(MoleculeIntegrityError::DuplicateParticipant { entity, atom });
        }
    }
    Ok(())
}

pub(crate) fn check_stereo_frame(
    entity: Entity,
    ligand_frame: &[StereoLigand],
) -> Result<(), MoleculeIntegrityError> {
    if ligand_frame.len() > MAX_DEGREE {
        return Err(MoleculeIntegrityError::StereoFrameDegreeTooLarge {
            entity,
            degree: ligand_frame.len(),
            maximum: MAX_DEGREE,
        });
    }

    let mut seen = HashSet::new();
    for &ligand in ligand_frame {
        if !seen.insert(ligand) {
            return Err(MoleculeIntegrityError::DuplicateStereoLigand { entity, ligand });
        }
    }
    Ok(())
}

pub(crate) fn check_stereo_atom_entry(
    entity: Entity,
    site: AtomId,
    ligand_frame: &[StereoLigand],
    attributes: &super::super::stereo::StereoAtomForm,
) -> Result<(), MoleculeIntegrityError> {
    check_stereo_frame(entity, ligand_frame)?;
    check_unique_participants(
        entity,
        iter::once(site).chain(
            ligand_frame
                .iter()
                .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
                .map(|ligand| ligand.atom_id),
        ),
    )?;
    check_stereo_atom(entity, ligand_frame.len(), attributes)
}

pub(crate) fn check_stereo_bond_entry(
    entity: Entity,
    ligand_frame: &[StereoLigand],
    attributes: &super::super::stereo::StereoBondForm,
) -> Result<(), MoleculeIntegrityError> {
    check_stereo_frame(entity, ligand_frame)?;
    check_unique_participants(
        entity,
        ligand_frame
            .iter()
            .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
            .map(|ligand| ligand.atom_id),
    )?;
    check_stereo_bond(entity, ligand_frame.len(), attributes)
}

pub(crate) fn check_stereo_atom_kind(
    entity: Entity,
    kind: StereoKind,
) -> Result<(), MoleculeIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Atom)
}

pub(crate) fn check_stereo_bond_kind(
    entity: Entity,
    kind: StereoKind,
) -> Result<(), MoleculeIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Bond)
}

pub(super) fn require_reference(
    contains: &dyn Fn(Entity) -> bool,
    entity: Entity,
) -> Result<(), MoleculeIntegrityError> {
    if contains(entity) {
        Ok(())
    } else {
        Err(MoleculeIntegrityError::InvalidReference { entity })
    }
}

pub(super) fn require_references(
    contains: &dyn Fn(Entity) -> bool,
    entities: impl IntoIterator<Item = Entity>,
) -> Result<(), MoleculeIntegrityError> {
    for entity in entities {
        require_reference(contains, entity)?;
    }
    Ok(())
}

fn check_electron_count_length(
    entity: Entity,
    participants: usize,
    electrons: &ElectronCountsForm,
) -> Result<(), MoleculeIntegrityError> {
    if let ElectronCountsForm::Lit(counts) = electrons {
        if counts.len() != participants {
            return Err(MoleculeIntegrityError::ElectronCountLengthMismatch {
                entity,
                participants,
                electron_counts: counts.len(),
            });
        }
    }
    Ok(())
}

fn check_stereo_atom(
    entity: Entity,
    ligand_count: usize,
    attributes: &super::super::stereo::StereoAtomForm,
) -> Result<(), MoleculeIntegrityError> {
    check_configuration_site_kind(entity, &attributes.configuration, StereoSite::Atom)?;
    check_configuration(entity, ligand_count, &attributes.configuration)?;
    for constraint in attributes.constraints.iter() {
        check_stereo_atom_constraint(entity, ligand_count, constraint)?;
    }
    Ok(())
}

fn check_stereo_bond(
    entity: Entity,
    ligand_count: usize,
    attributes: &super::super::stereo::StereoBondForm,
) -> Result<(), MoleculeIntegrityError> {
    check_configuration_site_kind(entity, &attributes.configuration, StereoSite::Bond)?;
    check_configuration(entity, ligand_count, &attributes.configuration)?;
    for constraint in attributes.constraints.iter() {
        check_stereo_bond_constraint(entity, ligand_count, constraint)?;
    }
    Ok(())
}

fn check_molecule_constraint(
    molecule: &Molecule,
    constraint: &Constraint,
) -> Result<(), MoleculeIntegrityError> {
    match constraint {
        Constraint::StereoAtom(id, kind, constraint) => {
            let entity = Entity::StereoAtom(*id);
            let ligand_count = molecule.stereo_atom(*id).ligand_count();
            check_stereo_site_kind(entity, *kind, StereoSite::Atom)?;
            check_stereo_frame_arity(entity, ligand_count, *kind)?;
            check_stereo_atom_constraint(entity, ligand_count, constraint)
        }
        Constraint::StereoBond(id, kind, constraint) => {
            let entity = Entity::StereoBond(*id);
            let ligand_count = molecule.stereo_bond(*id).ligand_count();
            check_stereo_site_kind(entity, *kind, StereoSite::Bond)?;
            check_stereo_frame_arity(entity, ligand_count, *kind)?;
            check_stereo_bond_constraint(entity, ligand_count, constraint)
        }
        Constraint::And(constraints) | Constraint::Or(constraints) => {
            for constraint in constraints {
                check_molecule_constraint(molecule, constraint)?;
            }
            Ok(())
        }
        Constraint::Not(constraint) => check_molecule_constraint(molecule, constraint),
        _ => Ok(()),
    }
}

/// Which site a stereo entry sits on. Local to the admissibility check; the entity id already
/// carries the distinction everywhere else.
#[derive(Clone, Copy)]
enum StereoSite {
    Atom,
    Bond,
}

/// A stereo kind describes a coordination geometry, and a geometry belongs to an atom or to a bond.
/// Arity cannot separate them: `Tetrahedral`, `CisTrans`, `Axial`, and `SquarePlanar` all have
/// degree 4. `Axial` is admissible on both, since axial chirality arises at an allene's central
/// atom and about an atropisomeric biaryl bond.
///
/// Matched exhaustively so that a new stereo kind must decide its site here.
fn check_configuration_site_kind(
    entity: Entity,
    configuration: &StereoConfigurationForm,
    site: StereoSite,
) -> Result<(), MoleculeIntegrityError> {
    let Some(kind) = configuration.kind() else {
        return Ok(());
    };
    check_stereo_site_kind(entity, kind, site)
}

fn check_stereo_site_kind(
    entity: Entity,
    kind: StereoKind,
    site: StereoSite,
) -> Result<(), MoleculeIntegrityError> {
    let admissible = match (site, kind) {
        (
            StereoSite::Atom,
            StereoKind::Tetrahedral
            | StereoKind::SquarePlanar
            | StereoKind::TrigonalBipyramidal
            | StereoKind::Octahedral
            | StereoKind::Axial,
        ) => true,
        (StereoSite::Atom, StereoKind::CisTrans) => false,
        (StereoSite::Bond, StereoKind::CisTrans | StereoKind::Axial) => true,
        (
            StereoSite::Bond,
            StereoKind::Tetrahedral
            | StereoKind::SquarePlanar
            | StereoKind::TrigonalBipyramidal
            | StereoKind::Octahedral,
        ) => false,
    };
    if admissible {
        Ok(())
    } else {
        Err(MoleculeIntegrityError::StereoKindSiteMismatch { entity, kind })
    }
}

fn check_stereo_frame_arity(
    entity: Entity,
    ligand_count: usize,
    kind: StereoKind,
) -> Result<(), MoleculeIntegrityError> {
    if ligand_count != kind.degree() {
        return Err(MoleculeIntegrityError::StereoLigandArity {
            entity,
            kind,
            expected: kind.degree(),
            actual: ligand_count,
        });
    }
    Ok(())
}

fn check_stereo_atom_constraint(
    entity: Entity,
    ligand_count: usize,
    constraint: &StereoAtomConstraintForm,
) -> Result<(), MoleculeIntegrityError> {
    match constraint {
        StereoAtomConstraintForm::LigandSymmetry(value) => {
            check_permutation(entity, ligand_count, value.permutation.permutation.0)
        }
        StereoAtomConstraintForm::Fluxionality(value) => {
            check_permutation(entity, ligand_count, value.permutation.0)
        }
        StereoAtomConstraintForm::Topicity(value) => check_pair(entity, ligand_count, value.pair),
        StereoAtomConstraintForm::Stereogenicity(_) => Ok(()),
    }
}

fn check_stereo_bond_constraint(
    entity: Entity,
    ligand_count: usize,
    constraint: &StereoBondConstraintForm,
) -> Result<(), MoleculeIntegrityError> {
    match constraint {
        StereoBondConstraintForm::LigandSymmetry(value) => {
            check_permutation(entity, ligand_count, value.permutation.permutation.0)
        }
        StereoBondConstraintForm::Fluxionality(value) => {
            check_permutation(entity, ligand_count, value.permutation.0)
        }
        StereoBondConstraintForm::Topicity(value) => check_pair(entity, ligand_count, value.pair),
        StereoBondConstraintForm::Stereogenicity(_) => Ok(()),
    }
}

fn check_configuration(
    entity: Entity,
    ligand_count: usize,
    configuration: &StereoConfigurationForm,
) -> Result<(), MoleculeIntegrityError> {
    let StereoConfigurationForm::Kinded(kind, coset) = configuration else {
        return Ok(());
    };
    check_stereo_frame_arity(entity, ligand_count, *kind)?;
    check_coset(entity, *kind, coset)?;
    Ok(())
}

fn check_coset(
    entity: Entity,
    kind: StereoKind,
    coset: &StereoCoset,
) -> Result<(), MoleculeIntegrityError> {
    match coset {
        StereoCoset::Undetermined => Ok(()),
        StereoCoset::Lit(value) => check_coset_index(entity, kind, *value),
        StereoCoset::LitSet(values) => {
            for &value in values {
                check_coset_index(entity, kind, value)?;
            }
            Ok(())
        }
        StereoCoset::Term(term) => check_term(entity, kind, term),
    }
}

fn check_term(
    entity: Entity,
    kind: StereoKind,
    term: &StereoTerm,
) -> Result<(), MoleculeIntegrityError> {
    match term {
        StereoTerm::Var(value) => {
            if let Some(domain) = &value.1 {
                for &coset in domain {
                    check_coset_index(entity, kind, coset)?;
                }
            }
            Ok(())
        }
        StereoTerm::Lit(value) => check_coset_index(entity, kind, *value),
        StereoTerm::LitSet(values) => {
            for &value in values {
                check_coset_index(entity, kind, value)?;
            }
            Ok(())
        }
        StereoTerm::Swap(inner) | StereoTerm::Mirror(inner) => check_term(entity, kind, inner),
        StereoTerm::Apply(inner, permutation) => {
            check_permutation(entity, kind.degree(), *permutation)?;
            check_term(entity, kind, inner)
        }
    }
}

fn check_coset_index(
    entity: Entity,
    kind: StereoKind,
    coset: u32,
) -> Result<(), MoleculeIntegrityError> {
    if coset as usize >= kind.count() {
        Err(MoleculeIntegrityError::StereoCosetOutOfRange {
            entity,
            kind,
            coset,
            count: kind.count(),
        })
    } else {
        Ok(())
    }
}

fn check_permutation(
    entity: Entity,
    expected: usize,
    permutation: Permutation,
) -> Result<(), MoleculeIntegrityError> {
    if permutation.degree() != expected {
        Err(MoleculeIntegrityError::StereoPermutationDegree {
            entity,
            expected,
            actual: permutation.degree(),
        })
    } else {
        Ok(())
    }
}

fn check_pair(
    entity: Entity,
    degree: usize,
    pair: StereoLigandPair,
) -> Result<(), MoleculeIntegrityError> {
    for position in [pair.first(), pair.second()] {
        if position.index() >= degree {
            return Err(MoleculeIntegrityError::StereoLigandPositionOutOfRange {
                entity,
                position: position.index(),
                degree,
            });
        }
    }
    Ok(())
}
