//! Representation-integrity checks for [`Molecule`].

use std::collections::HashSet;

use smallvec::SmallVec;
use thiserror::Error;
use umol_graph_core::NodeId;

use super::super::constraint::{Constraint, MoleculeConstraint, RelationalConstraint};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::id::{AtomId, BondId};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::stereo::integrity::{
    check_stereo_atom_constraint_on_frame, check_stereo_atom_entry,
    check_stereo_bond_constraint_on_frame, check_stereo_bond_entry, StereoIntegrityError,
};
use super::super::stereo::StereoKind;
use super::{Molecule, MoleculeEntries};

/// Failure of the representation contract required to interpret a [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeIntegrityError {
    #[error("molecule references unavailable {entity}")]
    InvalidReference { entity: Entity },
    #[error(
        "{entity}: electron-count vector has length {electron_counts}, expected {participants}"
    )]
    ElectronCountLengthMismatch {
        entity: Entity,
        participants: usize,
        electron_counts: usize,
    },
    #[error("{entity}: participant atom {atom:?} is duplicated")]
    DuplicateAtom { entity: Entity, atom: AtomId },
    #[error("bond: parallel bonds on atoms {atoms:?}")]
    ParallelBonds { atoms: [AtomId; 2] },
    #[error("dative bonds: identical acceptor {acceptor:?} and donor set {donors:?}")]
    IdenticalDativeBonds {
        acceptor: AtomId,
        donors: Vec<AtomId>,
    },
    #[error("noncovalent bond: parallel bonds on atoms {atoms:?}")]
    ParallelNoncovalentBonds { atoms: [AtomId; 2] },
    #[error("aromatic systems: overlap on atom {atom:?}")]
    AromaticSystemsOverlap { atom: AtomId },
    #[error("multicenter bonds: identical participant set {atoms:?}")]
    IdenticalMulticenterBonds { atoms: Vec<AtomId> },
    #[error("stereo atom: duplicate site {atom:?}")]
    DuplicateStereoAtomSites { atom: AtomId },
    #[error("stereo bond: duplicate site {bond:?}")]
    DuplicateStereoBondSites { bond: BondId },
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

impl From<StereoIntegrityError> for MoleculeIntegrityError {
    fn from(error: StereoIntegrityError) -> Self {
        match error {
            StereoIntegrityError::DuplicateAtom { entity, atom } => {
                Self::DuplicateAtom { entity, atom }
            }
            StereoIntegrityError::DuplicateStereoLigand { entity, ligand } => {
                Self::DuplicateStereoLigand { entity, ligand }
            }
            StereoIntegrityError::StereoFrameDegreeTooLarge {
                entity,
                degree,
                maximum,
            } => Self::StereoFrameDegreeTooLarge {
                entity,
                degree,
                maximum,
            },
            StereoIntegrityError::StereoKindSiteMismatch { entity, kind } => {
                Self::StereoKindSiteMismatch { entity, kind }
            }
            StereoIntegrityError::StereoLigandArity {
                entity,
                kind,
                expected,
                actual,
            } => Self::StereoLigandArity {
                entity,
                kind,
                expected,
                actual,
            },
            StereoIntegrityError::StereoCosetOutOfRange {
                entity,
                kind,
                coset,
                count,
            } => Self::StereoCosetOutOfRange {
                entity,
                kind,
                coset,
                count,
            },
            StereoIntegrityError::StereoPermutationDegree {
                entity,
                expected,
                actual,
            } => Self::StereoPermutationDegree {
                entity,
                expected,
                actual,
            },
            StereoIntegrityError::StereoLigandPositionOutOfRange {
                entity,
                position,
                degree,
            } => Self::StereoLigandPositionOutOfRange {
                entity,
                position,
                degree,
            },
        }
    }
}

impl Molecule {
    /// Check the representation invariants required to interpret this molecule.
    ///
    /// This checks stored references, parallel collection shapes, the fixed relation semantics of
    /// every entity kind, and kind-dependent stereo domains. It does not check chemistry or
    /// constraint satisfaction.
    pub(crate) fn check_integrity(&self) -> Result<(), MoleculeIntegrityError> {
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

        for view in self.bonds().iter() {
            let entity = Entity::Bond(view.id);
            let atoms = view.atom_ids();
            check_unique_pair(entity, atoms)?;
            let pair = unordered_pair(atoms);
            let upper = NodeId::from(pair[1]);
            let neighbors = self.graph.neighbors(NodeId::from(pair[0]));
            let position = neighbors.partition_point(|neighbor| neighbor.node < upper);
            if neighbors
                .get(position + 1)
                .is_some_and(|neighbor| neighbor.node == upper)
            {
                return Err(MoleculeIntegrityError::ParallelBonds { atoms: pair });
            }
        }

        check_dative_bonds(self, &contains)?;

        check_aromatic_systems(self, &contains)?;

        check_multicenter_bonds(self, &contains)?;

        {
            let mut noncovalent_pairs: SmallVec<[u64; 16]> =
                SmallVec::with_capacity(self.noncovalent_bonds.count());
            for view in self.noncovalent_bonds().iter() {
                let entity = Entity::NoncovalentBond(view.id);
                let atoms = view.atom_ids();
                require_references(&contains, atoms.into_iter().map(Entity::Atom))?;
                check_unique_pair(entity, atoms)?;
                let pair = unordered_pair(atoms);
                noncovalent_pairs.push((u64::from(pair[0].0) << 32) | u64::from(pair[1].0));
            }
            noncovalent_pairs.sort_unstable();
            for keys in noncovalent_pairs.windows(2) {
                if keys[0] == keys[1] {
                    let key = keys[0];
                    let pair = [AtomId((key >> 32) as u32), AtomId(key as u32)];
                    return Err(MoleculeIntegrityError::ParallelNoncovalentBonds { atoms: pair });
                }
            }
        }

        {
            let mut stereo_atom_sites = HashSet::with_capacity(self.stereo_atoms.count());
            for view in self.stereo_atoms().iter() {
                let entity = Entity::StereoAtom(view.id);
                let site = view.site_id();
                let ligand_frame = self.stereo_atoms.ligands(view.id);
                require_reference(&contains, Entity::Atom(site))?;
                require_references(
                    &contains,
                    ligand_frame
                        .iter()
                        .map(|ligand| Entity::Atom(ligand.atom_id)),
                )?;
                check_stereo_atom_entry(entity, site, ligand_frame, view.attributes)?;
                if !stereo_atom_sites.insert(site) {
                    return Err(MoleculeIntegrityError::DuplicateStereoAtomSites { atom: site });
                }
                if ligand_frame
                    .iter()
                    .any(|&ligand| !ligand_matches_site(self, ligand, site, None))
                {
                    return Err(MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity });
                }
            }
        }

        for view in self.stereo_bonds().iter() {
            let entity = Entity::StereoBond(view.id);
            let site = view.site_id();
            let ligand_frame = self.stereo_bonds.ligands(view.id);
            require_reference(&contains, Entity::Bond(site))?;
            require_references(
                &contains,
                ligand_frame
                    .iter()
                    .map(|ligand| Entity::Atom(ligand.atom_id)),
            )?;
            check_stereo_bond_entry(entity, ligand_frame, view.attributes)?;
            if self
                .stereo_bonds
                .incident_to_bond_ids(site)
                .next()
                .is_some_and(|id| id < view.id)
            {
                return Err(MoleculeIntegrityError::DuplicateStereoBondSites { bond: site });
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
            check_constraint_references(constraint, &contains)
                .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })?;
            check_molecule_constraint(self, constraint)?;
        }
        Ok(())
    }
}

fn check_dative_bonds(
    molecule: &Molecule,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), MoleculeIntegrityError> {
    if molecule.atoms.len() <= 128 {
        let mut identities = HashSet::with_capacity(molecule.dative_bonds.count());
        for view in molecule.dative_bonds().iter() {
            let entity = Entity::DativeBond(view.id);
            let acceptor = view.acceptor_id();
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            let mut donors = [0_u64; 2];
            for atom in view.donor_ids() {
                let (word, bit) = atom_bit(atom);
                if donors[word] & bit != 0 {
                    return Err(MoleculeIntegrityError::DuplicateAtom { entity, atom });
                }
                donors[word] |= bit;
            }
            let (word, bit) = atom_bit(acceptor);
            if donors[word] & bit != 0 {
                return Err(MoleculeIntegrityError::DuplicateAtom {
                    entity,
                    atom: acceptor,
                });
            }
            if !identities.insert((acceptor, donors)) {
                let mut donors: Vec<_> = view.donor_ids().collect();
                donors.sort_unstable();
                return Err(MoleculeIntegrityError::IdenticalDativeBonds { acceptor, donors });
            }
        }
    } else {
        let mut identities = HashSet::with_capacity(molecule.dative_bonds.count());
        for view in molecule.dative_bonds().iter() {
            let entity = Entity::DativeBond(view.id);
            let acceptor = view.acceptor_id();
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            let mut donors: Vec<_> = view.donor_ids().collect();
            donors.sort_unstable();
            if donors.windows(2).any(|pair| pair[0] == pair[1]) {
                check_unique_participants(entity, view.donor_ids())?;
            }
            if donors.binary_search(&acceptor).is_ok() {
                return Err(MoleculeIntegrityError::DuplicateAtom {
                    entity,
                    atom: acceptor,
                });
            }
            if !identities.insert((acceptor, donors)) {
                let mut donors: Vec<_> = view.donor_ids().collect();
                donors.sort_unstable();
                return Err(MoleculeIntegrityError::IdenticalDativeBonds { acceptor, donors });
            }
        }
    }
    Ok(())
}

fn check_aromatic_systems(
    molecule: &Molecule,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), MoleculeIntegrityError> {
    if molecule.atoms.len() <= 128 {
        let mut membership = [0_u64; 2];
        for view in molecule.aromatic_systems().iter() {
            let entity = Entity::AromaticSystem(view.id);
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            let mut row = [0_u64; 2];
            for atom in view.atom_ids() {
                let (word, bit) = atom_bit(atom);
                if row[word] & bit != 0 {
                    return Err(MoleculeIntegrityError::DuplicateAtom { entity, atom });
                }
                row[word] |= bit;
            }
            if row[0] & membership[0] != 0 || row[1] & membership[1] != 0 {
                for atom in view.atom_ids() {
                    let (word, bit) = atom_bit(atom);
                    if membership[word] & bit != 0 {
                        return Err(MoleculeIntegrityError::AromaticSystemsOverlap { atom });
                    }
                }
            }
            membership[0] |= row[0];
            membership[1] |= row[1];
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
    } else {
        let mut membership = HashSet::new();
        for view in molecule.aromatic_systems().iter() {
            let entity = Entity::AromaticSystem(view.id);
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            {
                let mut row: Vec<_> = view.atom_ids().collect();
                row.sort_unstable();
                if row.windows(2).any(|pair| pair[0] == pair[1]) {
                    check_unique_participants(entity, view.atom_ids())?;
                }
            }
            for atom in view.atom_ids() {
                if !membership.insert(atom) {
                    return Err(MoleculeIntegrityError::AromaticSystemsOverlap { atom });
                }
            }
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
    }
    Ok(())
}

fn check_multicenter_bonds(
    molecule: &Molecule,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), MoleculeIntegrityError> {
    if molecule.atoms.len() <= 128 {
        let mut identities = HashSet::with_capacity(molecule.multicenter_bonds.count());
        for view in molecule.multicenter_bonds().iter() {
            let entity = Entity::MulticenterBond(view.id);
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            let mut atoms = [0_u64; 2];
            for atom in view.atom_ids() {
                let (word, bit) = atom_bit(atom);
                if atoms[word] & bit != 0 {
                    return Err(MoleculeIntegrityError::DuplicateAtom { entity, atom });
                }
                atoms[word] |= bit;
            }
            if !identities.insert(atoms) {
                return Err(MoleculeIntegrityError::IdenticalMulticenterBonds {
                    atoms: view.atom_ids().collect(),
                });
            }
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
    } else {
        let mut identities = HashSet::with_capacity(molecule.multicenter_bonds.count());
        for view in molecule.multicenter_bonds().iter() {
            let entity = Entity::MulticenterBond(view.id);
            require_references(contains, view.atom_ids().map(Entity::Atom))?;
            let mut atoms: Vec<_> = view.atom_ids().collect();
            atoms.sort_unstable();
            if atoms.windows(2).any(|pair| pair[0] == pair[1]) {
                check_unique_participants(entity, view.atom_ids())?;
            }
            if !identities.insert(atoms) {
                return Err(MoleculeIntegrityError::IdenticalMulticenterBonds {
                    atoms: view.atom_ids().collect(),
                });
            }
            check_electron_count_length(
                entity,
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
    }
    Ok(())
}

fn atom_bit(atom: AtomId) -> (usize, u64) {
    (atom.index() / 64, 1_u64 << (atom.index() % 64))
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
        Err(MoleculeIntegrityError::DuplicateAtom {
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
            return Err(MoleculeIntegrityError::DuplicateAtom { entity, atom });
        }
    }
    Ok(())
}

fn require_reference(
    contains: &impl Fn(Entity) -> bool,
    entity: Entity,
) -> Result<(), MoleculeIntegrityError> {
    if contains(entity) {
        Ok(())
    } else {
        Err(MoleculeIntegrityError::InvalidReference { entity })
    }
}

fn require_references(
    contains: &impl Fn(Entity) -> bool,
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

/// Check stereo wrappers after the constraint's entity references have been checked.
pub(crate) fn check_molecule_constraint(
    molecule: &Molecule,
    constraint: &Constraint,
) -> Result<(), MoleculeIntegrityError> {
    match constraint {
        Constraint::StereoAtom(id, kind, constraint) => {
            let entity = Entity::StereoAtom(*id);
            let ligand_count = molecule.stereo_atom(*id).ligand_count();
            check_stereo_atom_constraint_on_frame(entity, ligand_count, *kind, constraint)
                .map_err(Into::into)
        }
        Constraint::StereoBond(id, kind, constraint) => {
            let entity = Entity::StereoBond(*id);
            let ligand_count = molecule.stereo_bond(*id).ligand_count();
            check_stereo_bond_constraint_on_frame(entity, ligand_count, *kind, constraint)
                .map_err(Into::into)
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

pub(crate) fn check_entry_references(
    entries: &MoleculeEntries,
) -> Result<(), MoleculeIntegrityError> {
    check_entry_references_inner(entries)
        .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })
}

fn check_entry_references_inner(entries: &MoleculeEntries) -> Result<(), Entity> {
    let contains = |entity| match entity {
        Entity::Atom(id) => id.index() < entries.atoms.len(),
        Entity::Bond(id) => id.index() < entries.bonds.len(),
        Entity::DativeBond(id) => id.index() < entries.dative.len(),
        Entity::AromaticSystem(id) => id.index() < entries.aromatic.len(),
        Entity::MulticenterBond(id) => id.index() < entries.multicenter.len(),
        Entity::NoncovalentBond(id) => id.index() < entries.noncovalent.len(),
        Entity::StereoAtom(id) => id.index() < entries.stereo_atoms.len(),
        Entity::StereoBond(id) => id.index() < entries.stereo_bonds.len(),
    };

    for &(first, second, _) in &entries.bonds {
        check_reference(&contains, Entity::Atom(first))?;
        check_reference(&contains, Entity::Atom(second))?;
    }
    for (donors, acceptor, _) in &entries.dative {
        check_references(&contains, donors.iter().copied().map(Entity::Atom))?;
        check_reference(&contains, Entity::Atom(*acceptor))?;
    }
    for (atoms, _) in &entries.aromatic {
        check_references(&contains, atoms.iter().copied().map(Entity::Atom))?;
    }
    for (atoms, _) in &entries.multicenter {
        check_references(&contains, atoms.iter().copied().map(Entity::Atom))?;
    }
    for (atoms, _) in &entries.noncovalent {
        check_references(&contains, atoms.iter().copied().map(Entity::Atom))?;
    }
    for (site, ligands, _) in &entries.stereo_atoms {
        check_reference(&contains, Entity::Atom(*site))?;
        check_references(
            &contains,
            ligands.iter().map(|ligand| Entity::Atom(ligand.atom_id)),
        )?;
    }
    for (site, ligands, _) in &entries.stereo_bonds {
        check_reference(&contains, Entity::Bond(*site))?;
        check_references(
            &contains,
            ligands.iter().map(|ligand| Entity::Atom(ligand.atom_id)),
        )?;
    }
    for constraint in entries.constraints.iter() {
        check_constraint_references(constraint, &contains)?;
    }
    Ok(())
}

fn check_reference(contains: &impl Fn(Entity) -> bool, entity: Entity) -> Result<(), Entity> {
    if contains(entity) {
        Ok(())
    } else {
        Err(entity)
    }
}

fn check_references(
    contains: &impl Fn(Entity) -> bool,
    entities: impl IntoIterator<Item = Entity>,
) -> Result<(), Entity> {
    for entity in entities {
        check_reference(contains, entity)?;
    }
    Ok(())
}

pub(crate) fn check_constraint_references(
    constraint: &Constraint,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        Constraint::Atom(id, _) => check_reference(contains, Entity::Atom(*id)),
        Constraint::Bond(id, _) => check_reference(contains, Entity::Bond(*id)),
        Constraint::DativeBond(id, _) => check_reference(contains, Entity::DativeBond(*id)),
        Constraint::AromaticSystem(id, _) => check_reference(contains, Entity::AromaticSystem(*id)),
        Constraint::MulticenterBond(id, _) => {
            check_reference(contains, Entity::MulticenterBond(*id))
        }
        Constraint::NoncovalentBond(id, _) => {
            check_reference(contains, Entity::NoncovalentBond(*id))
        }
        Constraint::StereoAtom(id, _, _) => check_reference(contains, Entity::StereoAtom(*id)),
        Constraint::StereoBond(id, _, _) => check_reference(contains, Entity::StereoBond(*id)),
        Constraint::Relational(constraint) => {
            check_relational_constraint_references(constraint, contains)
        }
        Constraint::Molecule(constraint) => {
            check_molecule_constraint_references(constraint, contains)
        }
        Constraint::And(constraints) | Constraint::Or(constraints) => {
            for constraint in constraints {
                check_constraint_references(constraint, contains)?;
            }
            Ok(())
        }
        Constraint::Not(constraint) => check_constraint_references(constraint, contains),
    }
}

fn check_relational_constraint_references(
    constraint: &RelationalConstraint,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        RelationalConstraint::DativeBondDonors { bond, atoms }
        | RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
            check_reference(contains, Entity::DativeBond(*bond))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::DativeBondDonor { bond, atom }
        | RelationalConstraint::DativeBondAcceptor { bond, atom } => {
            check_reference(contains, Entity::DativeBond(*bond))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::DativeBondAllDonors { bond, .. }
        | RelationalConstraint::DativeBondAnyDonor { bond, .. }
        | RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => {
            check_reference(contains, Entity::DativeBond(*bond))
        }
        RelationalConstraint::DativeBondParallels { dative, parallel } => {
            check_reference(contains, Entity::DativeBond(*dative))?;
            check_reference(contains, Entity::Bond(*parallel))
        }
        RelationalConstraint::AromaticSystemAtoms { system, atoms }
        | RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
            check_reference(contains, Entity::AromaticSystem(*system))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::AromaticSystemContains { system, atom } => {
            check_reference(contains, Entity::AromaticSystem(*system))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::AromaticSystemAllAtoms { system, .. }
        | RelationalConstraint::AromaticSystemAnyAtom { system, .. } => {
            check_reference(contains, Entity::AromaticSystem(*system))
        }
        RelationalConstraint::MulticenterBondAtoms { bond, atoms }
        | RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
            check_reference(contains, Entity::MulticenterBond(*bond))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::MulticenterBondContains { bond, atom } => {
            check_reference(contains, Entity::MulticenterBond(*bond))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::MulticenterBondAllAtoms { bond, .. }
        | RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => {
            check_reference(contains, Entity::MulticenterBond(*bond))
        }
        RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
            check_reference(contains, Entity::NoncovalentBond(*bond))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::NoncovalentBondContains { bond, atom } => {
            check_reference(contains, Entity::NoncovalentBond(*bond))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
            check_reference(contains, Entity::NoncovalentBond(*bond))
        }
        RelationalConstraint::StereoAtomSite { stereo_atom, atom }
        | RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
            check_reference(contains, Entity::StereoAtom(*stereo_atom))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
            check_reference(contains, Entity::StereoAtom(*stereo_atom))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::StereoAtomAllLigands { stereo_atom, .. }
        | RelationalConstraint::StereoAtomAnyLigand { stereo_atom, .. } => {
            check_reference(contains, Entity::StereoAtom(*stereo_atom))
        }
        RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
            check_reference(contains, Entity::StereoBond(*stereo_bond))?;
            check_reference(contains, Entity::Bond(*bond))
        }
        RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
            check_reference(contains, Entity::StereoBond(*stereo_bond))?;
            check_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
            check_reference(contains, Entity::StereoBond(*stereo_bond))?;
            check_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::StereoBondAllLigands { stereo_bond, .. }
        | RelationalConstraint::StereoBondAnyLigand { stereo_bond, .. } => {
            check_reference(contains, Entity::StereoBond(*stereo_bond))
        }
    }
}

fn check_molecule_constraint_references(
    constraint: &MoleculeConstraint,
    contains: &impl Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        MoleculeConstraint::ChargeSum { atoms, .. }
        | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
        | MoleculeConstraint::Connected { atoms } => {
            check_references(contains, atoms.iter().flatten().copied().map(Entity::Atom))
        }
        MoleculeConstraint::BondOrderSum { bonds, .. } => {
            check_references(contains, bonds.iter().flatten().copied().map(Entity::Bond))
        }
    }
}
