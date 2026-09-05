//! A partial bijection between two `Molecule` id spaces, per entity kind.
//!
//! The atom part is a `Correspondence<AtomId>`; the seven other entity kinds each carry a
//! `Correspondence` over their entity id. Valueless — pairing only; adding values and a direction
//! lifts it to a reaction span.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use umol_graph_core::{Compaction, Correspondence, CorrespondenceComposeError};

use super::compact::MoleculeCompaction;
use super::entity::{Entity, EntityKind};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::Molecule;
#[cfg(test)]
use super::molecule::MoleculeEntries;
use super::remap::MoleculeRemapping;

/// A per-entity partial bijection between two molecules: atoms, bonds, and the six overlay kinds.
/// The matched/unmatched reads of each component are those of its `Correspondence`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCorrespondence {
    atoms: Correspondence<AtomId>,
    bonds: Correspondence<BondId>,
    dative_bonds: Correspondence<DativeBondId>,
    aromatic_systems: Correspondence<AromaticSystemId>,
    multicenter_bonds: Correspondence<MulticenterBondId>,
    noncovalent_bonds: Correspondence<NoncovalentBondId>,
    stereo_atoms: Correspondence<StereoAtomId>,
    stereo_bonds: Correspondence<StereoBondId>,
}

impl From<&MoleculeRemapping> for MoleculeCorrespondence {
    /// Preserve all eight permutations and their source and target counts.
    fn from(remapping: &MoleculeRemapping) -> Self {
        let nodes = remapping.graph().nodes();
        let edges = remapping.graph().edges();
        Self::new(
            Correspondence::from_images(
                &(0..nodes.len())
                    .map(|idx| remapping.map_atom(AtomId::from(idx)))
                    .collect::<Vec<_>>(),
                nodes.len(),
            ),
            Correspondence::from_images(
                &(0..edges.len())
                    .map(|idx| remapping.map_bond(BondId::from(idx)))
                    .collect::<Vec<_>>(),
                edges.len(),
            ),
            remapping.dative_bonds().into(),
            remapping.aromatic_systems().into(),
            remapping.multicenter_bonds().into(),
            remapping.noncovalent_bonds().into(),
            remapping.stereo_atoms().into(),
            remapping.stereo_bonds().into(),
        )
    }
}

impl From<&MoleculeCompaction> for MoleculeCorrespondence {
    /// Preserve all eight source/result counts and order-preserving survivor pairs.
    fn from(compaction: &MoleculeCompaction) -> Self {
        let nodes = compaction.graph().nodes();
        let edges = compaction.graph().edges();
        Self::new(
            Correspondence::new(
                (0..nodes.source_count())
                    .filter_map(|idx| {
                        let id = AtomId::from(idx);
                        compaction.compact_atom(id).map(|image| (id, image))
                    })
                    .collect(),
                nodes.source_count(),
                nodes.result_count(),
            )
            .expect("compaction preserves an injective survivor mapping"),
            Correspondence::new(
                (0..edges.source_count())
                    .filter_map(|idx| {
                        let id = BondId::from(idx);
                        compaction.compact_bond(id).map(|image| (id, image))
                    })
                    .collect(),
                edges.source_count(),
                edges.result_count(),
            )
            .expect("compaction preserves an injective survivor mapping"),
            compaction.dative_bonds().into(),
            compaction.aromatic_systems().into(),
            compaction.multicenter_bonds().into(),
            compaction.noncovalent_bonds().into(),
            compaction.stereo_atoms().into(),
            compaction.stereo_bonds().into(),
        )
    }
}

impl MoleculeCorrespondence {
    /// A correspondence from its eight per-entity-kind correspondences (the fully materialized form).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        atoms: Correspondence<AtomId>,
        bonds: Correspondence<BondId>,
        dative_bonds: Correspondence<DativeBondId>,
        aromatic_systems: Correspondence<AromaticSystemId>,
        multicenter_bonds: Correspondence<MulticenterBondId>,
        noncovalent_bonds: Correspondence<NoncovalentBondId>,
        stereo_atoms: Correspondence<StereoAtomId>,
        stereo_bonds: Correspondence<StereoBondId>,
    ) -> Self {
        Self {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        }
    }

    /// A correspondence between empty molecules: no pairs and zero counts on both sides
    /// for every entity kind.
    pub const fn empty() -> Self {
        Self {
            atoms: Correspondence::empty(),
            bonds: Correspondence::empty(),
            dative_bonds: Correspondence::empty(),
            aromatic_systems: Correspondence::empty(),
            multicenter_bonds: Correspondence::empty(),
            noncovalent_bonds: Correspondence::empty(),
            stereo_atoms: Correspondence::empty(),
            stereo_bonds: Correspondence::empty(),
        }
    }

    /// Derive the full per-entity correspondence between `lhs` and `rhs` from their atom
    /// correspondence. Bonds are the induced edge correspondence; each overlay's lhs entities are
    /// matched to an rhs entity by their atom constituents mapped through `atoms`. An entity whose
    /// constituents are not all matched is unmatched. Stereo entities with different determined
    /// geometry kinds are also unmatched, so a geometry change is represented as removal and addition.
    ///
    /// Returns `None` when the atom correspondence is not compatible with the supplied molecule
    /// pair, or
    /// when entity incidence does not induce a unique right partner.
    pub fn induce(lhs: &Molecule, rhs: &Molecule, atoms: Correspondence<AtomId>) -> Option<Self> {
        if atoms.left_count() != lhs.atoms().count() || atoms.right_count() != rhs.atoms().count() {
            return None;
        }

        let bonds = induced_bonds(lhs, rhs, &atoms)?;
        let dative_bonds = induced_dative_bonds(lhs, rhs, &atoms)?;
        let aromatic_systems = induced_aromatic_systems(lhs, rhs, &atoms)?;
        let multicenter_bonds = induced_multicenter_bonds(lhs, rhs, &atoms)?;
        let noncovalent_bonds = induced_noncovalent_bonds(lhs, rhs, &atoms)?;

        let stereo_atoms = induce_by_key(
            lhs.stereo_atoms().iter().filter_map(|stereo| {
                let (Some(site), Some(ligands)) = (
                    map_atom(&atoms, stereo.site_id()),
                    map_ligands(&atoms, stereo.ligand_frame()),
                ) else {
                    return None;
                };
                Some((stereo.id, (site, sorted_ligands(ligands))))
            }),
            lhs.stereo_atoms().count(),
            rhs.stereo_atoms().iter().map(|stereo| {
                (
                    stereo.id,
                    (stereo.site_id(), sorted_ligands(stereo.ligand_frame())),
                )
            }),
            rhs.stereo_atoms().count(),
        )?;
        let stereo_atoms = Correspondence::new(
            stereo_atoms
                .matched_pairs()
                .iter()
                .copied()
                .filter(|&(left, right)| {
                    lhs.stereo_atom(left)
                        .attributes
                        .configuration
                        .kind()
                        .zip(rhs.stereo_atom(right).attributes.configuration.kind())
                        .is_none_or(|(left, right)| left == right)
                })
                .collect(),
            lhs.stereo_atoms().count(),
            rhs.stereo_atoms().count(),
        )
        .expect("filtering an induced correspondence preserves its partial bijection");

        let stereo_bonds = induce_by_key(
            lhs.stereo_bonds().iter().filter_map(|stereo| {
                let (Some(site), Some(ligands)) = (
                    bonds.right_of(stereo.site_id()),
                    map_ligands(&atoms, stereo.ligand_frame()),
                ) else {
                    return None;
                };
                Some((stereo.id, (site, sorted_ligands(ligands))))
            }),
            lhs.stereo_bonds().count(),
            rhs.stereo_bonds().iter().map(|stereo| {
                (
                    stereo.id,
                    (stereo.site_id(), sorted_ligands(stereo.ligand_frame())),
                )
            }),
            rhs.stereo_bonds().count(),
        )?;
        let stereo_bonds = Correspondence::new(
            stereo_bonds
                .matched_pairs()
                .iter()
                .copied()
                .filter(|&(left, right)| {
                    lhs.stereo_bond(left)
                        .attributes
                        .configuration
                        .kind()
                        .zip(rhs.stereo_bond(right).attributes.configuration.kind())
                        .is_none_or(|(left, right)| left == right)
                })
                .collect(),
            lhs.stereo_bonds().count(),
            rhs.stereo_bonds().count(),
        )
        .expect("filtering an induced correspondence preserves its partial bijection");

        Some(Self::new(
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        ))
    }

    /// Append unmatched ids to one right-hand entity domain, retaining all pair vectors.
    pub fn extend_right(mut self, kind: EntityKind, count: usize) -> Self {
        match kind {
            EntityKind::Atom => self.atoms = self.atoms.extend_right(count),
            EntityKind::Bond => self.bonds = self.bonds.extend_right(count),
            EntityKind::DativeBond => self.dative_bonds = self.dative_bonds.extend_right(count),
            EntityKind::AromaticSystem => {
                self.aromatic_systems = self.aromatic_systems.extend_right(count)
            }
            EntityKind::MulticenterBond => {
                self.multicenter_bonds = self.multicenter_bonds.extend_right(count)
            }
            EntityKind::NoncovalentBond => {
                self.noncovalent_bonds = self.noncovalent_bonds.extend_right(count)
            }
            EntityKind::StereoAtom => self.stereo_atoms = self.stereo_atoms.extend_right(count),
            EntityKind::StereoBond => self.stereo_bonds = self.stereo_bonds.extend_right(count),
        }
        self
    }

    /// Compact every right-hand entity domain, reusing its pair vector.
    ///
    /// Pairs whose right entity is removed are discarded.
    /// Only the removed atom/bond id lists are copied to adapt their graph index types.
    ///
    /// # Errors
    ///
    /// Returns the first source-count mismatch in entity-kind order. Consumes the receiver.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to composition with the compaction's correspondence.
    pub fn compact_right(
        self,
        compaction: &MoleculeCompaction,
    ) -> Result<Self, MoleculeCorrespondenceComposeError> {
        let atoms = Compaction::new(
            compaction.graph().nodes().source_count(),
            compaction
                .graph()
                .nodes()
                .removed()
                .iter()
                .copied()
                .map(AtomId::from)
                .collect(),
        )
        .expect("typed atom ids preserve the graph compaction");
        let bonds = Compaction::new(
            compaction.graph().edges().source_count(),
            compaction
                .graph()
                .edges()
                .removed()
                .iter()
                .copied()
                .map(BondId::from)
                .collect(),
        )
        .expect("typed bond ids preserve the graph compaction");
        Ok(Self::new(
            self.atoms.compact_right(&atoms).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Atom,
                    source,
                }
            })?,
            self.bonds.compact_right(&bonds).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Bond,
                    source,
                }
            })?,
            self.dative_bonds
                .compact_right(compaction.dative_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::DativeBond,
                    source,
                })?,
            self.aromatic_systems
                .compact_right(compaction.aromatic_systems())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::AromaticSystem,
                    source,
                })?,
            self.multicenter_bonds
                .compact_right(compaction.multicenter_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::MulticenterBond,
                    source,
                })?,
            self.noncovalent_bonds
                .compact_right(compaction.noncovalent_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::NoncovalentBond,
                    source,
                })?,
            self.stereo_atoms
                .compact_right(compaction.stereo_atoms())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoAtom,
                    source,
                })?,
            self.stereo_bonds
                .compact_right(compaction.stereo_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoBond,
                    source,
                })?,
        ))
    }

    /// Expand every right-hand entity domain through the inverse compaction.
    ///
    /// Restored ids remain unmatched; discarded pairs are not recreated.
    /// Only the removed atom/bond id lists are copied to adapt their graph index types.
    ///
    /// # Errors
    ///
    /// Returns the first result-count mismatch in entity-kind order. Consumes the receiver.
    ///
    /// # Semantic properties
    ///
    /// Equivalent to composition with the reversed compaction correspondence.
    pub fn uncompact_right(
        self,
        compaction: &MoleculeCompaction,
    ) -> Result<Self, MoleculeCorrespondenceComposeError> {
        let atoms = Compaction::new(
            compaction.graph().nodes().source_count(),
            compaction
                .graph()
                .nodes()
                .removed()
                .iter()
                .copied()
                .map(AtomId::from)
                .collect(),
        )
        .expect("typed atom ids preserve the graph compaction");
        let bonds = Compaction::new(
            compaction.graph().edges().source_count(),
            compaction
                .graph()
                .edges()
                .removed()
                .iter()
                .copied()
                .map(BondId::from)
                .collect(),
        )
        .expect("typed bond ids preserve the graph compaction");
        Ok(Self::new(
            self.atoms.uncompact_right(&atoms).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Atom,
                    source,
                }
            })?,
            self.bonds.uncompact_right(&bonds).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Bond,
                    source,
                }
            })?,
            self.dative_bonds
                .uncompact_right(compaction.dative_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::DativeBond,
                    source,
                })?,
            self.aromatic_systems
                .uncompact_right(compaction.aromatic_systems())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::AromaticSystem,
                    source,
                })?,
            self.multicenter_bonds
                .uncompact_right(compaction.multicenter_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::MulticenterBond,
                    source,
                })?,
            self.noncovalent_bonds
                .uncompact_right(compaction.noncovalent_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::NoncovalentBond,
                    source,
                })?,
            self.stereo_atoms
                .uncompact_right(compaction.stereo_atoms())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoAtom,
                    source,
                })?,
            self.stereo_bonds
                .uncompact_right(compaction.stereo_bonds())
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoBond,
                    source,
                })?,
        ))
    }

    /// Relational composition, per entity kind: `self` (lhs↔middle) followed by `other`
    /// (middle↔rhs), yielding a lhs↔rhs correspondence.
    ///
    /// # Errors
    /// Returns the first entity kind with incompatible intermediate counts, in constructor order.
    /// Equal counts do not establish intermediate molecule identity.
    pub fn compose(
        &self,
        other: &MoleculeCorrespondence,
    ) -> Result<MoleculeCorrespondence, MoleculeCorrespondenceComposeError> {
        Ok(MoleculeCorrespondence::new(
            self.atoms.compose(&other.atoms).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Atom,
                    source,
                }
            })?,
            self.bonds.compose(&other.bonds).map_err(|source| {
                MoleculeCorrespondenceComposeError {
                    kind: EntityKind::Bond,
                    source,
                }
            })?,
            self.dative_bonds
                .compose(&other.dative_bonds)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::DativeBond,
                    source,
                })?,
            self.aromatic_systems
                .compose(&other.aromatic_systems)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::AromaticSystem,
                    source,
                })?,
            self.multicenter_bonds
                .compose(&other.multicenter_bonds)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::MulticenterBond,
                    source,
                })?,
            self.noncovalent_bonds
                .compose(&other.noncovalent_bonds)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::NoncovalentBond,
                    source,
                })?,
            self.stereo_atoms
                .compose(&other.stereo_atoms)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoAtom,
                    source,
                })?,
            self.stereo_bonds
                .compose(&other.stereo_bonds)
                .map_err(|source| MoleculeCorrespondenceComposeError {
                    kind: EntityKind::StereoBond,
                    source,
                })?,
        ))
    }

    /// Compose molecule correspondences in iteration order. Returns `Ok(None)` for an empty input
    /// and `Ok(Some(value))` for a singleton.
    ///
    /// # Errors
    /// Returns the first composition error in iteration order.
    pub fn compose_all(
        correspondences: impl IntoIterator<Item = Self>,
    ) -> Result<Option<Self>, MoleculeCorrespondenceComposeError> {
        let mut correspondences = correspondences.into_iter();
        let Some(first) = correspondences.next() else {
            return Ok(None);
        };
        correspondences
            .try_fold(first, |left, right| left.compose(&right))
            .map(Some)
    }

    /// The inverse correspondence (rhs↔lhs), per entity kind: each component's `reverse`.
    pub fn reverse(&self) -> MoleculeCorrespondence {
        MoleculeCorrespondence::new(
            self.atoms.reverse(),
            self.bonds.reverse(),
            self.dative_bonds.reverse(),
            self.aromatic_systems.reverse(),
            self.multicenter_bonds.reverse(),
            self.noncovalent_bonds.reverse(),
            self.stereo_atoms.reverse(),
            self.stereo_bonds.reverse(),
        )
    }

    /// The rhs entity matched to `left`, preserving its entity kind.
    pub fn right_of(&self, left: Entity) -> Option<Entity> {
        match left {
            Entity::Atom(id) => self.atoms.right_of(id).map(Entity::Atom),
            Entity::Bond(id) => self.bonds.right_of(id).map(Entity::Bond),
            Entity::DativeBond(id) => self.dative_bonds.right_of(id).map(Entity::DativeBond),
            Entity::AromaticSystem(id) => self
                .aromatic_systems
                .right_of(id)
                .map(Entity::AromaticSystem),
            Entity::MulticenterBond(id) => self
                .multicenter_bonds
                .right_of(id)
                .map(Entity::MulticenterBond),
            Entity::NoncovalentBond(id) => self
                .noncovalent_bonds
                .right_of(id)
                .map(Entity::NoncovalentBond),
            Entity::StereoAtom(id) => self.stereo_atoms.right_of(id).map(Entity::StereoAtom),
            Entity::StereoBond(id) => self.stereo_bonds.right_of(id).map(Entity::StereoBond),
        }
    }

    /// The lhs entity matched to `right`, preserving its entity kind.
    pub fn left_of(&self, right: Entity) -> Option<Entity> {
        match right {
            Entity::Atom(id) => self.atoms.left_of(id).map(Entity::Atom),
            Entity::Bond(id) => self.bonds.left_of(id).map(Entity::Bond),
            Entity::DativeBond(id) => self.dative_bonds.left_of(id).map(Entity::DativeBond),
            Entity::AromaticSystem(id) => self
                .aromatic_systems
                .left_of(id)
                .map(Entity::AromaticSystem),
            Entity::MulticenterBond(id) => self
                .multicenter_bonds
                .left_of(id)
                .map(Entity::MulticenterBond),
            Entity::NoncovalentBond(id) => self
                .noncovalent_bonds
                .left_of(id)
                .map(Entity::NoncovalentBond),
            Entity::StereoAtom(id) => self.stereo_atoms.left_of(id).map(Entity::StereoAtom),
            Entity::StereoBond(id) => self.stereo_bonds.left_of(id).map(Entity::StereoBond),
        }
    }

    /// Whether every id of all eight entity kinds on the left is matched.
    pub fn is_total_on_left(&self) -> bool {
        self.atoms.is_total_on_left()
            && self.bonds.is_total_on_left()
            && self.dative_bonds.is_total_on_left()
            && self.aromatic_systems.is_total_on_left()
            && self.multicenter_bonds.is_total_on_left()
            && self.noncovalent_bonds.is_total_on_left()
            && self.stereo_atoms.is_total_on_left()
            && self.stereo_bonds.is_total_on_left()
    }

    /// Whether every id of all eight entity kinds on the right is matched.
    pub fn is_total_on_right(&self) -> bool {
        self.atoms.is_total_on_right()
            && self.bonds.is_total_on_right()
            && self.dative_bonds.is_total_on_right()
            && self.aromatic_systems.is_total_on_right()
            && self.multicenter_bonds.is_total_on_right()
            && self.noncovalent_bonds.is_total_on_right()
            && self.stereo_atoms.is_total_on_right()
            && self.stereo_bonds.is_total_on_right()
    }

    /// Whether every id of all eight entity kinds is matched on both sides.
    pub fn is_total(&self) -> bool {
        self.is_total_on_left() && self.is_total_on_right()
    }

    /// Whether this correspondence actually relates `lhs` to `rhs`: every component's declared counts
    /// match the molecules, and every matched pair's participants map onto its counterpart's under
    /// the atom correspondence.
    ///
    /// A structural property of the correspondence, asked once. It is not a value comparison — no
    /// payload is read — so an operation that compares values under a correspondence establishes
    /// this first and then only compares values.
    pub fn is_compatible(&self, lhs: &Molecule, rhs: &Molecule) -> bool {
        let correspondence = self;
        let counts = [
            (
                correspondence.atoms().left_count(),
                lhs.atoms().count(),
                correspondence.atoms().right_count(),
                rhs.atoms().count(),
            ),
            (
                correspondence.bonds().left_count(),
                lhs.bonds().count(),
                correspondence.bonds().right_count(),
                rhs.bonds().count(),
            ),
            (
                correspondence.dative_bonds().left_count(),
                lhs.dative_bonds().count(),
                correspondence.dative_bonds().right_count(),
                rhs.dative_bonds().count(),
            ),
            (
                correspondence.aromatic_systems().left_count(),
                lhs.aromatic_systems().count(),
                correspondence.aromatic_systems().right_count(),
                rhs.aromatic_systems().count(),
            ),
            (
                correspondence.multicenter_bonds().left_count(),
                lhs.multicenter_bonds().count(),
                correspondence.multicenter_bonds().right_count(),
                rhs.multicenter_bonds().count(),
            ),
            (
                correspondence.noncovalent_bonds().left_count(),
                lhs.noncovalent_bonds().count(),
                correspondence.noncovalent_bonds().right_count(),
                rhs.noncovalent_bonds().count(),
            ),
            (
                correspondence.stereo_atoms().left_count(),
                lhs.stereo_atoms().count(),
                correspondence.stereo_atoms().right_count(),
                rhs.stereo_atoms().count(),
            ),
            (
                correspondence.stereo_bonds().left_count(),
                lhs.stereo_bonds().count(),
                correspondence.stereo_bonds().right_count(),
                rhs.stereo_bonds().count(),
            ),
        ];
        if counts.into_iter().any(
            |(declared_left, actual_left, declared_right, actual_right)| {
                declared_left != actual_left || declared_right != actual_right
            },
        ) {
            return false;
        }

        let atoms = correspondence.atoms();
        let same_atom_set = |left: Vec<AtomId>, mut right: Vec<AtomId>| {
            let Some(mut mapped): Option<Vec<_>> =
                left.into_iter().map(|atom| atoms.right_of(atom)).collect()
            else {
                return false;
            };
            mapped.sort_unstable();
            right.sort_unstable();
            mapped == right
        };
        let same_ligand_set = |left: Vec<StereoLigand>, mut right: Vec<StereoLigand>| {
            let Some(mut mapped): Option<Vec<_>> = left
                .into_iter()
                .map(|ligand| {
                    atoms
                        .right_of(ligand.atom_id)
                        .map(|atom| StereoLigand::new(atom, ligand.kind))
                })
                .collect()
            else {
                return false;
            };
            mapped.sort_unstable();
            right.sort_unstable();
            mapped == right
        };

        if !correspondence
            .bonds()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                same_atom_set(
                    lhs.bond(left).atom_ids().to_vec(),
                    rhs.bond(right).atom_ids().to_vec(),
                )
            })
        {
            return false;
        }
        if !correspondence
            .dative_bonds()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                let lhs = lhs.dative_bond(left);
                let rhs = rhs.dative_bond(right);
                atoms.right_of(lhs.acceptor_id()) == Some(rhs.acceptor_id())
                    && same_atom_set(lhs.donor_ids().collect(), rhs.donor_ids().collect())
            })
        {
            return false;
        }
        if !correspondence
            .aromatic_systems()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                same_atom_set(
                    lhs.aromatic_system(left).atom_ids().collect(),
                    rhs.aromatic_system(right).atom_ids().collect(),
                )
            })
        {
            return false;
        }
        if !correspondence
            .multicenter_bonds()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                same_atom_set(
                    lhs.multicenter_bond(left).atom_ids().collect(),
                    rhs.multicenter_bond(right).atom_ids().collect(),
                )
            })
        {
            return false;
        }
        if !correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                same_atom_set(
                    lhs.noncovalent_bond(left).atom_ids().to_vec(),
                    rhs.noncovalent_bond(right).atom_ids().to_vec(),
                )
            })
        {
            return false;
        }
        if !correspondence
            .stereo_atoms()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                let lhs = lhs.stereo_atom(left);
                let rhs = rhs.stereo_atom(right);
                atoms.right_of(lhs.site_id()) == Some(rhs.site_id())
                    && same_ligand_set(lhs.ligand_frame(), rhs.ligand_frame())
            })
        {
            return false;
        }
        correspondence
            .stereo_bonds()
            .matched_pairs()
            .iter()
            .all(|&(left, right)| {
                let lhs = lhs.stereo_bond(left);
                let rhs = rhs.stereo_bond(right);
                correspondence.bonds().right_of(lhs.site_id()) == Some(rhs.site_id())
                    && same_ligand_set(lhs.ligand_frame(), rhs.ligand_frame())
            })
    }

    /// The atom correspondence — the spine the other entity correspondences are induced from.
    pub fn atoms(&self) -> &Correspondence<AtomId> {
        &self.atoms
    }

    /// The bond correspondence.
    pub fn bonds(&self) -> &Correspondence<BondId> {
        &self.bonds
    }

    /// The dative-bond correspondence.
    pub fn dative_bonds(&self) -> &Correspondence<DativeBondId> {
        &self.dative_bonds
    }

    /// The aromatic-system correspondence.
    pub fn aromatic_systems(&self) -> &Correspondence<AromaticSystemId> {
        &self.aromatic_systems
    }

    /// The multicenter-bond correspondence.
    pub fn multicenter_bonds(&self) -> &Correspondence<MulticenterBondId> {
        &self.multicenter_bonds
    }

    /// The noncovalent-bond correspondence.
    pub fn noncovalent_bonds(&self) -> &Correspondence<NoncovalentBondId> {
        &self.noncovalent_bonds
    }

    /// The stereo-atom correspondence.
    pub fn stereo_atoms(&self) -> &Correspondence<StereoAtomId> {
        &self.stereo_atoms
    }

    /// The stereo-bond correspondence.
    pub fn stereo_bonds(&self) -> &Correspondence<StereoBondId> {
        &self.stereo_bonds
    }
}

/// A molecule correspondence component has incompatible intermediate counts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCorrespondenceComposeError {
    pub kind: EntityKind,
    pub source: CorrespondenceComposeError,
}

impl Display for MoleculeCorrespondenceComposeError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.kind, self.source)
    }
}

impl Error for MoleculeCorrespondenceComposeError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

/// The bond correspondence induced by an atom correspondence: the two molecular graphs' edge
/// correspondence under `atoms`.
pub(crate) fn induced_bonds(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<BondId>> {
    induce_by_key(
        left.bonds().iter().filter_map(|bond| {
            let [first, second] = bond.atom_ids();
            Some((
                bond.id,
                ordered_pair(atoms.right_of(first)?, atoms.right_of(second)?),
            ))
        }),
        left.bonds().count(),
        right
            .bonds()
            .iter()
            .map(|bond| (bond.id, ordered_pair_from(bond.atom_ids()))),
        right.bonds().count(),
    )
}

/// The dative-bond correspondence induced by an atom correspondence: each left dative bond whose
/// acceptor and donors are all matched with the right dative bond over the same roles.
pub(crate) fn induced_dative_bonds(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<DativeBondId>> {
    induce_by_key(
        left.dative_bonds().iter().filter_map(|dative| {
            let (Some(acceptor), Some(donors)) = (
                map_atom(atoms, dative.acceptor_id()),
                map_atoms(atoms, dative.donor_ids()),
            ) else {
                return None;
            };
            Some((dative.id, (acceptor, sorted_atoms(donors))))
        }),
        left.dative_bonds().count(),
        right.dative_bonds().iter().map(|dative| {
            (
                dative.id,
                (dative.acceptor_id(), sorted_atoms(dative.donor_ids())),
            )
        }),
        right.dative_bonds().count(),
    )
}

/// The aromatic-system correspondence induced by an atom correspondence: each left system whose
/// atoms are all matched with the right system over the same atom set.
pub(crate) fn induced_aromatic_systems(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<AromaticSystemId>> {
    induce_by_key(
        left.aromatic_systems().iter().filter_map(|aromatic| {
            Some((
                aromatic.id,
                sorted_atoms(map_atoms(atoms, aromatic.atom_ids())?),
            ))
        }),
        left.aromatic_systems().count(),
        right
            .aromatic_systems()
            .iter()
            .map(|aromatic| (aromatic.id, sorted_atoms(aromatic.atom_ids()))),
        right.aromatic_systems().count(),
    )
}

/// The multicenter-bond correspondence induced by an atom correspondence: each left bond whose
/// atoms are all matched with the right bond over the same atom set.
pub(crate) fn induced_multicenter_bonds(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<MulticenterBondId>> {
    induce_by_key(
        left.multicenter_bonds().iter().filter_map(|multicenter| {
            Some((
                multicenter.id,
                sorted_atoms(map_atoms(atoms, multicenter.atom_ids())?),
            ))
        }),
        left.multicenter_bonds().count(),
        right
            .multicenter_bonds()
            .iter()
            .map(|multicenter| (multicenter.id, sorted_atoms(multicenter.atom_ids()))),
        right.multicenter_bonds().count(),
    )
}

/// The noncovalent-bond correspondence induced by an atom correspondence: each left bond whose two
/// atoms are both matched with the right bond over the same atom pair.
pub(crate) fn induced_noncovalent_bonds(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<NoncovalentBondId>> {
    induce_by_key(
        left.noncovalent_bonds().iter().filter_map(|noncovalent| {
            let [first, second] = noncovalent.atom_ids();
            Some((
                noncovalent.id,
                ordered_pair(map_atom(atoms, first)?, map_atom(atoms, second)?),
            ))
        }),
        left.noncovalent_bonds().count(),
        right
            .noncovalent_bonds()
            .iter()
            .map(|noncovalent| (noncovalent.id, ordered_pair_from(noncovalent.atom_ids()))),
        right.noncovalent_bonds().count(),
    )
}

/// The rhs partner of a lhs atom under the atom correspondence, if matched.
pub(crate) fn map_atom(atoms: &Correspondence<AtomId>, atom: AtomId) -> Option<AtomId> {
    atoms.right_of(atom)
}

/// The rhs partners of a set of lhs atoms, or `None` if any is unmatched.
fn map_atoms(
    atoms: &Correspondence<AtomId>,
    lhs: impl IntoIterator<Item = AtomId>,
) -> Option<Vec<AtomId>> {
    lhs.into_iter().map(|atom| map_atom(atoms, atom)).collect()
}

/// The rhs-frame ligands (each ligand's atom mapped, its kind kept), or `None` if any ligand's
/// atom is unmatched.
pub(crate) fn map_ligands(
    atoms: &Correspondence<AtomId>,
    ligands: Vec<StereoLigand>,
) -> Option<Vec<StereoLigand>> {
    ligands
        .into_iter()
        .map(|ligand| {
            map_atom(atoms, ligand.atom_id).map(|atom| StereoLigand::new(atom, ligand.kind))
        })
        .collect()
}

fn induce_by_key<Id, Key>(
    left: impl IntoIterator<Item = (Id, Key)>,
    left_count: usize,
    right: impl IntoIterator<Item = (Id, Key)>,
    right_count: usize,
) -> Option<Correspondence<Id>>
where
    Id: Copy + Ord + From<usize>,
    Key: Ord,
{
    let mut right_by_key = BTreeMap::new();
    for (id, key) in right {
        right_by_key
            .entry(key)
            .and_modify(|entry| *entry = None)
            .or_insert(Some(id));
    }

    let mut matched_pairs = Vec::new();
    for (left, key) in left {
        match right_by_key.get(&key) {
            Some(Some(right)) => matched_pairs.push((left, *right)),
            Some(None) => return None,
            None => {}
        }
    }
    Correspondence::new(matched_pairs, left_count, right_count).ok()
}

fn ordered_pair<Id: Ord>(first: Id, second: Id) -> [Id; 2] {
    if first <= second {
        [first, second]
    } else {
        [second, first]
    }
}

fn ordered_pair_from<Id: Ord>([first, second]: [Id; 2]) -> [Id; 2] {
    ordered_pair(first, second)
}

fn sorted_atoms(atoms: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
    let mut atoms: Vec<_> = atoms.into_iter().collect();
    atoms.sort_unstable();
    atoms
}

fn sorted_ligands(ligands: impl IntoIterator<Item = StereoLigand>) -> Vec<StereoLigand> {
    let mut ligands: Vec<_> = ligands.into_iter().collect();
    ligands.sort_unstable();
    ligands
}

#[cfg(test)]
mod tests {

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{Compaction, EdgeId, GraphCompaction, NodeId};

    use super::*;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::Constraints;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::NoncovalentBondForm;
    use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoKind};

    #[fixture]
    fn cascade_molecule() -> Molecule {
        let mut editor = Molecule::default().edit();
        for _ in 0..12 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        for base in [0, 6] {
            for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (0, 4)] {
                editor.add_bond(AtomId(base + a), AtomId(base + b), BondForm::from_order(1));
            }
            editor.add_dative_bond(
                vec![AtomId(base)],
                AtomId(base + 1),
                DativeBondForm::from_order(1),
            );
            editor.add_aromatic_system(
                vec![AtomId(base), AtomId(base + 1), AtomId(base + 2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 1]),
            );
            editor.add_multicenter_bond(
                vec![AtomId(base), AtomId(base + 1), AtomId(base + 3)],
                MulticenterBondForm::from_electrons(vec![1, 0, 1]),
            );
            editor.add_noncovalent_bond(
                [AtomId(base), AtomId(base + 5)],
                NoncovalentBondForm::default(),
            );
            editor.add_stereo_atom(
                AtomId(base),
                [1, 2, 3, 4]
                    .map(|idx| StereoLigand::new(AtomId(base + idx), StereoLigandKind::Atom))
                    .to_vec(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            );
            editor.add_stereo_bond(
                BondId(base),
                [2, 3, 4, 5]
                    .map(|idx| StereoLigand::new(AtomId(base + idx), StereoLigandKind::Atom))
                    .to_vec(),
                StereoBondForm::default(),
            );
        }
        editor.build()
    }

    #[rstest]
    #[case::first_component(0, 6, 1)]
    #[case::last_component(6, 0, 0)]
    fn test_molecule_correspondence_from_compaction_cascade(
        cascade_molecule: Molecule,
        #[case] removed_start: u32,
        #[case] survivor_start: u32,
        #[case] surviving_overlay: u32,
    ) {
        let source = cascade_molecule;
        let mut editor = source.clone().edit();
        let removed = (removed_start..removed_start + 6)
            .map(AtomId)
            .collect::<Vec<_>>();
        let compaction = editor.tracked_remove(&removed, &[]);
        let result = editor.build();
        let mut plain = source.edit();
        plain.remove(&removed, &[]);
        assert_eq!(plain.build(), result);
        let expected_compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(12, (removed_start..removed_start + 6).map(NodeId).collect())
                    .unwrap(),
                Compaction::new(12, (removed_start..removed_start + 6).map(EdgeId).collect())
                    .unwrap(),
            ),
            Compaction::new(2, vec![DativeBondId(1 - surviving_overlay)]).unwrap(),
            Compaction::new(2, vec![AromaticSystemId(1 - surviving_overlay)]).unwrap(),
            Compaction::new(2, vec![MulticenterBondId(1 - surviving_overlay)]).unwrap(),
            Compaction::new(2, vec![NoncovalentBondId(1 - surviving_overlay)]).unwrap(),
            Compaction::new(2, vec![StereoAtomId(1 - surviving_overlay)]).unwrap(),
            Compaction::new(2, vec![StereoBondId(1 - surviving_overlay)]).unwrap(),
        );
        assert_eq!(compaction, expected_compaction);
        let correspondence = MoleculeCorrespondence::from(&compaction);
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (0, 4)]
                .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                .to_vec(),
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 1]),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![1, 0, 1]),
            )],
            noncovalent: vec![([AtomId(0), AtomId(5)], NoncovalentBondForm::default())],
            stereo_atoms: vec![(
                AtomId(0),
                [1, 2, 3, 4]
                    .map(|idx| StereoLigand::new(AtomId(idx), StereoLigandKind::Atom))
                    .to_vec(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            )],
            stereo_bonds: vec![(
                BondId(0),
                [2, 3, 4, 5]
                    .map(|idx| StereoLigand::new(AtomId(idx), StereoLigandKind::Atom))
                    .to_vec(),
                StereoBondForm::default(),
            )],
            constraints: Constraints::default(),
        });
        assert_eq!(result, expected);
        let expected_correspondence = MoleculeCorrespondence::new(
            Correspondence::new(
                (0..6)
                    .map(|idx| (AtomId(survivor_start + idx), AtomId(idx)))
                    .collect(),
                source.atoms().count(),
                result.atoms().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..6)
                    .map(|idx| (BondId(survivor_start + idx), BondId(idx)))
                    .collect(),
                source.bonds().count(),
                result.bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(DativeBondId(surviving_overlay), DativeBondId(0))],
                source.dative_bonds().count(),
                result.dative_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(AromaticSystemId(surviving_overlay), AromaticSystemId(0))],
                source.aromatic_systems().count(),
                result.aromatic_systems().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(MulticenterBondId(surviving_overlay), MulticenterBondId(0))],
                source.multicenter_bonds().count(),
                result.multicenter_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(NoncovalentBondId(surviving_overlay), NoncovalentBondId(0))],
                source.noncovalent_bonds().count(),
                result.noncovalent_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(StereoAtomId(surviving_overlay), StereoAtomId(0))],
                source.stereo_atoms().count(),
                result.stereo_atoms().count(),
            )
            .unwrap(),
            Correspondence::new(
                vec![(StereoBondId(surviving_overlay), StereoBondId(0))],
                source.stereo_bonds().count(),
                result.stereo_bonds().count(),
            )
            .unwrap(),
        );
        assert_eq!(correspondence, expected_correspondence);
    }

    #[rstest]
    #[case::empty(false, true)]
    #[case::identity(true, false)]
    #[case::full_removal(true, true)]
    fn test_molecule_correspondence_from_compaction_boundary(
        cascade_molecule: Molecule,
        #[case] populated: bool,
        #[case] remove_all: bool,
    ) {
        let source = if populated {
            cascade_molecule
        } else {
            Molecule::default()
        };
        let mut editor = source.clone().edit();
        let removed = if remove_all {
            (0..source.atoms().count()).map(AtomId::from).collect()
        } else {
            Vec::new()
        };
        let compaction = editor.tracked_remove(&removed, &[]);
        let result = editor.build();
        assert_eq!(
            result,
            if remove_all {
                Molecule::default()
            } else {
                source.clone()
            }
        );
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(
                (0..result.atoms().count())
                    .map(|idx| (AtomId::from(idx), AtomId::from(idx)))
                    .collect(),
                source.atoms().count(),
                result.atoms().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.bonds().count())
                    .map(|idx| (BondId::from(idx), BondId::from(idx)))
                    .collect(),
                source.bonds().count(),
                result.bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.dative_bonds().count())
                    .map(|idx| (DativeBondId::from(idx), DativeBondId::from(idx)))
                    .collect(),
                source.dative_bonds().count(),
                result.dative_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.aromatic_systems().count())
                    .map(|idx| (AromaticSystemId::from(idx), AromaticSystemId::from(idx)))
                    .collect(),
                source.aromatic_systems().count(),
                result.aromatic_systems().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.multicenter_bonds().count())
                    .map(|idx| (MulticenterBondId::from(idx), MulticenterBondId::from(idx)))
                    .collect(),
                source.multicenter_bonds().count(),
                result.multicenter_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.noncovalent_bonds().count())
                    .map(|idx| (NoncovalentBondId::from(idx), NoncovalentBondId::from(idx)))
                    .collect(),
                source.noncovalent_bonds().count(),
                result.noncovalent_bonds().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.stereo_atoms().count())
                    .map(|idx| (StereoAtomId::from(idx), StereoAtomId::from(idx)))
                    .collect(),
                source.stereo_atoms().count(),
                result.stereo_atoms().count(),
            )
            .unwrap(),
            Correspondence::new(
                (0..result.stereo_bonds().count())
                    .map(|idx| (StereoBondId::from(idx), StereoBondId::from(idx)))
                    .collect(),
                source.stereo_bonds().count(),
                result.stereo_bonds().count(),
            )
            .unwrap(),
        );
        assert_eq!(MoleculeCorrespondence::from(&compaction), expected);
    }

    #[fixture]
    fn correspondence() -> MoleculeCorrespondence {
        // distinct pairs per entity kind so a mis-wired accessor is caught.
        MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(1))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(BondId(0), BondId(2))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(DativeBondId(0), DativeBondId(3))], 2, 4)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(4))], 2, 5)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(5))], 2, 6)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(6))], 2, 7)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(7))], 2, 8)
                .expect("correspondence producer preserves partial-bijection invariants"),
            Correspondence::new(vec![(StereoBondId(0), StereoBondId(8))], 2, 9)
                .expect("correspondence producer preserves partial-bijection invariants"),
        )
    }

    #[fixture]
    fn update_correspondence() -> MoleculeCorrespondence {
        MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(2)), (AtomId(2), AtomId(0))], 3, 4)
                .unwrap(),
            Correspondence::new(vec![(BondId(0), BondId(2)), (BondId(2), BondId(0))], 3, 5)
                .unwrap(),
            Correspondence::new(
                vec![
                    (DativeBondId(0), DativeBondId(2)),
                    (DativeBondId(2), DativeBondId(0)),
                ],
                3,
                6,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (AromaticSystemId(0), AromaticSystemId(2)),
                    (AromaticSystemId(2), AromaticSystemId(0)),
                ],
                3,
                7,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(0), MulticenterBondId(2)),
                    (MulticenterBondId(2), MulticenterBondId(0)),
                ],
                3,
                8,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (NoncovalentBondId(0), NoncovalentBondId(2)),
                    (NoncovalentBondId(2), NoncovalentBondId(0)),
                ],
                3,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(2)),
                    (StereoAtomId(2), StereoAtomId(0)),
                ],
                3,
                10,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoBondId(0), StereoBondId(2)),
                    (StereoBondId(2), StereoBondId(0)),
                ],
                3,
                11,
            )
            .unwrap(),
        )
    }

    #[rstest]
    #[case::atoms(EntityKind::Atom)]
    #[case::bonds(EntityKind::Bond)]
    #[case::dative_bonds(EntityKind::DativeBond)]
    #[case::aromatic_systems(EntityKind::AromaticSystem)]
    #[case::multicenter_bonds(EntityKind::MulticenterBond)]
    #[case::noncovalent_bonds(EntityKind::NoncovalentBond)]
    #[case::stereo_atoms(EntityKind::StereoAtom)]
    #[case::stereo_bonds(EntityKind::StereoBond)]
    fn test_molecule_correspondence_extend_right(
        update_correspondence: MoleculeCorrespondence,
        #[case] kind: EntityKind,
    ) {
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(
                vec![(AtomId(0), AtomId(2)), (AtomId(2), AtomId(0))],
                3,
                4 + if kind == EntityKind::Atom { 2 } else { 0 },
            )
            .unwrap(),
            Correspondence::new(
                vec![(BondId(0), BondId(2)), (BondId(2), BondId(0))],
                3,
                5 + if kind == EntityKind::Bond { 2 } else { 0 },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (DativeBondId(0), DativeBondId(2)),
                    (DativeBondId(2), DativeBondId(0)),
                ],
                3,
                6 + if kind == EntityKind::DativeBond { 2 } else { 0 },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (AromaticSystemId(0), AromaticSystemId(2)),
                    (AromaticSystemId(2), AromaticSystemId(0)),
                ],
                3,
                7 + if kind == EntityKind::AromaticSystem {
                    2
                } else {
                    0
                },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(0), MulticenterBondId(2)),
                    (MulticenterBondId(2), MulticenterBondId(0)),
                ],
                3,
                8 + if kind == EntityKind::MulticenterBond {
                    2
                } else {
                    0
                },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (NoncovalentBondId(0), NoncovalentBondId(2)),
                    (NoncovalentBondId(2), NoncovalentBondId(0)),
                ],
                3,
                9 + if kind == EntityKind::NoncovalentBond {
                    2
                } else {
                    0
                },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(2)),
                    (StereoAtomId(2), StereoAtomId(0)),
                ],
                3,
                10 + if kind == EntityKind::StereoAtom { 2 } else { 0 },
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoBondId(0), StereoBondId(2)),
                    (StereoBondId(2), StereoBondId(0)),
                ],
                3,
                11 + if kind == EntityKind::StereoBond { 2 } else { 0 },
            )
            .unwrap(),
        );
        let atoms_ptr = update_correspondence.atoms().matched_pairs().as_ptr();
        let bonds_ptr = update_correspondence.bonds().matched_pairs().as_ptr();
        let dative_bonds_ptr = update_correspondence
            .dative_bonds()
            .matched_pairs()
            .as_ptr();
        let aromatic_systems_ptr = update_correspondence
            .aromatic_systems()
            .matched_pairs()
            .as_ptr();
        let multicenter_bonds_ptr = update_correspondence
            .multicenter_bonds()
            .matched_pairs()
            .as_ptr();
        let noncovalent_bonds_ptr = update_correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .as_ptr();
        let stereo_atoms_ptr = update_correspondence
            .stereo_atoms()
            .matched_pairs()
            .as_ptr();
        let stereo_bonds_ptr = update_correspondence
            .stereo_bonds()
            .matched_pairs()
            .as_ptr();
        let result = update_correspondence.extend_right(kind, 2);
        assert_eq!(result.atoms().matched_pairs().as_ptr(), atoms_ptr);
        assert_eq!(result.bonds().matched_pairs().as_ptr(), bonds_ptr);
        assert_eq!(
            result.dative_bonds().matched_pairs().as_ptr(),
            dative_bonds_ptr
        );
        assert_eq!(
            result.aromatic_systems().matched_pairs().as_ptr(),
            aromatic_systems_ptr
        );
        assert_eq!(
            result.multicenter_bonds().matched_pairs().as_ptr(),
            multicenter_bonds_ptr
        );
        assert_eq!(
            result.noncovalent_bonds().matched_pairs().as_ptr(),
            noncovalent_bonds_ptr
        );
        assert_eq!(
            result.stereo_atoms().matched_pairs().as_ptr(),
            stereo_atoms_ptr
        );
        assert_eq!(
            result.stereo_bonds().matched_pairs().as_ptr(),
            stereo_bonds_ptr
        );
        assert_eq!(result, expected);
    }

    #[rstest]
    fn test_molecule_correspondence_compact_right(update_correspondence: MoleculeCorrespondence) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(4, vec![NodeId(1)]).unwrap(),
                Compaction::new(5, vec![EdgeId(1)]).unwrap(),
            ),
            Compaction::new(6, vec![DativeBondId(1)]).unwrap(),
            Compaction::new(7, vec![AromaticSystemId(1)]).unwrap(),
            Compaction::new(8, vec![MulticenterBondId(1)]).unwrap(),
            Compaction::new(9, vec![NoncovalentBondId(1)]).unwrap(),
            Compaction::new(10, vec![StereoAtomId(1)]).unwrap(),
            Compaction::new(11, vec![StereoBondId(1)]).unwrap(),
        );
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(1)), (AtomId(2), AtomId(0))], 3, 3)
                .unwrap(),
            Correspondence::new(vec![(BondId(0), BondId(1)), (BondId(2), BondId(0))], 3, 4)
                .unwrap(),
            Correspondence::new(
                vec![
                    (DativeBondId(0), DativeBondId(1)),
                    (DativeBondId(2), DativeBondId(0)),
                ],
                3,
                5,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (AromaticSystemId(0), AromaticSystemId(1)),
                    (AromaticSystemId(2), AromaticSystemId(0)),
                ],
                3,
                6,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(0), MulticenterBondId(1)),
                    (MulticenterBondId(2), MulticenterBondId(0)),
                ],
                3,
                7,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (NoncovalentBondId(0), NoncovalentBondId(1)),
                    (NoncovalentBondId(2), NoncovalentBondId(0)),
                ],
                3,
                8,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(1)),
                    (StereoAtomId(2), StereoAtomId(0)),
                ],
                3,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoBondId(0), StereoBondId(1)),
                    (StereoBondId(2), StereoBondId(0)),
                ],
                3,
                10,
            )
            .unwrap(),
        );
        let atoms_ptr = update_correspondence.atoms().matched_pairs().as_ptr();
        let bonds_ptr = update_correspondence.bonds().matched_pairs().as_ptr();
        let dative_bonds_ptr = update_correspondence
            .dative_bonds()
            .matched_pairs()
            .as_ptr();
        let aromatic_systems_ptr = update_correspondence
            .aromatic_systems()
            .matched_pairs()
            .as_ptr();
        let multicenter_bonds_ptr = update_correspondence
            .multicenter_bonds()
            .matched_pairs()
            .as_ptr();
        let noncovalent_bonds_ptr = update_correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .as_ptr();
        let stereo_atoms_ptr = update_correspondence
            .stereo_atoms()
            .matched_pairs()
            .as_ptr();
        let stereo_bonds_ptr = update_correspondence
            .stereo_bonds()
            .matched_pairs()
            .as_ptr();
        let result = update_correspondence.compact_right(&compaction).unwrap();
        assert_eq!(result.atoms().matched_pairs().as_ptr(), atoms_ptr);
        assert_eq!(result.bonds().matched_pairs().as_ptr(), bonds_ptr);
        assert_eq!(
            result.dative_bonds().matched_pairs().as_ptr(),
            dative_bonds_ptr
        );
        assert_eq!(
            result.aromatic_systems().matched_pairs().as_ptr(),
            aromatic_systems_ptr
        );
        assert_eq!(
            result.multicenter_bonds().matched_pairs().as_ptr(),
            multicenter_bonds_ptr
        );
        assert_eq!(
            result.noncovalent_bonds().matched_pairs().as_ptr(),
            noncovalent_bonds_ptr
        );
        assert_eq!(
            result.stereo_atoms().matched_pairs().as_ptr(),
            stereo_atoms_ptr
        );
        assert_eq!(
            result.stereo_bonds().matched_pairs().as_ptr(),
            stereo_bonds_ptr
        );
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case::atoms(EntityKind::Atom, 4)]
    #[case::bonds(EntityKind::Bond, 5)]
    #[case::dative_bonds(EntityKind::DativeBond, 6)]
    #[case::aromatic_systems(EntityKind::AromaticSystem, 7)]
    #[case::multicenter_bonds(EntityKind::MulticenterBond, 8)]
    #[case::noncovalent_bonds(EntityKind::NoncovalentBond, 9)]
    #[case::stereo_atoms(EntityKind::StereoAtom, 10)]
    #[case::stereo_bonds(EntityKind::StereoBond, 11)]
    fn test_molecule_correspondence_compact_right_error(
        update_correspondence: MoleculeCorrespondence,
        #[case] kind: EntityKind,
        #[case] right_count: usize,
    ) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(4 + usize::from(kind == EntityKind::Atom)),
                Compaction::identity(5 + usize::from(kind == EntityKind::Bond)),
            ),
            Compaction::identity(6 + usize::from(kind == EntityKind::DativeBond)),
            Compaction::identity(7 + usize::from(kind == EntityKind::AromaticSystem)),
            Compaction::identity(8 + usize::from(kind == EntityKind::MulticenterBond)),
            Compaction::identity(9 + usize::from(kind == EntityKind::NoncovalentBond)),
            Compaction::identity(10 + usize::from(kind == EntityKind::StereoAtom)),
            Compaction::identity(11 + usize::from(kind == EntityKind::StereoBond)),
        );
        assert_eq!(
            update_correspondence.compact_right(&compaction),
            Err(MoleculeCorrespondenceComposeError {
                kind,
                source: CorrespondenceComposeError {
                    right_count,
                    next_left_count: right_count + 1
                },
            })
        );
    }

    #[rstest]
    fn test_molecule_correspondence_uncompact_right(update_correspondence: MoleculeCorrespondence) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(5, vec![NodeId(1)]).unwrap(),
                Compaction::new(6, vec![EdgeId(1)]).unwrap(),
            ),
            Compaction::new(7, vec![DativeBondId(1)]).unwrap(),
            Compaction::new(8, vec![AromaticSystemId(1)]).unwrap(),
            Compaction::new(9, vec![MulticenterBondId(1)]).unwrap(),
            Compaction::new(10, vec![NoncovalentBondId(1)]).unwrap(),
            Compaction::new(11, vec![StereoAtomId(1)]).unwrap(),
            Compaction::new(12, vec![StereoBondId(1)]).unwrap(),
        );
        let expected = MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(3)), (AtomId(2), AtomId(0))], 3, 5)
                .unwrap(),
            Correspondence::new(vec![(BondId(0), BondId(3)), (BondId(2), BondId(0))], 3, 6)
                .unwrap(),
            Correspondence::new(
                vec![
                    (DativeBondId(0), DativeBondId(3)),
                    (DativeBondId(2), DativeBondId(0)),
                ],
                3,
                7,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (AromaticSystemId(0), AromaticSystemId(3)),
                    (AromaticSystemId(2), AromaticSystemId(0)),
                ],
                3,
                8,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (MulticenterBondId(0), MulticenterBondId(3)),
                    (MulticenterBondId(2), MulticenterBondId(0)),
                ],
                3,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (NoncovalentBondId(0), NoncovalentBondId(3)),
                    (NoncovalentBondId(2), NoncovalentBondId(0)),
                ],
                3,
                10,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(3)),
                    (StereoAtomId(2), StereoAtomId(0)),
                ],
                3,
                11,
            )
            .unwrap(),
            Correspondence::new(
                vec![
                    (StereoBondId(0), StereoBondId(3)),
                    (StereoBondId(2), StereoBondId(0)),
                ],
                3,
                12,
            )
            .unwrap(),
        );
        let atoms_ptr = update_correspondence.atoms().matched_pairs().as_ptr();
        let bonds_ptr = update_correspondence.bonds().matched_pairs().as_ptr();
        let dative_bonds_ptr = update_correspondence
            .dative_bonds()
            .matched_pairs()
            .as_ptr();
        let aromatic_systems_ptr = update_correspondence
            .aromatic_systems()
            .matched_pairs()
            .as_ptr();
        let multicenter_bonds_ptr = update_correspondence
            .multicenter_bonds()
            .matched_pairs()
            .as_ptr();
        let noncovalent_bonds_ptr = update_correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .as_ptr();
        let stereo_atoms_ptr = update_correspondence
            .stereo_atoms()
            .matched_pairs()
            .as_ptr();
        let stereo_bonds_ptr = update_correspondence
            .stereo_bonds()
            .matched_pairs()
            .as_ptr();
        let result = update_correspondence.uncompact_right(&compaction).unwrap();
        assert_eq!(result.atoms().matched_pairs().as_ptr(), atoms_ptr);
        assert_eq!(result.bonds().matched_pairs().as_ptr(), bonds_ptr);
        assert_eq!(
            result.dative_bonds().matched_pairs().as_ptr(),
            dative_bonds_ptr
        );
        assert_eq!(
            result.aromatic_systems().matched_pairs().as_ptr(),
            aromatic_systems_ptr
        );
        assert_eq!(
            result.multicenter_bonds().matched_pairs().as_ptr(),
            multicenter_bonds_ptr
        );
        assert_eq!(
            result.noncovalent_bonds().matched_pairs().as_ptr(),
            noncovalent_bonds_ptr
        );
        assert_eq!(
            result.stereo_atoms().matched_pairs().as_ptr(),
            stereo_atoms_ptr
        );
        assert_eq!(
            result.stereo_bonds().matched_pairs().as_ptr(),
            stereo_bonds_ptr
        );
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case::atoms(EntityKind::Atom, 4)]
    #[case::bonds(EntityKind::Bond, 5)]
    #[case::dative_bonds(EntityKind::DativeBond, 6)]
    #[case::aromatic_systems(EntityKind::AromaticSystem, 7)]
    #[case::multicenter_bonds(EntityKind::MulticenterBond, 8)]
    #[case::noncovalent_bonds(EntityKind::NoncovalentBond, 9)]
    #[case::stereo_atoms(EntityKind::StereoAtom, 10)]
    #[case::stereo_bonds(EntityKind::StereoBond, 11)]
    fn test_molecule_correspondence_uncompact_right_error(
        update_correspondence: MoleculeCorrespondence,
        #[case] kind: EntityKind,
        #[case] right_count: usize,
    ) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(4 + usize::from(kind == EntityKind::Atom)),
                Compaction::identity(5 + usize::from(kind == EntityKind::Bond)),
            ),
            Compaction::identity(6 + usize::from(kind == EntityKind::DativeBond)),
            Compaction::identity(7 + usize::from(kind == EntityKind::AromaticSystem)),
            Compaction::identity(8 + usize::from(kind == EntityKind::MulticenterBond)),
            Compaction::identity(9 + usize::from(kind == EntityKind::NoncovalentBond)),
            Compaction::identity(10 + usize::from(kind == EntityKind::StereoAtom)),
            Compaction::identity(11 + usize::from(kind == EntityKind::StereoBond)),
        );
        assert_eq!(
            update_correspondence.uncompact_right(&compaction),
            Err(MoleculeCorrespondenceComposeError {
                kind,
                source: CorrespondenceComposeError {
                    right_count,
                    next_left_count: right_count + 1
                },
            })
        );
    }

    #[rstest]
    fn test_molecule_correspondence_extend_right_identity(
        update_correspondence: MoleculeCorrespondence,
    ) {
        for kind in [
            EntityKind::Atom,
            EntityKind::Bond,
            EntityKind::DativeBond,
            EntityKind::AromaticSystem,
            EntityKind::MulticenterBond,
            EntityKind::NoncovalentBond,
            EntityKind::StereoAtom,
            EntityKind::StereoBond,
        ] {
            assert_eq!(
                update_correspondence.clone().extend_right(kind, 0),
                update_correspondence
            );
        }
    }

    #[rstest]
    fn test_molecule_correspondence_compact_right_identity(
        update_correspondence: MoleculeCorrespondence,
    ) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(Compaction::identity(4), Compaction::identity(5)),
            Compaction::identity(6),
            Compaction::identity(7),
            Compaction::identity(8),
            Compaction::identity(9),
            Compaction::identity(10),
            Compaction::identity(11),
        );
        assert_eq!(
            update_correspondence.clone().compact_right(&compaction),
            Ok(update_correspondence)
        );
        let empty = MoleculeCorrespondence::empty();
        assert_eq!(
            empty.clone().compact_right(&MoleculeCompaction::empty()),
            Ok(empty)
        );
    }

    #[rstest]
    fn test_molecule_correspondence_uncompact_right_identity(
        update_correspondence: MoleculeCorrespondence,
    ) {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(Compaction::identity(4), Compaction::identity(5)),
            Compaction::identity(6),
            Compaction::identity(7),
            Compaction::identity(8),
            Compaction::identity(9),
            Compaction::identity(10),
            Compaction::identity(11),
        );
        assert_eq!(
            update_correspondence.clone().uncompact_right(&compaction),
            Ok(update_correspondence)
        );
        let empty = MoleculeCorrespondence::empty();
        assert_eq!(
            empty.clone().uncompact_right(&MoleculeCompaction::empty()),
            Ok(empty)
        );
    }

    #[fixture]
    fn composition_chain() -> [MoleculeCorrespondence; 3] {
        [
            MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(0), AtomId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(BondId(0), BondId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(DativeBondId(0), DativeBondId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
            MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(1), AtomId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(BondId(1), BondId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(DativeBondId(1), DativeBondId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(AromaticSystemId(1), AromaticSystemId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(MulticenterBondId(1), MulticenterBondId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(NoncovalentBondId(1), NoncovalentBondId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoAtomId(1), StereoAtomId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoBondId(1), StereoBondId(2))], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
            MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(2), AtomId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(BondId(2), BondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(DativeBondId(2), DativeBondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(AromaticSystemId(2), AromaticSystemId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(MulticenterBondId(2), MulticenterBondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(NoncovalentBondId(2), NoncovalentBondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoAtomId(2), StereoAtomId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoBondId(2), StereoBondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
        ]
    }

    #[rstest]
    fn test_molecule_correspondence_empty() {
        let expected = MoleculeCorrespondence {
            atoms: Correspondence::empty(),
            bonds: Correspondence::empty(),
            dative_bonds: Correspondence::empty(),
            aromatic_systems: Correspondence::empty(),
            multicenter_bonds: Correspondence::empty(),
            noncovalent_bonds: Correspondence::empty(),
            stereo_atoms: Correspondence::empty(),
            stereo_bonds: Correspondence::empty(),
        };

        assert_eq!(MoleculeCorrespondence::empty(), expected);
    }

    #[rstest]
    fn test_molecule_correspondence_accessors(correspondence: MoleculeCorrespondence) {
        assert_eq!(
            correspondence.atoms().matched_pairs(),
            &[(AtomId(0), AtomId(1))]
        );
        assert_eq!(
            correspondence.bonds().matched_pairs(),
            &[(BondId(0), BondId(2))]
        );
        assert_eq!(
            correspondence.dative_bonds().matched_pairs(),
            &[(DativeBondId(0), DativeBondId(3))]
        );
        assert_eq!(
            correspondence.aromatic_systems().matched_pairs(),
            &[(AromaticSystemId(0), AromaticSystemId(4))]
        );
        assert_eq!(
            correspondence.multicenter_bonds().matched_pairs(),
            &[(MulticenterBondId(0), MulticenterBondId(5))]
        );
        assert_eq!(
            correspondence.noncovalent_bonds().matched_pairs(),
            &[(NoncovalentBondId(0), NoncovalentBondId(6))]
        );
        assert_eq!(
            correspondence.stereo_atoms().matched_pairs(),
            &[(StereoAtomId(0), StereoAtomId(7))]
        );
        assert_eq!(
            correspondence.stereo_bonds().matched_pairs(),
            &[(StereoBondId(0), StereoBondId(8))]
        );
    }

    #[rstest]
    fn test_molecule_correspondence_induce() {
        // lhs C-C-C with a dative (donor 2 → acceptor 1); rhs adds a fourth atom + bond.
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        let rhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        let atoms = Correspondence::new(
            vec![
                (AtomId(0), AtomId(0)),
                (AtomId(1), AtomId(1)),
                (AtomId(2), AtomId(2)),
            ],
            3,
            4,
        )
        .expect("correspondence producer preserves partial-bijection invariants");

        let c = MoleculeCorrespondence::induce(&lhs, &rhs, atoms)
            .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            c.atoms().matched_pairs(),
            &[
                (AtomId(0), AtomId(0)),
                (AtomId(1), AtomId(1)),
                (AtomId(2), AtomId(2))
            ]
        );
        assert_eq!(
            c.bonds().matched_pairs(),
            &[(BondId(0), BondId(0)), (BondId(1), BondId(1))]
        );
        assert_eq!(c.bonds().right_unmatched(), vec![BondId(2)]);
        assert_eq!(
            c.dative_bonds().matched_pairs(),
            &[(DativeBondId(0), DativeBondId(0))]
        );
    }

    #[rstest]
    #[case::tetrahedral(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        vec![(StereoAtomId(0), StereoAtomId(0))],
    )]
    #[case::square_planar(
        StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
        StereoAtomForm::new(StereoKind::SquarePlanar, 1u32),
        vec![(StereoAtomId(0), StereoAtomId(0))],
    )]
    #[case::different_kinds(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
        vec![],
    )]
    #[case::undetermined(
        StereoAtomForm::default(),
        StereoAtomForm::default(),
        vec![(StereoAtomId(0), StereoAtomId(0))],
    )]
    #[case::one_undetermined(
        StereoAtomForm::default(),
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        vec![(StereoAtomId(0), StereoAtomId(0))],
    )]
    fn test_molecule_correspondence_induce_stereo_atoms(
        #[case] left: StereoAtomForm,
        #[case] right: StereoAtomForm,
        #[case] pairs: Vec<(StereoAtomId, StereoAtomId)>,
    ) {
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: (1..=4)
                .map(|idx| (AtomId(0), AtomId(idx), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        };
        let ligands: Vec<_> = (1..=4)
            .map(|idx| StereoLigand::new(AtomId(idx), StereoLigandKind::Atom))
            .collect();
        let lhs = Molecule::from_entries(MoleculeEntries {
            stereo_atoms: vec![(AtomId(0), ligands.clone(), left)],
            ..entries.clone()
        });
        let rhs = Molecule::from_entries(MoleculeEntries {
            stereo_atoms: vec![(AtomId(0), ligands, right)],
            ..entries
        });
        let expected = MoleculeCorrespondence::new(
            Correspondence::identity(5),
            Correspondence::identity(4),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::new(pairs, 1, 1).unwrap(),
            Correspondence::empty(),
        );
        assert_eq!(
            MoleculeCorrespondence::induce(&lhs, &rhs, Correspondence::identity(5)),
            Some(expected.clone())
        );
        assert_eq!(
            MoleculeCorrespondence::induce(&rhs, &lhs, Correspondence::identity(5)),
            Some(expected.reverse())
        );
    }

    #[rstest]
    #[case::cis_trans(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        StereoBondForm::new(StereoKind::CisTrans, 1u32),
        vec![(StereoBondId(0), StereoBondId(0))],
    )]
    #[case::axial(
        StereoBondForm::new(StereoKind::Axial, 0u32),
        StereoBondForm::new(StereoKind::Axial, 1u32),
        vec![(StereoBondId(0), StereoBondId(0))],
    )]
    #[case::different_kinds(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        StereoBondForm::new(StereoKind::Axial, 0u32),
        vec![],
    )]
    #[case::undetermined(
        StereoBondForm::default(),
        StereoBondForm::default(),
        vec![(StereoBondId(0), StereoBondId(0))],
    )]
    #[case::one_undetermined(
        StereoBondForm::default(),
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        vec![(StereoBondId(0), StereoBondId(0))],
    )]
    fn test_molecule_correspondence_induce_stereo_bonds(
        #[case] left: StereoBondForm,
        #[case] right: StereoBondForm,
        #[case] pairs: Vec<(StereoBondId, StereoBondId)>,
    ) {
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
                .into_iter()
                .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        };
        let ligands: Vec<_> = (2..=5)
            .map(|idx| StereoLigand::new(AtomId(idx), StereoLigandKind::Atom))
            .collect();
        let lhs = Molecule::from_entries(MoleculeEntries {
            stereo_bonds: vec![(BondId(0), ligands.clone(), left)],
            ..entries.clone()
        });
        let rhs = Molecule::from_entries(MoleculeEntries {
            stereo_bonds: vec![(BondId(0), ligands, right)],
            ..entries
        });
        let expected = MoleculeCorrespondence::new(
            Correspondence::identity(6),
            Correspondence::identity(5),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::new(pairs, 1, 1).unwrap(),
        );
        assert_eq!(
            MoleculeCorrespondence::induce(&lhs, &rhs, Correspondence::identity(6)),
            Some(expected.clone())
        );
        assert_eq!(
            MoleculeCorrespondence::induce(&rhs, &lhs, Correspondence::identity(6)),
            Some(expected.reverse())
        );
    }

    #[rstest]
    #[case::left_count(1, 2, 0, 2)]
    #[case::right_count(1, 2, 1, 1)]
    fn test_molecule_correspondence_induce_dimension_error(
        #[case] left_atom_count: usize,
        #[case] right_atom_count: usize,
        #[case] declared_left_count: usize,
        #[case] declared_right_count: usize,
    ) {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); left_atom_count],
            ..Default::default()
        });
        let rhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); right_atom_count],
            ..Default::default()
        });
        let atoms = Correspondence::new(Vec::new(), declared_left_count, declared_right_count)
            .expect("the empty correspondence is a partial bijection");

        assert_eq!(MoleculeCorrespondence::induce(&lhs, &rhs, atoms), None);
    }

    #[rstest]
    fn test_molecule_correspondence_compose(composition_chain: [MoleculeCorrespondence; 3]) {
        let [left, right, _] = composition_chain;

        assert_eq!(
            left.compose(&right).unwrap(),
            MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(0), AtomId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(BondId(0), BondId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(DativeBondId(0), DativeBondId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoBondId(0), StereoBondId(2))], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ),
        );
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative(EntityKind::DativeBond)]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::multicenter(EntityKind::MulticenterBond)]
    #[case::noncovalent(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_correspondence_compose_error(
        composition_chain: [MoleculeCorrespondence; 3],
        #[case] kind: EntityKind,
    ) {
        let [left, mut right, _] = composition_chain;
        match kind {
            EntityKind::Atom => right.atoms = Correspondence::new(vec![], 5, 3).unwrap(),
            EntityKind::Bond => right.bonds = Correspondence::new(vec![], 5, 3).unwrap(),
            EntityKind::DativeBond => {
                right.dative_bonds = Correspondence::new(vec![], 5, 3).unwrap()
            }
            EntityKind::AromaticSystem => {
                right.aromatic_systems = Correspondence::new(vec![], 5, 3).unwrap()
            }
            EntityKind::MulticenterBond => {
                right.multicenter_bonds = Correspondence::new(vec![], 5, 3).unwrap()
            }
            EntityKind::NoncovalentBond => {
                right.noncovalent_bonds = Correspondence::new(vec![], 5, 3).unwrap()
            }
            EntityKind::StereoAtom => {
                right.stereo_atoms = Correspondence::new(vec![], 5, 3).unwrap()
            }
            EntityKind::StereoBond => {
                right.stereo_bonds = Correspondence::new(vec![], 5, 3).unwrap()
            }
        }
        let expected = MoleculeCorrespondenceComposeError {
            kind,
            source: CorrespondenceComposeError {
                right_count: 2,
                next_left_count: 5,
            },
        };
        assert_eq!(left.compose(&right), Err(expected.clone()));
        assert_eq!(
            MoleculeCorrespondence::compose_all([left, right]),
            Err(expected)
        );
    }

    #[rstest]
    #[case::empty(0)]
    #[case::singleton(1)]
    #[case::multiple(3)]
    fn test_molecule_correspondence_compose_all(
        composition_chain: [MoleculeCorrespondence; 3],
        #[case] count: usize,
    ) {
        let expected = match count {
            0 => None,
            1 => Some(composition_chain[0].clone()),
            3 => Some(MoleculeCorrespondence::new(
                Correspondence::new(vec![(AtomId(0), AtomId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(BondId(0), BondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(DativeBondId(0), DativeBondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                Correspondence::new(vec![(StereoBondId(0), StereoBondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            )),
            _ => unreachable!(),
        };

        assert_eq!(
            MoleculeCorrespondence::compose_all(composition_chain.into_iter().take(count)).unwrap(),
            expected,
        );
    }

    #[rstest]
    fn test_molecule_correspondence_reverse(correspondence: MoleculeCorrespondence) {
        let reversed = correspondence.reverse();
        assert_eq!(reversed.atoms().matched_pairs(), &[(AtomId(1), AtomId(0))]);
        assert_eq!(reversed.bonds().matched_pairs(), &[(BondId(2), BondId(0))]);
        assert_eq!(
            reversed.dative_bonds().matched_pairs(),
            &[(DativeBondId(3), DativeBondId(0))]
        );
        assert_eq!(
            reversed.aromatic_systems().matched_pairs(),
            &[(AromaticSystemId(4), AromaticSystemId(0))]
        );
        assert_eq!(
            reversed.multicenter_bonds().matched_pairs(),
            &[(MulticenterBondId(5), MulticenterBondId(0))]
        );
        assert_eq!(
            reversed.noncovalent_bonds().matched_pairs(),
            &[(NoncovalentBondId(6), NoncovalentBondId(0))]
        );
        assert_eq!(
            reversed.stereo_atoms().matched_pairs(),
            &[(StereoAtomId(7), StereoAtomId(0))]
        );
        assert_eq!(
            reversed.stereo_bonds().matched_pairs(),
            &[(StereoBondId(8), StereoBondId(0))]
        );
        // Counts swap too: rhs atom 0 becomes an unmatched lhs atom.
        assert_eq!(reversed.atoms().left_unmatched(), vec![AtomId(0)]);
    }

    #[rstest]
    #[case::atom_matched(Entity::Atom(AtomId(0)), Some(Entity::Atom(AtomId(1))))]
    #[case::atom_unmatched(Entity::Atom(AtomId(1)), None)]
    #[case::bond_matched(Entity::Bond(BondId(0)), Some(Entity::Bond(BondId(2))))]
    #[case::bond_unmatched(Entity::Bond(BondId(1)), None)]
    #[case::dative_bond_matched(
        Entity::DativeBond(DativeBondId(0)),
        Some(Entity::DativeBond(DativeBondId(3)))
    )]
    #[case::dative_bond_unmatched(Entity::DativeBond(DativeBondId(1)), None)]
    #[case::aromatic_system_matched(
        Entity::AromaticSystem(AromaticSystemId(0)),
        Some(Entity::AromaticSystem(AromaticSystemId(4)))
    )]
    #[case::aromatic_system_unmatched(Entity::AromaticSystem(AromaticSystemId(1)), None)]
    #[case::multicenter_bond_matched(
        Entity::MulticenterBond(MulticenterBondId(0)),
        Some(Entity::MulticenterBond(MulticenterBondId(5)))
    )]
    #[case::multicenter_bond_unmatched(Entity::MulticenterBond(MulticenterBondId(1)), None)]
    #[case::noncovalent_bond_matched(
        Entity::NoncovalentBond(NoncovalentBondId(0)),
        Some(Entity::NoncovalentBond(NoncovalentBondId(6)))
    )]
    #[case::noncovalent_bond_unmatched(Entity::NoncovalentBond(NoncovalentBondId(1)), None)]
    #[case::stereo_atom_matched(
        Entity::StereoAtom(StereoAtomId(0)),
        Some(Entity::StereoAtom(StereoAtomId(7)))
    )]
    #[case::stereo_atom_unmatched(Entity::StereoAtom(StereoAtomId(1)), None)]
    #[case::stereo_bond_matched(
        Entity::StereoBond(StereoBondId(0)),
        Some(Entity::StereoBond(StereoBondId(8)))
    )]
    #[case::stereo_bond_unmatched(Entity::StereoBond(StereoBondId(1)), None)]
    fn test_molecule_correspondence_right_of(
        correspondence: MoleculeCorrespondence,
        #[case] left: Entity,
        #[case] expected: Option<Entity>,
    ) {
        assert_eq!(correspondence.right_of(left), expected);
    }

    #[rstest]
    #[case::atom_matched(Entity::Atom(AtomId(1)), Some(Entity::Atom(AtomId(0))))]
    #[case::atom_unmatched(Entity::Atom(AtomId(0)), None)]
    #[case::bond_matched(Entity::Bond(BondId(2)), Some(Entity::Bond(BondId(0))))]
    #[case::bond_unmatched(Entity::Bond(BondId(1)), None)]
    #[case::dative_bond_matched(
        Entity::DativeBond(DativeBondId(3)),
        Some(Entity::DativeBond(DativeBondId(0)))
    )]
    #[case::dative_bond_unmatched(Entity::DativeBond(DativeBondId(1)), None)]
    #[case::aromatic_system_matched(
        Entity::AromaticSystem(AromaticSystemId(4)),
        Some(Entity::AromaticSystem(AromaticSystemId(0)))
    )]
    #[case::aromatic_system_unmatched(Entity::AromaticSystem(AromaticSystemId(1)), None)]
    #[case::multicenter_bond_matched(
        Entity::MulticenterBond(MulticenterBondId(5)),
        Some(Entity::MulticenterBond(MulticenterBondId(0)))
    )]
    #[case::multicenter_bond_unmatched(Entity::MulticenterBond(MulticenterBondId(1)), None)]
    #[case::noncovalent_bond_matched(
        Entity::NoncovalentBond(NoncovalentBondId(6)),
        Some(Entity::NoncovalentBond(NoncovalentBondId(0)))
    )]
    #[case::noncovalent_bond_unmatched(Entity::NoncovalentBond(NoncovalentBondId(1)), None)]
    #[case::stereo_atom_matched(
        Entity::StereoAtom(StereoAtomId(7)),
        Some(Entity::StereoAtom(StereoAtomId(0)))
    )]
    #[case::stereo_atom_unmatched(Entity::StereoAtom(StereoAtomId(1)), None)]
    #[case::stereo_bond_matched(
        Entity::StereoBond(StereoBondId(8)),
        Some(Entity::StereoBond(StereoBondId(0)))
    )]
    #[case::stereo_bond_unmatched(Entity::StereoBond(StereoBondId(1)), None)]
    fn test_molecule_correspondence_left_of(
        correspondence: MoleculeCorrespondence,
        #[case] right: Entity,
        #[case] expected: Option<Entity>,
    ) {
        assert_eq!(correspondence.left_of(right), expected);
    }

    #[rstest]
    fn test_molecule_correspondence_is_total() {
        let complete = MoleculeCorrespondence::new(
            Correspondence::from_images(&[AtomId(0)], 1),
            Correspondence::from_images(&[BondId(0)], 1),
            Correspondence::from_images(&[DativeBondId(0)], 1),
            Correspondence::from_images(&[AromaticSystemId(0)], 1),
            Correspondence::from_images(&[MulticenterBondId(0)], 1),
            Correspondence::from_images(&[NoncovalentBondId(0)], 1),
            Correspondence::from_images(&[StereoAtomId(0)], 1),
            Correspondence::from_images(&[StereoBondId(0)], 1),
        );
        let mut atom_unmatched = complete.clone();
        atom_unmatched.atoms = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut bond_unmatched = complete.clone();
        bond_unmatched.bonds = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut dative_unmatched = complete.clone();
        dative_unmatched.dative_bonds = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut aromatic_unmatched = complete.clone();
        aromatic_unmatched.aromatic_systems = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut multicenter_unmatched = complete.clone();
        multicenter_unmatched.multicenter_bonds = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut noncovalent_unmatched = complete.clone();
        noncovalent_unmatched.noncovalent_bonds = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut stereo_atom_unmatched = complete.clone();
        stereo_atom_unmatched.stereo_atoms = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");
        let mut stereo_bond_unmatched = complete.clone();
        stereo_bond_unmatched.stereo_bonds = Correspondence::new(Vec::new(), 1, 1)
            .expect("correspondence producer preserves partial-bijection invariants");

        assert_eq!(
            (
                complete.is_total_on_left(),
                complete.is_total_on_right(),
                complete.is_total(),
            ),
            (true, true, true),
        );
        assert_eq!(
            [
                (
                    atom_unmatched.is_total_on_left(),
                    atom_unmatched.is_total_on_right(),
                    atom_unmatched.is_total(),
                ),
                (
                    bond_unmatched.is_total_on_left(),
                    bond_unmatched.is_total_on_right(),
                    bond_unmatched.is_total(),
                ),
                (
                    dative_unmatched.is_total_on_left(),
                    dative_unmatched.is_total_on_right(),
                    dative_unmatched.is_total(),
                ),
                (
                    aromatic_unmatched.is_total_on_left(),
                    aromatic_unmatched.is_total_on_right(),
                    aromatic_unmatched.is_total(),
                ),
                (
                    multicenter_unmatched.is_total_on_left(),
                    multicenter_unmatched.is_total_on_right(),
                    multicenter_unmatched.is_total(),
                ),
                (
                    noncovalent_unmatched.is_total_on_left(),
                    noncovalent_unmatched.is_total_on_right(),
                    noncovalent_unmatched.is_total(),
                ),
                (
                    stereo_atom_unmatched.is_total_on_left(),
                    stereo_atom_unmatched.is_total_on_right(),
                    stereo_atom_unmatched.is_total(),
                ),
                (
                    stereo_bond_unmatched.is_total_on_left(),
                    stereo_bond_unmatched.is_total_on_right(),
                    stereo_bond_unmatched.is_total(),
                ),
            ],
            [(false, false, false); 8],
        );
    }
}
