//! A partial bijection between two `Molecule` id spaces, per entity family.
//!
//! The atom part is a `Correspondence<AtomId>`; the seven other families each carry a
//! `Correspondence` over their entity id. Valueless — pairing only; adding values and a direction
//! lifts it to a reaction span.

use std::collections::BTreeMap;

use umol_graph_core::Correspondence;

use super::entity::Entity;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::Molecule;
#[cfg(test)]
use super::molecule::MoleculeEntries;
use super::remap::IdRemapping;

/// A per-entity partial bijection between two molecules: atoms + bonds + the six overlay families.
/// The matched/unmatched reads of each family are those of its `Correspondence`.
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

impl MoleculeCorrespondence {
    /// A correspondence from its eight per-family correspondences (the fully materialized form).
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

    /// Derive the full per-entity correspondence between `lhs` and `rhs` from their atom
    /// correspondence. Bonds are the induced edge correspondence; each overlay's lhs entities are
    /// matched to an rhs entity by their atom constituents mapped through `atoms`. An entity whose
    /// constituents are not all matched is unmatched.
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

        let stereo_atoms = induce_family(
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

        let stereo_bonds = induce_family(
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

    /// Relational composition, per entity family: `self` (lhs↔middle) followed by `other`
    /// (middle↔rhs), yielding a lhs↔rhs correspondence.
    pub fn compose(&self, other: &MoleculeCorrespondence) -> MoleculeCorrespondence {
        MoleculeCorrespondence::new(
            self.atoms.compose(&other.atoms),
            self.bonds.compose(&other.bonds),
            self.dative_bonds.compose(&other.dative_bonds),
            self.aromatic_systems.compose(&other.aromatic_systems),
            self.multicenter_bonds.compose(&other.multicenter_bonds),
            self.noncovalent_bonds.compose(&other.noncovalent_bonds),
            self.stereo_atoms.compose(&other.stereo_atoms),
            self.stereo_bonds.compose(&other.stereo_bonds),
        )
    }

    /// Compose molecule correspondences in iteration order. Returns `None` for an empty input and
    /// the value itself for a singleton.
    pub fn compose_all(correspondences: impl IntoIterator<Item = Self>) -> Option<Self> {
        correspondences
            .into_iter()
            .reduce(|left, right| left.compose(&right))
    }

    /// The inverse correspondence (rhs↔lhs), per entity family: each family's `reverse`.
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

    /// Whether every id in all eight entity families on the left is matched.
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

    /// Whether every id in all eight entity families on the right is matched.
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

    /// Whether every id in all eight entity families is matched on both sides.
    pub fn is_total(&self) -> bool {
        self.is_total_on_left() && self.is_total_on_right()
    }

    /// Whether this correspondence actually relates `lhs` to `rhs`: every family's declared counts
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

    /// This correspondence as an [`IdRemapping`], or `None` unless every entity family is total on
    /// the left. Each left id then maps to its matched right id.
    pub fn to_remapping(&self) -> Option<IdRemapping> {
        if !self.is_total_on_left() {
            return None;
        }
        Some(IdRemapping::new(
            self.atoms.matched_pairs().iter().copied().collect(),
            self.bonds.matched_pairs().iter().copied().collect(),
            self.dative_bonds.matched_pairs().iter().copied().collect(),
            self.aromatic_systems
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            self.multicenter_bonds
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            self.noncovalent_bonds
                .matched_pairs()
                .iter()
                .copied()
                .collect(),
            self.stereo_atoms.matched_pairs().iter().copied().collect(),
            self.stereo_bonds.matched_pairs().iter().copied().collect(),
        ))
    }

    /// The atom correspondence — the spine the other families are induced from.
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

/// The bond correspondence induced by an atom correspondence: the two molecular graphs' edge
/// correspondence under `atoms`.
pub(crate) fn induced_bonds(
    left: &Molecule,
    right: &Molecule,
    atoms: &Correspondence<AtomId>,
) -> Option<Correspondence<BondId>> {
    induce_family(
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
    induce_family(
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
    induce_family(
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
    induce_family(
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
    induce_family(
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

fn induce_family<Id, Key>(
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
    use std::collections::HashMap;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;

    #[fixture]
    fn correspondence() -> MoleculeCorrespondence {
        // distinct pairs per family so a mis-wired accessor is caught.
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
            left.compose(&right),
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
            MoleculeCorrespondence::compose_all(composition_chain.into_iter().take(count)),
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

    #[rstest]
    fn test_molecule_correspondence_to_remapping() {
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&[AtomId(1), AtomId(0)], 3),
            Correspondence::from_images(&[BondId(1), BondId(0)], 3),
            Correspondence::from_images(&[DativeBondId(1), DativeBondId(0)], 3),
            Correspondence::from_images(&[AromaticSystemId(1), AromaticSystemId(0)], 3),
            Correspondence::from_images(&[MulticenterBondId(1), MulticenterBondId(0)], 3),
            Correspondence::from_images(&[NoncovalentBondId(1), NoncovalentBondId(0)], 3),
            Correspondence::from_images(&[StereoAtomId(1), StereoAtomId(0)], 3),
            Correspondence::from_images(&[StereoBondId(1), StereoBondId(0)], 3),
        );
        let expected = IdRemapping::new(
            HashMap::from([(AtomId(0), AtomId(1)), (AtomId(1), AtomId(0))]),
            HashMap::from([(BondId(0), BondId(1)), (BondId(1), BondId(0))]),
            HashMap::from([
                (DativeBondId(0), DativeBondId(1)),
                (DativeBondId(1), DativeBondId(0)),
            ]),
            HashMap::from([
                (AromaticSystemId(0), AromaticSystemId(1)),
                (AromaticSystemId(1), AromaticSystemId(0)),
            ]),
            HashMap::from([
                (MulticenterBondId(0), MulticenterBondId(1)),
                (MulticenterBondId(1), MulticenterBondId(0)),
            ]),
            HashMap::from([
                (NoncovalentBondId(0), NoncovalentBondId(1)),
                (NoncovalentBondId(1), NoncovalentBondId(0)),
            ]),
            HashMap::from([
                (StereoAtomId(0), StereoAtomId(1)),
                (StereoAtomId(1), StereoAtomId(0)),
            ]),
            HashMap::from([
                (StereoBondId(0), StereoBondId(1)),
                (StereoBondId(1), StereoBondId(0)),
            ]),
        );

        assert_eq!(
            (
                correspondence.is_total_on_left(),
                correspondence.is_total_on_right(),
                correspondence.to_remapping(),
            ),
            (true, false, Some(expected)),
        );
    }

    #[rstest]
    fn test_molecule_correspondence_to_remapping_partial(correspondence: MoleculeCorrespondence) {
        assert_eq!(correspondence.to_remapping(), None);
    }
}
