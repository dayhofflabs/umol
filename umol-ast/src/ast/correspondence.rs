//! A partial bijection between two `MoleculeAst` id spaces, per entity family.
//!
//! The atom part is a `Correspondence<AtomId>`; the seven other families each carry a
//! `Correspondence` over their entity id. Valueless — pairing only; adding values and a direction
//! lifts it to a reaction span.

use umol_graph_core::{Correspondence, NodeId};

use super::entity::Entity;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::MoleculeAst;
#[cfg(test)]
use super::molecule::MoleculeParts;
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
    /// Duplicate-incidence entities violate the entity-structure validator and do not have a unique
    /// induced pairing. To keep this infallible operation from panicking on such input, the first
    /// left entity found retains the right match and later collisions remain unmatched. No semantic
    /// correctness beyond a valid partial correspondence is promised for invalid input.
    pub fn induce(lhs: &MoleculeAst, rhs: &MoleculeAst, atoms: Correspondence<AtomId>) -> Self {
        let bonds = induced_bonds(lhs, rhs, &atoms);
        let dative_bonds = induced_dative_bonds(lhs, rhs, &atoms);
        let aromatic_systems = induced_aromatic_systems(lhs, rhs, &atoms);
        let multicenter_bonds = induced_multicenter_bonds(lhs, rhs, &atoms);
        let noncovalent_bonds = induced_noncovalent_bonds(lhs, rhs, &atoms);

        let stereo_atom = retain_unique_rights(
            lhs.stereo_atoms().iter().filter_map(|stereo| {
                let (Some(site), Some(ligands)) = (
                    map_atom(&atoms, stereo.site_id()),
                    map_ligands(&atoms, stereo.ligand_frame()),
                ) else {
                    return None;
                };
                rhs.stereo_atoms()
                    .of_id(site, &ligands)
                    .map(|right| (stereo.id, right))
            }),
            rhs.stereo_atoms().count(),
            StereoAtomId::index,
        );
        let stereo_atoms = Correspondence::new(
            stereo_atom,
            lhs.stereo_atoms().count(),
            rhs.stereo_atoms().count(),
        )
        .expect("correspondence producer preserves partial-bijection invariants");

        let stereo_bond = retain_unique_rights(
            lhs.stereo_bonds().iter().filter_map(|stereo| {
                let (Some(site), Some(ligands)) = (
                    bonds.right_of(stereo.site_id()),
                    map_ligands(&atoms, stereo.ligand_frame()),
                ) else {
                    return None;
                };
                rhs.stereo_bonds()
                    .of_id(site, &ligands)
                    .map(|right| (stereo.id, right))
            }),
            rhs.stereo_bonds().count(),
            StereoBondId::index,
        );
        let stereo_bonds = Correspondence::new(
            stereo_bond,
            lhs.stereo_bonds().count(),
            rhs.stereo_bonds().count(),
        )
        .expect("correspondence producer preserves partial-bijection invariants");

        Self::new(
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        )
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

    /// Whether every id in all eight entity families is matched on both sides.
    pub fn is_total(&self) -> bool {
        self.atoms.is_total()
            && self.bonds.is_total()
            && self.dative_bonds.is_total()
            && self.aromatic_systems.is_total()
            && self.multicenter_bonds.is_total()
            && self.noncovalent_bonds.is_total()
            && self.stereo_atoms.is_total()
            && self.stereo_bonds.is_total()
    }

    /// This correspondence as an [`IdRemapping`]. Requires every entity family to be total on the
    /// left: each left id maps to its matched right id.
    pub fn to_remapping(&self) -> IdRemapping {
        debug_assert!(
            self.atoms.matched_pair_count() == self.atoms.left_count()
                && self.bonds.matched_pair_count() == self.bonds.left_count()
                && self.dative_bonds.matched_pair_count() == self.dative_bonds.left_count()
                && self.aromatic_systems.matched_pair_count() == self.aromatic_systems.left_count()
                && self.multicenter_bonds.matched_pair_count()
                    == self.multicenter_bonds.left_count()
                && self.noncovalent_bonds.matched_pair_count()
                    == self.noncovalent_bonds.left_count()
                && self.stereo_atoms.matched_pair_count() == self.stereo_atoms.left_count()
                && self.stereo_bonds.matched_pair_count() == self.stereo_bonds.left_count(),
            "to_remapping requires every entity family to be total on the left",
        );
        IdRemapping::new(
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
        )
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
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<AtomId>,
) -> Correspondence<BondId> {
    let matched_pairs = left.bonds().iter().filter_map(|bond| {
        let [first, second] = bond.atom_ids();
        let (first, second) = (atoms.right_of(first)?, atoms.right_of(second)?);
        let right_bond = right
            .raw_graph()
            .find_edge(NodeId::from(first), NodeId::from(second))?;
        Some((bond.id, BondId::from(right_bond)))
    });
    Correspondence::new(
        retain_unique_rights(matched_pairs, right.bonds().count(), BondId::index),
        left.bonds().count(),
        right.bonds().count(),
    )
    .expect("correspondence producer preserves partial-bijection invariants")
}

/// The dative-bond correspondence induced by an atom correspondence: each left dative bond whose
/// acceptor and donors are all matched with the right dative bond over the same roles.
pub(crate) fn induced_dative_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<AtomId>,
) -> Correspondence<DativeBondId> {
    let matched_pairs = retain_unique_rights(
        left.dative_bonds().iter().filter_map(|dative| {
            let (Some(acceptor), Some(donors)) = (
                map_atom(atoms, dative.acceptor_id()),
                map_atoms(atoms, dative.donor_ids()),
            ) else {
                return None;
            };
            right
                .dative_bonds()
                .of_id(acceptor, &donors)
                .map(|right| (dative.id, right))
        }),
        right.dative_bonds().count(),
        DativeBondId::index,
    );
    Correspondence::new(
        matched_pairs,
        left.dative_bonds().count(),
        right.dative_bonds().count(),
    )
    .expect("correspondence producer preserves partial-bijection invariants")
}

/// The aromatic-system correspondence induced by an atom correspondence: each left system whose
/// atoms are all matched with the right system over the same atom set.
pub(crate) fn induced_aromatic_systems(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<AtomId>,
) -> Correspondence<AromaticSystemId> {
    let matched_pairs = retain_unique_rights(
        left.aromatic_systems().iter().filter_map(|aromatic| {
            let mapped = map_atoms(atoms, aromatic.atom_ids())?;
            right
                .aromatic_systems()
                .of_id(mapped)
                .map(|right| (aromatic.id, right))
        }),
        right.aromatic_systems().count(),
        AromaticSystemId::index,
    );
    Correspondence::new(
        matched_pairs,
        left.aromatic_systems().count(),
        right.aromatic_systems().count(),
    )
    .expect("correspondence producer preserves partial-bijection invariants")
}

/// The multicenter-bond correspondence induced by an atom correspondence: each left bond whose
/// atoms are all matched with the right bond over the same atom set.
pub(crate) fn induced_multicenter_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<AtomId>,
) -> Correspondence<MulticenterBondId> {
    let matched_pairs = retain_unique_rights(
        left.multicenter_bonds().iter().filter_map(|multicenter| {
            let mapped = map_atoms(atoms, multicenter.atom_ids())?;
            right
                .multicenter_bonds()
                .of_id(mapped)
                .map(|right| (multicenter.id, right))
        }),
        right.multicenter_bonds().count(),
        MulticenterBondId::index,
    );
    Correspondence::new(
        matched_pairs,
        left.multicenter_bonds().count(),
        right.multicenter_bonds().count(),
    )
    .expect("correspondence producer preserves partial-bijection invariants")
}

/// The noncovalent-bond correspondence induced by an atom correspondence: each left bond whose two
/// atoms are both matched with the right bond over the same atom pair.
pub(crate) fn induced_noncovalent_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<AtomId>,
) -> Correspondence<NoncovalentBondId> {
    let matched_pairs = retain_unique_rights(
        left.noncovalent_bonds().iter().filter_map(|noncovalent| {
            let [first, second] = noncovalent.atom_ids();
            let (Some(first), Some(second)) = (map_atom(atoms, first), map_atom(atoms, second))
            else {
                return None;
            };
            right
                .noncovalent_bonds()
                .of_id(first, second)
                .map(|right| (noncovalent.id, right))
        }),
        right.noncovalent_bonds().count(),
        NoncovalentBondId::index,
    );
    Correspondence::new(
        matched_pairs,
        left.noncovalent_bonds().count(),
        right.noncovalent_bonds().count(),
    )
    .expect("correspondence producer preserves partial-bijection invariants")
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

/// Preserve the first match to each right entity and drop later collisions. Collisions are possible
/// only for duplicate-incidence input rejected by the entity-structure validator; this guard keeps
/// induction non-panicking without assigning semantics to invalid input.
fn retain_unique_rights<Id>(
    pairs: impl IntoIterator<Item = (Id, Id)>,
    right_count: usize,
    index: impl Fn(Id) -> usize,
) -> Vec<(Id, Id)>
where
    Id: Copy,
{
    let mut used = vec![false; right_count];
    pairs
        .into_iter()
        .filter(|&(_, right)| {
            let used = &mut used[index(right)];
            if *used {
                false
            } else {
                *used = true;
                true
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::noncovalent::NoncovalentBondAst;

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
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondAst::from_order(1))],
            ..Default::default()
        });
        let rhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondAst::from_order(1))],
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

        let c = MoleculeCorrespondence::induce(&lhs, &rhs, atoms);

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
    fn test_molecule_correspondence_induce_partial() {
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 2],
            noncovalent: vec![
                (AtomId(0), AtomId(1), NoncovalentBondAst::default()),
                (AtomId(0), AtomId(1), NoncovalentBondAst::default()),
            ],
            ..Default::default()
        });
        let rhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 2],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::default())],
            ..Default::default()
        });

        let correspondence = MoleculeCorrespondence::induce(
            &lhs,
            &rhs,
            Correspondence::from_images(&[AtomId(0), AtomId(1)], 2),
        );

        assert_eq!(
            correspondence.noncovalent_bonds(),
            &Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(0))], 2, 1,)
                .expect("the guarded induced correspondence is a partial bijection"),
        );
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

        assert!(complete.is_total());
        assert_eq!(
            [
                atom_unmatched.is_total(),
                bond_unmatched.is_total(),
                dative_unmatched.is_total(),
                aromatic_unmatched.is_total(),
                multicenter_unmatched.is_total(),
                noncovalent_unmatched.is_total(),
                stereo_atom_unmatched.is_total(),
                stereo_bond_unmatched.is_total(),
            ],
            [false; 8],
        );
    }

    #[rstest]
    fn test_molecule_correspondence_to_remapping() {
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&[AtomId(1), AtomId(0)], 2),
            Correspondence::from_images(&[BondId(1), BondId(0)], 2),
            Correspondence::from_images(&[DativeBondId(1), DativeBondId(0)], 2),
            Correspondence::from_images(&[AromaticSystemId(1), AromaticSystemId(0)], 2),
            Correspondence::from_images(&[MulticenterBondId(1), MulticenterBondId(0)], 2),
            Correspondence::from_images(&[NoncovalentBondId(1), NoncovalentBondId(0)], 2),
            Correspondence::from_images(&[StereoAtomId(1), StereoAtomId(0)], 2),
            Correspondence::from_images(&[StereoBondId(1), StereoBondId(0)], 2),
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

        assert_eq!(correspondence.to_remapping(), expected);
    }
}
