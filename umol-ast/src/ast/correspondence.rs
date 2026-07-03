//! A partial bijection between two `MoleculeAst` id spaces, per entity family.
//!
//! The atom part is a node-level `Correspondence<NodeId>` (aligned with the molecular graph, so the
//! bond correspondence is its induced edge correspondence); the seven other families each carry a
//! `Correspondence` over their entity id. Valueless — pairing only; adding values and a direction
//! lifts it to a reaction span.

use umol_graph_core::{Correspondence, NodeId};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::MoleculeAst;

/// A per-entity partial bijection between two molecules: atoms + bonds + the six overlay families.
/// The mated/exposed reads of each family are those of its `Correspondence`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCorrespondence {
    atoms: Correspondence<NodeId>,
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
        atoms: Correspondence<NodeId>,
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
    /// constituents are not all mated is exposed (not matched).
    pub fn induce(lhs: &MoleculeAst, rhs: &MoleculeAst, atoms: Correspondence<NodeId>) -> Self {
        let bonds = induced_bonds(lhs, rhs, &atoms);
        let dative_bonds = induced_dative_bonds(lhs, rhs, &atoms);
        let aromatic_systems = induced_aromatic_systems(lhs, rhs, &atoms);
        let multicenter_bonds = induced_multicenter_bonds(lhs, rhs, &atoms);
        let noncovalent_bonds = induced_noncovalent_bonds(lhs, rhs, &atoms);

        let mut stereo_atom = Vec::new();
        for sp in lhs.stereo_atoms().iter() {
            let (Some(site), Some(ligands)) = (
                map_atom(&atoms, sp.site_id()),
                map_ligands(&atoms, sp.ligand_frame()),
            ) else {
                continue;
            };
            if let Some(id) = rhs.stereo_atoms().connecting_id(site, &ligands) {
                stereo_atom.push((sp.id, id));
            }
        }
        let stereo_atoms = Correspondence::new(
            stereo_atom,
            lhs.stereo_atoms().count(),
            rhs.stereo_atoms().count(),
        );

        let mut stereo_bond = Vec::new();
        for sp in lhs.stereo_bonds().iter() {
            let (Some(site), Some(ligands)) = (
                bonds.right_of(sp.site_id()),
                map_ligands(&atoms, sp.ligand_frame()),
            ) else {
                continue;
            };
            if let Some(id) = rhs.stereo_bonds().connecting_id(site, &ligands) {
                stereo_bond.push((sp.id, id));
            }
        }
        let stereo_bonds = Correspondence::new(
            stereo_bond,
            lhs.stereo_bonds().count(),
            rhs.stereo_bonds().count(),
        );

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

    /// The atom correspondence (node-level — the spine the other families are induced from).
    pub fn atoms(&self) -> &Correspondence<NodeId> {
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
    atoms: &Correspondence<NodeId>,
) -> Correspondence<BondId> {
    Correspondence::new(
        atoms
            .edge_mates(left.raw_graph(), right.raw_graph())
            .into_iter()
            .map(|(l, r)| (BondId::from(l), BondId::from(r)))
            .collect(),
        left.bonds().count(),
        right.bonds().count(),
    )
}

/// The dative-bond correspondence induced by an atom correspondence: each left dative bond whose
/// acceptor and donors are all mated pairs with the right dative bond over the same roles.
pub(crate) fn induced_dative_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<NodeId>,
) -> Correspondence<DativeBondId> {
    let mut mates = Vec::new();
    for d in left.dative_bonds().iter() {
        let (Some(acceptor), Some(donors)) = (
            map_atom(atoms, d.acceptor_id()),
            map_atoms(atoms, d.donor_ids()),
        ) else {
            continue;
        };
        if let Some(id) = right.dative_bonds().connecting_id(acceptor, &donors) {
            mates.push((d.id, id));
        }
    }
    Correspondence::new(
        mates,
        left.dative_bonds().count(),
        right.dative_bonds().count(),
    )
}

/// The aromatic-system correspondence induced by an atom correspondence: each left system whose
/// atoms are all mated pairs with the right system over the same atom set.
pub(crate) fn induced_aromatic_systems(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<NodeId>,
) -> Correspondence<AromaticSystemId> {
    let mut mates = Vec::new();
    for a in left.aromatic_systems().iter() {
        let Some(mapped) = map_atoms(atoms, a.atom_ids()) else {
            continue;
        };
        if let Some(id) = right.aromatic_systems().connecting_id(mapped) {
            mates.push((a.id, id));
        }
    }
    Correspondence::new(
        mates,
        left.aromatic_systems().count(),
        right.aromatic_systems().count(),
    )
}

/// The multicenter-bond correspondence induced by an atom correspondence: each left bond whose
/// atoms are all mated pairs with the right bond over the same atom set.
pub(crate) fn induced_multicenter_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<NodeId>,
) -> Correspondence<MulticenterBondId> {
    let mut mates = Vec::new();
    for m in left.multicenter_bonds().iter() {
        let Some(mapped) = map_atoms(atoms, m.atom_ids()) else {
            continue;
        };
        if let Some(id) = right.multicenter_bonds().connecting_id(mapped) {
            mates.push((m.id, id));
        }
    }
    Correspondence::new(
        mates,
        left.multicenter_bonds().count(),
        right.multicenter_bonds().count(),
    )
}

/// The noncovalent-bond correspondence induced by an atom correspondence: each left bond whose two
/// atoms are both mated pairs with the right bond over the same atom pair.
pub(crate) fn induced_noncovalent_bonds(
    left: &MoleculeAst,
    right: &MoleculeAst,
    atoms: &Correspondence<NodeId>,
) -> Correspondence<NoncovalentBondId> {
    let mut mates = Vec::new();
    for nc in left.noncovalent_bonds().iter() {
        let [a, b] = nc.atom_ids();
        let (Some(first), Some(second)) = (map_atom(atoms, a), map_atom(atoms, b)) else {
            continue;
        };
        if let Some(id) = right.noncovalent_bonds().connecting_id(first, second) {
            mates.push((nc.id, id));
        }
    }
    Correspondence::new(
        mates,
        left.noncovalent_bonds().count(),
        right.noncovalent_bonds().count(),
    )
}

/// The rhs partner of a lhs atom under the atom correspondence, if mated.
pub(crate) fn map_atom(atoms: &Correspondence<NodeId>, atom: AtomId) -> Option<AtomId> {
    atoms.right_of(NodeId::from(atom)).map(AtomId::from)
}

/// The rhs partners of a set of lhs atoms, or `None` if any is exposed.
fn map_atoms(
    atoms: &Correspondence<NodeId>,
    lhs: impl IntoIterator<Item = AtomId>,
) -> Option<Vec<AtomId>> {
    lhs.into_iter().map(|a| map_atom(atoms, a)).collect()
}

/// The rhs-frame ligands (each ligand's atom mapped, its kind kept), or `None` if any ligand's
/// atom is exposed.
pub(crate) fn map_ligands(
    atoms: &Correspondence<NodeId>,
    ligands: Vec<StereoLigand>,
) -> Option<Vec<StereoLigand>> {
    ligands
        .into_iter()
        .map(|l| map_atom(atoms, l.atom_id).map(|a| StereoLigand::new(a, l.kind)))
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::dative::DativeBondAst;

    #[fixture]
    fn correspondence() -> MoleculeCorrespondence {
        // distinct pairs per family so a mis-wired accessor is caught.
        MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(1))], 1, 2),
            Correspondence::new(vec![(BondId(0), BondId(2))], 1, 3),
            Correspondence::new(vec![(DativeBondId(0), DativeBondId(3))], 1, 4),
            Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(4))], 1, 5),
            Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(5))], 1, 6),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(6))], 1, 7),
            Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(7))], 1, 8),
            Correspondence::new(vec![(StereoBondId(0), StereoBondId(8))], 1, 9),
        )
    }

    #[rstest]
    fn test_molecule_correspondence_accessors(correspondence: MoleculeCorrespondence) {
        assert_eq!(correspondence.atoms().mates(), &[(NodeId(0), NodeId(1))]);
        assert_eq!(correspondence.bonds().mates(), &[(BondId(0), BondId(2))]);
        assert_eq!(
            correspondence.dative_bonds().mates(),
            &[(DativeBondId(0), DativeBondId(3))]
        );
        assert_eq!(
            correspondence.aromatic_systems().mates(),
            &[(AromaticSystemId(0), AromaticSystemId(4))]
        );
        assert_eq!(
            correspondence.multicenter_bonds().mates(),
            &[(MulticenterBondId(0), MulticenterBondId(5))]
        );
        assert_eq!(
            correspondence.noncovalent_bonds().mates(),
            &[(NoncovalentBondId(0), NoncovalentBondId(6))]
        );
        assert_eq!(
            correspondence.stereo_atoms().mates(),
            &[(StereoAtomId(0), StereoAtomId(7))]
        );
        assert_eq!(
            correspondence.stereo_bonds().mates(),
            &[(StereoBondId(0), StereoBondId(8))]
        );
    }

    #[rstest]
    fn test_molecule_correspondence_induce() {
        // lhs C-C-C with a dative (donor 2 → acceptor 1); rhs adds a fourth atom + bond.
        let lhs = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 3],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let rhs = MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 4],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let atoms = Correspondence::new(
            vec![
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2)),
            ],
            3,
            4,
        );

        let c = MoleculeCorrespondence::induce(&lhs, &rhs, atoms);

        assert_eq!(
            c.atoms().mates(),
            &[
                (NodeId(0), NodeId(0)),
                (NodeId(1), NodeId(1)),
                (NodeId(2), NodeId(2))
            ]
        );
        assert_eq!(
            c.bonds().mates(),
            &[(BondId(0), BondId(0)), (BondId(1), BondId(1))]
        );
        assert_eq!(c.bonds().right_exposed(), vec![BondId(2)]);
        assert_eq!(
            c.dative_bonds().mates(),
            &[(DativeBondId(0), DativeBondId(0))]
        );
    }

    #[rstest]
    fn test_molecule_correspondence_compose() {
        // Every family: self maps 0→1, other maps 1→2, so the composite maps 0→2.
        let ab = MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(0), NodeId(1))], 1, 2),
            Correspondence::new(vec![(BondId(0), BondId(1))], 1, 2),
            Correspondence::new(vec![(DativeBondId(0), DativeBondId(1))], 1, 2),
            Correspondence::new(vec![(AromaticSystemId(0), AromaticSystemId(1))], 1, 2),
            Correspondence::new(vec![(MulticenterBondId(0), MulticenterBondId(1))], 1, 2),
            Correspondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 1, 2),
            Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(1))], 1, 2),
            Correspondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2),
        );
        let bc = MoleculeCorrespondence::new(
            Correspondence::new(vec![(NodeId(1), NodeId(2))], 2, 3),
            Correspondence::new(vec![(BondId(1), BondId(2))], 2, 3),
            Correspondence::new(vec![(DativeBondId(1), DativeBondId(2))], 2, 3),
            Correspondence::new(vec![(AromaticSystemId(1), AromaticSystemId(2))], 2, 3),
            Correspondence::new(vec![(MulticenterBondId(1), MulticenterBondId(2))], 2, 3),
            Correspondence::new(vec![(NoncovalentBondId(1), NoncovalentBondId(2))], 2, 3),
            Correspondence::new(vec![(StereoAtomId(1), StereoAtomId(2))], 2, 3),
            Correspondence::new(vec![(StereoBondId(1), StereoBondId(2))], 2, 3),
        );

        let ac = ab.compose(&bc);

        assert_eq!(ac.atoms().mates(), &[(NodeId(0), NodeId(2))]);
        assert_eq!(ac.bonds().mates(), &[(BondId(0), BondId(2))]);
        assert_eq!(
            ac.dative_bonds().mates(),
            &[(DativeBondId(0), DativeBondId(2))]
        );
        assert_eq!(
            ac.aromatic_systems().mates(),
            &[(AromaticSystemId(0), AromaticSystemId(2))]
        );
        assert_eq!(
            ac.multicenter_bonds().mates(),
            &[(MulticenterBondId(0), MulticenterBondId(2))]
        );
        assert_eq!(
            ac.noncovalent_bonds().mates(),
            &[(NoncovalentBondId(0), NoncovalentBondId(2))]
        );
        assert_eq!(
            ac.stereo_atoms().mates(),
            &[(StereoAtomId(0), StereoAtomId(2))]
        );
        assert_eq!(
            ac.stereo_bonds().mates(),
            &[(StereoBondId(0), StereoBondId(2))]
        );
    }

    #[rstest]
    fn test_molecule_correspondence_reverse(correspondence: MoleculeCorrespondence) {
        let reversed = correspondence.reverse();
        assert_eq!(reversed.atoms().mates(), &[(NodeId(1), NodeId(0))]);
        assert_eq!(reversed.bonds().mates(), &[(BondId(2), BondId(0))]);
        assert_eq!(
            reversed.dative_bonds().mates(),
            &[(DativeBondId(3), DativeBondId(0))]
        );
        assert_eq!(
            reversed.aromatic_systems().mates(),
            &[(AromaticSystemId(4), AromaticSystemId(0))]
        );
        assert_eq!(
            reversed.multicenter_bonds().mates(),
            &[(MulticenterBondId(5), MulticenterBondId(0))]
        );
        assert_eq!(
            reversed.noncovalent_bonds().mates(),
            &[(NoncovalentBondId(6), NoncovalentBondId(0))]
        );
        assert_eq!(
            reversed.stereo_atoms().mates(),
            &[(StereoAtomId(7), StereoAtomId(0))]
        );
        assert_eq!(
            reversed.stereo_bonds().mates(),
            &[(StereoBondId(8), StereoBondId(0))]
        );
        // counts swap too: atoms went lhs_count 1 / rhs_count 2, so the new lhs id 0 is exposed.
        assert_eq!(reversed.atoms().left_exposed(), vec![NodeId(0)]);
    }
}
