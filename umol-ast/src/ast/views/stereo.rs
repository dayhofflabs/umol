//! Stereo atom and stereo bond views.

use std::iter;
use std::ops::Index;

use umol_graph_core::{EdgeId, FixedVarBirelationSet, NodeId, Ordered, RelationId};

use super::super::ids::{AtomId, BondId, StereoAtomId, StereoBondId};
use super::super::ligand::StereoLigand;
use super::super::molecule::MoleculeAst;
use super::super::stereo::{StereoAtomAst, StereoBondAst, StereoKind};
use super::super::traits::Lattice;

type StereoAtomSet =
    FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>;
type StereoBondSet =
    FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>;

/// Namespace accessor for stereo-atom views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct StereoAtomViews<'a> {
    molecule: &'a MoleculeAst,
    stereo_atoms: &'a StereoAtomSet,
}

impl<'a> StereoAtomViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, stereo_atoms: &'a StereoAtomSet) -> Self {
        Self {
            molecule,
            stereo_atoms,
        }
    }

    pub fn count(&self) -> usize {
        self.stereo_atoms.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = StereoAtomId> {
        self.stereo_atoms.relation_ids().map(StereoAtomId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = StereoAtomView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_atoms;
        set.relation_ids().map(move |rid| StereoAtomView {
            id: StereoAtomId::from(rid),
            site: set.participants_1(rid)[0],
            ligands: set.participants_2(rid),
            ast: set.data(rid),
            molecule,
        })
    }

    pub fn contains(&self, id: StereoAtomId) -> bool {
        self.stereo_atoms.contains(RelationId::from(id))
    }

    pub fn get(&self, id: StereoAtomId) -> Option<StereoAtomView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(StereoAtomView {
            id,
            site: self.stereo_atoms.participants_1(rid)[0],
            ligands: self.stereo_atoms.participants_2(rid),
            ast: self.stereo_atoms.data(rid),
            molecule: self.molecule,
        })
    }
}

impl<'a> Index<StereoAtomId> for StereoAtomViews<'a> {
    type Output = StereoAtomAst;
    fn index(&self, id: StereoAtomId) -> &StereoAtomAst {
        self.stereo_atoms.data(RelationId::from(id))
    }
}

/// Borrowed view of a stereo atom: the site atom, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoAtomView<'a> {
    pub id: StereoAtomId,
    site: NodeId,
    ligands: &'a [StereoLigand],
    pub ast: &'a StereoAtomAst,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoAtomView<'a> {
    /// The stereo site atom.
    pub fn site(&self) -> AtomId {
        AtomId::from(self.site)
    }

    /// View of the stereo site atom.
    pub fn site_atom(&self) -> super::atom::AtomView<'a> {
        self.molecule.atom(self.site())
    }

    /// The ordered ligands occupying the site's coordination positions.
    pub fn ligands(&self) -> &'a [StereoLigand] {
        self.ligands
    }

    /// The coordination-geometry kind.
    pub fn kind(&self) -> StereoKind {
        self.ast.kind
    }

    /// Site atom followed by the ligand atoms — the relation's full atom incidence.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let site = self.site();
        let ligands = self.ligands;
        iter::once(site).chain(ligands.iter().map(|l| l.atom()))
    }

    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }
}

/// Namespace accessor for stereo-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct StereoBondViews<'a> {
    molecule: &'a MoleculeAst,
    stereo_bonds: &'a StereoBondSet,
}

impl<'a> StereoBondViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, stereo_bonds: &'a StereoBondSet) -> Self {
        Self {
            molecule,
            stereo_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.stereo_bonds.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = StereoBondId> {
        self.stereo_bonds.relation_ids().map(StereoBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = StereoBondView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        set.relation_ids().map(move |rid| StereoBondView {
            id: StereoBondId::from(rid),
            site: set.participants_1(rid)[0],
            ligands: set.participants_2(rid),
            ast: set.data(rid),
            molecule,
        })
    }

    pub fn contains(&self, id: StereoBondId) -> bool {
        self.stereo_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: StereoBondId) -> Option<StereoBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(StereoBondView {
            id,
            site: self.stereo_bonds.participants_1(rid)[0],
            ligands: self.stereo_bonds.participants_2(rid),
            ast: self.stereo_bonds.data(rid),
            molecule: self.molecule,
        })
    }
}

impl<'a> Index<StereoBondId> for StereoBondViews<'a> {
    type Output = StereoBondAst;
    fn index(&self, id: StereoBondId) -> &StereoBondAst {
        self.stereo_bonds.data(RelationId::from(id))
    }
}

/// Borrowed view of a stereo bond: the site bond, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoBondView<'a> {
    pub id: StereoBondId,
    site: EdgeId,
    ligands: &'a [StereoLigand],
    pub ast: &'a StereoBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoBondView<'a> {
    /// The stereo site bond.
    pub fn site(&self) -> BondId {
        BondId::from(self.site)
    }

    /// The ordered ligands defining the bond's configuration.
    pub fn ligands(&self) -> &'a [StereoLigand] {
        self.ligands
    }

    /// The cis/trans kind.
    pub fn kind(&self) -> StereoKind {
        self.ast.kind
    }

    /// The site bond's two atoms followed by the ligand atoms — the relation's
    /// full atom incidence.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let [a, b] = self.molecule.bond(self.site()).atom_ids();
        let ligands = self.ligands;
        [a, b].into_iter().chain(ligands.iter().map(|l| l.atom()))
    }

    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }
}

// Builder-scope view bundles for stereo elements. `ligands` is a borrow into
// builder storage so old-state checks compare without cloning; callers clone
// only what they keep (the `ast`).

pub struct StereoAtomBuilderView<'a> {
    pub id: StereoAtomId,
    pub ast: &'a StereoAtomAst,
    pub site: AtomId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoBondBuilderView<'a> {
    pub id: StereoBondId,
    pub ast: &'a StereoBondAst,
    pub site: BondId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoAtomBuilderViewMut<'a> {
    pub id: StereoAtomId,
    pub ast: &'a mut StereoAtomAst,
    pub site: AtomId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoBondBuilderViewMut<'a> {
    pub id: StereoBondId,
    pub ast: &'a mut StereoBondAst,
    pub site: BondId,
    pub ligands: &'a [StereoLigand],
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::ids::{AtomId, BondId, StereoAtomId, StereoBondId};
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::stereo::{StereoAtomAst, StereoBondAst, StereoCosetAst, StereoKind};

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 6],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(2)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_stereo_atom_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atoms().count(), 1);
    }

    #[rstest]
    fn test_stereo_atom_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_atoms().ids().collect::<Vec<_>>(),
            vec![StereoAtomId(0)],
        );
    }

    #[rstest]
    #[case::present(StereoAtomId(0), true)]
    #[case::absent(StereoAtomId(99), false)]
    fn test_stereo_atom_views_contains(
        molecule: MoleculeAst,
        #[case] id: StereoAtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_atoms().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_atom_views_get(molecule: MoleculeAst) {
        let res = molecule.stereo_atoms().get(StereoAtomId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, StereoAtomId(0));
        assert_eq!(view.site(), AtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert_eq!(
            view.ligands(),
            &[
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
        );
        assert_eq!(
            view.ast,
            &StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_atom_views_get_none(molecule: MoleculeAst) {
        let res = molecule.stereo_atoms().get(StereoAtomId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
        );
    }

    #[rstest]
    fn test_stereo_atom_view_site_atom(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0)).site_atom();
        assert_eq!(view.id, AtomId(0));
        assert_eq!(view.ast, &AtomAst::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_views_index(molecule: MoleculeAst) {
        assert_eq!(
            &molecule.stereo_atoms()[StereoAtomId(0)],
            &StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bonds().count(), 1);
    }

    #[rstest]
    fn test_stereo_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_bonds().ids().collect::<Vec<_>>(),
            vec![StereoBondId(0)],
        );
    }

    #[rstest]
    #[case::present(StereoBondId(0), true)]
    #[case::absent(StereoBondId(99), false)]
    fn test_stereo_bond_views_contains(
        molecule: MoleculeAst,
        #[case] id: StereoBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_bond_views_get(molecule: MoleculeAst) {
        let res = molecule.stereo_bonds().get(StereoBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, StereoBondId(0));
        assert_eq!(view.site(), BondId(1));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert_eq!(
            view.ast,
            &StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_bond_views_get_none(molecule: MoleculeAst) {
        let res = molecule.stereo_bonds().get(StereoBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_bond(StereoBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
        );
    }

    #[rstest]
    fn test_stereo_bond_views_index(molecule: MoleculeAst) {
        assert_eq!(
            &molecule.stereo_bonds()[StereoBondId(0)],
            &StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        );
    }
}
