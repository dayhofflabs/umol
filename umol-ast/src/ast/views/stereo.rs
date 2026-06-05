//! Stereo atom and stereo bond views.

use std::collections::HashSet;
use std::iter;
use std::ops::Index;

use umol_graph_core::{EdgeId, FixedVarBirelationSet, NodeId, Ordered, RelationId};
use umol_perm::Permutation;

use super::super::ids::{AtomId, BondId, StereoAtomId, StereoBondId};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::molecule::MoleculeAst;
use super::super::stereo::{StereoAtomAst, StereoBondAst, StereoKind};
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::bond::BondView;
use super::ligand::StereoLigandView;
use crate::ast::{StereoAtomConstraints, StereoBondConstraints, StereoCosetAst};

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
    #[inline]
    /// The coordination-geometry kind.
    pub fn kind(&self) -> StereoKind {
        self.ast.kind
    }

    #[inline]
    /// The stereo coset.
    pub fn coset(&self) -> &'a StereoCosetAst {
        &self.ast.coset
    }

    #[inline]
    /// The stereo atom constraints.
    pub fn constraints(&self) -> &'a StereoAtomConstraints {
        &self.ast.constraints
    }

    /// ID of the stereo site atom.
    pub fn site_id(&self) -> AtomId {
        AtomId::from(self.site)
    }

    /// View of the stereo site atom.
    pub fn site(&self) -> AtomView<'a> {
        self.molecule.atom(self.site_id())
    }

    pub fn ligand_count(&self) -> usize {
        self.ligands.len()
    }

    /// The ordered ligands occupying the site's coordination positions.
    pub fn ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        let molecule = self.molecule;
        let ligands = self.ligands;
        ligands
            .iter()
            .map(move |ligand| StereoLigandView::new(*ligand, molecule))
    }

    pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
    }

    pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atom_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn atom_ligand_count(&self) -> usize {
        self.atom_ligands().count()
    }

    pub fn implicit_hydrogen_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::ImplicitHydrogen)
    }

    pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.implicit_hydrogen_ligands()
            .map(|ligand| ligand.atom_id())
    }

    pub fn implicit_hydrogen_count(&self) -> usize {
        self.implicit_hydrogen_ligands().count()
    }

    pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::LonePair)
    }

    pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.lone_pair_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn lone_pair_count(&self) -> usize {
        self.lone_pair_ligands().count()
    }

    pub fn permutation_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<Permutation> {
        permutation_for_ligands(self.ligands, ligands)
    }

    pub fn coset_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<StereoCosetAst> {
        let perm = self.permutation_for(ligands)?;
        Some(self.coset().apply_permutation(self.kind(), perm))
    }

    /// Site atom followed by the ligand atoms — the relation's full atom incidence.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let site = self.site_id();
        let ligands = self.ligands;
        iter::once(site).chain(ligands.iter().map(|l| l.atom_id))
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
    #[inline]
    /// The coordination-geometry kind.
    pub fn kind(&self) -> StereoKind {
        self.ast.kind
    }

    #[inline]
    /// The stereo coset.
    pub fn coset(&self) -> &'a StereoCosetAst {
        &self.ast.coset
    }

    #[inline]
    /// The stereo bond constraints.
    pub fn constraints(&self) -> &'a StereoBondConstraints {
        &self.ast.constraints
    }

    /// ID of the stereo site bond.
    pub fn site_id(&self) -> BondId {
        BondId::from(self.site)
    }

    /// View of the stereo site bond.
    pub fn site(&self) -> BondView<'a> {
        self.molecule.bond(self.site_id())
    }

    pub fn ligand_count(&self) -> usize {
        self.ligands.len()
    }

    /// The ordered ligands defining the bond's configuration.
    pub fn ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        let molecule = self.molecule;
        let ligands = self.ligands;
        ligands
            .iter()
            .map(move |ligand| StereoLigandView::new(*ligand, molecule))
    }

    pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
    }

    pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atom_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn atom_ligand_count(&self) -> usize {
        self.atom_ligands().count()
    }

    pub fn implicit_hydrogen_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::ImplicitHydrogen)
    }

    pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.implicit_hydrogen_ligands()
            .map(|ligand| ligand.atom_id())
    }

    pub fn implicit_hydrogen_count(&self) -> usize {
        self.implicit_hydrogen_ligands().count()
    }

    pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::LonePair)
    }

    pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.lone_pair_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn lone_pair_count(&self) -> usize {
        self.lone_pair_ligands().count()
    }

    pub fn permutation_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<Permutation> {
        permutation_for_ligands(self.ligands, ligands)
    }

    pub fn coset_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<StereoCosetAst> {
        let perm = self.permutation_for(ligands)?;
        Some(self.coset().apply_permutation(self.kind(), perm))
    }

    /// The site bond's two atoms followed by the ligand atoms — the relation's
    /// full atom incidence.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let [a, b] = self.site().atom_ids();
        let ligands = self.ligands;
        [a, b].into_iter().chain(ligands.iter().map(|l| l.atom_id))
    }

    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }
}

fn has_unique_ligands(ligands: &[StereoLigand]) -> bool {
    ligands.iter().copied().collect::<HashSet<_>>().len() == ligands.len()
}

fn permutation_for_ligands(
    current: &[StereoLigand],
    ligands: impl IntoIterator<Item = StereoLigand>,
) -> Option<Permutation> {
    let current: Vec<StereoLigand> = current.to_vec();
    let requested: Vec<StereoLigand> = ligands.into_iter().collect();
    if current.len() != requested.len()
        || !has_unique_ligands(&current)
        || !has_unique_ligands(&requested)
    {
        return None;
    }
    let current_set: HashSet<StereoLigand> = current.iter().copied().collect();
    let requested_set: HashSet<StereoLigand> = requested.iter().copied().collect();
    (current_set == requested_set).then(|| Permutation::between(&current, &requested))
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
    use umol_perm::Permutation;
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

    #[fixture]
    fn virtual_ligand_molecule() -> MoleculeAst {
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
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
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
        assert_eq!(view.site_id(), AtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert_eq!(
            view.ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(1)),
                (StereoLigandKind::Atom, AtomId(2)),
                (StereoLigandKind::Atom, AtomId(3)),
                (StereoLigandKind::Atom, AtomId(4)),
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
    fn test_stereo_atom_view_site_id(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_view_site(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0)).site();
        assert_eq!(view.id, AtomId(0));
        assert_eq!(view.ast, &AtomAst::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_coset(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_atom(StereoAtomId(0)).coset(),
            &StereoCosetAst::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_atom_view_ligand_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).ligand_count(), 4);
    }

    #[rstest]
    fn test_stereo_atom_view_ligands(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(1)),
                (StereoLigandKind::Atom, AtomId(2)),
                (StereoLigandKind::Atom, AtomId(3)),
                (StereoLigandKind::Atom, AtomId(4)),
            ],
        );
    }

    #[rstest]
    fn test_stereo_ligand_view_atom(molecule: MoleculeAst) {
        let ligand = molecule
            .stereo_atom(StereoAtomId(0))
            .ligands()
            .next()
            .unwrap();
        let atom = ligand.atom();
        assert_eq!(atom.id, AtomId(1));
        assert_eq!(atom.ast, &AtomAst::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligand_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligand_count(),
            2,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_implicit_hydrogen_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_lone_pair_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_permutation_for(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(2),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(3),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(
            view.permutation_for(ligands.clone()),
            Some(Permutation::identity(4)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(
            view.permutation_for(reordered),
            Some(Permutation::from_image(4, &[1, 0, 2, 3])),
        );
    }

    #[rstest]
    fn test_stereo_atom_view_permutation_for_none(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(2),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(3),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.permutation_for(ligands[..3].iter().copied()), None);
        assert_eq!(
            view.permutation_for([ligands[0], ligands[0], ligands[2], ligands[3]]),
            None,
        );
        assert_eq!(
            view.permutation_for([
                ligands[0],
                ligands[1],
                ligands[2],
                StereoLigand {
                    atom_id: AtomId(99),
                    kind: StereoLigandKind::Atom,
                },
            ]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_coset_for(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(2),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(3),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(
            view.coset_for(ligands.clone()),
            Some(StereoCosetAst::Lit(1)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(view.coset_for(reordered), Some(StereoCosetAst::Lit(0)));
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
        assert_eq!(view.site_id(), BondId(1));
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
    fn test_stereo_bond_view_site_id(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).site_id(), BondId(1));
    }

    #[rstest]
    fn test_stereo_bond_view_site(molecule: MoleculeAst) {
        let view = molecule.stereo_bond(StereoBondId(0)).site();
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_stereo_bond_view_coset(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_bond(StereoBondId(0)).coset(),
            &StereoCosetAst::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_bond_view_ligand_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).ligand_count(), 2);
    }

    #[rstest]
    fn test_stereo_bond_view_ligands(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_bond(StereoBondId(0))
                .ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(4)),
                (StereoLigandKind::Atom, AtomId(5)),
            ],
        );
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new(),
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligand_ids()
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new(),
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligand_count(),
            0,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_implicit_hydrogen_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_lone_pair_ligands(virtual_ligand_molecule: MoleculeAst) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(3)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(3)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_permutation_for(molecule: MoleculeAst) {
        let view = molecule.stereo_bond(StereoBondId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(5),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(
            view.permutation_for(ligands.clone()),
            Some(Permutation::identity(2)),
        );

        let reordered = vec![ligands[1], ligands[0]];
        assert_eq!(
            view.permutation_for(reordered),
            Some(Permutation::from_image(2, &[1, 0])),
        );
    }

    #[rstest]
    fn test_stereo_bond_view_permutation_for_none(molecule: MoleculeAst) {
        let view = molecule.stereo_bond(StereoBondId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(5),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.permutation_for(ligands[..1].iter().copied()), None);
        assert_eq!(view.permutation_for([ligands[0], ligands[0]]), None);
        assert_eq!(
            view.permutation_for([
                ligands[0],
                StereoLigand {
                    atom_id: AtomId(99),
                    kind: StereoLigandKind::Atom,
                },
            ]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_coset_for(molecule: MoleculeAst) {
        let view = molecule.stereo_bond(StereoBondId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(5),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(
            view.coset_for(ligands.clone()),
            Some(StereoCosetAst::Lit(1)),
        );

        let reordered = vec![ligands[1], ligands[0]];
        assert_eq!(view.coset_for(reordered), Some(StereoCosetAst::Lit(0)));
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
