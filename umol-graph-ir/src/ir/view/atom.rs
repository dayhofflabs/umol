//! Atom views: `AtomViews` namespace, `AtomView` / `AtomViewMut` AST bundles,
//! `AtomEditorView` / `AtomEditorViewMut` builder bundles.

use umol_chem::element::Element;
use umol_graph_core::NodeId;

use super::super::atom::{AtomForm, ElementForm, IsotopeMassForm};
use super::super::boolean::BooleanForm;
use super::super::constraint::AtomConstraintsForm;
use super::super::electrons::ElectronCountsForm;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId,
};
use super::super::molecule::Molecule;
use super::super::spin::UnpairedElectronsForm;
use super::super::stereo::{StereoKind, TetrahedralStereoForm};
use super::super::traits::Lattice;
use super::super::value::NumForm;
use super::aromatic::AromaticSystemView;
use super::dative::DativeBondView;
use super::multicenter::MulticenterBondView;
use super::neighbor::NeighborView;
use super::noncovalent::NoncovalentBondView;
use super::stereo::StereoAtomView;
use crate::ir::{AromaticValenceForm, AtomConstraintForm, MulticenterValenceForm};

/// Namespace accessor for atom views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    molecule: &'a Molecule,
    atoms: &'a [AtomForm],
}

impl<'a> AtomViews<'a> {
    pub(crate) fn new(molecule: &'a Molecule, atoms: &'a [AtomForm]) -> Self {
        Self { molecule, atoms }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_graph().node_count()
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AtomId> {
        self.molecule.raw_graph().node_ids().map(AtomId::from)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = AtomView<'a>> {
        let molecule = self.molecule;
        let atoms = self.atoms;
        let graph = molecule.raw_graph();
        graph.node_ids().map(move |id| AtomView {
            id: AtomId::from(id),
            ast: &atoms[id.index()],
            molecule,
        })
    }

    pub fn contains(&self, id: AtomId) -> bool {
        self.molecule.raw_graph().contains_node(NodeId::from(id))
    }

    pub fn get(&self, id: AtomId) -> Option<AtomView<'a>> {
        if !self.contains(id) {
            return None;
        }
        Some(AtomView {
            id,
            ast: &self.atoms[id.index()],
            molecule: self.molecule,
        })
    }
}

/// Borrowed view of an atom: index, underlying `AtomForm`, and the parent
/// `Molecule` for cross-relation chemistry methods.
///
/// Chemistry methods come in pairs: the topology-derived value (summed from
/// incident bonds / dative bonds / aromatic system / multicenter bonds) and
/// the matching local-constraint value carried in `data.constraints`. The
/// validator cross-checks the two when both are ground.
#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomForm,
    molecule: &'a Molecule,
}

impl<'a> AtomView<'a> {
    #[inline]
    pub fn element(&self) -> &'a ElementForm {
        &self.ast.element
    }

    #[inline]
    pub fn isotope_mass(&self) -> &'a IsotopeMassForm {
        &self.ast.isotope_mass
    }

    #[inline]
    pub fn charge(&self) -> &'a NumForm {
        &self.ast.charge
    }

    #[inline]
    pub fn implicit_hydrogens(&self) -> &'a NumForm {
        &self.ast.implicit_hydrogens
    }

    #[inline]
    pub fn lone_pairs(&self) -> &'a NumForm {
        &self.ast.lone_pairs
    }

    #[inline]
    pub fn unpaired_electrons(&self) -> &'a UnpairedElectronsForm {
        &self.ast.unpaired_electrons
    }

    #[inline]
    pub fn constraints(&self) -> &'a AtomConstraintsForm {
        &self.ast.constraints
    }

    /// Iterator over incident bonds and their neighbor atoms. Equivalent to
    /// `self.molecule.neighbors(self.id)` but exposed on the view so closures
    /// that take `&AtomView` (e.g. perception electron-counting) can inspect
    /// bonds without reaching back to the molecule.
    /// Incident neighbors, ordered by ascending neighbor atom id.
    pub fn neighbors(&self) -> impl ExactSizeIterator<Item = NeighborView<'a>> {
        self.molecule.neighbors(self.id)
    }

    /// Ids of incident bonds, in iteration order of `neighbors`.
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = BondId> + 'a {
        self.molecule.neighbors(self.id).map(|n| n.bond_id())
    }

    /// Localized valence: sum of incident `Bond.order` values. Returns
    /// `NumForm::Lit(n)` when every incident bond order is `Lit`; collapses
    /// to `Undetermined` if any bond order is non-`Lit`.
    pub fn valence(&self) -> NumForm {
        self.neighbors()
            .map(|n| n.bond().ast.order.clone())
            .fold(NumForm::Lit(0), |acc, order| acc + order)
    }

    /// Sum of `order` over incident binary dative bonds where this atom is the
    /// sole donor. Multi-donor entries are currently skipped: their per-atom
    /// projection is a stub pending the dative versus coordination/haptic
    /// entity split in discussion doc 117. Returns `NumForm::Lit(0)` when
    /// this atom donates to no single-donor dative bonds; collapses to
    /// `Undetermined` if any contributing dative's `order` is non-`Lit`.
    pub fn donated_pairs(&self) -> NumForm {
        let mut sum = NumForm::Lit(0);
        for view in self.dative_bonds() {
            let donor_ids: Vec<AtomId> = view.donor_ids().collect();
            // TODO(doc 117): define this projection after separating binary
            // dative bonds from coordination/haptic relations.
            if donor_ids.len() != 1 || donor_ids[0] != self.id {
                continue;
            }
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Sum of `order` over incident dative bonds where this atom is the
    /// acceptor. The contribution from multi-donor entries is provisional and
    /// must not define validation semantics before the dative versus
    /// coordination/haptic entity split in discussion doc 117. Returns
    /// `NumForm::Lit(0)` when this atom is not an acceptor; collapses to
    /// `Undetermined` if any contributing dative's `order` is non-`Lit`.
    pub fn accepted_pairs(&self) -> NumForm {
        let mut sum = NumForm::Lit(0);
        for view in self.dative_bonds() {
            if view.acceptor_id() != self.id {
                continue;
            }
            // TODO(doc 117): define the multi-donor acceptor projection after
            // separating binary dative bonds from coordination/haptic relations.
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Electron contribution from the aromatic system this atom belongs to.
    /// `NumForm::Lit(0)` if the atom is not in any aromatic system;
    /// `Undetermined` if the system's per-atom electron count is non-`Lit`.
    pub fn aromatic_valence(&self) -> NumForm {
        let Some(sys) = self.aromatic_system() else {
            return NumForm::Lit(0);
        };
        let Some(pos) = sys.atom_ids().position(|a| a == self.id) else {
            return NumForm::Undetermined;
        };
        match &sys.ast.electrons {
            ElectronCountsForm::Lit(counts) => counts
                .get(pos)
                .map(|&n| NumForm::Lit(n))
                .unwrap_or(NumForm::Undetermined),
            ElectronCountsForm::Undetermined => NumForm::Undetermined,
        }
    }

    /// Electrons gained from aromatic system this atom belongs to.
    pub fn aromatic_covalence(&self) -> NumForm {
        match self.aromatic_valence() {
            NumForm::Lit(1) => NumForm::Lit(1),
            NumForm::Lit(_) => NumForm::Lit(0),
            _ => NumForm::Undetermined,
        }
    }

    /// Count of multicenter co-participants across all incident multicenter
    /// bonds. Per the no-overlap structural rule these are not localized-
    /// bond neighbors. Always `Lit`.
    pub fn multicenter_degree(&self) -> NumForm {
        let count: usize = self
            .multicenter_bonds()
            .map(|mc| mc.atom_count().saturating_sub(1))
            .sum();
        NumForm::Lit(count as i64)
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `NumForm::Lit(0)` when not in any multicenter bond; collapses to
    /// `Undetermined` if any contribution is non-`Lit`.
    pub fn multicenter_valence(&self) -> NumForm {
        let mut sum = NumForm::Lit(0);
        for view in self.multicenter_bonds() {
            let Some(pos) = view.atom_ids().position(|a| a == self.id) else {
                return NumForm::Undetermined;
            };
            let term = match &view.ast.electrons {
                ElectronCountsForm::Lit(counts) => counts
                    .get(pos)
                    .map(|&n| NumForm::Lit(n))
                    .unwrap_or(NumForm::Undetermined),
                ElectronCountsForm::Undetermined => NumForm::Undetermined,
            };
            sum = sum + term;
        }
        sum
    }

    /// Count of incident localized bonds, each weighted 1. Always `Lit`.
    pub fn degree(&self) -> NumForm {
        NumForm::Lit(self.neighbors().count() as i64)
    }

    /// `degree` + `implicit_hydrogens` + `multicenter_degree`. Collapses to
    /// `Undetermined` if any term is non-`Lit`.
    pub fn total_degree(&self) -> NumForm {
        self.degree() + self.implicit_hydrogens() + self.multicenter_degree()
    }

    /// Count of incident localized bonds whose neighbor is not a literal
    /// hydrogen atom (Element::H). Always `Lit`; non-`Lit` neighbor
    /// elements count as heavy (i.e., not filtered out).
    pub fn heavy_atom_degree(&self) -> NumForm {
        let count = self
            .neighbors()
            .filter(|n| !matches!(n.atom().element(), ElementForm::Lit(Element::H)))
            .count();
        NumForm::Lit(count as i64)
    }

    /// `valence` over incident bonds whose neighbor is not a literal
    /// hydrogen. Collapses to `Undetermined` if any contributing bond order
    /// is non-`Lit`.
    pub fn heavy_atom_valence(&self) -> NumForm {
        self.neighbors()
            .filter(|n| !matches!(n.atom().element(), ElementForm::Lit(Element::H)))
            .map(|n| n.bond().order().clone())
            .fold(NumForm::Lit(0), |acc, order| acc + order)
    }

    /// Explicit hydrogens (incident neighbors with `Element::H`) plus
    /// `implicit_hydrogens`. Collapses to `Undetermined` if `implicit_hydrogens`
    /// is non-`Lit` (including `Normal`).
    pub fn total_hydrogens(&self) -> NumForm {
        let explicit = self
            .neighbors()
            .filter(|n| matches!(n.atom().element(), ElementForm::Lit(Element::H)))
            .count() as i64;
        NumForm::Lit(explicit) + self.implicit_hydrogens()
    }

    /// Full electron-sharing sum at this atom:
    /// `valence + implicit_hydrogens + aromatic_valence + multicenter_valence`.
    /// Diverges from SMARTS `v<n>` for aromatic lone-pair donors (pyrrole N,
    /// furan O) which contribute the donated pair via `aromatic_valence`.
    pub fn total_valence(&self) -> NumForm {
        self.valence()
            + self.implicit_hydrogens()
            + self.aromatic_valence()
            + self.multicenter_valence()
    }

    /// Covalence, count of electrons gained by atom from electron sharing.
    /// `valence + implicit_hydrogens + aromatic_covalence`.
    pub fn covalence(&self) -> NumForm {
        self.valence() + self.implicit_hydrogens() + self.aromatic_covalence()
    }

    pub fn is_in_dative_bond(&self) -> bool {
        self.molecule.dative_bonds().has_incident(self.id)
    }

    pub fn dative_bonds(&self) -> impl ExactSizeIterator<Item = DativeBondView<'a>> + 'a {
        self.molecule.dative_bonds().incident(self.id)
    }

    pub fn dative_bond_ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> + 'a {
        self.molecule.dative_bonds().incident_ids(self.id)
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.molecule.aromatic_systems().has_incident(self.id)
    }

    /// The aromatic system containing this atom, if any.
    pub fn aromatic_system(&self) -> Option<AromaticSystemView<'a>> {
        self.aromatic_system_id()
            .map(|id| self.molecule.aromatic_system(id))
    }

    pub fn aromatic_system_id(&self) -> Option<AromaticSystemId> {
        self.molecule
            .aromatic_systems()
            .incident_ids(self.id)
            .next()
    }

    pub fn is_in_multicenter_bond(&self) -> bool {
        self.molecule.multicenter_bonds().has_incident(self.id)
    }

    pub fn multicenter_bonds(&self) -> impl ExactSizeIterator<Item = MulticenterBondView<'a>> + 'a {
        self.molecule.multicenter_bonds().incident(self.id)
    }

    pub fn multicenter_bond_ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> + 'a {
        self.molecule.multicenter_bonds().incident_ids(self.id)
    }

    pub fn is_in_noncovalent_bond(&self) -> bool {
        self.molecule.noncovalent_bonds().has_incident(self.id)
    }

    pub fn noncovalent_bonds(&self) -> impl ExactSizeIterator<Item = NoncovalentBondView<'a>> + 'a {
        self.molecule.noncovalent_bonds().incident(self.id)
    }

    pub fn noncovalent_bond_ids(&self) -> impl ExactSizeIterator<Item = NoncovalentBondId> + 'a {
        self.molecule.noncovalent_bonds().incident_ids(self.id)
    }

    pub fn is_stereo_atom(&self) -> bool {
        self.stereo_atom().is_some()
    }

    pub fn stereo_atom_id(&self) -> Option<StereoAtomId> {
        self.stereo_atom().map(|s| s.id)
    }

    /// The stereo atom sited on this atom, if any — any coordination geometry. An
    /// atom is the site of at most one stereo atom.
    pub fn stereo_atom(&self) -> Option<StereoAtomView<'a>> {
        self.molecule.stereo_atoms().at(self.id)
    }

    /// Derive topological constraints from atom properties.
    /// Derives topology constraints for an atom. With `include_missing`, an absent
    /// overlay yields its definite negative (`NotAromatic` / `NotMulticenter` /
    /// `NotStereo`, zero dative pairs) — the fully-perceived reading used by
    /// substructure matching and conformance. Without it, the overlay-based fields
    /// (aromatic / multicenter / stereo, and the dative pair counts) are emitted
    /// only when the overlay is present — the pre-resolution reading, where an
    /// absent overlay is not yet perceived rather than known-absent. `valence`
    /// comes from localized bonds and is always emitted.
    pub fn derive_constraints(&self, include_missing: bool) -> AtomConstraintsForm {
        let mut constraints = AtomConstraintsForm::new();
        constraints.set(AtomConstraintForm::valence(self.valence()));

        if self.is_in_dative_bond() || include_missing {
            constraints.set(AtomConstraintForm::donated_pairs(self.donated_pairs()));
            constraints.set(AtomConstraintForm::accepted_pairs(self.accepted_pairs()));
        }

        if self.is_in_aromatic_system() {
            constraints.set(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::aromatic(self.aromatic_valence()),
            ));
        } else if self
            .neighbors()
            .any(|n| matches!(n.bond().constraints().aromatic(), BooleanForm::Lit(true)))
        {
            constraints.set(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::aromatic(NumForm::Undetermined),
            ));
        } else if include_missing {
            constraints.set(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::NotAromatic,
            ));
        }

        if self.is_in_multicenter_bond() {
            constraints.set(AtomConstraintForm::multicenter_valence(
                MulticenterValenceForm::multicenter(self.multicenter_valence()),
            ));
        } else if include_missing {
            constraints.set(AtomConstraintForm::multicenter_valence(
                MulticenterValenceForm::NotMulticenter,
            ));
        }

        if let Some(stereo) = self
            .stereo_atom()
            .filter(|s| s.kind() == StereoKind::Tetrahedral)
        {
            constraints.set(AtomConstraintForm::tetrahedral_stereo(
                TetrahedralStereoForm::stereo(stereo.coset().clone()),
            ));
        } else if include_missing {
            constraints.set(AtomConstraintForm::tetrahedral_stereo(
                TetrahedralStereoForm::NotStereo,
            ));
        }

        constraints
    }

    /// Is atom ground
    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }

    /// Is atom undetermined
    pub fn is_undetermined(&self) -> bool {
        self.ast.is_undetermined()
    }
}

/// Mutable borrowed view of an atom.
#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomForm,
}

// Editor-scope view bundles for atoms.

pub struct AtomEditorView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomForm,
}

pub struct AtomEditorViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomForm,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::assert_exact_size_by;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::{
        AromaticValenceForm, AtomConstraintForm, AtomConstraintsForm, MulticenterValenceForm,
    };
    use crate::ir::dative::DativeBondForm;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId,
    };
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::stereo::{StereoAtomForm, StereoCoset, StereoKind, TetrahedralStereoForm};
    use crate::ir::value::NumForm;
    use crate::mol_dsl;

    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    fn test_atom_views_count(molecule: Molecule) {
        assert_eq!(molecule.atoms().count(), 4);
    }

    #[rstest]
    fn test_atom_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().atoms().ids(), vec![], |id| id);
        assert_exact_size_by(
            molecule.atoms().ids(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            |id| id,
        );
    }

    #[rstest]
    fn test_atom_views_iter(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().atoms().iter(), vec![], |view| {
            (view.id, view.ast.clone())
        });
        assert_exact_size_by(
            molecule.atoms().iter(),
            vec![
                (AtomId(0), AtomForm::from_element(Element::C)),
                (AtomId(1), AtomForm::from_element(Element::C)),
                (AtomId(2), AtomForm::from_element(Element::N)),
                (AtomId(3), AtomForm::from_element(Element::O)),
            ],
            |view| (view.id, view.ast.clone()),
        );
    }

    #[rstest]
    #[case::present(AtomId(2), true)]
    #[case::absent(AtomId(999), false)]
    fn test_atom_views_contains(molecule: Molecule, #[case] id: AtomId, #[case] expected: bool) {
        assert_eq!(molecule.atoms().contains(id), expected);
    }

    #[rstest]
    fn test_atom_views_get(molecule: Molecule) {
        let res = molecule.atoms().get(AtomId(2));
        assert!(res.is_some());
        let atom = res.unwrap();
        assert_eq!(atom.id, AtomId(2));
        assert_eq!(atom.ast, &AtomForm::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_get_none(molecule: Molecule) {
        let res = molecule.atoms().get(AtomId(999));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_atom_view_neighbors(molecule: Molecule) {
        let isolated = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        assert_exact_size_by(isolated.atom(AtomId(0)).neighbors(), vec![], |neighbor| {
            neighbor.atom_id()
        });

        let view = molecule.atom(AtomId(1));
        assert_exact_size_by(
            view.neighbors(),
            vec![
                (BondId(0), AtomId(0), BondForm::from_order(1)),
                (BondId(1), AtomId(2), BondForm::from_order(2)),
            ],
            |neighbor| {
                (
                    neighbor.bond_id(),
                    neighbor.atom_id(),
                    neighbor.bond().ast.clone(),
                )
            },
        );
    }

    #[rstest]
    fn test_atom_view_bond_ids(molecule: Molecule) {
        let isolated = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        assert_exact_size_by(isolated.atom(AtomId(0)).bond_ids(), vec![], |id| id);
        assert_exact_size_by(
            molecule.atom(AtomId(1)).bond_ids(),
            vec![BondId(0), BondId(1)],
            |id| id,
        );
    }

    #[rstest]
    #[case::no_incident(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(3),
        NumForm::Lit(0),
    )]
    #[case::single(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(0),
        NumForm::Lit(1),
    )]
    #[case::three_around_center(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(1),
        NumForm::Lit(3),
    )]
    #[case::double(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(2),
        NumForm::Lit(2),
    )]
    #[case::undetermined_bond(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#),
        AtomId(0),
        NumForm::Undetermined,
    )]
    fn test_atom_view_valence(
        #[case] molecule: Molecule,
        #[case] center: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(molecule.atom(center).valence(), expected);
    }

    #[rstest]
    #[case::with_constraint(Some(AtomConstraintForm::valence(4)), Some(NumForm::Lit(4)))]
    #[case::absent(None, None)]
    fn test_atom_view_valence_constraint(
        #[case] constraint: Option<AtomConstraintForm>,
        #[case] expected: Option<NumForm>,
    ) {
        let mut atom = AtomForm::from_element(Element::C);
        if let Some(c) = constraint {
            atom.constraints.set(c);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().valence(),
            expected.as_ref()
        );
    }

    #[rstest]
    #[case::donor(AtomId(0), NumForm::Lit(1))]
    #[case::acceptor(AtomId(1), NumForm::Lit(0))]
    fn test_atom_view_donated_pairs(#[case] atom: AtomId, #[case] expected: NumForm) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        assert_eq!(molecule.atom(atom).donated_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_donated_pairs_constraint() {
        let mut atom = AtomForm::from_element(Element::N);
        atom.constraints.set(AtomConstraintForm::donated_pairs(1));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().donated_pairs(),
            Some(&NumForm::Lit(1)),
        );
    }

    #[rstest]
    #[case::donor(AtomId(0), NumForm::Lit(0))]
    #[case::acceptor(AtomId(1), NumForm::Lit(1))]
    fn test_atom_view_accepted_pairs(#[case] atom: AtomId, #[case] expected: NumForm) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        });
        assert_eq!(molecule.atom(atom).accepted_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_accepted_pairs_constraint() {
        let mut atom = AtomForm::from_element(Element::C);
        atom.constraints.set(AtomConstraintForm::accepted_pairs(2));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().accepted_pairs(),
            Some(&NumForm::Lit(2)),
        );
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_not_in_system() {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(molecule.atom(AtomId(0)).aromatic_valence(), NumForm::Lit(0));
    }

    #[rstest]
    #[case::not_in_system(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        NumForm::Lit(0),
    )]
    #[case::aromatic_one(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(
                vec![AtomId(0)],
                AromaticSystemForm::from_electrons(vec![1]),
            )],
            ..Default::default()
        }),
        NumForm::Lit(1),
    )]
    #[case::aromatic_zero(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(
                vec![AtomId(0)],
                AromaticSystemForm::from_electrons(vec![0]),
            )],
            ..Default::default()
        }),
        NumForm::Lit(0),
    )]
    #[case::aromatic_two(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(
                vec![AtomId(0)],
                AromaticSystemForm::from_electrons(vec![2]),
            )],
            ..Default::default()
        }),
        NumForm::Lit(0),
    )]
    #[case::undetermined(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::default())],
            ..Default::default()
        }),
        NumForm::Undetermined,
    )]
    fn test_atom_view_aromatic_covalence(#[case] molecule: Molecule, #[case] expected: NumForm) {
        assert_eq!(molecule.atom(AtomId(0)).aromatic_covalence(), expected);
    }

    #[rstest]
    #[case::in_system(AtomId(0), true)]
    #[case::not_in_system(AtomId(3), false)]
    fn test_atom_view_is_in_aromatic_system(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.atom(atom).is_in_aromatic_system(), expected);
    }

    #[rstest]
    #[case::participant(AtomId(0), Some(AromaticSystemId(0)))]
    #[case::not_participant(AtomId(3), None)]
    fn test_atom_view_aromatic_system(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let id = molecule.atom(atom).aromatic_system().map(|v| v.id);
        assert_eq!(id, expected);
    }

    #[rstest]
    #[case::donor(AtomId(2), vec![DativeBondId(0)])]
    #[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
    #[case::uninvolved(AtomId(0), vec![])]
    fn test_atom_view_dative_bonds(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<DativeBondId>,
    ) {
        assert_exact_size_by(
            molecule.atom(atom).dative_bonds(),
            expected.clone(),
            |view| view.id,
        );
        assert_exact_size_by(molecule.atom(atom).dative_bond_ids(), expected, |id| id);
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![MulticenterBondId(0)])]
    #[case::uninvolved(AtomId(3), vec![])]
    fn test_atom_view_multicenter_bonds(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<MulticenterBondId>,
    ) {
        assert_exact_size_by(
            molecule.atom(atom).multicenter_bonds(),
            expected.clone(),
            |view| view.id,
        );
        assert_exact_size_by(molecule.atom(atom).multicenter_bond_ids(), expected, |id| {
            id
        });
    }

    #[rstest]
    #[case::endpoint_0(AtomId(0), vec![NoncovalentBondId(0)])]
    #[case::endpoint_3(AtomId(3), vec![NoncovalentBondId(0)])]
    #[case::uninvolved(AtomId(1), vec![])]
    fn test_atom_view_noncovalent_bonds(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<NoncovalentBondId>,
    ) {
        assert_exact_size_by(
            molecule.atom(atom).noncovalent_bonds(),
            expected.clone(),
            |view| view.id,
        );
        assert_exact_size_by(molecule.atom(atom).noncovalent_bond_ids(), expected, |id| {
            id
        });
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 10],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(5), AtomId(6), BondForm::from_order(1)),
                (AtomId(5), AtomId(7), BondForm::from_order(1)),
                (AtomId(5), AtomId(8), BondForm::from_order(1)),
                (AtomId(5), AtomId(9), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
                (
                    AtomId(5),
                    vec![
                        StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(9), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(1)),
                ),
            ],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::tetrahedral_site(AtomId(0), true)]
    #[case::square_planar_site(AtomId(5), true)]
    #[case::ligand(AtomId(1), false)]
    fn test_atom_view_is_stereo_atom(
        stereo_molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(stereo_molecule.atom(atom).is_stereo_atom(), expected);
    }

    #[rstest]
    #[case::tetrahedral_site(AtomId(0), Some(StereoAtomId(0)))]
    #[case::square_planar_site(AtomId(5), Some(StereoAtomId(1)))]
    #[case::ligand(AtomId(1), None)]
    fn test_atom_view_stereo_atom_id(
        stereo_molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(stereo_molecule.atom(atom).stereo_atom_id(), expected);
    }

    #[rstest]
    fn test_atom_view_stereo_atom(stereo_molecule: Molecule) {
        // kind-generic: returns the sited stereo atom of any coordination geometry
        let tetrahedral = stereo_molecule.atom(AtomId(0)).stereo_atom().unwrap();
        assert_eq!(tetrahedral.id, StereoAtomId(0));
        assert_eq!(tetrahedral.kind(), StereoKind::Tetrahedral);
        let square_planar = stereo_molecule.atom(AtomId(5)).stereo_atom().unwrap();
        assert_eq!(square_planar.id, StereoAtomId(1));
        assert_eq!(square_planar.kind(), StereoKind::SquarePlanar);
        assert!(stereo_molecule.atom(AtomId(1)).stereo_atom().is_none());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_site(AtomId(0), AtomConstraintsForm::from_iter([
        AtomConstraintForm::valence(NumForm::Lit(4)),
        AtomConstraintForm::donated_pairs(NumForm::Lit(0)),
        AtomConstraintForm::accepted_pairs(NumForm::Lit(0)),
        AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::stereo(StereoCoset::Lit(1))),
    ]))]
    #[case::non_stereo_ligand(AtomId(1), AtomConstraintsForm::from_iter([
        AtomConstraintForm::valence(NumForm::Lit(1)),
        AtomConstraintForm::donated_pairs(NumForm::Lit(0)),
        AtomConstraintForm::accepted_pairs(NumForm::Lit(0)),
        AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
    ]))]
    #[case::square_planar_site(AtomId(5), AtomConstraintsForm::from_iter([
        AtomConstraintForm::valence(NumForm::Lit(4)),
        AtomConstraintForm::donated_pairs(NumForm::Lit(0)),
        AtomConstraintForm::accepted_pairs(NumForm::Lit(0)),
        AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
    ]))]
    fn test_atom_view_derive_constraints(
        stereo_molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: AtomConstraintsForm,
    ) {
        assert_eq!(stereo_molecule.atom(atom).derive_constraints(true), expected);
    }

    #[rstest]
    #[case::aromatic(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::default())],
            ..Default::default()
        }),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(0),
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(
                NumForm::Undetermined,
            )),
        ]),
    )]
    #[case::multicenter(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            ..Default::default()
        }),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::valence(0),
            AtomConstraintForm::multicenter_valence(MulticenterValenceForm::multicenter(
                NumForm::Undetermined,
            )),
        ]),
    )]
    fn test_atom_view_derive_constraints_partial(
        #[case] molecule: Molecule,
        #[case] expected: AtomConstraintsForm,
    ) {
        assert_eq!(molecule.atom(AtomId(0)).derive_constraints(false), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_constraint() {
        let mut atom = AtomForm::from_element(Element::C);
        atom.constraints.set(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::Aromatic(NumForm::Lit(1)),
        ));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().aromatic_valence(),
            Some(&AromaticValenceForm::Aromatic(NumForm::Lit(1))),
        );
    }

    #[rstest]
    #[case::single_bond(
        vec![(vec![AtomId(0), AtomId(1)], ElectronCountsForm::Lit(vec![2, 2]))],
        NumForm::Lit(2),
    )]
    #[case::two_bonds(
        vec![
            (vec![AtomId(0), AtomId(1)], ElectronCountsForm::Lit(vec![2, 2])),
            (vec![AtomId(0), AtomId(2)], ElectronCountsForm::Lit(vec![1, 1])),
        ],
        NumForm::Lit(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomId(0), AtomId(1)], ElectronCountsForm::Undetermined)],
        NumForm::Undetermined,
    )]
    fn test_atom_view_multicenter_valence(
        #[case] bonds: Vec<(Vec<AtomId>, ElectronCountsForm)>,
        #[case] expected: NumForm,
    ) {
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, MulticenterBondForm::new(electrons)))
            .collect();
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            multicenter,
            ..Default::default()
        });
        assert_eq!(molecule.atom(AtomId(0)).multicenter_valence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_valence_constraint() {
        let mut atom = AtomForm::from_element(Element::C);
        atom.constraints
            .set(AtomConstraintForm::multicenter_valence(
                MulticenterValenceForm::Multicenter(NumForm::Lit(2)),
            ));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().multicenter_valence(),
            Some(&MulticenterValenceForm::Multicenter(NumForm::Lit(2))),
        );
    }

    #[rstest]
    #[case::ethane_carbon(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#), AtomId(0), NumForm::Lit(1))]
    #[case::ethene_carbon(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#), AtomId(0), NumForm::Lit(1))]
    #[case::three_bonds(mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#), AtomId(0), NumForm::Lit(3))]
    fn test_atom_view_degree(
        #[case] mol: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(mol.atom(atom).degree(), expected);
    }

    #[rstest]
    #[case::lit(mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#), NumForm::Lit(4))]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C#h*"] :bonds []}"#), NumForm::Undetermined)]
    fn test_atom_view_total_degree(#[case] molecule: Molecule, #[case] expected: NumForm) {
        assert_eq!(molecule.atom(AtomId(0)).total_degree(), expected);
    }

    #[rstest]
    #[case::all_heavy(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        NumForm::Lit(2),
    )]
    #[case::one_h_neighbor(
        mol_dsl!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        NumForm::Lit(1),
    )]
    fn test_atom_view_heavy_atom_degree(
        #[case] mol: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_degree(), expected);
    }

    #[rstest]
    #[case::all_heavy(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "2"]]}"#),
        AtomId(0),
        NumForm::Lit(3),
    )]
    #[case::skips_h(
        mol_dsl!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "2"] [0 2 "1"]]}"#),
        AtomId(0),
        NumForm::Lit(2),
    )]
    fn test_atom_view_heavy_atom_valence(
        #[case] mol: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_valence(), expected);
    }

    #[rstest]
    #[case::implicit_only(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomId(0),
        NumForm::Lit(4),
    )]
    #[case::implicit_and_explicit(
        mol_dsl!(r#"{:atoms ["C#h2" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        NumForm::Lit(4),
    )]
    #[case::implicit_undetermined(
        mol_dsl!(r#"{:atoms ["C#h*"] :bonds []}"#),
        AtomId(0),
        NumForm::Undetermined,
    )]
    fn test_atom_view_total_hydrogens(
        #[case] mol: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(mol.atom(atom).total_hydrogens(), expected);
    }

    #[rstest]
    #[case::lit(mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#), NumForm::Lit(4))]
    #[case::undetemined(mol_dsl!(r#"{:atoms ["C#h*"] :bonds []}"#), NumForm::Undetermined)]
    fn test_atom_view_total_valence_sum_of_terms(
        #[case] molecule: Molecule,
        #[case] expected: NumForm,
    ) {
        assert_eq!(molecule.atom(AtomId(0)).total_valence(), expected);
    }

    #[rstest]
    #[case::ch4(mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#), NumForm::Lit(4))]
    #[case::undetermined_h(mol_dsl!(r#"{:atoms ["C#h*"] :bonds []}"#), NumForm::Undetermined)]
    fn test_atom_view_covalence_non_aromatic(
        #[case] molecule: Molecule,
        #[case] expected: NumForm,
    ) {
        assert_eq!(molecule.atom(AtomId(0)).covalence(), expected);
    }

    #[fixture]
    fn aromatic_ring() -> Molecule {
        // 3-membered C ring, each with 0 implicit H (valence 2 from two ring
        // bonds), aromatic system electrons [1, 2, 0].
        let carbon = AtomForm::from_element(Element::C).with_implicit_hydrogens(0_i64);
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![carbon.clone(), carbon.clone(), carbon],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(0), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 0]),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::standard(AtomId(0), NumForm::Lit(3))] // av=1 → +1
    #[case::donor(AtomId(1), NumForm::Lit(2))] // av=2 (donated pair) → +0
    #[case::acceptor(AtomId(2), NumForm::Lit(2))] // av=0 → +0
    fn test_atom_view_covalence_aromatic(
        aromatic_ring: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(aromatic_ring.atom(atom).covalence(), expected);
    }

    #[fixture]
    fn dative_pair() -> Molecule {
        // H₃N→BH₃: N (3 H) donates a pair to B (3 H). Localized valence plus
        // implicit hydrogens plus aromatic covalence is 3 for both; the dative
        // bond (donated on N, accepted on B) is excluded.
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N).with_implicit_hydrogens(3_i64),
                AtomForm::from_element(Element::B).with_implicit_hydrogens(3_i64),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::donor(AtomId(0), NumForm::Lit(3))] // donated pair excluded → v+h = 3
    #[case::acceptor(AtomId(1), NumForm::Lit(3))] // accepted pair excluded → v+h = 3
    fn test_atom_view_covalence_dative(
        dative_pair: Molecule,
        #[case] atom: AtomId,
        #[case] expected: NumForm,
    ) {
        assert_eq!(dative_pair.atom(atom).covalence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_degree() {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![2, 2, 2]),
            )],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).multicenter_degree(),
            NumForm::Lit(2),
        );
    }
}
