//! Atom views: `AtomViews` namespace, `AtomView` / `AtomViewMut` AST bundles,
//! `AtomBuilderView` / `AtomBuilderViewMut` builder bundles.

use std::ops::Index;

use umol_graph_core::NodeId;
use umol_shared::element::Element;

use super::super::atom::{AtomAst, ElementAst, IsotopeMassAst};
use super::super::constraint::{AtomConstraints, RingScope};
use super::super::electrons::ElectronCountsAst;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId,
};
use super::super::molecule::MoleculeAst;
use super::super::ring::{RingSet, RingView};
use super::super::spin::SpinStateAst;
use super::super::stereo::{StereoKind, TetrahedralStereoAst};
use super::super::traits::Lattice;
use super::super::value::ValueAst;
use super::aromatic::AromaticSystemView;
use super::dative::DativeBondView;
use super::multicenter::MulticenterBondView;
use super::neighbor::NeighborView;
use super::noncovalent::NoncovalentBondView;
use super::stereo::StereoAtomView;
use crate::ast::{AromaticValenceAst, AsLit, AtomConstraint, MulticenterValenceAst};

/// Namespace accessor for atom views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    molecule: &'a MoleculeAst,
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, atoms: &'a [AtomAst]) -> Self {
        Self { molecule, atoms }
    }

    pub fn count(&self) -> usize {
        self.molecule.raw_graph().node_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomId> {
        self.molecule.raw_graph().node_ids().map(AtomId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
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

impl<'a> Index<AtomId> for AtomViews<'a> {
    type Output = AtomAst;
    fn index(&self, id: AtomId) -> &AtomAst {
        &self.atoms[id.index()]
    }
}

/// Borrowed view of an atom: index, underlying `AtomAst`, and the parent
/// `MoleculeAst` for cross-relation chemistry methods.
///
/// Chemistry methods come in pairs: the topology-derived value (summed from
/// incident bonds / dative bonds / aromatic system / multicenter bonds) and
/// the matching local-constraint value carried in `data.constraints`. The
/// validator cross-checks the two when both are ground.
#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AtomView<'a> {
    #[inline]
    pub fn element(&self) -> &'a ElementAst {
        &self.ast.element
    }

    #[inline]
    pub fn isotope_mass(&self) -> &'a IsotopeMassAst {
        &self.ast.isotope_mass
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn implicit_hydrogens(&self) -> &'a ValueAst {
        &self.ast.implicit_hydrogens
    }

    #[inline]
    pub fn lone_pairs(&self) -> &'a ValueAst {
        &self.ast.lone_pairs
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a AtomConstraints {
        &self.ast.constraints
    }

    /// Iterator over incident bonds and their neighbor atoms. Equivalent to
    /// `self.molecule.neighbors(self.id)` but exposed on the view so closures
    /// that take `&AtomView` (e.g. perception electron-counting) can inspect
    /// bonds without reaching back to the molecule.
    /// Incident neighbors, ordered by ascending neighbor atom id.
    pub fn neighbors(&self) -> impl Iterator<Item = NeighborView<'a>> {
        self.molecule.neighbors(self.id)
    }

    /// Ids of incident bonds, in iteration order of `neighbors`.
    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule.neighbors(self.id).map(|n| n.bond_id())
    }

    /// Localized valence: sum of incident `Bond.order` values. Returns
    /// `ValueAst::Lit(n)` when every incident bond order is `Lit`; collapses
    /// to `Undetermined` if any bond order is non-`Lit`.
    pub fn valence(&self) -> ValueAst {
        self.neighbors()
            .map(|n| n.bond().ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Sum of `order` over incident dative bonds where this atom is the sole
    /// donor (multi-donor datives contribute nothing per individual donor —
    /// the donated pair is collective and has no well-defined per-atom
    /// share). Returns `ValueAst::Lit(0)` when this atom donates to no
    /// single-donor dative bonds; collapses to `Undetermined` if any
    /// contributing dative's `order` is non-`Lit`.
    pub fn donated_pairs(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.dative_bonds() {
            let donor_ids: Vec<AtomId> = view.donor_ids().collect();
            if donor_ids.len() != 1 || donor_ids[0] != self.id {
                continue;
            }
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Sum of `order` over incident dative bonds where this atom is the
    /// acceptor. Returns `ValueAst::Lit(0)` when this atom is not an
    /// acceptor; collapses to `Undetermined` if any contributing dative's
    /// `order` is non-`Lit`.
    pub fn accepted_pairs(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.dative_bonds() {
            if view.acceptor_id() != self.id {
                continue;
            }
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Electron contribution from the aromatic system this atom belongs to.
    /// `ValueAst::Lit(0)` if the atom is not in any aromatic system;
    /// `Undetermined` if the system's per-atom electron count is non-`Lit`.
    pub fn aromatic_valence(&self) -> ValueAst {
        let Some(sys) = self.aromatic_system() else {
            return ValueAst::Lit(0);
        };
        let Some(pos) = sys.atom_ids().position(|a| a == self.id) else {
            return ValueAst::Undetermined;
        };
        match &sys.ast.electrons {
            ElectronCountsAst::Lit(counts) => counts
                .get(pos)
                .map(|&n| ValueAst::Lit(n))
                .unwrap_or(ValueAst::Undetermined),
            ElectronCountsAst::Undetermined => ValueAst::Undetermined,
        }
    }

    /// Electrons gained from aromatic system this atom belongs to.
    pub fn aromatic_increment(&self) -> ValueAst {
        match self.aromatic_valence() {
            ValueAst::Lit(1) => ValueAst::Lit(1),
            ValueAst::Lit(_) => ValueAst::Lit(0),
            _ => ValueAst::Undetermined,
        }
    }

    /// Count of multicenter co-participants across all incident multicenter
    /// bonds. Per the no-overlap structural rule these are not localized-
    /// bond neighbors. Always `Lit`.
    pub fn multicenter_degree(&self) -> ValueAst {
        let count: usize = self
            .multicenter_bonds()
            .map(|mc| mc.atom_count().saturating_sub(1))
            .sum();
        ValueAst::Lit(count as i64)
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `ValueAst::Lit(0)` when not in any multicenter bond; collapses to
    /// `Undetermined` if any contribution is non-`Lit`.
    pub fn multicenter_valence(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.multicenter_bonds() {
            let Some(pos) = view.atom_ids().position(|a| a == self.id) else {
                return ValueAst::Undetermined;
            };
            let term = match &view.ast.electrons {
                ElectronCountsAst::Lit(counts) => counts
                    .get(pos)
                    .map(|&n| ValueAst::Lit(n))
                    .unwrap_or(ValueAst::Undetermined),
                ElectronCountsAst::Undetermined => ValueAst::Undetermined,
            };
            sum = sum + term;
        }
        sum
    }

    /// Count of incident localized bonds, each weighted 1. Always `Lit`.
    pub fn degree(&self) -> ValueAst {
        ValueAst::Lit(self.neighbors().count() as i64)
    }

    /// `degree` + `implicit_hydrogens` + `multicenter_degree`. Collapses to
    /// `Undetermined` if any term is non-`Lit`.
    pub fn total_degree(&self) -> ValueAst {
        self.degree() + self.implicit_hydrogens() + self.multicenter_degree()
    }

    /// Count of incident localized bonds whose neighbor is not a literal
    /// hydrogen atom (Element::H). Always `Lit`; non-`Lit` neighbor
    /// elements count as heavy (i.e., not filtered out).
    pub fn heavy_atom_degree(&self) -> ValueAst {
        let count = self
            .neighbors()
            .filter(|n| !matches!(n.atom().element(), ElementAst::Lit(Element::H)))
            .count();
        ValueAst::Lit(count as i64)
    }

    /// `valence` over incident bonds whose neighbor is not a literal
    /// hydrogen. Collapses to `Undetermined` if any contributing bond order
    /// is non-`Lit`.
    pub fn heavy_atom_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| !matches!(n.atom().element(), ElementAst::Lit(Element::H)))
            .map(|n| n.bond().order().clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Explicit hydrogens (incident neighbors with `Element::H`) plus
    /// `implicit_hydrogens`. Collapses to `Undetermined` if `implicit_hydrogens`
    /// is non-`Lit` (including `Normal`).
    pub fn total_hydrogens(&self) -> ValueAst {
        let explicit = self
            .neighbors()
            .filter(|n| matches!(n.atom().element(), ElementAst::Lit(Element::H)))
            .count() as i64;
        ValueAst::Lit(explicit) + self.implicit_hydrogens()
    }

    /// Full electron-sharing sum at this atom:
    /// `valence + implicit_hydrogens + aromatic_valence + multicenter_valence`.
    /// Diverges from SMARTS `v<n>` for aromatic lone-pair donors (pyrrole N,
    /// furan O) which contribute the donated pair via `aromatic_valence`.
    pub fn total_valence(&self) -> ValueAst {
        self.valence()
            + self.implicit_hydrogens()
            + self.aromatic_valence()
            + self.multicenter_valence()
    }

    /// Covalence, count of electrons gained by atom from electron sharing.
    /// `valence + implicit_hydrogens + aromatic_increment`.
    pub fn covalence(&self) -> ValueAst {
        self.valence() + self.implicit_hydrogens() + self.aromatic_increment()
    }

    pub fn is_in_dative_bond(&self) -> bool {
        self.molecule.dative_bonds().has_incident(self.id)
    }

    pub fn dative_bonds(&self) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        self.molecule.dative_bonds().incident(self.id)
    }

    pub fn dative_bond_ids(&self) -> impl Iterator<Item = DativeBondId> + 'a {
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

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = MulticenterBondView<'a>> + 'a {
        self.molecule.multicenter_bonds().incident(self.id)
    }

    pub fn multicenter_bond_ids(&self) -> impl Iterator<Item = MulticenterBondId> + 'a {
        self.molecule.multicenter_bonds().incident_ids(self.id)
    }

    pub fn is_in_noncovalent_bond(&self) -> bool {
        self.molecule.noncovalent_bonds().has_incident(self.id)
    }

    pub fn noncovalent_bonds(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> + 'a {
        self.molecule.noncovalent_bonds().incident(self.id)
    }

    pub fn noncovalent_bond_ids(&self) -> impl Iterator<Item = NoncovalentBondId> + 'a {
        self.molecule.noncovalent_bonds().incident_ids(self.id)
    }

    pub fn is_in_tetrahedral_stereo(&self) -> bool {
        self.tetrahedral_stereo().is_some()
    }

    pub fn tetrahedral_stereo_id(&self) -> Option<StereoAtomId> {
        self.tetrahedral_stereo().map(|s| s.id)
    }

    /// The tetrahedral stereo atom sited on this atom, if any. An atom is the
    /// site of at most one stereo atom; the kind filter selects the tetrahedral
    /// case from the other coordination geometries that share the relation.
    pub fn tetrahedral_stereo(&self) -> Option<StereoAtomView<'a>> {
        self.molecule
            .stereo_atoms()
            .coincident(self.id)
            .filter(|s| s.kind() == StereoKind::Tetrahedral)
    }

    /// True if this atom participates in any overlay relation (aromatic
    /// system, dative bond, multicenter bond, noncovalent bond, tetrahedral
    /// stereo). Mirror of `MoleculeAst::has_overlays` scoped to a single atom;
    /// useful as a pre-mutation predicate before structural removal.
    pub fn is_in_overlays(&self) -> bool {
        self.is_in_aromatic_system()
            || self.is_in_dative_bond()
            || self.is_in_multicenter_bond()
            || self.is_in_noncovalent_bond()
            || self.is_in_tetrahedral_stereo()
    }

    /// True if this atom belongs to any ring in the molecule's canonical
    /// ring set (Vismara relevant cycles, max ring size 22). Uses the
    /// molecule's cached canonical `RingSet`.
    pub fn is_in_ring(&self) -> bool {
        self.molecule.rings().contains_atom(self.id)
    }

    /// True if this atom appears in any ring of the supplied set.
    pub fn is_in_ring_from(&self, rings: &RingSet) -> bool {
        rings.contains_atom(self.id)
    }

    /// Rings containing this atom drawn from the molecule's canonical
    /// `RingSet` (Vismara relevant cycles, max ring size 22).
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.molecule
            .rings()
            .iter()
            .filter(move |v| v.atoms().contains(&id))
    }

    /// Rings from the supplied set that contain this atom.
    pub fn rings_from<'r>(&self, rings: &'r RingSet) -> impl Iterator<Item = RingView<'r>> + 'r {
        let id = self.id;
        rings.iter().filter(move |v| v.atoms().contains(&id))
    }

    /// Derived count of canonical rings containing this atom matching
    /// `scope` (`All` = any, `Size(s)` = size `s`). Always `Lit`.
    pub fn ring_membership(&self, scope: RingScope) -> ValueAst {
        let count = match scope {
            RingScope::All => self.rings().count(),
            RingScope::Size(s) => self.rings().filter(|r| r.len() == s as usize).count(),
        };
        ValueAst::Lit(count as i64)
    }

    /// `ring_membership(RingScope::All)`.
    pub fn ring_count(&self) -> ValueAst {
        self.ring_membership(RingScope::All)
    }

    /// `ring_membership(RingScope::Size(s))`.
    pub fn ring_size_count(&self, s: u8) -> ValueAst {
        self.ring_membership(RingScope::Size(s))
    }

    /// Smallest containing canonical ring size, or `None` if not in any
    /// ring. Chemistry-classification helper, not a constraint counterpart.
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.rings().map(|r| r.len()).min()
    }

    /// Count of incident bonds that participate in any canonical ring.
    /// Always `Lit`.
    pub fn ring_degree(&self) -> ValueAst {
        let count = self.neighbors().filter(|n| n.bond().is_in_ring()).count();
        ValueAst::Lit(count as i64)
    }

    /// Sum of bond orders of incident bonds that participate in any
    /// canonical ring. Collapses to `Undetermined` if any contributing
    /// bond's `order` is non-`Lit`.
    pub fn ring_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| n.bond().is_in_ring())
            .map(|n| n.bond().order().clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Derive topological constraints from atom properties.
    pub fn derive_constraints(&self) -> AtomConstraints {
        let valence = self.valence();
        let donated_pairs = self.donated_pairs();
        let accepted_pairs = self.accepted_pairs();
        let aromatic_valence = if self.is_in_aromatic_system() {
            AromaticValenceAst::aromatic(
                self.aromatic_valence()
                    .as_lit_expect("aromatic valence should be Lit"),
            )
        } else if self.neighbors().any(|n| n.bond().constraints().aromatic()) {
            AromaticValenceAst::aromatic(ValueAst::Undetermined)
        } else {
            AromaticValenceAst::NotAromatic
        };

        let multicenter_valence = if self.is_in_multicenter_bond() {
            MulticenterValenceAst::multicenter(
                self.multicenter_valence()
                    .as_lit_expect("multicenter valence should be Lit"),
            )
        } else {
            MulticenterValenceAst::NotMulticenter
        };

        let tetrahedral_stereo = match self.tetrahedral_stereo() {
            Some(stereo) => TetrahedralStereoAst::stereo(stereo.coset().clone()),
            None => TetrahedralStereoAst::NotStereo,
        };

        AtomConstraints::from_iter([
            AtomConstraint::valence(valence),
            AtomConstraint::donated_pairs(donated_pairs),
            AtomConstraint::accepted_pairs(accepted_pairs),
            AtomConstraint::aromatic_valence(aromatic_valence),
            AtomConstraint::multicenter_valence(multicenter_valence),
            AtomConstraint::tetrahedral_stereo(tetrahedral_stereo),
        ])
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
    pub ast: &'a mut AtomAst,
}

// Builder-scope view bundles for atoms.

pub struct AtomBuilderView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomAst,
}

pub struct AtomBuilderViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomAst,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraint, AtomConstraints, Constraints, MulticenterValenceAst,
        RingScope,
    };
    use crate::ast::dative::DativeBondAst;
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId,
    };
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::ring::RingFamily;
    use crate::ast::stereo::{StereoAtomAst, TetrahedralStereoAst, StereoCosetAst, StereoKind};
    use crate::ast::value::ValueAst;
    use crate::mol;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomId(0),
                AtomId(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[fixture]
    fn ring_with_chain() -> MoleculeAst {
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C); 7],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(4), BondAst::from_order(1)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
                (AtomId(5), AtomId(0), BondAst::from_order(1)),
                (AtomId(0), AtomId(6), BondAst::from_order(1)),
            ],
        )
    }

    #[rstest]
    fn test_atom_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.atoms().count(), 4);
    }

    #[rstest]
    fn test_atom_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.atoms().ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_atom_views_iter(molecule: MoleculeAst) {
        let views = molecule.atoms();
        let collected: Vec<(AtomId, AtomAst)> =
            views.iter().map(|v| (v.id, v.ast.clone())).collect();
        assert_eq!(
            collected,
            vec![
                (AtomId(0), AtomAst::from_element(Element::C)),
                (AtomId(1), AtomAst::from_element(Element::C)),
                (AtomId(2), AtomAst::from_element(Element::N)),
                (AtomId(3), AtomAst::from_element(Element::O)),
            ],
        );
    }

    #[rstest]
    #[case::present(AtomId(2), true)]
    #[case::absent(AtomId(999), false)]
    fn test_atom_views_contains(molecule: MoleculeAst, #[case] id: AtomId, #[case] expected: bool) {
        assert_eq!(molecule.atoms().contains(id), expected);
    }

    #[rstest]
    fn test_atom_views_get(molecule: MoleculeAst) {
        let res = molecule.atoms().get(AtomId(2));
        assert!(res.is_some());
        let atom = res.unwrap();
        assert_eq!(atom.id, AtomId(2));
        assert_eq!(atom.ast, &AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_get_none(molecule: MoleculeAst) {
        let res = molecule.atoms().get(AtomId(999));
        assert!(res.is_none());
    }

    #[rstest]
    fn test_atom_views_index(molecule: MoleculeAst) {
        let atom: &AtomAst = &molecule.atoms()[AtomId(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_view_neighbors(molecule: MoleculeAst) {
        let view = molecule.atom(AtomId(1));
        let collected: Vec<(BondId, AtomId, BondAst)> = view
            .neighbors()
            .map(|n| (n.bond_id(), n.atom_id(), n.bond().ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), AtomId(0), BondAst::from_order(1)),
                (BondId(1), AtomId(2), BondAst::from_order(2)),
            ],
        );
    }

    #[rstest]
    #[case::no_incident(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(3),
        ValueAst::Lit(0),
    )]
    #[case::single(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(0),
        ValueAst::Lit(1),
    )]
    #[case::three_around_center(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(1),
        ValueAst::Lit(3),
    )]
    #[case::double(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(2),
        ValueAst::Lit(2),
    )]
    #[case::undetermined_bond(
        mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#),
        AtomId(0),
        ValueAst::Undetermined,
    )]
    fn test_atom_view_valence(
        #[case] molecule: MoleculeAst,
        #[case] center: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(molecule.atom(center).valence(), expected);
    }

    #[rstest]
    #[case::with_constraint(Some(AtomConstraint::valence(4)), ValueAst::Lit(4))]
    #[case::absent(None, ValueAst::Undetermined)]
    fn test_atom_view_valence_constraint(
        #[case] constraint: Option<AtomConstraint>,
        #[case] expected: ValueAst,
    ) {
        let mut atom = AtomAst::from_element(Element::C);
        if let Some(c) = constraint {
            atom.constraints.add(c);
        }
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(molecule.atom(AtomId(0)).constraints().valence(), expected);
    }

    #[rstest]
    #[case::donor(AtomId(0), ValueAst::Lit(1))]
    #[case::acceptor(AtomId(1), ValueAst::Lit(0))]
    fn test_atom_view_donated_pairs(#[case] atom: AtomId, #[case] expected: ValueAst) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).donated_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_donated_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::N);
        atom.constraints.add(AtomConstraint::donated_pairs(1));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().donated_pairs(),
            ValueAst::Lit(1),
        );
    }

    #[rstest]
    #[case::donor(AtomId(0), ValueAst::Lit(0))]
    #[case::acceptor(AtomId(1), ValueAst::Lit(1))]
    fn test_atom_view_accepted_pairs(#[case] atom: AtomId, #[case] expected: ValueAst) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).accepted_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_accepted_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::accepted_pairs(2));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().accepted_pairs(),
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_not_in_system() {
        let molecule = mol!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).aromatic_valence(),
            ValueAst::Lit(0)
        );
    }

    #[rstest]
    #[case::in_system(AtomId(0), true)]
    #[case::not_in_system(AtomId(3), false)]
    fn test_atom_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.atom(atom).is_in_aromatic_system(), expected);
    }

    #[rstest]
    #[case::participant(AtomId(0), Some(AromaticSystemId(0)))]
    #[case::not_participant(AtomId(3), None)]
    fn test_atom_view_aromatic_system(
        molecule: MoleculeAst,
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
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<DativeBondId>,
    ) {
        let ids: Vec<DativeBondId> = molecule.atom(atom).dative_bonds().map(|v| v.id).collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![MulticenterBondId(0)])]
    #[case::uninvolved(AtomId(3), vec![])]
    fn test_atom_view_multicenter_bonds(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<MulticenterBondId>,
    ) {
        let ids: Vec<MulticenterBondId> = molecule
            .atom(atom)
            .multicenter_bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::endpoint_0(AtomId(0), vec![NoncovalentBondId(0)])]
    #[case::endpoint_3(AtomId(3), vec![NoncovalentBondId(0)])]
    #[case::uninvolved(AtomId(1), vec![])]
    fn test_atom_view_noncovalent_bonds(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<NoncovalentBondId>,
    ) {
        let ids: Vec<NoncovalentBondId> = molecule
            .atom(atom)
            .noncovalent_bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[fixture]
    fn stereo_molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 10],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(0), AtomId(2), BondAst::from_order(1)),
                (AtomId(0), AtomId(3), BondAst::from_order(1)),
                (AtomId(0), AtomId(4), BondAst::from_order(1)),
                (AtomId(5), AtomId(6), BondAst::from_order(1)),
                (AtomId(5), AtomId(7), BondAst::from_order(1)),
                (AtomId(5), AtomId(8), BondAst::from_order(1)),
                (AtomId(5), AtomId(9), BondAst::from_order(1)),
            ],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
                ),
                (
                    AtomId(5),
                    vec![
                        StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(9), StereoLigandKind::Atom),
                    ],
                    StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(1)),
                ),
            ],
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    #[case::tetrahedral_site(AtomId(0), true)]
    #[case::square_planar_site(AtomId(5), false)]
    #[case::ligand(AtomId(1), false)]
    fn test_atom_view_is_in_tetrahedral_stereo(
        stereo_molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(
            stereo_molecule.atom(atom).is_in_tetrahedral_stereo(),
            expected
        );
    }

    #[rstest]
    #[case::tetrahedral_site(AtomId(0), Some(StereoAtomId(0)))]
    #[case::square_planar_site(AtomId(5), None)]
    #[case::ligand(AtomId(1), None)]
    fn test_atom_view_tetrahedral_stereo_id(
        stereo_molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(stereo_molecule.atom(atom).tetrahedral_stereo_id(), expected);
    }

    #[rstest]
    fn test_atom_view_tetrahedral_stereo(stereo_molecule: MoleculeAst) {
        let view = stereo_molecule
            .atom(AtomId(0))
            .tetrahedral_stereo()
            .unwrap();
        assert_eq!(view.id, StereoAtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert!(stereo_molecule
            .atom(AtomId(5))
            .tetrahedral_stereo()
            .is_none());
    }

    #[rstest]
    #[case::ring_atom_0(AtomId(0), true)]
    #[case::ring_atom_3(AtomId(3), true)]
    #[case::ring_atom_5(AtomId(5), true)]
    #[case::chain_atom_6(AtomId(6), false)]
    fn test_atom_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(ring_with_chain.atom(atom).is_in_ring(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), true)]
    #[case::chain_atom(AtomId(6), false)]
    fn test_atom_view_is_in_ring_from(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        assert_eq!(ring_with_chain.atom(atom).is_in_ring_from(&rings), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), 1)]
    #[case::chain_atom(AtomId(6), 0)]
    fn test_atom_view_rings_from(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected_count: usize,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        let count = ring_with_chain.atom(atom).rings_from(&rings).count();
        assert_eq!(count, expected_count);
    }

    #[rstest]
    #[case::aromatic_and_multicenter(molecule(), AtomId(0), true)]
    #[case::aromatic_only_in_rich(molecule(), AtomId(1), true)]
    #[case::dative_donor(molecule(), AtomId(2), true)]
    #[case::dative_acceptor(molecule(), AtomId(3), true)]
    #[case::bare_atom_0(
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        AtomId(0),
        false,
    )]
    #[case::bare_atom_1(
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        AtomId(1),
        false,
    )]
    #[case::tetrahedral_stereo_site(stereo_molecule(), AtomId(0), true)]
    fn test_atom_view_is_in_overlays(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(mol.atom(atom).is_in_overlays(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_site(AtomId(0), AtomConstraints::from_iter([
        AtomConstraint::valence(ValueAst::Lit(4)),
        AtomConstraint::donated_pairs(ValueAst::Lit(0)),
        AtomConstraint::accepted_pairs(ValueAst::Lit(0)),
        AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraint::tetrahedral_stereo(TetrahedralStereoAst::stereo(StereoCosetAst::Lit(1))),
    ]))]
    #[case::non_stereo_ligand(AtomId(1), AtomConstraints::from_iter([
        AtomConstraint::valence(ValueAst::Lit(1)),
        AtomConstraint::donated_pairs(ValueAst::Lit(0)),
        AtomConstraint::accepted_pairs(ValueAst::Lit(0)),
        AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraint::tetrahedral_stereo(TetrahedralStereoAst::NotStereo),
    ]))]
    #[case::square_planar_site(AtomId(5), AtomConstraints::from_iter([
        AtomConstraint::valence(ValueAst::Lit(4)),
        AtomConstraint::donated_pairs(ValueAst::Lit(0)),
        AtomConstraint::accepted_pairs(ValueAst::Lit(0)),
        AtomConstraint::aromatic_valence(AromaticValenceAst::NotAromatic),
        AtomConstraint::multicenter_valence(MulticenterValenceAst::NotMulticenter),
        AtomConstraint::tetrahedral_stereo(TetrahedralStereoAst::NotStereo),
    ]))]
    fn test_atom_view_derive_constraints(
        stereo_molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: AtomConstraints,
    ) {
        assert_eq!(stereo_molecule.atom(atom).derive_constraints(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::aromatic_valence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        ));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().aromatic_valence(),
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        );
    }

    #[rstest]
    #[case::single_bond(
        vec![(vec![AtomId(0), AtomId(1)], ElectronCountsAst::Lit(vec![2, 2]))],
        ValueAst::Lit(2),
    )]
    #[case::two_bonds(
        vec![
            (vec![AtomId(0), AtomId(1)], ElectronCountsAst::Lit(vec![2, 2])),
            (vec![AtomId(0), AtomId(2)], ElectronCountsAst::Lit(vec![1, 1])),
        ],
        ValueAst::Lit(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomId(0), AtomId(1)], ElectronCountsAst::Undetermined)],
        ValueAst::Undetermined,
    )]
    fn test_atom_view_multicenter_valence(
        #[case] bonds: Vec<(Vec<AtomId>, ElectronCountsAst)>,
        #[case] expected: ValueAst,
    ) {
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, MulticenterBondAst::new(electrons)))
            .collect();
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![],
            vec![],
            multicenter,
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        assert_eq!(molecule.atom(AtomId(0)).multicenter_valence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::multicenter_valence(
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        ));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().multicenter_valence(),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        );
    }

    #[rstest]
    #[case::ethane_carbon(mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#), AtomId(0), ValueAst::Lit(1))]
    #[case::ethene_carbon(mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#), AtomId(0), ValueAst::Lit(1))]
    #[case::three_bonds(mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#), AtomId(0), ValueAst::Lit(3))]
    fn test_atom_view_degree(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).degree(), expected);
    }

    #[rstest]
    #[case::lit(mol!(r#"{:atoms ["C#h4"] :bonds []}"#), ValueAst::Lit(4))]
    #[case::undetermined(mol!(r#"{:atoms ["C#h*"] :bonds []}"#), ValueAst::Undetermined)]
    fn test_atom_view_total_degree(#[case] molecule: MoleculeAst, #[case] expected: ValueAst) {
        assert_eq!(molecule.atom(AtomId(0)).total_degree(), expected);
    }

    #[rstest]
    #[case::all_heavy(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(2),
    )]
    #[case::one_h_neighbor(
        mol!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(1),
    )]
    fn test_atom_view_heavy_atom_degree(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_degree(), expected);
    }

    #[rstest]
    #[case::all_heavy(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "2"]]}"#),
        AtomId(0),
        ValueAst::Lit(3),
    )]
    #[case::skips_h(
        mol!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "2"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(2),
    )]
    fn test_atom_view_heavy_atom_valence(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_valence(), expected);
    }

    #[rstest]
    #[case::implicit_only(
        mol!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    #[case::implicit_and_explicit(
        mol!(r#"{:atoms ["C#h2" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    #[case::implicit_undetermined(
        mol!(r#"{:atoms ["C#h*"] :bonds []}"#),
        AtomId(0),
        ValueAst::Undetermined,
    )]
    fn test_atom_view_total_hydrogens(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).total_hydrogens(), expected);
    }

    #[rstest]
    #[case::lit(mol!(r#"{:atoms ["C#h4"] :bonds []}"#), ValueAst::Lit(4))]
    #[case::undetemined(mol!(r#"{:atoms ["C#h*"] :bonds []}"#), ValueAst::Undetermined)]
    fn test_atom_view_total_valence_sum_of_terms(
        #[case] molecule: MoleculeAst,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(molecule.atom(AtomId(0)).total_valence(), expected);
    }

    #[rstest]
    #[case::ch4(mol!(r#"{:atoms ["C#h4"] :bonds []}"#), ValueAst::Lit(4))]
    #[case::undetermined_h(mol!(r#"{:atoms ["C#h*"] :bonds []}"#), ValueAst::Undetermined)]
    fn test_atom_view_covalence_non_aromatic(
        #[case] molecule: MoleculeAst,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(molecule.atom(AtomId(0)).covalence(), expected);
    }

    #[fixture]
    fn aromatic_ring() -> MoleculeAst {
        // 3-membered C ring, each with 0 implicit H (valence 2 from two ring
        // bonds), aromatic system electrons [1, 2, 0].
        let carbon = AtomAst::from_element(Element::C).with_implicit_hydrogens(0_i64);
        MoleculeAst::from_parts(
            vec![carbon.clone(), carbon.clone(), carbon],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(0), BondAst::from_order(1)),
            ],
            vec![],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::from_counts(vec![1, 2, 0]),
            )],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    #[case::standard(AtomId(0), ValueAst::Lit(3))] // av=1 → +1
    #[case::donor(AtomId(1), ValueAst::Lit(2))] // av=2 (donated pair) → +0
    #[case::acceptor(AtomId(2), ValueAst::Lit(2))] // av=0 → +0
    fn test_atom_view_covalence_aromatic(
        aromatic_ring: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(aromatic_ring.atom(atom).covalence(), expected);
    }

    #[fixture]
    fn dative_pair() -> MoleculeAst {
        // H₃N→BH₃: N (3 H) donates a pair to B (3 H). Covalence = v+h+ai for
        // both = 3; the dative bond (donated on N, accepted on B) is excluded.
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N).with_implicit_hydrogens(3_i64),
                AtomAst::from_element(Element::B).with_implicit_hydrogens(3_i64),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    #[case::donor(AtomId(0), ValueAst::Lit(3))] // donated pair excluded → v+h = 3
    #[case::acceptor(AtomId(1), ValueAst::Lit(3))] // accepted pair excluded → v+h = 3
    fn test_atom_view_covalence_dative(
        dative_pair: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(dative_pair.atom(atom).covalence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_degree() {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![],
            vec![],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::from_counts(vec![2, 2, 2]),
            )],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );
        assert_eq!(
            molecule.atom(AtomId(0)).multicenter_degree(),
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    #[case::all_ring_atom(AtomId(0), RingScope::All, ValueAst::Lit(1))]
    #[case::all_ring_atom_alt(AtomId(3), RingScope::All, ValueAst::Lit(1))]
    #[case::size_match(AtomId(0), RingScope::Size(6), ValueAst::Lit(1))]
    #[case::size_no_match(AtomId(0), RingScope::Size(5), ValueAst::Lit(0))]
    #[case::all_chain_atom(AtomId(6), RingScope::All, ValueAst::Lit(0))]
    #[case::size_chain_atom(AtomId(6), RingScope::Size(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_membership(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] scope: RingScope,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_membership(scope), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(1))]
    #[case::ring_atom_alt(AtomId(3), ValueAst::Lit(1))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_count(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_count(), expected);
    }

    #[rstest]
    #[case::size_match(AtomId(0), 6, ValueAst::Lit(1))]
    #[case::size_no_match(AtomId(0), 5, ValueAst::Lit(0))]
    #[case::chain_atom(AtomId(6), 6, ValueAst::Lit(0))]
    fn test_atom_view_ring_size_count(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] size: u8,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_size_count(size), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), Some(6))]
    #[case::chain_atom(AtomId(6), None)]
    fn test_atom_view_smallest_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(ring_with_chain.atom(atom).smallest_ring_size(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(2))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_degree(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_degree(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(2))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_valence(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_valence(), expected);
    }
}
