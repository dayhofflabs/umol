//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, FixedRelationSet, Graph, MatchingEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId, RelationId, ShortestCycleAlgorithm,
    SubgraphIsomorphismAlgorithm, VarRelationSet,
};
use umol_shared::element::Element;

use super::automorphism::AtomAutomorphism;
use super::matching::BondMatching;

use super::aromatic::AromaticSystemAst;
use super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraints, AtomConstraints, BondConstraints, DativeBondConstraints,
    MulticenterBondConstraints, NoncovalentBondConstraints,
};
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use super::rings::{RingSet, RingView};
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Namespace accessor for atom views on a `MoleculeAst`. Provides `count`,
/// `ids`, `iter`, `get`, and `Index` without burying them on `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    molecule: &'a MoleculeAst,
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, atoms: &'a [AtomAst]) -> Self {
        Self { molecule, atoms }
    }

    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomId> {
        (0..self.atoms.len() as u32).map(AtomId)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .enumerate()
            .map(move |(i, ast)| AtomView {
                id: AtomId(i as u32),
                ast,
                molecule,
            })
    }

    pub fn get(&self, id: AtomId) -> AtomView<'a> {
        AtomView {
            id,
            ast: &self.atoms[id.index()],
            molecule: self.molecule,
        }
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
    pub fn isotope_mass(&self) -> &'a IsotopeAst {
        &self.ast.isotope_mass
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn implicit_hydrogens(&self) -> &'a ImplicitHydrogensAst {
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
    pub fn neighbors(&self) -> impl Iterator<Item = NeighborView<'a>> {
        self.molecule.neighbors(self.id)
    }

    /// IDs of incident localized bonds, in iteration order of `neighbors`.
    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule.neighbors(self.id).map(|n| n.bond)
    }

    /// Localized valence: sum of incident `Bond.order` values. Returns
    /// `ValueAst::Lit(n)` when every incident bond order is `Lit`; collapses
    /// to `Undetermined` if any bond order is non-`Lit`.
    pub fn valence(&self) -> ValueAst {
        self.neighbors()
            .map(|n| n.ast.order.clone())
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
            if view.acceptor_id != self.id {
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
        sys.ast
            .electrons
            .get(pos)
            .cloned()
            .unwrap_or(ValueAst::Undetermined)
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.aromatic_system().is_some()
    }

    /// The aromatic system containing this atom, if any. Per-perception
    /// design an atom belongs to at most one aromatic system; this
    /// returns the first incident system.
    pub fn aromatic_system(&self) -> Option<AromaticSystemView<'a>> {
        self.aromatic_system_id()
            .map(|id| self.molecule.aromatic_system(id))
    }

    pub fn aromatic_system_id(&self) -> Option<AromaticSystemId> {
        self.molecule.aromatic_systems().incident_ids(self.id).next()
    }

    pub fn dative_bonds(&self) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        self.molecule.dative_bonds().incident(self.id)
    }

    pub fn dative_bond_ids(&self) -> impl Iterator<Item = DativeBondId> + 'a {
        self.molecule.dative_bonds().incident_ids(self.id)
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = MulticenterBondView<'a>> + 'a {
        self.molecule.multicenter_bonds().incident(self.id)
    }

    pub fn multicenter_bond_ids(&self) -> impl Iterator<Item = MulticenterBondId> + 'a {
        self.molecule.multicenter_bonds().incident_ids(self.id)
    }

    pub fn noncovalent_bonds(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> + 'a {
        self.molecule.noncovalent_bonds().incident(self.id)
    }

    pub fn noncovalent_bond_ids(&self) -> impl Iterator<Item = NoncovalentBondId> + 'a {
        self.molecule.noncovalent_bonds().incident_ids(self.id)
    }

    /// True if this atom participates in any of the four overlay relations
    /// (aromatic system, dative bond, multicenter bond, noncovalent bond).
    /// Mirror of `MoleculeAst::has_overlays` scoped to a single atom; useful
    /// as a pre-mutation predicate before structural removal.
    pub fn is_in_overlays(&self) -> bool {
        self.aromatic_system().is_some()
            || self.dative_bonds().next().is_some()
            || self.multicenter_bonds().next().is_some()
            || self.noncovalent_bonds().next().is_some()
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

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `ValueAst::Lit(0)` when not in any multicenter bond; collapses to
    /// `Undetermined` if any contribution is non-`Lit`.
    pub fn multicenter_valence(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.multicenter_bonds() {
            let Some(pos) = view.atom_ids().position(|a| a == self.id) else {
                return ValueAst::Undetermined;
            };
            let term = view
                .ast
                .electrons
                .get(pos)
                .cloned()
                .unwrap_or(ValueAst::Undetermined);
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
        self.degree()
            + ValueAst::from(self.ast.implicit_hydrogens.clone())
            + self.multicenter_degree()
    }

    /// Count of incident localized bonds whose neighbor is not a literal
    /// hydrogen atom (Element::H). Always `Lit`; non-`Lit` neighbor
    /// elements count as heavy (i.e., not filtered out).
    pub fn heavy_atom_degree(&self) -> ValueAst {
        let count = self
            .neighbors()
            .filter(|n| {
                !matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .count();
        ValueAst::Lit(count as i64)
    }

    /// `valence` over incident bonds whose neighbor is not a literal
    /// hydrogen. Collapses to `Undetermined` if any contributing bond order
    /// is non-`Lit`.
    pub fn heavy_atom_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| {
                !matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .map(|n| n.ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Explicit hydrogens (incident neighbors with `Element::H`) plus
    /// `implicit_hydrogens`. Collapses to `Undetermined` if `implicit_hydrogens`
    /// is non-`Lit` (including `Normal`).
    pub fn total_hydrogens(&self) -> ValueAst {
        let explicit = self
            .neighbors()
            .filter(|n| {
                matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .count() as i64;
        ValueAst::Lit(explicit) + ValueAst::from(self.ast.implicit_hydrogens.clone())
    }

    /// Full electron-sharing sum at this atom:
    /// `valence + implicit_hydrogens + aromatic_valence + multicenter_valence`.
    /// Diverges from SMARTS `v<n>` for aromatic lone-pair donors (pyrrole N,
    /// furan O) which contribute the donated pair via `aromatic_valence`.
    pub fn total_valence(&self) -> ValueAst {
        self.valence()
            + ValueAst::from(self.ast.implicit_hydrogens.clone())
            + self.aromatic_valence()
            + self.multicenter_valence()
    }

    /// Count of multicenter co-participants across all incident multicenter
    /// bonds. Per the no-overlap structural rule these are not localized-
    /// bond neighbors. Always `Lit`.
    pub fn multicenter_degree(&self) -> ValueAst {
        let count: usize = self
            .multicenter_bonds()
            .map(|mc| mc.ast.electrons.len().saturating_sub(1))
            .sum();
        ValueAst::Lit(count as i64)
    }

    /// Count of canonical rings (Vismara / max ring size 22) containing
    /// this atom. Always `Lit`.
    pub fn ring_count(&self) -> ValueAst {
        ValueAst::Lit(self.rings().count() as i64)
    }

    /// Sizes of canonical rings containing this atom, in iteration order.
    /// Multi-valued: an atom in fused rings yields one size per ring.
    pub fn ring_size(&self) -> impl Iterator<Item = usize> + 'a {
        self.rings().map(|r| r.len())
    }

    /// Smallest containing canonical ring size, or `None` if not in any
    /// ring. Chemistry-classification helper; not the constraint
    /// counterpart of `RingSize` (which uses interpretation B — "in some
    /// ring of size n").
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.ring_size().min()
    }

    /// Count of incident bonds that participate in any canonical ring.
    /// Always `Lit`.
    pub fn ring_degree(&self) -> ValueAst {
        let count = self
            .neighbors()
            .filter(|n| self.molecule.bond(n.bond).is_in_ring())
            .count();
        ValueAst::Lit(count as i64)
    }

    /// Sum of bond orders of incident bonds that participate in any
    /// canonical ring. Collapses to `Undetermined` if any contributing
    /// bond's `order` is non-`Lit`.
    pub fn ring_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| self.molecule.bond(n.bond).is_in_ring())
            .map(|n| n.ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }
}

/// Mutable borrowed view of an atom.
#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomAst,
}

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    molecule: &'a MoleculeAst,
    bonds: &'a [BondAst],
}

impl<'a> BondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, bonds: &'a [BondAst]) -> Self {
        Self { molecule, bonds }
    }

    pub fn count(&self) -> usize {
        self.bonds.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = BondId> {
        (0..self.bonds.len() as u32).map(BondId)
    }

    pub fn iter(&self) -> impl Iterator<Item = BondView<'a>> {
        let molecule = self.molecule;
        let bonds = self.bonds;
        let graph = molecule.raw_graph();
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                id: BondId::from(id),
                atoms: [AtomId::from(s), AtomId::from(t)],
                ast: &bonds[id.index()],
                molecule,
            }
        })
    }

    pub fn get(&self, id: BondId) -> BondView<'a> {
        let [s, t] = self.molecule.raw_graph().edge_endpoints(EdgeId::from(id));
        BondView {
            id,
            atoms: [AtomId::from(s), AtomId::from(t)],
            ast: &self.bonds[id.index()],
            molecule: self.molecule,
        }
    }

    /// ID of the bond between `a` and `b`, if any.
    pub fn connecting_id(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.molecule
            .raw_graph()
            .find_edge(NodeId::from(a), NodeId::from(b))
            .map(BondId::from)
    }

    /// View of the bond between `a` and `b`, if any.
    pub fn connecting(&self, a: AtomId, b: AtomId) -> Option<BondView<'a>> {
        self.connecting_id(a, b).map(|id| self.get(id))
    }

    /// IDs of bonds whose both endpoints lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<BondId> {
        let mut nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        nodes.sort_unstable();
        self.molecule
            .raw_graph()
            .induced_edges(&nodes)
            .map(BondId::from)
            .collect()
    }

    /// Views of bonds whose both endpoints lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<BondView<'a>> {
        self.induced_ids(atoms).into_iter().map(|id| self.get(id)).collect()
    }
}

impl<'a> Index<BondId> for BondViews<'a> {
    type Output = BondAst;
    fn index(&self, id: BondId) -> &BondAst {
        &self.bonds[id.index()]
    }
}

/// Borrowed view of a bond: its index, the two participating atoms, and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub id: BondId,
    atoms: [AtomId; 2],
    pub ast: &'a BondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> BondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a BondConstraints {
        &self.ast.constraints
    }

    /// The two atom indices incident to this bond.
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.atoms
    }

    /// Views of the two atoms incident to this bond.
    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        let [a, b] = self.atoms;
        [a, b].into_iter().map(move |id| molecule.atom(id))
    }

    /// The aromatic system this bond participates in, if any. A bond is in
    /// an aromatic system iff both endpoints belong to that system.
    pub fn aromatic_system(&self) -> Option<AromaticSystemView<'a>> {
        let [a, b] = self.atoms;
        self.molecule
            .aromatic_systems()
            .incident(a)
            .find(|sys| sys.atom_ids().any(|x| x == b))
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.aromatic_system().is_some()
    }

    /// True if this bond belongs to any ring in the molecule's canonical
    /// ring set (Vismara relevant cycles, max ring size 22). Uses the
    /// molecule's cached canonical `RingSet`.
    pub fn is_in_ring(&self) -> bool {
        self.molecule.rings().contains_bond(self.id)
    }

    /// True if this bond appears in any ring of the supplied set.
    pub fn is_in_ring_from(&self, rings: &RingSet) -> bool {
        rings.contains_bond(self.id)
    }

    /// Rings containing this bond drawn from the molecule's canonical
    /// `RingSet` (Vismara relevant cycles, max ring size 22).
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.molecule
            .rings()
            .iter()
            .filter(move |v| v.bonds().contains(&id))
    }

    /// Rings from the supplied set that contain this bond.
    pub fn rings_from<'r>(&self, rings: &'r RingSet) -> impl Iterator<Item = RingView<'r>> + 'r {
        let id = self.id;
        rings.iter().filter(move |v| v.bonds().contains(&id))
    }

    /// Count of canonical rings (Vismara / max ring size 22) containing
    /// this bond. Always `Lit`.
    pub fn ring_count(&self) -> ValueAst {
        ValueAst::Lit(self.rings().count() as i64)
    }

    /// Sizes of canonical rings containing this bond, in iteration order.
    /// Multi-valued: a bond shared between fused rings yields one size per
    /// ring.
    pub fn ring_size(&self) -> impl Iterator<Item = usize> + 'a {
        self.rings().map(|r| r.len())
    }
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub id: BondId,
    atoms: [AtomId; 2],
    pub ast: &'a mut BondAst,
}

impl<'a> BondViewMut<'a> {
    pub(super) fn new(id: BondId, atoms: [AtomId; 2], ast: &'a mut BondAst) -> Self {
        Self { id, atoms, ast }
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<DativeBondAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondId> {
        self.set.relation_ids().map(DativeBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let atoms = set.participants(rid);
            let ast = set.data(rid);
            let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
            DativeBondView {
                id: DativeBondId::from(rid),
                ast,
                acceptor_id,
                atoms,
                molecule,
            }
        })
    }

    pub fn get(&self, id: DativeBondId) -> DativeBondView<'a> {
        let rid = RelationId::from(id);
        let atoms = self.set.participants(rid);
        let ast = self.set.data(rid);
        let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
        DativeBondView {
            id,
            ast,
            acceptor_id,
            atoms,
            molecule: self.molecule,
        }
    }

    /// IDs of dative bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = DativeBondId> + 'a {
        self.set
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondId::from(rid))
    }

    /// Views of dative bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.set;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            let atoms = set.participants(rid);
            let ast = set.data(rid);
            let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
            DativeBondView { id, ast, acceptor_id, atoms, molecule }
        })
    }

    /// ID of the dative bond whose participant set equals `atoms`, if any.
    pub fn connecting_id(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<DativeBondId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> =
                self.set.participants(RelationId::from(id))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect();
            parts == target
        })
    }

    /// View of the dative bond whose participant set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<DativeBondView<'a>> {
        self.connecting_id(atoms).map(|id| self.get(id))
    }

    /// IDs of dative bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<DativeBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.set
            .relation_ids()
            .filter(|&rid| {
                self.set
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(DativeBondId::from)
            .collect()
    }

    /// Views of dative bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<DativeBondView<'a>> {
        self.induced_ids(atoms).into_iter().map(|id| self.get(id)).collect()
    }
}

impl<'a> Index<DativeBondId> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, id: DativeBondId) -> &DativeBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondAst`. Donor atoms and the full participant
/// set are reachable through `donors()` and `atoms()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub id: DativeBondId,
    pub acceptor_id: AtomId,
    atoms: &'a [NodeId],
    pub ast: &'a DativeBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn acceptor_slot(&self) -> u8 {
        self.ast.acceptor_slot
    }

    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn constraints(&self) -> &'a DativeBondConstraints {
        &self.ast.constraints
    }

    /// All atoms in this dative bond (donors + acceptor), sorted by `AtomId`.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    /// Views of all atoms in this dative bond (donors + acceptor).
    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    /// Donor atom ids (participants minus the acceptor slot).
    pub fn donor_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let acceptor_slot = self.ast.acceptor_slot as usize;
        self.atoms
            .iter()
            .enumerate()
            .filter(move |(i, _)| *i != acceptor_slot)
            .map(|(_, &n)| AtomId::from(n))
    }

    /// Donor atom views (participants minus the acceptor slot).
    pub fn donors(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.donor_ids().map(move |id| molecule.atom(id))
    }

    /// View of the acceptor atom.
    pub fn acceptor(&self) -> AtomView<'a> {
        self.molecule.atom(self.acceptor_id)
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<AromaticSystemAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemId> {
        self.set.relation_ids().map(AromaticSystemId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| AromaticSystemView {
            id: AromaticSystemId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    pub fn get(&self, id: AromaticSystemId) -> AromaticSystemView<'a> {
        let rid = RelationId::from(id);
        AromaticSystemView {
            id,
            ast: self.set.data(rid),
            atoms: self.set.participants(rid),
            molecule: self.molecule,
        }
    }

    /// IDs of aromatic systems incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = AromaticSystemId> + 'a {
        self.set
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| AromaticSystemId::from(rid))
    }

    /// Views of aromatic systems incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = AromaticSystemView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.set;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            AromaticSystemView { id, ast: set.data(rid), atoms: set.participants(rid), molecule }
        })
    }

    /// ID of the aromatic system whose atom set equals `atoms`, if any.
    pub fn connecting_id(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<AromaticSystemId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> =
                self.set.participants(RelationId::from(id))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect();
            parts == target
        })
    }

    /// View of the aromatic system whose atom set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<AromaticSystemView<'a>> {
        self.connecting_id(atoms).map(|id| self.get(id))
    }

    /// IDs of aromatic systems whose atoms all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<AromaticSystemId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.set
            .relation_ids()
            .filter(|&rid| {
                self.set
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(AromaticSystemId::from)
            .collect()
    }

    /// Views of aromatic systems whose atoms all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<AromaticSystemView<'a>> {
        self.induced_ids(atoms).into_iter().map(|id| self.get(id)).collect()
    }
}

impl<'a> Index<AromaticSystemId> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, id: AromaticSystemId) -> &AromaticSystemAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemAst`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub ast: &'a AromaticSystemAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AromaticSystemView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a [ValueAst] {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a AromaticSystemConstraints {
        &self.ast.constraints
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(BondId::from)
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondView<'a>> + 'a {
        let molecule = self.molecule;
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(move |edge| molecule.bond(BondId::from(edge)))
    }

    /// Sum of per-atom electron contributions on this aromatic system.
    /// `Lit(n)` when every entry is `Lit`; collapses to `Undetermined` if
    /// any entry is non-`Lit`.
    pub fn electron_count(&self) -> ValueAst {
        self.ast
            .electrons
            .iter()
            .cloned()
            .fold(ValueAst::Lit(0), |acc, e| acc + e)
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bond_ids().count()
    }

    /// Atom views for atoms in this system that also appear in `subset`.
    pub fn overlapping_atoms<'s>(
        &self,
        subset: &'s [AtomId],
    ) -> impl Iterator<Item = AtomView<'a>> + 's
    where
        'a: 's,
    {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(|&n| AtomId::from(n))
            .filter(move |a| subset.contains(a))
            .map(move |id| molecule.atom(id))
    }

    /// Bond views for bonds in this system that also appear in `subset`.
    pub fn overlapping_bonds<'s>(
        &self,
        subset: &'s [BondId],
    ) -> impl Iterator<Item = BondView<'a>> + 's
    where
        'a: 's,
    {
        let molecule = self.molecule;
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(BondId::from)
            .filter(move |b| subset.contains(b))
            .map(move |id| molecule.bond(id))
    }

    /// Rings from the molecule's canonical `RingSet` that share at least
    /// one atom with this aromatic system.
    pub fn overlapping_rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let atoms: Vec<AtomId> = self.atoms.iter().map(|&n| AtomId::from(n)).collect();
        self.molecule
            .rings()
            .iter()
            .filter(move |r| r.atoms().iter().any(|a| atoms.contains(a)))
    }
}

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<MulticenterBondAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = MulticenterBondId> {
        self.set.relation_ids().map(MulticenterBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| MulticenterBondView {
            id: MulticenterBondId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    pub fn get(&self, id: MulticenterBondId) -> MulticenterBondView<'a> {
        let rid = RelationId::from(id);
        MulticenterBondView {
            id,
            ast: self.set.data(rid),
            atoms: self.set.participants(rid),
            molecule: self.molecule,
        }
    }

    /// IDs of multicenter bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = MulticenterBondId> + 'a {
        self.set
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| MulticenterBondId::from(rid))
    }

    /// Views of multicenter bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = MulticenterBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.set;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            MulticenterBondView { id, ast: set.data(rid), atoms: set.participants(rid), molecule }
        })
    }

    /// ID of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn connecting_id(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<MulticenterBondId> {
        let target: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = target.iter().next()?;
        self.incident_ids(first).find(|&id| {
            let parts: HashSet<AtomId> =
                self.set.participants(RelationId::from(id))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect();
            parts == target
        })
    }

    /// View of the multicenter bond whose participant set equals `atoms`, if any.
    pub fn connecting(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<MulticenterBondView<'a>> {
        self.connecting_id(atoms).map(|id| self.get(id))
    }

    /// IDs of multicenter bonds whose participants all lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<MulticenterBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.set
            .relation_ids()
            .filter(|&rid| {
                self.set
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(MulticenterBondId::from)
            .collect()
    }

    /// Views of multicenter bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<MulticenterBondView<'a>> {
        self.induced_ids(atoms).into_iter().map(|id| self.get(id)).collect()
    }
}

impl<'a> Index<MulticenterBondId> for MulticenterBondViews<'a> {
    type Output = MulticenterBondAst;
    fn index(&self, id: MulticenterBondId) -> &MulticenterBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a multicenter bond: its index, member atoms via
/// `atoms()`, and underlying `MulticenterBondAst`.
#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub ast: &'a MulticenterBondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> MulticenterBondView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a [ValueAst] {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a MulticenterBondConstraints {
        &self.ast.constraints
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(move |&n| molecule.atom(AtomId::from(n)))
    }

    /// Sum of per-atom electron contributions on this multicenter bond.
    /// `Lit(n)` when every entry is `Lit`; collapses to `Undetermined` if
    /// any entry is non-`Lit`.
    pub fn electron_count(&self) -> ValueAst {
        self.ast
            .electrons
            .iter()
            .cloned()
            .fold(ValueAst::Lit(0), |acc, e| acc + e)
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Atom views for atoms in this multicenter bond that also appear in `subset`.
    pub fn overlapping_atoms<'s>(
        &self,
        subset: &'s [AtomId],
    ) -> impl Iterator<Item = AtomView<'a>> + 's
    where
        'a: 's,
    {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .map(|&n| AtomId::from(n))
            .filter(move |a| subset.contains(a))
            .map(move |id| molecule.atom(id))
    }
}

/// Namespace accessor for noncovalent-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(super) fn new(
        molecule: &'a MoleculeAst,
        set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
    ) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = NoncovalentBondId> {
        self.set.relation_ids().map(NoncovalentBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| NoncovalentBondView {
            id: NoncovalentBondId::from(rid),
            ast: set.data(rid),
            atoms: {
                let parts = set.participants(rid);
                [AtomId::from(parts[0]), AtomId::from(parts[1])]
            },
            molecule,
        })
    }

    pub fn get(&self, id: NoncovalentBondId) -> NoncovalentBondView<'a> {
        let rid = RelationId::from(id);
        let parts = self.set.participants(rid);
        NoncovalentBondView {
            id,
            ast: self.set.data(rid),
            atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
            molecule: self.molecule,
        }
    }

    /// IDs of noncovalent bonds incident on `atom`.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = NoncovalentBondId> + 'a {
        self.set
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| NoncovalentBondId::from(rid))
    }

    /// Views of noncovalent bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = NoncovalentBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.set;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            let parts = set.participants(rid);
            NoncovalentBondView {
                id,
                ast: set.data(rid),
                atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
                molecule,
            }
        })
    }

    /// ID of the noncovalent bond between `a` and `b`, if any.
    pub fn connecting_id(&self, a: AtomId, b: AtomId) -> Option<NoncovalentBondId> {
        self.incident_ids(a).find(|&id| {
            let parts = self.set.participants(RelationId::from(id));
            let x = AtomId::from(parts[0]);
            let y = AtomId::from(parts[1]);
            (x == a && y == b) || (x == b && y == a)
        })
    }

    /// View of the noncovalent bond between `a` and `b`, if any.
    pub fn connecting(&self, a: AtomId, b: AtomId) -> Option<NoncovalentBondView<'a>> {
        self.connecting_id(a, b).map(|id| self.get(id))
    }

    /// IDs of noncovalent bonds whose endpoints both lie in `atoms`.
    pub fn induced_ids(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.set
            .relation_ids()
            .filter(|&rid| {
                self.set
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(NoncovalentBondId::from)
            .collect()
    }

    /// Views of noncovalent bonds whose endpoints both lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondView<'a>> {
        self.induced_ids(atoms).into_iter().map(|id| self.get(id)).collect()
    }
}

impl<'a> Index<NoncovalentBondId> for NoncovalentBondViews<'a> {
    type Output = NoncovalentBondAst;
    fn index(&self, id: NoncovalentBondId) -> &NoncovalentBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a noncovalent bond: the two participating atoms plus data.
#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    pub id: NoncovalentBondId,
    atoms: [AtomId; 2],
    pub ast: &'a NoncovalentBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> NoncovalentBondView<'a> {
    #[inline]
    pub fn kind(&self) -> &'a NoncovalentBondKindAst {
        &self.ast.kind
    }

    #[inline]
    pub fn constraints(&self) -> &'a NoncovalentBondConstraints {
        &self.ast.constraints
    }

    /// The two atom ids in this noncovalent interaction.
    pub fn atom_ids(&self) -> [AtomId; 2] {
        self.atoms
    }

    /// Views of the two atoms in this noncovalent interaction.
    pub fn atoms(&self) -> [AtomView<'a>; 2] {
        let [a, b] = self.atoms;
        [self.molecule.atom(a), self.molecule.atom(b)]
    }
}

/// Neighbor-side view of a bond: the atom on the other end (`atom`), the
/// bond index, the bond data, and the parent `MoleculeAst` for navigation
/// to the neighbor's full atom view. Yielded by `MoleculeAst::neighbors`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub bond: BondId,
    pub atom: AtomId,
    pub ast: &'a BondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> NeighborView<'a> {
    pub(super) fn new(
        bond: BondId,
        atom: AtomId,
        ast: &'a BondAst,
        molecule: &'a MoleculeAst,
    ) -> Self {
        Self {
            bond,
            atom,
            ast,
            molecule,
        }
    }
}

/// AtomId/BondId-typed adapter over the underlying `Graph`. Holds the
/// pure-graph algorithms (connectivity, cycles, matchings, isomorphisms)
/// without exposing graph-core's `NodeId` / `EdgeId` types in the public
/// API. Construct via `MoleculeAst::graph()`.
#[derive(Clone, Copy)]
pub struct GraphView<'a> {
    graph: &'a Graph,
}

impl<'a> GraphView<'a> {
    pub(super) fn new(graph: &'a Graph) -> Self {
        Self { graph }
    }

    pub fn degree(&self, atom: AtomId) -> usize {
        self.graph.degree(NodeId::from(atom))
    }

    pub fn connected_components(&self, alg: ConnectedComponentsAlgorithm) -> Vec<Vec<AtomId>> {
        self.graph
            .connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn biconnected_components(
        &self,
        alg: BiconnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .biconnected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn shortest_cycle_through_bond(
        &self,
        bond: BondId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_edge(EdgeId::from(bond), alg)
    }

    pub fn shortest_cycle_through_atom(
        &self,
        atom: AtomId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_node(NodeId::from(atom), alg)
    }

    pub fn enumerate_cycles(
        &self,
        max_size: usize,
        alg: CycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .enumerate_cycles(max_size, alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn maximum_independent_set(&self, alg: MaxIndependentSetAlgorithm) -> Vec<AtomId> {
        self.graph
            .maximum_independent_set(alg)
            .into_iter()
            .map(AtomId::from)
            .collect()
    }

    pub fn maximum_matching(&self, alg: MaxMatchingAlgorithm) -> BondMatching {
        BondMatching(self.graph.maximum_matching(alg))
    }

    pub fn enumerate_perfect_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_perfect_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
    }

    pub fn enumerate_maximum_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_maximum_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
    }

    pub fn automorphisms<C: Ord + Copy>(
        &self,
        atom_color: impl Fn(AtomId) -> C,
        alg: AutomorphismAlgorithm,
    ) -> AtomAutomorphism {
        AtomAutomorphism(
            self.graph
                .automorphisms(|n| atom_color(AtomId::from(n)), alg),
        )
    }

    pub fn subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms(
                query.graph,
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn subgraph_isomorphisms_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms_at(
                query.graph,
                (NodeId::from(anchor.0), NodeId::from(anchor.1)),
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraint, Constraints, MulticenterValenceAst,
    };
    use crate::ast::dative::DativeBondAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::rings::RingFamily;
    use crate::ast::value::ValueAst;
    use crate::mol;

    /// 4-atom molecule with one of every relation kind:
    /// atoms C C N O; bonds 0-1 single, 1-2 double, 2-3 single;
    /// dative donor=2 → acceptor=3; aromatic system [0,1,2];
    /// multicenter bond [0,1,2]; noncovalent H-bond 0-3.
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
            Constraints::default(),
        )
    }

    // --- AtomViews ---

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
    fn test_atom_views_get(molecule: MoleculeAst) {
        let view = molecule.atoms().get(AtomId(2));
        assert_eq!(view.id, AtomId(2));
        assert_eq!(*view.ast, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_index(molecule: MoleculeAst) {
        let atom: &AtomAst = &molecule.atoms()[AtomId(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    // --- AtomView ---

    #[rstest]
    fn test_atom_view_neighbors(molecule: MoleculeAst) {
        let view = molecule.atom(AtomId(1));
        let collected: Vec<(BondId, AtomId, BondAst)> = view
            .neighbors()
            .map(|n| (n.bond, n.atom, n.ast.clone()))
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
        assert_eq!(molecule.atom(AtomId(0)).aromatic_valence(), ValueAst::Lit(0));
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
        let ids: Vec<MulticenterBondId> =
            molecule.atom(atom).multicenter_bonds().map(|v| v.id).collect();
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
        let ids: Vec<NoncovalentBondId> =
            molecule.atom(atom).noncovalent_bonds().map(|v| v.id).collect();
        assert_eq!(ids, expected);
    }

    /// Cyclohexane with one chain atom: 0-1-2-3-4-5-0 closing the ring, plus 0-6 dangling.
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
    fn test_atom_view_is_in_overlays(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(mol.atom(atom).is_in_overlays(), expected);
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
        vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)])],
        ValueAst::Lit(2),
    )]
    #[case::two_bonds(
        vec![
            (vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)]),
            (vec![AtomId(0), AtomId(2)], vec![ValueAst::Lit(1), ValueAst::Lit(1)]),
        ],
        ValueAst::Lit(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Undetermined, ValueAst::Lit(2)])],
        ValueAst::Undetermined,
    )]
    fn test_atom_view_multicenter_valence(
        #[case] bonds: Vec<(Vec<AtomId>, Vec<ValueAst>)>,
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
    fn test_atom_view_total_degree() {
        // Methane: 0 incident bonds in graph + implicit_h=4 + no multicenter.
        let molecule = mol!(r#"{:atoms ["C#h4"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_degree(),
            ValueAst::Lit(4),
        );
    }

    #[rstest]
    fn test_atom_view_total_degree_undetermined() {
        // implicit_hydrogens = Normal (placeholder) collapses to Undetermined.
        let molecule = mol!(r#"{:atoms ["C#h="] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_degree(),
            ValueAst::Undetermined,
        );
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
    #[case::implicit_normal_collapses(
        mol!(r#"{:atoms ["C#h="] :bonds []}"#),
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
    fn test_atom_view_total_valence_sum_of_terms() {
        // Methane with implicit_h=4: valence=0, implicit=4, aromatic=0,
        // multicenter=0 → total=4.
        let molecule = mol!(r#"{:atoms ["C#h4"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_valence(),
            ValueAst::Lit(4),
        );
    }

    #[rstest]
    fn test_atom_view_total_valence_implicit_normal_collapses() {
        let molecule = mol!(r#"{:atoms ["C#h="] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_valence(),
            ValueAst::Undetermined,
        );
    }

    #[rstest]
    fn test_atom_view_multicenter_degree() {
        // 3-atom multicenter bond: atom 0's multicenter_degree = co-participant
        // count = 2.
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
                MulticenterBondAst::new(vec![
                    ValueAst::Lit(2),
                    ValueAst::Lit(2),
                    ValueAst::Lit(2),
                ]),
            )],
            vec![],
            Constraints::default(),
        );
        assert_eq!(
            molecule.atom(AtomId(0)).multicenter_degree(),
            ValueAst::Lit(2),
        );
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
    #[case::ring_atom(AtomId(0), vec![6])]
    #[case::chain_atom(AtomId(6), vec![])]
    fn test_atom_view_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<usize>,
    ) {
        let sizes: Vec<_> = ring_with_chain.atom(atom).ring_size().collect();
        assert_eq!(sizes, expected);
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

    // --- BondViews ---

    #[rstest]
    fn test_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.bonds().count(), 3);
    }

    #[rstest]
    fn test_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.bonds().ids().collect::<Vec<_>>(),
            vec![BondId(0), BondId(1), BondId(2)],
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(BondId, [AtomId; 2], BondAst)> = molecule
            .bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), [AtomId(0), AtomId(1)], BondAst::from_order(1)),
                (BondId(1), [AtomId(1), AtomId(2)], BondAst::from_order(2)),
                (BondId(2), [AtomId(2), AtomId(3)], BondAst::from_order(1)),
            ],
        );
    }

    #[rstest]
    fn test_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.bonds().get(BondId(1));
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
        assert_eq!(*view.ast, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_views_index(molecule: MoleculeAst) {
        let bond: &BondAst = &molecule.bonds()[BondId(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    // --- BondView ---

    #[rstest]
    fn test_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule.bond(BondId(1)).atoms().map(|a| a.id).collect();
        assert_eq!(ids, vec![AtomId(1), AtomId(2)]);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), Some(AromaticSystemId(0)))]
    #[case::both_endpoints_aromatic_alt(BondId(1), Some(AromaticSystemId(0)))]
    #[case::one_endpoint_outside(BondId(2), None)]
    fn test_bond_view_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let id = molecule.bond(bond).aromatic_system().map(|v| v.id);
        assert_eq!(id, expected);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), true)]
    #[case::both_endpoints_aromatic_alt(BondId(1), true)]
    #[case::one_endpoint_outside(BondId(2), false)]
    fn test_bond_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.bond(bond).is_in_aromatic_system(), expected);
    }

    #[rstest]
    #[case::ring_bond_0_1(BondId(0), true)]
    #[case::ring_bond_5_0(BondId(5), true)]
    #[case::chain_bond_0_6(BondId(6), false)]
    fn test_bond_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(ring_with_chain.bond(bond).is_in_ring(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), true)]
    #[case::chain_bond(BondId(6), false)]
    fn test_bond_view_is_in_ring_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        assert_eq!(ring_with_chain.bond(bond).is_in_ring_from(&rings), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), 1)]
    #[case::chain_bond(BondId(6), 0)]
    fn test_bond_view_rings_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected_count: usize,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        let count = ring_with_chain.bond(bond).rings_from(&rings).count();
        assert_eq!(count, expected_count);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), ValueAst::Lit(1))]
    #[case::chain_bond(BondId(6), ValueAst::Lit(0))]
    fn test_bond_view_ring_count(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.bond(bond).ring_count(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), vec![6])]
    #[case::chain_bond(BondId(6), vec![])]
    fn test_bond_view_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Vec<usize>,
    ) {
        let sizes: Vec<_> = ring_with_chain.bond(bond).ring_size().collect();
        assert_eq!(sizes, expected);
    }

    // --- DativeBondViews ---

    #[rstest]
    fn test_dative_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bonds().count(), 1);
    }

    #[rstest]
    fn test_dative_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bonds().ids().collect::<Vec<_>>(),
            vec![DativeBondId(0)],
        );
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(DativeBondId, AtomId, DativeBondAst)> = molecule
            .dative_bonds()
            .iter()
            .map(|v| (v.id, v.acceptor_id, v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                DativeBondId(0),
                AtomId(3),
                DativeBondAst::from_order(1).with_acceptor_slot(1),
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.dative_bonds().get(DativeBondId(0));
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_index(molecule: MoleculeAst) {
        let dative: &DativeBondAst = &molecule.dative_bonds()[DativeBondId(0)];
        assert_eq!(dative.order, ValueAst::Lit(1));
    }

    // --- DativeBondView ---

    #[rstest]
    fn test_dative_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_donor_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor_id(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .donors()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2)]);
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor().id,
            AtomId(3),
        );
    }

    #[rstest]
    fn test_dative_bond_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).atom_count(), 2);
    }

    // --- AromaticSystemViews ---

    #[rstest]
    fn test_aromatic_system_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(AromaticSystemId, Vec<AtomId>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                AromaticSystemId(0),
                vec![AtomId(0), AtomId(1), AtomId(2)]
            )],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: MoleculeAst) {
        let view = molecule.aromatic_systems().get(AromaticSystemId(0));
        assert_eq!(view.id, AromaticSystemId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_index(molecule: MoleculeAst) {
        let _: &AromaticSystemAst = &molecule.aromatic_systems()[AromaticSystemId(0)];
    }

    // --- AromaticSystemView ---

    #[rstest]
    fn test_aromatic_system_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .bond_ids()
                .collect::<Vec<_>>(),
            vec![BondId(0), BondId(1)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: MoleculeAst) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![BondId(0), BondId(1)]);
    }

    #[rstest]
    fn test_aromatic_system_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).electron_count(),
            ValueAst::Lit(0),
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_system(AromaticSystemId(0)).atom_count(), 3);
    }

    #[rstest]
    fn test_aromatic_system_view_bond_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_system(AromaticSystemId(0)).bond_count(), 2);
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_aromatic_system_view_overlapping_atoms(
        molecule: MoleculeAst,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::one(vec![BondId(0)], vec![BondId(0)])]
    #[case::both(vec![BondId(0), BondId(1)], vec![BondId(0), BondId(1)])]
    #[case::other(vec![BondId(2)], vec![])]
    fn test_aromatic_system_view_overlapping_bonds(
        molecule: MoleculeAst,
        #[case] subset: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_bonds(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    fn test_aromatic_system_view_overlapping_rings(molecule: MoleculeAst) {
        let ids: Vec<usize> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_rings()
            .map(|r| r.len())
            .collect();
        assert_eq!(ids, Vec::<usize>::new());
    }

    // --- MulticenterBondViews ---

    #[rstest]
    fn test_multicenter_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bonds().ids().collect::<Vec<_>>(),
            vec![MulticenterBondId(0)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(MulticenterBondId, Vec<AtomId>)> = molecule
            .multicenter_bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                MulticenterBondId(0),
                vec![AtomId(0), AtomId(1), AtomId(2)],
            )],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.multicenter_bonds().get(MulticenterBondId(0));
        assert_eq!(view.id, MulticenterBondId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_index(molecule: MoleculeAst) {
        let _: &MulticenterBondAst = &molecule.multicenter_bonds()[MulticenterBondId(0)];
    }

    // --- MulticenterBondView ---

    #[rstest]
    fn test_multicenter_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_multicenter_bond_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bond(MulticenterBondId(0)).electron_count(),
            ValueAst::Lit(0),
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bond(MulticenterBondId(0)).atom_count(),
            3,
        );
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_multicenter_bond_view_overlapping_atoms(
        molecule: MoleculeAst,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    // --- NoncovalentBondViews ---

    #[rstest]
    fn test_noncovalent_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.noncovalent_bonds().count(), 1);
    }

    #[rstest]
    fn test_noncovalent_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.noncovalent_bonds().ids().collect::<Vec<_>>(),
            vec![NoncovalentBondId(0)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(NoncovalentBondId, [AtomId; 2], NoncovalentBondAst)> = molecule
            .noncovalent_bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(3)],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.noncovalent_bonds().get(NoncovalentBondId(0));
        assert_eq!(view.id, NoncovalentBondId(0));
        assert_eq!(view.atom_ids(), [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_index(molecule: MoleculeAst) {
        let _: &NoncovalentBondAst = &molecule.noncovalent_bonds()[NoncovalentBondId(0)];
    }

    // --- NoncovalentBondView ---

    #[rstest]
    fn test_noncovalent_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondId(0))
                .atom_ids(),
            [AtomId(0), AtomId(3)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_view_atoms(molecule: MoleculeAst) {
        let ids = molecule
            .noncovalent_bond(NoncovalentBondId(0))
            .atoms()
            .map(|v| v.id);
        assert_eq!(ids, [AtomId(0), AtomId(3)]);
    }

    // --- NeighborView ---

    #[rstest]
    fn test_neighbor_view_fields(molecule: MoleculeAst) {
        let collected: Vec<(BondId, AtomId, BondAst)> = molecule
            .neighbors(AtomId(2))
            .map(|n| (n.bond, n.atom, n.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(1), AtomId(1), BondAst::from_order(2)),
                (BondId(2), AtomId(3), BondAst::from_order(1)),
            ],
        );
    }
}
