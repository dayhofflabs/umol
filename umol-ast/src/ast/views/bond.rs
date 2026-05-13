//! Bond views: `BondViews` namespace, `BondView` / `BondViewMut` AST bundles,
//! `BondBuilderView` / `BondBuilderViewMut` builder bundles.

use std::ops::Index;

use umol_graph_core::{EdgeId, NodeId};

use super::super::bond::BondAst;
use super::super::constraint::BondConstraints;
use super::super::idx::{AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::super::rings::{RingSet, RingView};
use super::super::spin::SpinStateAst;
use super::super::value::ValueAst;
use super::aromatic_system::AromaticSystemView;
use super::atom::AtomView;

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    molecule: &'a MoleculeAst,
    bonds: &'a [BondAst],
}

impl<'a> BondViews<'a> {
    pub(in crate::ast) fn new(molecule: &'a MoleculeAst, bonds: &'a [BondAst]) -> Self {
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
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
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
    pub(in crate::ast) fn new(id: BondId, atoms: [AtomId; 2], ast: &'a mut BondAst) -> Self {
        Self { id, atoms, ast }
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

// Builder-scope view bundles for bonds.

pub struct BondBuilderView<'a> {
    pub id: BondId,
    pub ast: &'a BondAst,
    pub atoms: [AtomId; 2],
}

pub struct BondBuilderViewMut<'a> {
    pub id: BondId,
    pub ast: &'a mut BondAst,
    pub atoms: [AtomId; 2],
}
