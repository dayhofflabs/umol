//! Aromatic system views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, VarRelationSet};

use super::super::aromatic::AromaticSystemAst;
use super::super::constraint::AromaticSystemConstraints;
use super::super::idx::{AromaticSystemId, AtomId, BondId};
use super::super::molecule::MoleculeAst;
use super::super::rings::RingView;
use super::super::spin::SpinStateAst;
use super::super::value::ValueAst;
use super::atom::AtomView;
use super::bond::BondView;

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(in crate::ast) fn new(
        molecule: &'a MoleculeAst,
        set: &'a VarRelationSet<AromaticSystemAst>,
    ) -> Self {
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
            AromaticSystemView {
                id,
                ast: set.data(rid),
                atoms: set.participants(rid),
                molecule,
            }
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
            let parts: HashSet<AtomId> = self
                .set
                .participants(RelationId::from(id))
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
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
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

// Builder-scope view bundles for aromatic systems.

pub struct AromaticSystemBuilderView<'a> {
    pub id: AromaticSystemId,
    pub ast: &'a AromaticSystemAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> AromaticSystemBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct AromaticSystemBuilderViewMut<'a> {
    pub id: AromaticSystemId,
    pub ast: &'a mut AromaticSystemAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> AromaticSystemBuilderViewMut<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}
