//! Dative bond views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, VarRelationSet};

use super::super::constraint::DativeBondConstraints;
use super::super::dative::DativeBondAst;
use super::super::idx::{AtomId, DativeBondId};
use super::super::molecule::MoleculeAst;
use super::super::value::ValueAst;
use super::atom::AtomView;

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(in crate::ast) fn new(
        molecule: &'a MoleculeAst,
        set: &'a VarRelationSet<DativeBondAst>,
    ) -> Self {
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
            DativeBondView {
                id,
                ast,
                acceptor_id,
                atoms,
                molecule,
            }
        })
    }

    /// ID of the dative bond whose participant set equals `atoms`, if any.
    pub fn connecting_id(&self, atoms: impl IntoIterator<Item = AtomId>) -> Option<DativeBondId> {
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
            .filter(|&rid| self.set.participants(rid).iter().all(|p| set.contains(p)))
            .map(DativeBondId::from)
            .collect()
    }

    /// Views of dative bonds whose participants all lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<DativeBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
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

// Builder-scope view bundles for dative bonds.

pub struct DativeBondBuilderView<'a> {
    pub id: DativeBondId,
    pub ast: &'a DativeBondAst,
    pub(crate) atoms: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct DativeBondBuilderViewMut<'a> {
    pub id: DativeBondId,
    pub ast: &'a mut DativeBondAst,
    pub(crate) atoms: &'a [NodeId],
    pub acceptor_id: AtomId,
}

impl<'a> DativeBondBuilderViewMut<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}
