//! Noncovalent bond views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{FixedRelationSet, NodeId, RelationId};

use super::super::constraint::NoncovalentBondConstraints;
use super::super::idx::{AtomId, NoncovalentBondId};
use super::super::molecule::MoleculeAst;
use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use super::atom::AtomView;

/// Namespace accessor for noncovalent-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(in crate::ast) fn new(
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
            .filter(|&rid| self.set.participants(rid).iter().all(|p| set.contains(p)))
            .map(NoncovalentBondId::from)
            .collect()
    }

    /// Views of noncovalent bonds whose endpoints both lie in `atoms`.
    pub fn induced(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondView<'a>> {
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
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

// Builder-scope view bundles for noncovalent bonds.

pub struct NoncovalentBondBuilderView<'a> {
    pub id: NoncovalentBondId,
    pub ast: &'a NoncovalentBondAst,
    pub atoms: [AtomId; 2],
}

pub struct NoncovalentBondBuilderViewMut<'a> {
    pub id: NoncovalentBondId,
    pub ast: &'a mut NoncovalentBondAst,
    pub atoms: [AtomId; 2],
}
