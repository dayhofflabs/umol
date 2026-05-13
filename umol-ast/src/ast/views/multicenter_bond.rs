//! Multicenter bond views.

use std::collections::HashSet;
use std::ops::Index;

use umol_graph_core::{NodeId, RelationId, VarRelationSet};

use super::super::constraint::MulticenterBondConstraints;
use super::super::idx::{AtomId, MulticenterBondId};
use super::super::molecule::MoleculeAst;
use super::super::multicenter::MulticenterBondAst;
use super::super::spin::SpinStateAst;
use super::super::value::ValueAst;
use super::atom::AtomView;

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(in crate::ast) fn new(
        molecule: &'a MoleculeAst,
        set: &'a VarRelationSet<MulticenterBondAst>,
    ) -> Self {
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
            MulticenterBondView {
                id,
                ast: set.data(rid),
                atoms: set.participants(rid),
                molecule,
            }
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
            let parts: HashSet<AtomId> = self
                .set
                .participants(RelationId::from(id))
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
        self.induced_ids(atoms)
            .into_iter()
            .map(|id| self.get(id))
            .collect()
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

// Builder-scope view bundles for multicenter bonds.

pub struct MulticenterBondBuilderView<'a> {
    pub id: MulticenterBondId,
    pub ast: &'a MulticenterBondAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> MulticenterBondBuilderView<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

pub struct MulticenterBondBuilderViewMut<'a> {
    pub id: MulticenterBondId,
    pub ast: &'a mut MulticenterBondAst,
    pub(crate) atoms: &'a [NodeId],
}

impl<'a> MulticenterBondBuilderViewMut<'a> {
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}
