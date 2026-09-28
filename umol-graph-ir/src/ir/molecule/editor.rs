//! Structural editing for `Molecule`. The molecule itself only allows attribute
//! mutation; structural change (add atoms/bonds/relations, remove anything)
//! goes through `MoleculeEditor`.
//!
//! The editor owns a molecule draft. Mutable access uses its copy-on-write storage;
//! publication checks molecule integrity.

use std::mem;
use std::sync::Arc;

pub use transact::{Transaction, TransactionError};
use umol_graph_core::{Compaction, Correspondence, EdgeId, Graph, GraphCompaction, NodeId};
use umol_perm::{DynPermutation, Permutation};

use super::super::aromatic::{AromaticSystemForm, AromaticSystems};
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::compact::MoleculeCompaction;
use super::super::constraint::{Constraint, Constraints};
use super::super::correspondence::MoleculeCorrespondence;
use super::super::dative::{DativeBondForm, DativeBonds};
use super::super::edit::{
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond,
};
use super::super::entity::EntityKind;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::multicenter::{MulticenterBondForm, MulticenterBonds};
use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBonds};
use super::super::stereo::{StereoAtomForm, StereoAtoms, StereoBondForm, StereoBonds};
use super::super::traits::{FrameTransport, Normalize};
use super::super::view::{
    AromaticSystemEditorView, AromaticSystemEditorViewMut, AtomEditorView, AtomEditorViewMut,
    BondEditorView, BondEditorViewMut, DativeBondEditorView, DativeBondEditorViewMut,
    MulticenterBondEditorView, MulticenterBondEditorViewMut, NoncovalentBondEditorView,
    NoncovalentBondEditorViewMut, StereoAtomEditorView, StereoAtomEditorViewMut,
    StereoBondEditorView, StereoBondEditorViewMut,
};
use super::{Molecule, MoleculeIntegrityError};

mod transact;

/// Editor for structural and attribute changes to a `Molecule`.
///
/// Publication checks molecule integrity. Removal compacts surviving entity ids.
///
/// The session correspondence maps the editor's initial id spaces to its current ones.
/// Tracking stores only id pairs and counts, not a source molecule. Additions are right-unmatched;
/// removals discard pairs. Restoration expands the id spaces without recreating discarded pairs.
/// Attribute-only changes preserve the pairings.
///
/// Bulk additions move attributes and copy supplied atom/ligand slices into storage.
/// Mutation completes before return; the exact-size id iterators retain no receiver
/// or input borrow. These are direct mutations, without an undo journal or integrity
/// checks at addition.
///
/// # Semantic properties
///
/// Bulk additions preserve existing ids, attributes, and overlay frames. New ids are
/// contiguous in input order. Empty batches leave storage unchanged. Splitting a
/// batch into consecutive additions gives the same stored result. These properties
/// are exercised through publication in `tests/property/edit.rs`.
///
/// The session correspondence composes the id changes since editor creation. Discarding the
/// correspondence from a tracked publication gives the same molecule or integrity error as its
/// plain counterpart. Repeated snapshots without intervening edits are equal; later edits do not
/// change an earlier snapshot or its correspondence.
#[derive(Clone)]
pub struct MoleculeEditor {
    molecule: Molecule,
    correspondence: MoleculeCorrespondence,
}

impl MoleculeEditor {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        graph: Graph,
        atoms: Arc<Vec<AtomForm>>,
        bonds: Arc<Vec<BondForm>>,
        dative_bonds: DativeBonds,
        aromatic_systems: AromaticSystems,
        multicenter_bonds: MulticenterBonds,
        noncovalent_bonds: NoncovalentBonds,
        stereo_atoms: StereoAtoms,
        stereo_bonds: StereoBonds,
        constraints: Constraints,
    ) -> Self {
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::identity(atoms.len()),
            Correspondence::identity(bonds.len()),
            Correspondence::identity(dative_bonds.count()),
            Correspondence::identity(aromatic_systems.count()),
            Correspondence::identity(multicenter_bonds.count()),
            Correspondence::identity(noncovalent_bonds.count()),
            Correspondence::identity(stereo_atoms.count()),
            Correspondence::identity(stereo_bonds.count()),
        );
        Self {
            molecule: Molecule {
                graph,
                atoms,
                bonds,
                dative_bonds,
                aromatic_systems,
                multicenter_bonds,
                noncovalent_bonds,
                stereo_atoms,
                stereo_bonds,
                constraints,
            },
            correspondence,
        }
    }

    /// Append an atom directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. Use `transact` for checked
    /// atomic edits with rollback or `apply` for consuming application without an undo journal.
    pub fn add_atom(&mut self, atom: AtomForm) -> AtomId {
        let id = self.molecule.graph.add_node();
        Arc::make_mut(&mut self.molecule.atoms).push(atom);
        self.correspondence.extend_right(EntityKind::Atom, 1);
        AtomId::from(id)
    }

    /// Append atoms in input order and return their ids.
    pub fn add_atoms(
        &mut self,
        atoms: Vec<AtomForm>,
    ) -> impl ExactSizeIterator<Item = AtomId> + use<> {
        let ids = self.molecule.add_atoms(atoms);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::Atom, ids.len());
        }
        ids
    }

    /// Append a localized bond directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. It
    /// assumes `first` and `second` are valid atom ids in the current dense layout.
    pub fn add_bond(&mut self, first: AtomId, second: AtomId, bond: BondForm) -> BondId {
        let id = self
            .molecule
            .graph
            .add_edge(NodeId::from(first), NodeId::from(second));
        Arc::make_mut(&mut self.molecule.bonds).push(bond);
        self.correspondence.extend_right(EntityKind::Bond, 1);
        BondId::from(id)
    }

    /// Append localized bonds in input order and return their ids.
    ///
    /// Bond endpoints are stored in increasing atom-id order.
    ///
    /// # Panics
    ///
    /// Panics if an endpoint is outside the current atom space.
    pub fn add_bonds(
        &mut self,
        bonds: Vec<([AtomId; 2], BondForm)>,
    ) -> impl ExactSizeIterator<Item = BondId> + use<> {
        let ids = self.molecule.add_bonds(bonds);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::Bond, ids.len());
        }
        ids
    }

    /// Append a dative-bond overlay directly to the editor. The acceptor is factor 1; the donors
    /// are factor 2. Their supplied order is retained as the stored frame.
    pub fn add_dative_bond(
        &mut self,
        donors: &[AtomId],
        acceptor: AtomId,
        bond: DativeBondForm,
    ) -> DativeBondId {
        let id = self.molecule.dative_bonds.add(donors, acceptor, bond);
        self.correspondence.extend_right(EntityKind::DativeBond, 1);
        id
    }

    /// Append dative bonds in input order and return their ids.
    ///
    /// Donor and acceptor ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_dative_bonds(
        &mut self,
        entries: Vec<(&[AtomId], AtomId, DativeBondForm)>,
    ) -> impl ExactSizeIterator<Item = DativeBondId> + use<> {
        let ids = self.molecule.add_dative_bonds(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::DativeBond, ids.len());
        }
        ids
    }

    /// Append an aromatic-system overlay directly to the editor.
    pub fn add_aromatic_system(
        &mut self,
        atoms: &[AtomId],
        data: AromaticSystemForm,
    ) -> AromaticSystemId {
        let id = self.molecule.aromatic_systems.add(atoms, data);
        self.correspondence
            .extend_right(EntityKind::AromaticSystem, 1);
        id
    }

    /// Append aromatic systems in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_aromatic_systems(
        &mut self,
        entries: Vec<(&[AtomId], AromaticSystemForm)>,
    ) -> impl ExactSizeIterator<Item = AromaticSystemId> + use<> {
        let ids = self.molecule.add_aromatic_systems(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::AromaticSystem, ids.len());
        }
        ids
    }

    /// Append a multicenter-bond overlay directly to the editor.
    pub fn add_multicenter_bond(
        &mut self,
        atoms: &[AtomId],
        data: MulticenterBondForm,
    ) -> MulticenterBondId {
        let id = self.molecule.multicenter_bonds.add(atoms, data);
        self.correspondence
            .extend_right(EntityKind::MulticenterBond, 1);
        id
    }

    /// Append multicenter bonds in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_multicenter_bonds(
        &mut self,
        entries: Vec<(&[AtomId], MulticenterBondForm)>,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + use<> {
        let ids = self.molecule.add_multicenter_bonds(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::MulticenterBond, ids.len());
        }
        ids
    }

    /// Append a noncovalent-bond overlay directly to the editor.
    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomId; 2],
        bond: NoncovalentBondForm,
    ) -> NoncovalentBondId {
        let id = self.molecule.noncovalent_bonds.add(ends, bond);
        self.correspondence
            .extend_right(EntityKind::NoncovalentBond, 1);
        id
    }

    /// Append noncovalent bonds in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_noncovalent_bonds(
        &mut self,
        entries: Vec<([AtomId; 2], NoncovalentBondForm)>,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondId> + use<> {
        let ids = self.molecule.add_noncovalent_bonds(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::NoncovalentBond, ids.len());
        }
        ids
    }

    /// Append a stereo-atom overlay directly to the editor.
    pub fn add_stereo_atom(
        &mut self,
        site: AtomId,
        ligands: &[StereoLigand],
        attributes: StereoAtomForm,
    ) -> StereoAtomId {
        let id = self.molecule.stereo_atoms.add(site, ligands, attributes);
        self.correspondence.extend_right(EntityKind::StereoAtom, 1);
        id
    }

    /// Append stereo atoms in input order and return their ids.
    ///
    /// Site and ligand ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_stereo_atoms(
        &mut self,
        entries: Vec<(AtomId, &[StereoLigand], StereoAtomForm)>,
    ) -> impl ExactSizeIterator<Item = StereoAtomId> + use<> {
        let ids = self.molecule.add_stereo_atoms(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::StereoAtom, ids.len());
        }
        ids
    }

    /// Append a stereo-bond overlay directly to the editor.
    pub fn add_stereo_bond(
        &mut self,
        site: BondId,
        ligands: &[StereoLigand],
        attributes: StereoBondForm,
    ) -> StereoBondId {
        let id = self.molecule.stereo_bonds.add(site, ligands, attributes);
        self.correspondence.extend_right(EntityKind::StereoBond, 1);
        id
    }

    /// Append stereo bonds in input order and return their ids.
    ///
    /// Site and ligand ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_stereo_bonds(
        &mut self,
        entries: Vec<(BondId, &[StereoLigand], StereoBondForm)>,
    ) -> impl ExactSizeIterator<Item = StereoBondId> + use<> {
        let ids = self.molecule.add_stereo_bonds(entries);
        if ids.len() != 0 {
            self.correspondence
                .extend_right(EntityKind::StereoBond, ids.len());
        }
        ids
    }

    /// Add a molecule-level constraint (molecule-scope predicate or
    /// combinator). Unconditional per-entity constraints belong inline on the
    /// entity — use `atom_mut(id).attributes_mut().constraints.set(c)` etc.
    pub fn push_constraint(&mut self, c: Constraint) {
        self.molecule.constraints.push(c);
    }

    // -- Attribute access -----------------------------------------------------
    //
    // Mutable views edit entity data in place. Structural add/remove stays on
    // the editor itself because dense removal can compact many unrelated ids.

    pub fn atom(&self, id: AtomId) -> AtomEditorView<'_> {
        AtomEditorView::new(id, &self.molecule.atoms[id.index()])
    }

    pub fn atom_mut(&mut self, id: AtomId) -> AtomEditorViewMut<'_> {
        self.molecule.atom_view_mut(id)
    }

    pub fn bond(&self, id: BondId) -> BondEditorView<'_> {
        let endpoints = self.molecule.graph.edge_endpoints(EdgeId::from(id));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        BondEditorView::new(id, atoms, &self.molecule.bonds[id.index()])
    }

    pub fn bond_mut(&mut self, id: BondId) -> BondEditorViewMut<'_> {
        self.molecule.bond_view_mut(id)
    }

    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondEditorView<'_> {
        assert!(
            self.molecule.dative_bonds.contains(id),
            "invalid dative bond id"
        );
        DativeBondEditorView::new(&self.molecule.dative_bonds, id)
    }

    pub fn dative_bond_mut(&mut self, id: DativeBondId) -> DativeBondEditorViewMut<'_> {
        self.molecule.dative_bond_view_mut(id)
    }

    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemEditorView<'_> {
        assert!(
            self.molecule.aromatic_systems.contains(id),
            "invalid aromatic system id"
        );
        AromaticSystemEditorView::new(&self.molecule.aromatic_systems, id)
    }

    pub fn aromatic_system_mut(&mut self, id: AromaticSystemId) -> AromaticSystemEditorViewMut<'_> {
        self.molecule.aromatic_system_view_mut(id)
    }

    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondEditorView<'_> {
        assert!(
            self.molecule.multicenter_bonds.contains(id),
            "invalid multicenter bond id"
        );
        MulticenterBondEditorView::new(&self.molecule.multicenter_bonds, id)
    }

    pub fn multicenter_bond_mut(
        &mut self,
        id: MulticenterBondId,
    ) -> MulticenterBondEditorViewMut<'_> {
        self.molecule.multicenter_bond_view_mut(id)
    }

    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondEditorView<'_> {
        assert!(
            self.molecule.noncovalent_bonds.contains(id),
            "invalid noncovalent bond id"
        );
        NoncovalentBondEditorView::new(&self.molecule.noncovalent_bonds, id)
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        id: NoncovalentBondId,
    ) -> NoncovalentBondEditorViewMut<'_> {
        self.molecule.noncovalent_bond_view_mut(id)
    }

    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomEditorView<'_> {
        let set = &self.molecule.stereo_atoms;
        assert!(set.contains(id), "invalid stereo atom id");
        StereoAtomEditorView::new(set, id)
    }

    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondEditorView<'_> {
        let set = &self.molecule.stereo_bonds;
        assert!(set.contains(id), "invalid stereo bond id");
        StereoBondEditorView::new(set, id)
    }

    /// `true` iff noncovalent bond `id` structurally equals `(atoms, attributes)` — participants (unordered)
    /// and `attributes` up to normal form, `attributes` reindexed into the stored participant frame.
    pub(crate) fn noncovalent_bond_equiv(
        &self,
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        attributes: &NoncovalentBondForm,
    ) -> bool {
        let set = &self.molecule.noncovalent_bonds;
        let stored = set.atoms(id);
        set.is_coincident(id, atoms[0], atoms[1])
            && DynPermutation::between(&atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff aromatic system `id` structurally equals `(atoms, attributes)`.
    pub(crate) fn aromatic_system_equiv(
        &self,
        id: AromaticSystemId,
        atoms: &[AtomId],
        attributes: &AromaticSystemForm,
    ) -> bool {
        let set = &self.molecule.aromatic_systems;
        let stored: Vec<AtomId> = set.atoms(id).collect();
        set.is_coincident(id, atoms)
            && DynPermutation::between(atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff multicenter bond `id` structurally equals `(atoms, attributes)`.
    pub(crate) fn multicenter_bond_equiv(
        &self,
        id: MulticenterBondId,
        atoms: &[AtomId],
        attributes: &MulticenterBondForm,
    ) -> bool {
        let set = &self.molecule.multicenter_bonds;
        let stored: Vec<AtomId> = set.atoms(id).collect();
        set.is_coincident(id, atoms)
            && DynPermutation::between(atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff dative bond `id` structurally equals `(acceptor, donors, attributes)` — the acceptor
    /// (ordered, single) and donors (unordered) factors and `attributes` up to normal form.
    pub(crate) fn dative_bond_equiv(
        &self,
        id: DativeBondId,
        acceptor: AtomId,
        donors: &[AtomId],
        attributes: &DativeBondForm,
    ) -> bool {
        let set = &self.molecule.dative_bonds;
        let stored: Vec<AtomId> = set.donors(id).collect();
        set.is_coincident(id, acceptor, donors)
            && DynPermutation::between(donors, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff stereo atom `id` structurally equals `(site, ligands, attributes)`.
    pub(crate) fn stereo_atom_equiv(
        &self,
        id: StereoAtomId,
        site: AtomId,
        ligands: &[StereoLigand],
        attributes: &StereoAtomForm,
    ) -> bool {
        let set = &self.molecule.stereo_atoms;
        let stored = set.ligands(id);
        set.site(id) == site
            && Permutation::between(ligands, stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff stereo bond `id` structurally equals `(site, ligands, attributes)`.
    pub(crate) fn stereo_bond_equiv(
        &self,
        id: StereoBondId,
        site: BondId,
        ligands: &[StereoLigand],
        attributes: &StereoBondForm,
    ) -> bool {
        let set = &self.molecule.stereo_bonds;
        let stored = set.ligands(id);
        set.site(id) == site
            && Permutation::between(ligands, stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    pub fn stereo_atom_mut(&mut self, id: StereoAtomId) -> StereoAtomEditorViewMut<'_> {
        self.molecule.stereo_atom_view_mut(id)
    }

    pub fn stereo_bond_mut(&mut self, id: StereoBondId) -> StereoBondEditorViewMut<'_> {
        self.molecule.stereo_bond_view_mut(id)
    }

    pub fn constraints(&self) -> &Constraints {
        self.molecule.constraints()
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        self.molecule.constraints_mut()
    }

    pub fn atom_count(&self) -> usize {
        self.molecule.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.molecule.bonds.len()
    }

    pub fn dative_bond_count(&self) -> usize {
        self.molecule.dative_bonds.count()
    }

    pub fn aromatic_system_count(&self) -> usize {
        self.molecule.aromatic_systems.count()
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.molecule.multicenter_bonds.count()
    }

    pub fn noncovalent_bond_count(&self) -> usize {
        self.molecule.noncovalent_bonds.count()
    }

    pub fn stereo_atom_count(&self) -> usize {
        self.molecule.stereo_atoms.count()
    }

    pub fn stereo_bond_count(&self) -> usize {
        self.molecule.stereo_bonds.count()
    }

    // -- Relation removal -----------------------------------------------------

    /// Remove dative-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_dative_bonds(&mut self, ids: &[DativeBondId]) {
        self.tracked_remove_dative_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_dative_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_dative_bonds(&mut self, ids: &[DativeBondId]) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            self.molecule.dative_bonds.tracked_remove(ids),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    /// Remove aromatic-system overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_aromatic_systems(&mut self, ids: &[AromaticSystemId]) {
        self.tracked_remove_aromatic_systems(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_aromatic_systems`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_aromatic_systems(
        &mut self,
        ids: &[AromaticSystemId],
    ) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            self.molecule.aromatic_systems.tracked_remove(ids),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    /// Remove multicenter-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_multicenter_bonds(&mut self, ids: &[MulticenterBondId]) {
        self.tracked_remove_multicenter_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_multicenter_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_multicenter_bonds(
        &mut self,
        ids: &[MulticenterBondId],
    ) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            self.molecule.multicenter_bonds.tracked_remove(ids),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    /// Remove noncovalent-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_noncovalent_bonds(&mut self, ids: &[NoncovalentBondId]) {
        self.tracked_remove_noncovalent_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_noncovalent_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_noncovalent_bonds(
        &mut self,
        ids: &[NoncovalentBondId],
    ) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            self.molecule.noncovalent_bonds.tracked_remove(ids),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    /// Remove stereo-atom overlays directly from the editor.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_atoms(&mut self, ids: &[StereoAtomId]) {
        self.tracked_remove_stereo_atoms(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_stereo_atoms`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_stereo_atoms(&mut self, ids: &[StereoAtomId]) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            self.molecule.stereo_atoms.tracked_remove(ids),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    /// Remove stereo-bond overlays directly from the editor.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_bonds(&mut self, ids: &[StereoBondId]) {
        self.tracked_remove_stereo_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_stereo_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_stereo_bonds(&mut self, ids: &[StereoBondId]) -> MoleculeCompaction {
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            self.molecule.stereo_bonds.tracked_remove(ids),
        );
        self.molecule.constraints.compact(&compaction);
        self.correspondence
            .compact_right(&compaction)
            .expect("removal compaction describes the editor's current id spaces");
        compaction
    }

    // -- Topological removal --------------------------------------------------

    /// Remove topology directly from the editor, cascading dependent relations and constraints.
    ///
    /// This is the low-level dense topology-removal primitive. It removes the
    /// requested atoms and bonds, cascades relations whose participants were
    /// removed, and compacts molecule-level constraints. It does not build rollback
    /// data; checked transactions capture the removed payloads before calling
    /// this method.
    pub fn remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) {
        self.tracked_remove(atoms, bonds);
    }

    /// Remove topology and return the source-to-result compaction for all eight entity kinds.
    ///
    /// Leaves the same state as [`Self::remove`], including cascading relation and constraint
    /// removal. Every component retains the source count from before removal.
    pub fn tracked_remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) -> MoleculeCompaction {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let compaction = self.molecule.graph.tracked_remove_cascading(&nodes, &edges);

        let new_atoms = compaction.nodes().compact_vec(&self.molecule.atoms);
        let new_bonds = compaction.edges().compact_vec(&self.molecule.bonds);
        self.molecule.atoms = Arc::new(new_atoms);
        self.molecule.bonds = Arc::new(new_bonds);

        let (dative_bonds, removed_dative_bonds) =
            self.molecule.dative_bonds.tracked_compact(&compaction);
        self.molecule.dative_bonds = dative_bonds;
        let (aromatic_systems, removed_aromatic_systems) =
            self.molecule.aromatic_systems.tracked_compact(&compaction);
        self.molecule.aromatic_systems = aromatic_systems;
        let (multicenter_bonds, removed_multicenter_bonds) =
            self.molecule.multicenter_bonds.tracked_compact(&compaction);
        self.molecule.multicenter_bonds = multicenter_bonds;
        let (noncovalent_bonds, removed_noncovalent_bonds) =
            self.molecule.noncovalent_bonds.tracked_compact(&compaction);
        self.molecule.noncovalent_bonds = noncovalent_bonds;
        let (stereo_atoms, removed_stereo_atoms) =
            self.molecule.stereo_atoms.tracked_compact(&compaction);
        self.molecule.stereo_atoms = stereo_atoms;
        let (stereo_bonds, removed_stereo_bonds) =
            self.molecule.stereo_bonds.tracked_compact(&compaction);
        self.molecule.stereo_bonds = stereo_bonds;

        let id_compaction = MoleculeCompaction::new(
            compaction,
            removed_dative_bonds,
            removed_aromatic_systems,
            removed_multicenter_bonds,
            removed_noncovalent_bonds,
            removed_stereo_atoms,
            removed_stereo_bonds,
        );
        self.molecule.constraints.compact(&id_compaction);
        self.correspondence
            .compact_right(&id_compaction)
            .expect("removal compaction describes the editor's current id spaces");
        id_compaction
    }

    // -- Undo of additions ----------------------------------------------------

    fn remove_added_topology(&mut self, atoms: &[AddedAtom], bonds: &[AddedBond]) {
        let atom_ids: Vec<AtomId> = atoms.iter().map(|a| a.id).collect();
        let bond_ids: Vec<BondId> = bonds.iter().map(|b| b.id).collect();
        self.remove(&atom_ids, &bond_ids);
    }

    fn remove_added_dative_bond(&mut self, added: &AddedDativeBond) {
        self.remove_dative_bonds(&[added.id]);
    }

    fn remove_added_aromatic_system(&mut self, added: &AddedAromaticSystem) {
        self.remove_aromatic_systems(&[added.id]);
    }

    fn remove_added_multicenter_bond(&mut self, added: &AddedMulticenterBond) {
        self.remove_multicenter_bonds(&[added.id]);
    }

    fn remove_added_noncovalent_bond(&mut self, added: &AddedNoncovalentBond) {
        self.remove_noncovalent_bonds(&[added.id]);
    }

    fn remove_added_stereo_atom(&mut self, added: &AddedStereoAtom) {
        self.remove_stereo_atoms(&[added.id]);
    }

    fn remove_added_stereo_bond(&mut self, added: &AddedStereoBond) {
        self.remove_stereo_bonds(&[added.id]);
    }

    /// Materialize the editor's current state without consuming it, after checking molecule
    /// integrity.
    ///
    /// Subsequent editor changes are independent of the returned immutable snapshot.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError`] when the transient editor state cannot be published as a
    /// molecule.
    pub fn snapshot(&self) -> Result<Molecule, MoleculeIntegrityError> {
        self.molecule.check_integrity()?;
        Ok(self.molecule.clone())
    }

    /// Publish a reusable snapshot and the initial-to-current session correspondence.
    ///
    /// Subsequent edits do not change either returned value.
    ///
    /// # Errors
    ///
    /// Returns the same integrity error as [`Self::snapshot`].
    pub fn tracked_snapshot(
        &self,
    ) -> Result<(Molecule, MoleculeCorrespondence), MoleculeIntegrityError> {
        Ok((self.snapshot()?, self.correspondence.clone()))
    }

    /// Publish the editor's current state after checking molecule integrity.
    pub fn try_build(self) -> Result<Molecule, MoleculeIntegrityError> {
        self.molecule.check_integrity()?;
        Ok(self.molecule)
    }

    /// Consume the editor, publishing its molecule and initial-to-current session correspondence.
    ///
    /// Moves the accumulated id-pair vectors into the result without copying them.
    ///
    /// # Errors
    ///
    /// Returns the same integrity error as [`Self::try_build`]. The editor is consumed on failure.
    pub fn try_tracked_build(
        mut self,
    ) -> Result<(Molecule, MoleculeCorrespondence), MoleculeIntegrityError> {
        let correspondence =
            mem::replace(&mut self.correspondence, MoleculeCorrespondence::empty());
        Ok((self.try_build()?, correspondence))
    }

    /// Publish editor state whose molecule integrity is established by the producer.
    ///
    /// # Panics
    ///
    /// Panics when the editor does not contain a representation-integral molecule. Use
    /// [`Self::try_build`] for independently assembled or potentially conflicting edits.
    pub fn build(self) -> Molecule {
        self.try_build()
            .unwrap_or_else(|error| panic!("invalid molecule editor state: {error}"))
    }

    /// Consume an integrity-established editor, returning its molecule and session correspondence.
    ///
    /// # Panics
    ///
    /// Panics on the same integrity failure as [`Self::build`]. Use [`Self::try_tracked_build`]
    /// when the edits do not establish molecule integrity.
    pub fn tracked_build(self) -> (Molecule, MoleculeCorrespondence) {
        self.try_tracked_build()
            .unwrap_or_else(|error| panic!("invalid molecule editor state: {error}"))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_perm::Permutation;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::MoleculeConstraint;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::{RemovedAtom, RemovedBond, RemovedDativeBond};
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::molecule::MoleculeEntries;
    use crate::ir::noncovalent::NoncovalentBondKind;
    use crate::ir::num::NumForm;
    use crate::ir::stereo::StereoKind;
    use crate::mol_dsl;

    #[fixture]
    fn triatomic() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_atom(AtomForm::from_element(Element::O));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
        b
    }

    #[fixture]
    fn addition_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(2), AtomId(4), BondForm::from_order(2)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(4), AtomId(0), BondForm::from_order(2)),
            ],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_molecule_editor_add_atoms(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = editor.add_atoms(vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ]);

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(AtomId(6)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.atom_count(), 8);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.atoms.extend([
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![AtomId(7)]);
    }

    #[rstest]
    fn test_molecule_editor_add_atoms_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_atoms(vec![]);

        assert!(Arc::ptr_eq(&triatomic.molecule.atoms, &original.atoms));

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_bonds(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = editor.add_bonds(vec![
            ([AtomId(5), AtomId(1)], BondForm::from_order(3)),
            ([AtomId(3), AtomId(1)], BondForm::from_order(1)),
        ]);

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(BondId(6)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.bond_count(), 8);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.bonds.extend([
            (AtomId(5), AtomId(1), BondForm::from_order(3)),
            (AtomId(3), AtomId(1), BondForm::from_order(1)),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![BondId(7)]);
    }

    #[rstest]
    fn test_molecule_editor_add_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_bonds(vec![]);

        assert!(Arc::ptr_eq(&triatomic.molecule.bonds, &original.bonds));

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_dative_bond(mut triatomic: MoleculeEditor) {
        let attributes = DativeBondForm::from_order(1);
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            dative: vec![
                (vec![AtomId(0)], AtomId(1), attributes.clone()),
                (vec![AtomId(2)], AtomId(1), attributes.clone()),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_dative_bond(&[AtomId(0)], AtomId(1), attributes.clone());
        let second = triatomic.add_dative_bond(&[AtomId(2)], AtomId(1), attributes.clone());
        assert_eq!((first, second), (DativeBondId(0), DativeBondId(1)));
        assert_eq!(
            triatomic
                .molecule
                .dative_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![DativeBondId(0), DativeBondId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_dative_bonds(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .dative_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![DativeBondId(0)]
        );
        let mut remaining = entries;
        remaining.dative.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_dative_bonds(
            compaction.dative_bonds(),
            vec![(first, vec![AtomId(0)], AtomId(1), attributes)],
        );
        assert_eq!(
            triatomic
                .molecule
                .dative_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![DativeBondId(0), DativeBondId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_dative_bonds(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.dative = vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = {
            let first = [AtomId(3), AtomId(2)];
            let last = [AtomId(5)];
            editor.add_dative_bonds(vec![
                (&first, AtomId(4), DativeBondForm::from_order(2)),
                (&last, AtomId(0), DativeBondForm::from_order(3)),
            ])
        };

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(DativeBondId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.dative_bond_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.dative.extend([
            (
                vec![AtomId(3), AtomId(2)],
                AtomId(4),
                DativeBondForm::from_order(2),
            ),
            (vec![AtomId(5)], AtomId(0), DativeBondForm::from_order(3)),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![DativeBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_dative_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_dative_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_aromatic_system(mut triatomic: MoleculeEditor) {
        let attributes = AromaticSystemForm::default();
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            aromatic: vec![
                (vec![AtomId(0)], attributes.clone()),
                (vec![AtomId(1), AtomId(2)], attributes.clone()),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_aromatic_system(&[AtomId(0)], attributes.clone());
        let second = triatomic.add_aromatic_system(&[AtomId(1), AtomId(2)], attributes.clone());
        assert_eq!((first, second), (AromaticSystemId(0), AromaticSystemId(1)));
        assert_eq!(
            triatomic
                .molecule
                .aromatic_systems
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![AromaticSystemId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_aromatic_systems(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .aromatic_systems
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![AromaticSystemId(0)]
        );
        let mut remaining = entries;
        remaining.aromatic.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_aromatic_systems(
            compaction.aromatic_systems(),
            vec![(first, vec![AtomId(0)], attributes)],
        );
        assert_eq!(
            triatomic
                .molecule
                .aromatic_systems
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![AromaticSystemId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_aromatic_systems(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.aromatic = vec![(
            vec![AtomId(0), AtomId(1)],
            AromaticSystemForm::from_electrons(vec![1, 1]),
        )];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = {
            let first = [AtomId(3), AtomId(2)];
            let last = [AtomId(5), AtomId(4)];
            editor.add_aromatic_systems(vec![
                (&first, AromaticSystemForm::from_electrons(vec![1, 2])),
                (&last, AromaticSystemForm::from_electrons(vec![2, 1])),
            ])
        };

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(AromaticSystemId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.aromatic_system_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.aromatic.extend([
            (
                vec![AtomId(3), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(5), AtomId(4)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![AromaticSystemId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_aromatic_systems_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_aromatic_systems(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_multicenter_bond(mut triatomic: MoleculeEditor) {
        let attributes = MulticenterBondForm::default();
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            multicenter: vec![
                (vec![AtomId(0), AtomId(1)], attributes.clone()),
                (vec![AtomId(1), AtomId(2)], attributes.clone()),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_multicenter_bond(&[AtomId(0), AtomId(1)], attributes.clone());
        let second = triatomic.add_multicenter_bond(&[AtomId(1), AtomId(2)], attributes.clone());
        assert_eq!(
            (first, second),
            (MulticenterBondId(0), MulticenterBondId(1))
        );
        assert_eq!(
            triatomic
                .molecule
                .multicenter_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![MulticenterBondId(0), MulticenterBondId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_multicenter_bonds(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .multicenter_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![MulticenterBondId(0)]
        );
        let mut remaining = entries;
        remaining.multicenter.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_multicenter_bonds(
            compaction.multicenter_bonds(),
            vec![(first, vec![AtomId(0), AtomId(1)], attributes)],
        );
        assert_eq!(
            triatomic
                .molecule
                .multicenter_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![MulticenterBondId(0), MulticenterBondId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_multicenter_bonds(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.multicenter = vec![(
            vec![AtomId(0), AtomId(1)],
            MulticenterBondForm::from_electrons(vec![1, 1]),
        )];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = {
            let first = [AtomId(3), AtomId(2)];
            let last = [AtomId(5), AtomId(4)];
            editor.add_multicenter_bonds(vec![
                (&first, MulticenterBondForm::from_electrons(vec![1, 2])),
                (&last, MulticenterBondForm::from_electrons(vec![2, 1])),
            ])
        };

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(MulticenterBondId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.multicenter_bond_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.multicenter.extend([
            (
                vec![AtomId(3), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(5), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![MulticenterBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_multicenter_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_multicenter_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_noncovalent_bond(mut triatomic: MoleculeEditor) {
        let attributes = NoncovalentBondForm::default();
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            noncovalent: vec![
                ([AtomId(0), AtomId(1)], attributes.clone()),
                ([AtomId(1), AtomId(2)], attributes.clone()),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_noncovalent_bond([AtomId(0), AtomId(1)], attributes.clone());
        let second = triatomic.add_noncovalent_bond([AtomId(1), AtomId(2)], attributes.clone());
        assert_eq!(
            (first, second),
            (NoncovalentBondId(0), NoncovalentBondId(1))
        );
        assert_eq!(
            triatomic
                .molecule
                .noncovalent_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![NoncovalentBondId(0), NoncovalentBondId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_noncovalent_bonds(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .noncovalent_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![NoncovalentBondId(0)]
        );
        let mut remaining = entries;
        remaining.noncovalent.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_noncovalent_bonds(
            compaction.noncovalent_bonds(),
            vec![(first, [AtomId(0), AtomId(1)], attributes)],
        );
        assert_eq!(
            triatomic
                .molecule
                .noncovalent_bonds
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![NoncovalentBondId(0), NoncovalentBondId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_noncovalent_bonds(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.noncovalent =
            vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::default())];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = editor.add_noncovalent_bonds(vec![
            (
                [AtomId(3), AtomId(2)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            ([AtomId(5), AtomId(4)], NoncovalentBondForm::default()),
        ]);

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(NoncovalentBondId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.noncovalent_bond_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.noncovalent.extend([
            (
                [AtomId(3), AtomId(2)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            ([AtomId(5), AtomId(4)], NoncovalentBondForm::default()),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![NoncovalentBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_noncovalent_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_noncovalent_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_atom(mut triatomic: MoleculeEditor) {
        let attributes = StereoAtomForm::default();
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                    attributes.clone(),
                ),
                (
                    AtomId(1),
                    vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
                    attributes.clone(),
                ),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_stereo_atom(
            AtomId(0),
            &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes.clone(),
        );
        let second = triatomic.add_stereo_atom(
            AtomId(1),
            &[StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
            attributes.clone(),
        );
        assert_eq!((first, second), (StereoAtomId(0), StereoAtomId(1)));
        assert_eq!(
            triatomic
                .molecule
                .stereo_atoms
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoAtomId(0), StereoAtomId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_stereo_atoms(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .stereo_atoms
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoAtomId(0)]
        );
        let mut remaining = entries;
        remaining.stereo_atoms.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_stereo_atoms(
            compaction.stereo_atoms(),
            vec![(
                first,
                AtomId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                attributes,
            )],
        );
        assert_eq!(
            triatomic
                .molecule
                .stereo_atoms
                .incident_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoAtomId(0), StereoAtomId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_atoms(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.stereo_atoms = vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = {
            let first = [
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
            ];
            let last = [
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(4), StereoLigandKind::LonePair),
            ];
            editor.add_stereo_atoms(vec![
                (
                    AtomId(2),
                    &first,
                    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
                ),
                (
                    AtomId(4),
                    &last,
                    StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
                ),
            ])
        };

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(StereoAtomId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.stereo_atom_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.stereo_atoms.extend([
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
            ),
            (
                AtomId(4),
                vec![
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(4), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
            ),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![StereoAtomId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_atoms_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_stereo_atoms(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_bond(mut triatomic: MoleculeEditor) {
        let attributes = StereoBondForm::default();
        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            stereo_bonds: vec![
                (
                    BondId(0),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    attributes.clone(),
                ),
                (
                    BondId(1),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    ],
                    attributes.clone(),
                ),
            ],
            ..Default::default()
        };
        let expected = Molecule::from_entries(entries.clone());
        let first = triatomic.add_stereo_bond(
            BondId(0),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            attributes.clone(),
        );
        let second = triatomic.add_stereo_bond(
            BondId(1),
            &[
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
            ],
            attributes.clone(),
        );
        assert_eq!((first, second), (StereoBondId(0), StereoBondId(1)));
        assert_eq!(
            triatomic
                .molecule
                .stereo_bonds
                .incident_to_atom_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoBondId(0), StereoBondId(1)]
        );
        let snapshot = triatomic.snapshot().unwrap();
        assert_eq!(snapshot, expected);

        let compaction = triatomic.tracked_remove_stereo_bonds(&[first]);
        assert_eq!(
            triatomic
                .molecule
                .stereo_bonds
                .incident_to_atom_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoBondId(0)]
        );
        let mut remaining = entries;
        remaining.stereo_bonds.remove(0);
        assert_eq!(
            triatomic.snapshot().unwrap(),
            Molecule::from_entries(remaining)
        );
        assert_eq!(snapshot, expected);

        triatomic.molecule.restore_stereo_bonds(
            compaction.stereo_bonds(),
            vec![(
                first,
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                attributes,
            )],
        );
        assert_eq!(
            triatomic
                .molecule
                .stereo_bonds
                .incident_to_atom_ids(AtomId(1))
                .collect::<Vec<_>>(),
            vec![StereoBondId(0), StereoBondId(1)]
        );
        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_bonds(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        addition_entries.stereo_bonds = vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 0_u32),
        )];
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        drop(source);
        let mut ids = {
            let first = [
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
            ];
            let last = [
                StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(4), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(5), StereoLigandKind::LonePair),
            ];
            editor.add_stereo_bonds(vec![
                (
                    BondId(2),
                    &first,
                    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                ),
                (
                    BondId(4),
                    &last,
                    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
                ),
            ])
        };

        assert_eq!(ids.len(), 2);
        assert_eq!(ids.next(), Some(StereoBondId(1)));
        assert_eq!(ids.len(), 1);
        assert_eq!(editor.stereo_bond_count(), 3);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.stereo_bonds.extend([
            (
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            ),
            (
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(4), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(5), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 0_u32),
            ),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.try_build(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![StereoBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.snapshot().unwrap();
        let correspondence = triatomic.correspondence.clone();
        let mut ids = triatomic.add_stereo_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(
            triatomic.try_tracked_build(),
            Ok((original, correspondence))
        );
        assert_eq!(ids.next(), None);
    }

    /// Aromatic systems, where the alignment is genuinely used: `on_permutation` reindexes the
    /// electron counts and `is_permutation_invariant` is false for a determinate vector.
    ///
    /// Classes rather than generated inputs — the methods are crate-visible and the property target
    /// cannot reach them.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::different_counts(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 99], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], vec![10, 20], false)]
    fn test_molecule_editor_aromatic_system_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = AromaticSystemForm::from_electrons(electrons);

        assert_eq!(editor.aromatic_system_equiv(AromaticSystemId(0), &atoms, &offered), expected);
    }

    /// Multicenter bonds mirror aromatic systems: one frame-bearing factor with a position-indexed
    /// electron vector, so the alignment is used and the two readings agree.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    fn test_molecule_editor_multicenter_bond_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_multicenter_bond(
            &[AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = MulticenterBondForm::from_electrons(electrons);

        assert_eq!(editor.multicenter_bond_equiv(MulticenterBondId(0), &atoms, &offered), expected);
    }

    /// Noncovalent bonds and dative bonds carry frame-invariant payloads, so the alignment cannot
    /// change the answer and the two readings agree on identity of participants alone.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([AtomId(0), AtomId(1)], true)]
    #[case::reversed_frame([AtomId(1), AtomId(0)], true)]
    #[case::different_pair([AtomId(0), AtomId(2)], false)]
    fn test_molecule_editor_noncovalent_bond_equiv(
        #[case] atoms: [AtomId; 2],
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..3 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let offered = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);

        assert_eq!(editor.noncovalent_bond_equiv(NoncovalentBondId(0), atoms, &offered), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(AtomId(0), vec![AtomId(1), AtomId(2)], true)]
    #[case::reordered_donors(AtomId(0), vec![AtomId(2), AtomId(1)], true)]
    #[case::different_acceptor(AtomId(1), vec![AtomId(1), AtomId(2)], false)]
    #[case::different_donors(AtomId(0), vec![AtomId(1), AtomId(3)], false)]
    fn test_molecule_editor_dative_bond_equiv(
        #[case] acceptor: AtomId,
        #[case] donors: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_dative_bond(
            &[AtomId(1), AtomId(2)],
            AtomId(0),
            DativeBondForm::from_order(1),
        );
        let offered = DativeBondForm::from_order(1);

        assert_eq!(editor.dative_bond_equiv(DativeBondId(0), acceptor, &donors, &offered), expected);
    }

    /// The editor's structural equality must reject a different atom set even when the payload
    /// carries through unread — which is what an undetermined electron vector does.
    ///
    /// Without deriving the participant action, two undetermined forms compare equal even when
    /// the offered and stored atom sets differ.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::reordered_atoms(vec![AtomId(2), AtomId(0), AtomId(1)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], false)]
    fn test_molecule_editor_aromatic_system_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );

        assert_eq!(
            editor.aromatic_system_equiv(
                AromaticSystemId(0),
                &atoms,
                &AromaticSystemForm::default(),
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    fn test_molecule_editor_multicenter_bond_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_multicenter_bond(
            &[AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        );

        assert_eq!(
            editor.multicenter_bond_equiv(
                MulticenterBondId(0),
                &atoms,
                &MulticenterBondForm::default(),
            ),
            expected,
        );
    }

    /// A tetrahedral centre over four distinct ligands, stored in one frame.
    #[fixture]
    fn stereo_editor() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        for element in [Element::F, Element::Cl, Element::Br, Element::I] {
            b.add_atom(AtomForm::from_element(element));
        }
        for ligand in 1..=4 {
            b.add_bond(AtomId(0), AtomId(ligand), BondForm::from_order(1));
        }
        b.add_stereo_atom(
            AtomId(0),
            &(1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect::<Vec<_>>(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        );
        b
    }

    /// A coset is read against its ligand frame, so the offered configuration is restated into the
    /// stored frame before comparison and the same index under a transposed frame denotes the
    /// opposite arrangement.
    ///
    /// Classes rather than generated inputs — the methods are crate-visible and the property target
    /// cannot reach them.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([1, 2, 3, 4], 0, true)]
    #[case::stored_frame_other_coset([1, 2, 3, 4], 1, false)]
    #[case::transposed_frame_same_coset([2, 1, 3, 4], 0, false)]
    #[case::transposed_frame_other_coset([2, 1, 3, 4], 1, true)]
    #[case::multiset_differs([1, 2, 3, 5], 0, false)]
    fn test_molecule_editor_stereo_atom_equiv(
        stereo_editor: MoleculeEditor,
        #[case] ligands: [u32; 4],
        #[case] coset: u32,
        #[case] expected: bool,
    ) {
        let offered: Vec<StereoLigand> = ligands
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset);

        assert_eq!(
            stereo_editor.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &offered, &attributes),
            expected,
        );
    }

    /// The stereo equality check must transport the configuration into the stored ligand frame.
    ///
    /// A coset is read against a frame, so the same index under a swapped frame denotes the
    /// opposite arrangement. Presenting the stored entry's own configuration against a transposed
    /// frame therefore describes a different stereocentre, and the check must say so.
    ///
    #[rstest]
    fn test_molecule_editor_stereo_atom_equiv_reordered_frame(stereo_editor: MoleculeEditor) {
        let stored: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let configuration = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);

        assert!(
            stereo_editor.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &stored, &configuration),
            "the stored frame with its own configuration is equivalent to itself",
        );

        let transposed = Permutation::from_image(&[1, 0, 2, 3]);
        assert!(
            !stereo_editor.stereo_atom_equiv(
                StereoAtomId(0),
                AtomId(0),
                &transposed.act(&stored),
                &configuration,
            ),
            "coset 0 against a transposed frame is the opposite arrangement, not the stored one",
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(BondId(0), vec![2, 3, 4, 5], true)]
    #[case::within_endpoint(BondId(0), vec![3, 2, 4, 5], true)]
    #[case::endpoint_block_swap(BondId(0), vec![4, 5, 2, 3], true)]
    #[case::across_endpoints(BondId(0), vec![2, 4, 3, 5], false)]
    #[case::different_ligand(BondId(0), vec![2, 3, 4, 6], false)]
    #[case::different_site(BondId(1), vec![2, 3, 4, 5], false)]
    fn test_molecule_editor_stereo_bond_equiv(
        #[case] site: BondId,
        #[case] ligand_ids: Vec<u32>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..7 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        for (first, second) in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)] {
            editor.add_bond(AtomId(first), AtomId(second), BondForm::from_order(1));
        }
        editor.add_stereo_bond(
            BondId(0),
            &[2, 3, 4, 5].map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)),
            StereoBondForm::default(),
        );
        let ligands = ligand_ids
            .into_iter()
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();

        assert_eq!(
            editor.stereo_bond_equiv(
                StereoBondId(0),
                site,
                &ligands,
                &StereoBondForm::default(),
            ),
            expected,
        );
    }

    #[rstest]
    fn test_molecule_editor_tracked_remove_roundtrip(mut triatomic: MoleculeEditor) {
        let expected = triatomic.clone().build();
        let removed_atoms = vec![RemovedAtom {
            id: AtomId(1),
            attributes: triatomic.atom(AtomId(1)).attributes().clone(),
        }];
        let removed_bonds = vec![
            RemovedBond {
                id: BondId(0),
                endpoints: triatomic.bond(BondId(0)).atom_ids(),
                attributes: triatomic.bond(BondId(0)).attributes().clone(),
            },
            RemovedBond {
                id: BondId(1),
                endpoints: triatomic.bond(BondId(1)).atom_ids(),
                attributes: triatomic.bond(BondId(1)).attributes().clone(),
            },
        ];

        let compaction = triatomic.tracked_remove(&[AtomId(1)], &[]);
        triatomic
            .molecule
            .restore_topology(compaction.graph(), removed_atoms, removed_bonds);

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_remove_added_topology(mut triatomic: MoleculeEditor) {
        let expected = triatomic.clone().build();
        let added_atom = AddedAtom {
            id: triatomic.add_atom(AtomForm::from_element(Element::F)),
            attributes: AtomForm::from_element(Element::F),
        };
        let added_bond = AddedBond {
            id: triatomic.add_bond(AtomId(2), added_atom.id, BondForm::from_order(1)),
            endpoints: [AtomId(2), added_atom.id],
            attributes: BondForm::from_order(1),
        };

        triatomic.remove_added_topology(&[added_atom], &[added_bond]);

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_remove_dative_bonds_roundtrip() {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_dative_bond(&[AtomId(1)], AtomId(0), DativeBondForm::from_order(2));
        let expected = b.clone().build();
        let view = b.dative_bond(DativeBondId(0));
        let removed = RemovedDativeBond {
            id: DativeBondId(0),
            donors: view.donor_ids().collect(),
            acceptor: view.acceptor_id(),
            attributes: view.attributes().clone(),
        };

        b.remove_dative_bonds(&[DativeBondId(0)]);
        let undo = MoleculeCompaction::new(
            GraphCompaction::new(Compaction::identity(2), Compaction::empty()),
            Compaction::new(2, vec![removed.id])
                .expect("removed entities belong to the source table"),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
        )
        .undo_compaction();
        b.molecule.restore_dative_bonds(
            undo.forward().dative_bonds(),
            vec![(
                removed.id,
                removed.donors,
                removed.acceptor,
                removed.attributes,
            )],
        );

        assert_eq!(b.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_snapshot(mut triatomic: MoleculeEditor) {
        let snapshot = triatomic
            .snapshot()
            .expect("the editor contains an integral molecule");
        triatomic.add_atom(AtomForm::from_element(Element::F));

        assert_eq!(
            snapshot,
            mol_dsl!(r#"{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#)
        );
        assert_eq!(
            triatomic.build(),
            mol_dsl!(r#"{:atoms ["C" "N" "O" "F"] :bonds [[0 1 "1"] [1 2 "2"]]}"#)
        );
    }

    #[rstest]
    #[case::parallel_bond(MoleculeIntegrityError::ParallelBonds {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_snapshot_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] expected: MoleculeIntegrityError,
    ) {
        editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));

        assert_eq!(editor.snapshot(), Err(expected));
    }

    #[rstest]
    #[case::parallel_bond(AtomId(0), AtomId(1), MoleculeIntegrityError::ParallelBonds {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_try_build_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] first: AtomId,
        #[case] second: AtomId,
        #[case] expected: MoleculeIntegrityError,
    ) {
        editor.add_bond(first, second, BondForm::from_order(1));

        assert_eq!(editor.try_build(), Err(expected));
    }

    // `edit()` → `build()` reproduces the molecule including both stereo overlays.
    #[rstest]
    fn test_molecule_editor_build() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "F" "Cl"]
                :bonds [[0 1 "1"] [1 2 "2"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th1"}]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] [:h 2] [:lp 2]] :attrs "Ct1"}]}"#
        );
        assert_eq!(molecule.edit().build(), molecule);
    }

    #[rstest]
    fn test_molecule_editor_try_tracked_build_allocation() {
        let mut carbon = AtomForm::from_element(Element::C);
        carbon.charge = NumForm::Lit(1);

        let entries = MoleculeEntries {
            atoms: vec![
                carbon,
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(
                vec![AtomId(1), AtomId(2)],
                AtomId(3),
                DativeBondForm::from_order(1),
            )],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 0]),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![2, 1, 0]),
            )],
            noncovalent: vec![(
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1u32),
            )],
            constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(2)]),
            })),
        };
        let source = Molecule::from_entries(entries);
        let editor = source.edit();
        let atoms_ptr = editor.correspondence.atoms().matched_pairs().as_ptr();
        let bonds_ptr = editor.correspondence.bonds().matched_pairs().as_ptr();
        let dative_bonds_ptr = editor
            .correspondence
            .dative_bonds()
            .matched_pairs()
            .as_ptr();
        let aromatic_systems_ptr = editor
            .correspondence
            .aromatic_systems()
            .matched_pairs()
            .as_ptr();
        let multicenter_bonds_ptr = editor
            .correspondence
            .multicenter_bonds()
            .matched_pairs()
            .as_ptr();
        let noncovalent_bonds_ptr = editor
            .correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .as_ptr();
        let stereo_atoms_ptr = editor
            .correspondence
            .stereo_atoms()
            .matched_pairs()
            .as_ptr();
        let stereo_bonds_ptr = editor
            .correspondence
            .stereo_bonds()
            .matched_pairs()
            .as_ptr();
        let (result, witness) = editor.try_tracked_build().unwrap();
        assert_eq!(result, source);
        assert_eq!(witness.atoms().matched_pairs().as_ptr(), atoms_ptr);
        assert_eq!(witness.bonds().matched_pairs().as_ptr(), bonds_ptr);
        assert_eq!(
            witness.dative_bonds().matched_pairs().as_ptr(),
            dative_bonds_ptr
        );
        assert_eq!(
            witness.aromatic_systems().matched_pairs().as_ptr(),
            aromatic_systems_ptr
        );
        assert_eq!(
            witness.multicenter_bonds().matched_pairs().as_ptr(),
            multicenter_bonds_ptr
        );
        assert_eq!(
            witness.noncovalent_bonds().matched_pairs().as_ptr(),
            noncovalent_bonds_ptr
        );
        assert_eq!(
            witness.stereo_atoms().matched_pairs().as_ptr(),
            stereo_atoms_ptr
        );
        assert_eq!(
            witness.stereo_bonds().matched_pairs().as_ptr(),
            stereo_bonds_ptr
        );
    }

    // `remove` forward-compacts stereo-atom node refs: removing a non-participant
    // shifts the surviving site/ligand ids; removing the site drops the element.
    #[rstest]
    #[case::remaps_surviving(vec![AtomId(0)], vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]])]
    #[case::drops_on_site_removal(vec![AtomId(1)], vec![])]
    fn test_molecule_editor_remove_stereo_atom(
        #[case] remove_atoms: Vec<AtomId>,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "F" "Cl" "Br"]
                :bonds [[1 2 "1"] [1 3 "1"] [1 4 "1"]]
                :stereo-atoms [{:site 1 :ligands [2 3 4 [:h 1]] :attrs "Th1"}]}"#
        );
        let mut editor = molecule.edit();
        editor.remove(&remove_atoms, &[]);
        let surviving: Vec<Vec<AtomId>> = editor
            .build()
            .stereo_atoms()
            .iter()
            .map(|view| view.atom_ids().collect())
            .collect();
        assert_eq!(surviving, expected);
    }

    // `remove` forward-compacts the stereo-bond edge site: removing a non-site bond
    // shifts the surviving site; removing the site bond drops the element.
    #[rstest]
    #[case::remaps_surviving(vec![BondId(0)], vec![BondId(0)])]
    #[case::drops_on_site_removal(vec![BondId(1)], vec![])]
    fn test_molecule_editor_remove_stereo_bond(
        #[case] remove_bonds: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[4 5 "1"] [1 2 "2"] [0 1 "1"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        );
        let mut editor = molecule.edit();
        editor.remove(&[], &remove_bonds);
        let surviving: Vec<BondId> = editor
            .build()
            .stereo_bonds()
            .iter()
            .map(|view| view.site_id())
            .collect();
        assert_eq!(surviving, expected);
    }
}
