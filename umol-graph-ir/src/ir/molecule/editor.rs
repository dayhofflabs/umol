//! Owning structural and attribute editing for [`Molecule`].
//!
//! Direct mutation and edit batches change transient state without recording undo.
//! [`MoleculeEditor::probe`] lends a molecule after checking integrity;
//! [`MoleculeEditor::finish`] checks integrity and transfers ownership.

use std::sync::Arc;

use umol_graph_core::{Compaction, EdgeId, Graph, GraphCompaction, IdRange};

use super::super::aromatic::{AromaticSystemForm, AromaticSystems};
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::compact::MoleculeCompaction;
use super::super::constraint::Constraints;
use super::super::dative::{DativeBondForm, DativeBonds};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::multicenter::{MulticenterBondForm, MulticenterBonds};
use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBonds};
use super::super::stereo::{StereoAtomForm, StereoAtoms, StereoBondForm, StereoBonds};
use super::super::view::{
    AromaticSystemEditorView, AromaticSystemEditorViewMut, AtomEditorView, AtomEditorViewMut,
    BondEditorView, BondEditorViewMut, DativeBondEditorView, DativeBondEditorViewMut,
    MulticenterBondEditorView, MulticenterBondEditorViewMut, NoncovalentBondEditorView,
    NoncovalentBondEditorViewMut, StereoAtomEditorView, StereoAtomEditorViewMut,
    StereoBondEditorView, StereoBondEditorViewMut,
};
use super::{Molecule, MoleculeIntegrityError, RemovedOverlayIds, TransactionError};

mod transact;

/// Editor for structural and attribute changes to a `Molecule`.
///
/// Publication checks molecule integrity. Removal compacts surviving entity ids.
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
/// Probe and finish establish the same integrity contract as molecule construction. A failed
/// probe retains the editor for further changes. A failed finish drops its state.
#[derive(Clone)]
pub struct MoleculeEditor {
    molecule: Molecule,
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
        constraints: Arc<Constraints>,
    ) -> Self {
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
        }
    }

    /// Append an atom directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. Use Molecule::transact
    /// for checked atomic edits with rollback or apply for consuming application without an undo journal.
    pub fn add_atom(&mut self, atom: AtomForm) -> AtomId {
        self.molecule.add_atom(atom)
    }

    /// Append atoms in input order and return their ids.
    pub fn add_atoms(&mut self, atoms: Vec<AtomForm>) -> IdRange<AtomId> {
        self.molecule.add_atoms(atoms)
    }

    /// Append a localized bond directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. It
    /// assumes `first` and `second` are valid atom ids in the current dense layout.
    pub fn add_bond(&mut self, first: AtomId, second: AtomId, bond: BondForm) -> BondId {
        self.molecule.add_bond(first, second, bond)
    }

    /// Append localized bonds in input order and return their ids.
    ///
    /// Bond endpoints are stored in increasing atom-id order.
    ///
    /// # Panics
    ///
    /// Panics if an endpoint is outside the current atom space.
    pub fn add_bonds(&mut self, bonds: Vec<([AtomId; 2], BondForm)>) -> IdRange<BondId> {
        self.molecule.add_bonds(bonds)
    }

    /// Append atoms and localized bonds in input order and return their ids.
    ///
    /// Bond endpoints refer to the resulting atom space, including the appended atoms, and are
    /// stored in increasing atom-id order.
    ///
    /// # Panics
    ///
    /// Panics if an endpoint is outside the resulting atom space.
    pub fn add_topology(
        &mut self,
        atoms: Vec<AtomForm>,
        bonds: Vec<([AtomId; 2], BondForm)>,
    ) -> (IdRange<AtomId>, IdRange<BondId>) {
        self.molecule.add_topology(atoms, bonds)
    }

    /// Append a dative-bond overlay directly to the editor. The acceptor is factor 1; the donors
    /// are factor 2. Their supplied order is retained as the stored frame.
    pub fn add_dative_bond(
        &mut self,
        donors: &[AtomId],
        acceptor: AtomId,
        bond: DativeBondForm,
    ) -> DativeBondId {
        self.molecule.add_dative_bond(donors, acceptor, bond)
    }

    /// Append dative bonds in input order and return their ids.
    ///
    /// Donor and acceptor ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_dative_bonds(
        &mut self,
        entries: Vec<(&[AtomId], AtomId, DativeBondForm)>,
    ) -> IdRange<DativeBondId> {
        self.molecule.add_dative_bonds(entries)
    }

    /// Append an aromatic-system overlay directly to the editor.
    pub fn add_aromatic_system(
        &mut self,
        atoms: &[AtomId],
        data: AromaticSystemForm,
    ) -> AromaticSystemId {
        self.molecule.add_aromatic_system(atoms, data)
    }

    /// Append aromatic systems in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_aromatic_systems(
        &mut self,
        entries: Vec<(&[AtomId], AromaticSystemForm)>,
    ) -> IdRange<AromaticSystemId> {
        self.molecule.add_aromatic_systems(entries)
    }

    /// Append a multicenter-bond overlay directly to the editor.
    pub fn add_multicenter_bond(
        &mut self,
        atoms: &[AtomId],
        data: MulticenterBondForm,
    ) -> MulticenterBondId {
        self.molecule.add_multicenter_bond(atoms, data)
    }

    /// Append multicenter bonds in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_multicenter_bonds(
        &mut self,
        entries: Vec<(&[AtomId], MulticenterBondForm)>,
    ) -> IdRange<MulticenterBondId> {
        self.molecule.add_multicenter_bonds(entries)
    }

    /// Append a noncovalent-bond overlay directly to the editor.
    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomId; 2],
        bond: NoncovalentBondForm,
    ) -> NoncovalentBondId {
        self.molecule.add_noncovalent_bond(ends, bond)
    }

    /// Append noncovalent bonds in input order and return their ids.
    ///
    /// Atom ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_noncovalent_bonds(
        &mut self,
        entries: Vec<([AtomId; 2], NoncovalentBondForm)>,
    ) -> IdRange<NoncovalentBondId> {
        self.molecule.add_noncovalent_bonds(entries)
    }

    /// Append a stereo-atom overlay directly to the editor.
    pub fn add_stereo_atom(
        &mut self,
        site: AtomId,
        ligands: &[StereoLigand],
        attributes: StereoAtomForm,
    ) -> StereoAtomId {
        self.molecule.add_stereo_atom(site, ligands, attributes)
    }

    /// Append stereo atoms in input order and return their ids.
    ///
    /// Site and ligand ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_stereo_atoms(
        &mut self,
        entries: Vec<(AtomId, &[StereoLigand], StereoAtomForm)>,
    ) -> IdRange<StereoAtomId> {
        self.molecule.add_stereo_atoms(entries)
    }

    /// Append a stereo-bond overlay directly to the editor.
    pub fn add_stereo_bond(
        &mut self,
        site: BondId,
        ligands: &[StereoLigand],
        attributes: StereoBondForm,
    ) -> StereoBondId {
        self.molecule.add_stereo_bond(site, ligands, attributes)
    }

    /// Append stereo bonds in input order and return their ids.
    ///
    /// Site and ligand ids are not checked against the molecule here. Supplied
    /// frame order is retained; molecule integrity is checked at publication.
    pub fn add_stereo_bonds(
        &mut self,
        entries: Vec<(BondId, &[StereoLigand], StereoBondForm)>,
    ) -> IdRange<StereoBondId> {
        self.molecule.add_stereo_bonds(entries)
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
            self.molecule.tracked_remove_dative_bonds(ids),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.compact_constraints(&compaction);
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
            self.molecule.tracked_remove_aromatic_systems(ids),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.compact_constraints(&compaction);
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
            self.molecule.tracked_remove_multicenter_bonds(ids),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.compact_constraints(&compaction);
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
            self.molecule.tracked_remove_noncovalent_bonds(ids),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.compact_constraints(&compaction);
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
            self.molecule.tracked_remove_stereo_atoms(ids),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.molecule.compact_constraints(&compaction);
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
            self.molecule.tracked_remove_stereo_bonds(ids),
        );
        self.molecule.compact_constraints(&compaction);
        compaction
    }

    // -- Topological removal --------------------------------------------------

    /// Remove exactly the supplied atoms and bonds, transporting surviving references.
    ///
    /// Molecule constraints referring to removed entities are compacted as part of removal.
    /// Input order and repeated removal ids do not matter. Failure leaves the editor unchanged.
    ///
    /// # Errors
    ///
    /// Returns TransactionError::InvalidReference for unavailable stored participants or sites,
    /// or TransactionError::DanglingReference when a retained entity references a removed entity.
    /// These checks cover stored references, not reaction applicability or chemical validity.
    ///
    /// # Panics
    ///
    /// Panics if a supplied atom or bond id is outside the current topology.
    pub fn remove_topology(
        &mut self,
        atoms: &[AtomId],
        bonds: &[BondId],
    ) -> Result<(), TransactionError> {
        self.tracked_remove_topology(atoms, bonds).map(|_| ())
    }

    /// Exact removal with the source-to-result compaction for all eight entity kinds.
    ///
    /// Has the same mutation and prerequisites as [`Self::remove_topology`]. Each component retains
    /// its source count, including entity kinds with no removed rows.
    ///
    /// # Errors
    ///
    /// Returns the same reference errors as [`Self::remove_topology`], before mutation.
    ///
    /// # Panics
    ///
    /// Panics if a supplied atom or bond id is outside the current topology.
    pub fn tracked_remove_topology(
        &mut self,
        atoms: &[AtomId],
        bonds: &[BondId],
    ) -> Result<MoleculeCompaction, TransactionError> {
        let mut bonds = bonds.to_vec();
        bonds.sort_unstable();
        self.molecule.check_topology_removal(atoms, &bonds)?;
        let compaction =
            self.molecule
                .remove_selected_topology(atoms, &bonds, RemovedOverlayIds::default());
        self.molecule.compact_constraints(&compaction);
        Ok(compaction)
    }

    /// Remove topology, incident bonds, and overlays with unavailable stored participants or sites.
    ///
    /// Also removes draft overlays whose references were already outside the topology. Surviving
    /// references and molecule constraints are compacted once. Input order and repeated removal ids
    /// do not matter. This does not select additional overlays from reaction applicability rules.
    ///
    /// # Semantic properties
    ///
    /// On a published molecule, explicitly removing the dependent overlay rows and then exactly
    /// removing the complete atom/bond selection gives the same result, including participant frames
    /// and constraints. This equivalence is exercised in `tests/property/molecule/compaction.rs`.
    ///
    /// # Panics
    ///
    /// Panics if a supplied atom or bond id is outside the current topology.
    pub fn remove_topology_cascading(&mut self, atoms: &[AtomId], bonds: &[BondId]) {
        self.tracked_remove_topology_cascading(atoms, bonds);
    }

    /// Cascading removal with the source-to-result compaction for all eight entity kinds.
    ///
    /// Has the same mutation as [`Self::remove_topology_cascading`]. Each component retains its
    /// source count, including entity kinds with no removed rows.
    ///
    /// # Panics
    ///
    /// Panics if a supplied atom or bond id is outside the current topology.
    pub fn tracked_remove_topology_cascading(
        &mut self,
        atoms: &[AtomId],
        bonds: &[BondId],
    ) -> MoleculeCompaction {
        let mut atoms = atoms.to_vec();
        let mut bonds = bonds.to_vec();
        let overlays = self
            .molecule
            .select_topology_cascade(&mut atoms, &mut bonds);
        let compaction = self
            .molecule
            .remove_selected_topology(&atoms, &bonds, overlays);
        self.molecule.compact_constraints(&compaction);
        compaction
    }

    /// Borrow the current molecule after checking its integrity.
    ///
    /// The borrow prevents editing while it remains in use. No molecule copy is made.
    ///
    /// # Errors
    ///
    /// Returns MoleculeIntegrityError if the current state cannot be published as a molecule.
    /// Failure leaves the editor available for inspection and repair.
    pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError> {
        self.molecule.check_integrity()?;
        Ok(&self.molecule)
    }

    /// Consume the editor and publish its molecule after checking integrity.
    ///
    /// Moves the stored topology, attributes, overlays, and constraints into the result.
    ///
    /// # Errors
    ///
    /// Returns MoleculeIntegrityError if the current state cannot be published as a molecule.
    /// Failure drops the editor and its changes.
    pub fn finish(self) -> Result<Molecule, MoleculeIntegrityError> {
        self.molecule.check_integrity()?;
        Ok(self.molecule)
    }
}

#[cfg(test)]
mod tests {
    use std::slice;

    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::{Constraint, MoleculeConstraint};
    use crate::ir::dative::DativeBondForm;
    use crate::ir::edit::{
        AtomHandle, BondHandle, Edit, Edits, RemovedAtom, RemovedBond, RemovedDativeBond,
    };
    use crate::ir::entity::{Entity, EntityKind};
    use crate::ir::error::MoleculeApplyError;
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::molecule::apply::ApplicationState;
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![AtomId(7)]);
    }

    #[rstest]
    fn test_molecule_editor_add_atoms_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_atoms(vec![]);

        assert!(Arc::ptr_eq(&triatomic.molecule.atoms, &original.atoms));

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![BondId(7)]);
    }

    #[rstest]
    fn test_molecule_editor_add_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_bonds(vec![]);

        assert!(Arc::ptr_eq(&triatomic.molecule.bonds, &original.bonds));

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_molecule_editor_add_topology(
        mut addition_entries: MoleculeEntries,
        #[values(false, true)] shared: bool,
    ) {
        let source = Molecule::from_entries(addition_entries.clone());
        let original = shared.then(|| source.clone());
        let mut editor = source.edit();
        let (mut atom_ids, mut bond_ids) = editor.add_topology(
            vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            vec![
                ([AtomId(7), AtomId(1)], BondForm::from_order(1)),
                ([AtomId(6), AtomId(7)], BondForm::from_order(2)),
            ],
        );

        assert_eq!(atom_ids.len(), 2);
        assert_eq!(atom_ids.next(), Some(AtomId(6)));
        assert_eq!(bond_ids.len(), 2);
        assert_eq!(bond_ids.next(), Some(BondId(6)));
        assert_eq!(editor.atom_count(), 8);
        assert_eq!(editor.bond_count(), 8);
        if let Some(original) = original {
            assert_eq!(original, Molecule::from_entries(addition_entries.clone()));
        }
        addition_entries.atoms.extend([
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
        ]);
        addition_entries.bonds.extend([
            (AtomId(1), AtomId(7), BondForm::from_order(1)),
            (AtomId(6), AtomId(7), BondForm::from_order(2)),
        ]);
        let expected = Molecule::from_entries(addition_entries);
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(atom_ids.collect::<Vec<_>>(), vec![AtomId(7)]);
        assert_eq!(bond_ids.collect::<Vec<_>>(), vec![BondId(7)]);
    }

    #[rstest]
    fn test_molecule_editor_add_topology_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let (mut atom_ids, mut bond_ids) = triatomic.add_topology(vec![], vec![]);

        assert!(Arc::ptr_eq(&triatomic.molecule.atoms, &original.atoms));
        assert!(Arc::ptr_eq(&triatomic.molecule.bonds, &original.bonds));

        assert_eq!(atom_ids.len(), 0);
        assert_eq!(bond_ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
        assert_eq!(atom_ids.next(), None);
        assert_eq!(bond_ids.next(), None);
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![DativeBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_dative_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_dative_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![AromaticSystemId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_aromatic_systems_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_aromatic_systems(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![MulticenterBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_multicenter_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_multicenter_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![NoncovalentBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_noncovalent_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_noncovalent_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![StereoAtomId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_atoms_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_stereo_atoms(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
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
        let snapshot = triatomic.probe().unwrap().clone();
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
            triatomic.probe().unwrap().clone(),
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
        assert_eq!(triatomic.finish().unwrap(), expected);
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
        assert_eq!(editor.finish(), Ok(expected));
        assert_eq!(ids.collect::<Vec<_>>(), vec![StereoBondId(2)]);
    }

    #[rstest]
    fn test_molecule_editor_add_stereo_bonds_identity(mut triatomic: MoleculeEditor) {
        let original = triatomic.probe().unwrap().clone();
        let mut ids = triatomic.add_stereo_bonds(vec![]);

        assert_eq!(ids.len(), 0);
        assert_eq!(triatomic.finish(), Ok(original));
        assert_eq!(ids.next(), None);
    }

    #[derive(Clone, Copy)]
    enum RemovalMethod {
        Direct,
        Tracked,
        Apply,
        Journal,
    }

    #[rstest]
    #[case::bond(Entity::Bond(BondId(0)), false)]
    #[case::dative(Entity::DativeBond(DativeBondId(0)), false)]
    #[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)), false)]
    #[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)), false)]
    #[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)), false)]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)), false)]
    #[case::stereo_bond_ligand(Entity::StereoBond(StereoBondId(0)), false)]
    #[case::stereo_bond_site(Entity::StereoBond(StereoBondId(0)), true)]
    fn test_molecule_editor_remove_topology_error(
        #[case] dependent: Entity,
        #[case] site: bool,
        #[values(
            RemovalMethod::Direct,
            RemovalMethod::Tracked,
            RemovalMethod::Apply,
            RemovalMethod::Journal
        )]
        method: RemovalMethod,
    ) {
        let mut editor = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
            ..Default::default()
        })
        .edit();
        match dependent.kind() {
            EntityKind::Bond => {}
            EntityKind::DativeBond => {
                editor.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::default());
            }
            EntityKind::AromaticSystem => {
                editor.add_aromatic_system(&[AtomId(0), AtomId(1)], AromaticSystemForm::default());
            }
            EntityKind::MulticenterBond => {
                editor
                    .add_multicenter_bond(&[AtomId(0), AtomId(1)], MulticenterBondForm::default());
            }
            EntityKind::NoncovalentBond => {
                editor.add_noncovalent_bond([AtomId(0), AtomId(1)], NoncovalentBondForm::default());
            }
            EntityKind::StereoAtom => {
                editor.add_stereo_atom(AtomId(0), &[], StereoAtomForm::default());
            }
            EntityKind::StereoBond => {
                editor.add_stereo_bond(
                    BondId(0),
                    &[StereoLigand::new(
                        AtomId(0),
                        StereoLigandKind::ImplicitHydrogen,
                    )],
                    StereoBondForm::default(),
                );
            }
            EntityKind::Atom => unreachable!(),
        }
        let before = editor.molecule.clone();
        let atoms = if site { vec![] } else { vec![AtomId(0)] };
        let bonds = if dependent.kind() == EntityKind::Bond {
            vec![]
        } else {
            vec![BondId(0)]
        };
        let edit = Edit::RemoveTopology {
            atoms: atoms.iter().copied().map(AtomHandle::Id).collect(),
            bonds: bonds.iter().copied().map(BondHandle::Id).collect(),
        };
        let result = match method {
            RemovalMethod::Direct => editor
                .remove_topology(&atoms, &bonds)
                .map_err(MoleculeApplyError::from),
            RemovalMethod::Tracked => editor
                .tracked_remove_topology(&atoms, &bonds)
                .map(|_| ())
                .map_err(MoleculeApplyError::from),
            RemovalMethod::Apply => editor.apply(Edits::from_iter([edit])),
            RemovalMethod::Journal => {
                let mut state = ApplicationState::new(&editor.molecule);
                editor
                    .molecule
                    .apply_edit_with_undo(edit, &mut state)
                    .map(|_| ())
                    .map_err(MoleculeApplyError::from)
            }
        };
        assert_eq!(
            result,
            Err(TransactionError::DanglingReference {
                removed: if site {
                    Entity::Bond(BondId(0))
                } else {
                    Entity::Atom(AtomId(0))
                },
                dependent,
            }
            .into())
        );
        assert_eq!(editor.molecule, before);
    }

    #[rstest]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_editor_remove_topology_adjacency(
        #[case] kind: EntityKind,
        #[values(false, true)] cascading: bool,
    ) {
        let mut entries = MoleculeEntries {
            atoms: vec![AtomForm::default(); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::default()),
                (AtomId(0), AtomId(2), BondForm::default()),
            ],
            ..Default::default()
        };
        match kind {
            EntityKind::AromaticSystem => entries
                .aromatic
                .push((vec![AtomId(0), AtomId(2)], AromaticSystemForm::default())),
            EntityKind::StereoAtom => entries.stereo_atoms.push((
                AtomId(0),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
                StereoAtomForm::default(),
            )),
            EntityKind::StereoBond => entries.stereo_bonds.push((
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            )),
            _ => unreachable!(),
        }
        let source = Molecule::from_entries(entries);
        let mut expected = source.clone();
        expected.graph = Graph::new(3, &[[0, 1]]);
        expected.bonds = Arc::new(vec![BondForm::default()]);
        let mut editor = source.clone().edit();
        if cascading {
            editor.remove_topology_cascading(&[], &[BondId(1)]);
        } else {
            editor.remove_topology(&[], &[BondId(1)]).unwrap();
        }
        assert_eq!(editor.molecule, expected);
        editor.add_bond(AtomId(0), AtomId(2), BondForm::default());
        assert_eq!(editor.finish(), Ok(source));
    }

    #[rstest]
    #[case::dative(EntityKind::DativeBond, false)]
    #[case::aromatic(EntityKind::AromaticSystem, false)]
    #[case::multicenter(EntityKind::MulticenterBond, false)]
    #[case::noncovalent(EntityKind::NoncovalentBond, false)]
    #[case::stereo_atom_site(EntityKind::StereoAtom, true)]
    #[case::stereo_atom_ligand(EntityKind::StereoAtom, false)]
    #[case::stereo_bond_site(EntityKind::StereoBond, true)]
    #[case::stereo_bond_ligand(EntityKind::StereoBond, false)]
    fn test_molecule_editor_remove_topology_draft(
        #[case] kind: EntityKind,
        #[case] site: bool,
        #[values(
            RemovalMethod::Direct,
            RemovalMethod::Tracked,
            RemovalMethod::Apply,
            RemovalMethod::Journal
        )]
        method: RemovalMethod,
    ) {
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
            ..Default::default()
        });
        let mut editor = source.clone().edit();
        let invalid = AtomId(2);
        let ligand = StereoLigand::new(invalid, StereoLigandKind::Atom);
        match kind {
            EntityKind::DativeBond => {
                editor.add_dative_bond(&[invalid], AtomId(0), DativeBondForm::default());
            }
            EntityKind::AromaticSystem => {
                editor.add_aromatic_system(&[AtomId(0), invalid], AromaticSystemForm::default());
            }
            EntityKind::MulticenterBond => {
                editor.add_multicenter_bond(&[AtomId(0), invalid], MulticenterBondForm::default());
            }
            EntityKind::NoncovalentBond => {
                editor.add_noncovalent_bond([AtomId(0), invalid], NoncovalentBondForm::default());
            }
            EntityKind::StereoAtom => {
                editor.add_stereo_atom(
                    if site { invalid } else { AtomId(0) },
                    if site { &[] } else { slice::from_ref(&ligand) },
                    StereoAtomForm::default(),
                );
            }
            EntityKind::StereoBond => {
                editor.add_stereo_bond(
                    if site { BondId(1) } else { BondId(0) },
                    if site { &[] } else { slice::from_ref(&ligand) },
                    StereoBondForm::default(),
                );
            }
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        }
        let before = editor.molecule.clone();
        let edit = Edit::RemoveTopology {
            atoms: vec![],
            bonds: vec![],
        };
        let result = match method {
            RemovalMethod::Direct => editor
                .remove_topology(&[], &[])
                .map_err(MoleculeApplyError::from),
            RemovalMethod::Tracked => editor
                .tracked_remove_topology(&[], &[])
                .map(|_| ())
                .map_err(MoleculeApplyError::from),
            RemovalMethod::Apply => editor.apply(Edits::from_iter([edit])),
            RemovalMethod::Journal => {
                let mut state = ApplicationState::new(&editor.molecule);
                editor
                    .molecule
                    .apply_edit_with_undo(edit, &mut state)
                    .map(|_| ())
                    .map_err(MoleculeApplyError::from)
            }
        };
        assert_eq!(
            result,
            Err(TransactionError::InvalidReference {
                entity: if kind == EntityKind::StereoBond && site {
                    Entity::Bond(BondId(1))
                } else {
                    Entity::Atom(invalid)
                },
            }
            .into())
        );
        assert_eq!(editor.molecule, before);
        let edit = Edit::RemoveTopologyCascading {
            atoms: vec![],
            bonds: vec![],
        };
        match method {
            RemovalMethod::Direct => editor.remove_topology_cascading(&[], &[]),
            RemovalMethod::Tracked => {
                editor.tracked_remove_topology_cascading(&[], &[]);
            }
            RemovalMethod::Apply => editor.apply(Edits::from_iter([edit])).unwrap(),
            RemovalMethod::Journal => {
                let mut state = ApplicationState::new(&editor.molecule);
                let undo = editor
                    .molecule
                    .apply_edit_with_undo(edit, &mut state)
                    .unwrap()
                    .unwrap();
                assert_eq!(editor.molecule, source);
                editor.molecule.apply_undo(undo);
                assert_eq!(editor.molecule, before);
                return;
            }
        }
        assert_eq!(editor.finish(), Ok(source));
    }

    #[rstest]
    #[case::ordered(vec![AtomId(1)], vec![BondId(0), BondId(1)])]
    #[case::repeated(vec![AtomId(1), AtomId(1)], vec![BondId(1), BondId(0), BondId(1)])]
    fn test_molecule_editor_tracked_remove_topology_roundtrip(
        mut triatomic: MoleculeEditor,
        #[case] atoms: Vec<AtomId>,
        #[case] bonds: Vec<BondId>,
        #[values(false, true)] cascading: bool,
    ) {
        let expected = triatomic.clone().finish().unwrap();
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

        let compaction = if cascading {
            triatomic.tracked_remove_topology_cascading(&atoms, &bonds)
        } else {
            triatomic.tracked_remove_topology(&atoms, &bonds).unwrap()
        };
        assert_eq!(
            compaction.graph(),
            &GraphCompaction::new(
                Compaction::new(3, vec![umol_graph_core::NodeId(1)]).unwrap(),
                Compaction::new(2, vec![EdgeId(0), EdgeId(1)]).unwrap(),
            )
        );
        triatomic
            .molecule
            .restore_topology(compaction.graph(), removed_atoms, removed_bonds);

        assert_eq!(triatomic.finish().unwrap(), expected);
    }

    #[rstest]
    fn test_molecule_editor_remove_dative_bonds_roundtrip() {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_dative_bond(&[AtomId(1)], AtomId(0), DativeBondForm::from_order(2));
        let expected = b.clone().finish().unwrap();
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

        assert_eq!(b.finish().unwrap(), expected);
    }

    #[rstest]
    fn test_molecule_editor_probe(mut triatomic: MoleculeEditor) {
        let molecule = triatomic.probe().unwrap();
        assert_eq!(
            molecule,
            &mol_dsl!(r#"{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#)
        );
        triatomic.add_atom(AtomForm::from_element(Element::F));
        assert_eq!(
            triatomic.probe(),
            Ok(&mol_dsl!(
                r#"{:atoms ["C" "N" "O" "F"] :bonds [[0 1 "1"] [1 2 "2"]]}"#
            ))
        );
    }

    #[rstest]
    #[case::parallel_bond(MoleculeIntegrityError::ParallelBonds {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_probe_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let duplicate = editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));

        assert_eq!(editor.probe(), Err(expected));

        editor.remove_topology(&[], &[duplicate]).unwrap();
        assert_eq!(
            editor.finish(),
            Ok(mol_dsl!(
                r#"{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#
            ))
        );
    }

    #[rstest]
    fn test_molecule_editor_finish() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "F" "Cl"]
                :bonds [[0 1 "1"] [1 2 "2"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th1"}]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] [:h 2] [:lp 2]] :attrs "Ct1"}]}"#
        );
        let expected = molecule.clone();

        assert_eq!(molecule.edit().finish(), Ok(expected));
    }

    #[rstest]
    #[case::parallel_bond(AtomId(0), AtomId(1), MoleculeIntegrityError::ParallelBonds {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_finish_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] first: AtomId,
        #[case] second: AtomId,
        #[case] expected: MoleculeIntegrityError,
    ) {
        editor.add_bond(first, second, BondForm::from_order(1));

        let entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (first, second, BondForm::from_order(1)),
            ],
            ..Default::default()
        };
        assert_eq!(Molecule::try_from_entries(entries), Err(expected.clone()));
        assert_eq!(editor.finish(), Err(expected));
    }

    #[rstest]
    fn test_molecule_editor_finish_allocation() {
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
        let mut expected = Molecule::from_entries(entries.clone());
        expected.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(2);
        let source = Molecule::from_entries(entries);
        let atoms_ptr = source.atoms.as_ptr();
        let bonds_ptr = source.bonds.as_ptr();
        let constraints_ptr = source.constraints.iter().as_slice().as_ptr();

        let mut editor = source.edit();
        editor.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(2);
        let result = editor.finish().unwrap();

        assert_eq!(result, expected);
        assert_eq!(result.atoms.as_ptr(), atoms_ptr);
        assert_eq!(result.bonds.as_ptr(), bonds_ptr);
        assert_eq!(
            result.constraints.iter().as_slice().as_ptr(),
            constraints_ptr
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
        editor.remove_topology_cascading(&remove_atoms, &[]);
        let surviving: Vec<Vec<AtomId>> = editor
            .finish()
            .unwrap()
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
        editor.remove_topology_cascading(&[], &remove_bonds);
        let surviving: Vec<BondId> = editor
            .finish()
            .unwrap()
            .stereo_bonds()
            .iter()
            .map(|view| view.site_id())
            .collect();
        assert_eq!(surviving, expected);
    }
}
