//! Molecule graph IR.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::Arc;
use std::{iter, mem};

pub use build::MoleculeBuilder;
pub use editor::MoleculeEditor;
pub use fragment::{Fragment, Port, PortArg};
pub use integrity::MoleculeIntegrityError;
pub use spec::{AtomArg, MoleculeSpec, MoleculeSpecTerm};
use umol_graph_core::{
    Correspondence, EdgeId, Graph, NodeId, RelationParticipant, Remapping, UnionFind,
};

use super::aromatic::{reframe_aromatic_systems_with, AromaticSystemForm, AromaticSystems};
use super::atom::AtomForm;
use super::bond::BondForm;
use super::constraint::{
    Constraint, ConstraintFrameActionMap, Constraints, MoleculeConstraint, RelationalConstraint,
};
use super::correspondence::MoleculeCorrespondence;
use super::dative::{reframe_dative_bonds_with, DativeBondForm, DativeBonds};
use super::edit::{AtomHandle, BondHandle, Edits};
use super::entity::{Entity, EntityKind};
use super::error::{Contradiction, MoleculeApplyError};
use super::frame::OverlaysFrameAction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::multicenter::{reframe_multicenter_bonds_with, MulticenterBondForm, MulticenterBonds};
use super::noncovalent::{reframe_noncovalent_bonds_with, NoncovalentBondForm, NoncovalentBonds};
use super::remap::IdRemapping;
use super::ring::{RingConfig, RingModel, RingSet};
use super::stereo::{
    reframe_stereo_atoms_with, reframe_stereo_bonds_with, StereoAtomForm, StereoAtoms,
    StereoBondForm, StereoBonds,
};
use super::traits::{FrameTransport, Lattice, Normalize, Reframe};
use super::view::{
    AromaticSystemView, AromaticSystemViewMut, AromaticSystemViews, AtomView, AtomViewMut,
    AtomViews, BondView, BondViewMut, BondViews, DativeBondView, DativeBondViewMut,
    DativeBondViews, GraphView, MulticenterBondView, MulticenterBondViewMut, MulticenterBondViews,
    NeighborView, NoncovalentBondView, NoncovalentBondViewMut, NoncovalentBondViews, RingViews,
    StereoAtomView, StereoAtomViews, StereoBondView, StereoBondViews,
};

mod build;
mod editor;
mod fragment;
pub(crate) mod integrity;
mod pushout;
mod remapping;
pub mod spec;
pub(crate) mod transact;

/// Molecule graph IR: atom-bond topology, overlays (typed hyperedges), and constraints.
///
/// Per-entity data are `Arc`-shared (copy-on-write). The molecule itself only allows
/// attribute mutation; structural edits go through `MoleculeEditor` via [`Molecule::edit`].
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct Molecule {
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
}

/// Constructor input for [`Molecule::from_entries`].
#[derive(Debug, Default, Clone)]
pub struct MoleculeEntries {
    pub atoms: Vec<AtomForm>,
    pub bonds: Vec<(AtomId, AtomId, BondForm)>,
    pub dative: Vec<(Vec<AtomId>, AtomId, DativeBondForm)>,
    pub aromatic: Vec<(Vec<AtomId>, AromaticSystemForm)>,
    pub multicenter: Vec<(Vec<AtomId>, MulticenterBondForm)>,
    pub noncovalent: Vec<(AtomId, AtomId, NoncovalentBondForm)>,
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>,
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>,
    pub constraints: Constraints,
}

pub(super) fn validate_entry_references(
    entries: &MoleculeEntries,
) -> Result<(), MoleculeIntegrityError> {
    validate_entry_references_inner(entries)
        .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })
}

fn validate_entry_references_inner(entries: &MoleculeEntries) -> Result<(), Entity> {
    let contains = |entity| match entity {
        Entity::Atom(id) => id.index() < entries.atoms.len(),
        Entity::Bond(id) => id.index() < entries.bonds.len(),
        Entity::DativeBond(id) => id.index() < entries.dative.len(),
        Entity::AromaticSystem(id) => id.index() < entries.aromatic.len(),
        Entity::MulticenterBond(id) => id.index() < entries.multicenter.len(),
        Entity::NoncovalentBond(id) => id.index() < entries.noncovalent.len(),
        Entity::StereoAtom(id) => id.index() < entries.stereo_atoms.len(),
        Entity::StereoBond(id) => id.index() < entries.stereo_bonds.len(),
    };

    for &(first, second, _) in &entries.bonds {
        require_reference(&contains, Entity::Atom(first))?;
        require_reference(&contains, Entity::Atom(second))?;
    }
    for (donors, acceptor, _) in &entries.dative {
        require_references(&contains, donors.iter().copied().map(Entity::Atom))?;
        require_reference(&contains, Entity::Atom(*acceptor))?;
    }
    for (atoms, _) in &entries.aromatic {
        require_references(&contains, atoms.iter().copied().map(Entity::Atom))?;
    }
    for (atoms, _) in &entries.multicenter {
        require_references(&contains, atoms.iter().copied().map(Entity::Atom))?;
    }
    for &(first, second, _) in &entries.noncovalent {
        require_reference(&contains, Entity::Atom(first))?;
        require_reference(&contains, Entity::Atom(second))?;
    }
    for (site, ligands, _) in &entries.stereo_atoms {
        require_reference(&contains, Entity::Atom(*site))?;
        require_references(
            &contains,
            ligands.iter().map(|ligand| Entity::Atom(ligand.atom_id)),
        )?;
    }
    for (site, ligands, _) in &entries.stereo_bonds {
        require_reference(&contains, Entity::Bond(*site))?;
        require_references(
            &contains,
            ligands.iter().map(|ligand| Entity::Atom(ligand.atom_id)),
        )?;
    }
    for constraint in entries.constraints.iter() {
        validate_constraint_references(constraint, &contains)?;
    }
    Ok(())
}

fn require_reference(contains: &dyn Fn(Entity) -> bool, entity: Entity) -> Result<(), Entity> {
    if contains(entity) {
        Ok(())
    } else {
        Err(entity)
    }
}

fn require_references(
    contains: &dyn Fn(Entity) -> bool,
    entities: impl IntoIterator<Item = Entity>,
) -> Result<(), Entity> {
    for entity in entities {
        require_reference(contains, entity)?;
    }
    Ok(())
}

pub(super) fn validate_constraint_references(
    constraint: &Constraint,
    contains: &dyn Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        Constraint::Atom(id, _) => require_reference(contains, Entity::Atom(*id)),
        Constraint::Bond(id, _) => require_reference(contains, Entity::Bond(*id)),
        Constraint::DativeBond(id, _) => require_reference(contains, Entity::DativeBond(*id)),
        Constraint::AromaticSystem(id, _) => {
            require_reference(contains, Entity::AromaticSystem(*id))
        }
        Constraint::MulticenterBond(id, _) => {
            require_reference(contains, Entity::MulticenterBond(*id))
        }
        Constraint::NoncovalentBond(id, _) => {
            require_reference(contains, Entity::NoncovalentBond(*id))
        }
        Constraint::StereoAtom(id, _, _) => require_reference(contains, Entity::StereoAtom(*id)),
        Constraint::StereoBond(id, _, _) => require_reference(contains, Entity::StereoBond(*id)),
        Constraint::Relational(constraint) => {
            validate_relational_constraint_references(constraint, contains)
        }
        Constraint::Molecule(constraint) => {
            validate_molecule_constraint_references(constraint, contains)
        }
        Constraint::And(constraints) | Constraint::Or(constraints) => {
            for constraint in constraints {
                validate_constraint_references(constraint, contains)?;
            }
            Ok(())
        }
        Constraint::Not(constraint) => validate_constraint_references(constraint, contains),
    }
}

fn validate_relational_constraint_references(
    constraint: &RelationalConstraint,
    contains: &dyn Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        RelationalConstraint::DativeBondDonors { bond, atoms }
        | RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
            require_reference(contains, Entity::DativeBond(*bond))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::DativeBondDonor { bond, atom }
        | RelationalConstraint::DativeBondAcceptor { bond, atom } => {
            require_reference(contains, Entity::DativeBond(*bond))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::DativeBondAllDonors { bond, .. }
        | RelationalConstraint::DativeBondAnyDonor { bond, .. }
        | RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => {
            require_reference(contains, Entity::DativeBond(*bond))
        }
        RelationalConstraint::DativeBondParallels { dative, parallel } => {
            require_reference(contains, Entity::DativeBond(*dative))?;
            require_reference(contains, Entity::Bond(*parallel))
        }
        RelationalConstraint::AromaticSystemAtoms { system, atoms }
        | RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
            require_reference(contains, Entity::AromaticSystem(*system))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::AromaticSystemContains { system, atom } => {
            require_reference(contains, Entity::AromaticSystem(*system))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::AromaticSystemAllAtoms { system, .. }
        | RelationalConstraint::AromaticSystemAnyAtom { system, .. } => {
            require_reference(contains, Entity::AromaticSystem(*system))
        }
        RelationalConstraint::MulticenterBondAtoms { bond, atoms }
        | RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
            require_reference(contains, Entity::MulticenterBond(*bond))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::MulticenterBondContains { bond, atom } => {
            require_reference(contains, Entity::MulticenterBond(*bond))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::MulticenterBondAllAtoms { bond, .. }
        | RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => {
            require_reference(contains, Entity::MulticenterBond(*bond))
        }
        RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
            require_reference(contains, Entity::NoncovalentBond(*bond))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::NoncovalentBondContains { bond, atom } => {
            require_reference(contains, Entity::NoncovalentBond(*bond))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
            require_reference(contains, Entity::NoncovalentBond(*bond))
        }
        RelationalConstraint::StereoAtomSite { stereo_atom, atom }
        | RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
            require_reference(contains, Entity::StereoAtom(*stereo_atom))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
            require_reference(contains, Entity::StereoAtom(*stereo_atom))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::StereoAtomAllLigands { stereo_atom, .. }
        | RelationalConstraint::StereoAtomAnyLigand { stereo_atom, .. } => {
            require_reference(contains, Entity::StereoAtom(*stereo_atom))
        }
        RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
            require_reference(contains, Entity::StereoBond(*stereo_bond))?;
            require_reference(contains, Entity::Bond(*bond))
        }
        RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
            require_reference(contains, Entity::StereoBond(*stereo_bond))?;
            require_reference(contains, Entity::Atom(*atom))
        }
        RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
            require_reference(contains, Entity::StereoBond(*stereo_bond))?;
            require_references(contains, atoms.iter().copied().map(Entity::Atom))
        }
        RelationalConstraint::StereoBondAllLigands { stereo_bond, .. }
        | RelationalConstraint::StereoBondAnyLigand { stereo_bond, .. } => {
            require_reference(contains, Entity::StereoBond(*stereo_bond))
        }
    }
}

fn validate_molecule_constraint_references(
    constraint: &MoleculeConstraint,
    contains: &dyn Fn(Entity) -> bool,
) -> Result<(), Entity> {
    match constraint {
        MoleculeConstraint::ChargeSum { atoms, .. }
        | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
        | MoleculeConstraint::Connected { atoms } => {
            require_references(contains, atoms.iter().flatten().copied().map(Entity::Atom))
        }
        MoleculeConstraint::BondOrderSum { bonds, .. } => {
            require_references(contains, bonds.iter().flatten().copied().map(Entity::Bond))
        }
    }
}

impl Molecule {
    /// Empty molecule: zero atoms, zero bonds, zero overlays, zero constraints.
    pub fn new() -> Self {
        Self::default()
    }

    /// Start an empty `MoleculeEditor` for fluent / programmatic
    /// construction. Use [`Molecule::edit`] to start from an existing
    /// molecule.
    pub fn builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    /// Full structural constructor from a flat [`MoleculeEntries`]: every entity-type field is
    /// supplied directly. The topology-only case fills just `atoms` and `bonds`; relations and
    /// molecule-level constraints go in the remaining fields.
    ///
    /// # Panics
    ///
    /// Panics if an entry references an unavailable entity or otherwise violates molecule
    /// representation integrity. Use [`Self::try_from_entries`] for untrusted input.
    pub fn from_entries(entries: MoleculeEntries) -> Self {
        Self::try_from_entries(entries)
            .unwrap_or_else(|error| panic!("invalid molecule entries: {error}"))
    }

    /// Checked form of [`Self::from_entries`]. Validates molecule representation integrity,
    /// including the fixed simple-relation semantics of every entity kind, but does not enforce
    /// chemistry or constraint satisfiability.
    pub fn try_from_entries(entries: MoleculeEntries) -> Result<Self, MoleculeIntegrityError> {
        validate_entry_references(&entries)?;
        let MoleculeEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        } = entries;
        let node_count = atoms.len();
        let edges: Vec<[u32; 2]> = bonds
            .iter()
            .map(|(first, second, _)| [first.0, second.0])
            .collect();
        let bond_data: Vec<BondForm> = bonds.into_iter().map(|(_, _, d)| d).collect();
        let graph = Graph::new(node_count, &edges);

        Self::try_from_arcs(
            graph,
            Arc::new(atoms),
            Arc::new(bond_data),
            DativeBonds::new(dative),
            AromaticSystems::new(aromatic),
            MulticenterBonds::new(multicenter),
            NoncovalentBonds::new(noncovalent),
            StereoAtoms::new(stereo_atoms),
            StereoBonds::new(stereo_bonds),
            constraints,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn try_from_arcs(
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
    ) -> Result<Self, MoleculeIntegrityError> {
        let molecule = Self {
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
        };
        molecule.check_integrity()?;
        Ok(molecule)
    }

    /// AtomId/BondId-typed adapter exposing the pure-graph algorithms.
    pub fn graph(&self) -> GraphView<'_> {
        GraphView::new(&self.graph)
    }

    /// Raw underlying graph with `NodeId` / `EdgeId` types. Escape hatch
    /// for code that needs the graph-core API directly; use [`Self::graph`]
    /// for AtomId/BondId-typed access.
    pub fn raw_graph(&self) -> &Graph {
        &self.graph
    }

    /// Complete framed equality under a dense entity-id correspondence from `self` to `other`.
    ///
    /// The correspondence is checked as a complete dense remapping of `self`. The remapped
    /// molecule is then compared with `other` modulo participant frames. This returns `false`
    /// when the correspondence has the wrong source domain, is partial, is not bijective onto a
    /// dense target domain, or maps `self` to a molecule outside `other`'s framed-equality class.
    ///
    /// Under the identity correspondence this is exactly [`Reframe::framed_eq`]. Reversing
    /// a correspondence reverses the comparison, and sequential correspondences compose.
    /// For integrity-valid molecules whose complete canonicalizations both succeed,
    /// `canonical_eq` holds exactly when some admissible total correspondence makes this
    /// comparison true. Equality totalization for two intrinsic contradictions does not require
    /// such a witness.
    pub fn framed_eq_under(&self, other: &Self, correspondence: &MoleculeCorrespondence) -> bool {
        self.try_remap(correspondence)
            .is_some_and(|remapped| remapped.framed_eq(other))
    }

    /// Neighbors of `atom`, ordered by ascending neighbor atom id.
    pub fn neighbors(&self, atom: AtomId) -> impl ExactSizeIterator<Item = NeighborView<'_>> {
        self.graph
            .neighbors(NodeId::from(atom))
            .iter()
            .map(move |n| NeighborView::new(AtomId::from(n.node), BondId::from(n.edge), self))
    }

    // Slices have to store node/edge ids, otherwise need owned vectors.
    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(self, &self.atoms)
    }

    /// View of the atom with `id`.
    ///
    /// Panics if `id` is not an atom in this molecule. Use
    /// `self.atoms().get(id)` for checked lookup.
    pub fn atom(&self, id: AtomId) -> AtomView<'_> {
        self.atoms()
            .get(id)
            .expect("atom id must refer to an atom in this molecule")
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(self, &self.bonds)
    }

    /// View of the bond with `id`.
    ///
    /// Panics if `id` is not a bond in this molecule. Use
    /// `self.bonds().get(id)` for checked lookup.
    pub fn bond(&self, id: BondId) -> BondView<'_> {
        self.bonds()
            .get(id)
            .expect("bond id must refer to a bond in this molecule")
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(self, &self.dative_bonds)
    }

    /// View of the dative bond with `id`.
    ///
    /// Panics if `id` is not a dative bond in this molecule. Use
    /// `self.dative_bonds().get(id)` for checked lookup.
    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondView<'_> {
        self.dative_bonds()
            .get(id)
            .expect("dative bond id must refer to a dative bond in this molecule")
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(self, &self.aromatic_systems)
    }

    /// View of the aromatic system with `id`.
    ///
    /// Panics if `id` is not an aromatic system in this molecule. Use
    /// `self.aromatic_systems().get(id)` for checked lookup.
    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemView<'_> {
        self.aromatic_systems()
            .get(id)
            .expect("aromatic system id must refer to an aromatic system in this molecule")
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(self, &self.multicenter_bonds)
    }

    /// View of the multicenter bond with `id`.
    ///
    /// Panics if `id` is not a multicenter bond in this molecule. Use
    /// `self.multicenter_bonds().get(id)` for checked lookup.
    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondView<'_> {
        self.multicenter_bonds()
            .get(id)
            .expect("multicenter bond id must refer to a multicenter bond in this molecule")
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(self, &self.noncovalent_bonds)
    }

    /// View of the noncovalent bond with `id`.
    ///
    /// Panics if `id` is not a noncovalent bond in this molecule. Use
    /// `self.noncovalent_bonds().get(id)` for checked lookup.
    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondView<'_> {
        self.noncovalent_bonds()
            .get(id)
            .expect("noncovalent bond id must refer to a noncovalent bond in this molecule")
    }

    pub fn stereo_atoms(&self) -> StereoAtomViews<'_> {
        StereoAtomViews::new(self, &self.stereo_atoms)
    }

    /// View of the stereo atom with `id`.
    ///
    /// Panics if `id` is not a stereo atom in this molecule. Use
    /// `self.stereo_atoms().get(id)` for checked lookup.
    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomView<'_> {
        self.stereo_atoms()
            .get(id)
            .expect("stereo atom id must refer to a stereo atom in this molecule")
    }

    pub fn stereo_bonds(&self) -> StereoBondViews<'_> {
        StereoBondViews::new(self, &self.stereo_bonds)
    }

    /// View of the stereo bond with `id`.
    ///
    /// Panics if `id` is not a stereo bond in this molecule. Use
    /// `self.stereo_bonds().get(id)` for checked lookup.
    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondView<'_> {
        self.stereo_bonds()
            .get(id)
            .expect("stereo bond id must refer to a stereo bond in this molecule")
    }

    /// The subgraph induced by `atoms` (deduplicated, host order preserved), as an injective
    /// sub→host [`MoleculeCorrespondence`]: sub entity `i` maps to its host id. `extract` / `edits`
    /// materialize it. A bond/overlay is included iff all its participants are in `atoms`.
    pub fn induced_subgraph(&self, atoms: &[AtomId]) -> MoleculeCorrespondence {
        let mut host_atoms: Vec<AtomId> = Vec::with_capacity(atoms.len());
        let mut atom_set: HashSet<AtomId> = HashSet::with_capacity(atoms.len());
        for &a in atoms {
            if atom_set.insert(a) {
                host_atoms.push(a);
            }
        }

        let host_bonds: Vec<BondId> = self
            .bonds()
            .iter()
            .filter(|b| {
                let [a, b_end] = b.atom_ids();
                atom_set.contains(&a) && atom_set.contains(&b_end)
            })
            .map(|b| b.id)
            .collect();
        let host_dative_bonds: Vec<DativeBondId> = self
            .dative_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_aromatic_systems: Vec<AromaticSystemId> = self
            .aromatic_systems()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_multicenter_bonds: Vec<MulticenterBondId> = self
            .multicenter_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_noncovalent_bonds: Vec<NoncovalentBondId> = self
            .noncovalent_bonds()
            .iter()
            .filter(|v| {
                let [a, b] = v.atom_ids();
                atom_set.contains(&a) && atom_set.contains(&b)
            })
            .map(|v| v.id)
            .collect();
        let host_stereo_atoms: Vec<StereoAtomId> = self
            .stereo_atoms()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_stereo_bonds: Vec<StereoBondId> = self
            .stereo_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();

        MoleculeCorrespondence::new(
            Correspondence::from_images(&host_atoms, self.atoms().count()),
            Correspondence::from_images(&host_bonds, self.bonds().count()),
            Correspondence::from_images(&host_dative_bonds, self.dative_bonds().count()),
            Correspondence::from_images(&host_aromatic_systems, self.aromatic_systems().count()),
            Correspondence::from_images(&host_multicenter_bonds, self.multicenter_bonds().count()),
            Correspondence::from_images(&host_noncovalent_bonds, self.noncovalent_bonds().count()),
            Correspondence::from_images(&host_stereo_atoms, self.stereo_atoms().count()),
            Correspondence::from_images(&host_stereo_bonds, self.stereo_bonds().count()),
        )
    }

    /// Materialize an induced-subgraph correspondence `sub` (over `self` as host) as a standalone
    /// molecule: drop every host atom/bond absent from `sub`. Host order preserved, gaps compacted;
    /// overlay drops cascade through the builder.
    pub fn extract(&self, sub: &MoleculeCorrespondence) -> Molecule {
        let kept: HashSet<AtomId> = sub
            .atoms()
            .matched_pairs()
            .iter()
            .map(|&(_, host)| host)
            .collect();
        let remove_atoms: Vec<AtomId> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|a| !kept.contains(a))
            .collect();
        let remove_bonds: Vec<BondId> = self
            .bonds()
            .iter()
            .filter(|b| {
                let [a, b_end] = b.atom_ids();
                !kept.contains(&a) || !kept.contains(&b_end)
            })
            .map(|b| b.id)
            .collect();
        let mut builder = self.edit();
        builder.remove(&remove_atoms, &remove_bonds);
        builder.build()
    }

    /// Edits transforming `self` into the extracted subgraph `sub`: one `RemoveTopology` over the
    /// host atoms/bonds not in `sub` (empty when `sub` covers the whole molecule).
    pub fn edits(&self, sub: &MoleculeCorrespondence) -> Edits {
        let kept: HashSet<AtomId> = sub
            .atoms()
            .matched_pairs()
            .iter()
            .map(|&(_, host)| host)
            .collect();
        let kept_bonds: HashSet<BondId> = sub
            .bonds()
            .matched_pairs()
            .iter()
            .map(|&(_, host)| host)
            .collect();
        let removed_atoms: Vec<AtomHandle> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|a| !kept.contains(a))
            .map(AtomHandle::Id)
            .collect();
        let removed_bonds: Vec<BondHandle> = (0..self.bonds().count())
            .map(BondId::from)
            .filter(|b| !kept_bonds.contains(b))
            .map(BondHandle::Id)
            .collect();
        if removed_atoms.is_empty() && removed_bonds.is_empty() {
            return Edits::new();
        }
        let mut edits = Edits::new();
        edits.remove_topology(removed_atoms, removed_bonds);
        edits
    }

    /// Concreteness: all entities have ground inherent fields.
    pub fn is_concrete(&self) -> bool {
        self.atoms.iter().all(|atom| atom.is_concrete())
            && self.bonds.iter().all(|bond| bond.is_concrete())
            && self
                .dative_bonds
                .ids()
                .all(|id| self.dative_bonds.attributes(id).is_concrete())
            && self
                .aromatic_systems
                .ids()
                .all(|id| self.aromatic_systems.attributes(id).is_concrete())
            && self
                .multicenter_bonds
                .ids()
                .all(|id| self.multicenter_bonds.attributes(id).is_concrete())
            && self
                .noncovalent_bonds
                .ids()
                .all(|id| self.noncovalent_bonds.attributes(id).is_concrete())
            && self
                .stereo_atoms
                .ids()
                .all(|id| self.stereo_atoms.attributes(id).is_concrete())
            && self
                .stereo_bonds
                .ids()
                .all(|id| self.stereo_bonds.attributes(id).is_concrete())
    }

    /// Rings selected by `model` and computed using `config`.
    pub fn rings(&self, model: RingModel, config: RingConfig) -> RingViews<'_> {
        RingViews::new(self, RingSet::enumerate(&self.graph, model, config))
    }

    pub fn atom_mut(&mut self, id: AtomId) -> AtomViewMut<'_> {
        let attributes = &mut Arc::make_mut(&mut self.atoms)[id.index()];
        AtomViewMut { id, attributes }
    }

    /// Replace every atom with `f(atom)` in place (owned in, owned out — no
    /// `&mut AtomForm` escapes, so the container controls any re-interning).
    pub fn modify_atoms(&mut self, mut f: impl FnMut(AtomForm) -> AtomForm) {
        for atom in Arc::make_mut(&mut self.atoms).iter_mut() {
            *atom = f(mem::take(atom));
        }
    }

    pub fn bond_mut(&mut self, id: BondId) -> BondViewMut<'_> {
        let [s, t] = self.graph.edge_endpoints(id.into());
        let atoms = [AtomId::from(s), AtomId::from(t)];
        let attributes = &mut Arc::make_mut(&mut self.bonds)[id.index()];
        BondViewMut {
            id,
            atoms,
            attributes,
        }
    }

    /// Replace every bond with `f(bond)` in place.
    pub fn modify_bonds(&mut self, mut f: impl FnMut(BondForm) -> BondForm) {
        for bond in Arc::make_mut(&mut self.bonds).iter_mut() {
            *bond = f(mem::take(bond));
        }
    }

    pub fn dative_bond_mut(&mut self, id: DativeBondId) -> DativeBondViewMut<'_> {
        let acceptor = self.dative_bonds.acceptor(id);
        let donors = self.dative_bonds.donors(id).collect();
        let attributes = self.dative_bonds.attributes_mut(id);
        DativeBondViewMut {
            id,
            donors,
            acceptor,
            attributes,
        }
    }

    /// Replace every dative bond with `f(bond)` in place.
    pub fn modify_dative_bonds(&mut self, mut f: impl FnMut(DativeBondForm) -> DativeBondForm) {
        for dative_bond in self.dative_bonds.attributes_iter_mut() {
            *dative_bond = f(mem::take(dative_bond));
        }
    }

    pub(crate) fn aromatic_system_mut(
        &mut self,
        id: AromaticSystemId,
    ) -> AromaticSystemViewMut<'_> {
        let atoms = self.aromatic_systems.atoms(id).collect();
        let attributes = self.aromatic_systems.attributes_mut(id);
        AromaticSystemViewMut {
            id,
            atoms,
            attributes,
        }
    }

    /// Replace every aromatic system with `f(system)` in place.
    pub(crate) fn modify_aromatic_systems(
        &mut self,
        mut f: impl FnMut(AromaticSystemForm) -> AromaticSystemForm,
    ) {
        for aromatic_system in self.aromatic_systems.attributes_iter_mut() {
            *aromatic_system = f(mem::take(aromatic_system));
        }
    }

    /// Transactionally modify one aromatic-system form.
    ///
    /// The callback operates on a private candidate. The candidate replaces this molecule only if
    /// it still satisfies molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError::InvalidReference`] if `id` is unavailable, or the exact
    /// integrity error introduced by the callback. On error, this molecule is unchanged.
    pub fn try_modify_aromatic_system(
        &mut self,
        id: AromaticSystemId,
        f: impl FnOnce(&mut AromaticSystemForm),
    ) -> Result<(), MoleculeIntegrityError> {
        if !self.aromatic_systems.contains(id) {
            return Err(MoleculeIntegrityError::InvalidReference {
                entity: Entity::AromaticSystem(id),
            });
        }
        self.try_modify_checked(|candidate| f(candidate.aromatic_systems.attributes_mut(id)))
    }

    /// Transactionally modify every aromatic-system form.
    ///
    /// The callback operates on forms in a private candidate. The candidate replaces this molecule
    /// only if all modified forms still satisfy molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns the first molecule integrity error introduced by the callback. On error, this
    /// molecule is unchanged.
    pub fn try_modify_aromatic_systems(
        &mut self,
        mut f: impl FnMut(&mut AromaticSystemForm),
    ) -> Result<(), MoleculeIntegrityError> {
        self.try_modify_checked(|candidate| {
            for aromatic_system in candidate.aromatic_systems.attributes_iter_mut() {
                f(aromatic_system);
            }
        })
    }

    pub(crate) fn multicenter_bond_mut(
        &mut self,
        id: MulticenterBondId,
    ) -> MulticenterBondViewMut<'_> {
        let atoms = self.multicenter_bonds.atoms(id).collect();
        let attributes = self.multicenter_bonds.attributes_mut(id);
        MulticenterBondViewMut {
            id,
            atoms,
            attributes,
        }
    }

    /// Replace every multicenter bond with `f(bond)` in place.
    pub(crate) fn modify_multicenter_bonds(
        &mut self,
        mut f: impl FnMut(MulticenterBondForm) -> MulticenterBondForm,
    ) {
        for multicenter_bond in self.multicenter_bonds.attributes_iter_mut() {
            *multicenter_bond = f(mem::take(multicenter_bond));
        }
    }

    /// Transactionally modify one multicenter-bond form.
    ///
    /// The callback operates on a private candidate. The candidate replaces this molecule only if
    /// it still satisfies molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError::InvalidReference`] if `id` is unavailable, or the exact
    /// integrity error introduced by the callback. On error, this molecule is unchanged.
    pub fn try_modify_multicenter_bond(
        &mut self,
        id: MulticenterBondId,
        f: impl FnOnce(&mut MulticenterBondForm),
    ) -> Result<(), MoleculeIntegrityError> {
        if !self.multicenter_bonds.contains(id) {
            return Err(MoleculeIntegrityError::InvalidReference {
                entity: Entity::MulticenterBond(id),
            });
        }
        self.try_modify_checked(|candidate| f(candidate.multicenter_bonds.attributes_mut(id)))
    }

    /// Transactionally modify every multicenter-bond form.
    ///
    /// The callback operates on forms in a private candidate. The candidate replaces this molecule
    /// only if all modified forms still satisfy molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns the first molecule integrity error introduced by the callback. On error, this
    /// molecule is unchanged.
    pub fn try_modify_multicenter_bonds(
        &mut self,
        mut f: impl FnMut(&mut MulticenterBondForm),
    ) -> Result<(), MoleculeIntegrityError> {
        self.try_modify_checked(|candidate| {
            for multicenter_bond in candidate.multicenter_bonds.attributes_iter_mut() {
                f(multicenter_bond);
            }
        })
    }

    pub fn noncovalent_bond_mut(&mut self, id: NoncovalentBondId) -> NoncovalentBondViewMut<'_> {
        let atoms = self.noncovalent_bonds.atoms(id);
        let attributes = self.noncovalent_bonds.attributes_mut(id);
        NoncovalentBondViewMut {
            id,
            atoms,
            attributes,
        }
    }

    /// Replace every noncovalent bond with `f(bond)` in place.
    pub fn modify_noncovalent_bonds(
        &mut self,
        mut f: impl FnMut(NoncovalentBondForm) -> NoncovalentBondForm,
    ) {
        for noncovalent_bond in self.noncovalent_bonds.attributes_iter_mut() {
            *noncovalent_bond = f(mem::take(noncovalent_bond));
        }
    }

    pub(crate) fn stereo_atom_mut(&mut self, id: StereoAtomId) -> &mut StereoAtomForm {
        self.stereo_atoms.attributes_mut(id)
    }

    /// Replace every stereo atom with `f(stereo_atom)` in place.
    pub(crate) fn modify_stereo_atoms(
        &mut self,
        mut f: impl FnMut(StereoAtomForm) -> StereoAtomForm,
    ) {
        for stereo_atom in self.stereo_atoms.attributes_iter_mut() {
            *stereo_atom = f(mem::take(stereo_atom));
        }
    }

    /// Transactionally modify one stereo-atom form.
    ///
    /// The callback operates on a private candidate. The candidate replaces this molecule only if
    /// it still satisfies molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError::InvalidReference`] if `id` is unavailable, or the exact
    /// integrity error introduced by the callback. On error, this molecule is unchanged.
    pub fn try_modify_stereo_atom(
        &mut self,
        id: StereoAtomId,
        f: impl FnOnce(&mut StereoAtomForm),
    ) -> Result<(), MoleculeIntegrityError> {
        if !self.stereo_atoms.contains(id) {
            return Err(MoleculeIntegrityError::InvalidReference {
                entity: Entity::StereoAtom(id),
            });
        }
        self.try_modify_checked(|candidate| f(candidate.stereo_atoms.attributes_mut(id)))
    }

    /// Transactionally modify every stereo-atom form.
    ///
    /// The callback operates on forms in a private candidate. The candidate replaces this molecule
    /// only if all modified forms still satisfy molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns the first molecule integrity error introduced by the callback. On error, this
    /// molecule is unchanged.
    pub fn try_modify_stereo_atoms(
        &mut self,
        mut f: impl FnMut(&mut StereoAtomForm),
    ) -> Result<(), MoleculeIntegrityError> {
        self.try_modify_checked(|candidate| {
            for stereo_atom in candidate.stereo_atoms.attributes_iter_mut() {
                f(stereo_atom);
            }
        })
    }

    pub(crate) fn stereo_bond_mut(&mut self, id: StereoBondId) -> &mut StereoBondForm {
        self.stereo_bonds.attributes_mut(id)
    }

    /// Replace every stereo bond with `f(stereo_bond)` in place.
    pub(crate) fn modify_stereo_bonds(
        &mut self,
        mut f: impl FnMut(StereoBondForm) -> StereoBondForm,
    ) {
        for stereo_bond in self.stereo_bonds.attributes_iter_mut() {
            *stereo_bond = f(mem::take(stereo_bond));
        }
    }

    /// Transactionally modify one stereo-bond form.
    ///
    /// The callback operates on a private candidate. The candidate replaces this molecule only if
    /// it still satisfies molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError::InvalidReference`] if `id` is unavailable, or the exact
    /// integrity error introduced by the callback. On error, this molecule is unchanged.
    pub fn try_modify_stereo_bond(
        &mut self,
        id: StereoBondId,
        f: impl FnOnce(&mut StereoBondForm),
    ) -> Result<(), MoleculeIntegrityError> {
        if !self.stereo_bonds.contains(id) {
            return Err(MoleculeIntegrityError::InvalidReference {
                entity: Entity::StereoBond(id),
            });
        }
        self.try_modify_checked(|candidate| f(candidate.stereo_bonds.attributes_mut(id)))
    }

    /// Transactionally modify every stereo-bond form.
    ///
    /// The callback operates on forms in a private candidate. The candidate replaces this molecule
    /// only if all modified forms still satisfy molecule representation integrity.
    ///
    /// # Errors
    ///
    /// Returns the first molecule integrity error introduced by the callback. On error, this
    /// molecule is unchanged.
    pub fn try_modify_stereo_bonds(
        &mut self,
        mut f: impl FnMut(&mut StereoBondForm),
    ) -> Result<(), MoleculeIntegrityError> {
        self.try_modify_checked(|candidate| {
            for stereo_bond in candidate.stereo_bonds.attributes_iter_mut() {
                f(stereo_bond);
            }
        })
    }

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    #[cfg(test)]
    pub(crate) fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    /// Transactionally modify the molecule-level constraint tree.
    ///
    /// The callback operates on a private candidate. The candidate replaces this molecule only if
    /// all constraint references and stereo wrapper domains remain valid.
    ///
    /// # Errors
    ///
    /// Returns the first molecule integrity error introduced by the callback. On error, this
    /// molecule is unchanged.
    pub fn try_modify_constraints(
        &mut self,
        f: impl FnOnce(&mut Constraints),
    ) -> Result<(), MoleculeIntegrityError> {
        self.try_modify_checked(|candidate| f(&mut candidate.constraints))
    }

    fn try_modify_checked(
        &mut self,
        f: impl FnOnce(&mut Self),
    ) -> Result<(), MoleculeIntegrityError> {
        let mut candidate = self.clone();
        f(&mut candidate);
        candidate.check_integrity()?;
        *self = candidate;
        Ok(())
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// True if the molecule-scope `Constraints` tree is non-empty.
    /// Per-entity constraint stores are not consulted.
    pub fn has_constraints(&self) -> bool {
        !self.constraints.is_empty()
    }

    pub fn has_dative_bonds(&self) -> bool {
        self.dative_bonds.count() > 0
    }

    pub fn has_aromatic_systems(&self) -> bool {
        self.aromatic_systems.count() > 0
    }

    pub fn has_multicenter_bonds(&self) -> bool {
        self.multicenter_bonds.count() > 0
    }

    pub fn has_noncovalent_bonds(&self) -> bool {
        self.noncovalent_bonds.count() > 0
    }

    pub fn has_stereo_atoms(&self) -> bool {
        self.stereo_atoms.count() > 0
    }

    pub fn has_stereo_bonds(&self) -> bool {
        self.stereo_bonds.count() > 0
    }

    /// Drain every entity's inline `constraints` store into `self.constraints`
    /// as `Constraint::Atom` / `Bond` / `DativeBond` / `AromaticSystem` /
    /// `MulticenterBond` / `NoncovalentBond` / `StereoAtom` / `StereoBond`
    /// entries. The order of inserted
    /// entries in `self.constraints` is unspecified.
    pub fn lift_constraints(&mut self) {
        let atom_count = self.atoms().count();
        let bond_count = self.bonds().count();
        let dative_count = self.dative_bonds().count();
        let aromatic_count = self.aromatic_systems().count();
        let multicenter_count = self.multicenter_bonds().count();
        let noncovalent_count = self.noncovalent_bonds().count();
        let stereo_atom_count = self.stereo_atoms().count();
        let stereo_bond_count = self.stereo_bonds().count();

        let mut additions: Vec<Constraint> = Vec::new();
        for i in 0..atom_count {
            let id = AtomId::from(i);
            for c in self.atom_mut(id).attributes.constraints.take() {
                additions.push(Constraint::Atom(id, c));
            }
        }
        for i in 0..bond_count {
            let id = BondId::from(i);
            for c in self.bond_mut(id).attributes.constraints.take() {
                additions.push(Constraint::Bond(id, c));
            }
        }
        for i in 0..dative_count {
            let id = DativeBondId::from(i);
            for c in self.dative_bond_mut(id).attributes.constraints.take() {
                additions.push(Constraint::DativeBond(id, c));
            }
        }
        for i in 0..aromatic_count {
            let id = AromaticSystemId::from(i);
            for c in self.aromatic_system_mut(id).attributes.constraints.take() {
                additions.push(Constraint::AromaticSystem(id, c));
            }
        }
        for i in 0..multicenter_count {
            let id = MulticenterBondId::from(i);
            for c in self.multicenter_bond_mut(id).attributes.constraints.take() {
                additions.push(Constraint::MulticenterBond(id, c));
            }
        }
        for i in 0..noncovalent_count {
            let id = NoncovalentBondId::from(i);
            for c in self.noncovalent_bond_mut(id).attributes.constraints.take() {
                additions.push(Constraint::NoncovalentBond(id, c));
            }
        }
        for i in 0..stereo_atom_count {
            let id = StereoAtomId::from(i);
            let kind = self
                .stereo_atom_mut(id)
                .configuration
                .kind()
                .expect("molecule stereo atom has a concrete kind");
            for c in self.stereo_atom_mut(id).constraints.take() {
                additions.push(Constraint::StereoAtom(id, kind, c));
            }
        }
        for i in 0..stereo_bond_count {
            let id = StereoBondId::from(i);
            let kind = self
                .stereo_bond_mut(id)
                .configuration
                .kind()
                .expect("molecule stereo bond has a concrete kind");
            for c in self.stereo_bond_mut(id).constraints.take() {
                additions.push(Constraint::StereoBond(id, kind, c));
            }
        }
        for c in additions {
            self.constraints.push(c);
        }
    }

    /// Push every entity inline constraints from `self.constraints`
    /// into the targeted entity's inline `constraints` store, removing it
    /// from the molecule list. A leaf colliding with a stored entry of the
    /// same key combines by meet; a meet to `⊥` is a contradiction and the
    /// molecule is left unchanged. Combinator subtrees, `Relational`, and
    /// `Molecule` entries are left in place.
    pub fn inline_constraints(&mut self) -> Result<(), Contradiction> {
        // Validate every collision before mutating anything.
        let mut planned: Vec<Constraint> = Vec::with_capacity(self.constraints.iter().count());
        for c in self.constraints.iter() {
            let met = match c {
                Constraint::Atom(id, inner) => {
                    let met = match self.atom(*id).attributes.constraints.get(inner.key()) {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::Atom(*id, met)
                }
                Constraint::Bond(id, inner) => {
                    let met = match self.bond(*id).attributes.constraints.get(inner.key()) {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::Bond(*id, met)
                }
                Constraint::DativeBond(id, inner) => {
                    let met = match self
                        .dative_bond(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::DativeBond(*id, met)
                }
                Constraint::AromaticSystem(id, inner) => {
                    let met = match self
                        .aromatic_system(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::AromaticSystem(*id, met)
                }
                Constraint::MulticenterBond(id, inner) => {
                    let met = match self
                        .multicenter_bond(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::MulticenterBond(*id, met)
                }
                Constraint::NoncovalentBond(id, inner) => {
                    let met = match self
                        .noncovalent_bond(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::NoncovalentBond(*id, met)
                }
                Constraint::StereoAtom(id, kind, inner) => {
                    let met = match self
                        .stereo_atom(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::StereoAtom(*id, *kind, met)
                }
                Constraint::StereoBond(id, kind, inner) => {
                    let met = match self
                        .stereo_bond(*id)
                        .attributes
                        .constraints
                        .get(inner.key())
                    {
                        Some(existing) => existing.meet(inner).ok_or(Contradiction)?,
                        None => inner.clone(),
                    };
                    Constraint::StereoBond(*id, *kind, met)
                }
                c @ (Constraint::Relational(_)
                | Constraint::Molecule(_)
                | Constraint::And(_)
                | Constraint::Or(_)
                | Constraint::Not(_)) => c.clone(),
            };
            planned.push(met);
        }

        self.constraints.take();
        for c in planned {
            match c {
                Constraint::Atom(id, inner) => {
                    self.atom_mut(id).attributes.constraints.set(inner);
                }
                Constraint::Bond(id, inner) => {
                    self.bond_mut(id).attributes.constraints.set(inner);
                }
                Constraint::DativeBond(id, inner) => {
                    self.dative_bond_mut(id).attributes.constraints.set(inner);
                }
                Constraint::AromaticSystem(id, inner) => {
                    self.aromatic_system_mut(id)
                        .attributes
                        .constraints
                        .set(inner);
                }
                Constraint::MulticenterBond(id, inner) => {
                    self.multicenter_bond_mut(id)
                        .attributes
                        .constraints
                        .set(inner);
                }
                Constraint::NoncovalentBond(id, inner) => {
                    self.noncovalent_bond_mut(id)
                        .attributes
                        .constraints
                        .set(inner);
                }
                // The carried kind is dropped here; kind/degree consistency
                // against the element is the C4 validator's job.
                Constraint::StereoAtom(id, _kind, inner) => {
                    self.stereo_atom_mut(id).constraints.set(inner);
                }
                Constraint::StereoBond(id, _kind, inner) => {
                    self.stereo_bond_mut(id).constraints.set(inner);
                }
                c @ (Constraint::Relational(_)
                | Constraint::Molecule(_)
                | Constraint::And(_)
                | Constraint::Or(_)
                | Constraint::Not(_)) => self.constraints.push(c),
            }
        }
        Ok(())
    }

    pub fn edit(&self) -> MoleculeEditor {
        MoleculeEditor::from_parts(
            self.graph.clone(),
            Arc::clone(&self.atoms),
            Arc::clone(&self.bonds),
            self.dative_bonds.clone(),
            self.aromatic_systems.clone(),
            self.multicenter_bonds.clone(),
            self.noncovalent_bonds.clone(),
            self.stereo_atoms.clone(),
            self.stereo_bonds.clone(),
            self.constraints.clone(),
        )
    }

    /// Apply a checked edit batch to an immutable molecule, returning the modified molecule while
    /// leaving `self` unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeApplyError::Transaction`] when an edit handle, precondition, or shape is
    /// invalid for the evolving draft. Returns [`MoleculeApplyError::Integrity`] when the modified
    /// draft cannot be published as a representation-integral molecule.
    ///
    /// # Semantic properties
    ///
    /// On success, the returned molecule passes [`Molecule::check_integrity`]. On either failure,
    /// `self` remains unchanged and no partially modified molecule is returned.
    pub fn apply(&self, edits: Edits) -> Result<Molecule, MoleculeApplyError> {
        let editor = self.edit().apply(edits)?;
        Ok(editor.try_build()?)
    }

    /// Combine molecules by disjoint concatenation. Input order determines each entity kind's
    /// id ranges in the result. Returns one correspondence per input, in input order, mapping that
    /// molecule's ids into the combined molecule. Pure renumbering — no gluing, no chemistry.
    pub fn combine_all<'a>(
        molecules: impl IntoIterator<Item = &'a Molecule>,
    ) -> (Molecule, Vec<MoleculeCorrespondence>) {
        let molecules: Vec<&Molecule> = molecules.into_iter().collect();
        let atom_count = molecules.iter().map(|m| m.atoms().count()).sum();
        let bond_count = molecules.iter().map(|m| m.bonds().count()).sum();
        let dative_count = molecules.iter().map(|m| m.dative_bonds().count()).sum();
        let aromatic_count = molecules.iter().map(|m| m.aromatic_systems().count()).sum();
        let multicenter_count = molecules
            .iter()
            .map(|m| m.multicenter_bonds().count())
            .sum();
        let noncovalent_count = molecules
            .iter()
            .map(|m| m.noncovalent_bonds().count())
            .sum();
        let stereo_atom_count = molecules.iter().map(|m| m.stereo_atoms().count()).sum();
        let stereo_bond_count = molecules.iter().map(|m| m.stereo_bonds().count()).sum();

        let mut entries = MoleculeEntries {
            atoms: Vec::with_capacity(atom_count),
            bonds: Vec::with_capacity(bond_count),
            dative: Vec::with_capacity(dative_count),
            aromatic: Vec::with_capacity(aromatic_count),
            multicenter: Vec::with_capacity(multicenter_count),
            noncovalent: Vec::with_capacity(noncovalent_count),
            stereo_atoms: Vec::with_capacity(stereo_atom_count),
            stereo_bonds: Vec::with_capacity(stereo_bond_count),
            constraints: Constraints::new(),
        };
        let mut correspondences = Vec::with_capacity(molecules.len());
        let mut atom_offset = 0;
        let mut bond_offset = 0;
        let mut dative_offset = 0;
        let mut aromatic_offset = 0;
        let mut multicenter_offset = 0;
        let mut noncovalent_offset = 0;
        let mut stereo_atom_offset = 0;
        let mut stereo_bond_offset = 0;

        for molecule in molecules {
            let molecule_atom_count = molecule.atoms().count();
            let molecule_bond_count = molecule.bonds().count();
            let molecule_dative_count = molecule.dative_bonds().count();
            let molecule_aromatic_count = molecule.aromatic_systems().count();
            let molecule_multicenter_count = molecule.multicenter_bonds().count();
            let molecule_noncovalent_count = molecule.noncovalent_bonds().count();
            let molecule_stereo_atom_count = molecule.stereo_atoms().count();
            let molecule_stereo_bond_count = molecule.stereo_bonds().count();
            let shift_atom = |id: AtomId| AtomId(id.0 + atom_offset as u32);

            entries
                .atoms
                .extend(molecule.atoms().iter().map(|atom| atom.attributes.clone()));
            entries.bonds.extend(molecule.bonds().iter().map(|bond| {
                let [first, second] = bond.atom_ids();
                (
                    shift_atom(first),
                    shift_atom(second),
                    bond.attributes.clone(),
                )
            }));
            entries
                .dative
                .extend(molecule.dative_bonds().iter().map(|bond| {
                    (
                        bond.donors().map(|donor| shift_atom(donor.id)).collect(),
                        shift_atom(bond.acceptor_id()),
                        bond.attributes.clone(),
                    )
                }));
            entries
                .aromatic
                .extend(molecule.aromatic_systems().iter().map(|system| {
                    (
                        system.atom_ids().map(shift_atom).collect(),
                        system.attributes.clone(),
                    )
                }));
            entries
                .multicenter
                .extend(molecule.multicenter_bonds().iter().map(|bond| {
                    (
                        bond.atom_ids().map(shift_atom).collect(),
                        bond.attributes.clone(),
                    )
                }));
            entries
                .noncovalent
                .extend(molecule.noncovalent_bonds().iter().map(|bond| {
                    let [first, second] = bond.atom_ids();
                    (
                        shift_atom(first),
                        shift_atom(second),
                        bond.attributes.clone(),
                    )
                }));
            for id in molecule.stereo_atoms.ids() {
                let site = shift_atom(molecule.stereo_atoms.site(id));
                let ligands = molecule
                    .stereo_atoms
                    .ligands(id)
                    .iter()
                    .map(|ligand| StereoLigand::new(shift_atom(ligand.atom_id), ligand.kind))
                    .collect();
                entries.stereo_atoms.push((
                    site,
                    ligands,
                    molecule.stereo_atoms.attributes(id).clone(),
                ));
            }
            for id in molecule.stereo_bonds.ids() {
                let site = BondId(molecule.stereo_bonds.site(id).0 + bond_offset as u32);
                let ligands = molecule
                    .stereo_bonds
                    .ligands(id)
                    .iter()
                    .map(|ligand| StereoLigand::new(shift_atom(ligand.atom_id), ligand.kind))
                    .collect();
                entries.stereo_bonds.push((
                    site,
                    ligands,
                    molecule.stereo_bonds.attributes(id).clone(),
                ));
            }

            if !molecule.constraints.is_empty() {
                let remapping = IdRemapping::new(
                    offset_map(atom_offset, molecule_atom_count),
                    offset_map(bond_offset, molecule_bond_count),
                    offset_map(dative_offset, molecule_dative_count),
                    offset_map(aromatic_offset, molecule_aromatic_count),
                    offset_map(multicenter_offset, molecule_multicenter_count),
                    offset_map(noncovalent_offset, molecule_noncovalent_count),
                    offset_map(stereo_atom_offset, molecule_stereo_atom_count),
                    offset_map(stereo_bond_offset, molecule_stereo_bond_count),
                );
                for constraint in molecule.constraints.iter() {
                    entries
                        .constraints
                        .push(constraint.clone().remap(&remapping));
                }
            }

            correspondences.push(MoleculeCorrespondence::new(
                offset_correspondence(atom_offset, molecule_atom_count, atom_count),
                offset_correspondence(bond_offset, molecule_bond_count, bond_count),
                offset_correspondence(dative_offset, molecule_dative_count, dative_count),
                offset_correspondence(aromatic_offset, molecule_aromatic_count, aromatic_count),
                offset_correspondence(
                    multicenter_offset,
                    molecule_multicenter_count,
                    multicenter_count,
                ),
                offset_correspondence(
                    noncovalent_offset,
                    molecule_noncovalent_count,
                    noncovalent_count,
                ),
                offset_correspondence(
                    stereo_atom_offset,
                    molecule_stereo_atom_count,
                    stereo_atom_count,
                ),
                offset_correspondence(
                    stereo_bond_offset,
                    molecule_stereo_bond_count,
                    stereo_bond_count,
                ),
            ));

            atom_offset += molecule_atom_count;
            bond_offset += molecule_bond_count;
            dative_offset += molecule_dative_count;
            aromatic_offset += molecule_aromatic_count;
            multicenter_offset += molecule_multicenter_count;
            noncovalent_offset += molecule_noncovalent_count;
            stereo_atom_offset += molecule_stereo_atom_count;
            stereo_bond_offset += molecule_stereo_bond_count;
        }

        (Molecule::from_entries(entries), correspondences)
    }

    /// Combine `self` and `other` as a fresh molecule by disjoint concatenation. Returns the
    /// correspondence mapping `other` into the combined molecule.
    pub fn combine(&self, other: &Molecule) -> (Molecule, MoleculeCorrespondence) {
        let (combined, mut correspondences) = Self::combine_all([self, other]);
        let correspondence = correspondences
            .pop()
            .expect("two inputs produce two correspondences");
        (combined, correspondence)
    }

    /// Append `other` by disjoint concatenation. `self` keeps its ids as the prefix; `other`'s
    /// entities follow. Returns the correspondence mapping `other` into the combined molecule.
    pub fn combine_from(&mut self, other: &Molecule) -> MoleculeCorrespondence {
        let atom_offset = self.atoms().count();
        let bond_offset = self.bonds().count();
        let dative_offset = self.dative_bonds().count();
        let aromatic_offset = self.aromatic_systems().count();
        let multicenter_offset = self.multicenter_bonds().count();
        let noncovalent_offset = self.noncovalent_bonds().count();
        let stereo_atom_offset = self.stereo_atoms().count();
        let stereo_bond_offset = self.stereo_bonds().count();
        let shift_atom = |id: AtomId| AtomId(id.0 + atom_offset as u32);

        let Molecule {
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
        } = mem::take(self);
        let mut editor = MoleculeEditor::from_parts(
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
        );

        for atom in other.atoms().iter() {
            editor.add_atom(atom.attributes.clone());
        }
        for bond in other.bonds().iter() {
            let [first, second] = bond.atom_ids();
            editor.add_bond(
                shift_atom(first),
                shift_atom(second),
                bond.attributes.clone(),
            );
        }
        for bond in other.dative_bonds().iter() {
            editor.add_dative_bond(
                bond.donors().map(|donor| shift_atom(donor.id)).collect(),
                shift_atom(bond.acceptor_id()),
                bond.attributes.clone(),
            );
        }
        for system in other.aromatic_systems().iter() {
            editor.add_aromatic_system(
                system.atom_ids().map(shift_atom).collect(),
                system.attributes.clone(),
            );
        }
        for bond in other.multicenter_bonds().iter() {
            editor.add_multicenter_bond(
                bond.atom_ids().map(shift_atom).collect(),
                bond.attributes.clone(),
            );
        }
        for bond in other.noncovalent_bonds().iter() {
            let [first, second] = bond.atom_ids();
            editor.add_noncovalent_bond(
                [shift_atom(first), shift_atom(second)],
                bond.attributes.clone(),
            );
        }

        let ligand_remapping = Remapping::new(
            (0..other.atoms().count())
                .map(|index| NodeId((atom_offset + index) as u32))
                .collect(),
            (0..other.bonds().count())
                .map(|index| EdgeId((bond_offset + index) as u32))
                .collect(),
        );
        for id in other.stereo_atoms.ids() {
            let site = shift_atom(other.stereo_atoms.site(id));
            let ligands = other
                .stereo_atoms
                .ligands(id)
                .iter()
                .map(|ligand| ligand.remap(&ligand_remapping))
                .collect();
            editor.add_stereo_atom(site, ligands, other.stereo_atoms.attributes(id).clone());
        }
        for id in other.stereo_bonds.ids() {
            let site = BondId(other.stereo_bonds.site(id).0 + bond_offset as u32);
            let ligands = other
                .stereo_bonds
                .ligands(id)
                .iter()
                .map(|ligand| ligand.remap(&ligand_remapping))
                .collect();
            editor.add_stereo_bond(site, ligands, other.stereo_bonds.attributes(id).clone());
        }

        if !other.constraints.is_empty() {
            let remapping = IdRemapping::new(
                offset_map(atom_offset, other.atoms().count()),
                offset_map(bond_offset, other.bonds().count()),
                offset_map(dative_offset, other.dative_bonds().count()),
                offset_map(aromatic_offset, other.aromatic_systems().count()),
                offset_map(multicenter_offset, other.multicenter_bonds().count()),
                offset_map(noncovalent_offset, other.noncovalent_bonds().count()),
                offset_map(stereo_atom_offset, other.stereo_atoms().count()),
                offset_map(stereo_bond_offset, other.stereo_bonds().count()),
            );
            for constraint in other.constraints.iter() {
                editor
                    .constraints_mut()
                    .push(constraint.clone().remap(&remapping));
            }
        }
        *self = editor.build();

        MoleculeCorrespondence::new(
            offset_correspondence(
                atom_offset,
                other.atoms().count(),
                atom_offset + other.atoms().count(),
            ),
            offset_correspondence(
                bond_offset,
                other.bonds().count(),
                bond_offset + other.bonds().count(),
            ),
            offset_correspondence(
                dative_offset,
                other.dative_bonds().count(),
                dative_offset + other.dative_bonds().count(),
            ),
            offset_correspondence(
                aromatic_offset,
                other.aromatic_systems().count(),
                aromatic_offset + other.aromatic_systems().count(),
            ),
            offset_correspondence(
                multicenter_offset,
                other.multicenter_bonds().count(),
                multicenter_offset + other.multicenter_bonds().count(),
            ),
            offset_correspondence(
                noncovalent_offset,
                other.noncovalent_bonds().count(),
                noncovalent_offset + other.noncovalent_bonds().count(),
            ),
            offset_correspondence(
                stereo_atom_offset,
                other.stereo_atoms().count(),
                stereo_atom_offset + other.stereo_atoms().count(),
            ),
            offset_correspondence(
                stereo_bond_offset,
                other.stereo_bonds().count(),
                stereo_bond_offset + other.stereo_bonds().count(),
            ),
        )
    }

    /// Decompose into connected components — a conservative partition where every relation keeps its
    /// atoms in one component (a spanning overlay prevents the split rather than being dropped). Each
    /// component is a fresh, compactly-renumbered `Molecule` paired with the
    /// `MoleculeCorrespondence` mapping its ids back to `self`. Components are ordered by their lowest
    /// original atom.
    pub fn split(&self) -> Vec<(Molecule, MoleculeCorrespondence)> {
        let atom_count = self.atoms().count();
        let mut uf = UnionFind::new(atom_count);
        for bond in self.bonds().iter() {
            let [a, b] = bond.atom_ids();
            uf.union(a.index(), b.index());
        }
        for dative in self.dative_bonds().iter() {
            union_participants(&mut uf, dative.atom_ids());
        }
        for system in self.aromatic_systems().iter() {
            union_participants(&mut uf, system.atom_ids());
        }
        for bond in self.multicenter_bonds().iter() {
            union_participants(&mut uf, bond.atom_ids());
        }
        for bond in self.noncovalent_bonds().iter() {
            let [a, b] = bond.atom_ids();
            uf.union(a.index(), b.index());
        }
        for rid in self.stereo_atoms.ids() {
            let site = self.stereo_atoms.site(rid);
            for ligand in self.stereo_atoms.ligands(rid) {
                uf.union(site.index(), ligand.atom_id.index());
            }
        }
        for rid in self.stereo_bonds.ids() {
            let [a, b] = self.bond(BondId(self.stereo_bonds.site(rid).0)).atom_ids();
            uf.union(a.index(), b.index());
            for ligand in self.stereo_bonds.ligands(rid) {
                uf.union(a.index(), ligand.atom_id.index());
            }
        }
        for constraint in self.constraints.iter() {
            union_participants(&mut uf, self.constraint_atoms(constraint));
        }

        let mut atom_component = vec![0usize; atom_count];
        let mut atom_compact = vec![0u32; atom_count];
        let mut component_atoms: Vec<Vec<AtomId>> = Vec::new();
        let mut index_of_root: HashMap<usize, usize> = HashMap::new();
        for i in 0..atom_count {
            let root = uf.find(i);
            let component = *index_of_root.entry(root).or_insert_with(|| {
                component_atoms.push(Vec::new());
                component_atoms.len() - 1
            });
            atom_component[i] = component;
            atom_compact[i] = component_atoms[component].len() as u32;
            component_atoms[component].push(AtomId(i as u32));
        }
        let compact = |a: AtomId| AtomId(atom_compact[a.index()]);
        let component_of = |a: AtomId| atom_component[a.index()];
        let compaction = Remapping::new(
            (0..atom_count).map(|i| NodeId(atom_compact[i])).collect(),
            Vec::new(),
        );

        component_atoms
            .iter()
            .enumerate()
            .map(|(component, atoms)| {
                let mut editor = Molecule::new().edit();
                let mut atom_pairs = Vec::new();
                for atom in atoms {
                    let added = editor.add_atom(self.atom(*atom).attributes.clone());
                    atom_pairs.push((added, *atom));
                }
                let mut bond_pairs = Vec::new();
                let mut bond_compact: HashMap<BondId, BondId> = HashMap::new();
                for bond in self.bonds().iter() {
                    let [a, b] = bond.atom_ids();
                    if component_of(a) == component {
                        let new_bond =
                            editor.add_bond(compact(a), compact(b), bond.attributes.clone());
                        bond_pairs.push((new_bond, bond.id));
                        bond_compact.insert(bond.id, new_bond);
                    }
                }
                let mut dative_pairs = Vec::new();
                for dative in self.dative_bonds().iter() {
                    if component_of(dative.acceptor_id()) == component {
                        let donors = dative.donors().map(|d| compact(d.id)).collect();
                        let added = editor.add_dative_bond(
                            donors,
                            compact(dative.acceptor_id()),
                            dative.attributes.clone(),
                        );
                        dative_pairs.push((added, dative.id));
                    }
                }
                let mut aromatic_pairs = Vec::new();
                for system in self.aromatic_systems().iter() {
                    let members: Vec<AtomId> = system.atom_ids().collect();
                    if members
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        let added = editor.add_aromatic_system(
                            members.iter().map(|a| compact(*a)).collect(),
                            system.attributes.clone(),
                        );
                        aromatic_pairs.push((added, system.id));
                    }
                }
                let mut multicenter_pairs = Vec::new();
                for bond in self.multicenter_bonds().iter() {
                    let members: Vec<AtomId> = bond.atom_ids().collect();
                    if members
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        let added = editor.add_multicenter_bond(
                            members.iter().map(|a| compact(*a)).collect(),
                            bond.attributes.clone(),
                        );
                        multicenter_pairs.push((added, bond.id));
                    }
                }
                let mut noncovalent_pairs = Vec::new();
                for bond in self.noncovalent_bonds().iter() {
                    let [a, b] = bond.atom_ids();
                    if component_of(a) == component {
                        let added = editor.add_noncovalent_bond(
                            [compact(a), compact(b)],
                            bond.attributes.clone(),
                        );
                        noncovalent_pairs.push((added, bond.id));
                    }
                }
                let mut stereo_atom_pairs = Vec::new();
                for rid in self.stereo_atoms.ids() {
                    let site = self.stereo_atoms.site(rid);
                    if component_of(site) == component {
                        let ligands: Vec<StereoLigand> = self
                            .stereo_atoms
                            .ligands(rid)
                            .iter()
                            .map(|ligand| ligand.remap(&compaction))
                            .collect();
                        let added = editor.add_stereo_atom(
                            compact(site),
                            ligands,
                            self.stereo_atoms.attributes(rid).clone(),
                        );
                        stereo_atom_pairs.push((added, rid));
                    }
                }
                let mut stereo_bond_pairs = Vec::new();
                for rid in self.stereo_bonds.ids() {
                    let bond = BondId(self.stereo_bonds.site(rid).0);
                    let [a, _] = self.bond(bond).atom_ids();
                    if component_of(a) == component {
                        let ligands: Vec<StereoLigand> = self
                            .stereo_bonds
                            .ligands(rid)
                            .iter()
                            .map(|ligand| ligand.remap(&compaction))
                            .collect();
                        let added = editor.add_stereo_bond(
                            bond_compact[&bond],
                            ligands,
                            self.stereo_bonds.attributes(rid).clone(),
                        );
                        stereo_bond_pairs.push((added, rid));
                    }
                }
                let entities = editor.build();

                let correspondence = MoleculeCorrespondence::new(
                    Correspondence::new(atom_pairs, entities.atoms().count(), atom_count)
                        .expect("split assigns each component atom its original id"),
                    Correspondence::new(bond_pairs, entities.bonds().count(), self.bonds().count())
                        .expect("split assigns each component bond its original id"),
                    Correspondence::new(
                        dative_pairs,
                        entities.dative_bonds().count(),
                        self.dative_bonds().count(),
                    )
                    .expect("split assigns each component dative bond its original id"),
                    Correspondence::new(
                        aromatic_pairs,
                        entities.aromatic_systems().count(),
                        self.aromatic_systems().count(),
                    )
                    .expect("split assigns each component aromatic system its original id"),
                    Correspondence::new(
                        multicenter_pairs,
                        entities.multicenter_bonds().count(),
                        self.multicenter_bonds().count(),
                    )
                    .expect("split assigns each component multicenter bond its original id"),
                    Correspondence::new(
                        noncovalent_pairs,
                        entities.noncovalent_bonds().count(),
                        self.noncovalent_bonds().count(),
                    )
                    .expect("split assigns each component noncovalent bond its original id"),
                    Correspondence::new(
                        stereo_atom_pairs,
                        entities.stereo_atoms().count(),
                        self.stereo_atoms().count(),
                    )
                    .expect("split assigns each component stereo atom its original id"),
                    Correspondence::new(
                        stereo_bond_pairs,
                        entities.stereo_bonds().count(),
                        self.stereo_bonds().count(),
                    )
                    .expect("split assigns each component stereo bond its original id"),
                );

                // route each constraint to the component holding its atoms (conservative union
                // guarantees they share one), remapped to compact ids
                let remapping = idremapping_from_correspondence(&correspondence);
                let mut editor = entities.edit();
                for constraint in self.constraints.iter() {
                    if self
                        .constraint_atoms(constraint)
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        editor
                            .constraints_mut()
                            .push(constraint.clone().remap(&remapping));
                    }
                }
                (editor.build(), correspondence)
            })
            .collect()
    }

    /// Every atom a constraint transitively references — the conservative binding set `split` unions
    /// and the connectivity validator checks (mirrors the ids each `Constraint` arm carries, resolved
    /// to atoms; a `None`-scoped molecule constraint binds the whole molecule).
    pub fn constraint_atoms(&self, constraint: &Constraint) -> Vec<AtomId> {
        let all_atoms = || {
            (0..self.atoms().count() as u32)
                .map(AtomId)
                .collect::<Vec<_>>()
        };
        match constraint {
            Constraint::Atom(id, _) => vec![*id],
            Constraint::Bond(id, _) => self.bond(*id).atom_ids().to_vec(),
            Constraint::DativeBond(id, _) => self.dative_bond(*id).atom_ids().collect(),
            Constraint::AromaticSystem(id, _) => self.aromatic_system(*id).atom_ids().collect(),
            Constraint::MulticenterBond(id, _) => self.multicenter_bond(*id).atom_ids().collect(),
            Constraint::NoncovalentBond(id, _) => self.noncovalent_bond(*id).atom_ids().to_vec(),
            Constraint::StereoAtom(id, _, _) => {
                let view = self.stereo_atom(*id);
                iter::once(view.site_id())
                    .chain(view.ligands().map(|ligand| ligand.atom_id()))
                    .collect()
            }
            Constraint::StereoBond(id, _, _) => {
                let view = self.stereo_bond(*id);
                let [a, b] = self.bond(view.site_id()).atom_ids();
                [a, b]
                    .into_iter()
                    .chain(view.ligands().map(|ligand| ligand.atom_id()))
                    .collect()
            }
            Constraint::Relational(relational) => self.relational_constraint_atoms(relational),
            Constraint::Molecule(molecule) => match molecule {
                MoleculeConstraint::ChargeSum { atoms, .. }
                | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
                | MoleculeConstraint::Connected { atoms } => {
                    atoms.clone().unwrap_or_else(all_atoms)
                }
                MoleculeConstraint::BondOrderSum { bonds, .. } => match bonds {
                    Some(bonds) => bonds
                        .iter()
                        .flat_map(|b| self.bond(*b).atom_ids())
                        .collect(),
                    None => all_atoms(),
                },
            },
            Constraint::And(constraints) | Constraint::Or(constraints) => constraints
                .iter()
                .flat_map(|c| self.constraint_atoms(c))
                .collect(),
            Constraint::Not(constraint) => self.constraint_atoms(constraint),
        }
    }

    /// A relational constraint's atoms: its primary entity's atoms plus any explicit atom operands.
    fn relational_constraint_atoms(&self, constraint: &RelationalConstraint) -> Vec<AtomId> {
        let dative = |id, extra: &[AtomId]| {
            self.dative_bond(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let aromatic = |id, extra: &[AtomId]| {
            self.aromatic_system(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let multicenter = |id, extra: &[AtomId]| {
            self.multicenter_bond(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let noncovalent = |id, extra: &[AtomId]| {
            self.noncovalent_bond(id)
                .atom_ids()
                .into_iter()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let stereo_atom = |id, extra: &[AtomId]| {
            let view = self.stereo_atom(id);
            iter::once(view.site_id())
                .chain(view.ligands().map(|ligand| ligand.atom_id()))
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let stereo_bond = |id, extra: &[AtomId]| {
            let view = self.stereo_bond(id);
            let [a, b] = self.bond(view.site_id()).atom_ids();
            [a, b]
                .into_iter()
                .chain(view.ligands().map(|ligand| ligand.atom_id()))
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        match constraint {
            RelationalConstraint::DativeBondDonors { bond, atoms } => dative(*bond, atoms),
            RelationalConstraint::DativeBondDonor { bond, atom } => dative(*bond, &[*atom]),
            RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                dative(*bond, atoms)
            }
            RelationalConstraint::DativeBondAllDonors { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondAnyDonor { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondAcceptor { bond, atom } => dative(*bond, &[*atom]),
            RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondParallels {
                dative: id,
                parallel,
            } => dative(*id, &self.bond(*parallel).atom_ids()),
            RelationalConstraint::AromaticSystemAtoms { system, atoms } => aromatic(*system, atoms),
            RelationalConstraint::AromaticSystemContains { system, atom } => {
                aromatic(*system, &[*atom])
            }
            RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                aromatic(*system, atoms)
            }
            RelationalConstraint::AromaticSystemAllAtoms { system, .. } => aromatic(*system, &[]),
            RelationalConstraint::AromaticSystemAnyAtom { system, .. } => aromatic(*system, &[]),
            RelationalConstraint::MulticenterBondAtoms { bond, atoms } => multicenter(*bond, atoms),
            RelationalConstraint::MulticenterBondContains { bond, atom } => {
                multicenter(*bond, &[*atom])
            }
            RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                multicenter(*bond, atoms)
            }
            RelationalConstraint::MulticenterBondAllAtoms { bond, .. } => multicenter(*bond, &[]),
            RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => multicenter(*bond, &[]),
            RelationalConstraint::NoncovalentBondEnds { bond, atoms } => noncovalent(*bond, atoms),
            RelationalConstraint::NoncovalentBondContains { bond, atom } => {
                noncovalent(*bond, &[*atom])
            }
            RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
                noncovalent(*bond, &[])
            }
            RelationalConstraint::StereoAtomSite {
                stereo_atom: id,
                atom,
            } => stereo_atom(*id, &[*atom]),
            RelationalConstraint::StereoAtomContains {
                stereo_atom: id,
                atom,
            } => stereo_atom(*id, &[*atom]),
            RelationalConstraint::StereoAtomLigands {
                stereo_atom: id,
                atoms,
            } => stereo_atom(*id, atoms),
            RelationalConstraint::StereoAtomAllLigands {
                stereo_atom: id, ..
            } => stereo_atom(*id, &[]),
            RelationalConstraint::StereoAtomAnyLigand {
                stereo_atom: id, ..
            } => stereo_atom(*id, &[]),
            RelationalConstraint::StereoBondSite {
                stereo_bond: id,
                bond,
            } => stereo_bond(*id, &self.bond(*bond).atom_ids()),
            RelationalConstraint::StereoBondContains {
                stereo_bond: id,
                atom,
            } => stereo_bond(*id, &[*atom]),
            RelationalConstraint::StereoBondLigands {
                stereo_bond: id,
                atoms,
            } => stereo_bond(*id, atoms),
            RelationalConstraint::StereoBondAllLigands {
                stereo_bond: id, ..
            } => stereo_bond(*id, &[]),
            RelationalConstraint::StereoBondAnyLigand {
                stereo_bond: id, ..
            } => stereo_bond(*id, &[]),
        }
    }
}

impl Normalize for Molecule {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for attributes in Arc::make_mut(&mut self.atoms) {
            *attributes = mem::take(attributes).normalize()?;
        }
        for attributes in Arc::make_mut(&mut self.bonds) {
            *attributes = mem::take(attributes).normalize()?;
        }
        self.dative_bonds = mem::take(&mut self.dative_bonds).normalize()?;
        self.aromatic_systems = mem::take(&mut self.aromatic_systems).normalize()?;
        self.multicenter_bonds = mem::take(&mut self.multicenter_bonds).normalize()?;
        self.noncovalent_bonds = mem::take(&mut self.noncovalent_bonds).normalize()?;
        self.stereo_atoms = mem::take(&mut self.stereo_atoms).normalize()?;
        self.stereo_bonds = mem::take(&mut self.stereo_bonds).normalize()?;
        self.constraints = mem::take(&mut self.constraints).normalize()?;
        Ok(self)
    }
}

impl FrameTransport for Molecule {
    type Action = OverlaysFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        self.dative_bonds = self.dative_bonds.reframe_by(actions.dative_bonds())?;
        self.aromatic_systems = self
            .aromatic_systems
            .reframe_by(actions.aromatic_systems())?;
        self.multicenter_bonds = self
            .multicenter_bonds
            .reframe_by(actions.multicenter_bonds())?;
        self.noncovalent_bonds = self
            .noncovalent_bonds
            .reframe_by(actions.noncovalent_bonds())?;
        self.stereo_atoms = self.stereo_atoms.reframe_by(actions.stereo_atoms())?;
        self.stereo_bonds = self.stereo_bonds.reframe_by(actions.stereo_bonds())?;
        self.constraints = self.constraints.reframe_by(actions)?;
        Some(self)
    }
}

impl Reframe for Molecule {
    fn representative_action(&self) -> Self::Action {
        OverlaysFrameAction::new(
            self.dative_bonds.representative_action(),
            self.aromatic_systems.representative_action(),
            self.multicenter_bonds.representative_action(),
            self.noncovalent_bonds.representative_action(),
            self.stereo_atoms.representative_action(),
            self.stereo_bonds.representative_action(),
        )
    }

    fn reframe(mut self) -> Result<Self, Contradiction> {
        let action_domain = self.constraints.frame_action_domain();
        for attributes in Arc::make_mut(&mut self.atoms) {
            *attributes = mem::take(attributes).normalize()?;
        }
        for attributes in Arc::make_mut(&mut self.bonds) {
            *attributes = mem::take(attributes).normalize()?;
        }
        self.constraints = mem::take(&mut self.constraints).normalize()?;

        let mut actions = ConstraintFrameActionMap::default();
        self.dative_bonds = if action_domain.count(EntityKind::DativeBond) == 0 {
            mem::take(&mut self.dative_bonds).reframe()?
        } else {
            reframe_dative_bonds_with(mem::take(&mut self.dative_bonds), |id, action| {
                if action_domain.contains_dative_bond(id) {
                    actions.insert_dative_bond(id, action.clone());
                }
            })?
        };
        self.aromatic_systems = if action_domain.count(EntityKind::AromaticSystem) == 0 {
            mem::take(&mut self.aromatic_systems).reframe()?
        } else {
            reframe_aromatic_systems_with(mem::take(&mut self.aromatic_systems), |id, action| {
                if action_domain.contains_aromatic_system(id) {
                    actions.insert_aromatic_system(id, action.clone());
                }
            })?
        };
        self.multicenter_bonds = if action_domain.count(EntityKind::MulticenterBond) == 0 {
            mem::take(&mut self.multicenter_bonds).reframe()?
        } else {
            reframe_multicenter_bonds_with(mem::take(&mut self.multicenter_bonds), |id, action| {
                if action_domain.contains_multicenter_bond(id) {
                    actions.insert_multicenter_bond(id, action.clone());
                }
            })?
        };
        self.noncovalent_bonds = if action_domain.count(EntityKind::NoncovalentBond) == 0 {
            mem::take(&mut self.noncovalent_bonds).reframe()?
        } else {
            reframe_noncovalent_bonds_with(mem::take(&mut self.noncovalent_bonds), |id, action| {
                if action_domain.contains_noncovalent_bond(id) {
                    actions.insert_noncovalent_bond(id, action.clone());
                }
            })?
        };
        self.stereo_atoms = if action_domain.count(EntityKind::StereoAtom) == 0 {
            mem::take(&mut self.stereo_atoms).reframe()?
        } else {
            reframe_stereo_atoms_with(mem::take(&mut self.stereo_atoms), |id, action| {
                if action_domain.contains_stereo_atom(id) {
                    actions.insert_stereo_atom(id, action);
                }
            })?
        };
        self.stereo_bonds = if action_domain.count(EntityKind::StereoBond) == 0 {
            mem::take(&mut self.stereo_bonds).reframe()?
        } else {
            reframe_stereo_bonds_with(mem::take(&mut self.stereo_bonds), |id, action| {
                if action_domain.contains_stereo_bond(id) {
                    actions.insert_stereo_bond(id, action);
                }
            })?
        };
        self.constraints = self
            .constraints
            .reframe_by_actions(&actions)
            .map_err(|_| Contradiction)?
            .normalize()?;
        Ok(self)
    }
}

/// The correspondence mapping one input entity set's ids to their offset ids in a combined set.
fn offset_correspondence<Id: Copy + Ord + From<usize>>(
    offset: usize,
    input_count: usize,
    combined_count: usize,
) -> Correspondence<Id> {
    let images: Vec<Id> = (0..input_count).map(|k| Id::from(offset + k)).collect();
    Correspondence::from_images(&images, combined_count)
}

/// The per-entity-kind offset map used to remap constraints into a combined molecule.
fn offset_map<Id: Copy + Eq + Hash + From<usize>>(offset: usize, count: usize) -> HashMap<Id, Id> {
    (0..count)
        .map(|k| (Id::from(k), Id::from(offset + k)))
        .collect()
}

/// Union every atom of a relation into one component (all participants share the first's set).
fn union_participants(uf: &mut UnionFind, atoms: impl IntoIterator<Item = AtomId>) {
    let mut atoms = atoms.into_iter();
    if let Some(first) = atoms.next() {
        for atom in atoms {
            uf.union(first.index(), atom.index());
        }
    }
}

/// The per-entity-kind `original → compact` remapping a `split` component induces, read off its
/// `component → original` correspondence.
fn idremapping_from_correspondence(correspondence: &MoleculeCorrespondence) -> IdRemapping {
    IdRemapping::new(
        correspondence
            .atoms()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .bonds()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .dative_bonds()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .aromatic_systems()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .multicenter_bonds()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .noncovalent_bonds()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .stereo_atoms()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .stereo_bonds()
            .matched_pairs()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
    )
}

#[cfg(test)]
mod tests;
