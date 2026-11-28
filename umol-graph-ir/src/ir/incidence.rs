//! Incidence graph: localized bonds stay direct edges while overlays become pseudonodes.
//!
//! This is an internal graph-IR lowering onto `umol_graph_core::Graph` for symmetry analysis,
//! canonicalization, and `SubstructureMatchAlgorithm::Incidence`. It is not a bipartite Levi
//! graph: atoms occupy `0..atom_count`, included overlays occupy `atom_count..node_count`,
//! localized bonds occupy edges `0..bond_count`, and overlay–atom incidences occupy edges
//! `bond_count..edge_count`.

use bitflags::bitflags;
use umol_graph_core::{EdgeId, Graph, NodeId};

use super::entity::{Entity, EntityKind};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigandKind;
use super::molecule::Molecule;
use super::reaction_span::ReactionSpan;

/// Dense index of a node in an [`IncidenceGraph`] (atoms in `0..atom_count`, followed by
/// included overlay entities).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub(crate) struct IncidenceNodeId(pub(crate) u32);

impl IncidenceNodeId {
    #[inline]
    pub(crate) fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<usize> for IncidenceNodeId {
    #[inline]
    fn from(value: usize) -> Self {
        Self(value as u32)
    }
}

impl From<IncidenceNodeId> for usize {
    #[inline]
    fn from(id: IncidenceNodeId) -> Self {
        id.0 as usize
    }
}

impl From<NodeId> for IncidenceNodeId {
    #[inline]
    fn from(id: NodeId) -> Self {
        Self(id.0)
    }
}

impl From<IncidenceNodeId> for NodeId {
    #[inline]
    fn from(id: IncidenceNodeId) -> Self {
        Self(id.0)
    }
}

/// Dense index of an edge in an [`IncidenceGraph`] (localized bonds in `0..bond_count`, followed
/// by overlay–atom incidences in `bond_count..edge_count`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub(crate) struct IncidenceEdgeId(pub(crate) u32);

impl IncidenceEdgeId {
    #[inline]
    pub(crate) fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<usize> for IncidenceEdgeId {
    #[inline]
    fn from(value: usize) -> Self {
        Self(value as u32)
    }
}

impl From<IncidenceEdgeId> for usize {
    #[inline]
    fn from(id: IncidenceEdgeId) -> Self {
        id.0 as usize
    }
}

impl From<EdgeId> for IncidenceEdgeId {
    #[inline]
    fn from(id: EdgeId) -> Self {
        Self(id.0)
    }
}

impl From<IncidenceEdgeId> for EdgeId {
    #[inline]
    fn from(id: IncidenceEdgeId) -> Self {
        Self(id.0)
    }
}

/// Structural level represented by an [`IncidenceGraph`].
///
/// The variants form a nested hierarchy. [`Topology`](Self::Topology) contains
/// atoms and localized bonds (no overlay pseudonodes). [`Constitution`](Self::Constitution)
/// additionally contains dative bonds, aromatic systems, multicenter bonds, and noncovalent
/// bonds. [`Full`](Self::Full) additionally contains stereo atoms and stereo
/// bonds. Constraints do not contribute to any structural level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum IncidenceLevel {
    Topology,
    Constitution,
    Full,
}

bitflags! {
    /// Virtual stereo ligand kinds anchored at a stereo-atom site or stereo-bond endpoint.
    ///
    /// Each bit records whether the corresponding virtual [`StereoLigandKind`] occurs at that
    /// anchor.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub(crate) struct StereoVirtualLigands: u8 {
        const IMPLICIT_HYDROGEN = 1 << 0;
        const LONE_PAIR = 1 << 1;
        const VIRTUAL = 1 << 2;
    }
}

impl StereoVirtualLigands {
    fn add_kind(&mut self, kind: StereoLigandKind) {
        match kind {
            StereoLigandKind::Atom => {}
            StereoLigandKind::ImplicitHydrogen => self.insert(Self::IMPLICIT_HYDROGEN),
            StereoLigandKind::LonePair => self.insert(Self::LONE_PAIR),
            StereoLigandKind::Virtual => self.insert(Self::VIRTUAL),
        }
    }
}

/// The chemical role of an overlay–atom incidence edge in an [`IncidenceGraph`].
///
/// Localized bonds are direct atom–atom edges represented by [`IncidenceEdge::Bond`], not
/// variants of this enum. Stereo entities encode virtual ligands anchored at their site atom or
/// bond endpoints on [`StereoAnchor`](Self::StereoAnchor) edges and actual-atom ligands on
/// [`StereoLigand`](Self::StereoLigand) edges. Electron contributions and reaction-side
/// attributes remain in the owning entity forms; they are not incidence labels.
///
/// # Semantic properties
///
/// The total order follows the declaration order of the variants. Canonicalization
/// augments these roles with the owning forms' electron contributions when
/// constructing edge colors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Incidence {
    DativeDonor,
    DativeAcceptor,
    AromaticAtom,
    MulticenterAtom,
    NoncovalentEndpoint,
    StereoAnchor(StereoVirtualLigands),
    StereoLigand,
}

/// Classification of one edge in an [`IncidenceGraph`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum IncidenceEdge {
    /// A direct atom–atom edge corresponding to a localized bond.
    Bond(BondId),
    /// An overlay–atom edge carrying its participant role.
    Incidence(Incidence),
}

/// Lowered simple graph of a molecule or reaction span for symmetry, canonicalization, and
/// incidence substructure matching.
///
/// Atoms occupy `0..atom_count` and included overlays follow in contiguous [`EntityKind`] order
/// (`DativeBond`, `AromaticSystem`, `MulticenterBond`, `NoncovalentBond`, `StereoAtom`,
/// `StereoBond`). Localized bonds are direct edges in `0..bond_count`; overlay–atom incidences
/// occupy edges `bond_count..edge_count`.
#[derive(Clone, Debug)]
pub(crate) struct IncidenceGraph {
    graph: Graph,
    // End offsets for `[Atom, DativeBond, AromaticSystem, MulticenterBond, NoncovalentBond,
    // StereoAtom, StereoBond]`.
    offsets: [u32; 7],
    bond_count: u32,
    // Overlay–atom incidences for edges `bond_count..graph.edge_count()`.
    incidences: Vec<Incidence>,
}

impl IncidenceGraph {
    /// Underlying `umol_graph_core::Graph` (boundary for graph-core algorithms).
    #[inline]
    pub(crate) fn raw_graph(&self) -> &Graph {
        &self.graph
    }

    /// Number of atom nodes (`0..atom_count`).
    #[inline]
    pub(crate) fn atom_count(&self) -> usize {
        self.offsets[0] as usize
    }

    /// Total number of nodes (atoms + included overlays).
    #[inline]
    pub(crate) fn node_count(&self) -> usize {
        self.offsets[6] as usize
    }

    /// Number of localized bond edges (`0..bond_count`).
    #[inline]
    pub(crate) fn bond_count(&self) -> usize {
        self.bond_count as usize
    }

    /// Number of overlay–atom incidence edges (`bond_count..edge_count`).
    #[inline]
    pub(crate) fn incidence_count(&self) -> usize {
        self.incidences.len()
    }

    /// Total number of edges (`bond_count + incidence_count`).
    #[inline]
    pub(crate) fn edge_count(&self) -> usize {
        self.bond_count as usize + self.incidences.len()
    }

    /// Contiguous `[start, end)` node index range for `kind`.
    ///
    /// Localized bonds are edges rather than nodes and therefore have an empty range `(0, 0)`.
    #[inline]
    pub(crate) fn node_range(&self, kind: EntityKind) -> (usize, usize) {
        match kind {
            EntityKind::Atom => (0, self.offsets[0] as usize),
            EntityKind::Bond => (0, 0),
            EntityKind::DativeBond => (self.offsets[0] as usize, self.offsets[1] as usize),
            EntityKind::AromaticSystem => (self.offsets[1] as usize, self.offsets[2] as usize),
            EntityKind::MulticenterBond => (self.offsets[2] as usize, self.offsets[3] as usize),
            EntityKind::NoncovalentBond => (self.offsets[3] as usize, self.offsets[4] as usize),
            EntityKind::StereoAtom => (self.offsets[4] as usize, self.offsets[5] as usize),
            EntityKind::StereoBond => (self.offsets[5] as usize, self.offsets[6] as usize),
        }
    }

    /// Number of nodes of the given entity kind (`0` for [`EntityKind::Bond`]).
    #[inline]
    pub(crate) fn node_count_of(&self, kind: EntityKind) -> usize {
        let (start, end) = self.node_range(kind);
        end - start
    }

    /// Every node id in ascending order.
    #[inline]
    pub(crate) fn node_ids(&self) -> impl ExactSizeIterator<Item = IncidenceNodeId> {
        (0..self.offsets[6]).map(IncidenceNodeId)
    }

    /// The molecule entity a node represents (an atom or an included overlay, never a bond).
    ///
    /// # Panics
    ///
    /// Panics when `node` is out of range.
    #[inline]
    pub(crate) fn entity(&self, node: IncidenceNodeId) -> Entity {
        let n = node.0;
        if n < self.offsets[0] {
            return Entity::Atom(AtomId(n));
        }
        if n < self.offsets[1] {
            return Entity::DativeBond(DativeBondId(n - self.offsets[0]));
        }
        if n < self.offsets[2] {
            return Entity::AromaticSystem(AromaticSystemId(n - self.offsets[1]));
        }
        if n < self.offsets[3] {
            return Entity::MulticenterBond(MulticenterBondId(n - self.offsets[2]));
        }
        if n < self.offsets[4] {
            return Entity::NoncovalentBond(NoncovalentBondId(n - self.offsets[3]));
        }
        if n < self.offsets[5] {
            return Entity::StereoAtom(StereoAtomId(n - self.offsets[4]));
        }
        if n < self.offsets[6] {
            return Entity::StereoBond(StereoBondId(n - self.offsets[5]));
        }
        panic!("incidence node {n} out of range");
    }

    /// The node representing `entity`, or `None` for [`Entity::Bond`], kinds excluded by this
    /// graph's [`IncidenceLevel`], or out-of-range entity ids.
    #[inline]
    pub(crate) fn node_of(&self, entity: Entity) -> Option<IncidenceNodeId> {
        let (start, end) = self.node_range(entity.kind());
        let id = entity.kind_id();
        (id < end - start).then_some(IncidenceNodeId((start + id) as u32))
    }

    /// Classification of `edge` as either a localized bond or an overlay–atom incidence.
    ///
    /// # Panics
    ///
    /// Panics when `edge` is out of range.
    #[inline]
    pub(crate) fn edge(&self, edge: IncidenceEdgeId) -> IncidenceEdge {
        if edge.0 < self.bond_count {
            IncidenceEdge::Bond(BondId(edge.0))
        } else {
            IncidenceEdge::Incidence(self.incidences[(edge.0 - self.bond_count) as usize])
        }
    }

    /// The edge id corresponding to localized bond `bond`, or `None` if out of range.
    #[inline]
    pub(crate) fn bond_edge(&self, bond: BondId) -> Option<IncidenceEdgeId> {
        (bond.0 < self.bond_count).then_some(IncidenceEdgeId(bond.0))
    }

    /// Dense `0..incidence_count()` index of an overlay incidence edge, or `None` for localized
    /// bond edges or out-of-range edge ids.
    #[inline]
    pub(crate) fn overlay_incidence_index(&self, edge: IncidenceEdgeId) -> Option<usize> {
        if edge.0 < self.bond_count {
            return None;
        }
        let idx = (edge.0 - self.bond_count) as usize;
        (idx < self.incidences.len()).then_some(idx)
    }

    /// Endpoints of `edge` in this incidence graph's node space.
    #[inline]
    pub(crate) fn edge_endpoints(&self, edge: IncidenceEdgeId) -> [IncidenceNodeId; 2] {
        self.graph
            .edge_endpoints(EdgeId::from(edge))
            .map(IncidenceNodeId::from)
    }

    /// The `(overlay_entity, atom)` pair connected by an overlay incidence edge, or `None` for a
    /// localized bond edge or out-of-range edge id.
    #[inline]
    pub(crate) fn incidence_endpoints(&self, edge: IncidenceEdgeId) -> Option<(Entity, AtomId)> {
        self.overlay_incidence_index(edge)?;
        let [first, second] = self.edge_endpoints(edge);
        let (atom_node, overlay_node) = if first.0 < self.offsets[0] {
            (first, second)
        } else {
            (second, first)
        };
        Some((self.entity(overlay_node), AtomId(atom_node.0)))
    }

    /// All overlay–atom incidence edges (excluding localized bonds).
    pub(crate) fn incidences(
        &self,
    ) -> impl ExactSizeIterator<Item = (IncidenceEdgeId, Incidence)> + '_ {
        self.incidences
            .iter()
            .enumerate()
            .map(|(offset, &incidence)| {
                (IncidenceEdgeId(self.bond_count + offset as u32), incidence)
            })
    }
}

fn build_offsets(counts: [u32; 7]) -> [u32; 7] {
    let mut offsets = [0u32; 7];
    let mut total = 0u32;
    for (slot, count) in offsets.iter_mut().zip(counts) {
        total += count;
        *slot = total;
    }
    offsets
}

impl Molecule {
    /// Build the incidence graph at the selected structural level.
    ///
    /// Localized bonds remain direct atom–atom edges (`0..bond_count`). Included overlays become
    /// pseudonodes wired to their participant atoms: stereo atoms connect to their site atom via
    /// [`Incidence::StereoAnchor`] and to actual-atom ligands via [`Incidence::StereoLigand`];
    /// stereo bonds connect to both site-bond endpoint atoms via [`Incidence::StereoAnchor`] and
    /// to actual-atom ligands via [`Incidence::StereoLigand`]. Electron contributions remain in
    /// the owning forms and do not affect this graph.
    pub(crate) fn incidence_graph(&self, level: IncidenceLevel) -> IncidenceGraph {
        let atom_count = self.raw_graph().node_count() as u32;
        let bond_count = self.bonds().count() as u32;
        let constitution = matches!(level, IncidenceLevel::Constitution | IncidenceLevel::Full);
        let stereo = matches!(level, IncidenceLevel::Full);

        let counts = [
            atom_count,
            if constitution {
                self.dative_bonds().count() as u32
            } else {
                0
            },
            if constitution {
                self.aromatic_systems().count() as u32
            } else {
                0
            },
            if constitution {
                self.multicenter_bonds().count() as u32
            } else {
                0
            },
            if constitution {
                self.noncovalent_bonds().count() as u32
            } else {
                0
            },
            if stereo {
                self.stereo_atoms().count() as u32
            } else {
                0
            },
            if stereo {
                self.stereo_bonds().count() as u32
            } else {
                0
            },
        ];
        let offsets = build_offsets(counts);

        if offsets[6] == atom_count {
            return IncidenceGraph {
                graph: self.raw_graph().clone(),
                offsets,
                bond_count,
                incidences: Vec::new(),
            };
        }

        let mut overlay_edge_count = 0usize;
        if constitution {
            overlay_edge_count += self
                .dative_bonds()
                .iter()
                .map(|bond| bond.donor_ids().len() + 1)
                .sum::<usize>()
                + self
                    .aromatic_systems()
                    .iter()
                    .map(|system| system.atom_ids().len())
                    .sum::<usize>()
                + self
                    .multicenter_bonds()
                    .iter()
                    .map(|bond| bond.atom_ids().len())
                    .sum::<usize>()
                + 2 * self.noncovalent_bonds().count();
        }
        if stereo {
            overlay_edge_count += self
                .stereo_atoms()
                .iter()
                .map(|atom| {
                    1 + atom
                        .ligands()
                        .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
                        .count()
                })
                .sum::<usize>()
                + self
                    .stereo_bonds()
                    .iter()
                    .map(|bond| {
                        2 + bond
                            .ligands()
                            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
                            .count()
                    })
                    .sum::<usize>();
        }

        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count as usize + overlay_edge_count);
        let mut incidences = Vec::with_capacity(overlay_edge_count);

        for edge in self.raw_graph().edge_ids() {
            let [first, second] = self.raw_graph().edge_endpoints(edge);
            edges.push([first.0, second.0]);
        }

        let mut node = atom_count;
        let mut push_overlay = |entity: u32, participant: u32, incidence: Incidence| {
            edges.push([entity, participant]);
            incidences.push(incidence);
        };

        if constitution {
            for v in self.dative_bonds().iter() {
                for donor in v.donor_ids() {
                    push_overlay(node, donor.0, Incidence::DativeDonor);
                }
                push_overlay(node, v.acceptor_id().0, Incidence::DativeAcceptor);
                node += 1;
            }
            for v in self.aromatic_systems().iter() {
                for atom in v.atom_ids() {
                    push_overlay(node, atom.0, Incidence::AromaticAtom);
                }
                node += 1;
            }
            for v in self.multicenter_bonds().iter() {
                for atom in v.atom_ids() {
                    push_overlay(node, atom.0, Incidence::MulticenterAtom);
                }
                node += 1;
            }
            for v in self.noncovalent_bonds().iter() {
                for a in v.atom_ids() {
                    push_overlay(node, a.0, Incidence::NoncovalentEndpoint);
                }
                node += 1;
            }
        }

        if stereo {
            for v in self.stereo_atoms().iter() {
                let site = v.site_id();
                let mut virtual_ligands = StereoVirtualLigands::empty();
                for ligand in v.ligands() {
                    if ligand.kind() != StereoLigandKind::Atom {
                        debug_assert_eq!(ligand.atom_id(), site);
                        virtual_ligands.add_kind(ligand.kind());
                    }
                }
                push_overlay(node, site.0, Incidence::StereoAnchor(virtual_ligands));
                for ligand in v.ligands() {
                    if ligand.kind() == StereoLigandKind::Atom {
                        push_overlay(node, ligand.atom_id().0, Incidence::StereoLigand);
                    }
                }
                node += 1;
            }
            for v in self.stereo_bonds().iter() {
                let site = v.site_id();
                let [first, second] = self.bond(site).atom_ids();
                let mut first_virtual = StereoVirtualLigands::empty();
                let mut second_virtual = StereoVirtualLigands::empty();
                for ligand in v.ligands() {
                    if ligand.kind() != StereoLigandKind::Atom {
                        if ligand.atom_id() == first {
                            first_virtual.add_kind(ligand.kind());
                        } else {
                            debug_assert_eq!(ligand.atom_id(), second);
                            second_virtual.add_kind(ligand.kind());
                        }
                    }
                }
                push_overlay(node, first.0, Incidence::StereoAnchor(first_virtual));
                push_overlay(node, second.0, Incidence::StereoAnchor(second_virtual));
                for ligand in v.ligands() {
                    if ligand.kind() == StereoLigandKind::Atom {
                        push_overlay(node, ligand.atom_id().0, Incidence::StereoLigand);
                    }
                }
                node += 1;
            }
        }

        debug_assert_eq!(node, offsets[6]);
        let graph = Graph::new(node as usize, &edges);
        debug_assert_eq!(graph.edge_count(), bond_count as usize + incidences.len());
        IncidenceGraph {
            graph,
            offsets,
            bond_count,
            incidences,
        }
    }
}

impl ReactionSpan {
    /// Build the incidence graph of the union frame at the selected structural level.
    ///
    /// Localized union bonds remain direct atom–atom edges (`0..bond_count`). Included overlay
    /// spans become pseudonodes wired to their participant atoms.
    pub(crate) fn incidence_graph(&self, level: IncidenceLevel) -> IncidenceGraph {
        let atom_count = self.graph().node_count() as u32;
        let bond_count = self.bonds().len() as u32;
        let constitution = matches!(level, IncidenceLevel::Constitution | IncidenceLevel::Full);
        let stereo = matches!(level, IncidenceLevel::Full);

        let counts = [
            atom_count,
            if constitution {
                self.dative_bonds().count() as u32
            } else {
                0
            },
            if constitution {
                self.aromatic_systems().count() as u32
            } else {
                0
            },
            if constitution {
                self.multicenter_bonds().count() as u32
            } else {
                0
            },
            if constitution {
                self.noncovalent_bonds().count() as u32
            } else {
                0
            },
            if stereo {
                self.stereo_atoms().count() as u32
            } else {
                0
            },
            if stereo {
                self.stereo_bonds().count() as u32
            } else {
                0
            },
        ];
        let offsets = build_offsets(counts);

        if offsets[6] == atom_count {
            return IncidenceGraph {
                graph: self.graph().clone(),
                offsets,
                bond_count,
                incidences: Vec::new(),
            };
        }

        let mut edges = Vec::new();
        let mut incidences = Vec::new();
        for edge in self.graph().edge_ids() {
            let [first, second] = self.graph().edge_endpoints(edge);
            edges.push([first.0, second.0]);
        }

        let mut node = atom_count;
        let mut push_overlay = |entity: u32, participant: u32, incidence: Incidence| {
            edges.push([entity, participant]);
            incidences.push(incidence);
        };

        if constitution {
            for id in self.dative_bonds().ids() {
                for donor in self.dative_bonds().donors(id) {
                    push_overlay(node, donor.0, Incidence::DativeDonor);
                }
                push_overlay(
                    node,
                    self.dative_bonds().acceptor(id).0,
                    Incidence::DativeAcceptor,
                );
                node += 1;
            }
            for id in self.aromatic_systems().ids() {
                for atom in self.aromatic_systems().atoms(id) {
                    push_overlay(node, atom.0, Incidence::AromaticAtom);
                }
                node += 1;
            }
            for id in self.multicenter_bonds().ids() {
                for atom in self.multicenter_bonds().atoms(id) {
                    push_overlay(node, atom.0, Incidence::MulticenterAtom);
                }
                node += 1;
            }
            for id in self.noncovalent_bonds().ids() {
                for atom in self.noncovalent_bonds().atoms(id) {
                    push_overlay(node, atom.0, Incidence::NoncovalentEndpoint);
                }
                node += 1;
            }
        }

        if stereo {
            for id in self.stereo_atoms().ids() {
                let site = self.stereo_atoms().site(id);
                let mut virtual_ligands = StereoVirtualLigands::empty();
                for ligand in self.stereo_atoms().ligands(id) {
                    if ligand.kind != StereoLigandKind::Atom {
                        debug_assert_eq!(ligand.atom_id, site);
                        virtual_ligands.add_kind(ligand.kind);
                    }
                }
                push_overlay(node, site.0, Incidence::StereoAnchor(virtual_ligands));
                for ligand in self.stereo_atoms().ligands(id) {
                    if ligand.kind == StereoLigandKind::Atom {
                        push_overlay(node, ligand.atom_id.0, Incidence::StereoLigand);
                    }
                }
                node += 1;
            }
            for id in self.stereo_bonds().ids() {
                let site = self.stereo_bonds().site(id);
                let [first, second] = self.graph().edge_endpoints(EdgeId(site.0));
                let mut first_virtual = StereoVirtualLigands::empty();
                let mut second_virtual = StereoVirtualLigands::empty();
                for ligand in self.stereo_bonds().ligands(id) {
                    if ligand.kind != StereoLigandKind::Atom {
                        if ligand.atom_id.0 == first.0 {
                            first_virtual.add_kind(ligand.kind);
                        } else {
                            debug_assert_eq!(ligand.atom_id.0, second.0);
                            second_virtual.add_kind(ligand.kind);
                        }
                    }
                }
                push_overlay(node, first.0, Incidence::StereoAnchor(first_virtual));
                push_overlay(node, second.0, Incidence::StereoAnchor(second_virtual));
                for ligand in self.stereo_bonds().ligands(id) {
                    if ligand.kind == StereoLigandKind::Atom {
                        push_overlay(node, ligand.atom_id.0, Incidence::StereoLigand);
                    }
                }
                node += 1;
            }
        }

        debug_assert_eq!(node, offsets[6]);
        let graph = Graph::new(node as usize, &edges);
        debug_assert_eq!(graph.edge_count(), bond_count as usize + incidences.len());
        IncidenceGraph {
            graph,
            offsets,
            bond_count,
            incidences,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::MoleculeEntries;
    use crate::ir::multicenter::MulticenterBondForm;
    use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use crate::ir::reaction_span::ReactionSpanEntries;
    use crate::ir::span::EntitySpan;
    use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoCoset, StereoKind};

    // Six carbons; chain bonds 0-1-2-3 (BondId 0,1,2); a dative 0→3; an aromatic
    // system {0,1,2}; a multicenter {3,4,5}; a noncovalent 0···5; a stereo atom on
    // site 1; a stereo bond on site BondId(1).
    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(3), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 0, 2]),
            )],
            multicenter: vec![(
                vec![AtomId(3), AtomId(4), AtomId(5)],
                MulticenterBondForm::from_electrons(vec![2, 0, 1]),
            )],
            noncovalent: vec![(
                [AtomId(0), AtomId(5)],
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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::topology(
        IncidenceLevel::Topology,
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
        ],
    )]
    #[case::constitution(
        IncidenceLevel::Constitution,
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::DativeBond(DativeBondId(0)), Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(0)),
        ],
    )]
    #[case::full(
        IncidenceLevel::Full,
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::DativeBond(DativeBondId(0)), Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(0)),
            Entity::StereoAtom(StereoAtomId(0)), Entity::StereoBond(StereoBondId(0)),
        ],
    )]
    fn test_molecule_incidence_graph(
        molecule: Molecule,
        #[case] level: IncidenceLevel,
        #[case] expected: Vec<Entity>,
    ) {
        let inc = molecule.incidence_graph(level);
        assert_eq!(inc.node_count(), expected.len());
        assert_eq!(inc.raw_graph().node_count(), expected.len());
        assert_eq!(inc.bond_count(), 3);
        let got: Vec<Entity> = inc.node_ids().map(|node| inc.entity(node)).collect();
        assert_eq!(got, expected);
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(inc.node_of(e), Some(IncidenceNodeId(i as u32)));
        }
        assert_eq!(inc.node_of(Entity::Bond(BondId(0))), None);
        assert_eq!(inc.bond_edge(BondId(0)), Some(IncidenceEdgeId(0)));
        assert_eq!(inc.bond_edge(BondId(2)), Some(IncidenceEdgeId(2)));
        assert_eq!(inc.bond_edge(BondId(3)), None);
        if level == IncidenceLevel::Topology {
            assert_eq!(inc.incidence_count(), 0);
            assert_eq!(inc.raw_graph(), molecule.raw_graph());
            assert_eq!(inc.node_of(Entity::DativeBond(DativeBondId(0))), None);
        }
        if level != IncidenceLevel::Full {
            assert_eq!(inc.node_of(Entity::StereoAtom(StereoAtomId(0))), None);
            assert_eq!(inc.node_of(Entity::StereoBond(StereoBondId(0))), None);
        }
    }

    #[rstest]
    #[case::bond(IncidenceEdgeId(0), IncidenceEdge::Bond(BondId(0)))]
    #[case::dative_donor(IncidenceEdgeId(3), IncidenceEdge::Incidence(Incidence::DativeDonor))]
    #[case::dative_acceptor(
        IncidenceEdgeId(4),
        IncidenceEdge::Incidence(Incidence::DativeAcceptor)
    )]
    #[case::aromatic(IncidenceEdgeId(5), IncidenceEdge::Incidence(Incidence::AromaticAtom))]
    #[case::multicenter(
        IncidenceEdgeId(8),
        IncidenceEdge::Incidence(Incidence::MulticenterAtom)
    )]
    #[case::noncovalent(
        IncidenceEdgeId(11),
        IncidenceEdge::Incidence(Incidence::NoncovalentEndpoint)
    )]
    #[case::stereo_atom_anchor(
        IncidenceEdgeId(13),
        IncidenceEdge::Incidence(Incidence::StereoAnchor(
            StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
        ))
    )]
    #[case::stereo_ligand(IncidenceEdgeId(14), IncidenceEdge::Incidence(Incidence::StereoLigand))]
    #[case::stereo_bond_anchor(
        IncidenceEdgeId(16),
        IncidenceEdge::Incidence(Incidence::StereoAnchor(
            StereoVirtualLigands::IMPLICIT_HYDROGEN
        ))
    )]
    fn test_incidence_graph_edge(
        molecule: Molecule,
        #[case] edge: IncidenceEdgeId,
        #[case] expected: IncidenceEdge,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(incidence.edge(edge), expected);
    }

    #[rstest]
    #[case::aromatic(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::default())],
            ..Default::default()
        }),
        Incidence::AromaticAtom,
    )]
    #[case::multicenter(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            multicenter: vec![(vec![AtomId(0)], MulticenterBondForm::default())],
            ..Default::default()
        }),
        Incidence::MulticenterAtom,
    )]
    #[case::malformed_aromatic_electron_length(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1, 1])),
            )],
            ..Default::default()
        }),
        Incidence::AromaticAtom,
    )]
    fn test_incidence_graph_incidence_electrons(
        #[case] molecule: Molecule,
        #[case] expected: Incidence,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Constitution);
        assert_eq!(
            incidence.edge(IncidenceEdgeId(0)),
            IncidenceEdge::Incidence(expected)
        );
    }

    #[rstest]
    fn test_incidence_graph_edges(molecule: Molecule) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(incidence.edge_count(), 20);
        assert_eq!(incidence.bond_count(), 3);
        assert_eq!(incidence.incidences().len(), 17);
        let edges: Vec<_> = (0..incidence.edge_count() as u32)
            .map(IncidenceEdgeId)
            .map(|edge| (edge, incidence.edge(edge)))
            .collect();
        assert_eq!(
            edges,
            vec![
                (IncidenceEdgeId(0), IncidenceEdge::Bond(BondId(0))),
                (IncidenceEdgeId(1), IncidenceEdge::Bond(BondId(1))),
                (IncidenceEdgeId(2), IncidenceEdge::Bond(BondId(2))),
                (
                    IncidenceEdgeId(3),
                    IncidenceEdge::Incidence(Incidence::DativeDonor)
                ),
                (
                    IncidenceEdgeId(4),
                    IncidenceEdge::Incidence(Incidence::DativeAcceptor)
                ),
                (
                    IncidenceEdgeId(5),
                    IncidenceEdge::Incidence(Incidence::AromaticAtom)
                ),
                (
                    IncidenceEdgeId(6),
                    IncidenceEdge::Incidence(Incidence::AromaticAtom)
                ),
                (
                    IncidenceEdgeId(7),
                    IncidenceEdge::Incidence(Incidence::AromaticAtom)
                ),
                (
                    IncidenceEdgeId(8),
                    IncidenceEdge::Incidence(Incidence::MulticenterAtom)
                ),
                (
                    IncidenceEdgeId(9),
                    IncidenceEdge::Incidence(Incidence::MulticenterAtom)
                ),
                (
                    IncidenceEdgeId(10),
                    IncidenceEdge::Incidence(Incidence::MulticenterAtom)
                ),
                (
                    IncidenceEdgeId(11),
                    IncidenceEdge::Incidence(Incidence::NoncovalentEndpoint)
                ),
                (
                    IncidenceEdgeId(12),
                    IncidenceEdge::Incidence(Incidence::NoncovalentEndpoint)
                ),
                (
                    IncidenceEdgeId(13),
                    IncidenceEdge::Incidence(Incidence::StereoAnchor(
                        StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
                    )),
                ),
                (
                    IncidenceEdgeId(14),
                    IncidenceEdge::Incidence(Incidence::StereoLigand)
                ),
                (
                    IncidenceEdgeId(15),
                    IncidenceEdge::Incidence(Incidence::StereoLigand)
                ),
                (
                    IncidenceEdgeId(16),
                    IncidenceEdge::Incidence(Incidence::StereoAnchor(
                        StereoVirtualLigands::IMPLICIT_HYDROGEN,
                    )),
                ),
                (
                    IncidenceEdgeId(17),
                    IncidenceEdge::Incidence(Incidence::StereoAnchor(
                        StereoVirtualLigands::LONE_PAIR,
                    )),
                ),
                (
                    IncidenceEdgeId(18),
                    IncidenceEdge::Incidence(Incidence::StereoLigand)
                ),
                (
                    IncidenceEdgeId(19),
                    IncidenceEdge::Incidence(Incidence::StereoLigand)
                ),
            ],
        );
        assert_eq!(incidence.incidence_endpoints(IncidenceEdgeId(0)), None);
        assert_eq!(
            incidence.incidence_endpoints(IncidenceEdgeId(3)),
            Some((Entity::DativeBond(DativeBondId(0)), AtomId(0)))
        );
        assert_eq!(
            incidence.incidence_endpoints(IncidenceEdgeId(4)),
            Some((Entity::DativeBond(DativeBondId(0)), AtomId(3)))
        );
    }

    #[rstest]
    // Atoms connect directly through localized bonds plus incident overlays.
    #[case::dative(6, vec![0, 3])]
    #[case::aromatic(7, vec![0, 1, 2])]
    #[case::multicenter(8, vec![3, 4, 5])]
    #[case::noncovalent(9, vec![0, 5])]
    // Stereo nodes attach once to each anchor atom (folding virtual ligands into the
    // anchor edge) and once to each actual-atom ligand.
    #[case::stereo_atom(10, vec![0, 1, 2])]
    #[case::stereo_bond(11, vec![0, 1, 2, 3])]
    fn test_molecule_incidence_graph_neighbors(
        molecule: Molecule,
        #[case] node: u32,
        #[case] expected: Vec<u32>,
    ) {
        let inc = molecule.incidence_graph(IncidenceLevel::Full);
        let mut got: Vec<u32> = (0..inc.edge_count() as u32)
            .map(|index| inc.edge_endpoints(IncidenceEdgeId(index)))
            .filter_map(|[a, b]| {
                let (a, b) = (a.0, b.0);
                if a == node {
                    Some(b)
                } else if b == node {
                    Some(a)
                } else {
                    None
                }
            })
            .collect();
        got.sort_unstable();
        assert_eq!(got, expected);
    }

    #[rstest]
    #[case::stereo_atom_virtual_ligands(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            ([IncidenceNodeId(0), IncidenceNodeId(1)], IncidenceEdge::Bond(BondId(0))),
            ([IncidenceNodeId(0), IncidenceNodeId(2)], IncidenceEdge::Bond(BondId(1))),
            (
                [IncidenceNodeId(0), IncidenceNodeId(3)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(
                    StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
                )),
            ),
            ([IncidenceNodeId(1), IncidenceNodeId(3)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            ([IncidenceNodeId(2), IncidenceNodeId(3)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
        ],
    )]
    #[case::stereo_bond_endpoint_anchors(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
                (AtomId(1), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            ([IncidenceNodeId(0), IncidenceNodeId(1)], IncidenceEdge::Bond(BondId(0))),
            ([IncidenceNodeId(0), IncidenceNodeId(2)], IncidenceEdge::Bond(BondId(1))),
            ([IncidenceNodeId(1), IncidenceNodeId(3)], IncidenceEdge::Bond(BondId(2))),
            ([IncidenceNodeId(1), IncidenceNodeId(4)], IncidenceEdge::Bond(BondId(3))),
            (
                [IncidenceNodeId(0), IncidenceNodeId(5)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(StereoVirtualLigands::VIRTUAL)),
            ),
            (
                [IncidenceNodeId(1), IncidenceNodeId(5)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(StereoVirtualLigands::empty())),
            ),
            ([IncidenceNodeId(2), IncidenceNodeId(5)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            ([IncidenceNodeId(3), IncidenceNodeId(5)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            ([IncidenceNodeId(4), IncidenceNodeId(5)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
        ],
    )]
    fn test_molecule_incidence_graph_stereo(
        #[case] molecule: Molecule,
        #[case] expected: Vec<([IncidenceNodeId; 2], IncidenceEdge)>,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(
            (0..incidence.edge_count() as u32)
                .map(IncidenceEdgeId)
                .map(|edge| (incidence.edge_endpoints(edge), incidence.edge(edge)))
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::stereo_atom_and_bond(
        ReactionSpan::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::F)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::Cl)),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(2))),
                (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
                (AtomId(1), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                ],
                EntitySpan::Unchanged(StereoAtomForm::new(
                    StereoKind::Tetrahedral,
                    StereoCoset::Lit(0),
                )),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                EntitySpan::Unchanged(StereoBondForm::new(
                    StereoKind::CisTrans,
                    StereoCoset::Lit(1),
                )),
            )],
            ..Default::default()
        }),
        vec![
            ([IncidenceNodeId(0), IncidenceNodeId(1)], IncidenceEdge::Bond(BondId(0))),
            ([IncidenceNodeId(0), IncidenceNodeId(2)], IncidenceEdge::Bond(BondId(1))),
            ([IncidenceNodeId(1), IncidenceNodeId(3)], IncidenceEdge::Bond(BondId(2))),
            (
                [IncidenceNodeId(0), IncidenceNodeId(4)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(
                    StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::VIRTUAL,
                )),
            ),
            ([IncidenceNodeId(1), IncidenceNodeId(4)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            ([IncidenceNodeId(2), IncidenceNodeId(4)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            (
                [IncidenceNodeId(0), IncidenceNodeId(5)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(
                    StereoVirtualLigands::IMPLICIT_HYDROGEN,
                )),
            ),
            (
                [IncidenceNodeId(1), IncidenceNodeId(5)],
                IncidenceEdge::Incidence(Incidence::StereoAnchor(
                    StereoVirtualLigands::LONE_PAIR,
                )),
            ),
            ([IncidenceNodeId(2), IncidenceNodeId(5)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
            ([IncidenceNodeId(3), IncidenceNodeId(5)], IncidenceEdge::Incidence(Incidence::StereoLigand)),
        ],
    )]
    fn test_reaction_span_incidence_graph(
        #[case] span: ReactionSpan,
        #[case] expected: Vec<([IncidenceNodeId; 2], IncidenceEdge)>,
    ) {
        let incidence = span.incidence_graph(IncidenceLevel::Full);
        assert_eq!(
            (0..incidence.edge_count() as u32)
                .map(IncidenceEdgeId)
                .map(|edge| (incidence.edge_endpoints(edge), incidence.edge(edge)))
                .collect::<Vec<_>>(),
            expected,
        );
    }
}
