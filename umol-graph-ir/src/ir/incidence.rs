//! Incidence (Levi) graph: relations lifted to pseudonodes for symmetry analysis.

use bitflags::bitflags;
use strum::EnumCount;
use umol_graph_core::{EdgeId, Graph, NodeId};

use super::entity::{Entity, EntityKind};
use super::ligand::StereoLigandKind;
use super::molecule::Molecule;
use super::reaction_span::ReactionSpan;

/// Structural level represented by an [`IncidenceGraph`].
///
/// The variants form a nested hierarchy. [`Topology`](Self::Topology) contains
/// atoms and localized bonds. [`Constitution`](Self::Constitution) additionally
/// contains dative bonds, aromatic systems, multicenter bonds, and noncovalent
/// bonds. [`Full`](Self::Full) additionally contains stereo atoms and stereo
/// bonds. Constraints do not contribute to any structural level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IncidenceLevel {
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
    pub struct StereoVirtualLigands: u8 {
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

/// The chemical role of one edge in an [`IncidenceGraph`].
///
/// Every incidence edge is simple: each pair of entity nodes is connected by at most one edge.
/// Stereo entities encode virtual ligands anchored at their site atom or bond endpoints on the
/// corresponding [`StereoAnchor`](Self::StereoAnchor) edge, while actual-atom ligands use
/// [`StereoLigand`](Self::StereoLigand). Electron contributions and reaction-side attributes
/// remain in the owning entity forms; they are not incidence labels.
///
/// # Semantic properties
///
/// The total order follows the declaration order of the variants. Canonicalization
/// augments these roles with the owning forms' electron contributions when
/// constructing colors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Incidence {
    BondEndpoint,
    DativeDonor,
    DativeAcceptor,
    AromaticAtom,
    MulticenterAtom,
    NoncovalentEndpoint,
    StereoSite,
    StereoAnchor(StereoVirtualLigands),
    StereoLigand,
}

/// A molecule's incidence graph: one node per selected molecule entity, with
/// typed edges for the entity's participant occurrences. `entity(node)` recovers
/// the entity behind a node; atoms occupy `0..atom_count` (node index = atom
/// index), so localized bond `BondId(k)` is the node at `atom_count + k`.
#[derive(Clone, Debug)]
pub struct IncidenceGraph {
    graph: Graph,
    // Per-kind block sizes, indexed by `EntityKind as usize` (node-layout order).
    // An entity's id is its index within its block, so node↔entity is offset
    // arithmetic — no per-node table.
    entity_counts: [u32; EntityKind::COUNT],
    // Indexed by the corresponding `Graph` edge id.
    incidences: Vec<Incidence>,
}

impl IncidenceGraph {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Number of nodes of the given kind.
    pub fn entity_count(&self, kind: EntityKind) -> usize {
        self.entity_counts[kind as usize] as usize
    }

    /// The molecule entity a node represents.
    pub fn entity(&self, node: NodeId) -> Entity {
        let n = node.index();
        let mut offset = 0usize;
        for (block, &count) in self.entity_counts.iter().enumerate() {
            let end = offset + count as usize;
            if n < end {
                let kind = EntityKind::try_from(block as u8).expect("block is a valid kind index");
                return kind.with_id((n - offset) as u32);
            }
            offset = end;
        }
        panic!("incidence node {n} out of range");
    }

    /// The node representing `entity` — inverse of [`entity`](Self::entity).
    pub fn node_of(&self, entity: Entity) -> NodeId {
        let block = entity.kind() as usize;
        let offset: u32 = self.entity_counts[..block].iter().sum();
        NodeId(offset + entity.kind_id() as u32)
    }

    /// The typed participant occurrence represented by `edge`.
    ///
    /// # Panics
    ///
    /// Panics when `edge` is not in this incidence graph.
    pub fn incidence(&self, edge: EdgeId) -> &Incidence {
        &self.incidences[edge.index()]
    }

    /// All incidence edges and their typed participant occurrences.
    pub fn incidences(&self) -> impl ExactSizeIterator<Item = (EdgeId, &Incidence)> {
        self.graph.edge_ids().zip(self.incidences.iter())
    }
}

impl Molecule {
    /// Build the incidence (Levi) graph at the selected structural level. Localized
    /// bonds and included overlays become entity nodes wired to every participant
    /// occurrence. Stereo atoms connect to their site atom via [`Incidence::StereoAnchor`]
    /// (carrying anchored virtual ligands) and to actual-atom ligands via
    /// [`Incidence::StereoLigand`]. Stereo bonds connect to their site bond via
    /// [`Incidence::StereoSite`], to both bond endpoint atoms via [`Incidence::StereoAnchor`],
    /// and to actual-atom ligands via [`Incidence::StereoLigand`].
    /// Electron contributions remain in the owning forms and do not affect this graph.
    pub fn incidence_graph(&self, level: IncidenceLevel) -> IncidenceGraph {
        let atom_count = self.raw_graph().node_count();
        let constitution = matches!(level, IncidenceLevel::Constitution | IncidenceLevel::Full);
        let stereo = matches!(level, IncidenceLevel::Full);

        let entity_counts = [
            atom_count as u32,
            self.bonds().ids().count() as u32,
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

        // Entity nodes follow the atom block in the same fixed order as `entity_counts`,
        // each wired to its participant occurrences. Bonds come first, fixing
        // BondId(k) at atom_count + k — relied on by the stereo-bond site link.
        let mut incidence_count = 2 * self.bonds().count();
        if constitution {
            incidence_count += self
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
            incidence_count += self
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
                        3 + bond
                            .ligands()
                            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
                            .count()
                    })
                    .sum::<usize>();
        }
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(incidence_count);
        let mut incidences = Vec::with_capacity(incidence_count);
        let mut node = atom_count as u32;

        let mut push = |entity: u32, participant: u32, incidence: Incidence| {
            edges.push([entity, participant]);
            incidences.push(incidence);
        };

        for id in self.bonds().ids() {
            let [a, b] = self.bond(id).atom_ids();
            push(node, a.index() as u32, Incidence::BondEndpoint);
            push(node, b.index() as u32, Incidence::BondEndpoint);
            node += 1;
        }

        if constitution {
            for v in self.dative_bonds().iter() {
                for donor in v.donor_ids() {
                    push(node, donor.index() as u32, Incidence::DativeDonor);
                }
                push(
                    node,
                    v.acceptor_id().index() as u32,
                    Incidence::DativeAcceptor,
                );
                node += 1;
            }
            for v in self.aromatic_systems().iter() {
                for atom in v.atom_ids() {
                    push(node, atom.index() as u32, Incidence::AromaticAtom);
                }
                node += 1;
            }
            for v in self.multicenter_bonds().iter() {
                for atom in v.atom_ids() {
                    push(node, atom.index() as u32, Incidence::MulticenterAtom);
                }
                node += 1;
            }
            for v in self.noncovalent_bonds().iter() {
                for a in v.atom_ids() {
                    push(node, a.index() as u32, Incidence::NoncovalentEndpoint);
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
                push(
                    node,
                    site.index() as u32,
                    Incidence::StereoAnchor(virtual_ligands),
                );
                for ligand in v.ligands() {
                    if ligand.kind() == StereoLigandKind::Atom {
                        push(
                            node,
                            ligand.atom_id().index() as u32,
                            Incidence::StereoLigand,
                        );
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
                push(
                    node,
                    atom_count as u32 + site.index() as u32,
                    Incidence::StereoSite,
                );
                push(
                    node,
                    first.index() as u32,
                    Incidence::StereoAnchor(first_virtual),
                );
                push(
                    node,
                    second.index() as u32,
                    Incidence::StereoAnchor(second_virtual),
                );
                for ligand in v.ligands() {
                    if ligand.kind() == StereoLigandKind::Atom {
                        push(
                            node,
                            ligand.atom_id().index() as u32,
                            Incidence::StereoLigand,
                        );
                    }
                }
                node += 1;
            }
        }

        let graph = Graph::new(node as usize, &edges);
        debug_assert_eq!(graph.edge_count(), incidences.len());
        IncidenceGraph {
            graph,
            entity_counts,
            incidences,
        }
    }
}

impl ReactionSpan {
    /// Build the incidence graph of the union frame at the selected structural level.
    /// Edges identify bond endpoints, dative donors and acceptors, overlay atoms,
    /// stereo sites, stereo anchors (carrying virtual ligands), and actual-atom stereo ligands.
    pub fn incidence_graph(&self, level: IncidenceLevel) -> IncidenceGraph {
        let atom_count = self.graph().node_count();
        let constitution = matches!(level, IncidenceLevel::Constitution | IncidenceLevel::Full);
        let stereo = matches!(level, IncidenceLevel::Full);
        let entity_counts = [
            atom_count as u32,
            self.bonds().len() as u32,
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

        let mut edges = Vec::new();
        let mut incidences = Vec::new();
        let mut node = atom_count as u32;
        let mut push = |entity: u32, participant: u32, incidence: Incidence| {
            edges.push([entity, participant]);
            incidences.push(incidence);
        };

        for edge in self.graph().edge_ids() {
            let [first, second] = self.graph().edge_endpoints(edge);
            push(node, first.0, Incidence::BondEndpoint);
            push(node, second.0, Incidence::BondEndpoint);
            node += 1;
        }

        if constitution {
            for id in self.dative_bonds().ids() {
                for donor in self.dative_bonds().donors(id) {
                    push(node, donor.0, Incidence::DativeDonor);
                }
                push(
                    node,
                    self.dative_bonds().acceptor(id).0,
                    Incidence::DativeAcceptor,
                );
                node += 1;
            }
            for id in self.aromatic_systems().ids() {
                for atom in self.aromatic_systems().atoms(id) {
                    push(node, atom.0, Incidence::AromaticAtom);
                }
                node += 1;
            }
            for id in self.multicenter_bonds().ids() {
                for atom in self.multicenter_bonds().atoms(id) {
                    push(node, atom.0, Incidence::MulticenterAtom);
                }
                node += 1;
            }
            for id in self.noncovalent_bonds().ids() {
                for atom in self.noncovalent_bonds().atoms(id) {
                    push(node, atom.0, Incidence::NoncovalentEndpoint);
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
                push(node, site.0, Incidence::StereoAnchor(virtual_ligands));
                for ligand in self.stereo_atoms().ligands(id) {
                    if ligand.kind == StereoLigandKind::Atom {
                        push(node, ligand.atom_id.0, Incidence::StereoLigand);
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
                push(node, atom_count as u32 + site.0, Incidence::StereoSite);
                push(node, first.0, Incidence::StereoAnchor(first_virtual));
                push(node, second.0, Incidence::StereoAnchor(second_virtual));
                for ligand in self.stereo_bonds().ligands(id) {
                    if ligand.kind == StereoLigandKind::Atom {
                        push(node, ligand.atom_id.0, Incidence::StereoLigand);
                    }
                }
                node += 1;
            }
        }

        let graph = Graph::new(node as usize, &edges);
        IncidenceGraph {
            graph,
            entity_counts,
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
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
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
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
        ],
    )]
    #[case::constitution(
        IncidenceLevel::Constitution,
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
            Entity::DativeBond(DativeBondId(0)), Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(0)),
        ],
    )]
    #[case::full(
        IncidenceLevel::Full,
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
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
        assert_eq!(inc.graph().node_count(), expected.len());
        let got: Vec<Entity> = (0..expected.len())
            .map(|i| inc.entity(NodeId(i as u32)))
            .collect();
        assert_eq!(got, expected);
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(inc.node_of(e), NodeId(i as u32));
        }
    }

    #[rstest]
    #[case::bond(EdgeId(0), Incidence::BondEndpoint)]
    #[case::dative_donor(EdgeId(6), Incidence::DativeDonor)]
    #[case::dative_acceptor(EdgeId(7), Incidence::DativeAcceptor)]
    #[case::aromatic(EdgeId(8), Incidence::AromaticAtom)]
    #[case::multicenter(EdgeId(11), Incidence::MulticenterAtom)]
    #[case::noncovalent(EdgeId(14), Incidence::NoncovalentEndpoint)]
    #[case::stereo_atom_anchor(
        EdgeId(16),
        Incidence::StereoAnchor(
            StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR
        )
    )]
    #[case::stereo_ligand(EdgeId(17), Incidence::StereoLigand)]
    #[case::stereo_bond_site(EdgeId(19), Incidence::StereoSite)]
    #[case::stereo_bond_anchor(
        EdgeId(20),
        Incidence::StereoAnchor(StereoVirtualLigands::IMPLICIT_HYDROGEN)
    )]
    fn test_incidence_graph_incidence(
        molecule: Molecule,
        #[case] edge: EdgeId,
        #[case] expected: Incidence,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(incidence.incidence(edge), &expected);
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
    fn test_incidence_graph_incidence_electrons(
        #[case] molecule: Molecule,
        #[case] expected: Incidence,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Constitution);
        assert_eq!(incidence.incidence(EdgeId(0)), &expected);
    }

    #[rstest]
    fn test_incidence_graph_incidences(molecule: Molecule) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        let iterator = incidence.incidences();
        assert_eq!(iterator.len(), 24);
        assert_eq!(
            iterator
                .map(|(edge, &value)| (edge, value))
                .collect::<Vec<_>>(),
            vec![
                (EdgeId(0), Incidence::BondEndpoint),
                (EdgeId(1), Incidence::BondEndpoint),
                (EdgeId(2), Incidence::BondEndpoint),
                (EdgeId(3), Incidence::BondEndpoint),
                (EdgeId(4), Incidence::BondEndpoint),
                (EdgeId(5), Incidence::BondEndpoint),
                (EdgeId(6), Incidence::DativeDonor),
                (EdgeId(7), Incidence::DativeAcceptor),
                (EdgeId(8), Incidence::AromaticAtom),
                (EdgeId(9), Incidence::AromaticAtom),
                (EdgeId(10), Incidence::AromaticAtom),
                (EdgeId(11), Incidence::MulticenterAtom),
                (EdgeId(12), Incidence::MulticenterAtom),
                (EdgeId(13), Incidence::MulticenterAtom),
                (EdgeId(14), Incidence::NoncovalentEndpoint),
                (EdgeId(15), Incidence::NoncovalentEndpoint),
                (
                    EdgeId(16),
                    Incidence::StereoAnchor(
                        StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
                    ),
                ),
                (EdgeId(17), Incidence::StereoLigand),
                (EdgeId(18), Incidence::StereoLigand),
                (EdgeId(19), Incidence::StereoSite),
                (
                    EdgeId(20),
                    Incidence::StereoAnchor(StereoVirtualLigands::IMPLICIT_HYDROGEN),
                ),
                (
                    EdgeId(21),
                    Incidence::StereoAnchor(StereoVirtualLigands::LONE_PAIR),
                ),
                (EdgeId(22), Incidence::StereoLigand),
                (EdgeId(23), Incidence::StereoLigand),
            ],
        );
    }

    #[rstest]
    // Localized bonds are pseudonodes wired to both endpoints (not atom-atom edges).
    #[case::bond(6, vec![0, 1])]
    // Overlays wire to all participant atoms.
    #[case::dative(9, vec![0, 3])]
    #[case::aromatic(10, vec![0, 1, 2])]
    #[case::multicenter(11, vec![3, 4, 5])]
    #[case::noncovalent(12, vec![0, 5])]
    // Stereo nodes attach once to each anchor atom (folding virtual ligands into the
    // anchor edge) and once to each actual-atom ligand; stereo bonds also attach to
    // their site bond node.
    #[case::stereo_atom(13, vec![0, 1, 2])]
    #[case::stereo_bond(14, vec![0, 1, 2, 3, 7])]
    fn test_molecule_incidence_graph_neighbors(
        molecule: Molecule,
        #[case] node: u32,
        #[case] expected: Vec<u32>,
    ) {
        let inc = molecule.incidence_graph(IncidenceLevel::Full);
        let graph = inc.graph();
        let mut got: Vec<u32> = graph
            .edge_ids()
            .map(|e| graph.edge_endpoints(e))
            .filter_map(|[a, b]| {
                let (a, b) = (a.index() as u32, b.index() as u32);
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
            ([NodeId(0), NodeId(3)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(3)], Incidence::BondEndpoint),
            ([NodeId(0), NodeId(4)], Incidence::BondEndpoint),
            ([NodeId(2), NodeId(4)], Incidence::BondEndpoint),
            (
                [NodeId(0), NodeId(5)],
                Incidence::StereoAnchor(
                    StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
                ),
            ),
            ([NodeId(1), NodeId(5)], Incidence::StereoLigand),
            ([NodeId(2), NodeId(5)], Incidence::StereoLigand),
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
            ([NodeId(0), NodeId(5)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(5)], Incidence::BondEndpoint),
            ([NodeId(0), NodeId(6)], Incidence::BondEndpoint),
            ([NodeId(2), NodeId(6)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(7)], Incidence::BondEndpoint),
            ([NodeId(3), NodeId(7)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(8)], Incidence::BondEndpoint),
            ([NodeId(4), NodeId(8)], Incidence::BondEndpoint),
            ([NodeId(5), NodeId(9)], Incidence::StereoSite),
            (
                [NodeId(0), NodeId(9)],
                Incidence::StereoAnchor(StereoVirtualLigands::VIRTUAL),
            ),
            (
                [NodeId(1), NodeId(9)],
                Incidence::StereoAnchor(StereoVirtualLigands::empty()),
            ),
            ([NodeId(2), NodeId(9)], Incidence::StereoLigand),
            ([NodeId(3), NodeId(9)], Incidence::StereoLigand),
            ([NodeId(4), NodeId(9)], Incidence::StereoLigand),
        ],
    )]
    fn test_molecule_incidence_graph_stereo(
        #[case] molecule: Molecule,
        #[case] expected: Vec<([NodeId; 2], Incidence)>,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(
            incidence
                .incidences()
                .map(|(edge, &value)| (incidence.graph().edge_endpoints(edge), value))
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
            ([NodeId(0), NodeId(4)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(4)], Incidence::BondEndpoint),
            ([NodeId(0), NodeId(5)], Incidence::BondEndpoint),
            ([NodeId(2), NodeId(5)], Incidence::BondEndpoint),
            ([NodeId(1), NodeId(6)], Incidence::BondEndpoint),
            ([NodeId(3), NodeId(6)], Incidence::BondEndpoint),
            (
                [NodeId(0), NodeId(7)],
                Incidence::StereoAnchor(
                    StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::VIRTUAL,
                ),
            ),
            ([NodeId(1), NodeId(7)], Incidence::StereoLigand),
            ([NodeId(2), NodeId(7)], Incidence::StereoLigand),
            ([NodeId(4), NodeId(8)], Incidence::StereoSite),
            (
                [NodeId(0), NodeId(8)],
                Incidence::StereoAnchor(StereoVirtualLigands::IMPLICIT_HYDROGEN),
            ),
            (
                [NodeId(1), NodeId(8)],
                Incidence::StereoAnchor(StereoVirtualLigands::LONE_PAIR),
            ),
            ([NodeId(2), NodeId(8)], Incidence::StereoLigand),
            ([NodeId(3), NodeId(8)], Incidence::StereoLigand),
        ],
    )]
    fn test_reaction_span_incidence_graph(
        #[case] span: ReactionSpan,
        #[case] expected: Vec<([NodeId; 2], Incidence)>,
    ) {
        let incidence = span.incidence_graph(IncidenceLevel::Full);
        assert_eq!(
            incidence
                .incidences()
                .map(|(edge, &value)| (incidence.graph().edge_endpoints(edge), value))
                .collect::<Vec<_>>(),
            expected,
        );
    }
}
