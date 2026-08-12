//! Incidence (Levi) graph: relations lifted to pseudonodes for symmetry analysis.

use std::cmp::Ordering;

use strum::EnumCount;
use umol_graph_core::{EdgeId, Graph, NodeId};

use super::delta::EntitySpan;
use super::electrons::ElectronCountsForm;
use super::entity::{Entity, EntityKind};
use super::ligand::StereoLigandKind;
use super::molecule::Molecule;
use super::num::NumForm;
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

/// The complete meaning of one edge in an [`IncidenceGraph`].
///
/// Each value describes one participant occurrence. Parallel occurrences remain
/// separate graph edges and therefore have separate `Incidence` values.
///
/// # Semantic properties
///
/// The total order uses the frozen aggregate-canonicalization schema rather
/// than enum declaration order. For normalized values, it agrees with the
/// typed incidence keys used to form canonicalization classes.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Incidence {
    BondEndpoint,
    DativeDonor,
    DativeAcceptor,
    AromaticParticipant(NumForm),
    AromaticParticipantSpan(EntitySpan<NumForm>),
    MulticenterParticipant(NumForm),
    MulticenterParticipantSpan(EntitySpan<NumForm>),
    NoncovalentEndpoint,
    StereoSite,
    StereoLigand(StereoLigandKind),
}

impl Incidence {
    fn position(&self) -> u8 {
        match self {
            Self::BondEndpoint => 0,
            Self::DativeDonor => 1,
            Self::DativeAcceptor => 2,
            Self::AromaticParticipant(_) | Self::AromaticParticipantSpan(_) => 3,
            Self::MulticenterParticipant(_) | Self::MulticenterParticipantSpan(_) => 4,
            Self::NoncovalentEndpoint => 5,
            Self::StereoSite => 6,
            Self::StereoLigand(_) => 7,
        }
    }
}

impl Ord for Incidence {
    fn cmp(&self, other: &Self) -> Ordering {
        self.position()
            .cmp(&other.position())
            .then_with(|| match (self, other) {
                (Self::AromaticParticipant(lhs), Self::AromaticParticipant(rhs))
                | (Self::MulticenterParticipant(lhs), Self::MulticenterParticipant(rhs)) => {
                    lhs.cmp(rhs)
                }
                (Self::AromaticParticipantSpan(lhs), Self::AromaticParticipantSpan(rhs))
                | (Self::MulticenterParticipantSpan(lhs), Self::MulticenterParticipantSpan(rhs)) => {
                    entity_span_cmp(lhs, rhs)
                }
                (Self::AromaticParticipant(_), Self::AromaticParticipantSpan(_))
                | (Self::MulticenterParticipant(_), Self::MulticenterParticipantSpan(_)) => {
                    Ordering::Less
                }
                (Self::AromaticParticipantSpan(_), Self::AromaticParticipant(_))
                | (Self::MulticenterParticipantSpan(_), Self::MulticenterParticipant(_)) => {
                    Ordering::Greater
                }
                (Self::StereoLigand(lhs), Self::StereoLigand(rhs)) => lhs.cmp(rhs),
                _ => Ordering::Equal,
            })
    }
}

fn entity_span_cmp<T: Ord>(left: &EntitySpan<T>, right: &EntitySpan<T>) -> Ordering {
    let position = |span: &EntitySpan<T>| match span {
        EntitySpan::Unchanged(_) => 0,
        EntitySpan::Added(_) => 1,
        EntitySpan::Removed(_) => 2,
        EntitySpan::Modified { .. } => 3,
    };
    position(left)
        .cmp(&position(right))
        .then_with(|| match (left, right) {
            (EntitySpan::Unchanged(left), EntitySpan::Unchanged(right))
            | (EntitySpan::Added(left), EntitySpan::Added(right))
            | (EntitySpan::Removed(left), EntitySpan::Removed(right)) => left.cmp(right),
            (
                EntitySpan::Modified {
                    lhs: left_lhs,
                    rhs: left_rhs,
                },
                EntitySpan::Modified {
                    lhs: right_lhs,
                    rhs: right_rhs,
                },
            ) => left_lhs
                .cmp(right_lhs)
                .then_with(|| left_rhs.cmp(right_rhs)),
            _ => Ordering::Equal,
        })
}

impl PartialOrd for Incidence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
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
        NodeId(offset + entity.id_index() as u32)
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
    /// occurrence. Stereo elements attach to their site and every ligand-bearing
    /// atom; the incidence type distinguishes the site and ligand roles.
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
        let mut edges: Vec<[u32; 2]> = Vec::new();
        let mut incidences = Vec::new();
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
                for (position, atom) in v.atom_ids().enumerate() {
                    let electrons = match v.electrons() {
                        ElectronCountsForm::Undetermined => NumForm::Undetermined,
                        ElectronCountsForm::Lit(counts) => NumForm::Lit(counts[position]),
                    };
                    push(
                        node,
                        atom.index() as u32,
                        Incidence::AromaticParticipant(electrons),
                    );
                }
                node += 1;
            }
            for v in self.multicenter_bonds().iter() {
                for (position, atom) in v.atom_ids().enumerate() {
                    let electrons = match v.electrons() {
                        ElectronCountsForm::Undetermined => NumForm::Undetermined,
                        ElectronCountsForm::Lit(counts) => NumForm::Lit(counts[position]),
                    };
                    push(
                        node,
                        atom.index() as u32,
                        Incidence::MulticenterParticipant(electrons),
                    );
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
                push(node, v.site_id().index() as u32, Incidence::StereoSite);
                for ligand in v.ligands() {
                    push(
                        node,
                        ligand.atom_id().index() as u32,
                        Incidence::StereoLigand(ligand.kind()),
                    );
                }
                node += 1;
            }
            for v in self.stereo_bonds().iter() {
                push(
                    node,
                    atom_count as u32 + v.site_id().index() as u32,
                    Incidence::StereoSite,
                );
                for ligand in v.ligands() {
                    push(
                        node,
                        ligand.atom_id().index() as u32,
                        Incidence::StereoLigand(ligand.kind()),
                    );
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
    /// Entity-span tags and both sides of positional electron-count values remain available to
    /// canonicalization through the typed entity and incidence values.
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
            for id in self.dative_bonds().relation_ids() {
                for &donor in self.dative_bonds().participants_2(id) {
                    push(node, donor.0, Incidence::DativeDonor);
                }
                push(
                    node,
                    self.dative_bonds().participants_1(id)[0].0,
                    Incidence::DativeAcceptor,
                );
                node += 1;
            }
            for id in self.aromatic_systems().relation_ids() {
                for (position, &atom) in self.aromatic_systems().participants(id).iter().enumerate()
                {
                    push(
                        node,
                        atom.0,
                        Incidence::AromaticParticipantSpan(electron_span(
                            self.aromatic_systems().data(id),
                            position,
                            |value| &value.electrons,
                        )),
                    );
                }
                node += 1;
            }
            for id in self.multicenter_bonds().relation_ids() {
                for (position, &atom) in
                    self.multicenter_bonds().participants(id).iter().enumerate()
                {
                    push(
                        node,
                        atom.0,
                        Incidence::MulticenterParticipantSpan(electron_span(
                            self.multicenter_bonds().data(id),
                            position,
                            |value| &value.electrons,
                        )),
                    );
                }
                node += 1;
            }
            for id in self.noncovalent_bonds().relation_ids() {
                for &atom in self.noncovalent_bonds().participants(id) {
                    push(node, atom.0, Incidence::NoncovalentEndpoint);
                }
                node += 1;
            }
        }

        if stereo {
            for id in self.stereo_atoms().relation_ids() {
                push(
                    node,
                    self.stereo_atoms().participants_1(id)[0].0,
                    Incidence::StereoSite,
                );
                for ligand in self.stereo_atoms().participants_2(id) {
                    push(node, ligand.atom_id.0, Incidence::StereoLigand(ligand.kind));
                }
                node += 1;
            }
            for id in self.stereo_bonds().relation_ids() {
                push(
                    node,
                    atom_count as u32 + self.stereo_bonds().participants_1(id)[0].0,
                    Incidence::StereoSite,
                );
                for ligand in self.stereo_bonds().participants_2(id) {
                    push(node, ligand.atom_id.0, Incidence::StereoLigand(ligand.kind));
                }
                node += 1;
            }
        }

        IncidenceGraph {
            graph: Graph::new(node as usize, &edges),
            entity_counts,
            incidences,
        }
    }
}

fn electron_span<T>(
    span: &EntitySpan<T>,
    position: usize,
    electrons: impl Fn(&T) -> &ElectronCountsForm,
) -> EntitySpan<NumForm> {
    let at = |value: &T| match electrons(value) {
        ElectronCountsForm::Undetermined => NumForm::Undetermined,
        ElectronCountsForm::Lit(counts) => NumForm::Lit(counts[position]),
    };
    match span {
        EntitySpan::Unchanged(value) => EntitySpan::Unchanged(at(value)),
        EntitySpan::Added(value) => EntitySpan::Added(at(value)),
        EntitySpan::Removed(value) => EntitySpan::Removed(at(value)),
        EntitySpan::Modified { lhs, rhs } => EntitySpan::Modified {
            lhs: at(lhs),
            rhs: at(rhs),
        },
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
                AtomId(0),
                AtomId(5),
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
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
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
    #[case::aromatic(EdgeId(8), Incidence::AromaticParticipant(NumForm::Lit(1)))]
    #[case::multicenter(EdgeId(11), Incidence::MulticenterParticipant(NumForm::Lit(2)))]
    #[case::noncovalent(EdgeId(14), Incidence::NoncovalentEndpoint)]
    #[case::stereo_site(EdgeId(16), Incidence::StereoSite)]
    #[case::stereo_ligand(
        EdgeId(19),
        Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen)
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
        Incidence::AromaticParticipant(NumForm::Undetermined),
    )]
    #[case::multicenter(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            multicenter: vec![(vec![AtomId(0)], MulticenterBondForm::default())],
            ..Default::default()
        }),
        Incidence::MulticenterParticipant(NumForm::Undetermined),
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
        assert_eq!(iterator.len(), 26);
        assert_eq!(
            iterator
                .map(|(edge, value)| (edge, value.clone()))
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
                (EdgeId(8), Incidence::AromaticParticipant(NumForm::Lit(1)),),
                (EdgeId(9), Incidence::AromaticParticipant(NumForm::Lit(0)),),
                (EdgeId(10), Incidence::AromaticParticipant(NumForm::Lit(2)),),
                (
                    EdgeId(11),
                    Incidence::MulticenterParticipant(NumForm::Lit(2)),
                ),
                (
                    EdgeId(12),
                    Incidence::MulticenterParticipant(NumForm::Lit(0)),
                ),
                (
                    EdgeId(13),
                    Incidence::MulticenterParticipant(NumForm::Lit(1)),
                ),
                (EdgeId(14), Incidence::NoncovalentEndpoint),
                (EdgeId(15), Incidence::NoncovalentEndpoint),
                (EdgeId(16), Incidence::StereoSite),
                (EdgeId(17), Incidence::StereoLigand(StereoLigandKind::Atom),),
                (EdgeId(18), Incidence::StereoLigand(StereoLigandKind::Atom),),
                (
                    EdgeId(19),
                    Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
                ),
                (
                    EdgeId(20),
                    Incidence::StereoLigand(StereoLigandKind::LonePair),
                ),
                (EdgeId(21), Incidence::StereoSite),
                (EdgeId(22), Incidence::StereoLigand(StereoLigandKind::Atom),),
                (EdgeId(23), Incidence::StereoLigand(StereoLigandKind::Atom),),
                (
                    EdgeId(24),
                    Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
                ),
                (
                    EdgeId(25),
                    Incidence::StereoLigand(StereoLigandKind::LonePair),
                ),
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
    // Stereo nodes attach to their site and to every ligand-bearing atom. Site
    // and virtual-ligand occurrences on the same atom remain parallel edges.
    #[case::stereo_atom(13, vec![0, 1, 1, 1, 2])]
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
    #[case::repeated_virtual_ligand_anchor(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            ([NodeId(0), NodeId(3)], Incidence::StereoSite),
            (
                [NodeId(1), NodeId(3)],
                Incidence::StereoLigand(StereoLigandKind::Atom),
            ),
            (
                [NodeId(2), NodeId(3)],
                Incidence::StereoLigand(StereoLigandKind::Atom),
            ),
            (
                [NodeId(0), NodeId(3)],
                Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            ),
            (
                [NodeId(0), NodeId(3)],
                Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            ),
        ],
    )]
    fn test_molecule_incidence_graph_parallel(
        #[case] molecule: Molecule,
        #[case] expected: Vec<([NodeId; 2], Incidence)>,
    ) {
        let incidence = molecule.incidence_graph(IncidenceLevel::Full);
        assert_eq!(
            incidence
                .incidences()
                .map(|(edge, value)| (incidence.graph().edge_endpoints(edge), value.clone()))
                .collect::<Vec<_>>(),
            expected,
        );
    }
}
