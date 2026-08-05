//! Incidence (Levi) graph: relations lifted to pseudonodes for symmetry analysis.

use bitflags::bitflags;
use strum::EnumCount;
use umol_graph_core::{Graph, NodeId};

use super::entity::{Entity, EntityKind};
use super::molecule::MoleculeAst;
#[cfg(test)]
use super::molecule::MoleculeEntries;

bitflags! {
    /// Which relation kinds become pseudonodes in [`MoleculeAst::incidence_graph`].
    /// Atoms and localized bonds are always present (the base topology); these
    /// toggle the rest.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct IncidenceNodeSelection: u8 {
        const OVERLAYS = 1 << 0; // dative, aromatic, multicenter, noncovalent
        const STEREO = 1 << 1; // stereo atoms and stereo bonds
    }
}

impl IncidenceNodeSelection {
    /// Atoms + localized bonds only.
    pub fn topological() -> Self {
        Self::empty()
    }

    /// Topological + all overlays — the full constitution.
    pub fn constitution() -> Self {
        Self::OVERLAYS
    }

    /// Constitution + stereo elements.
    pub fn full() -> Self {
        Self::OVERLAYS | Self::STEREO
    }
}

/// A molecule's incidence graph: atoms plus a pseudonode per selected relation,
/// each wired to its participants. `entity(node)` recovers the entity behind a
/// node; atoms occupy `0..atom_count` (node index = atom index), so localized
/// bond `BondId(k)` is the pseudonode at `atom_count + k`.
///
/// Bond direction is not encoded structurally: the coloring separates the
/// endpoints of a directed bond (a dative donor and acceptor are never
/// automorphism-equivalent), and the direction itself is retained in the AST.
#[derive(Clone, Debug)]
pub struct IncidenceGraph {
    graph: Graph,
    // Per-kind block sizes, indexed by `EntityKind as usize` (node-layout order).
    // An entity's id is its index within its block, so node↔entity is offset
    // arithmetic — no per-node table.
    entity_counts: [u32; EntityKind::COUNT],
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
}

impl MoleculeAst {
    /// Build the incidence (Levi) graph over the selected relation kinds. Localized
    /// bonds and overlays become pseudonodes wired to their participant atoms;
    /// stereo elements attach to their site only (an atom, or the site bond's
    /// pseudonode) — the ligand topology is already present via the bonds, so the
    /// only new information a stereo node carries is its site and (at colour time)
    /// its stereo label.
    pub fn incidence_graph(&self, selection: IncidenceNodeSelection) -> IncidenceGraph {
        let atom_count = self.raw_graph().node_count();
        let overlays = selection.contains(IncidenceNodeSelection::OVERLAYS);
        let stereo = selection.contains(IncidenceNodeSelection::STEREO);

        let entity_counts = [
            atom_count as u32,
            self.bonds().ids().count() as u32,
            if overlays {
                self.dative_bonds().count() as u32
            } else {
                0
            },
            if overlays {
                self.aromatic_systems().count() as u32
            } else {
                0
            },
            if overlays {
                self.multicenter_bonds().count() as u32
            } else {
                0
            },
            if overlays {
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

        // Pseudonodes follow the atom block in the same fixed order as `entity_counts`,
        // each wired to its participant atoms. Bonds come first, fixing BondId(k)
        // at atom_count + k — relied on by the stereo-bond site link. Stereo nodes
        // attach to their site only (atom, or the site bond's pseudonode).
        let mut edges: Vec<[u32; 2]> = Vec::new();
        let mut node = atom_count as u32;

        for id in self.bonds().ids() {
            let [a, b] = self.bond(id).atom_ids();
            edges.push([node, a.index() as u32]);
            edges.push([node, b.index() as u32]);
            node += 1;
        }

        if overlays {
            for v in self.dative_bonds().iter() {
                for a in v.atom_ids() {
                    edges.push([node, a.index() as u32]);
                }
                node += 1;
            }
            for v in self.aromatic_systems().iter() {
                for a in v.atom_ids() {
                    edges.push([node, a.index() as u32]);
                }
                node += 1;
            }
            for v in self.multicenter_bonds().iter() {
                for a in v.atom_ids() {
                    edges.push([node, a.index() as u32]);
                }
                node += 1;
            }
            for v in self.noncovalent_bonds().iter() {
                for a in v.atom_ids() {
                    edges.push([node, a.index() as u32]);
                }
                node += 1;
            }
        }

        if stereo {
            for v in self.stereo_atoms().iter() {
                edges.push([node, v.site_id().index() as u32]);
                node += 1;
            }
            for v in self.stereo_bonds().iter() {
                edges.push([node, atom_count as u32 + v.site_id().index() as u32]);
                node += 1;
            }
        }

        let graph = Graph::new(node as usize, &edges);
        IncidenceGraph {
            graph,
            entity_counts,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoAtomId, StereoBondId,
    };
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::stereo::{StereoAtomAst, StereoBondAst, StereoCoset, StereoKind};

    // Six carbons; chain bonds 0-1-2-3 (BondId 0,1,2); a dative 0→3; an aromatic
    // system {0,1,2}; a multicenter {3,4,5}; a noncovalent 0···5; a stereo atom on
    // site 1; a stereo bond on site BondId(1).
    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(3), DativeBondAst::from_order(1))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            multicenter: vec![(
                vec![AtomId(3), AtomId(4), AtomId(5)],
                MulticenterBondAst::default(),
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(5),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[rstest]
    #[case::topological(
        IncidenceNodeSelection::topological(),
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
        ],
    )]
    #[case::constitution(
        IncidenceNodeSelection::constitution(),
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
            Entity::DativeBond(DativeBondId(0)), Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(0)),
        ],
    )]
    #[case::full(
        IncidenceNodeSelection::full(),
        vec![
            Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), Entity::Atom(AtomId(2)),
            Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)),
            Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), Entity::Bond(BondId(2)),
            Entity::DativeBond(DativeBondId(0)), Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(0)),
            Entity::StereoAtom(StereoAtomId(0)), Entity::StereoBond(StereoBondId(0)),
        ],
    )]
    fn test_molecule_ast_incidence_graph(
        molecule: MoleculeAst,
        #[case] selection: IncidenceNodeSelection,
        #[case] expected: Vec<Entity>,
    ) {
        let inc = molecule.incidence_graph(selection);
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
    // Localized bonds are pseudonodes wired to both endpoints (not atom-atom edges).
    #[case::bond(6, vec![0, 1])]
    // Overlays wire to all participant atoms.
    #[case::dative(9, vec![0, 3])]
    #[case::aromatic(10, vec![0, 1, 2])]
    #[case::multicenter(11, vec![3, 4, 5])]
    #[case::noncovalent(12, vec![0, 5])]
    // Stereo nodes attach to their site only: the stereo atom to site atom 1, the
    // stereo bond to the pseudonode of its site BondId(1) = node 6 + 1 = 7.
    #[case::stereo_atom(13, vec![1])]
    #[case::stereo_bond(14, vec![7])]
    fn test_molecule_ast_incidence_graph_neighbors(
        molecule: MoleculeAst,
        #[case] node: u32,
        #[case] expected: Vec<u32>,
    ) {
        let inc = molecule.incidence_graph(IncidenceNodeSelection::full());
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
}
