//! Incidence (Levi) graph: relations lifted to pseudonodes for symmetry analysis.

use bitflags::bitflags;
use umol_graph_core::{Graph, NodeId};

use super::entity::Entity;
use super::ids::AtomId;
use super::molecule::MoleculeAst;

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
    node_entity: Vec<Entity>,
}

impl IncidenceGraph {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn entity(&self, node: NodeId) -> Entity {
        self.node_entity[node.index()]
    }

    pub fn entities(&self) -> &[Entity] {
        &self.node_entity
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
        let mut node_entity: Vec<Entity> =
            (0..atom_count).map(|i| Entity::Atom(AtomId::from(i))).collect();
        let mut edges: Vec<[u32; 2]> = Vec::new();

        let mut add = |entity: Entity, neighbors: &[u32]| {
            let node = node_entity.len() as u32;
            node_entity.push(entity);
            for &n in neighbors {
                edges.push([node, n]);
            }
        };

        // Localized bonds always: one pseudonode per bond, in BondId order, wired
        // to both endpoints. This fixes the bond pseudonode of BondId(k) at
        // atom_count + k, relied on by the stereo-bond site link below.
        for id in self.bonds().ids() {
            let [a, b] = self.bond(id).atom_ids();
            add(Entity::Bond(id), &[a.index() as u32, b.index() as u32]);
        }

        if selection.contains(IncidenceNodeSelection::OVERLAYS) {
            for v in self.dative_bonds().iter() {
                let atoms: Vec<u32> = v.atom_ids().map(|a| a.index() as u32).collect();
                add(Entity::DativeBond(v.id), &atoms);
            }
            for v in self.aromatic_systems().iter() {
                let atoms: Vec<u32> = v.atom_ids().map(|a| a.index() as u32).collect();
                add(Entity::AromaticSystem(v.id), &atoms);
            }
            for v in self.multicenter_bonds().iter() {
                let atoms: Vec<u32> = v.atom_ids().map(|a| a.index() as u32).collect();
                add(Entity::MulticenterBond(v.id), &atoms);
            }
            for v in self.noncovalent_bonds().iter() {
                let atoms: Vec<u32> = v.atom_ids().into_iter().map(|a| a.index() as u32).collect();
                add(Entity::NoncovalentBond(v.id), &atoms);
            }
        }

        if selection.contains(IncidenceNodeSelection::STEREO) {
            for v in self.stereo_atoms().iter() {
                add(Entity::StereoAtom(v.id), &[v.site_id().index() as u32]);
            }
            for v in self.stereo_bonds().iter() {
                let bond_node = (atom_count + v.site_id().index()) as u32;
                add(Entity::StereoBond(v.id), &[bond_node]);
            }
        }

        let graph = Graph::new(node_entity.len(), &edges);
        IncidenceGraph { graph, node_entity }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::ids::{
        AromaticSystemId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
        StereoBondId,
    };
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::stereo::{StereoAtomAst, StereoBondAst, StereoCosetAst, StereoKind};

    // Six carbons; chain bonds 0-1-2-3 (BondId 0,1,2); a dative 0→3; an aromatic
    // system {0,1,2}; a multicenter {3,4,5}; a noncovalent 0···5; a stereo atom on
    // site 1; a stereo bond on site BondId(1).
    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 6],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(0)], AtomId(3), DativeBondAst::from_order(1))],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomId(3), AtomId(4), AtomId(5)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomId(0),
                AtomId(5),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
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
        assert_eq!(inc.entities(), expected.as_slice());
        assert_eq!(inc.graph().node_count(), expected.len());
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
