//! Operation-local SMILES traversal in atom-index order.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::iter;
use std::ops::{ControlFlow, Range};

use umol_graph_core::{
    visit_depth_first, DepthFirstEvent, EdgeId, Neighbor as GraphNeighbor, NodeId,
};

use crate::table_ir::{Molecule, Neighbor};

/// Preorder atom visits and ring encounters over one unchanged molecular table.
/// No chemical fields are interpreted and no connectivity is stored in the table.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct Traversal {
    pub(super) atoms: Vec<AtomVisit>,
    pub(super) rings: Vec<RingVisit>,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) struct AtomVisit {
    pub(super) atom: u32,
    pub(super) parent: Option<Neighbor>,
    /// Exclusive preorder boundary, so a child's subtree can be skipped in constant time.
    pub(super) subtree_end: usize,
    pub(super) rings: Range<usize>,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) struct RingVisit {
    pub(super) neighbor: Neighbor,
    pub(super) label: usize,
    pub(super) opening: bool,
}

impl Traversal {
    /// Visits all atoms by DFS with roots and neighbors in atom-index order.
    ///
    /// Roots have no parent. Children occur in preorder; all but the last are branches.
    /// Ring encounters precede children and follow neighbor atom order at each atom. Labels use the
    /// smallest available positive integer and become reusable after the closing atom.
    /// Actual stereo-ligand encounter order is parent, rings, then children; virtual-ligand
    /// placement and configuration transport belong to atom formatting.
    ///
    /// Correct molecular output requires in-range endpoints and simple connectivity. This
    /// traversal does not certify an open table's integrity or format representability.
    /// Malformed endpoint rows inherit AtomNeighbors' omission behavior without panicking;
    /// loops and parallel bonds remain present for the format boundary to reject.
    pub(super) fn new(molecule: &Molecule) -> Self {
        let neighbors = molecule.atom_neighbors();
        let mut atoms = Vec::with_capacity(molecule.atoms.len());
        let mut positions = vec![0; molecule.atoms.len()];
        let _: ControlFlow<()> = visit_depth_first(
            molecule.atoms.len(),
            molecule.bonds.len(),
            (0..molecule.atoms.len()).map(|atom| NodeId(atom as u32)),
            |node| {
                neighbors
                    .neighbors(node.0)
                    .iter()
                    .map(|neighbor| GraphNeighbor {
                        node: NodeId(neighbor.atom),
                        edge: EdgeId(neighbor.bond),
                    })
            },
            |event| {
                match event {
                    DepthFirstEvent::Discover { node, parent } => {
                        positions[node.index()] = atoms.len();
                        atoms.push(AtomVisit {
                            atom: node.0,
                            parent: parent.map(|parent| Neighbor {
                                atom: parent.node.0,
                                bond: parent.edge.0,
                            }),
                            subtree_end: 0,
                            rings: 0..0,
                        });
                    }
                    DepthFirstEvent::Finish { node } => {
                        atoms[positions[node.index()]].subtree_end = atoms.len();
                    }
                    DepthFirstEvent::NonTreeEdge { .. } | DepthFirstEvent::FinishTree { .. } => {}
                }
                ControlFlow::Continue(())
            },
        );

        let mut rings = Vec::new();
        let mut active = HashMap::new();
        let mut available = BinaryHeap::new();
        let mut next_label = 1;
        for position in 0..atoms.len() {
            let start = rings.len();
            for &neighbor in neighbors.neighbors(atoms[position].atom) {
                let other = positions[neighbor.atom as usize];
                if atoms[position]
                    .parent
                    .is_some_and(|p| p.bond == neighbor.bond)
                    || atoms[other].parent.is_some_and(|p| p.bond == neighbor.bond)
                {
                    continue;
                }
                let opening = position <= other;
                let label = if opening {
                    let label = available.pop().map_or_else(
                        || {
                            let label = next_label;
                            next_label += 1;
                            label
                        },
                        |Reverse(label)| label,
                    );
                    active.insert(neighbor.bond, label);
                    label
                } else {
                    active.remove(&neighbor.bond).expect("earlier ring opening")
                };
                rings.push(RingVisit {
                    neighbor,
                    label,
                    opening,
                });
                if position == other {
                    active.remove(&neighbor.bond);
                    rings.push(RingVisit {
                        neighbor,
                        label,
                        opening: false,
                    });
                }
            }
            atoms[position].rings = start..rings.len();
            for ring in &rings[start..] {
                if !ring.opening {
                    available.push(Reverse(ring.label));
                }
            }
        }
        Self { atoms, rings }
    }

    /// Direct children in output order, skipping descendants via subtree boundaries.
    pub(super) fn children(&self, position: usize) -> impl Iterator<Item = &AtomVisit> {
        let end = self.atoms.get(position).map_or(0, |atom| atom.subtree_end);
        iter::successors(
            position.checked_add(1).filter(|&next| next < end),
            move |&index| {
                let next = self.atoms[index].subtree_end;
                (next < end).then_some(next)
            },
        )
        .map(|index| &self.atoms[index])
    }

    /// Actual neighbors in lexical encounter order for stereo-frame transport.
    pub(super) fn neighbors(&self, position: usize) -> impl Iterator<Item = Neighbor> + '_ {
        self.atoms.get(position).into_iter().flat_map(move |atom| {
            atom.parent
                .into_iter()
                .chain(
                    self.rings[atom.rings.clone()]
                        .iter()
                        .map(|ring| ring.neighbor),
                )
                .chain(self.children(position).filter_map(|child| {
                    child.parent.map(|parent| Neighbor {
                        atom: child.atom,
                        bond: parent.bond,
                    })
                }))
        })
    }

    /// Bond row and directed endpoints in first lexical encounter order.
    /// Ring bonds occur at opening, oriented toward the closing atom.
    pub(super) fn bonds(&self) -> impl Iterator<Item = (u32, [u32; 2])> + '_ {
        self.atoms.iter().flat_map(|atom| {
            atom.parent
                .map(|parent| (parent.bond, [parent.atom, atom.atom]))
                .into_iter()
                .chain(
                    self.rings[atom.rings.clone()]
                        .iter()
                        .filter(|ring| ring.opening)
                        .map(|ring| (ring.neighbor.bond, [atom.atom, ring.neighbor.atom])),
                )
        })
    }
}

#[cfg(test)]
mod tests;
