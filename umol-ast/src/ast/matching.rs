//! Typed matching result wrapper over `umol_graph_core::Matching`.

use umol_graph_core::{Matching, NodeId};

use super::idx::{AtomIdx, BondIdx};

/// Bond-level wrapper over `umol_graph_core::Matching`. Exposes matched
/// bonds and matched-atom membership in terms of `AtomIdx` / `BondIdx`.
#[derive(Clone, Debug)]
pub struct BondMatching(pub(crate) Matching);

impl BondMatching {
    pub fn bonds(&self) -> impl Iterator<Item = BondIdx> + '_ {
        self.0.edges().iter().map(|&e| BondIdx::from(e))
    }

    pub fn size(&self) -> usize {
        self.0.size()
    }

    pub fn is_perfect(&self, atom_count: usize) -> bool {
        self.0.is_perfect(atom_count)
    }

    pub fn mate(&self, atom: AtomIdx) -> Option<AtomIdx> {
        self.0.mate(NodeId::from(atom)).map(AtomIdx::from)
    }

    pub fn is_matched(&self, atom: AtomIdx) -> bool {
        self.0.is_matched(NodeId::from(atom))
    }
}
