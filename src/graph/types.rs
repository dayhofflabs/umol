// General types for the graph module

use petgraph::graph::{EdgeIndex, NodeIndex};
use std::fmt::{self, Display};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomIndex(NodeIndex);

impl AtomIndex {
    pub fn new(idx: usize) -> Self {
        Self(NodeIndex::new(idx))
    }

    pub fn index(&self) -> NodeIndex {
        self.0
    }
}


impl From<NodeIndex> for AtomIndex {
    fn from(idx: NodeIndex) -> Self {
        Self(idx)
    }
}

impl From<AtomIndex> for NodeIndex {
    fn from(idx: AtomIndex) -> Self {
        idx.0
    }
}

impl Display for AtomIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0.index())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BondIndex(EdgeIndex);

impl BondIndex {
    pub fn new(idx: EdgeIndex) -> Self {
        Self(idx)
    }

    pub fn index(&self) -> EdgeIndex {
        self.0
    }
}

impl From<EdgeIndex> for BondIndex {
    fn from(idx: EdgeIndex) -> Self {
        Self(idx)
    }
}

impl From<BondIndex> for EdgeIndex {
    fn from(idx: BondIndex) -> Self {
        idx.0
    }
}

impl Display for BondIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0.index())
    }
}
