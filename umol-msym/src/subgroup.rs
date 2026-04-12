//! Subgroups of point groups.

use std::ptr;

use crate::point_group::PointGroup;
use crate::types::SchoenfliesSymbol;

#[derive(Debug, Clone)]
pub(crate) struct SubgroupData {
    pub(crate) symbol: SchoenfliesSymbol,
    pub(crate) name: String,
    pub(crate) order: usize,
    /// How many inequivalent embeddings of this subgroup type exist in the parent.
    pub(crate) multiplicity: usize,
}

/// A specific subgroup embedding within a parent point group. Holds the parent
/// group reference, an opaque index for FFI, and the subgroup's operation data.
#[derive(Debug, Clone)]
pub struct Subgroup {
    parent: &'static PointGroup,
    index: usize,
    data: SubgroupData,
}

impl Subgroup {
    pub(crate) fn new(parent: &'static PointGroup, index: usize, data: SubgroupData) -> Self {
        Self {
            parent,
            index,
            data,
        }
    }

    pub fn parent(&self) -> &'static PointGroup {
        self.parent
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn symbol(&self) -> SchoenfliesSymbol {
        self.data.symbol
    }

    pub fn name(&self) -> &str {
        &self.data.name
    }

    pub fn order(&self) -> usize {
        self.data.order
    }

    /// How many inequivalent embeddings of this subgroup type exist in the parent.
    /// Returns 1 when this is the only embedding of its type.
    pub fn multiplicity(&self) -> usize {
        self.data.multiplicity
    }
}

impl PartialEq for Subgroup {
    fn eq(&self, other: &Self) -> bool {
        ptr::eq(self.parent, other.parent) && self.index == other.index
    }
}

impl Eq for Subgroup {}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::types::SchoenfliesSymbol;

    #[rstest]
    fn test_subgroup_accessors() {
        let parent = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let sg = Subgroup::new(
            parent,
            3,
            SubgroupData {
                symbol: SchoenfliesSymbol::Cn(2),
                name: "C2".into(),
                order: 2,
                multiplicity: 1,
            },
        );
        assert_eq!(sg.parent().symbol(), SchoenfliesSymbol::Cnv(2));
        assert_eq!(sg.symbol(), SchoenfliesSymbol::Cn(2));
        assert_eq!(sg.name(), "C2");
        assert_eq!(sg.order(), 2);
        assert_eq!(sg.multiplicity(), 1);
        assert_eq!(sg.index(), 3);
    }
}
