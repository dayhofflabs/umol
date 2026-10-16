//! Dense entity ids and contiguous id ranges.

use std::iter::FusedIterator;
use std::marker::PhantomData;
use std::ops::{Add, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub u32);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeId(pub u32);

impl NodeId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<NodeId> for usize {
    fn from(id: NodeId) -> Self {
        id.0 as usize
    }
}

impl From<usize> for NodeId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

/// Dense positional offset, for renumbering a dense id space.
impl Add<usize> for NodeId {
    type Output = Self;

    fn add(self, offset: usize) -> Self {
        Self(self.0 + offset as u32)
    }
}

impl Sub<usize> for NodeId {
    type Output = Self;

    fn sub(self, offset: usize) -> Self {
        Self(self.0 - offset as u32)
    }
}

impl EdgeId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<EdgeId> for usize {
    fn from(id: EdgeId) -> Self {
        id.0 as usize
    }
}

impl From<usize> for EdgeId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

/// Dense positional offset, for renumbering a dense id space.
impl Add<usize> for EdgeId {
    type Output = Self;

    fn add(self, offset: usize) -> Self {
        Self(self.0 + offset as u32)
    }
}

impl Sub<usize> for EdgeId {
    type Output = Self;

    fn sub(self, offset: usize) -> Self {
        Self(self.0 - offset as u32)
    }
}

/// Dense row index within a relation set.
///
/// This value carries no set identity or range validation. Indexed operations check it
/// against their receiver; conversion from usize narrows to u32.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    /// Return the row index as usize.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<RelationId> for usize {
    fn from(id: RelationId) -> Self {
        id.0 as usize
    }
}

impl From<usize> for RelationId {
    fn from(index: usize) -> Self {
        Self(index as u32)
    }
}

impl Add<usize> for RelationId {
    type Output = Self;

    fn add(self, offset: usize) -> Self {
        Self(self.0 + offset as u32)
    }
}

impl Sub<usize> for RelationId {
    type Output = Self;

    fn sub(self, offset: usize) -> Self {
        Self(self.0 - offset as u32)
    }
}

/// The contiguous ids of entities appended in one batch, in order. Positions `start..end`
/// convert to ids on iteration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IdRange<Id> {
    start: usize,
    end: usize,
    id: PhantomData<fn() -> Id>,
}

impl<Id> IdRange<Id> {
    /// The ids at positions `start..end`.
    ///
    /// # Panics
    ///
    /// Panics if `start` exceeds `end`.
    pub fn new(start: usize, end: usize) -> Self {
        assert!(start <= end, "id range start {start} exceeds its end {end}");
        Self {
            start,
            end,
            id: PhantomData,
        }
    }
}

impl<Id: From<usize>> Iterator for IdRange<Id> {
    type Item = Id;

    fn next(&mut self) -> Option<Id> {
        (self.start < self.end).then(|| {
            self.start += 1;
            Id::from(self.start - 1)
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.end - self.start;
        (len, Some(len))
    }
}

impl<Id: From<usize>> DoubleEndedIterator for IdRange<Id> {
    fn next_back(&mut self) -> Option<Id> {
        (self.start < self.end).then(|| {
            self.end -= 1;
            Id::from(self.end)
        })
    }
}

impl<Id: From<usize>> ExactSizeIterator for IdRange<Id> {}

impl<Id: From<usize>> FusedIterator for IdRange<Id> {}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::empty(IdRange::new(3, 3), vec![])]
    #[case::ids(IdRange::new(1, 4), vec![NodeId(1), NodeId(2), NodeId(3)])]
    fn test_id_range_next(#[case] range: IdRange<NodeId>, #[case] expected: Vec<NodeId>) {
        assert_eq!(range.len(), expected.len());
        assert_eq!(range.collect::<Vec<_>>(), expected);
    }

    #[rstest]
    fn test_id_range_next_back() {
        assert_eq!(
            IdRange::<NodeId>::new(1, 4).rev().collect::<Vec<_>>(),
            vec![NodeId(3), NodeId(2), NodeId(1)],
        );
    }

    #[rstest]
    #[should_panic(expected = "id range start 2 exceeds its end 1")]
    fn test_id_range_new_error() {
        IdRange::<NodeId>::new(2, 1);
    }
}
