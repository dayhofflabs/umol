//! General-purpose vocabulary shared across the crate's engines.

/// Direction of a sort key component: `Ascending` is the field's natural
/// order, `Descending` its reverse.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SortingDirection {
    Ascending,
    Descending,
}
