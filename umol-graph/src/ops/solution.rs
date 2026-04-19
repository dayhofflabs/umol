//! Three-valued outcome of a constraint-system operation.

/// Outcome of resolution, validation, or any other unification-shaped operation.
///
/// `Determined` and `Underdetermined` both carry a payload; `Contradictory` does not.
/// Callers that treat both non-contradictory outcomes as success use [`into_result`];
/// callers that demand a fully-determined answer use [`determined`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solution<T> {
    Determined(T),
    Underdetermined(T),
    Contradictory,
}

/// Residual produced when [`Solution::into_result`] collapses to an error.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Contradiction;

impl<T> Solution<T> {
    pub fn is_determined(&self) -> bool {
        matches!(self, Self::Determined(_))
    }

    pub fn is_underdetermined(&self) -> bool {
        matches!(self, Self::Underdetermined(_))
    }

    pub fn is_contradictory(&self) -> bool {
        matches!(self, Self::Contradictory)
    }

    /// Collapse to `Err(Contradiction)` on contradictory, else `Ok(value)`.
    /// Loses the determined/underdetermined distinction.
    pub fn into_result(self) -> Result<T, Contradiction> {
        match self {
            Self::Determined(v) | Self::Underdetermined(v) => Ok(v),
            Self::Contradictory => Err(Contradiction),
        }
    }

    /// `Some(v)` only if `Determined`; `None` otherwise.
    pub fn determined(self) -> Option<T> {
        match self {
            Self::Determined(v) => Some(v),
            Self::Underdetermined(_) | Self::Contradictory => None,
        }
    }

    pub fn map<U, F: FnOnce(T) -> U>(self, f: F) -> Solution<U> {
        match self {
            Self::Determined(v) => Solution::Determined(f(v)),
            Self::Underdetermined(v) => Solution::Underdetermined(f(v)),
            Self::Contradictory => Solution::Contradictory,
        }
    }
}
