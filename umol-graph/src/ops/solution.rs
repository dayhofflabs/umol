//! Three-valued outcome of an engine pass.
//!
//! `Solution<T, C>` distinguishes a fully-determined result, a
//! partially-resolved one (the engine ran but the AST is still not ground),
//! and a chemistry-level contradiction with a typed diagnostic payload `C`.
//!
//! Engine setup or parameter-table errors travel separately in `Result<_, _>`
//! and never collapse into `Solution`.

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Solution<T, C> {
    Determined(T),
    Underdetermined(T),
    Contradictory(C),
}

impl<T, C> Solution<T, C> {
    pub fn is_determined(&self) -> bool {
        matches!(self, Self::Determined(_))
    }

    pub fn is_underdetermined(&self) -> bool {
        matches!(self, Self::Underdetermined(_))
    }

    pub fn is_contradictory(&self) -> bool {
        matches!(self, Self::Contradictory(_))
    }

    /// Borrow the success payload (Determined or Underdetermined).
    pub fn data(&self) -> Option<&T> {
        match self {
            Self::Determined(v) | Self::Underdetermined(v) => Some(v),
            Self::Contradictory(_) => None,
        }
    }

    pub fn contradiction(&self) -> Option<&C> {
        match self {
            Self::Contradictory(c) => Some(c),
            Self::Determined(_) | Self::Underdetermined(_) => None,
        }
    }

    /// Extract the determined value only; `None` for the other two variants.
    pub fn into_determined(self) -> Option<T> {
        match self {
            Self::Determined(v) => Some(v),
            Self::Underdetermined(_) | Self::Contradictory(_) => None,
        }
    }

    /// Extract the success payload (Determined or Underdetermined).
    pub fn into_data(self) -> Option<T> {
        match self {
            Self::Determined(v) | Self::Underdetermined(v) => Some(v),
            Self::Contradictory(_) => None,
        }
    }

    pub fn into_contradiction(self) -> Option<C> {
        match self {
            Self::Contradictory(c) => Some(c),
            Self::Determined(_) | Self::Underdetermined(_) => None,
        }
    }

    /// Transform the success payload type. Contradiction passes through.
    pub fn map<U, F>(self, f: F) -> Solution<U, C>
    where
        F: FnOnce(T) -> U,
    {
        match self {
            Self::Determined(v) => Solution::Determined(f(v)),
            Self::Underdetermined(v) => Solution::Underdetermined(f(v)),
            Self::Contradictory(c) => Solution::Contradictory(c),
        }
    }

    /// Transform the contradiction payload type. Used by composite engines
    /// to wrap a sub-engine's contradiction in their union enum.
    pub fn map_contradiction<D, F>(self, f: F) -> Solution<T, D>
    where
        F: FnOnce(C) -> D,
    {
        match self {
            Self::Determined(v) => Solution::Determined(v),
            Self::Underdetermined(v) => Solution::Underdetermined(v),
            Self::Contradictory(c) => Solution::Contradictory(f(c)),
        }
    }

    /// Validator-style mapping: Determined and Underdetermined are both
    /// successful observations; Contradictory is the only failure. The
    /// payload is discarded.
    pub fn into_observation(self) -> Result<(), C> {
        match self {
            Self::Determined(_) | Self::Underdetermined(_) => Ok(()),
            Self::Contradictory(c) => Err(c),
        }
    }

    /// Transformer-style mapping: only Determined is successful; both
    /// Underdetermined and Contradictory map to `Err`. The caller supplies
    /// the error value used for the Underdetermined case.
    pub fn into_decisive<E>(self, on_underdetermined: E) -> Result<T, E>
    where
        C: Into<E>,
    {
        match self {
            Self::Determined(v) => Ok(v),
            Self::Underdetermined(_) => Err(on_underdetermined),
            Self::Contradictory(c) => Err(c.into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct Payload(i32);

    #[derive(Clone, Debug, PartialEq, Eq)]
    enum Mismatch {
        Reason(&'static str),
    }

    #[fixture]
    fn determined() -> Solution<Payload, Mismatch> {
        Solution::Determined(Payload(7))
    }

    #[fixture]
    fn underdetermined() -> Solution<Payload, Mismatch> {
        Solution::Underdetermined(Payload(3))
    }

    #[fixture]
    fn contradictory() -> Solution<Payload, Mismatch> {
        Solution::Contradictory(Mismatch::Reason("nope"))
    }

    #[rstest]
    #[case::determined(determined(), true, false, false)]
    #[case::underdetermined(underdetermined(), false, true, false)]
    #[case::contradictory(contradictory(), false, false, true)]
    fn test_solution_predicates(
        #[case] s: Solution<Payload, Mismatch>,
        #[case] det: bool,
        #[case] und: bool,
        #[case] con: bool,
    ) {
        assert_eq!(s.is_determined(), det);
        assert_eq!(s.is_underdetermined(), und);
        assert_eq!(s.is_contradictory(), con);
    }

    #[rstest]
    fn test_solution_data_determined(determined: Solution<Payload, Mismatch>) {
        assert_eq!(determined.data(), Some(&Payload(7)));
        assert_eq!(determined.contradiction(), None);
    }

    #[rstest]
    fn test_solution_data_underdetermined(underdetermined: Solution<Payload, Mismatch>) {
        assert_eq!(underdetermined.data(), Some(&Payload(3)));
        assert_eq!(underdetermined.contradiction(), None);
    }

    #[rstest]
    fn test_solution_data_contradictory(contradictory: Solution<Payload, Mismatch>) {
        assert_eq!(contradictory.data(), None);
        assert_eq!(contradictory.contradiction(), Some(&Mismatch::Reason("nope")));
    }

    #[rstest]
    fn test_solution_into_determined_only_unwraps_determined(
        determined: Solution<Payload, Mismatch>,
        underdetermined: Solution<Payload, Mismatch>,
        contradictory: Solution<Payload, Mismatch>,
    ) {
        assert_eq!(determined.into_determined(), Some(Payload(7)));
        assert_eq!(underdetermined.into_determined(), None);
        assert_eq!(contradictory.into_determined(), None);
    }

    #[rstest]
    fn test_solution_into_data_unwraps_both_success_variants(
        determined: Solution<Payload, Mismatch>,
        underdetermined: Solution<Payload, Mismatch>,
        contradictory: Solution<Payload, Mismatch>,
    ) {
        assert_eq!(determined.into_data(), Some(Payload(7)));
        assert_eq!(underdetermined.into_data(), Some(Payload(3)));
        assert_eq!(contradictory.into_data(), None);
    }

    #[rstest]
    fn test_solution_into_contradiction_only_unwraps_contradictory(
        determined: Solution<Payload, Mismatch>,
        underdetermined: Solution<Payload, Mismatch>,
        contradictory: Solution<Payload, Mismatch>,
    ) {
        assert_eq!(determined.into_contradiction(), None);
        assert_eq!(underdetermined.into_contradiction(), None);
        assert_eq!(
            contradictory.into_contradiction(),
            Some(Mismatch::Reason("nope")),
        );
    }

    #[rstest]
    fn test_solution_map_determined(determined: Solution<Payload, Mismatch>) {
        let mapped = determined.map(|Payload(n)| Payload(n + 1));
        assert_eq!(mapped, Solution::Determined(Payload(8)));
    }

    #[rstest]
    fn test_solution_map_underdetermined(underdetermined: Solution<Payload, Mismatch>) {
        let mapped = underdetermined.map(|Payload(n)| Payload(n * 2));
        assert_eq!(mapped, Solution::Underdetermined(Payload(6)));
    }

    #[rstest]
    fn test_solution_map_contradictory_passes_through(
        contradictory: Solution<Payload, Mismatch>,
    ) {
        let mapped = contradictory.map(|Payload(n)| Payload(n + 100));
        assert_eq!(mapped, Solution::Contradictory(Mismatch::Reason("nope")));
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    enum Wrapped {
        FromMismatch(Mismatch),
    }

    #[rstest]
    fn test_solution_map_contradiction_determined_passes_through(
        determined: Solution<Payload, Mismatch>,
    ) {
        let mapped = determined.map_contradiction(Wrapped::FromMismatch);
        assert_eq!(mapped, Solution::Determined(Payload(7)));
    }

    #[rstest]
    fn test_solution_map_contradiction_underdetermined_passes_through(
        underdetermined: Solution<Payload, Mismatch>,
    ) {
        let mapped = underdetermined.map_contradiction(Wrapped::FromMismatch);
        assert_eq!(mapped, Solution::Underdetermined(Payload(3)));
    }

    #[rstest]
    fn test_solution_map_contradiction_wraps_contradiction(
        contradictory: Solution<Payload, Mismatch>,
    ) {
        let mapped = contradictory.map_contradiction(Wrapped::FromMismatch);
        assert_eq!(
            mapped,
            Solution::Contradictory(Wrapped::FromMismatch(Mismatch::Reason("nope"))),
        );
    }

    #[rstest]
    fn test_solution_into_observation(
        determined: Solution<Payload, Mismatch>,
        underdetermined: Solution<Payload, Mismatch>,
        contradictory: Solution<Payload, Mismatch>,
    ) {
        assert_eq!(determined.into_observation(), Ok(()));
        assert_eq!(underdetermined.into_observation(), Ok(()));
        assert_eq!(
            contradictory.into_observation(),
            Err(Mismatch::Reason("nope")),
        );
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    enum DecisiveError {
        Undetermined,
        Mismatch(Mismatch),
    }

    impl From<Mismatch> for DecisiveError {
        fn from(m: Mismatch) -> Self {
            DecisiveError::Mismatch(m)
        }
    }

    #[rstest]
    fn test_solution_into_decisive(
        determined: Solution<Payload, Mismatch>,
        underdetermined: Solution<Payload, Mismatch>,
        contradictory: Solution<Payload, Mismatch>,
    ) {
        assert_eq!(
            determined.into_decisive::<DecisiveError>(DecisiveError::Undetermined),
            Ok(Payload(7)),
        );
        assert_eq!(
            underdetermined.into_decisive::<DecisiveError>(DecisiveError::Undetermined),
            Err(DecisiveError::Undetermined),
        );
        assert_eq!(
            contradictory.into_decisive::<DecisiveError>(DecisiveError::Undetermined),
            Err(DecisiveError::Mismatch(Mismatch::Reason("nope"))),
        );
    }
}
