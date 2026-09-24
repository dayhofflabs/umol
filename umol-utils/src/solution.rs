//! Three-valued outcome of an engine pass.
//!
//! `Solution<T, C, U = T>` distinguishes a fully-determined result carrying `T`,
//! an underdetermined outcome carrying `U`, and a chemistry-level contradiction
//! with a typed diagnostic payload `C`. The two-parameter form uses `T` for both
//! non-contradictory outcomes.
//!
//! Engine setup or parameter-table errors travel separately in `Result<_, _>`
//! and never collapse into `Solution`.

/// Semantic outcome with separate determined, underdetermined, and contradiction payloads.
///
/// `U` defaults to `T`; the two-parameter form retains one payload type for both
/// non-contradictory outcomes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Solution<T, C, U = T> {
    Determined(T),
    Underdetermined(U),
    Contradictory(C),
}

impl<T, C, U> Solution<T, C, U> {
    pub fn is_determined(&self) -> bool {
        matches!(self, Self::Determined(_))
    }

    pub fn is_underdetermined(&self) -> bool {
        matches!(self, Self::Underdetermined(_))
    }

    pub fn is_contradictory(&self) -> bool {
        matches!(self, Self::Contradictory(_))
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

    pub fn into_contradiction(self) -> Option<C> {
        match self {
            Self::Contradictory(c) => Some(c),
            Self::Determined(_) | Self::Underdetermined(_) => None,
        }
    }

    /// Transform the contradiction payload type. Used by composite engines
    /// to wrap a sub-engine's contradiction in their union enum.
    pub fn map_contradiction<D, F>(self, f: F) -> Solution<T, D, U>
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

impl<T, C> Solution<T, C> {
    /// Borrow the success payload (Determined or Underdetermined).
    pub fn data(&self) -> Option<&T> {
        match self {
            Self::Determined(v) | Self::Underdetermined(v) => Some(v),
            Self::Contradictory(_) => None,
        }
    }

    /// Extract the success payload (Determined or Underdetermined).
    pub fn into_data(self) -> Option<T> {
        match self {
            Self::Determined(v) | Self::Underdetermined(v) => Some(v),
            Self::Contradictory(_) => None,
        }
    }

    /// Transform the success payload type. Contradiction passes through.
    pub fn map<V, F>(self, f: F) -> Solution<V, C>
    where
        F: FnOnce(T) -> V,
    {
        match self {
            Self::Determined(v) => Solution::Determined(f(v)),
            Self::Underdetermined(v) => Solution::Underdetermined(f(v)),
            Self::Contradictory(c) => Solution::Contradictory(c),
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
        assert_eq!(
            contradictory.contradiction(),
            Some(&Mismatch::Reason("nope"))
        );
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
    fn test_solution_map_contradictory_passes_through(contradictory: Solution<Payload, Mismatch>) {
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

    #[rstest]
    #[case::determined(Solution::Determined(7), (true, false, false))]
    #[case::underdetermined(Solution::Underdetermined(()), (false, true, false))]
    #[case::contradictory(Solution::Contradictory("nope"), (false, false, true))]
    fn test_solution_predicates_distinct(
        #[case] solution: Solution<i32, &'static str, ()>,
        #[case] expected: (bool, bool, bool),
    ) {
        assert_eq!(
            (
                solution.is_determined(),
                solution.is_underdetermined(),
                solution.is_contradictory(),
            ),
            expected
        );
    }

    #[rstest]
    #[case::determined(Solution::Determined(7), Some(7))]
    #[case::underdetermined(Solution::Underdetermined(()), None)]
    #[case::contradictory(Solution::Contradictory("nope"), None)]
    fn test_solution_into_determined_distinct(
        #[case] solution: Solution<i32, &'static str, ()>,
        #[case] expected: Option<i32>,
    ) {
        assert_eq!(solution.into_determined(), expected);
    }

    #[rstest]
    #[case::determined(Solution::Determined(7), Solution::Determined(7))]
    #[case::underdetermined(Solution::Underdetermined(()), Solution::Underdetermined(()))]
    #[case::contradictory(
        Solution::Contradictory("nope"),
        Solution::Contradictory(String::from("nope"))
    )]
    fn test_solution_map_contradiction_distinct(
        #[case] solution: Solution<i32, &'static str, ()>,
        #[case] expected: Solution<i32, String, ()>,
    ) {
        assert_eq!(solution.map_contradiction(String::from), expected);
    }

    #[rstest]
    #[case::determined(Solution::Determined(7), Ok(()))]
    #[case::underdetermined(Solution::Underdetermined(()), Ok(()))]
    #[case::contradictory(Solution::Contradictory("nope"), Err("nope"))]
    fn test_solution_into_observation_distinct(
        #[case] solution: Solution<i32, &'static str, ()>,
        #[case] expected: Result<(), &'static str>,
    ) {
        assert_eq!(solution.into_observation(), expected);
    }

    #[rstest]
    #[case::determined(Solution::Determined(7), Ok(7))]
    #[case::underdetermined(Solution::Underdetermined(()), Err(String::from("underdetermined")))]
    #[case::contradictory(Solution::Contradictory("nope"), Err(String::from("nope")))]
    fn test_solution_into_decisive_distinct(
        #[case] solution: Solution<i32, &'static str, ()>,
        #[case] expected: Result<i32, String>,
    ) {
        assert_eq!(
            solution.into_decisive(String::from("underdetermined")),
            expected
        );
    }
}
