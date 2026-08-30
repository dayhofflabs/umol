use rstest::rstest;
use umol_perm::{Permutation, PermutationError, MAX_DEGREE};

#[rstest]
fn test_max_degree() {
    let maximum_image = (0..MAX_DEGREE).collect::<Vec<_>>();
    let excessive_image = (0..=MAX_DEGREE).collect::<Vec<_>>();

    assert_eq!(
        Permutation::try_from(maximum_image.as_slice()).map(Permutation::degree),
        Ok(MAX_DEGREE),
    );
    assert_eq!(
        Permutation::try_from(excessive_image.as_slice()),
        Err(PermutationError::ImageTooLong {
            length: MAX_DEGREE + 1,
            maximum: MAX_DEGREE,
        }),
    );
}
