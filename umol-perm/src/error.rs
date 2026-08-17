//! Errors from parsing and constructing permutation values.

use std::error::Error;
use std::fmt;

/// A malformed runtime permutation image or cycle representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PermutationError {
    /// The one-line image does not fit the fixed permutation representation.
    ImageTooLong { length: usize, maximum: usize },
    /// A one-line image value lies outside the image's degree.
    ImageValueOutOfRange {
        position: usize,
        value: usize,
        degree: usize,
    },
    /// A one-line image contains the same value more than once.
    DuplicateImageValue { value: usize },
    /// A cycle point lies outside the explicitly supplied degree.
    CyclePointOutOfRange {
        cycle: usize,
        position: usize,
        point: usize,
        degree: usize,
    },
    /// A point occurs more than once in the cycle representation.
    DuplicateCyclePoint { point: usize },
}

impl fmt::Display for PermutationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ImageTooLong { length, maximum } => {
                write!(
                    f,
                    "permutation image length {length} exceeds maximum {maximum}"
                )
            }
            Self::ImageValueOutOfRange {
                position,
                value,
                degree,
            } => write!(
                f,
                "permutation image value {value} at position {position} is outside 0..{degree}"
            ),
            Self::DuplicateImageValue { value } => {
                write!(f, "permutation image value {value} occurs more than once")
            }
            Self::CyclePointOutOfRange {
                cycle,
                position,
                point,
                degree,
            } => write!(
                f,
                "cycle point {point} at cycle {cycle}, position {position} is outside 0..{degree}"
            ),
            Self::DuplicateCyclePoint { point } => {
                write!(f, "cycle point {point} occurs more than once")
            }
        }
    }
}

impl Error for PermutationError {}

/// Invalid text supplied to `ClassKey` through [`FromStr::from_str`](std::str::FromStr::from_str).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ParseClassKeyError {
    /// The input is neither a fixed class key nor a recognized class family.
    UnknownClassKey { input: String },
    /// A recognized class family has a missing or malformed degree.
    InvalidDegree { input: String },
    /// The parsed degree does not fit the fixed permutation representation.
    DegreeTooLarge { degree: usize, maximum: usize },
}

impl fmt::Display for ParseClassKeyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownClassKey { input } => write!(f, "unknown class key: {input}"),
            Self::InvalidDegree { input } => {
                write!(f, "invalid degree in class key: {input}")
            }
            Self::DegreeTooLarge { degree, maximum } => {
                write!(f, "class key degree {degree} exceeds maximum {maximum}")
            }
        }
    }
}

impl Error for ParseClassKeyError {}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::{ParseClassKeyError, PermutationError};

    #[rstest]
    #[case::image_too_long(
        PermutationError::ImageTooLong { length: 7, maximum: 6 },
        "permutation image length 7 exceeds maximum 6",
    )]
    #[case::image_value_out_of_range(
        PermutationError::ImageValueOutOfRange { position: 2, value: 4, degree: 4 },
        "permutation image value 4 at position 2 is outside 0..4",
    )]
    #[case::duplicate_image_value(
        PermutationError::DuplicateImageValue { value: 1 },
        "permutation image value 1 occurs more than once",
    )]
    #[case::cycle_point_out_of_range(
        PermutationError::CyclePointOutOfRange { cycle: 1, position: 2, point: 5, degree: 5 },
        "cycle point 5 at cycle 1, position 2 is outside 0..5",
    )]
    #[case::duplicate_cycle_point(
        PermutationError::DuplicateCyclePoint { point: 3 },
        "cycle point 3 occurs more than once",
    )]
    fn test_permutation_error_display(#[case] error: PermutationError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    #[case::unknown_class_key(
        ParseClassKeyError::UnknownClassKey { input: "Xyz3".to_string() },
        "unknown class key: Xyz3",
    )]
    #[case::invalid_degree(
        ParseClassKeyError::InvalidDegree { input: "Sym".to_string() },
        "invalid degree in class key: Sym",
    )]
    #[case::degree_too_large(
        ParseClassKeyError::DegreeTooLarge { degree: 7, maximum: 6 },
        "class key degree 7 exceeds maximum 6",
    )]
    fn test_parse_class_key_error_display(
        #[case] error: ParseClassKeyError,
        #[case] expected: &str,
    ) {
        assert_eq!(error.to_string(), expected);
    }
}
