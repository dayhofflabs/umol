//! Arbitrary-degree permutations.

use crate::PermutationError;

/// A permutation with runtime degree, in one-line notation.
///
/// Unlike [`Permutation`](crate::Permutation), this action has no fixed maximum degree and is not
/// `Copy`. Its direction is `new[i] = old[action[i]]`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DynPermutation {
    image: Box<[usize]>,
}

impl DynPermutation {
    /// The identity action of `degree`.
    pub fn identity(degree: usize) -> Self {
        Self {
            image: (0..degree).collect(),
        }
    }

    /// The number of positions moved by this action.
    pub fn degree(&self) -> usize {
        self.image.len()
    }

    /// The one-line image of this action.
    pub fn image(&self) -> &[usize] {
        &self.image
    }

    /// Reorder `items` by this action, returning `None` when their length differs from its degree.
    pub fn act<T: Clone>(&self, items: &[T]) -> Option<Vec<T>> {
        (items.len() == self.degree()).then(|| {
            self.image
                .iter()
                .map(|&position| items[position].clone())
                .collect()
        })
    }

    /// Function composition `self ∘ other`.
    ///
    /// Returns `None` when the actions have different degrees.
    pub fn compose(&self, other: &Self) -> Option<Self> {
        (self.degree() == other.degree()).then(|| Self {
            image: other
                .image
                .iter()
                .map(|&position| self.image[position])
                .collect(),
        })
    }

    /// The inverse action.
    pub fn inverse(&self) -> Self {
        let mut image = vec![0; self.degree()];
        for (position, &source) in self.image.iter().enumerate() {
            image[source] = position;
        }
        Self {
            image: image.into_boxed_slice(),
        }
    }

    /// The unique action carrying `from` to `to`, such that `action.act(from) == Some(to)`.
    ///
    /// Returns `None` when the slices have different lengths, are not reorderings of the same
    /// values, or repeated equal values make the action non-unique.
    pub fn between<T: Eq>(from: &[T], to: &[T]) -> Option<Self> {
        if from.len() != to.len() {
            return None;
        }
        let mut image = Vec::with_capacity(from.len());
        let mut used = vec![false; from.len()];
        for target in to {
            let source = from.iter().position(|value| value == target)?;
            if used[source] {
                return None;
            }
            used[source] = true;
            image.push(source);
        }
        Some(Self {
            image: image.into_boxed_slice(),
        })
    }
}

impl TryFrom<&[usize]> for DynPermutation {
    type Error = PermutationError;

    fn try_from(image: &[usize]) -> Result<Self, Self::Error> {
        let degree = image.len();
        let mut seen = vec![false; degree];
        for (position, &value) in image.iter().enumerate() {
            if value >= degree {
                return Err(PermutationError::ImageValueOutOfRange {
                    position,
                    value,
                    degree,
                });
            }
            if seen[value] {
                return Err(PermutationError::DuplicateImageValue { value });
            }
            seen[value] = true;
        }
        Ok(Self {
            image: image.into(),
        })
    }
}

impl TryFrom<Vec<usize>> for DynPermutation {
    type Error = PermutationError;

    fn try_from(image: Vec<usize>) -> Result<Self, Self::Error> {
        Self::try_from(image.as_slice()).map(|_| Self {
            image: image.into_boxed_slice(),
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::empty(vec![], vec![])]
    #[case::nonidentity(vec![2, 0, 1], vec![2, 0, 1])]
    #[case::above_stereo_limit(
        vec![6, 5, 4, 3, 2, 1, 0],
        vec![6, 5, 4, 3, 2, 1, 0],
    )]
    fn test_dyn_permutation_try_from(#[case] image: Vec<usize>, #[case] expected: Vec<usize>) {
        assert_eq!(
            DynPermutation::try_from(image)
                .expect("case is a permutation")
                .image(),
            expected,
        );
    }

    #[rstest]
    #[case::out_of_range(
        vec![0, 3, 1],
        PermutationError::ImageValueOutOfRange {
            position: 1,
            value: 3,
            degree: 3,
        },
    )]
    #[case::duplicate(
        vec![0, 1, 1],
        PermutationError::DuplicateImageValue { value: 1 },
    )]
    fn test_dyn_permutation_try_from_error(
        #[case] image: Vec<usize>,
        #[case] expected: PermutationError,
    ) {
        assert_eq!(DynPermutation::try_from(image), Err(expected));
    }

    #[rstest]
    #[case::empty(0, Vec::<u8>::new())]
    #[case::populated(8, (0_u8..8).collect())]
    fn test_dyn_permutation_identity(#[case] degree: usize, #[case] values: Vec<u8>) {
        let identity = DynPermutation::identity(degree);
        assert_eq!(identity.degree(), degree);
        assert_eq!(identity.act(&values), Some(values));
    }

    #[rstest]
    #[case::compatible(
        vec![2, 0, 3, 1],
        vec![10, 11, 12, 13],
        Some(vec![12, 10, 13, 11]),
    )]
    #[case::short(vec![1, 0, 2], vec![10, 11], None)]
    fn test_dyn_permutation_act(
        #[case] image: Vec<usize>,
        #[case] values: Vec<u8>,
        #[case] expected: Option<Vec<u8>>,
    ) {
        let action = DynPermutation::try_from(image).expect("case is a permutation");
        assert_eq!(action.act(&values), expected);
    }

    #[rstest]
    #[case::compatible(
        vec![1, 2, 0],
        vec![2, 1, 0],
        Some(vec![0, 2, 1]),
    )]
    #[case::degree(vec![1, 0], vec![1, 2, 0], None)]
    fn test_dyn_permutation_compose(
        #[case] left: Vec<usize>,
        #[case] right: Vec<usize>,
        #[case] expected: Option<Vec<usize>>,
    ) {
        let left = DynPermutation::try_from(left).expect("case left is a permutation");
        let right = DynPermutation::try_from(right).expect("case right is a permutation");
        assert_eq!(
            left.compose(&right)
                .map(|permutation| permutation.image().to_vec()),
            expected,
        );
    }

    #[rstest]
    #[case::cycle(vec![1, 2, 0], vec![2, 0, 1])]
    #[case::transpositions(vec![1, 0, 3, 2], vec![1, 0, 3, 2])]
    fn test_dyn_permutation_inverse(#[case] image: Vec<usize>, #[case] expected: Vec<usize>) {
        let action = DynPermutation::try_from(image).expect("case is a permutation");
        assert_eq!(action.inverse().image(), expected);
        assert_eq!(
            action.compose(&action.inverse()),
            Some(DynPermutation::identity(action.degree())),
        );
    }

    #[rstest]
    #[case::nonidentity(
        &['a', 'b', 'c'],
        &['c', 'a', 'b'],
        Some(vec![2, 0, 1]),
    )]
    #[case::identity(&['a', 'b'], &['a', 'b'], Some(vec![0, 1]))]
    #[case::length(&['a', 'b'], &['a'], None)]
    #[case::membership(&['a', 'b'], &['a', 'c'], None)]
    #[case::repetition(&['a', 'a'], &['a', 'a'], None)]
    fn test_dyn_permutation_between(
        #[case] from: &[char],
        #[case] to: &[char],
        #[case] expected: Option<Vec<usize>>,
    ) {
        assert_eq!(
            DynPermutation::between(from, to).map(|permutation| permutation.image().to_vec()),
            expected,
        );
    }
}
