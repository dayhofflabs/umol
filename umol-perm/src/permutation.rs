//! A single permutation in one-line notation, degree ≤ 6, `Copy`.
//!
//! Conventions:
//! - one-line image: `apply(i)` = σ(i); slots beyond `degree` are identity.
//! - `act(items)[i] = items[apply(i)]` — reorders a sequence by σ.
//! - `compose(other)` = σ ∘ τ, i.e. `apply(other.apply(i))`.
//! - `between(from, to)` = the τ with `act(τ, from) == to`.
//! - `rank`/`unrank` are the Lehmer order — an internal canonical numbering for
//!   coset representatives, NOT the OpenSMILES arrangement index.

const MAX_DEGREE: usize = 6;

/// A permutation of `0..degree` (`degree ≤ 6`) in one-line notation.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Permutation {
    image: [u8; MAX_DEGREE],
    degree: u8,
}

impl Permutation {
    /// The identity permutation of the given degree.
    pub fn identity(degree: usize) -> Self {
        assert!(degree <= MAX_DEGREE);
        let mut image = [0u8; MAX_DEGREE];
        for (i, slot) in image.iter_mut().enumerate() {
            *slot = i as u8;
        }
        Self {
            image,
            degree: degree as u8,
        }
    }

    /// Build from a one-line image; panics unless `image` is a bijection of
    /// `0..degree`.
    pub fn from_image(degree: usize, image: &[u8]) -> Self {
        assert!(degree <= MAX_DEGREE);
        assert_eq!(image.len(), degree);
        let mut seen = [false; MAX_DEGREE];
        for &x in image {
            let x = x as usize;
            assert!(x < degree, "image entry out of range");
            assert!(!seen[x], "image is not a bijection");
            seen[x] = true;
        }
        let mut full = [0u8; MAX_DEGREE];
        for (i, slot) in full.iter_mut().enumerate() {
            *slot = if i < degree { image[i] } else { i as u8 };
        }
        Self {
            image: full,
            degree: degree as u8,
        }
    }

    pub fn degree(self) -> usize {
        self.degree as usize
    }

    /// σ(i).
    pub fn apply(self, i: usize) -> usize {
        self.image[i] as usize
    }

    /// Reorder `items` by σ: `act(items)[i] = items[σ(i)]`.
    pub fn act<T: Copy>(self, items: &[T]) -> Vec<T> {
        (0..self.degree()).map(|i| items[self.apply(i)]).collect()
    }

    /// Function composition σ ∘ τ: `compose(τ).apply(i) == σ(τ(i))`.
    pub fn compose(self, other: Self) -> Self {
        debug_assert_eq!(self.degree, other.degree);
        let mut image = [0u8; MAX_DEGREE];
        for (i, slot) in image.iter_mut().enumerate() {
            *slot = self.image[other.image[i] as usize];
        }
        Self {
            image,
            degree: self.degree,
        }
    }

    /// σ⁻¹.
    pub fn inverse(self) -> Self {
        let mut image = [0u8; MAX_DEGREE];
        for i in 0..MAX_DEGREE {
            image[self.image[i] as usize] = i as u8;
        }
        Self {
            image,
            degree: self.degree,
        }
    }

    /// Parity: `+1` even, `-1` odd.
    pub fn sign(self) -> i8 {
        let d = self.degree();
        let mut inversions = 0usize;
        for i in 0..d {
            for j in (i + 1)..d {
                if self.image[i] > self.image[j] {
                    inversions += 1;
                }
            }
        }
        if inversions % 2 == 0 {
            1
        } else {
            -1
        }
    }

    /// The τ that relabels `from` into `to`: `act(τ, from) == to`. Panics if the
    /// two slices are not orderings of the same set.
    pub fn between<T: Eq>(from: &[T], to: &[T]) -> Self {
        assert_eq!(from.len(), to.len());
        let degree = from.len();
        let mut image = vec![0u8; degree];
        for (i, target) in to.iter().enumerate() {
            let pos = from
                .iter()
                .position(|source| source == target)
                .expect("between: slices are not orderings of the same set");
            image[i] = pos as u8;
        }
        Self::from_image(degree, &image)
    }

    /// Lehmer rank in `0..degree!` — the internal canonical numbering.
    pub fn rank(self) -> usize {
        let d = self.degree();
        let mut rank = 0;
        let mut factorial = 1;
        for i in (0..d).rev() {
            let mut smaller = 0;
            for j in (i + 1)..d {
                if self.image[j] < self.image[i] {
                    smaller += 1;
                }
            }
            rank += smaller * factorial;
            factorial *= d - i;
        }
        rank
    }

    /// The `rank`-th permutation of the given degree in Lehmer order.
    pub fn unrank(degree: usize, mut rank: usize) -> Self {
        let mut available: Vec<u8> = (0..degree as u8).collect();
        let mut image = vec![0u8; degree];
        let mut factorial: usize = (1..=degree).product();
        for slot in image.iter_mut() {
            factorial /= available.len();
            let idx = rank / factorial;
            rank %= factorial;
            *slot = available.remove(idx);
        }
        Self::from_image(degree, &image)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_permutation_identity() {
        let p = Permutation::identity(4);
        assert_eq!(p.act(&[10, 11, 12, 13]), vec![10, 11, 12, 13]);
        assert_eq!(p.sign(), 1);
        assert_eq!(p.rank(), 0);
    }

    #[rstest]
    fn test_permutation_apply() {
        let p = Permutation::from_image(4, &[2, 0, 3, 1]);
        assert_eq!(p.apply(0), 2);
        assert_eq!(p.apply(1), 0);
        assert_eq!(p.apply(3), 1);
    }

    #[rstest]
    #[case::duplicate(vec![0, 1, 1])]
    #[case::out_of_range(vec![0, 1, 3])]
    #[should_panic]
    fn test_permutation_from_image_error(#[case] image: Vec<u8>) {
        Permutation::from_image(image.len(), &image);
    }

    #[rstest]
    fn test_permutation_act() {
        let p = Permutation::from_image(3, &[2, 0, 1]);
        assert_eq!(p.act(&['a', 'b', 'c']), vec!['c', 'a', 'b']);
    }

    #[rstest]
    fn test_permutation_compose() {
        let sigma = Permutation::from_image(3, &[1, 2, 0]);
        let tau = Permutation::from_image(3, &[2, 1, 0]);
        assert_eq!(sigma.compose(tau), Permutation::from_image(3, &[0, 2, 1]));
    }

    #[rstest]
    fn test_permutation_inverse() {
        let p = Permutation::from_image(3, &[1, 2, 0]);
        assert_eq!(p.inverse(), Permutation::from_image(3, &[2, 0, 1]));
        assert_eq!(p.compose(p.inverse()), Permutation::identity(3));
    }

    #[rstest]
    #[case::identity(Permutation::identity(3), 1)]
    #[case::transposition(Permutation::from_image(3, &[1, 0, 2]), -1)]
    #[case::three_cycle(Permutation::from_image(3, &[1, 2, 0]), 1)]
    fn test_permutation_sign(#[case] p: Permutation, #[case] expected: i8) {
        assert_eq!(p.sign(), expected);
    }

    #[rstest]
    fn test_permutation_between() {
        let tau = Permutation::between(&['a', 'b', 'c'], &['c', 'a', 'b']);
        assert_eq!(tau, Permutation::from_image(3, &[2, 0, 1]));
        assert_eq!(tau.act(&['a', 'b', 'c']), vec!['c', 'a', 'b']);
    }

    #[rstest]
    #[case::degree_2(2)]
    #[case::degree_3(3)]
    #[case::degree_4(4)]
    #[case::degree_5(5)]
    fn test_permutation_rank_unrank(#[case] degree: usize) {
        let count: usize = (1..=degree).product();
        for i in 0..count {
            assert_eq!(Permutation::unrank(degree, i).rank(), i);
        }
    }
}
