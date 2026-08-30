//! A single permutation in one-line notation, degree ≤ 6, `Copy`.
//!
//! Conventions:
//! - one-line image: `apply(i)` = σ(i); slots beyond `degree` are identity.
//! - `act(items)[i] = items[apply(i)]` — reorders a sequence by σ.
//! - `compose(other)` = σ ∘ τ, i.e. `apply(other.apply(i))`.
//! - `between(from, to)` = the τ with `act(τ, from) == to`.
//! - `rank`/`unrank` are the Lehmer order — an internal canonical numbering for
//!   coset representatives, NOT the OpenSMILES arrangement index.

use std::fmt;

use crate::error::PermutationError;

/// Maximum degree supported by [`Permutation`] and permutation-backed class keys.
///
/// Checked constructors reject inputs above this representation limit; asserted constructors
/// panic. The limit is not a chemistry-model restriction.
pub const MAX_DEGREE: usize = 6;

/// A permutation of `0..degree` (`degree <= MAX_DEGREE`) in one-line notation.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Permutation {
    image: [u32; MAX_DEGREE],
    degree: u32,
}

impl Permutation {
    /// The identity permutation of the given degree.
    pub fn identity(degree: usize) -> Self {
        assert!(degree <= MAX_DEGREE);
        let mut image = [0u32; MAX_DEGREE];
        for (i, target) in image.iter_mut().enumerate() {
            *target = i as u32;
        }
        Self {
            image,
            degree: degree as u32,
        }
    }

    /// Build from a one-line image; panics unless `image` is a bijection of
    /// `0..image.len()`.
    pub fn from_image(image: &[usize]) -> Self {
        Self::try_from(image).expect("invalid permutation image")
    }

    pub fn degree(self) -> usize {
        self.degree as usize
    }

    /// σ(i). Panics if `i >= degree`.
    pub fn apply(self, i: usize) -> usize {
        assert!(i < self.degree(), "apply index out of range");
        self.image[i] as usize
    }

    /// Reorder `items` by σ: `act(items)[i] = items[σ(i)]`. Panics unless
    /// `items.len() == degree`.
    pub fn act<T: Clone>(self, items: &[T]) -> Vec<T> {
        assert!(
            items.len() == self.degree(),
            "act slice length must equal degree"
        );
        (0..self.degree())
            .map(|i| items[self.apply(i)].clone())
            .collect()
    }

    /// Function composition σ ∘ τ: `compose(τ).apply(i) == σ(τ(i))`. Panics on a
    /// degree mismatch.
    pub fn compose(self, other: Self) -> Self {
        assert_eq!(self.degree, other.degree);
        let mut image = [0u32; MAX_DEGREE];
        for (i, target) in image.iter_mut().enumerate() {
            *target = self.image[other.image[i] as usize];
        }
        Self {
            image,
            degree: self.degree,
        }
    }

    /// σ⁻¹.
    pub fn inverse(self) -> Self {
        let mut image = [0u32; MAX_DEGREE];
        for i in 0..MAX_DEGREE {
            image[self.image[i] as usize] = i as u32;
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
        if inversions.is_multiple_of(2) {
            1
        } else {
            -1
        }
    }

    /// The unique τ that relabels `from` into `to`: `act(τ, from) == to`.
    /// Returns `None` unless exactly one such permutation exists. Repeated equal
    /// values therefore return `None` when their occurrences can be exchanged.
    /// Panics when the common length exceeds the fixed representation maximum.
    pub fn between<T: Eq>(from: &[T], to: &[T]) -> Option<Self> {
        if from.len() != to.len() {
            return None;
        }
        let degree = from.len();
        assert!(degree <= MAX_DEGREE);
        let mut image = vec![0usize; degree];
        let mut used = [false; MAX_DEGREE];
        for (i, target) in to.iter().enumerate() {
            let pos = from.iter().position(|source| source == target)?;
            if used[pos] {
                return None;
            }
            used[pos] = true;
            image[i] = pos;
        }
        Some(Self::from_image(&image))
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

    /// The `rank`-th permutation of the given degree in Lehmer order. Panics unless
    /// `degree <= MAX_DEGREE` and `rank < degree!`.
    pub fn unrank(degree: usize, mut rank: usize) -> Self {
        assert!(degree <= MAX_DEGREE);
        let mut factorial: usize = (1..=degree).product();
        assert!(rank < factorial, "rank out of range");
        let mut available: Vec<usize> = (0..degree).collect();
        let mut image = vec![0usize; degree];
        for target in image.iter_mut() {
            factorial /= available.len();
            let idx = rank / factorial;
            rank %= factorial;
            *target = available.remove(idx);
        }
        Self::from_image(&image)
    }

    /// Build from disjoint cycles: each `[c0,…,ck]` sets `σ(c0)=c1, …, σ(ck)=c0`;
    /// unlisted points are fixed; `[]` is the identity. Returns an error when a
    /// point is out of range or occurs more than once. Panics if `degree` exceeds
    /// the fixed representation maximum.
    pub fn from_cycles(degree: usize, cycles: &[Vec<usize>]) -> Result<Self, PermutationError> {
        assert!(degree <= MAX_DEGREE);
        let mut seen = [false; MAX_DEGREE];
        for (cycle_index, cycle) in cycles.iter().enumerate() {
            for (position, &point) in cycle.iter().enumerate() {
                if point >= degree {
                    return Err(PermutationError::CyclePointOutOfRange {
                        cycle: cycle_index,
                        position,
                        point,
                        degree,
                    });
                }
                if seen[point] {
                    return Err(PermutationError::DuplicateCyclePoint { point });
                }
                seen[point] = true;
            }
        }
        let mut image: Vec<usize> = (0..degree).collect();
        for cycle in cycles {
            let len = cycle.len();
            for (w, &point) in cycle.iter().enumerate() {
                image[point] = cycle[(w + 1) % len];
            }
        }
        Ok(Self::from_image(&image))
    }

    /// Disjoint-cycle decomposition, fixed points dropped. Canonical: each cycle
    /// starts at its least element, cycles ordered by least element; identity → `[]`.
    pub fn cycles(self) -> Vec<Vec<usize>> {
        let degree = self.degree();
        let mut visited = [false; MAX_DEGREE];
        let mut cycles = Vec::new();
        for start in 0..degree {
            if visited[start] || self.apply(start) == start {
                visited[start] = true;
                continue;
            }
            let mut cycle = Vec::new();
            let mut point = start;
            while !visited[point] {
                visited[point] = true;
                cycle.push(point);
                point = self.apply(point);
            }
            cycles.push(cycle);
        }
        cycles
    }
}

impl TryFrom<&[usize]> for Permutation {
    type Error = PermutationError;

    fn try_from(image: &[usize]) -> Result<Self, Self::Error> {
        let degree = image.len();
        if degree > MAX_DEGREE {
            return Err(PermutationError::ImageTooLong {
                length: degree,
                maximum: MAX_DEGREE,
            });
        }
        let mut seen = [false; MAX_DEGREE];
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
        let mut full = [0u32; MAX_DEGREE];
        for (position, target) in full.iter_mut().enumerate() {
            *target = if position < degree {
                image[position] as u32
            } else {
                position as u32
            };
        }
        Ok(Self {
            image: full,
            degree: degree as u32,
        })
    }
}

impl fmt::Display for Permutation {
    /// Product of disjoint cycles, 0-indexed and comma-separated; identity → `()`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let cycles = self.cycles();
        if cycles.is_empty() {
            return write!(f, "()");
        }
        for cycle in cycles {
            write!(f, "(")?;
            for (i, point) in cycle.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{point}")?;
            }
            write!(f, ")")?;
        }
        Ok(())
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
    #[case::duplicate(vec![0, 1, 1])]
    #[case::out_of_range(vec![0, 1, 3])]
    #[should_panic]
    fn test_permutation_from_image_error(#[case] image: Vec<usize>) {
        Permutation::from_image(&image);
    }

    #[rstest]
    fn test_permutation_apply() {
        let p = Permutation::from_image(&[2, 0, 3, 1]);
        assert_eq!(p.apply(0), 2);
        assert_eq!(p.apply(1), 0);
        assert_eq!(p.apply(3), 1);
    }

    #[rstest]
    fn test_permutation_act() {
        let p = Permutation::from_image(&[2, 0, 1]);
        assert_eq!(
            p.act(&[
                String::from("alpha"),
                String::from("beta"),
                String::from("gamma"),
            ]),
            vec![
                String::from("gamma"),
                String::from("alpha"),
                String::from("beta"),
            ],
        );
    }

    #[rstest]
    #[case::short(&[10, 11])]
    #[case::long(&[10, 11, 12, 13])]
    #[should_panic(expected = "act slice length must equal degree")]
    fn test_permutation_act_error(#[case] items: &[u32]) {
        Permutation::from_image(&[2, 0, 1]).act(items);
    }

    #[rstest]
    fn test_permutation_compose() {
        let sigma = Permutation::from_image(&[1, 2, 0]);
        let tau = Permutation::from_image(&[2, 1, 0]);
        assert_eq!(sigma.compose(tau), Permutation::from_image(&[0, 2, 1]));
    }

    #[rstest]
    fn test_permutation_inverse() {
        let p = Permutation::from_image(&[1, 2, 0]);
        assert_eq!(p.inverse(), Permutation::from_image(&[2, 0, 1]));
        assert_eq!(p.compose(p.inverse()), Permutation::identity(3));
    }

    #[rstest]
    #[case::identity(Permutation::identity(3), 1)]
    #[case::transposition(Permutation::from_image(&[1, 0, 2]), -1)]
    #[case::three_cycle(Permutation::from_image(&[1, 2, 0]), 1)]
    fn test_permutation_sign(#[case] p: Permutation, #[case] expected: i8) {
        assert_eq!(p.sign(), expected);
    }

    #[rstest]
    fn test_permutation_between() {
        let tau = Permutation::between(&['a', 'b', 'c'], &['c', 'a', 'b']);
        assert_eq!(tau, Some(Permutation::from_image(&[2, 0, 1])));
    }

    #[rstest]
    #[case::length(&['a', 'b'], &['a'])]
    #[case::source_repetition(&['a', 'a'], &['a', 'b'])]
    #[case::target_repetition(&['a', 'b'], &['a', 'a'])]
    #[case::ambiguous_repetition(&['a', 'a'], &['a', 'a'])]
    #[case::membership(&['a', 'b'], &['a', 'c'])]
    fn test_permutation_between_error(#[case] from: &[char], #[case] to: &[char]) {
        assert_eq!(Permutation::between(from, to), None);
    }

    #[rstest]
    #[should_panic]
    fn test_permutation_between_degree_error() {
        Permutation::between(&[0, 1, 2, 3, 4, 5, 6], &[0, 1, 2, 3, 4, 5, 6]);
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

    #[rstest]
    #[case::single_cycle(3, vec![vec![0, 1, 2]], Permutation::from_image(&[1, 2, 0]))]
    #[case::double_transposition(4, vec![vec![0, 1], vec![2, 3]], Permutation::from_image(&[1, 0, 3, 2]))]
    #[case::partial(4, vec![vec![0, 1, 2]], Permutation::from_image(&[1, 2, 0, 3]))]
    #[case::empty(4, vec![], Permutation::identity(4))]
    #[case::empty_cycle(4, vec![vec![]], Permutation::identity(4))]
    fn test_permutation_from_cycles(
        #[case] degree: usize,
        #[case] cycles: Vec<Vec<usize>>,
        #[case] expected: Permutation,
    ) {
        assert_eq!(Permutation::from_cycles(degree, &cycles), Ok(expected));
    }

    #[rstest]
    #[case::overlap(
        3,
        vec![vec![0, 1], vec![1, 2]],
        PermutationError::DuplicateCyclePoint { point: 1 },
    )]
    #[case::repeated(
        3,
        vec![vec![0, 1, 0]],
        PermutationError::DuplicateCyclePoint { point: 0 },
    )]
    #[case::out_of_range(
        3,
        vec![vec![0, 3]],
        PermutationError::CyclePointOutOfRange {
            cycle: 0,
            position: 1,
            point: 3,
            degree: 3,
        },
    )]
    fn test_permutation_from_cycles_error(
        #[case] degree: usize,
        #[case] cycles: Vec<Vec<usize>>,
        #[case] expected: PermutationError,
    ) {
        assert_eq!(Permutation::from_cycles(degree, &cycles), Err(expected));
    }

    #[rstest]
    #[case::identity(Permutation::identity(4), vec![])]
    #[case::single_cycle(Permutation::from_image(&[1, 2, 0, 3]), vec![vec![0, 1, 2]])]
    #[case::double_transposition(Permutation::from_image(&[1, 0, 3, 2]), vec![vec![0, 1], vec![2, 3]])]
    #[case::interleaved(Permutation::from_image(&[2, 3, 0, 1]), vec![vec![0, 2], vec![1, 3]])]
    fn test_permutation_cycles(#[case] p: Permutation, #[case] expected: Vec<Vec<usize>>) {
        assert_eq!(p.cycles(), expected);
    }

    #[rstest]
    #[case::valid(
        &[2, 0, 1],
        Ok(Permutation { image: [2, 0, 1, 3, 4, 5], degree: 3 }),
    )]
    #[case::identity(
        &[],
        Ok(Permutation { image: [0, 1, 2, 3, 4, 5], degree: 0 }),
    )]
    #[case::too_long(
        &[0, 1, 2, 3, 4, 5, 6],
        Err(PermutationError::ImageTooLong { length: 7, maximum: MAX_DEGREE }),
    )]
    #[case::out_of_range(
        &[0, 1, 3],
        Err(PermutationError::ImageValueOutOfRange { position: 2, value: 3, degree: 3 }),
    )]
    #[case::duplicate(
        &[0, 1, 1],
        Err(PermutationError::DuplicateImageValue { value: 1 }),
    )]
    fn test_permutation_try_from(
        #[case] image: &[usize],
        #[case] expected: Result<Permutation, PermutationError>,
    ) {
        assert_eq!(Permutation::try_from(image), expected);
    }

    #[rstest]
    #[case::identity(Permutation::identity(3), "()")]
    #[case::three_cycle(Permutation::from_image(&[1, 2, 0]), "(0,1,2)")]
    #[case::double_transposition(Permutation::from_image(&[1, 0, 3, 2]), "(0,1)(2,3)")]
    fn test_permutation_display(#[case] p: Permutation, #[case] expected: &str) {
        assert_eq!(p.to_string(), expected);
    }
}
