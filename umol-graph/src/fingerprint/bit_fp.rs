//! Folded bit fingerprint: a fixed-width bit array over `width` buckets.

use bitvec::order::Lsb0;
use bitvec::vec::BitVec;

use super::feature_set::FeatureSet;
use super::featurizer::FingerprintError;

/// A fixed-width bit fingerprint: bit `i` (for `i` in `0..width`) is a bucket.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BitFp {
    bits: BitVec<u64, Lsb0>,
}

impl BitFp {
    fn zeros(width: usize) -> Self {
        Self {
            bits: BitVec::repeat(false, width),
        }
    }

    fn set(&mut self, bit: usize) {
        self.bits.set(bit, true);
    }

    pub fn width(&self) -> usize {
        self.bits.len()
    }

    /// Return whether `bit` is set, or `None` when it is outside the fingerprint.
    pub fn get(&self, bit: usize) -> Option<bool> {
        self.bits.get(bit).map(|value| *value)
    }

    pub fn count_ones(&self) -> usize {
        self.bits.count_ones()
    }

    /// Tanimoto over equal-width bit sets; two empty fingerprints give `0.0`.
    ///
    /// Returns [`FingerprintError::WidthMismatch`] when the widths differ.
    pub fn tanimoto(&self, other: &Self) -> Result<f64, FingerprintError> {
        if self.bits.len() != other.bits.len() {
            return Err(FingerprintError::WidthMismatch {
                left: self.bits.len(),
                right: other.bits.len(),
            });
        }
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        let intersection: u32 = a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum();
        let union: u32 = a.iter().zip(b).map(|(x, y)| (x | y).count_ones()).sum();
        if union == 0 {
            Ok(0.0)
        } else {
            Ok(f64::from(intersection) / f64::from(union))
        }
    }

    /// Sørensen–Dice over equal-width bit sets; two empty fingerprints give `0.0`.
    ///
    /// Returns [`FingerprintError::WidthMismatch`] when the widths differ.
    pub fn dice(&self, other: &Self) -> Result<f64, FingerprintError> {
        if self.bits.len() != other.bits.len() {
            return Err(FingerprintError::WidthMismatch {
                left: self.bits.len(),
                right: other.bits.len(),
            });
        }
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        let intersection: u32 = a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum();
        let total = self.count_ones() + other.count_ones();
        if total == 0 {
            Ok(0.0)
        } else {
            Ok(2.0 * f64::from(intersection) / total as f64)
        }
    }

    /// Every set bit of `self` is set in `other` — `query.is_subset(target)`.
    ///
    /// Returns [`FingerprintError::WidthMismatch`] when the widths differ.
    pub fn is_subset(&self, other: &Self) -> Result<bool, FingerprintError> {
        if self.bits.len() != other.bits.len() {
            return Err(FingerprintError::WidthMismatch {
                left: self.bits.len(),
                right: other.bits.len(),
            });
        }
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        Ok(a.iter().zip(b).all(|(x, y)| x & !y == 0))
    }
}

impl<Id> FeatureSet<Id>
where
    Id: Clone + Copy + Ord,
    u128: From<Id>,
{
    /// Fold to a fixed-width [`BitFp`]: bit `id % width` set for each identifier.
    ///
    /// Returns [`FingerprintError::ZeroWidth`] when `width` is zero.
    pub fn fold(&self, width: usize) -> Result<BitFp, FingerprintError> {
        if width == 0 {
            return Err(FingerprintError::ZeroWidth);
        }
        let mut bits = BitFp::zeros(width);
        for &id in self.ids() {
            bits.set((u128::from(id) % width as u128) as usize);
        }
        Ok(bits)
    }
}

#[cfg(test)]
mod tests {
    use bitvec::bitvec;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::populated(vec![1u64, 7], 8, 8)]
    #[case::empty(vec![], 4, 4)]
    fn test_bit_fp_width(#[case] ids: Vec<u64>, #[case] width: usize, #[case] expected: usize) {
        assert_eq!(
            FeatureSet::from_features(ids).fold(width).unwrap().width(),
            expected
        );
    }

    #[rstest]
    #[case::last_valid(vec![7u64], 8, 7, Some(true))]
    #[case::unset(vec![7u64], 8, 6, Some(false))]
    #[case::first_invalid(vec![7u64], 8, 8, None)]
    #[case::empty(vec![], 8, 0, Some(false))]
    fn test_bit_fp_get(
        #[case] ids: Vec<u64>,
        #[case] width: usize,
        #[case] bit: usize,
        #[case] expected: Option<bool>,
    ) {
        assert_eq!(
            FeatureSet::from_features(ids).fold(width).unwrap().get(bit),
            expected
        );
    }

    #[rstest]
    #[case::populated(vec![1u64, 2, 5, 9], 8, 3)]
    #[case::empty(vec![], 8, 0)]
    fn test_bit_fp_count_ones(
        #[case] ids: Vec<u64>,
        #[case] width: usize,
        #[case] expected: usize,
    ) {
        assert_eq!(
            FeatureSet::from_features(ids)
                .fold(width)
                .unwrap()
                .count_ones(),
            expected
        );
    }

    #[rstest]
    #[case::partial(vec![1u64, 2, 3], vec![2, 3, 4], 0.5)]
    #[case::identical(vec![1u64, 2, 3], vec![1, 2, 3], 1.0)]
    #[case::disjoint(vec![1u64, 2], vec![3, 4], 0.0)]
    #[case::both_empty(vec![], vec![], 0.0)]
    fn test_bit_fp_tanimoto(#[case] a: Vec<u64>, #[case] b: Vec<u64>, #[case] expected: f64) {
        let a = FeatureSet::from_features(a).fold(8).unwrap();
        let b = FeatureSet::from_features(b).fold(8).unwrap();
        assert_eq!(a.tanimoto(&b), Ok(expected));
    }

    #[rstest]
    #[case::unequal_widths(vec![1u64], 8, vec![1u64], 4)]
    fn test_bit_fp_tanimoto_error(
        #[case] a: Vec<u64>,
        #[case] a_width: usize,
        #[case] b: Vec<u64>,
        #[case] b_width: usize,
    ) {
        let a = FeatureSet::from_features(a).fold(a_width).unwrap();
        let b = FeatureSet::from_features(b).fold(b_width).unwrap();
        assert_eq!(
            a.tanimoto(&b),
            Err(FingerprintError::WidthMismatch {
                left: a_width,
                right: b_width,
            })
        );
    }

    #[rstest]
    #[case::partial(vec![1u64, 2, 3], vec![2, 3, 4], 4.0 / 6.0)]
    #[case::identical(vec![1u64, 2, 3], vec![1, 2, 3], 1.0)]
    #[case::both_empty(vec![], vec![], 0.0)]
    fn test_bit_fp_dice(#[case] a: Vec<u64>, #[case] b: Vec<u64>, #[case] expected: f64) {
        let a = FeatureSet::from_features(a).fold(8).unwrap();
        let b = FeatureSet::from_features(b).fold(8).unwrap();
        assert_eq!(a.dice(&b), Ok(expected));
    }

    #[rstest]
    #[case::unequal_widths(vec![1u64], 8, vec![1u64], 4)]
    fn test_bit_fp_dice_error(
        #[case] a: Vec<u64>,
        #[case] a_width: usize,
        #[case] b: Vec<u64>,
        #[case] b_width: usize,
    ) {
        let a = FeatureSet::from_features(a).fold(a_width).unwrap();
        let b = FeatureSet::from_features(b).fold(b_width).unwrap();
        assert_eq!(
            a.dice(&b),
            Err(FingerprintError::WidthMismatch {
                left: a_width,
                right: b_width,
            })
        );
    }

    #[rstest]
    #[case::proper_subset(vec![1u64, 2], vec![1u64, 2, 3], true)]
    #[case::equal(vec![1u64, 2, 3], vec![1u64, 2, 3], true)]
    #[case::empty_query(vec![], vec![1u64, 2], true)]
    #[case::missing_bit(vec![1u64, 4], vec![1u64, 2, 3], false)]
    fn test_bit_fp_is_subset(
        #[case] query: Vec<u64>,
        #[case] target: Vec<u64>,
        #[case] expected: bool,
    ) {
        let query = FeatureSet::from_features(query).fold(8).unwrap();
        let target = FeatureSet::from_features(target).fold(8).unwrap();
        assert_eq!(query.is_subset(&target), Ok(expected));
    }

    #[rstest]
    #[case::unequal_widths(vec![1u64], 8, vec![1u64], 4)]
    fn test_bit_fp_is_subset_error(
        #[case] query: Vec<u64>,
        #[case] query_width: usize,
        #[case] target: Vec<u64>,
        #[case] target_width: usize,
    ) {
        let query = FeatureSet::from_features(query).fold(query_width).unwrap();
        let target = FeatureSet::from_features(target)
            .fold(target_width)
            .unwrap();
        assert_eq!(
            query.is_subset(&target),
            Err(FingerprintError::WidthMismatch {
                left: query_width,
                right: target_width,
            })
        );
    }

    #[rstest]
    #[case::width32_collision(
        FeatureSet::from_features([1u32, 2, 5, 9]).fold(8),
        BitFp {
            bits: bitvec![u64, Lsb0; 0, 1, 1, 0, 0, 1, 0, 0],
        }
    )]
    #[case::width64_empty(
        FeatureSet::from_features(Vec::<u64>::new()).fold(4),
        BitFp {
            bits: bitvec![u64, Lsb0; 0; 4],
        }
    )]
    #[case::width128_high_identifier(
        FeatureSet::from_features([(u64::MAX as u128) + 2]).fold(8),
        BitFp {
            bits: bitvec![u64, Lsb0; 0, 1, 0, 0, 0, 0, 0, 0],
        }
    )]
    fn test_feature_set_fold(
        #[case] actual: Result<BitFp, FingerprintError>,
        #[case] expected: BitFp,
    ) {
        assert_eq!(actual, Ok(expected));
    }

    #[rstest]
    #[case::width32(FeatureSet::from_features([1u32, 2]).fold(0))]
    #[case::width64(FeatureSet::from_features([1u64, 2]).fold(0))]
    #[case::width128(FeatureSet::from_features([1u128, 2]).fold(0))]
    fn test_feature_set_fold_error(#[case] actual: Result<BitFp, FingerprintError>) {
        assert_eq!(actual, Err(FingerprintError::ZeroWidth));
    }
}
