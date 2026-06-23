//! Folded bit fingerprint: a fixed-width bit array over `width` buckets.

use bitvec::order::Lsb0;
use bitvec::vec::BitVec;

use super::feature_set::FeatureSet;

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

    pub fn get(&self, bit: usize) -> bool {
        self.bits[bit]
    }

    pub fn count_ones(&self) -> usize {
        self.bits.count_ones()
    }

    /// Tanimoto over the bit sets; two empty fingerprints give `0.0`.
    pub fn tanimoto(&self, other: &Self) -> f64 {
        assert_eq!(self.bits.len(), other.bits.len(), "BitFp width mismatch");
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        let intersection: u32 = a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum();
        let union: u32 = a.iter().zip(b).map(|(x, y)| (x | y).count_ones()).sum();
        if union == 0 {
            0.0
        } else {
            f64::from(intersection) / f64::from(union)
        }
    }

    /// Sørensen–Dice over the bit sets; two empty fingerprints give `0.0`.
    pub fn dice(&self, other: &Self) -> f64 {
        assert_eq!(self.bits.len(), other.bits.len(), "BitFp width mismatch");
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        let intersection: u32 = a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum();
        let total = self.count_ones() + other.count_ones();
        if total == 0 {
            0.0
        } else {
            2.0 * f64::from(intersection) / total as f64
        }
    }

    /// Every set bit of `self` is set in `other` — `query.is_subset(target)`.
    pub fn is_subset(&self, other: &Self) -> bool {
        assert_eq!(self.bits.len(), other.bits.len(), "BitFp width mismatch");
        let (a, b) = (self.bits.as_raw_slice(), other.bits.as_raw_slice());
        a.iter().zip(b).all(|(x, y)| x & !y == 0)
    }
}

impl FeatureSet<u64> {
    /// Fold to a fixed-width [`BitFp`]: bit `id % width` set for each identifier.
    pub fn fold(&self, width: usize) -> BitFp {
        assert!(width > 0, "fold width must be positive");
        let mut bits = BitFp::zeros(width);
        for &id in self.ids() {
            bits.set((id % width as u64) as usize);
        }
        bits
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn bits(width: usize, positions: &[usize]) -> BitFp {
        let mut fp = BitFp::zeros(width);
        for &position in positions {
            fp.set(position);
        }
        fp
    }

    #[rstest]
    #[case::collision_and_distinct(vec![1u64, 5, 2049, 2], 2048, vec![1, 2, 5])]
    #[case::no_collision(vec![0u64, 63, 64], 128, vec![0, 63, 64])]
    fn test_feature_set_fold(
        #[case] ids: Vec<u64>,
        #[case] width: usize,
        #[case] set_bits: Vec<usize>,
    ) {
        let folded = FeatureSet::from_features(ids).fold(width);
        assert_eq!(folded.width(), width);
        assert_eq!(folded.count_ones(), set_bits.len());
        for bit in set_bits {
            assert!(folded.get(bit));
        }
    }

    #[rstest]
    #[case::partial(&[1, 2, 3], &[2, 3, 4], 0.5)]
    #[case::identical(&[1, 2, 3], &[1, 2, 3], 1.0)]
    #[case::disjoint(&[1, 2], &[3, 4], 0.0)]
    #[case::both_empty(&[], &[], 0.0)]
    fn test_bit_fp_tanimoto(#[case] a: &[usize], #[case] b: &[usize], #[case] expected: f64) {
        assert!((bits(8, a).tanimoto(&bits(8, b)) - expected).abs() < 1e-12);
    }

    #[rstest]
    #[case::partial(&[1, 2, 3], &[2, 3, 4], 4.0 / 6.0)]
    #[case::identical(&[1, 2, 3], &[1, 2, 3], 1.0)]
    #[case::both_empty(&[], &[], 0.0)]
    fn test_bit_fp_dice(#[case] a: &[usize], #[case] b: &[usize], #[case] expected: f64) {
        assert!((bits(8, a).dice(&bits(8, b)) - expected).abs() < 1e-12);
    }

    #[rstest]
    #[case::proper_subset(&[1, 2], &[1, 2, 3], true)]
    #[case::equal(&[1, 2, 3], &[1, 2, 3], true)]
    #[case::empty_query(&[], &[1, 2], true)]
    #[case::missing_bit(&[1, 4], &[1, 2, 3], false)]
    fn test_bit_fp_is_subset(#[case] query: &[usize], #[case] target: &[usize], #[case] expected: bool) {
        assert_eq!(bits(8, query).is_subset(&bits(8, target)), expected);
    }
}
