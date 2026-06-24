//! Unfolded binary feature set: a sorted, duplicate-free vector of identifiers.

use std::cmp::Ordering;

/// A set of feature identifiers, sorted and duplicate-free. `Id` is the hash
/// width (`u32` / `u64` / `u128`) for hashed fingerprints, or a canonical
/// structural key (`Vec<u8>`) for the unhashed substructure fingerprint.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureSet<Id> {
    ids: Vec<Id>,
}

impl<Id: Clone + Ord> FeatureSet<Id> {
    /// Wrap an already sorted, duplicate-free identifier vector (e.g.
    /// [`umol_graph_core::Refinement::features`]) — the zero-copy fast path.
    pub fn from_sorted_unique(ids: Vec<Id>) -> Self {
        debug_assert!(
            ids.windows(2).all(|w| w[0] < w[1]),
            "ids must be sorted and duplicate-free"
        );
        Self { ids }
    }

    /// Build from arbitrary identifiers, sorting and de-duplicating in place.
    pub fn from_features(ids: impl IntoIterator<Item = Id>) -> Self {
        let mut ids: Vec<Id> = ids.into_iter().collect();
        ids.sort_unstable();
        ids.dedup();
        Self { ids }
    }

    pub fn ids(&self) -> &[Id] {
        &self.ids
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Identifiers present in both sets, by linear merge of the sorted runs.
    fn intersection_size(&self, other: &Self) -> usize {
        let (mut i, mut j, mut shared) = (0, 0, 0);
        while i < self.ids.len() && j < other.ids.len() {
            match self.ids[i].cmp(&other.ids[j]) {
                Ordering::Less => i += 1,
                Ordering::Greater => j += 1,
                Ordering::Equal => {
                    shared += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
        shared
    }

    /// Jaccard / Tanimoto similarity; two empty sets give `0.0`.
    pub fn tanimoto(&self, other: &Self) -> f64 {
        let shared = self.intersection_size(other);
        let union = self.ids.len() + other.ids.len() - shared;
        if union == 0 {
            0.0
        } else {
            shared as f64 / union as f64
        }
    }

    /// Sørensen–Dice similarity; two empty sets give `0.0`.
    pub fn dice(&self, other: &Self) -> f64 {
        let total = self.ids.len() + other.ids.len();
        if total == 0 {
            0.0
        } else {
            2.0 * self.intersection_size(other) as f64 / total as f64
        }
    }

    /// Every id in `self` occurs in `other` — `query.is_subset(target)` for
    /// prescreening. Binary-searches each id into the shrinking tail of `other`:
    /// O(|self|·log|other|), which beats a linear merge when the query is much
    /// smaller than the target.
    pub fn is_subset(&self, other: &Self) -> bool {
        let mut tail = other.ids.as_slice();
        for id in &self.ids {
            match tail.binary_search(id) {
                Ok(pos) => tail = &tail[pos + 1..],
                Err(_) => return false,
            }
        }
        true
    }
}

/// A counted multiset of feature identifiers: sorted by `Id`, each with its
/// occurrence count. The counted counterpart of [`FeatureSet`] — the count of an
/// identifier is its multiplicity in the featurizer's pre-dedup output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountedFeatureSet<Id> {
    entries: Vec<(Id, u32)>,
}

impl<Id: Clone + Ord> CountedFeatureSet<Id> {
    /// Count occurrences of each identifier, sorted by identifier.
    pub fn from_features(ids: impl IntoIterator<Item = Id>) -> Self {
        let mut ids: Vec<Id> = ids.into_iter().collect();
        ids.sort_unstable();
        let mut entries: Vec<(Id, u32)> = Vec::new();
        for id in ids {
            match entries.last_mut() {
                Some((last, count)) if *last == id => *count += 1,
                _ => entries.push((id, 1)),
            }
        }
        Self { entries }
    }

    /// Build from `(identifier, count)` pairs with unique identifiers (e.g.
    /// [`umol_graph_core::Refinement::counts`]), sorting by identifier.
    pub fn from_counts(counts: impl IntoIterator<Item = (Id, u32)>) -> Self {
        let mut entries: Vec<(Id, u32)> = counts.into_iter().collect();
        entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        debug_assert!(
            entries.windows(2).all(|w| w[0].0 < w[1].0),
            "counts must have unique identifiers"
        );
        Self { entries }
    }

    /// The `(identifier, count)` pairs, sorted by identifier.
    pub fn entries(&self) -> &[(Id, u32)] {
        &self.entries
    }

    /// Count of `id`, or 0 if absent.
    pub fn count(&self, id: &Id) -> u32 {
        self.entries
            .binary_search_by(|(entry, _)| entry.cmp(id))
            .map(|i| self.entries[i].1)
            .unwrap_or(0)
    }

    /// Number of distinct identifiers.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// A signed difference of two counted feature sets: sorted by `Id`, each with a
/// nonzero signed count. Used for the reaction difference fingerprint (product
/// counts minus reactant counts).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SignedFeatureSet<Id> {
    entries: Vec<(Id, i32)>,
}

impl<Id: Clone + Ord> SignedFeatureSet<Id> {
    /// `minuend − subtrahend` per identifier (signed), dropping zeros. Linear merge
    /// of the two sorted counted sets.
    pub fn difference(minuend: &CountedFeatureSet<Id>, subtrahend: &CountedFeatureSet<Id>) -> Self {
        let (a, b) = (minuend.entries(), subtrahend.entries());
        let mut entries: Vec<(Id, i32)> = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].0.cmp(&b[j].0) {
                Ordering::Less => {
                    entries.push((a[i].0.clone(), a[i].1 as i32));
                    i += 1;
                }
                Ordering::Greater => {
                    entries.push((b[j].0.clone(), -(b[j].1 as i32)));
                    j += 1;
                }
                Ordering::Equal => {
                    let delta = a[i].1 as i32 - b[j].1 as i32;
                    if delta != 0 {
                        entries.push((a[i].0.clone(), delta));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        entries.extend(a[i..].iter().map(|(id, c)| (id.clone(), *c as i32)));
        entries.extend(b[j..].iter().map(|(id, c)| (id.clone(), -(*c as i32))));
        Self { entries }
    }

    /// The `(identifier, signed count)` pairs, sorted by identifier.
    pub fn entries(&self) -> &[(Id, i32)] {
        &self.entries
    }

    /// Signed count of `id`, or 0 if absent.
    pub fn count(&self, id: &Id) -> i32 {
        self.entries
            .binary_search_by(|(entry, _)| entry.cmp(id))
            .map(|i| self.entries[i].1)
            .unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn set(ids: &[u64]) -> FeatureSet<u64> {
        FeatureSet::from_features(ids.iter().copied())
    }

    #[rstest]
    #[case::sorts_and_dedups(vec![3u64, 1, 2, 3, 1], vec![1, 2, 3])]
    #[case::already_sorted(vec![1u64, 2, 3], vec![1, 2, 3])]
    #[case::empty(vec![], vec![])]
    fn test_feature_set_from_features(#[case] input: Vec<u64>, #[case] expected: Vec<u64>) {
        assert_eq!(FeatureSet::from_features(input).ids(), expected.as_slice());
    }

    #[rstest]
    #[case(vec![1u64, 5, 9])]
    fn test_feature_set_from_sorted_unique(#[case] ids: Vec<u64>) {
        assert_eq!(FeatureSet::from_sorted_unique(ids.clone()).ids(), ids.as_slice());
    }

    #[rstest]
    #[case::partial(vec![1, 2, 3], vec![2, 3, 4], 0.5)]
    #[case::identical(vec![1, 2, 3], vec![1, 2, 3], 1.0)]
    #[case::disjoint(vec![1, 2], vec![3, 4], 0.0)]
    #[case::both_empty(vec![], vec![], 0.0)]
    fn test_feature_set_tanimoto(#[case] a: Vec<u64>, #[case] b: Vec<u64>, #[case] expected: f64) {
        assert!((set(&a).tanimoto(&set(&b)) - expected).abs() < 1e-12);
    }

    #[rstest]
    #[case::partial(vec![1, 2, 3], vec![2, 3, 4], 4.0 / 6.0)]
    #[case::identical(vec![1, 2, 3], vec![1, 2, 3], 1.0)]
    #[case::disjoint(vec![1, 2], vec![3, 4], 0.0)]
    #[case::both_empty(vec![], vec![], 0.0)]
    fn test_feature_set_dice(#[case] a: Vec<u64>, #[case] b: Vec<u64>, #[case] expected: f64) {
        assert!((set(&a).dice(&set(&b)) - expected).abs() < 1e-12);
    }

    #[rstest]
    #[case::proper_subset(vec![1, 2], vec![1, 2, 3], true)]
    #[case::equal(vec![1, 2, 3], vec![1, 2, 3], true)]
    #[case::empty_query(vec![], vec![1, 2], true)]
    #[case::missing_id(vec![1, 4], vec![1, 2, 3], false)]
    #[case::superset(vec![1, 2, 3], vec![1, 2], false)]
    fn test_feature_set_is_subset(
        #[case] query: Vec<u64>,
        #[case] target: Vec<u64>,
        #[case] expected: bool,
    ) {
        assert_eq!(set(&query).is_subset(&set(&target)), expected);
    }

    #[rstest]
    #[case::counts_and_sorts(vec![3u64, 1, 2, 3, 1], vec![(1, 2), (2, 1), (3, 2)])]
    #[case::all_distinct(vec![1u64, 2, 3], vec![(1, 1), (2, 1), (3, 1)])]
    #[case::empty(vec![], vec![])]
    fn test_counted_feature_set_from_features(
        #[case] input: Vec<u64>,
        #[case] expected: Vec<(u64, u32)>,
    ) {
        assert_eq!(
            CountedFeatureSet::from_features(input).entries(),
            expected.as_slice()
        );
    }

    #[rstest]
    fn test_signed_feature_set_difference() {
        let minuend = CountedFeatureSet::from_features(vec![1u64, 1, 2, 3]);
        let subtrahend = CountedFeatureSet::from_features(vec![2u64, 3, 3, 4]);
        // 1: 2−0, 2: 1−1=0 (dropped), 3: 1−2=−1, 4: 0−1=−1.
        assert_eq!(
            SignedFeatureSet::difference(&minuend, &subtrahend).entries(),
            &[(1u64, 2), (3, -1), (4, -1)]
        );
    }
}
