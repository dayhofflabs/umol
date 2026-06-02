//! The coset space Sₙ/R: orderings partitioned by a proper-rotation group R.
//!
//! Equivalence is by the left coset σR = `{ σ ∘ r : r ∈ R }`; the canonical
//! representative of a coset is its min-rank element. This is the algebraic
//! layer — which orderings are equivalent — not yet a numbering (that is the
//! OpenSMILES index, added on top).

use crate::group::PermutationGroup;
use crate::permutation::Permutation;

/// A coset space Sₙ/R for a proper-rotation group R.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CosetSpace {
    group: PermutationGroup,
}

impl CosetSpace {
    pub fn new(group: PermutationGroup) -> Self {
        Self { group }
    }

    pub fn degree(&self) -> usize {
        self.group.degree()
    }

    pub fn group(&self) -> &PermutationGroup {
        &self.group
    }

    /// The number of cosets, `n! / |R|`.
    pub fn count(&self) -> usize {
        let total: usize = (1..=self.group.degree()).product();
        total / self.group.order()
    }

    /// The canonical representative of σ's coset: the min-rank element of σR.
    pub fn coset_rep(&self, perm: Permutation) -> Permutation {
        self.group
            .elements()
            .iter()
            .map(|&r| perm.compose(r))
            .min()
            .expect("R contains the identity")
    }

    /// The canonical representatives of all cosets, sorted (length = `count`).
    pub fn coset_reps(&self) -> Vec<Permutation> {
        let degree = self.degree();
        let total: usize = (1..=degree).product();
        let mut reps: Vec<Permutation> = (0..total)
            .map(|i| self.coset_rep(Permutation::unrank(degree, i)))
            .collect();
        reps.sort();
        reps.dedup();
        reps
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::tetrahedral(PermutationGroup::alternating(4), 2)]
    #[case::square_planar(PermutationGroup::dihedral(4), 3)]
    #[case::cyclic(PermutationGroup::cyclic(4), 6)]
    #[case::whole_group(PermutationGroup::symmetric(3), 1)]
    fn test_coset_space_count(#[case] group: PermutationGroup, #[case] count: usize) {
        let space = CosetSpace::new(group);
        assert_eq!(space.count(), count);
        assert_eq!(space.coset_reps().len(), count);
    }

    #[rstest]
    fn test_coset_space_coset_rep() {
        let space = CosetSpace::new(PermutationGroup::alternating(4));
        let identity = Permutation::identity(4);
        let even = Permutation::from_image(4, &[1, 2, 0, 3]);
        let odd = Permutation::from_image(4, &[1, 0, 2, 3]);
        assert_eq!(space.coset_rep(identity), space.coset_rep(even));
        assert_ne!(space.coset_rep(identity), space.coset_rep(odd));
    }
}
