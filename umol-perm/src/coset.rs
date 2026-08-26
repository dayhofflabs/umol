//! The coset space R\P for a parent group P (Sₙ unless restricted), with a per-class numbering.
//!
//! A configuration is a ligand→position map; two configurations are the same
//! arrangement when related by a proper rotation acting on the positions, i.e.
//! by left multiplication `r∘σ`. The class is therefore the right coset Rσ, and
//! the canonical representative is its min-rank element. On top of that algebra
//! sits a `Decomposition` that numbers the cosets — `CanonicalRank` (the generic
//! Lehmer-min ordering, which is the parity bit for a 2-coset space) or a
//! geometry that reproduces the OpenSMILES arrangement number.

use std::collections::HashMap;

use crate::group::PermutationGroup;
use crate::permutation::Permutation;

/// How a space numbers its cosets.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Decomposition {
    CanonicalRank,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

/// A coset space P/R for a proper-rotation group R inside a parent group P, with a fixed
/// numbering. P is the group of realizable arrangements: Sₙ for the geometry classes (any
/// ligand may take any position), or a partition subgroup for cis/trans (substituents are
/// bonded to fixed sp² carbons).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CosetSpace {
    parent: PermutationGroup,
    is_partitioned: bool,
    group: PermutationGroup,
    numbering: HashMap<Permutation, u32>,
    representatives: Vec<Permutation>,
    improper: Permutation,
}

impl CosetSpace {
    pub(crate) fn new(
        parent: PermutationGroup,
        is_partitioned: bool,
        group: PermutationGroup,
        decomposition: Decomposition,
        improper: Permutation,
    ) -> Self {
        assert!(
            group.elements().iter().all(|g| parent.contains(g)),
            "coset group R is not a subgroup of the parent P"
        );
        let ordered = decomposition.representatives(&parent, &group);
        let mut numbering = HashMap::new();
        let mut representatives = Vec::with_capacity(ordered.len());
        for rep in ordered {
            let canonical = coset_rep(&group, rep);
            let previous = numbering.insert(canonical, representatives.len() as u32);
            assert!(
                previous.is_none(),
                "decomposition assigns one coset two numbers"
            );
            representatives.push(rep);
        }
        assert_eq!(
            numbering.len(),
            parent.order() / group.order(),
            "decomposition does not cover every coset"
        );
        Self {
            parent,
            is_partitioned,
            group,
            numbering,
            representatives,
            improper,
        }
    }

    pub fn degree(&self) -> usize {
        self.group.degree()
    }

    pub fn group(&self) -> &PermutationGroup {
        &self.group
    }

    /// The orientation-reversing generator (the improper/mirror operation). The
    /// identity for achiral classes.
    pub fn improper(&self) -> Permutation {
        self.improper
    }

    /// The number of cosets, `n! / |R|`.
    pub fn count(&self) -> usize {
        self.representatives.len()
    }

    /// The canonical representative of σ's coset (the min-rank element of Rσ), or
    /// `None` if σ's degree does not match this space.
    pub fn coset_rep(&self, permutation: Permutation) -> Option<Permutation> {
        (permutation.degree() == self.degree()).then(|| coset_rep(&self.group, permutation))
    }

    /// The arrangement number of σ's coset, or `None` if σ is not in the parent
    /// group (wrong degree, or a coset this space does not number).
    pub fn index(&self, permutation: Permutation) -> Option<u32> {
        self.numbering.get(&self.coset_rep(permutation)?).copied()
    }

    /// A representative permutation for arrangement number `index`, or `None` if
    /// `index >= count`.
    pub fn unindex(&self, index: u32) -> Option<Permutation> {
        self.representatives.get(index as usize).copied()
    }

    /// The arrangement number of configuration `index` after the neighbor list is
    /// relabeled by `permutation` (`permutation(i)` is the original-list position of
    /// the relabeled list's i-th neighbor), or `None` if `index >= count` or
    /// `permutation` is not in the parent group. Carries a parsed `@`-number into
    /// umol's incidence order.
    pub fn reindex(&self, index: u32, permutation: Permutation) -> Option<u32> {
        if !self.parent.contains(&permutation) {
            return None;
        }
        self.index(self.unindex(index)?.compose(permutation))
    }

    /// The enantiomeric coset, or `None` if `index >= count`.
    pub fn enantiomer(&self, index: u32) -> Option<u32> {
        self.index(self.unindex(index)?.compose(self.improper))
    }

    /// Chiral iff the improper generator moves some coset.
    pub fn is_chiral(&self) -> bool {
        (0..self.count() as u32).any(|i| self.enantiomer(i).expect("0..count is in range") != i)
    }

    /// Whether `permutation` is an allowed action in this space's parent group.
    ///
    /// A permutation of another degree is not allowed.
    pub fn allows(&self, permutation: Permutation) -> bool {
        self.parent.contains(&permutation)
    }

    /// An allowed permutation taking `frame` to its lexicographically least presentation.
    ///
    /// Returns `None` when the frame length differs from the space degree. Unrestricted parent
    /// groups sort the complete frame. The restricted four-position parent sorts each two-position
    /// block and then orders the blocks.
    ///
    /// # Semantic properties
    ///
    /// The returned permutation is accepted by [`Self::allows`]. Applying the operation again to
    /// the resulting frame returns the identity permutation.
    pub fn normalizer<T: Ord>(&self, frame: &[T]) -> Option<Permutation> {
        if frame.len() != self.degree() {
            return None;
        }

        let mut image: Vec<usize> = (0..self.degree()).collect();
        if !self.is_partitioned {
            image.sort_by(|&left, &right| frame[left].cmp(&frame[right]).then(left.cmp(&right)));
        } else {
            debug_assert_eq!(self.degree(), 4);
            image[..2]
                .sort_by(|&left, &right| frame[left].cmp(&frame[right]).then(left.cmp(&right)));
            image[2..]
                .sort_by(|&left, &right| frame[left].cmp(&frame[right]).then(left.cmp(&right)));
            if [&frame[image[0]], &frame[image[1]]] > [&frame[image[2]], &frame[image[3]]] {
                image.swap(0, 2);
                image.swap(1, 3);
            }
        }
        Some(Permutation::from_image(&image))
    }

    /// The canonical (minimum) orbit representative per coset under the group
    /// generated by `generators` acting on the right, or `None` if a generator
    /// does not lie in the parent group.
    pub fn orbit_reps(&self, generators: &[Permutation]) -> Option<Vec<u32>> {
        if generators
            .iter()
            .any(|generator| !self.parent.contains(generator))
        {
            return None;
        }
        fn root(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        let count = self.count() as u32;
        let mut parent: Vec<u32> = (0..count).collect();
        for i in 0..count {
            let rep = self.unindex(i).expect("0..count is in range");
            for &generator in generators {
                let j = self.index(rep.compose(generator))?;
                let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                if ri != rj {
                    parent[ri.max(rj) as usize] = ri.min(rj);
                }
            }
        }
        Some((0..count).map(|i| root(&mut parent, i)).collect())
    }

    /// The observable coset under the supergroup generated by the proper-rotation group and
    /// `supergroup_generators`, or `None` if `index >= count` or a generator does not lie in the
    /// parent group.
    pub fn observable_coset(
        &self,
        index: u32,
        supergroup_generators: &[Permutation],
    ) -> Option<u32> {
        self.orbit_reps(supergroup_generators)?
            .get(index as usize)
            .copied()
    }
}

fn coset_rep(group: &PermutationGroup, permutation: Permutation) -> Permutation {
    group
        .elements()
        .iter()
        .map(|&r| r.compose(permutation))
        .min()
        .expect("R contains the identity")
}

impl Decomposition {
    /// One representative permutation per arrangement number, in numbering order.
    fn representatives(
        self,
        parent: &PermutationGroup,
        group: &PermutationGroup,
    ) -> Vec<Permutation> {
        match self {
            Decomposition::CanonicalRank => {
                let mut reps: Vec<Permutation> = parent
                    .elements()
                    .iter()
                    .map(|&permutation| coset_rep(group, permutation))
                    .collect();
                reps.sort();
                reps.dedup();
                reps
            }
            Decomposition::SquarePlanar => square_planar_reps(),
            Decomposition::TrigonalBipyramidal => trigonal_bipyramidal_reps(),
            Decomposition::Octahedral => octahedral_reps(),
        }
    }
}

/// `@SP1`/`@SP2`/`@SP3` = the U/4/Z path shapes (OpenSMILES §3.8.5); one
/// representative per shape, taken from the spec's enumeration (U = 1234,
/// 4 = 2413 with two diagonal steps, Z = 2314 with one).
fn square_planar_reps() -> Vec<Permutation> {
    vec![
        Permutation::from_image(&[0, 1, 2, 3]),
        Permutation::from_image(&[1, 3, 0, 2]),
        Permutation::from_image(&[1, 2, 0, 3]),
    ]
}

/// `@TB1`..`@TB20` (OpenSMILES §3.8.6): an ordered axial pair `(from, towards)`
/// at vertices `0,4`, with the three equatorial atoms anticlockwise (`@`) or
/// clockwise (`@@`) viewed `from → towards`.
fn trigonal_bipyramidal_reps() -> Vec<Permutation> {
    const AXES: [(u8, u8, bool); 20] = [
        (0, 4, true),
        (0, 4, false),
        (0, 3, true),
        (0, 3, false),
        (0, 2, true),
        (0, 2, false),
        (0, 1, true),
        (0, 1, false),
        (1, 4, true),
        (1, 3, true),
        (1, 4, false),
        (1, 3, false),
        (1, 2, true),
        (1, 2, false),
        (2, 4, true),
        (2, 3, true),
        (3, 4, true),
        (3, 4, false),
        (2, 3, false),
        (2, 4, false),
    ];
    AXES.iter()
        .map(|&(from, towards, anticlockwise)| {
            let mut image = [0usize; 5];
            image[from as usize] = 0;
            image[towards as usize] = 4;
            let plane: Vec<u8> = (0..5u8).filter(|&l| l != from && l != towards).collect();
            let vertices: [u8; 3] = if anticlockwise { [1, 2, 3] } else { [3, 2, 1] };
            for (&label, &vertex) in plane.iter().zip(vertices.iter()) {
                image[label as usize] = usize::from(vertex);
            }
            Permutation::from_image(&image)
        })
        .collect()
}

/// `@OH1`..`@OH30` (OpenSMILES §3.8.7): the axis runs from ligand `a` (vertex 0)
/// to `towards` (vertex 5); the remaining four ligands fill the equatorial
/// square (vertices 1–4) in one of six classes — each path shape (U with no
/// diagonal step, Z with one, 4 with two) in its two windings. The `@` and `@@`
/// windings are the two C₄ orbits of the shape (not an ordering and its
/// reverse, which can lie in the same orbit); the reps below are taken from the
/// spec's shape enumeration.
fn octahedral_reps() -> Vec<Permutation> {
    const U_ANTI: [u8; 4] = [1, 2, 3, 4];
    const U_CW: [u8; 4] = [4, 3, 2, 1];
    const FOUR_ANTI: [u8; 4] = [2, 4, 1, 3];
    const FOUR_CW: [u8; 4] = [1, 3, 2, 4];
    const Z_ANTI: [u8; 4] = [2, 3, 1, 4];
    const Z_CW: [u8; 4] = [2, 1, 3, 4];
    const ARRANGEMENTS: [(u8, [u8; 4]); 30] = [
        (5, U_ANTI),
        (5, U_CW),
        (4, U_ANTI),
        (5, Z_ANTI),
        (4, Z_ANTI),
        (3, U_ANTI),
        (3, Z_ANTI),
        (5, FOUR_CW),
        (4, FOUR_CW),
        (5, FOUR_ANTI),
        (4, FOUR_ANTI),
        (3, FOUR_CW),
        (3, FOUR_ANTI),
        (5, Z_CW),
        (4, Z_CW),
        (4, U_CW),
        (3, Z_CW),
        (3, U_CW),
        (2, U_ANTI),
        (2, Z_ANTI),
        (2, FOUR_CW),
        (2, FOUR_ANTI),
        (2, Z_CW),
        (2, U_CW),
        (1, U_ANTI),
        (1, Z_ANTI),
        (1, FOUR_CW),
        (1, FOUR_ANTI),
        (1, Z_CW),
        (1, U_CW),
    ];
    ARRANGEMENTS
        .iter()
        .map(|&(towards, equatorial)| {
            let mut image = [0usize; 6];
            image[towards as usize] = 5;
            let plane: Vec<u8> = (0..6u8).filter(|&l| l != 0 && l != towards).collect();
            for (&label, &vertex) in plane.iter().zip(equatorial.iter()) {
                image[label as usize] = usize::from(vertex);
            }
            Permutation::from_image(&image)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::class::ClassKey;

    #[rstest]
    #[case::tetrahedral(
        ClassKey::Tetrahedral,
        Permutation::from_image(&[1, 0, 2, 3]),
        true,
    )]
    #[case::cis_trans(
        ClassKey::CisTrans,
        Permutation::from_image(&[1, 0, 2, 3]),
        true,
    )]
    #[case::axial(
        ClassKey::Axial,
        Permutation::from_image(&[2, 3, 0, 1]),
        true,
    )]
    #[case::square_planar(
        ClassKey::SquarePlanar,
        Permutation::from_image(&[2, 0, 3, 1]),
        true,
    )]
    #[case::trigonal_bipyramidal(
        ClassKey::TrigonalBipyramidal,
        Permutation::from_image(&[4, 1, 2, 3, 0]),
        true,
    )]
    #[case::octahedral(
        ClassKey::Octahedral,
        Permutation::from_image(&[5, 1, 2, 3, 4, 0]),
        true,
    )]
    #[case::restricted_parent(
        ClassKey::CisTrans,
        Permutation::from_image(&[2, 1, 0, 3]),
        false,
    )]
    #[case::wrong_degree(ClassKey::Tetrahedral, Permutation::identity(3), false)]
    fn test_coset_space_allows(
        #[case] key: ClassKey,
        #[case] permutation: Permutation,
        #[case] expected: bool,
    ) {
        assert_eq!(key.space().allows(permutation), expected);
    }

    #[rstest]
    #[case::tetrahedral(
        ClassKey::Tetrahedral,
        vec![3, 1, 4, 2],
        Permutation::from_image(&[1, 3, 0, 2]),
    )]
    #[case::cis_trans(
        ClassKey::CisTrans,
        vec![4, 2, 3, 1],
        Permutation::from_image(&[3, 2, 1, 0]),
    )]
    #[case::axial(
        ClassKey::Axial,
        vec![3, 1, 4, 2],
        Permutation::from_image(&[1, 0, 3, 2]),
    )]
    #[case::square_planar(
        ClassKey::SquarePlanar,
        vec![3, 1, 4, 2],
        Permutation::from_image(&[1, 3, 0, 2]),
    )]
    #[case::trigonal_bipyramidal(
        ClassKey::TrigonalBipyramidal,
        vec![4, 0, 3, 1, 2],
        Permutation::from_image(&[1, 3, 4, 2, 0]),
    )]
    #[case::octahedral(
        ClassKey::Octahedral,
        vec![5, 3, 1, 4, 0, 2],
        Permutation::from_image(&[4, 2, 5, 1, 3, 0]),
    )]
    #[case::unrestricted_repeated(
        ClassKey::Tetrahedral,
        vec![2, 1, 1, 2],
        Permutation::from_image(&[1, 2, 0, 3]),
    )]
    #[case::restricted_equal_blocks(
        ClassKey::CisTrans,
        vec![2, 1, 1, 2],
        Permutation::from_image(&[1, 0, 2, 3]),
    )]
    fn test_coset_space_normalizer(
        #[case] key: ClassKey,
        #[case] frame: Vec<i32>,
        #[case] expected: Permutation,
    ) {
        let space = key.space();
        let normalizer = space.normalizer(&frame).unwrap();
        assert_eq!(normalizer, expected);
        assert!(space.allows(normalizer));
        assert_eq!(
            space.normalizer(&normalizer.act(&frame)),
            Some(Permutation::identity(space.degree()))
        );
    }

    #[rstest]
    #[case::short(ClassKey::Tetrahedral, vec![0, 1, 2])]
    #[case::long(ClassKey::CisTrans, vec![0, 1, 2, 3, 4])]
    fn test_coset_space_normalizer_degree(#[case] key: ClassKey, #[case] frame: Vec<i32>) {
        assert_eq!(key.space().normalizer(&frame), None);
    }

    #[rstest]
    #[case::tetrahedral(PermutationGroup::alternating(4), Decomposition::CanonicalRank, 2)]
    #[case::square_planar(PermutationGroup::dihedral(4), Decomposition::SquarePlanar, 3)]
    fn test_coset_space_count(
        #[case] group: PermutationGroup,
        #[case] decomposition: Decomposition,
        #[case] count: usize,
    ) {
        let degree = group.degree();
        let parent = PermutationGroup::symmetric(degree);
        assert_eq!(
            CosetSpace::new(
                parent,
                false,
                group,
                decomposition,
                Permutation::identity(degree),
            )
            .count(),
            count
        );
    }

    #[rstest]
    fn test_coset_space_coset_rep() {
        let space = CosetSpace::new(
            PermutationGroup::symmetric(4),
            false,
            PermutationGroup::alternating(4),
            Decomposition::CanonicalRank,
            Permutation::identity(4),
        );
        let even = Permutation::from_image(&[1, 2, 0, 3]);
        let odd = Permutation::from_image(&[1, 0, 2, 3]);
        assert_eq!(
            space.coset_rep(Permutation::identity(4)),
            space.coset_rep(even)
        );
        assert_ne!(
            space.coset_rep(Permutation::identity(4)),
            space.coset_rep(odd)
        );
    }

    #[rstest]
    #[case::canonical(PermutationGroup::alternating(4), Decomposition::CanonicalRank)]
    #[case::square_planar(PermutationGroup::dihedral(4), Decomposition::SquarePlanar)]
    fn test_coset_space_index_unindex(
        #[case] group: PermutationGroup,
        #[case] decomposition: Decomposition,
    ) {
        let degree = group.degree();
        let parent = PermutationGroup::symmetric(degree);
        let space = CosetSpace::new(
            parent,
            false,
            group,
            decomposition,
            Permutation::identity(degree),
        );
        for n in 0..space.count() as u32 {
            assert_eq!(space.index(space.unindex(n).unwrap()), Some(n));
        }
    }

    #[rstest]
    #[case::u_shape([0, 1, 2, 3], 0)]
    #[case::four_shape([1, 3, 0, 2], 1)]
    #[case::z_shape([1, 2, 0, 3], 2)]
    fn test_coset_space_index_square_planar(#[case] image: [usize; 4], #[case] expected: u32) {
        let space = CosetSpace::new(
            PermutationGroup::symmetric(4),
            false,
            PermutationGroup::dihedral(4),
            Decomposition::SquarePlanar,
            Permutation::identity(4),
        );
        assert_eq!(space.index(Permutation::from_image(&image)), Some(expected));
    }

    #[rstest]
    fn test_coset_space_cis_trans() {
        // parent D₄ = S₂ ≀ S₂ (within-side swaps + carbon swap); R = Klein four V; 8/4 = 2 cosets.
        let parent = PermutationGroup::generate(
            4,
            &[
                Permutation::from_image(&[1, 0, 2, 3]),
                Permutation::from_image(&[0, 1, 3, 2]),
                Permutation::from_image(&[2, 3, 0, 1]),
            ],
        );
        let group = PermutationGroup::generate(
            4,
            &[
                Permutation::from_image(&[1, 0, 3, 2]),
                Permutation::from_image(&[2, 3, 0, 1]),
            ],
        );
        let space = CosetSpace::new(
            parent,
            true,
            group,
            Decomposition::CanonicalRank,
            Permutation::identity(4),
        );
        assert_eq!(space.count(), 2);
        // The within-side swap (0 1) is the cis↔trans flip: a different coset from identity.
        let within_side_swap = Permutation::from_image(&[1, 0, 2, 3]);
        assert_ne!(
            space.coset_rep(Permutation::identity(4)),
            space.coset_rep(within_side_swap)
        );
        // The carbon swap (0 2)(1 3) ∈ R: same coset (writing the bond with either carbon first).
        let carbon_swap = Permutation::from_image(&[2, 3, 0, 1]);
        assert_eq!(
            space.coset_rep(Permutation::identity(4)),
            space.coset_rep(carbon_swap)
        );
    }

    #[rstest]
    #[case::coset_0(0, 1)]
    #[case::coset_1(1, 0)]
    fn test_coset_space_enantiomer(#[case] index: u32, #[case] expected: u32) {
        let space = CosetSpace::new(
            PermutationGroup::symmetric(4),
            false,
            PermutationGroup::alternating(4),
            Decomposition::CanonicalRank,
            Permutation::from_image(&[1, 0, 2, 3]),
        );
        assert_eq!(space.enantiomer(index), Some(expected));
    }

    #[rstest]
    #[case::improper_swap(Permutation::from_image(&[1, 0, 2, 3]), true)]
    #[case::identity(Permutation::identity(4), false)]
    fn test_coset_space_is_chiral(#[case] improper: Permutation, #[case] expected: bool) {
        let space = CosetSpace::new(
            PermutationGroup::symmetric(4),
            false,
            PermutationGroup::alternating(4),
            Decomposition::CanonicalRank,
            improper,
        );
        assert_eq!(space.is_chiral(), expected);
    }

    #[rstest]
    #[case::no_generators(ClassKey::Tetrahedral, vec![], Some(vec![0, 1]))]
    #[case::odd_swap(
        ClassKey::Tetrahedral,
        vec![Permutation::from_image(&[1, 0, 2, 3])],
        Some(vec![0, 0]),
    )]
    #[case::even_keeps(
        ClassKey::Tetrahedral,
        vec![Permutation::from_image(&[1, 2, 0, 3])],
        Some(vec![0, 1]),
    )]
    #[case::outside_parent(
        ClassKey::CisTrans,
        vec![Permutation::from_image(&[1, 2, 0, 3])],
        None,
    )]
    #[case::wrong_degree(
        ClassKey::Tetrahedral,
        vec![Permutation::identity(3)],
        None,
    )]
    fn test_coset_space_orbit_reps(
        #[case] key: ClassKey,
        #[case] generators: Vec<Permutation>,
        #[case] expected: Option<Vec<u32>>,
    ) {
        assert_eq!(key.space().orbit_reps(&generators), expected);
    }

    #[rstest]
    #[case::no_generators(ClassKey::Tetrahedral, 1, vec![], Some(1))]
    #[case::additional_generator(
        ClassKey::Tetrahedral,
        1,
        vec![Permutation::from_image(&[1, 0, 2, 3])],
        Some(0),
    )]
    #[case::index_out_of_range(ClassKey::Tetrahedral, 2, vec![], None)]
    #[case::outside_parent(
        ClassKey::CisTrans,
        0,
        vec![Permutation::from_image(&[1, 2, 0, 3])],
        None,
    )]
    #[case::wrong_degree(
        ClassKey::Tetrahedral,
        0,
        vec![Permutation::identity(3)],
        None,
    )]
    fn test_coset_space_observable_coset(
        #[case] key: ClassKey,
        #[case] index: u32,
        #[case] supergroup_generators: Vec<Permutation>,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(
            key.space().observable_coset(index, &supergroup_generators),
            expected
        );
    }
}
