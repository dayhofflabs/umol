//! Stereo class keys and the interned coset-space registry.
//!
//! `space(key)` builds each `CosetSpace` once and leaks it for `'static`,
//! mirroring umol-msym's point-group registry. The geometry classes pin the
//! proper-rotation group as a permutation group acting on the ligand positions.

use std::collections::HashMap;
use std::str::FromStr;
use std::sync::{LazyLock, Mutex};
use std::{fmt, ptr};

use crate::coset::{CosetSpace, Decomposition};
use crate::group::PermutationGroup;
use crate::permutation::Permutation;

/// A stereo class — either a named permutation-group family (with its degree)
/// or a coordination geometry whose proper-rotation group fixes the cosets.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ClassKey {
    Symmetric(u8),
    Alternating(u8),
    Cyclic(u8),
    Dihedral(u8),
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl ClassKey {
    fn build(self) -> CosetSpace {
        // Each arm gives (parent P, symmetry R, decomposition). P is the group of realizable
        // arrangements: Sₙ for the geometry classes, the partition group D₄ for cis/trans.
        let (parent, group, decomposition) = match self {
            ClassKey::Symmetric(n) => (
                PermutationGroup::symmetric(n as usize),
                PermutationGroup::symmetric(n as usize),
                Decomposition::CanonicalRank,
            ),
            ClassKey::Alternating(n) => (
                PermutationGroup::symmetric(n as usize),
                PermutationGroup::alternating(n as usize),
                Decomposition::CanonicalRank,
            ),
            ClassKey::Cyclic(n) => (
                PermutationGroup::symmetric(n as usize),
                PermutationGroup::cyclic(n as usize),
                Decomposition::CanonicalRank,
            ),
            ClassKey::Dihedral(n) => (
                PermutationGroup::symmetric(n as usize),
                PermutationGroup::dihedral(n as usize),
                Decomposition::CanonicalRank,
            ),
            ClassKey::Tetrahedral => (
                PermutationGroup::symmetric(4),
                PermutationGroup::alternating(4),
                Decomposition::CanonicalRank,
            ),
            // Substituents are bonded to fixed sp² carbons and the bond may be written with
            // either carbon first, so the parent is S₂ ≀ S₂ = D₄ — within-side swaps (0 1), (2 3)
            // and the carbon swap (0 2)(1 3). R is the Klein four V (face flip (0 1)(2 3) and
            // carbon swap), giving 8/4 = 2 cosets (cis, trans).
            // Axial (allene / biaryl …) shares this coset space but is chiral: its
            // improper generator swaps the two cosets (see the `improper` match below).
            ClassKey::CisTrans | ClassKey::Axial => (
                PermutationGroup::generate(
                    4,
                    &[
                        Permutation::from_image(4, &[1, 0, 2, 3]),
                        Permutation::from_image(4, &[0, 1, 3, 2]),
                        Permutation::from_image(4, &[2, 3, 0, 1]),
                    ],
                ),
                PermutationGroup::generate(
                    4,
                    &[
                        Permutation::from_image(4, &[1, 0, 3, 2]),
                        Permutation::from_image(4, &[2, 3, 0, 1]),
                    ],
                ),
                Decomposition::CanonicalRank,
            ),
            ClassKey::SquarePlanar => (
                PermutationGroup::symmetric(4),
                PermutationGroup::dihedral(4),
                Decomposition::SquarePlanar,
            ),
            // 2 axial (0,4) + 3 equatorial (1,2,3); C₃ cycles equatorial,
            // C₂ swaps the axial pair and two equatorial.
            ClassKey::TrigonalBipyramidal => (
                PermutationGroup::symmetric(5),
                PermutationGroup::generate(
                    5,
                    &[
                        Permutation::from_image(5, &[0, 2, 3, 1, 4]),
                        Permutation::from_image(5, &[4, 3, 2, 1, 0]),
                    ],
                ),
                Decomposition::TrigonalBipyramidal,
            ),
            // 6 octahedron vertices (0,5 = axial; 1,2,3,4 = equatorial square);
            // C₄ about the 0–5 axis, C₃ about a body diagonal.
            ClassKey::Octahedral => (
                PermutationGroup::symmetric(6),
                PermutationGroup::generate(
                    6,
                    &[
                        Permutation::from_image(6, &[0, 2, 3, 4, 1, 5]),
                        Permutation::from_image(6, &[1, 2, 0, 4, 5, 3]),
                    ],
                ),
                Decomposition::Octahedral,
            ),
        };
        // The orientation-reversing generator.
        // TH σ_d (0 1); Axial the cis↔trans flip (0 1); TB/OH σ_h swapping the axial pair.
        // Achiral classes (CisTrans, SP, and the abstract families) take the identity,
        // so `is_chiral` is false for them.
        // NOTE: σ_h gives the correct `is_chiral` for TB/OH; the exact enantiomer *pairing* is
        // class-geometry data still to be verified against the OpenSMILES @↔@@ numbering.
        let improper = match self {
            ClassKey::Tetrahedral | ClassKey::Axial => Permutation::from_image(4, &[1, 0, 2, 3]),
            ClassKey::TrigonalBipyramidal => Permutation::from_image(5, &[4, 1, 2, 3, 0]),
            ClassKey::Octahedral => Permutation::from_image(6, &[5, 1, 2, 3, 4, 0]),
            ClassKey::CisTrans | ClassKey::SquarePlanar => Permutation::identity(4),
            ClassKey::Symmetric(n)
            | ClassKey::Alternating(n)
            | ClassKey::Cyclic(n)
            | ClassKey::Dihedral(n) => Permutation::identity(n as usize),
        };
        CosetSpace::new(parent, group, decomposition, improper)
    }
}

impl fmt::Display for ClassKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClassKey::Symmetric(n) => write!(f, "Sym{n}"),
            ClassKey::Alternating(n) => write!(f, "Alt{n}"),
            ClassKey::Cyclic(n) => write!(f, "Cyc{n}"),
            ClassKey::Dihedral(n) => write!(f, "Dih{n}"),
            ClassKey::Tetrahedral => write!(f, "TH"),
            ClassKey::CisTrans => write!(f, "CT"),
            ClassKey::Axial => write!(f, "AX"),
            ClassKey::SquarePlanar => write!(f, "SP"),
            ClassKey::TrigonalBipyramidal => write!(f, "TB"),
            ClassKey::Octahedral => write!(f, "OH"),
        }
    }
}

impl FromStr for ClassKey {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "TH" => return Ok(ClassKey::Tetrahedral),
            "CT" => return Ok(ClassKey::CisTrans),
            "AX" => return Ok(ClassKey::Axial),
            "SP" => return Ok(ClassKey::SquarePlanar),
            "TB" => return Ok(ClassKey::TrigonalBipyramidal),
            "OH" => return Ok(ClassKey::Octahedral),
            _ => {}
        }
        let split = s
            .find(|c: char| c.is_ascii_digit())
            .ok_or_else(|| format!("unknown class key: {s}"))?;
        let (prefix, degree) = s.split_at(split);
        let n: u8 = degree
            .parse()
            .map_err(|_| format!("bad degree in class key: {s}"))?;
        match prefix {
            "Sym" => Ok(ClassKey::Symmetric(n)),
            "Alt" => Ok(ClassKey::Alternating(n)),
            "Cyc" => Ok(ClassKey::Cyclic(n)),
            "Dih" => Ok(ClassKey::Dihedral(n)),
            _ => Err(format!("unknown class family: {prefix}")),
        }
    }
}

static REGISTRY: LazyLock<Mutex<HashMap<ClassKey, &'static CosetSpace>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

/// The interned coset space for `key`, built once and leaked for `'static`.
pub fn space(key: ClassKey) -> &'static CosetSpace {
    let mut registry = REGISTRY.lock().expect("coset-space registry poisoned");
    if let Some(&interned) = registry.get(&key) {
        return interned;
    }
    let leaked: &'static CosetSpace = Box::leak(Box::new(key.build()));
    registry.insert(key, leaked);
    leaked
}

/// A configuration: an index into the cosets of an interned space. Identity is
/// the interned space pointer plus the index.
#[derive(Clone, Copy, Debug)]
pub struct Coset {
    space: &'static CosetSpace,
    index: u32,
}

impl Coset {
    pub fn new(key: ClassKey, index: u32) -> Self {
        Self {
            space: space(key),
            index,
        }
    }

    pub fn index(self) -> u32 {
        self.index
    }

    pub fn space(self) -> &'static CosetSpace {
        self.space
    }
}

impl PartialEq for Coset {
    fn eq(&self, other: &Self) -> bool {
        ptr::eq(self.space, other.space) && self.index == other.index
    }
}

impl Eq for Coset {}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::tetrahedral(ClassKey::Tetrahedral, 2)]
    #[case::cis_trans(ClassKey::CisTrans, 2)]
    #[case::axial(ClassKey::Axial, 2)]
    #[case::square_planar(ClassKey::SquarePlanar, 3)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, 20)]
    #[case::octahedral(ClassKey::Octahedral, 30)]
    fn test_space_count(#[case] key: ClassKey, #[case] count: usize) {
        assert_eq!(space(key).count(), count);
    }

    #[rstest]
    #[case::tetrahedral(ClassKey::Tetrahedral, true)]
    #[case::cis_trans(ClassKey::CisTrans, false)]
    #[case::axial(ClassKey::Axial, true)]
    #[case::square_planar(ClassKey::SquarePlanar, false)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, true)]
    #[case::octahedral(ClassKey::Octahedral, true)]
    fn test_space_is_chiral(#[case] key: ClassKey, #[case] expected: bool) {
        assert_eq!(space(key).is_chiral(), expected);
    }

    #[rstest]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, 5, 6)]
    #[case::octahedral(ClassKey::Octahedral, 6, 24)]
    fn test_geometry_group_order(
        #[case] key: ClassKey,
        #[case] degree: usize,
        #[case] order: usize,
    ) {
        let group = space(key).group();
        assert_eq!(group.degree(), degree);
        assert_eq!(group.order(), order);
    }

    #[rstest]
    #[case::tetrahedral(ClassKey::Tetrahedral)]
    #[case::cis_trans(ClassKey::CisTrans)]
    #[case::square_planar(ClassKey::SquarePlanar)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal)]
    #[case::octahedral(ClassKey::Octahedral)]
    fn test_space_index_unindex(#[case] key: ClassKey) {
        let space = space(key);
        for n in 0..space.count() as u32 {
            assert_eq!(space.index(space.unindex(n)), n);
        }
    }

    // Each case is an "all equivalent" SMILES pair from the OpenSMILES spec: the
    // same structure written with two neighbor orders and the two `@`-numbers
    // (0-based) they require.
    #[rstest]
    #[case::tb1_tb2(ClassKey::TrigonalBipyramidal, 0, vec!["S", "F", "Cl", "Br", "N"], vec!["S", "Br", "Cl", "F", "N"], 1)]
    #[case::tb5_tb10(ClassKey::TrigonalBipyramidal, 4, vec!["S", "F", "N", "Cl", "Br"], vec!["F", "S", "Cl", "N", "Br"], 9)]
    #[case::tb15_tb20(ClassKey::TrigonalBipyramidal, 14, vec!["F", "Cl", "S", "Br", "N"], vec!["Br", "Cl", "S", "F", "N"], 19)]
    #[case::oh1_oh2(ClassKey::Octahedral, 0, vec!["C", "F", "Cl", "Br", "I", "S"], vec!["F", "S", "I", "C", "Cl", "Br"], 1)]
    #[case::oh5_oh9(ClassKey::Octahedral, 4, vec!["S", "F", "I", "Cl", "C", "Br"], vec!["Br", "C", "S", "Cl", "F", "I"], 8)]
    #[case::oh12_oh15(ClassKey::Octahedral, 11, vec!["Br", "Cl", "I", "F", "S", "C"], vec!["Cl", "C", "Br", "F", "I", "S"], 14)]
    #[case::oh19_oh27(ClassKey::Octahedral, 18, vec!["Cl", "C", "I", "F", "S", "Br"], vec!["I", "Cl", "Br", "F", "S", "C"], 26)]
    fn test_space_reindex(
        #[case] key: ClassKey,
        #[case] from_index: u32,
        #[case] order: Vec<&str>,
        #[case] relabeled: Vec<&str>,
        #[case] to_index: u32,
    ) {
        let relabeling = Permutation::between(&order, &relabeled);
        assert_eq!(space(key).reindex(from_index, relabeling), to_index);
    }

    #[rstest]
    fn test_space_interned() {
        assert!(ptr::eq(
            space(ClassKey::Octahedral),
            space(ClassKey::Octahedral)
        ));
    }

    #[rstest]
    #[case::symmetric(ClassKey::Symmetric(4), "Sym4")]
    #[case::dihedral(ClassKey::Dihedral(5), "Dih5")]
    #[case::tetrahedral(ClassKey::Tetrahedral, "TH")]
    #[case::axial(ClassKey::Axial, "AX")]
    #[case::octahedral(ClassKey::Octahedral, "OH")]
    fn test_class_key_display(#[case] key: ClassKey, #[case] text: &str) {
        assert_eq!(key.to_string(), text);
        assert_eq!(ClassKey::from_str(text), Ok(key));
    }

    #[rstest]
    #[case::unknown("Xyz3")]
    #[case::no_degree("Sym")]
    fn test_class_key_from_str_error(#[case] text: &str) {
        assert!(ClassKey::from_str(text).is_err());
    }
}
