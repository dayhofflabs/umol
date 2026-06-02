//! Stereo class keys and the interned coset-space registry.
//!
//! `space(key)` builds each `CosetSpace` once and leaks it for `'static`,
//! mirroring umol-msym's point-group registry. The geometry classes pin the
//! proper-rotation group as a permutation group acting on the ligand positions.

use std::collections::HashMap;
use std::str::FromStr;
use std::sync::{LazyLock, Mutex};
use std::{fmt, ptr};

use crate::coset::CosetSpace;
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
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl ClassKey {
    fn build(self) -> CosetSpace {
        let group = match self {
            ClassKey::Symmetric(n) => PermutationGroup::symmetric(n as usize),
            ClassKey::Alternating(n) => PermutationGroup::alternating(n as usize),
            ClassKey::Cyclic(n) => PermutationGroup::cyclic(n as usize),
            ClassKey::Dihedral(n) => PermutationGroup::dihedral(n as usize),
            ClassKey::Tetrahedral => PermutationGroup::alternating(4),
            ClassKey::CisTrans => PermutationGroup::alternating(2),
            ClassKey::SquarePlanar => PermutationGroup::dihedral(4),
            // 2 apical (0,1) + 3 equatorial (2,3,4); C₃ cycles equatorial,
            // C₂ swaps apical and two equatorial.
            ClassKey::TrigonalBipyramidal => PermutationGroup::generate(
                5,
                &[
                    Permutation::from_image(5, &[0, 1, 3, 4, 2]),
                    Permutation::from_image(5, &[1, 0, 2, 4, 3]),
                ],
            ),
            // 6 octahedron vertices (0,1=±x; 2,3=±y; 4,5=±z); C₄ about z,
            // C₃ about the (1,1,1) body diagonal.
            ClassKey::Octahedral => PermutationGroup::generate(
                6,
                &[
                    Permutation::from_image(6, &[2, 3, 1, 0, 4, 5]),
                    Permutation::from_image(6, &[2, 3, 4, 5, 0, 1]),
                ],
            ),
        };
        CosetSpace::new(group)
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
    #[case::square_planar(ClassKey::SquarePlanar, 3)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, 20)]
    #[case::octahedral(ClassKey::Octahedral, 30)]
    fn test_space_count(#[case] key: ClassKey, #[case] count: usize) {
        assert_eq!(space(key).count(), count);
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
