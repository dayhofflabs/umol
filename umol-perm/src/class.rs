//! Stereo class keys and the interned coset-space registry.
//!
//! `ClassKey::space` builds each `CosetSpace` once and leaks it for `'static`,
//! mirroring umol-msym's point-group registry. The geometry classes pin the
//! proper-rotation group as a permutation group acting on the ligand positions.

use std::collections::HashMap;
use std::fmt;
use std::ops::RangeInclusive;
use std::str::FromStr;
use std::sync::{LazyLock, Mutex};

use crate::coset::{CosetSpace, Decomposition};
use crate::error::ParseClassKeyError;
use crate::group::PermutationGroup;
use crate::permutation::{Permutation, MAX_DEGREE};

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
        let is_partitioned = matches!(self, ClassKey::CisTrans | ClassKey::Axial);
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
                        Permutation::from_image(&[1, 0, 2, 3]),
                        Permutation::from_image(&[0, 1, 3, 2]),
                        Permutation::from_image(&[2, 3, 0, 1]),
                    ],
                ),
                PermutationGroup::generate(
                    4,
                    &[
                        Permutation::from_image(&[1, 0, 3, 2]),
                        Permutation::from_image(&[2, 3, 0, 1]),
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
                        Permutation::from_image(&[0, 2, 3, 1, 4]),
                        Permutation::from_image(&[4, 3, 2, 1, 0]),
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
                        Permutation::from_image(&[0, 2, 3, 4, 1, 5]),
                        Permutation::from_image(&[1, 2, 0, 4, 5, 3]),
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
            ClassKey::Tetrahedral | ClassKey::Axial => Permutation::from_image(&[1, 0, 2, 3]),
            ClassKey::TrigonalBipyramidal => Permutation::from_image(&[4, 1, 2, 3, 0]),
            ClassKey::Octahedral => Permutation::from_image(&[5, 1, 2, 3, 4, 0]),
            ClassKey::CisTrans | ClassKey::SquarePlanar => Permutation::identity(4),
            ClassKey::Symmetric(n)
            | ClassKey::Alternating(n)
            | ClassKey::Cyclic(n)
            | ClassKey::Dihedral(n) => Permutation::identity(n as usize),
        };
        CosetSpace::new(parent, is_partitioned, group, decomposition, improper)
    }

    /// The interned coset space for this class, built once and leaked for
    /// `'static`.
    ///
    /// # Establishes
    ///
    /// The returned reference is canonical for its key: every call with an equal key returns the
    /// identical pointer, so pointer identity coincides with key identity. `Coset` equality
    /// relies on this.
    ///
    /// # Panics
    ///
    /// Panics when the key's degree lies outside [`Self::degree_domain`]. Parsing rejects such
    /// keys, so the panic is reachable only from directly constructed ones. A build panic poisons
    /// the shared registry: every later `space` call on any key then panics.
    pub fn space(self) -> &'static CosetSpace {
        let mut registry = REGISTRY.lock().expect("coset-space registry poisoned");
        if let Some(&interned) = registry.get(&self) {
            return interned;
        }
        let interned: &'static CosetSpace = Box::leak(Box::new(self.build()));
        registry.insert(self, interned);
        interned
    }

    /// The supported degree domain; a key whose degree lies outside it has no
    /// coset space.
    pub fn degree_domain(self) -> RangeInclusive<u8> {
        match self {
            ClassKey::Symmetric(_)
            | ClassKey::Alternating(_)
            | ClassKey::Cyclic(_)
            | ClassKey::Dihedral(_) => 1..=(MAX_DEGREE as u8),
            ClassKey::Tetrahedral
            | ClassKey::CisTrans
            | ClassKey::Axial
            | ClassKey::SquarePlanar => 4..=4,
            ClassKey::TrigonalBipyramidal => 5..=5,
            ClassKey::Octahedral => 6..=6,
        }
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

/// Parses the rendered key text, e.g. `"TH"` or `"Sym4"`.
///
/// # Errors
///
/// Returns `UnknownClassKey` for an unrecognized prefix, `InvalidDegree` for a
/// missing, malformed, or out-of-domain degree, and `DegreeTooLarge` for a
/// degree above `MAX_DEGREE`.
///
/// # Semantic properties
///
/// Parsing inverts rendering: `key.to_string().parse() == Ok(key)` for every
/// key inside its degree domain.
impl FromStr for ClassKey {
    type Err = ParseClassKeyError;

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
        enum Family {
            Symmetric,
            Alternating,
            Cyclic,
            Dihedral,
        }
        let (degree, family) = if let Some(degree) = s.strip_prefix("Sym") {
            (degree, Family::Symmetric)
        } else if let Some(degree) = s.strip_prefix("Alt") {
            (degree, Family::Alternating)
        } else if let Some(degree) = s.strip_prefix("Cyc") {
            (degree, Family::Cyclic)
        } else if let Some(degree) = s.strip_prefix("Dih") {
            (degree, Family::Dihedral)
        } else {
            return Err(ParseClassKeyError::UnknownClassKey {
                input: s.to_string(),
            });
        };
        let degree = degree
            .parse::<usize>()
            .map_err(|_| ParseClassKeyError::InvalidDegree {
                input: s.to_string(),
            })?;
        if degree > MAX_DEGREE {
            return Err(ParseClassKeyError::DegreeTooLarge {
                degree,
                maximum: MAX_DEGREE,
            });
        }
        let key = match family {
            Family::Symmetric => ClassKey::Symmetric(degree as u8),
            Family::Alternating => ClassKey::Alternating(degree as u8),
            Family::Cyclic => ClassKey::Cyclic(degree as u8),
            Family::Dihedral => ClassKey::Dihedral(degree as u8),
        };
        if !key.degree_domain().contains(&(degree as u8)) {
            return Err(ParseClassKeyError::InvalidDegree {
                input: s.to_string(),
            });
        }
        Ok(key)
    }
}

static REGISTRY: LazyLock<Mutex<HashMap<ClassKey, &'static CosetSpace>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

#[cfg(test)]
mod tests {
    use std::ptr;

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
    fn test_class_key_space(#[case] key: ClassKey, #[case] count: usize) {
        assert_eq!(key.space().count(), count);
    }

    #[rstest]
    #[case::tetrahedral(ClassKey::Tetrahedral, true)]
    #[case::cis_trans(ClassKey::CisTrans, false)]
    #[case::axial(ClassKey::Axial, true)]
    #[case::square_planar(ClassKey::SquarePlanar, false)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, true)]
    #[case::octahedral(ClassKey::Octahedral, true)]
    fn test_class_key_space_chirality(#[case] key: ClassKey, #[case] expected: bool) {
        assert_eq!(key.space().is_chiral(), expected);
    }

    #[rstest]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal, 5, 6)]
    #[case::octahedral(ClassKey::Octahedral, 6, 24)]
    fn test_class_key_space_group(
        #[case] key: ClassKey,
        #[case] degree: usize,
        #[case] order: usize,
    ) {
        let group = key.space().group();
        assert_eq!(group.degree(), degree);
        assert_eq!(group.order(), order);
    }

    #[rstest]
    #[case::tetrahedral(ClassKey::Tetrahedral)]
    #[case::cis_trans(ClassKey::CisTrans)]
    #[case::square_planar(ClassKey::SquarePlanar)]
    #[case::trigonal_bipyramidal(ClassKey::TrigonalBipyramidal)]
    #[case::octahedral(ClassKey::Octahedral)]
    fn test_class_key_space_index_roundtrip(#[case] key: ClassKey) {
        let space = key.space();
        for n in 0..space.count() as u32 {
            assert_eq!(space.index(space.unindex(n).unwrap()), Some(n));
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
    fn test_class_key_space_reindex(
        #[case] key: ClassKey,
        #[case] from_index: u32,
        #[case] order: Vec<&str>,
        #[case] relabeled: Vec<&str>,
        #[case] to_index: u32,
    ) {
        let relabeling = Permutation::between(&order, &relabeled)
            .expect("fixed frames are orderings of the same ligands");
        assert_eq!(key.space().reindex(from_index, relabeling), Some(to_index));
    }

    #[rstest]
    fn test_class_key_space_interning() {
        assert!(ptr::eq(
            ClassKey::Octahedral.space(),
            ClassKey::Octahedral.space()
        ));
    }

    #[rstest]
    fn test_class_key_degree_domain() {
        let family_keys = (0..=u8::MAX).flat_map(|degree| {
            [
                (ClassKey::Symmetric(degree), degree),
                (ClassKey::Alternating(degree), degree),
                (ClassKey::Cyclic(degree), degree),
                (ClassKey::Dihedral(degree), degree),
            ]
        });
        let geometry_keys = [
            (ClassKey::Tetrahedral, 4),
            (ClassKey::CisTrans, 4),
            (ClassKey::Axial, 4),
            (ClassKey::SquarePlanar, 4),
            (ClassKey::TrigonalBipyramidal, 5),
            (ClassKey::Octahedral, 6),
        ];
        for (key, degree) in family_keys.chain(geometry_keys) {
            let in_domain = key.degree_domain().contains(&degree);
            let parsed = key.to_string().parse::<ClassKey>();
            assert_eq!(parsed.is_ok(), in_domain, "{key}");
            if in_domain {
                assert_eq!(parsed, Ok(key), "{key}");
                assert_eq!(key.space().degree(), degree as usize, "{key}");
            }
        }
    }

    #[rstest]
    #[case::symmetric(ClassKey::Symmetric(4), "Sym4")]
    #[case::alternating(ClassKey::Alternating(4), "Alt4")]
    #[case::cyclic(ClassKey::Cyclic(5), "Cyc5")]
    #[case::dihedral(ClassKey::Dihedral(5), "Dih5")]
    #[case::tetrahedral(ClassKey::Tetrahedral, "TH")]
    #[case::axial(ClassKey::Axial, "AX")]
    #[case::octahedral(ClassKey::Octahedral, "OH")]
    fn test_class_key_display_roundtrip(#[case] key: ClassKey, #[case] text: &str) {
        assert_eq!(key.to_string(), text);
        assert_eq!(ClassKey::from_str(text), Ok(key));
    }

    #[rstest]
    #[case::unknown(
        "Xyz3",
        ParseClassKeyError::UnknownClassKey { input: "Xyz3".to_string() },
    )]
    #[case::no_degree(
        "Sym",
        ParseClassKeyError::InvalidDegree { input: "Sym".to_string() },
    )]
    #[case::malformed_degree(
        "Altfour",
        ParseClassKeyError::InvalidDegree { input: "Altfour".to_string() },
    )]
    #[case::degree_too_large(
        "Cyc7",
        ParseClassKeyError::DegreeTooLarge { degree: 7, maximum: MAX_DEGREE },
    )]
    #[case::zero_degree_cyclic(
        "Cyc0",
        ParseClassKeyError::InvalidDegree { input: "Cyc0".to_string() },
    )]
    #[case::zero_degree_symmetric(
        "Sym0",
        ParseClassKeyError::InvalidDegree { input: "Sym0".to_string() },
    )]
    #[case::degree_exceeds_u8(
        "Sym256",
        ParseClassKeyError::DegreeTooLarge { degree: 256, maximum: MAX_DEGREE },
    )]
    fn test_class_key_from_str_error(#[case] text: &str, #[case] expected: ParseClassKeyError) {
        assert_eq!(ClassKey::from_str(text), Err(expected));
    }
}
