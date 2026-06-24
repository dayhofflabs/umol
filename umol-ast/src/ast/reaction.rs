//! Reaction rule AST: L ← K → R as two molecule ASTs plus atom and stereo maps.
//!
//! Interim, SMIRKS semantics (partial atom map, injective in both directions); the
//! full reaction redesign (concrete reactions, non-injective templates, overlay
//! rewrite) is separate. See discussion doc 127.

use umol_perm::Permutation;

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::molecule::MoleculeAst;

/// A reaction as a double-pushout rewrite — homoiconic, holding a concrete
/// reaction or a rule.
///
/// `lhs` and `rhs` are the left- and right-hand-side molecule patterns. `atom_map`
/// pairs the interface K: `(Some, Some)` is preserved, `(Some, None)` is deleted,
/// `(None, Some)` is created. `stereo_atoms` / `stereo_bonds` carry how each stereo
/// element transforms across the rule (see [`StereoAtomCorrespondence`]). Bonds and
/// overlay relations in K are inferred from the topology of L and R over the mapped
/// atoms.
#[derive(Clone, Debug)]
pub struct ReactionAst {
    pub lhs: MoleculeAst,
    pub rhs: MoleculeAst,
    pub atom_map: Vec<(Option<AtomId>, Option<AtomId>)>,
    pub stereo_atoms: Vec<StereoAtomCorrespondence>,
    pub stereo_bonds: Vec<StereoBondCorrespondence>,
}

/// Generates an L→R correspondence type for one stereo-element kind. `lhs`/`rhs`
/// give presence (a missing side is create / destroy); when both are present,
/// `frame` is the permutation of the ordered ligand frame (including virtual
/// ligands) that carries the L coset to the R coset through the geometry's coset
/// algebra. A `None` frame means the frame is preserved (identity).
macro_rules! stereo_correspondence {
    ($(#[$meta:meta])* $name:ident, $id:ty) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name {
            pub lhs: Option<$id>,
            pub rhs: Option<$id>,
            pub frame: Option<Permutation>,
        }

        impl $name {
            /// Stereo element preserved unchanged across the rule.
            pub fn preserve(lhs: $id, rhs: $id) -> Self {
                Self { lhs: Some(lhs), rhs: Some(rhs), frame: None }
            }

            /// Stereo element created on the right-hand side (no L preimage).
            pub fn create(rhs: $id) -> Self {
                Self { lhs: None, rhs: Some(rhs), frame: None }
            }

            /// Stereo element destroyed from the left-hand side (no R image).
            pub fn destroy(lhs: $id) -> Self {
                Self { lhs: Some(lhs), rhs: None, frame: None }
            }

            /// Stereo element whose ligand frame is transposed across the rule;
            /// `frame` is the transposition taking the L frame to the R frame.
            pub fn swap(lhs: $id, rhs: $id, frame: Permutation) -> Self {
                Self { lhs: Some(lhs), rhs: Some(rhs), frame: Some(frame) }
            }
        }
    };
}

stereo_correspondence!(
    /// L→R correspondence of a tetrahedral (atom) stereo element.
    StereoAtomCorrespondence,
    StereoAtomId
);
stereo_correspondence!(
    /// L→R correspondence of a cis/trans (bond) stereo element.
    StereoBondCorrespondence,
    StereoBondId
);

/// Assignment mapping L-entity ids to target (G) entity ids. Produced by
/// substructure matching; consumed by `apply_rule`.
#[derive(Clone, Debug, Default)]
pub struct Assignment {
    pub atoms: Vec<(AtomId, AtomId)>,
    pub bonds: Vec<(BondId, BondId)>,
    pub dative_bonds: Vec<(DativeBondId, DativeBondId)>,
    pub aromatic_systems: Vec<(AromaticSystemId, AromaticSystemId)>,
    pub multicenter_bonds: Vec<(MulticenterBondId, MulticenterBondId)>,
    pub noncovalent_bonds: Vec<(NoncovalentBondId, NoncovalentBondId)>,
    pub stereo_atoms: Vec<(StereoAtomId, StereoAtomId)>,
    pub stereo_bonds: Vec<(StereoBondId, StereoBondId)>,
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::preserve(
        StereoAtomCorrespondence::preserve(StereoAtomId(0), StereoAtomId(1)),
        StereoAtomCorrespondence { lhs: Some(StereoAtomId(0)), rhs: Some(StereoAtomId(1)), frame: None }
    )]
    #[case::create(
        StereoAtomCorrespondence::create(StereoAtomId(2)),
        StereoAtomCorrespondence { lhs: None, rhs: Some(StereoAtomId(2)), frame: None }
    )]
    #[case::destroy(
        StereoAtomCorrespondence::destroy(StereoAtomId(3)),
        StereoAtomCorrespondence { lhs: Some(StereoAtomId(3)), rhs: None, frame: None }
    )]
    #[case::swap(
        StereoAtomCorrespondence::swap(
            StereoAtomId(0),
            StereoAtomId(1),
            Permutation::from_image(4, &[1, 0, 2, 3]),
        ),
        StereoAtomCorrespondence {
            lhs: Some(StereoAtomId(0)),
            rhs: Some(StereoAtomId(1)),
            frame: Some(Permutation::from_image(4, &[1, 0, 2, 3])),
        }
    )]
    fn test_stereo_atom_correspondence(
        #[case] actual: StereoAtomCorrespondence,
        #[case] expected: StereoAtomCorrespondence,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::preserve(
        StereoBondCorrespondence::preserve(StereoBondId(0), StereoBondId(1)),
        StereoBondCorrespondence { lhs: Some(StereoBondId(0)), rhs: Some(StereoBondId(1)), frame: None }
    )]
    #[case::create(
        StereoBondCorrespondence::create(StereoBondId(2)),
        StereoBondCorrespondence { lhs: None, rhs: Some(StereoBondId(2)), frame: None }
    )]
    #[case::destroy(
        StereoBondCorrespondence::destroy(StereoBondId(3)),
        StereoBondCorrespondence { lhs: Some(StereoBondId(3)), rhs: None, frame: None }
    )]
    #[case::swap(
        StereoBondCorrespondence::swap(
            StereoBondId(0),
            StereoBondId(1),
            Permutation::from_image(4, &[1, 0, 2, 3]),
        ),
        StereoBondCorrespondence {
            lhs: Some(StereoBondId(0)),
            rhs: Some(StereoBondId(1)),
            frame: Some(Permutation::from_image(4, &[1, 0, 2, 3])),
        }
    )]
    fn test_stereo_bond_correspondence(
        #[case] actual: StereoBondCorrespondence,
        #[case] expected: StereoBondCorrespondence,
    ) {
        assert_eq!(actual, expected);
    }
}
