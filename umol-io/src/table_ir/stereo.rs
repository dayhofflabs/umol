//! Stereochemistry records and metadata for TableIR.

#[cfg_attr(not(test), expect(dead_code))]
pub(crate) mod derive;

/// Whether the molecule's stereo descriptors fix the absolute configuration or
/// only the relative one. Populated from format-specific flags:
/// - CTFile counts chiral flag (`ccc`)
/// - CXSMILES enhanced stereo markers (`a:` and `r`)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConfigurationScope {
    Absolute,
    Relative,
}

/// An explicit tetrahedral configuration in an ordered ligand frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoAtom {
    pub atom: u32,
    pub ligands: Vec<StereoLigand>,
    pub winding: Winding,
}

/// An actual neighbor or a virtual ligand borne by a stereo atom's site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoLigand {
    Atom(u32),
    ImplicitHydrogen,
    LonePair,
}

/// Winding of the last three tetrahedral ligands with the first toward the viewer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Winding {
    Clockwise,
    CounterClockwise,
}

/// A stereo assertion at a bond in the owning TableIR molecule.
///
/// This is an open record. The first consumer requiring a meaningful assertion checks the
/// site index, reference incidence, and consistency with other assertions in that molecule.
/// Construction does not check molecular context or chemical validity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoBond {
    /// Index in the owning molecule's bond table; transported when that table is reordered.
    pub bond: u32,
    pub configuration: BondConfiguration,
}

/// A site-only indefinite assertion or a definite configuration in an actual-atom frame.
///
/// Absence of a stereo-bond record means no assertion; it is distinct from [`Self::Either`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondConfiguration {
    /// Explicitly indefinite configuration, corresponding to CTfile Either, without references.
    Either,
    /// A definite relation between one actual substituent at each ordered site endpoint.
    Framed {
        /// Atom indices in the owning molecule. Position 0 is a substituent of the site's
        /// [`AtomPair::first`](super::AtomPair::first) endpoint; position 1 belongs to
        /// [`AtomPair::second`](super::AtomPair::second). Neither reference is the other endpoint.
        ///
        /// Renumbering atoms preserves reference identity and exchanges these slots if the
        /// ordered endpoints reverse. Choosing minimum-index references is a producer policy,
        /// not a requirement on independently supplied frames.
        references: [u32; 2],
        relation: BondRelation,
    },
}

/// The relative sides of two reference substituents at a stereo bond.
///
/// This is frame-relative, independent of CIP ranking, drawing orientation, and SMILES traversal.
/// In complete endpoint ligand blocks `[a, b]; [c, d]` with references `[a, c]`, SameSide is
/// [`ClassKey::CisTrans`](umol_perm::ClassKey::CisTrans) coset 0 and OppositeSide is coset 1.
/// The correspondence does not depend on enum discriminants.
///
/// # Semantic properties
///
/// Selecting the complementary substituent at one endpoint flips the relation. Selecting both
/// complements preserves it. Exchanging both endpoints and their reference slots preserves it.
/// Moving a ligand between endpoint blocks is not a frame change of the same site incidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondRelation {
    SameSide,
    OppositeSide,
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{FrameTransport, StereoBondForm, StereoKind};
    use umol_perm::{ClassKey, Permutation};

    use super::{BondConfiguration, BondRelation, StereoBond};
    use crate::table_ir::AtomPair;

    #[rstest]
    #[case::identity(
        [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4],
        AtomPair::new(1, 4), 2, [0, 3], [0, 2, 3, 5],
    )]
    #[case::bond_order(
        [0, 1, 2, 3, 4, 5], [2, 4, 1, 0, 3],
        AtomPair::new(1, 4), 1, [0, 3], [0, 2, 3, 5],
    )]
    #[case::reference_order(
        [2, 1, 0, 5, 4, 3], [0, 1, 2, 3, 4],
        AtomPair::new(1, 4), 2, [2, 5], [2, 0, 5, 3],
    )]
    #[case::endpoint_order(
        [0, 4, 2, 3, 1, 5], [0, 1, 2, 3, 4],
        AtomPair::new(1, 4), 2, [3, 0], [3, 5, 0, 2],
    )]
    #[case::atoms_and_bonds(
        [5, 4, 3, 2, 1, 0], [2, 4, 1, 0, 3],
        AtomPair::new(1, 4), 1, [2, 5], [2, 0, 5, 3],
    )]
    #[case::new_endpoints(
        [1, 5, 4, 0, 2, 3], [4, 3, 0, 2, 1],
        AtomPair::new(2, 5), 0, [0, 1], [0, 3, 1, 4],
    )]
    fn test_stereo_bond_frame_remapping(
        #[case] atom_image: [usize; 6],
        #[case] bond_image: [usize; 5],
        #[case] expected_endpoints: AtomPair,
        #[case] expected_bond: u32,
        #[case] expected_references: [u32; 2],
        #[case] expected_frame: [u32; 4],
        #[values(BondConfiguration::Either,
            BondConfiguration::Framed { references: [0, 3], relation: BondRelation::SameSide },
            BondConfiguration::Framed { references: [0, 3], relation: BondRelation::OppositeSide })]
        configuration: BondConfiguration,
    ) {
        let atoms = Permutation::from_image(&atom_image);
        let bonds = Permutation::from_image(&bond_image);
        let source = StereoBond {
            bond: 2,
            configuration,
        };
        let first = atoms.apply(1) as u32;
        let second = atoms.apply(4) as u32;
        assert_eq!(AtomPair::new(first, second), expected_endpoints);
        let configuration = match source.configuration {
            BondConfiguration::Either => BondConfiguration::Either,
            BondConfiguration::Framed {
                references,
                relation,
            } => {
                let mut references = references.map(|atom| atoms.apply(atom as usize) as u32);
                if first > second {
                    references.swap(0, 1);
                }
                BondConfiguration::Framed {
                    references,
                    relation,
                }
            }
        };
        let expected_configuration = match source.configuration {
            BondConfiguration::Either => BondConfiguration::Either,
            BondConfiguration::Framed { relation, .. } => BondConfiguration::Framed {
                references: expected_references,
                relation,
            },
        };
        assert_eq!(
            StereoBond {
                bond: bonds.apply(source.bond as usize) as u32,
                configuration
            },
            StereoBond {
                bond: expected_bond,
                configuration: expected_configuration
            },
        );

        let renamed = [0, 2, 3, 5].map(|atom| atoms.apply(atom) as u32);
        let action = Permutation::between(&renamed, &expected_frame).unwrap();
        assert_eq!([expected_frame[0], expected_frame[2]], expected_references);
        for coset in [0, 1] {
            assert_eq!(
                ClassKey::CisTrans.space().reindex(coset, action),
                Some(coset)
            );
        }
    }

    #[rstest]
    #[case::identity([0, 2, 3, 5], [0, 3], [0, 1])]
    #[case::first_reference([2, 0, 3, 5], [2, 3], [1, 0])]
    #[case::second_reference([0, 2, 5, 3], [0, 5], [1, 0])]
    #[case::both_references([2, 0, 5, 3], [2, 5], [0, 1])]
    #[case::endpoints([3, 5, 0, 2], [3, 0], [0, 1])]
    #[case::endpoints_first_reference([3, 5, 2, 0], [3, 2], [1, 0])]
    #[case::endpoints_second_reference([5, 3, 0, 2], [5, 0], [1, 0])]
    #[case::endpoints_both_references([5, 3, 2, 0], [5, 2], [0, 1])]
    fn test_bond_relation_frame(
        #[case] target: [u32; 4],
        #[case] references: [u32; 2],
        #[case] expected_cosets: [u32; 2],
        #[values((BondRelation::SameSide, 0), (BondRelation::OppositeSide, 1))] (relation, coset): (
            BondRelation,
            u32,
        ),
    ) {
        let arrangement = match relation {
            BondRelation::SameSide => [0, 2, 3, 5],
            BondRelation::OppositeSide => [0, 2, 5, 3],
        };
        assert_eq!(
            ClassKey::CisTrans
                .space()
                .index(Permutation::between(&arrangement, &[0, 2, 3, 5]).unwrap()),
            Some(coset),
        );
        let action = Permutation::between(&[0, 2, 3, 5], &target).unwrap();
        let expected_coset = expected_cosets[coset as usize];
        assert_eq!(
            ClassKey::CisTrans.space().reindex(coset, action),
            Some(expected_coset),
        );
        assert_eq!(
            StereoBondForm::new(StereoKind::CisTrans, coset).reframe_by(&action),
            Some(StereoBondForm::new(StereoKind::CisTrans, expected_coset)),
        );
        assert_eq!([target[0], target[2]], references);

        // Independent side labels for the four source ligands give the reference relation.
        let sides = match relation {
            BondRelation::SameSide => [1, -1, 1, -1],
            BondRelation::OppositeSide => [1, -1, -1, 1],
        };
        let same_side = sides[action.apply(0)] == sides[action.apply(2)];
        assert_eq!(expected_coset, if same_side { 0 } else { 1 });
    }

    #[rstest]
    #[case::cross_endpoint([0, 3, 2, 5])]
    #[case::interleaved([0, 3, 5, 2])]
    fn test_bond_relation_frame_error(#[case] target: [u32; 4], #[values(0, 1)] coset: u32) {
        let action = Permutation::between(&[0, 2, 3, 5], &target).unwrap();
        assert_eq!(ClassKey::CisTrans.space().reindex(coset, action), None);
        assert_eq!(
            StereoBondForm::new(StereoKind::CisTrans, coset).reframe_by(&action),
            None,
        );
    }
}
