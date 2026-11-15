//! Molecule views over owned or borrowed graph-symmetry data.

use std::borrow::Cow;

use super::super::id::AtomId;
use super::super::molecule::{GraphSymmetry, Molecule};

/// Graph-symmetry readings over owned or borrowed auxiliary data and a borrowed molecule.
///
/// The data must describe the molecule for meaningful results. Pairing them does not check
/// compatibility; mismatched context does not panic but has no specified semantic result.
///
/// # Semantic properties
///
/// Owned and borrowed views over compatible data give the same answers as the molecule
/// queries. Every query reads the supplied molecule and preserves its exact representation.
#[derive(Debug)]
pub struct GraphSymmetryView<'a> {
    molecule: &'a Molecule,
    symmetry: Cow<'a, GraphSymmetry>,
}

impl<'a> GraphSymmetryView<'a> {
    pub(crate) fn new(molecule: &'a Molecule, symmetry: GraphSymmetry) -> Self {
        Self {
            molecule,
            symmetry: Cow::Owned(symmetry),
        }
    }

    /// Borrow existing auxiliary data and its molecule without copying or recomputing either.
    ///
    /// The data must describe the supplied molecule. This relationship is not validated;
    /// mismatched context does not panic but has no specified semantic result.
    pub fn from_graph_symmetry(symmetry: &'a GraphSymmetry, molecule: &'a Molecule) -> Self {
        Self {
            molecule,
            symmetry: Cow::Borrowed(symmetry),
        }
    }

    /// Extract owned auxiliary data, cloning it when this view borrows it.
    pub fn into_graph_symmetry(self) -> GraphSymmetry {
        self.symmetry.into_owned()
    }

    /// Whether the molecule has no remapping to its mirror image.
    pub fn is_chiral(&self) -> bool {
        self.molecule.is_chiral(&self.symmetry)
    }

    /// Whether a preserving remapping sends `a` to `b`; false for unavailable ids.
    pub fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.molecule.same_proper_orbit(&self.symmetry, a, b)
    }

    /// Whether a preserving or mirror remapping sends `a` to `b`; false for unavailable ids.
    pub fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.molecule.same_star_orbit(&self.symmetry, a, b)
    }

    /// The atom-id-ordered proper orbit of `a`, empty for an unavailable id.
    pub fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId> {
        self.molecule.proper_orbit_of(&self.symmetry, a)
    }

    /// The atom-id-ordered proper-or-mirror orbit of `a`, empty for an unavailable id.
    pub fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId> {
        self.molecule.star_orbit_of(&self.symmetry, a)
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::ir::GraphSymmetryConfig;
    use crate::{mol_dsl, mol_dsl_concrete};

    #[fixture]
    fn molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{
            :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                           {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
        }"#
        )
    }

    #[rstest]
    #[case::different(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#))]
    #[case::larger(mol_dsl!(r#"{
        :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 8 "1"]]
    }"#))]
    #[case::smaller(Molecule::default())]
    fn test_graph_symmetry_view_context(molecule: Molecule, #[case] supplied: Molecule) {
        let symmetry = molecule
            .graph_symmetry(&GraphSymmetryConfig::default())
            .into_graph_symmetry();
        let view = GraphSymmetryView::from_graph_symmetry(&symmetry, &supplied);
        let _ = view.is_chiral();
        for atom in supplied.atoms().ids().chain([AtomId(9)]) {
            let _ = view.proper_orbit_of(atom);
            let _ = view.star_orbit_of(atom);
            let _ = view.same_proper_orbit(atom, AtomId(9));
            let _ = view.same_star_orbit(atom, AtomId(9));
        }
    }

    #[rstest]
    fn test_graph_symmetry_view_into_graph_symmetry(molecule: Molecule) {
        let symmetry = {
            let source = molecule.clone();
            source
                .graph_symmetry(&GraphSymmetryConfig::default())
                .into_graph_symmetry()
        };
        let borrowed = GraphSymmetryView::from_graph_symmetry(&symmetry, &molecule);
        let cloned = borrowed.into_graph_symmetry();
        let moved = molecule
            .graph_symmetry(&GraphSymmetryConfig::default())
            .into_graph_symmetry();
        for symmetry in [&symmetry, &cloned, &moved] {
            assert!(!molecule.is_chiral(symmetry));
            assert_eq!(
                molecule.proper_orbit_of(symmetry, AtomId(0)),
                vec![AtomId(0)]
            );
            assert_eq!(
                molecule.star_orbit_of(symmetry, AtomId(0)),
                vec![AtomId(0), AtomId(1)]
            );
        }
    }

    #[rstest]
    #[case::benzene(mol_dsl!(r#"{
        :atoms ["C" "C" "C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
    }"#))]
    fn test_graph_symmetry_view_is_chiral(#[case] input: Molecule) {
        assert!(!input
            .graph_symmetry(&GraphSymmetryConfig::default())
            .is_chiral());
    }
}
