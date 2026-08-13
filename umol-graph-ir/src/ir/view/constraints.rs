//! Constraint reading of an atom, reached by accessor chaining:
//! `molecule.atom(id).constraints()`.
//!
//! The bare read API is inherited from the stored container and keeps its
//! meanings — typed getters, `iter`, `is_empty` all read the asserted side,
//! which is what constraints store. The keyed core names both sides
//! explicitly: `asserted(key)` reads storage, `derived(key)` reads present
//! relations (vacuous on absence), and `derived_complete(key)` adds the
//! closure — under the caller's claim that the relation set is complete,
//! absence of a resolution-written overlay closes to its definite negative.
//! Topology keys read identically under both derived accessors; only the
//! overlay keys have an absence cell for the closure to fill. Mutation never
//! routes through this view; it belongs to the stored container.

use super::super::constraint::{
    AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
    MulticenterValenceForm,
};
use super::super::id::AtomId;
use super::super::molecule::Molecule;
use super::super::num::NumForm;
use super::super::ring::RingSet;
use super::super::stereo::TetrahedralStereoForm;
use super::super::traits::Lattice;
use super::atom::{asserted_constraints, derived_constraint};

/// Constraint reading of one atom: the asserted side under the container's
/// read API, and both sides under the keyed accessors.
#[derive(Clone, Copy, Debug)]
pub struct AtomConstraintsView<'a> {
    molecule: &'a Molecule,
    atom: AtomId,
    rings: Option<&'a RingSet>,
}

impl<'a> AtomConstraintsView<'a> {
    pub(crate) fn new(molecule: &'a Molecule, atom: AtomId) -> Self {
        Self {
            molecule,
            atom,
            rings: None,
        }
    }

    /// Attach ring context for the ring keys (`#R`, `#x`, `#y`).
    pub fn with_rings(mut self, rings: &'a RingSet) -> Self {
        self.rings = Some(rings);
        self
    }

    /// The stored side of `key`; absence is the vacuous constraint.
    pub fn asserted(&self, key: AtomConstraintKey) -> Option<&'a AtomConstraintForm> {
        asserted_constraints(self.molecule, self.atom).get(key)
    }

    /// The derived side of `key`, obtained by projection from present
    /// relations only; vacuous on absence.
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived(&self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        derived_constraint(self.molecule, self.atom, self.rings, key, false)
    }

    /// The derived side of `key` under the closure: absence of a
    /// resolution-written overlay yields its definite negative. Positive
    /// incidence and the topology keys agree with [`Self::derived`].
    ///
    /// # Panics
    ///
    /// A ring key without ring context ([`Self::with_rings`]) is a caller
    /// error — the caller scanning keys decides whether to build the ring
    /// set.
    pub fn derived_complete(&self, key: AtomConstraintKey) -> Option<AtomConstraintForm> {
        derived_constraint(self.molecule, self.atom, self.rings, key, true)
    }

    /// Whether this atom's constraint reading satisfies `pattern`: every
    /// pattern entry is refined by the meet of the asserted and
    /// [`Self::derived_complete`] sides at its key — the query-against-host
    /// reading. An internally conflicted key (the sides meet to `⊥`)
    /// satisfies nothing. Evaluation is driven by the pattern's keys; an
    /// empty pattern is satisfied.
    ///
    /// # Panics
    ///
    /// A ring key in `pattern` without ring context ([`Self::with_rings`]) is
    /// a caller error.
    pub fn satisfies(&self, pattern: &AtomConstraintsForm) -> bool {
        pattern.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived_complete(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => entry.as_undetermined(),
            };
            host.satisfies(entry)
        })
    }

    /// Whether `other` is compatible with this atom's constraint reading:
    /// for every key of `other`, a meet with the asserted and
    /// [`Self::derived`] sides exists — the narrowing-admissibility reading.
    /// A key on which this atom carries nothing constrains nothing; an
    /// internally conflicted key (the sides meet to `⊥`) is compatible with
    /// nothing. Evaluation is driven by `other`'s keys; an empty `other` is
    /// compatible.
    ///
    /// # Panics
    ///
    /// A ring key in `other` without ring context ([`Self::with_rings`]) is a
    /// caller error.
    pub fn is_compatible(&self, other: &AtomConstraintsForm) -> bool {
        other.iter().all(|entry| {
            let key = entry.key();
            let host = match (self.asserted(key), self.derived(key)) {
                (Some(asserted), Some(derived)) => match asserted.meet(&derived) {
                    Some(host) => host,
                    None => return false,
                },
                (Some(asserted), None) => asserted.clone(),
                (None, Some(derived)) => derived,
                (None, None) => return true,
            };
            entry.is_compatible(&host)
        })
    }

    // The stored container's read API, inherited with its meanings intact:
    // every accessor below reads the asserted side.

    pub fn valence(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).valence()
    }

    pub fn aromatic_valence(&self) -> Option<&'a AromaticValenceForm> {
        asserted_constraints(self.molecule, self.atom).aromatic_valence()
    }

    pub fn multicenter_valence(&self) -> Option<&'a MulticenterValenceForm> {
        asserted_constraints(self.molecule, self.atom).multicenter_valence()
    }

    pub fn tetrahedral_stereo(&self) -> Option<&'a TetrahedralStereoForm> {
        asserted_constraints(self.molecule, self.atom).tetrahedral_stereo()
    }

    pub fn degree(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).degree()
    }

    pub fn total_degree(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).total_degree()
    }

    pub fn total_valence(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).total_valence()
    }

    pub fn ring_degree(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).ring_degree()
    }

    pub fn ring_valence(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).ring_valence()
    }

    pub fn total_hydrogens(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).total_hydrogens()
    }

    pub fn donated_pairs(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).donated_pairs()
    }

    pub fn accepted_pairs(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).accepted_pairs()
    }

    pub fn ring_count(&self) -> Option<&'a NumForm> {
        asserted_constraints(self.molecule, self.atom).ring_count()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = &'a AtomConstraintForm> {
        asserted_constraints(self.molecule, self.atom).iter()
    }

    pub fn is_empty(&self) -> bool {
        asserted_constraints(self.molecule, self.atom).is_empty()
    }

    pub fn len(&self) -> usize {
        asserted_constraints(self.molecule, self.atom).len()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use crate::ir::aromatic::AromaticSystemForm;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::{
        AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
        MulticenterValenceForm, RingScope,
    };
    use crate::ir::dative::DativeBondForm;
    use crate::ir::id::AtomId;
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries};
    use crate::ir::num::NumForm;
    use crate::ir::ring::{RingConfig, RingModel};
    use crate::ir::stereo::{StereoAtomForm, StereoCoset, StereoKind, TetrahedralStereoForm};
    use crate::mol_dsl;

    #[rstest]
    #[case::present(
        Some(AtomConstraintForm::valence(4)),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(4))
    )]
    #[case::absent(None, AtomConstraintKey::Degree, None)]
    fn test_atom_constraints_view_asserted(
        #[case] set: Option<AtomConstraintForm>,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let mut atom = AtomForm::from_element(Element::C);
        if let Some(constraint) = set {
            atom.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().asserted(key),
            expected.as_ref()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintKey::Valence,
        Some(AtomConstraintForm::valence(2)),
    )]
    #[case::degree(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintKey::Degree,
        Some(AtomConstraintForm::degree(1)),
    )]
    #[case::total_degree(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalDegree,
        Some(AtomConstraintForm::total_degree(4)),
    )]
    #[case::total_hydrogens(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalHydrogens,
        Some(AtomConstraintForm::total_hydrogens(4)),
    )]
    #[case::total_valence(
        mol_dsl!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomConstraintKey::TotalValence,
        Some(AtomConstraintForm::total_valence(4)),
    )]
    #[case::donated_pairs_absent(
        mol_dsl!(r#"{:atoms ["N"] :bonds []}"#),
        AtomConstraintKey::DonatedPairs,
        None,
    )]
    #[case::donated_pairs_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::C),
            ],
            dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
            ..Default::default()
        }),
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(1)),
    )]
    #[case::aromatic_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::AromaticValence,
        None,
    )]
    #[case::aromatic_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            aromatic: vec![(vec![AtomId(0)], AromaticSystemForm::from_electrons(vec![1]))],
            ..Default::default()
        }),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::aromatic(NumForm::Lit(1)),
        )),
    )]
    #[case::aromatic_bond_flag(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(
            AromaticValenceForm::aromatic(NumForm::Undetermined),
        )),
    )]
    #[case::multicenter_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::MulticenterValence,
        None,
    )]
    #[case::tetrahedral_absent(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        AtomConstraintKey::TetrahedralStereo,
        None,
    )]
    #[case::tetrahedral_present(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 5],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        }),
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::stereo(StereoCoset::Lit(1)),
        )),
    )]
    fn test_atom_constraints_view_derived(
        #[case] molecule: Molecule,
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().derived(key),
            expected
        );
    }

    #[rstest]
    #[case::valence_agrees(AtomConstraintKey::Valence, Some(AtomConstraintForm::valence(0)))]
    #[case::donated_pairs_closed(
        AtomConstraintKey::DonatedPairs,
        Some(AtomConstraintForm::donated_pairs(0))
    )]
    #[case::accepted_pairs_closed(
        AtomConstraintKey::AcceptedPairs,
        Some(AtomConstraintForm::accepted_pairs(0))
    )]
    #[case::aromatic_closed(
        AtomConstraintKey::AromaticValence,
        Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic))
    )]
    #[case::multicenter_closed(
        AtomConstraintKey::MulticenterValence,
        Some(AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter,))
    )]
    #[case::tetrahedral_closed(
        AtomConstraintKey::TetrahedralStereo,
        Some(AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo))
    )]
    fn test_atom_constraints_view_derived_complete(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().derived_complete(key),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_derived_complete_bond_flag() {
        // The closure does not override positive Kekulé-bond evidence.
        let molecule = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#);
        assert_eq!(
            molecule
                .atom(AtomId(0))
                .constraints()
                .derived_complete(AtomConstraintKey::AromaticValence),
            Some(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::aromatic(NumForm::Undetermined),
            )),
        );
    }

    #[rstest]
    #[case::ring_degree(
        AtomConstraintKey::RingDegree,
        Some(AtomConstraintForm::ring_degree(2))
    )]
    #[case::ring_valence(
        AtomConstraintKey::RingValence,
        Some(AtomConstraintForm::ring_valence(2))
    )]
    #[case::membership_all(
        AtomConstraintKey::RingMembership(RingScope::All),
        Some(AtomConstraintForm::ring_membership(RingScope::All, 1))
    )]
    #[case::membership_size_match(
        AtomConstraintKey::RingMembership(RingScope::Size(6)),
        Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1))
    )]
    #[case::membership_size_no_match(
        AtomConstraintKey::RingMembership(RingScope::Size(5)),
        Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 0))
    )]
    fn test_atom_constraints_view_derived_ring(
        #[case] key: AtomConstraintKey,
        #[case] expected: Option<AtomConstraintForm>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let rings = molecule
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        assert_eq!(
            molecule
                .atom(AtomId(0))
                .constraints()
                .with_rings(&rings)
                .derived(key),
            expected
        );
    }

    #[rstest]
    #[should_panic(expected = "ring constraint key requires ring context")]
    fn test_atom_constraints_view_derived_error() {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        molecule
            .atom(AtomId(0))
            .constraints()
            .derived(AtomConstraintKey::RingDegree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::new(),
        true,
    )]
    #[case::valence_match(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        true,
    )]
    #[case::valence_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]),
        false,
    )]
    #[case::vacuous_entry(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(NumForm::Undetermined)]),
        true,
    )]
    #[case::closure_not_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        true,
    )]
    #[case::closure_aromatic_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::Lit(1))),
        ]),
        false,
    )]
    #[case::conflicted_host_mismatching_pattern(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        false,
    )]
    #[case::conflicted_host_matching_pattern(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        false,
    )]
    fn test_atom_constraints_view_satisfies(
        #[case] molecule: Molecule,
        #[case] pattern: AtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().satisfies(&pattern),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_satisfies_ring() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let rings = molecule
            .rings(RingModel::default(), RingConfig::default())
            .into_ring_set();
        let pattern = AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(
            RingScope::Size(6),
            1,
        )]);
        assert!(molecule
            .atom(AtomId(0))
            .constraints()
            .with_rings(&rings)
            .satisfies(&pattern));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::new(),
        true,
    )]
    #[case::shared_key_compatible(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(2)]),
        true,
    )]
    #[case::shared_key_conflict(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(3)]),
        false,
    )]
    #[case::absent_overlay_skipped(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        true,
    )]
    #[case::absent_stereo_skipped(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
        ]),
        true,
    )]
    #[case::asserted_conflict(
        mol_dsl!(r#"{:atoms ["C#a1"] :bonds []}"#),
        AtomConstraintsForm::from_iter([
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic),
        ]),
        false,
    )]
    #[case::conflicted_host(
        mol_dsl!(r#"{:atoms ["C#v4" "C"] :bonds [[0 1 "2"]]}"#),
        AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4)]),
        false,
    )]
    fn test_atom_constraints_view_is_compatible(
        #[case] molecule: Molecule,
        #[case] other: AtomConstraintsForm,
        #[case] expected: bool,
    ) {
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().is_compatible(&other),
            expected
        );
    }

    #[rstest]
    fn test_atom_constraints_view_inherited_getters() {
        // The container read API is inherited with its meanings intact: every
        // getter reads the asserted side.
        let mut atom = AtomForm::from_element(Element::C);
        for constraint in [
            AtomConstraintForm::valence(4),
            AtomConstraintForm::donated_pairs(1),
            AtomConstraintForm::accepted_pairs(2),
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::Lit(1))),
            AtomConstraintForm::multicenter_valence(MulticenterValenceForm::NotMulticenter),
            AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::NotStereo),
            AtomConstraintForm::degree(3),
            AtomConstraintForm::total_degree(4),
            AtomConstraintForm::total_valence(4),
            AtomConstraintForm::ring_degree(2),
            AtomConstraintForm::ring_valence(2),
            AtomConstraintForm::total_hydrogens(1),
            AtomConstraintForm::ring_membership(RingScope::All, 1),
        ] {
            atom.constraints.set(constraint);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom],
            ..Default::default()
        });
        let view = molecule.atom(AtomId(0)).constraints();
        assert_eq!(view.valence(), Some(&NumForm::Lit(4)));
        assert_eq!(view.donated_pairs(), Some(&NumForm::Lit(1)));
        assert_eq!(view.accepted_pairs(), Some(&NumForm::Lit(2)));
        assert_eq!(
            view.aromatic_valence(),
            Some(&AromaticValenceForm::Aromatic(NumForm::Lit(1)))
        );
        assert_eq!(
            view.multicenter_valence(),
            Some(&MulticenterValenceForm::NotMulticenter)
        );
        assert_eq!(
            view.tetrahedral_stereo(),
            Some(&TetrahedralStereoForm::NotStereo)
        );
        assert_eq!(view.degree(), Some(&NumForm::Lit(3)));
        assert_eq!(view.total_degree(), Some(&NumForm::Lit(4)));
        assert_eq!(view.total_valence(), Some(&NumForm::Lit(4)));
        assert_eq!(view.ring_degree(), Some(&NumForm::Lit(2)));
        assert_eq!(view.ring_valence(), Some(&NumForm::Lit(2)));
        assert_eq!(view.total_hydrogens(), Some(&NumForm::Lit(1)));
        assert_eq!(view.ring_count(), Some(&NumForm::Lit(1)));
        assert_eq!(view.len(), 13);
        assert!(!view.is_empty());
        assert_eq!(view.iter().count(), 13);
    }
}
