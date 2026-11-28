//! Substructure matching: occurrences of a pattern `Molecule` within a host
//! `Molecule`, each an injective pattern→host [`MoleculeCorrespondence`]. The receiver is the
//! pattern, so it parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the atom-bond topology then post-verifies
//! overlays; `Incidence` matches the overlay-augmented incidence graph for
//! connectivity carried by overlays.

mod prepared;

use std::ops::ControlFlow;

pub use prepared::{PreparedSubstructureHost, PreparedSubstructurePattern};
use umol_graph_core::{
    AutomorphismAlgorithm, Correspondence, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};
use umol_perm::{DynPermutation, Permutation};

use super::atom::AtomForm;
use super::bond::BondForm;
use super::correspondence::{
    induced_aromatic_systems, induced_bonds, induced_dative_bonds, induced_multicenter_bonds,
    induced_noncovalent_bonds, map_atom, map_ligands, MoleculeCorrespondence,
};
use super::electrons::ElectronCountsForm;
use super::id::AtomId;
use super::molecule::Molecule;
use super::traits::{FrameTransport, Lattice};

/// Strategy for `substructure_matches`, each composing over a
/// [`SubgraphIsomorphismAlgorithm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubstructureMatchAlgorithm {
    /// Match the atom-bond topology, then post-verify the N-ary / special
    /// overlays against the atom correspondence.
    GraphAndOverlays,
    /// Match the overlay-augmented incidence graph (localized bonds as direct
    /// edges, overlays as pseudonodes) — useful when pattern connectivity is
    /// carried by overlays.
    Incidence,
}

/// Algorithms and symmetry-reduction options used to enumerate substructure matches.
///
/// This type deliberately has no default at the graph-IR layer: every graph
/// algorithm selection remains explicit at the call site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureMatchConfig {
    /// Strategy used to represent and match molecule entities.
    pub match_algorithm: SubstructureMatchAlgorithm,
    /// Algorithm used to enumerate embeddings of the selected graph representation.
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    /// Algorithm used to derive relevant-ring constraints requested by the pattern.
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    /// Algorithm used to compute proper symmetry groups when pattern or host
    /// symmetry reduction is enabled.
    pub automorphism_algorithm: AutomorphismAlgorithm,
    /// Emit one canonical representative per pattern-symmetry orbit instead of
    /// all pattern-symmetric embeddings.
    pub reduce_pattern_symmetry: bool,
    /// Emit one canonical representative per host-symmetry orbit instead of
    /// all host-symmetric embeddings.
    pub reduce_host_symmetry: bool,
}

impl Molecule {
    /// Visits each occurrence of `self` (the pattern) within `host` as an injective
    /// pattern→host [`MoleculeCorrespondence`] until traversal completes or the
    /// visitor returns [`ControlFlow::Break`]. Pattern entity predicates are
    /// evaluated per key against the host's constraint reading under the closure
    /// (`constraints().satisfies`); ring keys use the fixed Relevant projection
    /// through size 22. Traversal is deterministic for a fixed representation,
    /// but its order is not a canonical ordering contract. A matched aromatic or
    /// multicenter overlay whose literal electron-count length disagrees with its
    /// atom list rejects that candidate, on either the pattern or host side.
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    ///
    /// # Semantic properties
    ///
    /// - [`SubstructureMatchAlgorithm::GraphAndOverlays`] and
    ///   [`SubstructureMatchAlgorithm::Incidence`] emit the same match set (and, when symmetry
    ///   reduction is enabled, the same proper symmetry-orbit set and count).
    /// - Within a strategy, the emitted match set (including canonical symmetry-orbit
    ///   representatives) is invariant across [`SubgraphIsomorphismAlgorithm`] and
    ///   [`AutomorphismAlgorithm`] variants.
    /// - [`Self::substructure_matches`], [`Self::visit_substructure_embeddings`],
    ///   [`Self::has_substructure_match`], [`Self::first_substructure_match`],
    ///   [`Self::count_substructure_matches`], and [`PreparedSubstructurePattern`] agree with
    ///   this operation for every `(pattern, host, config)`.
    pub fn visit_substructure_matches<B, F>(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
        visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    {
        self.assert_matchable_pattern();
        if !passes_quick_molecule_prescreen(self, host) {
            return ControlFlow::Continue(());
        }
        PreparedSubstructurePattern::new(self, config).visit_matches(host, visitor)
    }

    /// Visits each occurrence of `self` (the pattern) within `host` as a borrowed
    /// slice of matched host [`AtomId`]s indexed by pattern [`AtomId`], without
    /// materializing a [`MoleculeCorrespondence`], until traversal completes or the
    /// visitor returns [`ControlFlow::Break`].
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn visit_substructure_embeddings<B, F>(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
        visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(&[AtomId]) -> ControlFlow<B>,
    {
        self.assert_matchable_pattern();
        if !passes_quick_molecule_prescreen(self, host) {
            return ControlFlow::Continue(());
        }
        PreparedSubstructurePattern::new(self, config).visit_embeddings(host, visitor)
    }

    /// Occurrences of `self` (the pattern) within `host`, collected from
    /// [`Molecule::visit_substructure_matches`].
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn substructure_matches(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
    ) -> Vec<MoleculeCorrespondence> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> = self.visit_substructure_matches(host, config, |correspondence| {
            occurrences.push(correspondence);
            ControlFlow::Continue(())
        });
        occurrences
    }

    /// Returns whether `self` (the pattern) has at least one occurrence within `host`.
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn has_substructure_match(&self, host: &Molecule, config: SubstructureMatchConfig) -> bool {
        self.assert_matchable_pattern();
        if !passes_quick_molecule_prescreen(self, host) {
            return false;
        }
        PreparedSubstructurePattern::new(self, config).has_match(host)
    }

    /// Returns the first occurrence of `self` (the pattern) within `host` as a
    /// [`MoleculeCorrespondence`], or `None` if no match exists.
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn first_substructure_match(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
    ) -> Option<MoleculeCorrespondence> {
        self.assert_matchable_pattern();
        if !passes_quick_molecule_prescreen(self, host) {
            return None;
        }
        PreparedSubstructurePattern::new(self, config).first_match(host)
    }

    /// Counts every occurrence of `self` (the pattern) within `host` without
    /// materializing correspondences.
    ///
    /// # Panics
    ///
    /// Panics if `self` carries molecule constraints, which substructure matching
    /// does not evaluate.
    pub fn count_substructure_matches(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
    ) -> u64 {
        self.assert_matchable_pattern();
        if !passes_quick_molecule_prescreen(self, host) {
            return 0;
        }
        PreparedSubstructurePattern::new(self, config).count_matches(host)
    }

    fn assert_matchable_pattern(&self) {
        assert!(
            self.constraints().is_empty(),
            "substructure matching does not evaluate molecule constraints"
        );
    }

    /// Post-verify a topology occurrence's overlays against the atom correspondence, returning the
    /// injective pattern→host [`MoleculeCorrespondence`] or `None` if any pattern overlay has no
    /// matching host overlay. Each N-ary / special overlay is located by **exact participant set**
    /// via the per-kind inducer (which already checks dative donor/acceptor roles); the pattern
    /// overlay's predicate is then required to match the located host overlay's, and every pattern
    /// overlay must be matched. Stereo overlays are matched by the bespoke coset filter.
    fn verify_overlays(
        &self,
        host: &Molecule,
        atoms: Correspondence<AtomId>,
    ) -> Option<MoleculeCorrespondence> {
        let pattern = self;
        let bonds = induced_bonds(pattern, host, &atoms)?;

        let dative_bonds = if pattern.dative_bonds().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.dative_bonds().count())
                .expect("empty correspondence is valid")
        } else if host.dative_bonds().count() < pattern.dative_bonds().count() {
            return None;
        } else {
            let dative_bonds = induced_dative_bonds(pattern, host, &atoms)?;
            if dative_bonds.matched_pair_count() != pattern.dative_bonds().count() {
                return None;
            }
            for &(p, h) in dative_bonds.matched_pairs() {
                if !pattern
                    .dative_bond(p)
                    .attributes()
                    .matches(host.dative_bond(h).attributes())
                {
                    return None;
                }
            }
            dative_bonds
        };

        let aromatic_systems = if pattern.aromatic_systems().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.aromatic_systems().count())
                .expect("empty correspondence is valid")
        } else if host.aromatic_systems().count() < pattern.aromatic_systems().count() {
            return None;
        } else {
            let aromatic_systems = induced_aromatic_systems(pattern, host, &atoms)?;
            if aromatic_systems.matched_pair_count() != pattern.aromatic_systems().count() {
                return None;
            }
            for &(p, h) in aromatic_systems.matched_pairs() {
                let p_view = pattern.aromatic_system(p);
                let h_view = host.aromatic_system(h);
                if matches!(p_view.electrons(), ElectronCountsForm::Lit(counts)
                    if counts.len() != p_view.atom_ids().len())
                    || matches!(h_view.electrons(), ElectronCountsForm::Lit(counts)
                        if counts.len() != h_view.atom_ids().len())
                {
                    return None;
                }
                let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
                // The host frame named in the pattern's own atom ids, so the pattern form can be
                // restated into it before the two are compared.
                let host_frame: Vec<AtomId> = h_view
                    .atom_ids()
                    .map(|host_atom| {
                        atoms
                            .left_of(host_atom)
                            .expect("host overlay atom is matched")
                    })
                    .collect();
                let action = DynPermutation::between(&pat_atoms, &host_frame)?;
                if !p_view
                    .attributes()
                    .clone()
                    .reframe_by(&action)?
                    .matches(h_view.attributes())
                {
                    return None;
                }
            }
            aromatic_systems
        };

        let multicenter_bonds = if pattern.multicenter_bonds().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.multicenter_bonds().count())
                .expect("empty correspondence is valid")
        } else if host.multicenter_bonds().count() < pattern.multicenter_bonds().count() {
            return None;
        } else {
            let multicenter_bonds = induced_multicenter_bonds(pattern, host, &atoms)?;
            if multicenter_bonds.matched_pair_count() != pattern.multicenter_bonds().count() {
                return None;
            }
            for &(p, h) in multicenter_bonds.matched_pairs() {
                let p_view = pattern.multicenter_bond(p);
                let h_view = host.multicenter_bond(h);
                if matches!(p_view.electrons(), ElectronCountsForm::Lit(counts)
                    if counts.len() != p_view.atom_ids().len())
                    || matches!(h_view.electrons(), ElectronCountsForm::Lit(counts)
                        if counts.len() != h_view.atom_ids().len())
                {
                    return None;
                }
                let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
                // The host frame named in the pattern's own atom ids, so the pattern form can be
                // restated into it before the two are compared.
                let host_frame: Vec<AtomId> = h_view
                    .atom_ids()
                    .map(|host_atom| {
                        atoms
                            .left_of(host_atom)
                            .expect("host overlay atom is matched")
                    })
                    .collect();
                let action = DynPermutation::between(&pat_atoms, &host_frame)?;
                if !p_view
                    .attributes()
                    .clone()
                    .reframe_by(&action)?
                    .matches(h_view.attributes())
                {
                    return None;
                }
            }
            multicenter_bonds
        };

        let noncovalent_bonds = if pattern.noncovalent_bonds().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.noncovalent_bonds().count())
                .expect("empty correspondence is valid")
        } else if host.noncovalent_bonds().count() < pattern.noncovalent_bonds().count() {
            return None;
        } else {
            let noncovalent_bonds = induced_noncovalent_bonds(pattern, host, &atoms)?;
            if noncovalent_bonds.matched_pair_count() != pattern.noncovalent_bonds().count() {
                return None;
            }
            for &(p, h) in noncovalent_bonds.matched_pairs() {
                if !pattern
                    .noncovalent_bond(p)
                    .attributes()
                    .matches(host.noncovalent_bond(h).attributes())
                {
                    return None;
                }
            }
            noncovalent_bonds
        };

        // Stereo: a pattern stereo overlay matches iff the corresponding host site bears a stereo
        // element whose form, restated into the pattern's ligand frame (via the atom
        // correspondence), the pattern form admits. The comparison is over the whole form, so the
        // kind, the configuration and every frame-relative constraint are one question; an
        // `Undetermined` pattern configuration admits both handednesses.
        let stereo_atoms = if pattern.stereo_atoms().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.stereo_atoms().count())
                .expect("empty correspondence is valid")
        } else if host.stereo_atoms().count() < pattern.stereo_atoms().count() {
            return None;
        } else {
            let mut stereo_atom = Vec::with_capacity(pattern.stereo_atoms().count());
            for sp in pattern.stereo_atoms().iter() {
                let host_site =
                    map_atom(&atoms, sp.site_id()).expect("a matched pattern atom is matched");
                let sh = host.atom(host_site).stereo_atom()?;
                let frame = map_ligands(&atoms, sp.ligand_ids())
                    .expect("matched pattern ligands are matched");
                let action = Permutation::between(sh.ligand_ids(), &frame)?;
                let restated = sh.attributes().clone().reframe_by(&action)?;
                sp.attributes().matches(&restated).then_some(())?;
                stereo_atom.push((sp.id(), sh.id()));
            }
            Correspondence::new(
                stereo_atom,
                pattern.stereo_atoms().count(),
                host.stereo_atoms().count(),
            )
            .ok()?
        };

        let stereo_bonds = if pattern.stereo_bonds().count() == 0 {
            Correspondence::new(Vec::new(), 0, host.stereo_bonds().count())
                .expect("empty correspondence is valid")
        } else if host.stereo_bonds().count() < pattern.stereo_bonds().count() {
            return None;
        } else {
            let mut stereo_bond = Vec::with_capacity(pattern.stereo_bonds().count());
            for sp in pattern.stereo_bonds().iter() {
                let host_site = bonds
                    .right_of(sp.site_id())
                    .expect("a matched pattern bond is matched");
                let sh = host.bond(host_site).stereo_bond()?;
                let frame = map_ligands(&atoms, sp.ligand_ids())
                    .expect("matched pattern ligands are matched");
                let action = Permutation::between(sh.ligand_ids(), &frame)?;
                let restated = sh.attributes().clone().reframe_by(&action)?;
                sp.attributes().matches(&restated).then_some(())?;
                stereo_bond.push((sp.id(), sh.id()));
            }
            Correspondence::new(
                stereo_bond,
                pattern.stereo_bonds().count(),
                host.stereo_bonds().count(),
            )
            .ok()?
        };

        Some(MoleculeCorrespondence::new(
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        ))
    }
}

/// Cheap necessary cardinality, max-degree, and overlay-count check before preparing a
/// pattern for a one-shot `Molecule` match call.
fn passes_quick_molecule_prescreen(pattern: &Molecule, host: &Molecule) -> bool {
    pattern.atoms().count() <= host.atoms().count()
        && pattern.bonds().count() <= host.bonds().count()
        && pattern.raw_graph().max_degree() <= host.raw_graph().max_degree()
        && pattern.dative_bonds().count() <= host.dative_bonds().count()
        && pattern.aromatic_systems().count() <= host.aromatic_systems().count()
        && pattern.multicenter_bonds().count() <= host.multicenter_bonds().count()
        && pattern.noncovalent_bonds().count() <= host.noncovalent_bonds().count()
        && pattern.stereo_atoms().count() <= host.stereo_atoms().count()
        && pattern.stereo_bonds().count() <= host.stereo_bonds().count()
}

/// Inherent-field match: each pattern field admits the host's. Constraints are
/// compared separately through the host's constraint reading; the exhaustive
/// destructure makes a new `AtomForm` field a compile error here.
fn atom_fields_match(pattern: &AtomForm, host: &AtomForm) -> bool {
    let AtomForm {
        element,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired_electrons,
        constraints: _,
    } = pattern;
    element.matches(&host.element)
        && isotope_mass.matches(&host.isotope_mass)
        && charge.matches(&host.charge)
        && implicit_hydrogens.matches(&host.implicit_hydrogens)
        && lone_pairs.matches(&host.lone_pairs)
        && unpaired_electrons.matches(&host.unpaired_electrons)
}

/// Inherent-field match for bonds; see [`atom_fields_match`].
fn bond_fields_match(pattern: &BondForm, host: &BondForm) -> bool {
    let BondForm {
        order,
        charge,
        unpaired_electrons,
        constraints: _,
    } = pattern;
    order.matches(&host.order)
        && charge.matches(&host.charge)
        && unpaired_electrons.matches(&host.unpaired_electrons)
}

/// `pattern_form` matches `host_form` for an overlay whose payload is position-indexed by member
/// (aromatic / multicenter electron counts). The two overlays store their members in their own
/// participant order and `matches` compares the count vector whole, so the pattern payload is first
/// reindexed into the host's member order (via the atom correspondence) with
/// [`RelationData::on_permutation`].
#[cfg(test)]
mod tests {
    use std::ops::ControlFlow;

    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
    use umol_graph_core::{
        AutomorphismAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
    };

    use super::super::id::AtomId;
    use super::super::molecule::Molecule;
    use super::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
    use super::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
    use crate::ir::{
        AromaticSystemForm, AtomForm, ElectronCountsForm, EntityKind, MoleculeEntries,
        MulticenterBondForm,
    };
    use crate::mol_dsl;

    const SUBISO_ALGS: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];

    const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

    #[rstest]
    #[case::topology(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::noncovalent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]] :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#),
        vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]]
    )]
    #[case::no_match(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]]}"#),
        vec![]
    )]
    fn test_molecule_visit_substructure_matches(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let mut occurrences: Vec<Vec<AtomId>> = Vec::new();
                let flow: ControlFlow<()> =
                    pattern.visit_substructure_matches(&host, config, |correspondence| {
                        occurrences.push(
                            correspondence
                                .atoms()
                                .matched_pairs()
                                .iter()
                                .map(|&(_, host)| host)
                                .collect(),
                        );
                        ControlFlow::Continue(())
                    });
                assert_eq!(flow, ControlFlow::Continue(()), "{strategy:?}/{subiso:?}");
                occurrences.sort();
                assert_eq!(occurrences, expected, "{strategy:?}/{subiso:?}");
            }
        }
    }

    #[rstest]
    #[case::topology(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    fn test_molecule_visit_substructure_matches_termination(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let first = pattern.visit_substructure_matches(&host, config, |correspondence| {
                    ControlFlow::Break(
                        correspondence
                            .atoms()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, host)| host)
                            .collect::<Vec<_>>(),
                    )
                });
                let ControlFlow::Break(atoms) = first else {
                    panic!("expected Break on first occurrence: {strategy:?}/{subiso:?}");
                };
                assert!(
                    expected.contains(&atoms),
                    "{strategy:?}/{subiso:?}: invalid occurrence {atoms:?}"
                );
            }
        }
    }

    #[rstest]
    #[case::one_molecule_constraint(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#)
    )]
    #[case::two_molecule_constraints(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}} {:charge-sum {:sum 0}}]}"#)
    )]
    #[should_panic(expected = "substructure matching does not evaluate molecule constraints")]
    fn test_molecule_visit_substructure_matches_error(
        #[case] pattern: Molecule,
        #[values(GraphAndOverlays, Incidence)] strategy: SubstructureMatchAlgorithm,
    ) {
        let host = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let config = SubstructureMatchConfig {
            match_algorithm: strategy,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        let _: ControlFlow<()> =
            pattern.visit_substructure_matches(&host, config, |_| ControlFlow::Continue(()));
    }

    #[rstest]
    #[case::topology(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::atom_relevant_ring_membership(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#R3"] :bonds []}"#),
        vec![
            vec![AtomId(0)],
            vec![AtomId(1)],
            vec![AtomId(2)],
            vec![AtomId(3)],
        ]
    )]
    #[case::atom_simple_ring_membership(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#R6"] :bonds []}"#),
        vec![]
    )]
    #[case::atom_ring_degree_and_valence(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#x3#y3"] :bonds []}"#),
        vec![
            vec![AtomId(0)],
            vec![AtomId(1)],
            vec![AtomId(2)],
            vec![AtomId(3)],
        ]
    )]
    #[case::bond_relevant_ring_membership(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#R2"]]}"#),
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(0), AtomId(2)],
            vec![AtomId(0), AtomId(3)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(1), AtomId(3)],
            vec![AtomId(2), AtomId(0)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
            vec![AtomId(3), AtomId(0)],
            vec![AtomId(3), AtomId(1)],
            vec![AtomId(3), AtomId(2)],
        ]
    )]
    #[case::bond_simple_ring_membership(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#R4"]]}"#),
        vec![]
    )]
    #[case::atom_valence(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "2"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#v3"] :bonds []}"#),
        vec![vec![AtomId(1)]]
    )]
    #[case::atom_degree(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#D2"] :bonds []}"#),
        vec![vec![AtomId(1)]]
    )]
    #[case::atom_total_degree(
        mol_dsl!(r#"{:atoms ["C #h3" "C #h1"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#X4"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_total_hydrogens(
        mol_dsl!(r#"{:atoms ["C #h1" "H" "C #h0"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#H2"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_total_valence(
        mol_dsl!(r#"{:atoms ["C #h1" "C #h3"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#V4"] :bonds []}"#),
        vec![vec![AtomId(1)]]
    )]
    #[case::atom_aromatic_valence(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        mol_dsl!(r#"{:atoms ["C#a2"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_not_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#a!"] :bonds []}"#),
        vec![vec![AtomId(0)], vec![AtomId(1)]]
    )]
    #[case::atom_aromatic_bond_evidence(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        mol_dsl!(r#"{:atoms ["C#a+"] :bonds []}"#),
        vec![vec![AtomId(0)], vec![AtomId(1)]]
    )]
    #[case::atom_not_aromatic_bond_evidence(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        mol_dsl!(r#"{:atoms ["C#a!"] :bonds []}"#),
        vec![]
    )]
    #[case::atom_multicenter_valence(
        mol_dsl!(r#"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#),
        mol_dsl!(r#"{:atoms ["B#m2"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_not_multicenter(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#m!"] :bonds []}"#),
        vec![vec![AtomId(0)], vec![AtomId(1)]]
    )]
    #[case::atom_donated_pairs(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N#d1"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_accepted_pairs(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["B#t1"] :bonds []}"#),
        vec![vec![AtomId(1)]]
    )]
    #[case::atom_stereo(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C#T+"] :bonds []}"#),
        vec![vec![AtomId(0)]]
    )]
    #[case::atom_not_stereo(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#T!"] :bonds []}"#),
        vec![vec![AtomId(0)], vec![AtomId(1)]]
    )]
    #[case::bond_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,2]"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*#a"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::bond_aromatic_asserted_without_overlay(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*#a"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::bond_not_cis_trans(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#C!"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::noncovalent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]] :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#),
        vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]]
    )]
    #[case::noncovalent_absent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#),
        vec![]
    )]
    #[case::dative(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        vec![vec![AtomId(0), AtomId(1)]]
    )]
    #[case::dative_roles_swapped(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        vec![]
    )]
    #[case::stereo_chiral(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_enantiomer(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        vec![]
    )]
    #[case::stereo_agnostic_in_r(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_agnostic_in_s(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_reframed(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        vec![vec![AtomId(0), AtomId(3), AtomId(2), AtomId(1)]]
    )]
    #[case::stereo_reframed_enantiomer(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        vec![]
    )]
    #[case::stereo_constraint_agrees(
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)="}]}"##),
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)="}]}"##),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_constraint_contradicts(
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)/"}]}"##),
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)="}]}"##),
        vec![]
    )]
    #[case::stereo_constraint_absent_in_host(
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"##),
        mol_dsl!(r##"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)="}]}"##),
        vec![]
    )]
    #[case::stereo_bond(
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"] [2 0 "1"] [2 1 "1"] [3 4 "1"] [3 5 "1"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :attrs "Ct1"}]}"#),
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"] [2 0 "1"] [2 1 "1"] [3 4 "1"] [3 5 "1"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :attrs "Ct1"}]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)]]
    )]
    fn test_molecule_substructure_matches(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let mut occurrences: Vec<Vec<AtomId>> = pattern
                    .substructure_matches(
                        &host,
                        SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                            reduce_pattern_symmetry: false,
                            reduce_host_symmetry: false,
                        },
                    )
                    .iter()
                    .map(|c| {
                        c.atoms()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, host)| host)
                            .collect()
                    })
                    .collect();
                occurrences.sort();
                assert_eq!(occurrences, expected, "{strategy:?}/{subiso:?}");
            }
        }
    }

    #[rstest]
    #[case::ordered(
        [Element::C, Element::N, Element::O],
        vec![AtomId(0), AtomId(1), AtomId(2)],
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)]],
    )]
    #[case::reordered(
        [Element::C, Element::N, Element::O],
        vec![AtomId(2), AtomId(0), AtomId(1)],
        ElectronCountsForm::Lit(vec![0, 2, 1]),
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)]],
    )]
    #[case::different_association(
        [Element::C, Element::N, Element::O],
        vec![AtomId(2), AtomId(0), AtomId(1)],
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        vec![],
    )]
    #[case::undetermined_pattern(
        [Element::C, Element::N, Element::O],
        vec![AtomId(2), AtomId(0), AtomId(1)],
        ElectronCountsForm::Undetermined,
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)]],
    )]
    #[case::undetermined_host(
        [Element::C, Element::N, Element::O],
        vec![AtomId(2), AtomId(0), AtomId(1)],
        ElectronCountsForm::Lit(vec![0, 2, 1]),
        ElectronCountsForm::Undetermined,
        vec![],
    )]
    #[case::symmetric_nonuniform(
        [Element::C; 3],
        vec![AtomId(2), AtomId(0), AtomId(1)],
        ElectronCountsForm::Lit(vec![0, 2, 1]),
        ElectronCountsForm::Lit(vec![2, 1, 0]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)]],
    )]
    #[case::symmetric_repeated_counts(
        [Element::C; 3],
        vec![AtomId(0), AtomId(1), AtomId(2)],
        ElectronCountsForm::Lit(vec![1, 1, 0]),
        ElectronCountsForm::Lit(vec![1, 1, 0]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)],
             vec![AtomId(1), AtomId(0), AtomId(2)]],
    )]
    #[case::symmetric_uniform(
        [Element::C; 3],
        vec![AtomId(0), AtomId(1), AtomId(2)],
        ElectronCountsForm::Lit(vec![1, 1, 1]),
        ElectronCountsForm::Lit(vec![1, 1, 1]),
        vec![vec![AtomId(0), AtomId(1), AtomId(2)],
             vec![AtomId(0), AtomId(2), AtomId(1)],
             vec![AtomId(1), AtomId(0), AtomId(2)],
             vec![AtomId(1), AtomId(2), AtomId(0)],
             vec![AtomId(2), AtomId(0), AtomId(1)],
             vec![AtomId(2), AtomId(1), AtomId(0)]],
    )]
    #[case::symmetric_undetermined(
        [Element::C; 3],
        vec![AtomId(0), AtomId(1), AtomId(2)],
        ElectronCountsForm::Undetermined,
        ElectronCountsForm::Undetermined,
        vec![vec![AtomId(0), AtomId(1), AtomId(2)],
             vec![AtomId(0), AtomId(2), AtomId(1)],
             vec![AtomId(1), AtomId(0), AtomId(2)],
             vec![AtomId(1), AtomId(2), AtomId(0)],
             vec![AtomId(2), AtomId(0), AtomId(1)],
             vec![AtomId(2), AtomId(1), AtomId(0)]],
    )]
    fn test_molecule_substructure_matches_electrons(
        #[values(EntityKind::AromaticSystem, EntityKind::MulticenterBond)] kind: EntityKind,
        #[case] elements: [Element; 3],
        #[case] atoms: Vec<AtomId>,
        #[case] pattern_counts: ElectronCountsForm,
        #[case] host_counts: ElectronCountsForm,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let [pattern, host] = [
            (atoms, pattern_counts),
            (vec![AtomId(0), AtomId(1), AtomId(2)], host_counts),
        ]
        .map(|(atoms, electrons)| {
            let mut entries = MoleculeEntries {
                atoms: elements.map(AtomForm::from_element).to_vec(),
                ..Default::default()
            };
            match kind {
                EntityKind::AromaticSystem => {
                    entries
                        .aromatic
                        .push((atoms, AromaticSystemForm::new(electrons)));
                }
                EntityKind::MulticenterBond => {
                    entries
                        .multicenter
                        .push((atoms, MulticenterBondForm::new(electrons)));
                }
                _ => unreachable!(),
            }
            Molecule::from_entries(entries)
        });
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let mut occurrences: Vec<Vec<AtomId>> = pattern
                    .substructure_matches(
                        &host,
                        SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                            reduce_pattern_symmetry: false,
                            reduce_host_symmetry: false,
                        },
                    )
                    .iter()
                    .map(|correspondence| {
                        correspondence
                            .atoms()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, host)| host)
                            .collect()
                    })
                    .collect();
                occurrences.sort();
                assert_eq!(occurrences, expected, "{strategy:?}/{subiso:?}");
            }
        }
    }

    #[rstest]
    #[case::overlay_free_pattern_on_overlay_host(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "N" "O"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}] :dative-bonds [{:acceptor 3 :donors [4] :attrs "1"}] :noncovalent-bonds [{:atoms [0 4] :attrs "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        6,
        (0, 0, 0, 0, 0, 0),
        (1, 1, 0, 1, 0, 0),
    )]
    #[case::overlay_pattern_on_matching_overlay_host(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "N" "O"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}] :dative-bonds [{:acceptor 3 :donors [4] :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]}"#),
        6,
        (0, 1, 0, 0, 0, 0),
        (1, 1, 0, 0, 0, 0),
    )]
    #[case::overlay_pattern_on_overlay_free_host(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]}"#),
        0,
        (0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0),
    )]
    #[case::overlay_free_pattern_on_stereo_and_multicenter_host(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br" "B" "B" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :multicenter-bonds [{:atoms [4 5 6] :attrs "[2,0,0]"}] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "F"] :bonds [[0 1 "1"]]}"#),
        1,
        (0, 0, 0, 0, 0, 0),
        (0, 0, 1, 0, 1, 0),
    )]
    fn test_molecule_substructure_matches_overlays(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected_count: usize,
        #[case] expected_matched_overlays: (usize, usize, usize, usize, usize, usize),
        #[case] expected_right_counts: (usize, usize, usize, usize, usize, usize),
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let matches = pattern.substructure_matches(
                    &host,
                    SubstructureMatchConfig {
                        match_algorithm: strategy,
                        subgraph_isomorphism_algorithm: subiso,
                        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                        automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                        reduce_pattern_symmetry: false,
                        reduce_host_symmetry: false,
                    },
                );
                assert_eq!(matches.len(), expected_count, "{strategy:?}/{subiso:?}");
                for correspondence in &matches {
                    assert_eq!(
                        correspondence.bonds().matched_pair_count(),
                        pattern.bonds().count(),
                    );
                    assert_eq!(correspondence.bonds().left_count(), pattern.bonds().count());
                    assert_eq!(correspondence.bonds().right_count(), host.bonds().count());
                    assert_eq!(
                        (
                            correspondence.dative_bonds().matched_pair_count(),
                            correspondence.aromatic_systems().matched_pair_count(),
                            correspondence.multicenter_bonds().matched_pair_count(),
                            correspondence.noncovalent_bonds().matched_pair_count(),
                            correspondence.stereo_atoms().matched_pair_count(),
                            correspondence.stereo_bonds().matched_pair_count(),
                        ),
                        expected_matched_overlays,
                    );
                    assert_eq!(
                        (
                            correspondence.dative_bonds().right_count(),
                            correspondence.aromatic_systems().right_count(),
                            correspondence.multicenter_bonds().right_count(),
                            correspondence.noncovalent_bonds().right_count(),
                            correspondence.stereo_atoms().right_count(),
                            correspondence.stereo_bonds().right_count(),
                        ),
                        expected_right_counts,
                    );
                }
            }
        }
    }

    #[rstest]
    #[case::one_molecule_constraint(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#)
    )]
    #[should_panic(expected = "substructure matching does not evaluate molecule constraints")]
    fn test_molecule_substructure_matches_error(#[case] pattern: Molecule) {
        let host = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let config = SubstructureMatchConfig {
            match_algorithm: GraphAndOverlays,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        let _ = pattern.substructure_matches(&host, config);
    }

    #[rstest]
    #[case::topology(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::overlay_aromatic(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![AtomId(0), AtomId(2), AtomId(1)],
        ]
    )]
    #[case::no_match_larger_pattern(
        mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        vec![]
    )]
    fn test_molecule_visit_substructure_embeddings(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let config = SubstructureMatchConfig {
                    match_algorithm: strategy,
                    subgraph_isomorphism_algorithm: subiso,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
                    reduce_pattern_symmetry: false,
                    reduce_host_symmetry: false,
                };
                let mut visited = Vec::new();
                let flow: ControlFlow<()> =
                    pattern.visit_substructure_embeddings(&host, config, |embedding| {
                        visited.push(embedding.to_vec());
                        ControlFlow::Continue(())
                    });
                assert_eq!(flow, ControlFlow::Continue(()), "{strategy:?}/{subiso:?}");
                visited.sort();
                assert_eq!(visited, expected, "{strategy:?}/{subiso:?}");

                if !expected.is_empty() {
                    let early = pattern.visit_substructure_embeddings(&host, config, |embedding| {
                        ControlFlow::Break(embedding.to_vec())
                    });
                    let ControlFlow::Break(first) = early else {
                        panic!("expected Break on first embedding: {strategy:?}/{subiso:?}");
                    };
                    assert!(expected.contains(&first), "{strategy:?}/{subiso:?}");
                }
            }
        }
    }

    #[rstest]
    #[case::topology_match(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        (2, 1, 2, 1)
    )]
    #[case::cyclohexane_c3_symmetry(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        (12, 6, 1, 1)
    )]
    #[case::stereo_enantiomer_no_match(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]}"#),
        (0, 0, 0, 0)
    )]
    #[case::tetrahedral_methane_th0(
        mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
        (12, 1, 1, 1)
    )]
    #[case::larger_pattern_rejected(
        mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        (0, 0, 0, 0)
    )]
    fn test_molecule_substructure_match_modes(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] expected_counts: (u64, u64, u64, u64),
    ) {
        let modes = [
            (false, false, expected_counts.0),
            (true, false, expected_counts.1),
            (false, true, expected_counts.2),
            (true, true, expected_counts.3),
        ];
        for (reduce_pattern_symmetry, reduce_host_symmetry, expected) in modes {
            for strategy in STRATEGIES {
                for subiso in SUBISO_ALGS {
                    for automorphism_algorithm in
                        [AutomorphismAlgorithm::TrailIr, AutomorphismAlgorithm::Nauty]
                    {
                        let config = SubstructureMatchConfig {
                            match_algorithm: strategy,
                            subgraph_isomorphism_algorithm: subiso,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                            automorphism_algorithm,
                            reduce_pattern_symmetry,
                            reduce_host_symmetry,
                        };
                        let all = pattern.substructure_matches(&host, config);
                        let count = pattern.count_substructure_matches(&host, config);
                        let has = pattern.has_substructure_match(&host, config);
                        let first = pattern.first_substructure_match(&host, config);

                        assert_eq!(
                            count, expected,
                            "{strategy:?}/{subiso:?}/{automorphism_algorithm:?}"
                        );
                        assert_eq!(all.len() as u64, expected);
                        assert_eq!(has, expected > 0);
                        assert_eq!(first.is_some(), expected > 0);
                        if let Some(first_match) = first {
                            assert!(all.contains(&first_match));
                        }
                    }
                }
            }
        }
    }

    #[rstest]
    #[case::visit_embeddings(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        0
    )]
    #[case::has_match(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        1
    )]
    #[case::first_match(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        2
    )]
    #[case::count_matches(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        3
    )]
    #[should_panic(expected = "substructure matching does not evaluate molecule constraints")]
    fn test_molecule_substructure_match_modes_error(
        #[case] host: Molecule,
        #[case] pattern: Molecule,
        #[case] mode: u8,
    ) {
        let config = SubstructureMatchConfig {
            match_algorithm: GraphAndOverlays,
            subgraph_isomorphism_algorithm: BitsetRi,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        };
        match mode {
            0 => {
                let _: ControlFlow<()> =
                    pattern.visit_substructure_embeddings(&host, config, |_| {
                        ControlFlow::Continue(())
                    });
            }
            1 => {
                let _ = pattern.has_substructure_match(&host, config);
            }
            2 => {
                let _ = pattern.first_substructure_match(&host, config);
            }
            3 => {
                let _ = pattern.count_substructure_matches(&host, config);
            }
            _ => unreachable!(),
        }
    }
}
