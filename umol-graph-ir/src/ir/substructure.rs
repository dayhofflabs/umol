//! Substructure matching: occurrences of a pattern `Molecule` within a host
//! `Molecule`, each an injective pattern→host [`MoleculeCorrespondence`]. The receiver is the
//! pattern, so it parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the localized atom-bond skeleton then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use std::borrow::Cow;
use std::ops::ControlFlow;

use umol_graph_core::{
    Correspondence, ParticipantPosition, RelationData, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};

use super::atom::AtomForm;
use super::bond::BondForm;
use super::constraint::{AtomConstraintForm, BondConstraintForm, RingScope};
use super::correspondence::{
    induced_aromatic_systems, induced_bonds, induced_dative_bonds, induced_multicenter_bonds,
    induced_noncovalent_bonds, map_atom, map_ligands, MoleculeCorrespondence,
};
use super::entity::Entity;
use super::id::{AtomId, BondId};
use super::incidence::IncidenceNodeSelection;
use super::molecule::Molecule;
use super::ring::{RingConfig, RingModel, RingSetKind};
use super::stereo::coset_matches;
use super::traits::Lattice;

/// Strategy for `substructure_matches`, each composing over a
/// [`SubgraphIsomorphismAlgorithm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubstructureMatchAlgorithm {
    /// Match the localized atom-bond skeleton, then post-verify the N-ary / special
    /// overlays against the atom correspondence.
    GraphAndOverlays,
    /// Match the incidence (Levi) graph — true hypergraph matching for connectivity
    /// carried solely by overlays.
    Incidence,
}

/// Algorithms used to enumerate substructure matches.
///
/// This type deliberately has no default at the AST layer: every graph
/// algorithm selection remains explicit at the call site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureMatchConfig {
    /// Strategy used to represent and match molecule entities.
    pub match_algorithm: SubstructureMatchAlgorithm,
    /// Algorithm used to enumerate embeddings of the selected graph representation.
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    /// Algorithm used to derive relevant-ring constraints requested by the pattern.
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}

impl Molecule {
    /// Visits each occurrence of `self` (the pattern) within `host` as an injective
    /// pattern→host [`MoleculeCorrespondence`] until traversal completes or the
    /// visitor returns [`ControlFlow::Break`]. Pattern predicates are evaluated as
    /// `pattern.matches(host)` against the host atom/bond augmented with its derived
    /// topological constraints. Traversal is deterministic for a fixed
    /// representation, but its order is not a canonical ordering contract.
    pub fn visit_substructure_matches<B, F>(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
        mut visitor: F,
    ) -> ControlFlow<B>
    where
        F: FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    {
        match config.match_algorithm {
            SubstructureMatchAlgorithm::GraphAndOverlays => self
                .visit_substructure_matches_graph_and_overlays(
                    host,
                    config.subgraph_isomorphism_algorithm,
                    config.relevant_cycle_algorithm,
                    &mut visitor,
                ),
            SubstructureMatchAlgorithm::Incidence => self.visit_substructure_matches_incidence(
                host,
                config.subgraph_isomorphism_algorithm,
                config.relevant_cycle_algorithm,
                &mut visitor,
            ),
        }
    }

    /// Occurrences of `self` (the pattern) within `host`, collected from
    /// [`Molecule::visit_substructure_matches`].
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

    /// Host atom/bond match-targets: each stored entity with the topological constraints requested
    /// anywhere in the pattern folded in (last-wins). Ring constraints use the fixed Relevant
    /// projection through size 22. An unconstrained pattern never consults derived constraints, so
    /// element/bond patterns over SMILES-raised hosts skip the work.
    fn host_match_targets<'h>(
        &self,
        host: &'h Molecule,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> (Vec<Cow<'h, AtomForm>>, Vec<Cow<'h, BondForm>>) {
        let derive_atoms = self
            .atoms()
            .iter()
            .any(|a| !a.attributes.constraints.is_empty());
        let mut atom_ring_scopes = Vec::new();
        let mut derive_ring_degree = false;
        let mut derive_ring_valence = false;
        for constraint in self
            .atoms()
            .iter()
            .flat_map(|atom| atom.constraints().iter())
        {
            match constraint {
                AtomConstraintForm::RingDegree(_) => derive_ring_degree = true,
                AtomConstraintForm::RingValence(_) => derive_ring_valence = true,
                AtomConstraintForm::RingMembership(membership) => {
                    atom_ring_scopes.push(membership.scope);
                }
                _ => {}
            }
        }
        atom_ring_scopes.sort_unstable();
        atom_ring_scopes.dedup();

        let mut bond_ring_scopes: Vec<RingScope> = self
            .bonds()
            .iter()
            .flat_map(|bond| bond.constraints().iter())
            .filter_map(|constraint| match constraint {
                BondConstraintForm::RingMembership(membership) => Some(membership.scope),
                _ => None,
            })
            .collect();
        bond_ring_scopes.sort_unstable();
        bond_ring_scopes.dedup();

        let derive_rings = derive_ring_degree
            || derive_ring_valence
            || !atom_ring_scopes.is_empty()
            || !bond_ring_scopes.is_empty();
        let rings = derive_rings.then(|| {
            host.rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: 22,
                },
                RingConfig {
                    relevant_cycle_algorithm,
                    ..RingConfig::default()
                },
            )
        });

        let host_atoms = host
            .atoms()
            .iter()
            .map(|a| {
                if derive_atoms {
                    let mut constraints = a.derive_constraints(true);
                    if let Some(rings) = rings.as_ref() {
                        let ring = rings.atom(a.id);
                        if derive_ring_degree {
                            constraints.set(AtomConstraintForm::ring_degree(ring.ring_degree()));
                        }
                        if derive_ring_valence {
                            constraints.set(AtomConstraintForm::ring_valence(ring.ring_valence()));
                        }
                        for &scope in &atom_ring_scopes {
                            constraints.set(AtomConstraintForm::ring_membership(
                                scope,
                                ring.ring_membership(scope),
                            ));
                        }
                    }
                    Cow::Owned(a.attributes.clone().with_constraints(constraints))
                } else {
                    Cow::Borrowed(a.attributes)
                }
            })
            .collect();
        let derive_bonds = self
            .bonds()
            .iter()
            .any(|b| !b.attributes.constraints.is_empty());
        let host_bonds = host
            .bonds()
            .iter()
            .map(|b| {
                if derive_bonds {
                    let mut constraints = b.derive_constraints(true);
                    if let Some(rings) = rings.as_ref() {
                        let ring = rings.bond(b.id);
                        for &scope in &bond_ring_scopes {
                            constraints.set(BondConstraintForm::ring_membership(
                                scope,
                                ring.ring_membership(scope),
                            ));
                        }
                    }
                    Cow::Owned(b.attributes.clone().with_constraints(constraints))
                } else {
                    Cow::Borrowed(b.attributes)
                }
            })
            .collect();
        (host_atoms, host_bonds)
    }

    fn visit_substructure_matches_graph_and_overlays<B>(
        &self,
        host: &Molecule,
        subiso: SubgraphIsomorphismAlgorithm,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
        visitor: &mut impl FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return ControlFlow::Continue(());
        }
        let (host_atoms, host_bonds) = pattern.host_match_targets(host, relevant_cycle_algorithm);

        host.raw_graph().visit_subgraph_isomorphisms(
            pattern.raw_graph(),
            &mut |query_node, host_node| {
                pattern
                    .atom(AtomId::from(query_node))
                    .attributes
                    .matches(&host_atoms[host_node.index()])
            },
            &mut |query_edge, host_edge| {
                pattern
                    .bond(BondId::from(query_edge))
                    .attributes
                    .matches(&host_bonds[host_edge.index()])
            },
            subiso,
            |embedding| {
                let atoms = Correspondence::new(
                    embedding
                        .iter()
                        .enumerate()
                        .map(|(pattern_atom, &host_node)| {
                            (AtomId(pattern_atom as u32), AtomId::from(host_node))
                        })
                        .collect(),
                    embedding.len(),
                    host.atoms().count(),
                )
                .expect("subgraph match preserves atom correspondence invariants");
                match pattern.verify_overlays(host, atoms) {
                    Some(correspondence) => visitor(correspondence),
                    None => ControlFlow::Continue(()),
                }
            },
        )
    }

    /// Match on the incidence (Levi) graph: relations become pseudonodes wired to
    /// their participant atoms, so overlay-only connectivity (a 3c-2e bond, an H-bond
    /// that is the sole link) constrains placement — the case `GraphAndOverlays`
    /// degrades on. The Levi subiso supplies only the atom correspondence; the same
    /// exact `verify_overlays` then filters and builds the embedding, so this returns
    /// the identical match set as `GraphAndOverlays`.
    fn visit_substructure_matches_incidence<B>(
        &self,
        host: &Molecule,
        subiso: SubgraphIsomorphismAlgorithm,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
        visitor: &mut impl FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    ) -> ControlFlow<B> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return ControlFlow::Continue(());
        }
        let selection = IncidenceNodeSelection::constitution();
        let pattern_levi = pattern.incidence_graph(selection);
        let host_levi = host.incidence_graph(selection);
        let (host_atoms, host_bonds) = pattern.host_match_targets(host, relevant_cycle_algorithm);
        let atom_count = pattern.atoms().count();

        host_levi.graph().visit_subgraph_isomorphisms(
            pattern_levi.graph(),
            // Atoms/bonds carry their predicates; overlay pseudonodes match by
            // kind only (the exact AST/participation check is `verify_overlays`).
            &mut |pq, hq| match (pattern_levi.entity(pq), host_levi.entity(hq)) {
                (Entity::Atom(pa), Entity::Atom(ha)) => {
                    pattern.atom(pa).attributes.matches(&host_atoms[ha.index()])
                }
                (Entity::Bond(pb), Entity::Bond(hb)) => {
                    pattern.bond(pb).attributes.matches(&host_bonds[hb.index()])
                }
                (pe, he) => pe.kind() == he.kind(),
            },
            &mut |_, _| true,
            subiso,
            |embedding| {
                let atoms = Correspondence::new(
                    embedding[..atom_count]
                        .iter()
                        .enumerate()
                        .map(|(a, &host_node)| match host_levi.entity(host_node) {
                            Entity::Atom(id) => (AtomId(a as u32), id),
                            _ => unreachable!("a pattern atom node maps to a host atom node"),
                        })
                        .collect(),
                    pattern.atoms().count(),
                    host.atoms().count(),
                )
                .expect("correspondence producer preserves partial-bijection invariants");
                match pattern.verify_overlays(host, atoms) {
                    Some(correspondence) => visitor(correspondence),
                    None => ControlFlow::Continue(()),
                }
            },
        )
    }

    /// Post-verify a skeleton occurrence's overlays against the atom correspondence, returning the
    /// injective pattern→host [`MoleculeCorrespondence`] or `None` if any pattern overlay has no
    /// matching host overlay. Each N-ary / special overlay is located by **exact participant set**
    /// via the per-family inducer (which already checks dative donor/acceptor roles); the pattern
    /// overlay's predicate is then required to match the located host overlay's, and every pattern
    /// overlay must be matched. Stereo overlays are matched by the bespoke coset filter.
    fn verify_overlays(
        &self,
        host: &Molecule,
        atoms: Correspondence<AtomId>,
    ) -> Option<MoleculeCorrespondence> {
        let pattern = self;
        let bonds = induced_bonds(pattern, host, &atoms)?;

        let dative_bonds = induced_dative_bonds(pattern, host, &atoms)?;
        if dative_bonds.matched_pair_count() != pattern.dative_bonds().count() {
            return None;
        }
        for &(p, h) in dative_bonds.matched_pairs() {
            if !pattern
                .dative_bond(p)
                .attributes
                .matches(host.dative_bond(h).attributes)
            {
                return None;
            }
        }

        let aromatic_systems = induced_aromatic_systems(pattern, host, &atoms)?;
        if aromatic_systems.matched_pair_count() != pattern.aromatic_systems().count() {
            return None;
        }
        for &(p, h) in aromatic_systems.matched_pairs() {
            let p_view = pattern.aromatic_system(p);
            let h_view = host.aromatic_system(h);
            let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
            let host_atoms: Vec<AtomId> = h_view.atom_ids().collect();
            if !overlay_matches(
                p_view.attributes,
                h_view.attributes,
                &pat_atoms,
                &host_atoms,
                &atoms,
            ) {
                return None;
            }
        }

        let multicenter_bonds = induced_multicenter_bonds(pattern, host, &atoms)?;
        if multicenter_bonds.matched_pair_count() != pattern.multicenter_bonds().count() {
            return None;
        }
        for &(p, h) in multicenter_bonds.matched_pairs() {
            let p_view = pattern.multicenter_bond(p);
            let h_view = host.multicenter_bond(h);
            let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
            let host_atoms: Vec<AtomId> = h_view.atom_ids().collect();
            if !overlay_matches(
                p_view.attributes,
                h_view.attributes,
                &pat_atoms,
                &host_atoms,
                &atoms,
            ) {
                return None;
            }
        }

        let noncovalent_bonds = induced_noncovalent_bonds(pattern, host, &atoms)?;
        if noncovalent_bonds.matched_pair_count() != pattern.noncovalent_bonds().count() {
            return None;
        }
        for &(p, h) in noncovalent_bonds.matched_pairs() {
            if !pattern
                .noncovalent_bond(p)
                .attributes
                .matches(host.noncovalent_bond(h).attributes)
            {
                return None;
            }
        }

        // Stereo: a pattern stereo overlay matches iff the corresponding host site
        // bears a stereo element of the same class whose coset, reindexed from the
        // host ligand frame into the pattern's frame (via the atom correspondence),
        // is admitted by the pattern coset. An `Undetermined` pattern coset admits
        // both handednesses. TODO: a pattern that asserts stereo via `#T`/`#C` atom
        // /bond constraints rather than a `:stereo-atoms`/`:stereo-bonds` overlay is
        // not handled here — that needs the pattern run through stereo perception
        // (but not grounding, so no valence resolution).
        let mut stereo_atom = Vec::new();
        for sp in pattern.stereo_atoms().iter() {
            let host_site =
                map_atom(&atoms, sp.site_id()).expect("a matched pattern atom is matched");
            // `incident` returns stereo atoms where `host_site` is the site *or* a ligand; select the
            // one it is the site of (≤1 by the site-uniqueness invariant), not merely the first.
            let sh = host
                .stereo_atoms()
                .incident(host_site)
                .find(|sh| sh.site_id() == host_site)?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame = map_ligands(&atoms, sp.ligand_frame())
                .expect("matched pattern ligands are matched");
            let host_coset = sh.coset_for(frame)?;
            if !coset_matches(sp.coset(), &host_coset, sp.kind()) {
                return None;
            }
            stereo_atom.push((sp.id, sh.id));
        }
        let stereo_atoms = Correspondence::new(
            stereo_atom,
            pattern.stereo_atoms().count(),
            host.stereo_atoms().count(),
        )
        .ok()?;

        let mut stereo_bond = Vec::new();
        for sp in pattern.stereo_bonds().iter() {
            let host_site = bonds
                .right_of(sp.site_id())
                .expect("a matched pattern bond is matched");
            let sh = host.bond(host_site).stereo_bond()?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame = map_ligands(&atoms, sp.ligand_frame())
                .expect("matched pattern ligands are matched");
            let host_coset = sh.coset_for(frame)?;
            if !coset_matches(sp.coset(), &host_coset, sp.kind()) {
                return None;
            }
            stereo_bond.push((sp.id, sh.id));
        }
        let stereo_bonds = Correspondence::new(
            stereo_bond,
            pattern.stereo_bonds().count(),
            host.stereo_bonds().count(),
        )
        .ok()?;

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

/// `pattern_ast` matches `host_ast` for an overlay whose payload is position-indexed by member
/// (aromatic / multicenter electron counts). The two overlays store their members in their own
/// participant order and `matches` compares the count vector whole, so the pattern payload is first
/// reindexed into the host's member order (via the atom correspondence) with
/// [`RelationData::on_permutation`].
fn overlay_matches<D: Lattice + RelationData>(
    pattern_ast: &D,
    host_ast: &D,
    pattern_atoms: &[AtomId],
    host_atoms: &[AtomId],
    atoms: &Correspondence<AtomId>,
) -> bool {
    if pattern_ast.is_permutation_invariant() {
        return pattern_ast.matches(host_ast);
    }
    let order: Vec<ParticipantPosition> = host_atoms
        .iter()
        .map(|&host_atom| {
            let pattern_atom = atoms
                .left_of(host_atom)
                .expect("host overlay atom is matched");
            ParticipantPosition(
                pattern_atoms
                    .iter()
                    .position(|&a| a == pattern_atom)
                    .expect("host member maps to a pattern member") as u32,
            )
        })
        .collect();
    let mut probe = pattern_ast.clone();
    probe.on_permutation(&order);
    probe.matches(host_ast)
}

#[cfg(test)]
mod tests {
    use std::ops::ControlFlow;

    use rstest::rstest;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{
        ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
    };
    use umol_graph_core::{
        RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
        ARCMATCH_DEFAULT_PATH_LENGTH,
    };

    use super::super::id::AtomId;
    use super::super::molecule::Molecule;
    use super::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
    use super::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
    use crate::mol_dsl;

    const SUBISO_ALGS: [SubgraphIsomorphismAlgorithm; 6] = [
        Vf2,
        Ullmann,
        Ri,
        ArcMatch {
            path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
        },
        Vf2Rdkit,
        RayKirsch,
    ];

    const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

    #[rstest]
    #[case::skeleton(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::noncovalent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]] :noncovalent-bonds [{:atoms [0 2] :type "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#),
        vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]]
    )]
    #[case::no_match(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"]]}"#),
        vec![]
    )]
    fn test_molecule_ast_visit_substructure_matches(
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
    #[case::skeleton(
        mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    fn test_molecule_ast_visit_substructure_matches_termination(
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
    #[case::skeleton(
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
    #[case::noncovalent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]] :noncovalent-bonds [{:atoms [0 2] :type "Hbd"}]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#),
        vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]]
    )]
    #[case::noncovalent_absent(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#),
        vec![]
    )]
    #[case::dative(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#),
        vec![vec![AtomId(0), AtomId(1)]]
    )]
    #[case::dative_roles_swapped(
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [1] :acceptor 0 :type "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#),
        vec![]
    )]
    #[case::stereo_chiral(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_enantiomer(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        vec![]
    )]
    #[case::stereo_agnostic_in_r(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_agnostic_in_s(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_reframed(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        vec![vec![AtomId(0), AtomId(3), AtomId(2), AtomId(1)]]
    )]
    #[case::stereo_reframed_enantiomer(
        mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol_dsl!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        vec![]
    )]
    #[case::stereo_bond(
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :type "Ct1"}]}"#),
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :type "Ct1"}]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)]]
    )]
    fn test_molecule_ast_substructure_matches(
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
}
