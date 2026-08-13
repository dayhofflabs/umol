//! Substructure matching: occurrences of a pattern `Molecule` within a host
//! `Molecule`, each an injective pattern→host [`MoleculeCorrespondence`]. The receiver is the
//! pattern, so it parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the atom-bond topology then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use std::ops::ControlFlow;

use thiserror::Error;
use umol_graph_core::{
    Correspondence, ParticipantPosition, RelationData, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};

use super::atom::AtomForm;
use super::bond::BondForm;
use super::constraint::{AtomConstraintForm, BondConstraintForm};
use super::correspondence::{
    induced_aromatic_systems, induced_bonds, induced_dative_bonds, induced_multicenter_bonds,
    induced_noncovalent_bonds, map_atom, map_ligands, MoleculeCorrespondence,
};
use super::entity::Entity;
use super::id::{AtomId, BondId};
use super::incidence::{Incidence, IncidenceLevel};
use super::molecule::Molecule;
use super::ring::{RingConfig, RingModel, RingSet, RingSetKind};
use super::stereo::coset_matches;
use super::traits::Lattice;

/// Strategy for `substructure_matches`, each composing over a
/// [`SubgraphIsomorphismAlgorithm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubstructureMatchAlgorithm {
    /// Match the atom-bond topology, then post-verify the N-ary / special
    /// overlays against the atom correspondence.
    GraphAndOverlays,
    /// Match the incidence (Levi) graph — true hypergraph matching for connectivity
    /// carried solely by overlays.
    Incidence,
}

/// Algorithms used to enumerate substructure matches.
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
}

/// Matching rejected the pattern: it uses a construct matching does not
/// evaluate. Sound but incomplete — a query carrying an unevaluated construct
/// fails loudly instead of returning silently weakened results.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum SubstructureMatchError {
    /// The pattern's molecule-scope `Constraints` list (relational and
    /// molecule leaves, combinators) is not evaluated by matching.
    #[error(
        "pattern carries {count} molecule-scope constraint(s), which matching does not evaluate"
    )]
    MoleculeScopeConstraints { count: usize },
}

impl Molecule {
    /// Visits each occurrence of `self` (the pattern) within `host` as an injective
    /// pattern→host [`MoleculeCorrespondence`] until traversal completes or the
    /// visitor returns [`ControlFlow::Break`]. Pattern entity predicates are
    /// evaluated per key against the host's constraint reading under the closure
    /// (`constraints().satisfies`); ring keys use the fixed Relevant projection
    /// through size 22. Traversal is deterministic for a fixed representation,
    /// but its order is not a canonical ordering contract.
    ///
    /// # Errors
    ///
    /// A pattern carrying molecule-scope constraints is rejected: matching does
    /// not evaluate them, and failing loudly beats silently weakened results.
    pub fn visit_substructure_matches<B, F>(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
        mut visitor: F,
    ) -> Result<ControlFlow<B>, SubstructureMatchError>
    where
        F: FnMut(MoleculeCorrespondence) -> ControlFlow<B>,
    {
        if !self.constraints().is_empty() {
            return Err(SubstructureMatchError::MoleculeScopeConstraints {
                count: self.constraints().len(),
            });
        }
        Ok(match config.match_algorithm {
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
        })
    }

    /// Occurrences of `self` (the pattern) within `host`, collected from
    /// [`Molecule::visit_substructure_matches`].
    ///
    /// # Errors
    ///
    /// A pattern carrying molecule-scope constraints is rejected: matching does
    /// not evaluate them, and failing loudly beats silently weakened results.
    pub fn substructure_matches(
        &self,
        host: &Molecule,
        config: SubstructureMatchConfig,
    ) -> Result<Vec<MoleculeCorrespondence>, SubstructureMatchError> {
        let mut occurrences = Vec::new();
        let _: ControlFlow<()> =
            self.visit_substructure_matches(host, config, |correspondence| {
                occurrences.push(correspondence);
                ControlFlow::Continue(())
            })?;
        Ok(occurrences)
    }

    /// The host ring set, built once per match run iff a pattern constraint
    /// carries a ring key, with the fixed Relevant projection through size 22.
    fn host_ring_context(
        &self,
        host: &Molecule,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Option<RingSet> {
        let needs_rings = self
            .atoms()
            .iter()
            .flat_map(|atom| atom.attributes.constraints.iter())
            .any(|constraint| {
                matches!(
                    constraint,
                    AtomConstraintForm::RingDegree(_)
                        | AtomConstraintForm::RingValence(_)
                        | AtomConstraintForm::RingMembership(_)
                )
            })
            || self
                .bonds()
                .iter()
                .flat_map(|bond| bond.attributes.constraints.iter())
                .any(|constraint| matches!(constraint, BondConstraintForm::RingMembership(_)));
        needs_rings.then(|| {
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
            .into_ring_set()
        })
    }

    /// Per pattern×host entity constraint-satisfaction tables, evaluated once
    /// per run through the hosts' constraint readings under the closure; `None`
    /// when no pattern entity of that family carries constraints (the common
    /// element/bond-pattern case, which then skips all derivation).
    fn constraint_tables(
        &self,
        host: &Molecule,
        rings: Option<&RingSet>,
    ) -> (Option<Vec<bool>>, Option<Vec<bool>>) {
        let host_atom_count = host.atoms().count();
        let atoms = self
            .atoms()
            .iter()
            .any(|a| !a.attributes.constraints.is_empty())
            .then(|| {
                let mut table = vec![true; self.atoms().count() * host_atom_count];
                for p in self.atoms().iter() {
                    if p.attributes.constraints.is_empty() {
                        continue;
                    }
                    for h in host.atoms().iter() {
                        let mut reading = h.constraints();
                        if let Some(rings) = rings {
                            reading = reading.with_rings(rings);
                        }
                        table[p.id.index() * host_atom_count + h.id.index()] =
                            reading.satisfies(&p.attributes.constraints);
                    }
                }
                table
            });
        let host_bond_count = host.bonds().count();
        let bonds = self
            .bonds()
            .iter()
            .any(|b| !b.attributes.constraints.is_empty())
            .then(|| {
                let mut table = vec![true; self.bonds().count() * host_bond_count];
                for p in self.bonds().iter() {
                    if p.attributes.constraints.is_empty() {
                        continue;
                    }
                    for h in host.bonds().iter() {
                        let mut reading = h.constraints();
                        if let Some(rings) = rings {
                            reading = reading.with_rings(rings);
                        }
                        table[p.id.index() * host_bond_count + h.id.index()] =
                            reading.satisfies(&p.attributes.constraints);
                    }
                }
                table
            });
        (atoms, bonds)
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
        let rings = pattern.host_ring_context(host, relevant_cycle_algorithm);
        let (atom_table, bond_table) = pattern.constraint_tables(host, rings.as_ref());
        let host_atom_count = host.atoms().count();
        let host_bond_count = host.bonds().count();

        host.raw_graph().visit_subgraph_isomorphisms(
            pattern.raw_graph(),
            &mut |query_node, host_node| {
                let pa = AtomId::from(query_node);
                let ha = AtomId::from(host_node);
                atom_fields_match(pattern.atom(pa).attributes, host.atom(ha).attributes)
                    && atom_table
                        .as_ref()
                        .is_none_or(|table| table[pa.index() * host_atom_count + ha.index()])
            },
            &mut |query_edge, host_edge| {
                let pb = BondId::from(query_edge);
                let hb = BondId::from(host_edge);
                bond_fields_match(pattern.bond(pb).attributes, host.bond(hb).attributes)
                    && bond_table
                        .as_ref()
                        .is_none_or(|table| table[pb.index() * host_bond_count + hb.index()])
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
        let pattern_levi = pattern.incidence_graph(IncidenceLevel::Constitution);
        let host_levi = host.incidence_graph(IncidenceLevel::Constitution);
        let rings = pattern.host_ring_context(host, relevant_cycle_algorithm);
        let (atom_table, bond_table) = pattern.constraint_tables(host, rings.as_ref());
        let host_atom_count = host.atoms().count();
        let host_bond_count = host.bonds().count();
        let atom_count = pattern.atoms().count();

        host_levi.graph().visit_subgraph_isomorphisms(
            pattern_levi.graph(),
            // Atoms/bonds carry their predicates; overlay pseudonodes match by
            // kind only (the exact form/participation check is `verify_overlays`).
            &mut |pq, hq| match (pattern_levi.entity(pq), host_levi.entity(hq)) {
                (Entity::Atom(pa), Entity::Atom(ha)) => {
                    atom_fields_match(pattern.atom(pa).attributes, host.atom(ha).attributes)
                        && atom_table
                            .as_ref()
                            .is_none_or(|table| table[pa.index() * host_atom_count + ha.index()])
                }
                (Entity::Bond(pb), Entity::Bond(hb)) => {
                    bond_fields_match(pattern.bond(pb).attributes, host.bond(hb).attributes)
                        && bond_table
                            .as_ref()
                            .is_none_or(|table| table[pb.index() * host_bond_count + hb.index()])
                }
                (pe, he) => pe.kind() == he.kind(),
            },
            &mut |pattern_edge, host_edge| match (
                pattern_levi.incidence(pattern_edge),
                host_levi.incidence(host_edge),
            ) {
                (Incidence::AromaticParticipant(pattern), Incidence::AromaticParticipant(host))
                | (
                    Incidence::MulticenterParticipant(pattern),
                    Incidence::MulticenterParticipant(host),
                ) => pattern.matches(host),
                (pattern, host) => pattern == host,
            },
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

    /// Post-verify a topology occurrence's overlays against the atom correspondence, returning the
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
fn overlay_matches<D: Lattice + RelationData>(
    pattern_form: &D,
    host_form: &D,
    pattern_atoms: &[AtomId],
    host_atoms: &[AtomId],
    atoms: &Correspondence<AtomId>,
) -> bool {
    if pattern_form.is_permutation_invariant() {
        return pattern_form.matches(host_form);
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
    let mut probe = pattern_form.clone();
    probe.on_permutation(&order);
    probe.matches(host_form)
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
    use super::{SubstructureMatchAlgorithm, SubstructureMatchConfig, SubstructureMatchError};
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
                };
                let mut occurrences: Vec<Vec<AtomId>> = Vec::new();
                let flow: ControlFlow<()> = pattern
                    .visit_substructure_matches(&host, config, |correspondence| {
                        occurrences.push(
                            correspondence
                                .atoms()
                                .matched_pairs()
                                .iter()
                                .map(|&(_, host)| host)
                                .collect(),
                        );
                        ControlFlow::Continue(())
                    })
                    .unwrap();
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
                };
                let first = pattern
                    .visit_substructure_matches(&host, config, |correspondence| {
                        ControlFlow::Break(
                            correspondence
                                .atoms()
                                .matched_pairs()
                                .iter()
                                .map(|&(_, host)| host)
                                .collect::<Vec<_>>(),
                        )
                    })
                    .unwrap();
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
    #[case::one_molecule_scope_constraint(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        SubstructureMatchError::MoleculeScopeConstraints { count: 1 }
    )]
    #[case::two_molecule_scope_constraints(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}} {:charge-sum {:sum 0}}]}"#),
        SubstructureMatchError::MoleculeScopeConstraints { count: 2 }
    )]
    fn test_molecule_visit_substructure_matches_error(
        #[case] pattern: Molecule,
        #[case] expected: SubstructureMatchError,
    ) {
        let host = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        for strategy in STRATEGIES {
            let config = SubstructureMatchConfig {
                match_algorithm: strategy,
                subgraph_isomorphism_algorithm: Vf2,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            };
            let result = pattern
                .visit_substructure_matches(&host, config, |_| ControlFlow::Continue::<()>(()));
            assert_eq!(result, Err(expected), "{strategy:?}");
        }
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
    #[case::atom_aromatic_kekule_flag(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        mol_dsl!(r#"{:atoms ["C#a+"] :bonds []}"#),
        vec![vec![AtomId(0)], vec![AtomId(1)]]
    )]
    #[case::atom_not_aromatic_kekule_flag(
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
    #[case::bond_aromatic_asserted_conflict(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*#a"]]}"#),
        vec![]
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
    #[case::stereo_bond(
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :attrs "Ct1"}]}"#),
        mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :attrs "Ct1"}]}"#),
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
                        },
                    )
                    .unwrap()
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
    #[case::one_molecule_scope_constraint(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#),
        SubstructureMatchError::MoleculeScopeConstraints { count: 1 }
    )]
    fn test_molecule_substructure_matches_error(
        #[case] pattern: Molecule,
        #[case] expected: SubstructureMatchError,
    ) {
        let host = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let config = SubstructureMatchConfig {
            match_algorithm: GraphAndOverlays,
            subgraph_isomorphism_algorithm: Vf2,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        };
        assert_eq!(pattern.substructure_matches(&host, config), Err(expected));
    }
}
