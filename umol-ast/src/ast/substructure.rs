//! Substructure matching: occurrences of a pattern `MoleculeAst` within a host
//! `MoleculeAst`, each an injective pattern→host [`MoleculeCorrespondence`]. The receiver is the
//! pattern, so it parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the localized atom-bond skeleton then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use std::borrow::Cow;

use umol_graph_core::{
    Correspondence, NodeId, ParticipantPosition, RelationData, SubgraphIsomorphismAlgorithm,
};

use super::atom::AtomAst;
use super::bond::BondAst;
use super::correspondence::{
    induced_aromatic_systems, induced_bonds, induced_dative_bonds, induced_multicenter_bonds,
    induced_noncovalent_bonds, map_atom, map_ligands, MoleculeCorrespondence,
};
use super::entity::Entity;
use super::id::{AtomId, BondId};
use super::incidence::IncidenceNodeSelection;
use super::molecule::MoleculeAst;
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

impl MoleculeAst {
    /// Occurrences of `self` (the pattern) within `host`, one injective pattern→host
    /// [`MoleculeCorrespondence`] per occurrence. Pattern predicates are evaluated as
    /// `pattern.matches(host)` against the host atom/bond augmented with its derived topological
    /// constraints.
    pub fn substructure_matches(
        &self,
        host: &MoleculeAst,
        strategy: SubstructureMatchAlgorithm,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeCorrespondence> {
        match strategy {
            SubstructureMatchAlgorithm::GraphAndOverlays => {
                self.substructure_matches_graph_and_overlays(host, subiso)
            }
            SubstructureMatchAlgorithm::Incidence => {
                self.substructure_matches_incidence(host, subiso)
            }
        }
    }

    /// Host atom/bond match-targets: each stored entity with its derived
    /// topological constraints folded in (last-wins) — but only for the entity kind
    /// the *pattern* constrains. An unconstrained pattern never consults the host's
    /// derived constraints (empty pattern constraints match any target), so deriving
    /// them is wasted; element/bond patterns over SMILES-raised hosts skip it.
    fn host_match_targets<'h>(
        &self,
        host: &'h MoleculeAst,
    ) -> (Vec<Cow<'h, AtomAst>>, Vec<Cow<'h, BondAst>>) {
        let derive_atoms = self.atoms().iter().any(|a| !a.ast.constraints.is_empty());
        let host_atoms = host
            .atoms()
            .iter()
            .map(|a| {
                if derive_atoms {
                    Cow::Owned(a.ast.clone().with_constraints(a.derive_constraints(true)))
                } else {
                    Cow::Borrowed(a.ast)
                }
            })
            .collect();
        let derive_bonds = self.bonds().iter().any(|b| !b.ast.constraints.is_empty());
        let host_bonds = host
            .bonds()
            .iter()
            .map(|b| {
                if derive_bonds {
                    Cow::Owned(b.ast.clone().with_constraints(b.derive_constraints(true)))
                } else {
                    Cow::Borrowed(b.ast)
                }
            })
            .collect();
        (host_atoms, host_bonds)
    }

    fn substructure_matches_graph_and_overlays(
        &self,
        host: &MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeCorrespondence> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return Vec::new();
        }
        let (host_atoms, host_bonds) = pattern.host_match_targets(host);

        host.raw_graph()
            .subgraph_isomorphisms(
                pattern.raw_graph(),
                &mut |query_node, host_node| {
                    pattern
                        .atom(AtomId::from(query_node))
                        .ast
                        .matches(&host_atoms[host_node.index()])
                },
                &mut |query_edge, host_edge| {
                    pattern
                        .bond(BondId::from(query_edge))
                        .ast
                        .matches(&host_bonds[host_edge.index()])
                },
                subiso,
            )
            .into_iter()
            .filter_map(|atoms| pattern.verify_overlays(host, atoms))
            .collect()
    }

    /// Match on the incidence (Levi) graph: relations become pseudonodes wired to
    /// their participant atoms, so overlay-only connectivity (a 3c-2e bond, an H-bond
    /// that is the sole link) constrains placement — the case `GraphAndOverlays`
    /// degrades on. The Levi subiso supplies only the atom correspondence; the same
    /// exact `verify_overlays` then filters and builds the embedding, so this returns
    /// the identical match set as `GraphAndOverlays`.
    fn substructure_matches_incidence(
        &self,
        host: &MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeCorrespondence> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return Vec::new();
        }
        let selection = IncidenceNodeSelection::constitution();
        let pattern_levi = pattern.incidence_graph(selection);
        let host_levi = host.incidence_graph(selection);
        let (host_atoms, host_bonds) = pattern.host_match_targets(host);
        let atom_count = pattern.atoms().count();

        host_levi
            .graph()
            .subgraph_isomorphisms(
                pattern_levi.graph(),
                // Atoms/bonds carry their predicates; overlay pseudonodes match by
                // kind only (the exact AST/participation check is `verify_overlays`).
                &mut |pq, hq| match (pattern_levi.entity(pq), host_levi.entity(hq)) {
                    (Entity::Atom(pa), Entity::Atom(ha)) => {
                        pattern.atom(pa).ast.matches(&host_atoms[ha.index()])
                    }
                    (Entity::Bond(pb), Entity::Bond(hb)) => {
                        pattern.bond(pb).ast.matches(&host_bonds[hb.index()])
                    }
                    (pe, he) => pe.kind() == he.kind(),
                },
                &mut |_, _| true,
                subiso,
            )
            .into_iter()
            .filter_map(|levi_match| {
                let atoms = Correspondence::new(
                    (0..atom_count as u32)
                        .map(|a| {
                            let host_node = levi_match
                                .right_of(NodeId(a))
                                .expect("a pattern atom node is mated");
                            match host_levi.entity(host_node) {
                                Entity::Atom(id) => (NodeId(a), NodeId::from(id)),
                                _ => unreachable!("a pattern atom node maps to a host atom node"),
                            }
                        })
                        .collect(),
                    pattern.atoms().count(),
                    host.atoms().count(),
                );
                pattern.verify_overlays(host, atoms)
            })
            .collect()
    }

    /// Post-verify a skeleton occurrence's overlays against the atom correspondence, returning the
    /// injective pattern→host [`MoleculeCorrespondence`] or `None` if any pattern overlay has no
    /// matching host overlay. Each N-ary / special overlay is located by **exact participant set**
    /// via the per-family inducer (which already checks dative donor/acceptor roles); the pattern
    /// overlay's predicate is then required to match the located host overlay's, and every pattern
    /// overlay must be mated. Stereo overlays are matched by the bespoke coset filter.
    fn verify_overlays(
        &self,
        host: &MoleculeAst,
        atoms: Correspondence<NodeId>,
    ) -> Option<MoleculeCorrespondence> {
        let pattern = self;
        let bonds = induced_bonds(pattern, host, &atoms);

        let dative_bonds = induced_dative_bonds(pattern, host, &atoms);
        if dative_bonds.mate_count() != pattern.dative_bonds().count() {
            return None;
        }
        for &(p, h) in dative_bonds.mates() {
            if !pattern.dative_bond(p).ast.matches(host.dative_bond(h).ast) {
                return None;
            }
        }

        let aromatic_systems = induced_aromatic_systems(pattern, host, &atoms);
        if aromatic_systems.mate_count() != pattern.aromatic_systems().count() {
            return None;
        }
        for &(p, h) in aromatic_systems.mates() {
            let p_view = pattern.aromatic_system(p);
            let h_view = host.aromatic_system(h);
            let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
            let host_atoms: Vec<AtomId> = h_view.atom_ids().collect();
            if !overlay_matches(p_view.ast, h_view.ast, &pat_atoms, &host_atoms, &atoms) {
                return None;
            }
        }

        let multicenter_bonds = induced_multicenter_bonds(pattern, host, &atoms);
        if multicenter_bonds.mate_count() != pattern.multicenter_bonds().count() {
            return None;
        }
        for &(p, h) in multicenter_bonds.mates() {
            let p_view = pattern.multicenter_bond(p);
            let h_view = host.multicenter_bond(h);
            let pat_atoms: Vec<AtomId> = p_view.atom_ids().collect();
            let host_atoms: Vec<AtomId> = h_view.atom_ids().collect();
            if !overlay_matches(p_view.ast, h_view.ast, &pat_atoms, &host_atoms, &atoms) {
                return None;
            }
        }

        let noncovalent_bonds = induced_noncovalent_bonds(pattern, host, &atoms);
        if noncovalent_bonds.mate_count() != pattern.noncovalent_bonds().count() {
            return None;
        }
        for &(p, h) in noncovalent_bonds.mates() {
            if !pattern
                .noncovalent_bond(p)
                .ast
                .matches(host.noncovalent_bond(h).ast)
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
                map_atom(&atoms, sp.site_id()).expect("a matched pattern atom is mated");
            // `incident` returns stereo atoms where `host_site` is the site *or* a ligand; select the
            // one it is the site of (≤1 by the site-uniqueness invariant), not merely the first.
            let sh = host
                .stereo_atoms()
                .incident(host_site)
                .find(|sh| sh.site_id() == host_site)?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame =
                map_ligands(&atoms, sp.ligand_frame()).expect("matched pattern ligands are mated");
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
        );

        let mut stereo_bond = Vec::new();
        for sp in pattern.stereo_bonds().iter() {
            let host_site = bonds
                .right_of(sp.site_id())
                .expect("a matched pattern bond is mated");
            let sh = host.bond(host_site).stereo_bond()?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame =
                map_ligands(&atoms, sp.ligand_frame()).expect("matched pattern ligands are mated");
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
        );

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
    atoms: &Correspondence<NodeId>,
) -> bool {
    if pattern_ast.is_permutation_invariant() {
        return pattern_ast.matches(host_ast);
    }
    let order: Vec<ParticipantPosition> = host_atoms
        .iter()
        .map(|&host_atom| {
            let pattern_atom = AtomId::from(
                atoms
                    .left_of(NodeId::from(host_atom))
                    .expect("host overlay atom is matched"),
            );
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
    use rstest::rstest;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{
        ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
    };
    use umol_graph_core::{SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH};

    use super::super::id::AtomId;
    use super::super::molecule::MoleculeAst;
    use super::SubstructureMatchAlgorithm;
    use super::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
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
    fn test_molecule_ast_substructure_matches(
        #[case] host: MoleculeAst,
        #[case] pattern: MoleculeAst,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        for strategy in STRATEGIES {
            for subiso in SUBISO_ALGS {
                let mut occurrences: Vec<Vec<AtomId>> = pattern
                    .substructure_matches(&host, strategy, subiso)
                    .iter()
                    .map(|c| {
                        c.atoms()
                            .mates()
                            .iter()
                            .map(|&(_, host)| AtomId::from(host))
                            .collect()
                    })
                    .collect();
                occurrences.sort();
                assert_eq!(occurrences, expected, "{strategy:?}/{subiso:?}");
            }
        }
    }
}
