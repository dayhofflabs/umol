//! Substructure matching: occurrences of a pattern `MoleculeAst` within a host
//! `MoleculeAst`, as [`MoleculeEmbedding`]s. The receiver is the pattern, so it
//! parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the localized atom-bond skeleton then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use std::collections::HashSet;

use umol_graph_core::SubgraphIsomorphismAlgorithm;

use super::atom::AtomAst;
use super::bond::BondAst;
use super::embedding::MoleculeEmbedding;
use super::id::AtomId;
use super::molecule::MoleculeAst;
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
    /// Occurrences of `self` (the pattern) within `host`, one [`MoleculeEmbedding`]
    /// per occurrence. Pattern predicates are evaluated as `pattern.matches(host)`
    /// against the host atom/bond augmented with its derived topological
    /// constraints.
    pub fn substructure_matches<'h>(
        &self,
        host: &'h MoleculeAst,
        strategy: SubstructureMatchAlgorithm,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
        match strategy {
            SubstructureMatchAlgorithm::GraphAndOverlays => {
                self.substructure_matches_graph_and_overlays(host, subiso)
            }
            SubstructureMatchAlgorithm::Incidence => {
                unimplemented!("Incidence substructure strategy")
            }
        }
    }

    fn substructure_matches_graph_and_overlays<'h>(
        &self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return Vec::new();
        }
        // Host match-targets: stored constraints with derived topological ones folded
        // in (last-wins). Precomputed once, not per candidate.
        let host_atoms: Vec<AtomAst> = host
            .atoms()
            .iter()
            .map(|a| a.ast.clone().with_constraints(a.derive_constraints()))
            .collect();
        let host_bonds: Vec<BondAst> = host
            .bonds()
            .iter()
            .map(|b| b.ast.clone().with_constraints(b.derive_constraints()))
            .collect();

        host.graph()
            .subgraph_isomorphisms(
                &pattern.graph(),
                &mut |query_atom, host_atom| {
                    pattern.atom(query_atom).ast.matches(&host_atoms[host_atom.index()])
                },
                &mut |query_bond, host_bond| {
                    pattern.bond(query_bond).ast.matches(&host_bonds[host_bond.index()])
                },
                subiso,
            )
            .into_iter()
            .filter_map(|atom_map| pattern.verify_overlays(host, atom_map))
            .collect()
    }

    /// Post-verify a skeleton occurrence's overlays against the atom correspondence,
    /// returning the enriched embedding or `None` if any pattern overlay has no
    /// matching host overlay. Each N-ary / special overlay is located by **exact
    /// participant set** (`connecting`); dative donor/acceptor roles are checked
    /// explicitly. Stereo overlays are deferred to the coset post-filter.
    fn verify_overlays<'h>(
        &self,
        host: &'h MoleculeAst,
        atom_map: Vec<AtomId>,
    ) -> Option<MoleculeEmbedding<'h>> {
        let pattern = self;

        let mut dative = Vec::new();
        for d in pattern.dative_bonds().iter() {
            let host_dative = host
                .dative_bonds()
                .connecting(d.atom_ids().map(|a| atom_map[a.index()]))?;
            if atom_map[d.acceptor_id().index()] != host_dative.acceptor_id() {
                return None;
            }
            let pattern_donors: HashSet<AtomId> =
                d.donor_ids().map(|a| atom_map[a.index()]).collect();
            let host_donors: HashSet<AtomId> = host_dative.donor_ids().collect();
            if pattern_donors != host_donors || !d.ast.matches(host_dative.ast) {
                return None;
            }
            dative.push(host_dative.id);
        }

        let mut noncovalent = Vec::new();
        for nc in pattern.noncovalent_bonds().iter() {
            let [a, b] = nc.atom_ids();
            let host_nc = host
                .noncovalent_bonds()
                .connecting(atom_map[a.index()], atom_map[b.index()])?;
            if !nc.ast.matches(host_nc.ast) {
                return None;
            }
            noncovalent.push(host_nc.id);
        }

        let mut aromatic = Vec::new();
        for ar in pattern.aromatic_systems().iter() {
            let host_ar = host
                .aromatic_systems()
                .connecting(ar.atom_ids().map(|a| atom_map[a.index()]))?;
            if !ar.ast.matches(host_ar.ast) {
                return None;
            }
            aromatic.push(host_ar.id);
        }

        let mut multicenter = Vec::new();
        for mc in pattern.multicenter_bonds().iter() {
            let host_mc = host
                .multicenter_bonds()
                .connecting(mc.atom_ids().map(|a| atom_map[a.index()]))?;
            if !mc.ast.matches(host_mc.ast) {
                return None;
            }
            multicenter.push(host_mc.id);
        }

        Some(MoleculeEmbedding::from_match(
            host,
            pattern,
            atom_map,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            Vec::new(),
            Vec::new(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{
        ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
    };
    use umol_graph_core::{SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH};

    use crate::mol;

    use super::super::id::AtomId;
    use super::SubstructureMatchAlgorithm::GraphAndOverlays;

    const SUBISO_ALGS: [SubgraphIsomorphismAlgorithm; 6] = [
        Vf2,
        Ullmann,
        Ri,
        ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH },
        Vf2Rdkit,
        RayKirsch,
    ];

    // Pattern C-C in host C-C-O: matches the C-C edge in both orientations, never the
    // C-O edge. The same occurrence set under every subgraph-isomorphism algorithm.
    #[rstest]
    fn test_molecule_ast_substructure_matches() {
        let host = mol!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let pattern = mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        for subiso in SUBISO_ALGS {
            let mut maps: Vec<Vec<AtomId>> = pattern
                .substructure_matches(&host, GraphAndOverlays, subiso)
                .iter()
                .map(|e| e.host_atoms().to_vec())
                .collect();
            maps.sort();
            assert_eq!(
                maps,
                vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]],
                "subiso {subiso:?}"
            );
        }
    }

    // Overlay verification: a pattern carrying a hydrogen-bond noncovalent overlay
    // (its two atoms otherwise unbonded) matches only host atom pairs that actually
    // carry that overlay, and the embedding records the host overlay id. Absent the
    // host overlay there is no match.
    #[rstest]
    fn test_molecule_ast_substructure_matches_overlay() {
        // host: C-C-C chain with a hydrogen bond between atoms 0 and 2.
        let host = mol!(
            r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]
                :noncovalent-bonds [{:a 0 :b 2 :type "Hbd"}]}"#
        );
        // same skeleton, without the noncovalent overlay.
        let host_bare = mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        // pattern: two carbons joined only by a hydrogen bond.
        let pattern = mol!(
            r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:a 0 :b 1 :type "Hbd"}]}"#
        );
        for subiso in SUBISO_ALGS {
            let embeddings = pattern.substructure_matches(&host, GraphAndOverlays, subiso);
            let mut maps: Vec<Vec<AtomId>> =
                embeddings.iter().map(|e| e.host_atoms().to_vec()).collect();
            maps.sort();
            assert_eq!(
                maps,
                vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]],
                "subiso {subiso:?}"
            );
            assert!(
                embeddings
                    .iter()
                    .all(|e| e.host_noncovalent_bonds().len() == 1),
                "subiso {subiso:?}"
            );
            assert!(
                pattern
                    .substructure_matches(&host_bare, GraphAndOverlays, subiso)
                    .is_empty(),
                "subiso {subiso:?}"
            );
        }
    }
}
