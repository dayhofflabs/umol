//! Substructure matching: occurrences of a pattern `MoleculeAst` within a host
//! `MoleculeAst`, as [`MoleculeEmbedding`]s. The receiver is the pattern, so it
//! parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the localized atom-bond skeleton then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use std::borrow::Cow;
use std::collections::HashSet;

use umol_graph_core::{NodeId, SubgraphIsomorphismAlgorithm};

use super::atom::AtomAst;
use super::bond::BondAst;
use super::embedding::MoleculeEmbedding;
use super::entity::Entity;
use super::id::{AtomId, BondId};
use super::incidence::IncidenceNodeSelection;
use super::ligand::StereoLigand;
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

    fn substructure_matches_graph_and_overlays<'h>(
        &self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return Vec::new();
        }
        let (host_atoms, host_bonds) = pattern.host_match_targets(host);

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

    /// Match on the incidence (Levi) graph: relations become pseudonodes wired to
    /// their participant atoms, so overlay-only connectivity (a 3c-2e bond, an H-bond
    /// that is the sole link) constrains placement — the case `GraphAndOverlays`
    /// degrades on. The Levi subiso supplies only the atom correspondence; the same
    /// exact `verify_overlays` then filters and builds the embedding, so this returns
    /// the identical match set as `GraphAndOverlays`.
    fn substructure_matches_incidence<'h>(
        &self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
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
            .filter_map(|m| {
                let atom_map: Vec<AtomId> = (0..atom_count)
                    .map(|a| match host_levi.entity(NodeId(m[a] as u32)) {
                        Entity::Atom(id) => id,
                        _ => unreachable!("a pattern atom node maps to a host atom node"),
                    })
                    .collect();
                pattern.verify_overlays(host, atom_map)
            })
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

        // Stereo: a pattern stereo overlay matches iff the corresponding host site
        // bears a stereo element of the same class whose coset, reindexed from the
        // host ligand frame into the pattern's frame (via the atom correspondence),
        // is admitted by the pattern coset. An `Undetermined` pattern coset admits
        // both handednesses. TODO: a pattern that asserts stereo via `#T`/`#C` atom
        // /bond constraints rather than a `:stereo-atoms`/`:stereo-bonds` overlay is
        // not handled here — that needs the pattern run through stereo perception
        // (but not grounding, so no valence resolution).
        let mut stereo_atoms = Vec::new();
        for sp in pattern.stereo_atoms().iter() {
            let host_atom = atom_map[sp.site_id().index()];
            let sh = host.stereo_atoms().incident(host_atom).next()?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame = sp
                .ligand_frame()
                .into_iter()
                .map(|l| StereoLigand { atom_id: atom_map[l.atom_id.index()], kind: l.kind });
            let host_coset = sh.coset_for(frame)?;
            if !coset_matches(sp.coset(), &host_coset, sp.kind()) {
                return None;
            }
            stereo_atoms.push(sh.id);
        }

        let mut stereo_bonds = Vec::new();
        for sp in pattern.stereo_bonds().iter() {
            let [a, b] = pattern.bond(sp.site_id()).atom_ids();
            let host_edge = host
                .raw_graph()
                .find_edge(
                    NodeId::from(atom_map[a.index()]),
                    NodeId::from(atom_map[b.index()]),
                )
                .expect("a matched query bond maps to a host bond");
            let sh = host.bond(BondId::from(host_edge)).cis_trans_stereo()?;
            if sp.kind() != sh.kind() {
                return None;
            }
            let frame = sp
                .ligand_frame()
                .into_iter()
                .map(|l| StereoLigand { atom_id: atom_map[l.atom_id.index()], kind: l.kind });
            let host_coset = sh.coset_for(frame)?;
            if !coset_matches(sp.coset(), &host_coset, sp.kind()) {
                return None;
            }
            stereo_bonds.push(sh.id);
        }

        Some(MoleculeEmbedding::from_match(
            host,
            pattern,
            atom_map,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
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
    use super::super::molecule::MoleculeAst;
    use super::SubstructureMatchAlgorithm;
    use super::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};

    const SUBISO_ALGS: [SubgraphIsomorphismAlgorithm; 6] = [
        Vf2,
        Ullmann,
        Ri,
        ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH },
        Vf2Rdkit,
        RayKirsch,
    ];

    const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

    #[rstest]
    #[case::skeleton(
        mol!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]]
    )]
    #[case::noncovalent(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]] :noncovalent-bonds [{:atoms [0 2] :type "Hbd"}]}"#),
        mol!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#),
        vec![vec![AtomId(0), AtomId(2)], vec![AtomId(2), AtomId(0)]]
    )]
    #[case::noncovalent_absent(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        mol!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#),
        vec![]
    )]
    #[case::dative(
        mol!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1"}]}"#),
        mol!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1"}]}"#),
        vec![vec![AtomId(0), AtomId(1)]]
    )]
    #[case::dative_roles_swapped(
        mol!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donor 1 :acceptor 0 :type "1"}]}"#),
        mol!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donor 0 :acceptor 1 :type "1"}]}"#),
        vec![]
    )]
    #[case::stereo_chiral(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_enantiomer(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        vec![]
    )]
    #[case::stereo_agnostic_in_r(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_agnostic_in_s(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]]
    )]
    #[case::stereo_reframed(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th0"}]}"#),
        vec![vec![AtomId(0), AtomId(3), AtomId(2), AtomId(1)]]
    )]
    #[case::stereo_reframed_enantiomer(
        mol!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
        mol!(r#"{:atoms ["C #h1" "Br" "Cl" "F"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#),
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
                    .map(|e| e.host_atoms().to_vec())
                    .collect();
                occurrences.sort();
                assert_eq!(occurrences, expected, "{strategy:?}/{subiso:?}");
            }
        }
    }

}
