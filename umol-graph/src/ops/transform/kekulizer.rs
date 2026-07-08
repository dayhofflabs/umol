//! Kekulizer: aromatic-system form → Kekulé bond orders.
//!
//! For each aromatic system in the input, finds a perfect matching on the
//! induced subgraph (atoms in the system, bonds between them). Matched bonds
//! become order 2; unmatched bonds become order 1. The aromatic constraint is
//! removed from the system's bonds and from the system's atoms; the system
//! entry itself is removed at the end.
//!
//! The matching algorithm is configurable via [`KekulizationModel`]; the
//! atom processing order — which controls determinism — is fixed at
//! construction time and is the caller's responsibility (e.g., from a
//! nauty/Traces canonical labeling).
//!
//! TODO: Expand to charged systems.

use std::collections::{HashMap, HashSet};

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemId, AtomConstraintKey, AtomId, BondConstraintKey, BondId, MoleculeAst, ValueAst,
};
use umol_graph_core::PerfectMatchingAlgorithm;

use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum KekulizerError {
    #[error("no perfect matching exists for aromatic system {0:?}")]
    NoMatching(AromaticSystemId),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KekulizationModel {
    pub algorithm: PerfectMatchingAlgorithm,
}

impl KekulizationModel {
    pub fn new(algorithm: PerfectMatchingAlgorithm) -> Self {
        Self { algorithm }
    }
}

impl Default for KekulizationModel {
    fn default() -> Self {
        Self::new(PerfectMatchingAlgorithm::BacktrackingDfs)
    }
}

#[derive(Clone, Debug)]
pub struct Kekulizer {
    model: KekulizationModel,
    node_order: Vec<AtomId>,
}

impl Kekulizer {
    pub fn new(model: KekulizationModel, node_order: Vec<AtomId>) -> Self {
        Self { model, node_order }
    }
}

impl Transformer for Kekulizer {
    type Error = KekulizerError;

    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), KekulizerError> {
        if ast.aromatic_systems().count() == 0 {
            return Ok(());
        }

        // Plan the per-system matching against an immutable AST snapshot, then
        // apply the bond-order writes and structural cleanup in passes that
        // require &mut.
        let plans = self.plan_systems(ast)?;

        // Pass 1: bond-order writes and Aromatic-constraint stripping.
        for plan in &plans {
            for &bid in &plan.matched_bonds {
                let bond = ast.bond_mut(bid).ast;
                bond.order = ValueAst::Lit(2);
                bond.constraints.remove(BondConstraintKey::Aromatic);
            }
            for &bid in &plan.unmatched_bonds {
                let bond = ast.bond_mut(bid).ast;
                bond.order = ValueAst::Lit(1);
                bond.constraints.remove(BondConstraintKey::Aromatic);
            }
            for &aidx in &plan.atoms {
                let atom = ast.atom_mut(aidx).ast;
                atom.constraints.remove(AtomConstraintKey::AromaticValence);
            }
        }

        // Pass 2: drop the aromatic system entries via the builder.
        let to_remove: Vec<AromaticSystemId> = plans.iter().map(|p| p.system_idx).collect();
        let mut builder = ast.edit();
        builder.remove_aromatic_systems(&to_remove);
        *ast = builder.build();

        Ok(())
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a> {
        Box::new(self.transform(ast).ok().into_iter())
    }
}

struct SystemPlan {
    system_idx: AromaticSystemId,
    atoms: Vec<AtomId>,
    matched_bonds: Vec<BondId>,
    unmatched_bonds: Vec<BondId>,
}

impl Kekulizer {
    /// Build the per-system matching plan against an immutable AST.
    fn plan_systems(&self, ast: &MoleculeAst) -> Result<Vec<SystemPlan>, KekulizerError> {
        let mut plans = Vec::with_capacity(ast.aromatic_systems().count());
        for view in ast.aromatic_systems().iter() {
            let system_idx = view.id;
            let atoms: Vec<AtomId> = view.atom_ids().collect();
            let bonds: Vec<BondId> = view.bond_ids().collect();

            let subgraph = ast.induced_subgraph(&atoms);
            let extracted = ast.extract(&subgraph);

            // extracted has atoms in host-id order; build host→sub map.
            let mut sorted_host: Vec<AtomId> = subgraph
                .atoms()
                .mates()
                .iter()
                .map(|&(_, host)| AtomId::from(host))
                .collect();
            sorted_host.sort_unstable();
            let host_to_sub: HashMap<AtomId, AtomId> = sorted_host
                .iter()
                .enumerate()
                .map(|(i, &h)| (h, AtomId(i as u32)))
                .collect();

            let sub_order: Vec<AtomId> = self
                .node_order
                .iter()
                .copied()
                .filter_map(|a| host_to_sub.get(&a).copied())
                .collect();

            // Bonds in extracted preserve host-bond-id order, matching the sub→host bond images.
            let host_bonds: Vec<BondId> = subgraph
                .bonds()
                .mates()
                .iter()
                .map(|&(_, host)| host)
                .collect();
            let matched: HashSet<BondId> = extracted
                .graph()
                .perfect_matching(&sub_order, self.model.algorithm)
                .ok_or(KekulizerError::NoMatching(system_idx))?
                .bonds()
                .map(|sub| host_bonds[sub.index()])
                .collect();
            let (matched_bonds, unmatched_bonds): (Vec<BondId>, Vec<BondId>) =
                bonds.iter().copied().partition(|b| matched.contains(b));

            plans.push(SystemPlan {
                system_idx,
                atoms,
                matched_bonds,
                unmatched_bonds,
            });
        }
        Ok(plans)
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{AromaticSystemId, AtomId, MoleculeAst};
    use umol_ast::mol_ground;

    use super::*;

    #[rstest]
    #[case::benzene(
        mol_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]}"#),
        mol_ground!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
    )]
    fn test_kekulizer_transform_into(#[case] input: MoleculeAst, #[case] expected: MoleculeAst) {
        let mut ast = input;
        Kekulizer::new(KekulizationModel::default(), (0..6).map(AtomId).collect())
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast, expected);
    }

    #[rstest]
    #[case::kekule_benzene( mol_ground!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#))]
    fn test_kekulizer_transform_into_identity(#[case] input: MoleculeAst) {
        let mut ast = input.clone();
        Kekulizer::new(KekulizationModel::default(), (0..6).map(AtomId).collect())
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast, input);
    }

    #[rstest]
    #[case::no_matching( mol_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]"}]}"#))]
    fn test_kekulizer_transform_into_error(#[case] input: MoleculeAst) {
        let mut ast = input;
        let result = Kekulizer::new(KekulizationModel::default(), (0..5).map(AtomId).collect())
            .transform_into(&mut ast);
        assert_eq!(result, Err(KekulizerError::NoMatching(AromaticSystemId(0))));
    }
}
