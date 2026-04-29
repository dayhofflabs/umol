//! Kekulize: aromatic-system form → Kekulé bond orders.
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

use std::collections::{HashMap, HashSet};

use umol_ast::ast::{
    AromaticSystemIdx, AtomConstraintKind, AtomIdx, BondConstraintKind, BondIdx, MoleculeAst,
    ValueAst,
};
use umol_graph_core::{EdgeId, Graph, NodeId, PerfectMatchingAlgorithm};

use crate::ops::transform::{Transformer, TransformerError};

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
pub struct Kekulize {
    model: KekulizationModel,
    node_order: Vec<AtomIdx>,
}

impl Kekulize {
    pub fn new(model: KekulizationModel, node_order: Vec<AtomIdx>) -> Self {
        Self { model, node_order }
    }
}

impl Transformer for Kekulize {
    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), TransformerError> {
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
                let bond = ast.bond_mut(bid).data;
                bond.order = ValueAst::Lit(2);
                bond.constraints
                    .retain(|c| c.kind() != BondConstraintKind::Aromatic);
            }
            for &bid in &plan.unmatched_bonds {
                let bond = ast.bond_mut(bid).data;
                bond.order = ValueAst::Lit(1);
                bond.constraints
                    .retain(|c| c.kind() != BondConstraintKind::Aromatic);
            }
            for &aidx in &plan.atoms {
                let atom = ast.atom_mut(aidx).data;
                atom.constraints.remove(AtomConstraintKind::AromaticValence);
            }
        }

        // Pass 2: drop the aromatic system entries via the builder.
        let to_remove: Vec<AromaticSystemIdx> = plans.iter().map(|p| p.system_idx).collect();
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
    system_idx: AromaticSystemIdx,
    atoms: Vec<AtomIdx>,
    matched_bonds: Vec<BondIdx>,
    unmatched_bonds: Vec<BondIdx>,
}

impl Kekulize {
    /// Build the per-system matching plan against an immutable AST.
    fn plan_systems(&self, ast: &MoleculeAst) -> Result<Vec<SystemPlan>, TransformerError> {
        let mut plans = Vec::with_capacity(ast.aromatic_systems().count());
        for view in ast.aromatic_systems().iter() {
            let system_idx = view.idx;
            let atoms: Vec<AtomIdx> = view.atoms().collect();
            let bonds: Vec<BondIdx> = view.bonds().collect();

            let atom_to_node: HashMap<AtomIdx, u32> = atoms
                .iter()
                .enumerate()
                .map(|(i, &a)| (a, i as u32))
                .collect();
            let local_edges: Vec<[u32; 2]> = bonds
                .iter()
                .map(|&bid| {
                    let bond = ast.bond(bid);
                    [atom_to_node[&bond.src], atom_to_node[&bond.tgt]]
                })
                .collect();
            let subgraph = Graph::new(atoms.len(), &local_edges);

            let local_order: Vec<NodeId> = self
                .node_order
                .iter()
                .filter_map(|aidx| atom_to_node.get(aidx).map(|&n| NodeId(n)))
                .collect();

            let matching = subgraph
                .perfect_matching(&local_order, self.model.algorithm)
                .ok_or(TransformerError::KekulizeNoMatching(system_idx))?;

            let matched_local: HashSet<EdgeId> = matching.edges().iter().copied().collect();
            let mut matched_bonds = Vec::new();
            let mut unmatched_bonds = Vec::new();
            for (i, &bid) in bonds.iter().enumerate() {
                if matched_local.contains(&EdgeId(i as u32)) {
                    matched_bonds.push(bid);
                } else {
                    unmatched_bonds.push(bid);
                }
            }

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
    use umol_ast::ast::{
        AromaticSystemAst, AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst,
        BondConstraint, Constraints, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;

    fn aromatic_carbon(pi: i64) -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.spin = SpinStateAst::closed_shell();
        atom.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(pi)),
        ));
        atom
    }

    fn benzene_aromatic() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| aromatic_carbon(1)).collect();
        let mut bonds = Vec::new();
        for i in 0..6 {
            let mut bond = BondAst::from_order(1);
            bond.constraints.add(BondConstraint::Aromatic);
            bonds.push((AtomIdx(i), AtomIdx((i + 1) % 6), bond));
        }
        let system = AromaticSystemAst::new(
            vec![ValueAst::Lit(1); 6],
            ValueAst::Lit(0),
            SpinStateAst::closed_shell(),
        );
        let aromatic_systems = vec![((0..6).map(AtomIdx).collect(), system)];
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            aromatic_systems,
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    fn ascending(n: u32) -> Vec<AtomIdx> {
        (0..n).map(AtomIdx).collect()
    }

    #[rstest]
    fn test_kekulize_benzene_assigns_alternating_orders() {
        let mut ast = benzene_aromatic();
        Kekulize::new(KekulizationModel::default(), ascending(6))
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 0);

        let orders: Vec<i64> = ast
            .bonds()
            .iter()
            .map(|view| match view.data.order {
                ValueAst::Lit(n) => n,
                _ => -1,
            })
            .collect();
        assert_eq!(orders.len(), 6);
        let doubles = orders.iter().filter(|&&o| o == 2).count();
        let singles = orders.iter().filter(|&&o| o == 1).count();
        assert_eq!(doubles, 3, "exactly 3 double bonds in a Kekulé benzene");
        assert_eq!(singles, 3, "exactly 3 single bonds in a Kekulé benzene");

        for view in ast.bonds().iter() {
            let has_aromatic = view
                .data
                .constraints
                .iter()
                .any(|c| c.kind() == BondConstraintKind::Aromatic);
            assert!(!has_aromatic, "Aromatic bond constraint must be cleared");
        }

        for view in ast.atoms().iter() {
            let has_av = view
                .data
                .constraints
                .iter()
                .any(|c| c.kind() == AtomConstraintKind::AromaticValence);
            assert!(!has_av, "AromaticValence constraint must be cleared");
        }
    }

    #[rstest]
    fn test_kekulize_no_aromatic_systems_is_noop() {
        let mut ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C); 2],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let original = ast.clone();
        Kekulize::new(KekulizationModel::default(), ascending(2))
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast, original);
    }

    #[rstest]
    fn test_kekulize_pentagon_returns_no_matching_error() {
        // Hypothetical aromatic system over an odd ring (5 atoms): no perfect
        // matching exists, so kekulize must error.
        let atoms: Vec<AtomAst> = (0..5).map(|_| aromatic_carbon(1)).collect();
        let mut bonds = Vec::new();
        for i in 0..5 {
            let mut bond = BondAst::from_order(1);
            bond.constraints.add(BondConstraint::Aromatic);
            bonds.push((AtomIdx(i), AtomIdx((i + 1) % 5), bond));
        }
        let system = AromaticSystemAst::new(
            vec![ValueAst::Lit(1); 5],
            ValueAst::Lit(0),
            SpinStateAst::closed_shell(),
        );
        let aromatic_systems = vec![((0..5).map(AtomIdx).collect(), system)];
        let mut ast = MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            aromatic_systems,
            vec![],
            vec![],
            Constraints::default(),
        );

        let result = Kekulize::new(KekulizationModel::default(), ascending(5))
            .transform_into(&mut ast);
        assert!(matches!(
            result,
            Err(TransformerError::KekulizeNoMatching(_))
        ));
    }
}
