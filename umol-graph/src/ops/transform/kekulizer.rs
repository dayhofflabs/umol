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

use std::collections::{HashMap, HashSet};

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemId, AtomConstraintKind, AtomId, BondConstraintKind, BondId, MoleculeAst, ValueAst,
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
                bond.constraints
                    .retain(|c| c.kind() != BondConstraintKind::Aromatic);
            }
            for &bid in &plan.unmatched_bonds {
                let bond = ast.bond_mut(bid).ast;
                bond.order = ValueAst::Lit(1);
                bond.constraints
                    .retain(|c| c.kind() != BondConstraintKind::Aromatic);
            }
            for &aidx in &plan.atoms {
                let atom = ast.atom_mut(aidx).ast;
                atom.constraints.remove(AtomConstraintKind::AromaticValence);
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
    use umol_ast::ast::{
        AromaticSystemAst, AromaticValenceAst, AtomAst, AtomConstraint, AtomId, BondAst,
        BondConstraint, BooleanAst, Constraints, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_chem::element::Element;

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
            bond.constraints
                .add(BondConstraint::Aromatic(BooleanAst::Lit(true)));
            bonds.push((AtomId(i), AtomId((i + 1) % 6), bond));
        }
        let system = AromaticSystemAst::from_electrons(vec![1; 6]);
        let aromatic_systems = vec![((0..6).map(AtomId).collect(), system)];
        MoleculeAst::from_parts(
            atoms,
            bonds,
            vec![],
            aromatic_systems,
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    fn ascending(n: u32) -> Vec<AtomId> {
        (0..n).map(AtomId).collect()
    }

    #[rstest]
    fn test_kekulize_benzene_assigns_alternating_orders() {
        let mut ast = benzene_aromatic();
        Kekulizer::new(KekulizationModel::default(), ascending(6))
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast.aromatic_systems().count(), 0);

        let orders: Vec<i64> = ast
            .bonds()
            .iter()
            .map(|view| match view.ast.order {
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
                .ast
                .constraints
                .iter()
                .any(|c| c.kind() == BondConstraintKind::Aromatic);
            assert!(!has_aromatic, "Aromatic bond constraint must be cleared");
        }

        for view in ast.atoms().iter() {
            let has_av = view
                .ast
                .constraints
                .iter()
                .any(|c| c.kind() == AtomConstraintKind::AromaticValence);
            assert!(!has_av, "AromaticValence constraint must be cleared");
        }
    }

    #[rstest]
    fn test_kekulize_no_aromatic_systems_is_noop() {
        let mut ast = MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C); 2],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );
        let original = ast.clone();
        Kekulizer::new(KekulizationModel::default(), ascending(2))
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
            bond.constraints
                .add(BondConstraint::Aromatic(BooleanAst::Lit(true)));
            bonds.push((AtomId(i), AtomId((i + 1) % 5), bond));
        }
        let system = AromaticSystemAst::from_electrons(vec![1; 5]);
        let aromatic_systems = vec![((0..5).map(AtomId).collect(), system)];
        let mut ast = MoleculeAst::from_parts(
            atoms,
            bonds,
            vec![],
            aromatic_systems,
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        );

        let result =
            Kekulizer::new(KekulizationModel::default(), ascending(5)).transform_into(&mut ast);
        assert!(matches!(result, Err(KekulizerError::NoMatching(_))));
    }
}
