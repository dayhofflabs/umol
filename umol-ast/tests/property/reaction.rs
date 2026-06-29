//! Property tests for reaction application and composition. Generates valid localized reactions:
//! deltas consistent with `lhs` (`ModifyField` with the lhs value as `old`), appended atoms, and
//! DPO-valid deletions (a removed atom takes all its incident bonds), so `apply` stays
//! dangling-free. Reactions drive the public surface only.

use proptest::bool::weighted;
use proptest::prelude::*;
use umol_ast::ast::{AtomDelta, BondDelta, CompositionScope, Delta, Deltas, ReactionAst};
use umol_graph_core::{EdgeId, SubgraphIsomorphismAlgorithm};

use crate::strategies::*;

const ALG: SubgraphIsomorphismAlgorithm = SubgraphIsomorphismAlgorithm::Vf2;

/// A small localized molecule: 1–4 element atoms over a simple edge set, bond orders 1–3.
fn simple_molecule_strategy() -> impl Strategy<Value = MoleculeAst> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                prop::collection::vec(
                    element_strategy().prop_map(AtomAst::from_element),
                    atom_count,
                ),
                edge_set_strategy(atom_count),
            )
        })
        .prop_flat_map(|(atoms, edges)| {
            let orders = prop::collection::vec(1u8..=3, edges.len());
            (Just(atoms), Just(edges), orders)
        })
        .prop_map(|(atoms, edges, orders)| {
            let bonds = edges
                .iter()
                .zip(orders)
                .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondAst::from_order(order)))
                .collect();
            MoleculeAst::from_atoms_and_bonds(atoms, bonds)
        })
}

/// A valid reaction over a generated `lhs`: DPO-valid atom deletions (each removed atom takes its
/// incident bonds), per-surviving-atom optional charge change and per-surviving-bond optional
/// order change (the `old` read from `lhs`, so apply's precondition holds), plus up to two new
/// atoms bonded to the lowest survivor. No dangling by construction.
fn reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    simple_molecule_strategy()
        .prop_flat_map(|lhs| {
            let atom_count = lhs.atoms().count();
            let bond_count = lhs.bonds().count();
            (
                Just(lhs),
                prop::collection::vec(weighted(0.25), atom_count),
                prop::collection::vec(prop::option::of(-2i64..=2), atom_count),
                prop::collection::vec(prop::option::of(1i64..=3), bond_count),
                prop::collection::vec(element_strategy(), 0..=2),
            )
        })
        .prop_map(|(lhs, removals, charges, orders, additions)| {
            build_reaction(lhs, removals, charges, orders, additions)
        })
}

fn build_reaction(
    lhs: MoleculeAst,
    removals: Vec<bool>,
    charges: Vec<Option<i64>>,
    orders: Vec<Option<i64>>,
    additions: Vec<Element>,
) -> ReactionAst {
    let atom_count = lhs.atoms().count();
    let bond_count = lhs.bonds().count();
    let removed_atoms: HashSet<AtomId> = removals
        .iter()
        .enumerate()
        .filter(|&(_, &remove)| remove)
        .map(|(index, _)| AtomId(index as u32))
        .collect();
    // A removed atom takes all its incident bonds with it (DPO-valid; apply never dangles).
    let mut removed_bonds: HashSet<BondId> = HashSet::new();
    for j in 0..bond_count as u32 {
        let [x, y] = lhs.raw_graph().edge_endpoints(EdgeId(j));
        if removed_atoms.contains(&AtomId::from(x)) || removed_atoms.contains(&AtomId::from(y)) {
            removed_bonds.insert(BondId(j));
        }
    }

    let mut deltas: Vec<Delta> = Vec::new();
    for &id in &removed_atoms {
        deltas.push(Delta::Atom(AtomDelta::Remove {
            id,
            ast: lhs.atom(id).ast.clone(),
        }));
    }
    for &id in &removed_bonds {
        let [x, y] = lhs.raw_graph().edge_endpoints(EdgeId(id.0));
        deltas.push(Delta::Bond(BondDelta::Remove {
            id,
            atoms: [AtomId::from(x), AtomId::from(y)],
            ast: lhs.bond(id).ast.clone(),
        }));
    }
    for (index, new_charge) in charges.into_iter().enumerate() {
        let id = AtomId(index as u32);
        if removed_atoms.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.atom(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::Atom(AtomDelta::ModifyField {
                id,
                change: AtomFieldChange::Charge { old, new },
            }));
        }
    }
    for (index, new_order) in orders.into_iter().enumerate() {
        let id = BondId(index as u32);
        if removed_bonds.contains(&id) {
            continue;
        }
        let Some(order) = new_order else { continue };
        let old = lhs.bond(id).ast.order.clone();
        let new = ValueAst::Lit(order);
        if old != new {
            deltas.push(Delta::Bond(BondDelta::ModifyField {
                id,
                change: BondFieldChange::Order { old, new },
            }));
        }
    }
    // Append atoms bonded to the lowest surviving atom (isolated if every atom is removed).
    let anchor = (0..atom_count as u32)
        .map(AtomId)
        .find(|id| !removed_atoms.contains(id));
    for (offset, element) in additions.into_iter().enumerate() {
        let atom = AtomId((atom_count + offset) as u32);
        deltas.push(Delta::Atom(AtomDelta::Add {
            id: atom,
            ast: AtomAst::from_element(element),
        }));
        if let Some(anchor) = anchor {
            deltas.push(Delta::Bond(BondDelta::Add {
                id: BondId((bond_count + offset) as u32),
                atoms: [anchor, atom],
                ast: BondAst::from_order(1),
            }));
        }
    }
    ReactionAst::new(lhs, Deltas::from_iter(deltas))
}

proptest! {
    /// Applying a reaction at the identity occurrence of its own `lhs` reproduces the span's
    /// `right()` — the `transact`-apply path agrees with the span projection.
    #[test]
    fn test_reaction_ast_apply_reproduces_right(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let right = span.right();
            prop_assert!(reaction.apply(&reaction.lhs, ALG).any(|product| product == right));
        }
    }

    /// `reverse` swaps the span's sides. The reverse reaction's reactant is *exactly* the forward
    /// product. Its product is the forward reactant only up to atom renumbering (re-added atoms
    /// append rather than reoccupy their original ids), so structurally we check the reconstructed
    /// size — the exact value is covered by the unit tests on fixed frames.
    #[test]
    fn test_reaction_ast_reverse_swaps_sides(reaction in reaction_strategy()) {
        if let (Ok(span), Ok(reverse)) = (reaction.to_reaction_span(), reaction.reverse()) {
            if let Ok(reverse_span) = reverse.to_reaction_span() {
                prop_assert_eq!(reverse_span.left(), span.right());
                let forward_reactant = span.left();
                let reverse_product = reverse_span.right();
                prop_assert_eq!(
                    reverse_product.atoms().count(),
                    forward_reactant.atoms().count()
                );
                prop_assert_eq!(
                    reverse_product.bonds().count(),
                    forward_reactant.bonds().count()
                );
            }
        }
    }

    /// Every composite is itself a valid reaction: applying it at its own `lhs` reproduces its
    /// `right()`. Catches frame-algebra errors in the composite construction.
    #[test]
    fn test_reaction_ast_compose_well_formed(
        a in reaction_strategy(),
        b in reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            if let Ok(span) = composite.to_reaction_span() {
                let right = span.right();
                prop_assert!(composite.apply(&composite.lhs, ALG).any(|product| product == right));
            }
        }
    }

    /// Soundness: every product of a composite applied to `A`'s reactant is also a product of
    /// applying B after A — `compose` invents no reactions.
    #[test]
    fn test_reaction_ast_compose_sound(
        a in reaction_strategy(),
        b in reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composites = a.compose(&b, CompositionScope::Full);
        let composed: Vec<MoleculeAst> = composites
            .iter()
            .flat_map(|composite| composite.apply(&host, ALG))
            .collect();

        let intermediates: Vec<MoleculeAst> = a.apply(&host, ALG).collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(b.apply(intermediate, ALG));
        }

        for product in &composed {
            prop_assert!(sequential.contains(product));
        }
    }

    /// The reaction round-trips through the EDN surface: render → parse reaches a
    /// fixpoint, exercising the atom/bond add / remove / modify-field delta ops
    /// (`ReactionAst::to_edn` then `from_edn`, twice, must agree).
    #[test]
    fn test_reaction_ast_edn_roundtrip_stable(reaction in reaction_strategy()) {
        let once = ReactionAst::from_edn(&reaction.to_edn())
            .map_err(|e| TestCaseError::fail(format!("first reparse failed: {e}")))?;
        let twice = ReactionAst::from_edn(&once.to_edn())
            .map_err(|e| TestCaseError::fail(format!("second reparse failed: {e}")))?;
        prop_assert_eq!(once, twice);
    }
}
