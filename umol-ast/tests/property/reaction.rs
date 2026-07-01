//! Property tests for reaction application and composition. Generates valid localized reactions:
//! deltas consistent with `lhs` (`ModifyField` with the lhs value as `old`), appended atoms, and
//! DPO-valid deletions (a removed atom takes all its incident bonds), so `apply` stays
//! dangling-free. Reactions drive the public surface only.

use proptest::bool::weighted;
use proptest::prelude::*;
use umol_ast::ast::{
    AromaticSystemDelta, AtomDelta, BondDelta, CompositionScope, DativeBondDelta, Delta, Deltas,
    DpoValidator, MulticenterBondDelta, NoncovalentBondDelta, ReactionAst,
};
use umol_graph_core::{EdgeId, SubgraphIsomorphismAlgorithm};
use umol_utils::solution::Solution;

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

fn reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(simple_molecule_strategy())
}

/// A localized molecule with DAMN overlays (dative / aromatic / multicenter / noncovalent) — no
/// stereo (stereo reaction deltas are I6, and compose bails on a stereo lhs) and no molecule
/// constraints (orthogonal). 1–4 atoms; overlays generated as in `molecule_ast_strategy`, scoped.
fn overlay_molecule_strategy() -> impl Strategy<Value = MoleculeAst> {
    (1usize..=4)
        .prop_flat_map(|atom_count| {
            (
                Just(atom_count),
                prop::collection::vec(
                    element_strategy().prop_map(AtomAst::from_element),
                    atom_count,
                ),
                edge_set_strategy(atom_count),
            )
        })
        .prop_flat_map(|(atom_count, atoms, edges)| {
            let orders = prop::collection::vec(1u8..=3, edges.len());
            let datives = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    dative_bond_strategy(),
                ),
                0..=1,
            );
            let aromatics = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), aromatic_system_ast_for(n))
                    },
                ),
                0..=1,
            );
            let multicenters = prop::collection::vec(
                distinct_atoms_strategy(atom_count, 3, 4.min(atom_count.max(3))).prop_flat_map(
                    |atoms| {
                        let n = atoms.len();
                        (Just(atoms), multicenter_bond_ast_for(n))
                    },
                ),
                0..=1,
            );
            let noncovalents = prop::collection::vec(
                (
                    distinct_atoms_strategy(atom_count, 2, 2),
                    noncovalent_bond_ast_strategy(),
                ),
                0..=1,
            );
            (
                Just(atoms),
                Just(edges),
                orders,
                datives,
                aromatics,
                multicenters,
                noncovalents,
            )
        })
        .prop_map(
            |(atoms, edges, orders, datives, aromatics, multicenters, noncovalents)| {
                let bonds = edges
                    .iter()
                    .zip(orders)
                    .map(|(&[a, b], order)| (AtomId(a), AtomId(b), BondAst::from_order(order)))
                    .collect();
                let dative = datives
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((vec![*a], *b, data)),
                        _ => None,
                    })
                    .collect();
                let aromatic = aromatics
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let multicenter = multicenters
                    .into_iter()
                    .filter(|(atoms, _)| atoms.len() >= 3)
                    .collect();
                let noncovalent = noncovalents
                    .into_iter()
                    .filter_map(|(atoms, data)| match atoms.as_slice() {
                        [a, b] if a != b => Some((*a, *b, data)),
                        _ => None,
                    })
                    .collect();
                MoleculeAst::from_parts(
                    atoms,
                    bonds,
                    dative,
                    aromatic,
                    multicenter,
                    noncovalent,
                    vec![],
                    vec![],
                    Constraints::new(),
                )
            },
        )
}

/// A reaction whose `lhs` carries DAMN overlays — exercises overlay carry, correspondence, and
/// co-deletion through compose.
fn overlay_reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(overlay_molecule_strategy())
}

/// A valid reaction over any generated `lhs`: DPO-valid atom deletions (each removed atom takes its
/// incident bonds and overlays), per-surviving-atom optional charge change and per-surviving-bond
/// optional order change (the `old` read from `lhs`, so apply's precondition holds), plus up to two
/// new atoms bonded to the lowest survivor. No dangling by construction.
fn reaction_over(molecule: impl Strategy<Value = MoleculeAst>) -> impl Strategy<Value = ReactionAst> {
    molecule
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
    // A removed atom also takes its incident overlays (DPO-valid; apply never dangles on overlays).
    let mut removed_dative: HashSet<DativeBondId> = HashSet::new();
    let mut removed_aromatic: HashSet<AromaticSystemId> = HashSet::new();
    let mut removed_multicenter: HashSet<MulticenterBondId> = HashSet::new();
    let mut removed_noncovalent: HashSet<NoncovalentBondId> = HashSet::new();
    for &id in &removed_atoms {
        let view = lhs.atom(id);
        removed_dative.extend(view.dative_bond_ids());
        if let Some(system) = view.aromatic_system_id() {
            removed_aromatic.insert(system);
        }
        removed_multicenter.extend(view.multicenter_bond_ids());
        removed_noncovalent.extend(view.noncovalent_bond_ids());
    }
    for &id in &removed_dative {
        let view = lhs.dative_bond(id);
        deltas.push(Delta::DativeBond(DativeBondDelta::Remove {
            id,
            donors: view.donor_ids().collect(),
            acceptor: view.acceptor_id(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_aromatic {
        let view = lhs.aromatic_system(id);
        deltas.push(Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_multicenter {
        let view = lhs.multicenter_bond(id);
        deltas.push(Delta::MulticenterBond(MulticenterBondDelta::Remove {
            id,
            atoms: view.atom_ids().collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_noncovalent {
        let view = lhs.noncovalent_bond(id);
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
            id,
            atoms: view.atom_ids(),
            ast: view.ast.clone(),
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

// Overlay-bearing reactions (DAMN lhs, DPO-valid): the compose properties over the overlay carry /
// correspondence / co-deletion machinery. `overlay_reaction_strategy` subsumes the atom/bond case
// (overlay counts are 0..=1). Completeness (Full: sequential ⊆ composed) is a separate follow-on.
proptest! {
    /// Isolation probe: a plain overlay reaction's `apply` at its own `lhs` reproduces its
    /// `right()`. If this fails, the discrepancy is in apply-vs-span for overlays, not compose.
    #[test]
    fn test_reaction_ast_apply_reproduces_right_overlay(reaction in overlay_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let right = span.right();
            prop_assert!(reaction.apply(&reaction.lhs, ALG).any(|product| product == right));
        }
    }

    /// Every composite of two overlay reactions is a valid reaction: applying it at its own `lhs`
    /// reproduces its `right()`. Catches overlay frame-algebra errors in the composite, and (the
    /// reason it once failed) `apply_at` removing multiple same-kind overlays: composites routinely
    /// remove ≥2 overlays of one kind, which the pre-batching single-id lowering mishandled.
    #[test]
    fn test_reaction_ast_compose_well_formed_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            if let Ok(span) = composite.to_reaction_span() {
                let right = span.right();
                prop_assert!(composite.apply(&composite.lhs, ALG).any(|product| product == right));
            }
        }
    }

    /// Soundness with overlays: every product of a composite applied to A's reactant is also a
    /// product of applying B after A — compose invents no reactions, overlays included.
    #[test]
    fn test_reaction_ast_compose_sound_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
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

    /// P1 completeness (`Full`): every sequential product (B applied after A) is also some
    /// composite's product. Together with `compose_sound_overlay` (composed ⊆ seq) this is set
    /// equality at `host = a.lhs`.
    ///
    /// IGNORED — completeness needs three coordinated compose fixes, none landed:
    /// (1) **monomorphism overlaps** — compose enumerates *induced* common subgraphs but `apply` is a
    ///     monomorphism, so an overlap where R_A carries a bond/overlay L_B lacks (context in A's
    ///     product) is dropped, missing the sequential product B would produce over that structure;
    /// (2) **`meet` interface** — `lhs_c`'s overlap entity is A's lhs only, so B's specificity there
    ///     is lost unless (3) preserves it; it must be `meet(A-lhs, B-lhs)`;
    /// (3) **delta rebasing** — B's overlap-entity delta carries B's lhs old-state, but in the
    ///     composite it acts on A's product state, so the `Remove` ast / `ModifyField` `old` must be
    ///     R_A's value for `canonicalize` to fold A's modify with B's op.
    /// (3) alone is unsound (it drops B's constraints without (2)); they must land together.
    #[test]
    #[ignore = "compose completeness needs monomorphism overlaps + meet interface + delta rebasing together"]
    fn test_reaction_ast_compose_complete_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|composite| composite.apply(&host, ALG))
            .collect();

        let intermediates: Vec<MoleculeAst> = a.apply(&host, ALG).collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(b.apply(intermediate, ALG));
        }

        for product in &sequential {
            prop_assert!(
                composed.contains(product),
                "sequential product missing from composed set (P1 completeness)"
            );
        }
    }

    /// Every composite is DPO-valid — no deleted atom leaves a dangling bond or overlay. Confirms
    /// the compose during-check yields dangling-free composites (via the tier-2 `DpoValidator`).
    #[test]
    fn test_reaction_ast_compose_dangling_free(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            prop_assert_eq!(
                DpoValidator.validate_reaction(&composite).unwrap(),
                Solution::Determined(())
            );
        }
    }

    /// Reaction ↔ span roundtrip fidelity: recovering the reaction from a span and re-materializing
    /// reproduces the span (`to_reaction` then `to_reaction_span` is the identity on spans),
    /// exercising the overlay `EntitySpan` columns and the span→delta recovery in both directions.
    #[test]
    fn test_reaction_ast_span_roundtrip(reaction in overlay_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            if let Ok(rebuilt) = span.to_reaction().to_reaction_span() {
                prop_assert_eq!(rebuilt, span);
            }
        }
    }

    /// `RcAnchored` is a sound filter: every RC-anchored composite is also a `Full` composite.
    #[test]
    fn test_reaction_ast_compose_rc_anchored_subset(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let full = a.compose(&b, CompositionScope::Full);
        for composite in a.compose(&b, CompositionScope::RcAnchored) {
            prop_assert!(full.contains(&composite));
        }
    }

    /// P4 — determinism: `compose` returns the identical `Vec` on repeated calls, and is invariant
    /// under pre-canonicalizing the inputs (compose canonicalizes the deltas itself).
    #[test]
    fn test_reaction_ast_compose_determinism(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        prop_assert_eq!(
            a.compose(&b, CompositionScope::Full),
            a.compose(&b, CompositionScope::Full)
        );
        prop_assert_eq!(
            a.compose(&b, CompositionScope::RcAnchored),
            a.compose(&b, CompositionScope::RcAnchored)
        );
        if let (Ok(ac), Ok(bc)) = (a.clone().canonicalize(), b.clone().canonicalize()) {
            prop_assert_eq!(
                a.compose(&b, CompositionScope::Full),
                ac.compose(&bc, CompositionScope::Full)
            );
        }
    }

    /// P3 — every composite's deltas are in canonical normal form.
    #[test]
    fn test_reaction_ast_compose_canonical_deltas(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for c in a.compose(&b, CompositionScope::Full) {
            let canonical = c
                .deltas
                .clone()
                .canonicalize()
                .map_err(|e| TestCaseError::fail(format!("composite deltas not canonical: {e:?}")))?;
            prop_assert_eq!(canonical, c.deltas);
        }
    }

    /// P6 — no parallel overlays: within each kind a composite's overlays have distinct participant
    /// sets, so correspondence reuses an id and never duplicates (spec §4.1).
    #[test]
    fn test_reaction_ast_compose_distinct_overlays(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for c in a.compose(&b, CompositionScope::Full) {
            let m = &c.lhs;

            let mut dative: Vec<(Vec<AtomId>, AtomId)> = m
                .dative_bonds()
                .iter()
                .map(|x| {
                    let mut donors: Vec<AtomId> = x.donor_ids().collect();
                    donors.sort();
                    (donors, x.acceptor_id())
                })
                .collect();
            let dative_count = dative.len();
            dative.sort();
            dative.dedup();
            prop_assert_eq!(dative.len(), dative_count, "duplicate dative bonds");

            let mut aromatic: Vec<Vec<AtomId>> = m
                .aromatic_systems()
                .iter()
                .map(|x| {
                    let mut v: Vec<AtomId> = x.atom_ids().collect();
                    v.sort();
                    v
                })
                .collect();
            let aromatic_count = aromatic.len();
            aromatic.sort();
            aromatic.dedup();
            prop_assert_eq!(aromatic.len(), aromatic_count, "duplicate aromatic systems");

            let mut multicenter: Vec<Vec<AtomId>> = m
                .multicenter_bonds()
                .iter()
                .map(|x| {
                    let mut v: Vec<AtomId> = x.atom_ids().collect();
                    v.sort();
                    v
                })
                .collect();
            let multicenter_count = multicenter.len();
            multicenter.sort();
            multicenter.dedup();
            prop_assert_eq!(
                multicenter.len(),
                multicenter_count,
                "duplicate multicenter bonds"
            );

            let mut noncovalent: Vec<[AtomId; 2]> = m
                .noncovalent_bonds()
                .iter()
                .map(|x| {
                    let mut p = x.atom_ids();
                    p.sort();
                    p
                })
                .collect();
            let noncovalent_count = noncovalent.len();
            noncovalent.sort();
            noncovalent.dedup();
            prop_assert_eq!(
                noncovalent.len(),
                noncovalent_count,
                "duplicate noncovalent bonds"
            );
        }
    }
}
