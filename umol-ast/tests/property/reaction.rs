//! Property tests for reaction application and composition. Generates valid localized reactions:
//! deltas consistent with `lhs` (`ModifyField` with the lhs value as `old`), appended atoms, and
//! DPO-valid deletions (a removed atom takes all its incident bonds), so `apply` stays
//! dangling-free. Reactions drive the public surface only.
//!
//! ## Comparison contract
//!
//! | Operation | Property-test use |
//! | --- | --- |
//! | `==` | Exact stored representation in one fixed ID and participant frame. |
//! | entity `canonical_eq` | Semantic equality of one entity AST; no topology or incidence. |
//! | relation `equiv` / `equiv_under` | Canonical entity-data equality in the same or a supplied participant frame. |
//! | molecule `equiv` / `equiv_under` | Complete semantic equality of topology, incidence, entity data, and constraints in the same or a supplied ID/participant frame. |
//!
//! Application results below are complete molecules in the host's existing frame, so they use
//! `MoleculeAst::equiv`. Properties comparing a derived molecule through a non-identity match use
//! `MoleculeAst::equiv_under`; exact `==` remains appropriate only when representation identity is
//! itself the invariant.

use proptest::bool::weighted;
use proptest::prelude::*;
use umol_ast::ast::{
    AromaticSystemDelta, AromaticSystemFieldChange, AtomDelta, BondDelta, Canonicalize,
    CompositionScope, DativeBondConstraintAst, DativeBondDelta, DativeBondFieldChange, Delta,
    Deltas, DpoValidator, MoleculeParts, MulticenterBondDelta, MulticenterBondFieldChange,
    NoncovalentBondAst, NoncovalentBondDelta, NoncovalentBondId, NoncovalentBondKind,
    NoncovalentBondKindAst, ReactionAst, ReactionSpanAst, StereoAtomAst, StereoAtomDelta,
    StereoAtomFieldChange, StereoBondAst, StereoBondDelta, StereoBondFieldChange,
    StereoConfigurationAst, StereoKind, StereoLigand,
};
use umol_graph_core::{EdgeId, SubgraphIsomorphismAlgorithm};
use umol_perm::Permutation;
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
            MoleculeAst::from_parts(MoleculeParts {
                atoms,
                bonds,
                ..Default::default()
            })
        })
}

fn reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(simple_molecule_strategy())
}

/// A localized molecule with DAMN overlays (dative / aromatic / multicenter / noncovalent) plus
/// stereo (tetrahedral atoms / cis-trans bonds) and no molecule constraints (orthogonal). 1–4 atoms;
/// overlays generated as in `molecule_ast_strategy`, scoped.
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
            // A tetrahedral stereo atom: a site atom plus four ligands. Real atoms fill the first
            // slots (ids need not be graph neighbors — tier-1 only requires the kind's ligand
            // count); virtual implicit-H / lone-pair fills pad to `degree == 4`, all bearing the
            // site atom. 0..=1 so many molecules have none.
            let stereo_atoms = stereo_atom_overlay_strategy(atom_count);
            // A cis/trans stereo bond: a bond as site plus two ligand atoms (padded with virtual
            // fills to `degree == 4`). Requires a bond to name as site.
            let stereo_bonds = if edges.is_empty() {
                Just(Vec::new()).boxed()
            } else {
                stereo_bond_overlay_strategy(atom_count, edges.len())
            };
            (
                Just(atoms),
                Just(edges),
                orders,
                datives,
                aromatics,
                multicenters,
                noncovalents,
                stereo_atoms,
                stereo_bonds,
            )
        })
        .prop_map(
            |(
                atoms,
                edges,
                orders,
                datives,
                aromatics,
                multicenters,
                noncovalents,
                stereo_atoms,
                stereo_bonds,
            )| {
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
                MoleculeAst::from_parts(MoleculeParts {
                    atoms,
                    bonds,
                    dative,
                    aromatic,
                    multicenter,
                    noncovalent,
                    stereo_atoms,
                    stereo_bonds,
                    constraints: Constraints::new(),
                })
            },
        )
}

/// Cosets valid for `kind`: `Undetermined` or an in-range `Lit` index (`0..kind.count()`). Relative
/// reaction ops (`Swap` / `Mirror` / `Apply`) act on the coset through the kind's algebra, which
/// panics on an out-of-range index — so unlike the generic `stereo_coset_strategy`, indices are
/// bounded by the kind's coset count.
fn stereo_coset_for_kind(kind: StereoKind) -> impl Strategy<Value = StereoCosetAst> {
    let count = kind.count() as u32;
    prop_oneof![
        Just(StereoCosetAst::Undetermined),
        (0..count).prop_map(StereoCosetAst::Lit),
    ]
}

fn aromatic_system_update_for(atom_count: usize) -> impl Strategy<Value = AromaticSystemUpdate> {
    (
        prop::option::of(prop_oneof![
            Just(ElectronCountsAst::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsAst::Lit),
        ]),
        prop::option::of(value_basic(-2..=2)),
        spin_state_update_strategy(),
        aromatic_system_update_constraints_strategy(),
    )
        .prop_map(
            |(electrons, charge, spin, constraints)| AromaticSystemUpdate {
                electrons,
                charge,
                spin,
                constraints,
            },
        )
}

fn multicenter_bond_update_for(atom_count: usize) -> impl Strategy<Value = MulticenterBondUpdate> {
    (
        prop::option::of(prop_oneof![
            Just(ElectronCountsAst::Undetermined),
            prop::collection::vec(0i64..=2, atom_count).prop_map(ElectronCountsAst::Lit),
        ]),
        prop::option::of(value_basic(-2..=2)),
        spin_state_update_strategy(),
        multicenter_bond_update_constraints_strategy(),
    )
        .prop_map(
            |(electrons, charge, spin, constraints)| MulticenterBondUpdate {
                electrons,
                charge,
                spin,
                constraints,
            },
        )
}

fn stereo_atom_application_update_strategy() -> impl Strategy<Value = StereoAtomUpdate> {
    (
        prop_oneof![
            Just(StereoConfigurationUpdate::Unchanged),
            Just(StereoConfigurationUpdate::Undetermined),
            prop::option::of(stereo_coset_for_kind(StereoKind::Tetrahedral)).prop_map(|coset| {
                StereoConfigurationUpdate::Kinded {
                    kind: StereoKind::Tetrahedral,
                    coset,
                }
            }),
        ],
        stereo_atom_update_constraints_strategy(StereoKind::Tetrahedral),
    )
        .prop_map(|(configuration, constraints)| StereoAtomUpdate {
            configuration,
            constraints,
        })
}

fn stereo_bond_application_update_strategy() -> impl Strategy<Value = StereoBondUpdate> {
    (
        prop_oneof![
            Just(StereoConfigurationUpdate::Unchanged),
            Just(StereoConfigurationUpdate::Undetermined),
            prop::option::of(stereo_coset_for_kind(StereoKind::CisTrans)).prop_map(|coset| {
                StereoConfigurationUpdate::Kinded {
                    kind: StereoKind::CisTrans,
                    coset,
                }
            }),
        ],
        stereo_bond_update_constraints_strategy(StereoKind::CisTrans),
    )
        .prop_map(|(configuration, constraints)| StereoBondUpdate {
            configuration,
            constraints,
        })
}

/// A `degree`-length ligand frame of *distinct* `StereoLigand`s over `atom_count` atoms. The overlay
/// matcher (`permutation_for_ligands`) rejects a non-unique frame, so `apply` finds no identity
/// match — hence ligands must be unique. Real-atom ligands come first (distinct atoms); virtual
/// implicit-H / lone-pair fills pad by distinct `(atom, kind)` pairs. A frame of `degree` unique
/// ligands needs `atom_count * 3 >= degree`, so callers gate on `atom_count`.
fn unique_ligand_frame(
    atom_count: usize,
    degree: usize,
) -> impl Strategy<Value = Vec<StereoLigand>> {
    let pool: Vec<StereoLigand> = (0..atom_count as u32)
        .flat_map(|a| {
            [
                StereoLigandKind::Atom,
                StereoLigandKind::ImplicitHydrogen,
                StereoLigandKind::LonePair,
            ]
            .into_iter()
            .map(move |kind| StereoLigand::new(AtomId(a), kind))
        })
        .collect();
    Just(pool).prop_shuffle().prop_map(move |mut pool| {
        pool.truncate(degree);
        pool
    })
}

/// 0..=1 tetrahedral stereo atoms over an `atom_count`-atom molecule (needs `atom_count >= 2` for a
/// `degree`-length unique ligand frame). Site is any atom; ligands are distinct real/virtual ligands
/// whose atoms need not be graph neighbors (tier-1 requires only the ligand count for the kind).
fn stereo_atom_overlay_strategy(
    atom_count: usize,
) -> BoxedStrategy<Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>> {
    let degree = StereoKind::Tetrahedral.degree();
    if atom_count * 3 < degree {
        return Just(Vec::new()).boxed();
    }
    prop::collection::vec(
        (
            0..atom_count as u32,
            unique_ligand_frame(atom_count, degree),
            stereo_coset_for_kind(StereoKind::Tetrahedral),
        ),
        0..=1,
    )
    .prop_map(move |entries| {
        entries
            .into_iter()
            .map(|(site, ligands, coset)| {
                let ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
                (AtomId(site), ligands, ast)
            })
            .collect()
    })
    .boxed()
}

/// 0..=1 cis/trans stereo bonds (needs `atom_count >= 2` for a `degree`-length unique frame). Site is
/// any bond; ligands are distinct real/virtual ligands (their atoms need not be double-bond termini).
fn stereo_bond_overlay_strategy(
    atom_count: usize,
    bond_count: usize,
) -> BoxedStrategy<Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>> {
    let degree = StereoKind::CisTrans.degree();
    if bond_count == 0 || atom_count * 3 < degree {
        return Just(Vec::new()).boxed();
    }
    prop::collection::vec(
        (
            0..bond_count as u32,
            unique_ligand_frame(atom_count, degree),
            stereo_coset_for_kind(StereoKind::CisTrans),
        ),
        0..=1,
    )
    .prop_map(move |entries| {
        entries
            .into_iter()
            .map(|(site, ligands, coset)| {
                let ast = StereoBondAst::new(StereoKind::CisTrans, coset);
                (BondId(site), ligands, ast)
            })
            .collect()
    })
    .boxed()
}

/// A reaction whose `lhs` carries DAMN overlays — exercises overlay carry, correspondence, and
/// co-deletion through compose.
fn overlay_reaction_strategy() -> impl Strategy<Value = ReactionAst> {
    reaction_over(overlay_molecule_strategy())
}

/// An optional edit to a surviving stereo entity. The relative ops (`Swap` / `Mirror` / `Apply`)
/// resolve `old` from the matched host coset at apply, carrying no pre-state; `SetCoset` becomes a
/// `ModifyField { Configuration }` whose `old` is read from `lhs`, so apply's precondition holds.
#[derive(Clone, Debug)]
enum StereoOp {
    Swap,
    Mirror,
    Apply(Permutation),
    SetCoset(StereoCosetAst),
}

/// Per-surviving-stereo-entity optional op: `Swap` / `Mirror` use the kind's in-group generators,
/// `SetCoset` is bounded to the kind's in-range cosets, and `Apply` a permutation in the kind's
/// parent group. The coset algebra rejects out-of-group permutations (`reindex` → `None`, which
/// `act` unwraps), so `Apply` only draws arbitrary permutations for kinds whose parent is the full
/// symmetric group (Tetrahedral); other kinds omit it and lean on the in-group `Swap` / `Mirror`.
fn stereo_op_strategy(kind: StereoKind) -> impl Strategy<Value = Option<StereoOp>> {
    let base = prop_oneof![
        Just(StereoOp::Swap),
        Just(StereoOp::Mirror),
        stereo_coset_for_kind(kind).prop_map(StereoOp::SetCoset),
    ]
    .boxed();
    let ops = if kind == StereoKind::Tetrahedral {
        prop_oneof![
            base,
            permutation_strategy(kind.degree()).prop_map(StereoOp::Apply),
        ]
        .boxed()
    } else {
        base
    };
    prop::option::weighted(0.5, ops)
}

/// A valid reaction over any generated `lhs`: DPO-valid atom deletions (each removed atom takes its
/// incident bonds, overlays, and stereo entities), per-surviving-entity optional field / relative
/// edits (the absolute `old` read from `lhs`, so apply's precondition holds), plus up to two new
/// atoms bonded to the lowest survivor. No dangling by construction.
fn reaction_over(
    molecule: impl Strategy<Value = MoleculeAst>,
) -> impl Strategy<Value = ReactionAst> {
    molecule
        .prop_flat_map(|lhs| {
            let atom_count = lhs.atoms().count();
            let bond_count = lhs.bonds().count();
            let dative_count = lhs.dative_bonds().count();
            let aromatic_count = lhs.aromatic_systems().count();
            let multicenter_count = lhs.multicenter_bonds().count();
            let stereo_atom_count = lhs.stereo_atoms().count();
            let stereo_bond_count = lhs.stereo_bonds().count();
            (
                Just(lhs),
                prop::collection::vec(weighted(0.25), atom_count),
                prop::collection::vec(prop::option::of(-2i64..=2), atom_count),
                prop::collection::vec(prop::option::of(1i64..=3), bond_count),
                prop::collection::vec(element_strategy(), 0..=2),
                (
                    // Overlay `ModifyField` on survivors: dative order, aromatic / multicenter charge.
                    prop::collection::vec(prop::option::of(1i64..=3), dative_count),
                    prop::collection::vec(prop::option::of(-2i64..=2), aromatic_count),
                    prop::collection::vec(prop::option::of(-2i64..=2), multicenter_count),
                    // Add an `Aromatic` constraint to a surviving dative (guarded on absence).
                    prop::collection::vec(weighted(0.3), dative_count),
                    // Add a noncovalent overlay between the two newly-added atoms.
                    weighted(0.4),
                ),
                (
                    prop::collection::vec(
                        stereo_op_strategy(StereoKind::Tetrahedral),
                        stereo_atom_count,
                    ),
                    prop::collection::vec(
                        stereo_op_strategy(StereoKind::CisTrans),
                        stereo_bond_count,
                    ),
                ),
            )
        })
        .prop_map(
            |(lhs, removals, charges, orders, additions, overlay_ops, stereo_ops)| {
                build_reaction(
                    lhs,
                    removals,
                    charges,
                    orders,
                    additions,
                    overlay_ops,
                    stereo_ops,
                )
            },
        )
}

/// Per-entity overlay `ModifyField` / `Add` / `ModifyConstraint` randomness: dative orders, aromatic
/// charges, multicenter charges, dative-Aromatic-constraint flags, and the add-noncovalent flag.
type OverlayOps = (
    Vec<Option<i64>>,
    Vec<Option<i64>>,
    Vec<Option<i64>>,
    Vec<bool>,
    bool,
);

/// Per-stereo-entity optional op: stereo atoms, then stereo bonds.
type StereoOps = (Vec<Option<StereoOp>>, Vec<Option<StereoOp>>);

fn build_reaction(
    lhs: MoleculeAst,
    removals: Vec<bool>,
    charges: Vec<Option<i64>>,
    orders: Vec<Option<i64>>,
    additions: Vec<Element>,
    overlay_ops: OverlayOps,
    stereo_ops: StereoOps,
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
    // A removed atom also takes its incident stereo entities (site OR ligand incidence), else
    // apply / span / DpoValidator dangle. `incident_ids` covers both.
    let mut removed_stereo_atom: HashSet<StereoAtomId> = HashSet::new();
    let mut removed_stereo_bond: HashSet<StereoBondId> = HashSet::new();
    for &id in &removed_atoms {
        removed_stereo_atom.extend(lhs.stereo_atoms().incident_ids(id));
        removed_stereo_bond.extend(lhs.stereo_bonds().incident_ids(id));
    }
    for &id in &removed_stereo_atom {
        let view = lhs.stereo_atom(id);
        deltas.push(Delta::StereoAtom(StereoAtomDelta::Remove {
            id,
            site: view.site_id(),
            ligands: view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect(),
            ast: view.ast.clone(),
        }));
    }
    for &id in &removed_stereo_bond {
        let view = lhs.stereo_bond(id);
        deltas.push(Delta::StereoBond(StereoBondDelta::Remove {
            id,
            site: view.site_id(),
            ligands: view
                .ligands()
                .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                .collect(),
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
    // Part A — overlay `ModifyField` on survivors: read `old` from `lhs`, emit only when it changes.
    let (
        dative_orders,
        aromatic_charges,
        multicenter_charges,
        dative_aromatic_flags,
        add_noncovalent,
    ) = overlay_ops;
    for (index, new_order) in dative_orders.into_iter().enumerate() {
        let id = DativeBondId(index as u32);
        if removed_dative.contains(&id) {
            continue;
        }
        let Some(order) = new_order else { continue };
        let old = lhs.dative_bond(id).ast.order.clone();
        let new = ValueAst::Lit(order);
        if old != new {
            deltas.push(Delta::DativeBond(DativeBondDelta::ModifyField {
                id,
                change: DativeBondFieldChange::Order { old, new },
            }));
        }
    }
    for (index, new_charge) in aromatic_charges.into_iter().enumerate() {
        let id = AromaticSystemId(index as u32);
        if removed_aromatic.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.aromatic_system(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                id,
                change: AromaticSystemFieldChange::Charge { old, new },
            }));
        }
    }
    for (index, new_charge) in multicenter_charges.into_iter().enumerate() {
        let id = MulticenterBondId(index as u32);
        if removed_multicenter.contains(&id) {
            continue;
        }
        let Some(charge) = new_charge else { continue };
        let old = lhs.multicenter_bond(id).ast.charge.clone();
        let new = ValueAst::Lit(charge);
        if old != new {
            deltas.push(Delta::MulticenterBond(MulticenterBondDelta::ModifyField {
                id,
                change: MulticenterBondFieldChange::Charge { old, new },
            }));
        }
    }
    // Part A — add an `Aromatic` constraint to a surviving dative, guarded on its absence (apply's
    // `old: None` precondition requires no existing constraint under that key).
    for (index, add) in dative_aromatic_flags.into_iter().enumerate() {
        let id = DativeBondId(index as u32);
        if !add || removed_dative.contains(&id) {
            continue;
        }
        let has_aromatic = lhs
            .dative_bond(id)
            .ast
            .constraints
            .iter()
            .any(|c| matches!(c, DativeBondConstraintAst::Aromatic(_)));
        if has_aromatic {
            continue;
        }
        deltas.push(Delta::DativeBond(DativeBondDelta::ModifyConstraint {
            id,
            old: None,
            new: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        }));
    }
    // Part B — stereo edits on survivors. Relative ops resolve `old` from the host at apply;
    // `SetCoset` reads `old` from `lhs`. Every op is emitted only when it *changes* the entity's
    // configuration value: a value no-op (a relative op on an `Undetermined` coset, a `Mirror` on an
    // achiral kind, a stabilizer permutation, or a `SetCoset` to the current value) would materialize
    // a spurious `Modified { X, X }` span state that `to_reaction` diffs back to empty, breaking the
    // span roundtrip.
    let (stereo_atom_ops, stereo_bond_ops) = stereo_ops;
    for (index, op) in stereo_atom_ops.into_iter().enumerate() {
        let id = StereoAtomId(index as u32);
        if removed_stereo_atom.contains(&id) {
            continue;
        }
        let Some(op) = op else { continue };
        let kind = lhs.stereo_atom(id).kind();
        let old = lhs.stereo_atom(id).ast.configuration.clone();
        let (new, delta) = match &op {
            StereoOp::Swap => (old.swap(), StereoAtomDelta::Swap { id, kind }),
            StereoOp::Mirror => (old.mirror(), StereoAtomDelta::Mirror { id, kind }),
            StereoOp::Apply(permutation) => (
                old.apply(*permutation),
                StereoAtomDelta::Apply {
                    id,
                    kind,
                    permutation: *permutation,
                },
            ),
            StereoOp::SetCoset(coset) => {
                let new = StereoConfigurationAst::kinded(kind, coset.clone());
                (
                    new.clone(),
                    StereoAtomDelta::ModifyField {
                        id,
                        change: StereoAtomFieldChange::Configuration {
                            old: old.clone(),
                            new,
                        },
                    },
                )
            }
        };
        if new != old {
            deltas.push(Delta::StereoAtom(delta));
        }
    }
    for (index, op) in stereo_bond_ops.into_iter().enumerate() {
        let id = StereoBondId(index as u32);
        if removed_stereo_bond.contains(&id) {
            continue;
        }
        let Some(op) = op else { continue };
        let kind = lhs.stereo_bond(id).kind();
        let old = lhs.stereo_bond(id).ast.configuration.clone();
        let (new, delta) = match &op {
            StereoOp::Swap => (old.swap(), StereoBondDelta::Swap { id, kind }),
            StereoOp::Mirror => (old.mirror(), StereoBondDelta::Mirror { id, kind }),
            StereoOp::Apply(permutation) => (
                old.apply(*permutation),
                StereoBondDelta::Apply {
                    id,
                    kind,
                    permutation: *permutation,
                },
            ),
            StereoOp::SetCoset(coset) => {
                let new = StereoConfigurationAst::kinded(kind, coset.clone());
                (
                    new.clone(),
                    StereoBondDelta::ModifyField {
                        id,
                        change: StereoBondFieldChange::Configuration {
                            old: old.clone(),
                            new,
                        },
                    },
                )
            }
        };
        if new != old {
            deltas.push(Delta::StereoBond(delta));
        }
    }
    // Append atoms bonded to the lowest surviving atom (isolated if every atom is removed).
    let anchor = (0..atom_count as u32)
        .map(AtomId)
        .find(|id| !removed_atoms.contains(id));
    let mut added_atom_ids: Vec<AtomId> = Vec::new();
    for (offset, element) in additions.into_iter().enumerate() {
        let atom = AtomId((atom_count + offset) as u32);
        added_atom_ids.push(atom);
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
    // Part A — overlay `Add`: a noncovalent bond between the two newly-added atoms (both created in
    // this reaction, so no dangling). Ids append past the lhs noncovalent count.
    if add_noncovalent && added_atom_ids.len() >= 2 {
        deltas.push(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(lhs.noncovalent_bonds().count() as u32),
            atoms: [added_atom_ids[0], added_atom_ids[1]],
            ast: NoncovalentBondAst {
                kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::VanDerWaals),
                constraints: Default::default(),
            },
        }));
    }
    ReactionAst::new(lhs, Deltas::from_iter(deltas))
}

proptest! {
    /// A pattern-relative atom update lowers against the matched host atom, including independent
    /// spin leaves and keyed constraint set / replace / remove operations.
    #[test]
    fn test_reaction_ast_apply_atom_update(
        host_atom in atom_ast_strategy(),
        update in atom_update_strategy(),
    ) {
        let pattern_atom = AtomAst::default();
        let effective_update = pattern_atom.difference_to(&pattern_atom.update(&update));
        let expected_atom = host_atom.update(&effective_update).canonicalize().unwrap();
        let atom_deltas = AtomDelta::for_update(AtomId(0), &pattern_atom, &effective_update);
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::default()],
                ..Default::default()
            }),
            Deltas::from_iter(atom_deltas.into_iter().map(Delta::Atom)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![host_atom],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![expected_atom],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative localized-bond update lowers against the matched host bond.
    #[test]
    fn test_reaction_ast_apply_bond_update(
        host_bond in bond_ast_strategy(),
        update in bond_update_strategy(),
    ) {
        let pattern_bond = BondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let bond_deltas = BondDelta::for_update(BondId(0), &pattern_bond, &effective_update);
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(bond_deltas.into_iter().map(Delta::Bond)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative dative-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_dative_bond_update(
        host_bond in dative_bond_strategy(),
        update in dative_bond_update_strategy(),
    ) {
        let pattern_bond = DativeBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let dative_deltas = DativeBondDelta::for_update(
            DativeBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(dative_deltas.into_iter().map(Delta::DativeBond)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            dative: vec![(vec![AtomId(0)], AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            dative: vec![(vec![AtomId(0)], AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative aromatic-system update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_aromatic_system_update(
        mut host_system in aromatic_system_ast_for(3),
        update in aromatic_system_update_for(3),
    ) {
        host_system.spin = SpinStateAst::from((2_u8, 3_u8));
        let pattern_system = AromaticSystemAst::default();
        let effective_update = pattern_system.difference_to(&pattern_system.update(&update));
        let expected_system = host_system.update(&effective_update).canonicalize().unwrap();
        let aromatic_deltas = AromaticSystemDelta::for_update(
            AromaticSystemId(0),
            &pattern_system,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::O),
                ],
                aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], AromaticSystemAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(aromatic_deltas.into_iter().map(Delta::AromaticSystem)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], host_system)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], expected_system)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }
}

proptest! {
    /// A pattern-relative multicenter-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_multicenter_bond_update(
        mut host_bond in multicenter_bond_ast_for(3),
        update in multicenter_bond_update_for(3),
    ) {
        host_bond.spin = SpinStateAst::from((2_u8, 3_u8));
        let pattern_bond = MulticenterBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let multicenter_deltas = MulticenterBondDelta::for_update(
            MulticenterBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::O),
                ],
                multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(multicenter_deltas.into_iter().map(Delta::MulticenterBond)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative noncovalent-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_noncovalent_bond_update(
        host_bond in noncovalent_bond_ast_strategy(),
        update in noncovalent_bond_update_strategy(),
    ) {
        let pattern_bond = NoncovalentBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let noncovalent_deltas = NoncovalentBondDelta::for_update(
            NoncovalentBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(noncovalent_deltas.into_iter().map(Delta::NoncovalentBond)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative stereo-atom update lowers against the matched host configuration and
    /// keyed constraints.
    #[test]
    fn test_reaction_ast_apply_stereo_atom_update(
        host_coset in stereo_coset_for_kind(StereoKind::Tetrahedral),
        host_constraints in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
        update in stereo_atom_application_update_strategy(),
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
        let pattern_atom = StereoAtomAst::new(
            StereoKind::Tetrahedral,
            StereoCosetAst::Undetermined,
        );
        let host_atom = StereoAtomAst::new(StereoKind::Tetrahedral, host_coset)
            .with_constraints(host_constraints);
        let effective_update = pattern_atom.difference_to(&pattern_atom.update(&update));
        let expected_atom = host_atom.update(&effective_update).canonicalize().unwrap();
        let stereo_atom_deltas =
            StereoAtomDelta::for_update(StereoAtomId(0), &pattern_atom, &effective_update);
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::F),
        ];
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: atoms.clone(),
                stereo_atoms: vec![(AtomId(0), ligands.clone(), pattern_atom.clone())],
                ..Default::default()
            }),
            Deltas::from_iter(stereo_atom_deltas.into_iter().map(Delta::StereoAtom)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            stereo_atoms: vec![(AtomId(0), ligands.clone(), host_atom)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            stereo_atoms: vec![(AtomId(0), ligands, expected_atom)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative stereo-bond update lowers against the matched host configuration and
    /// keyed constraints.
    #[test]
    fn test_reaction_ast_apply_stereo_bond_update(
        host_coset in stereo_coset_for_kind(StereoKind::CisTrans),
        host_constraints in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        update in stereo_bond_application_update_strategy(),
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
        let pattern_bond = StereoBondAst::new(
            StereoKind::CisTrans,
            StereoCosetAst::Undetermined,
        );
        let host_bond = StereoBondAst::new(StereoKind::CisTrans, host_coset)
            .with_constraints(host_constraints);
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let stereo_bond_deltas =
            StereoBondDelta::for_update(StereoBondId(0), &pattern_bond, &effective_update);
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::F),
        ];
        let bonds = vec![(AtomId(0), AtomId(1), BondAst::from_order(2))];
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: atoms.clone(),
                bonds: bonds.clone(),
                stereo_bonds: vec![(BondId(0), ligands.clone(), pattern_bond.clone())],
                ..Default::default()
            }),
            Deltas::from_iter(stereo_bond_deltas.into_iter().map(Delta::StereoBond)),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), ligands.clone(), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(BondId(0), ligands, expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }
}

proptest! {
    /// Applying a reaction at the identity occurrence of its own `lhs` reproduces the span's
    /// `right()` — the `transact`-apply path agrees with the span projection.
    #[test]
    fn test_reaction_ast_apply_reproduces_right(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let right = span.rhs();
            prop_assert!(reaction
                .apply(&reaction.lhs, ALG)
                .unwrap()
                .any(|derivation| derivation.unwrap().rhs() == &right));
        }
    }

    /// Cross-validate the two span constructions: the direct `superimpose` (Strategy A) reproduces
    /// the span the delta path (`to_reaction_span`) builds. Recover `(L, R, C)` from the delta-path
    /// span and reassemble; a mismatch flags a diff-completeness or frame gap between the paths.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, span);
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
                prop_assert_eq!(reverse_span.lhs(), span.rhs());
                let forward_reactant = span.lhs();
                let reverse_product = reverse_span.rhs();
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
                let right = span.rhs();
                prop_assert!(composite
                    .apply(&composite.lhs, ALG)
                    .unwrap()
                    .any(|derivation| derivation.unwrap().rhs() == &right));
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
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
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
            let right = span.rhs();
            prop_assert!(reaction
                .apply(&reaction.lhs, ALG)
                .unwrap()
                .any(|derivation| derivation.unwrap().rhs() == &right));
        }
    }

    /// Cross-validate the two span constructions with overlays present: the direct `superimpose`
    /// reassembles the delta-path span across all overlay families, not just atoms/bonds.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path_overlay(
        reaction in overlay_reaction_strategy(),
    ) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, span);
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
                let right = span.rhs();
                prop_assert!(composite
                    .apply(&composite.lhs, ALG)
                    .unwrap()
                    .any(|derivation| derivation.unwrap().rhs() == &right));
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
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
        }

        for product in &composed {
            prop_assert!(sequential.contains(product));
        }
    }

    /// P1 completeness (`Full`): every sequential product (B applied after A) is also some
    /// composite's product. Together with `compose_sound_overlay` (composed ⊆ seq) this is set
    /// equality at `host = a.lhs`. Covers stereo: the reactants carry stereo overlays and the deltas
    /// stereo ops, glued and applied across ligand frames by `meet_pushout` / `apply_at`.
    #[test]
    fn test_reaction_ast_compose_complete_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
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
                DpoValidator
                    .validate_reaction(&composite.lhs, &composite.deltas)
                    .unwrap(),
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

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_atom_modification(
        old in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let new = 1 - old;
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let rule_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let old_ast = StereoAtomAst::new(StereoKind::Tetrahedral, old);
        let new_ast = StereoAtomAst::new(StereoKind::Tetrahedral, new);
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), rule_frame.clone(), old_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: old_ast.configuration,
                    new: new_ast.configuration,
                },
            })]),
        );
        let host_frame = permutation.act(&rule_frame);
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                host_frame.clone(),
                StereoAtomAst::new(StereoKind::Tetrahedral, old).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                host_frame,
                StereoAtomAst::new(StereoKind::Tetrahedral, new).apply(permutation),
            )],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_atom_removal(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let rule_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let rule_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), rule_frame.clone(), rule_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: rule_frame.clone(),
                ast: rule_ast,
            })]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                permutation.act(&rule_frame),
                StereoAtomAst::new(StereoKind::Tetrahedral, coset).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_bond_modification(
        old in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let new = 1 - old;
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let rule_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let old_ast = StereoBondAst::new(StereoKind::CisTrans, old);
        let new_ast = StereoBondAst::new(StereoKind::CisTrans, new);
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), rule_frame.clone(), old_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::ModifyField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: old_ast.configuration,
                    new: new_ast.configuration,
                },
            })]),
        );
        let host_frame = permutation.act(&rule_frame);
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(
                BondId(0),
                host_frame.clone(),
                StereoBondAst::new(StereoKind::CisTrans, old).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                host_frame,
                StereoBondAst::new(StereoKind::CisTrans, new).apply(permutation),
            )],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_bond_removal(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let rule_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let rule_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let lhs = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), rule_frame.clone(), rule_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                id: StereoBondId(0),
                site: BondId(0),
                ligands: rule_frame.clone(),
                ast: rule_ast,
            })]),
        );
        let host = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(
                BondId(0),
                permutation.act(&rule_frame),
                StereoBondAst::new(StereoKind::CisTrans, coset).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(&host, ALG)
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }
}
