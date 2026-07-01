//! Sequential reaction composition (A;B): the reaction whose application equals applying A then
//! B. Per overlap of A's product `R_A` with B's reactant `L_B`, the composite is built in one
//! id space and `canonicalize`d; overlaps with no `B.apply(A.apply(H))` witness (the DPO gluing
//! conditions) are rejected.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;

use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm, EdgeId, FactorOrdering, NodeId, Unordered,
};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::Constraints;
use super::dative::DativeBondAst;
use super::delta::{
    remap_delta, AromaticSystemDelta, AtomDelta, BondDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta,
};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::reaction::ReactionAst;
use super::remap::IdRemapping;
use super::traits::{Canonicalize, Lattice};

/// Which overlaps `compose` keeps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompositionScope {
    /// Keep only overlaps touching A's reaction center (sequential composites that actually
    /// chain A into B).
    RcAnchored,
    /// Keep every overlap, including the empty one (the free rule-algebra sum).
    Full,
}

impl ReactionAst {
    /// Sequential composites of `self` (A) then `other` (B): one per admissible overlap of A's
    /// product with B's reactant. `compose(A,B).apply(H)` equals `B.apply(A.apply(H))`.
    pub fn compose(&self, other: &ReactionAst, scope: CompositionScope) -> Vec<ReactionAst> {
        compose_all(self, other, scope).unwrap_or_default()
    }
}

fn created_atom_ids(deltas: &Deltas) -> Vec<AtomId> {
    let mut ids: Vec<AtomId> = deltas
        .iter()
        .filter_map(|d| match d {
            Delta::Atom(AtomDelta::Add { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    ids.sort();
    ids
}

fn created_bond_ids(deltas: &Deltas) -> Vec<BondId> {
    let mut ids: Vec<BondId> = deltas
        .iter()
        .filter_map(|d| match d {
            Delta::Bond(BondDelta::Add { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    ids.sort();
    ids
}

/// Assign composite overlay ids for one overlay kind over an overlap, in four ordered classes:
/// (1) lhs_A, (2) L_B context, (3) A-created, (4) B-created. Classes 1 and 2 live in `lhs_c`, in
/// push order. `a_lhs`/`lb` carry, per overlay, its composite correspondence key and its prebuilt
/// `lhs_c` entry; `lb` entries are in L_B id space, so `relabel` maps a context one into the
/// composite frame (re-aligning any positional payload). An L_B overlay whose participants all lie
/// in the overlap corresponds — by key — to an A-side overlay present in R_A (non-removed lhs_A +
/// A-created), reusing its id; a required overlap overlay with no match returns `None` (skip the
/// overlap). Keys are unique per overlay (spec §4.1). The key (`K`) and entry (`E`) are opaque, so a
/// bond-anchored kind (stereo bond) keys on a bond + atoms and the relabel routes its site through
/// the bond frame.
/// One overlay kind's composite placement: A-id→composite and B-id→composite maps, plus the kind's
/// `lhs_c` entries (classes 1 and 2, in id order).
type OverlayPlacement<I, E> = (HashMap<I, I>, HashMap<I, I>, Vec<E>);

fn place_overlays<I, K, E>(
    a_lhs: Vec<(I, K, E)>,
    a_removed: &HashSet<I>,
    a_created: &[(I, K)],
    lb: Vec<(I, K, E, bool)>,
    b_created: &[I],
    mut relabel: impl FnMut(E) -> E,
) -> Option<OverlayPlacement<I, E>>
where
    I: Copy + Eq + Hash + From<usize>,
    K: Eq + Hash + Clone,
{
    let mut da: HashMap<I, I> = HashMap::new();
    let mut db: HashMap<I, I> = HashMap::new();
    let mut lc: Vec<E> = Vec::new();
    let mut index: HashMap<K, I> = HashMap::new();

    // Class 1: lhs_A overlays, carried as-is; survivors seed the correspondence index.
    for (id, key, entry) in a_lhs {
        let cid = I::from(lc.len());
        da.insert(id, cid);
        if !a_removed.contains(&id) {
            index.insert(key, cid);
        }
        lc.push(entry);
    }

    // Class 2: L_B context overlays, relabeled into the composite frame; overlap-region deferred.
    let mut overlap_region: Vec<(I, K)> = Vec::new();
    for (id, key, entry, in_overlap) in lb {
        if in_overlap {
            overlap_region.push((id, key));
        } else {
            db.insert(id, I::from(lc.len()));
            lc.push(relabel(entry));
        }
    }
    let class12 = lc.len();

    // Class 3: A-created (also correspondents). Class 4: B-created.
    for (rank, (id, key)) in a_created.iter().enumerate() {
        let cid = I::from(class12 + rank);
        da.insert(*id, cid);
        index.insert(key.clone(), cid);
    }
    let class123 = class12 + a_created.len();
    for (rank, id) in b_created.iter().enumerate() {
        db.insert(*id, I::from(class123 + rank));
    }

    // Each overlap-region L_B overlay reuses its R_A correspondent's id.
    for (id, key) in overlap_region {
        db.insert(id, *index.get(&key)?);
    }
    Some((da, db, lc))
}

fn compose_all(
    a: &ReactionAst,
    b: &ReactionAst,
    scope: CompositionScope,
) -> Option<Vec<ReactionAst>> {
    let da = a.deltas.clone().canonicalize().ok()?;
    let db = b.deltas.clone().canonicalize().ok()?;
    // Stereo overlays arrive in I6; bail if either reactant carries one.
    if a.lhs.has_stereo_atoms()
        || a.lhs.has_stereo_bonds()
        || b.lhs.has_stereo_atoms()
        || b.lhs.has_stereo_bonds()
    {
        return None;
    }
    let span_a = a.to_reaction_span().ok()?;
    let r_a = span_a.right();
    let l_b = &b.lhs;
    let n_a = a.lhs.atoms().count();
    let m_a = a.lhs.bonds().count();
    let n_b = l_b.atoms().count();
    let m_b = l_b.bonds().count();

    // R_A id ⇒ A id: the span's union id space is `lhs_A` in place then A-created appended, and
    // `right()` keeps the survivors in that order, so the k-th survivor's A-id index is the
    // A id of R_A atom/bond k.
    let ra_atom_a_id: Vec<usize> = span_a
        .atoms()
        .iter()
        .enumerate()
        .filter(|(_, change)| change.right().is_some())
        .map(|(a_id, _)| a_id)
        .collect();
    let mut ra_bond_a_id: Vec<usize> = Vec::new();
    for (a_id, change) in span_a.bonds().iter().enumerate() {
        if change.right().is_none() {
            continue;
        }
        let [x, y] = span_a.graph().edge_endpoints(EdgeId(a_id as u32));
        if span_a.atoms()[x.index()].right().is_some()
            && span_a.atoms()[y.index()].right().is_some()
        {
            ra_bond_a_id.push(a_id);
        }
    }
    let a_created_atoms = span_a.atoms().len() - n_a;
    let a_created_bonds = span_a.bonds().len() - m_a;

    // A-created delta id ⇒ A-id rank (the span appends them sorted by id).
    let a_atom_rank: HashMap<AtomId, usize> = created_atom_ids(&da)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let a_bond_rank: HashMap<BondId, usize> = created_bond_ids(&da)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();

    // Reaction center, projected onto R_A: A-created elements (all) plus `lhs_A` elements A
    // modifies. A id of a changed element → which R_A atoms/bonds it is.
    let mut rc_a_atom_ids: HashSet<usize> = HashSet::new();
    let mut rc_a_bond_ids: HashSet<usize> = HashSet::new();
    for delta in da.iter() {
        match delta {
            Delta::Atom(AtomDelta::ModifyField { id, .. })
            | Delta::Atom(AtomDelta::ModifyConstraint { id, .. }) => {
                rc_a_atom_ids.insert(id.index());
            }
            Delta::Atom(AtomDelta::Add { id, .. }) => {
                rc_a_atom_ids.insert(n_a + a_atom_rank[id]);
            }
            Delta::Bond(BondDelta::ModifyField { id, .. })
            | Delta::Bond(BondDelta::ModifyConstraint { id, .. }) => {
                rc_a_bond_ids.insert(id.index());
            }
            Delta::Bond(BondDelta::Add { id, .. }) => {
                rc_a_bond_ids.insert(m_a + a_bond_rank[id]);
            }
            _ => {}
        }
    }
    // Any overlay delta touches its participant atoms, which join A's reaction center (so an
    // overlay-only A edit anchors an overlap). `Add`/`Remove` carry their atoms; the `Modify`
    // variants name only the overlay, so look its atoms up in `lhs_A`.
    for delta in da.iter() {
        let atoms: Vec<AtomId> = match delta {
            Delta::DativeBond(
                DativeBondDelta::Add {
                    donors, acceptor, ..
                }
                | DativeBondDelta::Remove {
                    donors, acceptor, ..
                },
            ) => donors.iter().copied().chain([*acceptor]).collect(),
            Delta::DativeBond(
                DativeBondDelta::ModifyField { id, .. }
                | DativeBondDelta::ModifyConstraint { id, .. },
            ) => a
                .lhs
                .dative_bonds()
                .get(*id)
                .map_or_else(Vec::new, |x| x.atom_ids().collect()),
            Delta::AromaticSystem(
                AromaticSystemDelta::Add { atoms, .. } | AromaticSystemDelta::Remove { atoms, .. },
            ) => atoms.clone(),
            Delta::AromaticSystem(
                AromaticSystemDelta::ModifyField { id, .. }
                | AromaticSystemDelta::ModifyConstraint { id, .. },
            ) => a
                .lhs
                .aromatic_systems()
                .get(*id)
                .map_or_else(Vec::new, |x| x.atom_ids().collect()),
            Delta::MulticenterBond(
                MulticenterBondDelta::Add { atoms, .. }
                | MulticenterBondDelta::Remove { atoms, .. },
            ) => atoms.clone(),
            Delta::MulticenterBond(
                MulticenterBondDelta::ModifyField { id, .. }
                | MulticenterBondDelta::ModifyConstraint { id, .. },
            ) => a
                .lhs
                .multicenter_bonds()
                .get(*id)
                .map_or_else(Vec::new, |x| x.atom_ids().collect()),
            Delta::NoncovalentBond(
                NoncovalentBondDelta::Add { atoms, .. }
                | NoncovalentBondDelta::Remove { atoms, .. },
            ) => atoms.to_vec(),
            Delta::NoncovalentBond(
                NoncovalentBondDelta::ModifyField { id, .. }
                | NoncovalentBondDelta::ModifyConstraint { id, .. },
            ) => a
                .lhs
                .noncovalent_bonds()
                .get(*id)
                .map_or_else(Vec::new, |x| x.atom_ids().to_vec()),
            _ => continue,
        };
        for atom in atoms {
            let index = a_atom_rank
                .get(&atom)
                .map_or(atom.index(), |&rank| n_a + rank);
            rc_a_atom_ids.insert(index);
        }
    }
    let rc_ra_atoms: HashSet<AtomId> = (0..ra_atom_a_id.len() as u32)
        .map(AtomId)
        .filter(|k| rc_a_atom_ids.contains(&ra_atom_a_id[k.index()]))
        .collect();
    let rc_ra_bonds: HashSet<BondId> = (0..ra_bond_a_id.len() as u32)
        .map(BondId)
        .filter(|k| rc_a_bond_ids.contains(&ra_bond_a_id[k.index()]))
        .collect();

    let mut node_match = |ra: NodeId, lb: NodeId| {
        r_a.atom(AtomId::from(ra))
            .ast
            .meet(l_b.atom(AtomId::from(lb)).ast)
            .is_some()
    };
    let mut edge_match = |re: EdgeId, le: EdgeId| {
        r_a.bond(BondId::from(re))
            .ast
            .meet(l_b.bond(BondId::from(le)).ast)
            .is_some()
    };
    let overlaps = r_a.raw_graph().enumerate_common_subgraphs(
        l_b.raw_graph(),
        &mut node_match,
        &mut edge_match,
        CommonSubgraphEnumerationAlgorithm::BronKerbosch,
    );

    let db_created_atom_rank: HashMap<AtomId, usize> = created_atom_ids(&db)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let db_created_bond_rank: HashMap<BondId, usize> = created_bond_ids(&db)
        .into_iter()
        .enumerate()
        .map(|(rank, id)| (id, rank))
        .collect();
    let db_removed_atoms: HashSet<AtomId> = db
        .iter()
        .filter_map(|d| match d {
            Delta::Atom(AtomDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    let db_removed_bonds: HashSet<BondId> = db
        .iter()
        .filter_map(|d| match d {
            Delta::Bond(BondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        })
        .collect();
    // Overlays B deletes, keyed by sorted L_B participant set (for the combined-frame overlay
    // dangling check). A deleted shared atom's R_A overlay must correspond to one of these.
    let sorted = |mut p: Vec<AtomId>| {
        p.sort();
        p
    };
    let b_removed_dative: HashSet<Vec<AtomId>> = db
        .iter()
        .filter_map(|d| match d {
            Delta::DativeBond(DativeBondDelta::Remove {
                donors, acceptor, ..
            }) => Some(sorted(donors.iter().copied().chain([*acceptor]).collect())),
            _ => None,
        })
        .collect();
    let b_removed_aromatic: HashSet<Vec<AtomId>> = db
        .iter()
        .filter_map(|d| match d {
            Delta::AromaticSystem(AromaticSystemDelta::Remove { atoms, .. }) => {
                Some(sorted(atoms.clone()))
            }
            _ => None,
        })
        .collect();
    let b_removed_multicenter: HashSet<Vec<AtomId>> = db
        .iter()
        .filter_map(|d| match d {
            Delta::MulticenterBond(MulticenterBondDelta::Remove { atoms, .. }) => {
                Some(sorted(atoms.clone()))
            }
            _ => None,
        })
        .collect();
    let b_removed_noncovalent: HashSet<Vec<AtomId>> = db
        .iter()
        .filter_map(|d| match d {
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { atoms, .. }) => {
                Some(sorted(atoms.to_vec()))
            }
            _ => None,
        })
        .collect();

    let mut results = Vec::new();
    for overlap in overlaps {
        let mapping = overlap.mapping();
        let lb_to_ra: HashMap<AtomId, AtomId> = mapping
            .iter()
            .map(|&(ra, lb)| (AtomId::from(lb), AtomId::from(ra)))
            .collect();
        let ra_to_lb: HashMap<AtomId, AtomId> = mapping
            .iter()
            .map(|&(ra, lb)| (AtomId::from(ra), AtomId::from(lb)))
            .collect();
        let overlap_lb: HashSet<AtomId> = lb_to_ra.keys().copied().collect();
        let overlap_ra: HashSet<AtomId> = ra_to_lb.keys().copied().collect();

        if matches!(scope, CompositionScope::RcAnchored) {
            let touches_atom = overlap_ra.iter().any(|ra| rc_ra_atoms.contains(ra));
            let touches_bond = rc_ra_bonds.iter().any(|&rb| {
                let [x, y] = r_a.raw_graph().edge_endpoints(EdgeId(rb.0));
                overlap_ra.contains(&AtomId::from(x)) && overlap_ra.contains(&AtomId::from(y))
            });
            if !touches_atom && !touches_bond {
                continue;
            }
        }

        let is_ra_created = |ra: AtomId| ra_atom_a_id[ra.index()] >= n_a;

        // Extra (non-overlap) L_B atoms → composite class 2 (`n_a..n_a+e`).
        let lb_extra: Vec<AtomId> = (0..n_b as u32)
            .map(AtomId)
            .filter(|x| !overlap_lb.contains(x))
            .collect();
        let e = lb_extra.len();
        let lb_extra_comp: HashMap<AtomId, usize> = lb_extra
            .iter()
            .enumerate()
            .map(|(rank, &x)| (x, n_a + rank))
            .collect();

        // Classify L_B bonds; context bonds (incident to an extra atom) become composite class 2.
        // Pushout-complement: a context bond whose overlap endpoint is A-created has no place in
        // `lhs_C` → the overlap is inadmissible.
        let mut context_bonds: Vec<BondId> = Vec::new();
        let mut admissible = true;
        for j in 0..m_b as u32 {
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(j));
            let (u, v) = (AtomId::from(u), AtomId::from(v));
            let u_overlap = overlap_lb.contains(&u);
            let v_overlap = overlap_lb.contains(&v);
            if u_overlap && v_overlap {
                continue;
            }
            if (u_overlap && is_ra_created(lb_to_ra[&u]))
                || (v_overlap && is_ra_created(lb_to_ra[&v]))
            {
                admissible = false;
                break;
            }
            context_bonds.push(BondId(j));
        }
        if !admissible {
            continue;
        }
        let f = context_bonds.len();
        let context_comp: HashMap<BondId, usize> = context_bonds
            .iter()
            .enumerate()
            .map(|(rank, &x)| (x, m_a + rank))
            .collect();

        // Overlay pushout-complement: a context overlay (≥1 participant off the overlap) whose
        // overlap participant is A-created has no place in `lhs_c` (created atoms are not in `lhs_c`)
        // → the overlap is inadmissible, mirroring the context-bond check above.
        let overlay_context_inadmissible = |atoms: &[AtomId]| {
            atoms.iter().any(|a| !overlap_lb.contains(a))
                && atoms
                    .iter()
                    .any(|a| overlap_lb.contains(a) && is_ra_created(lb_to_ra[a]))
        };
        let inadmissible_overlay = l_b
            .dative_bonds()
            .iter()
            .any(|x| overlay_context_inadmissible(&x.atom_ids().collect::<Vec<_>>()))
            || l_b
                .aromatic_systems()
                .iter()
                .any(|x| overlay_context_inadmissible(&x.atom_ids().collect::<Vec<_>>()))
            || l_b
                .multicenter_bonds()
                .iter()
                .any(|x| overlay_context_inadmissible(&x.atom_ids().collect::<Vec<_>>()))
            || l_b
                .noncovalent_bonds()
                .iter()
                .any(|x| overlay_context_inadmissible(&x.atom_ids()));
        if inadmissible_overlay {
            continue;
        }

        // Composite-id-space dangling: if B deletes a shared (overlap) atom, every R_A bond or
        // overlay incident to its image must be one B also deletes; an A-product bond/overlay B
        // cannot see would dangle.
        let mut dangling = false;
        // An R_A overlay on a deleted shared atom dangles unless all its participants are in the
        // overlap and B removes the corresponding overlay (matched by sorted L_B participant set).
        let overlay_dangles = |parts: Vec<AtomId>, removed: &HashSet<Vec<AtomId>>| -> bool {
            let mut lb: Vec<AtomId> = Vec::with_capacity(parts.len());
            for p in parts {
                match ra_to_lb.get(&p) {
                    Some(&id) => lb.push(id),
                    None => return true,
                }
            }
            lb.sort();
            !removed.contains(&lb)
        };
        for &u in &db_removed_atoms {
            if !overlap_lb.contains(&u) {
                continue;
            }
            let ru = lb_to_ra[&u];
            for rb in r_a.atom(ru).bond_ids() {
                let [x, y] = r_a.raw_graph().edge_endpoints(EdgeId(rb.0));
                let other = if AtomId::from(x) == ru {
                    AtomId::from(y)
                } else {
                    AtomId::from(x)
                };
                if !overlap_ra.contains(&other) {
                    dangling = true;
                    break;
                }
                let w = ra_to_lb[&other];
                match l_b.raw_graph().find_edge(NodeId::from(u), NodeId::from(w)) {
                    Some(le) if db_removed_bonds.contains(&BondId::from(le)) => {}
                    _ => {
                        dangling = true;
                        break;
                    }
                }
            }
            if !dangling {
                for od in r_a.atom(ru).dative_bond_ids() {
                    if overlay_dangles(r_a.dative_bond(od).atom_ids().collect(), &b_removed_dative) {
                        dangling = true;
                        break;
                    }
                }
            }
            if !dangling {
                if let Some(oa) = r_a.atom(ru).aromatic_system_id() {
                    if overlay_dangles(
                        r_a.aromatic_system(oa).atom_ids().collect(),
                        &b_removed_aromatic,
                    ) {
                        dangling = true;
                    }
                }
            }
            if !dangling {
                for om in r_a.atom(ru).multicenter_bond_ids() {
                    if overlay_dangles(
                        r_a.multicenter_bond(om).atom_ids().collect(),
                        &b_removed_multicenter,
                    ) {
                        dangling = true;
                        break;
                    }
                }
            }
            if !dangling {
                for on in r_a.atom(ru).noncovalent_bond_ids() {
                    if overlay_dangles(
                        r_a.noncovalent_bond(on).atom_ids().to_vec(),
                        &b_removed_noncovalent,
                    ) {
                        dangling = true;
                        break;
                    }
                }
            }
            if dangling {
                break;
            }
        }
        if dangling {
            continue;
        }

        // Composite id maps. A id → composite: `lhs_A`/`lhs_A`-bonds keep their id; A-created
        // shift past the appended extras (atoms by `e`, bonds by `f`).
        let composite_atom_id = |a_id: usize| {
            if a_id < n_a {
                a_id
            } else {
                n_a + e + (a_id - n_a)
            }
        };
        let composite_bond_id = |a_id: usize| {
            if a_id < m_a {
                a_id
            } else {
                m_a + f + (a_id - m_a)
            }
        };
        let ra_atom_comp = |ra: AtomId| AtomId(composite_atom_id(ra_atom_a_id[ra.index()]) as u32);
        let ra_bond_comp = |rb: BondId| BondId(composite_bond_id(ra_bond_a_id[rb.index()]) as u32);

        let mut da_atom: HashMap<AtomId, AtomId> =
            (0..n_a as u32).map(|i| (AtomId(i), AtomId(i))).collect();
        for (&id, &rank) in &a_atom_rank {
            da_atom.insert(id, AtomId((n_a + e + rank) as u32));
        }
        let mut da_bond: HashMap<BondId, BondId> =
            (0..m_a as u32).map(|j| (BondId(j), BondId(j))).collect();
        for (&id, &rank) in &a_bond_rank {
            da_bond.insert(id, BondId((m_a + f + rank) as u32));
        }

        let b_atom_base = n_a + e + a_created_atoms;
        let b_bond_base = m_a + f + a_created_bonds;
        let mut db_atom: HashMap<AtomId, AtomId> = HashMap::new();
        for j in 0..n_b as u32 {
            let lb = AtomId(j);
            let comp = match lb_to_ra.get(&lb) {
                Some(&ra) => ra_atom_comp(ra),
                None => AtomId(lb_extra_comp[&lb] as u32),
            };
            db_atom.insert(lb, comp);
        }
        for (&id, &rank) in &db_created_atom_rank {
            db_atom.insert(id, AtomId((b_atom_base + rank) as u32));
        }
        let mut db_bond: HashMap<BondId, BondId> = HashMap::new();
        for j in 0..m_b as u32 {
            let lb = BondId(j);
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(j));
            let (u, v) = (AtomId::from(u), AtomId::from(v));
            let comp = if overlap_lb.contains(&u) && overlap_lb.contains(&v) {
                let re = r_a
                    .raw_graph()
                    .find_edge(NodeId::from(lb_to_ra[&u]), NodeId::from(lb_to_ra[&v]))
                    .expect("an induced overlap bond exists in R_A");
                ra_bond_comp(BondId::from(re))
            } else {
                BondId(context_comp[&lb] as u32)
            };
            db_bond.insert(lb, comp);
        }
        for (&id, &rank) in &db_created_bond_rank {
            db_bond.insert(id, BondId((b_bond_base + rank) as u32));
        }

        // lhs_C = lhs_A + extra L_B atoms, with lhs_A bonds + the L_B context bonds.
        let mut lc_atoms: Vec<AtomAst> = Vec::with_capacity(n_a + e);
        for i in 0..n_a as u32 {
            lc_atoms.push(a.lhs.atom(AtomId(i)).ast.clone());
        }
        for &x in &lb_extra {
            lc_atoms.push(l_b.atom(x).ast.clone());
        }
        let mut lc_bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::new();
        for j in 0..m_a as u32 {
            let [x, y] = a.lhs.raw_graph().edge_endpoints(EdgeId(j));
            lc_bonds.push((
                AtomId::from(x),
                AtomId::from(y),
                a.lhs.bond(BondId(j)).ast.clone(),
            ));
        }
        for &cb in &context_bonds {
            let [u, v] = l_b.raw_graph().edge_endpoints(EdgeId(cb.0));
            lc_bonds.push((
                db_atom[&AtomId::from(u)],
                db_atom[&AtomId::from(v)],
                l_b.bond(cb).ast.clone(),
            ));
        }

        // Overlay frame, DAMN order (stereo is I6): per kind, place the four-class ids and carry
        // class 1 (lhs_A, identity) + class 2 (L_B context, relabeled) into `lhs_c`; classes 3/4 are
        // delta-created. Aromatic/multicenter relabel re-sorts participants and permutes electrons.

        // Dative.
        let removed_d: HashSet<DativeBondId> = da
            .iter()
            .filter_map(|d| match d {
                Delta::DativeBond(DativeBondDelta::Remove { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        let a_lhs_d = a
            .lhs
            .dative_bonds()
            .iter()
            .map(|x| {
                let acceptor = x.acceptor_id();
                let donors: Vec<AtomId> = x.atom_ids().filter(|a| *a != acceptor).collect();
                let mut key: Vec<AtomId> = x.atom_ids().collect();
                key.sort();
                (x.id, key, (donors, acceptor, x.ast.clone()))
            })
            .collect();
        let mut a_created_d: Vec<(DativeBondId, Vec<AtomId>)> = da
            .iter()
            .filter_map(|d| match d {
                Delta::DativeBond(DativeBondDelta::Add {
                    id,
                    donors,
                    acceptor,
                    ..
                }) => {
                    let mut key: Vec<AtomId> = donors.iter().map(|a| da_atom[a]).collect();
                    key.push(da_atom[acceptor]);
                    key.sort();
                    Some((*id, key))
                }
                _ => None,
            })
            .collect();
        a_created_d.sort_by_key(|(id, _)| *id);
        let lb_d = l_b
            .dative_bonds()
            .iter()
            .map(|x| {
                let acceptor = x.acceptor_id();
                let atoms: Vec<AtomId> = x.atom_ids().collect();
                let in_overlap = atoms.iter().all(|a| overlap_lb.contains(a));
                let donors: Vec<AtomId> = atoms.iter().copied().filter(|a| *a != acceptor).collect();
                let mut key: Vec<AtomId> = atoms.iter().map(|a| db_atom[a]).collect();
                key.sort();
                (x.id, key, (donors, acceptor, x.ast.clone()), in_overlap)
            })
            .collect();
        let mut b_created_d: Vec<DativeBondId> = db
            .iter()
            .filter_map(|d| match d {
                Delta::DativeBond(DativeBondDelta::Add { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        b_created_d.sort();
        let Some((da_dative, db_dative, lc_dative)) = place_overlays(
            a_lhs_d,
            &removed_d,
            &a_created_d,
            lb_d,
            &b_created_d,
            |(donors, acceptor, ast): (Vec<AtomId>, AtomId, DativeBondAst)| {
                (
                    donors.iter().map(|a| db_atom[a]).collect(),
                    db_atom[&acceptor],
                    ast,
                )
            },
        ) else {
            continue;
        };

        // Aromatic.
        let removed_a: HashSet<AromaticSystemId> = da
            .iter()
            .filter_map(|d| match d {
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        let a_lhs_a = a
            .lhs
            .aromatic_systems()
            .iter()
            .map(|x| {
                let mut key: Vec<AtomId> = x.atom_ids().collect();
                key.sort();
                (x.id, key, (x.atom_ids().collect(), x.ast.clone()))
            })
            .collect();
        let mut a_created_a: Vec<(AromaticSystemId, Vec<AtomId>)> = da
            .iter()
            .filter_map(|d| match d {
                Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. }) => {
                    let mut key: Vec<AtomId> = atoms.iter().map(|a| da_atom[a]).collect();
                    key.sort();
                    Some((*id, key))
                }
                _ => None,
            })
            .collect();
        a_created_a.sort_by_key(|(id, _)| *id);
        let lb_a = l_b
            .aromatic_systems()
            .iter()
            .map(|x| {
                let atoms: Vec<AtomId> = x.atom_ids().collect();
                let in_overlap = atoms.iter().all(|a| overlap_lb.contains(a));
                let mut key: Vec<AtomId> = atoms.iter().map(|a| db_atom[a]).collect();
                key.sort();
                (x.id, key, (atoms, x.ast.clone()), in_overlap)
            })
            .collect();
        let mut b_created_a: Vec<AromaticSystemId> = db
            .iter()
            .filter_map(|d| match d {
                Delta::AromaticSystem(AromaticSystemDelta::Add { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        b_created_a.sort();
        let Some((da_aromatic, db_aromatic, lc_aromatic)) = place_overlays(
            a_lhs_a,
            &removed_a,
            &a_created_a,
            lb_a,
            &b_created_a,
            |(members, ast): (Vec<AtomId>, AromaticSystemAst)| {
                let mut members: Vec<AtomId> = members.iter().map(|a| db_atom[a]).collect();
                let order = Unordered::canonicalize_positions(&mut members);
                let mut ast = ast;
                ast.permute(&order);
                (members, ast)
            },
        ) else {
            continue;
        };

        // Multicenter.
        let removed_m: HashSet<MulticenterBondId> = da
            .iter()
            .filter_map(|d| match d {
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        let a_lhs_m = a
            .lhs
            .multicenter_bonds()
            .iter()
            .map(|x| {
                let mut key: Vec<AtomId> = x.atom_ids().collect();
                key.sort();
                (x.id, key, (x.atom_ids().collect(), x.ast.clone()))
            })
            .collect();
        let mut a_created_m: Vec<(MulticenterBondId, Vec<AtomId>)> = da
            .iter()
            .filter_map(|d| match d {
                Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. }) => {
                    let mut key: Vec<AtomId> = atoms.iter().map(|a| da_atom[a]).collect();
                    key.sort();
                    Some((*id, key))
                }
                _ => None,
            })
            .collect();
        a_created_m.sort_by_key(|(id, _)| *id);
        let lb_m = l_b
            .multicenter_bonds()
                .iter()
                .map(|x| {
                    let atoms: Vec<AtomId> = x.atom_ids().collect();
                    let in_overlap = atoms.iter().all(|a| overlap_lb.contains(a));
                    let mut key: Vec<AtomId> = atoms.iter().map(|a| db_atom[a]).collect();
                    key.sort();
                    (x.id, key, (atoms, x.ast.clone()), in_overlap)
                })
                .collect();
        let mut b_created_m: Vec<MulticenterBondId> = db
            .iter()
            .filter_map(|d| match d {
                Delta::MulticenterBond(MulticenterBondDelta::Add { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        b_created_m.sort();
        let Some((da_multicenter, db_multicenter, lc_multicenter)) = place_overlays(
            a_lhs_m,
            &removed_m,
            &a_created_m,
            lb_m,
            &b_created_m,
            |(members, ast): (Vec<AtomId>, MulticenterBondAst)| {
                let mut members: Vec<AtomId> = members.iter().map(|a| db_atom[a]).collect();
                let order = Unordered::canonicalize_positions(&mut members);
                let mut ast = ast;
                ast.permute(&order);
                (members, ast)
            },
        ) else {
            continue;
        };

        // Noncovalent.
        let removed_n: HashSet<NoncovalentBondId> = da
            .iter()
            .filter_map(|d| match d {
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        let a_lhs_n = a
            .lhs
            .noncovalent_bonds()
            .iter()
            .map(|x| {
                let [u, v] = x.atom_ids();
                let mut key = vec![u, v];
                key.sort();
                (x.id, key, (u, v, x.ast.clone()))
            })
            .collect();
        let mut a_created_n: Vec<(NoncovalentBondId, Vec<AtomId>)> = da
            .iter()
            .filter_map(|d| match d {
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, .. }) => {
                    let mut key: Vec<AtomId> = atoms.iter().map(|a| da_atom[a]).collect();
                    key.sort();
                    Some((*id, key))
                }
                _ => None,
            })
            .collect();
        a_created_n.sort_by_key(|(id, _)| *id);
        let lb_n = l_b
            .noncovalent_bonds()
                .iter()
                .map(|x| {
                    let [u, v] = x.atom_ids();
                    let in_overlap = overlap_lb.contains(&u) && overlap_lb.contains(&v);
                    let mut key = vec![db_atom[&u], db_atom[&v]];
                    key.sort();
                    (x.id, key, (u, v, x.ast.clone()), in_overlap)
                })
                .collect();
        let mut b_created_n: Vec<NoncovalentBondId> = db
            .iter()
            .filter_map(|d| match d {
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, .. }) => Some(*id),
                _ => None,
            })
            .collect();
        b_created_n.sort();
        let Some((da_noncovalent, db_noncovalent, lc_noncovalent)) = place_overlays(
            a_lhs_n,
            &removed_n,
            &a_created_n,
            lb_n,
            &b_created_n,
            |(u, v, ast): (AtomId, AtomId, NoncovalentBondAst)| (db_atom[&u], db_atom[&v], ast),
        ) else {
            continue;
        };

        let lhs_c = MoleculeAst::from_parts(
            lc_atoms,
            lc_bonds,
            lc_dative,
            lc_aromatic,
            lc_multicenter,
            lc_noncovalent,
            Vec::new(),
            Vec::new(),
            Constraints::new(),
        );

        // Stereo overlays have no deltas yet (I6).
        let da_map = IdRemapping::new(
            da_atom,
            da_bond,
            da_dative,
            da_aromatic,
            da_multicenter,
            da_noncovalent,
            HashMap::new(),
            HashMap::new(),
        );
        let db_map = IdRemapping::new(
            db_atom,
            db_bond,
            db_dative,
            db_aromatic,
            db_multicenter,
            db_noncovalent,
            HashMap::new(),
            HashMap::new(),
        );
        let mut deltas: Vec<Delta> = Vec::with_capacity(da.len() + db.len());
        for delta in da.iter() {
            deltas.push(remap_delta(delta.clone(), &da_map));
        }
        for delta in db.iter() {
            deltas.push(remap_delta(delta.clone(), &db_map));
        }
        let Ok(deltas) = Deltas::from_iter(deltas).canonicalize() else {
            continue;
        };
        results.push(ReactionAst::new(lhs_c, deltas));
    }
    Some(results)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::SubgraphIsomorphismAlgorithm;

    use super::super::constraint::Constraints;
    use super::super::edit::{BondFieldChange, NoncovalentBondFieldChange};
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
    use super::super::value::ValueAst;
    use super::*;

    // C-O order 1→2 then 2→3; the single overlap fuses to 1→3.
    #[rstest]
    #[case::fuse(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::Full,
        vec![ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A appends an O bonded to C (O is A-created); B raises that C-O 1→2. The composite appends the
    // O already at order 2 (create-then-modify fuses across the seam).
    #[case::created_atom(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::Full,
        vec![ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(2),
                }),
            ]),
        )]
    )]
    // A appends a C to N (R_A = N-C); B's reactant N-C-O maps the A-created C onto the middle atom,
    // whose bond to the extra O is a boundary bond on an A-created atom — unrealizable, rejected.
    #[case::inadmissible(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::N)], vec![]),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(1), AtomId(2), BondAst::from_order(1)),
                ],
            ),
            Deltas::new(),
        ),
        CompositionScope::Full,
        vec![]
    )]
    // A raises C-N 1→2 and adds a hydrogen bond across the pair (a created overlay); B raises 2→3.
    // The RC-anchored overlap fuses the bond to 1→3 and carries the noncovalent bond at id 0.
    #[case::overlay(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's lhs carries a hydrogen bond it never touches (only a covalent-order edit); B raises 2→3.
    // The composite carries the noncovalent bond (class ①) and fuses the order to 1→3.
    #[case::carry(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A removes its carried hydrogen bond. The composite carries it (class ①) and re-anchors A's
    // Remove delta onto composite noncovalent id 0.
    #[case::remove_carried(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // Both A's product and B's reactant carry the hydrogen bond on the overlap; B retypes it. The
    // overlap-region overlay corresponds (no fresh id), so B's modify re-anchors onto A's bond.
    #[case::correspondence(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // B's reactant requires a hydrogen bond on the overlap that A's product does not supply — the
    // overlap has no overlay correspondent, so it is skipped and compose yields nothing.
    #[case::required_absent(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
                vec![], vec![], vec![],
                vec![(
                    AtomId(0),
                    AtomId(1),
                    NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![]
    )]
    // A carries an aromatic system (a positional family) it never touches; the composite carries it
    // (class ①, identity participants) and fuses the order.
    #[case::aromatic_carry(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![],
                vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![],
                vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(3) },
            })]),
        )]
    )]
    // A's only edit modifies its own hydrogen bond (no atom/bond delta); B raises the order. The
    // overlay modify extends A's reaction center (via the lhs lookup), so the overlap is RC-anchored.
    #[case::rc_modify(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
                    id: NoncovalentBondId(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                    },
                }),
            ]),
        )]
    )]
    // A's only edit removes its own hydrogen bond (no atom/bond delta). The overlay remove (which
    // carries its atoms) anchors the overlap; the composite carries the bond and re-anchors the remove.
    #[case::rc_remove(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's only edit *creates* a hydrogen bond (no atom/bond delta). The overlay add extends A's
    // reaction center too, so the overlap is RC-anchored; the composite creates the bond at id 0.
    #[case::rc_add(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
                }),
            ]),
        )]
    )]
    // A's only edit removes its aromatic system (a positional kind, no atom/bond delta) — exercises
    // the aromatic RC arm and carry. The overlay remove anchors the overlap.
    #[case::rc_aromatic_remove(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![],
                vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                id: AromaticSystemId(0),
                atoms: vec![AtomId(0), AtomId(1)],
                ast: AromaticSystemAst::from_electrons(vec![1, 2]),
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                vec![],
                vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::AromaticSystem(AromaticSystemDelta::Remove {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                    ast: AromaticSystemAst::from_electrons(vec![1, 2]),
                }),
            ]),
        )]
    )]
    // Disjoint reactants (C-C, N-N — no matchable atom): the only overlap is empty, so `Full` is the
    // disjoint sum A ⊔ B — ids concatenated, both bond modifies relabeled (B's bond 0 → 1).
    #[case::disjoint_sum(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::Full,
        vec![ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::N),
                ],
                vec![
                    (AtomId(0), AtomId(1), BondAst::from_order(1)),
                    (AtomId(2), AtomId(3), BondAst::from_order(1)),
                ],
            ),
            Deltas::from_iter([
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
                Delta::Bond(BondDelta::ModifyField {
                    id: BondId(1),
                    change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
                }),
            ]),
        )]
    )]
    // The empty overlap misses A's reaction center, so `RcAnchored` drops it — no composite.
    #[case::disjoint_rc_anchored(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::N)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        ),
        CompositionScope::RcAnchored,
        vec![]
    )]
    fn test_reaction_ast_compose(
        #[case] a: ReactionAst,
        #[case] b: ReactionAst,
        #[case] scope: CompositionScope,
        #[case] expected: Vec<ReactionAst>,
    ) {
        assert_eq!(a.compose(&b, scope), expected);
    }

    #[rstest]
    fn test_reaction_ast_compose_apply_equivalence() {
        // compose(A,B).apply(H) == B.apply(A.apply(H)): C-O 1→2 then 2→3 on host C-O order 1.
        let a = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            })]),
        );
        let b = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(3) },
            })]),
        );
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        );

        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|c| c.apply(&host, SubgraphIsomorphismAlgorithm::Vf2))
            .collect();
        let sequential: Vec<MoleculeAst> = a
            .apply(&host, SubgraphIsomorphismAlgorithm::Vf2)
            .flat_map(|intermediate| {
                b.apply(&intermediate, SubgraphIsomorphismAlgorithm::Vf2)
                    .collect::<Vec<_>>()
            })
            .collect();

        let product = MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(3))],
        );
        assert_eq!(composed, vec![product.clone()]);
        assert_eq!(sequential, vec![product]);
    }

}
