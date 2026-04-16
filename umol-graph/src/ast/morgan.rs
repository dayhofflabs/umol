//! Morgan / ECFP fingerprints over MoleculeAst.
//!
//! Two implementations for benchmarking:
//! - `morgan_direct`: works directly on `MoleculeAst` fields (Option access)
//! - `morgan_view`: works on a pre-built `MorganTarget` (packed arrays)
//!
//! Algorithm: ECFP (Rogers & Hahn, J. Chem. Inf. Model. 2010, 50, 742-754).
//!
//! Initial atom invariants (Daylight invariants + ring flag):
//!   1. heavy-atom degree (count of non-H neighbors)
//!   2. heavy-atom valence (bond order sum to non-H neighbors)
//!   3. atomic number
//!   4. atomic mass (isotope mass, or 0 for natural)
//!   5. formal charge
//!   6. total hydrogen count (implicit + explicit)
//!   7. ring membership flag (0/1)
//!
//! Iterative update: hash array [iteration, current_id, bond_order_1,
//! neighbor_id_1, ...] with neighbors sorted by (bond_order, neighbor_id).
//!
//! Duplicate structure removal: each atom tracks its bond environment as a
//! bitset. If two features share the same bond set, the one from more
//! iterations is removed (tie: larger hash removed). Atoms whose environment
//! duplicates a previously seen one are marked dead and skipped.

use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::mem;

use fixedbitset::FixedBitSet;
use xxhash_rust::xxh3::{Xxh3, Xxh3DefaultBuilder};

use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::value_ast::ValueAst;

use crate::ast::{AtomIdx, BondIdx};
use crate::ast::molecule::MoleculeAst;

type XxHashMap<K, V> = HashMap<K, V, Xxh3DefaultBuilder>;

/// Sparse Morgan fingerprint: feature hash → count.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MorganFingerprint {
    pub features: XxHashMap<u32, u32>,
}

impl MorganFingerprint {
    fn new() -> Self {
        Self {
            features: XxHashMap::with_hasher(Xxh3DefaultBuilder),
        }
    }
}

// ---------------------------------------------------------------------------
// Hashing
// ---------------------------------------------------------------------------

fn morgan_hash(values: &[u32]) -> u32 {
    let mut h = Xxh3::with_seed(0);
    for v in values {
        v.hash(&mut h);
    }
    h.finish() as u32
}

// ---------------------------------------------------------------------------
// Bond environment bitset
// ---------------------------------------------------------------------------

/// Compact bitset over bond indices, used for duplicate detection.
#[derive(Clone, PartialEq, Eq, Hash)]
struct BondSet {
    words: Vec<u64>,
}

impl BondSet {
    fn new(bond_count: usize) -> Self {
        let word_count = bond_count.div_ceil(64);
        Self {
            words: vec![0; word_count],
        }
    }

    fn set(&mut self, idx: usize) {
        self.words[idx / 64] |= 1u64 << (idx % 64);
    }

    fn union_with(&mut self, other: &BondSet) {
        for (a, b) in self.words.iter_mut().zip(other.words.iter()) {
            *a |= *b;
        }
    }
}

// ---------------------------------------------------------------------------
// Atom invariant extraction from MoleculeAst
// ---------------------------------------------------------------------------

fn is_heavy(ast: &MoleculeAst, idx: AtomIdx) -> bool {
    !matches!(&ast.atom(idx).element, ElementAst::Lit(Element::H))
}

fn ast_atomic_number(ast: &MoleculeAst, idx: AtomIdx) -> u32 {
    match &ast.atom(idx).element {
        ElementAst::Lit(e) => e.atomic_number() as u32,
        _ => 0,
    }
}

fn ast_atomic_mass(ast: &MoleculeAst, idx: AtomIdx) -> u32 {
    match &ast.atom(idx).isotope_mass {
        IsotopeAst::Lit(m) => *m,
        _ => 0,
    }
}

fn ast_charge(ast: &MoleculeAst, idx: AtomIdx) -> i32 {
    match &ast.atom(idx).charge {
        ValueAst::Lit(n) => *n as i32,
        _ => 0,
    }
}

fn ast_h_count(ast: &MoleculeAst, idx: AtomIdx) -> u32 {
    match &ast.atom(idx).implicit_hydrogens {
        HydrogenAst::Value(ValueAst::Lit(n)) => *n as u32,
        _ => 0,
    }
}

fn bond_order(ast: &MoleculeAst, bond_idx: BondIdx) -> u32 {
    match &ast.bond(bond_idx).order {
        ValueAst::Lit(n) => *n as u32,
        _ => 1,
    }
}

struct AtomInvariants {
    heavy_degree: u32,
    heavy_valence: u32,
    h_count: u32,
}

fn compute_atom_invariants(ast: &MoleculeAst) -> Vec<AtomInvariants> {
    let n = ast.atom_count();
    let mut result: Vec<AtomInvariants> = (0..n)
        .map(|_| AtomInvariants {
            heavy_degree: 0,
            heavy_valence: 0,
            h_count: 0,
        })
        .collect();

    for (bi, source, target, _) in ast.bonds() {
        let order = bond_order(ast, bi);
        if is_heavy(ast, target) {
            result[source.index()].heavy_degree += 1;
            result[source.index()].heavy_valence += order;
        } else {
            result[source.index()].h_count += 1;
        }
        if is_heavy(ast, source) {
            result[target.index()].heavy_degree += 1;
            result[target.index()].heavy_valence += order;
        } else {
            result[target.index()].h_count += 1;
        }
    }

    // Add implicit H count
    for (i, entry) in result.iter_mut().enumerate().take(n) {
        entry.h_count += ast_h_count(ast, AtomIdx::from_usize(i));
    }

    result
}

/// Compute ring membership flags via DFS cycle detection.
fn compute_ring_flags(ast: &MoleculeAst) -> Vec<bool> {
    let n = ast.atom_count();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, source, target, _) in ast.bonds() {
        adj[source.index()].push((target.index(), bi.index()));
        adj[target.index()].push((source.index(), bi.index()));
    }

    let mut in_ring = vec![false; n];
    let mut visited = vec![false; n];
    let mut parent_edge = vec![usize::MAX; n];
    let mut stack: Vec<(usize, usize)> = Vec::new();
    let mut path: Vec<usize> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        stack.push((start, 0));
        path.push(start);

        while let Some((node, cursor)) = stack.last_mut() {
            if *cursor < adj[*node].len() {
                let (neighbor, edge_idx) = adj[*node][*cursor];
                *cursor += 1;
                if edge_idx == parent_edge[*node] {
                    continue;
                }
                if visited[neighbor] {
                    if let Some(pos) = path.iter().rposition(|&p| p == neighbor) {
                        for &ring_atom in &path[pos..] {
                            in_ring[ring_atom] = true;
                        }
                    }
                } else {
                    visited[neighbor] = true;
                    parent_edge[neighbor] = edge_idx;
                    stack.push((neighbor, 0));
                    path.push(neighbor);
                }
            } else {
                stack.pop();
                path.pop();
            }
        }
    }
    in_ring
}

// ---------------------------------------------------------------------------
// Shared ECFP core
// ---------------------------------------------------------------------------

fn initial_identifier(
    atomic_number: u32,
    atomic_mass: u32,
    charge: i32,
    inv: &AtomInvariants,
    in_ring: bool,
) -> u32 {
    morgan_hash(&[
        inv.heavy_degree,
        inv.heavy_valence,
        atomic_number,
        atomic_mass,
        (charge + 4) as u32, // shift to avoid negative
        inv.h_count,
        in_ring as u32,
    ])
}

/// ECFP core loop. Generic over how atom/bond data is accessed.
///
/// `adj[i]` is a slice of `(neighbor_atom_idx, bond_vec_idx)`.
/// `bond_orders[bond_vec_idx]` gives the integer bond order.
fn ecfp_loop(
    n: usize,
    bond_count: usize,
    initial_ids: &[u32],
    adj: &[&[(usize, usize)]],
    bond_orders: &[u32],
    radius: usize,
) -> MorganFingerprint {
    let mut fp = MorganFingerprint::new();
    if n == 0 {
        return fp;
    }

    let mut identifiers = initial_ids.to_vec();
    let mut dead = vec![false; n];

    // Track bond environments per atom
    let mut atom_envs: Vec<BondSet> = (0..n).map(|_| BondSet::new(bond_count)).collect();
    // Global set of already-seen environments
    let mut seen_envs: HashSet<BondSet> = HashSet::new();

    // Radius 0: add initial features (no duplicates possible at radius 0)
    for &id in &identifiers {
        *fp.features.entry(id).or_insert(0) += 1;
    }

    let mut new_ids = vec![0u32; n];
    let mut round_envs: Vec<BondSet> = (0..n).map(|_| BondSet::new(bond_count)).collect();
    let mut neighbor_pairs: Vec<(u32, u32)> = Vec::new();

    for layer in 0..radius {
        // Compute new identifiers and environments
        let mut round_features: Vec<(BondSet, u32, usize)> = Vec::new(); // (env, id, atom)

        for i in 0..n {
            if dead[i] {
                new_ids[i] = identifiers[i];
                round_envs[i] = atom_envs[i].clone();
                continue;
            }

            if adj[i].is_empty() {
                dead[i] = true;
                new_ids[i] = identifiers[i];
                round_envs[i] = atom_envs[i].clone();
                continue;
            }

            // Build update array: [iteration, current_id, bond_order, neighbor_id, ...]
            neighbor_pairs.clear();
            let mut new_env = atom_envs[i].clone();

            for &(neighbor, bond_idx) in adj[i] {
                neighbor_pairs.push((bond_orders[bond_idx], identifiers[neighbor]));
                new_env.set(bond_idx);
                new_env.union_with(&atom_envs[neighbor]);
            }
            neighbor_pairs.sort_unstable();

            let mut hash_input: Vec<u32> = Vec::with_capacity(2 + neighbor_pairs.len() * 2);
            hash_input.push(layer as u32);
            hash_input.push(identifiers[i]);
            for &(bo, nid) in &neighbor_pairs {
                hash_input.push(bo);
                hash_input.push(nid);
            }
            new_ids[i] = morgan_hash(&hash_input);
            round_envs[i] = new_env.clone();

            round_features.push((new_env, new_ids[i], i));
        }

        // Sort round features for deterministic duplicate resolution
        round_features.sort_by(|a, b| a.0.words.cmp(&b.0.words).then(a.1.cmp(&b.1)));

        // Duplicate removal: check each feature's bond set against seen environments
        for (env, id, atom_idx) in &round_features {
            if dead[*atom_idx] {
                continue;
            }
            if seen_envs.contains(env) {
                dead[*atom_idx] = true;
            } else {
                *fp.features.entry(*id).or_insert(0) += 1;
                seen_envs.insert(env.clone());
            }
        }

        // Advance state
        identifiers.copy_from_slice(&new_ids);
        mem::swap(&mut atom_envs, &mut round_envs);
    }

    fp
}

// ---------------------------------------------------------------------------
// Implementation A: direct over MoleculeAst
// ---------------------------------------------------------------------------

/// Morgan fingerprint directly over MoleculeAst fields.
pub fn morgan_direct(ast: &MoleculeAst, radius: usize) -> MorganFingerprint {
    let n = ast.atom_count();
    if n == 0 {
        return MorganFingerprint::new();
    }

    let ring_flags = compute_ring_flags(ast);
    let atom_invs = compute_atom_invariants(ast);

    let initial_ids: Vec<u32> = (0..n)
        .map(|i| {
            let idx = AtomIdx::from_usize(i);
            initial_identifier(
                ast_atomic_number(ast, idx),
                ast_atomic_mass(ast, idx),
                ast_charge(ast, idx),
                &atom_invs[i],
                ring_flags[i],
            )
        })
        .collect();

    // Build adjacency: (neighbor_atom, bond_vec_index)
    let mut adj_vecs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, source, target, _) in ast.bonds() {
        adj_vecs[source.index()].push((target.index(), bi.index()));
        adj_vecs[target.index()].push((source.index(), bi.index()));
    }
    let adj_slices: Vec<&[(usize, usize)]> = adj_vecs.iter().map(|v| v.as_slice()).collect();

    let bond_orders: Vec<u32> = (0..ast.bond_count())
        .map(|bi| bond_order(ast, BondIdx::from_usize(bi)))
        .collect();

    ecfp_loop(n, ast.bond_count(), &initial_ids, &adj_slices, &bond_orders, radius)
}

// ---------------------------------------------------------------------------
// Implementation B: MorganTarget packed view
// ---------------------------------------------------------------------------

/// Precomputed packed view of a MoleculeAst for Morgan fingerprinting.
pub struct MorganTarget {
    atom_count: usize,
    bond_count: usize,
    initial_ids: Vec<u32>,
    /// Per-atom adjacency: (neighbor_atom, bond_vec_index).
    /// Flattened with CSR offsets.
    adj: Vec<(usize, usize)>,
    offsets: Vec<usize>,
    bond_orders: Vec<u32>,
}

impl MorganTarget {
    pub fn new(ast: &MoleculeAst) -> Self {
        let n = ast.atom_count();
        let m = ast.bond_count();

        let ring_flags = compute_ring_flags(ast);
        let atom_invs = compute_atom_invariants(ast);

        let initial_ids: Vec<u32> = (0..n)
            .map(|i| {
                let idx = AtomIdx::from_usize(i);
                initial_identifier(
                    ast_atomic_number(ast, idx),
                    ast_atomic_mass(ast, idx),
                    ast_charge(ast, idx),
                    &atom_invs[i],
                    ring_flags[i],
                )
            })
            .collect();

        let bond_orders: Vec<u32> = (0..m)
            .map(|bi| bond_order(ast, BondIdx::from_usize(bi)))
            .collect();

        // Build CSR adjacency
        let mut adj_lists: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
        for (bi, source, target, _) in ast.bonds() {
            adj_lists[source.index()].push((target.index(), bi.index()));
            adj_lists[target.index()].push((source.index(), bi.index()));
        }

        let total: usize = adj_lists.iter().map(|a| a.len()).sum();
        let mut adj = Vec::with_capacity(total);
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        for list in &adj_lists {
            adj.extend_from_slice(list);
            offsets.push(adj.len());
        }

        Self {
            atom_count: n,
            bond_count: m,
            initial_ids,
            adj,
            offsets,
            bond_orders,
        }
    }
}

/// Morgan fingerprint over a precomputed MorganTarget view.
pub fn morgan_view(target: &MorganTarget, radius: usize) -> MorganFingerprint {
    let n = target.atom_count;
    if n == 0 {
        return MorganFingerprint::new();
    }

    let adj_slices: Vec<&[(usize, usize)]> = (0..n)
        .map(|i| &target.adj[target.offsets[i]..target.offsets[i + 1]])
        .collect();

    ecfp_loop(
        n,
        target.bond_count,
        &target.initial_ids,
        &adj_slices,
        &target.bond_orders,
        radius,
    )
}

// ---------------------------------------------------------------------------
// Implementation C: optimized MorganTarget
// ---------------------------------------------------------------------------

/// Optimized packed view with pre-sorted adjacency and FixedBitSet environments.
pub struct MorganTargetOpt {
    atom_count: usize,
    bond_count: usize,
    initial_ids: Vec<u32>,
    /// CSR adjacency, pre-sorted by (bond_order, 0) per atom.
    /// Entries: (neighbor_atom, bond_vec_index, bond_order).
    adj: Vec<(u32, u32, u32)>,
    offsets: Vec<u32>,
}

impl MorganTargetOpt {
    pub fn new(ast: &MoleculeAst) -> Self {
        let n = ast.atom_count();
        let m = ast.bond_count();

        let ring_flags = compute_ring_flags(ast);
        let atom_invs = compute_atom_invariants(ast);

        let initial_ids: Vec<u32> = (0..n)
            .map(|i| {
                let idx = AtomIdx::from_usize(i);
                initial_identifier(
                    ast_atomic_number(ast, idx),
                    ast_atomic_mass(ast, idx),
                    ast_charge(ast, idx),
                    &atom_invs[i],
                    ring_flags[i],
                )
            })
            .collect();

        let bond_orders: Vec<u32> = (0..m)
            .map(|bi| bond_order(ast, BondIdx::from_usize(bi)))
            .collect();

        // Build per-atom adjacency lists, pre-sorted by bond_order
        let mut adj_lists: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); n];
        for (bi, source, target, _) in ast.bonds() {
            let bo = bond_orders[bi.index()];
            adj_lists[source.index()].push((target.0, bi.0, bo));
            adj_lists[target.index()].push((source.0, bi.0, bo));
        }
        for list in &mut adj_lists {
            list.sort_unstable_by_key(|&(_, _, bo)| bo);
        }

        // Flatten to CSR
        let total: usize = adj_lists.iter().map(|a| a.len()).sum();
        let mut adj = Vec::with_capacity(total);
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0u32);
        for list in &adj_lists {
            adj.extend_from_slice(list);
            offsets.push(adj.len() as u32);
        }

        Self {
            atom_count: n,
            bond_count: m,
            initial_ids,
            adj,
            offsets,
        }
    }
}

/// u64-based ECFP loop for molecules with ≤ 64 bonds.
fn ecfp_opt_small(target: &MorganTargetOpt, radius: usize) -> MorganFingerprint {
    let n = target.atom_count;
    let mut fp = MorganFingerprint::new();
    let mut identifiers = target.initial_ids.clone();
    let mut dead = vec![false; n];

    let mut atom_envs = vec![0u64; n];
    let mut round_envs = vec![0u64; n];
    let mut seen_envs: HashSet<u64> = HashSet::new();

    for &id in &identifiers {
        *fp.features.entry(id).or_insert(0) += 1;
    }

    let mut new_ids = vec![0u32; n];
    let mut neighbor_pairs: Vec<(u32, u32)> = Vec::with_capacity(16);
    let mut hash_input: Vec<u32> = Vec::with_capacity(32);
    let mut round_features: Vec<(usize, u32)> = Vec::with_capacity(n);

    for layer in 0..radius {
        round_features.clear();

        for i in 0..n {
            if dead[i] {
                new_ids[i] = identifiers[i];
                round_envs[i] = atom_envs[i];
                continue;
            }

            let start = target.offsets[i] as usize;
            let end = target.offsets[i + 1] as usize;
            if start == end {
                dead[i] = true;
                new_ids[i] = identifiers[i];
                round_envs[i] = atom_envs[i];
                continue;
            }

            let mut env = atom_envs[i];
            neighbor_pairs.clear();
            for &(neighbor, bond_idx, bo) in &target.adj[start..end] {
                env |= 1u64 << bond_idx;
                env |= atom_envs[neighbor as usize];
                neighbor_pairs.push((bo, identifiers[neighbor as usize]));
            }
            round_envs[i] = env;

            neighbor_pairs.sort_unstable();

            hash_input.clear();
            hash_input.push(layer as u32);
            hash_input.push(identifiers[i]);
            for &(bo, nid) in &neighbor_pairs {
                hash_input.push(bo);
                hash_input.push(nid);
            }
            new_ids[i] = morgan_hash(&hash_input);

            round_features.push((i, new_ids[i]));
        }

        round_features.sort_unstable_by(|a, b| {
            round_envs[a.0].cmp(&round_envs[b.0]).then(a.1.cmp(&b.1))
        });

        for &(atom_idx, id) in &round_features {
            if dead[atom_idx] {
                continue;
            }
            if seen_envs.contains(&round_envs[atom_idx]) {
                dead[atom_idx] = true;
            } else {
                *fp.features.entry(id).or_insert(0) += 1;
                seen_envs.insert(round_envs[atom_idx]);
            }
        }

        identifiers.copy_from_slice(&new_ids);
        mem::swap(&mut atom_envs, &mut round_envs);
    }

    fp
}

/// Optimized Morgan fingerprint with u64 fast path for ≤ 64 bonds,
/// FixedBitSet fallback otherwise.
pub fn morgan_view_opt(target: &MorganTargetOpt, radius: usize) -> MorganFingerprint {
    let n = target.atom_count;
    if n == 0 {
        return MorganFingerprint::new();
    }
    if target.bond_count <= 64 {
        return ecfp_opt_small(target, radius);
    }
    let m = target.bond_count;

    let mut fp = MorganFingerprint::new();
    let mut identifiers = target.initial_ids.clone();
    let mut dead = vec![false; n];

    // FixedBitSet environments — one per atom, plus a scratch set
    let mut atom_envs: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(m)).collect();
    let mut round_envs: Vec<FixedBitSet> =
        (0..n).map(|_| FixedBitSet::with_capacity(m)).collect();
    let mut seen_envs: HashSet<FixedBitSet> = HashSet::new();

    // Radius 0
    for &id in &identifiers {
        *fp.features.entry(id).or_insert(0) += 1;
    }

    let mut new_ids = vec![0u32; n];
    // Preallocated scratch buffers
    let mut neighbor_pairs: Vec<(u32, u32)> = Vec::with_capacity(16);
    let mut hash_input: Vec<u32> = Vec::with_capacity(32);
    let mut round_features: Vec<(usize, u32)> = Vec::with_capacity(n);

    for layer in 0..radius {
        round_features.clear();

        for i in 0..n {
            if dead[i] {
                new_ids[i] = identifiers[i];
                round_envs[i].clone_from(&atom_envs[i]);
                continue;
            }

            let start = target.offsets[i] as usize;
            let end = target.offsets[i + 1] as usize;
            if start == end {
                dead[i] = true;
                new_ids[i] = identifiers[i];
                round_envs[i].clone_from(&atom_envs[i]);
                continue;
            }

            // Build environment: union of previous env + attachment bonds + neighbor envs
            round_envs[i].clone_from(&atom_envs[i]);
            neighbor_pairs.clear();
            for &(neighbor, bond_idx, bo) in &target.adj[start..end] {
                round_envs[i].insert(bond_idx as usize);
                round_envs[i].union_with(&atom_envs[neighbor as usize]);
                neighbor_pairs.push((bo, identifiers[neighbor as usize]));
            }

            // Adjacency is pre-sorted by bond_order; stable-sort by full pair
            // to get (bond_order, neighbor_id) order. Since bond_order is already
            // sorted, this only needs to reorder within same-bond-order groups.
            neighbor_pairs.sort_unstable();

            // Hash: [layer, current_id, bo1, nid1, bo2, nid2, ...]
            hash_input.clear();
            hash_input.push(layer as u32);
            hash_input.push(identifiers[i]);
            for &(bo, nid) in &neighbor_pairs {
                hash_input.push(bo);
                hash_input.push(nid);
            }
            new_ids[i] = morgan_hash(&hash_input);

            round_features.push((i, new_ids[i]));
        }

        // Sort for deterministic duplicate resolution (same env → smaller hash wins)
        round_features.sort_unstable_by(|a, b| {
            round_envs[a.0]
                .as_slice()
                .cmp(round_envs[b.0].as_slice())
                .then(a.1.cmp(&b.1))
        });

        // Duplicate removal
        for &(atom_idx, id) in &round_features {
            if dead[atom_idx] {
                continue;
            }
            if seen_envs.contains(&round_envs[atom_idx]) {
                dead[atom_idx] = true;
            } else {
                *fp.features.entry(id).or_insert(0) += 1;
                seen_envs.insert(round_envs[atom_idx].clone());
            }
        }

        identifiers.copy_from_slice(&new_ids);
        mem::swap(&mut atom_envs, &mut round_envs);
    }

    fp
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;

    fn carbon(h: u8) -> AtomAst {
        AtomAst {
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(h as i64)),
            ..AtomAst::from_element(Element::C)
        }
    }

    fn methane() -> MoleculeAst {
        MoleculeAst::new(vec![carbon(4)], vec![], vec![], vec![], vec![], vec![], vec![])
    }

    fn ethane() -> MoleculeAst {
        MoleculeAst::new(
            vec![carbon(3), carbon(3)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![], vec![],
        )
    }

    fn ethanol() -> MoleculeAst {
        MoleculeAst::new(
            vec![
                carbon(3),
                carbon(2),
                AtomAst {
                    implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(1)),
                    ..AtomAst::from_element(Element::O)
                },
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            ],
            vec![], vec![], vec![], vec![], vec![],
        )
    }

    fn cyclohexane() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| carbon(2)).collect();
        let bonds: Vec<(AtomIdx, AtomIdx, BondAst)> = (0..6)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], vec![])
    }

    fn formaldehyde() -> MoleculeAst {
        MoleculeAst::new(
            vec![carbon(2), AtomAst::from_element(Element::O)],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![], vec![], vec![], vec![], vec![],
        )
    }

    #[rstest]
    #[case::methane(methane())]
    #[case::ethane(ethane())]
    #[case::ethanol(ethanol())]
    #[case::cyclohexane(cyclohexane())]
    #[case::formaldehyde(formaldehyde())]
    fn test_morgan_direct_eq_view(#[case] ast: MoleculeAst) {
        for radius in 0..=3 {
            let direct = morgan_direct(&ast, radius);
            let target = MorganTarget::new(&ast);
            let view = morgan_view(&target, radius);
            assert_eq!(direct.features, view.features, "radius {radius}");
        }
    }

    #[rstest]
    #[case::methane(methane())]
    #[case::ethane(ethane())]
    #[case::ethanol(ethanol())]
    #[case::cyclohexane(cyclohexane())]
    #[case::formaldehyde(formaldehyde())]
    fn test_morgan_direct_eq_view_opt(#[case] ast: MoleculeAst) {
        for radius in 0..=3 {
            let direct = morgan_direct(&ast, radius);
            let target = MorganTargetOpt::new(&ast);
            let opt = morgan_view_opt(&target, radius);
            assert_eq!(direct.features, opt.features, "radius {radius}");
        }
    }

    #[test]
    fn test_morgan_empty() {
        let ast = MoleculeAst::default();
        let fp = morgan_direct(&ast, 2);
        assert!(fp.features.is_empty());
    }

    #[test]
    fn test_morgan_single_atom() {
        let fp = morgan_direct(&methane(), 2);
        // No neighbors → same id each radius, but duplicate removal keeps only radius 0
        assert_eq!(fp.features.len(), 1);
    }

    #[test]
    fn test_morgan_ethane_symmetry() {
        // Two equivalent C atoms: at radius 0 they share one identifier (count=2).
        // After radius 1, both see the same environment → duplicate removal
        // keeps only one, so feature counts should be even at radius 0.
        let fp = morgan_direct(&ethane(), 0);
        for &count in fp.features.values() {
            assert_eq!(count % 2, 0);
        }
    }

    #[test]
    fn test_morgan_cyclohexane_ring_flag() {
        let ring_flags = compute_ring_flags(&cyclohexane());
        assert!(ring_flags.iter().all(|&f| f));
    }

    #[test]
    fn test_morgan_chain_no_ring() {
        let ring_flags = compute_ring_flags(&ethanol());
        assert!(ring_flags.iter().all(|&f| !f));
    }

    #[test]
    fn test_morgan_formaldehyde_degree_vs_valence() {
        // C in H2C=O: heavy_degree=1, heavy_valence=2 (double bond)
        let ast = formaldehyde();
        let invs = compute_atom_invariants(&ast);
        assert_eq!(invs[0].heavy_degree, 1);
        assert_eq!(invs[0].heavy_valence, 2);
        assert_eq!(invs[0].h_count, 2);
        // O: heavy_degree=1, heavy_valence=2, h_count=0
        assert_eq!(invs[1].heavy_degree, 1);
        assert_eq!(invs[1].heavy_valence, 2);
        assert_eq!(invs[1].h_count, 0);
    }

    #[test]
    fn test_morgan_duplicate_removal() {
        // Ethane at radius 1: both C atoms see identical bond environments
        // (one C-C bond). The second should be removed as duplicate.
        let fp_r0 = morgan_direct(&ethane(), 0);
        let fp_r1 = morgan_direct(&ethane(), 1);
        // Radius 0 has 1 unique feature (both C atoms identical)
        assert_eq!(fp_r0.features.len(), 1);
        // Radius 1 should have 2 features total: the radius-0 feature
        // plus one radius-1 feature (the duplicate is removed)
        assert_eq!(fp_r1.features.len(), 2);
    }

    #[test]
    fn test_morgan_radius_increases_features() {
        let ast = ethanol();
        let fp0 = morgan_direct(&ast, 0);
        let fp2 = morgan_direct(&ast, 2);
        assert!(fp2.features.len() >= fp0.features.len());
    }
}
