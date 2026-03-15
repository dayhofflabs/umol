//! Utility functions for GraphIR molecule.

use std::collections::{HashMap, HashSet};

use super::AtomIndex;

pub(crate) fn biconnected_components(
    atoms: impl Iterator<Item = AtomIndex>,
    adj: HashMap<AtomIndex, Vec<AtomIndex>>,
) -> Vec<Vec<AtomIndex>> {
    let atoms: Vec<AtomIndex> = atoms.collect();
    if atoms.is_empty() {
        return Vec::new();
    }

    let mut disc: HashMap<AtomIndex, u32> = HashMap::new();
    let mut low: HashMap<AtomIndex, u32> = HashMap::new();
    let mut timer: u32 = 0;
    let mut edge_stack: Vec<(AtomIndex, AtomIndex)> = Vec::new();
    let mut components: Vec<Vec<AtomIndex>> = Vec::new();

    for &start in &atoms {
        if disc.contains_key(&start) {
            continue;
        }
        biconnected_components_dfs(
            start,
            Option::<AtomIndex>::None,
            &adj,
            &mut disc,
            &mut low,
            &mut timer,
            &mut edge_stack,
            &mut components,
        );
    }

    components
}

fn biconnected_components_dfs(
    u: AtomIndex,
    parent: Option<AtomIndex>,
    adj: &HashMap<AtomIndex, Vec<AtomIndex>>,
    disc: &mut HashMap<AtomIndex, u32>,
    low: &mut HashMap<AtomIndex, u32>,
    timer: &mut u32,
    edge_stack: &mut Vec<(AtomIndex, AtomIndex)>,
    components: &mut Vec<Vec<AtomIndex>>,
) {
    // TODO: rewrite iteratively to avoid potential stack overflow on large/deep graphs.
    disc.insert(u, *timer);
    low.insert(u, *timer);
    *timer += 1;
    let mut child_count = 0u32;

    let neighbors: Vec<AtomIndex> = adj.get(&u).cloned().unwrap_or_default();
    for v in neighbors {
        if !disc.contains_key(&v) {
            child_count += 1;
            edge_stack.push((u, v));
            biconnected_components_dfs(v, Some(u), adj, disc, low, timer, edge_stack, components);

            let low_v = low[&v];
            let low_u = low[&u];
            if low_v < low_u {
                low.insert(u, low_v);
            }

            let is_articulation = match parent {
                None => child_count > 1,
                Some(_) => low_v >= disc[&u],
            };
            if is_articulation {
                let mut component_atoms: HashSet<AtomIndex> = HashSet::new();
                loop {
                    if let Some((a, b)) = edge_stack.pop() {
                        component_atoms.insert(a);
                        component_atoms.insert(b);
                        if (a == u && b == v) || (a == v && b == u) {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                if component_atoms.len() >= 3 {
                    let mut atoms_vec: Vec<AtomIndex> = component_atoms.into_iter().collect();
                    atoms_vec.sort_unstable();
                    components.push(atoms_vec);
                }
            }
        } else if Some(v) != parent && disc[&v] < disc[&u] {
            edge_stack.push((u, v));
            let disc_v = disc[&v];
            let low_u = low[&u];
            if disc_v < low_u {
                low.insert(u, disc_v);
            }
        }
    }

    if parent.is_none() && !edge_stack.is_empty() {
        let mut component_atoms: HashSet<AtomIndex> = HashSet::new();
        while let Some((a, b)) = edge_stack.pop() {
            component_atoms.insert(a);
            component_atoms.insert(b);
        }
        if component_atoms.len() >= 3 {
            let mut atoms_vec: Vec<AtomIndex> = component_atoms.into_iter().collect();
            atoms_vec.sort_unstable();
            components.push(atoms_vec);
        }
    }
}

pub(crate) fn enumerate_rings(
    adj: &HashMap<AtomIndex, Vec<AtomIndex>>,
    max_ring_size: usize,
) -> Vec<Vec<AtomIndex>> {
    if max_ring_size < 3 || adj.len() < 3 {
        return Vec::new();
    }

    let mut sorted_atoms: Vec<AtomIndex> = adj.keys().copied().collect();
    sorted_atoms.sort_unstable();

    let mut raw_rings: Vec<Vec<AtomIndex>> = Vec::new();

    for &start in &sorted_atoms {
        let mut path: Vec<AtomIndex> = vec![start];
        let mut visited: HashSet<AtomIndex> = HashSet::new();
        visited.insert(start);
        find_rings_dfs(
            start,
            start,
            &mut path,
            &mut visited,
            max_ring_size,
            adj,
            &mut raw_rings,
        );
    }

    let mut seen: HashSet<Vec<AtomIndex>> = HashSet::new();
    let mut result: Vec<Vec<AtomIndex>> = Vec::new();
    for ring in raw_rings {
        let normalized = normalize_ring(&ring);
        if seen.insert(normalized.clone()) {
            result.push(normalized);
        }
    }
    result.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    result
}

fn find_rings_dfs(
    start: AtomIndex,
    current: AtomIndex,
    path: &mut Vec<AtomIndex>,
    visited: &mut HashSet<AtomIndex>,
    max_size: usize,
    adj: &HashMap<AtomIndex, Vec<AtomIndex>>,
    rings: &mut Vec<Vec<AtomIndex>>,
) {
    if path.len() > max_size {
        return;
    }
    let neighbors = match adj.get(&current) {
        Some(n) => n,
        None => return,
    };
    for &next in neighbors {
        if next == start && path.len() >= 3 {
            rings.push(path.clone());
        } else if next.index() > start.index() && !visited.contains(&next) && path.len() < max_size
        {
            visited.insert(next);
            path.push(next);
            find_rings_dfs(start, next, path, visited, max_size, adj, rings);
            path.pop();
            visited.remove(&next);
        }
    }
}

fn normalize_ring(ring: &[AtomIndex]) -> Vec<AtomIndex> {
    let n = ring.len();
    debug_assert!(n >= 3);
    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, idx)| idx)
        .unwrap()
        .0;
    let mut rotated: Vec<AtomIndex> = Vec::with_capacity(n);
    for i in 0..n {
        rotated.push(ring[(min_pos + i) % n]);
    }
    if n > 1 && rotated[1] > rotated[n - 1] {
        rotated[1..].reverse();
    }
    rotated
}
