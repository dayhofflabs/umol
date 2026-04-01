use std::collections::HashSet;

pub fn enumerate_simple_cycles(
    node_count: usize,
    adj: &[Vec<usize>],
    max_cycle_size: usize,
) -> Vec<Vec<usize>> {
    if max_cycle_size < 3 || node_count < 3 {
        return Vec::new();
    }

    let mut raw_cycles: Vec<Vec<usize>> = Vec::new();
    for start in 0..node_count {
        let mut path = vec![start];
        let mut visited = HashSet::new();
        visited.insert(start);
        dfs_cycles(
            start,
            start,
            &mut path,
            &mut visited,
            max_cycle_size,
            adj,
            &mut raw_cycles,
        );
    }

    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut result = Vec::new();
    for cycle in raw_cycles {
        let normalized = normalize_cycle(&cycle);
        if seen.insert(normalized.clone()) {
            result.push(normalized);
        }
    }

    result.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    result
}

fn dfs_cycles(
    start: usize,
    current: usize,
    path: &mut Vec<usize>,
    visited: &mut HashSet<usize>,
    max_cycle_size: usize,
    adj: &[Vec<usize>],
    cycles: &mut Vec<Vec<usize>>,
) {
    if path.len() > max_cycle_size {
        return;
    }

    for &next in &adj[current] {
        if next >= adj.len() {
            continue;
        }
        if next == start && path.len() >= 3 {
            cycles.push(path.clone());
        } else if next > start && !visited.contains(&next) && path.len() < max_cycle_size {
            visited.insert(next);
            path.push(next);
            dfs_cycles(start, next, path, visited, max_cycle_size, adj, cycles);
            path.pop();
            visited.remove(&next);
        }
    }
}

fn normalize_cycle(cycle: &[usize]) -> Vec<usize> {
    let n = cycle.len();
    debug_assert!(n >= 3);

    let min_pos = cycle
        .iter()
        .enumerate()
        .min_by_key(|&(_, idx)| idx)
        .expect("non-empty cycle")
        .0;

    let mut rotated = Vec::with_capacity(n);
    for i in 0..n {
        rotated.push(cycle[(min_pos + i) % n]);
    }

    if n > 1 && rotated[1] > rotated[n - 1] {
        rotated[1..].reverse();
    }

    rotated
}

#[cfg(test)]
mod tests {
    use super::enumerate_simple_cycles;

    #[test]
    fn simple_ring() {
        let adj = vec![
            vec![1, 5],
            vec![0, 2],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
            vec![4, 0],
        ];
        let cycles = enumerate_simple_cycles(6, &adj, 6);
        assert_eq!(cycles, vec![vec![0, 1, 2, 3, 4, 5]]);
    }

    #[test]
    fn fused_rings() {
        // naphthalene skeleton
        let adj = vec![
            vec![1, 5],
            vec![0, 2],
            vec![1, 3],
            vec![2, 4, 6],
            vec![3, 5, 9],
            vec![0, 4],
            vec![3, 7],
            vec![6, 8],
            vec![7, 9],
            vec![4, 8],
        ];
        let cycles = enumerate_simple_cycles(10, &adj, 10);
        assert_eq!(cycles.len(), 3);
    }

    #[test]
    fn max_size_cutoff() {
        let adj = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        assert!(enumerate_simple_cycles(5, &adj, 4).is_empty());
        assert_eq!(
            enumerate_simple_cycles(5, &adj, 5),
            vec![vec![0, 1, 2, 3, 4]]
        );
    }

    #[test]
    fn deterministic_ordering() {
        let adj = vec![
            vec![1, 2],
            vec![0, 2],
            vec![0, 1, 3, 4],
            vec![2, 4],
            vec![2, 3],
        ];
        let cycles = enumerate_simple_cycles(5, &adj, 5);
        assert_eq!(cycles, vec![vec![0, 1, 2], vec![2, 3, 4]]);
    }
}
