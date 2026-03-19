pub fn maximum_independent_set(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    if n == 0 {
        return Vec::new();
    }

    let mut available = vec![true; n];
    let mut current = Vec::new();
    let mut best = Vec::new();
    branch(adj, &mut available, &mut current, &mut best);
    best.sort_unstable();
    best
}

fn branch(
    adj: &[Vec<usize>],
    available: &mut [bool],
    current: &mut Vec<usize>,
    best: &mut Vec<usize>,
) {
    let remaining = available.iter().filter(|&&v| v).count();
    if current.len() + remaining <= best.len() {
        return;
    }

    let Some(node) = next_available(available) else {
        let mut candidate = current.clone();
        candidate.sort_unstable();
        if candidate.len() > best.len() || (candidate.len() == best.len() && candidate < *best) {
            *best = candidate;
        }
        return;
    };

    // Include branch.
    let mut changed = Vec::new();
    if available[node] {
        changed.push(node);
        available[node] = false;
    }
    for &nbr in &adj[node] {
        if nbr < available.len() && available[nbr] {
            changed.push(nbr);
            available[nbr] = false;
        }
    }
    current.push(node);
    branch(adj, available, current, best);
    current.pop();
    for idx in changed {
        available[idx] = true;
    }

    // Exclude branch.
    available[node] = false;
    branch(adj, available, current, best);
    available[node] = true;
}

fn next_available(available: &[bool]) -> Option<usize> {
    available.iter().position(|&v| v)
}

#[cfg(test)]
mod tests {
    use super::maximum_independent_set;

    #[test]
    fn empty() {
        assert!(maximum_independent_set(&[]).is_empty());
    }

    #[test]
    fn clique() {
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        assert_eq!(maximum_independent_set(&adj), vec![0]);
    }

    #[test]
    fn path_graph() {
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        assert_eq!(maximum_independent_set(&adj), vec![0, 2]);
    }

    #[test]
    fn cycle_graph() {
        let adj = vec![vec![1, 3], vec![0, 2], vec![1, 3], vec![2, 0]];
        assert_eq!(maximum_independent_set(&adj), vec![0, 2]);
    }

    #[test]
    fn deterministic_tie_break() {
        // Two maximum sets of size 2: {0,3} and {1,3}; pick lexicographically smaller.
        let adj = vec![vec![1], vec![0, 2], vec![1], vec![]];
        assert_eq!(maximum_independent_set(&adj), vec![0, 2, 3]);
    }
}
