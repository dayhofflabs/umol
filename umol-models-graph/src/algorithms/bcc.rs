use std::collections::HashSet;

pub fn biconnected_components(node_count: usize, adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    if node_count == 0 {
        return Vec::new();
    }

    let mut disc: Vec<Option<u32>> = vec![None; node_count];
    let mut low: Vec<u32> = vec![0; node_count];
    let mut timer: u32 = 0;
    let mut edge_stack: Vec<(usize, usize)> = Vec::new();
    let mut components: Vec<Vec<usize>> = Vec::new();

    for start in 0..node_count {
        if disc[start].is_some() {
            continue;
        }
        bcc_dfs(
            start,
            None,
            adj,
            &mut disc,
            &mut low,
            &mut timer,
            &mut edge_stack,
            &mut components,
        );
    }

    components
}

fn bcc_dfs(
    u: usize,
    parent: Option<usize>,
    adj: &[Vec<usize>],
    disc: &mut [Option<u32>],
    low: &mut [u32],
    timer: &mut u32,
    edge_stack: &mut Vec<(usize, usize)>,
    components: &mut Vec<Vec<usize>>,
) {
    disc[u] = Some(*timer);
    low[u] = *timer;
    *timer += 1;

    let mut child_count = 0u32;
    for &v in &adj[u] {
        if v >= adj.len() {
            continue;
        }
        if disc[v].is_none() {
            child_count += 1;
            edge_stack.push((u, v));
            bcc_dfs(v, Some(u), adj, disc, low, timer, edge_stack, components);

            if low[v] < low[u] {
                low[u] = low[v];
            }

            let is_articulation = match parent {
                None => child_count > 1,
                Some(_) => low[v] >= disc[u].expect("u must be discovered"),
            };
            if is_articulation {
                let mut component_nodes = HashSet::new();
                while let Some((a, b)) = edge_stack.pop() {
                    component_nodes.insert(a);
                    component_nodes.insert(b);
                    if (a == u && b == v) || (a == v && b == u) {
                        break;
                    }
                }
                if component_nodes.len() >= 3 {
                    let mut component: Vec<usize> = component_nodes.into_iter().collect();
                    component.sort_unstable();
                    components.push(component);
                }
            }
        } else if Some(v) != parent
            && disc[v].expect("v must be discovered") < disc[u].expect("u discovered")
        {
            edge_stack.push((u, v));
            let disc_v = disc[v].expect("v must be discovered");
            if disc_v < low[u] {
                low[u] = disc_v;
            }
        }
    }

    if parent.is_none() && !edge_stack.is_empty() {
        let mut component_nodes = HashSet::new();
        while let Some((a, b)) = edge_stack.pop() {
            component_nodes.insert(a);
            component_nodes.insert(b);
        }
        if component_nodes.len() >= 3 {
            let mut component: Vec<usize> = component_nodes.into_iter().collect();
            component.sort_unstable();
            components.push(component);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::biconnected_components;

    #[test]
    fn empty_graph() {
        assert!(biconnected_components(0, &[]).is_empty());
    }

    #[test]
    fn chain_has_no_ring_components() {
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        assert!(biconnected_components(4, &adj).is_empty());
    }

    #[test]
    fn single_cycle_component() {
        let adj = vec![vec![1, 3], vec![0, 2], vec![1, 3], vec![2, 0]];
        let components = biconnected_components(4, &adj);
        assert_eq!(components, vec![vec![0, 1, 2, 3]]);
    }

    #[test]
    fn articulation_splits_components() {
        // two cycles sharing one articulation node (2)
        let adj = vec![
            vec![1, 2],
            vec![0, 2],
            vec![0, 1, 3, 4],
            vec![2, 4],
            vec![2, 3],
        ];
        let mut components = biconnected_components(5, &adj);
        components.sort();
        assert_eq!(components, vec![vec![0, 1, 2], vec![2, 3, 4]]);
    }

    #[test]
    fn disconnected_graph() {
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1], vec![4], vec![3]];
        let components = biconnected_components(5, &adj);
        assert_eq!(components, vec![vec![0, 1, 2]]);
    }
}
