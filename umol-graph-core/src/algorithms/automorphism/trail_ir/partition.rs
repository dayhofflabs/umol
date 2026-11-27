//! Zero-copy edge-colored CSR view and trail-backed ordered vertex partition.

use crate::graph::Graph;

#[inline]
pub(super) fn common_prefix_len(a: &[u32], b: &[u32]) -> usize {
    let limit = a.len().min(b.len());
    let mut i = 0;
    while i < limit && a[i] == b[i] {
        i += 1;
    }
    i
}

#[derive(Clone, Copy)]
pub(super) struct EdgeColoredView<'a> {
    pub(super) graph: &'a Graph,
    pub(super) edge_colors: &'a [u32],
    pub(super) uniform_edges: bool,
}

impl<'a> EdgeColoredView<'a> {
    pub(super) fn new(graph: &'a Graph, edge_colors: &'a [u32]) -> Self {
        debug_assert_eq!(edge_colors.len(), graph.edge_count());
        let uniform_edges = edge_colors
            .first()
            .is_none_or(|&first| edge_colors.iter().all(|&c| c == first));
        Self {
            graph,
            edge_colors,
            uniform_edges,
        }
    }

    #[inline]
    pub(super) fn node_count(&self) -> usize {
        self.graph.node_count()
    }
}

#[derive(Clone, Copy, Debug)]
struct TrailEntry {
    start: u32,
    old_len: u32,
    added_cells: u32,
}

#[derive(Clone, Debug, Default)]
pub(super) struct Partition {
    pub(super) lab: Vec<u32>,
    pub(super) inv_lab: Vec<u32>,
    pub(super) cell_start: Vec<u32>,
    pub(super) cell_len: Vec<u32>,
    pub(super) num_cells: usize,
    trail: Vec<TrailEntry>,
}

impl Partition {
    pub(super) fn reset_from_colors(&mut self, node_colors: &[u32], initial_cells: &mut Vec<u32>) {
        let n = node_colors.len();
        self.lab.clear();
        self.lab.extend(0..n as u32);
        self.inv_lab.resize(n, 0);
        self.cell_start.resize(n, 0);
        self.cell_len.resize(n, 0);
        self.trail.clear();
        initial_cells.clear();

        if n == 0 {
            self.num_cells = 0;
            return;
        }

        if n <= 24 {
            for i in 1..n {
                let v = self.lab[i];
                let cv = node_colors[v as usize];
                let mut j = i;
                while j > 0 && node_colors[self.lab[j - 1] as usize] > cv {
                    self.lab[j] = self.lab[j - 1];
                    j -= 1;
                }
                self.lab[j] = v;
            }
        } else {
            self.lab.sort_unstable_by_key(|&v| node_colors[v as usize]);
        }

        let mut start = 0usize;
        let mut cells = 0usize;
        while start < n {
            let color = node_colors[self.lab[start] as usize];
            let mut end = start + 1;
            while end < n && node_colors[self.lab[end] as usize] == color {
                end += 1;
            }
            let len = (end - start) as u32;
            let start_u32 = start as u32;
            self.cell_len[start] = len;
            for (pos, &v) in self.lab[start..end].iter().enumerate() {
                let vi = v as usize;
                self.cell_start[vi] = start_u32;
                self.inv_lab[vi] = (start + pos) as u32;
            }
            initial_cells.push(start_u32);
            cells += 1;
            start = end;
        }
        self.num_cells = cells;
    }

    #[inline]
    pub(super) fn swap_positions(&mut self, i: usize, j: usize) {
        if i != j {
            let vi = self.lab[i];
            let vj = self.lab[j];
            self.lab[i] = vj;
            self.lab[j] = vi;
            self.inv_lab[vj as usize] = i as u32;
            self.inv_lab[vi as usize] = j as u32;
        }
    }

    #[inline]
    pub(super) fn refresh_inv_range(&mut self, start: usize, end: usize) {
        for pos in start..end {
            self.inv_lab[self.lab[pos] as usize] = pos as u32;
        }
    }

    #[inline]
    pub(super) fn is_discrete(&self) -> bool {
        self.num_cells == self.lab.len()
    }

    #[inline]
    pub(super) fn trail_mark(&self) -> usize {
        self.trail.len()
    }

    #[inline]
    pub(super) fn record_split(&mut self, start: u32, old_len: u32, added_cells: u32) {
        debug_assert!(added_cells > 0);
        self.trail.push(TrailEntry {
            start,
            old_len,
            added_cells,
        });
        self.num_cells += added_cells as usize;
    }

    pub(super) fn individualize(&mut self, start: usize, vertex: u32) {
        let len = self.cell_len[start] as usize;
        debug_assert!(len > 1);
        let pos = self.inv_lab[vertex as usize] as usize;
        debug_assert!(pos >= start && pos < start + len);
        self.swap_positions(start, pos);

        self.cell_len[start] = 1;
        self.cell_len[start + 1] = (len - 1) as u32;
        let rem_start = (start + 1) as u32;
        for &u in &self.lab[start + 1..start + len] {
            self.cell_start[u as usize] = rem_start;
        }
        self.record_split(start as u32, len as u32, 1);
    }

    pub(super) fn collapse_cell_to_singletons(&mut self, start: usize) {
        let len = self.cell_len[start] as usize;
        debug_assert!(len > 1);
        for pos in start..start + len {
            self.cell_len[pos] = 1;
            self.cell_start[self.lab[pos] as usize] = pos as u32;
        }
        self.record_split(start as u32, len as u32, (len - 1) as u32);
    }

    pub(super) fn unwind(&mut self, mark: usize) {
        while self.trail.len() > mark {
            let entry = self.trail.pop().expect("trail non-empty above mark");
            let start = entry.start as usize;
            let end = start + entry.old_len as usize;
            self.cell_len[start] = entry.old_len;
            for &v in &self.lab[start..end] {
                self.cell_start[v as usize] = entry.start;
            }
            self.num_cells -= entry.added_cells as usize;
        }
    }

    #[inline]
    pub(super) fn cell_slice(&self, start: usize) -> &[u32] {
        &self.lab[start..start + self.cell_len[start] as usize]
    }
}
