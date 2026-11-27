//! Non-Uniform Component Recursion (NUCR) for edge-colored partitions.

use std::cmp::Ordering;

use super::partition::{common_prefix_len, EdgeColoredView, Partition};
use super::preprocess::{collapse_twin_cells_in_ranges, TwinCell};
use super::refine::{refine, RefineWorkspace};
use crate::id::NodeId;

#[derive(Clone, Debug)]
pub(super) struct NonUniformComponent {
    pub(super) total_verts: usize,
    pub(super) first_start: u32,
    pub(super) single_uniform_twin: bool,
    pub(super) cells: Vec<(u32, u32)>,
}

#[derive(Clone, Debug, Default)]
pub(super) struct NucWorkspace {
    parent: Vec<u32>,
    internal_uniform: Vec<bool>,
    root_to_comp: Vec<usize>,
    pub(super) components: Vec<NonUniformComponent>,
}

impl NucWorkspace {
    #[inline]
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let p = self.parent[x as usize];
            let g = self.parent[p as usize];
            self.parent[x as usize] = g;
            x = p;
        }
        x
    }

    #[inline]
    fn union(&mut self, a: u32, b: u32, num_comps: &mut usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            if ra < rb {
                self.parent[rb as usize] = ra;
            } else {
                self.parent[ra as usize] = rb;
            }
            *num_comps -= 1;
        }
    }
}

pub(super) fn find_non_uniform_components(
    view: EdgeColoredView<'_>,
    partition: &Partition,
    ws: &mut RefineWorkspace,
    nuc_ws: &mut NucWorkspace,
    min_non_singleton_cells: usize,
) -> bool {
    let n = view.node_count();
    let min_needed = min_non_singleton_cells.max(2);
    if n - partition.num_cells < min_needed {
        return false;
    }

    ws.ensure_capacity(n);
    ws.tc_start.clear();
    ws.tc_size.clear();

    let mut i = 0usize;
    let mut nnt = 0u32;
    while i < n {
        let len = partition.cell_len[i] as usize;
        if len > 1 {
            ws.tc_start.push(i as u32);
            ws.tc_size.push(len as u32);
            for &v in &partition.lab[i..i + len] {
                ws.tc_nntcell[v as usize] = nnt;
            }
            nnt += 1;
        } else {
            ws.tc_nntcell[partition.lab[i] as usize] = u32::MAX;
        }
        i += len;
    }

    let nnt_usize = nnt as usize;
    if nnt_usize < min_needed {
        return false;
    }

    nuc_ws.parent.clear();
    nuc_ws.parent.extend(0..nnt);
    nuc_ws.internal_uniform.clear();
    nuc_ws.internal_uniform.resize(nnt_usize, true);

    ws.tc_hits.clear();
    ws.tc_hits.resize(nnt_usize, 0);
    ws.tc_first_wt.clear();
    ws.tc_first_wt.resize(nnt_usize, 0);
    ws.tc_mixed.clear();
    ws.tc_mixed.resize(nnt_usize, 0);

    let mut num_comps = nnt_usize;

    for ci in 0..nnt_usize {
        let rep = partition.lab[ws.tc_start[ci] as usize];
        for nb in view.graph.neighbors(NodeId(rep)) {
            let k = ws.tc_nntcell[nb.node.index()];
            if k != u32::MAX {
                let ki = k as usize;
                let wt = view.edge_colors[nb.edge.index()];
                if ws.tc_hits[ki] == 0 {
                    ws.tc_first_wt[ki] = wt;
                } else if ws.tc_first_wt[ki] != wt {
                    ws.tc_mixed[ki] = 1;
                }
                ws.tc_hits[ki] += 1;
            }
        }

        for nb in view.graph.neighbors(NodeId(rep)) {
            let k = ws.tc_nntcell[nb.node.index()];
            if k != u32::MAX {
                let ki = k as usize;
                let h = ws.tc_hits[ki];
                if h > 0 {
                    let mixed = ws.tc_mixed[ki] != 0;
                    if ki == ci {
                        let uniform_self = h == ws.tc_size[ki] - 1 && !mixed;
                        nuc_ws.internal_uniform[ci] = uniform_self;
                    } else if h < ws.tc_size[ki] || mixed {
                        nuc_ws.union(ci as u32, ki as u32, &mut num_comps);
                    }
                    ws.tc_hits[ki] = 0;
                    ws.tc_mixed[ki] = 0;
                }
            }
        }

        if num_comps == 1 {
            return false;
        }
    }

    nuc_ws.root_to_comp.clear();
    nuc_ws.root_to_comp.resize(nnt_usize, usize::MAX);
    nuc_ws.components.clear();

    for ci in 0..nnt_usize {
        let root = nuc_ws.find(ci as u32) as usize;
        let idx = if nuc_ws.root_to_comp[root] == usize::MAX {
            let next = nuc_ws.components.len();
            nuc_ws.root_to_comp[root] = next;
            nuc_ws.components.push(NonUniformComponent {
                total_verts: 0,
                first_start: ws.tc_start[ci],
                single_uniform_twin: nuc_ws.internal_uniform[ci],
                cells: Vec::new(),
            });
            next
        } else {
            nuc_ws.root_to_comp[root]
        };
        let comp = &mut nuc_ws.components[idx];
        comp.total_verts += ws.tc_size[ci] as usize;
        comp.cells.push((ws.tc_start[ci], ws.tc_size[ci]));
        if comp.cells.len() > 1 {
            comp.single_uniform_twin = false;
        }
    }

    nuc_ws
        .components
        .sort_unstable_by_key(|c| (c.total_verts, c.first_start));
    true
}

#[derive(Clone, Debug, Default)]
struct CompCanonicalGraph {
    offsets: Vec<usize>,
    entries: Vec<(u32, u32)>,
}

impl CompCanonicalGraph {
    fn build_from(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &Partition,
        cells: &[(u32, u32)],
        total_verts: usize,
    ) {
        self.offsets.clear();
        self.offsets.reserve(total_verts + 1);
        self.entries.clear();
        self.offsets.push(0);

        for &(start, len) in cells {
            let s = start as usize;
            let e = s + len as usize;
            for &v in &partition.lab[s..e] {
                let row_start = self.entries.len();
                for nb in view.graph.neighbors(NodeId(v)) {
                    let nb_cs = partition.cell_start[nb.node.index()] as usize;
                    if partition.cell_len[nb_cs] == 1 {
                        let col = partition.inv_lab[nb.node.index()];
                        let wt = view.edge_colors[nb.edge.index()];
                        self.entries.push((col, wt));
                    }
                }
                self.entries[row_start..].sort_unstable();
                self.offsets.push(self.entries.len());
            }
        }
    }

    fn compare_with(
        &self,
        view: EdgeColoredView<'_>,
        partition: &Partition,
        cells: &[(u32, u32)],
        row_buf: &mut Vec<(u32, u32)>,
    ) -> Ordering {
        let mut row = 0usize;
        for &(start, len) in cells {
            let s = start as usize;
            let e = s + len as usize;
            for &v in &partition.lab[s..e] {
                let best_row = &self.entries[self.offsets[row]..self.offsets[row + 1]];
                row_buf.clear();
                for nb in view.graph.neighbors(NodeId(v)) {
                    let nb_cs = partition.cell_start[nb.node.index()] as usize;
                    if partition.cell_len[nb_cs] == 1 {
                        row_buf.push((
                            partition.inv_lab[nb.node.index()],
                            view.edge_colors[nb.edge.index()],
                        ));
                    }
                }
                match row_buf.len().cmp(&best_row.len()) {
                    Ordering::Equal => {}
                    ord => return ord,
                }
                row_buf.sort_unstable();
                match row_buf.as_slice().cmp(best_row) {
                    Ordering::Equal => {}
                    ord => return ord,
                }
                row += 1;
            }
        }
        Ordering::Equal
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct ComponentSolver {
    target_num_cells: usize,
    total_verts: usize,
    emit_symmetry: bool,
    path_fixed: Vec<u32>,
    path_codes: Vec<u64>,
    path_tcs: Vec<u32>,
    cand_stack: Vec<u32>,
    explored_stack: Vec<u32>,
    first_lab: Vec<u32>,
    first_codes: Vec<u64>,
    first_tcs: Vec<u32>,
    first_path_fixed: Vec<u32>,
    pub(super) best_lab: Vec<u32>,
    best_codes: Vec<u64>,
    best_tcs: Vec<u32>,
    best_path_fixed: Vec<u32>,
    best_graph: CompCanonicalGraph,
    best_version: u64,
    has_first: bool,
    orbits: Vec<u32>,
    pub(super) generators: Vec<Vec<(u32, u32)>>,
    pub(super) group_factors: Vec<u64>,
    twin_buf: Vec<TwinCell>,
    row_buf: Vec<(u32, u32)>,
    perm_map: Vec<u32>,
    aut_mark: Vec<u32>,
    aut_wmark: Vec<u32>,
    aut_epoch: u32,
}

impl ComponentSolver {
    pub(super) fn solve(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &mut Partition,
        ws: &mut RefineWorkspace,
        cells: &[(u32, u32)],
        total_verts: usize,
        emit_symmetry: bool,
    ) {
        let n = view.node_count();
        self.total_verts = total_verts;
        self.target_num_cells = partition.num_cells + (total_verts - cells.len());
        self.emit_symmetry = emit_symmetry;
        self.has_first = false;
        self.best_version = 0;
        self.path_fixed.clear();
        self.path_codes.clear();
        self.path_codes.push(0);
        self.path_tcs.clear();
        self.cand_stack.clear();
        self.explored_stack.clear();
        self.first_codes.clear();
        self.first_tcs.clear();
        self.first_path_fixed.clear();
        self.best_codes.clear();
        self.best_tcs.clear();
        self.best_path_fixed.clear();
        self.generators.clear();
        self.group_factors.clear();

        if self.orbits.len() < n {
            let old_len = self.orbits.len();
            self.orbits.resize(n, 0);
            self.perm_map.resize(n, 0);
            for v in old_len..n {
                self.perm_map[v] = v as u32;
            }
            self.aut_mark.resize(n, 0);
            self.aut_wmark.resize(n, 0);
        }
        for &(start, len) in cells {
            let s = start as usize;
            for &v in &partition.lab[s..s + len as usize] {
                self.orbits[v as usize] = v;
            }
        }

        let entry_mark = partition.trail_mark();
        let _ = self.dfs(view, partition, ws, cells, 0, true, true, Ordering::Equal);
        partition.unwind(entry_mark);
    }

    #[inline]
    fn orbit_find(&mut self, mut x: u32) -> u32 {
        while self.orbits[x as usize] != x {
            let p = self.orbits[x as usize];
            let g = self.orbits[p as usize];
            self.orbits[x as usize] = g;
            x = p;
        }
        x
    }

    fn orbit_union(&mut self, a: u32, b: u32) {
        let ra = self.orbit_find(a);
        let rb = self.orbit_find(b);
        if ra != rb {
            if ra < rb {
                self.orbits[rb as usize] = ra;
            } else {
                self.orbits[ra as usize] = rb;
            }
        }
    }

    fn select_target_cell(&self, partition: &Partition, cells: &[(u32, u32)]) -> u32 {
        let mut best_start = u32::MAX;
        let mut best_len = usize::MAX;
        for &(start, len) in cells {
            let mut i = start as usize;
            let end = i + len as usize;
            while i < end {
                let clen = partition.cell_len[i] as usize;
                if clen > 1 && clen < best_len {
                    best_len = clen;
                    best_start = i as u32;
                    if clen == 2 {
                        return best_start;
                    }
                }
                i += clen;
            }
        }
        best_start
    }

    fn record_component_automorphism(
        &mut self,
        from_first: bool,
        partition: &Partition,
        cells: &[(u32, u32)],
    ) {
        let mut moves = Vec::new();
        let mut idx = 0usize;
        for &(start, len) in cells {
            let s = start as usize;
            for &dst in &partition.lab[s..s + len as usize] {
                let src = if from_first {
                    self.first_lab[idx]
                } else {
                    self.best_lab[idx]
                };
                idx += 1;
                if src != dst {
                    self.orbit_union(src, dst);
                    moves.push((src, dst));
                }
            }
        }
        if !moves.is_empty() && self.emit_symmetry {
            self.generators.push(moves);
        }
    }

    fn is_component_automorphism(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &Partition,
        cells: &[(u32, u32)],
    ) -> bool {
        let mut idx = 0usize;
        for &(start, len) in cells {
            let s = start as usize;
            for &dst in &partition.lab[s..s + len as usize] {
                let src = self.first_lab[idx];
                self.perm_map[src as usize] = dst;
                idx += 1;
            }
        }

        let mut ok = true;
        'check: for &src in &self.first_lab {
            let dst = self.perm_map[src as usize];
            let deg = view.graph.degree(NodeId(src));
            if view.graph.degree(NodeId(dst)) != deg {
                ok = false;
                break;
            }
            if src == dst && deg == 0 {
                continue;
            }
            self.aut_epoch = self.aut_epoch.wrapping_add(1);
            if self.aut_epoch == 0 {
                self.aut_mark.fill(0);
                self.aut_epoch = 1;
            }
            let epoch = self.aut_epoch;
            for nb in view.graph.neighbors(NodeId(src)) {
                let nb_cs = partition.cell_start[nb.node.index()] as usize;
                if partition.cell_len[nb_cs] == 1 {
                    let img = self.perm_map[nb.node.index()] as usize;
                    self.aut_mark[img] = epoch;
                    self.aut_wmark[img] = view.edge_colors[nb.edge.index()];
                }
            }
            for nb in view.graph.neighbors(NodeId(dst)) {
                let nb_cs = partition.cell_start[nb.node.index()] as usize;
                if partition.cell_len[nb_cs] == 1 {
                    let u = nb.node.index();
                    if self.aut_mark[u] != epoch
                        || self.aut_wmark[u] != view.edge_colors[nb.edge.index()]
                    {
                        ok = false;
                        break 'check;
                    }
                }
            }
        }

        for &src in &self.first_lab {
            self.perm_map[src as usize] = src;
        }
        ok
    }

    fn apply_local_twin_collapse(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &mut Partition,
        cells: &[(u32, u32)],
        emit_symmetry: bool,
    ) {
        collapse_twin_cells_in_ranges(view, partition, cells, &mut self.twin_buf);
        if !emit_symmetry || self.twin_buf.is_empty() {
            return;
        }
        let count = self.twin_buf.len();
        for idx in 0..count {
            let tc = self.twin_buf[idx];
            let first = partition.lab[tc.start];
            for pos in tc.start + 1..tc.start + tc.len {
                self.orbit_union(first, partition.lab[pos]);
            }
            for factor in 2..=tc.len as u64 {
                self.group_factors.push(factor);
            }
            let v0 = partition.lab[tc.start];
            let v1 = partition.lab[tc.start + 1];
            self.generators.push(vec![(v0, v1), (v1, v0)]);
            if tc.len >= 3 {
                let mut cyc = Vec::with_capacity(tc.len);
                for offset in 0..tc.len {
                    let src = partition.lab[tc.start + offset];
                    let dst = partition.lab[tc.start + ((offset + 1) % tc.len)];
                    cyc.push((src, dst));
                }
                self.generators.push(cyc);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn dfs(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &mut Partition,
        ws: &mut RefineWorkspace,
        cells: &[(u32, u32)],
        depth: usize,
        on_first_path: bool,
        eq_first: bool,
        cmp_best: Ordering,
    ) -> Option<usize> {
        if depth > 0 && partition.num_cells != self.target_num_cells {
            self.apply_local_twin_collapse(
                view,
                partition,
                cells,
                on_first_path && self.emit_symmetry,
            );
        }
        if partition.num_cells == self.target_num_cells {
            return self.process_leaf(view, partition, cells, eq_first, cmp_best);
        }

        let tc = self.select_target_cell(partition, cells);

        let mut eq_first = eq_first;
        let mut cmp_best = cmp_best;
        if self.has_first {
            if eq_first && (depth >= self.first_tcs.len() || tc != self.first_tcs[depth]) {
                eq_first = false;
            }
            if cmp_best == Ordering::Equal {
                cmp_best = if depth >= self.best_tcs.len() {
                    Ordering::Greater
                } else {
                    tc.cmp(&self.best_tcs[depth])
                };
            }
            if !eq_first && cmp_best == Ordering::Less {
                return None;
            }
        }

        self.path_tcs.push(tc);

        let cand_base = self.cand_stack.len();
        self.cand_stack
            .extend_from_slice(partition.cell_slice(tc as usize));
        let cand_len = self.cand_stack.len() - cand_base;
        let first_vertex = self.cand_stack[cand_base];
        let exp_base = self.explored_stack.len();
        let mut early_return: Option<usize> = None;

        for cand_idx in 0..cand_len {
            let v = self.cand_stack[cand_base + cand_idx];
            let is_first_child = cand_idx == 0;

            if !is_first_child && on_first_path {
                let rv = self.orbit_find(v);
                let mut already_in_orbit = false;
                for exp_idx in exp_base..self.explored_stack.len() {
                    let prev = self.explored_stack[exp_idx];
                    if self.orbit_find(prev) == rv {
                        already_in_orbit = true;
                        break;
                    }
                }
                if already_in_orbit {
                    continue;
                }
            }

            let best_ver_before = self.best_version;
            let mark = partition.trail_mark();
            partition.individualize(tc as usize, v);
            self.path_fixed.push(v);

            let active_single = [tc];
            let child_code = refine(view, partition, &active_single, ws);
            self.path_codes.push(child_code);
            let next_depth = depth + 1;

            let child_on_first = on_first_path && is_first_child;
            let mut child_eq_first = eq_first;
            let mut child_cmp_best = cmp_best;
            let mut prune = false;

            if self.has_first {
                if child_eq_first
                    && (next_depth >= self.first_codes.len()
                        || child_code != self.first_codes[next_depth])
                {
                    child_eq_first = false;
                }
                if child_cmp_best == Ordering::Equal {
                    debug_assert!(next_depth < self.best_codes.len());
                    child_cmp_best = child_code.cmp(&self.best_codes[next_depth]);
                }
                if !child_eq_first && child_cmp_best == Ordering::Less {
                    prune = true;
                }
            }

            let backjump = if prune {
                None
            } else {
                self.dfs(
                    view,
                    partition,
                    ws,
                    cells,
                    next_depth,
                    child_on_first,
                    child_eq_first,
                    child_cmp_best,
                )
            };

            self.path_codes.pop();
            self.path_fixed.pop();
            partition.unwind(mark);
            self.explored_stack.push(v);

            if self.best_version != best_ver_before {
                cmp_best = Ordering::Equal;
            }

            if let Some(target_depth) = backjump {
                if target_depth < depth {
                    early_return = Some(target_depth);
                    break;
                }
            }
        }

        self.explored_stack.truncate(exp_base);
        self.path_tcs.pop();

        if let Some(target_depth) = early_return {
            self.cand_stack.truncate(cand_base);
            return Some(target_depth);
        }

        if on_first_path && self.emit_symmetry {
            let rep = self.orbit_find(first_vertex);
            let mut orbit_size = 0u64;
            for idx in cand_base..cand_base + cand_len {
                let c = self.cand_stack[idx];
                if self.orbit_find(c) == rep {
                    orbit_size += 1;
                }
            }
            if orbit_size > 1 {
                self.group_factors.push(orbit_size);
            }
        }

        self.cand_stack.truncate(cand_base);
        None
    }

    fn copy_component_lab(partition: &Partition, cells: &[(u32, u32)], dst: &mut Vec<u32>) {
        dst.clear();
        for &(start, len) in cells {
            let s = start as usize;
            dst.extend_from_slice(&partition.lab[s..s + len as usize]);
        }
    }

    fn process_leaf(
        &mut self,
        view: EdgeColoredView<'_>,
        partition: &Partition,
        cells: &[(u32, u32)],
        eq_first: bool,
        cmp_best: Ordering,
    ) -> Option<usize> {
        if !self.has_first {
            self.has_first = true;
            Self::copy_component_lab(partition, cells, &mut self.first_lab);
            self.first_codes.clear();
            self.first_codes.extend_from_slice(&self.path_codes);
            self.first_tcs.clear();
            self.first_tcs.extend_from_slice(&self.path_tcs);
            self.first_path_fixed.clear();
            self.first_path_fixed.extend_from_slice(&self.path_fixed);

            self.best_lab.clear();
            self.best_lab.extend_from_slice(&self.first_lab);
            self.best_codes.clear();
            self.best_codes.extend_from_slice(&self.path_codes);
            self.best_tcs.clear();
            self.best_tcs.extend_from_slice(&self.path_tcs);
            self.best_path_fixed.clear();
            self.best_path_fixed.extend_from_slice(&self.path_fixed);
            self.best_graph
                .build_from(view, partition, cells, self.total_verts);
            self.best_version += 1;
            return None;
        }

        if eq_first
            && self.path_fixed.len() == self.first_path_fixed.len()
            && self.is_component_automorphism(view, partition, cells)
        {
            self.record_component_automorphism(true, partition, cells);
            return Some(common_prefix_len(&self.path_fixed, &self.first_path_fixed));
        }

        let leaf_ord = match cmp_best {
            Ordering::Less => return None,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => match self.path_fixed.len().cmp(&self.best_path_fixed.len()) {
                Ordering::Equal => {
                    self.best_graph
                        .compare_with(view, partition, cells, &mut self.row_buf)
                }
                ord => ord,
            },
        };

        match leaf_ord {
            Ordering::Less => None,
            Ordering::Equal => {
                self.record_component_automorphism(false, partition, cells);
                Some(common_prefix_len(&self.path_fixed, &self.best_path_fixed))
            }
            Ordering::Greater => {
                Self::copy_component_lab(partition, cells, &mut self.best_lab);
                self.best_codes.clear();
                self.best_codes.extend_from_slice(&self.path_codes);
                self.best_tcs.clear();
                self.best_tcs.extend_from_slice(&self.path_tcs);
                self.best_path_fixed.clear();
                self.best_path_fixed.extend_from_slice(&self.path_fixed);
                self.best_graph
                    .build_from(view, partition, cells, self.total_verts);
                self.best_version += 1;
                None
            }
        }
    }
}
