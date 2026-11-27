//! Trail-based DFS individualization--refinement engine for edge-colored graphs.

use std::cmp::Ordering;
use std::mem;

use super::super::{AutomorphismGroup, AutomorphismGroupOrder};
use super::components::{find_non_uniform_components, ComponentSolver, NucWorkspace};
use super::partition::{common_prefix_len, EdgeColoredView, Partition};
use super::preprocess::{collapse_twin_cells, TwinCell};
use super::refine::{refine, select_target_cell, RefineWorkspace};
use crate::graph::Graph;
use crate::id::NodeId;

pub(super) const NUCR_MIN_CELLS: usize = 8;

#[derive(Clone, Debug, Default)]
struct CanonicalGraph {
    offsets: Vec<usize>,
    entries: Vec<(u32, u32)>,
}

impl CanonicalGraph {
    fn build_from(&mut self, view: EdgeColoredView<'_>, lab: &[u32], inv: &[u32]) {
        let n = view.node_count();
        self.offsets.clear();
        self.offsets.reserve(n + 1);
        self.entries.clear();
        self.entries.reserve(2 * view.graph.edge_count());

        self.offsets.push(0);
        for &v in lab {
            let row_start = self.entries.len();
            for nb in view.graph.neighbors(NodeId(v)) {
                let col = inv[nb.node.index()];
                let wt = view.edge_colors[nb.edge.index()];
                self.entries.push((col, wt));
            }
            self.entries[row_start..].sort_unstable();
            self.offsets.push(self.entries.len());
        }
    }

    fn compare_with(
        &self,
        view: EdgeColoredView<'_>,
        lab: &[u32],
        inv: &[u32],
        row_buf: &mut Vec<(u32, u32)>,
    ) -> Ordering {
        for (row, &v) in lab.iter().enumerate() {
            let best_row = &self.entries[self.offsets[row]..self.offsets[row + 1]];
            let deg = view.graph.degree(NodeId(v));
            match deg.cmp(&best_row.len()) {
                Ordering::Equal => {}
                ord => return ord,
            }
            row_buf.clear();
            for nb in view.graph.neighbors(NodeId(v)) {
                row_buf.push((inv[nb.node.index()], view.edge_colors[nb.edge.index()]));
            }
            row_buf.sort_unstable();
            match row_buf.as_slice().cmp(best_row) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        Ordering::Equal
    }
}

#[derive(Clone, Debug, Default)]
pub(in super::super) struct TrailIrWorkspace {
    partition: Partition,
    ws: RefineWorkspace,
    initial_active: Vec<u32>,
    twin_buf: Vec<TwinCell>,
    nuc_ws: NucWorkspace,
    comp_solver: ComponentSolver,
    path_fixed: Vec<u32>,
    path_codes: Vec<u64>,
    path_tcs: Vec<u32>,
    cand_stack: Vec<u32>,
    explored_stack: Vec<u32>,
    first_leaf: Vec<u32>,
    first_codes: Vec<u64>,
    first_tcs: Vec<u32>,
    first_path_fixed: Vec<u32>,
    best_leaf: Vec<u32>,
    best_inv: Vec<u32>,
    best_codes: Vec<u64>,
    best_tcs: Vec<u32>,
    best_path_fixed: Vec<u32>,
    best_graph: CanonicalGraph,
    best_version: u64,
    has_first: bool,
    orbits: Vec<u32>,
    orbit_count: usize,
    generators: Vec<Vec<u32>>,
    group_exact: Option<u128>,
    group_mantissa: f64,
    group_exponent: i32,
    row_buf: Vec<(u32, u32)>,
    perm_buf: Vec<u32>,
    aut_mark: Vec<u32>,
    aut_wmark: Vec<u32>,
    aut_epoch: u32,
}

impl TrailIrWorkspace {
    pub(in super::super) fn solve(
        &mut self,
        graph: &Graph,
        node_colors: &[u32],
        edge_colors: &[u32],
        get_canon: bool,
    ) -> (AutomorphismGroup, Vec<NodeId>) {
        let n = graph.node_count();
        if n == 0 {
            return (
                AutomorphismGroup {
                    orbits: Vec::new(),
                    node_count: 0,
                    orbit_count: 0,
                    group_order: AutomorphismGroupOrder::Exact(1),
                    generators: Vec::new(),
                },
                Vec::new(),
            );
        }

        let view = EdgeColoredView::new(graph, edge_colors);
        self.partition
            .reset_from_colors(node_colors, &mut self.initial_active);
        let root_code = refine(
            view,
            &mut self.partition,
            &self.initial_active,
            &mut self.ws,
        );

        if self.partition.is_discrete() {
            let labels = if get_canon {
                self.partition.lab.iter().copied().map(NodeId).collect()
            } else {
                Vec::new()
            };
            return (
                AutomorphismGroup {
                    orbits: (0..n as u32).map(NodeId).collect(),
                    node_count: n,
                    orbit_count: n,
                    group_order: AutomorphismGroupOrder::Exact(1),
                    generators: Vec::new(),
                },
                labels,
            );
        }

        self.has_first = false;
        self.best_version = 0;
        self.path_fixed.clear();
        self.path_codes.clear();
        self.path_codes.push(root_code);
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
        self.orbits.clear();
        self.orbits.extend(0..n as u32);
        self.orbit_count = n;
        self.group_exact = Some(1);
        self.group_mantissa = 1.0;
        self.group_exponent = 0;

        self.apply_twin_collapse(view, true);
        if !self.partition.is_discrete() {
            self.try_apply_nucr(view, true);
        }

        if self.partition.is_discrete() {
            for v in 0..n {
                let root = self.orbit_find(v as u32);
                self.orbits[v] = root;
            }
            let group_order = match self.group_exact {
                Some(exact) => AutomorphismGroupOrder::Exact(exact),
                None => AutomorphismGroupOrder::from_scientific(
                    self.group_mantissa,
                    self.group_exponent,
                ),
            };
            let labels = if get_canon {
                self.partition.lab.iter().copied().map(NodeId).collect()
            } else {
                Vec::new()
            };
            return (
                AutomorphismGroup {
                    orbits: self.orbits.iter().copied().map(NodeId).collect(),
                    node_count: n,
                    orbit_count: self.orbit_count,
                    group_order,
                    generators: self
                        .generators
                        .iter()
                        .map(|g| g.iter().copied().map(NodeId).collect())
                        .collect(),
                },
                labels,
            );
        }

        self.perm_buf.resize(n, 0);
        if self.aut_mark.len() < n {
            self.aut_mark.resize(n, 0);
            self.aut_wmark.resize(n, 0);
        }

        let _ = self.dfs(view, 0, true, true, Ordering::Equal, get_canon);

        for v in 0..n {
            let root = self.orbit_find(v as u32);
            self.orbits[v] = root;
        }

        let group_order = match self.group_exact {
            Some(exact) => AutomorphismGroupOrder::Exact(exact),
            None => {
                AutomorphismGroupOrder::from_scientific(self.group_mantissa, self.group_exponent)
            }
        };
        let labels = if get_canon {
            self.best_leaf.iter().copied().map(NodeId).collect()
        } else {
            Vec::new()
        };

        (
            AutomorphismGroup {
                orbits: self.orbits.iter().copied().map(NodeId).collect(),
                node_count: n,
                orbit_count: self.orbit_count,
                group_order,
                generators: self
                    .generators
                    .iter()
                    .map(|g| g.iter().copied().map(NodeId).collect())
                    .collect(),
            },
            labels,
        )
    }

    fn multiply_group_order(&mut self, factor: u64) {
        if factor <= 1 {
            return;
        }
        if let Some(exact) = self.group_exact {
            self.group_exact = exact.checked_mul(factor as u128);
        }
        self.group_mantissa *= factor as f64;
        while self.group_mantissa >= 10.0 {
            self.group_mantissa /= 10.0;
            self.group_exponent += 1;
        }
    }

    #[inline]
    fn orbit_find(&mut self, mut x: u32) -> u32 {
        while self.orbits[x as usize] != x {
            let parent = self.orbits[x as usize];
            let grand = self.orbits[parent as usize];
            self.orbits[x as usize] = grand;
            x = parent;
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
            self.orbit_count -= 1;
        }
    }

    fn record_perm_buf_automorphism(&mut self) {
        let mut nontrivial = false;
        let n = self.perm_buf.len();
        for src in 0..n {
            let dst = self.perm_buf[src];
            if src as u32 != dst {
                nontrivial = true;
                self.orbit_union(src as u32, dst);
            }
        }
        if nontrivial {
            self.generators.push(self.perm_buf.clone());
        }
    }

    fn is_perm_buf_automorphism(&mut self, view: EdgeColoredView<'_>) -> bool {
        let n = view.node_count();
        for i in 0..n {
            let pi = self.perm_buf[i] as usize;
            let deg = view.graph.degree(NodeId(i as u32));
            if view.graph.degree(NodeId(pi as u32)) != deg {
                return false;
            }
            if i == pi && deg == 0 {
                continue;
            }
            self.aut_epoch = self.aut_epoch.wrapping_add(1);
            if self.aut_epoch == 0 {
                self.aut_mark.fill(0);
                self.aut_epoch = 1;
            }
            let epoch = self.aut_epoch;
            for nb in view.graph.neighbors(NodeId(i as u32)) {
                let img = self.perm_buf[nb.node.index()] as usize;
                self.aut_mark[img] = epoch;
                self.aut_wmark[img] = view.edge_colors[nb.edge.index()];
            }
            for nb in view.graph.neighbors(NodeId(pi as u32)) {
                let u = nb.node.index();
                if self.aut_mark[u] != epoch
                    || self.aut_wmark[u] != view.edge_colors[nb.edge.index()]
                {
                    return false;
                }
            }
        }
        true
    }

    fn emit_single_twin_cell_symmetry(&mut self, start: usize, len: usize) {
        let n = self.partition.lab.len();
        let first = self.partition.lab[start];
        for pos in start + 1..start + len {
            let v = self.partition.lab[pos];
            self.orbit_union(first, v);
        }
        for factor in 2..=len as u64 {
            self.multiply_group_order(factor);
        }
        let mut gen1: Vec<u32> = (0..n as u32).collect();
        let v0 = self.partition.lab[start];
        let v1 = self.partition.lab[start + 1];
        gen1[v0 as usize] = v1;
        gen1[v1 as usize] = v0;
        self.generators.push(gen1);

        if len >= 3 {
            let mut gen2: Vec<u32> = (0..n as u32).collect();
            for offset in 0..len {
                let src = self.partition.lab[start + offset];
                let dst = self.partition.lab[start + ((offset + 1) % len)];
                gen2[src as usize] = dst;
            }
            self.generators.push(gen2);
        }
    }

    fn apply_twin_collapse(&mut self, view: EdgeColoredView<'_>, emit_symmetry: bool) {
        collapse_twin_cells(view, &mut self.partition, &mut self.twin_buf);
        if !emit_symmetry || self.twin_buf.is_empty() {
            return;
        }
        let count = self.twin_buf.len();
        for idx in 0..count {
            let tc = self.twin_buf[idx];
            self.emit_single_twin_cell_symmetry(tc.start, tc.len);
        }
    }

    fn try_apply_nucr(&mut self, view: EdgeColoredView<'_>, emit_symmetry: bool) -> bool {
        if !find_non_uniform_components(
            view,
            &self.partition,
            &mut self.ws,
            &mut self.nuc_ws,
            NUCR_MIN_CELLS,
        ) {
            return false;
        }

        let n = view.node_count();
        let components = mem::take(&mut self.nuc_ws.components);
        let last_idx = components.len() - 1;

        for (idx, comp) in components.iter().enumerate() {
            if comp.single_uniform_twin {
                let (start, len) = comp.cells[0];
                let s = start as usize;
                let l = len as usize;
                self.partition.collapse_cell_to_singletons(s);
                if emit_symmetry {
                    self.emit_single_twin_cell_symmetry(s, l);
                }
                continue;
            }
            if idx == last_idx {
                break;
            }

            self.comp_solver.solve(
                view,
                &mut self.partition,
                &mut self.ws,
                &comp.cells,
                comp.total_verts,
                emit_symmetry,
            );

            let mut cursor = 0usize;
            for &(start, len) in &comp.cells {
                let s = start as usize;
                let l = len as usize;
                self.partition.lab[s..s + l]
                    .copy_from_slice(&self.comp_solver.best_lab[cursor..cursor + l]);
                self.partition.refresh_inv_range(s, s + l);
                self.partition.collapse_cell_to_singletons(s);
                cursor += l;
            }

            if emit_symmetry {
                let factor_count = self.comp_solver.group_factors.len();
                for fi in 0..factor_count {
                    let f = self.comp_solver.group_factors[fi];
                    self.multiply_group_order(f);
                }
                let gen_count = self.comp_solver.generators.len();
                for gi in 0..gen_count {
                    let mut perm: Vec<u32> = (0..n as u32).collect();
                    let move_count = self.comp_solver.generators[gi].len();
                    for mi in 0..move_count {
                        let (src, dst) = self.comp_solver.generators[gi][mi];
                        perm[src as usize] = dst;
                        self.orbit_union(src, dst);
                    }
                    self.generators.push(perm);
                }
            }
        }

        self.nuc_ws.components = components;
        true
    }

    fn dfs(
        &mut self,
        view: EdgeColoredView<'_>,
        depth: usize,
        on_first_path: bool,
        eq_first: bool,
        cmp_best: Ordering,
        get_canon: bool,
    ) -> Option<usize> {
        if depth > 0 && !self.partition.is_discrete() {
            self.apply_twin_collapse(view, on_first_path);
            if !self.partition.is_discrete() {
                self.try_apply_nucr(view, on_first_path);
            }
        }
        if self.partition.is_discrete() {
            return self.process_leaf(view, eq_first, cmp_best, get_canon);
        }

        let tc = select_target_cell(&self.partition) as u32;

        let mut eq_first = eq_first;
        let mut cmp_best = cmp_best;
        if self.has_first {
            if eq_first && (depth >= self.first_tcs.len() || tc != self.first_tcs[depth]) {
                eq_first = false;
            }
            if !get_canon {
                if !eq_first {
                    return None;
                }
            } else {
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
        }

        self.path_tcs.push(tc);

        let cand_base = self.cand_stack.len();
        self.cand_stack
            .extend_from_slice(self.partition.cell_slice(tc as usize));
        let cand_len = self.cand_stack.len() - cand_base;
        let first_vertex = self.cand_stack[cand_base];
        let exp_base = self.explored_stack.len();
        let mut early_return: Option<usize> = None;

        for cand_idx in 0..cand_len {
            let v = self.cand_stack[cand_base + cand_idx];
            let is_first_child = cand_idx == 0;

            if !is_first_child {
                if on_first_path {
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
                } else {
                    let explored_slice = &self.explored_stack[exp_base..];
                    let pruned = self.generators.iter().any(|gen_perm| {
                        self.path_fixed[..depth]
                            .iter()
                            .all(|&fixed| gen_perm[fixed as usize] == fixed)
                            && explored_slice.iter().any(|&prev| {
                                gen_perm[prev as usize] == v || gen_perm[v as usize] == prev
                            })
                    });
                    if pruned {
                        continue;
                    }
                }
            }

            let best_ver_before = self.best_version;
            let mark = self.partition.trail_mark();
            self.partition.individualize(tc as usize, v);
            self.path_fixed.push(v);

            let active_single = [tc];
            let child_code = refine(view, &mut self.partition, &active_single, &mut self.ws);
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
                if !get_canon {
                    if !child_eq_first {
                        prune = true;
                    }
                } else {
                    if child_cmp_best == Ordering::Equal {
                        debug_assert!(next_depth < self.best_codes.len());
                        child_cmp_best = child_code.cmp(&self.best_codes[next_depth]);
                    }
                    if !child_eq_first && child_cmp_best == Ordering::Less {
                        prune = true;
                    }
                }
            }

            let backjump = if prune {
                None
            } else {
                self.dfs(
                    view,
                    next_depth,
                    child_on_first,
                    child_eq_first,
                    child_cmp_best,
                    get_canon,
                )
            };

            self.path_codes.pop();
            self.path_fixed.pop();
            self.partition.unwind(mark);
            self.explored_stack.push(v);

            if get_canon && self.best_version != best_ver_before {
                debug_assert_eq!(&self.path_codes[..=depth], &self.best_codes[..=depth]);
                debug_assert_eq!(&self.path_tcs[..=depth], &self.best_tcs[..=depth]);
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

        if on_first_path {
            let rep = self.orbit_find(first_vertex);
            let mut orbit_size = 0u64;
            for idx in cand_base..cand_base + cand_len {
                let c = self.cand_stack[idx];
                if self.orbit_find(c) == rep {
                    orbit_size += 1;
                }
            }
            self.multiply_group_order(orbit_size);
        }

        self.cand_stack.truncate(cand_base);
        None
    }

    fn process_leaf(
        &mut self,
        view: EdgeColoredView<'_>,
        eq_first: bool,
        cmp_best: Ordering,
        get_canon: bool,
    ) -> Option<usize> {
        let n = view.node_count();

        if !self.has_first {
            self.has_first = true;
            self.first_leaf.clear();
            self.first_leaf.extend_from_slice(&self.partition.lab);
            self.first_codes.clear();
            self.first_codes.extend_from_slice(&self.path_codes);
            self.first_tcs.clear();
            self.first_tcs.extend_from_slice(&self.path_tcs);
            self.first_path_fixed.clear();
            self.first_path_fixed.extend_from_slice(&self.path_fixed);

            if get_canon {
                self.best_leaf.clear();
                self.best_leaf.extend_from_slice(&self.partition.lab);
                self.best_inv.clear();
                self.best_inv.extend_from_slice(&self.partition.inv_lab);
                self.best_codes.clear();
                self.best_codes.extend_from_slice(&self.path_codes);
                self.best_tcs.clear();
                self.best_tcs.extend_from_slice(&self.path_tcs);
                self.best_path_fixed.clear();
                self.best_path_fixed.extend_from_slice(&self.path_fixed);
                self.best_graph
                    .build_from(view, &self.best_leaf, &self.best_inv);
                self.best_version += 1;
            }
            return None;
        }

        if eq_first && self.path_fixed.len() == self.first_path_fixed.len() {
            for r in 0..n {
                self.perm_buf[self.first_leaf[r] as usize] = self.partition.lab[r];
            }
            if self.is_perm_buf_automorphism(view) {
                self.record_perm_buf_automorphism();
                return Some(common_prefix_len(&self.path_fixed, &self.first_path_fixed));
            }
        }

        if !get_canon {
            return None;
        }

        let leaf_ord = match cmp_best {
            Ordering::Less => return None,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => match self.path_fixed.len().cmp(&self.best_path_fixed.len()) {
                Ordering::Equal => self.best_graph.compare_with(
                    view,
                    &self.partition.lab,
                    &self.partition.inv_lab,
                    &mut self.row_buf,
                ),
                ord => ord,
            },
        };

        match leaf_ord {
            Ordering::Less => None,
            Ordering::Equal => {
                for r in 0..n {
                    self.perm_buf[self.best_leaf[r] as usize] = self.partition.lab[r];
                }
                debug_assert!(self.is_perm_buf_automorphism(view));
                self.record_perm_buf_automorphism();
                Some(common_prefix_len(&self.path_fixed, &self.best_path_fixed))
            }
            Ordering::Greater => {
                self.best_leaf.clear();
                self.best_leaf.extend_from_slice(&self.partition.lab);
                self.best_inv.clear();
                self.best_inv.extend_from_slice(&self.partition.inv_lab);
                self.best_codes.clear();
                self.best_codes.extend_from_slice(&self.path_codes);
                self.best_tcs.clear();
                self.best_tcs.extend_from_slice(&self.path_tcs);
                self.best_path_fixed.clear();
                self.best_path_fixed.extend_from_slice(&self.path_fixed);
                self.best_graph
                    .build_from(view, &self.best_leaf, &self.best_inv);
                self.best_version += 1;
                None
            }
        }
    }
}
