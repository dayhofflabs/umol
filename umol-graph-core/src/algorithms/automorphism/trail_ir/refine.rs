//! Hopcroft edge-colored partition refinement, invariant hashing, and target-cell selection.

use std::hash::{BuildHasher, Hasher};

use foldhash::fast::{FixedState, FoldHasher};

use super::partition::{EdgeColoredView, Partition};
use crate::id::NodeId;

const REFINEMENT_HASH_STATE: FixedState = FixedState::with_seed(0);

#[derive(Clone, Debug, Default)]
pub(super) struct RefineWorkspace {
    active_queue: Vec<u32>,
    in_active: Vec<u32>,
    active_epoch: u32,
    mark1: Vec<u32>,
    mark1_epoch: u32,
    mark2: Vec<u32>,
    mark2_epoch: u32,
    hits: Vec<u32>,
    cell_write: Vec<u32>,
    cell_hit_count: Vec<u32>,
    hit_cells: Vec<u32>,
    hit_verts: Vec<u32>,
    bucket0: Vec<u32>,
    bucket1: Vec<u32>,
    hit_edges: Vec<(u32, u32)>,
    pub(super) tc_start: Vec<u32>,
    pub(super) tc_size: Vec<u32>,
    pub(super) tc_nntcell: Vec<u32>,
    pub(super) tc_hits: Vec<u32>,
    pub(super) tc_first_wt: Vec<u32>,
    pub(super) tc_mixed: Vec<u8>,
}

impl RefineWorkspace {
    pub(super) fn ensure_capacity(&mut self, n: usize) {
        if self.in_active.len() < n {
            self.in_active.resize(n, 0);
            self.mark1.resize(n, 0);
            self.mark2.resize(n, 0);
            self.hits.resize(n, 0);
            self.cell_write.resize(n, 0);
            self.cell_hit_count.resize(n, 0);
            self.tc_nntcell.resize(n, 0);
        }
    }

    #[inline]
    fn reset_active(&mut self) {
        self.active_queue.clear();
        self.active_epoch = self.active_epoch.wrapping_add(1);
        if self.active_epoch == 0 {
            self.in_active.fill(0);
            self.active_epoch = 1;
        }
    }

    #[inline]
    fn push_active(&mut self, cell_start: u32) {
        let idx = cell_start as usize;
        if self.in_active[idx] != self.active_epoch {
            self.in_active[idx] = self.active_epoch;
            self.active_queue.push(cell_start);
        }
    }

    #[inline]
    fn is_active(&self, cell_start: u32) -> bool {
        self.in_active[cell_start as usize] == self.active_epoch
    }

    #[inline]
    fn clear_active_flag(&mut self, cell_start: u32) {
        self.in_active[cell_start as usize] = 0;
    }

    #[inline]
    fn next_mark1(&mut self) -> u32 {
        self.mark1_epoch = self.mark1_epoch.wrapping_add(1);
        if self.mark1_epoch == 0 {
            self.mark1.fill(0);
            self.mark1_epoch = 1;
        }
        self.mark1_epoch
    }

    #[inline]
    fn next_mark2(&mut self) -> u32 {
        self.mark2_epoch = self.mark2_epoch.wrapping_add(1);
        if self.mark2_epoch == 0 {
            self.mark2.fill(0);
            self.mark2_epoch = 1;
        }
        self.mark2_epoch
    }
}

#[derive(Clone, Copy)]
enum HitTargetSource {
    Bucket0,
    Bucket1,
    HitEdges(usize, usize),
}

#[inline]
fn sort_u32_asc(slice: &mut [u32]) {
    let n = slice.len();
    if n <= 1 {
        return;
    }
    if n <= 16 {
        for i in 1..n {
            let x = slice[i];
            let mut j = i;
            while j > 0 && slice[j - 1] > x {
                slice[j] = slice[j - 1];
                j -= 1;
            }
            slice[j] = x;
        }
    } else {
        slice.sort_unstable();
    }
}

pub(super) fn refine(
    view: EdgeColoredView<'_>,
    partition: &mut Partition,
    initial_active: &[u32],
    ws: &mut RefineWorkspace,
) -> u64 {
    let n = view.node_count();
    ws.ensure_capacity(n);
    ws.reset_active();

    for &cell in initial_active {
        ws.push_active(cell);
    }

    let mut hasher = REFINEMENT_HASH_STATE.build_hasher();
    hasher.write_usize(partition.num_cells);

    while !ws.active_queue.is_empty() && partition.num_cells < n {
        let qlen = ws.active_queue.len();
        let search_start = qlen.saturating_sub(10);
        let mut chosen_idx = qlen - 1;
        for idx in (search_start..qlen).rev() {
            if partition.cell_len[ws.active_queue[idx] as usize] == 1 {
                chosen_idx = idx;
                break;
            }
        }
        let isplit = ws.active_queue.swap_remove(chosen_idx);
        ws.clear_active_flag(isplit);

        let split_start = isplit as usize;
        let split_len = partition.cell_len[split_start] as usize;
        hasher.write_u32(isplit);

        if split_len == 1 {
            let splitv = partition.lab[split_start];
            let m1 = ws.next_mark1();
            ws.hit_cells.clear();

            if view.uniform_edges {
                for nb in view.graph.neighbors(NodeId(splitv)) {
                    let u = nb.node.index();
                    let cs = partition.cell_start[u] as usize;
                    let clen = partition.cell_len[cs] as usize;
                    if clen > 1 {
                        if ws.mark1[cs] != m1 {
                            ws.mark1[cs] = m1;
                            ws.hit_cells.push(cs as u32);
                            ws.cell_write[cs] = (cs + clen) as u32;
                        }
                        ws.cell_write[cs] -= 1;
                        let w = ws.cell_write[cs] as usize;
                        let pos = partition.inv_lab[u] as usize;
                        partition.swap_positions(pos, w);
                    }
                }
                sort_u32_asc(&mut ws.hit_cells);
                hasher.write_usize(ws.hit_cells.len());

                for idx in 0..ws.hit_cells.len() {
                    let v1 = ws.hit_cells[idx] as usize;
                    let len = partition.cell_len[v1] as usize;
                    let end = v1 + len;
                    let write = ws.cell_write[v1] as usize;
                    if write == v1 || write == end {
                        continue;
                    }

                    let left_len = (write - v1) as u32;
                    let right_len = (end - write) as u32;
                    let right_start = write as u32;

                    partition.cell_len[v1] = left_len;
                    partition.cell_len[write] = right_len;
                    for &u in &partition.lab[write..end] {
                        partition.cell_start[u as usize] = right_start;
                    }
                    partition.record_split(v1 as u32, len as u32, 1);

                    hasher.write_usize(v1);
                    hasher.write_u32(right_start);

                    if ws.is_active(v1 as u32) {
                        ws.push_active(right_start);
                    } else if left_len <= right_len {
                        ws.push_active(v1 as u32);
                    } else {
                        ws.push_active(right_start);
                    }
                }
            } else {
                for nb in view.graph.neighbors(NodeId(splitv)) {
                    let u = nb.node.index();
                    let cs = partition.cell_start[u] as usize;
                    let clen = partition.cell_len[cs] as usize;
                    if clen > 1 {
                        if ws.mark1[cs] != m1 {
                            ws.mark1[cs] = m1;
                            ws.hit_cells.push(cs as u32);
                            ws.cell_write[cs] = (cs + clen) as u32;
                        }
                        ws.cell_write[cs] -= 1;
                        let w = ws.cell_write[cs] as usize;
                        let pos = partition.inv_lab[u] as usize;
                        partition.swap_positions(pos, w);
                        ws.hits[u] = view.edge_colors[nb.edge.index()];
                    }
                }
                sort_u32_asc(&mut ws.hit_cells);
                hasher.write_usize(ws.hit_cells.len());

                for idx in 0..ws.hit_cells.len() {
                    let v1 = ws.hit_cells[idx] as usize;
                    split_cell_sparse_hits(partition, ws, v1, &mut hasher);
                }
            }
        } else if view.uniform_edges {
            let m1 = ws.next_mark1();
            let m2 = ws.next_mark2();
            ws.hit_cells.clear();
            ws.hit_verts.clear();
            for pos in split_start..split_start + split_len {
                let sv = partition.lab[pos];
                for nb in view.graph.neighbors(NodeId(sv)) {
                    let u = nb.node.index();
                    let cs = partition.cell_start[u] as usize;
                    if partition.cell_len[cs] > 1 {
                        if ws.mark1[cs] != m1 {
                            ws.mark1[cs] = m1;
                            ws.hit_cells.push(cs as u32);
                            ws.cell_hit_count[cs] = 0;
                        }
                        if ws.mark2[u] != m2 {
                            ws.mark2[u] = m2;
                            ws.hits[u] = 1;
                            ws.cell_hit_count[cs] += 1;
                            ws.hit_verts.push(u as u32);
                        } else {
                            ws.hits[u] += 1;
                        }
                    }
                }
            }

            if ws.hit_cells.is_empty() {
                continue;
            }
            finish_nonsingleton_round(partition, ws, &mut hasher);
        } else {
            ws.bucket0.clear();
            ws.bucket1.clear();
            ws.hit_edges.clear();
            let mut w0: Option<u32> = None;
            let mut w1: Option<u32> = None;
            let mut three_or_more = false;

            for pos in split_start..split_start + split_len {
                let sv = partition.lab[pos];
                for nb in view.graph.neighbors(NodeId(sv)) {
                    let u = nb.node.index();
                    let cs = partition.cell_start[u] as usize;
                    if partition.cell_len[cs] > 1 {
                        let wt = view.edge_colors[nb.edge.index()];
                        if three_or_more {
                            ws.hit_edges.push((wt, u as u32));
                        } else if let Some(c0) = w0 {
                            if wt == c0 {
                                ws.bucket0.push(u as u32);
                            } else if let Some(c1) = w1 {
                                if wt == c1 {
                                    ws.bucket1.push(u as u32);
                                } else {
                                    three_or_more = true;
                                    ws.hit_edges
                                        .reserve(ws.bucket0.len() + ws.bucket1.len() + 16);
                                    for &v in &ws.bucket0 {
                                        ws.hit_edges.push((c0, v));
                                    }
                                    for &v in &ws.bucket1 {
                                        ws.hit_edges.push((c1, v));
                                    }
                                    ws.hit_edges.push((wt, u as u32));
                                }
                            } else {
                                w1 = Some(wt);
                                ws.bucket1.push(u as u32);
                            }
                        } else {
                            w0 = Some(wt);
                            ws.bucket0.push(u as u32);
                        }
                    }
                }
            }

            let Some(c0) = w0 else {
                continue;
            };

            if !three_or_more {
                match w1 {
                    None => {
                        hasher.write_u32(c0);
                        refine_hit_targets(partition, ws, HitTargetSource::Bucket0, &mut hasher);
                    }
                    Some(c1) => {
                        if c0 < c1 {
                            hasher.write_u32(c0);
                            refine_hit_targets(
                                partition,
                                ws,
                                HitTargetSource::Bucket0,
                                &mut hasher,
                            );
                            if partition.num_cells < n {
                                hasher.write_u32(c1);
                                refine_hit_targets(
                                    partition,
                                    ws,
                                    HitTargetSource::Bucket1,
                                    &mut hasher,
                                );
                            }
                        } else {
                            hasher.write_u32(c1);
                            refine_hit_targets(
                                partition,
                                ws,
                                HitTargetSource::Bucket1,
                                &mut hasher,
                            );
                            if partition.num_cells < n {
                                hasher.write_u32(c0);
                                refine_hit_targets(
                                    partition,
                                    ws,
                                    HitTargetSource::Bucket0,
                                    &mut hasher,
                                );
                            }
                        }
                    }
                }
            } else {
                ws.hit_edges.sort_unstable_by_key(|&(wt, _)| wt);
                let total = ws.hit_edges.len();
                let mut gstart = 0usize;
                while gstart < total && partition.num_cells < n {
                    let cur_wt = ws.hit_edges[gstart].0;
                    let mut gend = gstart + 1;
                    while gend < total && ws.hit_edges[gend].0 == cur_wt {
                        gend += 1;
                    }
                    hasher.write_u32(cur_wt);
                    refine_hit_targets(
                        partition,
                        ws,
                        HitTargetSource::HitEdges(gstart, gend),
                        &mut hasher,
                    );
                    gstart = gend;
                }
            }
        }
    }

    hasher.write_usize(partition.num_cells);
    hasher.finish()
}

fn refine_hit_targets(
    partition: &mut Partition,
    ws: &mut RefineWorkspace,
    source: HitTargetSource,
    hasher: &mut FoldHasher,
) {
    let m1 = ws.next_mark1();
    let m2 = ws.next_mark2();
    ws.hit_cells.clear();
    ws.hit_verts.clear();

    let (start, end) = match source {
        HitTargetSource::Bucket0 => (0, ws.bucket0.len()),
        HitTargetSource::Bucket1 => (0, ws.bucket1.len()),
        HitTargetSource::HitEdges(start, end) => (start, end),
    };

    for idx in start..end {
        let u = match source {
            HitTargetSource::Bucket0 => ws.bucket0[idx] as usize,
            HitTargetSource::Bucket1 => ws.bucket1[idx] as usize,
            HitTargetSource::HitEdges(_, _) => ws.hit_edges[idx].1 as usize,
        };
        let cs = partition.cell_start[u] as usize;
        if partition.cell_len[cs] > 1 {
            if ws.mark1[cs] != m1 {
                ws.mark1[cs] = m1;
                ws.hit_cells.push(cs as u32);
                ws.cell_hit_count[cs] = 0;
            }
            if ws.mark2[u] != m2 {
                ws.mark2[u] = m2;
                ws.hits[u] = 1;
                ws.cell_hit_count[cs] += 1;
                ws.hit_verts.push(u as u32);
            } else {
                ws.hits[u] += 1;
            }
        }
    }

    if ws.hit_cells.is_empty() {
        return;
    }
    finish_nonsingleton_round(partition, ws, hasher);
}

fn finish_nonsingleton_round(
    partition: &mut Partition,
    ws: &mut RefineWorkspace,
    hasher: &mut FoldHasher,
) {
    let mut has_sparse = false;
    for &cs_u32 in &ws.hit_cells {
        let cs = cs_u32 as usize;
        let clen = partition.cell_len[cs];
        if ws.cell_hit_count[cs] < clen {
            ws.cell_write[cs] = cs_u32 + clen;
            has_sparse = true;
        } else {
            ws.cell_write[cs] = cs_u32;
        }
    }

    if has_sparse {
        for idx in 0..ws.hit_verts.len() {
            let u = ws.hit_verts[idx] as usize;
            let cs = partition.cell_start[u] as usize;
            if ws.cell_hit_count[cs] < partition.cell_len[cs] {
                ws.cell_write[cs] -= 1;
                let w = ws.cell_write[cs] as usize;
                let pos = partition.inv_lab[u] as usize;
                partition.swap_positions(pos, w);
            }
        }
    }

    sort_u32_asc(&mut ws.hit_cells);
    hasher.write_usize(ws.hit_cells.len());
    for idx in 0..ws.hit_cells.len() {
        let v1 = ws.hit_cells[idx] as usize;
        split_cell_sparse_hits(partition, ws, v1, hasher);
    }
}

fn split_cell_sparse_hits(
    partition: &mut Partition,
    ws: &mut RefineWorkspace,
    v1: usize,
    hasher: &mut FoldHasher,
) {
    let len = partition.cell_len[v1] as usize;
    let end = v1 + len;
    let write = ws.cell_write[v1] as usize;
    let has_unhit = write > v1;

    let first_val = ws.hits[partition.lab[write] as usize];
    let mut hit_all_equal = true;
    for &u in &partition.lab[write + 1..end] {
        if ws.hits[u as usize] != first_val {
            hit_all_equal = false;
            break;
        }
    }

    if !has_unhit && hit_all_equal {
        return;
    }

    let was_active = ws.is_active(v1 as u32);

    if has_unhit && hit_all_equal {
        let left_len = (write - v1) as u32;
        let right_len = (end - write) as u32;
        let right_start = write as u32;

        partition.cell_len[v1] = left_len;
        partition.cell_len[write] = right_len;
        for &u in &partition.lab[write..end] {
            partition.cell_start[u as usize] = right_start;
        }
        partition.record_split(v1 as u32, len as u32, 1);

        hasher.write_usize(v1);
        hasher.write_u32(left_len);
        hasher.write_u32(first_val);
        hasher.write_u32(right_len);

        if was_active {
            ws.push_active(right_start);
        } else if left_len <= right_len {
            ws.push_active(v1 as u32);
        } else {
            ws.push_active(right_start);
        }
        return;
    }

    let hit_slice = &mut partition.lab[write..end];
    if hit_slice.len() <= 16 {
        for i in 1..hit_slice.len() {
            let lv = hit_slice[i];
            let kv = ws.hits[lv as usize];
            let mut j = i;
            while j > 0 && ws.hits[hit_slice[j - 1] as usize] > kv {
                hit_slice[j] = hit_slice[j - 1];
                j -= 1;
            }
            hit_slice[j] = lv;
        }
    } else {
        hit_slice.sort_unstable_by_key(|&v| ws.hits[v as usize]);
    }
    partition.refresh_inv_range(write, end);

    let mut biggest_start = v1 as u32;
    let mut biggest_len = 0u32;
    let mut added_cells = 0u32;
    let mut first_sub = true;

    hasher.write_usize(v1);

    if has_unhit {
        let unhit_len = (write - v1) as u32;
        partition.cell_len[v1] = unhit_len;
        hasher.write_u32(unhit_len);
        biggest_len = unhit_len;
        biggest_start = v1 as u32;
        if !was_active {
            ws.push_active(v1 as u32);
        }
        first_sub = false;
    }

    let mut sub_start = write;
    while sub_start < end {
        let val = ws.hits[partition.lab[sub_start] as usize];
        let mut sub_end = sub_start + 1;
        while sub_end < end && ws.hits[partition.lab[sub_end] as usize] == val {
            sub_end += 1;
        }
        let sub_len = (sub_end - sub_start) as u32;
        let sub_start_u32 = sub_start as u32;

        partition.cell_len[sub_start] = sub_len;
        if !first_sub {
            for &u in &partition.lab[sub_start..sub_end] {
                partition.cell_start[u as usize] = sub_start_u32;
            }
            added_cells += 1;
        }
        hasher.write_u32(val);
        hasher.write_u32(sub_len);

        if sub_len > biggest_len {
            biggest_len = sub_len;
            biggest_start = sub_start_u32;
        }

        if was_active {
            if !first_sub {
                ws.push_active(sub_start_u32);
            }
        } else {
            ws.push_active(sub_start_u32);
        }

        first_sub = false;
        sub_start = sub_end;
    }

    if !was_active {
        ws.clear_active_flag(biggest_start);
        if let Some(pos) = ws.active_queue.iter().rposition(|&c| c == biggest_start) {
            ws.active_queue.swap_remove(pos);
        }
    }

    partition.record_split(v1 as u32, len as u32, added_cells);
}

pub(super) fn select_target_cell(partition: &Partition) -> usize {
    let n = partition.lab.len();
    let mut best_start = n;
    let mut best_len = usize::MAX;
    let mut i = 0usize;
    while i < n {
        let len = partition.cell_len[i] as usize;
        if len > 1 && len < best_len {
            best_len = len;
            best_start = i;
            if len == 2 {
                return best_start;
            }
        }
        i += len;
    }
    best_start
}
