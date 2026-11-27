//! Exact false-twin and monochromatic true-twin cell reduction for equitable partitions.

use super::partition::{EdgeColoredView, Partition};
use crate::id::NodeId;

#[derive(Clone, Copy, Debug)]
pub(super) struct TwinCell {
    pub(super) start: usize,
    pub(super) len: usize,
}

pub(super) fn collapse_twin_cells(
    view: EdgeColoredView<'_>,
    partition: &mut Partition,
    out: &mut Vec<TwinCell>,
) {
    out.clear();
    let n = view.node_count();
    let mut i = 0usize;
    while i < n {
        let len = partition.cell_len[i] as usize;
        if len > 1 && is_pure_twin_cell(view, partition, i, len) {
            out.push(TwinCell { start: i, len });
        }
        i += len;
    }

    for &tc in out.iter() {
        partition.collapse_cell_to_singletons(tc.start);
    }
}

pub(super) fn collapse_twin_cells_in_ranges(
    view: EdgeColoredView<'_>,
    partition: &mut Partition,
    ranges: &[(u32, u32)],
    out: &mut Vec<TwinCell>,
) {
    out.clear();
    for &(rstart, rlen) in ranges {
        let mut i = rstart as usize;
        let end = i + rlen as usize;
        while i < end {
            let len = partition.cell_len[i] as usize;
            if len > 1 && is_pure_twin_cell(view, partition, i, len) {
                out.push(TwinCell { start: i, len });
            }
            i += len;
        }
    }

    for &tc in out.iter() {
        partition.collapse_cell_to_singletons(tc.start);
    }
}

#[inline]
fn is_pure_twin_cell(
    view: EdgeColoredView<'_>,
    partition: &Partition,
    start: usize,
    len: usize,
) -> bool {
    let rep = partition.lab[start];
    let start_u32 = start as u32;
    let mut internal_count = 0usize;
    let mut internal_wt: Option<u32> = None;

    for nb in view.graph.neighbors(NodeId(rep)) {
        let cs = partition.cell_start[nb.node.index()];
        if partition.cell_len[cs as usize] == 1 {
            continue;
        }
        if cs != start_u32 {
            return false;
        }
        let wt = view.edge_colors[nb.edge.index()];
        match internal_wt {
            None => internal_wt = Some(wt),
            Some(first) if first != wt => return false,
            _ => {}
        }
        internal_count += 1;
    }

    internal_count == 0 || internal_count == len - 1
}
