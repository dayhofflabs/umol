//! Shared bit-parallel primitives for `BitsetRi` and `Glasgow`.

use std::cmp::Reverse;
use std::mem;
use std::ops::{BitAnd, BitOr, Not, Shl};

use super::predicate::SubgraphIsomorphismPredicate;
use super::prepared::{PreparedSubgraphQuery, SubgraphIsomorphismScratch};
use super::symmetry::TargetSymmetry;
use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

/// Maximum target node count for lazily allocating a dense `n2 * words` wide adjacency matrix in
/// scratch; larger sparse targets use `Graph::neighbors` / `Graph::find_edge` directly to avoid
/// quadratic memory.
pub(super) const WIDE_LAZY_ADJ_LIMIT: usize = 4096;

pub(super) trait WordMask:
    Copy
    + Eq
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + Not<Output = Self>
    + Shl<usize, Output = Self>
{
    const ZERO: Self;

    fn bit(index: usize) -> Self;
    fn all_below(count: usize) -> Self;
    fn range_mask(min_v: usize, max_v: usize) -> Self;
    fn is_zero(self) -> bool;
    fn count_ones(self) -> u32;
    fn trailing_zeros(self) -> u32;
    fn clear_lowest(self) -> Self;
    fn test(self, index: usize) -> bool;
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self;
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self;
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self>;
    fn take_trail(scratch: &mut SubgraphIsomorphismScratch) -> Vec<Self>;
    fn restore_trail(scratch: &mut SubgraphIsomorphismScratch, trail: Vec<Self>);
}

impl WordMask for u64 {
    const ZERO: Self = 0;

    #[inline(always)]
    fn bit(index: usize) -> Self {
        1u64 << index
    }

    #[inline(always)]
    fn all_below(count: usize) -> Self {
        if count == 0 {
            0
        } else if count >= 64 {
            u64::MAX
        } else {
            (1u64 << count) - 1
        }
    }

    #[inline(always)]
    fn range_mask(min_v: usize, max_v: usize) -> Self {
        if min_v >= max_v || min_v >= 64 {
            0
        } else {
            Self::all_below(max_v.min(64)) & (u64::MAX << min_v)
        }
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == 0
    }

    #[inline(always)]
    fn count_ones(self) -> u32 {
        u64::count_ones(self)
    }

    #[inline(always)]
    fn trailing_zeros(self) -> u32 {
        u64::trailing_zeros(self)
    }

    #[inline(always)]
    fn clear_lowest(self) -> Self {
        self & self.wrapping_sub(1)
    }

    #[inline(always)]
    fn test(self, index: usize) -> bool {
        ((self >> index) & 1) != 0
    }

    #[inline(always)]
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self {
        target_sym.orbit_reps_u64()
    }

    #[inline(always)]
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self {
        pred.node_candidates_u64(query_node, eligible)
    }

    #[inline(always)]
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self> {
        pred.edge_neighbors_u64(query_edge, target_node)
    }

    #[inline]
    fn take_trail(scratch: &mut SubgraphIsomorphismScratch) -> Vec<Self> {
        mem::take(&mut scratch.trail_u64)
    }

    #[inline]
    fn restore_trail(scratch: &mut SubgraphIsomorphismScratch, trail: Vec<Self>) {
        scratch.trail_u64 = trail;
    }
}

impl WordMask for u128 {
    const ZERO: Self = 0;

    #[inline(always)]
    fn bit(index: usize) -> Self {
        1u128 << index
    }

    #[inline(always)]
    fn all_below(count: usize) -> Self {
        if count == 0 {
            0
        } else if count >= 128 {
            u128::MAX
        } else {
            (1u128 << count) - 1
        }
    }

    #[inline(always)]
    fn range_mask(min_v: usize, max_v: usize) -> Self {
        if min_v >= max_v || min_v >= 128 {
            0
        } else {
            Self::all_below(max_v.min(128)) & (u128::MAX << min_v)
        }
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == 0
    }

    #[inline(always)]
    fn count_ones(self) -> u32 {
        u128::count_ones(self)
    }

    #[inline(always)]
    fn trailing_zeros(self) -> u32 {
        u128::trailing_zeros(self)
    }

    #[inline(always)]
    fn clear_lowest(self) -> Self {
        self & self.wrapping_sub(1)
    }

    #[inline(always)]
    fn test(self, index: usize) -> bool {
        ((self >> index) & 1) != 0
    }

    #[inline(always)]
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self {
        target_sym.orbit_reps_u128()
    }

    #[inline(always)]
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self {
        pred.node_candidates_u128(query_node, eligible)
    }

    #[inline(always)]
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self> {
        pred.edge_neighbors_u128(query_edge, target_node)
    }

    #[inline]
    fn take_trail(scratch: &mut SubgraphIsomorphismScratch) -> Vec<Self> {
        mem::take(&mut scratch.trail_u128)
    }

    #[inline]
    fn restore_trail(scratch: &mut SubgraphIsomorphismScratch, trail: Vec<Self>) {
        scratch.trail_u128 = trail;
    }
}

pub(super) trait WordAdj<W: WordMask> {
    fn row(&mut self, u: usize) -> W;
}

pub(super) struct PrecomputedWordAdj<'a, W>(pub(super) &'a [W]);

impl<W: WordMask> WordAdj<W> for PrecomputedWordAdj<'_, W> {
    #[inline(always)]
    fn row(&mut self, u: usize) -> W {
        self.0[u]
    }
}

pub(super) struct LazyWordAdj<'a, W> {
    pub(super) target: &'a Graph,
    pub(super) adj: &'a mut [W],
    pub(super) built: W,
}

impl<W: WordMask> WordAdj<W> for LazyWordAdj<'_, W> {
    #[inline(always)]
    fn row(&mut self, u: usize) -> W {
        let bit = W::bit(u);
        if (self.built & bit).is_zero() {
            let mut mask = W::ZERO;
            for nb in self.target.neighbors(NodeId(u as u32)) {
                mask = mask | W::bit(nb.node.index());
            }
            self.adj[u] = mask;
            self.built = self.built | bit;
            mask
        } else {
            self.adj[u]
        }
    }
}

pub(super) enum WideAdj<'a> {
    Precomputed(&'a [u64]),
    Lazy {
        adj: &'a mut [u64],
        built: &'a mut [u64],
    },
    None,
}

impl WideAdj<'_> {
    #[inline(always)]
    pub(super) fn ensure_row(&mut self, target: &Graph, u: usize, words: usize) {
        if let Self::Lazy { adj, built } = self {
            let bw = u >> 6;
            let bm = 1u64 << (u & 63);
            if (built[bw] & bm) == 0 {
                let base = u * words;
                adj[base..base + words].fill(0);
                for nb in target.neighbors(NodeId(u as u32)) {
                    let v = nb.node.index();
                    adj[base + (v >> 6)] |= 1u64 << (v & 63);
                }
                built[bw] |= bm;
            }
        }
    }

    #[inline(always)]
    pub(super) fn has_edge_bit(&self, u: usize, vw: usize, vb: u64, words: usize) -> Option<bool> {
        match self {
            Self::Precomputed(adj) => Some((adj[u * words + vw] & vb) != 0),
            Self::Lazy { adj, .. } => Some((adj[u * words + vw] & vb) != 0),
            Self::None => None,
        }
    }

    /// Return the bitset row of target neighbors of `u`, materializing into `fallback_row` when
    /// `Self::None` is active.
    #[inline]
    pub(super) fn row_into<'r>(
        &'r mut self,
        target: &Graph,
        u: usize,
        words: usize,
        fallback_row: &'r mut [u64],
    ) -> &'r [u64] {
        match self {
            Self::Precomputed(adj) => &adj[u * words..(u + 1) * words],
            Self::Lazy { adj, built } => {
                let bw = u >> 6;
                let bm = 1u64 << (u & 63);
                let base = u * words;
                if (built[bw] & bm) == 0 {
                    adj[base..base + words].fill(0);
                    for nb in target.neighbors(NodeId(u as u32)) {
                        let v = nb.node.index();
                        adj[base + (v >> 6)] |= 1u64 << (v & 63);
                    }
                    built[bw] |= bm;
                }
                &adj[base..base + words]
            }
            Self::None => {
                fallback_row.fill(0);
                for nb in target.neighbors(NodeId(u as u32)) {
                    let v = nb.node.index();
                    fallback_row[v >> 6] |= 1u64 << (v & 63);
                }
                fallback_row
            }
        }
    }
}

/// O(n1 + n2) necessary degree-sequence domination filter for direct calls.
#[inline]
pub(super) fn degree_sequence_ok(query: &Graph, target: &Graph) -> bool {
    let mut q_hist = [0u32; 16];
    let mut t_hist = [0u32; 16];
    for q in query.node_ids() {
        q_hist[query.degree(q).min(15)] += 1;
    }
    for t in target.node_ids() {
        t_hist[target.degree(t).min(15)] += 1;
    }
    let mut q_ge = 0u32;
    let mut t_ge = 0u32;
    for d in (1..16).rev() {
        q_ge += q_hist[d];
        t_ge += t_hist[d];
        if t_ge < q_ge {
            return false;
        }
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn init_word_candidates_by_class<W: WordMask, P: SubgraphIsomorphismPredicate>(
    n1: usize,
    q_deg: &[usize],
    q_by_deg: &[usize],
    anchor: Option<(NodeId, NodeId)>,
    anchor_clear: W,
    root_orbit_filter: Option<(usize, W)>,
    mut eligible_for_deg: impl FnMut(usize) -> W,
    pred: &mut P,
    cand: &mut [W],
    mut cand_cnt: Option<&mut [usize]>,
    class_min_deg: &mut [usize],
    class_count: &mut [usize],
    class_computed: &mut [bool],
    class_masks: &mut [W],
) -> bool {
    for &q in &q_by_deg[..n1] {
        if anchor.is_some_and(|(qa, _)| q == qa.index()) {
            continue;
        }
        let c = pred.query_node_class(NodeId(q as u32));
        debug_assert!(c < class_min_deg.len());
        class_min_deg[c] = class_min_deg[c].min(q_deg[q]);
        class_count[c] += 1;
    }

    for &q in &q_by_deg[..n1] {
        if let Some((qa, ta)) = anchor {
            if q == qa.index() {
                cand[q] = W::bit(ta.index());
                if let Some(cnts) = cand_cnt.as_deref_mut() {
                    cnts[q] = 1;
                }
                continue;
            }
        }
        let c = pred.query_node_class(NodeId(q as u32));
        if !class_computed[c] {
            let mut eligible = eligible_for_deg(class_min_deg[c]) & !anchor_clear;
            if class_count[c] == 1 {
                if let Some((p0, orbit_reps)) = root_orbit_filter {
                    if q == p0 {
                        eligible = eligible & orbit_reps;
                    }
                }
            }
            class_masks[c] = if eligible.is_zero() {
                W::ZERO
            } else {
                W::node_candidates(pred, NodeId(q as u32), eligible)
            };
            class_computed[c] = true;
        }
        let mut mask = class_masks[c] & eligible_for_deg(q_deg[q]);
        if let Some((p0, orbit_reps)) = root_orbit_filter {
            if q == p0 {
                mask = mask & orbit_reps;
            }
        }
        if mask.is_zero() {
            return false;
        }
        cand[q] = mask;
        if let Some(cnts) = cand_cnt.as_deref_mut() {
            cnts[q] = mask.count_ones() as usize;
        }
    }
    true
}

/// Populate per-node candidate bitmasks for `n2 <= 128` using descending-degree fail-fast and
/// optional query-node class deduplication. Returns `false` if any query node has an empty domain.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_word_candidates<W: WordMask, const N: usize, P: SubgraphIsomorphismPredicate>(
    n1: usize,
    q_deg: &[usize],
    q_by_deg: &[usize],
    anchor: Option<(NodeId, NodeId)>,
    root_orbit_filter: Option<(usize, W)>,
    mut eligible_for_deg: impl FnMut(usize) -> W,
    pred: &mut P,
    cand: &mut [W],
    mut cand_cnt: Option<&mut [usize]>,
) -> bool {
    let anchor_clear = anchor.map_or(W::ZERO, |(_, ta)| W::bit(ta.index()));

    if let Some(num_classes) = pred.query_node_class_count() {
        if num_classes <= N {
            let mut class_min_deg = [usize::MAX; N];
            let mut class_count = [0usize; N];
            let mut class_computed = [false; N];
            let mut class_masks = [W::ZERO; N];
            return init_word_candidates_by_class::<W, P>(
                n1,
                q_deg,
                q_by_deg,
                anchor,
                anchor_clear,
                root_orbit_filter,
                eligible_for_deg,
                pred,
                cand,
                cand_cnt,
                &mut class_min_deg[..num_classes],
                &mut class_count[..num_classes],
                &mut class_computed[..num_classes],
                &mut class_masks[..num_classes],
            );
        }
        let mut class_min_deg = vec![usize::MAX; num_classes];
        let mut class_count = vec![0usize; num_classes];
        let mut class_computed = vec![false; num_classes];
        let mut class_masks = vec![W::ZERO; num_classes];
        return init_word_candidates_by_class::<W, P>(
            n1,
            q_deg,
            q_by_deg,
            anchor,
            anchor_clear,
            root_orbit_filter,
            eligible_for_deg,
            pred,
            cand,
            cand_cnt,
            &mut class_min_deg,
            &mut class_count,
            &mut class_computed,
            &mut class_masks,
        );
    }

    for &q in &q_by_deg[..n1] {
        if let Some((qa, ta)) = anchor {
            if q == qa.index() {
                cand[q] = W::bit(ta.index());
                if let Some(cnts) = cand_cnt.as_deref_mut() {
                    cnts[q] = 1;
                }
                continue;
            }
        }
        let mut eligible = eligible_for_deg(q_deg[q]) & !anchor_clear;
        if let Some((p0, orbit_reps)) = root_orbit_filter {
            if q == p0 {
                eligible = eligible & orbit_reps;
            }
        }
        if eligible.is_zero() {
            return false;
        }
        let mask = W::node_candidates(pred, NodeId(q as u32), eligible);
        if mask.is_zero() {
            return false;
        }
        cand[q] = mask;
        if let Some(cnts) = cand_cnt.as_deref_mut() {
            cnts[q] = mask.count_ones() as usize;
        }
    }
    true
}

/// Populate `scratch.cand_wide` (`n1 * words`) and `scratch.cand_counts` (`n1`) for `n2 > 128`
/// using descending-degree fail-fast and optional query-node class deduplication.
#[allow(clippy::too_many_arguments)]
pub(super) fn init_wide_candidates<P: SubgraphIsomorphismPredicate>(
    target: &Graph,
    q_deg: &[usize],
    prepared_query: Option<&PreparedSubgraphQuery>,
    prepared_deg_ge: Option<(&[u64], usize)>,
    anchor: Option<(NodeId, NodeId)>,
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
) -> bool {
    let n1 = q_deg.len();
    let n2 = target.node_count();
    let words = prepared_deg_ge.map_or_else(|| n2.div_ceil(64), |(_, w)| w);

    let mut q_by_deg: Vec<usize> = (0..n1).collect();
    q_by_deg.sort_unstable_by_key(|&q| Reverse(q_deg[q]));

    let root_orbit_filter = if anchor.is_none() {
        target_sym.filter(|sym| !sym.is_trivial()).map(|sym| {
            let p0 = prepared_query
                .expect("target symmetry is paired with prepared query")
                .priority_order()[0] as usize;
            (p0, sym.orbit_reps_words())
        })
    } else {
        None
    };

    scratch.cand_wide.clear();
    scratch.cand_wide.resize(n1 * words, 0);
    scratch.cand_counts.clear();
    scratch.cand_counts.resize(n1, 0);

    let t_deg: Vec<usize> = if prepared_deg_ge.is_none() {
        target.node_ids().map(|t| target.degree(t)).collect()
    } else {
        Vec::new()
    };
    let mut eligible_buf = if prepared_deg_ge.is_none() || anchor.is_some() {
        vec![0u64; words]
    } else {
        Vec::new()
    };

    let fill_eligible =
        |qd: usize, clear_anchor: Option<usize>, orbit_mask: Option<&[u64]>, buf: &mut [u64]| {
            if let Some((deg_ge, _)) = prepared_deg_ge {
                buf.copy_from_slice(&deg_ge[qd * words..(qd + 1) * words]);
            } else {
                buf.fill(0);
                for (t, &td) in t_deg.iter().enumerate() {
                    if td >= qd {
                        buf[t >> 6] |= 1u64 << (t & 63);
                    }
                }
            }
            if let Some(ta) = clear_anchor {
                buf[ta >> 6] &= !(1u64 << (ta & 63));
            }
            if let Some(orbits) = orbit_mask {
                for (dst, &m) in buf.iter_mut().zip(orbits.iter()) {
                    *dst &= m;
                }
            }
        };

    let anchor_ta = anchor.map(|(_, ta)| ta.index());
    if let Some(num_classes) = pred.query_node_class_count() {
        let mut class_min_deg = vec![usize::MAX; num_classes];
        let mut class_count = vec![0usize; num_classes];
        for &q in &q_by_deg {
            if anchor.is_some_and(|(qa, _)| q == qa.index()) {
                continue;
            }
            let c = pred.query_node_class(NodeId(q as u32));
            debug_assert!(c < num_classes);
            class_min_deg[c] = class_min_deg[c].min(q_deg[q]);
            class_count[c] += 1;
        }
        let mut class_computed = vec![false; num_classes];
        scratch.class_wide.clear();
        scratch.class_wide.resize(num_classes * words, 0);

        for &q in &q_by_deg {
            let base = q * words;
            if let Some((qa, ta)) = anchor {
                if q == qa.index() {
                    let t = ta.index();
                    scratch.cand_wide[base + (t >> 6)] = 1u64 << (t & 63);
                    scratch.cand_counts[q] = 1;
                    continue;
                }
            }
            let c = pred.query_node_class(NodeId(q as u32));
            let cbase = c * words;
            if !class_computed[c] {
                let solo_orbit = root_orbit_filter
                    .filter(|&(p0, _)| class_count[c] == 1 && q == p0)
                    .map(|(_, mask)| mask);
                if let (Some((deg_ge, _)), None, None) = (prepared_deg_ge, anchor_ta, solo_orbit) {
                    let elig = &deg_ge[class_min_deg[c] * words..(class_min_deg[c] + 1) * words];
                    pred.fill_node_candidates_words(
                        NodeId(q as u32),
                        elig,
                        &mut scratch.class_wide[cbase..cbase + words],
                    );
                } else {
                    if eligible_buf.len() != words {
                        eligible_buf.resize(words, 0);
                    }
                    fill_eligible(class_min_deg[c], anchor_ta, solo_orbit, &mut eligible_buf);
                    pred.fill_node_candidates_words(
                        NodeId(q as u32),
                        &eligible_buf,
                        &mut scratch.class_wide[cbase..cbase + words],
                    );
                }
                class_computed[c] = true;
            }

            let qd = q_deg[q];
            let orbit_slice = root_orbit_filter
                .filter(|&(p0, _)| q == p0)
                .map(|(_, mask)| mask);
            let mut cnt = 0usize;
            if let Some((deg_ge, _)) = prepared_deg_ge {
                let deg_slice = &deg_ge[qd * words..(qd + 1) * words];
                for w in 0..words {
                    let mut word = scratch.class_wide[cbase + w] & deg_slice[w];
                    if let Some(orbits) = orbit_slice {
                        word &= orbits[w];
                    }
                    scratch.cand_wide[base + w] = word;
                    cnt += word.count_ones() as usize;
                }
            } else {
                fill_eligible(qd, anchor_ta, orbit_slice, &mut eligible_buf);
                for (w, &elig_word) in eligible_buf.iter().enumerate().take(words) {
                    let word = scratch.class_wide[cbase + w] & elig_word;
                    scratch.cand_wide[base + w] = word;
                    cnt += word.count_ones() as usize;
                }
            }
            if cnt == 0 {
                return false;
            }
            scratch.cand_counts[q] = cnt;
        }
    } else {
        for &q in &q_by_deg {
            let base = q * words;
            if let Some((qa, ta)) = anchor {
                if q == qa.index() {
                    let t = ta.index();
                    scratch.cand_wide[base + (t >> 6)] = 1u64 << (t & 63);
                    scratch.cand_counts[q] = 1;
                    continue;
                }
            }
            let qd = q_deg[q];
            let orbit_slice = root_orbit_filter
                .filter(|&(p0, _)| q == p0)
                .map(|(_, mask)| mask);
            if let (Some((deg_ge, _)), None, None) = (prepared_deg_ge, anchor_ta, orbit_slice) {
                let elig = &deg_ge[qd * words..(qd + 1) * words];
                pred.fill_node_candidates_words(
                    NodeId(q as u32),
                    elig,
                    &mut scratch.cand_wide[base..base + words],
                );
            } else {
                if eligible_buf.len() != words {
                    eligible_buf.resize(words, 0);
                }
                fill_eligible(qd, anchor_ta, orbit_slice, &mut eligible_buf);
                pred.fill_node_candidates_words(
                    NodeId(q as u32),
                    &eligible_buf,
                    &mut scratch.cand_wide[base..base + words],
                );
            }
            let mut cnt = 0usize;
            for &word in &scratch.cand_wide[base..base + words] {
                cnt += word.count_ones() as usize;
            }
            if cnt == 0 {
                return false;
            }
            scratch.cand_counts[q] = cnt;
        }
    }
    true
}
