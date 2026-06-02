//! Relation sets: N-ary relations over typed participants (`RelationParticipant`
//! — a `NodeId` or `EdgeId`), each carrying a shared union incidence index (a
//! node index and an edge index) routed from every participant's `refs()`.
//!
//! `FixedRelationSet<P, O, D, N>` stores relations of compile-time-known arity,
//! `VarRelationSet<P, O, D>` stores variable-arity relations.
//! Participants are typed `P` (`RelationParticipant`);
//! the factor ordering `O` (`Unordered`/`Ordered`) controls canonicalization.
//!
//! `FixedFixedBirelationSet`, `FixedVarBirelationSet`, and `VarVarBirelationSet`
//! relate two factors, each with its own participant type, ordering, and arity.
//! The union incidence spans both factors, so a relation is reachable from any
//! of its participants regardless of id-space.

use std::hash::Hash;
use std::marker::PhantomData;

use crate::graph::{EdgeId, NodeId, Remapping};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Canonicalization of a relation factor's participants, applied on
/// construction and after a remap relabels them. `Unordered` sorts (membership
/// is the datum); `Ordered` preserves input order (position is the datum).
pub trait FactorOrdering {
    fn canonicalize<P: Ord>(participants: &mut [P]);
}

/// Set-valued factor: participants are canonicalized by sorting.
#[derive(Clone, Copy, Debug)]
pub struct Unordered;

/// Positional factor: canonicalization is a no-op, input order is preserved.
#[derive(Clone, Copy, Debug)]
pub struct Ordered;

impl FactorOrdering for Unordered {
    fn canonicalize<P: Ord>(participants: &mut [P]) {
        participants.sort_unstable();
    }
}

impl FactorOrdering for Ordered {
    fn canonicalize<P: Ord>(_participants: &mut [P]) {}
}

/// The id-space contents of a participant, surfaced for the incidence index.
/// At most one ref per space today (a node or an edge); a future port type
/// could fill both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ParticipantRefs {
    pub node: Option<NodeId>,
    pub edge: Option<EdgeId>,
}

/// A value that can occupy a relation factor: routes through a `Remapping` in
/// both directions and exposes its node/edge refs for incidence. One impl per
/// concrete id type — dispatch is static, since a factor is homogeneous.
pub trait RelationParticipant: Copy + Ord + Hash {
    fn remap(self, remapping: &Remapping) -> Option<Self>;
    fn unmap(self, remapping: &Remapping) -> Self;
    fn refs(self) -> ParticipantRefs;
}

impl RelationParticipant for NodeId {
    fn remap(self, remapping: &Remapping) -> Option<Self> {
        remapping.map_node(self)
    }

    fn unmap(self, remapping: &Remapping) -> Self {
        remapping.unmap_node(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: Some(self),
            edge: None,
        }
    }
}

impl RelationParticipant for EdgeId {
    fn remap(self, remapping: &Remapping) -> Option<Self> {
        remapping.map_edge(self)
    }

    fn unmap(self, remapping: &Remapping) -> Self {
        remapping.unmap_edge(self)
    }

    fn refs(self) -> ParticipantRefs {
        ParticipantRefs {
            node: None,
            edge: Some(self),
        }
    }
}

/// Union incidence index: a node → relations pair and an edge → relations pair,
/// each a sorted `(keys, rels)` slice. Participants self-route via `refs()`, so a
/// set with only node participants leaves the edge half empty, and vice versa.
#[derive(Clone, Debug, Default)]
struct Incidence {
    node_keys: Vec<NodeId>,
    node_rels: Vec<RelationId>,
    edge_keys: Vec<EdgeId>,
    edge_rels: Vec<RelationId>,
}

impl Incidence {
    /// `fill(i, out)` pushes every participant's `refs()` for relation `i`.
    fn build(
        relation_count: usize,
        mut fill: impl FnMut(usize, &mut Vec<ParticipantRefs>),
    ) -> Self {
        let mut node_entries: Vec<(NodeId, RelationId)> = Vec::new();
        let mut edge_entries: Vec<(EdgeId, RelationId)> = Vec::new();
        let mut refs: Vec<ParticipantRefs> = Vec::new();
        let mut nodes: Vec<NodeId> = Vec::new();
        let mut edges: Vec<EdgeId> = Vec::new();
        for i in 0..relation_count {
            let rid = RelationId(i as u32);
            refs.clear();
            fill(i, &mut refs);
            nodes.clear();
            edges.clear();
            for r in &refs {
                if let Some(node) = r.node {
                    nodes.push(node);
                }
                if let Some(edge) = r.edge {
                    edges.push(edge);
                }
            }
            nodes.sort_unstable();
            nodes.dedup();
            edges.sort_unstable();
            edges.dedup();
            node_entries.extend(nodes.iter().map(|&n| (n, rid)));
            edge_entries.extend(edges.iter().map(|&e| (e, rid)));
        }
        node_entries.sort_by_key(|&(k, _)| k);
        edge_entries.sort_by_key(|&(k, _)| k);
        Self {
            node_keys: node_entries.iter().map(|&(k, _)| k).collect(),
            node_rels: node_entries.iter().map(|&(_, r)| r).collect(),
            edge_keys: edge_entries.iter().map(|&(k, _)| k).collect(),
            edge_rels: edge_entries.iter().map(|&(_, r)| r).collect(),
        }
    }

    fn incident(&self, node: NodeId) -> &[RelationId] {
        let start = self.node_keys.partition_point(|n| *n < node);
        let end = start + self.node_keys[start..].partition_point(|n| *n <= node);
        &self.node_rels[start..end]
    }

    fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        let start = self.edge_keys.partition_point(|e| *e < edge);
        let end = start + self.edge_keys[start..].partition_point(|e| *e <= edge);
        &self.edge_rels[start..end]
    }

    fn has_incident(&self, node: NodeId) -> bool {
        self.node_keys.binary_search(&node).is_ok()
    }

    fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.edge_keys.binary_search(&edge).is_ok()
    }
}

/// Fixed-arity relation set. Each relation connects exactly N nodes.
///
/// Flat CSR storage: participants are `Vec<[NodeId; N]>`, incidence is
/// a flat array with offset table. No heap allocations per node or
/// per relation.
#[derive(Clone, Debug)]
pub struct FixedRelationSet<P, O, D, const N: usize> {
    participants: Vec<[P; N]>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<O>,
}

impl<P: PartialEq, O, D: PartialEq, const N: usize> PartialEq for FixedRelationSet<P, O, D, N> {
    fn eq(&self, other: &Self) -> bool {
        self.participants == other.participants && self.data == other.data
    }
}

impl<P: Eq, O, D: Eq, const N: usize> Eq for FixedRelationSet<P, O, D, N> {}

impl<P: RelationParticipant, O: FactorOrdering, D, const N: usize> FixedRelationSet<P, O, D, N> {
    pub fn new(entries: Vec<([P; N], D)>) -> Self {
        let mut participants = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (mut p, d) in entries {
            O::canonicalize(&mut p);
            participants.push(p);
            data.push(d);
        }

        let incidence = Incidence::build(participants.len(), |i, out| {
            out.extend(participants[i].iter().map(|p| p.refs()));
        });

        Self {
            participants,
            data,
            incidence,
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn participants(&self, id: RelationId) -> &[P; N] {
        &self.participants[id.index()]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([P; N], D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let parts: [P; N] = parts?.try_into().ok()?;
                Some((parts, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<P, O, D, const N: usize> Default for FixedRelationSet<P, O, D, N> {
    fn default() -> Self {
        Self {
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Variable-arity relation set. Each relation connects an arbitrary
/// number of nodes. Participants are sorted by `NodeId` on construction.
///
/// Flat CSR storage: participant ranges via offset table, incidence
/// via a second offset table. No heap allocations per node or per
/// relation.
#[derive(Clone, Debug)]
pub struct VarRelationSet<P, O, D> {
    offsets: Vec<u32>,
    participants: Vec<P>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<O>,
}

impl<P: PartialEq, O, D: PartialEq> PartialEq for VarRelationSet<P, O, D> {
    fn eq(&self, other: &Self) -> bool {
        self.offsets == other.offsets
            && self.participants == other.participants
            && self.data == other.data
    }
}

impl<P: Eq, O, D: Eq> Eq for VarRelationSet<P, O, D> {}

impl<P: RelationParticipant, O: FactorOrdering, D> VarRelationSet<P, O, D> {
    pub fn new(entries: Vec<(Vec<P>, D)>) -> Self {
        let relation_count = entries.len();
        let mut offsets = Vec::with_capacity(relation_count + 1);
        offsets.push(0);

        let total_participants: usize = entries.iter().map(|(p, _)| p.len()).sum();
        let mut participants = Vec::with_capacity(total_participants);
        let mut data = Vec::with_capacity(relation_count);

        for (mut p, d) in entries {
            O::canonicalize(&mut p);
            participants.extend_from_slice(&p);
            offsets.push(participants.len() as u32);
            data.push(d);
        }

        let incidence = Incidence::build(relation_count, |i, out| {
            let start = offsets[i] as usize;
            let end = offsets[i + 1] as usize;
            out.extend(participants[start..end].iter().map(|p| p.refs()));
        });

        Self {
            offsets,
            participants,
            data,
            incidence,
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn participants(&self, id: RelationId) -> &[P] {
        let start = self.offsets[id.index()] as usize;
        let end = self.offsets[id.index() + 1] as usize;
        &self.participants[start..end]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<(Vec<P>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let parts: Option<Vec<P>> = self
                    .participants(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                Some((parts?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<P, O, D> Default for VarRelationSet<P, O, D> {
    fn default() -> Self {
        Self {
            offsets: vec![0],
            participants: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with two fixed-arity factors.
#[derive(Clone, Debug)]
pub struct FixedFixedBirelationSet<L1, O1, const N1: usize, L2, O2, const N2: usize, D> {
    factor_1: Vec<[L1; N1]>,
    factor_2: Vec<[L2; N2]>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> PartialEq
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.factor_1 == other.factor_1
            && self.factor_2 == other.factor_2
            && self.data == other.data
    }
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> Eq
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D>
    FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<([L1; N1], [L2; N2], D)>) -> Self {
        let relation_count = entries.len();
        let mut factor_1 = Vec::with_capacity(relation_count);
        let mut factor_2 = Vec::with_capacity(relation_count);
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            O2::canonicalize(&mut l2);
            factor_1.push(l1);
            factor_2.push(l2);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            out.extend(factor_1[i].iter().map(|p| p.refs()));
            out.extend(factor_2[i].iter().map(|p| p.refs()));
        });
        Self {
            factor_1,
            factor_2,
            data,
            incidence,
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn factor_1(&self, id: RelationId) -> &[L1; N1] {
        &self.factor_1[id.index()]
    }

    pub fn factor_2(&self, id: RelationId) -> &[L2; N2] {
        &self.factor_2[id.index()]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([L1; N1], [L2; N2], D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .factor_1(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let f1: [L1; N1] = f1?.try_into().ok()?;
                let f2: Option<Vec<L2>> = self
                    .factor_2(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let f2: [L2; N2] = f2?.try_into().ok()?;
                Some((f1, f2, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<L1, O1, const N1: usize, L2, O2, const N2: usize, D> Default
    for FixedFixedBirelationSet<L1, O1, N1, L2, O2, N2, D>
{
    fn default() -> Self {
        Self {
            factor_1: Vec::new(),
            factor_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with a fixed-arity factor 1 and a variable-arity factor 2. Each factor
/// has its own participant type and ordering; the union incidence spans both.
#[derive(Clone, Debug)]
pub struct FixedVarBirelationSet<L1, O1, const N1: usize, L2, O2, D> {
    factor_1: Vec<[L1; N1]>,
    f2_offsets: Vec<u32>,
    factor_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, const N1: usize, L2, O2, D> PartialEq for FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.factor_1 == other.factor_1
            && self.f2_offsets == other.f2_offsets
            && self.factor_2 == other.factor_2
            && self.data == other.data
    }
}

impl<L1, O1, const N1: usize, L2, O2, D> Eq for FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, const N1: usize, L2, O2, D> FixedVarBirelationSet<L1, O1, N1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<([L1; N1], Vec<L2>, D)>) -> Self {
        let relation_count = entries.len();
        let mut factor_1 = Vec::with_capacity(relation_count);
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut factor_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            factor_1.push(l1);
            O2::canonicalize(&mut l2);
            factor_2.extend_from_slice(&l2);
            f2_offsets.push(factor_2.len() as u32);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            out.extend(factor_1[i].iter().map(|p| p.refs()));
            let start = f2_offsets[i] as usize;
            let end = f2_offsets[i + 1] as usize;
            out.extend(factor_2[start..end].iter().map(|p| p.refs()));
        });
        Self {
            factor_1,
            f2_offsets,
            factor_2,
            data,
            incidence,
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn factor_1(&self, id: RelationId) -> &[L1; N1] {
        &self.factor_1[id.index()]
    }

    pub fn factor_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.factor_2[start..end]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<([L1; N1], Vec<L2>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .factor_1(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let f1: [L1; N1] = f1?.try_into().ok()?;
                let f2: Option<Vec<L2>> = self
                    .factor_2(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                Some((f1, f2?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<L1, O1, const N1: usize, L2, O2, D> Default for FixedVarBirelationSet<L1, O1, N1, L2, O2, D> {
    fn default() -> Self {
        Self {
            factor_1: Vec::new(),
            f2_offsets: vec![0],
            factor_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

/// Birelation with two variable-arity factors.
#[derive(Clone, Debug)]
pub struct VarVarBirelationSet<L1, O1, L2, O2, D> {
    f1_offsets: Vec<u32>,
    factor_1: Vec<L1>,
    f2_offsets: Vec<u32>,
    factor_2: Vec<L2>,
    data: Vec<D>,
    incidence: Incidence,
    _ordering: PhantomData<(O1, O2)>,
}

impl<L1, O1, L2, O2, D> PartialEq for VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: PartialEq,
    L2: PartialEq,
    D: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.f1_offsets == other.f1_offsets
            && self.factor_1 == other.factor_1
            && self.f2_offsets == other.f2_offsets
            && self.factor_2 == other.factor_2
            && self.data == other.data
    }
}

impl<L1, O1, L2, O2, D> Eq for VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: Eq,
    L2: Eq,
    D: Eq,
{
}

impl<L1, O1, L2, O2, D> VarVarBirelationSet<L1, O1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
{
    pub fn new(entries: Vec<(Vec<L1>, Vec<L2>, D)>) -> Self {
        let relation_count = entries.len();
        let mut f1_offsets = Vec::with_capacity(relation_count + 1);
        f1_offsets.push(0);
        let mut factor_1 = Vec::new();
        let mut f2_offsets = Vec::with_capacity(relation_count + 1);
        f2_offsets.push(0);
        let mut factor_2 = Vec::new();
        let mut data = Vec::with_capacity(relation_count);
        for (mut l1, mut l2, d) in entries {
            O1::canonicalize(&mut l1);
            factor_1.extend_from_slice(&l1);
            f1_offsets.push(factor_1.len() as u32);
            O2::canonicalize(&mut l2);
            factor_2.extend_from_slice(&l2);
            f2_offsets.push(factor_2.len() as u32);
            data.push(d);
        }
        let incidence = Incidence::build(relation_count, |i, out| {
            let s1 = f1_offsets[i] as usize;
            let e1 = f1_offsets[i + 1] as usize;
            out.extend(factor_1[s1..e1].iter().map(|p| p.refs()));
            let s2 = f2_offsets[i] as usize;
            let e2 = f2_offsets[i + 1] as usize;
            out.extend(factor_2[s2..e2].iter().map(|p| p.refs()));
        });
        Self {
            f1_offsets,
            factor_1,
            f2_offsets,
            factor_2,
            data,
            incidence,
            _ordering: PhantomData,
        }
    }

    pub fn relation_count(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self, id: RelationId) -> &D {
        &self.data[id.index()]
    }

    pub fn data_mut(&mut self, id: RelationId) -> &mut D {
        &mut self.data[id.index()]
    }

    pub fn data_iter_mut(&mut self) -> impl Iterator<Item = &mut D> {
        self.data.iter_mut()
    }

    pub fn factor_1(&self, id: RelationId) -> &[L1] {
        let start = self.f1_offsets[id.index()] as usize;
        let end = self.f1_offsets[id.index() + 1] as usize;
        &self.factor_1[start..end]
    }

    pub fn factor_2(&self, id: RelationId) -> &[L2] {
        let start = self.f2_offsets[id.index()] as usize;
        let end = self.f2_offsets[id.index() + 1] as usize;
        &self.factor_2[start..end]
    }

    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        self.incidence.incident(node)
    }

    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId] {
        self.incidence.incident_edge(edge)
    }

    pub fn has_incident(&self, node: NodeId) -> bool {
        self.incidence.has_incident(node)
    }

    pub fn has_incident_edge(&self, edge: EdgeId) -> bool {
        self.incidence.has_incident_edge(edge)
    }

    pub fn contains(&self, id: RelationId) -> bool {
        id.index() < self.data.len()
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> {
        (0..self.data.len() as u32).map(RelationId)
    }

    pub fn apply_remapping(&self, remapping: &Remapping) -> Self
    where
        D: Clone,
    {
        let entries: Vec<(Vec<L1>, Vec<L2>, D)> = (0..self.relation_count())
            .filter_map(|i| {
                let rid = RelationId(i as u32);
                let f1: Option<Vec<L1>> = self
                    .factor_1(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                let f2: Option<Vec<L2>> = self
                    .factor_2(rid)
                    .iter()
                    .map(|&p| p.remap(remapping))
                    .collect();
                Some((f1?, f2?, self.data(rid).clone()))
            })
            .collect();
        Self::new(entries)
    }
}

impl<L1, O1, L2, O2, D> Default for VarVarBirelationSet<L1, O1, L2, O2, D> {
    fn default() -> Self {
        Self {
            f1_offsets: vec![0],
            factor_1: Vec::new(),
            f2_offsets: vec![0],
            factor_2: Vec::new(),
            data: Vec::new(),
            incidence: Incidence::default(),
            _ordering: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[rstest]
    #[case::already_sorted(vec![0, 1, 2], vec![0, 1, 2])]
    #[case::reversed(vec![2, 1, 0], vec![0, 1, 2])]
    #[case::shuffled(vec![2, 0, 3, 1], vec![0, 1, 2, 3])]
    fn test_unordered_canonicalize(#[case] mut input: Vec<i32>, #[case] expected: Vec<i32>) {
        Unordered::canonicalize(&mut input);
        assert_eq!(input, expected);
    }

    #[rstest]
    #[case::already_sorted(vec![0, 1, 2])]
    #[case::reversed(vec![2, 1, 0])]
    #[case::shuffled(vec![2, 0, 3, 1])]
    fn test_ordered_canonicalize(#[case] input: Vec<i32>) {
        let mut actual = input.clone();
        Ordered::canonicalize(&mut actual);
        assert_eq!(actual, input);
    }

    #[rstest]
    #[case::before_removed(NodeId(0), Some(NodeId(0)))]
    #[case::removed(NodeId(1), None)]
    #[case::after_removed(NodeId(2), Some(NodeId(1)))]
    fn test_node_id_remap(#[case] id: NodeId, #[case] expected: Option<NodeId>) {
        let remapping = Remapping::new(vec![1], vec![]);
        assert_eq!(id.remap(&remapping), expected);
    }

    #[rstest]
    #[case::before_gap(NodeId(0), NodeId(0))]
    #[case::after_gap(NodeId(1), NodeId(2))]
    fn test_node_id_unmap(#[case] id: NodeId, #[case] expected: NodeId) {
        let remapping = Remapping::new(vec![1], vec![]);
        assert_eq!(id.unmap(&remapping), expected);
    }

    #[rstest]
    fn test_node_id_refs() {
        assert_eq!(
            NodeId(3).refs(),
            ParticipantRefs {
                node: Some(NodeId(3)),
                edge: None,
            }
        );
    }

    #[rstest]
    #[case::removed(EdgeId(0), None)]
    #[case::after_removed(EdgeId(2), Some(EdgeId(1)))]
    fn test_edge_id_remap(#[case] id: EdgeId, #[case] expected: Option<EdgeId>) {
        let remapping = Remapping::new(vec![], vec![0]);
        assert_eq!(id.remap(&remapping), expected);
    }

    #[rstest]
    #[case::before_gap(EdgeId(0), EdgeId(0))]
    #[case::after_gap(EdgeId(1), EdgeId(2))]
    fn test_edge_id_unmap(#[case] id: EdgeId, #[case] expected: EdgeId) {
        let remapping = Remapping::new(vec![], vec![1]);
        assert_eq!(id.unmap(&remapping), expected);
    }

    #[rstest]
    fn test_edge_id_refs() {
        assert_eq!(
            EdgeId(2).refs(),
            ParticipantRefs {
                node: None,
                edge: Some(EdgeId(2)),
            }
        );
    }

    #[rstest]
    fn test_fixed_relation_set_new() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], "dative"), ([n(1), n(2)], "noncov")]);
        assert_eq!(rs.relation_count(), 2);
        assert_eq!(rs.data(RelationId(0)), &"dative");
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(2)]);
    }

    #[rstest]
    fn test_fixed_relation_set_data_mut() {
        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_relation_set_data_iter_mut() {
        let mut rs: FixedRelationSet<NodeId, Unordered, i32, 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], 1),
            ([n(1), n(2)], 2),
            ([n(2), n(3)], 3),
        ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[rstest]
    fn test_fixed_relation_set_participants_sorted() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(3)]);
    }

    #[rstest]
    fn test_fixed_relation_set_participants_ordered() {
        let rs: FixedRelationSet<NodeId, Ordered, &str, 2> =
            FixedRelationSet::new(vec![([n(2), n(0)], "a"), ([n(3), n(1)], "b")]);
        assert_eq!(rs.participants(RelationId(0)), &[n(2), n(0)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(3), n(1)]);
    }

    #[rstest]
    fn test_fixed_relation_set_incidence() {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> = FixedRelationSet::new(vec![
            ([n(0), n(1)], ()),
            ([n(0), n(2)], ()),
            ([n(2), n(3)], ()),
        ]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident(n(1)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(1), RelationId(2)]);
        assert_eq!(rs.incident(n(3)), &[RelationId(2)]);
        assert!(rs.has_incident(n(0)));
        assert!(!rs.has_incident(n(5)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::last(RelationId(1), true)]
    #[case::out_of_range(RelationId(2), false)]
    fn test_fixed_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_relation_set_relation_ids() {
        let rs: FixedRelationSet<NodeId, Unordered, (), 2> =
            FixedRelationSet::new(vec![([n(0), n(1)], ()), ([n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_relation_set_apply_remapping() {
        let rs: FixedRelationSet<NodeId, Unordered, &str, 2> =
            FixedRelationSet::new(vec![([n(0), n(2)], "keep"), ([n(1), n(3)], "drop")]);
        let remapping = Remapping::new(vec![1], vec![]);
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_relation_set_default() {
        let rs = FixedRelationSet::<NodeId, Unordered, (), 2>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_var_relation_set_new() {
        let rs: VarRelationSet<NodeId, Unordered, &str> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2), n(3), n(4), n(5)], "benzene")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.data(RelationId(0)), &"benzene");
        assert_eq!(
            rs.participants(RelationId(0)),
            &[n(0), n(1), n(2), n(3), n(4), n(5)]
        );
    }

    #[rstest]
    fn test_var_relation_set_data_mut() {
        let mut rs: VarRelationSet<NodeId, Unordered, i32> =
            VarRelationSet::new(vec![(vec![n(0), n(1), n(2)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_var_relation_set_data_iter_mut() {
        let mut rs: VarRelationSet<NodeId, Unordered, i32> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], 1),
            (vec![n(2), n(3), n(4)], 2),
            (vec![n(5)], 3),
        ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
        assert_eq!(rs.data(RelationId(2)), &30);
    }

    #[rstest]
    fn test_var_relation_set_participants_sorted() {
        let rs: VarRelationSet<NodeId, Unordered, ()> = VarRelationSet::new(vec![
            (vec![n(5), n(2), n(0), n(3)], ()),
            (vec![n(4), n(1)], ()),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(2), n(3), n(5)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(1), n(4)]);
    }

    #[rstest]
    fn test_var_relation_set_participants_ordered() {
        let rs: VarRelationSet<NodeId, Ordered, ()> = VarRelationSet::new(vec![
            (vec![n(5), n(2), n(0), n(3)], ()),
            (vec![n(4), n(1)], ()),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(5), n(2), n(0), n(3)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(4), n(1)]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(4)), &[RelationId(1)]);
    }

    #[rstest]
    fn test_var_relation_set_variable_arity() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(1)], "pair"),
            (vec![n(2), n(3), n(4), n(5)], "quad"),
        ]);
        assert_eq!(rs.participants(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.participants(RelationId(1)), &[n(2), n(3), n(4), n(5)]);
    }

    #[rstest]
    fn test_var_relation_set_incidence() {
        let rs: VarRelationSet<NodeId, Unordered, ()> = VarRelationSet::new(vec![
            (vec![n(0), n(1), n(2)], ()),
            (vec![n(2), n(3), n(4)], ()),
        ]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident(n(4)), &[RelationId(1)]);
        assert!(rs.has_incident(n(0)));
        assert!(!rs.has_incident(n(7)));
    }

    #[rstest]
    fn test_var_relation_set_edge_incidence() {
        let rs: VarRelationSet<EdgeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![EdgeId(0), EdgeId(2)], "a"),
            (vec![EdgeId(1), EdgeId(2)], "b"),
        ]);
        assert_eq!(rs.incident_edge(EdgeId(2)), &[RelationId(0), RelationId(1)]);
        assert_eq!(rs.incident_edge(EdgeId(0)), &[RelationId(0)]);
        assert!(rs.has_incident_edge(EdgeId(2)));
        assert!(!rs.has_incident_edge(EdgeId(5)));
        assert!(rs.incident(n(0)).is_empty());
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_var_relation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ())]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_relation_set_relation_ids() {
        let rs: VarRelationSet<NodeId, Unordered, ()> =
            VarRelationSet::new(vec![(vec![n(0), n(1)], ()), (vec![n(1), n(2)], ())]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_relation_set_apply_remapping() {
        let rs: VarRelationSet<NodeId, Unordered, &str> = VarRelationSet::new(vec![
            (vec![n(0), n(2), n(4)], "keep"),
            (vec![n(1), n(3)], "drop"),
        ]);
        let remapping = Remapping::new(vec![1], vec![]);
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.participants(RelationId(0)), &[n(0), n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_relation_set_default() {
        let rs = VarRelationSet::<NodeId, Unordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_new() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(2), n(1)], "x")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.factor_1(RelationId(0)), &[n(0)]);
        assert_eq!(rs.factor_2(RelationId(0)), &[n(1), n(2)]);
        assert_eq!(rs.data(RelationId(0)), &"x");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_data_mut() {
        let mut rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_data_iter_mut() {
        let mut rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, i32> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], 1), ([n(2)], [n(3)], 2)]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_incidence() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, EdgeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [EdgeId(7)], "x")]);
        assert_eq!(rs.incident(n(0)), &[RelationId(0)]);
        assert_eq!(rs.incident_edge(EdgeId(7)), &[RelationId(0)]);
        assert!(rs.has_incident(n(0)));
        assert!(rs.has_incident_edge(EdgeId(7)));
        assert!(!rs.has_incident(n(5)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_fixed_fixed_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "x")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_relation_ids() {
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 1, &str> =
            FixedFixedBirelationSet::new(vec![([n(0)], [n(1)], "a"), ([n(2)], [n(3)], "b")]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_apply_remapping() {
        // dropped relation loses a factor-1 participant
        let rs: FixedFixedBirelationSet<NodeId, Unordered, 1, NodeId, Unordered, 2, &str> =
            FixedFixedBirelationSet::new(vec![
                ([n(0)], [n(2), n(4)], "keep"),
                ([n(1)], [n(5), n(6)], "drop"),
            ]);
        let remapping = Remapping::new(vec![1], vec![]);
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.factor_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.factor_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_fixed_birelation_set_default() {
        let rs =
            FixedFixedBirelationSet::<NodeId, Unordered, 1, NodeId, Unordered, 1, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    fn test_fixed_var_birelation_set_new() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2), n(3)], "ct")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.factor_1(RelationId(0)), &[EdgeId(0)]);
        assert_eq!(rs.factor_2(RelationId(0)), &[n(1), n(2), n(3)]);
        assert_eq!(rs.data(RelationId(0)), &"ct");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_data_mut() {
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, i32> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_data_iter_mut() {
        let mut rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, i32> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], 1),
                ([EdgeId(1)], vec![n(2)], 2),
            ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_incidence() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1), n(2)], "ct")]);
        assert_eq!(rs.incident_edge(EdgeId(0)), &[RelationId(0)]);
        assert_eq!(rs.incident(n(2)), &[RelationId(0)]);
        assert!(rs.has_incident_edge(EdgeId(0)));
        assert!(rs.has_incident(n(1)));
        assert!(rs.incident(n(0)).is_empty());
        assert!(!rs.has_incident(n(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_fixed_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![([EdgeId(0)], vec![n(1)], "ct")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_relation_ids() {
        let rs: FixedVarBirelationSet<EdgeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([EdgeId(0)], vec![n(1)], "a"),
                ([EdgeId(1)], vec![n(2)], "b"),
            ]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_fixed_var_birelation_set_apply_remapping() {
        let rs: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Ordered, &str> =
            FixedVarBirelationSet::new(vec![
                ([n(0)], vec![n(2), n(4)], "keep"),
                ([n(5)], vec![n(1), n(3)], "drop"),
            ]);
        let remapping = Remapping::new(vec![1], vec![]);
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.factor_1(RelationId(0)), &[n(0)]);
        assert_eq!(out.factor_2(RelationId(0)), &[n(1), n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_fixed_var_birelation_set_default() {
        let rs = FixedVarBirelationSet::<EdgeId, Ordered, 1, NodeId, Ordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    fn test_var_var_birelation_set_new() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.factor_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(rs.factor_2(RelationId(0)), &[EdgeId(5)]);
        assert_eq!(rs.data(RelationId(0)), &"y");
    }

    #[rstest]
    fn test_var_var_birelation_set_data_mut() {
        let mut rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, i32> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], 1)]);
        *rs.data_mut(RelationId(0)) = 99;
        assert_eq!(rs.data(RelationId(0)), &99);
    }

    #[rstest]
    fn test_var_var_birelation_set_data_iter_mut() {
        let mut rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, i32> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], 1),
                (vec![n(2)], vec![EdgeId(3)], 2),
            ]);
        for d in rs.data_iter_mut() {
            *d *= 10;
        }
        assert_eq!(rs.data(RelationId(0)), &10);
        assert_eq!(rs.data(RelationId(1)), &20);
    }

    #[rstest]
    fn test_var_var_birelation_set_incidence() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0), n(1)], vec![EdgeId(5)], "y")]);
        assert_eq!(rs.incident(n(1)), &[RelationId(0)]);
        assert_eq!(rs.incident_edge(EdgeId(5)), &[RelationId(0)]);
        assert!(rs.has_incident(n(0)));
        assert!(rs.has_incident_edge(EdgeId(5)));
        assert!(!rs.has_incident_edge(EdgeId(0)));
    }

    #[rstest]
    #[case::first(RelationId(0), true)]
    #[case::out_of_range(RelationId(1), false)]
    fn test_var_var_birelation_set_contains(#[case] id: RelationId, #[case] expected: bool) {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![(vec![n(0)], vec![EdgeId(1)], "y")]);
        assert_eq!(rs.contains(id), expected);
    }

    #[rstest]
    fn test_var_var_birelation_set_relation_ids() {
        let rs: VarVarBirelationSet<NodeId, Unordered, EdgeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0)], vec![EdgeId(1)], "a"),
                (vec![n(2)], vec![EdgeId(3)], "b"),
            ]);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![RelationId(0), RelationId(1)]);
    }

    #[rstest]
    fn test_var_var_birelation_set_apply_remapping() {
        // dropped relation loses a factor-2 participant
        let rs: VarVarBirelationSet<NodeId, Unordered, NodeId, Unordered, &str> =
            VarVarBirelationSet::new(vec![
                (vec![n(0), n(2)], vec![n(4)], "keep"),
                (vec![n(5)], vec![n(1)], "drop"),
            ]);
        let remapping = Remapping::new(vec![1], vec![]);
        let out = rs.apply_remapping(&remapping);
        assert_eq!(out.relation_count(), 1);
        assert_eq!(out.factor_1(RelationId(0)), &[n(0), n(1)]);
        assert_eq!(out.factor_2(RelationId(0)), &[n(3)]);
        assert_eq!(out.data(RelationId(0)), &"keep");
    }

    #[rstest]
    fn test_var_var_birelation_set_default() {
        let rs = VarVarBirelationSet::<NodeId, Unordered, EdgeId, Unordered, ()>::default();
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.has_incident(n(0)));
    }
}
