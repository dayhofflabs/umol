use crate::graph::NodeId;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RelationId(pub u32);

impl RelationId {
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RelationData<R> {
    participants: Vec<NodeId>,
    data: R,
}

impl<R: PartialEq> PartialEq for RelationSet<R> {
    fn eq(&self, other: &Self) -> bool {
        self.relation_count == other.relation_count && self.relations == other.relations
    }
}

impl<R: Eq> Eq for RelationSet<R> {}

/// A typed set of relations (hyperedges) over a shared `NodeId` space.
///
/// Each relation connects one or more nodes and carries data of type `R`.
/// Per-node incidence lists enable O(incident) removal when a node is
/// deleted from the parent graph.
#[derive(Clone, Debug)]
pub struct RelationSet<R> {
    relations: Vec<Option<RelationData<R>>>,
    incidence: Vec<Vec<RelationId>>,
    relation_count: usize,
    free: Vec<RelationId>,
}

impl<R> RelationSet<R> {
    pub fn new() -> Self {
        Self {
            relations: Vec::new(),
            incidence: Vec::new(),
            relation_count: 0,
            free: Vec::new(),
        }
    }

    /// Ensure incidence tracking covers at least `node_bound` node slots.
    /// Call this when the parent graph grows.
    pub fn ensure_node_bound(&mut self, node_bound: usize) {
        if self.incidence.len() < node_bound {
            self.incidence.resize_with(node_bound, Vec::new);
        }
    }

    pub fn add(&mut self, participants: Vec<NodeId>, data: R) -> RelationId {
        let id = if let Some(id) = self.free.pop() {
            self.relations[id.index()] = Some(RelationData {
                participants: participants.clone(),
                data,
            });
            id
        } else {
            let id = RelationId(self.relations.len() as u32);
            self.relations.push(Some(RelationData {
                participants: participants.clone(),
                data,
            }));
            id
        };

        for &node in &participants {
            let idx = node.index();
            if idx >= self.incidence.len() {
                self.incidence.resize_with(idx + 1, Vec::new);
            }
            self.incidence[idx].push(id);
        }

        self.relation_count += 1;
        id
    }

    pub fn remove(&mut self, id: RelationId) -> Option<R> {
        let rel = self.relations.get_mut(id.index())?.take()?;
        for &node in &rel.participants {
            if node.index() < self.incidence.len() {
                self.incidence[node.index()].retain(|&r| r != id);
            }
        }
        self.free.push(id);
        self.relation_count -= 1;
        Some(rel.data)
    }

    /// Remove all relations that reference `node`. Returns the count removed.
    pub fn remove_participant(&mut self, node: NodeId) -> usize {
        if node.index() >= self.incidence.len() {
            return 0;
        }
        let incident: Vec<RelationId> = self.incidence[node.index()].drain(..).collect();
        let mut removed = 0;
        for rel_id in incident {
            if let Some(rel) = self.relations[rel_id.index()].take() {
                for &other in &rel.participants {
                    if other != node && other.index() < self.incidence.len() {
                        self.incidence[other.index()].retain(|&r| r != rel_id);
                    }
                }
                self.free.push(rel_id);
                self.relation_count -= 1;
                removed += 1;
            }
        }
        removed
    }

    pub fn contains(&self, id: RelationId) -> bool {
        self.relations
            .get(id.index())
            .is_some_and(|r| r.is_some())
    }

    pub fn data(&self, id: RelationId) -> Option<&R> {
        self.relations.get(id.index())?.as_ref().map(|r| &r.data)
    }

    pub fn data_mut(&mut self, id: RelationId) -> Option<&mut R> {
        self.relations
            .get_mut(id.index())?
            .as_mut()
            .map(|r| &mut r.data)
    }

    pub fn participants(&self, id: RelationId) -> Option<&[NodeId]> {
        self.relations
            .get(id.index())?
            .as_ref()
            .map(|r| r.participants.as_slice())
    }

    /// Relations incident to a node.
    pub fn incident(&self, node: NodeId) -> &[RelationId] {
        if node.index() < self.incidence.len() {
            &self.incidence[node.index()]
        } else {
            &[]
        }
    }

    /// Whether a node participates in any relation in this set.
    pub fn has_incident(&self, node: NodeId) -> bool {
        !self.incident(node).is_empty()
    }

    pub fn relation_count(&self) -> usize {
        self.relation_count
    }

    pub fn relation_ids(&self) -> impl Iterator<Item = RelationId> + '_ {
        self.relations
            .iter()
            .enumerate()
            .filter_map(|(i, r)| r.as_ref().map(|_| RelationId(i as u32)))
    }
}

impl<R> Default for RelationSet<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    #[test]
    fn test_relation_set_add() {
        let mut rs = RelationSet::<&str>::new();
        let r = rs.add(vec![n(0), n(1)], "dative");
        assert_eq!(rs.relation_count(), 1);
        assert_eq!(rs.data(r), Some(&"dative"));
        assert_eq!(rs.participants(r), Some(vec![n(0), n(1)].as_slice()));
    }

    #[test]
    fn test_relation_set_remove() {
        let mut rs = RelationSet::<()>::new();
        let r = rs.add(vec![n(0), n(1), n(2)], ());
        assert_eq!(rs.remove(r), Some(()));
        assert_eq!(rs.relation_count(), 0);
        assert!(!rs.contains(r));
        assert!(!rs.has_incident(n(0)));
        assert!(!rs.has_incident(n(1)));
        assert!(!rs.has_incident(n(2)));
    }

    #[test]
    fn test_relation_set_remove_participant() {
        let mut rs = RelationSet::<&str>::new();
        let r1 = rs.add(vec![n(0), n(1)], "bond_a");
        let r2 = rs.add(vec![n(1), n(2)], "bond_b");
        let r3 = rs.add(vec![n(2), n(3)], "bond_c");

        let removed = rs.remove_participant(n(1));
        assert_eq!(removed, 2);
        assert_eq!(rs.relation_count(), 1);
        assert!(!rs.contains(r1));
        assert!(!rs.contains(r2));
        assert!(rs.contains(r3));
        // n(0) was co-participant in r1, should have no incident relations left
        assert!(!rs.has_incident(n(0)));
        // n(2) was in r2 (removed) and r3 (kept)
        assert_eq!(rs.incident(n(2)), &[RelationId(2)]);
    }

    #[test]
    fn test_relation_set_nary() {
        let mut rs = RelationSet::<()>::new();
        let r = rs.add(vec![n(0), n(1), n(2), n(3), n(4), n(5)], ());
        assert!(rs.has_incident(n(3)));
        assert_eq!(rs.participants(r).unwrap().len(), 6);

        rs.remove_participant(n(2));
        assert_eq!(rs.relation_count(), 0);
        // All co-participants cleaned
        for i in 0..6 {
            assert!(!rs.has_incident(n(i)));
        }
    }

    #[test]
    fn test_relation_set_free_list_reuse() {
        let mut rs = RelationSet::<i32>::new();
        let r1 = rs.add(vec![n(0)], 10);
        rs.remove(r1);
        let r2 = rs.add(vec![n(1)], 20);
        assert_eq!(r2, RelationId(0));
        assert_eq!(rs.data(r2), Some(&20));
    }

    #[test]
    fn test_relation_set_incident() {
        let mut rs = RelationSet::<()>::new();
        let r1 = rs.add(vec![n(0), n(1)], ());
        let r2 = rs.add(vec![n(0), n(2)], ());
        let _r3 = rs.add(vec![n(3), n(4)], ());

        let mut inc: Vec<RelationId> = rs.incident(n(0)).to_vec();
        inc.sort();
        assert_eq!(inc, vec![r1, r2]);

        assert!(rs.incident(n(3)).len() == 1);
        assert!(rs.incident(n(5)).is_empty());
    }

    #[test]
    fn test_relation_set_data_mut() {
        let mut rs = RelationSet::<i32>::new();
        let r = rs.add(vec![n(0)], 1);
        *rs.data_mut(r).unwrap() = 99;
        assert_eq!(rs.data(r), Some(&99));
    }

    #[test]
    fn test_relation_set_remove_nonexistent() {
        let mut rs = RelationSet::<()>::new();
        assert_eq!(rs.remove(RelationId(0)), None);
        assert_eq!(rs.remove_participant(n(0)), 0);
    }

    #[test]
    fn test_relation_set_relation_ids() {
        let mut rs = RelationSet::<()>::new();
        let r1 = rs.add(vec![n(0)], ());
        let r2 = rs.add(vec![n(1)], ());
        let r3 = rs.add(vec![n(2)], ());
        rs.remove(r2);
        let ids: Vec<RelationId> = rs.relation_ids().collect();
        assert_eq!(ids, vec![r1, r3]);
    }
}
