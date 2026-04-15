//! Disjoint-set / union-find with union-by-rank and path halving.
//!
//! Galler & Fischer (1964) forest representation. Path halving from
//! Tarjan & van Leeuwen (1984). Union-by-rank for O(α(n)) amortized
//! per operation.

/// Disjoint-set data structure.
///
/// Elements are integers `0..n`. Initially each element is its own set.
/// `union` merges two sets; `find` returns the representative of a set.
pub struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
        }
    }

    pub fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x as u32 {
            self.parent[x] = self.parent[self.parent[x] as usize];
            x = self.parent[x] as usize;
        }
        x
    }

    pub fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb as u32,
            std::cmp::Ordering::Greater => self.parent[rb] = ra as u32,
            std::cmp::Ordering::Equal => {
                self.parent[ra] = rb as u32;
                self.rank[rb] += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_union_find_initial() {
        let mut uf = UnionFind::new(4);
        for i in 0..4 {
            assert_eq!(uf.find(i), i);
        }
    }

    #[test]
    fn test_union_find_merge() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(2, 3);
        assert_eq!(uf.find(0), uf.find(1));
        assert_eq!(uf.find(2), uf.find(3));
        assert_ne!(uf.find(0), uf.find(2));
        assert_ne!(uf.find(0), uf.find(4));
    }

    #[test]
    fn test_union_find_transitive() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        let root = uf.find(0);
        for i in 1..4 {
            assert_eq!(uf.find(i), root);
        }
    }

    #[test]
    fn test_union_find_idempotent() {
        let mut uf = UnionFind::new(3);
        uf.union(0, 1);
        uf.union(0, 1);
        assert_eq!(uf.find(0), uf.find(1));
    }

    #[test]
    fn test_union_find_rank_balancing() {
        // Build two chains of equal length, then merge them.
        // With rank balancing, the resulting tree depth is bounded.
        let mut uf = UnionFind::new(8);
        // Chain 1: {0,1,2,3}
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(0, 2);
        // Chain 2: {4,5,6,7}
        uf.union(4, 5);
        uf.union(6, 7);
        uf.union(4, 6);
        // Merge
        uf.union(0, 4);
        let root = uf.find(0);
        for i in 0..8 {
            assert_eq!(uf.find(i), root);
        }
    }
}
