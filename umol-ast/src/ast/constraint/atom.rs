//! Per-atom constraints.

use std::mem::replace;

use smallvec::SmallVec;
use strum::{EnumCount, EnumDiscriminants, EnumIter};

use super::super::remap::IdxRemapping;
use super::super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash, EnumCount, EnumIter))]
#[repr(u8)]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceAst),
    MulticenterValence(MulticenterValenceAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    Connectivity(ValueAst),
    RingConnectivity(ValueAst),
    TotalHydrogens(ValueAst),
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl AtomConstraint {
    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }

    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Valence(v)
            | Self::DonatedPairs(v)
            | Self::AcceptedPairs(v)
            | Self::Degree(v)
            | Self::Connectivity(v)
            | Self::RingConnectivity(v)
            | Self::TotalHydrogens(v)
            | Self::RingCount(v)
            | Self::RingSize(v) => v.is_undetermined(),
            Self::AromaticValence(c) => c.is_undetermined(),
            Self::MulticenterValence(c) => c.is_undetermined(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum AromaticValenceAst {
    #[default]
    Undetermined,
    NotAromatic,
    Aromatic(ValueAst),
}

impl AromaticValenceAst {
    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum MulticenterValenceAst {
    #[default]
    Undetermined,
    NotMulticenter,
    Multicenter(ValueAst),
}

impl MulticenterValenceAst {
    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }
}

/// Per-atom constraints: at most one entry per [`AtomConstraintKind`].
/// Stored kind-sorted in an inline-capacity-2 `SmallVec`; the common cases
/// after resolution (0–2 constraints) never touch the heap.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomConstraints {
    entries: SmallVec<[AtomConstraint; 2]>,
}

impl AtomConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn contains(&self, kind: AtomConstraintKind) -> bool {
        self.find(kind).is_ok()
    }

    pub fn get(&self, kind: AtomConstraintKind) -> Option<&AtomConstraint> {
        self.find(kind).ok().map(|i| &self.entries[i])
    }

    pub fn get_mut(&mut self, kind: AtomConstraintKind) -> Option<&mut AtomConstraint> {
        match self.find(kind) {
            Ok(i) => Some(&mut self.entries[i]),
            Err(_) => None,
        }
    }

    /// Insert a constraint at its kind's sorted position, returning the
    /// previous entry of the same kind if any. Every `AtomConstraintKind` is
    /// single-valued per atom, so `add` always replaces same-kind entries
    /// (last-wins).
    pub fn add(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        match self.find(constraint.kind()) {
            Ok(i) => Some(replace(&mut self.entries[i], constraint)),
            Err(i) => {
                self.entries.insert(i, constraint);
                None
            }
        }
    }

    pub fn retain(&mut self, mut f: impl FnMut(&AtomConstraint) -> bool) {
        self.entries.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.entries.clear();
    }

    pub fn remove(&mut self, kind: AtomConstraintKind) -> Option<AtomConstraint> {
        self.find(kind).ok().map(|i| self.entries.remove(i))
    }

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.entries.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut AtomConstraint> {
        self.entries.iter_mut()
    }

    /// No-op: no `AtomConstraint` variant carries an entity index.
    pub fn remap(self, _remap: &IdxRemapping) -> Self {
        self
    }

    fn find(&self, kind: AtomConstraintKind) -> Result<usize, usize> {
        self.entries
            .binary_search_by_key(&(kind as u8), |c| c.kind() as u8)
    }
}

impl FromIterator<AtomConstraint> for AtomConstraints {
    fn from_iter<I: IntoIterator<Item = AtomConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.add(c);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::Valence(ValueAst::Lit(4)), AtomConstraintKind::Valence)]
    #[case::aromatic_valence(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic), AtomConstraintKind::AromaticValence)]
    #[case::multicenter_valence(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined), AtomConstraintKind::MulticenterValence)]
    #[case::donated_pairs(AtomConstraint::DonatedPairs(ValueAst::Lit(1)), AtomConstraintKind::DonatedPairs)]
    #[case::accepted_pairs(AtomConstraint::AcceptedPairs(ValueAst::Lit(2)), AtomConstraintKind::AcceptedPairs)]
    #[case::degree(AtomConstraint::Degree(ValueAst::Lit(3)), AtomConstraintKind::Degree)]
    #[case::connectivity(AtomConstraint::Connectivity(ValueAst::Lit(4)), AtomConstraintKind::Connectivity)]
    #[case::ring_connectivity(AtomConstraint::RingConnectivity(ValueAst::Lit(2)), AtomConstraintKind::RingConnectivity)]
    #[case::total_hydrogens(AtomConstraint::TotalHydrogens(ValueAst::Lit(3)), AtomConstraintKind::TotalHydrogens)]
    #[case::ring_count(AtomConstraint::RingCount(ValueAst::Lit(1)), AtomConstraintKind::RingCount)]
    #[case::ring_size(AtomConstraint::RingSize(ValueAst::Lit(6)), AtomConstraintKind::RingSize)]
    fn test_atom_constraint_kind(
        #[case] constraint: AtomConstraint,
        #[case] expected: AtomConstraintKind,
    ) {
        assert_eq!(constraint.kind(), expected);
    }

    #[rstest]
    fn test_atom_constraints_new() {
        let cs = AtomConstraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.iter().collect::<Vec<_>>(), Vec::<&AtomConstraint>::new());
    }

    #[rstest]
    fn test_atom_constraints_set_inserts_and_returns_none() {
        let mut cs = AtomConstraints::new();
        let prev = cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        assert_eq!(prev, None);
        assert!(cs.contains(AtomConstraintKind::Valence));
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(4)))
        );
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_atom_constraints_set_replaces_same_kind() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(3)));
        let prev = cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        assert_eq!(prev, Some(AtomConstraint::Valence(ValueAst::Lit(3))));
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(4)))
        );
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_atom_constraints_set_different_kinds_coexist() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.add(AtomConstraint::Degree(ValueAst::Lit(3)));
        cs.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::NotAromatic,
        ));
        assert_eq!(cs.len(), 3);
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(4)))
        );
        assert_eq!(
            cs.get(AtomConstraintKind::Degree),
            Some(&AtomConstraint::Degree(ValueAst::Lit(3)))
        );
        assert_eq!(
            cs.get(AtomConstraintKind::AromaticValence),
            Some(&AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic
            ))
        );
    }

    #[rstest]
    fn test_atom_constraints_contains_absent() {
        let cs = AtomConstraints::new();
        assert!(!cs.contains(AtomConstraintKind::Valence));
        assert_eq!(cs.get(AtomConstraintKind::Valence), None);
    }

    #[rstest]
    fn test_atom_constraints_remove_present() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        let removed = cs.remove(AtomConstraintKind::Valence);
        assert_eq!(removed, Some(AtomConstraint::Valence(ValueAst::Lit(4))));
        assert!(!cs.contains(AtomConstraintKind::Valence));
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_atom_constraints_remove_absent() {
        let mut cs = AtomConstraints::new();
        assert_eq!(cs.remove(AtomConstraintKind::Valence), None);
    }

    #[rstest]
    fn test_atom_constraints_get_mut() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(3)));
        let slot = cs.get_mut(AtomConstraintKind::Valence).unwrap();
        *slot = AtomConstraint::Valence(ValueAst::Lit(5));
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(5)))
        );
    }

    #[rstest]
    fn test_atom_constraints_get_mut_absent() {
        let mut cs = AtomConstraints::new();
        assert!(cs.get_mut(AtomConstraintKind::Valence).is_none());
    }

    #[rstest]
    fn test_atom_constraints_iter_skips_empty_slots() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.add(AtomConstraint::RingSize(ValueAst::Lit(6)));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                AtomConstraint::Valence(ValueAst::Lit(4)),
                AtomConstraint::RingSize(ValueAst::Lit(6)),
            ]
        );
    }

    #[rstest]
    fn test_atom_constraints_iter_mut_allows_mutation() {
        let mut cs = AtomConstraints::new();
        cs.add(AtomConstraint::Valence(ValueAst::Lit(3)));
        cs.add(AtomConstraint::Degree(ValueAst::Lit(2)));
        for c in cs.iter_mut() {
            if let AtomConstraint::Valence(v) = c {
                *v = ValueAst::Lit(7);
            }
        }
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(7)))
        );
        assert_eq!(
            cs.get(AtomConstraintKind::Degree),
            Some(&AtomConstraint::Degree(ValueAst::Lit(2)))
        );
    }

    #[rstest]
    fn test_atom_constraints_from_iter_deduplicates_same_kind() {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(3)),
            AtomConstraint::Valence(ValueAst::Lit(4)),
        ]);
        assert_eq!(cs.len(), 1);
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(4)))
        );
    }

    #[rstest]
    fn test_atom_constraints_from_iter_empty() {
        let cs = AtomConstraints::from_iter([]);
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_atom_constraints_from_iter_preserves_distinct_kinds() {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(4)),
            AtomConstraint::Degree(ValueAst::Lit(3)),
            AtomConstraint::RingCount(ValueAst::Lit(2)),
        ]);
        assert_eq!(cs.len(), 3);
        assert_eq!(
            cs.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(4)))
        );
        assert_eq!(
            cs.get(AtomConstraintKind::Degree),
            Some(&AtomConstraint::Degree(ValueAst::Lit(3)))
        );
        assert_eq!(
            cs.get(AtomConstraintKind::RingCount),
            Some(&AtomConstraint::RingCount(ValueAst::Lit(2)))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence_lit(AtomConstraint::Valence(ValueAst::Lit(4)), false)]
    #[case::valence_undetermined(AtomConstraint::Valence(ValueAst::Undetermined), true)]
    #[case::degree_undetermined(AtomConstraint::Degree(ValueAst::Undetermined), true)]
    #[case::ring_size_undetermined(AtomConstraint::RingSize(ValueAst::Undetermined), true)]
    #[case::aromatic_undetermined(AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined), true)]
    #[case::aromatic_not_aromatic(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic), false)]
    #[case::aromatic_with_value(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1))), false)]
    #[case::multicenter_undetermined(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined), true)]
    #[case::multicenter_not(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter), false)]
    #[case::multicenter_with_value(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1))), false)]
    fn test_atom_constraint_is_undetermined(
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, true)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, false)]
    #[case::aromatic_lit(AromaticValenceAst::Aromatic(ValueAst::Lit(1)), false)]
    #[case::aromatic_inner_undetermined(
        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
        false
    )]
    fn test_aromatic_valence_ast_is_undetermined(
        #[case] v: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, true)]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, false)]
    #[case::multicenter_lit(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)), false)]
    fn test_multicenter_valence_ast_is_undetermined(
        #[case] v: MulticenterValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(v.is_undetermined(), expected);
    }

    #[rstest]
    fn test_atom_constraints_retain_keeps_matching() {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(4)),
            AtomConstraint::Degree(ValueAst::Lit(3)),
            AtomConstraint::RingCount(ValueAst::Lit(2)),
        ]);
        cs.retain(|c| matches!(c, AtomConstraint::Valence(_) | AtomConstraint::RingCount(_)));
        assert_eq!(cs.len(), 2);
        assert!(cs.contains(AtomConstraintKind::Valence));
        assert!(cs.contains(AtomConstraintKind::RingCount));
        assert!(!cs.contains(AtomConstraintKind::Degree));
    }

    #[rstest]
    fn test_atom_constraints_retain_predicate_false_empties() {
        let mut cs = AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))]);
        cs.retain(|_| false);
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_atom_constraints_clear_removes_all() {
        let mut cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(4)),
            AtomConstraint::Degree(ValueAst::Lit(3)),
        ]);
        cs.clear();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[rstest]
    fn test_atom_constraints_remap_is_noop() {
        let cs = AtomConstraints::from_iter([
            AtomConstraint::Valence(ValueAst::Lit(4)),
            AtomConstraint::Degree(ValueAst::Lit(3)),
        ]);
        let remap = IdxRemapping::new(
            umol_graph_core::Remapping {
                removed_nodes: vec![0, 1, 2],
                removed_edges: vec![0],
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        let after = cs.clone().remap(&remap);
        assert_eq!(after, cs);
    }
}
