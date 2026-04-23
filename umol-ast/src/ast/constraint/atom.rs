//! Per-atom constraints.

use strum::{EnumCount, EnumDiscriminants};

use super::super::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash, EnumCount))]
#[repr(u8)]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceConstraint),
    MulticenterValence(MulticenterValenceConstraint),
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
pub enum AromaticValenceConstraint {
    #[default]
    Undetermined,
    NotAromatic,
    Aromatic(ValueAst),
}

impl AromaticValenceConstraint {
    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum MulticenterValenceConstraint {
    #[default]
    Undetermined,
    NotMulticenter,
    Multicenter(ValueAst),
}

impl MulticenterValenceConstraint {
    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }
}

/// Per-atom constraint slotmap. Fixed-size array indexed by
/// [`AtomConstraintKind`]; each slot holds at most one constraint of that
/// kind. O(1) access and update by kind.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomConstraints {
    slots: [Option<AtomConstraint>; AtomConstraintKind::COUNT],
}

impl Default for AtomConstraints {
    fn default() -> Self {
        Self {
            slots: [const { None }; AtomConstraintKind::COUNT],
        }
    }
}

impl AtomConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.iter().all(Option::is_none)
    }

    pub fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.is_some()).count()
    }

    pub fn contains(&self, kind: AtomConstraintKind) -> bool {
        self.slots[kind as usize].is_some()
    }

    pub fn get(&self, kind: AtomConstraintKind) -> Option<&AtomConstraint> {
        self.slots[kind as usize].as_ref()
    }

    pub fn get_mut(&mut self, kind: AtomConstraintKind) -> Option<&mut AtomConstraint> {
        self.slots[kind as usize].as_mut()
    }

    /// Insert a constraint in its kind's slot, returning the previous occupant.
    pub fn set(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        let slot = &mut self.slots[constraint.kind() as usize];
        slot.replace(constraint)
    }

    pub fn remove(&mut self, kind: AtomConstraintKind) -> Option<AtomConstraint> {
        self.slots[kind as usize].take()
    }

    pub fn remove_undetermined(&mut self) {
        let kinds_to_strip: Vec<AtomConstraintKind> = self
            .iter()
            .filter(|c| c.is_undetermined())
            .map(|c| c.kind())
            .collect();
        for k in kinds_to_strip {
            self.remove(k);
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.slots.iter().filter_map(Option::as_ref)
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut AtomConstraint> {
        self.slots.iter_mut().filter_map(Option::as_mut)
    }

}

impl FromIterator<AtomConstraint> for AtomConstraints {
    fn from_iter<I: IntoIterator<Item = AtomConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
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
    #[case::aromatic(AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic), AtomConstraintKind::AromaticValence)]
    #[case::multicenter(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Undetermined), AtomConstraintKind::MulticenterValence)]
    #[case::donated(AtomConstraint::DonatedPairs(ValueAst::Lit(1)), AtomConstraintKind::DonatedPairs)]
    #[case::accepted(AtomConstraint::AcceptedPairs(ValueAst::Lit(2)), AtomConstraintKind::AcceptedPairs)]
    #[case::degree(AtomConstraint::Degree(ValueAst::Lit(3)), AtomConstraintKind::Degree)]
    #[case::connectivity(AtomConstraint::Connectivity(ValueAst::Lit(4)), AtomConstraintKind::Connectivity)]
    #[case::ring_conn(AtomConstraint::RingConnectivity(ValueAst::Lit(2)), AtomConstraintKind::RingConnectivity)]
    #[case::total_h(AtomConstraint::TotalHydrogens(ValueAst::Lit(3)), AtomConstraintKind::TotalHydrogens)]
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
        let prev = cs.set(AtomConstraint::Valence(ValueAst::Lit(4)));
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(3)));
        let prev = cs.set(AtomConstraint::Valence(ValueAst::Lit(4)));
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.set(AtomConstraint::Degree(ValueAst::Lit(3)));
        cs.set(AtomConstraint::AromaticValence(
            AromaticValenceConstraint::NotAromatic,
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
                AromaticValenceConstraint::NotAromatic
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(4)));
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(3)));
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(4)));
        cs.set(AtomConstraint::RingSize(ValueAst::Lit(6)));
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
        cs.set(AtomConstraint::Valence(ValueAst::Lit(3)));
        cs.set(AtomConstraint::Degree(ValueAst::Lit(2)));
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
}
