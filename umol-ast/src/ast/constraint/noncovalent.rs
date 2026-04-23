//! Per-noncovalent-bond constraints.

use std::mem;
use std::slice::Iter;

use super::atom::AtomConstraint;
use crate::ast::idx::AtomIdx;
use crate::ast::remap::IdxRemapping;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraint {
    Ends([AtomIdx; 2]),
    Contains(AtomIdx),
    EndsSatisfy([Box<AtomConstraint>; 2]),
}

impl NoncovalentBondConstraint {
    /// Single-valued per noncovalent bond: `Ends` (the two endpoints of the
    /// interaction). Multi-valued: `Contains`, `EndsSatisfy` (each is an
    /// independent endpoint-predicate pair that AND together with others).
    pub fn is_unique(&self) -> bool {
        matches!(self, Self::Ends(_))
    }

    /// Topology references (`Ends`, `Contains`) are never undetermined.
    /// `EndsSatisfy` is undetermined iff both endpoint constraints are.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Ends(_) | Self::Contains(_) => false,
            Self::EndsSatisfy([a, b]) => a.is_undetermined() && b.is_undetermined(),
        }
    }

    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::Ends([a, b]) => {
                let a = remap.atom(a)?;
                let b = remap.atom(b)?;
                Some(Self::Ends([a, b]))
            }
            Self::Contains(a) => remap.atom(a).map(Self::Contains),
            Self::EndsSatisfy(cs) => Some(Self::EndsSatisfy(cs)),
        }
    }
}

/// Per-noncovalent-bond constraint container. Enforces the per-variant
/// cardinality policy in [`NoncovalentBondConstraint::is_unique`] on insert.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBondConstraints(Vec<NoncovalentBondConstraint>);

impl NoncovalentBondConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[NoncovalentBondConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, NoncovalentBondConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: NoncovalentBondConstraint) -> Option<NoncovalentBondConstraint> {
        if c.is_unique() {
            if let Some(pos) = self
                .0
                .iter()
                .position(|e| mem::discriminant(e) == mem::discriminant(&c))
            {
                return Some(mem::replace(&mut self.0[pos], c));
            }
        }
        self.0.push(c);
        None
    }

    pub fn retain(&mut self, mut f: impl FnMut(&NoncovalentBondConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn remap(self, remap: &IdxRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl FromIterator<NoncovalentBondConstraint> for NoncovalentBondConstraints {
    fn from_iter<I: IntoIterator<Item = NoncovalentBondConstraint>>(iter: I) -> Self {
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
    use crate::ast::value::ValueAst;

    fn idx_remapping(removed_nodes: Vec<u32>) -> IdxRemapping {
        IdxRemapping::new(
            umol_graph_core::Remapping {
                removed_nodes,
                removed_edges: Vec::new(),
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    fn atom_v(n: i64) -> Box<AtomConstraint> {
        Box::new(AtomConstraint::Valence(ValueAst::Lit(n)))
    }

    fn atom_undet() -> Box<AtomConstraint> {
        Box::new(AtomConstraint::Valence(ValueAst::Undetermined))
    }

    #[rstest]
    #[case::ends(NoncovalentBondConstraint::Ends([AtomIdx(0), AtomIdx(1)]), true)]
    #[case::contains(NoncovalentBondConstraint::Contains(AtomIdx(0)), false)]
    #[case::ends_satisfy(NoncovalentBondConstraint::EndsSatisfy([atom_v(3), atom_v(4)]), false)]
    fn test_noncovalent_bond_constraint_is_unique(
        #[case] c: NoncovalentBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rstest]
    #[case::ends(NoncovalentBondConstraint::Ends([AtomIdx(0), AtomIdx(1)]), false)]
    #[case::contains(NoncovalentBondConstraint::Contains(AtomIdx(0)), false)]
    #[case::ends_satisfy_both_lit(NoncovalentBondConstraint::EndsSatisfy([atom_v(3), atom_v(4)]), false)]
    #[case::ends_satisfy_one_undetermined(NoncovalentBondConstraint::EndsSatisfy([atom_undet(), atom_v(4)]), false)]
    #[case::ends_satisfy_both_undetermined(NoncovalentBondConstraint::EndsSatisfy([atom_undet(), atom_undet()]), true)]
    fn test_noncovalent_bond_constraint_is_undetermined(
        #[case] c: NoncovalentBondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_add_unique_replaces() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.add(NoncovalentBondConstraint::Ends([AtomIdx(0), AtomIdx(1)]));
        let prev = cs.add(NoncovalentBondConstraint::Ends([AtomIdx(2), AtomIdx(3)]));
        assert_eq!(
            prev,
            Some(NoncovalentBondConstraint::Ends([AtomIdx(0), AtomIdx(1)]))
        );
        assert_eq!(
            cs.as_slice(),
            &[NoncovalentBondConstraint::Ends([AtomIdx(2), AtomIdx(3)])]
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_add_multi_appends() {
        let mut cs = NoncovalentBondConstraints::new();
        cs.add(NoncovalentBondConstraint::Contains(AtomIdx(0)));
        cs.add(NoncovalentBondConstraint::Contains(AtomIdx(1)));
        cs.add(NoncovalentBondConstraint::EndsSatisfy([
            atom_v(3),
            atom_v(4),
        ]));
        cs.add(NoncovalentBondConstraint::EndsSatisfy([
            atom_v(5),
            atom_v(6),
        ]));
        assert_eq!(cs.len(), 4);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_retain() {
        let mut cs = NoncovalentBondConstraints::from_iter([
            NoncovalentBondConstraint::Contains(AtomIdx(0)),
            NoncovalentBondConstraint::Ends([AtomIdx(1), AtomIdx(2)]),
        ]);
        cs.retain(|c| matches!(c, NoncovalentBondConstraint::Ends(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_clear() {
        let mut cs = NoncovalentBondConstraints::from_iter([NoncovalentBondConstraint::Contains(
            AtomIdx(0),
        )]);
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_remap_shifts_atom_refs() {
        let cs = NoncovalentBondConstraints::from_iter([
            NoncovalentBondConstraint::Ends([AtomIdx(2), AtomIdx(3)]),
            NoncovalentBondConstraint::Contains(AtomIdx(4)),
        ]);
        let remap = idx_remapping(vec![1]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[
                NoncovalentBondConstraint::Ends([AtomIdx(1), AtomIdx(2)]),
                NoncovalentBondConstraint::Contains(AtomIdx(3)),
            ]
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_remap_drops_when_ends_removed() {
        let cs = NoncovalentBondConstraints::from_iter([
            NoncovalentBondConstraint::Ends([AtomIdx(1), AtomIdx(3)]),
            NoncovalentBondConstraint::Contains(AtomIdx(0)),
        ]);
        let remap = idx_remapping(vec![1]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[NoncovalentBondConstraint::Contains(AtomIdx(0))]
        );
    }
}
