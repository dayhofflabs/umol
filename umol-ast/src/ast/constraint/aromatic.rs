//! Per-aromatic-system constraints.

use std::mem;
use std::slice::Iter;

use super::atom::AtomConstraint;
use crate::ast::idx::AtomIdx;
use crate::ast::remap::IdxRemapping;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraint {
    Atoms(Vec<AtomIdx>),
    Contains(AtomIdx),
    ContainsAll(Vec<AtomIdx>),
    AllAtoms(Box<AtomConstraint>),
    AnyAtom(Box<AtomConstraint>),
}

impl AromaticSystemConstraint {
    /// Single-valued per aromatic system: `Atoms` (defines the atom set
    /// explicitly). Multi-valued: `Contains`, `ContainsAll`, `AllAtoms`,
    /// `AnyAtom`.
    pub fn is_unique(&self) -> bool {
        matches!(self, Self::Atoms(_))
    }

    /// Topology references (`Atoms`, `Contains`, `ContainsAll`) are never
    /// undetermined. `AllAtoms` / `AnyAtom` delegate to the inner atom
    /// constraint.
    pub fn is_undetermined(&self) -> bool {
        match self {
            Self::Atoms(_) | Self::Contains(_) | Self::ContainsAll(_) => false,
            Self::AllAtoms(c) | Self::AnyAtom(c) => c.is_undetermined(),
        }
    }

    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::Atoms(atoms) => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(Self::Atoms)
            }
            Self::Contains(a) => remap.atom(a).map(Self::Contains),
            Self::ContainsAll(atoms) => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(Self::ContainsAll)
            }
            Self::AllAtoms(c) => Some(Self::AllAtoms(c)),
            Self::AnyAtom(c) => Some(Self::AnyAtom(c)),
        }
    }
}

/// Per-aromatic-system constraint container. Enforces the per-variant
/// cardinality policy in [`AromaticSystemConstraint::is_unique`] on insert.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticSystemConstraints(Vec<AromaticSystemConstraint>);

impl AromaticSystemConstraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[AromaticSystemConstraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, AromaticSystemConstraint> {
        self.0.iter()
    }

    pub fn add(&mut self, c: AromaticSystemConstraint) -> Option<AromaticSystemConstraint> {
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

    pub fn retain(&mut self, mut f: impl FnMut(&AromaticSystemConstraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn remap(self, remap: &IdxRemapping) -> Self {
        Self(self.0.into_iter().filter_map(|c| c.remap(remap)).collect())
    }
}

impl FromIterator<AromaticSystemConstraint> for AromaticSystemConstraints {
    fn from_iter<I: IntoIterator<Item = AromaticSystemConstraint>>(iter: I) -> Self {
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

    #[rustfmt::skip]
    #[rstest]
    #[case::atoms(AromaticSystemConstraint::Atoms(vec![AtomIdx(0)]), true)]
    #[case::contains(AromaticSystemConstraint::Contains(AtomIdx(0)), false)]
    #[case::contains_all(AromaticSystemConstraint::ContainsAll(vec![AtomIdx(0), AtomIdx(1)]), false)]
    #[case::all_atoms(AromaticSystemConstraint::AllAtoms(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), false)]
    #[case::any_atom(AromaticSystemConstraint::AnyAtom(Box::new(AtomConstraint::Degree(ValueAst::Lit(2)))), false)]
    fn test_aromatic_system_constraint_is_unique(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_unique(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atoms(AromaticSystemConstraint::Atoms(vec![AtomIdx(0)]), false)]
    #[case::contains(AromaticSystemConstraint::Contains(AtomIdx(0)), false)]
    #[case::contains_all(AromaticSystemConstraint::ContainsAll(vec![AtomIdx(0)]), false)]
    #[case::all_atoms_lit(AromaticSystemConstraint::AllAtoms(Box::new(AtomConstraint::Valence(ValueAst::Lit(3)))), false)]
    #[case::all_atoms_undetermined(AromaticSystemConstraint::AllAtoms(Box::new(AtomConstraint::Valence(ValueAst::Undetermined))), true)]
    #[case::any_atom_undetermined(AromaticSystemConstraint::AnyAtom(Box::new(AtomConstraint::Degree(ValueAst::Undetermined))), true)]
    fn test_aromatic_system_constraint_is_undetermined(
        #[case] c: AromaticSystemConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_undetermined(), expected);
    }

    #[rstest]
    fn test_aromatic_system_constraints_add_unique_replaces() {
        let mut cs = AromaticSystemConstraints::new();
        cs.add(AromaticSystemConstraint::Atoms(vec![
            AtomIdx(0),
            AtomIdx(1),
        ]));
        let prev = cs.add(AromaticSystemConstraint::Atoms(vec![
            AtomIdx(2),
            AtomIdx(3),
        ]));
        assert_eq!(
            prev,
            Some(AromaticSystemConstraint::Atoms(vec![
                AtomIdx(0),
                AtomIdx(1)
            ]))
        );
        assert_eq!(
            cs.as_slice(),
            &[AromaticSystemConstraint::Atoms(vec![
                AtomIdx(2),
                AtomIdx(3)
            ])]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_add_multi_appends() {
        let mut cs = AromaticSystemConstraints::new();
        cs.add(AromaticSystemConstraint::Contains(AtomIdx(0)));
        cs.add(AromaticSystemConstraint::Contains(AtomIdx(1)));
        assert_eq!(
            cs.as_slice(),
            &[
                AromaticSystemConstraint::Contains(AtomIdx(0)),
                AromaticSystemConstraint::Contains(AtomIdx(1)),
            ]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_retain() {
        let mut cs = AromaticSystemConstraints::from_iter([
            AromaticSystemConstraint::Contains(AtomIdx(0)),
            AromaticSystemConstraint::Atoms(vec![AtomIdx(1)]),
        ]);
        cs.retain(|c| matches!(c, AromaticSystemConstraint::Atoms(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_aromatic_system_constraints_clear() {
        let mut cs =
            AromaticSystemConstraints::from_iter([AromaticSystemConstraint::Contains(AtomIdx(0))]);
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap_shifts_atom_refs() {
        let cs = AromaticSystemConstraints::from_iter([
            AromaticSystemConstraint::Contains(AtomIdx(3)),
            AromaticSystemConstraint::Atoms(vec![AtomIdx(2), AtomIdx(4)]),
        ]);
        let remap = idx_remapping(vec![1]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[
                AromaticSystemConstraint::Contains(AtomIdx(2)),
                AromaticSystemConstraint::Atoms(vec![AtomIdx(1), AtomIdx(3)]),
            ]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap_drops_removed_atom() {
        let cs = AromaticSystemConstraints::from_iter([
            AromaticSystemConstraint::Contains(AtomIdx(1)),
            AromaticSystemConstraint::Atoms(vec![AtomIdx(0)]),
        ]);
        let remap = idx_remapping(vec![1]);
        let after = cs.remap(&remap);
        assert_eq!(
            after.as_slice(),
            &[AromaticSystemConstraint::Atoms(vec![AtomIdx(0)])]
        );
    }

    #[rstest]
    fn test_aromatic_system_constraints_remap_drops_atoms_if_any_removed() {
        let cs = AromaticSystemConstraints::from_iter([AromaticSystemConstraint::Atoms(vec![
            AtomIdx(0),
            AtomIdx(3),
        ])]);
        let remap = idx_remapping(vec![3]);
        let after = cs.remap(&remap);
        assert!(after.is_empty());
    }
}
