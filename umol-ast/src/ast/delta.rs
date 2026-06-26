//! Resolved edit vocabulary: the `Delta` counterpart of the deferred `Edit`.
//!
//! A `Delta` is one resolved edit over a `MoleculeAst`, referencing entities by stable
//! ids in the molecule's own frame (no positional `New`). The vocabulary is closed
//! under inversion — every delta's inverse is another delta — so it needs no `Undo`
//! journal. Molecule-level and reuse-agnostic: reactions, base+delta storage, and
//! matched-pair transforms all build on it. Increment 1 covers the localized-topology
//! families (atoms, bonds, and constraints — inline per-entity and molecule-level);
//! overlay and stereo families follow.

use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{AtomConstraint, BondConstraint, Constraint};
use super::edit::{AtomFieldChange, BondFieldChange};
use super::id::{AtomId, BondId};

/// A resolved edit to a single atom. `SetConstraint` is a keyed old→new change of an
/// inline constraint: `(None, Some)` adds, `(Some, None)` removes, `(Some, Some)`
/// modifies (old and new sharing a key).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomDelta {
    Add {
        id: AtomId,
        ast: AtomAst,
    },
    Remove {
        id: AtomId,
        ast: AtomAst,
    },
    SetField {
        id: AtomId,
        change: AtomFieldChange,
    },
    SetConstraint {
        id: AtomId,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    },
}

impl AtomDelta {
    /// The inverse delta: `Add`↔`Remove`; `SetField` / `SetConstraint` swap old/new.
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, ast } => Self::Remove { id, ast },
            Self::Remove { id, ast } => Self::Add { id, ast },
            Self::SetField { id, change } => Self::SetField {
                id,
                change: change.inverse(),
            },
            Self::SetConstraint { id, old, new } => Self::SetConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single bond. Identified by `id` (the uniform per-family
/// identity); `Add`/`Remove` also carry `endpoints` — the structural payload (which
/// atoms) needed to create/restore the bond and invert without the host.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondDelta {
    Add {
        id: BondId,
        endpoints: [AtomId; 2],
        ast: BondAst,
    },
    Remove {
        id: BondId,
        endpoints: [AtomId; 2],
        ast: BondAst,
    },
    SetField {
        id: BondId,
        change: BondFieldChange,
    },
    SetConstraint {
        id: BondId,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    },
}

impl BondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, endpoints, ast } => Self::Remove { id, endpoints, ast },
            Self::Remove { id, endpoints, ast } => Self::Add { id, endpoints, ast },
            Self::SetField { id, change } => Self::SetField {
                id,
                change: change.inverse(),
            },
            Self::SetConstraint { id, old, new } => Self::SetConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved change to the molecule-level constraint set, as a set-diff.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConstraintDelta {
    Add(Constraint),
    Remove(Constraint),
}

impl ConstraintDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add(constraint) => Self::Remove(constraint),
            Self::Remove(constraint) => Self::Add(constraint),
        }
    }
}

/// One resolved edit across the localized-topology families.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Delta {
    Atom(AtomDelta),
    Bond(BondDelta),
    Constraint(ConstraintDelta),
}

impl Delta {
    /// The inverse delta. The vocabulary is closed under inversion, so this is total
    /// and lands back in `Delta`.
    pub fn inverse(self) -> Self {
        match self {
            Self::Atom(delta) => Self::Atom(delta.inverse()),
            Self::Bond(delta) => Self::Bond(delta.inverse()),
            Self::Constraint(delta) => Self::Constraint(delta.inverse()),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::constraint::MoleculeConstraint;
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    #[case::add_remove(
        AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) },
        AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) }
    )]
    #[case::set_field(
        AtomDelta::SetField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(1) },
        },
        AtomDelta::SetField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(1), new: ValueAst::Lit(0) },
        }
    )]
    #[case::set_constraint(
        AtomDelta::SetConstraint {
            id: AtomId(2),
            old: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
            new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
        },
        AtomDelta::SetConstraint {
            id: AtomId(2),
            old: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
            new: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
        }
    )]
    fn test_atom_delta_inverse(#[case] input: AtomDelta, #[case] expected: AtomDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::add_remove(
        BondDelta::Add {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        },
        BondDelta::Remove {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }
    )]
    #[case::set_field(
        BondDelta::SetField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
        },
        BondDelta::SetField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(1) },
        }
    )]
    #[case::set_constraint(
        BondDelta::SetConstraint {
            id: BondId(3),
            old: None,
            new: Some(BondConstraint::Aromatic),
        },
        BondDelta::SetConstraint {
            id: BondId(3),
            old: Some(BondConstraint::Aromatic),
            new: None,
        }
    )]
    fn test_bond_delta_inverse(#[case] input: BondDelta, #[case] expected: BondDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_constraint_delta_inverse() {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: ValueAst::Lit(0),
        });
        assert_eq!(
            ConstraintDelta::Add(constraint.clone()).inverse(),
            ConstraintDelta::Remove(constraint.clone()),
        );
        assert_eq!(
            ConstraintDelta::Add(constraint.clone()).inverse().inverse(),
            ConstraintDelta::Add(constraint),
        );
    }

    #[rstest]
    #[case::atom(
        Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
        Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) })
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Add {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }),
        Delta::Bond(BondDelta::Remove {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        })
    )]
    fn test_delta_inverse(#[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }
}
