//! Reaction AST: a left-hand-side molecule plus a resolved transformation (`Deltas`).
//!
//! Homoiconic — a molecule is the empty-deltas case, a rule is a pattern `lhs` plus
//! deltas, and applying a rule yields a concrete reaction of the same type. The atom
//! map, R-side, condensed (CGR) form, and reverse reaction are all *derived* from
//! `(lhs, deltas)` rather than stored (those derivations land with delta lowering).

use super::delta::Deltas;
use super::error::Contradiction;
use super::molecule::MoleculeAst;
use super::traits::Canonicalize;

/// A reaction as one full molecule state (`lhs`) plus one resolved delta (`deltas`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionAst {
    pub lhs: MoleculeAst,
    pub deltas: Deltas,
}

impl ReactionAst {
    pub fn new(lhs: MoleculeAst, deltas: Deltas) -> Self {
        Self { lhs, deltas }
    }
}

impl Canonicalize for ReactionAst {
    /// Value-level in a fixed atom frame: `deltas` are canonicalized (the #2 reduction);
    /// `lhs` is passed through (`MoleculeAst` has no whole-molecule canonical form — its
    /// equality is structural). Equality up to atom renumbering is a separate `umol-graph`
    /// operation.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(Self {
            lhs: self.lhs,
            deltas: self.deltas.canonicalize()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::super::delta::{AtomDelta, Delta};
    use super::super::edit::AtomFieldChange;
    use super::super::id::AtomId;
    use super::super::value::ValueAst;
    use super::*;

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::SetField {
            id: AtomId(id),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(old),
                new: ValueAst::Lit(new),
            },
        })
    }

    #[rstest]
    fn test_reaction_ast_canonicalize() {
        // The delta chain fuses; the lhs is passed through unchanged.
        let reaction = ReactionAst::new(
            MoleculeAst::default(),
            Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)]),
        );
        assert_eq!(
            reaction.canonicalize().unwrap(),
            ReactionAst::new(
                MoleculeAst::default(),
                Deltas::from_iter([charge_set(0, 0, 2)]),
            ),
        );
    }
}
