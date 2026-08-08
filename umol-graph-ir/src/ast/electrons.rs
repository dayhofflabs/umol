//! Electron-count leaf: the per-member-atom electron-count vector shared by
//! `AromaticSystemAst` and `MulticenterBondAst`.

use std::borrow::Cow;

use umol_graph_core::ParticipantPosition;

use super::error::{Contradiction, NoJoin};
use super::traits::{AsLit, Canonicalize, Lattice};

/// Per-position electron counts as one atomic lattice value: undetermined, or a
/// concrete vector. The vector is positional (cell = member atom), so it is
/// compared whole — never sorted, deduped, or matched cell-by-cell. `i64` to
/// mirror the other electron-count quantities (`charge`, the `#e` constraint).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ElectronCountsAst {
    #[default]
    Undetermined,
    Lit(Vec<i64>),
}

impl ElectronCountsAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(counts: Vec<i64>) -> Self {
        Self::Lit(counts)
    }

    /// Reorder the positional counts by `order` (`new[i] = old[order[i]]`) to
    /// track a participant reordering. `Undetermined` is unchanged. A length
    /// mismatch (a malformed count vector, rejected later by structure
    /// validation) is left untouched rather than reindexed.
    pub fn permute(&mut self, order: &[ParticipantPosition]) {
        if let Self::Lit(counts) = self {
            if order.len() != counts.len() {
                return;
            }
            let reordered: Vec<i64> = order.iter().map(|p| counts[p.index()]).collect();
            *counts = reordered;
        }
    }
}

impl From<Vec<i64>> for ElectronCountsAst {
    fn from(counts: Vec<i64>) -> Self {
        Self::Lit(counts)
    }
}

impl Canonicalize for ElectronCountsAst {
    /// Positional vector — both variants are already canonical (no sort/dedup).
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Borrowed(self))
    }
}

impl AsLit for ElectronCountsAst {
    type Lit = Vec<i64>;

    /// The concrete count vector, only when it is a literal.
    #[inline]
    fn as_lit(&self) -> Option<Vec<i64>> {
        match self {
            Self::Lit(counts) => Some(counts.clone()),
            Self::Undetermined => None,
        }
    }
}

impl Lattice for ElectronCountsAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    #[inline]
    fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// Atomic exact-match: the whole vector meets only an equal vector (both
    /// length and contents); any mismatch is `None`.
    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Lit(a), Self::Lit(b)) if a == b => Some(Self::Lit(a.clone())),
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        Ok(match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(a.clone()),
            _ => Self::Undetermined,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::triple(vec![1, 1, 1], ElectronCountsAst::Lit(vec![1, 1, 1]))]
    #[case::mixed(vec![2, 0, 2], ElectronCountsAst::Lit(vec![2, 0, 2]))]
    fn test_electron_counts_ast_from(
        #[case] counts: Vec<i64>,
        #[case] expected: ElectronCountsAst,
    ) {
        assert_eq!(ElectronCountsAst::from(counts), expected);
    }

    #[rstest]
    #[case::undetermined(ElectronCountsAst::Undetermined, ElectronCountsAst::Undetermined)]
    #[case::lit(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1, 1]))]
    fn test_electron_counts_ast_constructors(
        #[case] actual: ElectronCountsAst,
        #[case] expected: ElectronCountsAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::reorder(
        ElectronCountsAst::Lit(vec![10, 20, 30]),
        vec![ParticipantPosition(2), ParticipantPosition(0), ParticipantPosition(1)],
        ElectronCountsAst::Lit(vec![30, 10, 20]),
    )]
    #[case::undetermined(
        ElectronCountsAst::Undetermined,
        vec![ParticipantPosition(1), ParticipantPosition(0)],
        ElectronCountsAst::Undetermined,
    )]
    fn test_electron_counts_ast_permute(
        #[case] mut input: ElectronCountsAst,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: ElectronCountsAst,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rstest]
    #[case::undetermined(ElectronCountsAst::Undetermined)]
    #[case::lit(ElectronCountsAst::Lit(vec![1, 1, 1]))]
    fn test_electron_counts_ast_canonicalize_identity(#[case] input: ElectronCountsAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ElectronCountsAst::Lit(vec![2, 0, 2]), Some(vec![2, 0, 2]))]
    #[case::undetermined(ElectronCountsAst::Undetermined, None)]
    fn test_electron_counts_ast_as_lit(#[case] ast: ElectronCountsAst, #[case] expected: Option<Vec<i64>>) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::undetermined(ElectronCountsAst::Undetermined, true)]
    #[case::lit(ElectronCountsAst::Lit(vec![1, 1, 1]), false)]
    fn test_electron_counts_ast_is_undetermined(
        #[case] ast: ElectronCountsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElectronCountsAst::Undetermined, ElectronCountsAst::Lit(vec![1, 1, 1]), Some(ElectronCountsAst::Lit(vec![1, 1, 1])))]
    #[case::lit_und(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Undetermined, Some(ElectronCountsAst::Lit(vec![1, 1, 1])))]
    #[case::und_und(ElectronCountsAst::Undetermined, ElectronCountsAst::Undetermined, Some(ElectronCountsAst::Undetermined))]
    #[case::lit_lit_eq(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1, 1]), Some(ElectronCountsAst::Lit(vec![1, 1, 1])))]
    #[case::lit_lit_neq(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![2, 0, 2]), None)]
    #[case::lit_lit_len_mismatch(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1]), None)]
    fn test_electron_counts_ast_meet(
        #[case] a: ElectronCountsAst,
        #[case] b: ElectronCountsAst,
        #[case] expected: Option<ElectronCountsAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElectronCountsAst::Undetermined, ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Undetermined)]
    #[case::und_und(ElectronCountsAst::Undetermined, ElectronCountsAst::Undetermined, ElectronCountsAst::Undetermined)]
    #[case::lit_lit_eq(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1, 1]))]
    #[case::lit_lit_neq(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![2, 0, 2]), ElectronCountsAst::Undetermined)]
    #[case::lit_lit_len_mismatch(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1]), ElectronCountsAst::Undetermined)]
    fn test_electron_counts_ast_join(
        #[case] a: ElectronCountsAst,
        #[case] b: ElectronCountsAst,
        #[case] expected: ElectronCountsAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElectronCountsAst::Undetermined, ElectronCountsAst::Lit(vec![1, 1, 1]), true)]
    #[case::und_und(ElectronCountsAst::Undetermined, ElectronCountsAst::Undetermined, true)]
    #[case::lit_und(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Undetermined, false)]
    #[case::lit_lit_match(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1, 1]), true)]
    #[case::lit_lit_mismatch(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![2, 0, 2]), false)]
    #[case::lit_lit_len_mismatch(ElectronCountsAst::Lit(vec![1, 1, 1]), ElectronCountsAst::Lit(vec![1, 1]), false)]
    fn test_electron_counts_ast_matches(
        #[case] pattern: ElectronCountsAst,
        #[case] target: ElectronCountsAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
