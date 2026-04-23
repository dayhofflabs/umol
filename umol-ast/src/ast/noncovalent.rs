//! Noncovalent bond AST.

use super::constraint::NoncovalentBondConstraints;

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBondAst {
    pub kind: NoncovalentKindAst,
    pub constraints: NoncovalentBondConstraints,
}

impl NoncovalentBondAst {
    pub fn new(kind: NoncovalentKindAst) -> Self {
        Self {
            kind,
            constraints: NoncovalentBondConstraints::new(),
        }
    }

    pub fn from_kind(kind: NoncovalentKind) -> Self {
        Self::new(NoncovalentKindAst::Lit(kind))
    }

    pub fn is_ground(&self) -> bool {
        self.kind.is_ground()
    }

    pub fn matches(&self, target: &NoncovalentBondAst) -> bool {
        self.kind.matches(&target.kind)
    }
}

/// Noncovalent interaction kind expressions. Mirrors `ElementAst`:
/// wildcard, literal, set, bind, ref.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum NoncovalentKindAst {
    #[default]
    Undetermined,
    Lit(NoncovalentKind),
    Set(Vec<NoncovalentKind>),
    Bind {
        id: String,
        set: Vec<NoncovalentKind>,
    },
    Ref(String),
}

impl NoncovalentKindAst {
    pub fn new(kind: NoncovalentKind) -> Self {
        Self::Lit(kind)
    }

    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Lit(p), Self::Set(ts) | Self::Bind { set: ts, .. }) => ts.iter().all(|t| t == p),
            (Self::Set(ps) | Self::Bind { set: ps, .. }, Self::Lit(t)) => ps.contains(t),
            (
                Self::Set(ps) | Self::Bind { set: ps, .. },
                Self::Set(ts) | Self::Bind { set: ts, .. },
            ) => ts.iter().all(|t| ps.contains(t)),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentKind {
    HydrogenBond,
    HalogenBond,
    ChalcogenBond,
    Ionic,
    VanDerWaals,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), true)]
    #[case::undetermined(NoncovalentKindAst::Undetermined, false)]
    #[case::set(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]), false)]
    #[case::bind(NoncovalentKindAst::Bind { id: "k".into(), set: vec![NoncovalentKind::HydrogenBond] }, false)]
    #[case::reference(NoncovalentKindAst::Ref("k".into()), false)]
    fn test_noncovalent_kind_ast_is_ground(
        #[case] ast: NoncovalentKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(NoncovalentKindAst::Undetermined, NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), true)]
    #[case::undetermined_undetermined(NoncovalentKindAst::Undetermined, NoncovalentKindAst::Undetermined, true)]
    #[case::lit_undetermined(NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), NoncovalentKindAst::Undetermined, false)]
    #[case::lit_lit_match(NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), true)]
    #[case::lit_lit_mismatch(NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), NoncovalentKindAst::Lit(NoncovalentKind::Ionic), false)]
    #[case::set_lit_in(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]), NoncovalentKindAst::Lit(NoncovalentKind::Ionic), true)]
    #[case::set_lit_out(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond]), NoncovalentKindAst::Lit(NoncovalentKind::Ionic), false)]
    #[case::set_set_subset(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic, NoncovalentKind::VanDerWaals]),
        NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]), true)]
    #[case::set_set_superset(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond]),
        NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]), false)]
    #[case::bind_lit_match(NoncovalentKindAst::Bind { id: "k".into(), set: vec![NoncovalentKind::HydrogenBond] },
        NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), true)]
    #[case::ref_lit(NoncovalentKindAst::Ref("k".into()), NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond), false)]
    fn test_noncovalent_kind_ast_matches(
        #[case] pattern: NoncovalentKindAst,
        #[case] target: NoncovalentKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(NoncovalentKindAst::new(NoncovalentKind::HydrogenBond), NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond))]
    fn test_noncovalent_kind_ast_new(
        #[case] actual: NoncovalentKindAst,
        #[case] expected: NoncovalentKindAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default(), false)]
    #[case::ground_lit(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond), true)]
    fn test_noncovalent_bond_ast_is_ground(
        #[case] ast: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondAst::default(), NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentKind::Ionic), false)]
    fn test_noncovalent_bond_ast_matches(
        #[case] pattern: NoncovalentBondAst,
        #[case] target: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
