//! Noncovalent bond AST.

use super::constraint::NoncovalentBondConstraints;

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBondAst {
    pub kind: NoncovalentBondKindAst,
    pub constraints: NoncovalentBondConstraints,
}

impl NoncovalentBondAst {
    pub fn new(kind: NoncovalentBondKindAst) -> Self {
        Self {
            kind,
            constraints: NoncovalentBondConstraints::new(),
        }
    }

    pub fn from_kind(kind: NoncovalentBondKind) -> Self {
        Self::new(NoncovalentBondKindAst::Lit(kind))
    }

    pub fn is_ground(&self) -> bool {
        self.kind.is_ground()
    }

    pub fn matches(&self, target: &NoncovalentBondAst) -> bool {
        self.kind.matches(&target.kind)
    }

    /// Simplify every constraint's inner value in place. `kind` carries no
    /// `ValueAst`, so it is unchanged.
    pub fn simplify_values(&mut self) {
        self.constraints.simplify_each();
    }
}

/// Noncovalent interaction kind expressions. Mirrors `ElementAst`:
/// wildcard, literal, set, bind, ref.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKindAst {
    #[default]
    Undetermined,
    Lit(NoncovalentBondKind),
    Set(Vec<NoncovalentBondKind>),
    Bind {
        id: String,
        set: Vec<NoncovalentBondKind>,
    },
    Ref(String),
}

impl NoncovalentBondKindAst {
    pub fn new(kind: NoncovalentBondKind) -> Self {
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

/// Fundamental kind of a noncovalent interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKind {
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
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, false)]
    #[case::set(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), false)]
    #[case::bind(NoncovalentBondKindAst::Bind { id: "k".into(), set: vec![NoncovalentBondKind::HydrogenBond] }, false)]
    #[case::reference(NoncovalentBondKindAst::Ref("k".into()), false)]
    fn test_noncovalent_kind_ast_is_ground(
        #[case] ast: NoncovalentBondKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined_undetermined(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined, true)]
    #[case::lit_undetermined(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Undetermined, false)]
    #[case::lit_lit_match(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::lit_lit_mismatch(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), false)]
    #[case::set_lit_in(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), true)]
    #[case::set_lit_out(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond]), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), false)]
    #[case::set_set_subset(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic, NoncovalentBondKind::VanDerWaals]),
        NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), true)]
    #[case::set_set_superset(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond]),
        NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), false)]
    #[case::bind_lit_match(NoncovalentBondKindAst::Bind { id: "k".into(), set: vec![NoncovalentBondKind::HydrogenBond] },
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::ref_lit(NoncovalentBondKindAst::Ref("k".into()), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), false)]
    fn test_noncovalent_kind_ast_matches(
        #[case] pattern: NoncovalentBondKindAst,
        #[case] target: NoncovalentBondKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(NoncovalentBondKindAst::new(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_kind_ast_new(
        #[case] actual: NoncovalentBondKindAst,
        #[case] expected: NoncovalentBondKindAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default(), false)]
    #[case::ground_lit(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    fn test_noncovalent_bond_ast_is_ground(
        #[case] ast: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondAst::default(), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic), false)]
    fn test_noncovalent_bond_ast_matches(
        #[case] pattern: NoncovalentBondAst,
        #[case] target: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
