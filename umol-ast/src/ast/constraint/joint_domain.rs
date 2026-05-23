//! Joint-domain (relational) constraint over a tuple of atom-level fields.
use super::super::error::{Contradiction, JointDomainError};
use super::super::traits::Lattice;
use super::super::value::ValueAst;

/// Atom-level variable referenceable from a joint-domain constraint.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum JointVar {
    Charge,
    ImplicitHydrogens,
    LonePairs,
    UnpairedElectrons,
    Multiplicity,
    Valence,
    DonatedPairs,
    AcceptedPairs,
}

/// Concrete value occupying one slot of a `JointDomainAst` tuple. Only `Int`
/// today; new variants will land when `JointVar` gains `Element`,
/// `Isotope`, etc.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum JointValue {
    Int(i64),
}

impl From<i64> for JointValue {
    fn from(n: i64) -> Self {
        Self::Int(n)
    }
}

/// Relational constraint: `Undetermined` is the lattice top (no constraint);
/// `Domain(...)` asserts a finite set of admissible tuples.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum JointDomainAst {
    #[default]
    Undetermined,
    Domain(DomainState),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DomainState {
    vars: Vec<JointVar>,
    tuples: Vec<Vec<JointValue>>,
}

impl JointDomainAst {
    /// Construct a `Domain` from an integer-only tuple list. Validates
    /// arity and rejects empty vars / empty tuples / duplicate vars.
    /// Canonicalizes by sorting vars (and permuting each tuple to match),
    /// sorting tuples lexicographically, and dedup'ing the tuple list.
    pub fn from_ints(
        vars: Vec<JointVar>,
        tuples: Vec<Vec<i64>>,
    ) -> Result<Self, JointDomainError> {
        if vars.is_empty() {
            return Err(JointDomainError::TooFewVars(0));
        }
        let vars_len = vars.len();
        for (i, tuple) in tuples.iter().enumerate() {
            if tuple.len() != vars_len {
                return Err(JointDomainError::ArityMismatch {
                    tuple_index: i,
                    tuple_len: tuple.len(),
                    vars_len,
                });
            }
        }
        let permutation = sort_permutation(&vars)?;
        let sorted_vars: Vec<JointVar> = permutation.iter().map(|&i| vars[i]).collect();
        let mut sorted_tuples: Vec<Vec<JointValue>> = tuples
            .into_iter()
            .map(|tuple| {
                permutation
                    .iter()
                    .map(|&i| JointValue::Int(tuple[i]))
                    .collect()
            })
            .collect();
        sorted_tuples.sort();
        sorted_tuples.dedup();
        if sorted_tuples.is_empty() {
            return Err(JointDomainError::TooFewTuples(0));
        }
        Ok(Self::Domain(DomainState {
            vars: sorted_vars,
            tuples: sorted_tuples,
        }))
    }

    /// Variables this constraint binds, in canonical sorted order. `None`
    /// for the `Undetermined` top.
    pub fn vars(&self) -> Option<&[JointVar]> {
        match self {
            Self::Undetermined => None,
            Self::Domain(d) => Some(&d.vars),
        }
    }

    /// Admissible tuples over `vars`, in canonical sorted+dedup'd order.
    /// `None` for the `Undetermined` top.
    pub fn tuples(&self) -> Option<&[Vec<JointValue>]> {
        match self {
            Self::Undetermined => None,
            Self::Domain(d) => Some(&d.tuples),
        }
    }

    /// Re-normalize internal state of the `Domain` variant: sort and dedup
    /// the tuple list. Idempotent on values produced by `from_ints`; serves
    /// as a safety net for code paths that build a `Domain` via internal
    /// constructors (e.g., the relational meet path).
    pub fn simplify(mut self) -> Self {
        if let Self::Domain(d) = &mut self {
            d.tuples.sort();
            d.tuples.dedup();
        }
        self
    }

    /// Prune `Domain` tuples to those consistent with the current value of
    /// each var as supplied by `value_of`. `Err(Contradiction)` if pruning
    /// empties the tuple list. `Undetermined` returns itself unchanged.
    /// Used by per-entity `saturate` to project a JointDomain against the
    /// surrounding field state.
    pub fn project<F>(&self, value_of: F) -> Result<Self, Contradiction>
    where
        F: Fn(JointVar) -> ValueAst,
    {
        let Self::Domain(state) = self else {
            return Ok(self.clone());
        };
        let mut kept: Vec<Vec<JointValue>> = Vec::new();
        for tuple in &state.tuples {
            let mut admissible = true;
            for (var, value) in state.vars.iter().zip(tuple) {
                let JointValue::Int(n) = value;
                let field = value_of(*var);
                if field.meet(&ValueAst::Lit(*n)).is_none() {
                    admissible = false;
                    break;
                }
            }
            if admissible {
                kept.push(tuple.clone());
            }
        }
        if kept.is_empty() {
            return Err(Contradiction);
        }
        Ok(Self::Domain(DomainState {
            vars: state.vars.clone(),
            tuples: kept,
        }))
    }
}

impl Lattice for JointDomainAst {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// `Domain` is ground iff its tuple list has exactly one tuple — every
    /// var resolves to a single concrete value. `Undetermined` (top) is
    /// not ground.
    fn is_ground(&self) -> bool {
        matches!(self, Self::Domain(d) if d.tuples.len() == 1)
    }

    /// Relational meet (natural join): cartesian product on disjoint var
    /// sets, equijoin on shared, intersection on identical. `Undetermined`
    /// is the lattice top. Returns `None` only when the joined tuple set
    /// is empty (genuine contradiction).
    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Domain(s), Self::Domain(o)) => relational_meet(s, o),
        }
    }

    /// Join: `Undetermined` absorbs (top). Two `Domain` values join by
    /// projecting both onto their shared vars and unioning the projected
    /// tuple sets. If the shared var set is empty, the proper LUB is
    /// `Undetermined`.
    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Domain(s), Self::Domain(o)) => relational_join(s, o),
        }
    }

    /// `Undetermined` matches anything (wildcard). A `Domain` pattern
    /// matches a `Domain` target iff `pattern.vars ⊆ target.vars` and
    /// every target tuple, projected to `pattern.vars`, is in
    /// `pattern.tuples`. Pattern over a var that target does not constrain
    /// is stricter than target, and so does not match.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Domain(p), Self::Domain(t)) => {
                let Some(projection) = positional_projection(&p.vars, &t.vars) else {
                    return false;
                };
                t.tuples.iter().all(|row| {
                    p.tuples
                        .contains(&projection.iter().map(|&i| row[i].clone()).collect())
                })
            }
        }
    }
}

/// Relational meet (natural join) of two `Domain` states. Returns the
/// canonical result (sorted vars, sorted + dedup'd tuples) or `None` if the
/// joined set is empty.
fn relational_meet(s: &DomainState, o: &DomainState) -> Option<JointDomainAst> {
    let (union_vars, ls, lo) = merge_vars(&s.vars, &o.vars);
    let mut shared: Vec<(usize, usize)> = Vec::new();
    for (sl, ol) in ls.iter().zip(lo.iter()) {
        if let (Some(s_i), Some(o_i)) = (sl, ol) {
            shared.push((*s_i, *o_i));
        }
    }
    let mut tuples: Vec<Vec<JointValue>> = Vec::new();
    for t1 in &s.tuples {
        for t2 in &o.tuples {
            if shared.iter().all(|&(i, j)| t1[i] == t2[j]) {
                let mut row: Vec<JointValue> = Vec::with_capacity(union_vars.len());
                for u_i in 0..union_vars.len() {
                    if let Some(o_i) = lo[u_i] {
                        row.push(t2[o_i].clone());
                    } else {
                        row.push(t1[ls[u_i].expect("missing self projection")].clone());
                    }
                }
                tuples.push(row);
            }
        }
    }
    tuples.sort();
    tuples.dedup();
    if tuples.is_empty() {
        None
    } else {
        Some(JointDomainAst::Domain(DomainState {
            vars: union_vars,
            tuples,
        }))
    }
}

/// Relational join over shared vars: project both tuple-lists onto the
/// intersection of their var sets, union, dedup. If the shared set is
/// empty, returns `Undetermined` (the proper LUB).
fn relational_join(s: &DomainState, o: &DomainState) -> JointDomainAst {
    let shared_vars: Vec<JointVar> = s.vars.iter().filter(|v| o.vars.contains(v)).copied().collect();
    if shared_vars.is_empty() {
        return JointDomainAst::Undetermined;
    }
    let project_s: Vec<usize> = shared_vars
        .iter()
        .map(|sv| s.vars.iter().position(|v| v == sv).expect("shared in self"))
        .collect();
    let project_o: Vec<usize> = shared_vars
        .iter()
        .map(|sv| o.vars.iter().position(|v| v == sv).expect("shared in other"))
        .collect();
    let mut tuples: Vec<Vec<JointValue>> = Vec::new();
    for t in &s.tuples {
        tuples.push(project_s.iter().map(|&i| t[i].clone()).collect());
    }
    for t in &o.tuples {
        tuples.push(project_o.iter().map(|&i| t[i].clone()).collect());
    }
    tuples.sort();
    tuples.dedup();
    JointDomainAst::Domain(DomainState {
        vars: shared_vars,
        tuples,
    })
}

/// Merge two sorted var lists into the union (sorted). Returns the union
/// plus, for each union position, the indices into `vs` and `vo`
/// (`None` if absent).
fn merge_vars(
    vs: &[JointVar],
    vo: &[JointVar],
) -> (Vec<JointVar>, Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut union: Vec<JointVar> = Vec::new();
    let mut from_s: Vec<Option<usize>> = Vec::new();
    let mut from_o: Vec<Option<usize>> = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < vs.len() && j < vo.len() {
        match vs[i].cmp(&vo[j]) {
            std::cmp::Ordering::Less => {
                union.push(vs[i]);
                from_s.push(Some(i));
                from_o.push(None);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                union.push(vo[j]);
                from_s.push(None);
                from_o.push(Some(j));
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                union.push(vs[i]);
                from_s.push(Some(i));
                from_o.push(Some(j));
                i += 1;
                j += 1;
            }
        }
    }
    while i < vs.len() {
        union.push(vs[i]);
        from_s.push(Some(i));
        from_o.push(None);
        i += 1;
    }
    while j < vo.len() {
        union.push(vo[j]);
        from_s.push(None);
        from_o.push(Some(j));
        j += 1;
    }
    (union, from_s, from_o)
}

/// For `matches`: positional projection from target vars to pattern vars.
/// Returns `Some(indices)` where `indices[i]` is the position in target
/// of pattern var `i`. Returns `None` if any pattern var is missing from
/// target (i.e., pattern is over a var target does not constrain).
fn positional_projection(pattern: &[JointVar], target: &[JointVar]) -> Option<Vec<usize>> {
    pattern
        .iter()
        .map(|pv| target.iter().position(|tv| tv == pv))
        .collect()
}

/// Compute the permutation that sorts `vars` into ascending order, rejecting
/// duplicates encountered during the sort.
fn sort_permutation(vars: &[JointVar]) -> Result<Vec<usize>, JointDomainError> {
    let mut indices: Vec<usize> = (0..vars.len()).collect();
    indices.sort_by_key(|&i| vars[i]);
    for window in indices.windows(2) {
        if vars[window[0]] == vars[window[1]] {
            return Err(JointDomainError::DuplicateVar(vars[window[0]]));
        }
    }
    Ok(indices)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    fn test_joint_domain_ast_default_is_undetermined() {
        assert_eq!(JointDomainAst::default(), JointDomainAst::Undetermined);
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints() {
        let jd = JointDomainAst::from_ints(
            vec![JointVar::LonePairs, JointVar::UnpairedElectrons],
            vec![vec![3, 0], vec![1, 4]],
        )
        .unwrap();
        assert_eq!(
            jd.vars(),
            Some(&[JointVar::LonePairs, JointVar::UnpairedElectrons][..])
        );
        assert_eq!(
            jd.tuples(),
            Some(
                &[
                    vec![JointValue::Int(1), JointValue::Int(4)],
                    vec![JointValue::Int(3), JointValue::Int(0)],
                ][..]
            )
        );
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_single_var() {
        let jd =
            JointDomainAst::from_ints(vec![JointVar::Charge], vec![vec![0], vec![1], vec![2]])
                .unwrap();
        assert_eq!(jd.vars(), Some(&[JointVar::Charge][..]));
        assert_eq!(
            jd.tuples(),
            Some(
                &[
                    vec![JointValue::Int(0)],
                    vec![JointValue::Int(1)],
                    vec![JointValue::Int(2)],
                ][..]
            )
        );
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_single_tuple_is_ground() {
        let jd = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1]],
        )
        .unwrap();
        assert!(jd.is_ground());
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_canonicalizes_var_order() {
        let user_order = JointDomainAst::from_ints(
            vec![JointVar::UnpairedElectrons, JointVar::LonePairs],
            vec![vec![0, 3], vec![4, 1]],
        )
        .unwrap();
        let canonical = JointDomainAst::from_ints(
            vec![JointVar::LonePairs, JointVar::UnpairedElectrons],
            vec![vec![3, 0], vec![1, 4]],
        )
        .unwrap();
        assert_eq!(user_order, canonical);
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_dedups_tuples() {
        let jd = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2], vec![0, 1]],
        )
        .unwrap();
        assert_eq!(
            jd.tuples(),
            Some(
                &[
                    vec![JointValue::Int(0), JointValue::Int(1)],
                    vec![JointValue::Int(1), JointValue::Int(2)],
                ][..]
            )
        );
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_error_zero_vars() {
        let err =
            JointDomainAst::from_ints(Vec::<JointVar>::new(), vec![vec![0], vec![1]]).unwrap_err();
        assert_eq!(err, JointDomainError::TooFewVars(0));
    }

    #[rstest]
    #[case::zero(Vec::<Vec<i64>>::new(), JointDomainError::TooFewTuples(0))]
    fn test_joint_domain_ast_from_ints_error_zero_tuples(
        #[case] tuples: Vec<Vec<i64>>,
        #[case] expected: JointDomainError,
    ) {
        let err = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            tuples,
        )
        .unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_error_arity_mismatch() {
        let err = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2, 3]],
        )
        .unwrap_err();
        assert_eq!(
            err,
            JointDomainError::ArityMismatch {
                tuple_index: 1,
                tuple_len: 3,
                vars_len: 2,
            }
        );
    }

    #[rstest]
    fn test_joint_domain_ast_from_ints_error_duplicate_var() {
        let err = JointDomainAst::from_ints(
            vec![JointVar::Charge, JointVar::Charge],
            vec![vec![0, 0], vec![1, 1]],
        )
        .unwrap_err();
        assert_eq!(err, JointDomainError::DuplicateVar(JointVar::Charge));
    }

    fn jd(vars: Vec<JointVar>, tuples: Vec<Vec<i64>>) -> JointDomainAst {
        JointDomainAst::from_ints(vars, tuples).unwrap()
    }

    #[rstest]
    fn test_joint_domain_ast_is_undetermined() {
        assert!(JointDomainAst::Undetermined.is_undetermined());
        assert!(!jd(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2]],
        )
        .is_undetermined());
    }

    #[rstest]
    fn test_joint_domain_ast_is_ground() {
        assert!(!JointDomainAst::Undetermined.is_ground());
        assert!(jd(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1]],
        )
        .is_ground());
        assert!(!jd(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 1], vec![1, 2]],
        )
        .is_ground());
    }

    #[rstest]
    #[case::und_und(JointDomainAst::Undetermined, JointDomainAst::Undetermined, Some(JointDomainAst::Undetermined))]
    #[case::und_dom(
        JointDomainAst::Undetermined,
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        Some(jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]])),
    )]
    #[case::same_vars_intersect(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![1, 2], vec![2, 3]]),
        Some(jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![1, 2]])),
    )]
    #[case::same_vars_disjoint_is_none(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![1, 2]]),
        None,
    )]
    #[case::disjoint_vars_cartesian(
        jd(vec![JointVar::Charge], vec![vec![0], vec![1]]),
        jd(vec![JointVar::ImplicitHydrogens], vec![vec![2], vec![3]]),
        Some(jd(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens],
            vec![vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3]],
        )),
    )]
    #[case::overlap_equijoin(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
        jd(vec![JointVar::ImplicitHydrogens, JointVar::LonePairs], vec![vec![1, 3], vec![2, 4]]),
        Some(jd(
            vec![JointVar::Charge, JointVar::ImplicitHydrogens, JointVar::LonePairs],
            vec![vec![0, 1, 3], vec![1, 2, 4]],
        )),
    )]
    #[case::overlap_no_match(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        jd(vec![JointVar::ImplicitHydrogens, JointVar::LonePairs], vec![vec![2, 3]]),
        None,
    )]
    fn test_joint_domain_ast_meet(
        #[case] a: JointDomainAst,
        #[case] b: JointDomainAst,
        #[case] expected: Option<JointDomainAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::und_und(JointDomainAst::Undetermined, JointDomainAst::Undetermined, JointDomainAst::Undetermined)]
    #[case::und_dom(
        JointDomainAst::Undetermined,
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        JointDomainAst::Undetermined,
    )]
    #[case::disjoint_vars_top(
        jd(vec![JointVar::Charge], vec![vec![0]]),
        jd(vec![JointVar::ImplicitHydrogens], vec![vec![1]]),
        JointDomainAst::Undetermined,
    )]
    #[case::same_vars_union(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![1, 2]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
    )]
    #[case::overlap_project_to_shared(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
        jd(vec![JointVar::ImplicitHydrogens, JointVar::LonePairs], vec![vec![3, 4]]),
        jd(vec![JointVar::ImplicitHydrogens], vec![vec![1], vec![2], vec![3]]),
    )]
    fn test_joint_domain_ast_join(
        #[case] a: JointDomainAst,
        #[case] b: JointDomainAst,
        #[case] expected: JointDomainAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::und_und(JointDomainAst::Undetermined, JointDomainAst::Undetermined, true)]
    #[case::und_dom(
        JointDomainAst::Undetermined,
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        true,
    )]
    #[case::dom_und(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        JointDomainAst::Undetermined,
        false,
    )]
    #[case::pattern_superset_of_target(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        true,
    )]
    #[case::pattern_not_subset(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1], vec![1, 2]]),
        false,
    )]
    #[case::pattern_var_missing_from_target(
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 1]]),
        jd(vec![JointVar::Charge, JointVar::LonePairs], vec![vec![0, 1]]),
        false,
    )]
    #[case::pattern_fewer_vars_target_projects(
        jd(vec![JointVar::Charge], vec![vec![0], vec![1]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 5], vec![1, 6]]),
        true,
    )]
    #[case::pattern_fewer_vars_target_violates(
        jd(vec![JointVar::Charge], vec![vec![0]]),
        jd(vec![JointVar::Charge, JointVar::ImplicitHydrogens], vec![vec![0, 5], vec![1, 6]]),
        false,
    )]
    fn test_joint_domain_ast_matches(
        #[case] pattern: JointDomainAst,
        #[case] target: JointDomainAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
