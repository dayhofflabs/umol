//! Isotope completion and projection, independent of valence resolution.

use thiserror::Error;
use umol_graph_ir::ir::{AtomHandle, AtomId, Edits, IsotopeMassForm, Lattice, Molecule};
use umol_utils::solution::Solution;

use super::WildcardPolicy;

/// How isotope resolution handles an Undetermined isotope.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IsotopePolicy {
    /// Retain Undetermined; do not choose an isotope composition.
    #[default]
    Strict,
    /// Complete Undetermined with naturally occurring isotope composition.
    Natural,
}

/// Resolves and projects isotope fields using the same isotope and wildcard policies.
///
/// WildcardPolicy::Strict applies the isotope policy to every atom. WildcardPolicy::Admit
/// skips atoms whose element is undetermined during completion and default elision.
/// Projection requires ground isotopes on every atom. Only isotope fields can change.
/// Supplied masses, sets, and variables are preserved.
///
/// # Semantic properties
///
/// Resolution is idempotent and preserves every supplied non-Undetermined isotope.
/// It returns Determined when every participating atom has a ground isotope after
/// completion. resolve_into preserves the caller's molecule
/// on underdetermination or error. resolve consumes it on every outcome.
///
/// For every molecule with ground isotope fields, successful projection followed
/// by resolution under the same policies recovers the original molecule exactly.
/// This law does not require the molecule's other fields to be ground.
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input and policies.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IsotopeResolver {
    policy: IsotopePolicy,
    wildcards: WildcardPolicy,
}

/// Isotope defaulting has no semantic contradiction outcomes.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeContradiction {}

/// Isotope completion has no operational failure outcomes.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeError {}

/// Failures projecting a resolved isotope description.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeProjectError {
    #[error("atom {atom:?} has a non-ground isotope")]
    NonGroundIsotope { atom: AtomId },
}

impl IsotopeResolver {
    /// Constructs an isotope resolver with strict isotope and wildcard policies.
    pub fn new() -> Self {
        Self::default()
    }

    /// Constructs an isotope resolver with the selected isotope and wildcard policies.
    pub fn with_config(policy: IsotopePolicy, wildcards: WildcardPolicy) -> Self {
        Self { policy, wildcards }
    }

    /// Plans isotope defaults without mutating the molecule.
    ///
    /// WildcardPolicy::Admit skips wildcard atoms. On participating atoms, IsotopePolicy::Natural
    /// completes an Undetermined isotope with natural composition. IsotopePolicy::Strict makes no edits.
    /// A remaining non-ground isotope on a participating atom yields Underdetermined
    /// with any proposed defaults retained in the plan.
    pub fn plan_resolve(&self, molecule: &Molecule) -> Solution<Edits, IsotopeContradiction> {
        let mut isotope_masses = Vec::new();
        let mut determined = true;
        for atom in molecule.atoms().iter() {
            let attributes = atom.attributes();
            if self.wildcards == WildcardPolicy::Admit && attributes.is_wildcard() {
                continue;
            }
            if self.policy == IsotopePolicy::Natural && attributes.isotope_mass.is_undetermined() {
                isotope_masses.push((
                    AtomHandle::Id(atom.id()),
                    IsotopeMassForm::Undetermined,
                    IsotopeMassForm::Natural,
                ));
            } else if !attributes.isotope_mass.is_ground() {
                determined = false;
            }
        }
        let mut edits = Edits::new();
        if !isotope_masses.is_empty() {
            edits.modify_atom_isotope_masses(isotope_masses);
        }
        if determined {
            Solution::Determined(edits)
        } else {
            Solution::Underdetermined(edits)
        }
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [IsotopeResolver::resolve_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, IsotopeContradiction, ()>, IsotopeError> {
        let edits = match self.plan_resolve(&molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("isotope resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("isotope resolve plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Plans and applies a determined isotope completion without a recovery journal.
    ///
    /// Underdetermination leaves the entire molecule unchanged, including fields
    /// for which the plan proposed defaults.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), IsotopeContradiction>, IsotopeError> {
        let edits = match self.plan_resolve(molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(Solution::Determined(()))
    }

    /// Projects the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [IsotopeResolver::project_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    ///
    /// # Errors
    ///
    /// Returns the planning error described by [IsotopeResolver::plan_project].
    pub fn project(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, IsotopeContradiction, ()>, IsotopeProjectError> {
        let edits = self.plan_project(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("isotope project plan applies to its input");
        let molecule = editor
            .finish()
            .expect("isotope project plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Elides Natural isotope composition when this policy restores it on resolution.
    ///
    /// IsotopePolicy::Natural elides natural composition on participating atoms.
    /// Admitted wildcard atoms retain it. IsotopePolicy::Strict retains all isotope fields.
    /// Explicit masses are preserved and other fields may remain unresolved.
    /// Planning finishes before mutation and the plan applies without a recovery journal.
    ///
    /// # Errors
    ///
    /// Returns IsotopeProjectError::NonGroundIsotope for the first atom whose
    /// isotope is Undetermined, a set, or a variable. Every error leaves the input unchanged.
    pub fn project_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), IsotopeContradiction>, IsotopeProjectError> {
        let edits = self.plan_project(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(Solution::Determined(()))
    }

    /// Plans to elide natural isotope composition wherever resolution would restore it.
    ///
    /// Applying the returned edits to the same molecule produces the result of
    /// [IsotopeResolver::project]. The molecule remains unchanged while planning.
    ///
    /// # Errors
    ///
    /// Returns IsotopeProjectError::NonGroundIsotope for the first non-ground isotope,
    /// including one on an admitted wildcard atom.
    pub fn plan_project(&self, molecule: &Molecule) -> Result<Edits, IsotopeProjectError> {
        let mut isotope_masses = Vec::new();
        for atom in molecule.atoms().iter() {
            let attributes = atom.attributes();
            if !attributes.isotope_mass.is_ground() {
                return Err(IsotopeProjectError::NonGroundIsotope { atom: atom.id() });
            }
            if self.policy == IsotopePolicy::Natural
                && matches!(attributes.isotope_mass, IsotopeMassForm::Natural)
                && (self.wildcards == WildcardPolicy::Strict || !attributes.is_wildcard())
            {
                isotope_masses.push((
                    AtomHandle::Id(atom.id()),
                    IsotopeMassForm::Natural,
                    IsotopeMassForm::Undetermined,
                ));
            }
        }
        let mut edits = Edits::new();
        if !isotope_masses.is_empty() {
            edits.modify_atom_isotope_masses(isotope_masses);
        }
        Ok(edits)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AtomForm, Edit, MoleculeApplyError, MoleculeEntries, RollbackMethod, TransactionError,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    fn test_isotope_resolver_new() {
        assert_eq!(
            IsotopeResolver::new(),
            IsotopeResolver {
                policy: IsotopePolicy::Strict,
                wildcards: WildcardPolicy::Strict
            },
        );
    }

    #[rstest]
    fn test_isotope_resolver_with_config(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
    ) {
        assert_eq!(
            IsotopeResolver::with_config(policy, wildcards),
            IsotopeResolver { policy, wildcards },
        );
    }

    #[rstest]
    #[case::natural(
        IsotopePolicy::Natural,
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#h4" "C#i13" "N#i="]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomIsotopeMasses {
            changes: vec![(
                AtomHandle::Id(AtomId(0)),
                IsotopeMassForm::Undetermined,
                IsotopeMassForm::Natural,
            )],
        }])),
    )]
    #[case::strict(
        IsotopePolicy::Strict,
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C" "C#i13"]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::partial(
        IsotopePolicy::Natural,
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C" "C#i?mass" "O"]}"#),
        Solution::Underdetermined(Edits::from_iter([Edit::ModifyAtomIsotopeMasses {
            changes: vec![
                (AtomHandle::Id(AtomId(0)), IsotopeMassForm::Undetermined, IsotopeMassForm::Natural),
                (AtomHandle::Id(AtomId(2)), IsotopeMassForm::Undetermined, IsotopeMassForm::Natural),
            ],
        }])),
    )]
    #[case::strict_mixed(
        IsotopePolicy::Natural, WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["*" "C"]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomIsotopeMasses {
            changes: vec![
                (AtomHandle::Id(AtomId(0)), IsotopeMassForm::Undetermined, IsotopeMassForm::Natural),
                (AtomHandle::Id(AtomId(1)), IsotopeMassForm::Undetermined, IsotopeMassForm::Natural),
            ],
        }])),
    )]
    #[case::admit_mixed(
        IsotopePolicy::Natural, WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*" "C"]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomIsotopeMasses {
            changes: vec![(AtomHandle::Id(AtomId(1)), IsotopeMassForm::Undetermined, IsotopeMassForm::Natural)],
        }])),
    )]
    #[case::admit_wildcard(
        IsotopePolicy::Natural, WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*"]}"#), Solution::Determined(Edits::new()),
    )]
    #[case::admit_wildcard_strict_isotope(
        IsotopePolicy::Strict, WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*"]}"#), Solution::Determined(Edits::new()),
    )]
    #[case::admit_mass(
        IsotopePolicy::Natural, WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*#i1"]}"#), Solution::Determined(Edits::new()),
    )]
    #[case::admit_strict_isotope(
        IsotopePolicy::Strict, WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*" "C"]}"#), Solution::Underdetermined(Edits::new()),
    )]
    fn test_isotope_resolver_plan_resolve(
        #[case] policy: IsotopePolicy,
        #[case] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
        #[case] expected: Solution<Edits, IsotopeContradiction>,
    ) {
        assert_eq!(
            IsotopeResolver::with_config(policy, wildcards).plan_resolve(&molecule),
            expected
        );
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::ground(mol_dsl!(r#"{:atoms ["C#i=" "C#i13"]}"#))]
    fn test_isotope_resolver_plan_resolve_identity(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            IsotopeResolver::with_config(policy, wildcards).plan_resolve(&molecule),
            Solution::Determined(Edits::new()),
        );
    }

    #[rstest]
    #[case::set(IsotopeMassForm::lit_set([12, 13]))]
    #[case::singleton(IsotopeMassForm::lit_set([13]))]
    #[case::empty_set(IsotopeMassForm::lit_set([]))]
    #[case::variable(IsotopeMassForm::var("mass"))]
    #[case::restricted_variable(IsotopeMassForm::var_in("mass", [12, 13]))]
    fn test_isotope_resolver_plan_resolve_partial(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[case] isotope: IsotopeMassForm,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                isotope_mass: isotope,
                ..Default::default()
            }],
            ..Default::default()
        });
        assert_eq!(
            IsotopeResolver::with_config(policy, WildcardPolicy::Strict).plan_resolve(&molecule),
            Solution::Underdetermined(Edits::new()),
        );
    }

    #[rstest]
    fn test_isotope_resolver_plan_resolve_stale() {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let Solution::Determined(edits) =
            IsotopeResolver::with_config(IsotopePolicy::Natural, WildcardPolicy::Strict)
                .plan_resolve(&molecule)
        else {
            panic!("fixture must produce a determined edit plan");
        };
        molecule.atom_mut(AtomId(1)).attributes_mut().isotope_mass = IsotopeMassForm::Lit(18);
        let expected = molecule.clone();
        assert_eq!(
            molecule.clone().edit().apply(edits.clone()),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(
            molecule.transact([edits], RollbackMethod::Journal),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::partial_fields(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#i=#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
    )]
    #[case::dative(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["N#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N#i=#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    )]
    #[case::admit(
        WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*" "*#i1" "*#i{1,2}" "*#i?mass" "C"]}"#),
        mol_dsl!(r#"{:atoms ["*" "*#i1" "*#i{1,2}" "*#i?mass" "C#i="]}"#),
    )]
    #[case::symbolic_elements(
        WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["{C,N}" "?e"]}"#),
        mol_dsl!(r#"{:atoms ["{C,N}#i=" "?e#i="]}"#),
    )]
    fn test_isotope_resolver_resolve(
        #[case] wildcards: WildcardPolicy,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = IsotopeResolver::with_config(IsotopePolicy::Natural, wildcards);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Determined(())),
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new(), Solution::Determined(()))]
    #[case::ground(mol_dsl!(r#"{:atoms ["C#i=" "C#i13"]}"#), Solution::Determined(()))]
    #[case::variable(mol_dsl!(r#"{:atoms ["C#i?mass"]}"#), Solution::Underdetermined(()))]
    #[case::set(mol_dsl!(r#"{:atoms ["C#i{12,13}"]}"#), Solution::Underdetermined(()))]
    fn test_isotope_resolver_resolve_into_identity(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<(), IsotopeContradiction>,
    ) {
        let original = molecule.clone();
        assert_eq!(
            IsotopeResolver::with_config(policy, wildcards).resolve_into(&mut molecule),
            Ok(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::default(IsotopeResolver::default(), mol_dsl!(r#"{:atoms ["C"]}"#))]
    #[case::strict(IsotopeResolver::new(), mol_dsl!(r#"{:atoms ["C"]}"#))]
    #[case::natural_late_variable(IsotopeResolver::with_config(IsotopePolicy::Natural, WildcardPolicy::Strict), mol_dsl!(r#"{:atoms ["C" "O#i?mass"]}"#))]
    #[case::natural_early_variable(IsotopeResolver::with_config(IsotopePolicy::Natural, WildcardPolicy::Strict), mol_dsl!(r#"{:atoms ["C#i?mass" "O"]}"#))]
    fn test_isotope_resolver_resolve_partial(
        #[case] resolver: IsotopeResolver,
        #[case] mut molecule: Molecule,
    ) {
        let original = molecule.clone();
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::partial_fields(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["C#i=#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
    )]
    #[case::dative(
        WildcardPolicy::Strict,
        mol_dsl!(r#"{:atoms ["N#i=#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    )]
    #[case::admit(
        WildcardPolicy::Admit,
        mol_dsl!(r#"{:atoms ["*#i=" "*#i1" "C#i=" "C#i13"]}"#),
        mol_dsl!(r#"{:atoms ["*#i=" "*#i1" "C" "C#i13"]}"#),
    )]
    fn test_isotope_resolver_project(
        #[case] wildcards: WildcardPolicy,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = IsotopeResolver::with_config(IsotopePolicy::Natural, wildcards);
        assert_eq!(
            resolver.project(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            resolver.project_into(&mut molecule),
            Ok(Solution::Determined(())),
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::strict_natural(IsotopeResolver::new(), mol_dsl!(r#"{:atoms ["C#i=" "O#i18"]}"#))]
    #[case::default_natural(IsotopeResolver::default(), mol_dsl!(r#"{:atoms ["C#i=" "O#i18"]}"#))]
    #[case::strict_mass(IsotopeResolver::new(), mol_dsl!(r#"{:atoms ["C#i13"]}"#))]
    #[case::natural_mass(IsotopeResolver::with_config(IsotopePolicy::Natural, WildcardPolicy::Strict), mol_dsl!(r#"{:atoms ["C#i13"]}"#))]
    #[case::strict_empty(IsotopeResolver::new(), Molecule::new())]
    #[case::natural_empty(
        IsotopeResolver::with_config(IsotopePolicy::Natural, WildcardPolicy::Strict),
        Molecule::new()
    )]
    fn test_isotope_resolver_project_identity(
        #[case] resolver: IsotopeResolver,
        #[case] mut molecule: Molecule,
    ) {
        let original = molecule.clone();
        assert_eq!(
            resolver.project(molecule.clone()),
            Ok(Solution::Determined(original.clone()))
        );
        assert_eq!(
            resolver.project_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::undetermined(IsotopeMassForm::Undetermined)]
    #[case::set(IsotopeMassForm::lit_set([12, 13]))]
    #[case::singleton(IsotopeMassForm::lit_set([13]))]
    #[case::empty_set(IsotopeMassForm::lit_set([]))]
    #[case::variable(IsotopeMassForm::var("mass"))]
    #[case::restricted_variable(IsotopeMassForm::var_in("mass", [12, 13]))]
    fn test_isotope_resolver_project_error(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
        #[case] isotope: IsotopeMassForm,
    ) {
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm {
                    isotope_mass: IsotopeMassForm::Natural,
                    ..Default::default()
                },
                AtomForm {
                    isotope_mass: isotope,
                    ..Default::default()
                },
            ],
            ..Default::default()
        });
        let resolver = IsotopeResolver::with_config(policy, wildcards);
        let original = molecule.clone();
        assert_eq!(
            resolver.project(molecule.clone()),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
        assert_eq!(
            resolver.project_into(&mut molecule),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::natural(IsotopePolicy::Natural, Edits::from_iter([
        Edit::ModifyAtomIsotopeMasses {
            changes: vec![
                (AtomHandle::Id(AtomId(0)), IsotopeMassForm::Natural, IsotopeMassForm::Undetermined),
                (AtomHandle::Id(AtomId(2)), IsotopeMassForm::Natural, IsotopeMassForm::Undetermined),
            ],
        },
    ]))]
    #[case::strict(IsotopePolicy::Strict, Edits::new())]
    fn test_isotope_resolver_plan_project(#[case] policy: IsotopePolicy, #[case] expected: Edits) {
        let molecule = mol_dsl!(r#"{:atoms ["C#i=" "O#i18" "N#i="]}"#);
        assert_eq!(
            IsotopeResolver::with_config(policy, WildcardPolicy::Strict).plan_project(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C#i=" "O"]}"#))]
    #[case::set(mol_dsl!(r#"{:atoms ["C#i=" "O#i{16,18}"]}"#))]
    #[case::wildcard(mol_dsl!(r#"{:atoms ["C#i=" "*"]}"#))]
    fn test_isotope_resolver_plan_project_error(
        #[values(WildcardPolicy::Strict, WildcardPolicy::Admit)] wildcards: WildcardPolicy,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            IsotopeResolver::with_config(IsotopePolicy::Natural, wildcards).plan_project(&molecule),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
    }
}
