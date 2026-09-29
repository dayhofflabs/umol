//! Isotope completion and projection under one policy, independent of valence resolution.

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomHandle, AtomId, AtomUpdate, Edits, IsotopeMassForm, Lattice, Molecule, MoleculeApplyError,
};
use umol_utils::solution::Solution;

/// How isotope resolution handles an Undetermined isotope.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IsotopePolicy {
    /// Retain Undetermined; do not choose an isotope composition.
    #[default]
    Strict,
    /// Complete Undetermined with naturally occurring isotope composition.
    Natural,
}

/// Resolves and projects isotope fields using the same policy.
///
/// Only isotope fields participate. Other atom fields, topology, entities, and
/// constraints are preserved, whether determined or not. Supplied masses, sets,
/// and variables are neither normalized nor interpreted.
///
/// # Semantic properties
///
/// Resolution is idempotent and preserves every supplied non-Undetermined isotope.
/// It publishes changes only when every isotope is ground after completion.
/// resolve_into preserves the caller's molecule on underdetermination or error;
/// resolve consumes it on every outcome.
///
/// For every molecule with ground isotope fields, successful projection followed
/// by resolution under the same policy recovers the original molecule exactly.
/// This law does not require the molecule's other fields to be ground.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IsotopeResolver {
    policy: IsotopePolicy,
}

/// Isotope defaulting has no semantic contradiction outcomes.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeContradiction {}

/// Operational failures applying isotope completion.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeError {
    #[error(transparent)]
    Apply(#[from] MoleculeApplyError),
}

/// Failures projecting a resolved isotope description.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum IsotopeProjectError {
    #[error("atom {atom:?} has a non-ground isotope")]
    NonGroundIsotope { atom: AtomId },
    #[error(transparent)]
    Apply(#[from] MoleculeApplyError),
}

impl IsotopeResolver {
    /// Constructs an isotope resolver with the selected policy.
    pub fn new(policy: IsotopePolicy) -> Self {
        Self { policy }
    }

    /// Plans isotope defaults without mutating the molecule.
    ///
    /// Under Natural, only Undetermined becomes Natural. Strict makes no edits.
    /// Remaining non-ground isotopes yield Underdetermined with any proposed
    /// default edits retained in the plan.
    pub fn plan_resolve(&self, molecule: &Molecule) -> Solution<Edits, IsotopeContradiction> {
        let mut edits = Edits::new();
        let mut determined = true;
        for atom in molecule.atoms().iter() {
            if self.policy == IsotopePolicy::Natural
                && atom.attributes().isotope_mass.is_undetermined()
            {
                edits.update_atom(
                    AtomHandle::Id(atom.id()),
                    atom.attributes(),
                    &AtomUpdate {
                        isotope_mass: Some(IsotopeMassForm::Natural),
                        ..Default::default()
                    },
                );
            } else if !atom.attributes().isotope_mass.is_ground() {
                determined = false;
            }
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
    ///
    /// # Errors
    ///
    /// Returns IsotopeError::Apply if application or the final integrity check fails.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, IsotopeContradiction, ()>, IsotopeError> {
        let edits = match self.plan_resolve(&molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        let molecule = molecule
            .edit()
            .apply(edits)?
            .finish()
            .map_err(MoleculeApplyError::from)?;
        Ok(Solution::Determined(molecule))
    }

    /// Plans and atomically applies a determined isotope completion.
    ///
    /// Underdetermination leaves the entire molecule unchanged, including fields
    /// for which the plan proposed defaults.
    ///
    /// # Errors
    ///
    /// Returns IsotopeError::Apply if edit application or the final integrity check fails.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), IsotopeContradiction>, IsotopeError> {
        let edits = match self.plan_resolve(molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        molecule.transact([edits])?;
        Ok(Solution::Determined(()))
    }

    /// Projects the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [IsotopeResolver::project_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    ///
    /// # Errors
    ///
    /// Returns the planning or application error, including failure of the final integrity check.
    pub fn project(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, IsotopeContradiction, ()>, IsotopeProjectError> {
        let edits = self.plan_project(&molecule)?;
        let molecule = molecule
            .edit()
            .apply(edits)?
            .finish()
            .map_err(MoleculeApplyError::from)?;
        Ok(Solution::Determined(molecule))
    }

    /// Elides Natural isotope composition when this policy restores it on resolution.
    ///
    /// Natural becomes Undetermined only under Natural; Strict retains it. Both
    /// policies preserve explicit masses. Other fields may remain unresolved.
    /// Only a successful result publishes the projected molecule.
    ///
    /// # Errors
    ///
    /// Returns IsotopeProjectError::NonGroundIsotope for the first atom whose
    /// isotope is Undetermined, a set, or a variable, or Apply if edit application
    /// or the final integrity check fails. Every error preserves the caller's molecule
    /// under Molecule::normalized_eq.
    pub fn project_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), IsotopeContradiction>, IsotopeProjectError> {
        let edits = self.plan_project(molecule)?;
        molecule.transact([edits])?;
        Ok(Solution::Determined(()))
    }

    /// Plans isotope default elision under this policy without mutating the molecule.
    ///
    /// Applying the returned edits to the same molecule produces the result of
    /// [IsotopeResolver::project]. The molecule remains unchanged while planning.
    ///
    /// # Errors
    ///
    /// Returns IsotopeProjectError::NonGroundIsotope for the first non-ground isotope.
    pub fn plan_project(&self, molecule: &Molecule) -> Result<Edits, IsotopeProjectError> {
        let mut edits = Edits::new();
        for atom in molecule.atoms().iter() {
            if !atom.attributes().isotope_mass.is_ground() {
                return Err(IsotopeProjectError::NonGroundIsotope { atom: atom.id() });
            }
            if self.policy == IsotopePolicy::Natural
                && matches!(atom.attributes().isotope_mass, IsotopeMassForm::Natural)
            {
                edits.update_atom(
                    AtomHandle::Id(atom.id()),
                    atom.attributes(),
                    &AtomUpdate {
                        isotope_mass: Some(IsotopeMassForm::Undetermined),
                        ..Default::default()
                    },
                );
            }
        }
        Ok(edits)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AtomFieldChange, AtomForm, Edit, MoleculeEntries, TransactionError};
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::natural(
        IsotopePolicy::Natural,
        mol_dsl!(r#"{:atoms ["C#h4" "C#i13" "N#i="]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::IsotopeMass {
                old: IsotopeMassForm::Undetermined,
                new: IsotopeMassForm::Natural,
            },
        }])),
    )]
    #[case::strict(
        IsotopePolicy::Strict,
        mol_dsl!(r#"{:atoms ["C" "C#i13"]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::partial(
        IsotopePolicy::Natural,
        mol_dsl!(r#"{:atoms ["C" "C#i?mass" "O"]}"#),
        Solution::Underdetermined(Edits::from_iter([
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::IsotopeMass {
                    old: IsotopeMassForm::Undetermined,
                    new: IsotopeMassForm::Natural,
                },
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(2)),
                change: AtomFieldChange::IsotopeMass {
                    old: IsotopeMassForm::Undetermined,
                    new: IsotopeMassForm::Natural,
                },
            },
        ])),
    )]
    fn test_isotope_resolver_plan_resolve(
        #[case] policy: IsotopePolicy,
        #[case] molecule: Molecule,
        #[case] expected: Solution<Edits, IsotopeContradiction>,
    ) {
        assert_eq!(
            IsotopeResolver::new(policy).plan_resolve(&molecule),
            expected
        );
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::ground(mol_dsl!(r#"{:atoms ["C#i=" "C#i13"]}"#))]
    fn test_isotope_resolver_plan_resolve_identity(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] policy: IsotopePolicy,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            IsotopeResolver::new(policy).plan_resolve(&molecule),
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
            IsotopeResolver::new(policy).plan_resolve(&molecule),
            Solution::Underdetermined(Edits::new()),
        );
    }

    #[rstest]
    fn test_isotope_resolver_plan_resolve_stale() {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let Solution::Determined(edits) =
            IsotopeResolver::new(IsotopePolicy::Natural).plan_resolve(&molecule)
        else {
            panic!("fixture must produce a determined edit plan");
        };
        molecule.atom_mut(AtomId(1)).attributes_mut().isotope_mass = IsotopeMassForm::Lit(18);
        let expected = molecule.clone();
        assert_eq!(
            molecule.clone().edit().apply(edits.clone()).map(|_| ()),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(
            molecule.transact([edits]),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::partial_fields(
        mol_dsl!(r#"{:atoms ["C#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#i=#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
    )]
    #[case::dative(
        mol_dsl!(r#"{:atoms ["N#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N#i=#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    )]
    fn test_isotope_resolver_resolve(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            IsotopeResolver::new(IsotopePolicy::Natural).resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            IsotopeResolver::new(IsotopePolicy::Natural).resolve_into(&mut molecule),
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
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<(), IsotopeContradiction>,
    ) {
        let original = molecule.clone();
        assert_eq!(
            IsotopeResolver::new(policy).resolve_into(&mut molecule),
            Ok(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::default(IsotopeResolver::default(), mol_dsl!(r#"{:atoms ["C"]}"#))]
    #[case::strict(IsotopeResolver::new(IsotopePolicy::Strict), mol_dsl!(r#"{:atoms ["C"]}"#))]
    #[case::natural_late_variable(IsotopeResolver::new(IsotopePolicy::Natural), mol_dsl!(r#"{:atoms ["C" "O#i?mass"]}"#))]
    #[case::natural_early_variable(IsotopeResolver::new(IsotopePolicy::Natural), mol_dsl!(r#"{:atoms ["C#i?mass" "O"]}"#))]
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
        mol_dsl!(r#"{:atoms ["C#i=#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["C#c?charge#h{1,2}" "O#i18"] :bonds [[0 1 "1"]]}"#),
    )]
    #[case::dative(
        mol_dsl!(r#"{:atoms ["N#i=#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        mol_dsl!(r#"{:atoms ["N#c+#h3" "B#i11#c-"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    )]
    fn test_isotope_resolver_project(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            IsotopeResolver::new(IsotopePolicy::Natural).project(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            IsotopeResolver::new(IsotopePolicy::Natural).project_into(&mut molecule),
            Ok(Solution::Determined(())),
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::strict_natural(IsotopeResolver::new(IsotopePolicy::Strict), mol_dsl!(r#"{:atoms ["C#i=" "O#i18"]}"#))]
    #[case::default_natural(IsotopeResolver::default(), mol_dsl!(r#"{:atoms ["C#i=" "O#i18"]}"#))]
    #[case::strict_mass(IsotopeResolver::new(IsotopePolicy::Strict), mol_dsl!(r#"{:atoms ["C#i13"]}"#))]
    #[case::natural_mass(IsotopeResolver::new(IsotopePolicy::Natural), mol_dsl!(r#"{:atoms ["C#i13"]}"#))]
    #[case::strict_empty(IsotopeResolver::new(IsotopePolicy::Strict), Molecule::new())]
    #[case::natural_empty(IsotopeResolver::new(IsotopePolicy::Natural), Molecule::new())]
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
        let original = molecule.clone();
        assert_eq!(
            IsotopeResolver::new(policy).project(molecule.clone()),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
        assert_eq!(
            IsotopeResolver::new(policy).project_into(&mut molecule),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::natural(IsotopePolicy::Natural, Edits::from_iter([
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::IsotopeMass {
                old: IsotopeMassForm::Natural,
                new: IsotopeMassForm::Undetermined,
            },
        },
    ]))]
    #[case::strict(IsotopePolicy::Strict, Edits::new())]
    fn test_isotope_resolver_plan_project(#[case] policy: IsotopePolicy, #[case] expected: Edits) {
        let molecule = mol_dsl!(r#"{:atoms ["C#i=" "O#i18"]}"#);
        assert_eq!(
            IsotopeResolver::new(policy).plan_project(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C#i=" "O"]}"#))]
    #[case::set(mol_dsl!(r#"{:atoms ["C#i=" "O#i{16,18}"]}"#))]
    fn test_isotope_resolver_plan_project_error(#[case] molecule: Molecule) {
        assert_eq!(
            IsotopeResolver::new(IsotopePolicy::Natural).plan_project(&molecule),
            Err(IsotopeProjectError::NonGroundIsotope { atom: AtomId(1) }),
        );
    }
}
