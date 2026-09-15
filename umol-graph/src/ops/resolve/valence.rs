//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_graph_ir::ir::{AtomConstraintKey, AtomId, Molecule};
use umol_utils::solution::Solution;

use crate::ops::model::{ValenceCandidateSource, ValenceModel, ValenceTieBreak};
use crate::ops::resolve::ResolveState;
use crate::ops::valence::{AtomTypingError, AtomTypingValence, CountsError, CountsValence};
use crate::ops::validate::{
    DerivedKind, IncidenceConstraintInvariantsContradiction, IncidenceConstraintInvariantsValidator,
};

#[derive(Clone, Debug)]
pub enum ValenceResolver<'a> {
    AtomTyping(AtomTypingValence<'a>),
    Counts(CountsValence<'a>),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceContradiction {
    #[error(transparent)]
    Constraint(#[from] IncidenceConstraintInvariantsContradiction),
    #[error(transparent)]
    AtomTyping(#[from] AtomTypingError),
    #[error(transparent)]
    Counts(#[from] CountsError),
}

/// Operational failures of valence admission; currently uninhabited —
/// admission produces no edits and runs no transactions.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceError {}

/// Operational failures of valence projection; currently uninhabited.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceProjectError {}

impl<'a> ValenceResolver<'a> {
    pub fn new(model: &'a ValenceModel) -> Self {
        match &model.candidates {
            ValenceCandidateSource::AtomTyping { registry } => {
                Self::AtomTyping(AtomTypingValence::new(registry.as_ref()))
            }
            ValenceCandidateSource::Counts { table } => {
                Self::Counts(CountsValence::new(table.as_ref()))
            }
        }
    }

    /// Admission: the candidate sets of every atom under resolution, with
    /// the incidence-constraint invariants checked first; no edits are
    /// produced. The chemistry solution rides `Solution`; the operational
    /// channel is currently uninhabited.
    pub fn admit(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<ResolveState, ValenceContradiction>, ValenceError> {
        for atom in molecule.atoms().ids() {
            for key in [
                AtomConstraintKey::Valence,
                AtomConstraintKey::DonatedPairs,
                AtomConstraintKey::AcceptedPairs,
            ] {
                match IncidenceConstraintInvariantsValidator
                    .validate_molecule_atom_constraint(
                        molecule,
                        atom,
                        key,
                        DerivedKind::DerivedComplete,
                    )
                    .expect("atom id came from the molecule atom store")
                {
                    Solution::Determined(()) => {}
                    Solution::Underdetermined(()) => {
                        return Ok(Solution::Underdetermined(ResolveState::default()));
                    }
                    Solution::Contradictory(contradiction) => {
                        return Ok(Solution::Contradictory(contradiction.into()));
                    }
                }
            }
        }
        let completions = match self {
            Self::AtomTyping(resolver) => resolver
                .admit(molecule)
                .map_contradiction(ValenceContradiction::from),
            Self::Counts(resolver) => resolver
                .admit(molecule)
                .map_contradiction(ValenceContradiction::from),
        };
        Ok(completions.map(|completions| ResolveState {
            completions,
            systems: Vec::new(),
            tie_breaks: Vec::new(),
        }))
    }

    /// Preserve atom fields during valence projection.
    ///
    /// # Semantic properties
    ///
    /// Always returns Determined without reading or modifying the molecule, independently
    /// of the valence source and policy. Implicit-H counts and electron fields are preserved.
    pub fn project(
        &self,
        _molecule: &mut Molecule,
        _policy: ValenceTieBreak,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceProjectError> {
        Ok(Solution::Determined(()))
    }

    /// Infer implicit H using the selected valence source and policy, without changing the molecule.
    /// Returns None when the source cannot determine a count under that policy.
    pub(crate) fn infer_implicit_hydrogens(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        policy: ValenceTieBreak,
    ) -> Option<i64> {
        match self {
            Self::AtomTyping(resolver) => {
                resolver.infer_implicit_hydrogens(molecule, atom_id, policy)
            }
            Self::Counts(resolver) => resolver.infer_implicit_hydrogens(molecule, atom_id, policy),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use smallvec::smallvec;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomConstraintForm, AtomForm, AtomId};
    use umol_graph_ir::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::valence::{AtomCompletions, AtomTypeRegistry, ValenceTable};

    #[rstest]
    fn test_valence_resolver_new() {
        let counts = ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()));
        let atom_typing =
            ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                "C#c0#h4"
            )])));

        assert!(matches!(
            ValenceResolver::new(&counts),
            ValenceResolver::Counts(_)
        ));
        assert!(matches!(
            ValenceResolver::new(&atom_typing),
            ValenceResolver::AtomTyping(_)
        ));
    }

    #[rstest]
    #[case::counts(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        atom_dsl!("C#c0#h4#n0#u0#s#v0#a!"),
    )]
    #[case::atom_typing(
        ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#c0#h4#n0#u0#s#v0#a!"
        )]))),
        atom_dsl!("C#c0#h4#n0#u0#s#v0#a!"),
    )]
    fn test_valence_resolver_admit(#[case] model: ValenceModel, #[case] expected: AtomForm) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        let mut completions = AtomCompletions::new();
        completions.insert(AtomId(0), smallvec![expected]);
        assert_eq!(
            ValenceResolver::new(&model).admit(&molecule),
            Ok(Solution::Determined(ResolveState {
                completions,
                ..ResolveState::default()
            }))
        );
    }

    #[rstest]
    #[case::counts_contradictory(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["C#v1"]}"#),
        Solution::Contradictory(ValenceContradiction::Constraint(
            IncidenceConstraintInvariantsContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintForm::valence(1),
            },
        )),
    )]
    #[case::atom_typing_contradictory(
        ValenceModel::atom_typing(Cow::Borrowed(AtomTypeRegistry::default_registry())),
        mol_dsl!(r#"{:atoms ["C#v1"]}"#),
        Solution::Contradictory(ValenceContradiction::Constraint(
            IncidenceConstraintInvariantsContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintForm::valence(1),
            },
        )),
    )]
    #[case::counts_underdetermined(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "*"]]}"#),
        Solution::Underdetermined(ResolveState::default()),
    )]
    #[case::atom_typing_underdetermined(
        ValenceModel::atom_typing(Cow::Borrowed(AtomTypeRegistry::default_registry())),
        mol_dsl!(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "*"]]}"#),
        Solution::Underdetermined(ResolveState::default()),
    )]
    #[case::dative_pairs_contradictory(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["N#d0" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        Solution::Contradictory(ValenceContradiction::Constraint(
            IncidenceConstraintInvariantsContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintForm::donated_pairs(0),
            },
        )),
    )]
    #[case::accepted_pairs_contradictory(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["N" "B#t0"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
        Solution::Contradictory(ValenceContradiction::Constraint(
            IncidenceConstraintInvariantsContradiction::Atom {
                atom: AtomId(1),
                constraint: AtomConstraintForm::accepted_pairs(0),
            },
        )),
    )]
    #[case::multidonor_underdetermined(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["N#d1" "N" "B"] :dative-bonds [{:donors [0 1] :acceptor 2 :attrs "1"}]}"#),
        Solution::Underdetermined(ResolveState::default()),
    )]
    #[case::vacuous(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v*#a!"]}"#),
        Solution::Determined(ResolveState {
            completions: {
                let mut completions = AtomCompletions::new();
                completions.insert(
                    AtomId(0),
                    smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#a!")],
                );
                completions
            },
            systems: Vec::new(),
            tie_breaks: Vec::new(),
        }),
    )]
    fn test_valence_resolver_admit_constraints(
        #[case] model: ValenceModel,
        #[case] molecule: Molecule,
        #[case] expected: Solution<ResolveState, ValenceContradiction>,
    ) {
        assert_eq!(ValenceResolver::new(&model).admit(&molecule), Ok(expected));
    }

    #[rstest]
    #[case::counts(ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())))]
    #[case::atom_typing(ValenceModel::atom_typing(Cow::Borrowed(
        AtomTypeRegistry::default_registry()
    )))]
    fn test_valence_resolver_admit_partial(#[case] model: ValenceModel) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#);
        assert_eq!(
            ValenceResolver::new(&model).admit(&molecule),
            Ok(Solution::Underdetermined(ResolveState::default()))
        );
    }

    #[rstest]
    #[case::counts(
        ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
        mol_dsl!(r#"{:atoms ["C#c0#h4" "Fe#c0#h0#a+"]}"#),
        ValenceContradiction::Counts(CountsError::UndeterminedAromaticValence)
    )]
    #[case::atom_typing(
        ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0#h4")]))),
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0#h3"]}"#),
        ValenceContradiction::AtomTyping(AtomTypingError::NoMatch {
            atom_id: AtomId(1),
            element: Element::C,
            charge: Some(0),
        })
    )]
    fn test_valence_resolver_admit_error(
        #[case] model: ValenceModel,
        #[case] molecule: Molecule,
        #[case] expected: ValenceContradiction,
    ) {
        assert_eq!(
            ValenceResolver::new(&model).admit(&molecule),
            Ok(Solution::Contradictory(expected))
        );
    }

    #[rstest]
    #[case::methane(mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#))]
    #[case::fluorine(mol_dsl!(r#"{:atoms ["F#i=#c0#h1#n3#u0#s"]}"#))]
    #[case::zero_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h0#n0#u0#s#v4"]}"#))]
    #[case::aromatic_zero_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h0#n0#u0#s#v3#a1"]}"#))]
    #[case::lone_pairs(mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n2#u0#s"]}"#))]
    #[case::ethane(mol_dsl!(r#"{:atoms ["C#i=#c0#h3#n0#u0#s" "C#i=#c0#h3#n0#u0#s"] :bonds [[0 1 "1"]]}"#))]
    #[case::mixed(mol_dsl!(r#"{:atoms ["C#i=#c0#h4#u0#s" "C#i13#c0#h4#u0#s" "C#i=#c-#h3#u0#s" "C#i=#c0#h3#u1#s2"]}"#))]
    #[case::charge(mol_dsl!(r#"{:atoms ["C#i=#c+#h3#u0#s"]}"#))]
    #[case::radical(mol_dsl!(r#"{:atoms ["C#i=#c0#h3#u1#s2"]}"#))]
    #[case::open_shell_singlet(mol_dsl!(r#"{:atoms ["C#i=#c0#h2#u2#s"]}"#))]
    #[case::mass(mol_dsl!(r#"{:atoms ["C#i13#c0#h4#u0#s"]}"#))]
    #[case::mass_zero_h(mol_dsl!(r#"{:atoms ["C#i13#c0#h0#u0#s#v4"]}"#))]
    #[case::open_isotope(mol_dsl!(r#"{:atoms ["C#c0#h4#u0#s"]}"#))]
    #[case::isotope_set(mol_dsl!(r#"{:atoms ["C#i{12,13}#c0#h4#u0#s"]}"#))]
    #[case::isotope_variable(mol_dsl!(r#"{:atoms ["C#i?mass#c0#h4#u0#s"]}"#))]
    #[case::open_charge(mol_dsl!(r#"{:atoms ["C#i=#h4#u0#s"]}"#))]
    #[case::open_spin(mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#))]
    #[case::different_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h2#n1#u0#s"]}"#))]
    #[case::open_h(mol_dsl!(r#"{:atoms ["C#i=#c0#u0#s"]}"#))]
    #[case::nonliteral_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h{2,4}#u0#s"]}"#))]
    #[case::no_match(mol_dsl!(r#"{:atoms ["C#i=#c0#h0#u0#s#v5"]}"#))]
    #[case::open_aromatic(mol_dsl!(r#"{:atoms ["C#i=#c0#h1#u0#s#a+"]}"#))]
    #[case::aromatic_carbon(mol_dsl!(r#"{:atoms ["C#i=#c0#h1#n0#u0#s#v2#a1"]}"#))]
    #[case::aromatic_nitrogen(mol_dsl!(r#"{:atoms ["N#i=#c0#h1#n0#u0#s#v2#a2"]}"#))]
    #[case::tetrahedral_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h1#n0#u0#s#v3#T0"]}"#))]
    #[case::tetrahedral_zero_h(mol_dsl!(r#"{:atoms ["C#i=#c0#h0#n0#u0#s#v4#T1"]}"#))]
    #[case::tetrahedral_lone_pair(mol_dsl!(r#"{:atoms ["S#i=#c0#h0#n1#u0#s#v4#T0"]}"#))]
    fn test_valence_resolver_project_identity(
        #[case] mut molecule: Molecule,
        #[values(false, true)] typing: bool,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let model = if typing {
            ValenceModel::default()
        } else {
            ValenceModel::smiles()
        };
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, policy),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::carbon(
        r#"{
        :atoms ["C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"
                "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"
                "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#
    )]
    #[case::nitrogen(
        r#"{
        :atoms ["N#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"
                "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h1#n0#u0#s"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]
    }"#
    )]
    fn test_valence_resolver_project_aromatic_identity(
        #[case] input: &str,
        #[values(false, true)] typing: bool,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let model = if typing {
            ValenceModel::default()
        } else {
            ValenceModel::smiles()
        };
        let original = mol_dsl!(input);
        let mut molecule = original.clone();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, policy),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }
}
