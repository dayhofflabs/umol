//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_graph_ir::ir::{AtomConstraintKey, Molecule};
use umol_utils::solution::Solution;

use crate::ops::model::{ValenceCandidateSource, ValenceModel};
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

/// Operational failures of the valence phase; currently uninhabited — the
/// phase produces no edits and runs no transactions.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceError {}

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
    /// produced. The chemistry verdict rides `Solution`; the operational
    /// channel is currently uninhabited.
    pub fn admit(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<super::ResolveState, ValenceContradiction>, ValenceError> {
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
                        return Ok(Solution::Underdetermined(super::ResolveState::default()));
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
        Ok(completions.map(|completions| super::ResolveState {
            completions,
            systems: Vec::new(),
            tie_breaks: Vec::new(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use smallvec::smallvec;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomConstraintForm, AtomId};
    use umol_graph_ir::{atom_dsl, mol_dsl};

    use super::super::ResolveState;
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
    #[case::counts(ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())))]
    #[case::atom_typing(ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        )]))))]
    fn test_valence_resolver_admit(#[case] model: ValenceModel) {
        let molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        let mut completions = AtomCompletions::new();
        completions.insert(AtomId(0), smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#a!")]);
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
}
