//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintKey, AtomHandle, AtomUpdate, Edits, IsotopeMassForm,
    Molecule, NumForm, TetrahedralStereoForm, TransactionError,
};
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

/// Operational failures applying implicit-H projection.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceProjectError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

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

    /// Elide recoverable implicit H under the supplied valence policy.
    ///
    /// Only Natural isotope composition, zero charge, and zero unpaired electrons permit
    /// elision. A literal H count becomes Undetermined when inference selects that count;
    /// missing or ambiguous evidence leaves it unchanged. Isotope projection follows this
    /// operation so that Natural is still available for this external-format requirement.
    /// Tetrahedral assertions retain H, including zero, before any inference lookup.
    /// Aromatic atoms retain positive H counts, using asserted aromaticity when present
    /// and aromatic-system membership otherwise.
    ///
    /// # Semantic properties
    ///
    /// Projection is idempotent and changes only implicit-H fields. Stored lone pairs do
    /// not affect elision. Atom-typing row reordering or repetition preserves the result.
    /// Successful projection returns Determined even when some H counts remain unresolved.
    ///
    /// # Errors
    ///
    /// Returns Transaction if applying H edits fails. Errors preserve the molecule exactly.
    pub fn project(
        &self,
        molecule: &mut Molecule,
        policy: ValenceTieBreak,
    ) -> Result<Solution<(), ValenceContradiction>, ValenceProjectError> {
        let mut edits = Edits::new();
        for atom in molecule.atoms().iter() {
            if matches!(
                atom.constraints().tetrahedral_stereo(),
                Some(TetrahedralStereoForm::Stereo(_))
            ) || !matches!(atom.attributes.isotope_mass, IsotopeMassForm::Natural)
                || atom.charge().as_lit() != Some(0)
                || atom.unpaired_electrons().count.as_lit() != Some(0)
            {
                continue;
            }
            let Some(hydrogens) = atom.implicit_hydrogens().as_lit() else {
                continue;
            };
            if hydrogens > 0
                && match atom.constraints().aromatic_valence() {
                    Some(AromaticValenceForm::Aromatic(_)) => true,
                    None => atom.is_in_aromatic_system(),
                    _ => false,
                }
            {
                continue;
            }
            let inferred = match self {
                Self::AtomTyping(resolver) => {
                    resolver.infer_implicit_hydrogens(molecule, atom.id, policy)
                }
                Self::Counts(resolver) => {
                    resolver.infer_implicit_hydrogens(molecule, atom.id, policy)
                }
            };
            if inferred == Some(hydrogens) {
                edits.update_atom(
                    AtomHandle::Id(atom.id),
                    atom.attributes,
                    &AtomUpdate {
                        implicit_hydrogens: Some(NumForm::Undetermined),
                        ..Default::default()
                    },
                );
            }
        }
        if !edits.is_empty() {
            let mut editor = molecule.edit();
            editor.transact(edits)?;
            *molecule = editor.build();
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use smallvec::smallvec;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomConstraintForm, AtomForm, AtomId, MoleculeEntries};
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
    #[case::methane(
        ValenceTieBreak::MostSaturated,
        r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#,
        r#"{:atoms ["C#i=#c0#n0#u0#s"]}"#
    )]
    #[case::fluorine(
        ValenceTieBreak::Strict,
        r#"{:atoms ["F#i=#c0#h1#n3#u0#s"]}"#,
        r#"{:atoms ["F#i=#c0#n3#u0#s"]}"#
    )]
    #[case::zero_h(
        ValenceTieBreak::Strict,
        r#"{:atoms ["C#i=#c0#h0#n0#u0#s#v4"]}"#,
        r#"{:atoms ["C#i=#c0#n0#u0#s#v4"]}"#
    )]
    #[case::aromatic_zero_h(
        ValenceTieBreak::Strict,
        r#"{:atoms ["C#i=#c0#h0#n0#u0#s#v3#a1"]}"#,
        r#"{:atoms ["C#i=#c0#n0#u0#s#v3#a1"]}"#
    )]
    #[case::lone_pairs(
        ValenceTieBreak::MostSaturated,
        r#"{:atoms ["C#i=#c0#h4#n2#u0#s"]}"#,
        r#"{:atoms ["C#i=#c0#n2#u0#s"]}"#
    )]
    #[case::ethane(
        ValenceTieBreak::MostSaturated,
        r#"{:atoms ["C#i=#c0#h3#n0#u0#s" "C#i=#c0#h3#n0#u0#s"] :bonds [[0 1 "1"]]}"#,
        r#"{:atoms ["C#i=#c0#n0#u0#s" "C#i=#c0#n0#u0#s"] :bonds [[0 1 "1"]]}"#
    )]
    #[case::mixed(
        ValenceTieBreak::MostSaturated,
        r#"{:atoms ["C#i=#c0#h4#u0#s" "C#i13#c0#h4#u0#s" "C#i=#c-#h3#u0#s" "C#i=#c0#h3#u1#s2"]}"#,
        r#"{:atoms ["C#i=#c0#u0#s" "C#i13#c0#h4#u0#s" "C#i=#c-#h3#u0#s" "C#i=#c0#h3#u1#s2"]}"#
    )]
    fn test_valence_resolver_project(
        #[case] policy: ValenceTieBreak,
        #[case] input: &str,
        #[case] expected: &str,
        #[values(false, true)] typing: bool,
    ) {
        let model = if typing {
            ValenceModel::default()
        } else {
            ValenceModel::smiles()
        };
        let resolver = ValenceResolver::new(&model);
        let mut molecule = mol_dsl!(input);
        assert_eq!(
            resolver.project(&mut molecule, policy),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, mol_dsl!(expected));
        let projected = molecule.clone();
        assert_eq!(
            resolver.project(&mut molecule, policy),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, projected);
    }

    #[rstest]
    #[case::charge("C#i=#c+#h3#u0#s")]
    #[case::radical("C#i=#c0#h3#u1#s2")]
    #[case::open_shell_singlet("C#i=#c0#h2#u2#s")]
    #[case::mass("C#i13#c0#h4#u0#s")]
    #[case::mass_zero_h("C#i13#c0#h0#u0#s#v4")]
    #[case::open_isotope("C#c0#h4#u0#s")]
    #[case::isotope_set("C#i{12,13}#c0#h4#u0#s")]
    #[case::isotope_variable("C#i?mass#c0#h4#u0#s")]
    #[case::open_charge("C#i=#h4#u0#s")]
    #[case::open_spin("C#i=#c0#h4")]
    #[case::different_h("C#i=#c0#h2#n1#u0#s")]
    #[case::open_h("C#i=#c0#u0#s")]
    #[case::nonliteral_h("C#i=#c0#h{2,4}#u0#s")]
    #[case::no_match("C#i=#c0#h0#u0#s#v5")]
    #[case::open_aromatic("C#i=#c0#h1#u0#s#a+")]
    #[case::aromatic_carbon("C#i=#c0#h1#n0#u0#s#v2#a1")]
    #[case::aromatic_nitrogen("N#i=#c0#h1#n0#u0#s#v2#a2")]
    #[case::tetrahedral_h("C#i=#c0#h1#n0#u0#s#v3#T0")]
    #[case::tetrahedral_zero_h("C#i=#c0#h0#n0#u0#s#v4#T1")]
    #[case::tetrahedral_lone_pair("S#i=#c0#h0#n1#u0#s#v4#T0")]
    fn test_valence_resolver_project_identity(
        #[case] input: &str,
        #[values(false, true)] typing: bool,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let model = if typing {
            ValenceModel::default()
        } else {
            ValenceModel::smiles()
        };
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!(input)],
            ..Default::default()
        });
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
