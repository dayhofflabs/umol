//! Valence resolver. Dispatches between atom-typing and counts strategies
//! defined in [`crate::ops::valence`].

use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondId,
    Molecule, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use crate::ops::model::{ValenceCandidateSource, ValenceModel, ValenceTieBreak};
use crate::ops::resolve::ResolveState;
use crate::ops::valence::compare::compare_by_key;
use crate::ops::valence::{
    AtomTypingError, AtomTypingValence, CountsError, CountsValence, ResolveReport,
};
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

/// Failures to reduce atom fields to a reconstructible ordinary valence description.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceProjectError {
    #[error("atom {atom:?} has non-concrete fields")]
    NonConcreteAtom { atom: AtomId },
    #[error("bond {bond:?} has a non-literal order")]
    NonLiteralBondOrder { bond: BondId },
    #[error("aromatic systems require aromatic valence reconstruction")]
    AromaticSystems,
    #[error("atom {atom:?} has aromatic evidence outside ordinary valence reconstruction")]
    AromaticAtom { atom: AtomId },
    #[error("dative bonds are outside ordinary valence reconstruction")]
    DativeBonds,
    #[error("multicenter bonds are outside ordinary valence reconstruction")]
    MulticenterBonds,
    #[error("admission did not complete the fields of atom {atom:?}")]
    IncompleteAtom { atom: AtomId },
    #[error("the selected completion does not recover the fields of atom {atom:?}")]
    AtomMismatch { atom: AtomId },
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

    /// Reduce ordinary atom fields while establishing their recovery by this valence source.
    ///
    /// Retains element, isotope, charge, implicit H count, all bonds, entities, and assertions.
    /// Opens lone pairs and both unpaired-electron fields on a private copy, admits candidates,
    /// and applies `tie_break` to distinct completions. Candidate constraints remain solver
    /// evidence, as in forward resolution; only recovered inherent atom fields are compared.
    ///
    /// # Semantic properties
    ///
    /// A determined result publishes the reduced molecule only when the selected completion
    /// recovers every source atom field exactly. All other outcomes preserve the input exactly.
    /// This is a valence-phase guarantee, not a complete resolver or external-format roundtrip:
    /// stereo, other entities, and molecule-level assertions are preserved without interpretation.
    ///
    /// # Errors
    ///
    /// Rejects non-concrete source atoms, non-literal bond orders, aromatic evidence, dative or
    /// multicenter bonds, incomplete admissions, and a selected state different from the source.
    /// Admission contradictions and unresolved candidate plurality use `Solution`.
    pub fn project(
        &self,
        molecule: &mut Molecule,
        tie_break: ValenceTieBreak,
    ) -> Result<Solution<ResolveReport, ValenceContradiction>, ValenceProjectError> {
        if molecule.has_aromatic_systems() {
            return Err(ValenceProjectError::AromaticSystems);
        }
        if molecule.has_dative_bonds() {
            return Err(ValenceProjectError::DativeBonds);
        }
        if molecule.has_multicenter_bonds() {
            return Err(ValenceProjectError::MulticenterBonds);
        }
        for atom in molecule.atoms().iter() {
            if !atom.attributes.is_concrete() {
                return Err(ValenceProjectError::NonConcreteAtom { atom: atom.id });
            }
            if !matches!(
                atom.constraints()
                    .asserted_complete(AtomConstraintKey::AromaticValence),
                Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::NotAromatic
                ))
            ) {
                return Err(ValenceProjectError::AromaticAtom { atom: atom.id });
            }
        }
        for bond in molecule.bonds().iter() {
            if bond.order().as_lit().is_none() {
                return Err(ValenceProjectError::NonLiteralBondOrder { bond: bond.id });
            }
        }
        let mut projected = molecule.clone();
        projected.modify_atoms(|mut atom| {
            atom.lone_pairs = NumForm::Undetermined;
            atom.unpaired_electrons = UnpairedElectronsForm::default();
            atom
        });
        let state = match self.admit(&projected) {
            Ok(Solution::Determined(state)) => state,
            Ok(Solution::Underdetermined(state)) => {
                return Ok(Solution::Underdetermined(state.to_report()));
            }
            Ok(Solution::Contradictory(error)) => return Ok(Solution::Contradictory(error)),
            Err(error) => match error {},
        };
        let mut report = ResolveReport::default();
        for atom in molecule.atoms().iter() {
            let candidates = state
                .completions
                .get(atom.id)
                .ok_or(ValenceProjectError::IncompleteAtom { atom: atom.id })?;
            if candidates.iter().any(|candidate| !candidate.is_concrete()) {
                return Err(ValenceProjectError::IncompleteAtom { atom: atom.id });
            }
            let key = tie_break.key();
            let selected = if candidates.len() == 1 {
                Some(&candidates[0])
            } else if key.is_empty() {
                None
            } else {
                let best = candidates
                    .iter()
                    .max_by(|a, b| compare_by_key(key, a, b))
                    .expect("admission entries are nonempty");
                (candidates
                    .iter()
                    .filter(|candidate| compare_by_key(key, candidate, best).is_eq())
                    .count()
                    == 1)
                    .then_some(best)
            };
            let Some(selected) = selected else {
                report
                    .unresolved
                    .insert(atom.id, candidates.iter().cloned().collect());
                continue;
            };
            let AtomForm {
                element,
                isotope_mass,
                charge,
                implicit_hydrogens,
                lone_pairs,
                unpaired_electrons,
                constraints: _,
            } = selected;
            let source = atom.attributes;
            if element != &source.element
                || isotope_mass != &source.isotope_mass
                || charge != &source.charge
                || implicit_hydrogens != &source.implicit_hydrogens
                || lone_pairs != &source.lone_pairs
                || unpaired_electrons != &source.unpaired_electrons
            {
                return Err(ValenceProjectError::AtomMismatch { atom: atom.id });
            }
            if candidates.len() > 1 {
                report.tie_breaks.push(atom.id);
            }
        }
        if !report.unresolved.is_empty() {
            return Ok(Solution::Underdetermined(report));
        }
        *molecule = projected;
        Ok(Solution::Determined(report))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::rstest;
    use smallvec::smallvec;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomConstraintForm, AtomId, MoleculeEntries};
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::ChemistryModel;
    use crate::ops::resolve::Resolver;
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
    #[case::methane(atom_dsl!("C#i=#c0#h4#n0#u0#s"), atom_dsl!("C#i=#c0#h4"))]
    #[case::methyl(atom_dsl!("C#i=#c0#h3#n0#u1#s2"), atom_dsl!("C#i=#c0#h3"))]
    #[case::carbanion(atom_dsl!("C#i=#c-#h3#n1#u0#s"), atom_dsl!("C#i=#c-#h3"))]
    #[case::carbocation(atom_dsl!("C#i=#c+#h3#n0#u0#s"), atom_dsl!("C#i=#c+#h3"))]
    #[case::ammonia(atom_dsl!("N#i=#c0#h3#n1#u0#s"), atom_dsl!("N#i=#c0#h3"))]
    #[case::ammonium(atom_dsl!("N#i=#c+#h4#n0#u0#s"), atom_dsl!("N#i=#c+#h4"))]
    #[case::water(atom_dsl!("O#i=#c0#h2#n2#u0#s"), atom_dsl!("O#i=#c0#h2"))]
    #[case::hydroxide(atom_dsl!("O#i=#c-#h1#n3#u0#s"), atom_dsl!("O#i=#c-#h1"))]
    #[case::fluorane(atom_dsl!("F#i=#c0#h1#n3#u0#s"), atom_dsl!("F#i=#c0#h1"))]
    #[case::chloride(atom_dsl!("Cl#i=#c-#h0#n4#u0#s"), atom_dsl!("Cl#i=#c-#h0"))]
    #[case::phosphonium(atom_dsl!("P#i=#c+#h4#n0#u0#s"), atom_dsl!("P#i=#c+#h4"))]
    #[case::sulfane(atom_dsl!("S#i=#c0#h2#n2#u0#s"), atom_dsl!("S#i=#c0#h2"))]
    #[case::borane(atom_dsl!("B#i=#c0#h3#n0#u0#s"), atom_dsl!("B#i=#c0#h3"))]
    #[case::silane(atom_dsl!("Si#i=#c0#h4#n0#u0#s"), atom_dsl!("Si#i=#c0#h4"))]
    fn test_valence_resolver_project(
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
        #[case] source: AtomForm,
        #[case] expected: AtomForm,
    ) {
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![source],
            ..Default::default()
        });
        let original = molecule.clone();
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![expected],
            ..Default::default()
        });
        assert_eq!(
            ValenceResolver::new(&valence).project(&mut molecule, valence.tie_break),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, expected);
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::ethene(mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h2"] :bonds [[0 1 "2"]]}"#), mol_dsl!(r#"{:atoms ["C#i=#c0#h2" "C#i=#c0#h2"] :bonds [[0 1 "2#c0#u0#s"]]}"#))]
    #[case::explicit_hydrogen(mol_dsl_concrete!(r#"{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}"#), mol_dsl!(r#"{:atoms ["C#i=#c0#h3" "H#i=#c0#h0"] :bonds [[0 1 "1#c0#u0#s"]]}"#))]
    #[case::assertion(mol_dsl_concrete!(r#"{:atoms ["C#h4#D0"]}"#), mol_dsl!(r#"{:atoms ["C#i=#c0#h4#D0"]}"#))]
    fn test_valence_resolver_project_preservation(
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let model = ValenceModel::smiles();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, model.tie_break),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::non_concrete(mol_dsl!(r#"{:atoms ["C"]}"#), ValenceProjectError::NonConcreteAtom {atom:AtomId(0)})]
    #[case::pairing(mol_dsl_concrete!(r#"{:atoms ["C#h2#n0#u2#s3"]}"#), ValenceProjectError::AtomMismatch {atom:AtomId(0)})]
    #[case::late_mismatch(mol_dsl_concrete!(r#"{:atoms ["C#h4" "C#h2#n0#u2#s3"]}"#), ValenceProjectError::AtomMismatch {atom:AtomId(1)})]
    #[case::spin(mol_dsl_concrete!(r#"{:atoms ["C#h3#n0#u1#s1"]}"#), ValenceProjectError::AtomMismatch {atom:AtomId(0)})]
    #[case::aromatic(mol_dsl_concrete!(r#"{:atoms ["C#h4#a+"]}"#), ValenceProjectError::AromaticAtom {atom:AtomId(0)})]
    #[case::bond_order(mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "*"]]}"#), ValenceProjectError::NonLiteralBondOrder {bond:BondId(0)})]
    #[case::aromatic_bond(mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1#a"]]}"#), ValenceProjectError::AromaticAtom {atom:AtomId(0)})]
    #[case::dative(mol_dsl_concrete!(r#"{:atoms ["N#h3#n1" "B#h3"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#), ValenceProjectError::DativeBonds)]
    #[case::aromatic_system(mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h2"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]"}]}"#), ValenceProjectError::AromaticSystems)]
    #[case::multicenter(mol_dsl_concrete!(r#"{:atoms ["H" "B#h2" "H"] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,0,1]"}]}"#), ValenceProjectError::MulticenterBonds)]
    fn test_valence_resolver_project_error(
        #[case] mut molecule: Molecule,
        #[case] expected: ValenceProjectError,
    ) {
        let model = ValenceModel::smiles();
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, model.tie_break),
            Err(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_valence_resolver_project_contradiction() {
        let model = ValenceModel::smiles();
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h5"]}"#);
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, model.tie_break),
            Ok(Solution::Contradictory(ValenceContradiction::Counts(
                CountsError::NoMatch
            )))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::singlet(ValenceTieBreak::MostSaturated, atom_dsl!("C#i=#c0#h2#n1#u0#s"), Ok(Solution::Determined(ResolveReport {tie_breaks:vec![AtomId(0)], ..Default::default()})))]
    #[case::triplet(ValenceTieBreak::MostSaturated, atom_dsl!("C#i=#c0#h2#n0#u2#s3"), Err(ValenceProjectError::AtomMismatch {atom:AtomId(0)}))]
    #[case::strict(ValenceTieBreak::Strict, atom_dsl!("C#i=#c0#h2#n1#u0#s"), Ok(Solution::Underdetermined(ResolveReport {
        unresolved: AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#i=#c0#h2#n1#u0#s"), atom_dsl!("C#i=#c0#h2#n0#u2#s3")])]),
        tie_breaks:vec![],
    })))]
    fn test_valence_resolver_project_selection(
        #[case] tie_break: ValenceTieBreak,
        #[case] source: AtomForm,
        #[case] expected: Result<
            Solution<ResolveReport, ValenceContradiction>,
            ValenceProjectError,
        >,
    ) {
        let model = ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
            atom_dsl!("C#i*#c0#h2#n1#u0#s"),
            atom_dsl!("C#i*#c0#h2#n0#u2#s3"),
        ])));
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![source],
            ..Default::default()
        });
        let original = molecule.clone();
        let outcome = ValenceResolver::new(&model).project(&mut molecule, tie_break);
        assert_eq!(outcome, expected);
        if matches!(outcome, Ok(Solution::Determined(_))) {
            assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h2"]}"#));
            let model = ChemistryModel {
                valence: ValenceModel { tie_break, ..model },
                ..Default::default()
            };
            assert_eq!(
                Resolver::new(&model).resolve(&mut molecule),
                Ok(Solution::Determined(ResolveReport {
                    tie_breaks: vec![AtomId(0)],
                    ..Default::default()
                }))
            );
        }
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_valence_resolver_project_overlap(
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)]
        tie_break: ValenceTieBreak,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel {
                candidates: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
                    atom_dsl!("C#c0#h*#n0#u0#s"),
                    atom_dsl!("C#c0#h4#n0#u0#s"),
                ])))
                .candidates,
                tie_break,
            },
            ..Default::default()
        };
        let resolver = Resolver::new(&model);
        let mut molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#);
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        let original = molecule.clone();
        assert_eq!(
            resolver.valence.project(&mut molecule, tie_break),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#));
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::counts(
        ValenceModel::smiles(),
        Ok(Solution::Determined(ResolveReport::default()))
    )]
    #[case::default_registry(
        ValenceModel::default(),
        Ok(Solution::Determined(ResolveReport::default()))
    )]
    #[case::custom_registry(ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#i*#c0#h4#n0#u0#s")]))), Ok(Solution::Determined(ResolveReport::default())))]
    fn test_valence_resolver_project_isotope(
        #[case] valence: ValenceModel,
        #[case] expected: Result<
            Solution<ResolveReport, ValenceContradiction>,
            ValenceProjectError,
        >,
    ) {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let resolver = Resolver::new(&model);
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#i13#h4"]}"#);
        let original = molecule.clone();
        let outcome = resolver.valence.project(&mut molecule, resolver.tie_break);
        assert_eq!(outcome, expected);
        if matches!(outcome, Ok(Solution::Determined(_))) {
            assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i13#c0#h4"]}"#));
            assert_eq!(
                resolver.resolve(&mut molecule),
                Ok(Solution::Determined(ResolveReport::default()))
            );
        }
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_valence_resolver_project_incomplete() {
        let model =
            ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                "C#c0#h4"
            )])));
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#);
        let original = molecule.clone();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, model.tie_break),
            Err(ValenceProjectError::IncompleteAtom { atom: AtomId(0) })
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::low_spin(vec![atom_dsl!("C#i*#c0#h2#n0#u2#s1")], Ok(Solution::Determined(ResolveReport::default())))]
    #[case::high_spin(vec![atom_dsl!("C#i*#c0#h2#n0#u2#s3")], Err(ValenceProjectError::AtomMismatch {atom:AtomId(0)}))]
    #[case::spin_tie(vec![atom_dsl!("C#i*#c0#h2#n0#u2#s1"), atom_dsl!("C#i*#c0#h2#n0#u2#s3")], Ok(Solution::Underdetermined(ResolveReport {
        unresolved: AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#i=#c0#h2#n0#u2#s1"), atom_dsl!("C#i=#c0#h2#n0#u2#s3")])]),
        tie_breaks:vec![],
    })))]
    fn test_valence_resolver_project_spin(
        #[case] rows: Vec<AtomForm>,
        #[case] expected: Result<
            Solution<ResolveReport, ValenceContradiction>,
            ValenceProjectError,
        >,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms(rows)))
            },
            ..Default::default()
        };
        let resolver = Resolver::new(&model);
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h2#n0#u2#s1"]}"#);
        let original = molecule.clone();
        let outcome = resolver.valence.project(&mut molecule, resolver.tie_break);
        assert_eq!(outcome, expected);
        if matches!(outcome, Ok(Solution::Determined(_))) {
            assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h2"]}"#));
            assert_eq!(
                resolver.resolve(&mut molecule),
                Ok(Solution::Determined(ResolveReport::default()))
            );
        }
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    fn test_valence_resolver_project_identity(#[case] mut molecule: Molecule) {
        let original = molecule.clone();
        let model = ValenceModel::smiles();
        assert_eq!(
            ValenceResolver::new(&model).project(&mut molecule, model.tie_break),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }
}
