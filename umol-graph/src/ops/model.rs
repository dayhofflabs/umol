//! Chemistry model: top-level configuration consumed by the resolver and
//! validator engines.
//!
//! `ChemistryModel` wraps a `ValenceModel` (atom-typing or counts) and an
//! `AromaticityModel` (Hückel rule, HMO, or Clar). Engines read this; engines
//! and configs are kept as distinct types so multiple engine instances can
//! share one model.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, Lattice, MoleculeAst,
    MulticenterValenceAst, ValueAst,
};
use umol_shared::element::Element;

use crate::ops::valence::{
    AtomTypeRegistry, Mismatch, NormalValenceTable, ValenceEntry, ValenceInvariants, ValenceTable
};

#[derive(Debug, Clone)]
pub struct ChemistryModel {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
}

impl Default for ChemistryModel {
    fn default() -> Self {
        Self {
            valence: ValenceModel::AtomTyping {
                registry: AtomTypeRegistry::default_registry().clone(),
                normal_valence: NormalValenceTable::default_table().clone(),
            },
            aromaticity: AromaticityModel::daylight(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum ValenceModel {
    AtomTyping {
        registry: AtomTypeRegistry,
        normal_valence: NormalValenceTable,
    },
    Counts {
        table: ValenceTable,
        normal_valence: NormalValenceTable,
        allow_implicit_hydrogens: bool,
    },
}

impl ValenceModel {
    /// Universal per-atom electron-conservation invariants. Shared between
    /// `AtomTyping` and `Counts`: the equation is the same physics; the
    /// models differ only in how they narrow under-determined atoms.
    pub fn invariants(&self) -> ValenceInvariants {
        ValenceInvariants
    }

    /// Narrow every atom in `ast` per the model. Stops at the first atom
    /// that cannot be uniquely resolved.
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), ValenceResolveError> {
        for i in 0..ast.atoms().count() as u32 {
            self.resolve_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    /// Narrow `atom` per the model. AtomTyping matches against the registry;
    /// Counts enumerates the invariant equation, narrowing via min-unpaired
    /// then max-lone-pairs over the `covalence_set` / `aromatic_valence_set`
    /// trial space.
    pub fn resolve_atom(
        &self,
        ast: &mut MoleculeAst,
        atom: AtomId,
    ) -> Result<(), ValenceResolveError> {
        match self {
            Self::AtomTyping { registry, .. } => atom_typing_resolve_atom(registry, ast, atom),
            Self::Counts { table, .. } => counts_resolve_atom(table, ast, atom),
        }
    }

    /// Check every atom in `ast` per the model. Stops at the first atom
    /// that violates.
    pub fn validate(&self, ast: &MoleculeAst) -> Result<(), ValenceValidateError> {
        for i in 0..ast.atoms().count() as u32 {
            self.validate_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    /// Check `atom` per the model. AtomTyping requires at least one matching
    /// registry pattern; Counts evaluates the invariant equation.
    pub fn validate_atom(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
    ) -> Result<(), ValenceValidateError> {
        match self {
            Self::AtomTyping { registry, .. } => atom_typing_validate_atom(registry, ast, atom),
            Self::Counts { .. } => {
                ValenceInvariants::check_atom(ast, atom).map_err(ValenceValidateError::CountsMismatch)
            }
        }
    }
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceResolveError {
    #[error("atom {atom:?}: no atom-typing registry pattern matches (element {element}, charge {charge:?})")]
    AtomTypingNoMatch {
        atom: AtomId,
        element: Element,
        charge: ValueAst,
    },
    #[error("atom {atom:?}: {count} atom-typing registry patterns match (element {element}); expected exactly one")]
    AtomTypingAmbiguous {
        atom: AtomId,
        element: Element,
        count: usize,
    },
    #[error("atom {atom:?}: no candidate state satisfies the Counts model (element {element})")]
    CountsNoMatch { atom: AtomId, element: Element },
    #[error("atom {atom:?}: {count} candidate states survive min-unpaired and max-lone-pairs (element {element})")]
    CountsAmbiguous {
        atom: AtomId,
        element: Element,
        count: usize,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceValidateError {
    #[error("atom {atom:?}: no atom-typing registry pattern matches (element {element}, charge {charge:?})")]
    AtomTypingNoMatch {
        atom: AtomId,
        element: Element,
        charge: ValueAst,
    },
    #[error(transparent)]
    CountsMismatch(Mismatch),
}

fn counts_resolve_atom(
    table: &ValenceTable,
    ast: &mut MoleculeAst,
    atom_id: AtomId,
) -> Result<(), ValenceResolveError> {
    let view = ast.atom(atom_id);
    if view.ast.is_ground() {
        return Ok(());
    }
    let Some(element) = view.element().as_lit() else {
        return Ok(());
    };

    let entry = table.entry(element);
    let aromatic_trials = counts_aromatic_trials(&view, entry);
    let saved_atom = view.ast.clone();

    let mut all_candidates: Vec<AtomAst> = Vec::new();
    for av_trial in &aromatic_trials {
        {
            let atom_mut = ast.atom_mut(atom_id).ast;
            *atom_mut = saved_atom.clone();
            atom_mut
                .constraints
                .add(AtomConstraint::AromaticValence(av_trial.clone()));
        }
        all_candidates.extend(ValenceInvariants::solve_atom(ast, atom_id));
    }
    {
        let atom_mut = ast.atom_mut(atom_id).ast;
        *atom_mut = saved_atom;
    }

    if let Some(e) = entry {
        if !e.covalence_set.is_empty() {
            let topology_v = ast
                .atom(atom_id)
                .valence()
                .as_lit()
                .and_then(|n| u8::try_from(n).ok())
                .unwrap_or(0);
            all_candidates.retain(|c| {
                c.implicit_hydrogens
                    .as_lit()
                    .and_then(|h| u8::try_from(h).ok())
                    .is_some_and(|h| e.covalence_set.contains(&(topology_v + h)))
            });
        }
    }

    if let Some(min_u) = all_candidates
        .iter()
        .filter_map(|c| c.spin.unpaired.as_lit())
        .min()
    {
        all_candidates.retain(|c| c.spin.unpaired.as_lit() == Some(min_u));
    }
    if let Some(max_n) = all_candidates
        .iter()
        .filter_map(|c| c.lone_pairs.as_lit())
        .max()
    {
        all_candidates.retain(|c| c.lone_pairs.as_lit() == Some(max_n));
    }

    match all_candidates.len() {
        0 => Err(ValenceResolveError::CountsNoMatch {
            atom: atom_id,
            element,
        }),
        1 => {
            let cand = all_candidates.into_iter().next().unwrap();
            ast.atom_mut(atom_id).ast.narrow_from(&cand);
            Ok(())
        }
        n => Err(ValenceResolveError::CountsAmbiguous {
            atom: atom_id,
            element,
            count: n,
        }),
    }
}

fn counts_aromatic_trials(
    view: &AtomView<'_>,
    entry: Option<&ValenceEntry>,
) -> Vec<AromaticValenceAst> {
    let aromatic_constraint = view.constraints().aromatic_valence();
    let is_aromatic = view.is_in_aromatic_system()
        || AromaticValenceAst::aromatic(ValueAst::Undetermined).matches(&aromatic_constraint);
    if !is_aromatic {
        return vec![AromaticValenceAst::NotAromatic];
    }
    let av_set = entry
        .map(|e| e.aromatic_valence_set.as_slice())
        .unwrap_or(&[]);
    if av_set.is_empty() {
        return vec![
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
            AromaticValenceAst::Aromatic(ValueAst::Lit(2)),
        ];
    }
    av_set
        .iter()
        .map(|&v| AromaticValenceAst::Aromatic(ValueAst::Lit(v as i64)))
        .collect()
}

fn atom_typing_resolve_atom(
    registry: &AtomTypeRegistry,
    ast: &mut MoleculeAst,
    atom_id: AtomId,
) -> Result<(), ValenceResolveError> {
    let view = ast.atom(atom_id);
    if view.ast.is_ground() {
        return Ok(());
    }
    let Some(element) = view.element().as_lit() else {
        return Ok(());
    };
    let Some(prepared) = atom_typing_prepared(&view) else {
        return Ok(());
    };
    let charge_key = prepared.charge.as_lit().and_then(|n| i8::try_from(n).ok());
    let compatibles: Vec<&AtomAst> = registry
        .lookup(element, charge_key)
        .iter()
        .filter(|pat| prepared.matches(pat))
        .collect();
    match compatibles.len() {
        0 => Err(ValenceResolveError::AtomTypingNoMatch {
            atom: atom_id,
            element,
            charge: view.ast.charge.clone(),
        }),
        1 => {
            let cand = compatibles[0].clone();
            ast.atom_mut(atom_id).ast.narrow_from(&cand);
            Ok(())
        }
        n => Err(ValenceResolveError::AtomTypingAmbiguous {
            atom: atom_id,
            element,
            count: n,
        }),
    }
}

fn atom_typing_validate_atom(
    registry: &AtomTypeRegistry,
    ast: &MoleculeAst,
    atom_id: AtomId,
) -> Result<(), ValenceValidateError> {
    let view = ast.atom(atom_id);
    let Some(element) = view.element().as_lit() else {
        return Ok(());
    };
    let Some(prepared) = atom_typing_prepared(&view) else {
        return Ok(());
    };
    let charge_key = prepared.charge.as_lit().and_then(|n| i8::try_from(n).ok());
    let any_match = registry
        .lookup(element, charge_key)
        .iter()
        .any(|pat| prepared.matches(pat));
    if any_match {
        Ok(())
    } else {
        Err(ValenceValidateError::AtomTypingNoMatch {
            atom: atom_id,
            element,
            charge: view.ast.charge.clone(),
        })
    }
}

/// Project an atom into the shape registry patterns are written against:
/// add explicit constraints for the topology-derived valence, donated /
/// accepted pairs, aromatic membership, and multicenter sum so that
/// `prepared.matches(pattern)` filters via `AtomConstraints::matches`. No
/// implicit-hydrogen pre-fill — leftover ambiguity surfaces as
/// `AtomTypingAmbiguous`.
fn atom_typing_prepared(view: &AtomView<'_>) -> Option<AtomAst> {
    let valence = view.valence().as_lit().and_then(|n| u8::try_from(n).ok())?;
    let donated = view
        .donated_pairs()
        .as_lit()
        .and_then(|n| u8::try_from(n).ok())?;
    let accepted = view
        .accepted_pairs()
        .as_lit()
        .and_then(|n| u8::try_from(n).ok())?;
    let mut prepared = view.ast.clone();
    prepared
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(valence as i64)));
    prepared
        .constraints
        .add(AtomConstraint::DonatedPairs(ValueAst::Lit(donated as i64)));
    prepared
        .constraints
        .add(AtomConstraint::AcceptedPairs(ValueAst::Lit(accepted as i64)));
    if view.is_in_aromatic_system() {
        if let Some(pi) = view
            .aromatic_valence()
            .as_lit()
            .and_then(|n| u8::try_from(n).ok())
        {
            prepared
                .constraints
                .add(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(
                    ValueAst::Lit(pi as i64),
                )));
        }
    } else {
        let declared_aromatic = AromaticValenceAst::aromatic(ValueAst::Undetermined)
            .matches(&prepared.constraints.aromatic_valence());
        if !declared_aromatic {
            prepared.constraints.add(AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic,
            ));
        }
    }
    if let Some(mc) = view
        .multicenter_valence()
        .as_lit()
        .and_then(|n| u8::try_from(n).ok())
    {
        let mc_constraint = if mc == 0 {
            MulticenterValenceAst::NotMulticenter
        } else {
            MulticenterValenceAst::Multicenter(ValueAst::Lit(mc as i64))
        };
        prepared
            .constraints
            .add(AtomConstraint::MulticenterValence(mc_constraint));
    }
    Some(prepared)
}

#[derive(Debug, Clone)]
pub enum AromaticityModel {
    HueckelRule {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
    Hmo {
        scope: ElementScope,
        stabilization_threshold: f64,
    },
    Clar {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
}

impl AromaticityModel {
    /// Daylight (SMILES) aromaticity scope: C, N, O, S, Se, As.
    pub fn daylight() -> Self {
        Self::HueckelRule {
            scope: ElementScope::AllowList(vec![
                Element::C,
                Element::N,
                Element::O,
                Element::S,
                Element::Se,
                Element::As,
            ]),
            ring_limits: RingLimits::default(),
        }
    }

    /// MDL (MOL/SDF) aromaticity scope: C and N only, minimum ring size 6.
    pub fn mdl() -> Self {
        Self::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C, Element::N]),
            ring_limits: RingLimits {
                min_ring_size: 6,
                ..RingLimits::default()
            },
        }
    }

    /// Permissive aromaticity scope: any element.
    pub fn permissive() -> Self {
        Self::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        }
    }
}

/// Elements eligible for aromaticity perception.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementScope {
    Any,
    AllowList(Vec<Element>),
}

impl ElementScope {
    pub fn contains(&self, element: Element) -> bool {
        match self {
            Self::Any => true,
            Self::AllowList(list) => list.contains(&element),
        }
    }
}

/// Ring-size and fused-ring search bounds for ring-based aromaticity perception.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingLimits {
    pub min_ring_size: usize,
    pub max_ring_size: usize,
    pub include_fused: bool,
    pub max_fused_combination: usize,
    pub max_fused_search: usize,
}

impl Default for RingLimits {
    fn default() -> Self {
        Self {
            min_ring_size: 3,
            max_ring_size: 22,
            include_fused: true,
            max_fused_combination: 6,
            max_fused_search: 10_000,
        }
    }
}

/// Setup-time errors loading model data (TOML registries / valence tables).
/// Distinct from the per-engine `*Error` types that surface at resolve time.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ConfigError {
    #[error("invalid atom type registry: {0}")]
    InvalidAtomTypeRegistry(String),
    #[error("invalid valence table: {0}")]
    InvalidValenceTable(String),
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::{mol, mol_ground};
    use umol_shared::element::Element;

    use super::*;
    use crate::registry;

    #[rstest]
    #[case::counts_methane(
        ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            normal_valence: NormalValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        },
        mol!(r#"{:atoms ["C #c0 #n0"] :bonds []}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    #[case::atom_typing_methane(
        ValenceModel::AtomTyping {
            registry: registry!["C#c0#h4#n0#u0"],
            normal_valence: NormalValenceTable::default_table().clone(),
        },
        mol!(r#"{:atoms ["C #c0"] :bonds []}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    fn test_valence_model_resolve_atom(
        #[case] model: ValenceModel,
        #[case] mut ast: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected_implicit_h: ValueAst,
    ) {
        model.resolve_atom(&mut ast, atom).unwrap();
        assert_eq!(ast.atom(atom).ast.implicit_hydrogens, expected_implicit_h);
    }

    #[rstest]
    #[case::atom_typing_no_match(
        ValenceModel::AtomTyping {
            registry: registry!["C#c0#h4#n0#u0"],
            normal_valence: NormalValenceTable::default_table().clone(),
        },
        mol!(r#"{:atoms ["Cl #c0"] :bonds []}"#),
        AtomId(0),
        ValenceResolveError::AtomTypingNoMatch {
            atom: AtomId(0),
            element: Element::Cl,
            charge: ValueAst::Lit(0),
        },
    )]
    fn test_valence_model_resolve_atom_error(
        #[case] model: ValenceModel,
        #[case] mut ast: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValenceResolveError,
    ) {
        assert_eq!(model.resolve_atom(&mut ast, atom), Err(expected));
    }

    #[rstest]
    #[case::counts_methane(
        ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            normal_valence: NormalValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        },
        mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        AtomId(0),
    )]
    #[case::atom_typing_methane(
        ValenceModel::AtomTyping {
            registry: registry!["C#c0#h4#n0#u0"],
            normal_valence: NormalValenceTable::default_table().clone(),
        },
        mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        AtomId(0),
    )]
    fn test_valence_model_validate_atom(
        #[case] model: ValenceModel,
        #[case] ast: MoleculeAst,
        #[case] atom: AtomId,
    ) {
        model.validate_atom(&ast, atom).unwrap();
    }

    #[rstest]
    #[case::counts_orbital_mismatch(
        ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            normal_valence: NormalValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        },
        mol_ground!(r#"{:atoms ["C #h99"] :bonds []}"#),
        AtomId(0),
        ValenceValidateError::CountsMismatch(Mismatch::OrbitalCount {
            atom: AtomId(0),
            orbital_count: 198,
            electron_count: 103,
        }),
    )]
    #[case::atom_typing_no_match(
        ValenceModel::AtomTyping {
            registry: registry!["C#c0#h4#n0#u0"],
            normal_valence: NormalValenceTable::default_table().clone(),
        },
        mol_ground!(r#"{:atoms ["Cl #h0"] :bonds []}"#),
        AtomId(0),
        ValenceValidateError::AtomTypingNoMatch {
            atom: AtomId(0),
            element: Element::Cl,
            charge: ValueAst::Lit(0),
        },
    )]
    fn test_valence_model_validate_atom_error(
        #[case] model: ValenceModel,
        #[case] ast: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValenceValidateError,
    ) {
        assert_eq!(model.validate_atom(&ast, atom), Err(expected));
    }

    #[rstest]
    #[case::counts_ethane(
        ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            normal_valence: NormalValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        },
        mol!(r#"{:atoms ["C #c0 #n0" "C #c0 #n0"] :bonds [[0 1 "1"]]}"#),
        [ValueAst::Lit(3), ValueAst::Lit(3)],
    )]
    fn test_valence_model_resolve(
        #[case] model: ValenceModel,
        #[case] mut ast: MoleculeAst,
        #[case] expected_implicit_h: [ValueAst; 2],
    ) {
        model.resolve(&mut ast).unwrap();
        for (i, expected) in expected_implicit_h.iter().enumerate() {
            assert_eq!(&ast.atom(AtomId(i as u32)).ast.implicit_hydrogens, expected);
        }
    }

    #[rstest]
    #[case::counts_ethane(
        ValenceModel::Counts {
            table: ValenceTable::default_table().clone(),
            normal_valence: NormalValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        },
        mol_ground!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#),
    )]
    fn test_valence_model_validate(#[case] model: ValenceModel, #[case] ast: MoleculeAst) {
        model.validate(&ast).unwrap();
    }

    #[test]
    fn test_chemistry_model_default() {
        let model = ChemistryModel::default();
        assert!(matches!(model.valence, ValenceModel::AtomTyping { .. }));
        assert!(matches!(
            model.aromaticity,
            AromaticityModel::HueckelRule { .. }
        ));
    }

    #[rstest]
    #[case::any(ElementScope::Any, Element::U, true)]
    #[case::allow_match(ElementScope::AllowList(vec![Element::C]), Element::C, true)]
    #[case::allow_miss(ElementScope::AllowList(vec![Element::C]), Element::N, false)]
    fn test_element_scope_contains(
        #[case] scope: ElementScope,
        #[case] element: Element,
        #[case] expected: bool,
    ) {
        assert_eq!(scope.contains(element), expected);
    }

    #[test]
    fn test_aromaticity_model_daylight_scope() {
        match AromaticityModel::daylight() {
            AromaticityModel::HueckelRule { scope, .. } => {
                assert!(scope.contains(Element::C));
                assert!(scope.contains(Element::N));
                assert!(!scope.contains(Element::B));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }

    #[test]
    fn test_aromaticity_model_mdl_min_ring_size() {
        match AromaticityModel::mdl() {
            AromaticityModel::HueckelRule {
                ring_limits, scope, ..
            } => {
                assert_eq!(ring_limits.min_ring_size, 6);
                assert!(scope.contains(Element::N));
                assert!(!scope.contains(Element::O));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }

    #[test]
    fn test_aromaticity_model_permissive_scope() {
        match AromaticityModel::permissive() {
            AromaticityModel::HueckelRule { scope, .. } => {
                assert!(matches!(scope, ElementScope::Any));
            }
            other => panic!("expected HueckelRule, got {:?}", other),
        }
    }
}
