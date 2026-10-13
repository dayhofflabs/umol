//! Aromatize: Kekulé form → aromatic-system form.
//!
//! [`Aromatizer`] runs aromaticity perception against a Kekulé-form input —
//! atoms with explicit single/double bonds. Per-atom π contributions are derived
//! from localized bonds and lone pairs in the candidate ring context. Selected
//! contributions are transferred from those bonds and lone pairs into systems.
//! If the input molecule already carries one or more aromatic systems, this is a
//! no-op: re-aromatizing requires kekulizing first.

use std::{iter, mem};

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticSystemId, AromaticValence, AromaticValenceForm, AsLit, AtomConstraintForm,
    AtomConstraintKey, AtomHandle, AtomId, AtomUpdate, BondHandle, BondUpdate, Edits,
    ElectronCountsForm, ElementForm, Lattice, Molecule, MulticenterValenceForm, NumForm, RingSet,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    aromatic_valence_inconsistency, AromaticityConfig, AromaticityContradiction, AromaticityError,
    AromaticityInconsistency, AromaticityPerceiver,
};
use crate::ops::invariant::{ValenceInvariantsValidator, ValenceMismatch};
use crate::ops::model::AromaticityModel;
use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromatizeError {
    #[error("aromaticity setup: {0}")]
    Setup(#[from] AromaticityError),
    #[error("aromaticity contradiction: {0}")]
    Contradiction(#[from] AromaticityContradiction),
    #[error("aromatization input is underdetermined")]
    Underdetermined,
    #[error(
        "atom {atom:?} contradicts {constraint:?} in planned aromatic system {aromatic_system:?}"
    )]
    ContradictedAssertion {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("planned aromatic system {aromatic_system:?}: {mismatch}")]
    ValenceInvariant {
        aromatic_system: AromaticSystemId,
        mismatch: ValenceMismatch,
    },
}

/// Transfers localized π electrons into perceived aromatic systems.
///
/// Charge, hydrogens, spin, external bonds, and unrelated state are preserved.
/// A molecule already containing an aromatic system is unchanged. Fresh perception
/// need not recover the same system membership as an earlier Kekulé localization.
///
/// Required nonliteral contribution or transfer data produces
/// [`AromatizeError::Underdetermined`]; the molecule need not be wholly ground.
/// Planned valence contradictions are reported before mutation, so an in-place
/// failure leaves the input unchanged.
///
/// # Semantic properties
///
/// On inputs without stored systems, every added system agrees with its members' open evidence
/// and contains no bond asserting non-aromaticity. Non-aromatic atom evidence
/// excludes that atom from perception. An incompatible selected contribution or
/// an unclaimed positive atom or bond assertion is an aromaticity contradiction.
#[derive(Clone, Debug)]
pub struct Aromatizer {
    perception: AromaticityPerceiver,
    config: AromaticityConfig,
}

impl Aromatizer {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityConfig) -> Self {
        Self {
            perception: AromaticityPerceiver::new(model),
            config,
        }
    }

    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, AromatizeError> {
        let mut edits = Edits::new();
        if molecule.aromatic_systems().count() > 0 {
            return Ok(edits);
        }
        let rings = self.perception.candidate_rings(molecule, self.config);
        let mut contributions = Vec::with_capacity(if rings.count() == 0 {
            0
        } else {
            molecule.atoms().count()
        });
        let mut required_atoms = Vec::new();
        let mut mismatches = Vec::new();
        for atom in molecule.atoms().iter() {
            let evidence = atom
                .constraints()
                .evidenced(AtomConstraintKey::AromaticValence)
                .ok_or_else(|| {
                    AromaticityContradiction::from(aromatic_valence_inconsistency(
                        molecule,
                        atom.id(),
                    ))
                })?;
            if matches!(
                atom.attributes().constraints.aromatic_valence(),
                AromaticValenceForm::Aromatic(_)
            ) {
                required_atoms.push(atom.id());
            }
            if rings.count() == 0 {
                continue;
            }
            let contribution = if matches!(
                evidence,
                AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic)
            ) {
                AromaticValence::NotAromatic
            } else {
                match electrons_from_kekule(molecule, &rings, atom.id()) {
                    Err(AromatizeError::Underdetermined) if !rings.contains_atom(atom.id()) => {
                        AromaticValence::NotAromatic
                    }
                    result => result?,
                }
            };
            if let AromaticValence::Aromatic(electrons) = contribution {
                if !evidence.is_undetermined()
                    && evidence
                        .meet(&AtomConstraintForm::aromatic_valence(
                            AromaticValenceForm::aromatic(electrons),
                        ))
                        .is_none()
                {
                    mismatches.push(atom.id());
                }
            }
            contributions.push(contribution);
        }
        if rings.count() == 0 {
            if let Some(inconsistency) = self
                .perception
                .check_aromatic_systems(molecule, &required_atoms, &[])
                .into_iter()
                .next()
            {
                return Err(AromaticityContradiction::from(inconsistency).into());
            }
            return Ok(edits);
        }
        let derivation = self
            .perception
            .derive_from_rings(molecule, &rings, self.config, &required_atoms, |atom| {
                contributions[atom.index()]
            })?
            .into_decisive(AromatizeError::Underdetermined)?;
        if let Some(inconsistency) = derivation.inconsistencies.into_iter().next() {
            return Err(AromaticityContradiction::from(inconsistency).into());
        }
        for (index, (atoms, attributes)) in derivation.systems.into_iter().enumerate() {
            let aromatic_system = AromaticSystemId(index as u32);
            if !mismatches.is_empty() {
                if let Some(&atom) = atoms
                    .iter()
                    .find(|atom| mismatches.binary_search(atom).is_ok())
                {
                    return Err(AromaticityContradiction::from(
                        AromaticityInconsistency::AromaticValenceMismatch {
                            atom,
                            system: aromatic_system,
                        },
                    )
                    .into());
                }
            }
            let ElectronCountsForm::Lit(electrons) = &attributes.electrons else {
                unreachable!("perception supplies literal aromatic contributions");
            };
            for (&id, &electrons) in atoms.iter().zip(electrons) {
                let atom = molecule.atom(id);
                let mut planned = atom.attributes().clone();
                planned
                    .constraints
                    .remove(AtomConstraintKey::AromaticValence);
                let mut valence = 0;
                for neighbor in atom.neighbors() {
                    let bond = neighbor.bond();
                    let mut order = bond.order().as_lit().expect(
                        "Kekulé contribution recovery requires literal incident bond orders",
                    );
                    if order == 2 && atoms.binary_search(&neighbor.atom_id()).is_ok() {
                        order = 1;
                        if id < neighbor.atom_id() {
                            edits.update_bond(
                                BondHandle::Id(bond.id()),
                                bond.attributes(),
                                &BondUpdate {
                                    order: Some(NumForm::Lit(1)),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    valence += order;
                }
                if electrons == 2 {
                    let pairs = atom
                        .lone_pairs()
                        .as_lit()
                        .ok_or(AromatizeError::Underdetermined)?;
                    planned.lone_pairs = NumForm::Lit(pairs - 1);
                    edits.update_atom(
                        AtomHandle::Id(id),
                        atom.attributes(),
                        &AtomUpdate {
                            lone_pairs: Some(planned.lone_pairs.clone()),
                            ..Default::default()
                        },
                    );
                }
                for key in [
                    AtomConstraintKey::Valence,
                    AtomConstraintKey::AromaticValence,
                    AtomConstraintKey::MulticenterValence,
                    AtomConstraintKey::DonatedPairs,
                    AtomConstraintKey::AcceptedPairs,
                ] {
                    let derived = match key {
                        AtomConstraintKey::Valence => Some(AtomConstraintForm::valence(valence)),
                        AtomConstraintKey::AromaticValence => {
                            Some(AtomConstraintForm::aromatic_valence(
                                AromaticValenceForm::aromatic(electrons),
                            ))
                        }
                        _ => atom.constraints().derived_complete(key),
                    }
                    .ok_or(AromatizeError::ValenceInvariant {
                        aromatic_system,
                        mismatch: ValenceMismatch::ContradictoryConstraint { atom_id: id, key },
                    })?;
                    let value = planned
                        .constraints
                        .asserted(key)
                        .and_then(|asserted| asserted.meet(&derived))
                        .ok_or_else(|| AromatizeError::ContradictedAssertion {
                            aromatic_system,
                            atom: id,
                            constraint: planned
                                .constraints
                                .get(key)
                                .cloned()
                                .unwrap_or_else(|| derived.clone()),
                        })?;
                    match value {
                        AtomConstraintForm::DonatedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::AcceptedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::MulticenterValence(
                            MulticenterValenceForm::NotMulticenter,
                        ) => {
                            planned.constraints.remove(key);
                        }
                        _ => planned.constraints.set(value),
                    }
                }
                if let Solution::Contradictory(mismatch) = ValenceInvariantsValidator
                    .validate_atom(&planned)
                    .unwrap_or_else(|error| match error {})
                {
                    let mismatch = match mismatch {
                        ValenceMismatch::OrbitalCount {
                            orbital_count,
                            electron_count,
                            ..
                        } => ValenceMismatch::OrbitalCount {
                            atom_id: id,
                            orbital_count,
                            electron_count,
                        },
                        ValenceMismatch::ContradictoryConstraint { key, .. } => {
                            ValenceMismatch::ContradictoryConstraint { atom_id: id, key }
                        }
                    };
                    return Err(AromatizeError::ValenceInvariant {
                        aromatic_system,
                        mismatch,
                    });
                }
            }
            for edit in self
                .perception
                .plan_aromatic_system(molecule, atoms, attributes)
            {
                edits.push(edit);
            }
        }
        Ok(edits)
    }
}

impl Transformer for Aromatizer {
    type Error = AromatizeError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, AromatizeError> {
        let edits = self.plan_transform(&molecule)?;
        let editor = molecule
            .edit()
            .apply(edits)
            .expect("aromatization plan applies to its input");
        Ok(editor
            .finish()
            .expect("aromatization plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), AromatizeError> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            *molecule = mem::take(molecule)
                .apply(edits)
                .expect("aromatization plan applies and preserves molecule integrity");
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

/// Read localized contributions using the perception's candidate ring set.
/// Ring double bonds supply one electron; an external double bond remains
/// localized and leaves a carbon site with zero aromatic electrons.
/// Donors require an available pair; saturated and cumulated sites are excluded.
pub(crate) fn electrons_from_kekule(
    molecule: &Molecule,
    rings: &RingSet,
    atom: AtomId,
) -> Result<AromaticValence, AromatizeError> {
    let view = molecule.atom(atom);
    let ElementForm::Lit(element) = view.attributes().element else {
        return Err(AromatizeError::Underdetermined);
    };
    let mut double = None;
    for neighbor in view.neighbors() {
        match neighbor
            .bond()
            .order()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            1 => {}
            2 if double.is_none() => double = Some(neighbor.bond_id()),
            _ => return Ok(AromaticValence::NotAromatic),
        }
    }
    if let Some(double) = double {
        return Ok(
            if !rings.contains_atom(atom) || rings.contains_bond(double) {
                AromaticValence::Aromatic(1)
            } else if element == Element::C {
                AromaticValence::Aromatic(0)
            } else {
                AromaticValence::NotAromatic
            },
        );
    }
    let donor = match element {
        Element::C => match view
            .charge()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            1 => return Ok(AromaticValence::Aromatic(0)),
            -1 => true,
            _ => false,
        },
        Element::N | Element::O | Element::S | Element::Se | Element::P | Element::As => true,
        _ => false,
    };
    if donor {
        match view
            .lone_pairs()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            pairs if pairs > 0 => return Ok(AromaticValence::Aromatic(2)),
            _ => {}
        }
    }
    Ok(AromaticValence::NotAromatic)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{BondId, RingConfig};
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;
    use crate::ops::model::AromaticityRule;

    #[rstest]
    fn test_aromatizer_with_config() {
        let molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"]
                    [3 4 "1"] [4 5 "2"] [5 0 "1"]]
        }"#
        );
        let expected = Aromatizer::new(&AromaticityModel::daylight())
            .transform(molecule.clone())
            .unwrap();
        let configured = Aromatizer::with_config(
            &AromaticityModel::daylight(),
            AromaticityConfig {
                ring_config: RingConfig {
                    simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                },
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
            },
        )
        .transform(molecule);
        assert_eq!(configured, Ok(expected));
    }

    #[rstest]
    #[case::ring_double_bonds(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    #[case::nitrogen_donor(
        mol_dsl_concrete!(r#"{:atoms ["N#h#n" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    )]
    #[case::oxygen_donor(
        mol_dsl_concrete!(r#"{:atoms ["O#n2" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["O#n" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    )]
    #[case::sulfur_donor(
        mol_dsl_concrete!(r#"{:atoms ["S#n2" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["S#n" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    )]
    #[case::anionic_carbon_donor(
        mol_dsl_concrete!(r#"{:atoms ["C#c-#h#n" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#c-#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    )]
    #[case::external_carbonyl(
        mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#n2"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [5 6 "2"] [6 0 "1"] [0 7 "2"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#n2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"] [0 7 "2"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
    )]
    #[case::empty_carbon_orbital(
        mol_dsl_concrete!(r#"{:atoms ["C#c+#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [5 6 "2"] [6 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#c+#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
    )]
    #[case::aromatic_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a#D2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h#D2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    #[case::excluded_peripheral_bond(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C" "C" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a!"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                    [5 6 "1"] [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 4 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C" "C" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a!"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]
            :aromatic-systems [{:atoms [4 5 6 7 8 9] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    #[case::nonliteral_aromatic_assertion(
        {
            let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
                :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
            molecule.atom_mut(AtomId(0)).attributes_mut().constraints.set(
                AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::lit_set([0, 1])))
            );
            molecule
        },
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    #[case::unrelated_undetermined_fields(
        {
            let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
                :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
            molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens = NumForm::Undetermined;
            molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
            molecule
        },
        {
            let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
            molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens = NumForm::Undetermined;
            molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
            molecule
        },
    )]
    fn test_aromatizer_transform(#[case] molecule: Molecule, #[case] expected: Molecule) {
        let transformer = Aromatizer::new(&AromaticityModel::daylight());
        let mut in_place = molecule.clone();
        assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
        assert_eq!(in_place, expected);
        assert_eq!(
            transformer.transform_iter(&molecule).collect::<Vec<_>>(),
            vec![expected.clone()]
        );
        assert_eq!(transformer.transform(molecule), Ok(expected));
    }

    #[rstest]
    #[case::acyclic(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), AromaticityModel::daylight())]
    #[case::acyclic_donor(
        mol_dsl_concrete!(r#"{:atoms ["N#h3#n"]}"#),
        AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
    )]
    #[case::excluded_atom(mol_dsl_concrete!(r#"{
        :atoms ["C#h#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
    }"#), AromaticityModel::daylight())]
    #[case::excluded_undetermined_atom(mol_dsl_concrete!(r#"{
        :atoms ["*#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
    }"#), AromaticityModel::daylight())]
    #[case::saturated_ring(mol_dsl_concrete!(r#"{
        :atoms ["C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
    }"#), AromaticityModel::daylight())]
    #[case::unavailable_donor_pair(mol_dsl_concrete!(r#"{
        :atoms ["N#c+#h2" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]
    }"#), AromaticityModel::daylight())]
    #[case::cumulated_bonds(mol_dsl_concrete!(r#"{
        :atoms ["C" "C#h" "C#h2" "C#h2" "C#h2" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "2"]]
    }"#), AromaticityModel::daylight())]
    #[case::already_aromatic(mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#), AromaticityModel::daylight())]
    fn test_aromatizer_transform_identity(
        #[case] molecule: Molecule,
        #[case] model: AromaticityModel,
    ) {
        let transformer = Aromatizer::new(&model);
        let mut in_place = molecule.clone();
        assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
        assert_eq!(in_place, molecule);
        assert_eq!(transformer.transform(molecule.clone()), Ok(molecule));
    }

    #[rstest]
    #[case::aromatic_contribution(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
            atom: AtomId(0), system: AromaticSystemId(0),
        }.into()),
    )]
    #[case::later_aromatic_contribution(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                    "C#h" "C#h" "C#h#a2" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                    [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 10 "1"] [10 11 "2"] [11 6 "1"]]
        }"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
            atom: AtomId(8), system: AromaticSystemId(1),
        }.into()),
    )]
    #[case::nonliteral_aromatic_contribution(
        {
            let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
                :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
            molecule.atom_mut(AtomId(0)).attributes_mut().constraints.set(
                AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::lit_set([0, 2])))
            );
            molecule
        },
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
            atom: AtomId(0), system: AromaticSystemId(0),
        }.into()),
    )]
    #[case::unclaimed_ring_atom(
        mol_dsl_concrete!(r#"{:atoms ["C#h2#a" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceFailure {
            atom: AtomId(0),
        }.into()),
    )]
    #[case::unclaimed_ring_bond(
        mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticBondConstraintFailure {
            bond: BondId(0),
        }.into()),
    )]
    #[case::unclaimed_acyclic_atom(
        mol_dsl_concrete!(r#"{:atoms ["C#h4#a"]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceFailure {
            atom: AtomId(0),
        }.into()),
    )]
    #[case::unclaimed_acyclic_bond(
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1#a"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticBondConstraintFailure {
            bond: BondId(0),
        }.into()),
    )]
    #[case::conflicting_bond_evidence(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceBondMismatch {
            atom: AtomId(0), bond: BondId(0),
        }.into()),
    )]
    #[case::clar_heterocycle(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
        }"#),
        AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
        AromatizeError::Contradiction(AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".into()
        ))
    )]
    #[case::undetermined_donor_pair(
        {
            let mut molecule = mol_dsl_concrete!(r#"{:atoms ["N#h#n" "C#h" "C#h" "C#h" "C#h"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#);
            molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
            molecule
        },
        AromaticityModel::daylight(),
        AromatizeError::Underdetermined,
    )]
    #[case::undetermined_bond_order(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "*"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::Underdetermined,
    )]
    #[case::contradicted_valence(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h#v3" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        AromaticityModel::daylight(),
        AromatizeError::ContradictedAssertion {
            aromatic_system: AromaticSystemId(0),
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(3),
        },
    )]
    #[case::later_aromatic_system_mismatch(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                    "C#h" "C#h" "C#h#n" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                    [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 10 "1"] [10 11 "2"] [11 6 "1"]]
        }"#),
        AromaticityModel::daylight(),
        AromatizeError::ValenceInvariant {
            aromatic_system: AromaticSystemId(1),
            mismatch: ValenceMismatch::OrbitalCount {
                atom_id: AtomId(8), orbital_count: 10, electron_count: 8,
            },
        },
    )]
    fn test_aromatizer_transform_error(
        #[case] molecule: Molecule,
        #[case] model: AromaticityModel,
        #[case] expected: AromatizeError,
    ) {
        let transformer = Aromatizer::new(&model);
        let mut in_place = molecule.clone();
        assert_eq!(
            transformer.transform_into(&mut in_place),
            Err(expected.clone())
        );
        assert_eq!(in_place, molecule);
        assert_eq!(transformer.transform_iter(&molecule).next(), None);
        assert_eq!(transformer.transform(molecule), Err(expected));
    }
}
