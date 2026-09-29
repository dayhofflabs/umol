//! Aromatize: Kekulé form → aromatic-system form.
//!
//! [`Aromatizer`] runs aromaticity perception against a Kekulé-form input —
//! atoms with explicit single/double bonds and no aromatic hints. Per-atom π
//! contributions are derived from bond orders by [`electrons_from_kekule`]
//! rather than from the `AromaticValence` constraint that the resolver reads.
//! If the input molecule already carries one or more aromatic systems, this is a
//! no-op: re-aromatizing requires kekulizing first.

use std::iter;

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomHandle, AtomId, BondConstraintForm, BondHandle, BondUpdate, Edits, ElementForm, Molecule,
    NumForm,
};

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityError, AromaticityPerceiver,
};
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
}

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
        let systems = self
            .perception
            .find_systems(molecule, self.config, |atom| {
                electrons_from_kekule(molecule, atom)
            })?
            .into_decisive(AromatizeError::Underdetermined)?;
        for (atoms, attributes) in systems {
            let bond_ids = molecule.bonds().induced_ids(&atoms);
            edits.add_aromatic_system(atoms.into_iter().map(AtomHandle::Id).collect(), attributes);
            for id in bond_ids {
                edits.update_bond(
                    BondHandle::Id(id),
                    molecule.bond(id).attributes(),
                    &BondUpdate {
                        constraints: BondConstraintForm::aromatic(true).into(),
                        ..Default::default()
                    },
                );
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
            molecule
                .transact([edits])
                .expect("aromatization plan preserves molecule integrity");
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

/// Derive an atom's π contribution from a Kekulé bond-order layout.
///
/// - Exactly one incident double bond → 1 π electron (sp² atom on a single
///   π bond, e.g. benzene C, pyridine N).
/// - Zero incident double bonds, atom is an N/O/S/Se/P/As → 2 π electrons
///   (pyrrole-, furan-, thiophene-class heteroatom donating a lone pair).
/// - Zero incident double bonds, atom is C with charge `+1` → 0 π electrons
///   (sp² carbocation, empty p_z, e.g. tropylium C⁺).
/// - Anything else (sp³ C, two or more double bonds, undetermined data) →
///   `None`, marking the atom as not aromatic-eligible.
pub fn electrons_from_kekule(molecule: &Molecule, atom: AtomId) -> Option<u8> {
    let view = molecule.atom(atom);
    let ElementForm::Lit(element) = view.attributes().element else {
        return None;
    };
    let double_count = view
        .neighbors()
        .filter(|n| matches!(n.bond().attributes().order, NumForm::Lit(2)))
        .count();
    match double_count {
        1 => Some(1),
        0 => match element {
            Element::N | Element::O | Element::S | Element::Se | Element::P | Element::As => {
                Some(2)
            }
            Element::C if matches!(view.attributes().charge, NumForm::Lit(1)) => Some(0),
            _ => None,
        },
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AtomForm, AtomId, BondForm, Molecule, MoleculeEntries, RingConfig, UnpairedElectronsForm,
    };
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;
    use crate::ops::model::AromaticityRule;

    fn kekule_carbon() -> AtomForm {
        let mut atom = AtomForm::from_element(Element::C);
        atom.charge = NumForm::Lit(0);
        atom.unpaired_electrons = UnpairedElectronsForm::closed_shell();
        atom
    }

    fn benzene_kekule() -> Molecule {
        let atoms: Vec<AtomForm> = (0..6).map(|_| kekule_carbon()).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| {
                let order = if i % 2 == 0 { 2 } else { 1 };
                (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(order))
            })
            .collect();
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    #[rstest]
    fn test_aromatizer_with_config() {
        let molecule = mol_dsl_concrete!(
            r#"{
                :atoms ["C" "C" "C" "C" "C" "C"]
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
    #[case::benzene(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a"] [1 2 "1#a"] [2 3 "2#a"] [3 4 "1#a"] [4 5 "2#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    fn test_aromatizer_transform(#[case] molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            Aromatizer::new(&AromaticityModel::daylight()).transform(molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::acyclic(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#))]
    #[case::already_aromatic(mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#))]
    fn test_aromatizer_transform_identity(#[case] molecule: Molecule) {
        assert_eq!(
            Aromatizer::new(&AromaticityModel::daylight()).transform(molecule.clone()),
            Ok(molecule)
        );
    }

    #[rstest]
    #[case::clar_heterocycle(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h0#n" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 :double] [1 2 :single] [2 3 :double]
                    [3 4 :single] [4 5 :double] [5 0 :single]]
        }"#),
        AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
        AromatizeError::Contradiction(AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".into()
        ))
    )]
    fn test_aromatizer_transform_error(
        #[case] molecule: Molecule,
        #[case] model: AromaticityModel,
        #[case] expected: AromatizeError,
    ) {
        let transformer = Aromatizer::new(&model);
        let mut borrowed = molecule.clone();
        assert_eq!(
            transformer.transform_into(&mut borrowed),
            Err(expected.clone())
        );
        assert_eq!(borrowed, molecule);
        assert_eq!(transformer.transform_iter(&molecule).next(), None);
        assert_eq!(transformer.transform(molecule), Err(expected));
    }

    #[rstest]
    #[case::benzene(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a"] [1 2 "1#a"] [2 3 "2#a"] [3 4 "1#a"] [4 5 "2#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    fn test_aromatizer_transform_into(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        Aromatizer::new(&AromaticityModel::daylight())
            .transform_into(&mut molecule)
            .unwrap();
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_aromatizer_already_aromatic_is_noop() {
        let original = {
            let mut molecule = benzene_kekule();
            Aromatizer::new(&AromaticityModel::daylight())
                .transform_into(&mut molecule)
                .unwrap();
            molecule
        };
        let mut second = original.clone();
        Aromatizer::new(&AromaticityModel::daylight())
            .transform_into(&mut second)
            .unwrap();
        assert_eq!(original, second);
    }

    #[rstest]
    #[case::benzene(
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2#a"] [1 2 "1#a"] [2 3 "2#a"] [3 4 "1#a"] [4 5 "2#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    )]
    fn test_aromatizer_transform_iter(#[case] molecule: Molecule, #[case] expected: Molecule) {
        let original = molecule.clone();
        let transformer = Aromatizer::new(&AromaticityModel::daylight());
        let mut results = transformer.transform_iter(&molecule);
        assert_eq!(results.next(), Some(expected));
        assert_eq!(results.next(), None);
        assert_eq!(results.next(), None);
        assert_eq!(molecule, original);
    }
}
