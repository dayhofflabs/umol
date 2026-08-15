//! Delocalize atom-localized charge over homogeneous aromatic systems.

use std::convert::Infallible;
use std::iter;

use umol_graph_ir::ir::{
    AromaticSystemId, AromaticValenceForm, AtomConstraintForm, AtomId, ElectronCountsForm,
    ElementForm, Molecule, NumForm,
};

use crate::ops::transform::Transformer;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DelocalizeCharge;

#[derive(Clone, Debug, PartialEq, Eq)]
struct DelocalizationPlan {
    system: AromaticSystemId,
    charge: i64,
    electrons: Vec<i64>,
    atoms: Vec<(AtomId, i64)>,
}

impl DelocalizationPlan {
    fn derive(molecule: &Molecule, system: AromaticSystemId) -> Option<Self> {
        let view = molecule.aromatic_system(system);
        let atom_ids: Vec<AtomId> = view.atom_ids().collect();
        let (&first, rest) = atom_ids.split_first()?;
        let ElementForm::Lit(element) = molecule.atom(first).attributes.element else {
            return None;
        };
        if rest
            .iter()
            .any(|&atom| molecule.atom(atom).attributes.element != ElementForm::Lit(element))
        {
            return None;
        }

        let ElectronCountsForm::Lit(old_electrons) = &view.attributes.electrons else {
            return None;
        };
        if old_electrons.len() != atom_ids.len() {
            return None;
        }
        let NumForm::Lit(mut charge) = view.attributes.charge else {
            return None;
        };

        let mut electrons = Vec::with_capacity(atom_ids.len());
        let mut atoms = Vec::with_capacity(atom_ids.len());
        for atom_id in atom_ids {
            let atom = molecule.atom(atom_id);
            let NumForm::Lit(atom_charge) = atom.charge() else {
                return None;
            };
            let NumForm::Lit(degree) = atom.degree() else {
                return None;
            };
            let NumForm::Lit(implicit_hydrogens) = atom.implicit_hydrogens() else {
                return None;
            };
            let NumForm::Lit(lone_pairs) = atom.lone_pairs() else {
                return None;
            };
            let contribution = i64::from(element.valence_electrons())
                - degree
                - implicit_hydrogens
                - 2 * lone_pairs;
            charge += atom_charge;
            electrons.push(contribution);
            atoms.push((atom_id, contribution));
        }

        Some(Self {
            system,
            charge,
            electrons,
            atoms,
        })
    }

    fn apply(self, molecule: &mut Molecule) {
        for (atom_id, contribution) in self.atoms {
            let atom = &mut molecule.atom_mut(atom_id).attributes;
            atom.charge = NumForm::Lit(0);
            atom.constraints.set(AtomConstraintForm::AromaticValence(
                AromaticValenceForm::Aromatic(NumForm::Lit(contribution)),
            ));
        }
        let system = &mut molecule.aromatic_system_mut(self.system).attributes;
        system.charge = NumForm::Lit(self.charge);
        system.electrons = ElectronCountsForm::Lit(self.electrons);
    }
}

impl Transformer for DelocalizeCharge {
    type Error = Infallible;

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error> {
        let plans: Vec<DelocalizationPlan> = molecule
            .aromatic_systems()
            .ids()
            .filter_map(|system| DelocalizationPlan::derive(molecule, system))
            .collect();
        for plan in plans {
            plan.apply(molecule);
        }
        Ok(())
    }

    fn generate_all<'a>(
        &'a self,
        molecule: &'a Molecule,
    ) -> Box<dyn Iterator<Item = Molecule> + 'a> {
        let transformed = match self.transform(molecule) {
            Ok(transformed) => transformed,
            Err(never) => match never {},
        };
        Box::new(iter::once(transformed))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_graph_ir::{mol_dsl, mol_dsl_concrete};

    use super::*;

    #[rstest]
    #[case::cyclopentadienyl(
        mol_dsl_concrete!(r#"{
            :atoms ["C#c-#h#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]
        }"#)
    )]
    #[case::tropylium(
        mol_dsl_concrete!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                    [6 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                                :attrs "[0,1,1,1,1,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                    [6 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                                :attrs "[1,1,1,1,1,1,1]#c+"}]
        }"#)
    )]
    #[case::multiple_systems(
        mol_dsl_concrete!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"
                    "C#c+#h#a0" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
            :aromatic-systems [
                {:atoms [0 1 2] :attrs "[0,1,1]"}
                {:atoms [3 4 5] :attrs "[0,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
            :aromatic-systems [
                {:atoms [0 1 2] :attrs "[1,1,1]#c+"}
                {:atoms [3 4 5] :attrs "[1,1,1]#c+"}]
        }"#)
    )]
    fn test_delocalize_charge_transform(#[case] input: Molecule, #[case] expected: Molecule) {
        assert_eq!(DelocalizeCharge.transform(&input), Ok(expected));
    }

    #[rstest]
    #[case::already_delocalized(mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]
    }"#))]
    #[case::heterogeneous(mol_dsl_concrete!(r#"{
        :atoms ["B#c-#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#))]
    #[case::non_literal(mol_dsl!(r#"{
        :atoms ["C" "C" "C"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "*"}]
    }"#))]
    fn test_delocalize_charge_transform_identity(#[case] input: Molecule) {
        assert_eq!(DelocalizeCharge.transform(&input), Ok(input));
    }

    #[rstest]
    #[case(
        mol_dsl_concrete!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c+"}]
        }"#)
    )]
    fn test_delocalize_charge_generate_all(#[case] input: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            DelocalizeCharge.generate_all(&input).collect::<Vec<_>>(),
            vec![expected]
        );
    }
}
