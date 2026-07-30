//! Delocalize atom-localized charge over homogeneous aromatic systems.

use std::convert::Infallible;
use std::iter::once;

use umol_ast::ast::{
    AromaticSystemId, AromaticValenceAst, AtomConstraintAst, AtomId, ElectronCountsAst, ElementAst,
    MoleculeAst, ValueAst,
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
    fn derive(ast: &MoleculeAst, system: AromaticSystemId) -> Option<Self> {
        let view = ast.aromatic_system(system);
        let atom_ids: Vec<AtomId> = view.atom_ids().collect();
        let (&first, rest) = atom_ids.split_first()?;
        let ElementAst::Lit(element) = ast.atom(first).ast.element else {
            return None;
        };
        if rest
            .iter()
            .any(|&atom| ast.atom(atom).ast.element != ElementAst::Lit(element))
        {
            return None;
        }

        let ElectronCountsAst::Lit(old_electrons) = &view.ast.electrons else {
            return None;
        };
        if old_electrons.len() != atom_ids.len() {
            return None;
        }
        let ValueAst::Lit(mut charge) = view.ast.charge else {
            return None;
        };

        let mut electrons = Vec::with_capacity(atom_ids.len());
        let mut atoms = Vec::with_capacity(atom_ids.len());
        for atom_id in atom_ids {
            let atom = ast.atom(atom_id);
            let ValueAst::Lit(atom_charge) = atom.charge() else {
                return None;
            };
            let ValueAst::Lit(degree) = atom.degree() else {
                return None;
            };
            let ValueAst::Lit(implicit_hydrogens) = atom.implicit_hydrogens() else {
                return None;
            };
            let ValueAst::Lit(lone_pairs) = atom.lone_pairs() else {
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

    fn apply(self, ast: &mut MoleculeAst) {
        for (atom_id, contribution) in self.atoms {
            let atom = &mut ast.atom_mut(atom_id).ast;
            atom.charge = ValueAst::Lit(0);
            atom.constraints.set(AtomConstraintAst::AromaticValence(
                AromaticValenceAst::Aromatic(ValueAst::Lit(contribution)),
            ));
        }
        let system = &mut ast.aromatic_system_mut(self.system).ast;
        system.charge = ValueAst::Lit(self.charge);
        system.electrons = ElectronCountsAst::Lit(self.electrons);
    }
}

impl Transformer for DelocalizeCharge {
    type Error = Infallible;

    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), Self::Error> {
        let plans: Vec<DelocalizationPlan> = ast
            .aromatic_systems()
            .ids()
            .filter_map(|system| DelocalizationPlan::derive(ast, system))
            .collect();
        for plan in plans {
            plan.apply(ast);
        }
        Ok(())
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a> {
        let transformed = match self.transform(ast) {
            Ok(transformed) => transformed,
            Err(never) => match never {},
        };
        Box::new(once(transformed))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_ast::{mol_dsl, mol_dsl_ground};

    use super::*;

    #[rstest]
    #[case::cyclopentadienyl(
        mol_dsl_ground!(r#"{
            :atoms ["C#c-#h#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[2,1,1,1,1]"}]
        }"#),
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]#c-"}]
        }"#)
    )]
    #[case::tropylium(
        mol_dsl_ground!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                    [6 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                                :type "[0,1,1,1,1,1,1]"}]
        }"#),
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                    [6 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                                :type "[1,1,1,1,1,1,1]#c+"}]
        }"#)
    )]
    #[case::multiple_systems(
        mol_dsl_ground!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"
                    "C#c+#h#a0" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
            :aromatic-systems [
                {:atoms [0 1 2] :type "[0,1,1]"}
                {:atoms [3 4 5] :type "[0,1,1]"}]
        }"#),
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a"
                    "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
            :aromatic-systems [
                {:atoms [0 1 2] :type "[1,1,1]#c+"}
                {:atoms [3 4 5] :type "[1,1,1]#c+"}]
        }"#)
    )]
    fn test_delocalize_charge_transform(#[case] input: MoleculeAst, #[case] expected: MoleculeAst) {
        assert_eq!(DelocalizeCharge.transform(&input), Ok(expected));
    }

    #[rstest]
    #[case::already_delocalized(mol_dsl_ground!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]#c-"}]
    }"#))]
    #[case::heterogeneous(mol_dsl_ground!(r#"{
        :atoms ["B#c-#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
    }"#))]
    #[case::non_literal(mol_dsl!(r#"{
        :atoms ["C" "C" "C"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2] :type "*"}]
    }"#))]
    fn test_delocalize_charge_transform_identity(#[case] input: MoleculeAst) {
        assert_eq!(DelocalizeCharge.transform(&input), Ok(input));
    }

    #[rstest]
    #[case(
        mol_dsl_ground!(r#"{
            :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2] :type "[0,1,1]"}]
        }"#),
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]#c+"}]
        }"#)
    )]
    fn test_delocalize_charge_generate_all(
        #[case] input: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        assert_eq!(
            DelocalizeCharge.generate_all(&input).collect::<Vec<_>>(),
            vec![expected]
        );
    }
}
