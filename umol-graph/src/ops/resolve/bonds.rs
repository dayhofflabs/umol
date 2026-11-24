//! Covalent-bond resolver: fills `BondForm` charge and unpaired-electron defaults on bonds
//! whose corresponding inherent fields are still `Undetermined` after the
//! valence and aromaticity passes.

use thiserror::Error;
use umol_graph_ir::ir::{
    BondHandle, Edits, Lattice, Molecule, Normalize, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use super::WildcardPolicy;

/// Resolves localized-bond defaults from the molecule's evidence.
///
/// # Semantic properties
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input.
#[derive(Clone, Debug, Default)]
pub struct BondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsError {}

impl BondsResolver {
    pub fn new() -> Self {
        Self
    }

    /// Construct charge and unpaired-electron default edits without mutating `molecule`.
    pub fn plan_resolve(&self, molecule: &Molecule) -> Edits {
        self.plan_resolve_with_policy(molecule, None).0
    }

    /// Plan the edits and, when given an aggregate policy, check the fields they leave.
    pub(super) fn plan_resolve_with_policy(
        &self,
        molecule: &Molecule,
        wildcards: Option<WildcardPolicy>,
    ) -> (Edits, bool) {
        let mut admissible = true;
        let mut charges = Vec::new();
        let mut unpaired_electrons = Vec::new();
        for bond_id in molecule.bonds().ids() {
            let bond = molecule.bond(bond_id).attributes();
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            let default_charge = bond.charge.is_undetermined();
            if default_charge {
                charges.push((
                    BondHandle::Id(bond_id),
                    NumForm::Undetermined,
                    NumForm::Lit(0),
                ));
            }
            if selected_unpaired_electrons.is_undetermined() {
                selected_unpaired_electrons = UnpairedElectronsForm::closed_shell();
            } else {
                selected_unpaired_electrons.high_spin_complete();
            }
            if let Some(policy) = wildcards {
                admissible = admissible
                    && (bond.order.is_ground()
                        || (policy == WildcardPolicy::Admit && bond.is_wildcard()))
                    && (default_charge || bond.charge.is_ground())
                    && selected_unpaired_electrons.is_ground();
            }
            if !bond
                .unpaired_electrons
                .normalized_eq(&selected_unpaired_electrons)
            {
                unpaired_electrons.push((
                    BondHandle::Id(bond_id),
                    bond.unpaired_electrons.clone(),
                    selected_unpaired_electrons,
                ));
            }
        }
        let mut edits = Edits::new();
        if !charges.is_empty() {
            edits.modify_bond_charges(charges);
        }
        if !unpaired_electrons.is_empty() {
            edits.modify_bond_unpaired_electrons(unpaired_electrons);
        }
        if let Some(policy) = wildcards {
            admissible = admissible
                && molecule.dative_bonds().iter().all(|bond| {
                    bond.order().is_ground()
                        || (policy == WildcardPolicy::Admit && bond.attributes().is_wildcard())
                })
                && molecule.noncovalent_bonds().iter().all(|bond| {
                    bond.kind().is_ground()
                        || (policy == WildcardPolicy::Admit && bond.attributes().is_wildcard())
                });
        }
        (edits, admissible)
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [BondsResolver::resolve_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, BondsContradiction, ()>, BondsError> {
        let edits = self.plan_resolve(&molecule);
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("bond resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("bond resolve plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Plans and applies localized-bond defaults without a recovery journal.
    ///
    /// # Semantic properties
    ///
    /// Planning finishes before mutation. Only Determined changes the molecule;
    /// every other outcome leaves it unchanged.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), BondsContradiction>, BondsError> {
        let edits = self.plan_resolve(molecule);
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        BondId, Edit, Edits, MoleculeApplyError, RollbackMethod, TransactionError,
    };
    use umol_graph_ir::{bond_dsl, mol_dsl};

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        Edits::from_iter([
            Edit::ModifyBondCharges {
                changes: vec![(
                    BondHandle::Id(BondId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(0),
                )],
            },
            Edit::ModifyBondUnpairedElectrons {
                changes: vec![(
                    BondHandle::Id(BondId(0)),
                    UnpairedElectronsForm::default(),
                    UnpairedElectronsForm::closed_shell(),
                )],
            },
        ]),
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#s3"]]}"#),
        Edits::from_iter([Edit::ModifyBondUnpairedElectrons {
            changes: vec![(
                BondHandle::Id(BondId(0)),
                UnpairedElectronsForm {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(3),
                },
                UnpairedElectronsForm::from((2_u8, 3_u8)),
            )],
        }]),
    )]
    #[case::independent_fields(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1#c+#s3"] [2 3 "1#u0#s"]]}"#
        ),
        Edits::from_iter([
            Edit::ModifyBondCharges {
                changes: vec![
                    (
                        BondHandle::Id(BondId(0)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                    (
                        BondHandle::Id(BondId(2)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                ],
            },
            Edit::ModifyBondUnpairedElectrons {
                changes: vec![
                    (
                        BondHandle::Id(BondId(0)),
                        UnpairedElectronsForm::default(),
                        UnpairedElectronsForm::closed_shell(),
                    ),
                    (
                        BondHandle::Id(BondId(1)),
                        UnpairedElectronsForm {
                            count: NumForm::Undetermined,
                            multiplicity: NumForm::Lit(3),
                        },
                        UnpairedElectronsForm::from((2_u8, 3_u8)),
                    ),
                ],
            },
        ]),
    )]
    fn test_bonds_resolver_plan_resolve(#[case] molecule: Molecule, #[case] expected: Edits) {
        assert_eq!(BondsResolver::new().plan_resolve(&molecule), expected);
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#u2#s1"]]}"#))]
    fn test_bonds_resolver_plan_resolve_identity(#[case] molecule: Molecule) {
        assert_eq!(BondsResolver::new().plan_resolve(&molecule), Edits::new());
    }

    #[rstest]
    #[case::literal("1", "1#c0#u0#s1", true, true)]
    #[case::wildcard("*", "*#c0#u0#s1", false, true)]
    #[case::partial_spin("*#s3", "*#c0#u2#s3", false, true)]
    #[case::order_set("{1,2}", "{1,2}#c0#u0#s1", false, false)]
    #[case::order_variable("?order", "?order#c0#u0#s1", false, false)]
    #[case::charge_variable("*#c?charge", "*#c?charge#u0#s1", false, false)]
    #[case::spin_variable("1#u?count", "1#c0#u?count", false, false)]
    fn test_bonds_resolver_plan_resolve_with_policy(
        #[case] input: &str,
        #[case] resolved: &str,
        #[case] strict: bool,
        #[case] admit: bool,
    ) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "*"] [1 2 "1"]]}"#);
        let mut expected = molecule.clone();
        *molecule.bond_mut(BondId(0)).attributes_mut() = bond_dsl!(input);
        *expected.bond_mut(BondId(0)).attributes_mut() = bond_dsl!(resolved);
        *expected.bond_mut(BondId(1)).attributes_mut() = bond_dsl!("1#c0#u0#s1");
        let resolver = BondsResolver::new();
        let edits = resolver.plan_resolve(&molecule);
        assert_eq!(molecule.clone().apply(edits.clone()), Ok(expected));
        for (policy, admissible) in [
            (WildcardPolicy::Strict, strict),
            (WildcardPolicy::Admit, admit),
        ] {
            assert_eq!(
                resolver.plan_resolve_with_policy(&molecule, Some(policy)),
                (edits.clone(), admissible)
            );
        }
    }

    #[rstest]
    #[case::empty(Molecule::new(), true, true)]
    #[case::literal(mol_dsl!(r#"{:atoms ["N" "B" "O"]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
        :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#), true, true)]
    #[case::dative_wildcard(mol_dsl!(r#"{:atoms ["N" "B"]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#), false, true)]
    #[case::dative_variable(mol_dsl!(r#"{:atoms ["N" "B"]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "?order"}]}"#), false, false)]
    #[case::dative_set(mol_dsl!(r#"{:atoms ["N" "B"]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "{1,2}"}]}"#), false, false)]
    #[case::noncovalent_wildcard(mol_dsl!(r#"{:atoms ["N" "O"]
        :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]}"#), false, true)]
    fn test_bonds_resolver_plan_resolve_with_policy_retained(
        #[case] molecule: Molecule,
        #[case] strict: bool,
        #[case] admit: bool,
    ) {
        let resolver = BondsResolver::new();
        assert_eq!(resolver.plan_resolve(&molecule), Edits::new());
        for (policy, admissible) in [
            (WildcardPolicy::Strict, strict),
            (WildcardPolicy::Admit, admit),
        ] {
            assert_eq!(
                resolver.plan_resolve_with_policy(&molecule, Some(policy)),
                (Edits::new(), admissible)
            );
        }
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#s3"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c0#u2#s3"]]}"#)
    )]
    fn test_bonds_resolver_resolve(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            BondsResolver::new().resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            BondsResolver::new().resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_bonds_resolver_plan_resolve_stale(
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let edits = BondsResolver::new().plan_resolve(&molecule);
        molecule.bond_mut(BondId(1)).attributes_mut().charge = NumForm::Lit(9);
        let expected = molecule.clone();
        assert_eq!(
            molecule.clone().edit().apply(edits.clone()),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(
            molecule.transact([edits], rollback),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, expected);
    }
}
