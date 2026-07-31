//! The perceived-molecule boundary type. Bond perception over a geometric
//! `Molecule` yields per-atom elements, perceived bond orders, and the geometry's
//! total charge and multiplicity; the type lifts into a `MoleculeAst` via [`IntoAst`].
//! Resolution of the lifted AST (hydrogens, per-atom valence, aromaticity) is the
//! caller's job.

use umol_ast::ast::{
    AtomAst, AtomId, BondAst, Constraint, Constraints, ElementAst, IntoAst, MoleculeAst,
    MoleculeConstraint, MoleculeParts, UnpairedElectronsAst, ValueAst,
};
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_geometric::molecule::Molecule;

use crate::bond_perception::{perceive_bonds, BondPerceptionConfig};

/// A molecule perceived from 3D geometry: the per-atom elements, the perceived
/// bonds, and the geometry's total charge and multiplicity. This is the boundary between
/// the geometric model and the AST — it lifts into a `MoleculeAst` via [`IntoAst`].
#[derive(Clone, Debug, PartialEq)]
pub struct PerceivedMolecule {
    /// Element per atom, in geometric-atom order.
    pub elements: Vec<Element>,
    /// Total molecular charge carried over from the geometry.
    pub charge: i32,
    /// Total spin multiplicity carried over from the geometry.
    pub multiplicity: SpinMultiplicity,
    /// Perceived bonds as `(atom_i, atom_j, order)`.
    pub bonds: Vec<(usize, usize, u8)>,
    /// Whether bond perception satisfied every valence constraint.
    pub feasible: bool,
    /// Residual valence violation per atom (actual − target).
    pub valence_residuals: Vec<i32>,
}

impl PerceivedMolecule {
    /// Perceive bonds over `mol` and bundle the result with the geometry's
    /// elements, charge, and multiplicity into the boundary type.
    pub fn perceive(mol: &Molecule, config: &BondPerceptionConfig) -> Self {
        let result = perceive_bonds(mol, config);
        Self {
            elements: (0..mol.atom_count()).map(|i| mol.element(i)).collect(),
            charge: mol.charge(),
            multiplicity: mol.multiplicity(),
            bonds: result.bonds,
            feasible: result.feasible,
            valence_residuals: result.valence_residuals,
        }
    }
}

impl IntoAst<MoleculeAst> for PerceivedMolecule {
    type Ctx = ();

    /// Each atom carries only its element; hydrogens, per-atom charge, and per-atom
    /// spin are left undetermined for the caller's resolver. The total charge and
    /// spin become molecule-scope `ChargeSum` / `SpinSum` constraints over the whole
    /// molecule.
    fn into_ast(self, _ctx: &Self::Ctx) -> MoleculeAst {
        let atoms: Vec<AtomAst> = self
            .elements
            .iter()
            .map(|&element| AtomAst::new(ElementAst::Lit(element)))
            .collect();
        let bonds: Vec<(AtomId, AtomId, BondAst)> = self
            .bonds
            .iter()
            .map(|&(i, j, order)| {
                (
                    AtomId(i as u32),
                    AtomId(j as u32),
                    BondAst::new(ValueAst::Lit(i64::from(order))),
                )
            })
            .collect();
        let multiplicity = u8::from(self.multiplicity);
        let constraints = Constraints::from_iter([
            Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: ValueAst::Lit(i64::from(self.charge)),
            }),
            Constraint::Molecule(MoleculeConstraint::SpinSum {
                atoms: None,
                spin: UnpairedElectronsAst::from((multiplicity - 1, multiplicity)),
            }),
        ]);
        MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            constraints,
            ..Default::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_chem::element::Element::{C, H, O};
    use umol_chem::spin::SpinMultiplicity;
    use umol_geometric::molecule::Molecule;

    use super::*;

    #[rstest]
    fn test_perceived_molecule_perceive() {
        // Water: O-H 0.96 Å, H-O-H 104.5°; both O-H bonds single.
        let water = Molecule::from_cartesian_angstrom(
            vec![O, H, H],
            &[
                0.000, 0.000, 0.000, 0.960, 0.000, 0.000, -0.240, 0.930, 0.000,
            ],
            0,
            SpinMultiplicity::SINGLET,
        );
        let perceived = PerceivedMolecule::perceive(&water, &BondPerceptionConfig::default());
        assert_eq!(perceived.elements, vec![O, H, H]);
        assert_eq!(perceived.charge, 0);
        assert_eq!(perceived.multiplicity, SpinMultiplicity::SINGLET);
        assert_eq!(perceived.bonds, vec![(0, 1, 1), (0, 2, 1)]);
        assert!(perceived.feasible);
    }

    #[rstest]
    #[case::neutral_singlet(
        PerceivedMolecule {
            elements: vec![C, C],
            charge: 0,
            multiplicity: SpinMultiplicity::SINGLET,
            bonds: vec![(0, 1, 2)],
            feasible: true,
            valence_residuals: vec![0, 0],
        },
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::new(ElementAst::Lit(C)), AtomAst::new(ElementAst::Lit(C))],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::new(ValueAst::Lit(2)))],
            constraints: Constraints::from_iter([
                Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) }),
                Constraint::Molecule(MoleculeConstraint::SpinSum { atoms: None, spin: UnpairedElectronsAst::from((0u8, 1u8)) }),
            ]),
            ..Default::default()
        })
    )]
    #[case::anion_doublet(
        PerceivedMolecule {
            elements: vec![O],
            charge: -1,
            multiplicity: SpinMultiplicity::DOUBLET,
            bonds: vec![],
            feasible: true,
            valence_residuals: vec![0],
        },
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::new(ElementAst::Lit(O))],
            constraints: Constraints::from_iter([
                Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(-1) }),
                Constraint::Molecule(MoleculeConstraint::SpinSum { atoms: None, spin: UnpairedElectronsAst::from((1u8, 2u8)) }),
            ]),
            ..Default::default()
        })
    )]
    fn test_perceived_molecule_into_ast(
        #[case] perceived: PerceivedMolecule,
        #[case] expected: MoleculeAst,
    ) {
        assert_eq!(perceived.into_ast(&()), expected);
    }
}
