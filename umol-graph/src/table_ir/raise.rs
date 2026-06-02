//! TableIR → `umol_ast::MoleculeAst` raise.
//!
//! Implements `TryIntoAst<MoleculeAst> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! raise applies fixed IO ground semantics for resolution.

use std::collections::HashSet;

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomId, BondAst, BondConstraint, Constraints,
    DativeBondAst, ElementAst, IsotopeMassAst, MoleculeAst, MulticenterBondAst, NoncovalentBondAst,
    NoncovalentBondKind, SpinStateAst, TryIntoAst, ValueAst,
};

use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondNoncovalent as TableNoncovalent,
    BondOrder as TableBondOrder,
};
use crate::table_ir::Molecule as TableMolecule;

/// Failure modes for the TableIR → MoleculeAst raise.
///
/// Empty for now: every TableIR field maps mechanically to a MoleculeAst
/// equivalent. The variant slot exists so that future strictly-checked raise
/// paths (e.g., asserting bond endpoint validity, rejecting unsupported
/// `BondOrder::Any`) can land without changing the public signature.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {}

impl TryIntoAst<MoleculeAst> for &TableMolecule {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, ctx: &Self::Ctx) -> Result<MoleculeAst, RaiseError> {
        let mut atoms: Vec<AtomAst> = self
            .atoms
            .iter()
            .map(|a| a.try_into_ast(ctx))
            .collect::<Result<_, _>>()?;

        let mut bonds = Vec::new();
        let mut dative_bonds: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        for b in &self.bonds {
            let a_idx = AtomId(b.atoms.first());
            let b_idx = AtomId(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent_bonds.push((a_idx, b_idx, NoncovalentBondAst::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        bonds.push((a_idx, b_idx, b.try_into_ast(ctx)?));
                        continue;
                    }
                };
                let dative_bond = DativeBondAst::new(raise_bond_order(b.order));
                dative_bonds.push((vec![donor], acceptor, dative_bond));
            } else {
                bonds.push((a_idx, b_idx, b.try_into_ast(ctx)?));
            }
        }

        let multicenter_bond: Vec<(Vec<AtomId>, MulticenterBondAst)> = self
            .multicenter_bonds
            .iter()
            .map(|mc| {
                let mut seen = HashSet::new();
                let atoms: Vec<AtomId> = mc
                    .all_atoms()
                    .into_iter()
                    .filter(|a| seen.insert(*a))
                    .map(AtomId)
                    .collect();
                let n = atoms.len();
                (
                    atoms,
                    MulticenterBondAst::new(vec![ValueAst::Undetermined; n]),
                )
            })
            .collect();

        for (table_atom, atom) in self.atoms.iter().zip(atoms.iter_mut()) {
            match table_atom.aromatic {
                Some(true) => {
                    atom.constraints.add(AtomConstraint::AromaticValence(
                        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
                    ));
                }
                Some(false) => {
                    atom.constraints.add(AtomConstraint::AromaticValence(
                        AromaticValenceAst::NotAromatic,
                    ));
                }
                None => {}
            }

            if matches!(atom.isotope_mass, IsotopeMassAst::Undetermined) {
                atom.isotope_mass = IsotopeMassAst::Natural;
            }
            if matches!(atom.charge, ValueAst::Undetermined) {
                atom.charge = ValueAst::Lit(0);
            }
            if matches!(atom.spin.unpaired, ValueAst::Undetermined) {
                atom.spin.unpaired = ValueAst::Lit(0);
            }
            atom.constraints.retain(|c| !c.is_undetermined());
        }
        let constraints = Constraints::new();

        Ok(MoleculeAst::from_parts(
            atoms,
            bonds,
            dative_bonds,
            vec![],
            multicenter_bond,
            noncovalent_bonds,
            constraints,
        ))
    }
}

impl TryIntoAst<AtomAst> for &TableAtom {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, _ctx: &Self::Ctx) -> Result<AtomAst, RaiseError> {
        Ok(AtomAst {
            element: ElementAst::Lit(self.element),
            isotope_mass: match self.isotope_mass {
                Some(m) => IsotopeMassAst::Lit(m as i64),
                None => IsotopeMassAst::Undetermined,
            },
            charge: match self.charge {
                Some(c) => ValueAst::Lit(c as i64),
                None => ValueAst::Undetermined,
            },
            implicit_hydrogens: match self.implicit_hydrogens {
                Some(n) => ValueAst::Lit(n as i64),
                None => ValueAst::Undetermined,
            },
            lone_pairs: match self.lone_pairs {
                Some(n) => ValueAst::Lit(n as i64),
                None => ValueAst::Undetermined,
            },
            spin: match self.unpaired_electrons {
                Some(u) => SpinStateAst {
                    unpaired: ValueAst::Lit(u as i64),
                    multiplicity: ValueAst::Undetermined,
                },
                None => SpinStateAst::default(),
            },
            constraints: Default::default(),
        })
    }
}

impl TryIntoAst<BondAst> for &TableBond {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, _ctx: &Self::Ctx) -> Result<BondAst, RaiseError> {
        let mut bond = BondAst::new(raise_bond_order(self.order));
        bond.charge = match self.charge {
            Some(c) => ValueAst::Lit(c as i64),
            None => ValueAst::Undetermined,
        };
        if matches!(self.order, TableBondOrder::Aromatic) {
            bond.constraints.add(BondConstraint::Aromatic);
        }
        Ok(bond)
    }
}

fn raise_bond_order(order: TableBondOrder) -> ValueAst {
    match order {
        TableBondOrder::Zero => ValueAst::Lit(0),
        TableBondOrder::Single => ValueAst::Lit(1),
        TableBondOrder::Double => ValueAst::Lit(2),
        TableBondOrder::Triple => ValueAst::Lit(3),
        TableBondOrder::Quadruple => ValueAst::Lit(4),
        TableBondOrder::Quintuple => ValueAst::Lit(5),
        TableBondOrder::Sextuple => ValueAst::Lit(6),
        // Definite-aromatic: localized bond order is 1 by Kekulé convention;
        // the aromatic flag is added separately as `BondConstraint::Aromatic`.
        // Renders as `1#a`.
        TableBondOrder::Aromatic => ValueAst::Lit(1),
        // Fuzzy orders: no concrete bond order can be assigned; raise to
        // `Undetermined`. Aromatic-flag setting (where applicable) is left
        // off — the chemistry of these is too ambiguous for the raise.
        TableBondOrder::SingleOrDouble
        | TableBondOrder::SingleOrAromatic
        | TableBondOrder::DoubleOrAromatic
        | TableBondOrder::Any => ValueAst::Undetermined,
    }
}

fn noncovalent_kind(kind: TableNoncovalent) -> NoncovalentBondKind {
    match kind {
        TableNoncovalent::Hydrogen => NoncovalentBondKind::HydrogenBond,
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::ast::BondId;
    use umol_shared::element::Element;

    use super::*;
    use crate::io::ctfile::config::CtfileIoConfig;
    use crate::io::ctfile::{parse_mol_bytes_with, parse_mol_to_ast};
    use crate::io::smiles::parse_smiles_to_ast;
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, CountsModel, ElementScope, RingLimits, ValenceModel,
    };
    use crate::ops::valence::{CountsValence, ValenceTable};
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder as TableBondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    fn methane_table() -> TableMolecule {
        let mut atom = TableAtom::from_element(Element::C);
        atom.implicit_hydrogens = Some(4);
        let mut mol = TableMolecule::empty();
        mol.atoms.push(atom);
        mol
    }

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    const CARBON_H0_EXPLICIT_MOL: &str = "carbon-h0\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  1  0  0  0  0  0  0  0  0\nM  END\n";

    fn counts_model() -> CountsModel {
        CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        }
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_methane_table() {
        let mol = methane_table();
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        assert_eq!(ast.atoms().count(), 1);
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.element, ElementAst::Lit(Element::C));
        assert!(matches!(atom.implicit_hydrogens, ValueAst::Lit(4)));
        assert!(matches!(atom.charge, ValueAst::Lit(0)));
        assert!(matches!(atom.lone_pairs, ValueAst::Undetermined));
        assert!(matches!(atom.spin.unpaired, ValueAst::Lit(0)));
    }

    #[rstest]
    #[case::table_aromatic_none(None, AromaticValenceAst::Undetermined)]
    #[case::table_aromatic_false(Some(false), AromaticValenceAst::NotAromatic)]
    #[case::table_aromatic_true(Some(true), AromaticValenceAst::Aromatic(ValueAst::Undetermined))]
    fn test_table_molecule_try_into_ast(
        #[case] table_aromatic: Option<bool>,
        #[case] expected: AromaticValenceAst,
    ) {
        let mut mol = TableMolecule::empty();
        let mut atom = TableAtom::from_element(Element::C);
        atom.aromatic = table_aromatic;
        mol.atoms.push(atom);
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.constraints.aromatic_valence(),
            expected
        );
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_bond_order() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds.push(TableBond::new(0, 1, TableBondOrder::Double));
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert!(matches!(bond.order, ValueAst::Lit(2)));
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_aromatic_bond() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds
            .push(TableBond::new(0, 1, TableBondOrder::Aromatic));
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert!(matches!(bond.order, ValueAst::Lit(1)));
        assert!(bond
            .constraints
            .iter()
            .any(|c| matches!(c, BondConstraint::Aromatic)));
        for i in 0..2 {
            assert!(matches!(
                ast.atom(AtomId(i)).ast.constraints.aromatic_valence(),
                AromaticValenceAst::Undetermined
            ));
        }
    }

    #[rstest]
    #[case::methane(METHANE_MOL, "C#i=#c0#u0")]
    #[case::benzene(BENZENE_AROMATIC_MOL, "C#i=#c0#u0")]
    #[case::carbon_h0(CARBON_H0_EXPLICIT_MOL, "C#i=#c0#h0#u0")]
    fn test_parse_mol_to_ast(#[case] input: &str, #[case] expected_atom: &str) {
        let ast = parse_mol_to_ast(input).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.charge, ValueAst::Lit(0));
        assert!(matches!(
            atom.constraints.aromatic_valence(),
            AromaticValenceAst::Undetermined
        ));
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    #[case::methane(METHANE_MOL, 1, "C#i=#c0#h4#n0#u0#s#v0#a!")]
    #[case::benzene(BENZENE_AROMATIC_MOL, 6, "C#i=#c0#h#n0#u0#s#v2#a")]
    fn test_parse_mol_to_ast_counts_resolve(
        #[case] input: &str,
        #[case] atom_count: u32,
        #[case] expected_atom: &str,
    ) {
        let mut ast = parse_mol_to_ast(input).unwrap();
        CountsValence::new(&counts_model())
            .resolve(&mut ast)
            .unwrap();
        assert_eq!(ast.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(ast.atom(AtomId(i)).ast.to_string(), expected_atom);
        }
    }

    #[rstest]
    fn test_parse_smiles_to_ast_methane_counts_resolve() {
        let mut ast = parse_smiles_to_ast("C").unwrap();
        CountsValence::new(&counts_model())
            .resolve(&mut ast)
            .unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
    }

    #[rstest]
    #[case::organic("C", "C#i=#c0#u0#a!")]
    fn test_parse_smiles_to_ast(#[case] input: &str, #[case] expected_atom: &str) {
        let ast = parse_smiles_to_ast(input).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.charge, ValueAst::Lit(0));
        assert!(matches!(atom.implicit_hydrogens, ValueAst::Undetermined));
        assert!(matches!(
            atom.constraints.aromatic_valence(),
            AromaticValenceAst::NotAromatic
        ));
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    fn test_parse_mol_bytes_with_resolver_methane_determined() {
        let model = ChemistryModel {
            valence: ValenceModel::Counts(CountsModel {
                table: Cow::Borrowed(ValenceTable::default_table()),
            }),
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
        };
        let ast =
            parse_mol_bytes_with(METHANE_MOL.as_bytes(), &CtfileIoConfig::basic(), &model).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
    }
}
