//! TableIR → `umol_ast::MoleculeAst` raise.
//!
//! Implements `TryIntoAst<MoleculeAst> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! ground defaults (`GROUND_ATOM_DEFAULTS`) fill neutral fields for resolution.

use std::collections::HashSet;
use std::mem;

use strum::IntoEnumIterator;
use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraintKind, AtomConstraints, AtomId,
    BondAst, BondConstraint, Constraints, DativeBondAst, ElementAst, IsotopeAst, MoleculeAst,
    MulticenterBondAst, MulticenterValenceAst, NoncovalentBondAst, NoncovalentBondKind,
    SpinStateAst, TryIntoAst, ValueAst,
};
use umol_ast::dsl::{
    AromaticValenceDefault, AtomDefaults, IsotopeDefault, MulticenterValenceDefault,
    MultiplicityDefault, NumericDefault, UnpairedElectronsDefault,
};
use umol_shared::spin::SpinState;

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

fn ground_atom_defaults() -> AtomDefaults {
    AtomDefaults {
        isotope: IsotopeDefault::Natural,
        charge: NumericDefault::Zero,
        implicit_hydrogens: NumericDefault::Required,
        lone_pairs: NumericDefault::Required,
        unpaired_electrons: UnpairedElectronsDefault::Zero,
        multiplicity: MultiplicityDefault::Required,
        valence: NumericDefault::Required,
        donated_pairs: NumericDefault::Required,
        accepted_pairs: NumericDefault::Required,
        aromatic_valence: AromaticValenceDefault::Required,
        multicenter_valence: MulticenterValenceDefault::Required,
    }
}

impl TryIntoAst<MoleculeAst> for &TableMolecule {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, ctx: &Self::Ctx) -> Result<MoleculeAst, RaiseError> {
        let mut atoms: Vec<AtomAst> = self
            .atoms
            .iter()
            .map(|a| a.try_into_ast(ctx))
            .collect::<Result<_, _>>()?;

        let mut regular = Vec::new();
        let mut dative: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = Vec::new();
        let mut noncovalent = Vec::new();
        for b in &self.bonds {
            let a_idx = AtomId(b.atoms.first());
            let b_idx = AtomId(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent.push((a_idx, b_idx, NoncovalentBondAst::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        regular.push((a_idx, b_idx, b.try_into_ast(ctx)?));
                        continue;
                    }
                };
                let dative_bond = DativeBondAst::new(raise_bond_order(b.order));
                dative.push((vec![donor], acceptor, dative_bond));
            } else {
                regular.push((a_idx, b_idx, b.try_into_ast(ctx)?));
            }
        }

        let multicenter: Vec<(Vec<AtomId>, MulticenterBondAst)> = self
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

        let cfg = ground_atom_defaults();
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
            raise_ground_atom(atom, &cfg);
        }
        let constraints = Constraints::new();

        Ok(MoleculeAst::from_parts(
            atoms,
            regular,
            dative,
            vec![],
            multicenter,
            noncovalent,
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
                Some(m) => IsotopeAst::Lit(m as i64),
                None => IsotopeAst::Undetermined,
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
            spin: raise_atom_spin(self),
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
        bond.spin = raise_bond_spin(self);
        if matches!(self.order, TableBondOrder::Aromatic) {
            bond.constraints.add(BondConstraint::Aromatic);
        }
        Ok(bond)
    }
}

fn raise_atom_spin(atom: &TableAtom) -> SpinStateAst {
    match (atom.unpaired_electrons, atom.multiplicity) {
        (Some(u), Some(m)) => match SpinState::try_new(u, m) {
            Ok(s) => s.into(),
            Err(_) => SpinStateAst {
                unpaired: ValueAst::Lit(u as i64),
                multiplicity: ValueAst::Lit(u8::from(m) as i64),
            },
        },
        (Some(u), None) => SpinStateAst {
            unpaired: ValueAst::Lit(u as i64),
            multiplicity: ValueAst::Undetermined,
        },
        (None, Some(m)) => SpinStateAst {
            unpaired: ValueAst::Undetermined,
            multiplicity: ValueAst::Lit(u8::from(m) as i64),
        },
        (None, None) => SpinStateAst::default(),
    }
}

fn raise_bond_spin(bond: &TableBond) -> SpinStateAst {
    match (bond.unpaired_electrons, bond.multiplicity) {
        (Some(u), Some(m)) => match SpinState::try_new(u, m) {
            Ok(s) => s.into(),
            Err(_) => SpinStateAst {
                unpaired: ValueAst::Lit(u as i64),
                multiplicity: ValueAst::Lit(u8::from(m) as i64),
            },
        },
        (Some(u), None) => SpinStateAst {
            unpaired: ValueAst::Lit(u as i64),
            multiplicity: ValueAst::Undetermined,
        },
        (None, Some(m)) => SpinStateAst {
            unpaired: ValueAst::Undetermined,
            multiplicity: ValueAst::Lit(u8::from(m) as i64),
        },
        (None, None) => SpinStateAst::default(),
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

fn raise_ground_spin(
    spin: &mut SpinStateAst,
    u_mode: UnpairedElectronsDefault,
    m_mode: MultiplicityDefault,
) {
    let u = mem::replace(&mut spin.unpaired, ValueAst::Undetermined);
    let m = mem::replace(&mut spin.multiplicity, ValueAst::Undetermined);
    let resolved_u = if matches!(u, ValueAst::Undetermined) {
        match u_mode {
            UnpairedElectronsDefault::Zero => ValueAst::Lit(0),
            UnpairedElectronsDefault::Required => ValueAst::Undetermined,
            UnpairedElectronsDefault::Derived => match &m {
                ValueAst::Lit(mm) => ValueAst::Lit(mm - 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        u
    };
    let resolved_m = if matches!(m, ValueAst::Undetermined) {
        match m_mode {
            MultiplicityDefault::Required => ValueAst::Undetermined,
            MultiplicityDefault::Derived => match &resolved_u {
                ValueAst::Lit(uu) => ValueAst::Lit(uu + 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        m
    };
    spin.unpaired = resolved_u;
    spin.multiplicity = resolved_m;
}

fn raise_ground_atom(ast: &mut AtomAst, cfg: &AtomDefaults) {
    let AtomAst {
        element: _,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        spin,
        constraints,
    } = ast;

    if matches!(*isotope_mass, IsotopeAst::Undetermined) {
        *isotope_mass = match cfg.isotope {
            IsotopeDefault::Natural => IsotopeAst::Natural,
            IsotopeDefault::Required => IsotopeAst::Undetermined,
        };
    }
    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*implicit_hydrogens, ValueAst::Undetermined) {
        *implicit_hydrogens = match cfg.implicit_hydrogens {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*lone_pairs, ValueAst::Undetermined) {
        *lone_pairs = match cfg.lone_pairs {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_ground_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
    raise_ground_atom_constraints(constraints, cfg);
}

fn raise_ground_atom_constraints(constraints: &mut AtomConstraints, cfg: &AtomDefaults) {
    constraints.retain(|c| !c.is_undetermined());

    for kind in AtomConstraintKind::iter() {
        match kind {
            AtomConstraintKind::Valence => {
                if matches!(cfg.valence, NumericDefault::Zero) && !constraints.contains(kind) {
                    constraints.add(AtomConstraint::Valence(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::DonatedPairs => {
                if matches!(cfg.donated_pairs, NumericDefault::Zero) && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::DonatedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AcceptedPairs => {
                if matches!(cfg.accepted_pairs, NumericDefault::Zero) && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AromaticValence => {
                if !constraints.contains(kind) {
                    match cfg.aromatic_valence {
                        AromaticValenceDefault::NotAromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::NotAromatic,
                            ));
                        }
                        AromaticValenceDefault::Aromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::Aromatic(ValueAst::Undetermined),
                            ));
                        }
                        AromaticValenceDefault::Required => {}
                    }
                }
            }
            AtomConstraintKind::MulticenterValence => {
                if !constraints.contains(kind) {
                    match cfg.multicenter_valence {
                        MulticenterValenceDefault::NotMulticenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::NotMulticenter,
                            ));
                        }
                        MulticenterValenceDefault::Multicenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::Multicenter(ValueAst::Undetermined),
                            ));
                        }
                        MulticenterValenceDefault::Required => {}
                    }
                }
            }
            AtomConstraintKind::TotalValence
            | AtomConstraintKind::Degree
            | AtomConstraintKind::TotalDegree
            | AtomConstraintKind::RingDegree
            | AtomConstraintKind::RingValence
            | AtomConstraintKind::TotalHydrogens
            | AtomConstraintKind::RingCount
            | AtomConstraintKind::RingSize
            | AtomConstraintKind::JointDomain => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::BondId;
    use umol_shared::element::Element;

    use super::*;
    use crate::io::ctfile::config::CtfileIoConfig;
    use crate::io::ctfile::{parse_mol_bytes_with, parse_mol_to_ast};
    use crate::io::smiles::parse_smiles_to_ast;
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, ElementScope, RingLimits, ValenceModel,
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

    fn counts_valence() -> CountsValence {
        CountsValence::new(ValenceTable::default_table().clone())
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
        assert_eq!(ast.atom(AtomId(0)).ast.constraints.aromatic_valence(), expected);
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
        counts_valence().resolve(&mut ast).unwrap();
        assert_eq!(ast.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(ast.atom(AtomId(i)).ast.to_string(), expected_atom);
        }
    }

    #[rstest]
    fn test_parse_smiles_to_ast_methane_counts_resolve() {
        let mut ast = parse_smiles_to_ast("C").unwrap();
        counts_valence().resolve(&mut ast).unwrap();
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
            valence: ValenceModel::Counts {
                table: ValenceTable::default_table().clone(),
            },
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
