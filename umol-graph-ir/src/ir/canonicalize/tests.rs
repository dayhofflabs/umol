//! Unit and bounded-exhaustive tests for aggregate canonicalization.
//!
//! These tests exercise exact typed-key ordering, adapter projection, search invariance,
//! stereo-frame transport, canonical representatives, and aggregate failure boundaries.

use std::array;
use std::cell::Cell;
use std::cmp::Ordering;

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_edn::FromEdn;

use super::*;
use crate::ir::{
    ArithExpr, AromaticSystemConstraintForm, AromaticSystemDelta, AromaticSystemFieldChange,
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomDelta, AtomFieldChange, AtomForm,
    AtomId, BondConstraintForm, BondDelta, BondFieldChange, BondForm, BondId, BooleanForm,
    Constraint, ConstraintDelta, ConstraintSpan, DativeBondConstraintForm, DativeBondDelta,
    DativeBondFieldChange, DativeBondForm, DativeBondId, Deltas, EntitySpan, MoleculeConstraint,
    MoleculeEntries, MulticenterBondConstraintForm, MulticenterBondDelta,
    MulticenterBondFieldChange, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondConstraintForm, NoncovalentBondDelta, NoncovalentBondFieldChange,
    NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, NoncovalentBondKindForm, NumForm,
    ReactionSpanEntries, RelationalConstraint, StereoAtomConstraintForm, StereoAtomDelta,
    StereoAtomFieldChange, StereoAtomForm, StereoAtomId, StereoBondConstraintForm, StereoBondDelta,
    StereoBondFieldChange, StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, Stereogenicity,
    StereogenicityForm, Topicity, TopicityForm, TopicityRelationForm,
};
use crate::mol_dsl;

impl AutomorphismAdapter {
    pub(super) fn automorphisms(
        &self,
        algorithm: AutomorphismAlgorithm,
    ) -> AutomorphismAdapterOutput {
        let output = self
            .graph()
            .automorphisms(|node| self.color(node), algorithm);
        self.project_automorphisms(&output)
    }
}

#[fixture]
pub(super) fn initial_color_molecule() -> Molecule {
    let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(1),
        ArithExpr::Lit(2),
    ])));
    let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(0),
        ArithExpr::Lit(1),
    ])));
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::C)
                .with_charge(normalized_three)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
            AtomForm::from_element(Element::O).with_charge(3_i64),
            AtomForm::from_element(Element::C).with_charge(4_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
            (
                vec![AtomId(1)],
                AtomId(4),
                DativeBondForm::new(normalized_one),
            ),
            (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(1), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(2), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

#[fixture]
pub(super) fn canonicalize_context() -> CanonicalizeContext {
    CanonicalizeContext {
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

#[fixture]
pub(super) fn stereo_atom_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

#[fixture]
pub(super) fn stereo_bond_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        ..Default::default()
    })
}

#[fixture]
pub(super) fn stereo_constraint_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

/// Eight outer stereo atoms in two mirrored groups of four, each site bearing a distinct stereo
/// descriptor so partition refinement can separate them. Only three stereo kinds are admissible on
/// an atom at degree 4, so the fourth descriptor is a second `Tetrahedral` coset rather than a
/// fourth kind.
#[fixture]
pub(super) fn para_stereo_canonicalization_molecule() -> Molecule {
    let outer_ligands = vec![
        StereoLigand::new(AtomId(10), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(11), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(12), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(13), StereoLigandKind::Atom),
    ];

    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(AtomForm::from_element)
        .collect(),
        bonds: (2..=5)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .chain((6..=9).map(|ligand| (AtomId(1), AtomId(ligand), BondForm::from_order(1))))
            .chain((2..=9).flat_map(|site| {
                (10..=13).map(move |ligand| (AtomId(site), AtomId(ligand), BondForm::from_order(1)))
            }))
            .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [2, 3, 4, 5]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                [6, 8, 7, 9]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(2),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(3),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            (
                AtomId(4),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(5),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
            (
                AtomId(6),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(7),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            (
                AtomId(8),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(9),
                outer_ligands,
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
        ],
        ..Default::default()
    })
}

#[rstest]
#[case::topology_constitution(DescriptionLevel::Topology, DescriptionLevel::Constitution)]
#[case::constitution_structure(DescriptionLevel::Constitution, DescriptionLevel::Structure)]
#[case::structure_full(DescriptionLevel::Structure, DescriptionLevel::Full)]
fn test_description_level_order(#[case] lower: DescriptionLevel, #[case] upper: DescriptionLevel) {
    assert_eq!(lower.cmp(&upper), Ordering::Less);
    assert_eq!(lower.max(upper), upper);
}

#[rstest]
#[case::atom_add(
    Delta::Atom(AtomDelta::Add {
        id: AtomId(0),
        attributes: AtomForm::from_element(Element::C),
    }),
    DescriptionLevel::Topology,
)]
#[case::atom_remove(
    Delta::Atom(AtomDelta::Remove {
        id: AtomId(0),
        attributes: AtomForm::from_element(Element::C),
    }),
    DescriptionLevel::Topology,
)]
#[case::atom_modify_field(
    Delta::Atom(AtomDelta::ModifyField {
        id: AtomId(0),
        change: AtomFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Lit(1) },
    }),
    DescriptionLevel::Topology,
)]
#[case::atom_modify_constraint(
    Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: None, new: None }),
    DescriptionLevel::Topology,
)]
#[case::atom_inline_constraint(
    Delta::Atom(AtomDelta::Add {
        id: AtomId(0),
        attributes: AtomForm::from_element(Element::C)
            .with_constraint(AtomConstraintForm::valence(4)),
    }),
    DescriptionLevel::Topology,
)]
#[case::bond_add(
    Delta::Bond(BondDelta::Add {
        id: BondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: BondForm::from_order(1),
    }),
    DescriptionLevel::Topology,
)]
#[case::bond_remove(
    Delta::Bond(BondDelta::Remove {
        id: BondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: BondForm::from_order(1),
    }),
    DescriptionLevel::Topology,
)]
#[case::bond_modify_field(
    Delta::Bond(BondDelta::ModifyField {
        id: BondId(0),
        change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    }),
    DescriptionLevel::Topology,
)]
#[case::bond_modify_constraint(
    Delta::Bond(BondDelta::ModifyConstraint { id: BondId(0), old: None, new: None }),
    DescriptionLevel::Topology,
)]
#[case::bond_inline_constraint(
    Delta::Bond(BondDelta::Add {
        id: BondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: BondForm::from_order(1)
            .with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    }),
    DescriptionLevel::Topology,
)]
#[case::dative_add(
    Delta::DativeBond(DativeBondDelta::Add {
        id: DativeBondId(0),
        donors: vec![AtomId(0)],
        acceptor: AtomId(1),
        attributes: DativeBondForm::from_order(1),
    }),
    DescriptionLevel::Constitution,
)]
#[case::dative_remove(
    Delta::DativeBond(DativeBondDelta::Remove {
        id: DativeBondId(0),
        donors: vec![AtomId(0)],
        acceptor: AtomId(1),
        attributes: DativeBondForm::from_order(1),
    }),
    DescriptionLevel::Constitution,
)]
#[case::dative_modify_field(
    Delta::DativeBond(DativeBondDelta::ModifyField {
        id: DativeBondId(0),
        change: DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    }),
    DescriptionLevel::Constitution,
)]
#[case::dative_modify_constraint(
    Delta::DativeBond(DativeBondDelta::ModifyConstraint {
        id: DativeBondId(0),
        old: None,
        new: None,
    }),
    DescriptionLevel::Constitution,
)]
#[case::dative_inline_constraint(
    Delta::DativeBond(DativeBondDelta::Add {
        id: DativeBondId(0),
        donors: vec![AtomId(0)],
        acceptor: AtomId(1),
        attributes: DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    }),
    DescriptionLevel::Constitution,
)]
#[case::aromatic_add(
    Delta::AromaticSystem(AromaticSystemDelta::Add {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: AromaticSystemForm::from_electrons(vec![1, 1]),
    }),
    DescriptionLevel::Constitution,
)]
#[case::aromatic_remove(
    Delta::AromaticSystem(AromaticSystemDelta::Remove {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: AromaticSystemForm::from_electrons(vec![1, 1]),
    }),
    DescriptionLevel::Constitution,
)]
#[case::aromatic_modify_field(
    Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
        id: AromaticSystemId(0),
        change: AromaticSystemFieldChange::Charge {
            old: NumForm::Lit(0),
            new: NumForm::Lit(1),
        },
    }),
    DescriptionLevel::Constitution,
)]
#[case::aromatic_modify_constraint(
    Delta::AromaticSystem(AromaticSystemDelta::ModifyConstraint {
        id: AromaticSystemId(0),
        old: None,
        new: None,
    }),
    DescriptionLevel::Constitution,
)]
#[case::aromatic_inline_constraint(
    Delta::AromaticSystem(AromaticSystemDelta::Add {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: AromaticSystemForm::from_electrons(vec![1, 1])
            .with_constraint(AromaticSystemConstraintForm::electron_count(2)),
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_add(
    Delta::MulticenterBond(MulticenterBondDelta::Add {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: MulticenterBondForm::from_electrons(vec![1, 1]),
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_remove(
    Delta::MulticenterBond(MulticenterBondDelta::Remove {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: MulticenterBondForm::from_electrons(vec![1, 1]),
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_modify_field(
    Delta::MulticenterBond(MulticenterBondDelta::ModifyField {
        id: MulticenterBondId(0),
        change: MulticenterBondFieldChange::Charge {
            old: NumForm::Lit(0),
            new: NumForm::Lit(1),
        },
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_modify_constraint(
    Delta::MulticenterBond(MulticenterBondDelta::ModifyConstraint {
        id: MulticenterBondId(0),
        old: None,
        new: None,
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_inline_constraint(
    Delta::MulticenterBond(MulticenterBondDelta::Add {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)],
        attributes: MulticenterBondForm::from_electrons(vec![1, 1])
            .with_constraint(MulticenterBondConstraintForm::electron_count(2)),
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_add(
    Delta::NoncovalentBond(NoncovalentBondDelta::Add {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_remove(
    Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_modify_field(
    Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
        id: NoncovalentBondId(0),
        change: NoncovalentBondFieldChange::Kind {
            old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
        },
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_modify_constraint(
    Delta::NoncovalentBond(NoncovalentBondDelta::ModifyConstraint {
        id: NoncovalentBondId(0),
        old: None,
        new: None,
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_inline_constraint(
    Delta::NoncovalentBond(NoncovalentBondDelta::Add {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
    }),
    DescriptionLevel::Constitution,
)]
#[case::stereo_atom_add(
    Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: vec![],
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_atom_remove(
    Delta::StereoAtom(StereoAtomDelta::Remove {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: vec![],
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_atom_modify_field(
    Delta::StereoAtom(StereoAtomDelta::ModifyField {
        id: StereoAtomId(0),
        change: StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32),
            new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32),
        },
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_atom_modify_constraint(
    Delta::StereoAtom(StereoAtomDelta::ModifyConstraint {
        id: StereoAtomId(0),
        kind: Some(StereoKind::Tetrahedral),
        old: None,
        new: None,
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_atom_inline_constraint(
    Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: vec![],
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_bond_add(
    Delta::StereoBond(StereoBondDelta::Add {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: vec![],
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_bond_remove(
    Delta::StereoBond(StereoBondDelta::Remove {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: vec![],
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_bond_modify_field(
    Delta::StereoBond(StereoBondDelta::ModifyField {
        id: StereoBondId(0),
        change: StereoBondFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32),
            new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32),
        },
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_bond_modify_constraint(
    Delta::StereoBond(StereoBondDelta::ModifyConstraint {
        id: StereoBondId(0),
        kind: Some(StereoKind::CisTrans),
        old: None,
        new: None,
    }),
    DescriptionLevel::Structure,
)]
#[case::stereo_bond_inline_constraint(
    Delta::StereoBond(StereoBondDelta::Add {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: vec![],
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
    }),
    DescriptionLevel::Structure,
)]
#[case::constraint(
    Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
        MoleculeConstraint::Connected { atoms: None },
    ))),
    DescriptionLevel::Full,
)]
fn test_delta_canonicalize_level(#[case] delta: Delta, #[case] expected: DescriptionLevel) {
    assert_eq!(delta_canonicalize_level(&delta), expected);
    assert_eq!(delta_canonicalize_level(&delta.inverse()), expected);
}

#[rstest]
#[case::modified_rhs(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Modified {
            lhs: AtomForm::from_element(Element::C),
            rhs: AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        }],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
)]
fn test_reaction_span_canonicalize_level_modified(
    #[case] span: ReactionSpan,
    #[case] expected: DescriptionLevel,
) {
    assert_eq!(reaction_span_canonicalize_level(&span), expected);
}

#[rstest]
#[case::dative_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::B),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::B),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    }),
)]
#[case::stereo_atom(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        noncovalent: vec![([AtomId(1), AtomId(2)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        noncovalent: vec![([AtomId(1), AtomId(2)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
        )],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        ..Default::default()
    }),
)]
#[case::atom_constraint(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        ..Default::default()
    }),
)]
fn test_canonical_eq_additional_entity(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
) {
    let left_reaction = Reaction::new(left.clone(), Deltas::new());
    let right_reaction = Reaction::new(right.clone(), Deltas::new());
    let left_span = left_reaction
        .clone()
        .to_reaction_span()
        .expect("identity reaction materializes");
    let right_span = right_reaction
        .clone()
        .to_reaction_span()
        .expect("identity reaction materializes");

    assert_eq!(
        (
            left.canonical_eq(&right, &canonicalize_context),
            left_reaction.canonical_eq(&right_reaction, &canonicalize_context),
            left_span.canonical_eq(&right_span, &canonicalize_context),
        ),
        (false, false, false),
    );
    assert_eq!(
        (
            right.canonical_eq(&left, &canonicalize_context),
            right_reaction.canonical_eq(&left_reaction, &canonicalize_context),
            right_span.canonical_eq(&left_span, &canonicalize_context),
        ),
        (false, false, false),
    );
}

#[rstest]
fn test_reaction_span_canonicalize_contradiction(canonicalize_context: CanonicalizeContext) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set(Vec::<i64>::new())),
        )],
        ..Default::default()
    });

    assert_eq!(
        span.canonicalize(&canonicalize_context),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
#[case::atom_charge(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_absent(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_pair(
    r#"{:atoms ["C#c+" "C#c+" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C#c+" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_absent(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_pair(
    r#"{:atoms ["C#v4" "C#v4" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C#v4" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "2"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1#a"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::dative_bond(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [1] :acceptor 4 :attrs "1"}]}"#,
    true
)]
#[case::dative_bond_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 2 :attrs "1"}]}"#,
    false
)]
#[case::dative_bond_role(
    r#"{:atoms ["N" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["N" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [3] :acceptor 0 :attrs "1"}]}"#,
    false
)]
#[case::aromatic_electron_counts(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}]}"#,
    true
)]
#[case::aromatic_electron_counts_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
    false
)]
#[case::aromatic_electron_counts_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}]}"#,
    false
)]
#[case::aromatic_system_charge(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#c-"}]}"#,
    true
)]
#[case::aromatic_system_charge_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#e2"}]}"#,
    true
)]
#[case::aromatic_system_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::multicenter_electron_counts(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[0,0,2]"}]}"#,
    true
)]
#[case::multicenter_electron_counts_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}]}"#,
    false
)]
#[case::noncovalent_bond(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [1 4] :attrs "Hbd"}]}"#,
    true
)]
#[case::noncovalent_bond_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#,
    false
)]
#[case::stereo_configuration(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    true
)]
#[case::stereo_configuration_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_stereogenicity(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}]}"#,
    true
)]
#[case::stereo_stereogenicity_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}]}"#,
    true
)]
#[case::stereo_topicity_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::molecule_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_opposite(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [1 4] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_adjacent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3 2] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_count(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::molecule_constraint_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::molecule_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
fn test_molecule_canonicalize_placement(
    canonicalize_context: CanonicalizeContext,
    #[case] lhs: &str,
    #[case] rhs: &str,
    #[case] expected: bool,
) {
    let [lhs, rhs] =
        [lhs, rhs].map(|source| Molecule::from_edn_str(source).expect("case molecule parses"));

    assert_eq!(
        lhs.clone().canonicalize(&canonicalize_context)
            == rhs.clone().canonicalize(&canonicalize_context),
        expected
    );
    assert_eq!(lhs.canonical_eq(&rhs, &canonicalize_context), expected);
    if expected {
        assert_eq!(
            lhs.canonical_hash(&canonicalize_context),
            rhs.canonical_hash(&canonicalize_context),
        );
    }
}

#[rstest]
fn test_molecule_canonicalize_stereo_frame(
    stereo_constraint_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    let mut editor = stereo_constraint_canonicalization_molecule.edit();
    editor
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints = constraint.clone().into();
    *editor.constraints_mut() =
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
    let source = editor.finish().unwrap();
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        reframed.canonicalize(&canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical.clone().canonicalize(&canonicalize_context),
        Ok(canonical),
    );
}

#[rstest]
#[case::atom_charge(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C).with_charge(NumForm::lit_set(Vec::<i64>::new()))],
    ..Default::default()
}))]
#[case::atom_constraint(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C)
        .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())))],
    ..Default::default()
}))]
#[case::dative_order(Molecule::from_entries(MoleculeEntries {
    atoms: vec![
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::B),
    ],
    dative: vec![(
        vec![AtomId(0)],
        AtomId(1),
        DativeBondForm::new(NumForm::lit_set(Vec::<i64>::new())),
    )],
    ..Default::default()
}))]
#[case::stereo_coset(Molecule::from_entries(MoleculeEntries {
    atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
        .map(AtomForm::from_element)
        .into(),
    bonds: (1..=4)
        .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
        .collect(),
    stereo_atoms: vec![(
        AtomId(0),
        (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
    )],
    ..Default::default()
}))]
#[case::electron_count_length(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 3],
    aromatic: vec![(
        vec![AtomId(0), AtomId(1), AtomId(2)],
        AromaticSystemForm::from_electrons(vec![2, 1]),
    )],
    ..Default::default()
}))]
fn test_molecule_canonicalize_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
) {
    assert_eq!(
        molecule.canonicalize(&canonicalize_context),
        Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
fn test_molecule_canonical_eq(
    initial_color_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let renumbered = initial_color_molecule
        .clone()
        .remap(&reverse_remapping(&initial_color_molecule));

    assert!(initial_color_molecule.canonical_eq(&renumbered, &canonicalize_context));
}

#[rstest]
#[case::ethanol_ethylamine(
    mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    false
)]
#[case::hexane_2_methylpentane(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    false
)]
#[case::bromochlorofluoroiodomethane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::bromochlorofluoroiodomethane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
    ),
    false
)]
#[case::decalin_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[9 8 "1"] [8 7 "1"] [7 6 "1"] [6 5 "1"] [5 4 "1"] [4 9 "1"]
                    [4 3 "1"] [3 2 "1"] [2 1 "1"] [1 0 "1"] [0 5 "1"]]}"#
    ),
    true
)]
#[case::cyclohexane_two_cyclopropanes(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"] [3 4 "1"] [4 5 "1"] [5 3 "1"]]}"#),
    false
)]
#[case::cuneane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[7 6 "1"] [6 5 "1"] [5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"]
                    [1 0 "1"] [0 7 "1"] [7 3 "1"] [6 4 "1"] [5 1 "1"] [2 0 "1"]]}"#
    ),
    true
)]
#[case::cubane_cuneane(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 4 "1"] [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    false
)]
#[case::meso_dichlorobutane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::meso_and_chiral_dichlorobutane(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    false
)]
#[case::dichlorobutane_enantiomers(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    false
)]
#[case::charge_sum_rotated(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(3), AtomId(4)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    true
)]
#[case::charge_sum_other_atoms(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(2)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    false
)]
fn test_molecule_canonical_eq_compounds(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    assert_eq!(right.canonical_eq(&left, &canonicalize_context), expected);
}

#[rstest]
#[case::same_frame([4, 2, 1, 3], 0, true)]
#[case::odd_frame([2, 4, 1, 3], 1, true)]
#[case::even_frame([2, 1, 4, 3], 0, true)]
#[case::opposite_odd_frame([2, 4, 1, 3], 0, false)]
#[case::opposite_even_frame([2, 1, 4, 3], 1, false)]
fn test_molecule_canonical_eq_stereo_atom(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] ligand_ids: [u32; 4],
    #[case] coset: u32,
    #[case] expected: bool,
    #[values(false, true)] relabel: bool,
) {
    let mut elements = [Element::C, Element::F, Element::Cl, Element::Br, Element::I];
    if relabel {
        elements.swap(1, 2);
    }
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: elements.into_iter().map(AtomForm::from_element).collect(),
        bonds: (1..=4)
            .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            ligand_ids
                .into_iter()
                .map(|atom| {
                    let atom = match (relabel, atom) {
                        (true, 1) => 2,
                        (true, 2) => 1,
                        _ => atom,
                    };
                    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
                })
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, coset),
        )],
        ..Default::default()
    });
    let source = stereo_atom_canonicalization_molecule;
    let canonical = source.clone().canonicalize(&canonicalize_context).unwrap();
    let other_canonical = other.clone().canonicalize(&canonicalize_context).unwrap();

    if expected {
        assert_eq!(canonical, other_canonical);
    } else {
        assert_ne!(canonical, other_canonical);
    }
    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::same_frame([3, 2, 5, 4], 1, true)]
#[case::one_end([2, 3, 5, 4], 0, true)]
#[case::both_ends([2, 3, 4, 5], 1, true)]
#[case::exchanged_ends([5, 4, 3, 2], 1, true)]
#[case::opposite([2, 3, 5, 4], 1, false)]
fn test_molecule_canonical_eq_stereo_bond(
    stereo_bond_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] ligand_ids: [u32; 4],
    #[case] coset: u32,
    #[case] expected: bool,
    #[values(false, true)] relabel: bool,
) {
    let mut elements = [
        Element::C,
        Element::C,
        Element::F,
        Element::Cl,
        Element::Br,
        Element::I,
    ];
    if relabel {
        elements.swap(2, 3);
    }
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: elements.into_iter().map(AtomForm::from_element).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            ligand_ids
                .into_iter()
                .map(|atom| {
                    let atom = match (relabel, atom) {
                        (true, 2) => 3,
                        (true, 3) => 2,
                        _ => atom,
                    };
                    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
                })
                .collect(),
            StereoBondForm::new(StereoKind::CisTrans, coset),
        )],
        ..Default::default()
    });
    let source = stereo_bond_canonicalization_molecule;
    let canonical = source.clone().canonicalize(&canonicalize_context).unwrap();
    let other_canonical = other.clone().canonicalize(&canonicalize_context).unwrap();

    if expected {
        assert_eq!(canonical, other_canonical);
    } else {
        assert_ne!(canonical, other_canonical);
    }
    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::renumbering(
    mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["N" "C"] :bonds [[0 1 "1"]]}"#),
    true
)]
#[case::atom_form(
    mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["F" "C"] :bonds [[0 1 "1"]]}"#),
    false
)]
#[case::ammonia_borane_renumbered(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    mol_dsl!(r#"{:atoms ["B" "N"] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#),
    true
)]
#[case::methoxyamine_methylhydroxylamine(
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 2 "1"] [2 1 "1"]]}"#),
    false
)]
#[case::bromochlorofluoroiodomethane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::bromochlorofluoroiodomethane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
    ),
    false
)]
fn test_molecule_canonical_eq_remapping_witness(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert!(left.clone().canonicalize(&canonicalize_context).is_ok());
    assert!(right.clone().canonicalize(&canonicalize_context).is_ok());
    assert_eq!(explicitly_dense_equivalent(&left, &right), expected);
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::hydrogen(StereoLigandKind::ImplicitHydrogen, false)]
#[case::lone_pair(StereoLigandKind::LonePair, false)]
#[case::virtual_ligand(StereoLigandKind::Virtual, true)]
fn test_molecule_canonical_eq_virtual_ligand(
    canonicalize_context: CanonicalizeContext,
    #[case] kind: StereoLigandKind,
    #[case] expected: bool,
) {
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    });
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: [Element::Br, Element::C, Element::F, Element::Cl]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: vec![
            (AtomId(1), AtomId(0), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), kind),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        ..Default::default()
    });

    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::one_fails(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_identical(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_equal_normal_forms(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_different_frames(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(0), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different_counts(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([])),
            AtomForm::from_element(Element::O),
        ],
        ..Default::default()
    }),
    false
)]
fn test_molecule_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::atom_charge_unchanged(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_unchanged_absent(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_modified(
    r#"{:atoms [{:modify ["C" "C#c+"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" {:modify ["C" "C#c+"]} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_modified_reversed(
    r#"{:atoms [{:modify ["C" "C#c+"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms [{:modify ["C#c+" "C"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [3 6 "1"]}]}"#,
    true
)]
#[case::atom_charge_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_charge_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [3 6 "1"]}]}"#,
    true
)]
#[case::atom_charge_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_constraint_unchanged(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_unchanged_absent(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_modified(
    r#"{:atoms [{:modify ["C" "C#v4"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" {:modify ["C" "C#v4"]} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_modified_reversed(
    r#"{:atoms [{:modify ["C" "C#v4"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms [{:modify ["C#v4" "C"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [3 6 "1"]}]}"#,
    true
)]
#[case::atom_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [3 6 "1"]}]}"#,
    true
)]
#[case::atom_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::bond_order_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "2"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:modify [3 4 ["1" "2"]]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "2"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["2" "1"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "2"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "2"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1#a"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "1#a"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:modify [3 4 ["1" "1#a"]]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "1#a"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1#a" "1"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "1#a"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "1#a"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_added(
    r#"{:atoms [{:add "C"} "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:add [0 5 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" {:add "C"} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] {:add [2 3 "1"]} {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_added_removed(
    r#"{:atoms [{:add "C"} "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:add [0 5 "1"]}]}"#,
    r#"{:atoms [{:remove "C"} "C" "C" "C" "C" "C"]
          :bonds
          [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:remove [0 5 "1"]}]}"#,
    false
)]
#[case::bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} {:add [1 2 "1"]} [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::dative_bond_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [1] :acceptor 4 :attrs "1"}]}"#,
    true
)]
#[case::dative_bond_unchanged_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 2 :attrs "1"}]}"#,
    false
)]
#[case::dative_bond_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["1" "2"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [1] :acceptor 4 :attrs ["1" "2"]}}]}"#,
    true
)]
#[case::dative_bond_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["1" "2"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["2" "1"]}}]}"#,
    false
)]
#[case::dative_bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [1] :acceptor 4 :attrs "1"}}]}"#,
    true
)]
#[case::dative_bond_added_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 2 :attrs "1"}}]}"#,
    false
)]
#[case::dative_bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [1] :acceptor 4 :attrs "1"}}]}"#,
    true
)]
#[case::dative_bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    false
)]
#[case::aromatic_electron_counts_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}]}"#,
    true
)]
#[case::aromatic_electron_counts_unchanged_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}]}"#,
    false
)]
#[case::aromatic_electron_counts_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[1,1,1,1,1,1]" "[2,1,1,1,1,1]"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[1,1,1,1,1,1]" "[1,1,1,2,1,1]"]}}]}"#,
    true
)]
#[case::aromatic_electron_counts_modified_shift(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[2,1,1,1,1,1]" "[1,2,1,1,1,1]"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[2,1,1,1,1,1]" "[1,1,1,2,1,1]"]}}]}"#,
    false
)]
#[case::aromatic_electron_counts_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}}]}"#,
    true
)]
#[case::aromatic_electron_counts_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}}]}"#,
    false
)]
#[case::aromatic_electron_counts_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}}]}"#,
    true
)]
#[case::aromatic_electron_counts_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    false
)]
#[case::aromatic_system_charge_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#c-"}]}"#,
    true
)]
#[case::aromatic_system_charge_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#c-"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:modify {:atoms [3 4 5] :attrs ["*" "*#c-"]}}]}"#,
    true
)]
#[case::aromatic_system_charge_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#c-"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*#c-" "*"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:add {:atoms [3 4 5] :attrs "*#c-"}}]}"#,
    true
)]
#[case::aromatic_system_charge_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:remove {:atoms [3 4 5] :attrs "*#c-"}}]}"#,
    true
)]
#[case::aromatic_system_charge_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#e2"}]}"#,
    true
)]
#[case::aromatic_system_constraint_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#e2"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:modify {:atoms [3 4 5] :attrs ["*" "*#e2"]}}]}"#,
    true
)]
#[case::aromatic_system_constraint_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#e2"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*#e2" "*"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:add {:atoms [3 4 5] :attrs "*#e2"}}]}"#,
    true
)]
#[case::aromatic_system_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:remove {:atoms [3 4 5] :attrs "*#e2"}}]}"#,
    true
)]
#[case::aromatic_system_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::multicenter_electron_counts_unchanged(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[0,0,2]"}]}"#,
    true
)]
#[case::multicenter_electron_counts_unchanged_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}]}"#,
    false
)]
#[case::multicenter_electron_counts_modified(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[0,2,0]"]}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[0,0,2]" "[2,0,0]"]}}]}"#,
    true
)]
#[case::multicenter_electron_counts_modified_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[0,2,0]"]}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[1,1,0]"]}}]}"#,
    false
)]
#[case::multicenter_electron_counts_added(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[0,0,2]"}}]}"#,
    true
)]
#[case::multicenter_electron_counts_added_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[1,1,0]"}}]}"#,
    false
)]
#[case::multicenter_electron_counts_removed(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[0,0,2]"}}]}"#,
    true
)]
#[case::multicenter_electron_counts_removed_added(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    false
)]
#[case::noncovalent_bond_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [1 4] :attrs "Hbd"}]}"#,
    true
)]
#[case::noncovalent_bond_unchanged_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#,
    false
)]
#[case::noncovalent_bond_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Hbd" "Xbd"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [1 4] :attrs ["Hbd" "Xbd"]}}]}"#,
    true
)]
#[case::noncovalent_bond_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Hbd" "Xbd"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Xbd" "Hbd"]}}]}"#,
    false
)]
#[case::noncovalent_bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [1 4] :attrs "Hbd"}}]}"#,
    true
)]
#[case::noncovalent_bond_added_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 2] :attrs "Hbd"}}]}"#,
    false
)]
#[case::noncovalent_bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [1 4] :attrs "Hbd"}}]}"#,
    true
)]
#[case::noncovalent_bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    false
)]
#[case::stereo_configuration_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    true
)]
#[case::stereo_configuration_unchanged_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_configuration_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th1"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th1"]}}]}"#,
    true
)]
#[case::stereo_configuration_modified_reversed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th1"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th1" "Th0"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_configuration_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0"}}]}"#,
    true
)]
#[case::stereo_configuration_added_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th1"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_configuration_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0"}}]}"#,
    true
)]
#[case::stereo_configuration_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_stereogenicity_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}]}"#,
    true
)]
#[case::stereo_stereogenicity_unchanged_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#g/"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th0#g/"]}}]}"#,
    true
)]
#[case::stereo_stereogenicity_modified_reversed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#g/"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0#g/" "Th0"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}}]}"#,
    true
)]
#[case::stereo_stereogenicity_added_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}}]}"#,
    true
)]
#[case::stereo_stereogenicity_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}]}"#,
    true
)]
#[case::stereo_topicity_unchanged_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,1)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th0#o(0,1)'"]}}]}"#,
    true
)]
#[case::stereo_topicity_modified_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,1)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,2)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}}]}"#,
    true
)]
#[case::stereo_topicity_added_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}}]}"#,
    true
)]
#[case::stereo_topicity_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3] :sum 0}}]}"#,
    true
)]
#[case::constraint_unchanged_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [3] :sum 0}}}]}"#,
    true
)]
#[case::constraint_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0 1] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0 3] :sum 0}}}]}"#,
    false
)]
#[case::constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [3] :sum 0}}}]}"#,
    true
)]
#[case::constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    false
)]
fn test_reaction_span_canonicalize_placement(
    canonicalize_context: CanonicalizeContext,
    #[case] left: &str,
    #[case] right: &str,
    #[case] expected: bool,
) {
    let [left, right] =
        [left, right].map(|source| ReactionSpan::from_edn_str(source).expect("case span parses"));

    assert_eq!(
        left.clone().canonicalize(&canonicalize_context)
            == right.clone().canonicalize(&canonicalize_context),
        expected
    );
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    if expected {
        assert_eq!(
            left.canonical_hash(&canonicalize_context),
            right.canonical_hash(&canonicalize_context),
        );
    }
}

fn reaction_canonicalization_fixture() -> Reaction {
    let atom = AtomForm::from_element(Element::C).with_charge(0_i64);
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom.clone(), atom],
            ..Default::default()
        }),
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    )
}

#[rustfmt::skip]
#[rstest]
#[case::molecule_subset(ReactionSpan::from_entries(ReactionSpanEntries {
    atoms: vec![
        EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::N)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
    ],
    constraints: vec![ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
        sum: NumForm::Lit(0),
    }))],
    ..Default::default()
}))]
#[case::relational_atom_list(ReactionSpan::from_entries(ReactionSpanEntries {
    atoms: vec![
        EntitySpan::Unchanged(AtomForm::from_element(Element::Fe)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::N)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
    ],
    dative: vec![(vec![AtomId(1), AtomId(2)], AtomId(0), EntitySpan::Unchanged(DativeBondForm::from_order(1)))],
    constraints: vec![ConstraintSpan::Unchanged(Constraint::Relational(RelationalConstraint::DativeBondDonors {
        bond: DativeBondId(0),
        atoms: vec![AtomId(1), AtomId(2)],
    }))],
    ..Default::default()
}))]
fn test_reaction_span_canonical_eq_renumbering(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
) {
    let images = span.entity_counts().map(|count| (0..count).rev().collect::<Vec<_>>());
    let renumbered = span.clone().remap(&molecule_remapping(&images));

    assert!(span.canonical_eq(&renumbered, &canonicalize_context));
    assert!(span.to_reaction().canonical_eq(&renumbered.to_reaction(), &canonicalize_context));
}

#[rstest]
#[case::one_fails(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_identical(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_equal_normal_forms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Modified { lhs: AtomForm::from_element(Element::C), rhs: AtomForm::from_element(Element::C) },
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_different_frames(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(0), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different_counts(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([]))),
            EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
        ],
        ..Default::default()
    }),
    false
)]
fn test_reaction_span_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: ReactionSpan,
    #[case] right: ReactionSpan,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::contradiction(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
                ..Default::default()
            }),
            [Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })].into_iter().collect(),
        ),
        ReactionCanonicalizeError::Contradiction(Contradiction),
    )]
fn test_reaction_canonicalize_error(
    canonicalize_context: CanonicalizeContext,
    #[case] reaction: Reaction,
    #[case] expected: ReactionCanonicalizeError,
) {
    assert_eq!(reaction.canonicalize(&canonicalize_context), Err(expected));
}

#[rstest]
#[case::stereo_atom_configuration(ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1))))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([]))),
        )],
        ..Default::default()
}))]
fn test_reaction_span_canonical_hash_error(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
) {
    assert_eq!(
        span.canonical_hash(&canonicalize_context),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
fn test_reaction_canonical_eq(canonicalize_context: CanonicalizeContext) {
    let source = reaction_canonicalization_fixture();
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed reaction canonicalizes");

    assert!(source.canonical_eq(&canonical, &canonicalize_context));
}

#[rstest]
#[case::one_fails(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    false
)]
#[case::both_fail_identical(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    true
)]
#[case::both_fail_equal_reframed_forms(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O).with_charge(0_i64),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
                AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1]),
            )],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O).with_charge(0_i64),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(0)],
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 2]),
            )],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    true
)]
#[case::both_fail_different(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(3),
            },
        })]),
    ),
    false
)]
fn test_reaction_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Reaction,
    #[case] right: Reaction,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

pub(super) fn permutations(count: usize) -> Vec<Vec<usize>> {
    fn visit(values: &mut [usize], position: usize, output: &mut Vec<Vec<usize>>) {
        if position == values.len() {
            output.push(values.to_vec());
            return;
        }
        for next in position..values.len() {
            values.swap(position, next);
            visit(values, position + 1, output);
            values.swap(position, next);
        }
    }

    let mut values = (0..count).collect::<Vec<_>>();
    let mut output = Vec::new();
    visit(&mut values, 0, &mut output);
    output
}

pub(super) fn explicitly_dense_equivalent(left: &Molecule, right: &Molecule) -> bool {
    fn visit(
        entity_kind: usize,
        permutations: &[Vec<Vec<usize>>; 8],
        images: &mut [Vec<usize>; 8],
        left: &Molecule,
        right: &Molecule,
    ) -> bool {
        if entity_kind == images.len() {
            return left.framed_eq_under(right, &molecule_remapping(images));
        }

        permutations[entity_kind].iter().any(|permutation| {
            images[entity_kind].clone_from(permutation);
            visit(entity_kind + 1, permutations, images, left, right)
        })
    }

    let left_counts = left.entity_counts();
    if left_counts != right.entity_counts() {
        return false;
    }
    let permutations = left_counts.map(permutations);
    visit(
        0,
        &permutations,
        &mut array::from_fn(|_| Vec::new()),
        left,
        right,
    )
}

pub(super) fn reverse_remapping(molecule: &Molecule) -> MoleculeRemapping {
    let images = molecule
        .entity_counts()
        .map(|count| (0..count).rev().collect::<Vec<_>>());
    molecule_remapping(&images)
}

fn direct_graph_adapter(source: &Graph) -> AutomorphismAdapter {
    AutomorphismAdapter {
        graph: source.clone(),
        colors: vec![AutomorphismColor::Entity(0); source.node_count()],
        node_sources: source.node_ids().map(SubdivisionNodeSource::Node).collect(),
        source_node_count: source.node_count(),
    }
}

#[rstest]
fn test_incidence_cmp() {
    let incidences = [
        Incidence::BondEndpoint,
        Incidence::DativeDonor,
        Incidence::DativeAcceptor,
        Incidence::AromaticAtom,
        Incidence::MulticenterAtom,
        Incidence::NoncovalentEndpoint,
        Incidence::StereoSite,
        Incidence::StereoLigand(StereoLigandKind::Atom),
        Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
        Incidence::StereoLigand(StereoLigandKind::LonePair),
    ];

    for pair in incidences.windows(2) {
        assert_eq!(pair[0].cmp(&pair[1]), Ordering::Less);
    }
}

#[rstest]
#[case::empty(Vec::new(), 2, vec![0, 1], vec![0, 1], 2)]
#[case::distinct(
    vec![(20, 0), (10, 1), (30, 2)],
    0,
    vec![0, 0, 0],
    vec![1, 0, 2],
    3
)]
#[case::equal(
    vec![(20, 0), (10, 1), (20, 2)],
    0,
    vec![0, 0, 0],
    vec![1, 0, 1],
    2
)]
#[case::first_above_zero(
    vec![(20, 2), (10, 3), (20, 4)],
    2,
    vec![0, 1, 0, 0, 0],
    vec![0, 1, 3, 2, 3],
    4
)]
fn test_rank_colors(
    #[case] colors: Vec<(i32, usize)>,
    #[case] first: u32,
    #[case] mut ranks: Vec<u32>,
    #[case] expected_ranks: Vec<u32>,
    #[case] expected_next: u32,
) {
    assert_eq!(rank_colors(colors, first, &mut ranks), expected_next);
    assert_eq!(ranks, expected_ranks);
}

#[rstest]
#[case::empty(vec![], vec![], Some(vec![]))]
#[case::reordered(vec![0, 0, 1], vec![0, 1, 0], Some(vec![(0, 0), (2, 1)]))]
#[case::offset(vec![5, 6, 6], vec![6, 5, 6], Some(vec![(0, 1), (1, 0)]))]
#[case::different_counts(vec![0, 0, 1], vec![0, 1, 1], None)]
#[case::different_lengths(vec![0], vec![0, 0], None)]
#[case::integer_above(vec![0, 1], vec![0, 2], None)]
#[case::integer_below(vec![1, 2], vec![0, 1], None)]
fn test_color_representatives(
    #[case] source: Vec<u32>,
    #[case] target: Vec<u32>,
    #[case] expected: Option<Vec<(usize, usize)>>,
) {
    assert_eq!(color_representatives(&source, &target), expected);
}

#[rstest]
#[case::empty(
    4,
    Vec::new(),
    vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
)]
#[case::one_generator(
    4,
    vec![vec![NodeId(1), NodeId(0), NodeId(2), NodeId(3)]],
    vec![NodeId(0), NodeId(0), NodeId(2), NodeId(3)],
)]
#[case::generated_subgroup(
    4,
    vec![
        vec![NodeId(1), NodeId(0), NodeId(2), NodeId(3)],
        vec![NodeId(0), NodeId(2), NodeId(1), NodeId(3)],
    ],
    vec![NodeId(0), NodeId(0), NodeId(0), NodeId(3)],
)]
fn test_source_orbits_from_generators(
    #[case] source_node_count: usize,
    #[case] generators: Vec<Vec<NodeId>>,
    #[case] expected: Vec<NodeId>,
) {
    assert_eq!(
        source_orbits_from_generators(source_node_count, &generators),
        expected,
    );
}

#[rstest]
#[case::ordered_classes(
        vec![2_u32, 0, 2, 1, 0],
        OrderedPartition {
            cells: vec![
                vec![NodeId(1), NodeId(4)],
                vec![NodeId(3)],
                vec![NodeId(0), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_from_descriptors(
    #[case] descriptors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(OrderedPartition::from_descriptors(&descriptors), expected,);
}

#[rstest]
#[case::empty(Vec::new(), OrderedPartition { cells: Vec::new() })]
#[case::repeated(
    vec![2, 0, 2, 1, 0],
    OrderedPartition {
        cells: vec![
            vec![NodeId(1), NodeId(4)],
            vec![NodeId(3)],
            vec![NodeId(0), NodeId(2)],
        ],
    },
)]
#[case::ordered_distinct(
    vec![0, 1, 2],
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
)]
fn test_ordered_partition_from_color_ranks(
    #[case] colors: Vec<u32>,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(OrderedPartition::from_color_ranks(colors), expected);
}

#[rstest]
#[case::empty(
    OrderedPartition { cells: Vec::new() },
    Graph::new(0, &[]),
    OrderedPartition { cells: Vec::new() },
)]
#[case::discrete(
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
    Graph::new(3, &[[0, 1], [1, 2]]),
    OrderedPartition {
        cells: vec![vec![NodeId(0)], vec![NodeId(1)], vec![NodeId(2)]],
    },
)]
#[case::unsplit(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]],
    },
    Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)]],
    },
)]
#[case::multiply_split(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4)]],
    },
    Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4]]),
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(1), NodeId(3)],
            vec![NodeId(0), NodeId(4)],
        ],
    },
)]
#[case::disconnected(
    OrderedPartition {
        cells: vec![vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4)]],
    },
    Graph::new(5, &[[0, 1], [2, 3]]),
    OrderedPartition {
        cells: vec![
            vec![NodeId(0), NodeId(1), NodeId(2), NodeId(3)],
            vec![NodeId(4)],
        ],
    },
)]
#[case::individualized(
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(0), NodeId(1), NodeId(3)],
        ],
    },
    Graph::new(4, &[[0, 1], [0, 2], [0, 3]]),
    OrderedPartition {
        cells: vec![
            vec![NodeId(2)],
            vec![NodeId(0)],
            vec![NodeId(1), NodeId(3)],
        ],
    },
)]
fn test_ordered_partition_refine(
    #[case] partition: OrderedPartition,
    #[case] graph: Graph,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(partition.clone().refine(&graph), expected);
    assert_eq!(partition.refine_by_active_cells(&graph), expected);
}

#[rstest]
#[case::naphthalene(
    10,
    &[
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
    ],
)]
#[case::disconnected_rings(
    12,
    &[
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6),
    ],
)]
#[case::long_path(
    77,
    &(0_u32..76)
        .map(|first| (first, first + 1))
        .collect::<Vec<_>>(),
)]
fn test_ordered_partition_refine_active_cell_equivalence(
    #[case] atom_count: usize,
    #[case] edges: &[(u32, u32)],
) {
    let graph = Graph::new(
        atom_count,
        &edges
            .iter()
            .map(|&(first, second)| [first, second])
            .collect::<Vec<_>>(),
    );
    let partition = OrderedPartition::from_descriptors(&vec![0_u32; atom_count]);

    assert_eq!(
        partition.clone().refine(&graph),
        partition.clone().refine_by_all_cells(&graph),
    );
    assert_eq!(
        partition.clone().refine_by_active_cells(&graph),
        partition.refine_by_all_cells(&graph),
    );
}

#[rstest]
#[case::first_cell(
        OrderedPartition {
            cells: vec![vec![NodeId(0), NodeId(3)], vec![NodeId(1), NodeId(2)]],
        },
        0,
        NodeId(3),
        OrderedPartition {
            cells: vec![
                vec![NodeId(3)],
                vec![NodeId(0)],
                vec![NodeId(1), NodeId(2)],
            ],
        },
    )]
fn test_ordered_partition_individualize(
    #[case] partition: OrderedPartition,
    #[case] cell_index: usize,
    #[case] node: NodeId,
    #[case] expected: OrderedPartition,
) {
    assert_eq!(partition.individualize(cell_index, node), expected);
}

#[rstest]
#[case::singleton(
    1,
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![],
        refinement_calls: 1,
        branch_order_calls: 0,
        backend_calls: 0,
        visited_leaves: 1,
        leaf_comparisons: 0,
        orbit_pruned_branches: 0,
    },
)]
#[case::orbit_pruned(
    2,
    CanonicalSearchStats {
        initial_residual_cell_sizes: vec![2],
        refinement_calls: 2,
        branch_order_calls: 1,
        backend_calls: 1,
        visited_leaves: 1,
        leaf_comparisons: 0,
        orbit_pruned_branches: 1,
    },
)]
fn test_canonical_search_stats(#[case] node_count: usize, #[case] expected: CanonicalSearchStats) {
    let source = Graph::new(node_count, &[]);
    let adapter = direct_graph_adapter(&source);
    let leaf_candidate = |order: &[NodeId]| CanonicalCandidate {
        key: order.to_vec(),
        entity_order: order.to_vec(),
    };
    let actual = canonical_search(
        &adapter,
        &adapter.colors,
        AutomorphismAlgorithm::Nauty,
        &leaf_candidate,
    );

    assert_eq!(actual.stats, expected);
}

#[rstest]
#[case::unbuilt(false, None, 0)]
#[case::compared(
    true,
    Some(vec![NodeId(1), NodeId(0)]),
    2,
)]
fn test_deferred_key_into_built(
    #[case] compared: bool,
    #[case] expected: Option<Vec<NodeId>>,
    #[case] expected_builds: usize,
) {
    let builds = Cell::new(0_usize);
    let build = |order: &[NodeId]| {
        builds.set(builds.get() + 1);
        order.to_vec()
    };
    let key = DeferredKey::new(&[NodeId(1), NodeId(0)], &build);
    if compared {
        assert_eq!(
            key.cmp(&DeferredKey::new(&[NodeId(0), NodeId(1)], &build)),
            Ordering::Greater,
        );
    }

    assert_eq!(key.into_built(), expected);
    assert_eq!(builds.get(), expected_builds);
}

#[rstest]
#[case::less(vec![NodeId(0), NodeId(1)], vec![NodeId(1), NodeId(0)], Ordering::Less)]
#[case::equal(vec![NodeId(0), NodeId(1)], vec![NodeId(0), NodeId(1)], Ordering::Equal)]
#[case::greater(vec![NodeId(1), NodeId(0)], vec![NodeId(0), NodeId(1)], Ordering::Greater)]
fn test_deferred_key_cmp(
    #[case] left: Vec<NodeId>,
    #[case] right: Vec<NodeId>,
    #[case] expected: Ordering,
) {
    let builds = Cell::new(0_usize);
    let build = |order: &[NodeId]| {
        builds.set(builds.get() + 1);
        order.to_vec()
    };
    let left = DeferredKey::new(&left, &build);
    let right = DeferredKey::new(&right, &build);

    assert_eq!(left.cmp(&right), expected);
    assert_eq!(left.cmp(&right), expected);
    assert_eq!(builds.get(), 2);
}

#[rstest]
#[case::contradiction(
    MoleculeCanonicalizeError::from(Contradiction),
    MoleculeCanonicalizeError::Contradiction(Contradiction)
)]
fn test_molecule_canonicalize_error_from(
    #[case] actual: MoleculeCanonicalizeError,
    #[case] expected: MoleculeCanonicalizeError,
) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case::contradiction(
    ReactionSpanCanonicalizeError::from(Contradiction),
    ReactionSpanCanonicalizeError::Contradiction(Contradiction)
)]
fn test_reaction_span_canonicalize_error_from(
    #[case] actual: ReactionSpanCanonicalizeError,
    #[case] expected: ReactionSpanCanonicalizeError,
) {
    assert_eq!(actual, expected);
}

#[rstest]
#[case::contradiction(
    ReactionCanonicalizeError::from(Contradiction),
    ReactionCanonicalizeError::Contradiction(Contradiction)
)]
fn test_reaction_canonicalize_error_from(
    #[case] actual: ReactionCanonicalizeError,
    #[case] expected: ReactionCanonicalizeError,
) {
    assert_eq!(actual, expected);
}
