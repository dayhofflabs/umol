//! Unit and bounded-exhaustive tests for aggregate canonicalization.
//!
//! These tests exercise exact typed-key ordering, adapter projection, search invariance,
//! stereo-frame transport, canonical representatives, and aggregate failure boundaries.

use std::array;
use std::cell::Cell;
use std::cmp::Ordering;

use rstest::{fixture, rstest};
use umol_chem::element::Element;

use super::*;
use crate::ir::{
    AromaticSystemConstraintForm, AromaticSystemDelta, AromaticSystemFieldChange,
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomDelta, AtomFieldChange, AtomForm,
    AtomId, BondConstraintForm, BondDelta, BondFieldChange, BondForm, BondId, BooleanForm,
    Constraint, ConstraintDelta, DativeBondConstraintForm, DativeBondDelta, DativeBondFieldChange,
    DativeBondForm, DativeBondId, Deltas, EntitySpan, Incidence, MoleculeConstraint,
    MoleculeEntries, MulticenterBondConstraintForm, MulticenterBondDelta,
    MulticenterBondFieldChange, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondConstraintForm, NoncovalentBondDelta, NoncovalentBondFieldChange,
    NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, NoncovalentBondKindForm, NumForm,
    ReactionSpanEntries, StereoAtomConstraintForm, StereoAtomDelta, StereoAtomFieldChange,
    StereoAtomForm, StereoAtomId, StereoBondConstraintForm, StereoBondDelta, StereoBondFieldChange,
    StereoBondForm, StereoBondId, StereoConfigurationForm, StereoKind, StereoLigand,
    StereoLigandKind, StereoVirtualLigands, Stereogenicity, StereogenicityForm,
};

#[fixture]
pub(crate) fn canonicalize_context() -> CanonicalizeContext {
    CanonicalizeContext {
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
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

pub(crate) fn permutations(count: usize) -> Vec<Vec<usize>> {
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

pub(crate) fn explicitly_dense_equivalent(
    counts: [usize; 8],
    framed_eq_under: impl Fn(&MoleculeRemapping) -> bool,
) -> bool {
    fn visit(
        entity_kind: usize,
        permutations: &[Vec<Vec<usize>>; 8],
        images: &mut [Vec<usize>; 8],
        framed_eq_under: &impl Fn(&MoleculeRemapping) -> bool,
    ) -> bool {
        if entity_kind == images.len() {
            return framed_eq_under(&molecule_remapping(images));
        }

        permutations[entity_kind].iter().any(|permutation| {
            images[entity_kind].clone_from(permutation);
            visit(entity_kind + 1, permutations, images, framed_eq_under)
        })
    }

    let permutations = counts.map(permutations);
    visit(
        0,
        &permutations,
        &mut array::from_fn(|_| Vec::new()),
        &framed_eq_under,
    )
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
        Incidence::StereoAnchor(StereoVirtualLigands::empty()),
        Incidence::StereoAnchor(StereoVirtualLigands::IMPLICIT_HYDROGEN),
        Incidence::StereoAnchor(StereoVirtualLigands::LONE_PAIR),
        Incidence::StereoAnchor(
            StereoVirtualLigands::IMPLICIT_HYDROGEN | StereoVirtualLigands::LONE_PAIR,
        ),
        Incidence::StereoAnchor(StereoVirtualLigands::VIRTUAL),
        Incidence::StereoLigand,
    ];

    for pair in incidences.windows(2) {
        assert_eq!(pair[0].cmp(&pair[1]), Ordering::Less);
    }
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
