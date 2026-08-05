use std::collections::HashSet;
use std::sync::Arc;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, BipartiteMaximumMatchingAlgorithm,
    ConnectedComponentsAlgorithm, Correspondence, EdgeId, GeneralMaximumMatchingAlgorithm,
    MatchingEnumerationAlgorithm, MaximumIndependentSetAlgorithm, NodeId, NonBipartiteGraphError,
    RelevantCycleEnumerationAlgorithm, ShortestCycleAlgorithm, SimpleCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
};

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::{AtomAst, ElementAst, IsotopeMassAst};
use super::super::bond::BondAst;
use super::super::boolean::BooleanAst;
use super::super::constraint::{
    AromaticSystemConstraintAst, AtomConstraintAst, AtomConstraintsAst, BondConstraintAst,
    BondConstraintsAst, Constraint, Constraints, DativeBondConstraintAst, DativeBondConstraintsAst,
    MoleculeConstraint, MulticenterBondConstraintAst, NoncovalentBondConstraintAst,
    RelationalConstraint, RingScope, StereoAtomConstraintAst, StereoBondConstraintAst,
    StereogenicityAst, SubPatternAnchor,
};
use super::super::correspondence::MoleculeCorrespondence;
use super::super::dative::DativeBondAst;
use super::super::edit::{AtomFieldChange, AtomHandle, BondHandle, Edit, Edits};
use super::super::electrons::ElectronCountsAst;
use super::super::entity::Entity;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
use super::super::ring::{RingConfig, RingModel, RingSetKind};
use super::super::spin::UnpairedElectronsAst;
use super::super::stereo::{StereoAtomAst, StereoBondAst, StereoCoset, StereoKind};
use super::super::value::ValueAst;
use super::{MoleculeAst, MoleculeEntries, MoleculeEntriesError, TransactionError};
use crate::{mol_dsl, mol_dsl_ground};

fn ground_atom() -> AtomAst {
    let mut a = AtomAst::from_element(Element::C);
    a.isotope_mass = IsotopeMassAst::Natural;
    a.charge = ValueAst::Lit(0);
    a.implicit_hydrogens = ValueAst::Lit(4);
    a.lone_pairs = ValueAst::Lit(0);
    a.unpaired_electrons = UnpairedElectronsAst::from((0_u8, 1_u8));
    a
}

fn constraints_with_molecule(c: Constraint) -> Constraints {
    let mut out = Constraints::new();
    out.push(c);
    out
}

#[rstest]
fn test_molecule_ast_new() {
    let m = MoleculeAst::new();
    assert_eq!(m.atoms().count(), 0);
    assert_eq!(m.bonds().count(), 0);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
    assert_eq!(m.constraints().len(), 0);
}

#[rstest]
fn test_molecule_ast_default_equals_new() {
    assert_eq!(MoleculeAst::default(), MoleculeAst::new());
}

#[rstest]
fn test_molecule_ast_from_entries() {
    let atoms = vec![
        AtomAst::from_element(Element::C),
        AtomAst::from_element(Element::O),
    ];
    let bonds = vec![(AtomId(0), AtomId(1), BondAst::from_order(1))];
    let m = MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    });
    assert_eq!(m.atoms().count(), 2);
    assert_eq!(m.bonds().count(), 1);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
    assert_eq!(m.atom(AtomId(0)).ast.element, ElementAst::Lit(Element::C));
    assert_eq!(m.atom(AtomId(1)).ast.element, ElementAst::Lit(Element::O));
    assert_eq!(m.bond(BondId(0)).ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_builder() {
    assert_eq!(MoleculeAst::builder().build(), MoleculeAst::new());
}

#[rstest]
#[case::empty(MoleculeAst::default(), true)]
#[case::ground_atom(
    mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
    true,
)]
#[case::wildcard_element(
    mol_dsl!(r#"{:atoms ["*"] :bonds []}"#),
    false,
)]
#[case::wildcard_bond(
    mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "*"]]}"#),
    false,
)]
#[case::ground_atom_with_undetermined_constraint(
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![ground_atom()],
        constraints: constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![]),
            sum: ValueAst::Undetermined,
        })),
        ..Default::default()
    }),
    true,
)]
#[case::stereo_atom_ground_coset(
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![ground_atom()],
        stereo_atoms: vec![(AtomId(0), vec![], StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))],
        constraints: Constraints::new(),
        ..Default::default()
    }),
    true,
)]
#[case::stereo_atom_undetermined_coset(
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![ground_atom()],
        stereo_atoms: vec![(AtomId(0), vec![], StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))],
        constraints: Constraints::new(),
        ..Default::default()
    }),
    false,
)]
fn test_molecule_ast_is_ground(#[case] ast: MoleculeAst, #[case] expected: bool) {
    assert_eq!(ast.is_ground(), expected);
}

#[rstest]
#[case::hub(AtomId(0), vec![(AtomId(1), BondId(0)), (AtomId(2), BondId(1))])]
#[case::leaf_o(AtomId(1), vec![(AtomId(0), BondId(0))])]
#[case::leaf_n(AtomId(2), vec![(AtomId(0), BondId(1))])]
#[case::isolated(AtomId(3), vec![])]
fn test_molecule_ast_neighbors(#[case] atom: AtomId, #[case] expected: Vec<(AtomId, BondId)>) {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::C),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(0), AtomId(2), BondAst::from_order(2)),
        ],
        ..Default::default()
    });
    let mut neighbors = ast.neighbors(atom);
    assert_eq!(neighbors.len(), expected.len());
    assert_eq!(
        neighbors.size_hint(),
        (expected.len(), Some(expected.len())),
    );
    while let Some(expected_neighbor) = expected.get(expected.len() - neighbors.len()) {
        let previous = neighbors.len();
        assert_eq!(
            neighbors.next().map(|n| (n.atom_id(), n.bond_id())),
            Some(*expected_neighbor),
        );
        let remaining = neighbors.len();
        assert_eq!(remaining, previous - 1);
        assert_eq!(neighbors.size_hint(), (remaining, Some(remaining)));
    }
    assert_eq!(neighbors.next().map(|n| (n.atom_id(), n.bond_id())), None,);
    assert_eq!(neighbors.len(), 0);
}

#[rstest]
fn test_molecule_editor_add_aromatic_system() {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ..Default::default()
    });
    let mut b = ast.edit();
    let id = b.add_aromatic_system(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default());
    let new_ast = b.build();
    assert_eq!(id, AromaticSystemId(0));
    let new_atoms: Vec<AtomId> = new_ast
        .aromatic_system(AromaticSystemId(0))
        .atom_ids()
        .collect();
    assert_eq!(new_atoms, vec![AtomId(0), AtomId(1)]);
    assert_eq!(
        new_ast.aromatic_systems().ids().collect::<Vec<_>>(),
        vec![AromaticSystemId(0)]
    );
    assert_eq!(
        ast.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
}

#[fixture]
fn rich_molecule() -> MoleculeAst {
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(2)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        dative: vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::default(),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondAst::default(),
        )],
        noncovalent: vec![(
            AtomId(0),
            AtomId(3),
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        ..Default::default()
    })
}

#[fixture]
fn equiv_molecule_entries() -> MoleculeEntries {
    let mut carbon = AtomAst::from_element(Element::C);
    carbon.charge = ValueAst::Lit(1);

    MoleculeEntries {
        atoms: vec![
            carbon,
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(2)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        dative: vec![(
            vec![AtomId(1), AtomId(2)],
            AtomId(3),
            DativeBondAst::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::from_electrons(vec![1, 2, 0]),
        )],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondAst::from_electrons(vec![2, 1, 0]),
        )],
        noncovalent: vec![(
            AtomId(0),
            AtomId(3),
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondAst::new(StereoKind::CisTrans, 1u32),
        )],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(2)]),
            },
        )),
    }
}

#[rstest]
#[case::bond_endpoint(
    |entries: &mut MoleculeEntries| entries.bonds[0].0 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::dative_donor(
    |entries: &mut MoleculeEntries| entries.dative[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::dative_acceptor(
    |entries: &mut MoleculeEntries| entries.dative[0].1 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::aromatic_participant(
    |entries: &mut MoleculeEntries| entries.aromatic[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::multicenter_participant(
    |entries: &mut MoleculeEntries| entries.multicenter[0].0[0] = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::noncovalent_endpoint(
    |entries: &mut MoleculeEntries| entries.noncovalent[0].0 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_atom_site(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].0 = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_atom_ligand(
    |entries: &mut MoleculeEntries| entries.stereo_atoms[0].1[0].atom_id = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
#[case::stereo_bond_site(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].0 = BondId(3),
    Entity::Bond(BondId(3)),
)]
#[case::stereo_bond_ligand(
    |entries: &mut MoleculeEntries| entries.stereo_bonds[0].1[0].atom_id = AtomId(4),
    Entity::Atom(AtomId(4)),
)]
fn test_molecule_ast_try_from_entries_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] invalidate: fn(&mut MoleculeEntries),
    #[case] entity: Entity,
) {
    invalidate(&mut entries);

    assert_eq!(
        MoleculeAst::try_from_entries(entries),
        Err(MoleculeEntriesError::InvalidReference { entity }),
    );
}

#[rstest]
#[case::atom(Entity::Atom(AtomId(4)))]
#[case::bond(Entity::Bond(BondId(3)))]
#[case::dative_bond(Entity::DativeBond(DativeBondId(1)))]
#[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(1)))]
#[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(1)))]
#[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(1)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(1)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(1)))]
fn test_molecule_ast_try_from_entries_constraint_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] entity: Entity,
) {
    let constraint = match entity {
        Entity::Atom(id) => Constraint::Atom(id, AtomConstraintAst::valence(ValueAst::Lit(4))),
        Entity::Bond(id) => Constraint::Bond(id, BondConstraintAst::aromatic(false)),
        Entity::DativeBond(id) => {
            Constraint::DativeBond(id, DativeBondConstraintAst::aromatic(false))
        }
        Entity::AromaticSystem(id) => Constraint::AromaticSystem(
            id,
            AromaticSystemConstraintAst::electron_count(ValueAst::Lit(6)),
        ),
        Entity::MulticenterBond(id) => Constraint::MulticenterBond(
            id,
            MulticenterBondConstraintAst::electron_count(ValueAst::Lit(2)),
        ),
        Entity::NoncovalentBond(id) => {
            Constraint::NoncovalentBond(id, NoncovalentBondConstraintAst::intramolecular(true))
        }
        Entity::StereoAtom(id) => Constraint::StereoAtom(
            id,
            StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
        Entity::StereoBond(id) => Constraint::StereoBond(
            id,
            StereoKind::CisTrans,
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
    };
    entries.constraints = Constraint::Not(Box::new(constraint)).into();

    assert_eq!(
        MoleculeAst::try_from_entries(entries),
        Err(MoleculeEntriesError::InvalidReference { entity }),
    );
}

#[rstest]
#[case::relational_atom(
    Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }),
    Entity::Atom(AtomId(4)),
)]
#[case::relational_bond(
    Constraint::Relational(RelationalConstraint::DativeBondParallels {
        dative: DativeBondId(0),
        parallel: BondId(3),
    }),
    Entity::Bond(BondId(3)),
)]
#[case::atom_subset(
    Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(4)]),
    }),
    Entity::Atom(AtomId(4)),
)]
#[case::bond_subset(
    Constraint::Molecule(MoleculeConstraint::BondOrderSum {
        bonds: Some(vec![BondId(3)]),
        sum: ValueAst::Lit(1),
    }),
    Entity::Bond(BondId(3)),
)]
fn test_molecule_ast_try_from_entries_molecule_constraint_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] constraint: Constraint,
    #[case] entity: Entity,
) {
    entries.constraints = constraint.into();

    assert_eq!(
        MoleculeAst::try_from_entries(entries),
        Err(MoleculeEntriesError::InvalidReference { entity }),
    );
}

#[rstest]
#[case::target(
    {
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomId(4), AtomId(0));
        anchor
    },
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::default()],
        ..Default::default()
    }),
    Entity::Atom(AtomId(4)),
)]
#[case::pattern(
    {
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomId(0), AtomId(1));
        anchor
    },
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::default()],
        ..Default::default()
    }),
    Entity::Atom(AtomId(1)),
)]
fn test_molecule_ast_try_from_entries_subpattern_error(
    #[from(equiv_molecule_entries)] mut entries: MoleculeEntries,
    #[case] anchor: SubPatternAnchor,
    #[case] pattern: MoleculeAst,
    #[case] entity: Entity,
) {
    entries.constraints = Constraint::Molecule(MoleculeConstraint::SubPattern {
        anchor,
        pattern: Box::new(pattern),
    })
    .into();

    assert_eq!(
        MoleculeAst::try_from_entries(entries),
        Err(MoleculeEntriesError::InvalidReference { entity }),
    );
}

#[rstest]
#[should_panic(
    expected = "invalid molecule entries: molecule entries reference unavailable atom 1"
)]
fn test_molecule_ast_from_entries_error() {
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::default()],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::default())],
        ..Default::default()
    });
}

#[fixture]
fn equiv_under_molecules(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) -> (MoleculeAst, MoleculeAst, MoleculeCorrespondence) {
    let atom_images = [AtomId(2), AtomId(3), AtomId(0), AtomId(1)];
    let map_atom = |id: AtomId| atom_images[id.index()];

    let mut right_atoms = vec![AtomAst::default(); entries.atoms.len()];
    for (index, atom) in entries.atoms.iter().cloned().enumerate() {
        right_atoms[map_atom(AtomId(index as u32)).index()] = atom;
    }
    let right_bonds = entries
        .bonds
        .iter()
        .cloned()
        .map(|(first, second, ast)| (map_atom(first), map_atom(second), ast))
        .collect();
    let right_dative = entries
        .dative
        .iter()
        .cloned()
        .map(|(donors, acceptor, ast)| {
            (
                donors.into_iter().map(map_atom).collect(),
                map_atom(acceptor),
                ast,
            )
        })
        .collect();
    let right_aromatic = entries
        .aromatic
        .iter()
        .cloned()
        .map(|(atoms, ast)| (atoms.into_iter().map(map_atom).collect(), ast))
        .collect();
    let right_multicenter = entries
        .multicenter
        .iter()
        .cloned()
        .map(|(atoms, ast)| (atoms.into_iter().map(map_atom).collect(), ast))
        .collect();
    let right_noncovalent = entries
        .noncovalent
        .iter()
        .cloned()
        .map(|(first, second, ast)| (map_atom(first), map_atom(second), ast))
        .collect();
    let right_stereo_atoms = entries
        .stereo_atoms
        .iter()
        .cloned()
        .map(|(site, ligands, ast)| {
            (
                map_atom(site),
                ligands
                    .into_iter()
                    .map(|ligand| StereoLigand::new(map_atom(ligand.atom_id), ligand.kind))
                    .collect(),
                ast,
            )
        })
        .collect();
    let right_stereo_bonds = entries
        .stereo_bonds
        .iter()
        .cloned()
        .map(|(site, ligands, ast)| {
            (
                site,
                ligands
                    .into_iter()
                    .map(|ligand| StereoLigand::new(map_atom(ligand.atom_id), ligand.kind))
                    .collect(),
                ast,
            )
        })
        .collect();

    let left = MoleculeAst::from_entries(entries);
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: right_atoms,
        bonds: right_bonds,
        dative: right_dative,
        aromatic: right_aromatic,
        multicenter: right_multicenter,
        noncovalent: right_noncovalent,
        stereo_atoms: right_stereo_atoms,
        stereo_bonds: right_stereo_bonds,
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(2), AtomId(0)]),
            },
        )),
    });
    let atom_correspondence = Correspondence::from_images(&atom_images, atom_images.len());
    let correspondence = MoleculeCorrespondence::induce(&left, &right, atom_correspondence);

    (left, right, correspondence)
}

#[rstest]
fn test_molecule_ast_equiv_entity_data(#[from(equiv_molecule_entries)] entries: MoleculeEntries) {
    let base = MoleculeAst::from_entries(entries.clone());

    let mut canonical_encoding = entries.clone();
    canonical_encoding.atoms[0].charge = ValueAst::lit_set([1]);
    let canonical_encoding = MoleculeAst::from_entries(canonical_encoding);
    assert_ne!(base, canonical_encoding);
    assert!(base.equiv(&canonical_encoding));

    let mut differences = Vec::new();

    let mut atom = entries.clone();
    atom.atoms[0].charge = ValueAst::Lit(2);
    differences.push(MoleculeAst::from_entries(atom));

    let mut bond = entries.clone();
    bond.bonds[0].2.order = ValueAst::Lit(2);
    differences.push(MoleculeAst::from_entries(bond));

    let mut dative = entries.clone();
    dative.dative[0].2.order = ValueAst::Lit(2);
    differences.push(MoleculeAst::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic[0].1.electrons = ElectronCountsAst::Lit(vec![2, 0, 1]);
    differences.push(MoleculeAst::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter[0].1.electrons = ElectronCountsAst::Lit(vec![2, 0, 0]);
    differences.push(MoleculeAst::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent[0].2.kind = NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    differences.push(MoleculeAst::from_entries(noncovalent));

    let mut stereo_atom = entries.clone();
    stereo_atom.stereo_atoms[0].2 = StereoAtomAst::new(StereoKind::Tetrahedral, 0u32);
    differences.push(MoleculeAst::from_entries(stereo_atom));

    let mut stereo_bond = entries.clone();
    stereo_bond.stereo_bonds[0].2 = StereoBondAst::new(StereoKind::CisTrans, 0u32);
    differences.push(MoleculeAst::from_entries(stereo_bond));

    let mut constraint = entries;
    constraint.constraints =
        constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2)]),
        }));
    differences.push(MoleculeAst::from_entries(constraint));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_ast_equiv_relation_frames(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = MoleculeAst::from_entries(entries.clone());
    let mut differences = Vec::new();

    let mut dative = entries.clone();
    dative.dative[0].0 = vec![AtomId(0), AtomId(2)];
    differences.push(MoleculeAst::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic[0].0 = vec![AtomId(0), AtomId(1), AtomId(3)];
    differences.push(MoleculeAst::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter[0].0 = vec![AtomId(0), AtomId(1), AtomId(3)];
    differences.push(MoleculeAst::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent[0].1 = AtomId(2);
    differences.push(MoleculeAst::from_entries(noncovalent));

    let mut stereo_atom_site = entries.clone();
    stereo_atom_site.stereo_atoms[0].0 = AtomId(2);
    differences.push(MoleculeAst::from_entries(stereo_atom_site));

    let mut stereo_atom_ligand = entries.clone();
    stereo_atom_ligand.stereo_atoms[0].1[2] = StereoLigand::new(AtomId(1), StereoLigandKind::Atom);
    differences.push(MoleculeAst::from_entries(stereo_atom_ligand));

    let mut stereo_bond_site = entries.clone();
    stereo_bond_site.stereo_bonds[0].0 = BondId(2);
    differences.push(MoleculeAst::from_entries(stereo_bond_site));

    let mut stereo_bond_ligand = entries;
    stereo_bond_ligand.stereo_bonds[0].1[1] = StereoLigand::new(AtomId(2), StereoLigandKind::Atom);
    differences.push(MoleculeAst::from_entries(stereo_bond_ligand));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 8],
    );
}

#[rstest]
fn test_molecule_ast_equiv_structure_and_counts(
    #[from(equiv_molecule_entries)] entries: MoleculeEntries,
) {
    let base = MoleculeAst::from_entries(entries.clone());
    let mut differences = Vec::new();

    let mut topology = entries.clone();
    topology.bonds[2].1 = AtomId(1);
    differences.push(MoleculeAst::from_entries(topology));

    let mut atoms = entries.clone();
    atoms.atoms.push(AtomAst::from_element(Element::F));
    differences.push(MoleculeAst::from_entries(atoms));

    let mut bonds = entries.clone();
    bonds
        .bonds
        .push((AtomId(0), AtomId(3), BondAst::from_order(1)));
    differences.push(MoleculeAst::from_entries(bonds));

    let mut dative = entries.clone();
    dative.dative.pop();
    differences.push(MoleculeAst::from_entries(dative));

    let mut aromatic = entries.clone();
    aromatic.aromatic.pop();
    differences.push(MoleculeAst::from_entries(aromatic));

    let mut multicenter = entries.clone();
    multicenter.multicenter.pop();
    differences.push(MoleculeAst::from_entries(multicenter));

    let mut noncovalent = entries.clone();
    noncovalent.noncovalent.pop();
    differences.push(MoleculeAst::from_entries(noncovalent));

    let mut stereo_atom = entries.clone();
    stereo_atom.stereo_atoms.pop();
    differences.push(MoleculeAst::from_entries(stereo_atom));

    let mut stereo_bond = entries;
    stereo_bond.stereo_bonds.pop();
    differences.push(MoleculeAst::from_entries(stereo_bond));

    assert_eq!(
        differences
            .iter()
            .map(|other| base.equiv(other))
            .collect::<Vec<_>>(),
        vec![false; 9],
    );
}

#[rstest]
fn test_molecule_ast_equiv_under_non_identity(
    #[from(equiv_under_molecules)] case: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;

    assert!(correspondence.is_total());
    assert!(!left.equiv(&right));
    assert!(left.equiv_under(&right, &correspondence));
    assert!(right.equiv_under(&left, &correspondence.reverse()));
}

#[rstest]
fn test_molecule_ast_equiv_under_rejects_partial_correspondence(
    #[from(equiv_under_molecules)] case: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;
    let partial = MoleculeCorrespondence::new(
        Correspondence::new(
            vec![
                (AtomId(0), AtomId(2)),
                (AtomId(1), AtomId(3)),
                (AtomId(2), AtomId(0)),
            ],
            4,
            4,
        )
        .expect("correspondence producer preserves partial-bijection invariants"),
        correspondence.bonds().clone(),
        correspondence.dative_bonds().clone(),
        correspondence.aromatic_systems().clone(),
        correspondence.multicenter_bonds().clone(),
        correspondence.noncovalent_bonds().clone(),
        correspondence.stereo_atoms().clone(),
        correspondence.stereo_bonds().clone(),
    );

    assert!(!left.equiv_under(&right, &partial));
}

#[rstest]
fn test_molecule_ast_equiv_under_rejects_inconsistent_correspondence(
    #[from(equiv_under_molecules)] case: (MoleculeAst, MoleculeAst, MoleculeCorrespondence),
) {
    let (left, right, correspondence) = case;
    let inconsistent = MoleculeCorrespondence::new(
        correspondence.atoms().clone(),
        Correspondence::from_images(&[BondId(1), BondId(0), BondId(2)], 3),
        correspondence.dative_bonds().clone(),
        correspondence.aromatic_systems().clone(),
        correspondence.multicenter_bonds().clone(),
        correspondence.noncovalent_bonds().clone(),
        correspondence.stereo_atoms().clone(),
        correspondence.stereo_bonds().clone(),
    );

    assert!(inconsistent.is_total());
    assert!(!left.equiv_under(&right, &inconsistent));
}

#[rstest]
#[case::c_c(BondId(0), AtomId(0), AtomId(1), ValueAst::Lit(1))]
#[case::c_n(BondId(1), AtomId(1), AtomId(2), ValueAst::Lit(2))]
#[case::n_o(BondId(2), AtomId(2), AtomId(3), ValueAst::Lit(1))]
fn test_molecule_ast_bond(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] id: BondId,
    #[case] first: AtomId,
    #[case] second: AtomId,
    #[case] order: ValueAst,
) {
    let bv = ast.bond(id);
    assert_eq!(bv.id, id);
    assert_eq!(bv.atom_ids()[0], first);
    assert_eq!(bv.atom_ids()[1], second);
    assert_eq!(bv.ast.order, order);
}

#[rstest]
fn test_molecule_ast_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(BondId, AtomId, AtomId, ValueAst)> = ast
        .bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (BondId(0), AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (BondId(1), AtomId(1), AtomId(2), ValueAst::Lit(2)),
            (BondId(2), AtomId(2), AtomId(3), ValueAst::Lit(1)),
        ]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let dv = ast.dative_bond(DativeBondId(0));
    assert_eq!(dv.id, DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(3));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(2)]);
    assert_eq!(
        dv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(2), AtomId(3)]
    );
    assert_eq!(dv.ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(DativeBondId, Vec<AtomId>, AtomId)> = ast
        .dative_bonds()
        .iter()
        .map(|v| (v.id, v.donor_ids().collect(), v.acceptor_id()))
        .collect();
    assert_eq!(
        projected,
        vec![(DativeBondId(0), vec![AtomId(2)], AtomId(3))]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    let av = ast.aromatic_system(AromaticSystemId(0));
    assert_eq!(av.id, AromaticSystemId(0));
    assert_eq!(
        av.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
    assert_eq!(
        av.bond_ids().collect::<Vec<_>>(),
        vec![BondId(0), BondId(1)]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(AromaticSystemId, Vec<AtomId>, Vec<BondId>)> = ast
        .aromatic_systems()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect(), v.bond_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(
            AromaticSystemId(0),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![BondId(0), BondId(1)],
        )]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mv = ast.multicenter_bond(MulticenterBondId(0));
    assert_eq!(mv.id, MulticenterBondId(0));
    assert_eq!(
        mv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(MulticenterBondId, Vec<AtomId>)> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)],)]
    );
}

#[rstest]
fn test_molecule_ast_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let nv = ast.noncovalent_bond(NoncovalentBondId(0));
    assert_eq!(nv.id, NoncovalentBondId(0));
    assert_eq!(nv.atom_ids(), [AtomId(0), AtomId(3)]);
}

#[rstest]
fn test_molecule_ast_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(NoncovalentBondId, [AtomId; 2])> = ast
        .noncovalent_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids()))
        .collect();
    assert_eq!(
        projected,
        vec![(NoncovalentBondId(0), [AtomId(0), AtomId(3)])]
    );
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_of_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(ast.bonds().of_id(a, b), expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
#[case::wrong_donor(AtomId(3), vec![AtomId(1)], None)]
fn test_dative_bond_views_of_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(ast.dative_bonds().of_id(acceptor, &donors), expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_of_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(ast.noncovalent_bonds().of_id(a, b), expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let inc: Vec<_> = ast.dative_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let inc: Vec<_> = ast.aromatic_systems().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let inc: Vec<_> = ast.multicenter_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let inc: Vec<_> = ast.noncovalent_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
#[case::disjoint(vec![AtomId(0), AtomId(1)], vec![])]
fn test_dative_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    assert_eq!(ast.dative_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_aromatic_system_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    assert_eq!(ast.aromatic_systems().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_multicenter_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    assert_eq!(ast.multicenter_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(1), AtomId(2)], vec![])]
fn test_noncovalent_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    assert_eq!(ast.noncovalent_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_of(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(ast.bonds().of(a, b).map(|v| v.id), expected);
}

#[rstest]
#[case::pair(vec![AtomId(0), AtomId(1)], vec![BondId(0)])]
#[case::triangle(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1)])]
#[case::singleton(vec![AtomId(0)], vec![])]
fn test_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<BondId>,
) {
    let mut got: Vec<BondId> = ast
        .bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    got.sort_unstable();
    assert_eq!(got, expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = ast.dative_bonds().incident(atom).map(|v| v.id).collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
fn test_dative_bond_views_of(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(
        ast.dative_bonds().of(acceptor, &donors).map(|v| v.id),
        expected
    );
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
fn test_dative_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = ast
        .dative_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = ast
        .aromatic_systems()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_aromatic_system_views_of(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(ast.aromatic_systems().of(atoms).map(|v| v.id), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_aromatic_system_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = ast
        .aromatic_systems()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = ast
        .multicenter_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_multicenter_bond_views_of(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(ast.multicenter_bonds().of(atoms).map(|v| v.id), expected,);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_multicenter_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = ast
        .multicenter_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = ast
        .noncovalent_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_of(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(ast.noncovalent_bonds().of(a, b).map(|v| v.id), expected,);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_noncovalent_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = ast
        .noncovalent_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::atom_0(AtomId(0), Element::C)]
#[case::atom_1(AtomId(1), Element::C)]
#[case::atom_2(AtomId(2), Element::N)]
#[case::atom_3(AtomId(3), Element::O)]
fn test_molecule_ast_atom(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] id: AtomId,
    #[case] element: Element,
) {
    let av = ast.atom(id);
    assert_eq!(av.id, id);
    assert_eq!(av.ast.element, ElementAst::Lit(element));
}

#[rstest]
fn test_molecule_ast_is_empty() {
    assert!(MoleculeAst::default().is_empty());
}

#[rstest]
fn test_molecule_ast_is_empty_rich(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(!ast.is_empty());
}

#[rstest]
fn test_molecule_ast_has_constraints_empty() {
    assert!(!MoleculeAst::default().has_constraints());
}

#[rstest]
fn test_molecule_ast_has_constraints_rich(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(!ast.has_constraints());
}

#[rstest]
fn test_molecule_ast_has_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_dative_bonds());
}

#[rstest]
fn test_molecule_ast_has_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_aromatic_systems());
}

#[rstest]
fn test_molecule_ast_has_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_multicenter_bonds());
}

#[rstest]
fn test_molecule_ast_has_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_noncovalent_bonds());
}

#[rstest]
fn test_molecule_ast_has_stereo_atoms() {
    let ast = mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#
    );
    assert!(ast.has_stereo_atoms());
    assert!(!ast.has_stereo_bonds());
}

#[rstest]
fn test_molecule_ast_has_stereo_bonds() {
    let ast = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#
    );
    assert!(ast.has_stereo_bonds());
    assert!(!ast.has_stereo_atoms());
}

#[rstest]
fn test_molecule_ast_atoms(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(AtomId, ElementAst)> = ast
        .atoms()
        .iter()
        .map(|v| (v.id, v.ast.element.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (AtomId(0), ElementAst::Lit(Element::C)),
            (AtomId(1), ElementAst::Lit(Element::C)),
            (AtomId(2), ElementAst::Lit(Element::N)),
            (AtomId(3), ElementAst::Lit(Element::O)),
        ]
    );
}

#[test]
fn test_bond_views_induced_ids() {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
        ],
        ..Default::default()
    });
    let bonds = ast.bonds().induced_ids(&[AtomId(0), AtomId(1)]);
    assert_eq!(bonds, vec![BondId(0)]);

    let mut all = ast.bonds().induced_ids(&[AtomId(0), AtomId(1), AtomId(2)]);
    all.sort_unstable();
    assert_eq!(all, vec![BondId(0), BondId(1), BondId(2)]);
}

fn chain(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n.saturating_sub(1))
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId((i + 1) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn ring(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n)
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId(((i + 1) % n) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn two_components() -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); 4];
    let bonds = vec![
        (AtomId(0), AtomId(1), BondAst::from_order(1)),
        (AtomId(2), AtomId(3), BondAst::from_order(1)),
    ];
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

#[rstest]
#[case::isolated(chain(1), AtomId(0), 0)]
#[case::chain_end(chain(3), AtomId(0), 1)]
#[case::chain_mid(chain(3), AtomId(1), 2)]
#[case::ring_vertex(ring(6), AtomId(0), 2)]
fn test_molecule_ast_degree(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: usize,
) {
    assert_eq!(ast.graph().degree(atom), expected);
}

#[rstest]
#[case::single(chain(3), 1)]
#[case::two(two_components(), 2)]
#[case::empty(MoleculeAst::default(), 0)]
fn test_molecule_ast_connected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let cc = ast
        .graph()
        .connected_components(ConnectedComponentsAlgorithm::Bfs);
    assert_eq!(cc.len(), expected);
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain(chain(5), 0)]
fn test_molecule_ast_biconnected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let bcc = ast
        .graph()
        .biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
    assert_eq!(bcc.len(), expected);
}

#[rstest]
#[case::ring_bond(ring(6), BondId(0), Some(6))]
#[case::chain_bond(chain(3), BondId(0), None)]
fn test_molecule_ast_shortest_cycle_through_bond(
    #[case] ast: MoleculeAst,
    #[case] bond: BondId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.graph()
            .shortest_cycle_through_bond(bond, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::ring_atom(ring(6), AtomId(0), Some(6))]
#[case::chain_atom(chain(3), AtomId(1), None)]
fn test_molecule_ast_shortest_cycle_through_atom(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.graph()
            .shortest_cycle_through_atom(atom, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
#[case::chain(chain(5), 10, vec![])]
#[case::empty(MoleculeAst::default(), 10, vec![])]
fn test_graph_view_enumerate_simple_cycles(
    #[case] ast: MoleculeAst,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    assert_eq!(
        ast.graph()
            .enumerate_simple_cycles(max_size, SimpleCycleEnumerationAlgorithm::ReadTarjan),
        expected
    );
}

#[rstest]
#[case::hexagon(
    ring(6),
    6,
    vec![vec![
        AtomId(0),
        AtomId(1),
        AtomId(2),
        AtomId(3),
        AtomId(4),
        AtomId(5),
    ]],
)]
#[case::hexagon_cutoff(ring(6), 5, vec![])]
#[case::chain(chain(5), 10, vec![])]
#[case::empty(MoleculeAst::default(), 10, vec![])]
fn test_graph_view_enumerate_relevant_cycles(
    #[case] ast: MoleculeAst,
    #[case] max_size: usize,
    #[case] expected: Vec<Vec<AtomId>>,
) {
    assert_eq!(
        ast.graph()
            .enumerate_relevant_cycles(max_size, RelevantCycleEnumerationAlgorithm::Vismara),
        expected
    );
}

#[rstest]
#[case::triangle(ring(3), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_maximum_independent_set(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let mis = ast
        .graph()
        .maximum_independent_set(MaximumIndependentSetAlgorithm::BranchAndBound);
    assert_eq!(mis.len(), expected);
}

#[rstest]
#[case::ring_6(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 4 :single] [4 5 :single] [5 0 :single]]}"#),
    vec![BondId(0), BondId(2), BondId(4)],
)]
fn test_graph_view_bipartite_maximum_matching(
    #[case] ast: MoleculeAst,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = ast.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        ast.graph()
            .bipartite_maximum_matching(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
            )
            .expect("case graph is bipartite")
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::triangle(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#))]
fn test_graph_view_bipartite_maximum_matching_error(#[case] ast: MoleculeAst) {
    let node_order: Vec<AtomId> = ast.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        ast.graph()
            .bipartite_maximum_matching(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
            )
            .unwrap_err(),
        NonBipartiteGraphError,
    );
}

#[rstest]
#[case::chain_4(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single]]}"#),
    vec![BondId(0), BondId(2)],
)]
#[case::triangle(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#),
    vec![BondId(0)],
)]
#[case::single(mol_dsl!(r#"{:atoms ["C"]}"#), vec![])]
fn test_graph_view_general_maximum_matching(
    #[case] ast: MoleculeAst,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = ast.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        ast.graph()
            .general_maximum_matching(&node_order, GeneralMaximumMatchingAlgorithm::Edmonds)
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::bipartite(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [2 3 :single] [3 0 :single]]}"#),
    vec![BondId(0), BondId(2)],
)]
#[case::general(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 :single] [1 2 :single] [0 2 :single]]}"#),
    vec![BondId(0)],
)]
fn test_graph_view_bipartite_maximum_matching_or_general(
    #[case] ast: MoleculeAst,
    #[case] expected: Vec<BondId>,
) {
    let node_order: Vec<AtomId> = ast.atoms().iter().map(|atom| atom.id).collect();
    assert_eq!(
        ast.graph()
            .bipartite_maximum_matching_or_general(
                &node_order,
                BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                GeneralMaximumMatchingAlgorithm::Edmonds,
            )
            .bonds()
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::ring_6(ring(6), 2)]
fn test_molecule_ast_enumerate_perfect_matchings(
    #[case] ast: MoleculeAst,
    #[case] expected: usize,
) {
    let ms = ast
        .graph()
        .enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound);
    assert_eq!(ms.len(), expected);
    for m in &ms {
        assert!(m.is_perfect(ast.atoms().count()));
    }
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_automorphisms(#[case] ast: MoleculeAst, #[case] expected_orbits: usize) {
    let auto = ast
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert_eq!(auto.orbit_count(), expected_orbits);
    assert_eq!(auto.atom_count(), ast.atoms().count());
}

#[test]
fn test_atom_automorphism_same_orbit() {
    let ast = ring(6);
    let auto = ast
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert!(auto.same_orbit(AtomId(0), AtomId(3)));
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().subgraph_isomorphisms(
        &query.graph(),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(0), AtomId(5)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
            vec![AtomId(3), AtomId(2)],
            vec![AtomId(3), AtomId(4)],
            vec![AtomId(4), AtomId(3)],
            vec![AtomId(4), AtomId(5)],
            vec![AtomId(5), AtomId(0)],
            vec![AtomId(5), AtomId(4)],
        ]
    );
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().subgraph_isomorphisms_at(
        &query.graph(),
        (AtomId(0), AtomId(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(5)],]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    let extracted = ast.extract(&sub);
    let atom_elements: Vec<_> = extracted
        .atoms()
        .iter()
        .map(|v| v.ast.element.clone())
        .collect();
    assert_eq!(
        atom_elements,
        vec![
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::N),
        ]
    );
    let bonds: Vec<(AtomId, AtomId, ValueAst)> = extracted
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (AtomId(1), AtomId(2), ValueAst::Lit(2)),
        ]
    );
    assert_eq!(
        sub.atoms().matched_pairs(),
        &[
            (AtomId(0), AtomId(0)),
            (AtomId(1), AtomId(1)),
            (AtomId(2), AtomId(2))
        ]
    );
    assert_eq!(
        sub.bonds().matched_pairs(),
        &[(BondId(0), BondId(0)), (BondId(1), BondId(1))]
    );
    assert_eq!(
        sub.aromatic_systems().matched_pairs(),
        &[(AromaticSystemId(0), AromaticSystemId(0))]
    );
    assert_eq!(
        sub.multicenter_bonds().matched_pairs(),
        &[(MulticenterBondId(0), MulticenterBondId(0))]
    );
    assert_eq!(
        sub.dative_bonds().matched_pairs(),
        &[] as &[(DativeBondId, DativeBondId)]
    );
    assert_eq!(
        sub.noncovalent_bonds().matched_pairs(),
        &[] as &[(NoncovalentBondId, NoncovalentBondId)]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph_preserves_dative(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(2), AtomId(3)]);
    assert_eq!(
        sub.atoms().matched_pairs(),
        &[(AtomId(0), AtomId(2)), (AtomId(1), AtomId(3))]
    );
    assert_eq!(
        sub.dative_bonds().matched_pairs(),
        &[(DativeBondId(0), DativeBondId(0))]
    );
    let extracted = ast.extract(&sub);
    let dv = extracted.dative_bond(DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(1));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(0)]);
    assert_eq!(dv.ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_edits(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    assert_eq!(
        ast.edits(&sub),
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(3))],
            bonds: vec![BondHandle::Id(BondId(2))],
        }])
    );
}

#[rstest]
fn test_molecule_ast_edits_identity(#[from(rich_molecule)] ast: MoleculeAst) {
    let atom_ids: Vec<AtomId> = ast.atoms().iter().map(|v| v.id).collect();
    let sub = ast.induced_subgraph(&atom_ids);
    assert_eq!(ast.edits(&sub), Edits::new());
}

#[rstest]
#[case::add_atom(
    mol_dsl!(r#"{:atoms ["C"]}"#),
    Edits::from_iter([Edit::AddAtoms {
        atoms: vec![AtomAst::from_element(Element::N)],
    }]),
    mol_dsl!(r#"{:atoms ["C" "N"]}"#),
)]
fn test_molecule_ast_apply(
    #[case] molecule: MoleculeAst,
    #[case] edits: Edits,
    #[case] expected: MoleculeAst,
) {
    let original = molecule.clone();

    assert_eq!(molecule.apply(edits), Ok(expected));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::stale_after_add(
    mol_dsl!(r#"{:atoms ["C"]}"#),
    Edits::from_iter([
        Edit::AddAtoms {
            atoms: vec![AtomAst::from_element(Element::N)],
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        },
    ]),
    TransactionError::OldStateMismatch,
)]
fn test_molecule_ast_apply_error(
    #[case] molecule: MoleculeAst,
    #[case] edits: Edits,
    #[case] expected: TransactionError,
) {
    let original = molecule.clone();

    assert_eq!(molecule.apply(edits), Err(expected));
    assert_eq!(molecule, original);
}

#[rstest]
fn test_molecule_ast_extract(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1)]);
    let extracted = ast.extract(&sub);
    assert_eq!(extracted.atoms().count(), 2);
}

#[rstest]
fn test_molecule_editor_remove_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    let result = b.build();
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.atoms().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]
    );
    assert_eq!(
        result.bonds().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![BondId(0), BondId(1), BondId(2)]
    );
}

#[rstest]
fn test_molecule_editor_remove_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_dative_bonds(&[DativeBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_remove_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_multicenter_bonds(&[MulticenterBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bonds().ids().collect::<Vec<_>>(),
        Vec::<MulticenterBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_remove_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
fn test_molecule_editor_atom_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomId(0)).ast.element = ElementAst::Lit(Element::N);
    let result = b.build();
    assert_eq!(
        result.atom(AtomId(0)).ast.element,
        ElementAst::Lit(Element::N)
    );
    assert_eq!(ast.atom(AtomId(0)).ast.element, ElementAst::Lit(Element::C));
}

#[rstest]
fn test_molecule_editor_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.bond_mut(BondId(0)).ast.order = ValueAst::Lit(3);
    let result = b.build();
    assert_eq!(result.bond(BondId(0)).ast.order, ValueAst::Lit(3));
    assert_eq!(ast.bond(BondId(0)).ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_editor_atom_constraint_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomId(0))
        .ast
        .constraints
        .set(AtomConstraintAst::Degree(ValueAst::Lit(2)));
    let result = b.build();
    assert_eq!(
        result.atom(AtomId(0)).ast.constraints,
        AtomConstraintsAst::from_iter([AtomConstraintAst::Degree(ValueAst::Lit(2))])
    );
    assert!(ast.atom(AtomId(0)).ast.constraints.is_empty());
}

#[rstest]
fn test_molecule_editor_add_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_dative_bond(vec![AtomId(1)], AtomId(0), DativeBondAst::from_order(1));
    let result = b.build();
    assert_eq!(id, DativeBondId(1));
    let view = result.dative_bond(id);
    assert_eq!(view.acceptor_id(), AtomId(0));
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![AtomId(1)]);
}

#[rstest]
fn test_molecule_editor_add_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_multicenter_bond(
        vec![AtomId(1), AtomId(2), AtomId(3)],
        MulticenterBondAst::default(),
    );
    let result = b.build();
    assert_eq!(id, MulticenterBondId(1));
    let atoms: Vec<AtomId> = result.multicenter_bond(id).atom_ids().collect();
    assert_eq!(atoms, vec![AtomId(1), AtomId(2), AtomId(3)]);
}

#[rstest]
fn test_molecule_editor_add_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_noncovalent_bond(
        [AtomId(1), AtomId(2)],
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let result = b.build();
    assert_eq!(id, NoncovalentBondId(1));
    let view = result.noncovalent_bond(id);
    assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
}

#[rstest]
fn test_molecule_editor_push_constraint_and_constraints_mut(
    #[from(rich_molecule)] ast: MoleculeAst,
) {
    let mut b = ast.edit();
    b.push_constraint(Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    }));
    b.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: ValueAst::Lit(0),
        }));
    let result = b.build();
    assert_eq!(result.constraints().len(), 2);
}

#[rstest]
fn test_molecule_editor_dative_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.dative_bond_mut(DativeBondId(0)).ast.constraints.set(
        DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1),
    );
    let result = b.build();
    assert!(!result
        .dative_bond(DativeBondId(0))
        .ast
        .constraints
        .is_empty());
    assert!(ast.dative_bond(DativeBondId(0)).ast.constraints.is_empty());
}

#[rstest]
fn test_molecule_editor_aromatic_system_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.aromatic_system_mut(AromaticSystemId(0)).ast.charge = ValueAst::Lit(0);
    let result = b.build();
    assert_eq!(
        result.aromatic_system(AromaticSystemId(0)).ast.charge,
        ValueAst::Lit(0)
    );
}

#[rstest]
fn test_molecule_editor_multicenter_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.multicenter_bond_mut(MulticenterBondId(0)).ast.electrons =
        ElectronCountsAst::Lit(vec![1, 1, 0]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bond(MulticenterBondId(0)).ast.electrons,
        ElectronCountsAst::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_editor_noncovalent_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.noncovalent_bond_mut(NoncovalentBondId(0)).ast.kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bond(NoncovalentBondId(0)).ast.kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
    );
}

#[rstest]
fn test_molecule_editor_remove_empty_is_noop(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_dative_bonds(&[]);
    b.remove_aromatic_systems(&[]);
    b.remove_multicenter_bonds(&[]);
    b.remove_noncovalent_bonds(&[]);
    let result = b.build();
    assert_eq!(result.dative_bonds().count(), 1);
    assert_eq!(result.aromatic_systems().count(), 1);
    assert_eq!(result.multicenter_bonds().count(), 1);
    assert_eq!(result.noncovalent_bonds().count(), 1);
}

#[rstest]
#[case::hexagon(
    ring(6),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 6,
    },
    vec![(
        vec![
            AtomId(0),
            AtomId(1),
            AtomId(2),
            AtomId(3),
            AtomId(4),
            AtomId(5),
        ],
        vec![
            BondId(0),
            BondId(1),
            BondId(2),
            BondId(3),
            BondId(4),
            BondId(5),
        ],
    )],
)]
#[case::hexagon_cutoff(
    ring(6),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 5,
    },
    vec![],
)]
#[case::chain(
    chain(5),
    RingModel {
        kind: RingSetKind::Relevant,
        max_ring_size: 10,
    },
    vec![],
)]
#[case::empty(MoleculeAst::default(), RingModel::default(), vec![])]
fn test_molecule_ast_rings(
    #[case] ast: MoleculeAst,
    #[case] model: RingModel,
    #[case] expected: Vec<(Vec<AtomId>, Vec<BondId>)>,
) {
    let rings = ast
        .rings(model, RingConfig::default())
        .iter()
        .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
        .collect::<Vec<_>>();
    assert_eq!(rings, expected);
}

#[rstest]
#[case::simple(
    RingSetKind::Simple,
    vec![
        vec![BondId(0), BondId(3), BondId(1)],
        vec![BondId(0), BondId(3), BondId(5), BondId(2)],
        vec![BondId(0), BondId(4), BondId(2)],
        vec![BondId(0), BondId(4), BondId(5), BondId(1)],
        vec![BondId(1), BondId(3), BondId(4), BondId(2)],
        vec![BondId(1), BondId(5), BondId(2)],
        vec![BondId(3), BondId(5), BondId(4)],
    ],
)]
#[case::relevant(
    RingSetKind::Relevant,
    vec![
        vec![BondId(0), BondId(3), BondId(1)],
        vec![BondId(0), BondId(4), BondId(2)],
        vec![BondId(1), BondId(5), BondId(2)],
        vec![BondId(3), BondId(5), BondId(4)],
    ],
)]
fn test_molecule_ast_rings_kind(#[case] kind: RingSetKind, #[case] mut expected: Vec<Vec<BondId>>) {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(1)),
            (AtomId(1), AtomId(3), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        ..Default::default()
    });
    let mut actual = ast
        .rings(
            RingModel {
                kind,
                max_ring_size: 4,
            },
            RingConfig::default(),
        )
        .iter()
        .map(|ring| ring.bonds().to_vec())
        .collect::<Vec<_>>();
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected);
}

#[rstest]
fn test_molecule_ast_rings_parallel_bond_identity() {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(1)),
            (AtomId(2), AtomId(0), BondAst::from_order(1)),
        ],
        ..Default::default()
    });
    let mut actual = ast
        .rings(
            RingModel {
                kind: RingSetKind::Simple,
                max_ring_size: 3,
            },
            RingConfig::default(),
        )
        .iter()
        .map(|ring| ring.bonds().to_vec())
        .collect::<Vec<_>>();
    actual.sort_unstable();
    assert_eq!(
        actual,
        vec![
            vec![BondId(0), BondId(2), BondId(3)],
            vec![BondId(1), BondId(2), BondId(3)],
        ]
    );
}

#[rstest]
#[case::self_loop(1, vec![[0, 0]])]
#[case::parallel_pair(2, vec![[0, 1], [0, 1]])]
fn test_molecule_ast_rings_invalid(#[case] atom_count: usize, #[case] edges: Vec<[u32; 2]>) {
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C); atom_count],
        bonds: edges
            .into_iter()
            .map(|[first, second]| (AtomId(first), AtomId(second), BondAst::from_order(1)))
            .collect(),
        ..Default::default()
    });
    let rings = ast.rings(
        RingModel {
            kind: RingSetKind::Simple,
            max_ring_size: usize::MAX,
        },
        RingConfig::default(),
    );
    assert_eq!(
        rings
            .iter()
            .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
            .collect::<Vec<_>>(),
        vec![]
    );
}

#[rstest]
fn test_molecule_editor_add_and_remove(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let new_a = b.add_atom(AtomAst::from_element(Element::Br));
    b.add_bond(AtomId(0), new_a, BondAst::from_order(1));
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    let _compaction = b.remove(&[AtomId(3)], &[BondId(2)]);
    let result = b.build();
    let atoms: Vec<Element> = result
        .atoms()
        .iter()
        .map(|v| match v.ast.element {
            ElementAst::Lit(e) => e,
            _ => panic!("non-ground element in editor result"),
        })
        .collect();
    assert_eq!(atoms, vec![Element::C, Element::C, Element::N, Element::Br]);
    let bonds: Vec<(AtomId, AtomId, ValueAst)> = result
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (AtomId(1), AtomId(2), ValueAst::Lit(2)),
            (AtomId(0), AtomId(3), ValueAst::Lit(1)),
        ]
    );
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
#[case::donor_below_acceptor(AtomId(0), AtomId(1))]
#[case::donor_above_acceptor(AtomId(1), AtomId(0))]
fn test_molecule_ast_dative_acceptor_donor(#[case] donor: AtomId, #[case] acceptor: AtomId) {
    let atoms = vec![ground_atom(), ground_atom()];
    let ast = MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        dative: vec![(vec![donor], acceptor, DativeBondAst::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let view = ast.dative_bond(DativeBondId(0));
    assert_eq!(view.acceptor_id(), acceptor);
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![donor]);
}

#[rstest]
fn test_molecule_ast_eq_canonical_across_bond_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let bond = BondAst {
        order: ValueAst::Lit(1),
        charge: ValueAst::Lit(0),
        unpaired_electrons: UnpairedElectronsAst::closed_shell(),
        constraints: BondConstraintsAst::new(),
    };
    let forward = MoleculeAst::from_entries(MoleculeEntries {
        atoms: atoms_a,
        bonds: vec![(AtomId(0), AtomId(1), bond.clone())],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let reverse = MoleculeAst::from_entries(MoleculeEntries {
        atoms: atoms_b,
        bonds: vec![(AtomId(1), AtomId(0), bond)],
        constraints: Constraints::new(),
        ..Default::default()
    });
    assert_eq!(forward, reverse);
}

#[rstest]
fn test_molecule_ast_eq_canonical_across_dative_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let forward = MoleculeAst::from_entries(MoleculeEntries {
        atoms: atoms_a,
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let reverse = MoleculeAst::from_entries(MoleculeEntries {
        atoms: atoms_b,
        dative: vec![(vec![AtomId(1)], AtomId(0), DativeBondAst::from_order(1))],
        constraints: Constraints::new(),
        ..Default::default()
    });
    assert_ne!(
        forward, reverse,
        "acceptor identity is part of dative bond; swapping donor/acceptor should differ"
    );
}

#[rstest]
fn test_molecule_ast_raw_graph(#[from(rich_molecule)] ast: MoleculeAst) {
    let g = ast.raw_graph();
    assert_eq!(g.node_count(), 4);
    assert_eq!(g.edge_count(), 3);
    assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
#[case::disjoint(
    HashSet::from([AtomId(3)]),
    None,
)]
fn test_aromatic_system_views_of_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(ast.aromatic_systems().of_id(atoms), expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
fn test_multicenter_bond_views_of_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(ast.multicenter_bonds().of_id(atoms), expected);
}

#[rstest]
fn test_molecule_ast_enumerate_maximum_matchings() {
    let ast = ring(4);
    let mut ms: Vec<Vec<(AtomId, AtomId)>> = ast
        .graph()
        .enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound)
        .into_iter()
        .map(|m| {
            let mut pairs: Vec<_> = (0..ast.atoms().count())
                .map(AtomId::from)
                .filter_map(|a| m.mate(a).filter(|b| a < *b).map(|b| (a, b)))
                .collect();
            pairs.sort_unstable();
            pairs
        })
        .collect();
    ms.sort_unstable();
    assert_eq!(
        ms,
        vec![
            vec![(AtomId(0), AtomId(1)), (AtomId(2), AtomId(3))],
            vec![(AtomId(0), AtomId(3)), (AtomId(1), AtomId(2))],
        ]
    );
}

#[rstest]
fn test_molecule_ast_index_atom(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast.atom(AtomId(2)).ast.element, ElementAst::Lit(Element::N));
}

#[rstest]
fn test_molecule_ast_index_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast.bond(BondId(1)).ast.order, ValueAst::Lit(2));
}

#[rstest]
fn test_molecule_ast_index_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast.dative_bond(DativeBondId(0)).ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_index_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast.aromatic_system(AromaticSystemId(0)).ast.electrons,
        ElectronCountsAst::Undetermined
    );
}

#[rstest]
fn test_molecule_ast_index_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast.multicenter_bond(MulticenterBondId(0)).ast.electrons,
        ElectronCountsAst::Undetermined
    );
}

#[rstest]
fn test_molecule_ast_index_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast.noncovalent_bond(NoncovalentBondId(0)).ast.kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)
    );
}

#[rstest]
fn test_molecule_ast_modify_atoms(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.modify_atoms(|mut a| {
        a.charge = ValueAst::Lit(1);
        a
    });
    let charges: Vec<ValueAst> = ast.atoms().iter().map(|v| v.ast.charge.clone()).collect();
    assert_eq!(
        charges,
        vec![
            ValueAst::Lit(1),
            ValueAst::Lit(1),
            ValueAst::Lit(1),
            ValueAst::Lit(1),
        ]
    );
}

#[rstest]
fn test_molecule_ast_modify_bonds(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.modify_bonds(|mut b| {
        b.order = ValueAst::Lit(1);
        b
    });
    let orders: Vec<ValueAst> = ast.bonds().iter().map(|v| v.ast.order.clone()).collect();
    assert_eq!(
        orders,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(1)]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.dative_bond_mut(DativeBondId(0)).ast.constraints.set(
        DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
    );
    assert_eq!(
        ast.dative_bond(DativeBondId(0)).ast.constraints,
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(
            RingScope::Size(6),
            1
        )])
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.aromatic_system_mut(AromaticSystemId(0)).ast.electrons = ElectronCountsAst::Lit(vec![1; 3]);
    assert_eq!(
        ast.aromatic_system(AromaticSystemId(0)).ast.electrons,
        ElectronCountsAst::Lit(vec![1, 1, 1]),
    );
}

#[rstest]
fn test_molecule_ast_modify_aromatic_systems(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.modify_aromatic_systems(|mut a| {
        a.electrons = ElectronCountsAst::Lit(vec![1; 3]);
        a
    });
    let electrons: Vec<ElectronCountsAst> = ast
        .aromatic_systems()
        .iter()
        .map(|v| v.ast.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsAst::Lit(vec![1; 3])]);
}

#[rstest]
fn test_molecule_ast_multicenter_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.multicenter_bond_mut(MulticenterBondId(0)).ast.electrons =
        ElectronCountsAst::Lit(vec![1, 1, 0]);
    assert_eq!(
        ast.multicenter_bond(MulticenterBondId(0)).ast.electrons,
        ElectronCountsAst::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_ast_modify_multicenter_bonds(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.modify_multicenter_bonds(|mut m| {
        m.electrons = ElectronCountsAst::Lit(vec![1, 1, 0]);
        m
    });
    let electrons: Vec<ElectronCountsAst> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| v.ast.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsAst::Lit(vec![1, 1, 0])],);
}

#[rstest]
fn test_molecule_ast_noncovalent_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.noncovalent_bond_mut(NoncovalentBondId(0)).ast.kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    assert_eq!(
        ast.noncovalent_bond(NoncovalentBondId(0)).ast.kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)
    );
}

// -- lift_constraints / inline_constraints ---------------------

/// Set-equality assertion for constraint vecs: order is unspecified, so the
/// test compares as multisets via sort+eq.
fn assert_same_constraints(a: &Constraints, b: &Constraints) {
    let mut x: Vec<&Constraint> = a.iter().collect();
    let mut y: Vec<&Constraint> = b.iter().collect();
    x.sort_by_key(|c| format!("{c:?}"));
    y.sort_by_key(|c| format!("{c:?}"));
    assert_eq!(x, y);
}

#[rstest]
fn test_molecule_ast_lift_constraints_empty() {
    let mut ast = MoleculeAst::default();
    ast.lift_constraints();
    assert!(ast.constraints().is_empty());
}

#[rstest]
fn test_molecule_ast_lift_constraints_drains_inline_stores(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .set(AtomConstraintAst::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomId(2))
        .ast
        .constraints
        .set(AtomConstraintAst::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondId(0))
        .ast
        .constraints
        .set(BondConstraintAst::Aromatic(BooleanAst::Lit(true)));
    ast.dative_bond_mut(DativeBondId(0)).ast.constraints.set(
        DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)),
    );

    ast.lift_constraints();

    assert!(ast.atom(AtomId(0)).ast.constraints.is_empty());
    assert!(ast.atom(AtomId(2)).ast.constraints.is_empty());
    assert!(ast.bond(BondId(0)).ast.constraints.is_empty());
    assert!(ast.dative_bond(DativeBondId(0)).ast.constraints.is_empty());

    let mut expected = Constraints::new();
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));
    expected.push(Constraint::Atom(
        AtomId(2),
        AtomConstraintAst::Degree(ValueAst::Lit(3)),
    ));
    expected.push(Constraint::Bond(
        BondId(0),
        BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
    ));
    expected.push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_lift_constraints_appends_to_existing(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let prior = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    ast.constraints_mut().push(prior.clone());
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .set(AtomConstraintAst::Valence(ValueAst::Lit(4)));

    ast.lift_constraints();

    let mut expected = Constraints::new();
    expected.push(prior);
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_inline_constraints_drains_top_level_leaves(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));
    ast.constraints_mut().push(Constraint::Bond(
        BondId(0),
        BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
    ));
    ast.constraints_mut().push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1),
    ));

    ast.inline_constraints();

    assert!(ast.constraints().is_empty());
    assert_eq!(
        ast.atom(AtomId(0)).ast.constraints,
        AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::Lit(4))])
    );
    assert_eq!(
        ast.bond(BondId(0)).ast.constraints,
        BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanAst::Lit(true))])
    );
    assert_eq!(
        ast.dative_bond(DativeBondId(0)).ast.constraints,
        DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(
            RingScope::Size(5),
            1
        )])
    );
}

#[rstest]
fn test_molecule_ast_inline_constraints_last_wins_on_collision(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(3)),
    ));
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    // Only one Valence survives; with two competing inserts of the same kind,
    // exactly one wins (which one is unspecified). Verify count and kind.
    assert_eq!(ast.atom(AtomId(0)).ast.constraints.len(), 1);
    let v = ast
        .atom(AtomId(0))
        .ast
        .constraints
        .iter()
        .next()
        .unwrap()
        .clone();
    assert!(matches!(v, AtomConstraintAst::Valence(_)));
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_combinator_nested(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let leaf = Constraint::Atom(AtomId(0), AtomConstraintAst::Valence(ValueAst::Lit(4)));
    let nested = Constraint::And(vec![
        leaf.clone(),
        Constraint::Bond(
            BondId(0),
            BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
        ),
    ]);
    ast.constraints_mut().push(nested.clone());

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(nested);
    assert_same_constraints(ast.constraints(), &expected);
    assert!(ast.atom(AtomId(0)).ast.constraints.is_empty());
    assert!(ast.bond(BondId(0)).ast.constraints.is_empty());
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_relational_and_molecule(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let rel = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    let mol = Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    });
    ast.constraints_mut().push(rel.clone());
    ast.constraints_mut().push(mol.clone());
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(rel);
    expected.push(mol);
    assert_same_constraints(ast.constraints(), &expected);
    assert_eq!(
        ast.atom(AtomId(0)).ast.constraints,
        AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::Lit(4))])
    );
}

#[rstest]
fn test_molecule_ast_lift_then_inline_roundtrips_inline_state(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .set(AtomConstraintAst::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .set(AtomConstraintAst::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondId(0))
        .ast
        .constraints
        .set(BondConstraintAst::Aromatic(BooleanAst::Lit(true)));
    ast.dative_bond_mut(DativeBondId(0)).ast.constraints.set(
        DativeBondConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)),
    );

    let original = ast.clone();

    ast.lift_constraints();
    assert!(ast.atom(AtomId(0)).ast.constraints.is_empty());
    ast.inline_constraints();

    assert_eq!(ast, original);
}

#[rstest]
#[case::empty(Vec::new(), MoleculeAst::new(), Vec::new())]
#[case::singleton(
    vec![MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    })],
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    }),
    vec![vec![(AtomId(0), AtomId(0))]],
)]
#[case::multiple(
    vec![
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        }),
        MoleculeAst::new(),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::O),
                AtomAst::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        }),
    ],
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![(AtomId(1), AtomId(2), BondAst::from_order(2))],
        ..Default::default()
    }),
    vec![
        vec![(AtomId(0), AtomId(0))],
        vec![],
        vec![(AtomId(0), AtomId(1)), (AtomId(1), AtomId(2))],
    ],
)]
fn test_molecule_ast_combine_all(
    #[case] molecules: Vec<MoleculeAst>,
    #[case] expected: MoleculeAst,
    #[case] expected_atom_matched_pairs: Vec<Vec<(AtomId, AtomId)>>,
) {
    let (combined, correspondences) = MoleculeAst::combine_all(&molecules);

    assert_eq!(combined, expected);
    assert_eq!(
        correspondences
            .iter()
            .map(|correspondence| correspondence.atoms().matched_pairs().to_vec())
            .collect::<Vec<_>>(),
        expected_atom_matched_pairs,
    );
    for (molecule, correspondence) in molecules.iter().zip(&correspondences) {
        assert_eq!(combined.extract(correspondence), *molecule);
        assert_eq!(
            correspondence.atoms().right_count(),
            combined.atoms().count()
        );
        assert_eq!(
            correspondence.bonds().right_count(),
            combined.bonds().count()
        );
        assert_eq!(
            correspondence.dative_bonds().right_count(),
            combined.dative_bonds().count()
        );
        assert_eq!(
            correspondence.aromatic_systems().right_count(),
            combined.aromatic_systems().count()
        );
        assert_eq!(
            correspondence.multicenter_bonds().right_count(),
            combined.multicenter_bonds().count()
        );
        assert_eq!(
            correspondence.noncovalent_bonds().right_count(),
            combined.noncovalent_bonds().count()
        );
        assert_eq!(
            correspondence.stereo_atoms().right_count(),
            combined.stereo_atoms().count()
        );
        assert_eq!(
            correspondence.stereo_bonds().right_count(),
            combined.stereo_bonds().count()
        );
    }
}

#[rstest]
fn test_molecule_ast_combine() {
    let left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
        ..Default::default()
    });
    let (union, correspondence) = left.combine(&right);

    assert_eq!(union.atoms().count(), 4);
    assert_eq!(union.bonds().count(), 2);
    assert_eq!(union.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    assert_eq!(union.bond(BondId(1)).atom_ids(), [AtomId(2), AtomId(3)]);
    assert_eq!(union.bond(BondId(1)).ast, &BondAst::from_order(2));
    // right's ids map to their offset union ids; left's are the prefix (unchanged)
    assert_eq!(correspondence.atoms().right_of(AtomId(0)), Some(AtomId(2)));
    assert_eq!(correspondence.atoms().right_of(AtomId(1)), Some(AtomId(3)));
    assert_eq!(correspondence.bonds().right_of(BondId(0)), Some(BondId(1)));
}

#[rstest]
fn test_molecule_ast_combine_from() {
    let mut left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ..Default::default()
    });
    let correspondence = left.combine_from(&right);

    assert_eq!(left.atoms().count(), 3);
    assert_eq!(left.bond(BondId(0)).atom_ids(), [AtomId(1), AtomId(2)]);
    assert_eq!(correspondence.atoms().right_of(AtomId(0)), Some(AtomId(1)));
    assert_eq!(correspondence.atoms().right_of(AtomId(1)), Some(AtomId(2)));
}

#[rstest]
fn test_molecule_ast_combine_from_storage() {
    let mut left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ..Default::default()
    });
    Arc::get_mut(&mut left.atoms).unwrap().reserve(2);
    Arc::get_mut(&mut left.bonds).unwrap().reserve(1);
    let atom_storage = left.atoms.as_ptr();
    let bond_storage = left.bonds.as_ptr();
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::F),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
        ..Default::default()
    });

    left.combine_from(&right);

    assert_eq!(left.atoms.as_ptr(), atom_storage);
    assert_eq!(left.bonds.as_ptr(), bond_storage);
}

#[rstest]
fn test_molecule_ast_combine_overlay() {
    let left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1)],
            AromaticSystemAst::from_electrons(vec![1, 1]),
        )],
        ..Default::default()
    });
    let (union, correspondence) = left.combine(&right);

    assert_eq!(union.aromatic_systems().count(), 1);
    // right's overlay over its atoms [0, 1] shifts by left's one atom
    assert_eq!(
        union
            .aromatic_system(AromaticSystemId(0))
            .atom_ids()
            .collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2)]
    );
    assert_eq!(
        correspondence
            .aromatic_systems()
            .right_of(AromaticSystemId(0)),
        Some(AromaticSystemId(0))
    );
}

#[rstest]
fn test_molecule_ast_combine_stereo() {
    let left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            ],
            StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let (union, _) = left.combine(&right);

    assert_eq!(union.stereo_atoms().count(), 1);
    let stereo = union.stereo_atoms().iter().next().unwrap();
    // right's site (atom 0) and ligands (atoms 1, 2) shift by left's one atom
    assert_eq!(stereo.site_id(), AtomId(1));
    assert_eq!(
        stereo.ligands().map(|l| l.atom_id()).collect::<Vec<_>>(),
        vec![AtomId(2), AtomId(3)]
    );
}

#[rstest]
fn test_molecule_ast_combine_constraint() {
    let left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
                sum: ValueAst::Lit(0),
            },
        )),
        ..Default::default()
    });
    let (union, _) = left.combine(&right);

    // right's constraint over atoms [0, 1] is remapped to [1, 2] in the union
    let expected = Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(1), AtomId(2)]),
        sum: ValueAst::Lit(0),
    });
    assert_eq!(
        union.constraints.iter().collect::<Vec<_>>(),
        vec![&expected]
    );
}

#[rstest]
fn test_molecule_ast_split() {
    // two disconnected bonds → two components
    let mol = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(2)),
        ],
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 2);
    let (first, first_corr) = &components[0];
    assert_eq!(first.atoms().count(), 2);
    assert_eq!(first.bond(BondId(0)).ast, &BondAst::from_order(1));
    assert_eq!(first_corr.atoms().right_of(AtomId(0)), Some(AtomId(0)));
    assert_eq!(first_corr.atoms().right_of(AtomId(1)), Some(AtomId(1)));
    let (second, second_corr) = &components[1];
    assert_eq!(second.bond(BondId(0)).ast, &BondAst::from_order(2));
    assert_eq!(second_corr.atoms().right_of(AtomId(0)), Some(AtomId(2)));
    assert_eq!(second_corr.atoms().right_of(AtomId(1)), Some(AtomId(3)));
}

#[rstest]
fn test_molecule_ast_split_overlay_binds() {
    // two disconnected bonds, but an aromatic system over {1, 2} keeps all four atoms in one component
    let mol = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(2)],
            AromaticSystemAst::from_electrons(vec![1, 1]),
        )],
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 1);
    assert_eq!(components[0].0.atoms().count(), 4);
}

#[rstest]
fn test_molecule_ast_combine_split_roundtrip() {
    let left = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ..Default::default()
    });
    let right = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
        ..Default::default()
    });
    let (union, _) = left.combine(&right);
    let components = union.split();

    assert_eq!(components.len(), 2);
    assert_eq!(components[0].0, left);
    assert_eq!(components[1].0, right);
}

#[rstest]
fn test_molecule_ast_split_stereo() {
    // a stereo atom binds its site + ligands into one component, separate from a lone bond
    let mol = MoleculeAst::from_entries(MoleculeEntries {
        atoms: (0..7).map(|_| AtomAst::from_element(Element::C)).collect(),
        bonds: vec![(AtomId(5), AtomId(6), BondAst::from_order(1))],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        )],
        constraints: Constraints::new(),
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 2);
    let (bound, _) = &components[0];
    assert_eq!(bound.atoms().count(), 5);
    assert_eq!(bound.stereo_atoms().count(), 1);
    let stereo = bound.stereo_atoms().iter().next().unwrap();
    assert_eq!(stereo.site_id(), AtomId(0));
    assert_eq!(
        stereo.ligands().map(|l| l.atom_id()).collect::<Vec<_>>(),
        vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)]
    );
    let (lone, _) = &components[1];
    assert_eq!(lone.atoms().count(), 2);
    assert_eq!(lone.stereo_atoms().count(), 0);
    assert_eq!(lone.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
}

#[rstest]
fn test_molecule_ast_split_constraint_binds() {
    // two disconnected bonds, but a ChargeSum over {1, 2} binds all four atoms into one component
    let mol = MoleculeAst::from_entries(MoleculeEntries {
        atoms: (0..4).map(|_| AtomAst::from_element(Element::C)).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
                sum: ValueAst::Lit(0),
            },
        )),
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 1);
    assert_eq!(
        components[0].0.constraints.iter().collect::<Vec<_>>(),
        vec![&Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(1), AtomId(2)]),
            sum: ValueAst::Lit(0),
        })]
    );
}

#[rstest]
fn test_molecule_ast_split_constraint_routed() {
    // a constraint over the second component's atoms routes there, remapped to compact ids
    let mol = MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(2), AtomId(3), BondAst::from_order(2)),
        ],
        constraints: constraints_with_molecule(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(2), AtomId(3)]),
                sum: ValueAst::Lit(0),
            },
        )),
        ..Default::default()
    });
    let components = mol.split();

    assert_eq!(components.len(), 2);
    assert!(components[0].0.constraints.is_empty());
    assert_eq!(
        components[1].0.constraints.iter().collect::<Vec<_>>(),
        vec![&Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: ValueAst::Lit(0),
        })]
    );
}
