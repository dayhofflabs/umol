use proptest::prelude::*;
use rstest::rstest;
// The DSL ref family (`Structural` + `resolve`), distinct from the `ast::edit` handles.
use umol_ast::ast::MoleculeCorrespondence;
use umol_ast::dsl::{
    AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef, NoncovalentBondRef,
    StereoAtomRef, StereoBondRef,
};
use umol_graph_core::{Correspondence, EdgeId, GraphCorrespondence, NodeId};

use crate::strategies::*;

fn identity_correspondence(ast: &MoleculeAst) -> MoleculeCorrespondence {
    fn identity<Id>(count: usize) -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        let images: Vec<Id> = (0..count).map(Id::from).collect();
        Correspondence::from_images(&images, count)
    }

    MoleculeCorrespondence::new(
        identity::<NodeId>(ast.atoms().count()),
        identity::<BondId>(ast.bonds().count()),
        identity::<DativeBondId>(ast.dative_bonds().count()),
        identity::<AromaticSystemId>(ast.aromatic_systems().count()),
        identity::<MulticenterBondId>(ast.multicenter_bonds().count()),
        identity::<NoncovalentBondId>(ast.noncovalent_bonds().count()),
        identity::<StereoAtomId>(ast.stereo_atoms().count()),
        identity::<StereoBondId>(ast.stereo_bonds().count()),
    )
}

proptest! {
    #[test]
    fn test_molecule_ast_equiv_reflexive(ast in molecule_ast_with_constraints_strategy()) {
        prop_assert!(ast.equiv(&ast));
    }

    #[test]
    fn test_molecule_ast_equiv_symmetric(
        left in molecule_ast_with_constraints_strategy(),
        right in molecule_ast_with_constraints_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), right.equiv(&left));
    }

    #[test]
    fn test_molecule_ast_equiv_agrees_with_equality_for_canonical_asts(
        left in molecule_ast_strategy(),
        right in molecule_ast_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), left == right);
    }

    #[test]
    fn test_molecule_ast_equiv_under_identity_reduces_to_equiv(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let correspondence = identity_correspondence(&ast);
        let mut other = ast.clone();
        if other.atoms().count() > 0 {
            other.atom_mut(AtomId(0)).ast.charge = ValueAst::Lit(99);
        }

        prop_assert_eq!(
            ast.equiv_under(&other, &correspondence),
            ast.equiv(&other),
        );
    }

    #[test]
    fn test_molecule_ast_equiv_under_symmetric_under_reverse(
        atoms in prop::collection::vec(atom_ast_strategy(), 0..=5),
        change_mapped_atom in any::<bool>(),
    ) {
        let count = atoms.len();
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            ..Default::default()
        });
        let mut right = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.into_iter().rev().collect(),
            ..Default::default()
        });
        if change_mapped_atom && count > 0 {
            right.atom_mut(AtomId((count - 1) as u32)).ast.charge = ValueAst::Lit(99);
        }
        let images: Vec<NodeId> = (0..count).rev().map(NodeId::from).collect();
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&images, count),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
        );

        let forward = left.equiv_under(&right, &correspondence);
        let reverse = right.equiv_under(&left, &correspondence.reverse());
        prop_assert_eq!(forward, reverse);
        prop_assert_eq!(forward, !change_mapped_atom || count == 0);
    }

    #[test]
    fn test_molecule_dsl_to_edn_from_edn_tree_roundtrip(dsl in molecule_dsl_strategy()) {
        let edn = dsl.to_edn();
        let parsed = MoleculeDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_to_edn_from_edn_str_roundtrip(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_molecule_dsl_streaming_matches_tree(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let via_stream = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming: {e}")))?;
        let tree = read_string(&rendered)
            .map_err(|e| TestCaseError::fail(format!("edn parse: {e}")))?;
        let via_tree = MoleculeDsl::from_edn(&tree)
            .map_err(|e| TestCaseError::fail(format!("tree: {e}")))?;
        prop_assert_eq!(via_stream, via_tree);
    }

    /// Direct `MoleculeAst::ToEdn` / `FromEdn` round-trips are the identity.
    /// Refs render as positional integers (no id keywords); the AST carries
    /// no metadata, so canonical EDN parses back to an equal AST.
    #[test]
    fn test_molecule_ast_to_edn_from_edn_tree_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let edn = ast.to_edn();
        let parsed = MoleculeAst::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("tree parse failed: {e}")))?;
        prop_assert_eq!(ast, parsed);
    }

    #[test]
    fn test_molecule_ast_to_edn_from_edn_str_roundtrip(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let rendered = ast.to_edn().to_string();
        let parsed = MoleculeAst::from_edn_str(&rendered)
            .map_err(|e| TestCaseError::fail(format!("streaming parse failed: {e}\nrendered: {rendered}")))?;
        prop_assert_eq!(ast, parsed);
    }

    #[test]
    fn test_molecule_dsl_double_render_is_stable(dsl in molecule_dsl_strategy()) {
        let s1 = dsl.to_edn().to_string();
        let d1 = MoleculeDsl::from_edn_str(&s1)
            .map_err(|e| TestCaseError::fail(format!("first parse: {e}")))?;
        let s2 = d1.to_edn().to_string();
        prop_assert_eq!(s1, s2);
    }
}

/// When `MoleculeMetadata` records an id for an entity, refs in molecule constraints
/// render as the keyword `:id` rather than the positional integer. Rendered
/// EDN must carry the keyword form, never the integer index, and must
/// roundtrip back through both the tree and streaming parsers.
#[rstest]
fn test_constraint_ref_uses_keyword_when_metadata_id_present() {
    let atoms = vec![AtomAst::default(), AtomAst::default()];
    let mut cs = Constraints::new();
    cs.push(Constraint::Atom(
        AtomId(0),
        AtomConstraintAst::Valence(ValueAst::Lit(4)),
    ));
    let ast = MoleculeAst::from_parts(MoleculeParts {
        atoms,
        constraints: cs,
        ..Default::default()
    });

    let mut metadata = MoleculeMetadata::new();
    metadata.set_atom_keyword(AtomId(0), "carbon".to_string());

    let dsl = MoleculeDsl::from_parts(ast, metadata);
    let rendered = dsl.to_edn().to_string();

    assert!(
        rendered.contains(":carbon"),
        "expected :carbon in rendered output: {rendered}",
    );
    assert!(
        !rendered.contains("[0 {:valence"),
        "rendered output must not use positional ref when id is present: {rendered}",
    );

    let via_tree = MoleculeDsl::from_edn(&dsl.to_edn()).expect("tree parse");
    assert_eq!(dsl, via_tree, "tree roundtrip");
    let via_stream = MoleculeDsl::from_edn_str(&rendered).expect("streaming parse");
    assert_eq!(dsl, via_stream, "streaming roundtrip");
}

/// Whether `keys[i]` occurs exactly once — the entity's participant set names it unambiguously.
fn is_unique<T: PartialEq>(keys: &[T], i: usize) -> bool {
    keys.iter().filter(|k| **k == keys[i]).count() == 1
}

/// Atom indices of a two-atom entity, sorted (the participant key is unordered).
fn sorted_pair(a: usize, b: usize) -> [usize; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

/// The atom-index set of a multi-atom entity (sorted, deduplicated).
fn atom_set(atoms: impl Iterator<Item = AtomId>) -> Vec<usize> {
    let mut set: Vec<usize> = atoms.map(|a| a.index()).collect();
    set.sort_unstable();
    set.dedup();
    set
}

proptest! {
    /// A structural ref (an entity named by its constituent atom/bond refs) resolves to the same id
    /// as the positional ref, for every non-atom entity whose participant set is unique — the seven
    /// `resolve_*_structural` paths cross-checked against positional resolution.
    #[test]
    fn test_structural_ref_resolves_like_positional(ast in molecule_ast_strategy()) {
        let ns = MoleculeNamespace::from_ast(&ast);

        let bond_keys: Vec<[usize; 2]> = ast
            .bonds()
            .iter()
            .map(|v| {
                let [a, b] = v.atom_ids();
                sorted_pair(a.index(), b.index())
            })
            .collect();
        for (i, view) in ast.bonds().iter().enumerate() {
            if !is_unique(&bond_keys, i) {
                continue;
            }
            let [a, b] = view.atom_ids();
            let structural =
                BondRef::Structural([AtomRef::Index(a.index()), AtomRef::Index(b.index())]);
            prop_assert_eq!(structural.resolve(&ns), Ok(BondId(i as u32)));
            prop_assert_eq!(BondRef::Index(i).resolve(&ns), Ok(BondId(i as u32)));
        }

        let dative_keys: Vec<(Vec<usize>, usize)> = ast
            .dative_bonds()
            .iter()
            .map(|v| (atom_set(v.donor_ids()), v.acceptor_id().index()))
            .collect();
        for (i, view) in ast.dative_bonds().iter().enumerate() {
            if !is_unique(&dative_keys, i) {
                continue;
            }
            let donors = view.donor_ids().map(|a| AtomRef::Index(a.index())).collect();
            let acceptor = AtomRef::Index(view.acceptor_id().index());
            let structural = DativeBondRef::Structural(DativeBondParticipants { donors, acceptor });
            prop_assert_eq!(structural.resolve(&ns), Ok(DativeBondId(i as u32)));
        }

        let aromatic_keys: Vec<Vec<usize>> = ast
            .aromatic_systems()
            .iter()
            .map(|v| atom_set(v.atom_ids()))
            .collect();
        for (i, view) in ast.aromatic_systems().iter().enumerate() {
            if !is_unique(&aromatic_keys, i) {
                continue;
            }
            let atoms = view.atom_ids().map(|a| AtomRef::Index(a.index())).collect();
            let structural = AromaticSystemRef::Structural(atoms);
            prop_assert_eq!(structural.resolve(&ns), Ok(AromaticSystemId(i as u32)));
        }

        let multicenter_keys: Vec<Vec<usize>> = ast
            .multicenter_bonds()
            .iter()
            .map(|v| atom_set(v.atom_ids()))
            .collect();
        for (i, view) in ast.multicenter_bonds().iter().enumerate() {
            if !is_unique(&multicenter_keys, i) {
                continue;
            }
            let atoms = view.atom_ids().map(|a| AtomRef::Index(a.index())).collect();
            let structural = MulticenterBondRef::Structural(atoms);
            prop_assert_eq!(structural.resolve(&ns), Ok(MulticenterBondId(i as u32)));
        }

        let noncovalent_keys: Vec<[usize; 2]> = ast
            .noncovalent_bonds()
            .iter()
            .map(|v| {
                let [a, b] = v.atom_ids();
                sorted_pair(a.index(), b.index())
            })
            .collect();
        for (i, view) in ast.noncovalent_bonds().iter().enumerate() {
            if !is_unique(&noncovalent_keys, i) {
                continue;
            }
            let [a, b] = view.atom_ids();
            let structural =
                NoncovalentBondRef::Structural([AtomRef::Index(a.index()), AtomRef::Index(b.index())]);
            prop_assert_eq!(structural.resolve(&ns), Ok(NoncovalentBondId(i as u32)));
        }

        // Stereo sites are distinct (one element per site), so keys never collide.
        for (i, view) in ast.stereo_atoms().iter().enumerate() {
            let site = AtomRef::Index(view.site_id().index());
            let ligands = view
                .ligands()
                .map(|l| StereoLigandRef {
                    kind: l.kind(),
                    atom: AtomRef::Index(l.atom_id().index()),
                })
                .collect();
            let structural = StereoAtomRef::Structural(StereoAtomParticipants { site, ligands });
            prop_assert_eq!(structural.resolve(&ns), Ok(StereoAtomId(i as u32)));
        }

        for (i, view) in ast.stereo_bonds().iter().enumerate() {
            let site = BondRef::Index(view.site_id().index());
            let ligands = view
                .ligands()
                .map(|l| StereoLigandRef {
                    kind: l.kind(),
                    atom: AtomRef::Index(l.atom_id().index()),
                })
                .collect();
            let structural = StereoBondRef::Structural(StereoBondParticipants { site, ligands });
            prop_assert_eq!(structural.resolve(&ns), Ok(StereoBondId(i as u32)));
        }
    }

    /// A structural ref over a set that is not any entity's participant set fails to resolve
    /// (`InvalidRef`), never silently hitting a wrong id — one guaranteed-miss perturbation per kind.
    #[test]
    fn test_structural_ref_wrong_participants_error(ast in molecule_ast_strategy()) {
        let ns = MoleculeNamespace::from_ast(&ast);

        // A two-atom entity: a self-pair is never a bond (endpoints are distinct).
        for view in ast.bonds().iter() {
            let [a, _] = view.atom_ids();
            let wrong = BondRef::Structural([AtomRef::Index(a.index()), AtomRef::Index(a.index())]);
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "bond", .. })
            ), "structural ref over wrong participants must not resolve");
        }
        for view in ast.noncovalent_bonds().iter() {
            let [a, _] = view.atom_ids();
            let wrong =
                NoncovalentBondRef::Structural([AtomRef::Index(a.index()), AtomRef::Index(a.index())]);
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "noncovalent-bond", .. })
            ), "structural ref over wrong participants must not resolve");
        }

        // A single-atom set is never a multi-atom entity (aromatic/multicenter need ≥ 3 atoms).
        for view in ast.aromatic_systems().iter() {
            let first = view.atom_ids().next().expect("aromatic systems are non-empty");
            let wrong = AromaticSystemRef::Structural(vec![AtomRef::Index(first.index())]);
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "aromatic-system", .. })
            ), "structural ref over wrong participants must not resolve");
        }
        for view in ast.multicenter_bonds().iter() {
            let first = view.atom_ids().next().expect("multicenter bonds are non-empty");
            let wrong = MulticenterBondRef::Structural(vec![AtomRef::Index(first.index())]);
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "multicenter-bond", .. })
            ), "structural ref over wrong participants must not resolve");
        }

        // The acceptor coinciding with a donor is never a real dative bond.
        for view in ast.dative_bonds().iter() {
            let donors = view.donor_ids().map(|a| AtomRef::Index(a.index())).collect();
            let acceptor = AtomRef::Index(view.donor_ids().next().expect("≥ 1 donor").index());
            let wrong = DativeBondRef::Structural(DativeBondParticipants { donors, acceptor });
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "dative-bond", .. })
            ), "structural ref over wrong participants must not resolve");
        }

        // Repeating a ligand changes the multiset, which no real stereo element carries.
        for view in ast.stereo_atoms().iter() {
            let site = AtomRef::Index(view.site_id().index());
            let mut ligands: Vec<StereoLigandRef> = view
                .ligands()
                .map(|l| StereoLigandRef {
                    kind: l.kind(),
                    atom: AtomRef::Index(l.atom_id().index()),
                })
                .collect();
            let Some(first) = ligands.first().cloned() else {
                continue;
            };
            ligands.push(first);
            let wrong = StereoAtomRef::Structural(StereoAtomParticipants { site, ligands });
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "stereo-atom", .. })
            ), "structural ref over wrong participants must not resolve");
        }
        for view in ast.stereo_bonds().iter() {
            let site = BondRef::Index(view.site_id().index());
            let mut ligands: Vec<StereoLigandRef> = view
                .ligands()
                .map(|l| StereoLigandRef {
                    kind: l.kind(),
                    atom: AtomRef::Index(l.atom_id().index()),
                })
                .collect();
            let Some(first) = ligands.first().cloned() else {
                continue;
            };
            ligands.push(first);
            let wrong = StereoBondRef::Structural(StereoBondParticipants { site, ligands });
            prop_assert!(matches!(
                wrong.resolve(&ns),
                Err(ParseError::InvalidRef { kind: "stereo-bond", .. })
            ), "structural ref over wrong participants must not resolve");
        }
    }

    #[test]
    fn test_molecule_ast_meet_pushout_reframes_stereo_atom(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                permutation.act(&left_frame),
                left_ast.apply(permutation),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..4u32).map(EdgeId).collect::<Vec<_>>(), 4),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap).map(|pushout| pushout.object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_ast_meet_pushout_rejects_changed_stereo_atom_ligand(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(5), StereoLigandKind::Atom);
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_atoms: vec![(AtomId(0), right_frame, left_ast.apply(permutation))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..4u32).map(EdgeId).collect::<Vec<_>>(), 4),
        );

        prop_assert!(left.meet_pushout(&right, &overlap).is_none());
    }

    #[test]
    fn test_molecule_ast_meet_pushout_reframes_stereo_bond(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                permutation.act(&left_frame),
                left_ast.apply(permutation),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap).map(|pushout| pushout.object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_ast_meet_pushout_rejects_changed_stereo_bond_ligand(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(6), StereoLigandKind::Atom);
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(BondId(0), right_frame, left_ast.apply(permutation))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..7u32).map(NodeId).collect::<Vec<_>>(), 7),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert!(left.meet_pushout(&right, &overlap).is_none());
    }
}
