//! Molecule DSL reference-resolution properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_ast::dsl::{
    AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef, NoncovalentBondRef,
    StereoAtomRef, StereoBondRef,
};

use crate::strategies::*;

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
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
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
}
