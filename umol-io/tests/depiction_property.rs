//! Semantic properties of tetrahedral depiction.
//!
//! These properties exercise two independent frames: the stored graph-IR ligand order and the
//! geometric display frame. Reframing a literal tetrahedral form into any permutation of the same
//! distinct ligands must not change its depiction. Once a wedge is emitted, applying the inverse
//! TableIR winding convention to its kind and the displayed coordinates must recover the requested
//! tetrahedral coset. The geometry domain consists of a nondegenerate four-ligand star under finite
//! translations, positive integer scale, quarter-turn rotations, and reflection.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_geometric_core::{signed_volume, Point2D, Point3D};
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, Entity, FrameTransport, Molecule, MoleculeEntries, StereoAtomForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
};
use umol_io::depiction::molecule::depict;
use umol_io::depiction::{Depiction, DepictionItem, DepictionReference, WedgeItem, WedgeKind};
use umol_io::layout::MoleculeLayout;
use umol_perm::Permutation;

fn display_ligands() -> Vec<StereoLigand> {
    (1..=4)
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect()
}

fn molecule(stored_ligands: Vec<StereoLigand>, attributes: StereoAtomForm) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: (1..=4)
            .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(AtomId(0), stored_ligands, attributes)],
        ..Default::default()
    })
}

fn reframed_molecule(coset: u32, permutation_rank: usize) -> Molecule {
    let display_ligands = display_ligands();
    let stored_ligands = Permutation::unrank(4, permutation_rank).act(&display_ligands);
    let action = Permutation::between(&display_ligands, &stored_ligands)
        .expect("generated ligands are a permutation of the display frame");
    let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset)
        .reframe_by(&action)
        .expect("every degree-four permutation acts on a tetrahedral form");
    molecule(stored_ligands, attributes)
}

fn base_positions() -> Vec<Point2D> {
    vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(0.0, 1.0),
        Point2D::new(-1.0, 0.0),
        Point2D::new(0.0, -1.0),
    ]
}

fn transformed_layout(
    quarter_turns: u8,
    reflected: bool,
    scale: u8,
    translation: [i8; 2],
) -> MoleculeLayout {
    let positions = base_positions()
        .into_iter()
        .map(|point| {
            let mut x = point.x;
            let mut y = point.y;
            if reflected {
                x = -x;
            }
            for _ in 0..quarter_turns {
                (x, y) = (-y, x);
            }
            Point2D::new(
                x * f64::from(scale) + f64::from(translation[0]),
                y * f64::from(scale) + f64::from(translation[1]),
            )
        })
        .collect();
    MoleculeLayout::try_new(positions).expect("generated coordinates are finite")
}

fn wedge(depiction: &Depiction) -> &WedgeItem {
    let mut wedges = depiction.items().iter().filter_map(|item| match item {
        DepictionItem::Wedge(wedge) => Some(wedge),
        _ => None,
    });
    let wedge = wedges
        .next()
        .expect("literal tetrahedral fixture has a wedge");
    assert!(wedges.next().is_none(), "fixture has one stereo atom");
    wedge
}

fn recovered_coset(molecule: &Molecule, depiction: &Depiction, layout: &MoleculeLayout) -> u32 {
    let wedge = wedge(depiction);
    let wedged_atom = wedge
        .references
        .iter()
        .find_map(|reference| match reference {
            DepictionReference::Molecule(Entity::Bond(bond)) => molecule
                .bond(*bond)
                .atom_ids()
                .into_iter()
                .find(|&atom| atom != AtomId(0)),
            _ => None,
        })
        .expect("a tetrahedral wedge carries its display-bond reference");
    let points = (1..=4)
        .map(|atom| {
            let atom = AtomId(atom);
            let point = *layout
                .position(atom)
                .expect("ligand is in the layout frame");
            let z = if atom == wedged_atom {
                match wedge.kind {
                    WedgeKind::Solid => 1.0,
                    WedgeKind::Hashed => -1.0,
                }
            } else {
                0.0
            };
            Point3D::new(point.x, point.y, z)
        })
        .collect::<Vec<_>>();
    let [first, second, third, fourth] = points.as_slice() else {
        unreachable!("tetrahedral fixture has four ligands")
    };
    let volume = signed_volume(*first, *second, *third, *fourth);
    assert!(volume.is_finite() && volume != 0.0);
    u32::from(volume >= 0.0)
}

proptest! {
    #![proptest_config(Config {
        cases: 256,
        failure_persistence: Some(Box::new(FileFailurePersistence::WithSource(
            "proptest-regressions",
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_tetrahedral_depiction_is_invariant_under_stored_ligand_reframing(
        coset in 0u32..2,
        permutation_rank in 0usize..24,
    ) {
        let layout = MoleculeLayout::try_new(base_positions()).expect("base geometry is finite");
        let baseline = molecule(
            display_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, coset),
        );
        let reframed = reframed_molecule(coset, permutation_rank);

        prop_assert_eq!(
            depict(&reframed, &layout),
            depict(&baseline, &layout),
        );
        prop_assert_eq!(
            reframed
                .stereo_atoms()
                .iter()
                .next()
                .expect("fixture has one stereo atom")
                .coset_for(display_ligands()),
            Some(StereoCoset::Lit(coset)),
        );
    }

    #[test]
    fn test_tetrahedral_wedge_recovers_display_coset(
        coset in 0u32..2,
        permutation_rank in 0usize..24,
        quarter_turns in 0u8..4,
        reflected in any::<bool>(),
        scale in 1u8..5,
        translation in prop::array::uniform2(-8i8..=8),
    ) {
        let molecule = reframed_molecule(coset, permutation_rank);
        let layout = transformed_layout(quarter_turns, reflected, scale, translation);
        let depiction = depict(&molecule, &layout).expect("generated geometry is nondegenerate");

        prop_assert_eq!(recovered_coset(&molecule, &depiction, &layout), coset);
    }
}
