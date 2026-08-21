//! Graph-IR projection for the CoordGen molecule-layout backend.

use umol_coordgen_sys::{generate_coordinates, Bond, CoordgenError};
use umol_graph_ir::ir::{AsLit, Molecule};

use super::{MoleculeLayout, Point2D};

// CoordGen's fixed working scale is defined by BONDLENGTH in sketcherMinimizerMaths.h.
const COORDGEN_BOND_LENGTH: f64 = 50.0;

pub(crate) fn layout(molecule: &Molecule) -> Result<MoleculeLayout, CoordgenError> {
    let atomic_numbers = molecule
        .atoms()
        .iter()
        .map(|atom| {
            atom.element()
                .as_lit()
                .map(|element| u16::from(element.atomic_number()))
                .unwrap_or(0)
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let [atom_0, atom_1] = bond.atom_ids();
            Bond {
                atom_0: atom_0.index(),
                atom_1: atom_1.index(),
                order: match bond.order().as_lit() {
                    Some(order @ 1..=3) => order as u8,
                    _ => 1,
                },
            }
        })
        .collect::<Vec<_>>();

    let positions = generate_coordinates(&atomic_numbers, &bonds)?
        .into_iter()
        .map(|point| {
            Point2D::new(
                point.x / COORDGEN_BOND_LENGTH,
                point.y / COORDGEN_BOND_LENGTH,
            )
        })
        .collect();

    Ok(MoleculeLayout::try_new(positions)
        .expect("scaling finite CoordGen points by its finite bond length preserves finiteness"))
}
