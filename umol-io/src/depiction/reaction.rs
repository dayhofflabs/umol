//! Indexed-side depiction construction from two molecules, layouts, and an atom correspondence.

use std::any::Any;

use thiserror::Error;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomId, Entity, Molecule};
use umol_utils::error::UmolError;

use super::{molecule, ArrowItem, Depiction, DepictionItem, DepictionReference, TextItem};
use crate::layout::{MoleculeLayout, MoleculeLayoutError, Point2D};

const ARROW_HALF_LENGTH: f64 = 1.0;
const SIDE_ARROW_GAP: f64 = 1.0;
const MAP_INDEX_OFFSET: Point2D = Point2D::new(0.35, 0.35);

/// Constructs an indexed-side reaction depiction from independently laid-out molecular sides.
///
/// Each side is translated without rotation or regeneration so that it is vertically centered and
/// separated from a fixed horizontal reaction arrow. Matching atom pairs receive equal, zero-based
/// display indices in correspondence left-id order. Unmatched atoms receive no index. The indices
/// are depiction-local labels, not atom ids.
///
/// This operation does not validate chemistry, induce a full molecule correspondence, align
/// conserved substructures, or canonicalize either side.
///
/// # Establishes
///
/// The result uses lhs/rhs reaction references for every lowered molecular item, places the arrow
/// between the translated side layouts, and represents each matched atom on both sides with the
/// same [`DepictionReference::CorrespondencePair`] index.
///
/// # Errors
///
/// Returns [`ReactionDepictionError`] if either layout or either atom-correspondence frame size
/// disagrees with its supplied molecule.
///
/// # Semantic properties
///
/// Side-item order is preserved. Correspondence-pair indices depend only on the correspondence's
/// left-id order, so identical inputs produce structurally equal depictions. Each supplied layout
/// is changed only by one translation.
pub fn depict(
    lhs: &Molecule,
    lhs_layout: &MoleculeLayout,
    rhs: &Molecule,
    rhs_layout: &MoleculeLayout,
    atom_correspondence: &Correspondence<AtomId>,
) -> Result<Depiction, ReactionDepictionError> {
    let lhs_depiction = molecule::depict(lhs, lhs_layout).map_err(lhs_layout_error)?;
    let rhs_depiction = molecule::depict(rhs, rhs_layout).map_err(rhs_layout_error)?;

    let lhs_atom_count = lhs.atoms().count();
    if atom_correspondence.left_count() != lhs_atom_count {
        return Err(ReactionDepictionError::LhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: lhs_atom_count,
            correspondence_atom_count: atom_correspondence.left_count(),
        });
    }
    let rhs_atom_count = rhs.atoms().count();
    if atom_correspondence.right_count() != rhs_atom_count {
        return Err(ReactionDepictionError::RhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: rhs_atom_count,
            correspondence_atom_count: atom_correspondence.right_count(),
        });
    }

    let lhs_offset = side_offset(lhs_layout, ReactionSide::Lhs);
    let rhs_offset = side_offset(rhs_layout, ReactionSide::Rhs);
    let mut items = lhs_depiction
        .items()
        .iter()
        .cloned()
        .map(|item| translate_item(item, lhs_offset, ReactionSide::Lhs))
        .collect::<Vec<_>>();

    items.extend(index_items(
        lhs_layout,
        lhs_offset,
        atom_correspondence
            .matched_pairs()
            .iter()
            .enumerate()
            .map(|(index, &(left, _))| (index, left)),
        ReactionSide::Lhs,
    ));
    items.push(DepictionItem::Arrow(ArrowItem {
        start: Point2D::new(-ARROW_HALF_LENGTH, 0.0),
        end: Point2D::new(ARROW_HALF_LENGTH, 0.0),
        references: Vec::new(),
    }));
    items.extend(
        rhs_depiction
            .items()
            .iter()
            .cloned()
            .map(|item| translate_item(item, rhs_offset, ReactionSide::Rhs)),
    );
    items.extend(index_items(
        rhs_layout,
        rhs_offset,
        atom_correspondence
            .matched_pairs()
            .iter()
            .enumerate()
            .map(|(index, &(_, right))| (index, right)),
        ReactionSide::Rhs,
    ));

    Ok(Depiction::from_items(items))
}

/// Failures caused by mismatched molecular, layout, or atom-correspondence frames.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ReactionDepictionError {
    #[error(
        "lhs molecule atom count {molecule_atom_count} does not match layout atom count {layout_atom_count}"
    )]
    LhsLayoutFrameSizeMismatch {
        molecule_atom_count: usize,
        layout_atom_count: usize,
    },
    #[error(
        "rhs molecule atom count {molecule_atom_count} does not match layout atom count {layout_atom_count}"
    )]
    RhsLayoutFrameSizeMismatch {
        molecule_atom_count: usize,
        layout_atom_count: usize,
    },
    #[error(
        "lhs molecule atom count {molecule_atom_count} does not match correspondence atom count {correspondence_atom_count}"
    )]
    LhsCorrespondenceFrameSizeMismatch {
        molecule_atom_count: usize,
        correspondence_atom_count: usize,
    },
    #[error(
        "rhs molecule atom count {molecule_atom_count} does not match correspondence atom count {correspondence_atom_count}"
    )]
    RhsCorrespondenceFrameSizeMismatch {
        molecule_atom_count: usize,
        correspondence_atom_count: usize,
    },
}

impl UmolError for ReactionDepictionError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Clone, Copy)]
enum ReactionSide {
    Lhs,
    Rhs,
}

fn lhs_layout_error(error: MoleculeLayoutError) -> ReactionDepictionError {
    match error {
        MoleculeLayoutError::FrameSizeMismatch {
            molecule_atom_count,
            layout_atom_count,
        } => ReactionDepictionError::LhsLayoutFrameSizeMismatch {
            molecule_atom_count,
            layout_atom_count,
        },
        MoleculeLayoutError::NonFinitePosition { .. }
        | MoleculeLayoutError::AtomOutOfFrame { .. } => {
            unreachable!("molecule depiction checks only layout frame agreement")
        }
    }
}

fn rhs_layout_error(error: MoleculeLayoutError) -> ReactionDepictionError {
    match error {
        MoleculeLayoutError::FrameSizeMismatch {
            molecule_atom_count,
            layout_atom_count,
        } => ReactionDepictionError::RhsLayoutFrameSizeMismatch {
            molecule_atom_count,
            layout_atom_count,
        },
        MoleculeLayoutError::NonFinitePosition { .. }
        | MoleculeLayoutError::AtomOutOfFrame { .. } => {
            unreachable!("molecule depiction checks only layout frame agreement")
        }
    }
}

fn side_offset(layout: &MoleculeLayout, side: ReactionSide) -> Point2D {
    let Some(first) = layout.positions().first() else {
        return Point2D::new(0.0, 0.0);
    };
    let mut min = *first;
    let mut max = *first;
    for &position in &layout.positions()[1..] {
        min.x = min.x.min(position.x);
        min.y = min.y.min(position.y);
        max.x = max.x.max(position.x);
        max.y = max.y.max(position.y);
    }
    let x = match side {
        ReactionSide::Lhs => -ARROW_HALF_LENGTH - SIDE_ARROW_GAP - max.x,
        ReactionSide::Rhs => ARROW_HALF_LENGTH + SIDE_ARROW_GAP - min.x,
    };
    Point2D::new(x, -(min.y + max.y) / 2.0)
}

fn translate_item(mut item: DepictionItem, offset: Point2D, side: ReactionSide) -> DepictionItem {
    let references = match &mut item {
        DepictionItem::Atom(item) => {
            item.position = translate(item.position, offset);
            &mut item.references
        }
        DepictionItem::Bond(item) => {
            item.start = translate(item.start, offset);
            item.end = translate(item.end, offset);
            &mut item.references
        }
        DepictionItem::Text(item) => {
            item.position = translate(item.position, offset);
            &mut item.references
        }
        DepictionItem::Marker(item) => {
            item.position = translate(item.position, offset);
            &mut item.references
        }
        DepictionItem::Arrow(item) => {
            item.start = translate(item.start, offset);
            item.end = translate(item.end, offset);
            &mut item.references
        }
    };
    for reference in references {
        if let DepictionReference::Molecule(entity) = *reference {
            *reference = reaction_reference(side, entity);
        }
    }
    item
}

fn index_items(
    layout: &MoleculeLayout,
    offset: Point2D,
    indexed_atoms: impl IntoIterator<Item = (usize, AtomId)>,
    side: ReactionSide,
) -> Vec<DepictionItem> {
    indexed_atoms
        .into_iter()
        .map(|(index, atom)| {
            let position = *layout
                .position(atom)
                .expect("correspondence frame agreement establishes every atom position");
            DepictionItem::Text(TextItem {
                position: translate(translate(position, offset), MAP_INDEX_OFFSET),
                text: index.to_string(),
                references: vec![
                    reaction_reference(side, Entity::Atom(atom)),
                    DepictionReference::CorrespondencePair(index as u32),
                ],
            })
        })
        .collect()
}

fn reaction_reference(side: ReactionSide, entity: Entity) -> DepictionReference {
    match side {
        ReactionSide::Lhs => DepictionReference::ReactionLhs(entity),
        ReactionSide::Rhs => DepictionReference::ReactionRhs(entity),
    }
}

fn translate(point: Point2D, offset: Point2D) -> Point2D {
    Point2D::new(point.x + offset.x, point.y + offset.y)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::Correspondence;
    use umol_graph_ir::ir::{AtomId, BondId, Entity};
    use umol_graph_ir::mol_dsl;

    use super::{depict, ReactionDepictionError};
    use crate::depiction::{
        ArrowItem, AtomItem, BondItem, Bounds, DepictionItem, DepictionReference, TextItem,
    };
    use crate::layout::{MoleculeLayout, Point2D};

    #[rstest]
    fn test_depict() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let rhs = mol_dsl!(r#"{:atoms ["O" "C"] :bonds [[0 1 "1"]]}"#);
        let lhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let rhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let correspondence =
            Correspondence::new(vec![(AtomId(1), AtomId(0)), (AtomId(0), AtomId(1))], 2, 2)
                .unwrap();

        let depiction = depict(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        assert_eq!(
            depiction.items(),
            [
                DepictionItem::Bond(BondItem {
                    start: Point2D::new(-3.0, 0.0),
                    end: Point2D::new(-2.0, 0.0),
                    line_count: 1,
                    references: vec![DepictionReference::ReactionLhs(Entity::Bond(BondId(0),))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(-3.0, 0.0),
                    label: "C".to_owned(),
                    references: vec![DepictionReference::ReactionLhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(-2.0, 0.0),
                    label: "O".to_owned(),
                    references: vec![DepictionReference::ReactionLhs(Entity::Atom(AtomId(1)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-2.65, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-1.65, 0.35),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                }),
                DepictionItem::Arrow(ArrowItem {
                    start: Point2D::new(-1.0, 0.0),
                    end: Point2D::new(1.0, 0.0),
                    references: Vec::new(),
                }),
                DepictionItem::Bond(BondItem {
                    start: Point2D::new(2.0, 0.0),
                    end: Point2D::new(3.0, 0.0),
                    line_count: 1,
                    references: vec![DepictionReference::ReactionRhs(Entity::Bond(BondId(0),))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(2.0, 0.0),
                    label: "O".to_owned(),
                    references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(3.0, 0.0),
                    label: "C".to_owned(),
                    references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(1)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(3.35, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(2.35, 0.35),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                }),
            ]
        );
        assert_eq!(
            depiction.bounds(),
            Some(&Bounds {
                min: Point2D::new(-3.0, 0.0),
                max: Point2D::new(3.35, 0.35),
            })
        );
        assert_eq!(
            depict(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence,),
            Ok(depiction)
        );
    }

    #[rstest]
    fn test_depict_partial() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "F" "O"] :bonds []}"#);
        let rhs = mol_dsl!(r#"{:atoms ["O" "C"] :bonds [[0 1 "1"]]}"#);
        let lhs_layout = MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
        ])
        .unwrap();
        let rhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let correspondence =
            Correspondence::new(vec![(AtomId(2), AtomId(0)), (AtomId(0), AtomId(1))], 3, 2)
                .unwrap();

        let depiction = depict(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();
        let index_items = depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Text(item) => Some(item.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(
            index_items,
            [
                TextItem {
                    position: Point2D::new(-3.65, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                },
                TextItem {
                    position: Point2D::new(-1.65, 0.35),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(2))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                },
                TextItem {
                    position: Point2D::new(3.35, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                },
                TextItem {
                    position: Point2D::new(2.35, 0.35),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                },
            ]
        );
    }

    #[rstest]
    fn test_depict_projection_omission() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "*"] :bonds []}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let lhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(10.0, 0.0)]).unwrap();
        let rhs_layout = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap();
        let correspondence = Correspondence::new(vec![(AtomId(1), AtomId(0))], 2, 1).unwrap();

        let depiction = depict(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        assert_eq!(
            depiction.items(),
            [
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(-12.0, 0.0),
                    label: "C".to_owned(),
                    references: vec![DepictionReference::ReactionLhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-1.65, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Arrow(ArrowItem {
                    start: Point2D::new(-1.0, 0.0),
                    end: Point2D::new(1.0, 0.0),
                    references: Vec::new(),
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(2.0, 0.0),
                    label: "C".to_owned(),
                    references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(2.35, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
            ]
        );
    }

    #[rstest]
    #[case::lhs_layout(
        MoleculeLayout::try_new(Vec::new()).unwrap(),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        Correspondence::new(vec![(AtomId(0), AtomId(0))], 1, 1).unwrap(),
        ReactionDepictionError::LhsLayoutFrameSizeMismatch {
            molecule_atom_count: 1,
            layout_atom_count: 0,
        }
    )]
    #[case::rhs_layout(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(Vec::new()).unwrap(),
        Correspondence::new(vec![(AtomId(0), AtomId(0))], 1, 1).unwrap(),
        ReactionDepictionError::RhsLayoutFrameSizeMismatch {
            molecule_atom_count: 1,
            layout_atom_count: 0,
        }
    )]
    #[case::lhs_correspondence(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        Correspondence::new(Vec::new(), 0, 1).unwrap(),
        ReactionDepictionError::LhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 1,
            correspondence_atom_count: 0,
        }
    )]
    #[case::rhs_correspondence(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        Correspondence::new(Vec::new(), 1, 0).unwrap(),
        ReactionDepictionError::RhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 1,
            correspondence_atom_count: 0,
        }
    )]
    fn test_depict_error(
        #[case] lhs_layout: MoleculeLayout,
        #[case] rhs_layout: MoleculeLayout,
        #[case] correspondence: Correspondence<AtomId>,
        #[case] expected: ReactionDepictionError,
    ) {
        let lhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);

        assert_eq!(
            depict(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence,),
            Err(expected)
        );
    }
}
