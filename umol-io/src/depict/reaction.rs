//! Indexed-side depiction from two molecules, layouts, and an atom correspondence.

use std::any::Any;
use std::f64::consts::{PI, SQRT_2, TAU};

use thiserror::Error;
use umol_geometric_core::Point2D;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomId, Entity, Molecule};
#[cfg(feature = "coordgen")]
use umol_graph_ir::ir::{Contradiction, Reaction};
use umol_utils::error::UmolError;

#[cfg(feature = "coordgen")]
use super::Depict;
use super::{molecule, ArrowItem, Depiction, DepictionItem, DepictionReference, TextItem};
use crate::depict::molecule::MoleculeDepictionError;
use crate::layout::MoleculeLayout;
#[cfg(feature = "coordgen")]
use crate::layout::{layout_molecule, MoleculeLayoutAlgorithm};

const ARROW_HALF_LENGTH: f64 = 0.75;
const SIDE_ARROW_GAP: f64 = 1.0;
const MAP_INDEX_COMPONENT_OFFSET: f64 = 0.35;
const MAP_INDEX_DISTANCE: f64 = MAP_INDEX_COMPONENT_OFFSET * SQRT_2;
const MAP_INDEX_FALLBACK_OFFSET: Point2D =
    Point2D::new(MAP_INDEX_COMPONENT_OFFSET, MAP_INDEX_COMPONENT_OFFSET);

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
/// Returns [`DepictFromSidesError`] if either molecular depiction fails or either
/// atom-correspondence frame size disagrees with its supplied molecule.
///
/// # Semantic properties
///
/// Side-item order is preserved. Correspondence-pair indices depend only on the correspondence's
/// left-id order, so identical inputs produce structurally equal depictions. Each supplied layout
/// is changed only by one translation. Index positions are selected atom-locally away from
/// incident bond directions and do not alter either molecular layout.
pub fn depict_from_sides(
    lhs: &Molecule,
    lhs_layout: &MoleculeLayout,
    rhs: &Molecule,
    rhs_layout: &MoleculeLayout,
    atom_correspondence: &Correspondence<AtomId>,
) -> Result<Depiction, DepictFromSidesError> {
    let lhs_depiction =
        molecule::depict(lhs, lhs_layout).map_err(DepictFromSidesError::LhsDepiction)?;
    let rhs_depiction =
        molecule::depict(rhs, rhs_layout).map_err(DepictFromSidesError::RhsDepiction)?;

    let lhs_atom_count = lhs.atoms().count();
    if atom_correspondence.left_count() != lhs_atom_count {
        return Err(DepictFromSidesError::LhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: lhs_atom_count,
            correspondence_atom_count: atom_correspondence.left_count(),
        });
    }
    let rhs_atom_count = rhs.atoms().count();
    if atom_correspondence.right_count() != rhs_atom_count {
        return Err(DepictFromSidesError::RhsCorrespondenceFrameSizeMismatch {
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
        lhs,
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
        rhs,
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

/// Generates both molecular layouts and composes an indexed-side reaction depiction.
///
/// The selected molecule-layout algorithm is applied independently to each side. The atom
/// correspondence is passed only to [`depict_from_sides`], where it determines the displayed map
/// indices; it does not influence either generated layout.
///
/// # Establishes
///
/// The result is structurally equal to generating each side with [`layout_molecule`] and passing
/// those layouts, the molecules, and the correspondence to [`depict_from_sides`].
///
/// # Errors
///
/// Returns [`DepictFromSidesError::LhsDepiction`] or
/// [`DepictFromSidesError::RhsDepiction`] if layout generation or molecule depiction fails for the
/// corresponding side. Returns a frame-size error if the atom correspondence does not agree with
/// the molecular sides.
///
/// # Semantic properties
///
/// Correspondence-pair indices retain the zero-based, left-id ordering defined by
/// [`depict_from_sides`].
#[cfg(feature = "coordgen")]
pub fn depict_from_sides_with(
    lhs: &Molecule,
    rhs: &Molecule,
    atom_correspondence: &Correspondence<AtomId>,
    layout_algorithm: MoleculeLayoutAlgorithm,
) -> Result<Depiction, DepictFromSidesError> {
    let lhs_layout = layout_molecule(lhs, layout_algorithm)
        .map_err(MoleculeDepictionError::Layout)
        .map_err(DepictFromSidesError::LhsDepiction)?;
    let rhs_layout = layout_molecule(rhs, layout_algorithm)
        .map_err(MoleculeDepictionError::Layout)
        .map_err(DepictFromSidesError::RhsDepiction)?;

    depict_from_sides(lhs, &lhs_layout, rhs, &rhs_layout, atom_correspondence)
}

#[cfg(feature = "coordgen")]
impl Depict for Reaction {
    type Error = ReactionDepictionError;

    fn depict_with(
        &self,
        layout_algorithm: MoleculeLayoutAlgorithm,
    ) -> Result<Depiction, Self::Error> {
        let span = self
            .to_reaction_span()
            .map_err(ReactionDepictionError::Materialization)?;
        let lhs = span.lhs();
        let rhs = span.rhs();
        let correspondence = span.correspondence();
        let lhs_layout = layout_molecule(&lhs, layout_algorithm)
            .map_err(MoleculeDepictionError::Layout)
            .map_err(ReactionDepictionError::LhsDepiction)?;
        let rhs_layout = layout_molecule(&rhs, layout_algorithm)
            .map_err(MoleculeDepictionError::Layout)
            .map_err(ReactionDepictionError::RhsDepiction)?;

        depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, correspondence.atoms()).map_err(
            |error| match error {
                DepictFromSidesError::LhsDepiction(error) => {
                    ReactionDepictionError::LhsDepiction(error)
                }
                DepictFromSidesError::RhsDepiction(error) => {
                    ReactionDepictionError::RhsDepiction(error)
                }
                DepictFromSidesError::LhsCorrespondenceFrameSizeMismatch { .. }
                | DepictFromSidesError::RhsCorrespondenceFrameSizeMismatch { .. } => {
                    unreachable!(
                        "reaction materialization establishes matching correspondence frames"
                    )
                }
            },
        )
    }
}

/// Failures while depicting independently supplied molecular sides.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum DepictFromSidesError {
    #[error("lhs depiction: {0}")]
    LhsDepiction(#[source] MoleculeDepictionError),
    #[error("rhs depiction: {0}")]
    RhsDepiction(#[source] MoleculeDepictionError),
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

impl UmolError for DepictFromSidesError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Failures while depicting a [`Reaction`].
#[cfg(feature = "coordgen")]
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ReactionDepictionError {
    #[error("reaction materialization: {0}")]
    Materialization(#[source] Contradiction),
    #[error("lhs depiction: {0}")]
    LhsDepiction(#[source] MoleculeDepictionError),
    #[error("rhs depiction: {0}")]
    RhsDepiction(#[source] MoleculeDepictionError),
}

#[cfg(feature = "coordgen")]
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
        DepictionItem::Wedge(item) => {
            item.tip = translate(item.tip, offset);
            item.base = translate(item.base, offset);
            &mut item.references
        }
        DepictionItem::DashedContour(item) => {
            for point in &mut item.points {
                *point = translate(*point, offset);
            }
            &mut item.references
        }
        DepictionItem::Text(item) => {
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
    molecule: &Molecule,
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
                position: translate(
                    translate(position, offset),
                    mapping_index_offset(molecule, layout, atom),
                ),
                text: index.to_string(),
                references: vec![
                    reaction_reference(side, Entity::Atom(atom)),
                    DepictionReference::CorrespondencePair(index as u32),
                ],
            })
        })
        .collect()
}

fn mapping_index_offset(molecule: &Molecule, layout: &MoleculeLayout, atom: AtomId) -> Point2D {
    let origin = *layout
        .position(atom)
        .expect("molecule/layout frame agreement establishes the atom position");
    let mut angles = molecule
        .neighbors(atom)
        .filter_map(|neighbor| {
            let position = *layout
                .position(neighbor.atom_id())
                .expect("molecule/layout frame agreement establishes every neighbor position");
            let dx = position.x - origin.x;
            let dy = position.y - origin.y;
            let length = dx.hypot(dy);
            (length > f64::EPSILON).then(|| {
                let angle = dy.atan2(dx);
                if angle < 0.0 {
                    angle + TAU
                } else {
                    angle
                }
            })
        })
        .collect::<Vec<_>>();

    let angle = match angles.len() {
        0 => return MAP_INDEX_FALLBACK_OFFSET,
        1 => angles[0] + PI,
        _ => {
            angles.sort_by(f64::total_cmp);
            let mut largest_start = angles[0];
            let mut largest_extent = 0.0;
            for index in 0..angles.len() {
                let start = angles[index];
                let end = if index + 1 == angles.len() {
                    angles[0] + TAU
                } else {
                    angles[index + 1]
                };
                let extent = end - start;
                if extent > largest_extent {
                    largest_start = start;
                    largest_extent = extent;
                }
            }
            largest_start + largest_extent / 2.0
        }
    };
    let x = MAP_INDEX_DISTANCE * angle.cos();
    let y = MAP_INDEX_DISTANCE * angle.sin();
    Point2D::new(
        if x.abs() < f64::EPSILON { 0.0 } else { x },
        if y.abs() < f64::EPSILON { 0.0 } else { y },
    )
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
    use float_cmp::approx_eq;
    use rstest::rstest;
    use umol_geometric_core::Point2D;
    use umol_graph_core::Correspondence;
    use umol_graph_ir::ir::{AtomId, BondId, Entity, Molecule, StereoAtomId};
    #[cfg(feature = "coordgen")]
    use umol_graph_ir::ir::{
        BondDelta, BondFieldChange, Contradiction, Delta, Deltas, NumForm, Reaction,
    };
    use umol_graph_ir::mol_dsl;

    #[cfg(feature = "coordgen")]
    use super::depict_from_sides_with;
    #[cfg(feature = "coordgen")]
    use super::ReactionDepictionError;
    use super::{
        depict_from_sides, mapping_index_offset, translate_item, DepictFromSidesError, ReactionSide,
    };
    use crate::depict::molecule::MoleculeDepictionError;
    #[cfg(feature = "coordgen")]
    use crate::depict::Depict;
    use crate::depict::{
        ArrowItem, AtomItem, AtomLabel, BondItem, Bounds, DashedContourItem, DepictionItem,
        DepictionReference, TextItem, WedgeItem, WedgeKind,
    };
    #[cfg(feature = "coordgen")]
    use crate::layout::{layout_molecule, MoleculeLayoutAlgorithm};
    use crate::layout::{MoleculeLayout, MoleculeLayoutError};

    #[rstest]
    fn test_translate_item_geometry() {
        let wedge = DepictionItem::Wedge(WedgeItem {
            tip: Point2D::new(1.0, 2.0),
            base: Point2D::new(3.0, 4.0),
            kind: WedgeKind::Hashed,
            references: vec![
                DepictionReference::Molecule(Entity::Bond(BondId(2))),
                DepictionReference::CorrespondencePair(4),
            ],
        });
        let contour = DepictionItem::DashedContour(DashedContourItem {
            points: vec![Point2D::new(-1.0, 0.0), Point2D::new(2.0, 3.0)],
            closed: true,
            references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
        });

        assert_eq!(
            translate_item(wedge, Point2D::new(5.0, -2.0), ReactionSide::Lhs),
            DepictionItem::Wedge(WedgeItem {
                tip: Point2D::new(6.0, 0.0),
                base: Point2D::new(8.0, 2.0),
                kind: WedgeKind::Hashed,
                references: vec![
                    DepictionReference::ReactionLhs(Entity::Bond(BondId(2))),
                    DepictionReference::CorrespondencePair(4),
                ],
            })
        );
        assert_eq!(
            translate_item(contour, Point2D::new(-3.0, 4.0), ReactionSide::Rhs),
            DepictionItem::DashedContour(DashedContourItem {
                points: vec![Point2D::new(-4.0, 4.0), Point2D::new(-1.0, 7.0)],
                closed: true,
                references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(1)))],
            })
        );
    }

    #[rstest]
    #[case::isolated(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        AtomId(0),
        Point2D::new(0.35, 0.35)
    )]
    #[case::degree_one(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap(),
        AtomId(0),
        Point2D::new(-0.4949747468305833, 0.0)
    )]
    #[case::degree_two(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap(),
        AtomId(0),
        Point2D::new(-0.35, -0.35)
    )]
    #[case::collinear(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(-1.0, 0.0),
        ])
        .unwrap(),
        AtomId(0),
        Point2D::new(0.0, 0.4949747468305833)
    )]
    #[case::tied(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#
        ),
        MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
            Point2D::new(-1.0, 0.0),
            Point2D::new(0.0, -1.0),
        ])
        .unwrap(),
        AtomId(0),
        Point2D::new(0.35, 0.35)
    )]
    fn test_mapping_index_offset(
        #[case] molecule: Molecule,
        #[case] layout: MoleculeLayout,
        #[case] atom: AtomId,
        #[case] expected: Point2D,
    ) {
        let actual = mapping_index_offset(&molecule, &layout, atom);

        assert!(approx_eq!(f64, actual.x, expected.x, epsilon = 1e-12));
        assert!(approx_eq!(f64, actual.y, expected.y, epsilon = 1e-12));
    }

    #[rstest]
    fn test_depict_from_sides() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let rhs = mol_dsl!(r#"{:atoms ["O" "C"] :bonds [[0 1 "1"]]}"#);
        let lhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let rhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let correspondence =
            Correspondence::new(vec![(AtomId(1), AtomId(0)), (AtomId(0), AtomId(1))], 2, 2)
                .unwrap();

        let depiction =
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        assert_eq!(
            depiction.items(),
            [
                DepictionItem::Bond(BondItem {
                    start: Point2D::new(-2.75, 0.0),
                    end: Point2D::new(-1.75, 0.0),
                    line_count: 1,
                    references: vec![DepictionReference::ReactionLhs(Entity::Bond(BondId(0),))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(-1.75, 0.0),
                    label: AtomLabel {
                        base: "O".to_owned(),
                        left_superscript: None,
                        right_subscript: None,
                        right_superscript: None,
                    },
                    references: vec![DepictionReference::ReactionLhs(Entity::Atom(AtomId(1)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-3.244974746830583, 0.0),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-1.2550252531694168, 0.0),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                }),
                DepictionItem::Arrow(ArrowItem {
                    start: Point2D::new(-0.75, 0.0),
                    end: Point2D::new(0.75, 0.0),
                    references: Vec::new(),
                }),
                DepictionItem::Bond(BondItem {
                    start: Point2D::new(1.75, 0.0),
                    end: Point2D::new(2.75, 0.0),
                    line_count: 1,
                    references: vec![DepictionReference::ReactionRhs(Entity::Bond(BondId(0),))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(1.75, 0.0),
                    label: AtomLabel {
                        base: "O".to_owned(),
                        left_superscript: None,
                        right_subscript: None,
                        right_superscript: None,
                    },
                    references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(3.244974746830583, 0.0),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(1.2550252531694168, 0.0),
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
                min: Point2D::new(-3.244974746830583, 0.0),
                max: Point2D::new(3.244974746830583, 0.0),
            })
        );
        assert_eq!(
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence,),
            Ok(depiction)
        );
    }

    #[rstest]
    fn test_depict_from_sides_partial() {
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

        let depiction =
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();
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
                    position: Point2D::new(-3.4, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                },
                TextItem {
                    position: Point2D::new(-1.4, 0.35),
                    text: "1".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(2))),
                        DepictionReference::CorrespondencePair(1),
                    ],
                },
                TextItem {
                    position: Point2D::new(3.244974746830583, 0.0),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                },
                TextItem {
                    position: Point2D::new(1.2550252531694168, 0.0),
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
    fn test_depict_from_sides_tetrahedral_wedge() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let layout = MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
            Point2D::new(-1.0, 0.0),
            Point2D::new(0.0, -1.0),
        ])
        .unwrap();
        let correspondence = Correspondence::new(
            (0..5)
                .map(|index| (AtomId::from(index), AtomId::from(index)))
                .collect(),
            5,
            5,
        )
        .unwrap();

        let depiction =
            depict_from_sides(&molecule, &layout, &molecule, &layout, &correspondence).unwrap();
        let wedges = depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Wedge(wedge) => Some(wedge.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(
            wedges,
            [
                WedgeItem {
                    tip: Point2D::new(-2.75, 0.0),
                    base: Point2D::new(-1.75, 0.0),
                    kind: WedgeKind::Solid,
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Bond(BondId(0))),
                        DepictionReference::ReactionLhs(Entity::StereoAtom(StereoAtomId(0))),
                    ],
                },
                WedgeItem {
                    tip: Point2D::new(2.75, 0.0),
                    base: Point2D::new(3.75, 0.0),
                    kind: WedgeKind::Solid,
                    references: vec![
                        DepictionReference::ReactionRhs(Entity::Bond(BondId(0))),
                        DepictionReference::ReactionRhs(Entity::StereoAtom(StereoAtomId(0))),
                    ],
                },
            ]
        );
    }

    #[rstest]
    fn test_depict_from_sides_projection() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "*"] :bonds []}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let lhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(10.0, 0.0)]).unwrap();
        let rhs_layout = MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap();
        let correspondence = Correspondence::new(vec![(AtomId(1), AtomId(0))], 2, 1).unwrap();

        let depiction =
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        assert_eq!(
            depiction.items(),
            [
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(-11.75, 0.0),
                    label: AtomLabel {
                        base: "C".to_owned(),
                        left_superscript: None,
                        right_subscript: None,
                        right_superscript: None,
                    },
                    references: vec![DepictionReference::ReactionLhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(-1.4, 0.35),
                    text: "0".to_owned(),
                    references: vec![
                        DepictionReference::ReactionLhs(Entity::Atom(AtomId(1))),
                        DepictionReference::CorrespondencePair(0),
                    ],
                }),
                DepictionItem::Arrow(ArrowItem {
                    start: Point2D::new(-0.75, 0.0),
                    end: Point2D::new(0.75, 0.0),
                    references: Vec::new(),
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(1.75, 0.0),
                    label: AtomLabel {
                        base: "C".to_owned(),
                        left_superscript: None,
                        right_subscript: None,
                        right_superscript: None,
                    },
                    references: vec![DepictionReference::ReactionRhs(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Text(TextItem {
                    position: Point2D::new(2.1, 0.35),
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
        DepictFromSidesError::LhsDepiction(MoleculeDepictionError::LayoutFrame(
            MoleculeLayoutError::FrameSizeMismatch {
                molecule_atom_count: 1,
                layout_atom_count: 0,
            }
        ))
    )]
    #[case::rhs_layout(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(Vec::new()).unwrap(),
        Correspondence::new(vec![(AtomId(0), AtomId(0))], 1, 1).unwrap(),
        DepictFromSidesError::RhsDepiction(MoleculeDepictionError::LayoutFrame(
            MoleculeLayoutError::FrameSizeMismatch {
                molecule_atom_count: 1,
                layout_atom_count: 0,
            }
        ))
    )]
    #[case::lhs_correspondence(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        Correspondence::new(Vec::new(), 0, 1).unwrap(),
        DepictFromSidesError::LhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 1,
            correspondence_atom_count: 0,
        }
    )]
    #[case::rhs_correspondence(
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
        Correspondence::new(Vec::new(), 1, 0).unwrap(),
        DepictFromSidesError::RhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 1,
            correspondence_atom_count: 0,
        }
    )]
    fn test_depict_from_sides_error(
        #[case] lhs_layout: MoleculeLayout,
        #[case] rhs_layout: MoleculeLayout,
        #[case] correspondence: Correspondence<AtomId>,
        #[case] expected: DepictFromSidesError,
    ) {
        let lhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);

        assert_eq!(
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence,),
            Err(expected)
        );
    }

    #[rstest]
    #[case::lhs(true)]
    #[case::rhs(false)]
    fn test_depict_from_sides_depiction_error(#[case] lhs_fails: bool) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let valid = MoleculeLayout::try_new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
            Point2D::new(-1.0, 0.0),
            Point2D::new(0.0, -1.0),
        ])
        .unwrap();
        let degenerate = MoleculeLayout::try_new(
            (0..5)
                .map(|index| Point2D::new(index as f64, 0.0))
                .collect(),
        )
        .unwrap();
        let (lhs_layout, rhs_layout, expected) = if lhs_fails {
            (
                &degenerate,
                &valid,
                DepictFromSidesError::LhsDepiction(MoleculeDepictionError::TetrahedralGeometry {
                    stereo_atom: StereoAtomId(0),
                }),
            )
        } else {
            (
                &valid,
                &degenerate,
                DepictFromSidesError::RhsDepiction(MoleculeDepictionError::TetrahedralGeometry {
                    stereo_atom: StereoAtomId(0),
                }),
            )
        };
        let correspondence = Correspondence::new(Vec::new(), 5, 5).unwrap();

        assert_eq!(
            depict_from_sides(
                &molecule,
                lhs_layout,
                &molecule,
                rhs_layout,
                &correspondence,
            ),
            Err(expected)
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
    fn test_depict_from_sides_with(#[case] algorithm: MoleculeLayoutAlgorithm) {
        let lhs = mol_dsl!(r#"{:atoms ["C" "F" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let rhs = mol_dsl!(r#"{:atoms ["O" "C" "F"] :bonds [[0 2 "1"] [2 1 "1"]]}"#);
        let correspondence =
            Correspondence::new(vec![(AtomId(2), AtomId(0)), (AtomId(0), AtomId(1))], 3, 3)
                .unwrap();
        let lhs_layout = layout_molecule(&lhs, algorithm).unwrap();
        let rhs_layout = layout_molecule(&rhs, algorithm).unwrap();
        let expected =
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        let actual = depict_from_sides_with(&lhs, &rhs, &correspondence, algorithm).unwrap();
        let index_references = actual
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Text(item) => Some(item.references.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(actual, expected);
        assert_eq!(
            index_references,
            [
                vec![
                    DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                    DepictionReference::CorrespondencePair(0),
                ],
                vec![
                    DepictionReference::ReactionLhs(Entity::Atom(AtomId(2))),
                    DepictionReference::CorrespondencePair(1),
                ],
                vec![
                    DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                    DepictionReference::CorrespondencePair(0),
                ],
                vec![
                    DepictionReference::ReactionRhs(Entity::Atom(AtomId(0))),
                    DepictionReference::CorrespondencePair(1),
                ],
            ]
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::lhs_correspondence(
        1,
        2,
        DepictFromSidesError::LhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 2,
            correspondence_atom_count: 1,
        }
    )]
    #[case::rhs_correspondence(
        2,
        1,
        DepictFromSidesError::RhsCorrespondenceFrameSizeMismatch {
            molecule_atom_count: 2,
            correspondence_atom_count: 1,
        }
    )]
    fn test_depict_from_sides_with_error(
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: DepictFromSidesError,
    ) {
        let lhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let correspondence = Correspondence::new(Vec::new(), left_count, right_count).unwrap();

        assert_eq!(
            depict_from_sides_with(
                &lhs,
                &rhs,
                &correspondence,
                MoleculeLayoutAlgorithm::CoordGen,
            ),
            Err(expected)
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict() {
        let reaction = Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        );

        assert_eq!(
            reaction.depict(),
            reaction.depict_with(MoleculeLayoutAlgorithm::CoordGen)
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_error() {
        let reaction = Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(3),
                },
            })]),
        );

        assert_eq!(
            reaction.depict(),
            reaction.depict_with(MoleculeLayoutAlgorithm::CoordGen)
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
    fn test_reaction_depict_with(#[case] algorithm: MoleculeLayoutAlgorithm) {
        let reaction = Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        );
        let span = reaction.to_reaction_span().unwrap();
        let lhs = span.lhs();
        let rhs = span.rhs();
        let correspondence = span.correspondence();
        let expected =
            depict_from_sides_with(&lhs, &rhs, correspondence.atoms(), algorithm).unwrap();

        let depiction = reaction.depict_with(algorithm).unwrap();
        let line_counts = depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Bond(item) => Some(item.line_count),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(depiction, expected);
        assert_eq!(line_counts, [1, 2]);
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_with_error() {
        let reaction = Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(3),
                },
            })]),
        );

        assert_eq!(
            reaction.depict_with(MoleculeLayoutAlgorithm::CoordGen),
            Err(ReactionDepictionError::Materialization(Contradiction))
        );
    }
}
