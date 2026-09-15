//! Reaction-owned depiction and indexed-side composition.

use std::any::Any;
use std::f64::consts::{PI, SQRT_2, TAU};

use thiserror::Error;
use umol_geometric_core::Point2D;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomId, Entity, Molecule};
#[cfg(feature = "coordgen")]
use umol_graph_ir::ir::{Contradiction, Reaction};
use umol_utils::error::UmolError;

use super::{molecule, ArrowItem, Depiction, DepictionItem, DepictionReference, TextItem};
#[cfg(feature = "coordgen")]
use super::{Depict, DepictConfig};
use crate::depict::molecule::MoleculeDepictionError;
#[cfg(feature = "coordgen")]
use crate::layout::{check_arrow, ReactionLayoutError};
use crate::layout::{MoleculeLayout, ReactionLayout};

const MAP_INDEX_COMPONENT_OFFSET: f64 = 0.35;
const MAP_INDEX_DISTANCE: f64 = MAP_INDEX_COMPONENT_OFFSET * SQRT_2;
const MAP_INDEX_FALLBACK_OFFSET: Point2D =
    Point2D::new(MAP_INDEX_COMPONENT_OFFSET, MAP_INDEX_COMPONENT_OFFSET);

fn compose_sides(
    lhs: &Molecule,
    lhs_depiction: Depiction,
    rhs: &Molecule,
    rhs_depiction: Depiction,
    layout: &ReactionLayout,
    atom_correspondence: &Correspondence<AtomId>,
) -> Depiction {
    let mut items = lhs_depiction
        .items
        .into_iter()
        .map(|item| reference_side(item, ReactionSide::Lhs))
        .collect::<Vec<_>>();

    items.extend(index_items(
        lhs,
        layout.lhs(),
        atom_correspondence
            .matched_pairs()
            .iter()
            .enumerate()
            .map(|(index, &(left, _))| (index, left)),
        ReactionSide::Lhs,
    ));
    items.push(DepictionItem::Arrow(ArrowItem {
        start: layout.arrow_start(),
        end: layout.arrow_end(),
        references: Vec::new(),
    }));
    items.extend(
        rhs_depiction
            .items
            .into_iter()
            .map(|item| reference_side(item, ReactionSide::Rhs)),
    );
    items.extend(index_items(
        rhs,
        layout.rhs(),
        atom_correspondence
            .matched_pairs()
            .iter()
            .enumerate()
            .map(|(index, &(_, right))| (index, right)),
        ReactionSide::Rhs,
    ));

    Depiction::from_items(items)
}

#[cfg(feature = "coordgen")]
impl Depict for Reaction {
    type Layout = ReactionLayout;
    type Error = ReactionDepictionError;

    fn layout_with(&self, config: &DepictConfig) -> Result<Self::Layout, Self::Error> {
        let span = self
            .to_reaction_span()
            .map_err(ReactionDepictionError::Materialization)?;
        let lhs = span
            .lhs()
            .layout_with(config)
            .map_err(ReactionDepictionError::LhsDepiction)?;
        let rhs = span
            .rhs()
            .layout_with(config)
            .map_err(ReactionDepictionError::RhsDepiction)?;
        ReactionLayout::arrange(lhs, rhs).map_err(ReactionDepictionError::Layout)
    }

    fn verify_layout(&self, layout: &Self::Layout) -> Result<(), Self::Error> {
        let span = self
            .to_reaction_span()
            .map_err(ReactionDepictionError::Materialization)?;
        verify_sides(&span.lhs(), &span.rhs(), layout)
    }

    fn depict_layout(&self, layout: &Self::Layout) -> Result<Depiction, Self::Error> {
        let span = self
            .to_reaction_span()
            .map_err(ReactionDepictionError::Materialization)?;
        let lhs = span.lhs();
        let rhs = span.rhs();
        verify_sides(&lhs, &rhs, layout)?;
        let lhs_depiction =
            molecule::depict(&lhs, layout.lhs()).map_err(ReactionDepictionError::LhsDepiction)?;
        let rhs_depiction =
            molecule::depict(&rhs, layout.rhs()).map_err(ReactionDepictionError::RhsDepiction)?;

        Ok(compose_sides(
            &lhs,
            lhs_depiction,
            &rhs,
            rhs_depiction,
            layout,
            span.correspondence().atoms(),
        ))
    }
}

#[cfg(feature = "coordgen")]
fn verify_sides(
    lhs: &Molecule,
    rhs: &Molecule,
    layout: &ReactionLayout,
) -> Result<(), ReactionDepictionError> {
    lhs.verify_layout(layout.lhs())
        .map_err(ReactionDepictionError::LhsDepiction)?;
    rhs.verify_layout(layout.rhs())
        .map_err(ReactionDepictionError::RhsDepiction)?;
    check_arrow(layout.arrow_start(), layout.arrow_end()).map_err(ReactionDepictionError::Layout)
}

/// Failures while laying out, verifying, or depicting a [`Reaction`].
#[cfg(feature = "coordgen")]
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ReactionDepictionError {
    /// The reaction deltas could not be materialized into a two-sided reaction span.
    #[error("reaction materialization: {0}")]
    Materialization(#[source] Contradiction),
    /// Layout, verification, or depiction of the materialized left-hand side failed.
    #[error("lhs depiction: {0}")]
    LhsDepiction(#[source] MoleculeDepictionError),
    /// Layout, verification, or depiction of the materialized right-hand side failed.
    #[error("rhs depiction: {0}")]
    RhsDepiction(#[source] MoleculeDepictionError),
    /// The reaction layout could not be arranged, or its arrow cannot be drawn.
    #[error("reaction layout: {0}")]
    Layout(#[source] ReactionLayoutError),
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

fn reference_side(mut item: DepictionItem, side: ReactionSide) -> DepictionItem {
    let references = match &mut item {
        DepictionItem::Atom(item) => &mut item.references,
        DepictionItem::Bond(item) => &mut item.references,
        DepictionItem::Wedge(item) => &mut item.references,
        DepictionItem::DashedContour(item) => &mut item.references,
        DepictionItem::Text(item) => &mut item.references,
        DepictionItem::Arrow(item) => &mut item.references,
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
                position: translate(position, mapping_index_offset(molecule, layout, atom)),
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
    use umol_graph_ir::ir::{AtomId, BondId, Entity, Molecule};
    #[cfg(feature = "coordgen")]
    use umol_graph_ir::ir::{
        BondDelta, BondFieldChange, Contradiction, Delta, Deltas, NumForm, Reaction, StereoAtomId,
    };
    use umol_graph_ir::mol_dsl;

    #[cfg(feature = "coordgen")]
    use super::ReactionDepictionError;
    use super::{
        mapping_index_offset, reference_side, ArrowItem, DepictionItem, DepictionReference,
        ReactionSide,
    };
    #[cfg(feature = "coordgen")]
    use crate::depict::molecule::MoleculeDepictionError;
    use crate::depict::{DashedContourItem, WedgeItem, WedgeKind};
    #[cfg(feature = "coordgen")]
    use crate::depict::{Depict, DepictConfig};
    use crate::layout::MoleculeLayout;
    #[cfg(feature = "coordgen")]
    use crate::layout::{MoleculeLayoutError, ReactionLayout, ReactionLayoutError};

    #[rstest]
    fn test_reference_side() {
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
            reference_side(wedge, ReactionSide::Lhs),
            DepictionItem::Wedge(WedgeItem {
                tip: Point2D::new(1.0, 2.0),
                base: Point2D::new(3.0, 4.0),
                kind: WedgeKind::Hashed,
                references: vec![
                    DepictionReference::ReactionLhs(Entity::Bond(BondId(2))),
                    DepictionReference::CorrespondencePair(4),
                ],
            })
        );
        assert_eq!(
            reference_side(contour, ReactionSide::Rhs),
            DepictionItem::DashedContour(DashedContourItem {
                points: vec![Point2D::new(-1.0, 0.0), Point2D::new(2.0, 3.0)],
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

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict() {
        let reaction = bond_order_reaction();

        assert_eq!(
            reaction.depict().unwrap().render_svg(),
            reaction
                .depict_with(&DepictConfig::default())
                .unwrap()
                .render_svg()
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_items() {
        let depiction = bond_order_reaction().depict().unwrap();
        let bonds = depiction
            .items
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Bond(bond) => Some((bond.line_count, bond.references.clone())),
                _ => None,
            })
            .collect::<Vec<_>>();
        let mapping_references = depiction
            .items
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Text(text)
                    if text.references.iter().any(|reference| {
                        matches!(reference, DepictionReference::CorrespondencePair(_))
                    }) =>
                {
                    Some(text.references.clone())
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        let arrows = depiction
            .items
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Arrow(arrow) => Some(arrow.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(
            bonds,
            [
                (
                    1,
                    vec![DepictionReference::ReactionLhs(Entity::Bond(BondId(0)))]
                ),
                (
                    2,
                    vec![DepictionReference::ReactionRhs(Entity::Bond(BondId(0)))]
                ),
            ]
        );
        assert_eq!(
            mapping_references,
            [
                vec![
                    DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
                    DepictionReference::CorrespondencePair(0),
                ],
                vec![
                    DepictionReference::ReactionLhs(Entity::Atom(AtomId(1))),
                    DepictionReference::CorrespondencePair(1),
                ],
                vec![
                    DepictionReference::ReactionRhs(Entity::Atom(AtomId(0))),
                    DepictionReference::CorrespondencePair(0),
                ],
                vec![
                    DepictionReference::ReactionRhs(Entity::Atom(AtomId(1))),
                    DepictionReference::CorrespondencePair(1),
                ],
            ]
        );
        assert_eq!(
            arrows,
            [ArrowItem {
                start: Point2D::new(-0.75, 0.0),
                end: Point2D::new(0.75, 0.0),
                references: Vec::new(),
            }]
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_stereo() {
        let reaction = Reaction::new(
            mol_dsl!(
                r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
            ),
            Deltas::new(),
        );

        let depiction = reaction.depict().unwrap();
        let wedges = depiction
            .items
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Wedge(wedge) => Some(wedge),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(wedges.len(), 2);
        assert_eq!(wedges[0].kind, wedges[1].kind);
        assert_eq!(
            wedges[0].references,
            [
                DepictionReference::ReactionLhs(Entity::Bond(BondId(0))),
                DepictionReference::ReactionLhs(Entity::StereoAtom(StereoAtomId(0))),
            ]
        );
        assert_eq!(
            wedges[1].references,
            [
                DepictionReference::ReactionRhs(Entity::Bond(BondId(0))),
                DepictionReference::ReactionRhs(Entity::StereoAtom(StereoAtomId(0))),
            ]
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_materialization_error() {
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
            reaction.depict().err(),
            Some(ReactionDepictionError::Materialization(Contradiction))
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::lhs(2, 1, true)]
    #[case::rhs(1, 2, false)]
    fn test_reaction_depict_side_error(
        #[case] old_order: i64,
        #[case] new_order: i64,
        #[case] lhs: bool,
    ) {
        let reaction = Reaction::new(
            mol_dsl!(&format!(
                r#"{{:atoms ["C" "F" "Cl" "Br" "I"]
                         :bonds [[0 1 "{old_order}"] [0 2 "2"] [0 3 "2"] [0 4 "2"]]
                         :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "Th0"}}]}}"#
            )),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(old_order),
                    new: NumForm::Lit(new_order),
                },
            })]),
        );
        let expected = MoleculeDepictionError::TetrahedralGeometry {
            stereo_atom: StereoAtomId(0),
        };

        assert_eq!(
            reaction.depict().err(),
            Some(if lhs {
                ReactionDepictionError::LhsDepiction(expected)
            } else {
                ReactionDepictionError::RhsDepiction(expected)
            })
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_layout_equals_depict() {
        for reaction in [
            bond_order_reaction(),
            stereo_reaction(),
            aromatic_reaction(),
        ] {
            let layout = reaction.layout().unwrap();
            let supplied = reaction.depict_layout(&layout).unwrap();
            let generated = reaction.depict().unwrap();

            assert_eq!(supplied.items, generated.items);
            assert_eq!(supplied.render_svg(), generated.render_svg());
        }
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_layout_arrow() {
        let reaction = bond_order_reaction();
        let layout = reaction.layout().unwrap();

        assert_eq!(layout.arrow_start(), Point2D::new(-0.75, 0.0));
        assert_eq!(layout.arrow_end(), Point2D::new(0.75, 0.0));
        assert_eq!(reaction.verify_layout(&layout), Ok(()));
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_layout_moves_arrow() {
        let reaction = bond_order_reaction();
        let mut layout = reaction.layout().unwrap();
        layout
            .set_arrow(Point2D::new(-1.0, 0.5), Point2D::new(1.0, 0.5))
            .unwrap();

        let depiction = reaction.depict_layout(&layout).unwrap();
        let arrows = depiction
            .items
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Arrow(arrow) => Some(arrow.clone()),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(
            arrows,
            [ArrowItem {
                start: Point2D::new(-1.0, 0.5),
                end: Point2D::new(1.0, 0.5),
                references: Vec::new(),
            }]
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::lhs(true)]
    #[case::rhs(false)]
    fn test_reaction_verify_layout_side_error(#[case] lhs: bool) {
        let reaction = bond_order_reaction();
        let mut layout = reaction.layout().unwrap();
        let side = if lhs {
            layout.lhs_mut()
        } else {
            layout.rhs_mut()
        };
        let coincident = *side.position(AtomId(1)).unwrap();
        side.set_position(AtomId(0), coincident).unwrap();
        let expected = MoleculeDepictionError::NonFiniteGeometry {
            entity: Entity::Bond(BondId(0)),
        };

        assert_eq!(
            reaction.verify_layout(&layout),
            Err(if lhs {
                ReactionDepictionError::LhsDepiction(expected.clone())
            } else {
                ReactionDepictionError::RhsDepiction(expected.clone())
            })
        );
        assert_eq!(
            reaction.depict_layout(&layout).err(),
            Some(if lhs {
                ReactionDepictionError::LhsDepiction(expected)
            } else {
                ReactionDepictionError::RhsDepiction(expected)
            })
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_depict_layout_frame_error() {
        let reaction = bond_order_reaction();
        let generated = reaction.layout().unwrap();
        let layout = ReactionLayout::try_new(
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0)]).unwrap(),
            generated.rhs().clone(),
            generated.arrow_start(),
            generated.arrow_end(),
        )
        .unwrap();

        assert_eq!(
            reaction.depict_layout(&layout).err(),
            Some(ReactionDepictionError::LhsDepiction(
                MoleculeDepictionError::LayoutFrame(MoleculeLayoutError::FrameSizeMismatch {
                    molecule_atom_count: 2,
                    layout_atom_count: 1,
                })
            ))
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_layout_materialization_error() {
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
        let layout = ReactionLayout::arrange(
            MoleculeLayout::try_new(Vec::new()).unwrap(),
            MoleculeLayout::try_new(Vec::new()).unwrap(),
        )
        .unwrap();

        assert_eq!(
            reaction.layout().err(),
            Some(ReactionDepictionError::Materialization(Contradiction))
        );
        assert_eq!(
            reaction.verify_layout(&layout),
            Err(ReactionDepictionError::Materialization(Contradiction))
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_reaction_layout_error_wraps_arrow() {
        let error = ReactionDepictionError::Layout(ReactionLayoutError::DegenerateArrow {
            position: Point2D::new(0.0, 0.0),
        });

        assert_eq!(
            error.to_string(),
            "reaction layout: reaction arrow starts and ends at Point2D { x: 0.0, y: 0.0 }"
        );
    }

    #[cfg(feature = "coordgen")]
    fn stereo_reaction() -> Reaction {
        Reaction::new(
            mol_dsl!(
                r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
            ),
            Deltas::new(),
        )
    }

    #[cfg(feature = "coordgen")]
    fn aromatic_reaction() -> Reaction {
        Reaction::new(
            mol_dsl!(
                r#"{:atoms ["C" "C" "C" "C" "C" "C" "O#h1" "N#h2"]
                    :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                            [0 6 "1"] [3 7 "1"]]
                    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#c+"}]}"#
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(6),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        )
    }

    #[cfg(feature = "coordgen")]
    fn bond_order_reaction() -> Reaction {
        Reaction::new(
            mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        )
    }
}
