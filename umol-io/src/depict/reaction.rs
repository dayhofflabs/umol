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
use crate::layout::layout_molecule;
use crate::layout::MoleculeLayout;

const ARROW_HALF_LENGTH: f64 = 0.75;
const SIDE_ARROW_GAP: f64 = 1.0;
const MAP_INDEX_COMPONENT_OFFSET: f64 = 0.35;
const MAP_INDEX_DISTANCE: f64 = MAP_INDEX_COMPONENT_OFFSET * SQRT_2;
const MAP_INDEX_FALLBACK_OFFSET: Point2D =
    Point2D::new(MAP_INDEX_COMPONENT_OFFSET, MAP_INDEX_COMPONENT_OFFSET);

fn compose_sides(
    lhs: &Molecule,
    lhs_layout: &MoleculeLayout,
    lhs_depiction: Depiction,
    rhs: &Molecule,
    rhs_layout: &MoleculeLayout,
    rhs_depiction: Depiction,
    atom_correspondence: &Correspondence<AtomId>,
) -> Depiction {
    let lhs_offset = side_offset(lhs_layout, ReactionSide::Lhs);
    let rhs_offset = side_offset(rhs_layout, ReactionSide::Rhs);
    let mut items = lhs_depiction
        .items
        .into_iter()
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
            .items
            .into_iter()
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

    Depiction::from_items(items)
}

#[cfg(feature = "coordgen")]
impl Depict for Reaction {
    type Error = ReactionDepictionError;

    fn depict_with(&self, config: &DepictConfig) -> Result<Depiction, Self::Error> {
        let span = self
            .to_reaction_span()
            .map_err(ReactionDepictionError::Materialization)?;
        let lhs = span.lhs();
        let rhs = span.rhs();
        let correspondence = span.correspondence();
        let lhs_layout = layout_molecule(&lhs, config.layout_algorithm)
            .map_err(MoleculeDepictionError::Layout)
            .map_err(ReactionDepictionError::LhsDepiction)?;
        let rhs_layout = layout_molecule(&rhs, config.layout_algorithm)
            .map_err(MoleculeDepictionError::Layout)
            .map_err(ReactionDepictionError::RhsDepiction)?;
        let lhs_depiction =
            molecule::depict(&lhs, &lhs_layout).map_err(ReactionDepictionError::LhsDepiction)?;
        let rhs_depiction =
            molecule::depict(&rhs, &rhs_layout).map_err(ReactionDepictionError::RhsDepiction)?;

        Ok(compose_sides(
            &lhs,
            &lhs_layout,
            lhs_depiction,
            &rhs,
            &rhs_layout,
            rhs_depiction,
            correspondence.atoms(),
        ))
    }
}

/// Failures while depicting a [`Reaction`].
#[cfg(feature = "coordgen")]
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ReactionDepictionError {
    /// The reaction deltas could not be materialized into a two-sided reaction span.
    #[error("reaction materialization: {0}")]
    Materialization(#[source] Contradiction),
    /// Layout or depiction of the materialized left-hand side failed.
    #[error("lhs depiction: {0}")]
    LhsDepiction(#[source] MoleculeDepictionError),
    /// Layout or depiction of the materialized right-hand side failed.
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
    use umol_graph_ir::ir::{AtomId, BondId, Entity, Molecule};
    #[cfg(feature = "coordgen")]
    use umol_graph_ir::ir::{
        BondDelta, BondFieldChange, Contradiction, Delta, Deltas, NumForm, Reaction, StereoAtomId,
    };
    use umol_graph_ir::mol_dsl;

    #[cfg(feature = "coordgen")]
    use super::ReactionDepictionError;
    use super::{
        mapping_index_offset, translate_item, ArrowItem, DepictionItem, DepictionReference,
        ReactionSide,
    };
    #[cfg(feature = "coordgen")]
    use crate::depict::molecule::MoleculeDepictionError;
    use crate::depict::{DashedContourItem, WedgeItem, WedgeKind};
    #[cfg(feature = "coordgen")]
    use crate::depict::{Depict, DepictConfig};
    use crate::layout::MoleculeLayout;

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
