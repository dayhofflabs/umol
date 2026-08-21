//! Format-neutral drawing items for molecular and reaction depictions.

pub mod molecule;

use umol_graph_ir::ir::Entity;

use crate::layout::Point2D;

/// An ordered, format-neutral molecular drawing scene.
///
/// Item order is drawing order. Bounds cover item anchors and segment endpoints; they do
/// not include renderer-dependent glyph extents. A depiction is issued by a graph-IR lowering or
/// composition operation rather than assembled through a public aggregate constructor.
#[derive(Clone, Debug, PartialEq)]
pub struct Depiction {
    items: Vec<DepictionItem>,
    bounds: Option<Bounds>,
}

impl Depiction {
    fn from_items(items: Vec<DepictionItem>) -> Self {
        let bounds = Bounds::from_items(&items);
        Self { items, bounds }
    }

    /// Drawing items in their rendering order.
    pub fn items(&self) -> &[DepictionItem] {
        &self.items
    }

    /// Anchor and endpoint bounds, or `None` when the scene has no items.
    pub fn bounds(&self) -> Option<&Bounds> {
        self.bounds.as_ref()
    }
}

/// One format-neutral item in a depiction.
#[derive(Clone, Debug, PartialEq)]
pub enum DepictionItem {
    /// A positioned atom label.
    Atom(AtomItem),
    /// A localized bond segment.
    Bond(BondItem),
    /// Free text.
    Text(TextItem),
    /// An aromatic or stereo mark.
    Marker(MarkerItem),
    /// A directed arrow.
    Arrow(ArrowItem),
}

impl DepictionItem {
    /// Structured links from this item to the graph or reaction data it depicts.
    pub fn references(&self) -> &[DepictionReference] {
        match self {
            Self::Atom(item) => &item.references,
            Self::Bond(item) => &item.references,
            Self::Text(item) => &item.references,
            Self::Marker(item) => &item.references,
            Self::Arrow(item) => &item.references,
        }
    }
}

/// A positioned atom label.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomItem {
    /// Center of the atom label.
    pub position: Point2D,
    /// Display text selected by depiction construction.
    pub label: String,
    /// Structured source references carried into rendered output.
    pub references: Vec<DepictionReference>,
}

/// A localized bond drawn between two points.
#[derive(Clone, Debug, PartialEq)]
pub struct BondItem {
    /// First line endpoint.
    pub start: Point2D,
    /// Second line endpoint.
    pub end: Point2D,
    /// Number of parallel lines selected by depiction construction.
    pub line_count: u8,
    /// Structured source references carried into rendered output.
    pub references: Vec<DepictionReference>,
}

/// Free text placed at one scene position.
#[derive(Clone, Debug, PartialEq)]
pub struct TextItem {
    /// Center of the text anchor.
    pub position: Point2D,
    /// Exact display text.
    pub text: String,
    /// Structured source references carried into rendered output.
    pub references: Vec<DepictionReference>,
}

/// A small semantic mark placed at one scene position.
#[derive(Clone, Debug, PartialEq)]
pub struct MarkerItem {
    /// Center of the marker.
    pub position: Point2D,
    /// Meaning of the visible marker.
    pub kind: MarkerKind,
    /// Structured source references carried into rendered output.
    pub references: Vec<DepictionReference>,
}

/// Marker meanings supported by the first depiction projection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MarkerKind {
    /// Marks aromatic participation.
    Aromatic,
    /// Marks a stereo site without prescribing wedge or hatch geometry.
    Stereo,
}

/// A directed reaction or annotation arrow.
#[derive(Clone, Debug, PartialEq)]
pub struct ArrowItem {
    /// Tail of the arrow.
    pub start: Point2D,
    /// Tip of the arrow.
    pub end: Point2D,
    /// Structured source references carried into rendered output.
    pub references: Vec<DepictionReference>,
}

/// A structured link from presentation data to its graph or reaction source.
///
/// Entity ids remain local to the named molecular frame. Correspondence-pair and delta values are
/// zero-based positions in their respective ordered sequences, not persistent identifiers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DepictionReference {
    /// An entity in a standalone molecule depiction.
    Molecule(Entity),
    /// An entity on the lhs of a reaction depiction.
    ReactionLhs(Entity),
    /// An entity on the rhs of a reaction depiction.
    ReactionRhs(Entity),
    /// The zero-based displayed index of a correspondence pair.
    CorrespondencePair(u32),
    /// The zero-based ordinal of a reaction delta.
    Delta(u32),
}

/// Axis-aligned bounds over item anchors and segment endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds {
    /// Component-wise minimum anchor coordinates.
    pub min: Point2D,
    /// Component-wise maximum anchor coordinates.
    pub max: Point2D,
}

impl Bounds {
    fn from_items(items: &[DepictionItem]) -> Option<Self> {
        let mut bounds = None;
        for item in items {
            match item {
                DepictionItem::Atom(item) => include_point(&mut bounds, item.position),
                DepictionItem::Bond(item) => {
                    include_point(&mut bounds, item.start);
                    include_point(&mut bounds, item.end);
                }
                DepictionItem::Text(item) => include_point(&mut bounds, item.position),
                DepictionItem::Marker(item) => include_point(&mut bounds, item.position),
                DepictionItem::Arrow(item) => {
                    include_point(&mut bounds, item.start);
                    include_point(&mut bounds, item.end);
                }
            }
        }
        bounds
    }

    fn include(&mut self, point: Point2D) {
        self.min.x = self.min.x.min(point.x);
        self.min.y = self.min.y.min(point.y);
        self.max.x = self.max.x.max(point.x);
        self.max.y = self.max.y.max(point.y);
    }
}

fn include_point(bounds: &mut Option<Bounds>, point: Point2D) {
    match bounds {
        Some(bounds) => bounds.include(point),
        None => {
            *bounds = Some(Bounds {
                min: point,
                max: point,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, Entity, StereoAtomId};

    use super::*;

    #[rstest]
    #[case::empty(Vec::new(), None)]
    #[case::all_item_kinds(
        vec![
            DepictionItem::Atom(AtomItem {
                position: Point2D::new(2.0, 3.0),
                label: "C".to_owned(),
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(0)))],
            }),
            DepictionItem::Bond(BondItem {
                start: Point2D::new(-1.0, 1.0),
                end: Point2D::new(4.0, -2.0),
                line_count: 2,
                references: vec![DepictionReference::Molecule(Entity::Bond(BondId(0)))],
            }),
            DepictionItem::Text(TextItem {
                position: Point2D::new(0.0, 5.0),
                text: "0".to_owned(),
                references: vec![DepictionReference::CorrespondencePair(0)],
            }),
            DepictionItem::Marker(MarkerItem {
                position: Point2D::new(-3.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                    AromaticSystemId(0),
                ))],
            }),
            DepictionItem::Arrow(ArrowItem {
                start: Point2D::new(1.0, -4.0),
                end: Point2D::new(6.0, 2.0),
                references: Vec::new(),
            }),
        ],
        Some(Bounds {
            min: Point2D::new(-3.0, -4.0),
            max: Point2D::new(6.0, 5.0),
        })
    )]
    fn test_depiction_from_items(
        #[case] items: Vec<DepictionItem>,
        #[case] expected_bounds: Option<Bounds>,
    ) {
        let depiction = Depiction::from_items(items.clone());

        assert_eq!(
            depiction,
            Depiction {
                items,
                bounds: expected_bounds,
            }
        );
    }

    #[rstest]
    #[case::empty(Vec::new())]
    #[case::ordered(vec![
        DepictionItem::Text(TextItem {
            position: Point2D::new(1.0, 2.0),
            text: "first".to_owned(),
            references: Vec::new(),
        }),
        DepictionItem::Text(TextItem {
            position: Point2D::new(3.0, 4.0),
            text: "second".to_owned(),
            references: Vec::new(),
        }),
    ])]
    fn test_depiction_items(#[case] items: Vec<DepictionItem>) {
        let depiction = Depiction::from_items(items.clone());

        assert_eq!(depiction.items(), items);
    }

    #[rstest]
    #[case::empty(Vec::new(), None)]
    #[case::present(
        vec![DepictionItem::Atom(AtomItem {
            position: Point2D::new(-1.0, 2.0),
            label: "N".to_owned(),
            references: Vec::new(),
        })],
        Some(Bounds {
            min: Point2D::new(-1.0, 2.0),
            max: Point2D::new(-1.0, 2.0),
        })
    )]
    fn test_depiction_bounds(#[case] items: Vec<DepictionItem>, #[case] expected: Option<Bounds>) {
        let depiction = Depiction::from_items(items);

        assert_eq!(depiction.bounds().copied(), expected);
    }

    #[rstest]
    #[case::atom(
        DepictionItem::Atom(AtomItem {
            position: Point2D::new(0.0, 0.0),
            label: "C".to_owned(),
            references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
        }),
        vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))]
    )]
    #[case::bond(
        DepictionItem::Bond(BondItem {
            start: Point2D::new(0.0, 0.0),
            end: Point2D::new(1.0, 0.0),
            line_count: 1,
            references: vec![DepictionReference::ReactionLhs(Entity::Bond(BondId(2)))],
        }),
        vec![DepictionReference::ReactionLhs(Entity::Bond(BondId(2)))]
    )]
    #[case::text(
        DepictionItem::Text(TextItem {
            position: Point2D::new(0.5, 1.0),
            text: "3".to_owned(),
            references: vec![DepictionReference::CorrespondencePair(3)],
        }),
        vec![DepictionReference::CorrespondencePair(3)]
    )]
    #[case::marker(
        DepictionItem::Marker(MarkerItem {
            position: Point2D::new(2.0, 1.0),
            kind: MarkerKind::Stereo,
            references: vec![
                DepictionReference::ReactionRhs(Entity::StereoAtom(StereoAtomId(0))),
                DepictionReference::Delta(4),
            ],
        }),
        vec![
            DepictionReference::ReactionRhs(Entity::StereoAtom(StereoAtomId(0))),
            DepictionReference::Delta(4),
        ]
    )]
    #[case::arrow(
        DepictionItem::Arrow(ArrowItem {
            start: Point2D::new(-1.0, 0.0),
            end: Point2D::new(1.0, 0.0),
            references: Vec::new(),
        }),
        Vec::new()
    )]
    fn test_depiction_item_references(
        #[case] item: DepictionItem,
        #[case] expected: Vec<DepictionReference>,
    ) {
        assert_eq!(item.references(), expected);
    }
}
