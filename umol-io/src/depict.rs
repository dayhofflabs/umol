//! Format-neutral molecule and reaction depiction.

pub(crate) mod molecule;
mod reaction;

pub use molecule::MoleculeDepictionError;
pub use reaction::ReactionDepictionError;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::Entity;

#[cfg(feature = "coordgen")]
use crate::layout::MoleculeLayoutAlgorithm;
use crate::svg;

/// Operational configuration for molecule and reaction depiction.
///
/// The initial configuration selects the layout algorithm. CoordGen is currently the only
/// available variant.
#[cfg(feature = "coordgen")]
pub struct DepictConfig {
    /// Algorithm used to generate the molecular layout.
    pub layout_algorithm: MoleculeLayoutAlgorithm,
}

#[cfg(feature = "coordgen")]
impl Default for DepictConfig {
    fn default() -> Self {
        Self {
            layout_algorithm: MoleculeLayoutAlgorithm::CoordGen,
        }
    }
}

/// Constructs a format-neutral depiction using default or explicitly configured operations.
#[cfg(feature = "coordgen")]
pub trait Depict {
    /// Failure produced while laying out or depicting this value.
    type Error;

    /// Constructs the depiction with [`DepictConfig::default`].
    ///
    /// # Errors
    ///
    /// Returns [`Self::Error`] when the default layout or depiction operation cannot produce the
    /// result.
    fn depict(&self) -> Result<Depiction, Self::Error> {
        self.depict_with(&DepictConfig::default())
    }

    /// Constructs the depiction with `config`.
    ///
    /// # Errors
    ///
    /// Returns [`Self::Error`] when the configured layout or depiction operation cannot produce
    /// the result.
    fn depict_with(&self, config: &DepictConfig) -> Result<Depiction, Self::Error>;
}

/// An ordered, format-neutral molecular drawing scene.
///
/// Item order is drawing order. Bounds cover item anchors and segment endpoints; they do
/// not include renderer-dependent glyph extents. A depiction is issued by a graph-IR lowering or
/// composition operation rather than assembled through a public aggregate constructor.
pub struct Depiction {
    items: Vec<DepictionItem>,
    bounds: Option<Bounds>,
}

impl Depiction {
    fn from_items(items: Vec<DepictionItem>) -> Self {
        let bounds = Bounds::from_items(&items);
        Self { items, bounds }
    }

    pub(crate) fn items(&self) -> &[DepictionItem] {
        &self.items
    }

    pub(crate) fn bounds(&self) -> Option<&Bounds> {
        self.bounds.as_ref()
    }

    /// Renders this depiction as a complete SVG document fragment.
    ///
    /// Item order is preserved. Coordinates are converted from the depiction's y-up convention to
    /// SVG's y-down convention. Molecular strokes are masked beneath estimated atom-label bounds,
    /// and structured source references are encoded in `data-umol-references` attributes.
    pub fn render_svg(&self) -> String {
        svg::render(self)
    }
}

/// One format-neutral item in a depiction.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum DepictionItem {
    /// A positioned atom label.
    Atom(AtomItem),
    /// A localized bond segment.
    Bond(BondItem),
    /// A localized stereochemical wedge.
    Wedge(WedgeItem),
    /// A dashed contour through an ordered sequence of points.
    DashedContour(DashedContourItem),
    /// Free text.
    Text(TextItem),
    /// A directed arrow.
    Arrow(ArrowItem),
}

impl DepictionItem {
    pub(crate) fn references(&self) -> &[DepictionReference] {
        match self {
            Self::Atom(item) => &item.references,
            Self::Bond(item) => &item.references,
            Self::Wedge(item) => &item.references,
            Self::DashedContour(item) => &item.references,
            Self::Text(item) => &item.references,
            Self::Arrow(item) => &item.references,
        }
    }
}

/// Typographic content of an atom label.
///
/// Each optional script is positioned relative to `base`. The carrier describes presentation
/// structure only; it does not require nonempty text or assign chemical meaning to a slot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AtomLabel {
    /// Text on the ordinary baseline.
    pub(crate) base: String,
    /// Optional superscript preceding the base.
    pub(crate) left_superscript: Option<String>,
    /// Optional subscript following the base.
    pub(crate) right_subscript: Option<String>,
    /// Optional superscript following the base and any subscript.
    pub(crate) right_superscript: Option<String>,
}

/// A positioned atom label.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AtomItem {
    /// Center of the atom label.
    pub(crate) position: Point2D,
    /// Structured display text selected by depiction construction.
    pub(crate) label: AtomLabel,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// A localized bond drawn between two points.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BondItem {
    /// First line endpoint.
    pub(crate) start: Point2D,
    /// Second line endpoint.
    pub(crate) end: Point2D,
    /// Number of parallel lines selected by depiction construction.
    pub(crate) line_count: u8,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// A stereochemical wedge between a narrow tip and a wider base.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct WedgeItem {
    /// Narrow endpoint, located at the stereocenter in an issued depiction.
    pub(crate) tip: Point2D,
    /// Center of the wide endpoint at the ligand.
    pub(crate) base: Point2D,
    /// Visible wedge treatment.
    pub(crate) kind: WedgeKind,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// Visible treatments for a stereochemical wedge.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum WedgeKind {
    /// A filled triangular wedge.
    Solid,
    /// A wedge represented by transverse hash marks.
    Hashed,
}

/// A dashed contour through an ordered sequence of scene points.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DashedContourItem {
    /// Points visited by the contour in drawing order.
    pub(crate) points: Vec<Point2D>,
    /// Whether the last point is joined back to the first.
    pub(crate) closed: bool,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// Free text placed at one scene position.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TextItem {
    /// Center of the text anchor.
    pub(crate) position: Point2D,
    /// Exact display text.
    pub(crate) text: String,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// A directed reaction or annotation arrow.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ArrowItem {
    /// Tail of the arrow.
    pub(crate) start: Point2D,
    /// Tip of the arrow.
    pub(crate) end: Point2D,
    /// Structured source references carried into rendered output.
    pub(crate) references: Vec<DepictionReference>,
}

/// A structured link from presentation data to its graph or reaction source.
///
/// Entity ids remain local to the named molecular frame. Correspondence-pair and delta values are
/// zero-based positions in their respective ordered sequences, not persistent identifiers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum DepictionReference {
    /// An entity in a standalone molecule depiction.
    Molecule(Entity),
    /// An entity on the lhs of a reaction depiction.
    ReactionLhs(Entity),
    /// An entity on the rhs of a reaction depiction.
    ReactionRhs(Entity),
    /// The zero-based displayed index of a correspondence pair.
    CorrespondencePair(u32),
}

/// Axis-aligned bounds over item anchors and segment endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Bounds {
    /// Component-wise minimum anchor coordinates.
    pub(crate) min: Point2D,
    /// Component-wise maximum anchor coordinates.
    pub(crate) max: Point2D,
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
                DepictionItem::Wedge(item) => {
                    include_point(&mut bounds, item.tip);
                    include_point(&mut bounds, item.base);
                }
                DepictionItem::DashedContour(item) => {
                    for &point in &item.points {
                        include_point(&mut bounds, point);
                    }
                }
                DepictionItem::Text(item) => include_point(&mut bounds, item.position),
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
    use roxmltree::Document;
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AromaticSystemId, AtomId, BondId, DativeBondId, Entity, MulticenterBondId,
        NoncovalentBondId, StereoAtomId, StereoBondId,
    };

    use super::*;
    use crate::svg::render;

    #[rstest]
    #[case::empty(Vec::new(), None)]
    #[case::all_item_kinds(
        vec![
            DepictionItem::Atom(AtomItem {
                position: Point2D::new(2.0, 3.0),
                label: AtomLabel {
                    base: "C".to_owned(),
                    left_superscript: None,
                    right_subscript: None,
                    right_superscript: None,
                },
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(0)))],
            }),
            DepictionItem::Bond(BondItem {
                start: Point2D::new(-1.0, 1.0),
                end: Point2D::new(4.0, -2.0),
                line_count: 2,
                references: vec![DepictionReference::Molecule(Entity::Bond(BondId(0)))],
            }),
            DepictionItem::Wedge(WedgeItem {
                tip: Point2D::new(-4.0, -5.0),
                base: Point2D::new(7.0, 2.0),
                kind: WedgeKind::Solid,
                references: vec![DepictionReference::Molecule(Entity::StereoAtom(
                    StereoAtomId(0),
                ))],
            }),
            DepictionItem::DashedContour(DashedContourItem {
                points: vec![Point2D::new(-6.0, 1.0), Point2D::new(0.0, 8.0)],
                closed: true,
                references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                    AromaticSystemId(0),
                ))],
            }),
            DepictionItem::Text(TextItem {
                position: Point2D::new(0.0, 5.0),
                text: "0".to_owned(),
                references: vec![DepictionReference::CorrespondencePair(0)],
            }),
            DepictionItem::Arrow(ArrowItem {
                start: Point2D::new(1.0, -4.0),
                end: Point2D::new(6.0, 2.0),
                references: Vec::new(),
            }),
        ],
        Some(Bounds {
            min: Point2D::new(-6.0, -5.0),
            max: Point2D::new(7.0, 8.0),
        })
    )]
    fn test_depiction_from_items(
        #[case] items: Vec<DepictionItem>,
        #[case] expected_bounds: Option<Bounds>,
    ) {
        let depiction = Depiction::from_items(items.clone());

        assert_eq!(depiction.items, items);
        assert_eq!(depiction.bounds, expected_bounds);
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
            label: AtomLabel {
                base: "N".to_owned(),
                left_superscript: None,
                right_subscript: None,
                right_superscript: None,
            },
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
            label: AtomLabel {
                base: "C".to_owned(),
                left_superscript: None,
                right_subscript: None,
                right_superscript: None,
            },
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
    #[case::wedge(
        DepictionItem::Wedge(WedgeItem {
            tip: Point2D::new(0.0, 0.0),
            base: Point2D::new(1.0, 0.0),
            kind: WedgeKind::Hashed,
            references: vec![
                DepictionReference::Molecule(Entity::Bond(BondId(2))),
                DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(1))),
            ],
        }),
        vec![
            DepictionReference::Molecule(Entity::Bond(BondId(2))),
            DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(1))),
        ]
    )]
    #[case::dashed_contour(
        DepictionItem::DashedContour(DashedContourItem {
            points: vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)],
            closed: false,
            references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(2),
            ))],
        }),
        vec![DepictionReference::Molecule(Entity::AromaticSystem(
            AromaticSystemId(2),
        ))]
    )]
    #[case::text(
        DepictionItem::Text(TextItem {
            position: Point2D::new(0.5, 1.0),
            text: "3".to_owned(),
            references: vec![DepictionReference::CorrespondencePair(3)],
        }),
        vec![DepictionReference::CorrespondencePair(3)]
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

    #[rstest]
    fn test_depiction_render_svg() {
        let text = "C<&>\"'";
        let depiction = Depiction::from_items(vec![
            DepictionItem::Atom(AtomItem {
                position: Point2D::new(1.0, 2.0),
                label: AtomLabel {
                    base: text.to_owned(),
                    left_superscript: None,
                    right_subscript: None,
                    right_superscript: None,
                },
                references: vec![
                    DepictionReference::Molecule(Entity::Atom(AtomId(0))),
                    DepictionReference::ReactionLhs(Entity::Bond(BondId(1))),
                    DepictionReference::ReactionRhs(Entity::DativeBond(DativeBondId(2))),
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(3))),
                    DepictionReference::ReactionLhs(Entity::MulticenterBond(MulticenterBondId(4))),
                    DepictionReference::ReactionRhs(Entity::NoncovalentBond(NoncovalentBondId(5))),
                    DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(6))),
                    DepictionReference::ReactionLhs(Entity::StereoBond(StereoBondId(7))),
                    DepictionReference::CorrespondencePair(8),
                ],
            }),
            DepictionItem::Bond(BondItem {
                start: Point2D::new(-2.0, -1.0),
                end: Point2D::new(2.0, -1.0),
                line_count: 2,
                references: Vec::new(),
            }),
            DepictionItem::Text(TextItem {
                position: Point2D::new(0.0, 3.0),
                text: text.to_owned(),
                references: vec![DepictionReference::CorrespondencePair(8)],
            }),
            DepictionItem::Arrow(ArrowItem {
                start: Point2D::new(-3.0, -2.0),
                end: Point2D::new(3.0, -2.0),
                references: Vec::new(),
            }),
        ]);

        let svg = depiction.render_svg();

        assert_eq!(svg, render(&depiction));

        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let groups = root
            .children()
            .filter(|child| child.has_tag_name("g"))
            .collect::<Vec<_>>();

        assert_eq!(root.attribute("viewBox"), Some("-3.5 -3.5 7 6"));
        assert_eq!(
            groups
                .iter()
                .map(|group| group.attribute("data-umol-item").unwrap())
                .collect::<Vec<_>>(),
            ["atom", "bond", "text", "arrow"]
        );
        assert_eq!(
            groups[0].attribute("data-umol-references"),
            Some(
                "molecule/atom/0 reaction-lhs/bond/1 reaction-rhs/dative-bond/2 \
                 molecule/aromatic-system/3 reaction-lhs/multicenter-bond/4 \
                 reaction-rhs/noncovalent-bond/5 molecule/stereo-atom/6 \
                 reaction-lhs/stereo-bond/7 correspondence-pair/8"
            )
        );
        assert_eq!(
            groups[0]
                .first_element_child()
                .unwrap()
                .first_element_child()
                .unwrap()
                .text(),
            Some(text)
        );
        assert_eq!(groups[1].attribute("data-umol-references"), None);
        assert_eq!(
            groups[1].attribute("mask"),
            Some("url(#umol-atom-label-mask)")
        );
        let bond_lines: Vec<_> = groups[1]
            .children()
            .filter(|child| child.is_element())
            .collect();
        assert_eq!(bond_lines.len(), 2);
        assert_eq!(bond_lines[0].attribute("x1"), Some("-2"));
        assert_eq!(bond_lines[0].attribute("y1"), Some("1.06"));
        assert_eq!(bond_lines[0].attribute("x2"), Some("2"));
        assert_eq!(bond_lines[0].attribute("y2"), Some("1.06"));
        assert_eq!(bond_lines[1].attribute("x1"), Some("-2"));
        assert_eq!(bond_lines[1].attribute("y1"), Some("0.94"));
        assert_eq!(bond_lines[1].attribute("x2"), Some("2"));
        assert_eq!(bond_lines[1].attribute("y2"), Some("0.94"));
        assert_eq!(groups[2].first_element_child().unwrap().text(), Some(text));
        assert_eq!(
            groups[3]
                .children()
                .filter(|child| child.is_element())
                .map(|child| child.tag_name().name())
                .collect::<Vec<_>>(),
            ["line", "polygon"]
        );
        assert_eq!(groups[3].attribute("data-umol-references"), None);
        assert_eq!(svg.matches("C&lt;&amp;&gt;&quot;&apos;").count(), 2);
    }
}
