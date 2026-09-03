//! Deterministic SVG rendering of format-neutral depictions.

use std::fmt::Write;

use umol_geometric_core::Point2D;
use umol_graph_ir::ir::Entity;

use crate::depiction::{
    ArrowItem, AtomItem, BondItem, Bounds, DashedContourItem, Depiction, DepictionItem,
    DepictionReference, MarkerItem, MarkerKind, TextItem, WedgeItem, WedgeKind,
};

const VIEW_MARGIN: f64 = 0.5;
const BOND_GAP: f64 = 0.12;
const BOND_WIDTH: f64 = 0.06;
const WEDGE_BASE_HALF_WIDTH: f64 = 0.1;
const WEDGE_HASH_COUNT: usize = 5;
const DASHED_CONTOUR_WIDTH: f64 = 0.04;
const DASHED_CONTOUR_PATTERN: &str = "0.12 0.1";
const TEXT_SIZE: f64 = 0.45;
const ATOM_LABEL_MASK_ID: &str = "umol-atom-label-mask";
const ATOM_LABEL_MASK_EXPANSION: f64 = 0.3;
const AROMATIC_MARKER_RADIUS: f64 = 0.16;
const STEREO_MARKER_RADIUS: f64 = 0.19;
const MARKER_WIDTH: f64 = 0.04;
const ARROW_HEAD_LENGTH: f64 = 0.24;
const ARROW_HEAD_HALF_WIDTH: f64 = 0.11;

/// Renders a [`Depiction`] as a complete SVG document fragment.
///
/// Depiction item order is preserved among item groups. When atom labels are present, one leading
/// definition masks molecular strokes beneath duplicates of the visible glyphs without painting a
/// page-background color. Coordinates are converted from the depiction's y-up convention to SVG's
/// y-down convention. The view box extends the depiction's anchor bounds by half a nominal bond
/// length on each side; an empty depiction uses a centered one-by-one view box. Structured
/// references are encoded in the `data-umol-references` attribute as ordered, space-separated path
/// tokens: molecular entity references use `molecule`, `reaction-lhs`, or `reaction-rhs` followed
/// by entity kind and id; `correspondence-pair` and `delta` references use their kind followed by
/// their zero-based position.
pub fn render(depiction: &Depiction) -> String {
    let mut output = String::new();
    let view_box = SvgViewBox::from_bounds(depiction.bounds());

    output.push_str(r#"<svg xmlns="http://www.w3.org/2000/svg" class="umol-depiction" viewBox=""#);
    write_number(&mut output, view_box.x);
    output.push(' ');
    write_number(&mut output, view_box.y);
    output.push(' ');
    write_number(&mut output, view_box.width);
    output.push(' ');
    write_number(&mut output, view_box.height);
    output.push_str("\">\n");

    let has_atom_mask = depiction
        .items()
        .iter()
        .any(|item| matches!(item, DepictionItem::Atom(_)));
    if has_atom_mask {
        render_atom_mask(&mut output, depiction, view_box);
    }

    for item in depiction.items() {
        render_item(
            &mut output,
            item,
            has_atom_mask && item_uses_atom_mask(item),
        );
    }

    output.push_str("</svg>");
    output
}

#[derive(Clone, Copy)]
struct SvgViewBox {
    x: f64,
    y: f64,
    width: f64,
    height: f64,
}

impl SvgViewBox {
    fn from_bounds(bounds: Option<&Bounds>) -> Self {
        match bounds {
            Some(bounds) => Self {
                x: bounds.min.x - VIEW_MARGIN,
                y: -bounds.max.y - VIEW_MARGIN,
                width: bounds.max.x - bounds.min.x + 2.0 * VIEW_MARGIN,
                height: bounds.max.y - bounds.min.y + 2.0 * VIEW_MARGIN,
            },
            None => Self {
                x: -VIEW_MARGIN,
                y: -VIEW_MARGIN,
                width: 2.0 * VIEW_MARGIN,
                height: 2.0 * VIEW_MARGIN,
            },
        }
    }
}

fn render_atom_mask(output: &mut String, depiction: &Depiction, view_box: SvgViewBox) {
    output.push_str(
        r#"<defs><mask id="umol-atom-label-mask" maskUnits="userSpaceOnUse" maskContentUnits="userSpaceOnUse" mask-type="luminance" x=""#,
    );
    write_number(output, view_box.x);
    output.push_str(r#"" y=""#);
    write_number(output, view_box.y);
    output.push_str(r#"" width=""#);
    write_number(output, view_box.width);
    output.push_str(r#"" height=""#);
    write_number(output, view_box.height);
    output.push_str(r#""><rect x=""#);
    write_number(output, view_box.x);
    output.push_str(r#"" y=""#);
    write_number(output, view_box.y);
    output.push_str(r#"" width=""#);
    write_number(output, view_box.width);
    output.push_str(r#"" height=""#);
    write_number(output, view_box.height);
    output.push_str(r#"" fill="white"/>"#);
    for item in depiction.items() {
        if let DepictionItem::Atom(atom) = item {
            render_atom_glyph(
                output,
                "umol-atom-mask",
                atom,
                "black",
                Some(2.0 * ATOM_LABEL_MASK_EXPANSION),
            );
        }
    }
    output.push_str("</mask></defs>\n");
}

fn render_item(output: &mut String, item: &DepictionItem, atom_mask: bool) {
    output.push_str("<g data-umol-item=\"");
    output.push_str(item_kind(item));
    output.push('"');
    if let DepictionItem::Marker(marker) = item {
        output.push_str(" data-umol-marker=\"");
        output.push_str(marker_kind(marker.kind));
        output.push('"');
    }
    write_references(output, item.references());
    if atom_mask {
        write!(output, r##" mask="url(#{ATOM_LABEL_MASK_ID})""##)
            .expect("writing to a String cannot fail");
    }
    output.push('>');

    match item {
        DepictionItem::Atom(atom) => render_atom(output, atom),
        DepictionItem::Bond(bond) => render_bond(output, bond),
        DepictionItem::Wedge(wedge) => render_wedge(output, wedge),
        DepictionItem::DashedContour(contour) => render_dashed_contour(output, contour),
        DepictionItem::Text(text) => render_text(output, text),
        DepictionItem::Marker(marker) => render_marker(output, marker),
        DepictionItem::Arrow(arrow) => render_arrow(output, arrow),
    }

    output.push_str("</g>\n");
}

fn item_uses_atom_mask(item: &DepictionItem) -> bool {
    matches!(
        item,
        DepictionItem::Bond(_) | DepictionItem::Wedge(_) | DepictionItem::DashedContour(_)
    )
}

fn item_kind(item: &DepictionItem) -> &'static str {
    match item {
        DepictionItem::Atom(_) => "atom",
        DepictionItem::Bond(_) => "bond",
        DepictionItem::Wedge(_) => "wedge",
        DepictionItem::DashedContour(_) => "dashed-contour",
        DepictionItem::Text(_) => "text",
        DepictionItem::Marker(_) => "marker",
        DepictionItem::Arrow(_) => "arrow",
    }
}

fn render_atom(output: &mut String, atom: &AtomItem) {
    render_atom_glyph(output, "umol-atom", atom, "currentColor", None);
}

fn render_atom_glyph(
    output: &mut String,
    class: &str,
    atom: &AtomItem,
    fill: &str,
    stroke_width: Option<f64>,
) {
    output.push_str("<text class=\"");
    output.push_str(class);
    output.push('"');
    write_point_attributes(output, atom.position);
    write!(
        output,
        r#" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="{TEXT_SIZE}" fill="{fill}""#
    )
    .expect("writing to a String cannot fail");
    if let Some(stroke_width) = stroke_width {
        output.push_str(r#" stroke="black" stroke-width=""#);
        write_number(output, stroke_width);
        output.push_str(r#"" stroke-linecap="round" stroke-linejoin="round""#);
    }
    output.push('>');
    write_escaped_text(output, &atom.label);
    output.push_str("</text>");
}

fn render_bond(output: &mut String, bond: &BondItem) {
    let dx = bond.end.x - bond.start.x;
    let dy = bond.end.y - bond.start.y;
    let length = dx.hypot(dy);
    let perpendicular = if length == 0.0 {
        Point2D::new(0.0, 0.0)
    } else {
        Point2D::new(-dy / length, dx / length)
    };
    let center = (f64::from(bond.line_count) - 1.0) / 2.0;

    for index in 0..bond.line_count {
        let distance = (f64::from(index) - center) * BOND_GAP;
        let offset = Point2D::new(perpendicular.x * distance, perpendicular.y * distance);
        render_line(
            output,
            "umol-bond",
            Point2D::new(bond.start.x + offset.x, bond.start.y + offset.y),
            Point2D::new(bond.end.x + offset.x, bond.end.y + offset.y),
        );
    }
}

fn render_wedge(output: &mut String, wedge: &WedgeItem) {
    let dx = wedge.base.x - wedge.tip.x;
    let dy = wedge.base.y - wedge.tip.y;
    let length = dx.hypot(dy);
    let perpendicular = if length == 0.0 {
        Point2D::new(0.0, 0.0)
    } else {
        Point2D::new(-dy / length, dx / length)
    };

    match wedge.kind {
        WedgeKind::Solid => {
            let first = Point2D::new(
                wedge.base.x + WEDGE_BASE_HALF_WIDTH * perpendicular.x,
                wedge.base.y + WEDGE_BASE_HALF_WIDTH * perpendicular.y,
            );
            let second = Point2D::new(
                wedge.base.x - WEDGE_BASE_HALF_WIDTH * perpendicular.x,
                wedge.base.y - WEDGE_BASE_HALF_WIDTH * perpendicular.y,
            );
            output.push_str(r#"<polygon class="umol-wedge-solid" points=""#);
            write_svg_point(output, wedge.tip);
            output.push(' ');
            write_svg_point(output, first);
            output.push(' ');
            write_svg_point(output, second);
            output.push_str(r#"" fill="currentColor"/>"#);
        }
        WedgeKind::Hashed => {
            for index in 1..=WEDGE_HASH_COUNT {
                let numerator = index as f64;
                let denominator = WEDGE_HASH_COUNT as f64;
                let center = Point2D::new(
                    wedge.tip.x + dx * numerator / denominator,
                    wedge.tip.y + dy * numerator / denominator,
                );
                let half_width = (WEDGE_BASE_HALF_WIDTH / denominator) * numerator;
                let first = Point2D::new(
                    center.x + half_width * perpendicular.x,
                    center.y + half_width * perpendicular.y,
                );
                let second = Point2D::new(
                    center.x - half_width * perpendicular.x,
                    center.y - half_width * perpendicular.y,
                );
                render_line(output, "umol-wedge-hash", first, second);
            }
        }
    }
}

fn render_dashed_contour(output: &mut String, contour: &DashedContourItem) {
    output.push_str(r#"<path class="umol-dashed-contour" d=""#);
    if let Some((first, rest)) = contour.points.split_first() {
        output.push('M');
        write_svg_point(output, *first);
        for &point in rest {
            output.push_str(" L");
            write_svg_point(output, point);
        }
        if contour.closed {
            output.push_str(" Z");
        }
    }
    write!(
        output,
        r#"" fill="none" stroke="currentColor" stroke-width="{DASHED_CONTOUR_WIDTH}" stroke-dasharray="{DASHED_CONTOUR_PATTERN}" stroke-linecap="round" stroke-linejoin="round"/>"#
    )
    .expect("writing to a String cannot fail");
}

fn render_text(output: &mut String, text: &TextItem) {
    output.push_str(r#"<text class="umol-text""#);
    write_point_attributes(output, text.position);
    write!(
        output,
        r#" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="{TEXT_SIZE}" fill="currentColor">"#
    )
    .expect("writing to a String cannot fail");
    write_escaped_text(output, &text.text);
    output.push_str("</text>");
}

fn render_marker(output: &mut String, marker: &MarkerItem) {
    let (class, radius, dash) = match marker.kind {
        MarkerKind::Aromatic => (
            "umol-aromatic-marker",
            AROMATIC_MARKER_RADIUS,
            Some("0.05 0.05"),
        ),
        MarkerKind::Stereo => ("umol-stereo-marker", STEREO_MARKER_RADIUS, None),
    };
    output.push_str("<circle class=\"");
    output.push_str(class);
    output.push('"');
    write_point_attributes_with_names(output, marker.position, "cx", "cy");
    write!(
        output,
        r#" r="{radius}" fill="none" stroke="currentColor" stroke-width="{MARKER_WIDTH}""#
    )
    .expect("writing to a String cannot fail");
    if let Some(dash) = dash {
        write!(output, r#" stroke-dasharray="{dash}""#).expect("writing to a String cannot fail");
    }
    output.push_str("/>");
}

fn render_arrow(output: &mut String, arrow: &ArrowItem) {
    render_line(output, "umol-arrow-shaft", arrow.start, arrow.end);

    let dx = arrow.end.x - arrow.start.x;
    let dy = arrow.end.y - arrow.start.y;
    let length = dx.hypot(dy);
    if length == 0.0 {
        return;
    }
    let along = Point2D::new(dx / length, dy / length);
    let perpendicular = Point2D::new(-along.y, along.x);
    let base = Point2D::new(
        arrow.end.x - ARROW_HEAD_LENGTH * along.x,
        arrow.end.y - ARROW_HEAD_LENGTH * along.y,
    );
    let first = Point2D::new(
        base.x + ARROW_HEAD_HALF_WIDTH * perpendicular.x,
        base.y + ARROW_HEAD_HALF_WIDTH * perpendicular.y,
    );
    let second = Point2D::new(
        base.x - ARROW_HEAD_HALF_WIDTH * perpendicular.x,
        base.y - ARROW_HEAD_HALF_WIDTH * perpendicular.y,
    );

    output.push_str(r#"<polygon class="umol-arrow-head" points=""#);
    write_svg_point(output, arrow.end);
    output.push(' ');
    write_svg_point(output, first);
    output.push(' ');
    write_svg_point(output, second);
    output.push_str(r#"" fill="currentColor"/>"#);
}

fn render_line(output: &mut String, class: &str, start: Point2D, end: Point2D) {
    output.push_str("<line class=\"");
    output.push_str(class);
    output.push('"');
    write_coordinate_attribute(output, "x1", start.x);
    write_coordinate_attribute(output, "y1", -start.y);
    write_coordinate_attribute(output, "x2", end.x);
    write_coordinate_attribute(output, "y2", -end.y);
    write!(
        output,
        r#" fill="none" stroke="currentColor" stroke-width="{BOND_WIDTH}" stroke-linecap="round""#
    )
    .expect("writing to a String cannot fail");
    output.push_str("/>");
}

fn write_point_attributes(output: &mut String, point: Point2D) {
    write_point_attributes_with_names(output, point, "x", "y");
}

fn write_point_attributes_with_names(
    output: &mut String,
    point: Point2D,
    x_name: &str,
    y_name: &str,
) {
    write_coordinate_attribute(output, x_name, point.x);
    write_coordinate_attribute(output, y_name, -point.y);
}

fn write_coordinate_attribute(output: &mut String, name: &str, value: f64) {
    output.push(' ');
    output.push_str(name);
    output.push_str("=\"");
    write_number(output, value);
    output.push('"');
}

fn write_svg_point(output: &mut String, point: Point2D) {
    write_number(output, point.x);
    output.push(',');
    write_number(output, -point.y);
}

fn write_number(output: &mut String, value: f64) {
    if value == 0.0 {
        output.push('0');
    } else {
        write!(output, "{value}").expect("writing to a String cannot fail");
    }
}

fn write_escaped_text(output: &mut String, text: &str) {
    for character in text.chars() {
        match character {
            '&' => output.push_str("&amp;"),
            '<' => output.push_str("&lt;"),
            '>' => output.push_str("&gt;"),
            '"' => output.push_str("&quot;"),
            '\'' => output.push_str("&apos;"),
            character => output.push(character),
        }
    }
}

fn write_references(output: &mut String, references: &[DepictionReference]) {
    if references.is_empty() {
        return;
    }
    output.push_str(" data-umol-references=\"");
    for (index, reference) in references.iter().enumerate() {
        if index != 0 {
            output.push(' ');
        }
        write_reference(output, *reference);
    }
    output.push('"');
}

fn write_reference(output: &mut String, reference: DepictionReference) {
    match reference {
        DepictionReference::Molecule(entity) => write_entity_reference(output, "molecule", entity),
        DepictionReference::ReactionLhs(entity) => {
            write_entity_reference(output, "reaction-lhs", entity);
        }
        DepictionReference::ReactionRhs(entity) => {
            write_entity_reference(output, "reaction-rhs", entity);
        }
        DepictionReference::CorrespondencePair(index) => {
            write!(output, "correspondence-pair/{index}").expect("writing to a String cannot fail");
        }
        DepictionReference::Delta(ordinal) => {
            write!(output, "delta/{ordinal}").expect("writing to a String cannot fail");
        }
    }
}

fn write_entity_reference(output: &mut String, frame: &str, entity: Entity) {
    output.push_str(frame);
    output.push('/');
    match entity {
        Entity::Atom(id) => write!(output, "atom/{}", id.index()),
        Entity::Bond(id) => write!(output, "bond/{}", id.index()),
        Entity::DativeBond(id) => write!(output, "dative-bond/{}", id.index()),
        Entity::AromaticSystem(id) => write!(output, "aromatic-system/{}", id.index()),
        Entity::MulticenterBond(id) => write!(output, "multicenter-bond/{}", id.index()),
        Entity::NoncovalentBond(id) => write!(output, "noncovalent-bond/{}", id.index()),
        Entity::StereoAtom(id) => write!(output, "stereo-atom/{}", id.index()),
        Entity::StereoBond(id) => write!(output, "stereo-bond/{}", id.index()),
    }
    .expect("writing to a String cannot fail");
}

fn marker_kind(kind: MarkerKind) -> &'static str {
    match kind {
        MarkerKind::Aromatic => "aromatic",
        MarkerKind::Stereo => "stereo",
    }
}

#[cfg(test)]
mod tests {
    use roxmltree::Document;
    use rstest::rstest;
    use umol_graph_core::Correspondence;
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, Entity, Molecule};
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::depiction::molecule::depict;
    use crate::depiction::reaction::depict_from_sides;
    use crate::layout::MoleculeLayout;

    #[rstest]
    fn test_render_atom_mask() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#);
        let layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)]).unwrap();
        let depiction = depict(&molecule, &layout).unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let mask = root
            .descendants()
            .find(|node| node.has_tag_name("mask"))
            .unwrap();
        let mask_children = mask
            .children()
            .filter(|node| node.is_element())
            .collect::<Vec<_>>();
        let groups = root
            .children()
            .filter(|node| node.has_tag_name("g"))
            .collect::<Vec<_>>();

        assert_eq!(mask.attribute("id"), Some(ATOM_LABEL_MASK_ID));
        assert_eq!(mask.attribute("maskUnits"), Some("userSpaceOnUse"));
        assert_eq!(mask.attribute("maskContentUnits"), Some("userSpaceOnUse"));
        assert_eq!(mask.attribute("mask-type"), Some("luminance"));
        assert_eq!(
            mask_children
                .iter()
                .map(|node| node.tag_name().name())
                .collect::<Vec<_>>(),
            ["rect", "text"]
        );
        assert_eq!(mask_children[0].attribute("fill"), Some("white"));
        assert_eq!(mask_children[1].attribute("class"), Some("umol-atom-mask"));
        assert_eq!(mask_children[1].attribute("fill"), Some("black"));
        assert_eq!(mask_children[1].attribute("stroke"), Some("black"));
        assert_eq!(mask_children[1].attribute("stroke-width"), Some("0.6"));
        assert_eq!(mask_children[1].text(), Some("N"));
        assert_eq!(
            groups
                .iter()
                .map(|group| group.attribute("mask"))
                .collect::<Vec<_>>(),
            [Some("url(#umol-atom-label-mask)"), None]
        );
        assert_eq!(
            groups[0].attribute("data-umol-references"),
            Some("molecule/bond/0")
        );
        assert_eq!(
            groups[1].attribute("data-umol-references"),
            Some("molecule/atom/1")
        );
        assert_eq!(
            groups[0].first_element_child().unwrap().attribute("stroke"),
            Some("currentColor")
        );
        assert_eq!(
            groups[1].first_element_child().unwrap().attribute("fill"),
            Some("currentColor")
        );
        assert_eq!(
            root.children()
                .filter(|node| node.has_tag_name("rect"))
                .collect::<Vec<_>>(),
            []
        );
    }

    #[rstest]
    #[case::visible(
        "umol-atom",
        "currentColor",
        None,
        r#"<text class="umol-atom" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.45" fill="currentColor">&lt;&amp;&gt;&quot;&apos;</text>"#
    )]
    #[case::mask(
        "umol-atom-mask",
        "black",
        Some(0.6),
        r#"<text class="umol-atom-mask" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.45" fill="black" stroke="black" stroke-width="0.6" stroke-linecap="round" stroke-linejoin="round">&lt;&amp;&gt;&quot;&apos;</text>"#
    )]
    fn test_render_atom_glyph(
        #[case] class: &str,
        #[case] fill: &str,
        #[case] stroke_width: Option<f64>,
        #[case] expected: &str,
    ) {
        let atom = AtomItem {
            position: Point2D::new(1.0, -2.0),
            label: "<&>\"'".to_owned(),
            references: Vec::new(),
        };
        let mut output = String::new();

        render_atom_glyph(&mut output, class, &atom, fill, stroke_width);

        assert_eq!(output, expected);
    }

    #[rstest]
    #[case::solid(
        WedgeKind::Solid,
        r#"<polygon class="umol-wedge-solid" points="0,0 5,-0.1 5,0.1" fill="currentColor"/>"#
    )]
    #[case::hashed(
        WedgeKind::Hashed,
        concat!(
            r#"<line class="umol-wedge-hash" x1="1" y1="-0.02" x2="1" y2="0.02" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="2" y1="-0.04" x2="2" y2="0.04" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="3" y1="-0.06" x2="3" y2="0.06" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="4" y1="-0.08" x2="4" y2="0.08" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="5" y1="-0.1" x2="5" y2="0.1" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
        )
    )]
    fn test_render_wedge(#[case] kind: WedgeKind, #[case] expected_glyph: &str) {
        let item = DepictionItem::Wedge(WedgeItem {
            tip: Point2D::new(0.0, 0.0),
            base: Point2D::new(5.0, 0.0),
            kind,
            references: vec![DepictionReference::Molecule(Entity::Bond(BondId(3)))],
        });
        let mut output = String::new();

        render_item(&mut output, &item, true);

        assert_eq!(
            output,
            format!(
                r#"<g data-umol-item="wedge" data-umol-references="molecule/bond/3" mask="url(#umol-atom-label-mask)">{expected_glyph}</g>
"#
            )
        );
    }

    #[rstest]
    #[case::open(
        false,
        r#"<g data-umol-item="dashed-contour" data-umol-references="molecule/aromatic-system/2" mask="url(#umol-atom-label-mask)"><path class="umol-dashed-contour" d="M-1,-2 L0,1 L3,0" fill="none" stroke="currentColor" stroke-width="0.04" stroke-dasharray="0.12 0.1" stroke-linecap="round" stroke-linejoin="round"/></g>
"#
    )]
    #[case::closed(
        true,
        r#"<g data-umol-item="dashed-contour" data-umol-references="molecule/aromatic-system/2" mask="url(#umol-atom-label-mask)"><path class="umol-dashed-contour" d="M-1,-2 L0,1 L3,0 Z" fill="none" stroke="currentColor" stroke-width="0.04" stroke-dasharray="0.12 0.1" stroke-linecap="round" stroke-linejoin="round"/></g>
"#
    )]
    fn test_render_dashed_contour(#[case] closed: bool, #[case] expected: &str) {
        let item = DepictionItem::DashedContour(DashedContourItem {
            points: vec![
                Point2D::new(-1.0, 2.0),
                Point2D::new(0.0, -1.0),
                Point2D::new(3.0, 0.0),
            ],
            closed,
            references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(2),
            ))],
        });
        let mut output = String::new();

        render_item(&mut output, &item, true);

        assert_eq!(output, expected);
    }

    #[rstest]
    fn test_render_empty_depiction() {
        let molecule = Molecule::new();
        let layout = MoleculeLayout::try_new(Vec::new()).unwrap();
        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            render(&depiction),
            r#"<svg xmlns="http://www.w3.org/2000/svg" class="umol-depiction" viewBox="-0.5 -0.5 1 1">
</svg>"#
        );
    }

    #[rstest]
    fn test_render_molecule() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let layout =
            MoleculeLayout::try_new(vec![Point2D::new(-1.0, 2.0), Point2D::new(3.0, -4.0)])
                .unwrap();
        let depiction = depict(&molecule, &layout).unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let children = root
            .children()
            .filter(|child| child.is_element())
            .collect::<Vec<_>>();
        let groups = children
            .iter()
            .copied()
            .filter(|child| child.has_tag_name("g"))
            .collect::<Vec<_>>();

        assert_eq!(root.attribute("viewBox"), Some("-1.5 -2.5 5 7"));
        assert_eq!(
            children
                .iter()
                .map(|child| child.tag_name().name())
                .collect::<Vec<_>>(),
            ["defs", "g", "g"]
        );
        assert_eq!(
            groups
                .iter()
                .map(|group| group.attribute("data-umol-item").unwrap())
                .collect::<Vec<_>>(),
            ["bond", "atom"]
        );
        assert_eq!(
            groups[0].attribute("data-umol-references"),
            Some("molecule/bond/0")
        );
        assert_eq!(
            groups[1].attribute("data-umol-references"),
            Some("molecule/atom/1")
        );
        assert_eq!(
            groups[0].attribute("mask"),
            Some("url(#umol-atom-label-mask)")
        );

        let line = groups[0].first_element_child().unwrap();
        assert_eq!(line.tag_name().name(), "line");
        assert_eq!(line.attribute("x1"), Some("-1"));
        assert_eq!(line.attribute("y1"), Some("-2"));
        assert_eq!(line.attribute("x2"), Some("3"));
        assert_eq!(line.attribute("y2"), Some("4"));

        let oxygen = groups[1].first_element_child().unwrap();
        assert_eq!(oxygen.attribute("x"), Some("3"));
        assert_eq!(oxygen.attribute("y"), Some("4"));
        assert_eq!(oxygen.text(), Some("O"));
    }

    #[rstest]
    fn test_render_reaction() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let rhs = mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#);
        let lhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(-1.0, 0.0), Point2D::new(0.0, 0.0)]).unwrap();
        let rhs_layout =
            MoleculeLayout::try_new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();
        let correspondence =
            Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 2)
                .unwrap();
        let depiction =
            depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence).unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let mask_labels = root
            .descendants()
            .filter(|node| node.attribute("class") == Some("umol-atom-mask"))
            .map(|node| node.text().unwrap())
            .collect::<Vec<_>>();
        let bond_groups = root
            .children()
            .filter(|node| node.attribute("data-umol-item") == Some("bond"))
            .collect::<Vec<_>>();
        let arrow = root
            .children()
            .find(|node| node.attribute("data-umol-item") == Some("arrow"))
            .unwrap();

        assert_eq!(mask_labels, ["O", "N"]);
        assert_eq!(
            bond_groups
                .iter()
                .map(|group| group.attribute("data-umol-references"))
                .collect::<Vec<_>>(),
            [Some("reaction-lhs/bond/0"), Some("reaction-rhs/bond/0")]
        );
        assert_eq!(
            bond_groups
                .iter()
                .map(|group| group.attribute("mask"))
                .collect::<Vec<_>>(),
            [
                Some("url(#umol-atom-label-mask)"),
                Some("url(#umol-atom-label-mask)"),
            ]
        );
        assert_eq!(arrow.attribute("mask"), None);
    }
}
