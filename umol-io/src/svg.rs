//! Deterministic SVG rendering of format-neutral depictions.

use std::fmt::Write;

use umol_geometric_core::Point2D;
use umol_graph_ir::ir::Entity;

use crate::depiction::{
    ArrowItem, AtomItem, BondItem, Bounds, Depiction, DepictionItem, DepictionReference,
    MarkerItem, MarkerKind, TextItem,
};

const VIEW_MARGIN: f64 = 0.5;
const BOND_GAP: f64 = 0.12;
const BOND_WIDTH: f64 = 0.06;
const TEXT_SIZE: f64 = 0.45;
const AROMATIC_MARKER_RADIUS: f64 = 0.16;
const STEREO_MARKER_RADIUS: f64 = 0.19;
const MARKER_WIDTH: f64 = 0.04;
const ARROW_HEAD_LENGTH: f64 = 0.24;
const ARROW_HEAD_HALF_WIDTH: f64 = 0.11;

/// Renders a [`Depiction`] as a complete SVG document fragment.
///
/// Depiction item order is preserved as direct child order. Coordinates are converted from the
/// depiction's y-up convention to SVG's y-down convention. The view box extends the depiction's
/// anchor bounds by half a nominal bond length on each side; an empty depiction uses a centered
/// one-by-one view box. Structured references are encoded in the `data-umol-references` attribute
/// as ordered, space-separated path tokens: molecular entity references use `molecule`,
/// `reaction-lhs`, or `reaction-rhs` followed by entity kind and id; `correspondence-pair` and
/// `delta` references use their kind followed by their zero-based position.
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

    for item in depiction.items() {
        render_item(&mut output, item);
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

fn render_item(output: &mut String, item: &DepictionItem) {
    output.push_str("<g data-umol-item=\"");
    output.push_str(item_kind(item));
    output.push('"');
    if let DepictionItem::Marker(marker) = item {
        output.push_str(" data-umol-marker=\"");
        output.push_str(marker_kind(marker.kind));
        output.push('"');
    }
    write_references(output, item.references());
    output.push('>');

    match item {
        DepictionItem::Atom(atom) => render_atom(output, atom),
        DepictionItem::Bond(bond) => render_bond(output, bond),
        DepictionItem::Text(text) => render_text(output, text),
        DepictionItem::Marker(marker) => render_marker(output, marker),
        DepictionItem::Arrow(arrow) => render_arrow(output, arrow),
    }

    output.push_str("</g>\n");
}

fn item_kind(item: &DepictionItem) -> &'static str {
    match item {
        DepictionItem::Atom(_) => "atom",
        DepictionItem::Bond(_) => "bond",
        DepictionItem::Text(_) => "text",
        DepictionItem::Marker(_) => "marker",
        DepictionItem::Arrow(_) => "arrow",
    }
}

fn render_atom(output: &mut String, atom: &AtomItem) {
    output.push_str(r#"<text class="umol-atom""#);
    write_point_attributes(output, atom.position);
    write!(
        output,
        r#" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="{TEXT_SIZE}" fill="currentColor">"#
    )
    .expect("writing to a String cannot fail");
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
    use umol_graph_ir::ir::Molecule;
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::depiction::molecule::depict;
    use crate::layout::MoleculeLayout;

    #[test]
    fn test_render_empty_depiction() {
        let molecule = Molecule::new();
        let layout = MoleculeLayout::try_new(Vec::new()).unwrap();
        let depiction = depict(&molecule, &layout).unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();

        assert_eq!(root.tag_name().name(), "svg");
        assert_eq!(
            root.tag_name().namespace(),
            Some("http://www.w3.org/2000/svg")
        );
        assert_eq!(root.attribute("class"), Some("umol-depiction"));
        assert_eq!(root.attribute("viewBox"), Some("-0.5 -0.5 1 1"));
        assert_eq!(
            root.children().filter(|child| child.is_element()).count(),
            0
        );
    }

    #[test]
    fn test_render_preserves_molecule_item_order_and_coordinates() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let layout =
            MoleculeLayout::try_new(vec![Point2D::new(-1.0, 2.0), Point2D::new(3.0, -4.0)])
                .unwrap();
        let depiction = depict(&molecule, &layout).unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let groups: Vec<_> = root.children().filter(|child| child.is_element()).collect();

        assert_eq!(root.attribute("viewBox"), Some("-1.5 -2.5 5 7"));
        assert_eq!(
            groups
                .iter()
                .map(|group| group.attribute("data-umol-item").unwrap())
                .collect::<Vec<_>>(),
            ["bond", "atom", "atom"]
        );
        assert_eq!(
            groups[0].attribute("data-umol-references"),
            Some("molecule/bond/0")
        );
        assert_eq!(
            groups[1].attribute("data-umol-references"),
            Some("molecule/atom/0")
        );
        assert_eq!(
            groups[2].attribute("data-umol-references"),
            Some("molecule/atom/1")
        );

        let line = groups[0].first_element_child().unwrap();
        assert_eq!(line.tag_name().name(), "line");
        assert_eq!(line.attribute("x1"), Some("-1"));
        assert_eq!(line.attribute("y1"), Some("-2"));
        assert_eq!(line.attribute("x2"), Some("3"));
        assert_eq!(line.attribute("y2"), Some("4"));

        let carbon = groups[1].first_element_child().unwrap();
        let oxygen = groups[2].first_element_child().unwrap();
        assert_eq!(carbon.attribute("x"), Some("-1"));
        assert_eq!(carbon.attribute("y"), Some("-2"));
        assert_eq!(carbon.text(), Some("C"));
        assert_eq!(oxygen.attribute("x"), Some("3"));
        assert_eq!(oxygen.attribute("y"), Some("4"));
        assert_eq!(oxygen.text(), Some("O"));
    }
}
