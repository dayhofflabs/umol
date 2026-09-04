//! Deterministic SVG rendering of format-neutral depictions.

use std::fmt::Write;

use umol_geometric_core::Point2D;
use umol_graph_ir::ir::Entity;

use crate::depict::{
    ArrowItem, AtomItem, AtomLabel, BondItem, Bounds, DashedContourItem, Depiction, DepictionItem,
    DepictionReference, TextItem, WedgeItem, WedgeKind,
};

const VIEW_MARGIN: f64 = 0.5;
const BOND_GAP: f64 = 0.12;
const BOND_WIDTH: f64 = 0.06;
const WEDGE_TIP_HALF_WIDTH: f64 = 0.03;
const WEDGE_BASE_HALF_WIDTH: f64 = 0.1;
const WEDGE_HASH_COUNT: usize = 8;
const DASHED_CONTOUR_WIDTH: f64 = 0.04;
const DASHED_CONTOUR_PATTERN: &str = "0.12 0.1";
const TEXT_SIZE: f64 = 0.45;
const MAPPING_INDEX_TEXT_SIZE: f64 = TEXT_SIZE * 0.85;
const SCRIPT_TEXT_SIZE: f64 = 0.315;
const ATOM_LABEL_MASK_ID: &str = "umol-atom-label-mask";
const ATOM_LABEL_CHARACTER_ADVANCE: f64 = 0.36;
const ATOM_LABEL_SCRIPT_CHARACTER_ADVANCE: f64 = 0.252;
const ATOM_LABEL_HORIZONTAL_CLEARANCE: f64 = 0.08;
const ATOM_LABEL_VERTICAL_CLEARANCE: f64 = 0.08;
const ATOM_LABEL_BASE_HALF_HEIGHT: f64 = 0.2475;
const ATOM_LABEL_SCRIPT_HALF_HEIGHT: f64 = 0.17325;
const ATOM_LABEL_SUPERSCRIPT_RISE: f64 = 0.1575;
const ATOM_LABEL_SUBSCRIPT_DROP: f64 = 0.1125;
const ARROW_HEAD_LENGTH: f64 = 0.24;
const ARROW_HEAD_HALF_WIDTH: f64 = 0.11;

/// Renders a [`Depiction`] as a complete SVG document fragment.
///
/// Depiction item order is preserved among item groups. When atom labels are present, one leading
/// definition masks molecular strokes beneath continuous conservative label rectangles without
/// painting a page-background color. Coordinates are converted from the depiction's y-up
/// convention to SVG's y-down convention. The view box covers both depiction anchors and estimated
/// atom-label extents, then adds half a nominal bond length on each side; an empty depiction uses a
/// centered one-by-one view box. Structured
/// references are encoded in the `data-umol-references` attribute as ordered, space-separated path
/// tokens: molecular entity references use `molecule`, `reaction-lhs`, or `reaction-rhs` followed
/// by entity kind and id; `correspondence-pair` and `delta` references use their kind followed by
/// their zero-based position.
pub(crate) fn render(depiction: &Depiction) -> String {
    let mut output = String::new();
    let view_box = SvgViewBox::from_depiction(depiction);

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
    fn from_depiction(depiction: &Depiction) -> Self {
        let mut extents = depiction.bounds().map(SvgExtents::from_depiction_bounds);
        for item in depiction.items() {
            if let DepictionItem::Atom(atom) = item {
                let label_box = atom_label_box(atom);
                include_svg_box(&mut extents, label_box);
            }
        }

        match extents {
            Some(extents) => Self {
                x: extents.min_x - VIEW_MARGIN,
                y: extents.min_y - VIEW_MARGIN,
                width: extents.max_x - extents.min_x + 2.0 * VIEW_MARGIN,
                height: extents.max_y - extents.min_y + 2.0 * VIEW_MARGIN,
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

#[derive(Clone, Copy)]
struct SvgExtents {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
}

impl SvgExtents {
    fn from_depiction_bounds(bounds: &Bounds) -> Self {
        Self {
            min_x: bounds.min.x,
            min_y: -bounds.max.y,
            max_x: bounds.max.x,
            max_y: -bounds.min.y,
        }
    }

    fn include_box(&mut self, label_box: SvgBox) {
        self.min_x = self.min_x.min(label_box.x);
        self.min_y = self.min_y.min(label_box.y);
        self.max_x = self.max_x.max(label_box.x + label_box.width);
        self.max_y = self.max_y.max(label_box.y + label_box.height);
    }
}

#[derive(Clone, Copy)]
struct SvgBox {
    x: f64,
    y: f64,
    width: f64,
    height: f64,
}

fn include_svg_box(extents: &mut Option<SvgExtents>, label_box: SvgBox) {
    match extents {
        Some(extents) => extents.include_box(label_box),
        None => {
            *extents = Some(SvgExtents {
                min_x: label_box.x,
                min_y: label_box.y,
                max_x: label_box.x + label_box.width,
                max_y: label_box.y + label_box.height,
            });
        }
    }
}

fn atom_label_box(atom: &AtomItem) -> SvgBox {
    let label = &atom.label;
    let base_width = character_count(&label.base) * ATOM_LABEL_CHARACTER_ADVANCE;
    let script_width = [
        &label.left_superscript,
        &label.right_subscript,
        &label.right_superscript,
    ]
    .into_iter()
    .flatten()
    .map(|text| character_count(text) * ATOM_LABEL_SCRIPT_CHARACTER_ADVANCE)
    .sum::<f64>();
    let width = base_width + script_width + 2.0 * ATOM_LABEL_HORIZONTAL_CLEARANCE;
    let has_superscript = label.left_superscript.is_some() || label.right_superscript.is_some();
    let top = if has_superscript {
        (ATOM_LABEL_SUPERSCRIPT_RISE + ATOM_LABEL_SCRIPT_HALF_HEIGHT)
            .max(ATOM_LABEL_BASE_HALF_HEIGHT)
    } else {
        ATOM_LABEL_BASE_HALF_HEIGHT
    } + ATOM_LABEL_VERTICAL_CLEARANCE;
    let bottom = if label.right_subscript.is_some() {
        (ATOM_LABEL_SUBSCRIPT_DROP + ATOM_LABEL_SCRIPT_HALF_HEIGHT).max(ATOM_LABEL_BASE_HALF_HEIGHT)
    } else {
        ATOM_LABEL_BASE_HALF_HEIGHT
    } + ATOM_LABEL_VERTICAL_CLEARANCE;

    SvgBox {
        x: atom.position.x - width / 2.0,
        y: -atom.position.y - top,
        width,
        height: top + bottom,
    }
}

fn character_count(text: &str) -> f64 {
    text.chars().count() as f64
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
            let label_box = atom_label_box(atom);
            output.push_str(r#"<rect class="umol-atom-mask" x=""#);
            write_number(output, label_box.x);
            output.push_str(r#"" y=""#);
            write_number(output, label_box.y);
            output.push_str(r#"" width=""#);
            write_number(output, label_box.width);
            output.push_str(r#"" height=""#);
            write_number(output, label_box.height);
            output.push_str(r#"" fill="black"/>"#);
        }
    }
    output.push_str("</mask></defs>\n");
}

fn render_item(output: &mut String, item: &DepictionItem, atom_mask: bool) {
    output.push_str("<g data-umol-item=\"");
    output.push_str(item_kind(item));
    output.push('"');
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

    let AtomLabel {
        base,
        left_superscript,
        right_subscript,
        right_superscript,
    } = &atom.label;
    if let Some(text) = left_superscript {
        render_atom_script(output, "umol-atom-left-superscript", "super", text);
    }
    output.push_str(r#"<tspan class="umol-atom-base">"#);
    write_escaped_text(output, base);
    output.push_str("</tspan>");
    if let Some(text) = right_subscript {
        render_atom_script(output, "umol-atom-right-subscript", "sub", text);
    }
    if let Some(text) = right_superscript {
        render_atom_script(output, "umol-atom-right-superscript", "super", text);
    }
    output.push_str("</text>");
}

fn render_atom_script(output: &mut String, class: &str, baseline_shift: &str, text: &str) {
    write!(
        output,
        r#"<tspan class="{class}" font-size="{SCRIPT_TEXT_SIZE}" baseline-shift="{baseline_shift}">"#
    )
    .expect("writing to a String cannot fail");
    write_escaped_text(output, text);
    output.push_str("</tspan>");
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
            let tip_first = Point2D::new(
                wedge.tip.x + WEDGE_TIP_HALF_WIDTH * perpendicular.x,
                wedge.tip.y + WEDGE_TIP_HALF_WIDTH * perpendicular.y,
            );
            let base_first = Point2D::new(
                wedge.base.x + WEDGE_BASE_HALF_WIDTH * perpendicular.x,
                wedge.base.y + WEDGE_BASE_HALF_WIDTH * perpendicular.y,
            );
            let base_second = Point2D::new(
                wedge.base.x - WEDGE_BASE_HALF_WIDTH * perpendicular.x,
                wedge.base.y - WEDGE_BASE_HALF_WIDTH * perpendicular.y,
            );
            let tip_second = Point2D::new(
                wedge.tip.x - WEDGE_TIP_HALF_WIDTH * perpendicular.x,
                wedge.tip.y - WEDGE_TIP_HALF_WIDTH * perpendicular.y,
            );
            output.push_str(r#"<polygon class="umol-wedge-solid" points=""#);
            write_svg_point(output, tip_first);
            output.push(' ');
            write_svg_point(output, base_first);
            output.push(' ');
            write_svg_point(output, base_second);
            output.push(' ');
            write_svg_point(output, tip_second);
            output.push_str(r#"" fill="currentColor"/>"#);
        }
        WedgeKind::Hashed => {
            for index in 0..WEDGE_HASH_COUNT {
                let numerator = index as f64;
                let denominator = (WEDGE_HASH_COUNT - 1) as f64;
                let center = Point2D::new(
                    wedge.tip.x + dx * numerator / denominator,
                    wedge.tip.y + dy * numerator / denominator,
                );
                let half_width = WEDGE_TIP_HALF_WIDTH
                    + (WEDGE_BASE_HALF_WIDTH - WEDGE_TIP_HALF_WIDTH) * numerator / denominator;
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
    let text_size = if text
        .references
        .iter()
        .any(|reference| matches!(reference, DepictionReference::CorrespondencePair(_)))
    {
        MAPPING_INDEX_TEXT_SIZE
    } else {
        TEXT_SIZE
    };
    output.push_str(r#"<text class="umol-text""#);
    write_point_attributes(output, text.position);
    write!(
        output,
        r#" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="{text_size}" fill="currentColor">"#
    )
    .expect("writing to a String cannot fail");
    write_escaped_text(output, &text.text);
    output.push_str("</text>");
}

fn render_arrow(output: &mut String, arrow: &ArrowItem) {
    let dx = arrow.end.x - arrow.start.x;
    let dy = arrow.end.y - arrow.start.y;
    let length = dx.hypot(dy);
    if length == 0.0 {
        render_line(output, "umol-arrow-shaft", arrow.start, arrow.end);
        return;
    }
    let along = Point2D::new(dx / length, dy / length);
    let perpendicular = Point2D::new(-along.y, along.x);
    let base = Point2D::new(
        arrow.end.x - ARROW_HEAD_LENGTH * along.x,
        arrow.end.y - ARROW_HEAD_LENGTH * along.y,
    );
    render_line(output, "umol-arrow-shaft", arrow.start, base);
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

#[cfg(test)]
mod tests {
    use roxmltree::Document;
    use rstest::rstest;
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, Entity, Molecule};
    #[cfg(feature = "coordgen")]
    use umol_graph_ir::ir::{Deltas, Reaction};
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::depict::molecule::depict;
    #[cfg(feature = "coordgen")]
    use crate::depict::Depict;
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
            ["rect", "rect"]
        );
        assert_eq!(mask_children[0].attribute("fill"), Some("white"));
        assert_eq!(mask_children[1].attribute("class"), Some("umol-atom-mask"));
        assert_eq!(mask_children[1].attribute("x"), Some("1.74"));
        assert_eq!(mask_children[1].attribute("y"), Some("-0.3275"));
        assert_eq!(mask_children[1].attribute("width"), Some("0.52"));
        assert_eq!(mask_children[1].attribute("height"), Some("0.655"));
        assert_eq!(mask_children[1].attribute("fill"), Some("black"));
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
    #[case::base(
        AtomLabel {
            base: "<&>\"'".to_owned(),
            left_superscript: None,
            right_subscript: None,
            right_superscript: None,
        },
        r#"<text class="umol-atom" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.45" fill="currentColor"><tspan class="umol-atom-base">&lt;&amp;&gt;&quot;&apos;</tspan></text>"#
    )]
    #[case::structured(
        AtomLabel {
            base: "CH".to_owned(),
            left_superscript: Some("13".to_owned()),
            right_subscript: Some("2".to_owned()),
            right_superscript: Some("+••".to_owned()),
        },
        concat!(
            r#"<text class="umol-atom" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.45" fill="currentColor">"#,
            r#"<tspan class="umol-atom-left-superscript" font-size="0.315" baseline-shift="super">13</tspan>"#,
            r#"<tspan class="umol-atom-base">CH</tspan>"#,
            r#"<tspan class="umol-atom-right-subscript" font-size="0.315" baseline-shift="sub">2</tspan>"#,
            r#"<tspan class="umol-atom-right-superscript" font-size="0.315" baseline-shift="super">+••</tspan>"#,
            "</text>",
        )
    )]
    fn test_render_atom(#[case] label: AtomLabel, #[case] expected: &str) {
        let atom = AtomItem {
            position: Point2D::new(1.0, -2.0),
            label,
            references: Vec::new(),
        };
        let mut output = String::new();

        render_atom(&mut output, &atom);

        assert_eq!(output, expected);
    }

    #[rstest]
    #[case::solid(
        WedgeKind::Solid,
        r#"<polygon class="umol-wedge-solid" points="0,-0.03 7,-0.1 7,0.1 0,0.03" fill="currentColor"/>"#
    )]
    #[case::hashed(
        WedgeKind::Hashed,
        concat!(
            r#"<line class="umol-wedge-hash" x1="0" y1="-0.03" x2="0" y2="0.03" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="1" y1="-0.04" x2="1" y2="0.04" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="2" y1="-0.05" x2="2" y2="0.05" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="3" y1="-0.06" x2="3" y2="0.06" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="4" y1="-0.07" x2="4" y2="0.07" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="5" y1="-0.08" x2="5" y2="0.08" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="6" y1="-0.09" x2="6" y2="0.09" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<line class="umol-wedge-hash" x1="7" y1="-0.1" x2="7" y2="0.1" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
        )
    )]
    fn test_render_wedge(#[case] kind: WedgeKind, #[case] expected_glyph: &str) {
        let item = DepictionItem::Wedge(WedgeItem {
            tip: Point2D::new(0.0, 0.0),
            base: Point2D::new(7.0, 0.0),
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
    #[case::annotation(
        Vec::new(),
        r#"<text class="umol-text" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.45" fill="currentColor">&lt;&amp;</text>"#
    )]
    #[case::mapping_index(
        vec![
            DepictionReference::ReactionLhs(Entity::Atom(AtomId(0))),
            DepictionReference::CorrespondencePair(3),
        ],
        r#"<text class="umol-text" x="1" y="2" text-anchor="middle" dominant-baseline="central" font-family="sans-serif" font-size="0.3825" fill="currentColor">&lt;&amp;</text>"#
    )]
    fn test_render_text(#[case] references: Vec<DepictionReference>, #[case] expected: &str) {
        let text = TextItem {
            position: Point2D::new(1.0, -2.0),
            text: "<&".to_owned(),
            references,
        };
        let mut output = String::new();

        render_text(&mut output, &text);

        assert_eq!(output, expected);
    }

    #[rstest]
    #[case::horizontal(
        ArrowItem {
            start: Point2D::new(-0.75, 0.0),
            end: Point2D::new(0.75, 0.0),
            references: Vec::new(),
        },
        concat!(
            r#"<line class="umol-arrow-shaft" x1="-0.75" y1="0" x2="0.51" y2="0" fill="none" stroke="currentColor" stroke-width="0.06" stroke-linecap="round"/>"#,
            r#"<polygon class="umol-arrow-head" points="0.75,0 0.51,-0.11 0.51,0.11" fill="currentColor"/>"#,
        )
    )]
    fn test_render_arrow(#[case] arrow: ArrowItem, #[case] expected: &str) {
        let mut output = String::new();

        render_arrow(&mut output, &arrow);

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

        assert_eq!(root.attribute("viewBox"), Some("-1.5 -2.5 5.26 7.3275"));
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
        let oxygen_base = oxygen.first_element_child().unwrap();
        assert_eq!(oxygen_base.attribute("class"), Some("umol-atom-base"));
        assert_eq!(oxygen_base.text(), Some("O"));
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_render_reaction() {
        let lhs = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let depiction = Reaction::new(lhs, Deltas::new()).depict().unwrap();

        let svg = render(&depiction);
        let document = Document::parse(&svg).unwrap();
        let root = document.root_element();
        let mask_boxes = root
            .descendants()
            .filter(|node| node.attribute("class") == Some("umol-atom-mask"))
            .collect::<Vec<_>>();
        let bond_groups = root
            .children()
            .filter(|node| node.attribute("data-umol-item") == Some("bond"))
            .collect::<Vec<_>>();
        let arrow = root
            .children()
            .find(|node| node.attribute("data-umol-item") == Some("arrow"))
            .unwrap();

        assert_eq!(mask_boxes.len(), 2);
        assert_eq!(mask_boxes[0].attribute("width"), Some("0.52"));
        assert_eq!(mask_boxes[1].attribute("width"), Some("0.52"));
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
