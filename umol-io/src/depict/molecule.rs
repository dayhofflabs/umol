//! Depiction construction from graph-IR molecules and supplied layouts.

use std::any::Any;
use std::fmt::Write;

use thiserror::Error;
use umol_chem::element::Element;
use umol_geometric_core::{complementary_direction, signed_volume, Point2D, Point3D};
use umol_graph_ir::ir::{
    AromaticSystemView, AsLit, AtomId, AtomView, BondId, Entity, IsotopeMass, Molecule,
    StereoAtomId, StereoAtomView, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
};
use umol_utils::error::UmolError;

use super::{
    AtomItem, AtomLabel, BondItem, DashedContourItem, Depiction, DepictionItem, DepictionReference,
    TextItem, WedgeItem, WedgeKind,
};
#[cfg(feature = "coordgen")]
use super::{Depict, DepictConfig};
use crate::layout::MoleculeLayout;
#[cfg(feature = "coordgen")]
use crate::layout::{layout_molecule, LayoutError};

const AROMATIC_CONTOUR_OFFSET: f64 = 0.18;
const AROMATIC_ANNOTATION_CLEARANCE: f64 = 0.28;
const AROMATIC_ANNOTATION_STEP: f64 = 0.35;
const AROMATIC_ANNOTATION_EXTERIOR_OFFSET: f64 = 0.35;
const GEOMETRY_EPSILON: f64 = 1.0e-9;
const MAX_CONTOUR_MITER: f64 = 4.0;

/// Constructs the first format-neutral depiction projection of `molecule` in `layout`.
///
/// Localized bonds and selected tetrahedral wedges are followed by visible atom labels, then one
/// trustworthy outer contour and any literal system annotation for each explicit aromatic system;
/// each group is ordered by graph-IR id. Carbon labels are omitted at non-isolated skeleton
/// vertices unless an isotope, charge, or radical count decorates the atom. Definite cis/trans
/// stereo is carried by the supplied coordinates. Nonliteral projected fields and unsupported
/// overlays, stereo kinds, local aromatic assertions, or constraints are omitted. The first
/// projection does not represent dative, multicenter, or noncovalent bonds, unprojected inherent
/// fields, or unsupported constraints. Crossed, degenerate, self-intersecting, and cage-like
/// aromatic projections receive no contour or system annotation.
///
/// # Errors
///
/// Returns [`MoleculeDepictionError::TetrahedralGeometry`] if a definite tetrahedral stereo atom
/// cannot be represented by a distinct, geometrically valid display wedge.
///
/// # Semantic properties
///
/// Every emitted tetrahedral wedge replaces its selected localized single bond. Reading the wedge
/// with the TableIR winding convention in the selected ligand frame reproduces the stored coset.
pub(crate) fn depict(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<Depiction, MoleculeDepictionError> {
    let wedges = tetrahedral_wedges(molecule, layout)?;
    let mut items = Vec::new();

    for bond in molecule.bonds().iter() {
        let Some(line_count) = bond
            .order()
            .as_lit()
            .and_then(|order| order.try_into().ok())
        else {
            continue;
        };
        let [first, second] = bond.atom_ids();
        let bond_reference = DepictionReference::Molecule(Entity::Bond(bond.id));
        if let Some(wedge) = wedges[bond.id.index()] {
            items.push(DepictionItem::Wedge(WedgeItem {
                tip: position(layout, wedge.tip),
                base: position(layout, wedge.base),
                kind: wedge.kind,
                references: vec![
                    bond_reference,
                    DepictionReference::Molecule(Entity::StereoAtom(wedge.stereo_atom)),
                ],
            }));
        } else {
            items.push(DepictionItem::Bond(BondItem {
                start: position(layout, first),
                end: position(layout, second),
                line_count,
                references: vec![bond_reference],
            }));
        }
    }

    for atom in molecule.atoms().iter() {
        let Some(label) = atom_label(atom) else {
            continue;
        };
        items.push(DepictionItem::Atom(AtomItem {
            position: position(layout, atom.id),
            label,
            references: vec![DepictionReference::Molecule(Entity::Atom(atom.id))],
        }));
    }

    for system in molecule.aromatic_systems().iter() {
        let Some(contour) = aromatic_contour(molecule, layout, system) else {
            continue;
        };
        let annotation = aromatic_annotation(molecule, layout, system, &contour.points);
        items.push(DepictionItem::DashedContour(contour));
        if let Some(annotation) = annotation {
            items.push(DepictionItem::Text(annotation));
        }
    }

    Ok(Depiction::from_items(items))
}

#[cfg(feature = "coordgen")]
impl Depict for Molecule {
    type Error = MoleculeDepictionError;

    fn depict_with(&self, config: &DepictConfig) -> Result<Depiction, Self::Error> {
        let layout = layout_molecule(self, config.layout_algorithm)?;
        depict(self, &layout)
    }
}

/// Failures while depicting a graph-IR [`Molecule`].
#[derive(Clone, Debug, Error, PartialEq)]
pub enum MoleculeDepictionError {
    /// The configured layout backend could not produce coordinates.
    #[cfg(feature = "coordgen")]
    #[error("layout: {0}")]
    Layout(#[from] LayoutError),
    /// A definite tetrahedral stereo atom could not be represented by a display wedge.
    #[error("tetrahedral geometry cannot establish a display wedge for stereo atom {stereo_atom}")]
    TetrahedralGeometry { stereo_atom: StereoAtomId },
}

impl UmolError for MoleculeDepictionError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[derive(Clone, Copy)]
struct SelectedWedge {
    stereo_atom: StereoAtomId,
    tip: AtomId,
    base: AtomId,
    kind: WedgeKind,
}

#[derive(Clone, Copy)]
struct WedgeCandidate {
    bond: BondId,
    base: AtomId,
    kind: WedgeKind,
}

struct TetrahedralCandidates {
    stereo_atom: StereoAtomId,
    site: AtomId,
    wedges: Vec<WedgeCandidate>,
}

fn tetrahedral_wedges(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<Vec<Option<SelectedWedge>>, MoleculeDepictionError> {
    let stereos = molecule
        .stereo_atoms()
        .iter()
        .filter(|stereo| {
            stereo
                .attributes
                .configuration
                .as_lit()
                .is_some_and(|configuration| configuration.kind == StereoKind::Tetrahedral)
        })
        .collect::<Vec<_>>();
    let mut tetrahedral_sites = vec![false; molecule.atoms().count()];
    for stereo in &stereos {
        tetrahedral_sites[stereo.site_id().index()] = true;
    }

    let candidates = stereos
        .into_iter()
        .map(|stereo| tetrahedral_candidates(molecule, layout, stereo, &tetrahedral_sites))
        .collect::<Vec<_>>();
    let mut bond_owners = vec![None; molecule.bonds().count()];
    for stereo_index in 0..candidates.len() {
        let mut visited_bonds = vec![false; molecule.bonds().count()];
        if !assign_distinct_wedge(
            stereo_index,
            &candidates,
            &mut bond_owners,
            &mut visited_bonds,
        ) {
            return Err(MoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: candidates[stereo_index].stereo_atom,
            });
        }
    }

    let mut selected = vec![None; molecule.bonds().count()];
    for (bond_index, owner) in bond_owners.into_iter().enumerate() {
        let Some(stereo_index) = owner else {
            continue;
        };
        let stereo = &candidates[stereo_index];
        let wedge = stereo
            .wedges
            .iter()
            .find(|wedge| wedge.bond.index() == bond_index)
            .expect("wedge assignment retains one candidate for its owning stereo atom");
        selected[bond_index] = Some(SelectedWedge {
            stereo_atom: stereo.stereo_atom,
            tip: stereo.site,
            base: wedge.base,
            kind: wedge.kind,
        });
    }
    Ok(selected)
}

fn tetrahedral_candidates(
    molecule: &Molecule,
    layout: &MoleculeLayout,
    stereo: StereoAtomView<'_>,
    tetrahedral_sites: &[bool],
) -> TetrahedralCandidates {
    let mut ligands = stereo
        .ligands()
        .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
        .collect::<Vec<_>>();
    ligands.sort_by_key(|ligand| {
        (
            ligand.kind != StereoLigandKind::Atom,
            ligand.atom_id,
            ligand.kind,
        )
    });
    let StereoCoset::Lit(coset) = stereo
        .coset_for(ligands.iter().copied())
        .expect("a closed stereo atom can be reframed over its stored distinct ligands")
    else {
        unreachable!("literal stereo configuration remains literal after reframing")
    };

    let mut wedges = ligands
        .iter()
        .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
        .filter_map(|ligand| {
            let bond = molecule.bonds().of(stereo.site_id(), ligand.atom_id)?;
            if bond.order().as_lit() != Some(1) {
                return None;
            }
            let kind = wedge_kind(layout, stereo.site_id(), &ligands, *ligand, coset)?;
            Some(WedgeCandidate {
                bond: bond.id,
                base: ligand.atom_id,
                kind,
            })
        })
        .collect::<Vec<_>>();
    wedges.sort_by_key(|wedge| (tetrahedral_sites[wedge.base.index()], wedge.bond.index()));

    TetrahedralCandidates {
        stereo_atom: stereo.id,
        site: stereo.site_id(),
        wedges,
    }
}

fn assign_distinct_wedge(
    stereo_index: usize,
    candidates: &[TetrahedralCandidates],
    bond_owners: &mut [Option<usize>],
    visited_bonds: &mut [bool],
) -> bool {
    for wedge in &candidates[stereo_index].wedges {
        let bond_index = wedge.bond.index();
        if visited_bonds[bond_index] {
            continue;
        }
        visited_bonds[bond_index] = true;
        let displaced = bond_owners[bond_index];
        if displaced.is_none_or(|owner| {
            assign_distinct_wedge(owner, candidates, bond_owners, visited_bonds)
        }) {
            bond_owners[bond_index] = Some(stereo_index);
            return true;
        }
    }
    false
}

fn wedge_kind(
    layout: &MoleculeLayout,
    site: AtomId,
    ligands: &[StereoLigand],
    wedged: StereoLigand,
    coset: u32,
) -> Option<WedgeKind> {
    [WedgeKind::Solid, WedgeKind::Hashed]
        .into_iter()
        .find(|&kind| wedge_coset(layout, site, ligands, wedged, kind) == Some(coset))
}

fn wedge_coset(
    layout: &MoleculeLayout,
    site: AtomId,
    ligands: &[StereoLigand],
    wedged: StereoLigand,
    kind: WedgeKind,
) -> Option<u32> {
    let center = point3(position(layout, site), 0.0);
    let actual_positions = ligands
        .iter()
        .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
        .map(|ligand| point3(position(layout, ligand.atom_id), 0.0))
        .collect::<Vec<_>>();
    let virtual_position = complementary_direction(center, &actual_positions);
    let wedge_z = match kind {
        WedgeKind::Solid => 1.0,
        WedgeKind::Hashed => -1.0,
    };
    let points = ligands
        .iter()
        .map(|ligand| match ligand.kind {
            StereoLigandKind::Atom if *ligand == wedged => {
                point3(position(layout, ligand.atom_id), wedge_z)
            }
            StereoLigandKind::Atom => point3(position(layout, ligand.atom_id), 0.0),
            StereoLigandKind::ImplicitHydrogen | StereoLigandKind::LonePair => virtual_position,
        })
        .collect::<Vec<_>>();
    let [first, second, third, fourth] = points.as_slice() else {
        return None;
    };
    let volume = signed_volume(*first, *second, *third, *fourth);
    if !volume.is_finite() || volume == 0.0 {
        return None;
    }
    Some(u32::from(volume >= 0.0))
}

fn point3(point: Point2D, z: f64) -> Point3D {
    Point3D::new(point.x, point.y, z)
}

fn atom_label(atom: AtomView<'_>) -> Option<AtomLabel> {
    let element = atom.element().as_lit()?;
    let isotope = match atom.isotope_mass().as_lit() {
        Some(IsotopeMass::MassNumber(mass)) => Some(mass),
        _ => None,
    };
    let charge = atom.charge().as_lit().filter(|&charge| charge != 0);
    let unpaired_electrons = atom
        .unpaired_electrons()
        .count
        .as_lit()
        .filter(|&count| count > 0);
    let isolated_carbon = element == Element::C && atom.neighbors().len() == 0;
    let decorated = isotope.is_some() || charge.is_some() || unpaired_electrons.is_some();
    if element == Element::C && !isolated_carbon && !decorated {
        return None;
    }

    let mut base = element.symbol().to_owned();
    let mut right_subscript = None;

    if !isolated_carbon || decorated {
        if let Some(hydrogens) = atom
            .implicit_hydrogens()
            .as_lit()
            .filter(|&count| count != 0)
        {
            base.push('H');
            if hydrogens != 1 {
                right_subscript = Some(hydrogens.to_string());
            }
        }
    }

    let mut right_superscript = String::new();
    append_charge_and_unpaired(&mut right_superscript, charge, unpaired_electrons);

    Some(AtomLabel {
        base,
        left_superscript: isotope.map(|mass| mass.to_string()),
        right_subscript,
        right_superscript: (!right_superscript.is_empty()).then_some(right_superscript),
    })
}

fn position(layout: &MoleculeLayout, atom: AtomId) -> Point2D {
    *layout
        .position(atom)
        .expect("frame agreement establishes every graph-IR atom position")
}

fn midpoint(first: Point2D, second: Point2D) -> Point2D {
    Point2D::new((first.x + second.x) / 2.0, (first.y + second.y) / 2.0)
}

fn aromatic_contour(
    molecule: &Molecule,
    layout: &MoleculeLayout,
    system: AromaticSystemView<'_>,
) -> Option<DashedContourItem> {
    let atoms = system.atom_ids().collect::<Vec<_>>();
    if atoms.len() < 3 {
        return None;
    }

    let mut member_indices = vec![None; molecule.atoms().count()];
    for (member_index, atom) in atoms.iter().copied().enumerate() {
        member_indices[atom.index()] = Some(member_index);
    }

    let mut edges = Vec::new();
    let mut adjacency = vec![Vec::new(); atoms.len()];
    for bond in system.bonds() {
        let [first, second] = bond.atom_ids();
        let first = member_indices[first.index()]?;
        let second = member_indices[second.index()]?;
        if first == second {
            return None;
        }
        edges.push((first, second));
        adjacency[first].push(second);
        adjacency[second].push(first);
    }

    if adjacency.iter().any(|neighbors| neighbors.len() < 2)
        || adjacency.iter().all(|neighbors| neighbors.len() != 2)
        || !is_connected(&adjacency)
    {
        return None;
    }

    let positions = atoms
        .iter()
        .map(|&atom| position(layout, atom))
        .collect::<Vec<_>>();
    if edges_cross(&edges, &positions) {
        return None;
    }
    for (site, neighbors) in adjacency.iter_mut().enumerate() {
        neighbors.sort_by(|&first, &second| {
            direction_angle(positions[site], positions[first])
                .total_cmp(&direction_angle(positions[site], positions[second]))
        });
        if degenerate_rotation(site, neighbors, &positions) {
            return None;
        }
    }

    let faces = walk_faces(&adjacency)?;
    let outer = faces
        .into_iter()
        .map(|face| {
            let area = signed_area_indices(&face, &positions);
            (area, face)
        })
        .filter(|(area, _)| *area < -GEOMETRY_EPSILON)
        .min_by(|(first, _), (second, _)| first.total_cmp(second))?
        .1;
    let outer = simplify_boundary(outer, &positions)?;
    let boundary = outer
        .iter()
        .map(|&member| positions[member])
        .collect::<Vec<_>>();
    let points = offset_clockwise_boundary(&boundary)?;

    Some(DashedContourItem {
        points,
        closed: true,
        references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
            system.id,
        ))],
    })
}

fn aromatic_annotation(
    molecule: &Molecule,
    layout: &MoleculeLayout,
    system: AromaticSystemView<'_>,
    contour: &[Point2D],
) -> Option<TextItem> {
    let charge = system.charge().as_lit().filter(|&charge| charge != 0);
    let unpaired_electrons = system
        .unpaired_electrons()
        .count
        .as_lit()
        .filter(|&count| count > 0);
    if charge.is_none() && unpaired_electrons.is_none() {
        return None;
    }

    let member_positions = system
        .atom_ids()
        .map(|atom| position(layout, atom))
        .collect::<Vec<_>>();
    let centroid = scale(
        member_positions
            .iter()
            .copied()
            .fold(Point2D::new(0.0, 0.0), add),
        1.0 / member_positions.len() as f64,
    );
    let position = annotation_candidates(centroid)
        .into_iter()
        .find(|&candidate| {
            point_in_polygon(candidate, contour)
                && annotation_position_is_clear(candidate, molecule, layout, contour)
        })
        .unwrap_or_else(|| {
            let max_x = contour
                .iter()
                .map(|point| point.x)
                .fold(f64::NEG_INFINITY, f64::max);
            let max_y = contour
                .iter()
                .map(|point| point.y)
                .fold(f64::NEG_INFINITY, f64::max);
            Point2D::new(
                max_x + AROMATIC_ANNOTATION_EXTERIOR_OFFSET,
                max_y + AROMATIC_ANNOTATION_EXTERIOR_OFFSET,
            )
        });

    let mut text = String::new();
    append_charge_and_unpaired(&mut text, charge, unpaired_electrons);
    Some(TextItem {
        position,
        text,
        references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
            system.id,
        ))],
    })
}

fn append_charge_and_unpaired(
    output: &mut String,
    charge: Option<i64>,
    unpaired_electrons: Option<i64>,
) {
    if let Some(charge) = charge {
        let magnitude = charge.unsigned_abs();
        if magnitude != 1 {
            write!(output, "{magnitude}").expect("writing to a String cannot fail");
        }
        output.push(if charge > 0 { '+' } else { '−' });
    }
    if let Some(count) = unpaired_electrons {
        for _ in 0..count {
            output.push('•');
        }
    }
}

fn annotation_candidates(centroid: Point2D) -> Vec<Point2D> {
    const OFFSETS: [(i8, i8); 17] = [
        (0, 0),
        (0, 1),
        (1, 0),
        (0, -1),
        (-1, 0),
        (1, 1),
        (1, -1),
        (-1, -1),
        (-1, 1),
        (0, 2),
        (2, 0),
        (0, -2),
        (-2, 0),
        (2, 2),
        (2, -2),
        (-2, -2),
        (-2, 2),
    ];
    OFFSETS
        .into_iter()
        .map(|(x, y)| {
            Point2D::new(
                centroid.x + f64::from(x) * AROMATIC_ANNOTATION_STEP,
                centroid.y + f64::from(y) * AROMATIC_ANNOTATION_STEP,
            )
        })
        .collect()
}

fn annotation_position_is_clear(
    candidate: Point2D,
    molecule: &Molecule,
    layout: &MoleculeLayout,
    contour: &[Point2D],
) -> bool {
    let clearance_squared = AROMATIC_ANNOTATION_CLEARANCE.powi(2);
    if molecule.atoms().iter().any(|atom| {
        squared_length(difference(candidate, position(layout, atom.id))) < clearance_squared
    }) {
        return false;
    }
    if molecule.bonds().iter().any(|bond| {
        let [first, second] = bond.atom_ids();
        squared_distance_to_segment(candidate, position(layout, first), position(layout, second))
            < clearance_squared
    }) {
        return false;
    }
    !(0..contour.len()).any(|index| {
        squared_distance_to_segment(
            candidate,
            contour[index],
            contour[(index + 1) % contour.len()],
        ) < clearance_squared
    })
}

fn squared_distance_to_segment(point: Point2D, start: Point2D, end: Point2D) -> f64 {
    let segment = difference(end, start);
    let length_squared = squared_length(segment);
    if length_squared <= GEOMETRY_EPSILON {
        return squared_length(difference(point, start));
    }
    let projection = (dot(difference(point, start), segment) / length_squared).clamp(0.0, 1.0);
    squared_length(difference(point, add(start, scale(segment, projection))))
}

fn is_connected(adjacency: &[Vec<usize>]) -> bool {
    let mut visited = vec![false; adjacency.len()];
    let mut pending = vec![0];
    visited[0] = true;
    while let Some(site) = pending.pop() {
        for &neighbor in &adjacency[site] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                pending.push(neighbor);
            }
        }
    }
    visited.into_iter().all(|is_visited| is_visited)
}

fn direction_angle(origin: Point2D, target: Point2D) -> f64 {
    (target.y - origin.y).atan2(target.x - origin.x)
}

fn degenerate_rotation(site: usize, neighbors: &[usize], positions: &[Point2D]) -> bool {
    neighbors.iter().enumerate().any(|(index, &first)| {
        let second = neighbors[(index + 1) % neighbors.len()];
        same_direction(positions[site], positions[first], positions[second])
    })
}

fn same_direction(origin: Point2D, first: Point2D, second: Point2D) -> bool {
    let first = difference(first, origin);
    let second = difference(second, origin);
    let first_length = squared_length(first);
    let second_length = squared_length(second);
    if first_length <= GEOMETRY_EPSILON || second_length <= GEOMETRY_EPSILON {
        return true;
    }
    cross(first, second).abs() <= GEOMETRY_EPSILON * first_length.sqrt() * second_length.sqrt()
        && dot(first, second) > 0.0
}

fn edges_cross(edges: &[(usize, usize)], positions: &[Point2D]) -> bool {
    edges.iter().enumerate().any(|(index, &(a, b))| {
        edges[index + 1..].iter().any(|&(c, d)| {
            a != c
                && a != d
                && b != c
                && b != d
                && segments_intersect(positions[a], positions[b], positions[c], positions[d])
        })
    })
}

fn segments_intersect(a: Point2D, b: Point2D, c: Point2D, d: Point2D) -> bool {
    let ab_c = orientation(a, b, c);
    let ab_d = orientation(a, b, d);
    let cd_a = orientation(c, d, a);
    let cd_b = orientation(c, d, b);
    (opposite_signs(ab_c, ab_d) && opposite_signs(cd_a, cd_b))
        || (ab_c.abs() <= GEOMETRY_EPSILON && point_on_segment(c, a, b))
        || (ab_d.abs() <= GEOMETRY_EPSILON && point_on_segment(d, a, b))
        || (cd_a.abs() <= GEOMETRY_EPSILON && point_on_segment(a, c, d))
        || (cd_b.abs() <= GEOMETRY_EPSILON && point_on_segment(b, c, d))
}

fn opposite_signs(first: f64, second: f64) -> bool {
    (first > GEOMETRY_EPSILON && second < -GEOMETRY_EPSILON)
        || (first < -GEOMETRY_EPSILON && second > GEOMETRY_EPSILON)
}

fn point_on_segment(point: Point2D, start: Point2D, end: Point2D) -> bool {
    point.x >= start.x.min(end.x) - GEOMETRY_EPSILON
        && point.x <= start.x.max(end.x) + GEOMETRY_EPSILON
        && point.y >= start.y.min(end.y) - GEOMETRY_EPSILON
        && point.y <= start.y.max(end.y) + GEOMETRY_EPSILON
}

fn orientation(first: Point2D, second: Point2D, third: Point2D) -> f64 {
    cross(difference(second, first), difference(third, first))
}

fn walk_faces(adjacency: &[Vec<usize>]) -> Option<Vec<Vec<usize>>> {
    let mut visited = adjacency
        .iter()
        .map(|neighbors| vec![false; neighbors.len()])
        .collect::<Vec<_>>();
    let directed_edge_count = adjacency.iter().map(Vec::len).sum::<usize>();
    let mut faces = Vec::new();

    for first in 0..adjacency.len() {
        for (second_position, &second) in adjacency[first].iter().enumerate() {
            if visited[first][second_position] {
                continue;
            }
            let start = (first, second);
            let mut edge = start;
            let mut face = Vec::new();
            loop {
                let outgoing = adjacency[edge.0]
                    .iter()
                    .position(|&neighbor| neighbor == edge.1)?;
                if visited[edge.0][outgoing] {
                    return None;
                }
                visited[edge.0][outgoing] = true;
                face.push(edge.0);
                let incoming = adjacency[edge.1]
                    .iter()
                    .position(|&neighbor| neighbor == edge.0)?;
                let next = adjacency[edge.1]
                    [(incoming + adjacency[edge.1].len() - 1) % adjacency[edge.1].len()];
                edge = (edge.1, next);
                if edge == start {
                    break;
                }
                if face.len() > directed_edge_count {
                    return None;
                }
            }
            let mut distinct = face.clone();
            distinct.sort_unstable();
            distinct.dedup();
            if face.len() < 3 || distinct.len() != face.len() {
                return None;
            }
            faces.push(face);
        }
    }
    Some(faces)
}

fn signed_area_indices(indices: &[usize], positions: &[Point2D]) -> f64 {
    let points = indices
        .iter()
        .map(|&index| positions[index])
        .collect::<Vec<_>>();
    signed_area(&points)
}

fn simplify_boundary(mut boundary: Vec<usize>, positions: &[Point2D]) -> Option<Vec<usize>> {
    loop {
        if boundary.len() < 3 {
            return None;
        }
        let removable = (0..boundary.len()).find(|&index| {
            let previous = positions[boundary[(index + boundary.len() - 1) % boundary.len()]];
            let current = positions[boundary[index]];
            let next = positions[boundary[(index + 1) % boundary.len()]];
            let incoming = difference(current, previous);
            let outgoing = difference(next, current);
            cross(incoming, outgoing).abs() <= GEOMETRY_EPSILON && dot(incoming, outgoing) > 0.0
        });
        let Some(index) = removable else {
            break;
        };
        boundary.remove(index);
    }

    if signed_area_indices(&boundary, positions) >= -GEOMETRY_EPSILON {
        return None;
    }
    let first = boundary
        .iter()
        .enumerate()
        .min_by_key(|(_, member)| **member)
        .map(|(index, _)| index)?;
    boundary.rotate_left(first);
    Some(boundary)
}

fn offset_clockwise_boundary(boundary: &[Point2D]) -> Option<Vec<Point2D>> {
    let mut points = Vec::with_capacity(boundary.len());
    for index in 0..boundary.len() {
        let previous = boundary[(index + boundary.len() - 1) % boundary.len()];
        let current = boundary[index];
        let next = boundary[(index + 1) % boundary.len()];
        let incoming = unit(difference(current, previous))?;
        let outgoing = unit(difference(next, current))?;
        let incoming_anchor = add(
            current,
            scale(right_normal(incoming), AROMATIC_CONTOUR_OFFSET),
        );
        let outgoing_anchor = add(
            current,
            scale(right_normal(outgoing), AROMATIC_CONTOUR_OFFSET),
        );
        let offset = line_intersection(incoming_anchor, incoming, outgoing_anchor, outgoing)?;
        if squared_length(difference(offset, current)).sqrt()
            > MAX_CONTOUR_MITER * AROMATIC_CONTOUR_OFFSET
        {
            return None;
        }
        points.push(offset);
    }

    if signed_area(&points) >= -GEOMETRY_EPSILON || polygon_self_intersects(&points) {
        return None;
    }
    for index in 0..points.len() {
        let midpoint = midpoint(points[index], points[(index + 1) % points.len()]);
        if !point_in_polygon(points[index], boundary) || !point_in_polygon(midpoint, boundary) {
            return None;
        }
    }
    Some(points)
}

fn polygon_self_intersects(points: &[Point2D]) -> bool {
    (0..points.len()).any(|first| {
        let first_next = (first + 1) % points.len();
        ((first + 1)..points.len()).any(|second| {
            let second_next = (second + 1) % points.len();
            first != second_next
                && first_next != second
                && segments_intersect(
                    points[first],
                    points[first_next],
                    points[second],
                    points[second_next],
                )
        })
    })
}

fn point_in_polygon(point: Point2D, polygon: &[Point2D]) -> bool {
    let mut inside = false;
    for index in 0..polygon.len() {
        let first = polygon[index];
        let second = polygon[(index + 1) % polygon.len()];
        if orientation(first, second, point).abs() <= GEOMETRY_EPSILON
            && point_on_segment(point, first, second)
        {
            return true;
        }
        if (first.y > point.y) != (second.y > point.y) {
            let x = (second.x - first.x) * (point.y - first.y) / (second.y - first.y) + first.x;
            if point.x < x {
                inside = !inside;
            }
        }
    }
    inside
}

fn signed_area(points: &[Point2D]) -> f64 {
    points
        .iter()
        .enumerate()
        .map(|(index, &point)| {
            let next = points[(index + 1) % points.len()];
            point.x * next.y - next.x * point.y
        })
        .sum::<f64>()
        / 2.0
}

fn line_intersection(
    first: Point2D,
    first_direction: Point2D,
    second: Point2D,
    second_direction: Point2D,
) -> Option<Point2D> {
    let denominator = cross(first_direction, second_direction);
    if denominator.abs() <= GEOMETRY_EPSILON {
        return None;
    }
    let distance = cross(difference(second, first), second_direction) / denominator;
    Some(add(first, scale(first_direction, distance)))
}

fn unit(vector: Point2D) -> Option<Point2D> {
    let length = squared_length(vector).sqrt();
    (length > GEOMETRY_EPSILON).then(|| scale(vector, 1.0 / length))
}

fn right_normal(vector: Point2D) -> Point2D {
    Point2D::new(vector.y, -vector.x)
}

fn difference(first: Point2D, second: Point2D) -> Point2D {
    Point2D::new(first.x - second.x, first.y - second.y)
}

fn add(first: Point2D, second: Point2D) -> Point2D {
    Point2D::new(first.x + second.x, first.y + second.y)
}

fn scale(vector: Point2D, factor: f64) -> Point2D {
    Point2D::new(vector.x * factor, vector.y * factor)
}

fn dot(first: Point2D, second: Point2D) -> f64 {
    first.x * second.x + first.y * second.y
}

fn cross(first: Point2D, second: Point2D) -> f64 {
    first.x * second.y - first.y * second.x
}

fn squared_length(vector: Point2D) -> f64 {
    dot(vector, vector)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AromaticSystemForm, AromaticSystemId, AtomForm, AtomId, BondForm, BondId,
        ElectronCountsForm, FrameTransport, MoleculeEntries, StereoAtomForm, StereoAtomId,
        TryIntoIr,
    };
    use umol_graph_ir::mol_dsl;
    use umol_perm::Permutation;

    use super::*;
    use crate::ctfile::parser::parse_mol_to_table_ir;
    use crate::depict::Bounds;
    #[cfg(feature = "coordgen")]
    use crate::depict::{Depict, DepictConfig};
    #[cfg(feature = "coordgen")]
    use crate::layout::{layout_molecule, MoleculeLayoutAlgorithm};

    fn layout(positions: &[[f64; 2]]) -> MoleculeLayout {
        MoleculeLayout::try_new(positions.iter().map(|&[x, y]| Point2D::new(x, y)).collect())
            .unwrap()
    }

    #[rstest]
    fn test_depict_literal_atom_and_bond_projection() {
        let molecule = mol_dsl!(r#"{:atoms ["C#i13#c+#h2#u2" "O#c-"] :bonds [[0 1 "2"]]}"#);
        let layout = layout(&[[0.0, 1.0], [2.0, -1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            depiction.items(),
            [
                DepictionItem::Bond(BondItem {
                    start: Point2D::new(0.0, 1.0),
                    end: Point2D::new(2.0, -1.0),
                    line_count: 2,
                    references: vec![DepictionReference::Molecule(Entity::Bond(BondId(0)))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(0.0, 1.0),
                    label: AtomLabel {
                        base: "CH".to_owned(),
                        left_superscript: Some("13".to_owned()),
                        right_subscript: Some("2".to_owned()),
                        right_superscript: Some("+••".to_owned()),
                    },
                    references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(2.0, -1.0),
                    label: AtomLabel {
                        base: "O".to_owned(),
                        left_superscript: None,
                        right_subscript: None,
                        right_superscript: Some("−".to_owned()),
                    },
                    references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
                }),
            ]
        );
        assert_eq!(
            depiction.bounds(),
            Some(&Bounds {
                min: Point2D::new(0.0, -1.0),
                max: Point2D::new(2.0, 1.0),
            })
        );
    }

    #[rstest]
    #[case::heteroatom(
        r#"{:atoms ["N#h2"] :bonds []}"#,
        Some(AtomLabel {
            base: "NH".to_owned(),
            left_superscript: None,
            right_subscript: Some("2".to_owned()),
            right_superscript: None,
        })
    )]
    #[case::isolated_carbon(
        r#"{:atoms ["C#h4"] :bonds []}"#,
        Some(AtomLabel {
            base: "C".to_owned(),
            left_superscript: None,
            right_subscript: None,
            right_superscript: None,
        })
    )]
    #[case::isolated_charged_carbon(
        r#"{:atoms ["C#c+#h4"] :bonds []}"#,
        Some(AtomLabel {
            base: "CH".to_owned(),
            left_superscript: None,
            right_subscript: Some("4".to_owned()),
            right_superscript: Some("+".to_owned()),
        })
    )]
    #[case::skeleton_carbon(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#, None)]
    #[case::skeleton_carbon_with_implicit_hydrogens(
        r#"{:atoms ["C#h3" "C"] :bonds [[0 1 "1"]]}"#,
        None
    )]
    #[case::isotopic_carbon(
        r#"{:atoms ["C#i13#h2" "C"] :bonds [[0 1 "1"]]}"#,
        Some(AtomLabel {
            base: "CH".to_owned(),
            left_superscript: Some("13".to_owned()),
            right_subscript: Some("2".to_owned()),
            right_superscript: None,
        })
    )]
    #[case::charged_carbon(
        r#"{:atoms ["C#c+" "C"] :bonds [[0 1 "1"]]}"#,
        Some(AtomLabel {
            base: "C".to_owned(),
            left_superscript: None,
            right_subscript: None,
            right_superscript: Some("+".to_owned()),
        })
    )]
    #[case::radical_carbon(
        r#"{:atoms ["C#u2" "C"] :bonds [[0 1 "1"]]}"#,
        Some(AtomLabel {
            base: "C".to_owned(),
            left_superscript: None,
            right_subscript: None,
            right_superscript: Some("••".to_owned()),
        })
    )]
    #[case::multi_digit_fields(
        r#"{:atoms ["N#i15#c-12#h12#u2"] :bonds []}"#,
        Some(AtomLabel {
            base: "NH".to_owned(),
            left_superscript: Some("15".to_owned()),
            right_subscript: Some("12".to_owned()),
            right_superscript: Some("12−••".to_owned()),
        })
    )]
    #[case::independently_nonliteral_fields(
        r#"{:atoms ["N#i*#c*#h*#u*#s3"] :bonds []}"#,
        Some(AtomLabel {
            base: "N".to_owned(),
            left_superscript: None,
            right_subscript: None,
            right_superscript: None,
        })
    )]
    fn test_atom_label(#[case] input: &str, #[case] expected: Option<AtomLabel>) {
        let molecule = mol_dsl!(input);

        assert_eq!(atom_label(molecule.atom(AtomId(0))), expected);
    }

    #[test]
    fn test_depict_projected_bond_order_changes_output() {
        let single = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let double = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0]]);

        assert_eq!(bond_line_counts(&depict(&single, &layout).unwrap()), [1]);
        assert_eq!(bond_line_counts(&depict(&double, &layout).unwrap()), [2]);
    }

    #[rstest]
    #[case::solid("Th0", WedgeKind::Solid)]
    #[case::hashed("Th1", WedgeKind::Hashed)]
    fn test_depict_tetrahedral_wedge(#[case] attributes: &str, #[case] kind: WedgeKind) {
        let molecule = mol_dsl!(&format!(
            r#"{{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "{attributes}"}}]}}"#
        ));
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();
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
            [WedgeItem {
                tip: Point2D::new(0.0, 0.0),
                base: Point2D::new(1.0, 0.0),
                kind,
                references: vec![
                    DepictionReference::Molecule(Entity::Bond(BondId(0))),
                    DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(0))),
                ],
            }]
        );
        assert_eq!(bond_line_counts(&depiction), [1, 1, 1]);
    }

    #[rstest]
    #[case::coset_zero(0)]
    #[case::coset_one(1)]
    fn test_depict_tetrahedral_reframing(#[case] coset: u32) {
        let display_ligands = (1..=4)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]);
        let baseline = Molecule::from_entries(MoleculeEntries {
            atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: (1..=4)
                .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                display_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, coset),
            )],
            ..Default::default()
        });
        let expected = depict(&baseline, &layout).unwrap();

        for permutation_rank in 0..24 {
            let stored_ligands = Permutation::unrank(4, permutation_rank).act(&display_ligands);
            let action = Permutation::between(&display_ligands, &stored_ligands).unwrap();
            let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset)
                .reframe_by(&action)
                .unwrap();
            let reframed = Molecule::from_entries(MoleculeEntries {
                atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                    .into_iter()
                    .map(AtomForm::from_element)
                    .collect(),
                bonds: (1..=4)
                    .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                    .collect(),
                stereo_atoms: vec![(AtomId(0), stored_ligands, attributes)],
                ..Default::default()
            });

            assert_eq!(
                reframed
                    .stereo_atoms()
                    .iter()
                    .next()
                    .unwrap()
                    .coset_for(display_ligands.iter().copied()),
                Some(StereoCoset::Lit(coset)),
            );
            assert_eq!(
                depict(&reframed, &layout).unwrap().items(),
                expected.items()
            );
        }
    }

    #[rstest]
    #[case::coset_zero(0)]
    #[case::coset_one(1)]
    fn test_depict_tetrahedral_geometry(#[case] coset: u32) {
        let display_ligands = (1..=4)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();
        let base_positions = [
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
            Point2D::new(-1.0, 0.0),
            Point2D::new(0.0, -1.0),
        ];

        for permutation_rank in 0..24 {
            let stored_ligands = Permutation::unrank(4, permutation_rank).act(&display_ligands);
            let action = Permutation::between(&display_ligands, &stored_ligands).unwrap();
            let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset)
                .reframe_by(&action)
                .unwrap();
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                    .into_iter()
                    .map(AtomForm::from_element)
                    .collect(),
                bonds: (1..=4)
                    .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                    .collect(),
                stereo_atoms: vec![(AtomId(0), stored_ligands, attributes)],
                ..Default::default()
            });

            for quarter_turns in 0..4 {
                for reflected in [false, true] {
                    let positions = base_positions
                        .iter()
                        .copied()
                        .map(|point| {
                            let mut x = if reflected { -point.x } else { point.x };
                            let mut y = point.y;
                            for _ in 0..quarter_turns {
                                (x, y) = (-y, x);
                            }
                            Point2D::new(x * 3.0 + 4.0, y * 3.0 - 5.0)
                        })
                        .collect();
                    let layout = MoleculeLayout::try_new(positions).unwrap();
                    let depiction = depict(&molecule, &layout).unwrap();
                    let emitted_wedges = wedges(&depiction);
                    let [wedge] = emitted_wedges.as_slice() else {
                        panic!("tetrahedral fixture must emit exactly one wedge");
                    };
                    let wedged_atom = wedge
                        .references
                        .iter()
                        .find_map(|reference| match reference {
                            DepictionReference::Molecule(Entity::Bond(bond)) => molecule
                                .bond(*bond)
                                .atom_ids()
                                .into_iter()
                                .find(|&atom| atom != AtomId(0)),
                            _ => None,
                        })
                        .unwrap();
                    let points = (1..=4)
                        .map(|atom| {
                            let atom = AtomId(atom);
                            let point = *layout.position(atom).unwrap();
                            let z = if atom == wedged_atom {
                                match wedge.kind {
                                    WedgeKind::Solid => 1.0,
                                    WedgeKind::Hashed => -1.0,
                                }
                            } else {
                                0.0
                            };
                            Point3D::new(point.x, point.y, z)
                        })
                        .collect::<Vec<_>>();
                    let [first, second, third, fourth] = points.as_slice() else {
                        unreachable!("tetrahedral fixture has four ligands")
                    };
                    let volume = signed_volume(*first, *second, *third, *fourth);

                    assert!(volume.is_finite() && volume != 0.0);
                    assert_eq!(u32::from(volume >= 0.0), coset);
                }
            }
        }
    }

    #[rstest]
    #[case::implicit_hydrogen("[:h 0]")]
    #[case::lone_pair("[:lp 0]")]
    fn test_depict_tetrahedral_virtual_ligand(#[case] virtual_ligand: &str) {
        let molecule = mol_dsl!(&format!(
            r#"{{:atoms ["C" "F" "Cl" "Br"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
                :stereo-atoms [{{:site 0 :ligands [1 2 3 {virtual_ligand}] :attrs "Th0"}}]}}"#
        ));
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            wedges(&depiction),
            [WedgeItem {
                tip: Point2D::new(0.0, 0.0),
                base: Point2D::new(1.0, 0.0),
                kind: WedgeKind::Solid,
                references: vec![
                    DepictionReference::Molecule(Entity::Bond(BondId(0))),
                    DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(0))),
                ],
            }]
        );
        assert_eq!(bond_line_counts(&depiction), [1, 1]);
    }

    #[rstest]
    fn test_depict_tetrahedral_ring_site() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "F" "Cl"]
                :bonds [[0 3 "1"] [0 1 "1"] [1 2 "1"] [2 0 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let layout = layout(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 0.866],
            [-1.0, 0.0],
            [0.0, -1.0],
        ]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(wedges(&depiction).len(), 1);
        assert_eq!(
            wedges(&depiction)[0].references,
            [
                DepictionReference::Molecule(Entity::Bond(BondId(0))),
                DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(0))),
            ]
        );
        assert_eq!(bond_line_counts(&depiction), [1, 1, 1, 1]);
    }

    #[rstest]
    fn test_depict_adjacent_tetrahedral_sites_use_distinct_bonds() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "F" "Cl" "Br" "F" "Cl" "Br"]
                :bonds [[0 1 "1"]
                        [0 2 "1"] [0 3 "1"] [0 4 "1"]
                        [1 5 "1"] [1 6 "1"] [1 7 "1"]]
                :stereo-atoms [
                    {:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                    {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]}"#
        );
        let layout = layout(&[
            [-0.5, 0.0],
            [0.5, 0.0],
            [-1.5, 0.0],
            [-0.5, 1.0],
            [-0.5, -1.0],
            [1.5, 0.0],
            [0.5, 1.0],
            [0.5, -1.0],
        ]);

        let depiction = depict(&molecule, &layout).unwrap();
        let wedge_references = wedges(&depiction)
            .into_iter()
            .map(|wedge| wedge.references)
            .collect::<Vec<_>>();

        assert_eq!(
            wedge_references,
            [
                vec![
                    DepictionReference::Molecule(Entity::Bond(BondId(1))),
                    DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(0))),
                ],
                vec![
                    DepictionReference::Molecule(Entity::Bond(BondId(4))),
                    DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(1))),
                ],
            ]
        );
        assert_eq!(bond_line_counts(&depiction), [1, 1, 1, 1, 1]);
    }

    #[rstest]
    fn test_depict_omits_nonliteral_projection_states() {
        let molecule = mol_dsl!(r#"{:atoms ["*#i13#c+#h2" "N#i*#c*#h*"] :bonds [[0 1 "*"]]}"#);
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            depiction.items(),
            [DepictionItem::Atom(AtomItem {
                position: Point2D::new(1.0, 0.0),
                label: AtomLabel {
                    base: "N".to_owned(),
                    left_superscript: None,
                    right_subscript: None,
                    right_superscript: None,
                },
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
            })]
        );
    }

    #[rstest]
    fn test_depict_single_aromatic_system_contour() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3] :attrs "*"}]}"#
        );
        let layout = layout(&[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        let contours = contours(&depiction);
        let [contour] = contours.as_slice() else {
            panic!("expected exactly one aromatic-system contour");
        };
        assert_points_close(
            &contour.points,
            &[
                Point2D::new(-0.82, -0.82),
                Point2D::new(-0.82, 0.82),
                Point2D::new(0.82, 0.82),
                Point2D::new(0.82, -0.82),
            ],
        );
        assert!(contour.closed);
        assert_eq!(
            contour.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))]
        );
    }

    #[rstest]
    fn test_depict_multiple_aromatic_system_contours() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]
                        [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 4 "1"]]
                :aromatic-systems [
                    {:atoms [0 1 2 3] :attrs "*"}
                    {:atoms [4 5 6 7] :attrs "*"}]}"#
        );
        let layout = layout(&[
            [-3.0, -1.0],
            [-1.0, -1.0],
            [-1.0, 1.0],
            [-3.0, 1.0],
            [1.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            [1.0, 1.0],
        ]);

        let contours = contours(&depict(&molecule, &layout).unwrap());
        let [first, second] = contours.as_slice() else {
            panic!("expected one contour for each aromatic system");
        };

        assert_points_close(
            &first.points,
            &[
                Point2D::new(-2.82, -0.82),
                Point2D::new(-2.82, 0.82),
                Point2D::new(-1.18, 0.82),
                Point2D::new(-1.18, -0.82),
            ],
        );
        assert_points_close(
            &second.points,
            &[
                Point2D::new(1.18, -0.82),
                Point2D::new(1.18, 0.82),
                Point2D::new(2.82, 0.82),
                Point2D::new(2.82, -0.82),
            ],
        );
        assert_eq!([first.closed, second.closed], [true, true]);
        assert_eq!(
            first.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))]
        );
        assert_eq!(
            second.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(1),
            ))]
        );
    }

    #[rstest]
    fn test_depict_concave_aromatic_system_contour() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "*"}]}"#
        );
        let layout = layout(&[
            [-2.0, -1.0],
            [2.0, -1.0],
            [2.0, 1.0],
            [0.0, 0.0],
            [-2.0, 1.0],
        ]);

        let contours = contours(&depict(&molecule, &layout).unwrap());
        let [contour] = contours.as_slice() else {
            panic!("expected exactly one concave aromatic-system contour");
        };

        assert_points_close(
            &contour.points,
            &[
                Point2D::new(-1.82, -0.82),
                Point2D::new(-1.82, 0.7087538820250189),
                Point2D::new(0.0, -0.20124611797498107),
                Point2D::new(1.82, 0.7087538820250189),
                Point2D::new(1.82, -0.82),
            ],
        );
        assert!(contour.closed);
        assert_eq!(
            contour.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))]
        );
    }

    #[rstest]
    #[case::crossed([[-1.0, -1.0], [1.0, 1.0], [-1.0, 1.0], [1.0, -1.0]])]
    #[case::degenerate([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])]
    fn test_depict_omits_untrustworthy_aromatic_contour(#[case] positions: [[f64; 2]; 4]) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3] :attrs "*"}]}"#
        );

        let depiction = depict(&molecule, &layout(&positions)).unwrap();

        assert_eq!(contours(&depiction), []);
    }

    #[rstest]
    fn test_depict_omits_aromatic_system_without_degree_two_member() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]
                        [1 2 "1"] [1 3 "1"] [2 3 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3] :attrs "*"}]}"#
        );
        let layout = layout(&[[0.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [0.0, 0.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(contours(&depiction), []);
    }

    #[rstest]
    fn test_depict_omits_c60_conformance_cage_contour() {
        let source = include_str!("../../tests/mol_parsing/data/molecule/scifinder/99685-96-8.mol");
        let table = parse_mol_to_table_ir(source).unwrap();
        let positions = table
            .positions
            .as_ref()
            .expect("the C60 conformance molecule carries coordinates")
            .iter()
            .map(|point| Point2D::new(point.x, point.y))
            .collect();
        let layout = MoleculeLayout::try_new(positions).unwrap();
        let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let mut editor = molecule.edit();
        editor.add_aromatic_system(
            (0..60).map(AtomId).collect(),
            AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(vec![1; 60]),
                ..Default::default()
            },
        );
        let molecule = editor.build();

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(molecule.atoms().count(), 60);
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).bond_count(),
            90
        );
        assert_eq!(contours(&depiction), []);
        assert_eq!(texts(&depiction), []);
    }

    #[rstest]
    #[case::charge_and_radicals("*#c-#u2#s3", Some("−••"))]
    #[case::literal_charge_nonliteral_radicals("*#c2#u*#s3", Some("2+"))]
    #[case::nonliteral_charge_literal_radicals("*#c*#u2#s*", Some("••"))]
    #[case::neutral_closed_shell("*#c0#u0#s1", None)]
    #[case::both_nonliteral("*#c*#u*#s*", None)]
    fn test_depict_aromatic_system_annotation_fields(
        #[case] attributes: &str,
        #[case] expected: Option<&str>,
    ) {
        let molecule = mol_dsl!(&format!(
            r#"{{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]
                :aromatic-systems [{{:atoms [0 1 2 3] :attrs "{attributes}"}}]}}"#
        ));
        let layout = layout(&[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();
        let expected = expected.map(|text| TextItem {
            position: Point2D::new(0.0, 0.0),
            text: text.to_owned(),
            references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))],
        });

        assert_eq!(texts(&depiction), expected.into_iter().collect::<Vec<_>>());
    }

    #[rstest]
    fn test_depict_moves_aromatic_annotation_off_internal_fusion_bond() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]
                        [5 0 "1"] [2 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"]
                        [9 3 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "*#c+"}]}"#
        );
        let layout = layout(&[
            [-1.732, 0.5],
            [-0.866, 1.0],
            [0.0, 0.5],
            [0.0, -0.5],
            [-0.866, -1.0],
            [-1.732, -0.5],
            [0.866, 1.0],
            [1.732, 0.5],
            [1.732, -0.5],
            [0.866, -1.0],
        ]);

        let depiction = depict(&molecule, &layout).unwrap();
        let contours = contours(&depiction);
        let [contour] = contours.as_slice() else {
            panic!("expected exactly one fused aromatic-system contour");
        };

        assert_points_close(
            &contour.points,
            &[
                Point2D::new(-1.552, 0.3960784758008787),
                Point2D::new(-0.866, 0.7921523788031881),
                Point2D::new(0.0, 0.2921523788031882),
                Point2D::new(0.866, 0.7921523788031881),
                Point2D::new(1.552, 0.3960784758008787),
                Point2D::new(1.552, -0.3960784758008787),
                Point2D::new(0.866, -0.7921523788031881),
                Point2D::new(0.0, -0.2921523788031882),
                Point2D::new(-0.866, -0.7921523788031881),
                Point2D::new(-1.552, -0.3960784758008787),
            ],
        );
        assert!(contour.closed);
        assert_eq!(
            contour.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))]
        );

        assert_eq!(
            texts(&depiction),
            [TextItem {
                position: Point2D::new(0.35, 0.0),
                text: "+".to_owned(),
                references: vec![DepictionReference::Molecule(Entity::AromaticSystem(
                    AromaticSystemId(0),
                ))],
            }]
        );
    }

    #[rstest]
    fn test_depict_uses_exterior_aromatic_annotation_fallback() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3] :attrs "*#c-"}]}"#
        );
        let layout = layout(&[[-0.2, -0.2], [0.2, -0.2], [0.2, 0.2], [-0.2, 0.2]]);

        let depiction = depict(&molecule, &layout).unwrap();
        let texts = texts(&depiction);
        let [annotation] = texts.as_slice() else {
            panic!("expected exactly one aromatic-system annotation");
        };

        assert_point_close(annotation.position, Point2D::new(0.37, 0.37));
        assert_eq!(annotation.text, "−");
        assert_eq!(
            annotation.references,
            [DepictionReference::Molecule(Entity::AromaticSystem(
                AromaticSystemId(0),
            ))]
        );
    }

    #[rstest]
    fn test_depict_local_aromatic_assertions_are_omitted() {
        let baseline = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let asserted = mol_dsl!(r#"{:atoms ["C#a+" "C#a+"] :bonds [[0 1 "1#a"]]}"#);
        let layout = layout(&[[0.0, 0.0], [2.0, 0.0]]);

        assert_ne!(baseline, asserted);
        assert_eq!(
            depict(&baseline, &layout).unwrap().items,
            depict(&asserted, &layout).unwrap().items
        );
    }

    #[rstest]
    #[case::undetermined("Th*")]
    #[case::set("Th{0,1}")]
    #[case::term("Th?configuration")]
    #[case::unsupported("Sp0")]
    fn test_depict_stereo_atom_omission(#[case] attributes: &str) {
        let molecule = mol_dsl!(&format!(
            r#"{{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{{:site 0 :ligands [1 2 3 4] :attrs "{attributes}"}}]}}"#
        ));
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert!(depiction
            .items()
            .iter()
            .all(|item| !matches!(item, DepictionItem::Wedge(_))));
        assert_eq!(bond_line_counts(&depiction), [1, 1, 1, 1]);
    }

    #[rstest]
    fn test_depict_cis_trans_stereo() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        );
        let layout = layout(&[[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, -1.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(bond_line_counts(&depiction), [1, 2, 1]);
    }

    #[rstest]
    #[case::lone_pairs(
        r#"{:atoms ["C"] :bonds []}"#,
        r#"{:atoms ["C#n1"] :bonds []}"#,
        vec![[0.0, 0.0]]
    )]
    #[case::dative(
        r#"{:atoms ["C" "N"] :bonds []}"#,
        r#"{:atoms ["C" "N"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs :single}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0]]
    )]
    #[case::multicenter(
        r#"{:atoms ["C" "C" "C"] :bonds []}"#,
        r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]
    )]
    #[case::noncovalent(
        r#"{:atoms ["N" "H"] :bonds []}"#,
        r#"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0]]
    )]
    #[case::constraint(
        r#"{:atoms ["C" "C"] :bonds []}"#,
        r#"{:atoms ["C" "C"] :bonds [] :constraints [{:connected {:atoms [0 1]}}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0]]
    )]
    fn test_depict_out_of_projection_difference_is_omitted(
        #[case] baseline: &str,
        #[case] changed: &str,
        #[case] positions: Vec<[f64; 2]>,
    ) {
        let baseline = mol_dsl!(baseline);
        let changed = mol_dsl!(changed);
        let layout = layout(&positions);

        assert_ne!(baseline, changed);
        assert_eq!(
            depict(&baseline, &layout).unwrap().items,
            depict(&changed, &layout).unwrap().items
        );
    }

    #[rstest]
    fn test_depict_tetrahedral_geometry_error() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]);

        assert_eq!(
            depict(&molecule, &layout).err(),
            Some(MoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: StereoAtomId(0),
            })
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    fn test_molecule_depict() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);

        assert_eq!(
            molecule.depict().unwrap().render_svg(),
            molecule
                .depict_with(&DepictConfig::default())
                .unwrap()
                .render_svg()
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
    fn test_molecule_depict_with(#[case] algorithm: MoleculeLayoutAlgorithm) {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);
        let layout = layout_molecule(&molecule, algorithm).unwrap();
        let expected = depict(&molecule, &layout).unwrap();
        let config = DepictConfig {
            layout_algorithm: algorithm,
        };

        assert_eq!(molecule.depict_with(&config).unwrap().items, expected.items);
    }

    fn bond_line_counts(depiction: &Depiction) -> Vec<u8> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Bond(bond) => Some(bond.line_count),
                _ => None,
            })
            .collect()
    }

    fn wedges(depiction: &Depiction) -> Vec<WedgeItem> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Wedge(wedge) => Some(wedge.clone()),
                _ => None,
            })
            .collect()
    }

    fn contours(depiction: &Depiction) -> Vec<DashedContourItem> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::DashedContour(contour) => Some(contour.clone()),
                _ => None,
            })
            .collect()
    }

    fn texts(depiction: &Depiction) -> Vec<TextItem> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Text(text) => Some(text.clone()),
                _ => None,
            })
            .collect()
    }

    fn assert_points_close(actual: &[Point2D], expected: &[Point2D]) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert!((actual.x - expected.x).abs() < 1.0e-12);
            assert!((actual.y - expected.y).abs() < 1.0e-12);
        }
    }

    fn assert_point_close(actual: Point2D, expected: Point2D) {
        assert!((actual.x - expected.x).abs() < 1.0e-12);
        assert!((actual.y - expected.y).abs() < 1.0e-12);
    }
}
