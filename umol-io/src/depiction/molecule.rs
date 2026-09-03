//! Depiction construction from graph-IR molecules and supplied layouts.

use std::any::Any;
use std::fmt::Write;

use thiserror::Error;
use umol_chem::element::Element;
use umol_geometric_core::{complementary_direction, signed_volume, Point2D, Point3D};
use umol_graph_ir::ir::{
    AsLit, AtomId, AtomView, BondId, Entity, IsotopeMass, Molecule, StereoAtomId, StereoAtomView,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
};
use umol_utils::error::UmolError;

#[cfg(feature = "coordgen")]
use super::Depict;
use super::{
    AtomItem, BondItem, Depiction, DepictionItem, DepictionReference, MarkerItem, MarkerKind,
    WedgeItem, WedgeKind,
};
#[cfg(feature = "coordgen")]
use crate::layout::{layout_molecule, LayoutError, MoleculeLayoutAlgorithm};
use crate::layout::{MoleculeLayout, MoleculeLayoutError};

/// Constructs the first format-neutral depiction projection of `molecule` in `layout`.
///
/// Localized bonds and selected tetrahedral wedges are followed by visible atom labels and
/// aromatic markers, with each group ordered by graph-IR id. Carbon labels are omitted at
/// non-isolated skeleton vertices unless an isotope, charge, or radical count decorates the atom.
/// Aromatic markers project both aromatic-system membership and definite aromatic atom or bond
/// constraints. Definite cis/trans stereo is carried by the supplied coordinates. Nonliteral
/// projected fields and unsupported overlays, stereo kinds, or constraints are omitted. The first
/// projection does not represent dative, multicenter, or noncovalent bonds, unprojected inherent
/// fields, or unsupported constraints.
///
/// # Errors
///
/// Returns [`MoleculeDepictionError::LayoutFrame`] if the molecule and layout do not use the same
/// dense atom frame, or [`MoleculeDepictionError::TetrahedralGeometry`] if a definite tetrahedral
/// stereo atom cannot be represented by a distinct, geometrically valid display wedge.
///
/// # Semantic properties
///
/// Every emitted tetrahedral wedge replaces its selected localized single bond. Reading the wedge
/// with the TableIR winding convention in the selected ligand frame reproduces the stored coset.
pub fn depict(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<Depiction, MoleculeDepictionError> {
    layout.check_frame(molecule)?;

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

    for atom in molecule.atoms().iter() {
        let system = atom.aromatic_system_id();
        let asserted = atom
            .constraints()
            .aromatic_valence()
            .is_some_and(|value| value.is_aromatic());
        if system.is_none() && !asserted {
            continue;
        }

        let mut references = Vec::with_capacity(2);
        if let Some(system) = system {
            references.push(DepictionReference::Molecule(Entity::AromaticSystem(system)));
        }
        references.push(DepictionReference::Molecule(Entity::Atom(atom.id)));
        items.push(DepictionItem::Marker(MarkerItem {
            position: position(layout, atom.id),
            kind: MarkerKind::Aromatic,
            references,
        }));
    }

    for bond in molecule.bonds().iter() {
        let system = bond.aromatic_system().map(|system| system.id);
        let asserted = bond.constraints().aromatic().as_lit() == Some(true);
        if system.is_none() && !asserted {
            continue;
        }

        let [first, second] = bond.atom_ids();
        let mut references = Vec::with_capacity(2);
        if let Some(system) = system {
            references.push(DepictionReference::Molecule(Entity::AromaticSystem(system)));
        }
        references.push(DepictionReference::Molecule(Entity::Bond(bond.id)));
        items.push(DepictionItem::Marker(MarkerItem {
            position: midpoint(position(layout, first), position(layout, second)),
            kind: MarkerKind::Aromatic,
            references,
        }));
    }

    Ok(Depiction::from_items(items))
}

#[cfg(feature = "coordgen")]
impl Depict for Molecule {
    type Error = MoleculeDepictionError;

    fn depict_with(
        &self,
        layout_algorithm: MoleculeLayoutAlgorithm,
    ) -> Result<Depiction, Self::Error> {
        let layout = layout_molecule(self, layout_algorithm)?;
        depict(self, &layout)
    }
}

/// Failures while depicting a graph-IR [`Molecule`].
#[derive(Clone, Debug, Error, PartialEq)]
pub enum MoleculeDepictionError {
    #[cfg(feature = "coordgen")]
    #[error("layout: {0}")]
    Layout(#[from] LayoutError),
    #[error("layout frame: {0}")]
    LayoutFrame(#[from] MoleculeLayoutError),
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

fn atom_label(atom: AtomView<'_>) -> Option<String> {
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

    let mut label = String::new();

    if let Some(mass) = isotope {
        write!(label, "{mass}").expect("writing to a String cannot fail");
    }
    label.push_str(element.symbol());

    if !isolated_carbon || decorated {
        if let Some(hydrogens) = atom
            .implicit_hydrogens()
            .as_lit()
            .filter(|&count| count != 0)
        {
            label.push('H');
            if hydrogens != 1 {
                write!(label, "{hydrogens}").expect("writing to a String cannot fail");
            }
        }
    }

    if let Some(charge) = charge {
        let magnitude = charge.unsigned_abs();
        if magnitude != 1 {
            write!(label, "{magnitude}").expect("writing to a String cannot fail");
        }
        label.push(if charge > 0 { '+' } else { '-' });
    }

    if let Some(count) = unpaired_electrons {
        for _ in 0..count {
            label.push('•');
        }
    }

    Some(label)
}

fn position(layout: &MoleculeLayout, atom: AtomId) -> Point2D {
    *layout
        .position(atom)
        .expect("frame agreement establishes every graph-IR atom position")
}

fn midpoint(first: Point2D, second: Point2D) -> Point2D {
    Point2D::new((first.x + second.x) / 2.0, (first.y + second.y) / 2.0)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, StereoAtomId};
    use umol_graph_ir::mol_dsl;

    use super::*;
    #[cfg(feature = "coordgen")]
    use crate::depiction::Depict;
    use crate::depiction::{Bounds, MarkerItem};
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
                    label: "13CH2+••".to_owned(),
                    references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(0)))],
                }),
                DepictionItem::Atom(AtomItem {
                    position: Point2D::new(2.0, -1.0),
                    label: "O-".to_owned(),
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
    #[case::heteroatom(r#"{:atoms ["N#h2"] :bonds []}"#, 0, Some("NH2"))]
    #[case::isolated_carbon(r#"{:atoms ["C#h4"] :bonds []}"#, 0, Some("C"))]
    #[case::isolated_charged_carbon(r#"{:atoms ["C#c+#h4"] :bonds []}"#, 0, Some("CH4+"))]
    #[case::skeleton_carbon(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#, 0, None)]
    #[case::skeleton_carbon_with_implicit_hydrogens(
        r#"{:atoms ["C#h3" "C"] :bonds [[0 1 "1"]]}"#,
        0,
        None
    )]
    #[case::isotopic_carbon(r#"{:atoms ["C#i13#h2" "C"] :bonds [[0 1 "1"]]}"#, 0, Some("13CH2"))]
    #[case::charged_carbon(r#"{:atoms ["C#c+" "C"] :bonds [[0 1 "1"]]}"#, 0, Some("C+"))]
    #[case::radical_carbon(r#"{:atoms ["C#u2" "C"] :bonds [[0 1 "1"]]}"#, 0, Some("C••"))]
    #[case::multi_digit_fields(
        r#"{:atoms ["N#i15#c-12#h12#u2"] :bonds []}"#,
        0,
        Some("15NH1212-••")
    )]
    #[case::independently_nonliteral_fields(
        r#"{:atoms ["N#i*#c*#h*#u*#s3"] :bonds []}"#,
        0,
        Some("N")
    )]
    fn test_atom_label(
        #[case] input: &str,
        #[case] atom_id: usize,
        #[case] expected: Option<&str>,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            atom_label(molecule.atom(AtomId::from(atom_id))).as_deref(),
            expected
        );
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

    #[test]
    fn test_depict_omits_nonliteral_projection_states() {
        let molecule = mol_dsl!(r#"{:atoms ["*#i13#c+#h2" "N#i*#c*#h*"] :bonds [[0 1 "*"]]}"#);
        let layout = layout(&[[0.0, 0.0], [1.0, 0.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            depiction.items(),
            [DepictionItem::Atom(AtomItem {
                position: Point2D::new(1.0, 0.0),
                label: "N".to_owned(),
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
            })]
        );
    }

    #[rstest]
    #[case::constraint_only(
        r#"{:atoms ["C#a+" "C#a+"] :bonds [[0 1 "1#a"]]}"#,
        vec![
            MarkerItem {
                position: Point2D::new(0.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(0)))],
            },
            MarkerItem {
                position: Point2D::new(2.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![DepictionReference::Molecule(Entity::Atom(AtomId(1)))],
            },
            MarkerItem {
                position: Point2D::new(1.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![DepictionReference::Molecule(Entity::Bond(BondId(0)))],
            },
        ]
    )]
    #[case::overlay_only(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "*"}]}"#,
        vec![
            MarkerItem {
                position: Point2D::new(0.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Atom(AtomId(0))),
                ],
            },
            MarkerItem {
                position: Point2D::new(2.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Atom(AtomId(1))),
                ],
            },
            MarkerItem {
                position: Point2D::new(1.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Bond(BondId(0))),
                ],
            },
        ]
    )]
    #[case::combined(
        r#"{:atoms ["C#a+" "C#a+"] :bonds [[0 1 "1#a"]] :aromatic-systems [{:atoms [0 1] :attrs "*"}]}"#,
        vec![
            MarkerItem {
                position: Point2D::new(0.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Atom(AtomId(0))),
                ],
            },
            MarkerItem {
                position: Point2D::new(2.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Atom(AtomId(1))),
                ],
            },
            MarkerItem {
                position: Point2D::new(1.0, 0.0),
                kind: MarkerKind::Aromatic,
                references: vec![
                    DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
                    DepictionReference::Molecule(Entity::Bond(BondId(0))),
                ],
            },
        ]
    )]
    #[case::nonaromatic(
        r#"{:atoms ["C#a!" "C#a!"] :bonds [[0 1 "1#a!"]]}"#,
        vec![]
    )]
    #[case::undetermined(
        r#"{:atoms ["C#a*" "C#a*"] :bonds [[0 1 "1#a*"]]}"#,
        vec![]
    )]
    fn test_depict_aromatic_projection(#[case] input: &str, #[case] expected: Vec<MarkerItem>) {
        let molecule = mol_dsl!(input);
        let layout = layout(&[[0.0, 0.0], [2.0, 0.0]]);

        assert_eq!(
            markers(&depict(&molecule, &layout).unwrap(), MarkerKind::Aromatic),
            expected
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

        assert!(depiction.items().iter().all(|item| !matches!(
            item,
            DepictionItem::Wedge(_)
                | DepictionItem::Marker(MarkerItem {
                    kind: MarkerKind::Stereo,
                    ..
                })
        )));
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
        assert!(markers(&depiction, MarkerKind::Stereo).is_empty());
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
        assert_eq!(depict(&baseline, &layout), depict(&changed, &layout));
    }

    #[rstest]
    fn test_depict_frame_size_mismatch() {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let layout = layout(&[]);

        assert_eq!(
            depict(&molecule, &layout),
            Err(MoleculeDepictionError::LayoutFrame(
                MoleculeLayoutError::FrameSizeMismatch {
                    molecule_atom_count: 1,
                    layout_atom_count: 0,
                }
            ))
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
            depict(&molecule, &layout),
            Err(MoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: StereoAtomId(0),
            })
        );
    }

    #[cfg(feature = "coordgen")]
    #[rstest]
    #[case::coordgen(MoleculeLayoutAlgorithm::CoordGen)]
    fn test_molecule_depict_with(#[case] algorithm: MoleculeLayoutAlgorithm) {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);
        let layout = layout_molecule(&molecule, algorithm).unwrap();
        let expected = depict(&molecule, &layout).unwrap();

        assert_eq!(molecule.depict_with(algorithm), Ok(expected));
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

    fn markers(depiction: &Depiction, kind: MarkerKind) -> Vec<MarkerItem> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Marker(marker) if marker.kind == kind => Some(marker.clone()),
                _ => None,
            })
            .collect()
    }
}
