//! Depiction construction from graph-IR molecules and supplied layouts.

use std::fmt::Write;

use umol_graph_ir::ir::{
    AsLit, AtomId, AtomView, Entity, IsotopeMass, Molecule, StereoAtomView, StereoBondView,
};

#[cfg(feature = "coordgen")]
use super::Depict;
use super::{
    AtomItem, BondItem, Depiction, DepictionItem, DepictionReference, MarkerItem, MarkerKind,
};
#[cfg(feature = "coordgen")]
use crate::layout::{layout_molecule, LayoutError, MoleculeLayoutAlgorithm};
use crate::layout::{MoleculeLayout, MoleculeLayoutError, Point2D};

/// Constructs the first format-neutral depiction projection of `molecule` in `layout`.
///
/// Localized bonds are followed by atom labels, aromatic markers, and stereo markers, with each
/// group ordered by graph-IR id. Nonliteral projected fields and unsupported overlays or
/// constraints are omitted. The first projection does not represent dative, multicenter, or
/// noncovalent bonds, unprojected inherent fields, or any constraint.
///
/// # Errors
///
/// Returns [`MoleculeLayoutError::FrameSizeMismatch`] if the molecule and layout do not use the
/// same dense atom frame.
pub fn depict(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<Depiction, MoleculeLayoutError> {
    layout.check_frame(molecule)?;

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
        items.push(DepictionItem::Bond(BondItem {
            start: position(layout, first),
            end: position(layout, second),
            line_count,
            references: vec![DepictionReference::Molecule(Entity::Bond(bond.id))],
        }));
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
        let system_reference = DepictionReference::Molecule(Entity::AromaticSystem(system.id));
        let mut atoms: Vec<_> = system.atoms().collect();
        atoms.sort_unstable_by_key(|atom| atom.id);
        for atom in atoms {
            items.push(DepictionItem::Marker(MarkerItem {
                position: position(layout, atom.id),
                kind: MarkerKind::Aromatic,
                references: vec![
                    system_reference,
                    DepictionReference::Molecule(Entity::Atom(atom.id)),
                ],
            }));
        }
        let mut bonds: Vec<_> = system.bonds().collect();
        bonds.sort_unstable_by_key(|bond| bond.id);
        for bond in bonds {
            let [first, second] = bond.atom_ids();
            items.push(DepictionItem::Marker(MarkerItem {
                position: midpoint(position(layout, first), position(layout, second)),
                kind: MarkerKind::Aromatic,
                references: vec![
                    system_reference,
                    DepictionReference::Molecule(Entity::Bond(bond.id)),
                ],
            }));
        }
    }

    items.extend(
        molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| stereo_atom_marker(layout, stereo)),
    );
    items.extend(
        molecule
            .stereo_bonds()
            .iter()
            .map(|stereo| stereo_bond_marker(layout, stereo)),
    );

    Ok(Depiction::from_items(items))
}

#[cfg(feature = "coordgen")]
impl Depict for Molecule {
    type Error = LayoutError;

    fn depict_with(
        &self,
        layout_algorithm: MoleculeLayoutAlgorithm,
    ) -> Result<Depiction, Self::Error> {
        let layout = layout_molecule(self, layout_algorithm)?;
        Ok(depict(self, &layout).expect("generated layout preserves the molecule atom frame"))
    }
}

fn atom_label(atom: AtomView<'_>) -> Option<String> {
    let element = atom.element().as_lit()?;
    let mut label = String::new();

    if let Some(IsotopeMass::MassNumber(mass)) = atom.isotope_mass().as_lit() {
        write!(label, "{mass}").expect("writing to a String cannot fail");
    }
    label.push_str(element.symbol());

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

    if let Some(charge) = atom.charge().as_lit().filter(|&charge| charge != 0) {
        let magnitude = charge.unsigned_abs();
        if magnitude != 1 {
            write!(label, "{magnitude}").expect("writing to a String cannot fail");
        }
        label.push(if charge > 0 { '+' } else { '-' });
    }

    Some(label)
}

fn stereo_atom_marker(layout: &MoleculeLayout, stereo: StereoAtomView<'_>) -> DepictionItem {
    let site = stereo.site_id();
    DepictionItem::Marker(MarkerItem {
        position: position(layout, site),
        kind: MarkerKind::Stereo,
        references: vec![
            DepictionReference::Molecule(Entity::StereoAtom(stereo.id)),
            DepictionReference::Molecule(Entity::Atom(site)),
        ],
    })
}

fn stereo_bond_marker(layout: &MoleculeLayout, stereo: StereoBondView<'_>) -> DepictionItem {
    let site = stereo.site();
    let [first, second] = site.atom_ids();
    DepictionItem::Marker(MarkerItem {
        position: midpoint(position(layout, first), position(layout, second)),
        kind: MarkerKind::Stereo,
        references: vec![
            DepictionReference::Molecule(Entity::StereoBond(stereo.id)),
            DepictionReference::Molecule(Entity::Bond(site.id)),
        ],
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, StereoAtomId, StereoBondId};
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

    #[test]
    fn test_depict_literal_atom_and_bond_projection() {
        let molecule = mol_dsl!(r#"{:atoms ["C#i13#c+#h2" "O#c-"] :bonds [[0 1 "2"]]}"#);
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
                    label: "13CH2+".to_owned(),
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
    #[case::element("C", "N", "C", "N")]
    #[case::isotope("C", "C#i13", "C", "13C")]
    #[case::charge("C", "C#c+", "C", "C+")]
    #[case::implicit_hydrogens("C", "C#h3", "C", "CH3")]
    fn test_depict_projected_atom_field_changes_output(
        #[case] before: &str,
        #[case] after: &str,
        #[case] before_label: &str,
        #[case] after_label: &str,
    ) {
        let before_input = format!(r#"{{:atoms ["{before}"] :bonds []}}"#);
        let after_input = format!(r#"{{:atoms ["{after}"] :bonds []}}"#);
        let before = mol_dsl!(before_input.as_str());
        let after = mol_dsl!(after_input.as_str());
        let layout = layout(&[[0.0, 0.0]]);

        assert_eq!(
            atom_labels(&depict(&before, &layout).unwrap()),
            [before_label]
        );
        assert_eq!(
            atom_labels(&depict(&after, &layout).unwrap()),
            [after_label]
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

    #[test]
    fn test_depict_aromatic_system_marks_members_and_induced_bonds() {
        let molecule = mol_dsl!(
            r##"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#e3"}]}"##
        );
        let layout = layout(&[[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]]);

        let depiction = depict(&molecule, &layout).unwrap();

        assert_eq!(
            markers(&depiction, MarkerKind::Aromatic),
            [
                MarkerItem {
                    position: Point2D::new(0.0, 0.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Atom(AtomId(0))),
                },
                MarkerItem {
                    position: Point2D::new(2.0, 0.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Atom(AtomId(1))),
                },
                MarkerItem {
                    position: Point2D::new(1.0, 2.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Atom(AtomId(2))),
                },
                MarkerItem {
                    position: Point2D::new(1.0, 0.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Bond(BondId(0))),
                },
                MarkerItem {
                    position: Point2D::new(1.5, 1.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Bond(BondId(1))),
                },
                MarkerItem {
                    position: Point2D::new(0.5, 1.0),
                    kind: MarkerKind::Aromatic,
                    references: aromatic_references(Entity::Bond(BondId(2))),
                },
            ]
        );
    }

    #[rstest]
    #[case::atom(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]],
        MarkerItem {
            position: Point2D::new(0.0, 0.0),
            kind: MarkerKind::Stereo,
            references: vec![
                DepictionReference::Molecule(Entity::StereoAtom(StereoAtomId(0))),
                DepictionReference::Molecule(Entity::Atom(AtomId(0))),
            ],
        }
    )]
    #[case::bond(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#,
        vec![[0.0, 0.0], [1.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
        MarkerItem {
            position: Point2D::new(2.0, 0.0),
            kind: MarkerKind::Stereo,
            references: vec![
                DepictionReference::Molecule(Entity::StereoBond(StereoBondId(0))),
                DepictionReference::Molecule(Entity::Bond(BondId(1))),
            ],
        }
    )]
    fn test_depict_stereo_entity_marks_site(
        #[case] input: &str,
        #[case] positions: Vec<[f64; 2]>,
        #[case] expected: MarkerItem,
    ) {
        let molecule = mol_dsl!(input);
        let layout = layout(&positions);

        assert_eq!(
            markers(&depict(&molecule, &layout).unwrap(), MarkerKind::Stereo),
            [expected]
        );
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

    #[test]
    fn test_depict_frame_size_mismatch() {
        let molecule = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let layout = layout(&[]);

        assert_eq!(
            depict(&molecule, &layout),
            Err(MoleculeLayoutError::FrameSizeMismatch {
                molecule_atom_count: 1,
                layout_atom_count: 0,
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

    fn atom_labels(depiction: &Depiction) -> Vec<&str> {
        depiction
            .items()
            .iter()
            .filter_map(|item| match item {
                DepictionItem::Atom(atom) => Some(atom.label.as_str()),
                _ => None,
            })
            .collect()
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

    fn aromatic_references(anchor: Entity) -> Vec<DepictionReference> {
        vec![
            DepictionReference::Molecule(Entity::AromaticSystem(AromaticSystemId(0))),
            DepictionReference::Molecule(anchor),
        ]
    }
}
