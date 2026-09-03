use rstest::rstest;
use umol_coordgen_sys::{generate_coordinates, Bond, CisTransBond, CoordgenError, SideRelation};

#[rstest]
#[case::empty_result(Vec::new(), Vec::new())]
fn test_generate_coordinates_empty_graph(
    #[case] atomic_numbers: Vec<u16>,
    #[case] bonds: Vec<Bond>,
) {
    let points = generate_coordinates(&atomic_numbers, &bonds, &[]).expect("coordinate generation");

    assert_eq!(points, Vec::new());
}

#[rstest]
#[case::finite_point(vec![6], Vec::new())]
fn test_generate_coordinates_one_atom(#[case] atomic_numbers: Vec<u16>, #[case] bonds: Vec<Bond>) {
    let points = generate_coordinates(&atomic_numbers, &bonds, &[]).expect("coordinate generation");

    assert_eq!(points.len(), 1);
    assert!(points[0].x.is_finite());
    assert!(points[0].y.is_finite());
}

#[rstest]
#[case::preserved_frame(
    vec![8, 7, 9, 6],
    vec![
        Bond {
            atom_0: 3,
            atom_1: 0,
            order: 1,
        },
        Bond {
            atom_0: 3,
            atom_1: 1,
            order: 2,
        },
        Bond {
            atom_0: 3,
            atom_1: 2,
            order: 1,
        },
    ]
)]
fn test_generate_coordinates_bonded_graph(
    #[case] atomic_numbers: Vec<u16>,
    #[case] bonds: Vec<Bond>,
) {
    let points = generate_coordinates(&atomic_numbers, &bonds, &[]).expect("coordinate generation");

    assert_eq!(points.len(), 4);
    assert!(points
        .iter()
        .all(|point| point.x.is_finite() && point.y.is_finite()));
    for terminal in 0..3 {
        let distance = (points[3].x - points[terminal].x).hypot(points[3].y - points[terminal].y);
        assert!((distance - 50.0).abs() < 0.1, "bond length was {distance}");
    }
    let terminal_distance = (points[0].x - points[1].x).hypot(points[0].y - points[1].y);
    assert!(terminal_distance > 50.0);
}

#[rstest]
#[case::same_side(SideRelation::SameSide)]
#[case::opposite_side(SideRelation::OppositeSide)]
fn test_generate_coordinates_cis_trans(#[case] relation: SideRelation) {
    let atomic_numbers = [6, 6, 9, 17];
    let bonds = [
        Bond {
            atom_0: 0,
            atom_1: 1,
            order: 2,
        },
        Bond {
            atom_0: 0,
            atom_1: 2,
            order: 1,
        },
        Bond {
            atom_0: 1,
            atom_1: 3,
            order: 1,
        },
    ];
    let cis_trans_bonds = [CisTransBond {
        bond: 0,
        first_ligand: 2,
        second_ligand: 3,
        relation,
    }];

    let points = generate_coordinates(&atomic_numbers, &bonds, &cis_trans_bonds)
        .expect("coordinate generation");

    assert_eq!(
        observed_relation(&points, &bonds[0], &cis_trans_bonds[0]),
        relation
    );
}

#[rstest]
fn test_generate_coordinates_cis_trans_endpoint_orientation() {
    let atomic_numbers = [6, 6, 9, 17];
    let bonds = [
        Bond {
            atom_0: 0,
            atom_1: 1,
            order: 2,
        },
        Bond {
            atom_0: 0,
            atom_1: 2,
            order: 1,
        },
        Bond {
            atom_0: 1,
            atom_1: 3,
            order: 1,
        },
    ];
    let reversed_bonds = [
        Bond {
            atom_0: 1,
            atom_1: 0,
            order: 2,
        },
        bonds[1],
        bonds[2],
    ];
    let cis_trans_bonds = [CisTransBond {
        bond: 0,
        first_ligand: 2,
        second_ligand: 3,
        relation: SideRelation::OppositeSide,
    }];
    let reversed_cis_trans_bonds = [CisTransBond {
        bond: 0,
        first_ligand: 3,
        second_ligand: 2,
        relation: SideRelation::OppositeSide,
    }];

    let points = generate_coordinates(&atomic_numbers, &bonds, &cis_trans_bonds)
        .expect("first coordinate generation");
    let repeated = generate_coordinates(&atomic_numbers, &bonds, &cis_trans_bonds)
        .expect("repeated coordinate generation");
    let reversed =
        generate_coordinates(&atomic_numbers, &reversed_bonds, &reversed_cis_trans_bonds)
            .expect("reversed coordinate generation");

    assert_eq!(points, repeated);
    assert_layouts_congruent(&points, &reversed);
    assert_eq!(
        observed_relation(&reversed, &reversed_bonds[0], &reversed_cis_trans_bonds[0]),
        SideRelation::OppositeSide
    );
}

#[rstest]
#[case::first_endpoint(
    Bond {
        atom_0: 2,
        atom_1: 0,
        order: 1,
    },
    2
)]
#[case::second_endpoint(
    Bond {
        atom_0: 0,
        atom_1: 2,
        order: 1,
    },
    2
)]
fn test_generate_coordinates_error(#[case] bond: Bond, #[case] atom_index: usize) {
    assert_eq!(
        generate_coordinates(&[6], &[bond], &[]),
        Err(CoordgenError::BondAtomOutOfBounds {
            bond_index: 0,
            atom_index,
            atom_count: 1,
        })
    );
}

#[derive(Clone, Copy, Debug)]
enum CisTransInputErrorCase {
    SiteOutOfBounds,
    SiteOrder,
    FirstLigandOutOfBounds,
    SecondLigandOutOfBounds,
    FirstLigandIsSiteAtom,
    SecondLigandIsSiteAtom,
    FirstLigandNotIncident,
    SecondLigandNotIncident,
    DuplicateSite,
}

#[rstest]
#[case::site_out_of_bounds(
    CisTransInputErrorCase::SiteOutOfBounds,
    CoordgenError::CisTransSiteOutOfBounds {
        cis_trans_index: 0,
        bond_index: 3,
        bond_count: 3,
    }
)]
#[case::site_order(
    CisTransInputErrorCase::SiteOrder,
    CoordgenError::CisTransSiteOrder {
        cis_trans_index: 0,
        bond_index: 0,
        order: 1,
    }
)]
#[case::first_ligand_out_of_bounds(
    CisTransInputErrorCase::FirstLigandOutOfBounds,
    CoordgenError::CisTransLigandOutOfBounds {
        cis_trans_index: 0,
        ligand_position: 0,
        atom_index: 4,
        atom_count: 4,
    }
)]
#[case::second_ligand_out_of_bounds(
    CisTransInputErrorCase::SecondLigandOutOfBounds,
    CoordgenError::CisTransLigandOutOfBounds {
        cis_trans_index: 0,
        ligand_position: 1,
        atom_index: 4,
        atom_count: 4,
    }
)]
#[case::first_ligand_is_site_atom(
    CisTransInputErrorCase::FirstLigandIsSiteAtom,
    CoordgenError::CisTransLigandIsSiteAtom {
        cis_trans_index: 0,
        ligand_position: 0,
        atom_index: 0,
    }
)]
#[case::second_ligand_is_site_atom(
    CisTransInputErrorCase::SecondLigandIsSiteAtom,
    CoordgenError::CisTransLigandIsSiteAtom {
        cis_trans_index: 0,
        ligand_position: 1,
        atom_index: 1,
    }
)]
#[case::first_ligand_not_incident(
    CisTransInputErrorCase::FirstLigandNotIncident,
    CoordgenError::CisTransLigandNotIncident {
        cis_trans_index: 0,
        ligand_position: 0,
        ligand_atom: 3,
        endpoint_atom: 0,
    }
)]
#[case::second_ligand_not_incident(
    CisTransInputErrorCase::SecondLigandNotIncident,
    CoordgenError::CisTransLigandNotIncident {
        cis_trans_index: 0,
        ligand_position: 1,
        ligand_atom: 2,
        endpoint_atom: 1,
    }
)]
#[case::duplicate_site(
    CisTransInputErrorCase::DuplicateSite,
    CoordgenError::DuplicateCisTransSite {
        first_cis_trans_index: 0,
        second_cis_trans_index: 1,
        bond_index: 0,
    }
)]
fn test_generate_coordinates_cis_trans_input_error(
    #[case] scenario: CisTransInputErrorCase,
    #[case] expected: CoordgenError,
) {
    let atomic_numbers = [6, 6, 6, 6];
    let mut bonds = vec![
        Bond {
            atom_0: 0,
            atom_1: 1,
            order: 2,
        },
        Bond {
            atom_0: 0,
            atom_1: 2,
            order: 1,
        },
        Bond {
            atom_0: 1,
            atom_1: 3,
            order: 1,
        },
    ];
    let mut cis_trans_bonds = vec![CisTransBond {
        bond: 0,
        first_ligand: 2,
        second_ligand: 3,
        relation: SideRelation::SameSide,
    }];

    match scenario {
        CisTransInputErrorCase::SiteOutOfBounds => cis_trans_bonds[0].bond = 3,
        CisTransInputErrorCase::SiteOrder => bonds[0].order = 1,
        CisTransInputErrorCase::FirstLigandOutOfBounds => {
            cis_trans_bonds[0].first_ligand = 4;
        }
        CisTransInputErrorCase::SecondLigandOutOfBounds => {
            cis_trans_bonds[0].second_ligand = 4;
        }
        CisTransInputErrorCase::FirstLigandIsSiteAtom => {
            cis_trans_bonds[0].first_ligand = 0;
        }
        CisTransInputErrorCase::SecondLigandIsSiteAtom => {
            cis_trans_bonds[0].second_ligand = 1;
        }
        CisTransInputErrorCase::FirstLigandNotIncident => {
            cis_trans_bonds[0].first_ligand = 3;
        }
        CisTransInputErrorCase::SecondLigandNotIncident => {
            cis_trans_bonds[0].second_ligand = 2;
        }
        CisTransInputErrorCase::DuplicateSite => cis_trans_bonds.push(cis_trans_bonds[0]),
    }

    assert_eq!(
        generate_coordinates(&atomic_numbers, &bonds, &cis_trans_bonds),
        Err(expected)
    );
}

fn observed_relation(
    points: &[umol_coordgen_sys::Point],
    site: &Bond,
    cis_trans: &CisTransBond,
) -> SideRelation {
    let first = signed_half_plane(
        points[site.atom_0],
        points[site.atom_1],
        points[cis_trans.first_ligand],
    );
    let second = signed_half_plane(
        points[site.atom_0],
        points[site.atom_1],
        points[cis_trans.second_ligand],
    );
    assert!(first.abs() > 1e-6);
    assert!(second.abs() > 1e-6);
    if first.is_sign_positive() == second.is_sign_positive() {
        SideRelation::SameSide
    } else {
        SideRelation::OppositeSide
    }
}

fn signed_half_plane(
    site_0: umol_coordgen_sys::Point,
    site_1: umol_coordgen_sys::Point,
    ligand: umol_coordgen_sys::Point,
) -> f64 {
    (site_1.x - site_0.x) * (ligand.y - site_0.y) - (site_1.y - site_0.y) * (ligand.x - site_0.x)
}

fn assert_layouts_congruent(
    first: &[umol_coordgen_sys::Point],
    second: &[umol_coordgen_sys::Point],
) {
    assert_eq!(first.len(), second.len());
    for atom_0 in 0..first.len() {
        for atom_1 in 0..first.len() {
            let first_distance =
                (first[atom_1].x - first[atom_0].x).hypot(first[atom_1].y - first[atom_0].y);
            let second_distance =
                (second[atom_1].x - second[atom_0].x).hypot(second[atom_1].y - second[atom_0].y);
            assert!(
                (first_distance - second_distance).abs() < 1e-3,
                "distance between atoms {atom_0} and {atom_1} changed from {first_distance} to {second_distance}"
            );
        }
    }
}
