//! Safe data boundary for umol's vendored CoordGen integration.
//!
//! Enabling the `native` feature compiles the vendored C++ source and exposes
//! [`generate_coordinates`]. The default feature set does not invoke a C++ compiler.

/// Vendored upstream CoordGen version.
pub const COORDGEN_VERSION: &str = "3.0.2";

#[cfg(feature = "native")]
mod native {
    use std::collections::HashMap;
    use std::error::Error;
    use std::fmt::{self, Display, Formatter};

    const OK: i32 = 0;
    const NULL_POINTER: i32 = 1;
    const ATOM_OUT_OF_BOUNDS: i32 = 2;
    const ALLOCATION_FAILED: i32 = 3;
    const BACKEND_EXCEPTION: i32 = 4;
    const CIS_TRANS_SITE_OUT_OF_BOUNDS: i32 = 5;
    const CIS_TRANS_LIGAND_OUT_OF_BOUNDS: i32 = 6;
    const INVALID_SIDE_RELATION: i32 = 7;
    const RELATIVE_SIDE_TOLERANCE: f64 = 1e-6;

    unsafe extern "C" {
        fn umol_coordgen_generate(
            atom_count: usize,
            atomic_numbers: *const u16,
            bond_count: usize,
            bonds: *const Bond,
            cis_trans_bond_count: usize,
            cis_trans_bonds: *const CisTransBond,
            points: *mut Point,
        ) -> i32;
    }

    /// Bond input for native coordinate generation.
    ///
    /// Atom positions are zero-based indices into the atomic-number slice supplied to
    /// [`generate_coordinates`].
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct Bond {
        pub atom_0: usize,
        pub atom_1: usize,
        pub order: u8,
    }

    /// Relative placement requested for two selected double-bond ligands.
    #[repr(u8)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    pub enum SideRelation {
        /// The selected ligands lie on the same side of the site-bond axis.
        SameSide = 0,
        /// The selected ligands lie on opposite sides of the site-bond axis.
        OppositeSide = 1,
    }

    /// Cis/trans input for one bond in a coordinate-generation request.
    ///
    /// Indices are interpreted in the atom and bond frames supplied to [`generate_coordinates`].
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct CisTransBond {
        /// Index of the double bond in the supplied bond slice.
        pub bond: usize,
        /// Substituent of the site's `atom_0` endpoint.
        pub first_ligand: usize,
        /// Substituent of the site's `atom_1` endpoint.
        pub second_ligand: usize,
        /// Requested relative placement of the selected substituents.
        pub relation: SideRelation,
    }

    /// Two-dimensional coordinate returned by CoordGen.
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Point {
        pub x: f64,
        pub y: f64,
    }

    /// Failure while validating coordinate-generation input or running CoordGen.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub enum CoordgenError {
        BondAtomOutOfBounds {
            bond_index: usize,
            atom_index: usize,
            atom_count: usize,
        },
        CisTransSiteOutOfBounds {
            cis_trans_index: usize,
            bond_index: usize,
            bond_count: usize,
        },
        CisTransSiteOrder {
            cis_trans_index: usize,
            bond_index: usize,
            order: u8,
        },
        CisTransLigandOutOfBounds {
            cis_trans_index: usize,
            ligand_position: usize,
            atom_index: usize,
            atom_count: usize,
        },
        CisTransLigandIsSiteAtom {
            cis_trans_index: usize,
            ligand_position: usize,
            atom_index: usize,
        },
        CisTransLigandNotIncident {
            cis_trans_index: usize,
            ligand_position: usize,
            ligand_atom: usize,
            endpoint_atom: usize,
        },
        DuplicateCisTransSite {
            first_cis_trans_index: usize,
            second_cis_trans_index: usize,
            bond_index: usize,
        },
        DegenerateCisTransGeometry {
            cis_trans_index: usize,
            bond_index: usize,
        },
        CisTransGeometryMismatch {
            cis_trans_index: usize,
            bond_index: usize,
            relation: SideRelation,
        },
        AllocationFailed,
        BackendException,
        NonFinitePoint {
            atom_index: usize,
        },
        UnexpectedStatus {
            status: i32,
        },
    }

    impl Display for CoordgenError {
        fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
            match self {
                Self::BondAtomOutOfBounds {
                    bond_index,
                    atom_index,
                    atom_count,
                } => write!(
                    formatter,
                    "bond {bond_index} references atom {atom_index} outside frame of size {atom_count}"
                ),
                Self::CisTransSiteOutOfBounds {
                    cis_trans_index,
                    bond_index,
                    bond_count,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} references bond {bond_index} outside frame of size {bond_count}"
                ),
                Self::CisTransSiteOrder {
                    cis_trans_index,
                    bond_index,
                    order,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} references bond {bond_index} with order {order}, expected 2"
                ),
                Self::CisTransLigandOutOfBounds {
                    cis_trans_index,
                    ligand_position,
                    atom_index,
                    atom_count,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} ligand {ligand_position} references atom {atom_index} outside frame of size {atom_count}"
                ),
                Self::CisTransLigandIsSiteAtom {
                    cis_trans_index,
                    ligand_position,
                    atom_index,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} ligand {ligand_position} is site atom {atom_index}"
                ),
                Self::CisTransLigandNotIncident {
                    cis_trans_index,
                    ligand_position,
                    ligand_atom,
                    endpoint_atom,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} ligand {ligand_position} atom {ligand_atom} is not bonded to site endpoint {endpoint_atom}"
                ),
                Self::DuplicateCisTransSite {
                    first_cis_trans_index,
                    second_cis_trans_index,
                    bond_index,
                } => write!(
                    formatter,
                    "cis/trans inputs {first_cis_trans_index} and {second_cis_trans_index} both reference site bond {bond_index}"
                ),
                Self::DegenerateCisTransGeometry {
                    cis_trans_index,
                    bond_index,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} on bond {bond_index} produced degenerate geometry"
                ),
                Self::CisTransGeometryMismatch {
                    cis_trans_index,
                    bond_index,
                    relation,
                } => write!(
                    formatter,
                    "cis/trans input {cis_trans_index} on bond {bond_index} did not produce {relation:?} geometry"
                ),
                Self::AllocationFailed => formatter.write_str("CoordGen allocation failed"),
                Self::BackendException => formatter.write_str("CoordGen raised an exception"),
                Self::NonFinitePoint { atom_index } => {
                    write!(formatter, "CoordGen returned a non-finite point for atom {atom_index}")
                }
                Self::UnexpectedStatus { status } => {
                    write!(formatter, "CoordGen returned unexpected status {status}")
                }
            }
        }
    }

    impl Error for CoordgenError {}

    /// Generates one point per input atom while preserving the input atom order.
    ///
    /// Atomic number zero denotes a generic layout atom. Empty input is valid. The backend's
    /// native bond-length scale is preserved; the graph-IR adapter owns normalization.
    ///
    /// # Errors
    ///
    /// Returns [`CoordgenError::BondAtomOutOfBounds`] for a bond endpoint outside the supplied
    /// atom frame. Native allocation failures, exceptions, and non-finite output are reported by
    /// their corresponding variants.
    pub fn generate_coordinates(
        atomic_numbers: &[u16],
        bonds: &[Bond],
        cis_trans_bonds: &[CisTransBond],
    ) -> Result<Vec<Point>, CoordgenError> {
        for (bond_index, bond) in bonds.iter().enumerate() {
            for atom_index in [bond.atom_0, bond.atom_1] {
                if atom_index >= atomic_numbers.len() {
                    return Err(CoordgenError::BondAtomOutOfBounds {
                        bond_index,
                        atom_index,
                        atom_count: atomic_numbers.len(),
                    });
                }
            }
        }
        validate_cis_trans_bonds(atomic_numbers.len(), bonds, cis_trans_bonds)?;

        let mut points = vec![Point { x: 0.0, y: 0.0 }; atomic_numbers.len()];
        let status = unsafe {
            umol_coordgen_generate(
                atomic_numbers.len(),
                atomic_numbers.as_ptr(),
                bonds.len(),
                bonds.as_ptr(),
                cis_trans_bonds.len(),
                cis_trans_bonds.as_ptr(),
                points.as_mut_ptr(),
            )
        };
        match status {
            OK => {}
            NULL_POINTER
            | ATOM_OUT_OF_BOUNDS
            | CIS_TRANS_SITE_OUT_OF_BOUNDS
            | CIS_TRANS_LIGAND_OUT_OF_BOUNDS
            | INVALID_SIDE_RELATION => {
                return Err(CoordgenError::UnexpectedStatus { status });
            }
            ALLOCATION_FAILED => return Err(CoordgenError::AllocationFailed),
            BACKEND_EXCEPTION => return Err(CoordgenError::BackendException),
            _ => return Err(CoordgenError::UnexpectedStatus { status }),
        }

        if let Some(atom_index) = points
            .iter()
            .position(|point| !point.x.is_finite() || !point.y.is_finite())
        {
            return Err(CoordgenError::NonFinitePoint { atom_index });
        }
        validate_cis_trans_geometry(&points, bonds, cis_trans_bonds)?;
        Ok(points)
    }

    fn validate_cis_trans_bonds(
        atom_count: usize,
        bonds: &[Bond],
        cis_trans_bonds: &[CisTransBond],
    ) -> Result<(), CoordgenError> {
        let mut sites = HashMap::new();
        for (cis_trans_index, cis_trans) in cis_trans_bonds.iter().enumerate() {
            let Some(site) = bonds.get(cis_trans.bond) else {
                return Err(CoordgenError::CisTransSiteOutOfBounds {
                    cis_trans_index,
                    bond_index: cis_trans.bond,
                    bond_count: bonds.len(),
                });
            };
            if site.order != 2 {
                return Err(CoordgenError::CisTransSiteOrder {
                    cis_trans_index,
                    bond_index: cis_trans.bond,
                    order: site.order,
                });
            }
            if let Some(first_cis_trans_index) = sites.insert(cis_trans.bond, cis_trans_index) {
                return Err(CoordgenError::DuplicateCisTransSite {
                    first_cis_trans_index,
                    second_cis_trans_index: cis_trans_index,
                    bond_index: cis_trans.bond,
                });
            }

            for (ligand_position, ligand_atom, endpoint_atom) in [
                (0, cis_trans.first_ligand, site.atom_0),
                (1, cis_trans.second_ligand, site.atom_1),
            ] {
                if ligand_atom >= atom_count {
                    return Err(CoordgenError::CisTransLigandOutOfBounds {
                        cis_trans_index,
                        ligand_position,
                        atom_index: ligand_atom,
                        atom_count,
                    });
                }
                if ligand_atom == site.atom_0 || ligand_atom == site.atom_1 {
                    return Err(CoordgenError::CisTransLigandIsSiteAtom {
                        cis_trans_index,
                        ligand_position,
                        atom_index: ligand_atom,
                    });
                }
                if !bonds.iter().any(|bond| {
                    [bond.atom_0, bond.atom_1] == [endpoint_atom, ligand_atom]
                        || [bond.atom_0, bond.atom_1] == [ligand_atom, endpoint_atom]
                }) {
                    return Err(CoordgenError::CisTransLigandNotIncident {
                        cis_trans_index,
                        ligand_position,
                        ligand_atom,
                        endpoint_atom,
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_cis_trans_geometry(
        points: &[Point],
        bonds: &[Bond],
        cis_trans_bonds: &[CisTransBond],
    ) -> Result<(), CoordgenError> {
        for (cis_trans_index, cis_trans) in cis_trans_bonds.iter().enumerate() {
            let site = bonds[cis_trans.bond];
            let Some(actual) = relative_side(
                points[site.atom_0],
                points[site.atom_1],
                points[cis_trans.first_ligand],
                points[cis_trans.second_ligand],
            ) else {
                return Err(CoordgenError::DegenerateCisTransGeometry {
                    cis_trans_index,
                    bond_index: cis_trans.bond,
                });
            };
            if actual != cis_trans.relation {
                return Err(CoordgenError::CisTransGeometryMismatch {
                    cis_trans_index,
                    bond_index: cis_trans.bond,
                    relation: cis_trans.relation,
                });
            }
        }
        Ok(())
    }

    fn relative_side(
        site_0: Point,
        site_1: Point,
        first_ligand: Point,
        second_ligand: Point,
    ) -> Option<SideRelation> {
        let axis_x = site_1.x - site_0.x;
        let axis_y = site_1.y - site_0.y;
        let axis_length = axis_x.hypot(axis_y);
        if axis_length == 0.0 || !axis_length.is_finite() {
            return None;
        }

        let first_x = first_ligand.x - site_0.x;
        let first_y = first_ligand.y - site_0.y;
        let second_x = second_ligand.x - site_1.x;
        let second_y = second_ligand.y - site_1.y;
        let first_side = half_plane(axis_x, axis_y, axis_length, first_x, first_y)?;
        let second_side = half_plane(axis_x, axis_y, axis_length, second_x, second_y)?;

        Some(if first_side == second_side {
            SideRelation::SameSide
        } else {
            SideRelation::OppositeSide
        })
    }

    fn half_plane(
        axis_x: f64,
        axis_y: f64,
        axis_length: f64,
        ligand_x: f64,
        ligand_y: f64,
    ) -> Option<bool> {
        let ligand_length = ligand_x.hypot(ligand_y);
        if ligand_length == 0.0 || !ligand_length.is_finite() {
            return None;
        }
        let cross = axis_x * ligand_y - axis_y * ligand_x;
        let tolerance = RELATIVE_SIDE_TOLERANCE * axis_length * ligand_length;
        (cross.abs() > tolerance).then_some(cross.is_sign_positive())
    }

    #[cfg(test)]
    mod tests {
        use rstest::rstest;

        use super::*;

        #[repr(C)]
        struct RawCisTransBond {
            bond: usize,
            first_ligand: usize,
            second_ligand: usize,
            relation: u8,
        }

        #[rstest]
        fn test_native_side_relation_discriminator() {
            let atomic_numbers = [6, 6, 6, 6];
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
            let cis_trans_bond = RawCisTransBond {
                bond: 0,
                first_ligand: 2,
                second_ligand: 3,
                relation: 2,
            };
            let mut points = [Point { x: 0.0, y: 0.0 }; 4];

            let status = unsafe {
                umol_coordgen_generate(
                    atomic_numbers.len(),
                    atomic_numbers.as_ptr(),
                    bonds.len(),
                    bonds.as_ptr(),
                    1,
                    (&raw const cis_trans_bond).cast(),
                    points.as_mut_ptr(),
                )
            };

            assert_eq!(status, INVALID_SIDE_RELATION);
        }

        #[rstest]
        #[case::same_side(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 1.0, y: 0.0 },
                Point { x: 0.0, y: 1.0 },
                Point { x: 1.0, y: 1.0 },
            ],
            SideRelation::SameSide,
            Ok(())
        )]
        #[case::opposite_side(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 1.0, y: 0.0 },
                Point { x: 0.0, y: 1.0 },
                Point { x: 1.0, y: -1.0 },
            ],
            SideRelation::OppositeSide,
            Ok(())
        )]
        #[case::mismatch(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 1.0, y: 0.0 },
                Point { x: 0.0, y: 1.0 },
                Point { x: 1.0, y: 1.0 },
            ],
            SideRelation::OppositeSide,
            Err(CoordgenError::CisTransGeometryMismatch {
                cis_trans_index: 0,
                bond_index: 0,
                relation: SideRelation::OppositeSide,
            })
        )]
        #[case::zero_length_site(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 0.0, y: 0.0 },
                Point { x: 0.0, y: 1.0 },
                Point { x: 1.0, y: 1.0 },
            ],
            SideRelation::SameSide,
            Err(CoordgenError::DegenerateCisTransGeometry {
                cis_trans_index: 0,
                bond_index: 0,
            })
        )]
        #[case::first_collinear_ligand(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 1.0, y: 0.0 },
                Point { x: -1.0, y: 0.0 },
                Point { x: 1.0, y: 1.0 },
            ],
            SideRelation::SameSide,
            Err(CoordgenError::DegenerateCisTransGeometry {
                cis_trans_index: 0,
                bond_index: 0,
            })
        )]
        #[case::second_collinear_ligand(
            [
                Point { x: 0.0, y: 0.0 },
                Point { x: 1.0, y: 0.0 },
                Point { x: 0.0, y: 1.0 },
                Point { x: 2.0, y: 0.0 },
            ],
            SideRelation::SameSide,
            Err(CoordgenError::DegenerateCisTransGeometry {
                cis_trans_index: 0,
                bond_index: 0,
            })
        )]
        fn test_validate_cis_trans_geometry(
            #[case] points: [Point; 4],
            #[case] relation: SideRelation,
            #[case] expected: Result<(), CoordgenError>,
        ) {
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

            assert_eq!(
                validate_cis_trans_geometry(&points, &bonds, &cis_trans_bonds),
                expected
            );
        }
    }
}

#[cfg(feature = "native")]
pub use native::{generate_coordinates, Bond, CisTransBond, CoordgenError, Point, SideRelation};
