//! Raw FFI bindings to libmsym.

#![allow(non_camel_case_types, non_upper_case_globals)]

use std::os::raw::{c_char, c_double, c_int, c_void};

// ---------------------------------------------------------------------------
// Opaque context
// ---------------------------------------------------------------------------

pub type msym_context = *mut c_void;

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

pub type msym_error_t = c_int;

pub const MSYM_SUCCESS: msym_error_t = 0;
pub const MSYM_INVALID_INPUT: msym_error_t = -1;
pub const MSYM_INVALID_CONTEXT: msym_error_t = -2;
pub const MSYM_INVALID_THRESHOLD: msym_error_t = -3;
pub const MSYM_INVALID_ELEMENTS: msym_error_t = -4;
pub const MSYM_INVALID_BASIS_FUNCTIONS: msym_error_t = -5;
pub const MSYM_INVALID_POINT_GROUP: msym_error_t = -6;
pub const MSYM_INVALID_EQUIVALENCE_SET: msym_error_t = -7;
pub const MSYM_INVALID_PERMUTATION: msym_error_t = -8;
pub const MSYM_INVALID_GEOMETRY: msym_error_t = -9;
pub const MSYM_INVALID_CHARACTER_TABLE: msym_error_t = -10;
pub const MSYM_INVALID_SUBSPACE: msym_error_t = -11;
pub const MSYM_INVALID_SUBGROUPS: msym_error_t = -12;
pub const MSYM_INVALID_AXES: msym_error_t = -13;
pub const MSYM_SYMMETRY_ERROR: msym_error_t = -14;
pub const MSYM_PERMUTATION_ERROR: msym_error_t = -15;
pub const MSYM_POINT_GROUP_ERROR: msym_error_t = -16;
pub const MSYM_SYMMETRIZATION_ERROR: msym_error_t = -17;
pub const MSYM_SUBSPACE_ERROR: msym_error_t = -18;
pub const MSYM_MEMORY_ERROR: msym_error_t = -128;

pub type msym_geometry_t = c_int;

pub const MSYM_GEOMETRY_UNKNOWN: msym_geometry_t = 0;
pub const MSYM_GEOMETRY_SPHERICAL: msym_geometry_t = 1;
pub const MSYM_GEOMETRY_LINEAR: msym_geometry_t = 2;
pub const MSYM_GEOMETRY_PLANAR_REGULAR: msym_geometry_t = 3;
pub const MSYM_GEOMETRY_PLANAR_IRREGULAR: msym_geometry_t = 4;
pub const MSYM_GEOMETRY_POLYHEDRAL_PROLATE: msym_geometry_t = 5;
pub const MSYM_GEOMETRY_POLYHEDRAL_OBLATE: msym_geometry_t = 6;
pub const MSYM_GEOMETRY_ASSYMETRIC: msym_geometry_t = 7;

pub type msym_symmetry_operation_type_t = c_int;

pub const MSYM_SYMMETRY_OPERATION_TYPE_IDENTITY: msym_symmetry_operation_type_t = 0;
pub const MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION: msym_symmetry_operation_type_t = 1;
pub const MSYM_SYMMETRY_OPERATION_TYPE_IMPROPER_ROTATION: msym_symmetry_operation_type_t = 2;
pub const MSYM_SYMMETRY_OPERATION_TYPE_REFLECTION: msym_symmetry_operation_type_t = 3;
pub const MSYM_SYMMETRY_OPERATION_TYPE_INVERSION: msym_symmetry_operation_type_t = 4;

pub type msym_symmetry_operation_orientation_t = c_int;

pub const MSYM_SYMMETRY_OPERATION_ORIENTATION_NONE: msym_symmetry_operation_orientation_t = 0;
pub const MSYM_SYMMETRY_OPERATION_ORIENTATION_HORIZONTAL: msym_symmetry_operation_orientation_t = 1;
pub const MSYM_SYMMETRY_OPERATION_ORIENTATION_VERTICAL: msym_symmetry_operation_orientation_t = 2;
pub const MSYM_SYMMETRY_OPERATION_ORIENTATION_DIHEDRAL: msym_symmetry_operation_orientation_t = 3;

pub type msym_point_group_type_t = c_int;

pub const MSYM_POINT_GROUP_TYPE_Kh: msym_point_group_type_t = 0;
pub const MSYM_POINT_GROUP_TYPE_K: msym_point_group_type_t = 1;
pub const MSYM_POINT_GROUP_TYPE_Ci: msym_point_group_type_t = 2;
pub const MSYM_POINT_GROUP_TYPE_Cs: msym_point_group_type_t = 3;
pub const MSYM_POINT_GROUP_TYPE_Cn: msym_point_group_type_t = 4;
pub const MSYM_POINT_GROUP_TYPE_Cnh: msym_point_group_type_t = 5;
pub const MSYM_POINT_GROUP_TYPE_Cnv: msym_point_group_type_t = 6;
pub const MSYM_POINT_GROUP_TYPE_Dn: msym_point_group_type_t = 7;
pub const MSYM_POINT_GROUP_TYPE_Dnh: msym_point_group_type_t = 8;
pub const MSYM_POINT_GROUP_TYPE_Dnd: msym_point_group_type_t = 9;
pub const MSYM_POINT_GROUP_TYPE_Sn: msym_point_group_type_t = 10;
pub const MSYM_POINT_GROUP_TYPE_T: msym_point_group_type_t = 11;
pub const MSYM_POINT_GROUP_TYPE_Td: msym_point_group_type_t = 12;
pub const MSYM_POINT_GROUP_TYPE_Th: msym_point_group_type_t = 13;
pub const MSYM_POINT_GROUP_TYPE_O: msym_point_group_type_t = 14;
pub const MSYM_POINT_GROUP_TYPE_Oh: msym_point_group_type_t = 15;
pub const MSYM_POINT_GROUP_TYPE_I: msym_point_group_type_t = 16;
pub const MSYM_POINT_GROUP_TYPE_Ih: msym_point_group_type_t = 17;

pub type msym_basis_type_t = c_int;

pub const MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC: msym_basis_type_t = 0;
pub const MSYM_BASIS_TYPE_CARTESIAN: msym_basis_type_t = 1;

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct msym_symmetry_operation_t {
    pub type_: msym_symmetry_operation_type_t,
    pub order: c_int,
    pub power: c_int,
    pub orientation: msym_symmetry_operation_orientation_t,
    pub v: [c_double; 3],
    pub cla: c_int,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct msym_thresholds_t {
    pub zero: c_double,
    pub geometry: c_double,
    pub angle: c_double,
    pub equivalence: c_double,
    pub eigfact: c_double,
    pub permutation: c_double,
    pub orthogonalization: c_double,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct msym_element_t {
    pub id: *mut c_void,
    pub m: c_double,
    pub v: [c_double; 3],
    pub n: c_int,
    pub name: [c_char; 4],
}

#[repr(C)]
pub struct msym_equivalence_set_t {
    pub elements: *mut *mut msym_element_t,
    pub err: c_double,
    pub length: c_int,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct msym_real_spherical_harmonic_t {
    pub n: c_int,
    pub l: c_int,
    pub m: c_int,
}

#[repr(C)]
pub struct msym_basis_function_t {
    pub id: *mut c_void,
    pub type_: msym_basis_type_t,
    pub element: *mut msym_element_t,
    pub f: msym_basis_function_union_t,
    pub name: [c_char; 8],
}

#[repr(C)]
#[derive(Clone, Copy)]
pub union msym_basis_function_union_t {
    pub rsh: msym_real_spherical_harmonic_t,
}

#[repr(C)]
pub struct msym_partner_function_t {
    pub i: c_int,
    pub d: c_int,
}

#[repr(C)]
pub struct msym_salc_t {
    pub d: c_int,
    pub fl: c_int,
    pub pf: *mut c_void, // double[d][fl]
    pub f: *mut *mut msym_basis_function_t,
}

#[repr(C)]
pub struct msym_subrepresentation_space_t {
    pub s: c_int,
    pub salcl: c_int,
    pub salc: *mut msym_salc_t,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct msym_symmetry_species_t {
    pub d: c_int,
    pub r: c_int,
    pub name: [c_char; 8],
}

#[repr(C)]
pub struct msym_character_table_t {
    pub d: c_int,
    pub classc: *mut c_int,
    pub sops: *mut *mut msym_symmetry_operation_t,
    pub s: *mut msym_symmetry_species_t,
    pub table: *mut c_void, // double[d][d]
}

#[repr(C)]
pub struct msym_subgroup_t {
    pub type_: msym_point_group_type_t,
    pub n: c_int,
    pub order: c_int,
    pub primary: *mut msym_symmetry_operation_t,
    pub sops: *mut *mut msym_symmetry_operation_t,
    pub generators: [*mut msym_subgroup_t; 2],
    pub name: [c_char; 8],
}

// ---------------------------------------------------------------------------
// Functions
// ---------------------------------------------------------------------------

unsafe extern "C" {
    // Context
    pub fn msymCreateContext() -> msym_context;
    pub fn msymReleaseContext(ctx: msym_context) -> msym_error_t;

    // Thresholds
    pub fn msymGetDefaultThresholds() -> *const msym_thresholds_t;
    pub fn msymSetThresholds(ctx: msym_context, thresholds: *const msym_thresholds_t)
        -> msym_error_t;
    pub fn msymGetThresholds(
        ctx: msym_context,
        thresholds: *mut *const msym_thresholds_t,
    ) -> msym_error_t;

    // Elements
    pub fn msymSetElements(
        ctx: msym_context,
        length: c_int,
        elements: *mut msym_element_t,
    ) -> msym_error_t;
    pub fn msymGetElements(
        ctx: msym_context,
        length: *mut c_int,
        elements: *mut *mut msym_element_t,
    ) -> msym_error_t;
    pub fn msymGenerateElements(
        ctx: msym_context,
        length: c_int,
        elements: *mut msym_element_t,
    ) -> msym_error_t;
    pub fn msymSymmetrizeElements(ctx: msym_context, err: *mut c_double) -> msym_error_t;
    pub fn msymApplyTranslation(
        ctx: msym_context,
        element: *mut msym_element_t,
        v: *mut c_double,
    ) -> msym_error_t;

    // Basis functions
    pub fn msymSetBasisFunctions(
        ctx: msym_context,
        length: c_int,
        basis: *mut msym_basis_function_t,
    ) -> msym_error_t;
    pub fn msymGetBasisFunctions(
        ctx: msym_context,
        length: *mut c_int,
        basis: *mut *mut msym_basis_function_t,
    ) -> msym_error_t;

    // Point group
    pub fn msymGetPointGroupType(
        ctx: msym_context,
        t: *mut msym_point_group_type_t,
        n: *mut c_int,
    ) -> msym_error_t;
    pub fn msymSetPointGroupByName(ctx: msym_context, name: *const c_char) -> msym_error_t;
    pub fn msymSetPointGroupByType(
        ctx: msym_context,
        type_: msym_point_group_type_t,
        n: c_int,
    ) -> msym_error_t;
    pub fn msymGetPointGroupName(
        ctx: msym_context,
        l: c_int,
        buf: *mut c_char,
    ) -> msym_error_t;
    pub fn msymGetSubgroups(
        ctx: msym_context,
        l: *mut c_int,
        subgroups: *mut *const msym_subgroup_t,
    ) -> msym_error_t;
    pub fn msymSelectSubgroup(
        ctx: msym_context,
        subgroup: *const msym_subgroup_t,
    ) -> msym_error_t;

    // Symmetry operations
    pub fn msymGetSymmetryOperations(
        ctx: msym_context,
        sopsl: *mut c_int,
        sops: *mut *const msym_symmetry_operation_t,
    ) -> msym_error_t;
    pub fn msymFindSymmetry(ctx: msym_context) -> msym_error_t;

    // Equivalence sets
    pub fn msymGetEquivalenceSets(
        ctx: msym_context,
        l: *mut c_int,
        es: *mut *const msym_equivalence_set_t,
    ) -> msym_error_t;
    pub fn msymGetEquivalenceSetByElement(
        ctx: msym_context,
        element: *mut msym_element_t,
        es: *mut *const msym_equivalence_set_t,
    ) -> msym_error_t;
    pub fn msymFindEquivalenceSets(ctx: msym_context) -> msym_error_t;
    pub fn msymFindEquivalenceSetPermutations(ctx: msym_context) -> msym_error_t;

    // Subrepresentation spaces & character tables
    pub fn msymGetSubrepresentationSpaces(
        ctx: msym_context,
        l: *mut c_int,
        srs: *mut *const msym_subrepresentation_space_t,
    ) -> msym_error_t;
    pub fn msymGetCharacterTable(
        ctx: msym_context,
        ct: *mut *const msym_character_table_t,
    ) -> msym_error_t;
    pub fn msymGenerateSubrepresentationSpaces(ctx: msym_context) -> msym_error_t;

    // SALCs
    pub fn msymSymmetrizeWavefunctions(
        ctx: msym_context,
        l: c_int,
        c: *mut c_double,
        species: *mut c_int,
        pf: *mut msym_partner_function_t,
    ) -> msym_error_t;
    pub fn msymGetSALCs(
        ctx: msym_context,
        l: c_int,
        c: *mut c_double,
        species: *mut c_int,
        pf: *mut msym_partner_function_t,
    ) -> msym_error_t;
    pub fn msymSymmetrySpeciesComponents(
        ctx: msym_context,
        wfl: c_int,
        wf: *mut c_double,
        sl: c_int,
        s: *mut c_double,
    ) -> msym_error_t;

    // Geometry & alignment
    pub fn msymGetCenterOfMass(ctx: msym_context, v: *mut c_double) -> msym_error_t;
    pub fn msymSetCenterOfMass(ctx: msym_context, v: *mut c_double) -> msym_error_t;
    pub fn msymGetRadius(ctx: msym_context, radius: *mut c_double) -> msym_error_t;
    pub fn msymGetGeometry(ctx: msym_context, geometry: *mut msym_geometry_t) -> msym_error_t;
    pub fn msymGetPrincipalMoments(ctx: msym_context, eigval: *mut c_double) -> msym_error_t;
    pub fn msymGetPrincipalAxes(ctx: msym_context, eigvec: *mut [c_double; 3]) -> msym_error_t;
    pub fn msymAlignAxes(ctx: msym_context) -> msym_error_t;
    pub fn msymGetAlignmentAxes(
        ctx: msym_context,
        primary: *mut c_double,
        secondary: *mut c_double,
    ) -> msym_error_t;
    pub fn msymSetAlignmentAxes(
        ctx: msym_context,
        primary: *mut c_double,
        secondary: *mut c_double,
    ) -> msym_error_t;
    pub fn msymGetAlignmentTransform(
        ctx: msym_context,
        transform: *mut [c_double; 3],
    ) -> msym_error_t;
    pub fn msymSetAlignmentTransform(
        ctx: msym_context,
        transform: *mut [c_double; 3],
    ) -> msym_error_t;

    // Error
    pub fn msymErrorString(error: msym_error_t) -> *const c_char;
    pub fn msymGetErrorDetails() -> *const c_char;
}
