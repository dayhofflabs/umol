use std::ffi::CStr;
use std::fmt;

use umol_msym_sys as ffi;

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Geometry {
    Unknown,
    Spherical,
    Linear,
    PlanarRegular,
    PlanarIrregular,
    PolyhedralProlate,
    PolyhedralOblate,
    Asymmetric,
}

impl From<ffi::msym_geometry_t> for Geometry {
    fn from(g: ffi::msym_geometry_t) -> Self {
        match g {
            ffi::MSYM_GEOMETRY_SPHERICAL => Self::Spherical,
            ffi::MSYM_GEOMETRY_LINEAR => Self::Linear,
            ffi::MSYM_GEOMETRY_PLANAR_REGULAR => Self::PlanarRegular,
            ffi::MSYM_GEOMETRY_PLANAR_IRREGULAR => Self::PlanarIrregular,
            ffi::MSYM_GEOMETRY_POLYHEDRAL_PROLATE => Self::PolyhedralProlate,
            ffi::MSYM_GEOMETRY_POLYHEDRAL_OBLATE => Self::PolyhedralOblate,
            ffi::MSYM_GEOMETRY_ASSYMETRIC => Self::Asymmetric,
            _ => Self::Unknown,
        }
    }
}

// ---------------------------------------------------------------------------
// Point group type
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PointGroupType {
    Kh,
    K,
    Ci,
    Cs,
    Cn,
    Cnh,
    Cnv,
    Dn,
    Dnh,
    Dnd,
    Sn,
    T,
    Td,
    Th,
    O,
    Oh,
    I,
    Ih,
}

impl From<ffi::msym_point_group_type_t> for PointGroupType {
    fn from(t: ffi::msym_point_group_type_t) -> Self {
        match t {
            ffi::MSYM_POINT_GROUP_TYPE_Kh => Self::Kh,
            ffi::MSYM_POINT_GROUP_TYPE_K => Self::K,
            ffi::MSYM_POINT_GROUP_TYPE_Ci => Self::Ci,
            ffi::MSYM_POINT_GROUP_TYPE_Cs => Self::Cs,
            ffi::MSYM_POINT_GROUP_TYPE_Cn => Self::Cn,
            ffi::MSYM_POINT_GROUP_TYPE_Cnh => Self::Cnh,
            ffi::MSYM_POINT_GROUP_TYPE_Cnv => Self::Cnv,
            ffi::MSYM_POINT_GROUP_TYPE_Dn => Self::Dn,
            ffi::MSYM_POINT_GROUP_TYPE_Dnh => Self::Dnh,
            ffi::MSYM_POINT_GROUP_TYPE_Dnd => Self::Dnd,
            ffi::MSYM_POINT_GROUP_TYPE_Sn => Self::Sn,
            ffi::MSYM_POINT_GROUP_TYPE_T => Self::T,
            ffi::MSYM_POINT_GROUP_TYPE_Td => Self::Td,
            ffi::MSYM_POINT_GROUP_TYPE_Th => Self::Th,
            ffi::MSYM_POINT_GROUP_TYPE_O => Self::O,
            ffi::MSYM_POINT_GROUP_TYPE_Oh => Self::Oh,
            ffi::MSYM_POINT_GROUP_TYPE_I => Self::I,
            ffi::MSYM_POINT_GROUP_TYPE_Ih => Self::Ih,
            _ => Self::Ci, // unreachable in practice
        }
    }
}

impl PointGroupType {
    pub(crate) fn to_ffi(self) -> ffi::msym_point_group_type_t {
        match self {
            Self::Kh => ffi::MSYM_POINT_GROUP_TYPE_Kh,
            Self::K => ffi::MSYM_POINT_GROUP_TYPE_K,
            Self::Ci => ffi::MSYM_POINT_GROUP_TYPE_Ci,
            Self::Cs => ffi::MSYM_POINT_GROUP_TYPE_Cs,
            Self::Cn => ffi::MSYM_POINT_GROUP_TYPE_Cn,
            Self::Cnh => ffi::MSYM_POINT_GROUP_TYPE_Cnh,
            Self::Cnv => ffi::MSYM_POINT_GROUP_TYPE_Cnv,
            Self::Dn => ffi::MSYM_POINT_GROUP_TYPE_Dn,
            Self::Dnh => ffi::MSYM_POINT_GROUP_TYPE_Dnh,
            Self::Dnd => ffi::MSYM_POINT_GROUP_TYPE_Dnd,
            Self::Sn => ffi::MSYM_POINT_GROUP_TYPE_Sn,
            Self::T => ffi::MSYM_POINT_GROUP_TYPE_T,
            Self::Td => ffi::MSYM_POINT_GROUP_TYPE_Td,
            Self::Th => ffi::MSYM_POINT_GROUP_TYPE_Th,
            Self::O => ffi::MSYM_POINT_GROUP_TYPE_O,
            Self::Oh => ffi::MSYM_POINT_GROUP_TYPE_Oh,
            Self::I => ffi::MSYM_POINT_GROUP_TYPE_I,
            Self::Ih => ffi::MSYM_POINT_GROUP_TYPE_Ih,
        }
    }
}

// ---------------------------------------------------------------------------
// Symmetry operation
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SymmetryOperationType {
    Identity,
    ProperRotation,
    ImproperRotation,
    Reflection,
    Inversion,
}

impl From<ffi::msym_symmetry_operation_type_t> for SymmetryOperationType {
    fn from(t: ffi::msym_symmetry_operation_type_t) -> Self {
        match t {
            ffi::MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION => Self::ProperRotation,
            ffi::MSYM_SYMMETRY_OPERATION_TYPE_IMPROPER_ROTATION => Self::ImproperRotation,
            ffi::MSYM_SYMMETRY_OPERATION_TYPE_REFLECTION => Self::Reflection,
            ffi::MSYM_SYMMETRY_OPERATION_TYPE_INVERSION => Self::Inversion,
            _ => Self::Identity,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SymmetryOperationOrientation {
    None,
    Horizontal,
    Vertical,
    Dihedral,
}

impl From<ffi::msym_symmetry_operation_orientation_t> for SymmetryOperationOrientation {
    fn from(o: ffi::msym_symmetry_operation_orientation_t) -> Self {
        match o {
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_HORIZONTAL => Self::Horizontal,
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_VERTICAL => Self::Vertical,
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_DIHEDRAL => Self::Dihedral,
            _ => Self::None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SymmetryOperation {
    pub type_: SymmetryOperationType,
    pub order: i32,
    pub power: i32,
    pub orientation: SymmetryOperationOrientation,
    pub vector: [f64; 3],
    pub class: i32,
}

impl From<&ffi::msym_symmetry_operation_t> for SymmetryOperation {
    fn from(op: &ffi::msym_symmetry_operation_t) -> Self {
        Self {
            type_: op.type_.into(),
            order: op.order,
            power: op.power,
            orientation: op.orientation.into(),
            vector: op.v,
            class: op.cla,
        }
    }
}

// ---------------------------------------------------------------------------
// Element (atom)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct SymmetryElement {
    pub atomic_number: i32,
    pub mass: f64,
    pub position: [f64; 3],
    pub name: String,
}

impl SymmetryElement {
    pub(crate) fn to_ffi(&self) -> ffi::msym_element_t {
        let mut name = [0i8; 4];
        for (i, &b) in self.name.as_bytes().iter().take(3).enumerate() {
            name[i] = b as i8;
        }
        ffi::msym_element_t {
            id: std::ptr::null_mut(),
            m: self.mass,
            v: self.position,
            n: self.atomic_number,
            name,
        }
    }

    pub(crate) fn from_ffi(e: &ffi::msym_element_t) -> Self {
        let name = unsafe { CStr::from_ptr(e.name.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        Self {
            atomic_number: e.n,
            mass: e.m,
            position: e.v,
            name,
        }
    }
}

// ---------------------------------------------------------------------------
// Thresholds
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub zero: f64,
    pub geometry: f64,
    pub angle: f64,
    pub equivalence: f64,
    pub eigfact: f64,
    pub permutation: f64,
    pub orthogonalization: f64,
}

impl Thresholds {
    pub fn defaults() -> Self {
        unsafe {
            let t = ffi::msymGetDefaultThresholds();
            Self::from_ffi(&*t)
        }
    }

    pub(crate) fn from_ffi(t: &ffi::msym_thresholds_t) -> Self {
        Self {
            zero: t.zero,
            geometry: t.geometry,
            angle: t.angle,
            equivalence: t.equivalence,
            eigfact: t.eigfact,
            permutation: t.permutation,
            orthogonalization: t.orthogonalization,
        }
    }

    pub(crate) fn to_ffi(&self) -> ffi::msym_thresholds_t {
        ffi::msym_thresholds_t {
            zero: self.zero,
            geometry: self.geometry,
            angle: self.angle,
            equivalence: self.equivalence,
            eigfact: self.eigfact,
            permutation: self.permutation,
            orthogonalization: self.orthogonalization,
        }
    }
}

// ---------------------------------------------------------------------------
// Irreducible representation
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Irrep {
    pub name: String,
    pub dimension: i32,
    pub index: usize,
}

impl fmt::Display for Irrep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.name)
    }
}

// ---------------------------------------------------------------------------
// Character table
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct CharacterTable {
    pub irreps: Vec<Irrep>,
    pub class_sizes: Vec<i32>,
    pub class_operations: Vec<SymmetryOperation>,
    pub characters: Vec<Vec<f64>>,
    pub order: usize,
}

impl CharacterTable {
    pub(crate) unsafe fn from_ffi(ct: &ffi::msym_character_table_t) -> Self {
        let d = ct.d as usize;

        let class_sizes = std::slice::from_raw_parts(ct.classc, d).to_vec();

        let class_operations: Vec<SymmetryOperation> = (0..d)
            .map(|i| SymmetryOperation::from(&**ct.sops.add(i)))
            .collect();

        let species = std::slice::from_raw_parts(ct.s, d);
        let irreps: Vec<Irrep> = species
            .iter()
            .enumerate()
            .map(|(i, s)| Irrep {
                name: CStr::from_ptr(s.name.as_ptr())
                    .to_string_lossy()
                    .into_owned(),
                dimension: s.d,
                index: i,
            })
            .collect();

        let table_ptr = ct.table as *const f64;
        let characters: Vec<Vec<f64>> = (0..d)
            .map(|i| std::slice::from_raw_parts(table_ptr.add(i * d), d).to_vec())
            .collect();

        let order = class_sizes.iter().sum::<i32>() as usize;

        Self {
            irreps,
            class_sizes,
            class_operations,
            characters,
            order,
        }
    }

    pub fn direct_product(&self, a: &Irrep, b: &Irrep) -> Vec<(Irrep, u32)> {
        let d = self.irreps.len();
        let h = self.order as f64;

        // χ_{a⊗b}(R) = χ_a(R) · χ_b(R)
        let product_chars: Vec<f64> = (0..d)
            .map(|c| self.characters[a.index][c] * self.characters[b.index][c])
            .collect();

        // Decompose: n_i = (1/h) Σ_c |C_c| · χ_i(c)* · χ_{a⊗b}(c)
        let mut result = Vec::new();
        for irrep in &self.irreps {
            let n: f64 = (0..d)
                .map(|c| {
                    self.class_sizes[c] as f64 * self.characters[irrep.index][c] * product_chars[c]
                })
                .sum::<f64>()
                / h;
            let n_rounded = n.round() as u32;
            if n_rounded > 0 {
                result.push((irrep.clone(), n_rounded));
            }
        }
        result
    }

    pub fn contains_totally_symmetric(&self, a: &Irrep, b: &Irrep, c: &Irrep) -> bool {
        let d = self.irreps.len();
        let h = self.order as f64;

        // n_{A1} in a⊗b⊗c = (1/h) Σ_c |C_c| · χ_a(c) · χ_b(c) · χ_c(c)
        let n: f64 = (0..d)
            .map(|cls| {
                self.class_sizes[cls] as f64
                    * self.characters[a.index][cls]
                    * self.characters[b.index][cls]
                    * self.characters[c.index][cls]
            })
            .sum::<f64>()
            / h;
        n.round() as u32 > 0
    }
}

// ---------------------------------------------------------------------------
// Equivalence set
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct EquivalenceSet {
    pub elements: Vec<SymmetryElement>,
    pub max_error: f64,
}

// ---------------------------------------------------------------------------
// Basis function
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BasisType {
    RealSphericalHarmonic,
    Cartesian,
}

#[derive(Debug, Clone)]
pub struct BasisFunction {
    pub type_: BasisType,
    pub element_index: usize,
    pub n: i32,
    pub l: i32,
    pub m: i32,
    pub name: String,
}

// ---------------------------------------------------------------------------
// SALC
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Salc {
    pub irrep: Irrep,
    pub dimension: i32,
    pub coefficients: Vec<Vec<f64>>,
    pub basis_indices: Vec<usize>,
}

// ---------------------------------------------------------------------------
// Subrepresentation space
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct SubrepresentationSpace {
    pub irrep: Irrep,
    pub salcs: Vec<Salc>,
}
