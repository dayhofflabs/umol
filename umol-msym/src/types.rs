use std::ffi::CStr;
use std::os::raw::c_int;
use std::{fmt, ptr};

use umol_msym_sys as ffi;

/// Structured Schoenflies symbol identifying a point group.
///
/// The rotational axis order is encoded on the variant for families that take it.
/// Cubic, icosahedral, and continuous groups carry no parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SchoenfliesLabel {
    Ci,
    Cs,
    Cn(u32),
    Cnh(u32),
    Cnv(u32),
    Sn(u32),
    Dn(u32),
    Dnh(u32),
    Dnd(u32),
    T,
    Td,
    Th,
    O,
    Oh,
    I,
    Ih,
    Coov,
    Dooh,
    K,
    Kh,
}

impl SchoenfliesLabel {
    pub(crate) fn from_ffi(t: ffi::msym_point_group_type_t, n: c_int) -> Self {
        let n = n as u32;
        match t {
            ffi::MSYM_POINT_GROUP_TYPE_Kh => Self::Kh,
            ffi::MSYM_POINT_GROUP_TYPE_K => Self::K,
            ffi::MSYM_POINT_GROUP_TYPE_Ci => Self::Ci,
            ffi::MSYM_POINT_GROUP_TYPE_Cs => Self::Cs,
            ffi::MSYM_POINT_GROUP_TYPE_Cn => Self::Cn(n),
            ffi::MSYM_POINT_GROUP_TYPE_Cnh => Self::Cnh(n),
            ffi::MSYM_POINT_GROUP_TYPE_Cnv if n == 0 => Self::Coov,
            ffi::MSYM_POINT_GROUP_TYPE_Cnv => Self::Cnv(n),
            ffi::MSYM_POINT_GROUP_TYPE_Dn => Self::Dn(n),
            ffi::MSYM_POINT_GROUP_TYPE_Dnh if n == 0 => Self::Dooh,
            ffi::MSYM_POINT_GROUP_TYPE_Dnh => Self::Dnh(n),
            ffi::MSYM_POINT_GROUP_TYPE_Dnd => Self::Dnd(n),
            ffi::MSYM_POINT_GROUP_TYPE_Sn => Self::Sn(n),
            ffi::MSYM_POINT_GROUP_TYPE_T => Self::T,
            ffi::MSYM_POINT_GROUP_TYPE_Td => Self::Td,
            ffi::MSYM_POINT_GROUP_TYPE_Th => Self::Th,
            ffi::MSYM_POINT_GROUP_TYPE_O => Self::O,
            ffi::MSYM_POINT_GROUP_TYPE_Oh => Self::Oh,
            ffi::MSYM_POINT_GROUP_TYPE_I => Self::I,
            ffi::MSYM_POINT_GROUP_TYPE_Ih => Self::Ih,
            _ => Self::Ci,
        }
    }

    pub(crate) fn to_ffi(self) -> (ffi::msym_point_group_type_t, c_int) {
        match self {
            Self::Kh => (ffi::MSYM_POINT_GROUP_TYPE_Kh, 0),
            Self::K => (ffi::MSYM_POINT_GROUP_TYPE_K, 0),
            Self::Ci => (ffi::MSYM_POINT_GROUP_TYPE_Ci, 0),
            Self::Cs => (ffi::MSYM_POINT_GROUP_TYPE_Cs, 0),
            Self::Cn(n) => (ffi::MSYM_POINT_GROUP_TYPE_Cn, n as c_int),
            Self::Cnh(n) => (ffi::MSYM_POINT_GROUP_TYPE_Cnh, n as c_int),
            Self::Cnv(n) => (ffi::MSYM_POINT_GROUP_TYPE_Cnv, n as c_int),
            Self::Sn(n) => (ffi::MSYM_POINT_GROUP_TYPE_Sn, n as c_int),
            Self::Dn(n) => (ffi::MSYM_POINT_GROUP_TYPE_Dn, n as c_int),
            Self::Dnh(n) => (ffi::MSYM_POINT_GROUP_TYPE_Dnh, n as c_int),
            Self::Dnd(n) => (ffi::MSYM_POINT_GROUP_TYPE_Dnd, n as c_int),
            Self::T => (ffi::MSYM_POINT_GROUP_TYPE_T, 3),
            Self::Td => (ffi::MSYM_POINT_GROUP_TYPE_Td, 3),
            Self::Th => (ffi::MSYM_POINT_GROUP_TYPE_Th, 3),
            Self::O => (ffi::MSYM_POINT_GROUP_TYPE_O, 4),
            Self::Oh => (ffi::MSYM_POINT_GROUP_TYPE_Oh, 4),
            Self::I => (ffi::MSYM_POINT_GROUP_TYPE_I, 5),
            Self::Ih => (ffi::MSYM_POINT_GROUP_TYPE_Ih, 5),
            Self::Coov => (ffi::MSYM_POINT_GROUP_TYPE_Cnv, 0),
            Self::Dooh => (ffi::MSYM_POINT_GROUP_TYPE_Dnh, 0),
        }
    }

    /// Parse a Schoenflies symbol string (e.g. "C2v", "Td", "D6h").
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "Ci" => return Some(Self::Ci),
            "Cs" => return Some(Self::Cs),
            "T" => return Some(Self::T),
            "Td" => return Some(Self::Td),
            "Th" => return Some(Self::Th),
            "O" => return Some(Self::O),
            "Oh" => return Some(Self::Oh),
            "I" => return Some(Self::I),
            "Ih" => return Some(Self::Ih),
            "Kh" => return Some(Self::Kh),
            "K" => return Some(Self::K),
            "C∞v" | "Coov" | "C0v" => return Some(Self::Coov),
            "D∞h" | "Dooh" | "D0h" => return Some(Self::Dooh),
            _ => {}
        }

        if let Some(rest) = s.strip_prefix('C') {
            parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Cn(n)),
                "h" => Some(Self::Cnh(n)),
                "v" => Some(Self::Cnv(n)),
                _ => None,
            })
        } else if let Some(rest) = s.strip_prefix('D') {
            parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Dn(n)),
                "h" => Some(Self::Dnh(n)),
                "d" => Some(Self::Dnd(n)),
                _ => None,
            })
        } else if let Some(rest) = s.strip_prefix('S') {
            parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Sn(n)),
                _ => None,
            })
        } else {
            None
        }
    }

}

fn parse_parametric<F>(s: &str, f: F) -> Option<SchoenfliesLabel>
where
    F: FnOnce(u32, &str) -> Option<SchoenfliesLabel>,
{
    let digit_end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    if digit_end == 0 {
        return None;
    }
    let n: u32 = s[..digit_end].parse().ok()?;
    let suffix = &s[digit_end..];
    f(n, suffix)
}

impl fmt::Display for SchoenfliesLabel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ci => write!(f, "Ci"),
            Self::Cs => write!(f, "Cs"),
            Self::Cn(n) => write!(f, "C{n}"),
            Self::Cnh(n) => write!(f, "C{n}h"),
            Self::Cnv(n) => write!(f, "C{n}v"),
            Self::Sn(n) => write!(f, "S{n}"),
            Self::Dn(n) => write!(f, "D{n}"),
            Self::Dnh(n) => write!(f, "D{n}h"),
            Self::Dnd(n) => write!(f, "D{n}d"),
            Self::T => write!(f, "T"),
            Self::Td => write!(f, "Td"),
            Self::Th => write!(f, "Th"),
            Self::O => write!(f, "O"),
            Self::Oh => write!(f, "Oh"),
            Self::I => write!(f, "I"),
            Self::Ih => write!(f, "Ih"),
            Self::Coov => write!(f, "C∞v"),
            Self::Dooh => write!(f, "D∞h"),
            Self::K => write!(f, "K"),
            Self::Kh => write!(f, "Kh"),
        }
    }
}

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SymmetryOpKind {
    Identity,
    ProperRotation,
    ImproperRotation,
    Reflection,
    Inversion,
}

impl From<ffi::msym_symmetry_operation_type_t> for SymmetryOpKind {
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
pub enum SymmetryOpOrientation {
    None,
    Horizontal,
    Vertical,
    Dihedral,
}

impl From<ffi::msym_symmetry_operation_orientation_t> for SymmetryOpOrientation {
    fn from(o: ffi::msym_symmetry_operation_orientation_t) -> Self {
        match o {
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_HORIZONTAL => Self::Horizontal,
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_VERTICAL => Self::Vertical,
            ffi::MSYM_SYMMETRY_OPERATION_ORIENTATION_DIHEDRAL => Self::Dihedral,
            _ => Self::None,
        }
    }
}

/// Per-operation abstract data stored on `PointGroup`, indexed by op position.
#[derive(Debug, Clone)]
pub(crate) struct OpData {
    pub kind: SymmetryOpKind,
    pub order: i32,
    pub power: i32,
    pub orientation: SymmetryOpOrientation,
    pub class: usize,
}

#[derive(Debug, Clone)]
pub struct SymmetryCenter {
    pub atomic_number: i32,
    pub mass: f64,
    pub position: [f64; 3],
    pub name: String,
}

impl SymmetryCenter {
    pub(crate) fn to_ffi(&self) -> ffi::msym_element_t {
        let mut name = [0i8; 4];
        for (i, &b) in self.name.as_bytes().iter().take(3).enumerate() {
            name[i] = b as i8;
        }
        ffi::msym_element_t {
            id: ptr::null_mut(),
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

#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub zero: f64,
    pub geometry: f64,
    pub angle: f64,
    pub equivalence: f64,
    pub jacobi: f64,
    pub permutation: f64,
    pub orthogonalization: f64,
    /// SVD singular value cutoff for symmetry coordinate projection.
    pub projection: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        unsafe {
            let t = ffi::msymGetDefaultThresholds();
            Self::from_ffi(&*t)
        }
    }
}

impl Thresholds {
    pub(crate) fn from_ffi(t: &ffi::msym_thresholds_t) -> Self {
        Self {
            zero: t.zero,
            geometry: t.geometry,
            angle: t.angle,
            equivalence: t.equivalence,
            jacobi: t.eigfact,
            permutation: t.permutation,
            orthogonalization: t.orthogonalization,
            projection: 1e-8,
        }
    }

    pub(crate) fn to_ffi(self) -> ffi::msym_thresholds_t {
        ffi::msym_thresholds_t {
            zero: self.zero,
            geometry: self.geometry,
            angle: self.angle,
            equivalence: self.equivalence,
            eigfact: self.jacobi,
            permutation: self.permutation,
            orthogonalization: self.orthogonalization,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IrrepData {
    pub symbol: String,
    pub dimension: i32,
    pub index: usize,
    /// Character values per conjugacy class. Empty for linear groups.
    pub characters: Vec<f64>,
    /// Angular momentum quantum number (λ) for linear group irreps.
    pub lambda: Option<u32>,
    /// σ_v parity for Σ irreps of linear groups: true = Σ+, false = Σ-.
    pub sigma_v: Option<bool>,
    /// Gerade/ungerade for D∞h irreps.
    pub gerade: Option<bool>,
}

#[derive(Debug, Clone)]
pub struct EquivalenceSet {
    pub centers: Vec<SymmetryCenter>,
    pub max_error: f64,
}
