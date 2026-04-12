//! Miscellanious type definitions.

use std::ffi::CStr;
use std::os::raw::c_int;
use std::str::FromStr;
use std::{fmt, ptr};

use nalgebra::Vector3;
use umol_msym_sys as ffi;

use crate::error::ParseError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SchoenfliesSymbol {
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

impl SchoenfliesSymbol {
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

    fn parse_parametric<F>(s: &str, f: F) -> Option<SchoenfliesSymbol>
    where
        F: FnOnce(u32, &str) -> Option<SchoenfliesSymbol>,
    {
        let digit_end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        if digit_end == 0 {
            return None;
        }
        let n: u32 = s[..digit_end].parse().ok()?;
        let suffix = &s[digit_end..];
        f(n, suffix)
    }
}

impl FromStr for SchoenfliesSymbol {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, ParseError> {
        match s {
            "Ci" => return Ok(Self::Ci),
            "Cs" => return Ok(Self::Cs),
            "T" => return Ok(Self::T),
            "Td" => return Ok(Self::Td),
            "Th" => return Ok(Self::Th),
            "O" => return Ok(Self::O),
            "Oh" => return Ok(Self::Oh),
            "I" => return Ok(Self::I),
            "Ih" => return Ok(Self::Ih),
            "Kh" => return Ok(Self::Kh),
            "K" => return Ok(Self::K),
            "C∞v" | "Coov" | "C0v" => return Ok(Self::Coov),
            "D∞h" | "Dooh" | "D0h" => return Ok(Self::Dooh),
            _ => {}
        }

        let result = if let Some(rest) = s.strip_prefix('C') {
            Self::parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Cn(n)),
                "h" => Some(Self::Cnh(n)),
                "v" => Some(Self::Cnv(n)),
                _ => None,
            })
        } else if let Some(rest) = s.strip_prefix('D') {
            Self::parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Dn(n)),
                "h" => Some(Self::Dnh(n)),
                "d" => Some(Self::Dnd(n)),
                _ => None,
            })
        } else if let Some(rest) = s.strip_prefix('S') {
            Self::parse_parametric(rest, |n, suffix| match suffix {
                "" => Some(Self::Sn(n)),
                _ => None,
            })
        } else {
            None
        };

        result.ok_or_else(|| ParseError(s.to_string()))
    }
}

impl fmt::Display for SchoenfliesSymbol {
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
pub enum MolecularShape {
    Unknown,
    Spherical,
    Linear,
    PlanarRegular,
    PlanarIrregular,
    PolyhedralProlate,
    PolyhedralOblate,
    Asymmetric,
}

impl From<ffi::msym_geometry_t> for MolecularShape {
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

#[derive(Debug, Clone)]
pub struct SymmetryCenter {
    pub atomic_number: i32,
    pub mass: f64,
    pub position: Vector3<f64>,
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
            v: self.position.into(),
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
            position: Vector3::from(e.v),
            name,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EquivalenceSet {
    pub centers: Vec<SymmetryCenter>,
}
