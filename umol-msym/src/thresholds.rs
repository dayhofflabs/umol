use umol_msym_sys as ffi;

pub(crate) const REDUCTION_INTEGRALITY: f64 = 0.01;
pub(crate) const COMPLEX_IRREP_NORM: f64 = 0.5;
pub(crate) const CHARACTER_DISPLAY_ROUNDING: f64 = 1e-6;
pub(crate) const DEFAULT_PROJECTION: f64 = 1e-8;

#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub zero: f64,
    pub geometry: f64,
    pub angle: f64,
    pub equivalence: f64,
    pub jacobi: f64,
    pub permutation: f64,
    pub orthogonalization: f64,
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
            projection: DEFAULT_PROJECTION,
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
