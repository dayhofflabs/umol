//! Symmetry operations and point group data.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::ffi::{CStr, CString};
use std::hash::{Hash, Hasher};
use std::os::raw::c_int;
use std::sync::{LazyLock, Mutex};
use std::{fmt, ptr, slice};

use nalgebra::{Matrix3, Vector3};
use umol_msym_sys::{
    self as ffi, MSYM_INVALID_CHARACTER_TABLE, MSYM_INVALID_INPUT, MSYM_MEMORY_ERROR,
};

use crate::error::{self, MsymError, ReductionError};
use crate::irrep::{Irrep, IrrepData, ReductionData};
use crate::linear;
use crate::thresholds::{CHARACTER_DISPLAY_ROUNDING, COMPLEX_IRREP_NORM, REDUCTION_INTEGRALITY};
use crate::types::SchoenfliesSymbol;

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

#[derive(Debug, Clone)]
pub(crate) struct SymmetryOpData {
    pub(crate) kind: SymmetryOpKind,
    pub(crate) order: i32,
    pub(crate) power: i32,
    pub(crate) orientation: SymmetryOpOrientation,
    pub(crate) class: usize,
}

static REGISTRY: LazyLock<Mutex<HashMap<SchoenfliesSymbol, &'static PointGroup>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

fn register(pg: PointGroup) -> &'static PointGroup {
    let label = pg.symbol;
    let mut map = REGISTRY.lock().unwrap();
    if let Some(&existing) = map.get(&label) {
        return existing;
    }
    let leaked: &'static PointGroup = Box::leak(Box::new(pg));
    map.insert(label, leaked);
    leaked
}

/// A symmetry operation handle: parent point group plus opaque positional index.
///
/// Equality and hashing use `(group, index)`, so σv and σv′ in C₂ᵥ compare
/// unequal even though their `(kind, order, power, orientation)` tuples agree.
#[derive(Clone, Copy)]
pub struct SymmetryOp {
    group: &'static PointGroup,
    index: usize,
}

impl SymmetryOp {
    pub fn group(&self) -> &'static PointGroup {
        self.group
    }

    pub fn index(&self) -> usize {
        self.index
    }

    fn data(&self) -> &SymmetryOpData {
        &self
            .group
            .finite
            .as_ref()
            .expect("SymmetryOp on linear group")
            .op_data[self.index]
    }

    pub fn kind(&self) -> SymmetryOpKind {
        self.data().kind
    }

    pub fn order(&self) -> i32 {
        self.data().order
    }

    pub fn power(&self) -> i32 {
        self.data().power
    }

    pub fn orientation(&self) -> SymmetryOpOrientation {
        self.data().orientation
    }

    pub fn class(&self) -> usize {
        self.data().class
    }

    pub fn is_proper(&self) -> bool {
        matches!(
            self.kind(),
            SymmetryOpKind::Identity | SymmetryOpKind::ProperRotation
        )
    }

    /// Character of this operation under the given irrep.
    pub fn character(&self, irrep: Irrep) -> f64 {
        debug_assert!(ptr::eq(irrep.group, self.group));
        irrep.data.characters()[self.class()]
    }
}

impl PartialEq for SymmetryOp {
    fn eq(&self, other: &Self) -> bool {
        ptr::eq(self.group, other.group) && self.index == other.index
    }
}

impl Eq for SymmetryOp {}

impl Hash for SymmetryOp {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (self.group as *const PointGroup).hash(state);
        self.index.hash(state);
    }
}

impl fmt::Debug for SymmetryOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SymmetryOp({}:{})", self.group.symbol, self.index)
    }
}

impl fmt::Display for SymmetryOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.data())
    }
}

impl fmt::Display for SymmetryOpData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind {
            SymmetryOpKind::Identity => return write!(f, "E"),
            SymmetryOpKind::Inversion => return write!(f, "i"),
            SymmetryOpKind::Reflection => {
                let s = match self.orientation {
                    SymmetryOpOrientation::Horizontal => "σh",
                    SymmetryOpOrientation::Vertical => "σv",
                    SymmetryOpOrientation::Dihedral => "σd",
                    SymmetryOpOrientation::None => "σ",
                };
                return write!(f, "{s}");
            }
            _ => {}
        }

        let letter = match self.kind {
            SymmetryOpKind::ProperRotation => "C",
            SymmetryOpKind::ImproperRotation => "S",
            _ => unreachable!(),
        };
        let prime = match self.orientation {
            SymmetryOpOrientation::Vertical => "'",
            SymmetryOpOrientation::Dihedral => "''",
            _ => "",
        };
        write!(f, "{letter}{}{prime}", self.order)?;
        if self.power != 1 {
            write!(f, "{}", superscript(self.power))?;
        }
        Ok(())
    }
}

fn superscript(n: i32) -> String {
    n.to_string()
        .chars()
        .map(|c| match c {
            '0' => '⁰',
            '1' => '¹',
            '2' => '²',
            '3' => '³',
            '4' => '⁴',
            '5' => '⁵',
            '6' => '⁶',
            '7' => '⁷',
            '8' => '⁸',
            '9' => '⁹',
            '-' => '⁻',
            _ => c,
        })
        .collect()
}

#[derive(Debug)]
pub(crate) struct FiniteGroupData {
    pub(crate) order: usize,
    pub(crate) op_data: Vec<SymmetryOpData>,
    #[allow(dead_code)]
    pub(crate) classes: Vec<Vec<usize>>,
    pub(crate) class_rep_indices: Vec<usize>,
    pub(crate) class_sizes: Vec<i32>,
    pub(crate) mul_table: Vec<Vec<usize>>,
    /// For each class `c`, the class containing the square of its representative.
    pub(crate) r_squared_class: Vec<usize>,
    /// Trace of ρ(R) per class, for the R³ representation.
    pub(crate) translation_chars: Vec<f64>,
    /// det·trace of ρ(R) per class, for the axial vector representation.
    pub(crate) rotation_chars: Vec<f64>,
    /// Sym² of the R³ representation, per class.
    pub(crate) quadratic_chars: Vec<f64>,
}

/// A molecular point group: orientation-independent algebraic data.
///
/// Point groups are `&'static` singletons cached in a process-global registry.
/// There is exactly one C₂ᵥ, one Td, etc. Access via macro (`group!(C2v)`) or
/// `PointGroup::parse("C2v")`.
#[derive(Debug)]
pub struct PointGroup {
    pub(crate) symbol: SchoenfliesSymbol,
    pub(crate) irreps: Vec<IrrepData>,
    pub(crate) finite: Option<FiniteGroupData>,
}

impl PointGroup {
    /// Construct point group from Schoenflies symbol
    pub fn parse(name: &str) -> Result<&'static PointGroup, MsymError> {
        let label: SchoenfliesSymbol = name.parse().map_err(|_| MsymError {
            code: MSYM_INVALID_INPUT,
            message: format!("cannot parse Schoenflies symbol '{name}'"),
        })?;
        Self::get_or_build(label)
    }

    /// Construct point group by parsed Schoenflies symbol
    pub fn from_symbol(label: SchoenfliesSymbol) -> Result<&'static PointGroup, MsymError> {
        Self::get_or_build(label)
    }

    fn get_or_build(label: SchoenfliesSymbol) -> Result<&'static PointGroup, MsymError> {
        {
            let map = REGISTRY.lock().unwrap();
            if let Some(&pg) = map.get(&label) {
                return Ok(pg);
            }
        }

        let pg = if matches!(label, SchoenfliesSymbol::Coov | SchoenfliesSymbol::Dooh) {
            Self::build_linear(label)
        } else if label == SchoenfliesSymbol::Cn(1) {
            Self::build_c1()
        } else if matches!(label, SchoenfliesSymbol::K | SchoenfliesSymbol::Kh) {
            return Err(MsymError {
                code: MSYM_INVALID_INPUT,
                message: format!("{label} is a continuous group and cannot be instantiated"),
            });
        } else {
            Self::build_finite(label)?
        };
        Ok(register(pg))
    }

    fn build_finite(label: SchoenfliesSymbol) -> Result<PointGroup, MsymError> {
        let name = label.to_string();
        let c_name = CString::new(name.as_str()).expect("no NUL in label");

        let ctx = unsafe { ffi::msymCreateContext() };
        if ctx.is_null() {
            return Err(MsymError {
                code: MSYM_MEMORY_ERROR,
                message: "failed to create msym context".into(),
            });
        }

        let result = build_finite_inner(ctx, label, &c_name);
        unsafe {
            ffi::msymReleaseContext(ctx);
        }
        result
    }

    fn build_c1() -> PointGroup {
        let op_data = vec![SymmetryOpData {
            kind: SymmetryOpKind::Identity,
            order: 1,
            power: 0,
            orientation: SymmetryOpOrientation::None,
            class: 0,
        }];
        PointGroup {
            symbol: SchoenfliesSymbol::Cn(1),
            irreps: vec![IrrepData {
                symbol: "A".into(),
                dimension: 1,
                reduction: ReductionData::Finite {
                    characters: vec![1.0],
                },
            }],
            finite: Some(FiniteGroupData {
                order: 1,
                op_data,
                classes: vec![vec![0]],
                class_rep_indices: vec![0],
                class_sizes: vec![1],
                mul_table: vec![vec![0]],
                r_squared_class: vec![0],
                translation_chars: vec![3.0],
                rotation_chars: vec![3.0],
                quadratic_chars: vec![6.0],
            }),
        }
    }

    fn build_linear(label: SchoenfliesSymbol) -> PointGroup {
        let is_dooh = label == SchoenfliesSymbol::Dooh;
        let lambda_symbols = ["Σ", "Π", "Δ", "Φ", "Γ", "H", "I"];
        let mut irrep_data = Vec::new();

        let gu_list: &[Option<bool>] = if is_dooh {
            &[Some(true), Some(false)]
        } else {
            &[None]
        };

        for &gerade in gu_list {
            let gu_suffix = match gerade {
                Some(true) => "g",
                Some(false) => "u",
                None => "",
            };

            for (lambda, &base) in lambda_symbols.iter().enumerate() {
                let lambda = lambda as u32;
                let dim = if lambda == 0 { 1 } else { 2 };

                if lambda == 0 {
                    for &sv in &[true, false] {
                        let sign = if sv { "+" } else { "-" };
                        irrep_data.push(IrrepData {
                            symbol: format!("{base}{sign}{gu_suffix}"),
                            dimension: dim,
                            reduction: ReductionData::Linear {
                                lambda,
                                sigma_v: Some(sv),
                                gerade,
                            },
                        });
                    }
                } else {
                    irrep_data.push(IrrepData {
                        symbol: format!("{base}{gu_suffix}"),
                        dimension: dim,
                        reduction: ReductionData::Linear {
                            lambda,
                            sigma_v: None,
                            gerade,
                        },
                    });
                }
            }
        }

        PointGroup {
            symbol: label,
            irreps: irrep_data,
            finite: None,
        }
    }
    pub fn symbol(&self) -> SchoenfliesSymbol {
        self.symbol
    }

    pub fn is_finite(&self) -> bool {
        self.finite.is_some()
    }

    pub fn is_abelian(&self) -> bool {
        matches!(
            self.symbol,
            SchoenfliesSymbol::Ci
                | SchoenfliesSymbol::Cs
                | SchoenfliesSymbol::Cn(_)
                | SchoenfliesSymbol::Cnh(_)
                | SchoenfliesSymbol::Sn(_)
                | SchoenfliesSymbol::Cnv(2)
                | SchoenfliesSymbol::Dn(2)
                | SchoenfliesSymbol::Dnh(2)
        )
    }

    pub fn is_cyclic(&self) -> bool {
        matches!(
            self.symbol,
            SchoenfliesSymbol::Ci
                | SchoenfliesSymbol::Cs
                | SchoenfliesSymbol::Cn(_)
                | SchoenfliesSymbol::Cnh(_)
                | SchoenfliesSymbol::Sn(_)
        )
    }

    pub fn is_dihedral(&self) -> bool {
        matches!(
            self.symbol,
            SchoenfliesSymbol::Dn(_)
                | SchoenfliesSymbol::Dnh(_)
                | SchoenfliesSymbol::Dnd(_)
                | SchoenfliesSymbol::Dooh
        )
    }

    pub fn is_cubic(&self) -> bool {
        matches!(
            self.symbol,
            SchoenfliesSymbol::T
                | SchoenfliesSymbol::Td
                | SchoenfliesSymbol::Th
                | SchoenfliesSymbol::O
                | SchoenfliesSymbol::Oh
        )
    }

    pub fn is_icosahedral(&self) -> bool {
        matches!(self.symbol, SchoenfliesSymbol::I | SchoenfliesSymbol::Ih)
    }

    pub fn is_simply_reducible(&self) -> bool {
        !self.is_icosahedral()
    }

    pub fn is_linear(&self) -> bool {
        matches!(
            self.symbol,
            SchoenfliesSymbol::Coov | SchoenfliesSymbol::Dooh
        )
    }

    pub fn is_spherical(&self) -> bool {
        matches!(self.symbol, SchoenfliesSymbol::K | SchoenfliesSymbol::Kh)
    }

    /// True if some irreps are complex conjugate pairs fused into real 2D representations.
    pub fn has_complex_irreps(&self) -> bool {
        let Some(f) = &self.finite else {
            return false;
        };
        let h = f.order as f64;
        self.irreps.iter().any(|ir| {
            if ir.dimension != 2 {
                return false;
            }
            let norm_sq: f64 = ir
                .characters()
                .iter()
                .zip(&f.class_sizes)
                .map(|(chi, &size)| size as f64 * chi * chi)
                .sum();
            (norm_sq - 2.0 * h).abs() < COMPLEX_IRREP_NORM
        })
    }

    /// True iff the group contains only proper operations (no reflections, inversions,
    /// or improper rotations). Chiral groups: Cn, Dn, T, O, I.
    pub fn is_chiral(&'static self) -> bool {
        self.ops().iter().all(|op| op.is_proper())
    }

    pub fn has_inversion(&self) -> bool {
        self.finite
            .as_ref()
            .map(|f| {
                f.op_data
                    .iter()
                    .any(|d| d.kind == SymmetryOpKind::Inversion)
            })
            .unwrap_or(false)
    }

    /// Order of the principal rotation axis. Returns 0 for linear groups.
    pub fn principal_axis_order(&self) -> u32 {
        match self.symbol {
            SchoenfliesSymbol::Cn(n)
            | SchoenfliesSymbol::Cnh(n)
            | SchoenfliesSymbol::Cnv(n)
            | SchoenfliesSymbol::Sn(n)
            | SchoenfliesSymbol::Dn(n)
            | SchoenfliesSymbol::Dnh(n)
            | SchoenfliesSymbol::Dnd(n) => n,
            SchoenfliesSymbol::Ci | SchoenfliesSymbol::Cs => 1,
            SchoenfliesSymbol::T | SchoenfliesSymbol::Td | SchoenfliesSymbol::Th => 3,
            SchoenfliesSymbol::O | SchoenfliesSymbol::Oh => 4,
            SchoenfliesSymbol::I | SchoenfliesSymbol::Ih => 5,
            SchoenfliesSymbol::K | SchoenfliesSymbol::Kh => 0,
            SchoenfliesSymbol::Coov | SchoenfliesSymbol::Dooh => 0,
        }
    }

    pub fn order(&self) -> usize {
        self.finite.as_ref().map(|f| f.order).unwrap_or(0)
    }

    pub fn class_sizes(&self) -> &[i32] {
        self.finite
            .as_ref()
            .map(|f| f.class_sizes.as_slice())
            .unwrap_or(&[])
    }

    /// All symmetry operations in positional order (one handle per group element).
    pub fn ops(&'static self) -> Vec<SymmetryOp> {
        match &self.finite {
            Some(f) => (0..f.op_data.len())
                .map(|i| SymmetryOp {
                    group: self,
                    index: i,
                })
                .collect(),
            None => Vec::new(),
        }
    }

    /// One representative per conjugacy class.
    /// `class_reps()[i]` corresponds to `class_sizes()[i]` and `Irrep::characters()[i]`.
    pub fn class_reps(&'static self) -> Vec<SymmetryOp> {
        match &self.finite {
            Some(f) => f
                .class_rep_indices
                .iter()
                .map(|&i| SymmetryOp {
                    group: self,
                    index: i,
                })
                .collect(),
            None => Vec::new(),
        }
    }

    /// Product of two ops in the group: returns the op `a · b`.
    pub fn multiply(&'static self, a: SymmetryOp, b: SymmetryOp) -> SymmetryOp {
        debug_assert!(ptr::eq(a.group, self));
        debug_assert!(ptr::eq(b.group, self));
        let finite = self.finite.as_ref().expect("multiply on linear group");
        SymmetryOp {
            group: self,
            index: finite.mul_table[a.index][b.index],
        }
    }

    pub fn irreps(&'static self) -> Vec<Irrep> {
        self.irreps
            .iter()
            .map(|d| Irrep {
                data: d,
                group: self,
            })
            .collect()
    }

    /// The totally symmetric irrep (A1, Ag, Σ+, etc.).
    pub fn totally_symmetric_irrep(&'static self) -> Irrep {
        Irrep {
            data: &self.irreps[0],
            group: self,
        }
    }

    pub fn irrep(&'static self, symbol: &str) -> Option<Irrep> {
        self.irreps
            .iter()
            .find(|d| d.symbol == symbol)
            .or_else(|| {
                // Use standard Mulliken notation (E in C3v) instead of libmsym (E1)
                let b = symbol.as_bytes();
                if b.len() >= 2
                    && b[0].is_ascii_uppercase()
                    && b[1] == b'1'
                    && (b.len() == 2 || !b[2].is_ascii_digit())
                {
                    let stripped = format!("{}{}", b[0] as char, &symbol[2..]);
                    self.irreps.iter().find(|d| d.symbol == stripped)
                } else {
                    None
                }
            })
            .map(|d| Irrep {
                data: d,
                group: self,
            })
    }

    /// Decompose the direct product a ⊗ b into irreps with multiplicities.
    pub fn direct_product(&'static self, a: Irrep, b: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        debug_assert!(ptr::eq(b.group, self));

        if self.finite.is_some() {
            let product_chars: Vec<f64> = a
                .characters()
                .iter()
                .zip(b.characters())
                .map(|(ca, cb)| ca * cb)
                .collect();
            self.reduce(&product_chars)
                .expect("valid product characters")
        } else {
            linear::direct_product(self, a, b)
        }
    }

    /// Symmetric square \[a²\]: decompose the symmetric part of a ⊗ a.
    pub fn symmetric_square(&'static self, a: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        let Some(f) = &self.finite else {
            return linear::symmetric_square(self, a);
        };
        let chars = a.characters();
        let sym_chars: Vec<f64> = (0..f.class_sizes.len())
            .map(|c| 0.5 * (chars[c] * chars[c] + chars[f.r_squared_class[c]]))
            .collect();
        self.reduce(&sym_chars)
            .expect("valid symmetric square characters")
    }

    /// Antisymmetric square {a²}: decompose the antisymmetric part of a ⊗ a.
    pub fn antisymmetric_square(&'static self, a: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        let Some(f) = &self.finite else {
            return linear::antisymmetric_square(self, a);
        };
        let chars = a.characters();
        let anti_chars: Vec<f64> = (0..f.class_sizes.len())
            .map(|c| 0.5 * (chars[c] * chars[c] - chars[f.r_squared_class[c]]))
            .collect();
        self.reduce(&anti_chars)
            .expect("valid antisymmetric square characters")
    }

    /// Reduce a representation (given by its class characters) into irreps with multiplicities.
    pub fn reduce(&'static self, characters: &[f64]) -> Result<Vec<(Irrep, u32)>, ReductionError> {
        let Some(f) = &self.finite else {
            return Err(ReductionError::InfiniteGroup);
        };
        let d = self.irreps.len();
        if characters.len() != d {
            return Err(ReductionError::DimensionMismatch {
                expected: d,
                got: characters.len(),
            });
        }
        let h = f.order as f64;

        let mut result = Vec::new();
        for ir_data in &self.irreps {
            let n: f64 = (0..d)
                .map(|c| f.class_sizes[c] as f64 * ir_data.characters()[c] * characters[c])
                .sum::<f64>()
                / h;
            if (n - n.round()).abs() > REDUCTION_INTEGRALITY {
                return Err(ReductionError::NonIntegralMultiplicity {
                    irrep: ir_data.symbol.clone(),
                    value: n,
                });
            }
            let n_rounded = n.round() as u32;
            if n_rounded > 0 {
                result.push((
                    Irrep {
                        data: ir_data,
                        group: self,
                    },
                    n_rounded,
                ));
            }
        }
        Ok(result)
    }

    /// Irreps spanned by the translational degrees of freedom (x, y, z).
    pub fn translation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.finite {
            Some(f) => self
                .reduce(&f.translation_chars)
                .expect("valid translation characters"),
            None => linear::translation_irreps(self),
        }
    }

    /// Irreps spanned by the rotational degrees of freedom (Rx, Ry, Rz).
    pub fn rotation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.finite {
            Some(f) => self
                .reduce(&f.rotation_chars)
                .expect("valid rotation characters"),
            None => linear::rotation_irreps(self),
        }
    }

    /// Irreps of the symmetric square of the vector representation.
    fn quadratic_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.finite {
            Some(f) => self
                .reduce(&f.quadratic_chars)
                .expect("valid quadratic characters"),
            None => linear::quadratic_irreps(self),
        }
    }

    pub fn electric_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.translation_irreps()
            .iter()
            .any(|(gamma_t, _)| self.contains_totally_symmetric(initial, *gamma_t, final_))
    }

    pub fn magnetic_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.rotation_irreps()
            .iter()
            .any(|(gamma_r, _)| self.contains_totally_symmetric(initial, *gamma_r, final_))
    }

    pub fn raman_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.quadratic_irreps()
            .iter()
            .any(|(gamma_q, _)| self.contains_totally_symmetric(initial, *gamma_q, final_))
    }

    pub fn electric_quadrupole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        self.raman_allowed(initial, final_)
    }

    pub fn contains_totally_symmetric(&'static self, a: Irrep, b: Irrep, c: Irrep) -> bool {
        debug_assert!(ptr::eq(a.group, self));
        debug_assert!(ptr::eq(b.group, self));
        debug_assert!(ptr::eq(c.group, self));

        let Some(f) = &self.finite else {
            return linear::contains_totally_symmetric(self, a, b, c);
        };
        let d = self.irreps.len();
        let h = f.order as f64;
        let n: f64 = (0..d)
            .map(|cls| {
                f.class_sizes[cls] as f64
                    * a.characters()[cls]
                    * b.characters()[cls]
                    * c.characters()[cls]
            })
            .sum::<f64>()
            / h;
        n.round() as u32 > 0
    }

    pub fn character_table(&'static self) -> CharacterTableDisplay {
        CharacterTableDisplay { group: self }
    }

    pub(crate) fn finite_data(&self) -> Option<&FiniteGroupData> {
        self.finite.as_ref()
    }
}

fn build_finite_inner(
    ctx: ffi::msym_context,
    label: SchoenfliesSymbol,
    c_name: &CStr,
) -> Result<PointGroup, MsymError> {
    error::check(unsafe { ffi::msymSetPointGroupByName(ctx, c_name.as_ptr()) })?;

    let mut len: c_int = 0;
    let mut sops_ptr: *const ffi::msym_symmetry_operation_t = ptr::null();
    error::check(unsafe { ffi::msymGetSymmetryOperations(ctx, &mut len, &mut sops_ptr) })?;
    let sops_slice = unsafe { slice::from_raw_parts(sops_ptr, len as usize) };

    let op_data: Vec<SymmetryOpData> = sops_slice
        .iter()
        .map(|sop| SymmetryOpData {
            kind: sop.type_.into(),
            order: sop.order,
            power: sop.power,
            orientation: sop.orientation.into(),
            class: sop.cla as usize,
        })
        .collect();
    let matrices: Vec<Matrix3<f64>> = sops_slice.iter().map(compute_op_matrix).collect();

    let order = op_data.len();
    let n_classes = op_data.iter().map(|d| d.class).max().unwrap_or(0) + 1;

    let mut classes: Vec<Vec<usize>> = vec![Vec::new(); n_classes];
    for (i, d) in op_data.iter().enumerate() {
        classes[d.class].push(i);
    }

    let class_rep_indices: Vec<usize> = classes.iter().map(|c| c[0]).collect();
    let class_sizes: Vec<i32> = classes.iter().map(|c| c.len() as i32).collect();

    let mul_table = build_mul_table(&matrices);

    let r_squared_class: Vec<usize> = class_rep_indices
        .iter()
        .map(|&rep| op_data[mul_table[rep][rep]].class)
        .collect();

    let translation_chars: Vec<f64> = class_rep_indices
        .iter()
        .map(|&i| matrices[i].trace())
        .collect();
    let rotation_chars: Vec<f64> = class_rep_indices
        .iter()
        .map(|&i| matrices[i].determinant() * matrices[i].trace())
        .collect();
    let quadratic_chars: Vec<f64> = class_rep_indices
        .iter()
        .map(|&i| {
            let m = &matrices[i];
            let tr = m.trace();
            let tr2 = (m * m).trace();
            (tr * tr + tr2) / 2.0
        })
        .collect();

    let mut ct_ptr: *const ffi::msym_character_table_t = ptr::null();
    error::check(unsafe { ffi::msymGetCharacterTable(ctx, &mut ct_ptr) })?;
    let mut irrep_data = unsafe { extract_irrep_data(&*ct_ptr)? };
    normalize_irrep_symbols(&mut irrep_data);

    Ok(PointGroup {
        symbol: label,
        irreps: irrep_data,
        finite: Some(FiniteGroupData {
            order,
            op_data,
            classes,
            class_rep_indices,
            class_sizes,
            mul_table,
            r_squared_class,
            translation_chars,
            rotation_chars,
            quadratic_chars,
        }),
    })
}

pub(crate) fn compute_op_matrix(sop: &ffi::msym_symmetry_operation_t) -> Matrix3<f64> {
    match sop.type_ {
        ffi::MSYM_SYMMETRY_OPERATION_TYPE_IDENTITY => Matrix3::identity(),
        ffi::MSYM_SYMMETRY_OPERATION_TYPE_INVERSION => -Matrix3::identity(),
        ffi::MSYM_SYMMETRY_OPERATION_TYPE_REFLECTION => reflection_matrix(sop.v),
        ffi::MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION => {
            let theta = 2.0 * PI * sop.power as f64 / sop.order as f64;
            rotation_matrix(sop.v, theta)
        }
        ffi::MSYM_SYMMETRY_OPERATION_TYPE_IMPROPER_ROTATION => {
            let theta = 2.0 * PI * sop.power as f64 / sop.order as f64;
            rotation_matrix(sop.v, theta) * reflection_matrix(sop.v)
        }
        _ => Matrix3::identity(),
    }
}

fn rotation_matrix(axis: [f64; 3], theta: f64) -> Matrix3<f64> {
    let axis = Vector3::from(axis).normalize();
    let (s, c) = theta.sin_cos();
    let one_c = 1.0 - c;
    let (x, y, z) = (axis.x, axis.y, axis.z);
    Matrix3::new(
        c + x * x * one_c,
        x * y * one_c - z * s,
        x * z * one_c + y * s,
        y * x * one_c + z * s,
        c + y * y * one_c,
        y * z * one_c - x * s,
        z * x * one_c - y * s,
        z * y * one_c + x * s,
        c + z * z * one_c,
    )
}

fn reflection_matrix(normal: [f64; 3]) -> Matrix3<f64> {
    let n = Vector3::from(normal).normalize();
    Matrix3::identity() - 2.0 * n * n.transpose()
}

fn build_mul_table(matrices: &[Matrix3<f64>]) -> Vec<Vec<usize>> {
    let n = matrices.len();
    let mut table = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in 0..n {
            let product = matrices[i] * matrices[j];
            let (best, _) = matrices
                .iter()
                .enumerate()
                .map(|(k, m)| {
                    let diff = product - m;
                    let dist: f64 = diff.iter().map(|x| x * x).sum();
                    (k, dist)
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap();
            table[i][j] = best;
        }
    }
    table
}

unsafe fn extract_irrep_data(
    ct: &ffi::msym_character_table_t,
) -> Result<Vec<IrrepData>, MsymError> {
    let d = ct.d as usize;
    if d == 0 {
        return Err(MsymError {
            code: MSYM_INVALID_CHARACTER_TABLE,
            message: "character table has zero dimension".into(),
        });
    }

    let species = unsafe { slice::from_raw_parts(ct.s, d) };
    let table: &[f64] = unsafe { slice::from_raw_parts(ct.table as *const f64, d * d) };

    let mut irrep_data = Vec::with_capacity(d);
    for i in 0..d {
        let symbol = unsafe { CStr::from_ptr(species[i].name.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        let characters: Vec<f64> = (0..d).map(|c| table[i * d + c]).collect();
        irrep_data.push(IrrepData {
            symbol,
            dimension: species[i].d,
            reduction: ReductionData::Finite { characters },
        });
    }
    Ok(irrep_data)
}

/// Strip redundant "1" subscripts from irrep symbols to match standard Mulliken notation.
/// libmsym always numbers E/A irreps (E1, A1) even when only one exists.
fn normalize_irrep_symbols(data: &mut [IrrepData]) {
    let mut max_sub: HashMap<u8, u32> = HashMap::new();
    for d in data.iter() {
        let b = d.symbol.as_bytes();
        if b.len() >= 2 && b[0].is_ascii_uppercase() && b[1].is_ascii_digit() {
            let end = b[1..]
                .iter()
                .position(|c| !c.is_ascii_digit())
                .map_or(b.len(), |i| i + 1);
            if let Ok(n) = d.symbol[1..end].parse::<u32>() {
                let entry = max_sub.entry(b[0]).or_insert(0);
                *entry = (*entry).max(n);
            }
        }
    }
    for d in data.iter_mut() {
        let b = d.symbol.as_bytes();
        if b.len() >= 2
            && b[1] == b'1'
            && (b.len() == 2 || !b[2].is_ascii_digit())
            && max_sub.get(&b[0]) == Some(&1)
        {
            d.symbol.remove(1);
        }
    }
}

pub struct CharacterTableDisplay {
    group: &'static PointGroup,
}

impl fmt::Display for CharacterTableDisplay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let g = self.group;
        let class_reps = g.class_reps();
        if class_reps.is_empty() {
            return write!(f, "{} (linear group, no finite character table)", g.symbol);
        }
        let class_sizes = g.class_sizes();
        let irreps = g.irreps();

        let col_headers: Vec<String> = class_reps
            .iter()
            .zip(class_sizes)
            .map(|(op, &size)| {
                if size == 1 {
                    op.to_string()
                } else {
                    format!("{size}{op}")
                }
            })
            .collect();

        let irrep_col_width = irreps.iter().map(|ir| ir.symbol().len()).max().unwrap_or(0);
        let label_width = g.symbol.to_string().len().max(irrep_col_width);

        let col_widths: Vec<usize> = col_headers
            .iter()
            .enumerate()
            .map(|(c, header)| {
                let max_char = irreps
                    .iter()
                    .map(|ir| format_character(ir.characters()[c]).len())
                    .max()
                    .unwrap_or(0);
                header.len().max(max_char)
            })
            .collect();

        write!(f, "{:width$} │", g.symbol, width = label_width)?;
        for (c, header) in col_headers.iter().enumerate() {
            write!(f, " {:>width$}", header, width = col_widths[c])?;
        }
        writeln!(f)?;

        write!(f, "{:─>width$}─┼", "", width = label_width)?;
        for &w in &col_widths {
            write!(f, "─{:─>width$}", "", width = w)?;
        }
        writeln!(f)?;

        for ir in &irreps {
            write!(f, "{:width$} │", ir.symbol(), width = label_width)?;
            for (c, &w) in col_widths.iter().enumerate() {
                write!(
                    f,
                    " {:>width$}",
                    format_character(ir.characters()[c]),
                    width = w
                )?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

fn format_character(value: f64) -> String {
    let rounded = value.round();
    if (value - rounded).abs() < CHARACTER_DISPLAY_ROUNDING {
        format!("{}", rounded as i64)
    } else {
        format!("{value:.4}")
    }
}

impl fmt::Display for PointGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_point_group_c1() {
        let g = group!(C1);
        assert_eq!(g.to_string(), "C1");
        assert_eq!(g.order(), 1);
        assert_eq!(g.ops().len(), 1);
        assert_eq!(g.symbol(), SchoenfliesSymbol::Cn(1));
        assert_eq!(g.irreps().len(), 1);
        assert_eq!(g.irreps()[0].symbol(), "A");
    }

    #[rstest]
    #[case::c2v("C2v", SchoenfliesSymbol::Cnv(2), 4)]
    #[case::td("Td", SchoenfliesSymbol::Td, 24)]
    #[case::oh("Oh", SchoenfliesSymbol::Oh, 48)]
    #[case::c3v("C3v", SchoenfliesSymbol::Cnv(3), 6)]
    #[case::d2h("D2h", SchoenfliesSymbol::Dnh(2), 8)]
    fn test_point_group_parse(
        #[case] name: &str,
        #[case] expected_label: SchoenfliesSymbol,
        #[case] expected_order: usize,
    ) {
        let g = PointGroup::parse(name).unwrap();
        assert_eq!(g.to_string(), name);
        assert_eq!(g.symbol(), expected_label);
        assert_eq!(g.order(), expected_order);
        assert_eq!(g.ops().len(), expected_order);
    }

    #[rstest]
    fn test_point_group_pointer_identity() {
        let a = group!(C2v);
        let b = PointGroup::parse("C2v").unwrap();
        assert!(ptr::eq(a, b));
    }

    #[rstest]
    fn test_group_macro() {
        assert_eq!(group!(Td).to_string(), "Td");
        assert_eq!(group!(Oh).to_string(), "Oh");
        assert_eq!(group!(Ih).to_string(), "Ih");
        assert_eq!(group!(D2h).to_string(), "D2h");
        assert_eq!(group!(Coov).to_string(), "C∞v");
        assert_eq!(group!(Dooh).to_string(), "D∞h");
    }

    #[rstest]
    #[case::c1("C1", true)]
    #[case::ci("Ci", true)]
    #[case::cs("Cs", true)]
    #[case::c2("C2", true)]
    #[case::c3("C3", true)]
    #[case::c2v("C2v", true)]
    #[case::c2h("C2h", true)]
    #[case::c3h("C3h", true)]
    #[case::d2("D2", true)]
    #[case::d2h("D2h", true)]
    #[case::s4("S4", true)]
    #[case::s6("S6", true)]
    #[case::c3v("C3v", false)]
    #[case::d3("D3", false)]
    #[case::d3h("D3h", false)]
    #[case::d2d("D2d", false)]
    #[case::td("Td", false)]
    #[case::oh("Oh", false)]
    #[case::ih("Ih", false)]
    fn test_point_group_is_abelian(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_abelian(), expected);
    }

    #[rstest]
    #[case::c1("C1", true)]
    #[case::c3("C3", true)]
    #[case::c3h("C3h", true)]
    #[case::s4("S4", true)]
    #[case::ci("Ci", true)]
    #[case::cs("Cs", true)]
    #[case::c2v("C2v", false)]
    #[case::d2("D2", false)]
    #[case::td("Td", false)]
    fn test_point_group_is_cyclic(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_cyclic(), expected);
    }

    #[rstest]
    #[case::c1("C1", false)]
    #[case::d2("D2", true)]
    #[case::d2h("D2h", true)]
    #[case::d3d("D3d", true)]
    #[case::dooh("Dooh", true)]
    #[case::td("Td", false)]
    fn test_point_group_is_dihedral(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_dihedral(), expected);
    }

    #[rstest]
    #[case::t("T", true)]
    #[case::td("Td", true)]
    #[case::th("Th", true)]
    #[case::o("O", true)]
    #[case::oh("Oh", true)]
    #[case::ih("Ih", false)]
    #[case::c2v("C2v", false)]
    #[case::d6h("D6h", false)]
    fn test_point_group_is_cubic(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_cubic(), expected);
    }

    #[rstest]
    #[case::t("T", false)]
    #[case::td("Td", false)]
    #[case::th("Th", false)]
    #[case::i("I", true)]
    #[case::ih("Ih", true)]
    #[case::c2v("C2v", false)]
    #[case::d6h("D6h", false)]
    fn test_point_group_is_icosahedral(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_icosahedral(), expected);
    }

    #[rstest]
    #[case::c1("C1", false)]
    #[case::coov("Coov", true)]
    #[case::dooh("Dooh", true)]
    fn test_point_group_is_linear(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_linear(), expected);
    }

    #[rstest]
    #[case::c3("C3", true)]
    #[case::c4("C4", true)]
    #[case::c5("C5", true)]
    #[case::c3h("C3h", true)]
    #[case::s4("S4", true)]
    #[case::s6("S6", true)]
    #[case::t("T", true)]
    #[case::th("Th", true)]
    #[case::c1("C1", false)]
    #[case::c2("C2", false)]
    #[case::c2v("C2v", false)]
    #[case::td("Td", false)]
    #[case::oh("Oh", false)]
    #[case::c3v("C3v", false)]
    #[case::d3h("D3h", false)]
    fn test_point_group_has_complex_irreps(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.has_complex_irreps(), expected);
    }

    #[rstest]
    #[case::ci("Ci", true)]
    #[case::c2h("C2h", true)]
    #[case::d2h("D2h", true)]
    #[case::oh("Oh", true)]
    #[case::ih("Ih", true)]
    #[case::th("Th", true)]
    #[case::c2v("C2v", false)]
    #[case::td("Td", false)]
    #[case::c3v("C3v", false)]
    #[case::d3("D3", false)]
    #[case::c1("C1", false)]
    fn test_point_group_has_inversion(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.has_inversion(), expected);
    }

    #[rstest]
    #[case::c1("C1", true)]
    #[case::c2("C2", true)]
    #[case::c3("C3", true)]
    #[case::d2("D2", true)]
    #[case::d3("D3", true)]
    #[case::t("T", true)]
    #[case::o("O", true)]
    #[case::i("I", true)]
    #[case::cs("Cs", false)]
    #[case::ci("Ci", false)]
    #[case::c2v("C2v", false)]
    #[case::c2h("C2h", false)]
    #[case::d2h("D2h", false)]
    #[case::d3h("D3h", false)]
    #[case::td("Td", false)]
    #[case::oh("Oh", false)]
    #[case::ih("Ih", false)]
    #[case::s4("S4", false)]
    fn test_point_group_is_chiral(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.is_chiral(), expected);
    }

    #[rstest]
    #[case::c1("C1", 1)]
    #[case::cs("Cs", 1)]
    #[case::ci("Ci", 1)]
    #[case::c2v("C2v", 2)]
    #[case::c3v("C3v", 3)]
    #[case::d6h("D6h", 6)]
    #[case::s4("S4", 4)]
    #[case::td("Td", 3)]
    #[case::oh("Oh", 4)]
    #[case::ih("Ih", 5)]
    fn test_point_group_principal_axis_order(#[case] group: &str, #[case] expected: u32) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.principal_axis_order(), expected);
    }

    #[rstest]
    fn test_point_group_irreps() {
        let g = group!(C2v);
        let irreps = g.irreps();
        assert_eq!(irreps.len(), 4);
        let symbols: Vec<&str> = irreps.iter().map(|ir| ir.symbol()).collect();
        assert!(symbols.contains(&"A1"));
        assert!(symbols.contains(&"A2"));
        assert!(symbols.contains(&"B1"));
        assert!(symbols.contains(&"B2"));
        for ir in &irreps {
            assert_eq!(ir.dimension(), 1);
            assert!(ptr::eq(ir.group(), g));
        }
    }

    #[rstest]
    fn test_point_group_irrep() {
        let g = group!(C2v);
        let a1 = g.irrep("A1").unwrap();
        assert_eq!(a1.symbol(), "A1");
        assert!(g.irrep("nonexistent").is_none());
    }

    #[rstest]
    #[case::c2v_b1_b2("C2v", "B1", "B2", &[("A2", 1)])]
    #[case::c2v_a1_b1("C2v", "A1", "B1", &[("B1", 1)])]
    #[case::td_e_t2("Td", "E", "T2", &[("T1", 1), ("T2", 1)])]
    #[case::oh_eg_t1u("Oh", "Eg", "T1u", &[("T1u", 1), ("T2u", 1)])]
    #[case::c3_a_e("C3", "A", "E", &[("E", 2)])]
    #[case::c3_e_e("C3", "E", "E", &[("A", 2), ("E", 2)])]
    #[case::c4_b_e("C4", "B", "E", &[("E", 2)])]
    #[case::c4_e_e("C4", "E", "E", &[("A", 2), ("B", 2)])]
    #[case::c3h_app_ep("C3h", "A''", "E'", &[("E''", 2)])]
    #[case::c3h_ep_ep("C3h", "E'", "E'", &[("A'", 2), ("E'", 2)])]
    #[case::c3h_ep_epp("C3h", "E'", "E''", &[("A''", 2), ("E''", 2)])]
    fn test_point_group_direct_product(
        #[case] group: &str,
        #[case] a: &str,
        #[case] b: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir_a = g.irrep(a).unwrap();
        let ir_b = g.irrep(b).unwrap();
        let product = g.direct_product(ir_a, ir_b);
        let actual: Vec<(&str, u32)> = product.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::c2v_b1("C2v", "B1", &[("A1", 1)], &[])]
    #[case::td_e("Td", "E", &[("A1", 1), ("E", 1)], &[("A2", 1)])]
    #[case::td_t2("Td", "T2", &[("A1", 1), ("E", 1), ("T2", 1)], &[("T1", 1)])]
    #[case::oh_t1u("Oh", "T1u", &[("A1g", 1), ("Eg", 1), ("T2g", 1)], &[("T1g", 1)])]
    #[case::ih_hg("Ih", "Hg", &[("Ag", 1), ("Gg", 1), ("Hg", 2)], &[("T1g", 1), ("T2g", 1), ("Gg", 1)])]
    fn test_point_group_symmetric_antisymmetric_square(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected_sym: &[(&str, u32)],
        #[case] expected_anti: &[(&str, u32)],
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();

        let sym = g.symmetric_square(ir);
        let sym_actual: Vec<(&str, u32)> = sym.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(sym_actual, expected_sym);

        let anti = g.antisymmetric_square(ir);
        let anti_actual: Vec<(&str, u32)> = anti.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(anti_actual, expected_anti);
    }

    #[rstest]
    fn test_point_group_contains_totally_symmetric() {
        let g = group!(C2v);
        let a1 = g.irrep("A1").unwrap();
        let b1 = g.irrep("B1").unwrap();
        let b2 = g.irrep("B2").unwrap();
        assert!(g.contains_totally_symmetric(a1, b1, b1));
        assert!(!g.contains_totally_symmetric(a1, b1, b2));
    }

    #[rstest]
    fn test_irrep_equality() {
        let g = group!(C2v);
        let a1_first = g.irrep("A1").unwrap();
        let a1_second = g.irrep("A1").unwrap();
        let b1 = g.irrep("B1").unwrap();
        assert_eq!(a1_first, a1_second);
        assert_ne!(a1_first, b1);
    }

    #[rstest]
    #[case::c2v_regular("C2v", &[4.0, 0.0, 0.0, 0.0], &[("A1", 1), ("A2", 1), ("B1", 1), ("B2", 1)])]
    #[case::c2v_totally_symmetric("C2v", &[1.0, 1.0, 1.0, 1.0], &[("A1", 1)])]
    #[case::c1_septuple("C1",  &[7.0], &[("A", 7)])]
    fn test_point_group_reduce(
        #[case] group: &str,
        #[case] characters: &[f64],
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::parse(group).unwrap();
        let result = g.reduce(characters).unwrap();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::c2v("C2v", &[("A1", 2), ("B2", 3)])]
    #[case::td("Td",  &[("A1", 1), ("E", 1), ("T2", 3)])]
    #[case::oh("Oh",  &[("Eg", 2), ("T1u", 1)])]
    fn test_point_group_reduce_roundtrip(#[case] group: &str, #[case] composition: &[(&str, u32)]) {
        let g = PointGroup::parse(group).unwrap();
        let n_classes = g.class_sizes().len();
        let mut chars = vec![0.0; n_classes];
        for &(sym, mult) in composition {
            let ir = g.irrep(sym).unwrap();
            for (c, ch) in ir.characters().iter().enumerate() {
                chars[c] += mult as f64 * ch;
            }
        }
        let result = g.reduce(&chars).unwrap();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        let expected: Vec<(&str, u32)> = composition.to_vec();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::c2v("C2v", &[("A1", 1), ("B1", 1), ("B2", 1)])]
    #[case::td("Td",  &[("T2", 1)])]
    #[case::c1("C1",  &[("A", 3)])]
    fn test_point_group_translation_irreps(#[case] group: &str, #[case] expected: &[(&str, u32)]) {
        let g = PointGroup::parse(group).unwrap();
        let result = g.translation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::c2v("C2v", &[("A2", 1), ("B1", 1), ("B2", 1)])]
    #[case::td("Td",  &[("T1", 1)])]
    #[case::c1("C1",  &[("A", 3)])]
    fn test_point_group_rotation_irreps(#[case] group: &str, #[case] expected: &[(&str, u32)]) {
        let g = PointGroup::parse(group).unwrap();
        let result = g.rotation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::identity(SymmetryOpKind::Identity, 1, 1, SymmetryOpOrientation::None, "E")]
    #[case::inversion(SymmetryOpKind::Inversion, 1, 1, SymmetryOpOrientation::None, "i")]
    #[case::sigma_h(SymmetryOpKind::Reflection, 1, 1, SymmetryOpOrientation::Horizontal, "σh")]
    #[case::sigma_v(SymmetryOpKind::Reflection, 1, 1, SymmetryOpOrientation::Vertical, "σv")]
    #[case::sigma_d(SymmetryOpKind::Reflection, 1, 1, SymmetryOpOrientation::Dihedral, "σd")]
    #[case::sigma_bare(SymmetryOpKind::Reflection, 1, 1, SymmetryOpOrientation::None, "σ")]
    #[case::c3(SymmetryOpKind::ProperRotation, 3, 1, SymmetryOpOrientation::None, "C3")]
    #[case::c3_squared(SymmetryOpKind::ProperRotation, 3, 2, SymmetryOpOrientation::None, "C3²")]
    #[case::c6_fifth(SymmetryOpKind::ProperRotation, 6, 5, SymmetryOpOrientation::None, "C6⁵")]
    #[case::s4(SymmetryOpKind::ImproperRotation, 4, 1, SymmetryOpOrientation::None, "S4")]
    #[case::s4_cubed(SymmetryOpKind::ImproperRotation, 4, 3, SymmetryOpOrientation::None, "S4³")]
    #[case::s10_seventh(SymmetryOpKind::ImproperRotation, 10, 7, SymmetryOpOrientation::None, "S10⁷")]
    #[case::c2_prime(SymmetryOpKind::ProperRotation, 2, 1, SymmetryOpOrientation::Vertical, "C2'")]
    #[case::c2_double_prime(SymmetryOpKind::ProperRotation, 2, 1, SymmetryOpOrientation::Dihedral, "C2''")]
    #[case::c3_prime_squared(SymmetryOpKind::ProperRotation, 3, 2, SymmetryOpOrientation::Vertical, "C3'²")]
    #[case::s4_prime(SymmetryOpKind::ImproperRotation, 4, 1, SymmetryOpOrientation::Vertical, "S4'")]
    fn test_op_data_display(
        #[case] kind: SymmetryOpKind,
        #[case] order: i32,
        #[case] power: i32,
        #[case] orientation: SymmetryOpOrientation,
        #[case] expected: &str,
    ) {
        let d = SymmetryOpData {
            kind,
            order,
            power,
            orientation,
            class: 0,
        };
        assert_eq!(d.to_string(), expected);
    }

    #[rstest]
    fn test_character_table_display_c2v() {
        let g = group!(C2v);
        let table = g.character_table().to_string();
        let expected = "\
C2v │ E C2  σv  σd
────┼─────────────
A1  │ 1  1   1   1
A2  │ 1  1  -1  -1
B1  │ 1 -1   1  -1
B2  │ 1 -1  -1   1
";
        assert_eq!(table, expected);
    }

    #[rstest]
    fn test_character_table_display_c3v() {
        let g = PointGroup::parse("C3v").unwrap();
        let table = g.character_table().to_string();
        let expected = "\
C3v │ E 2C3  3σv
────┼───────────
A1  │ 1   1    1
A2  │ 1   1   -1
E   │ 2  -1    0
";
        assert_eq!(table, expected);
    }

    #[rstest]
    #[case::c1("C1", "A")]
    #[case::c2v("C2v", "A1")]
    #[case::oh("Oh", "A1g")]
    #[case::td("Td", "A1")]
    fn test_point_group_totally_symmetric_irrep(#[case] group: &str, #[case] expected: &str) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(g.totally_symmetric_irrep().symbol(), expected);
    }

    #[rstest]
    fn test_irrep_gerade_finite() {
        let g = group!(Oh);
        assert_eq!(g.irrep("A1g").unwrap().gerade(), Some(true));
        assert_eq!(g.irrep("Eg").unwrap().gerade(), Some(true));
        assert_eq!(g.irrep("T1u").unwrap().gerade(), Some(false));
        assert_eq!(g.irrep("T2u").unwrap().gerade(), Some(false));
    }

    #[rstest]
    fn test_irrep_gerade_no_inversion() {
        let g = group!(C2v);
        for ir in g.irreps() {
            assert_eq!(ir.gerade(), None);
        }
    }

    #[rstest]
    fn test_symmetry_op_is_proper() {
        let g = group!(Oh);
        for op in g.ops() {
            assert_eq!(
                op.is_proper(),
                matches!(
                    op.kind(),
                    SymmetryOpKind::Identity | SymmetryOpKind::ProperRotation
                ),
            );
        }
    }

    #[rstest]
    fn test_symmetry_op_equality_via_handle_identity() {
        let g = group!(C2v);
        let ops = g.ops();
        // Two reflections in C2v share descriptive tuple but must compare unequal.
        let reflections: Vec<SymmetryOp> = ops
            .iter()
            .filter(|op| op.kind() == SymmetryOpKind::Reflection)
            .copied()
            .collect();
        assert_eq!(reflections.len(), 2);
        assert_ne!(reflections[0], reflections[1]);
        assert_ne!(reflections[0].index(), reflections[1].index());
    }

    #[rstest]
    fn test_symmetry_op_cross_group_inequality() {
        let c2v = group!(C2v);
        let d2h = group!(D2h);
        let c2v_inv = c2v
            .ops()
            .iter()
            .find(|op| op.kind() == SymmetryOpKind::Identity)
            .copied();
        let d2h_inv = d2h
            .ops()
            .iter()
            .find(|op| op.kind() == SymmetryOpKind::Identity)
            .copied();
        assert!(c2v_inv.is_some() && d2h_inv.is_some());
        assert_ne!(c2v_inv.unwrap(), d2h_inv.unwrap());
    }

    #[rstest]
    fn test_point_group_multiply_closure() {
        let g = group!(C2v);
        let ops = g.ops();
        // Closure: for every pair (a, b), multiply(a, b) returns an op in the group.
        for a in &ops {
            for b in &ops {
                let c = g.multiply(*a, *b);
                assert!(c.index() < ops.len());
                assert!(ptr::eq(c.group(), g));
            }
        }
    }

    #[rstest]
    #[case::c2v("C2v")]
    #[case::td("Td")]
    #[case::oh("Oh")]
    #[case::d3h("D3h")]
    #[case::c1("C1")]
    fn test_point_group_trans_rot_dimensions(#[case] group: &str) {
        let g = PointGroup::parse(group).unwrap();
        let trans_dim: u32 = g
            .translation_irreps()
            .iter()
            .map(|(ir, n)| ir.dimension() as u32 * n)
            .sum();
        let rot_dim: u32 = g
            .rotation_irreps()
            .iter()
            .map(|(ir, n)| ir.dimension() as u32 * n)
            .sum();
        assert_eq!(trans_dim, 3);
        assert_eq!(rot_dim, 3);
    }

    fn ir_active_symbols(g: &'static PointGroup) -> Vec<String> {
        let a1 = g.irreps()[0];
        g.irreps()
            .into_iter()
            .filter(|ir| g.electric_dipole_allowed(a1, *ir))
            .map(|ir| ir.symbol().to_owned())
            .collect()
    }

    fn raman_active_symbols(g: &'static PointGroup) -> Vec<String> {
        let a1 = g.irreps()[0];
        g.irreps()
            .into_iter()
            .filter(|ir| g.raman_allowed(a1, *ir))
            .map(|ir| ir.symbol().to_owned())
            .collect()
    }

    #[rstest]
    #[case::c1("C1",  &["A"],                  &["A"])]
    #[case::c2v("C2v", &["A1", "B1", "B2"],     &["A1", "A2", "B1", "B2"])]
    #[case::td("Td",  &["T2"],                 &["A1", "E", "T2"])]
    fn test_point_group_ir_active(
        #[case] group: &str,
        #[case] expected_ir: &[&str],
        #[case] expected_raman: &[&str],
    ) {
        let g = PointGroup::parse(group).unwrap();
        assert_eq!(ir_active_symbols(g), expected_ir);
        assert_eq!(raman_active_symbols(g), expected_raman);
    }
}
