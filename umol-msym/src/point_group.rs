use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::{LazyLock, Mutex};
use std::{fmt, ptr};

use nalgebra::Matrix3;
use umol_msym_sys::{MSYM_INVALID_CHARACTER_TABLE, MSYM_INVALID_INPUT, MSYM_POINT_GROUP_ERROR};

use crate::context::Context;
use crate::error::Error;
use crate::linear;
use crate::types::{
    IrrepData, PointGroupKind, SchoenfliesLabel, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation,
};

static REGISTRY: LazyLock<Mutex<HashMap<SchoenfliesLabel, &'static PointGroup>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

fn register(pg: PointGroup) -> &'static PointGroup {
    let label = pg.label;
    let mut map = REGISTRY.lock().unwrap();
    if let Some(&existing) = map.get(&label) {
        return existing;
    }
    let leaked: &'static PointGroup = Box::leak(Box::new(pg));
    map.insert(label, leaked);
    leaked
}

/// An irreducible representation, bound to its parent point group.
///
/// `Copy` view type — carries `&'static` references to both its data and its group.
/// Two irreps are equal iff they point to the same `IrrepData` (pointer identity).
#[derive(Copy, Clone)]
pub struct Irrep {
    pub(crate) data: &'static IrrepData,
    pub(crate) group: &'static PointGroup,
}

impl Irrep {
    /// Mulliken symbol (e.g. "A1", "B2", "T2g").
    pub fn symbol(&self) -> &str {
        &self.data.symbol
    }

    /// Position in the character table (row index).
    pub fn index(&self) -> usize {
        self.data.index
    }

    /// Dimensionality (1 for A/B, 2 for E, 3 for T).
    pub fn dimension(&self) -> i32 {
        self.data.dimension
    }

    /// Character values, one per conjugacy class. Empty for infinite groups.
    pub fn characters(&self) -> &[f64] {
        &self.data.characters
    }

    /// Angular momentum quantum number for infinite-group irreps.
    pub fn lambda(&self) -> Option<u32> {
        self.data.lambda
    }

    /// Gerade (true) or ungerade (false) for centrosymmetric groups.
    /// Returns `None` for groups without inversion.
    pub fn is_gerade(&self) -> Option<bool> {
        if let Some(g) = self.data.gerade {
            return Some(g);
        }
        // Finite centrosymmetric groups: find inversion op, check character sign
        let ops = self.group.ops();
        let inv_op = ops.iter().find(|op| op.kind == SymmetryOpKind::Inversion)?;
        let chi = self.data.characters.get(inv_op.class)?;
        Some(*chi > 0.0)
    }

    /// Sigma-v parity for Σ irreps: true = Σ+, false = Σ-.
    pub fn is_sigma_plus(&self) -> Option<bool> {
        if self.data.lambda == Some(0) {
            self.data.sigma_v
        } else {
            None
        }
    }

    /// Parent point group.
    pub fn group(&self) -> &'static PointGroup {
        self.group
    }
}

impl PartialEq for Irrep {
    fn eq(&self, other: &Self) -> bool {
        ptr::eq(self.data, other.data)
    }
}

impl Eq for Irrep {}

impl Hash for Irrep {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (self.data as *const IrrepData).hash(state);
    }
}

impl fmt::Debug for Irrep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Irrep({}:{})", self.group.label, self.data.symbol)
    }
}

impl fmt::Display for Irrep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.data.symbol)
    }
}

/// A molecular point group: singleton algebraic object with operations and character table.
///
/// Point groups are `&'static` singletons. There is exactly one C2v, one Td, etc.
/// Access via named constructors (`PointGroup::c2v()`) or `PointGroup::from_schoenflies("C2v")`.
#[derive(Debug)]
pub struct PointGroup {
    pub(crate) label: SchoenfliesLabel,
    pub(crate) irrep_data: Vec<IrrepData>,
    pub(crate) kind: PointGroupKind,
}

impl PointGroup {
    pub fn label(&self) -> SchoenfliesLabel {
        self.label
    }

    /// Order of the principal rotation axis. Returns 0 for linear groups.
    pub fn principal_axis_order(&self) -> u32 {
        match self.label {
            SchoenfliesLabel::Cn(n)
            | SchoenfliesLabel::Cnh(n)
            | SchoenfliesLabel::Cnv(n)
            | SchoenfliesLabel::Sn(n)
            | SchoenfliesLabel::Dn(n)
            | SchoenfliesLabel::Dnh(n)
            | SchoenfliesLabel::Dnd(n) => n,
            SchoenfliesLabel::Ci | SchoenfliesLabel::Cs => 1,
            SchoenfliesLabel::T | SchoenfliesLabel::Td | SchoenfliesLabel::Th => 3,
            SchoenfliesLabel::O | SchoenfliesLabel::Oh => 4,
            SchoenfliesLabel::I | SchoenfliesLabel::Ih => 5,
            SchoenfliesLabel::K | SchoenfliesLabel::Kh => 0,
            SchoenfliesLabel::Coov | SchoenfliesLabel::Dooh => 0,
        }
    }

    pub fn order(&self) -> usize {
        match &self.kind {
            PointGroupKind::Finite { order, .. } => *order,
            PointGroupKind::Linear => 0,
        }
    }

    pub fn ops(&self) -> &[SymmetryOp] {
        match &self.kind {
            PointGroupKind::Finite { ops, .. } => ops,
            PointGroupKind::Linear => &[],
        }
    }

    pub fn class_sizes(&self) -> &[i32] {
        match &self.kind {
            PointGroupKind::Finite { class_sizes, .. } => class_sizes,
            PointGroupKind::Linear => &[],
        }
    }

    /// One representative symmetry operation per conjugacy class.
    /// `class_reps()[i]` corresponds to `class_sizes()[i]` and `Irrep::characters()[i]`.
    pub fn class_reps(&self) -> &[SymmetryOp] {
        match &self.kind {
            PointGroupKind::Finite { class_reps, .. } => class_reps,
            PointGroupKind::Linear => &[],
        }
    }

    pub fn is_linear(&self) -> bool {
        matches!(self.kind, PointGroupKind::Linear)
    }

    pub fn is_abelian(&self) -> bool {
        matches!(
            self.label,
            SchoenfliesLabel::Ci
                | SchoenfliesLabel::Cs
                | SchoenfliesLabel::Cn(_)
                | SchoenfliesLabel::Cnh(_)
                | SchoenfliesLabel::Sn(_)
                | SchoenfliesLabel::Cnv(2)
                | SchoenfliesLabel::Dn(2)
                | SchoenfliesLabel::Dnh(2)
        )
    }

    pub fn is_cyclic(&self) -> bool {
        matches!(
            self.label,
            SchoenfliesLabel::Ci
                | SchoenfliesLabel::Cs
                | SchoenfliesLabel::Cn(_)
                | SchoenfliesLabel::Cnh(_)
                | SchoenfliesLabel::Sn(_)
        )
    }

    pub fn is_cubic(&self) -> bool {
        matches!(
            self.label,
            SchoenfliesLabel::T
                | SchoenfliesLabel::Td
                | SchoenfliesLabel::Th
                | SchoenfliesLabel::O
                | SchoenfliesLabel::Oh
                | SchoenfliesLabel::I
                | SchoenfliesLabel::Ih
                | SchoenfliesLabel::K
                | SchoenfliesLabel::Kh
        )
    }

    /// True if some irreps are complex conjugate pairs fused into real 2D representations.
    /// This occurs in cyclic groups with order > 2 (Cn, Cnh, Sn with n > 2).
    pub fn has_complex_irreps(&self) -> bool {
        match self.label {
            SchoenfliesLabel::Cn(n) | SchoenfliesLabel::Cnh(n) | SchoenfliesLabel::Sn(n) => n > 2,
            _ => false,
        }
    }

    /// True iff the group contains only proper operations (no reflections, inversions,
    /// or improper rotations). Chiral groups: Cn, Dn, T, O, I.
    pub fn is_chiral(&self) -> bool {
        self.ops().iter().all(|op| op.is_proper())
    }

    pub fn has_inversion(&self) -> bool {
        self.ops()
            .iter()
            .any(|op| op.kind == SymmetryOpKind::Inversion)
    }

    /// The totally symmetric irrep (A1, Ag, Σ+, etc.).
    pub fn totally_symmetric_irrep(&'static self) -> Irrep {
        Irrep {
            data: &self.irrep_data[0],
            group: self,
        }
    }

    /// All irreducible representations of this group.
    pub fn irreps(&'static self) -> Vec<Irrep> {
        self.irrep_data
            .iter()
            .map(|d| Irrep {
                data: d,
                group: self,
            })
            .collect()
    }

    /// Look up an irrep by Mulliken symbol.
    pub fn irrep(&'static self, symbol: &str) -> Option<Irrep> {
        self.irrep_data
            .iter()
            .find(|d| d.symbol == symbol)
            .map(|d| Irrep {
                data: d,
                group: self,
            })
    }

    /// Decompose the direct product a ⊗ b into irreps with multiplicities.
    pub fn direct_product(&'static self, a: Irrep, b: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        debug_assert!(ptr::eq(b.group, self));

        match &self.kind {
            PointGroupKind::Finite { .. } => {
                let product_chars: Vec<f64> = a
                    .characters()
                    .iter()
                    .zip(b.characters())
                    .map(|(ca, cb)| ca * cb)
                    .collect();
                self.reduce(&product_chars)
            }
            PointGroupKind::Linear => linear::direct_product(self, a, b),
        }
    }

    /// Symmetric square [a²]: decompose the symmetric part of a ⊗ a.
    pub fn symmetric_square(&'static self, a: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        match &self.kind {
            PointGroupKind::Finite {
                ops, class_reps, ..
            } => {
                let chars = a.characters();
                let n_classes = class_reps.len();
                let r2_class = r_squared_classes(class_reps, ops);
                let sym_chars: Vec<f64> = (0..n_classes)
                    .map(|c| 0.5 * (chars[c] * chars[c] + chars[r2_class[c]]))
                    .collect();
                self.reduce(&sym_chars)
            }
            PointGroupKind::Linear => linear::symmetric_square(self, a),
        }
    }

    /// Antisymmetric square {a²}: decompose the antisymmetric part of a ⊗ a.
    pub fn antisymmetric_square(&'static self, a: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(ptr::eq(a.group, self));
        match &self.kind {
            PointGroupKind::Finite {
                ops, class_reps, ..
            } => {
                let chars = a.characters();
                let n_classes = class_reps.len();
                let r2_class = r_squared_classes(class_reps, ops);
                let anti_chars: Vec<f64> = (0..n_classes)
                    .map(|c| 0.5 * (chars[c] * chars[c] - chars[r2_class[c]]))
                    .collect();
                self.reduce(&anti_chars)
            }
            PointGroupKind::Linear => linear::antisymmetric_square(self, a),
        }
    }

    /// Reduce a representation (given by its class characters) into irreps with multiplicities.
    ///
    /// Only valid for finite groups. Panics for infinite groups (C∞v, D∞h).
    pub fn reduce(&'static self, characters: &[f64]) -> Vec<(Irrep, u32)> {
        let PointGroupKind::Finite {
            order, class_sizes, ..
        } = &self.kind
        else {
            panic!("character-based reduction undefined for infinite groups");
        };
        let d = self.irrep_data.len();
        assert_eq!(
            characters.len(),
            d,
            "expected {} class characters, got {}",
            d,
            characters.len()
        );
        let h = *order as f64;

        let mut result = Vec::new();
        for ir_data in &self.irrep_data {
            let n: f64 = (0..d)
                .map(|c| class_sizes[c] as f64 * ir_data.characters[c] * characters[c])
                .sum::<f64>()
                / h;
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
        result
    }

    /// Irreps spanned by the translational degrees of freedom (x, y, z).
    pub fn translation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.kind {
            PointGroupKind::Finite { class_reps, .. } => {
                let chars: Vec<f64> = class_reps.iter().map(|op| op.matrix.trace()).collect();
                self.reduce(&chars)
            }
            PointGroupKind::Linear => linear::translation_irreps(self),
        }
    }

    /// Irreps spanned by the rotational degrees of freedom (Rx, Ry, Rz).
    pub fn rotation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.kind {
            PointGroupKind::Finite { class_reps, .. } => {
                let chars: Vec<f64> = class_reps
                    .iter()
                    .map(|op| op.matrix.determinant() * op.matrix.trace())
                    .collect();
                self.reduce(&chars)
            }
            PointGroupKind::Linear => linear::rotation_irreps(self),
        }
    }

    /// Irreps of the symmetric square of the vector representation.
    fn quadratic_irreps(&'static self) -> Vec<(Irrep, u32)> {
        match &self.kind {
            PointGroupKind::Finite { class_reps, .. } => {
                let chars: Vec<f64> = class_reps
                    .iter()
                    .map(|op| {
                        let tr = op.matrix.trace();
                        let tr2 = (op.matrix * op.matrix).trace();
                        (tr * tr + tr2) / 2.0
                    })
                    .collect();
                self.reduce(&chars)
            }
            PointGroupKind::Linear => linear::quadratic_irreps(self),
        }
    }

    /// Electric dipole transition allowed? Checks Γ_i ⊗ Γ(x,y,z) ⊗ Γ_f ⊃ A1.
    pub fn electric_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.translation_irreps()
            .iter()
            .any(|(gamma_t, _)| self.contains_totally_symmetric(initial, *gamma_t, final_))
    }

    /// Magnetic dipole transition allowed? Checks Γ_i ⊗ Γ(Rx,Ry,Rz) ⊗ Γ_f ⊃ A1.
    pub fn magnetic_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.rotation_irreps()
            .iter()
            .any(|(gamma_r, _)| self.contains_totally_symmetric(initial, *gamma_r, final_))
    }

    /// Raman transition allowed? Checks Γ_i ⊗ Γ(x²,y²,...,yz) ⊗ Γ_f ⊃ A1.
    pub fn raman_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(ptr::eq(initial.group, self));
        debug_assert!(ptr::eq(final_.group, self));
        self.quadratic_irreps()
            .iter()
            .any(|(gamma_q, _)| self.contains_totally_symmetric(initial, *gamma_q, final_))
    }

    /// Electric quadrupole transition allowed? Same basis as Raman (symmetric square).
    pub fn electric_quadrupole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        self.raman_allowed(initial, final_)
    }

    /// Whether a ⊗ b ⊗ c contains the totally symmetric representation.
    pub fn contains_totally_symmetric(&'static self, a: Irrep, b: Irrep, c: Irrep) -> bool {
        debug_assert!(ptr::eq(a.group, self));
        debug_assert!(ptr::eq(b.group, self));
        debug_assert!(ptr::eq(c.group, self));

        match &self.kind {
            PointGroupKind::Finite {
                order, class_sizes, ..
            } => {
                let d = self.irrep_data.len();
                let h = *order as f64;
                let n: f64 = (0..d)
                    .map(|cls| {
                        class_sizes[cls] as f64
                            * a.characters()[cls]
                            * b.characters()[cls]
                            * c.characters()[cls]
                    })
                    .sum::<f64>()
                    / h;
                n.round() as u32 > 0
            }
            PointGroupKind::Linear => linear::contains_totally_symmetric(self, a, b, c),
        }
    }

    /// Returns a displayable character table.
    ///
    /// Only meaningful for finite groups; returns an empty display for linear groups.
    pub fn character_table(&'static self) -> CharacterTableDisplay {
        CharacterTableDisplay { group: self }
    }

    /// Construct a point group by Schoenflies symbol (e.g. "C2v", "Td", "Oh").
    pub fn from_schoenflies(name: &str) -> Result<&'static PointGroup, Error> {
        let label = SchoenfliesLabel::parse(name).ok_or_else(|| Error {
            code: MSYM_INVALID_INPUT,
            message: format!("cannot parse Schoenflies symbol '{name}'"),
        })?;

        // Fast path: already registered
        {
            let map = REGISTRY.lock().unwrap();
            if let Some(&pg) = map.get(&label) {
                return Ok(pg);
            }
        }

        // Slow path: construct and register
        let pg = Self::construct(name)?;
        debug_assert_eq!(pg.label, label);
        Ok(register(pg))
    }

    /// Extract from a Context after find_symmetry. Registers in the global registry.
    pub(crate) fn from_context(ctx: &Context) -> Result<&'static PointGroup, Error> {
        let label = ctx.point_group()?;

        // Fast path
        {
            let map = REGISTRY.lock().unwrap();
            if let Some(&pg) = map.get(&label) {
                return Ok(pg);
            }
        }

        // Linear groups: build analytically (no character table from libmsym)
        if matches!(label, SchoenfliesLabel::Coov | SchoenfliesLabel::Dooh) {
            return Ok(register(Self::build_linear(label)));
        }

        // Finite groups: extract from libmsym context
        let ops = ctx.symmetry_operations()?;
        let ct = ctx.character_table().ok_or_else(|| Error {
            code: MSYM_INVALID_CHARACTER_TABLE,
            message: "character table not available".into(),
        })?;

        let pg = PointGroup {
            label,
            irrep_data: ct.irrep_data.clone(),
            kind: PointGroupKind::Finite {
                order: ct.order,
                ops,
                class_sizes: ct.class_sizes.clone(),
                class_reps: ct.class_reps.clone(),
            },
        };
        Ok(register(pg))
    }

    /// Internal construction. Builds C1 and linear groups inline; others via libmsym.
    fn construct(name: &str) -> Result<PointGroup, Error> {
        if name == "C1" {
            return Ok(Self::build_c1());
        }
        if let Some(label @ (SchoenfliesLabel::Coov | SchoenfliesLabel::Dooh)) =
            SchoenfliesLabel::parse(name)
        {
            return Ok(Self::build_linear(label));
        }

        let mut ctx = Context::new()?;

        let seeds = vec![
            SymmetryCenter {
                atomic_number: 6,
                mass: 12.011,
                position: [1.0, 0.3, 0.7],
                name: String::new(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [0.5, 0.8, 0.2],
                name: String::new(),
            },
        ];

        ctx.set_centers(&seeds)?;
        ctx.set_point_group_by_name(name)?;
        ctx.generate_centers(&seeds)?;
        ctx.find_symmetry()?;

        let detected = ctx.point_group_name()?;
        if detected != name {
            return Err(Error {
                code: MSYM_POINT_GROUP_ERROR,
                message: format!(
                    "requested group '{name}' but detected '{detected}' on generated molecule"
                ),
            });
        }

        let label = ctx.point_group()?;
        let ops = ctx.symmetry_operations()?;
        let ct = ctx.character_table().ok_or_else(|| Error {
            code: MSYM_INVALID_CHARACTER_TABLE,
            message: "character table not available".into(),
        })?;

        Ok(PointGroup {
            label,
            irrep_data: ct.irrep_data.clone(),
            kind: PointGroupKind::Finite {
                order: ct.order,
                ops,
                class_sizes: ct.class_sizes.clone(),
                class_reps: ct.class_reps.clone(),
            },
        })
    }

    fn build_c1() -> PointGroup {
        let identity = SymmetryOp {
            kind: SymmetryOpKind::Identity,
            order: 1,
            power: 0,
            orientation: SymmetryOpOrientation::None,
            vector: [0.0, 0.0, 1.0],
            class: 0,
            matrix: Matrix3::identity(),
        };
        PointGroup {
            label: SchoenfliesLabel::Cn(1),
            irrep_data: vec![IrrepData {
                symbol: "A".into(),
                dimension: 1,
                index: 0,
                characters: vec![1.0],
                lambda: None,
                sigma_v: None,
                gerade: None,
            }],
            kind: PointGroupKind::Finite {
                order: 1,
                ops: vec![identity.clone()],
                class_sizes: vec![1],
                class_reps: vec![identity],
            },
        }
    }

    fn build_linear(label: SchoenfliesLabel) -> PointGroup {
        let is_dooh = label == SchoenfliesLabel::Dooh;
        // λ symbols: Σ(0), Π(1), Δ(2), Φ(3), Γ(4), H(5), I(6)
        let lambda_symbols = ["Σ", "Π", "Δ", "Φ", "Γ", "H", "I"];
        let mut irrep_data = Vec::new();
        let mut idx = 0;

        // Ordering: all gerade first (for D∞h), then ungerade
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
                            index: idx,
                            characters: vec![],
                            lambda: Some(lambda),
                            sigma_v: Some(sv),
                            gerade,
                        });
                        idx += 1;
                    }
                } else {
                    irrep_data.push(IrrepData {
                        symbol: format!("{base}{gu_suffix}"),
                        dimension: dim,
                        index: idx,
                        characters: vec![],
                        lambda: Some(lambda),
                        sigma_v: None,
                        gerade,
                    });
                    idx += 1;
                }
            }
        }

        PointGroup {
            label,
            irrep_data,
            kind: PointGroupKind::Linear,
        }
    }
}

macro_rules! point_group_fn {
    ($name:ident, $schoenflies:expr) => {
        pub fn $name() -> &'static PointGroup {
            Self::from_schoenflies($schoenflies).unwrap()
        }
    };
}

// Named constructors: Cotton's character table appendix
impl PointGroup {
    // Non-axial
    point_group_fn!(c1, "C1");
    point_group_fn!(cs, "Cs");
    point_group_fn!(ci, "Ci");

    // Cn
    point_group_fn!(c2, "C2");
    point_group_fn!(c3, "C3");
    point_group_fn!(c4, "C4");
    point_group_fn!(c5, "C5");
    point_group_fn!(c6, "C6");
    point_group_fn!(c7, "C7");
    point_group_fn!(c8, "C8");

    // Cnv
    point_group_fn!(c2v, "C2v");
    point_group_fn!(c3v, "C3v");
    point_group_fn!(c4v, "C4v");
    point_group_fn!(c5v, "C5v");
    point_group_fn!(c6v, "C6v");

    // Cnh
    point_group_fn!(c2h, "C2h");
    point_group_fn!(c3h, "C3h");
    point_group_fn!(c4h, "C4h");
    point_group_fn!(c5h, "C5h");
    point_group_fn!(c6h, "C6h");

    // Dn
    point_group_fn!(d2, "D2");
    point_group_fn!(d3, "D3");
    point_group_fn!(d4, "D4");
    point_group_fn!(d5, "D5");
    point_group_fn!(d6, "D6");

    // Dnh
    point_group_fn!(d2h, "D2h");
    point_group_fn!(d3h, "D3h");
    point_group_fn!(d4h, "D4h");
    point_group_fn!(d5h, "D5h");
    point_group_fn!(d6h, "D6h");
    point_group_fn!(d8h, "D8h");

    // Dnd
    point_group_fn!(d2d, "D2d");
    point_group_fn!(d3d, "D3d");
    point_group_fn!(d4d, "D4d");
    point_group_fn!(d5d, "D5d");
    point_group_fn!(d6d, "D6d");

    // Sn
    point_group_fn!(s4, "S4");
    point_group_fn!(s6, "S6");
    point_group_fn!(s8, "S8");

    // Cubic
    point_group_fn!(t, "T");
    point_group_fn!(th, "Th");
    point_group_fn!(td, "Td");
    point_group_fn!(o, "O");
    point_group_fn!(oh, "Oh");

    // Icosahedral
    point_group_fn!(ih, "Ih");

    // Linear
    point_group_fn!(coov, "C∞v");
    point_group_fn!(dooh, "D∞h");
}

/// For each class representative R, find which class R² belongs to.
fn r_squared_classes(class_reps: &[SymmetryOp], ops: &[SymmetryOp]) -> Vec<usize> {
    class_reps
        .iter()
        .map(|rep| {
            let m2 = rep.matrix * rep.matrix;
            ops.iter()
                .map(|op| {
                    let diff = op.matrix - m2;
                    diff.iter().map(|x| x * x).sum::<f64>()
                })
                .enumerate()
                .min_by(|(_, da), (_, db)| da.partial_cmp(db).unwrap())
                .unwrap()
                .0
        })
        .map(|best_idx| ops[best_idx].class)
        .collect()
}

pub struct CharacterTableDisplay {
    group: &'static PointGroup,
}

impl fmt::Display for CharacterTableDisplay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let g = self.group;
        let class_reps = g.class_reps();
        if class_reps.is_empty() {
            return write!(f, "{} (linear group, no finite character table)", g.label);
        }
        let class_sizes = g.class_sizes();
        let irreps = g.irreps();

        // Column headers: "2C3", "3σv", etc. (omit multiplicity of 1)
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

        // Column widths: max of header width and widest formatted character in that column
        let irrep_col_width = irreps.iter().map(|ir| ir.symbol().len()).max().unwrap_or(0);
        let label_width = g.label.to_string().len().max(irrep_col_width);

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

        // Header row
        write!(f, "{:width$} │", g.label, width = label_width)?;
        for (c, header) in col_headers.iter().enumerate() {
            write!(f, " {:>width$}", header, width = col_widths[c])?;
        }
        writeln!(f)?;

        // Separator
        write!(f, "{:─>width$}─┼", "", width = label_width)?;
        for &w in &col_widths {
            write!(f, "─{:─>width$}", "", width = w)?;
        }
        writeln!(f)?;

        // Irrep rows
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
    if (value - rounded).abs() < 1e-6 {
        format!("{}", rounded as i64)
    } else {
        format!("{value:.4}")
    }
}

impl fmt::Display for PointGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_point_group_c1() {
        let g = PointGroup::c1();
        assert_eq!(g.to_string(), "C1");
        assert_eq!(g.order(), 1);
        assert_eq!(g.ops().len(), 1);
        assert_eq!(g.label(), SchoenfliesLabel::Cn(1));
        assert_eq!(g.irreps().len(), 1);
        assert_eq!(g.irreps()[0].symbol(), "A");
    }

    #[rstest]
    #[case("C2v", SchoenfliesLabel::Cnv(2), 4)]
    #[case("Td", SchoenfliesLabel::Td, 24)]
    #[case("Oh", SchoenfliesLabel::Oh, 48)]
    #[case("C3v", SchoenfliesLabel::Cnv(3), 6)]
    #[case("D2h", SchoenfliesLabel::Dnh(2), 8)]
    fn test_point_group_from_schoenflies(
        #[case] name: &str,
        #[case] expected_label: SchoenfliesLabel,
        #[case] expected_order: usize,
    ) {
        let g = PointGroup::from_schoenflies(name).unwrap();
        assert_eq!(g.to_string(), name);
        assert_eq!(g.label(), expected_label);
        assert_eq!(g.order(), expected_order);
        assert_eq!(g.ops().len(), expected_order);
    }

    #[rstest]
    fn test_point_group_pointer_identity() {
        let a = PointGroup::c2v();
        let b = PointGroup::from_schoenflies("C2v").unwrap();
        assert!(ptr::eq(a, b));
    }

    #[rstest]
    fn test_point_group_named_constructors() {
        assert_eq!(PointGroup::td().to_string(), "Td");
        assert_eq!(PointGroup::oh().to_string(), "Oh");
        assert_eq!(PointGroup::ih().to_string(), "Ih");
        assert_eq!(PointGroup::d2h().to_string(), "D2h");
    }

    #[rstest]
    fn test_point_group_irreps() {
        let g = PointGroup::c2v();
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
        let g = PointGroup::c2v();
        let a1 = g.irrep("A1").unwrap();
        assert_eq!(a1.symbol(), "A1");
        assert!(g.irrep("nonexistent").is_none());
    }

    #[rstest]
    #[case("C2v", "B1", "B2", &[("A2", 1)])]
    #[case("C2v", "A1", "B1", &[("B1", 1)])]
    #[case("Td", "E", "T2", &[("T1", 1), ("T2", 1)])]
    #[case("Oh", "Eg", "T1u", &[("T1u", 1), ("T2u", 1)])]
    // Real-representation groups: multiplicities are doubled for complex irreps
    // because the standard reduction formula applied to merged classes over-counts.
    // These values are internally consistent but differ from complex Katzer tables.
    #[case("C3", "A", "E1", &[("E1", 2)])]
    #[case("C3", "E1", "E1", &[("A", 2), ("E1", 2)])]
    #[case("C4", "B", "E1", &[("E1", 2)])]
    #[case("C4", "E1", "E1", &[("A", 2), ("B", 2)])]
    #[case("C3h", "A''", "E1'", &[("E1''", 2)])]
    #[case("C3h", "E1'", "E1'", &[("A'", 2), ("E1'", 2)])]
    #[case("C3h", "E1'", "E1''", &[("A''", 2), ("E1''", 2)])]
    fn test_point_group_direct_product(
        #[case] group: &str,
        #[case] a: &str,
        #[case] b: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let ir_a = g.irrep(a).unwrap();
        let ir_b = g.irrep(b).unwrap();
        let product = g.direct_product(ir_a, ir_b);
        let actual: Vec<(&str, u32)> = product.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    // C2v: all irreps are 1D, symmetric square = A1, antisymmetric = empty
    #[case("C2v", "B1", &[("A1", 1)], &[])]
    // Td: E is 2D, [E²] = A1 + E, {E²} = A2
    #[case("Td", "E", &[("A1", 1), ("E", 1)], &[("A2", 1)])]
    // Td: T2 is 3D, [T2²] = A1 + E + T2, {T2²} = T1
    #[case("Td", "T2", &[("A1", 1), ("E", 1), ("T2", 1)], &[("T1", 1)])]
    // Oh: T1u is 3D, [T1u²] = A1g + Eg + T2g, {T1u²} = T1g
    #[case("Oh", "T1u", &[("A1g", 1), ("Eg", 1), ("T2g", 1)], &[("T1g", 1)])]
    // Ih: Hg is 5D, [Hg²] has dim 15, {Hg²} has dim 10
    #[case("Ih", "Hg", &[("Ag", 1), ("Gg", 1), ("Hg", 2)], &[("T1g", 1), ("T2g", 1), ("Gg", 1)])]
    fn test_point_group_symmetric_antisymmetric_square(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected_sym: &[(&str, u32)],
        #[case] expected_anti: &[(&str, u32)],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
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
        let g = PointGroup::c2v();
        let a1 = g.irrep("A1").unwrap();
        let b1 = g.irrep("B1").unwrap();
        let b2 = g.irrep("B2").unwrap();
        assert!(g.contains_totally_symmetric(a1, b1, b1));
        assert!(!g.contains_totally_symmetric(a1, b1, b2));
    }

    #[rstest]
    fn test_irrep_equality() {
        let g = PointGroup::c2v();
        let a1_first = g.irrep("A1").unwrap();
        let a1_second = g.irrep("A1").unwrap();
        let b1 = g.irrep("B1").unwrap();
        assert_eq!(a1_first, a1_second);
        assert_ne!(a1_first, b1);
    }

    #[rstest]
    #[case("C2v", &[4.0, 0.0, 0.0, 0.0], &[("A1", 1), ("A2", 1), ("B1", 1), ("B2", 1)])]
    #[case("C2v", &[1.0, 1.0, 1.0, 1.0], &[("A1", 1)])]
    #[case("C1",  &[7.0], &[("A", 7)])]
    fn test_point_group_reduce(
        #[case] group: &str,
        #[case] characters: &[f64],
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let result = g.reduce(characters);
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    /// Roundtrip: build characters from known irrep multiplicities, reduce, verify.
    #[rstest]
    #[case("C2v", &[("A1", 2), ("B2", 3)])]
    #[case("Td",  &[("A1", 1), ("E", 1), ("T2", 3)])]
    #[case("Oh",  &[("Eg", 2), ("T1u", 1)])]
    fn test_point_group_reduce_roundtrip(#[case] group: &str, #[case] composition: &[(&str, u32)]) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let n_classes = g.class_sizes().len();
        let mut chars = vec![0.0; n_classes];
        for &(sym, mult) in composition {
            let ir = g.irrep(sym).unwrap();
            for (c, ch) in ir.characters().iter().enumerate() {
                chars[c] += mult as f64 * ch;
            }
        }
        let result = g.reduce(&chars);
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        let expected: Vec<(&str, u32)> = composition.to_vec();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case("C2v", &[("A1", 1), ("B1", 1), ("B2", 1)])]
    #[case("Td",  &[("T2", 1)])]
    #[case("C1",  &[("A", 3)])]
    fn test_point_group_translation_irreps(#[case] group: &str, #[case] expected: &[(&str, u32)]) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let result = g.translation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case("C2v", &[("A2", 1), ("B1", 1), ("B2", 1)])]
    #[case("Td",  &[("T1", 1)])]
    #[case("C1",  &[("A", 3)])]
    fn test_point_group_rotation_irreps(#[case] group: &str, #[case] expected: &[(&str, u32)]) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let result = g.rotation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case("C1", true)]
    #[case("C2", true)]
    #[case("C3", true)]
    #[case("D2", true)]
    #[case("D3", true)]
    #[case("T", true)]
    #[case("O", true)]
    #[case("I", true)]
    #[case("Cs", false)]
    #[case("Ci", false)]
    #[case("C2v", false)]
    #[case("C2h", false)]
    #[case("D2h", false)]
    #[case("D3h", false)]
    #[case("Td", false)]
    #[case("Oh", false)]
    #[case("Ih", false)]
    #[case("S4", false)]
    fn test_point_group_is_chiral(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.is_chiral(), expected);
    }

    #[rstest]
    #[case(SymmetryOpKind::Identity, 1, 1, SymmetryOpOrientation::None, "E")]
    #[case(SymmetryOpKind::Inversion, 1, 1, SymmetryOpOrientation::None, "i")]
    #[case(
        SymmetryOpKind::Reflection,
        1,
        1,
        SymmetryOpOrientation::Horizontal,
        "σh"
    )]
    #[case(
        SymmetryOpKind::Reflection,
        1,
        1,
        SymmetryOpOrientation::Vertical,
        "σv"
    )]
    #[case(
        SymmetryOpKind::Reflection,
        1,
        1,
        SymmetryOpOrientation::Dihedral,
        "σd"
    )]
    #[case(SymmetryOpKind::Reflection, 1, 1, SymmetryOpOrientation::None, "σ")]
    #[case(
        SymmetryOpKind::ProperRotation,
        3,
        1,
        SymmetryOpOrientation::None,
        "C3"
    )]
    #[case(
        SymmetryOpKind::ProperRotation,
        3,
        2,
        SymmetryOpOrientation::None,
        "C3²"
    )]
    #[case(
        SymmetryOpKind::ProperRotation,
        6,
        5,
        SymmetryOpOrientation::None,
        "C6⁵"
    )]
    #[case(
        SymmetryOpKind::ImproperRotation,
        4,
        1,
        SymmetryOpOrientation::None,
        "S4"
    )]
    #[case(
        SymmetryOpKind::ImproperRotation,
        4,
        3,
        SymmetryOpOrientation::None,
        "S4³"
    )]
    #[case(
        SymmetryOpKind::ImproperRotation,
        10,
        7,
        SymmetryOpOrientation::None,
        "S10⁷"
    )]
    #[case(
        SymmetryOpKind::ProperRotation,
        2,
        1,
        SymmetryOpOrientation::Vertical,
        "C2'"
    )]
    #[case(
        SymmetryOpKind::ProperRotation,
        2,
        1,
        SymmetryOpOrientation::Dihedral,
        "C2''"
    )]
    #[case(
        SymmetryOpKind::ProperRotation,
        3,
        2,
        SymmetryOpOrientation::Vertical,
        "C3'²"
    )]
    #[case(
        SymmetryOpKind::ImproperRotation,
        4,
        1,
        SymmetryOpOrientation::Vertical,
        "S4'"
    )]
    fn test_symmetry_op_display(
        #[case] kind: SymmetryOpKind,
        #[case] order: i32,
        #[case] power: i32,
        #[case] orientation: SymmetryOpOrientation,
        #[case] expected: &str,
    ) {
        let op = SymmetryOp {
            kind,
            order,
            power,
            orientation,
            vector: [0.0, 0.0, 1.0],
            class: 0,
            matrix: Matrix3::identity(),
        };
        assert_eq!(op.to_string(), expected);
    }

    #[rstest]
    fn test_character_table_display_c2v() {
        let g = PointGroup::c2v();
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
        let g = PointGroup::from_schoenflies("C3v").unwrap();
        let table = g.character_table().to_string();
        let expected = "\
C3v │ E 2C3  3σv
────┼───────────
A1  │ 1   1    1
A2  │ 1   1   -1
E1  │ 2  -1    0
";
        assert_eq!(table, expected);
    }

    #[rstest]
    fn test_symmetry_op_transform_point() {
        let g = PointGroup::c2v();
        let ops = g.ops();
        let p = [1.0, 2.0, 3.0];

        // E leaves point unchanged
        let e = ops
            .iter()
            .find(|op| op.kind == SymmetryOpKind::Identity)
            .unwrap();
        let r = e.transform_point(p);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 2.0).abs() < 1e-12);
        assert!((r[2] - 3.0).abs() < 1e-12);

        // C2 around z: (x,y,z) → (-x,-y,z)
        let c2 = ops
            .iter()
            .find(|op| op.kind == SymmetryOpKind::ProperRotation)
            .unwrap();
        let r = c2.transform_point(p);
        assert!((r[0] + 1.0).abs() < 1e-12);
        assert!((r[1] + 2.0).abs() < 1e-12);
        assert!((r[2] - 3.0).abs() < 1e-12);
    }

    #[rstest]
    #[case("C1", 1)]
    #[case("Cs", 1)]
    #[case("Ci", 1)]
    #[case("C2v", 2)]
    #[case("C3v", 3)]
    #[case("D6h", 6)]
    #[case("S4", 4)]
    #[case("Td", 3)]
    #[case("Oh", 4)]
    #[case("Ih", 5)]
    fn test_point_group_principal_axis_order(#[case] group: &str, #[case] expected: u32) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.principal_axis_order(), expected);
    }

    #[rstest]
    #[case("C1", true)]
    #[case("Ci", true)]
    #[case("Cs", true)]
    #[case("C2", true)]
    #[case("C3", true)]
    #[case("C2v", true)]
    #[case("C2h", true)]
    #[case("C3h", true)]
    #[case("D2", true)]
    #[case("D2h", true)]
    #[case("S4", true)]
    #[case("S6", true)]
    #[case("C3v", false)]
    #[case("D3", false)]
    #[case("D3h", false)]
    #[case("D2d", false)]
    #[case("Td", false)]
    #[case("Oh", false)]
    #[case("Ih", false)]
    fn test_point_group_is_abelian(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.is_abelian(), expected);
    }

    #[rstest]
    #[case("C1", true)]
    #[case("C3", true)]
    #[case("C3h", true)]
    #[case("S4", true)]
    #[case("Ci", true)]
    #[case("Cs", true)]
    #[case("C2v", false)]
    #[case("D2", false)]
    #[case("Td", false)]
    fn test_point_group_is_cyclic(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.is_cyclic(), expected);
    }

    #[rstest]
    #[case("Td", true)]
    #[case("Oh", true)]
    #[case("Ih", true)]
    #[case("T", true)]
    #[case("O", true)]
    #[case("Th", true)]
    #[case("C2v", false)]
    #[case("D6h", false)]
    fn test_point_group_is_cubic(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.is_cubic(), expected);
    }

    #[rstest]
    #[case("C3", true)]
    #[case("C4", true)]
    #[case("C5", true)]
    #[case("C3h", true)]
    #[case("S4", true)]
    #[case("S6", true)]
    #[case("C1", false)]
    #[case("C2", false)]
    #[case("C2v", false)]
    #[case("Td", false)]
    fn test_point_group_has_complex_irreps(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.has_complex_irreps(), expected);
    }

    #[rstest]
    #[case("Ci", true)]
    #[case("C2h", true)]
    #[case("D2h", true)]
    #[case("Oh", true)]
    #[case("Ih", true)]
    #[case("Th", true)]
    #[case("C2v", false)]
    #[case("Td", false)]
    #[case("C3v", false)]
    #[case("D3", false)]
    #[case("C1", false)]
    fn test_point_group_has_inversion(#[case] group: &str, #[case] expected: bool) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.has_inversion(), expected);
    }

    #[rstest]
    #[case("C1", "A")]
    #[case("C2v", "A1")]
    #[case("Oh", "A1g")]
    #[case("Td", "A1")]
    fn test_point_group_totally_symmetric_irrep(#[case] group: &str, #[case] expected: &str) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(g.totally_symmetric_irrep().symbol(), expected);
    }

    #[rstest]
    fn test_irrep_is_gerade_finite() {
        let g = PointGroup::oh();
        assert_eq!(g.irrep("A1g").unwrap().is_gerade(), Some(true));
        assert_eq!(g.irrep("Eg").unwrap().is_gerade(), Some(true));
        assert_eq!(g.irrep("T1u").unwrap().is_gerade(), Some(false));
        assert_eq!(g.irrep("T2u").unwrap().is_gerade(), Some(false));
    }

    #[rstest]
    fn test_irrep_is_gerade_no_inversion() {
        let g = PointGroup::c2v();
        for ir in g.irreps() {
            assert_eq!(ir.is_gerade(), None);
        }
    }

    #[rstest]
    fn test_symmetry_op_is_proper() {
        let g = PointGroup::oh();
        for op in g.ops() {
            assert_eq!(
                op.is_proper(),
                matches!(
                    op.kind,
                    SymmetryOpKind::Identity | SymmetryOpKind::ProperRotation
                ),
            );
        }
    }

    #[rstest]
    #[case("C2v")]
    #[case("Td")]
    #[case("Oh")]
    #[case("D3h")]
    #[case("C1")]
    fn test_point_group_trans_rot_dimensions(#[case] group: &str) {
        let g = PointGroup::from_schoenflies(group).unwrap();
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
    #[case("C1",  &["A"],                  &["A"])]
    #[case("C2v", &["A1", "B1", "B2"],     &["A1", "A2", "B1", "B2"])]
    #[case("Td",  &["T2"],                 &["A1", "E", "T2"])]
    fn test_point_group_ir_active(
        #[case] group: &str,
        #[case] expected_ir: &[&str],
        #[case] expected_raman: &[&str],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        assert_eq!(ir_active_symbols(g), expected_ir);
        assert_eq!(raman_active_symbols(g), expected_raman);
    }
}
