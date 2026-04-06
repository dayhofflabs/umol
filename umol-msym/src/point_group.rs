use std::collections::HashMap;
use std::fmt;
use std::sync::{LazyLock, Mutex};

use nalgebra::Matrix3;
use umol_msym_sys::{MSYM_INVALID_CHARACTER_TABLE, MSYM_INVALID_INPUT, MSYM_POINT_GROUP_ERROR};

use crate::context::Context;
use crate::error::Error;
use crate::types::{
    CharacterTable, IrrepData, SchoenfliesLabel, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation,
};

// ---------------------------------------------------------------------------
// Registry
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Irrep (view type)
// ---------------------------------------------------------------------------

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

    /// Dimensionality (1 for A/B, 2 for E, 3 for T).
    pub fn dimension(&self) -> i32 {
        self.data.dimension
    }

    /// Character values, one per conjugacy class.
    pub fn characters(&self) -> &[f64] {
        &self.data.characters
    }

    /// Parent point group.
    pub fn group(&self) -> &'static PointGroup {
        self.group
    }
}

impl PartialEq for Irrep {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.data, other.data)
    }
}

impl Eq for Irrep {}

impl std::hash::Hash for Irrep {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
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

// ---------------------------------------------------------------------------
// PointGroup
// ---------------------------------------------------------------------------

/// A molecular point group: singleton algebraic object with operations and character table.
///
/// Point groups are `&'static` singletons. There is exactly one C2v, one Td, etc.
/// Access via named constructors (`PointGroup::c2v()`) or `PointGroup::from_schoenflies("C2v")`.
#[derive(Debug)]
pub struct PointGroup {
    pub(crate) label: SchoenfliesLabel,
    pub(crate) order: usize,
    pub(crate) operations: Vec<SymmetryOp>,
    pub(crate) character_table: CharacterTable,
}

impl PointGroup {
    pub fn label(&self) -> SchoenfliesLabel {
        self.label
    }

    pub fn order(&self) -> usize {
        self.order
    }

    pub fn operations(&self) -> &[SymmetryOp] {
        &self.operations
    }

    pub fn class_sizes(&self) -> &[i32] {
        &self.character_table.class_sizes
    }

    // -----------------------------------------------------------------------
    // Irrep access
    // -----------------------------------------------------------------------

    /// All irreducible representations of this group.
    pub fn irreps(&'static self) -> Vec<Irrep> {
        self.character_table
            .irrep_data
            .iter()
            .map(|d| Irrep {
                data: d,
                group: self,
            })
            .collect()
    }

    /// Look up an irrep by Mulliken symbol.
    pub fn irrep(&'static self, symbol: &str) -> Option<Irrep> {
        self.character_table
            .irrep_data
            .iter()
            .find(|d| d.symbol == symbol)
            .map(|d| Irrep {
                data: d,
                group: self,
            })
    }

    // -----------------------------------------------------------------------
    // Algebraic methods
    // -----------------------------------------------------------------------

    /// Decompose the direct product a ⊗ b into irreps with multiplicities.
    pub fn direct_product(&'static self, a: Irrep, b: Irrep) -> Vec<(Irrep, u32)> {
        debug_assert!(std::ptr::eq(a.group, self));
        debug_assert!(std::ptr::eq(b.group, self));

        let product_chars: Vec<f64> = a
            .data
            .characters
            .iter()
            .zip(&b.data.characters)
            .map(|(ca, cb)| ca * cb)
            .collect();
        self.reduce(&product_chars)
    }

    /// Reduce a representation (given by its class characters) into irreps with multiplicities.
    ///
    /// `characters` must have one entry per conjugacy class, in the same order as the
    /// character table. Returns irreps with non-zero multiplicity.
    pub fn reduce(&'static self, characters: &[f64]) -> Vec<(Irrep, u32)> {
        let ct = &self.character_table;
        let d = ct.irrep_data.len();
        assert_eq!(
            characters.len(),
            d,
            "expected {} class characters, got {}",
            d,
            characters.len()
        );
        let h = ct.order as f64;

        let mut result = Vec::new();
        for ir_data in &ct.irrep_data {
            let n: f64 = (0..d)
                .map(|c| ct.class_sizes[c] as f64 * ir_data.characters[c] * characters[c])
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
    ///
    /// Characters computed from traces of the 3×3 operation matrices.
    pub fn translation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        let ct = &self.character_table;
        let chars: Vec<f64> = ct
            .class_operations
            .iter()
            .map(|op| op.matrix.trace())
            .collect();
        self.reduce(&chars)
    }

    /// Irreps spanned by the rotational degrees of freedom (Rx, Ry, Rz).
    ///
    /// Pseudovector representation: χ_rot(R) = det(M_R) · tr(M_R).
    pub fn rotation_irreps(&'static self) -> Vec<(Irrep, u32)> {
        let ct = &self.character_table;
        let chars: Vec<f64> = ct
            .class_operations
            .iter()
            .map(|op| op.matrix.determinant() * op.matrix.trace())
            .collect();
        self.reduce(&chars)
    }

    /// Irreps of the symmetric square of the vector representation (x², y², z², xy, xz, yz).
    ///
    /// Used for Raman and electric quadrupole selection rules.
    fn quadratic_irreps(&'static self) -> Vec<(Irrep, u32)> {
        let ct = &self.character_table;
        let chars: Vec<f64> = ct
            .class_operations
            .iter()
            .map(|op| {
                let tr = op.matrix.trace();
                let tr2 = (op.matrix * op.matrix).trace();
                (tr * tr + tr2) / 2.0
            })
            .collect();
        self.reduce(&chars)
    }

    /// Electric dipole transition allowed? Checks Γ_i ⊗ Γ(x,y,z) ⊗ Γ_f ⊃ A1.
    pub fn electric_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(std::ptr::eq(initial.group, self));
        debug_assert!(std::ptr::eq(final_.group, self));
        self.translation_irreps()
            .iter()
            .any(|(gamma_t, _)| self.contains_totally_symmetric(initial, *gamma_t, final_))
    }

    /// Magnetic dipole transition allowed? Checks Γ_i ⊗ Γ(Rx,Ry,Rz) ⊗ Γ_f ⊃ A1.
    pub fn magnetic_dipole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(std::ptr::eq(initial.group, self));
        debug_assert!(std::ptr::eq(final_.group, self));
        self.rotation_irreps()
            .iter()
            .any(|(gamma_r, _)| self.contains_totally_symmetric(initial, *gamma_r, final_))
    }

    /// Raman transition allowed? Checks Γ_i ⊗ Γ(x²,y²,...,yz) ⊗ Γ_f ⊃ A1.
    pub fn raman_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        debug_assert!(std::ptr::eq(initial.group, self));
        debug_assert!(std::ptr::eq(final_.group, self));
        self.quadratic_irreps()
            .iter()
            .any(|(gamma_q, _)| self.contains_totally_symmetric(initial, *gamma_q, final_))
    }

    /// Electric quadrupole transition allowed? Same basis as Raman (symmetric square).
    pub fn electric_quadrupole_allowed(&'static self, initial: Irrep, final_: Irrep) -> bool {
        self.raman_allowed(initial, final_)
    }

    /// Whether a ⊗ b ⊗ c contains the totally symmetric representation.
    pub fn contains_totally_symmetric(
        &'static self,
        a: Irrep,
        b: Irrep,
        c: Irrep,
    ) -> bool {
        debug_assert!(std::ptr::eq(a.group, self));
        debug_assert!(std::ptr::eq(b.group, self));
        debug_assert!(std::ptr::eq(c.group, self));

        let ct = &self.character_table;
        let d = ct.irrep_data.len();
        let h = ct.order as f64;

        let n: f64 = (0..d)
            .map(|cls| {
                ct.class_sizes[cls] as f64
                    * a.data.characters[cls]
                    * b.data.characters[cls]
                    * c.data.characters[cls]
            })
            .sum::<f64>()
            / h;
        n.round() as u32 > 0
    }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

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

        // Build and register
        let operations = ctx.symmetry_operations()?;
        let character_table = ctx
            .character_table()
            .ok_or_else(|| Error {
                code: MSYM_INVALID_CHARACTER_TABLE,
                message: "character table not available".into(),
            })?
            .clone();
        let order = character_table.order;

        let pg = PointGroup {
            label,
            order,
            operations,
            character_table,
        };
        Ok(register(pg))
    }

    /// Internal construction. Builds C1 inline; others via libmsym.
    fn construct(name: &str) -> Result<PointGroup, Error> {
        if name == "C1" {
            return Ok(Self::build_c1());
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
        let operations = ctx.symmetry_operations()?;
        let character_table = ctx
            .character_table()
            .ok_or_else(|| Error {
                code: MSYM_INVALID_CHARACTER_TABLE,
                message: "character table not available".into(),
            })?
            .clone();
        let order = character_table.order;

        Ok(PointGroup {
            label,
            order,
            operations,
            character_table,
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
            order: 1,
            operations: vec![identity.clone()],
            character_table: CharacterTable {
                irrep_data: vec![IrrepData {
                    symbol: "A".into(),
                    dimension: 1,
                    index: 0,
                    characters: vec![1.0],
                }],
                class_sizes: vec![1],
                class_operations: vec![identity],
                order: 1,
            },
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
        assert_eq!(g.operations().len(), 1);
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
        assert_eq!(g.operations().len(), expected_order);
    }

    #[rstest]
    fn test_point_group_pointer_identity() {
        let a = PointGroup::c2v();
        let b = PointGroup::from_schoenflies("C2v").unwrap();
        assert!(std::ptr::eq(a, b));
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
            assert!(std::ptr::eq(ir.group(), g));
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
    fn test_point_group_direct_product() {
        let g = PointGroup::c2v();
        let b1 = g.irrep("B1").unwrap();
        let b2 = g.irrep("B2").unwrap();
        let product = g.direct_product(b1, b2);
        assert_eq!(product.len(), 1);
        assert_eq!(product[0].0.symbol(), "A2");
        assert_eq!(product[0].1, 1);
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

    // --- reduce ---

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
    fn test_point_group_reduce_roundtrip(
        #[case] group: &str,
        #[case] composition: &[(&str, u32)],
    ) {
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

    // --- translation_irreps / rotation_irreps ---

    #[rstest]
    #[case("C2v", &[("A1", 1), ("B1", 1), ("B2", 1)])]
    #[case("Td",  &[("T2", 1)])]
    #[case("C1",  &[("A", 3)])]
    fn test_point_group_translation_irreps(
        #[case] group: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let result = g.translation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case("C2v", &[("A2", 1), ("B1", 1), ("B2", 1)])]
    #[case("Td",  &[("T1", 1)])]
    #[case("C1",  &[("A", 3)])]
    fn test_point_group_rotation_irreps(
        #[case] group: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::from_schoenflies(group).unwrap();
        let result = g.rotation_irreps();
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
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

    // --- selection rules ---

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
