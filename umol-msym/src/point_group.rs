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

        let ct = &self.character_table;
        let d = ct.irrep_data.len();
        let h = ct.order as f64;

        let product_chars: Vec<f64> = (0..d)
            .map(|c| a.data.characters[c] * b.data.characters[c])
            .collect();

        let mut result = Vec::new();
        for ir_data in &ct.irrep_data {
            let n: f64 = (0..d)
                .map(|c| ct.class_sizes[c] as f64 * ir_data.characters[c] * product_chars[c])
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

        ctx.set_elements(&seeds)?;
        ctx.set_point_group_by_name(name)?;
        ctx.generate_elements(&seeds)?;
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
}
