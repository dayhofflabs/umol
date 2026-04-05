use std::fmt;

use nalgebra::Matrix3;
use umol_msym_sys::{MSYM_INVALID_CHARACTER_TABLE, MSYM_POINT_GROUP_ERROR};

use crate::context::Context;
use crate::error::Error;
use crate::types::{
    CharacterTable, Irrep, PointGroupKind, SymmetryCenter, SymmetryOp, SymmetryOpKind,
    SymmetryOpOrientation,
};

/// A molecular point group: abstract algebraic object with operations and character table.
#[derive(Debug, Clone)]
pub struct PointGroup {
    pub kind: PointGroupKind,
    /// Fold number of the principal rotation (or improper rotation) axis.
    pub n: u32,
    pub name: String,
    /// Number of group elements (= sum of class sizes).
    pub order: usize,
    pub operations: Vec<SymmetryOp>,
    pub character_table: CharacterTable,
}

impl PointGroup {
    /// The trivial point group (identity only).
    pub fn c1() -> Self {
        let identity = SymmetryOp {
            kind: SymmetryOpKind::Identity,
            order: 1,
            power: 0,
            orientation: SymmetryOpOrientation::None,
            vector: [0.0, 0.0, 1.0],
            class: 0,
            matrix: Matrix3::identity(),
        };
        let irrep = Irrep {
            name: "A".into(),
            dimension: 1,
            index: 0,
        };
        Self {
            kind: PointGroupKind::Cn,
            n: 1,
            name: "C1".into(),
            order: 1,
            operations: vec![identity.clone()],
            character_table: CharacterTable {
                irreps: vec![irrep],
                class_sizes: vec![1],
                class_operations: vec![identity],
                characters: vec![vec![1.0]],
                order: 1,
            },
        }
    }

    /// Construct a point group by Schoenflies symbol (e.g. "C2v", "Td", "Oh").
    ///
    /// Uses libmsym internally: sets the group on a seed molecule and generates
    /// the full orbit so that find_symmetry recovers exactly the requested group.
    pub fn from_schoenflies(name: &str) -> Result<Self, Error> {
        if name == "C1" {
            return Ok(Self::c1());
        }

        let mut ctx = Context::new()?;

        // Seed: two atoms of different elements at generic positions (not on any
        // symmetry element). generate_elements creates their orbits under the
        // requested group, producing a molecule with exactly that symmetry.
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

        Self::from_context(&ctx)
    }

    /// Extract a PointGroup from a Context after find_symmetry has been called.
    pub fn from_context(ctx: &Context) -> Result<Self, Error> {
        let (kind, n) = ctx.point_group()?;
        let name = ctx.point_group_name()?;
        let operations = ctx.symmetry_operations()?;
        let character_table = ctx
            .character_table()
            .ok_or_else(|| Error {
                code: MSYM_INVALID_CHARACTER_TABLE,
                message: "character table not available".into(),
            })?
            .clone();
        let order = character_table.order;

        Ok(Self {
            kind,
            n: n as u32,
            name,
            order,
            operations,
            character_table,
        })
    }
}

impl fmt::Display for PointGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.name)
    }
}

impl PartialEq for PointGroup {
    fn eq(&self, other: &Self) -> bool {
        self.kind == other.kind && self.n == other.n
    }
}

impl Eq for PointGroup {}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_point_group_c1() {
        let g = PointGroup::c1();
        assert_eq!(g.name, "C1");
        assert_eq!(g.order, 1);
        assert_eq!(g.operations.len(), 1);
        assert_eq!(g.character_table.irreps.len(), 1);
        assert_eq!(g.character_table.irreps[0].name, "A");
    }

    #[rstest]
    #[case("C2v", 4)]
    #[case("Td", 24)]
    #[case("Oh", 48)]
    #[case("C3v", 6)]
    #[case("D2h", 8)]
    fn test_point_group_from_schoenflies(#[case] name: &str, #[case] expected_order: usize) {
        let g = PointGroup::from_schoenflies(name).unwrap();
        assert_eq!(g.name, name);
        assert_eq!(g.order, expected_order);
        assert_eq!(g.operations.len(), expected_order);
    }
}
