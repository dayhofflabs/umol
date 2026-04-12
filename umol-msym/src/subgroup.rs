use std::fmt;

use nalgebra::Matrix3;

use crate::error::ReductionError;
use crate::irrep::Irrep;
use crate::point_group::PointGroup;
use crate::types::SchoenfliesSymbol;

/// Subgroup descent transport: parent-group ops that survive into a target
/// subgroup, together with the libmsym subgroup-table index used to drive
/// `msymSelectSubgroup`.
#[derive(Debug, Clone)]
pub struct SubgroupData {
    pub symbol: SchoenfliesSymbol,
    pub name: String,
    pub order: usize,
    /// For each subgroup op, the 3×3 matrix in the parent coordinate frame
    /// and the parent-group class index.
    pub parent_ops: Vec<(Matrix3<f64>, usize)>,
    pub(crate) index: usize,
}

/// Correlation table mapping parent irreps to child irrep decompositions.
#[derive(Debug)]
pub struct CorrelationTable {
    pub parent: &'static PointGroup,
    pub child: &'static PointGroup,
    pub rows: Vec<Vec<(Irrep, u32)>>,
}

/// Build a correlation table from a class_map.
///
/// `class_map[i]` = the parent class index that child class `i` maps to.
pub fn correlation_table(
    parent: &'static PointGroup,
    child: &'static PointGroup,
    class_map: &[usize],
) -> Result<CorrelationTable, ReductionError> {
    let parent_irreps = parent.irreps();

    let mut rows = Vec::with_capacity(parent_irreps.len());
    for parent_irrep in &parent_irreps {
        let restricted: Vec<f64> = class_map
            .iter()
            .map(|&parent_class| parent_irrep.characters()[parent_class])
            .collect();
        let decomp = child.reduce(&restricted)?;
        rows.push(decomp);
    }

    Ok(CorrelationTable {
        parent,
        child,
        rows,
    })
}

impl fmt::Display for CorrelationTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let parent_irreps = self.parent.irreps();

        writeln!(f, "{} → {}", self.parent.symbol(), self.child.symbol())?;

        let label_width = parent_irreps
            .iter()
            .map(|ir| ir.symbol().len())
            .max()
            .unwrap_or(4);

        for (i, row) in self.rows.iter().enumerate() {
            write!(
                f,
                "{:>width$} │ ",
                parent_irreps[i].symbol(),
                width = label_width
            )?;
            let decomp: Vec<String> = row
                .iter()
                .map(|(ir, n)| {
                    if *n == 1 {
                        ir.symbol().to_string()
                    } else {
                        format!("{n}{}", ir.symbol())
                    }
                })
                .collect();
            writeln!(f, "{}", decomp.join(" + "))?;
        }

        Ok(())
    }
}
