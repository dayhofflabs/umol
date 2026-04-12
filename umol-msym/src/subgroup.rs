use std::fmt;

use nalgebra::Matrix3;

use crate::error::ReductionError;
use crate::irrep::Irrep;
use crate::point_group::PointGroup;
use crate::types::SchoenfliesSymbol;

#[derive(Debug, Clone)]
pub(crate) struct SubgroupData {
    pub(crate) symbol: SchoenfliesSymbol,
    pub(crate) name: String,
    pub(crate) order: usize,
    /// For each subgroup op, the 3×3 matrix in the parent coordinate frame
    /// and the parent-group class index.
    pub(crate) parent_ops: Vec<(Matrix3<f64>, usize)>,
    /// How many inequivalent embeddings of this subgroup type exist in the parent.
    pub(crate) multiplicity: usize,
}

/// A specific subgroup embedding within a parent point group. Holds the parent
/// group reference, an opaque index for FFI, and the subgroup's operation data.
#[derive(Debug, Clone)]
pub struct Subgroup {
    parent: &'static PointGroup,
    index: usize,
    data: SubgroupData,
}

impl Subgroup {
    pub(crate) fn new(
        parent: &'static PointGroup,
        index: usize,
        data: SubgroupData,
    ) -> Self {
        Self { parent, index, data }
    }

    pub fn parent(&self) -> &'static PointGroup {
        self.parent
    }

    pub fn symbol(&self) -> SchoenfliesSymbol {
        self.data.symbol
    }

    pub fn name(&self) -> &str {
        &self.data.name
    }

    pub fn order(&self) -> usize {
        self.data.order
    }

    /// The subgroup's operations as 3×3 matrices in the parent coordinate frame,
    /// paired with their parent-group class indices.
    pub fn parent_ops(&self) -> &[(Matrix3<f64>, usize)] {
        &self.data.parent_ops
    }

    /// How many inequivalent embeddings of this subgroup type exist in the parent.
    /// Returns 1 when this is the only embedding of its type.
    pub fn multiplicity(&self) -> usize {
        self.data.multiplicity
    }

    pub(crate) fn index(&self) -> usize {
        self.index
    }
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

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use nalgebra::{Rotation3, Vector3};
    use rstest::rstest;

    use super::*;
    use crate::types::SchoenfliesSymbol;

    fn cn(axis: Vector3<f64>, n: u32) -> Matrix3<f64> {
        Rotation3::new(axis.normalize() * (2.0 * PI / n as f64)).into_inner()
    }

    // --- Subgroup accessors ---

    #[rstest]
    fn test_subgroup_accessors() {
        let parent = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let c2_mat = cn(Vector3::z(), 2);
        let sg = Subgroup::new(
            parent,
            3,
            SubgroupData {
                symbol: SchoenfliesSymbol::Cn(2),
                name: "C2".into(),
                order: 2,
                parent_ops: vec![(Matrix3::identity(), 0), (c2_mat, 1)],
                multiplicity: 1,
            },
        );
        assert_eq!(sg.parent().symbol(), SchoenfliesSymbol::Cnv(2));
        assert_eq!(sg.symbol(), SchoenfliesSymbol::Cn(2));
        assert_eq!(sg.name(), "C2");
        assert_eq!(sg.order(), 2);
        assert_eq!(sg.parent_ops().len(), 2);
        assert_eq!(sg.multiplicity(), 1);
        assert_eq!(sg.index(), 3);
    }

    // --- correlation_table ---

    #[rstest]
    fn test_correlation_table_c2v_to_cs() {
        let parent = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let child = PointGroup::from_symbol(SchoenfliesSymbol::Cs).unwrap();

        // Find the parent class index containing σ operations.
        // Cs has 2 classes: E and σ. Map E→E (0→0), σ→one of the parent σ classes.
        let parent_irreps = parent.irreps();
        // In C2v: A1=(1,1,1,1), A2=(1,1,-1,-1), B1=(1,-1,1,-1), B2=(1,-1,-1,1)
        // Cs: A'=(1,1), A''=(1,-1)
        // Using σv (parent class 2): A1→A', A2→A'', B1→A', B2→A''
        // Using σv' (parent class 3): A1→A', A2→A'', B1→A'', B2→A'
        // Either is valid; pick class 2.
        let class_map = &[0, 2];

        let ct = correlation_table(parent, child, class_map).unwrap();
        assert_eq!(ct.rows.len(), parent_irreps.len());

        let child_irreps = child.irreps();
        let a_prime = &child_irreps[0];
        let a_double_prime = &child_irreps[1];

        // A1 → A', A2 → A'', B1 → A', B2 → A''
        assert_eq!(ct.rows[0].len(), 1);
        assert_eq!(ct.rows[0][0].0, *a_prime);
        assert_eq!(ct.rows[1].len(), 1);
        assert_eq!(ct.rows[1][0].0, *a_double_prime);
        assert_eq!(ct.rows[2].len(), 1);
        assert_eq!(ct.rows[2][0].0, *a_prime);
        assert_eq!(ct.rows[3].len(), 1);
        assert_eq!(ct.rows[3][0].0, *a_double_prime);
    }

    #[rstest]
    fn test_correlation_table_display() {
        let parent = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let child = PointGroup::from_symbol(SchoenfliesSymbol::Cs).unwrap();
        let ct = correlation_table(parent, child, &[0, 2]).unwrap();
        let s = ct.to_string();
        assert!(s.contains("C2v → Cs"));
        assert!(s.contains("A1 │"));
        assert!(s.contains("A2 │"));
    }
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
