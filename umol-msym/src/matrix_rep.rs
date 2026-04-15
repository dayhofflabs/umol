//! Symmetry group matrix representations.

use std::{fmt, ptr};

use nalgebra::{Matrix3, Vector3};

use crate::point_group::{PointGroup, SymmetryOp};

/// Per-molecule matrix realization of a point group.
///
/// `PointGroup` carries only abstract, orientation-independent data (multiplication
/// table, classes, characters). The 3×3 matrices and axis vectors live here, tied
/// to a specific coordinate frame produced by libmsym for a particular molecule
/// or a canonical name-based placement.
///
/// Slot `i` corresponds to `PointGroup::ops()[i]` (same canonical order).
#[derive(Debug, Clone)]
pub struct MatrixRep {
    group: &'static PointGroup,
    matrices: Vec<Matrix3<f64>>,
    axes: Vec<Vector3<f64>>,
}

impl MatrixRep {
    pub(crate) fn new(
        group: &'static PointGroup,
        matrices: Vec<Matrix3<f64>>,
        axes: Vec<Vector3<f64>>,
    ) -> Self {
        assert_eq!(
            matrices.len(),
            axes.len(),
            "MatrixRep: matrix/axis count mismatch",
        );
        if let Some(finite) = group.finite_data() {
            assert_eq!(
                finite.op_data.len(),
                matrices.len(),
                "MatrixRep op count mismatch: group has {}, got {}",
                finite.op_data.len(),
                matrices.len(),
            );
        }
        Self {
            group,
            matrices,
            axes,
        }
    }

    pub fn identity_only(group: &'static PointGroup) -> Self {
        Self {
            group,
            matrices: vec![Matrix3::identity()],
            axes: vec![Vector3::zeros()],
        }
    }

    pub fn group(&self) -> &'static PointGroup {
        self.group
    }

    pub fn order(&self) -> usize {
        self.matrices.len()
    }

    pub fn matrix(&self, op: SymmetryOp) -> &Matrix3<f64> {
        assert!(
            ptr::eq(op.group(), self.group),
            "SymmetryOp group does not match MatrixRep group",
        );
        &self.matrices[op.index()]
    }

    pub fn axis(&self, op: SymmetryOp) -> &Vector3<f64> {
        assert!(
            ptr::eq(op.group(), self.group),
            "SymmetryOp group does not match MatrixRep group",
        );
        &self.axes[op.index()]
    }

    pub fn transform_point(&self, op: SymmetryOp, p: Vector3<f64>) -> Vector3<f64> {
        self.matrix(op) * p
    }

    pub fn matrices(&self) -> &[Matrix3<f64>] {
        &self.matrices
    }

    pub fn axes(&self) -> &[Vector3<f64>] {
        &self.axes
    }

    pub fn iter(&self) -> impl Iterator<Item = (SymmetryOp, &Matrix3<f64>, &Vector3<f64>)> + '_ {
        self.group
            .ops()
            .into_iter()
            .zip(self.matrices.iter().zip(self.axes.iter()))
            .map(|(op, (m, a))| (op, m, a))
    }
}

impl fmt::Display for MatrixRep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "MatrixRep({}, order={})",
            self.group.symbol(),
            self.order()
        )
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;
    use rstest::rstest;

    use super::*;
    use crate::context::Context;
    use crate::point_group::SymmetryOpKind;
    use crate::types::SymmetryCenter;

    fn water_rep() -> MatrixRep {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&[
            SymmetryCenter {
                atomic_number: 8,
                mass: 15.999,
                position: Vector3::new(0.0, 0.0, 0.117_370_3),
                name: "O".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: Vector3::new(0.0, 0.757_160_4, -0.469_481_2),
                name: "H".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: Vector3::new(0.0, -0.757_160_4, -0.469_481_2),
                name: "H".into(),
            },
        ])
        .unwrap();
        ctx.find_symmetry().unwrap();
        ctx.matrix_rep().unwrap()
    }

    #[rstest]
    fn test_matrix_rep_identity_only() {
        let g = group!(C1);
        let rep = MatrixRep::identity_only(g);
        assert_eq!(rep.order(), 1);
        assert!(ptr::eq(rep.group(), g));
        assert_eq!(rep.matrices().len(), 1);
        assert_eq!(rep.axes().len(), 1);
        assert!((rep.matrices()[0] - Matrix3::identity()).norm() < 1e-12);
        assert!(rep.axes()[0].norm() < 1e-12);
    }

    #[rstest]
    fn test_matrix_rep_order() {
        let rep = water_rep();
        assert_eq!(rep.group().symbol().to_string(), "C2v");
        assert_eq!(rep.order(), 4);
    }

    #[rstest]
    fn test_matrix_rep_identity_op() {
        let rep = water_rep();
        let identity_op = rep
            .group()
            .ops()
            .into_iter()
            .find(|op| op.kind() == SymmetryOpKind::Identity)
            .unwrap();
        let m = rep.matrix(identity_op);
        assert!((m - Matrix3::identity()).norm() < 1e-12);
    }

    #[rstest]
    fn test_matrix_rep_matrices_orthogonal() {
        let rep = water_rep();
        for (op, m, _) in rep.iter() {
            let mmt = m * m.transpose();
            assert!(
                (mmt - Matrix3::identity()).norm() < 1e-10,
                "matrix for op {:?} is not orthogonal",
                op.kind()
            );
            assert!(
                (m.determinant().abs() - 1.0).abs() < 1e-10,
                "det != ±1 for op {:?}",
                op.kind()
            );
        }
    }

    #[rstest]
    fn test_matrix_rep_proper_det() {
        let rep = water_rep();
        for op in rep.group().ops() {
            let det = rep.matrix(op).determinant();
            if op.is_proper() {
                assert!((det - 1.0).abs() < 1e-10, "proper op has det != +1");
            } else {
                assert!((det + 1.0).abs() < 1e-10, "improper op has det != -1");
            }
        }
    }

    #[rstest]
    fn test_matrix_rep_transform_point() {
        let rep = water_rep();
        let p = Vector3::new(1.0, 2.0, 3.0);
        let identity_op = rep
            .group()
            .ops()
            .into_iter()
            .find(|op| op.kind() == SymmetryOpKind::Identity)
            .unwrap();
        let q = rep.transform_point(identity_op, p);
        assert!((q - p).norm() < 1e-12);
    }

    #[rstest]
    fn test_matrix_rep_iter() {
        let rep = water_rep();
        let items: Vec<_> = rep.iter().collect();
        assert_eq!(items.len(), 4);
        for (i, (op, m, a)) in items.iter().enumerate() {
            assert_eq!(op.index(), i);
            assert!(ptr::eq(*m, &rep.matrices()[i]));
            assert!(ptr::eq(*a, &rep.axes()[i]));
        }
    }

    #[rstest]
    fn test_matrix_rep_display() {
        let rep = water_rep();
        assert_eq!(format!("{}", rep), "MatrixRep(C2v, order=4)");

        let g = group!(C1);
        let rep = MatrixRep::identity_only(g);
        assert_eq!(format!("{}", rep), "MatrixRep(C1, order=1)");
    }

    #[rstest]
    #[should_panic(expected = "SymmetryOp group does not match MatrixRep group")]
    fn test_matrix_rep_matrix_error() {
        let rep = water_rep();
        let other = group!(D2h);
        let op = other.ops()[0];
        rep.matrix(op);
    }

    #[rstest]
    #[should_panic(expected = "SymmetryOp group does not match MatrixRep group")]
    fn test_matrix_rep_axis_error() {
        let rep = water_rep();
        let other = group!(D2h);
        let op = other.ops()[0];
        rep.axis(op);
    }
}
