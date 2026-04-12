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
            self.group.label(),
            self.order()
        )
    }
}
