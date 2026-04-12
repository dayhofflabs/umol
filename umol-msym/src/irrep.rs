//! Irrep data structures for umol-msym.

use std::hash::{Hash, Hasher};
use std::{fmt, ptr};

use crate::point_group::{PointGroup, SymmetryOpKind, SymmetryOpOrientation};
use crate::thresholds::COMPLEX_IRREP_NORM;

#[derive(Debug, Clone)]
pub(crate) enum ReductionData {
    Finite { characters: Vec<f64> },
    Linear {
        lambda: u32,
        sigma_v: Option<bool>,
        gerade: Option<bool>,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct IrrepData {
    pub(crate) symbol: String,
    pub(crate) dimension: i32,
    pub(crate) reduction: ReductionData,
}

impl IrrepData {
    pub(crate) fn characters(&self) -> &[f64] {
        match &self.reduction {
            ReductionData::Finite { characters } => characters,
            ReductionData::Linear { .. } => &[],
        }
    }
}

#[derive(Copy, Clone)]
pub struct Irrep {
    pub(crate) data: &'static IrrepData,
    pub(crate) group: &'static PointGroup,
}

impl Irrep {
    pub fn symbol(&self) -> &str {
        &self.data.symbol
    }

    pub fn dimension(&self) -> i32 {
        self.data.dimension
    }

    pub fn characters(&self) -> &[f64] {
        self.data.characters()
    }

    pub fn lambda(&self) -> Option<u32> {
        match &self.data.reduction {
            ReductionData::Linear { lambda, .. } => Some(*lambda),
            ReductionData::Finite { .. } => None,
        }
    }

    pub fn gerade(&self) -> Option<bool> {
        match &self.data.reduction {
            ReductionData::Linear { gerade, .. } => *gerade,
            ReductionData::Finite { characters } => {
                let finite = self.group.finite.as_ref()?;
                let inv_index = finite
                    .op_data
                    .iter()
                    .position(|d| d.kind == SymmetryOpKind::Inversion)?;
                let inv_class = finite.op_data[inv_index].class;
                let chi = characters.get(inv_class)?;
                Some(*chi > 0.0)
            }
        }
    }

    pub fn vertical_mirror_symmetric(&self) -> Option<bool> {
        match &self.data.reduction {
            ReductionData::Linear { lambda, sigma_v, .. } if *lambda == 0 => *sigma_v,
            _ => None,
        }
    }

    pub fn principal_axis_symmetric(&self) -> Option<bool> {
        if self.data.dimension != 1 {
            return None;
        }
        let finite = self.group.finite.as_ref()?;
        let n = self.group.principal_axis_order();
        let cn_index = finite.op_data.iter().position(|d| {
            d.kind == SymmetryOpKind::ProperRotation && d.order == n as i32 && d.power == 1
        })?;
        let cn_class = finite.op_data[cn_index].class;
        let chi = self.data.characters().get(cn_class)?;
        Some(*chi > 0.0)
    }

    pub fn horizontal_mirror_symmetric(&self) -> Option<bool> {
        if self.data.dimension != 1 {
            return None;
        }
        let finite = self.group.finite.as_ref()?;
        let sh_index = finite.op_data.iter().position(|d| {
            d.kind == SymmetryOpKind::Reflection
                && d.orientation == SymmetryOpOrientation::Horizontal
        })?;
        let sh_class = finite.op_data[sh_index].class;
        let chi = self.data.characters().get(sh_class)?;
        Some(*chi > 0.0)
    }

    pub fn totally_symmetric(&self) -> bool {
        match &self.data.reduction {
            ReductionData::Finite { characters } => characters.iter().all(|&c| (c - 1.0).abs() < 1e-10),
            ReductionData::Linear {
                lambda,
                sigma_v,
                gerade,
            } => *lambda == 0 && *sigma_v == Some(true) && *gerade != Some(false),
        }
    }

    pub fn complex(&self) -> bool {
        if self.data.dimension != 2 {
            return false;
        }
        let Some(finite) = self.group.finite.as_ref() else {
            return false;
        };
        let h = finite.order as f64;
        let norm_sq: f64 = self
            .data
            .characters()
            .iter()
            .zip(&finite.class_sizes)
            .map(|(chi, &size)| size as f64 * chi * chi)
            .sum();
        (norm_sq - 2.0 * h).abs() < COMPLEX_IRREP_NORM
    }

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
        write!(f, "Irrep({}:{})", self.group.symbol, self.data.symbol)
    }
}

impl fmt::Display for Irrep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.data.symbol)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use crate::point_group::PointGroup;

    #[rstest]
    #[case::c2v_a1("C2v", "A1", true)]
    #[case::c2v_b1("C2v", "B1", false)]
    #[case::oh_a1g("Oh", "A1g", true)]
    #[case::oh_t1u("Oh", "T1u", false)]
    fn test_irrep_totally_symmetric(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected: bool,
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();
        assert_eq!(ir.totally_symmetric(), expected);
    }

    #[rstest]
    fn test_irrep_totally_symmetric_linear() {
        let g = PointGroup::coov();
        assert!(g.irrep("Σ+").unwrap().totally_symmetric());
        assert!(!g.irrep("Σ-").unwrap().totally_symmetric());
        assert!(!g.irrep("Π").unwrap().totally_symmetric());

        let g = PointGroup::dooh();
        assert!(g.irrep("Σ+g").unwrap().totally_symmetric());
        assert!(!g.irrep("Σ+u").unwrap().totally_symmetric());
        assert!(!g.irrep("Σ-g").unwrap().totally_symmetric());
    }

    #[rstest]
    #[case::c2v_a1("C2v", "A1", Some(true))]
    #[case::c2v_a2("C2v", "A2", Some(true))]
    #[case::c2v_b1("C2v", "B1", Some(false))]
    #[case::c3v_a1("C3v", "A1", Some(true))]
    #[case::c3v_e("C3v", "E", None)]
    #[case::d4h_a1g("D4h", "A1g", Some(true))]
    #[case::d4h_b1g("D4h", "B1g", Some(false))]
    fn test_irrep_principal_axis_symmetric(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected: Option<bool>,
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();
        assert_eq!(ir.principal_axis_symmetric(), expected);
    }

    #[rstest]
    #[case::c2h_ag("C2h", "Ag", Some(true))]
    #[case::c2h_au("C2h", "Au", Some(false))]
    #[case::c2h_bg("C2h", "Bg", Some(false))]
    #[case::c2h_bu("C2h", "Bu", Some(true))]
    #[case::d6h_a1g("D6h", "A1g", Some(true))]
    #[case::d6h_a2u("D6h", "A2u", Some(false))]
    #[case::c2v_no_sigma_h("C2v", "A1", None)]
    fn test_irrep_horizontal_mirror_symmetric(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected: Option<bool>,
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();
        assert_eq!(ir.horizontal_mirror_symmetric(), expected);
    }

    #[rstest]
    #[case::oh_a1g("Oh", "A1g", Some(true))]
    #[case::oh_a2g("Oh", "A2g", Some(true))]
    #[case::oh_t1u("Oh", "T1u", Some(false))]
    #[case::d2h_a1g("D2h", "A1g", Some(true))]
    #[case::d2h_b1u("D2h", "B1u", Some(false))]
    #[case::c2v_no_inversion("C2v", "A1", None)]
    fn test_irrep_gerade(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected: Option<bool>,
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();
        assert_eq!(ir.gerade(), expected);
    }

    #[rstest]
    fn test_irrep_vertical_mirror_symmetric_linear() {
        let g = PointGroup::coov();
        assert_eq!(g.irrep("Σ+").unwrap().vertical_mirror_symmetric(), Some(true));
        assert_eq!(g.irrep("Σ-").unwrap().vertical_mirror_symmetric(), Some(false));
        assert_eq!(g.irrep("Π").unwrap().vertical_mirror_symmetric(), None);
    }

    #[rstest]
    #[case::c3_e("C3", "E", true)]
    #[case::c4_e("C4", "E", true)]
    #[case::c5_e1("C5", "E1", true)]
    #[case::c5_e2("C5", "E2", true)]
    #[case::c3v_e_real("C3v", "E", false)]
    #[case::c2v_a1_1d("C2v", "A1", false)]
    fn test_irrep_complex(
        #[case] group: &str,
        #[case] irrep: &str,
        #[case] expected: bool,
    ) {
        let g = PointGroup::parse(group).unwrap();
        let ir = g.irrep(irrep).unwrap();
        assert_eq!(ir.complex(), expected);
    }
}
