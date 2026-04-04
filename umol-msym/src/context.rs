use std::ffi::CString;
use std::os::raw::c_int;

use umol_msym_sys as ffi;

use crate::error::{self, Error};
use crate::types::*;

pub struct Context {
    ctx: ffi::msym_context,
    character_table: Option<CharacterTable>,
}

unsafe impl Send for Context {}

impl Drop for Context {
    fn drop(&mut self) {
        unsafe {
            ffi::msymReleaseContext(self.ctx);
        }
    }
}

impl Context {
    // -------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------

    pub fn new() -> Result<Self, Error> {
        let ctx = unsafe { ffi::msymCreateContext() };
        if ctx.is_null() {
            return Err(Error {
                code: ffi::MSYM_MEMORY_ERROR,
                message: "failed to create context".into(),
            });
        }
        Ok(Self {
            ctx,
            character_table: None,
        })
    }

    // -------------------------------------------------------------------
    // Thresholds
    // -------------------------------------------------------------------

    pub fn set_thresholds(&mut self, t: &Thresholds) -> Result<(), Error> {
        let ft = t.to_ffi();
        error::check(unsafe { ffi::msymSetThresholds(self.ctx, &ft) })
    }

    pub fn thresholds(&self) -> Result<Thresholds, Error> {
        let mut ptr = std::ptr::null();
        error::check(unsafe { ffi::msymGetThresholds(self.ctx, &mut ptr) })?;
        Ok(Thresholds::from_ffi(unsafe { &*ptr }))
    }

    // -------------------------------------------------------------------
    // Elements
    // -------------------------------------------------------------------

    pub fn set_elements(&mut self, elements: &[SymmetryElement]) -> Result<(), Error> {
        self.character_table = None;
        let mut ffi_elements: Vec<ffi::msym_element_t> =
            elements.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymSetElements(self.ctx, ffi_elements.len() as c_int, ffi_elements.as_mut_ptr())
        })
    }

    pub fn elements(&self) -> Result<Vec<SymmetryElement>, Error> {
        let mut len: c_int = 0;
        let mut ptr = std::ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { std::slice::from_raw_parts(ptr, len as usize) };
        Ok(slice.iter().map(SymmetryElement::from_ffi).collect())
    }

    // -------------------------------------------------------------------
    // Symmetry detection
    // -------------------------------------------------------------------

    pub fn find_symmetry(&mut self) -> Result<(), Error> {
        self.character_table = None;
        error::check(unsafe { ffi::msymFindSymmetry(self.ctx) })?;
        // Cache the character table
        let mut ct_ptr = std::ptr::null();
        error::check(unsafe { ffi::msymGetCharacterTable(self.ctx, &mut ct_ptr) })?;
        self.character_table = Some(unsafe { CharacterTable::from_ffi(&*ct_ptr) });
        Ok(())
    }

    // -------------------------------------------------------------------
    // Point group queries
    // -------------------------------------------------------------------

    pub fn point_group(&self) -> Result<(PointGroupType, i32), Error> {
        let mut pg_type: ffi::msym_point_group_type_t = 0;
        let mut n: c_int = 0;
        error::check(unsafe { ffi::msymGetPointGroupType(self.ctx, &mut pg_type, &mut n) })?;
        Ok((pg_type.into(), n))
    }

    pub fn point_group_name(&self) -> Result<String, Error> {
        let mut buf = [0i8; 16];
        error::check(unsafe { ffi::msymGetPointGroupName(self.ctx, 16, buf.as_mut_ptr()) })?;
        let name = unsafe { std::ffi::CStr::from_ptr(buf.as_ptr()) };
        Ok(name.to_string_lossy().into_owned())
    }

    pub fn set_point_group_by_name(&mut self, name: &str) -> Result<(), Error> {
        self.character_table = None;
        let cname = CString::new(name).map_err(|_| Error {
            code: ffi::MSYM_INVALID_INPUT,
            message: "invalid group name".into(),
        })?;
        error::check(unsafe { ffi::msymSetPointGroupByName(self.ctx, cname.as_ptr()) })
    }

    pub fn set_point_group(&mut self, type_: PointGroupType, n: i32) -> Result<(), Error> {
        self.character_table = None;
        error::check(unsafe { ffi::msymSetPointGroupByType(self.ctx, type_.to_ffi(), n) })
    }

    // -------------------------------------------------------------------
    // Symmetry operations
    // -------------------------------------------------------------------

    pub fn symmetry_operations(&self) -> Result<Vec<SymmetryOperation>, Error> {
        let mut len: c_int = 0;
        let mut ptr = std::ptr::null();
        error::check(unsafe { ffi::msymGetSymmetryOperations(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { std::slice::from_raw_parts(ptr, len as usize) };
        Ok(slice.iter().map(SymmetryOperation::from).collect())
    }

    // -------------------------------------------------------------------
    // Character table
    // -------------------------------------------------------------------

    pub fn character_table(&self) -> Option<&CharacterTable> {
        self.character_table.as_ref()
    }

    // -------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------

    pub fn geometry(&self) -> Result<Geometry, Error> {
        let mut g: ffi::msym_geometry_t = 0;
        error::check(unsafe { ffi::msymGetGeometry(self.ctx, &mut g) })?;
        Ok(g.into())
    }

    pub fn center_of_mass(&self) -> Result<[f64; 3], Error> {
        let mut v = [0.0f64; 3];
        error::check(unsafe { ffi::msymGetCenterOfMass(self.ctx, v.as_mut_ptr()) })?;
        Ok(v)
    }

    pub fn radius(&self) -> Result<f64, Error> {
        let mut r = 0.0f64;
        error::check(unsafe { ffi::msymGetRadius(self.ctx, &mut r) })?;
        Ok(r)
    }

    pub fn principal_moments(&self) -> Result<[f64; 3], Error> {
        let mut eigval = [0.0f64; 3];
        error::check(unsafe { ffi::msymGetPrincipalMoments(self.ctx, eigval.as_mut_ptr()) })?;
        Ok(eigval)
    }

    pub fn principal_axes(&self) -> Result<[[f64; 3]; 3], Error> {
        let mut eigvec = [[0.0f64; 3]; 3];
        error::check(unsafe { ffi::msymGetPrincipalAxes(self.ctx, eigvec.as_mut_ptr()) })?;
        Ok(eigvec)
    }

    // -------------------------------------------------------------------
    // Equivalence sets
    // -------------------------------------------------------------------

    pub fn equivalence_sets(&self) -> Result<Vec<EquivalenceSet>, Error> {
        let mut len: c_int = 0;
        let mut ptr = std::ptr::null();
        error::check(unsafe { ffi::msymGetEquivalenceSets(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { std::slice::from_raw_parts(ptr, len as usize) };
        Ok(slice
            .iter()
            .map(|es| {
                let elems = unsafe { std::slice::from_raw_parts(es.elements, es.length as usize) };
                EquivalenceSet {
                    elements: elems
                        .iter()
                        .map(|e| SymmetryElement::from_ffi(unsafe { &**e }))
                        .collect(),
                    max_error: es.err,
                }
            })
            .collect())
    }

    // -------------------------------------------------------------------
    // Symmetrization
    // -------------------------------------------------------------------

    pub fn symmetrize_elements(&mut self) -> Result<f64, Error> {
        self.character_table = None;
        let mut err = 0.0f64;
        error::check(unsafe { ffi::msymSymmetrizeElements(self.ctx, &mut err) })?;
        Ok(err)
    }

    // -------------------------------------------------------------------
    // Alignment
    // -------------------------------------------------------------------

    pub fn align_axes(&mut self) -> Result<(), Error> {
        error::check(unsafe { ffi::msymAlignAxes(self.ctx) })
    }

    pub fn alignment_transform(&self) -> Result<[[f64; 3]; 3], Error> {
        let mut transform = [[0.0f64; 3]; 3];
        error::check(unsafe {
            ffi::msymGetAlignmentTransform(self.ctx, transform.as_mut_ptr())
        })?;
        Ok(transform)
    }

    // -------------------------------------------------------------------
    // Element generation (from asymmetric unit)
    // -------------------------------------------------------------------

    pub fn generate_elements(&mut self, asymmetric_unit: &[SymmetryElement]) -> Result<(), Error> {
        self.character_table = None;
        let mut ffi_elements: Vec<ffi::msym_element_t> =
            asymmetric_unit.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymGenerateElements(
                self.ctx,
                ffi_elements.len() as c_int,
                ffi_elements.as_mut_ptr(),
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn water() -> Vec<SymmetryElement> {
        vec![
            SymmetryElement {
                atomic_number: 8,
                mass: 15.999,
                position: [0.0, 0.0, 0.117_370_3],
                name: "O".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [0.0, 0.757_160_4, -0.469_481_2],
                name: "H".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [0.0, -0.757_160_4, -0.469_481_2],
                name: "H".into(),
            },
        ]
    }

    fn methane() -> Vec<SymmetryElement> {
        vec![
            SymmetryElement {
                atomic_number: 6,
                mass: 12.011,
                position: [0.0, 0.0, 0.0],
                name: "C".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [0.629_118_5, 0.629_118_5, 0.629_118_5],
                name: "H".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [-0.629_118_5, -0.629_118_5, 0.629_118_5],
                name: "H".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [-0.629_118_5, 0.629_118_5, -0.629_118_5],
                name: "H".into(),
            },
            SymmetryElement {
                atomic_number: 1,
                mass: 1.008,
                position: [0.629_118_5, -0.629_118_5, -0.629_118_5],
                name: "H".into(),
            },
        ]
    }

    #[rstest]
    #[case(water(), PointGroupType::Cnv, 2, "C2v")]
    #[case(methane(), PointGroupType::Td, 3, "Td")]
    fn test_context_find_symmetry(
        #[case] elements: Vec<SymmetryElement>,
        #[case] expected_type: PointGroupType,
        #[case] expected_n: i32,
        #[case] expected_name: &str,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_elements(&elements).unwrap();
        ctx.find_symmetry().unwrap();

        let (pg_type, n) = ctx.point_group().unwrap();
        assert_eq!(pg_type, expected_type);
        assert_eq!(n, expected_n);
        assert_eq!(ctx.point_group_name().unwrap(), expected_name);
    }

    #[rstest]
    #[case(water(), 4)]   // C2v: E, C2, σv, σv'
    #[case(methane(), 24)] // Td: 24 operations
    fn test_context_symmetry_operations(
        #[case] elements: Vec<SymmetryElement>,
        #[case] expected_count: usize,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_elements(&elements).unwrap();
        ctx.find_symmetry().unwrap();
        assert_eq!(ctx.symmetry_operations().unwrap().len(), expected_count);
    }

    #[rstest]
    fn test_context_character_table() {
        let mut ctx = Context::new().unwrap();
        ctx.set_elements(&water()).unwrap();
        ctx.find_symmetry().unwrap();

        let ct = ctx.character_table().unwrap();

        // C2v has 4 irreps: A1, A2, B1, B2
        assert_eq!(ct.irreps.len(), 4);
        assert_eq!(ct.order, 4);

        let names: Vec<&str> = ct.irreps.iter().map(|ir| ir.name.as_str()).collect();
        assert!(names.contains(&"A1"));
        assert!(names.contains(&"A2"));
        assert!(names.contains(&"B1"));
        assert!(names.contains(&"B2"));

        // All irreps in C2v are 1-dimensional
        for irrep in &ct.irreps {
            assert_eq!(irrep.dimension, 1);
        }
    }

    #[rstest]
    fn test_character_table_direct_product() {
        let mut ctx = Context::new().unwrap();
        ctx.set_elements(&water()).unwrap();
        ctx.find_symmetry().unwrap();
        let ct = ctx.character_table().unwrap();

        let b1 = ct.irreps.iter().find(|ir| ir.name == "B1").unwrap();
        let b2 = ct.irreps.iter().find(|ir| ir.name == "B2").unwrap();

        // B1 ⊗ B2 = A2 in C2v
        let product = ct.direct_product(b1, b2);
        assert_eq!(product.len(), 1);
        assert_eq!(product[0].0.name, "A2");
        assert_eq!(product[0].1, 1);
    }

    #[rstest]
    fn test_character_table_contains_totally_symmetric() {
        let mut ctx = Context::new().unwrap();
        ctx.set_elements(&water()).unwrap();
        ctx.find_symmetry().unwrap();
        let ct = ctx.character_table().unwrap();

        let a1 = ct.irreps.iter().find(|ir| ir.name == "A1").unwrap();
        let b1 = ct.irreps.iter().find(|ir| ir.name == "B1").unwrap();
        let b2 = ct.irreps.iter().find(|ir| ir.name == "B2").unwrap();

        // A1 ⊗ B1 ⊗ B1 contains A1 (dipole-allowed if operator transforms as B1)
        assert!(ct.contains_totally_symmetric(a1, b1, b1));
        // A1 ⊗ B1 ⊗ B2 does not contain A1
        assert!(!ct.contains_totally_symmetric(a1, b1, b2));
    }
}
