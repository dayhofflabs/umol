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
    // Centers (atoms with positions and masses)
    // -------------------------------------------------------------------

    pub fn set_centers(&mut self, centers: &[SymmetryCenter]) -> Result<(), Error> {
        self.character_table = None;
        let mut ffi_elems: Vec<ffi::msym_element_t> =
            centers.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymSetElements(self.ctx, ffi_elems.len() as c_int, ffi_elems.as_mut_ptr())
        })
    }

    pub fn centers(&self) -> Result<Vec<SymmetryCenter>, Error> {
        let mut len: c_int = 0;
        let mut ptr = std::ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { std::slice::from_raw_parts(ptr, len as usize) };
        Ok(slice.iter().map(SymmetryCenter::from_ffi).collect())
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

    pub fn point_group(&self) -> Result<SchoenfliesLabel, Error> {
        let mut pg_type: ffi::msym_point_group_type_t = 0;
        let mut n: c_int = 0;
        error::check(unsafe { ffi::msymGetPointGroupType(self.ctx, &mut pg_type, &mut n) })?;
        Ok(SchoenfliesLabel::from_ffi(pg_type, n))
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

    pub fn set_point_group(&mut self, label: SchoenfliesLabel) -> Result<(), Error> {
        self.character_table = None;
        let (pg_type, n) = label.to_ffi();
        error::check(unsafe { ffi::msymSetPointGroupByType(self.ctx, pg_type, n) })
    }

    // -------------------------------------------------------------------
    // Symmetry operations
    // -------------------------------------------------------------------

    pub fn symmetry_operations(&self) -> Result<Vec<SymmetryOp>, Error> {
        let mut len: c_int = 0;
        let mut ptr = std::ptr::null();
        error::check(unsafe { ffi::msymGetSymmetryOperations(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { std::slice::from_raw_parts(ptr, len as usize) };
        Ok(slice.iter().map(SymmetryOp::from).collect())
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

    pub fn set_center_of_mass(&mut self, v: [f64; 3]) -> Result<(), Error> {
        let mut v = v;
        error::check(unsafe { ffi::msymSetCenterOfMass(self.ctx, v.as_mut_ptr()) })
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
                    centers: elems
                        .iter()
                        .map(|e| SymmetryCenter::from_ffi(unsafe { &**e }))
                        .collect(),
                    max_error: es.err,
                }
            })
            .collect())
    }

    // -------------------------------------------------------------------
    // Symmetrization
    // -------------------------------------------------------------------

    pub fn symmetrize_centers(&mut self) -> Result<f64, Error> {
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
    // Center generation (from asymmetric unit)
    // -------------------------------------------------------------------

    pub fn generate_centers(&mut self, asymmetric_unit: &[SymmetryCenter]) -> Result<(), Error> {
        self.character_table = None;
        let mut ffi_elems: Vec<ffi::msym_element_t> =
            asymmetric_unit.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymGenerateElements(
                self.ctx,
                ffi_elems.len() as c_int,
                ffi_elems.as_mut_ptr(),
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn water() -> Vec<SymmetryCenter> {
        vec![
            SymmetryCenter {
                atomic_number: 8,
                mass: 15.999,
                position: [0.0, 0.0, 0.117_370_3],
                name: "O".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [0.0, 0.757_160_4, -0.469_481_2],
                name: "H".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [0.0, -0.757_160_4, -0.469_481_2],
                name: "H".into(),
            },
        ]
    }

    fn methane() -> Vec<SymmetryCenter> {
        vec![
            SymmetryCenter {
                atomic_number: 6,
                mass: 12.011,
                position: [0.0, 0.0, 0.0],
                name: "C".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [0.629_118_5, 0.629_118_5, 0.629_118_5],
                name: "H".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [-0.629_118_5, -0.629_118_5, 0.629_118_5],
                name: "H".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [-0.629_118_5, 0.629_118_5, -0.629_118_5],
                name: "H".into(),
            },
            SymmetryCenter {
                atomic_number: 1,
                mass: 1.008,
                position: [0.629_118_5, -0.629_118_5, -0.629_118_5],
                name: "H".into(),
            },
        ]
    }

    #[rstest]
    #[case(water(), SchoenfliesLabel::Cnv(2), "C2v")]
    #[case(methane(), SchoenfliesLabel::Td, "Td")]
    fn test_context_find_symmetry(
        #[case] elements: Vec<SymmetryCenter>,
        #[case] expected_label: SchoenfliesLabel,
        #[case] expected_name: &str,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&elements).unwrap();
        ctx.find_symmetry().unwrap();

        assert_eq!(ctx.point_group().unwrap(), expected_label);
        assert_eq!(ctx.point_group_name().unwrap(), expected_name);
    }

    #[rstest]
    #[case(water(), 4)]   // C2v: E, C2, σv, σv'
    #[case(methane(), 24)] // Td: 24 operations
    fn test_context_symmetry_operations(
        #[case] elements: Vec<SymmetryCenter>,
        #[case] expected_count: usize,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&elements).unwrap();
        ctx.find_symmetry().unwrap();
        assert_eq!(ctx.symmetry_operations().unwrap().len(), expected_count);
    }

    #[rstest]
    fn test_context_character_table() {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&water()).unwrap();
        ctx.find_symmetry().unwrap();

        let ct = ctx.character_table().unwrap();

        // C2v has 4 irreps: A1, A2, B1, B2
        assert_eq!(ct.irrep_data.len(), 4);
        assert_eq!(ct.order, 4);

        let symbols: Vec<&str> = ct.irrep_data.iter().map(|ir| ir.symbol.as_str()).collect();
        assert!(symbols.contains(&"A1"));
        assert!(symbols.contains(&"A2"));
        assert!(symbols.contains(&"B1"));
        assert!(symbols.contains(&"B2"));

        for ir in &ct.irrep_data {
            assert_eq!(ir.dimension, 1);
        }
    }
}
