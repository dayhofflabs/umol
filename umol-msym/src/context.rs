//! Context wrapper for libmsym.

#![allow(dead_code)]

use std::ffi::{CStr, CString};
use std::os::raw::c_int;
use std::{ptr, slice};

use nalgebra::Vector3;
use umol_msym_sys as ffi;

use crate::basis::{BasisFunction, BasisKind, IrrepBasis, Salc};
use crate::error::{self, MsymError};
use crate::matrix_rep::MatrixRep;
use crate::point_group::{compute_op_matrix, PointGroup};
use crate::subgroup::{Subgroup, SubgroupData};
use crate::thresholds::Thresholds;
use crate::types::*;

pub(crate) struct Context {
    ctx: ffi::msym_context,
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
    pub(crate) fn new() -> Result<Self, MsymError> {
        let ctx = unsafe { ffi::msymCreateContext() };
        if ctx.is_null() {
            return Err(MsymError {
                code: ffi::MSYM_MEMORY_ERROR,
                message: "failed to create context".into(),
            });
        }
        Ok(Self { ctx })
    }

    /// Set numerical thresholds for symmetry detection and comparison. Optional — libmsym
    /// uses built-in defaults if not called.
    pub(crate) fn set_thresholds(&mut self, t: &Thresholds) -> Result<(), MsymError> {
        let ft = t.to_ffi();
        error::check(unsafe { ffi::msymSetThresholds(self.ctx, &ft) })
    }

    /// Read the active thresholds from the context.
    pub(crate) fn thresholds(&self) -> Result<Thresholds, MsymError> {
        let mut ptr = ptr::null();
        error::check(unsafe { ffi::msymGetThresholds(self.ctx, &mut ptr) })?;
        Ok(Thresholds::from_ffi(unsafe { &*ptr }))
    }

    /// Load a complete set of centers into the context. Used in two paths:
    /// `set_centers` → `find_symmetry` → `symmetrize_centers`, or
    /// `set_point_group` + `set_centers` → `symmetrize_centers`.
    pub(crate) fn set_centers(&mut self, centers: &[SymmetryCenter]) -> Result<(), MsymError> {
        let mut ffi_elems: Vec<ffi::msym_element_t> = centers.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymSetElements(self.ctx, ffi_elems.len() as c_int, ffi_elems.as_mut_ptr())
        })
    }

    /// Read back the current centers. Works after any of the three symmetrization paths.
    pub(crate) fn centers(&self) -> Result<Vec<SymmetryCenter>, MsymError> {
        let mut len: c_int = 0;
        let mut ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { slice::from_raw_parts(ptr, len as usize) };
        Ok(slice.iter().map(SymmetryCenter::from_ffi).collect())
    }

    /// Expand an asymmetric unit into the full set of centers under the current point
    /// group. Internally computes equivalence sets, permutations, and symmetrization.
    /// Requires `set_point_group`.
    pub(crate) fn generate_centers(
        &mut self,
        asymmetric_unit: &[SymmetryCenter],
    ) -> Result<(), MsymError> {
        let mut ffi_elems: Vec<ffi::msym_element_t> =
            asymmetric_unit.iter().map(|e| e.to_ffi()).collect();
        error::check(unsafe {
            ffi::msymGenerateElements(self.ctx, ffi_elems.len() as c_int, ffi_elems.as_mut_ptr())
        })
    }

    /// Snap centers to exact symmetry positions. Returns the RMS displacement.
    /// Requires a point group (from `find_symmetry` or `set_point_group`); lazily
    /// computes equivalence sets and permutations if not already cached.
    pub(crate) fn symmetrize_centers(&mut self) -> Result<f64, MsymError> {
        let mut err = 0.0f64;
        error::check(unsafe { ffi::msymSymmetrizeElements(self.ctx, &mut err) })?;
        Ok(err)
    }

    /// Displace one center by a vector while preserving the point group. The translation
    /// is projected onto the totally symmetric subspace: each symmetry operation is applied
    /// to `v`, the results are accumulated on equivalent centers, and the antisymmetric
    /// component is discarded. Requires a point group; lazily computes equivalence sets
    /// and permutations if not already cached.
    pub(crate) fn apply_translation(
        &mut self,
        center_index: usize,
        v: [f64; 3],
    ) -> Result<(), MsymError> {
        let mut len: c_int = 0;
        let mut ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut len, &mut ptr) })?;
        assert!(
            (center_index as c_int) < len,
            "center_index {} out of range ({})",
            center_index,
            len
        );
        let element = unsafe { ptr.add(center_index) };
        let mut v = v;
        error::check(unsafe { ffi::msymApplyTranslation(self.ctx, element, v.as_mut_ptr()) })
    }

    /// Register basis functions (AOs or displacement vectors) on the centers.
    /// Requires `set_centers` and `find_symmetry`. Needed before `generate_salcs` or `salcs`.
    pub(crate) fn set_basis_functions(&mut self, basis: &[BasisFunction]) -> Result<(), MsymError> {
        let mut elem_len: c_int = 0;
        let mut elem_ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut elem_len, &mut elem_ptr) })?;

        let mut ffi_basis: Vec<ffi::msym_basis_function_t> = basis
            .iter()
            .map(|bf| {
                assert!(
                    (bf.atom_index as c_int) < elem_len,
                    "atom_index {} out of range ({})",
                    bf.atom_index,
                    elem_len
                );
                let element = unsafe { elem_ptr.add(bf.atom_index) };
                let type_ = match bf.kind {
                    BasisKind::CartesianHarmonic => ffi::MSYM_BASIS_TYPE_CARTESIAN,
                    _ => ffi::MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC,
                };
                ffi::msym_basis_function_t {
                    id: ptr::null_mut(),
                    type_,
                    element,
                    f: ffi::msym_basis_function_union_t {
                        rsh: ffi::msym_real_spherical_harmonic_t {
                            n: bf.ffi_n(),
                            l: bf.l,
                            m: bf.m,
                        },
                    },
                    name: bf.libmsym_name_bytes(),
                }
            })
            .collect();

        error::check(unsafe {
            ffi::msymSetBasisFunctions(self.ctx, ffi_basis.len() as c_int, ffi_basis.as_mut_ptr())
        })
    }

    /// Read back basis functions from the context. Requires `set_basis_functions`.
    pub(crate) fn basis_functions(&self) -> Result<Vec<BasisFunction>, MsymError> {
        let mut elem_len: c_int = 0;
        let mut elem_ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut elem_len, &mut elem_ptr) })?;

        let mut len: c_int = 0;
        let mut ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetBasisFunctions(self.ctx, &mut len, &mut ptr) })?;
        let bfs = unsafe { slice::from_raw_parts(ptr, len as usize) };
        Ok(bfs
            .iter()
            .map(|bf| BasisFunction::from_ffi(bf, elem_ptr as *const _))
            .collect())
    }

    /// Return the detected or manually set point group as a Schoenflies symbol.
    pub(crate) fn point_group(&self) -> Result<SchoenfliesSymbol, MsymError> {
        let mut pg_type: ffi::msym_point_group_type_t = 0;
        let mut n: c_int = 0;
        error::check(unsafe { ffi::msymGetPointGroupType(self.ctx, &mut pg_type, &mut n) })?;
        Ok(SchoenfliesSymbol::from_ffi(pg_type, n))
    }

    /// Manually set the point group (skipping detection). Alternative to `find_symmetry`.
    pub(crate) fn set_point_group(&mut self, label: SchoenfliesSymbol) -> Result<(), MsymError> {
        let (pg_type, n) = label.to_ffi();
        error::check(unsafe { ffi::msymSetPointGroupByType(self.ctx, pg_type, n) })
    }

    /// Manually set the point group by Schoenflies name string (e.g. "C2v").
    pub(crate) fn set_point_group_by_name(&mut self, name: &str) -> Result<(), MsymError> {
        let cname = CString::new(name).map_err(|_| MsymError {
            code: ffi::MSYM_INVALID_INPUT,
            message: "invalid group name".into(),
        })?;
        error::check(unsafe { ffi::msymSetPointGroupByName(self.ctx, cname.as_ptr()) })
    }

    /// Return the point group name as a string. Requires `find_symmetry` or `set_point_group`.
    pub(crate) fn point_group_name(&self) -> Result<String, MsymError> {
        let mut buf = [0i8; 16];
        error::check(unsafe { ffi::msymGetPointGroupName(self.ctx, 16, buf.as_mut_ptr()) })?;
        let name = unsafe { CStr::from_ptr(buf.as_ptr()) };
        Ok(name.to_string_lossy().into_owned())
    }

    /// List subgroups of the current point group. Requires `find_symmetry`.
    pub(crate) fn subgroups(
        &self,
        parent: &'static PointGroup,
    ) -> Result<Vec<Subgroup>, MsymError> {
        let mut len: c_int = 0;
        let mut ptr = ptr::null();
        error::check(unsafe { ffi::msymGetSubgroups(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { slice::from_raw_parts(ptr, len as usize) };

        let symbols: Vec<SchoenfliesSymbol> = slice
            .iter()
            .map(|sg| SchoenfliesSymbol::from_ffi(sg.type_, sg.n))
            .collect();

        let mut result = Vec::with_capacity(len as usize);
        for (i, sg) in slice.iter().enumerate() {
            let symbol = symbols[i];
            let name = unsafe { CStr::from_ptr(sg.name.as_ptr()) }
                .to_string_lossy()
                .into_owned();
            let order = sg.order as usize;
            let multiplicity = symbols.iter().filter(|&&s| s == symbol).count();

            let sops_slice = unsafe { slice::from_raw_parts(sg.sops, order) };
            let parent_ops = sops_slice
                .iter()
                .map(|sop_ptr| {
                    let sop = unsafe { &**sop_ptr };
                    (compute_op_matrix(sop), sop.cla as usize)
                })
                .collect();

            result.push(Subgroup::new(
                parent,
                i,
                SubgroupData {
                    symbol,
                    name,
                    order,
                    parent_ops,
                    multiplicity,
                },
            ));
        }
        Ok(result)
    }

    /// Lower the context to a subgroup. Requires `subgroups` to obtain `sg`.
    pub(crate) fn select_subgroup(&mut self, sg: &Subgroup) -> Result<(), MsymError> {
        let mut len: c_int = 0;
        let mut ptr = ptr::null();
        error::check(unsafe { ffi::msymGetSubgroups(self.ctx, &mut len, &mut ptr) })?;

        let sg_ptr = unsafe { ptr.add(sg.index()) };
        error::check(unsafe { ffi::msymSelectSubgroup(self.ctx, sg_ptr) })
    }

    /// Build 3×3 matrix representations for each symmetry operation. Requires `find_symmetry`.
    pub(crate) fn symmetry_representation(
        &self,
        group: &'static PointGroup,
    ) -> Result<MatrixRep, MsymError> {
        let mut len: c_int = 0;
        let mut p = ptr::null();
        error::check(unsafe { ffi::msymGetSymmetryOperations(self.ctx, &mut len, &mut p) })?;
        let sops = unsafe { slice::from_raw_parts(p, len as usize) };
        let matrices = sops.iter().map(compute_op_matrix).collect();
        let axes = sops.iter().map(|s| Vector3::from(s.v)).collect();
        Ok(MatrixRep::new(group, matrices, axes))
    }

    /// Detect the point group from the loaded centers. Requires `set_centers`.
    pub(crate) fn find_symmetry(&mut self) -> Result<(), MsymError> {
        error::check(unsafe { ffi::msymFindSymmetry(self.ctx) })
    }

    /// Return all equivalence sets (orbits under the symmetry group).
    /// Requires `find_symmetry`.
    pub(crate) fn equivalence_sets(&self) -> Result<Vec<EquivalenceSet>, MsymError> {
        let mut len: c_int = 0;
        let mut ptr = ptr::null();
        error::check(unsafe { ffi::msymGetEquivalenceSets(self.ctx, &mut len, &mut ptr) })?;
        let slice = unsafe { slice::from_raw_parts(ptr, len as usize) };
        Ok(slice
            .iter()
            .map(|es| {
                let elems = unsafe { slice::from_raw_parts(es.elements, es.length as usize) };
                EquivalenceSet {
                    centers: elems
                        .iter()
                        .map(|e| SymmetryCenter::from_ffi(unsafe { &**e }))
                        .collect(),
                }
            })
            .collect())
    }

    /// Return the equivalence set containing a specific center.
    /// Requires `find_symmetry`.
    pub(crate) fn equivalence_set_for_center(
        &self,
        center_index: usize,
    ) -> Result<EquivalenceSet, MsymError> {
        let mut len: c_int = 0;
        let mut elem_ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetElements(self.ctx, &mut len, &mut elem_ptr) })?;
        assert!(
            (center_index as c_int) < len,
            "center_index {} out of range ({})",
            center_index,
            len
        );
        let element = unsafe { elem_ptr.add(center_index) };
        let mut es_ptr: *const ffi::msym_equivalence_set_t = ptr::null();
        error::check(unsafe {
            ffi::msymGetEquivalenceSetByElement(self.ctx, element, &mut es_ptr)
        })?;
        let es = unsafe { &*es_ptr };
        let elems = unsafe { slice::from_raw_parts(es.elements, es.length as usize) };
        Ok(EquivalenceSet {
            centers: elems
                .iter()
                .map(|e| SymmetryCenter::from_ffi(unsafe { &**e }))
                .collect(),
        })
    }

    /// Compute equivalence sets without running full symmetry detection.
    /// Requires `set_centers` and `set_point_group`.
    pub(crate) fn find_equivalence_sets(&mut self) -> Result<(), MsymError> {
        error::check(unsafe { ffi::msymFindEquivalenceSets(self.ctx) })
    }

    /// Compute how symmetry operations permute elements within equivalence sets.
    /// Requires `find_equivalence_sets` or `find_symmetry`.
    pub(crate) fn find_equivalence_set_permutations(&mut self) -> Result<(), MsymError> {
        error::check(unsafe { ffi::msymFindEquivalenceSetPermutations(self.ctx) })
    }

    /// Return SALCs grouped by irrep as structured data. Each SALC stores sparse
    /// coefficients referencing the basis function indices from `set_basis_functions`.
    /// Requires `generate_salcs`.
    pub(crate) fn salcs_by_irrep(
        &self,
        group: &'static PointGroup,
    ) -> Result<Vec<IrrepBasis>, MsymError> {
        let mut len: c_int = 0;
        let mut srs_ptr: *const ffi::msym_subrepresentation_space_t = ptr::null();
        error::check(unsafe {
            ffi::msymGetSubrepresentationSpaces(self.ctx, &mut len, &mut srs_ptr)
        })?;
        let spaces = unsafe { slice::from_raw_parts(srs_ptr, len as usize) };

        let mut bf_len: c_int = 0;
        let mut bf_ptr = ptr::null_mut();
        error::check(unsafe { ffi::msymGetBasisFunctions(self.ctx, &mut bf_len, &mut bf_ptr) })?;

        let irreps = group.irreps();
        let mut result = Vec::with_capacity(len as usize);
        for srs in spaces {
            let irrep = irreps[srs.s as usize];
            let salcs = unsafe { slice::from_raw_parts(srs.salc, srs.salcl as usize) };
            let mut salc_list = Vec::with_capacity(salcs.len());
            for salc in salcs {
                let fl = salc.fl as usize;
                let pf = salc.pf as *const f64;
                let coeffs = unsafe { slice::from_raw_parts(pf, fl) };
                let bf_ptrs = unsafe { slice::from_raw_parts(salc.f, fl) };
                let mut sparse = Vec::new();
                for j in 0..fl {
                    let c = coeffs[j];
                    if c.abs() > 1e-14 {
                        let bf_index = unsafe {
                            (bf_ptrs[j] as *const ffi::msym_basis_function_t)
                                .offset_from(bf_ptr as *const _)
                        } as usize;
                        sparse.push((bf_index, c));
                    }
                }
                salc_list.push(Salc {
                    coefficients: sparse,
                });
            }
            result.push(IrrepBasis {
                irrep,
                salcs: salc_list,
            });
        }
        Ok(result)
    }

    /// Compute SALC subspaces. Caches results for `salcs` and `salcs_by_irrep`.
    /// Requires `set_centers` → `find_symmetry` → `set_basis_functions`.
    pub(crate) fn generate_salcs(&mut self) -> Result<(), MsymError> {
        error::check(unsafe { ffi::msymGenerateSubrepresentationSpaces(self.ctx) })
    }

    /// Project MO coefficients onto SALC subspaces in-place. The coefficient matrix
    /// is rewritten so each row belongs to a single irrep; `species` is updated with
    /// the irrep index for each row. Requires `generate_salcs`.
    pub(crate) fn symmetrize_salcs(
        &mut self,
        basis_count: usize,
        coefficients: &mut [f64],
        species: &mut [i32],
    ) -> Result<(), MsymError> {
        let l = basis_count as c_int;
        let mut partner: Vec<ffi::msym_partner_function_t> = (0..l)
            .map(|_| ffi::msym_partner_function_t { i: 0, d: 0 })
            .collect();
        error::check(unsafe {
            ffi::msymSymmetrizeWavefunctions(
                self.ctx,
                l,
                coefficients.as_mut_ptr(),
                species.as_mut_ptr(),
                partner.as_mut_ptr(),
            )
        })
    }

    /// Return SALCs as a dense l×l row-major coefficient matrix and a species index
    /// vector mapping each row to an irrep in the character table.
    /// Requires `generate_salcs`.
    pub(crate) fn salcs(&self, basis_count: usize) -> Result<(Vec<f64>, Vec<i32>), MsymError> {
        let l = basis_count as c_int;
        let mut coefficients = vec![0.0f64; (l * l) as usize];
        let mut species = vec![0i32; l as usize];
        let mut partner: Vec<ffi::msym_partner_function_t> = (0..l)
            .map(|_| ffi::msym_partner_function_t { i: 0, d: 0 })
            .collect();

        error::check(unsafe {
            ffi::msymGetSALCs(
                self.ctx,
                l,
                coefficients.as_mut_ptr(),
                species.as_mut_ptr(),
                partner.as_mut_ptr(),
            )
        })?;

        Ok((coefficients, species))
    }

    /// Measure per-irrep projection norms of a single wavefunction vector.
    /// `components` must have length equal to the number of irreps. Read-only diagnostic.
    /// Requires `generate_salcs`.
    pub(crate) fn irrep_components(
        &self,
        wavefunction: &[f64],
        components: &mut [f64],
    ) -> Result<(), MsymError> {
        let mut wf = wavefunction.to_vec();
        error::check(unsafe {
            ffi::msymSymmetrySpeciesComponents(
                self.ctx,
                wf.len() as c_int,
                wf.as_mut_ptr(),
                components.len() as c_int,
                components.as_mut_ptr(),
            )
        })
    }

    /// Return the mass-weighted centroid. Requires `set_centers`.
    pub(crate) fn center_of_mass(&self) -> Result<[f64; 3], MsymError> {
        let mut v = [0.0f64; 3];
        error::check(unsafe { ffi::msymGetCenterOfMass(self.ctx, v.as_mut_ptr()) })?;
        Ok(v)
    }

    /// Override the center of mass used for alignment. Requires `set_centers`.
    pub(crate) fn set_center_of_mass(&mut self, v: [f64; 3]) -> Result<(), MsymError> {
        let mut v = v;
        error::check(unsafe { ffi::msymSetCenterOfMass(self.ctx, v.as_mut_ptr()) })
    }

    /// Return the maximum distance from the center of mass to any center.
    /// Requires `set_centers`.
    pub(crate) fn radius(&self) -> Result<f64, MsymError> {
        let mut r = 0.0f64;
        error::check(unsafe { ffi::msymGetRadius(self.ctx, &mut r) })?;
        Ok(r)
    }

    /// Classify the molecular geometry (linear, planar, polyhedral, etc.).
    /// Requires `set_centers`.
    pub(crate) fn molecular_shape(&self) -> Result<MolecularShape, MsymError> {
        let mut g: ffi::msym_geometry_t = 0;
        error::check(unsafe { ffi::msymGetGeometry(self.ctx, &mut g) })?;
        Ok(g.into())
    }

    /// Return the three principal moments of inertia. Requires `set_centers`.
    pub(crate) fn principal_moments(&self) -> Result<[f64; 3], MsymError> {
        let mut eigval = [0.0f64; 3];
        error::check(unsafe { ffi::msymGetPrincipalMoments(self.ctx, eigval.as_mut_ptr()) })?;
        Ok(eigval)
    }

    /// Return the three principal axes as row vectors. Requires `set_centers`.
    pub(crate) fn principal_axes(&self) -> Result<[[f64; 3]; 3], MsymError> {
        let mut eigvec = [[0.0f64; 3]; 3];
        error::check(unsafe { ffi::msymGetPrincipalAxes(self.ctx, eigvec.as_mut_ptr()) })?;
        Ok(eigvec)
    }

    /// Rotate centers so the principal axis is along z and the secondary along y.
    /// Requires `find_symmetry`.
    pub(crate) fn align_axes(&mut self) -> Result<(), MsymError> {
        error::check(unsafe { ffi::msymAlignAxes(self.ctx) })
    }

    /// Return the primary and secondary alignment axes. Requires `find_symmetry`.
    pub(crate) fn alignment_axes(&self) -> Result<([f64; 3], [f64; 3]), MsymError> {
        let mut primary = [0.0f64; 3];
        let mut secondary = [0.0f64; 3];
        error::check(unsafe {
            ffi::msymGetAlignmentAxes(self.ctx, primary.as_mut_ptr(), secondary.as_mut_ptr())
        })?;
        Ok((primary, secondary))
    }

    /// Override the primary and secondary alignment axes. Call before `align_axes`.
    pub(crate) fn set_alignment_axes(
        &mut self,
        primary: [f64; 3],
        secondary: [f64; 3],
    ) -> Result<(), MsymError> {
        let mut primary = primary;
        let mut secondary = secondary;
        error::check(unsafe {
            ffi::msymSetAlignmentAxes(self.ctx, primary.as_mut_ptr(), secondary.as_mut_ptr())
        })
    }

    /// Return the 3×3 rotation matrix applied during `align_axes`.
    /// Requires `find_symmetry`.
    pub(crate) fn alignment_transform(&self) -> Result<[[f64; 3]; 3], MsymError> {
        let mut transform = [[0.0f64; 3]; 3];
        error::check(unsafe { ffi::msymGetAlignmentTransform(self.ctx, transform.as_mut_ptr()) })?;
        Ok(transform)
    }

    /// Override the 3×3 alignment rotation matrix. Call before `align_axes`.
    pub(crate) fn set_alignment_transform(
        &mut self,
        transform: [[f64; 3]; 3],
    ) -> Result<(), MsymError> {
        let mut transform = transform;
        error::check(unsafe { ffi::msymSetAlignmentTransform(self.ctx, transform.as_mut_ptr()) })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use crate::basis::{BasisFunction, BasisKind, CartesianAxis};

    use super::*;

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

    fn water_ctx() -> Context {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&water()).unwrap();
        ctx.find_symmetry().unwrap();
        ctx
    }

    fn displacement_basis(atom_count: usize) -> Vec<BasisFunction> {
        (0..atom_count)
            .flat_map(|atom| {
                [(CartesianAxis::X, 1), (CartesianAxis::Y, -1), (CartesianAxis::Z, 0)]
                    .into_iter()
                    .map(move |(axis, m)| BasisFunction {
                        atom_index: atom,
                        kind: BasisKind::Displacement(axis),
                        shell_index: 0,
                        l: 1,
                        m,
                    })
            })
            .collect()
    }

    fn water_salc_ctx() -> (Context, &'static PointGroup) {
        let mut ctx = water_ctx();
        let group = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        ctx.set_basis_functions(&displacement_basis(3)).unwrap();
        ctx.generate_salcs().unwrap();
        (ctx, group)
    }

    fn sf6() -> Vec<SymmetryCenter> {
        vec![
            SymmetryCenter { atomic_number: 16, mass: 32.06, position: [0.0, 0.0, 0.0], name: "S".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [1.564, 0.0, 0.0], name: "F".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [-1.564, 0.0, 0.0], name: "F".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [0.0, 1.564, 0.0], name: "F".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [0.0, -1.564, 0.0], name: "F".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [0.0, 0.0, 1.564], name: "F".into() },
            SymmetryCenter { atomic_number: 9, mass: 18.998, position: [0.0, 0.0, -1.564], name: "F".into() },
        ]
    }

    #[rstest]
    #[case(water(), SchoenfliesSymbol::Cnv(2), SchoenfliesSymbol::Cs, 2)]
    #[case(water(), SchoenfliesSymbol::Cnv(2), SchoenfliesSymbol::Cn(2), 1)]
    #[case(sf6(), SchoenfliesSymbol::Oh, SchoenfliesSymbol::Cn(2), 9)]
    #[case(sf6(), SchoenfliesSymbol::Oh, SchoenfliesSymbol::Ci, 1)]
    fn test_subgroup_multiplicity(
        #[case] centers: Vec<SymmetryCenter>,
        #[case] parent_sym: SchoenfliesSymbol,
        #[case] child_sym: SchoenfliesSymbol,
        #[case] expected: usize,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&centers).unwrap();
        ctx.find_symmetry().unwrap();
        let group = PointGroup::from_symbol(parent_sym).unwrap();
        let sgs = ctx.subgroups(group).unwrap();
        let sg = sgs.iter().find(|s| s.symbol() == child_sym).unwrap();
        assert_eq!(sg.multiplicity(), expected);
    }

    // --- Thresholds ---

    #[rstest]
    fn test_context_thresholds() {
        let mut ctx = Context::new().unwrap();
        let t = Thresholds::default();
        ctx.set_thresholds(&t).unwrap();
        let t2 = ctx.thresholds().unwrap();
        assert!((t.geometry - t2.geometry).abs() < 1e-15);
        assert!((t.angle - t2.angle).abs() < 1e-15);
    }

    // --- Centers ---

    #[rstest]
    fn test_context_centers() {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&water()).unwrap();
        let back = ctx.centers().unwrap();
        assert_eq!(back.len(), 3);
        assert_eq!(back[0].atomic_number, 8);
        assert_eq!(back[1].atomic_number, 1);
    }

    #[rstest]
    fn test_context_symmetrize_centers() {
        let mut ctx = water_ctx();
        let rms = ctx.symmetrize_centers().unwrap();
        assert!(rms < 1e-6);
    }

    #[rstest]
    fn test_context_apply_translation() {
        let mut ctx = water_ctx();
        let before = ctx.centers().unwrap();
        ctx.apply_translation(0, [0.0, 0.0, 0.1]).unwrap();
        let after = ctx.centers().unwrap();
        let delta_sq: f64 = (0..3)
            .map(|i| (after[0].position[i] - before[0].position[i]).powi(2))
            .sum();
        assert!(delta_sq.sqrt() > 0.05);
    }

    // --- Basis functions ---

    #[rstest]
    fn test_context_basis_functions() {
        let mut ctx = water_ctx();
        let basis = displacement_basis(3);
        ctx.set_basis_functions(&basis).unwrap();
        let back = ctx.basis_functions().unwrap();
        assert_eq!(back.len(), 9);
        for (orig, read) in basis.iter().zip(back.iter()) {
            assert_eq!(orig.atom_index, read.atom_index);
            assert_eq!(orig.l, read.l);
            assert_eq!(orig.m, read.m);
        }
    }

    // --- Point group ---

    #[rstest]
    #[case(water(), SchoenfliesSymbol::Cnv(2), "C2v")]
    #[case(methane(), SchoenfliesSymbol::Td, "Td")]
    fn test_context_find_symmetry(
        #[case] elements: Vec<SymmetryCenter>,
        #[case] expected_label: SchoenfliesSymbol,
        #[case] expected_name: &str,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&elements).unwrap();
        ctx.find_symmetry().unwrap();
        assert_eq!(ctx.point_group().unwrap(), expected_label);
        assert_eq!(ctx.point_group_name().unwrap(), expected_name);
    }

    #[rstest]
    fn test_context_set_point_group() {
        let mut ctx = Context::new().unwrap();
        ctx.set_point_group(SchoenfliesSymbol::Cnv(3)).unwrap();
        assert_eq!(ctx.point_group().unwrap(), SchoenfliesSymbol::Cnv(3));
        assert_eq!(ctx.point_group_name().unwrap(), "C3v");
    }

    #[rstest]
    fn test_context_set_point_group_by_name() {
        let mut ctx = Context::new().unwrap();
        ctx.set_point_group_by_name("D6h").unwrap();
        assert_eq!(ctx.point_group().unwrap(), SchoenfliesSymbol::Dnh(6));
    }

    #[rstest]
    fn test_context_subgroups() {
        let ctx = water_ctx();
        let group = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let sgs = ctx.subgroups(group).unwrap();
        assert!(!sgs.is_empty());
        let names: Vec<&str> = sgs.iter().map(|s| s.name()).collect();
        assert!(names.contains(&"C2"));
        assert!(names.contains(&"Cs"));
    }

    #[rstest]
    fn test_context_select_subgroup() {
        let mut ctx = water_ctx();
        let group = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();
        let sgs = ctx.subgroups(group).unwrap();
        let c2 = sgs.iter().find(|s| s.name() == "C2").unwrap();
        ctx.select_subgroup(c2).unwrap();
        assert_eq!(ctx.point_group_name().unwrap(), "C2");
    }

    #[rstest]
    fn test_subgroup_indices_independent_of_orientation() {
        let group = PointGroup::from_symbol(SchoenfliesSymbol::Cnv(2)).unwrap();

        // Water in xz plane.
        let water_xz = vec![
            SymmetryCenter { atomic_number: 8, mass: 15.999, position: [0.0, 0.0, 0.117], name: "O".into() },
            SymmetryCenter { atomic_number: 1, mass: 1.008, position: [0.757, 0.0, -0.469], name: "H".into() },
            SymmetryCenter { atomic_number: 1, mass: 1.008, position: [-0.757, 0.0, -0.469], name: "H".into() },
        ];
        // Water in yz plane.
        let water_yz = vec![
            SymmetryCenter { atomic_number: 8, mass: 15.999, position: [0.0, 0.0, 0.117], name: "O".into() },
            SymmetryCenter { atomic_number: 1, mass: 1.008, position: [0.0, 0.757, -0.469], name: "H".into() },
            SymmetryCenter { atomic_number: 1, mass: 1.008, position: [0.0, -0.757, -0.469], name: "H".into() },
        ];

        let mut ctx_xz = Context::new().unwrap();
        ctx_xz.set_centers(&water_xz).unwrap();
        ctx_xz.find_symmetry().unwrap();
        let sgs_xz = ctx_xz.subgroups(group).unwrap();

        let mut ctx_yz = Context::new().unwrap();
        ctx_yz.set_centers(&water_yz).unwrap();
        ctx_yz.find_symmetry().unwrap();
        let sgs_yz = ctx_yz.subgroups(group).unwrap();

        assert_eq!(sgs_xz.len(), sgs_yz.len());
        for (a, b) in sgs_xz.iter().zip(sgs_yz.iter()) {
            assert_eq!(a.symbol(), b.symbol());
            assert_eq!(a.name(), b.name());
            assert_eq!(a.order(), b.order());
        }

        let cs_xz = sgs_xz.iter().find(|s| s.symbol() == SchoenfliesSymbol::Cs).unwrap();
        let cs_yz = sgs_yz.iter().find(|s| s.symbol() == SchoenfliesSymbol::Cs).unwrap();
        assert_eq!(cs_xz.index(), cs_yz.index());
    }

    // --- Symmetry operations ---

    #[rstest]
    #[case(water(), SchoenfliesSymbol::Cnv(2), 4)]
    #[case(methane(), SchoenfliesSymbol::Td, 24)]
    fn test_context_symmetry_representation(
        #[case] elements: Vec<SymmetryCenter>,
        #[case] label: SchoenfliesSymbol,
        #[case] expected_count: usize,
    ) {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&elements).unwrap();
        ctx.find_symmetry().unwrap();
        let group = PointGroup::from_symbol(label).unwrap();
        let rep = ctx.symmetry_representation(group).unwrap();
        assert_eq!(rep.order(), expected_count);
    }

    // --- Equivalence sets ---

    #[rstest]
    fn test_context_equivalence_sets() {
        let ctx = water_ctx();
        let sets = ctx.equivalence_sets().unwrap();
        assert_eq!(sets.len(), 2);
        let mut sizes: Vec<usize> = sets.iter().map(|s| s.centers.len()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![1, 2]);
    }

    #[rstest]
    #[case(0, 1)]
    #[case(1, 2)]
    fn test_context_equivalence_set_for_center(
        #[case] index: usize,
        #[case] expected_size: usize,
    ) {
        let ctx = water_ctx();
        let es = ctx.equivalence_set_for_center(index).unwrap();
        assert_eq!(es.centers.len(), expected_size);
    }

    #[rstest]
    fn test_context_find_equivalence_sets() {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&water()).unwrap();
        ctx.set_point_group(SchoenfliesSymbol::Cnv(2)).unwrap();
        ctx.find_equivalence_sets().unwrap();
        let sets = ctx.equivalence_sets().unwrap();
        assert_eq!(sets.len(), 2);
    }

    #[rstest]
    fn test_context_find_equivalence_set_permutations() {
        let mut ctx = water_ctx();
        ctx.find_equivalence_set_permutations().unwrap();
    }

    // --- SALCs ---

    #[rstest]
    fn test_context_generate_salcs() {
        let mut ctx = water_ctx();
        ctx.set_basis_functions(&displacement_basis(3)).unwrap();
        ctx.generate_salcs().unwrap();
    }

    #[rstest]
    fn test_context_salcs() {
        let (ctx, _) = water_salc_ctx();
        let (coeffs, species) = ctx.salcs(9).unwrap();
        assert_eq!(coeffs.len(), 81);
        assert_eq!(species.len(), 9);
    }

    #[rstest]
    fn test_context_salcs_by_irrep() {
        let (ctx, group) = water_salc_ctx();
        let irrep_bases = ctx.salcs_by_irrep(group).unwrap();
        let total: usize = irrep_bases.iter().map(|ib| ib.salcs.len()).sum();
        assert_eq!(total, 9);
    }

    #[rstest]
    fn test_context_symmetrize_salcs() {
        let (mut ctx, _) = water_salc_ctx();
        let (mut coeffs, mut species) = ctx.salcs(9).unwrap();
        ctx.symmetrize_salcs(9, &mut coeffs, &mut species).unwrap();
        for &s in &species {
            assert!(s >= 0);
        }
    }

    #[rstest]
    fn test_context_irrep_components() {
        let (ctx, group) = water_salc_ctx();
        let irrep_count = group.irreps().len();
        let mut components = vec![0.0f64; irrep_count];
        let wf = vec![1.0f64; 9];
        ctx.irrep_components(&wf, &mut components).unwrap();
        let total: f64 = components.iter().sum();
        assert!(total > 0.0);
    }

    // --- Geometry and alignment ---

    #[rstest]
    fn test_context_center_of_mass() {
        let ctx = water_ctx();
        let com = ctx.center_of_mass().unwrap();
        assert!(com[0].abs() < 0.5);
        assert!(com[1].abs() < 0.5);
    }

    #[rstest]
    fn test_context_set_center_of_mass() {
        let mut ctx = Context::new().unwrap();
        ctx.set_centers(&water()).unwrap();
        ctx.set_center_of_mass([1.0, 2.0, 3.0]).unwrap();
        let com = ctx.center_of_mass().unwrap();
        assert!((com[0] - 1.0).abs() < 1e-10);
        assert!((com[1] - 2.0).abs() < 1e-10);
        assert!((com[2] - 3.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_context_radius() {
        let ctx = water_ctx();
        let r = ctx.radius().unwrap();
        assert!(r > 0.0);
    }

    #[rstest]
    fn test_context_molecular_shape() {
        let ctx = water_ctx();
        let shape = ctx.molecular_shape().unwrap();
        assert_ne!(shape, MolecularShape::Linear);
        assert_ne!(shape, MolecularShape::Spherical);
    }

    #[rstest]
    fn test_context_principal_moments() {
        let ctx = water_ctx();
        let moments = ctx.principal_moments().unwrap();
        for m in &moments {
            assert!(*m >= 0.0);
        }
    }

    #[rstest]
    fn test_context_principal_axes() {
        let ctx = water_ctx();
        let axes = ctx.principal_axes().unwrap();
        for row in &axes {
            let norm_sq: f64 = row.iter().map(|x| x * x).sum();
            assert!((norm_sq - 1.0).abs() < 1e-10);
        }
    }

    #[rstest]
    fn test_context_align_axes() {
        let mut ctx = water_ctx();
        ctx.align_axes().unwrap();
    }

    #[rstest]
    fn test_context_alignment_axes() {
        let mut ctx = water_ctx();
        let (primary, secondary) = ctx.alignment_axes().unwrap();
        ctx.set_alignment_axes(primary, secondary).unwrap();
        let (p2, s2) = ctx.alignment_axes().unwrap();
        for i in 0..3 {
            assert!((p2[i] - primary[i]).abs() < 1e-10);
            assert!((s2[i] - secondary[i]).abs() < 1e-10);
        }
    }

    #[rstest]
    fn test_context_alignment_transform() {
        let ctx = water_ctx();
        let transform = ctx.alignment_transform().unwrap();
        for row in &transform {
            let norm_sq: f64 = row.iter().map(|x| x * x).sum();
            assert!((norm_sq - 1.0).abs() < 1e-10);
        }
    }

    #[rstest]
    fn test_context_set_alignment_transform() {
        let mut ctx = water_ctx();
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        ctx.set_alignment_transform(identity).unwrap();
        let t = ctx.alignment_transform().unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - identity[i][j]).abs() < 1e-10);
            }
        }
    }
}
