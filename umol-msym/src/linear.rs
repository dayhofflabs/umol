use crate::basis::{BasisFunction, IrrepBasis, Salc, SalcBasis};
use crate::point_group::{Irrep, PointGroup};
use crate::types::{SchoenfliesLabel, SymmetryCenter};

fn is_dooh(group: &PointGroup) -> bool {
    group.label() == SchoenfliesLabel::Dooh
}

fn linear_info(irrep: Irrep) -> (u32, Option<bool>, Option<bool>) {
    (
        irrep.lambda().expect("linear_info called on finite group irrep"),
        irrep.data.sigma_v,
        irrep.data.gerade,
    )
}

fn find_irrep(
    group: &'static PointGroup,
    lambda: u32,
    sigma_v: Option<bool>,
    gerade: Option<bool>,
) -> Option<Irrep> {
    group.irreps().into_iter().find(|ir| {
        let (l, sv, g) = linear_info(*ir);
        l == lambda && sv == sigma_v && g == gerade
    })
}

/// g⊗g=g, u⊗u=g, g⊗u=u
fn gu_product(g1: Option<bool>, g2: Option<bool>) -> Option<bool> {
    match (g1, g2) {
        (Some(a), Some(b)) => Some(a == b),
        (None, None) => None,
        _ => unreachable!(),
    }
}

pub fn direct_product(
    group: &'static PointGroup,
    a: Irrep,
    b: Irrep,
) -> Vec<(Irrep, u32)> {
    let (la, sva, ga) = linear_info(a);
    let (lb, svb, gb) = linear_info(b);
    let gu = gu_product(ga, gb);

    let mut result = Vec::new();

    match (la, lb) {
        (0, 0) => {
            // Σ ⊗ Σ: +⊗+=+, -⊗-=+, +⊗-=-
            let sv_product = Some(sva.unwrap() == svb.unwrap());
            if let Some(ir) = find_irrep(group, 0, sv_product, gu) {
                result.push((ir, 1));
            }
        }
        (0, _) => {
            // Σ± ⊗ Λ = Λ
            if let Some(ir) = find_irrep(group, lb, None, gu) {
                result.push((ir, 1));
            }
        }
        (_, 0) => {
            // Λ ⊗ Σ± = Λ
            if let Some(ir) = find_irrep(group, la, None, gu) {
                result.push((ir, 1));
            }
        }
        _ if la == lb => {
            // Λ ⊗ Λ (same λ) = Σ+ + Σ- + Λ_{2λ}
            if let Some(ir) = find_irrep(group, 0, Some(true), gu) {
                result.push((ir, 1));
            }
            if let Some(ir) = find_irrep(group, 0, Some(false), gu) {
                result.push((ir, 1));
            }
            if let Some(ir) = find_irrep(group, 2 * la, None, gu) {
                result.push((ir, 1));
            }
        }
        _ => {
            // Λ₁ ⊗ Λ₂ (different λ) = |λ₁-λ₂| + (λ₁+λ₂)
            let diff = la.abs_diff(lb);
            let sum = la + lb;
            // diff > 0 here since la ≠ lb and both > 0
            if let Some(ir) = find_irrep(group, diff, None, gu) {
                result.push((ir, 1));
            }
            if let Some(ir) = find_irrep(group, sum, None, gu) {
                result.push((ir, 1));
            }
        }
    }

    result
}

/// [a²]: For 1D irreps (Σ) → Σ+. For 2D irreps (λ>0) → Σ+ + irrep(2λ).
pub fn symmetric_square(group: &'static PointGroup, a: Irrep) -> Vec<(Irrep, u32)> {
    let (la, _sva, ga) = linear_info(a);
    let gu = gu_product(ga, ga);
    let mut result = Vec::new();
    if let Some(ir) = find_irrep(group, 0, Some(true), gu) {
        result.push((ir, 1));
    }
    if la > 0 {
        if let Some(ir) = find_irrep(group, 2 * la, None, gu) {
            result.push((ir, 1));
        }
    }
    result
}

/// {a²}: For 1D irreps (Σ) → empty. For 2D irreps (λ>0) → Σ-.
pub fn antisymmetric_square(group: &'static PointGroup, a: Irrep) -> Vec<(Irrep, u32)> {
    let (la, _sva, ga) = linear_info(a);
    if la == 0 {
        return Vec::new();
    }
    let gu = gu_product(ga, ga);
    let mut result = Vec::new();
    if let Some(ir) = find_irrep(group, 0, Some(false), gu) {
        result.push((ir, 1));
    }
    result
}

/// C∞v: Σ+ + Π. D∞h: Σ+u + Πu.
pub fn translation_irreps(group: &'static PointGroup) -> Vec<(Irrep, u32)> {
    let gu = if is_dooh(group) { Some(false) } else { None };
    let mut result = Vec::new();
    if let Some(ir) = find_irrep(group, 0, Some(true), gu) {
        result.push((ir, 1));
    }
    if let Some(ir) = find_irrep(group, 1, None, gu) {
        result.push((ir, 1));
    }
    result
}

/// C∞v: Π. D∞h: Πg. (Linear molecules: 2 rotational DOF, no Rz.)
pub fn rotation_irreps(group: &'static PointGroup) -> Vec<(Irrep, u32)> {
    let gu = if is_dooh(group) { Some(true) } else { None };
    let mut result = Vec::new();
    if let Some(ir) = find_irrep(group, 1, None, gu) {
        result.push((ir, 1));
    }
    result
}

/// Sym²(translation): C∞v: 2Σ+ + Π + Δ. D∞h: 2Σ+g + Πg + Δg.
pub fn quadratic_irreps(group: &'static PointGroup) -> Vec<(Irrep, u32)> {
    let gu = if is_dooh(group) { Some(true) } else { None };
    let mut result = Vec::new();
    if let Some(ir) = find_irrep(group, 0, Some(true), gu) {
        result.push((ir, 2));
    }
    if let Some(ir) = find_irrep(group, 1, None, gu) {
        result.push((ir, 1));
    }
    if let Some(ir) = find_irrep(group, 2, None, gu) {
        result.push((ir, 1));
    }
    result
}

/// a ⊗ b ⊗ c ⊃ Σ+ (or Σ+g)?
pub fn contains_totally_symmetric(
    group: &'static PointGroup,
    a: Irrep,
    b: Irrep,
    c: Irrep,
) -> bool {
    let ab = direct_product(group, a, b);
    ab.iter().any(|(ir, _)| *ir == c)
}

#[derive(Clone, Copy)]
enum AtomRole {
    Singleton,
    Positive(usize),
    Negative,
}

/// For D∞h: pair atoms at ±z. Returns (pairs, singletons).
fn find_pairs(centers: &[SymmetryCenter], equivalence: f64) -> (Vec<(usize, usize)>, Vec<usize>) {
    let n = centers.len();
    let mut used = vec![false; n];
    let mut pairs = Vec::new();
    let mut singletons = Vec::new();

    for i in 0..n {
        if used[i] {
            continue;
        }
        let zi = centers[i].position[2];

        if zi.abs() < equivalence {
            singletons.push(i);
            used[i] = true;
            continue;
        }

        let mut found = false;
        for j in (i + 1)..n {
            if used[j] || centers[j].atomic_number != centers[i].atomic_number {
                continue;
            }
            if (zi + centers[j].position[2]).abs() < equivalence {
                let (p, q) = if zi > 0.0 { (i, j) } else { (j, i) };
                pairs.push((p, q));
                used[i] = true;
                used[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            singletons.push(i);
            used[i] = true;
        }
    }

    (pairs, singletons)
}

pub fn compute_salcs(
    centers: &[SymmetryCenter],
    basis: &[BasisFunction],
    group: &'static PointGroup,
    equivalence: f64,
) -> SalcBasis {
    let dooh = is_dooh(group);
    let mut irrep_salcs: Vec<(Irrep, Vec<Salc>)> = Vec::new();

    let mut add_salc = |irrep: Irrep, salc: Salc| {
        if let Some(entry) = irrep_salcs.iter_mut().find(|(ir, _)| *ir == irrep) {
            entry.1.push(salc);
        } else {
            irrep_salcs.push((irrep, vec![salc]));
        }
    };

    if !dooh {
        // C∞v: each basis function is its own SALC
        for (i, bf) in basis.iter().enumerate() {
            let lambda = bf.m.unsigned_abs();
            let sigma_v = if lambda == 0 { Some(true) } else { None };
            let irrep = find_irrep(group, lambda, sigma_v, None)
                .unwrap_or_else(|| panic!("no irrep for λ={lambda}"));
            add_salc(irrep, Salc {
                coefficients: vec![(i, 1.0)],
            });
        }
    } else {
        // D∞h: form ±1/√2 combinations for paired atoms
        let (pairs, _singletons) = find_pairs(centers, equivalence);
        let inv2 = 1.0 / 2.0_f64.sqrt();

        let mut atom_role = vec![AtomRole::Singleton; centers.len()];
        for &(p, q) in &pairs {
            atom_role[p] = AtomRole::Positive(q);
            atom_role[q] = AtomRole::Negative;
        }

        let mut processed = vec![false; basis.len()];

        for (i, bf) in basis.iter().enumerate() {
            if processed[i] {
                continue;
            }

            let lambda = bf.m.unsigned_abs();
            let l_even = bf.l % 2 == 0;

            match atom_role[bf.atom_index] {
                AtomRole::Singleton => {
                    // Central atom: gerade if l even, ungerade if l odd
                    let gerade = Some(l_even);
                    let sigma_v = if lambda == 0 { Some(true) } else { None };
                    let irrep = find_irrep(group, lambda, sigma_v, gerade)
                        .unwrap_or_else(|| panic!("no irrep for λ={lambda}, g={gerade:?}"));
                    add_salc(irrep, Salc {
                        coefficients: vec![(i, 1.0)],
                    });
                    processed[i] = true;
                }
                AtomRole::Positive(partner_atom) => {
                    // Find corresponding basis function on partner atom
                    let j = basis
                        .iter()
                        .position(|bj| {
                            bj.atom_index == partner_atom
                                && bj.l == bf.l
                                && bj.m == bf.m
                                && bj.kind == bf.kind
                        })
                        .expect("partner basis function not found");

                    // Symmetric combination: (f_+ + f_-)/√2
                    // Gerade if l even (parity of spherical harmonic under inversion is (-1)^l)
                    let sym_gerade = Some(l_even);
                    let sigma_v = if lambda == 0 { Some(true) } else { None };
                    let sym_irrep = find_irrep(group, lambda, sigma_v, sym_gerade)
                        .unwrap_or_else(|| panic!("no irrep for λ={lambda}, g={sym_gerade:?}"));
                    add_salc(sym_irrep, Salc {
                        coefficients: vec![(i, inv2), (j, inv2)],
                    });

                    // Antisymmetric combination: (f_+ - f_-)/√2
                    let antisym_gerade = Some(!l_even);
                    let antisym_irrep =
                        find_irrep(group, lambda, sigma_v, antisym_gerade).unwrap_or_else(|| {
                            panic!("no irrep for λ={lambda}, g={antisym_gerade:?}")
                        });
                    add_salc(antisym_irrep, Salc {
                        coefficients: vec![(i, inv2), (j, -inv2)],
                    });

                    processed[i] = true;
                    processed[j] = true;
                }
                AtomRole::Negative => {
                    // Handled when positive partner is processed
                    continue;
                }
            }
        }
    }

    let irrep_bases: Vec<IrrepBasis> = irrep_salcs
        .into_iter()
        .map(|(irrep, salcs)| IrrepBasis { irrep, salcs })
        .collect();

    SalcBasis {
        basis_functions: basis.to_vec(),
        irreps: irrep_bases,
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use crate::point_group::PointGroup;

    #[rstest]
    fn test_coov_irreps() {
        let g = PointGroup::coov();
        assert_eq!(g.to_string(), "C∞v");
        assert_eq!(g.order(), 0);
        assert!(g.is_linear());
        assert_eq!(g.irreps().len(), 8);
        assert!(g.irrep("Σ+").is_some());
        assert!(g.irrep("Σ-").is_some());
        assert!(g.irrep("Π").is_some());
        assert!(g.irrep("Δ").is_some());
    }

    #[rstest]
    fn test_dooh_irreps() {
        let g = PointGroup::dooh();
        assert_eq!(g.to_string(), "D∞h");
        assert_eq!(g.order(), 0);
        assert!(g.is_linear());
        assert_eq!(g.irreps().len(), 16);
        assert!(g.irrep("Σ+g").is_some());
        assert!(g.irrep("Σ-u").is_some());
        assert!(g.irrep("Πg").is_some());
        assert!(g.irrep("Πu").is_some());
        assert!(g.irrep("Δg").is_some());
    }

    #[rstest]
    #[case("Σ+", "Σ+", &[("Σ+", 1)])]
    #[case("Σ+", "Σ-", &[("Σ-", 1)])]
    #[case("Σ-", "Σ-", &[("Σ+", 1)])]
    #[case("Σ+", "Π", &[("Π", 1)])]
    #[case("Σ-", "Π", &[("Π", 1)])]
    #[case("Π", "Π", &[("Σ+", 1), ("Σ-", 1), ("Δ", 1)])]
    #[case("Π", "Δ", &[("Π", 1), ("Φ", 1)])]
    #[case("Δ", "Δ", &[("Σ+", 1), ("Σ-", 1), ("Γ", 1)])]
    fn test_coov_direct_product(
        #[case] a: &str,
        #[case] b: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::coov();
        let ia = g.irrep(a).unwrap();
        let ib = g.irrep(b).unwrap();
        let result = g.direct_product(ia, ib);
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case("Πg", "Πg", &[("Σ+g", 1), ("Σ-g", 1), ("Δg", 1)])]
    #[case("Πu", "Πu", &[("Σ+g", 1), ("Σ-g", 1), ("Δg", 1)])]
    #[case("Πg", "Πu", &[("Σ+u", 1), ("Σ-u", 1), ("Δu", 1)])]
    #[case("Σ+u", "Πg", &[("Πu", 1)])]
    fn test_dooh_direct_product(
        #[case] a: &str,
        #[case] b: &str,
        #[case] expected: &[(&str, u32)],
    ) {
        let g = PointGroup::dooh();
        let ia = g.irrep(a).unwrap();
        let ib = g.irrep(b).unwrap();
        let result = g.direct_product(ia, ib);
        let actual: Vec<(&str, u32)> = result.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_coov_translation_rotation() {
        let g = PointGroup::coov();
        let ti = g.translation_irreps();
        let trans: Vec<(&str, u32)> = ti.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(trans, vec![("Σ+", 1), ("Π", 1)]);

        let ri = g.rotation_irreps();
        let rot: Vec<(&str, u32)> = ri.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(rot, vec![("Π", 1)]);
    }

    #[rstest]
    fn test_dooh_translation_rotation() {
        let g = PointGroup::dooh();
        let ti = g.translation_irreps();
        let trans: Vec<(&str, u32)> = ti.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(trans, vec![("Σ+u", 1), ("Πu", 1)]);

        let ri = g.rotation_irreps();
        let rot: Vec<(&str, u32)> = ri.iter().map(|(ir, n)| (ir.symbol(), *n)).collect();
        assert_eq!(rot, vec![("Πg", 1)]);
    }

    #[rstest]
    fn test_coov_selection_rules() {
        let g = PointGroup::coov();
        let sp = g.irrep("Σ+").unwrap();
        let sm = g.irrep("Σ-").unwrap();
        let pi = g.irrep("Π").unwrap();
        let de = g.irrep("Δ").unwrap();

        // Electric dipole: Σ+→Σ+ (via Σ+), Σ+→Π (via Π), Σ+→Δ forbidden
        assert!(g.electric_dipole_allowed(sp, sp));
        assert!(g.electric_dipole_allowed(sp, pi));
        assert!(!g.electric_dipole_allowed(sp, de));
        assert!(!g.electric_dipole_allowed(sp, sm));

        // Π→Π allowed (via Σ+)
        assert!(g.electric_dipole_allowed(pi, pi));
    }

    #[rstest]
    fn test_dooh_mutual_exclusion() {
        let g = PointGroup::dooh();
        let spg = g.irrep("Σ+g").unwrap();

        // D∞h has mutual exclusion: IR-active modes are Raman-inactive and vice versa
        // IR active: Σ+u, Πu (ungerade)
        // Raman active: Σ+g, Πg, Δg (gerade)
        let spu = g.irrep("Σ+u").unwrap();
        let pig = g.irrep("Πg").unwrap();
        let piu = g.irrep("Πu").unwrap();

        assert!(g.electric_dipole_allowed(spg, spu)); // IR active
        assert!(g.electric_dipole_allowed(spg, piu)); // IR active
        assert!(!g.electric_dipole_allowed(spg, pig)); // IR forbidden (gerade)

        assert!(g.raman_allowed(spg, pig)); // Raman active
        assert!(!g.raman_allowed(spg, piu)); // Raman forbidden (ungerade)
    }
}
