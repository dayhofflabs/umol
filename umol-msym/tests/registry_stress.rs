use std::thread;

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use umol_msym::PointGroup;

fn verify_group(pg: &'static PointGroup, name: &str, order: usize) {
    assert_eq!(pg.order(), order, "{name}: wrong order");

    let irreps = pg.irreps();
    assert!(!irreps.is_empty(), "{name}: no irreps");

    // Sum of dim² = group order (complex merged pairs contribute 2×1² not 1×2²)
    let dim_sq_sum: usize = irreps
        .iter()
        .map(|ir| {
            if ir.complex() {
                2
            } else {
                (ir.dimension() as usize).pow(2)
            }
        })
        .sum();
    assert_eq!(dim_sq_sum, order, "{name}: sum(dim²) != order");

    // Number of irreps = number of classes = number of class sizes
    assert_eq!(
        irreps.len(),
        pg.class_sizes().len(),
        "{name}: irrep count != class count"
    );

    // Totally symmetric irrep
    let ts = pg.totally_symmetric_irrep();
    assert!(ts.totally_symmetric(), "{name}: totally_symmetric() is false");

    // Row orthogonality
    let h = order as f64;
    let class_reps = pg.class_reps();
    for ir in &irreps {
        let norm_sq: f64 = class_reps
            .iter()
            .zip(pg.class_sizes())
            .map(|(op, &size)| {
                let chi = op.character(*ir);
                size as f64 * chi * chi
            })
            .sum();
        let expected = if ir.complex() { 2.0 * h } else { h };
        assert!(
            (norm_sq - expected).abs() < 0.01,
            "{name}: irrep {} norm² = {norm_sq}, expected {expected}",
            ir.symbol()
        );
    }
}

/// Parametric group families with order formulas.
const FAMILIES: &[(&str, fn(u32) -> (String, usize))] = &[
    ("Cn", |n| (format!("C{n}"), n as usize)),
    ("Cnv", |n| (format!("C{n}v"), 2 * n as usize)),
    ("Cnh", |n| (format!("C{n}h"), 2 * n as usize)),
    ("Dn", |n| (format!("D{n}"), 2 * n as usize)),
    ("Dnh", |n| (format!("D{n}h"), 4 * n as usize)),
    ("Dnd", |n| (format!("D{n}d"), 4 * n as usize)),
    ("Sn", |n| (format!("S{n}"), n as usize)),
];

/// Threads race to insert groups into the registry.
/// Each thread shuffles the same group list differently, maximizing insert contention.
#[test]
fn test_registry_concurrent_insertion() {
    let thread_count = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);

    let mut groups: Vec<(String, usize)> = vec![
        ("C1".into(), 1),
        ("Ci".into(), 2),
        ("Cs".into(), 2),
        ("T".into(), 12),
        ("Td".into(), 24),
        ("Th".into(), 24),
        ("O".into(), 24),
        ("Oh".into(), 48),
        ("I".into(), 60),
        ("Ih".into(), 120),
    ];
    for n in 2..=12u32 {
        for &(family, order_fn) in FAMILIES {
            let actual_n = if family == "Sn" { 2 * n } else { n };
            let (name, order) = order_fn(actual_n);
            groups.push((name, order));
        }
    }

    let handles: Vec<_> = (0..thread_count)
        .map(|t| {
            let mut groups = groups.clone();
            thread::spawn(move || {
                let mut rng = ChaCha8Rng::seed_from_u64(t as u64);
                groups.shuffle(&mut rng);
                for (name, order) in &groups {
                    let pg = PointGroup::parse(name).unwrap();
                    verify_group(pg, name, *order);
                }
            })
        })
        .collect();

    for h in handles {
        h.join().unwrap();
    }
}

/// Round-robin retrieval of common groups under contention.
#[test]
fn test_registry_concurrent_retrieval() {
    let thread_count = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let iterations = 1000;

    let groups: Vec<(&str, usize)> = vec![
        ("C1", 1),
        ("Ci", 2),
        ("Cs", 2),
        ("C2v", 4),
        ("C3v", 6),
        ("D2h", 8),
        ("D6h", 24),
        ("T", 12),
        ("Td", 24),
        ("Oh", 48),
        ("I", 60),
        ("Ih", 120),
    ];

    for &(name, _) in &groups {
        PointGroup::parse(name).unwrap();
    }

    let handles: Vec<_> = (0..thread_count)
        .map(|t| {
            let mut groups = groups.clone();
            thread::spawn(move || {
                let mut rng = ChaCha8Rng::seed_from_u64(t as u64);
                for _ in 0..iterations {
                    groups.shuffle(&mut rng);
                    let (name, order) = groups[0];
                    let pg = PointGroup::parse(name).unwrap();
                    verify_group(pg, name, order);
                }
            })
        })
        .collect();

    for h in handles {
        h.join().unwrap();
    }
}
