use nalgebra::{DMatrix, SymmetricEigen};

// Minimal Hückel MO computation utility over adjacency matrix
pub fn hmo_density_from_adjacency(edges: &[(usize, usize)], size: usize) -> DMatrix<f64> {
    if size == 0 {
        return DMatrix::<f64>::zeros(0, 0);
    }
    let mut a = DMatrix::<f64>::zeros(size, size);
    for &(u, v) in edges {
        if u < size && v < size && u != v {
            a[(u, v)] = 1.0;
            a[(v, u)] = 1.0;
        }
    }
    let se = SymmetricEigen::new(a);
    let q = se.eigenvectors; // columns = eigenvectors
    let n = size;
    let n_occ = (n + 1) / 2; // placeholder closed-shell
    let mut p = DMatrix::<f64>::zeros(n, n);
    for k in 0..n_occ {
        let col = q.column(k);
        for i in 0..n {
            for j in 0..n {
                p[(i, j)] += 2.0 * col[i] * col[j];
            }
        }
    }
    p
}


