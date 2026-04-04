//! Discrete optimization with linear equality constraints via Lagrangian relaxation.
//!
//! Maximizes Σ log f_i(x_i) subject to linear equality constraints Σ_k x_{jk} = v_j,
//! where each x_i takes integer values in a finite domain.

/// A variable with a finite set of log-likelihood values for each possible assignment.
pub struct Variable {
    /// Log-likelihood for each possible value (indexed by value).
    pub log_likelihoods: Vec<f64>,
    /// Constraint indices this variable participates in (index, coefficient).
    pub constraints: Vec<(usize, f64)>,
}

/// A linear equality constraint: Σ coeff_k * x_k = rhs.
pub struct Constraint {
    pub rhs: f64,
}

/// Configuration for the Lagrangian relaxation solver.
pub struct LagrangianConfig {
    pub max_iter: usize,
    pub step_scale: f64,
}

impl Default for LagrangianConfig {
    fn default() -> Self {
        Self {
            max_iter: 200,
            step_scale: 0.5,
        }
    }
}

/// Result of Lagrangian relaxation.
pub struct LagrangianResult {
    /// Optimal assignment for each variable.
    pub assignments: Vec<usize>,
    /// Final Lagrange multipliers.
    pub multipliers: Vec<f64>,
    /// Whether all constraints are satisfied.
    pub feasible: bool,
    /// Constraint residuals (actual - target).
    pub residuals: Vec<f64>,
}

/// Solve a discrete optimization problem with linear equality constraints.
///
/// Each variable x_i has a domain of size |log_likelihoods_i|. The objective
/// is to maximize Σ log_lik_i(x_i) subject to linear constraints.
///
/// If the Lagrangian solution is infeasible (degenerate case), falls back
/// to greedy primal recovery.
pub fn lagrangian_relaxation(
    variables: &[Variable],
    constraints: &[Constraint],
    config: &LagrangianConfig,
) -> LagrangianResult {
    let m = constraints.len();
    let mut lambda = vec![0.0_f64; m];
    let mut assignments = vec![0usize; variables.len()];

    for iter in 0..config.max_iter {
        let step = config.step_scale / (1.0 + iter as f64);

        // Solve subproblems: for each variable, pick the value maximizing
        // log f(x) - Σ_j lambda_j * coeff_j * x
        for (i, var) in variables.iter().enumerate() {
            let mut best_k = 0;
            let mut best_val = f64::NEG_INFINITY;
            let dual_cost: f64 = var.constraints.iter().map(|&(j, c)| lambda[j] * c).sum();
            for (k, &ll) in var.log_likelihoods.iter().enumerate() {
                let val = ll - dual_cost * k as f64;
                if val > best_val {
                    best_val = val;
                    best_k = k;
                }
            }
            assignments[i] = best_k;
        }

        // Compute subgradient and update multipliers
        let mut residuals = vec![0.0_f64; m];
        for (i, var) in variables.iter().enumerate() {
            for &(j, coeff) in &var.constraints {
                residuals[j] += coeff * assignments[i] as f64;
            }
        }

        let mut max_violation = 0.0_f64;
        for j in 0..m {
            residuals[j] -= constraints[j].rhs;
            max_violation = max_violation.max(residuals[j].abs());
            lambda[j] += step * residuals[j];
        }

        if max_violation < 0.5 {
            break;
        }
    }

    // Check feasibility
    let mut residuals = vec![0.0_f64; m];
    for (i, var) in variables.iter().enumerate() {
        for &(j, coeff) in &var.constraints {
            residuals[j] += coeff * assignments[i] as f64;
        }
    }
    for j in 0..m {
        residuals[j] -= constraints[j].rhs;
    }
    let feasible = residuals.iter().all(|&r| r.abs() < 0.5);

    if !feasible {
        greedy_recover(variables, constraints, &mut assignments);
        // Recompute residuals
        for r in residuals.iter_mut() {
            *r = 0.0;
        }
        for (i, var) in variables.iter().enumerate() {
            for &(j, coeff) in &var.constraints {
                residuals[j] += coeff * assignments[i] as f64;
            }
        }
        for j in 0..m {
            residuals[j] -= constraints[j].rhs;
        }
    }
    let feasible = residuals.iter().all(|&r| r.abs() < 0.5);

    LagrangianResult {
        assignments,
        multipliers: lambda,
        feasible,
        residuals,
    }
}

/// Greedy primal recovery for degenerate cases.
///
/// Starts from all-zero assignments and greedily increases variable values
/// by the best marginal log-likelihood gain, subject to constraint capacity.
fn greedy_recover(
    variables: &[Variable],
    constraints: &[Constraint],
    assignments: &mut [usize],
) {
    let m = constraints.len();
    let mut remaining = vec![0.0_f64; m];
    for j in 0..m {
        remaining[j] = constraints[j].rhs;
    }

    for a in assignments.iter_mut() {
        *a = 0;
    }

    loop {
        let mut best_idx = None;
        let mut best_score = f64::NEG_INFINITY;

        for (i, var) in variables.iter().enumerate() {
            let current = assignments[i];
            if current + 1 >= var.log_likelihoods.len() {
                continue;
            }
            // Check that all constraints have remaining capacity
            let can_increase = var.constraints.iter().all(|&(j, coeff)| {
                if coeff > 0.0 { remaining[j] >= coeff } else { true }
            });
            if !can_increase {
                continue;
            }
            let score = var.log_likelihoods[current + 1] - var.log_likelihoods[current];
            if score > best_score {
                best_score = score;
                best_idx = Some(i);
            }
        }

        match best_idx {
            Some(i) => {
                assignments[i] += 1;
                for &(j, coeff) in &variables[i].constraints {
                    remaining[j] -= coeff;
                }
            }
            None => break,
        }

        if remaining.iter().all(|&r| r <= 0.0) {
            break;
        }
    }
}
