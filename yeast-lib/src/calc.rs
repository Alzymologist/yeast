//! Data analysis tools
//!
//! # Yeast growth model
//!
//! Yeast concentration vs time growth curve is described as logistic function:
//!
//!
//! # Hypothesis 0
//!
//!
//!
//! # Assumptions
//!
//!

use multicalc::numerical_integration::{
    integrator::IntegratorSingleVariable, iterative_integration::SingleVariableSolver,
};
use probability::distribution::{Continuous, Gamma, Gaussian};

/// Probability to encounter given ln_n1 parameter value (logarithm of the final
/// cell councentration, CFU/m^3).
///
/// Using prior assumptions.
pub fn ln_n1_probability(ln_n1: f64) -> f64 {
    Gaussian::new(LN_N1_MEAN_PRIOR, LN_N1_VARIANCE_PRIOR).density(ln_n1)
}

/// Logarithm of the final cell councentration, CFU/m^3, prior assumption for
/// mean
pub const LN_N1_MEAN_PRIOR: f64 = 32.4;

/// Logarithm of the final cell councentration, CFU/m^3, prior assumption for
/// variance
pub const LN_N1_VARIANCE_PRIOR: f64 = 1.2;

pub const T_BEND_MEAN_PRIOR: f64 = 15.0 * 60.0;

pub const T_BEND_VARIANCE_PRIOR: f64 = 15.0 * 60.0;

pub fn t_bend_probability(t_bend: f64) -> f64 {
    Gaussian::new(T_BEND_MEAN_PRIOR, T_BEND_VARIANCE_PRIOR).density(t_bend)
}

pub const T_DOUBLING_MEAN_PRIOR: f64 = 3.0 * 60.0;

pub const T_DOUBLING_VARIANCE_PRIOR: f64 = 0.9 * 60.0;

pub fn t_doubling_probability(t_doubling: f64) -> f64 {
    Gaussian::new(T_DOUBLING_MEAN_PRIOR, T_DOUBLING_VARIANCE_PRIOR).density(t_doubling)
}

pub fn log_logistics(t: f64, ln_n1: f64, t_bend: f64, t_doubling: f64) -> f64 {
    ln_n1 - (1.0 + (-1.0 * (t - t_bend) / t_doubling).exp()).ln()
}

pub fn pr_log_logistics(
    t: f64,
    ln_n1: f64,
    t_bend: f64,
    t_doubling: f64,
    sigma: f64,
    y: f64,
    y_error: f64,
) -> f64 {
    let mu = log_logistics(t, ln_n1, t_bend, t_doubling);
    let gaussian_logistics = Gaussian::new(mu, sigma);
    let gaussian_y = Gaussian::new(y, y_error);
    let func = |arg: f64| -> f64 { gaussian_logistics.density(arg.ln()) * gaussian_y.density(arg) };

    let integrator = SingleVariableSolver::default();
    let integration_limit = [y - 3.0 * y_error, y + 3.0 * y_error];
    integrator.get_single(&func, &integration_limit).unwrap()
}

pub const PRECISION_SHAPE: f64 = 18.0;
pub const PRECISION_SCALE: f64 = 0.4;

pub fn precision_probability(precision: f64) -> f64 {
    Gamma::new(PRECISION_SHAPE, PRECISION_SCALE).density(precision)
}

/*
Uniformity procedure.
Detecting the situation then single uniformity run (typically consisting of 3-4 samples)

QA procedure.
Detecting the situation when yeast culture changed its properties (could be mutation, contamination, accidental sample swapping etc).
Samples after certain date and before certain date group in two distinct data clusters better than in a single one.
*/
