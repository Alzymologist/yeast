use serde::Serialize;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Uncertain {
    pub value: f64,
    pub error: f64,
}

impl Add for Uncertain {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            value: self.value + other.value,
            error: (self.error.powi(2) + other.error.powi(2)).sqrt(),
        }
    }
}

impl Add<f64> for Uncertain {
    type Output = Self;

    fn add(self, other: f64) -> Self {
        Self {
            value: self.value + other,
            error: self.error,
        }
    }
}

impl Add<Uncertain> for f64 {
    type Output = Uncertain;

    fn add(self, other: Uncertain) -> Uncertain {
        Uncertain {
            value: self + other.value,
            error: other.error,
        }
    }
}

impl Sub for Uncertain {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value - other.value,
            error: (self.error.powi(2) + other.error.powi(2)).sqrt(),
        }
    }
}

impl Sub<f64> for Uncertain {
    type Output = Self;

    fn sub(self, other: f64) -> Self {
        Self {
            value: self.value - other,
            error: self.error,
        }
    }
}

impl Sub<Uncertain> for f64 {
    type Output = Uncertain;

    fn sub(self, other: Uncertain) -> Uncertain {
        Uncertain {
            value: self - other.value,
            error: other.error,
        }
    }
}

impl Mul<f64> for Uncertain {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self {
            value: self.value * rhs,
            error: self.error * rhs,
        }
    }
}

impl Mul<Uncertain> for f64 {
    type Output = Uncertain;

    fn mul(self, rhs: Uncertain) -> Uncertain {
        Uncertain {
            value: self * rhs.value,
            error: self * rhs.error,
        }
    }
}

impl Mul<Uncertain> for Uncertain {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            error: ((self.error * rhs.value).powi(2) + (self.value * rhs.error).powi(2)).sqrt(),
        }
    }
}

impl Div<f64> for Uncertain {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            value: self.value / rhs,
            error: self.error / rhs,
        }
    }
}

impl Div<Uncertain> for Uncertain {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self {
            value: self.value / rhs.value,
            error: ((self.error / rhs.value).powi(2)
                + (self.value * rhs.error / (rhs.value).powi(2)).powi(2))
            .sqrt(),
        }
    }
}

impl Div<Uncertain> for f64 {
    type Output = Uncertain;

    fn div(self, rhs: Uncertain) -> Uncertain {
        Uncertain {
            value: self / rhs.value,
            error: self * rhs.error / (rhs.value).powi(2),
        }
    }
}

impl Uncertain {
    /// Mean and standard deviation calculation, a method to use for
    /// measurements with close distribution width, such as attenuation or cell
    /// doubling rate.
    ///
    /// p.54 in Data reduction and error analysis for physical sciences
    /// (Bevington, Robinson)
    pub fn mean_with_standard_deviation(set: &[Self]) -> Option<MeanWithStandardDeviation> {
        if set.is_empty() || set.len() == 1 {
            return None;
        }

        let n = set.len() as f64;

        let mut sum_values = 0.0f64;
        for uncertain in set.iter() {
            sum_values += uncertain.value;
        }
        let mean = sum_values / n;

        let mut sum_deviation_squares = 0.0f64;
        for uncertain in set.iter() {
            sum_deviation_squares += (mean - uncertain.value).powi(2);
        }
        let standard_deviation = (sum_deviation_squares / (n - 1.0)).sqrt();

        Some(MeanWithStandardDeviation {
            mean,
            standard_deviation,
        })
    }
}

/// Mean value with standard deviation
///
/// Used for set of values that are expected to have close error values
/// throughout the set and follow Normal distribution, such as real attenuation
/// for a particular strain.
#[derive(Debug, Copy, Clone, Serialize)]
pub struct MeanWithStandardDeviation {
    /// Mean of the sample data
    pub mean: f64,

    /// Standard deviation of the sample data
    pub standard_deviation: f64,
}
