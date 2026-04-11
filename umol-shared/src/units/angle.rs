//! Angle in radians.

use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AngleUnit {
    Radians,
    Degrees,
}

impl AngleUnit {
    pub const fn to_radians_factor(&self) -> f64 {
        const RADIANS_PER_DEGREE: f64 = std::f64::consts::PI / 180.0;
        match *self {
            AngleUnit::Radians => 1.0,
            AngleUnit::Degrees => RADIANS_PER_DEGREE,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Angle(f64);

impl Angle {
    pub const fn new(value: f64, unit: AngleUnit) -> Self {
        Self(value * unit.to_radians_factor())
    }

    pub const fn radians(v: f64) -> Self {
        Self::new(v, AngleUnit::Radians)
    }

    pub const fn degrees(v: f64) -> Self {
        Self::new(v, AngleUnit::Degrees)
    }

    pub const fn as_radians(self) -> f64 {
        self.0
    }

    pub const fn as_unit(self, unit: AngleUnit) -> f64 {
        self.0 / unit.to_radians_factor()
    }

    pub const fn as_degrees(self) -> f64 {
        self.as_unit(AngleUnit::Degrees)
    }

    pub fn sin(self) -> f64 {
        self.0.sin()
    }

    pub fn cos(self) -> f64 {
        self.0.cos()
    }

    pub fn tan(self) -> f64 {
        self.0.tan()
    }
}

impl Add for Angle {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Angle {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Angle {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Mul<f64> for Angle {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<Angle> for f64 {
    type Output = Angle;
    fn mul(self, rhs: Angle) -> Angle {
        Angle(self * rhs.0)
    }
}

impl Div<f64> for Angle {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl Div<Angle> for Angle {
    type Output = f64;
    fn div(self, rhs: Angle) -> f64 {
        self.0 / rhs.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::radians_identity(Angle::radians(1.0), 1.0)]
    #[case::degrees_to_radians(Angle::degrees(180.0), std::f64::consts::PI)]
    fn test_angle_as_radians(#[case] angle: Angle, #[case] expected: f64) {
        assert!((angle.as_radians() - expected).abs() < 1e-14);
    }

    #[rstest]
    #[case::radians(1.5, AngleUnit::Radians)]
    #[case::degrees(104.5, AngleUnit::Degrees)]
    fn test_angle_unit_roundtrip(#[case] value: f64, #[case] unit: AngleUnit) {
        let angle = Angle::new(value, unit);
        assert!((angle.as_unit(unit) - value).abs() < 1e-12);
    }

    #[rstest]
    fn test_angle_trig() {
        let right = Angle::degrees(90.0);
        assert!((right.sin() - 1.0).abs() < 1e-15);
        assert!(right.cos().abs() < 1e-15);
    }
}
