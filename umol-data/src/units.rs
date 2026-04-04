//! Physical units for umol.
//!
//! Internal unit system: atomic units (Bohr, radians).
//! Named constructors enforce unit conversion at API boundaries.

use std::ops::{Add, Div, Mul, Neg, Sub};

// CODATA 2018: 1 Bohr = 0.529177210903 Angstrom
const BOHR_PER_ANGSTROM: f64 = 1.8897259886;
const ANGSTROM_PER_BOHR: f64 = 1.0 / BOHR_PER_ANGSTROM;

/// Length in atomic units (Bohr).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Length(f64);

impl Length {
    pub const fn bohr(v: f64) -> Self {
        Self(v)
    }

    pub const fn angstrom(v: f64) -> Self {
        Self(v * BOHR_PER_ANGSTROM)
    }

    pub const fn picometer(v: f64) -> Self {
        Self(v * 0.01 * BOHR_PER_ANGSTROM)
    }

    pub const fn as_bohr(self) -> f64 {
        self.0
    }

    pub const fn as_angstrom(self) -> f64 {
        self.0 * ANGSTROM_PER_BOHR
    }
}

impl Add for Length {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Length {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Length {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Mul<f64> for Length {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<Length> for f64 {
    type Output = Length;
    fn mul(self, rhs: Length) -> Length {
        Length(self * rhs.0)
    }
}

impl Div<f64> for Length {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl Div<Length> for Length {
    type Output = f64;
    fn div(self, rhs: Length) -> f64 {
        self.0 / rhs.0
    }
}

const RADIANS_PER_DEGREE: f64 = std::f64::consts::PI / 180.0;
const DEGREES_PER_RADIAN: f64 = 180.0 / std::f64::consts::PI;

/// Angle in radians.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Angle(f64);

impl Angle {
    pub fn radians(v: f64) -> Self {
        Self(v)
    }

    pub fn degrees(v: f64) -> Self {
        Self(v * RADIANS_PER_DEGREE)
    }

    pub fn as_radians(self) -> f64 {
        self.0
    }

    pub fn as_degrees(self) -> f64 {
        self.0 * DEGREES_PER_RADIAN
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
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::bohr_identity(Length::bohr(1.0), 1.0)]
    #[case::angstrom_to_bohr(Length::angstrom(1.0), BOHR_PER_ANGSTROM)]
    #[case::picometer_to_bohr(Length::picometer(100.0), Length::angstrom(1.0).as_bohr())]
    fn test_length_as_bohr(#[case] length: Length, #[case] expected: f64) {
        assert!((length.as_bohr() - expected).abs() < 1e-10);
    }

    #[rstest]
    #[case::angstrom_roundtrip(1.54)]
    #[case::zero(0.0)]
    fn test_length_angstrom_roundtrip(#[case] v: f64) {
        assert!((Length::angstrom(v).as_angstrom() - v).abs() < 1e-14);
    }

    #[rstest]
    fn test_length_arithmetic() {
        let a = Length::bohr(3.0);
        let b = Length::bohr(1.0);
        assert_eq!((a + b).as_bohr(), 4.0);
        assert_eq!((a - b).as_bohr(), 2.0);
        assert_eq!((a * 2.0).as_bohr(), 6.0);
        assert_eq!((2.0 * a).as_bohr(), 6.0);
        assert_eq!((a / 3.0).as_bohr(), 1.0);
        assert_eq!(a / b, 3.0);
        assert_eq!((-a).as_bohr(), -3.0);
    }

    #[rstest]
    #[case::radians_identity(Angle::radians(1.0), 1.0)]
    #[case::degrees_to_radians(Angle::degrees(180.0), std::f64::consts::PI)]
    fn test_angle_as_radians(#[case] angle: Angle, #[case] expected: f64) {
        assert!((angle.as_radians() - expected).abs() < 1e-14);
    }

    #[rstest]
    #[case::degrees_roundtrip(104.5)]
    #[case::zero(0.0)]
    fn test_angle_degrees_roundtrip(#[case] v: f64) {
        assert!((Angle::degrees(v).as_degrees() - v).abs() < 1e-12);
    }

    #[rstest]
    fn test_angle_trig() {
        let right = Angle::degrees(90.0);
        assert!((right.sin() - 1.0).abs() < 1e-15);
        assert!(right.cos().abs() < 1e-15);
    }
}
