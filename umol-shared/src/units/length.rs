//! Length in atomic units (Bohr).

use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LengthUnit {
    Bohr,
    Angstrom,
    Picometer,
    Nanometer,
}

impl LengthUnit {
    // CODATA 2018: 1 Bohr = 0.529177210903 Angstrom
    pub const fn to_bohr_factor(&self) -> f64 {
        const BOHR_PER_ANGSTROM: f64 = 1.8897259886;
        match *self {
            LengthUnit::Bohr => 1.0,
            LengthUnit::Angstrom => BOHR_PER_ANGSTROM,
            LengthUnit::Picometer => BOHR_PER_ANGSTROM * 0.01,
            LengthUnit::Nanometer => BOHR_PER_ANGSTROM * 10.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Length(f64);

impl Length {
    pub const fn new(value: f64, unit: LengthUnit) -> Self {
        Self(value * unit.to_bohr_factor())
    }

    pub const fn bohr(v: f64) -> Self {
        Self::new(v, LengthUnit::Bohr)
    }

    pub const fn angstrom(v: f64) -> Self {
        Self::new(v, LengthUnit::Angstrom)
    }

    pub const fn picometer(v: f64) -> Self {
        Self::new(v, LengthUnit::Picometer)
    }

    pub const fn as_bohr(self) -> f64 {
        self.0
    }

    pub const fn as_unit(self, unit: LengthUnit) -> f64 {
        self.0 / unit.to_bohr_factor()
    }

    pub const fn as_angstrom(self) -> f64 {
        self.as_unit(LengthUnit::Angstrom)
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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::bohr_identity(Length::bohr(1.0), 1.0)]
    #[case::angstrom_to_bohr(Length::angstrom(1.0), 1.8897259886)]
    #[case::picometer_to_bohr(Length::picometer(100.0), Length::angstrom(1.0).as_bohr())]
    fn test_length_as_bohr(#[case] length: Length, #[case] expected: f64) {
        assert!((length.as_bohr() - expected).abs() < 1e-10);
    }

    #[rstest]
    #[case::angstrom(1.54, LengthUnit::Angstrom)]
    #[case::picometer(154.0, LengthUnit::Picometer)]
    fn test_length_unit_roundtrip(#[case] value: f64, #[case] unit: LengthUnit) {
        let length = Length::new(value, unit);
        assert!((length.as_unit(unit) - value).abs() < 1e-10);
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
}
