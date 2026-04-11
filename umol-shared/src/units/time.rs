//! Time in atomic units (seconds).

use std::ops::{Add, Div, Mul, Neg, Sub};

use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum TimeUnit {
    Yoctoseconds, // 1e-24 s
    Zeptoseconds, // 1e-21 s
    Attoseconds,  // 1e-18 s
    Femtoseconds, // 1e-15 s
    Picoseconds,  // 1e-12 s
    Nanoseconds,  // 1e-9 s
    Microseconds, // 1e-6 s
    Milliseconds, // 1e-3 s
    Seconds,
    Minutes,
    Hours,
    Days,
    Years,
    KiloYears,  // 1e3 y
    MegaYears,  // 1e6 y
    GigaYears,  // 1e9 y
    TeraYears,  // 1e12 y
    PetaYears,  // 1e15 y
    ExaYears,   // 1e18 y
    ZettaYears, // 1e21 y
    YottaYears, // 1e24 y
    ElectronVolts,
    KiloElectronVolts,
    MegaElectronVolts,
}

impl TimeUnit {
    pub const fn to_seconds_factor(&self) -> f64 {
        const S_PER_YEAR: f64 = 3.155_695_2e7; // seconds per mean tropical year
        const H_BAR: f64 = 6.582_119_569e-16; // eV·s

        match *self {
            TimeUnit::Yoctoseconds => 1e-24,
            TimeUnit::Zeptoseconds => 1e-21,
            TimeUnit::Attoseconds => 1e-18,
            TimeUnit::Femtoseconds => 1e-15,
            TimeUnit::Picoseconds => 1e-12,
            TimeUnit::Nanoseconds => 1e-9,
            TimeUnit::Microseconds => 1e-6,
            TimeUnit::Milliseconds => 1e-3,
            TimeUnit::Seconds => 1.0,
            TimeUnit::Minutes => 60.0,
            TimeUnit::Hours => 3600.0,
            TimeUnit::Days => 86400.0,
            TimeUnit::Years => S_PER_YEAR,
            TimeUnit::KiloYears => S_PER_YEAR * 1e3,
            TimeUnit::MegaYears => S_PER_YEAR * 1e6,
            TimeUnit::GigaYears => S_PER_YEAR * 1e9,
            TimeUnit::TeraYears => S_PER_YEAR * 1e12,
            TimeUnit::PetaYears => S_PER_YEAR * 1e15,
            TimeUnit::ExaYears => S_PER_YEAR * 1e18,
            TimeUnit::ZettaYears => S_PER_YEAR * 1e21,
            TimeUnit::YottaYears => S_PER_YEAR * 1e24,
            TimeUnit::ElectronVolts => H_BAR, // τ = ħ/Γ
            TimeUnit::KiloElectronVolts => H_BAR / 1e3,
            TimeUnit::MegaElectronVolts => H_BAR / 1e6,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Time(f64);

impl Time {
    pub const fn new(value: f64, unit: TimeUnit) -> Self {
        Self(value * unit.to_seconds_factor())
    }

    pub const fn seconds(v: f64) -> Self {
        Self(v)
    }

    pub const fn as_seconds(self) -> f64 {
        self.0
    }

    pub const fn as_unit(self, unit: TimeUnit) -> f64 {
        self.0 / unit.to_seconds_factor()
    }
}

impl Add for Time {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Time {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Time {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Mul<f64> for Time {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<Time> for f64 {
    type Output = Time;
    fn mul(self, rhs: Time) -> Time {
        Time(self * rhs.0)
    }
}

impl Div<f64> for Time {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl Div<Time> for Time {
    type Output = f64;
    fn div(self, rhs: Time) -> f64 {
        self.0 / rhs.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::seconds_identity(Time::seconds(1.0), 1.0)]
    #[case::years_to_seconds(Time::new(1.0, TimeUnit::Years), 3.155_695_2e7)]
    #[case::ms_to_seconds(Time::new(500.0, TimeUnit::Milliseconds), 0.5)]
    fn test_time_as_seconds(#[case] time: Time, #[case] expected: f64) {
        assert!((time.as_seconds() - expected).abs() < 1e-10);
    }

    #[rstest]
    #[case::years(12.32, TimeUnit::Years)]
    #[case::milliseconds(806.92, TimeUnit::Milliseconds)]
    fn test_time_unit_roundtrip(#[case] value: f64, #[case] unit: TimeUnit) {
        let time = Time::new(value, unit);
        assert!((time.as_unit(unit) - value).abs() < 1e-10);
    }

    #[rstest]
    fn test_time_arithmetic() {
        let a = Time::seconds(3.0);
        let b = Time::seconds(1.0);
        assert_eq!((a + b).as_seconds(), 4.0);
        assert_eq!((a - b).as_seconds(), 2.0);
        assert_eq!((a * 2.0).as_seconds(), 6.0);
        assert_eq!((2.0 * a).as_seconds(), 6.0);
        assert_eq!((a / 3.0).as_seconds(), 1.0);
        assert_eq!(a / b, 3.0);
        assert_eq!((-a).as_seconds(), -3.0);
    }
}
