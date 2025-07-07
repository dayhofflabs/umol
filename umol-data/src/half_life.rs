//! Time unnits and half-lives of isotopes for umol-data.

use serde::{Deserialize, Serialize};

/// Represents a unit of time for half-life values.
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
    pub fn to_seconds_factor(&self) -> f64 {
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

/// Represents the half-life of an isotope.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct HalfLife {
    pub value: f64,
    pub unit: TimeUnit,
}

impl HalfLife {
    /// Calculates the half-life in seconds.
    pub fn to_seconds(&self) -> f64 {
        self.value * self.unit.to_seconds_factor()
    }
}
