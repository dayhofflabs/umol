//! Stoichiometry model
//!
//! Simple model containing only atom counts

use map_macro::hash_set;
use nom::character::complete::one_of;
use nom::character::complete::u32 as nom_u32;
use nom::combinator::{map_res, opt, recognize};
use nom::multi::many0;
use nom::sequence::pair;
use nom::{error, Parser};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Display, Formatter};
use std::ops::{Add, Mul};
use std::str::FromStr;

use umol::error::SerializationError;
use umol::{property, Capability, Error, Model, Property, Result};
use umol_data::Element;

/// A stoichiometry model representing atom counts for elements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Stoichiometry {
    counts: HashMap<Element, u32>,
}

impl Stoichiometry {
    /// Create a new empty stoichiometry
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
        }
    }

    pub fn from_counts(counts: HashMap<Element, u32>) -> Self {
        Self { counts }
    }

    /// Get the count of a specific element
    fn get_count(&self, element: Element) -> u32 {
        *self.counts.get(&element).unwrap_or(&0)
    }
}

impl Model for Stoichiometry {
    type Data = Self;

    fn capabilities(&self) -> HashSet<Capability> {
        hash_set! {
            Capability::new("basic", "stoichiometry", 1)
        }
    }

    fn data(&self) -> &Self::Data {
        self
    }
}

impl Display for Stoichiometry {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut elements: Vec<_> = self.counts.iter().collect();

        // Sort elements: C first, H second, then alphabetically
        elements.sort_by(|(a, _), (b, _)| match (*a, *b) {
            (Element::C, _) => Ordering::Less,
            (_, Element::C) => Ordering::Greater,
            (Element::H, _) => Ordering::Less,
            (_, Element::H) => Ordering::Greater,
            _ => a.symbol().cmp(b.symbol()),
        });

        for (element, &count) in elements {
            if count == 1 {
                write!(f, "{}", element.symbol())?;
            } else {
                write!(f, "{}{}", element.symbol(), count)?;
            }
        }
        Ok(())
    }
}

/// Parse element symbol (uppercase + optional lowercase letters)
fn element_symbol<'a>() -> impl Parser<&'a str, Output = Element, Error = error::Error<&'a str>> {
    map_res(
        recognize((
            one_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
            opt(one_of("abcdefghijklmnopqrstuvwxyz")),
        )),
        |s: &str| Element::from_symbol(s).ok_or("Invalid element symbol"),
    )
}

fn stoichiometry<'a>() -> impl Parser<&'a str, Output = Stoichiometry, Error = error::Error<&'a str>>
{
    many0(pair(element_symbol(), opt(nom_u32))).map(|pairs| {
        let mut counts = HashMap::new();
        for (element, count) in pairs {
            *counts.entry(element).or_insert(0) += count.unwrap_or(1);
        }
        Stoichiometry::from_counts(counts)
    })
}

/// Parse full stoichiometry formula
fn parse_stoichiometry(input: &str) -> Result<Stoichiometry> {
    let (rest, stoichiometry) = stoichiometry()
        .parse(input)
        .map_err(|e| Error::Serialization(SerializationError::ParseError(e.to_string())))?;
    if !rest.is_empty() {
        return Err(Error::Serialization(SerializationError::ParseError(
            format!("Unexpected characters at end of formula: '{}'", rest),
        )));
    }
    Ok(stoichiometry)
}

impl FromStr for Stoichiometry {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        parse_stoichiometry(s)
    }
}

impl Add for Stoichiometry {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut counts = self.counts;
        for (element, count) in rhs.counts {
            *counts.entry(element).or_insert(0) += count;
        }
        Self { counts }
    }
}

impl Mul<u32> for Stoichiometry {
    type Output = Self;

    fn mul(self, rhs: u32) -> Self {
        Self {
            counts: self
                .counts
                .into_iter()
                .map(|(element, count)| (element, count * rhs))
                .collect(),
        }
    }
}

impl PartialEq for Stoichiometry {
    fn eq(&self, other: &Self) -> bool {
        self.counts == other.counts
    }
}

impl PartialOrd for Stoichiometry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mut all_elements: Vec<_> = self.counts.keys().chain(other.counts.keys()).collect();
        all_elements.sort();
        all_elements.dedup();

        let mut ordering = Ordering::Equal;
        for &element in all_elements {
            let self_count = self.get_count(element);
            let other_count = other.get_count(element);
            match self_count.cmp(&other_count) {
                Ordering::Less => {
                    if ordering == Ordering::Greater {
                        return None;
                    }
                    ordering = Ordering::Less;
                }
                Ordering::Greater => {
                    if ordering == Ordering::Less {
                        return None;
                    }
                    ordering = Ordering::Greater;
                }
                Ordering::Equal => {}
            }
        }
        Some(ordering)
    }
}

// Atom count property

struct AtomCount;

#[property]
impl Property<Stoichiometry> for AtomCount {
    type Value = u32;
    type Args = Element;

    fn name(&self) -> String {
        "atom_count".to_string()
    }

    fn compute(&self, model: &Stoichiometry, element: Self::Args) -> Result<Self::Value> {
        Ok(model.get_count(element))
    }
}

// Mass property

struct Mass;

#[property]
impl Property<Stoichiometry> for Mass {
    type Value = f64;
    type Args = ();

    fn name(&self) -> String {
        "mass".to_string()
    }

    fn compute(&self, model: &Stoichiometry, _args: Self::Args) -> Result<Self::Value> {
        Ok(model
            .counts
            .iter()
            .map(|(&element, &count)| element.mass() * (count as f64))
            .sum())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use map_macro::hash_map;
    use umol_data::e;

    #[test]
    fn test_stoichiometry_creation() {
        let s = Stoichiometry::new();
        assert_eq!(s.get_count(e!(C)), 0);
    }

    #[test]
    fn test_stoichiometry_from_counts() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(H) => 6,
            e!(O) => 1,
        });
        assert_eq!(s.get_count(e!(C)), 2);
        assert_eq!(s.get_count(e!(H)), 6);
        assert_eq!(s.get_count(e!(O)), 1);
    }

    #[test]
    fn test_stoichiometry_display() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(H) => 6,
            e!(O) => 1,
        });
        assert_eq!(format!("{}", s), "C2H6O");
    }

    #[test]
    fn test_stoichiometry_display_single_atom() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 1,
            e!(H) => 1,
        });

        assert_eq!(format!("{}", s), "CH");
    }

    #[test]
    fn test_stoichiometry_display_empty() {
        let s = Stoichiometry::new();
        assert_eq!(format!("{}", s), "");
    }

    #[test]
    fn test_stoichiometry_parsing() {
        // Test basic parsing
        let s: Stoichiometry = "C2H6O".parse().unwrap();
        assert_eq!(s.get_count(e!(C)), 2);
        assert_eq!(s.get_count(e!(H)), 6);
        assert_eq!(s.get_count(e!(O)), 1);

        // Test single atom without count
        let s: Stoichiometry = "CH".parse().unwrap();
        assert_eq!(s.get_count(e!(C)), 1);
        assert_eq!(s.get_count(e!(H)), 1);

        // Test multi-letter elements
        let s: Stoichiometry = "NaCl".parse().unwrap();
        assert_eq!(s.get_count(e!(Na)), 1);
        assert_eq!(s.get_count(e!(Cl)), 1);

        // Test empty string
        let s: Stoichiometry = "".parse().unwrap();
        assert_eq!(s.get_count(e!(C)), 0);
        assert_eq!(s.get_count(e!(H)), 0);
    }

    #[test]
    fn test_stoichiometry_parsing_errors() {
        // Test invalid element symbol
        assert!(matches!(
            "X".parse::<Stoichiometry>(),
            Err(Error::Serialization(SerializationError::ParseError(_)))
        ));

        // Test invalid count
        assert!(matches!(
            "C2H6O1.5".parse::<Stoichiometry>(),
            Err(Error::Serialization(SerializationError::ParseError(_)))
        ));

        // Test invalid character
        assert!(matches!(
            "C2H6O!".parse::<Stoichiometry>(),
            Err(Error::Serialization(SerializationError::ParseError(_)))
        ));

        // Test lowercase start
        assert!(matches!(
            "c2H6O".parse::<Stoichiometry>(),
            Err(Error::Serialization(SerializationError::ParseError(_)))
        ));
    }

    #[test]
    fn test_stoichiometry_serialization() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(H) => 6,
        });

        // Test serialization and deserialization
        let json = serde_json::to_string(&s).unwrap();
        let deserialized: Stoichiometry = serde_json::from_str(&json).unwrap();

        // Verify the counts match
        assert_eq!(deserialized.get_count(e!(C)), 2);
        assert_eq!(deserialized.get_count(e!(H)), 6);
        assert_eq!(deserialized.get_count(e!(O)), 0);
    }

    #[test]
    fn test_stoichiometry_addition() {
        let s1 = Stoichiometry::from_counts(hash_map! {
            e!(C) => 1,
            e!(H) => 4,
        });

        let s2 = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(O) => 1,
        });

        let sum = s1 + s2;
        assert_eq!(sum.get_count(e!(C)), 3);
        assert_eq!(sum.get_count(e!(H)), 4);
        assert_eq!(sum.get_count(e!(O)), 1);
    }

    #[test]
    fn test_stoichiometry_multiplication() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 1,
            e!(H) => 4,
        });

        let doubled = s * 2;
        assert_eq!(doubled.get_count(e!(C)), 2);
        assert_eq!(doubled.get_count(e!(H)), 8);
    }

    #[test]
    fn test_property_access() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(H) => 6,
        });

        // Test atom count property
        assert_eq!(s.atom_count(e!(C)).unwrap(), 2);
        assert_eq!(s.atom_count(e!(H)).unwrap(), 6);
        assert_eq!(s.atom_count(e!(O)).unwrap(), 0);

        // Test mass property
        let mass = s.mass().unwrap();
        let expected_mass = 2.0 * e!(C).mass() + 6.0 * e!(H).mass();
        assert!((mass - expected_mass).abs() < 1e-10);
    }

    #[test]
    fn test_property_capabilities() {
        let s = Stoichiometry::new();
        let caps = s.capabilities();

        assert!(caps.contains(&Capability::new("basic", "stoichiometry", 1)));
        assert_eq!(caps.len(), 1);
    }

    #[test]
    fn test_property_computation() {
        let s = Stoichiometry::from_counts(hash_map! {
            e!(C) => 2,
            e!(H) => 6,
        });

        assert_eq!(s.atom_count(e!(C)).unwrap(), 2);
        assert_eq!(s.atom_count(e!(H)).unwrap(), 6);
        assert_eq!(s.atom_count(e!(O)).unwrap(), 0);

        let mass = s.mass().unwrap();
        let expected_mass = 2.0 * e!(C).mass() + 6.0 * e!(H).mass();
        assert!((mass - expected_mass).abs() < 1e-10);
    }
}
