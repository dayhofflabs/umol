//! XYZ file format parser and formatter.
//!
//! Standard XYZ format: atom count on line 1, comment on line 2,
//! then one line per atom: `Symbol x y z` with coordinates in Angstroms.

use std::error::Error;
use std::fmt::{self, Write};
use std::str::FromStr;

use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_chem::units::length::Length;

use crate::molecule::Molecule;

/// Error type for XYZ parsing.
#[derive(Debug)]
pub enum ParseError {
    UnexpectedEof,
    InvalidAtomCount(String),
    AtomCountMismatch { expected: usize, found: usize },
    InvalidAtomLine { line: usize, message: String },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected end of input"),
            Self::InvalidAtomCount(s) => write!(f, "invalid atom count: {s}"),
            Self::AtomCountMismatch { expected, found } => {
                write!(f, "expected {expected} atoms, found {found}")
            }
            Self::InvalidAtomLine { line, message } => {
                write!(f, "line {line}: {message}")
            }
        }
    }
}

impl Error for ParseError {}

/// Parse an XYZ string into a `Molecule` with charge and multiplicity provided.
///
/// Coordinates are in Angstroms (converted to Bohr internally).
/// Charge and spin multiplicity default to 0 and singlet.
pub fn parse_xyz_with(
    input: &str,
    charge: i32,
    multiplicity: SpinMultiplicity,
) -> Result<(Molecule, String), ParseError> {
    let mut lines = input.lines();

    let count_line = lines.next().ok_or(ParseError::UnexpectedEof)?;
    let n: usize = count_line
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidAtomCount(count_line.trim().to_string()))?;

    let comment = lines.next().ok_or(ParseError::UnexpectedEof)?.to_string();

    let mut elements = Vec::with_capacity(n);
    let mut coords = Vec::with_capacity(3 * n);

    for (idx, line) in lines.enumerate() {
        if elements.len() == n {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(ParseError::InvalidAtomLine {
                line: idx + 3,
                message: format!("expected 4 fields, found {}", parts.len()),
            });
        }
        let elem = Element::from_str(parts[0]).map_err(|_| ParseError::InvalidAtomLine {
            line: idx + 3,
            message: format!("unknown element: {}", parts[0]),
        })?;
        for &coord_str in &parts[1..4] {
            let v: f64 = coord_str.parse().map_err(|_| ParseError::InvalidAtomLine {
                line: idx + 3,
                message: format!("invalid coordinate: {coord_str}"),
            })?;
            coords.push(v);
        }
        elements.push(elem);
    }

    if elements.len() != n {
        return Err(ParseError::AtomCountMismatch {
            expected: n,
            found: elements.len(),
        });
    }

    let mol = Molecule::from_cartesian_angstrom(elements, &coords, charge, multiplicity);
    Ok((mol, comment))
}

/// Parse an XYZ string into a `Molecule` with default charge and multiplicity.
pub fn parse_xyz(input: &str) -> Result<(Molecule, String), ParseError> {
    parse_xyz_with(input, 0, SpinMultiplicity::SINGLET)
}

/// Format a molecule as an XYZ string.
///
/// Coordinates are output in Angstroms.
pub fn format_xyz<W: Write>(mol: &Molecule, comment: &str, w: &mut W) -> fmt::Result {
    let n = mol.atom_count();
    writeln!(w, "{n}")?;
    writeln!(w, "{comment}")?;
    let m = mol.cartesian_coordinates();
    for i in 0..n {
        let x = Length::bohr(m[(0, i)]).as_angstrom();
        let y = Length::bohr(m[(1, i)]).as_angstrom();
        let z = Length::bohr(m[(2, i)]).as_angstrom();
        writeln!(
            w,
            "{:<2} {:>14.8} {:>14.8} {:>14.8}",
            mol.element(i).symbol(),
            x,
            y,
            z
        )?;
    }
    Ok(())
}

/// Convenience: format a molecule as an XYZ `String`.
pub fn to_xyz_string(mol: &Molecule, comment: &str) -> String {
    let mut s = String::new();
    format_xyz(mol, comment, &mut s).expect("String::write_str is infallible");
    s
}

#[cfg(test)]
mod tests {
    use float_cmp::approx_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;

    #[fixture]
    fn water_xyz() -> &'static str {
        return "\
3
water
O    0.00000000    0.00000000    0.00000000
H    0.96000000    0.00000000    0.00000000
H   -0.24000000    0.93000000    0.00000000
";
    }

    #[rstest]
    #[case::water(water_xyz(), 3, "water", vec![Element::O, Element::H, Element::H])]
    fn test_parse_xyz_with(
        #[case] input: &'static str,
        #[case] expected_atom_count: usize,
        #[case] expected_comment: &'static str,
        #[case] expected_elements: Vec<Element>,
    ) {
        let multiplicity = SpinMultiplicity::new(11).unwrap();
        let result = parse_xyz_with(input, 1, multiplicity);
        assert!(result.is_ok());
        let (mol, comment) = result.unwrap();
        assert_eq!(comment, expected_comment);
        assert_eq!(mol.atom_count(), expected_atom_count);
        assert_eq!(mol.charge(), 1);
        assert_eq!(mol.multiplicity(), multiplicity);
        let elements = (0..mol.atom_count())
            .map(|index: usize| mol.element(index))
            .collect::<Vec<_>>();
        assert_eq!(elements, expected_elements);
    }

    #[rstest]
    #[case::water(water_xyz(), 3, "water", vec![Element::O, Element::H, Element::H])]
    fn test_parse_xyz(
        #[case] input: &'static str,
        #[case] expected_atom_count: usize,
        #[case] expected_comment: &'static str,
        #[case] expected_elements: Vec<Element>,
    ) {
        let result = parse_xyz(input);
        assert!(result.is_ok());
        let (mol, comment) = result.unwrap();
        assert_eq!(comment, expected_comment);
        assert_eq!(mol.atom_count(), expected_atom_count);
        assert_eq!(mol.charge(), 0);
        assert_eq!(mol.multiplicity(), SpinMultiplicity::SINGLET);
        let elements = (0..mol.atom_count())
            .map(|index: usize| mol.element(index))
            .collect::<Vec<_>>();
        assert_eq!(elements, expected_elements);
    }

    #[rstest]
    #[case::water(water_xyz())]
    fn test_format_xyz_roundtrip(#[case] input: &'static str) {
        let (mol, comment) = parse_xyz(input).unwrap();
        let output = to_xyz_string(&mol, &comment);
        let (mol2, comment2) = parse_xyz(&output).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(comment, comment2);
        for i in 0..mol.atom_count() {
            assert_eq!(mol.element(i), mol2.element(i));
            approx_eq!(
                f64,
                mol.distance(0, i).as_bohr(),
                mol2.distance(0, i).as_bohr(),
                epsilon = 1e-6
            );
        }
    }

    #[rstest]
    #[case::missing_count("")]
    #[case::missing_comment("3")]
    #[case::bad_count("abc\ncomment\n")]
    #[case::too_few_atoms("3\ncomment\nH 0 0 0\n")]
    #[case::bad_element("1\ncomment\nXx 0 0 0\n")]
    #[case::bad_coord("1\ncomment\nH abc 0 0\n")]
    #[case::too_few_fields("1\ncomment\nH 0 0\n")]
    fn test_parse_xyz_error(#[case] input: &str) {
        assert!(parse_xyz(input).is_err());
    }
}
