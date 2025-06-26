//! Counts line parser for CTab files.

use nom::{
    bytes::take, character::complete::multispace0, combinator::{all_consuming, complete, map, verify}, error, sequence::{delimited, preceded, terminated}, Parser
};

use super::utils::fixed_width_int;

/// Parse counts line (39 characters wide)
/// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
///
/// *Values in the counts block*
/// ---------------------------------------------------------------
/// | Field   | Meaning                    | Values     | Notes   |
/// |---------|----------------------------|------------|---------|
/// | aaa     | number of atoms            | >0         | Generic |
/// | bbb     | number of bonds            | >=0        | Generic |
/// | lll     | number of atom lists       | 0..=30     | Generic |
/// | ccc     | chiral flag                | 0, 1       | Generic |
/// | sss     | number of stext entries    | >=0        | Generic |
/// | mmm     | number of properties lines | >=0        | Generic |
/// | vvvvvvv | version stamp              | V2000      | Generic |
/// ---------------------------------------------------------------
///
pub fn counts_input<'a>() -> impl Parser<&'a [u8], Output = Counts, Error = error::Error<&'a [u8]>>
{
    let atoms = fixed_width_int::<i32>(3);
    let bonds = fixed_width_int::<i32>(3);
    let atom_lists = fixed_width_int::<i32>(3);
    let chiral_flag = fixed_width_int::<i32>(3);
    let stext_entries = fixed_width_int::<i32>(3);
    let properties_lines = fixed_width_int::<i32>(3);
    let version = complete(verify(take(6usize), |s: &[u8]| {
        s == b" V2000" || s == b"V2000 "
    }));
    all_consuming(terminated(map(
        (
            atoms,
            bonds,
            atom_lists,
            preceded(take(3usize), chiral_flag),
            stext_entries,
            delimited(take(12usize), properties_lines, version),
        ),
        |(atoms, bonds, atom_lists, chiral_flag, stext_entries, properties_lines)| Counts {
                atoms,
                bonds,
                atom_lists,
                chiral_flag,
                stext_entries,
                properties_lines,
        },
    ), multispace0))
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Counts {
    atoms: i32,      // 'aaa' - number of atoms (max 255)
    bonds: i32,      // 'bbb' - number of bonds (max 255)
    atom_lists: i32, // 'lll' - number of atom lists (max 30)
    // fff is obsolete, skipping
    chiral_flag: i32,   // 'ccc' - chiral flag (0=not chiral, 1=chiral)
    stext_entries: i32, // 'sss' - number of stext entries
    // xxx is obsolete, skipping
    // rrr is obsolete, skipping
    // ppp is obsolete, skipping
    // iii is obsolete, skipping
    properties_lines: i32, // 'mmm' - number of additional properties lines
                           // 'vvvvvv' - version stamp (V2000), fixed string, skipping
}

impl Counts {
    pub fn atoms(&self) -> i32 {
        self.atoms
    }
    pub fn bonds(&self) -> i32 {
        self.bonds
    }
    pub fn atom_lists(&self) -> i32 {
        self.atom_lists
    }
    pub fn chiral_flag(&self) -> i32 {
        self.chiral_flag
    }
    pub fn stext_entries(&self) -> i32 {
        self.stext_entries
    }
    pub fn properties_lines(&self) -> i32 {
        self.properties_lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::{error::ErrorKind, Err};
    use rstest::rstest;

    #[rstest]
    // From CTab spec (Figure 3)
    #[case(b"  6  5  0  0  1                 3 V2000",
      Counts {atoms: 6, bonds: 5, atom_lists: 0, chiral_flag: 1, stext_entries: 0, properties_lines: 3})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0999 V2000",
      Counts {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0000 V2000    ",
      Counts {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 0})]
    #[case(b"  4  2  0  0  0  0  0  0  0  0999 V2000",
      Counts {atoms: 4, bonds: 2, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    #[case(b"  1  0  0  0  0  0            999 V2000",
      Counts {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    fn test_counts_input(#[case] input: &[u8], #[case] expected: Counts) {
        let (remaining, counts) = counts_input().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(counts, expected);
    }

    #[rstest]
    #[case(
        b"  4  2  0     0  0            999 V1000",
        "invalid version",
        ErrorKind::Verify
    )]
    #[case(b"  4  2  0     0  0            ", "too short", ErrorKind::Eof)]
    #[case(
        b" 1A  2  0     0  0            999 V2000",
        "non-numeric atom",
        ErrorKind::Eof
    )]
    #[case(
        b"  4 AA  0     0  0            999 V2000",
        "non-numeric bond",
        ErrorKind::Digit
    )]
    fn test_counts_input_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let res = counts_input().parse(input);
        assert!(res.is_err(), "{} should have failed", desc);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code),
        );
    }
}
