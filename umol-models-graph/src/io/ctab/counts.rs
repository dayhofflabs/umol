//! Counts line parser for CTab files.

use nom::{
    bytes::{complete::tag, take},
    character::complete::space0,
    combinator::{all_consuming, map},
    error,
    sequence::delimited,
    Parser,
};

use super::utils::fixed_width_int;

/// Parse counts line (39 characters wide)
/// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
///
/// *Values in the counts block*
/// -------------------------------------------------------------------
/// | Field   | Meaning                    | Values     | Notes       |
/// |---------|----------------------------|------------|-------------|
/// | aaa     | number of atoms            | >0         | *[Generic]* |
/// | bbb     | number of bonds            | >=0        | *[Generic]* |
/// | lll     | number of atom lists       | 0..=30     | *[Generic]* |
/// | ccc     | chiral flag                | 0, 1       | *[Generic]* |
/// | sss     | number of stext entries    | >=0        | *[Generic]* |
/// | mmm     | number of properties lines | >=0        | *[Generic]* |
/// | vvvvvvv | version stamp              | V2000      | *[Generic]* |
/// -------------------------------------------------------------------
///
pub(crate) fn counts_line<'a>(
) -> impl Parser<&'a [u8], Output = CountsLine, Error = error::Error<&'a [u8]>> {
    let atoms = fixed_width_int::<i32>(3);
    let bonds = fixed_width_int::<i32>(3);
    let atom_lists = fixed_width_int::<i32>(3);
    let chiral_flag = fixed_width_int::<i32>(3);
    let stext_entries = fixed_width_int::<i32>(3);
    let properties_lines = fixed_width_int::<i32>(3);
    let version = delimited(space0, tag("V2000"), space0);
    all_consuming(map(
        (
            atoms,
            bonds,
            atom_lists,
            take(3usize),
            chiral_flag,
            stext_entries,
            take(12usize),
            properties_lines,
            version,
        ),
        |(atoms, bonds, atom_lists, _, chiral_flag, stext_entries, _, properties_lines, _)| {
            CountsLine {
                atoms,
                bonds,
                atom_lists,
                chiral_flag,
                stext_entries,
                properties_lines,
            }
        },
    ))
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CountsLine {
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

#[cfg(test)]
mod tests {
    use super::*;
    use nom::error::ErrorKind;
    use nom::Err;
    use rstest::rstest;

    #[rstest]
    // From CTab spec (Figure 3)
    #[case(b"  6  5  0  0  1                 3 V2000",
      CountsLine {atoms: 6, bonds: 5, atom_lists: 0, chiral_flag: 1, stext_entries: 0, properties_lines: 3})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0999 V2000",
      CountsLine {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0000 V2000    ",
      CountsLine {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 0})]
    #[case(b"  4  2  0  0  0  0  0  0  0  0999 V2000",
      CountsLine {atoms: 4, bonds: 2, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    #[case(b"  1  0  0  0  0  0            999 V2000",
      CountsLine {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    fn test_counts_line(#[case] input: &[u8], #[case] expected: CountsLine) {
        let (remaining, counts) = counts_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(counts, expected);
    }

    #[rstest]
    #[case(b"  4  2  0     0  0            999 V1000", "invalid version", ErrorKind::Tag)]
    #[case(b"  4  2  0     0  0            ", "too short", ErrorKind::TakeWhileMN)]
    #[case(b" 1A  2  0     0  0            999 V2000", "non-numeric atom", ErrorKind::TakeWhileMN)]
    #[case(b"  4 AA  0     0  0            999 V2000", "non-numeric bond", ErrorKind::TakeWhileMN)]
    fn test_counts_line_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: ErrorKind) {
        let res = counts_line().parse(input);
        assert!(res.is_err(), "{}", desc);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code),
        );
    }
}
