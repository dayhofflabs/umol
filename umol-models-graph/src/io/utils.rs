//! IO utility traits and functions for MOL file parsing.

use fixed_width::LineBreak;
use num;
use std::error;
use std::io::{self, BufRead};
use std::iter::Iterator;

/// Enum representing the different line break types in a MOL file.
/// The V2000 format definition calls for CRLF line breaks, but most
/// implementations support all possible line break types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LineTerminator {
    CRLF,
    LF,
    CR,
    Unknown,
}

impl TryFrom<LineTerminator> for LineBreak {
    type Error = io::Error;

    fn try_from(line_break: LineTerminator) -> Result<Self, Self::Error> {
        match line_break {
            LineTerminator::CRLF => Ok(LineBreak::CRLF),
            LineTerminator::LF => Ok(LineBreak::Newline),
            _ => Err(io::Error::new(
                io::ErrorKind::Other,
                "Unsupported line break",
            )),
        }
    }
}

/// Detects the line break type in a MOL file.
pub(crate) fn detect_line_break(reader: &mut impl BufRead) -> Result<LineTerminator, io::Error> {
    let buffer = match reader.fill_buf() {
        Ok(buffer) if !buffer.is_empty() => buffer,
        _ => return Err(io::Error::new(io::ErrorKind::Other, "Empty file")),
    };
    for i in 0..buffer.len() - 1 {
        if buffer[i] == b'\r' && buffer[i + 1] == b'\n' {
            return Ok(LineTerminator::CRLF);
        } else if buffer[i] == b'\n' {
            return Ok(LineTerminator::LF);
        } else if buffer[i] == b'\r' {
            return Ok(LineTerminator::CR);
        }
    }
    Ok(LineTerminator::Unknown)
}

/// Auxiliary trait for items capable of being combined with other items.
/// Used to implement `combine_next_n` for iterators.
pub(crate) trait CombineableItem {
    fn combine(&mut self, other: &Self, separator: &str);
}

impl CombineableItem for String {
    fn combine(&mut self, other: &Self, separator: &str) {
        self.push_str(separator);
        self.push_str(other);
    }
}

impl CombineableItem for Vec<u8> {
    fn combine(&mut self, other: &Self, separator: &str) {
        self.extend_from_slice(separator.as_bytes());
        self.extend_from_slice(other);
    }
}

impl<C: CombineableItem> CombineableItem for Option<C> {
    fn combine(&mut self, other: &Self, separator: &str) {
        match (self, other) {
            (Some(s1), Some(s2)) => {
                s1.combine(s2, separator);
            }
            _ => {}
        }
    }
}

impl<C: CombineableItem, E: error::Error> CombineableItem for Result<C, E> {
    fn combine(&mut self, other: &Self, separator: &str) {
        match (self, other) {
            (Ok(s1), Ok(s2)) => {
                s1.combine(s2, separator);
            }
            _ => {}
        }
    }
}

/// Iterator adapter that combines the current item with the next n items
///  if the predicate returns `Some`.
pub(crate) trait CombineNextN: Iterator {
    fn combine_next_n<F, N>(self, predicate: F, separator: &str) -> CombineNextNIter<Self, F, N>
    where
        Self: Sized,
        F: FnMut(&Self::Item) -> Option<N>,
        N: num::Unsigned + Into<usize>,
        Self::Item: CombineableItem;
}

impl<I> CombineNextN for I
where
    I: Iterator,
    I::Item: CombineableItem,
{
    fn combine_next_n<F, N>(self, predicate: F, separator: &str) -> CombineNextNIter<Self, F, N>
    where
        F: FnMut(&Self::Item) -> Option<N>,
        N: num::Unsigned + Into<usize>,
    {
        CombineNextNIter {
            inner: self,
            predicate,
            separator: separator.to_string(),
            _phantom: std::marker::PhantomData,
        }
    }
}

pub(crate) struct CombineNextNIter<I, F, N>
where
    I: Iterator,
    I::Item: CombineableItem,
    F: FnMut(&I::Item) -> Option<N>,
    N: num::Unsigned + Into<usize>,
{
    inner: I,
    predicate: F,
    separator: String,
    _phantom: std::marker::PhantomData<N>,
}

impl<I, F, N> Iterator for CombineNextNIter<I, F, N>
where
    I: Iterator,
    I::Item: CombineableItem,
    F: FnMut(&I::Item) -> Option<N>,
    N: num::Unsigned + Into<usize>,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|item| {
            if let Some(count) = (self.predicate)(&item) {
                let mut result = item;

                for _ in 0..count.into() {
                    if let Some(next) = self.inner.next() {
                        result.combine(&next, &self.separator);
                    } else {
                        break;
                    }
                }

                result
            } else {
                item
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_line_break() {
        let mut reader = io::Cursor::new(b"Line 1\r\nLine 2\nLine 3\rLine 4\n");
        assert_eq!(
            detect_line_break(&mut reader).unwrap(),
            LineTerminator::CRLF
        );
        let mut reader = io::Cursor::new(b"Line 1\nLine 2\nLine 3\n");
        assert_eq!(detect_line_break(&mut reader).unwrap(), LineTerminator::LF);
        let mut reader = io::Cursor::new(b"Line 1\rLine 2\rLine 3\r");
        assert_eq!(detect_line_break(&mut reader).unwrap(), LineTerminator::CR);
        let mut reader = io::Cursor::new(b"Line 1\nLine 2\nLine 3\n");
        assert_eq!(detect_line_break(&mut reader).unwrap(), LineTerminator::LF);
        let mut reader = io::Cursor::new(b"");
        assert!(detect_line_break(&mut reader).is_err());
    }

    #[test]
    fn test_combine_next_n_string() {
        let lines = vec![
            "Line 1".to_string(),
            "Line 2 (join 2)".to_string(),
            "Line 3".to_string(),
            "Line 4".to_string(),
            "Line 5".to_string(),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if line.contains("join") {
                        Some(1u16)
                    } else {
                        None
                    }
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert_eq!(
            result,
            vec!["Line 1", "Line 2 (join 2)\nLine 3", "Line 4", "Line 5"]
        );
    }

    #[test]
    fn test_combine_next_n_option_string() {
        let lines = vec![
            Some("Line 1".to_string()),
            Some("Line 2 (join 2)".to_string()),
            Some("Line 3".to_string()),
            None,
            Some("Line 5".to_string()),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if let Some(line) = line.as_ref() {
                        if line.contains("join") {
                            return Some(1u16);
                        }
                    }
                    None
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert_eq!(
            result,
            vec![
                Some("Line 1".to_string()),
                Some("Line 2 (join 2)\nLine 3".to_string()),
                None,
                Some("Line 5".to_string()),
            ]
        );
    }

    #[test]
    fn test_combine_next_n_result_string() {
        let lines = vec![
            Ok("Line 1".to_string()),
            Ok("Line 2 (join 2)".to_string()),
            Ok("Line 3".to_string()),
            Err(std::io::Error::new(std::io::ErrorKind::Other, "Error")),
            Ok("Line 5".to_string()),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if let Ok(line) = line.as_ref() {
                        if line.contains("join") {
                            return Some(1u16);
                        }
                    }
                    None
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert!(matches!(result[0], Ok(_)) && result[0].as_ref().unwrap() == "Line 1");
        assert!(
            matches!(result[1], Ok(_)) && result[1].as_ref().unwrap() == "Line 2 (join 2)\nLine 3"
        );
        assert!(matches!(result[2], Err(_)));
        assert!(matches!(result[3], Ok(_)) && result[3].as_ref().unwrap() == "Line 5");
    }

    #[test]
    fn test_combine_next_n_vec_u8() {
        let lines = vec![
            b"Line 1".to_vec(),
            b"Line 2 \\".to_vec(),
            b"Line 3".to_vec(),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if line.contains(&b'\\') {
                        Some(1u16)
                    } else {
                        None
                    }
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert_eq!(
            result,
            vec![b"Line 1".to_vec(), b"Line 2 \\\nLine 3".to_vec()]
        );
    }

    #[test]
    fn test_combine_next_n_option_vec_u8() {
        let lines = vec![
            Some(b"Line 1".to_vec()),
            Some(b"Line 2 \\".to_vec()),
            Some(b"Line 3".to_vec()),
            None,
            Some(b"Line 5".to_vec()),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if line.is_some() && line.as_ref().unwrap().contains(&b'\\') {
                        Some(1u16)
                    } else {
                        None
                    }
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert_eq!(
            result,
            vec![
                Some(b"Line 1".to_vec()),
                Some(b"Line 2 \\\nLine 3".to_vec()),
                None,
                Some(b"Line 5".to_vec()),
            ]
        );
    }

    #[test]
    fn test_combine_next_n_result_vec_u8() {
        let lines = vec![
            Ok(b"Line 1".to_vec()),
            Ok(b"Line 2 \\".to_vec()),
            Ok(b"Line 3".to_vec()),
            Err(std::io::Error::new(std::io::ErrorKind::Other, "Error")),
            Ok(b"Line 5".to_vec()),
        ];
        let result = lines
            .into_iter()
            .combine_next_n(
                |line| {
                    if let Ok(line) = line.as_ref() {
                        if line.contains(&b'\\') {
                            return Some(1u16);
                        }
                    }
                    None
                },
                "\n",
            )
            .collect::<Vec<_>>();

        assert!(matches!(result[0], Ok(_)) && result[0].as_ref().unwrap() == b"Line 1");
        assert!(matches!(result[1], Ok(_)) && result[1].as_ref().unwrap() == b"Line 2 \\\nLine 3");
        assert!(matches!(result[2], Err(_)));
        assert!(matches!(result[3], Ok(_)) && result[3].as_ref().unwrap() == b"Line 5");
    }
}
