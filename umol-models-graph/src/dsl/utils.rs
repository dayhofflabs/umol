use std::fmt::Debug;

use nom::character::complete::{i8 as nom_i8, u32 as nom_u32, u8 as nom_u8};
use nom::error::Error as NomError;
use nom::Parser;

pub(crate) trait IntParser: Sized + Copy + PartialEq + Eq + PartialOrd + Ord + Debug {
    fn nom_parser<'inp>() -> impl Parser<&'inp str, Output = Self, Error = NomError<&'inp str>>;
}

impl IntParser for i8 {
    fn nom_parser<'inp>() -> impl Parser<&'inp str, Output = Self, Error = NomError<&'inp str>> {
        nom_i8
    }
}

impl IntParser for u8 {
    fn nom_parser<'inp>() -> impl Parser<&'inp str, Output = Self, Error = NomError<&'inp str>> {
        nom_u8
    }
}

impl IntParser for u32 {
    fn nom_parser<'inp>() -> impl Parser<&'inp str, Output = Self, Error = NomError<&'inp str>> {
        nom_u32
    }
}
