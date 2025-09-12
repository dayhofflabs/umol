//! Lexer for SMILES tokens

use std::fmt;
use std::num::ParseIntError;

use logos::{Logos, SpannedIter};

#[derive(Debug, Default, Clone, PartialEq)]
pub enum LexicalError {
    #[default]
    InvalidToken,
    InvalidNumber(ParseIntError),
}

impl From<ParseIntError> for LexicalError {
    fn from(error: ParseIntError) -> Self {
        LexicalError::InvalidNumber(error)
    }
}

#[derive(Logos, Debug, Clone, PartialEq)]
#[logos(error = LexicalError)]
pub enum Token {
    // Atoms
    #[token("Ac")]
    Ac,
    #[token("Ag")]
    Ag,
    #[token("Al")]
    Al,
    #[token("Am")]
    Am,
    #[token("Ar")]
    Ar,
    #[token("As")]
    As,
    #[token("At")]
    At,
    #[token("Au")]
    Au,
    #[token("B")]
    B,
    #[token("Ba")]
    Ba,
    #[token("Be")]
    Be,
    #[token("Bh")]
    Bh,
    #[token("Bi")]
    Bi,
    #[token("Bk")]
    Bk,
    #[token("Br")]
    Br,
    #[token("C")]
    C,
    #[token("Ca")]
    Ca,
    #[token("Cd")]
    Cd,
    #[token("Ce")]
    Ce,
    #[token("Cf")]
    Cf,
    #[token("Cl")]
    Cl,
    #[token("Cm")]
    Cm,
    #[token("Cn")]
    Cn,
    #[token("Co")]
    Co,
    #[token("Cr")]
    Cr,
    #[token("Cs")]
    Cs,
    #[token("Cu")]
    Cu,
    #[token("Db")]
    Db,
    #[token("Ds")]
    Ds,
    #[token("Dy")]
    Dy,
    #[token("Er")]
    Er,
    #[token("Es")]
    Es,
    #[token("Eu")]
    Eu,
    #[token("F")]
    F,
    #[token("Fe")]
    Fe,
    #[token("Fl")]
    Fl,
    #[token("Fm")]
    Fm,
    #[token("Fr")]
    Fr,
    #[token("Ga")]
    Ga,
    #[token("Gd")]
    Gd,
    #[token("Ge")]
    Ge,
    #[token("H")]
    H,
    #[token("He")]
    He,
    #[token("Hf")]
    Hf,
    #[token("Hg")]
    Hg,
    #[token("Ho")]
    Ho,
    #[token("Hs")]
    Hs,
    #[token("I")]
    I,
    #[token("In")]
    In,
    #[token("Ir")]
    Ir,
    #[token("K")]
    K,
    #[token("Kr")]
    Kr,
    #[token("La")]
    La,
    #[token("Li")]
    Li,
    #[token("Lr")]
    Lr,
    #[token("Lu")]
    Lu,
    #[token("Lv")]
    Lv,
    #[token("Mc")]
    Mc,
    #[token("Md")]
    Md,
    #[token("Mg")]
    Mg,
    #[token("Mn")]
    Mn,
    #[token("Mo")]
    Mo,
    #[token("Mt")]
    Mt,
    #[token("N")]
    N,
    #[token("Na")]
    Na,
    #[token("Nb")]
    Nb,
    #[token("Nd")]
    Nd,
    #[token("Ne")]
    Ne,
    #[token("Nh")]
    Nh,
    #[token("Ni")]
    Ni,
    #[token("No")]
    No,
    #[token("Np")]
    Np,
    #[token("O")]
    O,
    #[token("Og")]
    Og,
    #[token("Os")]
    Os,
    #[token("P")]
    P,
    #[token("Pa")]
    Pa,
    #[token("Pb")]
    Pb,
    #[token("Pd")]
    Pd,
    #[token("Pm")]
    Pm,
    #[token("Po")]
    Po,
    #[token("Pr")]
    Pr,
    #[token("Pt")]
    Pt,
    #[token("Pu")]
    Pu,
    #[token("Ra")]
    Ra,
    #[token("Rb")]
    Rb,
    #[token("Re")]
    Re,
    #[token("Rf")]
    Rf,
    #[token("Rg")]
    Rg,
    #[token("Rh")]
    Rh,
    #[token("Rn")]
    Rn,
    #[token("Ru")]
    Ru,
    #[token("S")]
    S,
    #[token("Sb")]
    Sb,
    #[token("Sc")]
    Sc,
    #[token("Se")]
    Se,
    #[token("Sg")]
    Sg,
    #[token("Si")]
    Si,
    #[token("Sm")]
    Sm,
    #[token("Sn")]
    Sn,
    #[token("Sr")]
    Sr,
    #[token("Ta")]
    Ta,
    #[token("Tb")]
    Tb,
    #[token("Tc")]
    Tc,
    #[token("Te")]
    Te,
    #[token("Th")]
    Th,
    #[token("Ti")]
    Ti,
    #[token("Tl")]
    Tl,
    #[token("Tm")]
    Tm,
    #[token("Ts")]
    Ts,
    #[token("U")]
    U,
    #[token("V")]
    V,
    #[token("W")]
    W,
    #[token("Xe")]
    Xe,
    #[token("Y")]
    Y,
    #[token("Yb")]
    Yb,
    #[token("Zn")]
    Zn,
    #[token("Zr")]
    Zr,

    // Aromatic atoms
    #[token("as")]
    AromAs,
    #[token("b")]
    AromB,
    #[token("c")]
    AromC,
    #[token("n")]
    AromN,
    #[token("o")]
    AromO,
    #[token("p")]
    AromP,
    #[token("s")]
    AromS,
    #[token("se")]
    AromSe,

    // Asterisk
    #[token("*")]
    Asterisk,

    // Brackets and parentheses
    #[token("[")]
    OpenBracket,
    #[token("]")]
    CloseBracket,
    #[token("(")]
    OpenParen,
    #[token(")")]
    CloseParen,

    // Numbers
    // Single decimal digit 0-9
    #[regex(r"[0-9]", |lex| lex.slice().parse::<u32>())]
    Digit(u32),
    // Percent-prefixed two-digit ring index (%10..%99)
    #[regex(r"%[1-9][0-9]", |lex| lex.slice()[1..].parse::<u32>())]
    Percent(u32),

    // Chirality flags
    #[token("@")]
    Clockwise,
    #[token("@@")]
    CounterClockwise,
    #[token("@TH")]
    Tetrahedral,
    #[token("@AL")]
    Allenal,
    #[token("@SP")]
    SquarePlanar,
    #[token("@TB")]
    TrigonalBipyramidal,
    #[token("@OH")]
    Octahedral,

    // Non-alphanumeric characters
    #[token("+")]
    Plus,
    #[token("-")]
    Minus, // Same as Dash
    #[token("++")]
    PlusTwo, // Deprecated
    #[token("--")]
    MinusTwo, // Deprecated
    #[token("=")]
    Equal,
    #[token("#")]
    Hash,
    #[token("$")]
    Dollar,
    #[token(":")]
    Colon,
    #[token("/")]
    Slash,
    #[token("\\")]
    Backslash,
    #[token(".")]
    Dot,

    // Terminator
    #[regex("[ \t\n\r]")]
    Stop,

    // Error
    Error,
}

impl fmt::Display for Token {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

// Lexer iterator type
pub type Spanned<Tok, Loc, Error> = Result<(Loc, Tok, Loc), Error>;

pub struct Lexer<'input> {
    token_stream: SpannedIter<'input, Token>,
}

impl<'input> Lexer<'input> {
    pub fn new(input: &'input str) -> Self {
        Self {
            token_stream: Token::lexer(input).spanned(),
        }
    }
}

impl<'input> Iterator for Lexer<'input> {
    type Item = Spanned<Token, usize, LexicalError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.token_stream.next().map(|(tok, span)| match tok {
            Ok(token) => Ok((span.start, token, span.end)),
            Err(_) => Ok((span.start, Token::Error, span.end)),
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case("C", vec![Token::C])]
    #[case("c", vec![Token::AromC])]
    #[case("CC", vec![Token::C, Token::C])]
    #[case("C.C", vec![Token::C, Token::Dot, Token::C])]
    #[case("C/C", vec![Token::C, Token::Slash, Token::C])]
    #[case("C\\C", vec![Token::C, Token::Backslash, Token::C])]
    #[case("C:C", vec![Token::C, Token::Colon, Token::C])]
    #[case("C=C", vec![Token::C, Token::Equal, Token::C])]
    #[case("C#C", vec![Token::C, Token::Hash, Token::C])]
    fn test_token(#[case] input: &str, #[case] expected: Vec<Token>) {
        let tokens = Token::lexer(input)
            .map(|t| t.ok())
            .collect::<Option<Vec<Token>>>();
        assert!(tokens.is_some(), "{} should have succeeded", input);
        let tokens = tokens.unwrap();
        assert_eq!(tokens, expected);
    }

    #[rstest]
    #[case("X", 1)]
    #[case("f", 1)]
    #[case(">>", 2)]
    #[case(",", 1)]
    fn test_token_invalid(#[case] input: &str, #[case] expected_count: usize) {
        let tokens = Token::lexer(input);
        let errors = tokens.map(|t| t.unwrap_err()).collect::<Vec<_>>();
        assert_eq!(errors.len(), expected_count);
    }

    #[rstest]
    #[case("C", vec![(0, Token::C, 1)])]
    #[case("c", vec![(0, Token::AromC, 1)])]
    #[case("CC", vec![(0, Token::C, 1), (1, Token::C, 2)])]
    #[case("C.C", vec![(0, Token::C, 1), (1, Token::Dot, 2), (2, Token::C, 3)])]
    #[case("C/C", vec![(0, Token::C, 1), (1, Token::Slash, 2), (2, Token::C, 3)])]
    #[case("C\\C", vec![(0, Token::C, 1), (1, Token::Backslash, 2), (2, Token::C, 3)])]
    #[case("C:C", vec![(0, Token::C, 1), (1, Token::Colon, 2), (2, Token::C, 3)])]
    #[case("C=C", vec![(0, Token::C, 1), (1, Token::Equal, 2), (2, Token::C, 3)])]
    #[case("C#C", vec![(0, Token::C, 1), (1, Token::Hash, 2), (2, Token::C, 3)])]
    fn test_lexer(#[case] input: &str, #[case] expected: Vec<(usize, Token, usize)>) {
        let lexer = Lexer::new(input);
        let tokens = lexer.map(|t| t.ok()).collect::<Option<Vec<_>>>();
        assert!(tokens.is_some(), "{} should have succeeded", input);
        let tokens = tokens.unwrap();
        assert_eq!(tokens, expected);
    }

    #[rstest]
    #[case("X", vec![(0, Token::Error, 1)])]
    #[case("f", vec![(0, Token::Error, 1)])]
    #[case(">>", vec![(0, Token::Error, 1), (1, Token::Error, 2)])]
    #[case(",", vec![(0, Token::Error, 1)])]
    fn test_lexer_invalid(#[case] input: &str, #[case] expected: Vec<(usize, Token, usize)>) {
        let lexer = Lexer::new(input);
        let errors = lexer.map(|t| t.unwrap()).collect::<Vec<_>>();
        assert_eq!(errors, expected);
    }
}
