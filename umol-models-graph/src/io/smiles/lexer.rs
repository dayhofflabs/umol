//! Lexer for SMILES tokens

use logos::Logos;

#[derive(Logos, Debug, PartialEq)]
enum Token {
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
    OpenParenthesis,
    #[token(")")]
    CloseParenthesis,

    // Numbers
    #[regex(r"[0-9]+", |lex| lex.slice().parse::<u32>().map_err(|_| ()))]
    Number(u32),

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
    #[token("%")]
    Percent,
    #[token(".")]
    Dot,

    // Terminator
    #[regex("[ \t\n\r]")]
    Terminator,
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case("C", vec![Token::C])]
    #[case("c", vec![Token::AromC])]
    #[case("CC", vec![Token::C, Token::C])]
    #[case("C.C", vec![Token::C, Token::Dot, Token::C])]
    #[case("C/C", vec![Token::C, Token::Slash, Token::C])]
    #[case("C\\C", vec![Token::C, Token::Backslash, Token::C])]
    #[case("C%C", vec![Token::C, Token::Percent, Token::C])]
    #[case("C:C", vec![Token::C, Token::Colon, Token::C])]
    #[case("C=C", vec![Token::C, Token::Equal, Token::C])]
    #[case("C#C", vec![Token::C, Token::Hash, Token::C])]
    fn test_lexer(#[case] input: &str, #[case] expected: Vec<Token>) {
        let tokens = Token::lexer(input)
            .map(|t| t.ok())
            .collect::<Option<Vec<Token>>>();
        assert!(tokens.is_some(), "{} should have succeeded", input);
        let tokens = tokens.unwrap();
        assert_eq!(tokens, expected);
    }
}
