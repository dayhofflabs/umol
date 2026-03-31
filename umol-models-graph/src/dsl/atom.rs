//! Atom-string DSL parser

use nom::character::complete::multispace0;
use nom::combinator::all_consuming;
use nom::multi::many0;
use nom::sequence::{delimited, pair, terminated};
use nom::{Err, IResult, Parser};
use umol_data::Element;

use super::error::ParseError;
use super::lowering::LowerAst;
use super::predicates::{
    atom_predicate, element_expr, AromaticExpr, AtomPredicate, ElementExpr, HydrogenExpr,
    IsotopeExpr,
};
use super::value::ValueAst;

/// Parsed atom-string AST
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomAst {
    pub element: ElementExpr,
    pub isotope_mass: Option<IsotopeExpr>,
    pub charge: Option<ValueAst>,
    pub implicit_hydrogens: Option<HydrogenExpr>,
    pub lone_pairs: Option<ValueAst>,
    pub unpaired_electrons: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
    pub valence: Option<ValueAst>,
    pub donated_pairs: Option<ValueAst>,
    pub accepted_pairs: Option<ValueAst>,
    pub aromatic_valence: Option<AromaticExpr>,
    pub multicenter_valence: Option<ValueAst>,
}

impl AtomAst {
    fn new(element: ElementExpr) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementExpr::Lit(element))
    }
}

impl LowerAst for AtomAst {
    type Config = AtomLowerConfig;
}

/// Isotope interpretation mode
#[derive(Clone, Debug, Default)]
pub enum IsotopeMode {
    #[default]
    Normal,
    Provided,
}

/// Charge interpretation mode
#[derive(Clone, Debug, Default)]
pub enum ChargeMode {
    Zero,
    #[default]
    Provided,
}

/// Implicit hydrogen interpretation mode
#[derive(Clone, Debug, Default)]
pub enum ImplicitHydrogenMode {
    Zero,
    Normal,
    #[default]
    Provided,
}

/// Aromatic interpretation mode
#[derive(Clone, Debug, Default)]
pub enum AromaticMode {
    None,
    Any,
    #[default]
    Provided,
}

/// Atom lowering configuration
#[derive(Clone, Debug, Default)]
pub struct AtomLowerConfig {
    pub isotope_mode: IsotopeMode,
    pub charge_mode: ChargeMode,
    pub implicit_h_mode: ImplicitHydrogenMode,
    pub aromatic_mode: AromaticMode,
}

/// Parse a complete atom-string
pub fn parse_atom_dsl(input: &str) -> Result<AtomAst, ParseError> {
    all_consuming(atom_dsl)
        .parse(input)
        .map(|(_, r)| r)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        })
}

/// Atom-string parser (does not require consuming all input)
pub fn atom_dsl(i: &str) -> IResult<&str, AtomAst, ParseError> {
    let (remaining, (element, preds)) = pair(
        delimited(multispace0, element_expr, multispace0),
        many0(terminated(atom_predicate, multispace0)),
    )
    .parse(i)?;

    let mut ast = AtomAst::new(element);
    update_atom_ast(&mut ast, preds).map_err(|e| Err::Error(e))?;
    Ok((remaining, ast))
}

fn update_atom_ast(ast: &mut AtomAst, preds: Vec<AtomPredicate>) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(v) => {
                if ast.isotope_mass.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#i".to_string()));
                }
                ast.isotope_mass = Some(v);
            }
            AtomPredicate::Charge(v) => {
                if ast.charge.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#c".to_string()));
                }
                ast.charge = Some(v);
            }
            AtomPredicate::ImplicitHydrogens(v) => {
                if ast.implicit_hydrogens.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#h".to_string()));
                }
                ast.implicit_hydrogens = Some(v);
            }
            AtomPredicate::LonePairs(v) => {
                if ast.lone_pairs.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#n".to_string()));
                }
                ast.lone_pairs = Some(v);
            }
            AtomPredicate::UnpairedElectrons(v) => {
                if ast.unpaired_electrons.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#u".to_string()));
                }
                ast.unpaired_electrons = Some(v);
            }
            AtomPredicate::Multiplicity(v) => {
                if ast.multiplicity.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#s".to_string()));
                }
                ast.multiplicity = Some(v);
            }
            AtomPredicate::Valence(v) => {
                if ast.valence.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#v".to_string()));
                }
                ast.valence = Some(v);
            }
            AtomPredicate::DonatedPairs(v) => {
                if ast.donated_pairs.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#d".to_string()));
                }
                ast.donated_pairs = Some(v);
            }
            AtomPredicate::AcceptedPairs(v) => {
                if ast.accepted_pairs.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#r".to_string()));
                }
                ast.accepted_pairs = Some(v);
            }
            AtomPredicate::AromaticValence(v) => {
                if ast.aromatic_valence.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#a".to_string()));
                }
                ast.aromatic_valence = Some(v);
            }
            AtomPredicate::MulticenterValence(v) => {
                if ast.multicenter_valence.is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#m".to_string()));
                }
                ast.multicenter_valence = Some(v);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use nom::Err;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::dsl::error::ParseError;
    use crate::dsl::value::{Expr, RelOp, ValueAst};

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomAst::new(ElementExpr::Lit(Element::C)))]
    #[case::iron("Fe", AtomAst::new(ElementExpr::Lit(Element::Fe)))]
    #[case::chlorine("Cl", AtomAst::new(ElementExpr::Lit(Element::Cl)))]
    #[case::whitespace("  C  ", AtomAst::new(ElementExpr::Lit(Element::C)))]
    #[case::wildcard("*", AtomAst::new(ElementExpr::Wildcard))]
    #[case::isotope("C#i12", AtomAst { isotope_mass: Some(IsotopeExpr::Lit(12)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::isotope_natural("C#i=", AtomAst { isotope_mass: Some(IsotopeExpr::Natural), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_pos("C#c+2", AtomAst { charge: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_neg("C#c-2", AtomAst { charge: Some(ValueAst::Lit(-2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_plus("C#c+", AtomAst { charge: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_minus("C#c-", AtomAst { charge: Some(ValueAst::Lit(-1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_zero("C#c0", AtomAst { charge: Some(ValueAst::Lit(0)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_count("C#h3", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_normal("C#h=", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Normal), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_wild("C#h*", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_omit("C#h", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::lone_pairs("O#n2", AtomAst { lone_pairs: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::O)) })]
    #[case::lone_pairs_omit("O#n", AtomAst { lone_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::O)) })]
    #[case::unpaired("C#u2", AtomAst { unpaired_electrons: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::unpaired_omit("C#u", AtomAst { unpaired_electrons: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::multiplicity("C#s3", AtomAst { multiplicity: Some(ValueAst::Lit(3)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::multiplicity_omit("C#s", AtomAst { multiplicity: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::valence("C#v4", AtomAst { valence: Some(ValueAst::Lit(4)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::donated_pairs("N#d1", AtomAst { donated_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::N)) })]
    #[case::accepted_pairs("B#r1", AtomAst { accepted_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::B)) })]
    #[case::arom_nonmember("C#a!", AtomAst { aromatic_valence: Some(AromaticExpr::None), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_wild("C#a*", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_zero("C#a0", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(0))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_one("C#a1", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_omit("C#a", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::multicenter("C#m2", AtomAst { multicenter_valence: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::ch3v1("C#h3#v1", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))), valence: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::nh2n1v3("N#h2#n1#v3", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(2))), lone_pairs: Some(ValueAst::Lit(1)), valence: Some(ValueAst::Lit(3)), ..AtomAst::new(ElementExpr::Lit(Element::N)) })]
    #[case::fe_highspin("Fe#c+2#u4#s5", AtomAst { charge: Some(ValueAst::Lit(2)), unpaired_electrons: Some(ValueAst::Lit(4)), multiplicity: Some(ValueAst::Lit(5)), ..AtomAst::new(ElementExpr::Lit(Element::Fe)) })]
    #[case::h_expr("C#h(?h >= 1)", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    fn test_parse_atom_dsl(#[case] input: &str, #[case] expected: AtomAst) {
        let result = atom_dsl(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, ast) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
        assert_eq!(ast, expected);
    }

    #[rstest]
    #[case::empty("", ParseError::InvalidElement("".to_string()))]
    #[case::no_element("#h3", ParseError::InvalidElement("#h3".to_string()))]
    #[case::unknown_pred("C#x", ParseError::UnknownAtomPredicate("#x".to_string()))]
    #[case::dup_h("C#h3#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_charge("C#c+#c-", ParseError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::dup_valence("C#v3#v4", ParseError::DuplicateAtomPredicate("#v".to_string()))]
    #[case::trailing("C#h3 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_atom_dsl_invalid(#[case] input: &str, #[case] expected: ParseError) {
        let result = atom_dsl(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        };
        assert_eq!(
            err, expected,
            "{input:?} should fail with {expected:?}, got {err:?}"
        );
    }
}
