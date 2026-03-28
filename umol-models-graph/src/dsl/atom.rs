//! Atom-string DSL parser — `spec/umol-dsl-spec.md` §7.3 / §7.4.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, multispace0, satisfy, u32 as nom_u32};
use nom::combinator::{all_consuming, map, recognize, value};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};
use umol_data::Element;

use super::error::ParseError;
use super::predicates::{atom_predicate, AtomPredicate};
use super::value::{op_char, parse_id, ValueAst};

/// Element position of an atom-string.
///
/// Categorical: supports literal, wildcard, set, bind, and reference.
/// No arithmetic; no numeric ordering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ElementExpr {
    Lit(Element),
    Wildcard,
    Set(Vec<Element>),
    Bind { id: String, set: Vec<Element> },
    Ref(String),
}

impl ElementExpr {
    pub fn new(element: Element) -> Self {
        Self::Lit(element)
    }
}

/// Isotope-mass position (`#i` payload).
///
/// Categorical like `ElementExpr`: indexed by mass number, no arithmetic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsotopeExpr {
    Lit(u32),
    Wildcard,
    Set(Vec<u32>),
    Bind { id: String, set: Vec<u32> },
    Ref(String),
}

/// `#h` payload — `=` (valence-model H) is not a `ValueAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum HydrogenExpr {
    /// `#h=` — implicit H count from the valence model.
    Normal,
    Value(ValueAst),
}

impl HydrogenExpr {
    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }
}

/// `#a` payload — `!` (non-member) is not a `ValueAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticExpr {
    /// `#a!` — atom is not a member of any aromatic system.
    None,
    Value(ValueAst),
}

/// Parsed atom-string AST.
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

fn element_literal(i: &str) -> IResult<&str, Element, ParseError> {
    let (rest, sym) = recognize(pair(
        satisfy(|c: char| c.is_ascii_uppercase()),
        many0(satisfy(|c: char| c.is_ascii_lowercase())),
    ))
    .parse(i)?;
    match Element::from_symbol(sym) {
        Some(el) => Ok((rest, el)),
        None => Err(Err::Error(ParseError::InvalidElement(sym.to_string()))),
    }
}

fn element_set(i: &str) -> IResult<&str, Vec<Element>, ParseError> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), element_literal),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn element_bind(i: &str) -> IResult<&str, (String, Vec<Element>), ParseError> {
    delimited(
        char('('),
        pair(
            delimited(multispace0, preceded(char('?'), parse_id), multispace0),
            preceded(
                pair(tag("::"), multispace0),
                terminated(element_set, multispace0),
            ),
        ),
        char(')'),
    )
    .parse(i)
}

fn element_ref(i: &str) -> IResult<&str, String, ParseError> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

pub(crate) fn element_expr(i: &str) -> IResult<&str, ElementExpr, ParseError> {
    alt((
        value(ElementExpr::Wildcard, char('*')),
        map(element_set, ElementExpr::Set),
        map(element_bind, |(id, set)| ElementExpr::Bind { id, set }),
        map(element_ref, ElementExpr::Ref),
        map(element_literal, ElementExpr::Lit),
    ))
    .parse(i)
    .map_err(|_| Err::Error(ParseError::InvalidElement(i.to_string())))
}

fn isotope_set(i: &str) -> IResult<&str, Vec<u32>, ParseError> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), nom_u32),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn isotope_bind(i: &str) -> IResult<&str, (String, Vec<u32>), ParseError> {
    delimited(
        char('('),
        pair(
            delimited(multispace0, preceded(char('?'), parse_id), multispace0),
            preceded(
                pair(tag("::"), multispace0),
                terminated(isotope_set, multispace0),
            ),
        ),
        char(')'),
    )
    .parse(i)
}

fn isotope_ref(i: &str) -> IResult<&str, String, ParseError> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

pub(crate) fn isotope_expr(i: &str) -> IResult<&str, IsotopeExpr, ParseError> {
    alt((
        value(IsotopeExpr::Wildcard, char('*')),
        map(isotope_set, IsotopeExpr::Set),
        map(isotope_bind, |(id, set)| IsotopeExpr::Bind { id, set }),
        map(isotope_ref, IsotopeExpr::Ref),
        map(nom_u32, IsotopeExpr::Lit),
    ))
    .parse(i)
}

/// Parse a complete atom-string.
pub fn parse_atom_dsl(input: &str) -> Result<AtomAst, ParseError> {
    all_consuming(atom_dsl)
        .parse(input)
        .map(|(_, r)| r)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        })
}

/// Atom-string parser (does not require consuming all input).
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
    #[case::carbon("C", ElementExpr::Lit(Element::C))]
    #[case::iron("Fe", ElementExpr::Lit(Element::Fe))]
    #[case::chlorine("Cl", ElementExpr::Lit(Element::Cl))]
    #[case::wildcard("*", ElementExpr::Wildcard)]
    #[case::set("{C,N,O}", ElementExpr::Set(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementExpr::Set(vec![Element::C, Element::N]))]
    #[case::bind("(?e :: {C,N})", ElementExpr::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })]
    #[case::ref_("(?e)", ElementExpr::Ref("e".to_string()))]
    fn test_element_expr(#[case] input: &str, #[case] expected: ElementExpr) {
        let result = element_expr(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, expr) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
        assert_eq!(expr, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::lowercase("c")]
    #[case::invalid("123")]
    #[case::unknown_element("Xx")]
    fn test_element_expr_invalid(#[case] input: &str) {
        let result = element_expr(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit("12", IsotopeExpr::Lit(12))]
    #[case::wildcard("*", IsotopeExpr::Wildcard)]
    #[case::set("{12,13,14}", IsotopeExpr::Set(vec![12, 13, 14]))]
    #[case::bind("(?m :: {12,13})", IsotopeExpr::Bind { id: "m".to_string(), set: vec![12, 13] })]
    #[case::ref_("(?m)", IsotopeExpr::Ref("m".to_string()))]
    fn test_isotope_expr(#[case] input: &str, #[case] expected: IsotopeExpr) {
        let result = isotope_expr(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, expr) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
        assert_eq!(expr, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomAst::new(ElementExpr::Lit(Element::C)))]
    #[case::iron("Fe", AtomAst::new(ElementExpr::Lit(Element::Fe)))]
    #[case::chlorine("Cl", AtomAst::new(ElementExpr::Lit(Element::Cl)))]
    #[case::whitespace("  C  ", AtomAst::new(ElementExpr::Lit(Element::C)))]
    #[case::wildcard("*", AtomAst::new(ElementExpr::Wildcard))]
    #[case::isotope("C#i12", AtomAst { isotope_mass: Some(IsotopeExpr::Lit(12)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
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
    fn test_atom_dsl(#[case] input: &str, #[case] expected: AtomAst) {
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
    fn test_atom_dsl_invalid(#[case] input: &str, #[case] expected: ParseError) {
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
