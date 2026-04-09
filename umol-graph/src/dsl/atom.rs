//! Atom-string DSL: parser, AST, and display

use std::fmt::{self, Display};
use std::str::FromStr;

use nom::character::complete::multispace0;
use nom::combinator::all_consuming;
use nom::multi::many0;
use nom::sequence::{delimited, pair, terminated};
use nom::{Err, IResult, Parser};
use serde::de::{Deserializer, Error as DeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::Element;

use super::ast::DslAst;
use super::config::AtomDslConfig;
use super::error::ParseError;
use super::predicates::{
    atom_predicate, element_expr, AromaticExpr, AtomPredicate, ElementExpr, HydrogenExpr,
    IsotopeExpr,
};
use super::value::ValueAst;

/// Parsed atom-string AST
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
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

impl DslAst for AtomAst {
    type Config = AtomDslConfig;
}

impl FromStr for AtomAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom_dsl(s)
    }
}

impl Display for AtomAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_element(f, &self.element)?;

        match &self.isotope_mass {
            None => {}
            Some(IsotopeExpr::Natural) => write!(f, "#i=")?,
            Some(IsotopeExpr::Lit(n)) => write!(f, "#i{}", n)?,
            Some(IsotopeExpr::Wildcard) => write!(f, "#i*")?,
            Some(IsotopeExpr::Set(ns)) => {
                write!(f, "#i{{")?;
                for (i, n) in ns.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}}")?;
            }
            Some(IsotopeExpr::Bind { id, set }) => {
                write!(f, "#i(?{} :: {{", id)?;
                for (i, n) in set.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}})")?;
            }
            Some(IsotopeExpr::Ref(id)) => write!(f, "#i(?{})", id)?,
        }

        match &self.charge {
            None | Some(ValueAst::Lit(0)) => {}
            Some(ValueAst::Lit(1)) => write!(f, "#c+")?,
            Some(ValueAst::Lit(-1)) => write!(f, "#c-")?,
            Some(ValueAst::Lit(n)) if *n > 0 => write!(f, "#c+{}", n)?,
            Some(ValueAst::Lit(n)) => write!(f, "#c{}", n)?,
            Some(ValueAst::Wildcard) => write!(f, "#c*")?,
            Some(v) => {
                write!(f, "#c")?;
                fmt_value(f, v)?;
            }
        }

        match &self.implicit_hydrogens {
            None | Some(HydrogenExpr::Value(ValueAst::Lit(0))) => {}
            Some(HydrogenExpr::Normal) => write!(f, "#h=")?,
            Some(HydrogenExpr::Value(ValueAst::Lit(1))) => write!(f, "#h")?,
            Some(HydrogenExpr::Value(ValueAst::Lit(n))) => write!(f, "#h{}", n)?,
            Some(HydrogenExpr::Value(ValueAst::Wildcard)) => write!(f, "#h*")?,
            Some(HydrogenExpr::Value(v)) => {
                write!(f, "#h")?;
                fmt_value(f, v)?;
            }
        }

        fmt_unsigned_field(f, "#n", &self.lone_pairs)?;
        fmt_unsigned_field(f, "#u", &self.unpaired_electrons)?;
        fmt_multiplicity(f, &self.multiplicity, &self.unpaired_electrons)?;
        fmt_unsigned_field(f, "#v", &self.valence)?;
        fmt_unsigned_field(f, "#d", &self.donated_pairs)?;
        fmt_unsigned_field(f, "#r", &self.accepted_pairs)?;

        match &self.aromatic_valence {
            None => {}
            Some(AromaticExpr::Unspecified) => write!(f, "#a?")?,
            Some(AromaticExpr::NotAromatic) => write!(f, "#a!")?,
            Some(AromaticExpr::Value(ValueAst::Lit(1))) => write!(f, "#a")?,
            Some(AromaticExpr::Value(ValueAst::Lit(n))) => write!(f, "#a{}", n)?,
            Some(AromaticExpr::Value(v)) => {
                write!(f, "#a")?;
                fmt_value(f, v)?;
            }
        }

        fmt_unsigned_field(f, "#m", &self.multicenter_valence)?;

        Ok(())
    }
}

impl Serialize for AtomAst {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for AtomAst {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        parse_atom_dsl(&s).map_err(DeError::custom)
    }
}

impl<'de> umol_edn::FromEdn<'de> for AtomAst {
    fn from_edn(edn: &umol_edn::Edn<'de>) -> Result<Self, umol_edn::DeError> {
        match edn {
            umol_edn::Edn::Str(s) => {
                parse_atom_dsl(s).map_err(|e| umol_edn::DeError::Custom(e.to_string()))
            }
            other => Err(umol_edn::DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl umol_edn::ToEdn for AtomAst {
    fn to_edn(&self) -> umol_edn::Edn<'_> {
        umol_edn::Edn::Str(std::borrow::Cow::Owned(self.to_string()))
    }
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
    update_atom_ast(&mut ast, preds).map_err(Err::Error)?;
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

fn fmt_element(f: &mut fmt::Formatter<'_>, expr: &ElementExpr) -> fmt::Result {
    match expr {
        ElementExpr::Lit(e) => write!(f, "{}", e),
        ElementExpr::Wildcard => write!(f, "*"),
        ElementExpr::Set(es) => {
            write!(f, "{{")?;
            for (i, e) in es.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", e)?;
            }
            write!(f, "}}")
        }
        ElementExpr::Bind { id, set } => {
            write!(f, "(?{} :: {{", id)?;
            for (i, e) in set.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", e)?;
            }
            write!(f, "}})")
        }
        ElementExpr::Ref(id) => write!(f, "(?{})", id),
    }
}

/// Suppress None and Lit(0); abbreviate Lit(1) to just the prefix.
fn fmt_unsigned_field(
    f: &mut fmt::Formatter<'_>,
    prefix: &str,
    v: &Option<ValueAst>,
) -> fmt::Result {
    match v {
        None | Some(ValueAst::Lit(0)) => Ok(()),
        Some(ValueAst::Lit(1)) => write!(f, "{}", prefix),
        Some(ValueAst::Lit(n)) => write!(f, "{}{}", prefix, n),
        Some(ValueAst::Wildcard) => write!(f, "{}*", prefix),
        Some(v) => {
            write!(f, "{}", prefix)?;
            fmt_value(f, v)
        }
    }
}

/// Suppress multiplicity when it equals unpaired_electrons + 1 (derivable default).
fn fmt_multiplicity(
    f: &mut fmt::Formatter<'_>,
    multiplicity: &Option<ValueAst>,
    unpaired: &Option<ValueAst>,
) -> fmt::Result {
    let m = match multiplicity {
        None => return Ok(()),
        Some(ValueAst::Lit(m)) => *m,
        Some(ValueAst::Wildcard) => return write!(f, "#s*"),
        Some(v) => {
            write!(f, "#s")?;
            return fmt_value(f, v);
        }
    };
    let u: i32 = match unpaired {
        Some(ValueAst::Lit(u)) => *u,
        None => 0,
        _ => -1, // non-literal: can't determine derivability, always print
    };
    if m as i32 == u + 1 {
        Ok(()) // derivable from unpaired, suppress
    } else if m == 1 {
        write!(f, "#s")
    } else {
        write!(f, "#s{}", m)
    }
}

fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Wildcard => write!(f, "*"),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => {
            write!(f, "{{")?;
            for (i, n) in s.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", n)?;
            }
            write!(f, "}}")
        }
        ValueAst::Expr(_) => write!(f, "<expr>"),
    }
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
    #[case::element_set("{C,N,O}", AtomAst::new(ElementExpr::Set(vec![Element::C, Element::N, Element::O])))]
    #[case::element_bind("(?e :: {C,N})", AtomAst::new(ElementExpr::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] }))]
    #[case::element_ref("(?e)", AtomAst::new(ElementExpr::Ref("e".to_string())))]
    #[case::isotope("C#i12", AtomAst { isotope_mass: Some(IsotopeExpr::Lit(12)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::isotope_natural("C#i=", AtomAst { isotope_mass: Some(IsotopeExpr::Natural), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_pos("C#c+2", AtomAst { charge: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_neg("C#c-2", AtomAst { charge: Some(ValueAst::Lit(-2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_plus("C#c+", AtomAst { charge: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_minus("C#c-", AtomAst { charge: Some(ValueAst::Lit(-1)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::charge_zero("C#c0", AtomAst { charge: Some(ValueAst::Lit(0)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_count("C#h3", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Lit(3))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_normal("C#h=", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Normal), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_wildcard("C#h*", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_bind("C#h(?h)", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Expr(Expr::Var("h".to_string())))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::h_set("N#h?h :: {2,3}", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![2, 3])))), ..AtomAst::new(ElementExpr::Lit(Element::N)) })]
    #[case::h_expr("C#h?h >= 1", AtomAst { implicit_hydrogens: Some(HydrogenExpr::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
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
    #[case::arom_unspecified("C#a?", AtomAst { aromatic_valence: Some(AromaticExpr::Unspecified), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_not_aromatic("C#a!", AtomAst { aromatic_valence: Some(AromaticExpr::NotAromatic), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_aromatic("C#a*", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_zero("C#a0", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(0))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_one("C#a1", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::arom_omit("C#a", AtomAst { aromatic_valence: Some(AromaticExpr::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
    #[case::multicenter("C#m2", AtomAst { multicenter_valence: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementExpr::Lit(Element::C)) })]
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
