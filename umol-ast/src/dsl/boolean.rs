//! Boolean DSL: parser, `Display`, EDN boundary.

use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::combinator::{alt, opt};
use winnow::Parser;

use super::edn_utils::eof_err;
use super::error::{PResult, ParseError};
use crate::ast::boolean::BooleanAst;
use crate::ast::traits::{FromAst, IntoAst};

/// Boundary type for [`BooleanAst`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BooleanDsl(pub BooleanAst);

impl FromAst<BooleanAst> for BooleanDsl {
    type Ctx = ();

    fn from_ast(ast: &BooleanAst, _ctx: &Self::Ctx) -> Self {
        Self(*ast)
    }
}

impl IntoAst<BooleanAst> for BooleanDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> BooleanAst {
        self.0
    }
}

/// Combinator: `!` → false, `*` → undetermined, `+`/(absent) → true.
pub(crate) fn boolean(i: &mut &str) -> PResult<BooleanDsl> {
    opt(alt((
        '!'.value(BooleanAst::Lit(false)),
        '*'.value(BooleanAst::Undetermined),
        '+'.value(BooleanAst::Lit(true)),
    )))
    .map(|sign| BooleanDsl(sign.unwrap_or(BooleanAst::Lit(true))))
    .parse_next(i)
}

pub fn parse_boolean(input: &str) -> Result<BooleanDsl, ParseError> {
    boolean.parse(input).map_err(|e| e.into_inner())
}

/// Streaming reader for the EDN form: `true` / `false` / `:undetermined` / a string (`"!"`/`"*"`).
pub(crate) fn read_boolean_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<BooleanDsl, EdnError> {
    let b = match de.peek_byte()?.ok_or_else(eof_err)? {
        b':' => {
            let name = de.read_keyword_name()?;
            if name.as_ref() == "undetermined" {
                BooleanAst::Undetermined
            } else {
                return Err(DeError::Custom(format!("unknown boolean keyword :{name}")).into());
            }
        }
        b'"' => {
            parse_boolean(de.read_string()?.as_ref())
                .map_err(|e| DeError::subgrammar("boolean", e))?
                .0
        }
        _ => match de.read_value_slice()? {
            "true" => BooleanAst::Lit(true),
            "false" => BooleanAst::Lit(false),
            other => return Err(DeError::Custom(format!("expected boolean, got {other}")).into()),
        },
    };
    Ok(BooleanDsl(b))
}

pub(crate) fn fmt_boolean(f: &mut fmt::Formatter<'_>, b: &BooleanAst) -> fmt::Result {
    match b {
        BooleanAst::Lit(true) => Ok(()),
        BooleanAst::Lit(false) => write!(f, "!"),
        BooleanAst::Undetermined => write!(f, "*"),
    }
}

impl Display for BooleanDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_boolean(f, &self.0)
    }
}

impl FromStr for BooleanDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_boolean(s)
    }
}

impl<'de> FromEdn<'de> for BooleanDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let b = match edn {
            Edn::Bool(b) => BooleanAst::Lit(*b),
            Edn::Keyword(k) if k.name() == "undetermined" => BooleanAst::Undetermined,
            Edn::Str(s) => {
                parse_boolean(s)
                    .map_err(|e| DeError::subgrammar("boolean", e))?
                    .0
            }
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "boolean (bool, :undetermined, or string)",
                    got: other.kind(),
                    path: Vec::new(),
                })
            }
        };
        Ok(Self(b))
    }
}

impl ToEdn for BooleanDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self.0 {
            BooleanAst::Lit(b) => Edn::Bool(b),
            BooleanAst::Undetermined => Edn::Keyword(EdnKeyword::owned("undetermined".to_string())),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;

    #[rstest]
    #[case::plus("+", BooleanDsl(BooleanAst::Lit(true)))]
    #[case::empty("", BooleanDsl(BooleanAst::Lit(true)))]
    #[case::bang("!", BooleanDsl(BooleanAst::Lit(false)))]
    #[case::star("*", BooleanDsl(BooleanAst::Undetermined))]
    fn test_parse_boolean(#[case] input: &str, #[case] expected: BooleanDsl) {
        assert_eq!(parse_boolean(input).unwrap(), expected);
    }

    #[rstest]
    #[case::trailing("!!")]
    #[case::unknown("x")]
    fn test_parse_boolean_error(#[case] input: &str) {
        assert!(parse_boolean(input).is_err());
    }

    #[rstest]
    #[case::truthy(BooleanDsl(BooleanAst::Lit(true)), "")]
    #[case::falsy(BooleanDsl(BooleanAst::Lit(false)), "!")]
    #[case::undetermined(BooleanDsl(BooleanAst::Undetermined), "*")]
    fn test_boolean_dsl_display(#[case] dsl: BooleanDsl, #[case] expected: &str) {
        assert_eq!(dsl.to_string(), expected);
    }

    #[rstest]
    #[case::bool_true("true", BooleanDsl(BooleanAst::Lit(true)))]
    #[case::bool_false("false", BooleanDsl(BooleanAst::Lit(false)))]
    #[case::keyword(":undetermined", BooleanDsl(BooleanAst::Undetermined))]
    #[case::string_bang(r##""!""##, BooleanDsl(BooleanAst::Lit(false)))]
    #[case::string_star(r##""*""##, BooleanDsl(BooleanAst::Undetermined))]
    fn test_boolean_dsl_from_edn(#[case] input: &str, #[case] expected: BooleanDsl) {
        assert_eq!(
            BooleanDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::truthy(BooleanDsl(BooleanAst::Lit(true)), Edn::Bool(true))]
    #[case::falsy(BooleanDsl(BooleanAst::Lit(false)), Edn::Bool(false))]
    #[case::undetermined(BooleanDsl(BooleanAst::Undetermined), Edn::Keyword(EdnKeyword::owned("undetermined".to_string())))]
    fn test_boolean_dsl_to_edn(#[case] dsl: BooleanDsl, #[case] expected: Edn<'static>) {
        assert_eq!(dsl.to_edn(), expected);
    }
}
