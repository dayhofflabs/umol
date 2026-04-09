//! Owned, serde-friendly EDN value.
//!
//! [`Value`] is an owned wrapper around [`Edn<'static>`] that carries the
//! full variant information of the internal EDN data model across the serde
//! boundary. Unlike `Edn<'a>`, which is borrow-parameterized to support
//! zero-copy parsing, `Value` has no lifetime and can be stored in struct
//! fields, sent across threads, and used as a dynamic-typed escape hatch
//! inside otherwise-typed deserialize paths.
//!
//! The EDN data model has variants (keyword, symbol, list, set, tagged,
//! bigint, bigdecimal) that do not correspond to serde's primitive set.
//! Under the `serde` feature, `Value` round-trips them losslessly through
//! the EDN serializer/deserializer while still producing sensible fallbacks
//! when the payload travels over a foreign format such as JSON or YAML.
//!
//! # Serialize
//!
//! The `Serialize` impl walks the inner `Edn` tree and, for each variant,
//! dispatches through the existing wrapper token (`KEYWORD_TOKEN`,
//! `SYMBOL_TOKEN`, `LIST_TOKEN`, `SET_TOKEN`, `TAGGED_TOKEN`,
//! `BIGINT_TOKEN`, `BIGDECIMAL_TOKEN`). Non-EDN serializers ignore the
//! tokens and see an ordinary string / sequence / tuple / map.
//!
//! # Deserialize
//!
//! The `Deserialize` impl calls `deserialize_newtype_struct(VALUE_TOKEN,
//! ValueVisitor)`. `EdnDeserializer` recognizes `VALUE_TOKEN` and, for
//! EDN-specific variants, hands a `ValueCarrier` deserializer to
//! `visit_newtype_struct`. The carrier presents itself through
//! `visit_enum` with variant name = token and payload = inner value,
//! which is the only visit method not used by any standard format when
//! responding to `deserialize_any`. Non-EDN deserializers default to
//! `visit_newtype_struct(self)`, and `ValueVisitor::visit_newtype_struct`
//! falls back to `deserialize_any(self)` so foreign formats degrade to
//! the lossy mapping (keyword/symbol → `Str`, list/set → `Vector`,
//! tagged → tuple).

use std::fmt;
use std::str::FromStr;

#[cfg(feature = "serde")]
use std::borrow::Cow;

#[cfg(feature = "serde")]
use serde::de::{
    self, DeserializeSeed, Deserializer, EnumAccess, MapAccess, SeqAccess, VariantAccess, Visitor,
};
#[cfg(feature = "serde")]
use serde::ser::{Serialize, SerializeMap, SerializeSeq, SerializeTupleStruct, Serializer};
#[cfg(feature = "serde")]
use serde::Deserialize;

use crate::config::ParseConfig;
use crate::edn::Edn;
use crate::error::EdnError;
use crate::native::{FromEdn, ToEdn};
use crate::reader::read_string_with;

#[cfg(feature = "serde")]
use crate::collections::{EdnMap, EdnSeq, EdnSet};
#[cfg(feature = "serde")]
use crate::de::EdnDeserializer;
#[cfg(feature = "serde")]
use crate::edn::{Keyword, Symbol};
#[cfg(feature = "serde")]
use crate::serde_tokens::{
    KEYWORD_TOKEN, LIST_TOKEN, SET_TOKEN, SYMBOL_TOKEN, TAGGED_TOKEN, VALUE_TOKEN,
};
#[cfg(all(feature = "serde", feature = "bignum"))]
use crate::serde_tokens::{BIGDECIMAL_TOKEN, BIGINT_TOKEN};

#[cfg(all(feature = "serde", feature = "bignum"))]
use crate::edn::EdnBigDecimal;

/// Owned EDN value. A lossless mirror of [`Edn<'static>`] for use where a
/// fully owned, lifetime-free value is required.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Value(pub(crate) Edn<'static>);

impl Value {
    /// Wrap an already-owned `Edn<'static>` value.
    pub fn new(edn: Edn<'static>) -> Self {
        Self(edn)
    }

    /// Parse an EDN string into a `Value` using the default config.
    ///
    /// Unknown tags are rejected with `InvalidTag`. Callers that need to
    /// preserve arbitrary tagged literals must supply a config with
    /// `allow_unknown_tags = true` via [`Value::parse_with`] or register
    /// the relevant tag readers.
    pub fn parse(input: &str) -> Result<Self, EdnError> {
        Self::parse_with(input, &ParseConfig::default())
    }

    /// Parse an EDN string into a `Value` using the supplied config.
    pub fn parse_with(input: &str, config: &ParseConfig) -> Result<Self, EdnError> {
        Ok(Self(read_string_with(input, config)?.into_owned()))
    }

    /// Borrow the inner `Edn<'static>`.
    pub fn as_edn(&self) -> &Edn<'static> {
        &self.0
    }

    /// Unwrap into the inner `Edn<'static>`.
    pub fn into_edn(self) -> Edn<'static> {
        self.0
    }
}

impl From<Edn<'static>> for Value {
    fn from(edn: Edn<'static>) -> Self {
        Self(edn)
    }
}

impl From<&Edn<'_>> for Value {
    fn from(edn: &Edn<'_>) -> Self {
        Self(edn.clone().into_owned())
    }
}

impl From<Value> for Edn<'static> {
    fn from(v: Value) -> Self {
        v.0
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for Value {
    type Err = EdnError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

impl<'de> FromEdn<'de> for Value {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError> {
        Ok(Self(edn.clone().into_owned()))
    }
}

impl ToEdn for Value {
    fn to_edn(&self) -> Edn<'_> {
        self.0.clone()
    }
}

#[cfg(feature = "serde")]
impl Serialize for Value {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        EdnRef(&self.0).serialize(serializer)
    }
}

/// Borrowing helper that serializes any `&Edn` through the same
/// variant-preserving token dispatch used by `Value`.
#[cfg(feature = "serde")]
struct EdnRef<'a>(&'a Edn<'a>);

#[cfg(feature = "serde")]
impl<'a> Serialize for EdnRef<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self.0 {
            Edn::Nil => serializer.serialize_unit(),
            Edn::Bool(b) => serializer.serialize_bool(*b),
            Edn::Int(i) => serializer.serialize_i64(*i),
            Edn::Float(f) => serializer.serialize_f64(*f),
            Edn::Char(c) => serializer.serialize_char(*c),
            Edn::Str(s) => serializer.serialize_str(s),
            Edn::Keyword(k) => serializer.serialize_newtype_struct(KEYWORD_TOKEN, k.as_str()),
            Edn::Symbol(s) => serializer.serialize_newtype_struct(SYMBOL_TOKEN, s.as_str()),
            Edn::Vector(v) => {
                let mut seq = serializer.serialize_seq(Some(v.len()))?;
                for item in v.iter() {
                    seq.serialize_element(&EdnRef(item))?;
                }
                seq.end()
            }
            Edn::List(v) => serializer.serialize_newtype_struct(LIST_TOKEN, &ListPayload(v)),
            Edn::Set(s) => serializer.serialize_newtype_struct(SET_TOKEN, &SetPayload(s)),
            Edn::Map(m) => {
                let mut map = serializer.serialize_map(Some(m.len()))?;
                for (k, v) in m.iter() {
                    map.serialize_entry(&EdnRef(k), &EdnRef(v))?;
                }
                map.end()
            }
            Edn::Tagged(tag, inner) => {
                let mut ts = serializer.serialize_tuple_struct(TAGGED_TOKEN, 2)?;
                ts.serialize_field(tag.as_ref())?;
                ts.serialize_field(&EdnRef(inner))?;
                ts.end()
            }
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => serializer.serialize_newtype_struct(BIGINT_TOKEN, &n.to_string()),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => {
                serializer.serialize_newtype_struct(BIGDECIMAL_TOKEN, &d.as_inner().to_string())
            }
        }
    }
}

#[cfg(feature = "serde")]
struct ListPayload<'a>(&'a EdnSeq<'a>);

#[cfg(feature = "serde")]
impl<'a> Serialize for ListPayload<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for item in self.0.iter() {
            seq.serialize_element(&EdnRef(item))?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
struct SetPayload<'a>(&'a EdnSet<'a>);

#[cfg(feature = "serde")]
impl<'a> Serialize for SetPayload<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for item in self.0.iter() {
            seq.serialize_element(&EdnRef(item))?;
        }
        seq.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for Value {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(VALUE_TOKEN, ValueVisitor)
    }
}

#[cfg(feature = "serde")]
pub(crate) struct ValueVisitor;

#[cfg(feature = "serde")]
impl<'de> Visitor<'de> for ValueVisitor {
    type Value = Value;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("any EDN value")
    }

    fn visit_unit<E: de::Error>(self) -> Result<Value, E> {
        Ok(Value(Edn::Nil))
    }

    fn visit_none<E: de::Error>(self) -> Result<Value, E> {
        Ok(Value(Edn::Nil))
    }

    fn visit_some<D: Deserializer<'de>>(self, deserializer: D) -> Result<Value, D::Error> {
        Value::deserialize(deserializer)
    }

    fn visit_bool<E: de::Error>(self, v: bool) -> Result<Value, E> {
        Ok(Value(Edn::Bool(v)))
    }

    fn visit_i64<E: de::Error>(self, v: i64) -> Result<Value, E> {
        Ok(Value(Edn::Int(v)))
    }

    fn visit_i128<E: de::Error>(self, v: i128) -> Result<Value, E> {
        i64::try_from(v)
            .map(|i| Value(Edn::Int(i)))
            .map_err(|_| E::custom(format!("i128 {v} out of range for Edn::Int")))
    }

    fn visit_u64<E: de::Error>(self, v: u64) -> Result<Value, E> {
        i64::try_from(v)
            .map(|i| Value(Edn::Int(i)))
            .map_err(|_| E::custom(format!("u64 {v} out of range for Edn::Int")))
    }

    fn visit_u128<E: de::Error>(self, v: u128) -> Result<Value, E> {
        i64::try_from(v)
            .map(|i| Value(Edn::Int(i)))
            .map_err(|_| E::custom(format!("u128 {v} out of range for Edn::Int")))
    }

    fn visit_f64<E: de::Error>(self, v: f64) -> Result<Value, E> {
        Ok(Value(Edn::Float(v)))
    }

    fn visit_char<E: de::Error>(self, v: char) -> Result<Value, E> {
        Ok(Value(Edn::Char(v)))
    }

    fn visit_str<E: de::Error>(self, v: &str) -> Result<Value, E> {
        Ok(Value(Edn::Str(Cow::Owned(v.to_string()))))
    }

    fn visit_borrowed_str<E: de::Error>(self, v: &'de str) -> Result<Value, E> {
        Ok(Value(Edn::Str(Cow::Owned(v.to_string()))))
    }

    fn visit_string<E: de::Error>(self, v: String) -> Result<Value, E> {
        Ok(Value(Edn::Str(Cow::Owned(v))))
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Value, A::Error> {
        let mut items: Vec<Edn<'static>> = Vec::with_capacity(seq.size_hint().unwrap_or(0));
        while let Some(element) = seq.next_element::<Value>()? {
            items.push(element.0);
        }
        Ok(Value(Edn::Vector(EdnSeq::from(items))))
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Value, A::Error> {
        let mut out: EdnMap<'static> = EdnMap::with_capacity(map.size_hint().unwrap_or(0));
        while let Some((k, v)) = map.next_entry::<Value, Value>()? {
            out.insert(k.0, v.0);
        }
        Ok(Value(Edn::Map(out)))
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<Value, D::Error> {
        deserializer.deserialize_any(self)
    }

    fn visit_enum<A: EnumAccess<'de>>(self, access: A) -> Result<Value, A::Error> {
        let (variant, payload): (String, A::Variant) = access.variant()?;
        match variant.as_str() {
            KEYWORD_TOKEN => {
                let s: String = payload.newtype_variant()?;
                Ok(Value(Edn::Keyword(Keyword::owned(s))))
            }
            SYMBOL_TOKEN => {
                let s: String = payload.newtype_variant()?;
                Ok(Value(Edn::Symbol(Symbol::owned(s))))
            }
            LIST_TOKEN => {
                let items: Vec<Value> = payload.newtype_variant()?;
                let owned: Vec<Edn<'static>> = items.into_iter().map(|v| v.0).collect();
                Ok(Value(Edn::List(EdnSeq::from(owned))))
            }
            SET_TOKEN => {
                let items: Vec<Value> = payload.newtype_variant()?;
                let mut set: EdnSet<'static> = EdnSet::new();
                for item in items {
                    set.insert(item.0);
                }
                Ok(Value(Edn::Set(set)))
            }
            TAGGED_TOKEN => {
                let (tag, inner): (String, Value) = payload.newtype_variant()?;
                Ok(Value(Edn::Tagged(Cow::Owned(tag), Box::new(inner.0))))
            }
            #[cfg(feature = "bignum")]
            BIGINT_TOKEN => {
                let s: String = payload.newtype_variant()?;
                let n = num_bigint::BigInt::from_str(&s).map_err(|e| {
                    <A::Error as de::Error>::custom(format!("invalid bigint {s:?}: {e}"))
                })?;
                Ok(Value(Edn::BigInt(n)))
            }
            #[cfg(feature = "bignum")]
            BIGDECIMAL_TOKEN => {
                let s: String = payload.newtype_variant()?;
                let d = bigdecimal::BigDecimal::from_str(&s).map_err(|e| {
                    <A::Error as de::Error>::custom(format!("invalid bigdecimal {s:?}: {e}"))
                })?;
                Ok(Value(Edn::BigDecimal(EdnBigDecimal::new(d))))
            }
            other => Err(<A::Error as de::Error>::custom(format!(
                "unexpected Value carrier variant: {other}"
            ))),
        }
    }
}

/// Payload variants carried from `EdnDeserializer` to `ValueVisitor` for
/// EDN-specific variants that do not map onto serde's primitive set.
#[cfg(feature = "serde")]
pub(crate) enum ValuePayload<'de> {
    Str(Cow<'de, str>),
    List(EdnSeq<'de>),
    Set(EdnSet<'de>),
    Tagged(Cow<'de, str>, Edn<'de>),
}

#[cfg(feature = "serde")]
pub(crate) struct ValueCarrier<'de> {
    variant: &'static str,
    payload: ValuePayload<'de>,
}

#[cfg(feature = "serde")]
impl<'de> ValueCarrier<'de> {
    pub(crate) fn new(variant: &'static str, payload: ValuePayload<'de>) -> Self {
        Self { variant, payload }
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserializer<'de> for ValueCarrier<'de> {
    type Error = EdnError;

    fn deserialize_any<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        visitor.visit_enum(ValueCarrierEnumAccess {
            variant: self.variant,
            payload: self.payload,
        })
    }

    serde::forward_to_deserialize_any! {
        bool i8 i16 i32 i64 i128 u8 u16 u32 u64 u128 f32 f64 char str string
        bytes byte_buf option unit unit_struct newtype_struct seq tuple
        tuple_struct map struct enum identifier ignored_any
    }
}

#[cfg(feature = "serde")]
struct ValueCarrierEnumAccess<'de> {
    variant: &'static str,
    payload: ValuePayload<'de>,
}

#[cfg(feature = "serde")]
impl<'de> EnumAccess<'de> for ValueCarrierEnumAccess<'de> {
    type Error = EdnError;
    type Variant = ValueCarrierVariantAccess<'de>;

    fn variant_seed<V: DeserializeSeed<'de>>(
        self,
        seed: V,
    ) -> Result<(V::Value, Self::Variant), Self::Error> {
        let variant_name =
            seed.deserialize(de::value::StrDeserializer::<EdnError>::new(self.variant))?;
        Ok((
            variant_name,
            ValueCarrierVariantAccess {
                payload: self.payload,
            },
        ))
    }
}

#[cfg(feature = "serde")]
struct ValueCarrierVariantAccess<'de> {
    payload: ValuePayload<'de>,
}

#[cfg(feature = "serde")]
impl<'de> VariantAccess<'de> for ValueCarrierVariantAccess<'de> {
    type Error = EdnError;

    fn unit_variant(self) -> Result<(), Self::Error> {
        Err(EdnError::Custom(
            "Value carrier does not support unit variants".into(),
        ))
    }

    fn newtype_variant_seed<T: DeserializeSeed<'de>>(
        self,
        seed: T,
    ) -> Result<T::Value, Self::Error> {
        match self.payload {
            ValuePayload::Str(s) => seed.deserialize(EdnDeserializer(Edn::Str(s))),
            ValuePayload::List(v) => seed.deserialize(EdnDeserializer(Edn::Vector(v))),
            ValuePayload::Set(s) => {
                let items: Vec<Edn<'de>> = s.into_iter().collect();
                seed.deserialize(EdnDeserializer(Edn::Vector(EdnSeq::from(items))))
            }
            ValuePayload::Tagged(tag, inner) => {
                let items = vec![Edn::Str(tag), inner];
                seed.deserialize(EdnDeserializer(Edn::Vector(EdnSeq::from(items))))
            }
        }
    }

    fn tuple_variant<V: Visitor<'de>>(
        self,
        _len: usize,
        _visitor: V,
    ) -> Result<V::Value, Self::Error> {
        Err(EdnError::Custom(
            "Value carrier does not support tuple variants".into(),
        ))
    }

    fn struct_variant<V: Visitor<'de>>(
        self,
        _fields: &'static [&'static str],
        _visitor: V,
    ) -> Result<V::Value, Self::Error> {
        Err(EdnError::Custom(
            "Value carrier does not support struct variants".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "serde")]
    use crate::{from_str, to_string};
    #[cfg(feature = "serde")]
    use rstest::rstest;

    #[test]
    fn test_value_parse() {
        assert_eq!(Value::parse("nil").unwrap().as_edn(), &Edn::Nil);
        assert_eq!(Value::parse("true").unwrap().as_edn(), &Edn::Bool(true));
        assert_eq!(Value::parse("7").unwrap().as_edn(), &Edn::Int(7));
    }

    #[test]
    fn test_value_parse_roundtrip() {
        let input = "{:name \"salt\" :atoms [:Na :Cl]}";
        let v = Value::parse(input).unwrap();
        let s = v.to_string();
        let v2 = Value::parse(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_parse_tagged() {
        let v = Value::parse("#inst \"2026-04-08\"").unwrap();
        let Edn::Tagged(tag, inner) = v.as_edn() else {
            panic!("expected tagged, got {:?}", v.as_edn());
        };
        assert_eq!(tag.as_ref(), "inst");
        assert_eq!(
            inner.as_ref(),
            &Edn::Str(std::borrow::Cow::Borrowed("2026-04-08"))
        );
    }

    #[test]
    fn test_value_from_edn() {
        let edn = crate::read_string("[1 2 3]").unwrap();
        let v = Value::from_edn(&edn).unwrap();
        assert_eq!(v.to_edn(), edn);
    }

    #[test]
    fn test_value_display() {
        let v = Value::parse(":foo").unwrap();
        assert_eq!(v.to_string(), ":foo");
    }

    #[test]
    fn test_value_from_str() {
        let v: Value = "(1 2 3)".parse().unwrap();
        assert!(matches!(v.as_edn(), Edn::List(_)));
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case("nil", Edn::Nil)]
    #[case("true", Edn::Bool(true))]
    #[case("false", Edn::Bool(false))]
    #[case("7", Edn::Int(7))]
    #[case("-3", Edn::Int(-3))]
    #[case("2.5", Edn::Float(2.5))]
    #[case(r#""hi""#, Edn::Str(Cow::Borrowed("hi")))]
    #[case(r#"\a"#, Edn::Char('a'))]
    fn test_value_roundtrip_primitives(#[case] input: &str, #[case] expected: Edn<'static>) {
        let v: Value = from_str(input).unwrap();
        assert_eq!(v.as_edn(), &expected);
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_keyword() {
        let v: Value = from_str(":name/foo").unwrap();
        assert!(matches!(v.as_edn(), Edn::Keyword(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, ":name/foo");
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[rstest]
    #[case("sym", "sym", |e: &Edn<'_>| matches!(e, Edn::Symbol(_)))]
    #[case("[1 2 3]", "[1 2 3]", |e: &Edn<'_>| matches!(e, Edn::Vector(_)))]
    fn test_value_roundtrip_shape(
        #[case] input: &str,
        #[case] expected: &str,
        #[case] kind: fn(&Edn<'_>) -> bool,
    ) {
        let v: Value = from_str(input).unwrap();
        assert!(kind(v.as_edn()));
        let s = to_string(&v).unwrap();
        assert_eq!(s, expected);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_list() {
        let v: Value = from_str("(1 2 3)").unwrap();
        assert!(matches!(v.as_edn(), Edn::List(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "(1 2 3)");
        let v2: Value = from_str(&s).unwrap();
        assert!(matches!(v2.as_edn(), Edn::List(_)));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_set() {
        let v: Value = from_str("#{1 2 3}").unwrap();
        assert!(matches!(v.as_edn(), Edn::Set(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert!(matches!(v2.as_edn(), Edn::Set(_)));
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_map() {
        let v: Value = from_str("{:a 1 :b 2}").unwrap();
        assert!(matches!(v.as_edn(), Edn::Map(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_tagged() {
        let v: Value = from_str("#inst \"2026-04-08\"").unwrap();
        let Edn::Tagged(tag, inner) = v.as_edn() else {
            panic!("expected tagged, got {:?}", v.as_edn());
        };
        assert_eq!(tag.as_ref(), "inst");
        assert_eq!(inner.as_ref(), &Edn::Str(Cow::Borrowed("2026-04-08")));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_roundtrip_nested() {
        let input = r#"{:name "salt" :atoms [:Na :Cl] :count 2}"#;
        let v: Value = from_str(input).unwrap();
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_serialize_in_struct() {
        use serde::{Deserialize, Serialize};
        #[derive(Serialize, Deserialize, PartialEq, Debug)]
        struct Wrapper {
            name: String,
            data: Value,
        }
        let input = r#"{:name "thing" :data (1 :two 3)}"#;
        let w: Wrapper = from_str(input).unwrap();
        assert_eq!(w.name, "thing");
        assert!(matches!(w.data.as_edn(), Edn::List(_)));
        let s = to_string(&w).unwrap();
        let w2: Wrapper = from_str(&s).unwrap();
        assert_eq!(w, w2);
    }

    #[cfg(all(feature = "serde", feature = "bignum"))]
    #[test]
    fn test_value_roundtrip_bigint() {
        let v: Value = from_str("123456789012345678901234567890N").unwrap();
        assert!(matches!(v.as_edn(), Edn::BigInt(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "123456789012345678901234567890N");
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(all(feature = "serde", feature = "bignum"))]
    #[test]
    fn test_value_roundtrip_bigdecimal() {
        let v: Value = from_str("3.14159265358979323846M").unwrap();
        assert!(matches!(v.as_edn(), Edn::BigDecimal(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_serialize_json() {
        let v = Value(Edn::Int(17));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "17");
        let back: Value = serde_json::from_str(&json).unwrap();
        assert_eq!(back, v);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_serialize_json_keyword() {
        let v = Value(Edn::Keyword(Keyword::owned("foo".to_string())));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "\"foo\"");
        let back: Value = serde_json::from_str(&json).unwrap();
        assert_eq!(back.as_edn(), &Edn::Str(Cow::Owned("foo".to_string())));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_serialize_json_list() {
        let v = Value(Edn::List(EdnSeq::from(vec![Edn::Int(1), Edn::Int(2)])));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "[1,2]");
        let back: Value = serde_json::from_str(&json).unwrap();
        assert!(matches!(back.as_edn(), Edn::Vector(_)));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_value_serialize_json_map() {
        let v: Value = from_str(r#"{"a" 1 "b" 2}"#).unwrap();
        let json = serde_json::to_string(&v).unwrap();
        let back: Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v, back);
    }
}
