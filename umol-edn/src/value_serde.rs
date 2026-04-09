//! Serde `Serialize`/`Deserialize` for [`Value`].
//!
//! The EDN data model has variants (keyword, symbol, list, set, tagged,
//! bigint, bigdecimal) that do not correspond to serde's primitive set.
//! `Value` round-trips them losslessly through the EDN serializer/
//! deserializer while still producing sensible fallbacks when the payload
//! travels over a foreign format such as JSON or YAML.
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

use std::borrow::Cow;
use std::fmt;

use serde::de::{
    self, DeserializeSeed, Deserializer, EnumAccess, MapAccess, SeqAccess, VariantAccess, Visitor,
};
use serde::ser::{
    Serialize, SerializeMap, SerializeSeq, SerializeTupleStruct, Serializer,
};
use serde::Deserialize;

use crate::collections::{EdnMap, EdnSeq, EdnSet};
use crate::de::EdnDeserializer;
use crate::edn::{Edn, Keyword, Symbol};
use crate::error::EdnError;
use crate::serde_tokens::{
    KEYWORD_TOKEN, LIST_TOKEN, SET_TOKEN, SYMBOL_TOKEN, TAGGED_TOKEN, VALUE_TOKEN,
};
#[cfg(feature = "bignum")]
use crate::serde_tokens::{BIGDECIMAL_TOKEN, BIGINT_TOKEN};
use crate::value::Value;

#[cfg(feature = "bignum")]
use crate::edn::EdnBigDecimal;
#[cfg(feature = "bignum")]
use std::str::FromStr;

// ---------------------------------------------------------------------------
// Serialize
// ---------------------------------------------------------------------------

impl Serialize for Value {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        EdnRef(&self.0).serialize(serializer)
    }
}

/// Borrowing helper that serializes any `&Edn` through the same
/// variant-preserving token dispatch used by `Value`.
struct EdnRef<'a>(&'a Edn<'a>);

impl<'a> Serialize for EdnRef<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self.0 {
            Edn::Nil => serializer.serialize_unit(),
            Edn::Bool(b) => serializer.serialize_bool(*b),
            Edn::Int(i) => serializer.serialize_i64(*i),
            Edn::Float(f) => serializer.serialize_f64(*f),
            Edn::Char(c) => serializer.serialize_char(*c),
            Edn::Str(s) => serializer.serialize_str(s),
            Edn::Keyword(k) => {
                serializer.serialize_newtype_struct(KEYWORD_TOKEN, k.as_str())
            }
            Edn::Symbol(s) => {
                serializer.serialize_newtype_struct(SYMBOL_TOKEN, s.as_str())
            }
            Edn::Vector(v) => {
                let mut seq = serializer.serialize_seq(Some(v.len()))?;
                for item in v.iter() {
                    seq.serialize_element(&EdnRef(item))?;
                }
                seq.end()
            }
            Edn::List(v) => {
                serializer.serialize_newtype_struct(LIST_TOKEN, &ListPayload(v))
            }
            Edn::Set(s) => {
                serializer.serialize_newtype_struct(SET_TOKEN, &SetPayload(s))
            }
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
            Edn::BigInt(n) => {
                serializer.serialize_newtype_struct(BIGINT_TOKEN, &n.to_string())
            }
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => {
                serializer
                    .serialize_newtype_struct(BIGDECIMAL_TOKEN, &d.as_inner().to_string())
            }
        }
    }
}

struct ListPayload<'a>(&'a EdnSeq<'a>);

impl<'a> Serialize for ListPayload<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for item in self.0.iter() {
            seq.serialize_element(&EdnRef(item))?;
        }
        seq.end()
    }
}

struct SetPayload<'a>(&'a EdnSet<'a>);

impl<'a> Serialize for SetPayload<'a> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.0.len()))?;
        for item in self.0.iter() {
            seq.serialize_element(&EdnRef(item))?;
        }
        seq.end()
    }
}

// ---------------------------------------------------------------------------
// Deserialize
// ---------------------------------------------------------------------------

impl<'de> Deserialize<'de> for Value {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_newtype_struct(VALUE_TOKEN, ValueVisitor)
    }
}

pub(crate) struct ValueVisitor;

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
        let mut out: EdnMap<'static> =
            EdnMap::with_capacity(map.size_hint().unwrap_or(0));
        while let Some((k, v)) = map.next_entry::<Value, Value>()? {
            out.insert(k.0, v.0);
        }
        Ok(Value(Edn::Map(out)))
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(
        self,
        deserializer: D,
    ) -> Result<Value, D::Error> {
        // Non-EDN fallback: recurse via `deserialize_any`. EDN itself never
        // reaches this arm — `EdnDeserializer::deserialize_newtype_struct`
        // routes EDN-specific variants through `visit_enum` on a carrier.
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
                let n = num_bigint::BigInt::from_str(&s)
                    .map_err(|e| <A::Error as de::Error>::custom(format!("invalid bigint {s:?}: {e}")))?;
                Ok(Value(Edn::BigInt(n)))
            }
            #[cfg(feature = "bignum")]
            BIGDECIMAL_TOKEN => {
                let s: String = payload.newtype_variant()?;
                let d = bigdecimal::BigDecimal::from_str(&s)
                    .map_err(|e| <A::Error as de::Error>::custom(format!("invalid bigdecimal {s:?}: {e}")))?;
                Ok(Value(Edn::BigDecimal(EdnBigDecimal::new(d))))
            }
            other => Err(<A::Error as de::Error>::custom(format!(
                "unexpected Value carrier variant: {other}"
            ))),
        }
    }
}

// ---------------------------------------------------------------------------
// ValueCarrier — the interop deserializer between `EdnDeserializer` and
// `ValueVisitor` for EDN-specific variants.
// ---------------------------------------------------------------------------

/// Payload variants carried from `EdnDeserializer` to `ValueVisitor`.
pub(crate) enum ValuePayload<'de> {
    /// String-backed payload (keyword name, symbol name, bigint/bigdec text).
    Str(Cow<'de, str>),
    /// Sequence payload (list items).
    List(EdnSeq<'de>),
    /// Set payload (converted to a seq of items).
    Set(EdnSet<'de>),
    /// Tagged literal — tag name and inner Edn.
    Tagged(Cow<'de, str>, Edn<'de>),
}

pub(crate) struct ValueCarrier<'de> {
    variant: &'static str,
    payload: ValuePayload<'de>,
}

impl<'de> ValueCarrier<'de> {
    pub(crate) fn new(variant: &'static str, payload: ValuePayload<'de>) -> Self {
        Self { variant, payload }
    }
}

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

struct ValueCarrierEnumAccess<'de> {
    variant: &'static str,
    payload: ValuePayload<'de>,
}

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

struct ValueCarrierVariantAccess<'de> {
    payload: ValuePayload<'de>,
}

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
                // Present as a 2-tuple (tag_str, inner_value).
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
    use crate::de::from_str;
    use crate::ser::to_string;
    use rstest::rstest;

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

    #[test]
    fn test_value_roundtrip_keyword() {
        let v: Value = from_str(":name/foo").unwrap();
        assert!(matches!(v.as_edn(), Edn::Keyword(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, ":name/foo");
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_roundtrip_symbol() {
        let v: Value = from_str("sym").unwrap();
        assert!(matches!(v.as_edn(), Edn::Symbol(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "sym");
    }

    #[test]
    fn test_value_roundtrip_vector() {
        let v: Value = from_str("[1 2 3]").unwrap();
        assert!(matches!(v.as_edn(), Edn::Vector(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "[1 2 3]");
    }

    #[test]
    fn test_value_roundtrip_list() {
        let v: Value = from_str("(1 2 3)").unwrap();
        assert!(matches!(v.as_edn(), Edn::List(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "(1 2 3)");
        let v2: Value = from_str(&s).unwrap();
        assert!(matches!(v2.as_edn(), Edn::List(_)));
    }

    #[test]
    fn test_value_roundtrip_set() {
        let v: Value = from_str("#{1 2 3}").unwrap();
        assert!(matches!(v.as_edn(), Edn::Set(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert!(matches!(v2.as_edn(), Edn::Set(_)));
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_roundtrip_map() {
        let v: Value = from_str("{:a 1 :b 2}").unwrap();
        assert!(matches!(v.as_edn(), Edn::Map(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_roundtrip_tagged() {
        let v: Value = from_str("#inst \"2026-04-08\"").unwrap();
        match v.as_edn() {
            Edn::Tagged(tag, inner) => {
                assert_eq!(tag.as_ref(), "inst");
                assert_eq!(inner.as_ref(), &Edn::Str(Cow::Borrowed("2026-04-08")));
            }
            other => panic!("expected tagged, got {other:?}"),
        }
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_roundtrip_nested() {
        let input = r#"{:name "salt" :atoms [:Na :Cl] :count 2}"#;
        let v: Value = from_str(input).unwrap();
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_in_struct_field() {
        #[derive(Serialize, Deserialize, PartialEq, Debug)]
        struct Wrapper {
            name: String,
            data: Value,
        }
        use serde::{Deserialize, Serialize};
        let input = r#"{:name "thing" :data (1 :two 3)}"#;
        let w: Wrapper = from_str(input).unwrap();
        assert_eq!(w.name, "thing");
        assert!(matches!(w.data.as_edn(), Edn::List(_)));
        let s = to_string(&w).unwrap();
        let w2: Wrapper = from_str(&s).unwrap();
        assert_eq!(w, w2);
    }

    #[cfg(feature = "bignum")]
    #[test]
    fn test_value_roundtrip_bigint() {
        let v: Value = from_str("123456789012345678901234567890N").unwrap();
        assert!(matches!(v.as_edn(), Edn::BigInt(_)));
        let s = to_string(&v).unwrap();
        assert_eq!(s, "123456789012345678901234567890N");
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[cfg(feature = "bignum")]
    #[test]
    fn test_value_roundtrip_bigdecimal() {
        let v: Value = from_str("3.14159265358979323846M").unwrap();
        assert!(matches!(v.as_edn(), Edn::BigDecimal(_)));
        let s = to_string(&v).unwrap();
        let v2: Value = from_str(&s).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_value_json_fallback_primitive() {
        let v = Value(Edn::Int(17));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "17");
        let back: Value = serde_json::from_str(&json).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn test_value_json_fallback_keyword_becomes_string() {
        let v = Value(Edn::Keyword(Keyword::owned("foo".to_string())));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "\"foo\"");
        let back: Value = serde_json::from_str(&json).unwrap();
        // Lossy: keyword → Str after JSON round-trip.
        assert_eq!(back.as_edn(), &Edn::Str(Cow::Owned("foo".to_string())));
    }

    #[test]
    fn test_value_json_fallback_list_becomes_vector() {
        let v = Value(Edn::List(EdnSeq::from(vec![
            Edn::Int(1),
            Edn::Int(2),
        ])));
        let json = serde_json::to_string(&v).unwrap();
        assert_eq!(json, "[1,2]");
        let back: Value = serde_json::from_str(&json).unwrap();
        assert!(matches!(back.as_edn(), Edn::Vector(_)));
    }

    #[test]
    fn test_value_json_fallback_map() {
        let v: Value = from_str(r#"{"a" 1 "b" 2}"#).unwrap();
        let json = serde_json::to_string(&v).unwrap();
        let back: Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v, back);
    }
}
