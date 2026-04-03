//! Serde `Deserializer` for `Edn` values.

use std::borrow::Cow;
use std::collections::hash_map::IntoIter as HashMapIntoIter;
use std::marker::PhantomData;

use serde::de::{
    self, Deserialize, DeserializeSeed, EnumAccess, MapAccess, SeqAccess, VariantAccess, Visitor,
};
use serde::Deserializer;

use crate::collections::{EdnMap, EdnSeq, EdnSeqIntoIter};
use crate::config::ParseConfig;
use crate::edn::Edn;
use crate::error::EdnError;
use crate::reader::Reader;
use crate::streaming::EdnStreamDeserializer;

/// Deserialize a Rust value from an EDN string.
///
/// Uses a streaming deserializer that parses directly into the target type
/// without building an intermediate `Edn` value tree.
pub fn from_str<'a, T: Deserialize<'a>>(s: &'a str) -> Result<T, EdnError> {
    from_str_with(s, crate::reader::default_config())
}

/// Deserialize a Rust value from an EDN string using a custom config.
pub fn from_str_with<'a, T: Deserialize<'a>>(s: &'a str, config: &ParseConfig) -> Result<T, EdnError> {
    let mut de = EdnStreamDeserializer::with_config(s, config);
    let val = T::deserialize(&mut de)?;
    de.expect_eof()?;
    Ok(val)
}

/// Deserialize a Rust value from a pre-parsed `Edn` value.
pub fn from_value<'a, T: Deserialize<'a>>(val: Edn<'a>) -> Result<T, EdnError> {
    T::deserialize(EdnDeserializer(val)).map_err(Into::into)
}

/// Streaming deserializer over multiple EDN values in a string.
pub struct StreamDeserializer<'a, T> {
    reader: Reader<'a>,
    _marker: PhantomData<T>,
}

impl<'a, T> StreamDeserializer<'a, T> {
    pub fn new(input: &'a str) -> Self {
        Self {
            reader: Reader::new(input),
            _marker: PhantomData,
        }
    }

    pub fn with_config(input: &'a str, config: ParseConfig) -> Self {
        Self {
            reader: Reader::with_config(input, config),
            _marker: PhantomData,
        }
    }
}

impl<'a, T: Deserialize<'a>> Iterator for StreamDeserializer<'a, T> {
    type Item = Result<T, EdnError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next().map(|result| {
            result.and_then(|edn| T::deserialize(EdnDeserializer(edn)).map_err(Into::into))
        })
    }
}

/// Deserializer wrapping an `Edn` value.
pub struct EdnDeserializer<'de>(pub Edn<'de>);

impl<'de> de::Deserializer<'de> for EdnDeserializer<'de> {
    type Error = EdnError;

    fn deserialize_any<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Nil => visitor.visit_unit(),
            Edn::Bool(b) => visitor.visit_bool(b),
            Edn::Int(i) => visitor.visit_i64(i),
            #[cfg(feature = "bignum")]
            Edn::BigInt(n) => {
                let mut s = n.to_string();
                s.push('N');
                visitor.visit_string(s)
            }
            Edn::Float(f) => visitor.visit_f64(f),
            #[cfg(feature = "bignum")]
            Edn::BigDecimal(d) => {
                let mut s = d.to_string();
                s.push('M');
                visitor.visit_string(s)
            }
            Edn::Char(c) => visitor.visit_char(c),
            Edn::Str(s) => match s {
                Cow::Borrowed(b) => visitor.visit_borrowed_str(b),
                Cow::Owned(o) => visitor.visit_string(o),
            },
            Edn::Keyword(k) => visitor.visit_string(k.as_str().to_string()),
            Edn::Symbol(s) => visitor.visit_string(s.as_str().to_string()),
            Edn::List(v) | Edn::Vector(v) => visitor.visit_seq(EdnSeqAccess::new(v)),
            Edn::Map(m) => visitor.visit_map(EdnMapAccess::new(m)),
            Edn::Set(s) => {
                let v: EdnSeq<'de> = s.into_iter().collect();
                visitor.visit_seq(EdnSeqAccess::new(v))
            }
            Edn::Tagged(_tag, inner) => EdnDeserializer(*inner).deserialize_any(visitor),
        }
    }

    fn deserialize_option<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match &self.0 {
            Edn::Nil => visitor.visit_none(),
            _ => visitor.visit_some(self),
        }
    }

    fn deserialize_newtype_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        visitor.visit_newtype_struct(self)
    }

    fn deserialize_enum<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _variants: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Keyword(k) => {
                visitor.visit_enum(de::value::StrDeserializer::<Self::Error>::new(k.as_str()))
            }
            Edn::Symbol(s) => {
                visitor.visit_enum(de::value::StrDeserializer::<Self::Error>::new(s.as_str()))
            }
            Edn::Str(ref s) => {
                visitor.visit_enum(de::value::StrDeserializer::<Self::Error>::new(s))
            }
            Edn::Tagged(tag, inner) => {
                visitor.visit_enum(EdnTaggedEnumAccess { tag, inner: *inner })
            }
            other => EdnDeserializer(other).deserialize_any(visitor),
        }
    }

    fn deserialize_map<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Map(m) => visitor.visit_map(EdnMapAccess::new(m)),
            Edn::Nil => visitor.visit_map(EdnMapAccess::new(EdnMap::new())),
            other => Err(EdnError::Custom(format!("expected map, got {other:?}"))),
        }
    }

    fn deserialize_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Map(m) => visitor.visit_map(EdnStructMapAccess::new(m)),
            Edn::Nil => visitor.visit_map(EdnStructMapAccess::new(EdnMap::new())),
            other => Err(EdnError::Custom(format!("expected map, got {other:?}"))),
        }
    }

    fn deserialize_unit_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_unit(visitor)
    }

    fn deserialize_tuple<V: Visitor<'de>>(
        self,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_seq(visitor)
    }

    fn deserialize_tuple_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        self.deserialize_seq(visitor)
    }

    fn deserialize_i8<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = i8::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for i8")))?;
                visitor.visit_i8(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_i16<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = i16::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for i16")))?;
                visitor.visit_i16(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_i32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = i32::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for i32")))?;
                visitor.visit_i32(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u8<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = u8::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for u8")))?;
                visitor.visit_u8(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u16<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = u16::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for u16")))?;
                visitor.visit_u16(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = u32::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for u32")))?;
                visitor.visit_u32(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u64<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Int(i) => {
                let v = u64::try_from(i)
                    .map_err(|_| EdnError::Custom(format!("{i} out of range for u64")))?;
                visitor.visit_u64(v)
            }
            other => Err(EdnError::Custom(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_f32<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Float(f) => visitor.visit_f32(f as f32),
            Edn::Int(i) => visitor.visit_f32(i as f32),
            other => Err(EdnError::Custom(format!("expected number, got {other:?}"))),
        }
    }

    fn deserialize_f64<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Float(f) => visitor.visit_f64(f),
            Edn::Int(i) => visitor.visit_f64(i as f64),
            other => Err(EdnError::Custom(format!("expected number, got {other:?}"))),
        }
    }

    fn deserialize_string<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.deserialize_str(visitor)
    }

    fn deserialize_identifier<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        match self.0 {
            Edn::Keyword(k) => visitor.visit_string(k.as_str().to_string()),
            Edn::Symbol(s) => visitor.visit_string(s.as_str().to_string()),
            Edn::Str(s) => match s {
                Cow::Borrowed(b) => visitor.visit_borrowed_str(b),
                Cow::Owned(o) => visitor.visit_string(o),
            },
            // Non-string keys: convert to string so serde handles them as
            // unknown field names, consistent with the streaming path.
            other => visitor.visit_string(other.to_string()),
        }
    }

    fn deserialize_bytes<V: Visitor<'de>>(self, _visitor: V) -> Result<V::Value, Self::Error> {
        Err(EdnError::Custom("bytes not supported".to_string()))
    }

    fn deserialize_byte_buf<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, Self::Error> {
        self.deserialize_bytes(visitor)
    }

    serde::forward_to_deserialize_any! {
        bool i64 char str unit ignored_any seq
    }
}

// --- SeqAccess ---

struct EdnSeqAccess<'de> {
    iter: EdnSeqIntoIter<'de>,
}

impl<'de> EdnSeqAccess<'de> {
    fn new(v: EdnSeq<'de>) -> Self {
        Self {
            iter: v.into_iter(),
        }
    }
}

impl<'de> SeqAccess<'de> for EdnSeqAccess<'de> {
    type Error = EdnError;

    fn next_element_seed<T: DeserializeSeed<'de>>(
        &mut self,
        seed: T,
    ) -> Result<Option<T::Value>, Self::Error> {
        match self.iter.next() {
            Some(edn) => seed.deserialize(EdnDeserializer(edn)).map(Some),
            None => Ok(None),
        }
    }
}

// --- MapAccess ---

struct EdnMapAccess<'de> {
    iter: HashMapIntoIter<Edn<'de>, Edn<'de>>,
    pending_value: Option<Edn<'de>>,
}

impl<'de> EdnMapAccess<'de> {
    fn new(map: EdnMap<'de>) -> Self {
        Self {
            iter: map.into_iter(),
            pending_value: None,
        }
    }
}

impl<'de> MapAccess<'de> for EdnMapAccess<'de> {
    type Error = EdnError;

    fn next_key_seed<K: DeserializeSeed<'de>>(
        &mut self,
        seed: K,
    ) -> Result<Option<K::Value>, Self::Error> {
        match self.iter.next() {
            Some((k, v)) => {
                self.pending_value = Some(v);
                seed.deserialize(EdnDeserializer(k)).map(Some)
            }
            None => Ok(None),
        }
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(
        &mut self,
        seed: V,
    ) -> Result<V::Value, Self::Error> {
        let v = self
            .pending_value
            .take()
            .ok_or_else(|| EdnError::Custom("next_value_seed called without preceding next_key_seed".into()))?;
        seed.deserialize(EdnDeserializer(v))
    }
}

// --- StructMap (filters to string-like keys) ---

struct EdnStructMapAccess<'de> {
    iter: HashMapIntoIter<Edn<'de>, Edn<'de>>,
    pending_value: Option<Edn<'de>>,
}

impl<'de> EdnStructMapAccess<'de> {
    fn new(map: EdnMap<'de>) -> Self {
        Self {
            iter: map.into_iter(),
            pending_value: None,
        }
    }
}

impl<'de> de::MapAccess<'de> for EdnStructMapAccess<'de> {
    type Error = EdnError;

    fn next_key_seed<K: DeserializeSeed<'de>>(
        &mut self,
        seed: K,
    ) -> Result<Option<K::Value>, Self::Error> {
        match self.iter.next() {
            Some((k, v)) => {
                self.pending_value = Some(v);
                seed.deserialize(EdnDeserializer(k)).map(Some)
            }
            None => Ok(None),
        }
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(
        &mut self,
        seed: V,
    ) -> Result<V::Value, Self::Error> {
        let v = self
            .pending_value
            .take()
            .ok_or_else(|| EdnError::Custom("next_value_seed called without preceding next_key_seed".into()))?;
        seed.deserialize(EdnDeserializer(v))
    }
}

// -- Tagged enum access (for #Variant value round-tripping) ------------------

struct EdnTaggedEnumAccess<'de> {
    tag: String,
    inner: Edn<'de>,
}

impl<'de> EnumAccess<'de> for EdnTaggedEnumAccess<'de> {
    type Error = EdnError;
    type Variant = EdnTaggedVariantAccess<'de>;

    fn variant_seed<V: DeserializeSeed<'de>>(
        self,
        seed: V,
    ) -> Result<(V::Value, Self::Variant), Self::Error> {
        let variant = seed.deserialize(de::value::StrDeserializer::<EdnError>::new(&self.tag))?;
        Ok((variant, EdnTaggedVariantAccess(self.inner)))
    }
}

struct EdnTaggedVariantAccess<'de>(Edn<'de>);

impl<'de> VariantAccess<'de> for EdnTaggedVariantAccess<'de> {
    type Error = EdnError;

    fn unit_variant(self) -> Result<(), Self::Error> {
        match self.0 {
            Edn::Nil => Ok(()),
            other => Err(EdnError::Custom(
                format!("unit variant expects nil payload, got {other:?}"),
            )),
        }
    }

    fn newtype_variant_seed<T: DeserializeSeed<'de>>(
        self,
        seed: T,
    ) -> Result<T::Value, Self::Error> {
        seed.deserialize(EdnDeserializer(self.0))
    }

    fn tuple_variant<V: Visitor<'de>>(
        self,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        EdnDeserializer(self.0).deserialize_seq(visitor)
    }

    fn struct_variant<V: Visitor<'de>>(
        self,
        _fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error> {
        EdnDeserializer(self.0).deserialize_map(visitor)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use rstest::rstest;
    use serde::{Deserialize, Serialize};

    use super::*;
    use crate::{from_str, read_string};

    #[rstest]
    #[case("12", 12i64)]
    #[case("-1", -1i64)]
    #[case("0", 0i64)]
    fn test_deserialize_i64(#[case] input: &str, #[case] expected: i64) {
        assert_eq!(from_str::<i64>(input).unwrap(), expected);
    }

    #[rstest]
    #[case("12", 12u32)]
    #[case("0", 0u32)]
    fn test_deserialize_u32(#[case] input: &str, #[case] expected: u32) {
        assert_eq!(from_str::<u32>(input).unwrap(), expected);
    }

    #[rstest]
    #[case("256")]
    #[case("-1")]
    fn test_deserialize_u8_error(#[case] input: &str) {
        assert!(from_str::<u8>(input).is_err());
    }

    #[rstest]
    #[case("3.14", 3.14f64)]
    #[case("12", 12.0f64)]
    #[case("-0.5", -0.5f64)]
    fn test_deserialize_f64(#[case] input: &str, #[case] expected: f64) {
        assert!((from_str::<f64>(input).unwrap() - expected).abs() < 1e-10);
    }

    #[rstest]
    #[case("true", true)]
    #[case("false", false)]
    fn test_deserialize_bool(#[case] input: &str, #[case] expected: bool) {
        assert_eq!(from_str::<bool>(input).unwrap(), expected);
    }

    #[test]
    fn test_deserialize_char() {
        assert_eq!(from_str::<char>(r"\a").unwrap(), 'a');
    }

    #[rstest]
    #[case(r#""hello""#, "hello")]
    #[case(r#""with \"quotes\"""#, "with \"quotes\"")]
    #[case(r#""line\nbreak""#, "line\nbreak")]
    fn test_deserialize_string(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(from_str::<String>(input).unwrap(), expected);
    }

    #[test]
    fn test_deserialize_keyword_as_string() {
        assert_eq!(from_str::<String>(":foo").unwrap(), "foo");
    }

    #[rstest]
    #[case("12", Some(12i64))]
    #[case("nil", None)]
    fn test_deserialize_option(#[case] input: &str, #[case] expected: Option<i64>) {
        assert_eq!(from_str::<Option<i64>>(input).unwrap(), expected);
    }

    #[rstest]
    #[case("[1 2 3]", vec![1, 2, 3])]
    #[case("(1 2 3)", vec![1, 2, 3])]
    #[case("[]", vec![])]
    fn test_deserialize_vec(#[case] input: &str, #[case] expected: Vec<i64>) {
        assert_eq!(from_str::<Vec<i64>>(input).unwrap(), expected);
    }

    #[derive(Debug, Deserialize, PartialEq)]
    struct Point {
        x: f64,
        y: f64,
    }

    #[test]
    fn test_deserialize_struct() {
        let input = "{:x 1.0 :y 2.0}";
        assert_eq!(from_str::<Point>(input).unwrap(), Point { x: 1.0, y: 2.0 });
    }

    #[test]
    fn test_deserialize_struct_from_nil() {
        assert!(from_str::<Point>("nil").is_err());
    }

    #[derive(Debug, Deserialize, PartialEq)]
    struct OptionalFields {
        name: String,
        age: Option<i64>,
    }

    #[test]
    fn test_deserialize_struct_optional_field() {
        let input = r#"{:name "Alice" :age 30}"#;
        assert_eq!(
            from_str::<OptionalFields>(input).unwrap(),
            OptionalFields {
                name: "Alice".into(),
                age: Some(30)
            },
        );
    }

    #[test]
    fn test_deserialize_hashmap() {
        let input = r#"{"a" 1 "b" 2}"#;
        let m: HashMap<String, i64> = from_str(input).unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m["a"], 1);
        assert_eq!(m["b"], 2);
    }

    #[derive(Debug, Deserialize, PartialEq)]
    #[serde(rename_all = "lowercase")]
    enum Color {
        Red,
        Green,
        Blue,
    }

    #[rstest]
    #[case(":red", Color::Red)]
    #[case(":green", Color::Green)]
    #[case(":blue", Color::Blue)]
    fn test_deserialize_enum(#[case] input: &str, #[case] expected: Color) {
        assert_eq!(from_str::<Color>(input).unwrap(), expected);
    }

    #[test]
    fn test_deserialize_tuple() {
        let input = "[1 2 3]";
        assert_eq!(from_str::<(i64, i64, i64)>(input).unwrap(), (1, 2, 3));
    }

    #[derive(Debug, Deserialize, PartialEq)]
    struct Nested {
        point: Point,
        label: String,
    }

    #[test]
    fn test_deserialize_nested() {
        let input = r#"{:point {:x 3.0 :y 4.0} :label "origin"}"#;
        assert_eq!(
            from_str::<Nested>(input).unwrap(),
            Nested {
                point: Point { x: 3.0, y: 4.0 },
                label: "origin".into(),
            },
        );
    }

    #[test]
    fn test_deserialize_tagged_unwraps() {
        // Unknown tags produce Tagged(...) which is stripped during deserialization.
        let val = read_string("#my/custom [1 2 3]").unwrap();
        let v: Vec<i64> = Vec::deserialize(EdnDeserializer(val)).unwrap();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_deserialize_unit() {
        assert_eq!(from_str::<()>("nil").unwrap(), ());
    }

    #[test]
    fn test_deserialize_newtype_struct() {
        #[derive(Debug, Deserialize, PartialEq)]
        struct Wrapper(i64);
        assert_eq!(from_str::<Wrapper>("12").unwrap(), Wrapper(12));
    }

    #[test]
    fn test_deserialize_struct_non_string_key() {
        // Non-string keys are stringified and treated as unknown fields by serde.
        let expected = Point { x: 1.0, y: 2.0 };

        let val = read_string("{:x 1.0 12 99 :y 2.0}").unwrap();
        assert_eq!(Point::deserialize(EdnDeserializer(val)).unwrap(), expected);
        assert_eq!(from_str::<Point>("{:x 1.0 12 99 :y 2.0}").unwrap(), expected);
    }

    #[test]
    fn test_deserialize_struct_non_string_key_denied() {
        #[derive(Debug, Deserialize, PartialEq)]
        #[serde(deny_unknown_fields)]
        struct Strict {
            x: f64,
            y: f64,
        }

        let val = read_string("{:x 1.0 12 99 :y 2.0}").unwrap();
        assert!(Strict::deserialize(EdnDeserializer(val)).is_err());
        assert!(from_str::<Strict>("{:x 1.0 12 99 :y 2.0}").is_err());
    }

    #[test]
    fn test_deserialize_f32() {
        assert!((from_str::<f32>("3.14").unwrap() - 3.14f32).abs() < 1e-5);
        assert!((from_str::<f32>("12").unwrap() - 12.0f32).abs() < 1e-5);
    }

    #[rstest]
    #[case("1 2 3", vec![1, 2, 3])]
    #[case("  1  2  3  ", vec![1, 2, 3])]
    #[case("", vec![])]
    #[case("7", vec![7])]
    fn test_stream_deserializer_i64(#[case] input: &str, #[case] expected: Vec<i64>) {
        let results: Result<Vec<i64>, _> = StreamDeserializer::new(input).collect();
        assert_eq!(results.unwrap(), expected);
    }

    #[test]
    fn test_stream_deserializer_mixed_types() {
        #[derive(Debug, Deserialize, PartialEq)]
        struct Record {
            name: String,
            value: i64,
        }
        let input = r#"{:name "a" :value 1} {:name "b" :value 2}"#;
        let results: Vec<Record> = StreamDeserializer::new(input)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(
            results,
            vec![
                Record {
                    name: "a".into(),
                    value: 1
                },
                Record {
                    name: "b".into(),
                    value: 2
                },
            ]
        );
    }

    #[test]
    fn test_stream_deserializer_error_propagation() {
        let input = "1 2 [invalid";
        let results: Vec<Result<i64, _>> = StreamDeserializer::new(input).collect();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].as_ref().unwrap(), &1);
        assert_eq!(results[1].as_ref().unwrap(), &2);
        assert!(results[2].is_err());
    }

    #[test]
    fn test_stream_deserializer_stops_after_error() {
        let input = "1 [invalid 3";
        let results: Vec<Result<i64, _>> = StreamDeserializer::new(input).collect();
        // Reader sets remaining to "" on error, so only 2 items
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].as_ref().unwrap(), &1);
        assert!(results[1].is_err());
    }

    // -- Enum round-tripping via tags ------------------------------------------

    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    enum Shape {
        Circle(f64),
        Rect(f64, f64),
        Named { width: f64, height: f64 },
    }

    #[rstest]
    #[case(Shape::Circle(3.14), "#Circle 3.14")]
    #[case(Shape::Rect(1.0, 2.0), "#Rect [1.0 2.0]")]
    #[case(Shape::Named { width: 5.0, height: 10.0 }, "#Named {:width 5.0 :height 10.0}")]
    fn test_enum_tagged_roundtrip(#[case] value: Shape, #[case] expected_edn: &str) {
        let serialized = crate::to_string(&value).unwrap();
        assert_eq!(serialized, expected_edn);
        let deserialized: Shape = from_str(&serialized).unwrap();
        assert_eq!(deserialized, value);
    }

    #[derive(Debug, Deserialize, PartialEq)]
    enum UnitEnum {
        Red,
    }

    #[test]
    fn test_deserialize_tagged_unit_variant() {
        let val: UnitEnum = from_str("#Red nil").unwrap();
        assert_eq!(val, UnitEnum::Red);
    }

    #[rstest]
    #[case("#Red 1")]
    #[case("#Red {:a 1}")]
    #[case("#Red [1 2]")]
    #[case("#Red \"hello\"")]
    fn test_deserialize_tagged_unit_variant_error(#[case] input: &str) {
        assert!(from_str::<UnitEnum>(input).is_err());
    }

}
