//! EDN serde wrappers.
//!
//! Custom `Deserializer` and `Serializer` for `clojure_reader::edn::Edn` values.
//! Atom and bond specs are plain strings — no EDN tagged literals.

use std::collections::BTreeMap;
use std::fmt;

use clojure_reader::edn::{self, Edn};
use clojure_reader::error::Error as ClojureReaderError;
use serde::de::{self, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use serde::{ser, Deserialize, Serialize};

use super::error::ParseError;

/// Thin wrapper around `clojure_reader::error::Error` that implements serde error traits.
#[derive(Debug)]
pub struct EdnError(pub ClojureReaderError);

impl EdnError {
    pub fn msg(s: impl fmt::Display) -> Self {
        EdnError(<ClojureReaderError as de::Error>::custom(s))
    }

    /// Convert to `ParseError`, extracting structured info from `Code` variants
    /// instead of stringifying the whole error.
    pub fn into_parse_error(self) -> ParseError {
        map_edn_error(self.0)
    }
}

impl fmt::Display for EdnError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for EdnError {}

impl de::Error for EdnError {
    fn custom<T: fmt::Display>(msg: T) -> Self {
        Self::msg(msg)
    }
}

impl ser::Error for EdnError {
    fn custom<T: fmt::Display>(msg: T) -> Self {
        Self::msg(msg)
    }
}

impl From<ClojureReaderError> for EdnError {
    fn from(e: ClojureReaderError) -> Self {
        EdnError(e)
    }
}

/// Map `ClojureReaderError` to a `ParseError` variant, extracting
/// structured info where possible and falling back to `EdnParse` for the rest.
pub fn map_edn_error(e: ClojureReaderError) -> ParseError {
    use clojure_reader::error::Code;
    match e.code {
        Code::Serde(msg) => ParseError::EdnParse(msg),
        Code::UnexpectedEOF => ParseError::Incomplete,
        code => ParseError::EdnParse(format!("{code:?}")),
    }
}

// Deserializer

/// EDN deserializer that maps `Edn` values to serde visitors.
pub struct EdnDeserializer<'de>(pub Edn<'de>);

impl<'de> de::Deserializer<'de> for EdnDeserializer<'de> {
    type Error = EdnError;

    fn deserialize_any<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Key(k) => visitor.visit_borrowed_str(k),
            Edn::Str(s) | Edn::Symbol(s) => visitor.visit_borrowed_str(s),
            Edn::Int(i) => visitor.visit_i64(i),
            Edn::Char(c) => visitor.visit_char(c),
            Edn::Bool(b) => visitor.visit_bool(b),
            Edn::Nil => visitor.visit_unit(),
            Edn::Vector(list) | Edn::List(list) => visitor.visit_seq(EdnSeq::new(list)),
            Edn::Map(map) => {
                if map.is_empty() {
                    visitor.visit_unit()
                } else {
                    visitor.visit_map(EdnMap::new(map))
                }
            }
            Edn::Set(set) => {
                let v: Vec<Edn<'de>> = set.into_iter().collect();
                visitor.visit_seq(EdnSeq::new(v))
            }
            other => Err(EdnError::msg(format!("unsupported EDN value: {other:?}"))),
        }
    }

    fn deserialize_option<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match &self.0 {
            Edn::Nil => visitor.visit_none(),
            _ => visitor.visit_some(self),
        }
    }

    fn deserialize_newtype_struct<V>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        visitor.visit_newtype_struct(self)
    }

    fn deserialize_enum<V>(
        self,
        _name: &'static str,
        _variants: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        // For string-based enums (like BondSpec keywords), try visit_enum with string
        match self.0 {
            Edn::Key(k) => {
                visitor.visit_enum(de::value::BorrowedStrDeserializer::<Self::Error>::new(k))
            }
            Edn::Str(s) | Edn::Symbol(s) => {
                visitor.visit_enum(de::value::BorrowedStrDeserializer::<Self::Error>::new(s))
            }
            other => EdnDeserializer(other).deserialize_any(visitor),
        }
    }

    fn deserialize_struct<V>(
        self,
        _name: &'static str,
        _fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_map(visitor)
    }

    fn deserialize_unit_struct<V>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_unit(visitor)
    }

    fn deserialize_tuple<V>(self, _len: usize, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_seq(visitor)
    }

    fn deserialize_tuple_struct<V>(
        self,
        _name: &'static str,
        _len: usize,
        visitor: V,
    ) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_seq(visitor)
    }

    fn deserialize_i8<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = i8::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for i8")))?;
                visitor.visit_i8(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_i16<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = i16::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for i16")))?;
                visitor.visit_i16(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_i32<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = i32::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for i32")))?;
                visitor.visit_i32(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u8<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = u8::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for u8")))?;
                visitor.visit_u8(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u16<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = u16::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for u16")))?;
                visitor.visit_u16(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u32<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = u32::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for u32")))?;
                visitor.visit_u32(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_u64<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        match self.0 {
            Edn::Int(i) => {
                let v = u64::try_from(i)
                    .map_err(|_| EdnError::msg(format!("integer {i} out of range for u64")))?;
                visitor.visit_u64(v)
            }
            other => Err(EdnError::msg(format!("expected integer, got {other:?}"))),
        }
    }

    fn deserialize_string<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_str(visitor)
    }

    fn deserialize_identifier<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_str(visitor)
    }

    fn deserialize_bytes<V>(self, _visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        Err(EdnError::msg("deserialize_bytes not supported"))
    }

    fn deserialize_byte_buf<V>(self, visitor: V) -> Result<V::Value, Self::Error>
    where
        V: Visitor<'de>,
    {
        self.deserialize_bytes(visitor)
    }

    serde::forward_to_deserialize_any! {
        bool i64 f32 f64 char str unit map ignored_any seq
    }
}

// -- SeqAccess wrapping elements in EdnDeserializer -------------------------

struct EdnSeq<'de> {
    iter: std::vec::IntoIter<Edn<'de>>,
}

impl<'de> EdnSeq<'de> {
    fn new(v: Vec<Edn<'de>>) -> Self {
        Self {
            iter: v.into_iter(),
        }
    }
}

impl<'de> SeqAccess<'de> for EdnSeq<'de> {
    type Error = EdnError;

    fn next_element_seed<T>(&mut self, seed: T) -> Result<Option<T::Value>, Self::Error>
    where
        T: DeserializeSeed<'de>,
    {
        match self.iter.next() {
            Some(edn) => seed.deserialize(EdnDeserializer(edn)).map(Some),
            None => Ok(None),
        }
    }
}

// -- MapAccess wrapping keys/values in EdnDeserializer ----------------------

struct EdnMap<'de> {
    iter: std::collections::btree_map::IntoIter<Edn<'de>, Edn<'de>>,
    pending_value: Option<Edn<'de>>,
}

impl<'de> EdnMap<'de> {
    fn new(map: BTreeMap<Edn<'de>, Edn<'de>>) -> Self {
        Self {
            iter: map.into_iter(),
            pending_value: None,
        }
    }
}

impl<'de> MapAccess<'de> for EdnMap<'de> {
    type Error = EdnError;

    fn next_key_seed<K>(&mut self, seed: K) -> Result<Option<K::Value>, Self::Error>
    where
        K: DeserializeSeed<'de>,
    {
        loop {
            match self.iter.next() {
                Some((k, v)) => match &k {
                    Edn::Key(_) | Edn::Symbol(_) | Edn::Str(_) => {
                        self.pending_value = Some(v);
                        return seed.deserialize(EdnDeserializer(k)).map(Some);
                    }
                    _ => continue, // skip non-string keys
                },
                None => return Ok(None),
            }
        }
    }

    fn next_value_seed<V>(&mut self, seed: V) -> Result<V::Value, Self::Error>
    where
        V: DeserializeSeed<'de>,
    {
        let v = self
            .pending_value
            .take()
            .expect("next_value_seed called without preceding next_key_seed");
        seed.deserialize(EdnDeserializer(v))
    }
}

/// EDN serializer.
pub struct EdnSerializer {
    output: String,
}

impl EdnSerializer {
    fn new() -> Self {
        Self {
            output: String::new(),
        }
    }
}

impl ser::Serializer for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;
    type SerializeSeq = Self;
    type SerializeTuple = Self;
    type SerializeTupleStruct = Self;
    type SerializeTupleVariant = Self;
    type SerializeMap = Self;
    type SerializeStruct = Self;
    type SerializeStructVariant = Self;

    fn serialize_bool(self, v: bool) -> Result<(), Self::Error> {
        self.output += if v { "true" } else { "false" };
        Ok(())
    }

    fn serialize_i8(self, v: i8) -> Result<(), Self::Error> {
        self.serialize_i64(i64::from(v))
    }

    fn serialize_i16(self, v: i16) -> Result<(), Self::Error> {
        self.serialize_i64(i64::from(v))
    }

    fn serialize_i32(self, v: i32) -> Result<(), Self::Error> {
        self.serialize_i64(i64::from(v))
    }

    fn serialize_i64(self, v: i64) -> Result<(), Self::Error> {
        self.output += &v.to_string();
        Ok(())
    }

    fn serialize_u8(self, v: u8) -> Result<(), Self::Error> {
        self.serialize_u64(u64::from(v))
    }

    fn serialize_u16(self, v: u16) -> Result<(), Self::Error> {
        self.serialize_u64(u64::from(v))
    }

    fn serialize_u32(self, v: u32) -> Result<(), Self::Error> {
        self.serialize_u64(u64::from(v))
    }

    fn serialize_u64(self, v: u64) -> Result<(), Self::Error> {
        self.output += &v.to_string();
        Ok(())
    }

    fn serialize_f32(self, v: f32) -> Result<(), Self::Error> {
        self.serialize_f64(f64::from(v))
    }

    fn serialize_f64(self, v: f64) -> Result<(), Self::Error> {
        self.output += &v.to_string();
        Ok(())
    }

    fn serialize_char(self, v: char) -> Result<(), Self::Error> {
        self.output.push('\\');
        self.output.push(v);
        Ok(())
    }

    fn serialize_str(self, v: &str) -> Result<(), Self::Error> {
        self.output += "\"";
        self.output += v;
        self.output += "\"";
        Ok(())
    }

    fn serialize_bytes(self, _v: &[u8]) -> Result<(), Self::Error> {
        Err(ser::Error::custom("serialize_bytes not supported"))
    }

    fn serialize_none(self) -> Result<(), Self::Error> {
        self.output += "nil";
        Ok(())
    }

    fn serialize_some<T: ?Sized + Serialize>(self, value: &T) -> Result<(), Self::Error> {
        value.serialize(self)
    }

    fn serialize_unit(self) -> Result<(), Self::Error> {
        self.output += "nil";
        Ok(())
    }

    fn serialize_unit_struct(self, _name: &'static str) -> Result<(), Self::Error> {
        self.serialize_unit()
    }

    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
    ) -> Result<(), Self::Error> {
        self.output += ":";
        self.output += variant;
        Ok(())
    }

    fn serialize_newtype_struct<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        value.serialize(self)
    }

    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        // For enum variants with data, emit as tagged
        self.output += "#";
        self.output += variant;
        self.output += " ";
        value.serialize(self)
    }

    fn serialize_seq(self, _len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        self.output += "[";
        Ok(self)
    }

    fn serialize_tuple(self, _len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        self.output += "[";
        Ok(self)
    }

    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleStruct, Self::Error> {
        self.serialize_tuple(len)
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleVariant, Self::Error> {
        self.output += "#";
        self.output += variant;
        self.output += " [";
        Ok(self)
    }

    fn serialize_map(self, _len: Option<usize>) -> Result<Self::SerializeMap, Self::Error> {
        self.output += "{";
        Ok(self)
    }

    fn serialize_struct(
        self,
        _name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStruct, Self::Error> {
        self.serialize_map(Some(len))
    }

    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeStructVariant, Self::Error> {
        self.output += "#";
        self.output += variant;
        self.output += " {";
        Ok(self)
    }
}

impl ser::SerializeSeq for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output += " ";
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "]";
        Ok(())
    }
}

impl ser::SerializeTuple for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output += " ";
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "]";
        Ok(())
    }
}

impl ser::SerializeTupleStruct for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output += " ";
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "]";
        Ok(())
    }
}

impl ser::SerializeTupleVariant for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output += " ";
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "]";
        Ok(())
    }
}

impl ser::SerializeMap for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_key<T: ?Sized + Serialize>(&mut self, key: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output += " ";
        }
        key.serialize(&mut **self)
    }

    fn serialize_value<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        self.output += " ";
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "}";
        Ok(())
    }
}

impl ser::SerializeStruct for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output += " ";
        }
        self.output += ":";
        self.output += key;
        self.output += " ";
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "}";
        Ok(())
    }
}

impl ser::SerializeStructVariant for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output += " ";
        }
        self.output += ":";
        self.output += key;
        self.output += " ";
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output += "}";
        Ok(())
    }
}

// Public API

/// Serialize a value to an EDN string.
pub fn edn_to_string<T: Serialize>(value: &T) -> Result<String, EdnError> {
    let mut serializer = EdnSerializer::new();
    value.serialize(&mut serializer)?;
    Ok(serializer.output)
}

/// Deserialize a value from an EDN string.
pub fn edn_from_str<'a, T: Deserialize<'a>>(s: &'a str) -> Result<T, ParseError> {
    let (top, rest) = edn::read(s).map_err(map_edn_error)?;
    let rest = rest.trim();
    if !rest.is_empty() {
        return Err(ParseError::EdnParse(format!(
            "unexpected trailing content: {rest}"
        )));
    }
    T::deserialize(EdnDeserializer(top)).map_err(EdnError::into_parse_error)
}

/// Helper for serializing an EDN keyword (`:name`).
pub struct EdnKeyword<'a>(pub &'a str);

impl Serialize for EdnKeyword<'_> {
    fn serialize<S: ser::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        // When used with EdnSerializer, we want `:name` output.
        // EdnSerializer.serialize_str emits `"name"` (quoted), so we need raw output.
        // For now, use a unit variant trick or just handle in the parent.
        serializer.serialize_str(self.0)
    }
}
