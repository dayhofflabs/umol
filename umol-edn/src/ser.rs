//! Serde `Serializer` for EDN. The compact path writes directly to a `String`;
//! the tree path builds an `Edn<'static>` so the size-aware formatter can lay
//! it out.
use std::borrow::Cow;
use std::fmt::Write;
#[cfg(feature = "bignum")]
use std::str::FromStr;

use ::serde::ser::{
    SerializeMap, SerializeSeq, SerializeStruct, SerializeStructVariant, SerializeTuple,
    SerializeTupleStruct, SerializeTupleVariant,
};
use ::serde::{Serialize, Serializer};
#[cfg(feature = "bignum")]
use bigdecimal::BigDecimal;

use crate::collections::{EdnMap, EdnSeq};
#[cfg(feature = "bignum")]
use crate::edn::EdnBigDecimal;
use crate::edn::{Edn, EdnKeyword};
use crate::error::{EdnError, SerError};
use crate::formatter::FormatConfig;
#[cfg(feature = "bignum")]
use crate::serde::{BIGDECIMAL_TOKEN, BIGINT_TOKEN};
use crate::serde::{KEYWORD_TOKEN, LIST_TOKEN, SET_TOKEN, SYMBOL_TOKEN, TAGGED_TOKEN};

/// Serialize a value to a compact EDN string.
pub fn to_string<T: Serialize>(value: &T) -> Result<String, EdnError> {
    let mut ser = EdnSerializer::new();
    value.serialize(&mut ser)?;
    Ok(ser.output)
}

/// Serialize a value into an `Edn<'static>` tree.
pub fn to_value<T: Serialize>(value: &T) -> Result<Edn<'static>, EdnError> {
    let mut ser = EdnTreeSerializer {
        keyword_mode: false,
        symbol_mode: false,
        set_mode: false,
        list_mode: false,
        #[cfg(feature = "bignum")]
        bigint_mode: false,
        #[cfg(feature = "bignum")]
        bigdec_mode: false,
    };
    value.serialize(&mut ser)
}

/// Serialize a value to a pretty-printed EDN string with default formatting.
pub fn to_string_pretty<T: Serialize>(value: &T) -> Result<String, EdnError> {
    to_string_with(value, &FormatConfig::default())
}

/// Serialize a value to an EDN string formatted with `fmt`. Builds an
/// `Edn<'static>` and lays it out with [`Edn::to_string_with`], so the
/// formatter's width and compaction settings apply uniformly.
pub fn to_string_with<T: Serialize>(value: &T, fmt: &FormatConfig) -> Result<String, EdnError> {
    let edn = to_value(value)?;
    Ok(edn.to_string_with(fmt))
}

/// Serializer that writes compact EDN into a `String`.
pub struct EdnSerializer {
    output: String,
    keyword_mode: bool,
    symbol_mode: bool,
    set_mode: bool,
    list_mode: bool,
    /// When set, the next `serialize_str` writes the inner string raw and
    /// appends `'N'` (EDN bigint suffix).
    #[cfg(feature = "bignum")]
    bigint_mode: bool,
    /// When set, the next `serialize_str` writes the inner string raw and
    /// appends `'M'` (EDN bigdecimal suffix).
    #[cfg(feature = "bignum")]
    bigdec_mode: bool,
    /// Tracks the active state of a tagged tuple-struct (`#tag value`) being
    /// serialized via `serialize_tuple_struct(TAGGED_TOKEN, 2)`.
    tagged_state: TaggedState,
    /// Stack of close delimiters for active sequence contexts. Each
    /// `serialize_seq`/`_tuple`/`_tuple_variant` pushes its close char and
    /// the matching `end()` pops. This lets `#{...}` (set), `(...)` (list),
    /// and `[...]` (vector) share a single `SerializeSeq` impl. The sentinel
    /// `'\0'` means "no close delimiter" (used for tagged tuple structs).
    close_stack: Vec<char>,
}

/// Position within a tagged tuple-struct serialization.
#[derive(Clone, Copy, PartialEq, Eq)]
enum TaggedState {
    None,
    ExpectingTag,
    ExpectingValue,
}

impl EdnSerializer {
    pub fn new() -> Self {
        Self {
            output: String::new(),
            keyword_mode: false,
            symbol_mode: false,
            set_mode: false,
            list_mode: false,
            #[cfg(feature = "bignum")]
            bigint_mode: false,
            #[cfg(feature = "bignum")]
            bigdec_mode: false,
            tagged_state: TaggedState::None,
            close_stack: Vec::new(),
        }
    }
}

impl Default for EdnSerializer {
    fn default() -> Self {
        Self::new()
    }
}

fn write_escaped_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(c),
        }
    }
    out.push('"');
}

impl Serializer for &mut EdnSerializer {
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
        let mut buf = itoa::Buffer::new();
        self.output.push_str(buf.format(v));
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
        let mut buf = itoa::Buffer::new();
        self.output.push_str(buf.format(v));
        Ok(())
    }

    fn serialize_f32(self, v: f32) -> Result<(), Self::Error> {
        self.serialize_f64(f64::from(v))
    }
    fn serialize_f64(self, v: f64) -> Result<(), Self::Error> {
        if v.is_nan() || v.is_infinite() {
            return Err(SerError::UnsupportedFloat("NaN or Infinity").into());
        }
        let mut buf = zmij::Buffer::new();
        let s = buf.format_finite(v);
        self.output.push_str(s);
        if !s.contains('.') && !s.contains('e') && !s.contains('E') {
            self.output.push_str(".0");
        }
        Ok(())
    }

    fn serialize_char(self, v: char) -> Result<(), Self::Error> {
        match v {
            '\n' => self.output.push_str("\\newline"),
            '\r' => self.output.push_str("\\return"),
            ' ' => self.output.push_str("\\space"),
            '\t' => self.output.push_str("\\tab"),
            c if (c as u32) < 0x20 || c == '\u{7F}' => {
                write!(self.output, "\\u{:04X}", c as u32).unwrap();
            }
            c => {
                self.output.push('\\');
                self.output.push(c);
            }
        }
        Ok(())
    }

    fn serialize_str(self, v: &str) -> Result<(), Self::Error> {
        if self.keyword_mode {
            self.keyword_mode = false;
            self.output.push(':');
            self.output += v;
        } else if self.symbol_mode {
            self.symbol_mode = false;
            self.output += v;
        } else {
            #[cfg(feature = "bignum")]
            if self.bigint_mode {
                self.bigint_mode = false;
                self.output += v;
                self.output.push('N');
                return Ok(());
            }
            #[cfg(feature = "bignum")]
            if self.bigdec_mode {
                self.bigdec_mode = false;
                self.output += v;
                self.output.push('M');
                return Ok(());
            }
            write_escaped_str(&mut self.output, v);
        }
        Ok(())
    }

    fn serialize_bytes(self, _v: &[u8]) -> Result<(), Self::Error> {
        Err(SerError::Unsupported("bytes").into())
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
        self.output.push(':');
        self.output += variant;
        Ok(())
    }

    fn serialize_newtype_struct<T: ?Sized + Serialize>(
        self,
        name: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        #[cfg(feature = "bignum")]
        {
            if name == BIGINT_TOKEN {
                self.bigint_mode = true;
                let result = value.serialize(&mut *self);
                self.bigint_mode = false;
                return result;
            }
            if name == BIGDECIMAL_TOKEN {
                self.bigdec_mode = true;
                let result = value.serialize(&mut *self);
                self.bigdec_mode = false;
                return result;
            }
        }
        match name {
            KEYWORD_TOKEN => {
                self.keyword_mode = true;
                let result = value.serialize(&mut *self);
                self.keyword_mode = false;
                result
            }
            SYMBOL_TOKEN => {
                self.symbol_mode = true;
                let result = value.serialize(&mut *self);
                self.symbol_mode = false;
                result
            }
            SET_TOKEN => {
                self.set_mode = true;
                let result = value.serialize(&mut *self);
                self.set_mode = false;
                result
            }
            LIST_TOKEN => {
                self.list_mode = true;
                let result = value.serialize(&mut *self);
                self.list_mode = false;
                result
            }
            _ => value.serialize(self),
        }
    }

    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        self.output.push('#');
        self.output += variant;
        self.output.push(' ');
        value.serialize(self)
    }

    fn serialize_seq(self, _len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        if self.set_mode {
            self.set_mode = false;
            self.output.push_str("#{");
            self.close_stack.push('}');
        } else if self.list_mode {
            self.list_mode = false;
            self.output.push('(');
            self.close_stack.push(')');
        } else {
            self.output.push('[');
            self.close_stack.push(']');
        }
        Ok(self)
    }

    fn serialize_tuple(self, _len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        self.output.push('[');
        self.close_stack.push(']');
        Ok(self)
    }

    fn serialize_tuple_struct(
        self,
        name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleStruct, Self::Error> {
        if name == TAGGED_TOKEN {
            self.output.push('#');
            self.tagged_state = TaggedState::ExpectingTag;
            self.close_stack.push('\0');
            Ok(self)
        } else {
            self.serialize_tuple(len)
        }
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleVariant, Self::Error> {
        self.output.push('#');
        self.output += variant;
        self.output += " [";
        self.close_stack.push(']');
        Ok(self)
    }

    fn serialize_map(self, _len: Option<usize>) -> Result<Self::SerializeMap, Self::Error> {
        self.output.push('{');
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
        self.output.push('#');
        self.output += variant;
        self.output += " {";
        Ok(self)
    }
}

/// Returns `true` if `output` is positioned at the start of a fresh sequence
/// (i.e. the previous char is the open delimiter of a vector, set, or list).
#[inline]
fn at_seq_start(output: &str) -> bool {
    matches!(
        output.as_bytes().last(),
        Some(b'[') | Some(b'{') | Some(b'(')
    )
}

impl SerializeSeq for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !at_seq_start(&self.output) {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        let close = self
            .close_stack
            .pop()
            .expect("SerializeSeq::end without matching open");
        self.output.push(close);
        Ok(())
    }
}

impl SerializeTuple for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !at_seq_start(&self.output) {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        let close = self
            .close_stack
            .pop()
            .expect("SerializeTuple::end without matching open");
        self.output.push(close);
        Ok(())
    }
}

impl SerializeTupleStruct for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        match self.tagged_state {
            TaggedState::ExpectingTag => {
                self.symbol_mode = true;
                value.serialize(&mut **self)?;
                self.symbol_mode = false;
                self.output.push(' ');
                self.tagged_state = TaggedState::ExpectingValue;
                Ok(())
            }
            TaggedState::ExpectingValue => {
                let result = value.serialize(&mut **self);
                self.tagged_state = TaggedState::None;
                result
            }
            TaggedState::None => {
                if !at_seq_start(&self.output) {
                    self.output.push(' ');
                }
                value.serialize(&mut **self)
            }
        }
    }

    fn end(self) -> Result<(), Self::Error> {
        let close = self
            .close_stack
            .pop()
            .expect("SerializeTupleStruct::end without matching open");
        if close != '\0' {
            self.output.push(close);
        }
        Ok(())
    }
}

impl SerializeTupleVariant for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !at_seq_start(&self.output) {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        let close = self
            .close_stack
            .pop()
            .expect("SerializeTupleVariant::end without matching open");
        self.output.push(close);
        Ok(())
    }
}

impl SerializeMap for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_key<T: ?Sized + Serialize>(&mut self, key: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output.push(' ');
        }
        key.serialize(&mut **self)
    }

    fn serialize_value<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        self.output.push(' ');
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push('}');
        Ok(())
    }
}

impl SerializeStruct for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output.push(' ');
        }
        self.output.push(':');
        self.output += key;
        self.output.push(' ');
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push('}');
        Ok(())
    }
}

impl SerializeStructVariant for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        if !self.output.ends_with('{') {
            self.output.push(' ');
        }
        self.output.push(':');
        self.output += key;
        self.output.push(' ');
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push('}');
        Ok(())
    }
}

/// Serializer that materializes a value as an `Edn<'static>` so the size-aware
/// formatter can lay it out. Used by `to_value`, `to_string_pretty`, and
/// `to_string_with`.
pub struct EdnTreeSerializer {
    keyword_mode: bool,
    symbol_mode: bool,
    set_mode: bool,
    list_mode: bool,
    #[cfg(feature = "bignum")]
    bigint_mode: bool,
    #[cfg(feature = "bignum")]
    bigdec_mode: bool,
}

/// Sequence shape requested by the next `serialize_seq` call. Driven by the
/// wrapper-token modes set on `EdnTreeSerializer`.
#[derive(Clone, Copy)]
enum SeqKind {
    Vector,
    Set,
    List,
    /// Tagged tuple struct: first element is the tag (serialized as a
    /// symbol), second is the inner value. Produces `Edn::Tagged`.
    Tagged,
}

fn nan_or_inf_error() -> EdnError {
    SerError::UnsupportedFloat("NaN or Infinity").into()
}

impl<'a> Serializer for &'a mut EdnTreeSerializer {
    type Ok = Edn<'static>;
    type Error = EdnError;
    type SerializeSeq = TreeSeqSerializer<'a>;
    type SerializeTuple = TreeSeqSerializer<'a>;
    type SerializeTupleStruct = TreeSeqSerializer<'a>;
    type SerializeTupleVariant = TreeVariantSeqSerializer<'a>;
    type SerializeMap = TreeMapSerializer<'a>;
    type SerializeStruct = TreeStructSerializer<'a>;
    type SerializeStructVariant = TreeVariantStructSerializer<'a>;

    fn serialize_bool(self, v: bool) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Bool(v))
    }

    fn serialize_i8(self, v: i8) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_i16(self, v: i16) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_i32(self, v: i32) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_i64(self, v: i64) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(v))
    }

    fn serialize_u8(self, v: u8) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_u16(self, v: u16) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_u32(self, v: u32) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Int(i64::from(v)))
    }
    fn serialize_u64(self, v: u64) -> Result<Self::Ok, Self::Error> {
        i64::try_from(v).map(Edn::Int).map_err(|_| {
            SerError::OutOfRange {
                value: v.to_string(),
                target: "i64",
            }
            .into()
        })
    }

    fn serialize_f32(self, v: f32) -> Result<Self::Ok, Self::Error> {
        self.serialize_f64(f64::from(v))
    }
    fn serialize_f64(self, v: f64) -> Result<Self::Ok, Self::Error> {
        if v.is_nan() || v.is_infinite() {
            return Err(nan_or_inf_error());
        }
        Ok(Edn::Float(v))
    }

    fn serialize_char(self, v: char) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Char(v))
    }

    fn serialize_str(self, v: &str) -> Result<Self::Ok, Self::Error> {
        if self.keyword_mode {
            self.keyword_mode = false;
            return Ok(Edn::Keyword(EdnKeyword::owned(v.to_string())));
        }
        if self.symbol_mode {
            self.symbol_mode = false;
            return Ok(Edn::Symbol(crate::edn::EdnSymbol::owned(v.to_string())));
        }
        #[cfg(feature = "bignum")]
        if self.bigint_mode {
            self.bigint_mode = false;
            return num_bigint::BigInt::from_str(v)
                .map(Edn::BigInt)
                .map_err(|e| SerError::Custom(format!("invalid bigint {v:?}: {e}")).into());
        }
        #[cfg(feature = "bignum")]
        if self.bigdec_mode {
            self.bigdec_mode = false;
            return BigDecimal::from_str(v)
                .map(|d| Edn::BigDecimal(EdnBigDecimal::new(d)))
                .map_err(|e| SerError::Custom(format!("invalid bigdecimal {v:?}: {e}")).into());
        }
        Ok(Edn::Str(Cow::Owned(v.to_string())))
    }

    fn serialize_bytes(self, _v: &[u8]) -> Result<Self::Ok, Self::Error> {
        Err(SerError::Unsupported("bytes").into())
    }

    fn serialize_none(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Nil)
    }

    fn serialize_some<T: ?Sized + Serialize>(self, value: &T) -> Result<Self::Ok, Self::Error> {
        value.serialize(self)
    }

    fn serialize_unit(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Nil)
    }

    fn serialize_unit_struct(self, _name: &'static str) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Nil)
    }

    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
    ) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Keyword(EdnKeyword::owned(variant.to_string())))
    }

    fn serialize_newtype_struct<T: ?Sized + Serialize>(
        self,
        name: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Self::Error> {
        #[cfg(feature = "bignum")]
        {
            if name == BIGINT_TOKEN {
                self.bigint_mode = true;
                let result = value.serialize(&mut *self);
                self.bigint_mode = false;
                return result;
            }
            if name == BIGDECIMAL_TOKEN {
                self.bigdec_mode = true;
                let result = value.serialize(&mut *self);
                self.bigdec_mode = false;
                return result;
            }
        }
        match name {
            KEYWORD_TOKEN => {
                self.keyword_mode = true;
                let result = value.serialize(&mut *self);
                self.keyword_mode = false;
                result
            }
            SYMBOL_TOKEN => {
                self.symbol_mode = true;
                let result = value.serialize(&mut *self);
                self.symbol_mode = false;
                result
            }
            SET_TOKEN => {
                self.set_mode = true;
                let result = value.serialize(&mut *self);
                self.set_mode = false;
                result
            }
            LIST_TOKEN => {
                self.list_mode = true;
                let result = value.serialize(&mut *self);
                self.list_mode = false;
                result
            }
            _ => value.serialize(self),
        }
    }

    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Self::Error> {
        let inner = value.serialize(&mut *self)?;
        Ok(Edn::Tagged(
            Cow::Owned(variant.to_string()),
            Box::new(inner),
        ))
    }

    fn serialize_seq(self, len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        let kind = if self.set_mode {
            self.set_mode = false;
            SeqKind::Set
        } else if self.list_mode {
            self.list_mode = false;
            SeqKind::List
        } else {
            SeqKind::Vector
        };
        Ok(TreeSeqSerializer {
            ser: self,
            items: Vec::with_capacity(len.unwrap_or(0)),
            kind,
        })
    }

    fn serialize_tuple(self, len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        Ok(TreeSeqSerializer {
            ser: self,
            items: Vec::with_capacity(len),
            kind: SeqKind::Vector,
        })
    }

    fn serialize_tuple_struct(
        self,
        name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleStruct, Self::Error> {
        let kind = if name == TAGGED_TOKEN {
            SeqKind::Tagged
        } else {
            SeqKind::Vector
        };
        Ok(TreeSeqSerializer {
            ser: self,
            items: Vec::with_capacity(len),
            kind,
        })
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleVariant, Self::Error> {
        Ok(TreeVariantSeqSerializer {
            ser: self,
            tag: variant,
            items: Vec::with_capacity(len),
        })
    }

    fn serialize_map(self, len: Option<usize>) -> Result<Self::SerializeMap, Self::Error> {
        Ok(TreeMapSerializer {
            ser: self,
            map: EdnMap::with_capacity(len.unwrap_or(0)),
            next_key: None,
        })
    }

    fn serialize_struct(
        self,
        _name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStruct, Self::Error> {
        Ok(TreeStructSerializer {
            ser: self,
            map: EdnMap::with_capacity(len),
        })
    }

    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _variant_index: u32,
        variant: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStructVariant, Self::Error> {
        Ok(TreeVariantStructSerializer {
            ser: self,
            tag: variant,
            map: EdnMap::with_capacity(len),
        })
    }
}

pub struct TreeSeqSerializer<'a> {
    ser: &'a mut EdnTreeSerializer,
    items: Vec<Edn<'static>>,
    kind: SeqKind,
}

impl TreeSeqSerializer<'_> {
    fn finish(self) -> Result<Edn<'static>, EdnError> {
        match self.kind {
            SeqKind::Vector => Ok(Edn::Vector(EdnSeq::from(self.items))),
            SeqKind::List => Ok(Edn::List(EdnSeq::from(self.items))),
            SeqKind::Set => {
                let mut set = crate::collections::EdnSet::new();
                for item in self.items {
                    set.insert(item);
                }
                Ok(Edn::Set(set))
            }
            SeqKind::Tagged => {
                let mut iter = self.items.into_iter();
                let tag = match iter.next() {
                    Some(Edn::Symbol(s)) => s.into_cow().into_owned(),
                    Some(other) => {
                        return Err(SerError::Custom(format!(
                            "tagged literal tag must be a symbol, got {other:?}"
                        ))
                        .into());
                    }
                    None => {
                        return Err(
                            SerError::Custom("tagged literal missing tag".to_string()).into()
                        );
                    }
                };
                let inner = iter
                    .next()
                    .ok_or_else(|| SerError::Custom("tagged literal missing value".to_string()))?;
                Ok(Edn::Tagged(Cow::Owned(tag), Box::new(inner)))
            }
        }
    }
}

impl SerializeSeq for TreeSeqSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        let v = value.serialize(&mut *self.ser)?;
        self.items.push(v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        self.finish()
    }
}

impl SerializeTuple for TreeSeqSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        let v = value.serialize(&mut *self.ser)?;
        self.items.push(v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        self.finish()
    }
}

impl SerializeTupleStruct for TreeSeqSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if matches!(self.kind, SeqKind::Tagged) && self.items.is_empty() {
            self.ser.symbol_mode = true;
            let v = value.serialize(&mut *self.ser)?;
            self.ser.symbol_mode = false;
            self.items.push(v);
        } else {
            let v = value.serialize(&mut *self.ser)?;
            self.items.push(v);
        }
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        self.finish()
    }
}

pub struct TreeVariantSeqSerializer<'a> {
    ser: &'a mut EdnTreeSerializer,
    tag: &'static str,
    items: Vec<Edn<'static>>,
}

impl SerializeTupleVariant for TreeVariantSeqSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        let v = value.serialize(&mut *self.ser)?;
        self.items.push(v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Tagged(
            Cow::Owned(self.tag.to_string()),
            Box::new(Edn::Vector(EdnSeq::from(self.items))),
        ))
    }
}

pub struct TreeMapSerializer<'a> {
    ser: &'a mut EdnTreeSerializer,
    map: EdnMap<'static>,
    next_key: Option<Edn<'static>>,
}

impl SerializeMap for TreeMapSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_key<T: ?Sized + Serialize>(&mut self, key: &T) -> Result<(), Self::Error> {
        self.next_key = Some(key.serialize(&mut *self.ser)?);
        Ok(())
    }

    fn serialize_value<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        let key = self
            .next_key
            .take()
            .expect("serialize_value called without serialize_key");
        let v = value.serialize(&mut *self.ser)?;
        self.map.insert(key, v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Map(self.map))
    }
}

pub struct TreeStructSerializer<'a> {
    ser: &'a mut EdnTreeSerializer,
    map: EdnMap<'static>,
}

impl SerializeStruct for TreeStructSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        let v = value.serialize(&mut *self.ser)?;
        self.map
            .insert(Edn::Keyword(EdnKeyword::owned(key.to_string())), v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Map(self.map))
    }
}

pub struct TreeVariantStructSerializer<'a> {
    ser: &'a mut EdnTreeSerializer,
    tag: &'static str,
    map: EdnMap<'static>,
}

impl SerializeStructVariant for TreeVariantStructSerializer<'_> {
    type Ok = Edn<'static>;
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        let v = value.serialize(&mut *self.ser)?;
        self.map
            .insert(Edn::Keyword(EdnKeyword::owned(key.to_string())), v);
        Ok(())
    }

    fn end(self) -> Result<Self::Ok, Self::Error> {
        Ok(Edn::Tagged(
            Cow::Owned(self.tag.to_string()),
            Box::new(Edn::Map(self.map)),
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use ::serde::Serialize;
    use rstest::rstest;

    use crate::edn::{Edn, EdnKeyword};
    use crate::reader::read_string;
    use crate::ser::{to_string, to_string_pretty, to_string_with, to_value, EdnSerializer};
    use crate::FormatConfig;

    #[rstest]
    #[case(true, "true")]
    #[case(false, "false")]
    fn test_serialize_bool(#[case] input: bool, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(12i64, "12")]
    #[case(-1i64, "-1")]
    #[case(0i64, "0")]
    fn test_serialize_i64(#[case] input: i64, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(12u64, "12")]
    #[case(0u64, "0")]
    fn test_serialize_u64(#[case] input: u64, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(3.14f64, "3.14")]
    #[case(1.0f64, "1.0")]
    #[case(12.0f64, "12.0")]
    fn test_serialize_f64(#[case] input: f64, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(f64::NAN)]
    #[case(f64::INFINITY)]
    #[case(f64::NEG_INFINITY)]
    fn test_serialize_f64_error(#[case] input: f64) {
        assert!(to_string(&input).is_err());
    }

    #[rstest]
    #[case("hello", r#""hello""#)]
    #[case("with \"quotes\"", r#""with \"quotes\"""#)]
    #[case("line\nbreak", r#""line\nbreak""#)]
    #[case("tab\there", r#""tab\there""#)]
    #[case("a\rb", r#""a\rb""#)]
    #[case(r"a\b", r#""a\\b""#)]
    fn test_serialize_string(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[test]
    fn test_serialize_none() {
        assert_eq!(to_string(&None::<i64>).unwrap(), "nil");
    }

    #[test]
    fn test_serialize_some() {
        assert_eq!(to_string(&Some(12i64)).unwrap(), "12");
    }

    #[test]
    fn test_serialize_unit() {
        assert_eq!(to_string(&()).unwrap(), "nil");
    }

    #[test]
    fn test_serialize_vec() {
        assert_eq!(to_string(&vec![1, 2, 3]).unwrap(), "[1 2 3]");
    }

    #[test]
    fn test_serialize_empty_vec() {
        assert_eq!(to_string(&Vec::<i64>::new()).unwrap(), "[]");
    }

    #[test]
    fn test_serialize_tuple() {
        assert_eq!(to_string(&(1, 2, 3)).unwrap(), "[1 2 3]");
    }

    #[derive(Serialize)]
    struct Point {
        x: f64,
        y: f64,
    }

    #[test]
    fn test_serialize_struct() {
        let p = Point { x: 1.0, y: 2.0 };
        assert_eq!(to_string(&p).unwrap(), "{:x 1.0 :y 2.0}");
    }

    #[derive(Serialize)]
    #[serde(rename_all = "lowercase")]
    enum Color {
        Red,
        Green,
        Blue,
    }

    #[rstest]
    #[case(Color::Red, ":red")]
    #[case(Color::Green, ":green")]
    #[case(Color::Blue, ":blue")]
    fn test_serialize_unit_variant(#[case] input: Color, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[derive(Serialize)]
    enum Shape {
        Circle(f64),
        Rect { w: f64, h: f64 },
    }

    #[test]
    fn test_serialize_newtype_variant() {
        assert_eq!(to_string(&Shape::Circle(5.0)).unwrap(), "#Circle 5.0");
    }

    #[test]
    fn test_serialize_struct_variant() {
        let s = Shape::Rect { w: 3.0, h: 4.0 };
        assert_eq!(to_string(&s).unwrap(), "#Rect {:w 3.0 :h 4.0}");
    }

    #[derive(Serialize)]
    struct Nested {
        point: Point,
        label: String,
    }

    #[test]
    fn test_serialize_nested() {
        let n = Nested {
            point: Point { x: 3.0, y: 4.0 },
            label: "origin".into(),
        };
        assert_eq!(
            to_string(&n).unwrap(),
            r#"{:point {:x 3.0 :y 4.0} :label "origin"}"#,
        );
    }

    #[test]
    fn test_serialize_newtype_struct() {
        #[derive(Serialize)]
        struct Wrapper(i64);
        assert_eq!(to_string(&Wrapper(12)).unwrap(), "12");
    }

    #[derive(Serialize)]
    enum Tagged {
        Pair(i64, i64),
    }

    #[test]
    fn test_serialize_tuple_variant() {
        assert_eq!(to_string(&Tagged::Pair(1, 2)).unwrap(), "#Pair [1 2]");
    }

    #[test]
    fn test_serialize_map() {
        let mut m = HashMap::new();
        m.insert("a", 1);
        m.insert("b", 2);
        let s = to_string(&m).unwrap();
        assert!(s == r#"{"a" 1 "b" 2}"# || s == r#"{"b" 2 "a" 1}"#);
    }

    #[rstest]
    #[case(7i8, "7")]
    #[case(-1i8, "-1")]
    fn test_serialize_i8(#[case] input: i8, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(300i16, "300")]
    #[case(-300i16, "-300")]
    fn test_serialize_i16(#[case] input: i16, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(100000i32, "100000")]
    fn test_serialize_i32(#[case] input: i32, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(7u8, "7")]
    #[case(255u8, "255")]
    fn test_serialize_u8(#[case] input: u8, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(1000u16, "1000")]
    fn test_serialize_u16(#[case] input: u16, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case(100000u32, "100000")]
    fn test_serialize_u32(#[case] input: u32, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[test]
    fn test_serialize_f32() {
        let s = to_string(&3.14f32).unwrap();
        assert!(s.starts_with("3.14"));
    }

    #[rstest]
    #[case('a', "\\a")]
    #[case('\n', "\\newline")]
    #[case('\r', "\\return")]
    #[case(' ', "\\space")]
    #[case('\t', "\\tab")]
    #[case('x', "\\x")]
    #[case('\x01', "\\u0001")]
    #[case('\x7F', "\\u007F")]
    fn test_serialize_char(#[case] input: char, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[test]
    fn test_serialize_bytes_error() {
        use ::serde::Serializer;
        let mut ser = EdnSerializer::new();
        assert!((&mut ser).serialize_bytes(b"data").is_err());
    }

    #[test]
    fn test_serialize_unit_struct() {
        #[derive(Serialize)]
        struct Marker;
        assert_eq!(to_string(&Marker).unwrap(), "nil");
    }

    #[test]
    fn test_serialize_tuple_struct() {
        #[derive(Serialize)]
        struct Pair(i64, i64);
        assert_eq!(to_string(&Pair(3, 4)).unwrap(), "[3 4]");
    }

    #[test]
    fn test_to_value_primitives() {
        assert_eq!(to_value(&true).unwrap(), Edn::Bool(true));
        assert_eq!(to_value(&12i64).unwrap(), Edn::Int(12));
        assert_eq!(to_value(&1.5f64).unwrap(), Edn::Float(1.5));
        assert_eq!(to_value(&None::<i64>).unwrap(), Edn::Nil);
    }

    #[test]
    fn test_to_value_struct_is_map() {
        let p = Point { x: 1.0, y: 2.0 };
        let edn = to_value(&p).unwrap();
        assert_eq!(edn.as_map().unwrap().len(), 2);
    }

    #[test]
    fn test_to_value_keyword_newtype() {
        let kw = EdnKeyword::new("foo");
        let edn = to_value(&kw).unwrap();
        assert_eq!(edn.as_keyword().unwrap().as_str(), "foo");
    }

    #[test]
    fn test_to_value_unit_variant_is_keyword() {
        let edn = to_value(&Color::Red).unwrap();
        assert_eq!(edn.as_keyword().unwrap().as_str(), "red");
    }

    #[test]
    fn test_to_value_newtype_variant_is_tagged() {
        let edn = to_value(&Shape::Circle(5.0)).unwrap();
        assert!(edn.is_tagged());
    }

    #[test]
    fn test_to_value_tuple_variant_is_tagged_vector() {
        let edn = to_value(&Tagged::Pair(1, 2)).unwrap();
        let (tag, inner) = match edn {
            Edn::Tagged(t, i) => (t, i),
            _ => panic!("expected tagged"),
        };
        assert_eq!(&*tag, "Pair");
        assert!(inner.is_vector());
    }

    #[test]
    fn test_to_string_pretty_short_struct_inline() {
        let p = Point { x: 1.0, y: 2.0 };
        assert_eq!(to_string_pretty(&p).unwrap(), "{:x 1.0 :y 2.0}");
    }

    #[test]
    fn test_to_string_pretty_short_nested_inline() {
        let n = Nested {
            point: Point { x: 3.0, y: 4.0 },
            label: "origin".into(),
        };
        assert_eq!(
            to_string_pretty(&n).unwrap(),
            r#"{:label "origin" :point {:x 3.0 :y 4.0}}"#,
        );
    }

    #[test]
    fn test_to_string_pretty_short_keyword_vec_inline() {
        let v = vec![EdnKeyword::new("a"), EdnKeyword::new("b")];
        assert_eq!(to_string_pretty(&v).unwrap(), "[:a :b]");
    }

    #[test]
    fn test_to_string_pretty_long_vec_wraps() {
        let v: Vec<i64> = (0..50).collect();
        let s = to_string_pretty(&v).unwrap();
        assert!(s.contains('\n'));
    }

    /// Pretty-printing serde values must converge with parsing the compact form
    /// and re-formatting the resulting `Edn` tree. This guards against the two
    /// pretty-print paths drifting apart.
    #[test]
    fn test_pretty_path_matches_parse_then_format() {
        #[derive(Serialize)]
        struct Mol {
            atoms: Vec<String>,
            bonds: Vec<(i64, i64, String)>,
            charge: i64,
        }
        let mol = Mol {
            atoms: vec!["C".into(), "O".into(), "H".into()],
            bonds: vec![(0, 1, "single".into()), (1, 2, "double".into())],
            charge: -1,
        };
        let fmt = FormatConfig {
            line_width: Some(30),
            ..Default::default()
        };

        let via_tree = to_string_with(&mol, &fmt).unwrap();

        let compact = to_string(&mol).unwrap();
        let parsed = read_string(&compact).unwrap();
        let via_parse = parsed.to_string_with(&fmt);

        assert_eq!(via_tree, via_parse);
    }
}
