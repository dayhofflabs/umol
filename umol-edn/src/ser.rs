//! Serde `Serializer` that writes compact EDN.

use serde::{ser, Serialize};

use crate::error::EdnError;

/// Serialize a value to a compact EDN string.
pub fn to_string<T: Serialize>(value: &T) -> Result<String, EdnError> {
    let mut ser = EdnSerializer::new();
    value.serialize(&mut ser)?;
    Ok(ser.output)
}

/// Serializer that writes EDN into a `String`.
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

fn write_escaped_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000C}' => out.push_str("\\f"),
            _ => out.push(c),
        }
    }
    out.push('"');
}

impl<'a> ser::Serializer for &'a mut EdnSerializer {
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
        if v.is_nan() {
            self.output += "##NaN";
        } else if v == f64::INFINITY {
            self.output += "##Inf";
        } else if v == f64::NEG_INFINITY {
            self.output += "##-Inf";
        } else {
            let s = v.to_string();
            self.output += &s;
            if !s.contains('.') && !s.contains('e') && !s.contains('E') {
                self.output += ".0";
            }
        }
        Ok(())
    }

    fn serialize_char(self, v: char) -> Result<(), Self::Error> {
        self.output.push('\\');
        self.output.push(v);
        Ok(())
    }

    fn serialize_str(self, v: &str) -> Result<(), Self::Error> {
        write_escaped_str(&mut self.output, v);
        Ok(())
    }

    fn serialize_bytes(self, _v: &[u8]) -> Result<(), Self::Error> {
        Err(EdnError::Custom("bytes not supported".to_string()))
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
        self.output.push('#');
        self.output += variant;
        self.output.push(' ');
        value.serialize(self)
    }

    fn serialize_seq(self, _len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        self.output.push('[');
        Ok(self)
    }

    fn serialize_tuple(self, _len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        self.output.push('[');
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
        self.output.push('#');
        self.output += variant;
        self.output += " [";
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

// --- Compound trait impls ---

impl ser::SerializeSeq for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push(']');
        Ok(())
    }
}

impl ser::SerializeTuple for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push(']');
        Ok(())
    }
}

impl ser::SerializeTupleStruct for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push(']');
        Ok(())
    }
}

impl ser::SerializeTupleVariant for &mut EdnSerializer {
    type Ok = ();
    type Error = EdnError;

    fn serialize_field<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), Self::Error> {
        if !self.output.ends_with('[') {
            self.output.push(' ');
        }
        value.serialize(&mut **self)
    }

    fn end(self) -> Result<(), Self::Error> {
        self.output.push(']');
        Ok(())
    }
}

impl ser::SerializeMap for &mut EdnSerializer {
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

impl ser::SerializeStruct for &mut EdnSerializer {
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

impl ser::SerializeStructVariant for &mut EdnSerializer {
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use serde::Serialize;

    use crate::to_string;

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
    #[case(f64::NAN, "##NaN")]
    #[case(f64::INFINITY, "##Inf")]
    #[case(f64::NEG_INFINITY, "##-Inf")]
    fn test_serialize_f64(#[case] input: f64, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[rstest]
    #[case("hello", r#""hello""#)]
    #[case("with \"quotes\"", r#""with \"quotes\"""#)]
    #[case("line\nbreak", r#""line\nbreak""#)]
    #[case("tab\there", r#""tab\there""#)]
    fn test_serialize_string(#[case] input: &str, #[case] expected: &str) {
        assert_eq!(to_string(&input).unwrap(), expected);
    }

    #[test]
    fn test_serialize_char() {
        assert_eq!(to_string(&'a').unwrap(), "\\a");
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
        use std::collections::HashMap;
        let mut m = HashMap::new();
        m.insert("a", 1);
        m.insert("b", 2);
        let s = to_string(&m).unwrap();
        assert!(s == r#"{"a" 1 "b" 2}"# || s == r#"{"b" 2 "a" 1}"#);
    }
}
