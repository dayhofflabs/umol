//! Tests for #[derive(FromEdn)] and #[derive(ToEdn)].

use std::collections::{HashMap, HashSet};

use rstest::rstest;
use umol_edn::{read_string, Edn, FromEdn, ToEdn};

// ---------------------------------------------------------------------------
// Unit enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
enum Color {
    Red,
    DarkBlue,
    LightGreen,
}

#[rstest]
#[case::red(":red", Color::Red)]
#[case::dark_blue(":dark-blue", Color::DarkBlue)]
#[case::light_green(":light-green", Color::LightGreen)]
fn test_unit_enum_from_edn(#[case] input: &str, #[case] expected: Color) {
    let tree = read_string(input).unwrap();
    assert_eq!(Color::from_edn(&tree).unwrap(), expected);
}

#[rstest]
#[case::red(Color::Red, ":red")]
#[case::dark_blue(Color::DarkBlue, ":dark-blue")]
#[case::light_green(Color::LightGreen, ":light-green")]
fn test_unit_enum_to_edn(#[case] value: Color, #[case] expected_edn: &str) {
    let edn = value.to_edn();
    assert_eq!(edn.to_string(), expected_edn);
}

#[test]
fn test_unit_enum_roundtrip() {
    for color in [Color::Red, Color::DarkBlue, Color::LightGreen] {
        let edn = color.to_edn();
        let back = Color::from_edn(&edn).unwrap();
        assert_eq!(color, back);
    }
}

#[test]
fn test_unit_enum_from_edn_str() {
    let color = Color::from_edn_str(":dark-blue").unwrap();
    assert_eq!(color, Color::DarkBlue);
}

#[test]
fn test_unit_enum_unknown_variant() {
    let tree = read_string(":magenta").unwrap();
    let err = Color::from_edn(&tree).unwrap_err();
    assert!(err.to_string().contains("unknown Color variant"));
}

#[test]
fn test_unit_enum_wrong_type() {
    let tree = read_string("42").unwrap();
    assert!(Color::from_edn(&tree).is_err());
}

// ---------------------------------------------------------------------------
// Unit enum with rename
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
enum Bond {
    #[edn(rename = "single")]
    Single,
    #[edn(rename = "double")]
    Double,
}

#[rstest]
#[case::single(":single", Bond::Single)]
#[case::double(":double", Bond::Double)]
fn test_unit_enum_rename(#[case] input: &str, #[case] expected: Bond) {
    let tree = read_string(input).unwrap();
    assert_eq!(Bond::from_edn(&tree).unwrap(), expected);
    assert_eq!(expected.to_edn().to_string(), input);
}

// ---------------------------------------------------------------------------
// Mixed enum (unit + newtype + tuple + struct variants)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, FromEdn, ToEdn)]
enum Value {
    Nil,
    Int(i64),
    Pair(i64, i64),
    Named { name: String, value: i64 },
}

#[rstest]
#[case::unit(":nil", Value::Nil)]
#[case::newtype("{:int 7}", Value::Int(7))]
#[case::tuple("{:pair [3 5]}", Value::Pair(3, 5))]
#[case::struct_variant("{:named {:name \"x\" :value 10}}", Value::Named { name: "x".into(), value: 10 })]
fn test_mixed_enum_from_edn(#[case] input: &str, #[case] expected: Value) {
    let tree = read_string(input).unwrap();
    assert_eq!(Value::from_edn(&tree).unwrap(), expected);
}

#[test]
fn test_mixed_enum_roundtrip() {
    let values = vec![
        Value::Nil,
        Value::Int(7),
        Value::Pair(3, 5),
        Value::Named {
            name: "x".into(),
            value: 10,
        },
    ];
    for v in values {
        let edn = v.to_edn();
        let back = Value::from_edn(&edn).unwrap();
        assert_eq!(v, back);
    }
}

// ---------------------------------------------------------------------------
// Struct with defaulted fields
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
struct Config {
    name: String,
    items: Vec<i64>,
    tags: HashSet<String>,
    metadata: HashMap<String, i64>,
    label: Option<String>,
}

#[test]
fn test_struct_defaults_missing_fields() {
    let input = "{:name \"test\"}";
    let tree = read_string(input).unwrap();
    let cfg = Config::from_edn(&tree).unwrap();
    assert_eq!(cfg.name, "test");
    assert!(cfg.items.is_empty());
    assert!(cfg.tags.is_empty());
    assert!(cfg.metadata.is_empty());
    assert_eq!(cfg.label, None);
}

#[test]
fn test_struct_defaults_all_present() {
    let input = r#"{:name "test" :items [1 2 3] :tags #{"a"} :metadata {"x" 5} :label "lbl"}"#;
    let tree = read_string(input).unwrap();
    let cfg = Config::from_edn(&tree).unwrap();
    assert_eq!(cfg.name, "test");
    assert_eq!(cfg.items, vec![1, 2, 3]);
    assert!(cfg.tags.contains("a"));
    assert_eq!(cfg.metadata.get("x"), Some(&5));
    assert_eq!(cfg.label, Some("lbl".into()));
}

#[test]
fn test_struct_missing_required_field() {
    let input = "{}";
    let tree = read_string(input).unwrap();
    let err = Config::from_edn(&tree).unwrap_err();
    assert!(err.to_string().contains("name"));
}

// ---------------------------------------------------------------------------
// #[edn(default)] attribute
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct Priority(u8);

impl<'de> FromEdn<'de> for Priority {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, umol_edn::DeError> {
        Ok(Priority(u8::from_edn(edn)?))
    }
}

impl ToEdn for Priority {
    fn to_edn(&self) -> Edn<'static> {
        self.0.to_edn()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, FromEdn)]
struct Task {
    title: String,
    #[edn(default)]
    priority: Priority,
}

#[test]
fn test_edn_default_attr_missing() {
    let input = r#"{:title "do stuff"}"#;
    let tree = read_string(input).unwrap();
    let task = Task::from_edn(&tree).unwrap();
    assert_eq!(task.title, "do stuff");
    assert_eq!(task.priority, Priority(0));
}

#[test]
fn test_edn_default_attr_present() {
    let input = r#"{:title "do stuff" :priority 3}"#;
    let tree = read_string(input).unwrap();
    let task = Task::from_edn(&tree).unwrap();
    assert_eq!(task.priority, Priority(3));
}

// ---------------------------------------------------------------------------
// Struct with fused streaming path
// ---------------------------------------------------------------------------

#[test]
fn test_struct_from_edn_str() {
    let input = r#"{:name "test" :items [1 2]}"#;
    let cfg = Config::from_edn_str(input).unwrap();
    assert_eq!(cfg.name, "test");
    assert_eq!(cfg.items, vec![1, 2]);
}

// ---------------------------------------------------------------------------
// Enum as struct field
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
struct Settings {
    color: Color,
    label: Option<String>,
}

#[test]
fn test_enum_in_struct() {
    let input = r#"{:color :dark-blue :label "foo"}"#;
    let tree = read_string(input).unwrap();
    let s = Settings::from_edn(&tree).unwrap();
    assert_eq!(s.color, Color::DarkBlue);
    assert_eq!(s.label, Some("foo".into()));

    let back = Settings::from_edn(&s.to_edn()).unwrap();
    assert_eq!(s, back);
}
