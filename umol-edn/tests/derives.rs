//! Tests for #[derive(FromEdn)] and #[derive(ToEdn)].

use std::collections::{HashMap, HashSet};

use rstest::rstest;
use umol_edn::{read_string, DeError, Edn, FromEdn, ToEdn};

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
    #[edn(default)]
    items: Vec<i64>,
    #[edn(default)]
    tags: HashSet<String>,
    #[edn(default)]
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
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
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
// Struct fused streaming path
// ---------------------------------------------------------------------------

#[test]
fn test_struct_from_edn_str() {
    let input = r#"{:name "test" :items [1 2]}"#;
    let cfg = Config::from_edn_str(input).unwrap();
    assert_eq!(cfg.name, "test");
    assert_eq!(cfg.items, vec![1, 2]);
}

#[test]
fn test_struct_from_edn_str_defaults() {
    let input = r#"{:name "bare"}"#;
    let cfg = Config::from_edn_str(input).unwrap();
    assert_eq!(cfg.name, "bare");
    assert!(cfg.items.is_empty());
    assert!(cfg.tags.is_empty());
    assert!(cfg.metadata.is_empty());
    assert_eq!(cfg.label, None);
}

#[test]
fn test_struct_from_edn_str_error() {
    let err = Config::from_edn_str("{}").unwrap_err();
    assert!(err.to_string().contains("name"));
}

// ---------------------------------------------------------------------------
// #[edn(default)] fused streaming path
// ---------------------------------------------------------------------------

#[test]
fn test_edn_default_attr_from_edn_str() {
    let task = Task::from_edn_str(r#"{:title "go"}"#).unwrap();
    assert_eq!(task.title, "go");
    assert_eq!(task.priority, Priority(0));

    let task = Task::from_edn_str(r#"{:title "go" :priority 5}"#).unwrap();
    assert_eq!(task.priority, Priority(5));
}

// ---------------------------------------------------------------------------
// Mixed enum fused streaming path
// ---------------------------------------------------------------------------

#[rstest]
#[case::unit(":nil", Value::Nil)]
#[case::newtype("{:int 7}", Value::Int(7))]
#[case::tuple("{:pair [3 5]}", Value::Pair(3, 5))]
#[case::struct_variant("{:named {:name \"x\" :value 10}}", Value::Named { name: "x".into(), value: 10 })]
fn test_mixed_enum_from_edn_str(#[case] input: &str, #[case] expected: Value) {
    assert_eq!(Value::from_edn_str(input).unwrap(), expected);
}

#[test]
fn test_mixed_enum_from_edn_str_roundtrip() {
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
        let edn_str = v.to_edn().to_string();
        let back = Value::from_edn_str(&edn_str).unwrap();
        assert_eq!(v, back);
    }
}

#[test]
fn test_mixed_enum_from_edn_str_unknown_variant() {
    let err = Value::from_edn_str(":bogus").unwrap_err();
    assert!(err.to_string().contains("unknown"));
}

#[test]
fn test_mixed_enum_from_edn_str_wrong_type() {
    assert!(Value::from_edn_str("123").is_err());
}

// ---------------------------------------------------------------------------
// Struct field rename
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
struct Renamed {
    #[edn(rename = "full-name")]
    name: String,
    #[edn(rename = "qty")]
    count: i64,
}

#[test]
fn test_struct_field_rename_from_edn() {
    let input = r#"{:full-name "x" :qty 3}"#;
    let tree = read_string(input).unwrap();
    let r = Renamed::from_edn(&tree).unwrap();
    assert_eq!(r.name, "x");
    assert_eq!(r.count, 3);
}

#[test]
fn test_struct_field_rename_to_edn() {
    let r = Renamed {
        name: "x".into(),
        count: 3,
    };
    let edn = r.to_edn().to_string();
    assert!(edn.contains(":full-name"));
    assert!(edn.contains(":qty"));
    assert!(!edn.contains(":name"));
    assert!(!edn.contains(":count"));
}

#[test]
fn test_struct_field_rename_from_edn_str() {
    let r = Renamed::from_edn_str(r#"{:full-name "y" :qty 7}"#).unwrap();
    assert_eq!(r.name, "y");
    assert_eq!(r.count, 7);
}

// ---------------------------------------------------------------------------
// Transparent newtype (tuple struct)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
#[edn(transparent)]
struct Score(i64);

#[test]
fn test_transparent_tuple_from_edn() {
    let tree = read_string("7").unwrap();
    assert_eq!(Score::from_edn(&tree).unwrap(), Score(7));
}

#[test]
fn test_transparent_tuple_to_edn() {
    assert_eq!(Score(7).to_edn().to_string(), "7");
}

#[test]
fn test_transparent_tuple_from_edn_str() {
    assert_eq!(Score::from_edn_str("7").unwrap(), Score(7));
}

#[test]
fn test_transparent_tuple_roundtrip() {
    let s = Score(13);
    let back = Score::from_edn(&s.to_edn()).unwrap();
    assert_eq!(s, back);
}

// ---------------------------------------------------------------------------
// Transparent newtype (named struct)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
#[edn(transparent)]
struct Label {
    value: String,
}

#[test]
fn test_transparent_named_from_edn() {
    let tree = read_string(r#""hello""#).unwrap();
    let l = Label::from_edn(&tree).unwrap();
    assert_eq!(l.value, "hello");
}

#[test]
fn test_transparent_named_to_edn() {
    let l = Label {
        value: "hello".into(),
    };
    assert_eq!(l.to_edn().to_string(), r#""hello""#);
}

#[test]
fn test_transparent_named_from_edn_str() {
    let l = Label::from_edn_str(r#""world""#).unwrap();
    assert_eq!(l.value, "world");
}

// ---------------------------------------------------------------------------
// deny_unknown_fields
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn)]
#[edn(deny_unknown_fields)]
struct Strict {
    name: String,
    count: i64,
}

#[test]
fn test_deny_unknown_fields_ok() {
    let input = r#"{:name "a" :count 1}"#;
    let tree = read_string(input).unwrap();
    let s = Strict::from_edn(&tree).unwrap();
    assert_eq!(s.name, "a");
    assert_eq!(s.count, 1);
}

#[test]
fn test_deny_unknown_fields_rejects_extra_from_edn() {
    let input = r#"{:name "a" :count 1 :extra true}"#;
    let tree = read_string(input).unwrap();
    let err = Strict::from_edn(&tree).unwrap_err();
    assert!(err.to_string().contains("extra"));
}

#[test]
fn test_deny_unknown_fields_ok_from_edn_str() {
    let s = Strict::from_edn_str(r#"{:name "b" :count 2}"#).unwrap();
    assert_eq!(s.name, "b");
    assert_eq!(s.count, 2);
}

#[test]
fn test_deny_unknown_fields_rejects_extra_from_edn_str() {
    let err = Strict::from_edn_str(r#"{:name "b" :count 2 :extra true}"#).unwrap_err();
    assert!(err.to_string().contains("extra"));
}

// ---------------------------------------------------------------------------
// Container default
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, FromEdn)]
#[edn(default)]
struct Opts {
    verbose: bool,
    count: i64,
    label: String,
}

#[test]
fn test_container_default_empty() {
    let tree = read_string("{}").unwrap();
    let o = Opts::from_edn(&tree).unwrap();
    assert!(!o.verbose);
    assert_eq!(o.count, 0);
    assert_eq!(o.label, "");
}

#[test]
fn test_container_default_partial() {
    let tree = read_string(r#"{:count 5}"#).unwrap();
    let o = Opts::from_edn(&tree).unwrap();
    assert_eq!(o.count, 5);
    assert_eq!(o.label, "");
}

#[test]
fn test_container_default_from_edn_str() {
    let o = Opts::from_edn_str(r#"{:verbose true}"#).unwrap();
    assert!(o.verbose);
    assert_eq!(o.count, 0);
}

// ---------------------------------------------------------------------------
// skip and skip_if
// ---------------------------------------------------------------------------

fn is_zero(v: &i64) -> bool {
    *v == 0
}

#[derive(Debug, Clone, PartialEq, Eq, FromEdn, ToEdn)]
struct Record {
    name: String,
    #[edn(skip)]
    cached: i64,
    #[edn(skip_if = "is_zero")]
    score: i64,
}

#[test]
fn test_skip_field_deser() {
    let tree = read_string(r#"{:name "a" :score 3}"#).unwrap();
    let r = Record::from_edn(&tree).unwrap();
    assert_eq!(r.name, "a");
    assert_eq!(r.cached, 0);
    assert_eq!(r.score, 3);
}

#[test]
fn test_skip_field_ignores_input() {
    // :cached appears in input but is ignored on deser
    let tree = read_string(r#"{:name "a" :cached 99 :score 3}"#).unwrap();
    let r = Record::from_edn(&tree).unwrap();
    assert_eq!(r.cached, 0);
}

#[test]
fn test_skip_field_ser() {
    let r = Record {
        name: "a".into(),
        cached: 99,
        score: 3,
    };
    let edn = r.to_edn().to_string();
    assert!(!edn.contains("cached"));
    assert!(edn.contains(":score"));
}

#[test]
fn test_skip_if_omits_zero() {
    let r = Record {
        name: "b".into(),
        cached: 0,
        score: 0,
    };
    let edn = r.to_edn().to_string();
    assert!(!edn.contains("cached"));
    assert!(!edn.contains("score"));
}

#[test]
fn test_skip_field_from_edn_str() {
    let r = Record::from_edn_str(r#"{:name "c" :score 7}"#).unwrap();
    assert_eq!(r.cached, 0);
    assert_eq!(r.score, 7);
}

// ---------------------------------------------------------------------------
// from_edn / from_edn_str parity
// ---------------------------------------------------------------------------

#[rstest]
#[case::red(":red")]
#[case::dark_blue(":dark-blue")]
#[case::light_green(":light-green")]
fn test_from_edn_str_agrees_with_from_edn_unit_enum(#[case] input: &str) {
    let tree = read_string(input).unwrap();
    let via_tree = Color::from_edn(&tree).unwrap();
    let via_str = Color::from_edn_str(input).unwrap();
    assert_eq!(via_tree, via_str);
}

#[rstest]
#[case::unit(":nil")]
#[case::newtype("{:int 7}")]
#[case::tuple("{:pair [3 5]}")]
#[case::struct_variant("{:named {:name \"x\" :value 10}}")]
fn test_from_edn_str_agrees_with_from_edn_mixed_enum(#[case] input: &str) {
    let tree = read_string(input).unwrap();
    let via_tree = Value::from_edn(&tree).unwrap();
    let via_str = Value::from_edn_str(input).unwrap();
    assert_eq!(via_tree, via_str);
}

#[rstest]
#[case::minimal(r#"{:name "a"}"#)]
#[case::with_items(r#"{:name "a" :items [1 2]}"#)]
#[case::all_fields(r#"{:name "a" :items [1] :tags #{"x"} :metadata {"k" 3} :label "l"}"#)]
fn test_from_edn_str_agrees_with_from_edn_struct(#[case] input: &str) {
    let tree = read_string(input).unwrap();
    let via_tree = Config::from_edn(&tree).unwrap();
    let via_str = Config::from_edn_str(input).unwrap();
    assert_eq!(via_tree, via_str);
}

#[rstest]
#[case::renamed(r#"{:full-name "y" :qty 7}"#)]
fn test_from_edn_str_agrees_with_from_edn_renamed(#[case] input: &str) {
    let tree = read_string(input).unwrap();
    let via_tree = Renamed::from_edn(&tree).unwrap();
    let via_str = Renamed::from_edn_str(input).unwrap();
    assert_eq!(via_tree, via_str);
}

// ---------------------------------------------------------------------------
// Roundtrip: value → to_edn → from_edn → value
// ---------------------------------------------------------------------------

#[test]
fn test_struct_roundtrip() {
    let cfg = Config {
        name: "test".into(),
        items: vec![1, 2],
        tags: ["a".into()].into_iter().collect(),
        metadata: [("k".into(), 3)].into_iter().collect(),
        label: Some("lbl".into()),
    };
    let back = Config::from_edn(&cfg.to_edn()).unwrap();
    assert_eq!(cfg, back);
}

#[test]
fn test_struct_roundtrip_defaults_inflate() {
    // Minimal input → deser → ser → deser produces same value
    let input = r#"{:name "x"}"#;
    let first = Config::from_edn_str(input).unwrap();
    let serialized = first.to_edn().to_string();
    let second = Config::from_edn_str(&serialized).unwrap();
    assert_eq!(first, second);
}

#[test]
fn test_renamed_struct_roundtrip() {
    let r = Renamed {
        name: "a".into(),
        count: 5,
    };
    let back = Renamed::from_edn(&r.to_edn()).unwrap();
    assert_eq!(r, back);
}

#[test]
fn test_transparent_roundtrip_via_str() {
    let s = Score(17);
    let text = s.to_edn().to_string();
    let back = Score::from_edn_str(&text).unwrap();
    assert_eq!(s, back);
}

#[test]
fn test_skip_if_roundtrip() {
    // Non-zero score survives roundtrip
    let r = Record {
        name: "a".into(),
        cached: 99,
        score: 5,
    };
    let back = Record::from_edn(&r.to_edn()).unwrap();
    assert_eq!(back.name, "a");
    assert_eq!(back.cached, 0); // skip → always default
    assert_eq!(back.score, 5);
}

#[test]
fn test_skip_if_without_default_is_lossy() {
    // score=0 is omitted by skip_if, but score is required on deser → error
    let r = Record {
        name: "b".into(),
        cached: 0,
        score: 0,
    };
    assert!(Record::from_edn(&r.to_edn()).is_err());
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
