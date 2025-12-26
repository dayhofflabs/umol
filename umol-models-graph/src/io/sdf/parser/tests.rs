use super::*;

#[test]
fn test_data_header() {
    let input = b"> <MELTING.POINT>\n";
    let result = data_header().parse(input);
    assert!(result.is_ok());
    let (_, field_name) = result.unwrap();
    assert_eq!(field_name, "MELTING.POINT");
}

#[test]
fn test_data_header_multiple_spaces() {
    let input = b">  <BOILING.POINT>\n";
    let result = data_header().parse(input);
    assert!(result.is_ok());
    let (_, field_name) = result.unwrap();
    assert_eq!(field_name, "BOILING.POINT");
}

#[test]
fn test_data_value() {
    let input = b"100.5\n\n";
    let result = data_value().parse(input);
    assert!(result.is_ok());
    let (_, value) = result.unwrap();
    assert_eq!(value, "100.5");
}

#[test]
fn test_data_value_multiline() {
    let input = b"Line 1\nLine 2\nLine 3\n\n";
    let result = data_value().parse(input);
    assert!(result.is_ok());
    let (_, value) = result.unwrap();
    assert_eq!(value, "Line 1\nLine 2\nLine 3");
}

#[test]
fn test_data_field() {
    let input = b"> <TEST_FIELD>\n100.5\n\n";
    let result = data_field().parse(input);
    assert!(result.is_ok());
    let (_, (field_name, field_value)) = result.unwrap();
    assert_eq!(field_name, "TEST_FIELD");
    assert_eq!(field_value, "100.5");
}

#[test]
fn test_sdf_delimiter() {
    let input = b"$$$$\n";
    let result = sdf_delimiter().parse(input);
    assert!(result.is_ok());
}

#[test]
fn test_sdf_delimiter_no_newline() {
    let input = b"$$$$";
    let result = sdf_delimiter().parse(input);
    assert!(result.is_ok());
}

#[test]
fn test_data_block() {
    let input = b"> <FIELD1>\nValue1\n\n> <FIELD2>\nValue2\n\n$$$$\n";
    let result = data_block().parse(input);
    assert!(result.is_ok());
    let (remaining, fields): (&[u8], IndexMap<String, String>) = result.unwrap();
    assert_eq!(fields.len(), 2);
    assert_eq!(fields.get("FIELD1"), Some(&"Value1".to_string()));
    assert_eq!(fields.get("FIELD2"), Some(&"Value2".to_string()));
    assert!(remaining.starts_with(b"$$$$"));
}
