//! Content-independent EDN serialization helpers shared across the DSL boundary types.

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer};

use super::refs::AtomRef;

pub(super) fn eof_err() -> EdnError {
    DeError::Custom("unexpected end of input".into()).into()
}

pub(super) fn two_atom_refs(
    atoms: Vec<AtomRef>,
    entry: &'static str,
) -> Result<[AtomRef; 2], DeError> {
    atoms.try_into().map_err(|v: Vec<AtomRef>| {
        DeError::Custom(format!(
            "{entry} :atoms must have exactly 2 atoms, got {}",
            v.len()
        ))
    })
}

pub(super) fn required_key<'e>(
    m: &'e EdnMap<'e>,
    key: &'static str,
    context: &'static str,
) -> Result<&'e Edn<'e>, DeError> {
    m.get_keyword(key).ok_or_else(|| DeError::MissingField {
        key: key.to_string(),
        path: vec![context.into()],
    })
}

pub(super) fn optional_id(m: &EdnMap<'_>) -> Result<Option<String>, DeError> {
    match m.get_keyword("id") {
        Some(Edn::Keyword(k)) => Ok(Some(k.name().to_string())),
        Some(other) => Err(DeError::TypeMismatch {
            expected: "keyword id",
            got: other.kind(),
            path: vec![":id".into()],
        }),
        None => Ok(None),
    }
}

/// The two elements of a `[left right]` vector.
pub(super) fn pair<'a, 'de>(
    edn: &'a Edn<'de>,
    context: &'static str,
) -> Result<(&'a Edn<'de>, &'a Edn<'de>), DeError> {
    if let Edn::Vector(v) = edn {
        if v.len() == 2 {
            return Ok((&v[0], &v[1]));
        }
    }
    Err(DeError::Custom(format!("{context} expects [left right]")))
}

pub(super) fn missing(key: &str, context: &'static str) -> EdnError {
    DeError::MissingField {
        key: key.to_string(),
        path: vec![context.into()],
    }
    .into()
}

pub(super) fn unexpected_byte_kind(b: u8) -> &'static str {
    match b {
        b'"' => "string",
        b':' => "keyword",
        b'[' => "vector",
        b'{' => "map",
        b'0'..=b'9' | b'-' | b'+' => "number",
        _ => "token",
    }
}

pub(crate) fn single_key_map(key: &str, value: Edn<'static>) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(Edn::Keyword(EdnKeyword::owned(key.into())), value);
    Edn::Map(m)
}

pub(super) fn read_vec<T>(
    de: &mut EdnStreamDeserializer<'_>,
    mut read_element: impl FnMut(&mut EdnStreamDeserializer<'_>) -> Result<T, EdnError>,
) -> Result<Vec<T>, EdnError> {
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            break;
        }
        out.push(read_element(de)?);
    }
    Ok(out)
}

pub(super) fn read_map(
    de: &mut EdnStreamDeserializer<'_>,
    mut on_entry: impl FnMut(&mut EdnStreamDeserializer<'_>, &str) -> Result<(), EdnError>,
) -> Result<(), EdnError> {
    de.consume_byte(b'{')?;
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?.into_owned();
        on_entry(de, key.as_str())?;
    }
    Ok(())
}

/// Consume `{:key value}` as a single-key map, returning the key and
/// leaving the stream positioned at the opening-map byte (caller has already
/// read the value). Errors if the map contains more than one key.
pub(super) fn read_single_key_map_header(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<String, EdnError> {
    de.consume_byte(b'{')?;
    Ok(de.read_keyword_name()?.into_owned())
}

pub(super) fn consume_single_key_map_close(
    de: &mut EdnStreamDeserializer<'_>,
    context: &'static str,
) -> Result<(), EdnError> {
    if !de.try_consume_byte(b'}')? {
        return Err(DeError::Custom(format!("{} must have exactly one key", context)).into());
    }
    Ok(())
}

pub(super) fn parse_vec<T>(
    edn: &Edn<'_>,
    context: &'static str,
    mut f: impl FnMut(&Edn<'_>) -> Result<T, DeError>,
) -> Result<Vec<T>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    v.iter().map(|e| f(e)).collect()
}

pub(super) fn parse_single_key_map<'a, 'de>(
    edn: &'a Edn<'de>,
    context: &'static str,
) -> Result<(&'a str, &'a Edn<'de>), DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "single-key map",
            got: edn.kind(),
            path: vec![context.into()],
        });
    };
    if m.len() != 1 {
        return Err(DeError::Custom(format!(
            "{context} must have exactly one key, got {}",
            m.len()
        )));
    }
    let (k, v) = m.iter().next().unwrap();
    let Edn::Keyword(key) = k else {
        return Err(DeError::TypeMismatch {
            expected: "keyword key",
            got: k.kind(),
            path: vec![context.into()],
        });
    };
    Ok((key.name(), v))
}
