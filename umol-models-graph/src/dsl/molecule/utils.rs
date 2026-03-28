//! EDN parsing utilities for molecule DSL.

use std::collections::BTreeMap;

use clojure_reader::edn::Edn;

use super::super::error::ParseError;

pub(super) fn map_get<'e>(map: &'e BTreeMap<Edn<'e>, Edn<'e>>, key: &str) -> Option<&'e Edn<'e>> {
    map.iter()
        .find(|(k, _)| matches!(k, Edn::Key(s) if *s == key))
        .map(|(_, v)| v)
}

pub(super) fn extract_map<'e>(
    edn: &'e Edn<'e>,
    ctx: &str,
) -> Result<&'e BTreeMap<Edn<'e>, Edn<'e>>, ParseError> {
    match edn {
        Edn::Map(m) => Ok(m),
        _ => Err(ParseError::InvalidMoleculeMap(format!(
            "expected EDN map for {ctx}"
        ))),
    }
}

pub(super) fn extract_label(edn: &Edn<'_>) -> Result<String, ParseError> {
    match edn {
        Edn::Key(s) => Ok((*s).to_string()),
        _ => Err(ParseError::InvalidMoleculeMap(
            "expected EDN keyword as label".to_string(),
        )),
    }
}

pub(super) fn extract_tagged_str<'e>(edn: &'e Edn<'e>, tag: &str) -> Result<&'e str, ParseError> {
    match edn {
        Edn::Tagged(t, v) if *t == tag => match v.as_ref() {
            Edn::Str(s) => Ok(s),
            _ => Err(ParseError::InvalidMoleculeMap(format!(
                "#{tag} value must be a string"
            ))),
        },
        _ => Err(ParseError::InvalidMoleculeMap(format!(
            "expected #{tag} tagged literal"
        ))),
    }
}

pub(super) fn extract_list<'e, T>(
    map: &'e BTreeMap<Edn<'e>, Edn<'e>>,
    key: &str,
    f: impl Fn(&'e Edn<'e>) -> Result<T, ParseError>,
) -> Result<Vec<T>, ParseError> {
    match map_get(map, key) {
        None => Ok(Vec::new()),
        Some(Edn::Vector(v)) => v.iter().map(f).collect(),
        Some(_) => Err(ParseError::InvalidMoleculeMap(format!(
            ":{key} must be a vector"
        ))),
    }
}
