use umol_data::{Element, SpinMultiplicity};

use crate::atom::{AromaticValence, IsotopeMass};
use crate::graph_ir::atom_pattern::AtomPattern;
use crate::graph_ir::atom::AtomError;

/// Parse a ground atom DSL string into an atom pattern
pub(super) fn parse_ground_atom_dsl(input: &str) -> Result<AtomPattern, AtomError> {
    let s = input.trim();
    let (element, mut index) = parse_element(s)?;
    let mut pattern = AtomPattern {
        element: Some(element),
        ..AtomPattern::default()
    };

    loop {
        while let Some(ch) = s[index..].chars().next() {
            if !ch.is_whitespace() {
                break;
            }
            index += ch.len_utf8();
        }
        let Some(ch) = s[index..].chars().next() else {
            break;
        };
        if ch != '#' {
            return Err(AtomError::UnexpectedTag(ch.to_string()));
        }
        index += '#'.len_utf8();

        let pred = s[index..]
            .chars()
            .next()
            .ok_or_else(|| AtomError::UnexpectedTag("#".to_string()))?;
        index += pred.len_utf8();

        while let Some(ch) = s[index..].chars().next() {
            if !ch.is_whitespace() {
                break;
            }
            index += ch.len_utf8();
        }
        let payload_start = index;
        while let Some(next) = s[index..].chars().next() {
            if next == '#' {
                break;
            }
            index += next.len_utf8();
        }
        let payload = s[payload_start..index].trim();

        match pred {
            'i' => {
                if pattern.isotope_mass.is_some() {
                    return Err(AtomError::DuplicateTag("#i".to_string()));
                }
                pattern.isotope_mass = Some(parse_isotope(payload)?);
            }
            'c' => {
                if pattern.charge.is_some() {
                    return Err(AtomError::DuplicateTag("#c".to_string()));
                }
                pattern.charge = Some(parse_charge(payload)?);
            }
            'h' => {
                if pattern.implicit_hydrogens.is_some() {
                    return Err(AtomError::DuplicateTag("#h".to_string()));
                }
                pattern.implicit_hydrogens = Some(parse_optional_u8(
                    payload,
                    |v| AtomError::InvalidImplicitHydrogens(v),
                )?);
            }
            'n' => {
                if pattern.lone_pairs.is_some() {
                    return Err(AtomError::DuplicateTag("#n".to_string()));
                }
                pattern.lone_pairs = Some(parse_optional_u8(payload, |v| AtomError::InvalidLonePairs(v))?);
            }
            'u' => {
                if pattern.unpaired_electrons.is_some() {
                    return Err(AtomError::DuplicateTag("#u".to_string()));
                }
                pattern.unpaired_electrons = Some(parse_optional_u8(
                    payload,
                    |v| AtomError::InvalidUnpairedElectrons(v),
                )?);
            }
            's' => {
                if pattern.multiplicity.is_some() {
                    return Err(AtomError::DuplicateTag("#s".to_string()));
                }
                let m = parse_optional_u8(payload, |v| AtomError::InvalidMultiplicity(v))?;
                pattern.multiplicity = Some(
                    SpinMultiplicity::from_multiplicity(m)
                        .ok_or_else(|| AtomError::InvalidMultiplicity(m.to_string()))?,
                );
            }
            'v' => {
                if pattern.valence.is_some() {
                    return Err(AtomError::DuplicateTag("#v".to_string()));
                }
                pattern.valence = Some(parse_optional_u8(payload, |v| AtomError::InvalidValence(v))?);
            }
            'd' => {
                if pattern.donated_pairs.is_some() {
                    return Err(AtomError::DuplicateTag("#d".to_string()));
                }
                pattern.donated_pairs =
                    Some(parse_optional_u8(payload, |v| AtomError::InvalidDonatedPairs(v))?);
            }
            'r' => {
                if pattern.accepted_pairs.is_some() {
                    return Err(AtomError::DuplicateTag("#r".to_string()));
                }
                pattern.accepted_pairs =
                    Some(parse_optional_u8(payload, |v| AtomError::InvalidAcceptedPairs(v))?);
            }
            'a' => {
                if pattern.aromatic_valence.is_some() {
                    return Err(AtomError::DuplicateTag("#a".to_string()));
                }
                pattern.aromatic_valence = Some(parse_aromatic(payload)?);
            }
            'm' => {
                if pattern.multicenter_valence.is_some() {
                    return Err(AtomError::DuplicateTag("#m".to_string()));
                }
                pattern.multicenter_valence = Some(parse_optional_u8(
                    payload,
                    |v| AtomError::InvalidMulticenterValence(v),
                )?);
            }
            _ => return Err(AtomError::UnexpectedTag(format!("#{}", pred))),
        }
    }

    Ok(pattern)
}

fn parse_element(s: &str) -> Result<(Element, usize), AtomError> {
    let mut chars = s.char_indices();
    let (_, first) = chars
        .next()
        .ok_or_else(|| AtomError::InvalidElement(s.to_string()))?;
    if !first.is_ascii_uppercase() {
        return Err(AtomError::InvalidElement(s.to_string()));
    }

    let mut end = first.len_utf8();
    for (i, c) in chars {
        if c.is_ascii_lowercase() {
            end = i + c.len_utf8();
        } else {
            break;
        }
    }

    let symbol = &s[..end];
    let element = Element::from_symbol(symbol)
        .ok_or_else(|| AtomError::InvalidElement(symbol.to_string()))?;
    Ok((element, end))
}

fn parse_isotope(payload: &str) -> Result<IsotopeMass, AtomError> {
    if payload.is_empty() {
        return Err(AtomError::InvalidTag("#i".to_string()));
    }
    payload
        .parse::<IsotopeMass>()
        .map_err(|_| AtomError::InvalidTag(format!("#i{}", payload)))
}

fn parse_charge(payload: &str) -> Result<i8, AtomError> {
    if payload.is_empty() {
        return Err(AtomError::InvalidCharge(payload.to_string()));
    }

    let signed: i16 = if let Some(rest) = payload.strip_prefix('+') {
        if rest.is_empty() {
            1
        } else {
            rest.parse::<i16>()
                .map_err(|_| AtomError::InvalidCharge(payload.to_string()))?
        }
    } else if let Some(rest) = payload.strip_prefix('-') {
        if rest.is_empty() {
            -1
        } else {
            -(rest
                .parse::<i16>()
                .map_err(|_| AtomError::InvalidCharge(payload.to_string()))?)
        }
    } else {
        payload
            .parse::<i16>()
            .map_err(|_| AtomError::InvalidCharge(payload.to_string()))?
    };

    i8::try_from(signed).map_err(|_| AtomError::InvalidCharge(payload.to_string()))
}

fn parse_aromatic(payload: &str) -> Result<AromaticValence, AtomError> {
    if payload.is_empty() {
        return Ok(AromaticValence::Valence(1));
    }
    payload
        .parse::<AromaticValence>()
        .map_err(|_| AtomError::InvalidAromaticValence(payload.to_string()))
}

fn parse_optional_u8<F>(payload: &str, err: F) -> Result<u8, AtomError>
where
    F: Fn(String) -> AtomError,
{
    if payload.is_empty() {
        return Ok(1);
    }
    payload.parse::<u8>().map_err(|_| err(payload.to_string()))
}
