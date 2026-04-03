//! Built-in tag readers and conversion functions for `#inst` and `#uuid`.

#[cfg(any(feature = "chrono", feature = "uuid"))]
use std::borrow::Cow;
#[cfg(feature = "chrono")]
use std::fmt::Display;

#[cfg(feature = "chrono")]
use chrono::{DateTime, TimeZone};
#[cfg(feature = "uuid")]
use uuid::Uuid;

#[cfg(any(feature = "chrono", feature = "uuid"))]
use crate::edn::Edn;
#[cfg(any(feature = "chrono", feature = "uuid"))]
use crate::error::EdnError;

// -- #inst (RFC 3339 instant) ------------------------------------------------

/// Parse-time reader: validates the string after `#inst` is valid RFC 3339.
#[cfg(feature = "chrono")]
pub(crate) fn read_inst(val: Edn) -> Result<Edn, EdnError> {
    match &val {
        Edn::Str(s) => {
            DateTime::parse_from_rfc3339(s).map_err(|e| EdnError::InvalidInst {
                reason: format!("\"{s}\": {e}"),
            })?;
            Ok(Edn::Tagged(Cow::Borrowed("inst"), Box::new(val)))
        }
        _ => Err(EdnError::InvalidInst {
            reason: "expected string after #inst".into(),
        }),
    }
}

/// Build an `Edn::Tagged("inst", ...)` from a chrono DateTime.
#[cfg(feature = "chrono")]
pub fn inst_to_edn<Tz: TimeZone>(dt: &DateTime<Tz>) -> Edn<'static>
where
    Tz::Offset: Display,
{
    Edn::Tagged(
        Cow::Borrowed("inst"),
        Box::new(Edn::Str(Cow::Owned(dt.to_rfc3339()))),
    )
}

// -- #uuid -------------------------------------------------------------------

/// Parse-time reader: validates the string after `#uuid` is a valid UUID.
#[cfg(feature = "uuid")]
pub(crate) fn read_uuid(val: Edn) -> Result<Edn, EdnError> {
    match &val {
        Edn::Str(s) => {
            Uuid::parse_str(s).map_err(|e| EdnError::InvalidUuid {
                reason: format!("\"{s}\": {e}"),
            })?;
            Ok(Edn::Tagged(Cow::Borrowed("uuid"), Box::new(val)))
        }
        _ => Err(EdnError::InvalidUuid {
            reason: "expected string after #uuid".into(),
        }),
    }
}

/// Build an `Edn::Tagged("uuid", ...)` from a uuid::Uuid.
#[cfg(feature = "uuid")]
pub fn uuid_to_edn(id: &uuid::Uuid) -> Edn<'static> {
    Edn::Tagged(
        Cow::Borrowed("uuid"),
        Box::new(Edn::Str(Cow::Owned(id.to_string()))),
    )
}
