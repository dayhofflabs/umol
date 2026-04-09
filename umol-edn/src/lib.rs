//! EDN (Extensible Data Notation) parser, formatter, and serde integration.

#[cfg(feature = "bignum")]
pub mod bigint_owned;
pub mod collections;
pub mod config;
mod display;
pub mod edn;
pub mod error;
pub mod formatter;
pub mod keyword_owned;
pub mod list_owned;
pub mod native;
pub mod parser;
pub mod reader;
pub mod set_owned;
pub mod streaming;
pub mod symbol_owned;
pub mod tagged_owned;
pub mod tags;
pub mod value;

#[cfg(all(feature = "serde", feature = "bignum"))]
pub mod bigdecimal_serde;
#[cfg(all(feature = "serde", feature = "bignum"))]
pub mod bigint_serde;
#[cfg(feature = "serde")]
pub mod de;
#[cfg(feature = "serde")]
pub mod keyword_serde;
#[cfg(feature = "serde")]
pub mod list_serde;
#[cfg(feature = "serde")]
pub mod ser;
#[cfg(feature = "serde")]
pub mod serde_tokens;
#[cfg(feature = "serde")]
pub mod set_serde;
#[cfg(feature = "serde")]
pub mod symbol_serde;
#[cfg(feature = "serde")]
pub mod tagged_serde;
#[cfg(feature = "serde")]
pub mod value_serde;

#[cfg(feature = "bignum")]
pub use bigdecimal::BigDecimal;
#[cfg(feature = "bignum")]
pub use bigint_owned::EdnBigInt;
pub use collections::{EdnKeyRef, EdnMap, EdnSeq, EdnSet};
pub use config::{DuplicateKeyPolicy, ParseConfig, TagFn, TagReaders};
#[cfg(feature = "bignum")]
pub use edn::EdnBigDecimal;
pub use edn::{Edn, Keyword, Symbol};
pub use error::EdnError;
pub use formatter::EdnFormatter;
pub use keyword_owned::EdnKeyword;
pub use list_owned::EdnList;
pub use native::{EdnMapHelper, FromEdn, ToEdn};
pub use set_owned::EdnHashSet;
pub use symbol_owned::EdnSymbol;
pub use tagged_owned::EdnTagged;
pub use value::Value;
#[cfg(feature = "serde")]
pub use formatter::{to_string_pretty, to_string_with};
#[cfg(feature = "bignum")]
pub use num_bigint::BigInt;
pub use reader::{read_all, read_all_with, read_string, read_string_with, Reader};

#[cfg(not(feature = "macros"))]
/// Construct an `Edn<'static>` value from an EDN string literal.
#[macro_export]
macro_rules! edn {
    ($s:literal) => {
        $crate::read_string($s).expect("invalid EDN in edn! macro")
    };
}

#[cfg(feature = "macros")]
pub use umol_edn_macros::{edn, FromEdn, ToEdn};
