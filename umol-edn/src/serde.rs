//! Serde integration for EDN.
//!
//! This module contains everything that requires the `serde` feature:
//! serialization/deserialization functions, the EDN serde adapters, and
//! wrapper types that carry EDN-only constructs (keyword, symbol, list,
//! set, tagged literal, bignum) through any serialization format.

// -- Wrapper types (available even without the `serde` feature) -----------

#[cfg(feature = "bignum")]
pub use crate::bigdecimal::EdnBigDecimal;
#[cfg(feature = "bignum")]
pub use crate::bigint::EdnBigInt;
// -- Functions and adapters (require the `serde` feature) -----------------
#[cfg(feature = "serde")]
pub use crate::de::{
    from_str, from_str_with, from_value, from_value_ref, EdnDeserializer,
    StreamDeserializer as EdnStreamDeserializer,
};
pub use crate::dyn_edn::DynEdn;
#[cfg(feature = "serde")]
pub use crate::formatter::{to_string_pretty, to_string_with};
pub use crate::list::EdnList;
#[cfg(feature = "serde")]
pub use crate::ser::{to_string, to_value, EdnSerializer};
pub use crate::set::EdnHashSet as EdnSet;
pub use crate::tagged::EdnTagged;

#[cfg(feature = "serde")]
pub(crate) const KEYWORD_TOKEN: &str = "$edn::keyword";
#[cfg(feature = "serde")]
pub(crate) const SYMBOL_TOKEN: &str = "$edn::symbol";
#[cfg(feature = "serde")]
pub(crate) const LIST_TOKEN: &str = "$edn::list";
#[cfg(feature = "serde")]
pub(crate) const SET_TOKEN: &str = "$edn::set";
#[cfg(feature = "serde")]
pub(crate) const TAGGED_TOKEN: &str = "$edn::tagged";
#[cfg(feature = "serde")]
pub(crate) const VALUE_TOKEN: &str = "$edn::value";
#[cfg(all(feature = "serde", feature = "bignum"))]
pub(crate) const BIGINT_TOKEN: &str = "$edn::bigint";
#[cfg(all(feature = "serde", feature = "bignum"))]
pub(crate) const BIGDECIMAL_TOKEN: &str = "$edn::bigdecimal";
