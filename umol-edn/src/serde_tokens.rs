//! Newtype-struct token strings used by the EDN serde compatibility layer
//! to distinguish wrapper types from their plain serde representations.
//!
//! Each wrapper type calls `serialize_newtype_struct(TOKEN, payload)` so that
//! the EDN serializers (`EdnSerializer`, `EdnTreeSerializer`) can recognize
//! the wrapper and emit the EDN-specific form. Non-EDN serializers simply
//! ignore the token name and serialize the payload transparently, which
//! gives wrapper types a sensible default representation in JSON, YAML, etc.
//!
//! On the deserialization side, wrapper types call
//! `deserialize_newtype_struct(TOKEN, visitor)`. `EdnDeserializer` matches
//! on the token and produces a visitor call appropriate to the wrapper.

pub const KEYWORD_TOKEN: &str = "$edn::keyword";
pub const SYMBOL_TOKEN: &str = "$edn::symbol";
pub const LIST_TOKEN: &str = "$edn::list";
pub const SET_TOKEN: &str = "$edn::set";
pub const TAGGED_TOKEN: &str = "$edn::tagged";
pub const VALUE_TOKEN: &str = "$edn::value";
#[cfg(feature = "bignum")]
pub const BIGINT_TOKEN: &str = "$edn::bigint";
#[cfg(feature = "bignum")]
pub const BIGDECIMAL_TOKEN: &str = "$edn::bigdecimal";

/// Returns `true` if `name` is one of the EDN wrapper newtype-struct tokens.
#[inline]
pub fn is_edn_wrapper_token(name: &str) -> bool {
    matches!(
        name,
        KEYWORD_TOKEN | SYMBOL_TOKEN | LIST_TOKEN | SET_TOKEN | TAGGED_TOKEN
    ) || {
        #[cfg(feature = "bignum")]
        {
            matches!(name, BIGINT_TOKEN | BIGDECIMAL_TOKEN)
        }
        #[cfg(not(feature = "bignum"))]
        {
            false
        }
    }
}
