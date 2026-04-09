//! Locks the crate-root import surface.
//!
//! Any future seal or rename that breaks these imports is a public API
//! change. Every symbol the crate exposes to downstream consumers must
//! be reachable from `umol_edn::` directly or via `umol_edn::serde::`.

// -- Root: native API (no feature flags required) -------------------------

#[cfg(feature = "chrono")]
#[allow(unused_imports)]
use umol_edn::inst_to_edn;
// -- serde module: functions and adapters (require `serde` feature) --------
#[allow(unused_imports)]
use umol_edn::serde::{
    from_str, from_str_with, from_value, from_value_ref, to_string, to_string_pretty,
    to_string_with, to_value, EdnDeserializer, EdnSerializer, EdnStreamDeserializer as SerdeStream,
};
// -- serde module: wrapper types (always available) -----------------------
#[allow(unused_imports)]
use umol_edn::serde::{DynEdn, EdnList, EdnSet as SerdeEdnSet, EdnTagged};
// -- Feature-gated re-exports ---------------------------------------------
#[cfg(feature = "bignum")]
#[allow(unused_imports)]
use umol_edn::serde::{EdnBigDecimal, EdnBigInt};
#[cfg(feature = "uuid")]
#[allow(unused_imports)]
use umol_edn::uuid_to_edn;
#[allow(unused_imports)]
use umol_edn::{
    edn, read_all, read_all_with, read_string, read_string_with, DeError, DuplicateKeyPolicy, Edn,
    EdnError, EdnKeyRef, EdnKeyword, EdnMap, EdnMapHelper, EdnSeq, EdnSet, EdnStreamDeserializer,
    EdnSymbol, FormatConfig, FromEdn, ParseConfig, ParseError, Reader, SerError, TagFn, TagReaders,
    ToEdn,
};

#[test]
fn test_import_paths_compile() {}
