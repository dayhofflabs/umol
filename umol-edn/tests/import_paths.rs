//! Locks the crate-root import surface.
//!
//! Any future seal or rename that breaks these imports is a public API
//! change. Every symbol the crate exposes to downstream consumers must
//! be reachable from `umol_edn::` directly — no submodule paths.

#[allow(unused_imports)]
use umol_edn::{
    edn, from_str, from_str_with, from_value, from_value_ref, read_all, read_all_with,
    read_string, read_string_with, to_string, to_string_pretty, to_string_with, to_value,
    DeError, DuplicateKeyPolicy, Edn, EdnDeserializer, EdnError, EdnHashSet, EdnKeyRef,
    EdnKeyword, EdnList, EdnMap, EdnMapHelper, EdnOwned, EdnSeq, EdnSerializer, EdnSet,
    EdnStreamDeserializer, EdnSymbol, EdnTagged, FormatConfig, FromEdn, Keyword, ParseConfig,
    ParseError, Reader, SerError, StreamDeserializer, Symbol, TagFn, TagReaders, ToEdn, Value,
};

#[cfg(feature = "bignum")]
#[allow(unused_imports)]
use umol_edn::{EdnBigDecimal, EdnBigInt};

#[cfg(feature = "chrono")]
#[allow(unused_imports)]
use umol_edn::inst_to_edn;

#[cfg(feature = "uuid")]
#[allow(unused_imports)]
use umol_edn::uuid_to_edn;

#[test]
fn test_import_paths_compile() {}
