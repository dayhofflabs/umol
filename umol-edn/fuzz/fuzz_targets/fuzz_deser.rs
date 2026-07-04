#![no_main]
use libfuzzer_sys::fuzz_target;
use std::collections::HashMap;
use umol_edn::serde::from_str;

fuzz_target!(|data: &str| {
    let _ = from_str::<Vec<i64>>(data);
    let _ = from_str::<HashMap<String, String>>(data);
    let _ = from_str::<i64>(data);
    let _ = from_str::<String>(data);
    let _ = from_str::<bool>(data);
    let _ = from_str::<f64>(data);
    let _ = from_str::<Vec<String>>(data);
    let _ = from_str::<Option<i64>>(data);
});
