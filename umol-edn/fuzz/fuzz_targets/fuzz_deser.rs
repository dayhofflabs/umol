#![no_main]
use libfuzzer_sys::fuzz_target;
use std::collections::HashMap;

fuzz_target!(|data: &str| {
    let _ = umol_edn::de::from_str::<Vec<i64>>(data);
    let _ = umol_edn::de::from_str::<HashMap<String, String>>(data);
    let _ = umol_edn::de::from_str::<i64>(data);
    let _ = umol_edn::de::from_str::<String>(data);
    let _ = umol_edn::de::from_str::<bool>(data);
    let _ = umol_edn::de::from_str::<f64>(data);
    let _ = umol_edn::de::from_str::<Vec<String>>(data);
    let _ = umol_edn::de::from_str::<Option<i64>>(data);
});
