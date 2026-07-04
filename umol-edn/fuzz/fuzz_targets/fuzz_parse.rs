#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_edn::read_string;

fuzz_target!(|data: &str| {
    let _ = read_string(data);
});
