#![no_main]
use libfuzzer_sys::fuzz_target;
use umol_edn::read_all;

fuzz_target!(|data: &str| {
    let _ = read_all(data);
});
