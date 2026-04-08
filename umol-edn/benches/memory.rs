//! Peak-memory bench for the two deserialization paths.
//!
//! Not a criterion bench — criterion measures time, not memory. This is a
//! standalone binary that installs a tracking global allocator, runs each
//! workload once, and reports peak bytes allocated.
//!
//! Question being answered: how much more peak memory does the tree-mediated
//! path (`read_string` → `from_value`) use vs. the direct streaming path
//! (`from_str`) for inputs of various sizes?

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicUsize, Ordering};

use serde::{Deserialize, Serialize};
use umol_edn::de::{from_str, from_value};
use umol_edn::read_string;

// ---------------------------------------------------------------------------
// Tracking allocator
// ---------------------------------------------------------------------------

struct TrackingAlloc;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static BASELINE: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for TrackingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            let n = LIVE.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            let mut peak = PEAK.load(Ordering::Relaxed);
            while n > peak {
                match PEAK.compare_exchange_weak(
                    peak,
                    n,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(cur) => peak = cur,
                }
            }
        }
        p
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout);
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }
}

#[global_allocator]
static ALLOCATOR: TrackingAlloc = TrackingAlloc;

fn reset_peak() {
    let live = LIVE.load(Ordering::Relaxed);
    BASELINE.store(live, Ordering::Relaxed);
    PEAK.store(live, Ordering::Relaxed);
}

fn peak_since_reset() -> usize {
    PEAK.load(Ordering::Relaxed).saturating_sub(BASELINE.load(Ordering::Relaxed))
}

// ---------------------------------------------------------------------------
// Workload
// ---------------------------------------------------------------------------

#[allow(dead_code)]
#[derive(Deserialize, Serialize)]
struct Record {
    atoms: Vec<String>,
    bonds: Vec<(String, String, String)>,
}

const MOLECULE_SMALL: &str = r#"{:atoms ["C" "O"] :bonds [["0" "1" "single"]]}"#;

fn build_input(bytes_target: usize) -> String {
    let unit = MOLECULE_SMALL;
    let count = bytes_target / (unit.len() + 1);
    let mut s = String::with_capacity(count * (unit.len() + 1));
    s.push('[');
    for i in 0..count {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(unit);
    }
    s.push(']');
    s
}

fn measure<F, R>(label: &str, f: F)
where
    F: FnOnce() -> R,
{
    // Force any lazy setup outside the measurement.
    reset_peak();
    let result = f();
    let peak = peak_since_reset();
    black_box(result);
    println!("  {label:40} peak = {:>10} bytes", peak);
}

fn run_sizes() {
    for &target in &[
        10 * 1024,        // 10 KB
        100 * 1024,       // 100 KB
        1024 * 1024,      // 1 MB
        10 * 1024 * 1024, // 10 MB
    ] {
        let input = build_input(target);
        println!("\ninput size: {:>10} bytes ({} records)", input.len(),
                 (input.len() / (MOLECULE_SMALL.len() + 1)).max(1));

        // Path A: direct serde streaming deserializer, never builds a full Edn tree.
        measure("direct  from_str::<Vec<Record>>", || {
            from_str::<Vec<Record>>(black_box(&input)).unwrap()
        });

        // Path B: eager tree construction, then from_value.
        measure("tree    read_string + from_value", || {
            let edn = read_string(black_box(&input)).unwrap();
            from_value::<Vec<Record>>(edn).unwrap()
        });

        // Path C: tree only (no deserialization) — isolates tree cost.
        measure("tree    read_string only", || {
            read_string(black_box(&input)).unwrap()
        });
    }
}

fn main() {
    println!("umol-edn peak-memory bench");
    println!("==========================");
    println!("Tracking allocator: peak bytes above baseline at end of each closure.");
    run_sizes();
}
