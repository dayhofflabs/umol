//! Standalone profiling target for the molecule DSL streaming parser.
//!
//! Runs `MoleculeDsl::from_edn_str` on the indole fixture in a tight loop
//! long enough to produce useful samples. Intended to be run under `samply`:
//!
//! ```
//! cargo build --release --example profile_indole -p umol-graph-ir
//! samply record ./target/release/examples/profile_indole
//! ```

use std::hint::black_box;

use umol_edn::FromEdn;
use umol_graph_ir::dsl::MoleculeDsl;

const MOL_INDOLE: &str = r##"{:atoms [[:n "N"] [:c2 "C"] [:c3 "C"] [:c3a "C"] [:c4 "C"] [:c5 "C"] [:c6 "C"] [:c7 "C"] [:c7a "C"]]
 :bonds [[:n :c2 "1"] [:c2 :c3 "1"] [:c3 :c3a "1"] [:c3a :c4 "1"] [:c4 :c5 "1"] [:c5 :c6 "1"] [:c6 :c7 "1"] [:c7 :c7a "1"] [:c7a :n "1"] [:c3a :c7a "1"]]
 :aromatic-systems [{:atoms [:n :c2 :c3 :c3a :c4 :c5 :c6 :c7 :c7a] :attrs "[2,1,1,1,1,1,1,1,1]#e10"}]}"##;

fn main() {
    const ITERATIONS: usize = 1_000_000;
    for _ in 0..ITERATIONS {
        let m = MoleculeDsl::from_edn_str(black_box(MOL_INDOLE)).unwrap();
        black_box(m);
    }
}
