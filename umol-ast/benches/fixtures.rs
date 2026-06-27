//! Shared molecule fixtures for the parsing and rendering benchmarks.
//!
//! Included via `#[path = "fixtures.rs"]` from each bench target; not a
//! standalone bench.

use std::sync::LazyLock;

pub const MOL_SMALL: &str = r##"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"##;

pub const MOL_BENZENE: &str = r##"{:atoms ["C" "C" "C" "C" "C" "C"]
 :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
 :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]#e6"}]}"##;

pub const MOL_INDOLE: &str = r##"{:atoms [[:n "N"] [:c2 "C"] [:c3 "C"] [:c3a "C"] [:c4 "C"] [:c5 "C"] [:c6 "C"] [:c7 "C"] [:c7a "C"]]
 :bonds [[:n :c2 "1"] [:c2 :c3 "1"] [:c3 :c3a "1"] [:c3a :c4 "1"] [:c4 :c5 "1"] [:c5 :c6 "1"] [:c6 :c7 "1"] [:c7 :c7a "1"] [:c7a :n "1"] [:c3a :c7a "1"]]
 :aromatic-systems [{:atoms [:n :c2 :c3 :c3a :c4 :c5 :c6 :c7 :c7a] :type "[2,1,1,1,1,1,1,1,1]#e10"}]}"##;

pub const MOL_DIBORANE: &str = r##"{:atoms ["B" "H" "B" "H" "H" "H" "H" "H"]
 :bonds [[0 4 "1"] [0 5 "1"] [2 6 "1"] [2 7 "1"]]
 :multicenter-bonds [{:atoms [0 1 2] :type "[1,0,1]#e2"}
               {:atoms [0 3 2] :type "[1,0,1]#e2"}]}"##;

pub const MOL_WITH_CONSTRAINTS: &str = r##"{:atoms [[:c1 "C"] [:c2 "C"] [:o "O"]]
 :bonds [{:id :b1 :atoms [:c1 :c2] :type "1"} {:id :b2 :atoms [:c2 :o] :type "1"}]
 :constraints [{:connected {:atoms [:c1 :c2 :o]}}
               {:bond-order-sum {:bonds [:b1 :b2] :sum 2}}
               {:not {:atom [:c1 {:valence 3}]}}]}"##;

const LARGE_N: usize = 100;

fn gen_large_molecule(id_every: Option<usize>) -> String {
    let mut atoms = String::from("[");
    for i in 0..LARGE_N {
        if i > 0 {
            atoms.push(' ');
        }
        match id_every {
            Some(every) if i.is_multiple_of(every) => {
                atoms.push_str(&format!("[:a{} \"C\"]", i));
            }
            _ => atoms.push_str("\"C\""),
        }
    }
    atoms.push(']');

    let mut bonds = String::from("[");
    for i in 0..LARGE_N - 1 {
        if i > 0 {
            bonds.push(' ');
        }
        bonds.push_str(&format!("[{} {} \"1\"]", i, i + 1));
    }
    bonds.push(']');

    format!("{{:atoms {} :bonds {}}}", atoms, bonds)
}

pub static MOL_LARGE_NO_IDS: LazyLock<String> = LazyLock::new(|| gen_large_molecule(None));
pub static MOL_LARGE_ALL_IDS: LazyLock<String> = LazyLock::new(|| gen_large_molecule(Some(1)));
pub static MOL_LARGE_PARTIAL_IDS: LazyLock<String> = LazyLock::new(|| gen_large_molecule(Some(10)));
