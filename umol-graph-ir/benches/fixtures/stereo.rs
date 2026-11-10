//! Stereo inputs shared by molecule parsing and rendering benchmarks.

pub const MOL_TETRAHEDRAL_EXPLICIT: &str = r#"{:atoms ["C" "F" "Cl" "Br" "I"]
 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#;

pub const MOL_CIS_TRANS_EXPLICIT: &str = r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
 :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
 :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#;
