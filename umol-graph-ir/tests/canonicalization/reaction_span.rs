//! Exact compatibility fixtures for canonical reaction spans.
//!
//! The cases freeze the selected union frame for lifecycle changes, every entity family, and a
//! constraint-only change. Generated relabeling and algebraic properties live in the property
//! target; these fixtures preserve the published canonical representation itself.

use rstest::rstest;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{Canonicalize, CanonicalizeContext, ReactionSpan};

const LIFECYCLE: &str = r#"{
    :atoms [{:add "O"} {:modify ["C" "N"]} {:remove "F"} "Cl"]
    :bonds [{:remove [2 3 :single]}
            {:add [0 1 :double]}
            {:modify [1 3 [:single :double]]}]
}"#;

const CONSTITUTION: &str = r#"{
    :atoms ["C" "N" "C" "C" "C" "B" "H" "B"]
    :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#R"}]
    :aromatic-systems [{:atoms [2 3 4] :attrs "*#e3"}]
    :multicenter-bonds [{:modify {:atoms [5 6 7] :attrs ["*#e2" "*#e4"]}}]
    :noncovalent-bonds [{:add {:atoms [0 7] :attrs "Hbd"}}]
}"#;

const STEREO: &str = r#"{
    :atoms ["C" "F" "Cl" "Br" "I" "C" "C"]
    :bonds [[5 6 :double]
            [0 1 :single] [0 2 :single] [0 3 :single] [0 4 :single]
            [5 1 :single] [5 2 :single] [6 3 :single] [6 4 :single]]
    :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th1"]}}]
    :stereo-bonds [{:site 0 :ligands [1 2 3 4] :attrs "Ct1"}]
}"#;

const CONSTRAINT: &str = r#"{:atoms ["C"] :constraints [{:add {:connected {}}}]}"#;

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

#[rstest]
#[case::lifecycle(LIFECYCLE, r#"{:atoms ["Cl" {:remove "F"} {:modify ["C" "N"]} {:add "O"}] :bonds [{:remove [0 1 :single]} {:modify [0 2 [:single :double]]} {:add [2 3 :double]}]}"#)]
#[case::constitution(CONSTITUTION, r#"{:aromatic-systems [{:atoms [4 5 6] :attrs "*#e3"}] :atoms ["H" "B" "B" "C" "C" "C" "C" "N"] :dative-bonds [{:acceptor 7 :attrs "1#R" :donors [3]}] :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["*#e2" "*#e4"]}}] :noncovalent-bonds [{:add {:atoms [1 3] :attrs "Hbd"}}]}"#)]
#[case::stereo(STEREO, r#"{:atoms ["C" "C" "C" "F" "Cl" "Br" "I"] :bonds [[0 3 :single] [1 3 :single] [0 4 :single] [1 4 :single] [0 5 :single] [2 5 :single] [0 6 :single] [2 6 :single] [1 2 :double]] :stereo-atoms [{:modify {:attrs [:ccw :cw] :ligands [3 4 5 6] :site 0}}] :stereo-bonds [{:attrs :e :ligands [3 4 5 6] :site 8}]}"#)]
#[case::constraint(CONSTRAINT, r#"{:atoms ["C"] :constraints [{:add {:connected {}}}]}"#)]
fn test_reaction_span_canonicalize(#[case] input: &str, #[case] expected: &str) {
    let span = input.parse::<ReactionSpan>().unwrap();
    let canonical = span.canonicalize(&context()).unwrap();

    assert_eq!(canonical.to_string(), expected);
}

#[rstest]
#[case::lifecycle(LIFECYCLE, r#"{:atoms ["Cl" {:remove "O"} {:modify ["N" "C"]} {:add "F"}] :bonds [{:remove [1 2 :double]} {:modify [0 2 [:double :single]]} {:add [0 3 :single]}]}"#)]
fn test_reaction_span_canonicalize_reversal(#[case] input: &str, #[case] expected: &str) {
    let span = input.parse::<ReactionSpan>().unwrap();
    let reversed = span
        .to_reaction()
        .reverse()
        .unwrap()
        .to_reaction_span()
        .unwrap()
        .canonicalize(&context())
        .unwrap();

    assert_eq!(reversed.to_string(), expected);
}
