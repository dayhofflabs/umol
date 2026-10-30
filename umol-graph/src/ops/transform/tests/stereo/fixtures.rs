//! Synthetic stereo inputs with explicit stored frames; no external runtime producer.

use rstest::fixture;
use umol_graph_ir::ir::{AtomId, ElementForm, Molecule};
use umol_graph_ir::mol_dsl_concrete;

#[fixture]
pub(super) fn tetrahedral() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
    }"#
    )
}

#[fixture]
pub(super) fn prochiral() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "F" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
pub(super) fn symmetric() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "H" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
pub(super) fn implicit_hydrogen() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]
    }"#
    )
}

#[fixture]
pub(super) fn lone_pair() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["P#n" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th*"}]
    }"#
    )
}

#[fixture]
pub(super) fn hydrogen_and_lone_pair() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["N#h#n" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 [:h 0] [:lp 0]] :attrs "Th*"}]
    }"#
    )
}

#[fixture]
pub(super) fn cis_trans() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct*"}]
    }"#
    )
}

#[fixture]
pub(super) fn imine() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "N#n" "F" "C#h3"]
        :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:lp 1]] :attrs "Ct*"}]
    }"#
    )
}

#[fixture]
pub(super) fn remote_stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "C" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th0"}
                       {:site 2 :ligands [0 8 9 10] :attrs "Th1"}]
    }"#
    )
}

#[fixture]
pub(super) fn small_ring() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [5 [:h 0] 2 [:h 1]] :attrs "Ct*"}]
    }"#
    )
}

#[fixture]
pub(super) fn square_planar() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
    }"#
    )
}

#[fixture]
pub(super) fn evidenced_stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#T1" "C#h" "F" "Cl" "Br" "C#h" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2#C0"] [5 6 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [0 [:h 1] 6 [:h 5]] :attrs "Ct*"}]
    }"#
    )
}

#[fixture]
pub(super) fn partial_stereo() -> Molecule {
    let mut molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C#h" "F" "Cl" "Br" "C#h" "Cl" "C"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2"] [5 6 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th{0,1}"}]
        :stereo-bonds [{:site 4 :ligands [0 [:h 1] 6 [:h 5]] :attrs "Ct?configuration"}]
    }"#
    );
    molecule.atom_mut(AtomId(7)).attributes_mut().element = ElementForm::Undetermined;
    molecule
}
