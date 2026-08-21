use std::env;
use std::path::PathBuf;

const SOURCES: [&str; 14] = [
    "CoordgenFragmentBuilder.cpp",
    "CoordgenFragmenter.cpp",
    "CoordgenMacrocycleBuilder.cpp",
    "CoordgenMinimizer.cpp",
    "CoordgenTemplates.cpp",
    "sketcherMinimizer.cpp",
    "sketcherMinimizerAtom.cpp",
    "sketcherMinimizerBond.cpp",
    "sketcherMinimizerFragment.cpp",
    "sketcherMinimizerMarchingSquares.cpp",
    "sketcherMinimizerMolecule.cpp",
    "sketcherMinimizerResidue.cpp",
    "sketcherMinimizerResidueInteraction.cpp",
    "sketcherMinimizerRing.cpp",
];

fn main() {
    println!("cargo:rerun-if-changed=coordgen");
    println!("cargo:rerun-if-changed=include/umol_coordgen.h");
    println!("cargo:rerun-if-changed=src/umol_coordgen.cpp");

    if env::var_os("CARGO_FEATURE_NATIVE").is_none() {
        return;
    }

    let source_dir = PathBuf::from("coordgen");
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .files(SOURCES.map(|source| source_dir.join(source)))
        .file("src/umol_coordgen.cpp")
        .include(&source_dir)
        .include("include")
        .define("STATIC_COORDGEN", None)
        .std("c++11")
        .warnings(false)
        .compile("umol_coordgen");
}
