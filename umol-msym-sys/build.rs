use std::path::PathBuf;

fn main() {
    let src = PathBuf::from("libmsym/src");

    let sources: Vec<PathBuf> = [
        "basis_function.c",
        "character_table.c",
        "context.c",
        "debug.c",
        "elements.c",
        "equivalence_set.c",
        "geometry.c",
        "linalg.c",
        "msym.c",
        "msym_error.c",
        "permutation.c",
        "point_group.c",
        "rsh.c",
        "subspace.c",
        "symmetrize.c",
        "symmetry.c",
        "symop.c",
    ]
    .iter()
    .map(|f| src.join(f))
    .collect();

    cc::Build::new()
        .files(&sources)
        .include(&src)
        .include("include")
        .std("c99")
        .define("MSYM_EXPORTS_BUILT_AS_STATIC", None)
        .warnings(false)
        .compile("msym");
}
