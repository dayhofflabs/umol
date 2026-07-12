use std::env;
use std::path::PathBuf;

const SOURCES: [&str; 5] = [
    "nauty.c",
    "nautil.c",
    "nausparse.c",
    "schreier.c",
    "naurng.c",
];

const HEADERS: [&str; 5] = [
    "nauty.h",
    "nausparse.h",
    "schreier.h",
    "naurng.h",
    "naututil.h",
];

// Namespace every global exported by the minimal nauty source closure. This
// allows migration tests to link the old bundled binding and this vendored copy
// into one process without one archive satisfying the other's references.
const PREFIXED_SYMBOLS: [&str; 71] = [
    "addgenerator",
    "addpermutation",
    "adjacencies_sg",
    "alloc_error",
    "aresame_sg",
    "breakout",
    "cheapautom_sg",
    "cleanup_sg",
    "comparelab_tr",
    "condaddgenerator",
    "copy_sg",
    "deleteunmarked",
    "dispatch_sparse",
    "distances_sg",
    "distvals",
    "doref",
    "dumpschreier",
    "expandschreier",
    "extra_autom",
    "extra_level",
    "filtercount",
    "findpermutation",
    "fmperm",
    "fmptn",
    "freeschreier",
    "getorbits",
    "getorbitsmin",
    "grouporder",
    "init_sg",
    "isautom_sg",
    "itos",
    "labelorg",
    "longprune",
    "maketargetcell",
    "multcount",
    "nausparse_check",
    "nausparse_freedyn",
    "nautil_check",
    "nautil_freedyn",
    "nauty",
    "nauty_check",
    "nauty_freedyn",
    "nauty_kill_request",
    "nauty_to_sg",
    "newgroup",
    "nextelement",
    "orbjoin",
    "permset",
    "pruneset",
    "put_sg",
    "putstring",
    "ran_init",
    "ran_init_2",
    "ran_init_time",
    "ran_nextran",
    "refine_sg",
    "schreier_check",
    "schreier_fails",
    "schreier_freedyn",
    "schreier_gens",
    "sg_to_nauty",
    "shortprune",
    "sortlists_sg",
    "sparsenauty",
    "targetcell_sg",
    "testcanlab_sg",
    "testcanlab_tr",
    "updatecan_sg",
    "updatecan_tr",
    "writegroupsize",
    "writeperm",
];

fn main() {
    let source_dir = PathBuf::from("nauty");
    let include_dir = PathBuf::from("include");
    let shim = PathBuf::from("src/umol_nauty.c");
    let word_size = target_word_size();

    println!("cargo:rerun-if-changed={}", shim.display());
    println!(
        "cargo:rerun-if-changed={}",
        include_dir.join("umol_nauty.h").display()
    );

    for file in SOURCES.iter().chain(HEADERS.iter()).chain(
        [
            "sorttemplates.c",
            "COPYRIGHT",
            "LICENSE-APACHE",
            "This_is_nauty2_9_3.txt",
        ]
        .iter(),
    ) {
        println!("cargo:rerun-if-changed={}", source_dir.join(file).display());
    }

    let mut shim_build = cc::Build::new();
    shim_build
        .file(&shim)
        .include(&include_dir)
        .include(&source_dir)
        .std("c11")
        .warnings(true)
        .extra_warnings(true);
    configure_nauty(&mut shim_build, word_size);
    shim_build.compile("umol_nauty_shim");

    let mut nauty_build = cc::Build::new();
    nauty_build
        .files(SOURCES.map(|file| source_dir.join(file)))
        .include(&source_dir)
        .std("c99")
        .warnings(false);
    configure_nauty(&mut nauty_build, word_size);
    nauty_build.compile("umol_nauty");
}

fn configure_nauty(build: &mut cc::Build, word_size: &'static str) {
    build
        .define("USE_TLS", None)
        .define("WORDSIZE", Some(word_size));
    for symbol in PREFIXED_SYMBOLS {
        build.define(symbol, Some(format!("umol_nauty_{symbol}").as_str()));
    }
}

fn target_word_size() -> &'static str {
    let target_os = env::var("CARGO_CFG_TARGET_OS").expect("target OS is set by Cargo");
    let pointer_width =
        env::var("CARGO_CFG_TARGET_POINTER_WIDTH").expect("target pointer width is set by Cargo");

    if target_os == "windows" || pointer_width == "32" {
        "32"
    } else {
        "64"
    }
}
