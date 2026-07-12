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

    cc::Build::new()
        .file(&shim)
        .include(&include_dir)
        .include(&source_dir)
        .define("USE_TLS", None)
        .define("WORDSIZE", Some(word_size))
        .std("c11")
        .warnings(true)
        .extra_warnings(true)
        .compile("umol_nauty_shim");

    cc::Build::new()
        .files(SOURCES.map(|file| source_dir.join(file)))
        .include(&source_dir)
        .define("USE_TLS", None)
        .define("WORDSIZE", Some(word_size))
        .std("c99")
        .warnings(false)
        .compile("umol_nauty");
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
