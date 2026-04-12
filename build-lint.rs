fn check_comment_art() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let src = std::path::Path::new(&manifest_dir).join("src");
    println!("cargo::rerun-if-changed={}", src.display());
    scan_dir(&src);
}

fn scan_dir(dir: &std::path::Path) {
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            scan_dir(&path);
        } else if path.extension().map_or(false, |e| e == "rs") {
            check_file(&path);
        }
    }
}

fn check_file(path: &std::path::Path) {
    let Ok(src) = std::fs::read_to_string(path) else { return };
    for (i, line) in src.lines().enumerate() {
        if is_comment_art(line) {
            println!(
                "cargo::warning=comment art at {}:{}: `{}`",
                path.display(),
                i + 1,
                line.trim(),
            );
        }
    }
}

fn is_comment_art(line: &str) -> bool {
    let trimmed = line.trim_start();
    let rest = if let Some(r) = trimmed.strip_prefix("///") {
        r
    } else if let Some(r) = trimmed.strip_prefix("//!") {
        r
    } else if let Some(r) = trimmed.strip_prefix("//") {
        r
    } else {
        return false;
    };
    let content = rest.trim();
    if content.len() < 4 {
        return false;
    }
    let first = match content.chars().next() {
        Some(c) if c == '-' || c == '=' => c,
        _ => return false,
    };
    content.chars().all(|c| c == first)
}
