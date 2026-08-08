fn main() {
    // `x` is referenced but never declared — the macro must reject this at compile time.
    let _ = umol_ast::mol! {
        (c: C) - (x),
    };
}
