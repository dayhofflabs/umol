use umol_graph_ir::mol;

fn main() {
    // `^name` port markers belong to `frag!`; `mol!` must reject them.
    let _ = mol!((c: C) - ^x);
}
