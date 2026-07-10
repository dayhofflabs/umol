use umol_ast::mol;

fn main() {
    // atom `x` and bond `x` collide: atom and bond labels share one namespace.
    let _ = mol!((x: C) -[ x: "1" ]- (o: O));
}
