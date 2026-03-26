umol is a framework for computation chemistry implemented in Rust.

Current goals of the umol project.

1. Representation of graph-based and geometric (3D) chemical objects.
2. Quantum chemistry-informed approach to structural elements with broad chemical compatibility.
3. Emphasis on standard algorithms, minimal number of ad hoc rules.
4. Lean on existing tools, do not reinvent algorithms. Use nom/winnow for parsing, petgraph for graph algorithms, etc.