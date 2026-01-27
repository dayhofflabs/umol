3. Provide GIR-native bond enumeration
Add an API on graph_ir::Molecule that yields all bonds between two atoms without redundant ordering (bond_pairs(atom_idx) or bond_between(a,b) returning iterator).
Refactor topology lint to use this method, removing the HashMap accumulation.
4. Normalize index types
Audit graph_ir::Molecule APIs and lint code to standardize on either usize or u32 for indices; update data structures so conversions disappear.
Adjust the parser → GIR conversion and lints so span lookups use the same index type end-to-end.
5. Standardize diagnostic numbering
Decide globally (and document) whether diagnostics expose 0-based or 1-based IDs.
Change the lints to follow that rule and drop all scattered +1 conversions.
Update existing tests/fixtures to match the chosen convention.
6. Surface spans without SIR fallback
Ensure GIR retains the link to original spans (store them during conversion) so topology lints can emit diagnostics without re-reading SIR data.
Adapt conversion-error reporting to reuse the same GIR-based code path: on failure, build a minimal temporary structure (or surface structured errors) that topology lint can run on or provide span info directly from the error.