- Check updated test cases
    moleculelike/mdldraw/18Na.mol → molecule/mdldraw/
    moleculelike/mdldraw/isochg.mol → molecule/mdldraw/
    moleculelike/mdldraw/radiso.mol → molecule/mdldraw/
    moleculelike/ketcher/chain-with-isotope-expected.mol → molecule/ketcher/
    moleculelike/ketcher/chain-with-isotope.mol → molecule/ketcher/
    moleculelike/indigo/t4_R_iso.mol → molecule/indigo/
    moleculelike/indigo/t4_S_iso.mol → molecule/indigo/
    Files moved from invalid/ to moleculelike/
    invalid/indigo/ind-295-biggy1v2000.mol → moleculelike/indigo/
    Contains Q query atom symbol, so correctly belongs in moleculelike/
    invalid/marvin/triglyceride.mol → moleculelike/marvin/
    Contains R1, R2, R3 R-group labels
- Merge ctab into ctfile, move header from mol to ctfile
- Check u8 -> u32, i8 -> i32, etc. conversions
  from umol_data:
    atomic number: u8
    isotope_mass: u32
    electrons: u8
    charge: i8
    unpaired_e: u8
    valence: u8
    class: u32
    span: u32
- Add flag for ignoring positions in CTab files. Should also think about "position-trivial" CTab files (all positions are exactly 0.0)
- Check which parsing functions take `&str` and `&[u8]`. The external representation should probably be like this: `parse_<FORMAT>` should take &str, `parse_<FORMAT>_bytes` should take `&[u8]`. Which type is used internally (most likely `&[u8]`) is irrelevant.
- Check if the formers semantics errors in the ctfile parser are correctly mapped onto diagnostics.
- Check if the number of atom lists is actually used in the code and if its typical usage agrees with expectation.

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