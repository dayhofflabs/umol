- [] Check updated test cases
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
- [] Check u8 -> u32, i8 -> i32, etc. conversions
  from umol_data:
    atomic number: u8
    isotope_mass: u32
    electrons: u8
    charge: i8
    unpaired_e: u8
    valence: u8
    class: u32
    span: u32
  - [] Create consistent naming for length-dependent atom and bond input parsers (should be the longest parsed string, the
    max length parser, 21 for bond, 69 for atom, allows whitespace padding at the end).
  - [] Remove extended_range from bond reacting center. Not sure what these values could mean.
  - [] ctfile::parser::bond: Construction of ExtendedBond
  ```

        let mut bond = ExtendedBond::with_order(order);
        if let Some((stereo_val, direction)) = stereo_direction {
            match order {
                BondOrder::Single => {
                    bond.direction = direction;
                }
                BondOrder::Double => {
                    bond.stereo = stereo_val;
                }
                _ => (),
            }
        }
        bond.topology = topology.flatten();
        bond.reacting_center = reacting_center.flatten();

        Ok((i, (first_atom, second_atom, bond)))

```
  -> why are the atom indices not included? Is with_order() a useful constructor? Seems weird that it only does have a job and I don't see how the SMILES parser would need it. Add a better constructor (all info is available at this point, may need to deconstruct stereo_direction, that's all). Should add an error if the bond order does not support stereo_direction (zero, triple, aromatic).
  Also check basic bond parser.
- [] Make parity tests for basic strict and extended strict parsers. Point of focus: in extended fields
  atom hhh, bbb, mmm, nnn, eee
  bond rrr, ccc
  only 0 or blanks should be accepted in the basic strict parser, 0 or other valid numeric balues should be accepted in the extended strict parser.
- [] Add atom_map_num tests to extended (and possibly basic) parser.
- [] Check if H count should be added to the basic parser. Check interaction with M  HYD property.
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