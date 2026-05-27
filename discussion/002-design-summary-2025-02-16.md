# Core Computational Chemistry Infrastructure Design

## Design goals

- Small core:
   1. structural data types
   2. structure generators (compositional and random)
   3. I/O
   4. structure converters
   5. substructure matching
   6. structure manipulation
   7. canonicalization
   8. chemical reaction
- Interpolate between cheminformatics, and computational/quantum chemistry.
  Interconversion between graph and 3D structures
- Different models of chemical structure allowed, e.g., discrete graphs, 3D structures with
  and without bonding information, fractional bond orders, with and without valence
  constraints, etc.
- Use well-defined algorithms, for example graph algorithms, instead of heuristic implementations
  to the extent possible
- Interface based: queryable capabilities of structure representations
- Static analysis, canonicalization, linting (with error codes)
- Reasonably efficient, scale of 10^4-10^6 structures, no very high-throughput
    optimization needed
- Extensible: New representations and features can be implemented as extensions or plugins

## Non-goals

- Coerce all chemical structures into a single internal representation. Explicit
  lossy conversion is preferable to unphysically unified model.
  Non-conformance to a unified model is a reality of chemistry, not a user error.
- Reproduce a significant portion of RDKit's functionality: fingerprints, descriptors, ML,
  structure decompositions and a lot more are outside the scope.
- Interface with C/C++ ABI of RDKit directly.
- Complete compatibility with existing formats, including ability to parse 100% of
    existing data resources without errors (SMILES, SMARTS, SDF/MOL, PDB)
- Include electronic structure, nuclear quantum effects

## Future goals

- Handling non-main group elements
- Integration with molecular structure manipulation: structure builders, optimization,
  conformer generation, MD, MC, transition state search
- Handling of ensembles of structures
- Higher throughput (> 10^6 structures)
- Principled aromaticity model, Clar structures-based (Green) or eigendecomposition-based?

## Minimal viable implementation

- Only discrete graph model
- Input/output to SMILES / SMARTS formats, potentially lossy
- Canonicalization
- Graph substructure search using established algorithm, try to avoid heuristics
- Graph reaction transformations using double push-out built on the graph substructure search
- Some aromaticity models