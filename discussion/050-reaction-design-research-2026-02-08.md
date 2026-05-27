# Reaction Parser Design Research

Research notes for Group 6 (Reactions) design, prior to implementation.

## 1. Non-Strict Atom Mapping in Practice

### Specs and Conventions

**SMIRKS (Daylight)** requires pairwise atom mapping: "every map class has exactly one reactant and one product atom" ([Daylight SMIRKS Tutorial](https://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html)). Stoichiometry is 1-1; non-unit stoichiometry requires repeating reactants/products. Unmapped atoms are created (product) or deleted (reactant).

**Reaction SMARTS** is broader than SMIRKS. The Daylight theory manual notes SMIRKS is a *restricted* subset of reaction SMARTS. Reaction SMARTS (query reactions) may allow more flexible mapping semantics; explicit documentation on duplicate map numbers is sparse.

### Major Data Sources

| Source | Scale | Mapping | Notes |
|--------|-------|---------|-------|
| **USPTO** | ~1M reactions (patents 1976–2016) | 1:1 via rxnmapper | Pipelines *remove* existing maps, then re-map. Output is atom-mapped SMILES. |
| **Pistachio** (NextMove) | 13M+ reactions | 1:1 via NameRxn | 71.5% mapped by NameRxn; fallback to Indigo. |
| **Rhea** | 18k reactions | ChEBI participants | Uses ChEBI IDs for participants; no explicit atom map in SMILES. RXN/RDF formats. |
| **ChEMBL** | Reactions present | Varies | Standard chemical DB; reaction handling less central. |

**Conclusion:** Main reaction datasets (USPTO, Pistachio) use or enforce 1:1 atom mapping. Mapping tools (rxnmapper, NameRxn, Indigo) produce pairwise maps. Non-strict mapping appears rare in curated reaction data; the generic case is more relevant for reaction *templates* (SMARTS) than for reaction *data* (SMILES).

### Toolkits

- **RDKit**: Uses Reaction SMARTS; parses SMIRKS. `react_idx` / `react_atom_idx` track origin. Implements transforms with pairwise mapping.
- **CDK**: `Mapping` class for reactant↔product atom correspondence; `ReactionManipulator` for mapped bonds/atoms.
- **Indigo**: Reaction support; used as fallback mapper in Pistachio.
- **OpenBabel**: Less detailed documentation on reaction mapping.

## 2. Semantic Gap: Variable Composition vs Query vs Markush

### Three Facets

1. **Variable composition** – Set of known structures (R-groups) or unknown from set.
2. **Query features** – SMARTS-like: "class of structures" (e.g. any halogen, any sp² C).
3. **Markush variation types** – ChemAxon / patent literature distinguish:
   - **Substituent variation** – R₁ = H, CH₃, OH, etc.
   - **Position variation** – Variable attachment points.
   - **Frequency variation** – Repeating units, multipliers.
   - **Homology variation** – Generic nodes (CHK, ARY, HET).

### How Others Handle This

**Query features** are well covered by SMARTS (RDKit, OpenBabel, CDK, Indigo). Reaction SMARTS combines query atoms with atom mapping.

**Variable composition** overlaps with R-groups (RG:, LOG: in CXSMILES). ChemAxon Marvin/JChem support R-group definitions, link nodes, repeating units. CXSMILES RG/LOG were removed from umol (code quality); reimplementation pending.

**Markush decomposition** is not well covered by SMARTS. ChemAxon's [Markush features](https://docs.chemaxon.com/display/docs/Homology+Groups+and+Markush+Structures) and [automatic Markush generation](https://chemaxon.com/blog/news/automatic-generation-of-markush-structures-from-specific-compounds) (2019, World Patent Information) address:
- Decomposition of compound series into Markush structures
- R-groups, atom/bond lists, link nodes, repeating units
- Position and homology variation

SMARTS expresses "query" (structural constraints) but not "variable composition" (enumerated set) or "frequency" (repeat counts). The WPI paper and patent tools treat these as separate layers.

**Implication:** A single `ExtendedMolecule` conflates:
- Query (SMARTS – structural class)
- Variable (R-groups – set of options)
- Markush (position, frequency, homology)

A `TemplateReaction` (reaction between template/query molecules) is semantically different from an `ExtendedReaction` (reaction between variable-composition molecules). Mixing them in one type complicates validation and downstream use.

## 3. Chemical Ontologies

### Rhea / ChEBI

- **Rhea**: Expert-curated biochemical reactions. Uses ChEBI for participants (balanced at pH 7.3).
- **ChEBI**: Chemical entity ontology; not reaction-specific.
- **Representation**: Participants as ChEBI IDs; no atom-level mapping in the ontology. Rhea provides RXN, RD files, SMILES (via RDKit) for structures.

### RSC RXNO

- **RXNO**: RSC name reaction ontology; ~500 reaction classes (e.g. Diels–Alder).
- Imports ChEBI.
- Classifies reaction *types*, not atom mappings.

### Summary

Ontologies describe *what* reacts (ChEBI) and *which reaction type* (RXNO), not atom-to-atom mapping. Atom mapping is a cheminformatics concern, not an ontological one. Rhea/ChEBI are useful for biocuration but do not define mapping semantics for SMILES/SMARTS.

## 4. Design Implications

1. **Strict mapping as default**: Aligns with SMIRKS, USPTO, Pistachio, and most tooling. `u32 -> (Option<u32>, Option<u32>)` is sufficient for the common case.

2. **Non-strict mapping**: Defer to `ExtendedReaction` or a separate "template" type if needed. Avoid complicating the primary `Reaction` type.

3. **Reaction vs TemplateReaction**: Consider separate types:
   - `Reaction`: `Molecule` components, strict mapping.
   - `ExtendedReaction`: `ExtendedMolecule` components, variable composition (R-groups, etc.).
   - `TemplateReaction`: Query molecules (SMARTS), for substructure replacement / reaction templates.

4. **Enforcement**: Validation (e.g. stoichiometry, injective mapping) can live in semantic layers; parser produces faithful IR. Strict vs non-strict can be enforced at conversion time.

## References

- [Daylight SMIRKS Tutorial](https://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html)
- [Daylight SMIRKS Theory](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)
- [ReactionUtils USPTO](https://molecularai.github.io/reaction_utils/uspto.html)
- [Pistachio (NextMove)](https://www.nextmovesoftware.com/pistachio.html)
- [Rhea](https://www.rhea-db.org/), [ChEBI](https://www.ebi.ac.uk/chebi/), [RXNO](https://www.ebi.ac.uk/ols4/ontologies/rxno)
- [ChemAxon Markush / Homology](https://docs.chemaxon.com/display/docs/Homology+Groups+and+Markush+Structures)
- [RDKit Reaction SMARTS discussion](https://github.com/rdkit/rdkit/discussions/5168)
