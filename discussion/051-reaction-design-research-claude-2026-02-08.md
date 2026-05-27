# Reaction representation design decisions for a Rust cheminformatics library

From Claude Opus 4.6 Research Mode

**SMIRKS and DPO graph rewriting impose identical bijectivity constraints on atom mapping, making a DPO span the natural internal representation for reaction rules in umol.** Real-world reaction databases almost universally contain non-bijective, partial, or erroneous atom mappings — meaning the library must parse permissively but represent strictly. No existing chemical ontology formally distinguishes the four key semantic categories (specific molecule, query pattern, Markush structure, reaction template), and among toolkits only ChemAxon treats Markush as a first-class entity. The design gap is significant: no Rust library handles reactions at all today.

---

## 1. Non-strict atom mappings dominate real-world databases

The most important empirical finding for umol's design is that **strict bijective atom mappings are the exception, not the rule**, in practice. Every major reaction dataset contains non-bijective mappings either by design or due to errors.

**USPTO (Lowe's extraction)** uses EPAM's Indigo toolkit for auto-mapping. The extraction pipeline requires all product atoms to be accounted for, but reactant atoms may remain unmapped (leaving groups, spectators). This makes the mapping injective from products→reactants at best, never a full bijection. A benchmarking study by Lin et al. (2021) tested 1,851 manually curated reactions against five mapping tools and found the best tool (RXNMapper) achieved only **83.74% accuracy**, implying ~16% of mappings contain errors — including non-bijective ones. The original Indigo mappings are described on Figshare as "wrong in many cases and hence should not be entirely relied on."

**RXNMapper** (Schwaller et al., *Science Advances* 2021) produces injective mappings by construction: its attention-guided algorithm zeroes out already-assigned attention weights (`zero_set_r=True`, `zero_set_p=True`), preventing duplicate assignments. But unmapped atoms (leaving groups, reagents) remain common. The output is an injection from product atoms to reactant atoms — **not a bijection**.

**NameRxn** (NextMove Software) was not originally designed for atom mapping and produces mechanistically-informed mappings only for its ~720 recognized reaction types. For unrecognized reactions (28.5% of Pistachio), it falls back to Indigo. **ChemAxon's automapper** offers three modes: "Complete" (attempts bijection), "Changing" (only reaction center atoms mapped), and "Matching" (atoms on both sides mapped). Only "Complete" mode aims for bijectivity.

**Reaxys and SciFinder/CAS do not systematically store atom mappings.** In both systems, atom mapping is primarily a query tool — users draw mapped atoms to specify reaction-center searches. The underlying databases store structures and conditions but not per-atom correspondences. This means manually-curated databases are effectively unmapped from a computational perspective.

**The Open Reaction Database (ORD)** stores atom mappings embedded within reaction SMILES strings but performs no bijectivity validation. The `validations.py` module checks SMILES parseability and schema compliance, not mapping correctness.

### Toolkit validation behavior is inconsistent

**RDKit** provides the strictest validation. Its `validate()` method reports duplicate atom mapping numbers as errors and unmapped atoms as warnings. Critically, **`RunReactants()` crashes with an invariant violation** (in `Reaction.cpp` line 173) when encountering duplicate atom class indices. RDKit effectively requires bijective mapping for reaction execution but tolerates partial mapping. The utility functions `RemoveUnmappedReactantTemplates()` and `RemoveUnmappedProductTemplates()` handle partially-mapped reactions by moving unmapped components to agents.

**CDK** stores atom mapping as integer properties on atoms (`CDKConstants.ATOM_ATOM_MAPPING`) with **no validation for bijectivity**. The mapping is simply a stored integer — any value, any duplication.

One intentional source of non-bijective mappings deserves special attention: **reaction SMARTS for searching**. The Daylight documentation shows examples like `C[C:2](=[O:1])[OH:1]>>C[C:2](=[O:1])OCC` where atom class `:1` appears on two oxygen atoms in the reactant — a deliberate non-bijective mapping expressing "either of these oxygens." This is a legitimate pattern-matching use case, distinct from SMIRKS transforms.

### Design implication

umol should define a strict internal `ReactionRule` type enforcing bijective atom mapping (the SMIRKS/DPO invariant), plus a separate `ReactionRecord` type that tolerates partial, non-bijective, or absent mappings for database interchange. Parsing from reaction SMILES should emit diagnostics classifying the mapping as bijective, injective, partial, or inconsistent.

---

## 2. Every major toolkit conflates molecules with queries — only ChemAxon treats Markush seriously

The central design tension across toolkits is whether to use a single molecule type (convenient but loses type safety) or separate types (safe but complicates the API). Existing approaches cluster into four patterns.

### The unified-type approach (RDKit)

RDKit uses a **single `ROMol` class** (with `RWMol` as the mutable variant) for both specific molecules and query molecules. When parsing SMILES, atoms are concrete `Atom` objects; when parsing SMARTS, atoms are `QueryAtom` objects (a subclass of `Atom` that overrides `Match()`). Both coexist in the same `ROMol` container. The PostgreSQL cartridge makes this explicit with separate `mol` and `qmol` types, but at the C++ level they share the same class.

This design means **nothing in the type system prevents computing molecular weight on a SMARTS pattern**, which would be meaningless. RDKit compensates with runtime checks and the `AdjustQueryProperties()` / `ReplaceAtomWithQueryAtom()` conversion functions.

R-group decomposition uses `RGroupDecomposition`, which takes `ROMol` cores and molecules and outputs dictionaries of `ROMol` fragments. R-groups are represented as molecules with dummy atoms at attachment points. **RDKit has no native Markush representation** — no position variation, frequency variation, or homology groups. SubstanceGroups (SRU polymer support) provide partial coverage for repeating units.

### The interface-based approach (CDK)

CDK separates `IAtomContainer` (molecules) from `IQueryAtomContainer` (query molecules) at the interface level. Modern CDK (≥2.2) uses an **`Expr` expression tree** for query features — a predicate tree with types like `IS_AROMATIC`, `ELEMENT`, `DEGREE`, `AND`, `OR`, `NOT`. This replaced an older design with dozens of specialized atom classes (`AliphaticAtom`, `AnyAtom`, `DegreeAtom`, etc.) — a combinatorial explosion that the expression tree elegantly solves.

CDK supports multiple SMARTS flavors through an enum: `FLAVOR_DAYLIGHT`, `FLAVOR_OECHEM`, `FLAVOR_CACTVS`, `FLAVOR_CDK_LEGACY`. This flexibility comes at the cost of more complex matching logic. CDK has **no dedicated Markush support**.

### The sibling-hierarchy approach (OEChem)

OEChem has the most explicitly designed type hierarchy: `OEMolBase` (abstract) branches into `OEGraphMol` (single conformer), `OEMol` (multi-conformer via `OEMCMolBase`), and `OEQMol` (query molecule via `OEQMolBase`). These are **siblings, not parent-child** — sharing a base interface but extending it differently. Functions accepting `OEMolBase` work on any molecule type; functions requiring query features take `OEQMolBase`.

OEChem tracks three distinct query origins (SMARTS, MDL query, converted-from-molecule) and uses different expression trees for each. Its design philosophy favors **free functions over methods**: molecules are "primarily data containers" and algorithms are external. **No Markush support.**

### The Markush-first approach (ChemAxon)

ChemAxon is unique in treating Markush structures as first-class entities. Its `RgMolecule` class stores a root scaffold (`MoleculeGraph`) plus R-group definitions, supporting all four Markush variation types:

- **Substituent variation**: Nested R-groups to arbitrary depth (up to R32767), multiple attachment points per R-atom
- **Position variation**: Variable attachment points implemented as multicenter group bonds
- **Frequency variation**: Repeating units with repetition ranges and link nodes
- **Homology variation**: Built-in groups (alkyl, aryl, alkenyl, cycloalkyl, etc.) as pseudo-atoms with structural-feature-based matching — **unique to ChemAxon**

ChemAxon's Markush Composer algorithm (Kovács et al., *World Patent Information* 2019, DOI: 10.1016/j.wpi.2019.03.006) automatically generates Markush structures from sets of specific compounds using MCS-based scaffold recognition and recursive R-group generation.

### The file-format approach (BIOVIA/MDL)

BIOVIA's approach centers on MDL's RG file format (V2000/V3000), which stores scaffold plus R-group member definitions. Pipeline Pilot provides dedicated components for RG file I/O and Markush library enumeration. The strength is interoperability with the vast installed base of MDL-format databases.

### Design recommendation for umol

A Rust enum discriminating between semantic categories maps naturally to the type system:

```rust
enum ChemicalEntity {
    Molecule(Molecule),           // ground-truth labeled graph
    QueryPattern(QueryMolecule),  // SMARTS-like predicate graph
    MarkushStructure(Markush),    // parameterized graph family
    ReactionRule(DPORule),        // graph transformation rule
}
```

CDK's `Expr` tree is the gold standard for query atom/bond predicates — a Rust enum-based expression tree would be natural and efficient. ChemAxon's `RgMolecule` shows Markush should be a **composition** (scaffold + R-group table), not an inheritance hierarchy.

---

## 3. No ontology formally captures the four semantic categories

The gap between chemical ontology work and computational representation needs is stark. **No existing ontology formally distinguishes between specific molecules, query patterns, Markush structures, and reaction templates** as separate ontological categories.

**ChEBI** comes closest to a molecule/class distinction through its naming convention: singular terms (e.g., "phenol," CHEBI:15882) denote specific compounds while plural terms (e.g., "phenols," CHEBI:33853) denote classes. Its three sub-ontologies — Molecular Structure, Role, and Subatomic Particle — classify entities by structure and function. But ChEBI's class definitions are manually curated natural-language descriptions, not SMARTS-like structural patterns that enable automated classification. ChEBI has no concepts for variable-composition structures or query patterns.

**CHEMINF** (Chemical Information Ontology) operates at a meta-level: it describes *information entities about* chemicals, not chemicals themselves. Its key contribution is the IAO `is_about` relationship — a SMILES string `is_about` a molecule. This cleanly separates the thing in the world from its representation. But CHEMINF defines no concept for "query pattern" or "reaction template" as distinct information types.

**InChI** represents specific molecular structures only, explicitly excluding "polymers, molecular class representations (Markush structures), and conformations." **RInChI** (Grethe et al., *J. Cheminformatics* 2018) extends InChI to specific reactions, separating reactants, products, and agents. Version 1.00 lacks atom mapping; v1.2 plans a `/MapAuxInfo` layer. Two experimental extensions are under development at the InChI Trust: **MarkInChI** (encoding Markush-like structures using Zz pseudo-atoms) and **VInChI** (compressing a set of related InChIs into a compact representation of their differences). Neither is part of the standard specification.

**CML** (Chemical Markup Language) has XML elements for `<molecule>` and `<reaction>` (with CMLReact extensions for reaction schemes), but **no elements for query patterns, Markush structures, or reaction templates** as distinct from specific reactions.

**IUPAC** defines "molecular entity," "chemical substance," and "chemical reaction" in the Gold Book, but has no published recommendations defining computational categories like "reaction template" or "molecular query pattern." The InChI/RInChI standards are IUPAC's primary computational outputs.

### The Daylight taxonomy remains the only complete formal framework

Only the Daylight SMILES/SMARTS/SMIRKS lineage provides a complete computational taxonomy:

| Category | Formalism | Graph-theoretic definition |
|---|---|---|
| Specific molecule | SMILES | Labeled graph G = (V, E) with atom/bond labels |
| Query pattern | SMARTS | Graph query with predicate functions, matched via subgraph isomorphism |
| Variable-composition structure | *(not formalized)* | Parameterized graph family G(R₁, R₂, ...) |
| Reaction template | SMIRKS | Graph transformation rule: pattern + rewrite + atom correspondence |

Markush structures remain the least formalized category — neither Daylight nor any ontology provides a complete formal treatment. For umol, this suggests defining Markush semantics carefully in the library's own type system, since no external standard provides guidance.

---

## 4. SMIRKS and DPO graph rewriting are isomorphic formalisms

The most actionable finding for umol's architecture is that **SMIRKS atom mapping semantics and DPO graph rewriting constraints are mathematically identical**, making DPO the natural internal representation.

### SMIRKS constraints verified

The Daylight specification states three key rules. First, "reactant and product atoms which are atom mapped must be mapped pairwise and be complete" — **bijective on the mapped subset**. Second, "any atoms which are not atom mapped are assumed to be added or deleted during the transformation" — unmapped reactant atoms are deleted, unmapped product atoms are created. Third, "stoichiometry is defined to be 1-1 for all atoms" — if non-unit stoichiometry is needed, reactants/products must be repeated.

Partial maps are allowed: not every atom needs mapping. But the mapped portion must be strictly bijective, and unmapped atoms have deterministic creation/destruction semantics. Bond expressions must be valid SMILES (not SMARTS queries), ensuring the transformation is deterministic.

The format hierarchy is strict: **Reaction SMILES ⊂ SMIRKS ⊂ Reaction SMARTS**. Reaction SMILES has optional, possibly non-bijective atom maps representing "equivalence classes of atoms within a reaction" (Daylight's words — explicitly accommodating ambiguity). Reaction SMARTS adds query features for searching. SMIRKS occupies the middle ground: it adds SMARTS pattern-matching on the reactant side while enforcing deterministic transformation semantics.

RDKit explicitly does not claim SMIRKS support. Greg Landrum has stated: "the RDKit does not claim to support SMIRKS. It supports something it calls 'Reaction SMARTS'." RDKit's `setImplicitPropertiesFlag` toggles between SMARTS semantics (copy unspecified properties from reactant) and RXN semantics (unspecified = default) — a concrete manifestation of the formal/legacy tension.

### DPO formalism maps directly to SMIRKS

A DPO rule is a span **L ← K → R** where L is the left-hand side (pattern to match), K is the interface (preserved structure), R is the right-hand side (replacement), and `l: K → L` and `r: K → R` are injective graph morphisms. The correspondence is exact:

- **K** (interface) = SMIRKS mapped atoms — the atoms preserved across the transformation, with bijective correspondence enforced by injectivity of both morphisms
- **L \ im(l)** = unmapped reactant atoms in SMIRKS — atoms deleted by the transformation
- **R \ im(r)** = unmapped product atoms in SMIRKS — atoms created by the transformation
- **The DPO dangling condition** = the chemical sense requirement that you can't remove an atom without removing its bonds

The MØD software (Andersen, Flamm, Merkle, Stadler; University of Southern Denmark) is the primary DPO chemistry implementation, written in C++ with Python bindings. Recent work by Phan et al. (*J. Chem. Inf. Model.* 2025) in their SynTemp system proves that the Imaginary Transition State (ITS) / Condensed Graph of Reaction (CGR) representation is formally equivalent to DPO rules, unifying two previously separate traditions.

### MCES provides the bridge from unmapped data to DPO rules

The Maximum Common Edge Subgraph (MCES) approach treats atom mapping as an optimization problem: find the subgraph common to reactant and product molecular graphs that maximizes preserved bonds (minimizing chemical edit distance). This is computed by building a modular product graph and finding maximum cliques. MCES is NP-complete but tractable for typical molecule sizes. For umol, MCES would be a separate module that takes unmapped reaction SMILES and produces DPO rules with inferred atom maps.

### Category-theoretic extensions add composability

Chemical graphs form an **adhesive category**, guaranteeing well-behaved pushout complements (Behr & Sobocinski, CSL 2018). This foundation enables formal rule composition — combining two reaction steps into one — as explored extensively by Andersen et al. The **sesqui-pushout** (SqPO) extension handles "deletion in unknown context" (relevant when molecular environment is underspecified). AlgebraicRewriting.jl (Julia) is the most complete implementation of category-theoretic rewriting, using **C-sets** (functors from a schema category to Set) as the data model — essentially typed, indexed tables. A Rust implementation using `BTreeMap`-based tables indexed by schema morphisms could provide an elegant foundation. **No Rust implementation of category-theoretic graph rewriting currently exists**, representing a significant opportunity.

### The Rust cheminformatics landscape is nearly empty

The most developed Rust cheminformatics library is **ChemCore** (Richard Apodaca), which handles SMILES reading/writing and basic molecular graph operations using the `gamma` graph crate. It achieves **3-4x faster SMILES parsing** than RDKit but has no reaction handling. Other Rust efforts include the `chem` crate (282 lines, minimal) and `rdkitcffi` (thin RDKit FFI wrapper). **No Rust library handles reactions, SMIRKS, atom mapping, or graph transformations.**

---

## 5. Practical architecture for umol's reaction representation

Synthesizing across all four topics, the following architecture emerges as well-grounded in both theory and practical compatibility.

**The internal reaction rule should be a DPO span.** Define types for `ReactionRule { left: MolecularPattern, interface: AtomMap, right: MolecularTemplate }` where the interface enforces bijectivity via a `BTreeMap<MapId, (LeftAtomIdx, RightAtomIdx)>` with uniqueness on both key and value. Atoms not in the interface are explicitly tagged `Created` or `Deleted`.

**Chemical entities should be an enum, not a class hierarchy.** Rust's enum-with-data is ideal for the four Daylight categories. The expression tree for query features (following CDK's `Expr` model) should be a first-class `enum AtomExpr { And(Box<AtomExpr>, Box<AtomExpr>), Or(...), Not(...), Element(u8), Aromatic, InRing, Degree(u8), ... }`.

**Parse permissively, represent strictly.** Legacy formats (reaction SMILES, RXN files, USPTO data) should parse into a `ReactionRecord` type that tolerates partial/non-bijective/absent mappings, with explicit `ValidationResult` diagnostics. A `ReactionRecord::to_rule()` method should attempt upgrade to a strict `ReactionRule`, failing explicitly when bijectivity cannot be established.

**SMIRKS is the natural serialization format.** The mapping between SMIRKS and DPO spans is essentially isomorphic. Implement SMIRKS as the canonical string representation, with reaction SMILES and RXN as lossy import/export formats.

**Markush needs its own type.** Following ChemAxon's lead, define `MarkushStructure` as a composition of scaffold + R-group table + variation metadata, rather than trying to shoehorn it into the molecule or query types. Among the four variation types (substituent, position, frequency, homology), substituent variation (R-groups from defined sets) is by far the most common and should be the first target.

**MCES as a service module.** Since most real-world data arrives without reliable atom maps, an MCES-based atom mapping inference module is essential for practical utility. This converts unmapped `ReactionRecord` data into `ReactionRule` objects with inferred DPO spans.

## Conclusion

The fundamental insight driving umol's design is that **formal graph-theoretic semantics (DPO/SMIRKS) and real-world chemical data (USPTO/ORD/Reaxys) occupy different quality regimes**, and the library must bridge them explicitly rather than papering over the gap. The DPO span `L ← K → R` is the mathematically canonical internal form — it captures exactly the SMIRKS bijectivity constraint, has well-understood composition properties, and maps cleanly to Rust's ownership and type system. The four semantic categories (molecule, query, Markush, reaction rule) should be distinct types because no existing ontology or toolkit has successfully unified them without introducing semantic confusion. RDKit's experience with its "split personality" between reaction SMARTS and RXN semantics is a cautionary tale. Finally, the absence of any Rust library handling reactions or graph transformations means umol has the rare opportunity to build on the strongest available theory (DPO, adhesive categories) from the ground up, rather than inheriting decades of legacy compromise.
