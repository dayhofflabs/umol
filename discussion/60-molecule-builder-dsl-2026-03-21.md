# Molecule / MoleculeBuilder DSL (sketch)

## Context

GraphIR tests and tooling need a **compact, unambiguous** way to build `MoleculeBuilder` / `Molecule` — covering **all** constituent parts: atoms, covalent bonds, dative bonds, multicenter bonds, noncovalent bonds, aromatic systems, molecular charge/spin — without encoding chemistry in SMILES-style aromatic atom/bond labels as the source of truth. Aromatic systems are the first case where the pain became acute, but the DSL must cover the full model.

Resolution conformance today uses **TOML** plus atom query tokens (`umol-models-graph/tests/resolution`). That can evolve: replace TOML with a dedicated DSL, or **embed** a DSL string inside TOML (or another container) if that keeps harness plumbing stable while the molecular notation becomes explicit.

## Design principles

1. **Align with GraphIR**  
   The DSL should correspond to **`MoleculeBuilder` / `Molecule` semantics** (and shared types such as bonds, atoms, aromatic systems) so there is no second hidden model. Prefer reusing existing tokens and concepts (e.g. atom type spec/query notation) where they already match GraphIR.

2. **Declarative**  
   Prefer stating **constituent pieces** (atoms, bonds, systems, attributes). If **sequences of operations** are needed (e.g. mutating a builder), they should still be explicit in the text, not implied by side effects elsewhere.

3. **Textual**  
   **ASCII** by default. **Selected Unicode** only where a single character carries clear, stable meaning; otherwise avoid decorative or ambiguous symbols.

4. **Layout**  
   **Multi-line is acceptable** (not SMILES-style single-line pressure). **Single-line forms are a plus** for table-driven tests and embedding in other formats.

5. **Diff-friendly vs compact**  
   These goals conflict. Mitigation: **whitespace-flexible** grammar — not line-oriented, **not indentation-sensitive** — so one logical document can be pretty-printed or minified without changing meaning.

6. **Serializable**  
   A canonical or normalizable textual form is desirable (ties to diff-friendliness, tooling, and optional embedding in TOML/JSON).

7. **Human-readable and reasonably compact**  
   Favor clarity and expressivity over extreme density; avoid SMILES-style encoding of aromaticity or implicit hydrogen as the only way to name a structure.

## Non-goals (for this DSL)

- **Not** a replacement for SMILES/MOL as **external interchange** at file boundaries; those remain practical inputs where already used.
- **Not** a notation that hides **which atoms belong to which aromatic system** or **per-atom π contribution** when those matter for tests or algorithms.

## Structural layout: flat edge list vs separate sections

Two candidate structures were compared using indole (9 atoms, 10 covalent bonds, 1 aromatic system) as a concrete example.

**Flat typed edge list** — all connectivity in one stream, bond kind as prefix/tag.
Aromatic systems are annotations, not edges, so the "one flat list" premise already breaks: you get edges + annotations regardless. Mixed bond types (covalent, dative, multicenter, noncovalent) must be distinguished inline.

**Separate sections per kind** — one section per bond kind, mirroring `MoleculeBuilder` fields.
Empty sections omitted. Each section independently scannable. Diffs localize to the section that changed. Matches the data model 1:1.

**Conclusion**: the flat list converges toward separate sections in practice (aromatic annotations force a break; bond type tags approximate section headers). Separate sections have the structural advantage: direct correspondence to `Molecule` internals, cleaner diffs, omit-when-empty. **Adopt separate-sections layout.**

## Repeat / multiplicity notation

Indole exposes a verbosity problem: six atoms are all `{CHv2a1}`. Repeating the full token six times is noisy and obscures the interesting structural differences. If the DSL is for *input* (builder side), most atoms may share a common pattern.

**Idea**: allow a multiplicity prefix or suffix on atom tokens to declare N atoms of the same spec, with indices assigned sequentially. This is purely syntactic sugar — the expanded form is canonical.

**Precedents outside chemistry**:

- **Regular expressions**: `a{6}` means six repetitions of `a`. Compact, well-understood, generalizes to ranges `{3,6}`.
- **Run-length encoding** (RLE): pixel/data compression — `6×W` or `W:6`. The simplest multiplicity notation.
- **Music notation**: repeat bars `‖: ... :‖` with a count, or `×4` on a section. Also *simile marks* (repeat previous bar).
- **Hardware description (Verilog/VHDL)**: bus declarations like `wire [5:0] d` declare 6 identical signals.
- **Format strings / printf**: `%6d` — the count prefix modifies the specifier.
- **Graph DSLs (DOT/Graphviz)**: not directly, but node attributes can be set for a group: `{ rank=same; a b c d e f }`.
- **CSS / array literals**: `repeat(6, 1fr)` in CSS grid; `[0; 6]` in Rust.

A natural fit for this DSL might be `6*{CHv2a1}` or `{CHv2a1}*6`, producing atom indices 0..5 (or wherever the allocation cursor is). The `*` is visually distinct, already means "repetition of" in regex-adjacent contexts, and does not collide with any atom spec token.

**Open question**: does repeat notation extend to bonds? E.g., a linear chain `0-1 1-2 2-3 3-4` could be expressed as a chain operator. This overlaps with adjacency-inline notations and needs further thought.

## Edge notation: binary infix + hyperedge brackets

Binary edges (covalent, dative) dominate; hyperedges (aromatic systems, multicenter bonds) are rarer but structurally essential. Precedents (Datalog predicates, factor graphs, hypergraph formats) converge on: **typed bracket groups** for hyperedges, **infix operators** for binary edges.

| Arity | Syntax | Example |
|-------|--------|---------|
| 2, covalent single | `a-b` | `0-1` |
| 2, covalent double | `a=b` | `0=1` |
| 2, dative | `a->b` | `0->9` |
| n, aromatic system | `pi[...]` | `pi[0:2 1:1 2:1 3:1]` |
| n, multicenter bond | `mc[...]` | `mc[0 1 2]` |

Binary infix is syntactic sugar for 2-element typed brackets (`0-1` ≡ `bond[0 1]`), but the sugar is justified by frequency. The `prefix[members]` pattern generalizes to any new hyperedge kind without syntax changes.

## Homoiconicity and relational structure

The atom spec/query DSL already exhibits proto-homoiconicity: `{CHv2a1}` (resolved datum) and `?{CH=}` (pattern with holes) share syntax. The difference is a prefix marker and wildcard tokens. Data and patterns use the same notation.

**Generalizing**: the separate-sections layout is a set of **named relations** (in the Datalog sense). Each section name corresponds to a `MoleculeBuilder` field; each entry is a fact/tuple in that relation.

```
atoms: {NHa2} {CHa1} {CHa1} {Ca1} {CHa1} {CHa1} {CHa1} {CHa1} {Ca1}
bonds: 0-1 1-2 2-3 3-4 4-5 5-6 6-7 7-8 8-0 3-8
pi:    [0:2 1:1 2:1 3:1 4:1 5:1 6:1 7:1 8:1]
```

**Consequences of this view**:

1. **Input and assertion share notation.** Build a molecule with the DSL; assert atom properties using the same tokens. No separate assertion language.
2. **Composability = relation union.** Two molecular fragments merge by concatenating their per-section entries — directly mirroring `MoleculeBuilder`'s `Vec` fields. Benzene fragment + pyrrole fragment + fusion bond = indole.
3. **Resolution = rule application.** Input relations use query-level tokens (`?{NH=}`); output relations use spec-level tokens (`{NHv2a2}`). Same syntax, different completeness. Resolution maps patterns to facts — which is what Datalog rules express.
4. **Embeddable in s-expression-family formats.** The structural correspondence (tagged data, nesting via brackets, whitespace separation) means the DSL can round-trip through EDN, JSON, or TOML without information loss, even if we use a custom surface syntax.

**Future possibility** (not in scope now): resolution rules in the same relational vocabulary.

```
atom(?x, C, H=, _) & sigma(?x, ?v) => atom(?x, C, H(4-?v), _)
```

Designing the data format as named relations of typed tuples keeps this door open at zero cost.

## Open decisions (later)

- Grammar: declarative-only vs mixed declarative + operation trace.
- How aromatic **systems** and **contributions** are delimited and named (indices vs labels).
- Relationship to existing atom **spec** / **query** mini-languages (reuse, subset, or prefix).
- Repeat notation: prefix vs suffix (`6*X` vs `X*6`), interaction with bond/chain shorthand.
- How far to push homoiconicity: data-only format vs. pattern/query support vs. rule expressiveness.
- Surface syntax: custom grammar vs. EDN-subset vs. other host format.

---

## Surface syntax: EDN (decided 2026-03-22)

The "surface syntax" open decision is resolved: **EDN**.

The rationale is not aesthetic. The separate-sections layout is a map of named vectors —
exactly what EDN represents natively. Tagged literals (`#atom`, `#query`) handle the
existing spec/query token distinction without a custom reader. The format is
whitespace-flexible and diff-friendly by construction. A Rust EDN parser exists (`edn-rs`),
so the implementation layer (Rust core) does not require a JVM.

The key constraint was that EDN should not be shoehorned — it fits only if the molecule
representation is naturally expressible as EDN data. Inspection of the GraphIR confirms
this: `Molecule` is five parallel collections (covalent graph, dative bonds, aromatic
systems, multicenter bonds, noncovalent bonds) that map 1:1 to EDN map keys. There is no
mismatch.

**Canonical form** for a resolved molecule (ground term):

```edn
{:atoms    [#atom "{NHv2a2}" #atom "{Cv2a1H}" ...]
 :bonds    [[0 1 :single] [1 2 :single] ...]
 :dative   []
 :aromatic [{:atoms [0 1 2 3 4 5 6 7 8] :electrons 10 :charge -1}]
 :mc       []
 :nc       []}
```

**Spec vs query distinction** via tagged literals:

- `#atom "{Cv4H2}"` — ground spec (resolved `AtomTypeSpec`); maps to `Atom`
- `#query "{C?v4H*}"` — pattern with wildcards; maps to `AtomTypeQuery`

The existing `AtomTypeSpec` string notation (`{Cv4H2+}`) is the natural payload for
`#atom`. The existing wildcard tokens (`H*`, `H=`) from `HydrogenConstraint` are the
natural payload for `#query`. The tagged literal wrapper is the only addition needed.

**Consequence**: spec and query objects share syntax modulo the tag prefix and wildcard
tokens. This is the homoiconicity property, now grounded in a concrete format.

## Unified term algebra: molecules, reactions, substructures (2026-03-22)

The broader motivation for the DSL is a **unified term algebra** in which molecules,
mixtures, substructure patterns, and reaction rules are all values of the same type,
composable with functions and logic. The racemate example:

```clojure
(mix (enantiomer-1) (invert-all (enantiomer-1)))
```

is a term in this algebra. `enantiomer-1` is an EDN molecule value; `invert-all` is a
function `Molecule → Molecule`; `mix` is a constructor `[Molecule] → Mixture`. A
substructure search on the racemate expands to a search over both components — no special
casing, just function application.

This is the reason EDN (rather than a purely custom syntax) is appropriate: EDN is the
natural notation for term algebras in the Lisp/Clojure ecosystem, and the Rust core can
consume it as a data format without requiring a JVM. The computation layer (functions,
transformations, queries) sits above the data layer (EDN-serialized molecules).

Reactions under this view are `{:lhs <pattern> :rhs <pattern>}` maps; substructure
queries are patterns; molecules are ground terms. Three concepts, one format.

## Charge and spin as computed properties (2026-03-22)

Confirmed by code review of `Molecule`, `Bond`, `AromaticSystem`, `MulticenterBond`:

**Charge** is a sum over features:

```
molecule.charge = Σ atom.charge
                + Σ bond.charge
                + Σ aromatic_system.charge
                + Σ multicenter_set.charge
```

This is not a stored field — it is derived. Storing it separately would create a
consistency obligation with no benefit. The examples that motivated this design:
- HO₃⁺: charge on atom
- (Br₂)⁺: charge on the bond (`Bond.charge = 1`)
- Cp⁻: charge on the aromatic system (`AromaticSystem.charge = -1`)
- (I₃)⁻: charge on a multicenter set

All are already representable in the struct. `Molecule::charge()` (currently `todo!()`)
is arithmetic over existing fields.

**Spin** is not a sum — it requires angular momentum coupling via CG coefficients.
`SpinState::is_compatible` (in `umol-data`) implements this correctly using sequential
coupling. The molecular `SpinState` must be either provided explicitly (from input
annotation) or validated against the space of achievable couplings. It is never silently
inferred. See `discussion/61-spin-state-builder-2026-03-22.md` for the builder design
that enforces this uniformly across all feature types.

**DSL implication**: charge and multiplicity annotations in the DSL belong on the features
that carry them (bonds, aromatic systems, multicenter sets), not on a top-level
`:charge`/`:spin` key. An optional top-level `:expected {:charge 0 :multiplicity :singlet}`
serves as a checksum assertion for tests, validated at parse time, not stored in the model.

## Repeat notation and bond chains: open question refined (2026-03-22)

The open question on repeat notation (`6*{CHv2a1}`) interacts with bond notation. If
atoms 0..5 are all `{CHv2a1}` and bonds `0-1 1-2 2-3 3-4 4-5` form a chain, the bond
section is still verbose even after atom compression. A chain shorthand like `0~5`
(atoms 0 through 5 connected as a chain) would address this, but must be designed
together with the atom repeat syntax so that index allocation is unambiguous.

This is the one area where EDN's data-model neutrality works against compactness: EDN
has no built-in range or chain concept. The shorthand would be a tagged literal
(`#chain [0 5]`) or a convention in the `:bonds` vector. Either way it is purely
syntactic sugar — the canonical form is the fully expanded vector.
