# Stereochemistry encoding (ordering-independent) — 2026-01-29

This note records a design direction for stereochemistry in umol with two goals:

1. **Ordering-independent semantics** suitable for GraphIR and graph algorithms.
2. A small, typed **TableIR-level “stereo interpretation” flag** to bridge CTFile + CXSMILES
   signals into the TableIR → GraphIR conversion.

The intent is to be **principled** (graph theory + group theory) and **extensible**, without
committing to a fully general “Platonic ideal” of stereochemistry.

---

## 0. Context: current state in umol (Jan 2026)

### 0.1 TableIR

- Local stereo is stored directly on carriers:
  - `Atom.chirality: Option<Chirality>` (SMILES `@`, `@@`, `@TH`, `@AL`, `@SP`, `@TB`, `@OH`; CTFile parity maps here too)
  - `Bond.wedge: Option<BondWedge>` (SMILES `/` `\`; CTFile wedge/dash)
  - `Bond.stereo: Option<BondStereo>` (cis/trans/either, mostly via CXSMILES `c:`/`t:` right now)
- CTFile **chiral flag** (counts line `ccc`) is currently stored as a string property:
  - `molecule.properties["chiral_flag"] = "true"|"false"` (same for `ExtendedMolecule`)
- CXSMILES enhanced stereo metadata is stored in `ExtendedMolecule.cx_data: Option<CxAnnotationData>`:
  - `stereo_mode: Option<StereoMode>` (`Absolute` or `Relative`)
  - `stereo_groups: HashMap<u32, StereoSet>` where each `StereoSet` is `{ atoms, mode }` with `mode` ∈ {`Correlated`, `Independent`}

### 0.2 GraphIR

GraphIR currently carries the same local flags (`Chirality`, `BondWedge`, `BondStereo`) because
it is mostly a graph-typed container + validators. It does **not** currently provide an
ordering-independent stereochemistry model.

---

## 1. Problem statement

SMILES and CTFile encode stereochemistry using **ordering-dependent constraints**:

- SMILES `@` / `@@` (and the extended `@THn`, `@SPn`, `@TBn`, `@OHn`, `@ALn`) are interpreted
  relative to an ordered list of attachments (the “ports” around a center) as defined by the
  linear notation’s parsing conventions.
- SMILES `/` and `\` encode relative orientation on adjacent single bonds and are interpreted
  across stereogenic bond systems (double bonds, cumulenes) using traversal conventions.
- CTFile wedges/dashes are drawing-derived constraints, and their semantics require a
  convention for how 2D drawing information maps to 3D parity.

For GraphIR, reaction templates, graph algorithms, and permutation-heavy transformations
(e.g., reindexing, canonicalization, subgraph extraction/replacement), we want a representation
whose semantics is:

- **Intrinsic to the graph** (invariant under mere atom renumbering),
- **Explicit about symmetry** (so equivalences and degeneracies are principled),
- and compatible with adding stereogenic phenomena not directly representable in SMILES/CTFile
  (axis/helix/planar/spiro/atropisomerism) as *derived* stereogenic elements.

This points to a model based on:

- **Graph theory**: stereogenic features refer to a carrier substructure and its attachment
  **ports** (half-edges).
- **Group theory**: stereochemical configurations are equivalence classes of port permutations
  under a symmetry group.

---

## 2. Graph-theoretic primitive: ports (half-edges)

### 2.1 Why ports, not “neighbor atoms”

Most stereo semantics is about **attachments**, not “neighbor atoms as an unordered set”.
Two reasons:

1. A single neighbor atom can be connected by multiple bonds (multi-edges), which are distinct
   attachment sites.
2. Future GraphIR designs may include fragments/links/ports explicitly (hypergraph-like IRs).

So we model local adjacency as **ports**:

- A **port** is an incident half-edge at a node (atom-like carrier).
- For ordinary graphs, each bond endpoint is one port.

At a minimum, GraphIR needs a stable identifier for each incident half-edge:

- Center ports: `(atom_id, bond_id, endpoint_side)`
- Bond-system ports: ports on the two end atoms adjacent to a stereogenic bond system
- More complex carriers (axes/helices/planes) will define their own port sets

This lets us talk about stereochemistry as constraints on permutations of ports.

### 2.2 Ports as an overlay (keep using an undirected graph)

The port concept does **not** require a new graph data structure. It can be implemented as an
overlay on top of an ordinary undirected atom–bond graph (e.g. `petgraph`’s `StableGraph`):

- The **graph** still models atoms (nodes) and bonds (edges), and existing algorithms operate on it
  unchanged (connectivity, shortest paths, ring finding, etc.).
- A **port table** is a separate indexed structure representing half-edges. A common minimal choice:
  - one port per bond endpoint
  - `PortId` maps to `(BondId, endpoint_side)` and therefore also to an incident `(AtomId, BondId)`

In other words: ports are “half-edges”, which are a standard graph concept; we are simply
making them explicit and indexable.

### 2.3 Port tables vs port order (what “preserve ordering metadata” means)

Having a port table (indexed half-edges) is **necessary** but not always **sufficient** for decoding
format-defined stereochemistry.

SMILES `@/@@` and `/` `\` are defined relative to an **input-defined ordering** of attachments
(ports). Therefore, to interpret those markers faithfully, we need not only the set of incident
ports but also an **ordered tuple** of ports per stereogenic element as implied by the input
(or an equivalent stored permutation).

Two common ways to obtain this:

- **Explicit port order in the representation** (e.g. a PortSMILES-like syntax that enumerates
  ports in a defined order). In that case, the serialization itself satisfies the “ordering metadata”.
  Example idea: `[C:p0:p1:p2:p3]` provides an ordered port tuple.
- **Record parse-time port order** when compiling SMILES/CTFile into the port table. This can be
  stored as a “rotation list” / local port order, e.g. `atom_id -> Vec<PortId>` (or per-stereo-element).

If we do not preserve the input-defined order, we are forced to invent one later (e.g. “sort by
atom index”), which tends to reproduce RDKit-like ad-hoc behavior.

---

## 3. Group-theoretic core: stereochemical configuration as a quotient of permutations

### 3.1 Abstract definition (centers)

For a stereogenic element with \(n\) attachment ports:

- Let \(P = \{p_1,\dots,p_n\}\) be the **ports** on the carrier.
- Let \(\Pi = \{1,\dots,n\}\) be the set of **idealized positions** in the element’s geometry.
- A **labeling** is a bijection \(f: \Pi \to P\), i.e. an element of \(S_n\) once a port ordering is chosen.
- Let \(R \le S_n\) be the **proper rotational symmetry group** of the idealized geometry acting on \(\Pi\).

Then two labelings \(f, g\) represent the same physical configuration if they differ by a rotation:

\[
g = f \circ r \quad \text{for some } r \in R
\]

Thus a stereochemical configuration is an equivalence class (orbit) under \(R\), commonly
represented as a **right coset** in \(S_n / R\).

This reframes “chirality arrangements” as **group quotients**, which is exactly what OpenSMILES
does for non-tetrahedral chirality (SP/TB/OH arrangement indices).

### 3.2 The OpenSMILES central-chirality counts emerge naturally

For distinct ligands (no extra molecular symmetries), the number of inequivalent labelings is:

\[
|S_n / R| = \frac{n!}{|R|}
\]

This matches OpenSMILES’ arrangement index ranges:

| Geometry | \(n\) | Proper rotation group \(R\) | \(|R|\) | \(n!/|R|\) | OpenSMILES arr range |
|---|---:|---|---:|---:|---|
| Tetrahedral (TH) | 4 | \(A_4\) | 12 | 2 | `@`/`@@` or `@TH1..2` |
| Square planar (SP) | 4 | \(D_4\) | 8 | 3 | `@SP1..3` |
| Trigonal bipyramidal (TB) | 5 | \(D_3\) | 6 | 20 | `@TB1..20` |
| Octahedral (OH) | 6 | \(O\) | 24 | 30 | `@OH1..30` |

This is a strong argument for encoding “central chirality” as a **coset element** rather than
as an ad-hoc enumeration.

### 3.3 Molecular symmetry (identical ligands, meso, etc.)

Even if the local geometry supports multiple cosets, the molecule may not:

- If two ports are equivalent under an automorphism of the molecular graph (including chemical
  identity constraints), some “different” assignments collapse.
- Meso cases are naturally explained as the presence of nontrivial automorphisms that map one
  local configuration to another, making the global structure achiral even if local centers exist.

The key design point: **local symmetry \(R\)** is only half the story; **molecular automorphisms**
act on ports too. An ordering-independent encoding should be compatible with both.

---

## 4. Extending beyond centers: stereogenic bond systems and other carriers

Central chirality is the cleanest case. For bond-based systems (E/Z, cumulenes) and
future “non-SMILES-native” phenomena (axes/helices/planar/spiro/atropisomerism), the same
principle still applies:

- Identify a **carrier substructure** (an edge, a path, a cycle, a rigid subgraph).
- Define a set of **ports** (attachment half-edges) whose relative arrangement defines the stereochemistry.
- Define the appropriate **symmetry group** acting on the idealized positions.
- Encode the configuration as an orbit/coset.

The hard part is not representation; it is **detection** (which subgraphs are stereogenic) and
**mapping** from format-specific encodings to the port+group representation.

For the purposes of umol’s near-term scope, it is sufficient that the representation can encode:

- What SMILES/CTFile can express now:
  - central chirality (`@`, `@SP`, `@TB`, `@OH`, `@AL`)
  - stereogenic double bonds/cumulenes (`/` `\`, CTFile bond stereo markers)
- And leave room to add derived stereogenic elements later without changing the model.

---

## 5. Proposed GraphIR-level representation (conceptual)

GraphIR is the right home for “ordering-independent stereochemistry” because:

- it already represents the molecule as a graph,
- it is where semantic validation and graph algorithms live,
- and it is the natural place to define ports and port permutations.

### 5.1 Two-layer concept: raw vs canonical

It is useful to distinguish:

1. **Raw stereo constraints** captured from input formats (ordering/drawing dependent).
   - Needed for faithful roundtripping and for diagnosing weird data.
2. **Canonical/order-independent stereo semantics** used for graph operations.

These can coexist; redundancy is acceptable if invariants are clear:

- TableIR can remain a “faithful parse” container.
- GraphIR can compute and store an ordering-independent `Stereochemistry` artifact derived from TableIR.

### 5.2 Stereochemistry as a set of stereogenic elements

At the conceptual level:

- A molecule has 0..N stereogenic elements.
- Each element is described by:
  - **carrier** (where it lives in the graph),
  - **ports** (what attachments matter),
  - **geometry / symmetry group** (how permutations are quotiented),
  - **state** (which orbit/coset, or “unspecified/either”),
  - optional **group membership** (correlation / mixture semantics).

This is enough structure to:

- transport stereo through atom/bond permutations,
- reason about equivalence under automorphisms,
- and attach richer derived stereogenic elements later.

---

## 6. Conversion: TableIR → GraphIR stereochemistry

The conversion step has to bridge **format-dependent constraints** to the ordering-independent model.

### 6.1 A key fact: `@/@@` needs an ordering reference

SMILES `@`/`@@` (and the extended center arrangements) are defined relative to an ordered list
of attachments (“ports around the atom”) as determined by the SMILES parse conventions.

TableIR currently stores only the *symbol* (`Chirality::{Clockwise, CounterClockwise, ...}`),
not the *ordered port list* used to interpret it.

Therefore, an ordering-independent interpretation requires one of:

- **(A) Preserve sufficient ordering metadata from parsing**: store the input-defined ordered tuple
  of ports for the stereogenic element (or equivalently store the permutation it implies). A port
  table alone is not enough unless its indexing is meaningful for stereo (e.g., PortSMILES explicitly
  enumerates ports, or the compiler records a local port order).
- **(B) Define a deterministic ordering rule on GraphIR and interpret `@` relative to that.**
  This is viable if we can reconstruct the SMILES port order or if we accept that `@` is
  only meaningful under umol’s own ordering convention.

For a principled design, (A) is preferable because it lets the semantic layer faithfully reflect
what the input format meant, and it makes conversions well-defined.

### 6.2 “Ports + symmetry” conversion sketch for centers

For each TableIR atom with `chirality != None`:

1. Identify the **carrier atom** and construct its **port set** \(P\) (incident half-edges).
2. Determine the element’s **geometry** (TH/SP/TB/OH/AL as currently parsed).
3. Obtain an ordered list of ports \( (p_{i_1},\dots,p_{i_n}) \) corresponding to the order used by
   the input format.
4. Compute the corresponding permutation \(p \in S_n\) that maps a chosen canonical port order
   to this ordered list.
5. Combine that with the parsed arrangement indicator (`@`/`@@` or `arr`) to pick the correct
   coset/orbit in \(S_n/R\).
6. Store a stereogenic element: `(carrier, ports, geometry, coset)`.

The essential point is that the **coset** is a good, symmetry-respecting “ordering-independent”
state. The ordering only appears in the conversion step, not in the model.

### 6.3 Bond stereo (`/` `\`) conversion sketch

SMILES bond stereo is inherently relational: it constrains a *system*, not an isolated bond.

For practical purposes, GraphIR can either:

- encode a **bond stereogenic element** directly (e.g. a Z2 state for a double bond system), or
- reduce it to a set of constraints on adjacent single-bond orientations and derive E/Z later.

The group-theoretic view still helps:

- A stereogenic double bond has two “sides” (a Z2-like state) once substituents are fixed.
- Molecular symmetries can collapse E and Z (e.g., identical substituents), which should be
  detected via graph automorphisms rather than special-case rules.

---

## 7. Back to TableIR: improving CTFile chiral flag + CXSMILES relative/absolute flags

Even if the “real” stereochemistry lives in GraphIR, the TableIR → GraphIR conversion benefits
from having a typed, non-string signal for “how to interpret stereo annotations”.

### 7.1 The common role of these flags

- CTFile counts **chiral flag** is a molecule-level statement about stereochemical intent.
- CXSMILES has molecule-level `a:` (absolute) and `r` (relative) and enhanced stereo groups
  (`o<n>`, `&<n>`) describing correlation/mixture semantics.

They are not the same mechanism, but they can feed one shared concept:

> A molecule-level **stereo interpretation context** that informs how local stereo markers
> (`@`, wedges, bond stereo) should be treated in the semantic layer.

### 7.2 Minimal TableIR-level structure (proposal)

Add a small, typed field to `Molecule` and `ExtendedMolecule` (exact naming TBD):

- `stereo_interpretation: Option<StereoInterpretation>`

Where:

- `StereoInterpretation::Absolute` means “stereo annotations specify an absolute stereoisomer
  (as drawn / single-enantiomer intent)”.
- `StereoInterpretation::Relative` means “stereo annotations specify only relative
  relationships (mixture/relative intent)”.

Population rules (initial, pragmatic):

- CTFile `chiral_flag == 1` ⇒ `Some(Absolute)`
- CXSMILES `a:` ⇒ `Some(Absolute)`
- CXSMILES `r` ⇒ `Some(Relative)`

Conflict handling can be explicit in the converter (e.g., “Relative overrides Absolute” or
“record both sources and treat as diagnostic-worthy”). The important part is avoiding
stringly-typed “chiral_flag” in conversion logic.

### 7.3 What about CXSMILES enhanced stereo groups?

The existing `CxAnnotationData.stereo_groups` (`StereoSet { atoms, mode }`) is already close to a
useful typed representation. If we keep Stereochemistry in GraphIR, TableIR can still expose:

- the molecule-level interpretation (`Absolute`/`Relative`)
- optional correlation constraints on sets of stereogenic elements

Those constraints are directly useful during GraphIR stereo construction because they specify
whether certain local Z2 choices are independent or correlated (and later, can generalize beyond Z2).

---

## 8. Design considerations and open questions

### 8.1 Canonicalization vs equivariance

There are two different “ordering independence” goals:

- **Equivariance**: if you permute atom indices and apply the same permutation to stereo data,
  the meaning is unchanged. This is what reaction templates and graph rewriting need.
- **Canonical invariance**: a unique representation for an unlabeled graph (canonical SMILES-like).

The proposed port+group model supports equivariance naturally. Canonical invariance requires a
canonical labeling algorithm (future).

### 8.2 Detection vs representation

Representing “helix chirality” or “atropisomerism” is not the hard part once ports are available.
The hard part is defining:

- which subgraphs are stereogenic under a chosen chemical model,
- what the ports are for that element,
- and how to handle flexibility/fluxionality (chemically real but format-invisible).

Those belong in the semantic/analysis layer, not in the parser.

### 8.3 Relation to CIP (R/S, E/Z naming)

The encoding here is about **stereochemical state** as a graph+group artifact, not about CIP naming.
CIP assignment is a derived, convention-heavy labeling that can be layered on later (and can be
made rigorous once the underlying state representation is sound).

---

## 9. Next step (future implementation, not in this note)

1. Introduce a typed `StereoInterpretation` field (TableIR) replacing string `"chiral_flag"` for
   conversion purposes.
2. Design a `GraphIR::Stereochemistry` container around stereogenic elements defined by ports and
   symmetry groups.
3. Extend `sir_to_gir` to **construct** this stereochemistry object (at least for TH/SP/TB/OH and E/Z),
   using whatever ordering metadata is available (spans, explicit port order, or both).

