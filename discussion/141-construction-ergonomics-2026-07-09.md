# 141 · Molecule construction ergonomics — survey and directions

Status: Active (exploratory design record)
Date: 2026-07-09
Relates: 060 (molecule-builder-dsl — the string path), 086 (edit API: join/split/`MoleculeParts`),
103 (ports/fragments — the operad substrate), 042 (relational representation), 140 (Python
entity bindings — eventual mirror). Prior inspiration: 008/017/019 (early builders), 009 (torch),
010 (diesel).

## Motivation

Programmatic, in-code molecule construction — build a molecule atom-by-atom / fragment-by-fragment
in code, not by parsing a SMILES/MOL string. "Everyone just reads a string" is a *symptom* of the
in-code builders being unpleasant, not evidence that construction doesn't matter. This doc surveys
what the field does, what the prior umol exploration missed, and unconventional directions; it is a
think-through before committing to the conventional path.

## Survey — `./materials/codes`

Three tiers on the reference-model axis (how a bond names its two atoms):

- **Index juggling (the norm): RDKit, OpenBabel, LillyMol.** `idx = AddAtom(...)`, then
  `AddBond(i, j)`. OpenBabel's is the classic footgun — `NewAtom()` returns an object but `AddBond`
  wants an integer, so you round-trip through `GetIdx()`.
- **String-only: MØD.** No atom-by-atom builder *at all* — construction is a GML/SMILES/DFS string.
  The most sophisticated reaction engine surveyed deliberately has no programmatic path.
- **Atom handles (the elegant minority):**
  - **mx** (Apodaca): `a = mol.addAtom("C"); mol.connect(a, b, 1)` — `connect`/`disconnect` verbs.
  - **Indigo**: `a.addBond(b, order)` — the bond verb lives *on the atom*.
  - **octet**: `AtomProxy` handles; `connectSingle(a, b)` — handles opaque, so an index can't leak.
  - **CDK**: `newBond(atomA, atomB, order)`.

The elegant answer is the ORM insight: *the reference is the object add-atom hands back, never an
index.* umol's Rust side already is this (`add_atom -> AtomId`); the pain was only the flat
constructors (`(AtomId(0), AtomId(1), …)` literals), now retired by `MoleculeParts` (086).

## Prior umol exploration was thin on the reference problem

- **009 (torch)** — API *organization* (`torch.nn` namespacing) + composition; explicitly "don't
  worry about the molecular example." Never touches atom-reference.
- **010 (diesel)** — type-level *query* composition (fluent chaining, wrapper types). A read DSL, not
  construction.
- **060** — put the answer in the *string* DSL (named-atom EDN labels). **019** — active-atom cursor,
  superseded. So the two "unconventional source" attempts chased fluency/organization, not the core
  ergonomic. The thin spot the user remembered is real.

## Unconventional sources — the distinct idea each brings

Handles are the ceiling of chemistry-world thinking. Above it:

| Source | Idea | What it offers a molecule builder |
|---|---|---|
| Graph DBs (Cypher) | the code is *isomorphic* to the graph — `CREATE (c:C)-[:B]->(o:O)` | a **visual literal**: `mol! { C=c -o- O }` reads like the structure, not like assembly steps |
| networkx | node key is a *meaningful value*, not a position — `add_edge("carbonyl_C", "amide_N")` | **named keys** resolved at build (`HashMap<Label, AtomId>`) — 060's names, in code, no string parse |
| ggplot2 / grammar of graphics | composition by `+`; each term self-contained; the spec is a value | **`+`-algebra of relation-bearing terms**; maps onto 042 (relation union) and the 086 `join` |
| CAD (build123d / OpenSCAD) | implicit build context + cursor + *selectors* (`select(aromatic)`) | **relative + predicate attachment** — the chain case where naming every atom is overkill (= SMILES chains) |
| Datalog | you *assert* relations, you don't build — `atom(c, carbon). bond(c, o).` | the declarative form of 042's relational model |
| parser combinators | fragments as first-class *composable values* — `benzene.fuse(pyridine, …)` | **fragment-first** assembly (the common real case) |
| octet / Apodaca | bonds as first-class verbs / *bonding systems* — `connectDouble` | **per-family construction verbs** (`single`, `aromatic_ring`, `dative`, `pi_system`) matching the de-privileged entities |

## The design axes (orthogonal — combine, don't pick one)

- **Reference model** — handle (AtomId/local var) · named key (label→id) · cursor ("current atom") ·
  none-needed (fragment algebra).
- **Composition** — imperative add-calls · `+`/`fuse`/`join` algebra over fragments (042 + 086).
- **Notation** — method calls · `+`-terms · a `mol!` *visual* macro (Cypher/SMILES-like, compile-time).

## The layered proposal

1. **Handle base** — `let c = b.atom(C); b.single(c, o)` (mx/Indigo; already have it). Per-family
   verbs (octet) rather than a generic `add_bond(kind)`.
2. **Fragment algebra** — fragments/molecules as values composed with `+` / `fuse` / `join`
   (ggplot + 042 + 086). *The biggest real win*: assembly-from-known-pieces is what people actually
   do, and it falls out of the relational core and the `join` already designed.
3. **A `mol!` visual macro** — Cypher/SMILES-like literals, compile-time, desugaring to the handle
   API — the in-code twin of 060 with no string parse, type-checked. *The second real win*: it's the
   only thing that beats a SMILES string on its own terms while staying in checked code.
4. **Optional cursor / selector sugar** for chains — only if the chain case earns it.

The rebellious conclusion: don't build yet-another-fluent-builder (the diesel/torch reflex). The
handle builder is fine and mostly done. Put effort where the field *hasn't* — **fragment-first
assembly** and the **`mol!` visual macro**. Those are where umol could be genuinely better, not just
less-bad than RDKit.

### `mol!` freed for visual literals

`mol!` today is the macro that reads the molecule-DSL string; **rename it** (candidate: `dsl!` /
`mol_dsl!` / `parse_mol!` — nomenclature TBD) to reserve `mol!` for the visual-literal macro above.

## Fragments as typed operads (the deep direction)

The fragment-first path is not just an ergonomic — it is a structural model worth taking seriously,
and it is the *construction side* of the port model in 103.

- A **fragment** is an operation with a **typed, ordered port interface** — its "arity" is its free
  ports, each typed by what it can bond to (valence, bond order, dative donor/acceptor, …). This is
  exactly 103's "fragment = subgraph with an ordered, typed port interface."
- Fragments form a **typed (coloured) operad**: the colours are **port types**, and **operadic
  composition** is *plugging one fragment's output port into another's input port* (γ), with the
  symmetric-group action = reordering ports. Building a molecule is composing fragment-operations
  until no free ports remain (a closed molecule) or some remain (a building block / open fragment).
- This unifies four things under one algebra: **construction** (compose fragments), the **port/entity
  model** (103 — fragments as primitives, atom|port|fragment ligand sum), **reactions** (a rule is an
  operation on fragments — the reaction-network goal), and **combinatorial libraries** (enumerate by
  composing a fragment set — 100k–1B nodes).
- Status: ambitious / research, gated like 103's `umol-port` spike. The ergonomic layers above stand
  on their own; the operad is the "do it right" target the spike would validate. Prototype the
  fragment-algebra signatures first (they are useful with or without the full operad), then see
  whether the port/operad substrate earns the base-graph rewrite.

## Open questions (engaged)

### (a) Are `AtomId`s sufficient as handles, or do we need named keys / aliases?

`AtomId` (u32 newtype) is the right **runtime handle** — it *is* the atom's identity, it's what
`add_atom` returns, it indexes storage directly. Keep it. But three *distinct* needs get conflated,
and only the first is a handle:

1. **Handle (identity of a placed atom)** → `AtomId`, bound to `let c = …`. Sufficient. (Caveat: an
   `AtomId` goes stale under a `Compaction`; during pure construction there's no removal, so it's
   stable. Edit-during-construction is the transactional layer's `New(n)`/keys, not a new id type.)
2. **Named key (readable cross-reference)** → a *build-scoped* label map (`HashMap<Label, AtomId>`),
   resolved at `finish`. Lets you write `bond("carbonyl_C", "amide_N")` — 060's names, in code. A
   construction-layer convenience, **not** a change to `AtomId`.
3. **Alias / template (a reusable atom *spec*, "simplifying atom definitions")** → `AtomAst` is
   already a `Clone` value, so a template is just a reused `AtomAst` (`let ar = AtomAst(C, aromatic);`
   place it six times). What's missing is a terse instantiation vocabulary and a small template
   library (common groups) — sugar, not a type.

The fragment/operad lift refines this: **within** a fragment you reference `AtomId`s (or fragment-
local names); **across** fragments the handle is a **port** (the typed interface), not a raw
`AtomId`. So `AtomId` is sufficient inside a fragment; the cross-fragment handle is a port. No new
atom-id type is needed — the additions are a build-time name map (2) and templates (3), plus ports
(from 103) as the composition-level handle.

### (d) Specifying relations for `+` efficiently

The trap is `+`-ing individual bonds: `mol + bond(a,b) + bond(b,c) + …` is verbose. The fix: `+`
layers **higher-level relation-bearing terms**, each internally specifying many relations compactly:

- `chain([a, b, c, d])` — a path of bonds (the cursor/SMILES idea, as a term).
- `ring([…])` — a cycle of bonds + closure.
- a **fragment / template** — a whole pre-wired internal structure in one term.
- per-family verbs — `aromatic([…])`, `dative(donors, acceptor)`, `pi_system([…])`.

So you rarely write bond-by-bond; you compose *terms*, and the efficiency is in the **term
vocabulary**, not a clever per-bond syntax. This is `+` = relation union (042) / operadic composition
(above); the terms are the operad's operations. When topology is the clearest expression, the `mol!`
visual macro draws it directly.

## Sketch — the four layers

The stack (details in the 2026-07-09 conversation): **L1 handle builder + per-family verbs**
(`b.atom(C) -> AtomId`; `single`/`double`/`triple`/`dative`; `chain`/`ring`; the overlay verbs
`aromatic_system`/`multicenter`/`noncovalent` — see *Implemented — L1*); **L2 `+`-spec** (layer
relation-*terms* — `atoms`/`chain`/`ring`/`dative`/overlay terms — onto one molecule, refs by name,
`+` = relation union); **L3 fragment/operad** (`Fragment { body, ports }`, `attach(port, other,
port)` = operadic γ, `+` = disjoint union, `close() -> MoleculeAst`, a fragment library); **L4
`mol!` visual macro** (Cypher/DOT/SMILES-flavored, desugars to L1/L2). Each lowers to the one below.

### Settled decisions (2026-07-09)

1. **`Fragment.body = MoleculeAst`** (reuse, simpler) — not a distinct lighter type.
2. **Ports now = attachment points** (`AtomId` + interface name + the `BondAst` formed on attach),
   lowering to real typed half-edges when the port model (103) lands — signatures unchanged.
3. **`+` is dual-use** — L2 relation-layering vs L3 disjoint union — disambiguated by operand type;
   fine, no rename.
4. **`mol!` = a `#[proc_macro]`** in `umol-ast-macros` (which already hosts the `Lattice`/
   `Canonicalize` derive proc-macros). The current `macro_rules!` DSL-string readers in
   `umol-ast/src/macros.rs` get a **`_dsl` suffix** to free the short names. Rename scheme (21
   macros, ~180 call sites — the first step of the construction slice):
   - **Remove** the metadata `mol_dsl!` (→ `MoleculeDsl`); its few sites move to
     `MoleculeDsl::from_str`.
   - Then suffix each reader: `mol→mol_dsl`, `mol_ground→mol_dsl_ground`, `atom→atom_dsl`,
     `atom_ground→atom_dsl_ground`, `partial_atom→partial_atom_dsl`, `bond→bond_dsl`,
     `bond_ground→bond_dsl_ground`, `partial_bond→partial_bond_dsl`, `dative→dative_dsl`
     (+`_ground`), `aromatic→aromatic_dsl` (+`_ground`), `multicenter→multicenter_dsl` (+`_ground`),
     `noncovalent→noncovalent_dsl` (+`_ground`), `stereo_atom→stereo_atom_dsl` (+`_ground`),
     `stereo_bond→stereo_bond_dsl` (+`_ground`). (`_ground` sits after `_dsl`.)
5. **`close()` leaves free ports as open valence**; resolution is a separate step. **Open:** an
   ergonomic way to assert **charges / unpaired electrons** at construction (not punted to
   resolution).
6. **Per-family bond verbs** (`single`/`double`/`triple`/`dative`), not a generic `bond(kind)`.
   **`dative` is not a bond subtype** — a separate family. Overlays: at **L1** they get primitive
   construction verbs (`aromatic_system`/`multicenter`/`noncovalent` — see *Implemented — L1*), since
   electron counts are usually known at build. **Still open** is how overlays participate in
   **fragments/ports** (L3) — fragments are atoms/bonds only for now; overlays-in-fragments needs the
   operad/port approach. (Superseded: `aromatic_ring` — dropped, see *Implemented — L1*.)

Note: the DSL already carries `:atom-aliases` (reusable atom specs) — the in-code analog of the
"template" answer to question (a).

### Sketch-pass resolutions (2026-07-09)

Sketched L1–L3 as a standalone compiling file (`todo!()` bodies); it typechecks. Four things it
settled:

- **Error handling (sharpens decision 5).** `mol!` (a proc-macro that sees every token) uses
  **`compile_error!`** for everything static — bad syntax, undeclared ident cross-refs, and
  statically-incompatible ports — the diesel lesson (push validation to compile time). The
  **runtime builder** (L2 names from a runtime `&str`, L3 `attach`) **panics** — a wrong name is a
  code bug, like out-of-bounds indexing, and it matches the existing `mol!(…).unwrap()` convention.
  `Result` only for construction from *external/untrusted* data.
- **L2 term naming.** Free-fn terms (`bond`, `dative`, `chain`, `ring`, `aromatic`) do **not** clash
  with the builder *methods* (`b.single`) or the `*_dsl!` macros — different namespaces (`bond` fn
  and `bond!` macro coexist). No `_term` suffix; a `build` prelude is an import convenience, not a
  necessity.
- **Heterogeneous atom specs.** The `atoms([…])` slice term stays **homogeneous** (generic
  `Into<AtomAst>`, no per-element `.into()`). **Mixed** specs (`Element` + `"C#h3"`) go through the
  **`mol!` macro**, where each atom decl is an independent expr — cleaner than eager `.into()`s.
- **Port colours = the bond lattice (sharpens decision 6 / the operad typing).** A `Port` carries a
  possibly-`Undetermined` `BondAst`; `attach` = **`BondAst::meet`** of the two specs (`BondAst`
  already impls `Lattice`, so it's free). An `Undetermined` (⊤) port **absorbs** the partner —
  `meet(Undetermined, double) = double` — so a port can be left order-open and inherit it. `meet =
  None` (⊥) → error (per the error rule above). Crucially this is the **same `meet` as
  `meet_pushout`** (reaction gluing): fragment-`attach` and reaction-gluing are one lattice-meet
  composition at different scales, and it generalizes to overlay-coloured ports via each family's
  `Lattice` when overlays enter fragments. The operad typing falls out of existing machinery.

Compiling reference sketch: `scratchpad/construction_sketch.rs`.

### Open items resolved (2026-07-09)

- **Charge / unpaired / multiplicity = the DSL defaults, reused.** This is `raise`'s job, so reuse
  `AtomDefaults` (and the per-family defaults) rather than inventing a construction-only mechanism.
  The builder/spec holds **partial** atoms + an `AtomDefaults`; **`build()` applies them** — exactly
  what DSL `raise` does, one machinery for both paths. Applied at `build` (not at `atom()`), so it's
  **order-independent** (a `+ charge(0)` term anywhere sets the molecule-wide default, per-atom specs
  override). Exposed at both layers: L1 `MoleculeBuilder::ground()` / `with_defaults(…)` /
  `default_charge(0)`; L2 `+ ground()` / `+ charge(0)` / `+ spin(singlet)` (desugars to the L1
  default). Prefer the named presets (`ground()` = neutral singlet, `zeroed()`) over spelling out
  `charge(0)+unpaired(0)+multiplicity(1)`; keep individual terms for partial control.
- **Overlays in fragments = n-ary port-generators.** An overlay participates by referencing **ports**
  (like a bond); a fragment can carry a **partially-formed overlay** whose participants are still free
  ports, closed on `attach`. Same port-wiring generalized to n-ary — colours are each family's
  `Lattice`, composition is the family `meet` (the identical primitive as bonds and `meet_pushout`).
  This folds directly into 103's open "can fragments subsume the overlays as primitives?" — in this
  algebra overlays and fragments are both boxes-with-ports; which one is primitive is the
  `umol-port` spike's call, not a new question.

### Underlying structure: a typed hypergraph

The compositional model is a **typed hypergraph** — Spivak's **operad of wiring diagrams** (the
cycle-allowing kind; algebras are **hypergraph categories**), *not* the tree operad. This is the
categorical semantics of the whole relation model, not just construction:

- **ports** = objects / colours;
- **bonds and overlays** = typed hyperedges — n-ary generators, arity + port-typing per family (bond
  = 2 symmetric legs, dative = n donor + 1 acceptor, aromatic/multicenter = n legs, stereo = site +
  ordered ligands);
- **composition** = port-wiring with per-family lattice-`meet` (⊤/`Undetermined` ports absorb the
  partner; ⊥ = incompatible → error);
- **disjoint union** = the monoidal product (086 `join`);
- **ring closure** = cyclic wiring — the move past a tree operad.

Construction (`attach`), reaction gluing (`meet_pushout`, 086), and the relation storage (042
relational core, 103 ports/fragments) are three faces of this one algebra. `ring_close` in the
sketch is the tell that we're in the wiring-diagram operad, not the tree one.

## Implemented — L1 `MoleculeBuilder` (2026-07-09)

`umol-ast/src/ast/molecule/build.rs`. Bare-verb-adds convention (every method adds/declares; lookup
is `MoleculeEditor`), reached via `MoleculeAst::builder()`, wrapping a `MoleculeEditor`. Surface:

| group | verbs |
|---|---|
| defaults context | `new()` · `ground()` |
| atoms | `atom(impl Into<AtomAst>) -> AtomId` — element, `AtomAst`, or `"C#h3"` string |
| bonds (single-entity) | `single` · `double` · `triple` · `dative(donors, acceptor)` |
| bond composites | `chain(atoms) -> Vec<BondId>` (path of singles) · `ring(atoms)` (+ closing single) |
| overlays (primitive) | `aromatic_system(atoms, electrons)` · `multicenter(atoms, electrons)` · `noncovalent(first, second, kind)` |
| finalize | `build() -> MoleculeAst` |

Three shifts vs. the original L1 sketch:

1. **Overlays are constructible at build.** Decision 6's "overlays open" was about overlays *in
   fragments/ports* (L3, still open); at L1 the overlay families get primitive verbs directly, since
   electron counts are usually known when building the molecule (benzene `[1;6]`, cyclopentadienyl
   `[1;5]`). They lower onto the editor's `add_aromatic_system` / `add_multicenter_bond` /
   `add_noncovalent_bond` via `from_electrons` / `from_kind`.
2. **`aromatic_ring` dropped.** It conflated three things (ring σ-bonds + per-atom aromatic flags +
   the overlay). Aromaticity is now expressed by the *atoms* passed (`atom("C#a")` →
   `AromaticValence(Aromatic)`) or by resolution; the ring is `ring(atoms)`; the overlay, when known,
   is `aromatic_system(atoms, electrons)`.
3. **`chain`/`ring` at L1** as the handle-wiring primitive (fold over `single`). The L2 `+`-spec
   `chain`/`ring` terms will be **create-atoms-then-delegate** wrappers over these — one wiring
   implementation, layered, so "chain" is not two parallel copies. The operations differ only by
   input: L1 wires *existing handles*; L2 *creates* atoms then wires.

Defaults: `ground()` applies `AtomAst::into_ground` per atom, preserving preset fields (a `+2` charge
survives; an unspecified field grounds to its default). Apply-at-add is order-independent here because
the mode is fixed at construction; the L2 `+ ground()`/`+ charge(0)` terms are the apply-at-build
case (order-independent across a spec).

Endpoints: the graph stores bond endpoints **normalized (min, max)** (undirected), so `ring`'s
closing bond `single(last, first)` reads back as `atom_ids() == [first, last]`.

## L2 `MoleculeSpec` design (settled 2026-07-09)

A `MoleculeSpec` is a value; free-fn **terms** return a `MoleculeSpecTerm` (name scoped for the coming
`ReactionSpec`); `spec + term -> spec`; `build() -> MoleculeAst`. Module `ast/molecule/spec.rs`
(sibling to `build.rs`); terms lower onto the L1 `MoleculeBuilder`.

### `AtomArg` — create-or-reference, unified

Every atom slot in every term is an `AtomArg`. What you write picks create vs reference by type — no
`&str` collision, so D1/D2/D3 are one mechanism and *every* term create-or-wires:

| write | means |
|---|---|
| `C` · `"C#h3"` · `AtomAst` | create a fresh **anonymous** atom |
| `("carbonyl", C)` tuple | create a fresh atom, **named** |
| `0` (bare integer) | reference existing atom by **position** (creation order; `From<i32>`+`From<u32>`, negative panics) |
| `name("carbonyl")` | reference existing atom by **name** |

Bare `&str` is a spec (matches L1 `atom("C")`); a by-name reference is the `name(…)` wrapper. So
`single(C, O)` mints a diatomic, `single(0, 1)` wires two existing, and
`ring([name("a"), name("b"), C, C])` closes a fused ring over two existing + two new atoms.

### Terms

| term | signature | lowers to |
|---|---|---|
| `atom` | `atom(impl Into<AtomArg>)` | introduce |
| `atoms` | `atoms(impl IntoIterator<Item = impl Into<AtomArg>>)` | introduce ×N |
| `single`/`double`/`triple` | `(AtomArg, AtomArg)` | L1 `single`/… |
| `bond` | `(AtomArg, AtomArg, impl Into<BondAst>)` | L1 `bond` (charge/spin/constraints) |
| `aromatic_bond` | `(AtomArg, AtomArg)` | L1 `aromatic_bond` (order 1 + `#a`) |
| `chain` | `(impl IntoIterator<Item = AtomArg>)` | L1 `chain` |
| `ring` | `(impl IntoIterator<Item = AtomArg>)` | L1 `ring` |
| `dative_bond` | `(donors, acceptor)` | L1 `dative_bond` |
| `aromatic_system` | `(atoms, electrons)` | L1 `aromatic_system` |
| `multicenter_bond` | `(atoms, electrons)` | L1 `multicenter_bond` |
| `noncovalent_bond` | `(AtomArg, AtomArg, NoncovalentBondKind)` | L1 `noncovalent_bond` |
| `ground` | `()` | build-time field fill |

Deferred (as at L1): `stereo_atom`/`stereo_bond` (stereo slice); molecule-level constraint terms.

### Decisions

- **Relations spelled out in full** — `dative_bond`, `multicenter_bond`, `noncovalent_bond`,
  `stereo_atom`, `stereo_bond`, `aromatic_system`, `aromatic_bond`; `single`/`double`/`triple` are
  bond *orders* (noun implied), `bond`/`chain`/`ring`/`atom(s)` generic/topology. Applied to L1 too —
  the earlier `dative`/`multicenter`/`noncovalent` were renamed.
- **`aromatic_bond` = order 1 + aromatic flag** (`1#a`), not undetermined order; resolution perceives
  the system. Not exclusive with `aromatic_system`.
- **`bond` exists at L1 and L2** — the only way to set a bond's `charge`/`spin`/constraints (which
  `BondAst` carries and the order verbs don't).
- **Named-create is a `(name, spec)` tuple**, not a `named(…)` verb.
- **Defaults = `ground()` only.** `zeroed` dropped (going away); no arbitrary-value field terms
  (`charge(-2)` is meaningless) — a field default is its canonical `0`/`1` or nothing. Per-field
  nullary terms omitted for now; revisit if granularity is needed.
- **Build order** fixed regardless of `+`-order: introduce (create + build the name/position map) →
  wire (resolve `AtomArg`s) → defaults (`ground` fill). So `+ ground()` anywhere is order-independent;
  positions are creation-order (the `+`-order of introduce terms).
- **`+` dual-use** with L3 disjoint union stays disambiguated by operand type (decision 3).

## L3 `Fragment` design + plan (2026-07-10) — unblocked (join/split done, 142)

`Fragment { body: MoleculeAst, ports: Vec<Port> }`; `Port { atom: AtomId, bond: BondAst, name:
Option<String> }`.

**Ports are a positional, colour-typed interface** — the operad arity: an ordered `Vec<Port>`, each
typed by its `bond: BondAst` (the colour). **Names are convenience, not typing** — an optional label to
address a port, exactly like the DSL's atom-aliases; compatibility is `meet` of the two bonds, never the
name. A port is addressed by `PortArg = Index(u32) | Name(String)` (`From<u32>`/`From<i32>` → index,
`From<&str>` → name — the `AtomArg` shape). Name lookup panics if not found or ambiguous (use the index).

Operations:
- **`attach(self, self_port, other, other_port) -> Fragment`** (operadic γ — the *only* thing that
  wires) = `join` the two bodies, add a bond between the two addressed port atoms with `BondAst =
  meet(p.bond, q.bond)`, drop the two consumed ports; remaining ports carry forward, `other`'s remapped
  through the join correspondence. Panics on a bad port ref or a ⊥-`meet` (both construction bugs, per
  the runtime-builder-panics convention).
- **`Fragment + Fragment -> Fragment`** (monoidal ⊗ — juxtaposition, forms **no** bond) = `join` alone,
  ports concatenated (`other`'s remapped). No port correlation is guessed. For genuinely multi-component
  structures — salts, mixtures, the two reactants of a bimolecular reaction (its main use). `+`
  disambiguated from L2 by operand type.
- **`close(self) -> MoleculeAst`** = the body; free ports become open valence (resolution fills later).

Construction: `Fragment::new(body)` + fluent infallible `with_port(name, atom, bond)` and an unnamed
variant. Overlay-as-port and a fragment library are deferred (the latter needs the 103 port spike).

### Impl plan — done (2026-07-10), `ast/molecule/fragment.rs`

- **S1** — `Port`, `PortArg` (+`From` impls), `Fragment` struct, `new`, `with_port` / unnamed variant,
  `close`. Tests: build a fragment with ports; `close` returns the body. Additive.
- **S2** — port-remap helper (map a port's atom through a join `MoleculeCorrespondence` —
  `corr.atoms().right_of`) + `impl Add<Fragment> for Fragment` (disjoint union; ports concatenated,
  `other`'s remapped). Tests: `a + b` body = `join`, all ports carried, `other`'s atoms remapped.
  `[dep: S1]` Additive.
- **S3** — `attach(self, impl Into<PortArg>, other, impl Into<PortArg>) -> Fragment` (join + `meet`-bond
  + drop the two ports, reusing the S2 remap helper). Tests: attach forms the `meet` bond, consumes the
  two ports, carries the rest; `#[should_panic]` for a missing port ref and for a ⊥-`meet`. `[dep: S2]`
  Additive.
- Register `mod fragment` + re-export `Fragment` / `Port` / `PortArg` at `molecule` and `ast`.

S1 (`Port`/`PortArg`/`Fragment`/`new`/`with_port`/`close`), S2 (`+` disjoint union + port remap), S3
(`attach` = join + `meet`-bond + drop ports) all landed; 10 tests. Nomenclature aligned with L2:
`AtomArg` variants `New`/`Index`/`Name`, `PortArg` variants `Index`/`Name`. Deferred (last): fragment
library, overlay-as-port (gated on the 103 spike).

## L4 `mol!` / `frag!` — visual-literal macros (design + plan, 2026-07-10)

Compile-checked **visual surface over L2**: parse a path-pattern literal, validate its structure at
expansion (undeclared reference → `compile_error!`), emit `MoleculeSpec` code. `mol!` → `MoleculeAst`;
`frag!` → `Fragment`. **Not a serialization/traversal** — a set of *local path-patterns joined by shared
names* (Cypher/Datalog), so there is no frontier or backedge machinery: names are the only long-range
reference (the traversal-notation analysis' "sparse long-range reference", done ergonomically — see
`materials/graph_traversal/`). Write procedurally, desugar declaratively. **Homoiconic for free**: the
same notation yields a ground molecule or a pattern — whichever the specs are (`*` / `Undetermined` /
constraints ⇒ pattern), since both are just `MoleculeAst`.

### Grammar (must lex as Rust tokens)

- **Atoms** `(name: elem)` — `elem` is a **quoteless element** (`C`, `Na`) for the common case, or a
  **quoted DSL spec** `"C#h3"` / `"*"` for anything with `#`, whitespace, constraints, or the
  undetermined-atom wildcard. First mention declares; bare `(name)` references (the non-tree edges).
- **Bonds** — order shortcuts `-` / `=` / `#`; general **Cypher-style** `-[ "spec" ]-` (optionally
  `-[name: "spec"]-`), spec = a DSL bond spec (order, `#a` aromatic, charge, spin, ring); the optional
  name enables bond references. Atoms `(…)` and bonds `-[…]-` are symmetric (shortcut or quoted spec).
- **Ports** (`frag!` only) — `*name` in an atom slot: a bond to `*name` declares a port on the other
  endpoint, colour read off the bond op (SMILES `*` / R-group / operad boundary). No separate syntax.
- **Overlays** — keyword statements in the same block (**not** nested macros), referencing named atoms:
  `aromatic [a, b, …] : "electrons"`, `dative [donors…] -> acceptor`, `multicenter [ … ] : …`,
  `noncovalent (a) ~ (b) : "kind"`, `stereo (site) [ligands…] : "kind"`. Each desugars to its L2 term.
- Statements comma/newline separated.

### Desugaring → L2

Atom decl → `atoms`/`atom` term (`(name, spec)` create); ref → `name(…)`; bond → `single`/`double`/
`triple`/`bond`; overlay → its term. `mol!` emits `(MoleculeSpec::new() + … ).build()`; `frag!` wraps in
`Fragment::new(…).with_port(…)`.

### Port finishing (refines L3 `close`)

The ground-vs-pattern duality of a fragment's free ports:
- **`close()`** — *ground*: every port must be paired; a free port is an error.
- **pattern close** (name TBD) — each free port becomes an undetermined `*` `AtomAst` bonded via the
  port's colour ("something attaches here").
(Supersedes decision 5's vaguer "open valence"; a lenient-ground "leave the valence open for resolution"
could be a third mode if wanted.)

### Compile-time checking — **decided: A** (2026-07-10)

The DSL parser is bound to `AtomAst` → `umol-ast` → `umol-ast-macros` (and `AtomAst` itself derives
`Lattice`/`Canonicalize` from `umol-ast-macros`), so the parser can never sit below `umol-ast-macros`;
relocating it doesn't help. The lever is `mol!`'s crate, not the parser's.

- **Structure / refs** — always validated at expansion → `compile_error!` (undeclared ref, duplicate
  name, malformed path). The primary "checked code" win, present in every option.
- **Bare elements** — `umol-chem` does not depend on `umol-ast-macros`, so bare `(c: C)` is
  compile-checked regardless (an invalid element makes the emitted `Element::…` path fail to compile; a
  friendlier macro-span error via a `umol-chem` dep is an easy refinement).
- **Rich quoted specs** (`"C#h3"`) — **A (chosen)**: `mol!`/`frag!` live in `umol-ast-macros`, emit
  `AtomAst::from(&str)`, so a bad spec panics at first runtime build. **B** (a separate `umol-mol-macro`
  crate depending on `umol-ast`) or **C** (an AST-free grammar crate shared with the real parser) add
  compile-time spec checking later; both are additive — same syntax, same codegen, only the validation
  step changes — so no rework.

### Impl plan (`umol-ast-macros`, `#[proc_macro]`)

- **S1 — done** (`umol-ast-macros/src/mol.rs`). `syn` parser for the path grammar (atoms
  `(name: elem | "spec")`, bare `(name)` refs, bonds `-`/`=`/`#`); structure validation → `compile_error!`
  (undeclared ref / duplicate decl, precise span); codegen emits `(MoleculeSpec::new() + atoms([…]) +
  … ).build()` via `::umol_ast::ast::spec::*` (atoms as strings, so bare `C` and `"C#h3"` both flow —
  runtime-parsed per decision A); `mol!` → `MoleculeAst`, re-exported at `umol_ast::mol`. Tests
  (`umol-ast/tests/`): three runtime builds (orders, quoted spec, `-[ "1#a" ]-` rich bond) + a `trybuild`
  compile-fail for an undeclared ref. Grammar complete for S1 (incl. `-[ "spec" ]-` → L2 `bond(…)`).
  Alongside: **`From<&str>` added to every entity AST** (bond/dative/aromatic/multicenter/noncovalent/
  stereo-atom/stereo-bond, mirroring `AtomAst`; stereo via the `stereo_element!` macro) so a quoted DSL
  spec becomes the entity via `Into` — the uniform string-spec surface the macro emits onto.
- **S2** — `frag!` + `*name` ports → `Fragment`; the two `close` modes (ground-error / pattern-`*`).
  `[dep: S1]`
- **S3** — overlay keyword statements (aromatic / dative / multicenter / noncovalent / stereo). `[dep: S1]`
- Deferrable: bond-spec compile-validation (the shared-parser-crate split).

`mol!` and `frag!` are split (revisit if they converge). Quoteless elements accepted; anything richer is
a quoted DSL spec.

## Next

- **L1 + L2 — done** (`build.rs`, `spec.rs`). See the *Implemented — L1* and *L2 `MoleculeSpec`*
  sections above.
- **Detour — join/split + connectivity validator — done** (doc 142). Unblocked L3.
- **L3 `Fragment` — done** (`fragment.rs`; S1→S2→S3, 10 tests).
- **L4 `mol!` / `frag!`** — design + staged plan above (S1→S2→S3); the last layer, ready to build.
- Overlays-in-*fragments* (L3) remains gated on the port/operad spike (103).
