# InChI 1.07: Rust migration vs vendor-and-wrap

Status: **Informational**
Date: 2026-08-09
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[190](190-enumeration-algorithm-candidates-2026-08-08.md)

## Scope

Findings from a review of the InChI checkout at `materials/codes/InChI`,
assessing two questions:

1. What it would cost to migrate the InChI 1.07 implementation to Rust while
   keeping the semantics of revision 1.07 faithfully, versus vendoring the C
   and wrapping it behind a Rust interface the way `umol-nauty-sys` and
   `umol-msym-sys` wrap nauty and libmsym.
2. What the actual improvement of either path is for the two stated purposes:
   natural integration (error handling in particular) at reaction-network call
   volumes (10⁸+ calls), and access to the genuinely innovative algorithms —
   canonicalization and tautomer handling — that are buried in the code.

Non-goals, fixed up front: no public fork, no evolution of InChI or InChIKey
construction, no overlays. RInChI is covered separately at the end.

All file references are relative to
`materials/codes/InChI/INCHI-1-SRC/` unless noted.

## The checkout is not a release

The working tree is the `dev` branch at `v1.07.5-261-g89f8dd4` — 261 commits
past the v1.07.5 tag (tagged 2026-02-17), HEAD dated 2026-08-03. The version
string `CURRENT_VER "1.07.5"` lives, of all places, in the bounds-checking
header `INCHI_BASE/src/bcf_s.h:52`. Two consequences:

- **Two unreleased, semantics-changing features are present**: Molecular
  Inorganics (revised bond-disconnection decision tree; a Ti bis-flavonoid
  case changes from disconnected `2C15H10O7.Ti` to connected `C30H18O14Ti`)
  and Enhanced Stereochemistry (emits an `InChI=1B` beta prefix). Both are
  runtime-gated and off by default, but one dev commit ("skip implicit H
  addition for metal atoms") changed default-path behavior and the bundled
  regression references were regenerated after it (2026-06-23).
- **A faithful target must be the v1.07.5 tag**, not this HEAD, and reference
  outputs must be regenerated at the tag.

One more tree-local surprise: `MAXVAL` in `inchi_api.h:105` was changed from
20 to 50 ("for testing purposes", per the inline comment). `MAXVAL` shapes
`inchi_Atom` (240 bytes here vs 120 stock), so the public ABI of this tree is
incompatible with any stock `libinchi.so`. A `-sys` crate must vendor and
compile the C; linking a system library is not an option.

What "1.07 semantics" means: the 1.07 series deliberately did not change the
algorithm. 1.07 is a hardening release (33 buffer overflows, 157 null derefs,
2,480 general fixes, all inline, 2,859 `djb-rwth` markers). The chemistry to
reproduce is 1.06's chemistry plus a small set of 1.07 fixes that sit on
semantic paths: nine rewritten math functions in the stereo-perception file
`ichister.c`, ~25 repaired conditionals across normalization and I/O files,
and FixedH/RecMet round-trip fixes. AuxInfo output is identical to 1.06 again
since 1.07.1 (1.07.0 briefly diverged and was reverted).

## Shape of the codebase

| Unit | LOC (.c + .h) |
| --- | --- |
| Whole repo C | ~225,400 |
| `INCHI_BASE/src` (the engine) | 175,237 |
| Compiled into `libinchi.so` | 177,146 (.c) |
| Compiled into the `inchi-1` CLI | 163,769 (.c) |
| `INCHI_API/libinchi` (API layer, incl. IXA) | 19,182 |
| Demos (not shipped) | 23,852 |

Subsystems of `INCHI_BASE`, with the parts a molecule→InChI path needs marked:

| Cluster | LOC | On mol→InChI path |
| --- | --- | --- |
| Normalization + tautomers (`ichi_bns.c` 12.5k, `ichitaut.c` 7.1k, `ichiqueu.c`, `ichiring.c`, headers) | ~23,800 | yes |
| Canonicalization (`ichican2.c` 6.9k, `ichicano.c`, `ichicans.c`, `ichimap1/2/4.c`, `ichisort.c`, headers) | 21,317 | yes |
| Stereo perception (`ichister.c`) | 4,995 | yes |
| Preprocessing (salt/metal disconnection, fix-ups; `strutil.c`, `runichi*.c` slices) | ~10,000 | yes |
| INChI assembly + serialization (`ichimake.c` 6.3k, `ichiprt1–3.c` 12.8k, key) | ~20,000 | yes |
| InChI string reader (`ichiread.c`) | 12,587 | no |
| InChI→structure reconstruction (`ichirvr1–7.c`) | 31,804 | no |
| Underivatization + ring-chain tautomerism (in `ichinorm.c`) | ~6,900 | no (optional features) |
| Option parsing (`ichiparm.c`, one 1,526-line function) | 3,196 | partially |

Legibility, measured rather than vibes: 1,319 top-level functions, median 43
LOC but p99 = 982; 45 functions exceed 500 lines. The largest on the forward
path are `MarkTautomerGroups` (1,997), `GetBaseCanonRanking` (1,806),
`Canon_INChI3` (1,354), `OutputINChI1` (1,030). The single largest function in
the tree, `FixFixedHRestoredStructure` at 6,202 lines, is on the reverse
(InChI→structure) path and out of scope. 2,665 `goto`s, of which ~72% are the
mechanical single-exit-cleanup idiom (`goto exit_function`) that translates
directly to `?`/RAII; ~600 are genuine jumps needing restructuring. 2,399
preprocessor conditionals; `mode.h` alone has 405 `#define`s, including
compile-time switches that change chemistry (`DISCONNECT_SALTS`,
`MOVE_CHARGES`, tautomer-rule gates). Hungarian naming is consistent
(`b`/`n`/`sz`/`p`). Twelve of the fifteen largest files carry no purpose
comment. Shipped optimization level is `-O1` in both build systems.

Global state: eliminated in 1.07.2 (moved into `CANON_GLOBALS`,
`INCHI_CLOCK`, `INPUT_PARMS` context structs threaded everywhere). The
structure→InChI direction is de-facto concurrently callable — the `mol2inchi`
demo runs a pthread pool over `MakeINCHIFromMolfileText`, and CI has a
multithreading test — with two benign-in-practice data races (a canonicalizer
bit-mask pair rewritten to the same constant on every call; a lazily
initialized halogen table) and one real corruption hazard on the
unmapped-error path (`ErrMsg` returns a pointer to a shared `static
char[64]`). The InChI-parsing direction is categorically not thread-safe: a
256 KB file-scope line buffer `szLine_i2i` in `readinch.c` is shared scratch.
No document in the tree states any thread-safety guarantee.

## The two innovative cores

These are the pieces worth writing down regardless of which path is taken.
Both are real algorithms with literature citations, buried under chemistry
plumbing.

### Tautomer and charge normalization is a balanced-flow problem

The core of `ichi_bns.c` implements Kocay–Stone balanced network flows
(cited in-code at `ichi_bns.c:10680`: J. Comb. Math. Comb. Comput. 19 (1995)
3–31). The construction:

- Each atom becomes a doubled vertex pair; each bond an edge with
  flow = bond order − 1. A vertex's st-edge capacity minus flow equals its
  count of unsatisfied valence "dots"; normalization is "eliminate dots by
  augmenting along alternating paths".
- Mobile hydrogens and mobile negative charges are literally flow units on
  edges from endpoint atoms to fictitious *t-group* vertices; movable positive
  charges to *c-group* vertices. Chemistry-specific veto hooks inside the
  search prevent H↔(−) exchange on non-acidic atoms.
- The search (`BalancedNetworkSearch`, ~2,700-LOC kernel including blossom
  contraction with union-find path compression) answers one query shape:
  "does an alternating path exist that moves this H/charge from A to B" —
  `bExistsAltPath` is the primitive every tautomer rule calls.
- Kekulization falls out of the same machinery: aromatic input bonds are
  forced single on entry and repaired by flow augmentation
  (`BnsAdjustFlowBondsRad`); a bond is marked "alternating" iff forcing its
  alternative order still admits a feasible flow (`BnsTestAndMarkAltBonds`).
  Alt-bond detection is a flow feasibility test, not ring perception.
- Proton normalization (the mobile-H layer content) is a fixed ladder —
  simple/hard H⁺ removal from N/P/O⁺, then acidic-proton removal or addition,
  each "hard" step phrased as flow between two t-groups — matching Technical
  Manual §5.1–5.2 step for step.
- The fixed-H layer is obtained by re-running the same entry point
  (`mark_alt_bonds_and_taut_groups`) with tautomer perception disabled
  (`t_group_info = NULL`); mobile vs fixed H is decided by the caller
  (`Create_INChI`), not inside the normalizer.

The tautomer *rules* sitting on top are the opposite of elegant: thirteen
recognized patterns (1,3-shift base rule, keto–enol, 1,5, pyridinol,
pyrazole, tropolone variants, and seven experimental `PT_*` rules), each an
~165-line copy-paste clone of the base block differing only in its endpoint
predicate and path mode. SMIRKS strings exist only as comments; nothing is
data-driven beyond a handful of small element tables. The cluster shares no
mutable state with canonicalization except `T_GROUP_INFO`, which
canonicalization reorders in place.

A caution for any extraction: the 31.8k-LOC InChI-reversing subsystem
(`ichirvr1–7.c`) is built on the same `BN_STRUCT`/`BN_DATA` types, so in the
C tree the flow engine serves both directions. A forward-only port severs
that dependency deliberately.

### Canonicalization is McKay 1981, transcribed literally

`CanonGraph` (`ichican2.c:3662`, 1,001 lines) cites and follows "Practical
Graph Isomorphism", Congressus Numerantium 30 (1981) 45–87 — the paper's step
numbers survive as C labels `L2:`…`L17:` and the variables keep the paper's
names (zeta, rho, theta, Omega, Phi). It is nauty's algorithm but not nauty's
code: naive sort-based 1-WL refinement to fixpoint (no worklist), first
non-singleton target cell (no smallest-cell heuristic), the three classic
prunings (automorphism mcr/fix sets, first-leaf, best-leaf CT comparison),
group order accumulated but generators discarded.

Two InChI-specific extensions carry the layer structure, and doc 186 already
cites the first:

- **Constrained-prefix mode** (`zb_rho_fix`): a `CanonGraph` call can be
  handed the connection-table certificate of a previous, coarser pass; the
  search is then constrained to reproduce that certificate on earlier layers
  and minimize only the new layer. `GetBaseCanonRanking` (1,806 lines) chains
  up to eight passes this way: skeleton → +non-tautomeric H → +isotopes →
  +t-groups (as extra vertices, directed edges) → +tautomer isotopes →
  fixed-H → fixed-H isotopes, each pass skipped when the previous orbits
  already separate the new colors. This is progressive refinement over a
  compact base graph — no bond-subdivision into edge vertices, confirming the
  benchmark question doc 186 poses.
- **Layered certificate comparison**: the leaf certificate carries parallel
  color arrays (H counts, fixed H, isotope keys) compared in a hard-coded
  six-layer order with memoized per-layer partial results.

The inner core is remarkably separable: `CanonGraph` plus the
partition/cell/nodeset/certificate machinery and refinement primitives is
**~5,600 LOC with zero references to atoms, elements, or parities** — its
entire chemistry interface is one adapter (`CreateNeighList`) and four color
arrays. What is not generic: the atom/t-group vertex split is baked into the
certificate layout, and the six comparison layers are hard-coded. Both are
expressible as parameters of a lifted module.

Stereo is a *separate second phase*, not part of the numbering loop: parities
are perceived from geometry before canonicalization (`ichister.c`), converted
to canonical parities via rank-permutation counting, and then a second,
independent backtracking search (`map_stereo_bonds4`/`map_stereo_atoms4`,
~7,000 LOC with the equivalence machinery) enumerates rank-preserving
renumberings minimizing the stereo CT, twice (once for the mirror image, to
decide the `m` layer). This search implicitly re-walks the automorphism group
that `CanonGraph` computed and threw away — a structural inefficiency worth
knowing about, not fixing (non-goal).

## API surface: what wrapping has to live with

The API is better than its reputation in one respect and worse in several.

**Input is not a string.** `inchi_Input` is a structured atom array
(`inchi_Atom`: element name `char[6]`, adjacency + bond types/2D-stereo up to
`MAXVAL`, implicit-H counts with a `-1` = "add automatically" convention,
isotope, radical, charge), plus optional `inchi_Stereo0D` records. The struct
path is complete — no molfile text round-trip anywhere
(`ExtractOneStructure` converts directly to internal atoms). Strings appear
in exactly four places: the options string, and the outputs
`szInChI`/`szAuxInfo`/`szMessage`/`szLog`.

**Options.** One string, tokenized per call (MSVC-style splitter, in-place,
max 32 tokens, silent truncation beyond that) and fed to the same 1,526-line
`ReadCommandLineParms` the CLI uses — 133 recognized literals,
case-insensitive. The prefix is `/` on Windows and `-` elsewhere;
`INCHI_ALT_OPT_PREFIX` is defined but dead. Failure modes: a *prefixed*
unknown option is reported only as text in `szLog` ("Unrecognized optionQ3")
with `inchi_Ret_OKAY` returned; an *unprefixed* token is silently treated as
an input file name and ignored entirely. `GetStdINCHI` does not reject
non-standard options — it silently scrubs them.

**Errors.** Return codes are six coarse values (OKAY/WARNING/ERROR/FATAL/
EOF/UNKNOWN; `inchi_Ret_BUSY` is a vestige never produced). All detail is
free text: fragments accumulated into a 256-byte buffer, joined with `"; "`,
truncated with `"..."`, e.g. `"Unknown element(s):"` + name, `"Unsupported in
this mode element '*'"`, `"Too many atoms [did you forget 'LargeMolecules'
switch?]"`. Severity is *derived from* numeric ranges of an internal `err`
value, not stored. `MakeINCHIFromMolfileText` clamps its negative codes to
OKAY on one path, so callers must additionally check `szInChI != NULL`.

**Sharp edges relevant to a wrapper.** `GetINCHI` dereferences the input
before any null check; the pseudoatom pre-check loops over `i` but reads
atom 0 only (`inchi_dll.c:301`), so a `*` at index ≥ 1 bypasses it;
`GetINCHIEx` *mutates the caller's array* (rewrites `"*"` → `"Zz"`), so the
Rust signature takes `&mut` or copies; element names are copied with an
unbounded `strcpy` into `char[6]`; `szAuxInfo` is an interior pointer into
the `szInChI` allocation and must never be freed separately; plain `GetINCHI`
cannot express polymers at all (only `GetINCHIEx`).

**The `*` freak-out, resolved.** `*`/`Zz` is the polymer pseudoatom. Without
`/Polymers*` or `/NPZz`, the struct API rejects it ("Unsupported in this mode
element '*'"; deeper in, "Invalid element(s):" / "Non-polymer-related Zz/star
atoms are not allowed"). There is no query-atom concept; any unknown symbol
is just "Unknown element(s)". For umol-constructed inputs this entire class
of failure is avoidable before the FFI boundary.

**IXA.** The "extensible" API (opaque handles, per-call status object holding
up to 50 typed-severity free-text messages) is the only part with a designed
error model, but it is a façade: internally it *builds an option string*,
re-parses it through the same parser, copies the molecule one extra time, and
caps atom degree at 20 regardless of `MAXVAL`. It adds work; it removes none.

### Per-call fixed cost — the 10⁸-call concern

The string-based exchange itself is not where the cost is. The molecule goes
in as structs; the InChI comes out as a string, which is irreducible — the
string is the product. The measured-by-reading fixed overhead per `GetINCHI`
call is:

1. `inchi_strbuf_init` → **a 256 KB zero-filled `calloc` on every call**
   (`inchi_dll.c:592`, size constants in `ichi_io.h:100`), freed before
   return. A 1 KB "smaller initial size" constant exists and is unused here.
2. Two more ~32 KB `calloc`s for the log and output streams (the log always
   receives a parameter banner — `PrintInputParms` writes "Generating
   standard InChI" even with no options). Growth policy is
   calloc-copy-free in 32 KB steps, quadratic in output length.
3. `strdup` + tokenization + the full 1,526-line option parser, per call.
4. ~1.7 KB of stack context memsets (`STRUCT_DATA`, `INPUT_PARMS`,
   `CANON_GLOBALS`, …), then a low-hundreds count of small mallocs through
   normalization/canonicalization (`AllocateCS` alone is ~25–30 `calloc`s,
   run twice per component). 347 allocation call sites are reachable on the
   generation path.

So ≈ 328 KB is allocated and zeroed per call regardless of molecule size —
tens of microseconds of memory bandwidth, comparable to or exceeding the
actual computation for a small molecule — plus parse and formatting work. At
10⁸ calls the zeroing alone is tens of terabytes of memory writes. **None of
this is avoidable through the public API**: there is no reusable context
object; the staged `INCHIGEN_*` API re-runs the option parser and banner too.
The presence of a bundled `tbbmalloc_proxy` option (Windows-only, commented
out) suggests allocator pressure is a known issue upstream.

Two caveats keep this honest. First, these are structural findings, not
measurements; a micro-benchmark (stock `GetINCHI` vs a bypassing shim on
representative molecules) is the obvious next step and is listed under open
items. Second, the InChI computation proper — normalization, two
canonicalization layer sets, stereo mapping, serialization — is the floor
under both paths; neither wrapping nor rewriting makes that free.

## Fit to umol's model

Independent of path, InChI as an identifier has semantics that differ from
umol's structural identity:

- Standard InChI **normalizes away** protonation states and mobile-H
  tautomers and **disconnects metals** (RecMet restores them only as a
  non-standard layer). For network node identity this is a decision, not a
  detail: InChIKey identity merges species that `MoleculeAst` distinguishes
  (doc 186's canonical form is exact by construction). Depending on the
  network semantics, tautomer-invariant node identity may be exactly right
  (species-level networks) or wrong (elementary-step networks).
- The input model covers the organic subset only: bond orders 1/2/3 (type 4
  "ALTERN" is documented "avoid by all means" — kekulize first, which umol
  can), no dative, multicenter, or noncovalent entities, atom count capped at
  1024 by default (32766 with `/LargeMolecules`, which switches the output to
  a beta prefix). Organometallic and multicenter chemistry — core scope for
  umol — is representable only through InChI's disconnection conventions.
- Doc 186 already positions InChI correctly for the internal question: an
  algorithmic precedent and comparison oracle, not the internal canonical
  form. An InChI/InChIKey emitter is interop and reporting surface.

## Path 1: vendor and wrap (the nauty/msym pattern)

What the existing pattern buys, transferred to libinchi:

- **Build**: `cc` in build.rs over the `INCHI_BASE/src` + `libinchi/src`
  closure (~107 files). Symbol-prefixing nauty-style (71 `-D` renames) does
  not scale to 1,319 functions; `objcopy --localize`/linker version scripts
  or simply accepting the static-link symbol space are the realistic options.
  Vendoring is mandatory anyway because of the ABI-shaping constants
  (`MAXVAL`, `MAX_ATOMS`).
- **Boundary**: umol constructs `inchi_Atom[]` from its own already-validated
  model, so the input-validation half of InChI's error surface (unknown
  elements, bond-to-nonexistent-atom, pseudoatoms, molfile syntax) is
  unreachable by construction. The remaining C-side errors map to a typed
  Rust enum from the coarse return code, with the free-text `szMessage`
  carried as payload — the same shape `NautyError` has. Known message
  strings ("Accepted unusual valence(s)", "Charges were rearranged") can be
  recognized, but that is string-matching on unversioned text; it works
  because the vendored copy is frozen.
- **Reentrancy**: forward-direction concurrent calling works today; the two
  benign races and the `ErrMsg` static are one-line fixes in a vendored copy.
  The InChI-parsing direction stays single-threaded or gets a lock; umol
  barely needs it.

The per-call fixed cost is the fork in this path:

- **Shallow wrap** (call `GetINCHI` as published): every pathology of the
  previous section is kept — 328 KB churn, option re-parse, banner
  formatting, per call. Integration improves; cost does not.
- **Deep shim** (the shim calls internal entry points): because globals are
  gone, a umol-owned C shim can pre-parse options once into a frozen
  `INPUT_PARMS` snapshot, hold per-thread reusable stream buffers, skip
  `PrintInputParms` and the log stream entirely, and enter at the
  `ExtractOneStructure`/`ProcessOneStructureEx` level. That removes
  essentially all fixed overhead while leaving the algorithm untouched. It is
  more shim than nauty's 210 lines — it re-plumbs the outer layer of
  `GetINCHI1` — but it patches *around* the engine, not inside it. Since
  upstream tracking is an explicit non-goal and the copy is frozen at
  v1.07.5, the usual objection to depending on internal symbols (drift) does
  not apply.

What wrapping does not deliver, at any depth: legibility. The algorithms stay
buried; the improvement on purpose (2) is zero. It also keeps ~177k LOC of
`-O1` C with a documented history of memory-safety findings in the process —
mitigated by the fact that umol feeds it only self-constructed structured
input, never hostile text.

## Path 2: faithful Rust reimplementation

Scope accounting for a molecule→InChI(+Key) port — the direction umol needs:

| Included | C LOC |
| --- | --- |
| Preprocessing (fix-ups, salt/metal disconnection) | ~10,000 |
| Normalization + tautomers (BNS) | ~23,800 |
| Canonicalization + stereo | ~26,300 |
| Assembly + serialization + InChIKey (SHA-224) | ~20,000 |
| Element data, utilities, option subset | ~8,000 |
| **Total** | **~88,000** |

Excluded by scoping to the forward direction: the InChI reader (12.6k), the
InChI→structure reconstruction (31.8k), underivatization/ring-chain (6.9k),
the API layer (19.2k), polymer support if Zz-chemistry stays out of scope,
and the CLI. That is roughly half the engine.

The Rust volume would be substantially below 88k: the seven `PT_*` tautomer
rules are clones of one block (a rule table in Rust), 72% of the gotos are
mechanical, the preprocessor layers (`TARGET_*` variants, dead `#if` code,
WinChI paths) vanish, and the four-way copies of stereo/isotope arrays in
`sp_ATOM` become generics. Against that: the port must reproduce, exactly,
several thousand lines of order-sensitive fixpoint loops (the normalization
main cycle), the eight-pass constrained canonicalization chain, and the
stereo mapping search — code where *the iteration order is the
specification*. Any "cleanup" that changes which of several equal-cost
outcomes is reached first changes output strings.

Fidelity is the actual cost driver, not volume. The verification
infrastructure in `INCHI-1-TEST` transfers well:

- **Byte-for-byte oracle**: SQLite references with full
  InChI + Key + AuxInfo + exit code for 4,190 structures (2,190 InChI legacy
  set + 2,000 mcule), driven over `ctypes` against any `cdylib` exporting
  `MakeINCHIFromMolfileText` and `GetINCHIKeyFromINCHI` — a Rust
  implementation drops in by exporting two symbols. References must be
  regenerated at the v1.07.5 tag (the bundled ones encode dev-branch
  behavior).
- **Invariance harness**: 10 atom-order permutations per structure, output
  must be identical — implementation-independent, no reference needed.
- **Scale**: 4,190 structures is thin. The 1.06 normalization fix affected
  ~1,700 of ~102M PubChem entries — 1 in 60,000. Divergences of that
  frequency are invisible until differential testing at PubChem scale, and
  the harness already has scale-out configs (`config_pubchem_*.py`) that
  download shards from NCBI. A faithful-port claim requires that run, with
  the vendored C as the oracle — which means Path 1's build exists inside
  Path 2 anyway, at least as a test fixture.
- The GoogleTest tier links internal C symbols by name and does not transfer.

What the rewrite delivers that no wrap can: typed errors end to end
(severity as types, not derived integer ranges; no free-text protocol); zero
fixed overhead by construction (contexts are values; buffers are reused;
options are an enum set resolved at compile time); the two core algorithms as
legible, documented Rust adjacent to umol's own canonicalization work; and
memory safety in the engine rather than at its boundary. What it costs: on
the order of a person-quarter to person-half-year of focused work — ~88k LOC
of C where perhaps a fifth is genuinely subtle, plus the PubChem-scale
differential campaign, which is compute and iteration time, not just code.
That estimate is structural, not a schedule.

## Decoupling the two purposes

The two stated purposes have different faithfulness requirements, and that
asymmetry is the main structural finding:

- **Purpose 1 (integration + call cost)** needs InChI-compatible *output*.
  Both paths deliver it; the deep shim delivers the cost profile without any
  reimplementation risk.
- **Purpose 2 (write down and reuse the innovative pieces)** does not need
  InChI compatibility at all. The balanced-flow normalization idea (mobile
  H/charge as flow, alt-bonds as flow feasibility) and the constrained-prefix
  layered canonicalization are *algorithm schemata* that could inform
  umol-native machinery (doc 186's progressive refinement question, resonance
  and protomer handling on `MoleculeAst`) under umol's own model — where doc
  186 explicitly warns against importing InChI's normalization semantics.
  Harvesting them is a documentation-plus-reimplementation effort against the
  ~2,700-LOC flow kernel and the ~5,600-LOC canonicalizer core, with no
  obligation to reproduce InChI's output, rule tables, or quirks.

A full faithful rewrite bundles these two deliverables; the
wrap-plus-harvest combination unbundles them. The bundled version's extra
value over the unbundled one is confined to: one toolchain instead of a C
dependency, typed errors inside the engine rather than at its boundary, and
the absence of unsafe C in-process. Its extra cost is the fidelity campaign.

## RInChI

No RInChI checkout exists under `materials/`; findings are from the v1.00
paper (Grethe, Blanke, Kraut, Goodman, J. Cheminform. 10:22, 2018;
`materials/representations/Grethe2018-RInChI.pdf`). Verification against the
reference implementation (InChI Trust download; C++ around libinchi, RXN/RD
input) needs a checkout — open item.

What RInChI is: a *string-level combinator* over standard InChI. Six layers —
version (`RInChI=1.00.1S`), two molecule-group layers (component InChIs with
prefixes stripped, sorted alphabetically within and between layers, separated
by `!` and `<>`), agents, a direction flag (`/d+`, `/d-`, `/d=`), and
no-structure counts (`/u2-0-1`). RAuxInfo mirrors the first four layers.
Three keys: Long (concatenated component InChIKeys), Short (fixed 63 chars,
hashed major/minor layers per group plus protonation-state letters), Web (47
chars, all components pooled, sorted, deduplicated, hashed — role-agnostic).
It performs no graph work of its own; all structure identity is delegated to
standard InChI **1.04 pinned** — the version letter is baked into the key
prefixes, so regenerating RInChIs over a 1.07 engine is technically trivial
and nominally out of spec at the same time.

Declared and observed limits: reactant/product roles live only in the
direction flag (alphabetic order decides which side is layer 2); v1.00 states
a no-duplicates rule but implements no duplicate check; agents are "present
at both ends", by design mechanism-free; stereo inherits standard InChI's
absolute-only handling, so racemic/relative reaction stereochemistry is
represented arbitrarily; and — decisive for umol — **atom correspondence is
discarded entirely**. A RInChI identifies a reaction as a (sorted multiset of
component InChIs, direction) pair. umol's planned reaction canonicalization
(doc 186 S11–S12, over `ReactionSpanAst` with correspondence intact) is
strictly finer; RInChI cannot serve as the internal reaction identity, only
as an interop/reporting identifier — a role it fits well precisely because
the literature uses component-level identity.

Cost assessment: given any molecule-level InChI provider (either path above),
RInChI generation is string sorting plus the InChIKey hash machinery
(SHA-224 truncation, letter encodings) — small, self-contained Rust with no
C beyond what the InChI provider already uses. Wrapping the official C++
implementation instead would add a C++ toolchain dependency to gain code
whose entire algorithmic content is the previous sentence. The "learn from
it" yield is the layer/key design itself (particularly the Web-key idea:
role-agnostic pooled identity for cross-database search), not the
implementation. Exact key compatibility with the reference implementation
would need byte-level verification against its outputs — checkout required.

## Cost / improvement summary

| Axis | Vendor + shallow wrap | Vendor + deep shim | Faithful Rust rewrite |
| --- | --- | --- | --- |
| Effort | ~nauty-sys scale + build plumbing | + re-plumbed outer layer in a umol-owned shim | ~88k C LOC forward-path port + PubChem-scale differential campaign |
| Error handling | typed at boundary; C free text as payload | same | typed end to end |
| Per-call fixed cost (10⁸ calls) | unchanged: ~328 KB zeroed + option parse per call | removed; floor = actual computation | removed by construction |
| Output fidelity to 1.07.5 | exact (it *is* 1.07.5) | exact (engine untouched) | claimed, must be proven; 1-in-60k tail risk |
| Algorithm legibility (purpose 2) | none | none | full, but bundled with fidelity work |
| Thread safety | forward path OK after 3 one-line fixes | same | inherent |
| Unsafe C in process | ~177k LOC at `-O1` | same | none |
| Upstream InChI revisions | re-vendor + regenerate references | re-vendor + re-check shim entry points | re-port deltas (non-goal, but a real asymmetry) |

The harvest of the two core algorithms (~2.7k flow kernel, ~5.6k
canonicalizer) into documentation and/or umol-native machinery is available
under every column and is independent of the identifier-generation choice.

## Open items

- Micro-benchmark: stock `GetINCHI` vs a deep-shim prototype on
  representative molecules, to convert the structural per-call findings into
  numbers before weighting the wrap columns.
- Re-pin the vendored tree to the v1.07.5 tag and regenerate the SQLite
  references there (the bundled ones encode post-tag dev behavior).
- Decide the option subset that constitutes "1.07 semantics" for umol:
  standard InChI only, or + `RecMet`/`FixedH` (the organometallic-relevant
  non-standard layers). This bounds both the shim surface and any rewrite
  scope.
- Obtain the official RInChI implementation checkout to verify the paper
  reading and key construction against code.
- Decide where the algorithm write-up lives (whitepaper rationale doc vs a
  dedicated discussion doc) once either is scheduled.
