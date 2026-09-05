# 205 — Atom-mapping test corpus

Status: Proposed
Date: 2026-08-20
Relates: [200](200-molecular-data-substrate-2026-08-19.md),
[201](201-molecular-data-first-steps-2026-08-19.md),
[203](203-atom-mapping-2026-08-19.md),
[207](207-reaction-network-spike-2026-08-24.md)

## Purpose

This document scopes the evidence and algorithm work needed for the next atom-
mapping iteration. It uses the external-representation and visualization path
established experimentally in doc 201, but it is not a staged implementation
plan and does not yet settle a public objective or configuration surface.

The work has two coupled goals:

1. improve annotation by combining automatic evidence--proved minimum-edit
   results where available, operation-issued classification mappings, external
   mapper output, and structural analyses--with high-value manual inspection;
   and
2. develop the mapping algorithm against the resulting annotated corpus.

The central artifact is an atom-mapping **test corpus**: structured reaction
inputs, neutral atom mappings and source evidence, exact results where
available, score decompositions, named equivalence relations, provenance,
bounded human judgments, and algorithm-run measurements. It should support
executable comparisons and regression tests, not merely a gallery of
interesting examples.

Manual review is the scarce resource. The expected primary campaign is at most
roughly one thousand cases; collaboration may extend that to a few thousand,
but no design may assume manual annotation of the approximately two hundred
thousand classification cases or repeated annotation campaigns after a
sequence of small algorithm changes. Automatic evidence must perform the
population-scale triage, and the anticipated algorithm and objective
improvements must be assembled before the main manual campaign.

The corpus must not be implemented as a private atom-mapping file format or a
one-off report generator. Atom mapping is the first demanding case study for
portable `Molecule`, `Reaction`, `ReactionSpan`, and correspondence records;
columnar storage; direct analytical queries; SVG export; and linked Jupyter
inspection. The generic substrate should emerge from real use here while
remaining useful to reaction-template libraries and reaction networks.

## Starting point from doc 203

Doc 203 delivered an exact, non-parametric structural baseline. Its objective
maximizes element-compatible atom matches and then minimizes unit localized-
bond additions, deletions, and order modifications. It enumerates every labeled
optimum and makes no mechanism claim.

The first RXNMapper comparison established two independent problems:

- the current branch-and-bound implementation proves small reactions but loses
  coverage rapidly as reaction size grows;
- the exact optimum of the current objective can be chemically wrong, as in
  the wildcard/hydrogen exchange in `Rh_12116`, and endpoint structures alone
  cannot identify the intended parent Diels--Alder map.

These must remain distinct in the next investigation. A faster solver for the
same objective changes execution but not semantics. New compatibility rules,
objective terms, anchors, or mechanism-informed selection define different
mapping semantics even if they reuse the same search machinery. Alternative
exact formulations can serve either purpose only after their contract is made
explicit.

The practical consequence is that the exact baseline is a reference only on
the subset for which it establishes an optimum. A time or resource limit on a
larger reaction produces an unknown result, not a weaker certificate and not
evidence that an external mapping is optimal. Refining the mapper only against
that tractable subset would bias the investigation toward small systems.
Manual annotation is useful for selected disagreements, but it cannot replace
reference evidence over a corpus containing hundreds of thousands of
reactions.

## One shared data path

The test corpus should exercise the same logical records in three physical
forms:

1. ordinary owned umol values used by the mapping and reaction operations;
2. Arrow-compatible columnar batches persisted through Parquet or a comparable
   analytical format;
3. notebook projections combining queryable tables with SVG molecular and
   reaction views.

The first implementation may use the version-scoped structural and canonical
identity proxies from doc 201. It must preserve source identifiers and
algorithm versions so all derived keys can be rebuilt when canonicalization or
the record schema changes.

The design should factor reusable chemistry records from experiment records.
`Molecule`, `Reaction`, `ReactionSpan`, and `Correspondence` values should not
be embedded in an atom-mapper-specific opaque payload when the same boundary
representation can serve template libraries or network derivations. Conversely,
run configuration, timing, solver statistics, confidence, and adjudication are
experiment facts rather than fields of the molecular model.

## Repository and distribution boundary

The corpus may remain a separate, non-published crate in the same Rust
workspace because that is a convenient place to reuse umol and develop the
experiment. This does not make corpus generation, source-specific schemas,
annotation workflows, query scripts, or notebooks part of the umol library.
They build on umol; they are not umol proper and must not shape its public
molecular API around the needs of this one study.

`umol-store` is different: faithful records and generic storage,
reconstruction, scan, and query capabilities are within umol's scope. Corpus
relations and corpus-specific operations remain on the consumer side of that
boundary.

The `umol-py::store` feature added by doc 201 is a temporary experimental
bridge. It may be used during this investigation, but it is not a durable home
for corpus or storage tooling and must not grow additional corpus-specific
surface. Before doc 201 is closed, that bridge must be removed or moved behind
a separately scoped distribution boundary. The core `umol` Python package must
not accumulate the corpus, database, query-engine, or annotation stack.

## Logical corpus contents

The following are required capabilities, not a settled table schema.

### Mapping inputs and molecular frames

Each case needs the reported source form and every normalized mapping input used
for comparison. A `MappingInput` identifies one ordered lhs/rhs pair of
`MoleculeRecord` keys together with its source, parsing and standardization
path, and umol version. Each molecule fixes one local entity-id frame. A
stored atom mapping is interpreted only against the pair identified by its
mapping input.

Raw and standardized readings of one source case are different mapping inputs.
When their relation is known, an alignment records separate lhs-to-lhs and
rhs-to-rhs correspondences. Comparisons of raw atom labels or edit scores are
valid only within one mapping input or after applying that explicit alignment.

Concrete reactions should use `Reaction` as the central umol representation.
`ReactionSpan` and condensed reaction graphs are derived analytical and visual
projections. The corpus must be able to retain both a source-provided atom map
and maps computed within umol without treating either as intrinsic to the
unmapped sides.

### Atom mappings, derived analyses, and runs

`AtomMapping` is the neutral chemical relation settled in doc 201: a complete
partial atom correspondence interpreted against one mapping input. It makes no
claim about optimality, chemical correctness, producer, or role. Source-
specific relations connect that value to an RXNMapper output, an operation-
issued network reconstruction, a manual assertion, or another concrete source
of evidence. Producer confidence remains on the producer output rather than
becoming an intrinsic property of the mapping.

Feasibility, edit decomposition, objective score, and equivalence-class
membership are derived from a mapping, its input, and an explicitly named
analysis or objective. They should remain computable on demand unless a
concrete query demonstrates that materializing one is useful. Exact molecular
symmetry, labeled correspondence equality, and induced-reaction canonical
equality remain distinct relations.

A mapping run names its mapping input, algorithm, version, configuration,
execution outcome, timing, and emitted atom mappings. Algorithm-specific
results retain only the facts that algorithm can establish, including
objective values, search statistics, proof status, and enumeration status.
There is no universal candidate-evaluation row or universal proof field
attached to every mapping or run.

For the first objective, feasibility means that an atom mapping has the mapping
input's carrier counts, is an in-range partial bijection, and pairs only
element-compatible atoms. It does not assert maximum cardinality, global
optimality, chemical correctness, or selection of a symmetry representative.
Proof status states whether the run established that no feasible correspondence
has a better score under the named objective. It is an operational global-
optimality claim, not a formal proof artifact. Enumeration completeness
separately states whether every labeled correspondence attaining that optimum
was emitted. An interrupted enumeration may therefore follow a proved optimum,
while a best observed score without a completed optimization remains unproved.

Complete labeled optima remain the semantic result of doc 203. Equivalence-
class representatives and counts are analytical projections, not replacements
for the underlying correspondences.

### Provenance and adjudication

External confidence, source versions, preprocessing, literature or database
identifiers, and precise record locations must survive ingestion. Curated
judgments should distinguish at least:

- known or strongly supported atom identity;
- symmetry-equivalent alternatives;
- mechanistically ambiguous cases;
- mappings contradicted by the available evidence; and
- unresolved cases.

An annotation must retain a reference to the evidence bundle presented for
review and must not silently become a hard constraint or a training label. The
reviewer does not manually recode those automatic evidence references. This
permits the same corpus to test exact search, compare approximate mappers, and
study alternative objectives without pretending that every source has equal
authority.

Manual annotations should state the supported identity claim, alternatives,
ambiguity, and evidence as directly as the source permits. Whether A0, a later
objective, or a particular implementation agrees with that claim is a derived
comparison and can be recomputed. The objective used to select a review queue
may change which cases receive attention, but it must not make the resulting
human judgment unusable after the algorithm changes.

The annotation interface may be validated on a small pilot before the main
campaign. That pilot is for checking targets, evidence links, and review
mechanics, not for running repeated rounds of manual review against successive
minor algorithm variants. The main campaign uses stratified and disagreement-
driven selection within the stated human-review budget.

#### Review operation and durable result

Each review item presents one existing `AtomMapping` as the proposed
correspondence for one mapping input. Other exact, external, or network-derived
mappings may appear beside it as evidence, but the manual operation remains
deliberately narrow:

- **confirm** the proposed correspondence; or
- **correct** the proposed correspondence by changing only the atom assignments
  that differ.

Correction is an interface operation, not a second patch-based corpus
representation. The interface starts from the complete proposed mapping and
lets the reviewer remove, replace, or add only changed assignments. It applies
all deviations together, checks the resulting carrier counts, id ranges, and
partial-bijection shape, and shows the resulting full correspondence before
saving it. The reviewer never has to re-enter every unchanged atom pair.

Every saved review persists the complete resulting atom correspondence and all
of its atom pairs, including when confirmation merely reproduces the proposal.
The review links that durable result to the mapping input, proposed mapping,
presented evidence bundle, and confirm-or-correct decision. The deviation list
may be discarded after materialization; it is not the semantic authority. A
corrected mapping is therefore usable by later algorithms without replaying a
notebook session or interpreting an editing command.

The first cut does not add reviewer accounts, assertion history, confidence
scales, or a general patch language. An item that is passed over remains
unreviewed and can be revisited later. More elaborate ambiguity or multiple-
reviewer semantics return only if the bounded campaign demonstrates the need.

The S6d storage slice in doc 201 does not define physical assessment or
synthesis tables before this annotation workflow exists. Their design is
deferred here so that actual review determines the required targets, evidence
references, treatment of multiple labeled mappings and named equivalence
relations, synthesized conclusions, notes, and read/write boundary. Empty
generic annotation tables are not placeholders for that design.

## Evidence sources

### Evidence roles and scale

The corpus should keep the role of each population explicit. Proved exact
results test optimization and objective behavior where they are tractable, and
their coverage must be reported by molecular size and structural class rather
than generalized to the whole corpus. Operation-issued reaction-network maps
and reconstructed endpoint maps provide scalable conditional references
without requiring a new exact search over the endpoints. External mapper
outputs remain challenge candidates and comparative evidence, not ground
truth. Human review is reserved for stratified samples and disagreement sets;
unreviewed large populations remain unadjudicated rather than silently
receiving manual or external labels.

Automatic evidence is responsible for reducing those populations to reviewable
strata. At minimum, selection should distinguish proved agreement, proved
objective disagreement, unproved search, path-stable and path-dependent network
mappings, aromatic cases, hydrogen or wildcard cases, symmetry multiplicity,
and explicit diagnostic families. Sampling must retain the population counts
behind every stratum so a curated subset is not mistaken for corpus-wide
prevalence.

### External reaction mappings

The saved Rhea/RXNMapper 0.4.0 results are the first population because doc 203
already established their provenance, parse, validity, and small-case exact
baseline. A normalized continuation may add the saved Indigo results, RDT
policies, an independently derived modular-product or clique formulation, and
other properly licensed mapper outputs. Raw label agreement is never the sole
comparison; score, correspondence shape, objective feasibility, molecular-
frame alignment, symmetry, mapping-input alignment, and induced reaction
equivalence are separate observations.

LocalMapper remains paper-only evidence. Its code, model, and generated
artifacts are outside this work.

### Reaction-network generation

umol reaction-network generation supplies atom-mapped reactions by
construction for a fixed rule set, stoichiometry, and application path. This
is valuable evidence unavailable to ordinary product/reactant-only corpora:
the mapping witness is part of the derivation rather than reconstructed after
the fact. The corpus should retain the rule-application witnesses, generated
products, composed endpoint correspondences, and aggregate path support so
that a mapper can be evaluated against both the endpoints and the operations
that connected them. Individual enumerated paths need not be separate stored
rows when a complete QRS graph and its transformation correspondences preserve
the derivation space.

These records are conditional references, not universal mechanistic ground
truth. They answer whether a mapping recovers the known derivation in the
generated system and expose when several derivations produce the same endpoint
reaction.

### Native QRS population

Doc 207 supplies the technical producer for a native quasireaction-subgraph
population: complete bounded network closure, QRS extraction, induced GraphML,
faithful transformation correspondences, heavy-atom mapping classes, support
aggregates, Parquet persistence, and per-job manifests. This removes the need
to reconstruct atom identity from the historical GraphML edges and is already
sufficient to feed the atom-mapping investigation over localized, non-stereo
structures.

The generated QRS population is not a curated collection of chemically
preferred reactions. Its distinctive value is completeness over a declared
finite domain: every state admitted by the rule catalog and molecular
composition participates, including oxidation states and reaction paths that
would be uncommon or deliberately excluded from ordinary reaction databases.
For each admitted network, the corpus-generation operation retains every
eligible QRS rather than selecting endpoints for apparent chemical interest.
Curated reactions remain diagnostics and later review cases; they do not define
the generated population.

Completeness has three separately reported boundaries:

- **network closure** means every state reachable under the declared catalog,
  seed, and closure bounds has been expanded;
- **QRS population completeness** means every eligible endpoint pair in that
  network has been considered under the declared endpoint predicate and path-
  length domain; and
- **path-evidence completeness** means every eligible path within the declared
  length and slack domain was enumerated rather than stopped at a per-QRS path
  bound.

A path-limited QRS may remain a useful corpus record, but its mapping and
support counts are bounded observations rather than complete evidence. These
three statuses must remain queryable independently so corpus size is not
mistaken for complete closure or complete mapping support.

Endpoint mapping evidence is set-valued. Full correspondences that differ only
in explicit-hydrogen placement project to the same heavy-atom map. Heavy maps
related by independent endpoint automorphisms form one mapping class, while
non-symmetry-equivalent classes remain separate even when their path evidence
is equal. Minimum path length and support multiplicity are evidence from the
generated system, not a claim that one class is the unique chemically correct
mapping.

The first population work should revisit the classification-corpus construction
at larger scale rather than stop with doc 207's calibration cases. It includes:

- complete normal-polarity populations across explicit compositions and
  oxidation-state domains, renewed seven- and eight-heavy-atom attempts, and
  positive-slack extraction where the expanded runtime now makes it feasible;
- extended and mixed rule populations that include non-normal-polarity paths,
  beginning with the implemented carbon--hydrogen and carbon--hydrogen--oxygen
  electron-allocation catalogs and then extending the element domain;
- aromatic mapping inputs derived consistently from generated Kekule structures
  through the existing aromaticity model and perception machinery, with the
  localized source and resulting aromatic-system representation related
  explicitly; and
- model-dependent stereo perception as dedicated implementation work, initially
  limited to tetrahedral atom and cis/trans bond stereo and following the
  established SMILES/MOL raise precedent.

#### Rule-catalog construction

The CH, CHO, and CHN electronic-state catalogs use the same **extended** rule
derivation and all three also have **mixed** constructions. Extended and mixed
catalogs do not share a state-domain formula. For a fully defined atom, the
counts invariant is

```text
B + 2n + u + c = V
```

where `B` is the sum of localized bond orders incident on the atom, `n` is its
lone-pair count, `u` its unpaired-electron count, `c` its formal charge, and
`V` the element's valence-electron count. Multiplicity is `s = u + 1` in these
catalogs. Catalog construction states every `c`, `n`, `u`, and `s` value
explicitly; resolution does not supply an assumed `n = 0`.

The extended construction is an exact whitelist of unmixed local electronic
states. It is not the set of tuples satisfying the mixed-state inequality. Its
carbon states are:

| Description | Atom DSL |
| --- | --- |
| ordinary closed shell | `C#c0#n0#u0#s` |
| cation | `C#c+#n0#u0#s` |
| anion | `C#c-#n1#u0#s` |
| neutral doublet radical | `C#c0#n0#u1#s2` |
| neutral singlet carbene | `C#c0#n1#u0#s` |
| neutral triplet carbene | `C#c0#n0#u2#s3` |

Its oxygen states are independently defined from oxygen's electron budget;
they are not obtained by adding one lone pair to every carbon state:

| Description | Atom DSL |
| --- | --- |
| ordinary closed shell | `O#c0#n2#u0#s` |
| cation | `O#c+#n1#u0#s` |
| anion | `O#c-#n3#u0#s` |
| neutral doublet radical | `O#c0#n2#u1#s2` |
| neutral singlet atom | `O#c0#n3#u0#s` |
| neutral triplet atom | `O#c0#n2#u2#s3` |

The extended nitrogen whitelist is:

| Description | Atom DSL |
| --- | --- |
| ordinary closed shell | `N#c0#n1#u0#s` |
| cation | `N#c+#n0#u0#s` |
| anion | `N#c-#n2#u0#s` |
| neutral doublet radical | `N#c0#n1#u1#s2` |
| neutral singlet nitrene | `N#c0#n2#u0#s` |
| neutral triplet nitrene | `N#c0#n1#u2#s3` |

Hydrogen retains the four states `H#c0#n0#u0#s`, `H#c+#n0#u0#s`,
`H#c-#n1#u0#s`, and `H#c0#n0#u1#s2`. The extended domain does not combine
charge with radical or singlet/triplet character; those combinations belong to
the mixed construction.

The mixed construction retains the same four hydrogen states and enumerates
carbon, oxygen, and nitrogen independently. Carbon has `|c| <= 2`, `u <= 2`,
`0 <= n <= 1`, and `|c| + u + n <= 2`. Oxygen has `|c| <= 2`, `u <= 2`,
`1 <= n <= 3`, and `|c| + u + (n - 2) <= 2`. Nitrogen has `|c| <= 2`,
`u <= 2`, `0 <= n <= 2`, and `|c| + u + (n - 1) <= 2`. The adjusted
lone-pair terms are literal differences, not clamped at zero. All three use
`s = u + 1`.

The inequality candidates are then filtered by their implied localized bond
budget `B = V - 2n - u - c`: retain only `B >= 0` and `B + n + u <= 4`.
This leaves ten carbon states, thirteen oxygen states, and fifteen nitrogen
states while retaining every extended state. The construction admits mixed
charged, radical, singlet, and triplet states without admitting second-row
states that require more than four occupied bonding, lone-pair, and singly
occupied orbitals. Relative to the extended nitrogen whitelist, mixed nitrogen
adds:

| Atom DSL | `B` |
| --- | ---: |
| `N#c+2#n0#u0#s` | 3 |
| `N#c+#n0#u1#s2` | 3 |
| `N#c+2#n0#u1#s2` | 2 |
| `N#c+#n0#u2#s3` | 2 |
| `N#c+#n1#u0#s` | 2 |
| `N#c+2#n1#u0#s` | 1 |
| `N#c+#n1#u1#s2` | 1 |
| `N#c+#n2#u0#s` | 0 |
| `N#c0#n2#u1#s2` | 0 |

For CH and CHO, enumerate C--C bonds through order three, C--H bonds at order
one, C--O bonds through order three, O--H bonds at order one, and O--O bonds
through order two. Extended and mixed CHN retain the same C/H families and add
C--N and N--N bonds through order three and N--H bonds at order one. A forward rule
decreases one localized bond order by one and allocates the released electron
pair as `2+0`, `0+2`, or `1+1`. Retain the decrement only when both resulting
endpoint states are in the selected domain and each endpoint obeys
`c_rhs = c_lhs + 1 - allocated_electrons`. Pair every retained decrement with
the exact structural inverse that increases the same bond order.

The extended and mixed catalogs also contain topology-neutral one- and
two-electron transfers. For an atom state define its nonbonding-electron count
as `e = 2n + u`. A `k`-electron donor transition, for `k` in `{1, 2}`, retains
the atom's bond budget and obeys `e_rhs = e_lhs - k` and
`c_rhs = c_lhs + k`; its reverse is an acceptor transition. Enumerate pairs of
complementary donor and acceptor transitions over the same element-pair
families as the bond rules, still excluding H--H. Homoelement families use
unordered transition pairs with repetition. Emit every resulting transfer and
its exact inverse.

The implemented catalog populations are:

| Catalog | Bond decrements | Bond entries | Redox pairs | Redox entries | Total entries |
| --- | ---: | ---: | ---: | ---: | ---: |
| extended CH | 54 | 108 | 9 | 18 | 126 |
| mixed CH | 208 | 416 | 42 | 84 | 500 |
| extended CHO | 144 | 288 | 14 | 28 | 316 |
| mixed CHO | 1,082 | 2,164 | 222 | 444 | 2,608 |
| extended CHN | 154 | 308 | 14 | 28 | 336 |
| mixed CHN | 1,603 | 3,206 | 241 | 482 | 3,688 |

Each larger catalog retains its corresponding implemented smaller-domain rules
exactly. The executable catalog audit derives the domains and redox transitions
independently and checks state membership, bond and orbital budgets, charge and
electron accounting, family completeness, name uniqueness, and reciprocal
graph edits.

The redox rule left-hand sides contain two atoms and no localized bond. The
current non-induced matcher can match such a pattern to any compatible pair of
distinct atoms, including a bonded pair or two atoms in the same connected
component. Constraining electron-transfer matching by bond absence or component
membership is separate matching and operation work; these catalogs do not
claim either constraint. The catalogs also do not by themselves settle a
population campaign or establish network, QRS-population, or path-evidence
completeness.

The first CHN catalog extends the historical normal-polarity model rather
than the extended electronic-state model. Its rules track formal-charge
polarity only: C--H cleavage assigns the pair to carbon, C--N and N--H cleavage
assign it to nitrogen, and the homonuclear C--C and N--N families retain one
canonical heterolytic orientation. The catalog covers C--C, C--N, and N--N
through triple bonds and C--H and N--H single bonds. It retains the eight
existing normal-polarity C/H entries exactly and adds eleven nitrogen-containing
forward decrements with their inverses, giving 15 forward decrements and 30
entries in `normal-polarity-nitrogen-rules.edn`. Its executable audit checks the
exact family/order/charge-transition set, total-charge conservation, one-order
bond changes, unique reciprocal links, structural inversion, and exact C/H
retention. Direct release-mode checks completed for the charge-only N2 and HCN
seeds: N2 produced 2 flasks and 3
transformations, while HCN produced 4 flasks and 6 transformations. Both
reversibility checks reported no missing inverse application.

The extended CHN check used a fully ground `H2C2N2` seed consisting of two
disconnected H--C triple-bond N molecules. Generation reached closure at 4,394
flasks and 48,749 transformations. The counts-based valence-invariant validator
accepted all 4,394 flasks, with no underdetermined or contradictory result.
Applying every transformation's declared inverse recovered its source in all
48,749 cases: 43,868 by identity and 4,881 through a nonidentity source
automorphism. These closures are operational witnesses, not population-
completeness evidence.

The latter expansions are active atom-mapping corpus work, not exclusions from
the investigation. They do not belong to doc 207 because they change the
generated chemical population or require new chemistry operations rather than
the reaction-network producer's storage and export boundary. Before this
document becomes an implementation plan, it must settle the concrete network
population, rule-catalog extensions, completeness requirements, artifact
partitioning and distribution, and the aromatic and stereo construction
contracts.

### Existing quasireaction subgraphs

The read-only corpus at
`/Users/dr/Projects/chemical-network/classification` contains an independent
source of small network-derived cases: 184 directed CHO networks and 194,778
nonempty reducible endpoint subgraphs with two to six non-hydrogen atoms and
shortest paths of length two, four, six, or eight. Its 32,260 Weisfeiler--Lehman
buckets are useful for stratified sampling, but they are not exact isomorphism
classes.

The archived subgraphs alone do not preserve a unique full atom mapping. Their
edges retain rule strings, while construction converted the directed network
to undirected endpoint subgraphs and did not retain the particular host match.
A usable conditional reference therefore requires returning to the original
directed network, applying the recorded rule to the source, retaining every
application that produces the target, and composing the resulting
correspondences along a path. The corpus can then distinguish path-stable,
path-dependent, and locally ambiguous endpoint maps rather than manufacturing
one answer.

This source may be read for the investigation but must not be modified.

### Curated diagnostic reactions

Small explicit cases remain necessary because large aggregate counts cannot
explain why an objective succeeds or fails. The initial set includes
`Rh_20309`, `Rh_63116`, `Rh_12116`, chemically dubious same-element
permutations such as `Rh_34527`, `Rh_10012`, `Rh_10044`, and the parent
Diels--Alder reaction. `Rh_10012` distinguishes the external oxygen provenance
from the lower structural edit score, while `Rh_10044` exposes aromatic versus
localized representation in both depiction and scoring. Each should be
represented as an ordinary corpus record with full mappings and visual
witnesses, not as prose disconnected from executable data.

## Query and visualization requirements

The columnar form should make the full evidence directly inspectable with
DuckDB, Polars, or equivalent tools before umol defines native chemical query
operators. Useful projections include:

- cases grouped by source, mapping input, preprocessing path, size, element
  composition, and exact-proof status;
- score gaps and component-wise disagreements between algorithms;
- labeled-optimum and equivalence-class multiplicities;
- runtime distributions, timeouts, search-node counts, and hard-tail slices;
- cases where confidence, structural optimality, and curated judgment disagree;
- network cases grouped by rule, path length, and mapping stability.

The notebook view should render lhs and rhs structures, atom correspondences,
bond changes, alternative maps, symmetry classes, and provenance together.
Static SVG is the baseline export. Structured element identifiers in the visual
payload should permit linked selection between a table row, atoms and bonds in
the drawing, score components, and the originating network or external record.

This is a scientific debugging interface, not a hand-picked abridged summary.
A user must be able to filter to a named evidence stratum, inspect any selected
member, record or revise a judgment with its evidence, and recover the machine-
readable records behind the drawing. The interface must make a bounded campaign
efficient; it is not an invitation to page manually through the complete
corpus.

The S6 notebook's source-id environment variable proves that one known case can
be reconstructed; it is not the annotation workflow. The review interface
opens a persistent, ordered review set and supports at least first-unreviewed,
previous, next, direct source-id lookup, evidence-stratum filters, and reviewed/
unreviewed progress counts. Saving a confirmation or correction persists it
immediately and advances to the next unreviewed item. Navigation without saving
leaves the item unreviewed rather than creating an implicit judgment.

## Literature and formulation gate

The useful size frontier of a strengthened native branch-and-bound search is
not yet known. It must not be inferred from the improvement list below or from
performance on the current small exact slice. Before selecting the search
implementation, review the mathematically explicit non-learned literature and
the available local implementations as algorithm-design evidence.

The review must distinguish four roles:

- an exact formulation with the same feasible set and objective can certify an
  optimum independently of the current implementation;
- a stronger exact native search can become another implementation of the same
  operation;
- an approximate or differently scored method can supply a checked incumbent
  or value ordering without supplying a proof; and
- a chemistry rule, reaction-center restriction, or different objective
  defines different mapping semantics rather than an optimization.

The initial literature set is:

- Heinonen et al.'s exact
  [A* minimum-edge-edit search](https://doi.org/10.1089/cmb.2009.0216), with
  particular attention to its admissible heuristics;
- Latendresse et al.'s
  [minimum weighted edit-distance MILP](https://doi.org/10.1021/ci3002217),
  which proved 87% of 7,501 MetaCyc instances within ten seconds under its own
  weighted objective;
- Flamm et al.'s
  [AltCyc and ILP2](https://doi.org/10.1007/978-3-319-40530-8_13)
  formulations,
  including their Rhea measurements. ILP2 systematically outperformed the
  search-tree method, but both showed a sharp easy/hard division;
- McSplit's compact bidomain partitioning and bounds from
  [McCreesh et al.](https://doi.org/10.24963/ijcai.2017/99). These are relevant
  branch-and-bound techniques, not a direct solution: induced maximum common
  subgraph preserves nonedges, while atom mapping must admit bond additions and
  deletions under a maximum-cardinality correspondence;
- the MCES/modular-product and maximum-clique construction in
  [AAM-Ising](https://doi.org/10.1021/acs.jcim.4c01871) and its local source.
  The preserved-edge clique, completion, symmetry reduction, hydrogen filter,
  and bond-order filter must be separated; the resulting lexicographic
  semantics are not assumed to equal A0 or A1;
- the recent partial-map results of
  [Gonzalez Laffitte et al.](https://doi.org/10.4230/LIPIcs.WABI.2025.12),
  which reduce completion to constrained isomorphism after a good partial map
  covers the reaction center. This is directly useful for supplied anchors and
  manual corrections, not for discovering the center from endpoints; and
- the 2025 SLAPMapper
  [preprint](https://doi.org/10.26434/chemrxiv-2025-hthwn) and
  [MIT-licensed source](https://github.com/shin1koda/slap-mapper). Its
  Weisfeiler--Lehman-like sequential linear assignments are an approximate,
  polynomial-time source of candidates and incumbents, not an optimality
  certificate.

Reaction Decoder Tool, ReactionMap, Indigo, and other MCS or chemistry-rule
systems remain useful for decomposition, candidate generation, ring handling,
and failure cases. Their output quality does not make their internal selection
an exact implementation of A0 or A1. Learned methods remain comparative
evidence rather than search-design inputs for this pass.

For each family, the review records the feasible set, objective, exact or
heuristic status, available bound or certificate, treatment of disconnected
and unbalanced inputs, aromaticity, symmetry and multiple optima, reported
scale, and reusable algorithmic ideas. The result is a semantic translation,
not a mapper leaderboard.

The old branch-and-bound implementation is not the only reference for larger
cases. Use an independent certificate ladder:

1. exhaustive enumeration and the current implementation agree on tiny cases;
2. a development-only MIP, CP-SAT, pseudo-Boolean, or equivalent formulation
   is translated exactly to A0 and then A1 and cross-checked on the shared
   tractable range;
3. on larger cases, a solver's matching bounds and completed proof establish
   the reference result even when the old implementation cannot finish; and
4. where no exact formulation closes the gap, report the incumbent and bound
   gap as unresolved. A network witness or external map can establish recovery
   of that correspondence but not optimality under the named objective.

The classification corpus is the cleaner search-engineering population because
its operation-issued mappings provide independent recovery evidence. Its
current `subgraphs.csv` verifies the no-stereo assumption but not the all-
Kekule assumption: none of 642,120 endpoint rows contains a SMILES `@`, `/`, or
`\` stereo marker, while 29,326 rows contain lowercase aromatic endpoint
notation, including 5,605 rows marked valid. Begin with the Kekule, non-stereo
stratum and retain the aromatic classification cases as a separate objective
and notation stratum. Rhea remains the necessary realistic-ingest population.
Proof coverage must be stratified by atom count, compatible-domain size,
components, symmetry, aromatic status, and objective version so that either
population cannot hide a hard class in the other.

This gate ends with a justified choice among a strengthened direct search, a
native alternative formulation, or a small exact portfolio. It also identifies
the independent certificate formulation and incumbent generator used during
development. No staged implementation plan should assume that the direct
branch-and-bound improvements alone cross the useful Rhea size frontier.

## Pre-annotation algorithm pass

After the literature and formulation gate selects the search approach, the
known algorithmic improvements should be implemented as one coherent pass
before the main manual campaign. This avoids spending the fixed review budget
on cases selected by an obviously temporary tractability frontier and avoids
asking reviewers to revisit small successive algorithm iterations.

Implementation variants may still be compared automatically while this pass
is developed. The restriction is on consuming the manual-review budget after
each variant, not on evidence-driven algorithm engineering.

Optimization and complete optimum enumeration should become separable
operations. Once the optimum score has been proved, an external or network-
witnessed mapping can be evaluated against that score without enumerating every
labeled optimum. Complete labeled enumeration remains available when the
mapping set itself is the required result, and proof status remains distinct
from enumeration completeness.

The fixed-semantics implementation pass includes:

- compute maximum compatible cardinality once rather than rediscovering it
  throughout the search;
- seed a valid incumbent from a maximum matching;
- order values by immediate edge reward without removing candidates;
- maintain stronger matching feasibility during assignment;
- strengthen the admissible remaining-edge bound;
- support a bounded corpus run that reports completed optimization, completed
  optimization with interrupted enumeration, or interrupted optimization
  without mislabeling the incumbent as an optimum; and
- instrument search nodes, domain reductions, bounds, incumbent changes, proof
  time, enumeration time, and output multiplicity.

Component decomposition and exact symmetry reduction with complete orbit
expansion remain candidates after the simpler evidence is available.
Development-only MIP, CP-SAT, partial-QAP, or weighted-clique formulations can
provide independent certificates and performance comparisons. A formulation
becomes another backend for doc 203's operation only if it preserves the exact
feasible set, objective, and all-labeled-optima contract.

The principal implementation result is increased proof coverage over a fixed,
stratified corpus, not a storage benchmark or a faster timing on already
trivial reactions. The corpus must report two results separately: the
population for which global optimality is established under fixed semantics,
and recovery of independently witnessed mappings on the larger network-derived
population. Improved runtime without increased proof coverage is not an exact-
baseline improvement, while agreement with a network witness does not by
itself prove that the witness minimizes the named structural objective.

## Objective investigation

Known objective defects must also be addressed, or assigned explicit exclusion
strata, before the main annotation campaign. The corpus should expose objective
components before selecting a successor to the A0 structural baseline.
Candidate observations include element-form preservation, explicit-hydrogen
changes, atom charge, radical and valence changes, bond-edit categories,
locality, and supplied partial maps. These are measurements first, not approved
weights or tie-breakers.

The pre-annotation design pass should select the objective changes justified by
the existing diagnostic cases and automatic evidence, then implement evaluation
and exact search for that named objective together. The review campaign records
human evidence independently of the selected objective, while every automatic
optimality claim names its objective version. Further objective research
remains possible, but the campaign will not be organized as repeated annotate--
retune cycles.

The key question is whether a proposed criterion merely resolves degeneracy
among A0 optima or must admit a mapping outside the A0 optimum set. The
wildcard case and parent Diels--Alder reaction already show why that distinction
matters: a selector over current minima cannot recover a supported map that the
current objective has excluded. Hard anchors, a revised objective, and a
mechanism- or rule-conditioned algorithm are therefore separate hypotheses to
test.

### Aromatic-system bootstrap

The first Rhea inspection exposes a representation-dependent A0 score. The
reaction-SMILES/TableIR path raises an aromatic bond as localized order one
plus a definite aromatic constraint, while A0 compares only localized bond
order and ignores the constraint and any aromatic-system overlay. A mapped
aromatic bond paired with a localized double bond is therefore counted as a
modification, whereas the same aromatic bond paired with a localized single
bond is not. Excluding every aromatic reaction would make the first comparison
unrepresentative; discounting aromatic transformations entirely would make the
objective uninformative on an important part of the corpus.

Call the first aggregate successor **A1**. A1 operates only after ingest and
resolution have produced supported, determined aromatic components. Ambiguous,
contradictory, or unsupported aromatic notation is an input-status stratum,
not a scoring fallback: A1 must report the whole mapping input as outside its
declared domain rather than silently evaluate that part with A0. The F420 case
from docs 174 and 194 demonstrates this distinction. Its fused aromatic system
can be resolved under a stated selection policy, but the source notation does
not uniquely determine the tautomer; an atom-mapping objective must not hide
that boundary decision.

The first aggregate approximation to test treats each supported aromatic
component as a reservoir of mobile double bonds. Under an initial bipartite,
perfect-matching assumption, define

```text
k(A) = size of a maximum matching of aromatic component A = |V(A)| / 2
```

This equals `ceil(|E(A)| / 2)` only for an isolated even cycle. Bipartiteness
alone does not justify the edge-count expression for fused systems: a
naphthalene graph has eleven aromatic edges but a maximum matching of size
five. The existing Hopcroft--Karp and Kekulizer machinery can calculate the
matching quantity directly.

For a candidate that maps the complete aromatic component onto the same bond
topology represented by localized orders on the other side, let `d` be the
number of order-two bonds in that image and contribute `abs(k(A) - d)` aromatic
modifications. Two mapped aromatic components contribute the difference of
their reservoir sizes, `abs(k(A_lhs) - k(A_rhs))`; ordinary localized bonds
retain the A0 order comparison.
Topology additions and deletions remain separate. Component splitting,
merging, partial component mapping, a missing perfect matching, and unsupported
localized orders are explicit exclusions from this first approximation rather
than silently assigned a score.

The aggregate contribution does not identify particular aromatic bonds as
modified and should be stored separately from ordinary localized-order
modifications. Under this proposal, `Rh_10044` contributes one aromatic
modification from the nicotinamide system and one ordinary alcohol-to-carbonyl
modification, rather than the current three localized-order modifications.

This score cannot be implemented by the current independent per-edge
equivalence callback: it depends on a mapped component and an aggregate count.
Evaluation and search must share the same versioned objective before an exact
result can be compared with externally supplied candidates. Implementing the
aggregate only as a report-time correction would create incomparable scores and
is not sufficient.

## Boundaries

This document does not yet:

- define a durable corpus schema or stable identifier;
- select a durable columnar layout, corpus distribution, or notebook frontend
  beyond the experimental path in doc 201;
- define a configurable public atom-mapping objective;
- select chemistry-derived weights, anchors, or narrowing rules;
- make an external mapper or solver a runtime dependency;
- equate network-derived mappings with experimental mechanisms;
- distribute corpus generation, annotation, or query tooling as part of `umol`
  or `umol-py`;
- schedule implementation stages;
- include geometry, 3D export, ontological services, or the broader template-
  library metadata model.

Doc 201 has supplied the first faithful records, columnar path, direct queries,
and notebook bridge. Before this document becomes an implementation plan, it
must settle the independently regenerable source-corpus and algorithm-run
boundary, annotation targets and evidence references, review and sampling
workflow, and the complete pre-annotation algorithm and objective pass. The
physical and public surfaces from doc 201 remain provisional while that
experience is collected.
