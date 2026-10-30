# 242 — Stereo addition, removal, and enumeration

Status: In Progress
Date: 2026-10-29
Relates: [166](166-molecule-ops-2026-07-27.md),
[240](240-validation-constraint-ops-2026-10-18.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md),
[Python API guide](../docs/development/python-api.md)

## Scope

Add StereoAdder, StereoRemover, and StereoEnumerator as separate Transformer
implementations. They share StereoPerception and the existing graph-symmetry
and coset machinery. Identification and configuration enumeration are separate
tasks: addition creates missing stereo entities; enumeration determines selected
configuration values on existing entities.

This record contains the agreed semantics, interfaces, and nomenclature.
S0's fixtures/references/baselines, S1's binary classification and
ligand-construction corrections, S2's selectors/configurations/shared queries,
S3's three transformer engines and absolute-assignment traversal, S4a's
properties, S4b's receiver methods, S4c's Python values and existing perception
queries, and S4d's Python operations are implemented. S5's documentation and
performance closeout remain future work.
External implementations under materials/codes may supply
algorithms, counterexamples, and comparison cases;
their public API organization does not determine umol's operation boundaries.

Variables and variable-based relative configurations are deferred. Initial
enumeration handles undetermined cosets, literal cosets, and literal coset sets.
The required discovery, classification correction, and enumeration algorithms
are tetrahedral and cis/trans only. Other geometries do not expand this work;
unrestricted removal can still strip existing records without interpreting them.
It performs no symmetry deduplication, canonical selection, geometry embedding,
or energetic ranking. Removing kindless forms from the IR is outside this work.
Retain the current para-stereo implementation and its limitations. Improving
para-stereo accuracy or convergence is follow-up work, not a prerequisite.

## Representation and shared perception

A stereo atom or stereo bond records an arrangement in an ordered ligand frame.
Its presence does not assert stereogenicity. Symmetric and prochiral entities
remain legitimate representations. The three operations use the existing
Stereogenicity values: Symmetric, Prochiral, and Stereogenic. Selection by those
values is an operation policy, not an additional stored assertion.

Keep the settled frame contract from
[214](214-aggregate-frame-semantics-2026-08-28.md): StereoLigand values
are pairwise-distinct, and frames are bounded by MAX_DEGREE. Equal virtual-H
occurrences do not acquire separate positions or identities. A use case needing
distinguishable hydrogen positions uses explicit atoms. Addition does not
implicitly unfold hydrogens or weaken frame integrity.

StereoPerception currently materializes positive atom/bond stereo evidence and
checks existing entities. Extend the same object with the identification and
classification required by these transformers. Preserve the evidence-driven
resolver contract: extending discovery must not make ordinary stereo resolution
start adding every unasserted candidate. Recovery policies and edit construction
remain outside perception; graph symmetry and coset algebra retain their
graph-IR and umol-perm owners.

The shared read-only queries are public methods on StereoPerception:

    pub fn recognize(
        &self,
        molecule: &Molecule,
        selection: &StereoSelection,
        automorphism_algorithm: AutomorphismAlgorithm,
        ring_config: RingConfig,
    ) -> Result<StereoDerivation, StereoPerceptionError>;

    pub fn classify(
        &self,
        molecule: &Molecule,
        sites: &BTreeSet<StereoSite>,
        automorphism_algorithm: AutomorphismAlgorithm,
    ) -> Result<BTreeMap<StereoSite, Stereogenicity>, StereoPerceptionError>;

recognize returns only missing selected entities. Positive evidence supplies
their cosets; otherwise the coset is Undetermined. Explicit negative evidence
excludes a site. Automatic chemistry eligibility applies to unasserted sites;
positive assertions use the existing structural frame rules. Existing entities
are retained, and contradictory evidence at their sites is reported using
StereoInconsistency. Recognition reads open evidence, so the absence of an
overlay does not imply a negative assertion.

When stereogenicities is restricted, recognition appends prospective entities
to a temporary molecule sharing the input's storage, calls classify in that
molecular context, and filters the additions. With an unrestricted class
selector it does no symmetry calculation. Empty selectors select nothing;
explicit topology ids are checked even when another selector is empty.

classify requires an existing entity at every supplied topology site. It
calculates one whole-molecule GraphSymmetry and a local projection per supplied
site. It does not impose automatic recognition eligibility or the model's
element scopes on stored arrangements; classification also works when their
producer kinds are disabled. It supports Th atoms and Ct bonds with Undetermined,
literal zero or one, or a nonempty finite subset of those values. Kindless
configurations, symbolic cosets, invalid coset domains, incompatible site/kind
combinations, and other geometries report unsupported input. Empty input returns
an empty map. Both queries retain the current para-stereo behavior and leave the
input unchanged.

StereoPerception's existing new, derive, derive_stereo_atom, derive_stereo_bond,
and skipped_stereo_bonds signatures remain unchanged. derive stays driven by
source evidence for resolver and validator callers. The frame builders remain
structural queries rather than automatic chemistry screens. No Molecule method,
perception configuration type, or prepared-result type is added in S2e.

StereoPerceptionError distinguishes InvalidSite { site }, MissingStereoEntity
{ site }, UnsupportedKind { kind }, UnsupportedConfiguration { site },
UnsupportedFluxionality { kind }, UnderdeterminedAtomField { atom, field }, and
UnderdeterminedBondOrder { bond }; field is the existing AtomFieldKind.
Inconsistency wraps the existing StereoInconsistency for incompatible stored
frame arity. Evidence failures remain in StereoDerivation.inconsistencies.
Recognition checks explicitly selected kinds, or enabled producer kinds when
the kind selector is unrestricted; only Th/Ct producers are implemented.
Fluxional processing is unsupported rather than silently ignored. Required
literal fields are checked only within the selected processing scope, without
requiring the whole molecule to be concrete.

Correct StereoSymmetry::is_stereogenic and its dependent classification for a
known Th/Ct kind with an undetermined coset. Both have two cosets: the carrier
is stereogenic when the local ligand-symmetry action does not identify those
two arrangements. Determine that action without requiring a stored literal
coset. This is distinct from an undetermined stereo kind. Preserve meaningful
classification of explicit prochiral and symmetric arrangements. No non-binary
classification redesign is required here.

Construct ligand lists from actual neighbors and inherent hydrogen/lone-pair
counts as specified by
[240, S13c](240-validation-constraint-ops-2026-10-18.md#s13--invariant-validation-and-stereo-ligand-construction).
That producer correction is an explicit dependency; the broader validators and
constraint transformations in 240 are independent work.

## Selection and configuration

The shared selection vocabulary is:

    pub enum StereoSite {
        Atom(AtomId),
        Bond(BondId),
    }

    pub struct StereoSelection {
        pub sites: Option<BTreeSet<StereoSite>>,
        pub kinds: Option<BTreeSet<StereoKind>>,
        pub stereogenicities: Option<BTreeSet<Stereogenicity>>,
    }

    impl StereoSelection {
        pub fn matches(
            &self,
            site: StereoSite,
            kind: Option<StereoKind>,
            stereogenicity: Option<Stereogenicity>,
        ) -> bool;
    }

    pub struct AddStereoConfig {
        pub selection: StereoSelection,
        pub automorphism_algorithm: AutomorphismAlgorithm,
        pub ring_config: RingConfig,
    }

    pub struct RemoveStereoConfig {
        pub selection: StereoSelection,
        pub mode: RemoveStereoMode,
        pub automorphism_algorithm: AutomorphismAlgorithm,
    }

    pub struct EnumerateStereoConfig {
        pub selection: StereoSelection,
        pub automorphism_algorithm: AutomorphismAlgorithm,
    }

None means unrestricted; an empty set selects nothing. Site, kind, and
stereogenicity restrictions intersect. Site ids refer to the input molecule's
topology, not its stereo-entity tables.
StereoSelection is an open configuration carrier. Its construction does not bind
ids to a molecule; the executing operation checks supplied ids before indexing
or mutation. Selection limits changes, while classification uses the molecular
context, including stereo outside the selected sites.

StereoSelection::default leaves all three selectors unrestricted. The matches
predicate tests only the supplied values; it neither checks ids against a
molecule nor calculates stereo. A missing kind or stereogenicity satisfies
only an unrestricted selector for that field.

Each operation uses one StereoSelection. Its configuration applies uniformly to
the selected sites, kinds, and stereogenicity classes. RemoveStereoConfig also
has one All/Undetermined mode for that call.

AddStereoConfig, RemoveStereoConfig, and EnumerateStereoConfig own their shared
selection and the operation's algorithm settings. Algorithm choice remains an
operation configuration, separate from StereoModel's chemistry choices. Keep
the existing model and supported-kind vocabulary. A requested producer or
classification that is unsupported must report that boundary, not silently skip
the selected kind. The adder's default kind selector contains Tetrahedral and
CisTrans. Initial discovery/classification/enumeration support is Th/Ct; a
selected unsupported algorithm reports that fact. Unrestricted All removal
requires no geometry-specific algorithm. Algorithm defaults use Nauty and, for
addition, RingConfig::default(), matching the existing operation defaults.
Enumeration's default selection is Th/Ct, all sites, and all stereogenicities.

Defaults and operation-specific settings are:

| Operation | Default | Additional setting |
| --- | --- | --- |
| Addition | Tetrahedral and CisTrans; Stereogenic only; all sites | Shared site/kind/stereogenicity selection |
| Removal | All configurations; all sites, kinds, and stereogenicities | All or Undetermined intersected with the same selectors |
| Enumeration | Selected existing open absolute configurations | Caller-supplied site order; no deduplication |

RemoveStereoMode has the variants All and Undetermined. Undetermined refers to
the coset, not to an unknown kind. Keep that distinction in the implementation,
errors, and tests; do not rename the variant UndeterminedCoset. A failed as_lit
query is not by itself the Undetermined predicate: literal sets and symbolic
forms retain their own meanings.

RemoveStereoMode::matches takes a StereoConfigurationForm and tests the mode
without calculating stereogenicity or resolving a coset. All accepts every
configuration; Undetermined accepts only a known kind with an explicit
StereoCoset::Undetermined. A configuration with an undetermined kind has no
coset and does not satisfy that mode.

Unrestricted All removal also handles kindless records without attempting a
symmetry query that requires a known kind. An explicit kind set selects only
records of those known kinds. Selection must not guess a kindless record's kind
or classify it by the open-coset fix.

## Transformer contracts

All three engines implement the existing Transformer trait. They plan every
rejection before publishing a result. transform consumes the input;
transform_into applies the same choice in place and preserves its entry state on
rejection. transform_iter borrows its source and yields independent molecules.
Retain the trait's current error contract: singular methods return their error;
the iterator yields molecules and omits rejected results. Enumeration reports
unsupported selected variable/relative configurations on its singular methods;
it does not treat them as independent unconstrained cosets. Unconditional
removal need not evaluate those configurations merely to remove them.

### StereoAdder

Add selected missing entities using shared perception. Preserve existing
entities, their frames, and configurations. Honor existing positive and negative
site evidence. A new entity receives its known kind and perceived ligand frame;
its coset remains undetermined unless evidence supplies a literal or restricted
value. Do not choose a determined isomer during addition.

New entities must have ligand lists agreeing with the site's neighbors and
hydrogen/lone-pair counts and satisfy representation integrity. Report incompatible
selected evidence or competing incompatible kind assignments before mutation.
Addition retains the input's #T/#C and unrelated assertions; its edit plan only
appends missing entities. AddStereoError::Perception wraps StereoPerceptionError;
AddStereoError::Inconsistency rejects a reported StereoInconsistency, with the
first diagnostic returned when perception reports several.
Successful addition is idempotent for the
same model and configuration. Its iterator produces its single transformed
result, or no result on rejection.

### StereoRemover

Remove only entities selected by StereoSelection and the All/Undetermined mode.
Non-stereogenic-only removal is expressed by selecting Symmetric and
Prochiral; a user may retain Prochiral by selecting Symmetric alone. Validate
explicit topology ids, then apply mode, site, and kind filters before requesting
classification. A class filter uses StereoPerception::classify once in the
unchanged whole molecular context. With no class filter, both modes omit
classification: All accepts any configuration, while Undetermined matches only
an explicitly undetermined coset of a known kind. Kindless, symbolic, and
other-geometry records can therefore be removed without interpretation where
the mode and selectors admit them.

Follow Kekulizer's assertion-removal behavior: removing a selected stereo entity
also clears its corresponding atom #T or bond #C assertion, with no extra switch.
Keep unrelated assertions everywhere, and retain #T/#C hints at sites where no
entity is removed. Use the existing graph-IR overlay removal and
constraint-compaction machinery for references to removed entities.

Removal preserves molecular topology and the retained entities' data, subject
to ordinary id compaction and reference transport. RemoveStereoError::Perception
wraps StereoPerceptionError. Invalid selection ids and failures of required
classification reject before publication; unsupported records excluded by the
earlier filters do not themselves require classification. Its iterator produces
its single transformed result, or no result on rejection. Without a class
filter, removal is idempotent for the same config; class-filtered removal can
change the context used by a later call.

### StereoEnumerator

Operate on existing entities. Preserve their sites, kinds, ligand frames, fixed
configurations, and unselected data. Enumerate selected undetermined cosets over
the kind's finite domain and selected literal sets over their admitted members.
Check supported input constraints; selecting a variable or a relative symbolic
configuration reports unsupported input. Do not infer variables' correlations
or treat expressions as independent free choices.

Check #T/#C through graph-IR evidence readings, retaining their transport from
the stored ligand frames. Check Th/Ct stereogenicity, topicity, and ligand-symmetry
assertions against molecular and local symmetry, including assertions on retained
entities. Evaluate affected logical constraints with other leaves handled by
the existing constraint evaluator. Retained open fields may leave a claim
undecided; they need not become determined. Reject a decided contradiction.
Unrelated standalone constraints are preserved without reevaluation.

Choose one result as follows:

1. Visit selected sites in the caller-supplied order.
2. Try admitted coset indices in ascending order in their stored ligand frames.
3. Reject inconsistent candidates and return the first complete consistent
   assignment, backtracking as needed for supported constraints.

transform and transform_into make that same choice. transform_iter emits every
consistent assignment in that order, including symmetry-equivalent results.
The first iterator result agrees with a successful singular result. If no
assignment survives, the singular methods reject and the iterator is empty.
If no selected field needs determination, a consistent input yields one unchanged
result. The order controls a reproducible choice in the supplied frame; it does
not promise invariance under renumbering or ligand reframing.

Validate the selected-site order before publication: all selected processing
sites must be covered without duplicates. As with Kekulizer's atom order,
selection and processing order are separate inputs; a selection set does not
silently supply a canonical order.

## Public interface inventory

| Public symbol | Role and boundary |
| --- | --- |
| StereoSite | Typed atom/bond topology site; no molecule provenance |
| StereoSelection | Shared open site/kind/stereogenicity selector |
| StereoSelection::default, StereoSelection::matches | Unrestricted default and intersection predicate over supplied values; no molecule lookup |
| StereoPerception::recognize | Public query for selected missing Th/Ct entities; preserves evidence and input; returns StereoDerivation |
| StereoPerception::classify | Public query for existing selected Th/Ct arrangements in the whole molecular context; independent of automatic chemistry eligibility |
| StereoPerceptionError | Shared query failures; reuses StereoSite, AtomFieldKind, StereoKind, and StereoInconsistency |
| AddStereoConfig, RemoveStereoConfig, EnumerateStereoConfig | Shared selection and explicit algorithm choices; removal also has one mode |
| AddStereoConfig::default | All sites; Th/Ct kinds; Stereogenic only; Nauty and RingConfig::default() |
| RemoveStereoMode | All and Undetermined variants; Undetermined tests the coset |
| RemoveStereoMode::matches | Tests the configuration against the mode without classification or resolution |
| RemoveStereoConfig::default | All; unrestricted sites/kinds/stereogenicities; Nauty |
| EnumerateStereoConfig::default | Th/Ct; all sites and stereogenicities; Nauty; processing order supplied separately |
| StereoAdder, StereoRemover, StereoEnumerator | Configured Transformer engines sharing StereoPerception |
| StereoAdder::new, StereoAdder::with_config | Infallible model snapshot and owned configuration; topology ids checked during execution |
| StereoAdder's Transformer implementation | Consuming, in-place, and lazy one-result operations; AddStereoError on singular rejection, no iterator result on rejection |
| StereoRemover::new, StereoRemover::with_config | Infallible model snapshot and owned configuration; topology ids checked during execution |
| StereoRemover's Transformer implementation | Consuming, in-place, and lazy one-result operations; RemoveStereoError on singular rejection, no iterator result on rejection |
| StereoEnumerator::new, StereoEnumerator::with_config | Infallible model snapshot, owned configuration, and separate site order; selected ids/order checked during execution |
| StereoEnumerator::transform_iter | Issues EnumerateStereoIter with an owned source snapshot and prepared traversal; rejection yields an empty iterator |
| StereoEnumerator's Transformer implementation | Consuming and in-place calls publish the first checked assignment; trait iteration delegates to the inherent owned iterator |
| EnumerateStereoIter | Owned operation-issued iterator, produced by StereoEnumerator::transform_iter |
| EnumerateStereoIter's Iterator and FusedIterator implementations | Independent Molecule items in supplied site order with ascending stored-frame values; exhaustion or rejection is permanent; no public constructor |
| AddStereoError, RemoveStereoError, EnumerateStereoError | Operation-specific rejection boundaries |
| AddStereoError::Perception, AddStereoError::Inconsistency | Existing shared query/evidence causes, with From conversions; no parallel diagnostic vocabulary |
| RemoveStereoError::Perception | Wraps StereoPerceptionError with a From conversion; no separate removal-specific diagnostic payload |
| EnumerateStereoError::Perception, EnumerateStereoError::Constraint | Wrap shared StereoPerceptionError and ConstraintInvariantsError causes with From conversions |
| EnumerateStereoError::InvalidSiteOrder | Reports missing selected sites and duplicate selected entries, ignoring unselected order entries |
| EnumerateStereoError::NoAssignment | No selected absolute assignment satisfies the checked constraints, including an empty selected domain |

StereoAdder and StereoRemover have new(model) and with_config(model, config)
constructors, following Aromatizer. StereoEnumerator additionally receives the
site order, following Kekulizer's explicit-order construction:

    impl StereoEnumerator {
        pub fn new(model: &StereoModel, site_order: Vec<StereoSite>) -> Self;
        pub fn with_config(
            model: &StereoModel,
            config: EnumerateStereoConfig,
            site_order: Vec<StereoSite>,
        ) -> Self;
    }

Constructors store configuration; contextual id/order checks belong to execution
against the molecule. Errors distinguish unavailable sites or invalid selected
orders, unsupported selected states, required but undetermined input,
incompatible evidence, and absence of a consistent assignment where applicable.
Reuse StereoInconsistency and existing graph-IR diagnostic vocabulary where they
describe the failure; do not create parallel errors for the same condition.

Transform supplies add_stereo / add_stereo_into and remove_stereo /
remove_stereo_into on Molecule. The engines remain the common consuming,
in-place, and iterator interface; enumeration is a separate transformer, not a
free-standing enumeration function attached to StereoAdder.

Python follows the established named-Molecule-method ownership contract for
addition and removal: the base method preserves its receiver and returns an
independent result; the _into method mutates and returns None. Bind shared
selectors and operation configurations faithfully. Python exposes enumeration as
Molecule.enumerate_stereo(model, site_order, *, config=None), returning the native
owned iterator without exposing a transformer class. It preserves native lazy
results, ordering, multiplicity, and exhaustion on rejection. next selects its
first consistent result. The native single-result and in-place methods remain
Rust-only; Python exposes enumeration through independent results.

### Native enumeration iterator

EnumerateStereoIter is an operation-issued concrete iterator over Molecule. It
owns one source snapshot through Molecule's shared storage, the prepared finite
domains and site order, and its traversal state. It has no public constructor
from independent parts. The source snapshot is captured when the iterator is
issued; candidate molecules are constructed lazily on next. Later changes to
the source or engine do not change that issued run.

StereoEnumerator has an inherent method with the same operation name:

    pub fn transform_iter(&self, molecule: &Molecule) -> EnumerateStereoIter;

Its Transformer::transform_iter implementation delegates to that method. This
preserves the common trait and its molecule-only/error-omitting iterator
contract while providing a concrete owned result for native and Python callers.
The iterator can outlive the engine and source. Singular transform and
transform_into use the same preparation and traversal, reporting their errors
and choosing the first surviving assignment. Their successful publication
reuses the consuming or borrowed receiver through the existing planning rules.

Python stores the native iterator directly, following the existing
ReactionApplicationIter and ReactionProductsIter bindings. Issuing the run
captures the source snapshot once; next does not reconstruct previous prefixes
or collect all results. The Rust operation owns this run state, as the existing
reaction iterators do.

No new molecule notation or stored stereo form is introduced. Existing kinds,
cosets, constraints, and frame transport remain the representation. Update
living operation and naming documentation to cover creation/removal of entities
and determination of their configuration fields; the current description of
all transformations as rewrites of determined representations is too narrow for
these operations.

## Algorithm

The source inspection used RDKit commit
78cb05c3db00d32c31c911ef06f4c8a67bc06b98, principally
materials/codes/rdkit/Code/GraphMol/FindStereo.cpp and ConjugHybrid.cpp, and CDK commit
34806d7a54ee4493d57b5c5d1b26a3d748da1cb1, principally
materials/codes/cdk/base/standard/src/main/java/org/openscience/cdk/stereo/Stereocenters.java.
The agreed algorithm below adapts those source cases to umol's representation.
Implementation progress and measurements are recorded in the plan below.

1. **Screen local environments.** Inspect topology and inherent element, charge,
   bond orders, hydrogen counts, and lone-pair counts independently of whether
   #T/#C is present. Use RDKit's isAtomPotentialTetrahedralCenter and CDK's
   createElements/supportedType as the source cases for the atom screen. Cover
   ordinary four-coordinate sites, three neighbors plus hydrogen, P/As with a
   lone-pair position, and the S/Se valence/charge cases such as sulfoxides.
   Treat saturated three-coordinate N, its conjugation, and ring/bridgehead
   cases explicitly; RDKit's ConjugHybrid.cpp supplies the conjugation screen
   used by its N rule. Translate the local valence cases to umol's explicit
   charges, bonds, and lone pairs, rather than importing a boundary format's
   valence shorthand. These automatic-discovery rules do not change what an
   asserted arrangement may represent or require removing existing records.

2. **Build Th frames.** Collect actual neighbors in the established
   order and add the permitted virtual hydrogen/lone-pair positions from the
   inherent counts. Require exactly four positions and the existing distinct
   ligand contract. Do not import a toolkit's duplicate-ligand or multiple-total-H
   rejection into this frame builder: explicit equivalent ligands remain useful
   for symmetric and prochiral entities. S13c supplies the shared
   frame mechanics; the local chemistry screen is additional discovery work.

3. **Build Ct candidates from double bonds.** Require a localized bond of literal
   order two. For each endpoint, exclude the opposite endpoint, collect every
   substituent, and complete exactly two positions using the permitted virtual
   ligand. Preserve the two endpoint blocks and the established frame order.
   Apply StereoModel's minimum-ring-size rule to every candidate, including
   unmarked bonds. The existing assertion-only skipped-bond scan is not enough.
   Include imine/oxime cases with an appropriate lone-pair position.

4. **Classify with molecular symmetry.** Both toolkits compare whole-molecule
   symmetry classes after local eligibility, rather than just neighbor elements.
   Reuse umol's graph-symmetry and local stabilizer machinery for the same job.
   On the binary Th/Ct space, test whether the projected ligand permutations
   exchange the two cosets; use the oriented group for topicity/prochirality.
   A missing literal coset does not supply a negative result. Retain eligible
   candidates of every classification until StereoSelection is applied. Keep
   the current para-stereo behavior without adding a convergence project.

5. **Apply operation policy.** Intersect the site, kind, and stereogenicity
   selectors. Addition retains existing entities and realizes compatible source
   evidence or leaves the new coset open. Removal applies All/Undetermined to
   that selected set and follows the settled assertion-removal rule. Resolver
   recovery and evidence-realization policies remain separate from discovery.

6. **Enumerate binary absolute assignments.** Prepare the allowed domain of each
   selected existing entity: a fixed coset, a literal subset of {0,1}, or both
   values for an undetermined coset. Traverse their product in the supplied site
   order with ascending values. Check stereo assertions affected by each
   complete assignment; keep unrelated open entities unchanged. Do not require
   the whole molecule to become determined. Emit every surviving assignment,
   preserve symmetric/prochiral entities, and perform no canonical-SMILES
   deduplication or geometry embedding. A cursor/domain stack gives bounded
   traversal state and permits stopping after the first result.

## Verification contract

Use exact fixtures, independent finite-domain references, properties, and
benchmarks from the start. In particular:

- A known-kind open coset does not force a negative stereogenicity result.
- Site, kind, and stereogenicity selectors compose by intersection. The selected
  removal mode applies uniformly across the selected kinds.
- In-place rejection preserves the complete input under structural equality;
  consuming and in-place successes agree under structural equality.
- Addition preserves existing configurations and is idempotent.
- Removal clears the selected entities and corresponding source assertions,
  while preserving unrelated data under the actual compaction.
- Enumeration is sound and complete against a bounded exhaustive reference,
  preserves fixed values, retains symmetry-equivalent assignments, and agrees
  with singular selection at its first result.
- Iterator construction does not materialize the product of all domains, early
  stopping does not enumerate the remainder, and yielded molecules are
  independently mutable.
- Relabeling and allowed frame transport preserve classification and the full
  set of assignments under the appropriate transported equivalence. Do not
  assert invariance of the frame-dependent first-result choice.

Benchmark setup, parsing, and independent expected-result construction stay
outside timed sections. Measure local eligibility/frame construction, molecular
symmetry, per-site stabilizers/classification, and operation planning and output
materialization separately. Include consuming/in-place execution with both
shared and uniquely owned storage. For enumeration, separate iterator setup,
first-result latency, bounded prefixes, and full drain; report emitted results
so output volume can be distinguished from per-assignment overhead.

Record molecule size, candidate and selected-site counts, ring structure, and
symmetry characteristics in the input manifest. Include restricted-selection
and unrestricted cases. Record symmetry/stabilizer call counts and enumeration
assignments attempted, rejected, and emitted where needed to explain costs.
Use benchmark-local instrumentation and profiles to establish bottlenecks;
choose optimizations from the measurements as each stage becomes available.

## Implementation plan

Every subitem below is additive (green), including behavior corrections whose
signatures remain unchanged. Each stage ends green. If a public signature
change proves necessary, identify it explicitly as breaking (red-to-green) and
include all of its callers in the same stage; do not let an internal route make
unrelated APIs fallible.

New graph implementation modules are grouped under
umol-graph/src/ops/transform/stereo.rs, with per-engine tests under
ops/transform/stereo/tests/. Use module.rs plus child files, not mod.rs. Shared
perception remains in ops/stereo.rs. Test additions follow test-writing and the
property-test guide; every engine's implementation item includes its exact
success, rejection, preservation, and iterator cases.

### S0 — Correctness references and baseline measurements

- **S0a — Stereo fixtures and finite reference domains** — additive (green).
  Modules: the stereo unit/property support in umol-graph-ir, and new stereo
  fixture/reference support beside umol-graph's transformer tests. [dep: none]
  Establish independently specified cases for tetrahedral and cis/trans sites,
  symmetric and prochiral arrangements, explicit hydrogen positions, remote
  stereo dependencies, small rings, known kinds with open cosets, and retained
  non-default kinds for preservation/removal. Include literal-set, variable,
  invalid-selection, and
  partial-input cases. Build a bounded exhaustive absolute-assignment oracle
  from the stated finite domains and constraints; it must not call the future
  production enumerator or construct expected outputs through the transformers.
  Preserve external comparison provenance without a runtime toolkit dependency.
  Keep this subitem green: install assertions against the corrected production
  classification with S1a, rather than introducing an ignored or failing test.

  **Completed — 2026-10-29.** The
  [graph corpus](../umol-graph/src/ops/transform/tests/stereo/fixtures.rs)
  contains 13 synthetic inputs covering the arrangements above, including
  source stereo assertions, retained square-planar stereo, finite/symbolic
  configurations, and unrelated partial state. Tests check exact stored sites,
  ligand frames, configurations, and integrity-preserving publication.
  The [assignment reference](../umol-graph/src/ops/transform/tests/stereo/reference.rs)
  exhausts binary assignments through eight selected sites, independently
  filters domains and supplied admissibility, and sorts by the requested order.
  The [IR reference](../umol-graph-ir/tests/property/stereo/reference.rs)
  checks all 96 Th/Ct coset/permutation combinations using inversion parity and
  Ct endpoint blocks. Open-coset classification coverage is supplied by S1a.
  Validation passed: 265 transformer unit tests, 29 feature-gated IR stereo
  tests, and Clippy for both crates' libraries/tests with proptest enabled and
  warnings denied.

- **S0b — Existing perception and symmetry baselines** — additive (green).
  Modules: umol-graph/benches/transform.rs and its existing benchmark setup;
  umol-graph-ir benchmark registration only if a separate IR measurement is
  necessary. [dep: S0a]
  Measure existing shared perception and graph/local symmetry over the fixture
  domains before algorithm changes. Keep construction and reference checking
  outside timing. Record the open-coset defect as baseline behavior without
  making it a correctness expectation. Include input-size and stereo-site-count
  series, both shared and unique storage where mutation is measured, and an
  unchanged-input control. Cover differing ring structures and molecular
  symmetries, recording the input characteristics above with each result. Do
  not time unavailable transformer stubs.

  **Completed — 2026-10-29.** The existing transform benchmark registers
  [stereo baselines](../umol-graph/benches/transform/stereo.rs) for 28 inputs and
  five separately timed operations: unchanged-input control, perception,
  whole-molecule symmetry, all local projections, and classification over
  precomputed local symmetries. Setup validates input integrity, retained input
  equality, positive source evidence, known-literal classifications, and ring
  sizes. Observed open-coset results are recorded without asserting them.
  All 140 benchmark cases pass test mode; benchmark Clippy passes with warnings
  denied. The 47 graph reference/fixture tests and two exhaustive IR reference
  tests also pass. These operations borrow their inputs; copying and storage
  detachment are measured with the engine stages.

S0 gate: the new reference support passes, fixtures construct through public
integrity-valid APIs, and the baseline benchmark target builds and runs in test
mode. Save the baseline command, input manifest, and measurements with the
implementation evidence. **Passed — 2026-10-29.**

#### S0 baseline measurements

Production revision: 45f007dbd1b43636644cc38fc62a9a70f05a29ab. Host: Apple M1 Pro
(10 cores, 32 GB), macOS 15.7.3 arm64. Compiler: rustc 1.99.1 (58cf8f9c3).
The bench profile enables optimization and debug information. Configuration:
Nauty, ConstitutionColoring::full(), para_stereo=false, max_iterations=16,
and the default stereo model's minimum ring size of eight. Criterion used
30 samples, 100 ms warm-up, 300 ms measurement, and 10,000 bootstrap resamples.

    cargo bench --locked -p umol-graph --bench transform -- stereo_baseline --save-baseline stereo-s0b-45f007db --sample-size 30 --warm-up-time 0.1 --measurement-time 0.3 --nresamples 10000 --noplot

The complete input domains are:

| Domain | Cases | Construction |
| --- | ---: | --- |
| S0a fixtures | 13 | The stored arrangements, explicit/virtual ligands, source assertions, partial state, and retained square-planar record in the shared corpus |
| Literal controls | 2 | The tetrahedral and cis/trans fixtures with coset zero |
| Atom-count series | 4 | One Th* site with F, Cl, Br, and an alkyl chain; 5, 12, 36, or 132 atoms |
| Site-count series | 5 | Fixed 80 atoms/64 bonds in 16 identical disconnected C(F)(Cl)(Br)I components; 1, 2, 4, 8, or 16 stored Th* sites |
| Ring series | 4 | One simple 4-, 6-, 8-, or 12-membered ring with one double bond, a positive #C assertion, and a stored Ct* entity |

Selected medians below are microseconds per iteration. Local projection and
classification process all stored sites. Parsing, construction, validation,
and manifest generation are outside timing; each phase is measured separately.
The unchanged-input control's median ranges from 0.522 to 0.533 ns.

| Input | Perception | Graph symmetry | Local projection | Classification |
| --- | ---: | ---: | ---: | ---: |
| Five-atom Th* | 0.223 | 2.959 | 1.151 | 0.358 |
| Five-atom Th0 | 0.226 | 3.092 | 1.160 | 0.255 |
| Six-atom Ct0 | 0.233 | 3.363 | 1.228 | 0.180 |
| Explicit methane, Th0 | 0.231 | 10.732 | 11.108 | 4.928 |
| 132 atoms, one Th* site | 1.932 | 74.698 | 26.888 | 0.358 |
| 80 atoms, one Th* site | 1.125 | 167.453 | 75.282 | 0.359 |
| 80 atoms, four Th* sites | 1.547 | 171.018 | 276.266 | 1.446 |
| 80 atoms, sixteen Th* sites | 3.537 | 210.334 | 1390.284 | 5.713 |
| Four-membered asserted ring | 3.747 | 3.996 | 1.892 | 0.343 |
| Twelve-membered asserted ring | 16.645 | 7.506 | 2.945 | 0.342 |

Local projection dominates the 80-atom/16-site case at 6.6 times the whole-molecule
symmetry cost. Its input already includes the precomputed GraphSymmetry; the
current projection computes a carrier stabilizer for every stored site.
Whole-molecule symmetry dominates the 132-atom acyclic case. Ring calculation
within perception is material on the asserted-ring series: perception is the
largest measured phase for the twelve-membered ring. These are measured
synthetic domains; the transformer and discovery costs are measured in S2/S3.

At the S0 baseline, the four-distinct-ligand Th* and Ct* fixtures classified as
Symmetric, while their literal controls classified as Stereogenic. S1a corrects
this undetermined-coset defect. The hydrogen-plus-lone-pair fixture produced
StereoAtomFailure at S0; S1b constructs its frame without an
inconsistency. Its perception timing therefore includes successful construction
instead of failure.

### S1 — Shared stereo foundations

- **S1a — Known-kind open-coset stereogenicity** — additive (green).
  Modules: umol-graph-ir/src/ir/symmetry.rs and the corresponding view and stereo
  property tests; affected graph validator/query callers. [dep: S0a, S0b]
  Correct is_stereogenic and stereogenicity so an undetermined coset does not
  manufacture a negative classification. Use the kind's finite coset action and
  ligand symmetry rather than inserting a chosen literal. Preserve the distinct
  unknown-kind boundary and current known-literal behavior where it is correct.
  Update the existing test that expects false solely for an open coset. Compare
  with independent finite action calculations on Th/Ct and retain
  prochiral/symmetric and malformed-attribute domains. Leave non-binary
  classification behavior outside this correction. Audit
  dependent classification and validation before closing the subitem.

  **Completed — 2026-10-29.** StereoSymmetry classifies known Th/Ct kinds from
  their two-coset action for undetermined cosets and valid nonempty finite coset
  sets. It preserves the supplied coset rather than choosing a literal. Empty
  or invalid finite domains, invalid literal indices, and incompatible actions
  retain the false result; non-binary open configurations retain their prior
  behavior.

  Binary local projection grades the queried Th carrier by permutation parity;
  the queried Ct carrier contributes no handedness. Other carriers retain their
  stored-state grading. This makes prochiral Th arrangements recognizable with
  open cosets and retains the actions identifying Ct cosets at non-stereogenic
  sites, including sites carrying literal configurations. The existing private
  grading helper receives the carrier iterator and accumulated orientation.
  Existing public Molecule and StereoSymmetry method signatures are preserved.

  View delegates use the corrected classification. StereoConformanceValidator
  retains its configuration-determinacy boundary: a known kind with an open
  coset still contributes Underdetermined. Tests cover that case separately
  from an undetermined kind, and cover a non-stereogenic literal Ct record.
  The consumer audit found no direct Python symmetry-query bindings to migrate.

  Verification passed: 11,736 graph-IR unit tests (three ignored), 2,683 graph
  unit tests, 45 feature-gated stereo tests, 683 resolution-conformance tests,
  all-target/all-feature Clippy for both crates, and all 140 benchmark cases in
  test mode. The independent reference exhausts 4,096 one-carrier inputs across
  both kinds, all four-label assignments, literal/open/finite cosets, and both
  refinement settings. It checks ligand-equality predictions, the independent
  binary action, and preservation of the supplied coset.

  The S1a phase comparison repeated 30 affected phase/input cases with S0's
  sampling settings. Longer confirmation runs used 50 samples, 300 ms warm-up,
  and one-second measurement for the large-chain and many-site cases.

  Selected medians in microseconds compare the same phase and input:

  | Phase and input | S0b | S1a |
  | --- | ---: | ---: |
  | Classification, five-atom Th* | 0.358 | 0.260 |
  | Classification, six-atom Ct* | 0.360 | 0.182 |
  | Classification, 80 atoms / 16 Th* sites | 5.713 | 4.147 |
  | Graph symmetry, 132 atoms / one Th* site (confirmation) | 74.698 | 74.514 |
  | Local projection, 132 atoms / one Th* site (confirmation) | 26.888 | 27.189 |
  | Graph symmetry, 80 atoms / 16 Th* sites (confirmation) | 210.334 | 210.275 |
  | Local projection, 80 atoms / 16 Th* sites (confirmation) | 1390.284 | 1402.543 |

  The open-coset rows now return Stereogenic on the distinct-ligand fixtures.
  Their classification costs therefore reflect corrected work. The longer
  confirmations place local-projection cost about one percent above S0 in the
  two larger domains; local projection remains the largest phase on the
  many-site input.

- **S1b — Stereo ligand construction** — additive (green).
  Modules: umol-graph/src/ops/stereo.rs and its resolver/conformance tests.
  [dep: S0a; external: 240/S13c]
  Reuse the existing S13c work item for stereo ligand construction.
  If it is still outstanding, implement and verify that item in its
  owning record before proceeding; do not maintain a second implementation or
  mark 240's independent validator work complete. Read zero/one hydrogen and
  lone-pair counts exactly, preserve frame order, and reject inputs whose
  required ligands cannot form a Th or Ct frame through the established
  perception outcomes.

  **Completed — 2026-10-29.** Implemented and verified under
  [240, S13c](240-validation-constraint-ops-2026-10-18.md#s13--invariant-validation-and-stereo-ligand-construction).
  Perception now uses all actual neighbors and exact zero/one
  hydrogen and lone-pair counts, including both virtual kinds together.
  Wrong arity, a repeated ligand across Ct endpoints, and unsupported or
  non-literal counts use the existing absence and inconsistency outcomes.
  Public interfaces and resolver failure policies are preserved. No standalone
  validator is part of this completion.

S1 gate: focused IR symmetry/stereo properties and graph perception, stereo
resolve, and stereo-conformance tests pass. Re-run the relevant baseline cases
once after these changes; investigate a measured change before optimizing it.

**Passed — 2026-10-29.** The graph suite passed 2,779 unit tests, 186 property
tests and 683 resolution-conformance tests. The feature-gated graph-IR property
target passed 99 stereo-filtered tests. Clippy passed for both crates with all
targets and features and warnings denied. All 140 S0 benchmark cases passed
in test mode. The stereo-discharge property retains its assertions and uses
explicit zero lone-pair counts for its realizable input frames.

All 28 perception cases were measured with S0's host, build, model, and sampling
settings. Valid Ct endpoint frames use inline storage. Selected medians below
are microseconds; the brackets give the S1b median's 95% bootstrap interval.

| Perception input | S0b | S1b [95% interval] |
| --- | ---: | ---: |
| Five-atom Th* | 0.223 | 0.224 [0.220, 0.228] |
| Th* with implicit H | 0.230 | 0.237 [0.235, 0.237] |
| Th* with lone pair | 0.229 | 0.235 [0.235, 0.237] |
| Th* with H and lone pair | 0.150 | 0.226 [0.225, 0.229] |
| Six-atom Ct* | 0.228 | 0.240 [0.239, 0.241] |
| Imine Ct* | 0.204 | 0.216 [0.215, 0.217] |
| 80 atoms / 16 Th* sites | 3.537 | 3.525 [3.511, 3.534] |

The small Ct and imine inputs each take about 12 ns more than S0. The
hydrogen-plus-lone-pair input now constructs a frame and has no
inconsistency, rather than yielding StereoAtomFailure; its timing measures a
different outcome. The sixteen-site median differs by less than one percent.

### S2 — Selectors, configurations, and shared identification

- **S2a — StereoSite and StereoSelection** — additive (green).
  Modules: umol-graph/src/ops/stereo.rs and its tests. [dep: S0a]
  Add the two shared public types and their documented open-carrier contract.
  Implement the common selection behavior with None/empty-set and intersection
  semantics. Site validation is contextual at operation execution. Test atom
  and bond sites, unrestricted and empty selections for all three fields,
  included/excluded kinds, all three stereogenicity classes, and unavailable ids
  without indexing panics. Test intersections of site, kind, and stereogenicity
  restrictions, including excluded kinds at otherwise eligible sites. Keep
  selection separate from the caller's processing order.

  **Completed — 2026-10-30.** Added StereoSite's Atom/Bond variants,
  StereoSelection's three public optional sets, its unrestricted Default, and
  the matches predicate under ops::stereo. Construction and matching accept
  unbound ids without indexing a molecule; executing operations check ids in
  their input context. Missing kind or stereogenicity values satisfy only
  unrestricted selectors for those fields.

  The 26 selector cases cover unrestricted and empty sets for every field,
  atom/bond identity, kind inclusion/exclusion, all three stereogenicity values,
  intersections, missing values, and unbound ids. Exact expected selections
  retain caller-supplied order. Validation passed: all 2,805 graph unit tests
  and all-target/all-feature Clippy with warnings denied. Strict rustdoc,
  formatting, and diff checks passed.

- **S2b — AddStereoConfig** — additive (green).
  Modules: umol-graph/src/ops/transform/stereo.rs and config tests. [dep: S2a]
  Implement AddStereoConfig with the shared StereoSelection and required
  operation-level algorithm settings. Constructors do not check ids against an
  absent molecule. Default the kind selector to Th/Ct, stereogenicity to
  Stereogenic, and sites to unrestricted. Verify kind include/exclude filters,
  their intersection with the other selectors, and explicit opt-in to
  Prochiral/Symmetric.

  **Completed — 2026-10-30.** Added AddStereoConfig under
  ops::transform::stereo and re-exported it from ops::transform. Its public
  surface is selection, automorphism_algorithm, ring_config, and Default.
  Default selects all sites, Tetrahedral/CisTrans kinds, and Stereogenic only,
  with Nauty and RingConfig::default(). Struct construction preserves supplied
  selectors without checking ids against a molecule.

  The 12 configuration cases check the exact default value, atom/bond site
  restrictions, kind inclusion/exclusion, selector intersections, explicit
  Prochiral/Symmetric selection, unrestricted stereogenicity, empty kinds,
  and unbound ids. All 2,817 graph unit tests and all-target/all-feature Clippy
  with warnings denied passed. Strict rustdoc, formatting, and diff checks
  passed.

- **S2c — RemoveStereoMode and RemoveStereoConfig** — additive (green).
  Modules: the same stereo transformation module and config tests. [dep: S2a]
  Add All/Undetermined and RemoveStereoConfig with the same shared
  site/kind/stereogenicity selector. Default to All and unrestricted selectors.
  Test each mode with Th/Ct selections and verify that one mode applies uniformly
  to all selected kinds. Cover excluded kinds at eligible sites, the distinction
  between an open coset and an unknown kind, and the unrestricted All route for
  kindless records. Configuration construction must not force an unnecessary
  symmetry calculation.

  **Completed — 2026-10-30.** Added RemoveStereoMode's All/Undetermined
  variants, its matches predicate on StereoConfigurationForm, and
  RemoveStereoConfig with public selection, mode, and automorphism_algorithm
  fields and Default. Both types are re-exported from ops::transform.
  Default uses All, unrestricted selectors, and Nauty. Construction and mode
  matching perform no molecule lookup, classification, or coset resolution.

  The 63 mode/configuration cases cover Th/Ct and SquarePlanar configuration values,
  literal and undetermined cosets, finite sets, symbolic terms, and unknown
  kinds. They check uniform mode application, selector intersections, kind
  exclusion at eligible sites, empty selections, and the unrestricted All
  route with missing kind/classification values. All 2,880 graph unit tests and
  all-target/all-feature Clippy with warnings denied passed. Strict rustdoc,
  formatting, and diff checks passed.

- **S2d — EnumerateStereoConfig** — additive (green).
  Modules: the same stereo transformation module and config tests. [dep: S2a]
  Add EnumerateStereoConfig with Default selecting Th/Ct, all sites, and all
  stereogenicity values, using Nauty. Reuse the shared
  site/kind/stereogenicity selector. No deduplication or relative-variable
  options are introduced. Test exact defaults, the shared selector's field
  overrides and intersections, empty selections, and construction with unbound
  ids. Keep processing order as the separate constructor argument above,
  following Kekulizer. Traversal-order tests belong to S3c; constructors and
  contextual order validation belong to the executing engine in S3d, where the
  selected processing sites are available. The native iterator producer and
  ownership contract above remain fixed.

  **Completed — 2026-10-30.** Added EnumerateStereoConfig's public selection
  and automorphism_algorithm fields and Default under ops::transform::stereo,
  re-exported from ops::transform. Construction is infallible and does not bind
  ids to a molecule. Its ten configuration cases cover the exact defaults,
  atom/bond site restrictions, kind selection, all stereogenicity values,
  selector intersections, empty kinds, and unbound ids. All 2,890 graph unit
  tests and all-target/all-feature Clippy with warnings denied passed. Strict
  rustdoc, formatting, and diff checks passed.

- **S2e — Shared identification and classification** — additive (green).
  Modules: umol-graph/src/ops/stereo.rs, with focused resolver and validator
  callers and tests. [dep: S1a, S1b, S2a]
  Extend StereoPerception to identify unasserted candidates and classify selected
  existing/candidate entities using the shared symmetry machinery. Reuse
  existing kind/coset/frame values and StereoInconsistency. Add public recognize
  and classify queries with StereoPerceptionError as specified above; local
  eligibility checks remain private. Compute context from the complete molecule, apply candidate
  restrictions without letting unrelated unsupported fields reject the chosen
  scope, and honor the model's element and ring rules. Preserve derive's
  evidence-driven resolver behavior and all source assertion/frame transport.
  Implement the local Th/Ct discovery screens described above separately from
  assertion realization. In particular, test unmarked double bonds, endpoint
  ligand lists against their neighbors and counts, and ring eligibility without
  relying on the assertion-only skipped-bond scan. Retain symmetric/prochiral
  candidates until selector application. Keep current para-stereo behavior and
  limits; improving them is outside this stage.
  Test explicit negatives, positive evidence, competing candidate kinds,
  already-present entities, kind coverage, and genuinely unsupported
  producers. Document actual supported producer coverage without silently
  ignoring enabled configuration. Add phase measurements for local candidate
  discovery, molecular symmetry, and per-site classification to the established
  benchmark target, including restricted-selection cases.

  **Completed — 2026-10-30.** Added public recognize and classify methods and
  StereoPerceptionError under ops::stereo. Automatic Th recognition covers
  four-coordinate arrangements, one virtual hydrogen, P/As lone-pair sites,
  S/Se valence/charge cases, and saturated pyramidal N in a three-membered ring
  or at a bridged-ring endpoint, subject to its local conjugation screen.
  The N screen checks adjacent pi bonds, main-group electron contributions,
  neutral hypervalence, and the heavier P/S-group exclusions used by the
  inspected reference algorithm. Existing asserted arrangements use the
  structural frame builders. Ct recognition checks unmarked double bonds,
  endpoint frames including imine/oxime lone pairs, and the ring threshold.
  Ring data is reused within a query, with a search bound at least the
  configured minimum ring size minus one. Symmetric and prochiral candidates
  remain available to selectors. Neither query changes resolver policy or
  adds a Molecule method.

  The 118 focused query cases cover recognition, literal/open/finite/symbolic
  evidence, negatives, existing entities and competing kinds, frame failures,
  selectors, ring thresholds including a 24-membered ring, model scope and
  disabled producers, required partial input, classification independent of
  automatic eligibility, and exact error payloads. Two generated properties
  compare Th/Ct recognition and classification with independent label-equality
  references over explicit H/F/Cl/Br/I ligands and binary coset domains.
  The full 3,008 graph unit tests, 188 graph property tests, 45 graph-IR stereo
  property tests, and 683 resolution conformance cases passed. All-target,
  all-feature Clippy for graph/graph-IR, strict graph rustdoc, formatting, and
  diff checks passed. The established benchmark target validates all phase
  inputs outside timing; measurements are recorded below.

  Property and resolution test roots now live at tests/property/suite.rs and
  tests/resolution/suite.rs, with explicit Cargo target paths and ordinary
  child module declarations. The conformance body moved unchanged to
  resolution_suite.rs; test names, snapshots, and saved property cases retain
  their existing locations. suite.rs avoids proptest treating a new main.rs as
  a different saved-case root. Rust/Python query bindings and their result and
  diagnostic values are included in S4c/S4d.

S2 gate: all new public configuration fields, constructors, defaults, visibility,
and failure boundaries are listed and reconciled with this record. Selector,
config, perception, resolver, and validator tests pass. No new DSL representation
is required; do not introduce a stored stereogenicity label as a shortcut.
Record the discovery/classification measurements before proceeding to engines.
**Passed — 2026-10-30.** Public query/configuration surfaces and their failure
boundaries are listed above; existing resolver and validator callers continue
using derive without signature or behavior changes.

#### S2 recognition and classification measurements

Host, compiler, optimization profile, Nauty, full constitution coloring,
para_stereo=false, and minimum ring size eight match S0b. The new
stereo_recognition group uses 25 named inputs from the S0 manifest with source
assertions and stored entities removed; the retained square-planar, symbolic
partial-state, and acyclic N hydrogen/lone-pair inputs are excluded from this
automatic-recognition domain. The 80-atom series explicitly selects 1, 2, 4, 8,
or 16 sites in the same 16-component topology. The ring series is unmarked.
Independent expected frames and configurations are checked before timing.

Each input has five phases: local candidates without a class filter; whole
molecular symmetry after prospective additions; local projections and
classification using that precomputed symmetry; the public classify query;
and recognize filtered to Stereogenic, including prospective construction.
The latter two query phases compute one molecular symmetry and one local
projection per candidate, or neither when there are no candidates.

    cargo bench --locked -p umol-graph --bench transform -- stereo_recognition --save-baseline stereo-s2e-final --sample-size 30 --warm-up-time 0.1 --measurement-time 0.3 --nresamples 10000 --noplot

All 125 phases were measured. The initially noisy 132-atom and 80-atom/one-site
and four-site controls were repeated with 50 samples, 300 ms warm-up, and one
second measurement; the table uses that isolated repeat for those inputs.

    cargo bench --locked -p umol-graph --bench transform -- 'stereo_recognition/.*/(sites_1|sites_4|atoms_132)$' --save-baseline stereo-s2e-controls --sample-size 50 --warm-up-time 0.3 --measurement-time 1 --nresamples 10000 --noplot

Medians below are microseconds per iteration. Parsing, expected-result
construction, and setup for the isolated symmetry/projection phases are
outside timing.

| Input | Local candidates | Graph symmetry | Local classification | classify | Filtered recognize |
| --- | ---: | ---: | ---: | ---: | ---: |
| Five-atom Th | 0.279 | 2.965 | 1.471 | 4.609 | 5.518 |
| Explicit methane, Th | 0.281 | 7.002 | 15.289 | 22.687 | 23.580 |
| Six-atom Ct | 0.891 | 3.232 | 1.460 | 4.981 | 6.370 |
| 132 atoms, one candidate | 4.202 | 74.032 | 27.562 | 102.217 | 107.818 |
| 80 atoms, one selected site | 0.498 | 168.360 | 75.854 | 247.772 | 248.750 |
| 80 atoms, four selected sites | 1.123 | 173.341 | 281.282 | 461.845 | 458.735 |
| 80 atoms, sixteen selected sites | 2.908 | 212.882 | 1425.375 | 1631.235 | 1659.985 |
| Twelve-membered ring | 17.072 | 7.626 | 3.060 | 11.023 | 28.949 |

Local recognition is small on the acyclic inputs. Molecular symmetry dominates
the 132-atom/one-candidate case; local projections and classification dominate the
80-atom/16-site case, at 6.7 times its whole-molecule symmetry cost. The
twelve-membered ring's local-candidate phase includes ring perception and
accounts for most of its filtered recognition time. These are separate phase
measurements, not alternative algorithms; they establish costs in the stated
synthetic domains. No optimization or para-stereo redesign is added here.

### S3 — Transformer engines

- **S3a — StereoAdder and AddStereoError** — additive (green).
  Modules: umol-graph/src/ops/transform/stereo.rs, its adder tests, re-exports in
  ops/transform.rs, and umol-graph/benches/transform.rs.
  [dep: S2b, S2e]
  Implement new, with_config, and all three Transformer methods with complete
  planning before mutation. Reuse source coset evidence and the perceived frame;
  add no arbitrary configuration or assertion of stereogenicity. Preserve
  existing entities and reject incompatible competing kind plans. Reuse existing
  diagnostic payloads in AddStereoError where possible. Test full outputs,
  selector intersections, positive/negative evidence, non-stereogenic opt-in,
  idempotence, unrelated partial fields, and exact unchanged input on
  every rejection. Verify the one-result iterator and independent storage.
  Add correctness-checked consuming/in-place benchmark cases with this engine.

  **Completed — 2026-10-30.** Added StereoAdder with new(model) and
  with_config(model, config), its three Transformer methods, and AddStereoError;
  the engine and error are exported from ops::transform and its stereo module.
  Each call uses shared recognize once, rejects its errors or first reported
  inconsistency, and prepares additions through Edits. Consuming publication
  uses an editor; in-place publication applies the same plan directly after
  successful planning. The lazy iterator follows the existing deterministic
  transformer implementation. No recovery copy, journal, or public planning
  surface is added. Model, selector, and algorithm interpretation remain with
  the existing configuration and perception objects.

  The 67 adder cases cover full Th/Ct outputs, virtual H and lone-pair frames,
  positive/negative and finite/symbolic evidence, retained frames and other
  kinds, remote stereo, selector intersections and empty sets, model scope and
  disabled kinds, symmetric/prochiral opt-in, idempotence, partial data, exact
  rejection errors and unchanged inputs, and independent iterator results.
  Cases with a valid early candidate and a later failed site verify that a
  rejection publishes no additions. Existing assertions and stored cosets are
  preserved.

  All 3,075 graph unit tests, 188 graph property tests, and 683 resolution
  conformance cases passed. All-target/all-feature Clippy for graph and
  graph-IR, strict graph rustdoc, formatting, and diff checks passed. All 20
  adder benchmark cases passed their complete-result checks for consuming and
  in-place calls with shared and uniquely owned input storage. The S4a
  transformer properties and S4b receiver methods remain planned work.

- **S3b — StereoRemover and RemoveStereoError** — additive (green).
  Modules: the same implementation, remover tests, re-exports, and benchmark
  target. [dep: S2c, S2e]
  Implement all three Transformer methods. Intersect the operation's mode with the
  shared site/kind/stereogenicity selectors; bypass classification when the
  class selector is unrestricted. Remove corresponding #T/#C assertions on
  removed sites, following Kekulizer, and delegate entity/reference compaction to
  graph IR. Test All and Undetermined separately from stereogenicity, Th/Ct
  include/exclude selectors, known kinds with open cosets, unrestricted kindless
  removal, literal sets, and
  symbolic attributes that can be removed without evaluation. Check retained
  frames/attributes, transported constraint references, assertions at untouched
  sites, no-op results, and exact failure preservation. Add consuming/in-place
  measurements and selected-versus-unrestricted controls.

  **Completed — 2026-10-30.** Added StereoRemover with new(model) and
  with_config(model, config), all three Transformer methods, and the
  RemoveStereoError::Perception wrapper. The engine and error are exported from
  ops::transform and its stereo module. Explicit selection ids are checked
  before empty-filter handling. Mode, site, and kind filters precede the
  optional shared classify call, which uses the original whole molecule.
  No class filter means no geometry, coset-domain, or groundness evaluation.
  Planning records the exact stored frames and attributes for bulk graph-IR
  removal, plus corresponding #T/#C assertion removals. The existing graph-IR
  operations handle id compaction and constraint transport. No Molecule methods
  or additional public planning surface are added.

  The 104 remover cases cover All and Undetermined, shared selector
  intersections and empty sets, Th/Ct and other kinds, kindless and symbolic
  records, literal and finite cosets, invalid and empty domains, partial inherent
  fields, symmetric/prochiral/stereogenic classes, remote stereo context,
  retained frames and attributes, corresponding and untouched assertions,
  transported atom/bond stereo constraint references, no-op results,
  idempotence without a class filter, exact errors and unchanged inputs, and
  independent iterator results. A later classification failure after a valid
  selected site publishes no removals.

  All 3,179 graph unit tests, 188 graph property tests, and 683 resolution
  conformance cases passed. All-target/all-feature Clippy for graph and
  graph-IR, strict graph rustdoc, formatting, and diff checks passed. All 20
  remover benchmark cases passed their full-result checks for consuming and
  in-place calls with shared and uniquely owned storage. The S4a transformer
  properties and S4b receiver methods remain planned work.

- **S3c — Absolute assignment enumeration** — additive (green).
  Modules: internal enumeration support in ops/transform/stereo.rs or a named
  child module, with operation-local reference tests. [dep: S2d, S2e]
  Implement bounded-state ordered traversal of admitted absolute coset domains.
  Preserve fixed values, enumerate literal-set members without broadening their
  domains, check supported constraints, and backtrack on inconsistent choices.
  Check stereo assertions affected by the assignment while preserving unselected
  open entities. Use stereo-specific checks rather than requiring the composite
  validator to report a determined whole molecule.
  Selected unimplemented variable/relative expressions report unsupported input.
  Reuse existing normalization when a constant expression already reduces to
  an admitted absolute value; do not implement a variable solver. Compare exact
  sequences with S0a's bounded exhaustive oracle, including empty and rejected
  domains, included/excluded kinds, and symmetry-equivalent assignments. Test
  different supplied orders over atom and bond sites and preserve their traversal
  order. Verify early stop without retaining the full Cartesian product.

  **Completed — 2026-10-30.** Added the private
  ops/transform/stereo/enumeration.rs module with AbsoluteStereoAssignments and
  StereoAssignmentChecks. The traversal takes the selected configurations in
  processing order, normalizes constant expressions through the existing
  kind-aware implementation, and prepares finite ascending domains. Literal
  cosets stay fixed; an empty absolute domain produces no assignments. Variables,
  including finite or empty variable domains, remain unsupported. Invalid
  indices, unsupported kinds, and incompatible site/kind combinations reuse
  StereoPerceptionError. Its UnsupportedConfiguration display text now applies
  to stereo processing generally, without adding another diagnostic type.

  The cursor stores one position and at most two values per selected site. It
  yields the ordered product incrementally, including the single empty
  assignment for empty work, and remains exhausted after completion. There is
  no symmetry deduplication or Cartesian-product collection. Constraint checks
  use the supplied candidate without changing its configuration or any retained
  data. Source #T/#C assertions use the existing evidence readings, including
  stored-frame transport and retained source restrictions in logical clauses.
  Th/Ct symmetry assertions use the same molecular/local symmetry machinery as
  perception, with one lazy molecular symmetry query and cached local results
  per candidate. Other leaves inside affected logical constraints use the
  existing model-independent evaluator. Every branch is evaluated so an
  unsupported operation cannot be hidden by a preceding true disjunct or false
  conjunct. Open retained fields are allowed; fluxional dynamics and the current
  para-stereo limitations are unchanged.

  The 221 cases cover exact ordered sequences against S0a's independent oracle
  through eight binary sites, mixed atom/bond orders, fixed and finite domains,
  rejected and empty domains, constant expressions and unsupported variables,
  unavailable sites and incompatible frames, frame parity, joint logical
  admissibility, retained source evidence and partial fields, inline and
  molecule-scope symmetry assertions, and retention of both cosets on symmetric
  and prochiral entities. A 128-site case checks the first two results and cursor
  state without materializing the product. The existing oracle's visibility is
  widened only within the transform test modules; its algorithm and assertions
  are unchanged.

  All 3,400 graph unit tests, 188 graph property tests, and 683 resolution
  conformance cases passed. All-target/all-feature Clippy for graph and graph-IR,
  strict graph rustdoc, formatting, and diff checks passed. S3c added no public
  symbols or Molecule methods. S3d connects the module to the public engine and
  removes its temporary non-test dead_code expectation. S4a's transformer
  properties remain future work.

- **S3d — StereoEnumerator and EnumerateStereoError** — additive (green).
  Modules: the stereo engine module, enumerator tests, ops/transform.rs
  re-exports, and the transform benchmark. [dep: S3c]
  Add new(model, site_order) and with_config(model, config, site_order), storing
  configuration and order without contextual validation. Check that new agrees
  with with_config using EnumerateStereoConfig::default(). At execution, project
  the caller's order onto selected processing sites, following Kekulizer;
  reject missing or duplicate selected sites and ignore entries outside that
  set. Test empty work, both site families, and excluded kinds that create no
  order requirements. Validate supplied selection ids before indexing or mutation.
  Expose the configured/order-bearing Transformer. Consuming and in-place
  methods select the first complete admissible assignment; the iterator uses
  that same traversal for all assignments. Add EnumerateStereoIter with owned
  source/configuration/traversal state, produced by the engine's inherent
  transform_iter; the trait implementation delegates to it. Follow the existing
  native reaction iterator ownership pattern. Test that an issued run survives
  the engine/source lifetime and is unaffected by subsequent source mutation.
  Reuse consuming storage for the
  singular result and plan rejection before in-place mutation. Test first-item
  agreement, complete ordered results, unchanged retained fields and frames,
  unsupported selected expressions, no solutions, no selected open fields,
  repeated exhaustion, source preservation, and independent outputs. Measure
  first-result latency, bounded prefixes, and complete enumeration separately;
  cover products with duplicate molecular isomers to confirm no hidden dedup.

  **Completed — 2026-10-30.** Added StereoEnumerator, EnumerateStereoIter,
  and EnumerateStereoError, exported from ops::transform and its stereo module.
  new(model, site_order) delegates to with_config with the default config.
  Constructors capture a model snapshot and own the configuration and order;
  they perform no contextual lookup. Execution validates explicit selection ids
  before empty-filter handling, applies site/kind filters, and requests shared
  classification only when needed. A temporary input normalizes constant
  selected expressions for classification without changing the source. The
  caller's order is projected onto the selected entities; missing and duplicate
  selected entries reject, while other entries are ignored. Existing entities
  are independent of automatic recognition eligibility, enabled discovery kinds,
  and the discovery ring threshold. Required fluxional processing remains
  unsupported.

  Consuming and in-place operations use the same prepared traversal to find the
  first checked candidate and its Edits plan. The probe is discarded before
  publication through the original receiver's editor or direct mutation;
  rejection leaves a borrowed receiver unchanged. Edits change only selected
  configurations. The inherent transform_iter captures an owned source snapshot
  and prepared domains/order/checks. next constructs candidates lazily, retains
  symmetry-equivalent assignments, and terminates permanently on exhaustion or
  rejection. The trait method delegates to it. The iterator has no public
  constructor or borrowed engine/source state. No Molecule methods or Python
  bindings are added in this item.

  The error boundary wraps shared perception and constraint-evaluation causes;
  InvalidSiteOrder supplies missing and duplicate site vectors, and NoAssignment
  covers empty or unsatisfied domains. The 105 cases cover full Th/Ct outputs,
  constructor agreement, selector intersections and empty sets, fixed/finite
  domains, constant expressions, unsupported variables/kinds, model snapshots
  and scope, source assertions in retained frames, joint logical constraints,
  partial and opaque unselected data, exact errors and unchanged inputs,
  both cosets on symmetric/prochiral entities, and native/trait iteration
  agreement. Ordered molecule sequences match S0a's exhaustive oracle through
  eight sites. An issued iterator survives engine/source destruction and source
  mutation; yielded results are independently mutable. The 128-site prefix case
  checks lazy output and first-result agreement without collecting the product.

  All 3,505 graph unit tests, 188 graph property tests, and 683 resolution
  conformance cases passed. All-target/all-feature Clippy for graph and graph-IR,
  strict graph rustdoc, formatting, and diff checks passed. The transform target's
  360 stereo benchmark paths passed in test mode, including all 55 new
  enumerator paths with full expected-result and integrity checks. Enumeration
  measurements and S4a's properties are recorded below. S4b's receiver methods
  and S4c/S4d's Python surface remain planned work.

S3 gate: the complete transformer unit group and the new stereo benchmark cases
pass in test mode. Audit all public exports and error variants, source assertion
handling, and producers/consumers before declaring any engine complete.
**Passed — 2026-10-30.** All three engines, their exports and failure boundaries,
source assertion handling, and the owned iterator's production and consumption
have been checked against the agreed contracts.

#### S3a addition measurements

Measured the S3a working tree on S0b's host and optimization profile with rustc
1.99.1, the default stereo model, and AddStereoConfig::default(): Th/Ct,
Stereogenic, all sites, Nauty, default ring algorithms, para_stereo=false, and
minimum ring size eight. Five inputs exercise a five-atom Th addition, a
six-atom Ct addition, a seven-atom input with positive Th/Ct assertions, one
already-present Th entity with a retained frame, and sixteen missing Th
entities in the same 80-atom/64-bond component series used in S2.

Each input is measured through consuming and in-place calls with shared and
uniquely owned storage. Shared setup clones the fixture; unique setup
constructs a fresh molecule from owned entries. Both setups and independent
expected-result construction are outside timing. Before timing, each of the
four paths is checked under full molecule equality, including source
assertions, and the expected result passes integrity checking. Every timed call
includes recognition, planning, and publication.

    cargo bench --locked -p umol-graph --bench transform -- stereo_adder --save-baseline stereo-s3a-adder --sample-size 30 --warm-up-time 0.1 --measurement-time 0.3 --nresamples 10000 --noplot

All 20 cases were measured. Medians are microseconds per call:

| Input | Consuming, shared | In-place, shared | Consuming, unique | In-place, unique |
| --- | ---: | ---: | ---: | ---: |
| Five-atom Th | 6.112 | 6.062 | 6.047 | 6.114 |
| Six-atom Ct | 6.972 | 6.949 | 7.045 | 7.015 |
| Positive Th/Ct evidence | 11.250 | 11.614 | 11.398 | 11.214 |
| Retained Th, no additions | 0.451 | 0.335 | 0.463 | 0.349 |
| 80 atoms, sixteen additions | 1666.205 | 1685.420 | 1667.052 | 1656.847 |

The sixteen-site total is consistent with S2's phase measurements, where local
projection and classification are the largest measured component. These totals
do not isolate planning and publication from perception; no optimization is
inferred from the small shared/unique differences in this synthetic domain.

#### S3b removal measurements

Measured on S0b's host and optimization profile with rustc 1.99.1, the default
stereo model, Nauty, and para_stereo=false. The nine-atom/seven-bond input has
one fixed Th entity, one explicitly undetermined Ct coset, corresponding #T/#C
assertions, and unrelated #D/#R assertions. All removes both entities;
Undetermined removes only Ct. Three controls share an 80-atom/64-bond input
containing sixteen disconnected C(F)(Cl)(Br)I components, each with Th* and
#T+. All removes sixteen entities; explicit site selection removes one; a
Stereogenic filter removes sixteen after classification. All and Stereogenic
produce the same full result on this input.

Each input is measured through consuming and in-place calls with shared and
uniquely owned storage. Shared setup clones the fixture; unique setup builds a
fresh molecule from owned entries. Setup and independent expected-result
construction are outside timing. Full molecule equality and result integrity
are checked before timing for all four paths. Timed calls include selection,
classification when requested, planning, and publication.

    cargo bench --locked -p umol-graph --bench transform -- stereo_remover --save-baseline stereo-s3b-remover --sample-size 30 --warm-up-time 0.1 --measurement-time 0.3 --nresamples 10000 --noplot

All 20 cases were measured. Medians are microseconds per call:

| Input and selection | Consuming, shared | In-place, shared | Consuming, unique | In-place, unique |
| --- | ---: | ---: | ---: | ---: |
| Th/Ct, All | 1.613 | 1.522 | 1.044 | 1.053 |
| Th/Ct, Undetermined | 0.967 | 0.908 | 0.661 | 0.648 |
| Sixteen Th, All | 5.047 | 4.997 | 3.807 | 3.772 |
| Sixteen Th, one selected site | 3.487 | 3.413 | 1.585 | 1.587 |
| Sixteen Th, Stereogenic | 1631.467 | 1634.271 | 1628.018 | 1624.589 |

On the matched sixteen-site input, requesting classification adds about
1.62 ms to a removal costing 3.8–5.0 µs without it. Restricting removal to one
site reduces the unclassified total to 1.6–3.5 µs. These controls distinguish
the cost of requested classification from raw selection/removal; they do not
separate symmetry projection, edit planning, and publication internally or
establish performance on a broader molecular domain.

#### S3d enumeration measurements

Measured on S0b's host and optimization profile with rustc 1.99.1, the default
stereo model, Nauty, and para_stereo=false. Each selected configuration is an
explicitly undetermined Th or Ct coset. The input manifest is:

| Input | Atoms/bonds | Selected entities | Assignment product | Emitted results |
| --- | ---: | --- | --- | --- |
| mixed | 9/7 | One Th and one Ct | 4 | 4 |
| evidenced | 9/7 | Same sites with #T0/#C1 assertions | 4 | 1, the last attempted assignment |
| symmetric_2 | 10/8 | Two Th entities on explicit C(H)(H)(H)H components | 4 | 4, including symmetry-equivalent assignments |
| all_8 | 40/32 | Eight Th entities on C(F)(Cl)(Br)I components | 256 | 256 |
| classified_8 | 40/32 | Same sites and input, Stereogenic filter | 256 | 256 |
| selected_1_of_128 | 640/512 | One of 128 existing Th entities | 2 | 2; other cosets remain undetermined |
| prefix_128 | 640/512 | All 128 Th entities | 2^128 | First four only; full collection is omitted |

All inputs are acyclic. The mixed inputs retain unrelated #D/#R assertions and
non-default stored frames. Expected configurations are specified independently
of the transformer, outside timing. Before timing, the four singular paths and
native prefix/full output sequences are checked under full molecule equality;
expected results pass integrity checking. Shared singular setup clones the
fixture, while unique setup reconstructs fresh owned entries. Both are excluded
from timing.

Iterator setup measures issuance and disposal of an unconsumed run. Other
iterator timings receive prepared runs with setup excluded: first requests one
result, prefix collects at most four, and full collects every result. The
candidate checks, output construction, and collection work are included.
These separate measurements are not an additive timing decomposition.

    cargo bench --locked -p umol-graph --bench transform -- stereo_enumerator --save-baseline stereo-s3d-enumerator --sample-size 30 --warm-up-time 0.1 --measurement-time 0.3 --nresamples 10000 --noplot

All 55 paths were measured. Medians are microseconds per call:

| Input | Consuming, shared | In-place, shared | Consuming, unique | In-place, unique |
| --- | ---: | ---: | ---: | ---: |
| mixed | 1.556 | 1.444 | 1.239 | 1.278 |
| evidenced | 4.785 | 4.561 | 4.316 | 4.269 |
| symmetric_2 | 1.303 | 1.211 | 1.095 | 1.131 |
| all_8 | 3.268 | 3.012 | 2.743 | 2.866 |
| classified_8 | 225.442 | 230.278 | 226.152 | 227.962 |
| selected_1_of_128 | 11.503 | 10.890 | 9.588 | 9.569 |
| prefix_128 | 47.182 | 43.094 | 40.111 | 42.314 |

| Input | Iterator setup | First after setup | Prefix after setup | Full after setup |
| --- | ---: | ---: | ---: | ---: |
| mixed | 0.460 | 0.514 | 2.018 | 2.046 |
| evidenced | 0.685 | 3.369 | 3.352 | 3.346 |
| symmetric_2 | 0.472 | 0.405 | 1.434 | 1.492 |
| all_8 | 1.390 | 0.856 | 3.075 | 204.864 |
| classified_8 | 226.815 | 1.009 | 3.131 | 202.640 |
| selected_1_of_128 | 1.227 | 1.903 | 3.766 | 3.780 |
| prefix_128 | 23.127 | 10.714 | 39.572 | Omitted |

The matched eight-site controls place the extra classification cost at
preparation: setup increases from 1.390 to 226.815 µs, while collecting all 256
prepared results takes about 203–205 µs in both controls. The first four results
from 128 selected sites cost 39.572 µs after setup; no product-sized preparation
or full result collection is performed. These synthetic inputs distinguish
selection/classification setup from lazy traversal and output volume. They do
not isolate every allocation, constraint check, or publication cost, or establish
performance for a broader molecular domain. No optimization is inferred from
the small differences between consuming/in-place or shared/unique paths.

### S4 — Properties, receiver methods, and Python

- **S4a — Stereo transformer properties** — additive (green).
  Modules: umol-graph/tests/property/transform/stereo.rs, its registration, and
  operation-local strategies/reference support. [dep: S3a, S3b, S3d]
  Check the stated preservation, idempotence, rejection, first-result, soundness,
  completeness, and independent-result properties. Generate known-kind open
  cosets, kind include/exclude selections, both removal modes, explicit
  prochiral/symmetric frames, and unrelated partial state; keep malformed
  attribute domains explicit. Compare transported classification and full result
  collections under atom relabeling and allowed ligand-frame actions, without
  asserting that the first chosen
  representation is frame invariant. Keep the bounded exhaustive comparison
  independent of the production traversal.

  **Completed — 2026-10-30.** Added 27 generative property cases with local
  strategies and independent references. Inputs contain zero to four explicit
  Th/Ct sites with H/F/Cl/Br/I terminal atoms, arbitrary stored frames, site, kind,
  and stereogenicity selections, source #T/#C assertions, and an unrelated Ne atom
  with symbolic charge. Addition and classification use explicitly undetermined
  cosets, with ligand classes supplied by terminal-element equality. Raw removal
  and enumeration also vary fixed cosets and finite sets, including empty sets.
  Unsupported kinds, kindless configurations, invalid indices, symbolic/relative
  cosets, and incompatible site/kind combinations are explicit cases on generated
  frames.

  Independently constructed entries check preservation, removal's record-id and
  constraint compaction, consuming/in-place agreement, failure without
  publication, and result independence. Idempotence is checked for addition and
  removal without a stereogenicity filter. An exhaustive binary-mask reference
  checks enumeration's ordered results, source-domain membership, source-frame
  assertions, joint selected constraints, completeness, multiplicity, and first
  result. It does not call perception, validation, or the production traversal.

  Atom relabeling transports ids, selectors, order, and implicit #T/#C assertion
  frames. Stored-frame transport uses producer-issued representative actions and
  nonempty absolute domains. Classification and represented addition/removal
  results are compared with their stated equality relations; enumeration compares
  full collections with multiplicity, without a first-choice frame-invariance
  assertion. Public rustdoc states these laws and their equality relations.

  Shared StereoPerception::classify now rejects Th on bonds and Ct on atoms with
  UnsupportedConfiguration, matching the supported site/kind boundary already
  enforced by enumeration. Class-filtered removal propagates that shared error;
  unrestricted All removal still strips those records. No public symbols,
  constructors, or Molecule methods were added.

  All 3,505 graph unit tests, 215 graph property cases, and 683 resolution
  conformance cases passed. All-target/all-feature graph and graph-IR Clippy,
  strict graph rustdoc, formatting, and diff checks passed.

- **S4b — Transform receiver interface** — additive (green).
  Modules: umol-graph/src/ops/transform.rs and its tests. [dep: S3a, S3b]
  Add add_stereo / add_stereo_into and remove_stereo / remove_stereo_into as
  direct delegates with the same model, selection, configuration, errors, and
  ownership semantics. Keep the separately configured StereoEnumerator's three
  methods as the enumeration interface; do not replace it with the discarded
  free function or merge it back into addition. Test receiver/engine agreement
  and sequential publication without an extra molecule clone or validation pass.

  **Completed — 2026-10-30.** Added the four agreed methods to Transform and its
  Molecule implementation. Each accepts an explicit StereoModel and operation
  configuration, constructs StereoAdder or StereoRemover, and directly forwards
  the owned or borrowed receiver. The receiver methods reuse the existing errors
  and add no molecule clone or validation pass. Enumeration remains on the
  separately configured StereoEnumerator.

  Eleven new receiver cases compare independently specified full molecules with
  direct transformer results and consuming/in-place receiver results. They cover
  model and selector forwarding, All/Undetermined removal, source assertion
  preservation and removal, exact perception/inconsistency errors, and unchanged
  receivers on rejection. A sequence of successful addition and removal followed
  by rejected addition preserves the preceding result. All 3,516 graph unit
  tests, 215 graph property cases, and 683 resolution conformance cases passed.
  All-target/all-feature graph and graph-IR Clippy, strict graph rustdoc,
  formatting, and diff checks passed.

- **S4c — Python configuration and selector values** — additive (green).
  Modules: umol-py/src/transform.rs, relevant stereo/algorithm wrappers,
  src/lib.rs registration, python/umol/__init__.py, and test_transform.py.
  [dep: S2b, S2c, S2d, S2e]
  Bind StereoSite, StereoSelection, RemoveStereoMode, and the three operation
  configurations. Bind StereoSelection.kinds as well as sites and
  stereogenicities, and bind StereoSelection::matches as a direct Rust delegate
  with optional kind and stereogenicity values. Bind RemoveStereoMode::matches
  using the existing StereoConfigurationForm wrapper. Preserve None/empty-set
  distinctions. Use
  from_rust/to_rust and the established immutable config conventions. Test exact values, selector
  intersections, both removal modes, imports, constructor failures, and attempted
  mutation. Reuse existing enum wrappers. Bind StereoPerception and the
  StereoDerivation and StereoInconsistency values needed by its public queries.
  Map StereoPerceptionError through the existing exception vocabulary, retaining
  its site, kind, field, and shared inconsistency information.

  **Completed — 2026-10-30.** Registered and exported StereoSite, StereoSelection,
  RemoveStereoMode, AddStereoConfig, RemoveStereoConfig, EnumerateStereoConfig,
  StereoPerception, StereoDerivation, and StereoInconsistency. Added the existing
  Rust AtomFieldKind enum to preserve typed atom-field diagnostics. The explicit
  Python export inventory includes all ten classes. Site, inconsistency, and
  atom-field values retain Rust equality, ordering, and hashing.

  Configuration constructors preserve Rust defaults and all shared selector
  fields, including kinds and None/empty-set distinctions. They borrow supplied
  selectors and copy their sets once into the retained configuration. Selection
  matching and removal-mode matching delegate directly to Rust. Immutable
  properties expose frozensets and the existing algorithm/ring wrappers.

  Bound StereoPerception construction, derive, skipped_stereo_bonds, and the two
  frame builders. Queries borrow the molecule; frame builders preserve optional
  results and return independently mutable forms. StereoDerivation preserves its
  four open fields, including independent construction, sequence order, and
  duplicates. Its properties return tuples and read-only form copies. Existing
  exceptions retain the Rust message and typed site, kind, atom/bond, field, or
  inconsistency payloads. The Python API guide records constructor syntax,
  immutability and copying, optional matching arguments, and the frame builders'
  bounds-check adaptation. recognize/classify and the transformation operations
  remain in S4d.

  Twenty-six new Rust cases and 68 new Python cases cover conversions, defaults,
  selector intersections, both removal modes, constructor errors, immutable
  observations, perception results, diagnostic fields, comparisons, and exports.
  On Python 3.13.16, the Rust binding suite passed 1,768 cases with two ignored;
  the freshly rebuilt extension passed 3,078 Python cases with two skipped.
  All-target/all-feature binding Clippy, strict binding rustdoc, Ruff, formatting,
  and diff checks passed.

- **S4d — Python operation bindings and iterator lifecycle** — additive (green).
  Modules: umol-py/src/molecule.rs, src/transform.rs, src/error.rs,
  src/stereo/perception.rs, and the transform/molecule/import tests. [dep: S3d, S4b, S4c]
  Bind the addition/removal receiver pairs with the established source-preserving
  base and mutating _into behavior. Bind Molecule.enumerate_stereo as the native
  lazy enumeration operation with the same configuration and ordering as Rust,
  without exporting a transformer class. Store EnumerateStereoIter directly
  using the same wrapper pattern as ReactionApplicationIter and ReactionProductsIter.
  The native ownership contract is established in S3d; this item binds it without another
  source snapshot, prefix reconstruction, or eager collection. Review the
  Rust/Python surface against python-api.md. Map errors into the existing exception
  vocabulary; retain Transformer iterator exhaustion on rejection. Check source
  and alias preservation, in-place failure, early stop, independent yielded
  molecules, and the selected iterator's interaction with later source mutation.
  Bind StereoPerception.recognize and classify as read-only Rust delegates,
  preserving explicit algorithm arguments, selectors, result fields, errors,
  and source preservation; reconcile the existing perception methods in the
  same public-surface inventory.

  **Completed — 2026-10-30.** Bound Molecule.add_stereo/remove_stereo and their
  _into counterparts, Molecule.enumerate_stereo, and
  StereoPerception.recognize/classify. The named receiver methods use the same
  model/configuration convention as the existing transformers. Reusable models
  and configurations supply settings without separate Python transformer
  classes. enumerate_stereo(model, site_order, *, config=None) constructs the
  native StereoEnumerator and returns its owned iterator. The operation retains
  Rust defaults and contextual checks, with rejection producing an empty
  iterator. Native single-result and in-place enumeration remain Rust-only;
  Python callers select the first result with next. No transformer class is
  registered or added to the export inventory.

  Copy-producing methods clone the Rust input once and preserve source aliases
  and views. In-place methods advance the existing invalidation counter before
  native execution and preserve the molecule on rejection. The iterator wrapper
  stores EnumerateStereoIter directly; no additional source snapshot, candidate
  collection, or prefix reconstruction is introduced. It retains native order,
  multiplicity, lazy independent outputs, and exhaustion on rejection, and can
  outlive both source and transformer. Its class remains unexported and has no
  independent constructor.

  Perception queries borrow the source and retain explicit algorithm arguments.
  classify collects typed iterable sites into the native set and returns an
  independent dict for Rust's owned map. Shared perception and inconsistency
  errors retain their existing exception classes and typed attributes for
  addition/removal and perception queries. Enumeration preserves the native
  iterator's exhaustion on rejection and has no exception adapter.
  The Python API guide records all call-syntax, ownership, collection, and error
  adaptations and inventories the complete perception surface.

  Seventy-two new Python cases cover defaults, filters,
  source/alias preservation, view invalidation, error payloads, explicit order,
  fixed/finite/undetermined domains, iterator lifetime and independence, symmetry
  multiplicity, and early stop after two results from 64 binary sites. On Python
  3.13.16, the Rust binding suite passed 1,768 cases with two ignored, and the
  freshly rebuilt extension passed 3,150 Python cases with two skipped.
  All-target/all-feature binding Clippy, strict binding rustdoc, stereo-test Ruff,
  formatting, and diff checks passed.

S4 gate: the feature-gated stereo/transformer properties and focused Python
binding tests pass. Rebuild the extension before Python tests. Record every
Rust/Python surface or ownership adaptation and its specific justification;
passing tests do not settle an unexplained difference.

**S4 gate passed — 2026-10-30.** In addition to the binding suites above, all
215 graph property cases and 683 resolution conformance cases passed with
proptest and conformance explicitly enabled. The surface review found no
unexplained Rust/Python differences in the implemented stereo operations.

### S5 — Documentation, performance evidence, and closeout

- **S5a — Public contracts and examples** — additive (green).
  Modules: transformer/perception/symmetry rustdoc; docs/development/data-types.md,
  nomenclature.md, and python-api.md; applicable Rust/Python usage documentation.
  [dep: S4a, S4b, S4d]
  Document the three operations, defaults, selector
  intersections, All/Undetermined meanings, source assertion removal, and
  explicit-order absolute enumeration without deduplication. State the relevant
  laws and equality relations next to the public APIs. Correct the narrow
  determined-representation description without turning constructors or ordinary
  conversion into transformations. Examples must use the implemented happy path.
  Source comments and rustdoc must not cite this discussion record.

- **S5b — Measurements and final verification** — additive (green).
  Modules: the established transform benchmark and this record's evidence and
  lifecycle sections. [dep: S5a]
  Compare the named baseline and final cases, report measured domains and
  first/prefix/full enumeration costs, identify the dominant phases with the
  supporting profiles/counters, and investigate any unresolved regressions.
  Run the gates below once over the completed related change. Review the full
  diff for scope, public visibility, constructors, consumer migration, shared
  errors, supported-kind coverage, DSL consequences, and stale lifecycle text.
  Update completed subitems and synchronize this record's header and status
  index only after all agreed work passes; do not close independent work in 240.

## Verification commands and completion boundary

Narrow checks during the relevant stages include:

    cargo test --locked -p umol-graph-ir --lib ir::symmetry
    cargo test --locked -p umol-graph-ir --features proptest --test property -- stereo
    cargo test --locked -p umol-graph --lib -- ops::stereo:: ops::resolve::stereo:: ops::validate::stereo:: ops::transform::
    cargo test --locked -p umol-graph --features proptest --test property -- transform::
    cargo test --locked -p umol-graph --features conformance --test resolution
    cargo bench --locked -p umol-graph --bench transform -- stereo --test
    cargo clippy --locked -p umol-graph-ir -p umol-graph --all-targets --all-features -- -D warnings

Apply python-build for every command that compiles or loads umol-py. Activate
umol-py/.venv in that same shell invocation and confirm Python 3.13 before the
first build. The final gate is:

    source umol-py/.venv/bin/activate && python --version
    source umol-py/.venv/bin/activate && cargo fmt --all
    source umol-py/.venv/bin/activate && cargo test --locked --workspace
    source umol-py/.venv/bin/activate && cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
    source umol-py/.venv/bin/activate && cargo test --locked -p umol-py --lib
    source umol-py/.venv/bin/activate && maturin develop --locked --manifest-path umol-py/Cargo.toml
    source umol-py/.venv/bin/activate && pytest -q umol-py/tests
    git diff --check

The explicit property and conformance runs above are part of the final gate;
workspace defaults do not replace them. Use the active repository formatting
toolchain specified by AGENTS.md when executing the format command.

Critical path: S0 -> S1 -> S2e -> S3 -> S4 -> S5. S2a-S2d can proceed after
their listed fixture dependencies; the adder and remover implementation items
are independent once shared perception is ready. S4c can follow its config
dependencies while Rust engine work continues.

There is no required optimization or deduplication stage after the core path.
Variable-based relative stereo, symmetry deduplication, canonical/energetic
selection, and additional geometry-model research remain deferred. The first
implementation preserves the agreed absolute-configuration contract.
