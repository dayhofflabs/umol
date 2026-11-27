# Operation semantics

This guide states shared rules once, then describes each operation's distinct behavior.
[Data type contracts](data-types.md) defines construction, ownership, and publication;
[representation integrity](integrity.md) defines valid storage and references. Chemistry models
select acceptance rules; operation configs select algorithms and execution limits.

## Shared contracts

Perception, validation, and planning leave their input unchanged. Transformations change only
the data identified by their contract; removing entities uses ordinary reference compaction.
Validators establish their named predicate and its prerequisites, not other validators' predicates.
Within that scope, an established contradiction takes precedence over uncertainty elsewhere.

### Outcomes

For operations using Result<Solution<T, C, U>, E>:

| Outcome | Meaning |
| --- | --- |
| Determined(T) | The named result is established, without a general concreteness guarantee. Validation passed its scoped checks; a derivation may still contain classified inconsistencies. |
| Underdetermined(U) | Required information remains unresolved. Consumers must inspect any established findings retained in the payload. |
| Contradictory(C) | A violation within the operation's scope is established. |
| Err(E) | The operation could not execute or accept its required input. |

### Planning and application

Transformers and individual resolvers establish all rejection conditions before mutation,
including applicability, integrity, and their stated preservation guarantees. A completed plan
applies and publishes infallibly on its unchanged input; failure there is an implementation defect.
Standalone in-place operations apply only determined plans. Composite Resolver operations may
apply phases before later rejection and use the rollback contract in
[data type contracts](data-types.md#resolution-projection-and-transformations).

## Perception and derivation

Perceivers apply no recovery policy or constraint discharge. Stereo #T/#C describe kind only;
neither they nor explicit entities assert stereogenicity. Stereo creation cutoffs restrict missing
bond creation, while stored entities are checked independently of those cutoffs.

| Operation | Contract |
| --- | --- |
| AromaticityPerceiver::derive | Combines atom assertions, incident bond evidence, and stored systems. Absent evidence excludes an atom; nonliteral contributions leave derivation unresolved. Positive requirements must be covered; stored systems must reappear with equal members and frame-aligned compatible attributes. |
| StereoPerceiver::derive | Positive #T/#C request missing entities with open cosets. Identifies Virtual as H/LP only when unique, preserving frame positions and attributes. Incomplete ligand evidence yields Underdetermined with established findings and unresolved stored rows; a missing unresolved frame supplies no row. An open coset does not itself cause underdetermination. |
| StereoPerceiver::recognize | Returns selected missing tetrahedral/CisTrans entities with open cosets. Positive assertions request the kind, negatives exclude it, and unasserted sites use automatic chemistry rules. Required nonliteral input or unsupported queries produce StereoPerceptionError; evidence failures remain in StereoDerivation. |
| StereoPerceiver::classify | Takes a site slice and supplied GraphSymmetry, coalesces duplicate requests and returns one StereogenicSiteKind per distinct requested existing site in topology-site order, including an empty vector for no requests. Supports every kind admitted by the graph-IR batch without interpreting stored cosets. Invalid sites, missing entities, incompatible kind/frame or wrong arity produce StereoPerceptionError. |
| AromaticityPerceiver::recognize | Takes ground AromaticValence values: NotAromatic excludes, Aromatic(0) permits participation. Returns disjoint systems with sorted members and aligned contributions; an empty result is valid. Negative bond assertions exclude systems containing both endpoints. Positive requirements and stored-system checks belong to derive. |

Stereo selectors intersect sites, geometry kinds and stereogenic_site_kinds: None is
unrestricted; an empty set selects nothing. The adder defaults to Constitutional and
Configurational sites. Classification uses the supplied whole-molecule auxiliary data and
returns a sorted vector keyed by StereoSite, distinct from the graph-IR batch's stereo-record
Entity keys. The method retains its perceiver receiver and model-bearing operation interface.

When site-kind filtering is active, recognition first collects all eligible prospective sites
with the existing records, classifies that context in a temporary molecule and then applies the
output selectors. Required recognition fields outside the selected output sites can therefore
cause underdetermination. Without that filter, recognition processes the selected sites and
geometries. Graph-derived ligand actions fix virtual positions, so distinguishable H positions
require explicit atoms.

Hückel counting needs contributions, not element identities: ElementScope::Any permits unknown
elements, while AllowList currently requires literal permitted elements. HMO needs known supported
elements; Clar requires carbon. Rule rejection uses Contradictory and missing parameters use Err.

### Graph symmetry and stereo readings

Molecule::graph_symmetry prepares GraphSymmetry and returns GraphSymmetryView. Queries on the
molecule take the auxiliary data explicitly; view methods use their retained data and borrowed
molecule. GraphSymmetryConfig selects the automorphism algorithm, include_constraints and
max_classification_trials. Defaults are TrailIr, false and 4096. All colors and checks read the
supplied forms; they do not copy or normalize the molecule. Entity constraints are included
only when selected, and molecule constraints are excluded. Encoding-incompatible fields use
singleton entity colors; this fallback does not diagnose or validate those fields.

is_chiral tests for a remapping to the mirror image. The proper-orbit queries use preserving
remappings; star-orbit queries also admit mirror remappings. Stereo views attach data with
with_symmetry and expose stereogenicity, topicity and ligand_symmetry. Stereogenicity is
Nonstereogenic, Prostereogenic or Stereogenic in the supplied stereo context. Topicity is
Homotopic, Enantiotopic, Diastereotopic or Distinct, with the focus record fixed. Other stereo
records are read as supplied. Stereo constraint views derive these keys as literal forms and
meet them with their own assertions for evidence. Kindless keyed readings remain vacuous;
kinded derivation requires auxiliary data.

stereogenic_site_kinds computes all kinded records together and returns the final classifications
in an owned vector, ordered by stereo-record ID. It omits kindless records and is independent of
stored coset assignments. Constitutional and Configurational sites are the literature's two
kinds of stereogenic unit; Prostereogenic and Nonstereogenic complete the classification.
The caller retains any results it wants to reuse. max_classification_trials bounds assignment
trials confirming one site's kind; exhaustion conservatively keeps the candidate as a unit,
so a bounded result can retain an unconfirmed site. It does not limit output length or the
number of enumerated stereoisomers.

AddStereoConfig, RemoveStereoConfig and EnumerateStereoConfig carry that classification limit.
The validator reads the three current-configuration stereo predicates and does not run the
site-kind batch, so StereoValidateConfig has only automorphism_algorithm and ring_config.
Validation establishes an asserted predicate when the required inputs suffice, even if an
unneeded coset is nonliteral; required unresolved evidence remains Underdetermined. This does
not promise general symbolic analysis. Kind support, element scope and ring limits belong
to the chemistry model rather than to the graph-symmetry config.

## Resolution

Resolution refines attributes and constructs compatible entities under the model. Tie-breaks select
among admissible completions; unresolved alternatives remain underdetermined. Standalone Determined
means an applicable plan, not complete validation or ground configurations.

### Refinement and conflicts

A refinement narrows admitted values while respecting stored information. Incompatible independently
meaningful descriptions are contradictory under every recovery policy. Perception distinguishes
unrealizable assertions, unrealizable entities, and mismatches; ambiguity is underdetermination.
Deleting an assertion or replacing an entity cannot serve as mismatch recovery.

### Discharge

Discharge removes an assertion only when the planned structure independently entails it. For
normalized assertion C and independent structural reading D:

| Meet | Action |
| --- | --- |
| C ∧ D is bottom | Reject the conflict. |
| C ∧ D = D | Discharge C. |
| Compatible, but C ∧ D differs from D | Retain C; it still adds information. |

D uses open derivation. None is bottom; Some(undetermined(key)) is top. Evidenced readings include
the assertion itself and cannot establish redundancy. Closing an absent relation to a negative
value cannot justify discharge. Admission may use complete evidence, but discharge may not.

The rule applies equally to new and retained entities. #T+/#C+ discharge beside an entity of the
same kind, regardless of coset; negative kind assertions conflict there. Without an entity, the
open reading stays unknown: a positive request requires construction or retention, and negatives
remain. Likewise, #d0/#t0 remain when no dative relation independently establishes them.

### Realizability failures and recovery

Recovery handles only the classified failure it names, never a mismatch:

| Action | Effect |
| --- | --- |
| Error | Return Contradictory with the diagnostic. |
| Keep | Retain the failed description. |
| Remove | Remove the identified failed assertion or entity. |

AromaticityFailurePolicy offers Error/Keep through aromatic_valence_failure and
aromatic_system_failure. StereoFailurePolicy offers all three actions through stereo_atom_failure
and stereo_bond_failure, each covering assertions and entities of its kind. Missing evidence for a
required reading yields Underdetermined under every recovery policy.

### Resolver scope

| Operation | Owned work and particular requirements |
| --- | --- |
| ValenceResolver::admit | Checks #v/#d/#t against complete derived readings on every atom. Registry and counts admission return candidates consistent with the supplied fields and constraints, preserving isotope information. Strict requires literal elements. Admit skips candidate generation for wildcard atoms. A required completion or assertion that depends on an open order is underdetermined. A symbolic element blocks admission under either policy. Established contradictions take precedence. |
| AromaticityResolver::plan_resolve | Creates systems from aromatic evidence while preserving atom fields and localized bonds. With Hückel and any element in scope, a wildcard member's #a1 can establish its contribution. A stored system keeps its electron vector. If the vector is open, member assertions supplying contributions stay in place. Missing required contributions or unresolved system charge or spin give an empty underdetermined plan. An element allow-list excludes wildcard elements. |
| StereoResolver::plan_resolve | Creates entities requested by #T+/#C+ when it can establish the frame. For a new tetrahedral frame, the plan creates `Th*` and discharges #T+. Strict yields Underdetermined with the planned edits. Admit yields Determined. The same policy applies to retained configurations. A retained kindless entity can stay open under Admit when its frame is established. A coset set or a required choice of kind stays underdetermined under either policy. Insufficient frame evidence gives an empty underdetermined plan after established conflicts are checked. Requests filtered by the creation cutoff remain. |
| IsotopeResolver::plan_resolve | WildcardPolicy::Admit preserves the isotope on every wildcard atom. On participating atoms, IsotopePolicy::Natural fills an undetermined isotope with natural composition and IsotopePolicy::Strict preserves it. Any remaining non-ground isotope on a participating atom yields Underdetermined with the proposed defaults. Standalone in-place resolution applies a determined plan. |
| BondsResolver::plan_resolve | Defaults unknown localized-bond charge to zero and wholly unknown spin to closed-shell singlet. Completes one missing spin component from a literal counterpart using the high-spin relation. Order, topology and constraints are preserved. Defaults apply even on an open order. In aggregate resolution the phase also checks localized and dative orders and noncovalent kinds against the wildcard policy. |
| MulticenterBondsResolver::plan_resolve | Checks #m against complete derived readings. Insufficient evidence yields an empty underdetermined plan. Charge and spin receive the same defaults as localized bonds. An open electron vector can survive aggregate resolution under Admit. |

Each phase discharges only its own assertions that the result independently entails, including
#v/#d/#t on already concrete atoms. Stereo preserves retained frames and configurations,
identifies uniquely determined virtual ligands, and discharges entailed #T/#C. Aromaticity
discharges entailed atom and bond assertions beside new or retained systems. Other compatible
assertions remain for validation.

### Resolver::resolve

Resolution runs isotope completion, joint valence and aromatic selection, stereo, then bond
defaults. Each phase accounts for its fields, including those on newly created entities.
An unresolved isotope or stereo configuration allows subsequent independent checks to establish
a contradiction. All four resolution methods derive their outcome from these phases.

ResolveConfig.wildcards defaults to WildcardPolicy::Strict, which requires ground inherent fields on every entity.
WildcardPolicy::Admit permits an undetermined principal field or component to stay open when the
requested work can finish. For example, with SMILES valence, `*C` resolves the carbon and bond
and retains the wildcard atom. In lenient SMILES, `C~C` remains underdetermined because its
hydrogen counts depend on bond order. For every assignment to the admitted values that full
resolution accepts under the same model, the fully resolved molecule refines the partial result.

Molecule::has_wildcards and Molecule::is_concrete are caller queries. Consuming methods return a
molecule only on Determined. In-place methods commit only on Determined and otherwise restore the
receiver under the selected rollback method. Reports record unresolved candidate sets and valence
tie-breaks. An admitted wildcard alone contributes no report entry. Conformance validation is a
separate operation. These chemistry outcomes assume the phase settings supplied at construction.
Deliberately inconsistent settings or altered intermediate results have an unspecified chemistry
outcome and must still be handled without a panic.

## Projection

Projection has the information loss and input domain stated below; representability in an external
format is a separate check.

| Operation | Contract |
| --- | --- |
| AromaticityResolver::plan_project | Replaces systems with member #a<n> and internal-bond #a assertions, meeting existing assertions. Requires one literal contribution per member, zero system charge, and closed-shell singlet spin; violations produce AromaticityProjectError. Preserves atom fields and localized bonds, but loses grouping and system-local constraints. No perception or selection; idempotent. |
| StereoResolver::plan_project | Reduces tetrahedral/CisTrans entities to #T+/#C+, discarding frames, configurations, and local constraints without interpreting them. Meets existing assertions; conflicting, unsupported, or unknown kinds produce StereoProjectError. Model-free and idempotent. Aggregate projection and export do not invoke it. |
| IsotopeResolver::plan_project | Requires a ground isotope on every atom, including wildcard atoms, or returns NonGroundIsotope with the atom ID. IsotopePolicy::Natural elides natural composition only on participating atoms. WildcardPolicy::Admit preserves natural composition on wildcard atoms. Explicit masses survive. Under IsotopePolicy::Strict projection preserves all ground isotopes. Isotope resolution under the same policies restores the input exactly. |
| ValenceResolver::project | Always returns Determined and preserves the molecule exactly. |

### Resolver::project

ProjectFlags select aromaticity → valence → isotope stages. H counts, atom charge/electron fields,
and explicit stereo entities remain. Localized and multicenter bonds must have literal zero charge
and closed-shell singlet spin even with empty flags. Success returns Determined(molecule), possibly
nonconcrete; prerequisite or phase failures produce ProjectError.

## Transformations

Transformer follows the planning/application contract above. Consuming and in-place forms agree;
rejected in-place transformations preserve the input. Iteration preserves its source and yields
independent results: deterministic operations yield one result or none on rejection; enumeration
uses the order below. Transform delegates named Molecule methods to these owners.

| Operation | Changes and conditions |
| --- | --- |
| StereoAdder | Uses recognize to add selected missing entities with open cosets. Preserves existing entities and scalar assertions; performs no discharge. Invalid selection, incomplete required input, or contradictory selected evidence reject planning. |
| StereoRemover | Removes selected entities and clears #T/#C only at their sites. RemoveStereoMode::All accepts any configuration; Undetermined requires a known kind with explicitly open coset, excluding sets and symbolic forms. Optional classification applies; creation cutoffs and automatic eligibility do not. |
| StereoEnumerator | Changes configurations of selected existing entities. Domains come from fixed/open cosets, finite sets, or supported constant expressions; variables, unsupported kinds reject preparation. The order must name each selected site exactly once; outside entries are ignored. Tries cosets ascending without symmetry deduplication; singular transformation takes the first consistent assignment or NoAssignment. |
| Aromatizer | Any stored system makes the whole operation a no-op. Otherwise perceives from localized bonds/donor pairs, lowers selected doubles to singles, transfers donor lone pairs, adds systems, and removes realized aromatic assertions. Required nonliteral transfer data, model/evidence failures, retained-assertion conflicts, or planned electron-count contradictions reject planning. Charge, H, spin, and external bonds remain. |
| Kekulizer | Localizes prescribed 0/1/2 contributions: donors regain lone pairs, acceptors retain local state, one-electron contributors require matching. Supports closed-shell systems; charge ±1 additionally requires all-one contributions and an exposed atom for charge localization. Uses the selected matching algorithm and each system atom once in the supplied order; outside entries are ignored. Rejects unsupported data, bad order, impossible matching, or count contradictions. Removes systems/aromatic assertions, retaining compatible valence assertions. Fresh aromatization need not recover the partition. |
| ChargeDelocalizer | For atom charge q and contribution e, writes e+q and zero atom charge; system charge gains Σq. Other occupations remain, preserving contributions minus system charge. Requires homogeneous literal elements, one literal contribution per member, and literal transfer charges. Other fields may be partial; ineligible systems are skipped. |
| HydrogenFolder | Folds default H on one ordinary single bond into a literal H count on a known non-H partner. Labels, nondefault state, non-stereo relations, stereo sites, and explicit molecular-constraint references protect H. Unknown partner counts skip folding. Negative eligible counts or duplicate implicit-H frame ligands reject planning. |
| HydrogenUnfolder | Expands positive literal H counts to appended default H/single bonds in bearer order and zeros those counts. Negative counts, nonliteral counts backing implicit-H stereo ligands, or count/frame disagreements reject planning; other nonliteral counts remain. Generic virtual ligands are not identified. |

Enumeration checks scalar, entity-local, and affected molecular constraints: false checks reject;
unrelated uncertainty is allowed. An assertion alone supplies no target. Inherent transform_iter
snapshots and prepares immediately, then produces candidates lazily; its owned iterator outlives
the source and transformer. No selected work yields one unchanged result if checks pass.

Hydrogen operations exchange actual/implicit-H ligands at fixed frame positions, preserving
configurations, stereo constraints, and #T/#C. They adjust localized valence/degree assertions;
total and ring assertions remain. Both are idempotent and conserve explicit-plus-implicit H on
their literal domain. Unfolding then folding is an exact inverse only when no explicit H was
initially foldable, counts are nonnegative literals, positive counts have known non-H bearers,
affected valence/degree assertions are literal or absent, and frames/assertions agree.
No general inverse is promised.

## Validation

General assertion/structure agreement belongs to ConstraintInvariantsValidator. A calculation using
combined evidence still rejects a contradictory input reading needed by its equation. None of these
operations completes fields, repairs data, selects a description, or discharges assertions.

| Operation | Predicate and required evidence |
| --- | --- |
| ValenceInvariantsValidator::validate | Electron conservation from complete evidenced localized/aromatic/multicenter valence and donated/accepted pairs. Assertions can supply counts for partial structure. Checks all required readings for bottom before requiring literal counts and atom state, then compares electrons with spin-orbital occupation. Unrelated constraints are outside its scope. |
| SpinInvariantsValidator::validate | Complete literal count/multiplicity pairs on atoms, localized bonds, aromatic systems, and multicenter bonds must form SpinState; incomplete pairs contribute uncertainty. Checks coupling payload validity, not coupling feasibility. Wholly unspecified coupling imposes no check; other unresolved coupling descriptions remain underdetermined. |

### Tier-3 conformance validators

Model-dependent predicates may differ between chemistry models. Integrity-valid input does not
establish tier-2 properties: each query establishes the prerequisites it needs and reports those
failures independently of model rejection.

### Domain checks

| Operation | Predicate |
| --- | --- |
| ConnectivityConformanceValidator::validate | Connectivity conditions enabled by ConnectivityModel. Components use localized bonds only; other relations do not join them. Overall connectivity and enabled relation/constraint containment use topology/references, so unresolved attributes do not prevent a verdict. Zero or one component is connected. |
| ValenceConformanceValidator::validate | Each atom is accepted by the selected registry/counts table. The classifier needs its ground atom state and literal complete contribution readings; unresolved required input contributes uncertainty. Other entities need not be concrete. |
| AromaticityConformanceValidator::validate | Aromaticity derivation has no classified failure or mismatch, including unrealized assertions and unreproduced stored systems. Missing model parameters use Err. |
| StereoConformanceValidator::validate | Kind agreement, realizability, arity/coset ranges, model ring limits, and asserted stereogenicity/topicity/ligand symmetry. Frame positions and permutation degrees are checked even with open cosets. Ring-relative E uses ring-continuation ligands, not a raw coset. Algorithms come from config; element scope, supported kinds and ring limits come from the model. This validator does not compute the site-kind batch. |

### Constraints and explicit entities

Stored supported stereo entities are checked with or without #T/#C; existence does not assert
stereogenicity. A positive assertion can be checked for realizability without storing an entity.
Model-independent incidence validation separately chooses Derived or DerivedComplete: #T+ without
an entity is undecided under the open reading and false under closure. Resolver discharge always
requires the independent open reading. Virtual-ligand identification for validation leaves
the stored frame unchanged; unresolved identification contributes uncertainty.

### Validator

| Method | Composition |
| --- | --- |
| Validator::validate_invariants | Constraint → valence-invariant → spin-invariant checks. Uses configured constraint algorithms and DerivedKind, defaulting to DerivedComplete. |
| Validator::validate_conformance | Connectivity → valence → aromaticity → stereo conformance; does not run invariant validation. |
| Validator::validate | Invariants, then conformance. |

Each pass continues after uncertainty to detect later contradictions; Determined requires all
participating predicates to be established. Operational errors propagate.

## External-format conversion

Boundary conversion preserves its source on every outcome. Raising establishes integrity, not
chemistry conformance; format syntax and chemical interpretation remain separate stages.

### TableIR raising

TryIntoIr copies representable boundary fields and applies fixed IO defaults, without chemistry
resolution: missing isotope → Natural, missing atom charge → zero; notation-owned H omission
supplies closed-shell count. Other missing occupations remain unresolved. Aromatic bonds become
single localized bonds plus assertions; fuzzy orders remain unknown. Explicit frames, MOL wedges,
and bond configurations become stereo entities with transported cosets; unnamed ligands remain
Virtual and raw MOL parity is ignored. Unsupported presentations or malformed entries reject
conversion. Standalone atom/bond conversion lacks topology context.

### Interpret::interpret

Raises Smiles/ReactionSmiles, then resolves under the supplied model and policy; only Determined
returns graph IR, including a partial result under Admit. Raising, contradiction, underdetermination, and execution failures stay distinct.
Reaction sides resolve independently; unique matched atom-map pairs form the correspondence,
unmatched atoms stay unmatched. Agents and repeated atoms in a map class are rejected.

### Convey::convey

Projects a private copy and converts it to the boundary, preserving atom/localized-bond
order and explicit stereo. Unsupported values, constraints, entities, generic virtual ligands,
or configurations reject conversion; #T/#C alone cannot supply stereo. H survives projection;
conversion may omit a count only when notation permits and same-model/policy inference recovers
it. Missing/ambiguous evidence and mapped atoms retain counts. Occupations remain in TableIR even
when ordinary SMILES cannot encode them independently. Rendering owns narrower syntax limits;
neither conversion nor rendering performs stereo recognition/resolution.

Reaction convey materializes/projects both sides, regenerates shared one-based labels for matched
pairs, and leaves unmatched atoms unlabeled. Agents and source-only metadata are unavailable from
graph IR; failures retain side context. Supported ingestion roundtrips under the same model/policy
use canonical equality, preserving correspondence but not spelling,
ring labels, or redundant slash markers.

### SMILES ingestion and export

ingest::ingest_smiles and ingest::ingest_reaction_smiles compose parse → interpret;
export::export_smiles and export::export_reaction_smiles compose convey → render. Defaults are
OpenSMILES, the SMILES valence preset, Natural isotope completion and Strict wildcard handling.
Explicit-config and byte-input variants retain the same composition. Unresolved ingestion stops
after parse → raise, retaining partial fields, constraints, and explicit stereo, with the same
reaction-map/agent restrictions.
Failure channels remain owned by their respective stages.
