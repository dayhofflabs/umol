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
| StereoPerceiver::classify | Returns one Symmetric/Prochiral/Stereogenic classification per requested existing site, including an empty map for no sites. Supports literal, open, or nonempty finite coset domains without choosing a coset. Invalid sites, missing entities, unsupported kind/configuration/fluxionality, or wrong arity produce StereoPerceptionError. |
| AromaticityPerceiver::recognize | Takes ground AromaticValence values: NotAromatic excludes, Aromatic(0) permits participation. Returns disjoint systems with sorted members and aligned contributions; an empty result is valid. Negative bond assertions exclude systems containing both endpoints. Positive requirements and stored-system checks belong to derive. |

Stereo selectors intersect: None is unrestricted; an empty set selects nothing. Classification
uses whole-molecule symmetry, the selected automorphism algorithm, para-stereo setting, and iteration
limit; recognition includes prospective additions when filtering by stereogenicity. Graph-derived
ligand actions fix virtual ligands, so distinguishable H positions require explicit atoms.

Hückel counting needs contributions, not element identities: ElementScope::Any permits unknown
elements, while AllowList currently requires literal permitted elements. HMO needs known supported
elements; Clar requires carbon. Rule rejection uses Contradictory and missing parameters use Err.

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
and stereo_bond_failure, each covering assertions and entities of its kind. Ambiguity remains
underdetermined under every policy.

### Resolver scope

| Operation | Owned work and particular requirements |
| --- | --- |
| ValenceResolver::admit | Checks #v/#d/#t against complete derived readings, then admits completions from the registry or counts table without selecting one. Candidates respect supplied fields, constraints, and unresolved isotopes. Nonliteral elements leave admission underdetermined. Resolver::plan_constitution applies selections and discharges entailed owned assertions, including on already concrete atoms. |
| AromaticityResolver::plan_resolve | Creates missing systems from aromatic evidence without valence admission; insufficient contributions yield an empty underdetermined plan. Preserves stored systems, atom fields, and localized bonds. Discharges entailed atom/bond aromatic assertions beside new and retained systems; other compatible assertions remain. |
| StereoResolver::plan_resolve | Creates requested entities with open cosets, identifies unique virtual kinds, and discharges entailed #T/#C. Preserves retained frames, configurations, and attributes. Requests filtered by the creation cutoff remain under every policy. Insufficient evidence yields an empty plan after established findings are checked; no configuration or stereogenicity is selected. |
| IsotopeResolver::plan_resolve | Strict preserves isotopes; Natural changes only Undetermined to Natural. Supplied masses, sets, and variables remain. Determined requires ground isotopes, not other fields. Underdetermined retains proposed defaults for aggregate resolution; standalone application waits for a determined plan. |
| BondsResolver::plan_resolve | Defaults unknown localized-bond charge to zero and wholly unknown spin to closed-shell singlet. Completes one missing spin component from a literal counterpart using the high-spin relation; other forms remain. Orders, topology, and constraints are preserved. |
| MulticenterBondsResolver::plan_resolve | First checks #m against complete derived readings; uncertainty alone yields an empty plan. Then defaults charge/spin as BondsResolver does, preserving contributions and participants. |

Each phase handles its owned keys; unrelated assertions remain for validation. Removing an entity
in a transformer is distinct from resolver discharge.

### Resolver::resolve

Runs isotope completion → joint valence/aromatic selection → stereo → localized/multicenter bond
defaults. Unresolved isotopes allow later phases to run. Determined requires Molecule::is_concrete,
including ground stereo configurations; resolution does not choose an open coset. Consuming
resolution returns a molecule only on Determined; conformance validation remains explicit.

## Projection

Projection has the information loss and input domain stated below; representability in an external
format is a separate check.

| Operation | Contract |
| --- | --- |
| AromaticityResolver::plan_project | Replaces systems with member #a<n> and internal-bond #a assertions, meeting existing assertions. Requires one literal contribution per member, zero system charge, and closed-shell singlet spin; violations produce AromaticityProjectError. Preserves atom fields and localized bonds, but loses grouping and system-local constraints. No perception or selection; idempotent. |
| StereoResolver::plan_project | Reduces tetrahedral/CisTrans entities to #T+/#C+, discarding frames, configurations, and local constraints without interpreting them. Meets existing assertions; conflicting, unsupported, or unknown kinds produce StereoProjectError. Model-free and idempotent. Aggregate projection and export do not invoke it. |
| IsotopeResolver::plan_project | Requires ground isotopes or returns NonGroundIsotope with the atom ID. Natural elides Natural composition to Undetermined; Strict preserves it; explicit masses remain. Same-policy resolution restores the input exactly. The projected result need not meet projection's own ground-input requirement. |
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
| StereoEnumerator | Changes configurations of selected existing entities. Domains come from fixed/open cosets, finite sets, or supported constant expressions; variables, unsupported kinds/fluxionality reject preparation. The order must name each selected site exactly once; outside entries are ignored. Tries cosets ascending without symmetry deduplication; singular transformation takes the first consistent assignment or NoAssignment. |
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
| StereoConformanceValidator::validate | Kind agreement, realizability, arity/coset ranges, model ring limits, and asserted stereogenicity/topicity/ligand symmetry. Frame positions and permutation degrees are checked even with open cosets. Ring-relative E uses ring-continuation ligands, not a raw coset. Algorithms/iterations come from config; chemistry limits and para-stereo acceptance come from the model. |

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
returns graph IR. Raising, contradiction, underdetermination, and execution failures stay distinct.
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
OpenSMILES, the SMILES valence preset, and Natural isotope policy; explicit-config and byte-input
variants retain the same composition. Unresolved ingestion stops after parse → raise, retaining
partial fields, constraints, and explicit stereo, with the same reaction-map/agent restrictions.
Failure channels remain owned by their respective stages.
