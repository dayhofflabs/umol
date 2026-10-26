# umol 0.10.0

Changes since 0.9.0. This release expands molecular editing, revises operation
ownership and transactions, and adds hydrogen transformations and unresolved
SMILES ingestion. It makes constraint readings consistent across normalization,
separates exact and cascading removal, and reduces matching allocations. It
includes breaking Rust and Python API changes.

## Molecular editing and ownership

- Expose direct mutation of entity attributes and constraints through Rust mutable
  views and Python properties. Expand structural edits and add batched field and
  constraint updates across entity kinds. Attribute mutation remains infallible;
  structural editing checks representation integrity before publication.
- Rust Molecule.apply consumes the molecule and Edits, returning the modified
  molecule on success. Failure returns an error and drops the consumed input.
  Python Molecule.apply returns a new molecule and preserves the source on both
  success and failure. The supplied Edits is consumed in either language.
- Molecule.apply_into replaces the receiver only on success. MoleculeEditor
  supports incremental structural construction, with finish checking and
  publishing the result. Ordinary edit and reaction application no longer need
  an undo journal.
- Separate Transaction and TrackedTransaction provide rollback where requested,
  with Clone and Journal strategies. Only tracked transactions maintain the
  correspondence returned at commit. The composite Resolver retains transactions
  for in-place resolution and projection.

## Constraint readings and ring views

- Entity constraint views distinguish asserted, derived, and evidenced readings.
  Asserted readings reflect stored assertions; derived readings use molecular
  relations; evidenced readings combine applicable evidence by meet. Results do
  not depend on whether the stored constraints have been normalized.
- Readings use None consistently for bottom, meaning contradictory evidence.
  Absent evidence produces an undetermined form. The _complete counterparts give
  closed readings for resolution-written overlays, including definite negative
  aromaticity and stereo information where appropriate.
- Extend Python access to constraint readings and ring data, including RingView,
  RingAtomView, RingBondView, RingRelation, RingConnection, and RingGraph.
  RingSet remains standalone perceived data; molecule-bound RingViews can reuse
  it without repeating ring perception.

## Transformations

- Add HydrogenFolder and HydrogenUnfolder for explicit and implicit hydrogen
  representations. Folding retains isotope-labelled and otherwise ineligible
  hydrogens. Both operations transport affected stereo frames and constraints;
  retained atoms preserve total valence, degree, and hydrogen readings.
- Add the Rust Transform extension trait with Molecule.aromatize, kekulize,
  delocalize_charge, fold_hydrogens, and unfold_hydrogens, plus _into counterparts.
  Python exposes the same named molecule methods.
- Rust value-returning transformations consume their input; Python counterparts
  preserve the source and return a new molecule. The _into methods mutate the
  receiver and leave it unchanged on rejection. Each operation publishes its own
  result; failure in a later operation does not undo earlier successful calls.

## Removal and reaction application

- Separate exact topology removal from explicitly cascading removal in editors
  and Edits. Exact removal rejects retained bonds and overlays that reference
  removed entities before mutation. Cascading removal also removes incident
  bonds and overlays whose stored participants or sites disappear.
- Check reaction dangling conditions before applying changes, including bond
  dependencies of aromatic systems and stereo entities. A match that would leave
  a retained dependent entity dangling is inapplicable; reaction application does
  not repair it by cascading removal.
- Align stereo identity and correspondence handling between application, realized
  reactions, and spans. Changes between different determined stereo kinds are
  represented as removal and addition in reaction spans.

## Resolution and chemistry corrections

- Individual resolvers, projectors, and transformers prepare complete changes
  before applying them. Resolution preserves the molecule-level constraint list
  instead of automatically placing or discharging those constraints.
- Make resolution reports optional through _with_report methods. Solution carries
  a molecule only for a determined result; undetermined outcomes no longer expose
  a partially resolved molecule.
- Correct heteroaromatic kekulization, aromaticity perception under non-aromatic
  assertions, and counts-based valence bounds. Stereo and aromatic-system
  comparisons account for participant ordering by reframing before comparing
  their values.

## Unresolved SMILES ingestion

- Rust adds ingest::ingest_smiles_unresolved and
  ingest::ingest_reaction_smiles_unresolved, with _with variants accepting an IO
  configuration. They parse and raise the format representation without running
  molecular resolution.
- Python Molecule.from_smiles and Reaction.from_reaction_smiles accept
  resolve=False. Resolution remains enabled by default. Unresolved ingestion
  retains format defaults and raised information, including stereo overlays;
  it leaves unspecified molecular properties unresolved.

## Storage, matching, and parsing

- Reduce matching allocations through reusable buffers and revised incidence
  construction. Extend shared storage for molecule constraints, reaction deltas,
  and reaction spans to reduce copying.
- Replace the mutex-protected stereo coset registry with lazy initialization per
  class key. Unused symmetry classes are not constructed at startup.
- Fix parsing and roundtrip failures found by fuzzing, including escaped token
  delimiters and excessive work on malformed input. EDN consistently recognizes
  Unicode White_Space plus comma, including form feed, and escapes delimiters
  when writing tokens.

## Migration from 0.9.0

- Update Rust callers for consuming molecule operations. Clone explicitly when
  retaining the source, or use the applicable _into method for in-place work.
  Python value-returning methods preserve their receiver, but transferred Edits
  are consumed even on failure; copy a batch explicitly if it must be reused.
- Use TrackedTransaction when commit must return a correspondence, and select
  the rollback strategy explicitly in Rust. Ordinary application does not need
  to be wrapped in a transaction.
- Replace remove_topology with remove_topology_cascading, or Edit::RemoveTopology
  with Edit::RemoveTopologyCascading, where the previous cascade is intended.
- Handle owned Option results from constraint readings: None is contradiction,
  not absence. Replace normalized with consuming normalize, or normalized_eq
  when only normalized equality is needed.
- Rename DelocalizeCharge to ChargeDelocalizer. Import Transform for the named
  Rust molecule transformation methods.
- Use resolution _with_report methods when diagnostics are needed. Python returns
  the report separately from Solution. Underdetermined outcomes carry no molecule.

umol remains an alpha library. Canonical representatives may change between
0.x releases and are not persistent identifiers.
