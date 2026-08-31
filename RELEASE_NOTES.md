# umol 0.7.0

umol 0.7.0 is the first breaking update after the 0.6.0 public alpha release.

## Breaking changes

- `Equiv`, `RelationEquiv`, `BiRelationEquiv`, and `Normalized` have been removed, and Python
  `equiv` methods are replaced by `normalized_eq`. Participant-frame transport and equality use
  `FrameTransport`, `Reframe`, and `framed_eq`.
- Public level-selected canonicalization has been removed from Rust and Python. `CanonicalizeLevel`,
  `canonicalize_by`, and `canonical_eq_by` are gone; complete canonicalization now selects the
  lowest sufficient description level internally.
- Relation sets now preserve the supplied participant frame instead of sorting participants.
  Ordering marker types and payload-permutation traits have been removed, and relation iteration,
  coincidence, permutation, remapping, and compaction APIs have been revised accordingly.
- `Reaction` is now a closed aggregate with private fields, and its deltas are validated on
  construction. Use `new` or checked `try_new`, the `lhs` and `deltas` accessors, and `into_parts`;
  aggregate integrity checks are no longer public post-construction validators.
- Stereo frames must contain pairwise-distinct ligands and no more than `MAX_DEGREE` ligands at
  every construction and raise boundary. Repeated virtual ligands are rejected, and
  `Permutation::between_all` has been removed.
- Reaction-span construction now preserves explicitly supplied equivalent `Modified` entries;
  normalization and canonicalization collapse them to `Unchanged`.
- `IdCompaction` has been replaced by `MoleculeCompaction`; graph-core now separates generic
  `Compaction` from graph-level `GraphCompaction` and `Remapping`.

## Performance improvements

- Complete canonicalization now selects the lowest sufficient description level internally and
  prunes symmetry-equivalent search branches. In retained topology-only benchmarks, operations
  that previously required millisecond-scale full searches now complete in tens of microseconds.
- A compact topology carrier, dense initial partitioning, and allocation-reduced partition
  refinement further reduce canonicalization work, including for larger molecular graphs.

## Other compatible changes

- Molecules, reactions, and reaction spans now implement the complete normalization, reframing,
  and canonicalization pipeline, including operation-issued frame actions and correspondences.
- `Canonicalize::canonical_hash` provides an id- and participant-frame-invariant hash. Graph,
  relation, and graph-IR aggregate types also implement structural hashing.
- Frame-relative overlay data, constraints, and deltas are transported consistently through
  pushout, superposition, canonicalization, and reaction application. Compatible reaction-removal
  frames are accepted and transported instead of being rejected for differing source order.
- Checked aromatic-system, multicenter-bond, and constraint mutation is transactional in Rust and
  Python: invalid updates report an integrity error without partially modifying the molecule.
- Dative-bond integrity now follows its factor structure: an acceptor may also be a donor, and
  distinct dative bonds may share participants as long as their complete identities differ.
- `DynPermutation` adds arbitrary-degree permutation actions, while `Permutation` remains the
  bounded representation used for stereochemistry and now exposes `MAX_DEGREE`.
- Public graph-IR entity-set aggregate and frame-action types provide kind-specific iteration,
  lookup, frame selection, and transport.
