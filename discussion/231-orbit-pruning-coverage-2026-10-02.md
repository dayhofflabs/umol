# 231 — Orbit pruning for reaction spans and atom-constrained searches

Status: Completed
Date: 2026-10-02
Relates: [186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[216](216-canonicalization-performance-2026-08-30.md)

## Purpose

Before this change, orbit pruning covered molecule search below the full level only:

| Search | Topology, constitution | Structure | Full |
| --- | --- | --- | --- |
| Molecule | every adapter generator | stereo-preserving generators (doc 208) | exhaustive |
| Reaction span | exhaustive | exhaustive | exhaustive |

Exhaustive search visits every leaf of the refined partition tree, and the leaf count grows with the
automorphism group. On the identity span of the retained `symmetry_heavy_radicals` case, one
canonicalization took 13.8 ms; with a valence constraint on every hydrogen, molecule
canonicalization of the same case took 14.7 ms.

## Decision

| Search | Topology, constitution | Structure | Full, atom-local constraints | Full, other constraints |
| --- | --- | --- | --- | --- |
| Molecule | unchanged | unchanged | stereo- and atom-constraint-preserving generators | exhaustive |
| Reaction span | every adapter generator | stereo-preserving generators | stereo- and atom-constraint-preserving generators | exhaustive |

### Reaction spans below the full level

Reaction-span search prunes as molecule search does. The structure-level filter applies the molecule
rule to both sides of a span: a generator is retained when it maps every stereo atom and bond to one
whose site, ligand frame, and configuration it preserves, on the left-hand side and on the
right-hand side. Spans without stereo entities return their generators unvisited.

### Atom-local constraints at the full level

A full-level leaf key is the structure key plus constraint blocks. Constraints lie outside the
adapter colors, which is why full-level search stayed exhaustive. A molecule or span has atom-local
constraints when it has no molecule constraints, no constraints on bonds or any other non-atom
entity, and no tetrahedral stereo atom constraints. The only added block then holds atom constraint
rows keyed by atom, and each row's value does not depend on how other entities are numbered.

Such searches prune at the full level. A generator is retained when it preserves stereo
configurations, as at the structure level, and maps every atom to one with the same constraints, on
each side of a span. It then preserves the structure key and the atom constraint block, and so the
complete key.

Tetrahedral stereo atom constraints are excluded because their value is a coset relative to a
ligand frame, so equal stored values do not show that a generator preserves them. Molecule
eligibility is decided on the normalized molecule.

### Selection is unchanged

Every pruned search uses node-id branch order, and orbit representatives are orbit minima. Suppose
the first minimum leaf of the exhaustive search lay below a pruned branch `w`. Some product of
retained generators maps `w` to a smaller candidate `v` of the same cell and maps that leaf to an
equal-key leaf below `v`, which the exhaustive search visits earlier. That contradicts the choice of
the first minimum leaf. The pruned search visits a subsequence of the exhaustive leaves, so it
selects the same leaf: the same entity order, canonical value, and remapping.

## Rejected alternatives

- **Backend canonical branch order for reaction spans**, as lower-level molecule searches use. It
  visits leaves in a different order, so among equal-key leaves it can select a different entity
  order and return a different remapping.
- **Pruning under other full-level constraints.** Molecule constraints, constraints on other
  entities, and frame-relative atom constraints need their own generator conditions. That work is
  deferred.

## Evidence

The `canonicalize` benchmark measures each part on public cases, single-threaded, with Nauty. Each
value is the mean of two Criterion medians, from runs interleaved across the four commits of the
change: the benchmark, then reaction spans below the full level, then atom-local reaction spans,
then atom-local molecules. Each table compares the commit before a part with the commit that adds
it; the other parts leave its cases unchanged.

`canonicalize/reaction_span` canonicalizes the identity span of every molecule in the complete
corpus. Reaction spans below the full level:

| Case | Before, µs | After, µs | Change |
| --- | ---: | ---: | --- |
| `ordinary_naphthalene` | 108.1 | 55.8 | −48.4% |
| `disconnected_rings` | 5,981 | 78.3 | 76.4× faster |
| `overlay_heavy` | 94.4 | 92.6 | −1.9% |
| `tetrahedral_stereo` | 35.5 | 35.6 | +0.1% |
| `cis_trans_stereo_bond` | 41.1 | 41.2 | +0.4% |
| `mixed_atom_and_bond_stereo` | 72.5 | 72.5 | 0.0% |
| `frame_relative_stereo_constraint` | 38.9 | 39.0 | +0.4% |
| `meso_dichlorobutane` | 69.3 | 73.5 | +6.2% |
| `para_stereo_trichloropentane` | 93.1 | 99.0 | +6.3% |
| `para_stereo_cascade` | 799.7 | 901.7 | +12.8% |
| `feature_free_connected` | 3,086 | 80.3 | 38.4× faster |
| `feature_free_disconnected` | 3,909 | 78.2 | 50.0× faster |
| `symmetry_heavy_radicals` | 13,834 | 44.0 | 314.6× faster |

The three slower cases are structure-level spans whose stereo configurations reject every adapter
generator. Each branch point there now makes a backend request that prunes nothing, as molecule
structure search does.

`canonicalize/atom_constrained` takes each retained scaling case with a valence constraint on every
hydrogen, so canonicalization runs at the full level while the adapter graph keeps its symmetry.
Atom-local reaction spans, on the identity spans of those cases:

| Case | Before, µs | After, µs | Change |
| --- | ---: | ---: | --- |
| `feature_free_connected` | 3,797 | 106.5 | 35.6× faster |
| `feature_free_disconnected` | 4,828 | 105.2 | 45.9× faster |
| `symmetry_heavy_radicals` | 18,081 | 58.4 | 309.4× faster |

Atom-local molecules:

| Case | Before, µs | After, µs | Change |
| --- | ---: | ---: | --- |
| `feature_free_connected` | 3,344 | 130.2 | 25.7× faster |
| `feature_free_disconnected` | 4,288 | 127.1 | 33.7× faster |
| `symmetry_heavy_radicals` | 14,668 | 66.4 | 220.9× faster |

Across the whole change, every case of `canonicalize/complete/without_para_stereo` stays within
−0.4% to +1.2% of its time before the change.

`test_reaction_span_canonical_search` compares pruned and unpruned reaction-span search: symmetric
spans at the topology and constitution levels, three stereocentres with equal and with mixed
configurations, and full-level spans with a constraint on every equivalent atom, on one of them, with
a tetrahedral stereo constraint, and with a molecule constraint. Without the stereo filter, pruning
selects a different candidate for the mixed configurations; without the atom constraint check, it
does so when only one atom is constrained.

`test_canonicalize_full_with_options` compares pruned and exhaustive full-level molecule
canonicalization on the same constraint shapes and on constrained stereocentres. Without the atom
constraint check, the result differs when one of three equivalent atoms is constrained; without the
stereo check, it differs for a constrained stereocentre whose four equivalent ligands admit only odd
adapter permutations.
