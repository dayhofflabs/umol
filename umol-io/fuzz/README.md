# SMILES parser fuzzing

`fuzz_parse_opensmiles` passes arbitrary bytes to `Smiles::parse_bytes` using the
OpenSMILES configuration. Successfully parsed values are raised to GraphIR with
`TryIntoIr`. Parse and raise errors are expected results; panics are failures.
The target does not resolve molecules or compare chemical/stereo semantics.

Run from the repository root, explicitly supplying the mutable corpus, checked-in
seeds, and dictionary:

```sh
cargo +nightly fuzz run --fuzz-dir umol-io/fuzz fuzz_parse_opensmiles \
  umol-io/fuzz/corpus/fuzz_parse_opensmiles \
  umol-io/fuzz/seeds/fuzz_parse_opensmiles -- \
  -dict=umol-io/fuzz/dictionaries/smiles.dict -max_total_time=3600
```

New inputs are written to the first corpus directory. The seed directory and
dictionary are not automatically included by a default `cargo fuzz run` command.
The time limit covers fuzz execution, not compilation. Keep the terminal output
when recording duration, execution count, coverage, and crash findings.

Replay the seeds without mutation:

```sh
cargo +nightly fuzz run --fuzz-dir umol-io/fuzz fuzz_parse_opensmiles \
  umol-io/fuzz/seeds/fuzz_parse_opensmiles -- \
  -dict=umol-io/fuzz/dictionaries/smiles.dict -runs=0
```

## Seed coverage

The 51 named seeds include both successful parse/raise paths and intentional
rejections. Seed coverage supplies entry points for mutation; it is not proof of
semantic correctness or exhaustive path coverage.

| Area | Representative seeds |
| --- | --- |
| Tetrahedral actual, implicit-H, explicit-H and LP participants | `tetra_four_explicit`, `tetra_bracket_h`, `tetra_explicit_h`, `tetra_lone_pair` |
| Tetrahedral ring encounters | `tetra_ring_open`, `tetra_ring_close`, `tetra_ring_explicit_h`, `tetra_ring_percent`, `tetra_multiple_closures` |
| Multiple stereo sites and components | `tetra_adjacent_centers`, `fused_ring_stereo`, `tetra_later_component` |
| Definite, partial and redundant cis/trans markers | `cis_trans_same`, `cis_trans_opposite`, `cis_trans_partial`, `cis_trans_four_substituents`, `cis_trans_redundant` |
| Shared direction assignments | `cis_trans_shared`, `cis_trans_shared_opposite`, `cis_trans_branched_shared`, `cis_trans_cyclic_shared` |
| Directional ring closures | `cis_trans_ring_open`, `cis_trans_ring_close`, `cis_trans_ring_both` |
| Ring topology and labels | `ring_c3`, `pct12`, `ring_label_reuse`, `ring_across_dot`, `bridged_ring`, `spiro_ring` |
| Aromatic rings and localized links | `aromatic6`, `aromatic_fused`, `aromatic_localized_link`, `aromatic_bracket_h` |
| Stereo and ring rejection paths | `cis_trans_partial_conflict`, `ring_direction_conflict`, `ring_unclosed_stereo`, `tetra_too_few_ligands`, `tetra_duplicate_virtual_h`, `parallel_bond_stereo`, `ring_self` |
| Wildcards and bracket syntax | Existing `wildcard_*` seeds |

Reaction parsing, extended/CX configuration, resolution, and stereochemical
equivalence are outside this target. Exact stereo expectations belong in unit
and property tests as well as seed inputs.

## Panic reporting

Let panics escape the target. `libfuzzer-sys` installs a panic hook that aborts
before unwinding so libFuzzer can report the crash. Its current hook also makes
an inner `catch_unwind` ineffective at hiding a panic; the redundant inner catch
has been removed to keep the target's intent explicit.
