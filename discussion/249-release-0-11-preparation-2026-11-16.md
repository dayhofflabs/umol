# 249 — 0.11 release preparation

Status: In Progress
Date: 2026-11-16
Relates: [237](237-release-0-10-scope-2026-10-15.md),
[246](246-stereo-derived-readings-2026-11-12.md),
[248](248-explicit-import-repair-2026-11-16.md),
[release guide](../docs/development/releases.md)

## Scope and author-owned work

Prepare 0.11 using v0.10.0 as the release baseline. Review the preceding release's
preparation changes, verify the whitepaper listings against the current API,
propose corrections to existing manuscript and README claims, replace the
release notes, update versions and record a reusable release procedure.

The author retains manuscript text and README edits. The proposals below have
not been applied to main.tex or README.md. Included listing corrections are
authorized; the author separately authorized the corresponding test-harness
updates. No Git mutation or publication is part of this preparation.

The explicit import-repair proposal in discussion 248 is outside 0.11. This
release retains Python 3.9 support, Rust 1.87, the existing ABI configuration and
Ubuntu 24.04 wheel-build/test runners. Moving the Python minimum to 3.11 and the
wheel runners to latest Ubuntu is deferred to the 0.12 cycle.

## Preparation and verification

- Updated the workspace version to 0.11.0, all 41 internal dependency requirements
  across 11 member manifests, and all 19 workspace package entries in Cargo.lock.
  Locked offline metadata agrees; external lockfile package entries are unchanged.
  Also updated 21 workspace package entries across the three checked-in fuzz
  lockfiles; locked offline metadata passes for those auxiliary workspaces too.
- Rebuilt the Python extension from this checkout with the repository Python
  3.13.18 environment. The installed distribution is umol-py 0.11.0.
- The full Python suite passed: 3,577 passed and two skipped.
- The whitepaper driver passed all 12 Rust listings (13 Rust test cases) and all
  30 Python/data test cases. It tests the external paper's included files, using
  local crate patches and the rebuilt Python extension.
- All 19 current README Python blocks and all 18 blocks in the proposed shorter
  README executed successfully. The composition examples received the first,
  second and host inputs that their surrounding text explicitly assumes.
- Verified that the release workflow's 19 publication targets match the workspace
  and occur in dependency order. No workflow or Python support-policy files were
  changed.
- Nightly formatting and git diff --check passed. The release workflow's full
  Rust CI, MSRV, documentation, packaging and installed-distribution gates remain
  required before publication; they were not rerun in full during this preparation.

The v0.10.0 release preparation was checked through its tag and the version,
README, whitepaper and release-note commits 1c88931e0, ee4d0f781, 14ea2710d and
eaccdb918. The same additional responsibilities apply here: publish only the
latest release notes, keep package/dependency versions synchronized, test the
extracted examples, refresh the distributed whitepaper after author edits, and
use the manually dispatched release workflow on an existing matching tag.
The living release guide now records that procedure.

## Whitepaper listing and harness corrections applied

The manuscript repository is ~/Documents/paper/umol. Its main.tex is unchanged.

- listings/chemistry-model.py and listings/chemistry-model.rs now put phosphorus
  patterns under localized valence three and fluorine patterns under one.
- listings/nitrogen-registry.toml retains its eight illustrative atom patterns,
  regrouped under localized valences zero through four instead of charge keys.
- listings/rust-dependencies.toml now requires the 0.11.0 crates.
- tests/test_data.py checks those localized-valence groups with registry.lookup.
  Three Rust harnesses no longer set the removed CanonicalizeContext.para_stereo
  field, and the chemistry-model harness imports ValenceContradiction from
  ops::valence.

No new-feature inventory is proposed for the whitepaper. The text corrections
below address material it already discusses: the version, registry indexing,
stereo-kind versus configuration readings, entity terminology, and canonical
comparison of descriptions containing constraints. The glossary correction also
avoids implying that storing a stereo entity asserts stereogenicity.

## Proposed manuscript text changes

The following compact diff is against main.tex and is for the author to merge.
It does not add sections for the new stereo operations or other release features.
A full-context patch is also available in scratch/0.11-release/whitepaper-proposal.patch;
the proposal itself is preserved below independently of that temporary file.

```diff
--- a/main.tex
+++ b/main.tex
@@ -721 +721 @@
-the time of writing is 0.10.0; \umol\ is in alpha testing, and the interfaces described here may still
+the time of writing is 0.11.0; \umol\ is in alpha testing, and the interfaces described here may still
@@ -885,2 +885,2 @@
-Stereo atom      & site atom, ordered ligands & class, configuration \\
-Stereo bond      & site bond, ordered ligands & class, configuration \\
+Stereo atom      & site atom, ordered ligands & kind, configuration \\
+Stereo bond      & site bond, ordered ligands & kind, configuration \\
@@ -936,2 +936,3 @@
-stereo site, respectively. Noncovalent bonds provide no derivable constraints at present and are
-omitted.}
+stereo site, respectively. These constraints describe the stereo kind only; the ligand frame and
+configuration remain on the explicit stereo entity. Noncovalent bonds provide no derivable
+constraints at present and are omitted.}
@@ -949 +950 @@
-stereo      &      &      &      &      &      & class, configuration \\
+stereo      &      &      &      &      &      & kind, configuration \\
@@ -992,3 +993,4 @@
-above, only inherent fields contribute to the entity labels (colors) used in the canonical labeling.
-Since canonical labeling is a computationally expensive procedure, canonicalization is lazy, that
-is, it is only applied on request. Two readings of ethanol that differ in the order of their atoms
+above, an entity's identity is defined by its inherent fields. Canonicalization of a complete
+molecular description also accounts for its asserted constraints and transports frame-dependent
+attributes with their participants. Since canonical labeling is a computationally expensive
+procedure, canonicalization is lazy, that is, it is only applied on request. Two readings of ethanol that differ in the order of their atoms
@@ -1031,2 +1033,2 @@
-enumerates the admissible states of each element and charge as atom strings in the notation of
-\Cref{sec:dsl}. Conformance to the valence model can utilize the attribute lattice: an atom must be
+enumerates the admissible states of each element and localized valence as atom strings in the
+notation of \Cref{sec:dsl}. Conformance to the valence model can utilize the attribute lattice: an atom must be
@@ -1410,2 +1412,3 @@
-\item[stereo atom, stereo bond] Overlays that hold the configuration of a stereocenter or stereogenic
-bond with respect to an ordered list of ligands. \Cref{prin:overlays}
+\item[stereo atom, stereo bond] Overlays that hold a stereo kind and configuration relative to a
+site atom or bond and an ordered ligand frame. Their presence does not assert stereogenicity.
+\Cref{prin:overlays}
@@ -1494 +1497 @@
-\texttt{\#T} & tetrahedral stereo configuration \\
+\texttt{\#T} & tetrahedral stereo kind \\
@@ -1520 +1523 @@
-\texttt{\#C} & cis/trans stereo configuration \\
+\texttt{\#C} & cis/trans stereo kind \\
```

## README review and proposed changes

Agreed: remove the 33-line stereo/symmetry addition made after v0.10.0. It
teaches selection defaults, batch classification, auxiliary-data attachment and
lifecycle details in the introductory README. Those contracts already belong in
the linked API guide, while the release notes cover the new operations and
migration. The proposal removes the block without expanding adjacent sections.

Keep the phosphorus/fluorine registry-key correction already present in README:
it is needed for the existing chemistry-model example to execute. Apart from
removing the stereo block, the compact diff below only bumps the release version,
the preceding-release baseline and the three Rust dependency examples. A
full-context patch is available in scratch/0.11-release/README-proposal.patch.

```diff
--- a/README.md
+++ b/README.md
@@ -28 +28 @@
-Version 0.10.0 is an alpha release; interfaces may still change. See the
+Version 0.11.0 is an alpha release; interfaces may still change. See the
@@ -30 +30 @@
-for changes since 0.9.0, including API migration guidance. The
+for changes since 0.10.0, including API migration guidance. The
@@ -257,33 +256,0 @@
-
-Stereo addition creates missing tetrahedral and cis/trans entities, selecting Constitutional
-and Configurational sites by default. Removal defaults to all stored stereo and clears corresponding `#T`/`#C`
-assertions. Enumeration determines existing configurations in the supplied sequence of sites,
-preserving symmetry-equivalent assignments:
-
-```python
-from umol import GraphSymmetryConfig, Molecule, MoleculeDefaults, StereoModel, StereoSite, Stereogenicity
-
-source = Molecule.parse(
-    '{:atoms ["C" "F" "Cl" "Br" "I"] '
-    ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}',
-    defaults=MoleculeDefaults.concrete(),
-)
-model = StereoModel.default()
-stereo = source.add_stereo(model)
-isomers = list(stereo.enumerate_stereo(model, sites=[StereoSite.Atom(0)]))
-assert len(isomers) == 2
-assert stereo.remove_stereo(model) == source
-
-symmetry = stereo.graph_symmetry(GraphSymmetryConfig()).into_graph_symmetry()
-record = stereo.stereo_atoms[0].with_symmetry(symmetry)
-assert record.stereogenicity() == Stereogenicity.Stereogenic
-site_kinds = stereo.stereogenic_site_kinds(symmetry)
-```
-
-The shared selection filters topology sites, geometries and stereogenic_site_kinds. This batch
-classification differs from stereogenicity in the current assignment. Nonstereogenic and
-prostereogenic arrangements are representable; their presence asserts no classification.
-Graph-symmetry readings preserve supplied forms and use auxiliary data appropriate to the
-unchanged molecule. The [Python API guide](docs/development/python-api.md#stereo-values-and-operation-configurations)
-describes selectors, undetermined-coset removal, iterator ownership and
-[attached readings](docs/development/python-api.md#graph-symmetry-and-attached-stereo-readings).
@@ -369,3 +336,3 @@
-umol-graph = "0.10.0"
-umol-graph-ir = "0.10.0"
-umol-io = { version = "0.10.0", features = ["depiction"] }
+umol-graph = "0.11.0"
+umol-graph-ir = "0.11.0"
+umol-io = { version = "0.11.0", features = ["depiction"] }
```

## Remaining author actions

Review and merge the manuscript and README proposals, then rebuild the paper and
refresh docs/umol-whitepaper.pdf. The current distributed PDF has not been
replaced with a document that still contains the old text and version.

After the release changes are merged, the author can create and push v0.11.0 and
dispatch the release workflow. Publication and its remote verification have not
occurred. This record stays In Progress while the author-owned release work is
outstanding.
