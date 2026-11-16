# Releases

## Scope and documentation

Use the immediately preceding release tag as the baseline. Review that release's
preparation changes and the subsequent implementation changes, including breaking
interfaces, DSL changes and chemistry-data changes. Open issues and unfinished
discussion documents do not become release requirements merely because they exist.
Record explicit inclusions, deferrals and author-owned work in the release record.

Check the external whitepaper's included listings against the current checkout,
using its test driver. Review the text for changed claims and terminology. Propose
only corrections to material already discussed there; the whitepaper is not a
feature inventory. Preserve the author's manuscript-editing boundary. After the
author merges text proposals, refresh the distributed whitepaper PDF.

Review README changes separately. Keep it an introduction with useful examples;
detailed operation contracts and migration material belong in the API guides and
release notes. Present proposed README edits before applying them.

Replace RELEASE_NOTES.md with notes for the new version only, describing changes
since the preceding release and necessary migration guidance. Show a compressed
draft before writing the full notes. The GitHub workflow publishes this entire
file through `--notes-file`; do not accumulate older release sections.

## Versions and verification

Bump the workspace version, all internal dependency requirements and workspace
package entries in the root Cargo.lock together. Update any other tracked lockfiles
that pin workspace packages. Update maintained dependency examples and release
references, including the whitepaper listings. Preserve historical
discussion versions and unrelated external dependency versions. The Python project
version is inherited from Cargo through maturin.

Inspect packaging and workflow configuration against the actual workspace members
and dependency order. Keep supported Python/Rust versions, ABI settings and runner
choices unchanged unless a compatibility-policy change is explicitly selected.

Use the repository Python environment for every PyO3 build or Python check:

```sh
source umol-py/.venv/bin/activate
python --version
cargo metadata --locked --offline --no-deps --format-version 1
maturin develop --locked --manifest-path umol-py/Cargo.toml
pytest -q umol-py/tests
```

Run the whitepaper's tests/run.py with the current checkout and rebuilt Python
extension. Check that the listings themselves are exercised. Run the relevant
local checks and record results; the release workflow additionally requires the
Rust CI and Python distribution gates in .github/workflows/ci.yml and python.yml,
including all-feature tests, documentation, MSRV and installed-artifact checks.
Review the complete release diff and run `git diff --check`.

## Publication

Git mutations remain author-controlled unless expressly delegated. After the
release changes and author-owned documentation are merged, create and push the
matching `vX.Y.Z` tag. Manually dispatch .github/workflows/release.yml with that
existing tag. It verifies the version, runs the gates, packages and publishes
workspace crates in dependency order, publishes Python distributions and creates
the GitHub release with those distributions attached.

Verify the workflow result and the published crate, Python and GitHub artifacts.
Preparing files or passing local tests is not evidence that publication occurred.
