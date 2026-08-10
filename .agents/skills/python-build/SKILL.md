---
name: python-build
description: MANDATORY — load and apply before compiling, checking, testing, linting, benchmarking, fuzzing, or installing `umol-py`, and before any workspace-wide Cargo command that includes `umol-py`. Also apply before running Python tests or rebuilding the native extension. Ensures every PyO3 compilation uses the repository's Python 3.13 virtual environment and prevents stale artifacts linked against another Python version.
---

# umol Python build environment

Activate `umol-py/.venv` in the same shell invocation as every command that can compile or load the
Python extension. Tool shells do not retain activation between calls.

Use this shape from the repository root:

```sh
source umol-py/.venv/bin/activate && cargo test -p umol-py --lib
source umol-py/.venv/bin/activate && cargo test --workspace
source umol-py/.venv/bin/activate && cargo clippy --workspace --all-targets -- -D warnings
source umol-py/.venv/bin/activate && maturin develop --manifest-path umol-py/Cargo.toml
source umol-py/.venv/bin/activate && pytest -q umol-py/tests
```

Before the first compilation in a turn, verify the interpreter when it is not already visible in
the command output:

```sh
source umol-py/.venv/bin/activate && python --version
```

The interpreter must be the Python 3.13 executable inside `umol-py/.venv`. Do not add ad hoc
`PYO3_*` environment variables.

If a Rust test binary reports a missing library for another Python version, such as
`libpython3.9.dylib`, treat it as stale Cargo state rather than a code failure. Clean only the Python
package, then rebuild with the venv active:

```sh
cargo clean -p umol-py
source umol-py/.venv/bin/activate && cargo test -p umol-py --lib
```

After changing binding code or public exports, run `maturin develop` before `pytest`; otherwise
Python may exercise a previously installed extension. Commands strictly limited to crates other
than `umol-py` do not require activation.
