---
name: python-build
description: MANDATORY — load and apply before compiling, checking, testing, linting, benchmarking, fuzzing, or installing `umol-py`, and before any workspace-wide Cargo command that includes `umol-py`. Also apply before running Python tests or rebuilding the native extension. Ensures every PyO3 compilation uses the repository's Python 3.13 virtual environment and prevents stale artifacts linked against another Python version.
---

# umol Python build environment

Read and follow [the canonical repository skill](../../../.agents/skills/python-build/SKILL.md)
in full before taking task actions. The `.agents` copy is authoritative.
