# 121 — Linux profiling setup (AWS EC2) for the matcher hot paths (2026-06-21)

## Why

macOS dtrace/SIP blocks sampling profilers (the `Some(42)` failures), and Docker-on-Mac
runs a VM kernel with no hardware PMU. The open performance question — how much of the
~2.3× substructure-match gap to RDKit is recoverable vs. inherent to the relational /
immutable `MoleculeAst` — needs the per-candidate inner loop measured, which is most
reliably done on Linux. Much of the upcoming perf work will need this.

**Workflow:** Claude Code and auth stay on the Mac. On the remote, just `git pull` the
repo and run the commands below by hand — no Claude Code on the remote.

## Instance choice

- **Flamegraphs (CPU time, "where does time go"):** any modern instance (e.g.
  `c7i.2xlarge`). `perf record -e cpu-clock` (software sampling) works on virtualized
  Nitro hosts.
- **Callgrind (deterministic instruction counts, "how much work per candidate"):** any
  instance — it's software instrumentation, no PMU. This is the star for the
  relational-overhead question. ~10–50× slowdown, so use a small corpus subset.
- **Hardware counters (`perf stat` cycles/instructions/cache, IPC):** require a
  **`.metal`** instance (e.g. `c7i.metal-24xl`). Virtualized Nitro instances do **not**
  expose the PMU — `perf` hardware events fail or read zero there, exactly like
  Docker-on-Mac.

For the immediate question, a plain `c7i` (or similar) with **callgrind** + **cpu-clock
flamegraphs** is enough. Reach for `.metal` only when IPC / cache-miss attribution is
needed.

## One-time setup

Amazon Linux 2023:

```
sudo dnf install -y git gcc clang clang-devel perf valgrind
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
```

Ubuntu 24.04:

```
sudo apt update && sudo apt install -y git build-essential clang libclang-dev valgrind linux-tools-common linux-tools-$(uname -r)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
```

Optional sampling tools (build flamegraphs from collected data without the browser):

```
cargo install flamegraph        # cargo flamegraph; uses perf under the hood
cargo install samply            # alternative sampler; serves a profile UI
```

perf permissions (sampling), once per boot:

```
sudo sysctl kernel.perf_event_paranoid=1
sudo sysctl kernel.kptr_restrict=0
```

### C/C++ build dependencies

`umol-graph-core` depends (non-optionally) on `nauty-Traces-sys` with the `bundled`
feature: it compiles vendored nauty/Traces C source (needs a C compiler — `gcc`) **and**
runs `bindgen` (needs **libclang** — the `clang` + `clang-devel`/`libclang-dev` packages
above). This is pulled in even when profiling `-p umol-graph`, because the automorphism
algorithm lives in `umol-graph-core`. If bindgen still can't find libclang:

```
export LIBCLANG_PATH=$(dirname "$(find /usr -name 'libclang.so*' 2>/dev/null | head -1)")
```

`umol-msym-sys` (the libmsym FFI — a git submodule, `git submodule update --init`) is
**not** on the `umol-graph` substructure path, so profiling `-p umol-graph` does not need
it; only a whole-workspace build does.

Portability is a known rough edge: the C deps (nauty via bundled-source + bindgen, libmsym
via submodule) require a C toolchain and libclang on every build host. Improving it —
feature-gating the C deps so the subiso/substructure core builds pure-Rust (nauty is used
only at runtime by `AutomorphismAlgorithm::Nauty`, never by matching), and/or checking in
pre-generated bindings to drop the bindgen/libclang build requirement — is future work,
tracked separately.

## Build with symbols

The workspace root already has a `profiling` profile (`inherits = "release"`,
`debug = 1`) — optimized with line tables. Use it for every profiling build:

```
cargo build --profile profiling -p umol-graph --example prof
```

## Ad-hoc profiling target

Not committed (keeps the crate clean). On the remote, create
`umol-graph/examples/prof.rs` with:

```rust
//! Ad-hoc profiling target for substructure matching. Not committed.
//! Args: <corpus_limit> <reps>. Callgrind: small + reps=1; flamegraph: larger.

use std::env;

use umol_ast::ast::SubstructureMatchAlgorithm::GraphAndOverlays;
use umol_ast::ast::{AtomAst, AtomId, BondAst, MoleculeAst, ValueAst};
use umol_graph::parse::parse_smiles;
use umol_graph_core::SubgraphIsomorphismAlgorithm::Vf2Rdkit;
use umol_shared::element::Element;
use walkdir::WalkDir;

fn load_corpus(limit: usize) -> Vec<MoleculeAst> {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../umol-io/tests/smiles_parsing/data/basic_opensmiles"
    );
    let mut molecules = Vec::new();
    for entry in WalkDir::new(dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "smiles"))
    {
        if molecules.len() >= limit {
            break;
        }
        let Ok(content) = std::fs::read_to_string(entry.path()) else {
            continue;
        };
        if let Some(s) = content.lines().nth(1) {
            if !s.is_empty() {
                if let Ok(m) = parse_smiles(s) {
                    molecules.push(m);
                }
            }
        }
    }
    molecules
}

fn branched() -> MoleculeAst {
    let c = || AtomAst::from_element(Element::C);
    let wb = || BondAst::new(ValueAst::Undetermined);
    MoleculeAst::from_atoms_and_bonds(
        vec![c(), c(), c(), c(), AtomAst::from_element(Element::N)],
        vec![
            (AtomId(0), AtomId(1), wb()),
            (AtomId(0), AtomId(2), wb()),
            (AtomId(2), AtomId(3), wb()),
            (AtomId(2), AtomId(4), wb()),
        ],
    )
}

fn main() {
    let mut args = env::args().skip(1);
    let limit: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(3000);
    let reps: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(20);
    let corpus = load_corpus(limit);
    let pat = branched();
    let mut total = 0usize;
    for _ in 0..reps {
        for target in &corpus {
            total += pat.substructure_matches(target, GraphAndOverlays, Vf2Rdkit).len();
        }
    }
    println!("matches={total} corpus={} reps={reps}", corpus.len());
}
```

## Callgrind — deterministic instruction counts (the relational-overhead question)

```
cargo build --profile profiling -p umol-graph --example prof
valgrind --tool=callgrind --callgrind-out-file=cg.out \
    ./target/profiling/examples/prof 300 1
callgrind_annotate --auto=yes cg.out | head -80
```

- Small corpus (300) and `reps=1` because callgrind is ~10–50× slower; instruction
  counts (`Ir`) are deterministic, so one pass is enough.
- `callgrind_annotate` gives per-function `Ir`. The number to extract: instructions per
  candidate inside `subgraph_isomorphisms` for `matches` / `host_match_targets` /
  `verify_overlays` / the `MoleculeEmbedding` build, vs. the raw search bookkeeping. That
  ratio is how much of the residual ~2× is the relational/immutable design vs. genuinely
  inherent.
- `--auto=yes` annotates source lines (needs the repo present, which it is).
- KCachegrind/`qcachegrind` can open `cg.out` for a call-tree UI if a desktop is handy.

## Flamegraph — CPU-time distribution

```
cargo flamegraph --profile profiling -p umol-graph --example prof -o fg.svg -- 3000 20
# or, manual perf + the inferno tool that `cargo install flamegraph` provides:
perf record -F 997 -g --call-graph dwarf -o perf.data -- ./target/profiling/examples/prof 3000 20
perf script -i perf.data | inferno-collapse-perf | inferno-flamegraph > fg.svg
```

Open `fg.svg` in a browser; invert (or read self-time) to confirm `meet`/`canonical` are
gone post-fix and to see what now dominates (field-wise `matches`, the host clone,
`verify_overlays`, embedding allocation).

## Optional: callgrind as a CI regression gate

`iai-callgrind` runs callgrind under a criterion-like harness for deterministic,
machine-independent instruction-count benchmarks — ideal for a perf-regression gate so a
future change can't silently re-introduce the allocation blow-up. Linux + valgrind only;
add later if perf-gating is wanted (`cargo install iai-callgrind-runner`).

## RDKit side (fair comparison)

The actual-RDKit baseline (`scripts/rdkit_substructure_baseline.py`, doc 104) runs in the
`rdkit-ref` micromamba env on the Mac. To compare on the same host, recreate the env on
the remote (`micromamba create -n rdkit-ref python rdkit`); otherwise the Mac numbers in
doc 104 stand as the reference. Note that baseline uses bare `GetSubstructMatches` (no
fingerprint screen) — for a screened comparison, add a `PatternFingerprint` pre-filter to
both sides.

## What this is meant to answer

1. Post-fix, is the per-candidate cost now field-wise `matches` CPU (relational generality)
   or something else (host clone / embedding / `verify_overlays`)? — callgrind `Ir`.
2. How many instructions/candidate vs. a hand-written fixed-struct comparison would take —
   bounds how much of the ~2× is recoverable vs. inherent to the immutable AST.
3. Whether the immutability (per-call clone, incidence rebuild) is a measurable fraction —
   guides whether a host-derivation cache or a mutable fast path is worth it.
