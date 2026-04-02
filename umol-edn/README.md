# umol-edn

Rust parser, formatter, and serde integration for
[EDN (Extensible Data Notation)](https://github.com/edn-format/edn).

Inspired by [edamame](https://github.com/borkdude/edamame) (Clojure),
[go-edn](https://github.com/go-edn/edn) (Go), and
[serde_json](https://github.com/serde-rs/json) (Rust).

## Dialects

umol-edn supports two dialects:

- **Edn** — strict interpretation of the [EDN spec](https://github.com/edn-format/edn).
  See [`spec/edn-spec.md`](spec/edn-spec.md) for the full specification
  including ambiguity resolutions.
- **Clojure** (default) — extends Edn with features present in Clojure's reader:
  `\b`/`\f` in strings, `\formfeed`/`\backspace` characters, `##NaN`/`##Inf`,
  octal string escapes, digit-start keywords, `#_` discard, `::` auto-resolve
  keywords, and lenient symbol rules.

## Usage

### Value-based parsing

Parse EDN into an `Edn` value tree, then inspect or traverse it:

```rust
use umol_edn::{read_string, read_string_with, Edn, ParseConfig, Dialect};

// Parse with default config (Clojure dialect)
let val = read_string("{:name \"Alice\" :age 30}").unwrap();
assert_eq!(val.get("name").unwrap().as_str(), Some("Alice"));

// Parse with strict Edn dialect
let config = ParseConfig { dialect: Dialect::Edn, ..Default::default() };
let val = read_string_with("[1 2 3]", &config).unwrap();
```

### Streaming serde deserialization

Deserialize EDN directly into Rust types without an intermediate value tree.
Requires the `serde` feature.

```rust
use serde::Deserialize;
use umol_edn::from_str;

#[derive(Deserialize)]
struct Person {
    name: String,
    age: u32,
}

let p: Person = from_str("{:name \"Alice\" :age 30}").unwrap();
```

### Serialization

```rust
use umol_edn::{to_string, to_string_pretty};

let edn = to_string(&vec![1, 2, 3]).unwrap();       // "[1 2 3]"
let edn = to_string_pretty(&vec![1, 2, 3]).unwrap(); // formatted with newlines
```

## Features

| Feature | Description |
|---------|-------------|
| `serde` | Serde `Serialize`/`Deserialize` support, streaming parser |
| `chrono` | `#inst` tag reader producing `chrono::DateTime<Utc>` |
| `uuid`  | `#uuid` tag reader producing `uuid::Uuid` |

## License

MIT OR Apache-2.0
