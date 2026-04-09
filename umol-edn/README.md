# umol-edn

Rust parser, formatter, and serde integration for
[EDN (Extensible Data Notation)](https://github.com/edn-format/edn).

Inspired by [edamame](https://github.com/borkdude/edamame) (Clojure),
[go-edn](https://github.com/go-edn/edn) (Go),
[toml-spanner](https://github.com/exrok/toml-spanner), and
[serde_json](https://github.com/serde-rs/json) (both Rust).

Implements the [EDN spec](https://github.com/edn-format/edn) strictly. See
[`spec/edn-spec.md`](spec/edn-spec.md) for the full specification including
ambiguity resolutions.

## Usage

### `edn!` macro

With the `macros` feature, construct `Edn` values using EDN syntax directly:

```rust
use umol_edn::edn;

let val = edn!({:name "Alice" :age 30});
assert_eq!(val.get_keyword("name").unwrap().as_str(), Some("Alice"));

let v = edn!([1 2 3 :keyword ns/sym #my/tag "value"]);
```

Without the `macros` feature, `edn!` accepts a string literal instead:

```rust
let val = edn!(r#"{:name "Alice" :age 30}"#);
```

### Value-based parsing

Parse EDN into an `Edn` value tree, then inspect or traverse it:

```rust
use umol_edn::{read_string, Edn};

let val = read_string("{:name \"Alice\" :age 30}").unwrap();
assert_eq!(val.get_keyword("name").unwrap().as_str(), Some("Alice"));
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

| Feature  | Description                                               |
| -------- | --------------------------------------------------------- |
| `serde`  | Serde `Serialize`/`Deserialize` support, streaming parser |
| `chrono` | `#inst` tag reader producing `chrono::DateTime<Utc>`      |
| `uuid`   | `#uuid` tag reader producing `uuid::Uuid`                 |
| `macros` | Proc-macro `edn!` with bare EDN syntax                    |
| `bignum` | `N`/`M` suffix for BigInt/BigDecimal                      |

## License

MIT OR Apache-2.0
