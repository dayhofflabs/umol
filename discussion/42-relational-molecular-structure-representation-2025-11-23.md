# Data structures for efficient linting of molecular IRs

# Prompt

I am curious if it would make sense to generate a "redundant" tree, where each atom and bond belongs to >= 1 fragments (e.g., functional group, tetrahedral center, ring, ...). are there precendents for "overcomplete" IR representations -- does not have to be the single source of truth, can be derived from the latter for the purpose of linting?

The underlying data model of CodeQL is explicitly **relational** and uses **Datalog**, a declarative logic programming language. This approach is fundamental to how it handles "overcomplete" or redundant data: it doesn't store a single tree but rather a collection of tables (relations) that describe the code.

# Response (Google Gemini 3 Pro via GitHub Copilot)

### The CodeQL Data Model: Relational & Logic-Based
In CodeQL, source code is treated as data. The "database" is literally a relational database containing tables extracted from the source code during the build process.

1.  **Relational Algebra (The "Tables"):**
    CodeQL extracts code into primitive tables. For example, in a Java project, you might have tables like:
    *   `expressions(id, type_id, parent_id)`
    *   `methods(id, name, class_id)`
    *   `calls(expression_id, target_method_id)`

    Crucially, this is **set-based**, not tree-based. An entity (like a variable) is just an ID that can appear in dozens of different tables simultaneously. This is the "overcomplete" representation you are looking for. A variable isn't "inside" a function in the way a node is inside a DOM tree; rather, the `variable_in_function` relation *links* them.

2.  **Datalog (The Query Language):**
    You write queries in QL (a Datalog dialect) to infer new relations from existing ones.
    *   *Logic:* "A method is 'publicly accessible' IF it is marked public OR if it is an interface method."
    *   This defines a new, virtual relation `publicly_accessible(method_id)` that is derived on the fly (or cached).

3.  **Data Flow as Relations:**
    Complex analyses like taint tracking (e.g., "does user input reach this SQL query?") are modeled as recursive queries over these relations.
    *   `flow(source, sink)` is a relation calculated by recursively joining `assignment(a, b)`, `call_argument(func, arg)`, etc.

### Why this fits `umol`
Using a relational/Datalog model for chemical structures avoids the "smooshing" problem of fingerprints. You retain 100% reversibility and precision because you are storing exact indices, not hashes.

*   **Table `Atoms`:** `{id: 1, symbol: "C"}`
*   **Table `Rings`:** `{ring_id: A, atom_id: 1}`, `{ring_id: A, atom_id: 2}`...
*   **Table `AromaticSets`:** `{set_id: X, atom_id: 1}` ...

Atom 1 is fully defined by the intersection of all tables containing `atom_id: 1`.

### Relevant Links & Inspiration
*   **CodeQL Overview**: [About CodeQL](https://codeql.github.com/docs/codeql-overview/about-codeql/) - Explains the database generation and query execution.
*   **Soufflé (Datalog Engine)**: [Soufflé Logic](https://souffle-lang.github.io/) - A high-performance Datalog engine often used for static analysis. It compiles Datalog relations into highly optimized C++ parallel code. This might be a great backend or inspiration for `umol`'s query engine.
*   **Rust's "Polonius"**: [Polonius Borrow Checker](https://github.com/rust-lang/polonius) - The Rust compiler team is moving towards a Datalog-based borrow checker. It models lifetimes and loans as relations (`origin_live_on_entry(origin, point)`).

For `umol`, you don't need a full SQL engine. A simple in-memory "Entity Component System" or a set of `map[string][]int` (Index Maps) as discussed previously effectively implements a lightweight relational database tailored for your graph.

# Prompt

Can you give some more background on ECS? How does this model differ from a relational? A set of flat indexes for specific fragment type makes sense but seems too restrictive for the case that one would want to create novel fragment types. Is that correct?

# Response

### ECS vs. Relational: The "Database of Objects"

To answer your question directly: **No, a set of flat indexes is not inherently restrictive if designed as a dynamic map.** The restrictiveness comes only if you hardcode the *types* of indexes (e.g., `struct { RingIndex, StereoIndex }`). If you implement it as `map[string]Set<ID>`, it becomes an infinitely extensible "Tagging System," which is effectively a dynamic ECS.

Here is the breakdown of ECS and how it relates to Relational models for your use case.

---

### 1. What is ECS (Entity-Component-System)?
ECS is an architectural pattern primarily used in high-performance game development (Unity DOTS, Bevy Engine) to solve the exact problem you have: **objects (Entities) need to be composed of arbitrary combinations of features (Components).**

*   **The Entity:** Just an Integer ID (e.g., `AtomID: 42`). It has no methods and no data. It is just a key.
*   **The Component:** A raw data bucket.
    *   `ElementComponent`: `{ Symbol: "C" }`
    *   `RingComponent`: `{ RingSize: 6, IsAromatic: true }`
    *   `ChiralComponent`: `{ Type: "R" }`
*   **The System (The Linter Rule):** Logic that runs on specific intersections.
    *   *AromaticitySystem:* "Iterate over all Entities that have `RingComponent` AND `ElementComponent`."

#### How it differs from Relational
They are conceptually cousins. You can think of ECS as a **high-performance, in-memory relational database** with one specific constraint: **The Primary Key for every table is always the Entity ID.**

| Feature | Relational (CodeQL / SQL) | ECS (Bevy / Unity) |
| :--- | :--- | :--- |
| **Relationship** | Arbitrary (Table A can link to Table B via any column) | **Star Schema**: Everything links back to the Entity ID. |
| **Storage** | Row-based (usually) | **Column-based** (Structure of Arrays). All `RingComponents` are stored in one contiguous array. |
| **Query** | `JOIN` tables on arbitrary keys | **Filter** entities by component presence (set intersection). |

### 2. Is strict ECS or Indexing "Too Restrictive"?
You asked if flat indexes are too restrictive for **novel fragment types** (e.g., a user defines a custom "Warhead" fragment for drug discovery).

**If you use a Static Schema:** Yes, it is restrictive.
```go
// Rigid struct - Bad for dynamic fragments
type AtomData struct {
    Rings  []RingInfo
    Stereo []StereoInfo
    // If a user wants to add "WarheadInfo", you have to recompile the code.
}
```

**If you use a Dynamic/Generic Schema:** No, it is extremely flexible.
This is where the "Tagging" or "Labeling" approach (a simplified, dynamic ECS) shines.

#### The Recommended Model: Dynamic Tagging
Instead of typed components for every possible chemical feature, you use a generic **Label** or **Tag** system. This is essentially a Relational model where the relation name is a string.

**Data Structure:**
```go
type LinterDB struct {
    // The "Entity" is just the index in these arrays/maps
    Atoms []Atom 
    
    // The "Dynamic Index" (The Overcomplete Representation)
    // Key: Fragment Type / Scope Name
    // Value: List of Atom IDs belonging to this fragment
    Fragments map[string][]int 
}
```

**Scenario: Adding a Novel Fragment Type**
1.  **User Config:** A user defines a new rule in YAML: "Flag all *nitro groups*."
2.  **Parser/Recognizer:** You run a SMARTS matcher for `[N+](=O)[O-]` against the graph.
3.  **Storage:** You don't need a new struct. You just write to the map:
    ```go
    // Dynamically add a new "column" to your database
    db.Fragments["group.nitro"] = []int{12, 13, 14}
    ```
4.  **Linting:** The rule queries the index:
    ```go
    // Rule checks intersection of 'group.nitro' and 'ring.aromatic'
    matches := Intersect(db.Fragments["group.nitro"], db.Fragments["ring.aromatic"])
    ```

### 3. Why this fits `umol` (The Verdict)
You should use a hybrid approach that borrows the **Querying** power of Relational/ECS but keeps the **Dynamic** nature of Maps.

*   **Don't use strict ECS structs** (like `type RingComponent struct...`) for high-level chemical features, because defining new struct types at runtime is hard/impossible in compiled languages like Go/Rust.
*   **Do use a "Fragment Index" Map.** This *is* your relational database.
    *   Table Name = Map Key (`"ring.benzene"`, `"stereo.tetrahedral"`)
    *   Rows = Map Value (`List<AtomID>`)

This gives you the "overcomplete" representation you wanted:
*   Atom 5 is in `db.Fragments["element.C"]`
*   Atom 5 is in `db.Fragments["ring.aromatic"]`
*   Atom 5 is in `db.Fragments["user_defined.warhead"]`

This is exactly how search engines (Inverted Indexes) and tag-based file systems work. It is maximally flexible for "novel types" because a type is just a string key.

# Prompt

For hash map based inverted indexes, no additional tooling is needed. Are there efficient Datalog implementations in Rust that one could use here, e.g., related to the new borrow checker?

# Response

Yes, there are highly efficient Datalog implementations in Rust that align perfectly with your "overcomplete" relational model.

The two most relevant libraries are **Crepe** and **Datafrog**. Both were born from or inspired by the Rust compiler team's efforts to rewrite the borrow checker (Polonius) using Datalog.

### 1. Crepe (Recommended)
**Crepe** is likely the best fit for `umol` because it allows you to write Datalog rules directly in Rust using a macro, generating high-performance code at compile time.

*   **How it works:** It compiles Datalog rules into efficient, specialized Rust iterators. It feels like writing SQL/Datalog inline.
*   **Why it fits:** You can populate the initial "Relations" (your recognizer output) and then run the engine to derive complex chemical properties.

```rust
use crepe::crepe;

crepe! {
    @input
    struct Atom(i32, &'static str); // (ID, Symbol)

    @input
    struct InRing(i32, i32); // (RingID, AtomID)

    @input
    struct IsAromatic(i32); // (RingID)

    @output
    struct AromaticAtom(i32); // Derived relation

    // Datalog Logic:
    // "An atom is aromatic if it is in a ring, and that ring is aromatic."
    AromaticAtom(a) <- Atom(a, _), InRing(r, a), IsAromatic(r);
}

fn main() {
    let mut runtime = Crepe::new();
    // Phase 1: Enrich (Populate tables)
    runtime.extend(&[Atom(1, "C"), Atom(2, "N")]); 
    runtime.extend(&[InRing(101, 1), InRing(101, 2)]);
    runtime.extend(&[IsAromatic(101)]);

    // Phase 2: Lint (Run query)
    let (aromatic_atoms,) = runtime.run();
    // aromatic_atoms contains {1, 2}
}
```

### 2. Datafrog
**Datafrog** is the lightweight engine powering "Polonius" (the new borrow checker). It is lower-level than Crepe and focuses on fast, incremental joins.

*   **Use Case:** If you need maximum performance and don't mind writing explicit "joins" in Rust code (building the query plan yourself) rather than using declarative syntax.
*   **Relevance:** It is designed for the exact problem of "static analysis on relations" where the relations are dense integers (like Atom IDs).

### 3. Differential Dataflow
If you need **incremental updates** (e.g., a user edits the SMILES string in a GUI and you only want to re-compute the minimal set of changes), **Differential Dataflow** is the industry standard. It is complex but allows your linter to be reactive.

### Summary Recommendation
*   **Use Crepe** for a clean, declarative implementation of your "Rules." It maps 1:1 with the "Tagging/Overcomplete" mental model.
*   **Data Model:**
    *   **Input Relations:** `Atom(id, type)`, `Bond(src, dst, type)`, `Fragment(frag_id, type)`, `InFragment(frag_id, atom_id)`.
    *   **Rules:** Written in Datalog to derive `Violation(atom_id, "reason")`.

This approach gives you the query power of a database with the zero-overhead performance of Rust iterators.