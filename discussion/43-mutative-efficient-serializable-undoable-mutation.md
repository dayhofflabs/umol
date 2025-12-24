# Mutative's approach to serializable, undoable JSON patching

# Prompt

Can you analyze Mutative (fast JSON mutation library), please, and check if the underlying mutation language, patch application methods, and data model could be a useful template for molecular manipulation? not trying to replace proper graph-derived methods (DPO or others) but looking for the engineering overlayer on top of that. This is for a molecular manipulation model, in which the molecules are stored as an attributed undirected graph. But also interested if this can be applied to the atom cloud (atom positions in 3D) data model of chemistry.

# Response

## Mutative's architecture as a template for molecular manipulation systems

**Mutative offers compelling engineering patterns—lazy proxying, structural sharing, and JSON Patch generation—that could inform molecular graph manipulation libraries, though the mapping requires significant adaptation for attributed undirected graphs and 3D coordinate systems.** The library's core insight—tracking changes via proxy interception and producing minimal patches—translates well conceptually to molecular editing, but existing domain-specific systems (SMIRKS, DPO/SPO) already solve some problems Mutative addresses differently. This analysis identifies which patterns transfer directly and where molecular-specific solutions remain necessary.

### Mutative's proxy-based draft system enables efficient change detection

Mutative implements immutable updates through a **three-phase architecture**: draft creation, mutation, and finalization. The system wraps the base state in a JavaScript `Proxy` that lazily creates nested proxies only when properties are accessed—a key optimization over eager proxy creation.

```javascript
const baseState = { atoms: [{ element: 'C', charge: 0 }], bonds: [] };
const nextState = create(baseState, (draft) => {
  draft.atoms[0].charge = 1;  // Proxy created here, not before
  draft.bonds.push({ from: 0, to: 1, order: 1 });
});
```

Each draft node internally tracks: the **original value** (base), a **modified flag**, a **shallow copy** (created on first mutation), and a **parent reference** (for path generation). During finalization, Mutative traverses only modified nodes, creating shallow copies along the path from changed leaf to root while preserving references to unchanged subtrees.

The **copy-on-write semantics** mean modifications are recorded rather than applied directly to the original. This enables both structural sharing and automatic patch generation. Mutative achieves **10-50x performance gains** over Immer primarily by: (1) not auto-freezing objects by default, (2) more aggressive lazy proxy creation, and (3) optimized finalization traversal.

### The patch system follows JSON Patch RFC 6902 with pragmatic extensions

Mutative generates patches that are **JSON Patch RFC 6902-compatible** when configured appropriately, supporting three operations: `add`, `remove`, and `replace`. The `move`, `copy`, and `test` operations from RFC 6902 are not generated.

```javascript
const [state, patches, inversePatches] = create(
  { atoms: [{ id: 0, element: 'C' }] },
  (draft) => {
    draft.atoms[0].element = 'N';
    draft.atoms.push({ id: 1, element: 'O' });
  },
  { enablePatches: true }
);

// patches:
// [
//   { op: "replace", path: ["atoms", 0, "element"], value: "N" },
//   { op: "add", path: ["atoms", 1], value: { id: 1, element: "O" } }
// ]
```

By default, paths are **arrays** (`["atoms", 0, "element"]`) rather than RFC 6902 strings (`"/atoms/0/element"`), and array truncation uses `length` property patches for performance. Full compliance requires `{ enablePatches: { pathAsArray: false, arrayLengthAssignment: false } }`.

**Inverse patches** are automatically generated alongside forward patches, enabling undo operations:

```javascript
const previousState = apply(state, inversePatches);  // Reverts changes
```

Patches compose through simple concatenation—`[...patches1, ...patches2]`—applied sequentially. There's no built-in patch optimization or compression.

### Structural sharing preserves identity for unchanged subgraphs

Mutative guarantees that **unchanged portions maintain reference equality** (`===`) with the original state. Only the path from a modified node to the root receives new object allocations:

```javascript
const base = { 
  molecule: { 
    atoms: [{ id: 0 }, { id: 1 }], 
    metadata: { name: 'ethane' } 
  }
};
const next = create(base, draft => { draft.molecule.atoms[0].charge = 1; });

next.molecule.atoms[0] !== base.molecule.atoms[0];  // Modified: new reference
next.molecule.atoms[1] === base.molecule.atoms[1];  // Unchanged: same reference  
next.molecule.metadata === base.molecule.metadata;  // Unchanged: same reference
```

This **path copying** strategy means memory overhead is proportional to modification depth, not total structure size. For a molecular graph with 10,000 atoms where one atom's charge changes, only ~log(n) new objects are created (the atom, its container, and ancestors to root).

The `original()` and `current()` helpers provide access to base and current state within the draft function:

```javascript
create(base, draft => {
  const before = original(draft.atoms[0]);  // Original atom reference
  const snapshot = current(draft.atoms);    // Current draft state snapshot
});
```

### Mapping to molecular graphs requires adapting the data model

Mutative operates on **directed, hierarchical JSON trees** where paths uniquely identify nodes. Molecular graphs are **attributed undirected graphs** with nodes (atoms) and edges (bonds) that don't naturally nest. This fundamental mismatch shapes the entire adaptation strategy.

#### Node operations translate cleanly with indexing considerations

Atom mutations map directly to object property changes:

```javascript
// Mutative-style atom editing
const molecule = {
  atoms: [
    { id: 0, element: 'C', charge: 0, position: [0, 0, 0] },
    { id: 1, element: 'N', charge: 0, position: [1.5, 0, 0] }
  ],
  bonds: [{ from: 0, to: 1, order: 1 }]
};

const [newMol, patches] = create(molecule, draft => {
  draft.atoms[0].element = 'O';      // Transmutation
  draft.atoms[0].charge = -1;        // Charge modification
  draft.atoms.push({ id: 2, element: 'H', charge: 0, position: [2.0, 1.0, 0] });
}, { enablePatches: true });
```

The challenge is **index stability**: array-index-based patches (`path: ["atoms", 2]`) break when atoms are inserted or removed earlier in the array. Solutions include:

- **ID-based addressing**: Store atoms in a Map keyed by stable IDs, generating patches like `["atoms", "atom_42", "charge"]`
- **Semantic paths**: Extend the patch format with molecular semantics: `{ op: "replace", atomId: 42, property: "charge", value: -1 }`

#### Edge operations require dual-keyed patch semantics

Bonds connect two atoms without a natural parent-child relationship. Representing bonds as an array with `from`/`to` properties creates indexing fragility. Alternative representations:

```javascript
// Option 1: Adjacency map (directed representation of undirected edges)
const molecule = {
  atoms: { /* ... */ },
  bonds: {
    "0-1": { order: 1, stereo: null },  // Canonical key: min(a,b)-max(a,b)
    "1-2": { order: 2, stereo: 'E' }
  }
};

// Option 2: Per-atom adjacency lists
const molecule = {
  atoms: {
    0: { element: 'C', neighbors: { 1: { order: 1 } } },
    1: { element: 'C', neighbors: { 0: { order: 1 }, 2: { order: 2 } } }
  }
};
```

The adjacency map approach produces cleaner patches:

```javascript
// Adding a bond
{ op: "add", path: ["bonds", "1-3"], value: { order: 1 } }

// Changing bond order  
{ op: "replace", path: ["bonds", "0-1", "order"], value: 2 }
```

#### Subgraph operations require composite patch sequences

Functional group substitutions—replacing -OH with -Cl, for instance—require atomic patches that aren't naturally expressed as single operations:

```javascript
// Replace hydroxyl with chlorine at atom 5
const patches = [
  { op: "remove", path: ["atoms", "H_bound_to_5"] },
  { op: "remove", path: ["bonds", "5-H"] },
  { op: "replace", path: ["atoms", "O_at_5", "element"], value: "Cl" },
  { op: "replace", path: ["bonds", "C-O", "order"], value: 1 }  // If needed
];
```

A molecular patch system could introduce **composite operations**:

```javascript
{ 
  op: "substitute", 
  target: { type: "subgraph", pattern: "[OH]" },
  replacement: { type: "subgraph", pattern: "[Cl]" },
  anchorAtom: 5
}
```

This mirrors SMIRKS but in a JSON-serializable, composable format.

#### Multi-step reactions compose as patch sequences with dependencies

For reaction sequences A → B → C, patches could compose:

```javascript
const reactionSequence = [
  { step: 1, name: "oxidation", patches: [...oxidationPatches] },
  { step: 2, name: "decarboxylation", patches: [...decarboxPatches] }
];

// Apply sequentially
let state = initialMolecule;
for (const step of reactionSequence) {
  state = apply(state, step.patches);
}

// Or generate combined inverse for full undo
const fullInverse = reactionSequence
  .slice().reverse()
  .flatMap(step => step.inversePatches);
```

The challenge is **patch validity across structural changes**: if step 1 removes an atom, step 2's patches referencing that atom become invalid. Solutions include:

- **Re-anchoring**: Patches reference stable IDs rather than indices
- **Validation**: Check patch applicability before application
- **Rebasing**: Like git rebase, transform patches when the base changes

### Coordinate systems benefit from selection-based operation semantics

For 3D coordinate manipulation, Mutative's individual-property-change model is inefficient. Translating 1000 atoms would generate 1000 × 3 = 3000 patches. Better abstractions:

#### Selection-based transformations

```javascript
// Proposed molecular coordinate patch format
{
  op: "transform",
  selection: { type: "residue", id: 42 },  // or atomIndices: [0, 1, 2, ...]
  transformation: {
    type: "translate",
    vector: [1.0, 0.0, 0.0]
  }
}

{
  op: "transform",
  selection: { type: "backbone" },
  transformation: {
    type: "rotate",
    axis: [0, 0, 1],
    angle: 45,
    center: [0, 0, 0]
  }
}
```

This is semantically closer to MDAnalysis transformations but in a serializable, invertible format.

#### Trajectory frame patches

For molecular dynamics trajectories, per-frame patches are impractical. Instead:

```javascript
// Frame-range operations
{
  op: "transform_frames",
  frameRange: [100, 200],
  selection: "protein",
  transformation: { type: "fit_to_reference", reference: frame0 }
}

// Or delta-encoded coordinates
{
  op: "update_trajectory",
  frameDeltas: {
    100: { atomDeltas: { 0: [0.1, -0.05, 0.02], ... } },
    101: { atomDeltas: { ... } }
  }
}
```

The structural sharing concept translates: unchanged frames reference the original data, while modified frames get new allocations.

### Comparison with domain-specific molecular systems reveals complementary strengths

| Aspect | Mutative | SMIRKS | DPO/SPO | RDKit |
|--------|----------|--------|---------|-------|
| **Change model** | Property patches | Before/after patterns | L←K→R morphisms | Direct mutation |
| **Undo support** | Built-in inverse patches | None | None | None (manual copy) |
| **Structural sharing** | Automatic | N/A | K interface defines | None |
| **Composition** | Patch concatenation | Pattern chaining | Rule composition | Method chaining |
| **Serialization** | JSON native | String notation | Graph encoding | Pickle/binary |
| **Matching** | Path-based (exact) | Substructure (fuzzy) | Isomorphism | Substructure |

**SMIRKS** excels at expressing chemical transformations compactly but lacks undo, produces ambiguous results requiring sanitization, and doesn't serialize well for transmission or storage. Mutative's patch model could complement SMIRKS by generating invertible patch sequences from SMIRKS-defined transformations.

**DPO/SPO** graph transformations are theoretically rigorous and handle the undirected graph case natively. The K interface (what's preserved) maps conceptually to Mutative's structural sharing. A molecular library could use DPO-style rules internally while exposing a Mutative-style patch API externally.

**RDKit's RWMol** provides low-level primitives but requires manual state preservation and careful index management. Wrapping RWMol operations in a Mutative-style draft/finalize pattern would add:
- Automatic structural sharing (no wasted copies of unchanged substructures)
- Patch generation for operation logging
- Inverse patches for undo

### Engineering patterns transferable to molecular manipulation

#### The draft/finalize pattern for transactional editing

```python
# Proposed API mirroring Mutative's pattern
class MoleculeEditor:
    def __init__(self, molecule):
        self.base = molecule
        self.draft = molecule.copy()  # Deep copy for drafting
        self.patches = []
        
    def __enter__(self):
        return DraftMolecule(self.draft, self.patches)
    
    def __exit__(self, *args):
        self.result = self.finalize()
        
    def finalize(self):
        # Apply structural sharing: only copy changed portions
        return self.draft.with_structural_sharing(self.base)

# Usage
with MoleculeEditor(mol) as draft:
    draft.atoms[0].charge = 1
    draft.add_bond(0, 5, order=2)
# editor.result has new state, editor.patches has changes
```

### Lazy proxy creation for large molecules

For molecules with 100,000+ atoms, eager copying is prohibitive. A proxy-based approach:

```python
class AtomProxy:
    def __init__(self, base_atom, change_tracker):
        self._base = base_atom
        self._changes = {}
        self._tracker = change_tracker
    
    def __setattr__(self, name, value):
        if name.startswith('_'):
            super().__setattr__(name, value)
        else:
            self._changes[name] = value
            self._tracker.mark_modified(self)
    
    def __getattr__(self, name):
        return self._changes.get(name, getattr(self._base, name))
```

#### Patch-based operation logging

Every molecular operation can generate a serializable patch:

```python
@generates_patch
def set_formal_charge(atom_id: int, charge: int):
    # Implementation
    ...
    return Patch(op="replace", path=["atoms", atom_id, "charge"], value=charge)
```

This enables: operation replay, undo/redo stacks, collaborative editing (operational transformation), and debugging/auditing.

#### Structural sharing for trajectory analysis

When analyzing trajectories, structural sharing prevents memory explosion:

```python
class TrajectoryFrame:
    def __init__(self, base_frame=None):
        self._base = base_frame
        self._modified_coords = {}  # atom_id -> new_coords
    
    def get_position(self, atom_id):
        if atom_id in self._modified_coords:
            return self._modified_coords[atom_id]
        return self._base.get_position(atom_id) if self._base else None
```

### Practical recommendations for implementation

**For a molecular manipulation library inspired by Mutative:**

1. **Adopt the draft/finalize transaction model** for all editing operations. This provides a clean API, automatic undo support, and operation batching.

2. **Use stable IDs rather than indices** in the internal data model. Patches referencing `atom_42` remain valid across insertions and deletions; patches referencing `atoms[3]` do not.

3. **Define domain-specific patch operations** beyond add/remove/replace: `substitute`, `fuse_rings`, `rotate_dihedral`, `transform_selection`. These compress common operations and preserve chemical semantics.

4. **Implement structural sharing at the subgraph level**. For molecular graphs, this means tracking which residues, chains, or functional groups are modified and sharing references to unchanged substructures.

5. **Generate inverse patches automatically** to enable undo without explicit state checkpointing.

6. **Consider a hybrid approach**: Use SMIRKS patterns for transformation *specification* but Mutative-style patches for transformation *execution and logging*. This combines SMIRKS's chemical expressiveness with patches' invertibility and serializability.

The key insight from Mutative is that **separating change specification (the draft) from change application (finalization)** enables powerful capabilities—undo, structural sharing, operation logging—that imperative APIs lack. This pattern transfers to molecular manipulation even though the underlying data structure (undirected attributed graphs vs. JSON trees) differs substantially.

# Prompt

Independent of the Mutative's approach now, how does the requirement for efficient serializable, composable, undoable mutations with structural sharing inform the design of the molecular internal representation? Again, applies both to molecular graphs and atom clouds. Ask clarifying questions before launching research tasks if necessary.

# Clarification

Good question—this is really about deriving representation from requirements rather than retrofitting an existing library's approach.

Before I reason through this systematically, a few clarifying questions that significantly affect the design:

**Scale & mutation patterns:**
1. What's the typical atom count range? (small molecules ~100, proteins ~10k, MD systems ~100k+)
2. Mutation frequency: rare discrete edits (reaction modeling) vs. continuous coordinate updates (trajectory manipulation) vs. both?

**Operational constraints:**
3. Is network transmission/sync a use case? (collaborative editing, client-server architecture, distributed simulation)
4. Do you need concurrent mutations from multiple "editors" that must merge?

**Implementation context:**
5. Rust-native, Python-native, or Rust core with Python bindings?
6. Should the representation support memory-mapping for out-of-core trajectories?

**Semantic requirements:**
7. For molecular graphs: do you need to preserve/track atom identity through reactions (atom mapping), or is structural equivalence sufficient?
8. For atom clouds: fixed topology with varying coordinates, or both topology and coordinates mutable?

The answers substantially change whether you'd want, say, a persistent HAMTs with content-addressing vs. a simpler copy-on-write tree vs. a log-structured approach with periodic compaction.

# Answers

1. Let's talk about the 100--1k atoms range, small to medium sized molecules, so no biomolecules
2. Mutations should describe reaction modeling, structure optimizations, so 100--1k events at most, likely much less.
3. Extreme memory or bandwidth efficiency, as for collaborative editing or distributed simulation (amazing idea!) are out of scope.
4. Concurrent mutation and conflict resolution (again, amazing thought!) are not required at the moment.
5. Rust core with Python bindings is the eventual goal. Rust-native is the short-term requirement.
6. Memory-mapping or other optimization techniques are out of scope. This should focus on data model design in the first place.
7. Yes, atom mapping is a must, pretty much all reactions rely on that, and retroactive atom mapping is possible (desirable) but currently of limited accuracy for complex transformations. Can't rely on post-hoc atom mapping, I'm afraid.
8. Both topology (implicit in atom cloud anyway) and coordinates are mutable. If atoms are swapped, atom index tracking (mapping) is a must, as for graph structures.
Please ask additional questions before proceeding if necessary.

# Clarification

Clear picture forming. A few more targeted questions:

**Coordinate handling:**
1. For structure optimizations producing many geometries—do you need the full trajectory of coordinates recoverable from patches, or just final state + ability to undo to initial? (Affects whether patches store deltas vs. absolute coordinates)

**Attribute model:**
2. What's the attribute surface area? Core chemistry (element, charge, bond order, stereochemistry) only, or extensible to computed properties (partial charges, force field types, orbital populations)? If extensible, are attributes typed or opaque blobs?

**Representation duality:**
3. Should a single `Molecule` support both graph queries (`neighbors_of(atom_id)`) and coordinate queries (`position_of(atom_id)`) through one unified structure, or are `MolecularGraph` and `AtomCloud` separate types that can be derived from each other?

**Identity scope:**
4. Are atom IDs globally unique (UUID-style) or locally unique within a molecule? The former enables tracking atoms across molecule splits/merges in reactions; the latter is simpler but loses identity when `A → B + C`.

# Answers

1. Let's assume that the full trajectory may be needed but that coalescing ops should be available. That's because for optimizations a full trajectory is rarely necessary but the approach could also be used for chain-of-states methods, where the intermediates are needed individually.
2. Let's assume that the atomic attribute set is fixed and no free-form properties are allowed. All refinement to this model are out of scope.
3. No, I am talking about two separate sets of data models (with conversion framework but that's not the point here): one for molecular graphs and one for atom clouds. No need to design of unification of manipulations of both types. Good catch, didn't specify it clearly enough before.
4. Good q. Reactions are a core application here so global atom IDs are probably better. No strong opinions on this point, can see both, depending on practicalities.

# Response

Got it. This is now a clean systems design problem. No research needed—let me reason through how each requirement constrains the representation.

---

## Core tension: Identity vs. Locality

Global atom IDs (for reaction tracking) push toward **map-based storage** (`HashMap<AtomId, Atom>`). But maps sacrifice iteration locality and make batch operations (transform all positions) slower than contiguous arrays.

For your scale (100–1k atoms), this tradeoff is acceptable. The identity requirement wins.

**Design implication:** Primary storage keyed by stable `AtomId`, with optional secondary indices (e.g., by element) built lazily.

---

## Undoability → Patches must be self-inverting

To generate inverse patches without access to base state, each patch must carry enough information to reverse itself:

```rust
enum AtomPatch {
    Add { id: AtomId, atom: Atom },
    Remove { id: AtomId, atom: Atom },  // atom data preserved for undo
    Modify { id: AtomId, attr: AtomAttr, old: AtomAttrValue, new: AtomAttrValue },
}
```

The `old` field is the key insight—without it, you need the base state to compute inverses.

**For bonds** (graph model), same pattern with canonical edge keys:

```rust
struct BondKey(AtomId, AtomId);  // invariant: .0 < .1

enum BondPatch {
    Add { key: BondKey, bond: Bond },
    Remove { key: BondKey, bond: Bond },
    Modify { key: BondKey, attr: BondAttr, old: BondAttrValue, new: BondAttrValue },
}
```

---

## Composability → Operation-level granularity, not snapshots

Patches should be individual operations, not diffs between snapshots. This enables:

- **Filtering**: Extract only bond-order changes from a reaction
- **Transformation**: Rebase patches when atom IDs change
- **Analysis**: "How many atoms were modified?"

A `MoleculeEdit` is then a sequence:

```rust
struct MolGraphEdit {
    patches: Vec<MolGraphPatch>,  // ordered
}
```

Composition is concatenation. The sequence is the source of truth.

---

## Structural sharing at this scale → Copy-on-write, not persistent structures

Persistent data structures (HAMT, RRB-vectors) give O(log n) sharing but add indirection, allocation overhead, and complexity. At 1k atoms:

| Approach | Cost of single-atom mutation |
|----------|------------------------------|
| Full clone | ~1000 entries copied |
| Persistent HashMap | ~10 nodes allocated |
| COW with Rc/Arc | Clone is free until mutation |

For 100–1k atoms with ≤1k mutations, **simple COW with `Arc<HashMap>`** is sufficient. Clone the map only when mutating:

```rust
struct MolGraph {
    atoms: Arc<HashMap<AtomId, Atom>>,
    bonds: Arc<HashMap<BondKey, Bond>>,
}

impl MolGraph {
    fn with_atom_modified(&self, id: AtomId, f: impl FnOnce(&mut Atom)) -> Self {
        let mut atoms = (*self.atoms).clone();  // COW
        f(atoms.get_mut(&id).unwrap());
        Self { atoms: Arc::new(atoms), bonds: Arc::clone(&self.bonds) }
    }
}
```

Multiple snapshots pointing to the same `Arc<HashMap>` share memory until divergence. This is structural sharing without persistent data structure complexity.

---

## Coordinates: The AoS vs. SoA question

For `AtomCloud`, you have two layouts:

**AoS (Array of Structures):**
```rust
struct AtomCloud {
    atoms: HashMap<AtomId, CloudAtom>,  // CloudAtom contains position
}
```
- Consistent with graph model
- Poor cache behavior for batch transforms

**SoA (Structure of Arrays):**
```rust
struct AtomCloud {
    ids: Vec<AtomId>,
    positions: Vec<[f64; 3]>,
    elements: Vec<Element>,
}
```
- Great for SIMD, batch ops
- Requires maintaining parallel arrays
- ID lookup needs auxiliary index

**Hybrid approach for your requirements:**

```rust
struct AtomCloud {
    index: HashMap<AtomId, usize>,        // ID → position in arrays
    ids: Vec<AtomId>,                      // dense, iteration-friendly
    positions: Vec<[f64; 3]>,
    elements: Vec<Element>,
}
```

Patches still reference `AtomId`, but application uses the index for O(1) position lookup.

---

## Coordinate patches: Granularity tradeoff

**Fine-grained** (one patch per atom move):
```rust
ModifyPosition { id: AtomId, old: [f64; 3], new: [f64; 3] }
```
- Fully undoable
- Verbose for optimizations (100 steps × 100 atoms = 10k patches)

**Semantic** (transformation-based):
```rust
TransformSelection { 
    ids: Vec<AtomId>, 
    transform: Affine3,
    old_positions: Vec<[f64; 3]>,  // for undo
}
```
- Compact
- Preserves intent ("rotate this fragment")
- `old_positions` makes it self-inverting

**Snapshot** (frame-level):
```rust
SetFrame { 
    positions: HashMap<AtomId, [f64; 3]>,
    // old positions implicit in previous frame
}
```
- Natural for trajectories
- Undo requires access to previous frame

For chain-of-states: use **snapshot patches** with explicit frame references. Coalescing an optimization trajectory means keeping only the `SetFrame` patches you want.

---

## Coalescing as patch transformation

Keep the canonical patch sequence as source of truth. Coalescing is a function:

```rust
fn coalesce(patches: &[Patch]) -> Vec<Patch>
```

Rules:
- `Add(id, A) → Modify(id, ...) → Modify(id, ...)` → `Add(id, A')` with final state
- `Add(id, _) → Remove(id, _)` → ε (empty)
- `Modify(id, attr, v0, v1) → Modify(id, attr, v1, v2)` → `Modify(id, attr, v0, v2)`

This preserves invertibility (coalesced patches still carry old values) while reducing size.

---

## Suggested representation (Rust sketch)

**Molecular Graph:**

```rust
type AtomId = u64;  // or Uuid, or newtype

struct MolGraph {
    atoms: Arc<HashMap<AtomId, Atom>>,
    bonds: Arc<HashMap<BondKey, Bond>>,
    next_id: AtomId,  // for generating new IDs
}

struct Atom {
    element: Element,
    charge: i8,
    isotope: Option<u16>,
    stereo: AtomStereo,
}

struct Bond {
    order: BondOrder,
    stereo: BondStereo,
}

enum MolGraphPatch {
    AddAtom { id: AtomId, atom: Atom },
    RemoveAtom { id: AtomId, atom: Atom },
    ModifyAtom { id: AtomId, old: Atom, new: Atom },  // whole-atom for simplicity
    AddBond { key: BondKey, bond: Bond },
    RemoveBond { key: BondKey, bond: Bond },
    ModifyBond { key: BondKey, old: Bond, new: Bond },
}
```

**Atom Cloud:**

```rust
struct AtomCloud {
    index: HashMap<AtomId, usize>,
    ids: Vec<AtomId>,
    elements: Vec<Element>,
    positions: Vec<[f64; 3]>,
}

enum CloudPatch {
    AddAtom { id: AtomId, element: Element, position: [f64; 3] },
    RemoveAtom { id: AtomId, element: Element, position: [f64; 3], former_index: usize },
    SetPosition { id: AtomId, old: [f64; 3], new: [f64; 3] },
    SetPositions { changes: Vec<(AtomId, [f64; 3], [f64; 3])> },  // batch: id, old, new
    // Optional semantic patches:
    Transform { ids: Vec<AtomId>, transform: Affine3, old_positions: Vec<[f64; 3]> },
}
```

---

## Key insight summary

| Requirement | Design consequence |
|-------------|-------------------|
| Atom mapping across reactions | Map-keyed storage with stable global IDs |
| Undoability | Patches carry `old` values |
| Composability | Operation-level patches, not state diffs |
| Structural sharing (at this scale) | `Arc<HashMap>` with COW, not HAMT |
| Serializability | Self-contained patches (no base state needed to interpret) |
| Coalescing | Transformation on patch sequence, preserves invertibility |
| Coordinate batch ops | Hybrid AoS/SoA with ID→index mapping |

# Prompt

That level of detail should suffice for now, just recording thoughts. One more q: this design makes it harder to do any graph-based operations, correct?

# Response

Yes, exactly. The `HashMap<BondKey, Bond>` representation optimizes for:

- O(1) bond lookup by atom pair
- O(1) bond insertion/removal
- Clean patch semantics (bond identity is the atom pair)

But `neighbors_of(atom_id)` becomes O(|E|)—scan all bonds. Same for degree, connected components, ring detection, subgraph matching, etc.

**Options:**

1. **Accept it.** At 100–1k atoms (~100–2k bonds), linear scans are fast enough for occasional graph queries. Simple, no sync issues.

2. **Secondary adjacency index.** Maintain `HashMap<AtomId, HashSet<AtomId>>` alongside bonds. Every `AddBond`/`RemoveBond` patch updates both. Index is derived, not serialized—rebuilt on load.

3. **Adjacency-primary storage.** Flip the design:
   ```rust
   struct MolGraph {
       atoms: HashMap<AtomId, Atom>,
       adjacency: HashMap<AtomId, HashMap<AtomId, Bond>>,
   }
   ```
   Neighbors are O(1), but bond data lives in one "direction" (need canonical access pattern or store twice).

4. **Lazy index.** Build adjacency on first graph query, invalidate on mutation, rebuild on next query. Good if graph ops are bursty.

For reaction modeling with occasional substructure queries, I'd lean toward **(2)** or **(4)**. The patch system remains clean (operates on the primary representation), and the index is a performance optimization that doesn't affect semantics.

The deeper point: **patch-friendly representations and traversal-friendly representations pull in opposite directions.** Patches want stable keys and minimal redundancy; graph algorithms want adjacency locality. You're paying for undo/compose/serialize somewhere.

# Prompt

Thinking of a proper graph storage, like petgraph, which includes adjacencies internally but exposes a clean graph API. Not looking to reinvent the wheel here.

# Response

Right. Then the design becomes a wrapper layer:

```rust
struct MolGraph {
    graph: StableGraph<Atom, Bond, Undirected>,  // petgraph does adjacency
    id_to_node: HashMap<AtomId, NodeIndex>,
    node_to_id: HashMap<NodeIndex, AtomId>,
    next_id: AtomId,
}
```

`StableGraph` specifically—regular `Graph` invalidates indices on removal.

**What you gain:** All of petgraph's algorithms (traversal, subgraph isomorphism via petgraph + vf2, connected components, etc.) work directly.

**What you lose:** Fine-grained structural sharing. `StableGraph` isn't `Arc`-splittable at the node level. Snapshots mean cloning the whole graph.

**Implications:**

| Concern | Consequence |
|---------|-------------|
| Structural sharing | Gone, or requires external snapshot storage (e.g., `Vec<Arc<MolGraph>>` for trajectory, each a full clone) |
| Patch application | Translates `AtomId` → `NodeIndex` via map, then mutates `StableGraph` |
| Serialization | Patches stay `AtomId`-based (stable); `NodeIndex` is internal, not serialized |
| Memory at 1k atoms | ~50–100 KB per snapshot; 1k snapshots ≈ 50–100 MB. Acceptable? |

If full-graph clones per snapshot are tolerable, this is the pragmatic path. petgraph is well-tested, and you're not maintaining graph internals.

The patch layer becomes purely a journaling/undo concern on top of a mutable graph, rather than the graph representation itself being patch-native. Cleaner separation, honestly.

# Prompt

Ok, that basically torpedoes the Mutative style design entirely. Thanks for the clarification. Will have to weigh the pros and cons.

# Response

Exactly. The tradeoff crystallizes to:

| | Patch-native (maps) | Graph-native (petgraph) |
|-|---------------------|-------------------------|
| Structural sharing | Yes | No (full clones) |
| Graph algorithms | Rebuild/reimplement | Free |
| Complexity | Higher (maintain invariants) | Lower (delegate to petgraph) |
| Memory (1k snapshots) | ~proportional to edits | ~100 MB |

For reaction modeling where you're likely storing tens of snapshots, not thousands, petgraph + journaling patches is probably the right call. You get proven graph infrastructure and the patches still give you undo/compose/serialize—just not COW memory sharing.

Mutative's design shines when the data is tree-shaped and snapshots are numerous. Molecular graphs are neither.