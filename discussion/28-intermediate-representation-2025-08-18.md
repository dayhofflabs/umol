# Summary

* Informational

# Prompt (Gemini 2.5 pro)

Ok, here is another question. Building the progressive parser and comparing the implementation with the actual parser seems to reveal a typical parsing pattern, in which the parsing procedure generates some intermediate representation, which can be processed in different ways. Here, it can be either used to build an in-memory molecule model or formatted for output. The accumulation is currently coupled to the parsing, and it's probably best to leave it as is. But I do wonder what a good pattern here would be, just as a thought exercise. I am primarily thinking of an object that could be used for semantic validation (in particular valence rules but possibly other aspects) *before* building the molecule. That would ensure that every Molecule or MoleculeStandard object is a valid molecule. RDKit uses the opposite approach, in which the Mol object may or may not be valid, depending on sanitization, and would randomly crash on some operations. The more I think about this, the more useful the intermediate representation appears to. If I understand it correctly, it is the equivalent to an AST in parsing program code. Now, for AST, is it necessary that the input code is valid? In other words, is the validation happening at the AST stage or prior to it? If there is an intermediate step, how is it usually structured? It also seems that this task is different from that handled by parser combinators. Is there some tooling for this task (for example from interpreter building toolkits) or is it so performance-critical that it is usually coded manually? I know that in yacc/bison and Antlr4, one typically has to call processing functions during parsing. Is this the only approach available? How are ASTs / IRs manipulated in some applications? Please give some background (I know very little about this) and reason through some ideas that could help us make MOL parsing more robust. No code is needed at this point.

# Response

That is an absolutely fantastic line of questioning. You've gone from the tactical details of parsing to the core architectural strategy of building robust scientific software. Your analogy to an AST (Abstract Syntax Tree) in compilers is exactly right, and it's the key to answering all your questions.

Let's break this down as a thought exercise.

### The AST and the "Two-Phase" Model

Your intuition is correct: the problem is that the current process combines two distinct tasks:
1.  **Syntactic Parsing:** Does this text *look* like a MOL file? Are the columns aligned? Are there the right number of lines? This is what `nom` excels at.
2.  **Semantic Validation:** Does this syntactically correct MOL file describe a *chemically plausible* molecule? Does Carbon have 8 bonds? Is this stereocenter defined correctly?

Compilers solved this problem decades ago. They almost never go directly from source code to machine code. They use a "two-phase" (or multi-phase) model, and the AST is the bridge between the phases.

**To answer your core question:** An AST represents the *structure* of the code, not its *correctness*. It is absolutely possible, and very common, to generate an AST for code that is syntactically perfect but semantically nonsensical (e.g., `let x: String = 5;` in Rust).

*   The **Parser** builds the AST. Its only job is to enforce the rules of the language's grammar. If the input doesn't follow the grammar, the parser fails.
*   The **Semantic Analyzer** (or "Validator" or "Type Checker") walks the AST. Its job is to enforce the *meaning* and *rules* of the language. This is where it would catch the `String = 5` error.

This separation is the key to robustness.

### A Pattern for Robust MOL Parsing: The IR

Let's apply this to our MOL parser. Instead of parsing directly to a `Molecule` object, we would introduce an Intermediate Representation (IR), which is our "MOL AST".

**Phase 1: Parsing to the Intermediate Representation (IR)**

The parser's only job is to create a "dumb" IR. This IR is a direct, one-to-one mapping of the data in the file. It's just a bag of structs.

```
// This is a conceptual struct, not real code
struct RawCtab {
    atoms: Vec<RawAtom>, // Contains just x,y,z, symbol, charge from file
    bonds: Vec<RawBond>, // Contains just atom1_idx, atom2_idx, type from file
    properties: Vec<PropertyEntries>,
}
```
The parsers we've built are already perfect for this. Their job is to successfully turn the text into this `RawCtab`. At this stage, we don't care if an atom has a valence of 99. We only care that the line was formatted correctly.

**Phase 2: Semantic Validation**

This is the new, crucial step. We create a `Validator` that takes the `RawCtab` as input.

```
// Conceptual code
struct Validator {
    // could hold valence rules, etc.
}

impl Validator {
    fn validate(&self, raw_ctab: &RawCtab) -> Result<ValidatedCtab, Vec<SemanticError>> {
        let mut errors = Vec::new();

        // Pass 1: Check atom valences
        for (i, atom) in raw_ctab.atoms.iter().enumerate() {
            let expected_valence = self.get_expected_valence(atom.symbol);
            let actual_valence = calculate_valence_from_bonds(i, &raw_ctab.bonds);
            if actual_valence > expected_valence {
                errors.push(SemanticError::InvalidValence { atom_index: i, ... });
            }
        }

        // Pass 2: Check stereochemistry, ring definitions, etc.
        // ...

        if errors.is_empty() {
            Ok(ValidatedCtab(raw_ctab)) // A newtype wrapper to prove it's valid
        } else {
            Err(errors)
        }
    }
}
```
This is a game-changer for diagnostics. Instead of getting a single, cryptic "parse failed," the user gets a list of specific chemical errors: `["Error: Atom 5 (Carbon) has valence of 5", "Error: Stereocenter at atom 12 is defined ambiguously"]`.

**Phase 3: Construction**

The `Molecule` and `MoleculeStandard` structs would now have a constructor that only accepts a `ValidatedCtab`.

```
impl Molecule {
    fn from_validated(validated_ctab: ValidatedCtab) -> Self {
        // ... build the petgraph graph and high-level Molecule struct ...
    }
}
```
This construction step is now guaranteed to succeed. Any `Molecule` object that exists in the program is, by definition, semantically valid. The RDKit problem of having an invalid `Mol` object that might crash is completely eliminated.

### Answering Your Specific Tooling Questions

*   **Parser Combinators vs. This:** Parser combinators (`nom`) are the perfect tool for **Phase 1**. The validation logic of **Phase 2** is entirely separate and is pure, domain-specific Rust code (chemistry rules). It's not a parsing task.
*   **Tooling (Yacc/Bison/Antlr):** You are right that these tools often encourage embedding semantic actions directly in the grammar rules (e.g., `NUMBER { check_valence($1); }`). This is the "one-phase" model. The "two-phase" model (Parse to IR, then Validate) is generally considered superior for complex systems because it's far more modular and testable. You can unit-test your validator completely independently of the file parsing logic. `nom` naturally encourages this better two-phase pattern because its parsers return data structures, not side effects.
*   **Performance:** This is not a performance concern. The text-parsing part (Phase 1) is the slow part, and our existing parsers are already optimized for that. The validation and construction phases happen on in-memory data structures and are incredibly fast.
*   **Manipulation:** The IR is incredibly powerful. The Validator is just one "pass" that can operate on it. You could write other passes:
    *   A **Sanitizer Pass:** Automatically fixes common errors (e.g., adds explicit hydrogens).
    *   A **Canonicalizer Pass:** Reorders atoms to produce a canonical representation.
    *   A **Feature Extractor Pass:** Calculates properties like ring counts or functional groups.

By separating the stages, you create a pipeline where these independent, testable components can operate on a common data structure. This is a hallmark of robust, extensible scientific software.

# Prompt

Ok, this does seem a big improvement compared to the current approach, which resembles RDKit's parse to molecule, then validate approach. There is also the tempation to generate the same AST from SMILES/SMARTS strings, which are the next feature to implement. I supoose there are two competing issues here. One is trying to unify the semantics of graph-based molecule representations, which would align well with the name and purpose of this crate (umol-models-graph). The other is wanting to keep the subtle differences between the semantics of the SMILES and MOL representations and their quirks apart. I am just a bit worried that trying to build the grand unified theory of graph-based molecule representations will represent neither format faithfully, creating the problems discussed in @02-design-summary-2025-02-16.md .  On the other hand, I want to be pragmatic and the advantage of having a unified AST that can be manipulated and validated is quite compelling. The MOL and SMILES/SMARTS are pretty sprawling and not well standardized formats so how would one go about unifying the two? Btw, I don't think I need to worry about any other representations here, so that's something. To go back to the RDKit example, imperfect as it may be, it defined its own mini-query language and seems to convert the input formats to this unified representations. If we assume that the advantages of the unified representation outweight the drawbacks of having to invent a "superformat" of the two, then we would need an AST (raw) format that can be produced from either MOL or SMILES/SMARTS and a validated molecule type (perhaps a standard and a general type again, the distniction could be useful). The question then arises how would one go about defining these representations under these circumstances?

# Response

You've hit the nail on the head. This is the central tension in designing cheminformatics toolkits: the pragmatic desire for a unified model versus the chemical reality of format-specific semantics. Your reference to `02-design-summary-2025-02-16.md` is perfect—the goal is to avoid coercing structures into a single representation that is "unphysically unified."

So, how do we get the best of both worlds? How do we build a unified pipeline without creating a "superformat" that is faithful to none?

The answer lies in making a clear distinction between the **raw Intermediate Representation (IR/AST)** and the final, **validated semantic `Molecule` model**. The IR doesn't have to be the "grand unified theory"; it just has to be a data structure capable of *holding all the information* from all supported formats, even the quirky, contradictory parts. The unification happens in the next step.

Here is a strategy for defining these representations.

### 1. The Raw IR: The Principle of Maximum Information

The primary rule for the raw IR is: **do not throw information away**. It should be a "superset" data structure that can losslessly represent the direct output of parsing either a MOL file or a SMILES string.

This means the IR structs will have a lot of `Option<T>` fields.

*   A **MOL parser** will populate fields like 3D coordinates. It will leave SMILES-specific fields (like explicit aromatic flags on atoms) as `None`. It will represent bonds as single/double/triple as specified in the file.
*   A **SMILES parser** will populate fields for atom-level aromaticity (from `c`, `n`, etc.) and bond-level aromaticity (from `:`). It will leave 3D coordinates as `None`.

Let's sketch this out conceptually:

```rust
// Conceptual sketch of a "superset" IR

struct RawAtom {
    // Core properties, common to all
    element: Element,
    formal_charge: i8,
    isotope: Option<u16>,
    
    // Properties often from MOL files
    position_3d: Option<(f64, f64, f64)>,
    
    // Properties often from SMILES/SMARTS
    is_explicitly_aromatic: bool, // was this a 'c' or 'C'?
    chiral_spec: Option<SmilesChirality>, // @ or @@
    
    // Query properties from SMARTS/MOL queries
    query_info: Option<AtomQuery>, // e.g., atom lists, wildcards
    rgroup_id: Option<u32>,
}

struct RawBond {
    atom_indices: (usize, usize),
    order: RawBondOrder, // Enum: Single, Double, Triple, Aromatic, Query
    stereo: Option<BondStereo>, // from '/' and '\' in SMILES
}

struct RawMoleculeIR {
    atoms: Vec<RawAtom>,
    bonds: Vec<RawBond>,
    // Could even have a field to track the source format
    source_format: InputFormat, // Enum: MOL, SMILES, SMARTS
}
```

With this approach, the parsing phase is purely syntactic. Its only job is to faithfully translate the source text into this rich, raw data structure. There is no loss of information and no attempt at unification yet.

### 2. The Validator: The Bridge to Semantic Unification

This is where the chemical intelligence lives. The validator takes a `RawMoleculeIR` and its job is to produce a `ValidatedMolecule`. This is the step where unification happens, and it's an *active, algorithmic process*, not a simple data mapping.

Consider how it would handle aromaticity for benzene:

*   **Input from SMILES (`c1ccccc1`):** The validator receives a `RawMoleculeIR` where `is_explicitly_aromatic` is `true` for all atoms and bond orders are `Aromatic`. It can trust this and pass it through.
*   **Input from a MOL file (Kekulé form):** The validator receives a `RawMoleculeIR` with alternating `Single` and `Double` bonds and `is_explicitly_aromatic` is `false`. Here, the validator must *run an algorithm*. It will detect the 6-membered ring with conjugated bonds, perceive it as aromatic, and *promote* the bonds in the `ValidatedMolecule` to `Aromatic`.

The validator is the component that embodies the "rules" of your toolkit. It decides what it means for a molecule to be "valid". This is where you can be opinionated and create the unified model. You could have different validators with different strictness levels.

### 3. The Final `Molecule` Models: Guaranteed Validity

Finally, the `MoleculeStandard` and `MoleculeGeneral` structs are constructed from the output of the validator.

*   `MoleculeStandard::from(validated_molecule)`: This constructor would be very strict. It would fail if the `validated_molecule` contained any query properties. Its existence guarantees that you are holding a representation of a single, specific, chemically plausible molecule.
*   `MoleculeGeneral::from(validated_molecule)`: This constructor would accept query properties, R-groups, etc. It provides a different set of guarantees.

### Summary of the Approach

1.  **Define a "Superset" Raw IR:** Design it by feature (chirality, aromaticity, coordinates, queries), not by format. Use `Option` liberally to accommodate features that only exist in some formats. The goal is a lossless syntactic representation.
2.  **Parse, Don't Interpret:** The parsers for MOL and SMILES simply populate the Raw IR. Their job ends there.
3.  **Unify in a Validator:** Create a separate validation stage that takes the Raw IR. This component contains the chemical logic. It perceives patterns (like aromaticity), checks rules (like valence), and produces a clean, semantically consistent output.
4.  **Construct Guaranteed-Valid Models:** The final user-facing `Molecule` types are built from the validator's output, ensuring that any such object in your system is in a known-good state.

This layered approach directly addresses your concerns. It allows you to build a powerful, unified manipulation and validation pipeline *without* losing the fidelity of the original formats during the initial, crucial parsing step. It avoids the "unphysically unified model" by making the unification an explicit, algorithmic step rather than an implicit side-effect of parsing.

# Prompt

Now, we don't have a SMILES/SMARTS parser at the moment. Your recommendation: Build a direct-to-structure SMILES/SMARTS now (after building the compliance suite for MOL) and, with that information, make the petit unified (opposite of grand unified) IR and validated molecules types? Or attempt to build the frist half of the IR and Molecule types based on MOL data alone (perhaps with some tentative SMILES/SMARTS info) and merge in the SMILES/SMARTS part as we go along?

# Response

That's the crucial strategic question. The trade-off is between making immediate architectural progress versus designing with complete information.

Here is my recommendation: **Build a direct-to-structure SMILES/SMARTS parser first, and *then* design the unified IR.**

Let's break down the reasoning for this "breadth-first" approach.

### The Case for Parsing SMILES First

1.  **Reduces Design Risk Drastically:** SMILES/SMARTS introduces concepts that have no direct equivalent in a standard MOL file.
    *   **Aromaticity:** The concept of atom-level aromaticity (`c`, `n`, `o`) is fundamental to SMILES but is an emergent property derived from a Kekulé structure in a MOL file.
    *   **Complex Queries:** SMARTS is a rich query language with features like recursive definitions (`$()`) and complex logical operators (`,`, `;`, `!`) that are far more expressive than MOL query features.
    *   **Implicit Hydrogens and Valence Model:** The SMILES valence model is implicit and context-dependent, which is very different from the explicit data in a MOL file.

    Attempting to design the "petit unified" IR based only on MOL data would force you to make assumptions about how to represent these SMILES-specific features. It is almost certain that these assumptions would be wrong in subtle but important ways, leading to a painful and costly refactoring of the IR—the very foundation of your new system—once you have a concrete SMILES implementation.

2.  **The "Throwaway" Work is Actually Requirements Gathering:** The initial SMILES-to-`Molecule` parser isn't truly throwaway work. It is the most effective form of requirements analysis. By building it, you will gain an intimate, implementation-level understanding of all the edge cases and semantic details. This deep knowledge is essential for designing a robust IR that is truly capable of representing both formats without loss of fidelity. The cost of refactoring this single parser's output is far lower than the cost of refactoring the core IR and all the components that depend on it.

3.  **Follows the "Concrete Before Abstract" Principle:** It's generally better to build two concrete implementations and then find the correct abstraction that unifies them, rather than creating an abstraction first and trying to force a new, complex implementation to fit into it.

### The Proposed Workflow

1.  **Solidify the Present:** Finish building the compliance suite for the current MOL parser. This gives you a solid, well-tested baseline for the existing functionality.
2.  **Explore the Future:** Implement the SMILES/SMARTS parser. The goal is feature-completeness, targeting the *current* `Molecule` and `MoleculeGeneral` structs. This is your deep dive into the second half of the problem domain.
3.  **Synthesize and Abstract:** With two working parsers, you now have all the information.
    *   Design the `RawMoleculeIR` that can losslessly represent the syntactic information from both MOL and SMILES/SMARTS.
    *   Design the `Validator` that contains the chemical logic to transform the `RawMoleculeIR` into a semantically consistent `ValidatedMolecule`.
    *   Refactor the `Molecule` types to be constructed from the `ValidatedMolecule`.
4.  **Refactor and Unify:**
    *   Update the MOL parser to target the new `RawMoleculeIR`.
    *   Update the SMILES parser to target the new `RawMoleculeIR`.

This approach is more methodical and significantly de-risks the most critical part of the process: the design of the core IR. It ensures that the foundation you build is based on a complete understanding of the structures it must support.