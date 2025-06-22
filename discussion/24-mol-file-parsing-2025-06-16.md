Here is a detailed specification for implementing the adaptive MOL file parser. This document is designed to be a precise guide for implementation.

***

## Specification: Adaptive Cheminformatics File Parser

### 1. Guiding Principles

The primary objective is to create a robust, adaptive, and specification-compliant parser for MDL MOL, SDF, and RXN files (V2000 and V3000). The implementation shall adhere to the following principles:

*   **Fidelity to Specification:** The MDL CTFile Formats document (`ctfile.pdf`) [1] is the ground truth. Deviations for usability or to handle common, non-standard conventions must be explicitly documented.
*   **Explicitness and Transparency:** The parser must never silently ignore recognized chemical features. If a feature is encountered but not fully supported by the current parsing profile, a clear warning must be issued to the user. The internal representation of parsed features must be programmatically accessible and unambiguous.
*   **Internal Consistency:** A unified internal data model must be used. For example, an atom list should have a single, consistent representation regardless of whether it was parsed from a V2000 `M ALS` line [1], a V3000 `TYPE="[...]" ` string [1], or constructed via an API.
*   **User Control and Adaptability:** The parser will be adaptive by default, intelligently handling file contents. However, it will also provide optional, explicit parsing profiles to allow for performance optimization and stricter validation.

### 2. Core Data Structures

The parser will populate the following Rust data structures. The implementation must adhere to these definitions precisely.

```rust
// --- HASHMAP AND ELEMENT/BOND DEFINITIONS (ASSUMED) ---
// use std::collections::HashMap;
// pub struct Element; /*... */
// pub struct Bond; /*... */
// pub enum AtomStereoParity {... }
// pub enum NamedIsotope { D, T }
// pub struct SGroup; /*... */
// pub struct Conformer; /*... */
// pub use petgraph::graph::StableGraph;
// pub use petgraph::Undirected;
// --- END ASSUMED DEFINITIONS ---

/// Defines the specification for an atom list query feature.
#
pub struct AtomListSpec {
    /// The list of elements included in (or excluded from) the list.
    pub elements: Vec<Element>,
    /// A boolean flag indicating if this is a "NOT" list (exclusion list).
    /// `true` if it is a NOT list, `false` otherwise.
    pub is_not_list: bool,
}

/// Enumerates the types of non-elemental, "atom-like" entities
/// that can be represented in a MOL file.
#
pub enum AtomLike {
    /// An atom list, defined by the associated AtomListSpec.
    AtomList(AtomListSpec),
    /// An unspecified or generic query atom, storing the symbol ('A', 'Q', '*').
    Unspecified(char),
    /// An explicit lone pair, parsed as a pseudo-atom.
    LonePair,
    /// An R-group attachment point, storing its integer index (e.g., 1 for R1).
    RGroup(usize),
}

/// Represents the fundamental identity of an atom node in the graph.
/// It can be a standard element, a named isotope, or a query-like entity.
#
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    AtomLike(AtomLike),
}

/// Represents an atom node in the molecular graph. It is designed to be
/// flexible enough to represent both standard elements and query features.
#
pub struct Atom {
    /// The fundamental identity of the atom.
    pub symbol: AtomSymbol,
    /// The formal charge. `None` if not applicable or unspecified.
    pub charge: Option<i8>,
    /// The specific integer mass of an isotope. `None` for natural abundance.
    pub isotope_mass: Option<u32>,
    /// Stereochemical parity, if applicable.
    pub stereo_parity: Option<AtomStereoParity>,
    /// Explicitly defined hydrogen count for query purposes.
    pub hydrogen_count: Option<u8>,
    /// Explicitly defined valence, overriding standard rules.
    pub valence: Option<u8>,
    /// Atom-atom mapping number, primarily for reactions.
    pub atom_map_num: Option<u32>,
    /// Radical state (e.g., monovalent, divalent).
    pub radical: Option<u8>,
    /// A key-value map for storing additional, non-standard properties.
    pub properties: HashMap<String, String>,
}

/// Represents a single, contiguous chemical entity.
pub struct Molecule {
    /// The graph representation of the molecule, with `Atom` nodes and `Bond` edges.
    pub graph: StableGraph<Atom, Bond, Undirected, usize>,
    /// A vector of 3D conformers.
    pub conformers: Vec<Conformer>,
    /// A key-value map for molecule-level properties (e.g., title, comments).
    pub properties: HashMap<String, String>,
    /// A vector of S-Groups defined in the molecule.
    pub sgroups: Vec<SGroup>,
    /// A flag indicating if any query features were parsed.
    pub is_query: bool,
}

/// Represents a chemical reaction, consisting of reactants and products.
pub struct Reaction {
    /// A vector of molecules representing the reactants.
    pub reactants: Vec<Molecule>,
    /// A vector of molecules representing the products.
    pub products: Vec<Molecule>,
    /// A key-value map for reaction-level properties.
    pub properties: HashMap<String, String>,
}
```

### 3. Adaptive Parsing Strategy and Profiles

The parser shall implement an adaptive strategy controlled by optional parsing profiles.

#### 3.1. Parsing Profiles

The public API for parsing must accept an optional `ParsingProfile` enum:

```rust
pub enum ParsingProfile {
    /// (Default) Adaptive parser. Attempts to parse all supported features
    /// and infers the nature of the file (e.g., query, generic).
    Rich,
    /// Optimized for speed. Parses only basic structural information
    /// (atoms, bonds, coordinates, charge, isotopes). Ignores all query
    /// features, S-Groups, and other advanced properties.
    FastGeneric,
    /// Specialized parser for RXN files. Expects a reaction file format
    /// and produces a `Reaction` object.
    Reaction,
}
```

#### 3.2. Profile Behavior

*   **`Profile::Rich` (Default):**
    *   This is the main adaptive parser. It must attempt to parse all features defined in this specification.
    *   It will populate a `Molecule` object.
    *   The `molecule.is_query` flag must be set to `true` if any query-specific atom symbol ('L', 'A', 'Q', '*'), query bond type, or query-specific `M` property (e.g., `M SUB`, `M RBC`) is successfully parsed.[1]
    *   If it encounters an `\$RXN` header, it should issue a warning suggesting the use of `Profile::Reaction` for correct parsing.

*   **`Profile::FastGeneric`:**
    *   This profile prioritizes performance. It will populate a `Molecule` object.
    *   **It MUST parse:** Atom coordinates, standard element symbols, charge, and isotope information (`dd`, `M CHG`, `M ISO`, `M RAD`, `MASS=`).[1]
    *   **It MUST ignore:** All `AtomLike` symbols ('L', 'A', 'Q', 'LP', 'R#'), query bonds, S-Groups, and all other properties not listed above. Ignored features must not produce warnings in this mode to maintain performance.

*   **`Profile::Reaction`:**
    *   This profile is exclusively for parsing `rxnfile` formats. It must produce a `Reaction` object.
    *   It must first check for the `\$RXN` header.[1] If absent, it must return an error.
    *   It will parse the reactant and product counts, then loop through the `\$MOL` blocks, parsing each into a separate `Molecule` object using the `Rich` profile's logic.
    *   It must correctly parse reaction-specific properties like atom-atom mapping (`aamap`) and reacting center status (`ccc`/`RXCTR`) and use them to populate the final `Reaction` object.[1]

### 4. Detailed Parsing Instructions for Atom Symbols

The following instructions apply to the **`Rich`** and **`Reaction`** profiles.

#### 4.1. Atom Lists ('L')

*   **V2000:**
    1.  When an atom symbol `L` is found in the atom block, identify it as an atom list.[1]
    2.  Search the properties block for the corresponding `M ALS` line for that atom's index.[1]
    3.  Parse the `M ALS` line to extract the list of element symbols and the exclusion flag ('T' for NOT list, 'F' for normal list).[1]
*   **V3000:**
    1.  Parse the `TYPE` field in the atom block.[1]
    2.  If the type is of the form `"[...]"`, parse the comma-separated element symbols within the brackets. Set `is_not_list` to `false`.
    3.  If the type is of the form `"NOT [...]"`, parse the elements and set `is_not_list` to `true`.
*   **Internal Representation:**
    *   Create an `Atom` struct.
    *   Set `atom.symbol` to `AtomSymbol::AtomLike(AtomLike::AtomList(spec))`.
    *   The `spec` (`AtomListSpec`) must contain the parsed `elements` and the `is_not_list` boolean.

#### 4.2. Generic / Unspecified Atoms ('A', 'Q', '*')

*   **Parsing:** Recognize the symbols `A`, `Q`, and `*` in the atom symbol field (V2000 `aaa` or V3000 `TYPE`).[1]
*   **Internal Representation:**
    *   Create an `Atom` struct.
    *   Set `atom.symbol` to `AtomSymbol::AtomLike(AtomLike::Unspecified(char_code))`, where `char_code` is the character ('A', 'Q', or '*').
*   **Semantic Interpretation (for downstream logic):**
    *   `'A'`: Any atom except Hydrogen.[2, 1]
    *   `'Q'`: Any heteroatom (any atom except Carbon and Hydrogen).[2, 1]
    *   `'*'`: In V3000, this means "any atom but C or H".[1] In V2000, if encountered, it should be treated as "any atom including hydrogen" ('AH' in ChemAxon terms [2]). This distinction must be documented.

#### 4.3. Lone Pairs ('LP')

*   **Guiding Principle:** The library will prioritize calculated lone pairs for chemical consistency. Parsing of explicit 'LP' atoms is an optional feature for compatibility with legacy files.
*   **Default Behavior:** The parser shall **not** recognize the 'LP' symbol. Lone pairs should be calculated post-parsing based on valence, formal charge, and bonding topology.
*   **Optional Behavior:** The parser shall expose a boolean configuration flag, e.g., `parse_explicit_lone_pairs` (defaulting to `false`).
    *   If `true`, the parser will recognize the `LP` symbol in the V2000 atom block.[1]
    *   **Internal Representation (if parsed):**
        *   Create an `Atom` struct with `symbol: AtomSymbol::AtomLike(AtomLike::LonePair)`.
        *   This `Atom` node must be added to the graph and connected via a zero-order bond to its parent atom (inferred from file context, typically the preceding atom or via S-Group information).
    *   **Warning:** If `parse_explicit_lone_pairs` is `true`, the parser must issue a warning if the file contains explicit 'LP' atoms, informing the user that these may conflict with standard valence calculations.
    *   **V3000:** The `LP` symbol is not part of the V3000 atom type specification and must be ignored, even if `parse_explicit_lone_pairs` is true.[1]

#### 4.4. R-Groups ('R#')

*   **Parsing:** Recognize the `R#` symbol (e.g., R1, R2, etc.) in the atom symbol field (V2000 `aaa` or V3000 `TYPE`).[1]
*   **Internal Representation:**
    *   Create an `Atom` struct.
    *   Set `atom.symbol` to `AtomSymbol::AtomLike(AtomLike::RGroup(index))`, where `index` is the parsed integer from the R# label.
*   **R-Group Definitions:**
    *   The parser must also process R-group definition blocks (`M RGP`, `M LOG` in V2000; `BEGIN RGROUP` in V3000) and RGfiles.[1]
    *   These definitions (the list of allowed substituents for each R-group index) are not stored in the `Atom` itself. They must be stored at the `Molecule` level, for example in a `HashMap<usize, Vec<Molecule>>` within a dedicated R-group definition structure associated with the `Molecule`.

#### 4.5. Named Isotopes ('D', 'T')

*   **Formal Parsing:** The primary method is to parse an atom with symbol 'H' and then apply isotopic mass information from:
    *   V2000: The mass difference `dd` field or the `M ISO` property line.[1]
    *   V3000: The `MASS=val` keyword.[1]
*   **Convenience Parsing (Default):** The parser must also recognize `D` and `T` directly in the atom symbol field as a common convention.[3]
*   **Internal Representation:**
    *   Regardless of parsing method, the internal representation must be consistent.
    *   Create an `Atom` struct.
    *   Set `atom.symbol` to `AtomSymbol::NamedIsotope(NamedIsotope::D)` or `AtomSymbol::NamedIsotope(NamedIsotope::T)`.
    *   The `atom.isotope_mass` field **must** be set to `Some(2)` for Deuterium and `Some(3)` for Tritium.
    *   Conceptually, the underlying element is Hydrogen. API functions retrieving the element type should reflect this.

### 5. Error and Warning Strategy

A robust feedback mechanism is critical for user trust and data integrity.

*   **Errors:** The parser must raise a fatal, descriptive error and halt processing if it encounters a structurally malformed file that makes further parsing impossible (e.g., incorrect counts line that mismatches the number of atom lines).
*   **Warnings:** The parser must issue non-fatal warnings for any of the following conditions:
    *   An unrecognized `M` property line is found in the properties block. The warning should state the line was ignored.
    *   A feature is encountered that does not match the current parsing profile (e.g., a query feature in `Profile::FastGeneric` if that profile is ever modified to produce warnings, or an `\$RXN` header in `Profile::Rich`).
    *   An ambiguous or potentially conflicting feature is parsed (e.g., an explicit 'LP' atom is found when `parse_explicit_lone_pairs` is enabled).
*   **Warning Content:** All warnings must be specific, identifying:
    1.  The feature or symbol that triggered the warning.
    2.  Its location in the source file (line number).
    3.  How the parser handled it (e.g., "Ignored", "Interpreted as...").