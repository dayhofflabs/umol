# Prompt

I'd like to work on creating a conformance suite for SMILES formats similar to those for MOL and SDF formats. It should not be limited to the OpenSMILES dialect but should serve a dual purpose:
1. Collect reference data from relevant sources:
a. SMILES inputs/outputs from drawing programs/viewers/editors/processing software
b. SMILES strings from major chemical, biochemical, pharma, tox, etc datases
c. SMILES strings from programmatic generators (rule-based and NN-based ones)
2. Establish an approximate classification of dialects (incl CXSMILES, ??)
We should also be looking for sources of SMARTS and extended SMARTS strings for future work.
Please make a plan for collecting a set of relevant sources first

# Response

Based on the existing MOL/SDF conformance suite structure and the materials already gathered, here's a plan for the SMILES conformance suite:

---

## SMILES Conformance Suite: Source Collection Plan

### 1. Cheminformatics Toolkits/Libraries (Producer + Consumer)

Sources that both generate and parse SMILES, useful for dialect fingerprinting:

| Source | Priority | SMILES Type | Notes |
|--------|----------|-------------|-------|
| RDKit | High | Canonical, Isomeric, CXSMILES | Already have test data in `materials/formats/smiles/rdkit-impl/test_data/` |
| CDK/Beam | High | OpenSMILES-adjacent | Have impl code in `materials/formats/smiles/cdk-impl/` |
| OpenBabel | High | Daylight-style | Have impl code in `materials/formats/smiles/openbabel-impl/` |
| ChemAxon | High | CXSMILES (origin) | Primary CXSMILES source |
| Indigo/EPAM | Medium | Own dialect | Already sourced for MOL suite |
| OEChem | Medium | Daylight-compatible | Commercial, limited access |
| LillyMol | Low | Own dialect | Already have in MOL suite |

### 2. Databases (Bulk Corpora)

| Source | Priority | Size | Access | SMILES Type |
|--------|----------|------|--------|-------------|
| PubChem | High | ~115M | FTP bulk download | Canonical + Isomeric |
| ChEMBL | High | ~2M | API + bulk | Canonical |
| ZINC | High | ~1B | Tranches available | Multiple formats |
| GDB-11/17 | High | 11M/166M | Already have GDB-11 in `materials/` | Exhaustive enumeration |
| ChEBI | Medium | ~130K | OBO/SDF available | Various |
| DrugBank | Medium | ~15K | Academic license | Various |
| COCONUT | Medium | ~400K | Open access | Natural products |
| HMDB | Medium | ~200K | Open access | Metabolites |
| KEGG | Medium | ~30K | Limited bulk | Pathways |
| Tox21 | Medium | ~8K | Open | Toxicology |
| LipidMaps | Low | ~50K | Open | Lipid-specific |
| CAS/SciFinder | Low | N/A | Commercial | Reference standard |

### 3. Drawing Programs/Viewers/Editors

| Source | Priority | Notes |
|--------|----------|-------|
| Ketcher (EPAM) | High | Open source, already in MOL suite |
| MarvinSketch | High | CXSMILES native |
| ChemDraw | Medium | Industry standard |
| ChemSketch | Medium | ACD/Labs |
| ChemDoodle | Low | Web-based |
| JSME | Low | Open source editor |

### 4. Programmatic Generators

**Rule-Based:**
| Source | Type | Notes |
|--------|------|-------|
| REINVENT | Generative | Randomized SMILES for training |
| GuacaMol | Benchmark | Standard benchmark corpus |
| MOSES | Benchmark | Molecular generation benchmark |
| Fragmentation-based | Enumeration | BRICS, RECAP outputs |

**NN-Based:**
| Source | Type | Notes |
|--------|------|-------|
| ChemBERTa/MolBERT | Tokenization | May produce unusual tokenizations |
| SELFIES→SMILES | Conversion | Guaranteed valid |
| MolGPT/ChemGPT | Generation | LLM-generated strings |
| JT-VAE | Latent space | Fragment-based outputs |
| MoleculeSTM | Multimodal | Text→SMILES outputs |
| DrugLLM outputs | Generation | Emerging source |

### 5. SMILES Dialect Classification

**Core Dialects:**
1. **Daylight SMILES** - Original spec (proprietary)
2. **OpenSMILES** - Community standardization attempt
3. **CXSMILES** - ChemAxon extended format (stereo groups, enhanced stereo, radicals, coordinates)

**Derived/Specialized:**
4. **Canonical SMILES** - Implementation-specific (RDKit ≠ CDK ≠ OEChem)
5. **Isomeric SMILES** - With stereochemistry
6. **Extended SMILES** - Various extensions (reactions, etc.)
7. **DeepSMILES** - ML-friendly modification
8. **SAFE SMILES** - Fragment-based encoding

**Classification Dimensions:**
- Aromaticity model (Daylight vs Hückel vs none)
- Stereochemistry encoding (slash vs @/@@ on vinylic)
- Extended features (radicals, atom classes, enhanced stereo)
- Ring closure syntax (%nn vs %nnn)
- Canonicalization algorithm

### 6. Documented Dialect Deviations from OpenSMILES

Major toolkits do NOT produce strict OpenSMILES. Known deviations:

**RDKit:**
- Accepts aromatic symbols like `te` (tellurium) which OpenSMILES disallows
- Accepts tetrahedral stereo on sulfoxides (3 substituents + lone pair) not defined in OpenSMILES
- `SmilesParserParams` allows disabling sanitization, accepting hypervalent atoms
- Accepts some ambiguous double-bond stereo encodings OpenSMILES rejects
- Supports subset of CXSMILES (coordinates, atom labels, radicals, etc.)

**CDK/Beam:**
- Supports CXSMILES extensions (atom values, labels, coordinates)
- Preserves hydrogen representation as-is rather than normalizing

**OpenBabel:**
- Claims OpenSMILES but adds extensions: radicals, "Universal SMILES"
- InChI-based canonicalization with tautomer/nitro normalization

**ChemAxon:**
- CXSMILES originator: trailing `|feature1,feature2,...|` syntax
- Atom labels, coordinates, R-groups, enhanced stereo, radicals, polymer/S-group info

**PubChem:**
- Input SMILES from heterogeneous sources (not controlled dialect)
- Standardization modifies ~44% of structures (de-aromatization, stereo normalization)

Sources: Depth-First blog (ChemCore vs RDKit comparison), ChemAxon CXSMILES docs,
RDKit Book, OpenBabel docs, PubChem standardization paper (J Cheminform 2018).


### 7. SMARTS and Extended SMARTS Sources

| Source | Type | Priority |
|--------|------|----------|
| Daylight SMARTS spec | Reference | High |
| RDKit patterns | Implementation | High |
| PAINS filters | Structural alerts | High |
| BRENK filters | Structural alerts | High |
| ChEMBL substructure patterns | Query patterns | Medium |
| OpenEye SMARTS | Extended patterns | Medium |
| MolPort filters | Commercial patterns | Low |

### 8. File Format

SMILES is single-line, so standard practice is newline-delimited `.smi` files.
Use a hybrid approach:

| Use Case | Format | Rationale |
|----------|--------|-----------|
| Bulk corpora | `.smi` (newline-delimited) | Standard format, efficient |
| Edge cases | `.smiles` (single input) | Descriptive filenames, snapshot-friendly |
| Invalid inputs | `.smiles` (single input) | Track specific failure modes |
| Toolkit extracts | Preserve original | If `.smi`, keep; if from code, use `.smiles` |

**Format details:**

```
# .smi files (bulk): SMILES<TAB>ID
CCO	CHEMBL545
c1ccccc1	CHEMBL277500

# .smiles files (edge cases): raw SMILES only, filename is identifier
# e.g., aromatic_tellurium.smiles contains: [te]1cccc1
```

### 9. Proposed Directory Structure

Following the MOL/SDF conformance suite pattern: raw inputs by source, then
automatic classification by parser configuration into `data/<config>/<source>`.

```
tests/smiles_parsing/
├── data_raw/                       # Original unmodified inputs by source
│   ├── rdkit/
│   │   └── rdkit_tests.smi         # extracted from test suite
│   ├── cdk/
│   ├── openbabel/
│   ├── indigo/
│   ├── chemaxon/
│   ├── pubchem/
│   │   └── pubchem_sample_10k.smi  # bulk download sample
│   ├── chembl/
│   ├── zinc/
│   ├── gdb/
│   ├── chebi/
│   ├── ketcher/
│   ├── marvin/
│   └── generators/
│       ├── reinvent/
│       └── selfies/
├── data/                           # Classified by parser config + source
│   ├── opensmiles_strict/          # Passes strict OpenSMILES parser
│   │   ├── rdkit/
│   │   │   └── rdkit_strict.smi    # subset from data_raw that passes
│   │   ├── pubchem/
│   │   └── edge/
│   │       ├── allene_stereo.smiles
│   │       └── ...
│   ├── opensmiles_lenient/         # Passes lenient parser
│   │   ├── rdkit/
│   │   └── ...
│   ├── cxsmiles/                   # Requires CXSMILES parser (|...| block)
│   │   ├── chemaxon/
│   │   ├── rdkit/
│   │   └── ...
│   └── invalid/                    # Fails all configurations
│       ├── rdkit/
│       ├── unclosed_ring.smiles
│       └── ...
└── snapshots/

tests/smarts_parsing/
├── data_raw/
│   ├── rdkit/
│   ├── daylight/
│   ├── pains/
│   └── brenk/
├── data/
│   ├── smarts_strict/
│   │   └── <source>/
│   ├── smarts_extended/
│   │   └── <source>/
│   └── invalid/
│       └── <source>/
└── snapshots/
```

**Parser configurations for classification:**

| Config | Description |
|--------|-------------|
| `opensmiles_strict` | Strict OpenSMILES spec: organic subset aromatics only, valid stereo templates |
| `opensmiles_lenient` | OpenSMILES + common extensions (extra aromatic atoms like `te`, relaxed stereo) |
| `cxsmiles` | OpenSMILES base + trailing `\|...\|` extension block |
| `invalid` | Fails all configurations |

**Classification dimensions to track per input:**
- Aromatic atom set used (organic subset vs extended)
- Stereo templates (TH/AL/SP/TB/OH)
- CXSMILES extensions present (coordinates, labels, enhanced stereo, radicals)
- Ring closure syntax (`%nn` only vs `%nnn`)
- Hydrogen handling (implicit vs explicit preserved)

### 10. Collection Priority Order

**Phase 1: Toolkit test suites (immediate)**
1. Extract SMILES from RDKit test suite (already have some)
2. Extract from CDK test suite
3. Extract from OpenBabel test suite
4. Extract from Indigo test suite

**Phase 2: Database samples (short-term)**
1. PubChem slice (10K diverse compounds)
2. ChEMBL slice (10K drug-like)
3. GDB-11 (already have)
4. ZINC slice (10K diverse)

**Phase 3: Edge cases and extensions (medium-term)**
1. CXSMILES from ChemAxon docs
2. Structural alerts (PAINS, BRENK)
3. Generator outputs

**Phase 4: Comprehensive coverage (long-term)**
1. Full database coverage
2. NN-generator outputs
3. Reaction SMILES

### 11. Implementation Status

**Completed (2026-01-23):**
- Classification binary (`classify_smiles_strings`) with sampling + individual file output
- Conformance test suite with insta snapshots (`tests/smiles_parsing/`)
- Molecule summary: sum formula, atom count, bond count
- Data structure: individual `.smiles` files (one SMILES per file)
- Sampling: 200 per source (random, seed 0)
- Sources extracted: rdkit, cdk, openbabel, indigo, zinc, gdb

**Current stats (opensmiles_strict parser):**
| Source | Total | Valid | Invalid | Valid % |
|--------|-------|-------|---------|---------|
| cdk | 200 | 199 | 1 | 99.5% |
| gdb | 200 | 200 | 0 | 100.0% |
| indigo | 200 | 200 | 0 | 100.0% |
| openbabel | 200 | 199 | 1 | 99.5% |
| rdkit | 200 | 195 | 5 | 97.5% |
| zinc | 200 | 200 | 0 | 100.0% |
| **Total** | **1200** | **1193** | **7** | **99.4%** |

**Next steps:**
- Add additional parser configs (opensmiles_lenient, cxsmiles)
- Expand molecule summary (rings, branches, topological params)
- Add more sources (pubchem, chembl, chebi)
- Extract CXSMILES examples from ChemAxon docs

---