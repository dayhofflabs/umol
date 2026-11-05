# Representing chemical structures

## Purpose

This guide describes how chemical structures are expressed in the umol representation and how the
result is validated. The normative grammar is the
[DSL specification](../../umol-graph-ir/spec/umol-dsl-spec.md).

## Representation and validity

umol does not impose a chemistry. The representation admits any structure expressible as a
molecular graph with overlays, and validity is assessed by the
[three tiers](../development/data-types.md#integrity-checks-and-validation-tiers) of integrity,
invariants and conformance. Only conformance depends on a chemistry model, and the default chemistry
model is one choice among others: one set of rules may reject diborane or ferrocene while another
accepts them.

## Graph and overlays

The molecular graph consists of atoms and bonds. An atom carries its element, isotope, charge,
implicit hydrogen count, lone pairs, unpaired electrons and spin multiplicity. A localized bond
carries its order, charge, unpaired electrons and spin multiplicity. Relations that an ordinary
edge cannot express are overlays:

| Relation | Field | Example |
|---|---|---|
| Dative bond | `:dative-bonds` | BH3·NH3, Ni(CO)4 |
| Haptic coordination | `:dative-bonds`, several `:donors` | η5-cyclopentadienyl, order 3 |
| Multicenter bond | `:multicenter-bonds` | Bridges of B2H6 |
| Aromatic system | `:aromatic-systems` | Furan, tropylium |
| Noncovalent bond | `:noncovalent-bonds` | Hydrogen bond of the water dimer |
| Stereo configuration | `:stereo-atoms`, `:stereo-bonds`, or inline `#T` and `#C` | L-alanine, E-2-butene |

The molecule-map examples use neutral atom-charge defaults and leave other omitted atom fields
undetermined; see [Full and partial specification](#full-and-partial-specification).

Charge (`#c`), unpaired electrons (`#u`) and spin multiplicity (`#s`) can also be assigned to an
aromatic system or multicenter bond. They need not all be assigned to individual atoms. For example,
this description of O2− places its charge and unpaired electron on the O–O bond:

```clojure
{:atoms ["O" "O"] :bonds [[0 1 "2#c-#u#s2"]]}
```

A dative bond is directional, and its order is the number of donated electron pairs:

```clojure
{:atoms ["Ni#u0" "C#c-" "O#c+" "C#c-" "O#c+" "C#c-" "O#c+" "C#c-" "O#c+"]
 :bonds [[1 2 "3"] [3 4 "3"] [5 6 "3"] [7 8 "3"]]
 :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}
                {:donors [3] :acceptor 0 :attrs "1"}
                {:donors [5] :acceptor 0 :attrs "1"}
                {:donors [7] :acceptor 0 :attrs "1"}]}
```

A multicenter bond records the electron contribution of each participating atom, in the order of
its `:atoms`:

```clojure
{:atoms ["B#h" "B#h" "H" "H" "H" "H"]
 :bonds [[0 4 "1"] [1 5 "1"]]
 :multicenter-bonds [{:atoms [0 2 1] :attrs "[1,1,0]#c0#u0#s1"}
                     {:atoms [0 3 1] :attrs "[0,1,1]#c0#u0#s1"}]}
```

A noncovalent bond records the kind of interaction between two atoms and does not contribute to
the covalent graph:

```clojure
{:atoms ["O" "O" "H" "H" "H" "H"]
 :bonds [[0 2 :single] [0 3 :single] [1 4 :single] [1 5 :single]]
 :noncovalent-bonds [{:atoms [2 1] :attrs "Hbd"}]}
```

Aromaticity is not a bond order. Membership is marked with `#a` on atoms and bonds, in which case
resolution determines the aromatic system, or the system is given explicitly in
`:aromatic-systems`, whose `:attrs` list the π contribution of each atom:

```clojure
{:atoms ["O#a+" "C#a+" "C#a+" "C#a+" "C#a+"]
 :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]}
```

Write the aromatic and Kekulé descriptions as separate inputs. The aromatic description uses the
localized bonds for the σ skeleton and the aromatic system for the π electrons. The Kekulé
description uses localized single and double bonds, without aromatic marks or an aromatic overlay.
Combining the Kekulé double bonds with an aromatic overlay for the same π electrons counts them
twice and fails the electron-count invariants.

## Full and partial specification

A field omitted from an input takes its value from the reader's defaults and may remain
undetermined. Resolution imputes undetermined fields under a chemistry model. It returns an
Underdetermined result if the input is insufficient to fix the structure and a Contradictory result
if the input is inconsistent.

For example, this resolution-conformance input defaults omitted atom charges to zero while leaving
hydrogen counts, lone pairs and unpaired electrons undetermined:

```clojure
{:config-overrides {:atom {:charge :zero}}
 :input {:atoms ["C" "C" "O"]
         :bonds [[0 1 "1"] [1 2 "1"]]}}
```

`:input` is the molecule DSL; `:config-overrides` belongs to the conformance-suite wrapper. It applies
the override to `MoleculeDefaults::default()`. Explicit values, such as `#c-` or `#h0`, still take
precedence. `MoleculeDefaults::concrete()` is a different choice: it also fills omitted H counts,
lone pairs and unpaired electrons with zero.

An input is therefore either fully specified, in which case its meaning does not depend on the
model, or partially specified, in which case the result does. In a partial specification, each
stated value should be a deliberate constraint. Values introduced incidentally by a conversion or
an external syntax cannot be distinguished from deliberate constraints and should be removed.

Some constraints determine the species under the default chemistry model:

| Constraint | Effect |
|---|---|
| `C#h0` on a divalent carbon | Carbene rather than CH2 |
| `N#h#a+` on one ring nitrogen | Selects the NH tautomer of a purine |
| `P#h0` on a trivalent phosphorus | P(III) rather than P(V) |

Whether a value constrains the result can be determined by removing it and resolving again. The
outcome applies only to the chemistry model used.

## Alternative descriptions

A compound frequently admits more than one representation:

| Compound | Descriptions |
|---|---|
| Metal halide | Covalent M–Cl, or `Cl#c-` as a donor to the metal |
| Carbonyl ligand | `C#c-` and `O#c+` with a triple bond, donating to the metal, or M=C=O |
| S2N2 | Aromatic ring, or singlet diradical with a multicenter bond |

The choice is one of description and should be recorded where it is significant. Alternatives that
are themselves of interest are kept as separate inputs. Rejection of a description by a chemistry
model concerns that model and that description; it does not bear on representability.

## Sources

Author new structures directly in the DSL from the chosen chemical description. Do not write
SMILES as an intermediate and convert it to DSL. Importing an existing SMILES record is a separate
source-transcription task. In large molecules, atom ids in place of indices keep the bonds legible:

```clojure
{:atoms [[:n "N"] [:ca "C"] [:c "C"] [:o1 "O"] [:o2 "O"] [:cb "C"]]
 :bonds [[:n :ca :single] [:ca :c :single] [:c :o1 :double] [:c :o2 :single] [:ca :cb :single]]}
```

Standard SMILES cannot express dative, multicenter, haptic or noncovalent relations. umol's extended
SMILES dialect supports two-center dative arrows, but remains a narrower representation than the
DSL. Its hydrogen counts reflect its syntax: bracket atoms always carry a count and plain atoms
never do, and brackets are required for charge, stereo and aromatic NH whether or not the count is
significant. When importing existing SMILES data, review the bonding description and the values
introduced by conversion against the intended structure.

## Validation

An input that fails an integrity check is malformed. A structure intended as a molecule that
violates an invariant contains an error and is corrected. Failure of conformance is a result under
that chemistry model and is reported with the model, the atom and the reason; it is not a reason to
adapt the input until it is accepted.

Acceptance is no evidence of correctness, since a chemistry model admits many structures other than
the intended one. The identity of a structure is established from sources independent of umol, such
as database records, primary literature or crystal structures. Textbooks supply material, and
apparent errors in them are corrected from such sources.

## Conformance inputs

The resolution conformance suite exercises the default chemistry model, and its inputs are written
relative to that model. Hydrogen counts are generally left undetermined and are specified only
where they alter the MostSaturated result under the AtomTyping or Counts candidate source.
