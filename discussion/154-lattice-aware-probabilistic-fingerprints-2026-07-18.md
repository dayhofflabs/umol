# Lattice-aware fingerprints for partially specified reactions

Status: **Exploratory design note**
Date: 2026-07-18
Relates: [113](113-ast-canonical-equality-and-lattice-2026-06-14.md), [126](126-fingerprint-facility-2026-06-22.md), [151](151-python-molecule-workflows-2026-07-13.md), [152](152-basic-molecule-wildcards-2026-07-18.md)

## Motivation: enzyme–reaction contrastive learning

The immediate use case is not substructure screening. Ordinary RDKit Morgan
fingerprints, combined with DRFP, have worked well as reaction inputs for
enzyme–reaction contrastive model training. The difficulty is that enzymatic
reaction records commonly contain `*`, while the ordinary circular
fingerprints require concrete atoms. Replacing each wildcard with `CH3` has
been a workable preprocessing rule, but it is an arbitrary imputation.

Methyl replacement asserts considerably more than the source reaction does:

- the attachment atom is carbon;
- it is terminal and connected by a particular bond;
- it carries the hydrogens and local valence of a methyl group;
- the unknown structure ends after that atom.

In many reaction records, `*` instead marks an unspecified continuation of the
molecular structure: an R group or other context outside the described
transformation. It is therefore not merely an atom whose element is missing.
Sampling another element in place of carbon would improve the imputation model
without addressing the unknown continuation.

The practical question is:

> Can a circular reaction representation remain approximately invariant when
> peripheral concrete substituents are replaced by lattice-valued wildcards,
> while retaining sensitivity to the transformed reaction center?

This is an ML representation problem. Exact recovery of the AST partial order
would be useful but is not required. The desired fingerprint may be lossy,
probabilistic, and frequency-weighted as long as it improves the contrastive
task and has a stable, reproducible contract.

## Target properties

A useful construction should provide the following behavior.

1. **Ground compatibility.** A fully specified reaction retains the ordinary
   Morgan information that is already useful to the model.
2. **Masking stability.** Replacing a peripheral substituent by `*` changes the
   representation less than replacing the reaction center or changing the
   transformation.
3. **Concrete overlap.** A wildcard environment shares features with compatible
   concrete environments. It must not be merely a new atom type whose circular
   hashes are unrelated to all concrete hashes.
4. **Reaction-side coherence.** When corresponding wildcards denote the same
   unspecified group on both sides, the representation preserves that coupling.
   Unchanged unknown context should cancel or remain neutral; a changed
   attachment environment should remain visible.
5. **Fixed-vector use.** Binary, count, or floating output remains suitable for
   the existing contrastive model input and reaction combinators.
6. **Reproducibility.** Projection rules, priors, sampling, and hashes are named
   and versioned when they affect output.

An exact no-false-negative substructure screen is a related application, not a
requirement of this use case. Exact matching remains available separately.

## What the wildcard denotes

The first data audit must distinguish several meanings that the same surface
`*` may currently carry:

- a single atom with an undetermined element but otherwise known local graph;
- an attachment point standing for an unspecified rooted substituent;
- a shared variable whose completion is the same on both reaction sides;
- unrelated unknowns that happen to use the same surface spelling;
- an unknown atom or group participating directly in the transformation.

These cases should not be collapsed accidentally. The AST lattice describes
undetermined attribute values exactly. An R-group continuation additionally
requires reaction-level correspondence and a convention about the omitted
subgraph. A fingerprint scheme can support both, but its input semantics must
say which interpretation applies.

## Exact order and why an ordinary hash cannot preserve it

Let `D` be a finite ground domain and let `A(x) ⊆ D` be the set of ground values
admitted by an AST attribute `x`. Use the umol lattice orientation:

```text
x ≤ y  iff  A(x) ⊆ A(y)
```

A literal is specific, `Undetermined` admits all of `D`, and a target matches a
pattern when `target ≤ pattern`.

Suppose an encoding `h` and comparison `f` satisfy

```text
f(h(x), h(y))  iff  x ≤ y.
```

Then `h` must be an order embedding: it must both preserve and reflect the
order. It must therefore be injective, because equality can be recovered by
checking the order in both directions. The powerset lattice of an `n`-element
domain contains `2^n` values, so an exact fixed-width encoding needs at least
`n` bits. An ordinary lossy `u64` hash cannot provide this contract for a larger
domain.

The direct exact encoding is the characteristic vector of the admitted set:

```text
h(x) = bits(A(x))
x ≤ y  iff  h(x) & !h(y) == 0
```

This is a lossless lattice representation rather than a hash. It preserves
intersection and union as bitwise AND and OR.

The complementary information encoding is useful for substructure screening:

```text
i(x) = bits(D \ A(x))
target ≤ pattern  implies  i(pattern) ⊆ i(target).
```

The general pattern contains less asserted information than the concrete
target. In this orientation, a wildcard naturally contributes no
element-specific information.

Base attributes do not require approximation merely because their declared
domain is large. Element alternatives fit in 128 bits. Numeric sets and ranges
can remain sparse or symbolic. The difficult growth occurs when several
attributes and several sites are combined into refinement neighborhoods.

## Exact lattice-aware refinement

Ordinary WL replaces a rooted neighborhood with an opaque color. A
lattice-aware variant can instead retain its recursive structure:

```text
c0(v) = atom_attributes(v)
c(r+1)(v) = (
    c(r)(v),
    multiset((bond_attributes(v, u), c(r)(u)) for each neighbor u)
)
```

The attribute leaves carry their existing lattice order. Recursive colors are
compared componentwise; unordered neighbor collections require a multiset
matching under that order. With a known graph correspondence this reduces to
pointwise comparison. Without one, the comparison includes a matching problem.

These colors can be canonicalized and interned. A conventional hash can index
them or reject obvious inequalities, but order comparison must still consult
the retained structure. Hashing the entire tuple through an avalanche hash
destroys the order.

An equivalent denotational construction lifts an ordinary ground refinement:

```text
C_r(x) = { c_r(g) | g is a ground refinement admitted by x }.
```

If `x ≤ y`, then `C_r(x) ⊆ C_r(y)` for any deterministic ground-color recipe.
The forward implication does not require a special hash. The converse requires
collision-free ground identifiers and a refinement description that separates
the relevant neighborhoods. Explicit `C_r` sets can grow exponentially as
uncertain choices combine, especially when constraints correlate choices at
different sites.

Even collision-free structured colors do not make 1-WL a complete molecular
order test. The ordinary 1-WL indistinguishability examples remain. A whole
molecule fingerprint also loses correspondence when it deduplicates colors
into a feature set. Exact molecular matching therefore remains the job of the
AST lattice and graph-matching facilities; WL can at most be an accelerator.

## Candidate representations

The main design choices are not equally ambitious.

| Construction | Useful property | Main limitation |
| --- | --- | --- |
| Replace `*` by `CH3` | Reuses ordinary Morgan unchanged | Invents a terminal methyl environment |
| Give `*` a dedicated atom invariant | Represents missingness explicitly | Wildcard colors remain unrelated to concrete colors |
| Omit every environment touched by `*` | Introduces no fictional chemistry | Removes potentially large neighborhoods around each wildcard |
| Emit circular features at several resolutions | Deterministic overlap between partial and concrete environments | Requires a small, frozen projection family |
| Average fingerprints over sampled completions | Direct probabilistic interpretation | Requires a completion model and additional computation |

The dedicated-token and omission baselines are worth retaining in an
experiment even though neither is a complete answer. They separate the value of
marking missingness from the value of lattice-aware overlap.

### Mask-aware multiresolution Morgan

The most direct first implementation is a Morgan-like featurizer that emits a
small family of projections for each circular environment instead of only its
fully typed identifier. Candidate levels are:

1. rooted topology and radius;
2. topology plus bond categories;
3. topology, bonds, and coarse atom classes;
4. every known atom attribute, with undetermined fields omitted;
5. the ordinary fully typed Morgan identifier when the environment is ground.

Each level has its own feature namespace. A concrete environment emits its
fully typed feature and the selected generalizations. A partial environment
emits only the levels justified by its known information. The overlap therefore
occurs in the fingerprint itself rather than depending entirely on the
downstream network to learn that a wildcard token resembles concrete atoms.

This must be a deliberately small projection family. Emitting every coarsening
of every attribute would recreate the full powerset and make feature volume
unmanageable. The useful projections should be frozen as part of the named
scheme and evaluated by ablation.

An R-group wildcard introduces a second boundary. Refinement must not walk into
a fictional methyl cap, but it should retain the known attachment bond and the
fact that structure continues beyond it. Environments entirely within the
known graph remain ordinary. Environments reaching the wildcard emit a boundary
projection. Generating every possible boundary projection for a concrete
molecule would again be combinatorial; the practical alternatives are a small
fixed family of cut projections or learned invariance from masked concrete
examples.

Count features are the natural substrate because projection levels can coexist
without losing multiplicity. They can subsequently be folded, concatenated by
level, or consumed as a floating vector. The ordinary exact Morgan channel
should remain available so ground reactions do not lose the representation that
already performs well.

### Expected Morgan fingerprints

Let `φ(G)` be the ordinary sparse or folded Morgan vector for a ground molecule
and let `X` denote a partially specified molecule. Given a distribution over
valid completions:

```text
Φ(X) = E[φ(G) | G is a completion of X].
```

For a bit fingerprint, each coordinate of `Φ(X)` is the probability that the
bit is set. For a count fingerprint, it is the expected count. The resulting
floating vector can be used directly by the contrastive model.

For a reaction, the expectation can feed the same representation shapes that
already work:

```text
Φ_difference(R) = E[φ(products) - φ(reactants)]
Φ_combined(R)   = E[concat(φ(reactants), φ(products))].
```

Corresponding wildcard variables must be sampled jointly across the two sides.
Otherwise the sampling process itself invents a change in an unchanged R group.
If a reaction-native feature depends jointly on both sides, compute that feature
on each coupled sample before averaging.

The completion distribution is the difficult part. If `*` denotes an R group,
the distribution is over rooted substituent environments, not merely elements.
A radius-limited fingerprint only needs the completion distribution out to its
radius, but even that distribution should respect attachment bond, atom state,
and chemical validity. Monte Carlo sampling is likely simpler than exact
distribution propagation for a first experiment.

### Masking as contrastive augmentation

Ground reactions provide supervision without requiring the native wildcard
records to have known completions:

1. identify peripheral bonds or substituents outside the reaction center;
2. replace selected rooted substituents by `*`;
3. treat the complete and masked reactions as positive views;
4. retain reactions with different transformed centers as negatives;
5. measure how performance changes with the mask's distance from the center and
   with the amount of removed structure.

This directly teaches the desired invariance. It also supplies an evaluation
corpus in which the hidden completion is known. Molecular contrastive learning
has already used atom masking and subgraph removal as positive augmentations;
the distinction here is that masking models a real reaction-data convention and
the encoder remains a circular fingerprint plus the existing contrastive model,
not necessarily a GNN.

The observed frequency imbalance can initially be learned from these examples
instead of being fully encoded in a hand-built prior. If multiresolution Morgan
plus masking augmentation is sufficient, the more expensive completion model
may not be needed.

## Secondary construction: frequency-weighted containment

The original lattice question still supplies a useful diagnostic or retrieval
score, but it is not the primary training representation. Let `μ` be a
probability measure over the ground domain and define

```text
S_μ(x, y) = μ(A(x) ∩ A(y)) / μ(A(x))
          = P[V ∈ A(y) | V ∈ A(x)].
```

This asks: after instantiating the more specific candidate `x`, how often does
`y` admit the result?

The score has the desired orientation:

- `x ≤ y` implies `S_μ(x, y) = 1`;
- if `μ` has full support, `S_μ(x, y) = 1` implies `x ≤ y`;
- the score is directional and is not generally symmetric;
- exclusions are penalized in proportion to their probability mass rather than
  their raw count.

For example, if `A(x) = {C, N, O}`, `A(y) = {C, N}`, and oxygen carries 0.02 of
the conditioned mass, then `S_μ(x, y) = 0.98`, rather than the uniform-set value
`2/3`.

Independently sampling one ground value from each operand and asking whether
the two literals are equal is not the same construction. Since distinct ground
values are incomparable, that score primarily measures agreement and
undesirably penalizes broad patterns. Directed conditional containment is the
closer relaxation of the lattice order.

This formula is not novel. It is a weighted version of the subsethood or
inclusion measures studied in fuzzy-set theory. Kosko related subsethood to
conditional probability, and later work axiomatized and generalized graded
subsethood. The relevance here is that the umol AST supplies a concrete,
canonical source for the admitted sets rather than beginning with arbitrary
fuzzy membership functions.

### The prior is part of the algorithm

The probability measure determines the meaning of the score and must therefore
be named, versioned, and reproducible. A useful initial form is

```text
μ = (1 - ε) μ_empirical + ε μ_base,
```

where `μ_empirical` captures observed chemistry and `μ_base` assigns every
valid value nonzero mass. Without the full-support component, score 1 means
containment only almost surely under the selected corpus; a semantically valid
but unobserved counterexample becomes invisible.

The useful distribution is joint rather than a product of independent field
marginals. Element, charge, valence, hydrogen count, spin, and bond environment
are strongly correlated. Independent sampling would create large numbers of
chemically invalid states. Reasonable levels of increasing complexity are:

1. empirical atom-state frequencies;
2. element-conditioned charge and valence frequencies;
3. bond-conditioned neighboring atom-state frequencies;
4. a model over complete valid rooted environments.

The first versions need not attempt the last level, but they should not describe
independent marginals as a chemical prior.

## Probabilistic WL and circular refinement

Let `P(. | x)` be a distribution over valid ground refinements admitted by an
uncertain rooted environment or molecule `x`. The most faithful extension of
the lattice score is

```text
S(x, y) = P[G is admitted by y | G ~ P(. | x)].
```

This definition does not depend on WL. It can be estimated by sampling ground
refinements of `x` and checking them against `y` with the ordinary AST matching
semantics. A WL or circular fingerprint can make the repeated checks cheaper.

There are three implementation families worth investigating.

### Sample ground refinements

Sample `k` valid ground assignments, run the existing frozen WL or circular
recipe on each, and retain a sketch of the resulting colors or feature sets.
The direct Monte Carlo estimator for containment has standard error

```text
sqrt(S(1 - S) / k).
```

The ground-color hash does not need lattice properties; it only needs to be
deterministic and have an acceptable collision rate. The approximation errors
then come from sampling, WL indistinguishability, and color/folding collisions,
which should be measured separately.

### Propagate distributions

Replace each scalar seed color with a distribution over ground colors and
propagate those distributions through the neighborhood. Exact propagation
requires products and convolutions over neighbor distributions and recreates
the combinatorial growth of explicit ground refinements. Factorization,
pruning, or sampling can control it, at the cost of modeling correlations only
approximately.

### Retain only entailed screening information

Emit only facts that every admitted instantiation entails. A wildcard emits no
element-specific fact. A composite circular feature whose identity depends on
an uncertain value must be generalized systematically or omitted. This gives a
weak but safe one-sided screen:

```text
exact match  implies  fingerprint screen accepts.
```

This is distinct from a calibrated probability score, but it is the easiest
way to support query structures without introducing false negatives. It also
matches the design direction of existing chemical substructure fingerprints.

## Sketches and hashes

Once approximation is allowed, the order semantics can live in the sampling or
feature construction rather than in the hash itself.

- Ordinary MinHash estimates unweighted Jaccard similarity, which is symmetric
  and is not itself the desired directed score.
- Consistent weighted sampling provides weighted MinHash sketches whose
  collision probability is weighted Jaccard. Given the two masses, an estimated
  intersection mass can be recovered and converted to directed containment.
- Asymmetric minwise hashing and containment-oriented MinHash variants directly
  address the bias that ordinary Jaccard/MinHash has when the operand sizes are
  different.
- Bloom-style OR sketches provide monotone containment screens with false
  positives. Frequent alternatives can receive dedicated bits while a rare
  tail shares buckets.
- Count sketches and random projections are suitable for distribution-vector
  similarity but do not naturally estimate directed containment.

A raw bit-containment ratio or Tanimoto value should not be called a probability
without calibration. If the public result is a number in `[0, 1]`, its defined
statistic, prior, sampling process, and uncertainty should be part of the named
scheme.

## Precedents

### Chemistry

The full reaction-learning construction does not appear to be a standard
molecular fingerprint, but several close pieces already exist.

1. **RDKit PatternFingerprint — query-safe information omission.** RDKit's
   pattern fingerprint is explicitly intended for substructure screening. When
   a matched feature contains a query atom or query bond, RDKit hashes only the
   fact that the generic pattern matched and omits the atom/bond typing detail.
   This is a direct precedent for retaining only information entailed by a
   query. It is a Boolean screen, not a learned reaction representation.

2. **RDKit LayeredFingerprint — several typing projections.** The layered
   fingerprint enumerates a subgraph and sets multiple bits using different atom
   and bond type definitions. This is close to the proposed multiresolution
   construction. It does not derive those projections from an AST lattice or
   address paired reaction-side wildcards.

3. **Molecular contrastive masking.** MolCLR uses atom masking and subgraph
   removal as molecular graph augmentations and maximizes agreement between the
   resulting views. This is a direct precedent for the proposed training
   experiment. The new point here is that masking corresponds to native
   wildcard-bearing reaction data and is applied around a circular-fingerprint
   encoder used for enzyme–reaction contrastive learning.

4. **Fuzzy pharmacophore fingerprints — graded feature assignment.** Bonachéra
   et al.'s fuzzy tricentric pharmacophore fingerprints map a molecular triplet
   gradually onto several related basis triplets instead of selecting one exact
   bin. This establishes fuzzy chemical fingerprint features and adapted
   similarities, but its fuzziness concerns geometric/pharmacophore resemblance,
   not solution-set inclusion of query attributes.

5. **CSI:FingerID — probabilistic structural-property vectors.** CSI:FingerID
   predicts posterior probabilities for molecular fingerprint properties from
   tandem mass spectra. Later scoring models dependencies among those properties
   with a Bayesian tree. This is a strong precedent for fingerprints whose
   components are probabilities and for treating feature dependence explicitly.
   Its uncertainty comes from structure inference from experimental data, not
   from a query AST lattice.

6. **Substructure fingerprinting generally.** Structural-key, layered, and
   pattern fingerprints rely on a one-sided implication: query bits must occur
   in every true target, while collisions and omitted detail may produce false
   positives. This remains relevant to the secondary safe-screen construction.

### Graph kernels and partially labeled graphs

1. The original WL subtree kernel turns refinement labels into feature vectors
   and compares their occurrence counts.
2. Propagation kernels operate on label distributions and explicitly support
   partially labeled and attributed graphs. This is close to propagating
   uncertain atom-color distributions, although the resulting kernel is a
   symmetric similarity rather than directed lattice containment.
3. Generalized and Wasserstein WL kernels replace strict equality of WL labels
   with graded distances between the corresponding unfolding trees or node
   feature distributions. They establish that fuzzy comparison can live inside
   the WL framework, but do not impose the umol partial order.

### Fuzzy sets and large-scale set search

1. Fuzzy subsethood supplies the mathematical form of directed graded
   containment and its relationship to conditional probability.
2. Consistent weighted sampling supplies compact sketches for weighted set
   similarity.
3. Asymmetric MinHash and containment MinHash supply retrieval schemes for
   directed containment rather than ordinary resemblance.

The ingredients are therefore established separately. The apparently uncommon
combination relevant to this use case is:

```text
partially specified, paired reaction structures
    + AST solution-set lattices and shared wildcard correspondence
    + multiresolution or expected circular features
    + masking invariance in enzyme–reaction contrastive training.
```

This conclusion comes from a focused precedent search, not a claim that no such
method exists anywhere in the chemical-fingerprint literature.

## Relationship to the current fingerprint work

The current ground-input contract for WL, ECFP, and Morgan remains coherent.
Those operations produce similarity fingerprints of concrete molecules, and
silently treating `*` as another element would violate the AST semantics.

Wildcard-aware reaction representation should be a separate experiment and,
if retained, a separate named featurizer or configuration. The existing methods
should not silently change behavior according to whether the input happens to
be ground.

The four contracts remain distinct:

- **ordinary similarity:** the current ground-molecule fingerprints.
- **partial-reaction ML representation:** multiresolution counts or expected
  floating features, optimized for masking stability and reaction-center
  sensitivity.
- **safe query screen:** Boolean, one-sided, no false negatives before the exact
  matcher.
- **graded containment:** directed score with a named prior and an error or
  confidence estimate.

The existing Morgan-plus-DRFP model is the empirical baseline. The new work
should measure whether DRFP remains useful rather than assuming in advance that
the wildcard-aware Morgan representation replaces it.

## Proposed design spike

The first investigation should answer the ML question before designing a broad
public API.

1. **Inventory native wildcard semantics.** Count wildcard-bearing reactions,
   their distance from mapped reaction centers, whether wildcard identities are
   shared across sides, and whether `*` denotes an atom or an omitted rooted
   substituent in each source.
2. **Build a controlled masking corpus.** From ground reactions, mask peripheral
   substituents at several distances from the center while retaining the known
   completion and the original enzyme association.
3. **Implement comparable representations.** Evaluate `CH3`, a dedicated
   wildcard token, omission of contaminated environments, and a small
   multiresolution Morgan projection family. Keep width, radius, reaction
   combination, and model architecture fixed.
4. **Measure representation behavior before training.** Compare complete–masked
   similarity, retrieval of the complete reaction among distractors, sensitivity
   to changed reaction centers, feature volume, collision rate, and runtime.
5. **Run the contrastive task.** Report the existing enzyme–reaction metrics for
   ground reactions, native wildcard reactions, and synthetically masked
   reactions separately. Confirm that wildcard gains do not degrade the ground
   subset.
6. **Ablate reaction combination.** Compare the new representation alone,
   ordinary Morgan plus DRFP, and the new representation plus DRFP.
7. **Escalate only if needed.** If deterministic multiresolution features and
   masking augmentation are insufficient, prototype expected Morgan vectors
   from coupled, radius-limited completion samples.

Open modeling questions remain: how to identify reaction-center distance
reliably, how shared wildcard variables are represented in source data, which
projection levels are useful, and whether native wildcard records have enough
ground analogues to validate completion assumptions.

## References

- RDKit, [Pattern Fingerprints](https://www.rdkit.org/docs/RDKit_Book.html#pattern-fingerprints), RDKit Book.
- Y. Wang et al., [Molecular Contrastive Learning of Representations via Graph Neural Networks](https://doi.org/10.1038/s42256-022-00447-x), *Nature Machine Intelligence* 4 (2022), 279–287.
- F. Bonachéra et al., [Fuzzy Tricentric Pharmacophore Fingerprints. 1. Topological Fuzzy Pharmacophore Triplets and Adapted Molecular Similarity Scoring Schemes](https://doi.org/10.1021/ci6002416), *J. Chem. Inf. Model.* 46 (2006), 2457–2477.
- K. Dührkop et al., [Searching molecular structure databases with tandem mass spectra using CSI:FingerID](https://doi.org/10.1073/pnas.1509788112), *PNAS* 112 (2015), 12580–12585.
- H. Shen et al., [Bayesian networks for mass spectrometric metabolite identification via molecular fingerprints](https://doi.org/10.1093/bioinformatics/bty245), *Bioinformatics* 34 (2018), i162–i171.
- N. Shervashidze et al., [Weisfeiler-Lehman Graph Kernels](https://www.jmlr.org/papers/v12/shervashidze11a.html), *JMLR* 12 (2011), 2539–2561.
- M. Neumann et al., [Propagation kernels: efficient graph kernels from propagated information](https://doi.org/10.1007/s10994-015-5517-9), *Machine Learning* 102 (2016), 209–245.
- M. Togninalli et al., [Wasserstein Weisfeiler-Lehman Graph Kernels](https://arxiv.org/abs/1906.01277), NeurIPS 2019.
- T. H. Schulz et al., [A generalized Weisfeiler-Lehman graph kernel](https://doi.org/10.1007/s10994-022-06131-w), *Machine Learning* 111 (2022), 2601–2629.
- B. Kosko, [Fuzziness vs. probability](https://doi.org/10.1080/03081079008935108), *International Journal of General Systems* 17 (1990), 211–240.
- V. R. Young, [Fuzzy subsethood](https://doi.org/10.1016/0165-0114(95)00045-3), *Fuzzy Sets and Systems* 77 (1996), 371–384.
- S. Ioffe, [Improved Consistent Sampling, Weighted Minhash and L1 Sketching](https://research.google/pubs/improved-consistent-sampling-weighted-minhash-and-l1-sketching/), ICDM 2010.
- A. Shrivastava and P. Li, [Asymmetric Minwise Hashing](https://arxiv.org/abs/1411.3787), WWW 2015.
- S. M. Koslicki and D. J. Zabeti, [Improving MinHash via the containment index with applications to metagenomic analysis](https://doi.org/10.1016/j.amc.2019.02.018), *Applied Mathematics and Computation* 354 (2019), 206–215.
