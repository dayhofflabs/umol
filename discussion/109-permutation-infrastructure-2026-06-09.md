# Permutation and Assignment Infrastructure Discussion

**Status:** Informational

## Prompt

The umol-perm crate has a very basic set of perm data structures and algorithms, mostly targeted towards
stereochemistry treatment. One thing that is related to permutations but not stereochemistry (not directly anyway)
is Hungarian algorithm and friends, which are really relevant to RMSD calculations. What infrastructure would be needed
to impl this algorithm? Looking beyond that, which specific algorithms would be useful for cheminformatics / molecular
structure? Not looking at molecular spectroscopy so the complete molecular symmetry group (Bunker style, incl
permutations and point group symmetry) is not an immediate goal although it would be good to know which primitives are
needed. Are there libraries or frameworks to emulate here?

## Response

The important distinction is:

  > The Hungarian algorithm returns a permutation, but it is an optimization algorithm over bipartite assignments, not a permutation-group algorithm.

  It should probably not live in umol-perm.

  ## Current Boundaries

  The existing crates already suggest sensible ownership:

  - umol-perm: small-degree stereochemical permutation and coset algebra.
  - umol-graph-core: graph matching, graph automorphisms, canonical labeling, subgraph isomorphism.
  - umol-geometric: molecular coordinates, point-group operations, and geometry-dependent algorithms.

  Currently, umol-perm/src/permutation.rs:15 is:

  - Copy;
  - limited to degree 6;
  - designed for enumerating small stereochemical groups.

  umol-perm/src/group.rs:9 stores every group element explicitly. That is appropriate for tetrahedral through octahedral stereochemistry, but unsuitable
  for molecular automorphism groups or atom assignments.

  Meanwhile, umol-graph-core/src/algorithms/auto.rs:22 already owns graph automorphism results.

  ## Immediate Hungarian Use

  There is already a concrete place where linear assignment would improve correctness.

  umol-geometric/src/molecule.rs:875 currently maps each transformed atom independently to its nearest same-element atom:

  for i in atoms {
      perm[i] = nearest_compatible_atom(i);
  }

  This does not guarantee a bijection. Two transformed atoms can independently select the same target atom.

  Instead, for each element or atom-type equivalence block, construct

  [
  C_{ij}=\left|R x_i-x_j\right|^2
  ]

  and solve the minimum-cost assignment. The resulting mapping is guaranteed to be bijective.

  This is probably the first useful application to implement.

  ## Assignment Infrastructure

  I would add a small, chemistry-independent assignment module, probably in a new crate such as umol-optimize or umol-combinatorial.

  Putting it in umol-graph-core is defensible because linear assignment is minimum-weight bipartite matching, but it does not naturally operate on the
  existing undirected Graph. A standalone optimization crate gives the cleaner dependency direction.

  A minimal API could look like:

  pub struct Assignment {
      /// row_to_col[i] is the column assigned to row i.
      row_to_col: Vec<Option<usize>>,
      total_cost: f64,
  }

  pub enum AssignmentAlgorithm {
      Hungarian,
      JonkerVolgenant,
  }

  pub fn linear_sum_assignment(
      rows: usize,
      cols: usize,
      cost: impl Fn(usize, usize) -> f64,
      algorithm: AssignmentAlgorithm,
  ) -> Result<Assignment, AssignmentError>;

  Required semantics:

  - rectangular matrices;
  - complete and partial assignments;
  - forbidden pairs;
  - minimization, with maximization optional;
  - deterministic tie-breaking;
  - rejection of NaN;
  - behavior for ±∞;
  - optional dense matrix convenience API;
  - row-to-column and column-to-row maps;
  - total cost;
  - explicit verification that the result is injective.

  Do not return umol_perm::Permutation: rectangular and partial assignments are not permutations.

  ### Algorithm Choice

  Implement classic Hungarian first:

  - straightforward;
  - (O(n^3));
  - adequate for typical molecular atom counts;
  - easy to validate.

  Then consider Jonker-Volgenant. SciPy uses a modified Jonker-Volgenant implementation and supports rectangular assignment directly.

  For sparse compatibility graphs, a sparse minimum-cost bipartite matching or minimum-cost-flow implementation is preferable to a dense cost matrix.

  ## RMSD Is More Than Assignment

  For two point sets (X) and (Y), permutation-invariant RMSD requires solving

  [
  \min_{R,t,\pi}
  \sum_i \left|R x_i+t-y_{\pi(i)}\right|^2.
  ]

  There are three related but distinct problems:

  ### Fixed correspondence

  Given (\pi), solve for (R,t):

  - center the coordinates;
  - calculate the covariance matrix;
  - solve using Kabsch/SVD or a quaternion method.

  This belongs in umol-geometric.

  ### Fixed transform

  Given (R,t), solve for (\pi):

  - construct the squared-distance cost matrix;
  - run linear assignment.

  This uses the proposed optimization crate.

  ### Joint correspondence and transform

  Alternating:

  1. Kabsch;
  2. Hungarian;
  3. repeat;

  is an ICP-like local optimization. It does not guarantee the globally minimal RMSD.

  For chemically identical molecules, topology should usually restrict the allowed correspondences. Arbitrary same-element Hungarian assignment can map
  atoms in chemically invalid ways.

  The strongest workflow is:

  1. obtain chemically valid mappings from graph isomorphism or automorphisms;
  2. perform Kabsch for each candidate or search over the automorphism group;
  3. select minimum RMSD.

  RDKit’s symmetry-aware RMSD similarly works from atom mappings and warns that enumerating all matching atom orders can cause combinatorial explosion.

  ## Needed Automorphism Extension

  umol-graph-core currently exposes:

  - vertex orbits;
  - canonical labeling;
  - automorphism-group order.

  It does not expose automorphism generators.

  For scalable symmetry-corrected RMSD, generators are more valuable than enumerated automorphisms. Nauty can produce generators, and Automorphism
  should eventually expose something like:

  pub struct Automorphism {
      // existing fields
      generators: Vec<Vec<NodeId>>,
  }

  Then a large-degree permutation representation can consume those generators.

  That representation should be separate from the current degree-6 Permutation, for example:

  pub struct DynPermutation {
      image: Vec<usize>,
  }

  The likely long-term split is:

  - Permutation: compact degree-(\leq6) stereo permutation;
  - DynPermutation: atom-scale permutation;
  - PermutationGroup: explicitly enumerated small group;
  - GeneratedPermutationGroup: generator-based large group.

  ## Useful Algorithms for Cheminformatics

  ### Highest Priority

   Algorithm                          Main applications                                                           Likely owner
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━
   Linear assignment                  atom correspondence, transformed-atom matching, reaction mapping scoring    optimization crate
  ─────────────────────────────────  ──────────────────────────────────────────────────────────────────────────  ─────────────────────────
   Kabsch/Umeyama                     RMSD, conformer alignment                                                   umol-geometric
  ─────────────────────────────────  ──────────────────────────────────────────────────────────────────────────  ─────────────────────────
   Weighted Kabsch                    mass-weighted or feature-weighted alignment                                 umol-geometric
  ─────────────────────────────────  ──────────────────────────────────────────────────────────────────────────  ─────────────────────────
   Graph isomorphism mappings         chemically valid atom correspondence                                        umol-graph-core
  ─────────────────────────────────  ──────────────────────────────────────────────────────────────────────────  ─────────────────────────
   Automorphism generators            symmetry-corrected RMSD, symmetry pruning                                   umol-graph-core
  ─────────────────────────────────  ──────────────────────────────────────────────────────────────────────────  ─────────────────────────
   Minimum-cost bipartite matching    sparse assignment and constrained mapping                                   optimization/graph core

  ### Valuable Next Tier

  - Murty’s algorithm / k-best assignments
    Useful when the lowest-cost geometric assignment is chemically questionable and several alternatives should be evaluated.

  - Bottleneck assignment
    Minimizes the maximum matched distance rather than the sum. Useful for tolerance-based symmetry detection and geometry equivalence.

  - Quadratic assignment approximations
    Useful when pairwise distances or connectivity must be preserved. Exact quadratic assignment is NP-hard.

  - Branch-and-bound over correspondences
    Useful for globally minimizing symmetry-aware RMSD.

  - Optimal transport / Sinkhorn methods
    Useful for fuzzy, partial, or unequal-size molecular comparisons. Usually approximate and not a replacement for exact atom mapping.

  - Canonical image under a permutation group
    Useful for canonical stereochemical representations, symmetry pruning, and avoiding enumeration of equivalent assignments.

  - Orbit and stabilizer calculations
    Useful for equivalent atoms, stereotopicity, symmetry-reduced enumeration, and conformer search.

  - Double-coset enumeration
    [
    H\backslash G/K
    ]
    Useful when both the source and target possess internal symmetry. This is the natural group-theoretic structure behind eliminating duplicate
    mappings between symmetric objects.

  ## Group-Theory Infrastructure Beyond Stereo

  If atom-scale generated groups become useful, the central primitive is a base and strong generating set, usually built using Schreier-Sims.

  That provides efficient:

  - group-order calculation;
  - membership testing;
  - point and tuple orbits;
  - pointwise and setwise stabilizers;
  - coset representatives;
  - transversal construction;
  - random group elements;
  - group intersection;
  - canonical-image search.

  These primitives are much more important than enumerating all group elements.

  They would support symmetry-aware RMSD without attempting to implement a complete Bunker-style molecular symmetry group.

  ## What Bunker-Style Work Would Additionally Need

  The current permutation concept only records atom relabeling. Full molecular symmetry-group treatment requires operations carrying both:

  [
  (\pi, A),
  ]

  where:

  - (\pi) permutes identical nuclei;
  - (A\in O(3)) is a spatial operation;
  - possibly an additional inversion/parity marker;
  - multiplication composes both consistently.

  For non-rigid molecules, one also needs feasible permutation-inversion operations generated by internal motions, not merely graph automorphisms or
  equilibrium point-group operations.

  Useful primitives would therefore include:

  - dynamic atom permutations;
  - matrix representations;
  - permutation-matrix operation pairs;
  - generated groups and stabilizers;
  - subgroup and coset calculations;
  - representations and characters, only if spectroscopy later becomes relevant.

  The existing umol-msym atom permutations and matrix representations already contain part of this information, but full molecular symmetry groups would
  require a deliberate combined-operation abstraction.

  ## Libraries and Frameworks to Emulate

  - SciPy linear_sum_assignment: good assignment API semantics; rectangular problems; modified Jonker-Volgenant.
  - LEMON: good conceptual separation of graph structures, optimization maps, and algorithms; includes minimum-cost matching and flow.
  - nauty/Traces: already used correctly in umol-graph-core; canonical labeling and automorphism-group discovery.
  - GAP: the model to emulate for generated permutation groups, stabilizers, orbits, cosets, and Schreier-Sims.
  - RDKit rdMolAlign: useful cautionary model for molecular RMSD APIs, mappings, symmetry handling, and combinatorial explosion.

  ## Recommended Roadmap

  1. Add a standalone dense rectangular linear-assignment implementation.
  2. Replace greedy nearest-neighbor mapping in compute_atom_permutations.
  3. Add Kabsch/weighted Kabsch and explicit RMSD result types to umol-geometric.
  4. Expose automorphism generators from umol-graph-core.
  5. Add a dynamic generator-based permutation-group layer only when an actual consumer requires stabilizers or symmetry-aware search.
  6. Add symmetry-corrected RMSD using graph-valid mappings or automorphism generators.
  7. Add k-best assignment and sparse minimum-cost matching when reaction mapping or large structures require them.

  The main architectural recommendation is to keep umol-perm narrowly focused on finite, small-degree stereochemical arrangement algebra. Assignment
  optimization, graph automorphism discovery, and geometric alignment interact with permutations, but they are separate algorithmic domains.

  Sources: SciPy linear assignment (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html), nauty/Traces
  (https://pallini.di.uniroma1.it/), RDKit molecular alignment (https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html), Schreier-Sims overview
  (https://en.wikipedia.org/wiki/Schreier%E2%80%93Sims_algorithm).


## Reply

 Would it make sense to adjust the current simple structure to the more scalable generator-based infra?

## Response

• Yes, but not by replacing the current enumerated representation outright.

  The scalable and enumerated representations solve different problems:

  - Stereo coset numbering needs deterministic enumeration of tiny groups.
  - Molecular automorphism work needs generator-based groups whose elements cannot be enumerated.
  - Both need compatible permutation conventions and group-action APIs.

  ## Recommended Structure

  Separate the concepts explicitly:

  /// Atom-scale or arbitrary-degree permutation.
  pub struct Permutation {
      image: Box<[usize]>,
  }

  /// Small group represented by every element.
  pub struct EnumeratedGroup {
      degree: usize,
      elements: Vec<Permutation>,
  }

  /// Potentially large group represented by generators.
  pub struct GeneratedGroup {
      degree: usize,
      generators: Vec<Permutation>,
      stabilizer_chain: Option<StabilizerChain>,
  }

  The current PermutationGroup should become EnumeratedGroup. Its behavior is specifically defined by storing every element, and that is useful behavior
  rather than an implementation deficiency.

  ## Shared Interface

  Both group representations can implement a limited common trait:

  pub trait PermutationGroup {
      fn degree(&self) -> usize;
      fn generators(&self) -> &[Permutation];
      fn order(&self) -> GroupOrder;
      fn contains(&self, permutation: &Permutation) -> bool;
      fn orbit(&self, point: usize) -> Vec<usize>;
  }

  Do not put elements() in the common trait. Enumeration is not generally feasible.

  impl EnumeratedGroup {
      pub fn elements(&self) -> &[Permutation];
  }

  impl GeneratedGroup {
      pub fn try_enumerate(&self, limit: usize)
          -> Result<EnumeratedGroup, EnumerationLimit>;
  }

  Likewise, avoid returning plain usize from scalable order() because automorphism-group orders can exceed machine integers:

  pub enum GroupOrder {
      Exact(BigUint),
      Unknown,
  }

  ## Preserve Generators Now

  The current umol-perm/src/group.rs:17 computes the closure and discards the supplied generators. I would change that first:

  pub struct EnumeratedGroup {
      degree: usize,
      generators: Vec<Permutation>,
      elements: Vec<Permutation>,
  }

  That gives:

  - a natural conversion into GeneratedGroup;
  - documented generating sets for stereo geometry groups;
  - a uniform interface with future automorphism generators;
  - no behavioral disruption to current coset code.

  Constructors such as symmetric, alternating, and dihedral should also provide compact generating sets rather than treating every element as a
  generator.

  ## Permutation Representation

  The current degree-six Copy permutation is valuable for stereo-heavy code, so replacing it everywhere with heap allocation would be a regression.

  Two reasonable designs exist.

  ### Distinct Small and Dynamic Types

  pub struct SmallPermutation {
      image: [u8; 6],
      degree: u8,
  }

  pub struct Permutation {
      image: Box<[usize]>,
  }

  This is the clearest design. Add conversions between them.

  ### Small-Vector Storage

  pub struct Permutation {
      image: SmallVec<[u32; 8]>,
  }

  This gives one public type and avoids allocation for stereo permutations, but it loses Copy. Existing stereo code appears to benefit from Copy, so I
  would prefer distinct types.

  ## Coset Spaces

  The current umol-perm/src/coset.rs:30 fundamentally requires enumerated groups:

  - it assigns every coset a stable index;
  - stores representatives;
  - finds minimum-rank representatives by testing every group element.

  Keep it explicitly tied to EnumeratedGroup.

  A scalable generated-group equivalent would be a different abstraction:

  pub struct GroupActionSpace {
      group: GeneratedGroup,
  }

  impl GroupActionSpace {
      pub fn same_orbit<T: GroupObject>(&self, a: &T, b: &T) -> bool;
      pub fn canonical_image<T: GroupObject>(&self, object: &T) -> T;
  }

  Canonical-image and orbit algorithms are search problems over a stabilizer chain. They should not be hidden behind the existing CosetSpace API.

  ## Ownership With umol-graph-core

  umol-graph-core should continue to discover graph automorphisms. It should expose their generators as atom-index mappings.

  umol-perm can then consume those mappings:

  let generators = graph
      .automorphisms(colors, Nauty)
      .generators()
      .map(Permutation::from_image);

  let group = GeneratedGroup::new(atom_count, generators);

  This keeps the dependency direction clean:

  umol-graph-core discovers graph symmetry
            ↓
  consumer converts mappings
            ↓
  umol-perm performs generic group calculations

  umol-graph-core should not depend on umol-perm merely to return generators.

  ## When Schreier-Sims Is Worth Implementing

  A GeneratedGroup containing only generators is immediately useful for applying generators and simple orbit traversal. A full stabilizer-chain
  implementation becomes worthwhile when you need:

  - exact group order;
  - membership testing;
  - pointwise stabilizers;
  - setwise stabilizers;
  - coset traversal;
  - canonical images;
  - symmetry-pruned RMSD search.

  Do not implement Schreier-Sims merely to modernize the current stereo code. Add it when exposing graph automorphism generators creates a concrete
  consumer.

  ## Suggested Migration

  1. Rename the current type to EnumeratedGroup.
  2. Preserve compact generators alongside enumerated elements.
  3. Keep the current small Copy permutation for stereo code.
  4. Add an arbitrary-degree Permutation type with the same composition convention.
  5. Add GeneratedGroup initially with generator application and orbit calculation.
  6. Expose automorphism generators from umol-graph-core.
  7. Add stabilizer chains and Schreier-Sims when symmetry-corrected RMSD or canonical-image search requires them.

  So: adjust the architecture now to acknowledge generated groups, but retain the current enumerated implementation as the optimized small-group
  backend. Replacing it with generator-based machinery would make the stereo use case more complicated without providing a practical benefit.

## Reply

How does one go about reconciling the graph, permutation/group, and optimization frameworks (linear alg or constrained opt, for example) that are all
  used to solve related (or even the same type of problem)? service crate by algorithmic domain, bundled at the umol-graph / umol-geometric level ?

## Response

• Yes: organize low-level crates by algorithmic domain, then compose them in chemistry-aware service layers such as umol-graph and umol-geometric.

  The frameworks should not be forced behind one universal abstraction. A graph, permutation group, and cost matrix describe different aspects of a
  problem:

  - Graph: which mappings are chemically feasible?
  - Permutation group: which feasible mappings are equivalent by symmetry?
  - Optimization: which non-equivalent feasible mapping has the best score?
  - Geometry: how is a mapping scored or aligned?

  ## Architectural Principle

  Separate three layers:

  1. Problem representation
  2. Generic solvers
  3. Chemistry-aware orchestration

  chemical problem
        |
        v
  construct feasible mappings, costs, and symmetries
        |
        v
  invoke graph / group / optimization / linear algebra solvers
        |
        v
  interpret result as atom mapping, RMSD, stereoisomer, etc.

  The same chemical operation may use several solver domains without belonging to any one of them.

  ## Proposed Crate Roles

  ### umol-graph-core

  Generic discrete graph structures and algorithms:

  - graph isomorphism and automorphisms;
  - canonical labeling;
  - subgraph isomorphism;
  - maximum common subgraphs;
  - cardinality matching;
  - graph refinement;
  - connectivity and cycles.

  It should answer questions about topology, not coordinates or chemical atom semantics.

  Graph automorphisms should expose generator mappings as ordinary index vectors.

  ### umol-perm

  Generic group actions and small stereochemical arrangement algebra:

  - small and dynamic permutations;
  - enumerated groups;
  - generated groups;
  - orbits and stabilizers;
  - cosets and double cosets;
  - canonical image under a group action;
  - stereo-class coset numbering.

  It consumes permutations discovered elsewhere. It should not discover graph automorphisms itself.

  ### New umol-optimize

  Generic discrete optimization algorithms:

  - linear sum assignment;
  - sparse minimum-cost bipartite matching;
  - k-best assignment;
  - minimum-cost flow;
  - bottleneck assignment;
  - possibly generic branch-and-bound infrastructure.

  It should know nothing about atoms, graphs, or Cartesian coordinates.

  An assignment result should not be a Permutation, because assignments may be rectangular or partial:

  pub struct Assignment {
      row_to_col: Vec<Option<usize>>,
      total_cost: f64,
  }

  A complete square assignment can then explicitly convert into a dynamic permutation.

  ### umol-geometric-core

  Keep this limited to lightweight geometric value types needed across layers:

  - points;
  - vectors;
  - orientations;
  - transformations;
  - planes.

  It probably should not acquire heavyweight solvers.

  ### umol-geometric

  Geometry-aware algorithms and orchestration:

  - Kabsch and weighted Kabsch;
  - RMSD;
  - point-set alignment;
  - symmetry-corrected RMSD;
  - point-group operation to atom-permutation assignment;
  - geometrical correspondence scoring.

  It can depend on umol-optimize and potentially umol-perm.

  ### umol-graph

  Chemistry-aware graph services:

  - molecular graph construction;
  - chemically colored automorphisms;
  - chemically constrained isomorphisms;
  - atom mapping;
  - reaction mapping;
  - topology-aware mapping candidates.

  It can depend on umol-graph-core, umol-perm, and umol-optimize.

  ## Dependency Direction

  A reasonable dependency graph is:

  umol-graph-core     umol-perm     umol-optimize     umol-geometric-core
         ^                ^              ^                    ^
         |                |              |                    |
         +---------- umol-graph          |                    |
                          ^              |                    |
                          +------- higher chemical services   |
                                         |                    |
                                    umol-geometric ------------+

  More concretely:

  umol-graph
    -> umol-graph-core
    -> umol-perm
    -> umol-optimize

  umol-geometric
    -> umol-geometric-core
    -> umol-optimize
    -> umol-perm
    -> possibly umol-graph-core or umol-graph

  Avoid making low-level solver crates depend on chemistry-level crates.

  ## Shared Interchange Types

  The domains should communicate through a few simple types rather than shared framework traits.

  ### Index Mapping

  pub struct Mapping {
      source_to_target: Vec<Option<usize>>,
  }

  This is the common result of graph matching and assignment.

  It may represent:

  - a partial graph match;
  - a rectangular assignment;
  - a complete atom correspondence.

  Conversions can enforce stronger invariants:

  impl TryFrom<Mapping> for DynPermutation { ... }

  ### Compatibility Relation

  A compatibility relation says which assignments are allowed:

  pub trait Compatibility {
      fn allowed(&self, source: usize, target: usize) -> bool;
  }

  Implementations may come from:

  - atom types;
  - graph refinement colors;
  - element identities;
  - user constraints;
  - symmetry orbits.

  For a dense assignment solver, forbidden pairs become infinite costs. For sparse matching, they become absent edges.

  ### Cost Oracle

  pub trait Cost {
      type Value;
      fn cost(&self, source: usize, target: usize) -> Self::Value;
  }

  Costs may be:

  - squared Cartesian distances;
  - graph-edit penalties;
  - atom-feature differences;
  - reaction-center penalties.

  Do not require costs to be stored as a matrix. Solvers may materialize one when needed.

  ### Group Action

  Group actions should transform mappings or candidate states:

  pub trait Action<T> {
      fn act(&self, permutation: &DynPermutation, value: &T) -> T;
  }

  The group layer can then identify equivalent candidates without understanding their chemical meaning.

  ## How The Frameworks Compose

  ### Symmetry-Corrected RMSD

  The full objective is:

  [
  \min_{\pi\in M,;R,t}
  \sum_i \left|R x_i+t-y_{\pi(i)}\right|^2,
  ]

  where (M) is the set of chemically valid mappings.

  Responsibilities:

  1. umol-graph determines valid atom correspondences or automorphism generators.
  2. umol-perm removes symmetry-equivalent mappings or searches group actions.
  3. umol-geometric computes optimal (R,t) with Kabsch.
  4. umol-optimize may solve assignment subproblems.
  5. umol-geometric owns the complete RMSD workflow and result.

  The orchestration belongs in umol-geometric because RMSD is fundamentally a geometric service.

  ### Point-Group Operation to Atom Permutation

  For a spatial operation (R):

  1. umol-geometric constructs costs:
     [
     C_{ij}=|Rx_i-x_j|^2.
     ]

  2. Element or atom-type compatibility forbids invalid pairs.
  3. umol-optimize solves the assignment.
  4. The complete result converts to DynPermutation.
  5. umol-geometric validates the maximum assignment distance against a tolerance.

  This directly replaces the current non-bijective nearest-neighbor procedure.

  ### Reaction Atom Mapping

  1. umol-graph generates compatibility constraints from topology and atom types.
  2. umol-optimize minimizes bond-change or feature costs.
  3. umol-perm removes equivalent mappings caused by reactant/product automorphisms.
  4. umol-graph interprets and validates the result.

  This workflow belongs in umol-graph, not umol-optimize.

  ## Multiple Formulations Of One Problem

  A problem may admit several formulations:

   Problem                   Possible formulation
  ━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   Atom correspondence       graph isomorphism
  ────────────────────────  ─────────────────────────────────────
   Atom correspondence       linear assignment
  ────────────────────────  ─────────────────────────────────────
   Atom correspondence       minimum-cost flow
  ────────────────────────  ─────────────────────────────────────
   Atom correspondence       quadratic assignment
  ────────────────────────  ─────────────────────────────────────
   Symmetry-aware mapping    group orbit search
  ────────────────────────  ─────────────────────────────────────
   RMSD correspondence       alternating Kabsch and assignment
  ────────────────────────  ─────────────────────────────────────
   Reaction mapping          constrained graph edit optimization

  Do not hide this behind a single solve_mapping() implementation. Expose the formulation as a strategy:

  pub enum AtomMappingAlgorithm {
      GraphIsomorphism,
      LinearAssignment,
      BranchAndBound,
  }

  Or, for more complicated workflows:

  pub struct AlignmentOptions {
      pub correspondence: CorrespondenceAlgorithm,
      pub superposition: SuperpositionAlgorithm,
      pub symmetry: SymmetryHandling,
  }

  This makes approximation, completeness, and complexity visible to callers.

  ## Groups As Symmetry Reduction

  Permutation groups are usually not competing with graph or optimization solvers. They wrap them by reducing redundant work.

  If (G_s) and (G_t) are source and target automorphism groups, equivalent mappings are related by:

  [
  \pi' = g_t,\pi,g_s^{-1}.
  ]

  The unique correspondence classes therefore have a double-coset-like structure:

  [
  G_t \backslash M / G_s.
  ]

  That observation supports:

  - pruning equivalent graph-isomorphism matches;
  - avoiding duplicate assignment solutions;
  - symmetry-aware k-best mappings;
  - reducing RMSD candidate enumeration.

  This is where a generated-group umol-perm layer becomes valuable.

  ## Avoid A Generic “Solver Framework”

  A universal solver trait tends to erase important distinctions:

  trait Solver<Problem> {
      type Solution;
      fn solve(&self, problem: &Problem) -> Self::Solution;
  }

  It provides little value because each domain has different:

  - feasibility conditions;
  - optimality guarantees;
  - error modes;
  - result multiplicity;
  - numerical tolerances;
  - completeness guarantees.

  Prefer narrow domain APIs and compose them in service code.

  What is worth standardizing across algorithms is metadata:

  pub enum Optimality {
      Proven,
      LocallyOptimal,
      Heuristic,
  }

  pub struct SearchLimits {
      pub max_solutions: Option<usize>,
      pub max_nodes: Option<usize>,
      pub time_limit: Option<Duration>,
  }

  ## Recommended Direction

  Use algorithmic-domain crates for reusable primitives and solvers, bundled by higher-level chemical service crates:

  - umol-graph-core: topology algorithms;
  - umol-perm: group-action algorithms;
  - umol-optimize: assignment and discrete optimization;
  - umol-geometric: geometric solvers and geometric workflows;
  - umol-graph: chemistry-aware graph workflows.

  The higher-level crate that owns the scientific question should orchestrate the lower-level solvers. Graphs define feasibility, optimization selects
  candidates, groups remove equivalent candidates, and geometry supplies spatial objectives.