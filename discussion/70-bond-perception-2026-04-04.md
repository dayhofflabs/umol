# Bond perception algorithm

# Query

Here is an optimization problem. Let's assume that we're given N integer
variables {b_i} with l_i <= b_i <= u_i and corresponding real valued functions
f_i. Maximize the product of f_i(b_i), i = 1, .., N subject to a set of M sum
rules in b_i : sum_k b_{jk} = v_j, where v_j are constant. Which approaches
exist to solve this problem?

# Response

Good problem. Since the objective is a product of separable functions and the
constraints are linear equality constraints on integers, you're in a nicely
structured space.

**Take the log.** The product $\prod_i f_i(b_i)$ becomes $\sum_i \log
f_i(b_i)$, so you have separable concave (assuming $f_i > 0$) maximization over
integers with linear equality constraints. This is the key structural
observation.

**Approaches, roughly ordered by practicality:**

**1. Integer programming with tabulated objectives.** Each $b_i$ lives in a
finite range $[l_i, u_i]$. Introduce binary variables $x_{i,k}$ for "$b_i = k$"
and set the objective coefficient to $\log f_i(k)$. The sum constraints become
linear in the $x_{i,k}$. Standard ILP solvers (HiGHS, SCIP, Gurobi) handle
this. Works well for moderate $N$ and ranges.

**2. Lagrangian relaxation / dual decomposition.** Dualize the $M$ sum
constraints with multipliers $\lambda_j$. The relaxed problem decomposes into
$N$ independent single-variable problems: $\max_{b_i} \log f_i(b_i) - \sum_j
\lambda_j [\text{coeff of } b_i \text{ in constraint } j] \cdot b_i$. Each
subproblem is trivial (enumerate the finite range). Update $\lambda$ by
subgradient or bundle methods. Often gives exact solutions when the relaxation
is tight (which it frequently is for this structure).

**3. Dynamic programming on the constraint structure.** If the constraints have
tree/chain structure (each $b_i$ appears in few constraints), you can do DP
over the constraint graph. Complexity scales with the product of the constraint
RHS values, so only practical when $v_j$ are small.

**4. LP relaxation + rounding.** Solve the binary formulation from (1) as an
LP. The constraint matrix is often totally unimodular or near-TU (depending on
the overlap structure of your sum rules), in which case the LP solution is
already integral.

**5. Belief propagation / message passing.** The factor graph has variable
nodes $b_i$ and factor nodes for each constraint. Max-product BP on discrete
variables. Exact on trees, good heuristic otherwise. Essentially the
probabilistic analog of (3).

**What matters for choosing:** How large are $N$, $M$, the ranges $u_i - l_i$,
and the overlap pattern (does each variable appear in one constraint or many)?
If the constraints partition the variables into independent groups, this
trivially decomposes. If the constraint matrix is sparse/structured, (2) or (3)
will dominate. For a general dense case with moderate size, (1) is the most
pragmatic.

The valence-constraint / integer-partition flavor of this (given your `umol`
work) suggests $M$ is small, ranges are tight, and variables cluster—Lagrangian
relaxation is probably your best bet there.


# Reply

If N is a few hundred and the ranges are at most 4, usually smaller?

# Response

With $N \sim$ hundreds and $|[l_i, u_i]| \leq 4$, the binary reformulation
(approach 1) is tiny. You get at most $\sim 4N$ binary variables (~few hundred)
and $M + N$ constraints. Any ILP solver will eat this in milliseconds.

But honestly, with ranges that small you might not even need a solver.
Lagrangian relaxation will converge in very few iterations since each
subproblem is just picking the best of ≤4 values, and the dual is
low-dimensional ($M$ multipliers). A simple subgradient loop in pure Rust would
be ~50 lines and zero dependencies.

If $M$ is also small (which I'd guess for valence constraints—maybe $M \leq
10$?), brute-force enumeration of the dual space on a coarse grid followed by
local refinement is also viable.

**My recommendation for `umol`:** Lagrangian dual with subgradient. No external
dependency, trivial to implement, and the structure (small ranges, few
constraints, separable objective) is exactly where it shines. Fall back to
HiGHS via the `highs` crate only if you hit cases where the duality gap doesn't
close.

# Reply

Ok you got me. It's just a bond order perception problem with valence constraints. ;-)

# Response

Right, so the variables are bond orders $b_i \in \{1,2,3\}$ (or $\{0,1,2,3\}$
in some contexts), the $f_i$ score how plausible each bond order is (from
distances, electronegativity, whatever), and the sum constraints are just
$\sum_{i \in \text{bonds}(a)} b_i = v_a$ for each atom's valence.

Clean problem. The constraint matrix is the incidence matrix of the molecular
graph, which is totally unimodular. So the LP relaxation of the binary
formulation is guaranteed integral—you don't even need ILP, just an LP solve
gives you the exact integer optimum.

That makes HiGHS via LP (not ILP) the pragmatic choice: one call, exact
solution, already in your dependency tree. The Lagrangian approach also works
but is now overkill given TU gives you integrality for free.

The only wrinkle is aromaticity—if you want fractional/resonance bond orders or
need to handle kekulization as a separate step. But for assigning integer bond
orders from 3D coordinates, LP on the incidence matrix is about as clean as it
gets.

# Reply 

Is unimodularity a function of the graph being connected? Was curious if I
could allow 0 bond orders.


# Response

TU of the incidence matrix holds for any bipartite graph, not just connected
ones. The classical result is:

- Incidence matrix of a bipartite graph → TU
- Incidence matrix of a non-bipartite graph → not TU (the odd cycle gives a ±2
  determinant)

Molecular graphs have odd cycles all the time (cyclopentadienyl, tropylium, any
odd ring), so strictly speaking you lose TU. In practice the LP relaxation is
still often integral due to the tight variable bounds, but the guarantee is
gone.

Allowing $b_i = 0$ doesn't affect TU—it's just widening the box constraints.
The TU property is about the constraint matrix, not the bounds. So $b_i \in
\{0,1,2,3\}$ vs $\{1,2,3\}$ is irrelevant to whether you get integrality.

But if you allow 0, you're effectively allowing bond deletion, which changes
the problem semantically—you'd want the $f_i$ to strongly penalize $b_i = 0$
for bonds that clearly exist in the connectivity graph. Otherwise the solver
    might "solve" a hard valence constraint by just dropping a bond.

Given odd cycles break TU, your real options are:

1. **LP and hope** — often works, check integrality of solution, fall back to
ILP if not
2. **ILP directly** — still trivial at your scale, negligible overhead over LP
3. **Lagrangian** — if you want zero dependencies

At $N \sim$ hundreds with range ≤ 4, just go ILP via HiGHS and stop thinking
about it.

# Reply

Is a solution guaranteed to exist? I don't think it can be unique.

# Response

**Existence:** Not guaranteed. The valence constraints may be infeasible for a
given connectivity—think a carbon with only one neighbor and no option for a
triple bond. Or more subtly, a global inconsistency where satisfying one atom's
valence forces a violation elsewhere. You need to handle the infeasible case
gracefully (which HiGHS will report).

**Uniqueness:** Definitely not in general. Kekulé structures are exactly the
non-unique case—benzene has two equivalent bond order assignments. Any molecule
with resonance structures gives multiple optima. Even without resonance, if the
$f_i$ are flat (e.g., bond distance is ambiguous between single and double),
you get degeneracy.

So you're really looking at an enumeration problem if you want all solutions,
which is much harder. Options:

- If you just want **one** optimal assignment: ILP, done.
- If you want to **count** or **enumerate** optimal solutions: that's #P-hard
  in general (counting perfect matchings / Kekulé structures is the classic
  example). For small molecules it's still tractable via backtracking with
  constraint propagation.
- If you want a **distribution** over near-optimal solutions: this is where the
  Lagrangian dual or belief propagation framing becomes more natural—you get
  marginal probabilities $P(b_i = k)$ essentially for free, which may be more
  useful than enumeration anyway.

For `umol`, I'd guess you want one good assignment plus maybe a flag indicating
degeneracy (multiple optima exist). You can detect degeneracy cheaply: solve
the LP relaxation and check if any $b_i$ is fractional or if the optimal
objective has multiplicity (perturb slightly and re-solve).

# Reply

Both are interesting but probabilities are cool.

# Response

Then you want the Boltzmann distribution over feasible assignments:

$$P(\mathbf{b}) \propto \prod_i f_i(b_i) = \exp\left(\sum_i \log
f_i(b_i)\right)$$

subject to the valence constraints. This is exactly a discrete graphical model
/ factor graph where the atom valence constraints are hard factors and the
$f_i$ are local potentials on edges.

**Belief propagation** is the natural fit. On the molecular graph:

- Messages on each bond edge, domain size ≤ 4, so messages are just 4-vectors
- Atom factor nodes enforce the valence sum constraint over incident bonds
- Converged beliefs give you $P(b_i = k)$ directly

For trees it's exact. Molecular graphs have cycles, so you need loopy BP—but
with tight constraints and small domains it typically converges fine. Junction
tree is exact but overkill unless your molecules have very high treewidth (they
rarely do for drug-like molecules, treewidth ≤ 3–4 typically).

**Implementation in Rust:** This is maybe 200 lines. The core loop is:

1. Each bond sends a message to each endpoint atom (a distribution over
$\{0,1,2,3\}$)
2. Each atom collects messages from all incident bonds, computes the
constrained marginals (enumerate all tuples summing to $v_a$, weight by
incoming messages), sends updated messages back
3. Repeat until convergence
4. Read off marginals

Step 2 is the only nontrivial part—but with degree ≤ 4 and bond orders ≤ 3, the
enumeration at each atom node is tiny. For a hexavalent atom with 6 neighbors
and range 4, that's $4^6 = 4096$ tuples to filter. In practice much less with
pruning.

Zero external dependencies, gives you the full posterior over bond orders, and
naturally handles degeneracy (Kekulé structures show up as $P(b_i=1) \approx
P(b_i=2) \approx 0.5$ on aromatic bonds). Feels like the right fit for `umol`.
