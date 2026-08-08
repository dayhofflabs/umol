# Graph structure of internal coordinates

Status: **Informational**
Date: 2026-08-08
Relates: [190](190-enumeration-algorithm-candidates-2026-08-08.md),
[158](158-ring-model-and-enumeration-2026-07-22.md),
[072](072-msym-integration-2026-04-04.md)

## Scope

Records findings on the graph-theoretic structure of molecular internal
coordinates: which subgraphs carry primitive coordinate types and why the list
terminates, how Cremer–Pople puckering generalizes to polycyclic systems, the
in-principle limits on singularity-free coordinates, and the formal accounting
for redundant coordinate systems built from overlapping trees. Reference
material for future geometric work; nothing here is scheduled.

## Trees: characteristic coordinates end at order 4

The trees of order 2–4 carry the classical primitive types: P2 the stretch,
P3 the bend, P4 the proper torsion, K1,3 the out-of-plane. Define the
*characteristic coordinate* of a tree T as a rotation–translation invariant of
its atom positions that is new modulo the invariants inherited from proper
subtrees. Then:

- P4 inherits 3 stretches + 2 bends = 5 invariants against 3k − 6 = 6; the
  sixth is the torsion.
- K1,3 inherits 3 stretches + 3 bends, generically complete but degenerate at
  planarity; the out-of-plane coordinate is the regularizing replacement.
- P5 inherits 4 + 3 + 2 = 9 = 3·5 − 6: exactly closed, no new coordinate.
  K1,4 closes with one redundancy.

The termination is a theorem, not an accident of counting. By the first
fundamental theorem of invariant theory for SO(3) (Weyl), every smooth
rotation–translation invariant of a point configuration is a function of
pairwise inner products of difference vectors (at most three atoms per
generator) and signed 3×3 determinants (four atoms). Primitive generators
therefore involve at most four points; every tree of order ≥ 5 contributes
only relations — Cayley–Menger / Gram rank-≤ 3 conditions — never a new
coordinate type. Orders 2–4 are the generators of the invariant ring; order
≥ 5 supplies the syzygies.

A graph-Laplacian eigenvector characterization of tree coordinates does not
fit here — Laplacian modes are collective linear combinations, a basis choice
over the primitives, not primitive types. The Laplacian belongs to the cycle
story below.

## Cycles: Cremer–Pople as the cycle-graph Laplacian, and its generalization

Cremer–Pople puckering is the cycle-graph special case of a Laplacian
construction: the CP modes are the Fourier modes of C_N, which are exactly the
graph-Laplacian eigenvectors of the N-cycle. The natural polycyclic
generalization is to expand the transverse displacement field of a ring system
in the Laplacian eigenmodes of the ring-system subgraph:

- The Laplacian commutes with every graph automorphism, so degenerate
  eigenspaces carry irreps of the automorphism group — the construction is
  symmetry-adapted by construction. Automorphisms come from the existing
  canonical-labeling machinery; the displacement-SALC machinery of the msym
  integration applies directly.
- CP's polar (q, φ) form is singular at planarity only because it is polar;
  the two Cartesian amplitudes of each degenerate eigenmode pair are linear
  functionals of positions and remove that singularity outright.
- Which cycles constitute the ring system is the relevant-cycle / unique
  ring family question of the ring model.

Established practice for fused systems is weaker: per-ring CP with
shared-atom constraints (redundant, awkward at fusion bonds), Hill–Reilly
triangular-tessellation coordinates, or endocyclic-torsion formulations
(Haasnoot). No canonical polycyclic CP analog appears to be established; the
Laplacian construction above is a candidate. Precedent at the macromolecular
scale: Gaussian and anisotropic network models use exactly Kirchhoff-matrix
eigenmodes of a contact graph as collective coordinates.

## Singularities: what is possible in principle

The obstruction is geometric. Shape space R^{3N}/SE(3) is an orbifold:
at collinear configurations the SO(3) action loses freeness and the quotient
is genuinely singular (Littlejohn–Reinsch; Kendall shape theory). A torsion
lives on S¹, which no single-valued smooth chart covers. Consequently:

- **Singularity-free + minimal + global: impossible.** No (3N − 6)-coordinate
  chart covers shape space; the collinear stratum is singular in the quotient
  itself, so no choice of functions rescues a chart there. The
  degenerate-linear-bend pairs of Wilson–Decius–Cross are chart switching,
  not a resolution.
- **Singularity-free + redundant: possible.** This is an embedding rather
  than a chart. The invariant generators are polynomial in Cartesians and
  hence globally smooth: inner products, and — the non-distance members —
  signed volumes. The signed volume of a K1,3 quadruple is a
  singularity-free out-of-plane coordinate; (cos φ, sin φ) pairs embed
  torsions smoothly. Redundancy is forced by topology, not a defect of any
  particular definition.

## Overlapping trees: existing systems and the formal bookkeeper

- The standard redundant internal coordinate set (Pulay et al.; Peng–Ayala–
  Schlegel–Frisch) is precisely the union of characteristic coordinates over
  all embedded order-≤ 4 subtrees of the bond graph. Delocalized internals
  (Baker–Kessi–Delley) are its SVD orthogonalization; natural internals
  (Fogarasi–Pulay) are its symmetry-local combinations. A Laplacian eigenmode
  basis is a candidate canonical combination scheme over the same primitives.
- A Z-matrix is a chart built on one spanning tree. Optimizers that swap
  Z-matrices at degeneracies informally run an atlas of overlapping
  spanning-tree charts; a formalized ensemble-of-spanning-trees construction
  does not appear in the literature.
- The independence and redundancy accounting for any such union is the 3D
  generic rigidity matroid: a coordinate set is independent exactly when its
  Wilson B-matrix rows are independent, which is rank in the rigidity
  matroid. The Molecular Conjecture (Tay–Whiteley; proved by Katoh–Tanigawa)
  establishes this framework for molecular body-hinge frameworks. Matroid
  basis enumeration — listed without chemical application in the enumeration
  candidates doc — acquires one here.

## Open threads

- Laplacian puckering for fused ring systems: transverse field, ring-system
  Laplacian eigenmodes, Cartesian amplitudes, automorphism-adapted
  degeneracies.
- Formalizing the spanning-tree atlas: chart-switching rules, overlap
  conditions, and its rigidity-matroid rank bookkeeping.
- Symmetry adaptation of the redundant primitive set through the existing
  SALC machinery.

## References

- I. Bahar, A. R. Atilgan, B. Erman. Direct evaluation of thermal
  fluctuations in proteins using a single-parameter harmonic potential.
  Fold. Des. 2 (1997) 173–181.
- J. Baker, A. Kessi, B. Delley. The generation and use of delocalized
  internal coordinates in geometry optimization. J. Chem. Phys. 105 (1996)
  192–212.
- L. M. Blumenthal. Theory and Applications of Distance Geometry. Oxford
  University Press, 1953.
- D. Cremer, J. A. Pople. A general definition of ring puckering
  coordinates. J. Am. Chem. Soc. 97 (1975) 1354–1358.
- G. Fogarasi, X. Zhou, P. W. Taylor, P. Pulay. The calculation of ab initio
  molecular geometries: efficient optimization by natural internal
  coordinates. J. Am. Chem. Soc. 114 (1992) 8191–8201.
- C. A. G. Haasnoot. The conformation of six-membered rings described by
  puckering coordinates derived from endocyclic torsion angles. J. Am.
  Chem. Soc. 114 (1992) 882–887.
- A. D. Hill, P. J. Reilly. Puckering coordinates of monocyclic rings by
  triangular decomposition. J. Chem. Inf. Model. 47 (2007) 1031–1035.
- D. G. Kendall, D. Barden, T. K. Carne, H. Le. Shape and Shape Theory.
  Wiley, 1999.
- N. Katoh, S. Tanigawa. A proof of the molecular conjecture. Discrete
  Comput. Geom. 45 (2011) 647–700.
- R. G. Littlejohn, M. Reinsch. Gauge fields in the separation of rotations
  and internal motions in the n-body problem. Rev. Mod. Phys. 69 (1997)
  213–276.
- C. Peng, P. Y. Ayala, H. B. Schlegel, M. J. Frisch. Using redundant
  internal coordinates to optimize equilibrium geometries and transition
  states. J. Comput. Chem. 17 (1996) 49–56.
- P. Pulay, G. Fogarasi, F. Pang, J. E. Boggs. Systematic ab initio gradient
  calculation of molecular geometries, force constants, and dipole moment
  derivatives. J. Am. Chem. Soc. 101 (1979) 2550–2560.
- T.-S. Tay, W. Whiteley. Recent advances in the generic rigidity of
  structures. Structural Topology 9 (1984) 31–38.
- H. Weyl. The Classical Groups: Their Invariants and Representations.
  Princeton University Press, 1946.
- E. B. Wilson, J. C. Decius, P. C. Cross. Molecular Vibrations.
  McGraw-Hill, 1955.
