# Unconventional stereochemistry types

Status: **Informational**

Codex GPT-5.5 high
 
Prompt: Apart from the tetrahedral, square-planar, trigonal-bipyramidal, octahedral, cis-trans, what other (rarer) arrangements can be found for stereogenic objects (atoms, bonds, fragments with typed ports)? There must be one for each (idealized) bonding geometry, which ones have been so far found in the literature? How are helical, planar (atropisomery, cyclophanes), axial symmetry fit into this model? Do websearches if needed.

Reply

• The key distinction emerging is that an idealized geometry supplies a port set and a symmetry group, but it is stereogenic only for typings that split into multiple symmetry orbits. I’m checking the less-common coordination polyhedra and the formal treatment of non-point stereogenic units now.


• Your “typed ports on an idealized object” model is broad enough to cover essentially all classical stereochemistry. The important correction is:

  > There is not a different abstract group for every bonding geometry.
  > A geometry is characterized by a particular permutation action of its spatial symmetry group on its ports.

  Different geometries may have isomorphic symmetry groups but different port orbits, embeddings, reflection behavior, and rearrangement pathways.

  ## General Model

  Let an idealized stereogenic object have:

  - a set (P) of (n) ports;
  - port types or attached fragments assigned to (P);
  - an orientation-preserving spatial symmetry group
    [
    R(X)\leq S_n
    ]
    acting on those ports.

  Two assignments are the same configuration when related by a proper spatial rotation in (R(X)).

  For (n) distinct port types:

  [
  \text{configurations}\cong S_n/R(X),
  \qquad
  N=\frac{n!}{|R(X)|}.
  ]

  This count includes enantiomers separately because reflections are not allowed when superimposing molecules.

  For repeated port types, the action is no longer free; the configurations must be counted using orbit enumeration, normally Burnside's lemma or Pólya enumeration.

  To identify enantiomeric pairs, extend (R(X)) to the full symmetry group (F(X)), including orientation-reversing operations. Two (R(X))-orbits merged by an element of (F(X)\setminus
  R(X)) form an enantiomeric pair.

  ## Established Coordination Geometries

  The IUPAC Red Book explicitly recognizes the following idealized coordination geometries through coordination number 9. These are all usable as typed-port stereochemical objects.

    CN    Geometry                             IUPAC symbol      Proper port-symmetry group    All ports distinct
  ━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━
     2    linear                               L-2                           (C_2), order 2                     1
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     2    angular                              A-2                           (C_2), order 2                     1
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     3    trigonal planar                      TP-3                          (D_3), order 6                     1
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     3    trigonal pyramidal (Ty)              TPY-3                         (C_3), order 3                     2
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     3    T-shaped (Ts)                        TS-3                          (C_2), order 2                     3
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     4    tetrahedral  (Th)                    T-4                   (T\cong A_4), order 12                     2
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     4    square planar (Sp)                   SP-4                          (D_4), order 8                     3
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     4    square-pyramidal arrangement (Py)    SPY-4                         (C_4) order 4                     6
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     4    seesaw (Sw)                          SS-4                          (C_2), order 2                    12
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     5    trigonal bipyramidal (Tb)            TBPY-5                        (D_3), order 6                    20
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     5    square pyramidal (Sy)                SPY-5                         (C_4), order 4                    30
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     6    octahedral (Oh)                      OC-6                  (O\cong S_4), order 24                    30
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     6    trigonal prismatic (Tr)              TPR-6                         (D_3), order 6                   120
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     7    pentagonal bipyramidal (Pb)          PBPY-7                       (D_5), order 10                   504
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     7    face-capped octahedral (Of)          OCF-7                         (C_3), order 3                  1680
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     7    square-face-capped trigonal prism (Pf) TPRS-7          approximately (C_2), order 2                  2520
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     8    cubic (Cu)                           CU-8                           (O), order 24                  1680
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     8    square antiprismatic (Sa)            SAPR-8                        (D_4), order 8                  5040
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     8    triangular dodecahedral (Dt)         DD-8                          (D_2), order 4                 10080
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     8    hexagonal bipyramidal (Xb)           HBPY-8                       (D_6), order 12                  3360
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     9    tricapped trigonal prismatic (Tt)    TPRS-9                        (D_3), order 6                 60480
  ─────  ───────────────────────────────────  ──────────────  ──────────────────────────────  ────────────────────
     9    heptagonal bipyramidal (Hb)          HBPY-9                       (D_7), order 14                 25920

  For example, a trigonal-prismatic object with six distinct ports has

  [
  \frac{6!}{|D_3|}

  \frac{720}{6}

  120
  ]

  configurations, compared with only

  [
  \frac{6!}{|O|}

  \frac{720}{24}

  30
  ]

  for octahedral geometry.

  The literature also treats trans-bicapped octahedra, several bicapped trigonal prisms, coordination numbers above 9, and distorted or irregular polyhedra. Eight-coordinate square-
  antiprismatic, dodecahedral, and bicapped-trigonal-prismatic isomer enumeration has explicitly been studied using Pólya theory.

  ## The Possible Symmetry Groups

  For finite port arrangements embedded in ordinary three-dimensional space, the proper rotational group must be one of:

  [
  C_n,\qquad D_n,\qquad T,\qquad O,\qquad I,
  ]

  or a subgroup of one of these.

  Thus there are infinitely many possible geometries, but only five general families of finite rotational symmetry groups:

  - (C_n): pyramids, directed rings, many capped structures;
  - (D_n): prisms, antiprisms, bipyramids, planar polygons;
  - (T): tetrahedral family;
  - (O): cube/octahedron family;
  - (I): icosahedron/dodecahedron family.

  Low-symmetry geometries may have only (C_1), meaning every permutation of distinct ports gives a separate configuration.

  ## Cis/Trans Versus Axial Chirality

  Both cis/trans double-bond geometry and axial chirality can be represented using four ports divided into two fixed pairs:

  [
  {{1,2},{3,4}}.
  ]

  The parent relabeling group preserving that partition is

  [
  S_2\wr S_2\cong D_8.
  ]

  Both often produce two configurations, but the two cases differ geometrically.

  ### Double bond

  The substituent pairs are coplanar. The two configurations are:

  [
  E \quad\text{and}\quad Z.
  ]

  These are generally diastereomeric arrangements and individually need not be chiral.

  ### Axial unit

  The two substituent pairs lie in different planes around an axis, as in allenes, cumulenes, biaryls, and some spiro compounds. The two configurations are opposite handednesses:

  [
  R_a/S_a
  \quad\text{or}\quad
  P/M.
  ]

  Thus the abstract two-element quotient alone does not distinguish (E/Z) from axial enantiomerism. The model must also retain the spatial embedding and record which symmetries reverse
  orientation.

  IUPAC describes an axial stereogenic unit as an elongated tetrahedron.

  ## Planar Chirality

  “Planar chirality” does not mean an ordinary planar polygon of ports. A planar polygon can normally be flipped by a proper three-dimensional rotation and is therefore not itself
  chiral.

  A planar stereogenic unit consists of:

  - a reference plane;
  - a directed or ordered path within that plane;
  - an out-of-plane component that distinguishes its two faces.

  Examples include:

  - substituted cyclophanes;
  - ansa compounds;
  - substituted ferrocenes;
  - some constrained cyclic systems;
  - mechanically planar chiral rotaxanes.

  IUPAC formally reduces assignment of (R_p/S_p) planar configuration to a derived tetrahedral stereogenic unit. In the more general port model, it is a typed, oriented plane together
  with a face-selecting component.

  ## Atropisomerism

  Atropisomerism is not itself a geometric class. It is a statement about dynamics:

  > Two conformations are atropisomers when interconversion by rotation is sufficiently hindered that they behave as isolable stereoisomers.

  Depending on the structure, atropisomerism can preserve:

  - axial chirality, especially biaryls;
  - planar chirality;
  - helical chirality;
  - ordinary conformational stereochemistry.

  In the mathematical model, atropisomerism means that an otherwise available continuous path between configurations is removed because its energy barrier is treated as prohibitive.

  ## Helical Chirality

  A helix is not naturally a finite local port arrangement. It is better represented as an oriented path, ordered sequence, or embedded curve.

  Its configuration space has two handed components:

  [
  {\text{right-handed},\text{left-handed}}
  \cong C_2,
  ]

  described by:

  [
  P/M
  \quad\text{or, in coordination chemistry, often}\quad
  \Delta/\Lambda.
  ]

  Examples include:

  - helicenes;
  - helical polymers;
  - propeller-shaped molecules;
  - tris-chelate octahedral complexes;
  - helically folded molecules.

  A helix can be approximated by finitely many ordered ports, but its handedness fundamentally comes from the orientation of a sequence through space, not merely from coloring vertices
  of one polyhedron.

  ## Further Uncommon Stereogenic Objects

  Beyond classical centers, axes, planes, and helices, the literature contains:

  - Conformationally stereogenic rings: chairs, boats, twisted rings, trans-cycloalkenes.
  - Chelate arrangements: (\Delta/\Lambda), ring conformations, fac/mer, and ligand wrapping.
  - Propeller chirality: coordinated aryl groups and molecular rotors.
  - Mechanical planar chirality: rotaxanes.
  - Mechanical axial chirality: catenanes.
  - Topological chirality: molecular knots and some links.
  - Co-conformational stereogenicity: stereochemistry determined by the relative positions of mechanically interlocked components.
  - High-coordinate polyhedral stereochemistry: capped prisms, antiprisms, dodecahedra, cubes, and larger coordination shells.
  - Cluster and cage stereochemistry: typed vertices, edges, faces, or entire subclusters rather than only atom-centered ports.

  Mechanical and topological chirality require extending the model from a single finite port set to embedded graphs, linked components, or knots.

  ## Practical Unified Model

  A sufficiently general stereochemical representation needs four layers:

  1. Combinatorial object: atoms, bonds, fragments, ports, rings, and their connectivity.
  2. Idealized embedding: tetrahedron, prism, axis, plane, helix, knot, etc.
  3. Symmetry action: proper and improper automorphisms acting on the typed components.
  4. Allowed dynamics: rotations, inversions, pseudorotations, ligand exchanges, and their barriers.

  The stereoisomers are then the connected components of the allowed configuration space after quotienting by proper spatial motions and chemically valid automorphisms.

  That formulation covers finite polyhedral stereochemistry, (E/Z), axial and planar chirality, helicity, atropisomerism, fluxional molecules, and mechanically interlocked structures.

  ### Sources

  - IUPAC Red Book: coordination geometries and polyhedral symbols (https://publications.iupac.org/books/rbook/Red_Book_2005.pdf)
  - IUPAC Blue Book P-9: stereogenic centers, axes, planes, and helicity (https://iupac.qmul.ac.uk/BlueBook/P9.html)
  - IUPAC Gold Book: axial chirality (https://goldbook.iupac.org/terms/view/A00547)
  - Pólya enumeration of eight-coordinate complexes (https://doi.org/10.1016/0277-5387(95)00292-Z)
  - Polytopal rearrangement model of stereoisomerization (https://pubmed.ncbi.nlm.nih.gov/39132020/)
  - Chirality in rotaxanes and catenanes (https://pmc.ncbi.nlm.nih.gov/articles/PMC6049620/)
