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


## Review of Stereogenic Bond Arrangements

IUPAC explicitly treats **axial chirality** as a non-planar arrangement of four groups in pairs about a chirality axis, exemplified by allenes and ortho-substituted biphenyls. ([iupac.qmul.ac.uk][1])

### Well-documented classes beyond ordinary alkene cis/trans

| Class                                                |                                        Bond/axis | Real documented stereogenic behavior?                                                           | Examples                                                                                                                                                                                                 |
| ---------------------------------------------------- | -----------------------------------------------: | ----------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Cumulenes / allenes**                              |                   C=C=C, longer cumulated chains | **Yes, classic axial chirality**                                                                | Substituted allenes `abC=C=Ccd` where each terminal carbon has two different substituents; allene chirality was predicted by van ’t Hoff and experimentally demonstrated later. ([Wikipedia][2])         |
| **Biaryls / heterobiaryls**                          |                        C(sp²)–C(sp²) single bond | **Yes, canonical atropisomerism**                                                               | BINOL, BINAP, 6,6′-dinitro-2,2′-diphenic acid; atropisomerism was first experimentally detected in substituted biphenyls by Christie and Kenner in 1922. ([Wikipedia][3])                                |
| **Non-biaryl C(sp²)–C(sp²)**                         | aryl–alkene, aryl–carbonyl, aryl–heteroaryl etc. | **Yes**                                                                                         | “Axially chiral alkenes” are an active class; non-biaryl atropisomers are widely reviewed. ([tetrahedron-chem.com][4])                                                                                   |
| **C–N single bonds**                                 |                  aryl–N, amide C–N, diarylamines | **Yes**                                                                                         | C–N atropisomeric amides, diarylamines, anilides; modern reviews summarize catalytic access to non-biaryl amide C–N atropisomers. ([PMC][5])                                                             |
| **Amide / peptide bonds**                            |              formally C–N single, partial double | **Yes, mostly cis/trans rotamers; sometimes isolable atropisomers**                             | Peptide/proline cis–trans isomerism; tertiary biaryl 2-amides with measurable rotational barriers. ([PMC][6])                                                                                            |
| **N–N single bonds**                                 |                                         N–N axis | **Yes, rarer**                                                                                  | Dixiamycin natural products were reported as N–N atropisomeric natural-product examples; newer atroposelective N–N examples exist. ([MDPI][7])                                                           |
| **C–B single bonds**                                 |                                         C–B axis | **Yes, rare but documented**                                                                    | Sterically hindered aminoarylboranes were reported as first non-biaryl systems with slow C–B rotation and HPLC-resolvable atropisomers. ([RSC Publishing][8])                                            |
| **C–P / C–S / C–O and related aryl–heteroatom axes** |                                  aryl–heteroatom | **Yes, depending on barrier**                                                                   | These fall under heteroatom-containing atropisomerism; C–N and C–P are especially developed in catalysis/ligand chemistry. ([Cell][9])                                                                   |
| **C(sp²)–C(sp³)**                                    |                                  aryl–alkyl axis | **Yes, but less common**                                                                        | Cannabidiol/naphthylcyclohexane derivatives with barriers around 91–95 kJ/mol; cyclazine-core C(sp²)–C(sp³) atropisomers with long calculated/observed configurational lifetimes. ([RSC Publishing][10]) |
| **C(sp³)–C(sp³)**                                    |                        hindered alkane-like axis | **Yes, but special / engineered**                                                               | Triptycene-derived systems and “ethanoid” systems are documented; usually requires very bulky, rigid frameworks. ([knowleslab.princeton.edu][11])                                                        |
| **C–I single bonds**                                 |                                         C–I axis | **Yes, very recent / exotic**                                                                   | Benziodoxole–aryl systems with C–I as sole chiral axis, barriers >30 kcal/mol and very long racemization half-lives. ([Tohoku University][12])                                                           |
| **Metal–metal quadruple bonds**                      |                          M≣≣M, e.g. Re–Re, Mo–Mo | **Torsionally constrained, but not usually treated like ordinary organic bond stereochemistry** | Quadruple bonds have a δ component that favors eclipsed arrangements and creates rotational barriers; Re₂Me₈²⁻ and Mo–Mo systems are classic structural/dynamic examples. ([Chemistry LibreTexts][13])   |

### Your specific proposals

**Hindered single bonds with “3 substituents on each side”**: yes, in principle and in practice. This is **C(sp³)–C(sp³) atropisomerism**. But ordinary substituted ethanes do not qualify because rotation is too fast. You need a rigid, overcrowded framework so the rotamers have lifetimes long enough to isolate or observe separately. Triptycene-derived systems are the classic answer here. ([knowleslab.princeton.edu][11])

**Biaryls with two substituents on each side**: yes, this is the standard textbook/medicinal-chemistry case. But “two substituents on each side” is not the real criterion; the criterion is that the two substituent paths on each side are distinguishable and rotation is slow enough. Ortho substituents are what usually raise the barrier.

**Quadruple bonds with “4 substituents on each side”**: chemically yes, metal–metal quadruple-bond complexes often look like two square-planar fragments joined face-to-face, so there are four ligands around each metal. The δ bond makes rotation non-free and favors particular torsion angles. But I would **not** put these in the same bucket as alkene E/Z or biaryl atropisomerism unless the ligand environment actually makes two configurationally distinct, isolable stereoisomers. The documented mainstream phenomenon is **restricted internal rotation / eclipsed-staggered stereochemistry**, not a general “quadruple-bond cis/trans” rule. ([Royal Society of Chemistry][14])

### Important edge cases

**Triple bonds** are usually not stereogenic as bonds because the local geometry is linear: there is no “same side/opposite side” relationship around the bond. But **cumulenes** behave differently: the stereogenic element is the entire cumulated π-system, not a single triple bond.

**Rings** can make formally single bonds stereogenic in a geometric sense: cis/trans cycloalkanes, fused rings, bridged bicyclics. But that is usually treated as **ring stereochemistry**, not as a freely rotatable bond becoming an atropisomeric axis.

**Amides** are a good borderline case. A normal amide C–N bond is formally single but has partial double-bond character, so cis/trans rotamers can be meaningful. In peptides this is biologically crucial, especially for X–Pro peptide bonds; in more hindered synthetic amides the C–N axis can become configurationally stable enough to isolate atropisomers. ([PMC][6])

### Practical rule for a cheminformatics stereochemistry model

A “stereogenic bond/axis” should not be limited to bond order. Better criteria:

1. There is an **axis or bond** with two distinguishable substituent paths on each side.
2. The relevant conformations are **not related by a fast low-barrier motion** on the timescale you care about.
3. The stereoisomers are either **isolable**, spectroscopically distinguishable, or biologically/materially meaningful.
4. The axis may be C–C, C–N, N–N, C–B, C–P, C–I, metal–metal, or a cumulated π-chain.

So the answer is: **yes, many bonds other than ordinary double bonds and biaryl bonds have documented stereogenic behavior.** The most robust/common are allenes/cumulenes, biaryls/heterobiaryls, hindered C–N amides/anilides, and selected C(sp²)–C(sp³) axes. C(sp³)–C(sp³), C–B, C–I, and metal–metal quadruple-bond cases are real but more specialized.

[1]: https://iupac.qmul.ac.uk/BlueBook/P9.html?utm_source=chatgpt.com "Blue Book P-9"
[2]: https://en.wikipedia.org/wiki/Allene?utm_source=chatgpt.com "Allene"
[3]: https://en.wikipedia.org/wiki/Atropisomer?utm_source=chatgpt.com "Atropisomer"
[4]: https://www.tetrahedron-chem.com/article/S2666-951X%2822%2900005-5/pdf?utm_source=chatgpt.com "Axially chiral alkenes - Tetrahedron Chem"
[5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9572367/?utm_source=chatgpt.com "Construction of Non-Biaryl Atropisomeric Amide Scaffolds ..."
[6]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6222500/?utm_source=chatgpt.com "Role of Steric Factors in Cis/Trans Isomerization - PMC - NIH"
[7]: https://www.mdpi.com/2073-8994/15/6/1261?utm_source=chatgpt.com "Atropselective Organocatalytic Synthesis of Chiral ..."
[8]: https://pubs.rsc.org/en/content/articlelanding/2020/ob/d0ob00421a?utm_source=chatgpt.com "Non-biaryl atropisomerism at the C–B bond in sterically ..."
[9]: https://www.cell.com/chem/fulltext/S2451-9294%2826%2900142-7?utm_source=chatgpt.com "Atropisomerism at the carbon–phosphorus bond: Chem"
[10]: https://pubs.rsc.org/en/content/articlelanding/2023/ob/d3ob01617j?utm_source=chatgpt.com "Atropisomerism about aryl–C(sp3) bonds: chemically ..."
[11]: https://knowleslab.princeton.edu/wp-content/uploads/presentations/2018-06-09-AJM_Atropisomerism_GM_website.pdf?utm_source=chatgpt.com "Atropisomers | The Knowles Group - Princeton University"
[12]: https://tohoku.elsevierpure.com/en/publications/stable-and-responsive-atropisomerism-around-a-carboniodine-bond/?utm_source=chatgpt.com "Stable and responsive atropisomerism around a carbon– ..."
[13]: https://chem.libretexts.org/Courses/Centre_College/CHE_332%3A_Inorganic_Chemistry/06%3A_Coordination_Chemistry-_Bonding/6.05%3A_Metal-Metal_Bonds?utm_source=chatgpt.com "Section 6.5: Metal-Metal Bonds"
[14]: https://books.rsc.org/books/edited-volume/1746/chapter/1119504/Recent-advances-in-the-chemistry-of-metal-metal?utm_source=chatgpt.com "Recent advances in the chemistry of metal–metal quadruple ..."


The physical object is a **torsional rotor system**:

[
\text{end A rotor} ;-; \text{bond axis} ;-; \text{end B rotor}
]

Each end has a set of substituent directions around the bond axis. The stereochemical question is: **what relative angular phases are stable, and are mirror-related phases superposable?**

A useful notation is **n:m rotor type**, where `n` is the number of distinguishable directions around one end of the bond axis and `m` around the other.

## 1. `2:2` coplanar rotor — alkene / imine / amide cis–trans

This is the familiar double-bond case.

Geometry:

```text
A        C
 \      /
  X == Y
 /      \
B        D
```

Both ends are approximately trigonal planar. The π bond forces the two substituent planes to be **coplanar**. The relevant torsion is locked near `0°` or `180°`.

This gives **E/Z or cis/trans**.

Usually not chiral, because the molecule has a stereogenic **plane-like arrangement**, not a helical axis. But it is still a stereogenic bond because exchanging substituents changes the stereoisomer.

Physical examples:

* C=C alkenes
* C=N imines / oximes
* N=N azo compounds / diazenes
* amide C–N partial double bond, especially tertiary amides and peptide X–Pro cis/trans states

So amides belong here if you mean the **amide C–N bond itself**. The geometry is not “free axial chirality”; it is closer to a degraded double bond.

## 2. `2:2` orthogonal rotor — allene / odd cumulene axial chirality

Same number of substituent directions, totally different geometry.

```text
A      C
 \    /
  C==C==C
 /    \
B      D
```

The two terminal substituent planes are approximately **perpendicular**, not coplanar.

This creates an **axis**, not a cis/trans plane. If both ends have distinguishable substituents, the object is chiral:

```text
front pair rotated clockwise relative to back pair  = P
front pair rotated counterclockwise relative to back pair = M
```

This is the cleanest example showing that “four substituents around a stereogenic bond/axis” can produce either:

* planar E/Z stereochemistry, or
* axial P/M chirality,

depending on the allowed torsional phase.

## 3. `2:2` oblique rotor — biaryl and general atropisomerism

Biaryls are not fundamentally special. They are just a `2:2` rotor where the stable torsion is neither coplanar nor orthogonal by electronic necessity, but set by sterics/electronics.

```text
ring plane A   twisted by φ   ring plane B
```

If `φ = 0°` or `180°`, the system may be planar and achiral.

If `φ ≠ 0°/180°`, then the two enantiomeric torsions are:

```text
+φ  and  -φ
```

Those are the `P/M` atropisomers.

This category includes much more than biaryls:

* aryl–aryl
* heteroaryl–aryl
* aryl–amide
* aryl–alkene
* aryl–C(sp³)
* aryl–N
* aryl–O
* aryl–S
* aryl–B
* aryl–P
* even aryl–I in exotic cases

Geometrically, they are all **two-bladed rotors with a frozen skew angle**.

## 4. `2:3` rotor — planar end against tetrahedral/pyramidal end

Now the two ends have different valences around the axis.

One side has two substituent directions, like an sp² carbon or aryl plane. The other has three directions, like a tetrahedral C, pyramidal N, P, S, etc.

```text
      a
      |
A -- axis -- X
      |    / \
      b   c   d
```

This gives a **blade-against-propeller** system.

Physical examples:

* C(sp²)–C(sp³) atropisomers
* aryl–cyclohexyl systems
* aryl–sulfoxide-like or aryl–phosphorus systems, depending on inversion barriers
* some natural-product atropisomers

The interesting thing here is that chirality can arise either from:

1. the bond-axis torsion itself,
2. a true stereocenter at the sp³ atom,
3. or coupling between the two.

So `2:3` systems often blur the line between **axis chirality** and **center chirality**.

If the sp³ end has three different substituents plus the axis bond, then that atom is already a stereocenter. But even if the end is not independently stereogenic, restricted rotation can make the whole assembly stereogenic.

## 5. `3:3` rotor — hindered C(sp³)–C(sp³), N–N, P–P, etc.

This is the ethane-like case, but frozen.

Each end is a three-bladed rotor around the bond axis:

```text
front end:     a, b, c
back end:      d, e, f
axis:          X—Y
```

The allowed torsional motifs are the familiar ethane motifs:

```text
eclipsed
staggered
gauche-like
anti-like, if labels distinguish directions
```

But if rotation is blocked, these become actual stereoisomeric states.

For ordinary ethane, nothing is configurationally stable. But in overcrowded systems, triptycenes, substituted hydrazines, hindered peroxides, biphosphines, etc., the torsional states can become long-lived.

This is not `E/Z` and not ordinary biaryl `P/M`. It is more like **locked conformational stereochemistry of two coupled propellers**.

Important subcases:

### `3:3` with equivalent threefold ends

If both ends are ideal `C3` rotors, then rotation by `120°` is identical. The torsional coordinate has threefold periodicity.

Ethane-like staggered and eclipsed structures may be achiral because they can have inversion, mirror, or improper rotation symmetry.

### `3:3` with labeled substituents

If the substituents are all different, the number of possible stereochemical states explodes. The bond axis can encode:

* helical sense,
* relative ordering of substituents,
* local stereocenters at each end,
* conformational diastereomerism.

So the physical system is richer than `2:2` axial chirality. It is better viewed as **permutation-plus-torsion stereochemistry**.

## 6. N–N single bonds: usually `2/3 : 2/3`, depending on how you count lone pairs

Hydrazine-like N–N systems are geometrically subtle.

Each nitrogen is often pyramidal:

```text
     R
     |
R — N — N — R
    LP   LP
```

Around the N–N axis, each N has:

* two real substituents,
* one lone-pair direction,
* and a pyramidal geometry.

So physically it is close to a `3:3` rotor, except one blade on each side is a lone pair.

That means N–N atropisomerism can be modeled as:

```text
2 real blades + 1 electronic blade
against
2 real blades + 1 electronic blade
```

The complication is **nitrogen inversion**. If inversion is fast, it destroys configurational stability. If inversion is slow or geometrically constrained, the N–N torsional states can become real stereoisomers.

So N–N systems are not merely “single-bond atropisomers”; they are **torsion–inversion coupled rotors**.

That is an important separate class.

## 7. Amides: two different systems depending on which bond you mean

### The amide C–N bond itself

This is a `2:2` coplanar rotor with partial double-bond character:

```text
O
||
C — N
```

The amide resonance forces planarity. The stereochemistry is cis/trans-like, not helical in the ordinary case.

So the amide C–N bond belongs with alkenes, imines, and azo compounds.

### An aryl–N amide axis

Example:

```text
Ar — N(COR)R
```

Now the stereogenic axis may be the **Ar–N bond**, not the amide C–N bond.

That is closer to `2:2` or `2:3` atropisomerism depending on the nitrogen geometry and substituent pattern.

So “amide atropisomerism” can mean two distinct physical things:

1. restricted rotation of the amide C–N partial double bond;
2. restricted rotation around an aryl–N or aryl–C amide-adjacent axis.

Those should not be conflated.

## 8. Quadruple bonds: `4:4` rotor

This is the most interesting one from your list.

A metal–metal quadruple bond is often approximately:

```text
      a           e
      |           |
b — M ==== M — f
      |           |
      c           g
      d           h
```

More abstractly, each metal end is a **four-bladed square rotor**:

```text
front square:  a b c d
back square:   e f g h
```

The δ bond makes the relative rotation electronically meaningful. Unlike a σ single bond, rotation changes overlap of the δ orbital.

The natural torsional states are:

```text
eclipsed square-prismatic       φ = 0°
staggered square-antiprismatic  φ = 45°
general twisted                 φ ≠ 0°, 45°
```

So this is a true `4:4` rotor system.

It is not analogous to alkene E/Z. It is more like a **four-bladed gear pair**.

If the ligand pattern is symmetric, the structure may be achiral despite having a twist. For example, ideal square antiprismatic arrangements can have improper rotations that make them achiral.

But if the four ligands on one or both ends are distinguishable, then the relative phase of the two squares can absolutely become stereogenic. You could get:

* square-prismatic stereoisomers,
* square-antiprismatic stereoisomers,
* helical `Δ/Λ`-like twisted forms,
* many label-permutation isomers.

So physically, quadruple bonds belong to the general family of **higher-order axial rotor stereochemistry**, specifically `4:4`.

## 9. General rotor taxonomy

Here is the broader classification I would use.

| Rotor type       | Geometry                         | Common chemistry                          | Stereochemical behavior                                          |
| ---------------- | -------------------------------- | ----------------------------------------- | ---------------------------------------------------------------- |
| `1:1`            | one substituent each side        | diatomic-like, simple azo with symmetry   | usually no stereogenicity unless lone pairs/coordination matter  |
| `1:2`            | terminal vs trigonal             | imines, oximes, heteroatom multiple bonds | E/Z possible if lone pair counts as direction                    |
| `2:2 planar`     | coplanar four-substituent system | alkenes, imines, amides, azo              | cis/trans, E/Z; usually achiral                                  |
| `2:2 orthogonal` | perpendicular terminal planes    | allenes, odd cumulenes                    | axial chirality, P/M                                             |
| `2:2 oblique`    | twisted two-blade rotors         | biaryls, heterobiaryls                    | atropisomerism, P/M                                              |
| `2:3`            | blade vs propeller               | aryl–alkyl, aryl–N/P/S                    | mixed axial/central/conformational chirality                     |
| `3:3`            | propeller vs propeller           | hindered C–C, N–N, P–P, peroxides         | frozen ethane-like rotamers; torsion/inversion coupling possible |
| `2:4`            | blade vs square                  | organometallic / coordination axes        | possible but less common; label-rich stereochemistry             |
| `3:4`            | propeller vs square              | coordination / cluster systems            | possible higher-order torsional stereochemistry                  |
| `4:4`            | square vs square                 | metal–metal quadruple bonds               | square-prismatic / antiprismatic / twisted rotor stereochemistry |

The key jump is this:

**Organic stereochemistry usually stops at `2:2`; physical stereochemistry does not.**

Once you think in terms of rotors, there is no reason to privilege double bonds and biaryls. They are just low-order cases.

## The deepest distinction: planar vs helical vs geared

You can reduce the examples to three geometric archetypes.

### A. Planar locking

```text
φ = 0° or 180°
```

Examples:

* alkenes
* imines
* azo compounds
* amide C–N bonds

Stereochemistry: **cis/trans**, usually achiral.

### B. Helical two-blade locking

```text
φ = +θ or -θ
```

Examples:

* allenes
* cumulenes
* biaryls
* non-biaryl atropisomers

Stereochemistry: **P/M axial chirality**.

### C. Geared multi-blade locking

```text
n:m rotor with n or m > 2
```

Examples:

* sp³–sp³ hindered bonds: `3:3`
* N–N with lone pairs: pseudo-`3:3`
* metal–metal quadruple bonds: `4:4`

Stereochemistry: **rotameric phase stereochemistry**, often combining helicity, permutation, and local inversion.

This third class is the one mainstream organic stereochemistry underdevelops. It is not exotic physically; it is just harder to name with the usual R/S and E/Z vocabulary.

## How I would encode it

For a molecular graph or stereochemistry engine, I would not start from named classes. I would define a stereogenic rotor as:

```text
(axis atoms A—B,
 substituent directions around A,
 substituent directions around B,
 allowed torsional phase set Φ,
 symmetry group acting on each end,
 inversion operations if present,
 barrier / lifetime metadata)
```

Then ordinary cases fall out:

```text
alkene E/Z:        2:2 rotor, Φ = {0, π}
allene P/M:        2:2 rotor, Φ = {+π/2, -π/2}
biaryl P/M:        2:2 rotor, Φ = {+θ, -θ}
amide cis/trans:   2:2 rotor, Φ = {0, π}, partial double bond
sp3-sp3:           3:3 rotor, Φ = staggered/eclipsed/gauche phases
N-N:               pseudo-3:3 rotor + inversion
quadruple bond:    4:4 rotor, Φ = square-prismatic / antiprismatic / twisted
```

That is the physically honest taxonomy. It also gives you a path to nonstandard stereochemistry without arguing with organic nomenclature.
