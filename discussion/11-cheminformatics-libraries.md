# Cheminformatics Libraries

## [Awesome Cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics)

## [MayaChemTools](http://www.mayachemtools.org)

MayaChemTools is a growing collection of Perl and Python scripts, modules, and classes to support a variety of day-to-day computational discovery needs.

The command line Python scripts based on Psi4 provide functionality for the following tasks:

* Calculation of single point energies
* Calculation of interaction energies
* Calculation of molecular properties and partial charges
* Performing structure minimization
* Generating molecular conformations
* Performing torsion scan
* Visualizing frontier molecular orbitals and dual descriptors
* Visualizing electrostatic potential on densities and molecular surfaces

The command line Python scripts based on AutoDock Vina provide functionality for the following tasks:

* Performing rigid and flexible docking
* Scoring molecules

## Python tools of interest

## [CGRtools](https://github.com/cimm-kzn/CGRtools)

* Paper: [CGRtools: Python Library for Molecule, Reaction, and Condensed Graph of Reaction Processing. Journal of Chemical Information and Modeling 2019 59 (6), 2516-2521.](https://doi.org/10.1021/acs.jcim.9b00102)

Condensed Graph of Reaction (CGR) -> Is this the same as bond shift graph?

## [chython](https://github.com/chython/chython)

* Fork of CGRtools

## C++

## [CDPKit](https://cdpkit.org)

CDPKit comes bundled with a set of ready-to-use command line tools and GUI applications that help to master recurring tasks in modern computer-aided drug design (CADD) workflows. After installation, the provided programs can be found in the Bin sub-folder of the CDPKit installation directory. Currently, the following application areas are covered:

* Substructure searching
* 3D structure generation
* Conformer ensemble generation
* Tautomer generation and standardization
* Stereoisomer enumeration

## [LillyMol](https://github.com/EliLillyCo/LillyMol)

* Own description
    LillyMol has some novel approaches to substructure searching, reaction enumeration and chemical similarity. These have been developed over many years, driven by the needs of Computational and Medicinal Chemists at Lilly and elsewhere.

* Molecule building seems to use a deferred approach (set_modified, _partially_built flags)
* Uses atom indices ("numbers") to manipulate molecules.

## Julia

## [MolecularGraph.jl](https://github.com/mojaie/MolecularGraph.jl)

* MIT license
* SDF, SMILES parser -> useful?

### Question:

In a separate conversation, we're have been working on an implementation of molecular structures in Rust. In preparation for the implementation, I've been collecting code samples from existing implementations for inspiration. Here is the code from the MolecularGraph.jl library. I don't know Julia at all, although I can guess a lot from Python/Perl experience. I'd like to ask you to analyze this implementation of the base MolGraph class , in particular with respect to its approach to storage and manipulation of the molecular graph structure. If I understand it correctly, it is fundamentally a mutable structure with cached properties, much like the RDKit implementation. Is that correct? What are the modes for building the structure from scratch, not via parsing a string representation like SMILES? How well-develeoped is that interface? What other key design patterns can you recorgnize in this implementation? In your estimation, is it a useful model for designing a modern molecule building interface, regardless of programming language? Let me know if I need to include additional files.

### Partial Response:

The MolGraph is indeed a mutable structure with cached properties:
```julia
struct MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::Dict{Symbol,Any}
    state::Dict{Symbol,Any}
    edge_rank::Dict{Edge{T},Int}
end
```
It uses:

A graph backbone (SimpleGraph) - likely from LightGraphs.jl
Dictionaries for vertex and edge properties
A state dictionary for caching and modification tracking
The :has_updates flag to trigger recalculation of derived properties

### Building from Scratch

There are several constructors:

1. Empty constructor:
```julia
MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph() = MolGraph{Int,Any,Any}()
```
2. From edge/property lists:
```julia
MolGraph{T,V,E}(edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...)
```
3. Incremental building:
```julia
add_vertex!(mol::MolGraph{T,V,E}, prop::V)
add_u_edge!(mol::MolGraph{T,V,E}, e::Edge, prop::E)
```
### Key Design Patterns

1. Type Parameterization: Generic over vertex index type T, vertex property type V, and edge property type E
2. Dispatch-based Event System:
```julia
dispatch!(mol, event) = mol.state[event](mol)
```
Uses a callback mechanism for initialization and updates
3. Lazy Evaluation with Caching:
```julia
get_state(mol, :has_updates) && dispatch!(mol, :updater)
```
Properties are recomputed only when needed
4. Dictionary-based Property Storage rather than parallel arrays
5. Graph Manipulation Methods that maintain property consistency

### Validation seems to be using callbacks

## Rust

## [molecule.rs](https://github.com/rapodaca/molecule.rs)

* MIT license
* Minimal implementation in Rust. Blog post: [Cheminformatics in Rust: Implementing a Minimal Molecule API](https://depth-first.com/articles/2020/05/11/cheminformatics-in-rust-implementing-a-minimal-molecule-api/)
* Asserts that the minimal molecular interface contains the following getter methods.

    - element. Returns the Element associated with an atom given its integer ID. Element provides methods to query atomic number, symbol, and other attributes.
    - isocomp. Returns an Isocomp ("isotopic composition") object associated with an atom given its integer ID. Isocomp represents the state of an atom's nucleus. At a minimum, it reports whether the composition is naturally occurring and iterates the set of mass numbers.
    - electrons. Returns the number of nonbonding electrons associated with an atom given its integer ID.
    - hydrogens. Returns the number of virtual hydrogens associated with an atom given its integer ID.
    - charge. Returns the formal charge of a member atom identified by integer id.
    - atomParity. Returns the parity of an atom identify by integer id. Atom parity, encoded as a three-member enumeration, expresses the tetrahedral arrangement of substituents as clockwise, counterclockwise, or undefined.
    - bondOrder. Returns the formal bond order between two atoms identified by their integer ids.
    - bondParity. Returns the parity of a bond between two atoms identified by integer ids. Bond parity, encoded as a three-member enumeration, expresses the conformation about a double bond as syn, anti, or undefined.

* Assumes that the molecule is predefined, seems to completely ignore manipulation and validation aspects of it.
* __Atom and Bond Parity__

    Tetrahedral stereochemical configuration is supported through the concept of atom parity, as defined by the atomParity method. This idea borrows from notions of atom parity present in both SMILES and the molfile format. View the axis defined by the bond from the neighbor with the __smallest ID__ to a central atom. If the remaining neighbors, ordered by ascending ID, sweep clockwise, atom parity is positive (e.g., return value "+1"). Otherwise, atom parity is negative (e.g., return value "-1"). Atoms without parity assignments return a null value (e.g., "0").

    E/Z double bond conformations are supported through the bond parity concept, as defined by the bondParity method. Identify the two neighbors of the double bond system having the lowest-valued IDs. If these neighbors appear on the same side of the double bond, parity is syn (e.g., "+1"), otherwise atom parity is anti (e.g., "-1"). A double bond without parity assignment returns a null value (e.g., "0").

    Atom and bond parities involving atoms bearing one virtual hydrogen can be computed by assigning it the lowest atom ID. In the case of atom parity, the axis defined by the virtual hydrogen and the central atom would be sighted. In the case of bond parity, hydrogen could be assigned the role of lowest-numbered neighbor ID.

    Should an atom or bond not fit the template, the corresponding parity must be returned as null (e.g., "0"). Non-template descriptions of atomic configuration or bond conformation can be realized via Rich Expression Constructs as described below.

    __Both atom and bond parity rely on stable atomic IDs__. After an atom is assigned an ID, it can never change in a way that affects relative ordering to neighbors. Fortunately, this constraint is easy to achieve.

* Claims that atoms and bonds are not unnecessary as concepts, can be replaced by atom indices or pairs of atom indices. Don't understand how he can make this big claim without having written necessary code. Seems stupid.

* Molecule extends graph. Fair enough but not natural; the terminology is different, which should create some strong friction.

## [chiral-data](https://github.com/chiral-data/)

* MIT license
* Seems to belong to a company? 4 years old, seems to be mostly wrapping OpenBabel. 
* Has a [SMILES writer](https://github.com/chiral-data/rust-chem/blob/main/src/smiles_writer.rs) of questionable quality?
* __Has graph symmetry perception and canonicalization algorithms [Link](https://github.com/chiral-data/rust-graph-symmetry).__

## [molrs](https://github.com/molrs/molrs-core)

* MIT license
* Pretty simple implementation, seems to base it on vectors of atoms, bonds, and rings.
* Bonds can be single, double, up, down, ... seems to be borrowed from RDKit.
* Has [SMILES parser](https://github.com/molrs/molrs-core/blob/main/src/molecule.rs) and some formatter.

## OCaml

## [consent](https://github.com/UnixJunkie/consent)

* GPL 3.0 license
* Paper: [Berenger, F., Zhang, K.Y.J. & Yamanishi, Y. Chemoinformatics and structural bioinformatics in OCaml. J Cheminform 11, 10 (2019).](https://doi.org/10.1186/s13321-019-0332-0)
* Mostly intro to OCaml, advantages of concise notation and mature language.
* Needs to interface via files to existing ecosystem

## [molenc](https://ocaml.org/p/molenc/16.13.0)

Molecular encoder/featurizer using rdkit and OCaml

## Scala

## [chemf](https://github.com/stefan-hoeck/chemf)

* GPL 3.0 license
* Paper: [Höck, S., Riedl, R. chemf: A purely functional chemistry toolkit. J Cheminform 4, 38 (2012)].(https://doi.org/10.1186/1758-2946-4-38).
* Immutable data structure for molecular graphs
* Complete SMILES parser in accordance with the OpenSMILES specification.
* Labeled graph representation for molecules, parametrized over the vertex and edge label types.
```scala
trait LGraph[+E,+V] {
    def graph: Graph
    def vLabel (v: Int): V
    def eLabel (e: Edge): E
}
```
* Implemented Functor, Foldable, Traverse for LGraph.
* Concrete implementation in private case class
```scala
private case class LgImpl[+E,+V] (
    graph: Graph,
    vertices: IndexedSeq[V],
    eMap: Map[Edge,E]
) extends LGraph[E,V] {
    def vLabel (v: Int) = vertices (v)
    def eLabel (e: Edge) = eMap (e)
}
```
* Element
```scala
sealed abstract class Element (
    val atomicNr: Int
)
object Element {
    case object H extends Element (1)
    case object He extends Element (2)
...
}
```
* Isotope
```scala
sealed trait Isotope {
    def element: Element
    def massNumber: Option[Int]
}
object Isotope {
    def apply (e: Element): Isotope = ...
    def apply (e: Element, mn: Int)
        : Isotope = ...
}
```
* Atom
```scala
case class Atom (
    isotope: Isotope,
    charge: Int,
    hydrogens: Int,
    stereo: Stereo
)
```
* Bond
```scala
sealed trait Bond
object Bond {
    case object Single extends Bond
    case object Double extends Bond
...
}
```
* Molecule
```scala
type Molecule = LGraph[Bond,Atom]
```
* Apart from using Isotope as (an unnecessary) intermediary between Atom and Molecule, very similar construction.
* Use `foldMap` to compute additive molecular properties (such as molecular mass) from atomic properties - nice but mostly academic. There are few properties like that.
* __Sum formula as a special molecule type__
```scala
type Formula = Map[Isotope,Int]
```
* __Functional SMILES parser__ -> useful
* SmilesAtom different from Atom
* Error handling using scalaz::Validation (similar to Either)


## Haskell

## [smiles](https://github.com/zmactep/smiles)

Provides parsing of OpenSMILES [spec](http://opensmiles.org/opensmiles.html) (SMILES and SMARTS?) using MegaParsec.

## [radium](https://github.com/klangner/radium)

Radium is Haskell library for the Chemistry. It has the following functionality:
* Periodic table with the element data.
* Readers and writers for the following formats: SMILES (examples), [Condensed](http://en.wikipedia.org/wiki/Structural_formula#Condensed_formulas).

## [ouch](https://github.com/odj/Ouch)

Over 15 year old. Uses Parsec for parsing SMILES.
