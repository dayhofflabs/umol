# Alternative stacks

## Summary

* Informational

## Implementation discussion

## Me

I'm thinking about embarking on an ambitious, possibly thankless project developing a new library for representing
molecular structures. We have had a conversation about it before. Here's you concise summary of that conversation. I
have written down some basic ideas about what the project should and shouldn't attempt to do and what the guiding
principles should be. I am a quantum chemist by background and have worked with cheminformatics tools for long enough
that I believe that the current state of the open-source software in the field of molecule building, manipulation,
standardization, and chemical reactions is woefully inadequate. The existing tools, predominantly, RDKit, but also
Openbabel, Indigo Toolkit, CDK are broken to various degrees and cannot support the ecosystem of AI-driven computational
chemistry as much as the field would want them to. So I think a new, well-thought out foundation is not merely a
nice-to-have but an absolute prerequisite for structuring the new methods of scientific discovery. I also think that the
obsession with "performance" has been very detrimental to the field and has produced ugly, bloated, and unmaintainable
software. I would love to see more principled approaches to programming, especially functional programming paradigms
implemented. For these reasons, my current thinking is to developed a Rust core library with a Python interface. My
biggest concerns are:

* Amount of work needed for even a minimal viable implementation
* State of Rust language for scientific and numerical computing
* Trying to do to much, especially regarding the integration of cheminformatics and computational chemistry
* Underestimating the inherent complexity of graph models of molecules and other fundamental concepts and an
  overly optimistic outlook on the feasibility of the project.

## Claude

Based on the materials provided and your background, I'll address your concerns directly while keeping in mind your
preference for functional programming approaches.

Regarding the amount of work for an MVP: Your core proposal is well-scoped. The minimal feature set (SMILES/SDF I/O,
SMARTS matching, reaction SMARTS) represents a solid foundation. The 6+2 month timeline seems optimistic but not
unrealistic if you maintain strict scope control. The key is starting with pure immutable data structures and building
complexity gradually through composition, which aligns well with both Rust's ownership model and functional programming
principles.

The state of Rust for scientific computing is maturing rapidly. The petgraph library provides solid foundations for
graph algorithms, though you'll likely need to extend it for chemistry-specific operations. The ndarray ecosystem, while
not as mature as NumPy, is sufficient for the core computational needs of molecular representation. The PyO3 bindings
are production-ready and should handle your Python interface requirements effectively.

Your concern about the integration of cheminformatics and computational chemistry is well-founded. The design document
wisely lists this as a future goal rather than a core requirement. I would suggest taking an even more conservative
approach: start with a pure graph representation and treat 3D structures as a separate concern entirely for the MVP.
This separation of concerns will let you build a solid foundation without getting bogged down in the complexities of
interconversion.

The complexity of molecular graph models is indeed often underestimated. However, your non-goals section shows good
awareness of these pitfalls. The explicit rejection of a unified internal representation is particularly important -
accepting lossy conversions between different models is more honest than pretending all chemical structures can fit one
paradigm.

Would you like me to elaborate on any of these points or discuss specific technical approaches to maintaining functional
purity while dealing with chemical graph representations?

## Me

Which alternative stacks (language + library + structuring approach) can you envision for such a project? Can you
contrast and compare them to the proposed Rust + Python approach? Important points of comparisons are, roughly in
decreasing order of importance:

* Clear mapping domain concepts onto implementation
* As much immutability as possible, with an eye to multithreaded execution
* Speed of development
* Speed of execution: it should be fast enough to be useful, not beat any records. RDKit is not very fast, so
  same order of magnitude is desirable

## Claude

From your criteria and background, I see three viable alternative stacks worth considering alongside Rust+Python:

OCaml + Python The OCaml ecosystem offers excellent support for building typed, immutable data structures with pattern
matching. The Jane Street Core library provides industrial-strength foundational tools. OCaml's module system enables
clear separation of domain concepts through functors. Its type system can encode complex domain constraints at compile
time. Python bindings are possible through pyml, though less mature than PyO3. Development speed would be comparable to
Rust for developers familiar with ML-family languages. Performance would be adequate, with good garbage collection
characteristics for scientific workloads.

Scala 3 (Pure) Scala 3's intersection types and union types enable precise modeling of chemical domain concepts. The
cats-effect ecosystem provides principled handling of effects and concurrency. Immutability comes naturally through case
classes and functional collections. Development could be faster than Rust due to the JVM ecosystem's maturity. Python
interop through py4j is stable but introduces runtime overhead. The JVM provides acceptable performance with some
garbage collection tuning required.

Haskell + Python Haskell enables perhaps the cleanest mapping of domain concepts through its powerful type system and
type classes. Libraries like containers and vector provide highly optimized immutable data structures. The strict
separation of pure and effectful code enforces good architecture. Development speed depends heavily on team experience
with advanced type-level programming. Python bindings through inline-python exist but are not industrial-strength.
Performance can match C++ when carefully optimized.

Comparing to Rust+Python: - Domain mapping: Rust's ownership system and enums work well for chemical structures, but
OCaml and Haskell offer cleaner abstractions through more sophisticated type systems - Immutability: All options support
it well, with Rust providing unique compile-time guarantees through borrowing - Development speed: Rust likely slower
initially due to fighting the borrow checker, but could pay off in reduced debugging time - Execution speed: All options
can meet the performance requirements, with Rust potentially offering the most predictable performance characteristics

The key advantage of Rust+Python remains its growing scientific computing ecosystem and excellent Python integration
through PyO3. While the alternatives offer some theoretical advantages in domain modeling, the practical benefits of
Rust's tooling and community momentum shouldn't be underestimated.

Would you like me to elaborate on any specific aspect of these alternatives?

## Me

I have some prior experience with Scala (although less so with Scala 3) but would very much like to avoid the JVM
ecosystem. I've also tried working with Haskell and was deeply impressed but my typelevel programming skills are still
rudimentary. I don't have any experience with Ocaml or any other ML-derived languages so that learning curve is probably
not worth taking, even if it sounds really interesting. Pure Python will likely not be fast enough even in the best of
cases (using a C++ graph library like igraph? Does it still exist? NetworkX is definitely not fast), correct? So Rust +
Python seems like a good compromise at this point, and has the momentum on its side.

## Claude

Your assessment aligns with my view of the tradeoffs. Pure Python would indeed struggle to meet performance requirements -
NetworkX remains quite slow, and while igraph exists and is faster, its Python bindings aren't ideal for building a
foundational chemistry library. The C++ bindings add complexity without providing the type safety and memory guarantees
that would make them worth the effort. Your experience with Scala but desire to avoid the JVM, combined with an interest
but not deep expertise in Haskell, makes Rust a pragmatic choice. Rust offers many of the benefits of functional
programming - immutability by default, powerful pattern matching, and a strong type system - while maintaining more
approachable learning curve than Haskell's advanced type-level programming. The momentum behind Rust in scientific
computing is significant. The development of libraries like PyO3, petgraph, and ndarray shows a clear trajectory toward
a mature scientific computing ecosystem. While these libraries may not yet match the sophistication of their
counterparts in more established languages, they provide a solid foundation that's actively improving. Your recognition
that the learning curve of OCaml, while potentially beneficial, might not justify the investment is sound. The practical
benefits of Rust's growing ecosystem and excellent documentation likely outweigh the theoretical advantages of a more
mathematically pure approach.

Would you like to explore specific aspects of how Rust's type system could model chemical concepts, particularly in
areas where you've found Haskell's type-level programming appealing?
