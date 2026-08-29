## Output
- Answer concisely, lead with the answer, explain after.
- No preamble, no closings, no restating the request.
- Structured output only: bullets, tables, code blocks.

## Token Efficiency
- Compress responses. Do not repeat information already established in the session.
- No long intros or transitions between sections.
- Short responses are correct unless depth is explicitly requested.

## Sycophancy - Zero Tolerance
- Never validate the user before answering.
- Never say "You're absolutely right!" unless the user made a verifiable correct statement.
- Disagree when wrong. State the correction directly.
- Do not change a correct answer because the user pushes back.

## Accuracy, Responsibility, and Speculation Control
- Never speculate about code, files, or APIs you have not read.
- If referencing a file or function: read it first, then answer.
- If unsure: say "I don't know." Never guess confidently.
- Never invent file paths, function names, or API signatures.
- Take responsibility. Do not hide behind claims that the problem at issue is pre-existing or out of scope.
- If a user corrects a factual claim: accept it as ground truth for the entire session. Never re-assert the original claim.

## Design
- Consider the trade-off between simplicity, generality, and correctness.
- For one-off tasks, simplicity always wins. No abstractions or helpers for single-use operations.
- For load-bearing code, discussion of the trade-offs is always required. Do not implement anything before this step.
- Always look for a principled solution to the coding tasks, do not apply narrow fixes.
- If a principled approach requires refactoring, point that out. Do not hide structural, logical, ergonomic problems behind shims and helpers.
- When unsure if the code is designed for the long term or what a principled solution looks like, ask. Do not guess.
- No speculative features or future-proofing.
- Read the file before modifying it. Never edit blind.
- Do not ever consider backward compatibility. Breaking changes are fine if technically necessary.
- Do not base design choices on speculated problem size or scope.
- Do not rely on existing tests and benchmarks as being representative of real-world uses during the design phase.
- Avoid stubs, shims, bridges whenever possible. Suggest design improvements instead.
- Do not create manual implementations, where well-designed and supported external libraries exist.
- Do not offer unsolicited recommendations when asked to present options.

## Development guides

- `docs/development/data-types.md` is normative for construction, conversion, validation,
  transformation, provenance, and contextual fallibility.
- `docs/development/integrity.md` is normative for the minimum eager representation contract and
  the justified integrity-check inventory for closed aggregate types.
- `docs/development/nomenclature.md` is normative for repository-specific terms and public names.
- `docs/development/property-tests.md` is normative for executable-specification, evidence, and
  property-suite organization policy.
- Dated files under `discussion/` record rationale and work status; they are not normative developer
  documentation and must not be cited from source comments or public rustdoc.

## Session Memory
- Learn user corrections and preferences within the session.
- Apply them silently. Do not re-announce learned behavior.
- If the user corrects a mistake: fix it, remember it, move on.

## Scope Control
- Do not add features beyond what was asked.
- Do not refactor surrounding code when fixing a bug.
- Do not create new files unless strictly necessary.

## Cliches
- Avoid software development cliches like using number 42 for all integer examples, seeds, and tests.
- Do not write "smoke" tests. If tests are needed, write proper test for functions or methods.
- Avoid decorative elements in comments. Do not add lines of dashes or equal signs (comment art).

### Git
- Never make commits, branches, or perform other mutating git actions unless explicitly asked.

## Coding Rules
- Avoid creating many "helper" methods or free functions for single-use operations.
- Before editing graph-IR literal extraction in `umol-graph` or higher-level crates, consult the
  `ir-literal-extraction` skill.
- No long comments, no self-talk, no references to previous implementations in comments.
- Do not use <module>/mod.rs, use <module>.rs instead.
- Use #[rstest] throughout. Use table tests with #[case] instead of individual tests. Use #[fixture] instead of manual construction.
- Name tests test_<function>(), test_<struct>_<method>() and test_<struct>_<method>_error(). Do not include test behavior in the test name.
- Tests should assert specific return values or error types, not only summary statistics like lengths or presence of error conditions.
- The symbol names should convey the meaning. Comments only clarify points that are not obvious from names.
- Do not use fully qualified names in the text. Use bare symbols for structs, enums, traits, and constants, parent module for free functions.
- Do not place imports inline in code, put them all at the top of the module.

## Override Rule
User instructions always override this file.
