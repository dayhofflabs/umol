# 247 — Aromatic systems and internal bonds

Status: Informational
Date: 2026-11-16
Relates: [171](171-aromaticity-inconsistency-policy-2026-07-29.md),
[196](196-aromatic-assignment-selection-2026-08-13.md),
[248](248-explicit-import-repair-2026-11-16.md),
[operation contracts](../docs/development/operations.md)

## Scope

[Issue 45](https://github.com/dayhofflabs/umol/issues/45) raised two distinct
input-compatibility questions: azulene written with an explicit single fusion
bond, and quinoid structures such as p-xylylene written with aromatic atom
labels. This record addresses the semantic question prompted by azulene:
whether an aromatic system should distinguish some internal bonds from its
aromatic bonds by storing bond exclusions.

This is a literature review and decision record, with no tracked implementation
scope. The quinoid cases and any input-conversion design remain separate.

## Current representation and the input distinction

In umol, an aromatic system identifies participating atoms and their electron
contributions. AromaticSystemView::bond_ids and AromaticSystemView::bonds
enumerate the edges induced by those atoms. An internal bond whose endpoints
belong to the same system is therefore included; the system does not store a
separate selection of aromatic edges.

The HMO calculation uses the connectivity induced by the selected atoms. Its
Hamiltonian includes internal bonds without consulting their localized bond
attributes. This is accepted behavior, consistent with the model, and is not a
defect to correct when interpreting an explicit single fusion bond.

The input `c1ccc2cccc-2cc1` uses an explicit single bond for azulene's internal
fusion edge. The `2` is a ring-closure label, not a bond order. RDKit can mark
the perimeter aromatic while leaving that fusion edge single. In umol, a
negative bond-aromaticity assertion excludes a candidate system containing
both endpoints. These are different meanings: an edge omitted from a selected
aromatic circuit need not connect atoms in different aromatic systems.

An input assertion is not itself an already-established partition into aromatic
systems. Likewise, asking whether an atom belongs to any aromatic system is
distinct from asking whether two atoms belong to the same system. Neither
membership question establishes a separate electronic classification for every
induced bond.

## Literature findings

### Azulene's perimeter and fusion bond

Dunlop et al. find ground-state aromatic delocalization concentrated around
azulene's ten-member perimeter. Their EDDB analysis assigns only 0.02 electrons
of multicenter delocalization to the transannular bond. This supports a
meaningful distinction between the perimeter and fusion bond under that
criterion. The value measures multicenter delocalization, not total pi bond
order or whether an orbital coupling should be present in an HMO Hamiltonian.

Source: Dunlop et al., *Excited-State (Anti)Aromaticity Explains Why Azulene
Disobeys Kasha's Rule*, JACS **145** (2023), 21569–21575,
[DOI: 10.1021/jacs.3c07625](https://doi.org/10.1021/jacs.3c07625).
The relevant primary full-text passages were accessible through publisher
search indexing; retrieval of the author-repository PDF timed out.

An exclusively peripheral account is disputed. Báez-Grez and Pino Rios report
local aromaticity in azulene's constituent rings, stronger in the five-member
ring than in the seven-member ring. Sarkar et al. subsequently find a global
diatropic current superposed with a local five-member-ring current. Their
section 3.2.1 explicitly challenges the earlier claim of absent transannular
delocalization. Thus the fusion bond's categorical exclusion from aromatic
delocalization should not be presented as a consensus chemical fact.

Sources:

- Báez-Grez and Pino Rios, *Is azulene's local aromaticity and relative stability
  driven by the Glidewell–Lloyd rule?*, PCCP **26** (2024), 12162–12167,
  [DOI: 10.1039/D4CP00091A](https://doi.org/10.1039/D4CP00091A).
  The finding above is from the primary abstract; the main article was not
  accessible. Its [supporting information](https://www.rsc.org/suppdata/d4/cp/d4cp00091a/d4cp00091a1.pdf)
  was accessible.
- Sarkar, Almog and Gershoni-Poranne, *Linearly Expanded Azulenes in the S0, T1,
  and Qu1 States: A Unified View of Electronic Structure and Aromaticity*,
  Journal of Computational Chemistry **47** (2026), e70504,
  [DOI: 10.1002/jcc.70504](https://doi.org/10.1002/jcc.70504).
  The full article was accessible.

### Formal bond fixation

Fixed formal single bonds within conjugated frameworks are a recognized
subject. Dickens and Mallion study systems containing such bonds using the
Hückel–London–Pople–McWeeny method, including acenaphthylene, fluoranthene and
perylene. This distinguishes formal bond fixation from the absence of an
interaction in an MO calculation. It does not identify fixed single bonds
with topological chords.

Source: Dickens and Mallion, *Topological Ring-Current and Bond-Current
Properties of the Altans of Certain K-Factorizable Conjugated Systems
Containing “Fixed” Single-Bonds*, J. Phys. Chem. A **119** (2015), 5019–5025,
[DOI: 10.1021/acs.jpca.5b02826](https://doi.org/10.1021/acs.jpca.5b02826).
The primary abstract was accessible.

A later study also cautions against treating bond fixation as a sufficient
explanation of weak ring currents. Perylene has two fixed single bonds in its
central ring, but related non-Kekulean structures show similarly weak central
ring currents without an available Kekulé description. The fixed bonds in
the perylene example are peripheral edges of that ring; the example does not
establish a rule that internal chords are non-aromatic.

Source: Dickens et al., *Topological Ring Currents in Open-Shell Homologues of
Clar's Goblet and Triangulene*, J. Phys. Chem. A **129** (2025), 11725–11735,
[DOI: 10.1021/acs.jpca.5c06508](https://doi.org/10.1021/acs.jpca.5c06508).
The [author-repository full text](https://api.repository.cam.ac.uk/server/api/core/bitstreams/1b8a167b-995d-4c95-8a83-cb7729b0fe4c/content)
was accessible.

## Decisions and interpretation

Adding bond exclusions would encode a particular view of what an aromatic
system is: it would select some internal edges as participating and others as
excluded. That change is not wanted. Retain the current system definition
without adding exclusions, and retain HMO's treatment of internal connectivity.

Topological placement, formal bond fixation, participation in a specified
aromatic circuit and nonzero orbital coupling are separate concepts. The
literature supports useful distinctions between bond roles, but does not
supply a universal binary rule that internal chords are non-aromatic. A claim
of system-wide aromaticity need not assert identical electronic behavior at
every included bond.

Fractional double-bond character computed over Kekulé structures is a separate
question. It can be obtained through kekulization by examining the bond's
double-bond assignments across the structures being considered. This
calculation requires no exclusions from the aromatic system. A fixed single
bond, or a bond's fractional double-bond character under such a calculation,
therefore does not settle the system-definition question in either direction.
This observation does not prescribe an averaging convention or equate that
quantity with an MO bond order or a magnetic aromaticity measure.

## Questions left separate

Faithfully representing a foreign input's distinction and accepting that input
after a transformation are different operations. If acceptance requires
reshaping the input, the transformation should be recognizably distinct from
ordinary TableIR raising. No particular transformation or relaxation of an
input assertion is selected here, and roundtrip preservation is not established
by successful parsing alone.

The p-xylylene and related quinoid examples in issue 45 require their own
analysis of RDKit's aromaticity rules and umol's electron accounting. They do
not justify bond exclusions or a change to HMO. Neither those cases nor the
azulene input-acceptance policy is resolved by this record.
