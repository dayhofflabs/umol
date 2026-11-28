import pytest

from umol import (
    AutomorphismAlgorithm,
    Molecule,
    PreparedSubstructureHost,
    PreparedSubstructurePattern,
    RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
    SubstructureSearchConfig,
)


def test_substructure_search_config_default():
    config = SubstructureSearchConfig()

    assert config == SubstructureSearchConfig.default()
    assert config.match_algorithm == SubstructureMatchAlgorithm.GraphAndOverlays()
    assert (
        config.subgraph_isomorphism_algorithm
        == SubgraphIsomorphismAlgorithm.BitsetRi()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert config.automorphism_algorithm == AutomorphismAlgorithm.TrailIr()
    assert config.reduce_pattern_symmetry is False
    assert config.reduce_host_symmetry is False
    assert repr(config) == (
        "SubstructureSearchConfig("
        "match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(), "
        "subgraph_isomorphism_algorithm="
        "SubgraphIsomorphismAlgorithm.BitsetRi(), "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara(), "
        "automorphism_algorithm=AutomorphismAlgorithm.TrailIr(), "
        "reduce_pattern_symmetry=False, "
        "reduce_host_symmetry=False)"
    )


@pytest.mark.parametrize(
    (
        "match_algorithm",
        "subgraph_isomorphism_algorithm",
        "automorphism_algorithm",
        "reduce_pattern_symmetry",
        "reduce_host_symmetry",
    ),
    [
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.BitsetRi(),
            AutomorphismAlgorithm.TrailIr(),
            True,
            False,
            id="graph-bitset-ri-pattern-sym",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.Incidence(),
            SubgraphIsomorphismAlgorithm.Glasgow(),
            AutomorphismAlgorithm.Nauty(),
            True,
            True,
            id="incidence-glasgow-both-sym",
        ),
    ],
)
def test_substructure_search_config_value(
    match_algorithm,
    subgraph_isomorphism_algorithm,
    automorphism_algorithm,
    reduce_pattern_symmetry,
    reduce_host_symmetry,
):
    config = SubstructureSearchConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subgraph_isomorphism_algorithm,
        relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara(),
        automorphism_algorithm=automorphism_algorithm,
        reduce_pattern_symmetry=reduce_pattern_symmetry,
        reduce_host_symmetry=reduce_host_symmetry,
    )

    assert config.match_algorithm == match_algorithm
    assert config.subgraph_isomorphism_algorithm == subgraph_isomorphism_algorithm
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert config.automorphism_algorithm == automorphism_algorithm
    assert config.reduce_pattern_symmetry is reduce_pattern_symmetry
    assert config.reduce_host_symmetry is reduce_host_symmetry
    assert repr(config) == (
        "SubstructureSearchConfig("
        f"match_algorithm={match_algorithm!r}, "
        f"subgraph_isomorphism_algorithm={subgraph_isomorphism_algorithm!r}, "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara(), "
        f"automorphism_algorithm={automorphism_algorithm!r}, "
        f"reduce_pattern_symmetry={reduce_pattern_symmetry!r}, "
        f"reduce_host_symmetry={reduce_host_symmetry!r})"
    )
    with pytest.raises(AttributeError):
        config.reduce_pattern_symmetry = False


def test_substructure_search_config_keyword_error():
    with pytest.raises(
        TypeError,
        match=(
            "^SubstructureSearchConfig.__new__\\(\\) takes 0 positional arguments "
            "but 2 were given$"
        ),
    ):
        SubstructureSearchConfig(
            SubstructureMatchAlgorithm.Incidence(),
            SubgraphIsomorphismAlgorithm.Glasgow(),
        )


def test_molecule_substructure_matches():
    pattern_source = '{:atoms ["C" "C"] :bonds [[0 1 "1"]]}'
    host_source = '{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}'
    pattern = Molecule.parse(pattern_source)
    host = Molecule.parse(host_source)
    pattern_before = Molecule.parse(pattern_source)
    host_before = Molecule.parse(host_source)

    matches = pattern.substructure_matches(host)

    assert isinstance(matches, list)
    assert [match.atoms.matched_pairs for match in matches] == [
        ((0, 0), (1, 1)),
        ((0, 1), (1, 0)),
    ]
    assert [match.bonds.matched_pairs for match in matches] == [((0, 0),), ((0, 0),)]
    assert [match.dative_bonds.matched_pairs for match in matches] == [(), ()]
    assert [match.aromatic_systems.matched_pairs for match in matches] == [(), ()]
    assert [match.multicenter_bonds.matched_pairs for match in matches] == [(), ()]
    assert [match.noncovalent_bonds.matched_pairs for match in matches] == [(), ()]
    assert [match.stereo_atoms.matched_pairs for match in matches] == [(), ()]
    assert [match.stereo_bonds.matched_pairs for match in matches] == [(), ()]
    assert pattern == pattern_before
    assert host == host_before

    with pytest.raises(AttributeError):
        matches[0].atoms.matched_pairs.append((1, 2))
    assert matches[0].atoms.matched_pairs == ((0, 0), (1, 1))


def test_molecule_substructure_matches_overlay():
    pattern = Molecule.parse(
        (
            '{:atoms ["N" "B"] :bonds [] '
            ':dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}'
        ),
    )
    host = Molecule.parse(
        (
            '{:atoms ["N" "B" "C"] :bonds [] '
            ':dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}'
        ),
    )

    matches = pattern.substructure_matches(
        host,
        config=SubstructureSearchConfig(
            match_algorithm=SubstructureMatchAlgorithm.Incidence(),
            subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Glasgow(),
        ),
    )

    assert len(matches) == 1
    match = matches[0]
    assert match.atoms.matched_pairs == ((0, 0), (1, 1))
    assert match.bonds.matched_pairs == ()
    assert match.dative_bonds.matched_pairs == ((0, 0),)
    assert match.aromatic_systems.matched_pairs == ()
    assert match.multicenter_bonds.matched_pairs == ()
    assert match.noncovalent_bonds.matched_pairs == ()
    assert match.stereo_atoms.matched_pairs == ()
    assert match.stereo_bonds.matched_pairs == ()


def test_molecule_substructure_matches_empty():
    pattern = Molecule.parse('{:atoms ["O"] :bonds []}')
    host = Molecule.parse('{:atoms ["C"] :bonds []}')

    assert pattern.substructure_matches(
        host,
        config=SubstructureSearchConfig(
            match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(),
            subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.BitsetRi(),
        ),
    ) == []


def test_molecule_substructure_matches_keyword_error():
    pattern = Molecule.parse('{:atoms ["C"] :bonds []}')
    host = Molecule.parse('{:atoms ["C"] :bonds []}')

    with pytest.raises(
        TypeError,
        match=(
            "^Molecule.substructure_matches\\(\\) takes 1 positional arguments "
            "but 2 were given$"
        ),
    ):
        pattern.substructure_matches(host, SubstructureSearchConfig())


@pytest.mark.parametrize(
    "method",
    [
        "has_substructure_match",
        "first_substructure_match",
        "count_substructure_matches",
    ],
)
def test_molecule_substructure_screening_keyword_error(method):
    pattern = Molecule.parse('{:atoms ["C"] :bonds []}')
    host = Molecule.parse('{:atoms ["C"] :bonds []}')

    with pytest.raises(TypeError):
        getattr(pattern, method)(host, SubstructureSearchConfig())


@pytest.mark.parametrize(
    ("reduce_pattern_symmetry", "reduce_host_symmetry", "expected_count"),
    [
        (False, False, 4),
        (True, False, 2),
        (False, True, 2),
        (True, True, 1),
    ],
)
def test_molecule_substructure_screening_and_symmetry(
    reduce_pattern_symmetry,
    reduce_host_symmetry,
    expected_count,
):
    butane = Molecule.parse(
        '{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]}'
    )
    chain3 = Molecule.parse('{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}')
    absent = Molecule.parse('{:atoms ["S"] :bonds []}')
    config = SubstructureSearchConfig(
        reduce_pattern_symmetry=reduce_pattern_symmetry,
        reduce_host_symmetry=reduce_host_symmetry,
    )

    matches = chain3.substructure_matches(butane, config=config)
    assert len(matches) == expected_count
    assert chain3.has_substructure_match(butane, config=config) is True
    assert chain3.first_substructure_match(butane, config=config) == matches[0]
    assert chain3.count_substructure_matches(butane, config=config) == expected_count

    assert absent.has_substructure_match(butane, config=config) is False
    assert absent.first_substructure_match(butane, config=config) is None
    assert absent.count_substructure_matches(butane, config=config) == 0


def test_prepared_substructure_host_contract():
    source = Molecule.parse('{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')
    expected = source.copy()
    prepared_host = PreparedSubstructureHost(source)

    # Mutating source after construction does not affect the prepared host.
    source.atoms[0].charge = 3
    snapshot = prepared_host.molecule()
    assert snapshot == expected
    assert snapshot is not source

    # Mutating the returned detached snapshot does not affect the prepared host.
    snapshot.atoms[0].charge = -2
    assert prepared_host.molecule() == expected

    assert repr(prepared_host) == "PreparedSubstructureHost(atoms=3, bonds=2)"
    with pytest.raises(AttributeError):
        prepared_host.extra = 1
    assert prepared_host == prepared_host
    assert prepared_host != PreparedSubstructureHost(expected)


@pytest.mark.parametrize(
    "strategy",
    [
        SubstructureMatchAlgorithm.GraphAndOverlays(),
        SubstructureMatchAlgorithm.Incidence(),
    ],
)
@pytest.mark.parametrize(
    "subiso",
    [
        SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
        SubgraphIsomorphismAlgorithm.RayKirsch(),
        SubgraphIsomorphismAlgorithm.BitsetRi(),
        SubgraphIsomorphismAlgorithm.Glasgow(),
    ],
)
@pytest.mark.parametrize(
    ("reduce_pattern_symmetry", "reduce_host_symmetry", "expected_count"),
    [
        (False, False, 12),
        (True, False, 6),
        (False, True, 1),
        (True, True, 1),
    ],
)
def test_prepared_substructure_pattern_parity(
    strategy,
    subiso,
    reduce_pattern_symmetry,
    reduce_host_symmetry,
    expected_count,
):
    cyclohexane = Molecule.parse(
        '{:atoms ["C" "C" "C" "C" "C" "C"] '
        ':bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}'
    )
    pattern = Molecule.parse('{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}')
    config = SubstructureSearchConfig(
        match_algorithm=strategy,
        subgraph_isomorphism_algorithm=subiso,
        reduce_pattern_symmetry=reduce_pattern_symmetry,
        reduce_host_symmetry=reduce_host_symmetry,
    )
    prepared_host = PreparedSubstructureHost(cyclohexane)
    prepared_pattern = PreparedSubstructurePattern(pattern, config=config)

    direct = pattern.substructure_matches(cyclohexane, config=config)
    from_mol = prepared_pattern.matches(cyclohexane)
    from_host = prepared_pattern.matches(prepared_host)

    assert len(direct) == expected_count
    assert from_mol == direct
    assert from_host == direct
    assert prepared_pattern.has_match(cyclohexane) is True
    assert prepared_pattern.has_match(prepared_host) is True
    assert prepared_pattern.first_match(cyclohexane) == direct[0]
    assert prepared_pattern.first_match(prepared_host) == direct[0]
    assert prepared_pattern.count_matches(cyclohexane) == expected_count
    assert prepared_pattern.count_matches(prepared_host) == expected_count


def test_prepared_substructure_pattern_contract_and_errors():
    source = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')
    expected = source.copy()
    prepared = PreparedSubstructurePattern(source)

    source.atoms[0].charge = 4
    snapshot = prepared.pattern()
    assert snapshot == expected
    snapshot.atoms[0].charge = -1
    assert prepared.pattern() == expected

    assert prepared.config == SubstructureSearchConfig.default()
    assert repr(prepared) == (
        f"PreparedSubstructurePattern(atoms=2, bonds=1, config={prepared.config!r})"
    )

    with pytest.raises(AttributeError):
        prepared.config = SubstructureSearchConfig()
    assert prepared == prepared
    assert prepared != PreparedSubstructurePattern(expected)
    with pytest.raises(TypeError):
        PreparedSubstructurePattern(source, SubstructureSearchConfig())
    with pytest.raises(TypeError):
        prepared.matches("not a host")
    with pytest.raises(TypeError):
        prepared.has_match(123)
    with pytest.raises(TypeError):
        prepared.first_match(object())
    with pytest.raises(TypeError):
        prepared.count_matches(None)
