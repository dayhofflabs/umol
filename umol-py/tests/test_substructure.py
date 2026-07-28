import pytest

from umol import (
    MoleculeAst,
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
        == SubgraphIsomorphismAlgorithm.Vf2Rdkit()
    )
    assert repr(config) == (
        "SubstructureSearchConfig("
        "match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(), "
        "subgraph_isomorphism_algorithm="
        "SubgraphIsomorphismAlgorithm.Vf2Rdkit())"
    )


@pytest.mark.parametrize(
    ("match_algorithm", "subgraph_isomorphism_algorithm"),
    [
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm.Vf2(),
            id="graph-vf2",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.Incidence(),
            SubgraphIsomorphismAlgorithm.Ullmann(),
            id="incidence-ullmann",
        ),
    ],
)
def test_substructure_search_config_value(
    match_algorithm,
    subgraph_isomorphism_algorithm,
):
    config = SubstructureSearchConfig(
        match_algorithm=match_algorithm,
        subgraph_isomorphism_algorithm=subgraph_isomorphism_algorithm,
    )

    assert config.match_algorithm == match_algorithm
    assert config.subgraph_isomorphism_algorithm == subgraph_isomorphism_algorithm
    assert repr(config) == (
        "SubstructureSearchConfig("
        f"match_algorithm={match_algorithm!r}, "
        f"subgraph_isomorphism_algorithm={subgraph_isomorphism_algorithm!r})"
    )


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
            SubgraphIsomorphismAlgorithm.Ullmann(),
        )


def test_molecule_ast_substructure_matches():
    pattern_source = '{:atoms ["C" "C"] :bonds [[0 1 "1"]]}'
    host_source = '{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}'
    pattern = MoleculeAst.parse(pattern_source)
    host = MoleculeAst.parse(host_source)
    pattern_before = MoleculeAst.parse(pattern_source)
    host_before = MoleculeAst.parse(host_source)

    matches = pattern.substructure_matches(host)

    assert isinstance(matches, list)
    assert [match.atoms.matched_pairs for match in matches] == [
        [(0, 0), (1, 1)],
        [(0, 1), (1, 0)],
    ]
    assert [match.bonds.matched_pairs for match in matches] == [[(0, 0)], [(0, 0)]]
    assert [match.dative_bonds.matched_pairs for match in matches] == [[], []]
    assert [match.aromatic_systems.matched_pairs for match in matches] == [[], []]
    assert [match.multicenter_bonds.matched_pairs for match in matches] == [[], []]
    assert [match.noncovalent_bonds.matched_pairs for match in matches] == [[], []]
    assert [match.stereo_atoms.matched_pairs for match in matches] == [[], []]
    assert [match.stereo_bonds.matched_pairs for match in matches] == [[], []]
    assert pattern == pattern_before
    assert host == host_before

    matches[0].atoms.matched_pairs.append((1, 2))
    assert matches[0].atoms.matched_pairs == [(0, 0), (1, 1)]


def test_molecule_ast_substructure_matches_overlay():
    pattern = MoleculeAst.parse(
        (
            '{:atoms ["N" "B"] :bonds [] '
            ':dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}'
        ),
    )
    host = MoleculeAst.parse(
        (
            '{:atoms ["N" "B" "C"] :bonds [] '
            ':dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}'
        ),
    )

    matches = pattern.substructure_matches(
        host,
        config=SubstructureSearchConfig(
            match_algorithm=SubstructureMatchAlgorithm.Incidence(),
            subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Ullmann(),
        ),
    )

    assert len(matches) == 1
    match = matches[0]
    assert match.atoms.matched_pairs == [(0, 0), (1, 1)]
    assert match.bonds.matched_pairs == []
    assert match.dative_bonds.matched_pairs == [(0, 0)]
    assert match.aromatic_systems.matched_pairs == []
    assert match.multicenter_bonds.matched_pairs == []
    assert match.noncovalent_bonds.matched_pairs == []
    assert match.stereo_atoms.matched_pairs == []
    assert match.stereo_bonds.matched_pairs == []


def test_molecule_ast_substructure_matches_empty():
    pattern = MoleculeAst.parse('{:atoms ["O"] :bonds []}')
    host = MoleculeAst.parse('{:atoms ["C"] :bonds []}')

    assert pattern.substructure_matches(
        host,
        config=SubstructureSearchConfig(
            match_algorithm=SubstructureMatchAlgorithm.GraphAndOverlays(),
            subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Vf2(),
        ),
    ) == []


def test_molecule_ast_substructure_matches_keyword_error():
    pattern = MoleculeAst.parse('{:atoms ["C"] :bonds []}')
    host = MoleculeAst.parse('{:atoms ["C"] :bonds []}')

    with pytest.raises(
        TypeError,
        match=(
            "^MoleculeAst.substructure_matches\\(\\) takes 1 positional arguments "
            "but 2 were given$"
        ),
    ):
        pattern.substructure_matches(host, SubstructureSearchConfig())
