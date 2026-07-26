import pytest

from umol import (
    AutomorphismAlgorithm,
    CommonSubgraphEnumerationAlgorithm,
    ConnectedComponentsAlgorithm,
    MaximumIndependentSetAlgorithm,
    RelevantCycleEnumerationAlgorithm,
    SimpleCycleEnumerationAlgorithm,
    SubgraphEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
)


@pytest.mark.parametrize(
    "algorithm,equal,expected_repr",
    [
        pytest.param(
            AutomorphismAlgorithm.Nauty(),
            AutomorphismAlgorithm.Nauty(),
            "AutomorphismAlgorithm.Nauty()",
            id="automorphism-nauty",
        ),
        pytest.param(
            ConnectedComponentsAlgorithm.Bfs(),
            ConnectedComponentsAlgorithm.Bfs(),
            "ConnectedComponentsAlgorithm.Bfs()",
            id="connected-components-bfs",
        ),
        pytest.param(
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
            "SimpleCycleEnumerationAlgorithm.ReadTarjan()",
            id="simple-cycle-enumeration-read-tarjan",
        ),
        pytest.param(
            RelevantCycleEnumerationAlgorithm.Vismara(),
            RelevantCycleEnumerationAlgorithm.Vismara(),
            "RelevantCycleEnumerationAlgorithm.Vismara()",
            id="relevant-cycle-enumeration-vismara",
        ),
        pytest.param(
            MaximumIndependentSetAlgorithm.BranchAndBound(),
            MaximumIndependentSetAlgorithm.BranchAndBound(),
            "MaximumIndependentSetAlgorithm.BranchAndBound()",
            id="maximum-independent-set-branch-and-bound",
        ),
        pytest.param(
            SubgraphEnumerationAlgorithm.Esu(),
            SubgraphEnumerationAlgorithm.Esu(),
            "SubgraphEnumerationAlgorithm.Esu()",
            id="subgraph-enumeration-esu",
        ),
    ],
)
def test_algorithm_value(algorithm, equal, expected_repr):
    assert algorithm == equal
    assert repr(algorithm) == expected_repr


@pytest.mark.parametrize(
    "algorithm,equal,unequal,expected_repr",
    [
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
            CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
            "CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking()",
            id="modular-product-backtracking",
        ),
        pytest.param(
            CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
            CommonSubgraphEnumerationAlgorithm.DirectBacktracking(),
            CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking(),
            "CommonSubgraphEnumerationAlgorithm.DirectBacktracking()",
            id="direct-backtracking",
        ),
    ],
)
def test_common_subgraph_enumeration_algorithm_value(
    algorithm, equal, unequal, expected_repr
):
    assert algorithm == equal
    assert algorithm != unequal
    assert repr(algorithm) == expected_repr


@pytest.mark.parametrize(
    "algorithm,equal,unequal,expected_repr",
    [
        pytest.param(
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            SubstructureMatchAlgorithm.Incidence(),
            "SubstructureMatchAlgorithm.GraphAndOverlays()",
            id="graph-and-overlays",
        ),
        pytest.param(
            SubstructureMatchAlgorithm.Incidence(),
            SubstructureMatchAlgorithm.Incidence(),
            SubstructureMatchAlgorithm.GraphAndOverlays(),
            "SubstructureMatchAlgorithm.Incidence()",
            id="incidence",
        ),
    ],
)
def test_substructure_match_algorithm_value(algorithm, equal, unequal, expected_repr):
    assert algorithm == equal
    assert algorithm != unequal
    assert repr(algorithm) == expected_repr


@pytest.mark.parametrize(
    "algorithm,equal,expected_repr",
    [
        pytest.param(
            SubgraphIsomorphismAlgorithm.Vf2(),
            SubgraphIsomorphismAlgorithm.Vf2(),
            "SubgraphIsomorphismAlgorithm.Vf2()",
            id="vf2",
        ),
        pytest.param(
            SubgraphIsomorphismAlgorithm.Ullmann(),
            SubgraphIsomorphismAlgorithm.Ullmann(),
            "SubgraphIsomorphismAlgorithm.Ullmann()",
            id="ullmann",
        ),
        pytest.param(
            SubgraphIsomorphismAlgorithm.Ri(),
            SubgraphIsomorphismAlgorithm.Ri(),
            "SubgraphIsomorphismAlgorithm.Ri()",
            id="ri",
        ),
        pytest.param(
            SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6),
            SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6),
            "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6)",
            id="arc-match",
        ),
        pytest.param(
            SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
            SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
            "SubgraphIsomorphismAlgorithm.Vf2Rdkit()",
            id="vf2-rdkit",
        ),
        pytest.param(
            SubgraphIsomorphismAlgorithm.RayKirsch(),
            SubgraphIsomorphismAlgorithm.RayKirsch(),
            "SubgraphIsomorphismAlgorithm.RayKirsch()",
            id="ray-kirsch",
        ),
    ],
)
def test_subgraph_isomorphism_algorithm_value(algorithm, equal, expected_repr):
    assert algorithm == equal
    assert repr(algorithm) == expected_repr
