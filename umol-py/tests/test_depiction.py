import xml.etree.ElementTree as ET

import pytest

from umol import _native as _native_module

if not hasattr(_native_module, "Svg"):
    pytest.skip("umol-py was built without the depiction feature", allow_module_level=True)

from umol import (  # noqa: E402
    BondDelta,
    BondFieldChange,
    ContradictionError,
    Delta,
    DepictConfig,
    Depiction,
    Deltas,
    Molecule,
    MoleculeLayoutAlgorithm,
    NumForm,
    Reaction,
    Svg,
)


def test_depict_config_new():
    config = DepictConfig()

    assert config == DepictConfig.default()
    assert config.layout_algorithm == MoleculeLayoutAlgorithm.CoordGen()
    assert repr(config) == "DepictConfig.default()"


def test_depiction_constructor_error():
    with pytest.raises(TypeError):
        Depiction()


def test_molecule_layout_algorithm():
    algorithm = MoleculeLayoutAlgorithm.CoordGen()

    assert algorithm == MoleculeLayoutAlgorithm.CoordGen()
    assert repr(algorithm) == "MoleculeLayoutAlgorithm.CoordGen()"


def test_svg_constructor_error():
    with pytest.raises(TypeError):
        Svg()


@pytest.mark.parametrize(
    ("source", "item_kind", "reference", "svg_class"),
    [
        pytest.param(
            '{:atoms ["C#i13#c+#h2#u2" "O#c-"] :bonds [[0 1 "2"]]}',
            "atom",
            "molecule/atom/0",
            "umol-atom-right-subscript",
            id="atom-labels",
        ),
        pytest.param(
            """{:atoms ["C" "F" "Cl" "Br" "I"]
                 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                 :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}""",
            "wedge",
            "molecule/stereo-atom/0",
            "umol-wedge-solid",
            id="tetrahedral-stereo",
        ),
        pytest.param(
            """{:atoms ["C" "C" "C" "C" "C" "C"]
                 :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                         [3 4 "1"] [4 5 "1"] [5 0 "1"]]
                 :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*"}]}""",
            "dashed-contour",
            "molecule/aromatic-system/0",
            "umol-dashed-contour",
            id="aromatic-system",
        ),
    ],
)
def test_molecule_depict_with(source, item_kind, reference, svg_class):
    svg = Molecule.parse(source).depict_with(MoleculeLayoutAlgorithm.CoordGen())
    text = svg._repr_svg_()
    root = ET.fromstring(text)
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))

    assert isinstance(svg, Svg)
    assert root.attrib["class"] == "umol-depiction"
    assert any(
        group.attrib.get("data-umol-item") == item_kind
        and reference in group.attrib.get("data-umol-references", "").split()
        for group in groups
    )
    assert f'class="{svg_class}"' in text
    assert 'data-umol-item="marker"' not in text


def test_molecule_depict_with_algorithm_error():
    with pytest.raises(TypeError):
        Molecule().depict_with()


def test_reaction_depict_with():
    svg = Reaction.parse(
        """{:lhs {:atoms ["C" "O"]
                   :bonds [{:id :co :atoms [0 1] :attrs "1"}]}
             :deltas [{:bond {:modify [:co "2"]}}
                      {:atom {:add [:n "N"]}}
                      {:bond {:add [0 :n "1"]}}]}"""
    ).depict_with(MoleculeLayoutAlgorithm.CoordGen())
    root = ET.fromstring(svg._repr_svg_())
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))
    mapping_groups = [
        group
        for group in groups
        if "correspondence-pair/" in group.attrib.get("data-umol-references", "")
    ]
    arrow = next(group for group in groups if group.attrib.get("data-umol-item") == "arrow")
    shaft, head = list(arrow)

    assert isinstance(svg, Svg)
    assert [group.attrib["data-umol-references"] for group in mapping_groups] == [
        "reaction-lhs/atom/0 correspondence-pair/0",
        "reaction-lhs/atom/1 correspondence-pair/1",
        "reaction-rhs/atom/0 correspondence-pair/0",
        "reaction-rhs/atom/1 correspondence-pair/1",
    ]
    assert [list(group)[0].attrib["font-size"] for group in mapping_groups] == [
        "0.3825",
        "0.3825",
        "0.3825",
        "0.3825",
    ]
    assert shaft.attrib["class"] == "umol-arrow-shaft"
    assert shaft.attrib["x1"] == "-0.75"
    assert shaft.attrib["x2"] == "0.51"
    assert head.attrib == {
        "class": "umol-arrow-head",
        "points": "0.75,0 0.51,-0.11 0.51,0.11",
        "fill": "currentColor",
    }


def test_reaction_depict_with_error():
    lhs = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')
    deltas = Deltas(
        [
            Delta.Bond(
                BondDelta.ModifyField(
                    id=0,
                    change=BondFieldChange.Order(
                        old=NumForm.Lit(2),
                        new=NumForm.Lit(3),
                    ),
                )
            )
        ]
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        Reaction(lhs, deltas).depict_with(MoleculeLayoutAlgorithm.CoordGen())
