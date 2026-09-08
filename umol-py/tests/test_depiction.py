import inspect
import xml.etree.ElementTree as ET

import pytest

from umol import _native as _native_module

if not hasattr(_native_module, "Depiction"):
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
    MoleculeLayout,
    MoleculeLayoutAlgorithm,
    NumForm,
    Reaction,
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


@pytest.mark.parametrize(
    ("callable_", "expected"),
    [
        (MoleculeLayout, "(positions)"),
        (Molecule.layout, "(self, /, *, algorithm=Ellipsis)"),
        (Molecule.depict_with_layout, "(self, /, layout)"),
    ],
)
def test_layout_signature(callable_, expected):
    assert str(inspect.signature(callable_)) == expected


def test_molecule_layout_new():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    assert layout.positions == [(0.0, 1.0), (2.0, -1.0)]
    assert len(layout) == 2
    assert layout == MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])
    assert layout != MoleculeLayout([(0.0, 1.0), (2.0, 1.0)])


@pytest.mark.parametrize(
    ("positions", "message"),
    [
        pytest.param(
            [(0.0, 0.0), (float("nan"), 1.0)],
            "atom 1 has non-finite position",
            id="nan",
        ),
        pytest.param(
            [(float("inf"), 0.0)],
            "atom 0 has non-finite position",
            id="infinity",
        ),
    ],
)
def test_molecule_layout_new_error(positions, message):
    with pytest.raises(ValueError, match=message):
        MoleculeLayout(positions)


def test_molecule_layout_position():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    assert layout.position(1) == (2.0, -1.0)
    with pytest.raises(IndexError, match="atom 2 is outside layout frame of size 2"):
        layout.position(2)


def test_molecule_layout_with_position():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    moved = layout.with_position(1, (3.5, 4.0))

    assert moved == MoleculeLayout([(0.0, 1.0), (3.5, 4.0)])
    assert layout == MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])


@pytest.mark.parametrize(
    ("atom_id", "position", "message"),
    [
        pytest.param(2, (0.0, 0.0), "atom 2 is outside layout frame of size 2", id="frame"),
        pytest.param(0, (float("nan"), 0.0), "atom 0 has non-finite position", id="nan"),
    ],
)
def test_molecule_layout_with_position_error(atom_id, position, message):
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    with pytest.raises(ValueError, match=message):
        layout.with_position(atom_id, position)


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
def test_molecule_depict(source, item_kind, reference, svg_class):
    depiction = Molecule.parse(source).depict()
    text = depiction.render_svg()
    root = ET.fromstring(text)
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))

    assert isinstance(depiction, Depiction)
    assert depiction._repr_svg_() == text
    assert root.attrib["class"] == "umol-depiction"
    assert any(
        group.attrib.get("data-umol-item") == item_kind
        and reference in group.attrib.get("data-umol-references", "").split()
        for group in groups
    )
    assert f'class="{svg_class}"' in text
    assert 'data-umol-item="marker"' not in text


def test_molecule_depict_with():
    molecule = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "2"]]}')

    assert (
        molecule.depict_with(DepictConfig()).render_svg()
        == molecule.depict().render_svg()
    )


def test_molecule_depict_write_svg(tmp_path):
    depiction = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "2"]]}').depict()
    output = tmp_path / "molecule.svg"

    output.write_text(depiction.render_svg())

    assert output.read_text() == depiction.render_svg()


def test_molecule_depict_with_config_error():
    with pytest.raises(TypeError):
        Molecule().depict_with()


def test_molecule_layout():
    molecule = Molecule.parse('{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')

    layout = molecule.layout()

    assert isinstance(layout, MoleculeLayout)
    assert len(layout) == 3
    assert layout == molecule.layout(algorithm=MoleculeLayoutAlgorithm.CoordGen())
    assert all(isinstance(coordinate, float) for x, y in layout.positions for coordinate in (x, y))


def test_molecule_layout_keyword_only_error():
    with pytest.raises(TypeError):
        Molecule().layout(MoleculeLayoutAlgorithm.CoordGen())


def test_molecule_depict_with_layout():
    molecule = Molecule.parse('{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')

    assert (
        molecule.depict_with_layout(molecule.layout()).render_svg()
        == molecule.depict().render_svg()
    )


def bond_lines(text):
    root = ET.fromstring(text)
    lines = {}
    for group in root.findall("{http://www.w3.org/2000/svg}g"):
        if group.attrib.get("data-umol-item") != "bond":
            continue
        (line,) = list(group)
        lines[group.attrib["data-umol-references"]] = {
            name: float(line.attrib[name]) for name in ("x1", "y1", "x2", "y2")
        }
    return lines


def test_molecule_depict_with_layout_roundtrip():
    molecule = Molecule.parse(
        '{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]}'
    )
    layout = molecule.layout()
    before = bond_lines(molecule.depict_with_layout(layout).render_svg())

    moved = layout.with_position(1, (2.5, -1.25))
    after = bond_lines(molecule.depict_with_layout(moved).render_svg())

    assert (after["molecule/bond/0"]["x2"], after["molecule/bond/0"]["y2"]) == (2.5, 1.25)
    assert (after["molecule/bond/1"]["x1"], after["molecule/bond/1"]["y1"]) == (2.5, 1.25)
    assert (after["molecule/bond/0"]["x1"], after["molecule/bond/0"]["y1"]) == (
        before["molecule/bond/0"]["x1"],
        before["molecule/bond/0"]["y1"],
    )
    assert (after["molecule/bond/1"]["x2"], after["molecule/bond/1"]["y2"]) == (
        before["molecule/bond/1"]["x2"],
        before["molecule/bond/1"]["y2"],
    )
    assert after["molecule/bond/2"] == before["molecule/bond/2"]


def test_molecule_depict_with_layout_error():
    molecule = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')

    with pytest.raises(
        ValueError,
        match="^layout frame: molecule atom count 2 does not match layout atom count 1$",
    ):
        molecule.depict_with_layout(MoleculeLayout([(0.0, 0.0)]))


def test_reaction_depict():
    depiction = Reaction.parse(
        """{:lhs {:atoms ["C" "O"]
                   :bonds [{:id :co :atoms [0 1] :attrs "1"}]}
             :deltas [{:bond {:modify [:co "2"]}}
                      {:atom {:add [:n "N"]}}
                      {:bond {:add [0 :n "1"]}}]}"""
    ).depict()
    text = depiction.render_svg()
    root = ET.fromstring(text)
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))
    mapping_groups = [
        group
        for group in groups
        if "correspondence-pair/" in group.attrib.get("data-umol-references", "")
    ]
    arrow = next(group for group in groups if group.attrib.get("data-umol-item") == "arrow")
    shaft, head = list(arrow)

    assert isinstance(depiction, Depiction)
    assert depiction._repr_svg_() == text
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


def test_reaction_depict_with():
    reaction = Reaction.parse(
        """{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"]]}
             :deltas []}"""
    )

    assert (
        reaction.depict_with(DepictConfig()).render_svg()
        == reaction.depict().render_svg()
    )


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
        Reaction(lhs, deltas).depict_with(DepictConfig())
