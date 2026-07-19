import re

import pytest

from umol import AtomAst, AtomTypeRegistry, Element, ValenceEntry, ValenceTable


def test_atom_type_registry_default():
    registry = AtomTypeRegistry.default()

    assert registry == AtomTypeRegistry.default()
    assert repr(registry) == "AtomTypeRegistry.default()"
    assert registry.content_hash_hex == f"{registry.content_hash:016x}"


def test_atom_type_registry_from_atoms():
    atom = AtomAst.parse("C#c0#v4")
    registry = AtomTypeRegistry.from_atoms(
        [atom, AtomAst.parse("C#c+#v3"), AtomAst.parse("O#c0#v2")]
    )

    assert registry == AtomTypeRegistry.from_atoms(
        [
            AtomAst.parse("C#c0#v4"),
            AtomAst.parse("C#c+#v3"),
            AtomAst.parse("O#c0#v2"),
        ]
    )
    assert registry != AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4")])

    atom.charge = 1
    assert registry.patterns_for_element_and_charge(Element("C"), 0) == [
        AtomAst.parse("C#c0#v4")
    ]


@pytest.mark.parametrize(
    ("atom", "message"),
    [
        (
            AtomAst.parse("*#c0"),
            "atom type registry entry 0 must have a literal element",
        ),
        (
            AtomAst.parse("C"),
            "atom type registry entry 0 must have a literal charge",
        ),
        (
            AtomAst.parse("C#c128"),
            "atom type registry entry 0 charge 128 is outside -128..=127",
        ),
    ],
)
def test_atom_type_registry_from_atoms_error(atom, message):
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        AtomTypeRegistry.from_atoms([atom])


def test_atom_type_registry_from_toml():
    registry = AtomTypeRegistry.from_toml(
        '[C]\n0 = ["C#c0#v4"]\n1 = ["C#c+#v3"]\n[O]\n0 = ["O#c0#v2"]'
    )

    assert registry.patterns_for_element(Element("C")) == [
        AtomAst.parse("C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!"),
        AtomAst.parse("C#i=#c+#h0#n0#u0#s#v3#d0#t0#a!#m!"),
    ]
    assert registry.patterns_for_element_and_charge(Element("O"), 0) == [
        AtomAst.parse("O#i=#c0#h0#n0#u0#s#v2#d0#t0#a!#m!")
    ]


def test_atom_type_registry_from_toml_error():
    with pytest.raises(
        ValueError,
        match="^invalid atom type registry: unknown element: X$",
    ):
        AtomTypeRegistry.from_toml('[X]\n0 = ["X#c0"]')


def test_atom_type_registry_patterns():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("C#c+#v3")]
    )
    patterns = registry.patterns_for_element(Element("C"))

    patterns[0].charge = -1
    patterns.pop()

    assert registry.patterns_for_element(Element("C")) == [
        AtomAst.parse("C#c0#v4"),
        AtomAst.parse("C#c+#v3"),
    ]


def test_atom_type_registry_repr():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("O#c0#v2")]
    )

    assert repr(registry) == (
        'AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4"), '
        'AtomAst.parse("O#c0#v2")])'
    )


def test_atom_type_registry_mutation():
    registry = AtomTypeRegistry.default()

    with pytest.raises(AttributeError):
        registry.content_hash = 0
    assert not hasattr(registry, "add")


@pytest.mark.parametrize(
    ("target_covalences", "aromatic_valences", "expected"),
    [
        ([], [], ValenceEntry()),
        (
            [6, 2, 4, 2],
            [4, 2],
            ValenceEntry(
                target_covalences=[2, 2, 4, 6], aromatic_valences=[2, 4]
            ),
        ),
    ],
)
def test_valence_entry_new(target_covalences, aromatic_valences, expected):
    assert (
        ValenceEntry(
            target_covalences=target_covalences,
            aromatic_valences=aromatic_valences,
        )
        == expected
    )


def test_valence_entry_properties():
    entry = ValenceEntry(target_covalences=[6, 2, 4], aromatic_valences=[4, 2])

    target_covalences = entry.target_covalences
    aromatic_valences = entry.aromatic_valences
    target_covalences.append(8)
    aromatic_valences.clear()

    assert entry.target_covalences == [2, 4, 6]
    assert entry.aromatic_valences == [2, 4]


@pytest.mark.parametrize(
    ("entry", "expected"),
    [
        (
            ValenceEntry(),
            "ValenceEntry(target_covalences=[], aromatic_valences=[])",
        ),
        (
            ValenceEntry(target_covalences=[6, 2, 4], aromatic_valences=[4, 2]),
            "ValenceEntry(target_covalences=[2, 4, 6], aromatic_valences=[2, 4])",
        ),
    ],
)
def test_valence_entry_repr(entry, expected):
    assert repr(entry) == expected


def test_valence_entry_mutation():
    entry = ValenceEntry(target_covalences=[2, 4], aromatic_valences=[2])

    with pytest.raises(AttributeError):
        entry.target_covalences = [6]


def test_valence_table_new():
    entries = {
        Element("C"): ValenceEntry(
            target_covalences=[6, 2, 4], aromatic_valences=[3, 2]
        ),
        Element("O"): ValenceEntry(target_covalences=[2]),
    }
    table = ValenceTable(entries=entries)

    entries.clear()

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4, 6], aromatic_valences=[2, 3]
    )
    assert table.entry(Element("O")) == ValenceEntry(target_covalences=[2])
    assert table.entry(Element("N")) is None


def test_valence_table_default():
    table = ValenceTable.default()

    assert table == ValenceTable.default()
    assert repr(table) == "ValenceTable.default()"
    assert table.content_hash_hex == f"{table.content_hash:016x}"


def test_valence_table_from_toml():
    table = ValenceTable.from_toml(
        "[C]\ntarget_covalences = [6, 2, 4]\naromatic_valences = [3, 2]"
    )

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4, 6], aromatic_valences=[2, 3]
    )
    assert table.entry(Element("O")) is None


def test_valence_table_from_toml_error():
    with pytest.raises(ValueError, match="^invalid valence table: unknown element: X$"):
        ValenceTable.from_toml("[X]\ntarget_covalences = [1]")


def test_valence_table_entry():
    table = ValenceTable(
        entries={
            Element("C"): ValenceEntry(
                target_covalences=[2, 4], aromatic_valences=[2]
            )
        }
    )
    entry = table.entry(Element("C"))
    target_covalences = entry.target_covalences
    aromatic_valences = entry.aromatic_valences

    target_covalences.append(6)
    aromatic_valences.clear()

    assert table.entry(Element("C")) == ValenceEntry(
        target_covalences=[2, 4], aromatic_valences=[2]
    )


@pytest.mark.parametrize(
    ("table", "expected"),
    [
        (ValenceTable(entries={}), "ValenceTable(entries={})"),
        (
            ValenceTable(
                entries={
                    Element("C"): ValenceEntry(target_covalences=[4, 2]),
                    Element("O"): ValenceEntry(target_covalences=[2]),
                }
            ),
            "ValenceTable(entries={Element('C'): "
            "ValenceEntry(target_covalences=[2, 4], aromatic_valences=[]), "
            "Element('O'): ValenceEntry(target_covalences=[2], aromatic_valences=[])})",
        ),
    ],
)
def test_valence_table_repr(table, expected):
    assert repr(table) == expected


def test_valence_table_mutation():
    table = ValenceTable(entries={})

    with pytest.raises(AttributeError):
        table.content_hash = 0
    assert not hasattr(table, "insert")
