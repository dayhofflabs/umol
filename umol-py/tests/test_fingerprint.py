import pytest

from umol import EcfpHashScheme, RefinementRounds, WlHashScheme


@pytest.mark.parametrize("rounds", [0, 3])
def test_refinement_rounds_fixed(rounds):
    value = RefinementRounds.Fixed(rounds=rounds)

    assert value.rounds == rounds
    assert value == RefinementRounds.Fixed(rounds=rounds)
    assert repr(value) == f"RefinementRounds.Fixed(rounds={rounds})"


def test_refinement_rounds_to_fixpoint():
    value = RefinementRounds.ToFixpoint()

    assert value == RefinementRounds.ToFixpoint()
    assert value != RefinementRounds.Fixed(rounds=0)
    assert repr(value) == "RefinementRounds.ToFixpoint()"
    with pytest.raises(AttributeError):
        value.rounds


def test_wl_hash_scheme():
    value = WlHashScheme.Xxh3SortedWidth64V1()

    assert value == WlHashScheme.Xxh3SortedWidth64V1()
    assert repr(value) == "WlHashScheme.Xxh3SortedWidth64V1()"


def test_ecfp_hash_scheme():
    value = EcfpHashScheme.Xxh3Width64V1()

    assert value == EcfpHashScheme.Xxh3Width64V1()
    assert repr(value) == "EcfpHashScheme.Xxh3Width64V1()"
