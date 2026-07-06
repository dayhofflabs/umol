from umol import (
    Permutation,
    StereoCosetAst,
    StereoTerm,
    TetrahedralStereoAst,
)


def test_permutation_image_degree():
    p = Permutation([1, 0, 2, 3])
    assert p.image() == [1, 0, 2, 3]
    assert p.degree == 4


def test_permutation_identity():
    assert Permutation.identity(4).image() == [0, 1, 2, 3]


def test_permutation_eq_hash():
    assert Permutation([1, 0, 2]) == Permutation([1, 0, 2])
    assert Permutation([1, 0, 2]) != Permutation([0, 1, 2])
    assert len({Permutation([1, 0, 2]), Permutation([1, 0, 2])}) == 1


def test_stereoterm_lit():
    match StereoTerm.Lit(1):
        case StereoTerm.Lit(n):
            assert n == 1
        case _:
            raise AssertionError


def test_stereoterm_apply():
    term = StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))
    match term:
        case StereoTerm.Apply(inner, perm):
            match inner:
                case StereoTerm.Lit(n):
                    assert n == 0
                case _:
                    raise AssertionError
            assert perm.image() == [1, 0, 2, 3]
        case _:
            raise AssertionError


def test_stereocosetast_term():
    match StereoCosetAst.Term(StereoTerm.Lit(2)):
        case StereoCosetAst.Term(inner):
            match inner:
                case StereoTerm.Lit(n):
                    assert n == 2
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_tetrahedralstereoast_stereo():
    match TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1)):
        case TetrahedralStereoAst.Stereo(coset):
            match coset:
                case StereoCosetAst.Lit(n):
                    assert n == 1
                case _:
                    raise AssertionError
        case _:
            raise AssertionError


def test_tetrahedralstereoast_not_stereo():
    match TetrahedralStereoAst.NotStereo():
        case TetrahedralStereoAst.NotStereo():
            pass
        case _:
            raise AssertionError
