from umol import (
    CisTransStereo,
    CisTransStereoAst,
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
    assert StereoTerm.Lit(1) == StereoTerm.Lit(1)


def test_stereoterm_apply():
    term = StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))
    assert term == StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))


def test_stereocosetast_term():
    assert StereoCosetAst.Term(StereoTerm.Lit(2)) == StereoCosetAst.Term(StereoTerm.Lit(2))


def test_tetrahedralstereoast_stereo():
    assert TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1)) == TetrahedralStereoAst.Stereo(
        StereoCosetAst.Lit(1)
    )


def test_tetrahedralstereoast_not_stereo():
    assert TetrahedralStereoAst.NotStereo() == TetrahedralStereoAst.NotStereo()


def test_stereocosetast_eq_hash_repr():
    assert StereoCosetAst.Lit(1) == StereoCosetAst.Lit(1)
    assert StereoCosetAst.Lit(1) != StereoCosetAst.Lit(0)
    assert len({StereoCosetAst.Lit(1), StereoCosetAst.Lit(1)}) == 1
    assert repr(StereoCosetAst.Lit(1)) == "StereoCosetAst.Lit(1)"


def test_stereoterm_eq_repr():
    assert StereoTerm.Lit(0) == StereoTerm.Lit(0)
    assert StereoTerm.Lit(0) != StereoTerm.Swap(StereoTerm.Lit(0))
    assert repr(StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))) == (
        "StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))"
    )


def test_tetrahedralstereoast_stereo_repr():
    assert repr(TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1))) == (
        "TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1))"
    )


def test_cistransstereoast_stereo():
    assert CisTransStereoAst.Stereo(StereoCosetAst.Lit(1)) == CisTransStereoAst.Stereo(
        StereoCosetAst.Lit(1)
    )


def test_cistransstereoast_not_stereo():
    assert CisTransStereoAst.NotStereo() == CisTransStereoAst.NotStereo()


def test_cistransstereoast_stereo_repr():
    assert repr(CisTransStereoAst.Stereo(StereoCosetAst.Lit(1))) == (
        "CisTransStereoAst.Stereo(StereoCosetAst.Lit(1))"
    )


def test_cistransstereo_enum():
    assert CisTransStereo.Z == CisTransStereo.Z
    assert CisTransStereo.Z != CisTransStereo.E
    assert len({CisTransStereo.E, CisTransStereo.E, CisTransStereo.Z}) == 2
