import umol


def test_import_umol():
    assert isinstance(umol.__version__, str)
    assert umol.__version__ != "0.0.0"


def test_native_present():
    assert hasattr(umol, "_native")
