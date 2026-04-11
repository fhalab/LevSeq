from levseq import seqfit


def test_seqfit_module_imports():
    assert callable(seqfit.gen_seqfitvis)
