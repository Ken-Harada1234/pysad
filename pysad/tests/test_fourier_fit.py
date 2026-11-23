"""MAIN の TRACK の FOURIER 解析 (pysad/fourier_fit.py = K64 の wjfit・wfres・wfbin・wiord)。
gs は MAIN の TRACK が SIGSEGV (tracka → tinip) で参照が無い。基準は K64 の Fortran 自体 (wjfit.f の tfstk に依る枝を外し、
wiord.f・wfres.f・wfbin.f と一緒に gfortran で組んだもの。wtune は素朴な DFT) の出力: data_fourier_ref.json
(文字は全行一致、係数 b は 1e-7 (|b| の下限 1e-6) の相対で一致、チューンは 1e-12)。"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, HERE)

from pysad import fourier_fit as ff

REF = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_fourier_ref.json')))


def _gen(lsp, np0, seed):
    n = np.arange(lsp)
    Z = []
    for i in range(np0):
        zz = np.zeros((lsp, 3), complex)
        a = [1.0 + 0.3 * i, 0.8, 0.05]
        nu0 = [0.2113 + 0.01 * i, 0.3071, 0.0123]
        for k in range(3):
            zz[:, k] = a[k] * np.exp(-2j * np.pi * nu0[k] * n)
        z0, z1, z2 = zz[:, 0].copy(), zz[:, 1].copy(), zz[:, 2].copy()
        zz[:, 0] = z0 * (1 + 0.05 * np.cos(2 * np.pi * 0.3071 * n) + 0.02 * np.conj(z1) * z1.real) + 0.01 * np.conj(z0) ** 2
        zz[:, 1] = z1 * (1 + 0.03 * np.cos(2 * np.pi * 0.2113 * n + 0.4))
        zz[:, 2] = z2 * (1 + 0.02 * np.cos(2 * np.pi * 0.2113 * n + 1.1))
        zz[:, 0] += 1e-3 * np.exp(-2j * np.pi * 0.07 * n)
        Z.append(zz)
    return np.array(Z)


def test_matches_k64_fortran():
    for key, ref in REF.items():
        lsp, np0, mo, nr, ad, rf, seed = ref['args']
        z = _gen(lsp, np0, seed)
        lines, recs = ff.fourier_report(z, np.zeros(np0, bool), mo, nr, bool(rf), ad)
        assert [l.rstrip() for l in lines] == [l.rstrip() for l in ref['text'][:len(lines)]], key
        for i, r in enumerate(recs):
            b = np.array(ref['b'][i])
            assert np.max(np.abs(r['b'] - b) / np.maximum(np.abs(b), 1e-6)) < 1e-7, key
            assert np.allclose(r['tune'], ref['tune'][i], atol=1e-12, rtol=0), key


def test_wiord_count_matches_lcoeff():
    for mo in range(1, 9):
        for rf in (False, True):
            nd = 3 if rf else 2
            io = ff.wiord(mo, rf)
            assert io.shape[1] == ff.lcoeff_of(mo, 0, rf)
            assert (io[:, :nd] == 0).all()


def test_pure_rotation_has_no_resonance_terms():
    lsp = 128
    n = np.arange(lsp)
    z = np.stack([1.5 * np.exp(-2j * np.pi * 27 / 128 * n), 0.7 * np.exp(-2j * np.pi * 41 / 128 * n), 0.0 * n + 0j], axis=1)
    lines, recs = ff.fourier_report(z[None], np.zeros(1, bool), 4, 0, False)
    r = recs[0]
    assert abs(r['b'][0] - 1.5 ** 2) < 1e-10 and abs(r['b'][1] - 0.7 ** 2) < 1e-10
    assert np.max(np.abs(r['b'][2:])) < 1e-10
    assert abs(r['tune'][0] - 27 / 128) < 2e-3 and abs(r['tune'][1] - 41 / 128) < 2e-3


def test_main_track_fourier_prints_analysis():
    import io
    from pysad.ffscmd import run_sad_string
    from pysad import track6 as t6
    s = open(os.path.join(HERE, 'bench_level27a_main.sad'), encoding='utf-8').read()
    head = s[:s.index('ON RFSW;')]
    body = chr(10).join(['ON RFSW;', 'OFF RAD;', 'EMIT USE=RING;', 'ON FOURIER;', 'MAXORDER=3;',
                      'TRACK USE=RING NPART=2 TURNS=63 NX=(1,3) NY=(1) NZ=(0);', 'STOP;', ''])
    t6.RNG.restore(17)
    buf = io.StringIO()
    run_sad_string(head + body, out=buf, cwd=HERE)
    out = buf.getvalue()
    assert '<Fourier analysis>' in out and 'Big 16 coefficients' in out and 'Tune =' in out
    assert out.count('Particle #') == 2
