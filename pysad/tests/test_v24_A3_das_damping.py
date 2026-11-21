"""v2.4.0: DynamicApertureSurvey の Damping->True・DampEnough (trackd の tpdamp)。参照は gs の ntloss (FODO リング 24 セル、6 GeV、RF 30 MV、
RAD ON の emit、追跡は NORAD、400 周、X = 1〜25 σ)。6 極なしでは損失は LOSSAMPL の振幅だけで決まり、gs と一致する (DampingRate は gs の値と 1e-11)。
6 極つき (K2 = 14・−20) の縁の点はカオス的で、gs と pysad で数周〜数百周違う点がある (縁の 4 点中 1〜3 点。損失の無い点・遠い点は一致)。"""
import io
import os
import re

from pysad.ffscmd import run_sad_string

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
FODO = re.sub(r'MOMENTUM\s*=\s*[^;]*;', 'MOMENTUM = 6GEV;', FODO).replace('VOLT=1.5E6', 'VOLT=3E7')


def _nt(opt):
    src = ('ON RAD COD RFSW RADCOD;\n' + FODO + '\nFFS;\nNPARA=1;\nUSE RING;\nCELL; CALC;\nFFS["cell calc emit"];\nNORAD; NOFLUC; LOSSAMPL=0.008;\n'
           'res = DynamicApertureSurvey[{{1,25},{1,1},{0}}, 400, Output->0' + opt + '];\nPrint["NT ", ResultOfDAPERT[[2,1,3]]];\nExit;\n')
    out = io.StringIO()
    run_sad_string(src, out=out)
    m = re.search(r'^NT \{([^}]*)\}', out.getvalue(), re.M | re.S)
    return [int(float(v)) for v in m.group(1).replace('\n', ' ').split(',')]


def test_damping_true_matches_gs_without_sextupoles(k64_mode):
    # gs: {400,400,400,400,400,1,1,1,0,0,…}。Damping->True・DampEnough も動き (エラーにならず)、同じ表になる
    ref = [400] * 5 + [1] * 3 + [0] * 43
    assert _nt(',Damping->True') == ref
    assert _nt(',Damping->True,DampEnough->0.01')[:8] == ref[:8]


def test_tcavi_rankick_ranphase_matches_gs(k64_mode):
    """v2.4.0: TCAVI の RANKICK・RANPHASE (tturn.f icTCAV: 1 周ごとに tgauss を 2 回、ak = K0 + RANKICK·g1、ph = PHI + RANPHASE·g2)。
    gs の 8 粒子・6 周の x・px と 1e-14 で一致 (3 GeV の 24 セル FODO、TCAVI L=0.2 K0=2e-4 FREQ=5e8 PHI=0.3 RANKICK=1e-4 RANPHASE=0.1、RFSW ON・NORAD)"""
    GX = [-0.00017057698673404373, -0.0003633211142255208, -0.0005553707885043064, -0.0007467099833424176, -0.0010323341295079725, -0.0012238859912803808, -0.0014147336433520822, -0.001604892044048937]
    GP = [-3.994622666589469e-05, -0.0001406553492524997, -0.00024151084501171277, -0.0003425124699587055, -0.00046934525500916735, -0.0005704425679040154, -0.000671679513972319, -0.0007730322038962854]
    f = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace('LINE RING = (ORIGIN RF 24*CELL);', 'TCAVI TC = (L=0.2 K0=2E-4 HARM=0 FREQ=5.0E8 PHI=0.3 RANKICK=1E-4 RANPHASE=0.1);\nLINE RING = (ORIGIN RF TC 24*CELL);')
    n = 8
    fmt = lambda v: '{' + ','.join(repr(float(a)) for a in v) + '}'
    xs = [1e-4 * i for i in range(n)]
    ys = [5e-5 * (i % 3) for i in range(n)]
    zs = [1e-3 * (i % 4) for i in range(n)]
    src = (f + '\nFFS;\nNORAD NOFLUC;\nUSE RING; CELL; CALC;\nbeam={%s,%s,%s,%s,%s,%s,%s};\n' % (fmt(xs), fmt([0] * n), fmt(ys), fmt([0] * n), fmt(zs), fmt([0] * n), fmt([1] * n))
           + 'r=TrackParticles[{1,beam},1,1,6];\nPrint["X ",r[[2,1]]];Print["P ",r[[2,2]]];\nExit;\n')
    out = io.StringIO()
    run_sad_string(src, out=out)
    import numpy as np
    for tag, ref in (('X', GX), ('P', GP)):
        m = re.search(r'^' + tag + r' \{([^}]*)\}', out.getvalue(), re.M | re.S)
        v = [float(a) for a in m.group(1).replace('\n', ' ').split(',')]
        assert np.allclose(v, ref, atol=1e-14, rtol=0)


def test_pspac_npara4_matches_gs(k64_mode):
    """v2.4.0: PSPAC の NPARA > 1。gs は粒子を NPARA 個 (npara = min(NPARA, np/9+1, ne·np/3000+1)) の子プロセスに分け、各プロセスは自分の
    粒子だけで電位を作る。4096 粒子・NPARA = 4 の rms が gs と 12 桁一致 (NPARA = 1 の値 .00123400976450 とは違う)"""
    import pytest
    from pysad.tests import test_v2_P4_pspac_k64 as t
    src = t._TPL % dict(nxy=128, dx='1e-4', nz=32, dz='0.005', pn=1, pc=0, lb='0.1', nturn=1)
    src = src.replace('NPARA = 1;', 'NPARA = 4;')
    out = io.StringIO()
    run_sad_string(src, out=out)
    m = re.search(r'^RMS\s+(.*)$', out.getvalue(), re.M)
    got = [float(v) for v in m.group(1).split()]
    ref = [.00123378738693, .00123555127440, .00207201468904, .00207304880925, 1.898841600457E-5, 3.786796913267E-4]
    for g, r in zip(got, ref):
        assert g == pytest.approx(r, rel=3e-12)
