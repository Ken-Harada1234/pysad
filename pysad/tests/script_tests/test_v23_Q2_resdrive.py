# -*- coding: utf-8 -*-
"""v2.3.0 Q2: 6 極の共鳴駆動項・ADTS・クロマティシティ (pysad/resdrive.py)。
基準は追跡 (pysad の追跡と NAFF)。小さな FODO リング (3 GeV、105.6 m、24 セル)、ASTRID2 (0.58 GeV、45.7 m)、
大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m、3563 要素)。追跡の基準値は互換モード k64 (conftest の既定) で取った。"""
import io
import math
import os

import numpy as np
import pytest

from pysad import resdrive as R
from pysad.parser import load_sad_file, load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BS_N = chr(10)
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'


def fodo(k_sf=4.0, k_sd=-6.0, extra=None, cell=None):
    src = FODO.replace(BLANK, 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k_sf, k_sd) + (extra or ''))
    if cell:
        src = src.replace('(SF DS QF D2 B D2 QD DS SD D1)', cell)
    return load_sad_string(src).build_lattice('RING')


def tw0(lat):
    r = lat.twiss_table()[0]
    return ((float(r.alfx), float(r.betx)), (float(r.alfy), float(r.bety)))


def test_parse_and_working_point():
    assert R.parse_term('h21000') == (2, 1, 0, 0, 0)
    assert R.parse_term('h11001') == (1, 1, 0, 0, 1)
    wp = R.working_point_check(6.696329, 4.906056)
    assert abs(wp['dist_x_half'] - 0.1963) < 1e-3 and abs(wp['dist_y_int'] - 0.0939) < 1e-3
    wp = R.working_point_check(6.5 + 1e-3, 4.1)
    assert any(t[0] == 'x half' for t in wp['near'])


def test_interior_optics_matches_table():
    """素子内の w(s)・η(s) (指数・多項式) が、素子の出口の Twiss の表と一致する (四極成分つきの厚い MULT と SEXT)"""
    lat = fodo(extra='\nMULT QFM=(L=0.3 K1=0.7 K2=9) QDM=(L=0.3 K1=-0.6 K2=-7);',
               cell='(SF DS QFM D2 B D2 QDM DS SD D1)')
    rows = lat.twiss_table()
    nl = R.collect_nonlinear(lat.elements, rows)
    kinds = {e.name: e.kind for e in nl}
    assert kinds['QFM'] == 'xfoc' and kinds['QDM'] == 'xdef' and kinds['SF'] == 'drift'
    for e in nl:
        F, rate, fx, fy = R._funcs_ep(e)
        assert abs(F['wx'].value(e.L, rate) - math.sqrt(e.opt1[0]) * np.exp(1j * e.opt1[2])) < 1e-9
        assert abs(F['wy'].value(e.L, rate) - math.sqrt(e.opt1[3]) * np.exp(1j * e.opt1[5])) < 1e-9
        assert abs(F['eta'].value(e.L, rate) - e.opt1[6]) < 1e-9
        f = R._funcs_num(e, np.array([e.L]))
        assert abs(f['wx'][0] - F['wx'].value(e.L, rate)) < 1e-9


def test_analytic_vs_numeric_integration_thick_mult():
    lat = fodo(extra='\nMULT QFM=(L=0.3 K1=0.7 K2=9) QDM=(L=0.3 K1=-0.6 K2=-7) OM=(L=0.3 K3=40);',
               cell='(SF DS QFM D2 B D2 QDM DS SD D1 OM)')
    a = R.adts(lat)
    n = R.adts(lat, method='numeric', nslice=60)
    for k in ('xx', 'xy', 'yy'):
        assert abs(a.sext[k] - n.sext[k]) < 2e-4 * abs(a.sext[k])
        assert abs(a.octu[k] - n.octu[k]) < 1e-4 * abs(a.octu[k])
    r1 = R.resonance_terms(lat)
    r2 = R.resonance_terms(lat, method='numeric', nslice=300)
    for k in r1.names:
        assert abs(r1.total[k] - r2.total[k]) < 1e-3 * (abs(r1.total[k]) + 1)


def test_thick_to_thin_limit_and_convergence():
    """厚い SEXT を n 個の薄レンズ列にした値は n を増やすと厚い式 (解析) に 1/n² で収束する。L→0 の極限は薄レンズ (MULT L=0) に一致"""
    lat = fodo()
    A = R.adts(lat)
    d = []
    for n in (2, 4, 8, 16):
        T = R.adts(R.slice_thin(lat, n))
        d.append(abs(T.sext['xx'] - A.sext['xx']))
    assert d[0] / d[1] > 3.5 and d[1] / d[2] > 3.5 and d[2] / d[3] > 3.5
    assert d[3] < 2e-3 * abs(A.sext['xx'])
    # 数値積分 (nslice) も同じ
    N = R.adts(lat, method='numeric', nslice=16)
    assert abs(N.sext['xx'] - A.sext['xx']) < 2e-3 * abs(A.sext['xx'])
    # 薄レンズ (L=0 の MULT) と、同じ積分強度の非常に短い SEXT: 長さ 0 の極限で連続
    thin = fodo(0.0, 0.0, extra='\nMULT SFT=(L=0 K2=0.8) SDT=(L=0 K2=-1.2);', cell='(SFT DS QF D2 B D2 QD DS SDT D1)')
    short = fodo(0.0, 0.0, extra='\nSEXT SFS=(L=1e-5 K2=0.8) SDS=(L=1e-5 K2=-1.2);', cell='(SFS DS QF D2 B D2 QD DS SDS D1)')
    t1, t2 = R.adts(thin), R.adts(short)
    for k in ('xx', 'xy', 'yy'):
        assert abs(t1.sext[k] - t2.sext[k]) < 3e-4 * abs(t1.sext[k])
    r1, r2 = R.resonance_terms(thin), R.resonance_terms(short)
    for k in r1.names:
        assert abs(r1.total[k] - r2.total[k]) < 1e-3 * (abs(r1.total[k]) + 1e-3)


def test_reference_shift():
    """基準位置を変えると h は位相の因子だけ変わる (|h| は同じ)。ADTS は基準によらない"""
    lat = fodo()
    r0 = R.resonance_terms(lat, ('h21000', 'h10110', 'h11001'))
    r1 = R.resonance_terms(lat, ('h21000', 'h10110', 'h11001'), ref=100)
    rows = lat.twiss_table()
    mx, my = rows[100].mux, rows[100].muy
    assert abs(r1.total['h21000'] - r0.total['h21000'] * np.exp(-1j * mx)) < 1e-9
    assert abs(r1.total['h10110'] - r0.total['h10110'] * np.exp(-1j * mx)) < 1e-9
    assert abs(r1.total['h11001'] - r0.total['h11001']) < 1e-9


def test_h_conjugate_pairs():
    lat = fodo()
    r = R.resonance_terms(lat, ('h21000', 'h12000', 'h10020', 'h01002'.replace('01002', '01200')))
    assert abs(r.total['h21000'] - np.conj(r.total['h12000'])) < 1e-9


def test_rdt_vs_one_turn_map():
    """ΔJ(1 周) の追跡 (6 極の強さ λ = ±1 の差の半分 = 1 次の部分) が −∂G/∂φ (G = Σ h ...) と一致する。δ の項 (分散) も入る"""
    from pysad.track6n import track_particles_fast
    terms = [(j, k, l, m, p) for j in range(4) for k in range(4) for l in range(4) for m in range(4) for p in range(4) if j + k + l + m + p == 3]
    res = R.resonance_terms(fodo(), terms, bend_cubic=False)         # 6 極の強さの奇数次だけを取り出す比較: 6 極だけの h (偏向磁石の項は R4 の試験)
    Jx, Jy, phx, phy, d = 1e-9, 2e-9, 0.7, 1.9, 1e-4
    px = py = 0.0
    for (j, k, l, m, p), nm in zip(res.terms, res.names):
        h = res.total[nm]
        amp = (2 * Jx) ** ((j + k) / 2) * (2 * Jy) ** ((l + m) / 2) * d ** p
        ph = np.exp(1j * ((j - k) * phx + (l - m) * phy))
        px += -(1j * (j - k) * h * amp * ph).real
        py += -(1j * (l - m) * h * amp * ph).real

    def dj(lam):
        lat = fodo(4 * lam, -6 * lam)
        r = lat.twiss_table()[0]
        ax, bx, ay, by, ex, epx = r.alfx, r.betx, r.alfy, r.bety, r.etx, r.etpx
        c = np.zeros((1, 6))
        c[0, 0] = math.sqrt(2 * Jx * bx) * math.cos(phx) + ex * d
        c[0, 1] = -math.sqrt(2 * Jx / bx) * (math.sin(phx) + ax * math.cos(phx)) + epx * d
        c[0, 2] = math.sqrt(2 * Jy * by) * math.cos(phy)
        c[0, 3] = -math.sqrt(2 * Jy / by) * (math.sin(phy) + ay * math.cos(phy))
        c[0, 5] = d
        o = track_particles_fast(lat.elements, c.copy(), 1, 3e9, 0.0, False, False, False)[0]
        J = lambda x, p_, a, b: (x * x + (a * x + b * p_) ** 2) / b / 2
        x, p_ = o[0] - ex * d, o[1] - epx * d
        x0, p0 = c[0, 0] - ex * d, c[0, 1] - epx * d
        return J(x, p_, ax, bx) - J(x0, p0, ax, bx), J(o[2], o[3], ay, by) - Jy
    a, b = dj(1), dj(-1)
    tx, ty = (a[0] - b[0]) / 2, (a[1] - b[1]) / 2
    assert abs(tx - px) < 0.02 * abs(px) and abs(ty - py) < 0.02 * abs(py)


def _adts_c(make, p, jb, nturn=512):
    """追跡の 6 極の 2 次の項 C = (f(+1) + f(−1) − 2 f(0))/2 (6 極の強さを ±1, 0 倍)"""
    f = {}
    for lam in (1, -1, 0):
        lat = make(lam)
        t = R.adts_tracking(lat.elements, p, tw0(lat), jref=(jb, jb), nturn=nturn)
        f[lam] = np.array([t['xx'], t['xy'], t['yy']])
    return (f[1] + f[-1] - 2 * f[0]) / 2


def test_adts_fodo_vs_tracking():
    C = _adts_c(lambda l: fodo(4 * l, -6 * l), 3e9, 1e-8)
    A = R.adts(fodo())
    ref = np.array([A.sext['xx'], A.sext['xy'], A.sext['yy']])
    assert np.all(np.abs(C - ref) < 1e-3 * np.abs(ref))


def test_octupole_first_order_formula_and_tracking():
    lat = fodo(0.0, 0.0, extra='\nMULT O1=(L=0 K3=50) O2=(L=0 K3=-30);', cell='(O1 DS QF D2 B D2 QD DS O2 D1)')
    rows = lat.twiss_table()
    b = [(rows[i].betx, rows[i].bety) for i, e in enumerate(lat.elements) if e.name in ('O1', 'O2')]
    # 薄い 8 極の式: ∂νx/∂Jx = Σ K3L βx²/(16π)、∂νx/∂Jy = −Σ K3L βxβy/(8π)、∂νy/∂Jy = Σ K3L βy²/(16π)  (24 セル)
    # ここでの b は行の出口側。薄レンズは入口と出口が同じ
    k = [50.0, -30.0]
    ex = 24 * sum(ki * bi[0] ** 2 for ki, bi in zip(k, b)) / (16 * math.pi)
    exy = -24 * sum(ki * bi[0] * bi[1] for ki, bi in zip(k, b)) / (8 * math.pi)
    eyy = 24 * sum(ki * bi[1] ** 2 for ki, bi in zip(k, b)) / (16 * math.pi)
    A = R.adts(lat)
    assert abs(A.octu['xx'] - ex) < 1e-9 * abs(ex) and abs(A.octu['xy'] - exy) < 1e-9 * abs(exy) and abs(A.octu['yy'] - eyy) < 1e-9 * abs(eyy)
    # 追跡: 8 極の符号を反転した 2 本の差の半分 (偶数次・8 極以外の項を消す) = 8 極の 1 次
    latm = fodo(0.0, 0.0, extra=BS_N + 'MULT O1=(L=0 K3=-50) O2=(L=0 K3=30);', cell='(O1 DS QF D2 B D2 QD DS O2 D1)')
    t = R.adts_tracking(lat.elements, 3e9, tw0(lat), jref=(5e-9, 5e-9), nturn=512)
    tm = R.adts_tracking(latm.elements, 3e9, tw0(latm), jref=(5e-9, 5e-9), nturn=512)
    for key in ('xx', 'xy', 'yy'):
        assert abs((t[key] - tm[key]) / 2 - A.octu[key]) < 2e-3 * abs(A.octu[key])


def test_chromaticity_sextupole_part():
    lat = fodo()
    lat0 = fodo(0.0, 0.0)
    sx, sy = R.chromaticity_sextupole(lat)
    c1, c0 = R.chromaticity_numeric(lat, order=2), R.chromaticity_numeric(lat0, order=2)
    assert abs((c1['d1'][0] - c0['d1'][0]) - sx) < 2e-3 * abs(sx)
    assert abs((c1['d1'][1] - c0['d1'][1]) - sy) < 2e-3 * abs(sy)
    # h11001/π = (1/4π) Σ K2L η β (薄い式) との対応: 6 極の h11001 は実数
    r = R.resonance_terms(lat, ('h11001', 'h00111'))
    assert abs(r.total['h11001'].imag) < 1e-12 and abs(r.total['h00111'].imag) < 1e-12


def test_tilted_and_skew_sextupole_kick_convention():
    """回転 (ROTATE) と SK2 の符号の向きが pysad の追跡の kick と一致する (解析の係数の前提)"""
    from pysad.track6n import track_particles_fast
    for decl, c in (('MULT S0=(L=0 K2=0.6 ROTATE=0.2);', 0.6 * complex(math.cos(0.6), math.sin(0.6))),
                    ('MULT S0=(L=0 K2=0.6 SK2=0.4);', complex(0.6, 0.4))):
        el = load_sad_string('MOMENTUM=3GEV;' + decl + ' LINE A=(S0);').build_lattice('A').elements[0]
        x, y = 1e-3, 2e-3
        o = track_particles_fast([el], np.array([[x, 0, y, 0, 0, 0]]), 1, 3e9, 0.0, False, False, False)[0]
        P = R.poly_coeffs(2, c)
        dx = sum(v * a * x ** (a - 1) * y ** b for (a, b), v in P.items() if a > 0)
        dy = sum(v * b * x ** a * y ** (b - 1) for (a, b), v in P.items() if b > 0)
        assert abs(o[1] + dx) < 1e-12 and abs(o[3] + dy) < 1e-12
    lat = fodo(0.0, 0.0, extra='\nMULT T1=(L=0 K2=0.8 ROTATE=0.3);', cell='(T1 DS QF D2 B D2 QD DS SD D1)')
    ra = R.collect_nonlinear(lat.elements, lat.twiss_table())
    assert abs(ra[0].strengths[2] - 0.8 * complex(math.cos(0.9), math.sin(0.9))) < 1e-12


def test_astrid2_adts_vs_tracking():
    path = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
    src = open(path, encoding='utf-8').read()
    k = {}
    sc = load_sad_string(src)
    lat1 = sc.build_lattice('RING')
    base = {e.name: float(e._raw['K2']) for e in lat1.elements if type(e).__name__ == 'Sextupole'}

    def make(l):
        lat = load_sad_string(src).build_lattice('RING')
        for e in lat.elements:
            if type(e).__name__ == 'Sextupole':
                e.K2 = e.K2 * l
        return lat
    C = _adts_c(make, 0.58e9, 2e-9)
    A = R.adts(lat1)
    ref = np.array([A.sext['xx'], A.sext['xy'], A.sext['yy']])
    assert np.all(np.abs(C - ref) < 3e-3 * np.abs(ref))


def test_large_ring_a_adts_against_recorded_tracking():
    """大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)。追跡 (NAFF、1024 周、J = 2e-10 m·rad) の記録値
    (example/14_nonlinear_dynamics/adts_rdt/ADTS_summary.md): 6 極 2 次 C = (86842, −75097, −1608)、8 極 1 次 B = (−36905, 32480, −7154)"""
    path = os.path.join(HERE, '..', 'kek-ls_v03_12.lat')
    if not os.path.exists(path):
        pytest.skip('kek-ls_v03_12.lat がない')
    lat = load_sad_file(path).build_lattice('RING')
    rows = lat.twiss_table()
    a = R.adts(lat, rows=rows)
    sc, oc = (86842.4, -75097.4, -1607.9), (-36905.0, 32480.0, -7154.2)
    for k, s, o in zip(('xx', 'xy', 'yy'), sc, oc):
        assert abs(a.sext[k] - s) < 3e-3 * max(abs(s), 5e3)
        assert abs(a.octu[k] - o) < 1e-3 * abs(o)
    assert a.imag_residual < 1e-6
    n = R.adts(lat, method='numeric', nslice=12, rows=rows)
    for k in ('xx', 'xy'):
        assert abs(n.sext[k] - a.sext[k]) < 2e-3 * abs(a.sext[k])
    # 作動点の整数・半整数共鳴: νx = 47.580034、νy = 17.620044 → 半整数まで 0.0800、0.1200
    assert abs(a.wp['dist_x_half'] - 0.0800) < 1e-3 and abs(a.wp['dist_y_half'] - 0.1200) < 1e-3
    r = R.resonance_terms(lat, rows=rows, bend_cubic=False)         # 6 極のみ (R4 の既定は偏向磁石の項つき: 17.811)
    assert abs(abs(r.total['h30000']) - 17.76) < 0.05
    rb = R.resonance_terms(lat, rows=rows)
    assert abs(abs(rb.total['h30000']) - 17.811) < 0.01


def test_ffs_functions():
    from pysad import ffscmd, compat
    old = compat._mode
    compat.set_mode('v2')
    try:
        src = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)') + '''
FFS;
USE RING; CELL; CALC;
Print["RDT ", ResonanceDrivingTerms[Terms->{"h30000","h21000"}, BendCubic->False]];
Print["ADTS ", AmplitudeTuneShift[]];
Print["CHR ", ChromaticityTerms[Order->2]];
Exit;
'''
        buf = io.StringIO()
        ffscmd.run_sad_string(src, out=buf)
    finally:
        compat.set_mode(old)
    txt = buf.getvalue()
    assert 'H30000->{-.37313' in txt and 'DNuxDJx->323.40955' in txt and 'SextupolePart->{17.4741' in txt
