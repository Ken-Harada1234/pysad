# -*- coding: utf-8 -*-
"""v2 の修正 B38: PHSROT (6 次元 Twiss による位相空間の回転)。K64 は行列がすべて 0 になり動作しない。

根拠:
  - M0 = N Rot N⁻¹ はシンプレクティックで、固有値は exp(±2πi PX)、exp(±2πi PY)、exp(∓2πi PZ)。
  - 結合・分散がなければ x の 2×2 は Courant-Snyder の式 [[cos+α sin, β sin], [-γ sin, cos-α sin]]。
  - リング (ST PHSROT) の CALC の Twiss・チューンが入力に一致。
  - 減衰 D = diag(d)、拡散 B = diag(b) の平衡: ε = b/(2d - d²) (1 回あたり x ← (1-d) x + √b ξ の不動点)。Emittance[] と、
    FLUC つきの追跡の統計 (独立な経路) がこの式に一致。
"""
import io
import math

import numpy as np
import pytest

from pysad import phsrot

PAR = dict(AX=0.5, BX=4.0, PX=0.31, AY=-0.3, BY=2.0, PY=0.22, AZ=0.1, BZ=10.0, PZ=0.01,
           R1=0.02, R2=0.1, R3=-0.01, R4=0.03, EX=0.2, EPX=0.03, EY=0.01, EPY=-0.002, ZX=0.001, ZPX=0.0005)
SM = np.zeros((6, 6))
for _i in range(3):
    SM[2 * _i, 2 * _i + 1] = 1.0
    SM[2 * _i + 1, 2 * _i] = -1.0


def test_b38_symplectic_and_tunes():
    m = phsrot.matrices(PAR)
    M = m['M0']
    assert np.abs(M.T @ SM @ M - SM).max() < 1e-13
    assert np.abs(m['N'] @ m['Ninv'] - np.eye(6)).max() < 1e-13
    ang = sorted(abs(np.angle(np.linalg.eigvals(M))) / (2 * math.pi))
    assert np.allclose(ang, sorted([0.31, 0.31, 0.22, 0.22, 0.01, 0.01]), atol=1e-12)


def test_b38_courant_snyder_block():
    a, b, mu = 0.5, 4.0, 2 * math.pi * 0.31
    M = phsrot.matrices(dict(AX=a, BX=b, PX=0.31, BY=1.0, BZ=1.0))['M0']
    c, s = math.cos(mu), math.sin(mu)
    exp = np.array([[c + a * s, b * s], [-(1 + a * a) / b * s, c - a * s]])
    assert np.abs(M[:2, :2] - exp).max() < 1e-14
    assert np.abs(M[2:, 2:] - np.eye(4)).max() < 1e-14


SRC = r'''MOMENTUM = 3GEV;
PHSROT PR = (AX=0.5 BX=4 PX=0.31 AY=-0.3 BY=2 PY=0.22 AZ=0 BZ=10 PZ=0.01
 D11=1e-4 D22=1e-4 D33=2e-4 D44=2e-4 D55=3e-4 D66=3e-4 B11=2e-12 B22=2e-12 B33=1e-14 B44=1e-14 B55=6e-9 B66=6e-9);
MARK ST = ();
LINE RING = (ST PR);
FFS;
USE RING;
CELL; CALC;
$FORM="20.14";
PageWidth=1999;
Print["NU ",Twiss["NX","$$$"]/2/Pi," ",Twiss["NY","$$$"]/2/Pi," ",Twiss["BX",1]," ",Twiss["AX",1]," ",Twiss["BY",1]," ",Twiss["AY",1]];
e=Emittance[];
em=Emittances/.e; dr=DampingRate/.e; tu=Tunes/.e;
Print["EM ",em[[1]]," ",em[[2]]," ",em[[3]]," ",dr[[1]]," ",dr[[2]]," ",dr[[3]]," ",tu[[1]]," ",tu[[2]]," ",tu[[3]]];
NORAD;
r=TrackParticles[{1,{{1e-3},{2e-4},{-1e-3},{1e-4},{1e-3},{1e-4},{1}}},1,1,1];
Print["TR ",r[[2,1,1]]," ",r[[2,2,1]]," ",r[[2,3,1]]," ",r[[2,4,1]]," ",r[[2,5,1]]," ",r[[2,6,1]]];
Exit;
'''


def _lines(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
            if l.split() and l.split()[0] in ('NU', 'EM', 'TR')}


def test_b38_ring_twiss_emittance_and_tracking(v2_mode):
    o = _lines(SRC)
    assert np.allclose(o['NU'], [0.31, 0.22, 4.0, 0.5, 2.0, -0.3], atol=1e-9)
    d = [1e-4, 2e-4, 3e-4]
    b = [2e-12, 1e-14, 6e-9]
    for k in range(3):
        assert abs(o['EM'][k] / (b[k] / (2 * d[k] - d[k]**2)) - 1.0) < 1e-6
        assert abs(o['EM'][3 + k] / (-d[k]) - 1.0) < 1e-6
    assert np.allclose(o['EM'][6:], [0.31, 0.22, -0.01], atol=1e-9)
    # 追跡 (NORAD): x ← M0 x。px は正準運動量 (TrackParticles の入出力)
    par = dict(AX=0.5, BX=4.0, PX=0.31, AY=-0.3, BY=2.0, PY=0.22, BZ=10.0, PZ=0.01)
    exp = phsrot.matrices(par)['M0'] @ np.array([1e-3, 2e-4, -1e-3, 1e-4, 1e-3, 1e-4])
    assert np.allclose(o['TR'], exp, rtol=0, atol=1e-15)


def test_b38_k64_mode_does_not_track(k64_mode):
    from pysad.elements import PhsRot
    with pytest.raises(NotImplementedError):
        phsrot.compile_phsrot(PhsRot(name='P', params=dict(BX=1.0)))


def test_b38_tracking_equilibrium_matches_formula(v2_mode):
    """FLUC つきの追跡の平衡分布 (独立な経路): <x²> = εx βx、ε = b/(2d - d²)"""
    from pysad import emit6d as _e
    from pysad.elements import PhsRot
    from pysad.track6n import make_state
    d, bx, by = 2e-2, 4.0, 2.0
    par = dict(AX=0.5, BX=bx, PX=0.31, AY=-0.3, BY=by, PY=0.22, BZ=10.0, PZ=0.013,
               B11=2e-10, B22=2e-10, B33=1e-12, B44=1e-12, B55=6e-9, B66=6e-9)
    for i in range(1, 7):
        par['D%d%d' % (i, i)] = d
    pr = phsrot.compile_phsrot(PhsRot(name='P', params=par))
    ph = _e._Phys(3.0e9, 100.0)
    n = 4000
    S = make_state(np.zeros((n, 6)), ph)
    rs = np.array([12345, 1], dtype=np.int64)
    for _ in range(400):
        phsrot.k_phsrot(S, n, *pr, True, True, rs, ph.p0, ph.h0)
    ex, ey, ez = (v / (2 * d - d * d) for v in (2e-10, 1e-12, 6e-9))
    assert abs(np.mean(S[:, 0]**2) / (ex * bx) - 1.0) < 0.08
    assert abs(np.mean(S[:, 2]**2) / (ey * by) - 1.0) < 0.08
    assert abs(np.mean(S[:, 4]**2) / (ez * 10.0) - 1.0) < 0.08
    assert abs(np.mean(S[:, 5]**2) / (ez / 10.0) - 1.0) < 0.08
    # 相関 <x px> = -α ε
    assert abs(np.mean(S[:, 0] * S[:, 1] * (1 + S[:, 5])) / (-0.5 * ex) - 1.0) < 0.15
