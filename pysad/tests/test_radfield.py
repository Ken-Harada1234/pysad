"""RadiationField / RadiationSpectrum (K64 tlinit.f)。RadiationSpectrum は gs と全桁一致 (rf2.sad)。
RadiationField は gs がこのビルドで SIGSEGV するので参照が取れない (docs/19 B49)。等速運動の電荷の場 (Liénard-Wiechert の速度場) と比べる。"""
import math

import numpy as np


def test_radiation_spectrum_matches_gs():
    from pysad.radfield import radiation_spectrum
    tau = [k * 1e-12 for k in range(60)]
    f = [math.exp(-((k - 30.0) / 6) ** 2) * math.cos(k * 0.4) for k in range(60)]
    k, c, s = radiation_spectrum(tau, f, [1e-3, 5e-3, 1e-3])
    assert abs(k[0] - 1256.63706144) < 1e-6 and abs(k[-1] - 6283.18530718) < 1e-6
    refc = [4.0000325e-12, -1.199126e-12, 2.5317349e-14, 1.1736258e-15, -9.067542e-18, 8.4260886e-21]
    refs = [-3.383378e-12, 2.2513645e-12, -2.374868e-13, 4.1272779e-15, -1.167625e-17, 4.9163675e-21]
    for a, b in zip(c, refc):
        assert abs(a - b) <= 2e-6 * abs(b) + 1e-25
    for a, b in zip(s, refs):
        assert abs(a - b) <= 2e-6 * abs(b) + 1e-25


def test_radiation_field_of_uniform_motion():
    """等速直線運動の電荷 (電子) の場は Liénard-Wiechert の速度場に一致する (内部の点で相対 6e-3: スプライン微分の離散化誤差)"""
    from pysad.radfield import radiation_field, _C, _ELEMCH
    beta = 0.8
    m = 200
    dt = 1e-10
    t = np.arange(m) * dt
    tr = np.stack([t, np.zeros(m), np.zeros(m), beta * _C * t], axis=1)
    obs = [0.3, 0.1, 2.0]
    f = radiation_field(tr, obs, charge=-1.0)
    i = 100
    zq = tr[i, 3]
    r = np.array(obs) - np.array([0.0, 0.0, zq])
    R = np.linalg.norm(r)
    n = r / R
    bvec = np.array([0.0, 0.0, beta])
    # 観測点で見た時刻 τ の場は、遅延位置 zq にある電荷が作る (R は遅延位置からの距離)
    q = -1.0 * _ELEMCH
    e_ref = q / (4 * math.pi * 8.8541878128e-12) * (n - bvec) * (1 - beta ** 2) / ((1 - n @ bvec) ** 3 * R ** 2)
    got = np.array([f[i, 1], f[i, 2], f[i, 3]])
    assert np.allclose(got, e_ref, rtol=6e-3, atol=1e-6 * np.abs(e_ref).max())
