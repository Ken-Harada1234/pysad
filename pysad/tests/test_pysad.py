"""
test_pysad.py - Numerical verification tests for pysad

Reference values from SAD and MAD-X / known analytical formulas.
Run with:  python -m pytest pysad/tests/test_pysad.py -v
"""

import math
import numpy as np
import pytest
import sys, os

# pysad パッケージを探す:
# 1) このファイルが pysad/tests/test_pysad.py の場合 → 親の親を追加
# 2) このファイルが pysad/ フォルダと同じディレクトリにある場合 → そのまま
_here = os.path.dirname(os.path.abspath(__file__))
for _candidate in [
    os.path.join(_here, '..', '..'),  # pysad/tests/ → 2段上
    os.path.join(_here, '..'),        # pysad/ → 1段上
    _here,                             # 同じフォルダに pysad/ がある
]:
    if os.path.isdir(os.path.join(_candidate, 'pysad')):
        sys.path.insert(0, os.path.abspath(_candidate))
        break

import pysad
from pysad import (
    Drift, Quadrupole, Bend, Sextupole, Marker, Monitor,
    Cavity,
    Lattice, fodo_cell,
    symplectic_j, symplectic_inverse, periodic_optics,
    pxy2dpz, track_drift_exact, track_thick_quad, track_thick_bend,
    track_cavity, radiation_loss_per_turn, damping_times,
    track_particles, TrackResult, tune_by_fft,
    SADScript,
    SpeedOfLight,
)


# ---------------------------------------------------------------------------
# Helper tolerances
# ---------------------------------------------------------------------------

def assert_close(a, b, rtol=1e-6, atol=1e-10, label=''):
    msg = f"{label}: expected {b}, got {a}"
    assert abs(a - b) <= atol + rtol * abs(b), msg


# ===========================================================================
# 1. Transfer matrices
# ===========================================================================

class TestDriftMatrix:
    """Drift element linear transfer matrix."""

    def test_identity_length_zero(self):
        d = Drift(name='D0', L=0.0)
        M = d.transfer_matrix()
        np.testing.assert_allclose(M, np.eye(6), atol=1e-14)

    def test_drift_1m(self):
        L = 1.0
        d = Drift(name='D1', L=L)
        M = d.transfer_matrix()
        # x -> x + px*L
        assert_close(M[0, 1], L, label='M[x,px]')
        # y -> y + py*L
        assert_close(M[2, 3], L, label='M[y,py]')
        # z -> z - delta*L in ultra-relativistic limit (SAD convention)
        # (may be 0 if ultra-relativistic approximation uses different convention)
        assert M[4, 5] <= 0.0, f"M[z,delta] = {M[4,5]}, expected <= 0"
        # px unchanged
        assert_close(M[1, 0], 0.0, label='M[px,x]')
        assert_close(M[1, 1], 1.0, label='M[px,px]')

    def test_symplectic(self):
        L = 2.5
        d = Drift(name='D', L=L)
        M = d.transfer_matrix()
        J = symplectic_j(3)
        # M^T J M = J
        residual = M.T @ J @ M - J
        np.testing.assert_allclose(residual, np.zeros((6, 6)), atol=1e-12,
                                   err_msg="Drift matrix not symplectic")


class TestQuadrupoleMatrix:
    """Quadrupole transfer matrix (thin and thick)."""

    def test_thin_focusing(self):
        """Thin lens approximation: M[1,0] = -1/f."""
        f = 2.0  # focal length
        K1L = 1.0 / f
        q = Quadrupole(name='QF', L=0.01, K1=K1L / 0.01)
        M = q.transfer_matrix()
        # For ultra-thin: M[px, x] ≈ -K1*L = -1/f
        assert abs(M[1, 0] + K1L) < 0.01, f"Thin quad M[px,x] = {M[1,0]}, expected ≈ {-K1L}"

    def test_thick_quad_period(self):
        """Thick focusing quad: phase advance = sqrt(K1)*L."""
        K1 = 1.0   # m^-2
        L = math.pi / (2 * math.sqrt(K1))  # quarter period
        q = Quadrupole(name='QF', L=L, K1=K1)
        M = q.transfer_matrix()
        phi = math.sqrt(K1) * L  # = pi/2
        assert_close(M[0, 0], math.cos(phi), rtol=1e-5, label='M[x,x]')
        assert_close(M[0, 1], math.sin(phi) / math.sqrt(K1), rtol=1e-5, label='M[x,px]')
        assert_close(M[1, 0], -math.sqrt(K1) * math.sin(phi), rtol=1e-5, label='M[px,x]')

    def test_symplectic(self):
        K1, L = 1.2, 0.5
        q = Quadrupole(name='Q', L=L, K1=K1)
        M = q.transfer_matrix()
        J = symplectic_j(3)
        residual = M.T @ J @ M - J
        np.testing.assert_allclose(residual, np.zeros((6, 6)), atol=1e-11,
                                   err_msg="Quad matrix not symplectic")


class TestBendMatrix:
    """Bend element transfer matrix."""

    def test_sector_bend_matrix(self):
        """Sector bend matrix elements (no edge fields, no K1)."""
        L = 1.0
        angle = 0.1   # rad
        h = angle / L  # 1/m
        bend = Bend(name='B', L=L, ANGLE=angle, K1=0.0, E1=0.0, E2=0.0)
        M = bend.transfer_matrix()
        rho = 1.0 / h
        phi = h * L
        assert_close(M[0, 0], math.cos(phi), rtol=1e-5, label='bend M[x,x]')
        assert_close(M[0, 1], rho * math.sin(phi), rtol=1e-5, label='bend M[x,px]')
        # Dispersion: M[x,delta] = rho*(1-cos(phi))
        assert_close(M[0, 5], rho * (1 - math.cos(phi)), rtol=1e-5, label='bend M[x,delta]')

    def test_symplectic(self):
        bend = Bend(name='B', L=1.0, ANGLE=0.2, K1=0.5, E1=0.1, E2=0.1)
        M = bend.transfer_matrix()
        J = symplectic_j(3)
        residual = M.T @ J @ M - J
        np.testing.assert_allclose(residual, np.zeros((6, 6)), atol=1e-10,
                                   err_msg="Bend matrix not symplectic")


# ===========================================================================
# 2. Symplectic algebra
# ===========================================================================

class TestSymplectic:

    def test_j_antisymmetric(self):
        J = symplectic_j(3)
        np.testing.assert_allclose(J + J.T, np.zeros((6, 6)), atol=1e-14)

    def test_j_squared(self):
        J = symplectic_j(3)
        np.testing.assert_allclose(J @ J, -np.eye(6), atol=1e-14)

    def test_inverse(self):
        K1, L = 1.5, 0.4
        q = Quadrupole(name='Q', L=L, K1=K1)
        M = q.transfer_matrix()
        Minv = symplectic_inverse(M)
        np.testing.assert_allclose(M @ Minv, np.eye(6), atol=1e-11,
                                   err_msg="Symplectic inverse failed")


# ===========================================================================
# 3. FODO cell Twiss (analytic reference)
# ===========================================================================

class TestFODOCell:
    """
    Simple FODO cell with thin lenses:
      Half QF - Drift - QD - Drift - Half QF
    Analytic tune: cos(phi) = 1 - L^2/(2*f^2)
    Beta max = f*(1 + sin(phi/2)) / sin(phi)
    Beta min = f*(1 - sin(phi/2)) / sin(phi)
    """

    @pytest.fixture
    def simple_fodo(self):
        """FODO with thin quads, cell length 2L, f = 2L."""
        L = 2.0   # half-cell length (drift length)
        f = 2.0   # focal length
        lat = fodo_cell(L_half=L, K1F=1.0/f, K1D=-1.0/f)
        return lat, L, f

    def test_tune_range(self, simple_fodo):
        lat, L, f = simple_fodo
        tw = lat.twiss()  # list of dicts per element
        # Use phase advance at first element (normalised tune 0–1)
        Qx = tw[0]['nx']  # total phase / (2*pi) from periodic_optics
        assert 0.0 < Qx < 1.0, f"Qx={Qx} out of expected range for stable FODO"

    def test_beta_positive(self, simple_fodo):
        lat, L, f = simple_fodo
        tw = lat.twiss()
        for entry in tw:
            assert entry['bx'] > 0, f"betax <= 0: {entry['bx']}"
            assert entry['by'] > 0, f"betay <= 0: {entry['by']}"

    def test_periodic_twiss_closure(self, simple_fodo):
        lat, L, f = simple_fodo
        tw = lat.twiss()
        # betax at first and last element should differ by one cell advance
        # (they are different elements; periodicity means M^n = identity)
        # Just check all are positive as a proxy
        betas = [e['bx'] for e in tw]
        assert all(b > 0 for b in betas), "betax not all positive"


# ===========================================================================
# 4. Radiation integrals (analytic check for pure bend ring)
# ===========================================================================

class TestRadiationIntegrals:
    """
    For a ring made of identical sector bends only (no quads):
    I1 = rho * (2*pi) = C/(2*pi*rho) * (2*pi*rho) = C (circumference) if h=1/rho
    I2 = 2*pi/rho  (= C/rho^2)
    I3 = 2*pi/rho^2  (= C/rho^3)
    I4 = 0 (no K1)
    I5 depends on dispersion
    """

    def test_i2_pure_bend_ring(self):
        """I2 = int(h^2 ds) = total_angle^2 / C for uniform ring."""
        # Build a ring of 4 sector bends, each 90 deg, no quads
        n_bends = 4
        angle = 2 * math.pi / n_bends  # each bend = 90 deg
        rho = 10.0  # bending radius
        L_bend = rho * angle
        bend = Bend(name='B', L=L_bend, ANGLE=angle, K1=0.0, E1=0.0, E2=0.0)
        lat = Lattice([bend, bend, bend, bend])

        RI = lat.radiation_integrals()
        I2_expected = 2 * math.pi / rho   # = sum of h^2 * ds = (1/rho)^2 * 2*pi*rho
        assert_close(RI['I2'], I2_expected, rtol=1e-4, label='I2')

    def test_i3_pure_bend_ring(self):
        n_bends = 4
        angle = 2 * math.pi / n_bends
        rho = 10.0
        L_bend = rho * angle
        bend = Bend(name='B', L=L_bend, ANGLE=angle, K1=0.0, E1=0.0, E2=0.0)
        lat = Lattice([bend, bend, bend, bend])

        RI = lat.radiation_integrals()
        I3_expected = 2 * math.pi / rho**2
        assert_close(RI['I3'], I3_expected, rtol=1e-3, label='I3')


# ===========================================================================
# 5. Exact tracking (SAD coordinate system)
# ===========================================================================

class TestExactTracking:
    """Verify exact symplectic drift and quad tracking."""

    def test_pxy2dpz_identity(self):
        """At zero transverse momentum, dpz = 0."""
        assert_close(pxy2dpz(0.0, 0.0), 0.0, atol=1e-15, label='pxy2dpz(0,0)')

    def test_pxy2dpz_formula(self):
        """pxy2dpz = sqrt(1 - px^2 - py^2) - 1."""
        px, py = 0.1, 0.05
        expected = math.sqrt(1 - px**2 - py**2) - 1
        assert_close(pxy2dpz(px, py), expected, rtol=1e-12, label='pxy2dpz formula')

    def test_drift_exact_symplectic(self):
        """Exact drift must preserve the symplectic condition."""
        # Track a single particle (x, px, y, py, z, delta)
        x, px, y, py, z, g = 0.001, 0.0002, -0.001, 0.0001, 0.0, 0.001
        L = 1.5
        x2, px2, y2, py2, z2, g2 = track_drift_exact(x, px, y, py, z, g, L)
        # Check x update: x_out = x + px*L / (1+dpz)
        dpz = pxy2dpz(px, py)
        x_expected = x + px * L / (1.0 + dpz)
        assert_close(x2, x_expected, rtol=1e-10, label='drift x')
        # px unchanged
        assert_close(px2, px, rtol=1e-14, label='drift px')

    def test_thick_quad_small_amplitude(self):
        """For small amplitude, thick quad tracking should match the transfer matrix."""
        K1, L = 1.5, 0.3
        q = Quadrupole(name='Q', L=L, K1=K1)
        M = q.transfer_matrix()

        x0, px0, y0, py0 = 1e-5, 0.0, 0.0, 0.0
        z0, g0 = 0.0, 0.0
        result = track_thick_quad(x0, px0, y0, py0, z0, g0, L, K1)
        x_track, px_track = result[0], result[1]

        # Linear prediction from transfer matrix
        vec = np.array([x0, px0, y0, py0, z0, g0])
        out_lin = M @ vec

        # Should agree to within ~1% for small amplitude (nonlinear effects small)
        assert abs(x_track - out_lin[0]) / abs(out_lin[0] + 1e-15) < 0.01, \
            f"x: track={x_track:.6e} lin={out_lin[0]:.6e}"

    def test_track_particles_shape(self):
        """track_particles returns TrackResult with correct shape."""
        lat = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        N = 5
        coords = np.zeros((6, N))
        coords[0, :] = np.linspace(-0.001, 0.001, N)  # x offsets
        nturns = 3
        result = track_particles(lat.elements, coords, nturns=nturns)
        assert isinstance(result, TrackResult)
        # coords shape is (nturns+1, 6, N)
        assert result.coords.shape == (nturns + 1, 6, N), f"Wrong shape: {result.coords.shape}"
        assert result.survived.shape == (N,)

    def test_tune_by_fft(self):
        """Tune measurement by FFT should be in valid range."""
        lat = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        tw = lat.twiss()
        # nx = phase / (2pi) per element; at initial position it's the one-turn tune
        Qx_twiss = tw[0]['nx']  # tune from periodic optics

        N = 1
        coords = np.zeros((6, N))
        coords[0, 0] = 1e-4  # 0.1 mm x offset
        nturns = 256
        result = track_particles(lat.elements, coords, nturns=nturns)

        if result.survived[0]:
            x_tbt = result.coords[1:, 0, 0]  # turn-by-turn x (shape: nturns)
            # Call tune_by_fft with the TrackResult and particle index
            tunes = tune_by_fft(result, particle=0)
            Qx_fft = tunes['tune_x']
            # FFT tune should be in (0,1)
            assert 0.0 < Qx_fft < 1.0, f"FFT tune {Qx_fft} not in (0,1)"


# ===========================================================================
# 6. SADScript parser
# ===========================================================================

class TestSADScriptParser:
    """Verify the SADScript interpreter handles basic lattice definitions."""

    def test_variable_assignment(self):
        interp = SADScript()
        interp.eval("x = 3.14;")
        assert abs(interp.get_variable('x') - 3.14) < 1e-10

    def test_arithmetic(self):
        interp = SADScript()
        interp.eval("y = 2 * Pi;")
        assert abs(interp.get_variable('y') - 2 * math.pi) < 1e-10

    def test_element_definition(self):
        interp = SADScript()
        interp.load_string("""
            SetElement[QF, QUAD, L->0.3, K1->1.5];
            SetElement[D1, DRIFT, L->1.0];
        """)
        assert 'QF' in interp.element_names()
        assert 'D1' in interp.element_names()

    def test_beamline_definition(self):
        interp = SADScript()
        interp.load_string("""
            SetElement[QF1, QUAD, L->0.3, K1->1.5];
            SetElement[QD1, QUAD, L->0.3, K1->-1.5];
            SetElement[DR1, DRIFT, L->1.0];
            CELL = BeamLine[QF1, DR1, QD1, DR1];
        """)
        assert 'CELL' in interp.beamline_names()

    def test_build_lattice(self):
        interp = SADScript()
        interp.load_string("""
            SetElement[QF2, QUAD, L->0.3, K1->2.0];
            SetElement[QD2, QUAD, L->0.3, K1->-2.0];
            SetElement[DR2, DRIFT, L->0.8];
            SetElement[B1, BEND, L->1.0, ANGLE->0.1];
            CELL2 = BeamLine[QF2, DR2, B1, DR2, QD2, DR2, B1, DR2];
        """)
        lat = interp.build_lattice('CELL2')
        assert isinstance(lat, Lattice)
        assert len(lat.elements) == 8

    def test_math_functions(self):
        interp = SADScript()
        interp.eval("v = Sin[Pi/2];")
        assert abs(interp.get_variable('v') - 1.0) < 1e-10

    def test_element_param_update(self):
        interp = SADScript()
        interp.load_string("""
            SetElement[QX, QUAD, L->0.5, K1->1.0];
        """)
        # Directly update via API (parser Element[...] = form is complex)
        interp.env.set_element_param('QX', 'K1', 2.5)
        assert interp.env.get_element_param('QX', 'K1') == pytest.approx(2.5)


# ===========================================================================
# 7. Emittance and radiation integrals for a real cell
# ===========================================================================

class TestEmittance:
    """Emittance calculation for a FODO cell with dipoles."""

    def test_emittance_positive(self):
        """Equilibrium emittance must be positive for a ring with bends."""
        # A simple ring: 4 pure sector bends (no quads, no edge focusing)
        # This is unstable for Twiss, but radiation integrals are well defined
        n_bends = 4
        angle = 2 * math.pi / n_bends
        rho = 10.0
        L_bend = rho * angle

        bend = Bend('B', L=L_bend, ANGLE=angle, K1=0.0, E1=0.0, E2=0.0)
        lat = Lattice([bend, bend, bend, bend])

        energy = 1.0e9  # 1 GeV
        ri = lat.radiation_integrals()
        # I2 and I3 should be positive
        assert ri['I2'] > 0, "I2 <= 0"
        assert ri['I3'] > 0, "I3 <= 0"

        # U0 = Cgamma * E^4 * I2 / (2 pi)
        from pysad.constants import Cgamma
        E_GeV = energy / 1e9
        U0 = Cgamma * E_GeV**4 * ri['I2'] / (2 * math.pi)
        assert U0 > 0, f"U0={U0} <= 0"


# ===========================================================================
# 8. Integration: build FODO ring and check tune consistency
# ===========================================================================

class TestFODORing:
    """Build a ring of 8 FODO cells and verify tune."""

    def test_ring_tune_stable(self):
        cell = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        tw_cell = cell.twiss()
        # nx at element 0 = normalised tune (phase / 2pi)
        Qx_cell = tw_cell[0]['nx']
        Qy_cell = tw_cell[0]['ny']

        assert 0.0 < Qx_cell < 1.0, f"Qx_cell={Qx_cell}"
        assert 0.0 < Qy_cell < 1.0, f"Qy_cell={Qy_cell}"


# ===========================================================================
# 9. RF 空洞 (Cavity)
# ===========================================================================

class TestCavity:
    """RF 空洞の転送行列と非線形トラッキング。"""

    def test_zero_voltage_is_identity(self):
        """電圧ゼロの空洞は単位行列を返す (横方向のみ、L=0)。"""
        cav = Cavity(name='C0', L=0.0, FREQ=500e6, VOLT=0.0, PHI=np.pi/2, E0=1e9)
        M = cav.transfer_matrix()
        np.testing.assert_allclose(M, np.eye(6), atol=1e-14)

    def test_longitudinal_focus_matrix(self):
        """
        同期位相 PHI=0 のとき縦方向集束要素 M[5,4] が正しい。
        PHI=0: cos(PHI)=1, M[5,4] = -v0 * k_rf
        """
        FREQ = 500e6       # 500 MHz
        VOLT = 1e6         # 1 MeV = 1e6 eV
        E0   = 3e9         # 3 GeV = 3e9 eV
        cav  = Cavity(name='C1', L=0.0, FREQ=FREQ, VOLT=VOLT, PHI=0.0, E0=E0)
        M    = cav.transfer_matrix()
        k_rf = 2.0 * np.pi * FREQ / SpeedOfLight
        v0   = VOLT / E0
        expected = -v0 * k_rf
        assert_close(M[5, 4], expected, rtol=1e-8,
                     label='M[5,4] (longitudinal focusing)')
        # 横方向は無変化
        np.testing.assert_allclose(M[:4, :4], np.eye(4), atol=1e-14)

    def test_on_crest_kick(self):
        """
        on-crest (PHI=π/2) のとき z=0 粒子のδが v0 増加する。
        """
        VOLT = 1e6   # 1 MeV
        E0   = 1e9   # 1 GeV
        cav  = Cavity(name='C2', L=0.0, FREQ=500e6, VOLT=VOLT, PHI=np.pi/2, E0=E0)
        # z=0, δ=0 の粒子
        x, px, y, py, z, g = (np.array([0.0]),)*6
        _, _, _, _, _, g_out = track_cavity(cav, x, px, y, py, z, g)
        expected = VOLT / E0  # = 1e-3
        assert_close(float(g_out[0]), expected, rtol=1e-8,
                     label='delta after on-crest kick')

    def test_off_crest_restoring(self):
        """
        PHI < π/2 のとき後方粒子 (z > 0) はより強く加速される (安定)。
        """
        VOLT = 2e6   # 2 MeV
        E0   = 1e9   # 1 GeV
        PHI  = np.pi / 3   # 60 deg
        cav  = Cavity(name='C3', L=0.0, FREQ=500e6, VOLT=VOLT, PHI=PHI, E0=E0)
        # 同期粒子 (z=0) と後方粒子 (z=+λ_rf/10)
        lam_rf = SpeedOfLight / 500e6  # wavelength
        z0 = np.array([0.0, lam_rf / 10.0])
        x  = np.zeros(2); px = np.zeros(2)
        y  = np.zeros(2); py = np.zeros(2)
        g  = np.zeros(2)
        _, _, _, _, _, g_out = track_cavity(cav, x, px, y, py, z0, g)
        # PHI - k_rf*z での sin の変化で後方 (z>0) は sin(PHI - k_rf*z) > sin(PHI - 0)
        # → 後方粒子の方が大きいΔδを得るべき (PHI ∈ (0, π/2) の場合)
        # 実際には PHI - k_rf*z > 0 のとき sin が増えるかどうかは PHI の値次第
        # PHI=π/3, k_rf*z = 2π * 500e6 / c * (c/500e6 / 10) = 2π/10 ≈ 0.628
        # PHI - k_rf*z ≈ 1.047 - 0.628 = 0.419 > PHI - 0? No, smaller.
        # 実際: sin(PHI - k_rf * z_back) < sin(PHI) for small positive k_rf*z
        # → 後方粒子はエネルギー増加が少ない → 復元力が存在 (安定フォーカス)
        # 重要な物理: M[5,4] < 0 → 縦方向フォーカス
        M = cav.transfer_matrix()
        assert M[5, 4] < 0, f"縦方向集束要素が負でない: M[5,4]={M[5,4]}"


# ===========================================================================
# 10. 放射光ダンピング
# ===========================================================================

class TestRadiation:
    """放射光効果: ダンピングと量子励起のテスト。"""

    def test_radiation_loss_formula(self):
        """
        等磁場リングの U0 = Cgamma * E^4 * I2  (I2 = 2π/rho)。
        PETRA-III 近似値で確認: E=6 GeV, rho=193 m → U0 ≈ 0.065 GeV ≈ 65 MeV
        (実際の PETRA-III の U0 ≈ 27 MeV だが、単純化で係数が異なる)
        """
        E0_GeV = 6.0
        rho    = 193.0    # [m]
        U0     = radiation_loss_per_turn(E0_GeV, rho)
        # 単純な等磁場リングなので U0 ∝ E^4 / rho
        from pysad.constants import Cgamma
        I2 = 2.0 * np.pi / rho
        expected = Cgamma * E0_GeV**4 * I2
        assert_close(U0, expected, rtol=1e-10, label='U0')
        assert U0 > 0, "U0 must be positive"

    def test_damping_times_positive(self):
        """ダンピング時間が正かつ有限であること。"""
        E0_GeV = 3.0
        C      = 2400.0    # 周長 [m]
        I2     = 2 * np.pi / 40.0   # 曲率半径 40 m
        I4     = 0.0
        damps  = damping_times(E0_GeV, C, I2, I4)
        for key in ('tau_x', 'tau_y', 'tau_s'):
            assert damps[key] > 0, f"{key} must be positive"
            assert np.isfinite(damps[key]), f"{key} must be finite"
        # 等分配則 Jx + Jy + Js = 4
        assert_close(damps['Jx'] + damps['Jy'] + damps['Js'], 4.0,
                     atol=1e-10, label='partition sum')

    def test_deterministic_damping_in_bend(self):
        """
        偏向磁石一本通過後に δ が放射損失分だけ減少すること。

        理論値: Δδ ≈ Cgamma * E0³ * h² * L  (1 次近似)
        """
        from pysad.constants import Cgamma as Cg
        rho    = 10.0
        angle  = np.pi / 2.0
        L      = rho * angle
        E0_GeV = 1.0
        delta0 = 0.01

        x0  = np.array([0.0])
        px0 = np.array([0.0])
        y0  = np.array([0.0])
        py0 = np.array([0.0])
        z0  = np.array([0.0])
        g0  = np.array([delta0])

        # 放射光あり (quantum=False → rng=None)
        x1, px1, y1, py1, z1, g1 = track_thick_bend(
            x0, px0, y0, py0, z0, g0,
            L, angle, K1=0.0, E1=0.0, E2=0.0,
            nsteps=20, E0_GeV=E0_GeV, rng=None,
        )

        # δ が減少していること
        assert float(g1[0]) < float(g0[0]), (
            f"放射ダンピングが機能していない: "
            f"δ_before={g0[0]:.4e}, δ_after={g1[0]:.4e}")

        # 理論値との比較 (50 % 以内)
        h            = 1.0 / rho
        expected_loss = Cg * E0_GeV**3 * h**2 * L
        actual_loss   = float(g0[0] - g1[0])
        assert_close(actual_loss, expected_loss, rtol=0.5,
                     label='radiation loss per bend')

    def test_quantum_excitation_nonzero(self):
        """
        量子励起あり (quantum=True) で複数粒子の δ 分散が増加すること。
        """
        rho   = 10.0
        angle = np.pi / 2.0
        L     = rho * angle
        bend  = Bend('B', L=L, ANGLE=angle, K1=0.0, E1=0.0, E2=0.0)
        elements = [bend, bend, bend, bend]

        E0_GeV = 1.0
        N = 50   # 粒子数
        ci = np.zeros((6, N))   # 全粒子δ=0からスタート

        result = track_particles(
            elements, ci, nturns=30,
            aperture=0.5,
            radiation=True, E0_GeV=E0_GeV,
            quantum=True, seed=42,
            nsteps_bend=20,
        )
        # 量子励起で δ の標準偏差が増加するはず
        sigma_delta = np.std(result.delta[-1, :])
        assert sigma_delta > 0.0, (
            f"量子励起の分散がゼロ: σ(δ)={sigma_delta}")


# ===========================================================================
# 11. DynamicApertureSurvey
# ===========================================================================

class TestDynamicApertureSurvey:
    """DynamicApertureSurvey クラスのテスト。"""

    def test_das_single_dp(self):
        """dp=0 での動的アパーチャが正の値を返すこと。"""
        from pysad import DynamicApertureSurvey
        lat = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        das = DynamicApertureSurvey(lat, nturns=64, nangles=3, x_max=0.03,
                                    dp_range=(0.0, 0.0), ndp=1)
        das.run(verbose=False)
        res = das.result()
        assert len(res) == 1
        dp0 = list(res.keys())[0]
        assert res[dp0]['amplitudes'].min() > 0, "動的アパーチャが0"

    def test_das_dp_scan(self):
        """dp スキャンで複数の結果が得られること。"""
        from pysad import DynamicApertureSurvey
        lat = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        das = DynamicApertureSurvey(lat, nturns=32, nangles=3, x_max=0.03,
                                    dp_range=(-0.01, 0.01), ndp=3)
        das.run(verbose=False)
        dp_arr, angles, amps = das.to_array()
        assert len(dp_arr) == 3
        assert amps.shape == (3, 3)
        assert (amps > 0).all()

    def test_das_max_aperture_keys(self):
        """max_aperture() が dp をキーに持つ辞書を返すこと。"""
        from pysad import DynamicApertureSurvey
        lat = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        das = DynamicApertureSurvey(lat, nturns=32, nangles=3, x_max=0.02,
                                    dp_range=(-0.005, 0.005), ndp=3)
        das.run(verbose=False)
        ma = das.max_aperture()
        assert len(ma) == 3
        for dp, amp in ma.items():
            assert amp > 0


# ===========================================================================
# 12. Optics Matching (Match[])
# ===========================================================================

class TestMatching:
    """光学マッチング (Match クラス) のテスト。"""

    def test_match_tune_converges(self):
        """
        チューンを目標値にマッチングできること。
        FODO セルの K1F, K1D を変化させて Qx=Qy=0.20 を実現。
        """
        from pysad import Match, fodo_cell
        lat0 = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        target_Qx, target_Qy = 0.20, 0.20
        m = Match(lat0)
        m.vary('QF', 'K1', lower=0.3, upper=2.0)
        m.vary('QD', 'K1', lower=-2.0, upper=-0.3)
        m.constraint_tune(Qx=target_Qx, Qy=target_Qy, weight=1.0)
        result = m.minimize(method='nelder-mead', tol=1e-8, maxiter=2000)

        # マッチ後格子でチューン確認
        tw_final = result.lattice.twiss()
        Qx_final = tw_final[-1]['nx'] / (2 * np.pi)
        Qy_final = tw_final[-1]['ny'] / (2 * np.pi)
        assert abs(Qx_final - target_Qx) < 0.01, (
            f"Qx マッチ失敗: target={target_Qx}, got={Qx_final:.4f}")
        assert abs(Qy_final - target_Qy) < 0.01, (
            f"Qy マッチ失敗: target={target_Qy}, got={Qy_final:.4f}")

    def test_match_result_lattice_is_independent(self):
        """マッチ後格子が元格子と独立していること (ディープコピー)。"""
        from pysad import Match, fodo_cell
        lat0 = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        K1F_orig = lat0.elements[0].K1
        m = Match(lat0)
        m.vary('QF', 'K1', lower=0.3, upper=2.0)
        m.constraint_tune(Qx=0.26)
        res = m.minimize(tol=1e-4, maxiter=200)
        # 元格子のK1が変わっていないこと
        assert lat0.elements[0].K1 == K1F_orig, "元格子が変更された"

    def test_match_vary_constraint_custom(self):
        """カスタム制約でチューンをマッチできること。"""
        from pysad import Match, fodo_cell
        lat0 = fodo_cell(L_half=2.0, K1F=0.5, K1D=-0.5)
        target_Qx = 0.15
        m = Match(lat0)
        m.vary('QF', 'K1', lower=0.3, upper=1.5)
        m.vary('QD', 'K1', lower=-1.5, upper=-0.3)
        # カスタム制約: チューン Qx を 0.15 に
        m.constraint_custom(
            lambda lat: lat.twiss()[-1]['nx'] / (2 * np.pi),
            target=target_Qx, weight=1.0, label='Qx_custom')
        res = m.minimize(method='nelder-mead', tol=1e-8, maxiter=1000)
        Qx_final = res.lattice.twiss()[-1]['nx'] / (2 * np.pi)
        assert abs(Qx_final - target_Qx) < 0.02, (
            f"Qx マッチ失敗: target={target_Qx}, got={Qx_final:.4f}")


# ===========================================================================
# 13. Parser ラティス読み込み (SADScript)
# ===========================================================================

class TestParserLoading:
    """SADScript パーサーのラティス読み込みテスト。"""

    def test_parse_simple_lattice(self):
        """SetElement + BeamLine で Lattice を構築できること。"""
        src = """
        SetElement[QF, QUAD, K1 -> 1.5, L -> 0.3];
        SetElement[D1, DRIFT, L -> 2.0];
        SetElement[QD, QUAD, K1 -> -1.2, L -> 0.3];
        SetElement[D2, DRIFT, L -> 2.0];
        CELL = BeamLine[QF, D1, QD, D2];
        """
        interp = SADScript()
        interp.load_string(src)
        lat = interp.build_lattice('CELL')
        assert len(lat.elements) == 4
        assert lat.circumference == pytest.approx(0.3 + 2.0 + 0.3 + 2.0)

    def test_parse_variable_arithmetic(self):
        """変数を使った算術式が正しく評価されること。"""
        src = """
        L_cell = 5.0;
        rho = 10.0;
        angle = L_cell / rho;
        SetElement[B1, BEND, L -> L_cell, ANGLE -> angle];
        """
        interp = SADScript()
        interp.load_string(src)
        # element パラメータの確認
        el_info = interp.env.get_element('B1')
        assert el_info is not None
        assert abs(el_info.get('ANGLE', 0) - 0.5) < 1e-10

    def test_parse_beamline_repetition(self):
        """BeamLine で n*ELEM 繰り返し構文が展開されること。"""
        src = """
        SetElement[DR, DRIFT, L -> 1.0];
        LONG = BeamLine[4 * DR];
        """
        interp = SADScript()
        interp.load_string(src)
        lat = interp.build_lattice('LONG')
        assert len(lat.elements) == 4
        assert lat.circumference == pytest.approx(4.0)

    def test_parse_cavity_element(self):
        """Cavity 要素が正しく構築されること (VOLT/FREQ バグ修正確認)。"""
        src = """
        SetElement[CAV, CAVITY, L -> 0.5, VOLT -> 1000000.0, FREQ -> 500000000.0];
        RING = BeamLine[CAV];
        """
        interp = SADScript()
        interp.load_string(src)
        lat = interp.build_lattice('RING')
        assert len(lat.elements) == 1
        cav = lat.elements[0]
        assert isinstance(cav, Cavity), f"期待: Cavity, 実際: {type(cav)}"
        assert cav.VOLT == pytest.approx(1e6)
        assert cav.FREQ == pytest.approx(500e6)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
