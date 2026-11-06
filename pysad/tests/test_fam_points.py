"""FITP の fam 点 (MatchingAmplitude / NFAMP: tffamsetup.f の有限振幅ファミリー) の自己検証 (docs/history/work_trpt_acc.md 8 節)。
SAD(WSL) は FITP>1 の GO で再起動するので参照値はなく、構造と整合性を確かめる。"""
import io
import math

import numpy as np

from pysad.ffscmd import run_sad_string

PI2 = 2.0 * math.pi


def _src(k2f=8.0, k2d=-12.0, amp='{{0, 1, 1}}', nfp=4, mode='CELL'):
    return f'''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5) DS = (L=0.2);
BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
SEXT SF = (L=0.1 K2={k2f}) SD = (L=0.1 K2={k2d});
MARK MK = (AX=-1.7 BX=6.3 AY=0.9 BY=3.1);
LINE CELL = (MK QF DS SF D B D QD DS SD D B D);
LINE RING = (8*CELL);
FFS;
USE RING;
{mode};
CALC;
EMITX = 1E-6; EMITY = 1E-7;
DP = 0.001;
fit $$$ NX 2.21 3 NY 1.31 3;
free QF K1 QD K1;
MatchingAmplitude = {amp};
NFAMP = {nfp};
'''


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def test_fam_points_definition():
    ffs, I, text = _run(_src() + 'calc;\n')
    m = ffs.matcher_()
    m._fam_define(1)
    fam = m._fam
    # {0, 1, 1} で NFAMP=4: x の 4 本 (kfam 1..4) と y の 4 本 (kfam -1..-4)、添字は nfr+1=2 から 2,-2,3,-3,4,-4,5,-5 (tffamsetup の nfa)
    assert sorted(fam) == [-5, -4, -3, -2, 2, 3, 4, 5]
    assert [fam[ix]['kfam'] for ix in (2, -2, 3, -3, 4, -4, 5, -5)] == [1, 2, 3, 4, -1, -2, -3, -4]
    assert all(f['jfam'] == 0 for f in fam.values())
    # dfam: 先頭 MARK の Twiss (CALC 後の周期解が書き戻されている) で x0 = sqrt(βx em)、位相 θ = 2π(i-1)/NFAMP
    r0 = ffs.row(1)
    em = 1.1e-6
    x0 = math.sqrt(r0.betx * em)
    px0 = x0 / r0.betx
    y0 = math.sqrt(r0.bety * em)
    py0 = y0 / r0.bety
    for i, ix in enumerate((2, -2, 3, -3)):
        c, s = math.cos(PI2 / 4 * i), math.sin(PI2 / 4 * i)
        exp = [c * x0, s * px0 - r0.alfx / r0.betx * c * x0, 0.0, 0.0]
        np.testing.assert_allclose(fam[ix]['dfam'], exp, atol=1e-15)
    for i, ix in enumerate((4, -4, 5, -5)):
        c, s = math.cos(PI2 / 4 * i), math.sin(PI2 / 4 * i)
        exp = [0.0, 0.0, c * y0, s * py0 - r0.alfy / r0.bety * c * y0]
        np.testing.assert_allclose(fam[ix]['dfam'], exp, atol=1e-15)
    # 評価点: 親 (idp=0) が範囲内なら fam 点も入る。軌道の条件は fam 点では評価しない
    assert sorted(m._eval_idps(1)) == [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
    assert m._eval_idps(0) == [0] + [ix for ix in fam]


def test_fam_points_need_nfr_and_positive_nfamp():
    ffs, I, text = _run(_src(nfp=0) + 'go;\n')
    assert ffs.matcher._fam == {}
    assert 'Non-positive NFAMP' in text
    # 1 点の条件 (nfr=0) では fam 点なし (dpm の範囲が空)
    ffs, I, text = _run(_src().replace('NX 2.21 3 NY 1.31 3', 'NX 2.21 NY 1.31') + 'go;\n')
    assert ffs.matcher._fam == {}
    ffs, I, text = _run(_src(amp='{{0, 1}}') + 'go;\n')
    assert ffs.matcher._fam == {} and 'MatchingAmplitude must be' in text


def test_fam_orbit_and_amplitude_dependent_tune():
    ffs, I, text = _run(_src() + 'go;\n')
    m = ffs.matcher
    rows0 = m._rows_idp[0]
    cod_prev_before = ffs.lat._cod_prev.copy()
    for ix, f in m._fam.items():
        rk = m._rows_idp[ix]
        # fam 点の入口の軌道 = 親 (idp=0) の閉軌道 + dfam、δ は親と同じ (0)
        exp = np.array([rows0[0].x_cod, rows0[0].px_cod, rows0[0].y_cod, rows0[0].py_cod]) + f['dfam']
        np.testing.assert_allclose([rk[0].x_cod, rk[0].px_cod, rk[0].y_cod, rk[0].py_cod], exp, atol=1e-18)
        assert rk[0].delta_cod == 0.0
        # 六極があるので振幅依存のチューン: 軸上と違う (x の振幅なら主に νx)
        if f['kfam'] > 0:
            assert abs(rk[-1].mux - rows0[-1].mux) > 1e-5
        # fam 点は常に安定扱い (optstat(ix)%stabx = fam .or. …)
        assert m.last_stab[ix] == (True, True)
    # 残差に fam 点が入る (軸上だけの残差と違う)。MatchingAmplitude なしは同じ条件で別の残差
    ffs2, I2, text2 = _run(_src(amp='{}') + 'go;\n')
    assert ffs2.matcher._fam == {}
    assert abs(m.residual.r - ffs2.matcher.residual.r) > 1e-9
    # fam 点の計算は閉軌道の初期値や先頭 MARK を変えない (qcod の codfnd=.true. では twmov しない)
    r = m._calc_optics_fixed_default(0.0, np.array([1e-3, 0, 0, 0, 0, 0.0]))
    assert r is not None and abs(r[0].x_cod - 1e-3) < 1e-18
    np.testing.assert_array_equal(ffs.lat._cod_prev, cod_prev_before)


def test_fam_linear_lattice_same_tune_and_go_matches():
    # 六極なし: ずれた軌道のまわりの一周行列は軸上と同じ → fam 点のチューンは軸上と一致し、GO は Matched.
    ffs, I, text = _run(_src(k2f=0.0, k2d=0.0, amp='{{0, 1, 1}, {0.001, 0.5, 0}}') + 'go;\n')
    m = ffs.matcher
    assert 'Matched.' in text
    fam = m._fam
    # 2 番目の項 (dp=0.001) は idp=1 の点のまわり (dpm(1) = 0.0005 ≤ 0.001 < dpm(2) = 0.001?) → dpm の範囲 [0.0005, 0.001) に入らない
    # (dpm(nfr+1) = dp(nfr) = 0.001 なので 0.001 は外) → x のみ 4 本が付かず、{0,1,1} の 8 本だけ
    assert sorted(fam) == [-5, -4, -3, -2, 2, 3, 4, 5]
    # 六極なしでも BEND/ドリフトの非線形で軌道がずれるとチューンは少し変わる (5e-4)。fam 点の周期解はその軌道のまわりの
    # 一周行列 (閉軌道の探索なし) から作る: 独立に一周行列を作って cos μ を比べる
    from pysad import lattice as L
    L.RFSW = False
    r0 = ffs.row(1)                      # GO 後の要素の値で計算し直す (_rows_idp は最後の微分の評価点のもの)
    for ix, f in fam.items():
        c0 = np.array([r0.x_cod, r0.px_cod, r0.y_cod, r0.py_cod, 0.0, 0.0])
        c0[:4] += f['dfam']
        rk = m._calc_optics_fixed_default(0.0, c0)
        M, _, _ = L._one_turn_cod_matrix(ffs.els, 0.0, c0)
        _R, A, B = L.periodic_r(M)          # y にずれた軌道は BEND で x-y 結合を生むので、正規モードのブロック A, B の跡で比べる
        assert abs(math.cos(rk[-1].mux) - 0.5 * float(np.trace(A))) < 1e-9
        assert abs(math.cos(rk[-1].muy) - 0.5 * float(np.trace(B))) < 1e-9
    # dp=0.0009 なら idp=1 の点に付く (親 jfam=1、dp は親と同じ)
    ffs, I, text = _run(_src(k2f=0.0, k2d=0.0, amp='{{0.0009, 0.5, 0}}') + 'calc;\n')
    m = ffs.matcher_()
    m._fam_define(1)
    assert sorted(m._fam) == [-3, -2, 2, 3] and all(f['jfam'] == 1 and f['kfam'] > 0 for f in m._fam.values())
    assert m._fam_dps({0: 0.0, 1: 0.001, -1: -0.001})[2] == 0.001


def test_fam_ins_propagates_from_displaced_orbit():
    ffs, I, text = _run(_src(mode='INS') + 'go;\n')
    m = ffs.matcher
    assert m._fam
    rows0 = m._rows_idp[0]
    for ix, f in m._fam.items():
        rk = m._rows_idp[ix]
        exp = np.array([rows0[0].x_cod, rows0[0].px_cod, rows0[0].y_cod, rows0[0].py_cod]) + f['dfam']
        np.testing.assert_allclose([rk[0].x_cod, rk[0].px_cod, rk[0].y_cod, rk[0].py_cod], exp, atol=1e-18)
        assert rk[0].betx == rows0[0].betx            # INS: 入口の Twiss は MARK のまま、軌道だけ振幅分ずれる


def test_fitvalue_receives_fam_id():
    src = _src(k2f=0.0, k2d=0.0) + '''
ids = {};
FitValue[comp_, fn_, {id_, dp_}, goal_, now_] := (ids = Append[ids, id]; goal);
go;
'''
    ffs, I, text = _run(src)
    ids = I.get('ids')
    assert ids and max(float(v) for v in ids) == 4 and min(float(v) for v in ids) == -4
    assert any(float(v) == 0 for v in ids)
