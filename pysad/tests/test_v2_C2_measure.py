# -*- coding: utf-8 -*-
"""v2 (担当 C2): TRCOD FIXCOD / MVCOD / STNCOD / ALIGN、ORBITJITTER、PUTMEASURE / GETMEASURE と、MEASURE の残りの確認
(DETR > 1 の分布、TWAKE / LWAKE ON の MEASURE、TRPT の扱い)。実装は pysad/measure.py、記録は docs/history/v2_changes_C2_measure.md。

K64 (gs) は MEASURE が無出力、TRCOD が SIGSEGV (B42) なので gs の参照は無い。根拠は次の独立な経路:
  - FIXCOD の基準: 基準を取った時点の CALC の軌道・分散そのもの。その後の差 Twiss["DDX"] は、薄いキック Δθ の応答の解析式
    Δx = −Δθ·√(β1 β2)·sin(ψ2 − ψ1) (β と ψ は CALC の Twiss)。
  - ALIGN: 区間の磁石を一様に (dx, dy) ずらした軌道 = ずれのない格子に入口の軌道 (−dx, −dy) を入れた軌道 + (dx, dy) (座標の平行移動)。
  - ORBITJITTER: 足された量 = xj × (同じ種から独立に作った乱数列)。
  - PUTMEASURE / GETMEASURE: 書いた値が読み戻る (7 桁)。fit の目標値は tfgetm.f の式。
  - DETR > 1 の分布: ttinit の線形写像 A から Σ0 = A·Aᵀ を作り、1 粒子の追跡で求めた転送行列 M で運んだ M Σ0 Mᵀ が、
    CALC がその位置まで運んだ Twiss (R1..R4 を含む) から作った Σ と一致する (結合のある輸送路。入口が DETR > 1 と DETR < 1 の両方)。
  - LWAKE / TWAKE ON の MEASURE: 一様なウェイク W の解析式。重心の Δδ = −e·N·W/(p0c)、Δpx = e·N·W1·<x>/(p0c)。
    階段状の縦ウェイクでは平均は半分、広がりは (e·N·W)²/12 だけ増える (z の分布が効いていること)。
k64 モードでは、これらのコマンドは従来どおり動かない (K64 の動作を変えない)。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import measure as ms
from pysad.ffscmd import run_sad_string
from pysad.sadrandom import SadRandom

ECHARGE = 1.602176634e-19


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def _loc(ffs, name):
    return int(ffs._ielmf(name)[0])


# ---------------------------------------------------------------------------- 格子
# 薄いステアリング ST (L=0) を持つ輸送路。入口の軌道は 0
TL = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 EX=0.2 EPX=0.01 %s) MM=();
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
BEND B1 = (L=1 ANGLE=0.05) ST=(L=0 ANGLE=0 K0=0);
LINE TL = (M0 QF D1 B1 D1 QD ST D1 MM B1 D1 QF D1 QD D1);
FFS;
USE TL;
INS; CALC;
'''
COD8 = ('EX', 'EPX', 'EY', 'EPY', 'DX', 'DPX', 'DY', 'DPY')


def _tw(ffs, key, k):
    return ffs.twiss_value(key, k)


# ============================================================================ TRCOD FIXCOD / MVCOD / STNCOD
def test_fixcod_stores_reference_and_twiss_stays_absolute(v2_mode):
    ffs, I, txt = _run(TL % 'DX=1E-3 DPX=-2E-4 DY=5E-4 DPY=1E-4' + 'TRCOD FIXCOD;\n')
    assert '  trcodf' in txt and txt.rstrip().endswith('$$$')
    n = ffs.nlat
    for k in (1, 5, n):
        for key in COD8:
            cur = _tw(ffs, key, k)
            assert _tw(ffs, 'R' + key, k) == cur                    # 基準 = 取った時点の値
            assert _tw(ffs, 'D' + key, k) == 0.0                    # 差は 0
    assert abs(_tw(ffs, 'DX', n)) > 1e-4                             # Twiss["DX"] は絶対値のまま (0 にならない)
    ref = ms.trcod_reference(ffs)
    assert ref.shape == (n, 8) and ms.COD8_KEYS == COD8
    assert ref[n - 1, 4] == _tw(ffs, 'DX', n) and ref[n - 1, 7] == _tw(ffs, 'DPY', n)   # DPY も入る (K64 は DPY を落とす)
    # 他の光学関数の基準は触らない (REFERENCE をしていなければ 0)
    assert _tw(ffs, 'RAX', n) == 0.0


def test_fixcod_difference_equals_analytic_kick_response(v2_mode):
    """基準を取った後でステアリングを Δθ 変える。Twiss["DDX"] (現在 − 基準) = −Δθ·√(β1 β2)·sin(ψ2 − ψ1)"""
    dth = 1e-6
    ffs, I, txt = _run(TL % '' + 'TRCOD FIXCOD;\nST K0 %g;\nCALC;\n' % dth)
    ks = _loc(ffs, 'ST')
    b1, n1 = _tw(ffs, 'BX', ks), _tw(ffs, 'NX', ks)
    a1 = _tw(ffs, 'AX', ks)
    worst = 0.0
    for name in ('MM', 'QF.2', 'QD.2', '$$$'):
        k = _loc(ffs, name)
        b2, n2, a2 = _tw(ffs, 'BX', k), _tw(ffs, 'NX', k), _tw(ffs, 'AX', k)
        dx = -dth * math.sqrt(b1 * b2) * math.sin(n2 - n1)
        dpx = -dth * math.sqrt(b1 / b2) * (math.cos(n2 - n1) - a2 * math.sin(n2 - n1))
        assert _tw(ffs, 'DDX', k) == pytest.approx(dx, rel=1e-5, abs=1e-13)
        assert _tw(ffs, 'DDPX', k) == pytest.approx(dpx, rel=1e-5, abs=1e-13)
        assert _tw(ffs, 'DDX', k) == _tw(ffs, 'DX', k) - _tw(ffs, 'RDX', k)
        worst = max(worst, abs(_tw(ffs, 'DDX', k)))
        assert _tw(ffs, 'DDY', k) == 0.0
    assert worst > 1e-6                                              # 意味のある大きさ
    assert a1 == a1
    # 上流は変わらない
    assert _tw(ffs, 'DDX', ks) == 0.0


def test_fixcod_difference_in_ring_equals_closed_orbit_formula(v2_mode):
    """リング (24 セルの FODO、2 GeV、周長 139.2 m): 基準を取った後の薄いキック Δθ による閉軌道の変化
    Twiss["DDX"] = −Δθ·√(β1 βk)·cos(|ψk − ψ1| − πν)/(2 sin πν)"""
    dth = 1e-5
    src = r'''
MOMENTUM = 2 GEV;
MARK M0 = () MID=();
DRIFT D1 = (L=0.5);
QUAD QF = (L=0.4 K1=0.32) QD=(L=0.4 K1=-0.32);
BEND B1 = (L=1.5 ANGLE=0.1308996938995747) ST=(L=0 ANGLE=0 K0=0);
CAVI CA = (VOLT=5E5 HARM=40);
LINE CELL1 = (QF D1 B1 D1 QD D1 B1 D1);
LINE RING = (M0 12*CELL1 MID ST 12*CELL1 CA);
FFS;
USE RING;
CELL; CALC;
TRCOD FIXCOD;
ST K0 %g;
CALC;
''' % dth
    ffs, I, txt = _run(src)
    k1, n = _loc(ffs, 'ST'), ffs.nlat
    mu = _tw(ffs, 'NX', n)                                           # 2πν
    assert abs(math.sin(mu / 2)) > 0.1                               # 整数・半整数共鳴から離れている
    for k in (1, 20, 50, 97, 150, n):
        ana = -dth * math.sqrt(_tw(ffs, 'BX', k1) * _tw(ffs, 'BX', k)) \
            * math.cos(abs(_tw(ffs, 'NX', k) - _tw(ffs, 'NX', k1)) - mu / 2) / (2 * math.sin(mu / 2))
        assert abs(ana) > 1e-6
        assert _tw(ffs, 'DDX', k) == pytest.approx(ana, rel=1e-4)
        assert _tw(ffs, 'RDX', k) == 0.0 and _tw(ffs, 'DDY', k) == 0.0


def test_fixcod_position_mvcod_stncod(v2_mode):
    src = (TL % 'DX=1E-3 DY=5E-4' + 'TRCOD FIXCOD MM;\nST K0 1E-4;\nCALC;\n')
    ffs, I, txt = _run(src)
    kmm, kq, n = _loc(ffs, 'MM'), _loc(ffs, 'QF.2'), ffs.nlat
    assert f'{kmm:12d} MM' in txt
    assert ffs._measp == kmm                                         # FIXCOD の位置は MEASURE / TRCOD の位置になる (K64 と同じ)
    assert _tw(ffs, 'RDX', kmm) != 0.0 and _tw(ffs, 'RDX', kmm + 1) == 0.0     # 位置より後ろの基準は入らない
    d_mm = _tw(ffs, 'DDX', kmm)
    assert abs(d_mm) > 5e-5
    # MVCOD: 1 点だけ基準を今の値にする (DX DPX DY DPY。EX は変えない)
    ex_ref = _tw(ffs, 'REX', kmm)
    ffs.run_text('TRCOD MVCOD MM;')
    assert _tw(ffs, 'DDX', kmm) == 0.0 and _tw(ffs, 'DDPX', kmm) == 0.0
    assert _tw(ffs, 'REX', kmm) == ex_ref
    assert _tw(ffs, 'RDX', kmm - 1) == pytest.approx(_tw(ffs, 'DX', kmm - 1) - _tw(ffs, 'DDX', kmm - 1))
    # MVCOD の位置の省略は MEASURE / TRCOD の位置。measp は MVCOD では変わらない
    ffs.run_text('TRCOD MVCOD QF.2;')
    assert ffs._measp == kmm and _tw(ffs, 'RDX', kq) == _tw(ffs, 'DX', kq)
    # REFERENCE は参照光学を全部上書きする → STNCOD で FIXCOD / MVCOD の基準 (1..MM) を戻す
    keep = [_tw(ffs, 'RDX', k) for k in range(1, kmm + 1)]
    ffs.run_text('REF;')
    assert _tw(ffs, 'DDX', kmm - 1) == 0.0 and _tw(ffs, 'RAX', n) != 0.0
    ffs.run_text('TRCOD STN;')
    assert [_tw(ffs, 'RDX', k) for k in range(1, kmm + 1)] == keep
    assert _tw(ffs, 'RAX', n) != 0.0                                 # 8 関数以外は触らない


def test_stncod_without_fixcod_and_fixcod_after_trcod(v2_mode):
    ffs, I, txt = _run(TL % 'DX=1E-3 EMITX=1E-16 EMITY=1E-16 DP=1E-5' + 'TRCOD STNCOD;\nNP=400; GAUSS;\nTRCOD;\nTRCOD FIXCOD;\n')
    assert 'no reference orbit is saved' in txt
    n = ffs.nlat
    # TRCOD が入れた重心と分散 (追跡) が基準になる。CALC し直すと差は追跡と CALC の差 (小さい)
    r_dx, r_ex = _tw(ffs, 'RDX', n), _tw(ffs, 'REX', n)
    assert r_dx == _tw(ffs, 'DX', n) and r_ex == _tw(ffs, 'EX', n)
    ffs.run_text('CALC;')
    assert abs(_tw(ffs, 'DDX', n)) < 1e-7 and abs(_tw(ffs, 'DDX', n)) > 0.0
    assert abs(_tw(ffs, 'DEX', n)) < 1e-3


def test_k64_mode_subcommands_unchanged(k64_mode, tmp_path, monkeypatch):
    monkeypatch.setenv('SADPOUT', str(tmp_path / 'p.txt'))
    ffs, I, txt = _run(TL % '' + 'TRCOD FIXCOD;\n')
    assert '未対応' in txt and 'trcodf' not in txt
    assert getattr(ffs, '_ref_rows', None) is None
    assert not ffs.is_word('PUTM') and not ffs.is_word('ORBITJ') and not ffs.is_word('GETM')
    assert not (tmp_path / 'p.txt').exists()


# ============================================================================ TRCOD ALIGN
AL = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 %s) MM=();
DRIFT D1 = (L=1.0);
QUAD QA = (L=0.3 K1=0.25 %s) QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
BEND B1 = (L=1 ANGLE=0.05);
SEXT SX = (L=0.2 K2=3);
OCT OC = (L=0.1 K3=50);
MULT MQ = (L=0.2 K1=0.1 K2=2);
LINE TL = (M0 QA D1 %s D1 QD SX D1 OC MQ D1 MM D1 QF D1);
FFS;
USE TL;
INS; CALC;
'''


@pytest.mark.parametrize('bend, tol', [('QF', 1e-13), ('B1', 2e-9)])
def test_align_equals_shifted_entrance_orbit(v2_mode, bend, tol):
    """bend='QF': 区間は QUAD SEXT OCT MULT だけ (平行移動と厳密に等しい)。bend='B1': BEND を含む (pysad の BEND の DX・DY の扱いで
    ずれの 2 次 dx²·θ ≈ 5e-10 の差が残る)"""
    dx, dy = 1e-4, -2e-4
    fa, I, ta = _run(AL % ('', 'DX=%g DY=%g' % (dx, dy), bend) + 'TRCOD ALIGN QA MM;\nCALC;\n')
    fb, I, tb = _run(AL % ('DX=%g DY=%g' % (-dx, -dy), '', bend))
    # 区間の磁石 (BEND QUAD SEXT OCT MULT) が先頭の QA と同じずれを持つ
    for name in ('QF.1' if bend == 'QF' else bend, 'QD', 'SX', 'OC', 'MQ'):
        el = fa.els[_loc(fa, name) - 1]
        assert (el.DX, el.DY) == (dx, dy), name
    assert fa.els[fa.nlat - 3].name == 'QF' and fa.els[fa.nlat - 3].DX == 0.0     # 区間の外 (最後の QF) は動かさない
    assert fa._align_last[:2] == (_loc(fa, 'QA'), _loc(fa, 'MM')) and fa._align_last[4] == 5
    k = _loc(fa, 'MM')
    got = [_tw(fa, key, k) for key in ('DX', 'DPX', 'DY', 'DPY')]
    ref = [_tw(fb, 'DX', k) + dx, _tw(fb, 'DPX', k), _tw(fb, 'DY', k) + dy, _tw(fb, 'DPY', k)]
    assert max(abs(v) for v in got) > 5e-5
    assert got == pytest.approx(ref, abs=tol)
    # ALIGN をしないと (QA だけがずれている) 軌道は違う
    fc, I, tc = _run(AL % ('', 'DX=%g DY=%g' % (dx, dy), bend))
    assert abs(_tw(fc, 'DX', k) - got[0]) > 1e-6


def test_align_order_groups_and_messages(v2_mode):
    src = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3) MM=();
DRIFT D1 = (L=1.0);
QUAD QA = (L=0.15 K1=0.12 DX=3E-4 DY=1E-4) QD=(L=0.3 K1=-0.25) QE=(L=0.3 K1=0.2);
MONI PM = ();
SOL S1 = (BZ=0 DX=7E-4);
LINE TL = (M0 QA QA D1 QD QD PM D1 MM QE D1 S1 D1);
FFS;
USE TL;
INS; CALC;
'''
    # 逆順 (入れ替える)。is は分割の先頭 QA.1、ie は QD.1 の群の末尾 QD.2 まで広がる
    ffs, I, txt = _run(src + 'TRCOD ALIGN QD QA.2;\n')
    assert ffs._align_last[:2] == (_loc(ffs, 'QA.1'), _loc(ffs, 'QD.2'))
    assert [ffs.els[_loc(ffs, n) - 1].DX for n in ('QA.2', 'QD.1', 'QD.2', 'QE')] == [3e-4, 3e-4, 3e-4, 0.0]
    assert ffs.els[_loc(ffs, 'QD.2') - 1].DY == 1e-4
    # 第 2 の要素が無い・見つからない: 終端まで (K64 と同じ)。SOL と MONI は動かさない
    ffs, I, txt = _run(src + 'TRCOD ALIGN QA;\n')
    assert 'aligned to the end of the beam line' in txt
    assert ffs.els[_loc(ffs, 'QE') - 1].DX == 3e-4 and ffs.els[_loc(ffs, 'S1') - 1].DX == 7e-4
    assert float(getattr(ffs.els[_loc(ffs, 'PM') - 1], 'DX', 0.0) or 0.0) == 0.0
    ffs, I, txt = _run(src + 'TRCOD ALIGN QA NOSUCH;\n')
    assert 'NOSUCH is not found' in txt and ffs.els[_loc(ffs, 'QE') - 1].DX == 3e-4
    # 先頭が見つからない: 何もしない。先頭が磁石でない: 0 に揃える (K64 と同じ)
    ffs, I, txt = _run(src + 'TRCOD ALIGN NOSUCH QE;\nTRCOD ALIGN;\n')
    assert 'undefined element NOSUCH' in txt and 'missing elements' in txt
    assert ffs.els[_loc(ffs, 'QD.1') - 1].DX == 0.0
    ffs, I, txt = _run(src + 'TRCOD ALIGN M0 QD;\n')
    assert [ffs.els[_loc(ffs, n) - 1].DX for n in ('QA.1', 'QA.2', 'QD.2', 'QE')] == [0.0, 0.0, 0.0, 0.0]


# ============================================================================ ORBITJITTER
def test_orbitjitter_adds_gauss_sequence(v2_mode):
    ffs, I, txt = _run(TL % 'DX=1E-3 DY=5E-4')
    n = ffs.nlat
    before = [(_tw(ffs, 'DX', k), _tw(ffs, 'DPX', k), _tw(ffs, 'DY', k), _tw(ffs, 'DPY', k)) for k in range(1, n + 1)]
    rng = ffs.rng
    ind = SadRandom(17)
    ind.seed, ind.filled, ind.next = rng.seed, rng.filled, rng.next       # 同じ状態から独立に乱数列を作る
    xj, yj = 1e-5, 3e-5
    ffs.run_text('ORBITJ %g %g;' % (xj, yj))
    gx, gy = [], []
    for k in range(n):
        gx.append(ind.gauss(1, ffs.gcut()))
        gy.append(ind.gauss(1, ffs.gcut()))
    for k in range(1, n + 1):
        b = before[k - 1]
        assert _tw(ffs, 'DX', k) - b[0] == pytest.approx(xj * gx[k - 1], abs=1e-18)
        assert _tw(ffs, 'DY', k) - b[2] == pytest.approx(yj * gy[k - 1], abs=1e-18)
        assert _tw(ffs, 'DPX', k) == b[1] and _tw(ffs, 'DPY', k) == b[3]
    assert np.std(gx) > 0.4                                          # 乱数が実際に入っている
    assert (rng.seed, rng.next) == (ind.seed, ind.next)             # 乱数列は 2·nlat 個進む
    # 次の CALC で元に戻る (K64 と同じ: Twiss を書き換えるだけ)
    ffs.run_text('CALC;')
    assert _tw(ffs, 'DX', n) == before[n - 1][0]


def test_orbitjitter_gcut_statistics_and_missing_parameter(v2_mode):
    src = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3);
DRIFT D1 = (L=0.01);
LINE TL = (M0 4000*D1);
FFS;
USE TL;
INS; CALC;
GCUT=1.5;
ORBITJITTER 2E-5 0;
'''
    ffs, I, txt = _run(src)
    v = np.array([_tw(ffs, 'DX', k) for k in range(1, ffs.nlat + 1)])
    y = np.array([_tw(ffs, 'DY', k) for k in range(1, ffs.nlat + 1)])
    assert np.abs(v).max() <= 1.5 * 2e-5 and np.abs(v).max() > 1.3 * 2e-5        # GCUT で切れる
    # 1.5σ で切った正規分布の rms = σ·√(1 − 2a·φ(a)/(2Φ(a) − 1))、a = 1.5 → 0.7437σ
    a = 1.5
    phi = math.exp(-a * a / 2) / math.sqrt(2 * math.pi)
    rms = math.sqrt(1 - 2 * a * phi / math.erf(a / math.sqrt(2)))
    assert v.std() == pytest.approx(rms * 2e-5, rel=0.04)
    assert abs(v.mean()) < 4 * rms * 2e-5 / math.sqrt(len(v))
    assert np.all(y == 0.0)
    ffs.run_text('ORBITJITTER 1E-5; CALC;')
    ffs2, I, txt2 = _run(src + 'ORBITJITTER 1E-5;\n')
    assert '?Missing parameter for ORBITJ_ITTER x y' in txt2


# ============================================================================ PUTMEASURE / GETMEASURE
PG = TL % 'DX=1E-3 DY=5E-4 EMITX=1E-9 EMITY=1E-9 DP=1E-4' + 'NP=300; GAUSS;\n'


def test_putmeasure_getmeasure_round_trip(v2_mode, tmp_path, monkeypatch, capfd):
    pf = tmp_path / 'meas.txt'
    monkeypatch.setenv('SADPOUT', str(pf))
    monkeypatch.setenv('SADPIN', str(pf))
    ffs, I, txt = _run(PG + 'MEASURE MM;\nPUTMEASURE;\nMEASURE $$$;\nPUTM;\n')
    try:
        sv_end = list(ffs._measure_sv)                              # xa ya xxa yya xya (終端)
        lines = pf.read_text().splitlines()
        assert len(lines) == 2 and all(len(ln.split(',')) == 5 and len(ln) == 15 * 5 + 4 for ln in lines)
        r1 = ffs.fn_ffs(['GETMEASURE'])
        r2 = ffs.fn_ffs(['GETM'])
        want = [sv_end[0], sv_end[1], sv_end[2], sv_end[4], sv_end[3]]          # ファイルの順は xa ya xxa xya yya
        assert r2 == pytest.approx(want, rel=2e-7)
        assert r1 != r2 and r1 == [float(t) for t in lines[0].split(',')]       # 1 行目は MM での値 (呼ぶたびに次の行)
        err = capfd.readouterr().err
        assert 'Putting data:' in err and 'data gotten from :' in err and 'opening pipe:' in err
        # データが尽きた: 文を出し、目標値は変えない
        out = io.StringIO()
        ffs.out = out
        ffs.run_text('GETM;')
        assert 'no more data' in out.getvalue() and ffs._getm_last is None
    finally:
        ms.close_measure_files(ffs)


def test_getmeasure_sets_fit_goals(v2_mode, tmp_path, monkeypatch):
    """tfgetm.f: 目標 = (xa1 − xa)/2·α + 今の目標·(1 − α)。MEASURE の位置 measp の 1 点の DX・DY の条件だけ"""
    pf = tmp_path / 'other.txt'
    xa1, ya1 = 2.5e-3, -1.5e-3
    pf.write_text('%15.7g,%15.7g,%15.7g,%15.7g,%15.7g\n' % (xa1, ya1, 1e-8, 0.0, 2e-9) * 2)
    monkeypatch.setenv('SADPIN', str(pf))
    ffs, I, txt = _run(PG + 'MEASURE MM;\nFIT MM; DX 1E-4; DY 2E-4; DPX 3E-5; BX 10;\nFIT QD; DX 4E-4;\nFIT MM $$$; DY 6E-4;\n')
    try:
        xa, ya = ffs._measure_sv[0], ffs._measure_sv[1]
        kmm = _loc(ffs, 'MM')

        def goal(kf, kp, kp1):
            c = [c for c in ffs.matcher_().conditions if c['kf'] == kf and c['kp'] == kp and c['kp1'] == kp1]
            assert len(c) == 1 and c[0]['fmin'] == c[0]['fmax']
            return c[0]['fmin']

        ffs.run_text('GETMEASURE 0.25;')
        assert goal('DX', kmm, kmm) == pytest.approx((xa1 - xa) * 0.5 * 0.25 + 1e-4 * 0.75, rel=1e-12)
        assert goal('DY', kmm, kmm) == pytest.approx((ya1 - ya) * 0.5 * 0.25 + 2e-4 * 0.75, rel=1e-12)
        assert goal('DPX', kmm, kmm) == 3e-5 and goal('BX', kmm, kmm) == 10.0     # 他の関数は変えない
        assert goal('DX', _loc(ffs, 'QD'), _loc(ffs, 'QD')) == 4e-4               # 他の位置は変えない
        assert goal('DY', kmm, ffs.nlat) == 6e-4                                  # 2 点の条件は変えない
        ffs.run_text('GETM;')                                                     # alpha の既定は 1
        assert goal('DX', kmm, kmm) == pytest.approx((xa1 - xa) * 0.5, rel=1e-12)
        assert goal('DY', kmm, kmm) == pytest.approx((ya1 - ya) * 0.5, rel=1e-12)
    finally:
        ms.close_measure_files(ffs)


def test_putmeasure_getmeasure_errors(v2_mode, tmp_path, monkeypatch):
    monkeypatch.delenv('SADPOUT', raising=False)
    monkeypatch.delenv('SADPIN', raising=False)
    ffs, I, txt = _run(PG + 'PUTM;\nGETM;\n')
    assert 'SADPOUT' in txt and 'is not set' in txt and 'SADPIN' in txt
    monkeypatch.setenv('SADPIN', str(tmp_path / 'nosuch.txt'))
    monkeypatch.setenv('SADPOUT', str(tmp_path / 'nodir' / 'x.txt'))
    bad = tmp_path / 'bad.txt'
    bad.write_text('1,2,3\n')
    out = io.StringIO()
    ffs.out = out
    ffs.run_text('GETM; PUTM;')
    assert 'unable to open' in out.getvalue() and out.getvalue().count('unable to open') == 2
    monkeypatch.setenv('SADPIN', str(bad))
    ffs.run_text('GETM;')
    ms.close_measure_files(ffs)
    assert '5 comma-separated numbers' in out.getvalue()


# ============================================================================ DETR > 1 の分布 (ttinit)
class _UnitRng:
    """ttinit に決まった乱数を渡す (gauss の呼び出し順: 軌道の揺らぎ 4 個、粒子ごとに z、dp、x、px、y、py)"""

    def __init__(self, seq):
        self.seq, self.i = list(seq), 0

    def gauss(self, n=1, gcut=1e300):
        v = self.seq[self.i] if self.i < len(self.seq) else 0.0
        self.i += 1
        return v

    def uniform0(self):
        raise AssertionError('GAUSS / TRGAUSS なので一様乱数は引かない')


def _ttinit_matrix(mk):
    """ttinit の線形写像 A (4×4): (x, px, y, py) = A·(gx, gpx, gy, gpy)。Σ = A·Aᵀ が設計の分散行列"""
    A = np.zeros((4, 4))
    for j in range(4):
        seq = [0.0] * 6 + [1.0 if k == j else 0.0 for k in range(4)]
        A[:, j] = ms.ttinit(mk, 1, _UnitRng(seq), gauss=True, trgauss=True, jitter=True)[0, :4]
    return A


CPL = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 R1=0.1 R2=0.3 R3=-0.05 R4=0.2 DETR=%s EMITX=1E-8 EMITY=2E-9 DP=0) P1=() P2=() P3=();
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25) QS=(L=0.3 K1=0.9 ROTATE=45) QT=(L=0.3 K1=-0.7 ROTATE=45);
LINE TL = (M0 QF D1 QS P1 D1 QD D1 QT P2 D1 QS D1 QF P3 D1 QD D1);
FFS;
USE TL;
INS; CALC; DP=0;
'''


@pytest.mark.parametrize('detr0', ['0', '1.5'])
def test_ttinit_detr_gt1_consistent_with_calc_coupled_twiss(v2_mode, detr0):
    ffs, I, txt = _run(CPL % detr0)
    mk = ms.mark_values(ffs)
    A0 = _ttinit_matrix(mk)
    S0 = A0 @ A0.T
    J = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], dtype=float)
    em = np.sort(np.abs(np.linalg.eigvals(S0 @ J)))[::2]
    assert em == pytest.approx([2e-9, 1e-8], rel=1e-10)              # 固有エミッタンスは EMITY、EMITX
    sig = np.sqrt(np.diag(S0))
    if detr0 == '1.5':
        # DETR > 1: x モード (EMITX, BX) は主に y 面に出る。R=0 なら σy = √(EMITX·BX)。ここは R≠0 なので大小関係だけ
        assert sig[2] > sig[0]
    else:
        assert sig[0] > sig[2]
    eps = 1e-7
    forms = []
    for name in ('P1', 'P2', 'P3', '$$$'):
        k = _loc(ffs, name)
        cols = [[0.0] * 5 for _ in range(6)] + [[1.0] * 5]
        for j in range(4):
            cols[j][j + 1] = eps
        r = ffs.fn_track([[1.0, cols], float(k)])[1]
        X = np.array(r[:6], dtype=float).T
        M = ((X[1:5, :4] - X[0, :4]) / eps).T                        # 1 粒子の追跡による 4×4 の転送行列
        assert abs(np.linalg.det(M) - 1.0) < 1e-9
        S = M @ S0 @ M.T
        row = ffs.row(k)
        mk2 = dict(mk, AX=row.alfx, BX=row.betx, AY=row.alfy, BY=row.bety, R1=row.r1, R2=row.r2, R3=row.r3, R4=row.r4)
        sc = np.sqrt(np.outer(np.diag(S), np.diag(S)))
        d = []
        for dd in (0.0, 1.5):                                        # その位置の Twiss を「通常」「入れ替え」の形で読む
            mk2['DETR'] = dd
            A1 = _ttinit_matrix(mk2)
            d.append(float(np.abs((A1 @ A1.T - S) / sc).max()))
        assert min(d) < 1e-9 and max(d) > 0.1, (name, d)             # どちらか一方の形で一致する
        forms.append(d.index(min(d)))
    assert set(forms) == {0, 1}                                      # 経路の途中に両方の形が現れる


# ============================================================================ LWAKE / TWAKE ON の MEASURE
WK = r'''
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 DX=1E-3 DY=-5E-4 EMITX=1E-9 EMITY=1E-9 DP=1E-4 SIGZ=3E-3) WK=() MM=();
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
LINE TL = (M0 WK MM QF D1 QD D1);
FFS;
USE TL;
INS; CALC;
NP=%d; GAUSS; NORAD; NOFLUC;
PBUNCH=1E10; NBUNCH=1;
WakeFunction[Longitudinal, "WK"] = %s;
WakeFunction[Transverse, "WK"] = %s;
r0=FFS["MEASURE QF"];
%s
r1=FFS["MEASURE QF"];
'''
ZERO = '{{-1,0},{0,0},{1,0},{2,0}}'


def _wake_run(np_, wl, wt, flags):
    ffs, I, txt = _run(WK % (np_, wl, wt, flags))
    r0 = ffs.I.run('r0', 'script')
    r1 = ffs.I.run('r1', 'script')
    return np.array(r0[1]), np.array(r0[2]), np.array(r1[1]), np.array(r1[2]), ffs


@pytest.mark.parametrize('np_', [600, 1500])
def test_measure_lwake_uniform_wake_centroid(v2_mode, np_):
    """一様な縦ウェイク W0 (V/C): 全粒子が同じ −e·N·W0/(p0c) を受ける。粒子数 NP に依らない"""
    w0, nb, p0 = 1e12, 1e10, 1e9
    c0, s0, c1, s1, ffs = _wake_run(np_, '{{-1,%g},{0,%g},{1,%g},{2,%g}}' % (w0, w0, w0, w0), ZERO, 'LWAKE;')
    d = c1 - c0
    assert d[5] == pytest.approx(-ECHARGE * nb * w0 / p0, rel=1e-9)
    assert np.abs(d[:5]).max() < 1e-15
    assert s1[5, 5] == pytest.approx(s0[5, 5], rel=1e-9)             # 広がりは変わらない
    assert not ffs.flags['TRPT']                                     # v2: MEASURE は TRPT を元に戻す


def test_measure_twake_uniform_wake_centroid(v2_mode):
    """一様な横ウェイク W1 (V/C/m): Δpx = e·N·W1·<x>/(p0c)、Δpy = e·N·W1·<y>/(p0c)"""
    w1, nb, p0 = 1e14, 1e10, 1e9
    c0, s0, c1, s1, ffs = _wake_run(1000, ZERO, '{{-1,%g},{0,%g},{1,%g},{2,%g}}' % (w1, w1, w1, w1), 'LWAKE; TWAKE;')
    d = c1 - c0
    f = ECHARGE * nb * w1 / p0
    assert abs(c0[0]) > 5e-4 and abs(c0[2]) > 2e-4
    assert d[1] == pytest.approx(f * c0[0], rel=1e-9)
    assert d[3] == pytest.approx(f * c0[2], rel=1e-9)
    assert max(abs(d[0]), abs(d[2]), abs(d[4]), abs(d[5])) < 1e-15
    # LWAKE・TWAKE が OFF なら効かない
    c0, s0, c1, s1, ffs = _wake_run(300, ZERO, '{{-1,%g},{0,%g},{1,%g},{2,%g}}' % (w1, w1, w1, w1), '')
    assert np.all(c1 == c0)


def test_measure_lwake_step_wake_uses_z_distribution(v2_mode):
    """階段状の縦ウェイク (片側だけ W0): 平均は −e·N·W0/(2 p0c)、各粒子の損失は前後の粒子の割合に比例するので、
    δ の分散は (e·N·W0/p0c)²/12 だけ増える (向きを逆にした二つの表の平均で、初期の δ との相関が消える)。
    K64 は wake ON で全粒子を z=0 にするので広がりは増えない。v2 は z の分布を作る"""
    w0, nb, p0, n = 1e12, 1e10, 1e9, 2000
    k = ECHARGE * nb * w0 / p0
    ta = '{{-1,%g},{-1E-6,%g},{1E-6,0},{1,0},{2,0}}' % (w0, w0)
    tb = '{{-1,0},{-1E-6,0},{1E-6,%g},{1,%g},{2,%g}}' % (w0, w0, w0)
    va = []
    for t in (ta, tb):
        c0, s0, c1, s1, ffs = _wake_run(n, t, ZERO, 'LWAKE;')
        assert (c1 - c0)[5] == pytest.approx(-k / 2, rel=5e-4)
        va.append(s1[5, 5] - s0[5, 5])
    assert va[0] * va[1] < 0 or min(va) > 0                          # 相関の符号は表の向きで変わる
    assert (va[0] + va[1]) / 2 == pytest.approx(k * k / 12, rel=0.03)
