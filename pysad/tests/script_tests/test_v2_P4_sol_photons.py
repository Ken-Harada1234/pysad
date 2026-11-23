# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 2: SOL 領域内の光子 (PhotonList、PHOTONS ON + RAD + FLUC の追跡)。

K64 (gs、2026-10-19 ビルド) は領域内の MULT・DRIFT の光子を記録する (3 GeV、BZ = 1 T、MULT K0 = 0.01: 光子の gx は MULT の範囲 0.5〜0.9 m に入り、
エネルギーの合計は古典値 285 eV の 1.1 倍 (統計の範囲))が、領域内の QUAD の光子は記録しない (RAD ON・BZ≠0 の領域内の四極の放射がない: B68 の疑い)。領域内の BEND は
偏向をほぼ無視する (B69)。pysad の k64 は領域内の MULT・BEND が未対応 (NotImplementedError)、DRIFT の光子は領域の始点の位置 (gx = 0) で記録する (以前の動作。位置は
意味がない)、QUAD は記録しない。v2 は領域内の DRIFT・QUAD・MULT の光子を、領域の外の同じ要素と同じ位置の情報 (`_ph_set`・刻みごとの `_ph_adv`) で記録する
(BEND は未対応: 領域の BEND は直線の座標で一様な双極磁場の蹴り (ANGLE + K0) を受けるので、ANGLE の円弧の幾何とは別の軌道になる)。

根拠 (gs は使わない): (1) 光子の位置が粒子の軌道上にある (MULT の双極 K0 = 0.01 の放物線 x = (K0/2L) s²、光子の座標は SAD の大域座標 (gx = 進行方向、gy = x、gz = y))、
(2) 光子のエネルギーの合計が古典的な放射損失 (C_γ E⁴ θ/(2πρ))・RAD NOFLUC の損失と合う (3 GeV、K0 = 0.01 (θ = 0.01 rad、ρ = 40 m): 285 eV)、
(3) 光子の要素番号 (l_track) が領域内の要素。使った条件: 3 GeV、BZ = 1 T (bz = 0.1 1/m)、MULT (L = 0.4、K0 = 0.01)・QUAD (L = 0.4、K1 = 0.8)、1500 粒子、FLUC の統計
(光子エネルギーの合計の標準誤差は数 %)。未確認: 傾いた・ずれた SOL の領域内の光子の位置 (geo_frac は要素の幾何をそのまま使う)、領域内の BEND・SEXT・OCT、
偏光 (ξ) の値、大きなリング。
"""
import io

import numpy as np
import pytest

_SRC = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=0.5);
MULT M = (L=0.4 K0=0.01);
QUAD Q = (L=0.4 K1=0.8);
SOL SA = (L=0 BZ=%(bz)s BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST=();
LINE TS = (ST SA D1 M D1 Q SB);
LINE T0 = (ST D1 M D1 Q);
LINE TSQ = (ST SA D1 Q D1 SB);
FFS;
NPARA=1;
USE %(line)s; INS; CALC;
RAD; FLUC; PHOTONS;
SeedRandom[5];
n=%(n)d;
x0=Table[%(x0)s,{i,n}]; px0=Table[%(x0)s,{i,n}]; y0=Table[0,{i,n}];
py0=Table[0,{i,n}]; z0=Table[0.0,{i,n}]; dp0=Table[0.0,{i,n}]; fl0=Table[1,{i,n}];
r=TrackParticles[{1,{x0,px0,y0,py0,z0,dp0,fl0}},1,1];
Exit;
'''


def _photons(line='TS', bz='1', n=1500, x0='0'):
    from pysad.ffscmd import run_sad_string
    ffs, I = run_sad_string(_SRC % dict(line=line, bz=bz, n=n, x0=x0), out=io.StringIO())
    return np.array(I.globals['PhotonList'], dtype=float).reshape(-1, 12) if len(I.globals.get('PhotonList', [])) else np.zeros((0, 12))


_E_MULT = 88.46e3 * 3.0 ** 4 * 0.01 / (2 * np.pi * 40.0)       # eV: C_γ E⁴ θ/(2π ρ)、ρ = L/θ = 40 m


def test_photons_in_sol_region_are_on_the_orbit_and_in_the_element_v2(v2_mode):
    a = _photons()
    assert len(a) > 1000
    m = a[:, 11] == 4                                           # MULT (ST = 1、SA = 2、D1 = 3、M = 4、D1 = 5、Q = 6)
    s = a[m, 1] - 0.5
    assert m.sum() > 500 and s.min() >= -1e-9 and s.max() <= 0.4 + 1e-9
    assert np.abs(a[m, 2] - 0.0125 * s ** 2).max() < 5e-7       # 放射点は軌道 x = (K0/2L) s² の上 (K0 > 0 は大域座標の +gy。最大 2e-3 に対し 2e-7 はソレノイドのらせんと反跳)
    assert np.abs(a[m, 3]).max() < 1e-4                         # y はほぼ 0 (ソレノイドのらせんの分だけ)
    for l, lo, hi in ((5, 0.9, 1.4), (6, 1.4, 1.8)):            # 領域内の DRIFT と QUAD (放射する)
        mm = a[:, 11] == l
        assert mm.sum() > 5 and a[mm, 1].min() >= lo - 1e-9 and a[mm, 1].max() <= hi + 1e-9
    assert set(a[:, 11]) <= {4.0, 5.0, 6.0}


def test_photon_energy_sum_in_sol_region_equals_classical_loss_v2(v2_mode):
    a = _photons()
    m = a[:, 11] == 4
    e_per = a[m, 0].sum() / 1500.0
    assert e_per == pytest.approx(_E_MULT, rel=0.12), (e_per, _E_MULT)           # 光子エネルギーの合計の統計誤差は約 5〜8%
    # 領域の外の同じ並び (bz = 0) の MULT: 同じ古典値と同じ精度で合う (領域内の光子の記録が領域の外と同じ式を使う)
    b = _photons(line='T0')
    mb = b[:, 11] == 3
    assert b[mb, 0].sum() / 1500.0 == pytest.approx(_E_MULT, rel=0.12)


def test_photons_from_quad_in_sol_region_are_recorded_in_v2(v2_mode):
    """x = 2 mm の粒子が領域内の QUAD (K1 = 0.8) で放射する光子: v2 は記録する (要素番号 4 = Q。位置は QUAD の範囲)。k64 は記録しない (下)"""
    a = _photons(line='TSQ', x0='2e-3', n=800)
    assert len(a) > 100 and set(a[:, 11]) <= {3.0, 4.0, 5.0} and a[a[:, 11] == 4][:, 1].min() >= 0.5 - 1e-9 and a[a[:, 11] == 4][:, 1].max() <= 0.9 + 1e-9


def test_k64_does_not_record_quad_photons_in_sol_region(k64_mode):
    """k64: 領域内の QUAD (要素番号 4) の光子は記録しない (今までどおり乱数だけ引く)。領域内の DRIFT の光子は、K64 の `tsetpcvt` を直訳した古い経路が
    領域の始点 (SOL、要素番号 2、gx = 0) の位置で記録する (今までの動作で、位置は意味がない。gs は別の古い位置の情報で記録する)"""
    a = _photons(line='TSQ', x0='2e-3', n=800)
    assert 4.0 not in set(a[:, 11]) and 3.0 not in set(a[:, 11])
