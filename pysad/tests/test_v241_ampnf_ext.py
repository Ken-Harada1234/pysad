# -*- coding: utf-8 -*-
"""v2.4.1: ampnf の限界を広げる拡張 (`ampnf_ext.py`・`ampnf.py` の `jorder`・`d2nu`・`keep`・`k0_nonlinear`)。
(1) J の次数: J² の項を足すと追跡のチューンとの差が減る (FODO、J¹ のみ → J²)、`d2nu` は対称で c20 の 2 倍、`nu(jorder=)` が級数を切る。
(2) 島 (1 本の共鳴の包絡): 非周期の FODO (24 セルのうち 12 セルに 6 極) の 6 νx = 5 で、理論の固定点の作用と島の幅が追跡の周期 6 の固定点 (O・X) と
    X から出る不安定多様体 (分離線) の作用の範囲に合う。24 回対称の FODO は 6νx = 5 が禁止で駆動項が 1e-3 倍以下 (島が作れない)。
(3) K0≠0 の磁石の非線形を他のオプションと独立に切る `k0_nonlinear`: ν0 は同じ、J¹ の差は磁石自身の非線形だけ。
条件: pysad k64、RF・放射なし、結合なし。"""
import os
import sys

import numpy as np
import pytest

from pysad import ampnf as A
from pysad import ampnf_ext as E
from pysad.parser import load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'normal_form'))
BASE = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
PAT = 'AABABBABAABBBAABABBAABBA'
OPT = dict(bend_cubic=True, kin=True, fringe=True, nslice_bend=8)


def fodo(kf=0.7, kd=-0.6, k2f=4.0, k2d=-6.0, pat=None):
    src = BASE.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g) SF0 = (L=0.2 K2=0) SD0 = (L=0.2 K2=0)' % (k2f, k2d))
    src = src.replace('K1=0.7', 'K1=%r' % kf).replace('K1=-0.6', 'K1=%r' % kd)
    if pat:
        src = src.replace('LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);', 'LINE CELLA = (SF DS QF D2 B D2 QD DS SD D1); LINE CELLB = (SF0 DS QF D2 B D2 QD DS SD0 D1);')
        src = src.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF %s);' % ' '.join('CELL' + c for c in pat))
    return load_sad_string(src).build_lattice('RING')


def test_jorder_truncation_d2nu_and_j2_improves_on_tracking():
    """nu(jorder=) は a+b ≤ jorder の項だけ。d2nu は [[2 c20, c11], [c11, 2 c02]]。追跡 (正規形の J から出発) との差は J² を足すと J¹ のみの 1/2.5 以下 (J = 8e-6・4e-6。残りは J³ と K64 の積分器の分)"""
    lat = fodo()
    r = A.amplitude_tune(lat, D=5, nslice=8, **OPT)
    jx, jy = 8e-6, 4e-6
    assert abs(r.nu(0, jx, jy, 0.0, jorder=1) - (r.nu0[0] + r.adts['xx'] * jx + r.adts['xy'] * jy)) < 1e-14
    assert abs(r.nu(0, jx, jy, 0.0, jorder=2) - r.nu(0, jx, jy, 0.0)) < 1e-14          # D = 5 は J² まで
    d2 = r.d2nu(0)
    assert d2[0, 0] == 2 * r.adts2['x20'] and d2[0, 1] == d2[1, 0] == r.adts2['x11'] and d2[1, 1] == 2 * r.adts2['x02']
    nus = np.array(A.tracked_tunes(lat, r, jx, jy, nturn=1024, momentum_eV=3e9))
    for pl in range(2):
        e1 = abs(nus[pl] - r.nu(pl, jx, jy, 0.0, jorder=1))
        e2 = abs(nus[pl] - r.nu(pl, jx, jy, 0.0, jorder=2))
        assert e2 < e1 / 2.5 and e2 < 6e-6


def test_island_6th_order_theory_matches_tracking():
    """非周期 FODO (12 セルに 6 極)、6νx = 5 の下側 δν = 0.004 (kf = 0.7103193987)。D = 5 の正規形から共鳴の単項式 ζ̄⁵ を残し、固定点の作用と島の幅 (振り子)
    を出す。追跡の周期 6 の固定点 (O・X) の作用の平均と 2%、X から出る不安定多様体の作用の範囲 [J_lo, J_hi] と 5% 以内、島の幅は J_r の 30% 以下 (振り子が成り立つ範囲)"""
    lat = fodo(kf=0.7103193987084756, pat=PAT)
    r = A.amplitude_tune(lat, D=5, nslice=8, **OPT)
    al = E.pure_alpha(6, 0)
    th = E.island_theory(r, E.resonance_normal_form(r.res, 0, al), 0, al)
    assert th['ok'] and th['pendulum_ok'] and th['k'][0] == -6 and th['p'] == -5
    tr = E.track_island(lat, r, 6, 0, 3e9, jr=th['Jr'][0], nturn_sep=300)
    assert tr['n_O'] == 6 and tr['n_X'] == 6 and not tr['escaped']
    jo = np.mean([j for _, j, _ in tr['O']])
    jx = np.mean([j for _, j, _ in tr['X']])
    assert abs(th['Jr'][0] / (0.5 * (jo + jx)) - 1) < 0.02
    assert abs(th['Jlo'] / tr['Jmin'] - 1) < 0.05 and abs(th['Jhi'] / tr['Jmax'] - 1) < 0.05


def test_symmetry_forbidden_resonance_has_no_drive():
    """24 回対称の FODO では 6νx = 5 (24 の倍数でない) は禁止: 駆動項 |r| は、同じ強さの 6 極を 12 セルに置いた非周期の配置の 1e-3 倍以下"""
    out = []
    for pat in (None, PAT):
        lat = fodo(kf=0.7103193987084756, pat=pat)
        r = A.amplitude_tune(lat, D=5, nslice=4)
        tb = E.drive_table(r.res, orders=(6,), planes=(0,))
        out.append(tb[(0, 6)]['abs'])
    assert out[0] < 1e-3 * out[1]


def test_island_theory_reports_unreachable_resonance():
    """作る振幅の範囲 (ν(J) の級数が jmax まで) に共鳴がなければ ok = False と理由を返す"""
    lat = fodo()
    r = A.amplitude_tune(lat, D=4, nslice=4)
    al = E.pure_alpha(3, 0)
    th = E.island_theory(r, E.resonance_normal_form(r.res, 0, al), 0, al, jmax=1e-7)
    assert not th['ok'] and 'reason' in th


def test_k0_nonlinear_switch_is_independent_of_other_options():
    """K0≠0 の磁石 (3 極の挿入光源) の非線形だけを切る: ν0 は同じ (線形の光学は残る)、J¹ の差は磁石自身の非線形の分 (xx・yy で 0.02〜0.1%)。
    既定 (None) は従来どおり: 3 つのオプションが全部 False のときだけ 1 次 (k0_nonlinear=False と同じ)、1 つでも True なら全部入り (True と同じ)"""
    import k0_ampnf as G
    lat = G.build(a=0.06)
    opt = dict(D=4, nslice=8, nslice_bend=4)
    full = A.amplitude_tune(lat, bend_cubic=True, kin=True, fringe=True, **opt)
    lin = A.amplitude_tune(lat, bend_cubic=True, kin=True, fringe=True, k0_nonlinear=False, **opt)
    assert abs(full.nu0[0] - lin.nu0[0]) < 1e-12 and abs(full.nu0[1] - lin.nu0[1]) < 1e-12
    dyy = full.adts['yy'] - lin.adts['yy']
    assert 0.02 < abs(dyy) < 1.0 and abs(dyy) / abs(full.adts['yy']) < 1e-2
    legacy_none = A.amplitude_tune(lat, **opt)
    forced_lin = A.amplitude_tune(lat, k0_nonlinear=False, **opt)
    assert abs(legacy_none.adts['yy'] - forced_lin.adts['yy']) < 1e-12
    legacy_bc = A.amplitude_tune(lat, bend_cubic=True, **opt)
    forced_full = A.amplitude_tune(lat, bend_cubic=True, k0_nonlinear=True, **opt)
    assert abs(legacy_bc.adts['yy'] - forced_full.adts['yy']) < 1e-12
