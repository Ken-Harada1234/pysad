# -*- coding: utf-8 -*-
"""v2.3.0 T1: elegant との比較 (局所運動量アクセプタンス・Touschek 寿命)。

elegant を使わない部分: SAD→elegant の変換 (要素の対応、K2 = −K1·ANGLE/L)、保存済みの elegant の結果 (example/15_lma_touschek/touschek_elegant/elegant_*) と pysad の結果 (Q1・T1 の表) の比較、
K1≠0 の BEND のオフモメンタムの光学 (解析式と elegant・gs の保存値)。elegant を使う部分 (WSL の elegant で ASTRID2 の光学を計算して pysad と比べる) は、elegant が無い環境では skip。
"""
import json
import math
import os
import re
import subprocess
import sys

import numpy as np
import pytest

from pysad import compat

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
EX = os.path.join(ROOT, 'example', '15_lma_touschek', 'touschek_elegant')
sys.path.insert(0, EX)


@pytest.fixture(autouse=True)
def _v2():
    old = compat._mode
    compat.set_mode('v2')
    yield
    compat.set_mode(old)


def _read_ele(name):
    rows = {}
    for ln in open(os.path.join(EX, 'elegant_%s.lma.txt' % name), encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        rows[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
    return rows


def _read_pysad(path):
    rows = {}
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        rows[int(t[0])] = (float(t[1]), float(t[3]), float(t[4]))
    return rows


def test_converter_astrid2(tmp_path):
    """変換: 要素の数・種類、複合機能の BEND の K2 = −K1·ANGLE/L、RF の位相、各要素の出口の MAXAMP"""
    import sad2ele as c
    path, st, ffs = c.convert('astrid2', str(tmp_path), 6166.703905, aper=0.010)
    txt = open(path, encoding='utf-8').read()
    defs = dict(re.findall(r'^(\w+): (.*)$', txt, re.M))
    line = txt.split('RING: LINE=(')[1].replace('&\n', '').rstrip(')\n')
    items = [x.strip() for x in line.split(',')]
    assert len(items) == 2 * len(st.els)                      # 要素ごとに MAXAMP が付く
    assert sum(1 for x in items if x == 'AP') == len(st.els)
    bh = defs['B_BH']
    kv = dict(re.findall(r'(\w+)=([-0-9.e+]+)', bh))
    assert abs(float(kv['K1']) + 1.6782) < 1e-12
    assert abs(float(kv['K2']) - 1.6782 * 0.523599 / 0.85) < 1e-9          # −K1·ANGLE/L
    assert kv['EDGE_ORDER'] == '2' and float(kv['E1']) == 0.0
    rf = dict(re.findall(r'(\w+)=([-0-9.e+]+)', defs['RF_CAVI']))
    assert abs(float(rf['VOLT']) - 5e4) < 1e-6 and abs(float(rf['FREQ']) - 104.95e6) < 1
    assert abs(float(rf['PHASE']) - (180 - math.degrees(math.asin(6166.703905 / 5e4)))) < 1e-6
    assert defs['M_PV'] == 'MARK'                                # 長さ 0・ANGLE 0 の BEND は恒等
    assert sum(1 for k in defs if k.startswith('B_')) == 1
    path2, _, _ = c.convert('astrid2', str(tmp_path), 6166.703905, sad_k1=False, tag='_exact')
    assert 'K2=0,' in open(path2, encoding='utf-8').read().split('B_BH:')[1].split('\n')[0]


def test_converter_kekls(tmp_path):
    import sad2ele as c
    path, st, ffs = c.convert('kekls', str(tmp_path), 296074.88, aper=0.010)
    txt = open(path, encoding='utf-8').read()
    defs = dict(re.findall(r'^(\w+): (.*)$', txt, re.M))
    for k in ('B_BPI6', 'B_BPI7'):
        kv = dict(re.findall(r'(\w+)=([-0-9.e+]+)', defs[k]))
        assert abs(float(kv['K2']) + float(kv['K1']) * float(kv['ANGLE']) / float(kv['L'])) < 1e-9 * abs(float(kv['K2']))
    kv = dict(re.findall(r'(\w+)=([-0-9.e+]+)', defs['B_BID5']))
    assert float(kv['K2']) == 0.0                              # 純粋な偏向磁石には K2 を足さない
    assert 'KOCT' in defs['O_OCFA'] and 'K3=-6000' in defs['O_OCFA']
    rf = dict(re.findall(r'(\w+)=([-0-9.e+]+)', defs['RF_RFC']))
    assert abs(float(rf['FREQ']) - 493092781.43) < 1.0


def test_lma_astrid2_saved_vs_pysad():
    """保存済みの elegant (全周 199 位置) と pysad (Q1、1 セル 25 位置): 同じ位置の δ± が elegant の最終刻み + pysad の tol の範囲で一致。
    分散の大きい位置はアパーチャが限界 (δ ≈ 0.83%)、直線部は RF バケット (δ_rf = 1.618%) で、どちらも相対差 0.5% 以内"""
    ele = _read_ele('astrid2')
    py = _read_pysad(os.path.join(EX, '..', 'lma', 'lma_astrid2.table.txt'))
    assert len(ele) == 199 and len(py) == 25
    for l, (s, dp, dm) in py.items():
        e = ele[l]
        assert abs(e[0] - s) < 1e-6
        assert abs(e[1] - dp) < 1.2e-4 and abs(e[2] - dm) < 1.2e-4
        assert abs(e[1] / dp - 1) < 0.005 and abs(e[2] / dm - 1) < 0.005
    dps = np.array([v[1] for v in py.values()])
    assert 0.0082 < dps.min() < 0.0084 and 0.0159 < dps.max() < 0.0162
    # elegant の全周: 6 セルが同じ (周期性。空洞の分だけ要素の番号がずれるので find_periodicity の starts を使う)
    from pysad.trksetup import TrackSetup, find_periodicity
    st = TrackSetup.from_source(open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read(), 'RING',
                                rad=False, rfsw=False)[0]
    per = find_periodicity(st.els, st.rows)
    assert per['nper'] == 6
    n = 0
    for k in range(2, 7):
        for l in range(per['starts'][0], per['starts'][1]):
            lk = per['starts'][k - 1] + (l - per['starts'][0])
            if lk in ele and l in ele:
                n += 1
                assert ele[lk][1] == ele[l][1] and ele[lk][2] == ele[l][2]
    assert n >= 150


def test_lma_kekls_saved_vs_pysad():
    """保存済みの elegant (1500 周、14 位置) と pysad (check_every = 1、491 周、52 位置) の kek-ls_v03_12.lat (3 GeV、20 セル): 共通の 11〜12 位置の δ± が
    elegant の刻み 7.8e-5 + pysad の tol 1e-4 の範囲で一致し、全部 RF バケット (δ_rf = 2.994% の 1.00〜1.01 倍) で決まる"""
    ele = _read_ele('kekls')
    py = _read_pysad(os.path.join(EX, 'pysad_kekls.lma.txt'))
    assert len(ele) == 14 and len(py) == 52
    common = sorted(set(ele) & set(py))
    assert len(common) >= 11
    bh = 0.02993977
    for l in common:
        assert abs(ele[l][0] - py[l][0]) < 1e-6
        assert abs(ele[l][1] - py[l][1]) < 1.3e-4 and abs(ele[l][2] - py[l][2]) < 1.3e-4
        assert abs(ele[l][1] / py[l][1] - 1) < 0.005 and abs(ele[l][2] / py[l][2] - 1) < 0.005
        for v in (ele[l][1], -ele[l][2], py[l][1], -py[l][2]):
            assert 0.99 * bh < v < 1.02 * bh                     # RF バケットが限界 (アパーチャではない)
    # elegant の周期性: 9 番目のセルの 3 点と最初のセルの差は刻み 1 つ以内
    others = sorted(set(ele) - set(common))
    assert 2 <= len(others) <= 3
    for l in others:
        assert abs(ele[l][1] - 0.030156) < 1.6e-4 and abs(ele[l][2] + 0.030156) < 1.6e-4


def test_touschek_saved_astrid2():
    """保存済みの寿命: 一様な受容 (表の点 δ = 0.4・0.6・1.0%、δ_rf) の elegant / pysad が ±5% 以内、局所の LMA は 5% 以内 (elegant の式は 2 つの積分法で 3e-5 一致)"""
    ref = json.load(open(os.path.join(EX, 'pysad_astrid2.ref.json')))
    ele = {}
    for ln in open(os.path.join(EX, 'elegant_astrid2.touschek.txt'), encoding='utf-8'):
        if ln.startswith('#'):
            continue
        lb, rest = ln.split(':', 1)
        m = re.search(r'tLifetime \(hour\) =\s*([-0-9.e+]+)', rest)
        ele[lb.strip()] = float(m.group(1)) * 3600.0
    for a in ('0.004', '0.006', '0.01', '0.01617565869'):
        t_py = ref['pysad_uniform'][a]
        key = [k for k in ele if k.startswith('uniform') and abs(float(k.split()[1]) - float(a)) < 6e-5]
        assert key, a
        r = ele[key[0]] / t_py
        assert 0.94 < r < 1.06, (a, r)
    assert abs(ele['elegant LMA'] / ele['pysad LMA'] - 1) < 0.01            # 同じ式で 2 つの LMA を比べる: LMA の違いは寿命で 1% 以内
    assert 0.93 < ele['pysad LMA'] / ref['pysad_local'] < 1.02              # pysad の局所の LMA の寿命 (Q1) と elegant の式の差


def test_bend_k1_offmomentum_formula():
    """K1≠0 の BEND (複合機能) のオフモメンタムの光学: 一様な磁場指数 n = 0.4 の 16 台の扇形磁石。pysad (K64 と同じ) の dν/dδ は、勾配の力に曲率の因子 (1 + h x) を
    掛けない模型の解析式 (−0.1721、+0.2108) に一致し、K1 が一定の勾配とする厳密な曲線座標の解析式 (−0.6025、+0.7379。elegant はこちらに 5 桁一致) とは違う。
    elegant の CSBEND に K2 = −h K1 を足すと pysad に一致する (elegant の値は example/15_lma_touschek/touschek_elegant/gs_bend_k1.out.txt に記録)"""
    from pysad import resdrive
    from pysad.trksetup import TrackSetup
    import gs_bend_k1 as g
    src = g.lattice()[0]
    st, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=1.0)
    d1 = resdrive.chromaticity_numeric(ffs.lat)['d1']
    n = g.NIDX
    nux, nuy = math.sqrt(1 - n), math.sqrt(n)
    sad = (-n * n / (2 * nux * (1 - n)), nuy * n / (2 * (1 - n)))
    exact = (-n * (1 + n) / (2 * nux * (1 - n)), nuy * (1 + n) / (2 * (1 - n)))
    assert abs(d1[0] - sad[0]) < 2e-3 and abs(d1[1] - sad[1]) < 2e-3
    assert abs(d1[0] - exact[0]) > 0.3 and abs(d1[1] - exact[1]) > 0.3


def _elegant_ok():
    try:
        p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'source ~/elegant/env.sh; which elegant'], capture_output=True, text=True,
                           timeout=60, env=dict(os.environ, MSYS_NO_PATHCONV='1'))
        return p.returncode == 0 and 'elegant' in p.stdout
    except Exception:
        return False


@pytest.mark.skipif(not _elegant_ok(), reason='WSL の elegant が無い')
def test_elegant_optics_astrid2(tmp_path):
    """elegant (WSL) で変換した ASTRID2 の光学が pysad と合う: チューン 1e-5、αc 2e-4、自然エミッタンス 1e-4、エネルギー拡がり 1e-4、U0 1e-4、色収差 0.05 以内"""
    import sad2ele as c
    ref = json.load(open(os.path.join(EX, 'pysad_astrid2.ref.json')))
    path, st, ffs = c.convert('astrid2', str(tmp_path), ref['U0_eV'], aper=0.010, maxamp=False)
    lte = os.path.basename(path)
    ele = ("&run_setup lattice = %s, use_beamline = RING, p_central_mev = 580.0, default_order = 3, rootname = t, &end\n&run_control &end\n"
           "&twiss_output filename = %%s.twi, radiation_integrals = 1, higher_order_chromaticity = 1, &end\n&bunched_beam n_particles_per_bunch = 1 &end\n&track &end\n") % lte
    open(os.path.join(str(tmp_path), 't.ele'), 'w').write(ele)

    def wp(p):
        p = os.path.abspath(p).replace(os.sep, '/')
        return '/mnt/' + p[0].lower() + p[2:]
    cmd = ('source ~/elegant/env.sh; rm -rf ~/t1/pytest_a2; mkdir -p ~/t1/pytest_a2; cd ~/t1/pytest_a2; cp "%s" "%s" .; elegant t.ele > log.txt 2>&1; '
           'sddsprintout t.twi -par=nux -par=nuy -par=dnux/dp -par=dnuy/dp -par=alphac -par=ex0 -par=Sdelta0 -par=U0') % (
        wp(path), wp(os.path.join(str(tmp_path), 't.ele')))
    p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=300,
                       env=dict(os.environ, MSYS_NO_PATHCONV='1'))
    vals = {}
    for k, v in re.findall(r'(nux|nuy|dnux/dp|dnuy/dp|ex0|U0)\s*\([^=]*\)\s*=\s*([-0-9.e+]+)', p.stdout):
        vals[k] = float(v)
    for k, v in re.findall(r'(alphac|Sdelta0)\s*=\s*([-0-9.e+]+)', p.stdout):
        vals[k] = float(v)
    assert abs(vals['nux'] - ref['tune'][0]) < 1e-5 and abs(vals['nuy'] - ref['tune'][1]) < 1e-5
    assert abs(vals['alphac'] / ref['alpha_c'] - 1) < 2e-4
    assert abs(vals['ex0'] / ref['natural']['EMITX'] - 1) < 1e-4
    assert abs(vals['Sdelta0'] / ref['natural']['SIGE'] - 1) < 1e-4
    assert abs(vals['dnux/dp'] - ref['chrom'][0]) < 0.05 and abs(vals['dnuy/dp'] - ref['chrom'][1]) < 0.05
    assert abs(vals['U0'] * 1e6 / ref['U0_eV'] - 1) < 1e-4


def test_touschek_saved_kekls():
    """kek-ls (3 GeV): 一様な受容 (表の点 0.4・0.6・1.0・2.0%) の elegant / pysad が ±4% 以内、δ_rf で +3% 以内、局所の LMA の寿命は pysad の式と +3% 以内"""
    ref = json.load(open(os.path.join(EX, 'pysad_kekls.ref.json')))
    ele = {}
    for ln in open(os.path.join(EX, 'elegant_kekls.touschek.txt'), encoding='utf-8'):
        if ln.startswith('#'):
            continue
        lb, rest = ln.split(':', 1)
        m = re.search(r'tLifetime \(hour\) =\s*([-0-9.e+]+)', rest)
        ele[lb.strip()] = float(m.group(1)) * 3600.0
    for a in ('0.004', '0.006', '0.01', '0.02'):
        key = [k for k in ele if k.startswith('uniform') and abs(float(k.split()[1]) - float(a)) < 6e-5]
        r = ele[key[0]] / ref['pysad_uniform'][a]
        assert 0.96 < r < 1.04, (a, r)
    k = [k for k in ele if k.startswith('uniform') and abs(float(k.split()[1]) - 0.0299) < 6e-5][0]
    assert abs(ele[k] / ref['pysad_uniform']['0.0299397749'] - 1) < 0.03
    assert abs(ele['elegant LMA'] / ref['pysad_local'] - 1) < 0.03
    assert abs(ele['elegant LMA'] / ele['pysad LMA'] - 1) < 0.01
    assert abs(ref['pysad_local'] / 26876.0 - 1) < 1e-3                      # Q1 の値 (TouschekLifetimeLMA)
