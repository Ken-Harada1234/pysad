# -*- coding: utf-8 -*-
"""v2.4.0 H (B260): K64 の qcod (周期閉軌道の探索) は収束の保証がない (大域化のない Newton 法、最大 10 回、途中で NaN なら即打ち切り)。
KEK-LS v03_68 (3224 要素) に KEKLS_correct_cod の誤差 (位置ずれ ae・強さ fe・回転 te、種 111) を 1.5 倍で与えると、k64 (= gs と同じ) は COD を見つけられず設計の光学のまま (COD ゼロ・設計チューン) で、
v2 は堅牢な解法 (信頼領域つきの Newton 法 + 線探索、不動点が不安定なら別の初期値で安定な不動点を探す) にフォールバックして COD を見つける。1 倍 (K64 が収束する) は v2 の値が従来のまま (フォールバックに入らない。k64 は gs の値 48.5762095890)。"""
import contextlib
import io
import os

from pysad import compat
from pysad.ffscmd import run_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
KL = os.path.join(ROOT, 'example', 'kek_ls')

BODY = '''
FFS;
$FORM = "15.10";
%(mod)s
ae = %(ae)g; fe = %(fe)g; te = %(te)g; GCUT = 1;
USE RING NOEXPAND;
cc = "cell calc noexpand";
FFS[cc];
GetMagnetNames[];
SeedRandom[111];
SetErrors[];
FFS[cc];
Print["TUNE ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];
cx = Twiss["DX","PM.*"]; cy = Twiss["DY","PM.*"];
Print["COD ", Max[Abs[cx]], " ", Max[Abs[cy]]];
stop;
'''


def _script(scale):
    lat = open(os.path.join(KL, 'lattice', 'kek_ls_v03_68.sad'), encoding='utf-8').read()
    mod = open(os.path.join(KL, 'modules', 'modules_errors.n'), encoding='utf-8').read()
    return lat + BODY % dict(mod=mod, ae=4e-5 * scale, fe=4e-4 * scale, te=1e-4 * scale)


def _run(scale, mode):
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            run_sad_string(_script(scale), cwd=os.path.join(KL, 'errors'))
    finally:
        compat.set_mode(old)
    t = [l.split() for l in buf.getvalue().split('\n') if l.startswith(('TUNE', 'COD'))]
    tune = [float(v) for v in t[-2][1:3]]
    cod = [float(v) for v in t[-1][1:3]]
    return tune, cod


def test_v2_robust_qcod_finds_cod_where_k64_fails():
    """1.5 倍: k64 は COD がゼロ (設計のチューン 48.5800・17.6201 のまま)、v2 は COD を見つける (最大 1.9・3.4 mm)。1 倍は従来の値のまま (K64 が収束する範囲は変わらない)"""
    tk, ck = _run(1.5, 'k64')
    assert max(ck) < 1e-10 and abs(tk[0] - 48.5800240787) < 1e-6
    tv, cv = _run(1.5, 'v2')
    # v2 の既定のチューンは、四極の本体の 4 次の合成 (B226) と、6 極 (L ≠ 0) の非線形フリンジの CALC への追加 (B291。COD が 1.9・3.4 mm 外れた 6 極に効く) で、
    # 従来の値 17.6923195779 (K64 の積分器・フリンジなし) から 5e-6 動いた。両方を戻す (方式 0・フリンジなし) と従来の値に戻る (下)。
    assert 1e-3 < cv[0] < 3e-3 and 2e-3 < cv[1] < 5e-3 and abs(tv[1] - 17.6923143372) < 1e-8
    from pysad import lattice, quadbody
    old_s, old_f = quadbody.set_scheme(0), lattice._thin_fringe_flag
    lattice._thin_fringe_flag = lambda *a: False
    try:
        tv0, cv0 = _run(1.5, 'v2')
    finally:
        quadbody.set_scheme(old_s)
        lattice._thin_fringe_flag = old_f
    assert abs(tv0[1] - 17.6923195779) < 1e-8 and abs(tv[1] - tv0[1]) > 1e-6
    t1k, c1k = _run(1.0, 'k64')
    t1v, c1v = _run(1.0, 'v2')
    assert abs(t1k[0] - 48.576209589) < 1e-8 and abs(t1v[0] - 48.5760353019) < 1e-8      # k64 = gs の値。v2 は B161 などで k64 と 2e-4 違う (qcod の経路は同じ)
    assert abs(c1v[0] - 1.3085717e-3) < 1e-9
    # v2 のチューンの従来の値 48.5760358865 (v2.4.0) からの 5.8e-7 の動きは 1.5 倍と同じ原因 (B226 の 4 次の合成と B291 の 6 極のフリンジ)。両方を戻すと従来の値と COD
    old_s, old_f = quadbody.set_scheme(0), lattice._thin_fringe_flag
    lattice._thin_fringe_flag = lambda *a: False
    try:
        t1v0, c1v0 = _run(1.0, 'v2')
    finally:
        quadbody.set_scheme(old_s)
        lattice._thin_fringe_flag = old_f
    assert abs(t1v0[0] - 48.5760358865) < 1e-8 and abs(c1v0[0] - 1.3085727e-3) < 1e-9
