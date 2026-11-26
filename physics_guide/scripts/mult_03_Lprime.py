# -*- coding: utf-8 -*-
r"""ANGLE != 0 の MULT の L' (フリンジ長さの補正) の二重適用の疑い (08_mult.md §6.4 項 4)。

K64 のソース (読み取り):
  tffs.f:1366-1371:  al = al - ((phi*fb1)^2 + (phi*fb2)^2)/al/48 * sin(w/2)/sin(phi/2)   [w = phi - psi1 - psi2]、cmp%value(p_L_MULT) = al
  tturn.f:761:       tmulti(..., cmp%value(p_L_MULT), ...)  -> tmulti.f:79 -> tmulta(al, ...)   (al はすでに L')
  tmulta.f:79-90:    aln = al/ndiv;  aln = aln - (phi*fb1)^2/al/48*sin(w/2)/sin(phi/2)/ndiv  (fb1、fb2 とも)
pysad (読み取りのみ。コードは編集しない): multa.compile_multa が tffs.f の補正を写し、multa.k_multa が tmulta.f の補正を写す。
ここでは、(1) 同じ式を数値で出して 1 回目 (L') と 2 回目 (tmulta) の補正量を並べる、
(2) compile_multa が実際に渡す al を確かめる、(3) 追跡で MULT と BEND (補正は 1 回だけ) を比べる。
"""
import io
import numpy as np
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6n_ext, multa, compat

L0, phi, fb = 1.0, 0.2, 0.2
src = ('MOMENTUM = 3GEV;\n%s MP=(L=%.15g %s);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n')


def build(kind, extra):
    ffs, _ = run_sad_string(src % (kind, L0, extra), out=io.StringIO())
    return ffs


def track(ffs, c):
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6n_ext.track_particles_ext(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


print("== (1) 補正量。L = %g m、ANGLE = %g、FB1 = FB2 = %g (FRINGE = 3)、E1 = E2 = 0 (w = phi、sin 比 = 1) ==" % (L0, phi, fb))
corr1 = ((phi * fb)**2 * 2) / L0 / 48.0
Lp = L0 - corr1
print("   1 回目 (tffs.f、L に対して):  L - L' = %.6e   -> L' = %.12f" % (corr1, Lp))
corr2 = ((phi * fb)**2 * 2) / Lp / 48.0            # スライスごとの引き算 (aln から /ndiv ずつ) の ndiv 個の合計
print("   2 回目 (tmulta.f、L' に対して。ndiv 個のスライスの合計は ndiv によらない): %.6e   -> 円弧の合計 = %.12f" % (corr2, Lp - corr2))
print("   2 回目の補正量は 1 回目とほぼ同じ (L' ≈ L のため)。")

print("\n== (2) compile_multa が tmulta に渡す al ==")
ex = 'ANGLE=%g E1=0 E2=0 FRINGE=3 FB1=%g FB2=%g DISFRIN=1' % (phi, fb, fb)
for K2 in (0.0, 40.0):
    el = [e for e in build('MULT', ex + (' K2=%g' % K2 if K2 else '')).lat.elements if e.name == 'MP'][0]
    kind_, par = multa.compile_multa(el)
    print("   K2 = %-4g al = %.12f (L' = %.12f と %s)" % (K2, par[0], Lp, '一致' if abs(par[0] - Lp) < 1e-12 else '不一致'))
print("   つまり k_multa は補正済みの L' を受け取り、さらに aln から引く (K64 の tturn → tmulti → tmulta の連鎖と同じ構造)。")

print("\n== (3) 追跡: 設計軌道からずれた粒子で、MULT (ANGLE、K2 = 0) と BEND (同じ値。補正は 1 回) を比べる。k64 (K64 の動作) と v2 (B331 を直した版) ==")
print("   (x = px = 0 の粒子は、追跡が基準軌道を引くので 0 のまま。FB1 = FB2 = 0 の対照は MULT と BEND が一致する。)")
for mode in ('k64', 'v2'):
    compat.set_mode(mode)
    print('  -- 互換モード %s' % mode)
    for lab, e in (('FB1 = FB2 = 0 (対照)', 'ANGLE=%g E1=0 E2=0 FRINGE=3 FB1=0 FB2=0 DISFRIN=1' % phi),
                   ('FB1 = FB2 = %g' % fb, ex)):
        print('   ' + lab)
        for c in ([1e-3, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0.01]):
            rm = track(build('MULT', e), c)
            rb = track(build('BEND', e), c)
            print('     初期値 %s: |MULT - BEND| の (x, px, z) = %.2e %.2e %.2e' % (c, abs(rm[0] - rb[0]), abs(rm[1] - rb[1]), abs(rm[4] - rb[4])))
compat.set_mode(None)
print("   k64 では FB ≠ 0 のときだけ差が出る (2 回目の補正による。gs でも確認: B331)。v2 は 2 回目の補正を省き、MULT が BEND と一致する。")
