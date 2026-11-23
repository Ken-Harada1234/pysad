# -*- coding: utf-8 -*-
"""ヘルプ (SADHelp.html) の BEND の式の確認 (05_bend.md 付録 A)。
 (a) 線形フリンジの V: 平文版と TeX 版で式が違う。どちらがヘルプ自身の x2・py2・z2 と合うか (exp(:V:) は V を時間 1 のハミルトニアンとする流れ)。
 (b) 本体の s: ヘルプの s = ANGLE ρb (asin - asin + ω) と、経路長 s = ρb p_h Θ の比較 (ANGLE× の有無)。
 (c) 出口の「面からの直進」の ψ2 の定義 (ヘルプは ANGLE×E2 で AE2 が無い)。 → 原文の読み (script では扱わない)
 (d) 非線形フリンジの py2 の 1/p1 (導出との比。数値は bend_03_fringe.out.txt (b))。
"""
import math
import sympy as sp
import _common as C

C.start(__file__)
x, y, px, py, p = sp.symbols('x y p_x p_y p', real=True)
F1, rb = sp.symbols('F1 rho_b', positive=True)
dxfr = F1**2 / (24 * rb)
dyfr = F1 / (6 * rb**2)
dyfra = sp.Rational(2, 3) / (F1 * rb**2)

V_ascii = -F1**2 / rb * px / p / 24 - F1 / rb**2 * y**2 / p / 12 + 1 / rb**2 / F1 * y**4 / p / 6
V_tex = -(F1**2 / (24 * rb) * px - F1 / (12 * rb) * y**2 + y**4 / (6 * F1 * rb**2)) / p
help_py = (dyfr - dyfra * y**2) * y / p       # q の単位 (py/p) では 1/p^2。正準なら 1/p
help_z = dxfr * px / p + (dyfr / 2 * y**2 - dyfra * y**4 / 4) / p**2   # q の単位の px を p px' と読んだ形

print('== (a) 線形フリンジ V (入口) ==')
for name, V in (('平文版', V_ascii), ('TeX 版', V_tex)):
    dpy = sp.simplify(-sp.diff(V, y))
    dx = sp.simplify(sp.diff(V, px))
    print(' %s: Δx = ∂V/∂px = %s、Δpy = -∂V/∂y = %s' % (name, dx, dpy))
    print('   ヘルプの py2 (正準の 1/p で読む) との差: %s' % sp.simplify(dpy - help_py))
    print('   次元 (F1・ρb・y [m]、px 無次元): 各項 → ' + ', '.join(
        '%s' % sp.simplify(term.subs({F1: sp.Symbol('m'), rb: sp.Symbol('m'), y: sp.Symbol('m')})) for term in sp.Add.make_args(sp.expand(V * p))))
print(' 平文版の ∂V/∂p (z の更新) と、ヘルプの z2 の括弧 [dxfr px + (dyfr - dyfra y²/2) y²/(2p)]/p の比較 (px を p·q と読む):')
dz = sp.simplify(sp.diff(V_ascii.subs(px, px), p))
print('   ∂V_ascii/∂p = %s' % dz)
print('   差 (px → q p, ヘルプ z2 の px1 を q と読む) = %s' % sp.simplify(dz.subs(px, 0) - help_z.subs(px, 0)))
print(' → 平文版 V はヘルプ自身の変換式と一致 (x と py)。TeX 版は y² の項の ρb のべき・符号、y⁴ の項の符号が違い、y² の項は次元も合わない (誤植)。')

print('\n== (b) 本体の s ==')
print(' 経路長 (面内) s = ρb p_h Θ (05 章 §5.3。括弧内 asin - asin + ω は円弧の角 Θ)。ヘルプの s は ANGLE × ρb × Θ で、ANGLE× が余計で p_h が無い (p_h = 1 の近軸では ρb Θ)。')
print(' 検証 (数値): bend_02_arc_faces.out.txt (d) の「ヘルプの本体の式」の欄。')

print('\n== (d) 非線形フリンジの py2 ==')
q = sp.symbols('q_x', real=True)
a_ = 1 - y**2 / (6 * rb**2)
py_help = -px * a_ * y / (p * rb * sp.sqrt(p**2 - px**2))
py_derived = -px * a_ * y / (rb * sp.sqrt(p**2 - px**2))       # 05 章 (5.33)
print(' ヘルプ / 導出 = %s' % sp.simplify(py_help / py_derived))
print(' → ヘルプは 1/p1 を余分に持つ (p1 = 1 では一致)。')

print('\n== (e) 本体の px2: TeX 版と平文版の違い ==')
print(' 平文版: ... + sin(ω) pz1 + cos(ω) px1 - x1/ρb sin(ω) (= -(x1/ρb) sin ω)。TeX 版: - x1/(ρb sin ω) (sin ω が分母)。')
print(' TeX 版は ω → 0 (ANGLE = ψ1 + ψ2) で発散し、平文版は x1 の係数が -sin ω/ρb → 0 (面から面への直進に戻る)。平文版を 05 章 (5.22) と同じ形として使う。')
