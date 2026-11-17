# -*- coding: utf-8 -*-
"""v2 の修正 B59: CSRSetupOY の zj (位置) と dpj (運動量)。K64 は dpjn を zjn と同じ記憶域に割り当て、zj に dp の値が入り dpj は 0。

根拠: csrjasign は作用 J_j の軌道を、折り返し点 (p=0) から半周期 (π/ω_j) を mphi 刻みの蛙跳びで追う。
  - z_{k+1} - z_k = p_k r56 dt1 (dt1 = π/ω_j/mphi)。
  - 半周期の間 p の符号は変わらず、両端 (折り返し点) で |p| が小さい。
  - z はバンチの範囲 (nr·dz 以内) の値。
"""
import numpy as np

from pysad.tests.test_csr import _run, _nums

BODY = '''nr=256; dz=2e-5;
zl=Table[0.02*Cos[0.3 i]/(1+i/10.), {i,1,nr}];
rho=CSRHaissin[zl,{dz,1e-8,-0.3,3e-4}];
s=CSRSetupOY[rho,zl,{dz,6,3,1e-8,0.05,-0.3,3}];
Print["OM ",s[[3]]];
Do[Print["ZJ ",Take[s[[4]],{48 (j-1)+1,48 j}]];Print["PJ ",Take[s[[5]],{48 (j-1)+1,48 j}]],{j,6}];
'''


def _get():
    lines = _run(BODY).splitlines()
    g = lambda t: np.array([_nums(l) for l in lines if l.startswith(t + ' ')])
    return g('OM')[0], g('ZJ'), g('PJ')


def test_b59_zj_and_dpj_are_separate(v2_mode):
    om, zj, pj = _get()
    r56 = 0.05
    assert np.all(np.abs(zj) < 256 * 2e-5)
    assert np.abs(pj).max() > 0.0
    for j in range(6):
        dt1 = np.pi / om[j] / 48
        assert np.allclose(np.diff(zj[j]), pj[j, :-1] * r56 * dt1, rtol=1e-9, atol=1e-18)
        assert np.all(pj[j] > 0) or np.all(pj[j] < 0)                 # 半周期の間 p の符号は同じ
        pm = np.abs(pj[j]).max()
        assert abs(pj[j, 0]) < 0.2 * pm and abs(pj[j, -1]) < 0.2 * pm    # 両端は折り返し点


def test_b59_k64_mode_overwrites_zj(k64_mode):
    om, zj, pj = _get()
    assert np.all(pj == 0.0)
