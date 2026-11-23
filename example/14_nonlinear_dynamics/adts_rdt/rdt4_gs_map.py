# -*- coding: utf-8 -*-
"""gs (SAD K64) の TrackParticles の 1 周の写像から 4 次の生成関数 f4 を取り出し (rdttrack.ring_map3(tracker=gs))、pysad の追跡の取り出し・解析 rdt4 と比べる。
使い方 (リポジトリ直下から)
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_gs_map.py fodo|astrid
手順: ring_map3 の有限差分に必要な点を記録 → 台本 (NORFSW;、1 周、6 極 ×1) を WSL の gs で 1 回走らせる → 同じ計算を gs の出力で再生。
gs は timeout 1800 s、1 本ずつ、終了後に sad1.exe が残っていないか確認する (残っていれば PID で kill)。"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import rdttrack as T, rdt2nd as Q, compat
compat.set_mode('k64')
import adts_bend as M
import adts_gs as G

SHOW = ['h22000', 'h11110', 'h00220', 'h31000', 'h40000', 'h20020', 'h20200', 'h11200', 'h20110', 'h00310', 'h00400']


class Recorder:
    def __init__(self):
        self.pts = []

    def __call__(self, Z):
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
        self.pts.append(Z.copy())
        return np.zeros_like(Z)


class Replay:
    def __init__(self, out):
        self.out, self.i = out, 0

    def __call__(self, Z):
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
        r = self.out[self.i:self.i + len(Z)]
        self.i += len(Z)
        return r


if __name__ == '__main__':
    key = sys.argv[1] if len(sys.argv) > 1 else 'fodo'
    S = M.setup(key, 1.0)
    rows = S['rows']
    amp, amp3 = 1e-5, 1e-4
    rec = Recorder()
    T.ring_map3(S['lat'].elements, S['p'], rows[0], amp=amp, amp3=amp3, orbit=np.zeros(4), tracker=rec)
    Zall = np.vstack(rec.pts)
    nps = len(Zall)
    x0 = np.zeros((nps, 6))
    x0[:, :4] = Zall
    script = G.make_script(S['src'], S['line'], S['p'], 1.0, 1.0, x0, 1, {}, {}, by_text=False)
    print('%s: %d 点を gs で 1 周追跡 (NORFSW、6 極 ×1)' % (key, nps), flush=True)
    out, f = G.run_gs(script, 1800)
    left = G.check_gs_left()
    print('gs 終了後の sad1.exe: %s' % (left or 'なし'))
    rr = G.parse(out, nps, 1)[0]
    if np.isnan(rr).any():
        print('gs の出力が不完全 (NaN %d 個)' % int(np.isnan(rr).sum()))
        print(out[-2000:])
        sys.exit(1)
    rm_gs = T.ring_map3(S['lat'].elements, S['p'], rows[0], amp=amp, amp3=amp3, orbit=np.zeros(4), tracker=Replay(rr))
    rm_py = T.ring_map3(S['lat'].elements, S['p'], rows[0], amp=amp, amp3=amp3, orbit=np.zeros(4))
    hg, ig = T.rdt4_from_map(rm_gs)
    hp, ip = T.rdt4_from_map(rm_py)
    ra = Q.rdt4(S['lat'], rows=rows)
    lines = ['R5 gs の 1 周の写像の 3 次のテイラー展開から取り出した 4 次の h4: %s (%.2f GeV、周長 %.1f m)。amp=%.0e、amp3=%.0e √m、NORFSW、6 極 ×1、gs の追跡の座標 %d 点 (S20.15)' % (key, S['p'] / 1e9, S['circ'], amp, amp3, nps)]
    lines.append('項      gs                      pysad k64               解析 rdt4                    |gs − pysad|/|gs|   |解析 − gs|/|gs|')
    for nm in SHOW:
        lines.append('%-7s %+.5f%+.5fj  %+.5f%+.5fj  %+.5f%+.5fj  %.2e  %.2e' % (nm, hg[nm].real, hg[nm].imag, hp[nm].real, hp[nm].imag,
                     ra.total[nm].real, ra.total[nm].imag, abs(hg[nm] - hp[nm]) / abs(hg[nm]), abs(ra.total[nm] - hg[nm]) / abs(hg[nm])))
    lines.append('線形の写像 R の gs − pysad の最大差 %.2e、2 次の部分 T の最大差 %.2e (|T| 最大 %.2e)' % (np.abs(rm_gs.R - rm_py.R).max(), np.abs(rm_gs.T - rm_py.T).max(), np.abs(rm_py.T).max()))
    lines.append('f4 のシンプレクティック残差: gs %.1e、pysad %.1e' % (ig['symp'], ip['symp']))
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(HERE, 'rdt4_gs_map_%s.gs.txt' % key), 'w', encoding='utf-8').write(txt + '\n')
