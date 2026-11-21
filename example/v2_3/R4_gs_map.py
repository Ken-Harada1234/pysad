# -*- coding: utf-8 -*-
"""v2.3.0 R4: gs (SAD K64) の TrackParticles の 1 周の写像から RDT を取り出し (pysad の追跡と同じ有限差分、rdttrack.ring_map(tracker=gs))、
pysad の追跡の取り出し・解析 resonance_terms と比べる。pysad の追跡が gs と 1e-9 以内で一致する (R2・R3) ことの、RDT 用の小振幅での直接の確認。
使い方 (リポジトリ直下から)
  py -3.12 example/v2_3/R4_gs_map.py fodo|astrid|kekls
手順: ring_map の有限差分に必要な点 (約 70 点) を記録 → 台本 (NORFSW;、1 周、6 極 ×1) を WSL の gs で 1 回走らせる → 同じ計算を gs の出力で再生。
gs は timeout 1800 s、1 本ずつ、終了後に sad1.exe が残っていないか確認する。"""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import resdrive as R, rdttrack as T, compat
compat.set_mode('k64')
import R3_adts_bend as M
import Q2_adts_gs as G

CUB5 = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')


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
    amp = 1e-5
    rec = Recorder()
    T.ring_map(S['lat'].elements, S['p'], rows[0], amp=amp, orbit=np.zeros(4), tracker=rec)
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
    rm_gs = T.ring_map(S['lat'].elements, S['p'], rows[0], amp=amp, orbit=np.zeros(4), tracker=Replay(rr))
    rm_py = T.ring_map(S['lat'].elements, S['p'], rows[0], amp=amp, orbit=np.zeros(4))
    hg, hp = T.rdt_from_map(rm_gs), T.rdt_from_map(rm_py)
    ra = R.resonance_terms(S['lat'], CUB5, rows=rows)
    lines = ['R4 gs の 1 周の写像から取り出した RDT: %s (%.2f GeV、周長 %.1f m)。ring_map の振幅 %.0e √m、NORFSW、6 極 ×1、gs の追跡の座標 %d 点 (S20.15)' % (key, S['p'] / 1e9, S['circ'], amp, nps)]
    lines.append('項      gs                      pysad k64               解析 R4 (6 極 + 偏向磁石)    |gs − pysad|/|gs|   |解析 − gs|/|gs|')
    for nm in CUB5:
        lines.append('%-7s %+.6f%+.6fj  %+.6f%+.6fj  %+.6f%+.6fj  %.2e  %.2e' % (nm, hg[nm].real, hg[nm].imag, hp[nm].real, hp[nm].imag,
                     ra.total[nm].real, ra.total[nm].imag, abs(hg[nm] - hp[nm]) / abs(hg[nm]), abs(ra.total[nm] - hg[nm]) / abs(hg[nm])))
    lines.append('線形の写像 R の gs − pysad の最大差 %.2e、2 次の部分 T の最大差 %.2e (|T| 最大 %.2e)' % (np.abs(rm_gs.R - rm_py.R).max(), np.abs(rm_gs.T - rm_py.T).max(), np.abs(rm_py.T).max()))
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(HERE, 'R4_gs_map_%s.gs.txt' % key), 'w', encoding='utf-8').write(txt + '\n')
