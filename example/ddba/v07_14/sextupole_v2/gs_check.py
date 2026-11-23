# -*- coding: utf-8 -*-
"""最終候補の gs (SAD の TrackParticles、WSL) での確認。pysad の高速追跡 (screening) と、同じ粗い格子の同じ点を gs と pysad で追跡して、生き残りの点ごとの一致と DA (光線ごとの最初に失われた振幅の手前) を比べる。
  py -3.12 gs_check.py NAME --rms 0.10 --seed 1 [--nturn 512] [--timeout 3600] [--nray 9] [--nu 24]
手順: (1) results/cands.json の K2L と、K2L のランダム誤差 (sxv2_common.sext_error_vector の同じ実現) で、ラティスのソースの K2 を置き換え、台ごとの誤差は `名前.番号 K2 値;` の文で入れる
(K2≠0 の 6 極 600 台)。(2) FFS; USE RING; CELL; CALC; NORFSW; NORAD; LOSSAMPL=0.05; (3) 光線 9 本 × 24 点 (規格化した振幅 u ≤ 4e-3、sextupole/sxt_track.da_rays と同じ格子。二分はしない) を 1 回の TrackParticles で追跡。
(4) 同じ点を pysad の追跡 (TrackSetup.track) で追跡して比べる。gs は WSL の local/wake/gsrun.sh (台本名は一意: 候補・誤差・シードを含む)、timeout は長め、終了後に sad1.exe が残っていないか確認する。
出力: results/gs_<name>_e<rms%>_s<seed>.json (既にあれば飛ばす)、.sad (台本)、.gsfull.txt (gs の出力)。"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--rms', type=float, default=0.0)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--nturn', type=int, default=512)
ap.add_argument('--timeout', type=int, default=3600)
ap.add_argument('--nray', type=int, default=9)
ap.add_argument('--nu', type=int, default=24)
ap.add_argument('--skip-pysad', action='store_true', help='同じ点の pysad 追跡を省く (gs だけ。DA は eval_*.json の二分の値と比べる)')
ap.add_argument('--cands', default=os.path.join(V.OUT, 'cands.json'))
args = ap.parse_args()
V.lower_priority()
src, lat, m, ks, xt = V.load()
cands = json.load(open(args.cands, encoding='utf-8'))
k = np.array(cands[args.name]['k'])
rms = 0.0 if args.seed == 0 else args.rms
tag = 'gs_%s_e%g_s%d' % (args.name, 100 * rms, args.seed)
path_json = os.path.join(V.OUT, tag + '.json')
if os.path.exists(path_json):
    print('skip', tag)
    sys.exit(0)
eps = None if args.seed == 0 else V.sext_error_vector(lat, rms, args.seed)

# ---- 格子 (pysad のリング先頭の β で規格化)
st = V.make_setup(src, m, k, eps)
row = st.rows[0]
bx, by = row[1], row[3]
th = np.deg2rad(np.linspace(0.0, 90.0, args.nray))
ugrid = np.linspace(4.0e-3 / args.nu, 4.0e-3, args.nu)
pts = []
for t in th:
    for u in ugrid:
        pts.append([u * math.cos(t) * math.sqrt(bx), 0.0, u * math.sin(t) * math.sqrt(by) + 1e-5, 0.0, 0.0, 0.0])
pts = np.array(pts)

# ---- pysad (screening)
t0 = time.time()
if args.skip_pysad:
    alive_py = None
    t_py = 0.0
else:
    r = st.track(1, pts, args.nturn)
    alive_py = (r['lost_turn'] == 0)
    t_py = time.time() - t0

# ---- gs の台本
s2 = V.S.replace_k2_in_source(src, m.names, k)
lines = [s2, 'FFS;', 'NORFSW; NORAD; NOFLUC;', 'USE RING; CELL; CALC;', 'LOSSAMPL = 0.05;']
# 台ごとの K2 の誤差: 成分ごとの代入 LINE["K2","名前.番号"] = 値 (CALC のあとに入れる。CALC は既定の EXPAND で成分ごとの値を消すので、この後は CALC を呼ばない。
# TrackParticles は成分ごとの値をそのまま使う: example/01_lattice_elements/FFS_USE_EXPAND.md)。K2 は K2L [1/m2] (SAD の SEXT の K2 は積分強さ)
if eps is not None:
    j = 0
    cnt = {}
    for e in lat.elements:
        if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0:
            cnt[e.name] = cnt.get(e.name, 0) + 1
            kf = float(k[m.names.index(e.name)])
            lines.append('LINE["K2","%s.%d"] = %.12g;' % (e.name, cnt[e.name], kf * (1.0 + float(eps[j]))))
            j += 1
    lines.append('Print["GSK2CHECK ", LINE["K2", {"%s.1", "%s.2", "%s.3"}]];' % ((m.names[0],) * 3))
xs = ','.join('%.10e' % v for v in pts[:, 0])
ys = ','.join('%.10e' % v for v in pts[:, 2])
n = len(pts)
lines += ['$FORM = "15.8";',
          'beam = {1, {{%s}, Table[0,{%d}], {%s}, Table[0,{%d}], Table[0,{%d}], Table[0,{%d}], Table[1,{%d}]}};' % (xs, n, ys, n, n, n, n),
          'r = TrackParticles[beam, 1, 1, %d];' % args.nturn,
          'Print["GSFLAGS ", r[[2,7]]];',
          'Print["GSX ", r[[2,1]]];',
          'Exit;']
sad_path = os.path.join(V.OUT, tag + '.sad')
open(sad_path, 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')

# ---- gs を走らせる (WSL。local/wake/gsrun.sh と同じ手順: 台本を TEMP にコピーして timeout つきで gs。
# gsrun.sh はパスに空白があると動かないので、同じことを直接書く)
tmp = os.environ.get('TEMP', os.path.join(os.path.expanduser('~'), 'AppData', 'Local', 'Temp'))
drive, rest = os.path.splitdrive(tmp)
wsl_tmp = '/mnt/%s%s' % (drive[0].lower(), rest.replace(chr(92), '/'))
name_in_tmp = '_gs_%s.sad' % tag
import shutil
shutil.copy(sad_path, os.path.join(tmp, name_in_tmp))
t0 = time.time()
print('gs 実行:', tag, flush=True)
gs_out = sad_path[:-4] + '.gsfull.txt'
with open(gs_out, 'wb') as fo:
    subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', 'cd %s && timeout %d /mnt/d/work/SAD-master/SAD-master/bin/gs %s < /dev/null' % (wsl_tmp, args.timeout, name_in_tmp)],
                   stdout=fo, stderr=subprocess.STDOUT, check=False)
t_gs = time.time() - t0
left = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', 'pgrep -a sad1 || true'], capture_output=True, text=True).stdout.strip()
print('gs 終了後の sad1 プロセス:', left if left else 'なし', flush=True)
txt = open(sad_path[:-4] + '.gsfull.txt', encoding='utf-8', errors='replace').read()
flags = None
i0 = txt.find('GSFLAGS {')
if i0 >= 0:
    i1 = txt.find('}', i0)
    body = txt[i0 + len('GSFLAGS'):i1].replace('{', ' ').replace(',', ' ')
    flags = np.array([float(v) for v in body.split()])
if flags is None or len(flags) != n:
    print('gs の出力が読めない (timeout か台本のエラー)。出力の末尾:')
    print(txt[-1500:])
    sys.exit(1)
alive_gs = flags > 0.5


def da_from(alive):
    a = alive.reshape(len(th), len(ugrid))
    u_da = []
    for i in range(len(th)):
        bad = np.where(~a[i])[0]
        j = bad[0] if len(bad) else len(ugrid)
        u_da.append(ugrid[j - 1] if j > 0 else 0.0)
    u_da = np.array(u_da)
    xs_ = u_da * np.cos(th) * math.sqrt(bx)
    ys_ = u_da * np.sin(th) * math.sqrt(by)
    px = np.concatenate([[0.0], xs_, [0.0]])
    py = np.concatenate([[0.0], ys_, [0.0]])
    area = 0.5 * abs(np.sum(px[:-1] * py[1:] - px[1:] * py[:-1]))
    return dict(u=u_da.tolist(), area_mm2=float(4 * area * 1e6), x_mm=float(xs_[0] * 1e3), y_mm=float(ys_[-1] * 1e3))


out = dict(name=args.name, rms=rms, seed=args.seed, nturn=args.nturn, nray=args.nray, nu=args.nu, t_gs=t_gs, t_pysad=t_py, n=int(n), gs=da_from(alive_gs),
           alive_gs=alive_gs.astype(int).tolist())
if alive_py is not None:
    out.update(agree=float(np.mean(alive_gs == alive_py)), n_disagree=int(np.sum(alive_gs != alive_py)), pysad=da_from(alive_py), alive_pysad=alive_py.astype(int).tolist())
json.dump(out, open(path_json, 'w', encoding='utf-8'), ensure_ascii=False)
if alive_py is not None:
    print('%s: gs %.0f s、pysad %.0f s | 生き残りの一致 %.4f (%d/%d 点が違う) | DA 面積 gs %.2f、pysad %.2f mm2 (9 光線 x %d 点、二分なし)' % (
        tag, t_gs, t_py, out['agree'], out['n_disagree'], n, out['gs']['area_mm2'], out['pysad']['area_mm2'], args.nu))
else:
    print('%s: gs %.0f s | DA 面積 gs %.2f mm2 (x %.2f mm、y %.2f mm。9 光線 x %d 点、二分なし)' % (tag, t_gs, out['gs']['area_mm2'], out['gs']['x_mm'], out['gs']['y_mm'], args.nu))
