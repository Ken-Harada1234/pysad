#!/usr/bin/env python3
"""
bench_level27a_pysad.py — Level 27a (MAIN レベルの TRACK コマンド) の pysad 側
参照: bench_level27a.result (SAD K64 で bench_level27a.sad を実行した値)

SAD K64 では MAIN の EMIT/TRACK が動かない (EMIT は何もせず、TRACK は emit の正規モード行列 r が未設定で SIGSEGV; docs/work_main_track.md)。
そこで bench_level27a.sad は ActTra.f / tracka.f の手順 (規格化振幅 → r·amp + codin、周の頭の規格化座標から 2J の和と tsmear、
FFT のピーク wtunem) を SADScript で書き、追跡だけを TrackParticles (tturn) に任せた版。
pysad 側は bench_level27a_main.sad (古典 SAD の MAIN の書き方: ON/OFF, EMIT USE=RING; TRACK USE=RING NPART= TURNS= NX= NY= NZ=;)
を run_sad.py と同じ経路で実行し、maintrack.run_track の全桁データ (smear・最終座標・TunesByTracking・EMITX 等) を SAD と比べる。
ケース: N1 = RFSW ON / RAD OFF、R1 = RAD ON / FLUC OFF、F1 = RAD ON / FLUC ON (SeedRandom[17])。
"""
import io
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad.ffscmd import run_sad_string
from pysad import track6 as _t6

NUM = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?'


def read_sad(name):
    s = open(os.path.join(HERE, name + '.result'), encoding='utf-8', errors='replace').read()
    s = re.sub(r',\s*\n\s*', ', ', s)
    return s


def nums(txt):
    return [float(v.replace('D', 'E').replace('d', 'e')) for v in re.findall(NUM, txt)]


def sad_case(s, tag, np0=6):
    d = {}
    m = re.search(r'^' + tag + r' EM \{(.*)\}', s, re.M)
    d['em'] = nums(m.group(1))
    m = re.search(r'^' + tag + r' COD \{(.*)\}', s, re.M)
    d['cod'] = nums(m.group(1))
    d['smear'] = {}
    for m in re.finditer(r'^' + tag + r' SMEAR\s+(' + NUM + r')\s+(' + NUM + r') \{(.*)\}', s, re.M):
        d['smear'][(int(float(m.group(1))), int(float(m.group(2))))] = nums(m.group(3))
    d['fin'] = {}
    for m in re.finditer(r'^' + tag + r' FIN P\s*(' + NUM + r') \{(.*)\}', s, re.M):
        d['fin'][int(float(m.group(1)))] = nums(m.group(2))
    d['lost'] = {}
    for m in re.finditer(r'^' + tag + r' LOST\s+(' + NUM + r')\s+(' + NUM + r')', s, re.M):
        d['lost'][int(float(m.group(1)))] = int(float(m.group(2)))
    d['tune'] = {}
    d['ampl'] = {}
    for key, nm in (('tune', 'TUNE'), ('ampl', 'AMPL')):
        for m in re.finditer(r'^' + tag + ' ' + nm + r' P\s*(' + NUM + r')\s+(' + NUM + r') \{(.*)\}', s, re.M):
            d[key][(int(float(m.group(1))), int(float(m.group(2))))] = nums(m.group(3))
    return d


def main(case='a'):
    name = 'bench_level27' + case
    tags = ('N1', 'R1', 'F1') if case == 'a' else ('N1', 'F1')
    sad = read_sad(name)
    src = open(os.path.join(HERE, name + '_main.sad'), encoding='utf-8').read()
    buf = io.StringIO()
    _t6.RNG.restore(17)
    ffs, I = run_sad_string(src, out=buf, cwd=HERE)
    text = buf.getvalue()
    open(os.path.join(HERE, name + '_pysad.out'), 'w', encoding='utf-8').write(text)
    print(f'pysad: TRACK {len(ffs.track_results)} 回, 出力 {len(text.splitlines())} 行 ({name}_pysad.out)')
    ok = True
    for tag, res in zip(tags, ffs.track_results):
        d = sad_case(sad, tag)
        print(f'--- {tag}')
        e_sad = d['em']
        e_py = [res['emx'], res['emz']]
        rel = [abs(a - b) / max(abs(b), 1e-300) for a, b in zip(e_py, e_sad[4:6])]
        print(f'  emx, emz   pysad {e_py[0]:.10e} {e_py[1]:.10e}  SAD {e_sad[4]:.10e} {e_sad[5]:.10e}  相対差 {max(rel):.2e}')
        dcod = max(abs(a - b) for a, b in zip(res['codin'], d['cod']))
        print(f'  codin 最大差 {dcod:.2e}')
        mxr = 0.0
        cnt = 0
        got = set()
        for (n, i, v) in res['smear']:
            ref = d['smear'].get((n, i))
            got.add((n, i))
            if ref is None:
                print(f'  SMEAR ({n},{i}) が SAD に無い')
                ok = False
                continue
            for a, b in zip(v, ref):
                if abs(b) > 1e-30:
                    mxr = max(mxr, abs(a - b) / abs(b))
            cnt += 1
        missing = [k for k in d['smear'] if k not in got]
        print(f'  smear {cnt} 行: 最大相対差 {mxr:.2e} (SAD にだけある {len(missing)} 行)')
        # 生きている粒子は 6 座標、落ちた粒子は x, y, z, dp だけ比べる (SAD の TrackParticles は落ちた粒子の px, py を
        # tconvm の対象 (生きている npa 個) から外すので px/py の値が違う。tapert の損失報告 (x, y, z) には px, py は出ない)
        dfin = 0.0
        for i in d['fin']:
            cols = [0, 1, 2, 3, 4, 5] if d['fin'][i][6] > 0 else [0, 2, 4, 5]
            dfin = max(dfin, float(np.max(np.abs(res['final'][i - 1][cols] - np.array(d['fin'][i][:6])[cols]))))
        print(f'  最終座標 (落ちた粒子は落ちた時点) 最大差 {dfin:.2e}')
        lost_py = {pid: t for (pid, t, l) in res['losses']}
        print(f'  損失: pysad {sorted(lost_py.items())}  SAD {sorted(d["lost"].items())}')
        dt = 0.0
        dar = 0.0
        for i, rec in enumerate(res['tunes'], start=1):
            for k in range(3):
                f_sad = np.array(d['tune'][(i, k + 1)])
                a_sad = np.array(d['ampl'][(i, k + 1)])
                f_py = np.array([p[0] for p in rec[1 + k]])
                a_py = np.array([p[1] for p in rec[1 + k]])
                dt = max(dt, float(np.max(np.abs(f_py - f_sad))))
                big = a_sad > 1e-30
                if np.any(big):
                    dar = max(dar, float(np.max(np.abs(a_py[big] - a_sad[big]) / a_sad[big])))
        print(f'  tune (8 本 × 3 面 × {len(res["tunes"])} 粒子) 最大差 {dt:.2e}, 振幅 最大相対差 {dar:.2e}')
        good = (max(rel) < 1e-8 and dcod < 1e-9 and mxr < 1e-6 and dfin < 1e-9 and dt < 1e-6 and dar < 1e-5
                and lost_py == d['lost'] and not missing)
        print('  ' + ('一致' if good else '不一致'))
        ok = ok and good
    print('Matched' if ok else 'MISMATCH')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'a')
