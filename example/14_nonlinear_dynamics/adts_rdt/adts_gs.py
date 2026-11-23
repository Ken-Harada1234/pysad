# -*- coding: utf-8 -*-
"""SAD (gs) の TrackParticles で ADTS を求め、pysad の解析値 (pysad/resdrive.py) と比べる。
使い方 (リポジトリ直下から):
  py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_gs.py fodo|astrid|kekls [--lam 1,-1,0] [--turns 512]
手順: ラティスを pysad で読む → 6 極の強さを λ 倍した台本 (gs 用) を作る → WSL の gs で追跡 (リング先頭から 1 周ずつ、RF・放射なし)
  → 座標列を NAFF (pysad/naff.py) で解析 → ADTS (中心差分、振幅 b と 2b の外挿)。λ = ±1, 0 の 3 本から 6 極の 2 次の項 C = (f(1)+f(−1)−2f(0))/2
  (6 極の 1 次 (奇数次) と 6 極以外の非線形 (偏向磁石・端部) との交差項を消す) を取り出す。
gs は timeout 付き (既定 1800 s)、終了後に WSL で sad1.exe が残っていないか確認する。"""
import math, os, subprocess, sys, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from pysad.parser import load_sad_file
from pysad import naff, resdrive, compat
compat.set_mode(os.environ.get("PYSAD_COMPAT", "k64"))   # gs (K64) と比べるので既定は k64

NORFSW = '--rfsw' not in sys.argv      # 既定は NORFSW (RF なし)。--rfsw で RF ON のまま追跡 (RF の有無の対照実験用)

RINGS = {
    'fodo': (os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), 'RING', 3e9,
             {'SF': 4.0, 'SD': -6.0}, 2e-8),
    'astrid': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING', 0.58e9, None, 2e-8),
    'kekls': (os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), 'RING', 3e9, None, 2e-10),
}


def read_ring(key):
    path, line, p, k2over, jb = RINGS[key]
    src = open(path, encoding='utf-8').read()
    if k2over:        # FODO は台本の K2=0 (SF, SD) を上書きして 6 極を入れる
        src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % (k2over['SF'], k2over['SD']))
    return src, line, p, jb


def particles(tw, jb, phx=0.0, phy=0.0):
    (ax, bx), (ay, by) = tw
    pts = []
    for sc in (1.0, 2.0):
        b = sc * jb
        h = b / 2
        pts += [(b + h, b), (b - h, b), (b, b + h), (b, b - h)]
    c = []
    for jx, jy in pts:
        c.append([math.sqrt(2 * jx * bx), -math.sqrt(2 * jx / bx) * ax, math.sqrt(2 * jy * by), -math.sqrt(2 * jy / by) * ay, 0.0, 0.0])
    return pts, np.array(c)


def scale_text(src, lam, mu):
    """ラティスの文字列の K2 = 数値・K3 = 数値 を λ・μ 倍する (K2・K3 が数値だけのラティス用)。大きなリング A の定義ファイルは同じ名前の 6 極を
    何度も再定義していて (ASTRID2 も末尾で SF・SD を再定義)、gs の FFS の `SDA K2 値;` は最後の定義にしか効かないので、文字列で直す"""
    import re
    num = r'(-?[0-9]+\.?[0-9]*(?:[eE][-+]?[0-9]+)?)'
    s = re.sub(r'(\bK2\s*=\s*)' + num, lambda m: m.group(1) + repr(float(m.group(2)) * lam), src)
    return re.sub(r'(\bK3\s*=\s*)' + num, lambda m: m.group(1) + repr(float(m.group(2)) * mu), s)


def make_script(src, line, p, lam, mu_oct, x0, nturn, sext_fams, oct_fams, by_text=False):
    nps = len(x0)
    if by_text:
        src, sext_fams, oct_fams = scale_text(src, lam, mu_oct), {}, {}
    # FFS に入ると RFSW は ON (track.f:146-149)。RF つきのラティスで NORFSW を入れないと、空洞が 1 周内の z (経路長の振幅依存) に応じて δ を与え、
    # 自然クロマティシティ × δ(J) の偽の ADTS が入る (6 極 ×0 で ASTRID2 は 94.9 対 83.2、大きなリング A は 3.5 万の差)
    L = [src, 'MOMENTUM = %r;' % p, 'FFS;'] + (['NORFSW;'] if NORFSW else []) + ['$FORM="S20.15";', 'USE %s;' % line]
    for nm, v in sext_fams.items():
        L.append('%s K2 %r;' % (nm, v * lam))
    for nm, v in oct_fams.items():
        L.append('%s K3 %r;' % (nm, v * mu_oct))
    cols = ['{' + ','.join(repr(float(x0[k][j])) for k in range(nps)) + '}' for j in range(6)]
    L.append('bm = {%s,{%s}};' % (','.join(cols), ','.join(['1'] * nps)))
    L.append('Do[ r = TrackParticles[{1,bm},1,1,1]; bm = Join[Take[r[[2]],4],{0*r[[2,5]]},Drop[r[[2]],5]];')
    L.append('  Do[ Print["P ",i," ",k," ",1e3*bm[[1,k]]," ",1e3*bm[[2,k]]," ",1e3*bm[[3,k]]," ",1e3*bm[[4,k]]], {k,%d}], {i,%d}];' % (nps, nturn))
    L.append('Exit;')
    return '\n'.join(L) + '\n'


def run_gs(script, timeout=1800):
    td = tempfile.mkdtemp(prefix='q2gs_')
    f = os.path.join(td, 'q2_%s.sad' % os.path.basename(td)[5:])      # gsrun.sh は一時フォルダに basename でコピーする。同時に走る別の実行と衝突しない名前にする
    open(f, 'w', encoding='utf-8').write(script)
    sh = os.path.join(ROOT, 'local', 'wake', 'gsrun.sh')
    subprocess.run(['bash', sh, f, str(timeout)], check=False)
    out = open(f[:-4] + '.gsfull.txt', encoding='utf-8', errors='replace').read()
    return out, f


def check_gs_left():
    r = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', 'ps -eo pid,etime,pcpu,cmd | grep sad1.exe | grep -v grep'],
                       capture_output=True, text=True)
    return r.stdout.strip()


def parse(out, nps, nturn):
    rec = np.full((nturn, nps, 4), np.nan)
    for ln in out.splitlines():
        if ln.startswith('P '):
            t = ln.split()
            rec[int(t[1]) - 1, int(t[2]) - 1] = [float(v) * 1e-3 for v in t[3:7]]
    return rec


def analyse(rec, tw, pts, jb):
    (ax, bx), (ay, by) = tw
    nus = []
    for k in range(rec.shape[1]):
        r6 = np.zeros((rec.shape[0], 6))
        r6[:, :4] = rec[:, k, :]
        (fx, fy, _), _a = naff.tunes_from_coords(r6, twiss=((ax, bx), (ay, by), (0.0, 1.0)))
        nus.append((fx, fy))
    nus = np.array(nus)
    fr = lambda a: (a + 0.5) % 1.0 - 0.5
    f = []
    for s in range(2):
        n = nus[4 * s:4 * s + 4]
        b = (s + 1) * jb
        h = b / 2
        f.append(np.array([fr(n[0][0] - n[1][0]) / (2 * h), fr(n[2][0] - n[3][0]) / (2 * h), fr(n[2][1] - n[3][1]) / (2 * h)]))
    return 2 * f[0] - f[1]        # (xx, xy, yy)


if __name__ == '__main__':
    key = sys.argv[1] if len(sys.argv) > 1 else 'fodo'
    cases = [(1.0, 1.0), (-1.0, 1.0), (0.0, 1.0)]          # (λ: 6 極の倍率, μ: 8 極の倍率)
    nturn = 512
    if '--cases' in sys.argv:        # 例 --cases 1:1,-1:1,0:1,0:-1,0:0
        cases = [tuple(float(v) for v in c.split(':')) for c in sys.argv[sys.argv.index('--cases') + 1].split(',')]
    if '--turns' in sys.argv:
        nturn = int(sys.argv[sys.argv.index('--turns') + 1])
    src, line, p, jb = read_ring(key)
    if '--jb' in sys.argv:
        jb = float(sys.argv[sys.argv.index('--jb') + 1])
    from pysad.parser import load_sad_string
    sc = load_sad_string(src)
    lat = sc.build_lattice(line)
    rows = lat.twiss_table()
    r0 = rows[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    sext, octu = {}, {}
    for e in lat.elements:
        tn = type(e).__name__
        if tn == 'Sextupole':
            sext[e.name] = float(e._raw['K2'])
        elif tn == 'Octupole':
            octu[e.name] = float(e._raw['K3'])
    pts, x0 = particles(tw, jb)
    res = {}
    for lam, mu in cases:
        scr = make_script(src, line, p, lam, mu, x0, nturn, sext, octu, by_text=True)
        out, f = run_gs(scr)
        rec = parse(out, len(x0), nturn)
        if np.isnan(rec).any():
            print('gs 出力が不完全: λ=%g μ=%g' % (lam, mu), out[-500:])
            continue
        res[(lam, mu)] = analyse(rec, tw, pts, jb)
        print('gs  λ=%+g μ=%+g  ADTS (xx, xy, yy) = %s' % (lam, mu, np.array2string(res[(lam, mu)], precision=6)))
        print('  gs の残り:', check_gs_left() or 'なし')
    A = resdrive.adts(lat, rows=rows)
    mu0 = cases[0][1]
    if all((l, mu0) in res for l in (1.0, -1.0, 0.0)):
        C = (res[(1.0, mu0)] + res[(-1.0, mu0)] - 2 * res[(0.0, mu0)]) / 2
        print('gs 6 極の 2 次の項 C   =', np.array2string(C, precision=6))
        print('解析 (6 極の 2 次)     =', np.array2string(np.array([A.sext['xx'], A.sext['xy'], A.sext['yy']]), precision=6))
    if (0.0, 1.0) in res and (0.0, -1.0) in res:
        B = (res[(0.0, 1.0)] - res[(0.0, -1.0)]) / 2
        print('gs 8 極の 1 次の項 B   =', np.array2string(B, precision=6))
        print('解析 (8 極の 1 次)     =', np.array2string(np.array([A.octu['xx'], A.octu['xy'], A.octu['yy']]), precision=6))
