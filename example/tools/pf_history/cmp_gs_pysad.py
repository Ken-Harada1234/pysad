# -*- coding: utf-8 -*-
"""4 つの例の gs の画面表示 (.gs.txt) と pysad (k64 モード) の画面表示を、計算結果の部分だけ比べる。
使い方: py -3.12 example/tools/pf_history/cmp_gs_pysad.py   (差があれば行を出す)"""
import os, sys, io, re
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE)); ROOT = os.path.dirname(EX)
sys.path.insert(0, ROOT)
CASES = ['1981_460nmrad', '1986_128nmrad', '1997_036nmrad', '2005_036nmrad']
KEYS = ('circumference', 'number of elements', 'tune nx ny', 'EMITX', 'MOMENTUM [eV]', 'bends =', 'CAVI VOLT', 'with the design', 'natural')
def results(txt):
    """Print の結果の行 (入力の復唱ではない行) と、emit の表の主な行を返す"""
    out = []
    for ln in txt.splitlines():
        s = ln.rstrip()
        if s.startswith(KEYS) or re.match(r'^\w+\s+# ', s) or re.match(r'^(Energy loss|Equilibrium|Bucket|Design momentum|Orbit dilation|Emittance|Bunch Length|Beam size)', s):
            out.append(s)
        elif re.match(r'^\s+X :\s', s) and ('Damp' in ''.join(out[-1:]) or True):
            pass
    # emit の表の減衰行
    for m in re.finditer(r'Damping time \(sec\):\n\s*X :\s*(.*)\n', txt):
        out.append('damping time ' + m.group(1))
    for m in re.finditer(r'Damping partition number:\n\s*X :\s*(.*)\n', txt):
        out.append('partition ' + m.group(1))
    return out
NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?')
def same(a, b, rtol=1e-8):
    """数値を取り出して相対 1e-8 で比べる。εy・ビームの傾き (0 の雑音) は除く。alpha の表示桁は gs が 5 桁、pysad が 3 桁なので 1e-3"""
    for pat in (r'Emittance Y\s+=\s*\S+ m', r'Beam tilt\s+=\s*\S+ rad'):
        a = re.sub(pat, '', a); b = re.sub(pat, '', b)
    ta = re.sub(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?', '#', a); tb = re.sub(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?', '#', b)
    if re.sub(r'\s+', ' ', ta) != re.sub(r'\s+', ' ', tb): return False
    fa = [float(x.replace('D', 'E')) for x in NUM.findall(a)]; fb = [float(x.replace('D', 'E')) for x in NUM.findall(b)]
    tol = 1e-3 if 'Momentum compact' in a else rtol
    return len(fa) == len(fb) and all(abs(x - y) <= tol * max(abs(x), abs(y), 1e-30) + 1e-14 for x, y in zip(fa, fb))
def run_pysad(case, mode):
    os.environ['PYSAD_COMPAT'] = mode
    from pysad import compat, ffscmd
    compat.set_mode(mode)
    src = open(os.path.join(EX, 'pf', 'optics', 'FFS_cell_calc_emit_PF_%s.sad' % case), encoding='utf-8').read()
    buf = io.StringIO(); ffscmd.run_sad_string(src, out=buf)
    return buf.getvalue()
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'k64'
    nd = 0
    for c in CASES:
        g = results(open(os.path.join(EX, 'pf', 'optics', 'FFS_cell_calc_emit_PF_%s.gs.txt' % c), encoding='utf-8', errors='replace').read())
        p = results(run_pysad(c, mode))
        d = [(a, b) for a, b in zip(g, p) if not same(a, b)]
        print(c, 'gs lines', len(g), 'pysad lines', len(p), 'differences', len(d) + abs(len(g) - len(p)))
        for a, b in d[:6]:
            print('  gs   ', a[:200]); print('  pysad', b[:200])
        nd += len(d) + abs(len(g) - len(p))
    print('total differences', nd)
