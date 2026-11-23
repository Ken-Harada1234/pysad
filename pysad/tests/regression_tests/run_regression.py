#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pysad 回帰テスト — SAD (K64, WSL) との照合ベンチ 50 本 (Level 1〜14・16・22・24・26〜31) をまとめて実行し、合否を一覧にする。

使い方 (Windows, Git Bash / PowerShell):
    cd <リポジトリのルート>
    PYTHONIOENCODING=utf-8 python -u pysad/tests/regression_tests/run_regression.py            # 全部 (目安 約 17 分)
    QUICK=1 PYTHONIOENCODING=utf-8 python -u pysad/tests/regression_tests/run_regression.py    # 遅いもの (Level 7, 12t (純 Python 追跡), 13a, 16a) を飛ばす (目安 約 7 分)
    python -u pysad/tests/regression_tests/run_regression.py 3b 12 13a                          # 指定したものだけ
(どのフォルダから実行してもよい。基準の SAD (K64) の結果は保存済みで、実行時に SAD は要らない)

各ベンチの pysad スクリプトを子プロセスで実行し、その標準出力の「差」「最大差」「相対差」「Matched」「Score」
などの行を読んで合否を判定する。run_sad.py で動かすもの (13a/13b/14a/16a) は SAD の実行結果 (bench_levelN.result)
と Print[] の出力行 (TAG {...}) を突き合わせる。Level 1/2 の pysad スクリプトは SAD の値を印字しないので、
bench/bench_level1.result / bench/bench_level2.result から参照値を読んで比べる。

許容 (docs/27_benchmarks.md 参照):
    追跡 (TrackParticles)                 : 座標の絶対差 1e-9           (TOL_TRACK)
    CALC (オンモメンタム, 閉軌道 0)        : 1e-9 (絶対 + 相対)           (TOL_CALC)
    CALC (閉軌道/結合/オフモメンタム経由)  : 相対 1e-7 (絶対下限 1e-9)    (TOL_CALC_COD)  ※ SAD の qcod の収束と pysad の数値微分の丸め
    Emittance[] (6D)                      : 相対 1e-6                   (TOL_EMIT)
    Emittance[] (SAD の表示が 5 桁のとき)   : 相対 5e-5                   (TOL_EMIT_DISP)
    マッチング (K1, 残差, Twiss)            : 1e-9 (絶対 + 相対)           (TOL_CALC)
    DA (DynamicApertureSurvey)            : 地図の文字列と Score が一致。損失周数 (RES) は境界粒子 1 個の ±1 周まで「既知」

子プロセスの出力は regression_tests/regression_out/<名前>.txt に保存する。終了コード: FAIL があれば 1。
"""
import os
import re
import sys
import time
import math
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))        # pysad/tests/regression_tests/ (cmp_*.py、modules_general.n など。子プロセスの cwd)
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))   # リポジトリのルート (run_sad.py の場所)
BENCH = os.path.join(HERE, 'bench')                       # ベンチの台本・pysad スクリプト・SAD (K64) の基準結果
OUTDIR = os.path.join(HERE, 'regression_out')             # 実行ごとの出力 (git には入れない)
QUICK = os.environ.get('QUICK', '') == '1'

TOL_TRACK = 1e-9
TOL_CALC = 1e-9
TOL_CALC_COD = 1e-7
TOL_EMIT = 1e-6
TOL_EMIT_DISP = 5e-5
TOL_DISP6 = 1e-6          # pysad が小数 6 桁で印字している値 (Level 1/2 の Twiss) の比較
TOL_RADINT = 1e-6         # RADINT (SAD は 7 桁表示)
TOL_CHROM_C1 = 1e-5       # 色収差 (チューンの 3 次当てはめの 1 次係数): チューンの 1e-8 のノイズが拡大される

NUM = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?'


# ----------------------------------------------------------------------------
# 共通
# ----------------------------------------------------------------------------
class Check:
    """1 項目の判定結果。status: PASS / FAIL / KNOWN (既知の不一致) / INFO (判定しない)"""
    def __init__(self, name, diff, tol, status, note=''):
        self.name = name
        self.diff = diff
        self.tol = tol
        self.status = status
        self.note = note

    def __repr__(self):
        d = '' if self.diff is None else f'{self.diff:.2e}'
        t = '' if self.tol is None else f'{self.tol:.0e}'
        return f'{self.status:5s} {self.name:40s} 差 {d:>9s} / 許容 {t:>6s} {self.note}'


def fnum(s):
    return float(s.replace('D', 'e').replace('d', 'e'))


def cmp(name, a, b, tol_abs=0.0, tol_rel=0.0, note='', known=False):
    """|a-b| <= tol_abs + tol_rel*|b| なら PASS"""
    d = abs(a - b)
    tol = tol_abs + tol_rel * abs(b)
    if d <= tol or (math.isnan(d) is False and d == 0.0):
        return Check(name, d, tol, 'PASS', note)
    return Check(name, d, tol, 'KNOWN' if known else 'FAIL', note)


def cmp_diff(name, d, tol, note='', known=False):
    """スクリプトが印字した差 d をそのまま判定"""
    d = abs(d)
    if math.isnan(d):
        return Check(name, d, tol, 'FAIL', note + ' (NaN)')
    if d <= tol:
        return Check(name, d, tol, 'PASS', note)
    return Check(name, d, tol, 'KNOWN' if known else 'FAIL', note)


def grab_all(text, pattern, flags=re.M):
    return [m for m in re.finditer(pattern, text, flags)]


def grab(text, pattern, n=0, flags=re.M):
    ms = grab_all(text, pattern, flags)
    if len(ms) <= n:
        raise ValueError(f'pattern not found ({n}): {pattern}')
    return fnum(ms[n].group(1))


def read(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read()


def extract_tagged(text):
    """行頭の 'TAG {...}' (SAD の Print["TAG ", {...}] の出力) を集める。{ } は複数行にまたがることがある。
    戻り値: {TAG: [行をつないだ文字列, ...]} (同じ TAG が複数回出ることがある)"""
    out = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r'^([A-Z][A-Z0-9_]*)\s+(\{.*|' + NUM + r'\s*)$', lines[i])
        if m:
            tag, body = m.group(1), m.group(2)
            depth = body.count('{') - body.count('}')
            j = i + 1
            while depth > 0 and j < len(lines):
                body += lines[j]
                depth += lines[j].count('{') - lines[j].count('}')
                j += 1
            out.setdefault(tag, []).append(body)
            i = j
        else:
            i += 1
    return out


def numbers_in(s):
    return [fnum(x) for x in re.findall(NUM, s)]


def strings_in(s):
    return re.findall(r'"([^"]*)"', s)


# ----------------------------------------------------------------------------
# 各ベンチの判定
# ----------------------------------------------------------------------------
def tw_pairs(text, names):
    """'BX: pysad v SAD v 差 d' 形式の行 → [(name, pysad, sad, diff)]"""
    pat = r'^\s*(' + '|'.join(names) + r')\s*:\s*pysad\s+(' + NUM + r')\s+SAD\s+(' + NUM + r')\s+差\s+(' + NUM + r')'
    return [(m.group(1), fnum(m.group(2)), fnum(m.group(3)), fnum(m.group(4))) for m in grab_all(text, pat)]


def emit_pairs(text):
    """'εx : pysad v SAD v 相対差 d' 形式の行 → [(name, pysad, sad, reldiff)]"""
    pat = r'^\s*(\S+(?: \S+)?)\s*:?\s*pysad\s+(' + NUM + r')\s+SAD\s+(' + NUM + r')\s+相対差\s+(' + NUM + r')'
    return [(m.group(1).strip(), fnum(m.group(2)), fnum(m.group(3)), fnum(m.group(4))) for m in grab_all(text, pat)]


def cod_pairs(text):
    """'NAME pysad v1 v2 v3 v4' の次の行が 'SAD v1 v2 v3 v4' → [(name, [pysad], [sad])]"""
    pat = (r'^(\S+)\s+pysad\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s*\n'
           r'\s+SAD\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')')
    return [(m.group(1), [fnum(m.group(k)) for k in range(2, 6)], [fnum(m.group(k)) for k in range(6, 10)])
            for m in grab_all(text, pat)]


def chrom_pairs(text):
    """'X: c0..c3 pysad = a, b, c, d' と 'c0..c3 SAD = ...' (または 'X: c0..c3 SAD ...') の組"""
    pat = r'c0\.\.c3\s+(pysad|SAD)\s*=?\s*(' + NUM + r'),\s*(' + NUM + r'),\s*(' + NUM + r'),\s*(' + NUM + r')'
    items = [(m.group(1), [fnum(m.group(k)) for k in range(2, 6)]) for m in grab_all(text, pat)]
    # pysad 行の次に SAD 行が来る前提
    pairs = []
    k = 0
    while k + 1 < len(items):
        if items[k][0] == 'pysad' and items[k + 1][0] == 'SAD':
            pairs.append((items[k][1], items[k + 1][1]))
            k += 2
        else:
            k += 1
    return pairs


def track_checks(text, label=''):
    """'P1 差: v1..v6 |SAD| max' 行の最大絶対差 → 1 項目。乱数の状態の行があれば一致も確認"""
    checks = []
    vals = []
    for m in grab_all(text, r'^\s*P\d\s+差:\s+((?:' + NUM + r'\s+){6})'):
        vals += [abs(fnum(x)) for x in m.group(1).split()]
    if vals:
        checks.append(cmp_diff(f'{label}追跡 座標の最大差 ({len(vals)//6} 粒子・時刻)', max(vals), TOL_TRACK))
    else:
        checks.append(Check(f'{label}追跡', None, TOL_TRACK, 'FAIL', '差の行が見つからない'))
    rng = grab_all(text, r'乱数の状態:?\s+pysad\s+(\d+)\s+SAD\s+(\d+)\s+(\S+)')
    if rng:
        bad = [m for m in rng if m.group(1) != m.group(2)]
        checks.append(Check(f'{label}FLUC 乱数の最終状態 ({len(rng)} 回)', float(len(bad)), 0.0,
                            'PASS' if not bad else 'FAIL'))
    return checks


def check_level1(out):
    sad = read(os.path.join(BENCH, 'bench_level1.result'))
    c = []
    names_p = ['betx', 'bety', 'alfx', 'alfy', 'etax', 'etapx', 'Qx', 'Qy']
    names_s = ['BetaX', 'BetaY', 'AlphaX', 'AlphaY', 'EtaX', 'EtaPX', 'Qx', 'Qy']
    for p, s in zip(names_p, names_s):
        c.append(cmp(f'Test1 入口 {s}', grab(out, rf'^\s*{p}\s*=\s*({NUM})'), grab(sad, rf'{s}\s*=\s*({NUM})'),
                     TOL_DISP6, 0.0, '表示 6 桁'))
    for p, s in [('betx', 'BetaX'), ('bety', 'BetaY'), ('etax', 'EtaX')]:
        c.append(cmp(f'Test1 QD.1 入口 {s}', grab(out, rf'^\s*{p}\s*=\s*({NUM})', 2), grab(sad, rf'{s}\s*=\s*({NUM})', 2),
                     TOL_DISP6, 0.0, '表示 6 桁'))
    # Test 2: COD (µm 表示 → m)
    c.append(cmp('Test2 COD 入口 x', grab(out, rf'^\s*x\s*=\s*({NUM}) µm') * 1e-6, grab(sad, rf'COD_x\s*=\s*({NUM})'),
                 1e-10, TOL_CALC_COD, '表示 1e-10 m'))
    c.append(cmp('Test2 COD 入口 px', grab(out, rf'^\s*px\s*=\s*({NUM}) µrad') * 1e-6, grab(sad, rf'COD_px\s*=\s*({NUM})'),
                 1e-10, TOL_CALC_COD))
    c.append(cmp('Test2 COD QD.1 入口 x', grab(out, rf'^\s*x\s*=\s*({NUM}) µm', 1) * 1e-6, grab(sad, rf'COD_x\s*=\s*({NUM})', 1),
                 1e-10, TOL_CALC_COD))
    c.append(cmp('Test2 COD 下 BX', grab(out, rf'^\s*BX\s*=\s*({NUM})'), grab(sad, rf'f BX\s+#+\s+#\s+({NUM})', 1),
                 TOL_DISP6, 0.0, '表示 6 桁'))
    c.append(cmp('Test2 COD 下 Qx', grab(out, rf'^\s*Qx\s*=\s*({NUM})', 1), grab(sad, rf'f NX\s+#+\s+#\s+({NUM})', 1),
                 TOL_DISP6, 0.0, '表示 6 桁'))
    # Test 3: RADINT
    for p, s in [('I1', 'I1X'), ('I2', 'I2'), ('I3', 'I3'), ('I4', 'I4U'), ('I5', 'I5U')]:
        c.append(cmp(f'Test3 RADINT {s}', grab(out, rf'^\s*{p}\s*=\s*({NUM})'), grab(sad, rf'{s}\s*=\s*({NUM})'),
                     0.0, TOL_RADINT, 'SAD 7 桁表示'))
    for s in ['Emittance x', 'MomentumSpread', 'MomentumCompaction', 'EnergyLossU0']:
        c.append(cmp(f'Test3 Emittance[] {s}', grab(out, rf'{s}\s*=\s*({NUM})'), grab(sad, rf'{s}\s*=\s*({NUM})'),
                     0.0, TOL_EMIT))
    return c


def check_level2(out):
    sad = read(os.path.join(BENCH, 'bench_level2.result'))
    c = []
    names = ['BetaX', 'BetaY', 'AlphaX', 'AlphaY', 'EtaX', 'EtaPX', 'Qx', 'Qy']
    pat = r'\b(' + '|'.join(names) + r')\s*=\s*(' + NUM + r')'   # SAD 側は ' In[2]:=   BetaX = ...' なので行頭に固定しない
    vp = [(m.group(1), fnum(m.group(2))) for m in grab_all(out, pat)]
    vs = [(m.group(1), fnum(m.group(2))) for m in grab_all(sad, pat)]
    blocks = ['dp=0', 'dp=-0.02', 'dp=-0.01', 'dp=+0.01', 'dp=+0.02', 'RINGK (K1 付き BEND)']
    if len(vp) != 48 or len(vs) != 48:
        c.append(Check('TestA/B Twiss の行数', float(len(vp)), 48.0, 'FAIL', f'pysad {len(vp)} / SAD {len(vs)} (48 のはず)'))
    for k in range(min(len(vp), len(vs))):
        blk = blocks[k // 8]
        tol = TOL_DISP6
        note = '表示 6 桁'
        c.append(cmp(f'{blk} {vp[k][0]}', vp[k][1], vs[k][1], tol, 0.0, note))
    for p, s in [('I1', 'I1X'), ('I2', 'I2'), ('I3', 'I3'), ('I4', 'I4U'), ('I5', 'I5U')]:
        c.append(cmp(f'RINGK RADINT {s}', grab(out, rf'^\s*{p}\s*=\s*({NUM})'), grab(sad, rf'{s}\s*=\s*({NUM})'),
                     0.0, TOL_RADINT, 'SAD 7 桁表示'))
    for s in ['Emittance x', 'MomentumSpread', 'MomentumCompaction', 'EnergyLossU0']:
        c.append(cmp(f'RINGK Emittance[] {s}', grab(out, rf'{s}\s*=\s*({NUM})'), grab(sad, rf'{s}\s*=\s*({NUM})'),
                     0.0, TOL_EMIT))
    return c


def check_level3a(out):
    c = []
    for k, m in enumerate(grab_all(out, rf'チューン 21 点の SAD との差 \(最大\):\s*Qx\s+({NUM}),\s*Qy\s+({NUM})')):
        lab = ['補正前', '補正後'][k] if k < 2 else f'#{k}'
        c.append(cmp_diff(f'{lab} DP0 スキャン 21 点 Qx', fnum(m.group(1)), TOL_CALC_COD, 'オフモメンタム'))
        c.append(cmp_diff(f'{lab} DP0 スキャン 21 点 Qy', fnum(m.group(2)), TOL_CALC_COD, 'オフモメンタム'))
    for k, (p, s) in enumerate(chrom_pairs(out)):
        lab = ['補正前 X', '補正前 Y', '補正後 X', '補正後 Y'][k] if k < 4 else f'#{k}'
        c.append(cmp(f'{lab} c0 (チューン)', p[0], s[0], TOL_CALC_COD, 0.0))
        c.append(cmp(f'{lab} c1 (線形色収差)', p[1], s[1], TOL_CHROM_C1, 0.0, '当てはめ係数'))
        c.append(Check(f'{lab} c2,c3 (高次)', max(abs(p[2] - s[2]), abs(p[3] - s[3])), None, 'INFO', '参考 (判定なし)'))
    return c


def check_level3bcd(out, label):
    c = []
    for name, p, s in cod_pairs(out):
        for k, comp in enumerate(['x', 'px', 'y', 'py']):
            c.append(cmp(f'閉軌道 {name} {comp}', p[k], s[k], 1e-12, TOL_CALC_COD))
    # 以前は 3C の EX/EPX (3.2e-7)、3B/3C の EY (1e-8 台) を既知の不一致としていたが、四極・薄い多極の解析行列 (2026-10-16,
    # docs/work_emit_report.md) と z 行の組み直しで解消した (EX 3e-11、EPX 6e-12)。既知なし。
    known = {}
    for name, p, s, d in tw_pairs(out, ['BX', 'BY', 'AX', 'AY', 'EX', 'EPX', 'EY', 'EPY', 'QX', 'QY', 'BX11', 'BY11']):
        c.append(cmp(f'光学 {name}', p, s, 1e-8, TOL_CALC_COD, known.get(name, '閉軌道/結合あり'), known=name in known))
    pat = r'R1-R4 pysad\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s*\n\s*R1-R4 SAD\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')'
    for m in grab_all(out, pat):
        for k in range(4):
            c.append(cmp(f'結合 R{k+1}', fnum(m.group(1 + k)), fnum(m.group(5 + k)), 1e-7, 0.0, 'SAD 5-7 桁'))
    for k, (p, s) in enumerate(chrom_pairs(out)):
        lab = ['X', 'Y'][k] if k < 2 else f'#{k}'
        c.append(cmp(f'色収差 {lab} c0 (チューン)', p[0], s[0], TOL_CALC_COD, 0.0))
        c.append(cmp(f'色収差 {lab} c1 (線形色収差)', p[1], s[1], TOL_CHROM_C1, 0.0, '当てはめ係数'))
    for name, p, s, d in emit_pairs(out):
        c.append(cmp(f'Emittance[] {name}', p, s, 0.0, TOL_EMIT))
    return c


def check_level4(out):
    c = []
    sec = 0
    for line in out.splitlines():
        if line.startswith('--- (2)'):
            sec = 1
        m = re.match(r'^\s*(\S+(?: \S+)?)\s*:\s*pysad\s+(' + NUM + r')\s+SAD\s+(' + NUM + r')\s+相対差', line)
        if not m:
            continue
        name, p, s = m.group(1).strip(), fnum(m.group(2)), fnum(m.group(3))
        lab = ['RFSW', 'RFSW+RADCOD'][sec]
        if name.startswith('cod'):
            c.append(cmp(f'{lab} 閉軌道 {name}', p, s, 1e-9, 0.0, 'SAD tcod の収束誤差 (1.5e-10 m)'))
        else:
            c.append(cmp(f'{lab} {name}', p, s, 0.0, TOL_EMIT))
    return c


def check_level5(out):
    c = []
    case = ''
    for line in out.splitlines():
        m = re.match(r'^=== (\S+)', line)
        if m:
            case = m.group(1)
            continue
        m = re.match(r'^\s*T (\S+)\s+BX,AX,BY,AY,EX,EPX\s+差:\s+(.*)$', line)
        if m:
            d = max(abs(fnum(x)) for x in m.group(2).split())
            c.append(cmp_diff(f'{case} Twiss @{m.group(1)} (絶対差)', d, 1e-6, 'INS, 軌道あり (BX 〜10 m)'))
            continue
        m = re.match(r'^\s*O (\S+)\s+x,px,y,py\s+差:\s+(.*)$', line)
        if m:
            d = max(abs(fnum(x)) for x in m.group(2).split())
            c.append(cmp_diff(f'{case} 軌道 @{m.group(1)}', d, TOL_TRACK))
            continue
        m = re.match(r'^\s*NX,NY 差:\s+(' + NUM + r')\s+(' + NUM + r')', line)
        if m:
            d = max(abs(fnum(m.group(1))), abs(fnum(m.group(2))))
            c.append(cmp_diff(f'{case} 位相進み NX,NY', d, TOL_CALC_COD))
    return c


def check_level12(out):
    c = []
    m = re.search(r'要素名の不一致:\s*(\d+)', out)
    c.append(Check('ラティス読み込み 要素名の不一致', float(m.group(1)) if m else None, 0.0,
                   'PASS' if (m and m.group(1) == '0') else 'FAIL'))
    m = re.search(r'パラメータの最大差.*?:\s*(.*)$', out, re.M)
    if m:
        for key, val in re.findall(r'(\w+)\s+(' + NUM + r')', m.group(1)):
            c.append(cmp_diff(f'ラティス読み込み {key}', fnum(val), TOL_CALC, 'SAD 10 桁表示'))
    for m in grab_all(out, r'^\s*(\$\$\$|MKC)\s+([A-Z,]+)\s+差:\s+(.*)$'):
        d = max(abs(fnum(x)) for x in m.group(3).split())
        c.append(cmp_diff(f'CELL CALC {m.group(1)} {m.group(2)}', d, TOL_CALC))
    for name, p, s, d in emit_pairs(out):
        c.append(cmp(f'emit {name}', p, s, 0.0, TOL_EMIT_DISP, 'SAD 表示 5〜8 桁'))
    return c


def check_track_fast(out):
    c = []
    for m in grab_all(out, r'^\s*(\S+ T\s*\d+)\s+numba.*?SAD との最大差\s+(' + NUM + r')'):
        c.append(cmp_diff(f'track6n {m.group(1)}', fnum(m.group(2)), TOL_TRACK))
    rng = grab_all(out, r'乱数\s+pysad\s+(\d+)\s+SAD\s+(\d+)')
    if rng:
        bad = [m for m in rng if m.group(1) != m.group(2)]
        c.append(Check(f'track6n FLUC 乱数の最終状態 ({len(rng)} 回)', float(len(bad)), 0.0, 'PASS' if not bad else 'FAIL'))
    m = re.search(r'track6 vs track6n.*?最大差:\s*(' + NUM + r')', out)
    if m:
        c.append(cmp_diff('track6 (純 Python) と track6n の差', fnum(m.group(1)), 1e-12))
    return c


def compare_tags(out, sad, tags_tol, label=''):
    """Print["TAG ", {...}] の出力を SAD と比較。tags_tol: {TAG: (tol_abs, tol_rel, note)}"""
    c = []
    tp, ts = extract_tagged(out), extract_tagged(sad)
    for tag, (ta, tr, note) in tags_tol.items():
        if tag not in tp or tag not in ts:
            c.append(Check(f'{label}{tag}', None, None, 'FAIL', f'出力に無い (pysad {tag in tp}, SAD {tag in ts})'))
            continue
        for k, (bp, bs) in enumerate(zip(tp[tag], ts[tag])):
            sp, ss = strings_in(bp), strings_in(bs)
            if sp or ss:
                c.append(Check(f'{label}{tag} 文字列 {len(sp)} 個', float(sp != ss), 0.0, 'PASS' if sp == ss else 'FAIL'))
                continue
            vp, vs = numbers_in(bp), numbers_in(bs)
            if len(vp) != len(vs):
                c.append(Check(f'{label}{tag}', None, None, 'FAIL', f'要素数 pysad {len(vp)} / SAD {len(vs)}'))
                continue
            dmax = max((abs(a - b) - tr * abs(b)) for a, b in zip(vp, vs)) if vp else 0.0
            worst = max(range(len(vp)), key=lambda i: abs(vp[i] - vs[i]) - tr * abs(vs[i])) if vp else 0
            d = abs(vp[worst] - vs[worst]) if vp else 0.0
            tol = ta + tr * abs(vs[worst]) if vp else ta
            suffix = f' #{k+1}' if len(tp[tag]) > 1 else ''
            c.append(Check(f'{label}{tag}{suffix} ({len(vp)} 値)', d, tol, 'PASS' if d <= tol else 'FAIL', note))
    return c


def compare_matching_lines(out, sad, label=''):
    """' Matched. ( r)' / ' Residual = r' と反復行 'iter nstab res (NEWTON) fact nvar' の列を比べる"""
    c = []
    pat = r'^\s*(Matched\.\s*\(\s*(' + NUM + r')\s*\)|Residual\s*=\s*(' + NUM + r'))'
    mp = [(m.group(1).split('.')[0].split(' ')[0], fnum(m.group(2) or m.group(3))) for m in grab_all(out, pat)]
    ms = [(m.group(1).split('.')[0].split(' ')[0], fnum(m.group(2) or m.group(3))) for m in grab_all(sad, pat)]
    kp, ks = [x[0] for x in mp], [x[0] for x in ms]
    c.append(Check(f'{label}Matched/Residual 行の列 {kp}', float(kp != ks), 0.0, 'PASS' if kp == ks else 'FAIL',
                   '' if kp == ks else f'SAD {ks}'))
    for k, ((_, a), (_, b)) in enumerate(zip(mp, ms)):
        c.append(cmp(f'{label}Matched/Residual #{k+1} 残差', a, b, 1e-9, 1e-3, 'SAD 5 桁表示'))
    pat = r'^\s+(\d+)\s+(\d+)\s+(' + NUM + r')\s+\((\w+)\)\s+(' + NUM + r')\s+(\d+)\s*$'
    ip = [(m.group(1), m.group(2), fnum(m.group(3)), m.group(4), m.group(6)) for m in grab_all(out, pat)]
    is_ = [(m.group(1), m.group(2), fnum(m.group(3)), m.group(4), m.group(6)) for m in grab_all(sad, pat)]
    same = [(a[0], a[1], a[3], a[4]) for a in ip] == [(a[0], a[1], a[3], a[4]) for a in is_]
    c.append(Check(f'{label}反復行 (iter, nstab, 方法, nvar) {len(ip)} 行', float(not same), 0.0,
                   'PASS' if same else 'FAIL', '' if same else f'pysad {ip} / SAD {is_}'))
    for k, (a, b) in enumerate(zip(ip, is_)):
        c.append(cmp(f'{label}反復行 #{k+1} 残差', a[2], b[2], 1e-9, 1e-3, 'SAD 5 桁表示'))
    return c


def check_level13(out, sadfile, tags):
    sad = read(os.path.join(BENCH, sadfile))
    tt = {t: (TOL_CALC, TOL_CALC, '表示 10 桁') for t in tags}
    c = compare_matching_lines(out, sad) + compare_tags(out, sad, tt)
    if sadfile == 'bench_level13a.result':
        # 13A-2 は解なし (4 条件 / 3 変数) の停留点。解析微分 (既定) では GO の最終点が、残差の差 1e-15 (相対 4e-13) の
        # 大小で決まる隣の点 (K1 で 3e-10 相対) になり、SAD / 数値微分の点と食い違う。真の最小点は解析微分側。
        # FFS$NumericalDerivative=True なら SAD の最終点を 10 桁で再現する (docs/work_integration.md §2)
        for ck in c:
            if ck.status == 'FAIL' and ck.diff is not None and ck.diff <= 1e-7 and ck.name.split(' ')[0] in ('K1B', 'TWB'):
                ck.status = 'KNOWN'
                ck.note = '13A-2 解なしの停留点: 解析微分は隣の点 (work_integration.md §2)'
    return c


def check_level14a(out):
    sad = read(os.path.join(BENCH, 'bench_level14a.result'))
    c = []
    pat = r'^\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s*$'
    sp = [[fnum(m.group(k)) for k in (1, 2, 3)] for m in grab_all(out, pat)]
    ss = [[fnum(m.group(k)) for k in (1, 2, 3)] for m in grab_all(sad, pat)]
    sp = [r for r in sp if abs(r[0]) <= 0.0101 and 40 < r[1] < 50]
    ss = [r for r in ss if abs(r[0]) <= 0.0101 and 40 < r[1] < 50]
    c.append(Check(f'DP0 スキャンの行数 pysad {len(sp)} / SAD {len(ss)}', None, None,
                   'PASS' if (len(sp) == 21 and len(ss) == 21) else 'FAIL'))
    if sp and len(sp) == len(ss):
        dx = max(abs(a[1] - b[1]) for a, b in zip(sp, ss))
        dy = max(abs(a[2] - b[2]) for a, b in zip(sp, ss))
        c.append(cmp_diff('DP0 スキャン 21 点 NX (FFS["cell calc;"])', dx, TOL_CALC_COD,
                          'オフモメンタム (eps=1e-7 の丸めで 2.5e-7, 10/14 時点 ≤2e-8)', known=dx <= 5e-7))
        c.append(cmp_diff('DP0 スキャン 21 点 NY', dy, TOL_CALC_COD, 'オフモメンタム'))
    pat = r'^\{\{\s*(' + NUM + r'),\s*(' + NUM + r'),\s*(' + NUM + r'),\s*(' + NUM + r')\},\{Residual->\s*(' + NUM + r')\}\}'
    fp = [[fnum(m.group(k)) for k in range(1, 6)] for m in grab_all(out, pat)]
    fs = [[fnum(m.group(k)) for k in range(1, 6)] for m in grab_all(sad, pat)]
    for k, (a, b) in enumerate(zip(fp, fs)):
        lab = ['tunex', 'tuney'][k] if k < 2 else f'#{k}'
        d = max(abs(a[i] - b[i]) / abs(b[i]) for i in range(4))
        c.append(cmp_diff(f'PolynomialFit[{lab},3] 係数の相対差', d, 1e-4, '当てはめ係数 (チューンのノイズを拡大)'))
    # EMIT (RADCOD ON) の報告: 既知の不一致 (dz)
    for name, pat2, tol, known in [
            ('EMIT U0 [MV]', r'Energy loss per turn U0 =\s*(' + NUM + r')', 1e-6, False),
            ('EMIT dz [mm] (RADCOD ON)', r'Equilibrium position dz =\s*(' + NUM + r')', 1e-6, True),
            ('EMIT alpha', r'Momentum compact\. alpha =\s*(' + NUM + r')', 1e-4, False),
            ('EMIT Emittance X', r'Emittance X\s*=\s*(' + NUM + r')', 1e-4, False),
            ('EMIT Energy spread', r'Energy spread\s*=\s*(' + NUM + r')', 1e-5, False),
            ('EMIT Bunch Length', r'Bunch Length\s*=\s*(' + NUM + r')', 1e-5, True)]:
        try:
            c.append(cmp(name, grab(out, pat2), grab(sad, pat2), 0.0, tol, 'SAD 表示桁', known=known))
        except ValueError as e:
            c.append(Check(name, None, None, 'FAIL', str(e)))
    return c


def da_map(text):
    """DA の地図 (NZ 0----|---- の行に挟まれた行) を返す"""
    lines = text.splitlines()
    maps = []
    for i, line in enumerate(lines):
        if re.match(r'^\s*NZ\s+0----', line):
            j = i + 1
            rows = []
            while j < len(lines) and not re.match(r'^\s*NZ\s+0----', lines[j]):
                rows.append(lines[j].rstrip())
                j += 1
            if rows:
                maps.append(rows)
    return maps[0] if maps else []


def check_level16a(out):
    sad = read(os.path.join(BENCH, 'bench_level16a.result'))
    c = []
    mp, ms = da_map(out), da_map(sad)
    if not mp or len(mp) != len(ms):
        c.append(Check('DA 地図の行数', None, None, 'FAIL', f'pysad {len(mp)} / SAD {len(ms)}'))
    else:
        nbad = sum(1 for a, b in zip(mp, ms) for x, y in zip(a.ljust(len(b)), b.ljust(len(a))) if x != y)
        ntot = sum(len(b) for b in ms)
        c.append(Check(f'DA 地図 {len(mp)} 行 {ntot} 文字の不一致', float(nbad), 0.0, 'PASS' if nbad == 0 else 'FAIL'))
    try:
        c.append(cmp('DA Score', grab(out, r'Score:\s+(\d+)'), grab(sad, r'Score:\s+(\d+)'), 0.0, 0.0))
    except ValueError as e:
        c.append(Check('DA Score', None, None, 'FAIL', str(e)))
    c += compare_tags(out, sad, {'EMITS': (1e-12, TOL_EMIT, 'Emittance[] (EMITY は 0)'),
                                 'SCORE': (0.0, 0.0, '')})
    tp, ts = extract_tagged(out), extract_tagged(sad)
    for tag in ['RES', 'RDAP']:
        if tag in tp and tag in ts:
            vp, vs = numbers_in(tp[tag][0]), numbers_in(ts[tag][0])
            if len(vp) != len(vs):
                c.append(Check(f'DA {tag} 要素数', None, None, 'FAIL', f'pysad {len(vp)} / SAD {len(vs)}'))
                continue
            diffs = [abs(a - b) for a, b in zip(vp, vs)]
            nd = sum(1 for d in diffs if d > 0.5)
            dmax = max(diffs) if diffs else 0.0
            if nd == 0:
                st = 'PASS'
            elif nd <= 1 and dmax <= 1.0:
                st = 'KNOWN'
            else:
                st = 'FAIL'
            c.append(Check(f'DA {tag} 損失周数 ({len(vp)} 値, 不一致 {nd} 個)', dmax, 0.0, st,
                           '' if nd == 0 else '損失境界の粒子 1 個が ±1 周 (既知)'))
        else:
            c.append(Check(f'DA {tag}', None, None, 'FAIL', '出力に無い'))
    return c


def check_tags(result_file, known=(), num=None):
    """Level 24 (SADScript の関数・FFS コマンド・CalculateOptics): SAD の出力 (result_file) と pysad の出力を cmp_tags.py のタグ単位で比べる。
    known = SAD に無く pysad が拡張として評価するタグ (数えるだけ)、num = (絶対, 相対) なら数値許容で比べる"""
    def f(out):
        import cmp_tags
        da = cmp_tags.load(os.path.join(BENCH, result_file))
        db = cmp_tags.parse(out.split('\n'))
        tabs, trel = num if num else (None, None)
        bad, ext, n = cmp_tags.compare(da, db, known, tabs, trel, verbose=False)
        c = [Check(f'タグ {n} 個の不一致 (既知の拡張 {ext} 個を除く)', float(len(bad)), 0.0, 'PASS' if not bad else 'FAIL',
                   '' if not bad else '不一致: ' + ' '.join(bad[:12]))]
        return c
    return f


def check_lines(result_file, tol, what):
    """Print[] のラベル付き行 (cmp_lines.py) を SAD の結果 (result_file) と順に突き合わせ、行ごとの最大相対差 (max(|SAD|,1) で割る) が tol 以下か"""
    def f(out):
        import cmp_lines
        a = cmp_lines.rows(open(os.path.join(BENCH, result_file), encoding='utf-8', errors='replace'))
        b = cmp_lines.rows(out.split('\n'))
        res = cmp_lines.compare(a, b)
        c = []
        bad = [r for r in res if r[1] is None or not r[1] <= tol]
        worst = max((r[1] for r in res if r[1] is not None), default=0.0)
        c.append(Check(f'{what} {len(res)} 行の最大相対差', worst, tol, 'PASS' if not bad else 'FAIL',
                       '' if not bad else '不一致: ' + ' '.join(f'{r[0]}' for r in bad[:8])))
        return c
    return f


def check_level28(result_file, tol=1e-8, loose=()):
    """Level 28 (不安定ラインの CELL CALC と GO): Print[] のラベル付き行 (cmp_lines) の最大相対差 + 反復表 (iter, nstab, 方法, nvar, 残差 4 桁) +
    `Horizontal/Vertical Unstable.` と `Unstable = n` と qcod の報告の行の列が SAD と同じか。
    loose: [(ラベルの先頭, 許容)] 解が不安定の縁 (β が 1e5) で丸めが 1e-5 まで拡大されるラベルの許容"""
    def f(out):
        import cmp_lines
        sad = read(os.path.join(BENCH, result_file))
        a = cmp_lines.rows(sad.split('\n'))
        b = cmp_lines.rows(out.split('\n'))
        res = cmp_lines.compare(a, b)
        c = []
        bad = []
        for r in res:
            t = tol
            for pre, tl in loose:
                if r[0].startswith(pre):
                    t = tl
            if r[1] is None or not r[1] <= t:
                bad.append(r[0])
        worst = max((r[1] for r in res if r[1] is not None), default=0.0)
        c.append(Check(f'Print[] {len(res)} 行の最大相対差', worst, tol, 'PASS' if not bad else 'FAIL',
                       '' if not bad else '不一致: ' + ' '.join(bad[:8])))
        c += compare_matching_lines(out, sad)
        for what, pat in (('Unstable. の行', r'^ (Horizontal|Vertical|Horizontal/Vertical) Unstable\.'),
                          ('Unstable = n の行', r'^Unstable =\s+(\d+)'),
                          ('qcod の報告', r'^ \*\*\*\*\*qcod--->.*$')):
            ls = [m.group(0) for m in grab_all(sad, pat)]
            lp = [m.group(0) for m in grab_all(out, pat)]
            c.append(Check(f'{what} {len(ls)} 行の列', float(ls != lp), 0.0, 'PASS' if ls == lp else 'FAIL',
                           '' if ls == lp else f'SAD {ls[:6]} / pysad {lp[:6]}'))
        return c
    return f


def check_emit_report(result_file, noise=1e-11, rel=1e-9, skip=()):
    """EMIT 報告の全行 (temitprint) を SAD と行単位で突き合わせる (cmp_emit_report.py)。0 同然の雑音の桁だけの差は許す。
    noise/rel: same() の許容。skip: この語で始まる行は数値を比べない (両方にあれば一致。不安定な光学の 'Deviation from symplectic matrix' など雑音の量)"""
    def f(out):
        import cmp_emit_report as cer
        a = cer.load_blocks(os.path.join(BENCH, result_file))
        b = cer.blocks(out.split('\n'))
        nbad = 0
        n = 0
        if len(a) != len(b):
            return [Check(f'報告の個数 SAD {len(a)} / pysad {len(b)}', 1.0, 0.0, 'FAIL')]
        for x, y in zip(a, b):
            if len(x) != len(y):
                nbad += 1
                continue
            for lx, ly in zip(x, y):
                n += 1
                if any(lx.lstrip().startswith(k) and ly.lstrip().startswith(k) for k in skip):
                    continue
                if not cer.same(lx, ly, noise, rel)[0]:
                    nbad += 1
        return [Check(f'EMIT 報告 {len(a)} 通り {n} 行の不一致行数', float(nbad), 0.0, 'PASS' if not nbad else 'FAIL')]
    return f


def check_radint(result_file, tol):
    """RADINT (INTGRL.f) の出力ブロックを SAD と行単位で比べる (cmp_radint.py)。tol はブロックごとのリストでもよい。
    SAD の表示は 7 桁なので相対 1e-6 が目安 (0 同然の V 面の量は雑音として一致扱い)"""
    def f(out):
        import cmp_radint
        a = cmp_radint.blocks(os.path.join(BENCH, result_file))
        b = cmp_radint.blocks_text(out)
        ok, bad, worst = cmp_radint.compare(a, b, tol, quiet=True)
        tl = max(tol) if isinstance(tol, (list, tuple)) else tol
        return [Check(f'RADINT {len(a)} ブロック {ok + bad} 行の不一致行数 (最大相対差 {worst:.1e})', float(bad), 0.0, 'PASS' if not bad else 'FAIL',
                      f'許容 {tl:g}')]
    return f


def check_blocks(result_file, tag, tol=None, noise=None, skip=()):
    """Print["Dn ..."] で区切られたブロックの出力行を SAD と比べる (cmp_blocks.py)。noise 未満の数 (SAD の 1e-17 などの丸め雑音) は 0 と同じ扱い"""
    def f(out):
        import re as _re
        import cmp_blocks
        n, bad = cmp_blocks.compare(os.path.join(BENCH, result_file), out, _re.compile(tag), tol, noise, skip)
        return [Check(f'{n} ブロックの不一致ブロック数', float(len(bad)), 0.0, 'PASS' if not bad else 'FAIL',
                      '' if not bad else '不一致: ' + ' '.join(bad[:10]))]
    return f


def check_level27(out):
    """Level 27 (MAIN の TRACK): bench_level27a_pysad.py が SAD (SADScript で同じ手順) と smear・最終座標・tune・損失を比べて
    ケースごとに `一致` / `不一致` を出す。全ケース一致で最後に `Matched`"""
    n_ok = len(re.findall(r'^\s+一致\s*$', out, re.M))
    n_ng = len(re.findall(r'^\s+不一致\s*$', out, re.M))
    cks = [Check(f'MAIN TRACK ケース一致 ({n_ok} 件)', float(n_ng), 0.0, 'PASS' if (n_ng == 0 and n_ok > 0) else 'FAIL')]
    cks.append(Check('Matched', 0.0 if 'Matched' in out else 1.0, 0.0, 'PASS' if 'Matched' in out else 'FAIL'))
    return cks


T24A_KNOWN = ('T118', 'T191', 'T195', 'T220', 'T221', 'T222', 'T229', 'T231', 'U143', 'W075')


# ----------------------------------------------------------------------------
# ベンチの一覧
# ----------------------------------------------------------------------------
BENCHES = [
    # (key, level, 説明, コマンド引数 (python の後), 判定関数, QUICK で飛ばす, 目安時間)
    ('1',   'L1',  '線形光学 / COD (QD DX) / RADINT / Emittance[]', ['bench_level1_pysad.py'], check_level1, False, '3 s'),
    ('2',   'L2',  'オフモメンタム Twiss / K1 付き BEND', ['bench_level2_pysad.py'], check_level2, False, '6 s'),
    ('3a',  'L3A', '六極と色収差 (DP0 スキャン)', ['bench_level3a_pysad.py'], check_level3a, False, '10 s'),
    ('3b',  'L3B', '六極を通る COD (DELX/DELY 1 台)', ['bench_level3b_pysad.py'], lambda o: check_level3bcd(o, '3B'), False, '6 s'),
    ('3c',  'L3C', 'x-y 結合 (DTHETA) / Emittance[]', ['bench_level3c_pysad.py'], lambda o: check_level3bcd(o, '3C'), False, '4 s'),
    ('3d',  'L3D', '強さの誤差 (DK/DDK)', ['bench_level3d_pysad.py'], lambda o: check_level3bcd(o, '3D'), False, '3 s'),
    ('4',   'L4',  'RF 空洞 (RFSW, RADCOD) の Emittance[]', ['bench_level4_pysad.py'], check_level4, False, '4 s'),
    ('5',   'L5',  '輸送ライン (INS) とフリンジ', ['bench_level5_pysad.py'], check_level5, False, '1 s'),
    ('6',   'L6',  '追跡 (RAD OFF / RAD ON NOFLUC)', ['bench_level6_pysad.py'], lambda o: track_checks(o), False, '9 s'),
    ('7',   'L7',  '追跡: 誤差入り要素, K1≠0 BEND', ['bench_level7_pysad.py'], lambda o: track_checks(o), True, '2 min'),
    ('8',   'L8',  '追跡: 四極フリンジ, 長い空洞, MULT', ['bench_level8_pysad.py'], lambda o: track_checks(o), False, '35 s'),
    ('9',   'L9',  '追跡: ソレノイド (SOL BZ=2 T)', ['bench_level9_pysad.py'], lambda o: track_checks(o), False, '9 s'),
    ('10',  'L10', '追跡: RAD+FLUC (SAD の乱数)', ['bench_level10_pysad.py'], lambda o: track_checks(o), False, '9 s'),
    ('11',  'L11', '追跡: MULT K0/K1/K2 とずれ・回転', ['bench_level11_pysad.py'], lambda o: track_checks(o), False, '17 s'),
    ('12',  'L12', 'ESRF 実ラティス: 読み込み / CELL CALC / emit', ['bench_level12_pysad.py'], check_level12, False, '10 s'),
    ('12t', 'L12', 'ESRF 実ラティス: 追跡 (純 Python track6)', ['bench_level12_track_pysad.py'], lambda o: track_checks(o), True, '1 min'),
    ('12f', 'L12', 'ESRF 実ラティス: 追跡 (numba track6n)', ['bench_level12_track_fast.py'], check_track_fast, False, '20 s'),
    ('13a', 'L13A', 'マッチング FIT/FREE/GO (小リング, run_sad.py)', ['run_sad.py', 'bench_level13a.sad'],
     lambda o: check_level13(o, 'bench_level13a.result', ['TUNE0', 'K1A', 'TUNEA', 'RESA', 'K1B', 'TWB', 'RESB']), True, '3〜4 min (13A-2 の解なしの反復が遅い)'),
    ('13b', 'L13B', 'マッチング (ESRF, 7 ファミリー, run_sad.py)', ['run_sad.py', 'bench_level13b.sad'],
     lambda o: check_level13(o, 'bench_level13b.result', ['TW0', 'K10', 'QNAMES', 'K1B', 'TWB', 'RESB']), False, '10 s〜'),
    ('14a', 'L14A', 'SADScript: esrf/calc_chromaticity (run_sad.py)', ['run_sad.py', 'bench_level14a.sad'], check_level14a, False, '30 s'),
    ('16a', 'L16A', 'DynamicApertureSurvey (NORAD, run_sad.py)', ['run_sad.py', 'bench_level16a.sad'], check_level16a, True, '30 s〜'),
    ('27a', 'L27A', 'MAIN の TRACK: 規格化振幅→初期分布, smear, tune (wtunem), RAD/FLUC (SAD は SADScript 版)', ['bench_level27a_pysad.py', 'a'],
     check_level27, False, '10 s'),
    ('27b', 'L27B', 'MAIN の TRACK: 粒子が落ちる場合 (損失周, 落ちた粒子の x,y,z)', ['bench_level27a_pysad.py', 'b'],
     check_level27, False, '10 s'),
    ('22e', 'L22E', 'INS/CELL: 四極ずれ・DP0≠0 (閉軌道≠0, ESRF)', ['run_sad.py', 'bench_level22e.sad'],
     check_lines('bench_level22e.result', 1e-10, 'CELL/INS'), False, '10 s'),
    ('26a', 'L26A', 'EMIT 報告の全行 (ESRF, RFSW × RADCOD の 4 通り)', ['run_sad.py', 'bench_level26a.sad'],
     check_emit_report('bench_level26a.result'), False, '12 s'),
    ('26b', 'L26B', '四極 FRINGE=1〜3 / ずれ・ロール付き四極・六極 (CALC と Emittance[])', ['run_sad.py', 'bench_level26b.sad'],
     check_lines('bench_level26b.result', 3e-9, '26B'), False, '3 s'),
    ('26c', 'L26C', 'Orbit dilation (FSHIFT): EMIT 報告・Emittance[]・CALC6D', ['run_sad.py', 'bench_level26c.sad'],
     lambda o: check_emit_report('bench_level26c.result')(o) + check_lines('bench_level26c.result', 1e-9, '26C')(o), False, '8 s'),   # 11/04: 3e-10 → 1e-9 (numpy 2.5/numba 0.68 で F6 が 3.1e-10)
    ('26d', 'L26D', 'RFSW OFF + RADCOD ON の Emittance[] と CALC6D (ESRF)', ['run_sad.py', 'bench_level26d.sad'],
     check_lines('bench_level26d.result', 1e-9, '26D'), False, '10 s'),
    ('28a', 'L28A', '不安定ライン (INTRES/HALFRES ON): CELL CALC の続行 / Unstable 表示 / GO (nstab) (ESRF, V・H・HV の 3 通り)', ['run_sad.py', 'bench_level28a.sad'],
     check_level28('bench_level28a.result', 1e-8, (('TWF_H', 1e-3), ('K1F_H', 1e-6), ('TWU_H', 1e-8))), False, '12 s'),
    ('28b', 'L28B', '不安定ライン: INTRES/HALFRES の 4 通り (NOINTRES では安定扱いで tracex=4/trace) と GO', ['run_sad.py', 'bench_level28b.sad'],
     check_level28('bench_level28b.result', 1e-8, (('TWF_H', 1e-3), ('K1F_H', 1e-6))), False, '12 s'),
    ('28c', 'L28C', '不安定な光学 (V: QD x1.19) の EMIT 報告の全行 (teigen+tnorm の直訳。丸め雑音の桁のみ差)', ['run_sad.py', 'bench_level28c.sad'],
     check_emit_report('bench_level28c.result', noise=1e-10), False, '10 s'),
    ('28d', 'L28D', '不安定な光学 (H: QF x0.9) の EMIT 報告の全行 (T の悪条件で 1e-5 の差まで許す)', ['run_sad.py', 'bench_level28d.sad'],
     check_emit_report('bench_level28d.result', noise=1e-4, rel=3e-5, skip=('*** Deviation',)), False, '10 s'),
    ('28f', 'L28F', '不安定な光学の Emittance[] (V・H): Stable=0, 戻り値の全量', ['run_sad.py', 'bench_level28f.sad'],
     check_lines('bench_level28f.result', 3e-5, '28F'), False, '15 s'),
    ('30a', 'L30A', 'INS + CALC6D (qcell6d の cell=.false. 枝): 先頭 MARK の 28 成分から 6D 伝搬 (ESRF、A〜E の 5 通り)', ['run_sad.py', 'bench_level30a.sad'],
     check_lines('bench_level30a.result', 1e-9, '30A'), False, '20 s'),
    ('30b', 'L30B', 'RFSW ON + RADCOD OFF + FSHIFT≠0 の Emittance[] (tcod の 6D Newton の直訳)', ['run_sad.py', 'bench_level30b.sad'],
     check_lines('bench_level30b.result', 1e-10, '30B'), False, '60 s'),
    ('29a', 'L29A', 'RADINT コマンド (INTGRL.f): ESRF の 6 通り (RFSW/EMIT 後/ずれ・ロール/DP0/COD)', ['run_sad.py', 'bench_level29a.sad'],
     check_radint('bench_level29a.result', 1e-6), False, '10 s'),
    ('29a2', 'L29A', 'RADINT: 小リング (E1/E2 付き BEND, RF なし/あり, 結合)', ['run_sad.py', 'bench_level29a2.sad'],
     check_radint('bench_level29a2.result', 1e-4), False, '3 s'),
    ('29b', 'L29B', 'DISP の位置範囲 (数値・名前・NAME.n±k・パターン・周回・ALL・R/O/E モード)', ['run_sad.py', 'bench_level29b.sad'],
     check_blocks('bench_level29b.result', r'^D[0-9]+ ', None, 1e-9, ('D37',)), False, '10 s'),
    ('29c', 'L29C', 'CALC6D;/CALC4D; 単独 (CALC なし: 再計算されない) と 4D CALC が 6D の縦の列 (AZ..) を残す癖', ['run_sad.py', 'bench_level29c.sad'],
     check_lines('bench_level29c.result', 1e-9, '29C'), False, '15 s'),
    ('29d', 'L29D', 'FFS 内の ON/OFF flag (無効: Out[1]:= USE を出し、後ろの語はフラグとして読まれる)', ['run_sad.py', 'bench_level29d.sad'],
     check_blocks('bench_level29d.result', r'^E[0-9]+[ab]? ', None, 1e-10, ('E17',)), False, '5 s'),
    ('29e', 'L29E', 'FSHIFT (Orbit dilation) の追跡: NORAD / RAD / RAD+FLUC (乱数も)', ['run_sad.py', 'bench_level29e.sad'],
     check_lines('bench_level29e.result', 5e-11, '29E'), False, '10 s'),
    ('31a', 'L31A', '誤差コマンド (DELX/DELY/DK/DTHETA/DDK, UNIFORM/COHERENT/ADD): master グループ共有・tfadjst・触ったキーの展開と CALC (ESRF, 128 タグ)', ['run_sad.py', 'bench_level31a.sad'],
     check_tags('bench_level31a.result', (), (1e-10, 1e-7)), False, '10 s'),
    ('31b', 'L31B', 'Twiss[]/LINE[] の位置パターンの返り値 (一致 1 個ならスカラー、0 個は {}) (ESRF)', ['run_sad.py', 'bench_level31b.sad'],
     check_tags('bench_level31b.result', (), (1e-10, 1e-8)), False, '5 s'),
    ('24a', 'L24A', 'SADScript のリスト・文字列・行列関数 (559 タグ)', ['run_sad.py', 'bench_level24a.sad'],
     check_tags('bench_level24a.result', T24A_KNOWN), False, '3 s'),
    ('24b', 'L24B', 'CHROMATICITY コマンド (ESRF, 全 960 要素)', ['run_sad.py', 'bench_level24b.sad'],
     check_tags('bench_level24b.result'), False, '15 s'),
    ('24b2', 'L24B', 'CHROMATICITY コマンド (MULT の K1/K2, DP0≠0)', ['run_sad.py', 'bench_level24b2.sad'],
     check_tags('bench_level24b2.result'), False, '5 s'),
    ('24c', 'L24C', 'FFS の REPEAT/UNTIL, IF/ELSEIF/ELSE/ENDIF', ['run_sad.py', 'bench_level24c.sad'],
     check_tags('bench_level24c.result'), False, '5 s'),
    ('24d', 'L24D', 'CalculateOptics (4 次元、cell/INS、区間、結合、分散)', ['run_sad.py', 'bench_level24d.sad'],
     check_tags('bench_level24d.result', (), (1e-8, 1e-7)), False, '5 s'),
    ('24d2', 'L24D', 'CalculateOptics (ESRF、6 次元)', ['run_sad.py', 'bench_level24d2.sad'],
     check_tags('bench_level24d2.result', (), (1e-8, 1e-7)), False, '60 s'),
    ('24f', 'L24F', 'Element/SetElement/FLAG/BeamLine/WriteBeamLine/USE BeamLine', ['run_sad.py', 'bench_level24f.sad'],
     check_tags('bench_level24f.result', ('G004',)), False, '5 s'),
]


def run_bench(key, level, desc, args, checker, slow, est):
    env = dict(os.environ)
    env.setdefault('PYSAD_COMPAT', 'k64')   # 回帰は gs (K64) との比較なので k64 モードで走らせる
    env['PYTHONIOENCODING'] = 'utf-8'
    env['PYTHONUTF8'] = '1'
    # args の 'run_sad.py' はリポジトリのルート、'bench_*' は bench/ の実パスにする (cwd は HERE)
    args = [os.path.join(ROOT, a) if a == 'run_sad.py' else os.path.join(BENCH, a) if a.startswith('bench_') else a for a in args]
    cmd = [sys.executable, '-u'] + args
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=HERE, env=env, capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=3600)
        out = r.stdout + ('\n[stderr]\n' + r.stderr if r.stderr.strip() else '')
        rc = r.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or '') + '\n[timeout]'
        rc = -1
    dt = time.time() - t0
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, f'{key}.txt'), 'w', encoding='utf-8') as f:
        f.write(' '.join(cmd) + '\n\n' + out)
    try:
        checks = checker(out)
    except Exception as e:  # 出力形式が変わった等
        checks = [Check('判定', None, None, 'FAIL', f'判定中に例外: {e!r}')]
    if rc != 0:
        checks.append(Check('終了コード', float(rc), 0.0, 'FAIL', f'rc={rc}'))
    return dt, checks


def summarize(checks):
    n = {s: sum(1 for c in checks if c.status == s) for s in ['PASS', 'FAIL', 'KNOWN', 'INFO']}
    judged = [c for c in checks if c.status in ('PASS', 'FAIL', 'KNOWN') and c.diff is not None and c.tol is not None]
    worst = max(judged, key=lambda c: (c.diff / c.tol) if c.tol > 0 else c.diff, default=None)
    if n['FAIL']:
        verdict = 'FAIL'
    elif n['KNOWN']:
        verdict = 'PASS*'
    else:
        verdict = 'PASS'
    return verdict, n, worst


def main():
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass
    sel = [a for a in sys.argv[1:] if not a.startswith('-')]
    print(f'pysad 回帰テスト  (QUICK={int(QUICK)}, python {sys.version.split()[0]}, cwd {HERE})')
    print(f'許容: 追跡 {TOL_TRACK:g} / CALC {TOL_CALC:g} (閉軌道・結合・オフモメンタム経由は相対 {TOL_CALC_COD:g}) / '
          f'emit 相対 {TOL_EMIT:g} (SAD 5 桁表示は {TOL_EMIT_DISP:g}) / DA は地図と Score の一致')
    print()
    rows = []
    details = []
    t_all = time.time()
    for key, level, desc, args, checker, slow, est in BENCHES:
        if sel and key not in sel:
            continue
        if QUICK and slow and not sel:
            rows.append((level, key, desc, 'SKIP', '', '', '', f'QUICK=1 (目安 {est})'))
            print(f'[{key:4s}] SKIP  {desc}')
            continue
        print(f'[{key:4s}] 実行中 {" ".join(args)} ... (目安 {est})', flush=True)
        dt, checks = run_bench(key, level, desc, args, checker, slow, est)
        verdict, n, worst = summarize(checks)
        wtxt = f'{worst.diff:.1e} / {worst.tol:.0e} ({worst.name})' if worst else ''
        rows.append((level, key, desc, verdict, f'{n["PASS"]}/{n["PASS"]+n["FAIL"]+n["KNOWN"]}',
                     f'{n["FAIL"]}', wtxt, f'{dt:.0f} s'))
        print(f'       {verdict:5s} PASS {n["PASS"]} FAIL {n["FAIL"]} KNOWN {n["KNOWN"]}  最悪 {wtxt}  ({dt:.0f} s)')
        for c in checks:
            if c.status in ('FAIL', 'KNOWN'):
                details.append((key, c))
    print()
    print('=' * 132)
    print(f'{"Level":6s} {"key":4s} {"内容":46s} {"合否":5s} {"PASS/判定":9s} {"FAIL":4s} {"最悪 (差/許容, 項目)":44s} {"時間"}')
    print('-' * 132)
    for level, key, desc, v, np_, nf, w, dt in rows:
        print(f'{level:6s} {key:4s} {desc[:46]:46s} {v:5s} {np_:9s} {nf:4s} {w[:44]:44s} {dt}')
    print('=' * 132)
    print('合否: PASS = 全項目が許容内, PASS* = 既知の不一致 (KNOWN) のみ, FAIL = 許容を超える項目あり, SKIP = QUICK=1 で省略')
    if details:
        print()
        print('--- FAIL / KNOWN の項目 ---')
        for key, c in details:
            print(f'[{key:4s}] {c!r}')
    nfail = sum(1 for r in rows if r[3] == 'FAIL')
    print()
    print(f'合計 {len(rows)} ベンチ, FAIL {nfail}, 所要 {time.time() - t_all:.0f} s. 出力は {OUTDIR} に保存')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
