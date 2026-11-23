#!/usr/bin/env python3
"""
bench_level21a_pysad.py — Level 21a: Emittance[OneTurnInformation->True, Matrix->True] の全規則と
LINE["SIG",l] (各要素入口の平衡ビーム行列 Σ(s)) を pysad.emit_extra.emittance_rules と照合する。
対応する SAD スクリプト: bench_level21a.sad (pysad/tests/kek-ls_v03_12.lat; CELL; RFSW; RAD; CALC6D; CALC)
使い方: PYTHONIOENCODING=utf-8 python -u bench_level21a_pysad.py
目標: 相対 1e-9 (絶対下限は量ごとに表示)
"""
import os
import re
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.emit_extra import emittance_rules, TWISS_FUNS, beam_size_at

s = open(os.path.join(HERE, 'bench_level21a.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def _parse(txt):
    """SAD の Print 出力 ({...} の入れ子) を Python のリストにする"""
    txt = txt.strip()
    txt = txt.replace('{', '[').replace('}', ']')
    txt = re.sub(r'(?<![\w.])True(?![\w])', '1.0', txt)
    txt = re.sub(r'(?<![\w.])False(?![\w])', '0.0', txt)
    txt = re.sub(r'(\d)E([+-]?\d)', r'\1e\2', txt)

    def tok(m):
        t = m.group(1).strip('"')
        try:
            return repr(float(t))          # SAD は '-.4199' のように先頭の 0 を省く
        except ValueError:
            return '"' + t + '"'
    txt = re.sub(r'([^\[\],\s]+)', tok, txt)
    import json
    return json.loads(txt)


def sad(label):
    """'CO100' のように英字+番号のときは SAD が番号を $FORM で印字する ('CO 100.0000000000000 {...}')"""
    m = re.match(r'^([A-Z]+)(\d+)$', label)
    if m and m.group(1) in ('CO', 'TM', 'SZ'):
        pat = r'^' + m.group(1) + r'\s+' + m.group(2) + r'\.0*\s+(.*)$'
    else:
        pat = r'^' + re.escape(label) + r'\s+(.*)$'
    m = re.search(pat, s, re.M)
    if m is None:
        raise KeyError(label)
    return _parse(m.group(1))


# --- pysad ---
sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
t0 = time.time()
e = emittance_rules(lat, 3.0e9, rfsw=True, radcod=False)
print(f'emittance_rules: {time.time() - t0:.1f} s')

worst = []


def cmp(tag, mine, ref, abs_floor=1e-12):
    mine = np.asarray(mine, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if mine.shape != ref.shape:
        print(f'{tag}: 形が違う pysad {mine.shape} SAD {ref.shape}')
        worst.append((tag, np.inf))
        return
    d = np.abs(mine - ref)
    # 平面ラティスなので y 成分は 0 の丸め誤差 (SAD 1e-21, pysad 1e-17 程度) しかない。
    # 配列の最大値の 1e-6 を下限にして、ノイズ成分が相対差を支配しないようにする (配列ノルムに対する差 ≤ 1e-6 x 相対差)
    scale = np.maximum(np.abs(ref), max(abs_floor, 1e-6 * float(np.abs(ref).max()) if ref.size else abs_floor))
    rel = d / scale
    i = np.unravel_index(np.argmax(rel), rel.shape) if rel.size else ()
    # 合否は「配列の最大値に対する最大絶対差」(ノルム相対差) で判定する (要素ごとの相対差は参考)
    nrel = float(d.max() / max(float(np.abs(ref).max()), abs_floor))
    print(f'{tag:8s} ノルム相対差 {nrel:.3e}  要素最大相対差 {rel.max():.3e} (絶対差 {d.max():.3e}, 位置 {i}, SAD {ref[i] if rel.size else ""})')
    worst.append((tag, nrel))


print('STB pysad', e['Stable'], 'SAD', sad('STB'))
print('REG pysad', e['Region'], 'SAD', sad('REG'))
cmp('TUN', e['Tunes'], sad('TUN'))
cmp('U0V', [e['EnergyLossU0'], e['RfVoltageVc'], e['EquilibriumPosition'], e['MomentumCompaction'],
            e['OrbitDilation'], e['BucketHeight'], e['HarmonicNumber']], sad('U0V'))
cmp('OEX', e['OrbitAtExit'], sad('OEX'), abs_floor=1.0)   # 閉軌道 ≈ 0 なので絶対差 [m]
cmp('DMP', e['DampingRate'], sad('DMP'))
cmp('EMT', e['Emittances'], sad('EMT'))
cmp('SIG', [e['MomentumSpread'], e['BunchLength']], sad('SIG'))
cmp('TSR', e['TuneShiftByRadiation'], sad('TSR'), abs_floor=1e-9)
cmp('NSP', [e['NominalSpinTune'], e['PolarizationTime'], e['SpinTune'], e['Polarization']], sad('NSP'))
pol = sad('POL')
cmp('POL', [list(e['Polarizations']), list(e['PolarizationVector'])], pol)
cmp('OEN', e['OrbitAtEntrance'], sad('OEN'), abs_floor=1.0)
cmp('OTM', e['OneTurnTransferMatrix'], sad('OTM'), abs_floor=1e-6)
cmp('ODM', e['OneTurnDampingMatrix'], sad('ODM'), abs_floor=1e-9)
cmp('NCO', e['NormalCoordinates'], sad('NCO'), abs_floor=1e-6)
etk = sad('ETK')
print('ETK SAD', etk)
print('ETK ok' if [k.strip('"') for k in etk] == list(TWISS_FUNS) else 'ETK 順が違う!')
etw = sad('ETW')
tw = e['ExtendedTwissParameters']
cmp('ETW', [tw[k] for k in TWISS_FUNS], etw, abs_floor=1e-9)
for k, a, b in zip(TWISS_FUNS, [tw[k] for k in TWISS_FUNS], etw):
    print(f'   {k:5s} pysad {a: .14e}  SAD {b: .14e}  diff {a - b: .3e}')
cmp('OTE', e['OneTurnExcitation'], sad('OTE'), abs_floor=1e-18)
cmp('EBM', e['EquilibriumBeamMatrix'], sad('EBM'), abs_floor=1e-18)
cmp('SPM', e['SpinPrecessionMatrix'], sad('SPM'), abs_floor=1e-9)
nl = sad('NLAT')
print('NLAT SAD', nl, 'pysad', len(e['ClosedOrbit']), len(e['TransferMatrices']), len(e['BeamMatrices']))
nlat = int(nl[0])
for k in (1, 2, 100, 500, 1000, 2000, 3000, nlat):
    cmp(f'CO{k}', e['ClosedOrbit'][k - 1], sad(f'CO{k}'), abs_floor=1.0)   # 絶対差 [m]
for k in (1, 2, 100, 500, 1000, 2000, 3000, nlat):
    cmp(f'TM{k}', e['TransferMatrices'][k - 1], sad(f'TM{k}'), abs_floor=1e-6)
for k in (1, 2, 100, 500, 1000, 2000, 3000, nlat):
    cmp(f'SZ{k}', e['BeamMatrices'][k - 1], sad(f'SZ{k}'), abs_floor=1e-18)
for k in (1, 1000):
    cmp(f'SZB{k}', e['BeamMatrices'][k - 1], sad(f'SZB{k}'), abs_floor=1e-18)
szx = [beam_size_at(e, 100, 'X'), beam_size_at(e, 100, 'XY'), beam_size_at(e, 100, 'PX'), beam_size_at(e, 100, 'DP')]
print('SZX100 pysad', szx, 'SAD', sad('SZX100'))
cmp('SZX100', szx, sad('SZX100'), abs_floor=1e-12)
# 自己検査: 一周後の Σ は平衡ビーム行列に戻るはず
d = np.abs(e['BeamMatrices'][-1] - e['EquilibriumBeamMatrix'])
print(f'自己検査 Σ(nlat) vs emitp: 最大絶対差 {d.max():.3e} (ノルム相対 {d.max() / np.abs(e["EquilibriumBeamMatrix"]).max():.3e})')
print('NAMES SAD', sad('NAMES'), 'pysad', [lat.elements[k - 1].name for k in (1, 2, 100, 500, 1000, 2000, 3000)])
w = max(worst, key=lambda t: t[1])
print(f'最大相対差 (全項目) {w[1]:.3e} @ {w[0]}')
print('OK (<=1e-9)' if w[1] <= 1e-9 else 'NG (>1e-9)')
