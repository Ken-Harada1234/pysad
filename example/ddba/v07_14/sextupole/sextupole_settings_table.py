# -*- coding: utf-8 -*-
"""六極の設定をいくつか選んで、RDT・ADTS・クロマティシティの数字を並べた表 (R7、最優先の依頼)。
使い方: py -3.12 sextupole_settings_table.py [--da] [--src-dir 元のラティス定義のフォルダ]    (--da: 追跡の DA も求める。1 設定に 1〜2 分)
出力: results/settings_table.json、results/settings_table.md、results/settings_table.out.txt、fig_settings.png
選ぶ設定 (六極の種類):
  (A) 別の研究者による最適化の結果 (全セルの種類が独立に自由な六極): v07_14 の 20231205 再調整版 (final)、v07_16 の再調整 (rev9)、第 1 次の再調整 (retune1)、
      DDBA_11 (v05_14 に入れた MH・LH の 10 ファミリーの結果。光学は v05_14)
  (B) 全セルを自由にしていない (ファミリーでまとめた、設計者自身の探索の) 台本の値: v07_11、v07_12 (v07_14 の光学に載せ替え)、v07_12 そのまま (DDBA_22)、
      v07_13 (元のフォルダがあるとき。チューンが違う)、DDBA_05 (v05_12。S1 = −S2 など)、DDBA_10 (v05_14)
  (C) 参考: 六極を全部 0 にしたもの (v07_14 の光学。DDBA_05・DDBA_10 の光学)
同じ光学 (v07_14) に六極だけ載せ替えたものと、元の光学のままのもの (表の「光学」の列) を分ける。
計算: resdrive.resonance_terms (1 次、δ の項。偏向磁石の 3 次の項込み)、rdt2nd.rdt4 (4 次、六極の 2 次込み)、resdrive.adts(fourth='r3'、optics='calc')、
resdrive.chromaticity_natural + 六極 + 偏向磁石の δ の項、chromaticity_numeric (CALC、order 3)。単位は docs/49 と同じ (h は (2J)^{n/2} で規格化、ADTS は 1/(m rad))。
"""
import argparse
import glob
import io
import json
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import rdt2nd as Q
from pysad import resdrive as R
from pysad import sxtopt as S
from pysad.parser import load_sad_file, load_sad_string

ap = argparse.ArgumentParser()
ap.add_argument('--da', action='store_true')
ap.add_argument('--src-dir', default=os.environ.get('V07_14_SRC_DIR'))
args = ap.parse_args()
LOG = io.StringIO()


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')


LATDIR = os.path.join(C.ROOT, 'example', 'ddba', 'lattice')
CUB = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
DEL = ('h11001', 'h00111', 'h20001', 'h00201', 'h10002')
H4 = ('h22000', 'h11110', 'h00220', 'h31000', 'h40000', 'h20020', 'h20200', 'h11200', 'h00310', 'h00400', 'h20110', 'h11020')


def momentum_of(src):
    m = re.search(r'MOMENTUM\s*=\s*([0-9.]+)\s*(GEV|MEV)', src, re.I)
    v = float(m.group(1)) * (1e9 if m.group(2).upper() == 'GEV' else 1e6)
    return v


def pattern(names, k):
    """セルの種類ごとの S1〜S5 の K2L の並びから、ファミリーでまとめたか全部独立かを判定して説明する"""
    groups = {}
    for n, v in zip(names, k):
        m = re.match(r'^S(\d)(\w*)$', n)
        if m:
            groups.setdefault(m.group(2), {})[int(m.group(1))] = float(v)
    desc, free = [], True
    for suf, d in groups.items():
        idx = sorted(d)
        rel = []
        unused = [i for i in range(1, 6) if i not in d]
        if unused:
            rel.append(','.join('S%d' % i for i in unused) + ' は 0 (未使用)')
            free = False
        zeros = [i for i in idx if abs(d[i]) < 1e-9]
        if zeros:
            rel.append('S' + '=S'.join(map(str, zeros)) + '=0')
            free = False
        nz = [i for i in idx if i not in zeros]
        for a_ in range(len(nz)):
            for b_ in range(a_ + 1, len(nz)):
                i, j = nz[a_], nz[b_]
                if abs(d[i] - d[j]) < 1e-6 * max(abs(d[i]), 1e-9):
                    rel.append('S%d=S%d' % (i, j))
                    free = False
                elif abs(d[i] + d[j]) < 1e-6 * max(abs(d[i]), 1e-9):
                    rel.append('S%d=-S%d' % (i, j))
                    free = False
        desc.append('%s: %s' % (suf or '-', ', '.join(rel) if rel else '5 個とも独立'))
    # 同じ強さを持つ別のセルの種類 (LH = MH など)
    return ('全ファミリー独立' if free else 'ファミリーでまとめた'), '; '.join(desc[:4]) + (' …' if len(desc) > 4 else '')


def analyze(lat, p0, with_calc=True):
    out = {}
    names, fam, k = S.sextupole_families(lat)
    sx_idx = {i for i, e in enumerate(lat.elements) if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0}
    r = R.resonance_terms(lat, CUB + DEL, momentum_eV=p0, bend_cubic=True, order2=False)
    out['h3'] = {t: abs(r.total[t]) for t in CUB}
    out['hdel'] = {t: abs(r.total[t]) for t in DEL}
    mask = np.array([i in sx_idx for i in r.elem_idx])
    canc, mx = {}, {}
    for it, t in enumerate(r.names):
        col = r.per_element[:, it]
        s_col = col[mask] if mask.any() else np.zeros(1, dtype=complex)
        tot = abs(col[mask].sum()) if mask.any() else 0.0
        den = float(np.sum(np.abs(s_col)))
        canc[t] = tot / den if den > 0 else float('nan')
        mx[t] = float(np.max(np.abs(s_col))) if mask.any() else 0.0
    out['canc'], out['elmax'] = canc, mx
    r4 = Q.rdt4(lat, order2=True, momentum_eV=p0)
    out['h4'] = {t: abs(r4.total[t]) for t in H4}
    a = R.adts(lat, momentum_eV=p0, fourth='r3', optics='calc')
    out['adts'] = [float(a.total['xx']), float(a.total['xy']), float(a.total['yy'])]
    out['adts_parts'] = {q: [float(getattr(a, q)[c]) for c in ('xx', 'xy', 'yy')] for q in ('sext', 'octu', 'kin', 'fringe', 'bend', 'cross')}
    nat = R.chromaticity_natural(lat)
    rs = R.resonance_terms(lat, ('h11001', 'h00111'), momentum_eV=p0, bend_cubic=True, order2=False)
    xb = (rs.total_bend['h11001'].real / math.pi, rs.total_bend['h00111'].real / math.pi)
    xs = (rs.total_sext['h11001'].real / math.pi, rs.total_sext['h00111'].real / math.pi)
    out['xi_nat'], out['xi_bend'], out['xi_sext'] = list(nat), list(xb), list(xs)
    out['xi_an'] = [nat[0] + xb[0] + xs[0], nat[1] + xb[1] + xs[1]]
    if with_calc:
        for dm in (2e-3, 5e-4):
            try:
                ch = R.chromaticity_numeric(lat, order=3, dmax=dm)
                out['xi_calc'], out['xi2_calc'] = list(ch['d1']), list(ch['d2'])
                out['nu'] = list(ch['nu0'])
                out['chrom_dmax'] = dm
                break
            except Exception as ex:
                out['xi_calc'] = None
                out['calc_err'] = str(ex)
    out['k2l_max'] = float(np.max(np.abs(k))) if len(k) else 0.0
    out['nfam'], out['nsext'] = len(names), sum(len(v) for v in fam.values())
    out['pattern'] = pattern(names, k) if len(k) else ('(六極なし)', '')
    return out


SETTINGS = []     # (名前, 種類 A/B/C, ラベル, 光学, src 文字列, ring, 運動量)
src_v = C.load_source()
sets = C.load_sets()
nm_fam = sets['family_names']
LAB = {'final': ('A', 'v07_14 の 20231205 再調整版 (最終)'), 'rev9': ('A', 'v07_16 の再調整 (rev9)'), 'retune1': ('A', '第 1 次の再調整 (2023-10)'),
       'v07_11': ('B', 'v07_11 の六極 (S1=S2=0、S3=S5 の 3 個)'), 'v07_12': ('B', 'v07_12 の六極 (同上、M-A だけ少し違う)')}
for key in ('final', 'rev9', 'retune1', 'v07_12', 'v07_11'):
    k = np.array(sets['sets'][key]['k'])
    SETTINGS.append((key, LAB[key][0], LAB[key][1], 'v07_14 (載せ替え)', S.replace_k2_in_source(src_v, nm_fam, k), 'RING', 5e9))
SETTINGS.append(('zero_v07_14', 'C', '六極を全部 0 (v07_14 の光学)', 'v07_14', S.replace_k2_in_source(src_v, nm_fam, np.zeros(len(nm_fam))), 'RING', 5e9))
native = [('DDBA_11', 'A', 'DDBA_11 (v05_14 に入れた MH・LH 10 ファミリーの結果。2023-08)', 'DDBA_11_v05_14_sx.sad'),
          ('DDBA_22', 'B', 'DDBA_22 = v07_12 そのまま (元の並び)', 'DDBA_22_v07_12.sad'),
          ('DDBA_05', 'B', 'DDBA_05 (v05_12、六極の総当たり探索の版)', 'DDBA_05_v05_12.sad'),
          ('DDBA_10', 'B', 'DDBA_10 (v05_14)', 'DDBA_10_v05_14.sad')]
for key, kind, lab, fn in native:
    path = os.path.join(LATDIR, fn)
    if os.path.exists(path):
        SETTINGS.append((key, kind, lab, '元の光学 (%s)' % fn.split('_', 2)[2].replace('.sad', ''), open(path, encoding='utf-8', errors='replace').read(), 'RING', None))
if args.src_dir and os.path.isdir(args.src_dir):
    p13 = os.path.join(args.src_dir, 'v07_13')
    if os.path.exists(p13):
        SETTINGS.append(('v07_13', 'B', 'v07_13 の六極 (元の光学。チューン 45.865/21.665 で違う)', '元の光学 (v07_13)', open(p13, encoding='utf-8', errors='replace').read(), 'RING', None))
for base in ('DDBA_05', 'DDBA_10'):
    for s_ in list(SETTINGS):
        if s_[0] == base:
            lat0 = load_sad_string(s_[4]).build_lattice(s_[5])
            names0, _, _ = S.sextupole_families(lat0)
            try:
                SETTINGS.append(('zero_' + base, 'C', '六極を全部 0 (%s の光学)' % base, s_[3], S.replace_k2_in_source(s_[4], names0, np.zeros(len(names0))), s_[5], None))
            except ValueError as ex:
                P('注意: %s の六極を 0 にできなかった (%s)' % (base, ex))

rows = {}
for name, kind, label, optics, src, ring, p0 in SETTINGS:
    lat = load_sad_string(src).build_lattice(ring)
    p0 = p0 or momentum_of(src)
    try:
        o = analyze(lat, p0)
    except Exception as ex:
        P('%s: 解析に失敗 (%s: %s)' % (name, type(ex).__name__, ex))
        continue
    o.update(name=name, kind=kind, label=label, optics=optics, p0_GeV=p0 / 1e9)
    rows[name] = o
    P('%-12s [%s] %s: ADTS = %s、ξ(CALC) = %s、|h30000| = %.2f' % (name, kind, label, np.round(o['adts'], 0).tolist(), np.round(o.get('xi_calc') or [np.nan] * 2, 3).tolist(), o['h3']['h30000']))

# 追跡の DA
if args.da:
    from sxt_track import da_rays, da_xdp
    for name, kind, label, optics, src, ring, p0 in SETTINGS:
        if name not in rows:
            continue
        path = os.path.join(C.OUT, 'track_%s.json' % name)
        if name in ('final', 'rev9', 'retune1', 'v07_11', 'v07_12') and os.path.exists(path):
            d = json.load(open(path, encoding='utf-8'))
            rows[name]['da'], rows[name]['xdp'] = d['da'], d['xdp']
            continue
        try:
            st = C.make_track_setup(src, rad=False, rfsw=False, alost=0.05, ring=ring)
            rows[name]['da'] = da_rays(st, 1, 512, False)
            rows[name]['xdp'] = da_xdp(st, 1, 512, False)
            P('%s: DA 面積 %.1f mm²、x_DA %.2f mm、y_DA %.2f mm、x-δ 幅 %.1f mm (リング先頭 βx=%.1f βy=%.1f)' % (
                name, rows[name]['da']['area_mm2'], rows[name]['da']['x_da_mm'], rows[name]['da']['y_da_mm'], rows[name]['xdp']['mean_width_mm'], *rows[name]['da']['beta']))
        except Exception as ex:
            P('%s: DA に失敗 (%s: %s)' % (name, type(ex).__name__, ex))
os.makedirs(C.OUT, exist_ok=True)
json.dump(rows, open(os.path.join(C.OUT, 'settings_table.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

# 表 (markdown)
ML = []
ML.append('### 表 1: 設定と六極の種類')
ML.append('| 名前 | 種類 | 設定 | 光学 | 運動量 GeV | ファミリー数 / 台数 | 全ファミリー独立か | K2L の並び (セルの種類ごと) | max|K2L| 1/m² |')
ML.append('|---|---|---|---|---|---|---|---|---|')
for n, o in rows.items():
    ML.append('| %s | %s | %s | %s | %g | %d / %d | %s | %s | %.1f |' % (n, o['kind'], o['label'], o['optics'], o['p0_GeV'], o['nfam'], o['nsext'], o['pattern'][0], o['pattern'][1], o['k2l_max']))
ML.append('')
ML.append('### 表 2: 1 次の RDT |h| (六極の 3 次、偏向磁石込み) と δ の項、要素ごとの寄与の相殺 (六極だけの Σ/Σ|寄与|、最大の 1 台の寄与)')
ML.append('| 名前 | h21000 | h30000 | h10110 | h10020 | h10200 | h11001 | h00111 | h20001 | h00201 | h10002 | 相殺 h21000 | h30000 | h10020 | 最大寄与 h21000 | h30000 | h10020 |')
ML.append('|' + '---|' * 17)
for n, o in rows.items():
    ML.append('| %s | %s | %s | %s | %s | %s | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g |' % (
        n, *['%.3g' % o['h3'][t] for t in CUB], *[o['hdel'][t] for t in DEL], o['canc']['h21000'], o['canc']['h30000'], o['canc']['h10020'],
        o['elmax']['h21000'], o['elmax']['h30000'], o['elmax']['h10020']))
ML.append('')
ML.append('### 表 3: 4 次の RDT |h| (写像の係数。六極の 2 次 + 八極 + 運動学 + フリンジ + 偏向磁石) と ADTS (fourth=\'r3\'、optics=\'calc\')')
ML.append('| 名前 | ' + ' | '.join(H4) + ' | ADTS ∂νx/∂Jx | ∂νx/∂Jy | ∂νy/∂Jy |')
ML.append('|' + '---|' * (len(H4) + 5))
for n, o in rows.items():
    ML.append('| %s | ' % n + ' | '.join('%.3g' % o['h4'][t] for t in H4) + ' | %.4g | %.4g | %.4g |' % tuple(o['adts']))
ML.append('')
ML.append('### 表 4: クロマティシティ (六極込みの ξ、自然・偏向磁石の δ の項・六極の内訳、CALC)')
ML.append('| 名前 | ξx 自然 | ξy 自然 | ξx 偏向磁石 | ξy 偏向磁石 | ξx 六極 | ξy 六極 | ξx 解析 | ξy 解析 | ξx CALC | ξy CALC | ∂²νx/∂δ² (CALC) | ∂²νy/∂δ² | チューン (CALC) |')
ML.append('|' + '---|' * 14)
for n, o in rows.items():
    c = o.get('xi_calc') or [float('nan')] * 2
    c2 = o.get('xi2_calc') or [float('nan')] * 2
    nu = o.get('nu') or [float('nan')] * 2
    ML.append('| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.0f | %.0f | %.3f, %.3f |' % (
        n, *o['xi_nat'], *o['xi_bend'], *o['xi_sext'], *o['xi_an'], *c, *c2, *nu))
if args.da:
    ML.append('')
    ML.append('### 表 5: 追跡の DA (簡略化: 512 周、RF・放射なし、δ=0、誤差なし、LOSSAMPL ±50 mm、リング先頭から。光学が違う設定は β が違うので DA の mm は直接比べられない)')
    ML.append('| 名前 | DA 面積 mm² | x_DA mm | y_DA mm | x-δ の幅 (平均) mm | β (x, y) 先頭 |')
    ML.append('|---|---|---|---|---|---|')
    for n, o in rows.items():
        if 'da' in o:
            ML.append('| %s | %.1f | %.2f | %.2f | %.1f | %.1f, %.1f |' % (n, o['da']['area_mm2'], o['da']['x_da_mm'], o['da']['y_da_mm'], o['xdp']['mean_width_mm'], *o['da']['beta']))
open(os.path.join(C.OUT, 'settings_table.md'), 'w', encoding='utf-8').write('\n'.join(ML) + '\n')
P('\n' + '\n'.join(ML))
open(os.path.join(C.OUT, 'settings_table.out.txt'), 'w', encoding='utf-8').write(LOG.getvalue())

# 図: 設定ごとの主な項 (v07_14 の光学の載せ替えだけ)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
    same = [n for n in ('zero_v07_14', 'v07_11', 'v07_12', 'retune1', 'final', 'rev9') if n in rows]
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (ttl, fn) in zip(axs, (('1 次の RDT |h| (5 項)', lambda o: [o['h3'][t] for t in CUB]),
                                   ('4 次の RDT |h|', lambda o: [o['h4'][t] for t in H4[:8]]),
                                   ('ADTS の絶対値 [1/(m rad)]', lambda o: [abs(v) for v in o['adts']]),
                                   ('ξ (CALC) の絶対値と |δ の項|', lambda o: [abs(o['xi_calc'][0]), abs(o['xi_calc'][1]), o['hdel']['h20001'], o['hdel']['h00201']]))):
        lab = {'1 次の RDT |h| (5 項)': CUB, '4 次の RDT |h|': H4[:8], 'ADTS の絶対値 [1/(m rad)]': ('xx', 'xy', 'yy'), 'ξ (CALC) の絶対値と |δ の項|': ('ξx', 'ξy', 'h20001', 'h00201')}[ttl]
        w = 0.8 / len(same)
        for i, n in enumerate(same):
            v = np.array(fn(rows[n]), dtype=float)
            ax.bar(np.arange(len(v)) + i * w, np.maximum(v, 1e-3), w, label=n)
        ax.set_yscale('log')
        ax.set_xticks(np.arange(len(lab)) + 0.4)
        ax.set_xticklabels(lab, rotation=60, fontsize=7)
        ax.set_title(ttl, fontsize=9)
    axs[0].legend(fontsize=7)
    fig.suptitle('v07_14 の光学での六極の組の比較 (zero_v07_14 = 六極なし)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'fig_settings.png'), dpi=120)
except Exception as ex:
    P('図に失敗:', ex)
