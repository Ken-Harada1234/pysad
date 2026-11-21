# -*- coding: utf-8 -*-
"""example/kek_ls の NAME.gs.txt と NAME.pysad.txt を、Print が出す行 (指定の接頭辞の行) の数値で比べる。
使い方 (リポジトリ直下から): py -3.12 example/tools/kek_ls_cmp.py [NAME ...]   (名前なしは全部)
判定: ○ 全数値が相対差 1e-5 以内 (基準は max(|a|,|b|,1e-3))、△ 数値の形は同じだが差がある、× 行数・NaN が違う・pysad が最後まで走らない。
戻り値の辞書は kek_ls_md.py が使う。"""
import glob, os, re, sys
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'kek_ls')
num = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?|NaN')

# 例の名前 → (フォルダ, 比べる行の接頭辞の正規表現)
SPEC = {
    'KEKLS_cell_calc_emit': ('lattice', r'(circumference|number of|tune nx|EMITX|SIGE|SIGZ|MOMENTUM)'),
    'KEKLS_cell_optics': ('lattice', r'(ARCA|MKS |MKC |MKSFA1|E )'),
    'KEKLS_match_cell_tune': ('design', r'(design|target|SF-SF|exmkc|K1:)'),
    'KEKLS_match_chrom': ('design', r'(before|after)'),
    'KEKLS_match_injcell': ('design', r'(ARCA end|design|before|iteration)'),
    'KEKLS_match_ring_geo': ('design', r'(design|start|result)'),
    'KEKLS_match_repeat': ('design', r'[ABCD] '),
    'KEKLS_scale_mini': ('design', r'(design|factor|   K1)'),
    'KEKLS_scale_ncell': ('design', r'(cells|   K1|   ring)'),
    'KEKLS_tune_survey_cell': ('design', r'POINT'),
    'KEKLS_ring_params': ('display', r'(circumference|tune|EMITX|SIGE|SIGZ|xi)'),
    'KEKLS_optics_table': ('display', r'(tune|E )'),
    'KEKLS_calc_H': ('display', r'(Jx|emittance|EMITX|ratio)'),
    'KEKLS_chromaticity': ('display', r'T[01] '),
    'KEKLS_amp_dep_ts': ('display', r'(design tune|AX|AY)'),
    'KEKLS_sextupole_tune_scan': ('display', r'(design|S )'),
    'KEKLS_scan_target': ('display', r'R '),
    'KEKLS_errors_cod': ('errors', r'(design|R )'),
    'KEKLS_correct_cod': ('errors', r'(=====|before|after|iteration|ICX|   COD max|   eta-beat)'),
    'KEKLS_cod_response': ('errors', r'(horizontal|Responsex|Responsey|steering)'),
}


def tagged(path, pat):
    out = []
    if not os.path.exists(path):
        return None
    rx = re.compile('^' + pat)
    for l in open(path, encoding='utf-8', errors='replace'):
        l = l.rstrip('\n')
        if rx.match(l):
            out.append(l)
    return out


def compare(name):
    d, pat = SPEC[name]
    base = os.path.join(EX, d, name)
    g = tagged(base + '.gs.txt', pat); p = tagged(base + '.pysad.txt', pat)
    if g is None or p is None:
        return {'mark': '-', 'note': 'no output', 'n': 0, 'maxrel': None}
    if len(g) != len(p) or not g:
        return {'mark': '×', 'note': 'lines gs %d / pysad %d' % (len(g), len(p)), 'n': len(g), 'maxrel': None}
    worst = 0.0; nn = 0; bad = 0
    for a, b in zip(g, p):
        ta = num.findall(a); tb = num.findall(b)
        if len(ta) != len(tb):
            bad += 1; continue
        for x, y in zip(ta, tb):
            if x == 'NaN' or y == 'NaN':
                if x != y: bad += 1
                continue
            fx, fy = float(x), float(y)
            if abs(fx) > 1e100 or abs(fy) > 1e100:
                if abs(fx) > 1e100 and abs(fy) > 1e100: continue
                bad += 1; continue
            nn += 1
            rel = abs(fx - fy) / max(abs(fx), abs(fy), 1e-3)
            if abs(fx - fy) < 1e-6:
                rel = 0.0   # 絶対差が 1e-6 未満は丸めの雑音
            worst = max(worst, rel)
    if bad:
        return {'mark': '×', 'note': 'token mismatch %d' % bad, 'n': len(g), 'maxrel': worst}
    return {'mark': '○' if worst < 1e-5 else '△', 'note': '', 'n': len(g), 'maxrel': worst}


if __name__ == '__main__':
    names = sys.argv[1:] or list(SPEC)
    for n in names:
        r = compare(n)
        print('%-28s %s  lines %-4s maxrel %s %s' % (n, r['mark'], r['n'], ('%.2e' % r['maxrel']) if r['maxrel'] is not None else '-', r['note']))
