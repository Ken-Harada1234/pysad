# -*- coding: utf-8 -*-
"""example/ddba/lattice/DDBA_NN_<台本名>.sad を作る (元の台本 + 整えたヘッダ)。
使い方: DDBA_SRC=<元の台本のあるフォルダ> py -3.12 example/tools/ddba/build_lattices.py <results_k64.json> <results_v2.json>
  results_*.json は run_all.py (src モード) の出力。ラティスの定義の中身は変えない。
  変えるのは (1) 先頭のヘッダ、(2) 改行を LF に、(3) コメントの中の氏名 (1 行。「! from <名前> san」の形の行) を中立な英語のコメントに。"""
import sys, os, json, re, textwrap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ddba_lib as L, ddba_info as D
PAPER = 'K. Harada et al., "Conceptual design of the Hybrid Ring with superconducting linac", J. Synchrotron Rad. 29, 118-124 (2022), DOI 10.1107/S1600577521012753'
SX_OTHER = (11, 14, 24, 25, 26, 27, 28)   # 六極の設定に、別の研究者による最適化の結果を含む版 (設計メモの記述による)
def wrap(prefix, text, w=48):
    lines = textwrap.wrap(text, width=w, break_long_words=True, break_on_hyphens=False) or ['']
    return [prefix + lines[0]] + ['!     ' + l for l in lines[1:]]
def header(nn, srcname, outname, res, res2):
    hidden = '*' in srcname or '[' in srcname
    title = outname if hidden else srcname
    srcname = '(元の台本名は氏名を含むため省略)' if hidden else srcname
    inf = D.I[nn]; sz = inf['SZ']
    h = ['! ' + '=' * 70, '! DDBA_%02d  %s' % (nn, title)]
    h += wrap('! 設計手順の段階: ', '%d. %s' % (inf['STAGE'], D.STAGES[inf['STAGE']]))
    h += wrap('! 元の台本: ', '%s (設計メモの添付 %s KB、元の台本 %d byte。%s)' % (srcname, sz[0], sz[1], sz[2]))
    h += ['! 設計メモのページの最終更新: %s' % (inf['DATE'] or '不明')]
    h += wrap('! 何の版か (設計メモの記述の要約): ', inf['SUM'])
    h += wrap('! 設計メモが書いている数値: ', inf['MEMO'])
    if res.get('C'):
        h += ['! pysad の cell calc emit (ON RAD FLUC COD RFSW RADCOD、リングのライン RING、k64 モード = gs と同じ):',
              '!     エネルギー %.1f GeV、周長 %.4f m、チューン (%.4f, %.4f)' % (res['P'] / 1e9, res['C'], res['nx'], res['ny']),
              '!     εx = %.4g m rad (k64)、%.4g m rad (v2 モード)、σE = %.4g、αc = %.4g、U0 = %.4g MV' % (res['ex'], res2['ex'], res['sige'], res['alpha'], res['U0']),
              '!     RF: 総電圧 %.3g MV、h = %d、周波数 %.6g MHz。偏向磁石 %d 台 (L > 0, |ANGLE| > 0.02)、逆曲げなどの小角 BEND %d 台、四極 %d 台、六極 %d 台、八極 %d 台' % (
                  res['V'] / 1e6, res['h'], res['f'] / 1e6, res['nb'], res['nr'], res['nq'], res['ns'], res['no'])]
    else:
        h += wrap('! pysad の cell calc emit: ', '走らなかった')
    if nn in SX_OTHER:
        h += wrap('! 六極の設定: ', 'pysad の作者とは別の研究者 (公開論文 ' + PAPER + ' の著者の一人) の手による最適化の結果を含む。この論文には六極の K2 の値や最適化の手法の記述はなく、論文の暫定設計 (3 GeV、周長 350 m、16 セル) はこの版とは別の版。論文の設定と一致する、同じ最適化の結果である、という意味ではない')
    h += ['! 元のヘッダ: 元の台本の先頭にコメントのヘッダはない (先頭は MOMENTUM の行)。',
          '! 変更: 改行を LF にした。ラティスの定義の中身は変えていない。', '! ' + '=' * 70, '']
    return h
if __name__ == '__main__':
    k64 = json.load(open(sys.argv[1])); v2 = json.load(open(sys.argv[2]))
    os.makedirs(L.LAT, exist_ok=True)
    for nn, srcname, outname in L.MAP:
        body = open(L.find_src(srcname), encoding='utf-8').read().replace('\r\n', '\n').replace('\r', '\n')
        lines = body.split('\n'); new = []
        for ln in lines:
            if re.match(r'\s*!\s*from \w+ san', ln):
                ln = '! sextupole set (optimization result by another researcher)'
            new.append(ln)
        txt = '\n'.join(header(nn, srcname, outname, k64[str(nn)], v2[str(nn)])) + '\n'.join(new)
        if not txt.endswith('\n'): txt += '\n'
        open(L.lat_path(nn), 'w', encoding='utf-8', newline='\n').write(txt)
        print(nn, os.path.basename(L.lat_path(nn)))
