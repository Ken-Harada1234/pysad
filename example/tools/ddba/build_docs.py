# -*- coding: utf-8 -*-
"""README.md・STAGES.md・ORIGINAL_MAP.md の表を作る (results_k64.json など → markdown)。
使い方: py -3.12 example/tools/ddba/build_docs.py  (表は example/ddba/_tables.md に出す。README の本文はここには入れない)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE)); D = os.path.join(EX, 'ddba')
sys.path.insert(0, HERE); import ddba_info as I, ddba_lib as L
R = json.load(open(os.path.join(D, 'results_k64.json'))); V = json.load(open(os.path.join(D, 'results_v2.json'))); C = json.load(open(os.path.join(D, 'chrom_k64.json')))
GS = {1: 5.2096e-10, 6: 1.891e-10, 15: 9.6006e-10, 22: 8.3032e-10, 25: 2.0794e-10, 28: 8.3219e-10}
def fn(nn):
    m = [m for m in L.MAP if m[0] == nn][0]; return 'DDBA_%02d_%s.sad' % (nn, m[2])
def SRC(n): return '(氏名を含むため省略)' if ('*' in n or '[' in n) else n
def short(s, n=60): return s if len(s) <= n else s[:n] + '…'
rows = []
for nn, src, out in L.MAP:
    r = R[str(nn)]; inf = I.I[nn]
    if r.get('C'):
        rows.append('| [%02d](lattice/%s) | %s | %d | %s | %.1f | %.3f | %.4f, %.4f | %.3g | %.3g | %.3g |' % (
            nn, fn(nn), inf['DATE'], inf['STAGE'], SRC(src), r['P'] / 1e9, r['C'], r['nx'], r['ny'], r['ex'] * 1e9, V[str(nn)]['ex'] * 1e9, r['alpha']))
ALL = '| DDBA_NN | 設計メモの日付 | 段階 | 元の台本 | E [GeV] | 周長 [m] | チューン (νx, νy) | εx k64 [nm rad] | εx v2 [nm rad] | αc |\n|---|---|---|---|---|---|---|---|---|---|\n' + '\n'.join(rows)
rows = []
for nn, src, out in L.MAP:
    r = R[str(nn)]
    c = C.get(str(nn))
    rows.append('| %02d | %.4g | %.3f | %.4g | %.3f | %.4g | %.3g | %s | %d / %d / %d / %d / %d | %.3g MV (%.3g MV)、h = %d、%.6g MHz |' % (
        nn, r['sige'], r['sigz'], r['alpha'], r['U0'], r['taux'] * 1e3, r['Jx'], ('%.1f, %.1f' % (c['nx'], c['ny'])) if c else '-',
        r['nb'], r['nr'], r['nq'], r['ns'], r['no'], r['V'] / 1e6, r['Vc'], r['h'], r['f'] / 1e6))
DET = '| DDBA_NN | σE | σz [mm] | αc | U0 [MV] | τx [ms] | Jx | 自然クロマティシティ (x, y) | 偏向 / 小角 BEND / 四極 / 六極 / 八極 の数 | RF 総電圧 (有効電圧)・h・周波数 |\n|---|---|---|---|---|---|---|---|---|---|\n' + '\n'.join(rows)
open(os.path.join(D, '_tables.md'), 'w', encoding='utf-8').write('## ALL\n' + ALL + '\n\n## DET\n' + DET + '\n')
# ORIGINAL_MAP.md
SXN = (11, 14, 24, 25, 26, 27, 28)
m = ['# ORIGINAL_MAP: DDBA_NN と元の台本の対応表', '',
     'DDBA_NN は、設計メモの版の通し番号 (1〜28)。設計メモの 1 ページが 1 つの版に対応する。この表は元の時系列の対応だけを残したもので、設計の流れは [README.md](README.md) (段階ごとの説明) に書いた。元の台本はファイル名だけで示す (置き場所の名前は書かない)。',
     '', '- 元の台本の確認の基準: 設計メモの添付 (名前と大きさ [KB]) と同名・同サイズ (1 KB = 1024 byte で小数 1 桁が合う) のもの。同名で複数あるときも、大きさが合うものを選んだ。',
     '- 台本ファイルの名前に氏名が入っているもの (DDBA_11、DDBA_27) は、元の台本名を「氏名を含むため省略」と書き、例のファイル名は氏名を除いた短い名前 (`v05_14_sx`、`v07_14_sx_20231205`) にした。台本の中の氏名のコメント 1 行 (DDBA_25・26) は英語の中立なコメントに置き換えた。それ以外は、ヘッダを足して改行を LF にしただけで、ラティスの定義は元のまま。',
     '- 設計メモのページに「添付された台本」がない (または新しい版がない) DDBA_NN は、そのページが使っている台本 (同じ設計メモの別のページの添付と同名・同サイズ) を割り当てた。',
     '- 六極の設定 (K2 の値) について: DDBA_11・14・24〜28 は、設計メモの記述によると、pysad の作者とは別の研究者 (公開論文 K. Harada et al., "Conceptual design of the Hybrid Ring with superconducting linac", J. Synchrotron Rad. 29, 118-124 (2022), DOI 10.1107/S1600577521012753 の著者の一人) の手による最適化の結果を含む。この論文には六極の K2 の値や最適化の手法の記述はなく (入射に必要な水平アパーチャを、通常の 16 回対称 DDBA リングの約 80% まで回復した、という記述のみ)、論文の暫定設計 (3 GeV、周長 350 m、16 セル) はこれらの版とは別の版である。論文の設定と一致する、同じ最適化の結果である、という意味ではない。他の版の六極は、設計メモの記述では設計者自身の探索 (S1 = -S2、S3 = S5 など) による。', '',
     '| DDBA_NN ← 元の台本 | 例のファイル | 設計メモの添付 [KB] | 元の台本 [byte] | 判定 | 設計メモの日付 | 設計メモの説明 (要約) | 六極の注記 |', '|---|---|---|---|---|---|---|---|']
for nn, src, out in L.MAP:
    inf = I.I[nn]; sz = inf['SZ']
    m.append('| DDBA_%02d ← %s | `%s` | %s | %d | %s | %s | %s | %s |' % (nn, '`%s`' % src if SRC(src) == src else SRC(src), fn(nn), sz[0], sz[1], sz[2], inf['DATE'], inf['SUM'], '別の研究者の手による最適化を含む' if nn in SXN else ''))
open(os.path.join(D, 'ORIGINAL_MAP.md'), 'w', encoding='utf-8', newline='\n').write('\n'.join(m) + '\n')
print('ok')
