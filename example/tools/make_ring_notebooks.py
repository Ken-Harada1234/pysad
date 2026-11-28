# -*- coding: utf-8 -*-
"""PF・ASTRID2・KEK-LS のラティスに紐付いた例 (example/pf・astrid・kek_ls の .sad) を、notebook (.ipynb) にする。
冒頭でラティスを `%%sad lattice` に貼り、そのあと各例のコマンドを、最初は `%%sad`、2 番目以降は `%%sad --continue` で順に流す。
コマンドの文は元の .sad の `FFS;` から `Exit;` までを機械的に取り出す (`$FORM`・`PageWidth`・`SAVE ALL`・`Exit` は除く)ので、元の例とずれない。

使い方 (notebook の環境は pysad_nb312 の venv。リポジトリの根を PYTHONPATH に置く):
    py -3.12 example/tools/make_ring_notebooks.py            # 3 つの .ipynb を (出力なしで) 書く
    py -3.12 example/tools/make_ring_notebooks.py --execute  # さらに本物のカーネルで実行して、出力つきで保存する
nbformat・nbclient・ipykernel が要る (`C:/Users/coding/venvs/pysad_nb312`)。
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.dirname(HERE)                      # example/
ROOT = os.path.dirname(EX)                      # リポジトリの根

DROP = re.compile(r'^\s*(\$FORM\s*=|PageWidth\s*=|Exit\s*;|SAVE\s+ALL\s*;)', re.I)


def read(path):
    return open(os.path.join(EX, path), encoding='utf-8', errors='replace').read().replace('\r\n', '\n')


def lattice_text(path):
    """ラティスファイルの中身 (FFS より前)"""
    t = read(path)
    m = re.search(r'^\s*FFS\s*;', t, re.M)
    return (t[:m.start()] if m else t).rstrip() + '\n'


def commands(path, start=None, fix=None):
    """スクリプトのコマンド部分: `FFS;` の行から `Exit;` まで (除く行を除いて)。start: この正規表現に最初に合う行から (その行を含む)。fix: [(旧, 新)]"""
    t = read(path)
    m = re.search(r'^\s*FFS\s*;', t, re.M)
    body = t[m.start():].split('\n')
    out = [ln for ln in body if not DROP.match(ln)]
    if start:
        k = next(i for i, ln in enumerate(out) if re.search(start, ln))
        out = out[k:]
    s = '\n'.join(out).strip('\n') + '\n'
    for a, b in (fix or []):
        assert a in s, a
        s = s.replace(a, b)
    return s


PLOT = '''import matplotlib.pyplot as plt
s  = sad.value('Twiss["S","*"]')
bx = sad.value('Twiss["BX","*"]'); by = sad.value('Twiss["BY","*"]'); ex = sad.value('Twiss["EX","*"]')
fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(s, bx, label='betax [m]'); ax.plot(s, by, label='betay [m]'); ax.plot(s, [100 * e for e in ex], label='etax x 100 [m]')
ax.set_xlabel('s [m]'); ax.legend(); fig.tight_layout(); plt.show()
print('要素の数', len(s), ' 周長 [m]', s[-1])'''


def ring_spec():
    pf = dict(
        out='pf/PF_in_notebook.ipynb',
        title='PF リング (2.5 GeV) を notebook で',
        intro=('PF リングの超低エミッタンスラティス (`example/pf/lattice/v9_16_asym.sad`。2.5 GeV、周長 187.035 m) を読み込み、'
               'example/pf の例のコマンドを順に流します。**最初にラティスを貼り、`%%sad` で最初の計算をして、`%%sad --continue` で続きを流します** '
               '(状態は前のセルから引き継がれます)。各例の説明と gs との比較は、それぞれの md にあります。'),
        lattice=('pf/lattice/v9_16_asym.sad', 'PF リング v9_16_asym (2.5 GeV)'),
        flags='ON RAD FLUC COD RFSW RADCOD;\n',
        steps=[
            ('1. セルの計算とエミッタンス (`optics/FFS_cell_calc_emit_PF.sad`)',
             '周長・チューン・エミッタンス・エネルギー拡がり・バンチ長を印字します。', commands('pf/optics/FFS_cell_calc_emit_PF.sad'), 'first'),
            ('2. β 関数と分散関数を描く', '`sad.value()` で全要素の値を取り出して、Python で描きます。', PLOT, 'python'),
            ('3. 放射積分 (`tune_radiation/calc_RadInt_PF.sad`)', '偏向磁石の中を 17 点で積分して、`cell calc emit` の値と比べます。',
             commands('pf/tune_radiation/calc_RadInt_PF.sad', start=r'^\s*bposi\s*='), 'continue'),
            ('4. 運動量依存チューンシフト (`tune_radiation/FFS_mom_dep_ts_PF.sad`)',
             '分散のずれ δ を ±2% まで振って、6 極あり・なし (`S* 0`) のチューンを求めます。**6 極を切るので、最後に流します** (切ったままになります)。',
             commands('pf/tune_radiation/FFS_mom_dep_ts_PF.sad', start=r'^\s*DP\s*=\s*0\s*;'), 'continue'),
        ],
    )
    return [pf]


def build_specs():
    specs = ring_spec()
    # 1) PF: 分散の確認 (RF を切るので、ON RFSW で戻す) と振幅依存チューンシフトを、運動量依存 (6 極を切る) より前に入れる
    pf = specs[0]
    pf['steps'].insert(2, ('3a. 分散の確認 (`optics/calc_dispersion_PF.sad`)',
                           '分散が閉軌道の δ 微分と一致するか、運動量圧縮係数が曲げ磁石の和と一致するかを確かめます。**RF を切って計算するので、最後に `ON RFSW` で戻します**。',
                           commands('pf/optics/calc_dispersion_PF.sad', start=r'^\s*pmex\s*=') + 'FFS["ON RFSW"];\nFFS["cell calc emit"];\n', 'continue'))
    pf['steps'].insert(4, ('3c. 振幅依存チューンシフト (`tune_radiation/FFS_amp_dep_ts_PF.sad`)',
                           'σ の何倍かの初期振幅から、`ins calc` で 1 周のチューンを求めて、楕円上で平均します (追跡ではなく線形でない写像の計算)。',
                           commands('pf/tune_radiation/FFS_amp_dep_ts_PF.sad', start=r'^\s*SetInitialParameters'), 'continue'))
    astrid = dict(
        out='astrid/ASTRID2_in_notebook.ipynb',
        title='ASTRID2 (0.58 GeV) を notebook で',
        intro=('ASTRID2 のラティス (`example/astrid/lattice/ASTRID2_v01_01.sad`。0.58 GeV。30° の複合機能偏向磁石 (K1≠0) を持つ) を読み込み、'
               'example/astrid の例のコマンドを順に流します。最後に、pysad 独自の解析式 (`AmplitudeTuneShift[]`・`ResonanceDrivingTerms[]`。gs には無い) も流します '
               '(解析の検証は docs/25 §11)。**最初にラティスを貼り、`%%sad` で最初の計算をして、`%%sad --continue` で続きを流します**。'),
        lattice=('astrid/lattice/ASTRID2_v01_01.sad', 'ASTRID2 v01_01 (0.58 GeV)'),
        flags='ON RAD COD RFSW RADCOD;\nOFF FLUC;\n',
        steps=[
            ('1. セルの計算とエミッタンス',
             '`USE RING; cell calc emit;` と主なパラメータの印字です (`tune_radiation/calc_RadInt_ASTRID2.sad` の冒頭と同じ計算)。',
             'FFS;\nUSE RING;\ncell calc emit;\nPrint["circumference [m] = ", LINE["S","$$$"]];\nPrint["tune nx ny = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi];\n'
             'Print["EMITX [m rad] = ", EMITX, "  SIGE = ", SIGE];\nPrint["MOMENTUM [eV] = ", MOMENTUM];\n', 'first'),
            ('2. β 関数と分散関数を描く', '`sad.value()` で全要素の値を取り出して、Python で描きます。', PLOT, 'python'),
            ('3. 放射積分 (`tune_radiation/calc_RadInt_ASTRID2.sad`)', '複合機能偏向磁石 (K1≠0) の中の積分 (I4 の勾配の項、H の積分) を、解析式で求めて `cell calc emit` と比べます。',
             commands('astrid/tune_radiation/calc_RadInt_ASTRID2.sad', start=r'^\s*bposi\s*='), 'continue'),
            ('4. 振幅依存チューンシフト (`tune_radiation/amp_dep_ts_ASTRID2.sad`)', '`ins calc` で 1 周のチューンを求めて、楕円上で平均します。',
             commands('astrid/tune_radiation/amp_dep_ts_ASTRID2.sad', start=r'^\s*AmplitudeDependentTuneShift'), 'continue'),
            ('5. pysad 独自の解析式: ADTS と共鳴駆動項 (gs には無い)',
             '`AmplitudeTuneShift[]` (解析) と `ResonanceDrivingTerms[]` を、ラティスの光学から求めます。偏向磁石の 3 次の項 (複合機能の本体・端) は、追跡の写像から取り出して入れてあります (docs/25 §8)。',
             'USE RING; cell calc emit; DP0 = 0; DP = 0;\nPrint["ADTS analytic ", AmplitudeTuneShift[]];\nPrint["RDT ", ResonanceDrivingTerms[Terms->{"h21000","h30000","h10110","h10020","h10200"}]];\n', 'continue'),
            ('6. 運動量依存チューンシフト (`tune_radiation/mom_dep_ts_ASTRID2.sad`)',
             'δ を −2〜+4% で振って、6 極あり・なし (`S* 0`) のチューンを求め、多項式でフィットします。**6 極を切るので、最後に流します**。',
             commands('astrid/tune_radiation/mom_dep_ts_ASTRID2.sad', start=r'^\s*DP\s*=\s*0\s*;'), 'continue'),
        ],
    )
    kek = dict(
        out='kek_ls/KEKLS_in_notebook.ipynb',
        title='KEK-LS (3 GeV、20 セル、ESRF 型 HMBA) を notebook で',
        intro=('KEK-LS v03_68 (`example/kek_ls/lattice/kek_ls_v03_68.sad`。3 GeV、周長 570.72 m、ESRF 型の HMBA。ESRF の Pantaleo Raimondi 氏と Simone Liuzzo 氏の協力で設計した) を読み込み、'
               'example/kek_ls の例のコマンドを順に流します。**最初にラティスを貼り、`%%sad` で最初の計算をして、`%%sad --continue` で続きを流します**。'
               'この notebook は `example/kek_ls/` で実行してください (振幅依存チューンシフトが `modules/modules_amp_dep_ts.n` を相対パスで読みます)。'),
        lattice=('kek_ls/lattice/kek_ls_v03_68.sad', 'KEK-LS v03_68 (3 GeV)'),
        flags='ON RAD FLUC COD RFSW RADCOD;\n',
        steps=[
            ('1. リングのパラメータ (`display/KEKLS_ring_params.sad`)',
             'RINGC (標準セル 20 個) と RING (入射セルつき) の周長・チューン・エミッタンス・自然クロマティシティとクロマティシティを印字します (6 極は一時的に切って戻します)。',
             commands('kek_ls/display/KEKLS_ring_params.sad'), 'first'),
            ('2. β 関数と分散関数を描く', 'RING の結果を `sad.value()` で取り出して、Python で描きます。', PLOT, 'python'),
            ('3. クロマティシティ (`display/KEKLS_chromaticity.sad`)', 'δ を ±3% まで振ったチューンを、6 極あり・なしで求めます (6 極の K2 は戻します)。',
             commands('kek_ls/display/KEKLS_chromaticity.sad', start=r'^\s*USE\s+RINGC'), 'continue'),
            ('4. 振幅依存チューンシフト (`display/KEKLS_amp_dep_ts.sad`)', 'σ の何倍かの初期振幅から、`ins calc` でチューンを求めて、楕円上で平均します。',
             commands('kek_ls/display/KEKLS_amp_dep_ts.sad', fix=[('"../modules/modules_amp_dep_ts.n"', '"modules/modules_amp_dep_ts.n"')]).replace('FFS;\n', '', 1), 'continue'),
        ],
    )
    return [specs[0], astrid, kek]


def make(spec, nbf):
    nb = nbf.v4.new_notebook()
    nb.metadata['kernelspec'] = {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'}
    c = nb.cells
    c.append(nbf.v4.new_markdown_cell('# ' + spec['title'] + '\n\n' + spec['intro'] +
                                      '\n\n出力は pysad (既定の `v2`) の実行結果です。K64 (gs) と同じ動作にするときは、最初に `PYSAD_COMPAT=k64` を環境変数に置きます。'
                                      '各例の md には、gs と pysad の出力の比較があります。このファイルは `example/tools/make_ring_notebooks.py` で作りました (元の .sad から機械的に取り出しています)。'))
    c.append(nbf.v4.new_code_cell('%load_ext pysad.notebook'))
    c.append(nbf.v4.new_markdown_cell('## 0. ラティスを読み込む: ' + spec['lattice'][1] + '\n\n`%%sad lattice` は、ラティスを保存します (まだ実行しません)。'))
    c.append(nbf.v4.new_code_cell('%%sad lattice\n' + lattice_text(spec['lattice'][0])))
    for i, (title, text, code, kind) in enumerate(spec['steps'], 1):
        title = '%d. ' % i + re.sub(r'^\d+[a-z]?\.\s*', '', title)
        if kind == 'continue' and code.startswith('FFS;\n'):
            code = code[len('FFS;\n'):]
        c.append(nbf.v4.new_markdown_cell('## ' + title + '\n\n' + text))
        if kind == 'python':
            c.append(nbf.v4.new_code_cell(code))
        elif kind == 'first':
            c.append(nbf.v4.new_code_cell('%%sad\n! 元の例のスクリプトと同じフラグ (FFS の前に置く)\n' + spec['flags'] + code))
        else:
            c.append(nbf.v4.new_code_cell('%%sad --continue\n' + code))
    return nb


def execute(path, nbf, cwd):
    """本物のカーネルで実行して、出力つきで保存する。cwd: カーネルの作業フォルダ (KEK-LS は modules/ を相対パスで読むので example/kek_ls/)"""
    from nbclient import NotebookClient
    nb = nbf.read(path, as_version=4)
    NotebookClient(nb, timeout=1800, kernel_name='python3', resources={'metadata': {'path': cwd}}).execute()
    nbf.write(nb, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--only', default='')
    ap.add_argument('--outdir', default='', help='書き出し先 (検証用。省略すると example/ の中の元の場所)')
    a = ap.parse_args()
    import nbformat as nbf
    for spec in build_specs():
        if a.only and a.only not in spec['out']:
            continue
        path = os.path.join(a.outdir or EX, os.path.basename(spec['out']) if a.outdir else spec['out'])
        os.makedirs(os.path.dirname(path), exist_ok=True)
        nbf.write(make(spec, nbf), path)
        print('wrote', spec['out'], flush=True)
        if a.execute:
            execute(path, nbf, os.path.join(EX, os.path.dirname(spec['out'])))
            print('executed', spec['out'], flush=True)


if __name__ == '__main__':
    main()
