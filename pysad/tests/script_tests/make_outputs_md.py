# -*- coding: utf-8 -*-
"""pysad/tests/script_tests/OUTPUTS.md を作る。scripts/ の NAME.sad・NAME.gs.txt・NAME.pysad.txt と試験関数の対応表。

使い方: py -3.12 pysad/tests/script_tests/make_outputs_md.py
NAME は `<試験ファイル名 (test_ を除く)>__<試験関数 (クラス名を含む)>[_2, _3, ...]`。
scripts/ のファイルは、試験を実行して run_sad_string に渡った台本と画面出力を取り出したもの
(同じ試験の中の台本が複数あるときは _2, _3, ...)。gs 側は同じ台本を WSL の gs で走らせた画面表示。
"""
import ast
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SC = os.path.join(HERE, 'scripts')


def first_doc(path):
    try:
        t = ast.parse(open(path, encoding='utf-8', errors='replace').read())
        d = ast.get_docstring(t) or ''
        return d.strip().splitlines()[0] if d.strip() else ''
    except Exception:
        return ''


def gs_state(p):
    if not os.path.exists(p):
        return 'なし'
    s = open(p, 'rb').read().decode('utf-8', 'replace')
    n = len([l for l in s.splitlines() if l.strip()])
    if 'Segmentation' in s or 'SIGSEGV' in s:
        return 'SIGSEGV (%d 行)' % n
    body = [l for l in s.splitlines() if l.strip() and not l.startswith(' ***')]
    if len(body) <= 3:
        return '出力なし (止まる・落ちる)'
    return '%d 行' % n


def main():
    sads = sorted(glob.glob(os.path.join(SC, '*.sad')))
    by_file = {}
    for p in sads:
        name = os.path.basename(p)[:-4]
        if '__' not in name:   # 試験関数に対応しない台本 (v241_* の長時間追跡・偏光 FLUC、v25_draw_gs_*) は表に入れない
            continue
        f, fn = name.split('__', 1)
        by_file.setdefault(f, []).append((name, fn))
    tests = sorted(glob.glob(os.path.join(HERE, 'test_*.py')))
    L = []
    A = L.append
    A('# pysad の試験が使う SAD スクリプトと、gs・pysad の画面表示')
    A('')
    A('最終更新: 2026-11-17')
    A('')
    A('pysad の試験 (`pysad/tests/script_tests/test_*.py`) が `run_sad_string` に渡す SAD スクリプト (試験の中の文字列の台本) を、そのまま '
      '`scripts/NAME.sad` に取り出し、同じ台本を pysad と gs で走らせた画面表示を `NAME.pysad.txt`・`NAME.gs.txt` として並べたものです。'
      'ファイル名は `<試験ファイル名 (test_ を除く)>__<試験関数 (クラス名を含む)>` で、同じ試験に台本が複数あるときは `_2`、`_3`…が付きます。'
      '表は `py -3.12 pysad/tests/script_tests/make_outputs_md.py` で作ります。')
    A('')
    A('## 読むときの注意')
    A('')
    A('- 試験は `scripts/` を読まない (台本は試験の文字列のまま)。`scripts/` は試験の実行から取り出した写しで、試験の判定は変えていない。台本を直したら取り出し直す。')
    A('- `.pysad.txt` は試験の実行中に台本が走ったときの画面出力 (標準出力 + 標準エラー)。`PYSAD_SAD_PACKAGES` を設定して全試験を通したときのもの。'
      '試験が途中でパラメータを変えて走らせるもの (`fallback` など) は、別の台本として `_2` 以降に入る。')
    A('- `.gs.txt` は同じ台本を WSL の `gs` (K64、2026-10-19 ビルド) にそのまま渡した画面表示 (`gs NAME.sad < /dev/null`、40 秒で打ち切り)。'
      '試験の台本は pysad 専用の整形 (`FFS;` の有無、`Exit;` ほか) のため、gs では走らないもの・出力が少ないものがある。gs の列にその状態を書いた。'
      '**gs の表示と pysad の表示の一致は、ここでは判定していない** (試験の判定は、試験に埋め込んだ SAD の参照値で行う)。')
    A('- 画像 (PNG): pysad の試験に画像を書き出すものは無い。描画・GUI の図は `example/11_graphics`・`example/12_gui_tkinter` ほかの PNG にある (`example/README.md`)。')
    A('- gs の `.gs.txt` に日付・時刻 (`Today:`)・一時ファイル名が入る。実行ごとに変わる。')
    A('')
    A('## 回帰ベンチ (`pysad/tests/regression_tests/`。コピーは作らない)')
    A('')
    A('`pysad/tests/regression_tests/run_regression.py` (50 ベンチ) の台本と SAD の参照結果は、`pysad/tests/regression_tests/bench/` に既にある (自動判定の対象外のものは `bench_extra/`)。')
    A('')
    A('| 何 | 場所 |')
    A('|---|---|')
    A('| SAD (gs) 用の台本 | `bench/bench_levelN.sad` (N = 1〜31 ほか) |')
    A('| SAD (gs) の画面表示 (参照結果) | `bench/bench_levelN.result` |')
    A('| pysad 用の台本 | `bench/bench_levelN_pysad.py` (一部は `.sad` を `run_sad.py` で走らせる) |')
    A('| pysad の画面表示 | `bench_extra/bench_levelN_pysad.result` (参考)、実行ごとの出力は `regression_tests/regression_out/<名前>.txt` |')
    A('| 判定と許容 | `regression_tests/run_regression.py`、`docs/27_benchmarks.md` |')
    A('')
    A('## 試験関数とスクリプト')
    A('')
    A('| 試験ファイル | 試験関数 | .sad | .gs.txt | .pysad.txt | gs の画面 |')
    A('|---|---|---|---|---|---|')
    nfile = 0
    for t in tests:
        tn = os.path.basename(t)
        key = tn[5:-3]
        if key not in by_file:
            continue
        nfile += 1
        for name, fn in by_file[key]:
            base = 'scripts/' + name
            A('| [`%s`](%s) | `%s` | [%s.sad](%s.sad) | [gs](%s.gs.txt) | [pysad](%s.pysad.txt) | %s |' % (
                tn, tn, fn, name, base, base, base, gs_state(os.path.join(SC, name + '.gs.txt'))))
    A('')
    A('## 対象外の試験 (SAD スクリプトを走らせない、または gs の参照が無い単体試験)')
    A('')
    A('台本を `run_sad_string` に渡さない試験ファイル。Python の関数・クラスを直接呼ぶ単体試験で、gs で走らせる台本がない。'
      '(`run_sad` や `load_sad_string`、`script_tests/*.sad` を使うものは下の表。)')
    A('')
    A('| 試験ファイル | 内容 (先頭の説明) |')
    A('|---|---|')
    for t in tests:
        tn = os.path.basename(t)
        if tn[5:-3] in by_file:
            continue
        A('| [`%s`](%s) | %s |' % (tn, tn, first_doc(t).replace('|', '/')))
    A('')
    A('## `pysad/tests/script_tests/` (と共有の `pysad/tests/`) に置いてある台本・ラティスの入力')
    A('')
    A('| ファイル | 使う試験 |')
    A('|---|---|')
    for f in sorted(glob.glob(os.path.join(HERE, '*.sad')) + glob.glob(os.path.join(HERE, '*.lat')) + glob.glob(os.path.join(HERE, '..', '*.lat'))):
        b = os.path.basename(f)
        users = sorted(os.path.basename(t) for t in tests if b in open(t, encoding='utf-8', errors='replace').read())
        A('| `%s` | %s |' % (b, ', '.join('`%s`' % u for u in users) or '(直接は参照されない)'))
    A('')
    A('`scripts/` の台本数: %d (試験ファイル %d 本)。' % (len(sads), nfile))
    A('')
    open(os.path.join(HERE, 'OUTPUTS.md'), 'w', encoding='utf-8', newline='\n').write('\n'.join(L))
    print(len(sads), nfile)


if __name__ == '__main__':
    main()
