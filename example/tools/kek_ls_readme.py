# -*- coding: utf-8 -*-
"""example/kek_ls の design・display・errors の README.md と、kek_ls/README.md の表を作る。
一致の判定 (○△×) と所要時間は、gs・pysad の出力と NAME.time.txt から自動で入れる。先に run_kekls.py と kek_ls_md.py を走らせる。
使い方 (リポジトリ直下から): py -3.12 example/tools/kek_ls_readme.py"""
import os, sys, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kek_ls_cmp as C
import kek_ls_md as M

EX = C.EX

SHORT = {
 'KEKLS_cell_calc_emit': ('RINGC・RING の周期解とエミッタンス', '-', 'lattice/kek_ls_v03_68.sad', '周長・チューン・εx・σδ'),
 'KEKLS_cell_optics': ('標準セルの配置・オプティクスの数表と図', 'disp_optics', '同上', '数表 160 行、PNG (`../tools/kek_ls_figs.py`)'),
 'KEKLS_match_cell_tune': ('セルのチューンのマッチング', 'matching2_01 前半', '同上', 'セルの位相進み・四極の K1'),
 'KEKLS_match_chrom': ('6 極のクロマティシティ調整', 'matching2_01 後半', '同上', 'ξ・6 極の K2'),
 'KEKLS_match_injcell': ('入射セルのマッチング', 'matching2_03_INJ', '同上', '入射点の Twiss・QS0〜QS5 の K1'),
 'KEKLS_match_ring_geo': ('入射セルを含めてリングを閉じる (ジオメトリ)', 'matching2_02_geo', '同上', 'DRID1INJ・GX・GY'),
 'KEKLS_match_repeat': ('マッチングの再現性の診断', '(新規)', '同上', 'K1・K2 (3 回の繰り返し)'),
 'KEKLS_scale_mini': ('セルの一様な縮小', 'matching23_v03_35_mini', '同上', 'セル長・位相進み・K1'),
 'KEKLS_scale_ncell': ('セル数を変える', 'matching21_v03_35', '同上', 'j セルのリングの εx・チューン・周長'),
 'KEKLS_tune_survey_cell': ('セルのチューンサーベイ', 'matching05 + modules_tune_sv_v03_68.n', '同上', 'ξ・εx・6 極の強さの 25 点の地図'),
 'KEKLS_ring_params': ('リングのパラメータ (εx・σδ・減衰・αc・U0・ξ)', 'calc_emit.sad 前半、disp、chromaticity_plot', '同上', 'emit の表・ξ'),
 'KEKLS_optics_table': ('リングのオプティクス・配置の数表と図', 'disp_optics、disp、draw_cad の入力', '同上', '要素ごとの表 3226 行、PNG'),
 'KEKLS_calc_H': ('H 関数の積分によるエミッタンス', 'calc_H.sad、calc_emit.sad 後半', '同上', 'εx (SAD と 2.3% 差)'),
 'KEKLS_chromaticity': ('チューンの運動量依存', 'chromaticity_plot.sad', '同上', 'ν(δ)・PNG'),
 'KEKLS_amp_dep_ts': ('振幅依存チューンシフト', 'amp_dep_ts*.sad', '同上 + modules/modules_amp_dep_ts.n', 'ν(振幅)・PNG'),
 'KEKLS_sextupole_tune_scan': ('6 極を変えたチューンシフト', 'dp_sv_sx.sad', '同上', 'Σ|Δν| の 5×5 の地図・PNG'),
 'KEKLS_scan_target': ('マッチングの目標値の振り (δ 依存・ADTS)', 'dp_sv.sad、dx_sv*.sad、dy_sv*.sad', '同上 + modules/modules_cell_match.n', 'ξ・ν(δ)・ADTS'),
 'KEKLS_errors_cod': ('誤差の種類ごとの COD・分散・チューン', 'print_all_error.sad ほか print_*', '同上 + modules/modules_errors.n', 'COD・分散・チューンのずれ'),
 'KEKLS_correct_cod': ('COD・分散の補正', 'correct_cod.sad', '同上 + modules/modules_general.n・modules_errors.n', '補正前後の COD・分散'),
 'KEKLS_cod_response': ('COD 応答行列と解析式の照合', 'write_cod_response.sad', '同上', '応答行列の代表値・解析式との差'),
}


def row(name, link=True):
    d = C.SPEC[name][0]
    r = C.compare(name)
    t = M.times(name, d)
    tm = ('gs %d s / pysad %d s' % (t['gs'], t['py'])) if t else '-'
    mark = r['mark']
    sh = SHORT[name]
    return d, '| [`%s`](%s.md) | %s | %s | %s | %s | %s |' % (name, name, sh[0], sh[1], sh[3], tm, mark)


HEAD = '| 例 | 内容 | 元のスクリプト | 出力 | 所要時間 | gs と pysad |\n|---|---|---|---|---|---|\n'

DESIGN = '''# design: KEK-LS の設計・マッチング

v03_68 のセル・リング・入射セルのマッチングと、セルの縮小・セル数の変更・チューンサーベイ。元のスクリプトを `kek_ls_v03_68.sad` で動くように整えた。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点 (`NAME.time.txt` は所要時間の記録)。

## 流れ

1. **セルのチューン** `KEKLS_match_cell_tune`: 標準セル `ARCA` の位相進みを (48.58, 17.62)/20 に合わせる (四極の K1)。
2. **6 極** `KEKLS_match_chrom`: リング `RINGC` の6 極 `SDH`・`SFH` でクロマティシティを 0 にする。
3. **入射セル** `KEKLS_match_injcell`: `ARCINJH` の四極 `QS0`〜`QS5` で、入射点の位相・ツイスを合わせる。
4. **リングを閉じる** `KEKLS_match_ring_geo`: 入射セルの長直線部 `DRID1INJ` で、ジオメトリを閉じる (GX = GY = 0)。
5. 設計の変形: `KEKLS_scale_mini` (全体の一様な縮小)、`KEKLS_scale_ncell` (セル数の変更)。
6. セルのチューンサーベイ `KEKLS_tune_survey_cell`、再現性の診断 `KEKLS_match_repeat`。

## 一覧 (gs と pysad の一致: ○ 一致、△ 差あり、× 不一致または未対応)

''' + '{ROWS}' + '''

マッチング後のラティスの保存: 各例は、結果を画面 (`type`・`Print`) に出す。元の台本の `SAVE ALL` はラティス全体を画面に書き出す (gs・pysad とも) ので、例では K1・K2 の値を `Print` する形にした。`SAVE ALL` の出力をファイルに取るには、台本の `SAVE ALL;` を FFS の最後に書いて標準出力を保存する (`FFS["SAVE ALL", ファイル]` の形はファイルが空になるので使えない)。

## 元の設定と例の設定

| 台本 | 元の設定 | 例の設定 | 完全な設定にするには | 見込みの所要時間 (完全) |
|---|---|---|---|---|
| `KEKLS_match_cell_tune` | `exmkc` 0.004 → 0 (0.001 刻み 5 点)、NPARA 32/16 | 0.004 → 0 (0.002 刻み 3 点)、NPARA 8 | `{exmkc,0.004,0.0,-0.001}` | gs 1 秒、pysad 数秒 |
| `KEKLS_match_injcell` | マッチングを 50 回繰り返す | 2 回 (1 回で収束) | `iterations = 50` | gs 数秒 |
| `KEKLS_scale_mini` | factor 1 → 0.5 を 0.01 刻み (51 点)、NPARA 8 | 0.1 刻み (6 点) | `{factor,1,0.5,-0.01}` | gs 約 10 秒、pysad 約 1 分 |
| `KEKLS_scale_ncell` | j = 20 → 2 (19 点)、NPARA 8 | j = 20, 18, 16, 12, 10 (5 点) | `names`・`ncs` を増やす | gs 数秒 |
| `KEKLS_tune_survey_cell` | tunex 46.05〜46.55 (6 点) × tuney 16.05〜19 (30 点) = 180 点、各点に DA (5000 周) | 48.30〜48.90 (5 点) × 17.30〜17.90 (5 点) = 25 点、DA なし | `tx0 = 46.05; dtx = 0.1; ntx = 6; ty0 = 16.05; dty = 0.1; nty = 30;` + `DAsv[]` | 25 点で gs 10 秒・pysad 約 160 秒から、DA なしで gs 約 1.2 分・pysad 約 19 分。DA は未測定 |

(見込みは、例の設定の所要時間からの比例の外挿。DA の時間は測っていない。元のスクリプトのファイルの更新日時は 2016 年で、当時の所要時間の記録は元のフォルダに無い。)

## 元の台本の問題と直したもの

- `matching2_02_geo`: `USE RINGI` (v03_68 に無い) → `RING`。
- `matching2_01`・`modules_tune_sv_v03_68.n`・`matching23`・`matching21`: `free BPI6 K1`・`free BPI7 K1` は v03_68 の名前 `BPI6H`・`BPI7H` に一致せず効かない。`matching23`・`matching21`・`tune_survey_cell` の例は `BPI6H`・`BPI7H` に直した (`freebend`)。`match_cell_tune` は元の挙動 (効かない) を既定に残した。
- `matching2_01`・`matching21`: FitFunction が `EMITX` を使うが、閉じていない ARCA では EMITX は NaN (`matching21` は `1^-9` の書き間違いも)。gs では何も効かず、pysad では fit が動かなくなる。この 2 本は FitFunction を除いた。
- `matching23`・`modules_tune_sv_v03_68.n`: v03_35 の `MKSFA2`・`SDA`・`SFA`・`QFA4A/B`・`ARCA1`・目標値が残っていた → v03_68 の名前・設計値に直した。
- `matching05`: `modules_tune_sv_v03_35.n` を読む。v03_68 用の `modules_tune_sv_v03_68.n` を読む台本は無い。
- 元の台本の `WriteGraphs[]` (`draw` コマンド)・`datas.txt` への書き出し・`SAVE ALL`・`stop; stop;` の後ろの旧い試行は、例では除いた (描画は PNG に別に、数表は画面に)。

## 要確認 (pysad)

- `FitFunction` を `fit` 条件と併用する `MatchingARCA` (`KEKLS_match_repeat`・`KEKLS_tune_survey_cell`・`KEKLS_scan_target`) で、pysad は gs と違う解に落ち、同じ入力の繰り返しでも結果が変わる (1 回目と 2 回目以降が違う)。NaN の FitFunction では fit が動かない。pysad のコードは変えていない (追跡に合わせる補正は足さない)。
- `GaussRandom[]` の `GCUT = 0` の扱い (gs は ≈ 0、pysad は打ち切りなし) は `../errors/` で確認 (台本は `GCUT = 1`)。
'''

DISPLAY = '''# display: KEK-LS のパラメータ表示

リング・オプティクス・チューンシフトの表示。元の `calc_emit.sad`・`disp*`・`chromaticity_plot.sad`・`calc_H*.sad`・`draw_cad.sad`・`amp_dep_ts*.sad`・`dp_sv*.sad`・`dx_sv*`・`dy_sv*`・`dp_sv_sx*.sad` を整えた。`Plot` 系 (Tk の窓) は使わず、台本の出力 (数表) を gs の出力から `tools/kek_ls_figs.py` が PNG にする (図は gs の出力から作った)。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点。

## 一覧 (gs と pysad の一致: ○ 一致、△ 差あり、× 不一致または未対応)

''' + '{ROWS}' + '''

図: [`KEKLS_ring_layout_optics.png`](KEKLS_ring_layout_optics.png) (リングの配置とオプティクス)、[`KEKLS_chromaticity.png`](KEKLS_chromaticity.png)、[`KEKLS_amp_dep_ts.png`](KEKLS_amp_dep_ts.png)、[`KEKLS_sextupole_tune_scan.png`](KEKLS_sextupole_tune_scan.png)。

![リングの配置とオプティクス](KEKLS_ring_layout_optics.png)

## 元の設定と例の設定

| 台本 | 元の設定 | 例の設定 | 完全な設定にするには | 見込みの所要時間 (完全) |
|---|---|---|---|---|
| `KEKLS_amp_dep_ts` | y の振幅 1〜401σ (20 刻み、21 点)、θ 10 点、NPARA 64 | x 1〜101σ (6 点)、y 1〜51σ (6 点) | `{nn,1,401,20}` (ただし 100σ 前後から先は NaN) | gs 数秒、pysad 約 2 分 |
| `KEKLS_chromaticity` | δ −3%〜+3% を 0.1% 刻み | 0.5% 刻み (13 点) | `dd = 0.001; nd = 60;` | gs 数秒、pysad 約 25 秒 |
| `KEKLS_sextupole_tune_scan` | dsf・dsd −3〜3 (0.2 刻み 31×31 = 961 点)、δ 13 点 | −3〜3 (1.5 刻み 5×5 = 25 点)、δ 7 点 | `ds0 = -3; dds = 0.2; nds = 31; dd = 0.005; nd = 12;` | gs 約 2.5 分、pysad 約 1 時間 (外挿) |
| `KEKLS_scan_target` | AY 1.08〜1.29 (0.01 刻み)、BY 3〜7 (0.1 刻み)、EX など。各点で ADTS 11 点 (x 1〜202σ、y 1〜401σ) | 3 つの目標 × −10%・0・+10% = 9 点、ADTS 61σ・31σ | `facs` を細かく | 1 点で gs 約 0.6 秒、pysad 約 13 秒。元の全部 (約 200 点) は pysad で約 45 分 (外挿) |

## 元のスクリプトのうち取り込んでいないもの

- `tune_survey_amp_dep_coef*`・`tune_survey_momentum`・`tune_survey_disp`・`tune_survey_linear*`: 別のリング用の `modules_tune_sv.n`・`modules_amp_dep_ts.n` の `MatchingWholeRing[]`・`ChromCorrect[]`・`LinearAmpDepTSCorrect[]` (v03_68 には無い) と、6 極 S1L・S2L… を使うため v03_68 には適用できない。ADTS・運動量依存の表示は、この例集の `KEKLS_amp_dep_ts`・`KEKLS_scan_target` で扱う。
- `draw_cad.sad`: 元の描画モジュール `drawmagvw3d.n` (元のフォルダに無い) が必要。
- `disp2`・`disp_arc_optics`・`disp_emit`・`disp_matching_tune`・`disp_matrix`: ID 入りラティス用の `modules_matching_noID.n`・`modules_write_correctdata.n` を使う (v03_68 に ID は無い)。
- `disp_b`・`disp_bunchlength`・`disp_eta`・`disp_NPARA_*`・`calc_H_2`・`dp_sv_R`・`dp_sv_2`・`dp_sv_v2_11`・`dx_sv_R`・`dx_sv_v2_11`・`dy_sv_R`・`dy_sv_v2_11`・`amp_dep_ts_{x,y,dp}`: 上の台本の版違い・断片 (別のラティスの名前や `stop;` で止まるもの)。
'''

ERRORS = '''# errors: KEK-LS の誤差入りの計算・補正

誤差 (位置ずれ・強さ・回転) を入れたリングの COD・分散・チューンのずれと、COD の補正。元の `correct_cod.sad`・`print_*.sad`・`write_cod_response.sad` を、短く動く形に整えた。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点。

## 誤差の与え方と補正の流れ

- 乱数の種: `SeedRandom[seed]`。例は seeds = {17, 111, 3} (補正の例は {17, 111})。元の台本は `SeedRandom[187]` と `Random[1000]` から作った 115 個の種。`GCUT = 1` (σ 単位の打ち切り。元の台本と同じ)。
- 位置ずれ ae = 40 μm: 四極・6 極の成分ごとの `LINE["DX"/"DY", 成分名] = ae GaussRandom[]` (`modules/modules_errors.n`)。元の `DELX`・`DELY` コマンドは、この K64 ビルドでは `cell calc` の閉軌道に効かない (DUMP には値が出るだけ。gs・pysad とも) ので使わない。
- 強さの誤差 fe = 4e-4 (`DK`)、回転の誤差 te = 100 μrad (`DTHETA`): FFS のコマンド (偏向磁石・四極・6 極に)。
- 成分ごとの値 (DX・DY・K0) は、`USE RING NOEXPAND` と `cell calc noexpand` でないと同名の要素の値に揃えられてしまう (`docs/11_sad_quirks.md` A56)。HD・VD (補正磁石) は同じ名前の繰り返しなので K0 は NOEXPAND が必要。
- 補正の流れ: (1) COD 応答行列を解析式で作る (`KEKLS_cod_response` で SAD の COD と照合)、(2) 誤差を与える、(3) 応答行列の特異値分解の疑似逆行列で K0 を求め、0.5 倍ずつ繰り返しかける、(4) COD・分散の rms と最大値、チューンの変化を表示する。

## 一覧 (gs と pysad の一致: ○ 一致、△ 差あり、× 不一致または未対応)

''' + '{ROWS}' + '''

## 元の設定と例の設定

| 台本 | 元の設定 | 例の設定 | 完全な設定にするには | 見込みの所要時間 (完全) |
|---|---|---|---|---|
| `KEKLS_errors_cod` | 種 115 個 | 種 3 個 × 誤差 4 通り | `seeds` を 115 個に、`nseed = 115` | 3 つの種 (各 4 通り) で gs 1 秒、pysad 6 秒 → 115 個で gs 約 40 秒、pysad 約 4 分 |
| `KEKLS_correct_cod` | 種 115 個 (k = 1〜115)、補正 40 回、毎回 DA (5000 周) | 種 2 個、補正 6 回、DA なし | `seeds`・`niter = 40`・`doDA = 1; DAturns = 5000` | 種 2 個で gs 1 秒、pysad 37 秒 (応答行列の準備に約 17 秒) → 115 個で gs 約 1 分、pysad 約 20 分。DA は未測定 |
| `KEKLS_cod_response` | 全ステアリング分を fort.92 に書く | 代表値と照合 3 個 | `j0` を 1〜ICX | pysad 約 9 秒 |

## 元のスクリプトのうち取り込んでいないもの

- `correct_cod_best`・`_worst`・`_noSX`・`_GCUT1`・`_v03_70`: `correct_cod.sad` の設定違い (monitor 名 `BPM.*`、DA の範囲、誤差の大きさ ae = 5e-5・fe = 5e-4、6 極を切る `S* 0`、チューンの窓 18〜18.5 など)。例の変数 (`ae`・`fe`・`te`・`doDA`、チューンの窓) を変えれば同じ形になる。
- `print_*.sad` 11 本: 同じ台本に表示だけ違う (`print_IEX_after`・`print_IEY_after`・`print_rep_after` は補正の固有値の数 IEX・IEY、`print_residual_cod`・`print_cod_kick_ratio` は補正後の COD・キック、`print_emittance`・`print_solution_ok`・`print_all_error_after`・`print_cod_form` は表示の形式)。誤差だけの `print_all_error`・`print_dk_effect`・`print_dt_error`・`print_cod_magnitude` の流れは `KEKLS_errors_cod`、補正は `KEKLS_correct_cod` に統合した。
- `modules_correctoptics.n`・`modules_write_correctdata.n`: `correct_cod.sad` は読み込んでいない (コメントアウトされた `Get` のみ)。光学補正 (ベータ関数の補正) の関数集で、v03_68 の補正には使っていない。`modules_correctoptics_pf.n` は PF 用。
'''

TOP_TABLE_HEAD = '| 例 | 元のスクリプト | 何をするか | 入力 | 出力 | 所要時間 (gs / pysad) | gs と pysad |\n|---|---|---|---|---|---|---|\n'


def top_rows():
    out = []
    for n in M.ORDER:
        d = C.SPEC[n][0]
        r = C.compare(n)
        t = M.times(n, d)
        tm = ('%d s / %d s' % (t['gs'], t['py'])) if t else '-'
        sh = SHORT[n]
        out.append('| [`%s/%s`](%s/%s.md) | %s | %s | %s | %s | %s | %s |' % (d, n, d, n, sh[1], sh[0], sh[2], sh[3], tm, r['mark']))
    return '\n'.join(out)


if __name__ == '__main__':
    for d, text in (('design', DESIGN), ('display', DISPLAY), ('errors', ERRORS)):
        rows = '\n'.join(row(n)[1] for n in M.ORDER if C.SPEC[n][0] == d)
        open(os.path.join(EX, d, 'README.md'), 'w', encoding='utf-8', newline='\n').write(text.replace('{ROWS}', HEAD + rows) + '')
        print(d, 'README')
    top = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kek_ls_src', 'README_top.md'), encoding='utf-8').read()
    open(os.path.join(EX, 'README.md'), 'w', encoding='utf-8', newline=chr(10)).write(top.replace('<!-- TABLE -->', TOP_TABLE_HEAD + top_rows()))
    print('README')
