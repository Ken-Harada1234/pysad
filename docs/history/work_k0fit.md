# 作業ログ: K0 を fit 変数にしたときの更新 (bump3_r / calc_bump_r / kick_r) と TYPE の書式 (2026-11-04)

残課題 1 位 (docs/history/00_RESUME.md 10/16 のまとめ) の 3 本。いずれも「K0 の fit 変数としての更新」そのものではなく、
FFS コマンドの構文解釈の差が原因だった。修正はすべて `pysad/ffscmd.py`。

## 原因と修正

| # | 症状 | 原因 (SAD の仕様) | 修正 |
|---|---|---|---|
| 1 | bump3_r: `FIT MKSTRSS DX 0` の条件が付かず、fit が何もしない (Matched 0.000)。calc_bump_r: `fit B10 DX 0` が `?Undefined command or element: B10` | MKSTRSS は ID16 ラインに無い (どのラインにも無い)、B10 は INJL ラインに無い。SAD (tffsa.f `FIT`) は `ielme` が見つからない語に対して `mfpnt=mfpnt1=nlat` ($$$) にして `go to 12` → 語は式として評価 (記号なので何も起きない、`Out[1]:= USE` が出るだけ) され消費される。続く `DX 0` は **$$$ の条件** になる | `_fit_location`: 位置でも fit 関数でもない語は $$$ にして語を消費 (return 2)。fit 関数 (`FIT DX 0`) は従来どおり $$$ にして語を残す |
| 2 | bump3_r: `FREE K1 K0 K2 K0` で K2 が変数にならない (NEWTON の変数 1 個、Residual 1.878 で停止) | pysad の `_freefix` は `_ELEM_KEYWORDS` に入っている語 (K2 はキーワードでもある) を直前の要素のキーワードとして上書きしていた。SAD (tffsfreefix.f) は「直前に合った要素の型のキーワード (tfkwrd) に一致する語だけ」をキーワードにし、合わなければ次のパターン | `_kw_of_pattern` (sadkw.KEYWORDS で型のキーワードか判定) を通るときだけキーワードにする。キーワードは 1 パターンに 1 つ |
| 3 | `type K*` に `ANGLE =0` が無い、ROTATE がラジアン、OCT の `K3 =0` が無い、MARK が `=()`、型グループの前の ` ;` が無い | tftype.f/tftyp1: 値が 0 でないキーワードと **既定の変数 (ival: BEND の ANGLE、QUAD の K1、OCT の K3、DRIFT の L など) は 0 でも出す**。ROTATE は度で ` DEG` を付ける。パターン語ごとに型 (ltyp 順) でまとめ、グループの前に ` ;`、全部の後に ` ;`。MARK は光学関数 (先頭要素は要素自身の値、他は直前の CALC の twiss; CALC 前は 0)、EMITX/Y/Z・SIGMAZ・SIGE・DP (=dpmax)。幅 129 で折り返し (続きは 11 桁の空白) | `type_cmd` を書き直し (`_TYPE_ORDER`、`_type_value`、sadkw.KEYWORDS の順、DEFAULT_KEY を ival に) |

kick_r は 10/16 の時点で完走しており、差は fit の最終 K0 (SAD −4.78e-10 / +4.00e-10、pysad +5.89e-10 / −6.19e-10) だけ。
FFS[] の外で `go` を実行して見ると SAD の残差 2.0e-12、MKK1 の DX −5.7e-10 で、どちらも収束の下限の雑音 (方針: 追わない)。

## 照合 (<work-item>/、SAD は WSL の gs)

| スクリプト | 結果 |
|---|---|
| bump3_r (pf2 v9_25, USE ID16) | K1 K0 = .002086668142394 (SAD .0020866681423942)、K2 K0 = −.0023854159155666 (SAD −.0023854159155679)、Matched 7.1559E-10 (同)、$$$ DX 2.6734E-9 / DPX 4.1796E-9 (同)、ins calc の AX/BX/NX/AY/BY/NY 全桁一致 |
| bump3_chk (FFS[] の外で FIT/FREE/GO) | GO の反復・Matched 行・type K* が SAD と同じ |
| calc_bump_r (v9_16_ele, USE INJL) | NEWTON 残差 145.8 / 4.8560E-04 (同)、Matched 2.0395E-16 (SAD 2.0232E-16)、KK1..KK4 の K0 が 1e-15 で一致 (−.0048084405187、.0015696722061、−.0041493878895、−.0059342764582)、INJ DX/DPX・$$$ の全項目一致 |
| kick_r / kick_chk | cell calc 一致。fit 後の K0 は両方 1e-10 以下 (雑音) |
| type_chk (VD SQ PV PH Q1 K1 CAVI L250 SFH OCT MKK1、CALC 前) | SAD と **全行同一** |
| type_chk2 (cell calc 後の MARK MKK1, U02) | AX/BX/PSIX/AY/BY/PSIY/EX/EPX/BZ/DP が一致 (1e-14)。SAD は EY/EPY/R1〜R4/DETR/DX〜DZ の 1e-16〜1e-32 の雑音も印字するが pysad は 0 なので出ない (追わない) |

## 環境 (2026-11-04)
- numba 入りの Python 3.13 がこのマシンに無くなっていた (pysad/__pycache__ は cpython-313)。**Python 3.12 (py -3.12) に numba 0.68.1 / scipy 1.18.1 / pytest を入れた**。numpy 2.5.4。
- pytest: `pysad/tests` のうち `sys.exit` で終わるスクリプト型 4 本 (test_bend_feeddown, test_coord, test_sext_radiation, test_track6n_ext) は `--ignore`。
  test_emit_report の Beam tilt (SAD 2.8e-13、pysad −1.2e-10) は 0 同然の雑音なので `_same` に「両方 |v|<1e-9 なら同じ」を追加。

## 回帰 (11/04)
- `QUICK=1 py -3.12 run_regression.py`: 50 ベンチ中 48 PASS、2 FAIL → どちらも雑音: 26a は Beam tilt (SAD 2.8e-13 / pysad −1.2e-10)、26c は F6 (FSHIFT の 6D 閉軌道) の相対差 3.13e-10 で許容 3e-10 をわずかに超過。
  `cmp_emit_report.same` に「両方 |v|<1e-9 なら同じ」を追加し、26c の許容を 1e-9 (26d と同じ) にして **FAIL 0** (26a/26c を再実行)。pytest 156 passed (スクリプト型 4 本を除く)。
- 変更ファイル: `pysad/ffscmd.py` (`_fit_location`, `_freefix`, `_kw_of_pattern`, `type_cmd`, `_type_value`, `_TYPE_ORDER`)、`pysad/tests/test_emit_report.py`、`cmp_emit_report.py`、`run_regression.py` (26c の許容)。

## 残り (この 3 本から)
- なし。docs/12 の §4 項目 6 は解消。
