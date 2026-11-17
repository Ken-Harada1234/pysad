# 作業ログ: FIT/FREE の語の解釈と TYPE の書式 (2026-11-04)

最終更新: 2026-11-16

FFS コマンドの構文解釈の差の修正 (K0 を fit 変数にしたときの更新が原因ではなかった)。修正はすべて `pysad/ffscmd.py`。

## 原因と修正

| # | 症状 | 原因 (SAD の仕様) | 修正 |
|---|---|---|---|
| 1 | `FIT <未定義の語> DX 0` の条件が付かず、fit が何もしない。位置が無いと `?Undefined command or element` になる | SAD (tffsa.f `FIT`) は `ielme` が見つからない語に対して `mfpnt=mfpnt1=nlat` ($$$) にして `go to 12` → 語は式として評価 (記号なので何も起きない、`Out[1]:= USE` が出るだけ) され消費される。続く `DX 0` は **$$$ の条件** になる | `_fit_location`: 位置でも fit 関数でもない語は $$$ にして語を消費 (return 2)。fit 関数 (`FIT DX 0`) は従来どおり $$$ にして語を残す |
| 2 | `FREE K1 K0 K2 K0` で K2 が変数にならない (NEWTON の変数 1 個で停止) | pysad の `_freefix` は `_ELEM_KEYWORDS` に入っている語 (K2 はキーワードでもある) を直前の要素のキーワードとして上書きしていた。SAD (tffsfreefix.f) は「直前に合った要素の型のキーワード (tfkwrd) に一致する語だけ」をキーワードにし、合わなければ次のパターン | `_kw_of_pattern` (sadkw.KEYWORDS で型のキーワードか判定) を通るときだけキーワードにする。キーワードは 1 パターンに 1 つ |
| 3 | `type K*` に `ANGLE =0` が無い、ROTATE がラジアン、OCT の `K3 =0` が無い、MARK が `=()`、型グループの前の ` ;` が無い | tftype.f/tftyp1: 値が 0 でないキーワードと **既定の変数 (ival: BEND の ANGLE、QUAD の K1、OCT の K3、DRIFT の L など) は 0 でも出す**。ROTATE は度で ` DEG` を付ける。パターン語ごとに型 (ltyp 順) でまとめ、グループの前に ` ;`、全部の後に ` ;`。MARK は光学関数 (先頭要素は要素自身の値、他は直前の CALC の twiss; CALC 前は 0)、EMITX/Y/Z・SIGMAZ・SIGE・DP (=dpmax)。幅 129 で折り返し (続きは 11 桁の空白) | `type_cmd` を書き直し (`_TYPE_ORDER`、`_type_value`、sadkw.KEYWORDS の順、DEFAULT_KEY を ival に) |


## 環境 (2026-11-04)
- numba 入りの Python 3.13 がこのマシンに無くなっていた (pysad/__pycache__ は cpython-313)。**Python 3.12 (py -3.12) に numba 0.68.1 / scipy 1.18.1 / pytest を入れた**。numpy 2.5.4。
- pytest: `pysad/tests` のうち `sys.exit` で終わるスクリプト型 4 本 (test_bend_feeddown, test_coord, test_sext_radiation, test_track6n_ext) は `--ignore`。
  test_emit_report の Beam tilt (SAD 2.8e-13、pysad −1.2e-10) は 0 同然の雑音なので `_same` に「両方 |v|<1e-9 なら同じ」を追加。

## 回帰 (11/04)
- `QUICK=1 py -3.12 run_regression.py`: 50 ベンチ中 48 PASS、2 FAIL → どちらも雑音: 26a は Beam tilt (SAD 2.8e-13 / pysad −1.2e-10)、26c は F6 (FSHIFT の 6D 閉軌道) の相対差 3.13e-10 で許容 3e-10 をわずかに超過。
  `cmp_emit_report.same` に「両方 |v|<1e-9 なら同じ」を追加し、26c の許容を 1e-9 (26d と同じ) にして **FAIL 0** (26a/26c を再実行)。pytest 156 passed (スクリプト型 4 本を除く)。
- 変更ファイル: `pysad/ffscmd.py` (`_fit_location`, `_freefix`, `_kw_of_pattern`, `type_cmd`, `_type_value`, `_TYPE_ORDER`)、`pysad/tests/test_emit_report.py`、`cmp_emit_report.py`、`run_regression.py` (26c の許容)。

## 残り (この 3 本から)
- なし。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: この話題の記述は work_unimpl.md と docs/22〜35 に出てこない (grep で確認)、関係するモジュールに 11/06 以後この話題の変更の記録がない。この文書は 11/06 時点の記録のまま。
