# 作業ログ: 幾何の fit (GX/GY/GZ/CHI)・DISP G/OG・ビームサイズ (SIGX) (2026-11-05)

残課題 2 の「DISP G モードと fit の GX/GY」「BEAMSIZE/Twiss["SIGX"]」。
ユーザー指示 (11/05): ファイルの読み書きは Python ネイティブでやるので SAD の OpenRead/Write の失敗の扱いは移植しない。
DISP G と GX/GY のマッチングは重要。SIGX/SIGY などは **EMIT の事前実行を条件に、エネルギー広がり・エミッタンス・Twiss からビームサイズを返す**
(SAD はグローバル EMITX/EMITY/SIGE をそのまま使うが、ここでは仕様変更)。ただし基準座標の値と物理座標への射影は SAD を反映する。

## 1. 幾何の fit 条件 GX / GY / GZ / CHI1 / CHI2 / CHI3

- SAD: tgfun.f (`GX..GZ` = geo(:,4,kp)、`CHI1..3` = tfchi(geo(:,:,kp)))、重みは tffscalc.f (GX.. は 1e-3·sqrt(ν/ε/L)、CHI は 1e3)、
  微分は tdgeo.f (DRIFT の L、BEND の L/ANGLE、ソレノイド内の DX/DY の解析微分)。
- pysad: `match6.Matcher._tgfun` が幾何の条件のとき `GeoResult` (pysad/geo.py の tfgeo の直訳) から値を返す。幾何は変数を反映した要素から
  `_tffscalc` のたびに作り直す (`Matcher.geo_fn` = `FFS._matcher_geo`: ORG 込み)。微分は tdgeo の解析微分ではなく **数値微分** (`_analytic_unsupported` が
  幾何の条件で中心差分に切り替える)。重みと tdfun (CHI の ±π の折り返し) は以前から入っていた。
- 直した点: FFS で要素の `L` を変えても `el._raw['L']` が更新されず、幾何 (geo.py の `_val` は `_raw` を先に読む) が変わらないので
  `free D L` の幾何の微分が 0 になっていた → `_set_value` で L/ROTATE/ANGLE/DX/DY/DZ を `_raw` にも書く。
- 条件表 (`print_fit_table`) は幾何の値を `fit_fun_value` で出す (Twiss["GX"] は K64 と同じくエラーのまま)。CHI の表示は度。
- `?No variable.` (変数なしの GO) のときは SAD と同じく `Matched.` の行だけで条件表を出さない (tffscalc.f:283 の error=.true.)。

## 2. DISP G_EOMETRY / OG_EOMETRY

tfdisp.f:322-370 の直訳 (`FFS._disp_geo`): ` Element` + 名前 18 桁、Gx Gy Gz (12.6)、s (12.6、tfgeo の pos)、Length (10.6)、Value (10.7、MARK/MONI は ` 0`)、
Chi1-3 (度、12.6)、#。OG は閉軌道 (DX DPX DY DPY) を乗せた幾何 (tforbitgeo)。ヘッダは 66 行ごと、反転要素は行頭 `-`。

## 3. ビームサイズ Twiss["SIG.."] / LINE["SIG.."] / BEAMSIZE

- SAD: tftwiss.f:270 で `SIG`/`SIZE` は tfline へ → tftwiss.f:982 `tfbeamkey` (`SIG` の後ろ: 空 → 6×6、X/PX/Y/PY/Z/DP → sqrt(Σii)、XY など → Σij)
  → tfsize/tfbeamfrac → **tfbeam.f** (CALC6D でないとき): 基準座標 (ノーマルモード) の Twiss (AX BX EX EPX AY BY EY EPY) と emx/emy
  (MINCOUP の下限 emmin=(emx+emy)·MINCOUP)、pfi = SIGE² (gauss=.true.: EMIT 後) からビーム行列を作り、qtent (R1-R4 の 4×4) で物理座標へ
  (tmulbs: Σ' = TΣTᵀ)。小数位置は qtwissfrac の Twiss で同じ式。CALC6D のときは temit の平衡ビーム行列。
- pysad (`FFS.beam_matrix_twiss`, `beam_size_value`): tfbeam の直訳。**EMIT を先に実行していないと SADError** (SAD はグローバルをそのまま使う)。
  EMIT が決めた EMITX/EMITY/SIGE (INTRA ならその収束値) と現在の Twiss から作るので、`emit; ...; cell calc; Twiss["SIGX",..]` は
  新しい Twiss と EMIT のエミッタンスの組み合わせになる。gammab の比 (加速) は 1。CALC6D のときは従来どおり 6D の平衡ビーム行列 (bench_level21a)。
- `BEAMSIZE` コマンド (tffsa.f BEAM_SIZE: tfsetparam + tfsize) は EMIT 済みかを確かめるだけ (値は使うときに作る)。EMIT 前なら `?BEAMSIZE: EMIT を先に実行してください`。

## 4. 検証

| スクリプト | 結果 |
|---|---|
| esrf matching_geo_hmba (`disp geo $$$`、fit gx/gy、free DRID1 L) | SAD と同じ出力 (Matched. 6.2796E-12 まで同じ、GY は 1e-13 の雑音) |
| pf2 matching_geo_ring_r (fit MKRC GX/GY、free BMA → No variable) | SAD と同じ (`?No variable.` → Matched. の行だけ、GX/GY は 1e-14 の雑音) |
| <work-perf>/geo_tiny1 (小さな輸送線: 回転 BEND 付き、disp geo/ogeo) | DISP G の全行が SAD と同一 (OG の 1e-33 の雑音を除く)、LINE["GX".."GCHI3"] が 1e-15 で一致 |
| <work-perf>/geo_tiny2 (fit M2 GX/GZ で free D2 L, D3 L → fit $$$ CHI1/CHI3 で free B3/B2 ANGLE) | 反復の表示 (3〜12 回目) が SAD と同一、D2 L 2.784721292 (1e-12)、角度は 15 桁一致 |
| pf2 beamsize_r (emit; beamsize; cell calc beamsize; Twiss["SIGX".."SIGPY"] の 146 行) | **146 行すべて SAD と同一** (1e-13 以下の雑音の列を除く) |

テスト: `pysad/tests/test_geo_fit_sig.py` (DISP G の行と LINE["GX"..]、幾何の fit の SAD 値、Twiss["SIG.."] の式と EMIT 前のエラー、BEAMSIZE)。
pytest: 250 passed (34 s)。`QUICK=1 run_regression.py`: FAIL 0 (50 ベンチ、213 s)。

## 5. 変更したファイル

- `pysad/match6.py`: `_GEO_KF`、`Matcher.geo_fn/_geo/_compute_geo`、`_tgfun` (幾何)、`_tffscalc` で幾何を作り直す、`no_variable`。
- `pysad/ffscmd.py`: `_matcher_geo`、`fit_fun_value`、`print_fit_table` (幾何・CHI の度)、`_disp_geo` (DISP G/OG)、`beam_size_value`/`beam_matrix_twiss`、
  `twiss_value` の SIG..、`line_value` の SIG.. (CALC6D で分岐)、BEAM_SIZE コマンド、`_set_value` の `_raw` 更新、GO 後の条件表 (No variable)。
- `pysad/tests/test_geo_fit_sig.py` (新規)。
- `<work-perf>/`: patch5_geo.py、patch6_sig.py (差し替えの記録)、geo_tiny1/2.sad (+ _full.sad、.sad.out = SAD、.pysad.out)、geo_hmba.pysad.out、geo_ring.pysad.out、beamsize_r.pysad.out。
