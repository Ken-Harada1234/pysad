# 12. ユーザー実スクリプト (esrf / pf2) の pysad 対応状況

最終更新: 2026-11-06 (現状のスナップショット)

目標: `ユーザーの esrf ディレクトリ` (ESRF 型 3GLS 設計、サブディレクトリ含む) と `ユーザーの pf2 ディレクトリ` (PF リング) の
スクリプトを `python run_sad.py lattice script` でそのまま動かすこと。ユーザーのディレクトリは読み取り専用に扱い
(スクリプトが datas.txt に追記する)、実行はすべて `<repo>\<work-esrf>` / `<work-pf2>` のコピーで行う。

- SAD 参照実行: `wsl -d Ubuntu -- bash -lc 'cd ".../<work-esrf>" && <SAD>/bin/gs X_full.sad > X.sad_r.out'`
  (`X_full.sad` = ラティス + スクリプトを cat で連結、`NPARA=…` は `NPARA = 1` に書き換え)
- pysad 実行: `PYTHONIOENCODING=utf-8 python -u run_sad.py pysad/tests/kek-ls_v03_12.lat <work-esrf>/X > <work-esrf>/X.pysad_r.out`
  (ラティスは esrf: `20cells/v03_12` のコピー `pysad/tests/kek-ls_v03_12.lat`、pf2: `<work-pf2>/v30_10`。Get[] はスクリプトのディレクトリ基準)
- 集計スクリプトは scratchpad (classify.py / ffscount.py / unsupported.py)。集計は正規表現ベースなので ±数 % の誤差がある。

## 1. ファイル分類

### 1.1 esrf (1398 ファイル、6 サブディレクトリ)

| ディレクトリ | 役割 | ラティス | スクリプト | modules_*.n | その他 |
|---|---|---|---|---|---|
| `esrf/` (トップ) | 現行の作業場所。v05_61 (DQBA 20 セル)、v09_01/02 (5 GeV) と約 100 本のスクリプト | 5 | 96 | 10 | PF_v9_21_simple, SPring8 (比較用ラティス), datas.txt |
| `16cells/` | 16 セル版の設計履歴 (v04_01〜v04_52 等 94 本) | 94 | 150 (+1 混在) | 11 | phase |
| `20cells/` | 20 セル版 (v01_01b〜v05_xx、matching05_* 等)。<work-esrf> の v03_12 はここのコピー | 144 | 128 (+4 混在) | 10 | *.vss, datas.txt |
| `20cells_HBR/` | 高輝度版 (v08_01〜v09_xx)。v08_20〜 はラティスとスクリプトが 1 ファイルに混在 (9 本) | 53 | 151 (+9 混在) | 13 | *.vss |
| `20cells_KEK-LS/` | KEK-LS 検討版 (v05_xx) | 68 | 101 | 10 | *.ps, *.vss |
| `20cells_old/` | 旧版 (v01/v02/v03、Octa 検討) | 173 | 116 (+7 混在) | 9 | MAD 出力 (.mad/.log/.print) |
| `simone/` | MAD-X / elegant 形式 (SAD ではない) | 0 | 0 | 0 | 4 |
| **合計** | | **537** | **742 (+21 混在)** | **63** | 約 35 |

スクリプトの命名: `disp*` (光学表示・派生計算)、`chrom*` (六極で色収差合わせ)、`matching*` (セル/リングの線形マッチング)、
`tune_survey*` (チューンサーベイ)、`amp_dep_ts*.sad` (振幅依存チューンシフト)、`dx_sv*/dy_sv*.sad` (マッチング条件を振るサーベイ)、
`correct_cod*.sad` (誤差 + COD 補正 + DA)、`calc_*.sad` (H 関数・色収差・エミッタンス・LMA・放射積分)、`print_*.sad`、`da_xy`、`touscheck*`、
`go_N.sad` (バッチ実行の束ね)、`draw_cad*.sad`/`drawmagvw3d.n` (作図)。modules_*.n はどのディレクトリにもほぼ同じ 9〜13 本が複製されている
(modules_general.n = WriteGraphs/DAsv、modules_amp_dep_ts.n、modules_correctoptics(_pf).n = β/η 応答と補正、modules_matching(_noID).n、
modules_tune_sv*.n、modules_write_correctdata.n、modules_measure_cod/eta.n (pf2)、modules_ID.n (pf2))。

### 1.2 pf2 (277 ファイル、サブディレクトリなし)

| 種別 | 数 | 内容 |
|---|---|---|
| ラティス | 136 (+7 混在) | v9_16〜v9_26 (現行 PF)、v10_01〜v10_19、v29_46、v30_0〜v30_51 (PF-II 4BA 14 セル)、v31_01〜v31_30、v32_01〜03、v40_01、loco_lattice_wig_on/off(_injN).sad、linac-pfbt.sad、ring_bt_lattice |
| スクリプト | 103 | disp*, chrom, matching* (v30_xx/v31_xx ごと), tune_sv*.sad, correct_cod*.sad, amp_dep_ts*, calc_*.sad, bump*, injmulti, kick, response, da_*, sad2ele*/sad2madx.sad/sad2at_python.sad (変換) |
| modules_*.n | 14 | esrf と同系 + modules_ID.n, modules_matching.n, modules_measure_cod.n, modules_measure_eta.n, modules_study.n |
| その他 | 24 | elegant (.lte/.ele/.twi/.mag/.mpr/.rf), pf.seqx, pflattice.py, .gif/.ps, datas.txt |

pf2 の USE 名: `RING`, `NN1`, `NCELL`, `QUARTER`, `STR`, `N2` などは v30_10 以降 (v30_10, v31_01〜) に定義があり、`<work-pf2>/v30_10` で全部解決できる。

## 2. 使用頻度と pysad の対応

対象: スクリプト 845 本 + modules 77 本 (esrf+pf2)。コメント (`!`, `(* *)`) は除外。
FFS コマンドは `FFS["..."]` 文字列内とベタ書きの両方を数え、大文字小文字を区別しない (回数 / ファイル数)。

### 2.1 FFS コマンド

| コマンド | 回数 | ファイル | pysad | コマンド | 回数 | ファイル | pysad |
|---|---|---|---|---|---|---|---|
| cell (calc/emit を含む) | 4778 | 876 | 対応 | fitp | 311 | 285 | 対応 (FITP>1 は遅い) |
| fit (+ NX/NY/BX/BY/AX/AY/EX/EPX … 条件) | 4480 | 428 | 対応 | coup | 252 | 210 | 対応 |
| stop | 2926 | 856 | 対応 | ins | 243 | 120 | 対応 (不安定域の扱いは §4) |
| free | 1727 | 426 | 対応 | codplot | 211 | 103 | 受理 (表示なし) |
| on/off (RAD FLUC COD RFSW RADCOD LOG …) | 1703 | 844 | 対応 | sfa/sda/sf/sd … (要素値設定・MAX/MIN) | 487 | 230 | 対応 |
| disp | 1447 | 683 | 対応 | fitfunction (FitFunction:=) | 172 | 114 | 対応 |
| use | 1424 | 842 | 対応 | gcut | 136 | 136 | 対応 |
| type | 1364 | 496 | 対応 | axi/bxi/ayi/byi/exi/epxi … (入口 Twiss) | 530 | 105 | 対応 |
| dp / dp0 | 921 / 368 | 525 / 172 | 対応 | mon / mcod / tcod | 108 / 89 / 89 | 108 / 89 / 89 | 受理のみ (SAD のこのビルドでも MON はエラー) |
| go | 849 | 426 | 対応 (多解問題で SAD と別解になることあり → §4) | dx/dy/dpx/dpy (初期軌道) | 214 | 90 | 対応 |
| dtheta / dk / delx / dely | 844 / 839 / 771 / 608 | 144 / 143 / 144 / 144 | 対応 | emit (単独) | 62 | 56 | 対応 |
| save | 758 | 395 | 対応 | fix | 54 | 18 | 対応 |
| npara (変数代入) | 662 | 620 | 無視 (並列なし) | draw | 27 | 17 | 無視 (警告表示) |
| nx/ny (fit 条件の短縮形) | 311 / 311 | 285 | 対応 | trpt / r1-r4 / reset / dump / geo / ncell | 各 ≤8 | ≤8 | 概ね対応 (geo は未確認) |

### 2.2 SADScript 関数 (組込み)

| 関数 | 回数 | ファイル | pysad | 関数 | 回数 | ファイル | pysad |
|---|---|---|---|---|---|---|---|
| FFS[] | 16805 | 760 | 対応 | SeedRandom | 192 | 106 | 対応 (SAD と同じ乱数列) |
| Twiss[] | 8973 | 792 | 対応 (`Twiss["NX","***"]` は 1 要素リストを返す → §4) | Sinh/Cosh/Tan/ArcCos | 281 | 40 | 対応 |
| LINE[] | 7138 | 607 | 対応 | Floor | 154 | 86 | 対応 |
| Print | 5670 | 845 | 対応 (書式は $FORM に従う) | CalculateOptics | 134 | 134 | 対応 (Monitor/Steer と組の定型が通る) |
| If / Do / Which / Switch / While | 3530 / 3296 / 45 / 25 / 4 | 529 / 742 / … | 対応 | Monitor / Steer | 132 / 132 | 132 / 132 | 対応 (定型の呼び出しが通る) |
| Length / Abs / Sqrt / Sin / Cos | 2382 / 2134 / 1789 / 1181 / 695 | | 対応 | Eigensystem / DiagonalMatrix / Sort | 109 / 109 / 109 | 109 | 対応 (Eigensystem は teigen の直訳。Sort[list, Greater] 未確認) |
| Write / WriteString / OpenAppend / OpenWrite / OpenRead / Read / Close | 1516 / 84 / 790 / 115 / 77 / 226 / 697 | 481 … | 対応 (ファイル入出力は Python ネイティブ。失敗時の SAD の挙動は移植しない) | Flatten / Join / Thread / Append / Transpose | 86 / 52 / 42 / 195 / 218 | | 対応 |
| Table / Range | 1469 / 246 | 265 / 240 | 対応 | Random / GaussRandom | 86 / 59 | 86 / 11 | 対応 |
| Get | 1081 | 719 | 対応 | Emittance[] | 64 | 35 | 対応 |
| DynamicApertureSurvey (ユーザー定義 DAsv[] 経由) | 558 (DAsv 537) | 272 (350) | 対応・**遅い** (5000 周 ×969 粒子 ≈ 50 分級) | Module / Head / Scan / Map / Select | 64 / 42 / 21 / 14 / 14 | ≤8 | 対応 (Scan は未確認) |
| FitValue (ユーザー定義の制約フック) | 536 | 197 | 対応 | TouschekLifetime | 34 | 6 | 対応 (Touschek 表は一致、INTRA 後の報告ブロックの差は §4) |
| PolynomialFit | 30 | 16 | 対応 (SAD の Fit.n をそのまま評価) | TrackParticles | 24 | 15 | 対応 |
| System / TkWait / ListPlot / Rectangle / Graphics / Show / Canvas / Text / Plot | 50 / 51 / 21 / 27 / 6 / 6 / 4 / 8 / 1 | ≤35 | TkWait は no-op、作図系は**未対応** | CaGet / CaPut (EPICS) / Sleep / SetElement / StringFill | 10 / 4 / 5 / 4 / 1 | ≤4 | 未対応 (実機制御・対象外) |

ユーザー定義関数で頻度の高いもの (modules_*.n): WriteGraphs 434/292、RMSWrite 420/143、CodCorr 325/109、MatchingDef 252/84、
CalculateCODResponse 203/103、RestoreMagnetValue 156/54、AmplitudeDependentTuneShift 135/132、SetInitialParameters 127/126、
MatchingARCA 127/84、InitMag 114/82、MatchingWholeRing 108/54、ChromCorrect 55/55、LinearAmpDepTSCorrect 90/30。

## 3. 実行結果 (pysad vs SAD)

86 本を `<work-userscripts>/<group>/<name>.sad` (前半ラティス + 後半コマンドの一体化ファイル。SAD の gs にそのまま渡せる) と `<name>.py` (同じ .sad を pysad で実行) にして、
同じフォルダで SAD (WSL の gs、1 本 300 s) と pysad (py -3.12、最低優先度、1 本 600 s) を走らせ、`cmp_out.py` で数値を比べた結果。
道具と間引きの規則は `<work-userscripts>/README.md`、実行の記録は `docs/history/work_recheck2.md`。

- 間引き (SAD/pysad 共通、動作テストのため): NPARA=1、DA 100 周・運動量点 2dp、誤差の実現 3/3/2、チューンサーベイ 3×3、DP0 走査 0.01、振幅走査 100 刻み、VOLT 4 点、PBUNCH 3 点。
  TkWait[] は削除、`FFS["draw …"]` は `Print["nodraw …"]` (WSL の SAD が止まる・落ちる)。
- 比較: 入力のエコー・バナー・In[]/Out[] を除き、数値を # に置き換えた骨格を difflib で揃え、対応する行の数値を比べる。|差| < 1e-9 と SAD の表示桁 1 ulp 以内は無視、
  相対差 1e-6 超を「不一致」と数える。`Momentum compact. alpha` は比べない (表示桁の差)。
- 判定 (`status_table.py`): rc=124 → 打ち切り、rc≠0 か Traceback → 例外、不一致 0 かつ骨格 ≥ 80 % → 数値一致、それ以外 → 一部一致。`notes.txt` で上書き
  (雑音だけの不一致は数値一致に、USE 先がラティスに無いものは SAD 側無効に)。表の列: SAD と pysad の実行秒、骨格の一致率、数値行 (比較 / 不一致)、最大相対差。

### esrf (52 本)

| スクリプト | SAD s | pysad s | 骨格 | 数値行 (比較 / 不一致) | 最大相対差 | 判定 | メモ |
|---|---|---|---|---|---|---|---|
| amp_dep_ts | 1 | 10 | 75% | 49 / 8 | 1.0e+00 | 一部一致 | 振幅走査の値は一致。最後の cell calc emit で SAD は不安定 (Twiss 1e200)、pysad は安定: 先頭 MARK の DX≈5 mm から始める閉軌道探索の収束の差 |
| amp_dep_ts_dp | 1 | 5 | 96% | 64 / 2 | 9.5e-01 | 一部一致 | dp=−0.045 で値が違い (48.047 と 47.989)、−0.04 で SAD は qcod が溢れて 1/cosμ の値 (2.38)、pysad は NaN。不安定域の縁。他の dp は一致 |
| amp_dep_ts_dz | 1 | 4 | 90% | 50 / 20 | 1.3e-04 | 一部一致 | dz 走査の値が 1.3e-4 違う (10/16 から既知の小さな数値差、原因未調査) |
| amp_dep_ts_q | 1 | 5 | 93% | 32 / 0 | 0.0e+00 | 数値一致 |  |
| amp_dep_ts_x | 1 | 30 | 89% | 58 / 18 | 1.3e+00 | 一部一致 | 同上 (i=201 は両方とも不安定) |
| amp_dep_ts_x_octa | 1 | 1 | 83% | 31 / 0 | 0.0e+00 | 数値一致 |  |
| amp_dep_ts_y | 2 | 10 | 75% | 49 / 8 | 1.0e+00 | 一部一致 | 同上 |
| calc_H | 1 | 2 | 95% | 95 / 0 | 0.0e+00 | 数値一致 |  |
| calc_H_2 | 0 | 1 | 5% | 0 / 0 | 0.0e+00 | SAD 側無効 | USE 先 (TEST2 等) がラティスに無く SAD も無意味な出力 |
| calc_LMA | 1 | 1 | 99% | 350 / 0 | 0.0e+00 | 数値一致 |  |
| calc_RadInt | 1 | 2 | 93% | 31 / 0 | 0.0e+00 | 数値一致 |  |
| calc_chromaticity | 1 | 3 | 95% | 55 / 0 | 1.0e-09 | 数値一致 |  |
| calc_emit | 1 | 2 | 93% | 33 / 0 | 0.0e+00 | 数値一致 |  |
| chrom | 1 | 8 | 89% | 51 / 20 | 1.3e-01 | 一部一致 | FITP>1 の GO で SAD が再起動する (比較は途中の値)。2 段目以降は多解 |
| chrom_b | 0 | 2 | 86% | 30 / 0 | 0.0e+00 | 数値一致 |  |
| chromaticity_plot | 1 | 3 | 94% | 43 / 0 | 0.0e+00 | 数値一致 |  |
| correct_cod | 2 | 11 | 83% | 32 / 0 | 0.0e+00 | 数値一致 |  |
| correct_cod_nomon | 6 | 23 | 48% | 88 / 47 | 1.0e+00 | 一部一致 | 誤差 (ae=5E-5, fe=5E-4, te=1E-4) を入れた後は SAD も pysad も閉軌道が見つからず (qcod Overflow 85 回ずつ) 不安定。SAD は最後に qtwiss も溢れて Twiss が 0、pysad は直前の値 (47.15) のまま。誤差なしの部分と流れは一致 |
| correct_cod_nomon_dq | 6 | 25 | 48% | 88 / 47 | 1.0e+00 | 一部一致 | 同上 |
| disp | 1 | 2 | 100% | 809 / 0 | 0.0e+00 | 数値一致 |  |
| disp2 | 47 | 11 | 85% | 68 / 1 | 1.5e-01 | 一部一致 | DA の地図 (RAD FLUC 100 周のカオス) |
| disp_L | 0 | 1 | 2% | 1 / 0 | 0.0e+00 | SAD 側無効 | USE 先 (ARCAL) がラティスに無い |
| disp_R56 | 2 | 21 | 100% | 2005 / 0 | 0.0e+00 | 数値一致 |  |
| disp_arc_optics | 0 | 1 | 86% | 19 / 0 | 0.0e+00 | 数値一致 |  |
| disp_b | 0 | 2 | 55% | 62 / 0 | 0.0e+00 | 数値一致 | SAD は FFS::error の後に式の形 (... in Do[...]) と定義の echo を出すので骨格 55%。数値 62 行は一致。Do を抜ける挙動も同じ (C25) |
| disp_bunchlength | 1 | 3 | 93% | 33 / 0 | 0.0e+00 | 数値一致 |  |
| disp_da | 39 | 9 | 85% | 67 / 0 | 0.0e+00 | 数値一致 |  |
| disp_emit | 1 | 2 | 92% | 30 / 0 | 0.0e+00 | 数値一致 |  |
| disp_eta | 1 | 1 | 4% | 1 / 0 | 0.0e+00 | SAD 側無効 | USE 先がラティスに無い |
| disp_matching_tune | 1 | 2 | 94% | 38 / 0 | 0.0e+00 | 数値一致 |  |
| disp_matrix | 1 | 2 | 95% | 200 / 0 | 0.0e+00 | 数値一致 |  |
| disp_optics | 1 | 1 | 2% | 1 / 0 | 0.0e+00 | SAD 側無効 | USE 先 (ARCINJ) がラティスに無い |
| disp_test | 178 | 18 | 82% | 77 / 3 | 3.2e-01 | 一部一致 | DA の地図と Score (RAD FLUC 100 周のカオス) |
| dx_sv | 0 | 1 | 49% | 19 / 0 | 0.0e+00 | 数値一致 | SAD の ... in の表示と、SAD 側のスクリプトの構文エラー (General::clquote) の分だけ骨格が低い。数値は一致 (不安定域は両方 NaN/1e100) |
| dy_sv | 0 | 2 | 78% | 20 / 0 | 0.0e+00 | 数値一致 | 同上 |
| go_1 | 8 | 1 | 0% | 0 / 0 | 0.0e+00 | 対象外 | System["go …"] は実行しない (SAD も sh: go: not found) |
| matching | 1 | 2 | 84% | 42 / 0 | 0.0e+00 | 数値一致 |  |
| matching21 | 0 | 1 | 79% | 23 / 0 | 0.0e+00 | 数値一致 | General::invset で Do を抜ける挙動が同じ。骨格 79% は ... in の表示の差 |
| matching3_01 | 1 | 10 | 81% | 191 / 24 | 9.4e-01 | 一部一致 | 2 段目以降の fit が別解 (多解)。1 段目は一致 |
| matching_geo_hmba | 1 | 1 | 79% | 15 / 0 | 0.0e+00 | 数値一致 | 骨格 79% は SAD の再起動行の分。数値は一致 |
| matching_test | 1 | 2 | 57% | 4 / 0 | 0.0e+00 | 数値一致 | 出力 7 行、数値は一致 |
| print_all_error | 1 | 4 | 93% | 32 / 0 | 0.0e+00 | 数値一致 |  |
| print_cod_magnitude | 4 | 25 | 97% | 109 / 10 | 1.0e+00 | 一部一致 | ae=5E-6〜5.5E-5 の 11 点 (実現 3 つずつ) のチューンのずれは一致。ae≥6E-5 は不安定 (NaN、チューンのずれ 0.3〜0.6) で SAD と pysad の壊れた値が違う (SAD は qcod 溢れ後の停留値 0) |
| print_emittance | 1 | 10 | 83% | 31 / 0 | 0.0e+00 | 数値一致 |  |
| tH | 1 | 1 | 79% | 10 / 0 | 0.0e+00 | 数値一致 | 出力が短く骨格 79% だが数値は全部一致 |
| tH2 | 1 | 1 | 73% | 8 / 0 | 0.0e+00 | 数値一致 | 同上 |
| touscheklife | 1 | 2 | 96% | 59 / 0 | 0.0e+00 | 数値一致 |  |
| tune_survey | 7 | 38 | 100% | 636 / 0 | 0.0e+00 | 数値一致 |  |
| tune_survey_disp | 6 | 38 | 100% | 725 / 0 | 0.0e+00 | 数値一致 |  |
| tune_survey_linear | 14 | 109 | 100% | 12793 / 0 | 0.0e+00 | 数値一致 |  |
| tune_survey_momentum | 0 | 1 | 0% | 0 / 0 | 0.0e+00 | 例外 | 入力ファイル (hsxdata) が無い (方針により追わない) |
| tw_frac | 1 | 1 | 96% | 64 / 0 | 0.0e+00 | 数値一致 |  |

### pf2 (34 本)

| スクリプト | SAD s | pysad s | 骨格 | 数値行 (比較 / 不一致) | 最大相対差 | 判定 | メモ |
|---|---|---|---|---|---|---|---|
| amp_dep_ts_dp_r | 1 | 3 | 95% | 72 / 0 | 0.0e+00 | 数値一致 |  |
| amp_dep_ts_x_r | 1 | 2 | 96% | 62 / 0 | 0.0e+00 | 数値一致 |  |
| beamsize_r | 1 | 3 | 99% | 370 / 0 | 0.0e+00 | 数値一致 |  |
| bump2_r | 1 | 2 | 95% | 45 / 0 | 0.0e+00 | 数値一致 |  |
| bump3_r | 0 | 1 | 94% | 42 / 0 | 0.0e+00 | 数値一致 |  |
| calc_LMA_r | 1 | 2 | 100% | 1086 / 0 | 0.0e+00 | 数値一致 |  |
| calc_bump_r | 0 | 1 | 89% | 44 / 0 | 0.0e+00 | 数値一致 |  |
| calc_chromaticity_r | 1 | 2 | 94% | 55 / 0 | 1.8e-10 | 数値一致 |  |
| calc_emit_r | 1 | 1 | 90% | 30 / 0 | 0.0e+00 | 数値一致 |  |
| chrom_r | 1 | 2 | 78% | 30 / 0 | 0.0e+00 | 数値一致 | SAD が FITP>1 の GO で再起動 (出力に Run time Environment が混じる) ので骨格 78%。出た数値は一致 |
| correct_cod_9_16_r | 1 | 1 | 87% | 68 / 1 | 4.3e-02 | 一部一致 | DA の地図 1 文字 (RAD FLUC 20 周、落ちる周の乱数のカオス)。COD 補正量は一致 |
| da_xy_r | 19 | 10 | 81% | 40 / 1 | 4.3e-02 | 一部一致 | DA の Score 110 と 115 (RAD FLUC 50 周のカオス、1〜3 %)。pysad は DA の表示行が多い |
| disp | 1 | 1 | 84% | 19 / 0 | 0.0e+00 | 数値一致 |  |
| disp_bt_r | 1 | 1 | 100% | 3549 / 0 | 0.0e+00 | 数値一致 |  |
| disp_bunchlength_r | 1 | 10 | 98% | 161 / 0 | 0.0e+00 | 数値一致 |  |
| disp_bx_by_ex_r | 1 | 1 | 100% | 1086 / 0 | 0.0e+00 | 数値一致 |  |
| disp_kick_effect_r | 1 | 2 | 92% | 38 / 0 | 0.0e+00 | 数値一致 |  |
| disp_q_K1_r | 1 | 1 | 93% | 44 / 0 | 0.0e+00 | 数値一致 |  |
| disp_response_r | 1 | 2 | 99% | 588 / 0 | 9.6e-07 | 数値一致 |  |
| error_ring_r | 1 | 1 | 100% | 928 / 4 | 2.0e+00 | 一部一致 | DP0=±0.025 の BY が MKMH.1 で負 (不安定域)。SAD は -0.330、pysad は +0.330 (符号)。安定域の誤差表 (924 行) は一致 |
| kick_r | 0 | 2 | 94% | 37 / 1 | 1.8e+00 | 数値一致 | K0 の残り 1e-10 rad (SAD -4.8E-10、pysad 5.9E-10): バンプの閉じ残りの雑音 |
| matching | 0 | 1 | 83% | 18 / 0 | 0.0e+00 | 数値一致 |  |
| matching02_r | 1 | 1 | 88% | 39 / 0 | 0.0e+00 | 数値一致 |  |
| matching_BMNCELL_r | 1 | 1 | 88% | 37 / 0 | 0.0e+00 | 数値一致 |  |
| matching_geo_ring_r | 1 | 2 | 89% | 33 / 0 | 0.0e+00 | 数値一致 |  |
| matching_v30_12_r | 1 | 2 | 96% | 176 / 18 | 9.3e-01 | 一部一致 | 2 段目の fit (K2) が別解。1 段目は一致 |
| matching_v31_23_r | 1 | 3 | 93% | 133 / 32 | 8.4e-01 | 一部一致 | 2 段目の GO 以降が別解 (多解)。1 段目・初回の emit は一致 |
| plot_inj_r | 1 | 1 | 91% | 34 / 0 | 0.0e+00 | 数値一致 |  |
| print_all_error_r | 1 | 2 | 92% | 30 / 0 | 0.0e+00 | 数値一致 |  |
| sad2ele_v10_16_r | 1 | 3 | 98% | 393 / 0 | 0.0e+00 | 数値一致 |  |
| sad2madx_r | 1 | 2 | 98% | 134 / 0 | 0.0e+00 | 数値一致 |  |
| touscheck_r | 1 | 18 | 98% | 174 / 0 | 0.0e+00 | 数値一致 |  |
| touscheklife_r | 1 | 31 | 99% | 345 / 0 | 0.0e+00 | 数値一致 |  |
| tune_sv_v30_24_1_r | 1 | 2 | 97% | 159 / 15 | 1.3e+00 | 一部一致 | 2 段目の反復が途中から分岐 (多解)。反復表の前半は一致 |

| 対象 | 本数 | 数値一致 | 一部一致 | 動くが不一致 | 例外 | 打ち切り | SAD 側無効 | 対象外 |
|---|---|---|---|---|---|---|---|---|
| esrf | 52 | 34 | 12 | 0 | 1 | 0 | 4 | 1 |
| pf2 | 34 | 28 | 6 | 0 | 0 | 0 | 0 | 0 |
| 合計 | 86 | 62 | 18 | 0 | 1 | 0 | 4 | 1 |

## 4. 残っている不一致・未対応 (種類別)

「本数」は 86 本のうちその原因が主因のスクリプト数。

| 種類 | 本数 | スクリプト | 内容 |
|---|---|---|---|
| fit の多解・FITP>1 | 5 | chrom, matching3_01, matching_v30_12_r, matching_v31_23_r, tune_sv_v30_24_1_r | 1 段目・初回の解は一致。2 段目以降の多解で別解に着く。SAD は FITP>1 の GO で再起動するので最終解の参照値がない |
| 不安定ラティス・閉軌道探索の差 | 8 | amp_dep_ts, amp_dep_ts_x, amp_dep_ts_y (最後の emit)、amp_dep_ts_dp (dp=−0.045, −0.04)、error_ring_r (MKMH.1 の BY)、correct_cod_nomon, correct_cod_nomon_dq (誤差を入れた後は両方とも閉軌道が見つからない)、print_cod_magnitude (ae≥6E-5) | SAD は qcod が溢れると壊れた値 (1e200、1/cosμ の値、qtwiss 溢れ後の 0) で続行、pysad は直前の値か NaN。安定な部分の値と流れは一致。先頭 MARK の DX≈5 mm から始める閉軌道探索は収束が分かれる |
| DA の乱数のカオス | 4 | disp2, disp_test, da_xy_r, correct_cod_9_16_r | RAD FLUC の DA の地図 1 文字・Score 1〜3 % (落ちる周の乱数)。NORAD の DA は一致 |
| 小さな数値差 | 1 | amp_dep_ts_dz (1.3e-4) | 原因未調査 |
| SAD 側無効・対象外 | 6 | calc_H_2, disp_L, disp_eta, disp_optics (USE 先がラティスに無い)、go_1 (System)、tune_survey_momentum (入力ファイルなし) | 比較しない |

解消したもの (11/06): 性能 (間引きで 600 s 打ち切りなし)、INTRA の EMIT 報告 (touscheck_r, touscheklife_r, disp_bunchlength_r, beamsize_r)、FFS[] のエラーで Do を抜ける (disp_b, chrom_b, matching21, dx_sv, dy_sv, amp_dep_ts_x_octa)、PBUNCH の既定値。

## 5. 補足 (ファイルの場所)

- `<work-userscripts>/esrf/`, `<work-userscripts>/pf2/`: 一体化した `<name>.sad` / `<name>.py`、modules_*.n、出力 `<name>.sad.out` / `<name>.py.out`、`sad_log.txt` / `py_log.txt`。
  `<work-userscripts>/manifest.txt` (元のラティス・スクリプトと効いた間引き)、`notes.txt` (判定の上書きとメモ)。
- 元のスクリプトの写し: `<work-esrf>/` (v03_12 = pysad/tests/kek-ls_v03_12.lat、correct_cod 系は v05_61)、`<work-pf2>/` (v30_10 ほか)。ユーザーのディレクトリ `ユーザーの SAD スクリプトのディレクトリ` は触らない。
- 10/16 の再点検の出力: `<work-recheck>/` (manifest.txt はそのまま `<work-userscripts>/build.py` の入力)。
- 小さな照合スクリプト: `<work-perf>/` (err_test.sad, err_test2.sad, use_test.sad, dx0_test.sad, codplot_test.sad など)。
