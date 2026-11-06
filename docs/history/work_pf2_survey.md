# 作業ログ: PF リング (pf2) スクリプト群の SAD / pysad 照合調査 (2026-10-15)

前担当が中断した調査の仕上げ。材料は `<work-pf2>/` の既存成果物 (`<name>_r`、`<name>_r_full.sad`、`.sad.out`、`.pysad.out`、`.sad_time`)。
コードは未変更、`ユーザーの pf2 ディレクトリ` は一覧と grep のみ (スクリプト実行なし)。新たに行った実行は pysad の小テスト 2 本 (後述) だけ。
pysad 出力は 12:53〜13:10 の一括実行 (`run_pysad_batch.sh`、1 本 300 s 制限) の結果で、その時点のコードの状態を反映する。

## 1. 要約

- 一括実行は `_r` 31 本 (+ 前回の chrom / disp / matching)。SAD は 29 本が `Good Bye` で終了し、異常終了は 2 本 (error_ring: `stop` 無しで入力終端、plot_inj: Tk 初期化失敗後に SIGSEGV)。
  ただし matching 系は SAD 側でも `???-Emittance[]-closed orbit not found` (NCELL は空洞なし) や `Missing value for S*` を出す。これはスクリプトの性質で、pysad の問題ではない。
- pysad は 31 本中、完走 13、例外 16、打ち切り 2 (calc_bump=300 s 超過、matching_BMNCELL=NaN 警告後に無出力)。
- 完走 13 本のうち数値まで SAD と一致したのは 4 本 (calc_emit、calc_LMA、disp_bx_by_ex、sad2madx。最大相対差は 1e-8 以下)。
  calc_chromaticity と amp_dep_ts_dp は例外までの部分が一致 (4e-8 以下)。
- 不一致の原因は大きく 3 つに集約される。(a) **BEND の K0 キックが COD/分散に効かない** (bump2/3、disp_response、disp_kick_effect の 4 本が「実行は通るが応答が 0」)。
  (b) **2 段目以降の fit (GO) の発散**。1 段目は SAD と NEWTON 残差が 5 桁一致するが、2 段目の初期残差が SAD の 10^6〜10^13 倍になり `optics calculation failed`。
  (c) 構文・関数の未実装 (`ON INTRA` の FFS 内使用、`Repeat/Until`、`CalculateOptics`、`Twiss["SIGX"]`、`GX` の fit 条件など)。

## 2. 実行結果の表

「一致」列は、SAD 出力から入力エコー (`_full.sad` に含まれる行) と cell calc 報告ブロックを除いた数値行を、数値を `#` に置換した骨格で difflib 整列し、
対応した行の数値を比較したもの (絶対差 1e-9 未満は無視)。報告ブロックは Real tune、Emittance、U0、alpha、dz を別途比較した。
値は最大相対差。pysad が例外で止まる場合は止まるまでの部分のみ。時間は pysad の実測 (SAD は全て 0.6〜4 s)。

| スクリプト (ラティス) | 目的 / 主な FFS・関数 | SAD | pysad (時間) | Print 一致度 |
|---|---|---|---|---|
| calc_emit (v30_10) | 光学・エミッタンス表示 `cell calc emit`、`LINE["NAME","B*"]` | 正常 | 正常 (12 s) | 一致。報告値は一致、σδ のみ表示桁の差 (.00134 と .00134253) |
| calc_LMA (v30_10) | 0.5 刻みで `Twiss["BX"/"EX"]` を Print (約 1060 行) | 正常 | 正常 (13 s) | 1057 行対応、最大 9e-13 |
| disp_bx_by_ex (v30_10) | 同様の BX/BY/EX 表 | 正常 | 正常 (10 s) | 1057 行対応、最大 7e-11 |
| calc_chromaticity (v30_10) | DP0 を ±1 % で振り tune を Print、`PolynomialFit` | 正常 (PolynomialFit は 1e10 級の係数で、スクリプト側の問題) | SADError 不安定ラティス: y 面の sin(μy) ≒ 0 (DP0=+0.005 で停止) | 停止までの 15 行は最大 1.6e-9。以降は欠落 |
| amp_dep_ts_dp (v30_10) | 上と同型の DP0 スキャン (`modules_amp_dep_ts.n`) | 正常 (SAD も ν_y が 7.0173→7.0448 と整数共鳴近傍をまたぐ) | 同じ SADError (80 s) | 37 行対応、最大 4e-8。DP0≧0.005 は欠落 |
| amp_dep_ts_x (v9_16) | `AmplitudeDependentTuneShift` (追跡由来の tune) | 正常 | 正常 (30 s) | **不一致**: νx 9.54507 と 9.52056、νy 5.24950 と 5.27059 (相対 2.6e-3 / 4.0e-3)。3 行とも同傾向。線形値 νx=9.54498 に近いのは SAD |
| disp_bunchlength (v30_10) | `VOLT` を振り `BucketHeight`、`BunchLength` (`/.Emittance[]`) | 正常 | 正常 (89 s) | 4 行。BucketHeight は 2e-10 で一致、BunchLength は 1.1e-4〜3.8e-4 ずれ (j=10: .0051159251 と .0051139968) |
| da_xy (v30_10) | `DynamicApertureSurvey` の角度スキャン、`Repeat … Until` | 正常 (`type s*` は Missing value、スクリプト側) | SADError Part: index 2 out of range (`Repeat`/`Until` が未定義で `res` が空) (15 s) | 最初の DA Score が 107 と 110 (2.7 %)。以降は欠落 |
| error_ring (v30_10) | `ERROR` で誤差を入れ COD/ETA RMS | **異常終了** (入力終端、`NaN,NaN,NaN,NaN` を出力) | TypeError: 'float' object is not subscriptable (`ffscmd.py` 1498 行 `_error_cmd` の `rng.gauss(1, gcut())[0]`) | 比較不能 |
| matching02 (NCELL) | `Fit NX`、`free Q*`、`free B K1`、`go` | 正常 (3 反復で Residual 4.072、未収束。Emittance[] は NaN) | 正常 (28 s) | GO ログの反復列は不一致 (最大 0.70)。最終 Residual 4.072 と 4.069、NX .20918 と .20932。同じ未収束解 |
| matching_v30_12 (v30_12) | NCELL で 2 段の `Fit NX/NY`、`go` | 正常 (2 段目は 8 反復で Residual 0.832) | RuntimeError optics calculation failed (`match6.py` 447 行、44 s) | 1 段目の残差 7.4750E-05 → 1.8141E-08 は完全一致。2 段目は初期残差 1.9e23 (SAD 4.0e6) で発散 |
| matching_v31_23 (v31_23) | 同型 (`DAsv` は定義のみ) | 正常 | 正常 (229 s) | 1 段目 Matched 2.6724E-11 一致。2 段目は別解 (U16 の BX 4.685 と 5.865、最大 1.03) |
| tune_sv_v30_24_1 (v30_24) | チューンサーベイ (2 重 Do、fit 付き) | 正常 | RuntimeError optics calculation failed (30 s) | 1 段目の残差列は一致 (1.7703E-06)。2 段目の初期残差 2.45e23 (SAD 2.6e10) |
| matching_BMNCELL (v30_33) | `LINE["ANGLE","BMA"]` を変え fit (`free BMA K1`) | 正常 (`LINE["ANGLE","BMB"]=` は invset エラー、Residual 4.9e-2) | **打ち切り** (1 段目 Matched 後 emit6d の overflow/NaN 警告のみ) | 数値行なし |
| matching_geo_ring (v9_16_sym) | `fit MKRC GX/GY`、`free BMA L` | 正常 (`go` が `No variable` エラー) | KeyError 'GX' (`match6.py` 361 行 `_tgfun`、1 s) | 数値行なし |
| disp_q_K1 (v31_24) | `LINE["K1",qlist[i]]` を Print | 正常 | SADError bad location Q1A (`{Q1A,…}` の裸の要素名シンボル) (9 s) | 数値行 0 |
| beamsize (v31_24) | `beamsize`、`Twiss["SIGX"/"SIGPX"…]` | 正常 | SADError Twiss: unknown function SIGX (9 s) | 表の行対応が取れず比較不能 (U* 行の並びが違う) |
| touscheck (v30_10) | `ON INTRA` + `Emittance[]` (Touschek 入力) | 正常 (CPU 1.5 s) | SADError unexpected token 'INTRA' (1 s) | 数値行 0 |
| touscheklife (v30_10) | 同上 | 正常 (CPU 3.1 s) | 同じ SADError | 数値行 0 |
| bump2 (v9_18) | `LINE["K0","BK3"]` を振り fit (`free BK1 K0`) | 正常 (DX −2.0e-8、DPX −1.5e-9) | 正常 (4 s) | **不一致**: K0 値の Print 以外、DX/DPX は 0 |
| bump3 (v9_25) | `USE ID16`、`ins calc`、`free K1 K0 K2 K0`、`type K*` | 正常 (K1 K0=.0020867、K2 K0=−.0023854) | 正常 (4 s) | `ins calc` 表は 4e-8 で一致。fit 結果の K0 は pysad に無い |
| kick (v9_25) | `K1 K0 0.0021;` (要素パラメータ代入文) と fit | 正常 | SADError unexpected token 'K0' (5 s) | 数値行 0 |
| disp_response (v9_25) | `LINE["K0",K1..K5]=0.04e-3` の PM 応答表 | 正常 | 正常 (4 s) | 534/559 行は形式一致だが DX 応答が全て 0 (SAD は PM013 で .2206 など) |
| disp_kick_effect (v9_16_inj) | `LINE["K0","KK1"]`、`ins calc` で INJ の DX | 正常 | 正常 (3 s) | 10 行とも DX が 0 (SAD は 11.49 など) |
| calc_bump (v9_16_ele) | `ins calc`、`fit INJ DX`、`free KK1..4 K0`、`go` | 正常 | **打ち切り** (300 s 超過)。`?Undefined command or element: B10` の後に無出力 | 数値行 0 |
| sad2madx (v9_16) | MAD-X 形式の書き出し (`LINE["…"]` 参照) | 正常 | 正常 (4 s) | 105/108 行対応、最大 4e-8。WIGWIG 行など 3 行が pysad に無い |
| sad2ele_v10_16 (v10_16_PV) | elegant lte 書き出し | 正常 | 正常 (6 s) | BEND 39 行が**不一致**。`LINE["E1"]` が SAD は .5 (係数)、pysad は .1122 (rad) |
| disp_bt (linac-pfbt.sad) | `USE LinacPFBT`、`TRPT`、`ins calc` | 正常 (7381 行) | `Missing beamline in USE` → SADError No beam line (1 s) | 数値行 0 |
| correct_cod_9_16 (v9_16) | ERROR + COD 補正 (`Eigensystem`)、DA | 正常 | NotImplementedError track6n: 長さ 0 の四極 (7 s) | 数値行 0 |
| print_all_error (v9_16) | `CalculateOptics`、COD/ETA 表 | 正常 | SADError CalculateOptics is not implemented (7 s) | 数値行 0 |
| plot_inj (v9_16_inj) | `GaussRandom` 1000 粒子の `TrackParticles`、作図 | **異常終了** (Tk 初期化失敗の後 SIGSEGV) | FloatingPointError 1 周目: track6: INJ で発散 [inf …] (4 s) | 比較不能 |
| chrom (v30_10, 前回) | `FITP`、DP=0.001、`FREE SD*/SF* K2` | 正常 (1 反復、残差 7.9e-8) | 900 s で打ち切り (GO 発散、NaN 警告) | 比較不能 |
| disp (NN1, 前回) | `disp $$$`、`draw` | 正常 | 正常 ($$$ の BX/NX/BY/NY が一致、`draw` は無視) | 一致 |
| matching (NCELL, 前回) | `Fit NX…go` | 正常 | 正常 (両者とも go の出力なし) | 一致 |

補足 (今回の追加確認、現行コードで実施)。
`FFS; ON INTRA;` は `unexpected token 'INTRA'` で失敗し、FFS の前に置いた `ON INTRA;` は通る (touscheck/touscheklife の原因)。
`LINE["K0","K1"]=0.0001` は値としては保持されるが、v9_25 の `cell calc` 後も PM013 の DX は 0 のままだった (K1 は ANGLE=0 の BEND)。

## 3. ユーザーの pf2 ディレクトリ の全ファイル (277 ファイル + サブディレクトリ 2018_10 が 1 つ)

| 種別 | 数 | 内容 |
|---|---|---|
| スクリプト | 103 | amp_dep_ts*(4)、beamsize、bump2/3、calc_*(7)、chrom、correct_cod*(7)、da_*(3)、disp*(15、うち disp_bt、disp_ringbt が BT 系)、error/error_ring、fit_kick、injmulti、kick、matching*(20)、plot_inj(2)、print_*(2)、response、sad2*(5)、test、touscheck、touscheklife、tune_sv*(20)、write_geometry9、draw_cad.sad、drawmagvw3d.n |
| modules_*.n | 15 | ID、amp_dep_ts(.n あり/なし)、amp_dept_ts (綴り違い、拡張子なし)、correctoptics、correctoptics_pf、general、matching、matching_noID、measure_cod、measure_eta、study、tune_sv、tune_sv_v03_35、write_correctdata |
| ラティス | 144 | v9_16〜v9_26_V (30)、v10_01〜v10_19 (22。v10_16_PV/_ocelot_cod を含む)、v02_11、v29_46、**v30_0〜v30_51_1GeV (47)**、**v31_01〜v31_30_NCB (31)**、**v32_01〜03**、v40_01、loco_lattice_wig_on/off(_inj2/3) (6)、linac-pfbt.sad、ring_bt_lattice |
| その他 | 15 | elegant (pf_v10_16.lte/.mag/.mpr/.rf/.twi、_ocelot_cod.lte、_out.lte、Twiss_pf.ele)、pf.seqx、pflattice.py、gif、a.vss、datas.txt、v31_03_STR.12.ps、v9_16_sym_STR.12.ps |

v30 系 (PF-II 4BA 14 セル) の主なもの: v30_10 (現行基準。RING/NN1/NCELL/N2/QUARTER/STR が全部解決する)、v30_12、v30_24、v30_33、v30_43 (空洞電圧違い 3 種)。
v31 系: v31_01〜v31_24 (NCELL 系の matching が v31_01/06/09/14/17/19/21/23 に対応)、v31_24 (beamsize、disp_q_K1)、v31_30_NCB。v32 は v32_01〜03 (matching_v32_01 のみ)。

## 4. スクリプトが USE する LINE と Get する .n (grep、コメント除外、103 本)

| USE した LINE | スクリプト数 | 主な使用者 |
|---|---|---|
| RING | 87 | ほぼ全部 (disp 系、correct_cod、bump、calc_*、tune_sv、matching の最終段) |
| NCELL | 35 | matching* 約 15 (v30_xx/v31_xx/v32_01、02、NCELLA)、tune_sv* 約 19 |
| N2 | 27 | matching_v30_xx/v31_xx/v32 (13)、tune_sv_hint/int/ko*/v30_24* (10)、disp 系 (4) |
| LCELL 5、RINGA 5、INJL 4、STR 3、NCELL2 3 | | beamsize、fit_kick、kick、disp_bump / tune_sv*A、calc_bump・disp_kick_effect、disp、matching_geo |
| ID16 2、NN1 2、QUARTER 2 | | bump3・disp_SQ、disp・matching_geo_NN1 |
| ARCA、EXBTINV、INJLINV、HBRING、LINACPFBT、RINGC、RINGH、NCELLA、ALL、TEST 各 1 | | calc_H、calc_bump、disp_kick_effect_K3K4、sad2ele_slitj、disp_bt、matching_BMNCELL、disp_RINGH、matching_NCELLA、disp_ringbt、disp_bump |

| Get した .n | スクリプト数 | 備考 |
|---|---|---|
| modules_general.n | 65 | WriteGraphs、DAsv。pysad で解決済み |
| modules_amp_dep_ts.n | 7 | amp_dep_ts*(4)、calc_H*(2)、calc_emit。実行確認済み |
| modules_write_correctdata.n | 2 | error、error_ring |
| modules_ID.n | 1 | da_tp |
| drawmagvw3d.n | 1 | draw_cad (Tk 描画) |
| modules_deal_coordinate.n | 1 | write_geometry9。**ディレクトリに存在しない** |
| `PFRBT:K1:ANGLE_SET_READ.VAL` など 8 件 | 2 | disp_bump、response。EPICS PV 名で、オフラインでは実行不能 |

## 5. 未対応機能の優先順位 (頻度 × 難度)

頻度は pf2 のスクリプト 103 本中で、その機能を使う本数 (grep、コメント除外)。難度は実装の見積もり。

| 順位 | 機能 | 頻度 | 難度 | 症状・メモ |
|---|---|---|---|---|
| 1 | BEND の K0 キック (`LINE["K0",x]=`、`FREE x K0`、`x K0 v;`) を COD/分散/fit に反映 | 26 本 | 低〜中 (推定) | bump2/3、kick、disp_response、disp_kick_effect、calc_bump、correct_cod*、plot_inj の基礎。現状は無言で 0 応答になり、誤った数値が出る点が最も危険 |
| 2 | 2 段目以降の fit の収束 (SAD と同じ NEWTON 軌跡)。ステップ制限・線探索、解析微分 (`work_analytic_deriv.md`) | 36 本 (matching* 16 + tune_sv* 20)、chrom も | 高 | 1 段目は SAD と一致する。2 段目で残差が 10^6〜10^13 倍から始まり `optics calculation failed`。高速化 (matching_v31_23 は 229 s) も同じ課題 |
| 3 | ERROR + COD 補正系: `CalculateOptics`、track6n の長さ 0 四極 (SQ)、`ERROR` の gauss 引数 | 14 本 (CalculateOptics 9) | 中 | correct_cod*(7)、print_*(2)、error*(2)、injmulti。pysad 側のバグ (TypeError) も 1 件 |
| 4 | DynamicApertureSurvey の一致度と `Repeat/Until` | 56 本が DAsv を定義、実使用は da_xy ほか数本。Repeat は 1 本 | 中 | DA Score が 2.7 % ずれる。`Repeat/Until` は構文のみで低難度 |
| 5 | `ON INTRA` の FFS 内受理 (touschek.py は別途存在) | 4 本 | 低 | 構文だけ。通れば touscheck/touscheklife/disp_bunchlength に進める。disp_bunchlength の BunchLength 1e-4 差の原因も要調査 |
| 6 | fit 条件・LINE での幾何量 `GX/GY` | 6 本 (matching_geo* 4、disp_RINGH、draw_cad) | 中 | `geo.py` はあるが fit 条件 `_ROW_ATTR` に無い。SAD 側も `go` で No variable になるものがある |
| 7 | `LINE["E1"/"E2"]` の返値を SAD の係数 (.5) に | 7 本 (sad2* 系の変換スクリプト) | 低 | 現状はラジアン値を返して elegant/MAD-X 出力が狂う |
| 8 | 個別: 裸の要素名シンボル、`USE LinacPFBT` の解決、`Twiss["SIGX"…]`、ADTS の追跡 tune の差 (0.02) | 各 1〜3 本 | 低〜中 | disp_q_K1、disp_bt/disp_ringbt、beamsize、amp_dep_ts*(4) |
| 対象外 | Tk/plot/EPICS 依存 (plot_inj の作図、draw_cad、response、disp_bump) | 約 5 本 | — | SAD も WSL では Tk 失敗で SIGSEGV |

着手順の提案: 1 (小さく効果大) → 5・7 (構文と返値) → 2 (本丸) → 3 → 4。
注意: 比較の整列は骨格一致ベースなので、表の「対応行数」は目安。数値の最終判断は各 `.out` を直接見ること。
