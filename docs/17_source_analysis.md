# 17. SAD ソース解析 総合リファレンスと pysad 対応表

最終更新: 2026-11-06 (現状のスナップショット)

対象 SAD は K64 版 (https://github.com/KatsOide/SAD)、コマンドの意味は CERN SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html) を先に見る。数値の根拠はローカルの `<SAD>\src`。行番号は 2026-10 時点。
関連: 01_overview (モジュール地図)、10_benchmarks (Level 1〜16)、14_fortran_status (Fortran ファイル別の状況)、00_RESUME (時系列)、work_*.md (各モジュールの作業ログ)。
「照合桁」は SAD (K64) の印字精度 (10 桁、追跡は 5e-11) との最大差で、ベンチ名を併記する。

---

## 1. 概要

SAD (Strategic Accelerator Design) は KEK 発の加速器光学・追跡コード。固定形式 Fortran (module/select case 等 F90 構文を多用) と一部 C/C++ のコアに、Mathematica 風の SADScript を持つ。

| 項目 | 実測 |
|---|---|
| `src/*.f` | 474 本、計 141,678 行 (`cat src/*.f \| wc -l`) |
| `src/*.c` / `src/sim/` | 38 本 / 42 ファイル (.c 20、.h 16 ほか) |
| `src/inc/` | 35 エントリ (.inc 33 + changelog + k64.readme)。大半は全行コメントアウト |
| `Packages/` | 82 エントリ (.n 73、.py 6 ほか) |
| `mk/` | .mk 19 本 + Makefile 2 本 |
| pysad (現行系統) | `pysad/*.py` 36,821 行 (非推奨の旧 API ffs/matching/tracking/optics 3,911 行を含む)。pytest 277 件 |

ファイル分類: 🔴A 物理コア / 🟡B 数学 / 🔵C インタープリタ / 🟢D グラフィックス / ⚫E ビルド / 🟣F SADScript パッケージ / 🟤G 外部ライブラリ。pysad の対象は A と B (C は `sadscript.py` + `ffscmd.py` で最小限を自前実装、D/E/G は不要)。

## 2. フォルダツリー (実在するもののみ)

```
<SAD>/
├ src/   *.f 474 / *.c 38 / *.h 9
│   ├ tftrack.f tturn.f track.f tturne.f tbend.f tquad.f tsol.f ... (追跡・行列)
│   ├ qcell.f qtwiss.f tffscalc.f temit.f INTGRL.f                  (CALC / EMIT / RADINT)
│   ├ tfstk.f tfetok.f tfeval.f tfestk.f tfeeval.f tfefun.f tfefun1.f (インタープリタ)
│   ├ calc.y/calc.c  SAD1 MAIN の数値式 (SADScript とは無関係)
│   ├ inc/  gdtoa/(51)  sim/(OS 抽象化 C 層)
├ Packages/ init.n Optics.n emit.n dapert.n Fit.n correction.n optimize.n ... Plots 系 / Tkinter 系 / CaSad 系 / sad.py
├ extensions/ RandomMT/  sad+/ (大見氏の DA, C++ 33 本)  Standard/ (Tkinter, Xlib)
├ contrib/ gdtoa BBBrem_dir libtai-0.60 tcl8.6.11 tk8.5.18 tk8.6.11 btkshim patches tools
├ mk/ (sad.obj.mk ほか)  config/ (GCC.spec)  sad.conf  Makefile
```

## 3. アーキテクチャ

### 3.1 6 層

| 層 | 内容 | pysad 側 |
|---|---|---|
| 1 SADScript (.n/.sad) | ユーザスクリプト・格子 | `sadscript.py` |
| 2 パーサ | tfetok (字句) + tfeval/tfestk (演算子順位)。式木はリスト。**バイトコードやスタックマシンではない**。calc.y は無関係 | `sadscript.py` (Pratt 構文) |
| 3 評価器 | tfeeval (式木の再帰評価)、tfefun/tfefun1 (組込関数、id>1000 は tfefun1 の computed goto) | `sadscript.py` / `ffscmd.py` の組込表、`sadlist.py` (リスト系組み込み) |
| 4 FFS | tffs / tffsa (コマンド分岐) / tffscalc | `ffscmd.py` |
| 5 物理コア Fortran | tturn/tturne/qtwiss/temit/... | 下記各モジュール |
| 6 C OS 層 | sim/ (メモリ、シグナル、動的リンク、X11/Tk 薄グルー) | 不要 |

### 3.2 NaN タギング

全 SADScript オブジェクトは 64 bit ディスクリプタ。実数以外は IEEE double の NaN 空間に型とアドレスを埋める (`tfstk.f:822-823`: ktfmask = z'fffe000000000000'、ktamask = z'0001ffffffffffff')。型は ktfnull / ktfother / ktfoper / ktfref / ktflist (= ktfobj) / ktfpat / ktfcx / ktfstring / ktfsymbol。要素・ビームラインは FFS 側データ構造で別管理 (要素コードは module maccode)。pysad は Python のネイティブ型で表すため不要。

### 3.3 32 bit 制約

ヒープ全体が 32 bit 前提 (MAP_32BIT) ということはない。アドレスは integer*8 (有効 49 bit)。残る 32 bit 制約は k64.readme の通りスタック索引とリスト/文字列長 (< 2**31)。`sim/unix_memory_.c` は既定 SAD_OFFSET_T=integer4 (MAP_32BIT は TRY_SAD_MAP_ADDR32 定義時のみ)、`unix_memory8_.c` は integer8 版の 22 行ラッパで、どちらも OBJSIM でビルドされ使用される。

### 3.4 ビルド (mk/ + config/)

- コンパイラフラグは `config/GCC.spec`: `-fno-second-underscore -fdollar-ok -fargument-alias -falign-functions=8 -fno-range-check` など。
- `mk/sad.obj.mk`: OBJMOD (tfstk, temit, tffs など 13 本)、OBJF、OBJRC、OBJAUTO (tturn, tturne, tbend, tquad, qtwiss, qcell, tcod, tftrack, track 等 127 本)、OBJ0〜6、OBJSIM、OBJGLUE。**OBJOBS (mbmp, pqcell, ttdr, pmeas 等 66 本) は定義のみでビルドされない**。
- `sad.conf`: USE_X11 / USE_TCLTK は YES、USE_EPICS は既定コメントアウト。`PySad.c` / `tfPyInter_.c` / `PyConfig.c` (Python 2 系 API) は mk/ に載らず既定ビルドに入らない。pysad はこのブリッジを使わない (Python 側が SAD を呼ぶ構成ではなく、SAD の計算を Python に直訳している)。

### 3.5 インタープリタチェーン

```
tfeval (tfetok でトークン化 → tfestk で式木) → tfeevalref (tfeeval) → tfefun (組込) → tfefun1 (ID 1001〜: FFS 系)
  → tftrack→tturn / temit / tfdapert / tffs ...
```

主要ファイル行数: tfstk.f 5,222、tfetok.f 640、tfeval.f 428、tfestk.f 975、tfeeval.f 1,059、tfefun.f 4,451、tfefun1.f 170、tfinitn.f 993、tffs.f 3,191、tffsa.f 1,838、tffscalc.f 629。`pqcell.f` は mbmp.f からのみ呼ばれ OBJOBS なので非ビルド。FFS の展開名は未確認 (「Flexible Function System」は根拠なし)。

---

## 4. 主要コマンド・関数と実装、pysad の対応

| コマンド / 関数 | SAD 実装 (呼出チェーン) | pysad 対応 | 照合桁 |
|---|---|---|---|
| `TrackParticles[]` | tfinitn:917 (ID 1005) → tfefun1:70 → tftrack:8 → tturn:36 → tturn1:153 | `ffscmd.FFS.fn_track` → `track6n.track_particles_fast` (numba) / `track6.track_particles` | ≤5e-11 (Level 6〜11)、FLUC の乱数状態一致 |
| FFS CALC (4D) | tffsa:544 → tffscalc:109 → qcell1 (qcell:27) → qcod/qtwiss1 | `lattice.Lattice.twiss`, `_qcod_periodic`, `_one_turn_cod_matrix` | 閉軌道 9-10 桁、Twiss ≤1e-10 (誤差・オフモメンタム付きは ≤3e-7、10_benchmarks §3.1) |
| FFS CALC6D | tffscalc → qcell1 → qcell6d (qcell.f:404) | `calc6d.calc6d` (`tfetwiss`, `etwiss2ri`, `_closed_orbit`) | 1.7e-10 (Level 18a/b/c) |
| FFS EMIT / `Emittance[]` | tffsa:988 → temitf → temit:78 → tcod → tturne → tecalc | `emit6d.emittance6d` / `tturne` / `tecalc`、規則一式は `emit_extra.emittance_rules` | 全量 ≤1.3e-11 (Level 16B)、OneTurn/Matrix 系 2.1e-10 (Level 21a) |
| FFS RADINT | tffsa:617 → INTGRL:11 | `radint.py` (直訳)。`lattice._radiation_integrals_bend` は BEND のみの独自実装 | 7 桁 (Level 29a/29a2) |
| FFS CHRO / 色収差 | tffsa:609 → tfchro.f | `CHROMATICITY` コマンドは対応 (18 の C-03)。色収差の値はチューンの DP0 スキャン + `PolynomialFit` (Fit.n) でも求める | チューン ≤2.5e-7 (Level 3A, 14A) |
| FFS GO (FIT) | tffsa:1207 → tffsmatch:15 (tfsolv, tgrad, fmincube) | `match6.Matcher.go` (数値微分)、`match6d.AnalyticMatcher` (tffsqu/qdcell 解析微分) | 反復数・残差・K1 一致 (Level 13)。解析版は SAD の反復列と K1 を 12 桁 |
| `CalculateOptics[]` | tfoptics.f (qcell1 の薄いラッパー) | `ffscmd.fn_calculate_optics` (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) | — |
| `DynamicApertureSurvey` | dapert.n → tfdapert:1 → trackd:37 → tturn、損失は tapert | `dapert.dynamic_aperture` (trackd/tapert の独自実装)、追跡は track6n | NORAD は地図・Score 完全一致 (Level 16A)。RAD+FLUC は損失周の乱数が未一致 |
| `TouschekLifetime` | emit.n → temit (INTRA ON: tintrb) の TouschekTable | `touschek.touschek_table` / `touschek_lifetime` | ≤7e-12 (Level 17a)、RADCOD ON は 1.7e-8 (17b) |
| ジオメトリ (GX/GY/GZ/CHI) | tfgeo:1 (tfgeo/tfgeo1) / tgrot.f (geolib) | `geo.compute_geo` | 4.8e-13 (Level 19a/b/c) |
| MAIN TRACK | ActTra:257 → track.f:1 | `maintrack.py` (TRACK コマンドの手順を再現。PHSPLOTS の描画は無視、FOURIE ON は未対応) | SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV のため、SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 |

---

## 5. ファイル分類表 (pysad 列つき)

### 5.1 粒子追跡 (非線形) 🔴A

| ファイル | ルーチン | 機能 | pysad: 関数 / 照合桁 |
|---|---|---|---|
| tftrack.f | tftrack:8 | TrackParticles 本体 | `track6.track_particles`, `track6n.track_particles_fast` / ≤5e-11 |
| tturn.f | tturn:36, tturn1:153 | 要素ごとの分岐 (select case L286-552) | `track6._tturn1`, `track6n.turn` / ≤5e-11 |
| track.f | track:1, twig:332 | MAIN TRACK、twig は空スタブ | `maintrack.py` (TRACK コマンドの手順を再現。SAD K64 の track.f 自体は tinip で SIGSEGV) |
| tbend.f | tbend:531, tbendcore, tbshift, tbrot | BEND 追跡 | `track6._tbend`, `_tbendcore`, `_tbshift`, `_tbrot` (n 版 `tbendcore1` ほか) / ≤5e-11 (Level 6〜9)、K1≠0 は `_tbendi`/`k_bendi` ≤1.7e-10 (z のみ) |
| tquad.f | tquad, tthin:118, ttfrin | QUAD / SEXT 以上の薄い多極 | `track6._tquad`, `_tthin_gen`, `_ttfrin*` (n 版 `k_quad`, `k_thin`) / ≤5e-11 |
| tsolqu.f | tsolqu, tsolqum, tsolqur | bz 付き厚肉四極 | `track6._tsolqu`, `_tsolqur_ens`, `_tsolqu_bz_ens` / ≤1.5e-10 (Level 9)。n 版は `track6n.tsolqu1`、`track6n_ext.tsolqu_bz_all` (track6 と 4e-17) |
| tsol.f | tsol:510 | ソレノイド領域 | `track6._tsol_region` (DRIFT/QUAD を含む領域のみ) / ≤1.5e-10。n 版 `track6n_ext.k_sol_*` / 4e-17 |
| tcav.f | tcav, tcavfrin | CAVI | `track6._tcav`, `_tcavfrin` (n 版 `k_cav`) / ≤5e-11 (長さ 0 と L≠0、Level 8) |
| tmulti.f | tmulti:37, tmulk | MULT | `track6._tmulti_ens` (ANGLE=0、K0/K1/K2+) / ≤5.9e-11 (Level 11)。n 版 `track6n_ext.k_mult` / 4e-17 |
| tdrift.f | tdrift_free, tdrift_solenoid | ドリフト | `track6._tdrift_free`, `_tdrift_solenoid` / ≤5e-11 |
| photons.f | tradk:1267, tradke:1362 | 放射 (追跡 / 行列側) | `track6._tradk`, `track6n.tradk_core1`, `sadrandom.drndsr/tdusrn` / ≤4.9e-11、FLUC の seed 一致 (Level 10) |
| ttcav.f | ttcav | TCAVI | **未実装** |
| tcoord.f | tcoord | COORD | **未実装** (CALC/EMIT 側 tcoorde も未) |
| temap.f / tins.f | temap, tins | MAP (ExternalMap) / INS (線形挿入写像) | **未実装** (INS の CALC は `lattice.ins_twiss` で別途対応、Level 5) |
| txwake.f / bbstrhl.f / undulator.f / phsrot.f / tspch_.c | — | WAKE / BEAMBEAM / UND / PHSROT / SPCH | **未実装** |
| tapert.f | tapert | 損失判定 (矩形 LOSSAMPL/LOSSDZ) | `track6n.k_apert`、`dapert` / 地図一致 |

### 5.2 線形光学・COD 🔴A

| ファイル | ルーチン | 機能 | pysad: 関数 / 照合桁 |
|---|---|---|---|
| qtwiss.f | qtwiss1:14, qcod:810, qquad:1737 | CALC の Twiss・4D COD | `lattice._qcod_periodic` (qcod 直訳)、`_transport_twiss` (qmat2twiss は独自、非結合の標準式) / 9-10 桁 |
| qcell.f | qcell1:27, qcell6d:404 | 周期解・チューン | `lattice.Lattice.twiss` (A/B/分散は直訳、qmdiag は numpy 固有値分解)、`calc6d.calc6d` / ≤1e-10、1.7e-10 |
| tquade.f / tsolque.f | tquade, tsolque:313 | QUAD の写像 (bz 付き) | `lattice._tquade_cod`, `_tsolque_cod`, `_tqfrie_cod`, `_tqlfre_cod`, `_tqente_cod`、`emit6d._tquade` / ≤1e-11 (Level 5, 8) |
| tbende.f / tbfrie.f | tbende:351, qbend:627, tbedge:88 | BEND 行列 | `lattice._bend_cod_update`, `_tbedge_cod`, `_tblfre_cod`、`emit6d._tbende`, `_tbendecorr`, `_tbendebody` / ≤1e-10 (K1≠0 は ≤3e-9)。ANGLE=0 の F1 分岐は未実装 |
| tthine.f | tthine | SEXT 等の行列 | `lattice._tthine_cod`, `_thin_offset_map` / ≤1e-10 |
| tmulte.f | tmulte | MULT 行列 | `lattice._mult_sequence`, `emit6d._tmulte` / Level 8/11 一致。ANGLE≠0 (tmultae) 未実装 |
| tcave.f | tcave, tcavfrie | CAVI 行列 | `emit6d._tcave`, `_tcavfrie_map` / Level 8 一致 |
| tsole.f / ttcave.f / tcoorde.f / tinse.f | — | SOL 行列 / TCAVI / COORD / INS | SOL・TCAVI・COORD は **未実装**。INS は `lattice.ins_twiss` |
| tturne.f | tturne:1, tturne1:334 | 6D 一周 (行列 + ビーム行列) | `emit6d.tturne` / ≤1.3e-11、`emit_extra._turn_record` |
| tcod.f | tcod:31 | 6D COD (Newton) | `calc6d._closed_orbit`、`emit6d.emittance6d` 内の割線法 (RADCOD OFF は tcodz 相当) / δ, trf0 一致 (Level 7) |
| tffsqu / qdcell.f | tffsqu, qdcell, qdtwis, qdquad, qdbend, qdmdia | 解析微分 | `match6d.tffsqu`, `qdcell`, `qdtwis`, `qdquad`, `qdbend`, `qdmdia` / 数値微分と 2e-8〜4e-7 |
| tdfun (qdcell.f) | tdfun | 残差 df | `match6._tdfun1`, `_trdfun` / 一致 |

### 5.3 エミッタンス・放射積分・ほか 🔴A

| ファイル | ルーチン | 機能 | pysad: 関数 / 照合桁 |
|---|---|---|---|
| temit.f | temit:78, tecalc:235 | 6D 平衡ビーム行列 (teigenc + tsolva/tsvd) | `emit6d.emittance6d`, `tecalc`, `_tsolva`, `_normal_form` / ≤1.3e-11 |
| photons.f | tradke:1362 | 行列側の放射 | `emit6d` 内 `_tquade` 等 (krad) / 同上 |
| INTGRL.f | INTGRL:11 | 放射積分 | `radint.py` (直訳、FFS RADINT に接続)、`lattice._radiation_integrals_bend`, `_symps` (BEND のみの別実装) / 7 桁 (Level 29a)。RF 関係 (INTGRL.f:249-273) は未実装 |
| tfdapert.f / trackd.f | tfdapert, trackd | DA | `dapert.dynamic_aperture` / 地図一致 |
| tffscalc.f / tffsmatch.f | tffscalc:12, tffsmatch:15 | 残差 / GO | `match6.Matcher` / 反復一致 |
| tfgeo.f / tgrot.f | tfgeo, tfgeo1, tfchi.. | ジオメトリ | `geo.compute_geo`, `tfchi`, `tfrotgeo`, `tgrot`, `tforbitgeo` / 4.8e-13 |
| tintrb (temit 経由) | touckl, touckm, toucke | IBS・Touschek | `touschek._tintrb`, `touschek_table` / ≤7e-12 |
| tfchro.f | — | CHRO | `CHROMATICITY` コマンドとして対応 (薄肉近似の表) |
| tfoptics.f | — | CalculateOptics | `ffscmd.fn_calculate_optics` (idim=3 は全周) |

### 5.4 数学 🟡B

| ファイル | ルーチン | 機能 | pysad |
|---|---|---|---|
| tsolvg.f | tsolvg, tsolva, tsvd | 最小二乗・SVD | `numpy.linalg` / `emit6d._tsolva` (SVD 細部の対応は未確認) |
| teigen.f | teigenc ほか | 固有値 | numpy (qmdiag も固有値分解) |
| matrix.f / tmultr.f / tinitr.f | 6×6 行列 | numpy (`match6d.tmultr45`) |
| spline.f | spline, splint | スプライン | `touschek.SadSpline`, `_spline_ddy` |
| gaussn.f / drndsr.f / random_plugin_sad.c | 乱数 | `sadrandom.py` (LCG、Box-Muller、drndsr、tdusrn) / seed 一致 |
| tffindroot.f | FindRoot / Fit | `sadscript.py` の Fit.n 部 (PolynomialFit) |

### 5.5 インタープリタ 🔵C / グラフィックス 🟢D / C ソース・sim・inc・Packages・extensions・contrib

要約:

- 🔵C: tfstk/tfetok/tfeval/tfestk/tfeeval/tfefun/tfefun1/tfinitn/tffs/tffsa/tffscalc/toplvl/terror/initb1。pysad は `sadscript.py` (字句・Pratt 構文・評価、Module/Do/Table/パターン/ファイル I/O、autofg 書式)、`sadlist.py` (MapThread, Cases, MatchQ 等のリスト系組み込み)、`parser.py` (MAIN 部)、`ffscmd.py` (tffsa の語の解釈順、tflag/terror/tftype) で必要分のみ再実装。
- 🟢D: tfTkInter_.c、tfXlib_.c、sad_tcltk.c、sad_xlib.c、geodrw.f、Packages の Plots 系。**pysad では不要** (描画は matplotlib: `plot_optics.py`)。
- C: calc.c は SAD1 の数式用 (byacc 生成)、PySad.c は Python 2 API で非ビルド、tspch_.c/spch.c は空間電荷、tgauss_.c/random_driver.c/random_plugin_sad.c は乱数、sim/ は OS 抽象化 (unix_memory 系、sad_api、dynl 系ほか)。
- inc/: 実コードを含むのは TENT/TEXIT/common/MACexpn/TFMACRO/TMACRO/TFSF/coroper 程度。要素コードの実定義は tfstk.f の module maccode (icDRFT=1, icBEND=2, icQUAD=4, icSEXT=6, icOCTU=8, icDECA=10, icDODECA=12, icUND=18, icWIG=19, icSOL=20, icST=21, icMULT=22, icCAVI=31, icTCAV=32, icMAP=33, icINS=34, icCOORD=35, icBEAM=36, icPROT=37, icSPCH=38, icMARK=41, icMONI=42, icAPRT=43, icLINE=100)。pysad の `ffscmd._TYPECODE` が SOL=20, MULT=22 を使う。
- Packages: init.n (NPARA=4 など既定値)、Optics.n (CALC6D, NormalCoordinates)、emit.n (Emittance[] の規則名、TouschekLifetime)、dapert.n、Fit.n、correction.n、optimize.n、extmap.n。pysad が参照しているのは init.n / emit.n / Fit.n / dapert.n / Optics.n の挙動。
- extensions: RandomMT (MT19937/SFMT、既定の乱数は "SAD" 合同法)、sad+ (大見氏の微分代数 C++)、Standard (Tkinter/Xlib)。pysad は使わない。
- contrib: gdtoa、BBBrem_dir (BBBrem.C)、libtai、tcl/tk ほか。不要。

---

## 6. コールスタックと pysad の関数名

### 6.1 TrackParticles

| SAD | pysad (純 Python: track6.py / numba: track6n.py) |
|---|---|
| `TrackParticles[...]` (tfinitn:917, ID 1005) | `ffscmd.FFS.fn_track` → `track6n.track_particles_fast` (MULT/SOL は `track6n_ext.track_particles_ext`) / `track6.track_particles` |
| tfefun1:70 → tftrack:8 | `track_particles` (粒子配列準備、RNG は `track6.RNG`) |
| tturn:36 → tturn1:153 | `track6._tturn1` / `track6n.turn` (要素外側・粒子内側) |
| 要素 select case | `track6._step_el` / `track6n.Prog` の kind 表 |
| 放射 tradk (photons) | `track6._tradk` → `_tradk_core`、`track6n.tradk_core1` |
| tapert (損失判定) | `track6n.k_apert` |

### 6.2 Emittance

| SAD | pysad |
|---|---|
| `Emittance[]` (emit.n) → temitf → temit:78 | `emit_extra.emittance_rules` → `emit6d.emittance6d` |
| tcod:31 (tturne の Newton 閉軌道) | `emittance6d` 内 (割線法/6D Newton)、`calc6d._closed_orbit` |
| tturne:1 / tturne1:334 | `emit6d.tturne` (`_tdrife0`, `_tquade`, `_tbende`, `_tcave`, `_tmulte`、tthine は lattice の関数) |
| tradke | 各 `_t***e` の krad 引数 |
| tecalc:235 → teigenc / tsolva | `emit6d.tecalc` / `_tsolva` / `_normal_form` |
| INTRA ON: tintrb | `touschek._tintrb`, `_tturne_intra` |

### 6.3 CALC (CAL)

| SAD | pysad |
|---|---|
| tffsa:544 → tffscalc:12 | `ffscmd` の CALC → `match6.Matcher._tffscalc` (残差) / `Lattice.twiss` |
| qcell1:27 | `Lattice.twiss`, `twiss_table` |
| qcod:810 | `lattice._qcod_periodic`, `_one_turn_cod_matrix` |
| qtwiss1:14 → qquad / qbend / qthin | 要素ごとの COD 写像 (`_tquade_cod`, `_bend_cod_update`, `_tthine_cod`) の数値微分ヤコビアン (`_num_jacobian`、eps=1e-7) |
| qmdiag | `Lattice.twiss` 内 (numpy 固有値分解) |
| qcell6d:404 (CALC6D) | `calc6d.calc6d` |

### 6.4 GO (FIT)

| SAD | pysad |
|---|---|
| tffsa:1207 → tffsmatch:15 | `match6.Matcher.go` |
| tfsolv / tgrad / fmincube | `Matcher` 内、`match6._fmincube`、tdfun は `_tdfun1`、重みは tweigh |
| tffsqu → qdcell → qdtwis → qdquad/qdbend/qdmdia | `match6d.tffsqu` → `qdcell` → `qdtwis` → `qdquad`/`qdbend`/`qdmdia` (`AnalyticMatcher`) |

### 6.5 DynamicApertureSurvey

| SAD | pysad |
|---|---|
| dapert.n → tfdapert:1 | `ffscmd.FFS.fn_dasurvey` → `dapert.dynamic_aperture` |
| trackd:37 (走査・境界探索・Score) | `dapert.dynamic_aperture` (独自実装。スロット順・DAPWIDTH 窓・62 進地図 `rad62`) |
| tturn (1 周) | `track6n.turn` |
| tapert:1 (損失) | `track6n.k_apert` |

---

## 7. tturn / tturne ディスパッチ表 (pysad 列)

### 7.1 tturn1 (tturn.f select case L286-552)

| コード | case | 呼出 | pysad |
|---|---|---|---|
| 1 | icDRFT | tdrift_free:294 | ✓ `_tdrift_free` (≤5e-11) |
| 2 | icBEND | tbend:327 | ✓ `_tbend` (K1≠0 `_tbendi`)、DX/DY/DROTATE も可 |
| 4 | icQUAD | tquad:355 → tsolqu | ✓ `_tquad` / `k_quad` |
| 6,8,10,12 | SEXT/OCTU/DECA/DODECA | tthin:380 | ✓ `_tthin_gen` / `k_thin` (SEXT 厚肉放射も一致) |
| 22 | icMULT | tmulti:395 | △ `_tmulti_ens` / `k_mult` (ANGLE=0、SK/F1/RF は未) |
| 20 | icSOL | tsol:410 | △ `_tsol_region` / `k_sol_*` (DRIFT/QUAD を含む領域のみ。内部 BEND/MULT/CAVI、DX/DY/DZ/CHI、F1 は未) |
| 31 | icCAVI | tcav:443 | ✓ `_tcav` / `k_cav` |
| 41,42 | MARK, MONI | — | ✓ (無作用) |
| 43 | icAPRT | tapert1:418 | △ `k_apert` (損失判定のみ) |
| 32 | icTCAV | ttcav:480 | ✗ |
| 33 | icMAP | temap:490 | ✗ |
| 34 | icINS | tins:502 | ✗ |
| 35 | icCOORD | tcoord:506 | ✗ |
| 36 | icBEAM | beambeam:498 | ✗ |
| 37 | icPROT | phsrot:517 | ✗ |
| 38 | icSPCH | tpspac:525 | ✗ |
| 18 | icUND | undulator:537 | ✗ |
| 19 | icWIG | twig:541 (空スタブ) | — (SAD も無作用) |
| 21 | icST | abort ("Use BEND with ANGLE=0") | — |
| — | wake / wspac | txwake:266 / twspac:274 | ✗ |

### 7.2 tturne1 (tturne.f select case L470-675)

| コード | case | 呼出 | pysad |
|---|---|---|---|
| 1 | icDRFT | tdrife0:475 | ✓ `emit6d._tdrife0` |
| 2 | icBEND | tbende:533 | ✓ `_tbende` (`_tbendecorr`, `_tbendebody`) |
| 4 | icQUAD | tquade:569 | ✓ `_tquade` |
| 6,8,10,12 | tthine:592 | tthine | ✓ (六極内のドリフトの dv 項を含む) |
| 20 | icSOL | tsole:597 | ✗ (emit6d は SOL 未対応) |
| 21 | icST | forcesf:602 | — |
| 22 | icMULT | tmulte:621 | △ `_tmulte` (ANGLE=0。ANGLE≠0 の tmultae は未) |
| 31 | icCAVI | tcave:635 | ✓ `_tcave`, `_tcavfrie_map` |
| 32 | icTCAV | ttcave:648 | ✗ |
| 33 | icMAP | temape:661 / qemap:657 | ✗ |
| 34 | icINS | tinse:665 | ✗ (CALC の INS は別実装) |
| 35 | icCOORD | tcoorde:668 | ✗ |
| — | default | 無作用 (WIG, UND, MARK, BEAMBEAM) | ✓ |

---

## 8. デッドコード

呼出元なし AND (Obsolete コメント OR 非ビルド) のもの: `tquads.f`, `tquase.f`, `tbrad.f` (tlinit の旧コピー)、`tmulbs.f` (ファイルのみ dead、関数 tmulbs は tffs.f:2761 に移設され現役)、`ttdr.f`、`trad.f` (現役は photons.f の tradk)、`tserad.f`, `tstrad.f`, `tqrad.f`, `trade.f` (現役は tradke)。`tmuld.f` は一部のみ dead (tmuld6 は photons.f:1494、temap.f:432 から現役)。`tbfrie.f` は現役 (tbedge が tbende/tsteee から)。OBJOBS 66 本ほか非ビルドの .f は個別判定未実施。pysad は dead を一切移植していない。

---

## 9. Python 化優先度

### 9.1 高 (コア): すべて完了

| 機能 | 状態 | pysad |
|---|---|---|
| 追跡: DRIFT/BEND/QUAD/SEXT/CAVI/MARK | ✓ | track6/track6n |
| 追跡: MULT (ANGLE=0)、SOL 領域 | △ 部分 | track6 / track6n_ext |
| qtwiss/qcell/qcod (CALC) | ✓ | lattice.py |
| tturne/temit/tecalc (Emittance[]) | ✓ (DRIFT/QUAD/BEND/SEXT/CAVI/MULT ANGLE=0/MARK) | emit6d.py, emit_extra.py |
| tffsmatch/qdcell (GO) | ✓ 数値微分 + 解析微分 | match6.py, match6d.py |
| INTGRL (RADINT) | ✓ 直訳 | radint.py |
| qcell6d (CALC6D) | ✓ | calc6d.py |
| trackd/DA | ✓ (独自実装、NORAD 一致) | dapert.py |
| tfgeo | ✓ | geo.py |
| touschek (INTRA) | ✓ | touschek.py |
| tfchro (CHROMATICITY)、tfoptics (CalculateOptics) | ✓ | ffscmd.py |
| MAIN の TRACK コマンド | ✓ (SADScript 版と照合) | maintrack.py |

### 9.2 中 (将来)

| 項目 | 状態 | 備考 |
|---|---|---|
| ttcav/ttcave, tcoord/tcoorde, temap, tins (追跡) | ✗ | 該当要素を使う格子で必要 |
| emit6d の SOL、tmultae (MULT ANGLE≠0)、MULT の SK/F1 | ✗ | |
| BEND の F1 フリンジ ANGLE=0 分岐、INTGRL の RF 関係 | ✗ | |
| qtwissfrac (端数位置の Twiss) の厳密化 | △ | 端数位置は補間 |
| CALC の解析行列の残り (MULT、ANGLE=0 BEND、F1 付き BEND、DROTATE≠0 BEND、空洞) | △ | 数値微分 (eps=1e-7)。四極・薄い多極・BEND 本体は解析行列 |
| 損失周の乱数消費 (RAD+FLUC の DA) | △ | 論理は再同期テストで SAD と一致。実寸 DAsv の最終 seed はカオスで一致しない |
| 追跡の並列粒子 fork (NPARA) | ✗ | NPARA=1 前提で一致 |

### 9.3 低・不要

🔵C インタープリタ (最小限を sadscript で自前実装)、🟢D グラフィックス (matplotlib)、⚫E ビルド、🟤G contrib、デッドコード 10 件。

### 9.4 旧 API (非推奨)

`ffs.py` / `matching.py` / `tracking.py` / `optics.py` は SAD 未照合の独自実装 (tracking.py の放射係数に疑い)。`__init__.py` の export は旧 API のまま。DA は `matching.py:507` の旧 DA ではなく `dapert.py` が現行。

---

## 10. pysad 現状

- 到達点: ベンチ bench_level1〜31 (run_regression.py の 50 ベンチ) が SAD と一致。線形光学 (CELL/INS/CALC6D/幾何/RADINT/TRPT の加速) ≤1e-10、
  Emittance ≤1e-11〜1e-9、Touschek ≤7e-12、追跡 ≤5e-11、マッチングは反復列と K1 を SAD と同じ桁で再現、DA の地図と Score (NORAD) が一致。
- 現行モジュール (`pysad/`、行数) は `15_implementation_summary.md` §2 (lattice 4038、ffscmd 3568、sadscript 2831、emit6d 2265、track6 2054、parser 1876、track6n 1826、
  elkern 1723、match6d 1468、match6 1421、sadlist 1167、touschek 849、track6n_ext 696、calc6d 660、radint 630、geo 556、trptacc 334、optics_fn 315 ほか)。
  calc6d / match6d / geo / emit_extra / track6n_ext / touschek はすべて ffscmd から呼ばれる。
- 実スクリプト: esrf 52 本 (数値一致 26、一部一致 13、例外 6、打ち切り 7)、pf2 34 本 (数値一致 22、一部一致 10、動くが不一致 2)。詳細は `12_user_scripts_status.md`。

## 11. 未実装・未解決の解説

1. **性能**: 誤差入り 6D emit × 多数の実現 (print_all_error)、DAsv 5000 周 × 969 粒子、CodCorr の全実現は 600 s を超える (SAD も DA 系は 300 s 超)。
2. **FITP>1 の最終解**: SAD (WSL) が GO で再起動するため参照値がない。1 段目の解は一致、2 段目以降の多解は別解に着くことがある。fam 点は自己検証のみ。
3. **FFS エラー時の Do の中断**: SAD は FFS::error で外側の Do を抜けるが pysad は続行する (chrom_b, disp_b, dx_sv)。
4. **不安定ラティス**: SAD の 1e199 級の値や NaN は再現しない (安定域は一致)。
5. **追跡**: 損失周の乱数消費は論理が一致するが、実寸 DAsv の最終 seed はカオスで一致しない (Score の差 1〜3 %)。
6. **要素未対応**: TCAVI、COORD、MAP、INS (追跡)、BEAMBEAM、PHSROT、SPCH、UND、WAKE。SOL の emit/幾何、MULT ANGLE≠0 (tmulta)、加速要素の K_n と VOLT の併用。
7. **解析的でない部分**: fit の幾何条件の微分 (tdgeo)、MULT・ANGLE=0 BEND・F1 付き BEND・空洞の CALC 行列 (数値微分)。
8. **MON/MCOD/TCOD**: この SAD ビルドでは何もしない。pysad も受理のみで一致。
9. 解析残り (「転記のみ・導出なし」「未検証」の一覧) は `docs/16_physics_guide.md` §22 を参照。RADINT と Emittance[] の 13% 差の主因 (16 の §11) は未特定。
