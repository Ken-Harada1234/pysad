# 32. SAD ソースの構造とコールスタック

この文書で分かること: K64 版 SAD のソースの構成 (フォルダ、6 層のアーキテクチャ、ビルド)、主要コマンドが Fortran のどのルーチンを通って実行されるか、pysad のどの関数がそれに対応するか、使われていないコード (デッドコード)。
前提: SAD (Strategic Accelerator Design。KEK 発の加速器光学・追跡コード) の基本的な使い方と、pysad の概要 (docs/01)。
読む順: §1 概要、§3 アーキテクチャ、§4 コマンドと実装の対応表、§5 コールスタックの順。§2 と §6 は必要なときに引く。

対象 SAD は K64 版 (https://github.com/KatsOide/SAD)。K64 は、このリポジトリが翻訳元にしている SAD のソース版である。コマンドの意味は CERN SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html) を先に見る。数値の根拠はローカルの `<SAD>\src`。ソースの行番号は 2026-10 時点、ファイル数・行数は 2026-11-06 に数えた値。

用語:
- **SADScript**: SAD の Mathematica 風のスクリプト言語 (拡張子 .n)。
- **FFS**: SAD のラティスに対して CALC (光学計算)・FIT (マッチング)・EMIT (エミッタンス)・USE などのコマンドを実行する部分。コマンド分岐は `tffsa.f`、本体は `tffs.f` と `tffscalc.f`。FFS の語の展開名は未確認 (「Flexible Function System」は根拠なし)。
- **pysad の「直訳」**: Fortran のルーチンを、構造・手順を変えずに Python に書き写したもの。

未対応・状況は docs/31、ファイル別の翻訳状況は docs/33 にある。この文書は構造とルーチンの対応だけを扱う。

---

## 1. 概要

SAD は、固定形式 Fortran (module / select case などの F90 構文を多用) と一部の C/C++ をコアに、Mathematica 風の SADScript を持つ。次の表は、ソースの規模と pysad の規模を並べたものである。

| 項目 | 実測 |
|---|---|
| `src/*.f` | 474 本、計 141,678 行 (`cat src/*.f \| wc -l`) |
| `src/*.c` / `src/sim/` | 38 本 / 42 ファイル (.c 20、.h 16 ほか) |
| `src/inc/` | 35 エントリ (.inc 33 + changelog + k64.readme)。大半は全行コメントアウト |
| `Packages/` | 82 エントリ (.n 73、.py 6 ほか) |
| `mk/` | .mk 19 本 + Makefile 2 本 |
| pysad (現行系統) | `pysad/*.py` 71 本、計 57,124 行 (2026-11-17 に数えた。非推奨の旧 API ffs/matching/tracking/optics の 3,911 行を含む)。試験の件数は docs/27 §5.1 の最新の通し (GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip) |

ファイルの分類 (この記号は docs/33 の表でも使う): 🔴A 物理コア / 🟡B 数学 / 🔵C インタープリタ / 🟢D グラフィックス / ⚫E ビルド / 🟣F SADScript パッケージ / 🟤G 外部ライブラリ。pysad の翻訳の対象は A と B。C は `sadscript.py` + `ffscmd.py` で最小限を自前実装し、D・E・G は不要としている。

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
├ extensions/ RandomMT/  sad+/ (DA (微分代数), C++ 33 本)  Standard/ (Tkinter, Xlib)
├ contrib/ gdtoa BBBrem_dir libtai-0.60 tcl8.6.11 tk8.5.18 tk8.6.11 btkshim patches tools
├ mk/ (sad.obj.mk ほか)  config/ (GCC.spec)  sad.conf  Makefile
```

## 3. アーキテクチャ

### 3.1 6 層

SAD は 6 層でできており、pysad は 1〜5 層を Python のモジュールに対応させ、6 層は不要としている。

| 層 | 内容 | pysad 側 |
|---|---|---|
| 1 SADScript (.n/.sad) | ユーザスクリプト・格子 | `sadscript.py` |
| 2 パーサ | tfetok (字句) + tfeval/tfestk (演算子順位)。式木はリスト。**バイトコードやスタックマシンではない**。calc.y は無関係 | `sadscript.py` (Pratt 構文) |
| 3 評価器 | tfeeval (式木の再帰評価)、tfefun/tfefun1 (組込関数、id>1000 は tfefun1 の computed goto) | `sadscript.py` / `ffscmd.py` の組込表、`sadlist.py` (リスト系組み込み) |
| 4 FFS | tffs / tffsa (コマンド分岐) / tffscalc | `ffscmd.py` |
| 5 物理コア Fortran | tturn/tturne/qtwiss/temit/... | docs/33 の各モジュール |
| 6 C OS 層 | sim/ (メモリ、シグナル、動的リンク、X11/Tk 薄グルー) | 不要 |

### 3.2 NaN タギング

全 SADScript オブジェクトは 64 bit ディスクリプタ。実数以外は IEEE double の NaN 空間に型とアドレスを埋める (`tfstk.f:822-823`: ktfmask = z'fffe000000000000'、ktamask = z'0001ffffffffffff')。型は ktfnull / ktfother / ktfoper / ktfref / ktflist (= ktfobj) / ktfpat / ktfcx / ktfstring / ktfsymbol。要素・ビームラインは FFS 側データ構造で別管理 (要素コードは module maccode)。pysad は Python のネイティブ型で表すため、この仕組みは不要。

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

主要ファイル行数: tfstk.f 5,222、tfetok.f 640、tfeval.f 428、tfestk.f 975、tfeeval.f 1,059、tfefun.f 4,451、tfefun1.f 170、tfinitn.f 993、tffs.f 3,191、tffsa.f 1,838、tffscalc.f 629。`pqcell.f` は mbmp.f からのみ呼ばれ OBJOBS なので非ビルド。

---

## 4. 主要コマンド・関数と実装、pysad の対応

次の表は、主要なコマンドが SAD でどのルーチンを通って実行され (呼出チェーン。`ファイル:行`)、pysad のどの関数が対応するかを示す。右端の列は、その対応の限界・注意 (gs との照合の桁は [history/32_source_analysis_validation.md](history/32_source_analysis_validation.md))。

| コマンド / 関数 | SAD 実装 (呼出チェーン) | pysad 対応 | 限界・注意 |
|---|---|---|---|
| `TrackParticles[]` | tfinitn:917 (ID 1005) → tfefun1:70 → tftrack:8 → tturn:36 → tturn1:153 | `ffscmd.FFS.fn_track` → `track6n.track_particles_fast` (numba) / `track6.track_particles`。POL かつ RAD ON は純 Python の `track6` (`sprot`・`bsi`)、RADLIGHT は `ffscmd._RadLight` (tlinit / tlstore / tlresult)、MULT ANGLE≠0 は `multa.py`、COORD・TCAVI・UND・BEAMBEAM・WAKE・MAP はそれぞれ `coord.py`・`tcav.py`・`und.py`・`beambeam.py`・`wake.py`・`elements.Map` | FLUC の乱数状態は SAD と一致する |
| FFS CALC (4D) | tffsa:544 → tffscalc:109 → qcell1 (qcell:27) → qcod/qtwiss1 | `lattice.Lattice.twiss`, `_qcod_periodic`, `_one_turn_cod_matrix` | 誤差・オフモメンタム付きは照合の桁が下がる (docs/27 §3.1) |
| FFS CALC6D | tffscalc → qcell1 → qcell6d (qcell.f:404) | `calc6d.calc6d` (`tfetwiss`, `etwiss2ri`, `_closed_orbit`) | |
| FFS EMIT / `Emittance[]` | tffsa:988 → temitf → temit:78 → tcod → tturne → tecalc | `emit6d.emittance6d` / `tturne` / `tecalc`、規則一式は `emit_extra.emittance_rules`。POL のスピンは `spin.py`、SOL 領域は `solemit.py`、COORD・TCAVI・UND・MULT ANGLE≠0・MAP は `coord.py`・`tcav.py`・`und.py`・`multae.py`・`elements.Map` | OneTurn/Matrix 系は照合の桁が別。POL の MULT・空洞・SOL は [`test_spin.py`](../pysad/tests/script_tests/test_spin.py) |
| FFS RADINT | tffsa:617 → INTGRL:11 | `radint.py` (直訳)。`lattice._radiation_integrals_bend` は BEND のみの独自実装 | |
| FFS CHRO / 色収差 | tffsa:609 → tfchro.f | `CHROMATICITY` コマンドに対応 (docs/33 §9)。色収差の値はチューンの DP0 スキャン + `PolynomialFit` (Fit.n) でも求める | |
| FFS GO (FIT) | tffsa:1207 → tffsmatch:15 (tfsolv, tgrad, fmincube) | `match6.Matcher.go` (数値微分)、`match6d.AnalyticMatcher` (tffsqu/qdcell 解析微分) | 反復数・残差・K1 は SAD と一致。解析版は SAD の反復列と K1 を再現 |
| `CalculateOptics[]` | tfoptics.f (qcell1 の薄いラッパー) | `ffscmd.fn_calculate_optics` (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) | |
| `DynamicApertureSurvey` | dapert.n → tfdapert:1 → trackd:37 → tturn、損失は tapert | `dapert.dynamic_aperture` (trackd/tapert の独自実装)、追跡は track6n | NORAD は地図・Score が完全一致。RAD+FLUC は損失周の乱数が未一致 |
| `TouschekLifetime` | emit.n → temit (INTRA ON: tintrb) の TouschekTable | `touschek.touschek_table` / `touschek_lifetime` | RADCOD ON は照合の桁が下がる |
| ジオメトリ (GX/GY/GZ/CHI) | tfgeo:1 (tfgeo/tfgeo1) / tgrot.f (geolib) | `geo.compute_geo` | |
| `SynchroBetaEmittance` | emit.n → tfemits / temits.f | `synchrobeta.py` | |
| `RadiationField` / `RadiationSpectrum` | tlinit.f の tlcalc ほか (RADLIGHT の軌道が入力) | `radfield.py` | RadiationField は gs が SIGSEGV (docs/30 B49) なので、等速運動の場との比較のみ。軌道は gs と一致 |
| CSR 系 (`CSRTrack` `CSRHaissin` `CSRSetupOY` `CSRMatrix` `CSRInit` `CSRConvert`) | csrtrack.f、csroy.f、csrtest.f (csrhelmat) | `csr.py` | CSROYMatrix は gs が wrongtype (B60) で未照合 |
| `GaussianCoulomb(Fitted)` / BBBrem | tspac.f (bbkick1)、tintrb.f (twspu)、tbbbrem.f | `gausscoul.py`、`bbbrem.py` | docs/06・14 |
| DELCOR | tcorr.f、tvcorr.f | `ffscmd._delcor` | docs/04 |
| MAIN TRACK | ActTra:257 → track.f:1 | `maintrack.py` (TRACK コマンドの手順を再現。PHSPLOTS の描画は無視、FOURIE ON は未対応) | SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV のため、SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 |

---

## 5. コールスタックと pysad の関数名

次の 5 つの表は、コマンドごとに SAD のコールスタックを上から順に並べ、各段に対応する pysad の関数を示す。

### 5.1 TrackParticles

| SAD | pysad (純 Python: track6.py / numba: track6n.py) |
|---|---|
| `TrackParticles[...]` (tfinitn:917, ID 1005) | `ffscmd.FFS.fn_track` → `track6n.track_particles_fast` (MULT/SOL は `track6n_ext.track_particles_ext`) / `track6.track_particles` |
| tfefun1:70 → tftrack:8 | `track_particles` (粒子配列準備、RNG は `track6.RNG`) |
| tturn:36 → tturn1:153 | `track6._tturn1` / `track6n.turn` (要素外側・粒子内側) |
| 要素 select case | `track6._step_el` / `track6n.Prog` の kind 表 |
| 放射 tradk (photons) | `track6._tradk` → `_tradk_core`、`track6n.tradk_core1` |
| tapert (損失判定) | `track6n.k_apert_gen`・`k_apert` |

### 5.2 Emittance

| SAD | pysad |
|---|---|
| `Emittance[]` (emit.n) → temitf → temit:78 | `emit_extra.emittance_rules` → `emit6d.emittance6d` |
| tcod:31 (tturne の Newton 閉軌道) | `emittance6d` 内 (割線法/6D Newton)、`calc6d._closed_orbit` |
| tturne:1 / tturne1:334 | `emit6d.tturne` (`_tdrife0`, `_tquade`, `_tbende`, `_tcave`, `_tmulte`、tthine は lattice の関数) |
| tradke | 各 `_t***e` の krad 引数 |
| tecalc:235 → teigenc / tsolva | `emit6d.tecalc` / `_tsolva` / `_normal_form` |
| INTRA ON: tintrb | `touschek._tintrb`, `_tturne_intra` |

### 5.3 CALC (CAL)

| SAD | pysad |
|---|---|
| tffsa:544 → tffscalc:12 | `ffscmd` の CALC → `match6.Matcher._tffscalc` (残差) / `Lattice.twiss` |
| qcell1:27 | `Lattice.twiss`, `twiss_table` |
| qcod:810 | `lattice._qcod_periodic`, `_one_turn_cod_matrix` |
| qtwiss1:14 → qquad / qbend / qthin | 要素ごとの COD 写像 (`_tquade_cod`, `_bend_cod_update`, `_tthine_cod`) の数値微分ヤコビアン (`_num_jacobian`、eps=1e-7) |
| qmdiag | `Lattice.twiss` 内 (numpy 固有値分解) |
| qcell6d:404 (CALC6D) | `calc6d.calc6d` |

### 5.4 GO (FIT)

| SAD | pysad |
|---|---|
| tffsa:1207 → tffsmatch:15 | `match6.Matcher.go` |
| tfsolv / tgrad / fmincube | `Matcher` 内、`match6._fmincube`、tdfun は `_tdfun1`、重みは tweigh |
| tffsqu → qdcell → qdtwis → qdquad/qdbend/qdmdia | `match6d.tffsqu` → `qdcell` → `qdtwis` → `qdquad`/`qdbend`/`qdmdia` (`AnalyticMatcher`) |

### 5.5 DynamicApertureSurvey

| SAD | pysad |
|---|---|
| dapert.n → tfdapert:1 | `ffscmd.FFS.fn_dasurvey` → `dapert.dynamic_aperture` |
| trackd:37 (走査・境界探索・Score) | `dapert.dynamic_aperture` (独自実装。スロット順・DAPWIDTH 窓・62 進地図 `rad62`) |
| tturn (1 周) | `track6n.turn` |
| tapert:1 (損失) | `track6n.k_apert_gen`・`k_apert_el` |

---

## 6. デッドコード

呼出元なし AND (Obsolete コメント OR 非ビルド) のもの: `tquads.f`, `tquase.f`, `tbrad.f` (tlinit の旧コピー)、`tmulbs.f` (ファイルのみ dead、関数 tmulbs は tffs.f:2761 に移設され現役)、`ttdr.f`、`trad.f` (現役は photons.f の tradk)、`tserad.f`, `tstrad.f`, `tqrad.f`, `trade.f` (現役は tradke)。`tmuld.f` は一部のみ dead (tmuld6 は photons.f:1494、temap.f:432 から現役)。`tbfrie.f` は現役 (tbedge が tbende/tsteee から)。OBJOBS 66 本ほか非ビルドの .f は個別判定未実施 (未確認)。pysad は dead を一切移植していない。

## 7. 関連文書

- 未対応・未確認の一覧: docs/31。
- ファイル別・ルーチン別・機能別の翻訳状況 (ファイル分類表、tturn / tturne のディスパッチ表を含む): docs/33。
- モジュール地図: docs/01。ベンチマークと回帰: docs/27。K64 のバグ目録: docs/30。
