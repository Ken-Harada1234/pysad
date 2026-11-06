# 01. pysad 全体像 — SAD (K64 版) の Python 直訳

最終更新: 2026-11-06 (現状のスナップショット)

> 場所: コード・文書は `<repo>\` (pysad 本体は `pysad\`、文書は `docs\`、旧文書は `<work-olddocs>\`)。
> 対象の SAD は KEK 本流ではなく **K64 版** (https://github.com/KatsOide/SAD)。コマンド・要素キーワードの意味は CERN の SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html、1.4.0.0k64 2025-07-14) か、ローカルの `<SAD>\Documents\SADHelp.html` (1.4.4.2k64 2026-10-19、使っている WSL ビルドと同じ版) を先に見る。**KEK 公式サイト (acc-physics.kek.jp) は別系統の SAD なので参照しない** (ユーザー指示)。数値の最終根拠はローカルのソース `<SAD>\src` と `Packages\`。
> 再開時はまず `docs/history/00_RESUME.md` (状態・残課題) を読む。

## 1. 目的と方針

- 目的: SAD の光学計算・エミッタンス・マッチング・粒子追跡・DA を Python で再現し、ユーザーの実スクリプト (`ユーザーの esrf ディレクトリ`, `ユーザーの pf2 ディレクトリ`) を `run_sad.py` でそのまま動かす。
- 方針: **再導出せず直訳**。SAD の Fortran の式・分割数・反復順序・乱数の引き方を写し、SAD の数値 (印字精度 1e-10〜1e-11) を再現することを一致の基準にする。SAD の癖 (バグに見えるものも) は再現する (→ `11_sad_quirks.md`)。
- 進め方: SAD スクリプトを下書き → ユーザーが確認・改修 → WSL で実行して参照値 (`bench_level*.result`) → pysad で直訳して一致させる。長い計算は進捗表示 + 発散で停止。返答・コメントは日本語。
- 対象外: Tcl/Tk・X11 描画、EPICS、Python ブリッジ (SAD 側の PySad.c は非コンパイル)。

## 2. SAD の構造

```
MAIN 層 (toplvl.f)   要素定義 / LINE / 数値代入 / FFS, TRACK, EMIT の起動。SADScript は使えない
   └ FFS USE=…;  → track.f (tfinitn で SADScript 初期化、RFSW=ON, CALC6D=OFF 強制) → tffs → tffsa (コマンドループ)
        FFS 層 (tffsa.f)  1 語の解釈順: tfprint (式に見える語は SADScript として評価) → tflag (フラグ) → tgetfv (FIT 条件) → terror (誤差) → tffile → 要素名 (tfgetlineps→tfgetv) → コマンド群 → tfgetv (要素値、失敗なら ?Undefined command) (詳細は 03_ffs_commands §1)
        SADScript 層      tfetok (字句) / tfeval (構文、優先順位は tfstk.f iprior) / tfefun* (組み込み関数)、Packages/*.n
物理計算   CAL/GO: tffsmatch → tffscalc → qcell1 (qcod → qmdiag → qtwiss)
           Emittance[]: emit.n → temit → tturne (tquade/tbende/tthine/tcave + tradke) → tecalc
           TrackParticles[]: tftrack → tturn → tturn1 (tdrift/tquad/tthin/tbend/tcav/tmulti/tsol + tradk)
           DynamicApertureSurvey: dapert.n → tfdapert → trackd → tturn
```

### 主要ソースと役割

| ファイル | 役割 |
|---|---|
| tffsa.f | FFS コマンドループ (USE/CELL/INS/CAL/GO/FIT/FREE/SAVE…)、`FFS[]` (tfffs) |
| tffscalc.f | 1 回分の光学計算 + 残差 (運動量点ループ、重み twfit、nstab) |
| tffsmatch.f | GO の反復本体 (Newton/DESCEND 切替、fuzz、tfsolv/tgrad、fmincube、tweigh) |
| qdcell.f | tdfun (残差 df)、tffsqu/qdcell (解析微分) |
| qcell.f | qcell1: CELL の周期解 (qcod → qmdiag、結合込み分散) と INS、qtwiss 呼び出し |
| qtwiss.f | Twiss 伝搬 (qtwiss1/qmat2twiss)、**qcod** (4D 閉軌道、Newton + SVD)、qquad/qthin (ファイル qquad.f は非コンパイル) |
| temit.f | temit/tecalc: 6D 平衡ビーム行列 (固有値分解 + SVD)、synchm (RF あり) |
| tturne.f | emit 用の放射込み 6D 一周 (tquade/tbende/tthine/tcave/tmulte、trf0 更新) |
| tturn.f | 追跡の一周 (tturn1 の要素分岐、tapert 損失判定) |
| tbend.f / tbendi.f / tbende.f | 追跡 BEND (tbendcore/tbendr、K1≠0 は tbendi) / emit 用 BEND 行列 (qbend も内部) |
| tquad.f / tquade.f | 追跡の四極 tquad・薄い多極 tthin・ttfrin / emit 用四極行列 |
| tsolqu.f / tsolque.f | 四極本体の分割写像 (追跡 tsolqu/tsolqur/tsolqum、emit tsolque)、ソレノイド bz 分岐 |
| tthine.f | emit 用の薄い多極 (L/6-2L/3-L/6) |
| tcav.f / tcave.f | 空洞 (追跡 / emit)、tcavfrin フリンジ |
| tmulti.f / tmulte.f | MULT (追跡 / emit、ANGLE=0) |
| photons.f | 放射: tradk→tradkn (平均) / tradkfn (FLUC、乱数)、tradke (emit 用、radcod) |
| trackd.f / tfdapert.f | DynamicAperture の粒子配置・境界探索・Score |
| tftwiss.f | Twiss[] / LINE[] / Element[] 関数 (値の取り出し。Twiss の計算本体ではない) |
| tgetfv.f / terror.f / tfgetv.f / tflag.f | FIT 条件 / 誤差コマンド DELX… / 要素値 "QF K1 v" / フラグ語 (tffs.f fname/sino) |
| autofg.f / tftype.f | 数値書式 ($FORM) / TYPE の書式 |
| tfstk.f | 中核モジュール (型コード icBEND=2 等、物理定数 macphys、演算子 opcode/iprior) |
| random_plugin_sad.c | 既定乱数 "SAD" (合同法、Box–Muller) |
| Packages/init.n, emit.n, Fit.n, dapert.n, Optics.n | 既定値 (NPARA=4 等)、Emittance[] の規則名、PolynomialFit、DynamicApertureSurvey、NormalCoordinates |

フラグ (RAD, RFSW, FLUC/DAMPONLY, CELL/INS, RADCOD, INTRES…) は tffs.f の `type flagset` と tflag.f。初期値は initb1.f の defflg と track.f のコピー。

## 3. pysad のモジュール地図

| モジュール | 内容 | 直訳元 |
|---|---|---|
| parser.py / sadkw.py / sadelem.py | MAIN 部 (要素定義、LINE、`-LINE`, `N*`, 再定義、LINE 名は書かれたとおり) の読み込み → `Lattice`。K1/K2/K3 は per-length に変換 (元の値は `_raw`)。Element/SetElement/BeamLine/WriteBeamLine | toplvl/doelem/doline/expln、tfelement |
| elements.py | 要素クラス (Drift, Bend, Quadrupole, Sextupole, Octupole, Multipole, Cavity, Solenoid, Marker…) と線形行列 | — |
| lattice.py / elkern.py / perfkern*.py / turnkern.py | CALC: 4D 閉軌道 `_qcod_periodic` (1 周を numba 1 回に)、要素ごとの COD 写像 (tquade/tbende/tthine/tqlfre/tblfre/tchge/tbrote/MULT/tsteee)、`twiss()`/`twiss_table()` (qcell1 結合版、不安定時の 1/cosμ 続行)、`ins_twiss()` (INS の初期軌道オフセット、TRPT の加速)、MARK OFFSET、`set_element`/`set_error`、要素カーネルの numba 化 | qtwiss.f, qcell.f, tquade.f, tbende.f, tthine.f, tchge.f |
| calc6d.py | CALC6D (qcell6d、CELL と INS)、NormalCoordinate6 (etwiss2ri) | qcell6d, temit |
| geo.py | 幾何 (GX/GY/GZ/CHI、GEO、ORG、DISP G/OG、fit の GX..) | tfgeo.f, geolib |
| radint.py | RADINT (放射積分 I1〜I5、結合・色収差) | INTGRL.f |
| trptacc.py | TRPT の加速 (CAVI と VOLT 付き MULT の gammab 分布、tcave/tmulte の acc 部分) | tfgeo1 tfacc, tcave.f, tmulte.f |
| emit6d.py / teigen.py / emit_extra.py | Emittance[] / EMIT: 6D 一周 (T + D、量子励起 B)、tradke、tcave/trf0、tecalc (synchm、固有値 teigen/tnorm の直訳)、close_z、RADCOD、報告の全行、EMITXE.. の下限、Emittance[] のオプション (OneTurnInformation/Matrix) | temit.f, tturne.f, photons.f tradke, tcave.f, tmulte.f, emit.n |
| touschek.py | INTRA (IBS 付き emit)、TouschekTable / TouschekLifetime | temit.f tintrb, emit.n |
| match6.py / match6d.py / tsvd.py | FIT/FREE/FIX/GO: tffsmatch (Newton/DESCEND、fuzz、fmincube、tweigh、tfsolv SVD)、tdfun、twfit、FITP>1 と fam 点 (MatchingAmplitude/NFAMP)、FitValue/FitWeight/FitFunction/VariableRange/VariableWeight、OpticsProlog/Epilog。既定は解析微分 (match6d、qdcell/tffsqu の直訳)、`FFS$NumericalDerivative` で数値微分 | tffsmatch.f, tffscalc.f, qdcell.f, tgetfv.f, tsolvg.f |
| track6.py / track6n.py / track6n_ext.py | TrackParticles: 純 Python 版 (要素外側・粒子内側、tdrift/tquad/tthin/tbend/tbendi/tsteer/tcav/tmulti/tsol、tradk/tradkfn) と numba 版 (粒子配列 S(n,10)、同じ演算順序、乱数も int64)。粒子損失 (tapert、alive フラグ)、MULT/SOL 領域のカーネル | tftrack.f, tturn.f, tquad.f, tsolqu.f, tbend.f, tcav.f, tmulti.f, photons.f |
| dapert.py | DynamicApertureSurvey (スロット順・tapert・Score・62 進地図) | dapert.n, tfdapert.f, trackd.f |
| maintrack.py | MAIN 層の TRACK コマンド (`TRACK USE= NPART= TURNS= …`) | trackm (K64 では動かない機能) |
| optics_fn.py | TransferMatrix / BeamMatrix / NormalCoordinates / OrbitGeo / GeoBase / SymplecticJ / SurvivedParticles / AccelerateParticles | Optics.n, emit.n, dapert.n (Python で書き直し) |
| sadrandom.py | SAD の乱数 (LCG、Box–Muller、drndsr/tdusrn、ParabolaRandom) | random_plugin_sad.c, tfRandom_.c, drndsr.f |
| sadscript.py / sadlist.py | SADScript インタプリタ (字句・Pratt 構文・評価、Module/Do/Table/純関数/パターン/REPEAT/UNTIL/IF、単位付き数値、`?FLAG`、fmt_real=autofg) とリスト・文字列関数 | tfetok.f, tfeval.f, tfstk.f iprior, autofg.f, functions.n, Misc.n |
| ffscmd.py | FFS コマンド層 (tffsa の解釈順)、USE/CELL/INS/CALC/CALC6D/GO/EMIT/TYPE/DISP/SAVE/RESET/RADINT/CHROMATICITY/誤差コマンド (DELX…)、Twiss[]/LINE[]/FFS[]/Element[]/Emittance[]/TrackParticles[]/CalculateOptics/SeedRandom[]、EMIT の結果変数、FFS$SHOW | tffsa.f, tftwiss.f, tflag.f, terror.f, tftype.f, tfdisp.f |
| diag_fit.py / diag_trace.py | 照合用の診断 (fit の反復表、追跡の途中値) | — |
| run_sad.py (<repo> 直下) | `py -3.12 run_sad.py lattice script …` = `gs` に複数ファイルを連結して渡すのと同じ | — |
| **非推奨** ffs.py / matching.py / tracking.py / optics.py | 独自実装 (SAD 未照合。tracking.py は放射係数に疑い)。新規コードは使わない。`__init__.py` の export は旧 API のまま残っている | — |

テスト: `pysad/tests/` (pytest 277 件。`--ignore` する 4 本は sys.exit で終わる旧式のスクリプト)。
ベンチ: `bench_level1〜31*.sad/.result` と `*_pysad.py/.result`、`run_regression.py` (50 ベンチ、`QUICK=1` で短縮)。
ユーザースクリプトの照合: `<work-esrf>/`、`<work-pf2>/` (SAD 出力 `*.sad.out`、pysad 出力 `*.pysad.out`)。

## 4. docs/ の案内

| 文書 | 内容 |
|---|---|
| history/00_RESUME.md | 作業メモ (時系列、残課題)。再開時に最初に読む。履歴はここと work_*.md にだけ書く |
| 01_overview.md | 本書 |
| 02〜09 | SAD の構造と各計算の直訳の解説: MAIN/LINE/要素 (02)、FFS コマンド (03)、SADScript (04)、CALC/COD/Twiss (05)、マッチング (06)、Emittance (07)、追跡 (08)、DynamicApertureSurvey (09) |
| 10_benchmarks.md | ベンチ (bench_level*) の一覧と現在の一致度 |
| 11_sad_quirks.md | SAD の癖・仕様・ビルド依存の一覧 (現象 / 根拠 / pysad での扱い) |
| 12_user_scripts_status.md | ユーザースクリプト 86 本の現在の状態と未対応項目 |
| 13_pysad_api.md | pysad の Python API |
| 14_fortran_status.md | SAD の Fortran ファイル別の翻訳状況 |
| 15_implementation_summary.md | 実装まとめ (モジュール行数、SAD 準拠状況、課題) |
| 16_physics_guide.md / 17_source_analysis.md / 18_translation_candidates.md | 物理の式 / SAD ソースの解析 / ユーザースクリプト × SAD 機能の翻訳候補 |
| 19_sad_bug_report.md | SAD K64 のバグ報告 (B01〜) |
| 13_sad_help_index.md / 14_pysad_vs_sadhelp.md / 15_sadhelp_purpose_policy.md / sadhelp_summary.html | SADHelp の一覧、pysad との対応、目的別分類と翻訳方針 (<work-perf>/sadhelp_*.py で生成) |
| work_*.md | 項目ごとの作業メモ (経緯・照合の記録) |
| html/ | `make_html.py` で docs/*.md から生成した HTML |

旧文書 (`<work-olddocs>/`) はソース照合済みの初期の解説で、pysad の状況は古い。現状は本書と上の各文書を見る。

## 5. 使い方

```
# SAD スクリプトをそのまま実行 (MAIN 部 + FFS 部)。numba 入りの Python 3.12 を使う
cd "<repo>"
py -3.12 run_sad.py <work-esrf>\pysad/tests/kek-ls_v03_12.lat <work-esrf>\matching_full.sad

# Python API (詳細は 13_pysad_api.md)
from pysad.parser import load_sad_file      # → SADScript (env に Lattice と LINE)
lat.twiss_table(delta=0.0)                   # CELL CALC (qcell1)。lat は Lattice
from pysad.emit6d import emittance6d; emittance6d(lat, 3e9, rfsw=True, radcod=False)
from pysad.match6 import Matcher; m = Matcher(lat); m.fit('NX', 47.6); m.free('Q*'); m.go()
from pysad.track6n import track_particles_fast   # TrackParticles の本体
from pysad.dapert import dynamic_aperture        # DynamicApertureSurvey の本体

# 回帰 (SAD 参照値と比較)
QUICK=1 py -3.12 run_regression.py
py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py --ignore=pysad/tests/test_track6n_ext.py

# SAD 本体 (参照値の取得): WSL Ubuntu
MSYS_NO_PATHCONV=1 wsl -d Ubuntu -- <SAD>/bin/gs file.sad > file.result   # NPARA=1; を入れる
```

## 6. 到達点と未対応

SAD と一致しているもの (詳細は 10_benchmarks.md、12_user_scripts_status.md):
- 線形光学 (CELL/INS、CALC6D、ROTATE/ずれ付きの全要素、ANGLE=0 BEND の K0/K1、TRPT の加速と MARK OFFSET、幾何 GX..、RADINT) ≤1e-10 (誤差・オフモメンタム付きは ≤1e-7)。
- Emittance[]/EMIT (RFSW × RADCOD の 4 通り、INTRA、不安定な光学、報告の全行) ≤1e-11〜1e-9、Touschek ≤7e-12。
- マッチング (解析微分、DESCEND、fuzz、tsvd、2 段目フィット、INS、フック) は反復列・残差・K1 を SAD と同じ桁で再現。
- 追跡 (NORAD/RAD/FLUC、粒子損失) ≤5e-11、DA の地図と Score、誤差の実現 (DELX.. の乱数列)。
- SADHelp の関数 94 中 46、特別変数 48 中 42、FFS コマンド 56 中 32 (無視 12) に対応 (14_pysad_vs_sadhelp.md)。

| 未対応・制限 | 備考 |
|---|---|
| 性能 | 誤差の実現 ×多数・DAsv 5000 周・チューンサーベイ 54×54 の本番サイズは長い (動作テストは <work-userscripts>/ の間引き版)。1 周の追跡は 176 ns/粒子・要素 |
| FITP>1 の最終解 | SAD (WSL) が GO で再起動するため参照値がない。1 段目・初回の解は一致、2 段目以降の多解は別解に着くことがある。fam 点は自己検証のみ |
| 不安定ラティス | SAD の 1e199 級の値や NaN は再現しない (安定域は一致) |
| 損失周の乱数 | 論理は SAD と同じだが、実寸 DAsv (5000 周) の最終 seed はカオスで一致しない (Score の差 1〜3 %) |
| 要素 | SOL の emit/幾何 (追跡のみ)、MULT の ANGLE≠0 (tmulta)、TCAVI、COORD、MAP、UND、BEAMBEAM、POL (偏極) |
| コマンド・関数 | MEASURE/TRCOD/DRAW/GEO 描画、SynchroBetaEmittance (SAD K64 に実体なし)、InitialOrbits/ElementValues (使用なし)、GaussianCoulomb/WakeFunction/RadiationField (優先度 低) |
| 対象外 | Tk/plot/EPICS/System、ファイル入出力の失敗の扱い (Python ネイティブで行う)。SAD 自身も WSL では draw/Tk で core dump |
| 意図的な差 | SADScript の `1GEV` は pysad では単位を掛ける (SAD は `1*GEV` = 1、11_sad_quirks C24) |

## 7. ユーザーの実スクリプトの対応状況

| 対象 | 本数 | 数値一致 | 一部一致 | 動くが不一致 | 例外 | 打ち切り | SAD 側無効 | 対象外 |
|---|---|---|---|---|---|---|---|---|
| esrf (`ユーザーの SAD スクリプトのディレクトリGSR\esrf`、<work-userscripts>/esrf) | 52 | 34 | 12 | 0 | 1 | 0 | 4 | 1 |
| pf2 (`ユーザーの pf2 ディレクトリ`、<work-userscripts>/pf2) | 34 | 28 | 6 | 0 | 0 | 0 | 0 | 0 |

例外 1 本は入力ファイルが無いもの、SAD 側無効 4 本は USE 先がラティスに無いもの、対象外 1 本は System[]。打ち切りなし (誤差の実現・DA・チューンサーベイは間引き)。
一部一致 18 本の内訳 (fit の多解・FITP>1、不安定域の縁、DA の乱数のカオス、誤差系の要調査 2 件、小さな数値差) は 12_user_scripts_status.md §4。
