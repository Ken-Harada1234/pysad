# 01. pysad の全体像と文書の案内

pysad は、加速器コード SAD (K64 版) の計算 (光学・エミッタンス・マッチング・粒子追跡・ダイナミックアパーチャ) を Python で再現し、SAD の台本をそのまま動かせるようにしたものである。

- **この文書で分かること**: pysad が何を計算できるか、SAD のどの部分がどのモジュールに対応するか、動かし方、どの目的のときにどの文書を読めばよいか (全 36 本の目次)。
- **前提**: 加速器物理 (光学・エミッタンス・追跡) の基本を知っていること。SAD と pysad は知らなくてよい。
- **読む順**: まず §1・§2 (目的と文書の選び方)、§3 (動かしてみる)。SAD のしくみは §4、pysad のモジュールの対応は §5。全体の目次は §8。
- **用語**: **SAD** = KEK の K. Oide らの加速器設計コード。**K64** = その GitHub 版 (https://github.com/KatsOide/SAD)。**gs** = K64 の実行ファイル。**MAIN 部** = 要素と LINE を定義する層。**FFS** = SAD の対話・台本の層 (`FFS;` 以降)。**SADScript** = SAD の台本言語。**v2 モード / k64 モード** = pysad の 2 つの動作 (§1)。

> 場所: コードは `<repo>\` (pysad 本体は `pysad\`、文書は `docs\`、図つきの HTML は `docs\html\`)。
> 対象の SAD は KEK 本流ではなく **K64 版**。コマンド・要素キーワードの意味は CERN の SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html、1.4.0.0k64 2025-07-14) か、ローカルの `<SAD>\Documents\SADHelp.html` (1.4.4.2k64 2026-10-19、使っている WSL ビルドと同じ版) を先に見る。**KEK 公式サイト (acc-physics.kek.jp) は別系統の SAD なので参照しない**。数値の最終根拠はローカルのソース `<SAD>\src` と `Packages\`。

## 1. 目的と方針

- 目的: SAD の光学計算・エミッタンス・マッチング・粒子追跡・DA を Python で計算し、SAD スクリプトを `run_sad.py` でそのまま動かす。
- 2 つの動作モード (互換モード [`pysad/compat.py`](../pysad/compat.py)):
  - **v2 モード (既定)**: K64 のバグ・物理の誤りを直し、機能を足した動作。利用者はふつうこちらを使う。
  - **k64 モード**: K64 に忠実な翻訳 (SAD の Fortran の式・分割数・反復順序・乱数の引き方を写し、バグに見える癖も再現する)。`PYSAD_COMPAT=k64` か `pysad.compat.set_mode('k64')` で切り替える。gs と数値を比べるときはこちら。既存の試験と回帰は k64 で走る。
  - 2 つの違いは [28_v2_vs_k64.md](28_v2_vs_k64.md)、K64 自身のバグ目録は [30_sad_bug_report.md](30_sad_bug_report.md)、SAD の癖は [29_sad_quirks.md](29_sad_quirks.md)。
- k64 モードの基準: 1e-10 以下の差、表示の桁、丸めカオス後の乱数 seed は追わない。「計算が完走し、物理量が SAD と有意に一致する」ことを基準にする。
- 翻訳の範囲: Fortran/C の物理計算は直接翻訳する。`Packages/*.n` は SADScript のまま pysad のインタプリタで実行するか、Python で書き直す。ファイル・リスト・GUI・描画は Python のライブラリで代替し翻訳しない ([36_sadhelp_purpose_policy.md](36_sadhelp_purpose_policy.md))。
- 実行環境: Python 3.12 (`py -3.12`) に numba / scipy / pytest。`PYTHONIOENCODING=utf-8`。参照値 (gs) は WSL Ubuntu の `bin/gs` で取る (`NPARA=1;` を入れる)。長い計算は優先度を下げて 1 本ずつ流す。
- 対象外: EPICS、Python ブリッジ (SAD 側の PySad.c は非コンパイル)。`System[]` は既定で実行しない (環境変数 `PYSAD_ALLOW_SYSTEM=1` で許可。[05_sadscript.md](05_sadscript.md))。Tcl/Tk の GUI と描画 (`ListPlot`・`Plot`・`OpticsPlot` など) は、SAD 自身の `Packages/*.n` を pysad のインタプリタで解釈して tkinter に表示する方式で対応している ([20_tcltk_gui.md](20_tcltk_gui.md)、[21_graphics.md](21_graphics.md))。見た目 (色・大きさ・フォント)・操作の動作は SAD と比べておらず、`Canvas3D*` は未対応。既知の差: 最初の `Window[]` の親が、gs では Toplevel、pysad ではルート窓。`Options[g]` は gs が全オプションを返し、pysad は `{}`。`OpticsPlot` の text 項目数は Tk のフォントによる `RemoveOverlap` の差で違う。`w=Window[]`・`e=Entry[w]` は `Widget[…]` 値になり、`e[Insert]`・`e[Get]`・`e[Delete]` は gs と同じく働く。
- 扱う要素: DRIFT / BEND / QUAD / SEXT / OCT / DECA / DODECA / MULT (ANGLE≠0 を含む) / SOL / CAVI / TCAVI / UND / COORD / MAP / BEAMBEAM / MARK / MONI / APERT / INS。PHSROT は K64 で動作しない (行列が 0、B38) ため、k64 モードでは追跡が NotImplementedError、v2 モードは意図された写像を実装した ([10_tracking.md](10_tracking.md) §9.7)。SPCH は K64 の `tturn.f:521` の `case (icSPCH)` が PSPAC フラグのとき `tpspac` (粒子のメッシュによる空間電荷) を呼ぶ要素で、**v2 モードでは SPAC・PSPAC・WSPAC が `TrackParticles` に入る** (k64 モードは入らない)。PROF は K64 のソースに現れず、parser が読むだけ。要素別の状況は [19_other_elements.md](19_other_elements.md) ほか各文書。
- 未対応・未確認の唯一の一覧は [31_unimplemented_and_unverified.md](31_unimplemented_and_unverified.md)。

## 2. 目的別: どの文書を読むか

初めての人は、01 (本書) → 02 (ラティスを書く) → 03 (FFS で計算を指示する) → 07 (光学) → 09 (エミッタンス) → 10 (追跡) の順に読むと、SAD の標準的な作業の流れ (ラティス定義 → 光学 → エミッタンス → 追跡) に沿う。次の表は、やりたいことから読む文書を引く索引である。

| やりたいこと | 読む文書 |
|---|---|
| ラティス (要素と LINE) を書く・読み込む | [02_lattice_main.md](02_lattice_main.md) |
| FFS のコマンド (USE・CELL・CALC・GO・TYPE・DISP など) を使う | [03_ffs_commands.md](03_ffs_commands.md)、小さなコマンドは [04_ffs_small_commands.md](04_ffs_small_commands.md) |
| SADScript (台本言語) を書く・特殊関数・Fit・SAD の Packages を使う | [05_sadscript.md](05_sadscript.md)、[06_sadscript_library.md](06_sadscript_library.md) |
| チューン・β 関数・分散・閉軌道・クロマティシティを求める | [07_calc_optics.md](07_calc_optics.md) |
| 光学を目標値に合わせる (マッチング) | [08_matching.md](08_matching.md) |
| エミッタンス・エネルギー広がり・放射減衰を求める | [09_emittance.md](09_emittance.md) |
| 粒子を追跡する・ダイナミックアパーチャを求める | [10_tracking.md](10_tracking.md)、[11_dynamic_aperture.md](11_dynamic_aperture.md) |
| Python から pysad を呼ぶ | [12_pysad_api.md](12_pysad_api.md) |
| ウェイク・ビーム・ビーム・CSR を入れる | [13_wake.md](13_wake.md)、[14_beambeam.md](14_beambeam.md)、[15_csr.md](15_csr.md) |
| スピン偏極・ソレノイド・放射光 (RADLIGHT) を扱う | [16_spin_pol.md](16_spin_pol.md)、[17_sol.md](17_sol.md)、[18_radlight.md](18_radlight.md) |
| COORD・TCAVI・UND・曲げつき MULT・MAP を使う | [19_other_elements.md](19_other_elements.md) |
| GUI・作図・CAD 図を使う | [20_tcltk_gui.md](20_tcltk_gui.md)、[21_graphics.md](21_graphics.md)、[22_lattice_dxf.md](22_lattice_dxf.md) |
| 計算の式 (物理・理論) を知りたい | [23_physics_guide.md](23_physics_guide.md)。ADTS・RDT・周波数マップは [24_adts_rdt_primer.md](24_adts_rdt_primer.md) → [25_resonance_terms_adts.md](25_resonance_terms_adts.md) → [26_frequency_map_lma.md](26_frequency_map_lma.md) |
| どの程度信頼できるか (検証・SAD との違い) を知りたい | [27_benchmarks.md](27_benchmarks.md)、[28_v2_vs_k64.md](28_v2_vs_k64.md)、[29_sad_quirks.md](29_sad_quirks.md)、[30_sad_bug_report.md](30_sad_bug_report.md) |
| 何ができないか・未確認かを知りたい | [31_unimplemented_and_unverified.md](31_unimplemented_and_unverified.md) |
| SAD のソースや翻訳の状況 (開発者向け) | [32_source_analysis.md](32_source_analysis.md)、[33_fortran_status.md](33_fortran_status.md)、[34_sad_help_index.md](34_sad_help_index.md)、[35_pysad_vs_sadhelp.md](35_pysad_vs_sadhelp.md)、[36_sadhelp_purpose_policy.md](36_sadhelp_purpose_policy.md) |

## 3. 動かしてみる

手順 (`<repo>` の直下で実行):

1. SAD の台本を実行する。ラティスの台本と計算の台本を続けて渡す (SAD の `gs` に複数ファイルを連結して渡すのと同じ)。

   ```
   cd "<repo>"
   set PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>
   py -3.12 run_sad.py example/00_lattice/FODO_ring.sad example/02_optics_cod/FFS_CALCULATE.sad
   ```

   24 セルの FODO リング (3 GeV、周長 105.6 m) の最後に `NX=  6.696329  NY=  4.906056` と出る。これは水平・垂直のチューン (νx, νy) = (6.696, 4.906) である。`PYSAD_SAD_PACKAGES` は SAD 自身の `Packages/*.n` (CorrectOrbit・GUI など) を使うときに要る。
2. Python から使う (詳細は [12_pysad_api.md](12_pysad_api.md))。

   ```python
   from pysad.parser import load_sad_file          # → SADScript (env に Lattice と LINE)
   lat.twiss_table(delta=0.0)                       # CELL CALC (qcell1)。lat は Lattice
   from pysad.emit6d import emittance6d; emittance6d(lat, 3e9, rfsw=True, radcod=False)
   from pysad.match6 import Matcher; m = Matcher(lat); m.fit('NX', 47.6); m.free('Q*'); m.go()
   from pysad.track6n import track_particles_fast   # TrackParticles の本体
   from pysad.dapert import dynamic_aperture        # DynamicApertureSurvey の本体
   ```
3. 試験と回帰を走らせる (§7)。
4. gs で参照値を取る (WSL Ubuntu): `MSYS_NO_PATHCONV=1 wsl -d Ubuntu -- <SAD>/bin/gs file.sad > file.result` (台本に `NPARA=1;` を入れる)。

## 4. SAD の構造

SAD は MAIN 層 → FFS 層 → SADScript 層の 3 層で、計算本体は FFS から呼ばれる。pysad も同じ層に分けて対応させている。

```
MAIN 層 (toplvl.f)   要素定義 / LINE / 数値代入 / FFS, TRACK, EMIT の起動。SADScript は使えない
   └ FFS USE=…;  → track.f (tfinitn で SADScript 初期化、RFSW=ON, CALC6D=OFF 強制) → tffs → tffsa (コマンドループ)
        FFS 層 (tffsa.f)  1 語の解釈順: tfprint (式に見える語は SADScript として評価) → tflag (フラグ) → tgetfv (FIT 条件) → terror (誤差) → tffile → 要素名 (tfgetlineps→tfgetv) → コマンド群 → tfgetv (要素値、失敗なら ?Undefined command) (詳細は 03_ffs_commands.md §1)
        SADScript 層      tfetok (字句) / tfeval (構文、優先順位は tfstk.f iprior) / tfefun* (組み込み関数)、Packages/*.n
物理計算   CAL/GO: tffsmatch → tffscalc → qcell1 (qcod → qmdiag → qtwiss)
           Emittance[]: emit.n → temit → tturne (tquade/tbende/tthine/tcave + tradke) → tecalc
           TrackParticles[]: tftrack → tturn → tturn1 (tdrift/tquad/tthin/tbend/tcav/tmulti/tsol + tradk)
           DynamicApertureSurvey: dapert.n → tfdapert → trackd → tturn
```

次の表は、SAD の計算で中心になるソースファイルと役割である (ソースを読むときの索引。構造の詳細は [32_source_analysis.md](32_source_analysis.md)、翻訳状況は [33_fortran_status.md](33_fortran_status.md))。

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
| tmulti.f / tmulte.f、tmulta.f | MULT (追跡 / emit、ANGLE=0 と ANGLE≠0) |
| photons.f | 放射: tradk→tradkn (平均) / tradkfn (FLUC、乱数)、tradke (emit 用、radcod) |
| trackd.f / tfdapert.f | DynamicAperture の粒子配置・境界探索・Score |
| tftwiss.f | Twiss[] / LINE[] / Element[] 関数 (値の取り出し。Twiss の計算本体ではない) |
| tgetfv.f / terror.f / tfgetv.f / tflag.f | FIT 条件 / 誤差コマンド DELX… / 要素値 "QF K1 v" / フラグ語 (tffs.f fname/sino) |
| autofg.f / tftype.f | 数値書式 (`$FORM`) / TYPE の書式 |
| tfstk.f | 中核モジュール (型コード icBEND=2 等、物理定数 macphys、演算子 opcode/iprior) |
| random_plugin_sad.c | 既定乱数 "SAD" (合同法、Box–Muller) |
| Packages/init.n, emit.n, Fit.n, dapert.n, Optics.n | 既定値 (NPARA=4 等)、Emittance[] の規則名、PolynomialFit、DynamicApertureSurvey、NormalCoordinates |


フラグ (RAD, RFSW, FLUC/DAMPONLY, CELL/INS, RADCOD, INTRES…) は tffs.f の `type flagset` と tflag.f。初期値は initb1.f の defflg と track.f のコピー。
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
| tmulti.f / tmulte.f、tmulta.f | MULT (追跡 / emit、ANGLE=0 と ANGLE≠0) |
| photons.f | 放射: tradk→tradkn (平均) / tradkfn (FLUC、乱数)、tradke (emit 用、radcod) |
| trackd.f / tfdapert.f | DynamicAperture の粒子配置・境界探索・Score |
| tftwiss.f | Twiss[] / LINE[] / Element[] 関数 (値の取り出し。Twiss の計算本体ではない) |
| tgetfv.f / terror.f / tfgetv.f / tflag.f | FIT 条件 / 誤差コマンド DELX… / 要素値 "QF K1 v" / フラグ語 (tffs.f fname/sino) |
| autofg.f / tftype.f | 数値書式 (`$FORM`) / TYPE の書式 |
| tfstk.f | 中核モジュール (型コード icBEND=2 等、物理定数 macphys、演算子 opcode/iprior) |
| random_plugin_sad.c | 既定乱数 "SAD" (合同法、Box–Muller) |
| Packages/init.n, emit.n, Fit.n, dapert.n, Optics.n | 既定値 (NPARA=4 等)、Emittance[] の規則名、PolynomialFit、DynamicApertureSurvey、NormalCoordinates |


フラグ (RAD, RFSW, FLUC/DAMPONLY, CELL/INS, RADCOD, INTRES…) は tffs.f の `type flagset` と tflag.f。初期値は initb1.f の defflg と track.f のコピー。

## 5. pysad のモジュール地図

pysad は `pysad/` の 102 本の Python ファイル (計 78235 行。2026-11-23 に実測) でできている。次の表は、SAD の機能ごとにどのモジュールが受け持つか (内容)、行数、翻訳元の SAD のファイル (直訳元) を並べたものである。表は上から、MAIN 部・光学・エミッタンス・マッチング・追跡・機能別・台本・FFS の順に近い。

| モジュール | 行数 | 内容 | 直訳元 |
|---|---:|---|---|
| parser.py / sadkw.py / sadelem.py | 2282 / 34 / 440 | MAIN 部 (要素定義、LINE、`-LINE`, `N*`, 再定義、LINE 名は書かれたとおり) の読み込み → `Lattice`。K1/K2/K3 は per-length に変換 (元の値は `_raw`)。Element/SetElement/BeamLine/WriteBeamLine | toplvl/doelem/doline/expln、tfelement |
| elements.py | 896 | 要素クラス (Drift, Bend, Quadrupole, Sextupole, Octupole, Multipole, Cavity, Solenoid, Marker…) と線形行列 | — |
| lattice.py / elkern.py / perfkern*.py (perfkern_bend.py を含む) / turnkern.py | 4644 / 1873 / 295 / 137 / 480 | CALC: 4D 閉軌道 `_qcod_periodic` (1 周を numba 1 回に)、要素ごとの COD 写像 (tquade/tbende/tthine/tqlfre/tblfre/tchge/tbrote/MULT/tsteee)、`twiss()`/`twiss_table()` (qcell1 結合版、不安定時の 1/cosμ 続行)、`ins_twiss()` (INS の初期軌道オフセット、TRPT の加速)、MARK OFFSET、`set_element`/`set_error`、要素カーネルの numba 化 | qtwiss.f, qcell.f, tquade.f, tbende.f, tthine.f, tchge.f |
| calc6d.py | 690 | CALC6D (qcell6d、CELL と INS)、NormalCoordinate6 (etwiss2ri) | qcell6d, temit |
| geo.py | 580 | 幾何 (GX/GY/GZ/CHI、GEO、ORG、DISP G/OG、fit の GX..) | tfgeo.f, geolib |
| radint.py | 831 | RADINT (放射積分 I1〜I5、結合・色収差) | INTGRL.f |
| wake.py / wakeopt.py | 308 / 200 | ウェイク場 (txwake、spsort。RF のある MULT は track6n_ext.mult_full)、光学の wake (WAKEOPT) と InitialOrbits | txwake.f, spkick.f, tffswake.f, tmulti.f tmultiacc |
| trptacc.py | 348 | TRPT の加速 (CAVI と VOLT 付き MULT の gammab 分布、tcave/tmulte の acc 部分) | tfgeo1 tfacc, tcave.f, tmulte.f |
| emit6d.py / teigen.py / emit_extra.py | 2845 / 490 / 358 | Emittance[] / EMIT: 6D 一周 (T + D、量子励起 B)、tradke、tcave/trf0、tecalc (synchm、固有値 teigen/tnorm の直訳)、close_z、RADCOD、報告の全行、EMITXE.. の下限、Emittance[] のオプション (OneTurnInformation/Matrix) | temit.f, tturne.f, photons.f tradke, tcave.f, tmulte.f, emit.n |
| solmap.py | 118 | SOL 領域の COD 写像と行列 (CALC/INS の Twiss。追跡のカーネルを数値微分。docs/17) | tsole.f tsole1 |
| solemit.py | 453 | SOL 領域の放射込みの 1 周 (Emittance[]・EMIT: bz 付き tradke、tsolque の刻み。docs/17) | tsole.f, tsolque.f, tqente.f, photons.f tradke |
| multa.py / multae.py | 206 / 224 | MULT ANGLE≠0 (追跡 / COD 写像・行列。docs/19) | tmulta.f (tmulta, tmultae) |
| touschek.py | 1032 | INTRA (IBS 付き emit)、TouschekTable / TouschekLifetime。`touschek_piwinski` (独立な Piwinski 型の式。v2 の `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定。`Formula->"K64"`・k64 モードは K64 の表。docs/26 §4.4、B162) | temit.f tintrb, emit.n |
| match6.py / match6d.py / tsvd.py | 1469 / 1665 / 404 | FIT/FREE/FIX/GO: tffsmatch (Newton/DESCEND、fuzz、fmincube、tweigh、tfsolv SVD)、tdfun、twfit、FITP>1 と fam 点 (MatchingAmplitude/NFAMP)、FitValue/FitWeight/FitFunction/VariableRange/VariableWeight、OpticsProlog/Epilog。既定は解析微分 (match6d、qdcell/tffsqu の直訳)、`FFS$NumericalDerivative` で数値微分 | tffsmatch.f, tffscalc.f, qdcell.f, tgetfv.f, tsolvg.f |
| track6.py / track6n.py / track6n_ext.py | 3484 / 2687 / 1385 | TrackParticles: 純 Python 版 (要素外側・粒子内側、tdrift/tquad/tthin/tbend/tbendi/tsteer/tcav/tmulti/tsol、tradk/tradkfn) と numba 版 (粒子配列 S(n,10)、同じ演算順序、乱数も int64)。粒子損失 (tapert、alive フラグ)、MULT/SOL 領域のカーネル。POL かつ RAD ON のときのスピン追跡 (`sprot`・`bsi`。純 Python 版のみで損失の処理なし、docs/16) | tftrack.f, tturn.f, tquad.f, tsolqu.f, tbend.f, tcav.f, tmulti.f, photons.f |
| dapert.py | 255 | DynamicApertureSurvey (スロット順・tapert・Score・62 進地図) | dapert.n, tfdapert.f, trackd.f |
| maintrack.py | 369 | MAIN 層の TRACK コマンド (`TRACK USE= NPART= TURNS= …`) | trackm (K64 では動かない機能) |
| optics_fn.py | 315 | TransferMatrix / BeamMatrix / NormalCoordinates / OrbitGeo / GeoBase / SymplecticJ / SurvivedParticles / AccelerateParticles | Optics.n, emit.n, dapert.n (Python で書き直し) |
| sgeo.py | 315 | SOL 領域の幾何 (tsgeo。k64 モードは領域内の BEND・MULT・CAVI と tdgeo の領域幾何が未対応。v2 モードは対応。docs/17) | tsgeo.f |
| gausscoul.py / bbbrem.py | 472 / 162 | GaussianCoulomb(U・Fitted) (bbkick1・twspu・twspfu。docs/06)、BBBrem1 (docs/14) | tspac.f, tintrb.f, bint.f, tbbbrem.f |
| beambeam.py | 739 | BEAMBEAM 要素 (強弱ビーム、ビームストラールング、ルミノシティ。追跡のみで Twiss には効かない。docs/14) | bbstrhl.f, bbtool.f, matrix.f |
| coord.py / tcav.py / und.py | 147 / 200 / 297 | COORD (追跡・6 次元写像・幾何、docs/19)、TCAVI (追跡 RAD FLUC まで、Emittance[] は L=0 のみ、docs/19)、UND (K64 は磁場を設定できず磁場なしのドリフト。v2 モードは `FBX`・`FBY` の磁場つきで追跡・Twiss・Emittance[] に入る、docs/19) | tcoord.f, tcoorde.f, qcoord.f, ttcav.f, ttcave.f, undulator.f |
| compat.py | 194 | **互換モード**。`mode()`・`set_mode()`・`v2()`・`k64()`。環境変数 `PYSAD_COMPAT=k64\|v2`。docs/12 §16.6、docs/28 | — (pysad 独自) |
| phsrot.py | 170 | PHSROT (位相空間の回転。v2 モードだけ。追跡・CALC・Emittance[]・EMIT。docs/10 §9.7) | phsrot.f, matrix.f |
| radtaper.py | 137 | RADTAPER (磁石の強さを閉軌道の局所の運動量に合わせる。Emittance[]・EMIT・追跡。両モード。docs/10 §9.7) | tcod.f, tturne.f, tturn.f の radtaper |
| spch.py | 1164 | 空間電荷 (SPAC・PSPAC・WSPAC。v2 モードだけ。`TrackParticles` に入る。docs/10 §9.7。Emittance[]・EMIT・CALC6D の WSPAC は `touschek.py`・`calc6d.py` に入っている: docs/09) | spkick.f, tspch_.c・spch.c・psn.c, tintrb.f twspac |
| specfn.py | 1132 | 特殊関数・数学関数 (Bessel・超幾何・ガンマ族・ゼータ・多重対数・D.n の特殊関数など。K64 の InverseErf の誤差も再現。docs/06) | gamma.f, hg.f, Packages/D.n |
| naff.py | 135 | NAFF 風のチューン抽出 (`naff`・`naff_real`・`tune_xp`・`tunes_from_coords`。放射を切った追跡の座標列からチューンを取る。MAIN の TRACK の FOURIER の代用。docs/12 §16.7、docs/28) | K64 の直訳ではない (wtunem・wjfit に代わる独立な実装) |
| specfn_cx.py | 254 | 複素の Zeta・HurwitzZeta・PolyLog・HypergeometricU (両モード。gs が誤る点でも真値。docs/06) | gamma.f の複素の部分に代わる独立な実装 |
| specfn_cx2.py | 314 | 複素の ExpIntegralE・ExpIntegralEi・PolyGamma・WhittakerM/W・ProductLog・CosIntegral・FresnelS/C・HarmonicNumber・GammaRegularized・Beta・Factorial・`Gamma[a,z]`・LegendreP・複素の次数の BesselJ (両モード。複素の引数が 1 つでもあれば `specfn.install` が切り替える。docs/06) | scipy の複素対応と級数・漸化式 |
| spch_k64.py | 287 | PSPAC の K64 の形状関数 (k64 モードだけ。CIC・「ほぼ最近接」・5 点差分 + 双 3 次補間・3 次スプライン) | spch.c・psn.c・tspch_.c |
| trksetup.py | 408 | 追跡の設定の共通部分 (`TrackSetup` = 閉軌道の求め直し・`track`・RF バケット、`find_periodicity`、`select_locations`、優先度を下げた並列。docs/26 §2) | K64 にない追加 (`track6n.py` を使う) |
| freqmap.py | 474 | 周波数マップ (1 つの位置の格子から NAFF でチューンの拡散 D、表・図、共鳴線の表。台本 `FrequencyMap[]`。docs/26 §3) | K64 にない追加 (naff.py・track6n.py を使う) |
| lma.py | 600 | 局所運動量アクセプタンス (位置ごとの δ±、Touschek への受け渡し `touschek_lifetime_local`。台本 `MomentumAcceptance[]`・`TouschekLifetimeLMA[]`。docs/26 §4) | K64 にない追加 (touschek.py の表を使う) |
| lmabucket.py | 224 | RF バケットの限界とアパーチャの限界の区別 (解析式の δ_rf、縦方向の 1 次元の写像、`limit_types`。docs/26 §4.5) | K64 にない追加 |
| fmapcheck.py | 310 | 周波数マップの信頼性の確認 (誤差の統計、警告 12 種、D の周回数の換算、6 極のフィードダウンの見積もり。docs/26 §3.6・§10) | K64 にない追加 |
| q1cmd.py | 243 | 周波数マップ・LMA・共鳴駆動項などの SAD の台本の関数の登録 (`ffscmd.py` の入口から `install`) | K64 にない追加 |
| resdrive.py | 1547 | 6 極の共鳴駆動項・ADTS・クロマティシティ (厚いレンズの解析式と数値積分。台本 `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]`。docs/25) | K64 にない追加 (運動学・端部のフリンジの ADTS。追跡 + NAFF + フィットの `adtsdirect.py`。偏向磁石の 3 次の項を追跡の写像から取り出す `bendnl.py`。共鳴駆動項にも偏向磁石の 3 次の項 (δ の項は分散による軌道のずれ) を入れ、追跡の 1 周の写像から RDT を取り出す `rdttrack.py` で検証) |
| fmapres.py | 348 | 周波数マップの D の大きい領域を、RDT・ADTS から予想した共鳴の領域と比べる (誤差の与え方・統計。docs/26 §3.5) | K64 にない追加 |
| fmappred.py | 354 | 周波数マップの共鳴線の予想の作り直し (3 次・4 次の h と非共鳴の正規形の補正、線形の ADTS、振り子の島の幅、連続のスコア、信頼度。docs/25 §12) | K64 にない追加 |
| rdt2nd.py・nform.py・rdttrack.py・bendnl.py・adtsdirect.py | 422 / 53 / 445 / 642 / 306 | 4 次の共鳴駆動項 (6 極の 2 次)・数値の正規形・追跡の 1 周の写像からの RDT の取り出し・偏向磁石の生成関数 (追跡の写像の 2 次のテイラー係数から)・追跡 + NAFF の ADTS (docs/25) | K64 にない追加 |
| sxtopt.py | 935 | [試作] 6 極ファミリーの K2L の多項式モデルと SLSQP の最適化 (RDT・ADTS の解析式を k の線形・二次形式で組む。検証用。台本の関数なし。[`example/ddba/v07_14/sextupole/`](../example/ddba/v07_14/sextupole/)) | K64 にない追加 |
| constants.py | 10 | 物理定数 (光速、電子の古典半径、電子・陽子の質量、Cγ、Cq など) | pysad の定数 |
| nbcache.py | 146 | numba のキャッシュの無効化 (`py -3.12 -m pysad.nbcache`。kernel を編集したら試験の前に 1 回走らせる。docs/12 §16.7) | — |
| fitroot.py / miscfn.py / fileops.py | 692 / 599 / 616 | Fit・FindRoot・NIntegrate・FitGaussian、DownhillSimplex、SingularValues、システム・ファイル関数、Parallelize (docs/06) | tffindroot.f, Fit.n, D.n, optimize.n, Matrix.n, FileSystem.n, shared.n |
| corr.py / pkgload.py | 293 / 148 | ResponseMatrix・LinearSolveConditioned・Master と、Packages/*.n (CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・D・rational など) の自動読み込み (docs/06) | pgrmat.f, pgsolvcond.f, pmicad.f, pgmast.f, correction.n ほか |
| ffsextra.py | 490 | FFS の VARIABLES・ATTRIBUTE・MATRIX・SHOW・VARY・RECOVER・REVERSE・STATUS・PRINT・EXECUTE・OUTPUT/APPEND/CLOSE・INPUT/READ (docs/04) | tfvars.f, tfattr.f, tfmat.f, tshow.f, tffsa.f |
| spin.py | 955 | POL (スピン): tradke のスピン部・tchge の srot 回転・spnorm・srequpol/spdepol・Emittance[]/EMIT の PolarizationTime/SpinTune/PolarizationVector (MULT・空洞・SOL まで。`Polarization` は K64 が NaN: B55。docs/16) | photons.f (tradke, tspin), tchge.f, temit.f |
| synchrobeta.py | 714 | SynchroBetaEmittance1 (シンクロ・ベータ結合のエミッタンス。運動量のずれごとの閉軌道と tturne を取り方位角モード展開を解く) | tfemit.f tfemits, temits.f ほか |
| measure.py | 1097 | MEASURE・TRCOD・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE (v2 モードだけ。K64 のソースの読み。gs との比較は無い: docs/03) と tlum / tluma (ルミノシティの推定。解析値との比較のみ) | trackb.f, tracke.f, trcoda.f, ttinit.f, ttstat.f, tbdecoup.f, tlum.f |
| csr.py | 1077 | CSRTrack・CSRHaissin・CSRSetupOY・CSROYMatrix・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRConvert (= tfcsrvec)・CSRInit と実FFT・spline1 (docs/15。ScanTrack・OYStability は未照合) | csrtrack.f, csroy.f, csrtest.f, tcftr.f, spline.f |
| tkcanvas.py | 434 | Tk のキャンバスの組み込み関数 (TkCanvasCreateItem、CanvasClipLine、CanvasSymbolDirect …)。ListPlot・Plot の描画 (docs/21) | tfTkInter_.c, tfcanvasclip.f |
| vecfn.py | 56 | Vector・ReleaseVector (ベクトル評価。Plots.n の PlotInterpolate) | — |
| radfield.py | 158 | RadiationField・RadiationSpectrum。RADLIGHT の軌道記録 (`TrackParticles` が `{結果, 軌道{t,x,y,z}}` を返す) は `ffscmd._RadLight` (v2 モードは 3 次元の軌道で BEND の中も記録。docs/18) | tlinit.f, tlstore, tlresult |
| sadrandom.py | 197 | SAD の乱数 (LCG、Box–Muller、drndsr/tdusrn、ParabolaRandom) | random_plugin_sad.c, tfRandom_.c, drndsr.f |
| sadscript.py / sadlist.py | 5604 / 1468 | SADScript インタプリタ (字句・Pratt 構文・評価、Module/Do/Table/純関数/パターン/REPEAT/UNTIL/IF、単位付き数値、`?FLAG`、fmt_real=autofg) とリスト・文字列関数 | tfetok.f, tfeval.f, tfstk.f iprior, autofg.f, functions.n, Misc.n |
| wake.py | 308 | 追跡へのウェイク場の適用 (WakeFunction、LWAKE/TWAKE。TrackParticles の TRPT のみ) ([13_wake.md](13_wake.md)) | txwake.f, tffssetupwake (tffswake.f), tcav.f / tturn.f の呼び出し箇所 |
| tclbridge.py / sadclass.py | 330 / 424 | Tcl/Tk の接点 (tkinter) と SAD の Class。SAD 自身の Tkinter.n / KBFrame.n を SADScript として解釈して GUI を動かす ([20_tcltk_gui.md](20_tcltk_gui.md)) | tfTkInter_.c, tfTclArg.f, tfXlib_.c, Tkinter.n, KBFrame.n, KBMainFrame.n, TabFrame.n, SADHelp の Class |
| ffscmd.py | 5785 | FFS コマンド層 (tffsa の解釈順)、USE/CELL/INS/CALC/CALC6D/GO/EMIT/TYPE/DISP/SAVE/RESET/RADINT/CHROMATICITY/誤差コマンド (DELX…)、Twiss[]/LINE[]/FFS[]/Element[]/Emittance[]/TrackParticles[]/CalculateOptics/SeedRandom[]、EMIT の結果変数、`FFS$SHOW` | tffsa.f, tftwiss.f, tflag.f, terror.f, tftype.f, tfdisp.f |
| diag_fit.py / diag_trace.py | 245 / 267 | 照合用の診断 (fit の反復表、追跡の途中値) | — |
| ampnf.py・ampnf_coupled.py・ampnf_ext.py・ampnf_k0.py | 1230 / 577 / 309 / 884 | 1 周の写像のべき級数 (TPSA) の Birkhoff 正規形で、振幅・δ に依るチューン (ADTS の J¹・J² の項、クロマティシティの全次数、振幅と δ の交差項) を求める (`ampnf`)。結合した光学 (ROTATE した磁石)・DX/DY/DZ/DROTATE の誤差・4 次の共鳴駆動項 (`ampnf_coupled`)、5・6 次の共鳴の駆動項の表と 1 本の共鳴の島の幅 (`ampnf_ext`)、K0≠0 の偏向磁石 (挿入光源・逆曲げ・MPW) の取り込み (`ampnf_k0`)。docs/25 | K64 にない追加 |
| fourier_fit.py | 320 | MAIN の TRACK の FOURIER 解析 (作用変数の多項式 + 共鳴項の最小二乗フィット。FFS からは呼べない。docs/31) | wjfit.f, wfres.f, wfbin.f, wiord.f |
| quadbody.py | 175 | 四極の本体の積分器の共通部 (K64 の Strang 分割と、v2 モードの方式の切り替え) | tsolqu.f, tsolque.f の bz=0 の枝 |
| radbeta.py | 171 | 放射の β の因子 (β < 1 の粒子: 陽子・イオン・低い γ の放射の補正係数。v2 モード) | K64 にない追加 (厳密な放射のスペクトルから決めた補正) |
| touschek_bunch.py | 211 | Touschek の損失率にバンチ内の δ の分布 (対の重心の δ) を入れる (v2 モードの既定。k64 は入れない) | K64 にない追加 |
| draw_mpl.py | 1155 | FFS の `DRAW` コマンドを matplotlib で実装 (CANVASDRAW.n の構文解析・OpticsPlot.n の点列と図形の作成の移植)。docs/21 | CANVASDRAW.n, OpticsPlot.n |
| mpl_display.py | 1326 | `$DisplayFunction` の `MatplotlibDrawer` (Plot・ListPlot・OpticsPlot・GeometryPlot などの Graphics を matplotlib に描く。環境変数 `PYSAD_DISPLAY=mpl` で切り替え)。docs/21 | Canvas.n の `CanvasDrawer` に代わる実装 |
| dxf_lattice.py | 813 | ラティスの平面図 (と立体の外形) をテキストの DXF で書き出す。docs/22 | drawmagvw3d.n の幾何の考え方を参考にした自前版 |
| run_sad.py (<repo> 直下) |  | `py -3.12 run_sad.py lattice script …` = `gs` に複数ファイルを連結して渡すのと同じ | — |
| **非推奨** ffs.py / matching.py / tracking.py / optics.py | 1599 / 602 / 1239 / 471 | 独自実装 (SAD 未照合。tracking.py は放射係数に疑い)。新規コードは使わない。`__init__.py` の export は旧 API のまま残っている | — |

行数は `wc -l pysad/*.py` の値。`__init__.py` の `from pysad import ...` は旧 API を export しているので、各モジュールを直接 import する (例 `from pysad.parser import load_sad_file`)。**非推奨**のモジュール (`ffs.py`・`matching.py`・`tracking.py`・`optics.py`) は独自実装で SAD と未照合 (`tracking.py` は放射係数に疑い)。新規コードでは使わない。

## 6. pysad で計算できること (領域別)

pysad が扱う領域と、その詳細を書いた文書の対応である。k64 モードの SAD との一致の度合い (桁) は [27_benchmarks.md](27_benchmarks.md)、v2 モードで値が変わる箇所は [28_v2_vs_k64.md](28_v2_vs_k64.md)、できないことは [31_unimplemented_and_unverified.md](31_unimplemented_and_unverified.md) を見る。

| 領域 | 内容 | 文書 |
|---|---|---|
| 光学 (CELL/INS の CALC) | 4D 閉軌道・Twiss・結合・オフモメンタム・不安定ラインの続行・TRPT の加速・MARK OFFSET・幾何 (GX..、GEO、ORG)・RADINT・CALC6D | 07、23 |
| 要素の写像 | 四極のフリンジ・薄い多極・BEND の K0/K1/ロール/F1・ANGLE=0 の BEND・MULT・ROTATE やずれつきの全要素 | 02、07、10 |
| エミッタンス | Emittance[]/EMIT (RFSW × RADCOD)・不安定な光学・INTRA/Touschek・SynchroBetaEmittance | 09、26 |
| マッチング | FIT/FREE/FIX/GO (解析微分・DESCEND・fuzz・2 段目のフィット・INS・フック)、FITP>1 | 08 |
| 追跡・DA | TrackParticles (NORAD/RAD/FLUC・粒子損失)・DynamicApertureSurvey・誤差コマンド (DELX など)・周波数マップ・局所運動量アクセプタンス | 10、11、26 |
| 共鳴・非線形 | 6 極の共鳴駆動項・ADTS・クロマティシティ・周波数マップの共鳴線の予想 | 24、25、26 |
| 集団効果 | WAKE (追跡・光学)・BEAMBEAM・CSR・空間電荷 (v2 モード) | 13、14、15、10 |
| スピン・ソレノイド・放射光 | POL (スピン偏極)・SOL 領域の Twiss/Emittance[]/追跡/幾何・RADLIGHT・RadiationField/Spectrum | 16、17、18 |
| その他の要素 | COORD・TCAVI・UND・曲げつき MULT・MAP・DELCOR | 19、04 |
| 台本・関数 | SADScript インタプリタ・特殊関数・Fit・FindRoot・NIntegrate・SAD の Packages (CorrectOrbit・MakeBump・OptimizeOptics など)・FFS の小さなコマンド | 05、06、04、03 |
| 描画・GUI | ListPlot・Plot・OpticsPlot・Tk の GUI・matplotlib の出口・DXF の CAD 図 | 20、21、22 |

MON / MCOD / TCOD はこの SAD ビルドでは何もしない (式として評価されるだけ) ので pysad も受理のみ。

## 7. 性能と試験

**性能** (2026-11-06 の測定値。再測定していない):

- 1 回の cell calc 0.03 s (不変時) / 0.145 s (変更後)、Twiss[] 0.6 ms、追跡 176 ns/粒子・要素 (numba)。
- 600 s 内に終わらない計算: 誤差入り 6D emit × 多数の実現 (300 実現 221 s、SAD 37 s)、DA の 5000 周 × 969 粒子 (約 50 分、SAD も 300 s 超)。残る改善点は 6D の tturne 1 周 (0.03 s) と DA 系。

**試験とベンチ**:

- テストは [`pysad/tests/`](../pysad/tests/) (`test_*.py`)。2026-11-21 の通しは、収集 1907 件、**1899 通過・8 skip・失敗 0** (その時点で `test_*.py` は 158 ファイルだった。以後ファイルが増えており、現在の件数は未再確認)。v2 の試験は `test_v2_*.py` (73 ファイル・1081 件) で `compat.set_mode('v2')` を明示する。その他は gs との比較なので k64 モードで走る ([`conftest.py`](../pysad/tests/script_tests/conftest.py))。
- 実行:

  ```
  py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/script_tests/test_bend_feeddown.py --ignore=pysad/tests/script_tests/test_coord.py --ignore=pysad/tests/script_tests/test_sext_radiation.py
  QUICK=1 py -3.12 -u pysad/tests/regression_tests/run_regression.py
  ```

  `--ignore` する 3 本 ([`test_bend_feeddown.py`](../pysad/tests/script_tests/test_bend_feeddown.py)・[`test_coord.py`](../pysad/tests/script_tests/test_coord.py)・[`test_sext_radiation.py`](../pysad/tests/script_tests/test_sext_radiation.py)) は sys.exit で終わる旧式のスクリプト。環境変数 `PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>` が無いと GUI 系は skip。`PYSAD_SLOW=1` で遅い試験 (CSR の MakeZL、約 7 分) も走る。GUI 系 (`test_plot`・`test_tclbridge`) は全体実行で数回に 1 回失敗することがある (単独なら通る)。
- 回帰 [`pysad/tests/regression_tests/run_regression.py`](../pysad/tests/regression_tests/run_regression.py) (k64 モード。説明は [`pysad/tests/regression_tests/README.md`](../pysad/tests/regression_tests/README.md)) は 50 ベンチ。`QUICK=1` は遅い 4 件 (7・12t・13a・16a) を SKIP して 46 件を実行する。2026-11-18 の `QUICK=1` の通しは **FAIL 0**。QUICK なしの全 50 ベンチ (205 s) の記録は 2026-11-16 夜 (`local/wake/verify_log.md`)。
- ベンチ: `pysad/tests/regression_tests/bench/` (自動判定の対象外は `bench_extra/`) の `bench_level1〜31*.sad/.result` (gs の出力) と `*_pysad.py/.result`。結果の一覧と一致度は [27_benchmarks.md](27_benchmarks.md)。
- 試験が走らせる SAD スクリプトと gs・pysad の画面表示は [`pysad/tests/script_tests/scripts/`](../pysad/tests/script_tests/scripts/) (`NAME.sad`・`NAME.gs.txt`・`NAME.pysad.txt`)、試験関数との対応は [`pysad/tests/script_tests/OUTPUTS.md`](../pysad/tests/script_tests/OUTPUTS.md)。

## 8. 文書の目次 (全 36 本)

次の表は docs/ の全文書の一覧である。図つきの HTML は、md と同じフォルダに同名の `.html` として置く (`tools/build_html.py` で作る)。

| 群 | 文書 | 内容 |
|---|---|---|
| 入門 | [01_overview.md](01_overview.md) | 本書 |
| 使い方・リファレンス | [02_lattice_main.md](02_lattice_main.md) | SAD の MAIN レベル (ラティス定義) と parser |
| | [03_ffs_commands.md](03_ffs_commands.md) | FFS のコマンドの解釈と ffscmd |
| | [04_ffs_small_commands.md](04_ffs_small_commands.md) | FFS の小さなコマンドと DELCOR |
| | [05_sadscript.md](05_sadscript.md) | SADScript 言語とインタプリタ |
| | [06_sadscript_library.md](06_sadscript_library.md) | 特殊関数・Fit・GaussianCoulomb などの関数と SAD の Packages |
| | [07_calc_optics.md](07_calc_optics.md) | CELL/INS の光学計算 (CALC) |
| | [08_matching.md](08_matching.md) | FIT/FREE/GO のマッチング |
| | [09_emittance.md](09_emittance.md) | Emittance[] と EMIT |
| | [10_tracking.md](10_tracking.md) | TrackParticles の追跡 |
| | [11_dynamic_aperture.md](11_dynamic_aperture.md) | DynamicApertureSurvey |
| | [12_pysad_api.md](12_pysad_api.md) | Python API リファレンス |
| 機能別 | [13_wake.md](13_wake.md) | ウェイク場 (追跡と光学) |
| | [14_beambeam.md](14_beambeam.md) | ビーム・ビーム、BSTRL、光子、ルミノシティ、BBBrem |
| | [15_csr.md](15_csr.md) | コヒーレント放射 (CSR) |
| | [16_spin_pol.md](16_spin_pol.md) | スピンと偏極 |
| | [17_sol.md](17_sol.md) | ソレノイド領域 |
| | [18_radlight.md](18_radlight.md) | RADLIGHT と放射場 |
| | [19_other_elements.md](19_other_elements.md) | COORD・TCAVI・UND・曲げつき MULT・MAP |
| | [20_tcltk_gui.md](20_tcltk_gui.md) | Tcl/Tk の GUI |
| | [21_graphics.md](21_graphics.md) | グラフィック描画 (Plot 系、Tk キャンバス、matplotlib) |
| | [22_lattice_dxf.md](22_lattice_dxf.md) | ラティスの CAD 図 (DXF) |
| 物理と理論 | [23_physics_guide.md](23_physics_guide.md) | K64 の計算の式 (座標・要素・閉軌道・放射積分・エミッタンス・追跡・マッチングなど) |
| | [24_adts_rdt_primer.md](24_adts_rdt_primer.md) | ADTS・RDT・周波数マップ・LMA の入門 (図入り) |
| | [25_resonance_terms_adts.md](25_resonance_terms_adts.md) | 共鳴駆動項・ADTS・クロマティシティ |
| | [26_frequency_map_lma.md](26_frequency_map_lma.md) | 周波数マップと局所運動量アクセプタンス (使い方の章を含む) |
| 検証・既知の問題・未実装 | [27_benchmarks.md](27_benchmarks.md) | ベンチマークと回帰・試験 (検証の入口) |
| | [28_v2_vs_k64.md](28_v2_vs_k64.md) | v2 モードと k64 モードの違い |
| | [29_sad_quirks.md](29_sad_quirks.md) | SAD の癖・仕様 |
| | [30_sad_bug_report.md](30_sad_bug_report.md) | K64 のバグ目録 |
| | [31_unimplemented_and_unverified.md](31_unimplemented_and_unverified.md) | 未対応・未確認の唯一の一覧 |
| 開発者向け | [32_source_analysis.md](32_source_analysis.md) | SAD ソースの構造、コールスタック、デッドコード |
| | [33_fortran_status.md](33_fortran_status.md) | SAD のソース (ファイル・ルーチン・機能) 対 pysad の翻訳状況 |
| | [34_sad_help_index.md](34_sad_help_index.md) | SADHelp の項目一覧 (機械抽出) |
| | [35_pysad_vs_sadhelp.md](35_pysad_vs_sadhelp.md) | SADHelp の項目と pysad の対応 |
| | [36_sadhelp_purpose_policy.md](36_sadhelp_purpose_policy.md) | 目的別の分類と翻訳方針 |

`docs/sadhelp_summary.html` は SADHelp の項目と pysad の対応のまとめ (34〜36 の内容)。経緯・作業の記録は `docs/history/` にある (本文は現在の状態だけを書く)。
