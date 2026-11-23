# 33. SAD のソース (ファイル・ルーチン・機能) と pysad の翻訳状況

この文書で分かること: SAD の Fortran ファイルごとに、ルーチン (ファイル:行)・役割・対応する pysad の関数・翻訳の状態・照合の桁 (§1〜§8)、機能別の対応状況 (§9)、pysad のコードに名前が出てこない src のファイル (§10)。
前提: SAD (K64 版) のソース構造 (docs/32)。互換モード (v2 と k64) の意味 (docs/28)。
読む順: §0 の凡例、探したいファイルの節 (§1 追跡、§2 線形光学、…)、機能名から探すなら §9。未対応・未確認の一覧は docs/31。

SAD は K64 版、ソースは `<SAD>\src`。表の「状態」は、特記がなければ **k64 モード** (K64 (gs) と同じ動作。gs との照合はこのモード) のもの。v2 モード (既定) で変わった動作は「v2 モード:」と書く (一覧は [28_v2_vs_k64.md](28_v2_vs_k64.md))。v2 の値は gs ではなく解析式・独立な数値計算で確かめたもので、gs と比較できないものがある (MEASURE・TRCOD、空間電荷の SPAC、領域内の SOL の要素など)。
「デッドコード・ビルド対象外」の判定は mk/sad.obj.mk の有効リストに無い、または call 0 件のもの。行数の多いデッドコードは節ごとに 1 行にまとめる。
`<SAD>\src` の *.f 474 本のうち、ファイル名・ルーチン名が `pysad/*.py` に 1 度も出てこないものの一覧は §10 に付けた (名前の grep だけの確認で、名前を変えて直訳した場合は漏れる)。
根拠: `01_overview.md`、`27_benchmarks.md`、`history/00_RESUME.md`、`history/work_*.md`、`pysad/*.py` の実コード、`history/work_unimpl.md`、`history/v2_changes_*.md`。
ルーチンの欄の `名前:数字` は SAD ソースの行番号 (2026-10 時点)。アーキテクチャとコールスタックは docs/32。

## 0. 凡例と現在の pysad モジュール

表の「状態」欄の言葉は次の 4 つである。

| 状態 | 意味 |
|---|---|
| 翻訳済み・照合済み (桁) | pysad に直訳があり、SAD の出力と比較済み。桁は最大誤差 (`27_benchmarks.md` の値) |
| 翻訳済み・未照合 | コードはあるが SAD との数値比較が未実施、または再現が不完全 (括弧内に理由) |
| 未翻訳 | pysad にコードがない (NotImplementedError で止まるものを含む) |
| 対象外 | デッドコード、非リンク、GUI/グラフィックス、EPICS など |

現行の対応モジュール: parser / sadelem / lattice / elkern / perfkern / perfkern_bend / turnkern / emit6d / teigen / match6 / match6d / tsvd / track6 / track6n / track6n_ext / dapert / sadrandom / sadscript / sadlist / sadclass / ffscmd / ffsextra / calc6d / touschek / geo / sgeo / emit_extra / radint / maintrack / trptacc / optics_fn、
要素・機能別のモジュール: coord / tcav / und / multa / multae / solmap / solemit / spin / csr / synchrobeta / wake / wakeopt / beambeam / gausscoul / bbbrem / measure / radfield / specfn / fitroot / miscfn / vecfn / fileops / corr / pkgload / tclbridge / tkcanvas、
v2 のモジュール: compat (互換モード) / phsrot / radtaper / spch (空間電荷)。`measure` は v2 で MEASURE・TRCOD ほかを実装している。
(`pysad/*.py` は 71 本・計 57124 行。`tracking.py` の ttcav 等への言及は旧 API の名残で、実装は `tcav.py`)
`ffs.py`, `matching.py`, `tracking.py`, `optics.py` は非推奨の旧 API (SAD 未照合の独自実装) で、表には数えない。
calc6d / match6d / touschek / geo / track6n_ext / emit_extra はいずれも `ffscmd.py` から呼ばれる (CALC6D、解析微分、INTRA/TouschekLifetime、LINE["GX"..]/GEO/ORG、MULT/SOL の numba 追跡、Emittance[] のオプション)。

試験の現状 (`PYSAD_SAD_PACKAGES` つき): `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/script_tests/test_bend_feeddown.py --ignore=pysad/tests/script_tests/test_coord.py --ignore=pysad/tests/script_tests/test_sext_radiation.py` を実行する (件数は docs/27 §5.1 の最新の通し)。互換モードは k64 が既定 ([`pysad/tests/script_tests/conftest.py`](../pysad/tests/script_tests/conftest.py))。v2 の試験 (`test_v2_*.py`) は `compat.set_mode('v2')` を明示する (本数・件数は未確認)。回帰 (`QUICK=1 pysad/tests/regression_tests/run_regression.py`) の結果は docs/27 §5.1 (k64 モード)。

照合桁の目安 (CALC 光学・Emittance[]・追跡・マッチング・DA) は [history/33_fortran_status_validation.md](history/33_fortran_status_validation.md)。表の「照合」欄には個別の桁を書く。

## 1. 粒子トラッキング

追跡 (TrackParticles) のほぼ全要素が翻訳済みで、SOL 領域内の一部と PHSROT・空間電荷はモードによる (次の表)。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| tftrack.f (tftrack:8)、tturn.f (tturn:36 / tturn1:153。要素ごとの分岐は select case L286-552) | 1 周のディスパッチ、tapert | track6.py `_tturn1`, track6n.py `turn`。TrackParticles 本体は `track6.track_particles`、`track6n.track_particles_fast` | 翻訳済み・照合済み | ≤5e-11 (L6〜12, 16) |
| tdrift.f | ドリフト (tdrift_free, tdrift_solenoid) | track6 `_tdrift_free`, `_tdrift_solenoid` | 翻訳済み・照合済み | ≤5e-11 (L6, L9) |
| tquad.f (tquad、tthin:118、ttfrin。tthin 同居) | 四極・薄い多極、ttfrin | track6 `_tquad`, `_tthin_gen`, `_ttfrin*` (n 版 `k_quad`, `k_thin`) | 翻訳済み・照合済み | ≤5e-11 (L6, L8) |
| tbend.f (tbend:531、tbendcore、tbshift、tbrot)、tbendi.f | BEND、K1≠0 の複合 (bendib) | track6 `_tbendcore`, `_tbendi` (ほか `_tbshift`, `_tbrot`。n 版 `tbendcore1` ほか、K1≠0 は `k_bendi`) | 翻訳済み・照合済み | ≤5e-11、z のみ 1.7e-10 (L7) |
| tbfrie.f の tbedge, tqfrie, tqlfre, tqente | フリンジ | track6 `_tsetfringep`, `_tqlfr`; lattice.py | 翻訳済み・照合済み | ≤5e-11 (L8)、emit/CALC 側も一致 |
| tcav.f (tcav、tcavfrin) | RF 空洞 (tcavfrin 含む) | track6 `_tcav` (n 版 `k_cav`) | 翻訳済み・照合済み | ≤5e-11 (L6, L8-2) |
| tmulti.f (tmulti:37、tmulk) | MULT (ANGLE=0, K0..K21, ずれ・回転) | track6 `_tmulti_ens`、numba 版 track6n_ext.py | 翻訳済み・照合済み (track6 経由) | ≤5.9e-11 (L11)。ext は track6 と 4e-17、seed 一致 |
| tsolqu.f (tsolqu、tsolqum、tsolqur)、tsol.f (tsol:510) | ソレノイド領域 (bz 分岐、tsolqur) | track6 `_tsol_region`, `_tsolqu_bz_*`、track6n_ext.py (ほか `_tsolqur_ens`。n 版 `track6n.tsolqu1`、`track6n_ext.tsolqu_bz_all` (track6 と 4e-17)、`k_sol_*` (4e-17)) | 翻訳済み・照合済み (領域内の DRIFT・QUAD・SOL)。v2 モード: 領域内の BEND・MULT・CAVI・APERT、ずれ・回転 (trots)・F1 (tradks)、領域内の空洞の TRPT 加速・wake・AUTOPHI も実装 (ローレンツ力の数値積分との一致で、gs との比較ではない。docs/17)。v2 は tsolqur の bzp (B08)・tsolqum の dv 項 (B14) も直している。K64 は SOL (BOUND≠0) の長さを周長に入れない (B63) ので、追跡の `circ` も `emit6d._emit_circ` を使う | ≤1.5e-10 (L9)。EPS≠1 は制限あり |
| photons.f (tradk:1267) の tphoton | 放射キック (平均・FLUC) | track6 `_tradk`, sadrandom.py (ほか `track6n.tradk_core1`、`sadrandom.drndsr/tdusrn`) | 翻訳済み・照合済み | ≤4.9e-11、乱数状態一致 (L10) |
| tapert.f (tapert、tapert1) | 損失判定・APERT 要素 (tapert1) | track6n.py `k_apert_gen`・`k_apert`・`k_apert_el`、track6.py `_apert_pass`、dapert.py、`TrackRun.l_track0` | 翻訳済み・照合済み (APERT 要素: 楕円・矩形・DX DY・ROTATE・全部 0・複数周・JDPX JDPY DP が gs と一致、`tests/test_apert.py`。SOL 領域内の APERT は k64 モードでは未対応 (v2 モードは領域内の APERT で落ちる)。TRPT・RAD つきは未確認。JDPX の乱数は粒子 12 個以上で gs と不一致) | DA の地図一致 (L16A)。ただし損失周の乱数消費は不一致 (RAD+FLUC)。K64 は落とす印に `kptbl(i,4)=l_track` を書くので、セッション最初の TrackParticles (l_track=0) の最初の tapert では最初から LOSSAMPL の外の粒子が落ちない。k64 モードはこの癖を再現し gs (1・2 回目) と全桁一致 (B33 の出力座標の差は出ない、docs/29 A52)。v2 モード: 最初の tapert で落とす。LOSSMAP (戻りに落ちた要素・周の 2 行) は両モードで gs の例と一致 (docs/10 §9.7) |
| tsteer.f | BEND ANGLE=0, K0≠0, L≠0 | track6.py `_tsteer_ens` (F1/FB1/FB2 フリンジ付きは NotImplementedError) | 翻訳済み・照合済み | bendk1_* (history/work_bend_k0.md, history/work_trpt_acc.md 10 節) |
| tmulta.f (tmulta, tmultae) | MULT ANGLE≠0 | multa.py (追跡)、multae.py (Twiss・Emittance[]) | 翻訳済み・照合済み (docs/19) | 追跡は相対 0、Twiss 1e-13、Emittance[] 1e-11。BZ・RF・F1・SK0 は未対応 |
| ttcav.f, ttcave.f | TCAVI (クラブ空洞) の追跡と 6 次元写像 | pysad/tcav.py | 翻訳済み・照合済み (追跡 0、Emittance[] は L=0 のみ gs と一致。L≠0 は gs が壊れる) (ttcave の RFSW ON の CALC は未対応) | docs/19_other_elements.md |
| tcoord.f, tcoorde.f, qcoord.f | COORD 要素 | pysad/coord.py | 翻訳済み・照合済み (追跡 0、Twiss・幾何 1e-11、Emittance[] 5e-10) | docs/19_other_elements.md |
| temap.f (temap/temape/qemap)、tins.f (tins) | MAP 要素 (ExternalMap) | elements.Map、ffscmd._map_*、emit6d._tmap、track6n.K_MAP | TRACK・EMIT・OPTICS を翻訳し gs と照合済み (GEO・CompiledMap は未対応)。INS (tins/tinse) は `lattice.ins_twiss` で扱い、追跡・Emittance[] も gs と一致を確認した (2026-11-16。追加の実装は不要) | test_map.py |
| undulator.f | UND (追跡) | pysad/und.py | k64 モード: 翻訳済み・照合済み (追跡 0、Twiss・Emittance[] は恒等で一致)。K64 は磁場を設定できない (B37)。v2 モード: `FBX`・`FBY` の磁場つき (追跡・Twiss・Emittance[]・EMIT。dpds の二重掛け B78 も直している)。磁場なしでも光学にドリフトとして入る。track6 のスピン追跡も通る。解析式で確認し gs との比較は無い | docs/19_other_elements.md |
| bbstrhl.f (beambeam、bbinit)、bbtool.f、matrix.f (包絡行列) | BEAMBEAM | pysad/beambeam.py | 追跡 (track6n) を翻訳。BSTRL と光子 (PHOTONS、BEAMBEAM の分だけ)、LuminosityMonitor も対応。光学は恒等。傾き tilt の K64 の癖 (B34) は k64 モードでは直訳、v2 モードは (cos t, sin t) の回転に直している (残る差は docs/14) | docs/14_beambeam.md |
| photons.f の sprot、tradk のスピン部 | TrackParticles のスピン追跡 (POL かつ RAD ON のときだけ) | track6.py (`sprot`、`bsi`) | 翻訳済み・照合済み (gs と 4e-12 以下: 平面リング、ずれ・回転した四極・偏向・6 極、MULT、K1 のある BEND、SEXT・薄い MULT、SOL 内のドリフト)。k64 モードは損失 (flag) の処理なし (v2 モードは落とす。領域内の BEND・MULT・CAVI・ずれ・F1、UND、PHSROT のスピンも通る。docs/16)。SOL 内の四極のスピン追跡は gs の食い違い (docs/17_sol.md §7) のため未確認 | test_spin_track.py、docs/16_spin_pol.md |
| tlinit.f (tlinit/tlstore/tlresult) | RADLIGHT の軌道の記録 | ffscmd `_RadLight` | k64 モード: 翻訳済み・照合済み (TrackParticles の戻り `{結果, 軌道{t,x,y,z}}`、gs と差 0。K64 は BEND の中を記録せず横方向の座標も入れない)。v2 モード: 3 次元の軌道、BEND は円弧で中も `RADLIGHTDIV` 点記録 (B65。Liénard-Wiechert のピークと 0.25%、gs との比較は無い) | test_radlight.py、docs/18_radlight.md |
| tbbbrem.f | BBBrem1 関数 | bbbrem.py | 翻訳済み・gs と一致 (docs/06) | test_bbbrem.py |
| csrtrack.f, csroy.f, csrtest.f | CSR 関数 (CSR.n が使う 7 関数: CSRTrack CSRHaissin CSRSetupOY CSROYMatrix CSRMatrix CSRConvert CSRInit) | pysad/csr.py | 翻訳済み。CSRTrack・CSRHaissin は gs と 11 桁、CSRSetupOY は J・ω 1e-9、CSRMatrix (mode 22・33・23・32・110・220・330 と mode 1000 = csrhelmat)・CSRConvert (= K64 の tfcsrvec)・CSRInit は gs と相対 1e-10〜1e-13。CSROYMatrix は gs が wrongtype を返す (B60) ので未照合。`Packages/CSR.n` の MakeZL・SkinDepth・PathDifference などは SADScript のまま動かして gs と一致 (MakeZL は約 7 分、相対 1e-10)。ScanTrack・OYStability などの長い流れは未照合。`CSR.n` 全体は描画 (ListPlot) が要るので読み込まない | docs/15_csr.md、test_csr.py |
| track.f (track:1。twig:332 は空スタブ)、tracka.f、ActTra.f (ActTra:257)、ttinit.f | MAIN の TRACK コマンド・初期分布 | ffscmd `TrackParticles[]` (SADScript 関数経由)、MAIN の TRACK コマンドは maintrack.py | 翻訳済み・照合済み | TrackParticles は L6〜12。MAIN TRACK は SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 (SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV)。PHSPLOTS の描画は無視、FOURIE ON は未対応 |
| trackb.f, tracke.f, trcoda.f, ttinit.f, ttstat.f, tbdecoup.f, nalign.f, tfojit.f, tfgetm.f | MEASURE・TRCOD・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE、MAIN の `TRACK TURNS=-1` | `measure.py` | **v2 モードで翻訳** (K64 のソースの読み。TRCOD の FIXCOD/MVCOD/STNCOD/ALIGN も)。k64 モードは語を受理するのみ (無視)。**gs との比較は無い** (gs で TRCOD は SIGSEGV、MEASURE は出力なし: docs/30 B42)。解析式と pysad の独立な経路で確認。ttinit の読みの誤り B41 は v2 で直している。DRAW の中の粒子追跡 (tdrwdt.f) は K64 でも動かない死んだコードで翻訳しない (B105) | test_v2_C_measure.py、test_v2_C2_measure.py、docs/03 |
| tturn.f・tcod.f・tturne.f の radtaper | RADTAPER (磁石の強さを閉軌道の局所の運動量に合わせる) | `radtaper.py` | 翻訳済み・gs と照合 (両モード。CODCONV=1e-4 の gs と x 6 桁・U0 10 桁・チューン全桁)。K64 の tcod の反復の直訳ではなく同じ不動点を外側の反復で求める。4 次元の CALC のテーパーは K64 も無い | test_v2_C_track.py、test_v2_C2_optics.py |
| spkick.f、tspch_.c・spch.c・psn.c、tintrb.f twspac | 空間電荷 SPAC・PSPAC・WSPAC | `spch.py` | **v2 モードで翻訳** (k64 モードは入れない)。K64 の SPAC は力が 1/γ 倍に弱い (B89) ほか B90〜B98 の誤りがあり、v2 は直している。PSPAC は gs と 3〜4 桁、WSPAC は gs と 1% (Δν)、SPAC は gs の誤りのため解析式 (包絡方程式・Laslett の式) と比較。Emittance[]・EMIT・CALC6D・CALC (RFSW ON) の WSPAC は実装 (Emittance[]・EMIT は k64 も gs と全桁。K64 の最後のスライスのレンズが 8 倍強い誤り: B127)。追跡の WSPAC のビーム行列は空間電荷なしのまま | test_v2_C_spch.py、test_v2_C2_spch.py、docs/10 §9.7 |
| trackd.f (trackd:37) / tfdapert.f (tfdapert:1) | DynamicAperture 本体 | dapert.py | 翻訳済み・照合済み | 地図 153/153、Score 一致 (L16A、NORAD) |
| tquads, tquase, tbrad, tmulbs, ttdr, bb, hg 等 | Obsolete・非リンク | — | 対象外 | — |

### 1.1 tturn1 の要素ごとの分岐 (tturn.f select case L286-552)

追跡の 1 周のディスパッチ表で、✓ は翻訳済み、△ は一部、✗ は未実装を表す。

| コード | case | 呼出 | pysad |
|---|---|---|---|
| 1 | icDRFT | tdrift_free:294 | ✓ `_tdrift_free` (≤5e-11) |
| 2 | icBEND | tbend:327 | ✓ `_tbend` (K1≠0 `_tbendi`)、DX/DY/DROTATE も可 |
| 4 | icQUAD | tquad:355 → tsolqu | ✓ `_tquad` / `k_quad` |
| 6,8,10,12 | SEXT/OCTU/DECA/DODECA | tthin:380 | ✓ `_tthin_gen` / `k_thin` (SEXT 厚肉放射も一致) |
| 22 | icMULT | tmulti:395 | △ `_tmulti_ens` / `k_mult` (ANGLE=0。RF は `mult_full`)。ANGLE≠0 は ✓ `multa.py` (tmulta)。SK0・F1・TAPER・VOLT は未 |
| 20 | icSOL | tsol:410 | △ `_tsol_region` / `k_sol_*` (DRIFT/QUAD を含む領域のみ。内部 BEND/MULT/CAVI、DX/DY/DZ/CHI、F1 は未。領域内の BEND は gs が偏向をほぼ無視して合わず取り下げ: B69) |
| 31 | icCAVI | tcav:443 | ✓ `_tcav` / `k_cav` |
| 41,42 | MARK, MONI | — | ✓ (無作用) |
| 43 | icAPRT | tapert1:418 | ✓ `k_apert_el` (SOL 領域内の tapert1 は未翻訳) |
| 32 | icTCAV | ttcav:480 | ✓ `tcav.apply_tcav` (RANKICK・RANPHASE は未) |
| 33 | icMAP | temap:490 | ✓ `elements.Map`・`ffscmd._map_*` (ExternalMap の TRACK) |
| 34 | icINS | tins:502 | ✓ (gs と一致を確認。追加の実装は不要) |
| 35 | icCOORD | tcoord:506 | ✓ `coord.py` |
| 36 | icBEAM | beambeam:498 | ✓ `beambeam.py` |
| 37 | icPROT | phsrot:517 | ✗ (K64 で動作しない: B38。NotImplementedError) |
| 38 | icSPCH | tpspac:525 | ✗ |
| 18 | icUND | undulator:537 | ✓ `und.py` (磁場なしのドリフト: B37) |
| 19 | icWIG | twig:541 (空スタブ) | — (SAD も無作用) |
| 21 | icST | abort ("Use BEND with ANGLE=0") | — |
| — | wake / wspac | txwake:266 / twspac:274 | wake ✓ `wake.py`。wspac (空間電荷) ✗ |

### 1.2 tturne1 の要素ごとの分岐 (tturne.f select case L470-675)

Emittance[] 用の 6 次元 1 周のディスパッチ表。

| コード | case | 呼出 | pysad |
|---|---|---|---|
| 1 | icDRFT | tdrife0:475 | ✓ `emit6d._tdrife0` |
| 2 | icBEND | tbende:533 | ✓ `_tbende` (`_tbendecorr`, `_tbendebody`) |
| 4 | icQUAD | tquade:569 | ✓ `_tquade` |
| 6,8,10,12 | tthine:592 | tthine | ✓ (6 極内のドリフトの dv 項を含む) |
| 20 | icSOL | tsole:597 | △ `solemit.py` / `emit_extra._tturne1` の SOL 領域 (領域内は DRIFT・QUAD。POL のスピンも通る) |
| 21 | icST | forcesf:602 | — |
| 22 | icMULT | tmulte:621 | ✓ `_tmulte` (ANGLE=0)、`multae.py` (ANGLE≠0。BZ・RF・F1・SK0 は未) |
| 31 | icCAVI | tcave:635 | ✓ `_tcave`, `_tcavfrie_map` |
| 32 | icTCAV | ttcave:648 | ✓ `emit6d._ttcave` (L≠0 は gs が壊れる: B36) |
| 33 | icMAP | temape:661 / qemap:657 | ✓ `elements.Map` の EMIT・OPTICS フック |
| 34 | icINS | tinse:665 | ✓ (`lattice.ins_twiss`。gs と一致) |
| 35 | icCOORD | tcoorde:668 | ✓ `coord.py` |
| — | default | 無作用 (WIG, MARK)。UND・BEAMBEAM は ✓ (`und.py`・`beambeam.py`) | ✓ |

## 2. 線形光学

CALC・Emittance[] 用の線形・6 次元の写像は、MULT ANGLE≠0 と SOL 領域の一部を除いて翻訳済み・照合済みである。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| qtwiss.f (qtwiss1:14、qcod:810、qquad:1737、qmat2twiss) | Twiss 伝搬 | lattice.py `twiss_table` (`_transport_twiss`: qmat2twiss は独自の非結合の標準式) | 翻訳済み・照合済み | ≤1e-10 (L1,12,13)、結合・誤差経由 ≤3e-7 (L3, L5)。結合時の qmat2twiss は部分 |
| qtwiss.f の qcod:810 / tcod.f | 4D 閉軌道 | lattice.py `_qcod_periodic`, `find_cod` | 翻訳済み・照合済み | 9〜10 桁 (L3B) |
| qcell.f (qcell1:27、qcell6d:404) | CELL 周期解・INS | lattice.py `twiss`, `ins_twiss` (A/B/分散は直訳、qmdiag は numpy 固有値分解) | 翻訳済み・照合済み | CELL ≤1e-10。INS: 軌道 ≤1.3e-11、Twiss ≤2.6e-7 (L5)。初期軌道オフセット付き INS はチューン不一致 |
| qmdiag.f | 4x4 対角化 | lattice.py (numpy 固有値分解) | 翻訳済み・照合済み | ≤1e-10 |
| tquade.f、tsolque.f (tsolque:313)、tbende.f (tbende:351、qbend:627、tbedge:88)、tdrife.f、tmulte.f | 6D 転送行列 (emit 用) | emit6d.py。CALC 側の COD 写像は `lattice._tquade_cod`, `_tsolque_cod`, `_tqfrie_cod`, `_tqlfre_cod`, `_tqente_cod`, `_bend_cod_update`, `_tbedge_cod`, `_tblfre_cod`, `_mult_sequence`、emit 側は `emit6d._tquade`, `_tbende` (`_tbendecorr`, `_tbendebody`), `_tmulte` | 翻訳済み・照合済み。BEND の K1≠0 は ≤3e-9。ANGLE=0 の BEND の F1 分岐は未実装。MULT ANGLE≠0 (tmultae) は `multae.py` (Twiss 1e-13、Emittance[] 1e-11。BZ・RF・F1・SK0・FB1/FB2 の Twiss は未対応。docs/19)。tmulke の bsir0 の更新も写している (B64) | Emittance[] ≤1.3e-11 (L16B〜E) |
| tthine.f (tthine) | 薄レンズ行列 | emit6d.py (CALC 側 `lattice._tthine_cod`, `_thin_offset_map`) | 翻訳済み・照合済み。6 極内のドリフトの dv 項を含む | ≤1.3e-11 |
| tcave.f (tcave、tcavfrie) | 空洞行列 | emit6d.py (`_tcavfrie_map`) | 翻訳済み・照合済み | ≤1.1e-7 (L4、RADCOD ON) |
| tsole.f (tsole1)、tsolque.f、tqente.f (SOL 領域) | ソレノイド領域の写像・放射 | solmap.py (Twiss)、solemit.py (Emittance[]・EMIT) | 追跡カーネルの数値微分で翻訳。Twiss は領域内の β・α・分散・結合 1e-11、Emittance[] は εx 2e-12〜4e-9・εy 4e-10〜2e-8 で gs と照合済み。未対応 (両モード): CALC6D・Touschek。k64 モードで未対応: 領域内の MULT・BEND・CAVI、SOL の DX DY DZ CHI・F1 (Twiss・Emittance[]・追跡)。v2 モードはこれらを実装 (ローレンツ力の積分との一致。gs との比較は無い。Emittance[] の POL の量は領域内だけ K64 の式)。数値微分なので解析行列の gs と 1e-10 程度の差が残る | test_sol_optics.py、test_sol_emit.py、docs/17_sol.md |
| tchge.f, qchg.f | ずれ・回転座標変換 | lattice.py, emit6d.py、match6d.py `qchg` | 翻訳済み・照合済み | ≤1e-8 (L3B, L3C) |
| tturne.f (tturne:1、tturne1:334。要素ごとの分岐は select case L470-675) | emit 用 1 周ディスパッチ | emit6d.py `_turn` (ほか `tturne`、`emit_extra._turn_record`) | 翻訳済み・照合済み | ≤1.3e-11 |
| qdcell.f (qdcell、qdtwis、qdquad、qdbend、qdmdia、tdfun)、tffsqu | 解析微分 | match6d.py (`tffsqu`, `qdcell`, `qdtwis`, `qdquad`, `qdbend`, `qdmdia`。tdfun (残差 df) は `match6._tdfun1`, `_trdfun`) | 翻訳済み・照合済み | マッチングの反復列・K1 を 12 桁再現 (L22A)。数値微分とは 2e-8〜4e-7 |
| qcell.f の qcell6d:404 / temit plot 経路 | CALC6D | calc6d.py (`tfetwiss`, `etwiss2ri`, `_closed_orbit`) | 翻訳済み・照合済み | 1.7e-10 (L18A/B/C) |
| qins.f, tins.f, tinse.f | INS 要素 (領域の入口・出口を担う) | lattice.py の INS の扱い (`ins_twiss`) | ファイル名は pysad に出てこないが、INS を含む 3 セルのリング (BEND E1・E2 付き) の Emittance[] のチューン 1e-15、TrackParticles 1e-13 が gs と一致 (RFSW・RAD 込み)。追加の実装は不要と判断 (2026-11-16) | docs/history/work_unimpl.md |
| qquad, qthin, trad, twelm, pfcoup 等 | 非リンク・デッド | — | 対象外 | — |

## 3. エミッタンス

Emittance[]・RADINT・Touschek は翻訳済みで、SynchroBetaEmittance と POL の一部だけが一部照合である。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| temit.f (temit:78、tecalc:235、synchm) | 6D 平衡ビーム行列 | emit6d.py `emittance6d` (`tecalc`, `_tsolva`, `_normal_form`。teigenc + tsolva/tsvd) | 翻訳済み・照合済み | 全量 ≤1.3e-11 (L16)、RADCOD ON 1e-7 (L4) |
| INTGRL.f (INTGRL:11) | RADINT (I1〜I5, 結合, 色収差) | radint.py (FFS RADINT)、lattice `_radiation_integrals_bend` (BEND のみの別実装) | k64 モード: 翻訳済み・照合済み (K64 の癖: 最後の BEND の L/ANGLE を全 BEND に使う、6 極が直前の DL を使う、定数が単精度、CAVI は最後の値)。v2 モード: その BEND 自身の L/ANGLE、6 極は自分の L と K2、倍精度、VOLT は総和 (B06・B07・B12・B13)、量子寿命 (B72)、線形クロマティシティは dν/dδ (B73)、MULT・SOL が入らないときは表示 (B74)。RADINT の Jx が正しく、K64 の Emittance[] の Jx は端の角で誤る (B75。v2 は Emittance[] を直している)。MomentumCompaction の誤り (B77) も v2 で直している。MULT・SOL・BEND の K0 は積分に入らない。結合のあるリングの V 面 (I4V) は未確認。RF 関係 (INTGRL.f:249-273) は未実装 (モード別は未確認) | 7 桁 (L1、L29A)。v2: 解析式・自前の積分と 2e-5〜1e-13 (test_v2_A_radint.py) |
| temit.f の INTRA、tintrb.f (touckl、touckm、toucke)、tintraconv | IBS 付き emit / Touschek | touschek.py `touschek_table`, `touschek_lifetime` (`_tintrb`, `_tturne_intra`) | 翻訳済み・照合済み | 15 ケース ≤7e-12 (L17A)、RADCOD ON 1.7e-8 (L17B) |
| tfemit.f, emit.n | Emittance[] のオプション (OneTurnInformation, Matrix) | emit_extra.py | 翻訳済み・照合済み | 2.1e-10 (L21A)。SOL を含むと使えず `run_emit` の値に落ちる (Stable=0、HarmonicNumber=NaN)。軌道のある QUAD で Python 版 `_tquade` が 1e-5 ずれる (history/work_unimpl.md 不具合 1) |
| temitf.f | EMIT 報告 (temitprint ほか) | ffscmd (EMIT 報告) | 翻訳済み・照合済み (報告の行を gs と比較: bench_level26a、29c) | docs/09 |
| tfemit.f の tfemits、temits.f (temits1, tefsetup, tesolvd, tevdif, tesumb など) | SynchroBetaEmittance | synchrobeta.py、`ffscmd.fn_synchrobeta` | 翻訳済み・一部照合 (12 セルのリング + CAVI で emitx が相対 4e-7〜1e-5、emity は 1e-18 の雑音で符号も違う。24 セルの平面リング + CAVI VOLT=2e6 の 4 点で emitx 1.4e-7〜1.4e-6 (local/wake/sb11.sad)。結合のあるラティスは gs が SIGSEGV で比較不能、不安定なチューン、emiout、RFSW の組み合わせは未確認) | test_synchrobeta.py (許容 3e-5) |
| calpol、spnorm、srequpol、spdepol (photons.f の tradke のスピン部、tchge.f の srot) | POL (偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列) | spin.py、emit6d / emit_extra / elkern の spin 分岐 | 翻訳済み・照合済み (平面リングの PolarizationTime は全桁、縦ずれ・傾き・6 極のリングは 13 桁。MULT・空洞・SOL (ドリフトのある領域) まで)。`Polarization` は K64 が NaN (B55) で K64 の式のまま。SOL 領域内の MULT・BEND・CAVI のスピンは未確認 | test_spin.py、docs/16_spin_pol.md |

## 4. 放射

放射 (追跡・Emittance[] 用) は翻訳済み・照合済みである。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| photons.f の tradk:1267 / tradkn / tradkfn | 追跡の放射 | track6 `_tradk`, sadrandom.py | 翻訳済み・照合済み | ≤4.9e-11 (L10) |
| photons.f の tradke:1362 | emit 用放射 | emit6d.py (各 `_t***e` の krad 引数) | 翻訳済み・照合済み | ≤1.3e-11 |
| synradcl.f srfuncy | 光子エネルギー逆 CDF | track6 に内蔵 | 翻訳済み・照合済み | L10 の乱数状態一致 |
| tsdrad, tserad, tstrad, tqrad, trad | Obsolete | — | 対象外 | — |

## 5. COD 補正・FFS のコマンド類

COD と FFS の小さなコマンドは翻訳済みで、MON/MCOD/TCOD は SAD も何もしない。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| tcod.f (tcod:31) | 6D/4D COD、RADCOD | lattice.py, emit6d.py (6D Newton)。`calc6d._closed_orbit`、`emit6d.emittance6d` 内の割線法 (RADCOD OFF は tcodz 相当) | 翻訳済み・照合済み | 9〜10 桁 (L7 の δ, trf0 一致) |
| tffsmatch.f の変数の符号・重み (UNIPOL・RELW・BIPOL) | FIT のフラグ | ffscmd `matcher_` が match6 の bipol・absweit につなぐ | 翻訳済み・照合済み (UNIPOL は変数の符号を変えない、RELW は重みが変わる。fit1.sad の 3 通りで gs と全桁一致)。BIPOL の効果は FIT だけ。GAUSS・CTIME・CONV・EMIOUT はソースで効果を確認 (GAUSS・CTIME・CONV は MAIN の TRACK・MEASURE、CPU 時間の表示、収束の表示。EMIOUT は EMIT 報告の追加行で、pysad は spin の行だけ対応) | test_item4.py |
| tfdisp.f (DISP のモード GA・Z・P・A・B・D と REFERENCE・DREFERENCE) | 光学関数の表 | ffscmd `_disp_show` | 翻訳済み・照合済み (gs と見出し・桁が一致。DISP B の REF は gs も未実装。雑音 1e-16 の列は 0) |  |
| tffsa.f tfsetref、tftwiss.f の R・D 接頭辞 | REFERENCE、Twiss["RBX"/"DBX"/"RALL"/"DALL"] | ffscmd `_ref_rows`、`twiss_value`、`fn_twiss` | 翻訳済み・照合済み (全桁) |  |
| tfinit.f tffsrenumber/tfinimult | RENUMBER、LINE["MULT"] | ffscmd `_renumber_cmd`、`line_value` | 翻訳済み・照合済み (INDX 属性の初期値は未対応) |  |
| tffsa.f の USE/CALC/GO の EXPAND・NOEXPAND | 成分ごとの値を残す | ffscmd `use`、`calc`、`go`、`spread_pending` | 翻訳済み・照合済み (NOKEEPEXP は未対応) |  |
| DUMP コマンド | 要素の DK・DDK などの表 | ffscmd `_gwd` | 翻訳済み・照合済み (パターンは 20 文字まで、DDK は K0 のある要素だけ。誤差入りリングで値は gs と一致。DK・DDK 列の書式 (Fortran の 1p g12.5: 0 → `0.0000`、小さい値 → `4.22829E-04`) は gs と同じ) | test_ffsextra.py `test_dump_format_matches_gs` |
| MATRIX コマンド | 2 点間の転送行列 | ffsextra | 翻訳済み・照合済み (名前だけ (gs は NAME.1 として扱う) でも位置に指定でき、ラティスに無い名前を gs が読み飛ばす (to=`$$$`) のも再現) | test_ffsextra.py (mt1・mt2.sad) |
| LDSolve.n、行列演算 (Dot・Inverse) | Optics.n・LDSolve.n の線形代数 | pysad の式評価 | 翻訳済み・照合済み (`DiagonalMatrix` は複素数の列を受ける。行列の演算結果は要素ごとに虚部 0 の複素数を実数に戻す (gs と同じ)。LDSolve の終点と 3 通りの呼び出しが gs と一致。TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo も表示桁で一致)。PeriodicR は呼び出し方が不明で比較できず | test_evalfixes.py |
| OpenRead (開けないファイル) | 入力ファイルの読み | pysad の IO | gs と同じ (開けなければメッセージを出して入力行を打ち切り、Read・Close は wrongtype) | test_evalfixes.py |
| tcorr.f、tvcorr.f (DELCOR) | 相関付き誤差の補正 | ffscmd `_delcor` | 翻訳済み・照合済み (チューン 1e-14、軌道 1e-10、乱数列も同じ)。k64 モード: MARK は DETR・DPX に入る (K64 の添字のずれ B40)。CAVI は K64 のとおり HARM・PHI に書き gs と一致、OCT はチューンと QF の DX が gs と 12 桁で一致 (OCT 自体の DX は gs も pysad も 0 で、効果そのものの比較ではない)。corlength≠0 の枝は K64 が cv(n) を設定しない (B62。2 回目以降の corlength≠0 の DELCOR は再現できない)。v2 モード: MARK・CAVI も DX・DY に入れる (B40)、波数 −1 の成分を対称に置く (B62)、corlength≠0 の規格化を rms = delx に直す (B99。作者の文書による裏づけは無い推定)。v2 は統計・解析式で確認し、gs との比較は無い | docs/04_ffs_small_commands.md、test_delcor.py |
| tcftr.f (FFT) | 実 FFT | csr.py (CSRTrack 用)、ffscmd・maintrack | 翻訳済み (CSR 側は gs と照合) | docs/15 |
| trcoda.f (TRCOD) | 追跡による COD | measure.py (v2 モード)。k64 モードは ffscmd が受理のみ | v2 モードで翻訳 (gs で SIGSEGV: B42。gs との比較は無い。上の MEASURE の行) | docs/03 |
| terror.f (DELX/DK/DTHETA…) | 誤差コマンド | ffscmd `_error_cmd` | 翻訳済み・照合済み | L3C/3D、L7-2 |
| MON/MCOD/TCOD コマンド | COD 補正 | — | 対象外 (K64 の tffsa.f にこれらのコマンドは無いことをソースで確認 (2026-11-16)。FFS では式として評価されて何もしない) | — |
| tsvdm.f 等 | SVD | match6 の tfsolv が numpy で代替 | 翻訳済み・未照合 (numpy 流) | マッチング K1 ≤1e-10 |

## 6. FFS・マッチング

FFS のコマンドループとマッチングは翻訳済みで、FITP>1 と fam 点だけ参照値が無い。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| tffsa.f (CALC:544、CHRO:609、RADINT:617、EMIT:988、GO:1207 の各分岐) | コマンドループ | ffscmd.py | 翻訳済み・照合済み | L13, L14 |
| tffscalc.f (tffscalc:12)、tffsmatch.f (tffsmatch:15。tfsolv、tgrad、fmincube、tweigh) | 残差・GO (Newton/DESCEND, fuzz, fmincube) | match6.py (`_tffscalc`、`Matcher.go`、`_fmincube`、`_tdfun1`) | 翻訳済み・照合済み | 反復数・残差・K1 ≤1e-10 (L13A/B)。FITP>1 は参照値なし (L15A) |
| tffsqu / qdcell | 解析微分 | match6d.py | §2 参照 | L22A |
| tgetfv.f, tfgetv.f, tflag.f, tfsetv.f | FIT 条件・要素値・フラグ | ffscmd.py | 翻訳済み・照合済み | L13, L14 |
| tftwiss.f、tfchro.f (tfchro)、tfoptics.f (117 行。qcell1 の薄いラッパー) | Twiss[]/LINE[]/Element[]、色収差、CalculateOptics | ffscmd.py, lattice.py | Twiss[]/LINE[]: 照合済み (L14A、PolynomialFit 1e-7)。CHROMATICITY (tfchro)、CalculateOptics (tfoptics; 4 次元の cell/INS・区間・結合・分散、idim=3 は全周): 翻訳済み。CHROMATICITY は QUAD / SEXT / MULT ごとの薄肉近似の表で、色収差の値はチューンの DP0 スキャン + `PolynomialFit` (Fit.n) でも求められる | — |
| tfgeo.f (tfgeo:1、tfgeo1、tfchi)、tdgeo.f、tgrot.f | 幾何 (GX…、ORG) | geo.py (LINE["GX"] / fit GX.. / DISP G 経由)、match6d.py `tdgeo_term`。`tfchi`, `tfrotgeo`, `tgrot`, `tforbitgeo` | 翻訳済み・照合済み | 4.8e-13 (L19A/B/C)。tdgeo は DRIFT の L・BEND の ANGLE・L (幾何の条件だけの行) を解析微分に翻訳 (有限差分と 1e-5 以内、test_tdgeo.py)。四極などの DX・DY、Twiss の条件が同じ行にある列は中心差分のまま (半セル ARCAH の例で収束位置が gs と相対 5.5e-7 ずれる) |
| tfdapert.f (tfdapert:1)、trackd.f (trackd:37) | DynamicAperture[] | dapert.py (独自実装。スロット順・DAPWIDTH 窓・62 進地図 `rad62`。`ffscmd.FFS.fn_dasurvey` から呼ぶ) | 翻訳済み・照合済み | L16A |
| tffamsetup.f | FITP の有限振幅ファミリー (MatchingAmplitude/NFAMP) | match6.py `fam_setup` ほか | 翻訳済み・未照合 (SAD が FITP>1 の GO で再起動するため自己検証のみ) | test_fam_points.py |
| tffswake.f | ウェイク (`tffssetupwake` = 追跡用の表の取得、`tffswake` / `tffswakekick` = 光学) | ffscmd.py `_wake_table` `ini_optics`、wakeopt.py | 追跡・光学とも翻訳済み・照合済み (InitialOrbits = tffsinitialcond、tfshow の軌道列も) | test_wake.py、test_wake_optics.py |
| tmulti.f tmultiacc | RF のある MULT (TRPT か RFSW) | track6n_ext.py `mult_full` | 翻訳済み・照合済み (ANGLE=0、F1・FRINGE なし、AUTOPHI なし) | test_wake.py |
| tfstk, tfeval, tfetok 等 (SADScript) | インタプリタ | sadscript.py、sadlist.py (リスト系組み込み) | 翻訳済み・照合済み | L14A |

## 7. ビームライン定義・座標

MAIN 部とジオメトリは翻訳済みで、PHSROT だけモードで状態が違う。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| toplvl / doelem / doline / expln | MAIN 部 (LINE, 繰り返し, 再定義) | parser.py | 翻訳済み・照合済み | L12 (3563 要素) |
| talign.f, nalign.f | ずれ・回転の設定 | ffscmd, lattice | 翻訳済み・照合済み | L3B/3C |
| tfgeo.f | フロアー座標 | geo.py | §6 参照 | 4.8e-13 |
| tfloor.f | 数学関数 (cfunc) | 各所に個別実装 (xsincos 等) | 翻訳済み・未照合 (個別) | — |
| tsgeo.f | ソレノイドの幾何 | sgeo.py | 翻訳済み・照合済み (ずれ・傾きのある SOL の幾何)。k64 モード: 領域内の BEND・MULT・CAVI、領域内の L=0 / K1=0 の四極は未対応 (NotImplementedError)、tdgeo の領域幾何は未翻訳。v2 モード: 領域内の BEND・MULT・CAVI・K1=0 の四極も軌道を追い、領域内の変数の幾何 fit は中心差分 (設計軌道の粒子が出口で 0 に戻る 1e-12。gs との比較は無い) | test_sol_geo.py |
| phsrot.f | PHSROT (位相空間の回転) | phsrot.py | k64 モード: 未対応 (K64 で動作しない: 行列がすべて 0。追跡は NotImplementedError、光学は恒等。docs/30 B38)。v2 モード: 意図された写像 M = N·Rot·(1−D)·N⁻¹ + 拡散を実装 (追跡・CALC・Emittance[]・EMIT。固有値・シンプレクティック性・平衡エミッタンスで確認。gs との比較は無い。K64 の拡散の乱数の誤り B79 も直している)。スピンは動かさない | test_v2_C_phsrot.py |
| tpara, tscale 等 | 空・デッド | — | 対象外 | — |

## 8. 寿命・その他

寿命・ウェイク・空間電荷関数・特殊関数は翻訳済みで、MAIN の TRACK の Fourier 解析 (wjfit 系) だけ未翻訳である。

| Fortran (ルーチン:行) | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| Touschek (tintra 系、emit.n) | TouschekTable/Lifetime | touschek.py | §3 参照 | ≤7e-12 |
| tflifetrack.f | LifeTrack | — (SAD 側も空スタブ) | 対象外 | — |
| tlum.f, tluma.f | ルミノシティ (MEASURE の統計) | measure.py | 翻訳済み・未照合 (MEASURE が gs で動かないので解析値との比較だけ。ffscmd からは呼ばれず、test_measure.py のみ。BEAMBEAM の LuminosityMonitor は beambeam.py) | test_measure.py |
| txwake.f | ウェイクのキック (TrackParticles) | wake.py | 翻訳済み・照合済み (普通の要素、CAVI、SOL の領域の内側) | test_wake.py |
| tfwake.f、twinit.f 他 | 旧 WAKE コマンド | — (同じ経路で動く) | 移植不要: K64 の tffsa.f:972 は tfwake.f を呼ばず WAKECOMMAND[] (配布に無い) を評価する。gs はエラーを出して行の残りを捨て、pysad も同じ経路 (エラー文の表示だけ未再現) | docs/13 |
| wjfit、wfres、wfbin、wiord、wjdraw | MAIN の TRACK の Fourier 解析・共鳴の描画 | — | 未翻訳 (FFS から使えない。maintrack.py の `FOURIE ON` も未対応) | — |
| tspac.f (GaussianCoulomb) / tintrb.f (twspu、twspfu) | 空間電荷関数 | gausscoul.py | GaussianCoulomb・GaussianCoulombU・GaussianCoulombFitted を翻訳し gs と 1e-9 以内で一致。`tspac.f` の `tspac` は K64 のどこからも呼ばれない死んだコード。追跡の空間電荷 (SPAC・PSPAC・WSPAC) は v2 モードで `spch.py` に実装 (上の行) | test_gausscoul.py、docs/06 |
| tlinit.f (tlfield、tlcalc、tlspect) | RadiationField・RadiationSpectrum | radfield.py | 翻訳済み。RadiationField は gs が SIGSEGV (B49) なので等速運動の電荷の速度場と相対 6e-3 の確認のみ。RadiationSpectrum は gs と一致 (NPARA の fork 並列は逐次) | test_radfield.py、docs/18_radlight.md |
| gamma.f、hg.f、tfunc*.f (特殊関数)、tffindroot.f、D.n | 特殊関数・Fit・FindRoot・NIntegrate | specfn.py、fitroot.py、miscfn.py | 458 式が gs と一致 (K64 が誤値を返す B43〜B47、B53、B54 は真値)。Fit・FindRoot は数値微分で 5e-10〜1e-9。NIntegrate は複素数の被積分関数を扱い (gs と虚部 1e-15)、`LogGamma1` もある。複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)、複素引数の Bessel・超幾何・Airy (実数のパラメータ) と PolyLog の非整数次の大きい \|z\| (相対 1e-11、\|z\| ≤ 1e12) は両モードで実装。k64 と値が違う関数は InverseErf (B47) だけで、v2 は真値。複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta は両モードで実装 (`specfn_cx.py`。gs が誤る点でも真値: B118〜B120) | docs/06 |
| ActLie, ActGRA, tdlat, tmap, polpar 他 | 空スタブ・描画・デッド | — | 対象外 | — |

### 8.1 数学 🟡B

数学ルーチンは numpy・scipy で置き換えている。

| ファイル | ルーチン | 機能 | pysad |
|---|---|---|---|
| tsolvg.f | tsolvg, tsolva, tsvd | 最小二乗・SVD | `numpy.linalg` / `emit6d._tsolva` (SVD 細部の対応は未確認) |
| teigen.f | teigenc ほか | 固有値 | numpy (qmdiag も固有値分解) |
| matrix.f / tmultr.f / tinitr.f | — | 6×6 行列 | numpy (`match6d.tmultr45`) |
| spline.f | spline, splint | スプライン | `touschek.SadSpline`, `_spline_ddy` |
| gaussn.f / drndsr.f / random_plugin_sad.c | — | 乱数 | `sadrandom.py` (LCG、Box-Muller、drndsr、tdusrn) / seed 一致 |
| tffindroot.f | — | FindRoot / Fit | `fitroot.py` (tffit1・tfnewton の移植。FindRoot・Fit・FitGaussian・NIntegrate。gs と 5e-10〜1e-14)、PolynomialFit は `sadscript.py` で Fit.n を評価 (docs/06) |
| gamma.f / hg.f / D.n | — | 特殊関数 | `specfn.py` (scipy.special。458 式が gs と一致。K64 が誤値を返す式 (B43〜B47) は真値を返す。docs/06) |

### 8.2 インタープリタ 🔵C / グラフィックス 🟢D / C ソース・sim・inc・Packages・extensions・contrib

要約:

- 🔵C: tfstk/tfetok/tfeval/tfestk/tfeeval/tfefun/tfefun1/tfinitn/tffs/tffsa/tffscalc/toplvl/terror/initb1。pysad は `sadscript.py` (字句・Pratt 構文・評価、Module/Do/Table/パターン/ファイル I/O、autofg 書式)、`sadlist.py` (MapThread, Cases, MatchQ 等のリスト系組み込み)、`parser.py` (MAIN 部)、`ffscmd.py` (tffsa の語の解釈順、tflag/terror/tftype) で必要分のみ再実装。
- 🟢D: tfTkInter_.c、tfXlib_.c、sad_tcltk.c、sad_xlib.c、geodrw.f、Packages の Plots 系。pysad の描画は matplotlib (`plot_optics.py`) が基本だが、**例外あり**: GUI (Tkinter.n・KBFrame.n) と描画 (ListPlot・Plot・OpticsPlot など) は SAD 自身の .n を pysad のインタプリタで解釈し、tkinter のキャンバスに描く (`tclbridge.py`・`tkcanvas.py`・`pkgload.py`。docs/20・21)。tfTkInter_.c のキャンバス関数は `tkcanvas.py` が tkinter の Tcl に置き換える。matplotlib の描画 (`plot_optics.py`) は別に残る。
- C: calc.c は SAD1 の数式用 (byacc 生成)、PySad.c は Python 2 API で非ビルド、tspch_.c/spch.c は空間電荷、tgauss_.c/random_driver.c/random_plugin_sad.c は乱数、sim/ は OS 抽象化 (unix_memory 系、sad_api、dynl 系ほか)。
- inc/: 実コードを含むのは TENT/TEXIT/common/MACexpn/TFMACRO/TMACRO/TFSF/coroper 程度。要素コードの実定義は tfstk.f の module maccode (icDRFT=1, icBEND=2, icQUAD=4, icSEXT=6, icOCTU=8, icDECA=10, icDODECA=12, icUND=18, icWIG=19, icSOL=20, icST=21, icMULT=22, icCAVI=31, icTCAV=32, icMAP=33, icINS=34, icCOORD=35, icBEAM=36, icPROT=37, icSPCH=38, icMARK=41, icMONI=42, icAPRT=43, icLINE=100)。pysad の `ffscmd._TYPECODE` が SOL=20, MULT=22 を使う。
- Packages: init.n (NPARA=4 など既定値)、Optics.n (CALC6D, NormalCoordinates)、emit.n (Emittance[] の規則名、TouschekLifetime)、dapert.n、Fit.n、correction.n、optimize.n、extmap.n。さらに `pkgload.py` が correction.n・bump.n・corutility.n・optimize.n・D.n・Plots.n・ListPlot.n・CSR.n などを SADScript のまま読む (`PYSAD_SAD_PACKAGES` が要る。docs/06)。
- extensions: RandomMT (MT19937/SFMT、既定の乱数は "SAD" 合同法)、sad+ (微分代数の C++)、Standard (Tkinter/Xlib)。pysad は使わない。
- contrib: gdtoa、BBBrem_dir (BBBrem.C)、libtai、tcl/tk ほか。不要。

## 未対応・未確認の一覧

未翻訳・一部未対応・未照合の項目、優先度、未解決の差の一覧は docs/31 にある (未翻訳、一部未対応、モードで状態が違うもの、未照合、優先度、見送りの理由、未解決の差)。

## 9. 機能別の状況 (SAD の機能ごと)

目的: SAD-master (`src/*.f`、`Packages/*.n`、SADHelp) の機能ごとに、pysad での対応状況 (済み・見送りと理由) を並べる。照合は `bin/gs` (SAD の WSL ビルド) の出力との比較で、数値は各文書にある範囲 (小さなリング・試験の条件) のもの。実用のリングや強い非線形での確認は別である。
結論: FFS のコマンド、SADScript の関数、要素・物理モジュールのほぼすべてが対応済みで、未対応は FFSFlags の内容、MAP の GEO・CompiledMap、CanvasSymbol などに限られる (詳細は docs/31)。

### 9.1 FFS のコマンド (`tffsa.f` ほか)

| 機能 | SAD 側 | pysad の状況 |
|---|---|---|
| STOP / USE / VISIT / DISP / TYPE / SAVE / CALC / CAL / FIX / EXPAND / REJ / RESET / EXIT / ABORT / BYE | `tffsa.f` | 済み (`USE BeamLine[..]` も可。反転要素名 `-NAME` の差は docs/29) |
| FIT / FREE / GO / FITP / VariableRange / VariableWeight | `tffsmatch.f`、`qdcell.f` の tffsqu | 済み (解析微分 match6d.py が既定。多解問題の 2 段目以降は SAD と別解のことがある。FITP>1 は SAD 側が再起動を繰り返して落ちるので参照値なし) |
| DP / DP0 / GCUT / COUP / COUPLE | `tffsa.f` | 済み |
| DELX / DELY / DTHETA / DK / DDK / DL | `terror.f` | 済み (直訳、照合済み) |
| DUMP | `tffsa.f` | 済み (DK・DDK 列の書式 (Fortran の 1p g12.5) を `ffscmd._gwd` で合わせている。`test_dump_format_matches_gs`) |
| EMIT (コマンド/フラグ) | `temit.f` | 済み |
| 要素値設定 (`SDA K2 …`、`MAX`/`MIN`)・FITFUNCTION・入口 Twiss (`AXI` `BXI` … `EPXI`)・fit 条件 (`NX` `NY` `BX` …) | `tffsa.f` | 済み |
| UNTIL / REPEAT | `tffsa.f:255-290` (入力位置 `ipoint` を保存して巻き戻す) | 済み (最上位の FFS 文列) |
| IF / ELSEIF / ELSE / ENDIF | `tfif.f` | 済み (同上) |
| CHROMATICITY (CHRO) | `tfchro.f` (QUAD / SEXT / MULT ごとの薄肉近似) | 済み |
| MAIN の TRACK (`TRACK USE= NPART= TURNS= NX= NY= NZ=`) | `track.f`、`tracka.f`、`ttinit.f` | 済み ([`pysad/maintrack.py`](../pysad/maintrack.py)。K64 の MAIN の TRACK は SIGSEGV で動かず、SADScript で同じ手順を書いた版と一致)。PHSPLOTS の描画は無視、FOURIE ON は未対応 (docs/history/work_main_track.md) |
| PRINT / REF / ORG / SHOW / STATUS / MATRIX / BEAMSIZE / VARY / RECOVER / REVERSE | `tffsa.f` の各分岐、`tfmat.f` (MATRIX) | 済み (docs/04) |
| RENUMBER (要素番号の振り直し) | `tfinit.f` の tffsrenumber、tfinimult | 済み (`ffscmd._renumber_cmd`、`LINE["MULT"]` も。INDEX 属性が 2 以上なら番号をその値から振る。`LINE["INDEX",…]` も。gs と一致。両モード) |
| GEO・DRAW | `tfgeo.f` | 警告して無視 (描画・表示系) |
| RADINT | `INTGRL.f` (693 行) | 済み (`radint.py` を FFS RADINT に接続。6 通りと小リング 3 通りが SAD と 7 桁一致) |
| APERTURE / ZAP | `tfaprt.f`、`tfzap.f` | tapert (追跡側の損失判定) は翻訳済み。コマンドは未対応 (表示系) |
| DELCOR | `tcorr.f`、`tvcorr.f`、`tffsa.f:982` | 済み (k64 モードは gs と一致: CAVI は HARM・PHI に書く K64 のとおり、OCT も gs と一致。v2 モードは MARK・CAVI も DX・DY に書き、corlength≠0 の規格化は正しい式。docs/04) |
| MON / MCOD / TCOD | | 受理のみ (K64 ビルドでは SAD も何もしない) |
| MEASURE / PUTMEASURE / GETMEASURE / LOGMEASURE / ORBITJITTER / TRCOD (FIXCOD / MVCOD / STNCOD / ALIGN) | `tffsa.f:639-985`、`trcoda.f`、`trackb.f`、`tracke.f`、`ttinit.f`、`ttstat.f` | **v2 モードで実装** (K64 のソースの読み。gs は MEASURE が無出力・TRCOD が SIGSEGV (B42) なので gs との比較は無い)。k64 モードは受理のみ (警告して無視)。docs/03 |
| LTRACK / PKSEED / PUSHSEED / POPSEED / EXCGSEED / DELWAVE | `tffsa.f:639-985` | 未定義 (未確認) |

### 9.2 ON/OFF のフラグ

| フラグ | pysad の状況 |
|---|---|
| LOG / COD / CELL / RADCOD / RFSW / FLUC / RAD / RING / INS / TRPT / CODPLOT | 済み (INS の初期軌道オフセットも一致。CODPLOT は表示なし) |
| DAPERT | 済み (dapert.py、DA の地図一致。カオス的な領域の最終 seed は不一致) |
| INTRA | 済み (EMIT が IBS 付き、TouschekLifetime の表。`temit.f` の INTRA 経路) |
| CALC6D | 済み (CELL と INS、1.7e-10。`qcell.f` の qcell6d) |
| CTIME / ECHO / EMIOUT / GAUSS / CONV | フラグとして受理。効果はソースで確認 (2026-11-16): GAUSS は MAIN の TRACK・MEASURE の初期分布、CTIME は CPU 時間の表示、CONV は収束の判定・表示、EMIOUT は EMIT 報告の追加行 (pysad は spin の行だけ対応) |
| BIPOL / UNIPOL / RELW / ABSW | BIPOL は FIT だけ (tffsmatch.f)。UNIPOL・RELW は FIT の結果を変え、pysad も match6 につないで gs と全桁一致 (test_item4.py) |
| WSPAC / SPAC / SELFCOD / PSPAC | k64 モード: 受理のみ (空間電荷は追跡に入らない)。**v2 モード: `TrackParticles` に入る** (SPAC = DRIFT・MULT の RADIUS、PSPAC = SPCH 要素、WSPAC = ビーム行列のガウス分布の力、SELFCOD)。光学・Emittance への空間電荷は実装済み (Emittance[]・EMIT は両モード、CALC・CALC6D は v2 だけ。docs/09)。docs/10 §9.7。tspac.f の bbkick1 と tintrb.f の twspu は `GaussianCoulomb(Fitted)` として翻訳済み (docs/06) |

### 9.3 SADScript の関数

「組」は 組み込み (B、`src/*.f` の `itfunaloc`) / パッケージ (P、`Packages/*.n`) / ラッパー (W: パッケージ経由で組み込みを呼ぶ)。

| 関数 | 組 | SAD 側 | pysad の状況 |
|---|---|---|---|
| Twiss | B | `tfetwiss` ほか | 済み。キー NX NY BX BY EX AX AY DX DY EY EPX DPX S `*`、SIGX/SIGY/SIGPX/SIGPY/SIGZ/SIGDP/SIGE/SX/SY (emit_extra.py、2.1e-10) |
| LINE | B | | 済み。キー NAME K1 LENGTH K2 POSITION K0 L ANGLE S ELEMENT DIR DX DY TYPE TYPENAME DK VOLT ROTATE E1 E2 GX GY GZ CHI1-3 OFFSET。RHO/RHOX/RHOY/KX/SK0 は要確認 |
| FFS | B | | 済み |
| Emittance | W (emit.n) | `tfemit.f` (226 行)、`temit.f` | 済み。オプション OneTurnInformation・Matrix・OneTurnTransferMatrix・SaveEMIT・Output は emit_extra.py (2.1e-10)。偏極・SpinTune は 0 固定 |
| DynamicApertureSurvey | P (dapert.n) | | 済み (遅い。既定の Output / ExpandElementValues は ffscmd.py:880 で再現) |
| TrackParticles | B | `tftrack.f` | 済み |
| CalculateOptics | B | `tfoptics.f` (117 行。初期条件 28 個を行 lb に書き込み qcell1 を実行) | 済み (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) |
| Monitor / Steer | P (corutility.n) | `Master[]` (pgmast.f) と `LINE["type"/"element",…]` を呼ぶだけ | 済み (SAD の corutility.n を pysad のインタプリタで動かす。docs/06。FFS の MON とは別物) |
| ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine / ExpandBeamLine | B (tfbeamline.f:362-510) / P (beamline.n) | 逆向き要素は `-1*名前` | 済み (`USE BeamLine[..]` / `VISIT` / `BYE` も) |
| SetElement / Element | W (element.n → SetElement$ / Element$) | `tfbeamline.f:10-340` | 済み。`Element["KEYWORDS"/"DEFAULT"/"VALUE"/"POSITION"]`、型指定・新規定義とも。PROFILE, Saved->True, DK の連動を除く |
| FLAG | B (pgflag.f) | | 済み |
| GetMAIN | B (tffs.f tfmain) | | 済み (MAIN 形式ラティスの読込み) |
| FFSFlags | B (tfetok.f の tfflags) | | 登録名あり。内容は未対応 |
| CorrectOrbit / MakeBump | P (correction.n、bump.n) | `pgrmat.f` (応答行列)、`pgsolvcond.f` (条件付き SVD) | 済み (SAD の .n を動かす。ResponseMatrix・LinearSolveConditioned・Master は `corr.py` に直訳。補正量は gs と相対 1e-7。docs/06) |
| OptimizeOptics / DownhillSimplex | P (optimize.n) | | 済み (K1・チューンが gs と 1e-10。DownhillSimplex は `miscfn.py` で gs と 1e-13。docs/06) |
| TransferMatrix / BeamMatrix / NormalCoordinates / SymplecticInverse | P (Optics.n、emit.n) | | 済み (optics_fn.py、SAD と 1e-9) |
| PolynomialFit / Fit / FitEmit / FitGaussian / Spline / NIntegrate / Integrate / FindRoot / D | P (Fit.n、D.n) | | 済み (`fitroot.py`。gs と 5e-10〜1e-14。ベクトル独立変数は対応 (gs が誤値: B52)、複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)。複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta は両モードで実装 (`specfn_cx.py`。gs が誤る点は B118〜B120。docs/06)。NIntegrate は複素数の被積分関数も可。D・Integrate は D.n を動かす。docs/06) |
| TouschekLifetime / TouschekTable | P (emit.n:129-) | `temit.f` の intra | 済み (15 ケース ≤7e-12) |
| SynchroBetaEmittance | B (`tfemit.f:144` の tfemits) | `temits.f` (1325 行) | 済み (`synchrobeta.py`、gs と相対 1e-5。24 セルの平面リング + CAVI。24 セル FODO + skew 四極の結合のあるリングは gs も落ちず 7 桁一致。**大きなリング A では独立な値と 40 倍違い使えない (B128。v2 は警告)**。不安定なチューンは未確認) |
| SeedRandom / GaussRandom / Random | B | | 済み (SAD と同じ乱数列) |
| LDSolve ほか (ExpMatrix、rational.n、extmap.n、Misc.n) | P | | SAD 自身の `.n` を解釈する方式で動く (docs/06)。LDSolve は複素数・行列演算の虚部の扱いも含めて gs と一致。rational.n (Rational の和・積) と Couple.n の ElementValues (VISIT/BYE の退避・復元を除く) も gs と一致 (docs/06)。他は個別の照合なし |
| Polarization / SpinTune / PolarizationTime | | `temit.f:491,1188` の calpol、`polpar.f`、`photons.f`・`tturn.f` の spin 部 | 済み (Emittance[]・EMIT のスピン (MULT・空洞・SOL まで) と、POL かつ RAD ON の TrackParticles のスピン追跡。gs と 1e-11〜4e-12)。`Polarization` は K64 が NaN (B55)。k64 モードは NaN のまま、v2 モードは有限値 (線形のスピン軌道結合から。平面リング 92.376%。gs との照合は無く、FLUC つきのスピン追跡との比較のみ。docs/16) |
| PHOTONS (`PhotonList`) | | `tftrack.f:144,315`、`photons.f` の kphlist | 一部対応 (BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND の光子に対応 (k64 は gs と全列一致)。SOL 領域内・COORD・TCAVI・MAP のあるラティスの光子は未対応。RADLIGHT の軌道と RadiationField・RadiationSpectrum は docs/18) |
| LMA | | | SAD に組み込みはない |

### 9.4 SADScript 言語系の組み込み

[`pysad/sadlist.py`](../pysad/sadlist.py) に MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection, Union, Position, Count, MemberQ, FreeQ, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map, Fold, Nest, FixedPoint, Through ほかを実装済み (SAD 側は `tfinitn.f` の `itfunaloc` 登録と tfefun*.f / tfeeval.f のディスパッチ)。Total/Cross/Norm/Riffle/StringSplit は SAD にも無い (未評価のまま返る。gs で確認)。Sleep は no-op、Fourier / InverseFourier / BesselJ / BesselK は FFT / Bessel で実装。

### 9.5 パッケージ経由のものの注意

Emittance / Element / SetElement / Monitor / Steer / DynamicApertureSurvey / TouschekLifetime / PolynomialFit は SAD では SADScript パッケージ (`Packages/*.n`) の関数で、Fortran の組み込みではない。Python で個別に書き直しているものは、SAD のパッケージの既定値・オプションとの差が出やすい。`pkgload.py` が SAD の `Packages/*.n` を最初に使うときに読む (環境変数 `PYSAD_SAD_PACKAGES` が要る。docs/06)。

### 9.6 要素・物理モジュール

| 機能 | SAD 側 | pysad の状況 |
|---|---|---|
| BEND の K0 キック (steering) | `tsteer.f`、`tbende.f`、`tbend.f` | 済み (CALC の tsteee、追跡の tsteer) |
| 幾何 (LINE["GX"/"GY"/"GZ"/"CHI1".."CHI3"]、fit の GX) | `tfgeo.f` (411 行)、`tdgeo.f` | 済み (geo.py、4.8e-13。fit の幾何の微分は、DRIFT の L と BEND の ANGLE・L だけ解析微分 (tdgeo)、ほかは中心差分。SOL 領域内の四極の DX・DY の解析微分は `match6d.sol_dcod`。CHI3 の微分は K64 でも領域の終端のリセットと合わず fit は収束しない: B124) |
| MULT ANGLE≠0 | `tmulta.f` (433 行) | 済み (docs/19) |
| SOL | `tsole.f`、`tsolque.f`、`tsgeo.f` | Twiss・Emittance[]・EMIT・追跡 (k64 モードは領域内が DRIFT・QUAD のみ)、幾何は gs と一致 (docs/17)。v2 モードは領域内の BEND・MULT・CAVI・APERT、ずれて傾いたソレノイドと F1、領域内の空洞の TRPT の加速・wake・AUTOPHI も (ローレンツ力の積分との一致で gs との比較ではない)。K64 は SOL (BOUND≠0) の長さを周長に入れない (B63) ので同じにしている |
| TCAVI | `ttcav.f` (72 行) | 済み (docs/19)。L≠0 の Emittance[] が gs で壊れる: B36 |
| COORD | `tcoord.f` (82 行) | 済み (docs/19) |
| MAP | `temap.f` (625 行) | 済み (docs/19)。`ExternalMap["GEO"]` は未対応 (直進、長さ 0) |
| UND | `undulator.f` (165 行) | 済み (docs/19)。磁場は K64 に無い: B37。v2 モードは `FBX`・`FBY` の磁場つき (解析式で確認) |
| WAKE | `tffswake.f`、`txwake.f` | 追跡と光学 (WAKEOPT) が gs と一致 (docs/13)。旧式の WAKE コマンド (tfwake.f) は K64 が WAKECOMMAND.n (配布に無い) を呼ぶので k64 モードは同じ経路。v2 モードは `WAKE … ENDWAKE` の表を読む (B66)。wake の重み (B84) は v2 では正しい式 |
| BEAMBEAM | `bbstrhl.f` | 強弱ビームが gs と一致 (docs/14) |
| CSR | `csroy.f` ほか、CSR.n | CSRTrack・CSRHaissin・CSRSetupOY・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRInit・CSRConvert、MakeZL の 3 周波数 (docs/15)。ScanTrack・OYStability は未照合 |
| RADLIGHT・放射場 | | TrackParticles が軌道 {t,x,y,z} を返す (k64 モードは gs と差 0。v2 モードは 3 次元の軌道で BEND の中も記録: docs/18)。RadiationField・RadiationSpectrum を実装 (RadiationField は gs が落ちる: B49) (docs/18) |
| GaussianCoulomb(Fitted)・BBBrem | `tspac.f` の bbkick1、`tintrb.f` の twspu、`tbbbrem.f` | 直訳 (docs/06) |
| ルミノシティ | `tlum.f`、`tluma` | `measure.py`。MEASURE が gs で動かないので解析値との比較だけ |
| 空間電荷の追跡側 (SPCH 要素) | `spkick.f`、`tspch_.c`・`spch.c`・`psn.c`、`tintrb.f` の twspac、`tturn.f:521` | **v2 モードで実装** (`spch.py`。K64 の SPAC の誤りは B89〜B98。k64 モードは入れない)。`tspac.f` の `tspac` は K64 のどこからも呼ばれない死んだコード。PROF は K64 に無い。docs/10 §9.7 |

### 9.7 GUI・描画

Tkinter.n・Plots.n・OpticsPlot.n などを SAD の .n のまま tkinter で動かす (docs/20・21)。BeamPlot は未比較。CanvasSymbol・Canvas3D* は未着手。

## 10. `<SAD>\src` の *.f のうち pysad の中に名前が出てこないもの

`<SAD>\src\*.f` 474 本を `grep -w -i` で `pysad/*.py` と突き合わせた (2026-11-16)。名前が出ないだけで、別名のルーチンに直訳されている場合がある
(例: `tins`・`tinse` は `lattice.py` の INS の扱い、`talign`・`nalign` は ffscmd / lattice のずれ・回転)。
物理計算で名前が出ないもの: `qins tins tinse` (INS、§2。実装は別名で gs と一致)、`qmult qddrif` (qtwiss・qdcell 内の別実装に吸収: docs/29 D5)、`tfmap tmap` (MAP は別名)、
`tflifetrack twinit wfbin wfres wiord wjdraw` (§1・§8・docs/31 の未翻訳の表。`trackb` は v2 の `measure.py` の MEASURE の実装で名前が出るので、この一覧の再現 (2026-11-17) では出ない。名前の再数えはしていない: 要確認)、`tbrad tquads tquase tsdrad tserad tstrad tqrad trad trade` (Obsolete)、
`tconv` (追跡の変換は track6 内)、`tmuld tmultr` (行列の積は numpy)、`temitf` (EMIT 報告は ffscmd。ファイル名としては出ない)。
それ以外は、旧 MAIN のユーティリティ (`pmbump pwrite pundo pkill` など)、SADScript のインタプリタ内部 (`tfeval1 tfstk tfpart tfsort tfreplace itf*` など)、描画 (`geodrw phdrw pmbdrw tdplt twsdrw`)、ファイル・バッファ (`filaux filbuf flmgr getbuf`) で、
pysad の方針 (SADScript は Python のインタプリタ、ファイル・描画は Python のライブラリ。`36_sadhelp_purpose_policy.md`) により翻訳しない。
一覧の再現: `for f in <SAD>/src/*.f; do b=$(basename $f .f); grep -q -w -i -- "$b" pysad/*.py || echo $b; done`。
