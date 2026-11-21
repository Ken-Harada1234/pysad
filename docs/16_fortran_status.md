# 16. SAD Fortran ファイル別 pysad 翻訳状況

最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。内容は `history/work_unimpl.md` と `pysad/*.py` の docstring、v2 の作業の記録 `history/v2_changes_*.md` で裏を取った。§9 は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) と同じ状態)

**互換モード**: 表の「状態」は、特記がなければ **k64 モード** (K64 (gs) と同じ動作。gs との照合はこのモード) のもの。v2 モード (v2.0.0 から既定) で変わった動作は「v2 モード:」と書いた (一覧は [47_v2_changes.md](47_v2_changes.md))。v2 の値は gs ではなく解析式・独立な数値計算で確かめたもので、gs と比較できないものがある (MEASURE・TRCOD、空間電荷の SPAC、領域内の SOL の要素など)。

SAD は K64 版、ソースは `<SAD>\src`。
「デッドコード・ビルド対象外」の判定は mk/sad.obj.mk の有効リストに無い、または call 0 件のもの。行数の多いデッドコードは節ごとに 1 行にまとめる。
`<SAD>\src` の *.f 474 本のうち、ファイル名・ルーチン名が `pysad/*.py` に 1 度も出てこないものの一覧は §10 に付けた (名前の grep だけの確認で、名前を変えて直訳した場合は漏れる)。
根拠: `01_overview.md`、`10_benchmarks.md`、`history/00_RESUME.md`、`work_*.md`、`pysad/*.py` の実コード。

## 0. 凡例と現在の pysad モジュール

| 状態 | 意味 |
|---|---|
| 翻訳済み・照合済み (桁) | pysad に直訳があり、SAD の出力と比較済み。桁は最大誤差 (`10_benchmarks.md` の値) |
| 翻訳済み・未照合 | コードはあるが SAD との数値比較が未実施、または再現が不完全 (括弧内に理由) |
| 未翻訳 | pysad にコードがない (NotImplementedError で止まるものを含む) |
| 対象外 | デッドコード、非リンク、GUI/グラフィックス、EPICS など |

現行の対応モジュール: parser / sadelem / lattice / elkern / perfkern / perfkern_bend / turnkern / emit6d / teigen / match6 / match6d / tsvd / track6 / track6n / track6n_ext / dapert / sadrandom / sadscript / sadlist / sadclass / ffscmd / ffsextra / calc6d / touschek / geo / sgeo / emit_extra / radint / maintrack / trptacc / optics_fn、
追加分 (2026-11-13 以降の要素・機能別): coord / tcav / und / multa / multae / solmap / solemit / spin / csr / synchrobeta / wake / wakeopt / beambeam / gausscoul / bbbrem / measure / radfield / specfn / fitroot / miscfn / vecfn / fileops / corr / pkgload / tclbridge / tkcanvas、
v2.0.0 で追加: compat (互換モード) / phsrot / radtaper / spch (空間電荷)。`measure` は v2 で MEASURE・TRCOD ほかの実装を足した。
(`pysad/*.py` は 71 本・計 57124 行。`tracking.py` の ttcav 等への言及は旧 API の名残で、実装は `tcav.py`)
`ffs.py`, `matching.py`, `tracking.py`, `optics.py` は非推奨の旧 API (SAD 未照合の独自実装) で、表には数えない。
calc6d / match6d / touschek / geo / track6n_ext / emit_extra はいずれも `ffscmd.py` から呼ばれる (CALC6D、解析微分、INTRA/TouschekLifetime、LINE["GX"..]/GEO/ORG、MULT/SOL の numba 追跡、Emittance[] のオプション)。

試験の現状 (2026-11-17、`PYSAD_SAD_PACKAGES` つき): `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` で **937 通過・4 skip (941 件収集)**。互換モードは k64 が既定 ([`pysad/tests/conftest.py`](../pysad/tests/conftest.py))。v2 の試験 (`test_v2_*.py`) は 29 本・363 件で `compat.set_mode('v2')` を明示する。回帰 (`QUICK=1 run_regression.py`) は **50 ベンチ FAIL 0** (k64 モード)。

照合桁の目安 (10_benchmarks): CALC 光学 ≤1e-10 (誤差・オフモメンタム・結合経由は 1e-7〜3e-7)、Emittance[] ≤1.3e-11 (Level 3C の結合は 7.8e-7)、追跡 ≤5e-11 (Level 7/9/11 の z・trf0 は ≤1.7e-10)、マッチング K1 ≤1e-10、DA は地図・Score 一致。

## 1. 粒子トラッキング

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| tturn.f / tturn1 | 1 周のディスパッチ、tapert | track6.py `_tturn1`, track6n.py `turn` | 翻訳済み・照合済み | ≤5e-11 (L6〜12, 16) |
| tdrift.f | ドリフト (tdrift_free, tdrift_solenoid) | track6 `_tdrift_free`, `_tdrift_solenoid` | 翻訳済み・照合済み | ≤5e-11 (L6, L9) |
| tquad.f (tthin 同居) | 四極・薄い多極、ttfrin | track6 `_tquad`, `_tthin_gen`, `_ttfrin*` | 翻訳済み・照合済み | ≤5e-11 (L6, L8) |
| tbend.f, tbendi.f | BEND、K1≠0 の複合 (bendib) | track6 `_tbendcore`, `_tbendi` | 翻訳済み・照合済み | ≤5e-11、z のみ 1.7e-10 (L7) |
| tbfrie.f の tbedge, tqfrie, tqlfre, tqente | フリンジ | track6 `_tsetfringep`, `_tqlfr`; lattice.py | 翻訳済み・照合済み | ≤5e-11 (L8)、emit/CALC 側も一致 |
| tcav.f | RF 空洞 (tcavfrin 含む) | track6 `_tcav` | 翻訳済み・照合済み | ≤5e-11 (L6, L8-2) |
| tmulti.f | MULT (ANGLE=0, K0..K21, ずれ・回転) | track6 `_tmulti_ens`、numba 版 track6n_ext.py | 翻訳済み・照合済み (track6 経由) | ≤5.9e-11 (L11)。ext は track6 と 4e-17、seed 一致 |
| tsolqu.f / tsol.f | ソレノイド領域 (bz 分岐、tsolqur) | track6 `_tsol_region`, `_tsolqu_bz_*`、track6n_ext.py | 翻訳済み・照合済み (領域内の DRIFT・QUAD・SOL)。v2 モード: 領域内の BEND・MULT・CAVI・APERT、ずれ・回転 (trots)・F1 (tradks)、領域内の空洞の TRPT 加速・wake・AUTOPHI も実装 (ローレンツ力の数値積分との一致で、gs との比較ではない。docs/44 §7a)。v2 は tsolqur の bzp (B08)・tsolqum の dv 項 (B14) も直した | ≤1.5e-10 (L9)。EPS≠1 は制限あり |
| tphoton / photons.f の tradk | 放射キック (平均・FLUC) | track6 `_tradk`, sadrandom.py | 翻訳済み・照合済み | ≤4.9e-11、乱数状態一致 (L10) |
| tapert.f | 損失判定・APERT 要素 (tapert1) | track6n.py `k_apert_gen`・`k_apert`・`k_apert_el`、track6.py `_apert_pass`、dapert.py、`TrackRun.l_track0` | 翻訳済み・照合済み (APERT 要素: 楕円・矩形・DX DY・ROTATE・全部 0・複数周・JDPX JDPY DP が gs と一致、`tests/test_apert.py`。SOL 領域内の APERT は k64 モードでは未対応 (v2 モードは領域内の APERT で落ちる)。TRPT・RAD つきは未確認。JDPX の乱数は粒子 12 個以上で gs と不一致) | DA の地図一致 (L16A)。ただし損失周の乱数消費は不一致 (RAD+FLUC)。K64 は落とす印に `kptbl(i,4)=l_track` を書くので、セッション最初の TrackParticles (l_track=0) の最初の tapert では最初から LOSSAMPL の外の粒子が落ちない。k64 モードはこの癖を再現し gs (1・2 回目) と全桁一致 (B33 の出力座標の差は解消、docs/11 A52)。v2 モード: 最初の tapert で落とす。LOSSMAP (戻りに落ちた要素・周の 2 行) は両モードで gs の例と一致 (docs/08 §9.7) |
| tsteer.f | BEND ANGLE=0, K0≠0, L≠0 | track6.py `_tsteer_ens` (F1/FB1/FB2 フリンジ付きは NotImplementedError) | 翻訳済み・照合済み | bendk1_* (history/work_bend_k0.md, history/work_trpt_acc.md 10 節) |
| tmulta.f (tmulta, tmultae) | MULT ANGLE≠0 | multa.py (追跡)、multae.py (Twiss・Emittance[]) | 翻訳済み・照合済み (docs/28) | 追跡は相対 0、Twiss 1e-13、Emittance[] 1e-11。BZ・RF・F1・SK0 は未対応 |
| ttcav.f, ttcave.f | TCAVI (クラブ空洞) の追跡と 6 次元写像 | pysad/tcav.py | 翻訳済み・照合済み (追跡 0、Emittance[] は L=0 のみ gs と一致。L≠0 は gs が壊れる) | docs/26_tcav.md |
| tcoord.f, tcoorde.f, qcoord.f | COORD 要素 | pysad/coord.py | 翻訳済み・照合済み (追跡 0、Twiss・幾何 1e-11、Emittance[] 5e-10) | docs/25_coord.md |
| temap.f (temap/temape/qemap) | MAP 要素 (ExternalMap) | elements.Map、ffscmd._map_*、emit6d._tmap、track6n.K_MAP | TRACK・EMIT・OPTICS を翻訳し gs と照合済み (GEO・CompiledMap は未対応) | test_map.py |
| undulator.f | UND (追跡) | pysad/und.py | k64 モード: 翻訳済み・照合済み (追跡 0、Twiss・Emittance[] は恒等で一致)。K64 は磁場を設定できない (B37)。v2 モード: `FBX`・`FBY` の磁場つき (追跡・Twiss・Emittance[]・EMIT。dpds の二重掛け B78 も直した)。磁場なしでも光学にドリフトとして入る。track6 のスピン追跡も通る。解析式で確認し gs との比較は無い | docs/27_und.md |
| bbstrhl.f, bbtool.f, matrix.f (包絡行列) | BEAMBEAM | pysad/beambeam.py | 追跡 (track6n) を翻訳。BSTRL と光子 (PHOTONS、BEAMBEAM の分だけ)、LuminosityMonitor も対応。光学は恒等。傾き tilt の K64 の癖 (B34) は k64 モードでは直訳、v2 モードは (cos t, sin t) の回転に直した (残る差は docs/24) | docs/24_beambeam.md、docs/39_beambeam_guide.md |
| photons.f の sprot、tradk のスピン部 | TrackParticles のスピン追跡 (POL かつ RAD ON のときだけ) | track6.py (`sprot`、`bsi`) | 翻訳済み・照合済み (gs と 4e-12 以下: 平面リング、ずれ・回転した四極・偏向・6 極、MULT、K1 のある BEND、SEXT・薄い MULT、SOL 内のドリフト)。k64 モードは損失 (flag) の処理なし (v2 モードは落とす。領域内の BEND・MULT・CAVI・ずれ・F1、UND、PHSROT のスピンも通る。docs/43 §8a)。SOL 内の四極のスピン追跡は gs の食い違い (docs/44_sol_guide.md §7) のため未確認 | test_spin_track.py、docs/36、docs/43_pol_guide.md |
| tlinit.f (tlinit/tlstore/tlresult) | RADLIGHT の軌道の記録 | ffscmd `_RadLight` | k64 モード: 翻訳済み・照合済み (TrackParticles の戻り `{結果, 軌道{t,x,y,z}}`、gs と差 0。K64 は BEND の中を記録せず横方向の座標も入れない)。v2 モード: 3 次元の軌道、BEND は円弧で中も `RADLIGHTDIV` 点記録 (B65。Liénard-Wiechert のピークと 0.25%、gs との比較は無い) | test_radlight.py、docs/45_radlight_guide.md |
| tbbbrem.f | BBBrem1 関数 | bbbrem.py | 翻訳済み・gs と一致 (docs/32) | test_bbbrem.py |
| csrtrack.f, csroy.f, csrtest.f | CSR 関数 (CSR.n が使う 7 関数: CSRTrack CSRHaissin CSRSetupOY CSROYMatrix CSRMatrix CSRConvert CSRInit) | pysad/csr.py | 翻訳済み。CSRTrack・CSRHaissin は gs と 11 桁、CSRSetupOY は J・ω 1e-9、CSRMatrix (mode 22・33・23・32・110・220・330 と mode 1000 = csrhelmat)・CSRConvert (= K64 の tfcsrvec)・CSRInit は gs と相対 1e-10〜1e-13。CSROYMatrix は gs が wrongtype を返す (B60) ので未照合。`Packages/CSR.n` の MakeZL・SkinDepth・PathDifference などは SADScript のまま動かして gs と一致 (MakeZL は約 7 分、相対 1e-10)。ScanTrack・OYStability などの長い流れは未照合。`CSR.n` 全体は描画 (ListPlot) が要るので読み込まない | docs/36、docs/40_csr_guide.md、test_csr.py |
| track.f, tracka.f, ActTra.f, ttinit.f | MAIN の TRACK コマンド・初期分布 | ffscmd `TrackParticles[]` (SADScript 関数経由)、MAIN の TRACK コマンドは maintrack.py | 翻訳済み・照合済み | TrackParticles は L6〜12。MAIN TRACK は SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 (SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV)。PHSPLOTS の描画は無視、FOURIE ON は未対応 |
| trackb.f, tracke.f, trcoda.f, ttinit.f, ttstat.f, tbdecoup.f, nalign.f, tfojit.f, tfgetm.f | MEASURE・TRCOD・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE、MAIN の `TRACK TURNS=-1` | `measure.py` | **v2 モードで翻訳** (K64 のソースの読み。TRCOD の FIXCOD/MVCOD/STNCOD/ALIGN も)。k64 モードは語を受理するのみ (無視)。**gs との比較は無い** (gs で TRCOD は SIGSEGV、MEASURE は出力なし: docs/21 B42)。解析式と pysad の独立な経路で確認。ttinit の読みの誤り B41 は v2 で直した。DRAW の中の粒子追跡 (tdrwdt.f) は K64 でも動かない死んだコードで翻訳しない (B105) | test_v2_C_measure.py、test_v2_C2_measure.py、docs/03 |
| tturn.f・tcod.f・tturne.f の radtaper | RADTAPER (磁石の強さを閉軌道の局所の運動量に合わせる) | `radtaper.py` | 翻訳済み・gs と照合 (両モード。CODCONV=1e-4 の gs と x 6 桁・U0 10 桁・チューン全桁)。K64 の tcod の反復の直訳ではなく同じ不動点を外側の反復で求める。4 次元の CALC のテーパーは K64 も無い | test_v2_C_track.py、test_v2_C2_optics.py |
| spkick.f、tspch_.c・spch.c・psn.c、tintrb.f twspac | 空間電荷 SPAC・PSPAC・WSPAC | `spch.py` | **v2 モードで翻訳** (k64 モードは入れない)。K64 の SPAC は力が 1/γ 倍に弱い (B89) ほか B90〜B98 の誤りがあり、v2 は直した。PSPAC は gs と 3〜4 桁、WSPAC は gs と 1% (Δν)、SPAC は gs の誤りのため解析式 (包絡方程式・Laslett の式) と比較。Emittance[]・EMIT・CALC6D・CALC (RFSW ON) の WSPAC は v2.1.0 で実装 (Emittance[]・EMIT は k64 も gs と全桁。K64 の最後のスライスのレンズが 8 倍強い誤り: B127)。追跡の WSPAC のビーム行列は空間電荷なしのまま | test_v2_C_spch.py、test_v2_C2_spch.py、docs/08 §9.7 |
| trackd.f / tfdapert.f | DynamicAperture 本体 | dapert.py | 翻訳済み・照合済み | 地図 153/153、Score 一致 (L16A、NORAD) |
| tquads, tquase, tbrad, tmulbs, ttdr, bb, hg 等 | Obsolete・非リンク | — | 対象外 | — |

## 2. 線形光学

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| qtwiss.f (qtwiss1, qmat2twiss) | Twiss 伝搬 | lattice.py `twiss_table` | 翻訳済み・照合済み | ≤1e-10 (L1,12,13)、結合・誤差経由 ≤3e-7 (L3, L5)。結合時の qmat2twiss は部分 |
| qtwiss.f qcod / tcod.f | 4D 閉軌道 | lattice.py `_qcod_periodic`, `find_cod` | 翻訳済み・照合済み | 9〜10 桁 (L3B) |
| qcell.f (qcell1) | CELL 周期解・INS | lattice.py `twiss`, `ins_twiss` | 翻訳済み・照合済み | CELL ≤1e-10。INS: 軌道 ≤1.3e-11、Twiss ≤2.6e-7 (L5)。初期軌道オフセット付き INS はチューン不一致 |
| qmdiag.f | 4x4 対角化 | lattice.py | 翻訳済み・照合済み | ≤1e-10 |
| tquade.f, tbende.f, tdrife.f, tmulte.f | 6D 転送行列 (emit 用) | emit6d.py | 翻訳済み・照合済み | Emittance[] ≤1.3e-11 (L16B〜E) |
| tthine.f | 薄レンズ行列 | emit6d.py | 翻訳済み・照合済み | ≤1.3e-11 |
| tcave.f | 空洞行列 | emit6d.py | 翻訳済み・照合済み | ≤1.1e-7 (L4、RADCOD ON) |
| tsole.f (tsole1)、tsolque.f、tqente.f (SOL 領域) | ソレノイド領域の写像・放射 | solmap.py (Twiss)、solemit.py (Emittance[]・EMIT) | 追跡カーネルの数値微分で翻訳。Twiss は領域内の β・α・分散・結合 1e-11、Emittance[] は εx 2e-12〜4e-9・εy 4e-10〜2e-8 で gs と照合済み。未対応 (両モード): CALC6D・Touschek。k64 モードで未対応: 領域内の MULT・BEND・CAVI、SOL の DX DY DZ CHI・F1 (Twiss・Emittance[]・追跡)。v2 モードはこれらを実装 (ローレンツ力の積分との一致。gs との比較は無い。Emittance[] の POL の量は領域内だけ K64 の式)。数値微分なので解析行列の gs と 1e-10 程度の差が残る | test_sol_optics.py、test_sol_emit.py、docs/29、docs/44_sol_guide.md |
| tchge.f, qchg.f | ずれ・回転座標変換 | lattice.py, emit6d.py、match6d.py `qchg` | 翻訳済み・照合済み | ≤1e-8 (L3B, L3C) |
| tturne.f | emit 用 1 周ディスパッチ | emit6d.py `_turn` | 翻訳済み・照合済み | ≤1.3e-11 |
| qdcell.f, qdquad/qdbend/qdmdia/qdtwis.f, tffsqu | 解析微分 | match6d.py | 翻訳済み・照合済み | マッチングの反復列・K1 を 12 桁再現 (L22A)。数値微分とは 2e-8〜4e-7 |
| qcell6d / temit plot 経路 | CALC6D | calc6d.py | 翻訳済み・照合済み | 1.7e-10 (L18A/B/C) |
| qins.f, tins.f, tinse.f | INS 要素 (領域の入口・出口を担う) | lattice.py の INS の扱い (`ins_twiss`) | ファイル名は pysad に出てこないが、INS を含む 3 セルのリング (BEND E1・E2 付き) の Emittance[] のチューン 1e-15、TrackParticles 1e-13 が gs と一致 (RFSW・RAD 込み)。追加の実装は不要と判断 (2026-11-16) | docs/history/work_unimpl.md |
| qquad, qthin, trad, twelm, pfcoup 等 | 非リンク・デッド | — | 対象外 | — |

## 3. エミッタンス

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| temit.f (tecalc, synchm) | 6D 平衡ビーム行列 | emit6d.py `emittance6d` | 翻訳済み・照合済み | 全量 ≤1.3e-11 (L16)、RADCOD ON 1e-7 (L4) |
| INTGRL.f | RADINT (I1〜I5, 結合, 色収差) | radint.py (FFS RADINT)、lattice `_radiation_integrals_bend` (BEND のみの別実装) | k64 モード: 翻訳済み・照合済み (K64 の癖: 最後の BEND の L/ANGLE を全 BEND に使う、6 極が直前の DL を使う、定数が単精度、CAVI は最後の値)。v2 モード: その BEND 自身の L/ANGLE、6 極は自分の L と K2、倍精度、VOLT は総和 (B06・B07・B12・B13)、量子寿命 (B72)、線形クロマティシティは dν/dδ (B73)、MULT・SOL が入らないときは表示 (B74)。RADINT の Jx が正しく、K64 の Emittance[] の Jx は端の角で誤る (B75。v2 は Emittance[] を直した)。MomentumCompaction の誤り (B77) も v2 で直した。MULT・SOL・BEND の K0 は積分に入らない。結合のあるリングの V 面 (I4V) は未確認 | 7 桁 (L1、L29A)。v2: 解析式・自前の積分と 2e-5〜1e-13 (test_v2_A_radint.py) |
| temit.f の INTRA, tintraconv | IBS 付き emit / Touschek | touschek.py `touschek_table`, `touschek_lifetime` | 翻訳済み・照合済み | 15 ケース ≤7e-12 (L17A)、RADCOD ON 1.7e-8 (L17B) |
| tfemit.f, emit.n | Emittance[] のオプション (OneTurnInformation, Matrix) | emit_extra.py | 翻訳済み・照合済み | 2.1e-10 (L21A)。SOL を含むと使えず `run_emit` の値に落ちる (Stable=0、HarmonicNumber=NaN)。軌道のある QUAD で Python 版 `_tquade` が 1e-5 ずれる (history/work_unimpl.md 不具合 1) |
| temitf.f | EMIT 報告 (temitprint ほか) | ffscmd (EMIT 報告) | 翻訳済み・照合済み (報告の行を gs と比較: bench_level26a、29c) | docs/07 |
| tfemit.f の tfemits、temits.f (temits1, tefsetup, tesolvd, tevdif, tesumb など) | SynchroBetaEmittance | synchrobeta.py、`ffscmd.fn_synchrobeta` | 翻訳済み・一部照合 (12 セルのリング + CAVI で emitx が相対 4e-7〜1e-5、emity は 1e-18 の雑音で符号も違う。24 セルの平面リング + CAVI VOLT=2e6 の 4 点で emitx 1.4e-7〜1.4e-6 (local/wake/sb11.sad)。結合のあるラティスは gs が SIGSEGV で比較不能、不安定なチューン、emiout、RFSW の組み合わせは未確認) | test_synchrobeta.py (許容 3e-5) |
| calpol、spnorm、srequpol、spdepol (photons.f の tradke のスピン部、tchge.f の srot) | POL (偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列) | spin.py、emit6d / emit_extra / elkern の spin 分岐 | 翻訳済み・照合済み (平面リングの PolarizationTime は全桁、縦ずれ・傾き・6 極のリングは 13 桁。MULT・空洞・SOL (ドリフトのある領域) まで)。`Polarization` は K64 が NaN (B55) で K64 の式のまま。SOL 領域内の MULT・BEND・CAVI のスピンは未確認 | test_spin.py、docs/36、docs/43_pol_guide.md |

## 4. 放射

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| photons.f tradk / tradkn / tradkfn | 追跡の放射 | track6 `_tradk`, sadrandom.py | 翻訳済み・照合済み | ≤4.9e-11 (L10) |
| photons.f tradke | emit 用放射 | emit6d.py | 翻訳済み・照合済み | ≤1.3e-11 |
| synradcl.f srfuncy | 光子エネルギー逆 CDF | track6 に内蔵 | 翻訳済み・照合済み | L10 の乱数状態一致 |
| tsdrad, tserad, tstrad, tqrad, trad | Obsolete | — | 対象外 | — |

## 5. COD 補正

| Fortran | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| tcod.f | 6D/4D COD、RADCOD | lattice.py, emit6d.py (6D Newton) | 翻訳済み・照合済み | 9〜10 桁 |
| tffsmatch.f の変数の符号・重み (UNIPOL・RELW・BIPOL) | FIT のフラグ | ffscmd `matcher_` が match6 の bipol・absweit につなぐ | 翻訳済み・照合済み (UNIPOL は変数の符号を変えない、RELW は重みが変わる。fit1.sad の 3 通りで gs と全桁一致)。BIPOL の効果は FIT だけ。GAUSS・CTIME・CONV・EMIOUT はソースで効果を確認 (GAUSS・CTIME・CONV は MAIN の TRACK・MEASURE、CPU 時間の表示、収束の表示。EMIOUT は EMIT 報告の追加行で、pysad は spin の行だけ対応) | test_item4.py |
| tfdisp.f (DISP のモード GA・Z・P・A・B・D と REFERENCE・DREFERENCE) | 光学関数の表 | ffscmd `_disp_show` | 翻訳済み・照合済み (gs と見出し・桁が一致。DISP B の REF は gs も未実装。雑音 1e-16 の列は 0) |
| tffsa.f tfsetref、tftwiss.f の R・D 接頭辞 | REFERENCE、Twiss["RBX"/"DBX"/"RALL"/"DALL"] | ffscmd `_ref_rows`、`twiss_value`、`fn_twiss` | 翻訳済み・照合済み (全桁) |
| tfinit.f tffsrenumber/tfinimult | RENUMBER、LINE["MULT"] | ffscmd `_renumber_cmd`、`line_value` | 翻訳済み・照合済み (INDX 属性の初期値は未対応) |
| tffsa.f の USE/CALC/GO の EXPAND・NOEXPAND | 成分ごとの値を残す | ffscmd `use`、`calc`、`go`、`spread_pending` | 翻訳済み・照合済み (NOKEEPEXP は未対応) |
| DUMP コマンド | 要素の DK・DDK などの表 | ffscmd `_gwd` | 翻訳済み・照合済み (パターンは 20 文字まで、DDK は K0 のある要素だけ。誤差入りリングで値は gs と一致。DK・DDK 列の書式 (Fortran の 1p g12.5: 0 → `0.0000`、小さい値 → `4.22829E-04`) を gs に合わせた) | test_ffsextra.py `test_dump_format_matches_gs` |
| MATRIX コマンド | 2 点間の転送行列 | ffsextra | 翻訳済み・照合済み (名前だけ (gs は NAME.1 として扱う) でも位置に指定でき、ラティスに無い名前を gs が読み飛ばす (to=`$$$`) のも再現) | test_ffsextra.py (mt1・mt2.sad) |
| LDSolve.n、行列演算 (Dot・Inverse) | Optics.n・LDSolve.n の線形代数 | pysad の式評価 | 翻訳済み・照合済み (`DiagonalMatrix` は複素数の列を受ける。行列の演算結果は要素ごとに虚部 0 の複素数を実数に戻す (gs と同じ)。LDSolve の終点と 3 通りの呼び出しが gs と一致。TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo も表示桁で一致)。PeriodicR は呼び出し方が不明で比較できず | test_evalfixes.py |
| OpenRead (開けないファイル) | 入力ファイルの読み | pysad の IO | gs と同じ (開けなければメッセージを出して入力行を打ち切り、Read・Close は wrongtype) | test_evalfixes.py |
| tcorr.f、tvcorr.f (DELCOR) | 相関付き誤差の補正 | ffscmd `_delcor` | 翻訳済み・照合済み (チューン 1e-14、軌道 1e-10、乱数列も同じ)。k64 モード: MARK は DETR・DPX に入る (K64 の添字のずれ B40)。CAVI は K64 のとおり HARM・PHI に書き gs と一致、OCT はチューンと QF の DX が gs と 12 桁で一致 (OCT 自体の DX は gs も pysad も 0 で、効果そのものの比較ではない)。corlength≠0 の枝は K64 が cv(n) を設定しない (B62。2 回目以降の corlength≠0 の DELCOR は再現できない)。v2 モード: MARK・CAVI も DX・DY に入れる (B40)、波数 −1 の成分を対称に置く (B62)、corlength≠0 の規格化を rms = delx に直す (B99。作者の文書による裏づけは無い推定)。v2 は統計・解析式で確認し、gs との比較は無い | docs/31_delcor.md、test_delcor.py |
| tcftr.f (FFT) | 実 FFT | csr.py (CSRTrack 用)、ffscmd・maintrack | 翻訳済み (CSR 側は gs と照合) | docs/36 |
| trcoda.f (TRCOD) | 追跡による COD | measure.py (v2 モード)。k64 モードは ffscmd が受理のみ | v2 モードで翻訳 (gs で SIGSEGV: B42。gs との比較は無い。上の MEASURE の行) | docs/03 |
| terror.f (DELX/DK/DTHETA…) | 誤差コマンド | ffscmd `_error_cmd` | 翻訳済み・照合済み | L3C/3D、L7-2 |
| MON/MCOD/TCOD コマンド | COD 補正 | — | 対象外 (K64 の tffsa.f にこれらのコマンドは無いことをソースで確認 (2026-11-16)。FFS では式として評価されて何もしない) | — |
| tsvdm.f 等 | SVD | match6 の tfsolv が numpy で代替 | 翻訳済み・未照合 (numpy 流) | マッチング K1 ≤1e-10 |

## 6. FFS・マッチング

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| tffsa.f | コマンドループ | ffscmd.py | 翻訳済み・照合済み | L13, L14 |
| tffscalc.f, tffsmatch.f | 残差・GO (Newton/DESCEND, fuzz, fmincube) | match6.py | 翻訳済み・照合済み | 反復数・残差・K1 ≤1e-10 (L13A/B)。FITP>1 は参照値なし (L15A) |
| tffsqu / qdcell | 解析微分 | match6d.py | §2 参照 | L22A |
| tgetfv.f, tfgetv.f, tflag.f, tfsetv.f | FIT 条件・要素値・フラグ | ffscmd.py | 翻訳済み・照合済み | L13, L14 |
| tftwiss.f, tfchro.f, tfoptics.f | Twiss[]/LINE[]/Element[]、色収差、CalculateOptics | ffscmd.py, lattice.py | Twiss[]/LINE[]: 照合済み (L14A、PolynomialFit 1e-7)。CHROMATICITY (tfchro)、CalculateOptics (tfoptics; 4 次元の cell/INS・区間・結合・分散、idim=3 は全周): 翻訳済み | — |
| tfgeo.f, tdgeo.f, tgrot.f | 幾何 (GX…、ORG) | geo.py (LINE["GX"] / fit GX.. / DISP G 経由)、match6d.py `tdgeo_term` | 翻訳済み・照合済み | 4.8e-13 (L19A/B/C)。tdgeo は DRIFT の L・BEND の ANGLE・L (幾何の条件だけの行) を解析微分に翻訳 (有限差分と 1e-5 以内、test_tdgeo.py)。四極などの DX・DY、Twiss の条件が同じ行にある列は中心差分のまま (半セル ARCAH の例で収束位置が gs と相対 5.5e-7 ずれる) |
| tfdapert.f | DynamicAperture[] | dapert.py | 翻訳済み・照合済み | L16A |
| tffamsetup.f | FITP の有限振幅ファミリー (MatchingAmplitude/NFAMP) | match6.py `fam_setup` ほか | 翻訳済み・未照合 (SAD が FITP>1 の GO で再起動するため自己検証のみ) | test_fam_points.py |
| tffswake.f | ウェイク (`tffssetupwake` = 追跡用の表の取得、`tffswake` / `tffswakekick` = 光学) | ffscmd.py `_wake_table` `ini_optics`、wakeopt.py | 追跡・光学とも翻訳済み・照合済み (InitialOrbits = tffsinitialcond、tfshow の軌道列も) | test_wake.py、test_wake_optics.py |
| tmulti.f tmultiacc | RF のある MULT (TRPT か RFSW) | track6n_ext.py `mult_full` | 翻訳済み・照合済み (ANGLE=0、F1・FRINGE なし、AUTOPHI なし) | test_wake.py |
| tfstk, tfeval, tfetok 等 (SADScript) | インタプリタ | sadscript.py、sadlist.py (リスト系組み込み) | 翻訳済み・照合済み | L14A |

## 7. ビームライン定義・座標

| Fortran | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| toplvl / doelem / doline / expln | MAIN 部 (LINE, 繰り返し, 再定義) | parser.py | 翻訳済み・照合済み | L12 (3563 要素) |
| talign.f, nalign.f | ずれ・回転の設定 | ffscmd, lattice | 翻訳済み・照合済み | L3B/3C |
| tfgeo.f | フロアー座標 | geo.py | §6 参照 | 4.8e-13 |
| tfloor.f | 数学関数 (cfunc) | 各所に個別実装 (xsincos 等) | 翻訳済み・未照合 (個別) | — |
| tsgeo.f | ソレノイドの幾何 | sgeo.py | 翻訳済み・照合済み (ずれ・傾きのある SOL の幾何)。k64 モード: 領域内の BEND・MULT・CAVI、領域内の L=0 / K1=0 の四極は未対応 (NotImplementedError)、tdgeo の領域幾何は未翻訳。v2 モード: 領域内の BEND・MULT・CAVI・K1=0 の四極も軌道を追い、領域内の変数の幾何 fit は中心差分 (設計軌道の粒子が出口で 0 に戻る 1e-12。gs との比較は無い) | test_sol_geo.py |
| phsrot.f | PHSROT (位相空間の回転) | phsrot.py | k64 モード: 未対応 (K64 で動作しない: 行列がすべて 0。追跡は NotImplementedError、光学は恒等。docs/21 B38)。v2 モード: 意図された写像 M = N·Rot·(1−D)·N⁻¹ + 拡散を実装 (追跡・CALC・Emittance[]・EMIT。固有値・シンプレクティック性・平衡エミッタンスで確認。gs との比較は無い。K64 の拡散の乱数の誤り B79 も直した)。スピンは動かさない | test_v2_C_phsrot.py |
| tpara, tscale 等 | 空・デッド | — | 対象外 | — |

## 8. 寿命・その他

| Fortran | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| Touschek (tintra 系、emit.n) | TouschekTable/Lifetime | touschek.py | §3 参照 | ≤7e-12 |
| tflifetrack.f | LifeTrack | — (SAD 側も空スタブ) | 対象外 | — |
| tlum.f, tluma.f | ルミノシティ (MEASURE の統計) | measure.py | 翻訳済み・未照合 (MEASURE が gs で動かないので解析値との比較だけ。ffscmd からは呼ばれず、test_measure.py のみ。BEAMBEAM の LuminosityMonitor は beambeam.py) | test_measure.py |
| txwake.f | ウェイクのキック (TrackParticles) | wake.py | 翻訳済み・照合済み (普通の要素、CAVI、SOL の領域の内側) | test_wake.py |
| tfwake.f、twinit.f 他 | 旧 WAKE コマンド | — (同じ経路で動く) | 移植不要: K64 の tffsa.f:972 は tfwake.f を呼ばず WAKECOMMAND[] (配布に無い) を評価する。gs はエラーを出して行の残りを捨て、pysad も同じ経路 (エラー文の表示だけ未再現) | docs/23 |
| wjfit、wfres、wfbin、wiord、wjdraw | MAIN の TRACK の Fourier 解析・共鳴の描画 | — | 未翻訳 (FFS から使えない。maintrack.py の `FOURIE ON` も未対応) | — |
| tspac.f (GaussianCoulomb) / tintrb.f (twspu、twspfu) | 空間電荷関数 | gausscoul.py | GaussianCoulomb・GaussianCoulombU・GaussianCoulombFitted を翻訳し gs と 1e-9 以内で一致。`tspac.f` の `tspac` は K64 のどこからも呼ばれない死んだコード。追跡の空間電荷 (SPAC・PSPAC・WSPAC) は v2 モードで `spch.py` に実装 (上の行) | test_gausscoul.py、docs/32 |
| tlinit.f (tlfield、tlcalc、tlspect) | RadiationField・RadiationSpectrum | radfield.py | 翻訳済み。RadiationField は gs が SIGSEGV (B49) なので等速運動の電荷の速度場と相対 6e-3 の確認のみ。RadiationSpectrum は gs と一致 (NPARA の fork 並列は逐次) | test_radfield.py、docs/45_radlight_guide.md |
| gamma.f、hg.f、tfunc*.f (特殊関数)、tffindroot.f、D.n | 特殊関数・Fit・FindRoot・NIntegrate | specfn.py、fitroot.py、miscfn.py | 458 式が gs と一致 (K64 が誤値を返す B43〜B47、B53、B54 は真値)。Fit・FindRoot は数値微分で 5e-10〜1e-9。NIntegrate は複素数の被積分関数を扱い (gs と虚部 1e-15)、`LogGamma1` もある。複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)、複素引数の Bessel・超幾何・Airy (実数のパラメータ) と PolyLog の非整数次の大きい \|z\| (相対 1e-11、\|z\| ≤ 1e12) は両モードで実装。k64 と値が違う関数は InverseErf (B47) だけで、v2 は真値。複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で両モードで実装 (`specfn_cx.py`。gs が誤る点でも真値: B118〜B120) | docs/33 |
| ActLie, ActGRA, tdlat, tmap, polpar 他 | 空スタブ・描画・デッド | — | 対象外 | — |

## 9. 未翻訳・一部未対応・未照合の一覧 (2026-11-18、v2.1.0 の後)

最新の一覧は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) にあり、この節はその要約 (項目を直したら 2 つを一緒に更新する)。
各行の「モード」は **k64** (K64 (gs) と同じ動作)、**v2** (v2.0.0 から既定の、K64 のバグを直し機能を足した動作)、**両** (モードによらない)。v2 の機能の根拠は解析式・数値積分で、gs と比較できないものが多い (gs が動かない・K64 が誤っている)。
実装・解消済みのもの (APERT 要素の追跡 (SOL 領域内を除く)、FFS 内の `MOMENTUM=`・`[NO]EXPAND`・RENUMBER・DISP の各モード・REFERENCE・DUMP・OrbitGeo、MakeBump・Couple の ElementValues・Rational・NIntegrate の複素数・LogGamma1・ParametricPlot・ListContourPlot・OptimizeOptics・Entry/ListBox/Scale のメソッド、LOSSMAP・RADTAPER・MULT のスキュー成分・ResponseMatrix の TRPT・NOKEEPEXP・RENUMBER の INDEX・`BeginPackage`/`Begin` の文脈・`Shared`。**v2.1.0 で**: Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE と ANGLE=0 の BEND・MULT の F1、光学・Emittance への空間電荷 (WSPAC)、光子 (tphrec) の BEND・QUAD・MULT・SEXT・OCT・ドリフト・BEAMBEAM + BEND、tdgeo の領域内の四極の DX・DY の解析微分、CALC6D・RADCOD つきの CALC・INTRA と SOL、傾いたソレノイドの POL の量、POL + TRPT + 加速 (領域の外)、`MASS`・`CHARGE`、`NPARA=n;`、複素の PolyLog・HypergeometricU・Zeta、`Parallelize` (POSIX)、成分ごとの FREE、`$Context` の代入・`WithContext`、VISIT/BYE、BEND の F1 + RAD の追跡・CAVI の DX・DY・EMIT のあとの USE の差) は、ここでは未対応に入れない (46 §6)。

### 9.1 未翻訳 (理由つき)

| ファイル・機能 | モード | 理由 |
|---|---|---|
| wjfit、wfres、wfbin、wiord、wjdraw | 両 | MAIN の TRACK の Fourier 解析で、FFS から使えない (gs も SIGSEGV)。チューンの抽出は [`pysad/naff.py`](../pysad/naff.py) (v2.1.0) |
| 光子 (tphrec) の SOL 領域内の要素・COORD・TCAVI・MAP のあるラティス | 両 | BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND は対応。SOL 領域内は記録しない |
| TrackParticles の F1 (`tsteer`・`tmulti`) | 両 | ANGLE=0 の BEND・MULT。Emittance[]・EMIT・CALC は実装済み |
| CALC (4 次元) の RADTAPER | 両 | K64 の qtwiss にテーパーの式が無い (実装しない) |
| MON/MCOD/TCOD の機能、tdrwdt.f の粒子追跡 | 両 | K64 にコマンドが無い / 死んだコード (B105)。v2 は警告して読み飛ばす |
| tspac.f の `tspac` | 両 | K64 のどこからも呼ばれない死んだコード (追跡の空間電荷は v2 で `spch.py`) |
| `Out[n]:=` の gs の表示 | 両 | K64 の癖として写さない |
| 複素の特殊関数で gs が誤る点の k64 での再現 (B118〜B120) | k64 | K64 の gamma.f の複素の不完全ガンマなどの直訳が要る。k64 も真値を返す (要判断) |
| `DAPWIDTH` の窓 (既定 7)、K64 の `qwspac` (4D の CALC の空間電荷) の直訳 | 両 / v2 | trackd.f の追跡順の直訳 / gs が NaN を返すので検証できない |

### 9.1a モードで状態が違うもの (k64 モードでは未対応・v2 モードで実装、または結果が違う)

| ファイル・機能 | k64 モード | v2 モード |
|---|---|---|
| trackb.f、tracke.f、trcoda.f、ttinit.f、ttstat.f (MEASURE / TRCOD ほか) | 未対応 (gs で TRCOD が SIGSEGV、MEASURE は出力なし: B42) | 実装 (`measure.py`。K64 のソースの読み。**gs との比較は無い**) |
| tsole1 の領域内の MULT・BEND・CAVI (tmulte1)・APERT、SOL の DX DY DZ CHI・F1 (Twiss・Emittance[]・追跡)、領域の tsgeo・tdgeo | 未対応 (領域内の BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げた: B69) | 実装 (ローレンツ力の積分との一致。gs との比較は無い)。領域内の四極の DX・DY は両モードで gs と全桁 |
| 領域内の空洞の TRPT の加速・wake・AUTOPHI | 未対応 | 実装 (AUTOPHI は TRPT で gs と全桁一致、ほかは運動方程式の積分) |
| 空間電荷 (SPAC・PSPAC・SELFCOD の追跡。spkick.f・tspch_.c) | 追跡に入らない (K64 の SPAC は 1/γ 倍に弱い: B89 ほか) | 実装 (`spch.py`)。WSPAC の Emittance[]・EMIT は k64 も gs と全桁、CALC・CALC6D は v2 だけ |
| PHSROT | K64 で動作しない (B38)。追跡は NotImplementedError | 実装 (`phsrot.py`) |
| UND の磁場 | K64 が設定できない (B37) | `FBX`・`FBY` を実装 |
| Polarization (平衡偏極度) | K64 が NaN (B55)・率が約 1000 倍 (B125)・逆向きの BEND で相殺 (B116) | 有限値 (側帯の項つき。追跡の 0.43 倍の過小評価が残る) |
| Emittance[] の POL の量 (傾いたソレノイド) | K64 の式 (B130) | v2.1.0 で修正 (追跡と 1e-13〜5e-7) |
| RADINT の K0 の極・四極のずれ・BMAG (B115・B126・B107)、ANGLE=0 の BEND の端の角 (B108)、tradke の符号 (B111)、WSPAC のレンズ (B127)、K1≠0 の BEND の追跡の端 (B117)、ずれた BEND の CALC (B114) | K64 のまま | 修正 (47 の v2.1.0 の (a)) |
| System、Fork | `System` は実行せず、`Fork` は pid を返すだけ | `System` は許可があれば実行、`Fork` は POSIX で本当に分ける |
| 旧 WAKE コマンド (tfwake.f) | 移植しない (K64 は WAKECOMMAND.n を呼ぶ。B66) | `WAKE … ENDWAKE` の表を読む |
| 複素 FindRoot | TypeError (gs は壊れた値: B61) | 複素の Newton 法 |

### 9.2 一部未対応

| 機能 | モード | 未対応の部分 |
|---|---|---|
| MULT ANGLE≠0 (multa / multae) | 両 | Twiss の BZ・RF・F1・SK0・FB1/FB2 (追跡は FB1/FB2 まで対応)。MULT の端の FB1/FB2 の tblfre と AE1・AE2 (Emittance[]) |
| MAP | 両 | GEO・CompiledMap。TRACK は TrackParticles だけ (MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない)。MAP のスピン |
| POL (SOL 領域内) | k64 / v2 | k64: 領域内の MULT・BEND・CAVI のスピン、スピン追跡は損失の処理なし。v2: 領域内の APERT のスピン追跡 (NotImplementedError)、領域内の四極のずれの POL の追跡 (track6)、bz ≥ 3 1/m の QUAD・MULT のスピンは 1e-4 台、加速する CAVI の POL (Emittance[]) |
| POL + TRPT + 加速 (track6) | v2 | SOL 領域の外の CAVI だけ。領域内・RF つき MULT・TCAVI・リングの何周もの TRPT は SADError |
| Polarization | v2 | 4.2 GeV・12 MV で FLUC つきの追跡の 0.43 倍に過小評価 (軌道の非線形)。運動論的偏極なし。不安定なリングは線形の式。n が傾く・非平面のリングの緩和率は未確認 |
| TCAVI | 両 | Emittance[] は L=0 のみ gs と一致 (L≠0 は gs が壊れる: B36)。RANKICK・RANPHASE 未対応 |
| UND | k64 / v2 | k64: 磁場なしのドリフトのみ。v2: 端の磁極の整形なし。挿入光源のあるリングでは未確認 |
| PHSROT | v2 | スピンを動かさない |
| APERT | k64 / 両 | k64: SOL 領域内は NotImplementedError。両: JDPX・JDPY の再入射は粒子 12 個以上で gs と乱数列が違う (`NPARA=n;` で分け方を再現できる)。RAD・FLUC つき、TRPT、MAIN の TRACK の損失報告、DAPERT での APERT は gs と未照合 |
| LOSSMAP | 両 | MAP 要素・空間電荷が落とした粒子の要素番号は 0 |
| RF のある MULT (tmultiacc)・CAVI の V1 V20 V11 V02 | 両 | F1・FRINGE。WakeFunction のある CAVI・領域内の CAVI・Emittance[] の CAVI の V1 V20 V11 V02、ANGLE≠0 の MULT。AUTOPHI は k64 で NotImplementedError・v2 で実装 |
| 空間電荷 | v2 | ANGLE≠0・VOLT≠0 の MULT と SOL 領域内の MULT の SPAC、QUAD・BEND の中の SPAC、複数バンチ、MAIN の TRACK・DA・スピン追跡。追跡の WSPAC のビーム行列は空間電荷なしのまま。PSPAC の粗いメッシュは警告 (gs と違う)。ウェーク・放射との組み合わせは未確認 |
| MEASURE・TRCOD | v2 | gs との比較なし。DETR ≥ 1 (モードが入れ替わる) の MEASURE は未確認。FIFO での待ち合わせ未確認 |
| RADLIGHT | v2 | SOL 領域を含む区間、COORD による幾何の変化、曲げのある MULT の中の点 |
| RADTAPER | 両 | SOL 領域内の要素、RFSW OFF。RADTAPER の重みは「閉じる」を根拠にした選択 (B111) |
| SOL | 両 | BOUND=0 で BZ≠0 の SOL の Twiss。v2 は SEXT・COORD・BEAMBEAM が領域内に置けない (`SolRegionError`)、領域内の BEND の K1・DX・DY。CALC6D・RADCOD つきの CALC・INTRA と SOL は対応 (TouschekLifetime は未確認) |
| DELCOR | v2 | corlength≠0 の規格化 (B99) の意図は推定 (根拠は確定) |
| RADINT・Emittance[]・EMIT | v2 | RADINT は MULT・SOL・ANGLE≠0 の BEND の K0 が入らない。K0 の極の εx は 3.5%・Jx は 0.7% Emittance[] とずれる。SEXT・BEND のずれ。K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角 (AE)。v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う (未調査) |
| SynchroBetaEmittance | 両 | 大きなリング A で 40 倍違い収束しない (B128。v2 は警告) |
| Fit・FindRoot・特殊関数 | 両 | 複素の次数の PolyLog、Re a ≤ 0 の複素の a・b の HypergeometricU。複素 FindRoot の変数の範囲 (無視) |
| `ElementValues` (Couple.n)・FREE の成分 | 両 | フィット変数との連動 (`CouplingMatrix`・`CoupledKeys`)。先頭でない成分の既定のキーの変数は作らない |
| 文脈 (`BeginPackage`) | 両 | 文脈つきの `Names`、init.n が大域に作る名前の全部 (51 個しか調べていない) |
| 描画・GUI | 両 | CanvasSymbol、Canvas3D*、`UseBG->True` の文字背景。OpticsPlot の text 項目数 (gs 69、pysad 74。Tk のフォントによる)、BeamPlot (gs は SIGXX=0 で楕円なし)、最初の `Window[]` の親 (gs は Toplevel)、`Packages/CSR.n` 全体は読み込まない (docs/37)。Tk 9.0 の `cget` の単位つきの値 |

### 9.3 未照合 (コードはあるが gs と数値比較していない、または gs が参照を出さない)

| 機能 | モード | 理由 |
|---|---|---|
| v2 モードの機能全般 (MEASURE・TRCOD、SOL 領域内の要素、UND・PHSROT、空間電荷の SPAC、Polarization の有限値、RADLIGHT の 3 次元の軌道、v2.1.0 の修正の値) | v2 | 値の根拠は解析式・独立な数値計算。gs との比較は無い (gs が動かない・K64 が誤っている)。k64 モードの再現は gs と比べた |
| SP8 (大きなリング B。SPring-8 型、デザインレポートより再現。8 GeV、周長 1436 m、矩形の偏向磁石 88 台) | 両 | 数値だけ確認 (CALC・RADINT・Emittance[]・COD は k64 が gs と 1e-10)。U0/E が大きいリングの RADTAPER・tradke の符号、陽子の大きなリング、POL (不安定で NaN) は未確認 |
| FITP>1 の最終解、tffamsetup (fam 点) | 両 | FITP>1: pysad は独立な最小化と 5e-8 で一致、gs は最小値の手前で止まる。fam 点: 初期条件の不変量と8 極のキックのチューンは確認、有限振幅の周期解の条件そのものは未確認 |
| CSROYMatrix・OYStability・OYScan | 両 | 行列は手作りの入力で gs と 4e-20。gs は `CSRSetupOY` の出力で wrongtype (原因は zj の head の上書き: B60)。OYStability・OYScan は gs が完走しない |
| tbrote の srot 更新、tbendebody の bsi2 (B56・B58) | k64 | ソースの読みだけで、gs の個別確認なし (v2 は直した) |
| RadiationField | 両 | gs が SIGSEGV (B49)。Liénard-Wiechert の厳密な数値と合う (振幅・位相・スペクトル・全エネルギー) |
| measure.py (tlum・tluma) | 両 | gs で MEASURE が動かない。解析値のみ |
| tfloor.f の cfunc (xsincos 等) | 両 | 個別実装 |
| MASS・CHARGE を与えた粒子の Touschek・Intrabeam・RADINT、放射の β の因子 | 両 | 電子の古典半径のまま (B112)。陽子での確認なし |
| `Twiss["DETR"]` ≥ 1 | 両 | モードが入れ替わる INS・MEASURE は未確認 |
| 未解決の差: SOL 内の四極 + RAD (現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い。B68)、K3 の強い MULT リングの qcod (Emittance が NaN)、PSPAC の x と y の 16% の非対称 (再現しない)、v2 の Emittance[] (RADCOD) と gs の 1e-4〜1e-3 (意図した差: B111) | 両 | 原因の調査が未了または意図した差 (docs/11 §E、`history/work_unimpl.md`、docs/46 §4) |

### 9.4 優先度

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | Monitor/Steer/CalculateOptics は動く。誤差入り 6D emit と DA の反復が 600 s を超える |
| 中 | Polarization の追跡との差 (非線形の軌道)、v2 の Emittance[] と RADINT の矩形磁石の差、SynchroBeta の大きなリング、9.3 の未解決の差、v2 の要判断 (docs/47 v2.1.0 の (d)) | 原因の切り分け |
| 中 | TrackParticles の F1、MULT の AE1・AE2・FB1/FB2、SOL 領域内の光子 | K64 のソースの直訳 |
| 低 | wfres 系 | 必要になったら (NAFF で代用) |
| 不要 | MON/MCOD/TCOD の機能、Obsolete・非リンクのファイル | K64 で無効 |

## 10. `<SAD>\src` の *.f のうち pysad の中に名前が出てこないもの

`<SAD>\src\*.f` 474 本を `grep -w -i` で `pysad/*.py` と突き合わせた (2026-11-16)。名前が出ないだけで、別名のルーチンに直訳されている場合がある
(例: `tins`・`tinse` は `lattice.py` の INS の扱い、`talign`・`nalign` は ffscmd / lattice のずれ・回転)。
物理計算で名前が出ないもの: `qins tins tinse` (INS、§2。実装は別名で gs と一致)、`qmult qddrif` (qtwiss・qdcell 内の別実装に吸収: docs/11 D5)、`tfmap tmap` (MAP は別名)、
`tflifetrack twinit wfbin wfres wiord wjdraw` (§1・§8・上の 9.1。`trackb` は v2 の `measure.py` の MEASURE の実装で名前が出るようになったので、この一覧の再現 (2026-11-17) では出ない。名前の再数えはしていない: 要確認)、`tbrad tquads tquase tsdrad tserad tstrad tqrad trad trade` (Obsolete)、
`tconv` (追跡の変換は track6 内)、`tmuld tmultr` (行列の積は numpy)、`temitf` (EMIT 報告は ffscmd。ファイル名としては出ない)。
それ以外は、旧 MAIN のユーティリティ (`pmbump pwrite pundo pkill` など)、SADScript のインタプリタ内部 (`tfeval1 tfstk tfpart tfsort tfreplace itf*` など)、描画 (`geodrw phdrw pmbdrw tdplt twsdrw`)、ファイル・バッファ (`filaux filbuf flmgr getbuf`) で、
pysad の方針 (SADScript は Python のインタプリタ、ファイル・描画は Python のライブラリ。`15_sadhelp_purpose_policy.md`) により翻訳しない。
一覧の再現: `for f in <SAD>/src/*.f; do b=$(basename $f .f); grep -q -w -i -- "$b" pysad/*.py || echo $b; done`。
