# 16. SAD Fortran ファイル別 pysad 翻訳状況

最終更新: 2026-11-17 (現状のスナップショット。内容は `history/work_unimpl.md` と `pysad/*.py` の docstring で裏を取った。§9 は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) と同じ状態)

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
追加分 (2026-11-13 以降の要素・機能別): coord / tcav / und / multa / multae / solmap / solemit / spin / csr / synchrobeta / wake / wakeopt / beambeam / gausscoul / bbbrem / measure / radfield / specfn / fitroot / miscfn / vecfn / fileops / corr / pkgload / tclbridge / tkcanvas。
(`pysad/*.py` は 67 本。`tracking.py` の ttcav 等への言及は旧 API の名残で、実装は `tcav.py`)
`ffs.py`, `matching.py`, `tracking.py`, `optics.py` は非推奨の旧 API (SAD 未照合の独自実装) で、表には数えない。
calc6d / match6d / touschek / geo / track6n_ext / emit_extra はいずれも `ffscmd.py` から呼ばれる (CALC6D、解析微分、INTRA/TouschekLifetime、LINE["GX"..]/GEO/ORG、MULT/SOL の numba 追跡、Emittance[] のオプション)。

試験の現状 (2026-11-17、`PYSAD_SAD_PACKAGES` つき): `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` で **575 通過・3 skip (578 件収集)**。回帰は **50 ベンチをすべて実行して FAIL 0** (非 QUICK、205 s。QUICK で SKIP される 4 件 (7・12t・13a・16a) も実行して PASS: 5.0e-11、12t 3.5e-11、13a 4.3e-21、16a 1.4e-31)。全試験と `QUICK=1` の回帰 (50 ベンチ中 46 件、FAIL 0、125 s) は 2026-11-17 に実行。QUICK なしの記録は `local/wake/verify_log.md` (2026-11-16 夜。今回は再実行していない)。

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
| tsolqu.f / tsol.f | ソレノイド領域 (bz 分岐、tsolqur) | track6 `_tsol_region`, `_tsolqu_bz_*`、track6n_ext.py | 翻訳済み・照合済み | ≤1.5e-10 (L9)。EPS≠1 は制限あり |
| tphoton / photons.f の tradk | 放射キック (平均・FLUC) | track6 `_tradk`, sadrandom.py | 翻訳済み・照合済み | ≤4.9e-11、乱数状態一致 (L10) |
| tapert.f | 損失判定・APERT 要素 (tapert1) | track6n.py `k_apert_gen`・`k_apert`・`k_apert_el`、track6.py `_apert_pass`、dapert.py、`TrackRun.l_track0` | 翻訳済み・照合済み (APERT 要素: 楕円・矩形・DX DY・ROTATE・全部 0・複数周・JDPX JDPY DP が gs と一致、`tests/test_apert.py`。SOL 領域内の APERT・TRPT・RAD つきは未確認。JDPX の乱数は粒子 12 個以上で gs と不一致) | DA の地図一致 (L16A)。ただし損失周の乱数消費は不一致 (RAD+FLUC)。K64 は落とす印に `kptbl(i,4)=l_track` を書くので、セッション最初の TrackParticles (l_track=0) の最初の tapert では最初から LOSSAMPL の外の粒子が落ちない。この癖を再現し gs (1・2 回目) と全桁一致 (B33 の出力座標の差は解消、docs/11 A52) |
| tsteer.f | BEND ANGLE=0, K0≠0, L≠0 | track6.py `_tsteer_ens` (F1/FB1/FB2 フリンジ付きは NotImplementedError) | 翻訳済み・照合済み | bendk1_* (history/work_bend_k0.md, history/work_trpt_acc.md 10 節) |
| tmulta.f (tmulta, tmultae) | MULT ANGLE≠0 | multa.py (追跡)、multae.py (Twiss・Emittance[]) | 翻訳済み・照合済み (docs/28) | 追跡は相対 0、Twiss 1e-13、Emittance[] 1e-11。BZ・RF・F1・SK0 は未対応 |
| ttcav.f, ttcave.f | TCAVI (クラブ空洞) の追跡と 6 次元写像 | pysad/tcav.py | 翻訳済み・照合済み (追跡 0、Emittance[] は L=0 のみ gs と一致。L≠0 は gs が壊れる) | docs/26_tcav.md |
| tcoord.f, tcoorde.f, qcoord.f | COORD 要素 | pysad/coord.py | 翻訳済み・照合済み (追跡 0、Twiss・幾何 1e-11、Emittance[] 5e-10) | docs/25_coord.md |
| temap.f (temap/temape/qemap) | MAP 要素 (ExternalMap) | elements.Map、ffscmd._map_*、emit6d._tmap、track6n.K_MAP | TRACK・EMIT・OPTICS を翻訳し gs と照合済み (GEO・CompiledMap は未対応) | test_map.py |
| undulator.f | UND (追跡) | pysad/und.py | 翻訳済み・照合済み (追跡 0、Twiss・Emittance[] は恒等で一致)。K64 は磁場を設定できない (B37) | docs/27_und.md |
| bbstrhl.f, bbtool.f, matrix.f (包絡行列) | BEAMBEAM | pysad/beambeam.py | 追跡 (track6n) を翻訳。BSTRL と光子 (PHOTONS、BEAMBEAM の分だけ)、LuminosityMonitor も対応。光学は恒等。傾き tilt の K64 の癖 (B34) は直訳 (残る差は docs/24) | docs/24_beambeam.md、docs/39_beambeam_guide.md |
| photons.f の sprot、tradk のスピン部 | TrackParticles のスピン追跡 (POL かつ RAD ON のときだけ) | track6.py (`sprot`、`bsi`) | 翻訳済み・照合済み (gs と 4e-12 以下: 平面リング、ずれ・回転した四極・偏向・六極、MULT、K1 のある BEND、SEXT・薄い MULT、SOL 内のドリフト)。損失 (flag) の処理なし。SOL 内の四極のスピン追跡は gs の食い違い (docs/44_sol_guide.md §7) のため未確認 | test_spin_track.py、docs/36、docs/43_pol_guide.md |
| tlinit.f (tlinit/tlstore/tlresult) | RADLIGHT の軌道の記録 | ffscmd `_RadLight` | 翻訳済み・照合済み (TrackParticles の戻り `{結果, 軌道{t,x,y,z}}`、gs と差 0。K64 は BEND の中を記録しない) | test_radlight.py、docs/45_radlight_guide.md |
| tbbbrem.f | BBBrem1 関数 | bbbrem.py | 翻訳済み・gs と一致 (docs/32) | test_bbbrem.py |
| csrtrack.f, csroy.f, csrtest.f | CSR 関数 (CSR.n が使う 7 関数: CSRTrack CSRHaissin CSRSetupOY CSROYMatrix CSRMatrix CSRConvert CSRInit) | pysad/csr.py | 翻訳済み。CSRTrack・CSRHaissin は gs と 11 桁、CSRSetupOY は J・ω 1e-9、CSRMatrix (mode 22・33・23・32・110・220・330 と mode 1000 = csrhelmat)・CSRConvert (= K64 の tfcsrvec)・CSRInit は gs と相対 1e-10〜1e-13。CSROYMatrix は gs が wrongtype を返す (B60) ので未照合。`Packages/CSR.n` の MakeZL・SkinDepth・PathDifference などは SADScript のまま動かして gs と一致 (MakeZL は約 7 分、相対 1e-10)。ScanTrack・OYStability などの長い流れは未照合。`CSR.n` 全体は描画 (ListPlot) が要るので読み込まない | docs/36、docs/40_csr_guide.md、test_csr.py |
| track.f, tracka.f, ActTra.f, ttinit.f | MAIN の TRACK コマンド・初期分布 | ffscmd `TrackParticles[]` (SADScript 関数経由)、MAIN の TRACK コマンドは maintrack.py | 翻訳済み・照合済み | TrackParticles は L6〜12。MAIN TRACK は SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 (SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV)。PHSPLOTS の描画は無視、FOURIE ON は未対応 |
| trackb.f, tracke.f, trcoda.f, ttinit.f, ttstat.f | MEASURE / TRCOD | ffscmd が語を受理するのみ。ttstat から呼ばれる tlum・tluma だけ `measure.py` に翻訳 | 未翻訳 (gs で TRCOD は SIGSEGV、MEASURE は出力なし: docs/21 B42。参照が取れず、使うスクリプトも見つかっていない。ttinit の読みの誤り B41) | — |
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
| tsole.f (tsole1)、tsolque.f、tqente.f (SOL 領域) | ソレノイド領域の写像・放射 | solmap.py (Twiss)、solemit.py (Emittance[]・EMIT) | 追跡カーネルの数値微分で翻訳。Twiss は領域内の β・α・分散・結合 1e-11、Emittance[] は εx 2e-12〜4e-9・εy 4e-10〜2e-8 で gs と照合済み。未対応: CALC6D・Touschek、領域内の MULT・BEND・CAVI、SOL の DX DY DZ CHI・F1 (Twiss・Emittance[]・追跡)。数値微分なので解析行列の gs と 1e-10 程度の差が残る | test_sol_optics.py、test_sol_emit.py、docs/29、docs/44_sol_guide.md |
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
| INTGRL.f | RADINT (I1〜I5, 結合, 色収差) | radint.py (FFS RADINT)、lattice `_radiation_integrals_bend` (BEND のみの別実装) | 翻訳済み・照合済み | 7 桁 (L1、L29A) |
| temit.f の INTRA, tintraconv | IBS 付き emit / Touschek | touschek.py `touschek_table`, `touschek_lifetime` | 翻訳済み・照合済み | 15 ケース ≤7e-12 (L17A)、RADCOD ON 1.7e-8 (L17B) |
| tfemit.f, emit.n | Emittance[] のオプション (OneTurnInformation, Matrix) | emit_extra.py | 翻訳済み・照合済み | 2.1e-10 (L21A)。SOL を含むと使えず `run_emit` の値に落ちる (Stable=0、HarmonicNumber=NaN)。軌道のある QUAD で Python 版 `_tquade` が 1e-5 ずれる (history/work_unimpl.md 不具合 1) |
| temitf.f | EMIT 報告 (temitprint ほか) | ffscmd (EMIT 報告) | 翻訳済み・照合済み (報告の行を gs と比較: bench_level26a、29c) | docs/07 |
| tfemit.f の tfemits、temits.f (temits1, tefsetup, tesolvd, tevdif, tesumb など) | SynchroBetaEmittance | synchrobeta.py、`ffscmd.fn_synchrobeta` | 翻訳済み・一部照合 (12 セルのリング + CAVI で emitx が相対 4e-7〜1e-5、emity は 1e-18 の雑音で符号も違う。24 セルの平面リング + CAVI VOLT=2e6 の 4 点で emitx 1.4e-7〜1.4e-6 (local/wake/sb11.sad)。結合のあるラティスは gs が SIGSEGV で比較不能、不安定なチューン、emiout、RFSW の組み合わせは未確認) | test_synchrobeta.py (許容 3e-5) |
| calpol、spnorm、srequpol、spdepol (photons.f の tradke のスピン部、tchge.f の srot) | POL (偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列) | spin.py、emit6d / emit_extra / elkern の spin 分岐 | 翻訳済み・照合済み (平面リングの PolarizationTime は全桁、縦ずれ・傾き・六極のリングは 13 桁。MULT・空洞・SOL (ドリフトのある領域) まで)。`Polarization` は K64 が NaN (B55) で K64 の式のまま。SOL 領域内の MULT・BEND・CAVI のスピンは未確認 | test_spin.py、docs/36、docs/43_pol_guide.md |

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
| tcorr.f、tvcorr.f (DELCOR) | 相関付き誤差の補正 | ffscmd `_delcor` | 翻訳済み・照合済み (チューン 1e-14、軌道 1e-10、乱数列も同じ)。MARK は DETR・DPX に入る (K64 の添字のずれ B40)。CAVI は K64 のとおり HARM・PHI に書き gs と一致、OCT はチューンと QF の DX が gs と 12 桁で一致 (OCT 自体の DX は gs も pysad も 0 で、効果そのものの比較ではない)。corlength≠0 の枝は K64 が cv(n) を設定しない (B62。2 回目以降の corlength≠0 の DELCOR は再現できない) | docs/31_delcor.md、test_delcor.py |
| tcftr.f (FFT) | 実 FFT | csr.py (CSRTrack 用)、ffscmd・maintrack | 翻訳済み (CSR 側は gs と照合) | docs/36 |
| trcoda.f (TRCOD) | 追跡による COD | ffscmd が受理のみ | 未翻訳 (gs で SIGSEGV: B42) | — |
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
| tsgeo.f | ソレノイドの幾何 | sgeo.py | 翻訳済み・照合済み (ずれ・傾きのある SOL の幾何)。領域内の BEND・MULT・CAVI、領域内の L=0 / K1=0 の四極は未対応 (NotImplementedError)。tdgeo の領域幾何は未翻訳 | test_sol_geo.py |
| phsrot.f | PHSROT (位相空間の回転) | — | 対象外 (K64 で動作しない: 行列がすべて 0。docs/21 B38) | — |
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
| tspac.f (GaussianCoulomb) / tintrb.f (twspu、twspfu) | 空間電荷関数 | gausscoul.py | GaussianCoulomb・GaussianCoulombU・GaussianCoulombFitted を翻訳し gs と 1e-9 以内で一致。tspac.f の追跡側 (要素としての空間電荷キック) は未翻訳 | test_gausscoul.py、docs/32 |
| tlinit.f (tlfield、tlcalc、tlspect) | RadiationField・RadiationSpectrum | radfield.py | 翻訳済み。RadiationField は gs が SIGSEGV (B49) なので等速運動の電荷の速度場と相対 6e-3 の確認のみ。RadiationSpectrum は gs と一致 (NPARA の fork 並列は逐次) | test_radfield.py、docs/45_radlight_guide.md |
| gamma.f、hg.f、tfunc*.f (特殊関数)、tffindroot.f、D.n | 特殊関数・Fit・FindRoot・NIntegrate | specfn.py、fitroot.py、miscfn.py | 458 式が gs と一致 (K64 が誤値を返す B43〜B47、B53、B54 は真値)。Fit・FindRoot は数値微分で 5e-10〜1e-9。NIntegrate は複素数の被積分関数を扱い (gs と虚部 1e-15)、`LogGamma1` もある。複素 FindRoot (gs が壊れた値: B61)・複素引数の Bessel・超幾何・PolyLog の \|z\|≥150 の非整数次は未対応 | docs/33 |
| ActLie, ActGRA, tdlat, tmap, polpar 他 | 空スタブ・描画・デッド | — | 対象外 | — |

## 9. 未翻訳・一部未対応・未照合の一覧 (2026-11-17)

最新の一覧は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) にあり、この節はその要約 (項目を直したら 2 つを一緒に更新する)。
2026-11-17 までに対応・照合した APERT 要素の追跡 (SOL 領域内を除く)、FFS 内の `MOMENTUM=`・`[NO]EXPAND`・RENUMBER・DISP の各モード・REFERENCE・DUMP・OrbitGeo、MakeBump・Couple の ElementValues (VISIT/BYE を除く)・Rational・NIntegrate の複素数・LogGamma1・ParametricPlot・ListContourPlot・OptimizeOptics・Entry/ListBox/Scale のメソッドは、ここでは未対応に入れない。

### 9.1 未翻訳 (理由つき)

| ファイル・機能 | 理由 |
|---|---|
| trackb.f、tracke.f、trcoda.f、ttinit.f、ttstat.f (MEASURE / TRCOD) | gs で TRCOD が SIGSEGV、MEASURE は出力なし (B42) で参照が取れず、使うスクリプトも見つからない。tlum・tluma だけ翻訳 (`measure.py`) |
| wjfit、wfres、wfbin、wiord、wjdraw | MAIN の TRACK の Fourier 解析で、FFS から使えない |
| tsole1 の領域内の MULT・BEND・CAVI (tmulte1)、SOL の DX DY DZ CHI・F1 (Twiss・Emittance[]・追跡)、領域の tdgeo | 領域内の BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げた (B69)。MULT・CAVI は未着手 |
| tspac.f の追跡側 (SPCH 要素) | K64 の tturn.f:521 `case (icSPCH)` が PSPAC フラグのとき `tpspac` を呼ぶ。pysad は SPCH を読むだけで物理は未翻訳 (PSPAC・WSPAC・SPAC フラグも効果なし)。PROF は K64 のソースに無い。GaussianCoulomb 関数は翻訳済み |
| LOSSMAP | gs は TrackParticles の戻りに 2 行 (落ちた要素と周) を足すが、pysad は未実装 |
| RADTAPER | フラグを保持するだけで、磁石の値が変わらない |
| System、Fork、Shared、`BeginPackage`/`Begin` の文脈 | 登録だけの限定版・見送り (docs/46 §1) |
| tdgeo.f の四極などの DX・DY の解析微分 (`qdtrans`・`tftmatu` が要る) | 数値微分で代用。優先度は低い |
| Polarization (平衡偏極度) の有限値 | K64 が NaN (B55)。NaN のまま写す |
| PHSROT | K64 で動作しない (B38)。追跡は NotImplementedError |
| 光子 (PHOTONS) の tphrec | 他の要素の光子。BEAMBEAM の分だけ対応 |
| `NOKEEPEXP`、RENUMBER の INDX 属性の初期値、`Out[n]:=` の表示 | 未対応 |
| 旧 WAKE コマンド (tfwake.f) | 移植不要 (K64 は WAKECOMMAND.n を呼ぶ。B66) |

### 9.2 一部未対応

| 機能 | 未対応の部分 |
|---|---|
| MULT ANGLE≠0 (multa / multae) | Twiss の BZ・RF・F1・SK0・FB1/FB2 (追跡は FB1/FB2 まで対応) |
| MAP | GEO・CompiledMap。TRACK は TrackParticles だけ (MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない) |
| POL | SOL 領域内の MULT・BEND・CAVI のスピン。スピン追跡は損失の処理なし |
| TCAVI | Emittance[] は L=0 のみ gs と一致 (L≠0 は gs が壊れる: B36)。RANKICK・RANPHASE 未対応 |
| UND | K64 は磁場を設定できない (B37) ので磁場なしのドリフトのみ |
| APERT | SOL 領域内は NotImplementedError。JDPX・JDPY の再入射は粒子 12 個以上で gs と乱数列が違う。RAD・FLUC つき、TRPT、MAIN の TRACK の損失報告、DAPERT での APERT は gs と未照合 |
| RF のある MULT (tmultiacc) | F1・FRINGE・AUTOPHI なし |
| ResponseMatrix | TRPT (加速) の経路 |
| tsteer | F1/FB1/FB2 付きは NotImplementedError |
| Fit・FindRoot・特殊関数 | 複素 FindRoot (gs も壊れる: B61)、複素引数の Bessel・超幾何、PolyLog の \|z\|≥150 の非整数次 |
| `ElementValues` (Couple.n) | VISIT/BYE の退避・復元、フィット変数との連動 (`CouplingMatrix`・`CoupledKeys`) |
| 描画・GUI | CanvasSymbol、Canvas3D*、`UseBG->True` の文字背景。OpticsPlot の text 項目数 (gs 69、pysad 74。Tk のフォントによる)、BeamPlot (gs は SIGXX=0 で楕円なし)、最初の `Window[]` の親 (gs は Toplevel)、`Packages/CSR.n` 全体は読み込まない (docs/37) |

### 9.3 未照合 (コードはあるが gs と数値比較していない、または gs が参照を出さない)

| 機能 | 理由 |
|---|---|
| FITP>1 の最終解、tffamsetup (fam 点) | gs (WSL) が GO で再起動する (D1、B05)。自己検証のみ |
| CSROYMatrix | gs が wrongtype を返す (B60) |
| tbrote の srot 更新、tbendebody の bsi2 (B56・B58) | ソースの読みだけで、gs の個別確認なし |
| RadiationField | gs が SIGSEGV (B49)。速度場との比較のみ |
| measure.py (tlum・tluma) | gs で MEASURE が動かない。解析値のみ |
| SynchroBetaEmittance | 12 セルのリングで emitx 4e-7〜1e-5。24 セルの平面リング + CAVI VOLT=2e6 (AzimuthalModes->7) の 4 点で emitx が gs と相対 1.4e-7〜1.4e-6 (emity は雑音)。結合のあるラティスは gs が SIGSEGV で落ちて比較できず (pysad は NaN の ValueError で失敗)。不安定なチューンは未確認 |
| tfloor.f の cfunc (xsincos 等) | 個別実装 |
| `Optics.n` の PeriodicR、CorrectOrbit の Micado・CorrectLocal、`dapert.n`・`element.n`・`emit.n`・`beamline.n`・`functions.n`・`Misc.n` のファイル単位の比較 | 記録なし |
| 未解決の差: BEND の F1=0.1 FRINGE=3 + RAD の追跡 (相対 0.7)、L≠0 の CAVI + DX・DY (1e-3)、SOL 内の四極 + RAD (現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い。B68、local/wake/verify_log.md。ソースの箇所は未特定)、K3 の強い MULT リングの qcod (Emittance が NaN)、EMIT のあとに USE すると gs は EMIT の効果が追跡に残り pysad は捨てる、RADINT と EMIT の Jx の 2.4% の差、`BUCKET HIGHT` の 0.010045 と 0.010035 | 原因の調査が未了 (docs/11 §E、`history/work_unimpl.md`、docs/46 §4) |

### 9.4 優先度

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | Monitor/Steer/CalculateOptics は動く。誤差入り 6D emit と DA の反復が 600 s を超える |
| 中 | SOL 領域内の MULT・BEND・CAVI、CALC6D・Touschek | 要素の完全対応 |
| 中 | 9.3 の未解決の差 | 原因の切り分け |
| 低 | tdgeo の四極の微分、wfres 系、tspac の追跡 | 必要になったら |
| 不要 | MON/MCOD/TCOD、Obsolete・非リンクのファイル、旧 WAKE | K64 で無効 |

## 10. `<SAD>\src` の *.f のうち pysad の中に名前が出てこないもの

`<SAD>\src\*.f` 474 本を `grep -w -i` で `pysad/*.py` と突き合わせた (2026-11-16)。名前が出ないだけで、別名のルーチンに直訳されている場合がある
(例: `tins`・`tinse` は `lattice.py` の INS の扱い、`talign`・`nalign` は ffscmd / lattice のずれ・回転)。
物理計算で名前が出ないもの: `qins tins tinse` (INS、§2。実装は別名で gs と一致)、`qmult qddrif` (qtwiss・qdcell 内の別実装に吸収: docs/11 D5)、`tfmap tmap` (MAP は別名)、
`trackb tflifetrack twinit wfbin wfres wiord wjdraw` (§1・§8・上の 9.1)、`tbrad tquads tquase tsdrad tserad tstrad tqrad trad trade` (Obsolete)、
`tconv` (追跡の変換は track6 内)、`tmuld tmultr` (行列の積は numpy)、`temitf` (EMIT 報告は ffscmd。ファイル名としては出ない)。
それ以外は、旧 MAIN のユーティリティ (`pmbump pwrite pundo pkill` など)、SADScript のインタプリタ内部 (`tfeval1 tfstk tfpart tfsort tfreplace itf*` など)、描画 (`geodrw phdrw pmbdrw tdplt twsdrw`)、ファイル・バッファ (`filaux filbuf flmgr getbuf`) で、
pysad の方針 (SADScript は Python のインタプリタ、ファイル・描画は Python のライブラリ。`15_sadhelp_purpose_policy.md`) により翻訳しない。
一覧の再現: `for f in <SAD>/src/*.f; do b=$(basename $f .f); grep -q -w -i -- "$b" pysad/*.py || echo $b; done`。
