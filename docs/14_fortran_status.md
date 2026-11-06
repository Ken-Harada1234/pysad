# 14. SAD Fortran ファイル別 pysad 翻訳状況

最終更新: 2026-11-06 (現状のスナップショット)

SAD は K64 版、ソースは `<SAD>\src`。
「デッドコード・ビルド対象外」の判定は mk/sad.obj.mk の有効リストに無い、または call 0 件のもの。行数の多いデッドコードは節ごとに 1 行にまとめる。
根拠: `01_overview.md`, `10_benchmarks.md`, `history/00_RESUME.md`、`work_*.md`、`pysad/*.py` の実コード。

## 0. 凡例と現在の pysad モジュール

| 状態 | 意味 |
|---|---|
| 翻訳済み・照合済み (桁) | pysad に直訳があり、SAD の出力と比較済み。桁は最大誤差 (`10_benchmarks.md` の値) |
| 翻訳済み・未照合 | コードはあるが SAD との数値比較が未実施、または再現が不完全 (括弧内に理由) |
| 未翻訳 | pysad にコードがない (NotImplementedError で止まるものを含む) |
| 対象外 | デッドコード、非リンク、GUI/グラフィックス、EPICS など |

現行の対応モジュール: parser / sadelem / lattice / elkern / perfkern / turnkern / emit6d / teigen / match6 / match6d / tsvd / track6 / track6n / track6n_ext / dapert / sadrandom / sadscript / sadlist / ffscmd / calc6d / touschek / geo / emit_extra / radint / maintrack / trptacc / optics_fn。
`ffs.py`, `matching.py`, `tracking.py`, `optics.py` は非推奨の旧 API (SAD 未照合の独自実装) で、表には数えない。
calc6d / match6d / touschek / geo / track6n_ext / emit_extra はいずれも `ffscmd.py` から呼ばれる (CALC6D、解析微分、INTRA/TouschekLifetime、LINE["GX"..]/GEO/ORG、MULT/SOL の numba 追跡、Emittance[] のオプション)。

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
| tapert.f | 損失判定 | track6n.py `k_apert`、dapert.py | 翻訳済み・照合済み | DA の地図一致 (L16A)。ただし損失周の乱数消費は不一致 (RAD+FLUC) |
| tsteer.f | BEND ANGLE=0, K0≠0, L≠0 | track6.py `_tsteer_ens` (F1/FB1/FB2 フリンジ付きは NotImplementedError) | 翻訳済み・照合済み | bendk1_*/pf2 bump 系 (history/work_bend_k0.md, history/work_trpt_acc.md 10 節) |
| tmulta.f, tmultr.f | MULT ANGLE≠0 | — (NotImplementedError) | 未翻訳 | — |
| ttcav.f | TCAVI (クラブ空洞) 追跡 | elements.py にクラスのみ | 未翻訳 | — |
| tcoord.f | COORD 要素 | parser が受理のみ | 未翻訳 | — |
| temap.f | MAP 要素 | — | 未翻訳 | — |
| undulator.f | UND | — | 未翻訳 | — |
| bbstrhl.f, bbtool.f | BEAMBEAM | parser が受理のみ | 未翻訳 | — |
| tbbbrem.f, csrtrack.f, csroy.f, csrtest.f | BBBrem / CSR 関数 | — | 対象外 (衝突型・FEL 用) | — |
| track.f, tracka.f, ActTra.f, ttinit.f | MAIN の TRACK コマンド・初期分布 | ffscmd `TrackParticles[]` (SADScript 関数経由)、MAIN の TRACK コマンドは maintrack.py | 翻訳済み・照合済み | TrackParticles は L6〜12。MAIN TRACK は SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 (SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV)。PHSPLOTS の描画は無視、FOURIE ON は未対応 |
| trackb.f, tracke.f, trcoda.f | MEASURE / TRCOD | ffscmd が語を受理するのみ | 未翻訳 | — |
| trackd.f / tfdapert.f | DynamicAperture 本体 | dapert.py | 翻訳済み・照合済み | 地図 153/153、Score 一致 (L16A、NORAD) |
| tquads, tquase, tbrad, tmulbs, ttdr, bb, hg 等 | Obsolete・非リンク | — | 対象外 | — |

## 2. 線形光学

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| qtwiss.f (qtwiss1, qmat2twiss) | Twiss 伝搬 | lattice.py `twiss_table` | 翻訳済み・照合済み | ≤1e-10 (L1,12,13)、結合・誤差経由 ≤3e-7 (L3, L5)。結合時の qmat2twiss は部分 |
| qtwiss.f qcod / tcod.f | 4D 閉軌道 | lattice.py `_qcod_periodic`, `find_cod` | 翻訳済み・照合済み | 9〜10 桁 (L3B) |
| qcell.f (qcell1) | CELL 周期解・INS | lattice.py `twiss`, `ins_twiss` | 翻訳済み・照合済み | CELL ≤1e-10。INS: 軌道 ≤1.3e-11、Twiss ≤2.6e-7 (L5)。初期軌道オフセット付き INS はチューン不一致 (amp_dep_ts) |
| qmdiag.f | 4x4 対角化 | lattice.py | 翻訳済み・照合済み | ≤1e-10 |
| tquade.f, tbende.f, tdrife.f, tmulte.f | 6D 転送行列 (emit 用) | emit6d.py | 翻訳済み・照合済み | Emittance[] ≤1.3e-11 (L16B〜E) |
| tthine.f | 薄レンズ行列 | emit6d.py | 翻訳済み・照合済み | ≤1.3e-11 |
| tcave.f | 空洞行列 | emit6d.py | 翻訳済み・照合済み | ≤1.1e-7 (L4、RADCOD ON) |
| tsole.f, tsolque.f | ソレノイド行列 | emit6d.py に SOL なし | 未翻訳 | — |
| tchge.f, qchg.f | ずれ・回転座標変換 | lattice.py, emit6d.py、match6d.py `qchg` | 翻訳済み・照合済み | ≤1e-8 (L3B, L3C) |
| tturne.f | emit 用 1 周ディスパッチ | emit6d.py `_turn` | 翻訳済み・照合済み | ≤1.3e-11 |
| qdcell.f, qdquad/qdbend/qdmdia/qdtwis.f, tffsqu | 解析微分 | match6d.py | 翻訳済み・照合済み | esrf/matching の反復列・K1 を 12 桁再現 (L22A)。数値微分とは 2e-8〜4e-7 |
| qcell6d / temit plot 経路 | CALC6D | calc6d.py | 翻訳済み・照合済み | 1.7e-10 (L18A/B/C) |
| ttcave.f, tcoorde.f, qins.f の INS 行列 | TCAVI / COORD / INS 要素の行列 | — | 未翻訳 | — |
| qquad, qthin, trad, twelm, pfcoup 等 | 非リンク・デッド | — | 対象外 | — |

## 3. エミッタンス

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| temit.f (tecalc, synchm) | 6D 平衡ビーム行列 | emit6d.py `emittance6d` | 翻訳済み・照合済み | 全量 ≤1.3e-11 (L16)、RADCOD ON 1e-7 (L4) |
| INTGRL.f | RADINT (I1〜I5, 結合, 色収差) | radint.py (FFS RADINT)、lattice `_radiation_integrals_bend` (BEND のみの別実装) | 翻訳済み・照合済み | 7 桁 (L1、L29A) |
| temit.f の INTRA, tintraconv | IBS 付き emit / Touschek | touschek.py `touschek_table`, `touschek_lifetime` | 翻訳済み・照合済み | 15 ケース ≤7e-12 (L17A)、RADCOD ON 1.7e-8 (L17B) |
| tfemit.f, emit.n | Emittance[] のオプション (OneTurnInformation, Matrix) | emit_extra.py | 翻訳済み・照合済み | 2.1e-10 (L21A)。偏極・SpinTune は 0 固定 |
| temitf.f, temits.f | EMIT / SYNCHROB コマンド | ffscmd (EMIT 報告、SYNCHROB 語)、emit_extra | 翻訳済み・未照合 | — |
| calpol, spnorm (偏極) | POL | — | 未翻訳 | — |

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
| tcorr.f (DELCOR), tvcorr.f, tcftr.f | 相関付き誤差・FFT | — | 未翻訳 | — |
| trcoda.f (TRCOD) | 追跡による COD | ffscmd が受理のみ | 未翻訳 | — |
| terror.f (DELX/DK/DTHETA…) | 誤差コマンド | ffscmd `_error_cmd` | 翻訳済み・照合済み | L3C/3D、L7-2 |
| MON/MCOD/TCOD コマンド | COD 補正 | — | 対象外 (K64 ビルドでは何もしない) | — |
| tsvdm.f 等 | SVD | match6 の tfsolv が numpy で代替 | 翻訳済み・未照合 (numpy 流) | マッチング K1 ≤1e-10 |

## 6. FFS・マッチング

| Fortran | 役割 | pysad | 状態 | 照合 (Level) |
|---|---|---|---|---|
| tffsa.f | コマンドループ | ffscmd.py | 翻訳済み・照合済み | L13, L14 (実スクリプト) |
| tffscalc.f, tffsmatch.f | 残差・GO (Newton/DESCEND, fuzz, fmincube) | match6.py | 翻訳済み・照合済み | 反復数・残差・K1 ≤1e-10 (L13A/B)。FITP>1 は参照値なし (L15A) |
| tffsqu / qdcell | 解析微分 | match6d.py | §2 参照 | L22A |
| tgetfv.f, tfgetv.f, tflag.f, tfsetv.f | FIT 条件・要素値・フラグ | ffscmd.py | 翻訳済み・照合済み | L13, L14 |
| tftwiss.f, tfchro.f, tfoptics.f | Twiss[]/LINE[]/Element[]、色収差、CalculateOptics | ffscmd.py, lattice.py | Twiss[]/LINE[]: 照合済み (L14A、PolynomialFit 1e-7)。CHROMATICITY (tfchro)、CalculateOptics (tfoptics; 4 次元の cell/INS・区間・結合・分散、idim=3 は全周): 翻訳済み | — |
| tfgeo.f, tdgeo.f, tgrot.f | 幾何 (GX…、ORG) | geo.py (LINE["GX"] / fit GX.. / DISP G 経由) | 翻訳済み・照合済み | 4.8e-13 (L19A/B/C)。tdgeo (fit の GX.. の解析微分) は未翻訳で、幾何の条件は数値微分 |
| tfdapert.f | DynamicAperture[] | dapert.py | 翻訳済み・照合済み | L16A |
| tffamsetup.f | FITP の有限振幅ファミリー (MatchingAmplitude/NFAMP) | match6.py `fam_setup` ほか | 翻訳済み・未照合 (SAD が FITP>1 の GO で再起動するため自己検証のみ) | test_fam_points.py |
| tffswake.f | ウェイク | — | 未翻訳 | — |
| tfstk, tfeval, tfetok 等 (SADScript) | インタプリタ | sadscript.py、sadlist.py (リスト系組み込み) | 翻訳済み・照合済み (実スクリプト) | L14A (esrf calc_chromaticity) |

## 7. ビームライン定義・座標

| Fortran | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| toplvl / doelem / doline / expln | MAIN 部 (LINE, 繰り返し, 再定義) | parser.py | 翻訳済み・照合済み | L12 (3563 要素) |
| talign.f, nalign.f | ずれ・回転の設定 | ffscmd, lattice | 翻訳済み・照合済み | L3B/3C |
| tfgeo.f | フロアー座標 | geo.py | §6 参照 | 4.8e-13 |
| tfloor.f | 数学関数 (cfunc) | 各所に個別実装 (xsincos 等) | 翻訳済み・未照合 (個別) | — |
| tsgeo.f, phsrot.f | ソレノイド幾何・位相回転 | — | 未翻訳 | — |
| tpara, tscale 等 | 空・デッド | — | 対象外 | — |

## 8. 寿命・その他

| Fortran | 役割 | pysad | 状態 | 照合 |
|---|---|---|---|---|
| Touschek (tintra 系、emit.n) | TouschekTable/Lifetime | touschek.py | §3 参照 | ≤7e-12 |
| tflifetrack.f | LifeTrack | — (SAD 側も空スタブ) | 対象外 | — |
| tlum.f, tluma.f | ルミノシティ | — | 未翻訳 | — |
| txwake.f, twinit 他 | ウェイク | — | 未翻訳 | — |
| wfbin.f, wfres.f | 周波数解析・共鳴 | — | 未翻訳 | — |
| tspac.f (GaussianCoulomb) | 空間電荷関数 | — | 未翻訳 | — |
| ActLie, ActGRA, tdlat, tmap, polpar 他 | 空スタブ・描画・デッド | — | 対象外 | — |

## 9. 未翻訳ファイルの優先度

使用頻度は `12_user_scripts_status.md` (845 スクリプト) と `history/work_pf2_survey.md` (pf2 31 本) による。

| 優先 | ファイル | 理由 |
|---|---|---|
| 中 | correct_cod 系の全実現の完走 (性能) | Monitor/Steer/CalculateOptics/CodCorr は動く。誤差入り 6D emit と DAsv の反復が 600 s を超える |
| 中 | tdgeo.f (fit の GX/GY 微分) | pf2 matching_geo_ring |
| 中 | tcorr.f, tvcorr.f, tcftr.f (DELCOR) | 誤差・補正の実機スクリプト |
| 中 | tmulta.f (MULT ANGLE≠0)、tsole.f/tsolque.f (emit の SOL) | 要素の完全対応。実スクリプトでの頻度は低い |
| 低 | ttcav/ttcave, tcoord/tcoorde, temap, qins | tturn/tturne から有効に呼ばれるが esrf/pf2 では未使用 |
| 低 | tlum/tluma, txwake, wfres, tspac, phsrot, tsgeo | 衝突型・ウェイク・解析用 |
| 低 | undulator.f, bbstrhl.f | PF-II の ID は pf2 で MULT/BEND 扱い (要確認) |
| 不要 | CSR 系 (csroy 他)、BBBrem、MON/MCOD/TCOD、Obsolete・非リンクのファイル | FEL・衝突型専用または K64 で無効 |
