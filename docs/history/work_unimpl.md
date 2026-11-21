# 未実装要素の実装 (2026-11-14 以降)

最終更新: 2026-11-16 (整合: 同日、重複・古い記述に「(11/16 時点)」の注を足した。内容は消していない。夜: 未確認の項目を確かめ、解決したものと新しい差を下の要約に反映)

## 11/16 時点の未解決の差・未対応の一覧 (以下の節の要約。詳細は各節)

未解決の差 (pysad と gs が合わない。原因は確定していない):
1. BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると相対 0.7 (「POL の TrackParticles」節)。
2. L≠0 の CAVI に DX・DY で相対 1e-3 (同上)。
3. SOL (BZ≠0) 内の四極 + RAD の追跡で、現在の gs (10-19 ビルド) が保存済みの参照 (bench_level9、10-14) と別の値を返す (同上)。**11/16 夜に絞り込み: gs は RAD ON かつ BZ≠0 の領域内で四極の力を入れていない疑いが強い (四極を外した場合と 8 桁一致、BZ=0.01 と 0.2 でほぼ同じ値、2 つのビルドで 10 桁一致)。pysad は物理と保存済みの参照に合う。K64 ソースのどこかは未特定** (`work_session_20261116.md` §11.3)。
3b. Tkinter の Widget 値 (11/16 夜に発見): gs では `e=Entry[w]`・`w=Window[]` が `Widget[…]` を値に持ちメソッドが働くが、pysad は裸のシンボルのまま (ウィジェット自体は作られ、親が Root になる)。`Tkinter.n` の `w:=Widget[w,wid[master]]` が働く条件が未特定で未修正 (`work_gui_tclbridge.md`)。
4. K3 が強くずれの大きい MULT リングで qcod の閉軌道が gs と違う解になり Emittance が NaN (「POL の MULT」項、`local/wake/pm1.sad`)。
5. `emit6d._tquade` (Python 版) の δ 微分の欠落 (1e-5)。「後でまとめて直す一覧」1。
6. SOL 領域内の BEND は gs が偏向をほぼ無視して合わず取り下げ (「RADLIGHT」節の末尾)。

解決したもの (11/16 夜。`work_session_20261116.md` §11.2): B33 の出力座標の差 (原因はセッション最初の l_track=0 で最初の tapert が粒子を落とさない K64 の癖。`TrackRun.l_track0` で再現、gs と全桁一致)、UNIPOL・RELW (FIT にフラグがつながっていなかった。gs と全桁一致)、DUMP の DK・DDK の書式、MATRIX コマンド (名前だけの位置・ラティスに無い名前)、LDSolve の複素行列 (DiagonalMatrix・Dot・Inverse)、OpenRead の失敗の gs 互換。
新しく見つけた限界 (未修正): 結合リングの SynchroBetaEmittance (gs が SIGSEGV、pysad は NaN の ValueError)、RadiationField の複数区間 (区間の始めの点が重なり `x must be strictly increasing`。gs は SIGSEGV)、SPCH・LOSSMAP は未翻訳 (SPCH は `tturn.f:521` が PSPAC のとき `tpspac`)、`Options[g]` は pysad が {} を返す (gs は全オプション)。
確認済み (差なし): 11/16 夜の項目は `work_session_20261116.md` §11.1。QUICK が SKIP する 4 ベンチも実行して 50 ベンチすべて FAIL 0、試験は 520 通過・3 skip (523 件)。

未対応: SOL 領域内の MULT・BEND・CAVI、CALC6D・Touschek と SOL の組み合わせ、ExternalMap["GEO"]・CompiledMap、DELCOR の CAVI・OCT、PolyLog の |z|≥150 の非整数次、複素引数の Bessel・超幾何、ResponseMatrix の TRPT、BEAMBEAM 以外の要素の光子 (tphrec)、MessageName の関数形 (見送り)、BeamPlot・CanvasSymbol・Canvas3D*。
未翻訳 (理由つき): MEASURE/TRCOD (B42)、`wjfit/wfres/wfbin/wiord/wjdraw`、`tspac.f` の追跡側、tdgeo の四極などの変数、`Polarization` (B55)、PHSROT (B38)、旧式の WAKE (移植不要)。
未確認: SynchroBetaEmittance の実用のリング・不安定なチューン、CSR の長い流れ (ScanTrack・OYStability)、BBBrem の `Parallelize`、GUI の見た目と操作。

優先順位 (利用の多さと実装の小ささ。11/14 に立てた順。11/16 時点で 1〜5 と 7 の MAP、DELCOR は済み、6 PHSROT は K64 で動作せず): 1 COORD、2 TCAVI、3 UND、4 MULT ANGLE≠0 (済み)、5 SOL の emit・Twiss 行列 (Twiss は済み、emit の放射が残り)、6 PHSROT、7 MAP 以降。

## 済み

| 要素 | 結果 | 文書 |
| --- | --- | --- |
| COORD | 追跡・Twiss・幾何・Emittance[] が gs と一致 (`pysad/coord.py`) | docs/25_coord.md |
| TCAVI | 追跡 (RAD FLUC まで) が gs と一致。Emittance[] は L=0 のみ gs と一致、L≠0 は gs が壊れる (B36) | docs/26_tcav.md |
| UND | 追跡・Twiss・Emittance[] が gs と一致。K64 は磁場を設定できず、磁場なしのドリフト (B37) | docs/27_und.md |
| MULT ANGLE≠0 | 追跡 (11 通り、RAD FLUC 含む) は相対 0、Twiss 1e-13、Emittance[] 1e-11 で gs と一致 (`pysad/multa.py` `multae.py`)。BZ・RF・F1・SK0・FB1/FB2 の Twiss は未対応 | docs/28_mult_angle.md |
| DELCOR | tcorr/tvcorr が gs と一致 (チューン 1e-14、軌道 1e-10、乱数列も同じ)。MARK は DETR・DPX に入る (K64 の添字のずれ)。CAVI は未対応 (`ffscmd._delcor`) | docs/31_delcor.md |
| MAP | ExternalMap の TRACK・EMIT・OPTICS が gs と一致 (追跡 10 桁、Emittance[] 1e-11、CALC 1e-11)。GEO・CompiledMap は未対応 (`Map`、`ffscmd._map_*`) | docs/30_map.md |
| SOL の Twiss・Emittance[] | CALC / INS の Twiss が gs と一致 (チューン 1e-15、領域内の β・α・分散・結合 1e-11)。Emittance[]・EMIT は ε_x 2e-12〜4e-9、ε_y 4e-10〜2e-8、チューン 1e-11、軌道 1e-9 で一致 (`pysad/solmap.py` `solemit.py`)。CALC6D・Touschek・幾何 (tsgeo) は未対応 | docs/29_sol_optics.md |
| 特殊関数・数学関数 | 458 式が gs と一致 (specfn.py、基本関数の複素対応)。K64 が誤値を返す式 (B43〜B47) は真値 | docs/33_specfn_fit_misc.md |
| Fit / FindRoot / NIntegrate / FitGaussian / FitEmit | 数値微分で gs と 5e-10 (Fit は線形で 1e-9)。ベクトル独立変数・複素 FindRoot は未対応 (`fitroot.py`) | docs/33_specfn_fit_misc.md |
| DownhillSimplex ほか | optimize.n の Amoeba を移植 (gs と 1e-13)。StringToStream / TemporaryName / ListRandom / FixedPointList / SelectCases (`miscfn.py`、`sadlist.py`) | docs/33_specfn_fit_misc.md |
| 軌道補正・最適化パッケージ | CorrectOrbit / MakeBump / Monitor / Steer / OptimizeOptics など (correction.n、bump.n、corutility.n、optimize.n、rational.n、tableform.n) を SADScript のまま実行し gs と一致。ResponseMatrix・LinearSolveConditioned・Master・FFS["VAR"] を翻訳 (`corr.py`、`pkgload.py`) | docs/34_sad_packages.md |
| PHSROT | K64 で動作しない (行列がすべて 0、PX PY PZ が指定できない: B38)。追跡は NotImplementedError | docs/21 B38 |

## 見つけたこと

- 追跡の作り方: tmulta は tbend (k_bend0 / k_bendi) をスライスごとに呼ぶので、`k_bend0` に `ini` (放射の基準 pxr0 を設定するか) を足して再利用した。Twiss・Emittance[] は、スライスのパラメータをもった BEND を誤差なしの核 (`_tbende_core`、`_bend_core_cod_matrix`) に渡す。
- parser は MULT の `CHI2` `TAPER` を読んでいなかった (追加)。MULT の回転キーワードは `ROTATE` (`ROTATION` は gs で Undefined keyword)。
- ANGLE≠0 の MULT の CALC・Emittance[] は pysad が ANGLE を無視して黙って違う値 (水平だけ) を出していた (`lattice._el_has_offset` が設計軌道の近道に入れていた)。実装で解消。

- parser は COORD / TCAVI / UND / PHSROT / BEAMBEAM を黙ってドリフトにしていた。専用の要素クラスにし、`DX` などは誤差の属性にせず `params` に置いた。
- **SOL (BZ≠0) を含む CALC・EMIT・Emittance[] は pysad が SOL を無視して黙って違う値 (結合が 0) を出していた**。gs は結合を持つ。実装の前に、`SADError` で拒否するようにした (`ffscmd._check_sol_optics`)。
- TRPT の放射の分割 (ndivrad ×10) は BEAMBEAM の作業で直した (docs/24 §4)。

- SOL の Emittance[] (solemit): K64 の `dbeam(8)` が `beam(10)` を使う誤記らしい箇所 (B34 から数えて次の B39 とする) をそのまま写した。SOL に GEO=1 がないと K64 は S 座標が途切れ、MomentumCompaction が周長の誤りで 100 倍になる (pysad は GEO を見ないので影響を受けない)。
- `emit_extra.emittance_rules` (Emittance[] の Orbit・Matrix 規則) と `calc6d._turn` は Python 版の `_tquade` を使い、軌道のある QUAD で gs と 1e-5 (ε_y 1e-4) ずれる。SOL なしでも出る既存の差で、`emit6d.tturne` (elkern) は 1e-12 で合う。未調査。
- SOL つきの Emittance[] は emit_extra に入れず `run_emit` の値に落ちる (Stable=0、HarmonicNumber=NaN、Matrix 系の規則なし)。上の差を直せば emit_extra にもフックを入れられる。

## 次 (11/14 夜の記録。11/16 時点で、MEASURE/TRCOD 以外は下記の節のとおり済みまたは判断済み: tlum は `measure.py`、GaussianCoulomb は docs/32、CSR は docs/36、BBBrem は docs/32、tdgeo は DRIFT・BEND だけ解析微分)

- SOL の残り: CALC6D・Touschek、領域内の MULT・BEND・CAVI (tmulte1)、tsgeo/tdgeo の幾何 (GEO=1 のとき S 座標を領域でつなぐ)。Twiss と Emittance[] は済み。参考: gs は領域ごとに `tsole1` (tdrife0・tquade・tmulte1 に bz を渡す、tsconv、tsfrie) で、Twiss の各要素の出力も別の座標系。pysad は track6n_ext で追跡だけある。
- MEASURE / TRCOD (trackb / tracke / trcoda)、tlum、GaussianCoulomb (tspac)、CSR、BBBrem。


## 後でまとめて直す一覧 (自動で実装を進めた間に見つけたもの。直さずに記録だけした)

### pysad 側の不具合・改善点

1. **`emit6d._tquade` (Python 版) の解析ヤコビアンが、x がずれた軌道で δ の微分を落としている。**
   `emit_extra._turn_record` (Emittance[] の Orbit・Matrix・BeamMatrices の規則)・`calc6d._turn`・`touschek` が直接 `_tquade` を呼ぶので、
   DX/DY のある QUAD を含むリングで T の δ 列が 1e-5 ずれる (T[0,5] で 9.4e-6、T[4,0..3] で 1e-6〜9e-6)。ε_y で 1e-4、ε_x で 5e-6 相対。
   `emit6d.tturne` (elkern) は gs と 1e-12。`E._USE_ANALYTIC_JAC = False` にすると Python 版でも elkern と 5e-13 で合う (確認済み: `local/sol/h5.py`)。
   原因は `_jac_tqente` / `_jac_tsolque_body` の δ 依存 (集束が 1/(1+δ) で変わる項) が x ≠ 0 のとき抜けていること。
   直し方の候補: `_tquade` の軌道つきの経路を数値微分にする、または解析式に δ 項を足す。
2. `emit_extra.emittance_rules` が SOL を含むと使えず、`ffscmd.emittance` が `run_emit` の値に落ちる (Stable=0、HarmonicNumber=NaN、Orbit・Matrix 規則なし)。
   1 を直せば `solemit` を emit_extra / calc6d に入れて解消できる (フックの一行は docs/29 に書いた)。
3. SOL 領域の CALC (数値微分) は gs の解析行列と 1e-10 程度の差が残る。`tsolque` の解析行列を写せば消える。
4. USE 直後の `Matcher._calc_optics` で SOL 領域を annotate していなかった (直した: `_matcher_optics`)。同種の入口が他にもあるかもしれない (`calc` を経由しない光学計算)。
5. pysad は GEO がない SOL 領域で K64 のように S 座標を途切れさせない (警告だけ。K64 は周長が狂う)。tsgeo は翻訳した (`sgeo.py`、docs/29) が、領域内の BEND・MULT・CAVI は未対応で、
   ずれ・傾きのある SOL の追跡・Twiss・Emittance[] (tsconv の CHI1..3) も未対応。
6. MAP: ExternalMap["GEO"] と CompiledMap は未対応。TRACK は TrackParticles だけ (MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない)。
7. DELCOR: CAVI は未対応。OCT は確かめていない。
8. GaussianCoulombFitted (twspfu) は実装済み (docs/32)。BBBrem パッケージの `Parallelize` は動作未確認。
9. (11/16: 担当を分けて全文書を更新中。html は `make_html.py` で作り直す) ドキュメント: docs/html/*、docs/17_implementation_summary.md は (11/14 時点で) 新しいモジュール (solmap, solemit, multa, multae, tcav, und, coord, gausscoul, bbbrem) を反映していない。

### K64 側の不具合・疑い (docs/21 に B39, B40 を記録済み。以下を追加)

- B39 (疑い): photons.f tradke の dbeam(8) が beam(10) を使う。
- B40 (確定): tcorr.f が MARK の DETR・DPX に書く。CAVI は書く欄と読む欄が食い違う。
- B41 (確定、読解): ttinit.f (MEASURE の初期分布) の `dxj=rlist(kx+26)` … `dpyj=rlist(kx+29)` は ky_ZPX, ky_ZY, ky_ZPY, ky_DP (添字 26〜29) で、
  意図した ky_JDX..JDPY (34〜37) ではない。dpyj には DP (運動量幅) が入るので、初期軌道の py が DP·tgauss だけランダムにずれる (tgauss も 4 回余分に消費する)。
- B42 (確定): この環境の gs (K64 ビルド) で `TRCOD` は SIGSEGV で落ちる (CELL; CALC のあと `TRCOD;`)。`MEASURE` は (INS; CALC; あとでも) 何も出力しない。
  ともに gs と照合できないので trackb / tracke / trcoda / ttinit / ttstat は未翻訳のまま。
- K64 の `MARK` の DX を DELCOR で書く意図の添字ずれ (B40) と同じ型の誤りが他にもありうる (ky_ 添字を直接書いた箇所)。

### 未翻訳のまま残した物理計算 (規模順)

- POL / calpol: 実装済み (docs/36。Emittance[]・EMIT・TrackParticles の偏極時間・スピンチューン・偏極ベクトル・スピンは gs と一致。MULT・空洞・SOL (ドリフトのある領域) まで済み)。残り: SOL 領域内の MULT・BEND・CAVI (領域内の要素そのものが未対応)、`Polarization` (K64 が NaN: B55)。
- CSR: CSRTrack・CSRHaissin・CSRSetupOY・CSROYMatrix・CSRMatrix (mode 1000 を含む)・CSRConvert (K64 の tfcsrvec と同じ)・CSRInit は実装済み (docs/36)。`Packages/CSR.n` の MakeZL・SkinDepth・PathDifference・StandardKList・VariableMesh・ZThre などは gs と一致。ScanTrack・OYStability などの長い流れは未照合。
- tdgeo (fit の GX..): DRIFT・BEND の変数は解析微分、四極の DX・DY などは中心差分。tlum / tluma は `measure.py` に直訳 (MEASURE が gs で動かないので、解析値との比較だけ)。wfres / wfbin / wjfit は MAIN の TRACK (tracka.f) の Fourier 解析で、FFS からは使えない。旧式の WAKE は K64 では WAKECOMMAND.n (配布に無い) を呼ぶので移植不要 (docs/23)。
- MEASURE / TRCOD (tracke・ttinit・ttstat・trcoda・trackb): gs で参照が取れない (B42) うえ、使うスクリプトも見つからなかったので翻訳していない。
- 描画: BeamPlot の描画の gs との比較、CanvasSymbol、Canvas3D*、`UseBG->True` の文字背景 (docs/37。OpticsPlot・LatticePlot は枠・矩形・線が gs と一致、文字の数はフォント依存)。
- RadiationField・RadiationSpectrum と、その入力の RADLIGHT の軌道は実装済み (gs と一致: test_radlight.py)。

10. (済) 記号微分 `D[f,x]` は D.n を SADScript のまま実行して gs と一致 (59 式)。記号式の標準形 (`0*x` や `x^0` を簡約しない) は gs と少し違う (docs/34)。
11. 複素 FindRoot は gs が壊れた値を返す (B61)。FitGaussian の Plot オプションは実装済み (gs と同じ項目数)。Fit のベクトル独立変数は実装 (gs が誤値: B52)。
12. ResponseMatrix の TRPT (加速) 経路。LDSolve.n・Optics.n は読み込みだけ確認。
13. 複素引数の Bessel・超幾何、PolyLog の |z|>=150 の非整数次 (PolyLog の z>1・z<-1 の非整数次と HurwitzZeta の q<=0 は実装済み: docs/33)。
14. (`Skip` と `Read` の型 (Character・Real・Word の読み方、数でない語の Real は式の値) は gs に合わせた。test_evalfixes.py)。`MessageName[f, "usage"]` の関数形 (gs は "Undefined MessageName" を出す。見送り)、(`Short`・`Inequality`・`Fork` は gs に合わせた。`ToInputString` は一致)。
    (11/16 時点: GaussianCoulombFitted は docs/32、RadiationField・RadiationSpectrum と RADLIGHT の軌道記録は「RADLIGHT」節のとおり実装済み。この項の以前の記述は「未対応」だった。)
15. 前回までの pysad の誤り: スクリプトの `ElectronMass` が 510999.06 だった (K64 は 510998.95)。数の書式 (`1e6`、`1e+10`)、`Part`、`Head`、`VectorQ` の差は修正済み (docs/34)。

## SynchroBetaEmittance (2026-11-16 に移植した)

`pysad/synchrobeta.py` に K64 の `tfemits`・`temits`・`temits1`・`tefsetup`・`tesolvd`・`tevdif`・`tevdifj`・`tevdifj1`・`tespl`・`tesumb`・`teintp1`・`teintm`・`teintb` を直訳した
(`SynchroBetaEmittance1[nus0, nus1, dnus, mphi]` は `ffscmd.fn_synchrobeta`、`Packages/emit.n` の `SynchroBetaEmittance` が呼ぶ)。
運動量のずれ (2·ndims+1 点) ごとに `_tcod_radcod` (閉軌道) と `tturne` (T・D・B) を取り、方位角モード展開の連立方程式を解く。
12 セルのリング + CAVI (VOLT=2e6 HARM=64) で gs と比べた: nus=0.01、AzimuthalModes=9 は emitx が相対 6e-6、conv が 2e-4。
nus=0.005〜0.02、AzimuthalModes=5 の 4 点は emitx が相対 4e-7〜1e-5 (共鳴の構造は同じ)。emity は 1e-18 の雑音で符号も違う。差は閉軌道の収束の履歴と、
条件数 3e6 の連立方程式の解き方 (K64 の tsvd と numpy の SVD) から来ると考えている (solve と SVD の違いは 4e-10)。試験は `test_synchrobeta.py` (許容 3e-5)。
未確認: 実用の大きなリング、不安定なチューン、`emiout`、結合のあるラティス、RFSW の組み合わせ。

## INS 要素 (2026-11-16 に確認: 追加の実装は不要)

- INS を含む 3 セルのリング (BEND E1・E2 付き) で、`Emittance[]` のチューン (相対 1e-15)、`TrackParticles` (2 粒子・3 要素と 10 要素、相対 1e-13) が gs と一致する。
  RFSW・RAD を入れても同じ (`CAVI` VOLT=2e6 HARM=64、`Emittance[]` のチューン 1e-15、追跡 1e-13)。`qins`・`tins`・`tinse` が担う部分は pysad の既存の実装 (`INS` の領域の扱い) で gs と一致する。
- `TrackParticles` の呼び出し形式は `TrackParticles[{loc, {{x..},{px..},{y..},{py..},{z..},{dp..},{flg..}}}, 終点, 始点, 終点]`。

## 幾何の FIT (tdgeo) の状況 (2026-11-16)

- 幾何の条件 (`fit gx 0; fit gy 0` など) の 0 は **リングを閉じるための目標値**で、
  BT (輸送路) やセルでは目標を任意の値にする (`fit gx <値>`)。pysad の幾何の条件 (`match6.py` の GX GY GZ CHI1..3) は任意の目標値を受ける。
- `type` でドリフトを並べる表示は、pysad が 1 行ずつだったのを、gs の `twbuf` と同じ詰め方 (DRIFT だけ 1 行に並べる) に直した。
- 幾何の目標が実際に動く例は gs で作れていない: 小さなリングで `fit gx 1.0; fit gy 2.0; free (ドリフト X) L; free (ドリフト Y) L; go;` は gs でも pysad でも変数が動かない (残差 1e9〜1e12、同じ値)。
  不安定なリングでは gs が不安定の罰則を残差に入れず fit が動くが、pysad は入れるので結果が違う (不安定なリングの fit は比べない)。
  動く例が得られれば、`tdgeo` の微分の差が出るかを比べる (得られた例は gs と pysad で収束位置が 1e-6 の程度の差。差は fit の収束判定と微分の違い)。
- tdgeo.f (147 行) の中身: 幾何の量 (GX GY GZ CHI1..3) の変数 (DRIFT・BEND の L、BEND の ANGLE、四極・6 極などの DX・DY) に対する解析微分。ドリフトの L は `geo(j,3,k)` (向きの単位ベクトル)、
  BEND の ANGLE は sinc の式、四極などは `qdtrans`・`tftmatu` の写像を使う。pysad は数値微分で、その差は fit の収束位置の差に出るだけで、
  どちらも許容内で止まる。直訳は精度を上げる作業で、ドリフトと BEND の部分だけなら小さいが、四極などの部分は `qdtrans` の移植が要る。優先度は低い (必要になったら)。
- (10:41 に実装) `match6d.tdgeo_term` に tdgeo.f の直訳 (DRIFT の L、BEND の ANGLE (ROTATE 付きも)・L) を入れ、`tffsqu` の幾何の条件 (GX GY GZ CHI1..3) の行で使う。
  幾何の条件だけの行で、変数が DRIFT の L・BEND の ANGLE・L のときは解析微分、四極などの変数 (DX・DY は `qdtrans`・`tftmatu` が要る) と、Twiss の条件が同じ行にある列は中心差分のまま。
  DRIFT の L と BEND の ANGLE は幾何の有限差分と 1e-5 以内で一致し、BEND の L は K64 が出射の向き `geo(:,3,k+1)` で近似する式のまま (有限差分とは O(ANGLE^2) 違う)。試験 `test_tdgeo.py`。
  幾何の条件と Twiss の条件が同じ行にある列は数値微分のままで、収束位置の差は残る。


- 旧 WAKE コマンド (tfwake.f): K64 の tffsa.f:972 は tfwake.f を呼ばず WAKECOMMAND[] (init.n の AutoLoad、Packages/WAKECOMMAND.n は配布に無い) を評価する。gs は fileopen のエラーを出して行の残りを捨てる。pysad は同じ経路 (エラー文の表示だけ未再現)。tfwake.f の移植は不要と判断 (2026-11-16)。
- MessageName の関数形 (2026-11-16 確認): gs は `MessageName[f,"usage"]` を f 無視で `Undefined MessageName General::"usage"` とするが、pysad は `f::usage` と同じ表 (SAD の MessageName.n の定義) を引いて本文を返す。`f::usage` の形は一致。見送り (パーサが `::` の形と関数形を区別する必要がある)。
- PolyLog の非整数次数で |z|≥150: 未実装 (|z|<150 は一致)。
- POL の MULT (2026-11-16): K1・K2・K3 とずれのある MULT のリングは gs と 1e-11 で一致 (test_spin.py)。tmulke の bsir0 の更新も入れた (lattice.py の kick_cx、emit6d.py の _tmulte。この例では効かない)。未解決: K3=300 でずれ 1e-4 のように非線形が強いリング (local/wake/pm1.sad) は、pysad の CALC (qcod) の閉軌道が大きい値 (y=-0.033) になり Emittance が NaN。gs は小さい軌道で通る (gs の CALC は ZX=0 を返す)。qcod の初期値・後退の違いの可能性。空洞・SOL・追跡のスピンは未着手。
- POL の SOL (2026-11-16): `emit_extra._tturne1` に SOL 領域 (`_solmap`) を通す分岐を足し、SOL を含む Emittance[] で偏極時間・スピンチューン・偏極ベクトルが gs と一致 (test_spin.py の test_solenoid_ring_matches_gs)。K64 は SOL (BOUND≠0) の長さを周長に入れない (omega0・MomentumCompaction・偏極時間に効く。PLANAR + SOL で周長 45.6 のときと偏極時間が 15 桁一致) ので `emit6d._emit_circ` で同じにした。GEO のない SOL では K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分) が再現しない。SOL の領域内の MULT・BEND・CAVI のスピンは未確認。

## POL の TrackParticles (スピン追跡、2026-11-16)

- `track6.py` (Python の追跡) に粒子ごとのスピン (`_Particle.s`, `bsi`) と `photons.f` の `sprot` (tradk の中の歳差と放射偏極) を直訳した。
  POL かつ RAD が ON のときだけスピンが動く (K64 の tradk の中で計算するため、RAD OFF ではスピンは変わらない)。`fn_track` は POL ON・RAD ON で Python の追跡に切り替え、
  入力 9〜11 列 {x,px,y,py,z,dp, 角度, sy/|s|, flag}、戻りの 7・8 列は atan(sz, sx) と sy/|s| (tftrack.f の ktfresetparticles と tfsetparticles)。損失 (flag) の処理は高速追跡 (TrackRun) にしかないので、POL のスピン追跡では粒子が落ちない。
- bsi の入る場所: tsolqur (四極のスライスの両端)、tbendcore (入口・出口の y/rhob)、tbendi (K1 のある BEND の各スライス)、tmulti (tsolqu の ibsi と cx の項)、tdrift の ak0 付き、tdrift_solenoid (bz·al)、bz≠0 の tsolqur。tsolrot の chi3 のスピン回転 (四極の ROTATE)。
- gs と一致 (4e-12 以下、NOFLUC・FLUC とも): 平面リング、ずれ・回転した四極と偏向・6 極のリング、MULT (K1 K2 K3 とずれ)、K1 のある BEND、SEXT・薄い MULT (L=0)、SOL 領域内のドリフト (BZ=0.05)。`test_spin_track.py`。
- 追跡の途中で見つけた POL と無関係の差 (RAD ON、スピンなしでも出る。未調査):
  1. BEND に F1=0.1 FRINGE=3 を付けた追跡 (相対 0.7)。DX・DY・ROTATE だけなら 1e-13 で一致。
  2. L≠0 の CAVI に DX・DY (相対 1e-3)。
  3. SOL (BZ≠0) 領域内の四極 + RAD の追跡。bench_level9.sad (保存した gs の出力、10-14) は pysad と一致するが、現在の gs (10-19 ビルド) は同じ入力で別の値 (S R T1 P1 の x が -.0027867898 から -.0014539780) を返す。NORAD は両方で一致。保存した出力と現在の gs のどちらが正しいかは未確認。
- SOL の長さ: K64 は SOL (BOUND≠0) の長さを周長 (周回周波数・空洞の位相) に入れない。追跡の `circ` (fn_track) も `_emit_circ` を使うようにした。

## RADLIGHT (2026-11-16)

`TrackParticles` は RADLIGHT のとき RAD・TRPT を ON にし、戻りを `{{ls または ld, 座標}, 軌道}` にする (軌道は粒子ごとに {t, x, y, z} の 4 本の表)。K64 の tlinit / tlstore / tlresult を `ffscmd._RadLight` に直訳した。
gs と一致 (差 0、2 粒子・複数周・始点 ls≠1: test_radlight.py)。K64 の癖 (gs で確かめた): ① 記録するのは各区間 (tturn) の始めと、L≠0 の BEND 以外 (型 < CAVI) の要素の後だけで、BEND と曲げのある MULT の中は記録しない。
② tlinit は geo0 の軸 1・2 を設定せず (common の初期値 0)、横方向の座標は軌道に入らない。進行方向 (geo の z 軸) に沿った点だけが並ぶ。③ 位置は straight (dphi=0) の積み上げで、BEND で向きが変わらない。
RadiationField は実装済み (gs 照合は 458 式の特殊関数とは別に test_radfield 相当)。BEND を記録しない K64 の仕様のため、`example.oide.sad` の流れ (B を ANGLE 0 の STEER にする) でないと意味のある場にならない。
- SOL 領域内の BEND (直線、ANGLE か K0 のある): K64 の tsol.f は `tdrift(al, bzs, phiy, phix, …)` を呼ぶので、その直訳 (`track6._tdrift_kick` を bz≠0 に一般化) を試したが、gs は偏向をほぼ無視する (BZ=0.5 T の領域で ANGLE=0.1 でも粒子の位置・角度が 1e-5 しか変わらない。pysad の直訳は ANGLE·距離だけ曲がる)。
  原因は未確認 (tpara で K0/ANGLE が領域内の BEND に入らないなど K64 側の癖かもしれない)。gs と合わないので領域内の BEND は引き続き未対応 (NotImplementedError)。生の gs の出力は local/wake/ps23.sad〜ps28.sad。MULT・CAVI の領域内は未着手。
