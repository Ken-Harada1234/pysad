# 2026-11-16 の作業 (話題別)

最終更新: 2026-11-16

時系列ではなく話題別に書く。各項目に根拠のファイル (コード、試験、記録、gs の出力) を付ける。未解決の差の一覧は `work_unimpl.md`、現状の解説は `docs/34_spin_csr.md` ほか。
gs は WSL の K64 ビルド (`bin/gs`、2026-10-19 ビルド)。`local/wake/` は作業者のローカルの作業フォルダ (リポジトリに含まれない)。

## 1. POL のスピン追跡 (TrackParticles)

- 内容: `track6.py` (Python の追跡) に粒子ごとのスピン (`_Particle.s`、`bsi`) と `photons.f` の `sprot` (tradk の中の歳差と放射偏極) を直訳した。POL かつ RAD ON のときだけスピンが動く (K64 の tradk の中で計算するため。RAD OFF ではスピンは変わらない)。
  `ffscmd.fn_track` は POL ON・RAD ON で Python の追跡に切り替える。入力は 9〜11 列 {x,px,y,py,z,dp, 角度, sy/|s|, flag}、戻りの 7・8 列は atan(sz, sx) と sy/|s|。
- bsi を入れた場所: tsolqur (四極のスライスの両端)、tbendcore (入口・出口の y/rhob)、tbendi (K1 のある BEND の各スライス)、tmulti (tsolqu の ibsi と cx の項)、tdrift の ak0 付き、tdrift_solenoid (bz·al)、bz≠0 の tsolqur。tsolrot の chi3 のスピン回転 (四極の ROTATE)。
- 一致の度合い: gs と 4e-12 以下 (NOFLUC・FLUC とも)。平面リング、ずれ・回転した四極と偏向・六極のリング、MULT (K1 K2 K3 とずれ)、K1 のある BEND、SEXT・薄い MULT (L=0)、SOL 領域内のドリフト (BZ=0.05)。根拠: `pysad/tests/test_spin_track.py`。
- 限界・未確認: 損失 (flag) の処理は高速追跡 (TrackRun) にしかないので、POL のスピン追跡では粒子が落ちない。SOL 領域内の MULT・BEND・CAVI のスピンは未確認 (領域内の要素そのものが未対応)。
- Emittance・EMIT のスピン: MULT (K1・K2・K3 とずれのあるリングで gs と 1e-11。tmulke の bsir0 の更新も入れた: `lattice.py` の `kick_cx`、`emit6d.py` の `_tmulte`)、空洞、SOL (`emit_extra._tturne1` に SOL 領域の分岐) まで。偏極時間・スピンチューン・偏極ベクトルが gs と一致 (`test_spin.py`、SOL は `test_solenoid_ring_matches_gs`)。
- 根拠のファイル: `pysad/spin.py`、`pysad/track6.py`、`pysad/emit_extra.py`、`pysad/tests/test_spin.py`、`test_spin_track.py`、`work_unimpl.md` (「POL の MULT」「POL の SOL」「POL の TrackParticles」)、`docs/34_spin_csr.md`。
- 使った条件: 小さな試験ラティス (平面リング、BZ=0.05 の SOL など)。エネルギーは試験の既定 (試験ファイルを参照)。実用の大きなリングは未確認。`Polarization` は K64 が NaN (B55) で、pysad も NaN。

## 2. RADLIGHT

- 内容: `TrackParticles` は RADLIGHT のとき RAD・TRPT を ON にし、戻りを `{{ls または ld, 座標}, 軌道}` にする (軌道は粒子ごとに {t, x, y, z} の 4 本の表)。K64 の tlinit・tlstore・tlresult を `ffscmd._RadLight` に直訳した。
- 一致の度合い: gs と差 0 (2 粒子・複数周・始点 ls≠1)。根拠: `pysad/tests/test_radlight.py`。RadiationField・RadiationSpectrum は実装済み (`pysad/radfield.py`、`test_radfield.py`)。
- K64 の癖 (gs で確かめた): ① 記録するのは各区間 (tturn) の始めと、L≠0 の BEND 以外 (型 < CAVI) の要素の後だけで、BEND と曲げのある MULT の中は記録しない。② tlinit は geo0 の軸 1・2 を設定せず (common の初期値 0)、横方向の座標は軌道に入らない (進行方向の z 軸に沿った点だけが並ぶ)。③ 位置は straight (dphi=0) の積み上げで、BEND で向きが変わらない。
- 限界: BEND を記録しない K64 の仕様のため、`example.oide.sad` の流れ (B を ANGLE 0 の STEER にする) でないと意味のある場にならない。
- 根拠のファイル: `pysad/ffscmd.py` (`_RadLight`)、`test_radlight.py`、`work_unimpl.md` (「RADLIGHT」)、`docs/43_radlight_guide.md`。

## 3. MULT のスピンと SOL のスピン・周長

- MULT: §1 の Emittance・EMIT のスピンの項を参照 (K1・K2・K3 とずれのあるリングで gs と 1e-11、`test_spin.py`)。
- 未解決: K3=300、ずれ 1e-4 のような非線形が強い MULT リング (`local/wake/pm1.sad`) は、pysad の CALC (qcod) の閉軌道が大きい値 (y=-0.033) になり Emittance が NaN。gs は小さい軌道で通る (gs の CALC は ZX=0 を返す)。qcod の初期値・後退の違いと疑うが未確認。
- SOL の周長: K64 は SOL (BOUND≠0) の長さを周長 (周回周波数・空洞の位相・MomentumCompaction・偏極時間) に入れない。PLANAR + SOL で周長 45.6 のときと偏極時間が 15 桁一致する。`emit6d._emit_circ` で同じにし、追跡の `circ` (`fn_track`) にも使った。
  GEO のない SOL では K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分、MomentumCompaction が 100 倍) が、pysad は GEO を見ないので再現しない。
- SOL の領域内の BEND (直線、ANGLE か K0 のある): K64 の tsol.f は `tdrift(al, bzs, phiy, phix, …)` を呼ぶので直訳 (`track6._tdrift_kick` を bz≠0 に一般化) を試したが、gs は偏向をほぼ無視する (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない。pysad の直訳は ANGLE·距離だけ曲がる)。原因は未確認 (tpara で K0/ANGLE が領域内の BEND に入らないなど K64 側の癖かもしれない)。合わないので未対応 (NotImplementedError)。gs の出力は `local/wake/ps23.sad`〜`ps28.sad`。MULT・CAVI の領域内は未着手。
- 根拠のファイル: `pysad/emit6d.py` (`_emit_circ`)、`test_spin.py`、`work_unimpl.md`、`docs/27_sol_optics.md`、`docs/42_sol_guide.md`。

## 4. 旧 WAKE コマンド

- 旧式の WAKE コマンド (tfwake.f) は移植しないと判断した。K64 の `tffsa.f:972` は tfwake.f を呼ばず `WAKECOMMAND[]` (`init.n` の AutoLoad。`Packages/WAKECOMMAND.n` は配布に無い) を評価し、gs は fileopen のエラーを出して行の残りを捨てる。pysad も同じ経路 (エラー文の表示だけ未再現)。
- 根拠: K64 のソース `src/tffsa.f`、`Packages/init.n`、gs の出力、`work_unimpl.md`、`docs/21_wake.md`、`docs/36_wake_guide.md`。

## 5. tlum・tluma (measure.py)

- `pysad/measure.py` に K64 の tlum・tluma を直訳した。MEASURE が gs で動かない (B42) ので、gs との比較はできず、解析値との比較だけ (`pysad/tests/test_measure.py`)。MEASURE/TRCOD (`tracke`・`ttinit`・`ttstat`・`trcoda`・`trackb`) は gs の参照が取れず、使うスクリプトも見つからないので翻訳していない (docs/19 B41・B42)。
- 根拠のファイル: `pysad/measure.py`、`test_measure.py`、`work_unimpl.md`。

## 6. CSR の確認

- 実装 (gs と照合): CSRTrack・CSRHaissin・CSRSetupOY・CSROYMatrix・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRConvert (K64 の tfcsrvec と同じ)・CSRInit。`Packages/CSR.n` の MakeZL・SkinDepth・PathDifference・StandardKList・VariableMesh・ZThre などは gs と 1e-10〜1e-13 で一致。
- 限界・未確認: ScanTrack・OYStability などの長い流れは未照合。遅い試験 (MakeZL、約 7 分) は `PYSAD_SLOW=1` のときだけ走る。K64 の癖は docs/19 B55〜B60。
- 根拠のファイル: `pysad/csr.py`、`pysad/tests/test_csr.py`、`docs/34_spin_csr.md`、`docs/38_csr_guide.md`。

## 7. MOMENTUM の既定

- 初期の MOMENTUM は 30 GeV (K64 の `initb1.f`)。以前の文書は 1 GeV と書いていたので直した。エネルギーを書かないスクリプトで効く。
- 根拠: K64 の `src/initb1.f:1698` (`RsetGL('MOMENTUM',30.d9,idummy)`) と `pysad/ffscmd.py:376, 464` (既定 30e9)。gs での実測は未実施 (ソースの読解とコードの確認のみ)。

## 8. クラスのメンバー

- クラスのメンバーでない名前の `obj@name[...]` は式のまま (gs に合わせた)。評価器は `pysad/sadclass.py`・`sadscript.py`。GUI の試験は `PYSAD_SAD_PACKAGES` が無いと skip になるので、評価器を変えたら設定して全試験を回す (全体実行では GUI 系が数回に 1 回落ちる。単独なら通る)。
- 同日に gs に合わせた SADScript の小さな差: `Skip`・`Read` の型 (`test_evalfixes.py`)、数の書式、連鎖比較、`Null[]` の展開、With・Block の引数の局所化 (`work_unimpl.md`、`docs/32`・`35`)。`MessageName[f,"usage"]` の関数形は不一致のまま見送り。

## 9. そのほかの同日の追加

- `SynchroBetaEmittance` (`pysad/synchrobeta.py`): 12 セルのリング + CAVI (VOLT=2e6 HARM=64) で gs と emitx 相対 4e-7〜1e-5 (許容 3e-5: `test_synchrobeta.py`)。emity は 1e-18 の雑音。未確認: 実用の大きなリング、不安定なチューン、`emiout`、結合、RFSW の組み合わせ。
- tdgeo (`match6d.tdgeo_term`): DRIFT の L と BEND の ANGLE・L だけ解析微分 (他は中心差分)。`test_tdgeo.py`。Twiss の条件が同じ行にある fit は数値微分のままで、収束位置が gs と 1e-6 の程度違う。
- INS を含むリングの追跡・Emittance[] は gs と一致 (追加の実装は不要)。DELCOR の CAVI (B62)、PolyLog・HurwitzZeta の複素の領域 (docs/31)、CSR の MakeZL、GaussianCoulombFitted (docs/30)。
- 試験の全体 (昼の時点): 511 通過・3 skip (514 件収集)。回帰 `QUICK=1 run_regression.py` は FAIL 0 (遅い 4 件を SKIP して 46 件)。夜の修正後は 520 通過・3 skip (523 件収集)、50 ベンチすべて実行して FAIL 0 (§11)。

## 10. 見つけた未解決の差

1. BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)。DX・DY・ROTATE だけなら 1e-13 で一致。スピンなしでも出る (POL と無関係)。未調査。
2. L≠0 の CAVI に DX・DY を付けると相対 1e-3 ずれる。未調査。
3. SOL (BZ≠0) 領域内の四極 + RAD の追跡で、`bench_level9.sad` (保存した gs の出力、2026-10-14) は pysad と一致するが、現在の gs (2026-10-19 ビルド) は同じ入力で別の値 (S R T1 P1 の x が -.0027867898 から -.0014539780) を返す。NORAD は両方で一致。どちらが正しいか未確認。
4. SOL 領域内の BEND・MULT・CAVI は未対応 (§3)。
5. K3 が強くずれの大きい MULT リングで pysad の CALC (qcod) の閉軌道が gs と違う解になり Emittance が NaN (§3、`local/wake/pm1.sad`)。
6. 既存の差: `emit6d._tquade` (Python 版) が x のずれた軌道で δ の微分を落とす (1e-5)。`work_unimpl.md` 「後でまとめて直す一覧」1。

根拠: `work_unimpl.md`、`00_RESUME.md` §6・§8。夜の確認 (§11) で 3・5 の原因を絞り込み、B33 の差・UNIPOL・RELW などは解決した。

## 11. 夜の確認 (未確認の項目を確かめた、2026-11-16 夜)

文書の「未確認」「未調査」の項目を gs (WSL の K64 ビルド) と pysad で実際に動かして確かめた。記録は作業者のローカルの `local/wake/verify_log.md` (リポジトリに含まれない)。台本の名前 (v1.sad など) も `local/wake/` にある。
使った条件はどれも小さな試験ラティス (平面リング、12 セルや 24 セルなど)。実用の大きなリングは未確認。

### 11.1 確認できたもの (差なし)

- `Spline[{{x,y},..}]` (自然スプライン)・`Sort[list, Greater]`・`Scan` は gs と全桁一致。`Spline[data,0,0]` は gs も pysad も未評価のまま (v1.sad)。`Polarization2` は gs の `Emittance[]` が返さず、pysad も同じ (v2.sad)。
- POL のスピン (Emittance[]): ROTATE+DX+DY 付き空洞 (結合リング。偏極時間 13 桁、スピンチューン 1e-14、偏極ベクトル 1e-13: v3a.sad)、ANGLE≠0 の MULT (偏極時間 14 桁、スピンチューン 1e-14: v3b.sad)、BEND の DROTATE+DX (偏極時間 14 桁、偏極ベクトル 1e-14: v3c.sad)。
- BEAMBEAM を含むリングの CALC・Emittance (Twiss のチューン 1e-14、運動量収縮 4e-12。BEAMBEAM は Marker 扱い: v4.sad)。STURN・COUPLE (BEAMBEAM) は K64 のソースでキーワードの登録 (`initb1.f`) 以外に参照がなく、使われない。
- PhotonList の 12 列は `photons.f` の `tphotonlist` の定義と一致。`track6n` の「wake: SOL 領域を含む区間」の分岐は、Prog が SOL 領域でも 1 要素 1 kind なので当たらない防御。
- RADLIGHT の軌道: geo[0] から局所 z 軸の大域成分は (1,0,0) で、進行方向の距離は軌道の 2 行目 (大域 x) に出る。`docs/34` の「geo の z 軸」は正しい。RADLIGHT → RadiationField → FindPeak → RadiationSpectrum は pysad で最後まで通る (軌道 8 点、スペクトル {3,1001})。gs は RadiationField で SIGSEGV (B49) なので通し比較はできず、軌道の部分だけ gs と一致。
- 回帰 QUICK が SKIP していた 4 件 (7・12t・13a・16a) を実行し、すべて PASS (7: 5.0e-11、12t: 3.5e-11、13a: 4.3e-21、16a: 1.4e-31)。50 ベンチすべて実行して FAIL 0 (`reg_slow4.txt`)。
- SynchroBetaEmittance の大きなリング (24 セルの平面リング + CAVI VOLT=2e6、AzimuthalModes->7): emitx は gs と相対 1.4e-7〜1.4e-6 (4 点)。emity は雑音 (sb11.sad)。
- ListPlot・Plot のオプション (PlotRange、Scale、PlotLabel、GridLines、Frame、Legend、FontScale、Prolog・Epilog・Background、Plot->False、AspectRatio) の戻りの Graphics 式は gs と文字列で完全に一致 (10 例: pl2・pl3.sad。窓の見た目ではなく式の比較)。`Options[g]` は gs が全オプションを返し pysad は {} を返す (差のまま。pl1.sad)。
- Twiss と LINE[…] の端数位置 (3.5、4.25、4.5 (BEND の内側)、リスト指定): 平面の 12 セルのリングで 1e-14 一致 (fr1.sad)。`docs/05` の「端数位置の幾何は未確認」を確認済みにできる。
- DELCOR で OCT を対象に含める: チューンと QF の DX が 12 桁で一致 (OC 自体の DX は gs も pysad も 0 で、OCT の効果を見る強い比較ではない: dc1.sad)。COORD の対 (DX DY DZ CHI3 付き) を含むリングの FIT (NX・NY=0.2、QF K1・QD K1 を FREE、解析微分の経路) は gs と全桁一致 (cf1.sad)。COORD のパラメータを FIT 変数にする場合は未確認。
- `Optics.n` の TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo は gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し pysad は 0: op1.sad)。

### 11.2 差を見つけて pysad を直したもの

- **UNIPOL・RELW** (`docs/18` は「効果の有無は未確認」): gs は UNIPOL で変数の符号を変えず (`tffsmatch.f:669`)、RELW で重みが変わり、FIT の結果が変わる (fit1.sad: QF K1・QD K1 の 2 変数で NX・NY を 0.2 に fit。既定 K=−.0068645845/−.0057261776、UNIPOL K=.3289522375/−.3239879365、RELW K=−.0024047842/−.0107636499)。pysad は `match6` に bipol・absweit の実装はあるがフラグとつながらず、どれも既定の結果だった。`ffscmd.matcher_` でつなぎ、3 通りとも gs と全桁一致。試験: `test_item4.py` に 3 件。
- **DUMP**: 誤差入りのリングで値は一致したが、DK・DDK 列の書式が違った (gs は Fortran の `1p g12.5`: 0 → `0.0000`、小さい値 → `4.22829E-04`。pysad は `0`・`0.00042283`)。`ffscmd._gwd` を足して一致 (dm2.sad、`test_ffsextra.py` の `test_dump_format_matches_gs`)。
- **MATRIX コマンド**: 名前だけの引数 (gs は NAME.1 として扱う) を位置にできず、ラティスに無い名前は gs が読み飛ばす (to=$$$) のを再現していなかった。直した (mt1・mt2.sad、`test_ffsextra.py`)。
- **LDSolve.n (複素行列)**: `DiagonalMatrix` が複素数の列を受けず失敗していた (`float() argument must be … not 'complex'`)。あわせて、行列の演算 (Dot・Inverse) の結果は gs が要素ごとに虚部 0 の複素数を実数に戻すので、pysad も要素ごとに戻すようにした (以前は全要素が実のときだけ)。LDSolve の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (ld2・ev1.sad、`test_evalfixes.py` に 2 件)。`PeriodicR` は呼び方が不明で比較できず (gs も引数エラー)。
- **セッション最初の l_track=0** (B33 の差の原因。§11.3 と `work_wake.md`): K64 の tapert は落とす印に `kptbl(i,4)=l_track` を書く。l_track は tmacro のセッション変数で、セッション最初の TrackParticles では 0 (CALC をしても変わらない) なので、最初の tapert (tturn の最初の要素) では印が 0 = 「生きている」になり落ちない。座標は 3·LOSSAMPL に頭打ちされて追跡を続け、後の tapert で落ちる。2 回目の TrackParticles からは最初の要素で落ちる。pysad は常に最初の要素で落としていた。`TrackRun.l_track0`・`k_apert` の lcur==0 で再現し、gs (1 回目・2 回目の呼び出し) と直した SAD (b33.sad、15 値) と全桁一致 (lm2〜lm4.sad、`test_small_items.py` に 1 件)。`docs/36` の「B33 の出力座標が合わない」は解消。
- **OpenRead の失敗の gs 互換**: 入力ファイルが無いとき、失敗の出方が OpenRead の FileNotFoundError だった。gs と同じに、開けなければメッセージを出して入力行を打ち切り、`Read[h,…]`・`Close[h]` は wrongtype にした (rd3・rd4.sad、`test_evalfixes.py`)。

### 11.3 原因を絞り込んだもの (未修正)

- **Tkinter のウィジェットが Widget 値にならない** (`tk1〜tk5.sad`、`pv1.sad`): gs では `e=Entry[w]` のあと e は `Widget[e,Entry[w]]` を値に持ち、Insert・Get・Delete が働く (Entry の内容 AB→CD、Listbox 2 項目、`v=e[Get]`=AB)。pysad では e が裸のシンボルのままで、`e[Insert]=…`・`e[Get]` は何も呼ばない (Tcl 呼び出し 0 回、エラーもなし)。ウィジェット自体は作られる (winfo children に Entry・Listbox)。`w=Window[]` の直後も、gs は `Head[w]=Widget` だが pysad は Symbol で、ウィジェットの親が w でなく Root になる (winfo では '.' の直下。gs は Toplevel の下)。
  pysad は `Tk$DefCount=1`・`Root$WindowWidget=w`・`TkPathName[Root$Window]="."` まで進むが、`TkPathName[w]` が "" で、最初の Window (root) の `a$:=Widget[a$,w[master$]]` (`Tkinter.n` の DefineWidget 内) が w の値にならない。root の Window は toplevel を作らない。`f[x_]:=(x=5); f[aa]` は pysad も aa=5 で、パターン変数への代入は問題ない。
  玩具の再現 (With と `x:=Widget[x,…]`、UpValue の `(a_=W[opt___])^:=…`) は、UpValue の分岐までは gs と同じに動くが、`x:=Widget[x,…]` 自体は gs も pysad も再帰になり、実物がなぜ gs で再帰しないか (Widget の属性・評価の順) は未特定。pysad は再帰で RecursionError を出して落ちる (gs は General::deep のメッセージ)。Insert の位置に "End" (大文字) を渡すと gs は Tcl エラー 1501 (Tcl の添字は "end")。
- **SOL 領域内の四極 + RAD の gs の食い違い** (bench_level9 の差。b9a〜b9d.sad): 単一粒子 (粒子の並びの癖が効かない) の 1 周追跡 (2.5 GeV、RAD NOFLUC、emit のあと) で、SOL の BZ を 0 → 0.01 T に変えた。
  (1) 領域内がドリフトだけ (`SB D1 SE`): gs = pysad (BZ=0・0.01 とも 9 桁以上)。(2) 領域内が四極 (`SB QD SE`): BZ=0 は gs = pysad。BZ=0.01 で gs は x=6.33e-4、y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4、y=1.1e-8 (連続)。
  (3) `SB D1 QD SE` BZ=0.01: gs の結果が四極なしの `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていないように見える。BZ=0.01 と 0.2 でも gs の x はほぼ同じ (2.95e-5。BZ=0 の −2.23e-4 とは大きく違う) ので、BZ に連続でない。
  (4) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作り直した `~/sad_dbg`) は同じ値 (x=2.954816161E-5、10 桁一致) で、最適化の違いではなく現在の K64 ソースから作った gs の振る舞い。10-15 の古い実行ファイル (`obj/AMD64-Linux5/sad1.exe`) は symbol lookup error で動かず、保存済みの参照 (10-14) を出した版とは確認できない。
  結論: pysad が壊れているのではなく、現在の gs が SOL 領域内の四極 + RAD を物理的に誤って扱っている疑いが強い (pysad は物理と保存済みの参照に合う)。ソースのどこかは未特定。
- **B33 の差**: 元の K64 は SIGSEGV で参照にならず、直した SAD (`~/sad_dbg` の gs5) は最初から外にいる粒子を最初の tapert で外さず追跡を続ける (b33.sad: LOSSAMPL=1e-4、3 粒子 x0={1e-3,1e-6,2e-3}、NORAD、TRPT、2 周)。原因は §11.2 の l_track=0 の扱いで、解決した。

### 11.4 ソースの読解で確かめたもの (動かしてはいない)

- SPCH 要素: `tturn.f:521` の `case (icSPCH)` が PSPAC フラグのとき `tpspac` を呼ぶ。pysad は SPCH を読むだけで物理は未翻訳 (PSPAC・WSPAC・SPAC フラグも効果なし)。PROF は K64 のソースに現れない。LOSSMAP も未翻訳のまま。
- B29 (QUAD に RAD=1 で COD・エミッタンスが 0): K64 は `krad = (RAD_QUAD == 0) .and. al /= 0` (`tturn.f:360`、`tsol.f:677`、`tturne.f:572`) で、要素の RAD が 0 以外ならその要素の放射を止める。バグでなく仕様 (RAD=1 は「この要素は放射しない」)。
- MON コマンド (B30 相当): `tffsa.f` に MON・MCOD・TCOD のコマンドは無い (grep で 0 件)。FFS では式として評価されて何もしない。
- フラグ GAUSS (MAIN の TRACK・MEASURE の初期分布: `ttinit.f:139`、`track.f`)、CTIME (CPU 時間の表示)、CONV (収束の判定・表示: `tcod.f`・`qtwiss.f`・`temit.f`)、EMIT の EMIOUT (EMIT 報告の追加行: `temit.f:992` ほか。pysad は spin の行だけ対応)。BIPOL は `tffsmatch.f` だけ。
- `CSR.n` の TrackAnimate は ListPlot のあと `System["import …"]` (ImageMagick)・`System["mkdir …"]` を呼ぶので、Windows の pysad (System は実行しない) では画像が作れない。`Wait[mem]` は FreeMem (top) と Sleep を使い、gs はこの環境で top が無く FreeMem が OpenRead エラーになる (csr8.sad)。pysad の FreeMem は psutil で代替して 1.0 を返す。並列を増やしても速くなるかは未確認 (pysad は fork しない)。

### 11.5 確認できなかったもの・限界を発見したもの

- 複数区間 (複数周) の RADLIGHT の軌道は区間の始めの点が重なり (t が同じ点が 2 つ)、pysad の RadiationField は `x must be strictly increasing sequence` で失敗する (rf3.sad)。gs は RadiationField が SIGSEGV なので比較できない。
- 結合のあるラティス (QT・SF・DY のずれ入りの 8+4 セルのリング) の SynchroBetaEmittance: gs は SIGSEGV で落ちるので比較できない。pysad は `ValueError: cannot convert float NaN to integer` で失敗する (エラーの出し方が粗い。sb10.sad)。
- Tkinter の Widget 値と SOL 内の四極 + RAD の食い違いは §11.3 のとおり原因が確定していない。`PeriodicR` は呼び方が不明。COORD のパラメータを FIT 変数にする場合は未確認。
- 全体: 修正後の全試験は 520 通過・3 skip (523 件収集)、非 QUICK の回帰は 50 ベンチをすべて実行して FAIL 0 (205 s)。
