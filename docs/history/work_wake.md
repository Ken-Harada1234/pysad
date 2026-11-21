# WakeFunction (LWAKE / TWAKE) の実装記録

最終更新: 2026-11-16

作業日: 2026-11-14。現状の説明は `docs/23_wake.md`。

## 依頼
GUI の次に wake と beam-beam を実装する (実例の台本なし)。このメモは wake。根拠は SADHelp (CERN)、K64 のソース (`txwake.f`、`tffswake.f`、`tturn.f`、`tcav.f`)、WSL の gs での実測。

## 調べたこと
- 追跡の wake は `tftrack.f:120` で `(twake or lwake) and trpt` のときだけ。表は `tffssetupwake` が要素ごとに `WakeFunction[Longitudinal|Transverse, 要素名]` を評価して作る。数の 2 列のリストなら採用。
- 呼び出しは `tturn.f` (普通の要素: 前後で 1 回ずつ)、`tcav.f` (CAVI: 分割ごと)、`tmulti.f` (VOLT≠0 の MULT)。`tffswake.f` の本体は光学 (WAKEOPT) 用で追跡には使われない。`tfwake.f` は旧式の WAKE コマンド。
- `txwake` は z で並べて切片に分け、切片どうしの距離で表を線形補間する。

## 進め方
1. `pysad/wake.py` に `txwdefslice` と `txwake` を 1 始まりの添字のまま numba で写した。
2. `ffscmd.py` の `_wake_table` で表を取り、`track6n.py` の `turn` に前後のフックを足した。CAVI は `k_cav` を 3 段に分けた (`k_cav` 本体は触らない)。
3. gs で同じ台本を走らせて比べた (`local/wake/`: w1〜w5、c1〜c3、cmp.sh、exp.sh)。台本は先頭に MARK のあるラインで、`D` などの記号は gs では保護されているので別名にする。

## 見つけた差とその原因
- 最初は縦のウェイクだけ gs と合わなかった (横は 10 桁一致)。z が高い粒子は一致、低い粒子だけずれていたので、表の正の z 側の扱いを疑い、正の側だけに値がある表で 2 粒子の実験をしたところ gs は kick が 0 だった。
  原因は `twxiwp` の `lwl0=(ilist(1,iwpl-1)-2)/2` で、表の行数より 1 少ない。最後の行を捨てると全部一致した。
- TRPT の CAVI は wake なしでも gs とずれる。基準運動量が空洞で加速する (`vnominal`) のに、`TrackRun` は一定運動量のため。`turn` には `tcav.f` の `phis=0`、`offset=sin(dphis)` (TRPT) を入れたが、VOLT≠0 では合わない。VOLT=0 の CAVI は一致。

## 続き: 「未対応を対応させて」(2026-11-14)
- **TRPT の加速** (`vnominal`): `tcav.f` の全体を `track6n.cav_full` に写した。`turn` が基準運動量を持ち歩く (`_Phys.set_pgev`)。gs と 10 桁一致 (VOLT=1E5、3E7、厚い・薄い、フリンジあり・なし)。VOLT≠0 の CAVI の wake も一致。
  Cavity の compile の `pr` を 9 要素にした (`vnom`、`vq`)。`track6n_ext` の取り出しは `pr[:7]`。
- **SOL 領域**: 内側の要素と終端の SOL は要素ごとに前後の wake (tsol.f)。gs と一致。始点の SOL にウェイクがあると gs が SIGSEGV で落ちるので `SADError`。
  `WakeKick.apply` に回転 (TENT / TEXIT) を足した。
- (この時点では MULT の VOLT と WAKEOPT は未対応だった。次の節で対応した。)
- 試験を 6 本足した (計 12 本)。全体は 313 通過、6 skip。

## 続き: 「WAKE に必要な機能なら実装して。gs が落ちる以外を実装」(2026-11-14)
- **InitialOrbits** (tffsa.f `tffsinitialcond`) を `ffscmd.initial_orbits` に。6 個形式と 28 個形式。軌道は idp = nfr1..nfr (0 を除く) に入る。`match6.Matcher` は `ini` があると nfr・dp・光学を InitialOrbits から作る (`_optics_ini`)。
- **WAKEOPT** (`tffswake.f`) を `pysad/wakeopt.py` に直訳。軌道ごとに `ins_twiss` で伝搬し、ウェイクのある要素の前後で全軌道の z・x・y から dpx、dpy、ddp のキック (`tffswakekick`)。
  表の使い回し (`wbufxy`、`wzl` と `dsig = 0.1 sigg`) も写した。`FFS.calc(ini=True)` が GO / CALC の後の `Twiss[]` を中心の軌道にする。
- **tfshow の軌道の表** (`DP` / `Orbit ID` / `Res.` の 3 行と、軌道ごとの列) を `print_fit_table` に足した (InitialOrbits のときだけ)。残差の大きい軌道を選び、dp の昇順に並べる。
- gs と比べて分かったこと:
  - 奇数個の軌道で、縦のキックの代入 `ut1(mfitddp,:)=…-cp0*wbuf` が配列全体のため、要素の対応が 1 つずれる (K64 の癖、`_Kicker.shift`)。n=1〜5 の軌道の実験 (定数の縦 wake、z をずらす) で突き止めた。
  - 縦と横の表の片方しかない要素は gs が SIGSEGV (`SADError`)。
  - gs の中心軌道の z は δ の飛行時間 (`dvemit`) を含む。pysad の INS の z は含まない (既存の差。DP0≠0 の INS で再現)。
- **MULT の RF** (`tmulti.f` の `tmultiacc`) を `track6n_ext.mult_full` に。VOLT+DVOLT≠0 の MULT が (TRPT か RFSW) のとき。wake も分割ごとに呼ぶ。gs と 10 桁一致 (厚い・薄い、K0〜K3、DX・DY・ROTATE、DVOLT、DPHI、FREQ、フリンジ、TRPT の加速)。
  parser が MULT の DVOLT・DPHI・AUTOPHI と CAVI の AUTOPHI を落としていたのを直した。CAVI の位相に DPHI を足した (tcav.f:33)。
- **未確認だった項目の確認**: ROTATE≠0 の普通の要素は K64 の癖 (dy=ROTATE、theta=0) どおりで gs と一致。flag=0 の粒子は NP=全粒子数が最良 (wake なしと同じ差だけ)。同じ z の粒子は `spsort` を写して一致 (argsort では 2.4e-7 ずれた)。
  SOL の領域内で粒子が落ちる wake は gs が SIGSEGV。RAD ON は wake と無関係に gs と 1e-6 ずれる。
- 試験: `test_wake.py` 22 本、`test_wake_optics.py` 9 本 (計 31 本)。

## 続き: 元の SAD の SIGSEGV の原因と修正 (2026-11-14)
- 元の SAD (`<SAD>`) は書き換えず、WSL に複製 (`~/sad_dbg`) して `-O0 -g` で作り直した。手順の落とし穴:
  - `HostArch` が今は `AMD64-Linux6` を返すが、ビルド済みは `AMD64-Linux5`。`obj/AMD64-Linux5` で直接 `make exe`、`obj/.../config.mk` の SADDIR と `sad.conf` の `SAD_ROOT` を複製に向ける、`.*_done` の印を touch する (tcl/tk を作り直さない)。
  - `libsad.a` は Fortran の `.o` を中に持つ。`FOPT` を変えたら `*.o` のうち `src/*.f` があるものと `libsad.a`、`sad1.exe` を消す。`sad.conf` の `FOPT=` は 3 行ある (ubuntu 以外も)。
  - `bin/gs` の `SAD_ARCH` を `AMD64-Linux5` に固定した `gs_dbg` を作る。`setarch -R` + `addr2line -f -e sad1.exe (アドレス − 0x555555554000)` で行番号が出る (gdb は無い)。
- 原因は 3 つ。1) `tffswake.f:102` (表なしの添字 0)、2) `tturn.f:567` (SOL 領域で `aewak` が古いまま、`nwak` が表の外)、3) `tapert.f:184` (`l_track=0` の `idtypec`、wake・SOL と無関係)。詳細は `docs/23_wake.md` §8。
- 修正は `docs/history/sad_patch/wake_sigsegv_fixes.patch` (元の SAD に当てる差分)。直した SAD の出力は pysad と一致 (始点 SOL の wake 単独・領域内との併用で 0 差、光学の縦だけの表)。
- pysad: 領域始点の SOL の wake と、光学の片方だけの表を `SADError` にしていたのをやめた (直した SAD と同じ動作)。試験の参照は直した SAD の出力。
- 誤りの訂正: 「SOL 領域内で粒子が落ちると SIGSEGV」としたのは、台本の LOSSAMPL が初期の振幅より小さかった (最初の tapert で落ちる) ためで、SOL とは無関係だった。

## 試験
`pysad/tests/test_wake.py`、`test_wake_optics.py`。参照値は gs の出力 (10 桁)。

## 落とし穴
- gs に lattice 台本を渡すときは、先頭に `FFS;` を付けない (付けると SADScript として読まれ、`QUAD QF=(...)` が式になる)。SADScript (GUI) の台本は逆に `FFS;` が先頭に必要。
- 元の gs は、SOL 領域の始点にウェイク、光学の片方だけの表、TRPT で最初から LOSSAMPL の外の粒子、で SIGSEGV (§ 続き: 元の SAD の SIGSEGV の原因と修正)。
- TrackParticles の第 2 引数は終点の位置 (nlat 以下)。要素数より大きいと gs は 'Component number must be positive' (`TrackParticles[..., 12, 1]` が 8 要素のラインで失敗した)。
- gs に `Total` は無い (`Plus@@`)。
- bash の heredoc での python は文字が壊れる (`\` など)。`chr(92)` を使うか、ファイルに書く。

## 未解決
`docs/23_wake.md` §5 (VOLT≠0 の CAVI と MULT、SOL の中、光学の WAKEOPT、ROTATE のずれ、flag=0 の粒子、バンチ間隔の多い場合)。

## 次の作業
beam-beam (`bbstrhl.f` ほか)。

## 2026-11-16 の確認・追記

- 2026-11-16 追記: 旧式の WAKE コマンド (tfwake.f) は移植しない。K64 の `tffsa.f:972` は tfwake.f を呼ばず `WAKECOMMAND[]` (`init.n` の AutoLoad。`Packages/WAKECOMMAND.n` は配布に無い) を評価し、gs は fileopen のエラーを出して行の残りを捨てる。pysad も同じ経路 (エラー文の表示だけ未再現)。
- 2026-11-16 夜の追記 (B33 の差の解消): TRPT で最初から LOSSAMPL の外にいる粒子 (B33、元の gs は SIGSEGV で、直した SAD で比べていた) の出力座標が pysad と合わなかった原因が分かり、解消した。K64 の tapert は落とす印に `kptbl(i,4)=l_track` を書き、l_track はセッション最初の TrackParticles では 0 なので、最初の tapert では粒子が落ちず (座標は 3·LOSSAMPL に頭打ちで追跡を続け、後の tapert で落ちる)、2 回目の TrackParticles からは最初の要素で落ちる。pysad は常に最初の要素で落としていたので、`TrackRun.l_track0`・`k_apert` で再現した。直した SAD (`~/sad_dbg` の gs5) の 15 値と全桁一致 (b33.sad、lm2〜lm4.sad)、試験は `test_small_items.py` に 1 件。根拠: `work_session_20261116.md` §11.2・§11.3。
- 他の本文 (WakeFunction、WAKEOPT、SIGSEGV の修正) は 11/14 のまま。`wake.py`・`wakeopt.py` は 11/14 以後変更なし。wake の SOL 領域を含む区間の分岐 (`track6n`) は、Prog が SOL 領域でも 1 要素 1 kind なので当たらない防御と確認した (夜)。現状の解説は docs/23 と docs/38。
