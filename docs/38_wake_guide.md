# 38. WAKE (ウェイク場) 解説: WakeFunction、LWAKE / TWAKE、WAKEOPT

最終更新: 2026-11-17

初めて読む人向けの解説。詳細な検証表と SIGSEGV の修正は [`docs/23_wake.md`](23_wake.md)、経緯は [`docs/history/work_wake.md`](history/work_wake.md)。
対象コード: [`pysad/wake.py`](../pysad/wake.py) (txwake.f)、[`pysad/wakeopt.py`](../pysad/wakeopt.py) (tffswake.f)、[`pysad/ffscmd.py`](../pysad/ffscmd.py) (`_wake_table`、`_wake_opt_tables`、`fn_track` の wake 部、`_wake_cmd`)、
[`pysad/track6n.py`](../pysad/track6n.py) (`turn`、`TrackRun.set_wake`、CAVI の分割)、[`pysad/track6n_ext.py`](../pysad/track6n_ext.py) (`mult_full`)。

## 1. 何ができるか

SAD の「ウェイク場」は、要素ごとに与えた縦・横のウェイク関数で、バンチ内の粒子が互いに (先行する粒子から後続の粒子へ) 受けるキックを追跡に入れる機能。次の 2 通りがある。

| 用途 | 有効にする | 効く所 |
| --- | --- | --- |
| 粒子追跡 | `TRPT; LWAKE (縦); TWAKE (横)` | `TrackParticles` だけ (tftrack.f:120 で `(TWAKE or LWAKE) and TRPT`) |
| 光学 | 上に `WAKEOPT` と `InitialOrbits` | `CALC` / `GO` (INS のみ。tffsa.f:1186 `wake = (TWAKE|LWAKE) & TRPT & WAKEOPT`) |

どちらも TRPT (粒子運動の座標系) が前提。`LWAKE`・`TWAKE`・`TRPT`・`WAKEOPT` はスイッチで、`NO…` で OFF。

## 2. 使い方 (動く SAD スクリプト)

### 2.1 追跡 (`pysad/tests/test_wake.py` の W3_SRC。gs と 10 桁一致)

```
MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"]   = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"]  = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"]    = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];   px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];   py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];   fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
```

読み方:
- `WakeFunction[Longitudinal, 名前]` は縦 (z は m、w は V/C)、`WakeFunction[Transverse, 名前]` は横 (w は V/C/m)。表は `{{z1, w1}, ..., {zn, wn}}`、z の昇順。
- 名前は `elname` の形。同名の要素が複数あるときは `QF.1` のように番号付き (上の例は QF が 2 個あるので最初のものが `QF.1`)。パターン (`WakeFunction[Longitudinal, name_] := ...`) でもよい。
- 表が数の 2 列のリストでなければ、その要素にウェイクなし。要素に表がなくても他の要素に付ければよい。
- `PBUNCH` (バンチあたり粒子数)、`NBUNCH` (バンチ数)、`CHARGE` を使う。`TrackParticles` に渡す粒子数 NP が係数に入る。
- `TrackParticles[{1, {...}}, 17, 1]` の第 2 引数は終点の要素の位置 (要素数以下)。大きすぎると gs は `Component number must be positive` で失敗する。

試験では、このスクリプトの出力 (成分 1〜6 の総和、粒子 1・77・400 の値、生き残り数) を gs の値と相対 2e-8 以内で比べている (`test_wake_matches_sad`)。

### 2.2 光学 (`test_wake_optics.py` の A1_SRC。gs と一致)

```
(格子は 2.1 と同じ要素に WK を足したもの)
TRPT; INS; NORAD; NOFLUC;
LWAKE; TWAKE; WAKEOPT; PBUNCH=1E12; NBUNCH=1;
WakeFunction[Longitudinal, "WK"] = {{-0.2, 1e9}, {-0.05, 2e9}, {0, 3e9}, {0.05, 1e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK"]   = {{-0.1, 1e12}, {0.0, 3e12}, {0.05, 1e12}, {0.1, 0}, {0.3, 0}};
InitialOrbits = {{1e-3,1e-4,2e-3,0,-0.05,1e-3},{-1e-3,0,1e-3,2e-4,0.05,0}};
FITP 2;
FIT $$$ DX 0 DDP 0;
CALC;
Do[Print[k, " ", Twiss["DX",k], " ", Twiss["DPX",k], " ", Twiss["DY",k], " ", Twiss["DPY",k], " ", Twiss["DZ",k], " ", Twiss["DDP",k]], {k, 1, 18}];
```

`InitialOrbits` の 1 組は中心の軌道からのずれ `{x, px, y, py, z, dp}` (6 個) か、Twiss も与える 28 個形式 (`{AX, BX, NX, AY, BY, NY, EX, EPX, EY, EPY, R1..R4, DETR, DX, DPX, DY, DPY, DZ, DDP, AZ..ZPY}`)。
`CALC` の後の `Twiss[...]` は中心の軌道の値。他の軌道の値は表示 (`DP`、`Orbit ID`、`Res.` の 3 行と軌道ごとの列) と FitValue フックで読む。

## 3. 物理と K64 の式

### 3.1 追跡のキック (txwake.f)

1. 粒子を z の昇順に並べ (`spsort`)、`txwdefslice` で NBUNCH と z の間隔から切片に分ける。NBUNCH=1 なら 1 バンチを約 √np 個の切片に分ける。
2. 切片ごとに粒子数 `ws` と x・y・z の平均を取る。両端に仮の切片を足す。
3. 切片の組 (n, m) の距離 `dz = zs(m) − zs(n)` で表を線形補間。縦は `wz(n) += W(dz) ws(m)`、横は `wx(n) += Wt(dz) xs(m)`、`wy(n) += Wt(dz) ys(m)`。
4. 各粒子は両隣の切片の値を線形補間して、`dp`、`px`、`py`、`dv` を更新する。
5. 係数 `fw = |charge| e PBUNCH NBUNCH / amass / NP / 2`。要素にウェイクがあれば、要素の前に半分 (init=True)、後に残り半分 (tturn.f:256-269, 562-575)。要素のずれ DX は TENT / TEXIT で座標系を移して計算する。

### 3.2 光学のキック (tffswake.f)

軌道 (中心を含む n 本) を INS でウェイクのある要素の入口まで伝搬 → 全軌道の z・x・y から各軌道の dpx、dpy、ddp へキック (`tffswakekick`) → 要素を通す → 出口で再度キック → 次のウェイクまで伝搬。
追跡と違い前後で同じ大きさのキックを 1 回ずつ。縦は `dz = z_i − z_j` (境界は sigg = max(1e-5, SIGZ/np/NBUNCH) で滑らかにする)、横は `dz + 0.3989 sigg`。`cp0 = |charge e PBUNCH / np / (2 gammab amass)|`。
WAKEOPT・LWAKE・TWAKE のどれかが OFF か、ウェイクのある要素がなければ、各軌道は独立に伝搬する。

## 4. 要素の種類ごとの扱い

| 要素 | 追跡での扱い (K64 のソース) | pysad |
| --- | --- | --- |
| 普通の要素 (DRIFT、QUAD、BEND、MARK など) | tturn.f が前後で半分ずつ `txwake` | `track6n.turn` が前後に呼ぶ |
| CAVI | tcav.f が中で分割 (`ndiv`) ごとに呼ぶ。入口フリンジ → wake 前半 → ドリフトと電圧キック (分割ごとに wake) → wake 後半 → 出口フリンジ。dx、dy、theta は 0 を渡す | `k_cav_in` / `k_cav_mid` / `k_cav_out`、`cav_full`。`_wake_table` は dx=dy=theta=0 |
| MULT で VOLT+DVOLT≠0 (TRPT か RFSW のとき) | `tmultiacc` が CAVI と同じ形で分割ごとに wake | `track6n_ext.mult_full` (ANGLE=0、F1・FRINGE なしの MULT)。AUTOPHI は v2 モードだけ対応 |
| MULT で VOLT=0 | 普通の要素 | 普通の要素と同じ (gs と一致を確認) |
| SOL の領域内 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで) と終端の SOL | tsol.f が要素ごとに前後へ wake (dx=DX、dy=DY、theta=ROTATE)。領域内の CAVI・VOLT≠0 の MULT は tmultiacc の中で分割ごとに wake (ewak) | `_wake_table` が `in_sol` で判別して同じ形にする。**領域内の CAVI・RF つき MULT の wake は v2 モードだけ** (k64 モードは領域内の CAVI・MULT が未対応) |
| SOL 領域の始点の SOL | 元の K64 は SIGSEGV (後述 6 節) | 直した SAD と同じ動作 (普通の要素と同じ前後のキック) |
| BEAMBEAM | 要素にウェイクを付けられる | 縦・横・5 周で gs と一致 ([`docs/24_beambeam.md`](24_beambeam.md)) |

TRPT で VOLT≠0 の CAVI は、基準運動量が空洞で上がる (`vnominal`)。`cav_full` が `tcav.f` のとおり `ndiv` 分割と `p0`、`h0` の付け替えを行い、区間 (tturn) の終わりで元に戻す。wake の有無によらず効く。
光学の wake は、TRPT で加速要素があるときは未対応 (`SADError`)。

## 5. 旧 WAKE コマンド

Fortran には旧式の `WAKE` コマンド (`tfwake.f`) があるが、K64 の `tffsa.f` はこれを呼ばず、SADScript の `WAKECOMMAND[]` を評価する (`init.n` が `Packages/WAKECOMMAND.n` を AutoLoad)。
この配布には `WAKECOMMAND.n` が無い。gs で確かめた動作は、ファイルが無いと fileopen のエラーを出し、`WAKECOMMAND[]` が True にならないので行の残りを捨てる、というもの。
**k64 モード**の pysad の `_wake_cmd` は同じ経路を取る (`WAKECOMMAND.n` があればそれを評価する)。SADHelp にもない旧機能。

**v2 モード**は `tfwake.f` の構文を実装した (B66)。`WAKE` から `ENDWAKE` までの表を WakeFunction の表として登録する (CAVI だけ):

```
WAKE LONGITUDINAL RF
  -0.1   2e9          ! s [m]  w [V/C] (横 TRANSVERSE は V/C/m)
   0     5e8    0.02 -1e9
ENDWAKE;
```

- 表は行をまたいで書ける (`ENDWAKE` か `;` が来るまで続く)。1 行に何組でも、注釈 `!` 可、`ENDWAKE` の後の `;` は任意、同じ行の続き (`ENDWAKE; CALC;`) も実行する。大文字小文字・省略形 (`wake t RF`) 可。要素名はワイルドカード可 (合う CAVI すべてに登録)。値は数か式。点数の上限は無い (400 点で確認)。
- 追跡は `WakeFunction[…] = …` で同じ表を与えた場合と同じ結果 (gs の参照値と 2e-9、400 点・2 台の CAVI で 1e-12)。TRPT と LWAKE / TWAKE が ON の `TrackParticles` で効く。
- 使えない: 表の途中の `;` (エラー)、CAVI 以外の要素、k64 モード。`ENDWAKE` を書き忘れると次の `;` までを表として読んでエラーにする。
- 例 [`example/v2_C/C2w_wake_table.sad`](../example/v2_C/C2w_wake_table.sad)、試験 [`test_v2_C_wakecmd.py`](../pysad/tests/test_v2_C_wakecmd.py)・[`test_v2_C2_wakecmd.py`](../pysad/tests/test_v2_C2_wakecmd.py)。

## 5a. wake の重み (v2 モードで直した K64 の誤り)

K64 の `tcav.f`・`tmulti.f` は `fw=fw0*(wsn+ws(n-1)*.5d0)` で、空洞が 2 つ以上に分割されると wake の蹴りの合計が (1.5n − 0.5)/n 倍になる。**k64 モードは再現し、v2 モードは台形の重み fw0·(ws(n) + ws(n−1))/2 (合計 1) にする。** 位置によらない表では、エネルギーの変化が空洞の電圧 (= 分割数) によらなくなる (20 MV の 1 分割と 200 MV の 4 分割以上で 5e-3。k64 は 1.25〜1.5 倍)。分割が 1 のときは両モードで同じ。領域の外の空洞・RF つき MULT も。docs/23 §2c。

## 6. K64 の癖 (そのまま写したもの)

- 追跡では**表の最後の 1 行が使われない**。`twxiwp` の `lwl0 = (ブロック長 − 2)/2` が行数より 1 少ない (光学の `tffswake` は行数どおり)。最後の行だけ正の z に値がある表では、その側のウェイクが消える。
  gs の 2 粒子の実験で確認 (`test_wake_last_table_row_is_unused`)。使う行が 2 行未満になる表は、SAD は添字 0 を読むので、pysad はエラーにする。
- 普通の要素のウェイクは `dy = ROTATE` (実際の DY は捨てる)、theta は未初期化の局所変数 (実測 0)。ROTATE≠0 の QUAD (ROTATE=0.3、−0.2) で gs と一致。
- 粒子の整列 `spsort` は不安定なクイックソート。同じ z の粒子の順序が結果に出るので `wake._spsort` に直訳した (argsort では 2.4e-7 ずれた)。
- `dv` の更新で `h1` をキック前の `1+g` から作る。
- 光学の縦のキックは `ut1(mfitddp,:) = ... − cp0 wbuf` と配列全体への代入で、軌道の数が奇数だと要素の対応が 1 つずれる (軌道 s は s+1 のキックを受ける)。pysad も同じ (`_Kicker.shift`)。
- 元の K64 が SIGSEGV で落ちる 3 つの場合 ([`docs/21_sad_bug_report.md`](21_sad_bug_report.md) の B31-B33、修正の差分は [`docs/history/sad_patch/wake_sigsegv_fixes.patch`](history/sad_patch/wake_sigsegv_fixes.patch)):
  光学で縦と横の表の片方しかない要素 (`tffswake.f:102`)、SOL 領域の始点の SOL にウェイク (`tturn.f:567`)、TRPT で追跡の最初から LOSSAMPL の外にいる粒子 (`tapert.f:184`。wake と無関係)。
  pysad は前 2 つを直した SAD と同じ動作にしている。元の gs とは比べられない。

## 7. pysad の構成

- `ffscmd._wake_table`: 要素ごとに `WakeFunction[Longitudinal|Transverse, 要素名]` を評価し、`{位置: (縦, 横, dx, dy, theta)}` を作る (LWAKE・TWAKE が OFF の側は None)。CAVI と VOLT≠0 の MULT は dx=dy=theta=0、SOL 領域は DX・DY・ROTATE、普通の要素は dy=ROTATE。
- `ffscmd.fn_track`: `wake_on = (LWAKE or TWAKE) and TRPT` のとき表を作り、`kick_factor` で `fw` を出して `tr.set_wake(WakeKick, 表)`。高速追跡 (TrackRun) で扱えない要素の並びなら `SADError`。
- `wake.WakeKick.apply`: txwake の本体 (numba)。`track6n.turn` が要素の前後で呼ぶ。
- 光学: `ffscmd.initial_orbits`、`ini_optics`、`_wake_opt_tables` (LWAKE/TWAKE の ON/OFF はキック側で見る) から `wakeopt.wake_orbits`。`match6._optics_ini` が InitialOrbits を使う。

## 8. gs との一致の度合い

数値は [`docs/23_wake.md`](23_wake.md) §4・§7 と試験から。使った条件: 追跡はエネルギー 1.7 GeV 前後の小さなライン (INS、NORAD NOFLUC)、粒子 200〜400 個、PBUNCH=1E10〜1E12。実際のリングでの検証はしていない。

- 追跡: 成分ごとの総和と 3 粒子の値が 10 桁一致。ケースは縦だけ・横だけ・両方、2 か所のウェイク、NBUNCH=3、粒子損失 (LOSSAMPL=0.0025 で 107 個)、flag=0 の粒子、同じ z の粒子、ROTATE≠0、CAVI (VOLT=0 と加速あり、フリンジあり・なし)、RF のある MULT (K0〜K3、DX・DY・ROTATE、DVOLT、DPHI、FREQ)、MULT (VOLT=0)、SOL 領域。
- 光学: 条件の表は文字列で完全一致、中心の軌道 (`Twiss["DX"..]`) は相対 2e-8 以内 (z は絶対 3e-14 以内)。2〜5 軌道、NBUNCH=3 と表の使い回し、28 個形式、WAKEOPT なしを含む。
- 試験: `py -3.12 -m pytest pysad/tests/test_wake.py pysad/tests/test_wake_optics.py -q` で 32 通過 (本文書の作成時に確認)。SAD が無くても通る (参照値は試験に書き込み済み)。

## 9. 限界・未確認

- **RAD ON の追跡**は、wake の有無によらず gs と 1e-6 程度ずれる (wake なしでも 8.9e-7)。wake とは別の TrackParticles の放射の差で、wake 自体の RAD ON での一致は確認できていない。
- 光学の中心の軌道の z は、INS の `ins_twiss` が δ の飛行時間の項 (約 δL/γ²) を足さないため、gs と 1e-14〜1e-10 ずれる (wake と無関係。DP0≠0 の INS で再現)。
- CELL (周期系) と InitialOrbits の併用、TRPT で加速要素があるときの InitialOrbits は未対応 (`SADError`)。光学の `gammab` は基準運動量で一定。
- 光学の表は軌道が多いと列を選ぶ (`nn`)。極端な軌道数は未確認。
- MULT の AUTOPHI と CAVI の AUTOPHI は、k64 モードでは `NotImplementedError`、v2 モードでは実装した (docs/44)。ANGLE≠0 の MULT (tmulta)、MULT の F1・FRINGE は両モードで未対応。
- TRPT で最初から LOSSAMPL の外にいる粒子があると、元の gs は落ちる (B33)。直した SAD (`~/sad_dbg`) との出力座標の差は **解消** (2026-11-16 確認)。原因は K64 の `tapert` が落とす印に `kptbl(i,4)=l_track` を書くこと: `l_track` はセッション変数で、セッション最初の TrackParticles では 0 (CALC をしても変わらない)。最初の tapert では印が 0 =「生きている」になり落ちず (座標は 3·LOSSAMPL に頭打ちされて追跡を続け、後の tapert で落ちる)、2 回目の TrackParticles からは `l_track≠0` で最初の要素で落ちる。pysad は以前は常に最初の要素で落としていたので `TrackRun.l_track0`・`k_apert` の lcur==0 で再現し、gs の 1 回目・2 回目の呼び出しと直した SAD (LOSSAMPL=1e-4、3 粒子、NORAD、TRPT、2 周、15 値) と全桁一致 ([`test_small_items.py`](../pysad/tests/test_small_items.py))。元の K64 はここで落ちるので、元の K64 に対する基準は無い。
- `txwdefslice` で、バンチの間隔が NBUNCH を超えて空くときは Fortran が配列の外へ書く。pysad は読める値として扱うので結果が不定で、比較できない。
- 並列追跡 (NPARA) は wake があると 1 本になる (tftrack.f:133)。pysad は元から 1 本。
- 旧 WAKE コマンド (`tfwake.f`) は k64 モードでは移植していない (v2 モードは 5 節)。
- SOL 領域内の CAVI・RF つき MULT の wake は v2 モードで実装した (縦の wake によるエネルギー変化の合計が領域の外の同じ空洞と bz → 0 で 1e-6 で一致、横の蹴りの合計は 2% 以内)。**横の蹴りを gs と比べてはいない**。k64 モードでは領域内の CAVI・MULT そのものが未対応。TCAVI と wake の併用は未確認 ([`docs/26_tcav.md`](26_tcav.md))。
- `track6n.py` の `segment` には SOL 領域を含む区間で `NotImplementedError('wake: SOL 領域を含む区間')` を出す分岐が残るが、**当たらない防御**と確認した (2026-11-16)。Prog は SOL 領域でも 1 要素 1 kind (`compile_sol_region` が要素ごとに足す) なので、kinds の数と b-a は常に等しい。

## 10. 関連文書

- [`docs/23_wake.md`](23_wake.md): 検証表、SIGSEGV の調査と修正の詳細、WAKEOPT の仕様。
- [`docs/history/work_wake.md`](history/work_wake.md): 実装記録。[`docs/21_sad_bug_report.md`](21_sad_bug_report.md) (B31-B33)。
- [`docs/08_tracking.md`](08_tracking.md) (TrackParticles)、[`docs/24_beambeam.md`](24_beambeam.md) (BEAMBEAM と wake)、[`docs/26_tcav.md`](26_tcav.md)、[`docs/29_sol_optics.md`](29_sol_optics.md)、[`docs/11_sad_quirks.md`](11_sad_quirks.md)。
- 参照: SADHelp (CERN) の WakeFunction。
