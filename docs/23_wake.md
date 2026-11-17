# 23. ウェイク場 (WakeFunction、LWAKE / TWAKE) の追跡と光学 (WAKEOPT)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/work_wake.md`)

対象: `pysad/wake.py` (K64 の `txwake.f` の直訳、`spsort` を含む)、`pysad/ffscmd.py` の `_wake_table` (`tffswake.f` の `tffssetupwake`)、
`pysad/track6n.py` の `turn` / `TrackRun.set_wake` と CAVI の分割 (`k_cav_in` / `k_cav_mid` / `k_cav_out` / `cav_full`)、
`pysad/track6n_ext.py` の `mult_full` (RF のある MULT、`tmultiacc`)、
光学の wake: `pysad/wakeopt.py` (`tffswake.f` の直訳)、`pysad/ffscmd.py` の `initial_orbits` / `ini_optics`、`pysad/match6.py` の `_optics_ini`、
試験 `pysad/tests/test_wake.py` (23 件)、`pysad/tests/test_wake_optics.py` (9 件) (2026-11-16 の収集数)。解説: docs/38_wake_guide.md。

## 1. 使い方 (SAD の台本そのまま)

```
TRPT; INS; LWAKE; TWAKE;           (縦だけなら LWAKE、横だけなら TWAKE)
PBUNCH = 1E10; NBUNCH = 1;         (バンチあたり粒子数、バンチ数。CHARGE も使う)
WakeFunction[Longitudinal, "WK"] = {{z1, wl1}, ..., {zn, wln}};      z: m、w: V/C
WakeFunction[Transverse,   "WK"] = {{z1, wt1}, ..., {zn, wtn}};      w: V/C/m
r = TrackParticles[{1, {x, px, y, py, z, dp, flag}}, ld, nt];
```

- 効くのは TrackParticles だけ (SADHelp の記述どおり)。`TRPT` かつ `LWAKE` か `TWAKE` が ON のときだけ有効 (tftrack.f:120)。光学 (CALC / GO) は §7。
- 要素名は `elname` の形 (同名が複数なら `QF.1` のように番号付き)。名前ごとに表を `WakeFunction[…, 名前]` として評価し、数の 2 列のリストでなければ「その要素にウェイクなし」。
  パターン (`WakeFunction[Longitudinal, name_] := …`) でもよい。
- 要素 l にウェイクがあれば、要素の前に半分のキック (init=True) を、要素の後に残りの半分を与える (tturn.f:256-269, 562-575)。
  `fw = |charge| e PBUNCH NBUNCH / amass / NP × 1/2`、`NP` は TrackParticles に渡した粒子数 (flag=0 の粒子を含む)。

## 2. 計算の中身 (txwake.f)

1. 粒子を z の昇順に並べ (`spsort`、§3)、`txwdefslice` で z の間隔と NBUNCH から切片 (スライス) に分ける。NBUNCH=1 のときは 1 バンチを約 √np 個の切片に分ける。
2. 切片ごとに粒子数 `ws`、x・y・z の平均を取る。両端に仮の切片 (z の幅 / np だけ外側) を足す。
3. 切片の組 (n, m) について、距離 `dz = zs(m) − zs(n)` で表を線形補間し、縦は `wz(n) += W(dz) ws(m)`、横は `wx(n) += Wt(dz) xs(m)`、`wy(n) += Wt(dz) ys(m)`。dz=0 で表の最初の点が 0 以上のときは自己項に 1/2 を掛ける。
4. 各粒子は両隣の切片の値を線形補間して、`g` (dp)、`px`、`py`、`dv` を更新する。
5. 前後で z の並びを作り直すのは、前半 (init) か、前に並べたときの z の幅が 1e-5 (相対) 以上変わったときだけ。

要素のずれ DX は、要素の座標系 (x−DX) で横の wake を計算してから戻す (TENT / TEXIT)。

CAVI の中では tcav.f のとおり分けて呼ぶ: 入口のフリンジ → [wake 前半] → ドリフトと電圧キック (ndiv 分割ごとに wake) → [wake 後半] → 後半のドリフトと出口のフリンジ。dx、dy、theta は 0 を渡す。
VOLT≠0 の MULT (TRPT か RFSW のとき) も `tmultiacc` が同じ形で分割ごとに wake を呼ぶ (§2b)。

SOL の領域 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで) の内側の要素と終端の SOL は、tsol.f のとおり要素ごとに前後へ wake を与える (dx=DX、dy=DY、theta=ROTATE)。
普通の要素 (tturn.f) は dy に ROTATE を入れる癖がある (§3)。

## 2a. TRPT の加速 (CAVI の vnominal)

TRPT で VOLT≠0 の CAVI は、基準運動量が空洞で上がる (`vnominal = VOLT/amass |charge| sin(−PHI sign(charge))`)。
追跡では `tcav.f` のとおり、`ndiv` 分割 (重み `ws`、`dzn`)、空洞の後の `g`、`dv` の付け替え、基準運動量 (`p0`、`h0`、`crad`、`urad`) の更新をする
(`track6n.cav_full`、`emit6d._Phys.set_pgev`)。下流の要素は新しい `p0` で進む。区間 (tturn) の終わりで元の運動量に戻る (tturn.f:146)。
wake の有無によらず効く。`phis=0`、`offset=sin(dphis)` (TRPT) も入れた。位相は `(PHI + DPHI) charge` (tcav.f:33)。

## 2b. RF のある MULT (tmulti.f tmultiacc)

`VOLT + DVOLT ≠ 0` の MULT は、(TRPT か RFSW) のとき `tmulti1` が `tmultiacc` を呼ぶ。そうでなければ RF なしの `tmulti` (従来の `k_mult`)。
pysad は `track6n_ext.mult_full` (ANGLE=0、F1・FRINGE なしの MULT。AUTOPHI は未対応):

- CAVI と同じ `ndiv` 分割と `ws` (加速があるとき)。1 分割ごとに、wake → ソレノイド・四極 (`tsolqu`、K0 と K1) → 放射 (`tradk`) → n≥2 の多極キック → `dzn` と電圧キック。
- 電圧キックの定数と式は tcav と少し違う: `oneev = 1+1e-6`、`pmin = 1e-10` での下限、`tlim = 1e4`、`v20a = v02a = v (w (1/2p0 + 1/2pe))² / 4` (`vcorr`)、`v1a = v11a = 0`。
- フリンジ (`DISFRIN=0`): 入口で `tcavfrin` → `ttfrins`、出口で `ttfrins` の逆順 → `tcavfrin(−v)`。
- vnominal≠0 (TRPT) のとき、末尾で基準運動量を `pe` に付け替える (`cav_full` と同じ)。
- MULT の DVOLT、DPHI、AUTOPHI は parser が属性に付ける。

## 3. K64 の癖 (そのまま写した。gs で確認したもの)

- 追跡では**表の最後の 1 行が使われない**。`twxiwp` の `lwl0 = (ブロック長 − 2) / 2` が行数より 1 少ない (光学の `tffswake` は `(ブロック長 − 1) / 2` で行数どおり)。
  最後の行だけ正の z に値があるような表では、その側のウェイクが消える。2 粒子の実験と、最後の行だけ変えた比較で確認した (`test_wake_last_table_row_is_unused`)。
  使う行が 2 行未満になる表は、SAD では添字 0 を読むので、pysad ではエラーにする。
- `dv` の更新で、`h1` をキック前の `1+g` から作る。そのまま写した。
- 普通の要素の wake は、`tturn.f:263-264` が `dy = ROTATE` (実際の DY は捨てる)、`theta` は未初期化の局所変数 `rot` (実測 0)。pysad は dy に ROTATE (QUAD などは属性名 `ROTATION`) を入れ、theta=0。
  ROTATE=0.3、−0.2 の QUAD で gs と一致した (`rotate_on_wake_elements`)。領域内の SOL の theta=ROTATE は、横の wake が x・y に同じ係数なので結果に効かない。
- 粒子の z の整列は `spkick.f` の `spsort` (不安定なクイックソート)。同じ z の粒子の並びが結果に出るので、`wake._spsort` に直訳した (同じ z が 5 個ずつの試験で一致、`equal_z_ordering`)。
- 光学の縦のキック (`tffswakekick`) は `ut1(mfitddp,:) = ut0(mfitddp,:) − cp0 wbuf` と配列全体への代入。`wbuf` は nfam1..nfam の np 個で、軌道の数 n が奇数だと nfam1 = 1−nfam なので、
  要素の対応が 1 つずれる (軌道 s は軌道 s+1 のキックを受ける。最後の軌道は配列の外を読むが、実測では無視できる大きさ)。pysad も同じにしている (`_Kicker.shift`)。

## 4. 検証 (追跡)

SAD K64 (WSL の gs) と同じ台本を走らせて比べた。粒子 200〜400 個の成分ごとの総和と 3 粒子の値が、10 桁一致する。

| ケース | 内容 |
| --- | --- |
| 縦だけ、横だけ、縦と横 | マーカー 1 個、PBUNCH = 1E10 と 1E12 |
| 2 か所のウェイク | QUAD (DX≠0、同名 2 個なので `QF.1`) と MARK に別の表 |
| NBUNCH = 3 | 3 つのバンチ列 (z を 0.5 m ずつずらす) |
| 粒子損失 | LOSSAMPL = 0.0025 で 107 個が落ちる (np が途中で減っても一致) |
| flag=0 の粒子 | 7 の倍数を最初から落とす。NP = 全粒子数 (400) が最良。wake なしでも出る差 (初めから落ちた粒子の出力座標、3e-4) 以外の差はない |
| 同じ z の粒子 | 5 個ずつ同じ z (`spsort` を写して 0 差) |
| ROTATE≠0 の要素 | QF ROTATE=0.3、QD ROTATE=−0.2 (DX、DY も) |
| CAVI | VOLT=0 で L=0 と L=0.4 (フリンジあり)、DPHI≠0 |
| 加速のある CAVI | VOLT=1E5 (L=0)、3E7 (L=0.6、フリンジあり・なし)。wake なしと wake あり (厚い空洞と薄い空洞)。ndiv>1 になる |
| RF のある MULT | TRPT の加速、薄い MULT、K0・K1・K2・K3、DX・DY・ROTATE、DVOLT、DPHI、FREQ、フリンジ。wake なしと wake あり |
| MULT (VOLT=0) の wake | MULT に wake |
| SOL の領域 | 内側の QUAD (ROTATE=0.2) と終端の SOL にウェイク。wake なしの領域の追跡も一致 |
| 表の最後の行 | 最後の行を変えても結果が変わらない |

## 5. 未対応・未確認

- **元の SAD K64 が SIGSEGV で落ちる 3 つの場合**。原因を特定して SAD のソースを直し (§8、`history/sad_patch/wake_sigsegv_fixes.patch`)、pysad は直した SAD と同じ結果を出す。元の gs とは比べられない。
  1. SOL 領域の始点の SOL にウェイクがある。
  2. 光学の wake で、縦と横の表の片方しかない要素。
  3. TRPT で、追跡の最初の `tapert` の時点で LOSSAMPL の外にいる粒子がある (wake・SOL と無関係。§8)。この場合の出力座標の差は **解消済み** (2026-11-16 夜): K64 の `tapert` は落とす印に `kptbl(i,4)=l_track` を書き、`l_track` はセッション最初の TrackParticles では 0 (CALC をしても変わらない) なので、
     最初の `tapert` では印が 0 (「生きている」) になって落ちず、座標は 3·LOSSAMPL に頭打ちされたまま追跡が続き、後の `tapert` で落ちる。2 回目の TrackParticles からは `l_track≠0` で最初の要素で落ちる。
     pysad は `TrackRun.l_track0`・`k_apert` の lcur==0 で再現する。gs (1 回目・2 回目の呼び出し) と直した SAD (LOSSAMPL=1e-4、3 粒子、2 周、NORAD、TRPT、15 値) と全桁一致 (`test_small_items.py` に 1 件)。
     元の gs は SIGSEGV で落ちるので、比べた基準は直した SAD。
- **WAKE コマンドの旧式の入力 (`tfwake.f`) は移植しない**。K64 の `tffsa.f:972` は `tfwake.f` を呼ばず `WAKECOMMAND[]` (`init.n` の AutoLoad。`Packages/WAKECOMMAND.n` は配布にない) を評価する。
  gs は fileopen のエラーを出して行の残りを捨てる。pysad も同じ経路で、エラー文の表示だけ未再現 (2026-11-16 に確認。SADHelp にない旧機能)。
- `txwdefslice` で、バンチの間隔がバンチ数 NBUNCH を超えて空くときは Fortran が配列の外へ書く (`mb(nb+1)`)。pysad は読める値として扱う。結果が不定なので比較できない。
- RAD ON の追跡は、wake の有無によらず gs と 1e-6 程度ずれる (wake なしでも 8.9e-7、同じ量が wake ありでも出る)。wake とは別の、TrackParticles の放射 (TRPT) の差。
- MULT の AUTOPHI、CAVI の AUTOPHI (位置の z が要る)。`NotImplementedError`。MULT の F1・FRINGE。`mult_full` (RF のある MULT の wake 付き経路) は ANGLE=0 だけで、ANGLE≠0 の MULT そのものの追跡は docs/28 (`multa.py`) にあるが、wake・RF との併用は未確認。
- wake と COORD・TCAVI (docs/25・26)・UND・MAP の併用、SOL 領域内の TCAVI は確認していない。
- 並列追跡 (NPARA) は wake があると 1 本になる (tftrack.f:133)。pysad は元から 1 本。

## 6. 試験の実行

```
py -3.12 -m pytest pysad/tests/test_wake.py pysad/tests/test_wake_optics.py -q
```

SAD が無くても通る (参照値は gs の出力を試験に書き込んである)。

## 7. 光学の wake (WAKEOPT) と InitialOrbits

```
TRPT; INS; LWAKE; TWAKE; WAKEOPT;
WakeFunction[Longitudinal, "WK"] = {...};  WakeFunction[Transverse, "WK"] = {...};     (片方だけでも pysad は動く。元の SAD は落ちる: §8)
InitialOrbits = {{x1, px1, y1, py1, z1, dp1}, ...};      (6 個: 中心の軌道からのずれ)
InitialOrbits = {{AX, BX, NX, AY, BY, NY, EX, EPX, EY, EPY, R1..R4, DETR, DX, DPX, DY, DPY, DZ, DDP, AZ..ZPY}, ...};   (28 個: Twiss も与える)
CALC;   (または GO)
```

### InitialOrbits とは

光学計算 (CALC / GO) が、中心の軌道 (idp=0) のほかに、指定した初期条件の軌道を何本か持つ仕組み (tffsa.f の `tffsinitialcond`)。
軌道は中心の軌道 (閉軌道ではなく INS の入口の軌道) からのずれで与える。6 個ならずれだけ、28 個 (`ntwissfun`) なら Twiss も与える (BX>0 のとき使う。軌道は中心の軌道に足す)。
InitialOrbits があると、オフモメンタム整合と有限振幅整合は使わない。

- n 個の軌道は idp = nfr1..nfr (0 を除く) に順に入る。nfr=(n+1)/2、nfr1=−(n/2)。0 は中心。nfam1 は n が偶数なら −nfr、奇数なら 1−nfr。
- 表示 (CALC / GO の後): `DP`、`Orbit ID` (InitialOrbits の何番目か。中心は 0)、`Res.` (軌道ごとの残差) の 3 行と、各関数の軌道ごとの値 (`tfshow.f`)。列は残差の大きいものから nn 本を選び、dp の昇順に並べる。
- `Twiss[...]` は中心の軌道の値。他の軌道の値は表と、FitValue フック (条件が FITP でその軌道に掛かるとき) で読む。

### WAKEOPT のときの計算 (tffswake.f)

1. 全軌道 (中心を含む) を INS で、ウェイクのある要素の入口まで伝搬する (軌道ごとに `ins_twiss`。軌道の z、δ もそのまま追う)。
2. 全軌道の z・x・y から、各軌道の dpx、dpy、ddp へのキックを与える (`wakeopt._Kicker.kick`、`tffswakekick` の直訳)。
   縦は dz = z_i − z_j で表を引く (z_i − z_j ≥ −sigg の滑らかな境界、sigg = max(1e-5, SIGZ/np/NBUNCH))、横は dz + 0.3989 sigg。
   `cp0 = |charge e PBUNCH / np / (2 gammab amass)|`、np は軌道の数 (nfam−nfam1+1)。表が同じで z の配置がほぼ同じ (0.1 sigg 以内) のとき、前の値を使い回す条件も同じ。
3. 要素を伝搬し、要素の出口でもう一度キック (要素の前後で 1 回ずつ。`txwake` の前後の半分とは違い、光学は同じ大きさを 2 回)。
4. 次のウェイクまで伝搬する。

WAKEOPT が OFF、または LWAKE・TWAKE が OFF、またはウェイクのある要素がないときは、各軌道は独立に伝搬する (qcell1 の fam=.true. と同じ。中心の入口の Twiss と、中心の入口の軌道 + ずれ)。

### 検証 (光学)

gs と同じ台本の CALC の出力を比べた (`test_wake_optics.py`): 表は文字列で完全一致、中心の軌道 (`Twiss["DX"..]`) は相対 2e-8 以内 (z は絶対 3e-14 以内)。

| ケース | 内容 |
| --- | --- |
| 2 軌道 (縦+横) | `FIT $$$ DX DDP` の表。軌道ごとの AX BX NX .. 全部の列が一致 |
| 3 軌道 (縦だけ、横だけ、縦+横) | 奇数 (nfam1=1−nfam) の縦のキックの添字のずれ |
| 4、5 軌道 | 偶数・奇数、QUAD (DX・DY) の wake |
| NBUNCH=3、SIGZ=0.01、2 か所の wake | 表の使い回し (z がほぼ同じ) の分岐 |
| 28 個形式の軌道 | Twiss (BX>0) と NX の引き継ぎ |
| WAKEOPT なし | InitialOrbits だけ (独立な伝搬) |

### 既知の差と限界

- 軌道の z は、INS の `ins_twiss` が δ の飛行時間の項 (`dvemit`、約 δ L/γ²) を足さないため、gs と 1e-14〜1e-10 (δ と長さによる) ずれる。この差の影響は表と x, px の一致 (1e-9 台) に出る程度。
  wake と無関係な INS の軌道の差 (DP0≠0 の INS の Twiss["DZ"] で再現)。
- CELL (周期系) と InitialOrbits の併用、TRPT で加速要素があるときの InitialOrbits は未対応 (`SADError`。SADHelp も wake は INS だけとしている)。
- 光学の `gammab(le)` は基準運動量で一定 (加速要素なし)。
- gs の表は軌道が多いと列を選ぶ (`nn`)。pysad も同じ規則だが、極端な数 (nn<nfam−nfam1+1) は未確認。

## 8. 元の SAD K64 の SIGSEGV と修正

SAD を `-O0 -g` で作り直し (WSL に複製、`~/sad_dbg`)、落ちた場所を `addr2line` で特定した (`setarch -R` で PIE の基準を固定)。`-fcheck=bounds` は SAD の初期化 (`tfstk.f:487` の `klist(0)`) で落ちるため使えない。

| 症状 | 落ちる場所 | 原因 | 修正 |
| --- | --- | --- | --- |
| 光学の wake で表が片方だけ | `tffswake.f:102` | `iwl` か `iwt` が 0 (表なし) でも `nlw=(ilist(1,iwl-1)-1)/2` を計算して `ilist(1,-1)` を読む | 0 のとき、もう片方の表を指すダミーを使う (キックは `iwl /= 0` で守られている) |
| SOL 領域の始点にウェイク | `tturn.f:567` → `twxiwp` | `tsol` は領域全体を進めて戻る。その間 `sol` の分岐が `aewak` を古いままにして、領域の各要素で後半キックを重ねて呼び、`nwak` が表の数を超える。始点の後半キックを `tsol` の後で呼ぶため、領域内のウェイクの `nextwake` も古いままになる | `tsol` の前に始点の後半キックと `nwak` の更新を済ませ、`aewak` を偽にする |
| TRPT で最初から LOSSAMPL の外の粒子 | `tapert.f:184` | `idtypec(l_track)` を評価するが、最初の `tapert` は要素の前 (`l_track=0`) | `l_track>0` のときだけ評価する |

3 つ目は元の SAD の別の不具合で、wake も SOL も要らない (LOSSAMPL=0.0009 で、wake なしでも落ちる)。「SOL 領域内で粒子が落ちると落ちる」のではなく、LOSSAMPL=0.0009 の台本は初期の振幅がすでに外だった。
LOSSAMPL=0.002 (領域内・後で落ちる) は元の SAD でも落ちず、pysad と 10 桁一致する。

直した SAD の結果は、pysad の追跡と一致した: SOL 領域の始点の wake (単独と、領域内の wake との併用) で 0 差、光学の縦だけの表の中心軌道で一致。
pysad はこの 2 つを、直した SAD の動作として実装している (領域の始点の SOL は普通の要素と同じ前後のキック、光学は無い表の側をキックしない)。

解説: docs/38_wake_guide.md (WAKE)
